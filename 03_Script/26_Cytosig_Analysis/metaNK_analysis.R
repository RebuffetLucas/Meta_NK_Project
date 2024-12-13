library(Seurat)
library(dittoSeq)
library(biomaRt)
library(ggpubr)
library(scales)

setwd('~/Documents/Research/Niklas_projects/5.metaNK/data/')

##### Functions #####
biomaRT_convertID_hgnc <- function(geneIDs) {
  print('Get gene names and descriptions for Ensembl IDs...')
  
  dataset <- 'hsapiens_gene_ensembl'
  
  converted.df <- tryCatch(
    {
      print("Attempting useast mirror...")
      mart <- useDataset(dataset, useMart('ENSEMBL_MART_ENSEMBL', host = 'https://useast.ensembl.org'))
      getBM(filters='hgnc_symbol',
            attributes=c('hgnc_symbol','gene_biotype','chromosome_name'),
            values=geneIDs,
            mart=mart)
    },
    error = function(e) {
      tryCatch({
        print("Failed. Attempting uswest mirror...")
        mart <- useDataset(dataset, useMart('ENSEMBL_MART_ENSEMBL', host = 'https://uswest.ensembl.org'))
        getBM(filters='hgnc_symbol',
              attributes=c('hgnc_symbol','gene_biotype','chromosome_name'),
              values=geneIDs,
              mart=mart)
      },
      error = function(e2) {
        print("Failed. Attempting asia mirror...")
        mart <- useDataset(dataset, useMart('ENSEMBL_MART_ENSEMBL', host = 'https://asia.ensembl.org'))
        getBM(filters='hgnc_symbol',
              attributes=c('hgnc_symbol','gene_biotype','chromosome_name'),
              values=geneIDs,
              mart=mart)
      })
    }
  )
  return(converted.df)
}

bubbleZP <- function(df, celltype_col, value_col, replace = T, lower = -2, upper = 2) {
  celltypes <- unique(df[,celltype_col])
  out <- as.data.frame(matrix(nrow = length(celltypes), ncol = 3))
  colnames(out) <- c('celltype', 'z', 'p')
  out$celltype <- celltypes
  real_means <- c()
  pvalues <- c()
  
  for (aCelltype in celltypes) {
    print(aCelltype)
    real_mean <- mean(subset(df, get(celltype_col) == aCelltype)[,value_col])
    real_means <- c(real_means, real_mean)
    #pvalues <- c(pvalues, pnorm(real_mean, mean(df[,value_col]), sd(df[,value_col]), lower.tail=FALSE))
    pvalues <- c(pvalues, t.test(subset(df, get(celltype_col) == aCelltype)[,value_col], 
                                 subset(df, get(celltype_col) != aCelltype)[,value_col],
                                 alternative = 'greater')$p.value)
  }
  out$z <- real_means
  out$p <- pvalues
  out <- out %>% 
    mutate(logP = -log10(p),
           Z = rescale(z, to = c(0, 1))) %>% 
    mutate(logP = if_else(logP > 10, 10, logP))
  
  return(out)
}

##### 1. Cytosig analysis #####
## 1.1 Merge Seurat objects
# NK13.seu <- readRDS('NK1andNK3.rds')
# NK2.seu <- readRDS('NK_2.rds')
# NK2.seu$SecondClust <- NK2.seu$Secondclust
# NK2.seu$Secondclust <- NULL
# 
# NK.seu <- merge(NK13.seu, NK2.seu)
# NK0.seu <- readRDS('NK_All_AllGenes.rds')
NK.seu <- readRDS('PBMC_V2Chem_VF1.rds') # new clustering by Lucas and Erik

## 1.2 Prepare Cytosig data
# Filter genes, only keep protein-coding genes
counts <- GetAssayData(NK.seu, 'count')
genes <- rownames(counts)
# converted.df <- biomaRT_convertID_hgnc(genes)
# converted.good.df <-  subset(converted.df, ! hgnc_symbol == '' & 
#                                gene_biotype %in% c('protein_coding') &
#                                chromosome_name %in% c(1:22,'X','Y'))
# counts <- counts[converted.good.df$hgnc_symbol,] # remove non-protein-coding genes; retain 10234/11850 protein-coding genes
# counts <- counts[rowSums(counts) >= 100,] # 10230 final genes for analysis
# Convert to log2(TPM + 1) and mean centralize counts as recommended by the author (https://github.com/data2intelligence/CytoSig/issues/2)
size_factor <- 1E6/colSums(counts) # convert count to TPM
norm_counts <- sweep(counts, 2, size_factor, "*")
norm_counts <- log2(norm_counts + 1)
#meanCentral_counts <- norm_counts[rowMeans(norm_counts) > 0.1,] # mean TPM > 0.1 for each gene; further filtering out lowly expressed genes; 9340 genes in the final set!
background <- rowMeans(norm_counts)
meanCentral_counts <- sweep(norm_counts, 1, background, "-")
# Write table
set.seed(12345)
meanCentral_rand_counts <- meanCentral_counts[,sample(ncol(meanCentral_counts))] # shuffle cells first

numOfCellsPerRun <- 10000
for(i in 1:ceiling(ncol(meanCentral_rand_counts)/numOfCellsPerRun)) {
  print(i)
  temp_counts <- meanCentral_rand_counts[,(numOfCellsPerRun*(i-1)+1):(min(numOfCellsPerRun*i, ncol(meanCentral_rand_counts)))]
  write.table(temp_counts, paste0('new_seu_cytosig_mat_part', i, '.txt'), quote = F, sep = '\t')
}

## 1.2 Store Cytosig results (Z-score) to the seurat object
cytosig_Z.list <- list()
for (i in 1:4) {
  print(i)
  cytosig_Z.list[[i]] <- read.csv(paste0('new_cytosig_part',i,'.Zscore'), sep = '\t')
}

all_cells <- str_replace(rownames(NK.seu@meta.data), '-1', '.1')

cytosig_Z.all.res <- do.call(cbind, cytosig_Z.list)[,all_cells]

cytokines <- c('IL15','IL2','IL21','IL12','IL18','IL10','IL27','IFN1','IFNG','TGFB1','TNFA','PGE2','TRAIL')
for (cytokine in cytokines){
  print(cytokine)
  NK.seu[[paste0('cytosig_Z_',cytokine)]] <- t(cytosig_Z.all.res[cytokine,])
}

## 1.3 Calculate bubble Z and p values
bubble_z_p.res <- list()
for (cytokine in cytokines){
  print(cytokine)
  bubble_z_p.res[[cytokine]] <- bubbleZP(NK.seu@meta.data, celltype_col = 'FirstClust', value_col = paste0('cytosig_Z_',cytokine), lower = -3, upper = 3)
}

bubble.dat <- data.frame(bubble_z_p.res[['IL15']], cytokine = 'IL15')
for (cytokine in cytokines) {
  bubble.dat <- rbind(bubble.dat, data.frame(bubble_z_p.res[[cytokine]], cytokine = cytokine))
}
bubble.dat$cytokine <- factor(bubble.dat$cytokine, levels = cytokines)

## 1.4 Bubble plots
# ggballoonplot(bubble.dat, x = 'cytokine', y = 'celltype', color = 'Z', fill = 'Z', size = 'logP', ggtheme = theme_pubr()) +
#   scale_color_distiller(palette= "Spectral", limits = c(-max(abs(bubble.dat$Z)),max(abs(bubble.dat$Z)))) +
#   scale_fill_distiller(palette= "Spectral", limits = c(-max(abs(bubble.dat$Z)),max(abs(bubble.dat$Z)))) +
#   theme(legend.position='right',
#         strip.text = element_text(size = 12),
#         strip.background = element_blank(),
#         panel.border = element_rect(colour = "black", fill = NA))


ggballoonplot(bubble.dat, x = 'cytokine', y = 'celltype', color = 'z', fill = 'z', size = 'logP', ggtheme = theme_pubr()) +
  scale_color_gradientn(colors = colorRampPalette(colors = c('dodgerblue4','white','firebrick'))(12), limits = c(-max(abs(bubble.dat$z)),max(abs(bubble.dat$z)))) +
  scale_fill_gradientn(colors = colorRampPalette(colors = c('dodgerblue4','white','firebrick'))(12), limits = c(-max(abs(bubble.dat$z)),max(abs(bubble.dat$z)))) +
  theme(legend.position='right',
        strip.text = element_text(size = 12),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA)) 


ggballoonplot(bubble.dat, x = 'cytokine', y = 'celltype', color = 'Z', fill = 'Z', size = 'logP', ggtheme = theme_pubr()) +
  scale_color_gradientn(colors = colorRampPalette(colors = c('dodgerblue4','white','firebrick'))(12), limits = c(0,1)) +
  scale_fill_gradientn(colors = colorRampPalette(colors = c('dodgerblue4','white','firebrick'))(12), limits = c(0,1)) +
  theme(legend.position='right',
        strip.text = element_text(size = 12),
        strip.background = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA)) 

