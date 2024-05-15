#!/bin/bash


#singularity run -B /mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK3:/project /mnt/DOSI/PLATEFORMES/BIOINFORMATIQUE/03_WORKSPACE/rebuffet/docker_files/SCENIC/aertslab-pyscenic-0.12.0.sif \
#	pyscenic grn \
#       --num_workers 20 \
#        -o /project/05_Output/02_RegulatoryNetworkAnalysis/AllNK/expr_mat.adjacencies.tsv \
#        /project/05_Output/02_RegulatoryNetworkAnalysis/data/AllNK/AllNK.loom \
#        /project/01_Reference/TFList/allTFs_hg38.txt

#Replacing First step with Arboreto to prevent dask failure

singularity run -B /mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5:/project /mnt/DOSI/PLATEFORMES/BIOINFORMATIQUE/03_WORKSPACE/rebuffet/docker_files/SCENIC/aertslab-pyscenic-0.12.0.sif \
	arboreto_with_multiprocessing.py \
    	/project/05_Output/02_RegulatoryNetworkAnalysis/data/AllNK/AllNK.loom \
    	/project/01_Reference/TFList/allTFs_hg38.txt \
    	--method grnboost2 \
    	--output /project/05_Output/02_RegulatoryNetworkAnalysis/AllNK/expr_mat.adjacencies.tsv \
    	--num_workers 20 \
    	--seed 777

        
singularity run -B /mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5:/project /mnt/DOSI/PLATEFORMES/BIOINFORMATIQUE/03_WORKSPACE/rebuffet/docker_files/SCENIC/aertslab-pyscenic-0.12.0.sif \
    	pyscenic ctx \
        /project/05_Output/02_RegulatoryNetworkAnalysis/AllNK/expr_mat.adjacencies.tsv \
        /project/01_Reference/feather_database/hg19-tss-centered-5kb-7species.mc9nr.genes_vs_motifs.rankings.feather \
        /project/01_Reference/feather_database/hg19-tss-centered-10kb-7species.mc9nr.genes_vs_motifs.rankings.feather \
        --annotations_fname /project/01_Reference/motif/motifs-v9-nr.hgnc-m0.001-o0.0.tbl \
        --expression_mtx_fname /project/05_Output/02_RegulatoryNetworkAnalysis/data/AllNK/AllNK.loom \
        --mode "custom_multiprocessing" \
        --output /project/05_Output/02_RegulatoryNetworkAnalysis/AllNK/regulons.csv \
        --num_workers 20


singularity run -B /mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5:/project /mnt/DOSI/PLATEFORMES/BIOINFORMATIQUE/03_WORKSPACE/rebuffet/docker_files/SCENIC/aertslab-pyscenic-0.12.0.sif \
    pyscenic aucell \
        /project/05_Output/02_RegulatoryNetworkAnalysis/data/AllNK/AllNK.loom \
        /project/05_Output/02_RegulatoryNetworkAnalysis/AllNK/regulons.csv \
        -o /project/05_Output/02_RegulatoryNetworkAnalysis/AllNK/auc_mtx.loom \
        --num_workers 10
