


PBMC= readRDS("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/VF1/PBMC_V2Chem_VF1.rds")

DimPlot(PBMC, cols = palette)

NK3 = subset(PBMC, idents= "NK3")

p1 = VlnPlot(NK3, features= "KLRC2", group.by = "orig.ident")

VlnPlot(NK3, features= "KLRC2", group.by = "orig.ident")
VlnPlot(NK3, features= "IL32", group.by = "orig.ident")
VlnPlot(NK3, features= "CD52", group.by = "orig.ident")
VlnPlot(NK3, features= "CCL5", group.by = "orig.ident")


p2 = VlnPlot(PBMC, features= "KLRC2", group.by = "orig.ident")

png(file="/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/Final_Figures/VlnPlot_KLRC2_Across_Patient.png", width = 30, height = 15,  units = "cm", res=600 )
p1
dev.off()




#Look at 

#Extract only Romagnani Dataset
PBMC_Dataset4 = subset(PBMC, subset  = Dataset== "Dataset4")

#Separate CMVpos an neg

PBMC_Dataset4$orig.ident = droplevels(PBMC_Dataset4$orig.ident)
table(PBMC_Dataset4$orig.ident)

PBMC_Dataset4$CMVstatus = PBMC_Dataset4$orig.ident
levels(PBMC_Dataset4$CMVstatus) = c("CMVneg", "CMVneg", "CMVpos", "CMVpos", "CMVpos")

table(PBMC_Dataset4$CMVstatus)



NK3_Roma = subset(PBMC_Dataset4, idents= "NK3")

table(NK3_Roma$CMVstatus)

table(NK3_Roma$CMVstatus, NK3_Roma$nkg2c)

data_propo = table(PBMC_Dataset4$FirstClust, PBMC_Dataset4$nkg2c , PBMC_Dataset4$CMVstatus)



VlnPlot(PBMC_Dataset4 , group.by = "FirstClust", features= "KLRC2", split.by = "CMVstatus")



# Convert the 3D table to a long format dataframe
data_long <- as.data.frame.table(data_propo, responseName = "Count")

# Split the dataframe based on CMVstatus
data_CMVneg <- subset(data_long, Var3 == "CMVneg")
data_CMVpos <- subset(data_long, Var3 == "CMVpos")

# You might want to rename the columns for clarity
colnames(data_CMVneg) <- c("FirstClust", "NKG2Cstatus", "CMVstatus", "Count")
colnames(data_CMVpos) <- c("FirstClust", "NKG2Cstatus", "CMVstatus", "Count")

# Viewing the dataframes
print("CMV Negative:")
print(data_CMVneg)

print("CMV Positive:")
print(data_CMVpos)

# Calculate the proportions
data_long$Proportion <- with(data_long, ave(Count, Var1, Var3, FUN = function(x) x / sum(x)))

# Filter the data for NKG2Cpos only
data_NKG2Cpos <- subset(data_long, Var2 == "NKG2Cpos")

# Creating the bar graph
p6 = ggplot(data_NKG2Cpos, aes(x = Var1, y = Proportion, fill = Var3)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(title = "Comparison of Proportion of NKG2Cpos in each FirstClust for CMVpos and CMVneg",
       x = "FirstClust", y = "Proportion", fill= "HCMV Status") +
  theme_minimal()

p6

png(file="/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/01_GlobalHeterogeneity/Analysis_V2Chem/Final_Figures/BarGraph_NKG2C_CMVNegvsPos.png", width = 18, height = 12,  units = "cm", res=600 )
p6
dev.off()





