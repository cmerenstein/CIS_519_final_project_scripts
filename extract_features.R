library(Seurat)
library(tidyverse) # library(dplyr) might be enough
library(matrixStats)

# input arguments
data.dir <- "/path/to/data_directory" # e.g. "CART4/filtered_feature_bc_matrix/"
out.filepath <- "path/to/output.txt" # e.g. "2021-04-13_CART4_features.txt"

# read in data
data <- Read10X(data.dir)
sobj <-  CreateSeuratObject(counts = data$`Gene Expression`)
sobj[['Protein']] <- CreateAssayObject(counts = data$`Antibody Capture`)

# initialize output matrix
out.df <- sobj@meta.data %>% select(-orig.ident)

# remove sparsity to do calculations
protein_matrix <- as.matrix(sobj@assays$Protein@counts)
RNA_matrix <- as.matrix(sobj@assays$RNA@counts)

# store number of protein and gene features
nProtein <- dim(protein_matrix)[1]
nRNA <- dim(RNA_matrix)[1]

# extract features
out.df <- out.df %>% mutate(
  percFeature_Protein = nFeature_Protein / nProtein,
  average_Protein = colMeans(protein_matrix),
  median_Protein = colMedians(protein_matrix),
  variance_Protein= colVars(protein_matrix),
  percFeature_RNA= nFeature_Protein / nRNA,
  average_RNA= colMeans(RNA_matrix)
)

# write output
out.df %>% rownames_to_column("cell_ID") %>%
  write.table(out.filepath, sep="\t", quote=F, row.names=F)
