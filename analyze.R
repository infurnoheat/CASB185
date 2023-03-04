#install Sleuth and other dependencies
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()
BiocManager::install("devtools", force = TRUE)
BiocManager::install("pachterlab/sleuth", force = TRUE)
BiocManager::install("rhdf5", force = TRUE)
BiocManager::install("Rhdf5lib", force = TRUE)
BiocManager::install("biomaRt", force = TRUE)

#Loading Libraries
library(magrittr)
library(dplyr)
library(lazyeval)
library(ggplot2)
library(sleuth)
library(biomaRt)
library(rhdf5)
library(Rhdf5lib)

#Set WD
setwd("C:/Users/steve/CodeSpace/CASB185")

#Ignoring problematic H5 Files


#Setting directory to Kallisto Output Folder
base_dir <- "."
sample_id <- dir(file.path(base_dir, "kallisto_output"))
sample_id

#Setting directory to COVID data
kal_dirs <- sapply(sample_id, function(id) file.path(base_dir, "kallisto_output", id))
#kal_dirs_tsv <- paste0(kal_dirs, "/abundance.tsv")

list.files(kal_dirs)

#Loading Design Matrix
s2c <- read.table("C:/Users/steve/CodeSpace/CASB185/design.txt", header = TRUE, stringsAsFactors = FALSE)
s2c <- dplyr::mutate(s2c, path = kal_dirs)
s2c

#Loading Ensembl
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = 'ensembl.org')

t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name"), mart = mart)
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,ens_gene = ensembl_gene_id, ext_gene = external_gene_name)

#Running Sleuth
test1 <- H5Fopen("./kallisto_outputOLD/GSM4432378/abundance.h5")
test2 <- H5Fopen("C:/Users/steve/Downloads/abundance.h5")
test3 <- H5Fopen("C:/Users/steve/Downloads/moc_1_1_output/abundance.h5")
so <- sleuth_prep(s2c, ~ condition, target_mapping = t2g, num_cores = 1)
so <- sleuth_fit(so)
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')

sleuth_live(so)
