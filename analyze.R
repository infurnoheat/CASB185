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

#Setting directory to Kallisto Output Folder
base_dir <- "."
sample_id <- dir(file.path(base_dir, "kallisto_data"))
sample_id

#Setting directory to COVID data
kal_dirs <- sapply(sample_id, function(id) file.path(base_dir, "kallisto_data", id))

#Loading Design Matrix
s2c <- read.table("C:/Users/steve/CodeSpace/CASB185/design.txt", header = TRUE, stringsAsFactors = FALSE)
s2c <- dplyr::mutate(s2c, path = kal_dirs)
s2c

#Loading Ensembl
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = 'ensembl.org')

t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name"), mart = mart)
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,ens_gene = ensembl_gene_id, ext_gene = external_gene_name)

#Data Checking
# Get a list of paths to the H5 files
kal_files <- file.path(s2c$path, "abundance.h5")
# Use h5read function from rhdf5 to get the number of targets
n_targets <- sapply(kal_files, function(file) {
  ids <- rhdf5::h5read(file, "aux/ids")
  n_targets <- length(ids)
  n_targets
})
# Get the kallisto version
versions <- sapply(kal_files, function(file) {
  file_name <- basename(file)
  version <- rhdf5::h5read(file, "aux/kallisto_version")
  version
})
# give each entry the name of your sample
names(n_targets) <- names(versions) <- basename(s2c$path)
# if they are consistent, the unique function should just return one value
unique(n_targets)
unique(versions)
# if they are not consistent, you can inspect each entry to find the culprit(s)
n_targets
versions

#Running Sleuth
so <- sleuth_prep(s2c, ~ condition, target_mapping = t2g, num_cores = 1, verbose)
so <- sleuth_fit(so)
so <- sleuth_fit(so, ~1, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')

sleuth_live(so)
