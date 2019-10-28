# Get SNPs that are associated with white blood cell composition
library(tidyverse)

# Directories / filepaths
direc <- ""
genotype_path <- "" # directory where genotypes are stored per chromosome

# SNPLocations file
load(sprintf("%sdata/SNPLocations_cellcounts.RData", direc))

# Get SNPs per chromosome
chrs <- paste0("chr", 1:22, sep = "")

# Get SNPs
genotypes <- lapply(chrs, FUN = function(chr) {
	print(chr)
	load(sprintf("%s/%s_R20.3_hwe0.0001.RData", genotype_path, chr))
	genotypes <- genotypes[rownames(genotypes) %in% SNPLocations_cellcounts$POS_ID,]
	genotypes
})
genotypes <- do.call(rbind, genotypes)

# Save
save(genotypes, file = sprintf("%sdata/genotypes_cellcounts.RData", direc))
