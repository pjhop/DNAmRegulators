## Extract locations of cellcount-associated SNPs from published data (Roederer et al. and Orru et al.)

## Libraries
library(VariantAnnotation)
library(BiocParallel)
library(readr)
library(stringr)
library(dplyr)
library(readxl)
options("stringsAsFactors" = FALSE)

# Directories / filepaths
direc <- "" # Working directory
snplocations_path <- "" # .RData file

# Load SNPlocations file
load(snplocations_path)

# http://www.cell.com/cell/abstract/S0092-8674(15)00247-0
snps.roederer <- read_excel(sprintf("%sextdata/snps.cellcounts.roederer.xlsx", direc), skip=3) # These SNPS should also be in spel and csf
snps.roederer.spel <- read_excel(sprintf("%sextdata/snps.cellcounts.roederer.spel.xlsx", direc), skip=2)
snps.roederer.csf <- read_excel(sprintf("%sextdata/snps.cellcounts.roederer.csf.xlsx", direc), skip=2)

# http://www.cell.com/cell/fulltext/S0092-8674(13)01072-6
snps_orru <- read_excel(sprintf("%sextdata/snps.cellcounts.orru.xlsx", direc), skip=3)
snps_orru_lead <- read_excel(sprintf("%sextdata/snps.cellcounts.orru2.xlsx", direc), skip=3, sheet = 3)
snps_orru_notvalidated <- read_excel(sprintf("%sextdata/snps.cellcounts.orru2.xlsx", direc), skip=4, sheet=2)

# Format tables
snps_orru$`SNP (chromosome: position)` <- str_replace(snps_orru$`SNP (chromosome: position)`, "chr", "")
snps_orru_notvalidated$`topSNP (chr:position) /rsID` <- str_split(snps_orru_notvalidated$`topSNP (chr:position) /rsID`, pattern = "/", simplify = TRUE)[,1]
snps_orru_notvalidated$`topSNP (chr:position) /rsID` <- str_replace(snps_orru_notvalidated$`topSNP (chr:position) /rsID`, "chr", "")
snps_orru_notvalidated$`topSNP (chr:position) /rsID` <- str_trim(snps_orru_notvalidated$`topSNP (chr:position) /rsID`)
snps_orru_lead$POS_ID <- paste(snps_orru_lead$Chr, snps_orru_lead$Pos, sep = ":")

# Select these SNPs from SNPLocations file
SNPLocations_cellcounts <- SNPLocations %>% dplyr::filter(RS_ID %in% snps.roederer$Marker | RS_ID %in% snps.roederer.spel$Marker |
	                                               RS_ID %in% snps.roederer.csf$Marker | POS_ID %in% snps_orru$`SNP (chromosome: position)` | POS_ID %in% snps_orru_lead$POS_ID |
	                                               POS_ID %in% snps_orru_notvalidated$`topSNP (chr:position) /rsID`)

# Save
save(SNPLocations_cellcounts, file = sprintf("%sdata/SNPLocations_cellcounts.RData", direc))
