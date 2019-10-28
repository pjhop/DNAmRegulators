## Collect baconed test-statistics (corrected for nearby GIs + nearby SNPs known to influence white blood cell composition)

# Libraries
options('stringsAsFactors' = FALSE)
library(GenomicRanges)
library(tidyverse)
library(stringr)
library(BiocParallel)
library(BatchJobs)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(FDb.InfiniumMethylation.hg19)

# Directories / filepaths
direc <- "" # Working directory

# Load pairs
pairs <- read_tsv(sprintf("%soutput.genome.wide/ewas.corrected.nearby.gi/pairs.corrected.max.fval.5.txt", direc))
genes <- unique(pairs$grs)

# Files and tested genes
files <- list.files(sprintf("%soutput.genome.wide/ewas.corrected.nearby.gi.snps.cellcounts/output", direc),
  full.names = TRUE, pattern = "bacon")
genes <- list.files(sprintf("%soutput.genome.wide/ewas.corrected.nearby.gi.snps.cellcounts/output", direc), pattern = "bacon")
genes <- gsub(".bacon.txt.gz", "", genes)
genes <- genes[genes %in% pairs$grs]

# Get results corrected for cellcount SNPs
run <- function(gene, files, pairs) {
    require(readr)
    require(dplyr)
    require(magrittr)

    file <- files[grep(gene, files)]
    probes <- pairs %>% filter(grs == gene) %$% probe
    me <- read_tsv(file, col_names = TRUE, progress = FALSE)
    me <- me %>% filter(probe %in% probes)
    me
}

# Set up parallelization
conffile <- "config3.R"
BPPARAM <- BatchJobsParam(workers = 10, progressbar = FALSE, conffile = conffile)
register(BPPARAM)

# Results corrected for cellcount SNPs
pairs.snps.cellcounts  <- bplapply(genes, FUN = run, files = files, pairs = pairs)
pairs.snps.cellcounts  <- do.call(rbind, pairs.snps.cellcounts)

# gene annotation
ensembl <- read.table(sprintf("%s/annotations/ensembl.txt", direc), header = TRUE, sep = "\t")
ensembl <- subset(ensembl, !duplicated(ensembl_gene_id))
ensembl <- ensembl[, c("ensembl_gene_id", "external_gene_name", "start_position", "end_position", "chromosome_name")]

gene.locs <- GRanges(
	seqnames = paste0("chr", ensembl$chromosome_name),
	ranges = IRanges(
		start = ensembl$start_position,
		end = ensembl$end_position
	),
	gene = ensembl$ensembl_gene_id
)

# add information about causal genes
names(ensembl) <- c("EnsemblGeneID", "gene.grs", "start.grs",  "end.grs", "chr.grs")
pairs.snps.cellcounts <- merge(pairs.snps.cellcounts, ensembl, by.x = "grs", by.y = "EnsemblGeneID")

# add information on cpg-sites
hm450 <- data.frame(getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19))
hm450 <- hm450[,c("Name", "pos", "chr","UCSC_RefGene_Name")]
colnames(hm450) <- c("cpg.probe", "cpg.pos", "cpg.chr", "cpg.gene")
pairs.snps.cellcounts <- merge(pairs.snps.cellcounts, hm450, by.x = "probe", by.y = "cpg.probe")

# Save
write_tsv(pairs.snps.cellcounts, path = sprintf("%soutput.genome.wide/ewas.corrected.nearby.gi.snps.cellcounts/bonferroni.pairs.10.Mb.snps.cellcounts.txt", direc),
	       col_names = TRUE)
