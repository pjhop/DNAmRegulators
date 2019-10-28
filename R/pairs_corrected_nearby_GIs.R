## Collect baconed test-statistics (corrected for nearby GIs)

# Libraries
options("stringsAsFactors" = FALSE)
library(GenomicRanges)
library(tidyverse)
library(stringr)
library(BiocParallel)
library(BatchJobs)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(FDb.InfiniumMethylation.hg19)

# Directories / filepaths
direc <- "" # Working directory

# Load uncorrected results
pairs <- read_tsv(sprintf("%soutput.genome.wide/ewas/bonferroni.pairs.10.Mb.txt", direc))
pairs <- pairs %>% mutate(pair_id = paste(grs, probe, sep = "_"))
genes <- unique(pairs$grs)

# All test-statistics corrected for nearby GIs
files <- list.files(sprintf("%soutput.genome.wide/ewas.corrected.nearby.gi/output", direc), full.names = TRUE, pattern = "bacon")
files2 <- list.files(sprintf("%soutput.genome.wide/ewas.corrected.nearby.gi/output", direc), full.names = FALSE, pattern = "bacon")
files2 <- gsub("\\.bacon\\.txt\\.gz", "", files2)
# A few of the genes that are present in the uncorrected results are not present in the corrected results
# This is because these were excluded because of high correlation(s) with neighboring GI(s)
genes <- genes[genes %in% files2]

# For every pair, get the corrected p-value
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

# Run in parallel
conffile <- "config3.R"
BPPARAM <- BatchJobsParam(workers = 10, progressbar = FALSE, conffile = conffile)
register(BPPARAM)

## Corrected pairs
pairs.corrected  <- bplapply(genes, FUN = run, files = files, pairs = pairs)
pairs.corrected  <- do.call(rbind, pairs.corrected)

# Add ENSEMBL gene annotation
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

# Add information about genes
names(ensembl) <- c("EnsemblGeneID", "gene.grs", "start.grs",  "end.grs", "chr.grs")
pairs.corrected <- merge(pairs.corrected, ensembl, by.x = "grs", by.y = "EnsemblGeneID")

# Add information on cpg-sites
data(IlluminaHumanMethylation450kanno.ilmn12.hg19)
hm450 <- data.frame(getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19))
hm450 <- hm450[,c("Name", "pos", "chr","UCSC_RefGene_Name")]
colnames(hm450) <- c("cpg.probe", "cpg.pos", "cpg.chr", "cpg.gene")
pairs.corrected <- merge(pairs.corrected, hm450, by.x = "probe", by.y = "cpg.probe")

# Save corrected pairs
write_tsv(pairs.corrected, path = sprintf("%soutput.genome.wide/ewas.corrected.nearby.gi/pairs.10.Mb.corrected.txt", direc),
	       col_names = TRUE)
