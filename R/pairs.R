## Collect baconed test-statistics and select significant associations (Bonferroni)
## + exclude associations within <10Mb of the gene

# Libraries
options("stringsAsFactors" = FALSE)
library(BiocParallel)
library(GenomicRanges)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(FDb.InfiniumMethylation.hg19)
library(BatchJobs)

# Directories / filepaths
direc <- "" # Working directory

# Gene annotation
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

# get results (bacon-corrected), then returns significant pairs
run <- function(file, alpha) {
	require(readr)
    me <- read_tsv(file, col_names=TRUE)
    me <- subset(me, pvalue < alpha)
    return(me)
}

# Test-statistics for each gene (bacon-corrected)
files <- list.files(sprintf("%soutput.genome.wide/ewas/output", direc), pattern = "bacon", full.names = TRUE)

# Bonferroni cutoff
cpgs.tested <- read_lines(sprintf("%s/cpgs.tested.txt", direc))
alpha <- 0.05 / (length(cpgs.tested)  * length(files))

# Set up configuration for parallelization
conffile <- "config3.R"
BPPARAM <- BatchJobsParam(workers = 16, progressbar = FALSE, conffile = conffile)
register(BPPARAM)

# Run in parallel
pairs <- bplapply(files, FUN = run, alpha=alpha)
pairs <- do.call(rbind, pairs)

# add information about the genes
names(ensembl) <- c("EnsemblGeneID", "gene.grs", "start.grs",  "end.grs", "chr.grs")
pairs <- merge(pairs, ensembl, by.x = "grs", by.y = "EnsemblGeneID")

# add information on cpg-sites
hm450 <- data.frame(getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19))
hm450 <- hm450[,c("Name", "pos", "chr","UCSC_RefGene_Name")]
colnames(hm450) <- c("cpg.probe", "cpg.pos", "cpg.chr", "cpg.gene")
pairs <- merge(pairs, hm450, by.x = "probe", by.y = "cpg.probe")

# Save annotated results
write.table(pairs, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE,
    file = sprintf("%soutput.genome.wide/ewas/bonferroni.pairs.txt", direc))

gene.locs <- gene.locs[gene.locs$gene %in% pairs$grs,]


#################
## > 10 Mb
#################

# extract gene-CpG associations that are at least 10 Mb apart
fdb <- features(FDb.InfiniumMethylation.hg19)
fdb <- fdb[names(fdb) %in% pairs$probe,]

min.dist <- 10e6
overlaps <- findOverlaps(gene.locs, fdb, maxgap = min.dist)
overlaps <- data.frame(
	grs = gene.locs$gene[queryHits(overlaps)],
	probe = names(fdb)[subjectHits(overlaps)]
)
overlaps$pair <- with(overlaps, paste(grs, probe, sep = "_")) # assign each pair a unique ID

pairs <- subset(pairs, !paste(grs, probe, sep = "_") %in% overlaps$pair) # find those pairs that do not overlap

# Save annotated results (>10Mb)
write.table(pairs, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE,
    file = sprintf("%soutput.genome.wide/ewas/bonferroni.pairs.%s.Mb.txt", direc, min.dist/1e6))
