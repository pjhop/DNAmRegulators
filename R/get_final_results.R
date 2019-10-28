## Check sensitivity analyses:
# (1) adjusting for nearby SNPs known to influence cell type composition
# (2) adjusting for distant significant GIs associated with same CpG sites
## Remove associations that are insignificant in these analyses and save final results

# Libraries
options("stringsAsFactors" = FALSE)
library(tidyverse)
library(GenomicRanges)
library(magrittr)

# Directories / filepaths
direc <- "" # Working directory

# Load pairs
pairs <- read_tsv(sprintf("%soutput.genome.wide/ewas.corrected.nearby.gi/pairs.corrected.max.fval.5.txt", direc))

### Part1: Correction GIs on same chromosome that associate with same CpG(s) ###

# Load pairs of genes that influence same CpG(s)
gene_pairs <- read_tsv(sprintf("%sresultsv1/sensitivity/gene_pairs_distantgi.txt", direc))
genes <- unique(c(gene_pairs$gene1, gene_pairs$gene2))
pairs <- pairs %>% filter(grs %in% genes)

# Get test statistics
files <- list.files(sprintf("%soutput.genome.wide/ewas.corrected.distant.gi/output", direc), full.names = TRUE, pattern = "bacon")
run <- function(gene) {
    file <- files[grep(gene, files)]
    probes <- pairs %>% filter(grs == gene) %$% probe
    me <- read_tsv(file, col_names = TRUE, progress = FALSE)
    me <- me %>% filter(probe %in% probes)
    me
}

results <- lapply(genes, FUN = run)
results <- do.call(rbind, results)
results <- results %>% mutate(pair_id = paste(grs, probe, sep = "_"))

# Add corrected statistics
pairs <- results %>% rename(pvalue = "pvalue.corrected",
							tvalue = "tvalue.corrected") %>%
         select(pair_id, pvalue.corrected, tvalue.corrected) %>% right_join(pairs, by = "pair_id")

# Get number of genes tested to calculate alpha cutoff
min.F <- 10
F.values <- read.table(sprintf("%soutput.genome.wide/f.stats.txt", direc), header = TRUE, sep = "\t")
F.values.sig <- F.values %>% filter(f.stat > 10)

# Alpha cutoff
cpgs.tested <- read_lines(sprintf("%s/cpgs.tested.txt", direc))
alpha <- 0.05 / (length(cpgs.tested)  * length(files))

# Pairs that are significant after correction
pairs.sig <- pairs %>% filter(pvalue.corrected < alpha)

# Check pairs that are not significant (anymore)
pairs.notsig <- pairs %>% filter(!(pvalue.corrected < alpha))


# Genes that are predictive of other expression of a neighbouring gene associated with the same CpG
gene_pairs <- gene_pairs %>% filter(fval1 > 5 | fval2 > 5)
genes_fval_both <- gene_pairs %>% filter(fval1 > 5 & fval2 > 5)
genes_fval_both <- unique(c(genes_fval_both$gene1, genes_fval_both$gene2))
genes_fval_one <- gene_pairs %>% filter(!(fval1 > 5 & fval2 > 5))
genes_fval_one <- ifelse(genes_fval_one$fval1 > 5, genes_fval_one$gene1,
                     ifelse(genes_fval_one$fval2 > 5, genes_fval_one$gene2, NA))
genes_fval <- c(genes_fval_both, genes_fval_one)
pairs.sig <- pairs.sig %>% filter(!(grs %in% genes_fval))

# Pairs that should be removed
pairs.del.distant.gi <- pairs %>% filter(!(pair_id %in% pairs.sig$pair_id))
write_tsv(pairs.del.distant.gi, path = sprintf("%sresultsv1/sensitivity/pairs.deleted.distant.gi.txt", direc))

### Part2: Correction for SNPs that influence white blood cell composition ###

pairs <- read_tsv(sprintf("%soutput.genome.wide/ewas.corrected.nearby.gi/pairs.corrected.max.fval.5.txt", direc))

# Corrected pairs
pairs.corrected.snps.cellcounts <- read_tsv(sprintf("%soutput.genome.wide/ewas.corrected.nearby.gi.snps.cellcounts/bonferroni.pairs.10.Mb.snps.cellcounts.txt", direc))
pairs.corrected.snps.cellcounts <- mutate(pairs.corrected.snps.cellcounts, pair_id = paste(grs, probe, sep = "_"))
pairs.corrected.snps.cellcounts.sig <- pairs.corrected.snps.cellcounts %>% filter(pvalue < alpha)

# Pairs that should be removed (this includes 1 gene that was removed because it was strongly correlated with a cellcount-snp)
pairs.del.snps.cellcounts <- pairs %>% filter(!(pair_id %in% pairs.corrected.snps.cellcounts.sig$pair_id))
write_tsv(pairs.del.snps.cellcounts, path = sprintf("%sresultsv1/sensitivity/pairs.deleted.snps.cellcounts.txt", direc))

# Make plot of test statistics with and without correction

### Part3: Final dataset ###

pairs.final <- pairs %>% filter(!(pair_id %in% pairs.del.snps.cellcounts$pair_id), !(pair_id %in% pairs.del.distant.gi$pair_id))

# Write results
write_tsv(pairs.final, path = sprintf("%sresultsv1/pairs.final.txt", direc))
