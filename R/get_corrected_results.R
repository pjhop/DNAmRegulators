## Load test-statistics, corrected for nearby GI, and select significant GI-CpG associations (Bonferroni)
## Then test whether significant GIs are predictive of neighbouring significant genes (F > 5).
## If F>5 and the tested gene shares target CpGs with the neighbouring gene (at a gene-level threshold), the gene is excluded.

# Libraries
options('stringsAsFactors' = FALSE)
library(tidyverse)
library(stringr)
library(GenomicRanges)
library(magrittr)

# directories/filepaths
direc <- ""
rnaseq_path <- ""

## Load all pairs and F.values
min.F <- 10 # F-statistic from ANOVA should be at least this

# file with F-stats
F.values <- read_tsv(sprintf("%soutput.genome.wide/f.stats.txt", direc))
F.values.sig <- F.values %>% filter(f.stat > 10)
all.genes <- F.values$gene # All genes (regardless of F-statistic)

# Load pairs uncorrected for nearby GIs
pairs.uncorrected <- read_tsv(sprintf("%soutput.genome.wide/ewas/bonferroni.pairs.10.Mb.txt", direc))
pairs.uncorrected <- pairs.uncorrected %>% mutate(pair_id = paste(grs, probe, sep = "_"))

# Load pairs corrected for nearby GIs
pairs.corrected <- read_tsv(sprintf("%soutput.genome.wide/ewas.corrected.nearby.gi/pairs.10.Mb.corrected.txt", direc))
pairs.corrected <- pairs.corrected %>% mutate(pair_id = paste(grs, probe, sep = "_"))

# Genes not tested because high correlation with neighboring GIs
genes_lost <- pairs.uncorrected  %>% filter(!(grs %in% pairs.corrected$grs)) %$% unique(grs)

# Bonferroni threshold
cpgs.tested <- read_lines(sprintf("%s/cpgs.tested.txt", direc))
alpha <- 0.05 / (length(cpgs.tested) * nrow(F.values.sig))

# Pairs that are significant after correction for nearby GIs
pairs.corrected <- pairs.corrected %>% filter(pvalue < alpha)

# Phenotype files
cvrt.test <- read.table(sprintf("%scvrt.test.txt", direc), header = TRUE, sep = "\t")
cvrt.train <- read.table(sprintf("%scvrt.train.txt", direc), header = TRUE, sep = "\t")
cvrt <- rbind(cvrt.train, cvrt.test)

# Preprocessed RNAseq data
load(rnaseq_path)

# ensembl biomart annotation
ensembl <- read.table(sprintf("%s/annotations/ensembl.txt", direc), header = TRUE, sep = "\t")
ensembl <- subset(ensembl,
                  (gene_biotype == "protein_coding" | gene_biotype == "lincRNA") &
                    !duplicated(ensembl_gene_id) & # remove duplicated
                    chromosome_name %in% 1:22 # autosomes only
)

# gene locations
probes.locs <- GRanges(
  seqnames = paste0("chr", ensembl$chromosome_name),
  ranges = IRanges(
    start = ensembl$start_position,
    end = ensembl$end_position
  ),
  gene = ensembl$ensembl_gene_id
)
probes.locs <- probes.locs[probes.locs$gene %in% all.genes,]

# Rdata files of genetic instruments
rdata.files <- list.files(sprintf("%soutput.genome.wide/grs.5.pcs.1se_rnaseq/", direc), pattern = "Rdata", full.names = TRUE)

# get genetic instruments
grs <- lapply(all.genes, function(gene) {
  file <- grep(gene, rdata.files, value = TRUE)
  load(file)
  out <- cbind(instrument$grs)
  colnames(out) <- gene
  rownames(out) <- names(instrument$grs)
  return(out)
})

grs <- do.call(cbind, grs)
grs <- grs[match(cvrt$imputation_id, rownames(grs)), , drop = FALSE] # match sample IDs


## Function to calculate if a gene is predictive of the expression of neighbouring genes (as reflected by the F-statistic)
run_neighbours <- function(current.gene) {

  # get genetic instrument for current gene
  X <- grs[cvrt.test$imputation_id, current.gene]

  # observed expression for current gene
  y <- t(M[current.gene, cvrt.test$uuid, drop = FALSE])

  # get nearby genes
  corr.genes <- subsetByOverlaps(probes.locs, probes.locs[probes.locs$gene == current.gene], maxgap = 1e6)
  corr.genes <- setdiff(corr.genes$gene, current.gene) # Exclude current gene
  corr.genes <- corr.genes[corr.genes %in% all.genes]
  corr.genes_names <- corr.genes

  # Check probes associated with this gene at gene-level bonferroni threshold (0.05 / nr.cpgs.tested)
  pairs_ <- pairs.corrected.genelevel %>%
    dplyr::filter(grs == current.gene)
  # Check which other genes affect CpGs associated with the current gene
  pairs_ <- pairs.corrected.genelevel %>% dplyr::filter(probe %in% pairs_$probe)
  check_genes <- unique(pairs_$grs) # Consider those genes

  # if any genes are within 1Mb calculate added predictive power over nearby instruments and covariates
  if (length(corr.genes) > 0) {
    fstats <- unlist(lapply(corr.genes, function(gn) {

      # Reduced model: nearby instruments + covariates
      fit0 <- lm(M[gn,cvrt.test$uuid] ~  # fit with covariates and nearby genetic instruments
                   grs[cvrt.test$imputation_id, corr.genes] +
                   Sampling_Age +
                   Sex +
                   Neutrophils +
                   Lymphocytes +
                   Monocytes +
                   Eosinophils +
                   biobank_id + PC1 + PC2 + PC3 + PC4 + PC5, data = cvrt.test)

      # Full model: current GI + nearby instruments + covariates
      fit1 <- lm(M[gn,cvrt.test$uuid] ~ X + grs[cvrt.test$imputation_id, corr.genes] + # fit also includes genetic instrument currently under investigation
                   Sampling_Age +
                   Sex +
                   Neutrophils +
                   Lymphocytes +
                   Monocytes +
                   Eosinophils +
                   biobank_id + PC1 + PC2 + PC3 + PC4 + PC5, data = cvrt.test)
      f.stat <- anova(fit0, fit1)[2, "F"]
      f.stat
    }))

    names(fstats) <- corr.genes

    # Select genes that affect same CpG-sites (at gene-level threshold) or were not tested (because of high correlations)
    fstats <- fstats[names(fstats) %in% pairs_$grs | names(fstats) %in% genes_lost]

    # Get max F-value
    if (length(fstats) > 0) {
      max.f.stat <- max(fstats) # Check highest F-statistic
      return(max.f.stat)
    } else {
      return(NA)
    }

  } else { # if no nearby genes are present return NA
    return(NA)
  }}

genes <- unique(pairs.corrected$grs)

# Pairs at a gene-level threshold
pairs.corrected.genelevel <- read_tsv(sprintf("%soutput.genome.wide/ewas.corrected.nearby.gi/pairs.corrected.genelevelthreshold.txt", direc))
pairs.corrected.genelevel <- pairs.corrected.genelevel %>% dplyr::filter(grs %in% pairs.corrected$grs)

# Calculate F-values on neighboring genes
tested.genes <- c(genes, genes_lost)
f.vals.with.neighbours.corrected <- map_dbl(genes, run_neighbours) # Run above function for all genes
names(f.vals.with.neighbours.corrected) <- genes

# Add these F-values to pairs dataframe
pairs.corrected <- pairs.corrected %>% left_join(data.frame(grs = names(f.vals.with.neighbours.corrected),
                                                            f.stats.nearby.grs.neighbours = f.vals.with.neighbours.corrected), by = "grs")

write_tsv(pairs.corrected, path = sprintf("%soutput.genome.wide/ewas.corrected.nearby.gi/pairs.corrected.fval.txt", direc))


# Keep genes with maximal F-value of <5
pairs.corrected <- pairs.corrected %>% filter(is.na(f.stats.nearby.grs.neighbours) | f.stats.nearby.grs.neighbours < 5)

# Save
write_tsv(pairs.corrected, path = sprintf("%soutput.genome.wide/ewas.corrected.nearby.gi/pairs.corrected.max.fval.5.txt", direc),
          col_names = TRUE)
