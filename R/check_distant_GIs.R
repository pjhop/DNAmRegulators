## Check which significant GIs share a GI on the same chromosome affecting the same CpG-site(s)
## + check whether GIs are predictive (as reflected by the F-statistic) of the expression of those distant
## GIs affecting the same CpG(s)

# Libraries
options('stringsAsFactors' = FALSE)
library(tidyverse)
library(GenomicRanges)

# Directories / filepaths
direc <- "" # Working directory
rnaseq_path <- "" # .RData file

# Load pairs
pairs <- read_tsv(sprintf("%soutput.genome.wide/ewas.corrected.nearby.gi/pairs.corrected.max.fval.5.txt", direc))

# Get probes that associate with more than one index gene
pairs <- pairs %>%
  group_by(probe) %>%
  mutate(n = n()) %>%
  filter(n > 1) %>%
  arrange(desc(n)) %>%
  ungroup()

# Unique probes
probes <- unique(pairs$probe)

# Get ensembl annotation
ensembl <- read.table(sprintf("%s/annotations/ensembl.txt", direc), header = TRUE, sep = "\t")
ensembl <- subset(ensembl,
    (gene_biotype == "protein_coding" | gene_biotype == "lincRNA") &
    !duplicated(ensembl_gene_id) & # remove duplicated
    chromosome_name %in% 1:22 # autosomes only
)

probes.locs <- GRanges(
    seqnames = paste0("chr", ensembl$chromosome_name),
    ranges = IRanges(
    start = ensembl$start_position,
    end = ensembl$end_position
    ),
    gene = ensembl$ensembl_gene_id
)
names(probes.locs) <- probes.locs$gene

# For each probe check the distance(s) between its index genes
fun <- function(probe) {
	pairs.tmp <- pairs[pairs$probe == probe,]

  # If there are two index genes for a given probe
	if(nrow(pairs.tmp) == 2) {
		genes <- pairs.tmp$grs
		dist <- distance(probes.locs[probes.locs$gene == genes[1]], probes.locs[probes.locs$gene == genes[2]])
		df <- data.frame(gene1 = genes[1], gene2 = genes[2], distance = dist)
	} else {
		# Unique pairs:
		genes_df <- t(combn(pairs.tmp$grs,m = 2))
		probes.locs1 <- probes.locs[genes_df[,1]]
		probes.locs2 <- probes.locs[genes_df[,2]]
		dist <- distance(probes.locs1, probes.locs2)
		genes_df <- data.frame(gene1 = genes_df[,1], gene2 = genes_df[,2], distance = dist)
		genes_df
	}
}

# Get distances between genes that are associated with the same probe
gene_pairs <- map_df(probes, .f = fun)

# Only check genes that are on the same chromosome and >1Mb
gene_pairs <- gene_pairs %>%
  filter(!is.na(distance), distance > 1e6) # Only checked genes on same chromosome!
gene_pairs <- gene_pairs[!duplicated(t(apply(gene_pairs, 1, sort))),] # Remove doubles

# Pairs that will be checked
pairs_check <- pairs %>%
  filter(grs %in% unique(c(gene_pairs$gene1, gene_pairs$gene2)))

## Load all pairs and F.values
min.F <- 10 # F-statistic from ANOVA should be at least this

# file with F-stats
F.values <- read_tsv(sprintf("%soutput.genome.wide/f.stats.txt", direc))
F.values.sig <- F.values %>% filter(f.stat > 10)
all.genes <- F.values$gene

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

# Genetic instruments
rdata.files <- list.files(sprintf("%soutput.genome.wide/grs.5.pcs.1se_rnaseq", direc), pattern = "Rdata", full.names = TRUE)

# Get all genetic instruments
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
run_neighbours <- function(current.gene, target.gene) {

    # get genetic instrument for current gene
    X <- grs[cvrt.test$imputation_id, current.gene]

    # observed expression for current gene
    y <- t(M[current.gene, cvrt.test$uuid, drop = FALSE])

    # get nearby genes
    corr.genes <- subsetByOverlaps(probes.locs, probes.locs[probes.locs$gene == current.gene], maxgap = 1e6)
    corr.genes <- setdiff(corr.genes$gene, current.gene)
    corr.genes <- corr.genes[corr.genes %in% all.genes]
    corr.genes_names <- corr.genes

    #
    gns <- gene_pairs %>% filter(gene1 == current.gene | gene2 == current.gene)
    gns <- unique(c(gns$gene1, gns$gene2))
    gns <- setdiff(gns, current.gene)

    # calculate added predictive power
    if (length(corr.genes) > 0) {
                # Reduced model: nearby instruments + distant instruments + covariates
                fit0 <- lm(M[target.gene,cvrt.test$uuid] ~
                    grs[cvrt.test$imputation_id, corr.genes] + grs[cvrt.test$imputation_id, gns] +
                    Sampling_Age +
                    Sex +
                    Neutrophils +
                    Lymphocytes +
                    Monocytes +
                    Eosinophils +
                    biobank_id + PC1 + PC2 + PC3 + PC4 + PC5, data = cvrt.test)

                # Full model: current GI + nearby instruments + distant instruments + covariates
                fit1 <- lm(M[target.gene,cvrt.test$uuid] ~ X + grs[cvrt.test$imputation_id, corr.genes] +
                	grs[cvrt.test$imputation_id, gns]	+
                    Sampling_Age +
                    Sex +
                    Neutrophils +
                    Lymphocytes +
                    Monocytes +
                    Eosinophils +
                    biobank_id + PC1 + PC2 + PC3 + PC4 + PC5, data = cvrt.test)
                 f.stat <- anova(fit0, fit1)[2, "F"]
                 f.stat
            } else {
            	fit0 <- lm(M[target.gene,cvrt.test$uuid] ~  # fit with covariates and distant genetic instruments
                    grs[cvrt.test$imputation_id, gns] +
                    Sampling_Age +
                    Sex +
                    Neutrophils +
                    Lymphocytes +
                    Monocytes +
                    Eosinophils +
                    biobank_id + PC1 + PC2 + PC3 + PC4 + PC5, data = cvrt.test)

                # Full model: current GI + distant instruments + covariates
                fit1 <- lm(M[target.gene,cvrt.test$uuid] ~ X +
                	grs[cvrt.test$imputation_id, gns]	+
                    Sampling_Age +
                    Sex +
                    Neutrophils +
                    Lymphocytes +
                    Monocytes +
                    Eosinophils +
                    biobank_id + PC1 + PC2 + PC3 + PC4 + PC5, data = cvrt.test)

                 f.stat <- anova(fit0, fit1)[2, "F"]
                 f.stat
            }
   }

# Run above function for all pairs of genes
mp <- map2_dbl(gene_pairs$gene1, gene_pairs$gene2, .f = run_neighbours)
mp2 <- map2_dbl(gene_pairs$gene2, gene_pairs$gene1, .f = run_neighbours)
gene_pairs$fval1 <- mp
gene_pairs$fval2 <- mp2

# Add gene info
gene_pairs <- gene_pairs %>% left_join(ensembl[,c("ensembl_gene_id", "external_gene_name")], by = c("gene1" = "ensembl_gene_id"))
gene_pairs <- gene_pairs %>% left_join(ensembl[,c("ensembl_gene_id", "external_gene_name")], by = c("gene2" = "ensembl_gene_id"), suffix = c(".gene1", ".gene2"))

# Save
write_tsv(gene_pairs, sprintf("%sresultsv1/sensitivity/gene_pairs_distantgi.txt", direc))
