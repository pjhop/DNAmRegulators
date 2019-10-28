# Evaluate predictive ability of genetic instrument over covariates
# Predictive ability reflected by F-statistic (ANOVA)

## Libraries
options("stringsAsFactors" = FALSE)
library(BiocParallel)
library(matrixStats)
library(edgeR)
library(limma)
library(GenomicRanges)
library(MASS)
library(broom)
library(BatchJobs)

# Directories / filepaths
direc <- "" # Working directory
rnaseq_path <- "" # .RData file

###################################
## covariates
###################################

# covariates train and test set
cvrt.test <- read.table(sprintf("%scvrt.test.txt", direc), header = TRUE, sep = "\t")
cvrt.train <- read.table(sprintf("%scvrt.train.txt", direc), header = TRUE, sep = "\t")
cvrt <- rbind(cvrt.train, cvrt.test)

###################################
## RNA-seq data, gene annotation
###################################

# ensembl biomart annotation
ensembl <- read.table(sprintf("%s/annotations/ensembl.txt", direc), header = TRUE, sep = "\t")
ensembl <- subset(ensembl,
    (gene_biotype == "protein_coding" | gene_biotype == "lincRNA") &
    !duplicated(ensembl_gene_id) & # remove duplicated
    chromosome_name %in% 1:22 # autosomes only
)

# RNA-seq data
load(rnaseq_path)

# match sample IDs (RNA-seq IDs)
cvrt <- cvrt[match(colnames(M), cvrt$uuid), ]

###################################
## annotation in GRanges format
###################################

# gene locations
probes.locs <- GRanges(
    seqnames = paste0("chr", ensembl$chromosome_name),
    ranges = IRanges(
    start = ensembl$start_position,
    end = ensembl$end_position
    ),
    gene = ensembl$ensembl_gene_id
)
probes.locs <- probes.locs[probes.locs$gene %in% rownames(M)]


###################################
## run
###################################

# files
genes <- list.files(sprintf("%soutput.genome.wide/grs.5.pcs.1se_rnaseq", direc), pattern = "Rdata")
genes <- gsub("\\.Rdata", "", genes)
files <- list.files(sprintf("%soutput.genome.wide/grs.5.pcs.1se_rnaseq", direc), pattern = "Rdata", full.names = TRUE)


# function to run per gene
run <- function(current.gene, files, M, cvrt.test, probes.locs, genes) {
    require(broom)
    require(GenomicRanges)

    # function to calculate partial R-squared
    partial_R2 <- function(fit.full, fit.reduced) {
            dev.full <- glance(fit.full)$deviance
            dev.reduced <- glance(fit.reduced)$deviance
            part_r2 <- (dev.reduced - dev.full) / dev.reduced
            part_r2
    }

    # get genetic instrument for current gene
    file <- grep(current.gene, files, value = TRUE)
    load(file)
    grs <- instrument$grs[cvrt.test$imputation_id]
    # observed expression for current gene
    y <- t(M[current.gene, cvrt.test$uuid, drop = FALSE])

    ## Reduced + Full model

    fit0 <- lm(y ~ # fit with covariates only
        Sampling_Age +
        Sex +
        Neutrophils +
        Lymphocytes +
        Monocytes +
        Eosinophils +
        biobank_id +
        PC1 + PC2 + PC3 + PC4 + PC5 , data = cvrt.test)

    fit1 <- lm(y ~ grs + # fit also includes the genetic instrument currently under investigation
        Sampling_Age +
        Sex +
        Neutrophils +
        Lymphocytes +
        Monocytes +
        Eosinophils +
        biobank_id  +
        PC1 + PC2 + PC3 + PC4 + PC5 , data = cvrt.test)

    # F-statistic
    f.stat <- anova(fit0, fit1)[2, "F"]

     # Partial R-squared
    part.R2 <- partial_R2(fit1, fit0)

    return(data.frame(
        gene = current.gene,
        f.stat = f.stat,
        part.R2 = part.R2
    ))
}

# Run for all genes
f.stats <- lapply(genes, FUN = run, files = files, M = M, cvrt.test = cvrt.test, probes.locs = probes.locs, genes=genes) # run
f.stats <- do.call(rbind, f.stats) # aggregate

# Save f-statistics
write.table(f.stats, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE,
    file = sprintf("%soutput.genome.wide/f.stats.txt", direc))
