## Get results for each EWAS and inspect bias and inflation using bacon
## EWAS statistics after correction for inflation/bias using bacon are saved

# Libraries
options("stringsAsFactors" = FALSE)
library(BiocParallel)
library(bacon)
library(GenomicRanges)
library(BatchJobs)
library(readr)
library(dplyr)

# Directories / filepaths
direc <- "" # Working directory

# Covariates
cvrt <- read.table(sprintf("%scvrt.meth.txt", direc), header = TRUE, sep = "\t")

# nr.covariates
nr.covariates <- ncol(model.matrix(~ Sampling_Age +
    Sex +
    Neutrophils +
    Lymphocytes +
    Monocytes +
    Eosinophils +
    biobank_id +
    PC1 + PC2 + PC3 + PC4 + PC5,
    data = cvrt)
)

# Function that loads results and performs bacon
# For each gene the test-statistics are saved in the same directory with ".bacon" appended
run <- function(gene, direc, nr.covariates, cvrt) {
    require(GenomicRanges)
    require(bacon)
    require(readr)
    require(dplyr)

    # Filename
    file <- sprintf("%soutput.genome.wide/ewas/output/%s.txt.gz", direc, gene)

    # get results for this EWAS
    me <- read_tsv(file, col_names=TRUE)

    set.seed(1) # reproducibility

    # run bacon
    bc <- bacon(effectsizes = cbind(me$estimate), standarderrors = cbind(me$se))

    # get estimate, se, tvalue and pvalue
    me$estimate <- round(drop(es(bc)), digits = 6) # corrected regression coefficients
    me$se <- round(drop(se(bc)), digits = 6) # corrected standard errors
    me$tvalue <- drop(tstat(bc)) # corrected t-values
    me$pvalue <- 2*pt(abs(me$tvalue), nrow(cvrt) - nr.covariates - 1, lower.tail = FALSE) # corrected p-values

    write_tsv(me, gzfile(gsub("txt.gz", "bacon.txt.gz", file))) # same file, but just "bacon" added

    # return bias and inflation for this EWAS
    out <- c("bias" = bias(bc), "inflation" = inflation(bc))

    return(out)
}

# Genes (from ewas.R script)
genes <- list.files(sprintf("%soutput.genome.wide/ewas/output", direc), pattern = ".txt.gz")
genes <- gsub(".txt.gz", "", genes)

# Set up configuration for parallelization
conffile <- "config3.R"
BPPARAM <- BatchJobsParam(workers = 16, progressbar = FALSE, conffile = conffile)
register(BPPARAM)

# Run bacon for each gene
res <- bplapply(genes, FUN = run, direc = direc, nr.covariates = nr.covariates, cvrt = cvrt)
res <- do.call(rbind, res)

# Get bias and inflation for each gene
res <- as.data.frame(cbind(genes, round(res, digits = 4)))
names(res) <- c("gene", "bias", "inflation")

# Save bias and inflation
write.table(res, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t",
    file = sprintf("%soutput.genome.wide/ewas/bias_inflation_ewas.txt", direc))
