## Get results for each EWAS (corrected for nearby GIs) and inspect bias and inflation using bacon
## EWAS statistics after correction for inflation/bias using bacon are saved

# Libraries
options("stringsAsFactors" = FALSE)
library(BiocParallel)
library(bacon)
library(GenomicRanges)
library(BatchJobs)
library(readr)
library(dplyr)

# Working directory
direc <- ""

# Covariates
cvrt <- read.table(sprintf("%scvrt.meth.txt", direc), header = TRUE, sep = "\t")

# Get nr of covariates (these were saved in the ewas_corrected_nearby_GIs.R script)
genes <- gsub(".txt", "", list.files(sprintf("%soutput.genome.wide/ewas.corrected.nearby.gi/nr_covariates/", direc)))
nr.covariates <- unlist(lapply(genes, FUN = function(gene) {
    nr_covariate <- read.table(sprintf("%soutput.genome.wide/ewas.corrected.nearby.gi/nr_covariates/%s.txt", direc, gene))
    nr_covariate$V1
}))
names(nr.covariates) <- genes
save(nr.covariates, file = sprintf("%soutput.genome.wide/ewas.corrected.nearby.gi/nr.covariates.RData", direc))

# Function that loads results and performs bacon
# For each gene the test-statistics are saved in the same directory with ".bacon" appended
run <- function(gene, direc, nr.covariates, cvrt) {
    require(GenomicRanges)
    require(bacon)
    require(readr)
    require(dplyr)

    chrs <- paste0("chr",1:22)

    # Get results for current gene
    files <- lapply(chrs, FUN = function(chr) {
        file <- sprintf("%soutput.genome.wide/ewas.corrected.nearby.gi/output/%s.%s.txt", direc, chr, gene)
        me <- read_tsv(file, col_names=TRUE, progress=FALSE)
    })
    me <- bind_rows(files)

    set.seed(1) # reproducibility
    # run bacon
    bc <- bacon(effectsizes = cbind(me$estimate), standarderrors = cbind(me$se))

    me$estimate <- round(drop(es(bc)), digits = 6) # corrected regression coefficients
    me$se <- round(drop(se(bc)), digits = 6) # corrected standard errors
    me$tvalue <- drop(tstat(bc)) # corrected t-values
    me$pvalue <- 2*pt(abs(me$tvalue), nrow(cvrt) - nr.covariates[gene] - 1, lower.tail = FALSE) # corrected p-values

    write_tsv(me, gzfile(sprintf("%soutput.genome.wide/ewas.corrected.nearby.gi/output/%s.bacon.txt.gz", direc, gene)))

    # return bias and inflation for this EWAS
    out <- c("bias" = bias(bc), "inflation" = inflation(bc))
    return(out)
}

genes <- names(nr.covariates)

# Set up configuration for parallelization
conffile <- "config3.R"
BPPARAM <- BatchJobsParam(workers = 16, progressbar = FALSE, conffile = conffile)
register(BPPARAM)

# Run
res <- bplapply(genes, FUN = run, direc=direc, nr.covariates=nr.covariates, cvrt=cvrt)
res <- do.call(rbind, res)
res <- as.data.frame(cbind(genes, round(res, digits = 6)))
names(res) <- c("gene", "bias", "inflation")

# Write dataframe with bias and inflation per GI
write.table(res, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t",
    file = sprintf("%soutput.genome.wide/ewas.corrected.nearby.gi/bias_inflation_ewas.txt", direc))
