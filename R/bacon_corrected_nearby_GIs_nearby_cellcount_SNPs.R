## Get results for each EWAS (corrected for nearby GIs, and nearby SNPs that influence white blood cell composition)
## and inspect bias and inflation using bacon
## EWAS statistics after correction for inflation/bias using bacon are saved

options("stringsAsFactors" = FALSE)
library(BiocParallel)
library(bacon)
library(GenomicRanges)
library(BatchJobs)
library(readr)
library(dplyr)

# Working directory
direc <- "" # Working directory

# Phenotype file
cvrt <- read.table(sprintf("%scvrt.meth.txt", direc), header = TRUE, sep = "\t")

# Get nr of covariates (these were saved in the ewas_corrected_nearby_GIs_nearby_cellcount_SNPs.R script)
genes <- gsub(".txt", "", list.files(sprintf("%soutput.genome.wide/ewas.corrected.nearby.gi.snps.cellcounts/nr_covariates/", direc)))
nr.covariates <- unlist(lapply(genes, FUN = function(gene) {
    nr_covariate <- read.table(sprintf("%soutput.genome.wide/ewas.corrected.nearby.gi.snps.cellcounts/nr_covariates/%s.txt", direc, gene))
    nr_covariate$V1
}))
names(nr.covariates) <- genes
save(nr.covariates, file = sprintf("%soutput.genome.wide/ewas.corrected.nearby.gi.snps.cellcounts/nr.covariates.RData", direc))

run <- function(gene, direc, nr.covariates, cvrt) {
    require(GenomicRanges)
    require(bacon)
    require(readr)
    require(dplyr)

    chrs <- paste0("chr",1:22)

    # Get results for current gene
    files <- lapply(chrs, FUN = function(chr) {
        file <- sprintf("%soutput.genome.wide/ewas.corrected.nearby.gi.snps.cellcounts/output/%s.%s.txt", direc, chr, gene)
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

    write_tsv(me, gzfile(sprintf("%soutput.genome.wide/ewas.corrected.nearby.gi.snps.cellcounts/output/%s.bacon.txt.gz", direc, gene)))

    # return bias and inflation for this EWAS
    out <- c("bias" = bias(bc), "inflation" = inflation(bc))
    return(out)
}

# Set up configuration for parallelization
conffile <- "config3.R"
BPPARAM <- BatchJobsParam(workers = 16, progressbar = FALSE, conffile = conffile)
register(BPPARAM)

# Genes
genes <- list.files(sprintf("%soutput.genome.wide/ewas.corrected.nearby.gi.snps.cellcounts/output", direc))
genes <- gsub("chr[0-9]{1,}\\.|\\.txt", "", genes)
genes <- unique(genes)

# Run
res <- bplapply(genes, FUN = run, direc=direc, nr.covariates=nr.covariates, cvrt=cvrt)
res <- do.call(rbind, res)
res <- as.data.frame(cbind(genes, round(res, digits = 6)))
names(res) <- c("gene", "bias", "inflation")

# Save bias and inflation test-statistics
write.table(res, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t",
    file = sprintf("%soutput.genome.wide/ewas.corrected.nearby.gi.snps.cellcounts/bias_inflation_ewas.txt", direc))
