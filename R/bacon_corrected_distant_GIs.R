## Get results for each EWAS (corrected for nearby + distant GIs) and inspect bias and inflation using bacon
## EWAS statistics after correction for inflation/bias using bacon are saved

# Libraries
options('stringsAsFactors' = FALSE)
library(BiocParallel)
library(bacon)
library(GenomicRanges)
library(BatchJobs)
library(readr)
library(dplyr)

# Working directory
direc <- ""

# Phenotype file
cvrt <- read.table(sprintf("%scvrt.meth.txt", direc), header = TRUE, sep = "\t")

# Get nr of covariates
genes <- gsub(".txt", "", list.files(sprintf("%soutput.genome.wide/ewas.corrected.distant.gi/nr_covariates/", direc)))
nr.covariates <- unlist(lapply(genes, FUN = function(gene) {
    nr_covariate <- read.table(sprintf("%soutput.genome.wide/ewas.corrected.distant.gi/nr_covariates/%s.txt", direc, gene))
    nr_covariate$V1
}))
names(nr.covariates) <- genes
save(nr.covariates, file = sprintf("%soutput.genome.wide/ewas.corrected.distant.gi/nr.covariates.RData", direc))

run <- function(gene, direc, nr.covariates, cvrt) {
    require(GenomicRanges)
    require(bacon)
    require(readr)
    require(dplyr)

    chrs <- paste0("chr",1:22)

    # Get results for current gene
    files <- lapply(chrs, FUN = function(chr) {
        file <- sprintf("%soutput.genome.wide/ewas.corrected.distant.gi/output/%s.%s.txt", direc, chr, gene)
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

    write_tsv(me, gzfile(sprintf("%soutput.genome.wide/ewas.corrected.distant.gi/output/%s.bacon.txt.gz", direc, gene)))

    # return bias and inflation for this EWAS
    out <- c("bias" = bias(bc), "inflation" = inflation(bc))
    return(out)
}

# Genes
genes <- list.files(sprintf("%soutput.genome.wide/ewas.corrected.distant.gi/output", direc))
genes <- gsub("chr[0-9]{1,}\\.|\\.txt", "", genes)
genes <- unique(genes)

# Run
res <- lapply(genes, FUN = run, direc=direc, nr.covariates=nr.covariates, cvrt=cvrt)
res <- do.call(rbind, res)
res <- as.data.frame(cbind(genes, round(res, digits = 6)))
names(res) <- c("gene", "bias", "inflation")

# Save bias and inflation statistics
write.table(res, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t",
    file = sprintf("%soutput.genome.wide/ewas.corrected.distant.gi/bias_inflation_ewas.txt", direc))
