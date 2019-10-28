## For each genetic instrument, perform an EWAS

# Libraries
options("stringsAsFactors" = FALSE)
library(BiocParallel)
library(edgeR)
library(limma)
library(DBI)
library(GenomicRanges)
library(matrixStats)
library(readr)
library(dplyr)

# Collect arguments
args <- commandArgs(trailingOnly=TRUE)
N <- as.numeric(args[1]) # Total number of partitions (defined in 'submit_ewas.R' script)
i <- as.numeric(args[2]) # This partition

# Directories / filepaths
direc <- "" # Working directory
mvalues_path <- "" # .RData file

###############################
## select genes to test
###############################

min.F <- 10 # F-statistic from ANOVA should be at least this
# file with F-stats
F.values <- read.table(sprintf("%soutput.genome.wide/f.stats.txt", direc), header = TRUE, sep = "\t")

# Select genes with F > 10
genes <- F.values$gene[F.values$f.stat >= min.F] # Select genes with F > 10

###################################
## data
###################################

# covariates train and test set
cvrt <- read.table(sprintf("%scvrt.meth.txt", direc), header = TRUE, sep = "\t")

# Inverse rank transformed Mvalues
load(mvalues_path)

###############################
## Genetic instruments
###############################

rdata.files <- list.files(sprintf("%soutput.genome.wide/grs.5.pcs.1se_meth", direc), pattern = "Rdata", full.names = TRUE)

# get genetic instruments that were calculated in construct_GIs.R
grs <- bplapply(genes, FUN = function(gene) {
    file <- grep(gene, rdata.files, value = TRUE)
    load(file)
    out <- cbind(instrument$grs)
    colnames(out) <- gene
    rownames(out) <- names(instrument$grs)
    return(out)
}, BPPARAM = MulticoreParam(workers = 2, verbose = TRUE))

grs <- do.call(cbind, grs)
grs <- grs[match(cvrt$imputation_id, rownames(grs)), , drop = FALSE] # match sample IDs


###############################
## Run EWAS
###############################

# Transpose M-values
Y <- t(mvalues)
rm(mvalues); gc()

# model matrix
Z <- model.matrix(~ Sampling_Age +
    Sex +
    Neutrophils +
    Lymphocytes +
    Monocytes +
    Eosinophils +
    biobank_id +
    PC1 + PC2 + PC3 + PC4 + PC5,
    data = cvrt)

k <- ncol(Z) - 1
n <- nrow(Y)

# fast way of performing many linear models
# method explained here: http://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-14-166
U1 <- crossprod(Z, Y)
U2 <- solve(crossprod(Z), U1)
ytr <- Y - Z %*% U2

run <- function(current.gene) {
    X <- grs[cvrt$imputation_id, current.gene, drop = FALSE]
    U3 <- crossprod(Z, X)
    U4 <- solve(crossprod(Z), U3)
    Xtr <- X - Z %*% U4
    Xtr2 <- colSums(Xtr**2)

    b <- crossprod(ytr, Xtr)
    b <- b / matrix(Xtr2, nr = nrow(b), nc = ncol(b), byrow = TRUE)
    term1 <- colSums(ytr^2)
    term2 <- matrix(Xtr2, nc = ncol(b), nr = nrow(b), byrow = TRUE) * (b**2)
    term1 <- drop(term1)
    term2 <- drop(term2)
    sig <- (term1 - term2) / (n-k-2)
    err <- sqrt(sig * (1 / Xtr2))
    me <- data.frame(
        grs = current.gene, # gene currently investigated
        probe = names(err), # target cpgs
        estimate = round(drop(b), digits = 6), # regression coefficient
        se = round(drop(err), digits = 6) # standard error
    )

    write_tsv(me, col_names = TRUE,
        path = gzfile(sprintf("%soutput.genome.wide/ewas/output/%s.txt.gz", direc, current.gene)))
    return(NULL)
}

# Make N partitions, submit this script with the "submit_ewas.R" script
# specificying how many partitions should be used
nr <- floor(length(genes)/N)
ends <- cumsum(rep(nr,N))
ends[N] = length(genes)
begins <- dplyr::lag(ends) + 1
begins[1] <- 1
genes <- genes[begins[i]:ends[i]]

# Run for this partition
null <- lapply(genes, FUN = run)
