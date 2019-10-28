## Construct Genetic Instruments (GIs)

## Libraries
options("stringsAsFactors" = FALSE)
library(BiocParallel)
library(matrixStats)
library(edgeR)
library(limma)
library(GenomicRanges)
library(MASS)
library(glmnet)

# Chromosome (run in parallel over all chromosomes)
args = commandArgs(trailingOnly=TRUE)
chr = args[1]

# Directories / filepaths
direc <- "" # Working directory
snplocations_path <- "" # .RData file
rnaseq_path <- "" # .RData file
genotype_path <- "" # directory where genotypes are stored per chromosome

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

# Locations of SNPs
load(snplocations_path) # annotation for SNPs, object: SNPLocations

# RNAseq data
load(rnaseq_path)

# match sample IDs (RNA-seq IDs)
cvrt <- cvrt[match(colnames(M), cvrt$uuid), ]

###################################
## annotation in GRanges format
###################################

# Used later on to find the nearest SNPs per gene to predict expression with

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

# SNP locations
snps.locs <- GRanges(
    seqnames = paste0("chr", SNPLocations$CHROM),
    ranges = IRanges(
        start = SNPLocations$POS,
        end = SNPLocations$POS
    ),
    snp = SNPLocations$POS_ID
)

###################################
## run
###################################

nfolds <- 5 # K-fold cross-validation
type.measure <- "mse" # metric used to optimize in cross-validation
alpha <- 1 # mixing parameter in elastic net in glmnet package
# 1 corresponds to lasso; 0 corresponds ridge

# Select genes on current chromosome and find overlap with SNPs (<100kb)
genes <- ensembl$ensembl_gene_id[ensembl$chromosome_name == chr] # genes located on current chromosome
genes.tmp <- probes.locs[probes.locs$gene %in% genes] # only use annotation for those genes
overlaps <- findOverlaps(genes.tmp, snps.locs, maxgap = 1e5) # find nearby SNPs, per gene (max 100 Kb)

# SNPs used per gene in data frame
overlaps <- data.frame(
    gene = genes.tmp$gene[queryHits(overlaps)],
    snp = snps.locs$snp[subjectHits(overlaps)]
)

# genotype data per chromosome
load(sprintf("%s/chr%s_R20.3_hwe0.0001.RData", genotype_path, chr))
# extract only SNPs that are needed, saves memory
genotypes <- genotypes[, colnames(genotypes) %in% cvrt$imputation_id & !duplicated(colnames(genotypes))]
gc()
# Calculate MAFs, and select SNPs with MAF > 0.01
MAFs <- rowMeans(genotypes) / 2
MAFs <- ifelse(MAFs > 0.5, 1 - MAFs, MAFs)
genotypes <- genotypes[MAFs > 0.01,]

# extract only SNPs that are needed, saves memory
genotypes <- genotypes[rownames(genotypes) %in% overlaps$snp, , drop = FALSE]

# genes located on current chromosome, but now only unique values
genes <- unique(overlaps$gene); gc()

# loop over genes
# per gene, fit lasso on train data with nearby SNPs
out <- lapply(genes, function(current.gene) {

    # snps to use in prediction
    snps <- overlaps$snp[overlaps$gene == current.gene]
    snps <- intersect(snps, rownames(genotypes)) # snps should also be in data

    if (length(snps) == 0) {
        write.table(current.gene, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t",
            file = sprintf("%soutput.genome.wide/grs.5.pcs.1se_rnaseq/no.snps.available.chr.%s.%s.txt", direc, chr, current.gene))
        write.table(current.gene, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t",
            file = sprintf("%soutput.genome.wide/grs.5.pcs.1se_meth/no.snps.available.chr.%s.%s.txt", direc, chr, current.gene))
    } else {

    ## train

    # temporary object for SNP data, only train data
    geno <- t(genotypes[rownames(genotypes) %in% snps, cvrt.train$imputation_id, drop = FALSE])

    # temporary object for gene expression levels, should be matrix, only train data
    y <- t(M[current.gene, cvrt.train$uuid, drop = FALSE])

    # first part of model matrix.
    mm <- model.matrix(~ 0 + biobank_id + Sampling_Age + Sex + Neutrophils + Lymphocytes + Monocytes + Eosinophils +
                       PC1 + PC2 + PC3 + PC4 + PC5, data = cvrt.train)

    # full model matrix
    x <- cbind(mm, geno)

    # penalty factors. Covariates are unpenalized, only SNPs are penalized
    penalty.factors <- c(rep(0, ncol(mm)), rep(1, ncol(x)-ncol(mm)))

    set.seed(10) # reproducibility
    # determine folds
    foldid <- sample(rep(1:5, length = nrow(cvrt.train)))
    # fit model using 5-fold cross-validation
    cv.fit <- cv.glmnet(x = x, y = y, foldid = foldid, type.measure = type.measure, penalty.factor = penalty.factors)

    # get regression coefficients in model and corresponding SNP names
    beta <- drop(as.matrix(coef(cv.fit, s = "lambda.min")))
    beta <- beta[beta != 0]
    snps <- intersect(names(beta), snps)

    # second LASSO step:
    if (length(snps) > 0) {
        train.snps <- geno[cvrt.train$imputation_id,snps,drop=FALSE]

        if(length(snps) == 1) {
            model <- lm(M[current.gene, cvrt.train$uuid] ~ train.snps + Sampling_Age + Sex + biobank_id + Neutrophils +
                                                           Eosinophils + Monocytes + Lymphocytes +
                                                           PC1 + PC2 + PC3 + PC4 + PC5, data = cvrt.train)
            beta <- coef(model)[2:(ncol(train.snps)+1)]
            names(beta) <- colnames(train.snps)
        } else {
          mm <- model.matrix(~ 0 + biobank_id + Sampling_Age + Sex + Neutrophils + Lymphocytes + Monocytes + Eosinophils +
                             PC1 + PC2 + PC3 + PC4 + PC5, data = cvrt.train)
          x <- cbind(mm, train.snps)

          # Gene expression
          y <-  M[current.gene, cvrt.train$uuid]

          # Age etc. are unpenalized
          penalty.factors <- c(rep(0, ncol(mm)), rep(1, ncol(x)-ncol(mm)))

          # Fit model using 5-fold cross-validation
          foldid <- sample(rep(1:5, length = nrow(cvrt.train)))
          cv.fit <- cv.glmnet(x = x, y = y, foldid = foldid, type.measure = "mse", penalty.factor = penalty.factors)

          # Function from glmnet https://github.com/cran/glmnet
          # Adjusted so that at least 1 SNP is nonzero
          getmin <- function(lambda,cvm,cvsd, nzero, nr.covariates = ncol(mm)){
                        cvmin <- min(cvm, na.rm = TRUE)
                        idmin <- cvm <= cvmin
                        lambda.min <- max(lambda[idmin], na.rm = TRUE)
                        idmin <- match(lambda.min, lambda)
                        semin <- (cvm + cvsd)[idmin]
                        idmin <- cvm <= semin
                        nzero.1se <- nzero[idmin]
                        lambda.1se <- lambda[idmin]
                        if(max(nzero.1se) == nr.covariates) {
                          lambda.1se <- max(lambda[idmin], na.rm = TRUE)
                          } else {
                              lambda.1se <- lambda.1se[nzero.1se != nr.covariates]
                              lambda.1se <- max(lambda.1se, na.rm = TRUE)
                          }
                        list(lambda.min = lambda.min,lambda.1se = lambda.1se)
          }

          # Get maximum lambda within 1se of the minimum mse, with the constraint that at least one
          # SNP has a nonzero coefficient
          min <- getmin(cv.fit$lambda, cv.fit$cvm, cv.fit$cvsd, nzero = cv.fit$nzero, nr.covariates = ncol(mm))

          # Get SNPs with nonzero coefficient
          beta <- drop(as.matrix(coef(cv.fit, s = min$lambda.1se)))

          # Select non-zero coefficients
          beta <- beta[beta != 0]
          beta <- beta[intersect(names(beta), colnames(train.snps))]
        }

          geno <- t(genotypes[names(beta), cvrt$imputation_id,drop=FALSE])
          grs <- drop(geno %*% cbind(beta))
          instrument <- list(
              grs = grs,
              regression.coefs = beta,
              snps = names(beta)
          )
          # Save the instrument
          save(instrument, file = sprintf("%soutput.genome.wide/grs.5.pcs.1se_rnaseq/%s.Rdata", direc, current.gene))

    } else {
        # if none of the SNPs has any predictive ability, then write a file that says so
        write.table(current.gene, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t",
            file = sprintf("%soutput.genome.wide/grs.5.pcs.1se_rnaseq/no.snps.available.after.lasso.chr.%s.%s.txt", direc, chr, current.gene))
        write.table(current.gene, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t",
            file = sprintf("%soutput.genome.wide/grs.5.pcs.1se_meth/no.snps.available.after.lasso.chr.%s.%s.txt", direc, chr, current.gene))
        return(NULL)
    }}
})
