## Perform EWAS per genetic instrument (GI), while correcting for distant significant GIs (same chromosome)
# Note that this script runs per chromosome (chr.. should be supplied as argument to this script)
# The submit_ewas_corrected_distant_GIs.R script was used to run this script for each chromosome

options('stringsAsFactors' = FALSE)
library(BiocParallel)
library(edgeR)
library(limma)
library(DBI)
library(GenomicRanges)
library(matrixStats)
library(tidyverse)
args = commandArgs(trailingOnly=TRUE)

# Directories / filepaths
direc <- "" # Working directory
mvalues_path <- "" # # Directory containing m-values per chromosome

###############################
## select genes to test
###############################

min.F <- 10 # F-statistic from ANOVA should be at least this
# file with F-stats
F.values <- read.table(sprintf("%soutput.genome.wide/f.stats.txt", direc), header = TRUE, sep = "\t")
genes <- F.values$gene[F.values$f.stat > 10] # Genes that will be tested (F>10)
all.genes <- F.values$gene # All genes, regardless of F-statistic

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
probes.locs <- probes.locs[probes.locs$gene %in% all.genes]

rm(ensembl); gc()

# Load genes that associate with the same CpGs
gene_pairs <- read_tsv(sprintf("%sresultsv1/sensitivity/gene_pairs_distantgi.txt", direc))

###################################
## data
###################################

# Load covariates
cvrt <- read.table(sprintf("%scvrt.meth.txt", direc), header = TRUE, sep = "\t")

###############################
## Genetic Instruments
###############################

rdata.files <- list.files(sprintf("%soutput.genome.wide/grs.5.pcs.1se_meth", direc), pattern = "Rdata", full.names = TRUE)

# Get genetic instruments
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


###############################
## Run EWAS
###############################

genes <- unique(c(gene_pairs$gene1, gene_pairs$gene2))

# EWAS function
run <- function(current.gene) {
        require(GenomicRanges)
        require(readr)
        require(caret)

        # Get nearby genes
        corr.genes <- subsetByOverlaps(probes.locs, probes.locs[probes.locs$gene == current.gene], maxgap = 1e6)
        corr.genes <- setdiff(corr.genes$gene, current.gene)
        corr.genes <- grs[, corr.genes,drop=FALSE]

        # Genes on same chromosome (>1Mb) that are associated with same CpG-site
        gns <- gene_pairs %>% filter(gene1 == current.gene | gene2 == current.gene)
        gns <- unique(c(gns$gene1, gns$gene2))
        gns <- setdiff(gns, current.gene)
        gns <- grs[,gns, drop = FALSE]
        corr.genes <- cbind(corr.genes, gns)


        if(ncol(corr.genes) > 1) {
            cors <- cor(corr.genes)
            test <- findCorrelation(cors, cutoff = 0.95)
            if(length(test) > 0) {
                        corr.genes <- corr.genes[,-test, drop = FALSE]
                        } else {
                            corr.genes <- corr.genes
                        }

                # model matrix
        }

        Z <- model.matrix(~ corr.genes +
                    Sampling_Age +
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
        write_tsv(me,path = sprintf("%soutput.genome.wide/ewas.corrected.distant.gi/output/%s.%s.txt", direc, chr, current.gene), col_names=TRUE)
        nr.covariates <- ncol(Z)
        write.table(nr.covariates, file = sprintf("%soutput.genome.wide/ewas.corrected.distant.gi/nr_covariates/%s.txt", direc,
            current.gene), quote=FALSE, row.names=FALSE, col.names=FALSE)
}

# Run for current chromosome
chr <- args[1]

# Load M-values for the current chromosome
load(sprintf("%s/%s.mvalues.RData", mvalues_path, chr))
mvalues <- mvalues[,cvrt$uuid] # Match

# Run
Y <- t(mvalues)
rm(mvalues); gc()
null <- lapply(genes, FUN = run)
