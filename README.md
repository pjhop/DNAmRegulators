## Genome-wide identification of genes regulating DNA methylation using genetic anchors for causal inference

Overview of the scripts that we used for the main analyses.
Note that several scripts are based on work
by Rene Luijk (https://git.lumc.nl/r.luijk/DirectedGeneNetworks/tree/master/).

### Construct genetic instruments (GIs)
These scripts were used to construct genetic instruments (GIs) that are used
as proxies (causal anchors) for gene expression.

**construct_GIs.R:** Loads genotype and gene expression data and constructs genetic instruments (in the training set) for each
gene using genetic variants within 100kb of a gene's TSS or TES.
LASSO as implemented in the *glmnet* package is used to perform variable selection.

**fstats.R:** Calculate the predictive power (f-statistic and R<sup>2</sup>) of the genetic instruments in the test set.

### Identify *trans* DNA methylation effects

These scripts were used to identify gene-CpG associations, where we use
the genetic instruments as proxies for gene expression.
We first identify associations between each gene and all *trans* CpG-sites
without adjusting for LD/pleiotropy among neighbouring genes.
Then, we perform the analyses while correcting for LD/pleiotropy among neighbouring genes.
Finally, we run two sensitivity analyses: adjusting for nearby SNPs that influence white blood cell composition
and adjusting each GI for other (distant) significant GIs that influence the same CpG-site(s).

#### Not corrected for LD/Pleiotropy among neighboring genes
**ewas.R**: Performs an EWAS (CpG ~ GI + covar) for each gene.
Scripts by Sikorska *et al.* are used to perform fast EWASs.
(http://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-14-166)

**bacon.R**: Adjust each EWAS for bias/inflation using the *bacon* method.

**pairs.R**: Extract significant GI-CpG pairs (Bonferroni threshold) from bacon-corrected results.

#### Corrected for LD/Pleiotropy
**ewas_corrected_nearby_GIs.R**: Performs an EWAS (CpG ~ GI + nearby GIs (<1Mb) + covar) for each gene (this script expects a chromosome as argument, so that it can be run in parallel over all chromosomes).

**bacon_corrected_nearby_GIs.R**: Adjust each EWAS for bias/inflation using the *bacon* method.

**pairs_corrected_nearby_GIs.R**: Extract GI-CpG pairs from bacon-corrected results.

**get_corrected_results.R**: Loads corrected results, subsets GI-CpG pairs that remain significant and test whether
significant GIs are predictive of neighbouring significant genes (*F* > 5).
If *F*>5 and the tested gene shares target CpGs with the neighbouring gene (at a gene-level threshold), the gene is excluded.


#### Correct for nearby SNPs known to influence white blood cell composition

**get_SNPLocations_cellcounts.R**: Extract SNPs from the Roederer *et al.* and Orru *et al.* studies and check which SNPs are present in our data.

**extract_genotypes_cellcounts.R**: Extract cellcount-SNPs from genotype data.

**ewas_corrected_nearby_GIs_nearby_cellcount_SNPs.R**: Performs an EWAS (CpG ~ GI + nearby GIs + nearby cellcount SNPs + covar) for each gene.

**bacon_corrected_nearby_GIs_nearby_cellcount_SNPs.R**: Adjust each EWAS (adjusted for nearby SNPs known to influence cellcounts) for bias/inflation using the *bacon* method.

**pairs_corrected_nearby_GIs_nearby_cellcount_SNPs.R**: Extract GI-CpG pairs from bacon-corrected results.

#### Correct for distant (significant) GIs

**check_distant_GIs.R**
Check which significant GIs are associated with the same CpG and whether these GIs are predictive
of the expression of their respective genes.

**ewas_corrected_nearby_GIs_nearby_cellcount_SNPs.R**: Performs an EWAS (CpG ~ GI + nearby GIs + distant GIs + covar) for each gene.

**bacon_corrected_nearby_GIs_nearby_cellcount_SNPs.R**: Adjust each EWAS for bias/inflation using the *bacon* method.

**pairs_corrected_nearby_GIs_nearby_cellcount_SNPs.R**: Extract GI-CpG pairs from bacon-corrected results.

#### Final Results

**get_final_results.R:**
Remove gene-CpG pairs that become insignificant in the sensitivity analyses
(adding cellcount SNPs and/or adjusting for (distant) GIs that are associated with same CpG).
This results in the final dataset.

### Miscellaneous

#### Submitting jobs on SGE cluster

**submit_ewas.R**: Submit EWAS in partitions on SGE.  
**submit_ewas_corrected_nearby_GIs.R**: Submit EWAS (corrected for nearby GIs) per chromosome on SGE cluster.
