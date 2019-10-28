## Rscript for submitting array jobs

# Make memory requested dependent on number of CpGs on chromosome
direc <- ""
cpgs <- read.table(sprintf("%scpgs.tested.txt", direc))

# 450K annotation
library(FDb.InfiniumMethylation.hg19)
fdb <- features(FDb.InfiniumMethylation.hg19)
fdb <- fdb[as.character(cpgs$V1),]

# Nr cpgs per chr
nrcpgs_perchr <- table(seqnames(fdb))
nrcpgs_perchr <- nrcpgs_perchr[nrcpgs_perchr > 0]

# Submit jobs (1 job for each chromosome)
for(i in names(nrcpgs_perchr)) {
	if(nrcpgs_perchr[i] < 20000) {
		mem = 7
	} else if(nrcpgs_perchr[i] < 30000) {
		mem = 8
	} else if(nrcpgs_perchr[i] > 40000) {
		mem = 10
	}

	qsub_command <- sprintf('qsub -l h_rt=30:00:00 -l h_vmem=%1$sG -N %2$s template3.sh ewas_corrected_nearby_GIs %2$s', mem, i)

	# Submit
	system(qsub_command)
}
