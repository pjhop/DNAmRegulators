# Submit EWAS jobs
args <- commandArgs(trailingOnly=TRUE)
N = args[1]

for(n in 1:N) {
	qsub_command <- sprintf('qsub -l h_rt=18:00:00 -l h_vmem=60G -N ewas_%2$s template3.sh ewas %1$s %2$s', N, n)

	# Submit
	system(qsub_command)
}
