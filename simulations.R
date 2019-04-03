# setwd("/mnt/ubi/iferres/benchmark_pangenome_software")

library(simba)

ref <- "pan_genome_reference.fa"

#######################################################
## Simulate pangenomes at differents evolution times ##
#######################################################

# Simulation
dir.create('evo_times')
evo_times <- c(1e6, 5e6, 1e7, 5e7, 1e8, 5e8)
dout <- paste0('evo_times/', sub('[+]','_',as.character(evo_time)))

for (i in seq_along(evo_time)){
  set.seed(i)
	p <- simpg(ref, ne = evo_time[i], dir_out=dout[i])
	saveRDS(p, file = paste0(dout[i], '/sim_', basename(dout[i]) ,'.RDS'))
}


# Format


