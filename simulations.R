# setwd("/mnt/ubi/iferres/benchmark_pangenome_software")

library(simba)

ref <- "pan_genome_reference.fa"

#######################################################
## Simulate pangenomes at differents evolution times ##
#######################################################

# Simulation
dir.create('evo_times')
dir.create('evo_times/genes_fasta')
evo_times <- c(1e6, 5e6, 1e7, 5e7, 1e8, 5e8)
dout <- paste0('evo_times/genes_fasta/', sub('[+]','_',as.character(evo_times)))

replicates <- 1:5
seed <- 122
for (i in seq_along(evo_times)){
  for (j in replicates){
    sseed <- seed + j 
    set.seed(sseed)
    dir_out <- paste(dout[i], j, sep = '__')
    p <- simpg(ref, ne = evo_times[i], dir_out=dir_out)
    saveRDS(p, file = paste0(dir_out, '/sim_', basename(dir_out) ,'.RDS'))
  }
}


# Format


