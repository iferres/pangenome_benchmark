setwd('evo_times/')

library(seqinr)
fi <- list.files(path='.', 
                 pattern = 'gene\\d+[.]fasta$', 
                 full.names = TRUE,
                 recursive = TRUE)
dir.create('genomes_fasta')

for (i in seq_along(fi)){
  banam <- basename(fi[i])
  newdirpath <- sub('^[.]', './genomes_fasta', sub(banam, '',fi[i]))
  if (!dir.exists(newdirpath)){
    dir.create(newdirpath, recursive = TRUE)
  }
  
  rf <- read.fasta(fi[i],seqtype = 'DNA', as.string = TRUE)
  nn <- names(rf)
  gnm <- sub('_.+$', '', nn)
  he <- paste0('>', nn)
  sq <- unlist(rf)
  
  for (j in seq_along(gnm)){
    
    rr <- paste0(he[j], '\n', sq[j], '\n')
    geno <- paste0(newdirpath, gnm[j], '.fasta')
    cat(rr, file = geno, append = TRUE)
    
  }
  
}
