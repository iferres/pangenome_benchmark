setwd('/mnt/ubi/iferres/benchmark_pangenome_software/')
nwdir <- 'evo_times/genomes_gff/'
dir.create(nwdir)

require(seqinr)
ffn2gff <- function(ffn, dout = '.'){
  
  rf <- read.fasta(ffn, 
             seqtype = 'DNA', 
             as.string = TRUE, 
             forceDNAtolower = TRUE)
  
  li1 <- '##gff-version 3\n'
  
  contig <- 1:length(rf)
  len_sq <- sapply(rf, nchar)
  li2 <- paste0('##sequence-region contig.', contig, ' 1 ', len_sq, '\n')
  
  nn <- names(rf)
  li3 <- paste0('contig.', contig, 
                '\tProdigal:2.6\tCDS\t1\t',
                len_sq,
                '\t.\t+\t0\tID=',nn, 
                ';locus_tag=',nn, 
                ';gene=',nn,
                ';product=',nn,'\n') 
  li4 <- '##FASTA\n'
  li5 <- paste0('>contig.',contig,'\n',unlist(rf),'\n')
  
  o <- paste0(dout,basename(sub('[.]f\\w+$', '.gff', ffn)))
  cat(c(li1, li2, li3, li4, li5), file = o, sep = '')
  
  return(o)
}

fis <- list.files(path = 'evo_times', 
	pattern='genome\\d+[.]fasta$', 
	recursive=T, 
	full.names=T)
bn <- basename(dirname(fis))

for (i in 1:length(fis)){
	dout <- paste0(nwdir, bn[i], '/')
	if (!dir.exists(dout)) dir.create(dout)
	ffn2gff(fis[i], dout)
}




