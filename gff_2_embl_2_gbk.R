# from gff to genbank (and embl)
#> singularity pull docker:sangerpathogens/gff3toembl
gff3_2_embl <- function(gff3, 
                        authors = 'John', 
                        title = 'Some Title',
                        publication = 'Some Journal', 
                        genome_type = 'circular', 
                        classification = 'PROK',
                        organism = 'Organism', 
                        taxonid = 1234, 
                        accession = 'project',
                        description = 'description'){
  
  arg <- paste0("--authors ", "\'",authors,"\'", ' ',
                "--title ", "\'", title, "\'", ' ',
                "--publication ", "\'", publication, "\'", ' ',
                "--genome_type ","\'" , genome_type, "\'", ' ',
                "--classification ","\'" , classification, "\'", ' ',
                "--output_filename /dev/stdout", ' ',
                organism,  ' ',
                taxonid,  ' ',
                "\'", accession, "\'",  ' ',
                "\'", description, "\'", ' ',
                gff3)
  
  cmd <- 'singularity exec /mnt/ubi/iferres/singularity_images/gff3toembl_latest.sif gff3_to_embl'
  
  stout <- system(paste(cmd, arg), intern = T)
  
  deli <- grep('//', stout)
  
  fctr <- rep(seq_along(deli), c(deli[1], diff(deli)))
  
  chunks <- split(stout, fctr)
  
  # fields to add 
  DT1 <- 'DT   27-MAR-2019 (Rel. 01, Created)'
  DT2 <- 'DT   27-MAR-2019 (Rel. 02, Last updated, Version 1)'
  genomeid <- sub('[.]gff$', '', basename(gff3))
  OS <- paste0('OS   Fakeum bacterii (', genomeid, ')')
  
  chunks <- lapply(chunks, function(x){
    
    ln <- length(x)
    out <- vector('character', ln + 4 + 3)
    gpac <- rev(grep('^AC', x))[1]
    out[1:(gpac+1)] <- x[1:(gpac+1)]
    out[(gpac+2):(gpac+4)] <- c(DT1, DT2, 'XX')
    gpde <- rev(grep('^DE', x))[1]
    rg <- (gpac+2):(gpde+1)
    lim <- (gpac+5 + length(rg)-1)
    out[(gpac+5):lim] <- x[rg]
    ini <- lim+1
    lim <- ini + 1
    KW <- 'KW   function'
    out[ini:lim] <- c(KW, 'XX')
    ini <- lim+1
    lim <- ini+1
    out[ini:lim] <- c(OS, 'XX')
    ini <- lim+1
    lim <- ln + 4 + 3
    out[ini:lim] <- x[(gpde+2):length(x)]
    out
    
  })
  
  
  
  embl <- sub('[.]gff$', '.embl', gff3)
  if (file.exists(embl)){
    file.remove(embl)
  }
  cat(paste0('Writing embl at ', embl, '...'))
  lapply(chunks, cat, file = embl, sep = '\n', append=TRUE)
  cat(' DONE!\n')
  chunks
}

embl_2_gbk <- function(x, file){
  
  if (missing(file)) file <- stdout()
  #Line 1
  contig <- sub('[^contig\\d+]+', '', grep('^AC.+contig\\d+', x, value = TRUE))
  large <- regmatches(x[1], regexpr('\\d+ BP.', x[1]))
  large <- sub('BP[.]', 'bp', large)
  L1 <- paste('LOCUS', contig, large, 'DNA', 'linear', '01-APR-2019', sep='\t')
  
  #Line 2
  L2 <- paste('DEFINITION\tGenus species strain strain.')
  
  L3 <- paste('ACCESSION\t')
  
  L4 <- paste('VERSION\t')
  
  #Line 5
  KW <- sub('^KW[ ]+','',grep('^KW[ ]+', x, value = TRUE))
  L5 <- paste('KEYWORDS', KW, sep = '\t')
  
  #Line 6
  organism <- sub('^OS[ ]+','',grep('^OS[ ]+', x, value = TRUE))
  L6 <- paste('SOURCE', organism, sep = '\t')
  
  L7 <- paste('  ORGANISM', organism, sep = '\t')
  
  L8 <- paste('COMMENT', 'This is a fake GenBank file.')
  
  features <- gsub('FH|[ ]|Key','',grep('^FH[ ]+Key[ ]+', x, value = TRUE))
  L9 <- paste('FEATURES', features, sep = '\t')
  
  L10 <- sub('^FT','',grep('^FT', x, value = T))
  
  
  #L11: translation 
  
  sq <- grep('^SQ', x)
  L12 <- 'ORIGIN'
  
  dna <- x[(sq+1):(length(x)-1)]
  string <- paste0(gsub('[ ]|\\d', '', dna), collapse = '')
  trans <- paste0(seqinr::translate(strsplit(string, '')[[1]]), collapse = '')
  trans <- sub('[*]', '', trans)
  translation <- paste0('/translation="', trans, '"')
  splitInParts <- function(string, size){
    pat <- paste0('(?<=.{',size,'})')
    strsplit(string, pat, perl=TRUE)[[1]]
  }
  translation <- splitInParts(translation, 59)
  space <- sub('[^ ]+', '', L10)[2]
  L11 <- paste0(space, translation, sep = '')
  
  counts <- cumsum(c(1, rep(60, length(dna)-1)))
  mx <- max(nchar(counts))
  dna <- sub('[ ]+\\d+', '', dna)
  dna <- sub('^[ ]+','',dna)
  spp <- sapply((9-nchar(counts)), function(x){
    paste0(rep(' ', x), collapse = '')
  })
  L13 <- paste(paste(spp, counts, ' ', sep=''), dna, sep='')
  L14 <- x[length(x)]
  
  tofile <- c(L1, L2, L3, L4,
              L5, L6, L7, L8, 
              L9, L10, L11, L12, 
              L13, L14)
  
  cat(tofile, file = file,  sep = '\n', append = TRUE)
}



gff_2_embl_2_gbk <- function(gff){
  
  for (i in 1:length(gff)){
    embl_chunks <- gff3_2_embl(gff[i])
    gbk_out <- sub('gff$', 'gbk', gff[i])
    lapply(embl_chunks, embl_2_gbk, gbk_out)
    cat(paste0(gff[i], '\n'))
  }
  
}




