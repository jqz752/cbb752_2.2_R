setwd("~/University/2016 Spring/cbb752/project/rnaseq")
#setwd("/home2/qz93/cbb752")

#########################
##### main function #####
#########################

get.fpkm.tpm = function(input.sam, input.bed=NULL, input.gtf=NULL,
                        sam.num.header, mapq.thresh=NA, 
                        save.pileup=T, save.pileup.name='pileup.Rdata',
                        use.parallel=T){
  ##### input:
  # - input.sam/bed/gtf: filename of input sam/bed/gtf file
  # - sam.num.header: number of lines in head section of .sam; >=0
  # - mapq.thresh: threshold for MAPQ; >=0; applicable only if MAPQ!=255
  # - save.pileup, save.pileup.name: whether to save pileup as save.pileup.Rdata
  # - use.parallel: parallel computing via foreach
  # assumes that 3rd col (RNAME) in sam file contains chromosome number
  
  ##### output:
  
  ##### read in sam file #####
  # skip head section (lines beginning with @)
  # keep characters as if (as opposed to converting to factor)
  print('reading in sam file...')
  if (sam.num.header==0){
    sam = read.delim(input.sam, header=F, as.is=T)
  } else {
    sam = read.delim(input.sam, header=F, as.is=T, skip=sam.num.header)
  }
  
  ##### trim sam file #####
  ### keep only FLAG, RNAME, POS, MAPQ, and SEQ
  print('trimming sam file...')
  sam = sam[, c(2:5, 10)]
  colnames(sam) = c('flag', 'rname', 'pos', 'mapq', 'seq')
  
  ### filter by mapping quality (MAPQ)
  # only apply filter if MAPQ is available (255=unavail)
  if (!is.na(mapq.thresh)){
    sam = sam[(sam$mapq > mapq.thresh) | (sam$mapq != 255), ]
  }
  
  ### keep only reads with FLAG values 0 or 16
  sam = sam[(sam$flag==0) | (sam$flag==16), ]
  
  ##### compute pileup
  print('computing pileup based on sam file...')
  pileup = get.pileup.file(sam.df=sam, 
                           save=save.pileup, 
                           save.name=save.pileup.name,
                           use.parallel = use.parallel)
  
  return(0)
}

#########################################################
##### get pileup by chromosome for a given sam file #####
#########################################################

get.pileup.file = function(sam.df, save=F, 
                           save.name='pileup.Rdata', use.parallel=T){
  ##### input:
  # - data frame created after trimming sam file
  # assumes these columns exist: flag, rname, pos, mapq, seq
  # assumes that rname contains chromosome number
  # - save, save.name: save output as save.name.Rdata
  # - use.parallel: use parallel computing via foreach
  
  ##### output:
  # a list, each item is a pileup matrix for a chromosome
  # save.name.Rdata, if save is set to T
  
  ##### get number of chromosomes
  # format: chrx, were x is the chromosome number
  uniq.chr = names(table(sam.df$rname))
  #uniq.chr = sapply(uniq.chr, function(x){substring(x, first=4)})
  num.chr = length(uniq.chr)
  
  ##### calculate pileup for each chromosome
  # store result for each chromosome in a list
  if (!use.parallel) {
    # non-parallel
    pileup.file = list()
    for (i in 1:num.chr){
      print(paste('now computing pileup for', uniq.chr[i]))
      pileup.file = c(pileup.file, 
                      list(get.pileup.chr(sam.df, uniq.chr[i])))
    }
  } else {
    # parallel computing
    require(doMC)
    require(foreach)
    num.cores = detectCores() - 1
    registerDoMC(num.cores)
    
    pileup.file = foreach(i=1:num.chr) %dopar% {
      print(paste('now computing pileup for', uniq.chr[i]))
      get.pileup.chr(sam.df, uniq.chr[i])
    }
  }
  
  # assign name of chromosome to each list
  names(pileup.file) = uniq.chr
  
  ##### save?
  if (save){
    save(pileup.file, file=save.name)
  }
  
  return(pileup.file)
}

#######################################################################
##### get pileup for a given chromosome based on a given sam file #####
#######################################################################

get.pileup.chr = function(sam.df, chr){
  ##### input:
  # - sam.df: data frame created after trimming sam file
  # - chr: name of chromosome
  # assumes these columns exist: flag, rname, pos, mapq, seq
  # assumes that rname contains chromosome number
  
  ##### output:
  # pileup: a matrix with 2 cols
  # pos (position), and count (# reads)
  
  ##### get reads
  reads = sam.df[sam.df$rname==chr, ]
  
  ##### find smallest and largest POS
  # POS: 1-based leftmost mapping POSition of the first matching base
  pos.min = min(reads$pos)
  pos.max = max(reads$pos)
  
  ##### find length of the longest sequence starting at pos.max
  reads.pos.max.seq = reads[reads$pos==pos.max, 'seq']
  reads.pos.max.leng = sapply(reads.pos.max.seq, nchar)
  reads.pos.max.leng.longest = max(reads.pos.max.leng)
  
  ##### create pileup matrix
  # number of positions:
  # e.g. suppose that T is at pos.max, 
  # and that the longest read starting at T has lenght 8:
  # 12       18 19 20 21 22 23 24 25 
  # X  XXXXX T  X  X  X  X  X  X  X 
  # then, total no. of positions is 18-12+1 + (8-1) = 14
  pileup = matrix(0, ncol=2,
                  nrow=pos.max-pos.min+1+(reads.pos.max.leng.longest-1))
  colnames(pileup) = c('pos', 'count')
  pileup[, 'pos'] = pos.min:(pos.max+reads.pos.max.leng.longest-1)
  
  ##### fill in pileup matrix
  # increment of 1 in 'count' at each position covered by read
  for (i in 1:nrow(reads)){
    cur.read.start.pos = reads$pos[i]
    cur.read.leng = nchar(reads$seq[i])
    cur.read.end.pos = cur.read.start.pos + cur.read.leng - 1
    cur.read.idx = match(cur.read.start.pos:cur.read.end.pos, 
                         pileup[, 'pos'])
    pileup[cur.read.idx, 'count'] = pileup[cur.read.idx, 'count']+1
  }
  
  return(pileup)
}
