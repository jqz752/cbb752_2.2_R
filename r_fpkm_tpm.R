##### Julian Q. Zhou
##### compute FPKM/TPM based on .sam

#########################
##### main function #####
#########################

get.fpkm.tpm = function(input.sam, input.gtf,
                        sam.num.header, mapq.thresh=NA, 
                        save.pileup=T, save.pileup.name='pileup.Rdata',
                        use.parallel=T, count.est.mtd, count.verbose=F,
                        output.name){
  ##### input:
  # - input.sam/gtf: filename of input sam/gtf file (both 1-based)
  # - sam.num.header: number of lines in head section of .sam; >=0
  # - mapq.thresh: threshold for MAPQ; >=0; applicable only if MAPQ!=255
  # - save.pileup, save.pileup.name: whether to save pileup as save.pileup.Rdata
  # - use.parallel: parallel computing via foreach
  # - count.est.mtd: method to estimate read count: mean, median, min, max, quantile
  #   given as a list; e.g. list('mean'); list('quantile', 0.7)
  # - count.verbose: T/F; if T, print out messages when estimating counts
  # - output.name: output filename
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
  sam.nreads = nrow(sam)

  ##### compute pileup
  print('computing pileup based on sam file...')
  pileup = get.pileup.file(sam.df=sam, 
                           save=save.pileup, 
                           save.name=save.pileup.name,
                           use.parallel = use.parallel)
  
  ##### read in GTF
  print('reading in gtf file...')
  gtf = read.delim(input.gtf, header=F, as.is=T)
  
  # get gene id and gene name from V9 (attributes); drop ending semi-colon
  gtf.gene.id = c()
  gtf.gene.name = c()
  for (i in 1:nrow(gtf)){
    attributes = strsplit(gtf[i, 9], split=" ")[[1]]
    gtf.gene.id = c(gtf.gene.id, substr(attributes[2], start=1, stop=nchar(attributes[2])-1))
    gtf.gene.name = c(gtf.gene.name, substr(attributes[10], start=1, stop=nchar(attributes[10])-1))
  }
  
  gtf = cbind(gtf[, -9], gtf.gene.id, gtf.gene.name)
  rm(gtf.gene.id, gtf.gene.name, attributes, i)
  
  # keep only columns of interest
  gtf = gtf[, c(1,4,5,9,10)]
  colnames(gtf) = c('chr', 'start', 'end', 'id', 'name')
  
  ##### compute raw counts and calculate estimated counts
  gtf = cbind(gtf, counts = apply(gtf, 1, get.counts, 
                                  pileup.list = pileup.file,
                                  mtd = count.est.mtd,
                                  verbose = count.verbose))
  
  ##### compute rpkm/fpkm or tpm
  gtf = cbind(gtf, quantity = apply(gtf, 1, get.quantity, gtf.all=gtf, mtd='tpm', num.total.reads=sam.nreads)

  ##### export as tab-delimited file
  write.table(gtf, file=output.name, quote=F, row.names=F, sep="\t")
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

####################################################################
##### raw & estimated counts at each position for a given gene #####
####################################################################

get.counts = function(gtf.gene, pileup.list, mtd){
  # input:
  # - gtf.gene: a row from gtf (columns = chr, start, end, id, name)
  # - pileup.list: a list, each item is a pileup matrix for a chromosome
  # - mtd: method to estimate counts; a list of length 1 or 2
  # output:
  # a number, estimated no. of counts; or NA
  
  if (!is.null(pileup.list[[gtf.gene['chr']]])){
    # if pileup for the chromosome harboring current gene exists
    gene.idx.in.pileup = (pileup.list[[gtf.gene['chr']]][, 1]>=gtf.gene['start']) &
                         (pileup.list[[gtf.gene['chr']]][, 1]<=gtf.gene['end'])
    if (sum(gene.idx.in.pileup)>0){
      # if pileup for range corresponding to current gene exists
      raw.count = pileup.list[[gtf.gene['chr']]][gene.idx.in.pileup, 2]
      
      # compute estimated counts
      if (mtd[[1]]=='mean') {
        return(mean(raw.count))
      } else if (mtd[[1]]=='median') {
        return(median(raw.count))
      } else if (mtd[[1]]=='max'){
        return(max(raw.count))
      } else if (mtd[[1]]=='min'){
        return(min(raw.count))
      } else if (mtd[[1]]=='quantile'){
        return(quantile(raw.count, mtd[[2]]))
      }
      
    } else {
      print(paste('no coverage for gene', gtf.gene['id'], gtf.gene['name']))
      return(NA)
    }
  } else {
    print(paste('no data for the chromosome harboring gene', gtf.gene['id'], gtf.gene['name']))
    return(NA)
  }
}

#####################################
##### calculate RPKM/FPKM & TPM #####
#####################################

get.quantity = function(gtf.gene, gtf.all, mtd, num.total.reads){
  # input:
  # - gtf.gene: a row from gtf (columns = chr, start, end, id, name, counts)
  # - gtf.all: entire gtf (for calculating tpm)
  # - mtd: method to quantify, 'rpkm', 'fpkm', or 'tpm'
  # - num.total.reads: total number of reads in .sam
  # output:
  # a number, or NA
  if (is.na(gtf.gene['counts'])){
    return(NA)
  } else {
    # calculate gene length; used as effective length
    gene.length = abs(gtf.gene['end'] - gtf.gene['start'])
    
    if (mtd=='rpkm' | mtd=='fpkm') {
      # fpkm/rpkm
      rpkm = gtf.gene['counts'] / gene.length / num.total.reads * 10^9
      return(rpkm)

    } else if (mtd=='tpm') {
      # tpm
      tpm.norm.factor = get.tom.norm.factor(gtf.all)
      tpm = gtf.gene['counts'] / gene.length / tpm.norm.factor * 10^6
      return(tpm)
    }
  }
}

get.tpm.norm.factor = function(gtf.all){
  # input:
  # - gtf.all: all rows from gtf (columns = chr, start, end, id, name, counts)
  # output:
  # normalizatioin factor for tpm: sum_j{X_j/l_j}, where X_j = counts, l_j=eff.length
  
  gtf.all.non.na = gtf.all[!is.na(gtf.all$counts), ]
  gtf.all.non.na.length = abs(gtf.all.non.na$end - gtf.all.non.na$start)
  norm.factor = sum(gtf.all.non.na$counts / gtf.all.non.na.length, na.rm=T)
  return(norm.factor)
}