##### pre-process GTF from Gencode Reference V19 and keep only protein-coding genes

##### perform the following before running this script:

# download gencode v19 from 
# ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz
# unzipped gtf and rename as 'gencode19.gtf'

### keep only genes
gencode19 = read.delim('gencode19.gtf', header=F, skip=5)
gencode19 = gencode19[gencode19$V3=='gene',]
gencode19[, 9] = as.character(gencode19[, 9])

### keep only protein-coding genes
gene.type = vector()
for (i in 1:nrow(gencode19)){
  gene.type = c(gene.type, strsplit(gencode19[i, 9], split=' ')[[1]][6])
}

length(gene.type) == nrow(gencode19)

table(gene.type)

gencode19.prtn.coding = gencode19[gene.type=='protein_coding;', ]
dim(gencode19.prtn.coding) # 20345 x 9

### export as .Rdata and tab-delimited .gtf
save(gencode19.prtn.coding, file='sample_gencode19_prtn_coding.Rdata')
write.table(gencode19.prtn.coding, file='sample_gencode19_prtn_coding.gtf', quote=F, sep="\t", row.names=F)

##### proceed to manually delete first row of .gtf (V1 V2 ... V9)