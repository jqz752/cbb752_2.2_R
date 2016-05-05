#Quantifying RNA-seq Using RPKM/FPKM/TPM Based on SAM & GTF Files

*CBB752 Final Project 2.2, R card, by Julian Q Zhou*

## Objective

Given a .sam and a .gtf file, calculate either RPKM/FPKM or TPM as a measure of RNA-seq quantification for genes in .gtf based on reads from .sam

## Source code

Available [here](https://github.com/jqz752/cbb752_2.2_R)

* `r_preprocess_gtf.R`: for preprocessing a large comprehensive gtf file into a compact one containing only protein-coding genes

* `r_fpkm_tpm.R`: main script

## Sample input
* A .sam file that looks like below:

| 276945-1	| 16	| chr1	| 10560	| 255	| 36M	| *	| 0	| 0	| AAACGCAGCTCCGCCCTCGCGGTGCTCTCCGGGTCT	| IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	| XA:i:0	MD:Z:0C35	NM:i:1	YG:i:3	YC:Z:128,0,128 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 479920-1	| 16	| chr1	| 17390	| 255	| 36M	| *	| 0	| 0	| CAGGCAAGCTGACACCCGCTGTCCTGAGCCCATGTT	| IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	| XA:i:0	MD:Z:36	NM:i:0	YG:i:6	YC:Z:128,0,128 |
| 475014-1	| 0	| chr1	| 22401	| 255	| 16M	| *	| 0	| 0	| TCTGACAGGCGTACCA	| IIIIIIIIIIIIIIII	| XA:i:1	MD:Z:12G3	NM:i:1	YG:i:16	YC:Z:128,0,128 |


* A .gtf file that looks like below:

chr1 |	HAVANA |	gene	| 69091 | 70008	| . |	+	| .	| gene_id ENSG00000186092.4; transcript_id ENSG00000186092.4; gene_type protein_coding; gene_status KNOWN; gene_name OR4F5; transcript_type protein_coding; transcript_status KNOWN; transcript_name OR4F5; level 2; havana_gene OTTHUMG00000001094.1;
--- | --- | --- | --- | --- | --- | --- | --- | --- 
chr1 |	ENSEMBL|	gene	| 134901|	139379|	.	|-	|.	|gene_id ENSG00000237683.5; transcript_id ENSG00000237683.5; gene_type protein_coding; gene_status KNOWN; gene_name AL627309.1; transcript_type protein_coding; transcript_status KNOWN; transcript_name AL627309.1; level 3;
chr1 |	HAVANA	| gene	|367640	|368634	|.	|+	|.	|gene_id ENSG00000235249.1; transcript_id ENSG00000235249.1; gene_type protein_coding; gene_status KNOWN; gene_name OR4F29; transcript_type protein_coding; transcript_status KNOWN; transcript_name OR4F29; level 2; havana_gene OTTHUMG00000002860.1;


## Sample output
* A tab-delimited .txt file that looks like below:

chr |	start	| end |	id	| name	| counts |	tpm
--- | --- | --- | --- | --- | --- | --- 
chr1 |	110158726 |	110174673 |	ENSG00000116337.11 |	AMPD2 |	89.5790927065621 |	418.242573880954
chr1 |	110198703	| 110208118 |	ENSG00000168765.11 |	GSTM4 |	NA |	NA
chr1 |	110210644	| 110252171	| ENSG00000213366.8	 | GSTM2 |	148.055914670864 |	265.457974667923

## Usage
If applicable, before running `r_fpkm_tpm.R`, first run `r_preprocess_gtf.R` to convert a large comprehensive gtf file into a compact one containing only protein-coding genes.

Next, call the main function, `get.fpkm.tpm` from `r_fpkm_tpm.R` as follows:

`get.fpkm.tpm(input.sam, input.gtf, output.name, demo.mode = F,
              sam.num.header, mapq.thresh=NA, 
              save.pileup=T, save.pileup.name='pileup.Rdata', use.parallel=T, 
              count.est.mtd, count.verbose=F, quant.mtd)`

* `input.sam/gtf`: filename of input sam/gtf file (both 1-based).
* `output.name`: output filename.
* `demo.mode`: `T`/`F`; if `T`, for demonstration, skip processing input .sam and load processed .sam; load sample .gtf; skip computing pileup and load pre-stored pileup result.

* `sam.num.header`: number of lines in head section (`@...`) of .sam to skip; should be >=0.
* `mapq.thresh`: threshold for `MAPQ` in .sam; should be >=0; applicable only if `MAPQ!=255`.

* `save.pileup`, `save.pileup.name`: whether to save pileup as `save.pileup.Rdata`.
* `use.parallel`: whether to use parallel computing via the `foreach` package for computing pile up.

* `count.est.mtd`: method to estimate read counts of a gene based on raw read counts at each base position within the range of that range: `'mean'`, `'median'`, `'min'`, `'max'`, `'quantile'`; must be supplied as a list; e.g. `list('mean')`; `list('quantile', 0.7)`.
* `count.verbose`: `T`/`F`; if `T`, print out messages when estimating counts.
* `quant.mtd`: method to quantify RNA-seq; `'rpkm'`, `'fpkm'`, or `'tpm'`; `rpkm` and `fpkm` essentially give the same results.

This program assumes that the third column (`RNAME`) in the .sam file contains chromosome number. It also assumes that the reads in .sam are single-end reads. 

For each chromosome, it first computes the read depth, or pile-up, at each base position. Note that using parallel computing via `foreach` (`use.parallel=T`) to compute pile-up might not work for very large sam files as long vectors are not yet supported by `foreach`. Then, for each gene, based on the read depths of positions within the range of that gene, an estimated read count is calculated, using a method of the user's choice (`count.est.mtd`). Finally, RPKM/FPKM or TPM is calculated for each gene, depending on user's choice (`quant.mtd`), using the length of the protein-coding gene as effective length. `NA`s are reported as estimated counts and RPKM/FPKM/TPM for genes for which there is no read coverage or no sequencing data at all.

The output is a tab-delimited text file named after `output.name`, and contains 5 columns: chromosome number, start position, end position, gene ID, gene name, estimated read counts, and RPKM/FPKM/TPM (depending on `count.est.mtd`).

#### Example 1 (Demo mode)

In demo mode (`demo.mode=T`), as inputs, a pre-processed sample .sam file is loaded from `sample_Gm12878Cytosol_trimmed_sam.Rdata`; and a sample .gtf file is loaded from `sample_gencode19_prtn_coding.gtf`. To shorten the runtime, `sample_pileup.Rdata`, which stores the pile-up results computed in advance with exactly the same commands shown below, is also loaded so as to skip the waiting time of actually computing pile-up during demo mode.

As a result, two sample output files, `sample_tpm.txt` and `sample_rpkm.txt`, are produced, respectively, by running:

`get.fpkm.tpm(output.name='sample_tpm.txt', demo.mode=T, count.est.mtd=list('quantile', 0.75), quant.mtd='tpm')` , and

`get.fpkm.tpm(output.name='sample_rpkm.txt', demo.mode=T, count.est.mtd=list('quantile', 0.75), quant.mtd='rpkm')`.


#### Example 2 (Non-demo mode)

Here is an example running in non-demo mode (`demo.mode=F`), assuming that there are two input files, `Gm12878Cytosol.sam` and `sample_gencode19_prtn_coding.gtf`. 

`Gm12878Cytosol.sam` can be converted by `samtools` from `wgEncodeCshlShortRnaSeqGm12878CytosolShortAln.bam`, accessible from [here](http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeCshlShortRnaSeq/wgEncodeCshlShortRnaSeqGm12878CytosolShortAln.bam). This file is not provided here on Github due to its large size (>1 GB). Instead, a pre-processed version, `sample_Gm12878Cytosol_trimmed_sam.Rdata` is provided (see Example 1 (Demo mode)).

`sample_gencode19_prtn_coding.gtf` can be created by pre-processing `gencode.v19.annotation.gtf`, accessible at ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz, with `r_preprocess_gtf.R`. This is the same as the input gtf file used in Example 1 (Demo mode).

The output text file is the **same** as `sample_tpm.txt` produced in demo mode, since the input files used in demo mode and here are essentially the same, in addition to the same `count.est.mtd` and `quant.mtd` arguments. 

`> get.fpkm.tpm(input.sam='Gm12878Cytosol.sam', input.gtf='sample_gencode19_prtn_coding.gtf', output.name='sample_tpm.txt', demo.mode=F, sam.num.header=0, save.pileup.name='sample_pileup.Rdata', use.parallel=F, count.est.mtd=list('quantile', .75), quant.mtd='tpm')`

`[1] "reading in sam file..."`

`[1] "trimming sam file..."`

`[1] "computing pileup based on sam file (this might take a while)..."`

`[1] "now computing pileup for chr1"`

`[1] "now computing pileup for chr10"`

`[1] "now computing pileup for chr11"`

`[1] "now computing pileup for chr12"`

`[1] "now computing pileup for chr13"`

`[1] "now computing pileup for chr14"`

`[1] "now computing pileup for chr15"`

`[1] "now computing pileup for chr16"`

`[1] "now computing pileup for chr17"`

`[1] "now computing pileup for chr18"`

`[1] "now computing pileup for chr19"`

`[1] "now computing pileup for chr2"`

`[1] "now computing pileup for chr20"`

`[1] "now computing pileup for chr21"`

`[1] "now computing pileup for chr22"`

`[1] "now computing pileup for chr3"`

`[1] "now computing pileup for chr4"`

`[1] "now computing pileup for chr5"`

`[1] "now computing pileup for chr6"`

`[1] "now computing pileup for chr7"`

`[1] "now computing pileup for chr8"`

`[1] "now computing pileup for chr9"`

`[1] "now computing pileup for chrM"`

`[1] "now computing pileup for chrX"`

`[1] "saving pileup.file..."`

`[1] "reading in gtf file..."`

`[1] "computing raw and estimated counts for genes in gtf (this might take a while)..."`

`[1] "quantifying RNA-seq in tpm for genes in gtf..."`

`[1] "exporting output..."`
