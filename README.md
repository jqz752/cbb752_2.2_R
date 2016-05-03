# CBB752 Final Project 2.2, R card, by Julian Q Zhou

## Objective

Given a .sam and a .gtf file, calculate either RPKM/FPKM or TPM as a measure of RNA-seq quantification for genes in .gtf based on reads from .sam

## Source code

Available [here](https://github.com/jqz752/cbb752_2.2)

`r_preprocess_gtf.R`: for preprocessing a large comprehensive gtf file into a compact one containing only protein-coding genes

`r_fpkm_tpm.R`: main script

## Sample input
* A .sam file

| 276945-1	| 16	| chr1	| 10560	| 255	| 36M	| *	| 0	| 0	| AAACGCAGCTCCGCCCTCGCGGTGCTCTCCGGGTCT	| IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	| XA:i:0	MD:Z:0C35	NM:i:1	YG:i:3	YC:Z:128,0,128 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 479920-1	| 16	| chr1	| 17390	| 255	| 36M	| *	| 0	| 0	| CAGGCAAGCTGACACCCGCTGTCCTGAGCCCATGTT	| IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII	| XA:i:0	MD:Z:36	NM:i:0	YG:i:6	YC:Z:128,0,128 |
| 475014-1	| 0	| chr1	| 22401	| 255	| 16M	| *	| 0	| 0	| TCTGACAGGCGTACCA	| IIIIIIIIIIIIIIII	| XA:i:1	MD:Z:12G3	NM:i:1	YG:i:16	YC:Z:128,0,128 |


* A .gtf file

chr1 |	HAVANA |	gene	| 69091 | 70008	| . |	+	| .	| gene_id ENSG00000186092.4; transcript_id ENSG00000186092.4; gene_type protein_coding; gene_status KNOWN; gene_name OR4F5; transcript_type protein_coding; transcript_status KNOWN; transcript_name OR4F5; level 2; havana_gene OTTHUMG00000001094.1;
--- | --- | --- | --- | --- | --- | --- | --- | --- 
chr1 |	ENSEMBL|	gene	| 134901|	139379|	.	|-	|.	|gene_id ENSG00000237683.5; transcript_id ENSG00000237683.5; gene_type protein_coding; gene_status KNOWN; gene_name AL627309.1; transcript_type protein_coding; transcript_status KNOWN; transcript_name AL627309.1; level 3;
chr1 |	HAVANA	| gene	|367640	|368634	|.	|+	|.	|gene_id ENSG00000235249.1; transcript_id ENSG00000235249.1; gene_type protein_coding; gene_status KNOWN; gene_name OR4F29; transcript_type protein_coding; transcript_status KNOWN; transcript_name OR4F29; level 2; havana_gene OTTHUMG00000002860.1;

## Sample output
* A tab-delimited .txt file

chr |	start	| end |	id	| name	| counts |	tpm
--- | --- | --- | --- | --- | --- | --- 
chr1 |	110158726 |	110174673 |	ENSG00000116337.11 |	AMPD2 |	89.5790927065621 |	418.242573880954
chr1 |	110198703	| 110208118 |	ENSG00000168765.11 |	GSTM4 |	NA |	NA
chr1 |	110210644	| 110252171	| ENSG00000213366.8	 | GSTM2 |	148.055914670864 |	265.457974667923

## Usage
If applicable, before running `r_fpkm_tpm.R`, first run `r_preprocess_gtf.R` to convert a large comprehensive gtf file into a compact one containing only protein-coding genes.

Next, call the main function, `get.fpkm.tpm` from `r_fpkm_tpm.R` as follows:

`get.fpkm.tpm = function(input.sam, input.gtf,
                        sam.num.header, mapq.thresh=NA, 
                        save.pileup=T, save.pileup.name='pileup.Rdata',
                        use.parallel=T, count.est.mtd, count.verbose=F,
                        quant.mtd, output.name)`

* `input.sam/gtf`: filename of input sam/gtf file (both 1-based)
* `sam.num.header`: number of reads in head section of .sam; should be >=0; should exclude head section (`@...`)
* `mapq.thresh`: threshold for `MAPQ` in .sam; should be >=0; applicable only if `MAPQ!=255`
* `save.pileup`, `save.pileup.name`: whether to save pileup as `save.pileup.Rdata`
* `use.parallel`: whether to use parallel computing via the `foreach` package for computing pile up
* `count.est.mtd`: method to estimate read count: `'mean'`, `'median'`, `'min'`, `'max'`, `'quantile'`; must be supplied as a list; e.g. `list('mean')`; `list('quantile', 0.7)`
* `count.verbose`: `T`/`F`; if `T`, print out messages when estimating counts
* `quant.mtd`: method to quantify RNA-seq; `'rpkm'`, `'fpkm'`, or `'tpm'`
* `output.name`: output filename

This program assumes that the third column (`RNAME`) in the .sam file contains chromosome number.

`NA` will be reported as estimated counts and RPKM/FPKM/TPM for genes for which there is no read coverage or no sequencing data at all.
