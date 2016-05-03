# -*- coding: utf-8 -*-
"""
Created on Mon May  2 14:02:44 2016

@author: kevin
"""
#input is sample.sam
Samfname = "sample.sam"
with open(Samfname) as f:
    Sam = f.readlines()

lines = []

for i in range(len(Sam)):
    lines.append(Sam[i].split("\t", Sam[i].count("\t")))

    #input is sample.gtf
Gfffname = "sample.gtf"
with open(Gfffname) as f:
    Gff = f.readlines()

GFFlines = []

for i in range(len(Gff)):
    GFFlines.append(Gff[i].split("\t", Gff[i].count("\t")))
totalReads = 0
mappedReads = {}
GFFlengths = {}
gffCount = 0
for GFFline in GFFlines:
    gffCount +=1
    GFFrange = set(range(int(GFFline[2]),int(GFFline[3])))
    GFFlengths[str(gffCount)] = len(GFFrange)
    for i in range(2, len(lines)):
        totalReads += 1
        rangeOfLine = range(int(lines[i][3]), int(lines[i][3]) + len(lines[i][9]))
        if len(GFFrange.intersection(rangeOfLine))>0:
            if str(gffCount) in mappedReads:
                #mappedReads[lines[i][2]] += 1
                mappedReads[str(gffCount)] += 1
            else:
                #mappedReads[lines[i][2]] = 1
                mappedReads[str(gffCount)] = 1

#Calculate TPM:
#(X[i]/L[i])*(1/(sum(Xtot/Ltot)))*1e6
#breaking it down to steps:
#Take each gene read and divide by length (in kb) of gene
#Next take these normalized reads and sum them for each gene
#create scaling factor by dividing summed normalized reads by 1e6
#take your normalized reads from step 1 and divide by the scaling factor
#instructions from here: https://www.youtube.com/watch?v=TTUrtCY2k-w
#and here: https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/
GFFrpk = {}
ScalingFactor = 0#sum and divided by 1000000
for i in range(len(GFFlengths)):
    #calculate rpk by mappedRead/length
    GFFrpk[str(i+1)] = float(mappedReads[str(i+1)])/GFFlengths[str(i+1)]
    #get scaling factor by doing running total of the rpk
    ScalingFactor += (GFFrpk[str(i+1)])#get total rpk
    
#scale the total rpk by 1e6
ScalingFactor /= 1000000
tpm = {}
rpkm = {}
#calculate RPKM: 
#Get number of reads
#Calculate Number_of_reads/(length_of_gene/kilobase * total_number_reads/1e6)
#Instructions found here: https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/

for i in range(len(mappedReads)):
    tpm[str(i+1)] = mappedReads[str(i+1)]/ScalingFactor
    rpkm[str(i+1)] = float(mappedReads[str(i+1)])/((totalReads/1000000.) * (GFFlengths[str(i+1)]/1000.))

for i in range(len(mappedReads)):
    print("Gene: " + repr(GFFlines[i][0]) + ", Feature: " + repr(GFFlines[i][1]) +
    ", TPM: " +  repr(round(tpm[str(i+1)], 3)) + ", RPKM: " +  repr(round(rpkm[str(i+1)], 3)))
