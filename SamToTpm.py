# -*- coding: utf-8 -*-
"""
Created on Mon May  2 14:02:44 2016

@author: kevin
"""

Samfname = "sample.sam"
with open(Samfname) as f:
    Sam = f.readlines()

lines = []

for i in range(len(Sam)):
    lines.append(Sam[i].split("\t", Sam[i].count("\t")))

    
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

GFFrpk = {}
ScalingFactor = 0#sum and divided by 1000000
for i in range(len(GFFlengths)):
    GFFrpk[str(i+1)] = float(mappedReads[str(i+1)])/GFFlengths[str(i+1)]
    ScalingFactor += (GFFrpk[str(i+1)])#get total rpk

ScalingFactor /= 1000000#scale by 1000000
tpm = {}
rpkm = {}
for i in range(len(mappedReads)):
    tpm[str(i+1)] = mappedReads[str(i+1)]/ScalingFactor
    rpkm[str(i+1)] = float(mappedReads[str(i+1)])/((totalReads/1000000.) * (GFFlengths[str(i+1)]/1000.))
print("TPM: ", tpm, "RPKM: ", rpkm)
