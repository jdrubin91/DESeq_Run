__author__ = 'Jonathan Rubin'

import os

def run(bedfile):
    outfile = open(bedfile + ".count.header.bed",'w')
    with open(bedfile + ".count.bed") as F:
        line1 = F.readline()
        for i in range(len(line1.strip().split())):
            outfile.write("column " + str(i) + "\t")
        outfile.write(line1)
        for line in F:
            outfile.write(line)
            
    return bedfile + ".count.header.bed"