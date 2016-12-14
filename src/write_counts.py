__author__ = 'Jonathan Rubin'

#Writes a counts.sh script file and then runs it to get read coverage over
#specified intervals for each given bam file

import os

def run(bamdir,bamlist,scriptdir,bedfile,tempdir):
    outfile = open(scriptdir + '/counts.sh')
    outfile.write("### Run in desired queue\n")
    outfile.write("#PBS -q long8gb\n")
    outfile.write("### Use the bourne shell\n")
    outfile.write("#PBS -S /bin/sh\n")
    outfile.write("### Specify the number of nodes and processors for your job\n")
    outfile.write("#PBS -l nodes=1:ppn=1\n")
    outfile.write("#PBS -m ae\n")
    outfile.write("#PBS -M joru1876@colorado.edu\n")
    outfile.write("### Retrieve/use all modules loaded ###\n")
    outfile.wrtie("#PBS -V\n")
    outfile.write("/opt/bedtools/2.16.2/multiBamCov -bams ")
    for bamfile in bamlist:
        outfile.write(bamdir + bamfile + " ")
    outfile.write("-bed $bedfile >$bedfile.count.bed")
    outfile.close()
    
    os.system("qsub -v bedfile=" + bedfile + " -N refgenecount "+scriptdir + "/counts.sh > " + tempdir + "/Job_ID.txt")