__author__ = 'Jonathan Rubin'

#Runs DE-Seq to identify differentially transcribed regions on the genome.
#Requires:  1. Change bamdir to full path to bams (needs '/' at the end)
#           2. Change bamlist to have filenames for all bam files to be analyzed (comma separated)
#           3. Change bedfile to full path to bed file containing intervals to be analyzed
#           4. module load bedtools_2.16.2
#           5. Needs to be run on Vieques (in CU Boulder) because Vieques has DE-Seq installed

import os
import write_counts
import check_job
import add_header

#User-defined input
################################################################################
#Directory containing sorted bam files
bamdir = '/scratch/Users/joru1876/GROSeqRaw/flipped/bowtie2/sortedbam/'
#List of bam files in directory
bamlist = ['JDR_DMSO_SS102217_093015_CAGATC_L005_R1_001.flip.fastqbowtie2.sorted.bam','JDR_CA_SS102217_093015_ACTTGA_L005_R1_001.flip.fastqbowtie2.sorted.bam']
#Full path to bed file containing intervals to be analyzed
bedfile = '/scratch/Users/joru1876/GROSeqRaw/flipped/bowtie2/sortedbam/bedfiles/refGene.bed'
################################################################################

#Return parent directory
def parent_dir(directory):
    pathlist = directory.split('/')
    newdir = '/'.join(pathlist[0:len(pathlist)-1])
    
    return newdir

#Home directory
homedir = os.path.dirname(os.path.realpath(__file__))

#Scripts directory
scriptdir = parent_dir(homedir) + '/scripts'

#Temporary files directory
tempdir = parent_dir(homedir) + '/temp'

#Full path to temporary job ID file
job = tempdir + "/Job_ID.txt"

def run():
    write_counts.run(bamdir,bamlist,scriptdir,bedfile,tempdir)
    check_job.run(job,tempdir)
    
    headerbed = add_header.run(bedfile)
    
    return