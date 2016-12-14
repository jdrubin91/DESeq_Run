import sys


#to use
#python DESeq_Rscript.py <infile> <column1> <column2> <type_transcript> <condition1> <condition2>
#for example
#python DESeq_Rscript.py /projects/dowellde/groseq/data/set1/clipped_fastqM10/samfiles/sortedbamfiles/gffcoverage/compare_cov_fileless1neg_istead.txt  11 12 gffcoverage DMSO Nutlin
#Rscript file name and location will output
#after you make the R script run using /opt/R/2.15.1/bin/R CMD BATCH --no-save --no-restore <outfilename>


def create_DESeqRscript_replicates(infile="/projects/dowellde/groseq/data/replicates/gffcoverage/set1andset2.coverage.protein_coding",columns="[10, 15, 11, 14]", type_transcript="gffcoverage", conditions="['DMS0', 'DMSO', 'Nutlin', 'Nutlin']", condition1="DMSO", condition2="Nutlin",title_of_names_column="group"):
        """infile should be a compare coverage file. columns start counting at 1"""

        f = open(infile)
        headers = f.readline()
        headers = headers.strip("\n")
        headers = headers.split("\t")
        f.close()
        infile_dir = infile.split("/")[:-1]
        infile_dir = "/".join(infile_dir)+"/"
        infile_root = infile.split("/")[-1].strip(".txt")
	set_conditions = set(eval(conditions))
	set_conditions = list(set_conditions)
        outfile = infile_dir+infile_root+"."+condition1+condition2+type_transcript
        write_file = outfile+".R"
        print write_file
        wf = open(write_file ,"w")
        R_dump_file = outfile+".Rout"
        graph_file = outfile+".png"
        outfileallinputs = outfile+".res.txt"
        outfilesig = outfile+".resSig.txt"
        outfilesig_orderpval = outfile+".resSig_pvalue.txt"
        wf.write('sink("'+R_dump_file+'")\n')
        wf.write('library( DESeq )\n')
        wf.write('data <- read.delim("'+infile+r'", sep="\t", header=TRUE)'+"\n")#need to check that \t comes out like it should. Might write it wrong.
	columns_list = []
	columns = eval(columns)
	line = ", ".join(map(str,columns))
        wf.write('countsTable <- subset(data, select=c('+line+'))\n')
        wf.write('rownames(countsTable) <- data$'+title_of_names_column+'\n')
	conditions = eval(conditions)
        line = '", "'.join(conditions)
        wf.write('conds <- c("'+line+'")\n')
        wf.write('cds <- newCountDataSet( countsTable, conds )\n')
        wf.write('cds <- estimateSizeFactors( cds )\n')
        wf.write('sizeFactors(cds)\n')
        wf.write("cds <- estimateDispersions( cds )\n")
        wf.write('res <- nbinomTest( cds, "'+condition1+'", "'+condition2+'" )\n')
        wf.write('plotDE <- function( res ) plot(res$baseMean,res$log2FoldChange, log="x", pch=20, cex=.1, col = ifelse( res$padj < .1, "red", "black" ) )\n')
        wf.write("png('"+graph_file+"')\n")
        wf.write('plotDE( res )\n')
        wf.write('dev.off()\n')
        wf.write('resSig <- res[ res$padj < .1, ]\n')
        wf.write('write.table(res, file = "'+outfileallinputs+r'", append = FALSE, sep = "\t")'+"\n")
        wf.write('write.table(resSig, file = "'+outfilesig+r'", append = FALSE, sep = "\t")'+"\n")
        wf.write('write.table(resSig[ order(resSig$pval), ], file = "'+outfilesig_orderpval+r'", append = FALSE, sep = "\t")'+"\n")
        wf.write('sink()\n')



def create_DESeqRscript_no_replicates(infile="/projects/dowellde/groseq/data/set1/clipped_fastqM10/samfiles/sortedbamfiles/lncRNAs/compare_cov_fileless1neg_istead.txt",column1=11, column2=14, type_transcript="lncRNAs", condition1="DMS0", condition2="Nutlin", title_of_names_column="name", order_flip="N"):
	"""infile should be a compare coverage file. Start counting from 1."""

	f = open(infile)
	headers = f.readline()
	headers = headers.strip("\n")
	headers = headers.split("\t")
	f.close()
	infile_dir = infile.split("/")[:-1]
	infile_dir = "/".join(infile_dir)+"/"
	infile_root = infile.split("/")[-1].strip(".txt")
	headercondition1 = headers[column1-1]#adjust for the fact python starts counting with 0 and R with 1
	headercondition2 = headers[column2-1]#adjust for the fact python starts counting with 0 and R with 1
	if order_flip=="N":
		outfile = infile_dir+infile_root+"."+headercondition1+headercondition2+type_transcript
	else:
		outfile = infile_dir+infile_root+"."+headercondition2+headercondition1+type_transcript
	write_file = outfile+".R"
	print write_file
	wf = open(write_file ,"w")
	R_dump_file = outfile+".Rout"
	graph_file = outfile+".png"
	outfileallinputs = outfile+".res.txt"
	outfilesig = outfile+".resSig.txt"
	outfilesig_orderpval = outfile+".resSig_pvalue.txt"
	wf.write('sink("'+R_dump_file+'")\n')
	wf.write('library( DESeq )\n')
	wf.write('data <- read.delim("'+infile+r'", sep="\t", header=TRUE)'+"\n")#need to check that \t comes out like it should. Might write it wrong.
	wf.write('countsTable <- subset(data, select=c('+str(column1-1)+','+str(column2-1)+'))\n')
	wf.write('rownames(countsTable) <- data$'+title_of_names_column+'\n')
	wf.write('conds <- c("'+condition1+'", "'+condition2+'")\n')
	wf.write('cds <- newCountDataSet( countsTable, conds )\n')
	wf.write('cds <- estimateSizeFactors( cds )\n')
	wf.write('sizeFactors(cds)\n')
	wf.write("cds <- estimateDispersions( cds, method='blind', sharingMode='fit-only' )\n")
	if order_flip=="N":
		wf.write('res <- nbinomTest( cds, "'+condition1+'", "'+condition2+'" )\n')
	else:
		wf.write('res <- nbinomTest( cds, "'+condition2+'", "'+condition1+'" )\n')
	wf.write('plotDE <- function( res ) plot(res$baseMean,res$log2FoldChange, log="x", pch=20, cex=.1, col = ifelse( res$padj < .1, "red", "black" ) )\n')
	wf.write("png('"+graph_file+"')\n")
	wf.write('plotDE( res )\n')
	wf.write('dev.off()\n')
	wf.write('resSig <- res[ res$padj < .1, ]\n')
	wf.write('write.table(res, file = "'+outfileallinputs+r'", append = FALSE, sep = "\t")'+"\n")
	wf.write('write.table(resSig, file = "'+outfilesig+r'", append = FALSE, sep = "\t")'+"\n")
	wf.write('write.table(resSig[ order(resSig$pval), ], file = "'+outfilesig_orderpval+r'", append = FALSE, sep = "\t")'+"\n")
	wf.write('sink()\n')


def up_and_down(resSig_pvalue_file, split_names=1):
	f = open(resSig_pvalue_file)
	title = f.readline()
	wf1 = open(resSig_pvalue_file+".upname.txt", "w")
	wf2 = open(resSig_pvalue_file+".downname.txt", "w")
	for line in f:
		line = line.strip("\n")
		line = line.split("\t")
		name = line[1]
		if not name.startswith("NA"):
			if split_names==1:
				rootname = name.split(";")[1]
				rootname = rootname.split("=")[1]
			else:
				rootname = name
			fold_change = float(line[5])
			info = [rootname, str(fold_change)]
			if fold_change>1:
				wf1.write("\t".join(map(str,info))+"\n")
			else:
				wf2.write("\t".join(map(str,info))+"\n")
	f.close()
	wf1.close()
	wf2.close()

def gro_names(resfile="/projects/dowellLab/groseq/data/multiple_sets/remap/samfiles/sortedbamfiles/forpaper101712/doubles/gffcoverage/compare_cov_fileless1neg_instead.txt.protein_coding.DMSO2_3_depthNutlin2_3_depthgff.res.txt", split_names=1):
	f = open(resfile)
	title = f.readline()
	wf = open(resfile+".names.txt", "w")
	for line in f:
		line = line.strip("\n")
		line = line.split("\t")
		name = line[1]
		if not name.startswith("NA"):
			if split_names==1:
				rootname = name.split(";")[1]
				rootname = rootname.split("=")[1]
			else:
				rootname = name
			wf.write(rootname+"\n")
	f.close()
	wf.close()



if __name__=="__main__":
	if len(sys.argv)<3:
		print "python DESeq_Rscript.py <infile> <column1> <column2> <type_transcript> <condition1> <condition2>"
		print "column number as R would use it. So start counting at 1"
	else:
		infile=sys.argv[1]
		column1=int(sys.argv[2])
		column2=int(sys.argv[3])
		type_transcript=sys.argv[4]
		condition1=sys.argv[5]
		condition2=sys.argv[6]
		if len(sys.argv)>7:
			order_flip = sys.argv[7]
			create_DESeqRscript_no_replicates(infile=infile,column1=column1, column2=column2, type_transcript=type_transcript, condition1=condition1, condition2=condition2, order_flip=order_flip)
		else:
			create_DESeqRscript_no_replicates(infile=infile,column1=column1, column2=column2, type_transcript=type_transcript, condition1=condition1, condition2=condition2)	
