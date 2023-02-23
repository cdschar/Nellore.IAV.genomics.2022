#map and QC fastq files
source("bisTools.R");
homeDir = "~/pipeline/";
setwd(homeDir);

fqDir = homeDir;
fqFile = "RNAseq.sample.manifest.txt";
files = read.table(paste0(fqDir, fqFile), sep = "\t", header = T, as.is = T);
files = files[which(files$include == TRUE), ]

bowtieGenome = "~/bowtie/bowtie2.indices/hg19.noHap.ERCC/hg19.noHap.ERCC"

gtfFile = "~/seqTools/genomes/hg19/iGenomes.UCSC.hg19/genes.gtf";
tophatCmd = paste("tophat2 -p 8 -N 2 --max-multihits 1 --read-gap-length 1 --mate-std-dev 500 -G", gtfFile)

picardCmd = "java -jar ~/seqTools/picard/dist/picard.jar MarkDuplicates "

tophatBamFile = "accepted_hits.bam";
unMappedBamFile = "unmapped.bam"

baseDir = "~/data/"

bamSortExt = ".sort"

for (i in 1:dim(files)[1]) {	
	if (files$include[i]) {

		print(paste(files$sample[i], Sys.time()))
		start = Sys.time();

		#format and move fastq files
		mvFiles(files$fqMate1[i], files$dir[i], paste0(baseDir, files$sample[i], "/"))
		files$dir[i] = mvFiles(files$fqMate2[i], files$dir[i], paste0(baseDir, files$sample[i], "/"))
	
		#format and concatenate fastq files
		files$fqFile[i] = formatFastq(files$fqFile[i], files$dir[i], files$sample[i])
		files$fqMate1[i] = formatFastq(files$fqMate1[i], files$dir[i], paste0(files$sample[i], "_1"))
		files$fqMate2[i] = formatFastq(files$fqMate2[i], files$dir[i], paste0(files$sample[i], "_2"))

		#cut Nextera adapter sequences
		pelist = fqPECutadapt(files$fqMate1[i], files$fqMate2[i], files$dir[i])
		files$fqMate1[i] = pelist[[1]]
		files$fqMate2[i] = pelist[[2]]
		
		#system(paste0("fastqc ", files$dir[i],files$fqFile[i], " -o ", files$dir[i]))
		system(paste0("fastqc ", files$dir[i],files$fqMate1[i], " -o ", files$dir[i]))
		system(paste0("fastqc ", files$dir[i],files$fqMate2[i], " -o ", files$dir[i]))

		#single-end or paired-end tophat call
		if (!is.null(files$fqFile[i]) && !is.na(files$fqFile[i]) && files$fqFile[i] != "") {
			tophatCall = paste(tophatCmd, "-o", files$dir[i], bowtieGenome, paste0(files$dir[i], files$fqFile[i]))
			system(tophatCall)
			files$fqFile[i] = gzipFastq(paste0(files$dir[i], files$fqFile[i]));
		} else if (!is.null(files$fqMate1[i]) && !is.na(files$fqMate1[i]) && files$fqMate1[i] != "" && !is.null(files$fqMate2[i]) && !is.na(files$fqMate2[i]) && files$fqMate2[i] != "") {
			tophatCall = paste(tophatCmd, "-o", files$dir[i], bowtieGenome, paste0(files$dir[i], files$fqMate1[i], " ", files$dir[i], files$fqMate2[i]))
			system(tophatCall)
			files$fqMate1[i] = gzipFastq(paste0(files$dir[i], files$fqMate1[i]));
			files$fqMate2[i] = gzipFastq(paste0(files$dir[i], files$fqMate2[i]));
		}
	
		#sort bam
		files$bamFile[i] = sortBam(paste0(files$dir[i], tophatBamFile), bamSortFile = paste0(files$dir[i], files$sample[i], bamSortExt), delBam = T)

		#mark duplicates
		files$bamFile[i] = markDups(paste0(files$dir[i], files$bamFile[i]), picardCmd, delBam = T)

		#get counts	
		cts = getBamCts(paste0(files$dir[i], files$bamFile[i]), unMappedBamFile = paste0(files$dir[i], unMappedBamFile))

		files$unmapped.reads[i] = cts[1];
		files$mapped.reads[i] = cts[2];
		files$unique.reads[i] = cts[3];
		files$paired.reads[i] = cts[4];

		system(paste0("rm ", baseDir, files$sample[i], "/","flow_cell*gz"))
	
		end = Sys.time();
		files$timeTaken[i] = difftime(end, start, units= "hours");
		
		#update files manifest
		write.table(files, file = paste0(fqDir, fqFile), sep = "\t", row.names = F, quote = F);
	}
}


###########software versions #######
#bowtie v2.2.4
#TopHat v2.0.13
#FastQC v0.11.4
#R version 3.2.3

#samtools 1.5
#Using htslib 1.5



