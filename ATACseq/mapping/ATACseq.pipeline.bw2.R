#version 3.1

#R functions for processing data
source("bisTools.R");

#set pipeline directory
homeDir = "~/pipeline/";
setwd(homeDir);

#set the output directory for all the samples
baseDir = "~/data/"

#set genome specific options, choose one from hg19, hg38, mm9, mm10
bowtieGenome = select_ref_bw2( "hg38" )
macsGenome = select_ref_macs( "hg38" )

#chromosomes to exclude
exChr = c("chrM|chrY");

#read in the manifest file
filesDir = homeDir;
filesFile = "ATACseq.sample.manifest.txt";
files = read.table(paste0(filesDir, filesFile), sep = "\t", header = TRUE, as.is = T);

#mapping options
bowtieCmd = "bowtie2"; #bowtie2 v2.2.4
processes = 12;
bowtieOptionsPE = paste("-k 1 -X2000 -t --quiet --mm --threads", processes);
bowtieOptionsSE = paste("-k 1 -X2000 -t --quiet --mm --threads", processes);

#HOMER options #v4.8.2
tagCmd = "makeTagDirectory";
tagOptions = "-format sam";

#other tools
samCmd = "samtools";
picardCmd = "java -Xmx16G -jar ~/seqTools/picard/dist/picard.jar MarkDuplicates ";
macsCMD = "macs2";

bamSortExt = ".sort"

#adapter trimming options, uncomment one option
seq_platform = "nextera"  #Nextera Tn5 adapters
#seq_platform = "illumina" #TRUSEQ adapters

#go through each sample
for (i in 1:dim(files)[1]) {	
	
	#decide if sample should be processed
	if (files$include[i]) {

		print(paste(files$sample[i], Sys.time()))
		start = Sys.time();

		#move original fastq files to the sample directory
		mvFiles(files$fqMate1[i], files$dir[i], paste0(baseDir, files$sample[i], "/"))
		files$dir[i] = mvFiles(files$fqMate2[i], files$dir[i], paste0(baseDir, files$sample[i], "/"))
	
		#format and concatenate fastq files, keep original fastq names
		files$fqMate1_input[i] = files$fqMate1[i]
		files$fqMate2_input[i] = files$fqMate2[i]
		files$fqMate1[i] = formatFastq(files$fqMate1[i], files$dir[i], paste0(files$sample[i], "_1"), threads=processes)
		files$fqMate2[i] = formatFastq(files$fqMate2[i], files$dir[i], paste0(files$sample[i], "_2"), threads=processes)

		#cut 3 prime adapter sequences
		pelist = fqPECutadapt(files$fqMate1[i], files$fqMate2[i], files$dir[i], seq_platform, threads=processes)
		files$fqMate1[i] = pelist[[1]]
		files$fqMate2[i] = pelist[[2]]

		#run fastqc
		system(paste0("fastqc ", files$dir[i],files$fqMate1[i], " -o ", files$dir[i]))
		system(paste0("fastqc ", files$dir[i],files$fqMate2[i], " -o ", files$dir[i]))
		
		#single-end or paired-end bowtie call
		mapFile = paste0(files$dir[i], files$sample[i], ".sam")
		
		if (!is.null(files$fqFile[i]) && !is.na(files$fqFile[i]) && files$fqFile[i] != "") {
			bowtieCall = paste(bowtieCmd, bowtieOptionsSE, "-x", bowtieGenome, paste0(files$dir[i], files$fqFile[i]), "-S", mapFile)
			system(bowtieCall)
			files$fqFile[i] = gzipFastq(paste0(files$dir[i], files$fqFile[i]), threads=processes);
		
		} else if (!is.null(files$fqMate1[i]) && !is.na(files$fqMate1[i]) && files$fqMate1[i] != "" && !is.null(files$fqMate2[i]) && !is.na(files$fqMate2[i]) && files$fqMate2[i] != "") {
			
			bowtieCall = paste(bowtieCmd, bowtieOptionsPE, "-x", bowtieGenome, "-1", paste0(files$dir[i], files$fqMate1[i]), "-2", paste0(files$dir[i], files$fqMate2[i]), "-S", mapFile)
			system(bowtieCall)
			files$fqMate1[i] = gzipFastq(paste0(files$dir[i], files$fqMate1[i]), threads=processes);
			files$fqMate2[i] = gzipFastq(paste0(files$dir[i], files$fqMate2[i]), threads=processes);
		}
		
			
		#Convert sam to bam
		files$bamFile[i] = paste0(files$dir[i], files$sample[i], ".bam")
		system(paste(samCmd, "view -bS", mapFile, ">", files$bamFile[i]))
		
		#sort bam
		files$bamFile[i] = sortBam(files$bamFile[i], bamSortFile = paste0(files$dir[i], files$sample[i], bamSortExt), delBam = T, threads=processes)

		#mark duplicates
		files$bamFile[i] = markDups(paste0(files$dir[i], files$bamFile[i]), picardCmd, delBam = T)

		#make BigWig
		bamToBigWig(paste0(files$dir[i], files$bamFile[i]), isDup = FALSE, removeChrs = exChr)
		
		#call peaks with MACS2
		macsOutDir = paste0(files$dir[i], files$sample[i], ".MACS2.Peaks/");
		files$PeakFile[i] = paste0(macsOutDir, files$sample[i], "_peaks.narrowPeak");
		callMACS = paste(macsCMD, "callpeak -t", paste0(files$dir[i], files$bamFile[i]), "-f BAM -g", macsGenome, "-n", files$sample[i], "--outdir", macsOutDir);
		system(callMACS);
		
		#get mapping stats	
		cts = getBamCts(paste0(files$dir[i], files$bamFile[i]))	

		files$unmapped.reads[i] = cts[1];
		files$mapped.reads[i] = cts[2];
		files$unique.reads[i] = cts[3];
		files$paired.reads[i] = cts[4];
		
		#get peaks count
		if( file.exists(paste0(macsOutDir, files$sample[i], "_summits.bed")) ){
			peaks = read.table( paste0(macsOutDir, files$sample[i], "_summits.bed") )		
			files$MACS.peaks[i] = dim(peaks)[1]
		} else{
			files$MACS.peaks[i] = 0
		}
		
		#remove any remaining sam files
		system(paste0("rm ", baseDir, files$sample[i], "/*.sam"))

		#write the time taken to process
		end = Sys.time();
		files$timeTaken[i] = difftime(end, start, units= "hours");
	
		#update sample manifest file
		write.table(files, file = paste0(filesDir, filesFile), sep = "\t", row.names = F, quote = F);
	}
}

###
# software versions
###

#bowtie v1.1.1
#HOMER v4.8.2
#picard 1.127(6fd825a0b7c45834d5cc2eb92983322faa72e4d3_1418936317)
#samtools v1.5
#macs2 2.1.0.20140616

#R version 3.2.1 (2015-06-18)
#Platform: x86_64-apple-darwin13.4.0 (64-bit)
#Running under: OS X 10.10.5 (Yosemite)

#locale:
#[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

#attached base packages:
#[1] stats4    parallel  stats     graphics  grDevices utils     datasets 
#[8] methods   base     

#other attached packages:
# [1] BSgenome.Hsapiens.UCSC.hg19_1.4.0 data.table_1.9.6                 
# [3] chipseq_1.20.0                    ShortRead_1.28.0                 
# [5] GenomicAlignments_1.6.3           SummarizedExperiment_1.0.2       
# [7] Biobase_2.30.0                    BiocParallel_1.4.3               
# [9] Rsamtools_1.22.0                  BSgenome_1.38.0                  
#[11] rtracklayer_1.30.3                GenomicRanges_1.22.4             
#[13] GenomeInfoDb_1.6.3                Biostrings_2.38.4                
#[15] XVector_0.10.0                    IRanges_2.4.8                    
#[17] S4Vectors_0.8.11                  BiocGenerics_0.16.1              

#loaded via a namespace (and not attached):
# [1] zlibbioc_1.16.0      lattice_0.20-33      hwriter_1.3.2       
# [4] tools_3.2.1          grid_3.2.1           latticeExtra_0.6-28 
# [7] lambda.r_1.1.7       futile.logger_1.4.1  RColorBrewer_1.1-2  
#[10] futile.options_1.0.0 bitops_1.0-6         RCurl_1.95-4.8      
#[13] XML_3.98-1.4         chron_2.3-47        
