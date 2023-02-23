library("Biostrings");
library("BSgenome");
library("rtracklayer");
library("Rsamtools");
library("ShortRead");
library("chipseq")
library("data.table")
library("GenomeInfoDb")


# global constants
gSigDigits = 5;
gNormReads = 1e6

## bamToBigWig constants
gFragSize = NA;
gBinSize = 300;
gBwFile = "gsub('.bam', paste('.frag', paste(fragSizes, collapse = '_'), '.norm', format(normReads, scientific = FALSE), '.bw', sep = ''), bamFiles)";

gGenome = "BSgenome.Mmusculus.UCSC.mm9"

#### plotting constants
gMainSize = 1.5;
gLabSize = 1;
gAxisSize = 0.8;
gPch = 19;
gCex = 0.25
gHeight = 3;
gWidth = 3;

gSdLwd = 0.5;
gSdCol = rgb(0.5, 0.5, 0.5);

gRegLwd = 3;
gRegCol = rgb(0, 0, 0.5);
gRegCex = 2;

gMai = c(0.5, 0.5, 0.5, 0.25);
gMgp = c(1.5, 0.5, 0); 


#create dir and move files
mvFiles = function(files, dir, mvDir, split = c("\\|")) {
  
  #Debug	files = paste0(files$dir[i], strsplit(files$fqFile[i], "\\|")[[1]]); dir = paste0(files$dir[i], files$sample[i], "/");
  
  if (!is.null(files) && !is.na(files) && files != "") {	
    
    files = paste0(dir, strsplit(files, split)[[1]])
    
    if (!file.exists(mvDir)) dir.create(mvDir); #create new directory if necessary
    file.copy(files, paste0(mvDir, basename(files)))
    
    return(mvDir)
  }
}

formatFastq = function(fqFiles, dir, sampleName = c(NA), split = c("\\|"), threads=c(1)) {
  
  #Debug	fqFiles = files$fqMate1[i]; dir = files$dir[i]; sampleName = paste0(files$sample[i], "_1"); split = "\\|";
  #Debug  fqFiles = files$fqFile[i]; dir = files$dir[i]; sampleName = files$sample[i]; split = "\\|"; 
  
  if (!is.null(fqFiles) && !is.na(fqFiles) && fqFiles != "") {	
    
    fqFiles = paste0(dir, strsplit(fqFiles, split)[[1]])
    
    gzExt = ".gz";
    fqExt = ".fastq";
    
    print(paste("Formatting fastq files", Sys.time()))
    
    #if multiple fastq files are provided decompress and combine for ...
    if (length(fqFiles) > 1) {
      
      if (is.na(sampleName)) sampleName = gsub(fqExt, "", gsub(gzExt, "", basename(fqFiles[1])))
      
      fqFile = paste0(dir, sampleName, fqExt); #combined fastq File name
      
      for (i in 1:length(fqFiles)) if (grepl(paste0(gzExt, "$"), fqFiles[i])) system(paste("pigz -p", threads, "-cd", fqFiles[i], ">", gsub(gzExt, "", fqFiles[i]))); #decompress files
      
      system(paste("cat", paste(gsub(gzExt, "", fqFiles), collapse = " "), ">", fqFile))
      
      for (i in 1:length(fqFiles)) if (grepl(paste0(gzExt, "$"), fqFiles[i])) file.remove(gsub(gzExt, "", fqFiles[i])); #delete files after concatenated
      
    } else {
      
      if (grepl(paste0(gzExt, "$"), fqFiles)) {
        system(paste("pigz -p", threads, "-cd", fqFiles, ">", gsub(gzExt, "", fqFiles))); #decompress files
        fqFile = gsub(gzExt, "", fqFiles)
      } else fqFile = fqFiles;
    }
    return(as.character(basename(fqFile)))
    
  } else return(NULL)
  
}

#Function to cut adapter seqs for paired end reads
fqPECutadapt = function(fqMate1, fqMate2, dir, platform="nextera"){
  
  if( platform == "nextera"){
    adapter_seq = "CTGTCTCTTATA"
  } else if( platform == "illumina"){ 
    adapter_seq = "AGATCGGAAGAGC"
  } else if( platform == "small_rna"){
    adapter_seq = "TGGAATTCTCGG"
  }
  
  if( grepl(".gz", fqMate1) ){
    fqname = gsub("_1.fastq.gz", "", fqMate1)
    fout1 = paste0(fqname, "-trimmed-pair1.fastq.gz")
    fout2 = paste0(fqname, "-trimmed-pair2.fastq.gz")
  } else {
    fqname = gsub("_1.fastq", "", fqMate1)
    fout1 = paste0(fqname, "-trimmed-pair1.fastq")
    fout2 = paste0(fqname, "-trimmed-pair2.fastq")
  }
  
  system(paste0("skewer -x ", adapter_seq, " -q 3 -o ", dir,fqname, " ", dir,fqMate1, " ", dir,fqMate2))
  system(paste0("rm ", dir, fqMate1, " ", dir, fqMate2))
  return(list(fout1, fout2))
  
}

#check to see if fastq file is gzipped and if not compress and return name
gzipFastq = function(fqFile,  overwrite = c(F), rmIfExists = c(F), threads=8) {
  
  #Debug	fqFile = paste0(files$dir[i], files$fqMate1[i]);
  gz.ext = ".gz";
  fastq.ext = ".fastq"
  
  if (grepl(paste0(fastq.ext, "$"), fqFile) & (!file.exists(paste0(fqFile, gz.ext)) | overwrite)) { 
    system(paste("pigz -p", threads, fqFile))
    return(paste0(basename(fqFile), gz.ext))
  } else if (grepl(paste0(fastq.ext, "$"), fqFile) & file.exists(paste0(fqFile, gz.ext)) & rmIfExists) { 
    file.remove(fqFile)
    return(paste0(basename(fqFile), gz.ext))
  } else return(basename(fqFile))
}

#Sort bam file
sortBam = function(bamFile, bamSortFile = c(NA), delBam = c(F), threads = c(1)) {
  
  print(paste("sorting Bam", basename(bamFile), Sys.time()))
  
  idxExt = ".bai";
  bamExt = ".bam";
  bamSortExt = ".sort";
  
  sortCmd = "samtools sort "
  
  if (is.na(bamSortFile)) bamSortFile = gsub(bamExt, bamSortExt, bamFile)
  
  #sort BAM file
  system(paste0(sortCmd, "-o ", bamSortFile, ".bam ", bamFile, " -@ ", threads))
  
  #Remove unsorted bam file
  if (delBam) file.remove(bamFile)
  if (delBam & file.exists(paste0(bamFile, idxExt))) file.remove(paste0(bamFile, idxExt))
  
  #Create BAM index for fast access
  indexBam(paste0(bamSortFile, bamExt))		
  
  return(basename(paste0(bamSortFile, bamExt)))
  
}

#Use picard to mark duplicates in bam file 
markDups = function(bamFile, picardCmd, delBam = c(F)) {
  
  print(paste("Marking Duplicates", bamFile, Sys.time()))
  
  idxExt = ".bai";
  bamExt = ".bam";
  bamDupMarkExt = ".dupMark.bam";
  
  bamDupFile = gsub(paste0(bamExt, "$"), bamDupMarkExt, bamFile)
  
  #Mark duplicates using picard
  system(paste0(picardCmd, "INPUT=", bamFile, " OUTPUT=", bamDupFile, " METRICS_FILE=", bamDupFile, ".metrics"))
  
  #Create BAM index for fast access
  indexBam(bamDupFile)		
  
  #Remove BAM without duplicates marked
  if (delBam) file.remove(bamFile)
  if (delBam & file.exists(paste0(bamFile, idxExt))) file.remove(paste0(bamFile, idxExt))
  
  return(basename(bamDupFile))
}

#Get read counts for bam file
getBamCts = function(bamFile, unMappedBamFile = c(NA)) {
  
  print(paste("Getting Bam File Read Counts", basename(bamFile), Sys.time()))
  
  baiExt = ".bai";
  
  #create index for bamFile if it doesn't exist
  if (!file.exists(paste0(bamFile, baiExt))) indexBam(bamFile); 
  if (!is.na(unMappedBamFile)) {
    if (!file.exists(paste0(unMappedBamFile, baiExt))) indexBam(unMappedBamFile);
    unmappedCt = countBam(unMappedBamFile)$records;
  } else {
    sbp = ScanBamParam(flag = scanBamFlag(isUnmappedQuery = T))
    unmappedCt = countBam(bamFile, param = sbp)$records;			
  }
  
  sbp = ScanBamParam(flag = scanBamFlag(isUnmappedQuery = F))
  mappedCt = countBam(bamFile, param = sbp)$records;
  
  sbp = ScanBamParam(flag = scanBamFlag(isUnmappedQuery = F, isDuplicate = F))	
  uniqueCt = countBam(bamFile, param = sbp)$records;
  
  sbp = ScanBamParam(flag = scanBamFlag(isUnmappedQuery = F, isDuplicate = F, isPaired = T))	
  pairedCt = countBam(bamFile, param = sbp)$records;
  
  sbp = ScanBamParam(flag = scanBamFlag(isUnmappedQuery = F, isDuplicate = F, isPaired = T, isProperPair = T))	
  properPairedCt = countBam(bamFile, param = sbp)$records;
  
  return(c(unmappedCt, mappedCt, uniqueCt, pairedCt, properPairedCt))
  
}

#bamToBigWig function: converts a bam file into a bigWig coverage file.  Removes chromosomes that match the variable removeChrs by grepl.  Normalizes for read count .Only reads not in removeChrs are counted. The normReads variable sets the sequencing coverage. The default normalization is reads per million (e.g. normRead = 1e6).  
bamToBigWig = function(bamFiles, bwFile = c(gBwFile), fragSizes = c(gFragSize), sigDigits = c(gSigDigits), normReads = c(gNormReads), extCall = c(T), flag = scanBamFlag(), removeChrs = c(NA), ...) {
  
  #Debug
  #bamFiles = bamFile; bwFile = "test.bw"; fragSizes = gFragSize; sigDigits = gSigDigits; normReads = gNormReads; extCall = T; flag = scanBamFlag(); removeChrs = c("random");
  
  #Get genome version and seqinfo
  si = list(); for (i in 1:length(bamFiles)) si[[i]] = seqinfo(BamFile(bamFiles[i]))
  
  #if there are more bamFiles than fragSizes then recycle fragSizes
  if (!all(is.na(fragSizes)) & length(bamFiles) > length(fragSizes)) fragSizes = rep(fragSizes, ceiling(length(bamFiles) / length(fragSizes)));
  
  #Get readcounts for each bamFile
  #readCounts = list(); for (i in 1:length(bamFiles)) readCounts[[i]] = countBam(bamFiles[i], param = ScanBamParam(flag = flag))$records
  
  #Create bigwig and bedfile output file names	
  if (bwFile == gBwFile) bwFile = eval(parse(text = gBwFile))
  bgFile = paste0(bwFile, ".bedGraph")
  
  #get all chromosomes in bamFiles
  chrs = unique(unlist(lapply(si, seqnames)))
  
  #remove chromosomes that match removeChrs
  if (!is.na(removeChrs)) chrs = chrs[!grepl(removeChrs, chrs)] 
  
  #coverage and read count variables
  cv = RleList(compress = F); 
  ct = list(); for (bamFile in bamFiles) ct[[bamFile]] = 0
  
  #Recurse through each  chromosome
  for (chr in chrs) {
    print(paste(chr, Sys.time()))
    
    #Recursively build coverage for each bamFile 
    for (i in 1:length(bamFiles)) {
      
      #Check to make sure the bam file has reads on chromosome chr 
      if (chr %in% si[[i]]@seqnames) {
        
        #Set up chr specific bam files
        sbp = ScanBamParam(flag = flag, which = GRanges(seqnames = chr, ranges = IRanges(start = 0, end = si[[1]]@seqlengths[si[[1]]@seqnames == chr])))
        
        #Read in chromosome from each bam file
        bam = readGAlignments(as.character(bamFiles[i]), param = sbp);
        
        #Extend reads to fragment size
        if (!is.na(fragSizes[i])) bam = suppressWarnings(resize(granges(bam), width = fragSizes[i], fix = "start"));
        
        #cumulatively add up coverages and reads
        if (i == 1) cv[[chr]] = coverage(bam)[[chr]] else cv[[chr]] = cv[[chr]] + coverage(bam)[[chr]];
        
        #add up read counts
        ct[[bamFiles[i]]] = ct[[bamFiles[i]]] + length(bam)
        
        #clear memory
        rm(bam); gc();
      }
    }
  }
  
  #normalize coverage to reads
  for (chr in chrs) cv[[chr]] = round(cv[[chr]] * normReads / sum(unlist(ct)), sigDigits)
  
  #if external calll write to bedgraph then call Kent tools
  if(extCall) {
    export.bedGraph(cv, bgFile)
    writeSiChromSizes(si[[1]], "temp.chrom.size")
    system(paste("bedGraphToBigWig", bgFile, "temp.chrom.size", bwFile))
    file.remove(bgFile); file.remove("temp.chrom.size");
  } else {
    export.bw(rleCv, bwFile)
  }
}

