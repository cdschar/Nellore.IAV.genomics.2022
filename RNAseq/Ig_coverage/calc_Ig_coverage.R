#calculate Ig segment coverage

#set working directory
homeDir = "~/analysis/Ig/";
setwd(homeDir);

#load bistools suite
source("bisTools.R");

#set experiment name
experimentName = "TIV.Ig.Coverage"

#load manifest files
fqDir = "~/pipeline/";
fqFile = "RNAseq.sample.manifest.txt";
files = read.table(paste0(fqDir, fqFile), sep = "\t", header = TRUE);
files = files[files$include, ]

#set species file from IMGT
Ig = "IG.hg38.bed";

#load Ig annotation bed file into Granges object
bed = import.bed(Ig)

#set file names
ctsFile = paste0(experimentName, ".cts.txt")
rpmFile = paste0(experimentName, ".rpm.txt")
fpkmFile = paste0(experimentName, ".fpkm.txt");

#count reads in each Ig region for each sample
cts = data.frame(Segment = mcols(bed)[,"name"], loc = paste(seqnames(bed), start(bed), end(bed), sep = "_"), width = width(bed))

for (i in 1:dim(files)[1]) {

	print(paste(files$sample[i], Sys.time()))

	sbp = ScanBamParam(flag = scanBamFlag(isDuplicate = F), which = bed) 
	reads = readGAlignmentPairs(paste0(files$dir[i], files$bamFile[i]), param = sbp)

	#Get counts
	so = summarizeOverlaps(bed, reads, mode = "Union", singleEnd = F, ignore.strand = T)

	cts[match(paste(seqnames(rowRanges(so)), start(rowRanges(so)), end(rowRanges(so)), sep = "_"), cts$loc), as.character(files$sample[i])] = assays(so)$counts

}

#write raw counts file
write.table(cts, file = ctsFile, sep = "\t", quote = F, row.names = F)

#make RPM normalized counts file
geneRpm = cts
geneRpm = geneRpm[, grepl(paste(files$sample, collapse = "|"), colnames(geneRpm))]
geneRpm = t(t(geneRpm) * 1e6 / files$unique.reads)
geneRpm = cbind(cts[, !grepl(paste(files$sample, collapse = "|"), colnames(cts))], geneRpm)
write.table(geneRpm, file = rpmFile, sep = "\t", quote = F, row.names = F);

#make RPKM normalized counts file
geneRpkm = geneRpm
geneRpkm = geneRpkm[, grepl(paste(files$sample, collapse = "|"), colnames(geneRpkm))]
geneRpkm = (geneRpkm * 1000)/ cts$width
geneRpkm = cbind(cts[, !grepl(paste(files$sample, collapse = "|"), colnames(cts))], geneRpkm)
write.table(geneRpkm, file = fpkmFile, sep = "\t", quote = F, row.names = F);

