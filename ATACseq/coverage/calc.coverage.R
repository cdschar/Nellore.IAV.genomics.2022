#version 3.1

#R functions for processing data
source("bisTools.R");

#set name for output files
experimentName = "TIV.ATACseq"

#set genome
genome = "hg38"

#set working directory as coverage subfolder of the analysis folder
homeDir = "~/coverage/";
setwd(homeDir);

#read in the manifest file
filesDir = "~/pipeline/";
filesFile = "ATACseq.sample.manifest.txt";
files = read.table(paste0(filesDir, filesFile), sep = "\t", header = TRUE, as.is = T);
files = files[files$include,]

#chromosomes to exclude in analysis
exChr = c("chrM|chrY|random");

##########
# make union peak file from all samples with HOMER

mergePeakFiles = paste0(files$PeakFile, collapse = " ");
mergeFile = paste0(homeDir, experimentName, ".txt");
mergeCall = paste("mergePeaks", mergePeakFiles, ">", mergeFile);
system(mergeCall)

#convert to bed/remove peak file
bedFile = paste0(homeDir, experimentName, ".bed");
bedCall = paste("pos2bed.pl", mergeFile, ">", bedFile);
system(bedCall);
system(paste("rm", mergeFile));

###########
# calculate coverage from bam files

#output files
ctsFile = gsub(".bed", ".cts.txt", bedFile)
rpmFile = gsub(".bed", ".rpm.txt", bedFile)
rppmFile = gsub(".bed", ".rppm.txt", bedFile)
fripBarFile = gsub(".bed", ".FRiP.barplot.pdf", bedFile)

#load mergePeaks bed file into Granges object
bed = import.bed(paste0(bedFile), genome=genome)

#exclude ENCODE black listed regions
blacklist = import.bed( paste0("~/seqTools/genomes/", genome, "/", genome, ".blacklist.bed"), genome=genome)
bed = setdiff(bed, blacklist, ignore.strand=T)
bedFile2 = paste0(homeDir, experimentName, ".excluded.bed");
export.bed(bed, bedFile2)

#count reads in each peak for each file
cts = data.frame(peak = paste(seqnames(bed), start(bed), end(bed), sep = "_"), chr = paste(seqnames(bed)), start = paste(start(bed)), end = paste(end(bed)), strand = "+")

#remove peaks from certain chromosomes
cts = cts[!grepl(exChr, cts$peak), ]
cts = data.table(cts)

for (i in 1:dim(files)[1]) {

	print(paste(files$sample[i], Sys.time()))

	sbp = ScanBamParam(flag = scanBamFlag(isDuplicate = F), which = bed) 
	reads = readGAlignmentPairs(paste0(files$dir[i], files$bamFile[i]), param = sbp)

	#Get counts
	so = summarizeOverlaps(bed, reads, mode = "Union", singleEnd = F, ignore.strand = T)

	cts[match(paste(seqnames(rowRanges(so)), start(rowRanges(so)), end(rowRanges(so)), sep = "_"), cts$peak), as.character(files$sample[i])] = assays(so)$counts

}

write.table(cts, file = ctsFile, sep = "\t", quote = F, row.names = F)

###############
#annotate peaks
tempAnnFile = paste0(homeDir, "homer.annotated.txt");
homerAnn = paste("annotatePeaks.pl", ctsFile, genome, ">", tempAnnFile)
system(homerAnn);

#load annotated file
ann = read.table(file = tempAnnFile, sep = "\t", header = T, quote = "", comment.char = "")
colnames(ann)[1] = "peak"
#remove extra info
limit = c("Type|Nearest|Description|Score|Size|Chr|Strand|Start|End")
ann = subset(ann, , !grepl(limit, colnames(ann)))

#remove files
system(paste("rm", tempAnnFile))

################
#perform normalization

#make RPM normalized counts file and annotate
rpm = cts
rpm = subset(rpm, , files$sample)
colnames(rpm) = paste0(files$sample, ".rpm")
rpm = t(t(rpm) * 1e6 / files$unique.reads)
rpm = round(rpm, 3)
rpm = cbind(subset(cts, , !grepl(paste(files$sample, collapse = "|"), colnames(cts))), rpm)
rpm = merge(rpm, ann);
rpm = rpm[order(match(rpm$peak,cts$peak)),]
write.table(rpm, file = rpmFile, sep = "\t", quote = F, row.names = F);


#calculate the fraction of reads in peaks
files$frip = colSums(subset(cts, , files$sample)) / files$unique.reads
outFile = "ATACseq.sample.manifest.withFRiP.txt";
write.table(files, file = paste0(filesDir, outFile), sep = "\t", row.names = F, quote = F);

#plot fraction of reads in peaks
cairo_pdf(file = fripBarFile, height = 5, width = 5);
par(mai = c(1.5, 0.75, 0.75, 0.25), mgp = c(2, 0.5, 0), family = "Arial");
barplot(files$frip*100, las = 3, main = "FRiP Scores", ylab= "Fraction of reads in peaks (in percent)", names.arg = files$sample)
abline(h=5, lty = 2, col = "red");
dev.off();

############
#optional
#normalize for reads per peak million (RPPM) = RPM / FRIP and annotate
rppm = cts
rppm = subset(rppm, , files$sample)
colnames(rppm) = paste0(files$sample, ".rppm")
rppm = t(t(rppm) * 1e6 / (files$unique.reads * files$frip))
rppm = round(rppm, 3)
#Add peak column back on
rppm = cbind(subset(cts, , !grepl(paste(files$sample, collapse = "|"), colnames(cts))), rppm)
rppm = merge(rppm, ann);
rppm = rppm[order(match(rppm$peak,cts$peak)),]
write.table(rppm, file = rppmFile, sep = "\t", quote = F, row.names = F);

#############
#filter peaks by a group RPM cutoff
cutoff = 1
rpm = as.data.frame(rpm)
detected = apply(rpm[,grepl(paste(files$sample, collapse = "|"), colnames(rpm))], 1, function(x) filter_detected(x, files$group, cutoff) )

cts_detected = cts[detected,]
write.table(cts_detected, file = paste0(experimentName, ".cts.detected.RPM.", cutoff, ".exon.txt"), sep = "\t", quote = F, row.names = F)

rpm_detected = rpm[detected,]
write.table(rpm_detected, file = paste0(experimentName, ".rpm.detected.RPM.", cutoff, ".exon.txt"), sep = "\t", quote = F, row.names = F);

rppm_detected = rppm[detected,]
write.table(rppm_detected, file = paste0(experimentName, ".rppm.detected.RPM.", cutoff, ".exon.txt"), sep = "\t", quote = F, row.names = F);
