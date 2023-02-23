# calculate splicing for human hg38 Igh constant regions
#make custom txdb object using UCSC old annotation track that has Igh splice variants
library(SplicingGraphs)
#https://bioconductor.org/packages/release/bioc/html/SplicingGraphs.html
library(GenomicRanges)
library(GenomicFeatures)

#set working directory
homeDir = "~/analysis/spliceGraphs/"
setwd(homeDir)

#read in manifest files
fqDir = "~/pipeline/"
fqFile = "RNAseq.sample.manifest.txt"
files = read.table(paste0(fqDir, fqFile), sep = "\t", header = T)

#limit to comparisons of interest
files = files[files$day == 7, ]
files = files[which(files$group == "FcrL5pos_TETpos_IgDneg_MBC" | files$group == "FcrL5neg_TETpos_IgDneg_MBC"),]

#establish txdb object
#use ENSEMBLE database with alternative splicing for Igh genes
en = makeTxDbFromEnsembl(organism="Homo sapiens")
en

#Ensemble Igh gene annotations
IgSymb = c("IGHM", 
           "IGHG1", "IGHG2", "IGHG3", "IGHG4", 
           "IGHA1", "IGHA2", 
           "IGHE")
IgEnsemb = c("ENSG00000211899", 
             "ENSG00000211896", "ENSG00000211893", "ENSG00000211897", "ENSG00000211892", 
             "ENSG00000211895", "ENSG00000211890", 
             "ENSG00000211891")

#limit to chr14 to save time
isActiveSeq(en)[-match("14", names(isActiveSeq(en)))] <- FALSE
names(which(isActiveSeq(en)))

#create splicing graphs object
sg = SplicingGraphs(en)
sg

#limit to Ig genes
ig = sg[IgEnsemb]

#plot splice graphs for each Ig isotype
for(i in 1:length(IgEnsemb)){
  iso = IgEnsemb[i]
  file = paste0(IgSymb[i], ".pdf")

  pdf(file = file)
  plot(ig[IgEnsemb[i]])
  dev.off()
  }

#annotate reads to splice graph

#set bam file flags and region to read in
IgRegion = GRanges(seqnames="chr14", ranges=IRanges(start = 105579380, end = 105856761))
flag <- scanBamFlag(isSecondaryAlignment=FALSE,
                    isNotPassingQualityControls=FALSE,
                    isDuplicate=FALSE)
#set param 
param = ScanBamParam(flag = flag, which= IgRegion)

#loop through and calculate splicing
for(i in 1:dim(files)[1]){
  
  cat("Processing ", files$sample[i], " ... ", sep="")
  
  galp = readGAlignmentPairs(paste0(files$dir[i], files$bamFile[i]), use.names=TRUE, param=param)

  #change chromosome style to Ensembl
  seqlevelsStyle(galp) <- "Ensembl"

  #assign reads to splice graphs
  ig = assignReads(ig, galp, sample.name=files$sample[i])

  cat("OK\n")
}

#make data frame of counts  
ig_counts = countReads(ig)

#normalize counts by unique reads
ig_norm = ig_counts

for(i in 1:dim(files)[1]){
    
  norm = round((ig_counts[,2+i] * 1e6/files$unique.reads[i]), digits = 3)
  
  ig_norm[,2+i] = norm          
                
}

write.table(ig_norm, "SpliceGraph.Igh_rpm.txt", sep = "\t", row.names = F, quote = F)

#output R session
capture.output(sessionInfo(), file = paste0("Rsession.Info.", gsub("\\D", "", Sys.time()), ".txt"))

