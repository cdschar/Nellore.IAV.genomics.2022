library("edgeR");

#set working directory
homeDir = "~/analysis/diff/";
setwd(homeDir);

#set significance threshold
sigFDR = 0.05
sigFC = 1

#read in sample manifest
filesDir = "~/pipeline/";
filesFile = "ATACseq.sample.manifest.txt";
files = read.table(paste0(filesDir, filesFile), sep = "\t", header = TRUE, as.is = T);
#remove samples if necessary
files = files[files$include, ]

#read in raw and normalized data
covDir = "~/coverage/"

ctFile = "TIV.ATACseq.cts.detected.RPM.1.exon.txt";
cts = read.table(paste0(covDir, ctFile), sep = "\t", header = T, quote = "", comment.char = "")

rpmFile = "TIV.ATACseq.rppm.detected.RPM.1.exon.txt";
rpm = read.table(paste0(covDir, rpmFile), sep = "\t", header = T, quote = "", comment.char = "")

#reorder peaks in rpm file to match cts file
rpm = rpm[order(match(rpm$peak,cts$peak)),]

#establish all pairwise group combinations for comparisoins
grps = unique(files$group)
grpPairs = combn(grps, 2)

#establish edgeR DGE object
dge = DGEList(as.matrix(cts[, files$sample]), group = files$group)

#calculate dispersions/normalization factors
dge = calcNormFactors(dge)

#plotMeanVar(dge)
cairo_pdf(file = paste0("MDSplot.cellType.all.pdf"), height = 5, width = 5);
par(mai = c(1.5, 0.75, 0.75, 0.25), mgp = c(2, 0.5, 0), family = "Arial");
plotMDS.DGEList(dge, main = "MDS Plot for Count Data", labels = gsub("Sample\\_[0-9]+\\_", "", colnames(dge$counts)), cex = 0.3)
dev.off();

#######################################################################
## GLM Function

#set up differential matrix to append data to
diff = rpm

#set up design matrix
design = model.matrix(~0 + files$group + files$patient)
dge = estimateDisp(dge, design)
fit = glmFit(dge, design)
rownames(design) = colnames(dge)

#change design headers to simple strings
#this may have to be customized
colnames(design) = gsub("files\\$group", "", colnames(design))
colnames(design) = gsub("files\\$", "", colnames(design))

#loop through each comparison and append stats to diff matrix
for (i in 1:dim(grpPairs)[2]) {

	#set comparisons
	comp = paste0(grpPairs[2, i], ".v.", grpPairs[1, i]);
	print(paste(comp,  Sys.time()))
	
	#establish contrast matrix
	cont = makeContrasts(contrasts = paste0(grpPairs[2, i], "-", grpPairs[1, i]), levels = design)
	
	#create matrix for this comparisons and use glmLRT function for differential analysis
	lrt = glmLRT(fit, contrast = cont)
	
	#FDR correcting Pvalues
	lrt$table$fdr = p.adjust(lrt$table$PValue, method = "fdr")
  
  	#reduce number of digits of data
  	lrt$table$logFC = round(lrt$table$logFC, 3)
  	lrt$table$logCPM = round(lrt$table$logCPM, 3)
  
	#add header to the lrt$table slot that describes the comparison performed
	colnames(lrt$table) = paste0(comp, ".", names(lrt$table))
	
	#cbind diff test onto the diff matrix
	diff = cbind(diff, lrt$table)
	
	#Determine if a genes is significant by both fdr and log2FC > 1
	diff[[paste0(comp, ".sig")]] = (diff[[paste0(comp, ".fdr")]] < sigFDR & abs(diff[[paste0(comp, ".logFC")]]) >= sigFC)
}

#add column that tells if a gene is significant in any comparison
if( sum(grepl(".sig", colnames(diff)))>1 ){
	diff$sigAny <- apply(diff[,(grep(".sig"== "TRUE", diff))], 1, any)
} else{
	diff$sigAny <- diff[,grep(".sig"== "TRUE", diff)]
}

#remove LR columns and write to file
diff = diff[,-grep(".LR",names(diff))]
write.table(diff, file = paste0("diff.glm.ATACseq.txt"), sep = "\t", quote = F, row.names = F)

#filter only the significant genes
sig.diff = diff[which(diff$sigAny), ]
write.table(sig.diff, file = paste0("diff.significant.glm.ATACseq.txt"), sep = "\t", quote = F, row.names = F)	


#######################################################################
## Sig Stats function

cols = colnames(diff)[grepl(".sig", colnames(diff))]
stats = data.frame(Comp = gsub(".sig", "", cols))

for(i in 1:dim(stats)[1]){
  
  print(stats$Comp[i])
  print(table(diff[, cols[i]]))
  
  stats$notSig[i] = as.numeric(table(diff[, cols[i]])[1])
  stats$Sig[i] = as.numeric(table(diff[, cols[i]])[2])
  
}

write.table(stats, file = paste0("Sig.stats.glm.txt"), sep = "\t", quote = F, row.names = F);

