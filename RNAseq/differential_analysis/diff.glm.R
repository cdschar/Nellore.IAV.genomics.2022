library("edgeR");

experimentName = "TIV.secondary.all"

homeDir = "~/analysis/diff/"
setwd(homeDir);

#set significance threshold
sigFDR = 0.05
sigFC = 1

#read in manifest file
fqDir = "~/pipeline/";
fqFiles = "RNAseq.sample.manifest.txt";
files = read.table(paste(fqDir, fqFiles, sep = ""), sep = ",", header = TRUE, as.is = TRUE);
files$group = paste0(files$group, "_", files$day)

#read in raw and rpm normalized files
covDir = "~/coverage/";
cts = read.table(paste0(covDir, "TIV_secondary.geneCts.detected.RPM.3.exon.csv"), sep = ",", header = T)
rpkm = read.table(paste0(covDir, "TIV_secondary.geneRpkm.detected.RPM.3.exon.csv"), sep = ",", header = T);


#establish groups for comparisons
grps = unique(files$group)
grpPairs = combn(grps, 2)


##########################
# Establish edgeR object and filter genes

#establish edgeR DGE object
dge = DGEList(as.matrix(cts[, grepl(paste(files$sample, collapse = "|"), colnames(cts))]), group = files$group[match(colnames(cts[, grepl(paste(files$sample, collapse = "|"), colnames(cts))]), files$sample)])

#calculate normalization factor
dge = calcNormFactors(dge)

#plot MDS plot
cairo_pdf(file = paste0("MDSplot.", experimentName, ".pdf"), height = 5, width = 5);
par(mai = c(1.5, 0.75, 0.75, 0.25), mgp = c(2, 0.5, 0), family = "Arial");
plotMDS.DGEList(dge, main = "MDS Plot for Count Data", labels = gsub("Sample\\_[0-9]+\\_", "", colnames(dge$counts)), cex = 0.4)
dev.off();

#######################################################################
## GLM Function

#create design matrix and estimate dispersion
design = model.matrix(~0 + files$group + files$patient)
dge = estimateDisp(dge, design)

#call glmfit function
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
    
    #mark genes that are significant by FDR criteria or by FDR and logFC
    lrt$table$sig = (lrt$table$fdr < sigFDR & abs(lrt$table$logFC) >= sigFC)
    
    #add header to the lrt$table slot that describes the comparison performed
    colnames(lrt$table) = paste0(comp, ".", names(lrt$table))
    
    #cbind diff test onto the diff matrix
    rpkm[names(lrt$table)] = NA
    rpkm[rownames(lrt$table), names(lrt$table)] = lrt$table
}

#add column that tells if a gene is significant in any comparison
rpkm$sigAny <- apply(rpkm[,(grep(".sig"== "TRUE", rpkm))], 1, any)

#remove LR columns and write to file
rpkm = rpkm[,-grep(".LR",names(rpkm))]
write.table(rpkm, file = paste0(homeDir, "diff.glm.", experimentName, ".txt"), sep = "\t", quote = F, row.names = F)

#filter only the significant genes
sig.rpkm = rpkm[which(rpkm$sigAny == "TRUE"), ]
write.table(sig.rpkm, file = paste0("diff.significant.glm.", experimentName, ".txt"), sep = "\t", quote = F, row.names = F)

#################################################################
#write table that tells how many genes are significant for each comparison
cols = colnames(rpkm)[grepl(".sig", colnames(rpkm))]
stats = data.frame(Comp = gsub(".sig", "", cols))

for(i in 1:dim(stats)[1]){
    
    print(stats$Comp[i])
    print(table(rpkm[, cols[i]]))
    
    stats$notSig[i] = as.numeric(table(rpkm[, cols[i]])[1])
    stats$Sig[i] = as.numeric(table(rpkm[, cols[i]])[2])
    
}

write.table(stats, file = paste0("Sig.stats.glm.", experimentName, ".txt"), sep = "\t", quote = F, row.names = F);
