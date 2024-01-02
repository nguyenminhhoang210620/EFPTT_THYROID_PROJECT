library(limma)
library(minfi)
library(RColorBrewer)
library(missMethyl)
library(matrixStats)
library(umap)
library(stringr)
library(config)
library(GEOquery)
library(rtracklayer)
library(ggfortify)

data_dir="/g/data3/yo4/Cancer-Epigenetics/davban/Hoang/2022_Yao"
gse="GSE51090"
type="raw"
type="processed"


# SampleSheet.tsv
# Note 1. check this file to see column of interest (show groups comparision)
# Note 2. Make sure no space in this column name (defined in Note 1)  

# config file
# Change array type to 450K or EPIC 
# In specific case, change parameter - If not, keep it (threshold of DME database) 

configFile="/g/data3/yo4/Cancer-Epigenetics/davban/Hoang/project_1/config/GSE146003.config"
grCol="tissue:ch1"
grCol="disease state:ch1"
grCol="histology:ch1"
grcol="title source_name_ch1"
grCol="source_name_ch1"
grCol="title"
# args = commandArgs(trailingOnly = TRUE)
# # Input args
# data_dir = args[1]
# grCol = args[2]
# configFile = args[3]
# gse = args[4]
# type = args[5]

### load config
config <- config::get(file=configFile)
# get variables from config file
Pvalue <- config$Pvalue
xReactiveProbes <- config$xReactProbes
C <- config$C
lambda <- config$lambda
fdr <- config$fdr
genome1 <- config$genome1
genome2 <- config$genome2
array <- config$array
anno_dir <- config$anno_dir

Sys.setenv("VROOM_CONNECTION_SIZE"=131072*100)
options(timeout=6000)
set.seed(123)


raw <- paste0(data_dir, "/raw")
list.files(raw, recursive=TRUE)
targets <- read.table(paste0(raw, "/SampleSheet.tsv"), header=T, sep="\t", check.names=F, stringsAsFactors=F)
# make sure no space in group's names
targets$base <- str_replace_all(targets[,grCol], " ", ".")
targets$Basename <- paste0(raw, "/", targets$Slide, "_", targets$Array)
###################
for (i in 1:length(targets$supplementary_file)){
        targets$Slide[i]=strsplit(strsplit(as.character(targets$supplementary_file[i]), split="/")[[1]][9], split="_")[[1]][1]
}



QC <- paste0(data_dir, "/QC")
if (!dir.exists(QC)) {dir.create(QC)}

if (type=="raw") {
### loading data
# read in the raw data from the IDAT files
# rgSet <- read.metharray.exp(targets=targets, force=T)
setwd(raw)
system("gunzip -f *")
rgSet <- read.metharray(basenames=targets$Basename, force=T)
rgSet

sampleNames(rgSet) <- targets$title
rgSet
save(rgSet, file=paste0(raw, "/rgSet.RData"))
### Quality control
print("QC step")
# calculate p-values
print("detection P")
detP <- detectionP(rgSet)
head(detP)
write.table(detP, paste0(QC, "/detP.tsv"), sep="\t", quote=F, row.names=F)

# Normalization
print("Normalization using Quantile")
mSetSq <- preprocessQuantile(rgSet)
allProbes <- nrow(mSetSq)

# Filter based on p-values
detP <- detP[match(featureNames(mSetSq),rownames(detP)),]
keep <- rowSums(detP < Pvalue) == ncol(mSetSq)
table(keep)
mSetSqFlt <- mSetSq[keep,]
mSetSqFlt
gtPval <- data.frame(table(keep))[1,2]
ltPval <- data.frame(table(keep))[2,2]

rmPval=allProbes-gtPval

# pvalue
pdf(file = paste0(QC, "/Mean_detP.pdf"))
pal = brewer.pal(8,"Dark2")
barplot(colMeans(detP),col=pal[factor(targets$base)],las=2,cex.names=0.8,
        main="Mean detection p-values", xaxt='n')
abline(h=Pvalue, col="red")
dev.off()

}else {
load(paste0(raw, "/mSetSqFlt_processed.RData"))
sampleNames(mSetSqFlt) <- targets$title
gtPval <- 0
allProbes <- nrow(mSetSqFlt)
rmPval=allProbes-gtPval
}

# Probes with snps at CpGs
snpInfo <- getSnpInfo(mSetSqFlt)
mSetSqFlt <- dropLociWithSnps(mSetSqFlt)
mSetSqFlt
snpProbes <- rmPval - nrow(mSetSqFlt)

# exclude cross reactive probes
xReactiveProbes <- readRDS(xReactiveProbes)
d <- c()
for(i in which((featureNames(mSetSqFlt) %in% xReactiveProbes))) {d <- c(d, i)}
CrosReactiveProbes <- featureNames(mSetSqFlt)[d]
write.table(CrosReactiveProbes, paste0(QC, "/CrosReactiveProbes.tsv"), sep="\t", quote=F, col.names=NA) 
# remove Probes
if(!(featureNames(mSetSqFlt) %in% xReactiveProbes)) {print(featureNames(mSetSqFlt))}
keep <- !(featureNames(mSetSqFlt) %in% xReactiveProbes)
table(keep)
mSetSqFlt <- mSetSqFlt[keep,]
mSetSqFlt
ReactiveProbes <- data.frame(table(keep))[1,2]

# After filter Probes
FurtherUse <- allProbes-(ReactiveProbes + snpProbes + gtPval)

# filer pie chart
df <- t(data.frame(Pvalue.Filter=gtPval, SNP.Filter=(snpProbes + gtPval), CrosReactive.Filter=ReactiveProbes, FurtherUse=FurtherUse))
color <- brewer.pal(nrow(df), "Set2")
pdf(file = paste0(QC, "/Probe_filtering.pdf"))
pie_labels <- paste0(rownames(df), " = ", round(100 * df[,1]/sum(df[,1]), 2), "%")
par(mar=c(5.1, 4.1, 4.1, 6.5))
pie(df, main="Probe filtering", labels = pie_labels, col = color)
dev.off()
# UMAP
df <- getBeta(mSetSqFlt)
df <- df[complete.cases(df),]
tdf <- t(df)
Labels <- base::merge(tdf, targets, by.x=0, by.y="title")[,grCol]
tryCatch({
ump <- umap(as.matrix(tdf))
pdf(file = paste0(QC, "/UMAP.plot.pdf"))
pal = brewer.pal(8,"Dark2")
plot(ump$layout[,1], ump$layout[,2], col = pal[factor(Labels)], pch = 19, xlab="Dimension 1", ylab="Dimension 2")
legend("topright", legend=levels(factor(Labels)), text.col=pal, bty = "n")
dev.off()
}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
# PCA
pca_res <- prcomp(tdf, scale. = TRUE)
tdf=cbind(tdf,Labels)
pdf(file = paste0(QC, "/PCA.plot.pdf"))
autoplot(pca_res, data = tdf, colour = 'Labels')
dev.off()

# Density
pdf(file = paste0(QC, "/Dens.Plot.pdf"))
pal = brewer.pal(8,"Dark2")
densityPlot(getBeta(mSetSqFlt), sampGroups = targets$source_name_ch1, main="densityPlot_normalized")
dev.off()

# write processed data
methProbes <- paste0(data_dir, "/methProbes")
if (!dir.exists(methProbes)) {dir.create(methProbes)}
save(mSetSqFlt, file=paste0(methProbes, "/mSetSqFlt.RData"))

### get Beta_value
bVals <- as.data.frame(getBeta(mSetSqFlt))
bigTable <- as.data.frame(apply(bVals, 2, function(x) round(x, 2)))
probeID <- row.names(bigTable)
bigTable$probeID <- probeID
bigTable[1:5,1:5]

# coordination
# genome1 _38
probeID <- read.table('/g/data3/yo4/Cancer-Epigenetics/davban/Hoang/project_1/config/EPIC_CpG.tsv', header=T)
dt <- base::merge(bigTable, probeID, by="probeID")
bT.g1 <- dt[,c("CpG_chrm", "CpG_beg", "CpG_end", "probe_strand", "probeID", targets$title)]
bT.g1 <- bT.g1[order(bT.g1$CpG_chrm, bT.g1$CpG_beg),]
names(bT.g1)[1:5] <- c("chr", "start", "end", "strand", "probeID")
bT.g1[1:5,1:10]
write.table(bT.g1, paste0(methProbes, "/bigTable.", genome1,  ".tsv"), sep="\t", col.names=T, quote=F, row.names=F)
# genome2_19
probeID <- read.table('/g/data3/yo4/Cancer-Epigenetics/davban/Hoang/project_1/config/EPIC.hg19.manifest.tsv', header=T)
dt <- base::merge(bigTable, probeID, by="probeID")
bT.g2 <- dt[,c("CpG_chrm", "CpG_beg", "CpG_end", "probe_strand", "probeID", targets$title)]
bT.g2 <- bT.g2[order(bT.g2$CpG_chrm, bT.g2$CpG_beg),]
names(bT.g2)[1:5] <- c("chr", "start", "end", "strand", "probeID")
bT.g2[1:5,1:10]
write.table(bT.g2, paste0(methProbes, "/bigTable.", genome2,  ".tsv"), sep="\t", col.names=T, quote=F, row.names=F)
