# calling packages
library(minfi)
library(RColorBrewer)
library(missMethyl)
library(matrixStats)
library(Gviz)
library(DMRcate)
library(stringr)
library(DMRcatedata)
library(data.table)

# creating path 
data_dir  <- "/home/minhhoang/Documents/Thyroid cancer/TOTAL/Calling_DMP/First_calling_DMP"
# Loading file config 
configFile <- paste0(data_dir,"/Data/config/GSE146003.config")


### load config
config <- config::get(file=configFile)
Pvalue <- config$Pvalue
xReactiveProbes <- config$xReactProbes
C <- config$C
lambda <- config$lambda
fdr <- config$fdr
genome1 <- config$genome1
genome2 <- config$genome2
type <- config$type
array <- config$array
anno_dir <- config$anno_dir


# load processed expressionSet
data <- as.data.frame(fread(paste0(data_dir,"/Data/traindata.tsv")),header=T)
data <- data[!(data$subtype=="NG"),]
tdata <- as.data.frame(t(data[,-c(21254,21255)]))
colnames(tdata) <- data$title
tdata[1:5,1:5]
dim(tdata)
bVals <- tdata

### prepare design.cs
targets <- as.data.frame(fread(paste0(data_dir,'/Data/info_train.tsv'), header=T, sep="\t", check.names=F, stringsAsFactors=F,fill=T))
targets <- targets[!(targets$subtype=="NG"),]
info <- targets[,c("title","subtype")]
designTSV <- info
names(designTSV) <- c("title","groups")
write.table(designTSV, paste0(data_dir, "/Result/design.tsv"), sep="\t", quote=F, row.names=F)

### DMPs - DMRs calling


# possible comparisons
pair <- combn(unique(designTSV$groups),2)

for (i in 1:ncol(pair)) {
	print (i)
	gr <- pair[,i]
	print(gr)
	g1 <- gr[1]
	g2 <- gr[2]
	SubGr <- designTSV[designTSV$groups%in%gr,] 
	sN <- SubGr$title
	dt <- bVals[,sN]
	type <- ifelse(grepl(g1, SubGr$groups), "case", "control")
	type <- relevel(factor(type), "control")
	design <- model.matrix(~type)
	colnames(design) <- gsub("type", "", colnames(design))
	dt <- dt[complete.cases(dt),]
	myAnnotation <- cpg.annotate("array", as.matrix(dt), arraytype=array, analysis.type="differential", design=design, coef="case", fdr=fdr, what="Beta")
	cpSites <- as.data.frame(myAnnotation@ranges)
	DMP.query <- cpSites[cpSites$is.sig==TRUE,]
	if (nrow(DMP.query) ==0 ) {
  		print (paste0("Your contrast returned no individually significant probes under threshold fdr=", fdr, " Stop process pair ", g1, "_", g2))
		}else {
		DMR <- paste0(data_dir, "/DMR")
		if (!dir.exists(DMR)) {dir.create(DMR)}
		path <- paste0(data_dir, "/DMR/", paste0(g1, "-", g2))
		if (!dir.exists(path)) {dir.create(path)}
		}
	dmrcate.res <- dmrcate(myAnnotation, lambda=lambda, C=C, min.cpgs=2, betacutoff=0.2)
	DMRquery <- length(dmrcate.res@coord)
	if (DMRquery > 0) {
		saveRDS(dmrcate.res, file = paste0(path, "/DMR.rds"))
		dmrcate.res <- readRDS(paste0(path,"/DMR.rds"))
		# genome2
		dmrcate.ranges <- data.frame(extractRanges(dmrcate.res, genome=genome2))
		write.table(dmrcate.ranges, paste0(path, "/DMR.", genome2, ".tsv"))
		groups <- c(g1="magenta", g2="forestgreen")
		names(groups) <- c(g1, g2)
		cols <- groups[as.character(SubGr$groups)]
		n <- round(nrow(data.frame(dmrcate.ranges)))
		if (n>20) {
			print ("DMR gt 20 - plot first 20 DMRs ")
			pdf(paste0(path, "/DMR.", genome2, ".pdf"), width=20, height=20)
			for (idx in 1:20){
			print(paste0("DMR", idx))
			tryCatch({
			DMR.plot(ranges=dmrcate.ranges, dmr=idx, CpGs=as.matrix(dt), what="Beta", arraytype = array, phen.col=cols, genome=genome2)
			}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
			}
			dev.off()
		}else {
			print ("DMR lt 20 - plot all")
			pdf(paste0(path, "/DMR.", genome2, ".pdf"), width=20, height=20)
			for (idx in 1:n){
			print(paste0("DMR", idx))
			tryCatch({
			DMR.plot(ranges=dmrcate.ranges, dmr=idx, CpGs=as.matrix(dt), what="Beta", arraytype = array, phen.col=cols, genome=genome2)
			}, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
			}
			dev.off()
			}
	}else if (DMRquery==0){
			print("Your contrast returned no individually significant DMR. Stop process DMR calling, moving to extract DMPs process")
			}
	print("DMPs calling")
	DMPs <- cpSites[cpSites$is.sig==TRUE,]
	DMPs$probeID <- rownames(DMPs)
	DMP <- paste0(data_dir, "/DMP")
	if (!dir.exists(DMP)) {dir.create(DMP)}
	path <- paste0(data_dir, "/DMP/", paste0(g1, "-", g2))
	if (!dir.exists(path)) {dir.create(path)}
	Beta.g1 <- apply(dt[,SubGr[SubGr$groups==g1, "title"]], 1, mean)
	Beta.g2 <- apply(dt[,SubGr[SubGr$groups==g2, "title"]], 1, mean)
	delta <- data.frame(delta=Beta.g1 - Beta.g2)
	DMPs <- base::merge(DMPs, delta, by.x="probeID", by.y=0)
	pdf(file = paste0(path, "/delta_dist.pdf"))
	hist(abs(DMPs$delta), xlab="Beta value", main="Delta distribution")
	abline(v=0.2,col="red")
	dev.off()
	pdf(file = paste0(path, "/top.DMPs.pdf"))
	par(mfrow=c(2,2))
	sapply(head(DMPs, n=100)$probeID, function(cpg){
	plotCpg(dt, cpg=cpg, ylab="bVals", pheno = SubGr$groups)})
	dev.off()
	pdf(file = paste0(path, "/last.DMPs.pdf"))
	par(mfrow=c(2,2))
	sapply(tail(DMPs, n=100)$probeID, function(cpg){
	plotCpg(dt, cpg=cpg, ylab="bVals", pheno = SubGr$groups)})
	dev.off()
	SigDMPs <- DMPs
	# DMC
	# genome1
	probes <- read.table(paste0(data_dir,"/Data/config/HM27.hg38.manifest.tsv"), header=T, fill=TRUE)
	DMCprobes <- base::merge(SigDMPs, probes, by.x = "probeID", by.y = "probeID")
	DMCloc <- DMCprobes[,c("CpG_chrm", "CpG_beg", "CpG_end", "probe_strand", "probeID", "gene", "ind.fdr", "delta")]
	DMClocSorted <- DMCloc[order(DMCloc$CpG_chrm, DMCloc$CpG_beg),]
	names(DMClocSorted) <- c("chr", "start", "end", "strand", "probeID", "gene", "FDR", "mean.Delta")
	write.table(DMClocSorted, paste0(path, "/DMPs_filtered.", genome1, ".tsv"), quote=F, sep="\t", row.names=F)
	# genome2
	probes <- read.table(paste0(data_dir,"/Data/config/HM27.hg19.manifest.tsv"), header=T, fill=TRUE)
	DMCprobes <- base::merge(SigDMPs, probes, by.x = "probeID", by.y = "probeID")
	DMCloc <- DMCprobes[,c("CpG_chrm", "CpG_beg", "CpG_end", "probe_strand", "probeID", "gene", "ind.fdr", "delta")]
	DMClocSorted <- DMCloc[order(DMCloc$CpG_chrm, DMCloc$CpG_beg),]
	names(DMClocSorted) <- c("chr", "start", "end", "strand", "probeID", "gene", "FDR", "mean.Delta")
	write.table(DMClocSorted, paste0(path, "/DMPs_filtered.", genome2, ".tsv"), quote=F, sep="\t", row.names=F)
}

