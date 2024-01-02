library(GEOquery)
library(curl)
library(minfi)

data_dir="/g/data3/yo4/Cancer-Epigenetics/davban/Hoang/Thien"
gse="GSE152026"
type="processed"
type="raw"
# OR 
# type="processed"

# args=commandArgs(trailingOnly = TRUE)
# # Input args
# data_dir=args[1]
# gse=args[2]
# type=args[3]

Sys.setenv("VROOM_CONNECTION_SIZE"=131072*10000)
options(timeout=600000)

### download data
if (!dir.exists(data_dir)) {dir.create(data_dir)}
raw=paste0(data_dir, "/raw")
dir.create(raw)
gseInfo=getGEO(gse,GSEMatrix=TRUE, destdir=raw)
print("Be careful, columns name might different with different GSE. Check column detailedly before setting up 'type'")
colnames(pData(phenoData(gseInfo[[1]])))

tmp=pData(phenoData(gseInfo[[1]]))[,c(grep(":ch1", colnames(pData(phenoData(gseInfo[[1]])))), 1, 8, grep("supplementary_file", colnames(pData(phenoData(gseInfo[[1]])))))]
tmp
# remove white space 
tmp=as.data.frame(apply(tmp,2,function(x)gsub('\\s+', '_',x)))
tmp
# Array colums
for (i in 1:length(tmp$supplementary_file)){
	tmp$Array[i]=strsplit(strsplit(as.character(tmp$supplementary_file[i]), split="/")[[1]][9], split="_")[[1]][3]
}
# Slide columns
for (i in 1:length(tmp$supplementary_file)){
	tmp$Slide[i]=gsub("_R.*", "", strsplit(as.character(tmp$supplementary_file[i]), split="/")[[1]][9])
}
tmp$sample_acc=rownames(tmp)

head(tmp)

# write out
write.table(tmp, paste0(raw, "/SampleSheet.tsv"), quote=F, sep="\t", dec=",", row.names=F)
write.table(pData(phenoData(gseInfo[[1]])), paste0(raw, "/SampleSheet_full.tsv"), quote=F, sep="\t", dec=",", row.names=F)
write.table(pData(phenoData(gseInfo[[1]])), '/home/minhhoang/123/SampleSheet_full.tsv', quote=F, sep="\t", dec=",", row.names=F)
# download
if(type=="raw") {
url=tmp[,grepl("supplementary_file",names(tmp))]
for (i in 1:length(url[,1])) {
	down1=as.character(url[i,1])
	name1=strsplit(as.character(url[i,1]), "/")[[1]][9]
	name2=strsplit(as.character(url[i,2]), "/")[[1]][9]
	down2=as.character(url[i,2])
	download.file(down1, paste0(raw, "/", name1), method="wget")
	download.file(down2, paste0(raw, "/", name2), method="wget")
	}}else {
		print("You are using processed dataset - without idat files")
		mSetSqFlt=getGenomicRatioSetFromGEO(gse, path=raw)
		save(mSetSqFlt, file=paste0(raw, "/mSetSqFlt_processed.RData"))
}
  
