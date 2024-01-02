### library
library(TCGAbiolinks)
library(EDASeq) 
library(DT) 
library(SummarizedExperiment) 
library(RTCGAToolbox)
library(gaia)
library(GenomicRanges)
library(maftools)
library(EDASeq)
library(clusterProfiler)
library(ggpubr) 
library(ComplexHeatmap) 


setwd("/g/data3/yo4/Cancer-Epigenetics/davban/Hoang/TCGA")



#-=-=-=-=-=-=-=-=-=-Downloading data from TCGA data portal-=-=-=-=-=-=-=-=-=-=-=-=-
# Obs: The data in the legacy database has been aligned to hg19
# 1. DNA methylation
# a. query
query.met.gbm <- GDCquery(project = "TCGA-THCA", 
                          legacy = TRUE,
                          data.category = "DNA methylation",
                          platform = "Illumina Human Methylation 450") 
                          

# b. download
GDCdownload(query.met.gbm)

# c. prepare
met.gbm.450 <- GDCprepare(query = query.met.gbm, save = TRUE, save.filename = "gbmDNAmet450k.rda", summarizedExperiment = TRUE)


# Process data
library(data.table) 
a <- as.data.frame(fread("/g/data3/yo4/Cancer-Epigenetics/davban/Hoang/TCGA/DNAmethylationTCGA.tsv",header=T))
B=571
name_edit <- as.character(B)
name <- colnames(a[-572])
for (i in 1:length(name)) {
    name_edit[i] <-paste(strsplit(as.character(name[i]), split="-")[[1]][1:3], collapse = "-")
}
colnames(a)[1:571] <- name_edit

# clinical information 
TCGA_THCA <- GDCquery_clinic(project = "TCGA-THCA", type = "Clinical")
clect  <- TCGA_THCA$submitter_id

# lọc data có case 
data <- a[,c(clect,"probeID")]

# lọc 18 probe 
rownames(data) <-  data$probeID
data_1  <- data[c("cg20392764","cg01091565","cg02633817","cg21098323" ,"cg09119967" ,"cg14911395","cg21880903","cg14258236","cg11484872","cg18049750","cg08695223",#NIFTP
            "cg03001305","cg01796223", # FA
            "cg07447773","cg26607785","cg06454226","cg19282250",#FC
            "cg02192520"),]   #fvPTC
# lọc mẫu mong muốn 
aim <- subset(TCGA_THCA, primary_diagnosis=="Papillary carcinoma, follicular variant"|
                          primary_diagnosis=="Follicular adenocarcinoma, NOS"|
                          primary_diagnosis=="Follicular carcinoma, minimally invasive")

sample <- aim$submitter_id
data_2 <- data_1[,c(sample,"probeID")]
dim(data_2)

# save data #105 fvPTC, 1FA, 1FC 
path="/g/data3/yo4/Cancer-Epigenetics/davban/Hoang/TCGA"
write.table(data_2, paste0(path,"/data_process.tsv"), sep="\t", col.names=T, quote=F, row.names=F)
write.table(aim, paste0(path,"/clinical.tsv"),sep="\t", col.names=T, quote=F, row.names=F)



# processed  data
path="/home/minhhoang/Documents/Thyroid cancer/Data and analysis/one_again/non_removebatchefect/probe_differ_tumour_normal/data_without_normal/TCGA"
library(data.table)
data <- as.data.frame(fread("/home/minhhoang/Documents/Thyroid cancer/Data and analysis/one_again/non_removebatchefect/probe_differ_tumour_normal/data_without_normal/TCGA/data_process.tsv",header=T))
rownames(data) <- data$probeID
# không có probe cg21098323
info <- as.data.frame(fread("/home/minhhoang/Documents/Thyroid cancer/Data and analysis/one_again/non_removebatchefect/probe_differ_tumour_normal/data_without_normal/TCGA/clinical.tsv",header=T))
info <- info[,c("submitter_id","primary_diagnosis","race","gender","ajcc_pathologic_stage","age_at_index","updated_datetime","vital_status","year_of_birth","ajcc_pathologic_stage")]
info$subtype <- ifelse(info$primary_diagnosis=="Papillary carcinoma, follicular variant","fvPTC",ifelse(info$primary_diagnosis=="Follicular adenocarcinoma, NOS","FA","FC"))
write.table(info, paste0(path,"/SampleSheet.tsv"), sep="\t", col.names=T, quote=F, row.names=F )

#write.table(data, paste0(path,"/predict.tsv"), sep="\t", col.names=T, quote=F, row.names=F )
# check colnames 
colnames(data)[-108]==info$submitter_id

# do with data 

tdata <- as.data.frame(t(data[,-108]))
tdata$subtype <- info$subtype 

write.table(tdata, paste0(path,"/tdata.tsv"), sep="\t", col.names=T, quote=F, row.names=F )