# load data GSE121377
bigtable <- read.table('/g/data3/yo4/Cancer-Epigenetics/davban/Hoang/project_1/GSE121377/methProbes/bigTable.hg19.tsv', header=T)
targets <- read.table('/g/data3/yo4/Cancer-Epigenetics/davban/Hoang/project_1/GSE121377/raw/SampleSheet.tsv', header=T, sep="\t", check.names=F, stringsAsFactors=F) 
data_1 <- bigtable[,c("probeID","X10T._genomic_DNA_from_thyroid_tumor_tissue_3","X23T._genomic_DNA_from_thyroid_tumor_tissue_13","X27T._genomic_DNA_from_thyroid_tumor_tissue_16",
	                  "X34T._genomic_DNA_from_thyroid_tumor_tissue_20","X39T._genomic_DNA_from_thyroid_tumor_tissue_31","X40T._genomic_DNA_from_thyroid_tumor_tissue_33")]
names(targets)[1] <- "type"
a <- subset(targets[,3], targets$type=="NIFTP")
names(data_1)[2:7] <- a

#load data GSE97466 
bigtable <- read.table('/g/data3/yo4/Cancer-Epigenetics/davban/Hoang/project_1/GSE97466/methProbes/bigTable.hg19.tsv', header=T)
targets <- read.table('/g/data3/yo4/Cancer-Epigenetics/davban/Hoang/project_1/GSE97466/raw/SampleSheet.tsv', header=T, sep="\t", check.names=F, stringsAsFactors=F)
# 5 sample: follicular_adenoma, 3 samples: follicular_adenoma/Hürthle_cell
data_FA <- bigtable[, c("probeID","X201T_benign","X286T_benign","X294T_benign","X285T_benign", "X267T_benign", "X284T_benign","X287T_benign","X275T_benign")]
# 4 samples: follicullar_thyroid_cancer: 4 sample  minimally_invasive_follicular_carcinomas
data_FC <- bigtable[,c("X543T_tumor","X549T_tumor","X541T_tumor","X555T.A_tumor","X587T_tumor","X548T_tumor","X540T_tumor","X550T_tumor")]
# 6 samples NG 
data_NG <- bigtable[,c("X545N_benign","X314T_benign","X43T_benign","X48T_benign","X88T_benign","X93T_benign")]
# 7 sample fvPTC 
data_fvPTC <- bigtable[,c("X346T_tumor", "X347T_tumor", "X343T_tumor","X395T_tumor", "X473T_tumor","X333T_tumor","X328T_tumor")]
data_2 <- cbind(data_FA,data_FC,data_NG,data_fvPTC)
colnames(data_2)[2:30] <- gsub("X", "",colnames(data_2)[2:30])
colnames(data_2)[13] <- "555T-A_tumor"
#Không cần hàng này 
#name_col   <- c(rep("FA",8), rep("FC",8), rep("NG",6), rep("fvPTC",7), rep("non-neoplastic_adjacent_tissue",50))



#load data GSE53051
bigtable <- read.table('/g/data3/yo4/Cancer-Epigenetics/davban/Hoang/project_1/GSE53051/methProbes/bigTable.hg19.tsv', header=T)
targets <- read.table('/g/data3/yo4/Cancer-Epigenetics/davban/Hoang/project_1/GSE53051/raw/SampleSheet.tsv', header=T, sep="\t", check.names=F, stringsAsFactors=F)
data_3 <- bigtable[,c(1:5, grep("Thyroid",colnames(bigtable)))]
data_3<- data_3[,c(5, grep("FA",colnames(data_3)), grep("FC",colnames(data_3)), grep("FVPTC",colnames(data_3)))]

#không cần hàng này 
#names(data_3)[2:32] <- c(rep("FA",11),rep("FC",10),rep("fvPTC",10))

#load data GSE51090 
bigtable <- read.table('/g/data3/yo4/Cancer-Epigenetics/davban/Hoang/project_1/GSE51090/methProbes/bigTable.hg19.tsv', header=T)
targets <- read.table('/g/data3/yo4/Cancer-Epigenetics/davban/Hoang/project_1/GSE51090/raw/SampleSheet.tsv', header=T, sep="\t", check.names=F, stringsAsFactors=F)
data_4 <- bigtable[,-(1:4)]
# 5 fvPTC, 18 FA, 18 FC 
data_4 <- data_4[,c(1,49:66,67:84,44:48)]
# không cần hàng này 
#names(data_4)[2:42] <- c(rep("fvPTC",5),rep("FA",18),rep("FC",18))


# seach normal and merge together
# GSE121377_ 7 sample 
bigtable <- read.table('/g/data3/yo4/Cancer-Epigenetics/davban/Hoang/project_1/GSE121377/methProbes/bigTable.hg19.tsv', header=T)
targets <- read.table('/g/data3/yo4/Cancer-Epigenetics/davban/Hoang/project_1/GSE121377/raw/SampleSheet.tsv', header=T, sep="\t", check.names=F, stringsAsFactors=F)
names(targets)[1] <- "type"
targets_normal <- subset(targets, type=="Normal")
data_normal_1 <- bigtable[,c("probeID","X09N._genomic_DNA_from_thyroid_normal_tissue_1","X23N._genomic_DNA_from_thyroid_normal_tissue_12","X27N._genomic_DNA_from_thyroid_normal_tissue_15",
                          "X34N._genomic_DNA_from_thyroid_normal_tissue_19","X02N._genomic_DNA_from_thyroid_normal_tissue_24","X39N._genomic_DNA_from_thyroid_normal_tissue_30",
                          "X40N._genomic_DNA_from_thyroid_normal_tissue_32")]

a <- subset(targets[,3], targets$type=="Normal")
names(data_normal_1)[2:8] <- a
head(data_normal_1)

#GSE97466 - 50 sample 
bigtable <- read.table('/g/data3/yo4/Cancer-Epigenetics/davban/Hoang/project_1/GSE97466/methProbes/bigTable.hg19.tsv', header=T)
targets <- read.table('/g/data3/yo4/Cancer-Epigenetics/davban/Hoang/project_1/GSE97466/raw/SampleSheet.tsv', header=T, sep="\t", check.names=F, stringsAsFactors=F)
data_normal_2 <- bigtable[,c(5, grep("normal",colnames(bigtable)))]
colnames(data_normal_2)[2:51] <- gsub("X", "",colnames(data_normal_2)[2:51])
colnames(data_normal_2)[22] <- "586N-A_normal"
colnames(data_normal_2)[25] <- "583N-A_normal"

#GSE53051 - 12 samples
bigtable <- read.table('/g/data3/yo4/Cancer-Epigenetics/davban/Hoang/project_1/GSE53051/methProbes/bigTable.hg19.tsv', header=T)
targets <- read.table('/g/data3/yo4/Cancer-Epigenetics/davban/Hoang/project_1/GSE53051/raw/SampleSheet.tsv', header=T, sep="\t", check.names=F, stringsAsFactors=F)
data_3 <- bigtable[,c(1:5, grep("Thyroid",colnames(bigtable)))]
data_normal_3 <- data_3[,c(5, grep("Normal",colnames(data_3)))]

#GSE51090 - 8 sample 
bigtable <- read.table('/g/data3/yo4/Cancer-Epigenetics/davban/Hoang/project_1/GSE51090/methProbes/bigTable.hg19.tsv', header=T)
targets <- read.table('/g/data3/yo4/Cancer-Epigenetics/davban/Hoang/project_1/GSE51090/raw/SampleSheet.tsv', header=T, sep="\t", check.names=F, stringsAsFactors=F)
data_normal_4 <- bigtable[,c(5,grep("normal",colnames(bigtable)))]


# merge normal and all case 

path = '/g/data3/yo4/Cancer-Epigenetics/davban/Hoang/project_1/newdata'
df_list = list(data_1, data_normal_1, data_2, data_normal_2, data_3, data_normal_3, data_4, data_normal_4)
dt_19_27 <- Reduce(function(x, y) merge(x, y, all=TRUE, by='probeID'), df_list)
row.has.na <- apply(dt_19_27, 1, function(x){any(is.na(x))})
sum(row.has.na)
final.filtered <- dt_19_27[!row.has.na,]
probeID <- read.table('/g/data3/yo4/Cancer-Epigenetics/davban/Hoang/project_1/config/HM27.hg19.manifest.tsv', header=T)
dt <- merge(final.filtered, probeID, by="probeID")
colname_total  <- names(final.filtered)[2:185]
bT.g2 <- dt[,c("CpG_chrm", "CpG_beg", "CpG_end", "probe_strand", "probeID", colname_total)]
bT.g2 <- bT.g2[order(bT.g2$CpG_chrm, bT.g2$CpG_beg),]
names(bT.g2)[1:5] <- c("chr", "start", "end", "strand", "probeID")
dim(bT.g2)
bT.g2[1:5,1:10]
write.table(bT.g2, paste0(path, "/bigtable_merge_19_normal_full.tsv"), sep="\t", col.names=T, quote=F, row.names=F)

# check:
library(data.table)
data <- as.data.frame(fread("/g/data3/yo4/Cancer-Epigenetics/davban/Hoang/project_1/newdata/bigtable_merge_19_normal_full.tsv",header=T))
info <- read.table("/g/data3/yo4/Cancer-Epigenetics/davban/Hoang/project_1/newdata/info_merge_27k_full_tsv", header=T, sep="\t", check.names=F, stringsAsFactors=F)





#explore data 
library(ggplot2)
library(ggfortify)
path = '/g/data3/yo4/Cancer-Epigenetics/davban/Hoang/project_1/exploreDATA'
data <- read.table(paste0(path,'/bigtable_merge_19.tsv'), header=T)
info <- read.table(paste0(path,"info.tsv"), header=T)
data <- data[,-(1:4)]
colname_total  <- c(rep('NIFTP',6), rep("FA",8), rep("FC",8), rep("NG",6), rep("fvPTC",7), rep("FA",11),rep("FC",10),rep("fvPTC",10), rep("FA",18),rep("FC",18),rep("fvPTC",5))
data <- data[complete.cases(data),]
final_df <- as.data.frame(t(data[,-1]))
colnames(final_df) <- data$probeID
dim(final_df)
final_df[1:5,1:10]
res.pca <- prcomp(final_df[,c(1:21254)], scale=TRUE)
summary(res.pca)
final_df <- cbind(final_df,colname_total)
pdf(file = paste0(path, "/PCA.plot_1.pdf"))
autoplot(res.pca, data = final_df , colour = 'colname_total')
dev.off()

data <- read.table(paste0(path,'/bigtable_merge_19.tsv'), header=T)
info <- read.table(paste0(path,'/info.tsv'), header=T)
data <- data[,-(1:4)]
colname_total <- info$Platform
data <- data[complete.cases(data),]
final_df <- as.data.frame(t(data[,-1]))
colnames(final_df) <- data$probeID
res.pca <- prcomp(final_df[,c(1:21254)], scale=TRUE)
summary(res.pca)
final_df <- cbind(final_df,colname_total)
pdf(file = paste0(path, "/PCA.plot_1.pdf"))
autoplot(res.pca, data = final_df , colour = 'colname_total')
dev.off()




#vẽ trên máy tính 
data <- read.table('/home/minhhoang/Documents/Thyroid cancer/Data and analysis/MergeData_19_27k/bigtable_merge_19.tsv',header=T)
colname_total  <- c(rep('NIFTP',6), rep("FA",8), rep("FC",8), rep("NG",6), rep("fvPTC",7), rep("FA",11),rep("fvPTC",10),rep("FC",10), rep("fvPTC",5),rep("FA",18),rep("FC",18))
names(data)[2:108] <- colname_total 
library(reshape2)
tdata <- melt(data, id.vars="probeID")
ggplot(tdata,aes(x=variable,y=value,fill=variable)) + geom_boxplot()


# file inforrmation 
# GSE121377
targets_1 <- read.table('/g/data3/yo4/Cancer-Epigenetics/davban/Hoang/project_1/GSE121377/raw/SampleSheet.tsv', header=T, sep="\t", check.names=F, stringsAsFactors=F)
targets_2 <- read.table('/g/data3/yo4/Cancer-Epigenetics/davban/Hoang/project_1/GSE121377/raw/SampleSheet_full.tsv', header=T, sep="\t", check.names=F, stringsAsFactors=F)
names(targets_1)[1] <- "type"
names(targets_2)[33] <- "types"
info_1 <- subset(targets_1, type=="NIFTP")
info_1 <- info_1[,c(3,1)]
info_1$Age <- rep("NA",6)
info_1$Sex <- rep("NA",6)
info_1$Race <- rep("Korean",6)
info_1$Platform <- rep("850k",6)
info_1$Study <- rep("Park",6)

# GSE97466
targets_1 <- read.table('/g/data3/yo4/Cancer-Epigenetics/davban/Hoang/project_1/GSE97466/raw/SampleSheet.tsv', header=T, sep="\t", check.names=F, stringsAsFactors=F)
targets_2 <- read.table('/g/data3/yo4/Cancer-Epigenetics/davban/Hoang/project_1/GSE97466/raw/SampleSheet_full.tsv', header=T, sep="\t", check.names=F, stringsAsFactors=F)
names(targets_1)[3] <- "type"
names(targets_2)[39] <- "type"
# chưa có fvPTC
info_2 <- subset(targets_1, type=="follicular_adenoma"|type=="follicular_adenoma/Hürthle_cell"|type=="follicullar_thyroid_cancer"|type=="minimally_invasive_follicular_carcinomas"|type=="nodular_goiter")
info_subtype <- info_2[,c(6,3,1,2,4)]
info_subtype <- info_subtype[order(info_subtype$type),]
# thêm fvPTC 
info_fvPTC <- subset(targets_1, title=="328T_tumor"|title=="473T_tumor"|title=="347T_tumor"|title=="343T_tumor"|title=="395T_tumor"|title=="333T_tumor"|title=="346T_tumor")
info_fvPTC <- info_fvPTC[,c(6,3,1,2,4)]
# Kết hợp hai bảng 
info_2 <- rbind(info_subtype, info_fvPTC)
info_2$Platform <- rep("450k",29)
info_2$Study <- rep("Bisarro",29)
names(info_2)[3:5] <- c("Age","Sex","Race")

# GSE53051
targets_1 <- read.table('/g/data3/yo4/Cancer-Epigenetics/davban/Hoang/project_1/GSE53051/raw/SampleSheet.tsv', header=T, sep="\t", check.names=F, stringsAsFactors=F)
targets_2 <- read.table('/g/data3/yo4/Cancer-Epigenetics/davban/Hoang/project_1/GSE53051/raw/SampleSheet_full.tsv', header=T, sep="\t", check.names=F, stringsAsFactors=F)
names(targets_2)[47] <- "typecancer"
info_thyroid <- subset(targets_2, typecancer=="thyroid")
info_3 <- info_thyroid[,c("title","notes:ch1","age:ch1","Sex:ch1")]
names(info_3)[1:4] <- c("title","type","Age","Sex")
info_3 <- subset(info_3, type=="follicular.adenoma"|type=="follicular.variant.papillary.carcinoma"|type=="follicular.carcinoma")
# Fa, FC, fvPTC 
info_3 <- info_3[order(info_3$type),]
info_3$Race <- rep("NA",31)
info_3$Platform <- rep("450k",31)
info_3$Study <- rep("Timp", 31)

# GSE51090
targets_1 <- read.table('/g/data3/yo4/Cancer-Epigenetics/davban/Hoang/project_1/GSE51090/raw/SampleSheet.tsv', header=T, sep="\t", check.names=F, stringsAsFactors=F)
targets_2 <- read.table('/g/data3/yo4/Cancer-Epigenetics/davban/Hoang/project_1/GSE51090/raw/SampleSheet_full.tsv', header=T, sep="\t", check.names=F, stringsAsFactors=F,fill=T)
info_4 <- targets_2[,c("title","source_name_ch1")]
names(info_4) <- c("title","type")
info_4 <- subset(info_4, type=="FA"|type=="FAh"|type=="FTC"|type=="FTC+PDTC"|type=="fvPTC")
t <- gsub(" ", "_", info_4$title)
info_4$title <- t
# 62 = Fah, 82=FTC+PDTC
a <- c(rep("fvPTC",5),rep("FA",18),rep("FC",18))
info_4$type <- a
#FA, FC , fvPTC
info_4 <- info_4[order(info_4$type),]
info_4$Age <- rep("NA",41)
info_4$Sex <- rep("NA",41)
info_4$Race <- rep("NA",41)
info_4$Platform <- rep("270K",41)
info_4$Study <- rep("Mancikova",41)



#### infor about normal samples
#GSE121377
targets_1 <- read.table('/g/data3/yo4/Cancer-Epigenetics/davban/Hoang/project_1/GSE121377/raw/SampleSheet.tsv', header=T, sep="\t", check.names=F, stringsAsFactors=F)
targets_2 <- read.table('/g/data3/yo4/Cancer-Epigenetics/davban/Hoang/project_1/GSE121377/raw/SampleSheet_full.tsv', header=T, sep="\t", check.names=F, stringsAsFactors=F)
names(targets_1)[1] <- "type"
info_normal_1 <- subset(targets_1, type=="Normal")
info_normal_1 <- info_normal_1[,c(3,1)]
info_normal_1$Age <- rep("NA",7)
info_normal_1$Sex <- rep("NA",7)
info_normal_1$Race <- rep("Korean",7)
info_normal_1$Platform <- rep("850k",7)
info_normal_1$Study <- rep("Park",7)
# GSE97466
targets_1 <- read.table('/g/data3/yo4/Cancer-Epigenetics/davban/Hoang/project_1/GSE97466/raw/SampleSheet.tsv', header=T, sep="\t", check.names=F, stringsAsFactors=F)
targets_2 <- read.table('/g/data3/yo4/Cancer-Epigenetics/davban/Hoang/project_1/GSE97466/raw/SampleSheet_full.tsv', header=T, sep="\t", check.names=F, stringsAsFactors=F)
names(targets_1)[3] <- "type"
names(targets_2)[39] <- "type"
info_normal_2 <- subset(targets_1, type=="non-neoplastic_adjacent_tissue")
info_normal_2 <- info_normal_2[,c(6,3,1,2,4)]
info_normal_2$Platform <- rep("450k",50)
info_normal_2$Study <- rep("Bisarro",50)
names(info_normal_2)[3:5] <- c("Age","Sex","Race")
# GSE53051
targets_1 <- read.table('/g/data3/yo4/Cancer-Epigenetics/davban/Hoang/project_1/GSE53051/raw/SampleSheet.tsv', header=T, sep="\t", check.names=F, stringsAsFactors=F)
targets_2 <- read.table('/g/data3/yo4/Cancer-Epigenetics/davban/Hoang/project_1/GSE53051/raw/SampleSheet_full.tsv', header=T, sep="\t", check.names=F, stringsAsFactors=F)
names(targets_2)[47] <- "typecancer"
info_thyroid <- subset(targets_2, typecancer=="thyroid")
info_normal_3 <- info_thyroid[,c("title","notes:ch1","age:ch1","Sex:ch1")]
names(info_normal_3)[1:4] <- c("title","type","Age","Sex")
info_normal_3 <- subset(info_normal_3, type=="na")
info_normal_3$Race <- rep("NA",12)
info_normal_3$Platform <- rep("450k",12)
info_normal_3$Study <- rep("Timp", 12)
# GSE51090
targets_1 <- read.table('/g/data3/yo4/Cancer-Epigenetics/davban/Hoang/project_1/GSE51090/raw/SampleSheet.tsv', header=T, sep="\t", check.names=F, stringsAsFactors=F)
targets_2 <- read.table('/g/data3/yo4/Cancer-Epigenetics/davban/Hoang/project_1/GSE51090/raw/SampleSheet_full.tsv', header=T, sep="\t", check.names=F, stringsAsFactors=F,fill=T)
info_normal_4 <- targets_2[,c("title","source_name_ch1")]
names(info_normal_4) <- c("title","type")
info_normal_4 <- subset(info_normal_4, type=="NT")
t <- gsub(" ", "_", info_normal_4$title)
info_normal_4$title <- t
info_normal_4$Age <- rep("NA",8)
info_normal_4$Sex <- rep("NA",8)
info_normal_4$Race <- rep("NA",8)
info_normal_4$Platform <- rep("270K",8)
info_normal_4$Study <- rep("Mancikova",8)

#merge
infototal <- rbind(info_1, info_normal_1,info_2, info_normal_2,info_3, info_normal_3, info_4, info_normal_4)
a <- c(rep("NIFTP",6),rep("normal",7),
      rep("FA",8),rep("FC",8),rep("NG",6),rep("fvPTC",7),rep("normal",50),
      rep("FA",11),rep("FC",10),rep("fvPTC",10),rep("normal",12),
      rep("FA",18),rep("FC",18),rep("fvPTC",5), rep("normal",8))
infototal$subtype <- a
path = '/g/data3/yo4/Cancer-Epigenetics/davban/Hoang/project_1/newdata'
write.table(infototal, paste0(path, "/info_merge_27k_full_tsv"), sep="\t", col.names=T, quote=F, row.names=F)



# remove batch effect of data have 27k 
library(factoextra)
data <- read.table('/home/minhhoang/Documents/Thyroid cancer/Data and analysis/MergeData_19_27k/bigtable_merge_19.tsv',header=T)
data <- data[,-c(1:4)]
dim(data)
final_df <- as.data.frame(t(data[,-1]))
colnames(final_df) <- data$probeID
res.pca <- prcomp(final_df[,c(1:21254)], center=TRUE, scale=TRUE, retx=TRUE)
summary(res.pca)
fviz_eig(res.pca)
fviz_pca_ind(res.pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
             )
info <- read.table('/home/minhhoang/Documents/Thyroid cancer/Data and analysis/MergeData_19_27k/info.tsv',header=T)
final_df <- cbind(final_df,info$Study,info$Platform)
names(final_df)[21255:21256] <- c("Study","Platform")
group_platfrom <- as.factor(final_df$Platform)
fviz_pca_ind(res.pca,
             col.ind = group_platfrom, # color by groups
             palette = c("#00AFBB", "#E7B800", "#FC4E07"),
             addEllipses = TRUE, # Concentration ellipses
             ellipse.type = "confidence",
             legend.title = "Groups",
             repel = TRUE
             )
group_study <- as.factor(final_df$Study)
fviz_pca_ind(res.pca,
             col.ind = group_study, # color by groups
             addEllipses = TRUE, # Concentration ellipses
             ellipse.type = "confidence",
             legend.title = "Groups",
             repel = TRUE
             )
lmp <- function (modelobject) {
    if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
    f <- summary(modelobject)$fstatistic
    p <- pf(f[1],f[2],f[3],lower.tail=F)
    attributes(p) <- NULL
    return(p)
}

r = matrix(data=NA, nrow=10, ncol=2)
rownames(r) <- outcomes 
colnames(r) <- vars
p  = matrix(data=NA, nrow=10, ncol=2)
rownames(p) <- outcomes
colnames(p) <- vars
for (var in vars){
  for (outcome in outcomes){
  print(paste0(var, "..vs..", outcome))
  f <- as.formula(
      paste(outcome, paste(var, collapse = " + "), sep = " ~ ")
      )
  # do regression
    model.lm <- lm(f, data=dt)
    summary(model.lm)
    rsq <- summary(model.lm)$r.squared
	pvalue <- lmp(model.lm)
	print(paste0(round(sqrt(rsq), 4), "..vs..", pvalue))
	r[outcome,var] = round(sqrt(rsq), 4)
	p[outcome,var] = -log(pvalue,10)
	}
}
path = "/home/minhhoang/Documents/Thyroid cancer/Data and analysis/MergeData_19_27k"
write.table(r, paste0(path, "/r.tsv"), sep="\t", col.names=T, quote=F, row.names=F)
write.table(p, paste0(path, "/p.tsv"), sep="\t", col.names=T, quote=F, row.names=F)

library(ggplot2)
library(reshape2)
# correlation
r <- read.table("/home/minhhoang/Documents/Thyroid cancer/Data and analysis/MergeData_19_27k/r.tsv", header=T, sep="\t")
r$id <- rownames(r)
r <- melt(r)
colnames(r) <- c("id", "Feature", "R")
# p value
p <- read.table("/home/minhhoang/Documents/Thyroid cancer/Data and analysis/MergeData_19_27k/p.tsv", header=T, sep="\t")
p$id <- rownames(p)
p <- melt(p)
colnames(p) <- c("id", "Feature", "Pvalue")
# merge
dt <- merge(r, p, by=c("id", "Feature"))

# plot 1
ggplot(data=dt, aes(x=R, y=Pvalue, color=Feature, shape=Feature)) +
 geom_point(size=5) + theme(text = element_text(size=30))  

# plot 2
ggplot(data=dt, aes(x=R, y=Pvalue, color=Feature, shape=Feature)) +
 geom_point(size=5) + ylim(0, 30) + theme(text = element_text(size=30))


 # Remove batch efect PC1 PC2 PC3 
library(factoextra)
library(data.table)
data <- as.data.frame(fread('/home/minhhoang/Documents/Thyroid cancer/Data and analysis/MergeData_19_27k/bigtable_merge_19.tsv'))
data <- data[,-c(1:5)]
dim(data)
final_df <- as.data.frame(t(data[,-1]))
colnames(final_df) <- data$probeID
# Calculate mean
mu <- colMeans(final_df)
# Caculate PCA 
Xpca <- prcomp(final_df, center=TRUE, scale=TRUE, retx=TRUE)
# Remove batch efect
Xhat <- Xpca$x[,4:107] %*% t(Xpca$rotation[,4:107])
Xhat <- scale(Xhat, center=-mu, scale=FALSE)
# Caculate PCA 
Xpca_again <- prcomp(Xhat[,c(1:21254)], center=TRUE, scale=TRUE, retx=TRUE)
summary(Xpca_again)
fviz_eig(Xpca_again)
fviz_pca_ind(Xpca_again,
             col.ind = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE   
             )
info <- read.table('/home/minhhoang/Documents/Thyroid cancer/Data and analysis/MergeData_19_27k/info.tsv',header=T)
Xhat <- cbind(Xhat,Study=info$Study,Platform=info$Platform,Type=info$type)
group_type <- as.factor(Xhat$Type)
fviz_pca_ind(Xpca_again,
             col.ind = group_type, 
             addEllipses = TRUE,
             ellipse.type = "confidence",
             legend.title = "Groups",
             repel = TRUE
             )
a <- c(rep("NIFTP",6),rep("FA",8),rep("FC",8),rep("NG",6),rep("fvPTC",7),rep("FA",11),rep("FC",10),rep("fvPTC",10),rep("FA",18),rep("FC",18),rep("fvPTC",5))
final_df <- cbind(final_df,Type_1=a)
group_type_1 <- as.factor(final_df$Type_1)
fviz_pca_ind(Xpca_again,
             col.ind = group_type_1, 
             legend.title = "Groups",
             repel = TRUE
             )




#PCA NIFTP and fvPTC (6 NIFTP,22 fvPTC)
Xhat <- as.data.frame(Xhat)
Xhat$type <- a
Xhat[1:5,1:5]

data_NIFTP_fvPTC <- subset(Xhat,type=="NIFTP"|type=="fvPTC")

Xpca_again_2 <- prcomp(Xhat[,c(1:21254)], center=TRUE, scale=TRUE, retx=TRUE)
type_2 <- as.factor(data_NIFTP_fvPTC$type)
fviz_eig(Xpca_again_2)


fviz_pca_ind(Xpca_again_2,
             col.ind = type_2, 
             legend.title = "Groups",
             repel = TRUE
             )


