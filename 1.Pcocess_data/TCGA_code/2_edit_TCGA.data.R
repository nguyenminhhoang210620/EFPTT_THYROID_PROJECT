# dowload data from gadi
path="/g/data3/yo4/Cancer-Epigenetics/davban/Hoang/TCGA"
write.table(data, paste0(path,"/TCGA_data"), sep="\t", col.names=T, quote=F, row.names=F )
#
library(data.table)
data <- as.data.frame(fread('/home/minhhoang/Documents/Thyroid cancer/TCGA.DATA/TCGA_data', header=T))
tdata <- as.data.frame(t(data[-108]))
colnames(tdata) <- data$probeID 
# 13 probes## probe <- c("cg01796223",  "cg01091565", "cg04705866", "cg11484872", "cg26607785", "cg03409548","cg06454226","cg03001305",  "cg19282250", "cg15639951" ,"cg02192520" ,"cg08675585" ,"cg07447773")
# 8 probe obe### probe <- c("cg26607785","cg20387341","cg26521404","cg07447773","cg24623694","cg21280510","cg02510853","cg22467534") : FC 
# 4 probes ## probe <- c("cg18265887","cg08657449","cg02192520","cg17243643")
tdata <- (tdata[,c(probe)])
names <-  rownames(tdata)
tdata$name <- names 
colnames(tdata)[5] <- "id"

# info 
info <- as.data.frame(fread("/home/minhhoang/Documents/Thyroid cancer/data/TCGA.DATA/clinical.tsv",header=T))
# check 
info$submitter_id==tdata$id

# dignosis 
library(readxl)
final <- as.data.frame(read_excel('/home/minhhoang/Documents/Thyroid cancer/Data and analysis/one_again/non_removebatchefect/probe_differ_tumour_normal/data_without_normal/TCGA/Final.tsv'))
final <- final[,c(2,7)]
colnames(final) <- c("id", "final")

# merge 
data_TCGA <- merge(tdata, final, by="id")
data_TCGA <- data_TCGA[,-1]
# delete NIFTP and FA: have 29 cases
data_TCGA_1 <- data_TCGA[!(data_TCGA$final=="NIFTP"|data_TCGA$final=="FA"),]

# chuyển đổi đáu chấm thành dấu phẩy 
data_TCGA_1$cg01796223 <- gsub(",",".",data_TCGA_1$cg01796223)
data_TCGA_1$cg01091565 <- gsub(",",".",data_TCGA_1$cg01091565)
data_TCGA_1$cg04705866 <- gsub(",",".",data_TCGA_1$cg04705866)
data_TCGA_1$cg11484872 <- gsub(",",".",data_TCGA_1$cg11484872)
data_TCGA_1$cg26607785 <- gsub(",",".",data_TCGA_1$cg26607785)
data_TCGA_1$cg03409548 <- gsub(",",".",data_TCGA_1$cg03409548)
data_TCGA_1$cg06454226 <- gsub(",",".",data_TCGA_1$cg06454226)
data_TCGA_1$cg03001305 <- gsub(",",".",data_TCGA_1$cg03001305)
data_TCGA_1$cg19282250 <- gsub(",",".",data_TCGA_1$cg19282250)
data_TCGA_1$cg15639951 <- gsub(",",".",data_TCGA_1$cg15639951)
data_TCGA_1$cg02192520 <- gsub(",",".",data_TCGA_1$cg02192520)
data_TCGA_1$cg08675585 <- gsub(",",".",data_TCGA_1$cg08675585)
data_TCGA_1$cg07447773 <- gsub(",",".",data_TCGA_1$cg07447773) 

# fvPTC 
data_TCGA_1$cg18265887 <- gsub(",",".",data_TCGA_1$cg18265887)
data_TCGA_1$cg08657449 <- gsub(",",".",data_TCGA_1$cg08657449)
data_TCGA_1$cg02192520 <- gsub(",",".",data_TCGA_1$cg02192520)
data_TCGA_1$cg17243643 <- gsub(",",".",data_TCGA_1$cg17243643)



path = "/home/minhhoang/Documents/Thyroid cancer/Data and analysis/one_again/non_removebatchefect/Train_data/fvPTC-NON-fvPTC"
write.table(data_TCGA_1, paste0(path,"/TCGA_4PROBE.tsv"), sep="\t", col.names=T, quote=F, row.names=F )


# creat model to predict  with 13 probes in MUVR pPLS 
library(data.table)
traindata <- as.data.frame(fread("/home/minhhoang/Documents/Thyroid cancer/Data and analysis/one_again/non_removebatchefect/Train_data/traindata.tsv",header=T))
traindata1 <- traindata[,c(probe,"subtype")]
traindata1 <- subset(traindata1,subtype=="FA"|subtype=="fvPTC"|subtype=="FC"|subtype=="NIFTP")
traindata1$group <- ifelse(traindata1$subtype=="fvPTC","fvPTC","NON-fvPTC")
traindata1 <- traindata1[,-5]

# testdata = TCGA 
library(data.table)
testdata1 <- as.data.frame(fread('/home/minhhoang/Documents/Thyroid cancer/data/TCGA.DATA/TCGA_13PROBE.tsv',header=T))
colnames(testdata1)[5] <- "subtype"
testdata1$group <- ifelse(testdata1$subtype=="fvPTC","fvPTC","NON-fvPTC")
testdata1 <- testdata1[,-5]

# model Randomforesr  
library(data.table)
library(fastDummies)
library(ggplot2)
library(multiROC)
library(caret)
library(pROC)
library(GGally)
library(randomForest)
set.seed(123)
# Creat model 
traindata_rf <- randomForest(as.factor(subtype)~., data=traindata1, ntree=1000, proximity=TRUE)
#Call:
# randomForest(formula = as.factor(subtype) ~ ., data = traindata1,      ntree = 1000, proximity = TRUE) 
#               Type of random forest: classification
#                     Number of trees: 1000
#No. of variables tried at each split: 3

#        OOB estimate of  error rate: 22.22%
#Confusion matrix:
#      FA FC fvPTC NIFTP class.error
#FA    26  0     1     0  0.03703704
#FC     2 22     4     0  0.21428571
#fvPTC  5  4     4     0  0.69230769
#NIFTP  0  0     0     4  0.00000000

#Caculation parameters
testPred <- predict(traindata_rf, newdata=testdata1[1:13])
confusionMatrix(testPred, as.factor(testdata1$subtype), mode = "everything") # different levels 
table(testPred, as.factor(testdata1$subtype))
#testPred FC fvPTC
#   FA     2     1
#   FC     3     3
#   fvPTC  0    20
#   NIFTP  0     0

#caculation average AUC 
testPred <- predict(traindata_rf, newdata=testdata1, type="prob")
multiclass.roc(as.factor(testdata1$subtype), testPred)
#Call:
#multiclass.roc.default(response = as.factor(testdata1$subtype),     predictor = testPred)

#Data: multivariate predictor testPred with 2 levels of as.factor(testdata1$subtype): FC, fvPTC.
#Multi-class area under the curve: 0.75

# predict, output is probability 
testPred_rf <- as.data.frame(predict(traindata_rf, newdata=testdata1, type="prob"))
colnames(testPred_rf) <- paste(colnames(testPred_rf), "_pred_RandomForest")

# Create true label 
true_label <- dummy_cols(testdata1, select_columns = 'subtype')
true_label <- true_label[,c(15:16)]
colnames(true_label) <- c("FC _true", "fvPTC _true")
final_df <- cbind(true_label, testPred_rf)

### draw 
#ROC/AUC 
path="/home/minhhoang/Documents/Thyroid cancer/data/TCGA.DATA"
roc_res <- multi_roc(final_df, force_diag=T)
unlist(roc_res$AUC)
#RandomForest.FC  RandomForest.fvPTC   RandomForest.macro  RandomForest.micro 
#          0.5500000           0.9500000           0.7291667           0.8299643 
 
plot_roc_df <- as.data.frame(plot_roc_data(roc_res))
plot_roc_df_1 <- plot_roc_df[c(1:60),]
plot_roc_df_1$anno <-  c(rep("FC (AUC=0.55)",30),rep("fvPTC (AUC=0.95)",30))
                              
ggplot(plot_roc_df_1, aes(x = 1-Specificity, y=Sensitivity)) +
  geom_path(aes(color = anno, linetype = anno), size=1.5) +
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), 
                        colour='black', linetype = 'dotdash') +
  theme_bw() +
  theme(legend.justification=c(1, 0), legend.position=c(.95, .05),
                 legend.title=element_blank(), 
                 legend.background = element_rect(fill=NULL, size=0.5, 
                                                           linetype="longdash", colour ="black"))+
  theme(text = element_text(size = 18)) +
  labs(title = "RANDOMFOREST",
       subtitle = "Accuracy : 0.80")
ggsave(paste0(path,"/Rf_ROC_AUC_validation.pdf"), width = 11, height = 11)
# PR curve 
pr_res <- multi_pr(final_df, force_diag=T)
unlist(pr_res$AUC)
#RandomForest.FC  RandomForest.fvPTC   RandomForest.macro  RandomForest.micro 
#          0.2502157           0.9897157           0.6119580           0.8528491 
plot_pr_df <- plot_pr_data(pr_res)
plot_pr_df_1 <- plot_pr_df[c(1:60),]
plot_pr_df_1$anno <-  c(rep("FC (AUC=0.25)",30),rep("fvPTC (AUC=0.99)",30))
ggplot(plot_pr_df_1, aes(x=Recall, y=Precision)) + 
  geom_path(aes(color = anno), size=2) + 
  theme_bw() + 
  theme(legend.justification=c(1, 0), legend.position=c(.75, .00),
        legend.title=element_blank(), 
        legend.background = element_rect(fill=NULL, size=0.5, linetype="solid", colour ="black")) +
labs(title = "RANDOMFOREST",
       subtitle ="Accuracy : 0.80")
ggsave(paste0(path,"/Rf_PRcruve_validation.pdf"), width = 11, height = 11)

# creat table predict 
table <- data.frame(re_dig=testdata1$subtype ,predict=testPred)
write.table(table, paste0(path,"/table_predict.tsv"), sep="\t", col.names=T, quote=F, row.names=F )

#Aluvia plot 
library(GGally)
library(ggalluvial)
library(ggplot2)
library(data.table)


data <- as.data.frame(fread("/home/minhhoang/Documents/Thyroid cancer/TCGA.DATA/table_predict.tsv",header=T))
new_NIFTP <- c("NIFTP","NIFTP")
new_NIFTP_1 <- c("NIFTP","NIFTP")
new_FA <- c("FA","FA")
new_FA_1 <- c("FA","FA")
data_new <- rbind(data,new_NIFTP,new_NIFTP_1,new_FA,new_FA_1)

#edit data 
Samples <- rep(seq(1,33,1),2)
Diagnostics  <- c(rep("Label",33),rep("Predict",33))
Subtype <- c(data_new$re_dig ,data_new$predict)
table <- data.frame(Samples, Diagnostics, Subtype)


ggplot(table ,
       aes(x = Diagnostics, stratum = Subtype, alluvium = Samples,
           fill = Subtype, label = Subtype)) +
  scale_fill_brewer(type = "qual", palette = "Set2") +
  geom_flow(stat = "alluvium", lode.guidance = "frontback",
            color = "darkgray") +
  geom_stratum() +
  theme(text = element_text(size = 17), legend.position = "none") +
  #geom_text(stat = "flow",
  #          aes(label = after_stat(n),
  #              hjust = (after_stat(flow) == "to")),nudge_x=0.25)+ 
  geom_text(stat = "stratum", size = 6)+
    ggtitle("Original label cases and establishment of new diagnosis")
path="/home/minhhoang/Documents/Thyroid cancer/TCGA.DATA"    
ggsave(paste0(path,"/Rf_Validation_TCGA_1.pdf"), width = 8, height = 11)


ggplot(b, aes(x=Var1,y=Freq, fill=Gender))+ 
    geom_bar(stat="identity", position=position_dodge())+
    labs(x="Subtype", y="Frequency", title="Gender of Validation Data")+
    theme(text=element_text(size=20))+
    geom_text(aes(label=Freq), vjust=1.6,  color="white",
            position = position_dodge(0.9), size=5)
path="/home/minhhoang/Documents/Thyroid cancer/TCGA.DATA"    
ggsave(paste0(path,"/Gender.pdf"), width = 8, height = 11)



ggplot(a  , aes(x=age_at_index, color=Subtypes))+ 
    geom_density(size=2)+
    theme(text=element_text(size=20))+
    geom_vline(data=e, aes(xintercept=Mean, color=Subtype),
             linetype="dashed")+
    labs(title="Age distribution of validation data", x= "Age", y="Density")
path="/home/minhhoang/Documents/Thyroid cancer/TCGA.DATA"    
ggsave(paste0(path,"/Age.pdf"), width = 8, height = 11)
    
    