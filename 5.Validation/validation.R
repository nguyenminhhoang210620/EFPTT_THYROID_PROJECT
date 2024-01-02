####### VALIDATION ################
#Creat path 
path = '/home/minhhoang/Documents/Thyroid cancer/TOTAL/Validation'

# Loading two independent cohorts 
library(data.table)

#TCGA cohort 
TCGA <- as.data.frame(fread(paste0(path,"/Data/Validation_TCGA.tsv"),header=T))

#Juan cohort 
juan <- as.data.frame(fread(paste0(path,"/Data/validation_juan.tsv"),header=T))
juan <- juan[!juan$subtype=="NH",]

#input traning data 
traindata1 <- as.data.frame(fread(paste0(path,"/Data/traindata_model.tsv"),header=T))

#Loading packages 
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

####### validation with TCGA cohort ####################
testPred <- predict(traindata_rf , newdata=TCGA[,1:13])
TCGA_confusion_matrix <- as.matrix(table(testPred, TCGA$subtype))
write.table(TCGA_confusion_matrix, paste0(path,'/Result/TCGA/confusion_matrix.tsv'), sep = "\t", col.names = T, row.names = T, quote = F)


# Caculation averages AUC 
testPred <- predict(traindata_rf, newdata=TCGA[,1:13], type="prob")
AUC <- multiclass.roc(as.factor(TCGA$subtype), testPred)
saveRDS(AUC, file = paste0(path, '/Result/TCGA/AUC.RData'))

# predict, output is probability 
testPred_rf <- as.data.frame(predict(traindata_rf, newdata=TCGA[,1:13], type="prob"))
colnames(testPred_rf) <- paste(colnames(testPred_rf), "_pred_RandomForest")

# Create true label 
true_label <- dummy_cols(TCGA, select_columns = 'subtype')
true_label <- true_label[,c(15:16)]
colnames(true_label) <- c("FC _true", "fvPTC _true")
final_df <- cbind(true_label, testPred_rf)

### draw 
#ROC/AUC 
roc_res <- multi_roc(final_df, force_diag=T)
saveRDS(roc_res, file = paste0(path, '/Result/TCGA/roc_res.RData'))
unlist(roc_res$AUC)
RandomForest.FC  RandomForest.fvPTC   RandomForest.macro  RandomForest.micro 
          0.5500000           0.9500000           0.7291667           0.8299643 


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
  labs(title = "Validation_TCGA",
       subtitle = "Accuracy : 0.80")
ggsave(paste0(path,"/Result/TCGA/Rf_ROC_AUC.pdf"), width = 6, height = 6)

# PR curve 
pr_res <- multi_pr(final_df, force_diag=T)
saveRDS(pr_res, file = paste0(path, '/Result/TCGA/pr_res.RData'))
unlist(pr_res$AUC)
unlist(pr_res$AUC)
   RandomForest.FC  RandomForest.fvPTC   RandomForest.macro  RandomForest.micro 
          0.2502157           0.9897157           0.6119580           0.8528491 

plot_pr_df <- plot_pr_data(pr_res)
plot_pr_df_1 <- plot_pr_df[c(1:60),]
plot_pr_df_1$anno <-  c(rep("FC (AUC=0.25)",30),rep("fvPTC (AUC=0.99)",30))
ggplot(plot_pr_df_1, aes(x=Recall, y=Precision)) + 
  geom_path(aes(color = anno, linetype = anno), size=1.5) +
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), 
                        colour='black', linetype = 'dotdash') +
  theme_bw() +
  theme(legend.justification=c(1, 0), legend.position=c(.50, .05),
                 legend.title=element_blank(), 
                 legend.background = element_rect(fill=NULL, size=0.5, 
                                                           linetype="longdash", colour ="black"))+
  theme(text = element_text(size = 18)) +
  labs(title = "Validation_TCGA",
       subtitle = "Accuracy : 0.80")
ggsave(paste0(path,"/Result/TCGA/Rf_PRcruve.pdf"), width = 6, height = 6)


############ Validation juan cohort #########################
##### Try with 12 probes, cg11484872 was remove #####

set.seed(123)
# remove CpGs
traindata_12 <- traindata1[,-4]

# Creat model 
traindata_rf <- randomForest(as.factor(subtype)~., data=traindata_12, ntree=1000, proximity=TRUE)
saveRDS(traindata_rf, file = paste0(path, '/Result/juan/RF_model_12_probes.RData'))

#confusion matrix 
testPred <- predict(traindata_rf , newdata=juan[1:12])
juan_confusion_matrix <- as.matrix(table(testPred, juan$subtype))
write.table(juan_confusion_matrix, paste0(path,'/Result/juan/juan_confusion_matrix.tsv'), sep = "\t", col.names = T, row.names = T, quote = F)

# Caculation averages AUC 
testPred <- predict(traindata_rf, newdata=juan[,1:12], type="prob")
AUC <- multiclass.roc(as.factor(juan$subtype), testPred)
saveRDS(AUC, file = paste0(path, '/Result/juan/AUC.RData'))



# Create true label 
testPred_rf <- as.data.frame(predict(traindata_rf, newdata=juan, type="prob"))
colnames(testPred_rf) <- paste(colnames(testPred_rf), "_pred_RandomForest")
true_label <- dummy_cols(juan, select_columns = 'subtype')
true_label <- true_label[,c(15:17)]
colnames(true_label) <- c("FA _true","FC _true", "fvPTC _true")
final_df <- cbind(true_label, testPred_rf)

### draw 
#ROC/AUC 
roc_res <- multi_roc(final_df, force_diag=T)
saveRDS(roc_res, file = paste0(path, '/Result/juan/roc_res.RData'))
unlist(roc_res$AUC)
RandomForest.FA     RandomForest.FC  RandomForest.fvPTC   RandomForest.macro 
          0.8313725           0.7227273           0.8114286           0.7861459 
 RandomForest.micro 
          0.8007812 

plot_roc_df <- as.data.frame(plot_roc_data(roc_res))
plot_roc_df_1 <- plot_roc_df[c(1:99),]
plot_roc_df_1$anno <-  c(rep("FA (AUC=0.83)",33), rep("FC (AUC=0.72)",33),rep("fvPTC (AUC=0.81)",33))
                              
                              
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
  labs(title = "Validation_juan",
       subtitle = "Accuracy : 0.5625")
  ggsave(paste0(path,"/Result/juan/Rf_ROC_AUC.pdf"), width = 6, height = 6)


# PR curve 
pr_res <- multi_pr(final_df, force_diag=T)
saveRDS(, file = paste0(path, '/Result/juan/pr_res.RData'))
unlist(pr_res$AUC)

plot_pr_df <- plot_pr_data(pr_res)
plot_pr_df_1 <- plot_pr_df[c(1:99),]
plot_pr_df_1$anno <-  c(rep("FA (AUC=0.77)",33), rep("FC (AUC=0.59)",33),rep("fvPTC (AUC=0.61)",33))


ggplot(plot_pr_df_1, aes(x=Recall, y=Precision)) + 
  geom_path(aes(color = anno, linetype = anno), size=1.5) +
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), 
                        colour='black', linetype = 'dotdash') +
  theme_bw() +
  theme(legend.justification=c(1, 0), legend.position=c(.90, .05),
                 legend.title=element_blank(), 
                 legend.background = element_rect(fill=NULL, size=0.5, 
                                                           linetype="longdash", colour ="black"))+
  theme(text = element_text(size = 18)) +
  labs(title = "Validation_juan",
       subtitle = "Accuracy : 0.5625")
ggsave(paste0(path,"/Result/juan/Rf_PRcruve.pdf"), width = 6, height = 6)



