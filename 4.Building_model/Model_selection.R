# loading data 
library(data.table)
path = "/home/minhhoang/Documents/Thyroid cancer/TOTAL/Build_model"

traindata1 <- as.data.frame(fread(paste0(path,"/Data/traindata_model.tsv"), header=T))

testdata1 <- as.data.frame(fread(paste0(path,"/Data/testdata_model.tsv"), header=T))

#### Model 1: Using RandomForest algorithm
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
save(traindata_rf,file = paste0(path,"/Result/Random_forest/RF_model.RData"))

#Caculation parameters
testPred <- predict(traindata_rf, newdata=testdata1)
confu_matrix <- confusionMatrix(testPred, as.factor(testdata1$subtype), mode = "everything")

overal  <- data.frame(confu_matrix$overall)
class   <- data.frame(t(confu_matrix$byClass))
write.table(overal, paste0(path,"/Result/Random_forest/overal.tsv"), sep = "\t", row.names = F, col.names = T, quote = F)
write.table(class, paste0(path,"/Result/Random_forest/class.tsv"), sep = "\t", row.names = F, col.names = T, quote = F)


# Caculation averages AUC 
testPred <- predict(traindata_rf, newdata=testdata1, type="prob")
multi_AUC <- multiclass.roc(as.factor(testdata1$subtype), testPred)
save(multi_AUC,file = paste0(path,"/Result/Random_forest/AUC.RData"))


# predict, output is probability 
testPred_rf <- as.data.frame(predict(traindata_rf, newdata=testdata1, type="prob"))
colnames(testPred_rf) <- paste(colnames(testPred_rf), "_pred_RandomForest")

# Create true label 
true_label <- dummy_cols(testdata1, select_columns = 'subtype')
true_label <- true_label[,c(15:18)]
colnames(true_label) <- c("FA _true","FC _true", "fvPTC _true", "NIFTP _true")
final_df <- cbind(true_label, testPred_rf)

### draw 
#ROC/AUC 
roc_res <- multi_roc(final_df, force_diag=T)
save(roc_res,file = paste0(path,"/Result/Random_forest/roc_res.RData"))
unlist(roc_res$AUC)

plot_roc_df <- as.data.frame(plot_roc_data(roc_res))
plot_roc_df_1 <- plot_roc_df[c(1:120),]
plot_roc_df_1$anno <-  c(rep("FA (AUC=0.86)",30), rep("FC (AUC=0.65)",30),rep("fvPTC (AUC=0.82)",30),rep("NIFTP (AUC=1.00)",30))
                              
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
       subtitle = "Accuracy : 0.6897")
ggsave(paste0(path,"/Result/Random_forest/Rf_ROC_AUC.pdf"), width = 6, height = 6)

# PR curve 
pr_res <- multi_pr(final_df, force_diag=T)
save(pr_res,file = paste0(path,"/Result/Random_forest/pr_res.RData"))
unlist(pr_res$AUC)


plot_pr_df <- plot_pr_data(pr_res)
plot_pr_df_1 <- plot_pr_df[c(1:120),]
plot_pr_df_1$anno <-  c(rep("FA (AUC=0.64)",30), rep("FC (AUC=0.35)",30),rep("fvPTC (AUC=0.68)",30),rep("NIFTP (AUC=1.00)",30))
ggplot(plot_pr_df_1, aes(x=Recall, y=Precision)) + 
  geom_path(aes(color = anno, linetype = anno), size=1.5) +
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), 
                        colour='black', linetype = 'dotdash') +
  theme_bw() +
  theme(legend.justification=c(1, 0), legend.position=c(.75, .05),
                 legend.title=element_blank(), 
                 legend.background = element_rect(fill=NULL, size=0.5, 
                                                           linetype="longdash", colour ="black"))+
  theme(text = element_text(size = 18)) +
  labs(title = "RANDOMFOREST",
       subtitle = "Accuracy : 0.6897")
ggsave(paste0(path,"/Result/Random_forest/Rf_ROC_AUPRC.pdf"), width = 6, height = 6)


#######################
### model 2:  Using Boosting algorithm
library(gbm)
set.seed(123)
#creat model 
boost.boston<-gbm(subtype~., data = traindata1, distribution ="multinomial",n.trees=1000)
save(boost.boston,file = paste0(path,"/Result/Boosting/boosting_model.RData"))

# Predict 
hat.boost<-predict(boost.boston, newdata = testdata1, type="response")
hat.boost <- as.data.frame(hat.boost)
names(hat.boost) <-c("FA","FC","fvPTC","NIFTP")
predict <- names(hat.boost)[1:4][apply(hat.boost[,1:4], 1, which.max)]

# creat confusion matrix 
B_confusion_matrix <- confusionMatrix(as.factor(predict), as.factor(testdata1$subtype), mode = "everything")
overal  <- data.frame(B_confusion_matrix$overall)
class   <- data.frame(t(B_confusion_matrix$byClass))
write.table(overal, paste0(path,"/Result/Boosting/overal.tsv"), sep = "\t", row.names = F, col.names = T, quote = F)
write.table(class, paste0(path,"/Result/Boosting/class.tsv"), sep = "\t", row.names = F, col.names = T, quote = F)


# Caculation averages AUC 
multi_AUC <- multiclass.roc(as.factor(testdata1$subtype), hat.boost)
save(multi_AUC,file = paste0(path,"/Result/Boosting/AUC.RData"))


# Create true label 
hat.boost<-predict(boost.boston, newdata = testdata1, type="response")
colnames(hat.boost) <- paste(colnames(hat.boost), "_pred_Bososting")
true_label <- dummy_cols(testdata1, select_columns = 'subtype')
true_label <- true_label[,c(15:18)]
colnames(true_label) <- c("FA _true","FC _true", "fvPTC _true", "NIFTP _true")
final_df <- cbind(true_label, hat.boost)


#ROC/AUC 
roc_res <- multi_roc(final_df, force_diag=T)
save(roc_res,file = paste0(path,"/Result/Boosting/roc_res.RData"))
unlist(roc_res$AUC)


plot_roc_df <- as.data.frame(plot_roc_data(roc_res))
plot_roc_df_1 <- plot_roc_df[c(1:120),]
plot_roc_df_1$anno <-  c(rep("FA (AUC=0.78)",30), rep("FC (AUC= 0.66)",30),rep("fvPTC (AUC=0.81)",30),rep("NIFTP (AUC=1.00)",30))

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
  labs(title = "BOOSTING",
      subtitle ="Accuracy : 0.7241")
ggsave(paste0(path,"/Result/Boosting/Boosting _ROC_AUC.pdf"), width = 6, height = 6)


#multi pr curve 
pr_res <- multi_pr(final_df, force_diag=T)
save(pr_res,file = paste0(path,"/Result/Boosting/pr_res.RData"))
unlist(pr_res$AUC)

plot_pr_df <- plot_pr_data(pr_res)
plot_pr_df_1 <- plot_pr_df[c(1:120),]
plot_pr_df_1$anno <-  c(rep("FA (AUC=0.72)",30), rep("FC (AUC=0.36)",30),rep("fvPTC (AUC=0.57)",30),rep("NIFTP (AUC=1.00)",30))
ggplot(plot_pr_df_1, aes(x=Recall, y=Precision)) + 
  geom_path(aes(color = anno, linetype = anno), size=1.5) +
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), 
                        colour='black', linetype = 'dotdash') +
  theme_bw() +
  theme(legend.justification=c(1, 0), legend.position=c(.75, .05),
                 legend.title=element_blank(), 
                 legend.background = element_rect(fill=NULL, size=0.5, 
                                                           linetype="longdash", colour ="black"))+
  theme(text = element_text(size = 18)) +
  labs(title = "BOOSTING",
  subtitle ="Accuracy : 0.7241")
ggsave(paste0(path,"/Result/Boosting/Boosting_PRcruve.pdf"), width = 6, height = 6)

#############
# model 3: decision tree
library(tree)
library(rpart)
library(rpart.plot)

# creat model 
tree_1 <- rpart(subtype~.,data=traindata1, method="class")
save(tree_1,file = paste0(path,"/Result/Decision_tree/Decisiontree_model.RData"))

# creat plot decision tree
pdf(file = paste0(path,"/Result/Decision_tree/plot_tree.pdf"), width = 8, height = 6)
rpart.plot(tree_1 , extra = 0)
dev.off()

#creat confusion matrix 
tree.pred <- predict(tree_1, testdata1,type = "class")
decisiontree_matrix <- confusionMatrix(tree.pred, as.factor(testdata1$subtype), mode="everything")

overal  <- data.frame(decisiontree_matrix$overall)
class   <- data.frame(t(decisiontree_matrix$byClass))
write.table(overal, paste0(path,"/Result/Decision_tree/overal.tsv"), sep = "\t", row.names = F, col.names = T, quote = F)
write.table(class, paste0(path,"/Result/Decision_tree/class.tsv"), sep = "\t", row.names = F, col.names = T, quote = F)


# Caculation averages AUC 
testPred <- predict(tree_1, newdata=testdata1, type="prob")
AUC_dt <- multiclass.roc(as.factor(testdata1$subtype), testPred)
save(AUC_dt,file = paste0(path,"/Result/Decision_tree/AUC.RData"))


testPred_dt  <- as.data.frame(predict(tree_1, newdata=testdata1, type="prob"))
colnames(testPred_dt) <- paste(colnames(testPred_dt), "_pred_Decision_Tree")

# Create true label 
testPred_dt  <- as.data.frame(predict(tree_1, newdata=testdata1, type="prob"))
colnames(testPred_dt) <- paste(colnames(testPred_dt), "_pred_Decision_Tree")
true_label <- dummy_cols(testdata1, select_columns = 'subtype')
true_label <- true_label[,c(15:18)]
colnames(true_label) <- c("FA _true","FC _true", "fvPTC _true", "NIFTP _true")
final_df <- cbind(true_label, testPred_dt)

# Multi ROC/AUC
roc_res <- multi_roc(final_df, force_diag=T)
save(roc_res,file = paste0(path,"/Result/Decision_tree/roc_res.RData"))
unlist(roc_res$AUC)

plot_roc_df <- as.data.frame(plot_roc_data(roc_res))
plot_roc_df_1 <- plot_roc_df[c(1:120),]
plot_roc_df_1$anno <-  c(rep("FA (AUC=0.65)",30), rep("FC (AUC=0.48)",30),rep("fvPTC (AUC=0.60)",30),rep("NIFTP (AUC=0.87)",30))

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
  labs(title = "DECISION TREE",
       subtitle = "Accuracy : 0.4138")
ggsave(paste0(path,"/Result/Decision_tree/Decision_tree_ROC_AUC.pdf"), width = 6, height = 6)

# PR curve 
pr_res <- multi_pr(final_df, force_diag=T)
save(pr_res,file = paste0(path,"/Result/Decision_tree/pr_res.RData"))
unlist(pr_res$AUC)

plot_pr_df <- plot_pr_data(pr_res)
plot_pr_df_1 <- plot_pr_df[c(1:120),]
plot_pr_df_1$anno <-  c(rep("FA (AUC=0.50)",30), rep("FC (AUC=0.38)",30),rep("fvPTC (AUC=0.33)",30),rep("NIFTP (AUC=59)",30))
ggplot(plot_pr_df_1, aes(x=Recall, y=Precision)) + 
  geom_path(aes(color = anno, linetype = anno), size=1.5) +
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), 
                        colour='black', linetype = 'dotdash') +
  theme_bw() +
  theme(legend.justification=c(1, 0), legend.position=c(.95, .75),
                 legend.title=element_blank(), 
                 legend.background = element_rect(fill=NULL, size=0.5, 
                                                           linetype="longdash", colour ="black"))+
  theme(text = element_text(size = 18)) +
  labs(title = "DECISION TREE",
       subtitle ="Accuracy : 0.4138")
ggsave(paste0(path,"/Result/Decision_tree/Decision_tree_PRcruve.pdf"), width = 6, height = 6)

#####################
# model 4: SVM 
library(e1071)
library(data.table)
library(fastDummies)
library(ggplot2)
library(multiROC)
library(caret)
library(pROC)
library(GGally)

# data input  
x=traindata1[,1:13]
y=as.factor(traindata1$subtype)

#creat model 
model <- svm(x, y, probability = TRUE)
save(model,file = paste0(path,"/Result/SVM/SVM_model.RData"))
summary(model)

# Creat confusion matrix 
testPred <- predict(model, testdata1[,1:13])
svm_matrix  <- confusionMatrix(testPred, as.factor(testdata1$subtype), mode = "everything")

overal  <- data.frame(svm_matrix$overall)
class   <- data.frame(t(svm_matrix$byClass))
write.table(overal, paste0(path,"/Result/SVM/overal.tsv"), sep = "\t", row.names = F, col.names = T, quote = F)
write.table(class, paste0(path,"/Result/SVM/class.tsv"), sep = "\t", row.names = F, col.names = T, quote = F)

# Average AUC 
pred_prob <- predict(model, testdata1[,1:13], decision.values = TRUE, probability = TRUE)
AUC <- multiclass.roc(as.factor(testdata1$subtype), attr(pred_prob, "probabilities"))
save(AUC,file = paste0(path,"/Result/SVM/AUC.RData"))


testPred_rf <- as.data.frame(attr(pred_prob, "probabilities"))
colnames(testPred_rf) <- paste(colnames(testPred_rf), "_pred_SVM")

# Create true label 
testPred_rf <- as.data.frame(attr(pred_prob, "probabilities"))
colnames(testPred_rf) <- paste(colnames(testPred_rf), "_pred_SVM")
true_label <- dummy_cols(testdata1, select_columns = 'subtype')
true_label <- true_label[,c(15:18)]
colnames(true_label) <- c("FA _true","FC _true", "fvPTC _true", "NIFTP _true")
final_df <- cbind(true_label, testPred_rf)

### draw 
#ROC/AUC 
roc_res <- multi_roc(final_df, force_diag=T)
save(roc_res, file = paste0(path,"/Result/SVM/roc_res.RData"))
unlist(roc_res$AUC)
 
plot_roc_df <- as.data.frame(plot_roc_data(roc_res))
plot_roc_df_1 <- plot_roc_df[c(1:120),]
plot_roc_df_1$anno <-  c(rep("FA (AUC=0.78)",30), rep("FC (AUC=0.67)",30),rep("fvPTC (AUC=0.81)",30),rep("NIFTP (AUC=1.00)",30))
                              
ggplot(plot_roc_df_1, aes(x = 1-Specificity, y=Sensitivity)) +
  geom_path(aes(color = anno, linetype = anno), size=1.5) +
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), 
                        colour='black', linetype = 'dotdash') +
  theme_bw() +
  theme(legend.justification=c(1, 0), legend.position=c(.95, .05),
                 legend.title=element_blank(), 
                 legend.background = element_rect(fill=NULL, size=0.5, 
                                                           linetype="longdash", colour ="black"))+
  theme(text = element_text(size = 18))+
  labs(title = "Support Vector Machines (SVM)",
       subtitle = "Accuracy : 0.6552")
ggsave(paste0(path,"/Result/SVM/SVM_ROC_AUC.pdf"), width = 6, height = 6)

# PR curve 
pr_res <- multi_pr(final_df, force_diag=T)
save(pr_res, file = paste0(path,"/Result/SVM/pr_res.RData"))
unlist(pr_res$AUC)
  
plot_pr_df <- plot_pr_data(pr_res)
plot_pr_df_1 <- plot_pr_df[c(1:120),]
plot_pr_df_1$anno <-  c(rep("FA (AUC=0.57)",30), rep("FC (AUC=0.41)",30),rep("fvPTC (AUC=0.67)",30),rep("NIFTP (AUC=1.00)",30))
ggplot(plot_pr_df_1, aes(x=Recall, y=Precision)) + 
  geom_path(aes(color = anno, linetype = anno), size=1.5) +
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), 
                        colour='black', linetype = 'dotdash') +
  theme_bw() +
  theme(legend.justification=c(1, 0), legend.position=c(.75, .05),
                 legend.title=element_blank(), 
                 legend.background = element_rect(fill=NULL, size=0.5, 
                                                           linetype="longdash", colour ="black"))+
  theme(text = element_text(size = 18))+
  labs(title = "Support Vector Machines (SVM)",
       subtitle ="Accuracy : 0.6552")
ggsave(paste0(path,"/Result/SVM/SVM_PRcruve.pdf"), width = 6, height = 6)

