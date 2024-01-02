# cross_validation 
# number of subset probe 
probe_muvr_pls_9  <- c("cg01796223", "cg01091565",  "cg04705866", "cg11484872" , "cg26607785", "cg03409548", "cg06454226",  "cg19282250","cg03001305","subtype")
probe_muvr_pls_13 <- c("cg01796223",  "cg01091565", "cg04705866", "cg11484872", "cg26607785", "cg03409548","cg06454226","cg03001305",  "cg19282250", "cg15639951" ,"cg02192520" ,"cg08675585" ,"cg07447773","subtype")
probe_muvr_pls_20 <- c("cg01796223", "cg01091565", "cg04705866", "cg11484872", "cg26607785", "cg03409548", "cg06454226", "cg03001305", "cg19282250", "cg15639951", "cg02192520", "cg08675585", "cg10335112", "cg07447773" , "cg02719245", "cg13994177", "cg01344171", "cg04523589", "cg15888522", "cg08388746", "subtype")
probe_muvr_rf_9   <- c("cg01796223", "cg26607785", "cg11484872" , "cg15639951", "cg07447773", "cg03001305", "cg10335112", "cg02192520", "cg06454226","subtype")
probe_muvr_rf_13  <- c("cg01796223", "cg26607785", "cg11484872", "cg15639951", "cg07447773" , "cg03001305", "cg10335112", "cg02192520", "cg06454226", "cg19282250", "cg13994177", "cg22601215", "cg21098323","subtype")
probe_muvr_rf_20  <- c("cg01796223", "cg26607785", "cg11484872" , "cg15639951", "cg07447773", "cg03001305", "cg10335112", "cg02192520", "cg06454226", "cg19282250", "cg22601215", "cg13994177" ,"cg21098323", "cg08675585", "cg15782391", "cg04705866",  "cg24236938",  "cg10618882",  "cg14911395", "cg14258236","subtype")
probe_glmnet_16   <- c("cg10335112" ,"cg18575221" , "cg13994177" ,"cg01091565" ,"cg04705866",  "cg11484872" , "cg24236938" , "cg01344171" , "cg03001305" ,"cg15639951", "cg01796223",  "cg07447773" , "cg26607785" , "cg06454226" ,"cg19282250" , "cg02192520","subtype")
probe_feature_imp_13 <- c("cg01796223", "cg26607785", "cg07447773", "cg03001305", "cg15639951" ,"cg11484872", "cg02192520", "cg06454226", "cg10335112", "cg19282250","cg14911395", "cg19289461", "cg04705866","subtype")
probes_delta_fdr_18 <- c("cg20392764","cg01091565","cg02633817","cg09119967","cg14911395","cg21880903","cg14258236","cg11484872","cg18049750","cg08695223","cg03001305","cg01796223","cg07447773","cg26607785","cg06454226","cg19282250","cg02192520","subtype")

# creat list 
subset_probe <- list(probe_muvr_pls_9=probe_muvr_pls_9,probe_muvr_pls_13=probe_muvr_pls_13, probe_muvr_pls_20=probe_muvr_pls_20,
                        probe_muvr_rf_9=probe_muvr_rf_9,probe_muvr_rf_13=probe_muvr_rf_13,probe_muvr_rf_20=probe_muvr_rf_20,
                        probe_glmnet_16=probe_glmnet_16,
                        probe_feature_imp_13=probe_feature_imp_13,
                        probes_delta_fdr_18=probes_delta_fdr_18)

# data input 
library(data.table)
path = "/home/minhhoang/Documents/Thyroid cancer/TOTAL"
traindata <- as.data.frame(fread(paste0(path,'/Calling_DMP/First_calling_DMP/Data/traindata.tsv'),header=T))
traindata1 <- subset(traindata,subtype=="FA"|subtype=="fvPTC"|subtype=="FC"|subtype=="NIFTP")

# devind into groups 
set.seed(123)
# train: 70% ; test 30 % 
B = 100
sample_NIFTP <- c(1:4)
sample_FA <- c(5:11,23:31,44:54)
sample_FC <- c(12:18,32:37,55:69)
sample_fvPTC <- c(19:22,38:43,70:72)
test  = list()
for (k in 1:B ){
    a <- sample(sample_NIFTP,1,replace=F)
    b <- sample(sample_FA,8, replace=F)
    c <- sample(sample_FC,8, replace=F)
    d <- sample(sample_fvPTC,4, replace=F)
    e <- c(a,b,c,d )
    test[[k]] <- e
}


# result table 
res_all = data.frame(subset = character(),
                        AUC_average=numeric() ,
                        Accuracy_average=numeric())


for (i in 1:9) {
    print(names(subset_probe[i]))
    i_subset <- subset_probe[[i]]
    i_subset_names <- names(subset_probe[i])
    data <- traindata[,i_subset]
    for (j in test) {
        testdata <- data[j,]
        traindata <- data[-j,] 
        library(data.table)
        library(fastDummies)    
        library(ggplot2)
        library(multiROC)
        library(caret)
        library(pROC)
        library(GGally)
        library(randomForest)
        set.seed(123)
        traindata_rf <- randomForest(as.factor(subtype)~., data=traindata, ntree=1000, proximity=TRUE)
        testPred_1 <- predict(traindata_rf, newdata=testdata)
        Accuracy <- as.data.frame(confusionMatrix(testPred_1, as.factor(testdata$subtype), mode = "everything")$overal)[1,]
        testPred_2 <- predict(traindata_rf, newdata=testdata, type="prob")
        AUC <- as.data.frame(multiclass.roc(as.factor(testdata$subtype), testPred_2)$auc)[1,]
        #print(paste("AUC="," ",AUC,"Accuracy"," ",Accuracy))
        
    } 
    AUC_average <- mean(AUC)
    Accuracy_average <- mean(Accuracy)
    res_all[i,"AUC_average"] <- AUC_average
    res_all[i,"Accuracy_average"] <- Accuracy_average
    res_all[i,"subset"] <- i_subset_names
    print(paste0("AUC=",AUC_average, " ", "accuracy=",Accuracy_average))       
}

res_all
write.table(res_all, paste0(path, '/Feature_selection/Result/CV_result.tsv'), sep="\t", row.names=F, colnames=T, quote=T)


