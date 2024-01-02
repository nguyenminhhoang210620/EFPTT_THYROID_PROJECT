###### FEATURES SELECTION ###############
#### Using MUVR packages 
## Method: RF 
# loading packages 
library(MUVR)
library(data.table)
library(doParallel)
set.seed(123)

#creat path 
path <- "/home/minhhoang/Documents/Thyroid cancer/TOTAL/Feature_selection"
# loading data 
traindata <- as.data.frame(fread(paste0(path, "/Data/79_probes.tsv"),header=T))

#Create input for MUVR
y <- traindata$subtype
# Response
x <- traindata[,c(1:79)]

### RF (10 repetition)
nCore = detectCores()-1

# Number of MUVR repetitions 
nRep = 100 

# Number of outer cross-validation segments          
nOuter = 4    

# Proportion of variables kept per iteration            
varRatio = 0.8 

# Selected core modelling algorithm        
method = 'RF'    

# Set up parallel processing using doParallel
cl=makeCluster(nCore)   
registerDoParallel(cl)

# Perform modelling
classModel <- MUVR(X=x, Y=y, nRep=nRep, nOuter=nOuter, varRatio=varRatio, method=method)

classModel$miss
#min mid max 
#17  17  18 
classModel$nVar
#min mid max 
#  9  13  20

## model min: 9 probes 
min <- getVIP(classModel, model='min')
#         order       name      rank
#cg01796223     1 cg01796223  2.160000
#cg26607785     2 cg26607785  4.236667
#cg11484872     3 cg11484872 12.509167
#cg15639951     4 cg15639951 20.245000
#cg07447773     5 cg07447773 21.529167
#cg03001305     6 cg03001305 26.527500
#cg10335112     7 cg10335112 30.890833
#cg02192520     8 cg02192520 35.161667
#cg06454226     9 cg06454226 38.599167
probe_min_RF <- c("cg01796223", "cg26607785", "cg11484872" , "cg15639951", "cg07447773", "cg03001305", "cg10335112", "cg02192520", "cg06454226") 
saveRDS(probe_min_RF, paste0(path, "/Result/MUVR/probes_min_RF.rds"))

## model mid: 13 probes 
mid <- getVIP(classModel, model='mid')
#           order       name      rank
#cg01796223     1 cg01796223  1.860833
#cg26607785     2 cg26607785  4.200000
#cg11484872     3 cg11484872 10.457500
#cg15639951     4 cg15639951 15.545833
#cg07447773     5 cg07447773 18.510000
#cg03001305     6 cg03001305 20.869167
#cg10335112     7 cg10335112 26.409167
#cg02192520     8 cg02192520 29.108333
#cg06454226     9 cg06454226 31.510833
#cg19282250    10 cg19282250 41.549167
#cg13994177    11 cg13994177 51.170000
#cg22601215    12 cg22601215 51.905833
#cg21098323    13 cg21098323 54.816667
probe_mid_RF <- c("cg01796223", "cg26607785", "cg11484872", "cg15639951", "cg07447773" , "cg03001305", "cg10335112", "cg02192520", "cg06454226", "cg19282250", "cg13994177", "cg22601215", "cg21098323") 
saveRDS(probe_mid_RF, paste0(path, "/Result/MUVR/probes_mid_RF.rds"))

# model max: 20 probes 
max <- getVIP(classModel, model='max')
#       order       name      rank
#cg01796223     1 cg01796223  1.839167
#cg26607785     2 cg26607785  4.315833
#cg11484872     3 cg11484872 10.405833
#cg15639951     4 cg15639951 14.239167
#cg07447773     5 cg07447773 16.963333
#cg03001305     6 cg03001305 19.086667
#cg10335112     7 cg10335112 23.789167
#cg02192520     8 cg02192520 27.124167
#cg06454226     9 cg06454226 29.085000
#cg19282250    10 cg19282250 37.266667
#cg22601215    11 cg22601215 48.786667
#cg13994177    12 cg13994177 48.793333
#cg21098323    13 cg21098323 51.605833
#cg08675585    14 cg08675585 55.534167
#cg15782391    15 cg15782391 57.565000
#cg04705866    16 cg04705866 57.887500
#cg24236938    17 cg24236938 59.152500
#cg10618882    18 cg10618882 60.770833
#cg14911395    19 cg14911395 61.740833
#cg14258236    20 cg14258236 61.854167
probe_max_RF <- c("cg01796223", "cg26607785", "cg11484872" ,"cg15639951", "cg07447773", "cg03001305", "cg10335112", "cg02192520", "cg06454226", "cg19282250", "cg22601215", "cg13994177" ,"cg21098323", "cg08675585", "cg15782391", "cg04705866",  "cg24236938",  "cg10618882",  "cg14911395", "cg14258236")     
saveRDS(probe_max_RF, paste0(path, "/Result/MUVR/probes_max_RF.rds"))

#########################################
## Method: RF 

library(doParallel)
library(MUVR)
set.seed(123)
#data input 
#Create input for MUVR
y <- traindata$subtype
# Response
x <- traindata[,c(1:79)]

#Create input for MUVR
y <- traindata$subtype
# Response
x <- traindata[,c(1:79)]

### RF (10 repetition)
nCore = detectCores()-1

# Number of MUVR repetitions 
nRep = 100 

# Number of outer cross-validation segments          
nOuter = 4    

# Proportion of variables kept per iteration            
varRatio = 0.8 

# Selected core modelling algorithm        
method = 'PLS'    

# Set up parallel processing using doParallel
cl=makeCluster(nCore)   
registerDoParallel(cl)

# Perform modelling
classModel = MUVR(X=x, Y=y, nRep=nRep, nOuter=nOuter, varRatio=varRatio, method=method)

classModel$miss
#min mid max 
# 21  21  20
classModel$nVar
#min mid max 
#  9  13  20 

## model min: 9 probes 
min <- getVIP(classModel, model='min')
#            order       name      rank
#cg01796223     1 cg01796223  2.713333
#cg01091565     2 cg01091565  4.143333
#cg04705866     3 cg04705866  9.363333
#cg11484872     4 cg11484872 14.771667
#cg26607785     5 cg26607785 18.764167
#cg03409548     6 cg03409548 22.449167
#cg06454226     7 cg06454226 22.650000
#cg19282250     8 cg19282250 27.120000
#cg03001305     9 cg03001305 27.706667
probe_min_pls <- c("cg01796223", "cg01091565",  "cg04705866", "cg11484872" , "cg26607785", "cg03409548", "cg06454226",  "cg19282250","cg03001305") 
saveRDS(probe_min_pls, paste0(path, "/Result/MUVR/probes_min_pls.rds"))

## model mid: 13 probes 
mid  <- getVIP(classModel, model='mid')
#      order       name      rank
#cg01796223     1 cg01796223  2.823333
#cg01091565     2 cg01091565  4.535000
#cg04705866     3 cg04705866  8.656667
#cg11484872     4 cg11484872 12.122500
#cg26607785     5 cg26607785 15.660000
#cg03409548     6 cg03409548 16.639167
#cg06454226     7 cg06454226 18.187500
#cg03001305     8 cg03001305 22.242500
#cg19282250     9 cg19282250 22.958333
#cg15639951    10 cg15639951 24.288333
#cg02192520    11 cg02192520 27.537500
#cg08675585    12 cg08675585 30.230000
#cg07447773    13 cg07447773 36.248333
probe_mid_pls <- c("cg01796223",  "cg01091565", "cg04705866", "cg11484872", "cg26607785", "cg03409548","cg06454226","cg03001305",  "cg19282250", "cg15639951" ,"cg02192520" ,"cg08675585" ,"cg07447773")
saveRDS(probe_mid_pls, paste0(path, "/Result/MUVR/probes_mid_pls.rds"))

## model max: 20 probes 
max <- getVIP(classModel, model='max')
#            order       name      rank
#cg01796223     1 cg01796223  2.676667
#cg01091565     2 cg01091565  5.281667
#cg04705866     3 cg04705866  9.245833
#cg11484872     4 cg11484872 12.461667
#cg26607785     5 cg26607785 14.864167
#cg03409548     6 cg03409548 15.687500
#cg06454226     7 cg06454226 17.690000
#cg03001305     8 cg03001305 20.307500
#cg19282250     9 cg19282250 21.144167
#cg15639951    10 cg15639951 21.380833
#cg02192520    11 cg02192520 26.501667
#cg08675585    12 cg08675585 27.229167
#cg10335112    13 cg10335112 32.856667
#cg07447773    14 cg07447773 33.326667
#cg02719245    15 cg02719245 48.575833
#cg13994177    16 cg13994177 48.771667
#cg01344171    17 cg01344171 51.506667
#cg04523589    18 cg04523589 51.834167
#cg15888522    19 cg15888522 54.689167
#cg08388746    20 cg08388746 58.454167
probe_max_pls <- c("cg01796223", "cg01091565", "cg04705866", "cg11484872", "cg26607785", "cg03409548", "cg06454226", "cg03001305", "cg19282250", "cg15639951", "cg02192520", "cg08675585", "cg10335112", "cg07447773" , "cg02719245", "cg13994177", "cg01344171", "cg04523589", "cg15888522", "cg08388746")
saveRDS(probe_max_pls, paste0(path, "/Result/MUVR/probes_max_pls.rds"))

#########################
#### glmnet packaegs ##### 

## loading packages 
library(glmnet)
library(data.table)
set.seed(123)

#creat path 
path <- "/home/minhhoang/Documents/Thyroid cancer/TOTAL/Feature_selection"
# loading data 
traindata <- as.data.frame(fread(paste0(path, "/Data/79_probes.tsv"),header=T))

# input data 
x <- as.matrix(traindata[,c(1:79)])
y <- as.factor(traindata$subtype)

# creat model 
fit <- glmnet(x, y, family = "multinomial", type.multinomial = "grouped")

pdf(file = paste0(path, "/Result/glmnet/glmnet_plot_fit.pdf"), width = 4, height = 4)
plot(fit, xvar = "lambda", label = TRUE, type.coef = "2norm")
dev.off()


# creat model with Cross-validaition 
cvfit <- cv.glmnet(x, y, family = "multinomial", type.multinomial = "grouped")

pdf(file = paste0(path, "/Result/glmnet/glmnet_plot_cvfit.pdf"), width = 4, height = 4)
plot(cvfit)
dev.off()

### exact probes with lamda min 
# lamda min 
cvfit$lambda.min
#0.07595046

A  <- coef(cvfit, s = "lambda.min")
#FA 
B <- as.data.frame(as.matrix(A$FA))
B$probeID <- rownames(B)
colnames(B)[1] <- "lamda" 
B <- B[!(B$lamda==0),]
B$probeID 
# 16 probe with lamda min 
#[1] "(Intercept)" "cg10335112"  "cg18575221"  "cg13994177"  "cg01091565" 
#[6] "cg04705866"  "cg11484872"  "cg24236938"  "cg01344171"  "cg03001305" 
#[11] "cg15639951"  "cg01796223"  "cg07447773"  "cg26607785"  "cg06454226" 
#[16] "cg19282250"  "cg02192520" 

#FC 
C <- as.data.frame(as.matrix(A$FC))
C$probeID <- rownames(C)
colnames(C)[1] <- "lamda" 
C <- C[!(C$lamda==0),]
C$probeID 
#16 probe with lamda min
#[1] "(Intercept)" "cg10335112"  "cg18575221"  "cg13994177"  "cg01091565" 
#[6] "cg04705866"  "cg11484872"  "cg24236938"  "cg01344171"  "cg03001305" 
#[11] "cg15639951"  "cg01796223"  "cg07447773"  "cg26607785"  "cg06454226" 
#[16] "cg19282250"  "cg02192520" 

#NIFTP 
D <- as.data.frame(as.matrix(A$NIFTP))
D$probeID <- rownames(D)
colnames(D)[1] <- "lamda" 
D <- D[!(D$lamda==0),] 
D$probeID
#16 probe with lamda min
#[1] "(Intercept)" "cg10335112"  "cg18575221"  "cg13994177"  "cg01091565" 
#[6] "cg04705866"  "cg11484872"  "cg24236938"  "cg01344171"  "cg03001305" 
#[11] "cg15639951"  "cg01796223"  "cg07447773"  "cg26607785"  "cg06454226" 
#[16] "cg19282250"  "cg02192520"

#fvPTC 
E <- as.data.frame(as.matrix(A$fvPTC))
E$probeID <- rownames(E)
colnames(E)[1] <- "lamda" 
E <- E[!(E$lamda==0),] 
E$probeID
#16 probe with lamda min
#[1] "(Intercept)" "cg10335112"  "cg18575221"  "cg13994177"  "cg01091565" 
#[6] "cg04705866"  "cg11484872"  "cg24236938"  "cg01344171"  "cg03001305" 
#[11] "cg15639951"  "cg01796223"  "cg07447773"  "cg26607785"  "cg06454226" 
#[16] "cg19282250"  "cg02192520" 


probe_glmnet <- c("cg10335112" ,"cg18575221" , "cg13994177" ,"cg01091565" ,"cg04705866",  "cg11484872" , "cg24236938" , "cg01344171" , "cg03001305" ,"cg15639951", "cg01796223",  "cg07447773" , "cg26607785" , "cg06454226" ,"cg19282250" , "cg02192520" )
saveRDS(probe_glmnet, paste0(path, "/Result/glmnet/probes_glmnet.rds"))

######################################################
## importance feature in Randomforest  ##
# loading packages 
library(data.table)
library(randomForest)
set.seed(123)

#creat path 
path <- "/home/minhhoang/Documents/Thyroid cancer/TOTAL/Feature_selection"
# loading data 
traindata <- as.data.frame(fread(paste0(path, "/Data/79_probes.tsv"),header=T))

# Creat model 
traindata_rf <- randomForest(as.factor(subtype)~., data=traindata, ntree=1000, proximity=TRUE)

#plot importance plot 
pdf(file = paste0(path, "/Result/Importance_feature_RF/glmnet_plot_fit.pdf"), width = 4, height = 4)
randomForest::varImpPlot(traindata_rf, sort=TRUE, main="Variable Importance Plot")
dev.off()


### exact 13 probes with highest giini index 
impprobe <- as.data.frame(importance(traindata_rf))
impprobe$probeID <- rownames(impprobe)
impprobe1 <- impprobe[order(impprobe$MeanDecreaseGini,decreasing = TRUE),]
probes <- impprobe1[c(1:13),2]
saveRDS(probes, paste0(path, "/Result/Importance_feature_RF/probes_importance_RF.rds"))