library(data.table)
library(caret)
data <- as.data.frame(fread("/home/minhhoang/Documents/Thyroid cancer/Data and analysis/one_again/non_removebatchefect/bigtable_merge_19_normal_full.tsv",header=T))
info <- read.table("/home/minhhoang/Documents/Thyroid cancer/Data and analysis/one_again/non_removebatchefect/info_merge_27k_full_tsv", header=T, sep="\t", check.names=F, stringsAsFactors=F)
tdata <-as.data.frame(t(data[,-c(1:5)]))
colnames(tdata) <- data$probeID 
tdata$subtype <- info$subtype
dim(tdata)
tdata[1:5,1:5]
# devide data into test and train 
set.seed(123)
ind <- sample(2, nrow(tdata), replace=TRUE, prob=c(0.7,0.3))
table(ind)
traindata  <- tdata[ind==1,]
table(traindata[,21254])
#table(traindata[,21254])

    #FA     FC  fvPTC     NG  NIFTP normal 
    #27     28     13      4      4     55 
#dim(traindata)
#131 21254
traindata$title <- rownames(traindata)
write.table(traindata, "/home/minhhoang/Documents/Thyroid cancer/Data and analysis/one_again/non_removebatchefect/Train_data/traindata.tsv", quote=F, sep="\t", row.names=F)
testdata <- tdata[ind==2,]
table(testdata[,21254])
#table(testdata[,21254])

    #FA     FC  fvPTC     NG  NIFTP normal 
    #10      8      9      2      2     22 
#dim(testdata)
#53 21254
testdata$title <- rownames(testdata)
write.table(testdata, "/home/minhhoang/Documents/Thyroid cancer/Data and analysis/one_again/non_removebatchefect/Test_data/testdata.tsv", quote=F, sep="\t", row.names=F)

# devide info  into test and train 
# train info 
info_train <- subset(info, title%in%traindata$title)
write.table(info_train, "/home/minhhoang/Documents/Thyroid cancer/Data and analysis/one_again/non_removebatchefect/Train_data/info_train.tsv", quote=F, sep="\t", row.names=F)
info_test <- subset(info, title%in%testdata$title)
write.table(info_test, "/home/minhhoang/Documents/Thyroid cancer/Data and analysis/one_again/non_removebatchefect/Test_data/info_test.tsv", quote=F, sep="\t", row.names=F)