###### PLOT PCA ###########
# loading packages 
library(reshape2)
library(factoextra)
library(ggplot2)
library(dplyr)
library(ggforce)
library(data.table)

# creat path 
path = "/home/minhhoang/Documents/Thyroid cancer/TOTAL/Build_model"

#Loading data 
data1 <- as.data.frame(fread(paste0(path, "/Data/traindata_model.tsv")))

# transform data 
data2 <- as.data.frame(data1[,-14])

# Run PCA 
PC <- prcomp(data2, center = TRUE, scale = FALSE, retx = TRUE)

# plot 
pdf(file = paste0(path, '/Result/Images/fviz_plot.pdf'), width = 8, height = 5)
fviz_eig(PC,addlabels = TRUE, ncp=20)
dev.off()

# draw PCA 
loadings <- as.data.frame(PC$x)
loadings$subtype <- data1$subtype
loadings%>%
  ggplot(aes(x = PC1, y = PC2, color = subtype)) + geom_point(size=5) +
  xlab("PC1 (33.3%)") + ylab("PC2 (21.3%)") + labs(colour = "Subtype")+
  theme(legend.position="bottom") + theme(text = element_text(size = 20))+
  geom_text(aes(label = subtype ,size = 3, hjust=0, vjust=2),show.legend = FALSE)+
  stat_ellipse(geom="polygon", aes(fill = subtype),alpha = 0.2,
  show.legend = FALSE, level = 0.95)
ggsave(paste0(path,"/Result/Images/PCA_13_probes.pdf"), width = 11, height = 11)


######### t-SNE #############
# loading packages  
library(Rtsne)
library(tidyverse)
set.seed(123) 

# Run TSNE
train.sne <- Rtsne(data1[,-14],perplexity=round((nrow(data1)-2)/3))

# Show the objects in the 2D tsne representation
plot(train.sne$Y,col=factor(data1$subtype), asp=1)

# creat data.frame for ploting 
train.dt <- data.frame(train.sne$Y)
names(train.dt) <- c("dim1","dim2")
train.dt$subtype <- data1$subtype

train.dt %>%
  ggplot(aes(x = dim1, y = dim2, col= subtype))+ geom_point(size=5)+ 
  xlab("Dimension 1") + ylab("Dimension 2") + labs(col = "subtype") +
  theme(legend.position="bottom") + theme(text = element_text(size = 20))+
  geom_text(aes(label = subtype ,size = 3, hjust=0, vjust=2),show.legend = FALSE)+
  stat_ellipse(geom="polygon", aes(fill = subtype),alpha = 0.2,
  show.legend = FALSE, level = 0.95)
ggsave(paste0(path,"/Result/Images/t-SNE_13_probes.pdf"), width = 11, height = 11)

############## UMAP #######################
# Loading packages  
library(umap)
library(ggfortify)
library(ROCit)
library(ggpubr)
library(glmnet)
library(fivethirtyeight)
library(broom)
library(ggdist)

set.seed(123)

# Run UMAP  
train.umap <- umap(data1[,-14])

#plot(train.umap$layout,col=factor(data1$subtype), asp=1)

# creat data,frame for ploting 
train.dt <- data.frame(train.umap$layout)
names(train.dt) <- c("dim1","dim2")
train.dt$subtype <- data1$subtype

#ploting 
train.dt %>%
  ggplot(aes(x = dim1, y = dim2, color = subtype))+
  geom_point(size=5)+ xlab("Dimension 1") + ylab("Dimension 2") + labs(colour = "Class") +
  theme(legend.position="bottom") + theme(text = element_text(size = 20))+
  geom_text(aes(label = subtype ,size = 3, hjust=0, vjust=2),show.legend = FALSE)+
  stat_ellipse(geom="polygon", aes(fill = subtype),alpha = 0.2,
  show.legend = FALSE, level = 0.95)
ggsave(paste0(path,"/Result/Images/UMAP_13_probes.pdf"), width = 11, height = 11)

######### MDS plot ########
# Loading packages  
library(ggpubr)
library(glmnet)
library(fivethirtyeight)
library(broom)
library(ggdist)
set.seed(123)


#Run MDS 
dO <- dist(data1[,-14])
train.mds <- cmdscale(dO,eig=TRUE, k=2)

#plot(train.mds$points,col=factor(data1$subtype), asp=1)

#creat data.frame for ploting 
train.dt <- data.frame(train.mds$points)
names(train.dt) <- c("dim1","dim2")
train.dt$subtype <- data1$subtype

# PLotting 
train.dt %>%
  ggplot(aes(x = dim1, y = dim2, color = subtype))+
  geom_point(size=5)+ xlab("Dimension 1") + ylab("Dimension 2") + labs(colour = "Subtype") +
  theme(legend.position="bottom") + theme(text = element_text(size = 20))+
  geom_text(aes(label = subtype ,size = 3, hjust=0, vjust=2),show.legend = FALSE)+
  stat_ellipse(geom="polygon", aes(fill = subtype),alpha = 0.2,
  show.legend = FALSE, level = 0.95)
ggsave(paste0(path,"/Result/Images/MDS_13_probes.pdf"), width = 11, height = 11)