# loading packages 
library(ggplot2)
library(ggpubr)
library(reshape2)
library(data.table)

# creat path 
path <- "/home/minhhoang/Documents/Thyroid cancer/TOTAL/Build_model"

# loading data  
traindata1 <- as.data.frame(fread(paste0(path,"/Data/traindata_model.tsv"), header=T))

# draw for 13 CpGSS 
probes <- colnames(traindata1)[-14]
for (i in 1:length(probes)) {
traindata2<- traindata1[,c(probes[i],"subtype")]
mtraindata2<- melt(traindata2)
ggplot(mtraindata2, aes(x = subtype, y = value, fill = subtype)) +
  geom_boxplot(width = 0.15, lwd = 1, outlier.color = NA) + theme_classic() +
  stat_compare_means(aes(group = subtype), bracket.size=10, size=7)+
  ggdist::stat_halfeye(aes(fill = subtype), adjust = .5, width = .3,
                       justification = -.6,.width = 0, point_colour = NA) + 
  gghalves::geom_half_point(aes(color = subtype),side = "l", range_scale = .3, size = 3, alpha = .6) +
  labs(fill = "",color = "", title=probes[1]) + xlab("") + ylab("Beta values") +
  theme(text = element_text(size = 25), legend.position="NONE")
ggsave(paste0(path,"/Result/Images", probes[i], ".pdf"), width = 11.00, height = 11.00)
}