#https://cran.r-project.org/web/packages/ggalluvial/vignettes/ggalluvial.html
#Alluvial plot 
# Loading packages 
library(GGally)
library(ggalluvial)
library(ggplot2)
library(data.table)

#creat path 
path = "/home/minhhoang/Documents/Thyroid cancer/TOTAL/Build_model"

# load data 
data <- as.data.frame(fread(paste0(path, "/Result/Result_predict.tsv"),header=T))

#edit data 
Samples <- rep(seq(1,29,1),2)
Diagnostics  <- c(rep("Label",29),rep("Predict",29))
Subtype <- c(data$label,data$predict)
table <- data.frame(Samples, Diagnostics, Subtype)


#plot 
ggplot(table ,
       aes(x = Diagnostics, stratum = Subtype, alluvium = Samples,
           fill = Subtype, label = Subtype)) +
  scale_fill_brewer(type = "qual", palette = "Set2") +
  geom_flow(stat = "alluvium", lode.guidance = "frontback",
            color = "darkgray") +
  geom_stratum() +
  theme(text = element_text(size = 17),legend.position = "none") +
  geom_text(stat = "flow",
            aes(label = after_stat(n),
                hjust = (after_stat(flow) == "to")),nudge_x=0.22)+ geom_text(stat = "stratum", size = 7)+
    ggtitle("Original label cases and establishment of new diagnosis")

ggsave(paste0(path,"/Result/Images/Rf_Validation_testdata.pdf"), width = 8, height = 11)
