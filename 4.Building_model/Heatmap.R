############ CREAT HEATMAP ############
# Loading packages 
library(circlize)
library(ggplot2)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(data.table)

# crear path 
path = "/home/minhhoang/Documents/Thyroid cancer/TOTAL/Build_model"

#Loading data 
traindata1 <- as.data.frame(fread(paste0(path,"/Data/traindata_model.tsv")))

# edit data 
traindata1 <- traindata1[order(factor(traindata1$subtype, levels=c("FA" ,"FC" ,"fvPTC","NIFTP"))),]
hdata <- as.matrix(t(traindata1[,-14]))
colnames(hdata) <- traindata1$subtype 
dim(hdata)


# creat color 
col_fun = colorRamp2(c(0, 0.5, 1), c("#377EB8", "white", "#E41A1C"))
ha = HeatmapAnnotation(Subtypes = factor(traindata1$subtype, levels=c("FA" ,"FC" ,"fvPTC","NIFTP")),
                        col = list(Subtypes = c("FA" = "pink", "FC" = "#FFFFB3", "fvPTC" = "#377EB8", "NIFTP" = "green")),
                        show_legend = FALSE, simple_anno_size = unit(0.5, "cm"))


pdf(file = paste0(path, "/Result/Images/heatmap_13probes.pdf"), width = 11, height = 11)
Heatmap(hdata, name = "Methylation", col= col_fun,
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  column_split = factor(traindata1$subtype, levels = c("FA" ,"FC" ,"fvPTC","NIFTP")),
  row_km = 2,
  top_annotation = ha,
  show_column_names = FALSE,
  show_row_names = FALSE,
  show_row_dend = TRUE ,
  row_title = "DNA Methylation Probes (n=13)",
  column_title_side = "top",
  column_title_gp = gpar(fontsize = 16, fontface = "bold"),   
  row_title_gp = gpar(fontsize = 20, fontface = "bold"),
  row_names_gp = gpar(fontsize = 10),
  heatmap_legend_param = list(legend_height = unit(5, "cm"),title_position = "lefttop-rot"),
  column_gap = unit(2, "mm"),
  row_gap = unit(2, "mm"),
  border = TRUE)
dev.off()