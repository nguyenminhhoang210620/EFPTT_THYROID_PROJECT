# Loading packages  
library(circlize)
library(ggplot2)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(data.table)

# create path 
path = "/home/minhhoang/Documents/Thyroid cancer/TOTAL/Analysis/384_probes"

# input_status 
status_DMP <- as.data.frame(fread(paste0(path, "/Data/384_probes_DMP_status.tsv"))) 

# input DNA methylation data 
traindata <- as.data.frame(fread(paste0(path, "/Data/traindata.tsv"), header = T))
traindata1 <- traindata[,c(status_DMP$probeID,"subtype")]
traindata1 <- traindata1[!(traindata1$subtype=="NG"),]
traindata1 <- traindata1[order(factor(traindata1$subtype, levels=c("FC" ,"fvPTC" ,"NIFTP","FA","normal"))),]
hdata <- as.matrix(t(traindata1[,-385]))
colnames(hdata) <- traindata1$subtype 
dim(hdata)
rownames(hdata) <- status_DMP$status


# creat annot 
ha1 = rowAnnotation(Status=factor(rownames(hdata),levels=c("Hyper DMPs","Hypo DMPs")),
                        col=list(Status= c("Hyper DMPs"="#E7298A", "Hypo DMPs" = "#377EB8")),
                        show_legend = TRUE, simple_anno_size = unit(0.7, "cm"))
ha = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = 2:9),
     labels = c("FC" ,"fvPTC" ,"NIFTP","FA","normal"), 
     labels_gp = gpar(col = "white", fontsize = 10, fontface = "bold")))

# creat annot with DNA mehtylaion level 
meth_col_fun = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))

# draw plot 
pdf(file = paste0(path, "/Result/Heatmap_384.pdf"), width = 15, height = 11)
Heatmap(hdata, name="Methylation", col= meth_col_fun,
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  column_split=factor(colnames(hdata),levels=c("FC" ,"fvPTC" ,"NIFTP","FA","normal")),
  show_column_names = FALSE,
  show_row_names = FALSE,
  show_column_dend=FALSE,
  show_row_dend=TRUE,
  top_annotation = ha,
  right_annotation = ha1,
  row_split=factor(rownames(hdata)),
  column_title = "DNA METHYLATION BETWEEN EFPTT AND NORMAL",
  column_title_gp = gpar(fontsize = 16, fontface = "bold",title_position="bottom"),
  heatmap_legend_param = list(legend_height = unit(5, "cm"),title_position = "lefttop-rot"),
  row_title="DNA Methylation Probes (n=384)")
dev.off()
