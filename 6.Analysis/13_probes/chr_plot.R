#Loadding packages 
library(RIdeogram)
library(ggplot2)
library(data.table)
library(stringr)

# creat path 
path = "/home/minhhoang/Documents/Thyroid cancer/TOTAL/Analysis/13_probes"

# input data 
position <- as.data.frame(fread(paste0(path,"/Data/DMPs_13_probes.tsv")))

# select chr 
human_karyotype2 <- human_karyotype[c(1,5,6,7,9,12,14,15,17,18),]

position <- position[,c(1,2,3)]
position2 <- position[!duplicated(position), ]

colnames(position2) <- c("Chr","Start","End" )
position2$Type <- rep('CpGs',13)
position2$Shape <- rep('circle',13)
position2$color <- rep('33a02c',13)
position3 <- position2[,c("Type","Shape", "Chr", "Start", "End", "color")]
position3$Chr <- as.character(position3$Chr)

chr <- numeric()
for (i in 1:13){
	a <- as.numeric(paste(strsplit(position3$Chr, split = "")[[i]][-c(1:3)], collapse = ''))
	chr[i] <- a 
}
position3$Chr <- chr


ideogram(karyotype = human_karyotype2, overlaid = gene_density, label = position3, label_type = "marker")
convertSVG("chromosome.svg", device = "png")


