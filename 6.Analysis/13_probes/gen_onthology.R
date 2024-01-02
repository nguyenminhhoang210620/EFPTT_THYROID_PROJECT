library(missMethyl)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(limma)
library(data.table)
set.seed(123)

#creat path 
path = "/home/minhhoang/Documents/Thyroid cancer/TOTAL/Analysis/13_probes"

# input data 
data <- as.data.frame(fread(paste0(path, '/Data/traindata_model.tsv'), header = T))


# Call CpGs 
cpgorders <- colnames(data)[1:13]

# creat annot 
ann <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

# All CpG sites tested
allcpgs <- rownames(ann)

# GO testing ignoring bias
gst.bias <- gometh(sig.cpg = cpgorders, all.cpg = allcpgs, collection = "GO", 
                   array.type = "EPIC", prior.prob=FALSE, anno = ann)

table(gst.bias$P.DE<0.05)
table(gst.bias$FDR<0.05)

# filter with P.DE < 0.05
topgst.bias <- gst.bias[gst.bias$P.DE<0.05,]
dim(topgst.bias)
write.table(topgst.bias, paste0(path,"/Result/gen_ontholygy/ontology_GO.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)


# Significant genes to KEGG (Kyoto Encyclopedia of Genes and Genomes) output
kegg.siggenes <- gometh(sig.cpg = cpgorders, all.cpg = allcpgs, 
                        array.type = "EPIC", collection = "KEGG", anno = ann, sig.genes = TRUE)

table(kegg.siggenes$P.DE<0.05)
table(kegg.siggenes$FDR<0.05)

# FILTER WITH P.DE 
topkeggsig <- kegg.siggenes[kegg.siggenes$P.DE<0.05,]
dim(topkeggsig)
write.table(topkeggsig, paste0(path,"/Result/gen_ontholygy/ontology_KEGG.tsv.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)


################### - VISULIZATION GENE ONTHOLOGY  - #######################
### GO 
# N: number of genes in the GO or KEGG term
#DE: number of genes that are differentially methylated
#Ont:ontology that the GO term belongs to if testing GO pathways. "BP" - biological process, "CC" - cellular component, "MF" - molecular function.

# loading packages 
library(data.table)
library(ggplot2)
library(ggpubr)

# loadind data 
GO <- as.data.frame(fread(paste0(path, "/Result/gen_ontholygy/ontology_GO.tsv")))
table(GO$ONTOLOGY)
# BP  CC  MF 
#782   6  39 

BP <- GO[GO$ONTOLOGY=="BP",]
CC <- GO[GO$ONTOLOGY=="CC",]
MF <- GO[GO$ONTOLOGY=="MF",]

# save file 
BP <- BP[order(BP$FDR,decreasing=FALSE),]
CC <- CC[order(CC$FDR,decreasing=FALSE),]
MF <- MF[order(MF$FDR,decreasing=FALSE),]
write.table(BP, paste0(path,"/Result/gen_ontholygy/ontology_GO_BP.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
write.table(CC, paste0(path,"/Result/gen_ontholygy/ontology_GO_CC.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)
write.table(MF, paste0(path,"/Result/gen_ontholygy/ontology_GO_MF.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)