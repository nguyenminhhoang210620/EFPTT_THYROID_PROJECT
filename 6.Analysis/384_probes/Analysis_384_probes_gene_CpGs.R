##### Gene and CpGs analysis ############### 
library(data.table)
library(GenomicRanges)
library(wesanderson)
library(ggplot2)
library(data.table)
library(dplyr)
library(annotatr)
library(GenomicRanges)
library(annotatr)
library(Gviz)

# create path 
path <- "/home/minhhoang/Documents/Thyroid cancer/TOTAL/Analysis/384_probes"

# creat list files 
file <- c("FA_normal_probes.tsv", "FC_normal_probes.tsv", "fvPTC_normal_probes.tsv", "NIFTP_normal_probes.tsv")
folder <- c("FA", "FC","fvPTC", "NIFTP")

#craet annot 
annots = c('hg19_cpgs', 'hg19_basicgenes', 'hg19_genes_intergenic','hg19_genes_intronexonboundaries')
annotations = build_annotations(genome ='hg19', annotations = annots)


##### run annotation #######
for (i in 1:4){
    #loading files 
    data_probes <- as.data.frame(fread(paste0(path,"/Data/",file[i]), header = T))

    # transfrom 
    data_transform  <- makeGRangesFromDataFrame(data_probes , keep.extra.columns=TRUE)

    # run 
    dm_annotated = annotate_regions(regions = data_transform, annotations = annotations, ignore.strand = FALSE, quiet = FALSE)

    # as.data.frame 
    df_dm_annotated = data.frame(dm_annotated)

    # save files 
    write.table(df_dm_annotated, paste0(path,"/Result/",folder[i], "/", file[i]),sep = "\t", col.names = T, quote = F, row.names = F)
}

##### Creat file and exact following status of DMP ######
# creat list file form result above 
file_result  <- c("FA_normal_probes.tsv", "FC_normal_probes.tsv", "fvPTC_normal_probes.tsv", "NIFTP_normal_probes.tsv")

for (i in 1:4){
    print(file_result[i])
    # loading data 
    data <- as.data.frame(fread(paste0(path, "/Result/",folder[i], "/",  file_result[i]), header = T))

    # annnot gene location 
    data_1 <- data[grep("hg19_genes_", data$annot.type),]
    data_1$annot.type2 <- "Gene body"
    data_1["annot.type2"][data_1["annot.type"] == "hg19_genes_promoters"] <- "Promoter"
    data_1["annot.type2"][data_1["annot.type"] == "hg19_genes_intergenic"] <- "Intergenic"
    data_1["annot.type2"][data_1["annot.type"] == "hg19_genes_1to5kb"] <- "1-5kb gene"
    data_1["annot.type2"][data_1["annot.type"] == "hg19_genes_3UTRs"] <- "3UTR gene"
    data_1["annot.type2"][data_1["annot.type"] == "hg19_genes_5UTRs"] <- "5UTR gene"
    data_1["annot.type2"][data_1["annot.type"] == "hg19_genes_exons"] <- "Exon"
    data_1["annot.type2"][data_1["annot.type"] == "hg19_genes_intronexonboundaries"] <- "Intron Exon Boundary"
    data_1["annot.type2"][data_1["annot.type"] == "hg19_genes_introns"] <- "Intron"

    # exact with status (gene) 
    data_hyper_gene <- subset(data_1, data_1$status=="Hyper DMPs")
    data_hypo_gene <- subset(data_1, data_1$status=="Hypo DMPs")
    write.table(data_hyper_gene, paste0(path,"/Result/",folder[i], "/", "gene_hyper_", file_result[i]), sep = "\t", col.names = T, row.names = F, quote = F)
    write.table(data_hypo_gene, paste0(path,"/Result/",folder[i], "/", "gene_hypo_", file_result[i]), sep = "\t", col.names = T, row.names = F, quote = F)


    # annot CpGs location 
    data_2 <- data[grep("hg19_cpg_", data$annot.type),]
    data_2$annot.type <- gsub("hg19_cpg_", "CpG ", data_2$annot.type)


    # exact with status (CpGs)
    data_hyper_CpG <- subset(data_2, data_2$status=="Hyper DMPs")
    data_hypo_CpG <- subset(data_2, data_2$status=="Hypo DMPs")
    write.table(data_hyper_CpG, paste0(path,"/Result/",folder[i], "/", "CpGs_hyper_", file_result[i]), sep = "\t", col.names = T, row.names = F, quote = F)
    write.table(data_hypo_CpG, paste0(path,"/Result/", folder[i], "/","CpGs_hypo_", file_result[i]), sep = "\t", col.names = T, row.names = F, quote = F)
 }   


####### Creat backound from 21253 CpGs ###########
# loaling data 
df_dm_annotated <- as.data.frame(fread(paste0(path, "/Data/background_annotation.tsv"), header=T))

# exact background CpGs 
annottype.bg1 <- df_dm_annotated[grep("hg19_cpg_", df_dm_annotated$annot.type),]
annottype.bg1$annot.type <- gsub("hg19_cpg_", "CpG ", annottype.bg1$annot.type)

# save file 
write.table(annottype.bg1, paste0(path,"/Result/CpGs_background.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)

# Gene 
annottype.bg2 <- df_dm_annotated[grep("hg19_genes_", df_dm_annotated$annot.type),]
annottype.bg2$annot.type2 <- "Gene body"
annottype.bg2["annot.type2"][annottype.bg2["annot.type"] == "hg19_genes_promoters"] <- "Promoter"
annottype.bg2["annot.type2"][annottype.bg2["annot.type"] == "hg19_genes_intergenic"] <- "Intergenic"
annottype.bg2["annot.type2"][annottype.bg2["annot.type"] == "hg19_genes_1to5kb"] <- "1-5kb gene"
annottype.bg2["annot.type2"][annottype.bg2["annot.type"] == "hg19_genes_3UTRs"] <- "3UTR gene"
annottype.bg2["annot.type2"][annottype.bg2["annot.type"] == "hg19_genes_5UTRs"] <- "5UTR gene"
annottype.bg2["annot.type2"][annottype.bg2["annot.type"] == "hg19_genes_exons"] <- "Exon"
annottype.bg2["annot.type2"][annottype.bg2["annot.type"] == "hg19_genes_intronexonboundaries"] <- "Intron Exon Boundary"
annottype.bg2["annot.type2"][annottype.bg2["annot.type"] == "hg19_genes_introns"] <- "Intron"

# save file 
write.table(annottype.bg2, paste0(path,"/Result/Gene_background.tsv"), sep = "\t", col.names = T, row.names = F, quote = F)



############### Draw barplot ###################
# Gene location 
# creat list 
# loading packgaes 
library(dplyr)

for (i in 1:4) {
    print(folder[i])
    # gene background 
    background <- as.data.frame(fread(paste0(path, "/Result/Gene_background.tsv"), header = T))
    background$main <- ifelse(background$annot.type2=="Promoter","Promoter",ifelse(background$annot.type2=="Intergenic","Intergenic", "Gene Body"))
    background_gene <- as.data.frame(background %>% 
                        group_by(main) %>% 
                        summarise( percent = 100 * n() / nrow(background)))
    background_gene$diff <- "Genomic\nBackgorund"

    # hyper CpGs 
    tryCatch({
    hyper <- as.data.frame(fread(paste0(path,"/Result/", folder[i], "/gene_hyper_", folder[i], "_normal_probes.tsv"),header=T))
    hyper$main <- ifelse(hyper$annot.type2=="Promoter","Promoter",ifelse(hyper$annot.type2=="Intergenic","Intergenic", "Gene Body"))
    hyper_gene  <- as.data.frame(hyper %>% 
                        group_by(main) %>% 
                        summarise( percent = 100 * n() / nrow(hyper)))
    hyper_gene$diff <- "Hyper DMPs"
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

    # Hypo CpGs 
    tryCatch({
    hypo <- as.data.frame(fread(paste0(path,"/Result/", folder[i], "/gene_hypo_", folder[i], "_normal_probes.tsv"), header=T))
    hypo$main <- ifelse(hypo$annot.type2=="Promoter","Promoter",ifelse(hypo$annot.type2=="Intergenic","Intergenic", "Gene Body"))
    hypo_gene  <- as.data.frame(hypo %>% 
                        group_by(main) %>% 
                        summarise( percent = 100 * n() / nrow(hypo)))
    hypo_gene$diff <- "Hypo DMPs"
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

    # combind data 
    all.gene <- rbind(background_gene, hyper_gene, hypo_gene)

    #draw 
    tryCatch({
    library(ggplot2)
    ggplot(all.gene, aes(x=diff,y=percent, fill=main))+
            geom_bar(stat="identity", colour="black", size=2)+
            scale_fill_brewer(palette = "RdGy") +
            labs(x = NULL, y = "Percent (%)", fill = NULL, title= paste0(folder[i], "-", "Subtype")) +
            theme(text = element_text(size = 30)) +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_line(colour = "black"))

    ggsave(paste0(path, "/Result/", folder[i], "/", "gene_anotaion", ".pdf"), width = 11, height = 8)
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

}

------------------------------------------
# CpGs location 
# creat list 
# loading packgaes 
library(dplyr)

for (i in 1:4) {
    print(folder[i])
    # gene background 
    background <- as.data.frame(fread(paste0(path, "/Result/CpGs_background.tsv"), header = T))
    background$main <- ifelse(background$annot.type=="CpG islands","CpG islands",ifelse(background$annot.type=="CpG shores","CpG shores", "Non CpG Island"))
    background_CpG <- as.data.frame(background %>% 
                        group_by(main) %>% 
                        summarise( percent = 100 * n() / nrow(background)))
    background_CpG$diff <- "Genomic\nBackgorund"

    # hyper CpGs 
    tryCatch({
    hyper <- as.data.frame(fread(paste0(path,"/Result/", folder[i], "/CpGs_hyper_", folder[i], "_normal_probes.tsv"),header=T))
    hyper$main <- ifelse(hyper$annot.type=="CpG islands","CpG islands",ifelse(hyper$annot.type=="CpG shores","CpG shores", "Non CpG Island"))
    hyper_CpG  <- as.data.frame(hyper %>% 
                        group_by(main) %>% 
                        summarise( percent = 100 * n() / nrow(hyper)))
    hyper_CpG$diff <- "Hyper DMPs"
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})


    # Hypo CpGs 
    tryCatch({
    hypo <- as.data.frame(fread(paste0(path,"/Result/", folder[i], "/CpGs_hypo_", folder[i], "_normal_probes.tsv"), header=T))
    hypo$main <- ifelse(hypo$annot.type=="CpG islands","CpG islands",ifelse(hypo$annot.type=="CpG shores","CpG shores", "Non CpG Island"))
    hypo_CpG  <- as.data.frame(hypo %>% 
                        group_by(main) %>% 
                        summarise( percent = 100 * n() / nrow(hypo)))
    hypo_CpG$diff <- "Hypo DMPs"
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

    # combind data 
    all.CpG <- rbind(background_CpG, hyper_CpG, hypo_CpG)

    #draw 
    tryCatch({
    library(ggplot2)
    ggplot(all.CpG, aes(x=diff,y=percent, fill=main))+
            geom_bar(stat="identity", colour="black", size=2)+
            scale_fill_brewer(palette = "BuPu") +
            labs(x = NULL, y = "Percent (%)", fill = NULL, title= paste0(folder[i], "-", "Subtype")) +
            theme(text = element_text(size = 30)) +
            theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_line(colour = "black"))

    ggsave(paste0(path, "/Result/", folder[i], "/", "CpG_anotaion", ".pdf"), width = 11, height = 8)
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

}
