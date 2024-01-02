library(data.table)
library(dplyr)
########## -----GENE------ ############
# create path 
path = "/home/minhhoang/Documents/Thyroid cancer/TOTAL/Analysis/384_probes"

# Background 
background <- as.data.frame(fread(paste0(path, "/Result/Gene_background.tsv"), header = T))
background$main <- ifelse(background$annot.type2=="Promoter","Promoter",ifelse(background$annot.type2=="Intergenic","Intergenic", "Gene Body"))
#Gene Body Intergenic   Promoter 
#     40886         51      15100 
background_1 <- as.data.frame(table(background$main))
# other subtypes 

folder <- c("NIFTP", "FC", "fvPTC", "FA")
status <- c("hyper", "hypo")


for (i in 1:4){
    print(folder[i])
    for (j in 1:2){
        print(status[j])
        data <- as.data.frame(fread(paste0(path,"/Result/",folder[i], "/gene_", status[j], "_", folder[i], "_normal_probes.tsv"))) 
        data$main <- ifelse(data$annot.type2=="Promoter","Promoter",ifelse(data$annot.type2=="Intergenic","Intergenic", "Gene Body"))
        data_1 <- as.data.frame(table(data$main))
        target <- data_1$Freq 
        print(data_1)
    }
}

#Gene
FA <- c(14,0,2)
NIFTP_hyper <- c(19,0,7)
NIFTP_hypo <- c(302,3,106)
FC_hypo <- c(4,0,2)
FC_hyper <- c(211,0,86)
fvPTC_hyper <- c(8,0,0)
background <- c(40886,51,15100)

total <- cbind(background,FA,NIFTP_hyper,NIFTP_hypo,FC_hypo,FC_hyper,fvPTC_hyper)
rownames(total) <- background_1$Var1
tdt  <- as.data.frame(t(total))

res_table <- data.frame(
                        p_value1 = numeric(),
                        p_value2 = numeric(),
                        p_value3 = numeric(),
                        p_value4 = numeric(),
                        p_value5 = numeric(),
                        p_value6 = numeric(),
                        pair = character(),
                        stringsAsFactors=FALSE)

pair1 <- combn(unique(rownames(tdt)),2)
pair2 <- combn(unique(colnames(tdt)),2)


for (i in 1:6) {
    print(pair1[,i])
    sub <- pair1[,i]
    s1 <- sub[1]
    s2 <- sub[2] 
    colnames(res_table)[i] <- paste0(s1,"_",s2)
    tdte <- tdt[pair1[,i],]
    for (j in 1:ncol(tdte)) {
    	gr <- pair2[,j]
		g1 <- gr[1]
		g2 <- gr[2]
		res_table[j,7] <- paste0(g1,"_",g2)
		p.value <- fisher.test(tdte[,pair2[,j]])$p.value
	res_table[j,i] <- p.value    
	}
}

write.table(res_table,  paste0(path, "/Result/gene_test_result.tsv"), sep ="\t", quote=F, col.names=T, row.names=F)


res_table_f <- res_table[,c(7,1,2,3,4,5,6)]	

pair3 <- combn(unique(colnames(res_table_f)),2)

#FA 
dt <- res_table_f[,pair3[,1]]
dt <- arrange(dt, background_FA)
dt$adjust_p_value <- p.adjust(dt[,2], "BH")
# pair              background_FA  adjust_p_value
# Intron-Promoter  0.0006511933     0.01627983

dt <- res_table_f[,pair3[,2]]
dt <- arrange(dt, background_NIFTP_hyper)
dt$adjust_p_value <- p.adjust(dt[,2], "BH")
#                    pair background_NIFTP_hyper adjust_p_value
#Gene Body_Intergenic                      1              1
#Gene Body_Promoter                      1              1
#Intergenic_Promoter                      1              1
 

dt <- res_table_f[,pair3[,3]]
dt <- arrange(dt, background_NIFTP_hypo)
dt$adjust_p_value <- p.adjust(dt[,2], "BH")
#                 pair background_NIFTP_hypo adjust_p_value
#  Intergenic_Promoter           0.006750684     0.01128097
# Gene Body_Intergenic           0.007520647     0.01128097
#   Gene Body_Promoter           0.695206766     0.69520677



dt <- res_table_f[,pair3[,4]]
dt <- arrange(dt, background_FC_hypo)
dt$adjust_p_value <- p.adjust(dt[,2], "BH")
#                  pair background_FC_hypo adjust_p_value
#   Gene Body_Promoter          0.6638559              1
#   Gene Body_Intergenic          1.0000000              1
#   Intergenic_Promoter          1.0000000              1


dt <- res_table_f[,pair3[,5]]
dt <- arrange(dt, background_FC_hyper)
dt$adjust_p_value <- p.adjust(dt[,2], "BH")
#                   pair background_FC_hyper adjust_p_value
#   Gene Body_Promoter           0.4323151              1
#   Gene Body_Intergenic           1.0000000              1
#  Intergenic_Promoter           1.0000000              1


dt <- res_table_f[,pair3[,6]]
dt <- arrange(dt, background_fvPTC_hyper)
dt$adjust_p_value <- p.adjust(dt[,2], "BH")
#                  pair background_fvPTC_hyper adjust_p_value
#   Gene Body_Promoter                 0.1184         0.3552
# Gene Body_Intergenic                 1.0000         1.0000
#  Intergenic_Promoter                 1.0000         1.0000



####### -------------- CpGs-----------######### 
# Background 
library(data.table)
library(dplyr)
path = "/home/minhhoang/Documents/Thyroid cancer/TOTAL/Analysis/384_probes"
background <- as.data.frame(fread(paste0(path, "/Result/CpGs_background.tsv"), header = T))
background$main <- ifelse(background$annot.type=="CpG islands","CpG islands",ifelse(background$annot.type=="CpG shores","CpG shores", "Non CpG Island"))
background_1 <- as.data.frame(table(background$main))
#            Var1 Freq
#1    CpG islands 9147
#2     CpG shores 7246
#3 Non CpG Island 5202


# other subtypes 

folder <- c("NIFTP", "FC", "fvPTC", "FA")
status <- c("hyper", "hypo")


for (i in 1:4){
    print(folder[i])
    for (j in 1:2){
        print(status[j])
        data <- as.data.frame(fread(paste0(path,"/Result/",folder[i], "/CpGs_", status[j], "_", folder[i], "_normal_probes.tsv"))) 
        data$main <- ifelse(data$annot.type=="CpG islands","CpG islands",ifelse(data$annot.type=="CpG shores","CpG shores", "Non CpG Island"))
        data_1 <- as.data.frame(table(data$main))
        target <- data_1$Freq 
        print(data_1)
    }
}
 }

background <- c(9147, 7246, 5202)
NIFTP_hyper <- c(6,2,5)
NIFTP_hypo <- c(29,75,137)
FA_hyper <- c(0,1,5)
FC_hyper <- c(36,48,33)
FC_hypo <- c(0,2,1)
fvPTC_hyper <- c(0,1,2)
fvPTC_hypo <- c(0,0,2)

total <- cbind(background, NIFTP_hyper, NIFTP_hypo, FA_hyper, FC_hyper, FC_hypo, fvPTC_hyper, fvPTC_hypo)
rownames(total) <- background_1$Var1
mtotal <- as.data.frame(t(total))


res_table <- data.frame(
                        p_value1 = numeric(),
                        p_value2 = numeric(),
                        p_value3 = numeric(),
                        p_value4 = numeric(),
                        p_value5 = numeric(),
                        p_value6 = numeric(),
                        p_value7 = numeric(),
                        pair = character(),
                        stringsAsFactors=FALSE)

pair1 <- combn(unique(rownames(mtotal)),2)
pair2 <- combn(unique(colnames(mtotal)),2)



for (i in 1:7) {
    print(pair1[,i])
    sub <- pair1[,i]
    s1 <- sub[1]
    s2 <- sub[2] 
    colnames(res_table)[i] <- paste0(s1,"_",s2)
    tdte <- mtotal[pair1[,i],]
    for (j in 1:ncol(pair2)) {
        gr <- pair2[,j]
        g1 <- gr[1]
        g2 <- gr[2]
        res_table[j,8] <- paste0(g1,"_",g2)
        p.value <- fisher.test(tdte[,pair2[,j]])$p.value
    res_table[j,i] <- p.value    
    }
}
write.table(res_table,  paste0(path, "/Result/CpGs_test_result.tsv"), sep ="\t", quote=F, col.names=T, row.names=F)

res_table_f <- res_table[,c(8,1,2,3,4,5,6,7)] 

pair3 <- combn(unique(colnames(res_table_f)),2)
dt <- res_table_f[,pair3[,1]]
dt <- arrange(dt, background_NIFTP_hyper)
dt$adjust_p_value <- p.adjust(dt[,2], "BH")
#                       pair background_NIFTP_hyper adjust_p_value
#1  CpG shores_Non CpG Island              0.1373104      0.4119311
#2     CpG islands_CpG shores              0.4792761      0.5417056
#3 CpG islands_Non CpG Island              0.5417056      0.5417056


dt <- res_table_f[,pair3[,2]]
dt <- arrange(dt, background_NIFTP_hypo)
dt$adjust_p_value <- p.adjust(dt[,2], "BH")
#                         pair background_NIFTP_hypo adjust_p_value
#1 CpG islands_Non CpG Island          5.102863e-34   1.530859e-33
#2  CpG shores_Non CpG Island          3.576169e-11   5.364254e-11
#3     CpG islands_CpG shores          1.296622e-08   1.296622e-08


dt <- res_table_f[,pair3[,3]]
dt <- arrange(dt, background_FA_hyper)
dt$adjust_p_value <- p.adjust(dt[,2], "BH")
#                       pair background_FA_hyper adjust_p_value
#1 CpG islands_Non CpG Island         0.006273956     0.01882187
#2  CpG shores_Non CpG Island         0.088696869     0.13304530
#3     CpG islands_CpG shores         0.442051970     0.44205197


dt <- res_table_f[,pair3[,4]]
dt <- arrange(dt, background_FC_hyper)
dt$adjust_p_value <- p.adjust(dt[,2], "BH")
#                        pair background_FC_hyper adjust_p_value
#1     CpG islands_CpG shores          0.02035852     0.06107557
#2 CpG islands_Non CpG Island          0.05893126     0.08839689
#3  CpG shores_Non CpG Island          0.91023276     0.91023276


dt <- res_table_f[,pair3[,5]]
dt <- arrange(dt, background_FC_hypo)
dt$adjust_p_value <- p.adjust(dt[,2], "BH")
#                    pair background_FC_hypo adjust_p_value
#1     CpG islands_CpG shores          0.1954250      0.5438676
#2 CpG islands_Non CpG Island          0.3625784      0.5438676
#3  CpG shores_Non CpG Island          1.0000000      1.0000000

dt <- res_table_f[,pair3[,6]]
dt <- arrange(dt, background_fvPTC_hyper)
dt$adjust_p_value <- p.adjust(dt[,2], "BH")
#                       pair background_fvPTC_hyper adjust_p_value
#1 CpG islands_Non CpG Island              0.1314792      0.3944376
#2     CpG islands_CpG shores              0.4420520      0.5751779
#3  CpG shores_Non CpG Island              0.5751779      0.5751779


dt <- res_table_f[,pair3[,7]]
dt <- arrange(dt, background_fvPTC_hypo)
dt$adjust_p_value <- p.adjust(dt[,2], "BH")
#                       pair background_fvPTC_hypo adjust_p_value
#1 CpG islands_Non CpG Island             0.1314792      0.2620466
#2  CpG shores_Non CpG Island             0.1746977      0.2620466
#3     CpG islands_CpG shores             1.0000000      1.0000000
