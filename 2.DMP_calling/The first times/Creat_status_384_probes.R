###### creat file of 384 probes ######
# Create path 
path = "/home/minhhoang/Documents/Thyroid cancer/TOTAL/Calling_DMP/First_calling_DMP"
# loading packages 
library(data.table)
--------------------------------
### NIFTP vs normal ####
# loading 
NIFTP <- as.data.frame(fread(paste0(path, "/Result/DMP/NIFTP-normal/DMPs_filtered.hg19.tsv"),header=T))
dim(NIFTP)
# have 8085 probes 

summary(NIFTP$mean.Delta)

summary(NIFTP$FDR)
#   Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-0.38723 -0.08977 -0.04541 -0.06492 -0.03214  0.37382

# filter probes with mean.delta >= 0.2 and FDR < 0.01
NIFTP_f <- subset(NIFTP, abs(mean.Delta)>=0.2&FDR<0.01)
dim(NIFTP_f)
# have 284 probes 

# Determine status 
NIFTP_f$status <- ifelse(NIFTP_f$mean.Delta >0, "Hyper DMPs", "Hypo DMPs")
table(NIFTP_f$status)
#Hyper DMPs  Hypo DMPs 
#        18        266 
NIFTP_probe <- NIFTP_f$probeID

------------------------------------
### FC vs normal ### 
FC <- as.data.frame(fread(paste0(path, "/Result/DMP/normal-FC/DMPs_filtered.hg19.tsv"),header=T))

dim(FC)
# have 6902 probes 

summary(FC$mean.Delta)
#  Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-0.38817 -0.08256 -0.02097 -0.02747  0.02301  0.33238 

# filter probes with mean.delta >= 0.2 and FDR <0.01
FC_f <- subset(FC, abs(mean.Delta)>=0.2&FDR<0.01)

# creat status 
FC_f$status <- ifelse(FC_f$mean.Delta <0, "Hyper DMPs", "Hypo DMPs")
dim(FC_f)
# have  175 probes 

table(FC_f$status)
#Hyper DMPs  Hypo DMPs 
#      151        24
FC_probe <- FC_f$probeID

-------------------------------------
### FA vs normal ### 
FA <- as.data.frame(fread(paste0(path, "/Result/DMP/normal-FA/DMPs_filtered.hg19.tsv"),header=T))

dim(FA)
# have 5291 probes 

summary(FA$mean.Delta)
# 	Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-0.32446 -0.07897 -0.03511 -0.03035  0.02155  0.22276 

# filter probes with mean.delta >= 0.2 and FDR < 0.01
FA_f <- subset(FA, abs(mean.Delta)>=0.2&FDR<0.01)
FA_f$status <- ifelse(FA_f$mean.Delta <0, "Hyper DMPs", "Hypo DMPs")

dim(FA_f)
# have 36 probes 

table(FA_f$status)
#Hyper DMPs  Hypo DMPs 
#         32         4
FA_probe <- FA_f$probeID

-----------------------------------
# fvPTC vs normal 
fvPTC <- as.data.frame(fread(paste0(path, "/Result/DMP/normal-fvPTC/DMPs_filtered.hg19.tsv"),header=T))

dim(fvPTC)
# have 2158 probes 

summary(fvPTC$mean.Delta)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-0.24748 -0.09895 -0.06779 -0.04625  0.01709  0.26860 

#filter probes with mean.delta >= 0.2 and FDR < 0.01 
fvPTC_f <- subset(fvPTC, abs(mean.Delta)>=0.2&FDR<0.01)
fvPTC_f$status <- ifelse(fvPTC_f$mean.Delta < 0, "Hyper DMPs", "Hypo DMPs")

dim(fvPTC_f)
# have 32 probes 

table(fvPTC_f$status)
#Hyper DMPs  Hypo DMPs 
#        14         18 

fvPTC_probe <- fvPTC_f$probeID

######
#Check overlab probes 
NIFTP_probe <- NIFTP_f$probeID
FC_probe <- FC_f$probeID
FA_probe <- FA_f$probeID
fvPTC_probe <- fvPTC_f$probeID

# loading packages 
library(UpSetR)
library(ComplexHeatmap)
# fomat data by creating a list of probes 
probe <- list(NIFTP=NIFTP_probe,FA=FA_probe,FC=FC_probe,fvPTC=fvPTC_probe)
# creat the combination matrix
m = make_comb_mat(probe)
##A combination matrix with 4 sets and 12 combinations.
##  ranges of combination set size: c(2, 254).
##  mode for the combination size: distinct.
##  sets are on rows.

##Top 8 combination sets are:
##  NIFTP FA FC fvPTC code size
##      x             1000  254
##            x       0010  119
##         x  x       0110   21
##      x     x       1010    9
##            x     x 0011    9
##      x     x     x 1011    8
##      x           x 1001    6
##         x          0100    6

##Sets are:
##    set size
##  NIFTP  284
##     FA   36
##     FC  175
##  fvPTC   32
# draw upset plot 

pdf(file = paste0(path,'/Result/Upsetplot.pdf'), width = 10, height = 8)
UpSet(m, top_annotation = upset_top_annotation(m, add_numbers = TRUE),
    right_annotation = upset_right_annotation(m, add_numbers = TRUE))
dev.off()

# filter probe is specific for every subtypes => total 384 probes 
NIFTP_m <- extract_comb(m,"1000")      # 254 probes 
FA_m <- extract_comb(m,"0100")         # 6 probes 
FC_m <- extract_comb(m,"0010")         # 119 probes 
fvPTC_m <- extract_comb(m,"0001")      # 5 probes 


# data after filter 
NIFTP_f2 <- NIFTP_f[NIFTP_f$probeID%in%NIFTP_m,]
table(NIFTP_f2$status)
#Hyper DMPs  Hypo DMPs 
#        13        241 
write.table(NIFTP_f2, paste0(path,"/Result/NIFTP_normal_probes.tsv"),sep="\t", col.names=T, quote=F, row.names=F)

FC_f2 <- FC_f[FC_f$probeID%in%FC_m,]
table(FC_f2$status)
#Hyper DMPs  Hypo DMPs 
#        116        3
write.table(FC_f2, paste0(path,"/Result/FC_normal_probes.tsv"),sep="\t", col.names=T, quote=F, row.names=F)

FA_f2 <- FA_f[FA_f$probeID%in%FA_m,]
table(FA_f2$status)
#Hyper DMPs 
#       6 
write.table(FA_f2, paste0(path,"/Result/FA_normal_probes.tsv"),sep="\t", col.names=T, quote=F, row.names=F)

fvPTC_f2 <- fvPTC_f[fvPTC_f$probeID%in%fvPTC_m,]
table(fvPTC_f2$status)
#Hyper DMPs  Hypo DMPs 
#         3          2 
write.table(fvPTC_f2, paste0(path,"/Result/fvPTC_normal_probes.tsv"),sep="\t", col.names=T, quote=F, row.names=F)