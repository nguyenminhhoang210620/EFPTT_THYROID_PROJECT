##### Gene and CpGs analysis ############### 
library(data.table)
library(GenomicRanges)
library(wesanderson)
library(data.table)
library(annotatr)
library(GenomicRanges)
library(annotatr)

# create path 
path <- "/home/minhhoang/Documents/Thyroid cancer/TOTAL/Analysis/13_probes"

#craet annot 
annots = c('hg19_cpgs', 'hg19_basicgenes', 'hg19_genes_intergenic','hg19_genes_intronexonboundaries')
annotations = build_annotations(genome ='hg19', annotations = annots)

#loading files 
data_probes <- as.data.frame(fread(paste0(path,"/Data/DMPs_13_probes.tsv"), header = T))


# transfrom 
data_transform  <- makeGRangesFromDataFrame(data_probes , keep.extra.columns=TRUE)

# run 
dm_annotated = annotate_regions(regions = data_transform, annotations = annotations, ignore.strand = FALSE, quiet = FALSE)

# as.data.frame 
df_dm_annotated = data.frame(dm_annotated)

# save files 
write.table(df_dm_annotated, paste0(path,"/Result/gene_CpG_analysis.tsv"),sep = "\t", col.names = T, quote = F, row.names = F)
