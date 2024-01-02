# loadding packages 
library(wesanderson)
library(ggplot2)
library(data.table)
library(dplyr)
library(annotatr)
library(GenomicRanges)
library(annotatr)
library(Gviz)


# create path
path = "/home/minhhoang/Documents/Thyroid cancer/TOTAL/Analysis/13_probes/"


# DMP 
DMPs <- as.data.frame(fread(paste0(path, "/Data/DMPs_13_probes.tsv")))

## creat annot 
annots = c('hg19_cpgs', 'hg19_basicgenes', 'hg19_genes_intergenic','hg19_genes_intronexonboundaries')
annotations = build_annotations(genome ='hg19', annotations = annots)


# edit 
DMPs_1 <- DMPs[5,1:4]

#transform of data 
dm_regions <- makeGRangesFromDataFrame(DMPs_1, keep.extra.columns=TRUE)
      
# create annot 
dm_annotated = annotate_regions(regions = dm_regions,annotations = annotations, ignore.strand = FALSE, quiet = FALSE)

# make data.frame 
df_dm_annotated <- data.frame(dm_annotated)

# chromosome 
chr <- as.character(unique(seqnames(dm_annotated)))

# gene 
genome  <- "hg19"
gtrack <- GenomeAxisTrack()

# plot chr 
itrack <- IdeogramTrack(genome = gen, chromosome = chr)

#CpGs 
atrack <- (AnnotationTrack(dm_annotated[1], name = "cg11484872", fill = "red", background.title = "red"))
plotTracks(atrack)

# plot CpG region 
cpgisland <- dm_annotated@elementMetadata@listData[["annot"]][3]
cpgisland_1 <- as.data.frame(cpgisland)
atrack2 <- AnnotationTrack(cpgisland, name = "CpG Inter", background.title = "darkgreen", fill = "darkgreen")
plotTracks(atrack2)

#promoter 
gene_region <- dm_annotated@elementMetadata@listData[["annot"]][1]
gene_region_1 <- as.data.frame(gene_region)
atrack3 <- AnnotationTrack(cpgisland, name = "Gene Promoter", fill = "plum", background.title = "plum")
plotTracks(atrack3)


# calling cpG Region 
cpgisland <- dm_annotated@elementMetadata@listData[["annot"]][nrow(data.frame(dm_annotated))] 
cpgisland_1 <- as.data.frame(cpgisland)
df_dm_annotated = data.frame(dm_annotated)
txname <- df_dm_annotated$annot.tx_id[1:nrow(data.frame(dm_annotated))-1] 

# gene 
from <- 31543000
to <- 31546112
knownGenes <- UcscTrack(genome = "hg19", chromosome = "chr6", 
                        track = "knownGene", from = from, to = to,
                        trackType = "GeneRegionTrack", 
                        rstarts = "exonStarts", rends = "exonEnds", 
                        gene = "name", symbol = "name", 
                        transcript = "name", strand = "strand", 
                        fill = "#8282d2", name = "TNF",geneSymbol = T, background.title = "#8282d2" )


# CG _ CONTENT 
a <- as.data.frame(fread(paste0(path,'/Data/cpg_content .tsv')))
CG_content  <- makeGRangesFromDataFrame(a, keep.extra.columns=TRUE)
track <- DataTrack(CG_content, chromosome = "chr6", genome = "hg19", type = "hist", window = -1, windowSize = 150, 
                    fill.histogram = "black", col.histogram = "black",
                       ylim = c(20, 150), name = "GC Percent", background.title = "black" )


pdf(file = paste0(path, "/Result/Gviz_plot_2.pdf")
plotTracks(list(idxTrack, axTrack, knownGenes, atrack3, atrack2 ,track, atrack),from = 31543000, to = 31545500, sizes= c(2,3,5,4,4,4,3), col.border.title="#FFFEDB")
dev.off()



