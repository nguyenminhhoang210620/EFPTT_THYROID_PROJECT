# Loading packages 
library(data.table)
library(ggplot2)

# creat path 
path = "/home/minhhoang/Documents/Thyroid cancer/TOTAL/Analysis/13_probes"

# loadding data 
Kegg_13 <- as.data.frame(fread(paste0(path,"/Result/gen_ontholygy/ontology_KEGG.tsv.tsv")))


# filter pathways related to thyroid tumors 
Kegg_13_fill <- Kegg_13[c(4,5,6,12,13,14,16),] 
Kegg_13_fill$log10 <- -log10(Kegg_13_fill$P.DE)
Kegg_13_fill$rich_ratio <- Kegg_13_fill$DE/Kegg_13_fill$N

ggplot(Kegg_13_fill, aes(x =rich_ratio, y = Description)) +
		geom_segment(aes(yend = Description), xend = 0, colour = "grey50") +
  		geom_point(aes(colour = P.DE, size =N)) +
  		scale_color_gradientn(colours = rainbow(5))+
  		theme_bw() +
  		labs(x = 'Risk Factor', y = NULL, color = 'P-value',size = 'NO. Gene')+
  		theme(text=element_text(size = 20),
  			  axis.title = element_text(face = 'bold'),
  			  axis.text = element_text(face = 'bold', size = 15))
 ggsave(paste0(path,"/Result/gen_ontholygy/ontology_KEGG_13_top.pdf"))												


