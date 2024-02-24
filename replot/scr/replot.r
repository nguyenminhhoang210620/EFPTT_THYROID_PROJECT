# THIS SCRIPT TO VISUALIZATION PATHWAYS RELATED ###

# loading library
library(ggplot2)
library(dplyr)
#### ~~~~~~~~~~~~~~~~~~~~~ 13 and 38 PROBES go analysis ~~~~~~~~~~~~~~~~ ####
go_bp = read.delim("replot/data/13_probes/gen_ontholygy/ontology_GO_BP.tsv")
go_cc = read.delim("replot/data/13_probes/gen_ontholygy/ontology_GO_CC.tsv")
go_mf = read.delim("replot/data/13_probes/gen_ontholygy/ontology_GO_MF.tsv")


## sort 
sort_function = function(data){
    data = data[order(data$P.DE, decreasing = F),]
    data_collect = data[1:10,c("ONTOLOGY","TERM","P.DE")]
    return(data_collect)
}

go_bp_collect = sort_function(go_bp)
go_cc_collect = sort_function(go_cc)
go_mf_collect = sort_function(go_mf)

## data plot
dataPlot = do.call("rbind", list(go_bp_collect,go_cc_collect,go_mf_collect))
dataPlot = dataPlot[complete.cases(dataPlot),]
dataPlot$TERM = factor(dataPlot$TERM, levels=unique(dataPlot$TERM[order(dataPlot$ONTOLOGY,dataPlot$P.DE,dataPlot$TERM, decreasing = T)]), ordered=TRUE)


## ploting
ggplot(dataPlot, aes(x=-log10(P.DE), y= TERM, fill=ONTOLOGY)) + 
        geom_bar(stat="identity", width = 0.2)+
        theme_minimal()+
        labs(y="", fill="", x="-log10(p-value)",title="GO functional analysis (13 probes)")+
        scale_fill_manual(values=c("#FFA732","#C21292","#EF4040"))+
        theme(
            legend.position = c(.95, .95),
            legend.justification = c("right", "top"),
            legend.box.just = "right",
            legend.margin = margin(6, 6, 6, 6)
        )+
        theme(text=element_text(size=23))
ggsave("replot/result/GO_13_probes.pdf", width=20)
        


#### ~~~~~~~~~~~~~~ KEGG FUNCTION ANALYSIS ~~~~~~~~~~~~~ ####


# loading data Kegg
Kegg <- read.delim("replot/data/384_probes/gene_onthology/ontology_KEGG.tsv.tsv", sep = "\t")
Kegg <- Kegg[Kegg$P.DE <= 0.05, ]
Kegg <- Kegg[order(Kegg$P.DE), ]
Kegg$log10 <- -log10(Kegg$P.DE)
Kegg$rich_ratio <- Kegg$DE / Kegg$N
Kegg = Kegg[order(Kegg$P.DE),]
Kegg_collect = Kegg[1:15,]

# order the rich_ratio
Kegg_collect$Description <- with(Kegg_collect, reorder(Description, rich_ratio, na.rm = T))

ggplot(Kegg_collect, aes(x = rich_ratio, y = Description)) +
    geom_segment(aes(yend = Description), xend = 0, colour = "grey50") +
    geom_point(aes(colour = P.DE, size = N)) +
    scale_color_gradientn(colours = rainbow(5)) +
    scale_size(range = c(4,15))+
    theme_bw() +
    labs(x = "Risk Factor", y = NULL, color = "P-value", size = "NO. Gene", title="KEGG functional analysis (384 probes)") +
    theme(
        text = element_text(size = 25),
        axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold", size = 15)
    )
ggsave("replot/result/Kegg_384_probes.pdf", height = 12)
