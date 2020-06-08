#!/usr/bin/env Rscript

dir.create("int")

counts = readRDS("../data/counts-tumor.RDS")
metaInfo = readRDS("../data/cell_identity-tumor.RDS") 
weights = readRDS("../data/weights.RDS")

# load clustering result. 
k = 12
fit = readRDS("../data/bestFit.RDS")

library(plyr)
library(dplyr)
library(tidyr)
library(reshape2)
library(scales)
library(viridis)
library(ggpubr)
library(gridExtra)

source("poi_decom_gain/poi_decom.R", chdir = T)

# sample composition plots
# get the labels
ctl_labels <- fit$L
perm <- order(ctl_labels)

#' check the composition of each cluster
# barplots
cluster.colors <- setNames(c("turquoise", "blue", "brown", "yellow", "green", "red", "black", "pink", "magenta", "purple", "greenyellow", "tan"), levels(metaInfo$cluster) )
patient.colors <- setNames(c("yellow", "greenyellow", "magenta", "pink", "turquoise", "brown", "purple", "black", "red", "green", "blue"), levels(metaInfo$patient))

clusterCompo <- ddply(metaInfo, .(cluster, patient, treatmentPhase), nrow)
colnames(clusterCompo) <- c("cluster", "patient", "treatmentPhase",  "nCells")

pdf("sampleComposition-barplot.pdf", width = 26, height = 12)
ggplot(clusterCompo, aes(x = treatmentPhase, y = nCells, fill = cluster)) +
    geom_bar(position = "fill",stat = "identity") + theme(axis.line = element_line(colour = "black"),text = element_text(size=50),axis.text.x = element_text(angle=60, hjust=1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank()) +
    scale_y_continuous(labels = percent_format())+
    facet_grid(~ patient, scales = "free_x") +
    scale_fill_manual(values = cluster.colors) + 
    labs(y="percentage of cells", x = "treatment phase")
dev.off()

# make boxplot showing proportion changes
sampleCompo.prop <- clusterCompo %>% group_by(patient, treatmentPhase) %>% mutate(percent = nCells/sum(nCells)) %>% ungroup()
sampleCompo.prop <- as.data.frame(complete(sampleCompo.prop, cluster, patient, treatmentPhase, fill=list(nCells=0, percent=0))) 
sampleCompo.prop <- sampleCompo.prop[order(sampleCompo.prop$treatmentPhase), ]

pdf("pairComp-boxplot.pdf", height = 5, width = 6)
ggpaired(sampleCompo.prop, x = "treatmentPhase", y = "percent",
         palette = "jco",
	 point.size = 1.5,
	 ylab = "Proportion", 
         line.color = "gray33", line.size = 0.4,
         short.panel.labs = FALSE) +
	 facet_wrap(~ cluster, scale = "free_y") + 
	 theme(strip.background =element_rect(fill="white"))+
         stat_compare_means(label = "p.format", paired = TRUE)

dev.off()

# compute centroid for each gene in each cell
Z <- sapply(seq.int(ncol(counts)), function(i) poi_decom_control_centroid(fit$X %*% weights[, i, drop = F], counts[, i, drop = F], metaInfo$gain[i]))
colnames(Z) <- colnames(counts)
rownames(Z) <- rownames(counts)

Z.sqrt.scaled <- t(scale(t(sqrt(Z))))

# perform pair-wise DE 
source("poi_decom_gain/cluster.specific.genes.R", chdir = T)

clusterPairs = t(data.frame(combn(levels(metaInfo$cluster), 2), row.names = c("c1", "c2")))

clusterPair.degsStats = lapply(seq.int(nrow(clusterPairs)), function(x){
				c1 = clusterPairs[x, "c1"]
				c2 = clusterPairs[x, "c2"]
				degs.i = cluster.specific.genes(X = fit$X, D = weights, g = metaInfo$gain, Y = counts, which(metaInfo$cluster == c1), which(metaInfo$cluster == c2))
				degs.i$gene = rownames(degs.i)
				return(degs.i)
			}) 

saveRDS(clusterPair.degsStats, file = "int/clusterPair.degsStats.RDS")

clusterPair.degs.lik <- do.call("cbind", lapply(clusterPair.degsStats, function(x) x[, "lik"]))
rownames(clusterPair.degs.lik) <- rownames(counts)
saveRDS(clusterPair.degs.lik, file = "int/clusterPair.degs.lik.RDS")

clusterPair.degs.lik.rk <- apply(clusterPair.degs.lik, 2, function(x) rank(-x))
clusterPair.degs.p <- do.call("cbind", lapply(clusterPair.degsStats, function(x) x[, "p.adj"]))
clusterPair.degs.Z.p.e <- do.call("cbind", lapply(clusterPair.degsStats, function(x) x[, "Z.p.val.e"]))

clusterPair.degs.lik.flt <- clusterPair.degs.lik[rownames(clusterPair.degs.lik.rk)[rowSums((clusterPair.degs.lik.rk < 1000) * (clusterPair.degs.p < 0.01) * (clusterPair.degs.Z.p.e < 0.01)) > 0], ]

library(psych) 
clusterPair.lik.pearCoef <- cor(t(clusterPair.degs.lik.flt))
clusterPair.lik.pearPval <- corr.test(t(clusterPair.degs.lik.flt), adjust = "BH", ci = F) 

saveRDS(clusterPair.lik.pearCoef, file = "int/clusterPair.lik.pearCoef.flt.top1000.RDS")
saveRDS(clusterPair.lik.pearPval, file = "int/clusterPair.lik.pearPval.flt.top1000.RDS")

# use pearson correlation and take the pairs with abs(pearCoef) > 0.8 & pearPval < 0.01 
clusterPair.lik.nw <- clusterPair.lik.pearCoef
clusterPair.lik.nw[clusterPair.lik.pearCoef < 0.8 | clusterPair.lik.pearPval$p > 0.01] <- 0

clusterPair.lik.edges <- melt(clusterPair.lik.nw)
colnames(clusterPair.lik.edges) <- c("node1", "node2", "weight")
clusterPair.lik.edges <- subset(clusterPair.lik.edges, node1 != node2 & weight != 0)
write.table(clusterPair.lik.edges, "clusterPair_lik_gene_edges.csv", sep = "\t", row.names = F, col.names = F, quote = F)

# detect communities
library(igraph)
library(magrittr)
library(visNetwork)
library(data.table)

clusterPair.lik.graph <- graph_from_adjacency_matrix(clusterPair.lik.nw, weighted=TRUE, mode = "undirected", diag = F)
clusterPair.lik.graph <- simplify(clusterPair.lik.graph)

# walktrap.community
clusterPair.lik.communities.wt <- walktrap.community(clusterPair.lik.graph, weights = abs(E(clusterPair.lik.graph)$weight), steps = 3)

V(clusterPair.lik.graph)$community_wt <- clusterPair.lik.communities.wt$membership
saveRDS(clusterPair.lik.graph, file = "int/clusterPair.lik.graph.RDS")

# export the gene modules (based on walktrap with step = 3, only take the communities with >= 30 genes)
clusterPair.lik.communities.wt.s3.mbs <- table(clusterPair.lik.communities.wt$membership)
clusterPair.lik.communities.wt.s3.mbs <- clusterPair.lik.communities.wt.s3.mbs[order(clusterPair.lik.communities.wt.s3.mbs, decreasing = T)]
clusterPair.lik.communities.wt.s3.mbs.flt = clusterPair.lik.communities.wt.s3.mbs[clusterPair.lik.communities.wt.s3.mbs >= 30]

# significance test of each community based on node degrees
community.significance.test <- function(graph, vs, ...) {
    if (is.directed(graph)) stop("This method requires an undirected graph")
    subgraph <- induced.subgraph(graph, vs)
    in.degrees <- degree(subgraph)
    out.degrees <- degree(graph, vs) - in.degrees
    wilcox.test(in.degrees, out.degrees, ...)
}

community.p <- setNames(rep(NA, length(clusterPair.lik.communities.wt.s3.mbs.flt)), names(clusterPair.lik.communities.wt.s3.mbs.flt))
dir.create("comm2genes")
for (m in names(clusterPair.lik.communities.wt.s3.mbs.flt)){
    v.m = names(V(clusterPair.lik.graph))[which(V(clusterPair.lik.graph)$community_wt == m)]
    write.table( v.m, paste0("comm2genes/community_", m, ".csv"), row.names = F, col.names = F, quote = F)
    
    community.p[m] <- community.significance.test(clusterPair.lik.graph, v.m)$p.value 
}

clusterPair.lik.communities.wt.s3.mbs.flt = clusterPair.lik.communities.wt.s3.mbs.flt[community.p < 0.01]

# embedding of the communities (done in gene_nw.ipynb)
clusterPair.lik.communities.wt.s3.gc2genes.flt <- setNames(lapply(names(clusterPair.lik.communities.wt.s3.mbs.flt), function(m) names(V(clusterPair.lik.graph)[which(V(clusterPair.lik.graph)$community_wt == m)])), paste0("gc", names(clusterPair.lik.communities.wt.s3.mbs.flt)))

write.table(t(data.frame(lapply(clusterPair.lik.communities.wt.s3.gc2genes.flt, "length<-", max(lengths(clusterPair.lik.communities.wt.s3.gc2genes.flt))))), "comm2genes.csv", sep = "\t", row.names = T, col.names = F, quote = F) 
 
# K-core decomposition (for each community)
dir.create("comm2genes_maxCore")

clusterPair.lik.genes.maxCore <- list()
clusterPair.lik.hubness <- list()
clusterPair.lik.betweenness <- list()
clusterPair.lik.closeness <- list()
clusterPair.lik.centrality <- list()

for (m in names(clusterPair.lik.communities.wt.s3.mbs.flt)){
    clusterPair.lik.subgraph.m <- induced.subgraph(clusterPair.lik.graph, V(clusterPair.lik.graph)[which(V(clusterPair.lik.graph)$community_wt == m)] )
    m.coreness <- coreness(clusterPair.lik.subgraph.m)
    m.coreness <- m.coreness[order(m.coreness, decreasing = T)]
    m.genes.maxCore <- names(m.coreness)[which(m.coreness == max(m.coreness))]
    
    if (length(m.genes.maxCore) < 30){
	m.genes.maxCore <- names(m.coreness)[which(m.coreness >= m.coreness[30])]
    }

    clusterPair.lik.genes.maxCore[[m]] <- m.genes.maxCore
    write.table( m.genes.maxCore, paste0("comm2genes_maxCore/community_", m, ".csv"), row.names = F, col.names = F, quote = F)    

    # check hubness as well 
    m.hubness <- hub.score(clusterPair.lik.subgraph.m, weights = E(clusterPair.lik.subgraph.m)$weight)
    clusterPair.lik.hubness[[m]] <- m.hubness$vector[order(m.hubness$vector, decreasing = T)]
    
    m.centrality <- data.frame(row.names = names(m.coreness), gene = names(m.coreness), coreness = m.coreness, hubness = m.hubness$vector, stringsAsFactors = F)
    clusterPair.lik.centrality[[m]] <- m.centrality

}

names(clusterPair.lik.genes.maxCore) <- paste0("gcCore", names(clusterPair.lik.genes.maxCore))
saveRDS(clusterPair.lik.genes.maxCore, file = "int/clusterPair.lik.genes.maxCore.RDS")
saveRDS(clusterPair.lik.centrality, file = "int/clusterPair.lik.centrality.RDS")

# further decompose gcCores based on there overlappings with the significantly enriched gene-sets in cpdb 
oraRes.dir <- "/mnt/storageBig7/home/kaizhang/OvCa-SC/10XGC/180425_D00482_0211_ACCC3DANXX_cellRanger3.1/codes/data/ORA_cpdb/ORA_"

ora.ovl.mat <- list()
gcCore.gs.ovlp <- list() 
for (gc in names(clusterPair.lik.genes.maxCore)){
    ora.gc <- read.table(paste0(oraRes.dir, gc, ".csv"), sep = "\t", header = T, stringsAsFactors = F, comment.char = "")
    ora.gc <- subset(ora.gc, q.value < 0.05 & effective_size < 500)
    ora.gc$pathway <- gsub(" |, |-", "_", ora.gc$pathway) 
    ovlp.gc <- setNames(strsplit(ora.gc$members_input_overlap, "; "), ora.gc$pathway)
    gcCore.gs.ovlp[[gc]] <- ovlp.gc

    ovlp.gc.melt = melt(ovlp.gc, value.name = "gene")
    ovlp.gc.melt = ovlp.gc.melt[!duplicated(ovlp.gc.melt), ]
    colnames(ovlp.gc.melt) = c("gene", "gset")    

    ora.ovl.mat[[gc]] <- table(ovlp.gc.melt) 
}
saveRDS(gcCore.gs.ovlp, file = "int/gcCore.gs.ovlp.RDS")
saveRDS(ora.ovl.mat, file = "int/ora.ovl.mat.RDS") 

# get filtered clusterPair.lik.genes.maxCore, based on gene-sets in cpdb (overlapped genes with significantly enriched gene-sets(q.value < 0.05 & effective_size < 500))
clusterPair.lik.genes.maxCore.flt <- lapply(ora.ovl.mat, rownames) 
saveRDS(clusterPair.lik.genes.maxCore.flt, file = "int/clusterPair.lik.genes.maxCore.flt.RDS")

# block clustering for gcCores with number of genes > 20
ora.ovl.mat.large <- ora.ovl.mat[which(sapply(ora.ovl.mat, nrow) > 20 )]

library(blockcluster)

dir.create("gcCore_gsets_blockcluster")
gcCore.gs.bc <- list()
gcCore.gs.bc.maxICL <- list()
for (gc in names(ora.ovl.mat.large)){
    ora.ovl.mat.gc = ora.ovl.mat.large[[gc]]
    gc.gs.bc <- list()
    for (i in seq.int(min(nrow(ora.ovl.mat.gc), 10))){
        for (j in seq.int(min(ncol(ora.ovl.mat.gc), 10))){
            		
            set.seed(19901012)
            gc.gs.bc.i_j <- cocluster( as.matrix.data.frame(ora.ovl.mat.gc), "binary", semisupervised = FALSE, model = "pik_rhol_epsilonkl", nbcocluster = c(i, j) , strategy = coclusterStrategy(), nbCore = 6)
            gc.gs.bc[[paste0(i, "_", j)]] <- gc.gs.bc.i_j 

        }
    } 

    gcCore.gs.bc[[gc]] <- gc.gs.bc
    
    gc.gs.bc <- gc.gs.bc[c(1:5, 11:15, 21:25, 31:35, 41:45)]

    # get the clusterign with the max likelihood
    gc.gs.bc.maxICL <- gc.gs.bc[[names(which.max(lapply(gc.gs.bc, function(x) x@ICLvalue)))]]
    gcCore.gs.bc.maxICL[[gc]] <- gc.gs.bc.maxICL    

    pdf(paste0("gcCore_gsets_blockcluster/", gc, ".pdf"))
    plot(gc.gs.bc.maxICL,  asp = 0)
    dev.off()         
    
}
saveRDS(gcCore.gs.bc.maxICL, file = "int/gcCore.gs.bc.maxICL.RDS")
saveRDS(gcCore.gs.bc, file = "int/gcCore.gs.bc.RDS")

# further filter ora.ovl.mat.large based on annotations
clusterPair.lik.genes.maxCore.c <- clusterPair.lik.genes.maxCore.flt 
clusterPair.lik.genes.maxCore.c[["gcCore74"]] <- rownames(ora.ovl.mat.large[["gcCore74"]])[which(gcCore.gs.bc.maxICL[["gcCore74"]]@rowclass != 1)]

clusterPair.lik.genes.maxCore.c[["gcCore15"]] <- rownames(ora.ovl.mat.large[["gcCore15"]])[which(gcCore.gs.bc.maxICL[["gcCore15"]]@rowclass != 0)]

clusterPair.lik.genes.maxCore.c[["gcCore93"]] <- rownames(ora.ovl.mat.large[["gcCore93"]])[which(gcCore.gs.bc.maxICL[["gcCore93"]]@rowclass %in% c(0, 3, 2, 4))]
 
# renames gcCores based on their sizes 
gcNames = setNames(paste0("gc", seq.int(length(clusterPair.lik.genes.maxCore.c))), names(clusterPair.lik.genes.maxCore.c)[order(sapply(clusterPair.lik.genes.maxCore.c, length), decreasing = T)])

gene2gcCore.c = data.frame(gcCore = rep(names(clusterPair.lik.genes.maxCore.c), sapply(clusterPair.lik.genes.maxCore.c, length)), row.names = do.call("c", clusterPair.lik.genes.maxCore.c))

gene2gcCore.c$geneCommunity = factor(gcNames[match(gene2gcCore.c$gcCore, names(gcNames))], levels = gcNames[match(c("gcCore15", "gcCore9", "gcCore93", "gcCore74", "gcCore40", "gcCore21", "gcCore114", "gcCore49", "gcCore36", "gcCore119"), names(gcNames))] )

names(clusterPair.lik.genes.maxCore.c ) = gcNames[match(names(clusterPair.lik.genes.maxCore.c), names(gcNames))]

saveRDS(clusterPair.lik.genes.maxCore.c, file = "int/clusterPair.lik.genes.maxCore.c.RDS")

write.table(t(data.frame(lapply(clusterPair.lik.genes.maxCore.c, "length<-", max(lengths(clusterPair.lik.genes.maxCore.c))))), "comm2genes_maxCore-c.csv", sep = "\t", row.names = T, col.names = F, quote = F)


# make heatmap based on clusterPair.lik.genes.maxCore.c
cores.hm.c <- Z.sqrt.scaled[do.call("c", clusterPair.lik.genes.maxCore.c), order(ctl_labels)]

# this is only for visulization 
cores.hm.c[cores.hm.c > 2] = 2
cores.hm.c[cores.hm.c < -2] = -2

tars = c("HLA-DMA", "HLA-DMB", "HLA-DRA", "FOS", "JUN", "CDKN1A", "NFKBIA", "GADD45G", "ITGA2", "SMAD3", "COL1A2", "MMP1", "MUC15", "B3GNT7", "GCNT3", "MAD2L1", "BRCA2", "MCM3", "CHEK1", "RAD51", "BRCA1", "FANCI", "HMGB2", "CCL20", "CXCL1", "IL1R2", "IL1RN", "ROCK1", "TJP1", "ACIN1", "PSMB6", "PSMC5", "PSMA4", "NDUFB5", "UQCRC1", "IFIT1", "ISG15", "STAT2")

library(ComplexHeatmap)
jpeg("gcCores_flt-hm.jpeg",  width = 3200, height = 3000, quality = 100, res = 300)
Heatmap(cores.hm.c,
    show_row_names = F, 
    show_column_names = F, 
    top_annotation = HeatmapAnnotation(Patient = as.character(metaInfo$patient[match(colnames(cores.hm.c), rownames(metaInfo))]),
        col = list(Patient = patient.colors[levels(metaInfo$patient)]),
        simple_anno_size = unit(2.5, "mm")
        ),
    right_annotation = rowAnnotation(foo = anno_mark(at = match(tars, rownames(cores.hm.c)), labels = tars, labels_gp = gpar(fontsize = 11))),
    cluster_rows = F,
    cluster_columns = F,
    column_split = metaInfo$cluster[match(colnames(cores.hm.c), rownames(metaInfo))],
    row_split = gene2gcCore.c$geneCommunity[match(rownames(cores.hm.c), rownames(gene2gcCore.c))],
    heatmap_legend_param = list( title = "Expression", legend_height = unit(4, "cm"), labels_gp = gpar(fontsize = 11) ),
    column_gap = unit(1.4, "mm"),
    row_gap = unit(1.2, "mm")
)
dev.off()

#' score the cores of each community for each cell using ssgsea
# get columnwise relative ranks
unit.ranks <- function(X) {
    for (j in seq_len(ncol(X))) {
        m <- !is.na(X[, j])
        X[m, j] <- ( rank( X[m, j] ) - 1. ) / ( sum(m) - 1. )
    }
    return (X)
}

Z.r <- t(unit.ranks( t(Z) ))
saveRDS(Z.r, file = "int/Z.rank.RDS")

source("ssgsea/ssgsea.R", chdir = T)

ssgsea.gcCore.c.cw = lapply(seq.int(length(clusterPair.lik.genes.maxCore.c)), function(x) ssgsea(Z.r, clusterPair.lik.genes.maxCore.c[[x]], alpha = .75, compute.p = 1000L) )

saveRDS(ssgsea.gcCore.c.cw, file = "int/ssgsea.gcCore.c.cw.RDS")

ssgsea.gcCore.c.cw.es <- do.call("cbind", lapply(ssgsea.gcCore.c.cw, function(x) x[, "es"]))
rownames(ssgsea.gcCore.c.cw.es) <- colnames(Z)
colnames(ssgsea.gcCore.c.cw.es) <- names(clusterPair.lik.genes.maxCore.c)

ssgsea.gcCore.c.cw.es.melt = melt(ssgsea.gcCore.c.cw.es)
colnames(ssgsea.gcCore.c.cw.es.melt) <- c("cell", "gcCore", "ES") 
ssgsea.gcCore.c.cw.es.melt$cluster = metaInfo$poi_decom_cluster_k12[match(ssgsea.gcCore.c.cw.es.melt$cell, rownames(metaInfo))]

# violin plots of interval-primary pairs
treatmentPhase.colors = setNames(c("#00BFC4", "#F8766D"), c("primary", "interval"))

ssgsea.gcCore.c.cw.es.melt$treatmentPhase <- metaInfo$treatmentPhase[match(ssgsea.gcCore.c.cw.es.melt$cell, rownames(metaInfo))]
ssgsea.gcCore.c.cw.es.melt$patient <- metaInfo$patient[match(ssgsea.gcCore.c.cw.es.melt$cell, rownames(metaInfo))]

pdf("treatmentPhase_gcCore_c_violinPlots.pdf", width = 15, height = 15)
ggplot(ssgsea.gcCore.c.cw.es.melt, aes(y = ES, x = patient, fill = treatmentPhase)) +  geom_violin(trim = T, position = position_dodge()) + geom_point(shape=16, position = position_jitterdodge(seed = 1), size=0.1) + facet_wrap(~gcCore, nrow = 5, scales = "free") + theme(axis.line = element_line(colour = "black"), panel.grid.minor = element_blank(),panel.background = element_blank()) + scale_fill_manual(values = treatmentPhase.colors) 
dev.off()

pdf("treatmentPhase_gcCore6_violinPlots.pdf", width = 9.5, height = 3.5)
ggplot(subset(ssgsea.gcCore.c.cw.es.melt, gcCore == "gc6"), aes(y = ES, x = treatmentPhase, fill = treatmentPhase)) +  geom_violin(trim = T, position = position_dodge()) + geom_point(shape=16, position = position_jitterdodge(seed = 1), size=0.3) + facet_wrap(~patient, nrow = 1) + theme(axis.line = element_line(colour = "black"), axis.text.x = element_text(angle = 60, hjust = 1), panel.grid.minor = element_blank(),panel.background = element_blank(), text = element_text(size= 15)) + scale_fill_manual(values = treatmentPhase.colors) + labs(y="Stress score", x = "treatment phase") + stat_compare_means(label.y = 6500, label = "p.signif") 
dev.off()

# plot each gene community score in each cell cluster
ssgsea.gcCore.c.cw.es.melt$cluster = metaInfo$cluster[match(ssgsea.gcCore.c.cw.es.melt$cell, rownames(metaInfo))]

pdf("cellCluster_gcCore_flt_c_violinPlots.pdf", width = 9, height = 10)
ggplot(ssgsea.gcCore.c.cw.es.melt, aes(y = ES, x = cluster, fill = cluster)) +  geom_violin(trim = T) + geom_jitter(shape=16, position=position_jitter(0.1), size=0.1) + facet_wrap(~gcCore, nrow = 4, scales = "free") + theme(axis.line = element_line(colour = "black"), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.text.x = element_text(angle = 45), strip.text.x = element_text(size = 15))  + scale_fill_manual(values = cluster.colors) + labs(x="cluster", y="Enrichment score")
dev.off()


