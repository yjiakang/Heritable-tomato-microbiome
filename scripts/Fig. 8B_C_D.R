# Fig. 8B
library("dplyr")
library("ggplot2")
my_theme()
all_genes_top500 <- read.table("all_genes.txt", header = TRUE, row.names = 1) %>%
    mutate(sumfpkm = rowSums(.)) %>% 
    arrange(.,desc(sumfpkm)) %>%
    .[1:500,] %>% 
    select(-c("sumfpkm")) %>% 
    t %>%
    data.frame
    
pca_data <- prcomp(all_genes_top500, scale. = TRUE)
pca_df <- data.frame(pca_data$x) %>% 
    mutate(Treat = rep(c("CK", "Cra20", "KT2440"), each = 3))

percentage <- round(pca_data$sdev / sum(pca_data$sdev) * 100, 2)
percentage <- paste( colnames(pca_df), "(", paste( as.character(percentage), "%", ")", sep=""))

ggplot(pca_df, aes(PC1, PC2, color = Treat)) + 
    geom_point(size = 3) + 
    xlab(percentage[1]) + 
    ylab(percentage[2]) +
    scale_color_brewer(palette = "Set2")
ggsave("all_pca.pdf", height = 4, width = 5) # This is Fig. 8B


# Fig. 8C
library("dplyr")
library("magrittr")
library("ComplexHeatmap")
library("tibble")
library("RColorBrewer")
library("ggplot2")
my_theme()
all_genes <- read.table("all_genes.txt", header = TRUE) # fpkm

# DEGs resulted from DESeq2: |log2foldchange| > 1, padj < 0.05
Cra20 <- read.table("read_count.CK_vs_Cra20.DESeq2.DE_results.xls", sep = "\t", header = T) %>% 
    filter(abs(log2FoldChange.Cra20.CK.) > 1 & padj < 0.05) %>%
    mutate(Cra20_vs_CK = ifelse(log2FoldChange.Cra20.CK. > 0, "Up", "Down")) %>% 
    select(c("id","Cra20_vs_CK"))
    
KT2440 <- read.table("read_count.CK_vs_KT2440.DESeq2.DE_results.xls", sep = "\t", header = T) %>% 
    filter(abs(log2FoldChange.KT2440.CK.) > 1 & padj < 0.05) %>%
    mutate(KT2440_vs_CK = ifelse(log2FoldChange.KT2440.CK. > 0, "Up", "Down")) %>% 
    select(c("id","KT2440_vs_CK"))


all <- Cra20 %>% 
    full_join(KT2440) %>% 
    left_join(all_genes)

all[is.na(all)] <- "No"
## mean func
meanGp <- function(v) {
    require('magrittr')
    res <- v %>%
        split(rep(1 : 3, each = 3)) %>%
        sapply(mean, na.rm = TRUE)
    return(res)
}

all_for_cluster <- select(all, -contains("vs")) # select fpkm
rownames(all_for_cluster) <- all_for_cluster$id
all_for_cluster <- all_for_cluster[-1]

## sample name
sampleN <- c("CK","Cra20", "KT2440")

meanCount <- all_for_cluster %>%
    apply(1, meanGp) %>%
    t

colnames(meanCount) <- sampleN
## scale
scaleCount <- meanCount %>%
    t %>%
    scale %>%
    t
scaleCount %<>% .[complete.cases(.), ]

set.seed(123) # set seed to get consistent results
## choose cluster num
## 1. sum of squared error
wss <- (nrow(scaleCount) - 1) * sum(apply(scaleCount, 2, var))

for (i in 2:20) {
    wss[i] <- sum(kmeans(scaleCount,
                         centers=i,
                         algorithm = 'MacQueen')$withinss)
}

ggplot(tibble(k = 1:20, wss = wss), aes(k, wss)) +
    geom_point(colour = '#D55E00', size = 3) +
    geom_line(linetype = 'dashed') +
    xlab('Number of clusters') +
    ylab('Sum of squared error')

ggsave("Sum_of_squared_error.pdf", height = 3, width = 4)
## 2. Akaike information criterion
kmeansAIC = function(fit){
    m = ncol(fit$centers)
    n = length(fit$cluster)
    k = nrow(fit$centers)
    D = fit$tot.withinss
    return(D + 2*m*k)
}

aic <- numeric(20)
for (i in 1:20) {
    fit <- kmeans(x = scaleCount, centers = i, algorithm = 'MacQueen')
    aic[i] <- kmeansAIC(fit)
}

ggplot(tibble(k = 1:20, aic = aic), aes(k, wss)) +
    geom_point(colour = '#009E73', size = 3) +
    geom_line(linetype = 'dashed') +
    xlab('Number of clusters') +
    ylab('Akaike information criterion')

ggsave("Akaike_information_criterion.pdf", height = 3, width = 4)
# choose cluster 6

kclust6 <- kmeans(scaleCount, centers = 6, algorithm = "MacQueen", nstart = 1000, iter.max = 20)
cl <- as.data.frame(kclust6$cluster) %>% 
    rownames_to_column("id") %>% 
    set_colnames(c("id","cl"))

heat_all <- all %>% left_join(cl)


# export clustering result for further enrichment analysis
degs_cl <- heat_all %>%
    select(c("id","cl"))
# write.table(degs_cl, "./enrichment/degs_cl.txt", sep = "\t", quote = FALSE)


scaleC <- heat_all %>% 
    select(-contains("vs")) %>% 
    select(-c("id","cl")) %>% 
    t %>% 
    scale %>%
    t %>%
    as_tibble %>%
    bind_cols(heat_all %>% select(id, cl))

## col annotation
Treatment <- HeatmapAnnotation(Treatment = c(rep(c("CK", "Cra20", "KT2440"), each = 3)),
                             col = list(Treatment = c("CK" = "white", 
                                                        "Cra20" = "grey", "KT2440" = "grey50")),
                             gp = gpar(col = "black"))


## DEG annotation
deg <- heat_all %>% select(matches("vs"))

Heatmap(matrix = scaleC[1:9],
        name = 'Scaled Counts',
        row_split = scaleC$cl,
        row_gap = unit(2, "mm"),
        column_order = 1:9,
        # column_split = factor(rep(c("CK", "Cra20", "KT2440"), each = 3),
                              # levels = c("CK", "Cra20", "KT2440")),
        show_column_names = FALSE,
        col = colorRampPalette(rev(brewer.pal(n = 10, name = 'Spectral'))[c(-3,-8)])(100),
        top_annotation = Treatment,
        row_title_gp = gpar(fontsize = 10),
        use_raster = FALSE) +
        Heatmap(deg,
                col = c('Down' = '#00bbf9', 'No' = 'white', 'Up' = '#f94144'),
                column_names_gp = gpar(fontsize = 8),
                heatmap_legend_param = list(title = 'DEGs'),
                cluster_columns = FALSE,
                column_names_rot = 45,
                use_raster = FALSE)

# Fig. 8D
library("clusterProfiler")
library("magrittr")
library("tidyverse")
library("RColorBrewer")
library("AnnotationHub")
my_theme()

# > packageVersion("AnnotationHub")
# [1] ‘3,4,0’
# > packageVersion("clusterProfiler")
# [1] ‘4.4.1’

hub <- AnnotationHub()
# snapshotDate(): 2022-04-21
query(hub, "solanum")
# AnnotationHub with 8 records
# # snapshotDate(): 2022-04-21
# # $dataprovider: ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/, Inparanoid8, WikiPathways
# # $species: Solanum lycopersicum, Solanum tuberosum, Solanum pennellii, Solanum lycopersicum_...
# # $rdataclass: OrgDb, Inparanoid8Db, Tibble
# # additional mcols(): taxonomyid, genome, description, coordinate_1_based, maintainer,
# #   rdatadateadded, preparerclass, tags, rdatapath, sourceurl, sourcetype 
# # retrieve records with, e.g., 'object[["AH10593"]]' 
# 
# title                                             
# AH10593  | hom.Solanum_lycopersicum.inp8.sqlite              
# AH10606  | hom.Solanum_tuberosum.inp8.sqlite                 
# AH91815  | wikipathways_Solanum_lycopersicum_metabolites.rda 
# AH101115 | org.Solanum_pennellii.eg.sqlite                   
# AH101174 | org.Solanum_tuberosum.eg.sqlite                   
# AH101266 | org.Solanum_esculentum.eg.sqlite                  
# AH101267 | org.Solanum_lycopersicum.eg.sqlite                
# AH101268 | org.Solanum_lycopersicum_var._humboldtii.eg.sqlite

sly.db <- hub[["AH101267"]]
kmeansRes <- read.table("degs_cl.txt")

prefix <- 'kmeans6'
savepath <- "D:/phd/tomato_micr/RNASeq/Cra20_KT/cluster/enrichment/" # please change to your own work path

for (i in kmeansRes$cl %>% unique) {   
    ## KEGG
    kk2 <- enrichKEGG(gene = kmeansRes %>% filter(cl == i) %>% .$id %>%
                      bitr(.,"SYMBOL", "ENTREZID", sly.db) %>% dplyr::select("ENTREZID") %>%
                          unlist(),
                      organism = 'sly',
                      pvalueCutoff = 0.05)
    
    write.table(as.data.frame(kk2),
              paste0(prefix, '_cl', i, '_cp_KEGG.txt') %>% file.path(savepath, .),
              quote = FALSE,
              sep = "\t")
}

kallKEGG_input <- lapply(kmeansRes$cl %>% unique, function(x) {
    
    eachG <- kmeansRes %>% filter(cl == x) %>% 
        .$id %>% 
        bitr(.,"SYMBOL", "ENTREZID", sly.db) %>% 
        dplyr::select("ENTREZID") %>% 
        unlist()
    
    return(eachG)
    
}) %>%
    set_names(kmeansRes$cl %>% unique %>% paste0('cl', .))

# KEGG had been transfered to HTTPS from 2022.6.1, please refer to the following solution if needed.
options(clusterProfiler.download.method = "wininet")

kallKEGG <- compareCluster(geneCluster = kallKEGG_input,
                           fun = 'enrichKEGG',
                           organism = "sly",
                           pvalueCutoff = 0.05)

dotplot(kallKEGG)
ggsave('kmeans6_KEGGALL.pdf', width = 6, height = 4)

kallKEGG %>% 
    as.data.frame %>%
    write.table('kmeans6_KEGG.txt', quote = FALSE, sep = "\t")
# Adobe illustrator was used to adjust the figure legend "p.adjust" to "Padj"