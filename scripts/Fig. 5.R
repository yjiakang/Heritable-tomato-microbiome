library(DESeq2)
library(dplyr)
library(ggplot2)
library(tibble)
my_theme()
# package.version("DESeq2")
# [1] "1.28.1"

pathway <- read.delim("pab_keggpathway_unstratified.tsv",
                      row.names = 1, header = TRUE, sep = "\t") %>%
    round
design <- read.table("design.txt", header = TRUE, row.names = 1, sep = "\t")


# T1_HG64
design_T1_HG64 <- design[grep("T1", design$Treatment),] %>%
    .[-which(.$Treatment %in% c("HF12T1")),]
design_T1_HG64$Property <- factor(design_T1_HG64$Property, levels = c("Susceptible", "Resistant"))
pathway_T1_HG64 <- pathway[, rownames(design_T1_HG64)]
pathway_T1_HG64_dds <- DESeqDataSetFromMatrix(pathway_T1_HG64, design_T1_HG64, design = ~ Property)
keep <- rowSums(counts(pathway_T1_HG64_dds) >= 10) >= 3
pathway_T1_HG64_dds <- pathway_T1_HG64_dds[keep, ]
vs_pathway_T1_HG64_dds <- rlog(pathway_T1_HG64_dds, blind = FALSE) 

# DEGs analysis
T1_HG64_dds <- DESeq(pathway_T1_HG64_dds)
resultsNames(T1_HG64_dds)
# [1] "Intercept"                         "Property_Resistant_vs_Susceptible"
T1_HG64_dds_res <- results(T1_HG64_dds)
summary(T1_HG64_dds_res)
# out of 149 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 15, 10%
# LFC < 0 (down)     : 22, 15%
# outliers [1]       : 1, 0.67%
# low counts [2]     : 0, 0%
# (mean count < 6)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results
# plot
T1_HG64_dds_res_sig <- subset(T1_HG64_dds_res, (log2FoldChange > 1 | log2FoldChange < -1)) %>% 
    as.data.frame() %>%
    arrange(desc(log2FoldChange)) %>%
    replace(is.na(.),1) %>%
    rownames_to_column() %>% 
    rename(Category = rowname) %>%
    mutate(Signif = ifelse(.$padj <= 0.05, "Yes", "No"))

T1_HG64_dds_res_sig$Category <- factor(T1_HG64_dds_res_sig$Category, levels = T1_HG64_dds_res_sig$Category)
ggplot(T1_HG64_dds_res_sig, aes(x = Category, y = log2FoldChange, fill = Signif)) + 
    geom_point(size = 4, shape = 21, colour = "#40916c") +
    geom_hline(yintercept = 0, linetype = "dashed", size = 1) + # add ref line
    scale_fill_manual(values = c("transparent", "#40916c")) + # change color manually
    theme_bw() + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5), # center title
          axis.text = element_text(size = 8),
          axis.title = element_text(size = 12,face = "bold")) +
    ggtitle("HG64_T1 vs. Susceptible") + # add title
    xlab(NULL) + # remove xlab title
    ylab("log2 Fold Change") +
    # guides(fill = FALSE) + # remove legend 
    coord_flip()

ggsave("HG64_T1_vs_Susceptible.pdf", width = 5, height = 3) # save pdf



# T1_HF12
design_T1_HF12 <- design[grep("T1", design$Treatment),] %>%
    .[-which(.$Treatment %in% c("HG64T1")),]
design_T1_HF12$Property <- factor(design_T1_HF12$Property, levels = c("Susceptible", "Resistant"))
pathway_T1_HF12 <- pathway[, rownames(design_T1_HF12)]
pathway_T1_HF12_dds <- DESeqDataSetFromMatrix(pathway_T1_HF12, design_T1_HF12, design = ~ Property)
keep <- rowSums(counts(pathway_T1_HF12_dds) >= 10) >= 3
pathway_T1_HF12_dds <- pathway_T1_HF12_dds[keep, ]
vs_pathway_T1_HF12_dds <- rlog(pathway_T1_HF12_dds, blind = FALSE)

# DEGs analysis
T1_HF12_dds <- DESeq(pathway_T1_HF12_dds)
resultsNames(T1_HF12_dds)
# [1] "Intercept"                         "Property_Resistant_vs_Susceptible"
T1_HF12_dds_res <- results(T1_HF12_dds)
summary(T1_HF12_dds_res)
# out of 147 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 0, 0%
# LFC < 0 (down)     : 2, 1.4%
# outliers [1]       : 0, 0%
# low counts [2]     : 0, 0%
# (mean count < 5)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

T1_HF12_dds_res_sig <- subset(T1_HF12_dds_res, (log2FoldChange > 1 | log2FoldChange < -1)) %>% 
    as.data.frame() %>%
    arrange(desc(log2FoldChange)) %>%
    replace(is.na(.),1) %>%
    rownames_to_column() %>% 
    rename(Category = rowname) %>%
    mutate(Signif = ifelse(.$padj <= 0.05, "Yes", "No"))
T1_HF12_dds_res_sig$Category <- factor(T1_HF12_dds_res_sig$Category, levels = T1_HF12_dds_res_sig$Category)
ggplot(T1_HF12_dds_res_sig, aes(x = Category, y = log2FoldChange, fill = Signif)) + 
    geom_point(size = 4, shape = 21, colour = "#40916c") +
    geom_hline(yintercept = 0, linetype = "dashed", size = 1) + # add ref line
    scale_fill_manual(values = c("transparent", "#40916c")) + # change color manually
    theme_bw() + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5), # center title
          axis.text = element_text(size = 8),
          axis.title = element_text(size = 12,face = "bold")) +
    ggtitle("HF12_T1 vs. Susceptible") + # add title
    xlab(NULL) + # remove xlab title
    ylab("log2 Fold Change") +
    # guides(fill = FALSE) + # remove legend 
    coord_flip()
ggsave("HF12_T1_vs_Susceptible.pdf", width = 5, height = 3) # save pdf

# T2C_HG64
design_T2C_HG64 <- design[grep("T2C", design$Treatment),] %>%
    .[-which(.$Treatment == "HF12T2C"),]
design_T2C_HG64$Property <- factor(design_T2C_HG64$Property, levels = c("Susceptible", "Resistant"))
pathway_T2C_HG64 <- pathway[, rownames(design_T2C_HG64)]
pathway_T2C_HG64_dds <- DESeqDataSetFromMatrix(pathway_T2C_HG64, design_T2C_HG64, design = ~ Property)
keep <- rowSums(counts(pathway_T2C_HG64_dds) >= 10) >= 3
pathway_T2C_HG64_dds <- pathway_T2C_HG64_dds[keep, ]
vs_pathway_T2C_HG64_dds <- rlog(pathway_T2C_HG64_dds, blind = TRUE)

# DEGs
T2C_HG64_dds <- DESeq(pathway_T2C_HG64_dds)
resultsNames(T2C_HG64_dds)
# [1] "Intercept"                         "Property_Resistant_vs_Susceptible"
T2C_HG64_dds_res <- results(T2C_HG64_dds)
summary(T2C_HG64_dds_res)

T2C_HG64_dds_res_sig <- subset(T2C_HG64_dds_res, (log2FoldChange > 1 | log2FoldChange < -1)) %>% 
    as.data.frame() %>%
    arrange(desc(log2FoldChange)) %>%
    replace(is.na(.),1) %>%
    rownames_to_column() %>% 
    rename(Category = rowname) %>%
    mutate(Signif = ifelse(.$padj <= 0.05, "Yes", "No"))
T2C_HG64_dds_res_sig$Category <- factor(T2C_HG64_dds_res_sig$Category, levels = T2C_HG64_dds_res_sig$Category)
ggplot(T2C_HG64_dds_res_sig, aes(x = Category, y = log2FoldChange, fill = Signif)) + 
    geom_point(size = 4, shape = 21, colour = "#40916c") +
    geom_hline(yintercept = 0, linetype = "dashed", size = 1) + # add ref line
    scale_fill_manual(values = c("transparent", "#40916c")) + # change color manually
    theme_bw() + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5), # center title
          axis.text = element_text(size = 8),
          axis.title = element_text(size = 12,face = "bold")) +
    ggtitle("HG64_T2C vs. Susceptible") + # add title
    xlab(NULL) + # remove xlab title
    ylab("log2 Fold Change") +
    # guides(fill = FALSE) + # remove legend 
    coord_flip()
ggsave("HG64_T2C_vs_Susceptible.pdf", width = 5, height = 3) # save pdf

# T2C_HF12
design_T2C_HF12 <- design[grep("T2C", design$Treatment),] %>%
    .[-which(.$Treatment == "HG64T2C"),]
design_T2C_HF12$Property <- factor(design_T2C_HF12$Property, levels = c("Susceptible", "Resistant"))
pathway_T2C_HF12 <- pathway[, rownames(design_T2C_HF12)]
pathway_T2C_HF12_dds <- DESeqDataSetFromMatrix(pathway_T2C_HF12, design_T2C_HF12, design = ~ Property)
keep <- rowSums(counts(pathway_T2C_HF12_dds) >= 10) >= 3
pathway_T2C_HF12_dds <- pathway_T2C_HF12_dds[keep, ]
vs_pathway_T2C_HF12_dds <- rlog(pathway_T2C_HF12_dds, blind = TRUE)
plotPCA(vs_pathway_T2C_HF12_dds, intgroup = "Treatment")
# DEGs analysis
T2C_HF12_dds <- DESeq(pathway_T2C_HF12_dds)
resultsNames(T2C_HF12_dds)
# [1] "Intercept"                         "Property_Resistant_vs_Susceptible"
T2C_HF12_dds_res <- results(T2C_HF12_dds)
summary(T2C_HF12_dds_res)

T2C_HF12_dds_res_sig <- subset(T2C_HF12_dds_res, (log2FoldChange > 1 | log2FoldChange < -1)) %>% 
    as.data.frame() %>%
    arrange(desc(log2FoldChange)) %>%
    replace(is.na(.),1) %>%
    rownames_to_column() %>% 
    rename(Category = rowname) %>%
    mutate(Signif = ifelse(.$padj <= 0.05, "Yes", "No"))
T2C_HF12_dds_res_sig$Category <- factor(T2C_HF12_dds_res_sig$Category, levels = T2C_HF12_dds_res_sig$Category)
ggplot(T2C_HF12_dds_res_sig, aes(x = Category, y = log2FoldChange, fill = Signif)) + 
    geom_point(size = 4, shape = 21, colour = "#40916c") +
    geom_hline(yintercept = 0, linetype = "dashed", size = 1) + # add ref line
    scale_fill_manual(values = c("transparent", "#40916c")) + # change color manually
    theme_bw() + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5), # center title
          axis.text = element_text(size = 8),
          axis.title = element_text(size = 12,face = "bold")) +
    ggtitle("HF12_T2C vs. Susceptible") + # add title
    xlab(NULL) + # remove xlab title
    ylab("log2 Fold Change") +
    # guides(fill = FALSE) + # remove legend 
    coord_flip()
ggsave("HF12_T2C_vs_Susceptible.pdf", width = 5, height = 3) # save pdf

# Then use Adobe illustator to merge the plot to Fig. 5