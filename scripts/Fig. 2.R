# Fig. 2A
library(ggpubr)
library(car)
my_theme()
design = read.table("metadata.txt", header=T, row.names= 1, sep="\t")
alpha = read.table("alpha-diversity.tsv", header=T, row.names= 1, sep="\t")
index = merge(alpha, design,by="row.names")
# with(index,leveneTest(shannon,Treatment))
# Warning in leveneTest.default(shannon, Treatment): Treatment coerced to factor.
# Levene's Test for Homogeneity of Variance (center = median)
#       Df F value Pr(>F)
# group 19  1.1398 0.3525
#       40
# p>0.05: Meet the homogeneity of variance assumption

shapiro.test(alpha$shannon)

# Shapiro-Wilk normality test

# data:  alpha$shannon
# W = 0.81788, p-value = 3.879e-07


# x axis: Treatment
p <- ggboxplot(index, x = "Cultivar", y = "shannon",color = "Property", palette = "npg",
               facet.by = "Time", add = "jitter", ncol = 2, width = 0.6) +
    labs(y="Shannon") +
    rremove("xlab") + 
    stat_compare_means(label.y.npc = 0.2, label.x.npc = 0.2)

ggpar(p, ylim = c(2,12))


# Fig. 2B
# BiocManager::install("yjiakang/microVisu") 
# run the above line if microVisu not installed
library(microVisu)
fileDir <- "D:/phd/tomato_micr/16S/results_no_primer/dada2/R" # please choose your own data dir
meta <- paste0(fileDir, "/metadata.txt")
my_theme()
# weighted unifrac distance
weig <- paste0(fileDir, "/weighted_unifrac.tsv")
pcoaPlot(weig, meta, distType = "Weighted Unifrac", classForColor = "Time",
         classForShape = "Property", col = "Dark2")

