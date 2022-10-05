library(ggpubr)
library(car)
library(patchwork)
my_theme()
packageVersion("ggpubr")
# [1] ‘0.4.0’
sus <- read.table("A57_HG70.txt", header=T, row.names= 1, sep="\t")
resis <- read.table("HF12_HG64.txt", header=T, row.names= 1, sep="\t")

shapiro.test(sus$Disease_index)
# Shapiro-Wilk normality test
# 
# data:  sus$Disease_index
# W = 0.72731, p-value = 2.523e-09
shapiro.test(resis$Disease_index)
# Shapiro-Wilk normality test
# 
# data:  resis$Disease_index
# W = 0.69677, p-value = 4.849e-16

comparison_1 <- list(c("Sterile soil", "Natural soil"))
p1 <- ggboxplot(sus, x = "Treatment", y = "Disease_index", 
               color = "Treatment", palette = "npg",
               add = "jitter", add.params = list(fill="white", width = 0.1),
               facet.by = "Cultivar", ncol = 2, width = 0.8) +
    labs(y="Disease index") +
    rremove("xlab") +
    rremove("x.text") +
    rremove("x.ticks") +
    stat_compare_means(comparisons = comparison_1, method = "wilcox", 
                       method.args = list(exact = FALSE), label = "p.signif")
p1

comparison_2 <- list(c("Sterile soil", "Natural soil"), 
                     c("Sterile soil", "Natural soil extract"),
                     c("Natural soil", "Natural soil extract"))
p2 <- ggboxplot(resis, x = "Treatment", y = "Disease_index", 
                color = "Treatment", palette = "npg",
                add = "jitter", add.params = list(fill="white", width = 0.1),
                facet.by = "Cultivar", ncol = 2, width = 0.8) +
    labs(y="Disease index") +
    rremove("xlab") +
    rremove("x.text") +
    rremove("x.ticks") +
    rremove("ylab") +
    stat_compare_means(comparisons = comparison_2, method = "wilcox",
                       method.args = list(exact = FALSE), label = "p.signif")
p2
p3 <- p1 + p2
p3
ggsave("Fig1.pdf", p3, width = 9, height = 5.5)
# Then use Adobe illustrator to adjust the figure legend