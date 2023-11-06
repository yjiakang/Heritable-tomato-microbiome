function() {
    library(ggplot2)
    library(ggthemr)
    theme_set(theme_classic() +
                  theme(axis.text = element_text(size = 8, face = "bold"), 
                        axis.title = element_text(size = 12, face = "bold"), 
                        plot.margin = unit(rep(1, 4), "cm"),
                        panel.grid = element_blank()))
    }
