## Associations metabolomics and microbes

library(rio)
library(haven)
library(tidyverse)
library(tableone)
library(ggsci)

theme_Publication <- function(base_size=12, base_family="sans") {
    library(grid)
    library(ggthemes)
    (theme_foundation(base_size=base_size, base_family=base_family)
        + theme(plot.title = element_text(face = "bold",
                                          size = rel(0.8), hjust = 0.5),
                #family = 'Helvetica'
                text = element_text(),
                panel.background = element_rect(colour = NA),
                plot.background = element_rect(colour = NA),
                panel.border = element_rect(colour = NA),
                axis.title = element_text(face = "bold",size = rel(0.8)),
                axis.title.y = element_text(angle=90,vjust =2),
                axis.title.x = element_text(vjust = -0.2),
                axis.text = element_text(), 
                axis.line.x = element_line(colour="black"),
                axis.line.y = element_line(colour="black"),
                axis.ticks.x = element_line(),
                axis.ticks.y = element_line(),
                panel.grid.major = element_line(colour="#f0f0f0"),
                panel.grid.minor = element_blank(),
                legend.key = element_rect(colour = NA),
                legend.position = "right",
                # legend.direction = "horizontal",
                legend.key.size= unit(0.2, "cm"),
                legend.spacing  = unit(0, "cm"),
                # legend.title = element_text(face="italic"),
                plot.margin=unit(c(5,5,5,5),"mm"),
                strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                strip.text = element_text(face="bold"),
                plot.caption = element_text(face = "italic", size=rel(0.6))
        ))
} 








############### Associations with microbes ##################
mb <- readRDS("data/16S/phyloseq_paired16s.RDS")
mb <- as.data.frame(t(as(mb@otu_table, "matrix")))
head(mb)
mbsel <- mb[, c("ASV_211", "ASV_120")]
colnames(mbsel) <- c("Parabacteroides_distasonis", "Streptococcus_spp.")
mbsel <- mbsel %>% mutate(across(everything(.), ~log10(.x+1)))
mbsel$ID <- rownames(mbsel)

met <- readRDS("data/metabolomics/metabolomics_paired.RDS")
metba <- met[str_detect(rownames(met), "HELIBA"),!str_detect(colnames(met), "X-[0-9]*")]
metba <- as.data.frame(metba)
metba <- metba %>% mutate(across(everything(.), ~log10(.x+1)))
dim(metba)
metba$ID <- rownames(metba)

mbmet <- left_join(mbsel, metba, by = "ID")
dim(mbmet)

names(mbmet)[1:10]
res <- c()
mbmet$met <- NULL
for(a in 4:(ncol(mbmet))) {
    mbmet$met <- mbmet[[a]]
    print(a)
    print(names(mbmet)[a])
    cor1 <- cor.test(mbmet$met, y = mbmet[,1], method = "spearman")
    cor2 <- cor.test(mbmet$met, y = mbmet[,2], method = "spearman")
    row <- c(names(mbmet)[a], cor1$estimate[[1]], cor1$p.value, 
                                cor2$estimate[[1]], cor2$p.value)
    names(row) <- c("metabolite", "Parabacteroides_rho", "Parabacteroides_pvalue",
                                    "Streptococcus_rho", "Streptococcus_pvalue")
    res <- rbind(res, row)
}
res <- as.data.frame(res)
res <- res %>% mutate(across(c(2:5), as.numeric)) %>% 
    mutate(Parabacteroides_padj = p.adjust(Parabacteroides_pvalue, "fdr"),
           Streptococcus_padj = p.adjust(Streptococcus_pvalue, "fdr"))

res %>% filter(Parabacteroides_pvalue < 0.05 | Streptococcus_pvalue < 0.05) %>% nrow(.)
res %>% filter(Parabacteroides_pvalue < 0.05) %>% arrange(Parabacteroides_pvalue) %>% arrange(Parabacteroides_pvalue) %>% nrow(.)
res %>% filter(Parabacteroides_padj < 0.05 | Streptococcus_padj < 0.05) 

all(is.numeric(4:ncol(mbmet)))

para <- res %>% filter(Parabacteroides_pvalue < 0.01)
df <- mbmet
plotlist <- c()
for (a in para$metabolite) {
    print(a)
    df$met <- df[[a]]
    pl <- ggplot(df, mapping = aes(x = Parabacteroides_distasonis, y = met)) +
        geom_point(color = "royalblue4", alpha = 0.75) +
        geom_smooth(method = "lm", color = "firebrick3") +
        ggpubr::stat_cor(method = "spearman") +
        theme_Publication() +
        labs(x = "Parabacteroides distasonis (log10)", y = "concentration (log10(z-score))", 
             title = a)
    plotlist[[a]] <- pl
}
names(plotlist)[1:21]

ggpubr::ggarrange(plotlist = plotlist, nrow = 6, ncol = 4, labels = LETTERS[1:21])
ggsave("results/metabolomics/parabacteroides_corr.pdf", width = 15, height = 25)

strep <- res %>% filter(Streptococcus_padj < 0.05)
df <- mbmet
plotlist <- c()
for (a in strep$metabolite) {
    print(a)
    df$met <- df[[a]]
    pl <- ggplot(df, mapping = aes(x = Streptococcus_spp., y = met)) +
        geom_point(color = "royalblue4", alpha = 0.75) +
        geom_smooth(method = "lm", color = "firebrick3") +
        ggpubr::stat_cor(method = "spearman") +
        theme_Publication() +
        labs(x = "Streptococcus spp. (log10)", y = "concentration (log10(z-score))", 
             title = a)
    plotlist[[a]] <- pl
}
names(plotlist)

ggpubr::ggarrange(plotlist = plotlist, nrow = 1, ncol = 2, labels = LETTERS[1:2])
ggsave("results/metabolomics/strep_corr_padj.pdf", width = 8, height = 5)

strep <- res %>% filter(Streptococcus_pvalue < 0.01)
df <- mbmet
plotlist <- c()
for (a in strep$metabolite) {
    print(a)
    df$met <- df[[a]]
    pl <- ggplot(df, mapping = aes(x = Streptococcus_spp., y = met)) +
        geom_point(color = "royalblue4", alpha = 0.75) +
        geom_smooth(method = "lm", color = "firebrick3") +
        ggpubr::stat_cor(method = "spearman") +
        theme_Publication() +
        labs(x = "Streptococcus spp. (log10)", y = str_c("concentration (log10(z-score))"), 
             title = a)
    plotlist[[a]] <- pl
}
names(plotlist)

ggpubr::ggarrange(plotlist = plotlist, nrow = 9, ncol = 6, labels = c(LETTERS[1:26], 
                                                                      str_c("A", LETTERS[1:26]),
                                                                      str_c("B", LETTERS[1:2])))
ggsave("results/metabolomics/strep_corr.pdf", width = 15, height = 25)
