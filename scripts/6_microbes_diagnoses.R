## Associations microbes and new diagnoses

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

## Functions
log_group <- function(df, dfname, writetable = FALSE, figure = FALSE){
    theme_Publication <- function(base_size=12, base_family="sans") {
        library(grid)
        library(ggthemes)
        (theme_foundation(base_size=base_size, base_family=base_family)
            + theme(plot.title = element_text(face = "bold",
                                              size = rel(1.0), hjust = 0.5),
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
                    axis.ticks.x = element_line(),
                    axis.ticks.y = element_blank(),
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
    
    if(writetable == FALSE & figure == FALSE){
        print("No output set!")
    }    else{
        ## Variable selection
        dfsub <- df %>% select(CKD_group, Age, Sex, BMI, HT, CKDEPI, Ethnicity, FramRisk, tail(names(.), 20)) %>% 
            mutate(CKD_log=ifelse(CKD_group=="Controls", 0, 1))
        
        ## Models
        res_log <- c()
        for (i in c((ncol(dfsub)-19):ncol(dfsub)-1)){
            dfsub$met <- NULL    
            dfsub$met <- log10(dfsub[,i][[1]])
            m0 <- glm(CKD_log ~ scale(met), data = dfsub, family = "binomial")
            m1 <- glm(CKD_log ~ scale(met) + Age + Sex + BMI + HT, data=dfsub, family="binomial")
            
            metname <- colnames(dfsub)[i]
            m0 <- tidy(m0, conf.int=T)[2,]
            m1 <- tidy(m1, conf.int = T)[2,]
            
            resRow <- cbind(metname, exp(m0$estimate), exp(m0$conf.low), exp(m0$conf.high), m0$p.value,
                            exp(m1$estimate), exp(m1$conf.low), exp(m1$conf.high), m1$p.value)
            colnames(resRow) <- c("Metabolite", 
                                  "m0-est", "m0-l95", "m0-u95", "m0-p", 
                                  "m1-est", "m1-l95", "m1-u95", "m1-p")
            res_log <- rbind(res_log, resRow)
            dfsub$met <- NULL 
        }
        
        reslog <- as.data.frame(res_log)
        afronden2 <- function(x) return(as.numeric(format(round(x, 2),2)))
        afronden5 <- function(x) return(as.numeric(format(round(x, 5),5)))
        reslog2 <- reslog %>% 
            mutate_at(c(2:9), as.character) %>% 
            mutate_at(c(2:9), as.numeric) %>% 
            mutate_at(c(2:4, 6:8), afronden2) %>% 
            mutate_at(c(5,9), afronden5) %>% 
            mutate(
                `m0-q` = p.adjust(`m0-p`, 'fdr'),
                `m1-q` = p.adjust(`m1-p`, 'fdr')
            )
        
        ## Output
        if(figure == TRUE){
            labslist <- c("Unadjusted", "Adjusted")
            reslong <- reslog2 %>% 
                pivot_longer(c(2:9), names_to=c("model", "cat"), 
                             names_prefix="m", 
                             names_sep='-',
                             values_to="value") %>% 
                pivot_wider(names_from = cat, values_from = value) %>% 
                mutate(model = factor(model, levels = c("0", "1"), 
                                      labels = labslist),
                       Metabolite = factor(Metabolite, levels = colnames(dfsub)[9:28]),
                       Metabolite = fct_rev(Metabolite))
            
            ylab <- "OR for new diagnosis per log10 increase"
            colors <- c(pal_jco()(4)[1], pal_jco()(4)[4])
            pl <- ggplot(reslong, aes(x=Metabolite,y=est, color=model)) +
                geom_hline(yintercept = 1, color = "grey40") +
                geom_point(position=position_dodge(-0.5)) +
                geom_errorbar(aes(ymin=l95,ymax=u95,width=.3), position=position_dodge(-0.5)) +
                expand_limits(y=0)+
                scale_y_log10(breaks = c(0.5, 0.75, 1.0, 1.5, 2.0))+
                theme_Publication()+
                labs(title = "", x = "", y = ylab) +
                scale_color_manual(values = colors) +
                coord_flip()
            ggsave(pl, filename = str_c("results/logreg/CKD/", dfname, "_logreg.pdf"), device = "pdf", width = 8)
            ggsave(pl, filename = str_c("results/logreg/CKD/", dfname, "_logreg.svg"), device = "svg", width = 8)
        } 
        
        if(writetable == TRUE){
            openxlsx::write.xlsx(reslog2, file.path("results/logreg/CKD", str_c(dfname,"_logreg.xlsx")))
        }
    }
}

#### Data ####
mb <- readRDS("data/16S/phyloseq_paired16s.RDS")
mb <- as.data.frame(t(as(mb@otu_table, "matrix")))
head(mb)
mbsel <- mb[, c("ASV_211", "ASV_120")]
colnames(mbsel) <- c("Parabacteroides_distasonis", "Streptococcus_spp.")
mbsel <- mbsel %>% mutate(across(everything(.), ~log10(.x+1)))
mbsel$ID <- rownames(mbsel)

mbmet <- left_join(mbsel, metba, by = "ID")
dim(mbmet)

df <- readRDS("data/clinicaldata_long.RDS") %>% filter("baseline")

#### Associations with new diagnoses ####
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
