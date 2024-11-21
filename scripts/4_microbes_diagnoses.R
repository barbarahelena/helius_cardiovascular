## Associations microbes and new diagnoses

library(rio)
library(phyloseq)
library(haven)
library(tidyverse)
library(tableone)
library(ggsci)
library(broom)
library(patchwork)

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

theme_bar <- function() {
    (theme_Publication() +
        theme(plot.background = element_blank(), axis.text.x = element_blank(),
              axis.line.y = element_blank(), axis.ticks.y = element_blank(),
              axis.line.x = element_blank(), axis.ticks.x = element_blank(),
              panel.grid = element_blank()))
}

#### Data ####
df <- readRDS("data/clinicaldata_wide.RDS") 
mb <- readRDS("data/phyloseq_rarefied_cleaned.RDS")
tax <- readRDS("data/taxtable_rarefied_cleaned.RDS")

#### Preprocessing microbiome data ####
mb <- as.data.frame(t(as(mb@otu_table, "matrix")))
tk <- apply(mb, 2, function(x) sum(x > 5) > (0.2*length(x)))
mb <- mb[,tk]
head(mb)
mb <- mb %>% mutate(across(everything(.), ~log10(.x+1)))
mb$sampleID_baseline <- rownames(mb)
mbclin <- left_join(mb, df, by = "sampleID_baseline")
dim(mbclin)
head(mbclin)[1:5,1:5]

#### Associations with new diagnoses ####
names(mbclin)[1:ncol(mb)-1]
mbclin <- mbclin %>% mutate(
    across(c("DM_new", "HT_new", "MetSyn_new"), ~fct_recode(.x, "1" = "Yes", "0" = "No"))
)
res <- c()
for(a in 1:(ncol(mb)-1)) { # minus 1 because of sampleID_baseline var
    mbclin$asv <- mbclin[[a]]
    asvname <- colnames(mbclin)[a]
    print(asvname)
    # run models for each diagnosis while excluding participants with baseline diagnoses
    dm <- glm(DM_new ~ asv + Age_baseline, 
              data = mbclin %>% filter(DM_baseline == "No"), family = "binomial")
    ht <- glm(HT_new ~ asv  + Age_baseline, 
              data = mbclin %>% filter(HT_BPMed_baseline == "No"), family = "binomial")
    metsyn <- glm(MetSyn_new ~ asv + Age_baseline, 
                  data = mbclin %>% filter(MetSyn_baseline == "No"), family = "binomial")
    # extract estimates for variable sex
    dm <- tidy(dm, conf.int=TRUE, exponentiate = TRUE)[2,]
    ht <- tidy(ht, conf.int=TRUE, exponentiate = TRUE)[2,]
    metsyn <- tidy(metsyn, conf.int=TRUE, exponentiate = TRUE)[2,]
    # define rows
    row1 <- c(asvname, "diabetes", dm$estimate, dm$conf.low, dm$conf.high, dm$p.value)
    row2 <- c(asvname, "hypertension", ht$estimate, ht$conf.low, ht$conf.high, ht$p.value)
    row3 <- c(asvname, "metsyn", metsyn$estimate, metsyn$conf.low, metsyn$conf.high, metsyn$p.value)
    names(row1) <- c("ASV", "diagnosis", "OR", "lower", "upper", "pvalue")
    names(row2) <- c("ASV", "diagnosis", "OR", "lower", "upper", "pvalue")
    names(row3) <- c("ASV", "diagnosis", "OR", "lower", "upper", "pvalue")
    # bind
    res <- rbind(res, row1, row2, row3)
}
res <- as.data.frame(res)
res2 <- left_join(res, tax, by = 'ASV') %>% mutate(across(c(3:6), as.numeric))

dm <- res2 %>% filter(diagnosis == "diabetes") %>% arrange(pvalue) %>%
    mutate(padj = p.adjust(pvalue, "fdr"),
           Tax = as.factor(make.unique(Tax)),
           Tax = forcats::fct_rev(forcats::fct_inorder(Tax))
           ) %>% # calc FDR adj p value
    filter(padj < 0.05) # filter for FDR < 0.05
ht <- res2 %>% filter(diagnosis == "hypertension") %>% arrange(pvalue) %>%
    mutate(padj = p.adjust(pvalue, "fdr"),
           Tax = as.factor(make.unique(Tax)),
           Tax = forcats::fct_rev(forcats::fct_inorder(Tax))
           ) %>%  # calc FDR adj p value
    filter(padj < 0.05) # filter for FDR < 0.05
metsyn <- res2 %>% filter(diagnosis == "metsyn") %>% arrange(pvalue) %>%
    mutate(padj = p.adjust(pvalue, "fdr"),
           Tax = as.factor(make.unique(Tax)),
           Tax = forcats::fct_rev(forcats::fct_inorder(Tax))
           ) %>%  # calc FDR adj p value
    filter(padj < 0.05) # filter for FDR < 0.05
dim(dm)
dim(ht)
dim(metsyn)

summary(as.factor(dm$Family))
summary(as.factor(ht$Family))

coltab <- NULL
pal <- colorRampPalette(ggsci::pal_bmj()(9))
families <- c("Akkermansiaceae", "Anaerovoracaceae", "Bacteroidaceae", "Bifidobacteriaceae",
              "Christensenellaceae", "Clostridiaceae", "Desulfovibrionaceae", "Erysipelatoclostridiaceae",
              "Lachnospiraceae", "Marinifilaceae", "Monoglobaceae", "Oscillospiraceae", "Peptostreptococcaceae",
              "Rikenellaceae", "Ruminococcaceae", "Tannerellaceae", "Veillonellaceae"
              ) # to make colors in line with the registry analyses
coltab <- setNames(pal(length(families)), families)
families2 <- c("UCG-010", "Eggerthellaceae", "[Eubacterium] coprostanoligenes group", "Butyricicoccaceae",
               "Sutterellaceae", "Streptococcaceae", "Unknown") # not occuring in registry analyses
famtot <- c(families, families2)
pal2 <- colorRampPalette(ggsci::pal_futurama()(8)[c(1,3,4,7,8)])
coltab2 <- setNames(c(pal2(length(families2)-1), "grey80"), families2)
coltab <- coltab[which((names(coltab) %in% c(dm$Family, ht$Family, metsyn$Family)))]
coltab <- c(coltab, coltab2)

dm <- dm %>% mutate(Family = case_when(is.na(Family) ~ "Unknown", .default = Family),
                    Family = factor(Family, levels = famtot))
ht <- ht %>% mutate(Family = case_when(is.na(Family) ~ "Unknown", .default = Family),
                    Family = factor(Family, levels = famtot))
metsyn <- metsyn %>% mutate(Family = case_when(is.na(Family) ~ "Unknown", .default = Family),
                            Family = factor(Family, levels = famtot))

write.csv(dm, "results/diabetesassociations.csv")
write.csv(metsyn, "results/metsynassociationcs.csv")
write.csv(ht, "results/hypertensionassociations.csv")

(bar <- ggplot(dm, aes(y = Tax, x = 0, fill = Family)) +
        scale_fill_manual(values = coltab, drop = FALSE) +
        geom_tile(show.legend = TRUE) +
        labs(x = "", y = "", fill = "") +
        guides(fill = guide_legend(ncol = 1)) +
        theme_bar())

(pldm <- ggplot(dm, aes(x = OR, y = Tax)) +
            geom_vline(aes(xintercept = 1), color = "darkgrey") +
            geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.5,
                           color = pal_bmj()(1)[1]) +
            geom_point(color = pal_bmj()(1)[1]) +
            scale_x_continuous(n.breaks = 6, limits = c(0.25, 3.0)) +
            labs(title = "Diabetes", y = "", x = "OR per log10-increase ASV") +
            theme_Publication() +
            theme(axis.text.y = element_blank()))

(dmplot <- bar + plot_spacer() + pldm + 
    plot_layout(widths = c(0.1, -0.1, 0.7), axis_titles = "collect_y"))

(bar <- ggplot(metsyn, aes(y = Tax, x = 0, fill = Family)) +
    scale_fill_manual(values = coltab, drop = FALSE) +
    geom_tile(show.legend = TRUE) +
    labs(x = "", y = "", fill = "") +
    guides(fill = guide_legend(ncol = 1)) +
    theme_bar())

(plms <- ggplot(metsyn, aes(x = OR, y = Tax)) +
        geom_vline(aes(xintercept = 1), color = "darkgrey") +
        geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.5,
                       color = pal_bmj()(4)[4]) +
        geom_point(color = pal_bmj()(4)[4]) +
        scale_x_continuous(n.breaks = 6, limits = c(0.25, 3.0)) +
        labs(title = "MetSyn", y = "", x = "OR per log10-increase ASV") +
        theme_Publication() +
        theme(axis.text.y = element_blank()))

(msplot <- bar + plot_spacer() + plms + 
    plot_layout(widths = c(0.1, -0.1, 0.7), axis_titles = "collect_y"))
  
(bar <- ggplot(ht, aes(y = Tax, x = 0, fill = Family)) +
        scale_fill_manual(values = coltab, drop = FALSE) +
        geom_tile(show.legend = TRUE) +
        labs(x = "", y = "", fill = "") +
        guides(fill = guide_legend(ncol = 1)) +
        theme_bar())

(plht <- ggplot(ht, aes(x = OR, y = Tax)) +
        geom_vline(aes(xintercept = 1), color = "darkgrey") +
        geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.5,
                       color = pal_bmj()(3)[3]) +
        geom_point(color = pal_bmj()(3)[3]) +
        scale_x_continuous(n.breaks = 6, limits = c(0.25, 3.0)) +
        labs(title = "Hypertension", y = "", x = "OR per log10-increase ASV") +
        theme_Publication()+
        theme(axis.text.y = element_blank()))

(htplot <- bar + plot_spacer() + plht + 
    plot_layout(widths = c(0.1, -0.1, 0.7), axis_titles = "collect_y"))

dmplot / msplot / htplot +
    plot_layout(guides = "collect", nrow = 3, heights = c(0.45, 1.0, 1.1)) +
    plot_annotation(tag_levels = list(c("A", "", "B", "", "C", ""))) &
    theme(plot.tag = element_text(face = "bold"),
          legend.key.size= unit(0.4, "cm"),
          legend.text = element_text(size = rel(1.0)))

ggsave("results/diagnosispred_long.pdf", width = 14, height = 26)

#### Adjustment for PPI ####
#### Associations with new diagnoses ####
names(mbclin)[1:ncol(mb)-1]
summary(mbclin$PPI_baseline)
res <- c()
for(a in 1:(ncol(mb)-1)) { # minus 1 because of sampleID_baseline var
    mbclin$asv <- mbclin[[a]]
    asvname <- colnames(mbclin)[a]
    print(asvname)
    # run models for each diagnosis while excluding participants with baseline diagnoses
    dm <- glm(DM_new ~ asv + Age_baseline + PPI_baseline, 
              data = mbclin %>% filter(DM_baseline == "No"), family = "binomial")
    ht <- glm(HT_new ~ asv  + Age_baseline + PPI_baseline, 
              data = mbclin %>% filter(HT_BPMed_baseline == "No"), family = "binomial")
    metsyn <- glm(MetSyn_new ~ asv + Age_baseline + PPI_baseline, 
                  data = mbclin %>% filter(MetSyn_baseline == "No"), family = "binomial")
    # extract estimates for variable sex
    dm <- tidy(dm, conf.int=TRUE, exponentiate = TRUE)[2,]
    ht <- tidy(ht, conf.int=TRUE, exponentiate = TRUE)[2,]
    metsyn <- tidy(metsyn, conf.int=TRUE, exponentiate = TRUE)[2,]
    # define rows
    row1 <- c(asvname, "diabetes", dm$estimate, dm$conf.low, dm$conf.high, dm$p.value)
    row2 <- c(asvname, "hypertension", ht$estimate, ht$conf.low, ht$conf.high, ht$p.value)
    row3 <- c(asvname, "metsyn", metsyn$estimate, metsyn$conf.low, metsyn$conf.high, metsyn$p.value)
    names(row1) <- c("ASV", "diagnosis", "OR", "lower", "upper", "pvalue")
    names(row2) <- c("ASV", "diagnosis", "OR", "lower", "upper", "pvalue")
    names(row3) <- c("ASV", "diagnosis", "OR", "lower", "upper", "pvalue")
    # bind
    res <- rbind(res, row1, row2, row3)
}
res <- as.data.frame(res)
res2 <- left_join(res, tax, by = 'ASV') %>% mutate(across(c(3:6), as.numeric))

dm <- res2 %>% filter(diagnosis == "diabetes") %>% arrange(pvalue) %>%
    mutate(padj = p.adjust(pvalue, "fdr"),
           Tax = as.factor(make.unique(Tax)),
           Tax = forcats::fct_rev(forcats::fct_inorder(Tax))
    ) %>% # calc FDR adj p value
    filter(padj < 0.05) # filter for FDR < 0.05
ht <- res2 %>% filter(diagnosis == "hypertension") %>% arrange(pvalue) %>%
    mutate(padj = p.adjust(pvalue, "fdr"),
           Tax = as.factor(make.unique(Tax)),
           Tax = forcats::fct_rev(forcats::fct_inorder(Tax))
    ) %>%  # calc FDR adj p value
    filter(padj < 0.05) # filter for FDR < 0.05
metsyn <- res2 %>% filter(diagnosis == "metsyn") %>% arrange(pvalue) %>%
    mutate(padj = p.adjust(pvalue, "fdr"),
           Tax = as.factor(make.unique(Tax)),
           Tax = forcats::fct_rev(forcats::fct_inorder(Tax))
    ) %>%  # calc FDR adj p value
    filter(padj < 0.05) # filter for FDR < 0.05
dim(dm)
dim(ht)
dim(metsyn)
