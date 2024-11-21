## Alpha and beta diversity 
## Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

## Libraries
library(tidyverse)
library(phyloseq)
library(mixOmics)
library(vegan)

theme_Publication <- function(base_size=14, base_family="sans") {
    library(grid)
    library(ggthemes)
    library(stringr)
    (theme_foundation(base_size=base_size, base_family=base_family)
        + theme(plot.title = element_text(face = "bold",
                                          size = rel(1.0), hjust = 0.5),
                text = element_text(),
                panel.background = element_rect(colour = NA),
                plot.background = element_rect(colour = NA),
                panel.border = element_rect(colour = NA),
                axis.title = element_text(face = "bold",size = rel(0.8)),
                axis.title.y = element_text(angle=90, vjust =2),
                axis.title.x = element_text(vjust = -0.2),
                axis.text = element_text(size = rel(0.7)),
                axis.text.x = element_text(angle = 0), 
                axis.line = element_line(colour="black"),
                axis.ticks = element_line(),
                panel.grid.major = element_line(colour="#f0f0f0"),
                panel.grid.minor = element_blank(),
                legend.key = element_rect(colour = NA),
                legend.position = "bottom",
                # legend.direction = "horizontal",
                legend.key.size= unit(0.2, "cm"),
                legend.spacing  = unit(0, "cm"),
                # legend.title = element_text(face="italic"),
                plot.margin=unit(c(10,5,5,5),"mm"),
                strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                strip.text = element_text(face="bold"),
                plot.caption = element_text(size = rel(0.5), face = "italic")
        ))
    
} 

#### Output folder ####
resultsfolder <- "results/alphadiversity"
dir.create(resultsfolder, showWarnings = FALSE)
resultsfolder <- "results/ordination"
dir.create(resultsfolder, showWarnings = FALSE)

#### Load data ####
phydata <- readRDS("data/phyloseq_rarefied_cleaned.RDS")
df_new <- readRDS("data/clinicaldata_wide.RDS")
tab <- as.data.frame(t(as(phydata@otu_table, 'matrix')))
tab_matrix <- t(as(phydata@otu_table, 'matrix'))

## Diversity metrics
# Shannon plots
shannon <- vegan::diversity(tab, index = 'shannon')
df_shan <- data.frame(sampleID_baseline = names(shannon), shannon = shannon)
df_shan <- left_join(df_shan, df_new, by = "sampleID_baseline")

## Species richness
specrich <- specnumber(tab)
dfspec <- data.frame(sampleID_baseline = names(specrich), richness = specrich)
dfspec <- left_join(dfspec, df_shan, by = "sampleID_baseline")

## Faith's PD
faith <- picante::pd(samp = tab_matrix, tree = phydata@phy_tree)
dffai <- as.data.frame(faith)
dffai$sampleID_baseline <- rownames(faith)
dffai <- left_join(dffai, dfspec, by = "sampleID_baseline")

dftot <- dffai %>%  
    mutate(
        DM_new = fct_recode(DM_new, "No diabetes" = "No", "New-onset diabetes" = "Yes"),
        HT_new = fct_recode(HT_new, "No hypertension" = "No", "New-onset hypertension" = "Yes"),
        MetSyn_new = fct_recode(MetSyn_new, "No MetSyn" = "No", "New-onset MetSyn" = "Yes")
    )

#### Baseline alpha diversity ####
(shanpl1 <- ggplot(data = dftot %>% filter(!is.na(DM_new)), aes(x = DM_new, y = shannon)) +
                geom_violin(aes(alpha = DM_new), fill = pal_bmj()(3)[1], show.legend = FALSE) +
                scale_alpha_manual(values = c(0.4, 1.0)) +
                geom_boxplot(fill = "white", width = 0.2) +
                labs(y = "Shannon index", x = "", alpha = "", fill = "") +
                stat_compare_means(comparisons = list(c("No diabetes", "New-onset diabetes")),
                                   tip.length = 0, hide.ns = TRUE,
                                   label = "p.signif", method = "wilcox.test") +
                theme_Publication() + 
     theme(axis.text.x = element_text(hjust = 1, angle = 30)))
(richpl1 <- ggplot(data = dftot %>% filter(!is.na(DM_new)), aes(x = DM_new, y = richness)) +
        geom_violin(aes(alpha = DM_new), fill = pal_bmj()(3)[1], show.legend = FALSE) +
        scale_alpha_manual(values = c(0.4, 1.0)) +
        geom_boxplot(fill = "white", width = 0.2) +
        labs(y = "Richness", x = "", alpha = "", fill = "") +
        stat_compare_means(comparisons = list(c("No diabetes", "New-onset diabetes")),
                           tip.length = 0, hide.ns = TRUE,
                           label = "p.signif", method = "wilcox.test") +
        theme_Publication() + 
        theme(axis.text.x = element_text(hjust = 1, angle = 30)))
(pdpl1 <- ggplot(data = dftot %>% filter(!is.na(DM_new)), aes(x = DM_new, y = PD)) +
        geom_violin(aes(alpha = DM_new), fill = pal_bmj()(3)[1], show.legend = FALSE) +
        scale_alpha_manual(values = c(0.4, 1.0)) +
        geom_boxplot(fill = "white", width = 0.2) +
        labs(y = "Faith's phylogenetic diversity", x = "", alpha = "", fill = "") +
        stat_compare_means(comparisons = list(c("No diabetes", "New-onset diabetes")),
                           tip.length = 0, hide.ns = TRUE,
                           label = "p.signif", method = "wilcox.test") +
        theme_Publication() + 
        theme(axis.text.x = element_text(hjust = 1, angle = 30)))

(shanpl2 <- ggplot(data = dftot %>% filter(!is.na(MetSyn_new)), aes(x = MetSyn_new, y = shannon)) +
        geom_violin(aes(alpha = MetSyn_new), fill = pal_bmj()(4)[4], show.legend = FALSE) +
        scale_alpha_manual(values = c(0.4, 1.0)) +
        geom_boxplot(fill = "white", width = 0.2) +
        labs(y = "Shannon index", x = "", alpha = "", fill = "") +
        stat_compare_means(comparisons = list(c("No MetSyn", "New-onset MetSyn")),
                           tip.length = 0, hide.ns = TRUE,
                           label = "p.signif", method = "wilcox.test") +
        theme_Publication() + 
        theme(axis.text.x = element_text(hjust = 1, angle = 30)))
(richpl2 <- ggplot(data = dftot %>% filter(!is.na(MetSyn_new)), aes(x = MetSyn_new, y = richness)) +
        geom_violin(aes(alpha = MetSyn_new), fill = pal_bmj()(4)[4], show.legend = FALSE) +
        scale_alpha_manual(values = c(0.4, 1.0)) +
        geom_boxplot(fill = "white", width = 0.2) +
        labs(y = "Richness", x = "", alpha = "", fill = "") +
        stat_compare_means(comparisons = list(c("No MetSyn", "New-onset MetSyn")),
                           tip.length = 0, hide.ns = TRUE,
                           label = "p.signif", method = "wilcox.test") +
        theme_Publication() + 
        theme(axis.text.x = element_text(hjust = 1, angle = 30)))
(pdpl2 <- ggplot(data = dftot %>% filter(!is.na(MetSyn_new)), aes(x = MetSyn_new, y = PD)) +
        geom_violin(aes(alpha = MetSyn_new), fill = pal_bmj()(4)[4], show.legend = FALSE) +
        scale_alpha_manual(values = c(0.4, 1.0)) +
        geom_boxplot(fill = "white", width = 0.2) +
        labs(y = "Faith's phylogenetic diversity", x = "", alpha = "", fill = "") +
        stat_compare_means(comparisons = list(c("No MetSyn", "New-onset MetSyn")),
                           tip.length = 0, hide.ns = TRUE,
                           label = "p.signif", method = "wilcox.test") +
        theme_Publication() + 
        theme(axis.text.x = element_text(hjust = 1, angle = 30)))

(shanpl3 <- ggplot(data = dftot %>% filter(!is.na(HT_new)), aes(x = HT_new, y = shannon)) +
        geom_violin(aes(alpha = HT_new), fill = pal_bmj()(3)[3], show.legend = FALSE) +
        scale_alpha_manual(values = c(0.4, 1.0)) +
        geom_boxplot(fill = "white", width = 0.2) +
        labs(y = "Shannon index", x = "", alpha = "", fill = "") +
        stat_compare_means(comparisons = list(c("No hypertension", "New-onset hypertension")),
                           tip.length = 0, hide.ns = TRUE,
                           label = "p.signif", method = "wilcox.test") +
        theme_Publication()+ 
        theme(axis.text.x = element_text(hjust = 1, angle = 30)))
(richpl3 <- ggplot(data = dftot %>% filter(!is.na(HT_new)), aes(x = HT_new, y = richness)) +
        geom_violin(aes(alpha = HT_new), fill = pal_bmj()(3)[3], show.legend = FALSE) +
        scale_alpha_manual(values = c(0.4, 1.0)) +
        geom_boxplot(fill = "white", width = 0.2) +
        labs(y = "Richness", x = "", alpha = "", fill = "") +
        stat_compare_means(comparisons = list(c("No hypertension", "New-onset hypertension")),
                           tip.length = 0, hide.ns = TRUE,
                           label = "p.signif", method = "wilcox.test") +
        theme_Publication() + 
        theme(axis.text.x = element_text(hjust = 1, angle = 30)))
(pdpl3 <- ggplot(data = dftot %>% filter(!is.na(HT_new)), aes(x = HT_new, y = PD)) +
        geom_violin(aes(alpha = HT_new), fill = pal_bmj()(3)[3], show.legend = FALSE) +
        scale_alpha_manual(values = c(0.4, 1.0)) +
        geom_boxplot(fill = "white", width = 0.2) +
        labs(y = "Faith's phylogenetic diversity", x = "", alpha = "", fill = "") +
        stat_compare_means(comparisons = list(c("No hypertension", "New-onset hypertension")),
                           tip.length = 0, hide.ns = TRUE,
                           label = "p.signif", method = "wilcox.test") +
        theme_Publication()+ 
        theme(axis.text.x = element_text(hjust = 1, angle = 30)))


(shanplots <- ggarrange(shanpl1, richpl1, pdpl1, shanpl2, richpl2, pdpl2, shanpl3, richpl3, pdpl3,
                            nrow = 3, ncol = 3, labels = LETTERS[c(1:3,5:7,9:11)]))
ggsave("results/alphadiversity/shannon_baseline_newdiagnoses.pdf", width = 12, height = 14)


#### Bray Curtis ####
df_new <- readRDS("data/clinicaldata_wide.RDS")
phydata <- readRDS("data/phyloseq_rarefied_cleaned.RDS")

#### Diabetes ####
dfdm <- df_new %>% filter(!is.na(DM_new))
phydm <- prune_samples(sample_names(phydata) %in% dfdm$sampleID_baseline, phydata)
tabdm <- as.data.frame(t(as(phydm@otu_table, 'matrix')))

braydm <- vegan::vegdist(tabdm, method = 'bray') # takes about 15 min
pcoorddm <- ape::pcoa(braydm, correction = "cailliez") # takes about 20 min
expl_variance_braydm <- pcoorddm$values$Rel_corr_eig * 100
write_lines(expl_variance_braydm, "results/ordination/expl_var_bray_dm.csv")

dbray <- pcoorddm$vectors[, c('Axis.1', 'Axis.2')]
dbray <- as.data.frame(dbray)
dbray$sampleID_baseline <- rownames(dbray)
bray_dm <- dbray <- left_join(dbray, dfdm, by = 'sampleID_baseline') %>%  # add metadata / covariates
    dplyr::select(BrayPCo1 = `Axis.1`, BrayPCo2 = `Axis.2`, everything(.))
dfdm <- dfdm %>% filter(sampleID_baseline %in% attributes(braydm)[["Labels"]])
dfdm <- dfdm[order(match(dfdm$sampleID_baseline, attributes(braydm)[["Labels"]])),]
all(dfdm$sampleID_baseline == attributes(braydm)[["Labels"]]) # TRUE

resdm <- adonis2(braydm ~ DM_new, data = dfdm) # PERMANOVA

#### MetSyn ####
dfms <- df_new %>% filter(!is.na(MetSyn_new))
phyms <- prune_samples(sample_names(phydata) %in% dfms$sampleID_baseline, phydata)
tabms <- as.data.frame(t(as(phyms@otu_table, 'matrix')))

brayms <- vegan::vegdist(tabms, method = 'bray') # takes about 5-10 min
pcoordms <- ape::pcoa(brayms, correction = "cailliez") # takes about 15 min
expl_variance_brayms <- pcoordms$values$Rel_corr_eig * 100
write_lines(expl_variance_brayms, "results/ordination/expl_var_bray_ms.csv")

dbray <- pcoordms$vectors[, c('Axis.1', 'Axis.2')]
dbray <- as.data.frame(dbray)
dbray$sampleID_baseline <- rownames(dbray)
bray_ms <- dbray <- left_join(dbray, dfms, by = 'sampleID_baseline') %>%  # add metadata / covariates
    dplyr::select(BrayPCo1 = `Axis.1`, BrayPCo2 = `Axis.2`, everything(.))
dfms <- dfms %>% filter(sampleID_baseline %in% attributes(brayms)[["Labels"]])
dfms <- dfms[order(match(dfms$sampleID_baseline, attributes(brayms)[["Labels"]])),]
all(dfms$sampleID_baseline == attributes(brayms)[["Labels"]]) # TRUE

resms <- adonis2(brayms ~ MetSyn_new, data = dfms) # PERMANOVA

#### Hypertension ####
dfht <- df_new %>% filter(!is.na(HT_new))
phyht <- prune_samples(sample_names(phydata) %in% dfht$sampleID_baseline, phydata)
tabht <- as.data.frame(t(as(phyht@otu_table, 'matrix')))

brayht <- vegan::vegdist(tabht, method = 'bray') # takes about 15 min
pcoordht <- ape::pcoa(brayht, correction = "cailliez") # takes about 20 min
expl_variance_brayht <- pcoordht$values$Rel_corr_eig * 100
write_lines(expl_variance_brayht, "results/ordination/expl_var_bray_ht.csv")

dbray <- pcoordht$vectors[, c('Axis.1', 'Axis.2')]
dbray <- as.data.frame(dbray)
dbray$sampleID_baseline <- rownames(dbray)
bray_hypertension <- dbray <- left_join(dbray, dfht, by = 'sampleID_baseline') %>%  # add metadata / covariates
    dplyr::select(BrayPCo1 = `Axis.1`, BrayPCo2 = `Axis.2`, everything(.))
dfht <- dfht %>% filter(sampleID_baseline %in% attributes(brayht)[["Labels"]])
dfht <- dfht[order(match(dfht$sampleID_baseline, attributes(brayht)[["Labels"]])),]
all(dfht$sampleID_baseline == attributes(brayht)[["Labels"]]) # TRUE

resht <- adonis2(brayht ~ HT_new, data = dfht) # PERMANOVA

#### Plots bray-curtis PCoA ####
(brayplotdm <- bray_dm %>%
     mutate(DM_new = fct_recode(DM_new, "No diabetes" = "No", "New-onset diabetes" = "Yes"),
            BrayPCo2 = BrayPCo2 * -1) %>% 
     ggplot(aes(BrayPCo1, BrayPCo2)) +
     stat_ellipse(geom = "polygon", aes(color = DM_new, fill = DM_new), type = "norm",
                  alpha = 0.1) +
     geom_point(aes(color = DM_new), size = 1, alpha = 0.7) +
     ggtitle("New-onset diabetes") +
     xlab(paste0('PCo1 (', round(expl_variance_braydm[1], digits = 1),'%)')) +
     ylab(paste0('PCo2 (', round(expl_variance_braydm[2], digits = 1),'%)')) +
     scale_color_manual(values = c("grey80", pal_bmj()(1))) +
     scale_fill_manual(values = c("grey80", pal_bmj()(1)), guide = "none") +
     scale_alpha_manual(guide = "none") +
     theme_Publication() +
     labs(color = "", alpha = "") +
     annotate("text", x= Inf, y = Inf, hjust = 1, vjust = 1,
              label = str_c("PERMANOVA: p = ", resdm$`Pr(>F)`, ", r2 = ", 
                            format(round(resdm$R2[1],3), nsmall = 3))
     ))

(brayplotms <- bray_ms %>%
        mutate(MetSyn_new = fct_recode(MetSyn_new, "No MetSyn" = "No", "New-onset MetSyn" = "Yes")) %>% 
        ggplot(aes(BrayPCo1, BrayPCo2)) +
        stat_ellipse(geom = "polygon", aes(color = MetSyn_new, fill = MetSyn_new), type = "norm",
                     alpha = 0.1) +
        geom_point(aes(color = MetSyn_new), size = 1, alpha = 0.7) +
        ggtitle("New-onset MetSyn") +
        xlab(paste0('PCo1 (', round(expl_variance_brayms[1], digits = 1),'%)')) +
        ylab(paste0('PCo2 (', round(expl_variance_brayms[2], digits = 1),'%)')) +
        scale_color_manual(values = c("grey80", pal_bmj()(4)[4])) +
        scale_fill_manual(values = c("grey80", pal_bmj()(4)[4]), guide = "none") +
        scale_alpha_manual(guide = "none") +
        theme_Publication() +
        labs(color = "", alpha = "") +
        annotate("text", x= Inf, y = Inf, hjust = 1, vjust = 1,
                 label = str_c("PERMANOVA: p = ", resms$`Pr(>F)`, ", r2 = ", 
                               format(round(resms$R2[1],3), nsmall = 3))
        ))

(brayplotht <- bray_hypertension %>%
        mutate(HT_new = fct_recode(HT_new, "No hypertension" = "No", "New-onset hypertension" = "Yes"),
               BrayPCo1 = BrayPCo1 * -1) %>% 
        ggplot(aes(BrayPCo1, BrayPCo2)) +
        stat_ellipse(geom = "polygon", aes(color = HT_new, fill = HT_new), type = "norm",
                     alpha = 0.1) +
        geom_point(aes(color = HT_new), size = 1, alpha = 0.7) +
        ggtitle("New-onset hypertension") +
        xlab(paste0('PCo1 (', round(expl_variance_brayht[1], digits = 1),'%)')) +
        ylab(paste0('PCo2 (', round(expl_variance_brayht[2], digits = 1),'%)')) +
        scale_color_manual(values = c("grey80", pal_bmj()(3)[3])) +
        scale_fill_manual(values = c("grey80", pal_bmj()(3)[3]), guide = "none") +
        scale_alpha_manual(guide = "none") +
        theme_Publication() +
        labs(color = "", alpha = "") +
        annotate("text", x= Inf, y = Inf, hjust = 1, vjust = 1,
                 label = str_c("PERMANOVA: p = ", resht$`Pr(>F)`, ", r2 = ", 
                               format(round(resht$R2[1],3), nsmall = 3))
        ))

(brayplots <- ggarrange(brayplotdm, brayplotms, brayplotht, nrow = 3, ncol = 1, labels = c("D", "H", "L")))
ggsave("results/ordination/brayplots.pdf", width = 6, height = 18)

#### Arrange total plot ####
ggarrange(shanplots, brayplots, widths = c(0.6, 0.3))
ggsave("results/diversityplots_diagnoses.pdf", width = 12, height = 14)

#### Shannon stratified for ethnicity ####
dftot <- dftot %>% filter(EthnicityTot != "Other") %>% droplevels(.) %>% filter(!is.na(EthnicityTot))

(shaneth1 <- ggplot(data = dftot %>% filter(!is.na(DM_new)), aes(x = DM_new, y = shannon)) +
        geom_violin(aes(alpha = DM_new), fill = pal_bmj()(3)[1], show.legend = FALSE) +
        scale_alpha_manual(values = c(0.4, 1.0)) +
        geom_boxplot(fill = "white", width = 0.2) +
        labs(y = "Shannon index", x = "", alpha = "", fill = "") +
        facet_wrap(~EthnicityTot) +
        stat_compare_means(tip.length = 0, 
                           label = "p.signif", method = "wilcox.test") +
        theme_Publication() +
        theme(axis.text.x = element_text(hjust = 1, angle = 45)))

(shaneth2 <- ggplot(data = dftot %>% filter(!is.na(MetSyn_new)), aes(x = MetSyn_new, y = shannon)) +
        geom_violin(aes(alpha = MetSyn_new), fill = pal_bmj()(3)[2], show.legend = FALSE) +
        scale_alpha_manual(values = c(0.4, 1.0)) +
        geom_boxplot(fill = "white", width = 0.2) +
        labs(y = "Shannon index", x = "", alpha = "", fill = "") +
        facet_wrap(~EthnicityTot) +
        stat_compare_means(tip.length = 0, label = "p.signif", method = "wilcox.test") +
        theme_Publication() +
        theme(axis.text.x = element_text(hjust = 1, angle = 45)))

(shaneth3 <- ggplot(data = dftot %>% filter(!is.na(HT_new)), aes(x = HT_new, y = shannon)) +
        geom_violin(aes(alpha = HT_new), fill = pal_bmj()(3)[3], show.legend = FALSE) +
        scale_alpha_manual(values = c(0.4, 1.0)) +
        geom_boxplot(fill = "white", width = 0.2) +
        labs(y = "Shannon index", x = "", alpha = "", fill = "") +
        facet_wrap(~EthnicityTot) +
        stat_compare_means(tip.length = 0, label = "p.signif", method = "wilcox.test") +
        theme_Publication() +
        theme(axis.text.x = element_text(hjust = 1, angle = 45)))


ggarrange(shaneth1, shaneth2, shaneth3, nrow = 3, ncol = 1, labels = LETTERS[1:3])
ggsave("results/alphadiversity/shannon_ethnicity_newdiagnoses.pdf", width = 8, height = 24)