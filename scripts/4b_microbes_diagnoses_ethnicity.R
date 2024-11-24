## Associations microbes and new diagnoses

library(rio)
library(phyloseq)
library(tidyverse)
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

#### Data ####
df <- readRDS("data/clinicaldata_wide.RDS") %>% filter(SampleAB_baseline == "No")
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

#### Diabetes associations ####
dmpred <- read.csv("results/diabetesassociations.csv")
mbdm <- mb %>% dplyr::select(sampleID_baseline, dmpred$ASV)
mbdmclin <- left_join(mbdm, df, by = "sampleID_baseline")
dim(mbdmclin)
head(mbdmclin)[1:5,1:5]
names(mbdmclin)[1:ncol(mbdmclin)]
mbdmclin <- mbdmclin %>% mutate(
    across(c("DM_new", "HT_new", "MetSyn_new"), ~fct_recode(.x, "1" = "Yes", "0" = "No"))
)
resdm <- c()
length(names(mbdmclin)[2:ncol(mbdm)])
dim(dmpred)
ethsel <- c("Dutch", "South-Asian Surinamese", "African Surinamese")
for(a in 2:(ncol(mbdm))) {
    mbdmclin$asv <- mbdmclin[[a]]
    asvname <- colnames(mbdmclin)[a]
    print(asvname)
    # run models for each diagnosis while excluding participants with baseline diagnoses
    dm_dutch <- glm(DM_new ~ asv + Age_baseline, data = mbdmclin %>% 
                        filter(DM_baseline == "No" & EthnicityTot == "Dutch"), family = "binomial")
    dm_sas <- glm(DM_new ~ asv + Age_baseline, data = mbdmclin %>% 
                        filter(DM_baseline == "No" & EthnicityTot == "South-Asian Surinamese"), family = "binomial")
    dm_as <- glm(DM_new ~ asv + Age_baseline, data = mbdmclin %>% 
                      filter(DM_baseline == "No" & EthnicityTot == "African Surinamese"), family = "binomial")
    dm_ia <- glm(DM_new ~ asv + Age_baseline + asv*EthnicityTot, data = mbdmclin %>% 
                     filter(DM_baseline == "No" & EthnicityTot %in% ethsel), family = "binomial")

    # extract estimates for variable sex
    dm_dutch <- tidy(dm_dutch, conf.int=TRUE, exponentiate = TRUE)[2,]
    dm_sas <- tidy(dm_sas, conf.int=TRUE, exponentiate = TRUE)[2,]
    dm_as <- tidy(dm_as, conf.int=TRUE, exponentiate = TRUE)[2,]
    dm_int <- tidy(dm_ia)[6:7,]
    # define rows
    row1 <- c(asvname, "Dutch", dm_dutch$estimate, dm_dutch$conf.low, dm_dutch$conf.high, dm_dutch$p.value, "")
    row2 <- c(asvname, "South-Asian Surinamese", dm_sas$estimate, dm_sas$conf.low, dm_sas$conf.high, dm_sas$p.value, dm_int$p.value[1])
    row3 <- c(asvname, "African Surinamese", dm_as$estimate, dm_as$conf.low, dm_as$conf.high, dm_as$p.value, dm_int$p.value[2])
    names(row1) <- c("ASV", "Ethnicity", "OR", "lower", "upper", "pvalue", "interaction")
    names(row2) <- c("ASV", "Ethnicity", "OR", "lower", "upper", "pvalue", "interaction")
    names(row3) <- c("ASV", "Ethnicity", "OR", "lower", "upper", "pvalue", "interaction")
    # bind
    resdm <- rbind(resdm, row1, row2, row3)
}
resdm <- as.data.frame(resdm)

dm <- resdm %>%
    mutate(across(c(3:7), as.numeric)) %>% 
    pivot_wider(., id_cols = "ASV", 
                names_from = "Ethnicity", values_from = c("OR", "lower", "upper", 
                                                          "pvalue", "interaction")) %>% 
    left_join(., tax, by = "ASV") %>% 
    mutate(Tax = as.factor(make.unique(Tax)),
           Tax = forcats::fct_rev(forcats::fct_inorder(Tax)),
           interactionsig = case_when(`interaction_South-Asian Surinamese` < 0.05 |
                                       `interaction_African Surinamese` < 0.05 ~ "sig", 
                                      .default = "not sig")
    ) %>% 
    pivot_longer(., cols = 2:16, names_to = c("var", "Ethnicity"), 
                 names_sep = "_", values_to = "val") %>% 
    pivot_wider(., names_from = "var", values_from = "val") 

dm %>% filter(interactionsig == "sig")
interactions <- dm %>% filter(Ethnicity == "Dutch")

(pldm <- ggplot(dm, aes(x = OR, y = Tax, color = Ethnicity)) +
        geom_vline(aes(xintercept = 1), color = "darkgrey") +
        scale_color_bmj() +
        geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.5, position = position_dodge(0.5)) +
        geom_point(position = position_dodge(0.5)) +
        scale_x_continuous(n.breaks = 6, limits = c(0.0, 4.5)) +
        labs(title = "Diabetes", y = "", x = "OR per log10-increase ASV", color = "") +
        theme_Publication() +
        theme(axis.text.y = element_text(face = ifelse(rev(interactions$interactionsig) == "sig", "bold", "plain"))))
ggsave(pldm, filename = "results/diabetes_interactions.pdf", width = 10, height = 12)

#### MetSyn ####
mspred <- read.csv("results/metsynassociationcs.csv")
mbms <- mb %>% dplyr::select(sampleID_baseline, mspred$ASV)
mbmsclin <- left_join(mbms, df, by = "sampleID_baseline")
head(mbmsclin)[1:5,1:5]
mbmsclin <- mbmsclin %>% mutate(across(c("DM_new", "HT_new", "MetSyn_new"), ~fct_recode(.x, "1" = "Yes", "0" = "No")))
resms <- c()
length(names(mbmsclin)[2:ncol(mbms)])
dim(mspred)
ethsel <- c("Dutch", "South-Asian Surinamese", "African Surinamese")
for(a in 2:(ncol(mbms))) {
    mbmsclin$asv <- mbmsclin[[a]]
    asvname <- colnames(mbmsclin)[a]
    print(asvname)
    # run models for each diagnosis while excluding participants with baseline diagnoses
    ms_dutch <- glm(MetSyn_new ~ asv + Age_baseline, data = mbmsclin %>% 
                        filter(MetSyn_baseline == "No" & EthnicityTot == "Dutch"), family = "binomial")
    ms_sas <- glm(MetSyn_new ~ asv + Age_baseline, data = mbmsclin %>% 
                      filter(MetSyn_baseline == "No" & EthnicityTot == "South-Asian Surinamese"), family = "binomial")
    ms_as <- glm(MetSyn_new ~ asv + Age_baseline, data = mbmsclin %>% 
                     filter(MetSyn_baseline == "No" & EthnicityTot == "African Surinamese"), family = "binomial")
    ms_int <- glm(MetSyn_new ~ asv + Age_baseline + asv*EthnicityTot, data = mbmsclin %>% 
                     filter(MetSyn_baseline == "No" & EthnicityTot %in% ethsel), family = "binomial")
    
    # extract estimates for variable sex
    ms_dutch <- tidy(ms_dutch, conf.int=TRUE, exponentiate = TRUE)[2,]
    ms_sas <- tidy(ms_sas, conf.int=TRUE, exponentiate = TRUE)[2,]
    ms_as <- tidy(ms_as, conf.int=TRUE, exponentiate = TRUE)[2,]
    ms_int <- tidy(ms_int)[6:7,]
    # define rows
    row1 <- c(asvname, "Dutch", ms_dutch$estimate, ms_dutch$conf.low, ms_dutch$conf.high, ms_dutch$p.value, "")
    row2 <- c(asvname, "South-Asian Surinamese", ms_sas$estimate, ms_sas$conf.low, ms_sas$conf.high, ms_sas$p.value, ms_int$p.value[1])
    row3 <- c(asvname, "African Surinamese", ms_as$estimate, ms_as$conf.low, ms_as$conf.high, ms_as$p.value, ms_int$p.value[2])
    names(row1) <- c("ASV", "Ethnicity", "OR", "lower", "upper", "pvalue", "interaction")
    names(row2) <- c("ASV", "Ethnicity", "OR", "lower", "upper", "pvalue", "interaction")
    names(row3) <- c("ASV", "Ethnicity", "OR", "lower", "upper", "pvalue", "interaction")
    # bind
    resms <- rbind(resms, row1, row2, row3)
}
resms <- as.data.frame(resms)

ms <- resms %>%
    mutate(across(c(3:7), as.numeric)) %>% 
    pivot_wider(., id_cols = "ASV", 
                names_from = "Ethnicity", values_from = c("OR", "lower", "upper", 
                                                          "pvalue", "interaction")) %>% 
    left_join(., tax, by = "ASV") %>% 
    mutate(Tax = as.factor(make.unique(Tax)),
           Tax = forcats::fct_rev(forcats::fct_inorder(Tax)),
           interactionsig = case_when(`interaction_South-Asian Surinamese` < 0.05 |
                                          `interaction_African Surinamese` < 0.05 ~ "sig", 
                                      .default = "not sig")
    ) %>% 
    pivot_longer(., cols = 2:16, names_to = c("var", "Ethnicity"), 
                 names_sep = "_", values_to = "val") %>% 
    pivot_wider(., names_from = "var", values_from = "val") 

ms %>% filter(interactionsig == "sig")
interactions <- ms %>% filter(Ethnicity == "Dutch")

(plms <- ggplot(ms, aes(x = OR, y = Tax, color = Ethnicity)) +
        geom_vline(aes(xintercept = 1), color = "darkgrey") +
        scale_color_bmj() +
        geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.5, position = position_dodge(0.5)) +
        geom_point(position = position_dodge(0.5)) +
        scale_x_continuous(n.breaks = 6, limits = c(0.0, 4.5)) +
        labs(title = "MetSyn", y = "", x = "OR per log10-increase ASV", color = "") +
        theme_Publication() +
        theme(axis.text.y = element_text(face = ifelse(rev(interactions$interactionsig) == "sig", "bold", "plain"))))
ggsave(plms, filename = "results/metsyn_interactions.pdf", width = 10, height = 20)

#### Hypertension ####
htpred <- read.csv("results/hypertensionassociations.csv")
mbht <- mb %>% dplyr::select(sampleID_baseline, htpred$ASV)
mbhtclin <- left_join(mbht, df, by = "sampleID_baseline")
head(mbhtclin)[1:5,1:5]
mbhtclin <- mbhtclin %>% mutate(across(c("DM_new", "HT_new", "MetSyn_new"), ~fct_recode(.x, "1" = "Yes", "0" = "No")))
resht <- c()
length(names(mbhtclin)[2:ncol(mbht)])
dim(htpred)
ethsel <- c("Dutch", "South-Asian Surinamese", "African Surinamese")
for(a in 2:(ncol(mbht))) {
    mbhtclin$asv <- mbhtclin[[a]]
    asvname <- colnames(mbhtclin)[a]
    print(asvname)
    # run models for each diagnosis while excluding participants with baseline diagnoses
    ht_dutch <- glm(HT_new ~ asv + Age_baseline, data = mbhtclin %>% 
                        filter(HT_BPMed_baseline == "No" & EthnicityTot == "Dutch"), family = "binomial")
    ht_sas <- glm(HT_new ~ asv + Age_baseline, data = mbhtclin %>% 
                      filter(HT_BPMed_baseline == "No" & EthnicityTot == "South-Asian Surinamese"), family = "binomial")
    ht_as <- glm(HT_new ~ asv + Age_baseline, data = mbhtclin %>% 
                     filter(HT_BPMed_baseline == "No" & EthnicityTot == "African Surinamese"), family = "binomial")
    ht_int <- glm(HT_new ~ asv + Age_baseline + asv*EthnicityTot, data = mbhtclin %>% 
                      filter(HT_BPMed_baseline == "No" & EthnicityTot %in% ethsel), family = "binomial")
    
    # extract estimates for variable sex
    ht_dutch <- tidy(ht_dutch, conf.int=TRUE, exponentiate = TRUE)[2,]
    ht_sas <- tidy(ht_sas, conf.int=TRUE, exponentiate = TRUE)[2,]
    ht_as <- tidy(ht_as, conf.int=TRUE, exponentiate = TRUE)[2,]
    ht_int <- tidy(ht_int)[6:7,]
    # define rows
    row1 <- c(asvname, "Dutch", ht_dutch$estimate, ht_dutch$conf.low, ht_dutch$conf.high, ht_dutch$p.value, "")
    row2 <- c(asvname, "South-Asian Surinamese", ht_sas$estimate, ht_sas$conf.low, ht_sas$conf.high, ht_sas$p.value, ht_int$p.value[1])
    row3 <- c(asvname, "African Surinamese", ht_as$estimate, ht_as$conf.low, ht_as$conf.high, ht_as$p.value, ht_int$p.value[2])
    names(row1) <- c("ASV", "Ethnicity", "OR", "lower", "upper", "pvalue", "interaction")
    names(row2) <- c("ASV", "Ethnicity", "OR", "lower", "upper", "pvalue", "interaction")
    names(row3) <- c("ASV", "Ethnicity", "OR", "lower", "upper", "pvalue", "interaction")
    # bind
    resht <- rbind(resht, row1, row2, row3)
}
resht <- as.data.frame(resht)

ht <- resht %>%
    mutate(across(c(3:7), as.numeric)) %>% 
    pivot_wider(., id_cols = "ASV", 
                names_from = "Ethnicity", values_from = c("OR", "lower", "upper", 
                                                          "pvalue", "interaction")) %>% 
    left_join(., tax, by = "ASV") %>% 
    mutate(Tax = as.factor(make.unique(Tax)),
           Tax = forcats::fct_rev(forcats::fct_inorder(Tax)),
           interactionsig = case_when(`interaction_South-Asian Surinamese` < 0.05 |
                                          `interaction_African Surinamese` < 0.05 ~ "sig", 
                                      .default = "not sig")
    ) %>% 
    pivot_longer(., cols = 2:16, names_to = c("var", "Ethnicity"), 
                 names_sep = "_", values_to = "val") %>% 
    pivot_wider(., names_from = "var", values_from = "val")  %>% 
    mutate(intadj = p.adjust(interaction, method = "fdr"))

ht %>% filter(interactionsig == "sig")
interactions <- ht %>% filter(Ethnicity == "Dutch")

(plht <- ggplot(ht, aes(x = OR, y = Tax, color = Ethnicity)) +
        geom_vline(aes(xintercept = 1), color = "darkgrey") +
        scale_color_bmj() +
        geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.5, position = position_dodge(0.5)) +
        geom_point(position = position_dodge(0.5)) +
        scale_x_continuous(n.breaks = 6, limits = c(0.0, 5)) +
        labs(title = "Hypertension", y = "", x = "OR per log10-increase ASV", color = "") +
        theme_Publication() +
        theme(axis.text.y = element_text(face = ifelse(rev(interactions$interactionsig) == "sig", "bold", "plain"))))
ggsave(plht, filename = "results/hypertension_interactions.pdf", width = 10, height = 28)
