## Associations microbes and new diagnoses

library(rio)
library(phyloseq)
library(tidyverse)
library(ggsci)
library(broom)
library(patchwork)
library(compositions)  # for CLR transformation

theme_Publication <- function(base_size=12, base_family="sans") {
    library(grid)
    library(ggthemes)
    (theme_foundation(base_size=base_size, base_family=base_family)
        + theme(plot.title = element_text(face = "bold",
                                          size = rel(0.8), hjust = 0.5),
                #family = 'Helvetica'
                text = element_text(),
                panel.background = element_rect(colour = NA, fill = NA),
                plot.background = element_rect(colour = NA, fill = NA),
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
            #   axis.line.y = element_blank(), axis.ticks.y = element_blank(),
              axis.line.x = element_blank(), axis.ticks.x = element_blank(),
              panel.grid = element_blank()))
}

#### Data ####
df <- readRDS("data/clinicaldata_wide.RDS")
mb <- readRDS("data/phyloseq_rarefied_cleaned.RDS")
tax <- readRDS("data/taxtable_rarefied_cleaned.RDS")

#### Preprocessing microbiome data ####
mb <- as.data.frame(t(as(mb@otu_table, "matrix")))

# CLR (Centered Log-Ratio) transformation BEFORE filtering
# This preserves the compositional structure of the full dataset
mb_clr <- as.data.frame(clr(mb + 0.5))  # Add pseudocount of 0.5 to handle zeros
# mb_clr <- as.data.frame(log10(mb + 0.5))  # log10 transformation with pseudocount of 1

# Filter ASVs AFTER CLR transformation
tk <- apply(mb, 2, function(x) sum(x > 5) > (0.2*length(x)))
mb_clr <- mb_clr[,tk]  # Apply same filtering to CLR-transformed data

head(mb_clr)
mb <- mb_clr
mb$sampleID_baseline <- rownames(mb)
mbclin <- left_join(mb, df, by = "sampleID_baseline")
dim(mbclin)
head(mbclin)[1:5,1:5]

#### Associations with new diagnoses ####
names(mbclin)[1:ncol(mb)-1]
mbclin <- mbclin %>% mutate(
    across(c("DM_new", "HT_new", "Dyslip_new"), ~fct_recode(.x, "1" = "Yes", "0" = "No"))
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
    dys <- glm(Dyslip_new ~ asv + Age_baseline, 
                  data = mbclin %>% filter(Dyslipidemia_baseline == "No"), family = "binomial")
    # extract estimates for asv
    dm <- tidy(dm, conf.int=TRUE, exponentiate = TRUE)[2,]
    ht <- tidy(ht, conf.int=TRUE, exponentiate = TRUE)[2,]
    dys <- tidy(dys, conf.int=TRUE, exponentiate = TRUE)[2,]
    # define rows
    row1 <- c(asvname, "diabetes", "model 1", dm$estimate, dm$conf.low, dm$conf.high, dm$p.value)
    row2 <- c(asvname, "hypertension", "model 1", ht$estimate, ht$conf.low, ht$conf.high, ht$p.value)
    row3 <- c(asvname, "dyslipidemia","model 1", dys$estimate, dys$conf.low, dys$conf.high, dys$p.value)
    names(row1) <- c("ASV", "diagnosis", "model", "OR", "lower", "upper", "pvalue")
    names(row2) <- c("ASV", "diagnosis", "model","OR", "lower", "upper", "pvalue")
    names(row3) <- c("ASV", "diagnosis","model", "OR", "lower", "upper", "pvalue")
    # bind
    res <- rbind(res, row1, row2, row3)

    # adjusted models
    dm2 <- glm(DM_new ~ asv + Age_baseline + Sex +CurrentSmoking_baseline + AlcBin_baseline + BMI_baseline, 
              data = mbclin %>% filter(DM_baseline == "No"), family = "binomial")
    ht2 <- glm(HT_new ~ asv  + Age_baseline + Sex + CurrentSmoking_baseline + AlcBin_baseline + BMI_baseline, 
              data = mbclin %>% filter(HT_BPMed_baseline == "No"), family = "binomial")
    dys2 <- glm(Dyslip_new ~ asv + Age_baseline + Sex + CurrentSmoking_baseline + AlcBin_baseline + BMI_baseline, 
                  data = mbclin %>% filter(Dyslipidemia_baseline == "No"), family = "binomial")
    # extract estimates for asv
    dm2 <- tidy(dm2, conf.int=TRUE, exponentiate = TRUE)[2,]
    ht2 <- tidy(ht2, conf.int=TRUE, exponentiate = TRUE)[2,]
    dys2 <- tidy(dys2, conf.int=TRUE, exponentiate = TRUE)[2,]
    # define rows
    row1 <- c(asvname, "diabetes", "model 2", dm2$estimate, dm2$conf.low, dm2$conf.high, dm2$p.value)
    row2 <- c(asvname, "hypertension","model 2",  ht2$estimate, ht2$conf.low, ht2$conf.high, ht2$p.value)
    row3 <- c(asvname, "dyslipidemia", "model 2", dys2$estimate, dys2$conf.low, dys2$conf.high, dys2$p.value)
    names(row1) <- c("ASV", "diagnosis", "model", "OR", "lower", "upper", "pvalue")
    names(row2) <- c("ASV", "diagnosis", "model","OR", "lower", "upper", "pvalue")
    names(row3) <- c("ASV", "diagnosis","model", "OR", "lower", "upper", "pvalue")
    # bind
    res <- rbind(res, row1, row2, row3)
}
res <- as.data.frame(res)
res2 <- left_join(res, tax, by = 'ASV') %>% mutate(across(c(4:7), as.numeric),
                                           model = as.factor(model),
                                           diagnosis = as.factor(diagnosis),
                                           model = fct_recode(model, 
                                                "Adjusted for age, sex, BMI, smoking, alcohol" = "model 2",
                                                "Adjusted for age" = "model 1")
                                           ) 

asvdm <- res2 |> filter(diagnosis == "diabetes" & model == "Adjusted for age" & pvalue < 0.05) |> pull(ASV)
dm <- res2 %>% filter(diagnosis == "diabetes" & ASV %in% asvdm) %>%
    group_by(model) |> 
    mutate(padj = p.adjust(pvalue, "fdr"))
asvdm <- dm |> filter(model == "Adjusted for age" & padj < 0.05) |> pull(ASV)
dm <- dm |> filter(ASV %in% asvdm) |> # filter for FDR < 0.05                  
        mutate(Tax = as.factor(make.unique(Tax)),
            Tax = fct_reorder(Tax, ifelse(model == "Adjusted for age", OR, NA), .na_rm = TRUE, .desc = TRUE),
            psig = case_when(padj < 0.05 ~ "yes", TRUE ~ "no")
           )

asvht <- res2 |> filter(diagnosis == "hypertension" & model == "Adjusted for age" & pvalue < 0.05) |> pull(ASV)
ht <- res2 %>% filter(diagnosis == "hypertension" & ASV %in% asvht) %>%
    group_by(model) |> 
    mutate(padj = p.adjust(pvalue, "fdr"))
asvht <- ht |> filter(model == "Adjusted for age" & padj < 0.05) |> pull(ASV)
ht <- ht |> filter(ASV %in% asvht) |> # filter for FDR < 0.05
    mutate(Tax = as.factor(make.unique(Tax)),
           Tax = fct_reorder(Tax, ifelse(model == "Adjusted for age", OR, NA), .na_rm = TRUE, .desc = TRUE),
           psig = case_when(padj < 0.05 ~ "yes", TRUE ~ "no")
           )
head(ht)

asvdys <- res2 |> filter(diagnosis == "dyslipidemia" & model == "Adjusted for age" & pvalue < 0.05) |> pull(ASV)
dyslipidemia <- res2 %>% filter(diagnosis == "dyslipidemia"& ASV %in% asvht) %>%
     group_by(model) |> 
    mutate(padj = p.adjust(pvalue, "fdr")) # calc FDR adj p value
asvdys <- dyslipidemia |> filter(model == "Adjusted for age" & padj < 0.05) |> pull(ASV)
dyslipidemia <- dyslipidemia |> filter(ASV %in% asvdys) |> # filter for FDR < 0.05
    mutate(Tax = as.factor(make.unique(Tax)),
           Tax = fct_reorder(Tax, ifelse(model == "Adjusted for age", OR, NA), .na_rm = TRUE, .desc = TRUE),
           psig = case_when(padj < 0.05 ~ "yes", TRUE ~ "no")
           )

dim(dm)
dim(ht)
dim(dyslipidemia)

summary(as.factor(dm$Family))
summary(as.factor(dyslipidemia$Family))
summary(as.factor(ht$Family))
dm %>% filter(OR > 1)
ht %>% filter(OR > 1)
ht %>% filter(OR < 1)
dyslipidemia %>% filter(OR > 1)
summary(dyslipidemia$ASV %in% dm$ASV)
summary(dm$ASV %in% dyslipidemia$ASV)

coltab <- NULL
pal <- colorRampPalette(ggsci::pal_bmj()(9))
families <- c("Akkermansiaceae", "Anaerovoracaceae", "Bacteroidaceae", "Bifidobacteriaceae",
              "Christensenellaceae", "Clostridiaceae", "Desulfovibrionaceae", "Erysipelatoclostridiaceae",
              "Lachnospiraceae", "Marinifilaceae", "Monoglobaceae", "Oscillospiraceae", "Peptostreptococcaceae",
              "Rikenellaceae", "Ruminococcaceae", "Tannerellaceae", "Veillonellaceae"
              ) # to make colors in line with the registry analyses
coltab <- setNames(pal(length(families)), families)
families2 <- c("UCG-010", "Eggerthellaceae", "[Eubacterium] coprostanoligenes group", "Butyricicoccaceae",
               "Sutterellaceae", "Streptococcaceae", "Erysipelotrichaceae", "Pasteurellaceae", "Unknown") # not occuring in registry analyses
# check if Pasteurellaceae is not in registry analyses
famtot <- c(families, families2)
pal2 <- colorRampPalette(ggsci::pal_futurama()(8)[c(1,3,4,7,8)])
coltab2 <- setNames(c(pal2(length(families2)-3), "yellow1", "purple", "grey80"), families2)
coltab <- coltab[which((names(coltab) %in% c(dm$Family, ht$Family, dyslipidemia$Family)))]
coltab <- c(coltab, coltab2)

dm <- dm %>% mutate(Family = case_when(is.na(Family) ~ "Unknown", .default = Family),
                    Family = factor(Family, levels = famtot))
ht <- ht %>% mutate(Family = case_when(is.na(Family) ~ "Unknown", .default = Family),
                    Family = factor(Family, levels = famtot))
dyslipidemia <- dyslipidemia %>% mutate(Family = case_when(is.na(Family) ~ "Unknown", .default = Family),
                            Family = factor(Family, levels = famtot))

write.csv(dm, "results/diabetesassociations_smoking_alcohol_bmi.csv")
write.csv(dyslipidemia, "results/dyslipidemiaassociations_smoking_alcohol_bmi.csv")
write.csv(ht, "results/hypertensionassociations_smoking_alcohol_bmi.csv")

(bar <- ggplot(dm, aes(y = Tax, x = 0, fill = Family)) +
        scale_fill_manual(values = coltab, drop = FALSE) +
        geom_tile(show.legend = TRUE) +
        labs(x = "", y = "", fill = "") +
        guides(fill = guide_legend(ncol = 1)) +
        theme_bar())

(pldm <- ggplot(dm, aes(x = OR, y = Tax)) +
            geom_vline(aes(xintercept = 1), color = "darkgrey") +
            geom_errorbarh(aes(xmin = lower, xmax = upper, color = model), width = 0.5, 
                position = position_dodge(-0.8)) +
            geom_point(aes(color = model, shape = psig), position = position_dodge(-0.8), size = 2.5) +
            scale_color_manual(values = c(pal_bmj()(7)[7], "grey60")) +
            scale_shape_manual(values = c(21, 19)) +
            scale_x_continuous(n.breaks = 6, limits = c(0.50, 1.50)) +
            labs(title = "Diabetes", y = "", x = "OR per CLR-increase in ASV", shape = "padj < 0.05") +
            theme_Publication() +
            theme(axis.text.y = element_blank()))

(dmplot <- bar + plot_spacer() + pldm + 
    plot_layout(widths = c(0.1, -0.075, 0.7), axis_titles = "collect_y"))

(bar <- ggplot(dyslipidemia, aes(y = Tax, x = 0, fill = Family)) +
    scale_fill_manual(values = coltab, drop = FALSE) +
    geom_tile(show.legend = TRUE) +
    labs(x = "", y = "", fill = "") +
    guides(fill = guide_legend(ncol = 1)) +
    theme_bar())

(pldl <- ggplot(dyslipidemia, aes(x = OR, y = Tax)) +
        geom_vline(aes(xintercept = 1), color = "darkgrey") +
        geom_errorbarh(aes(xmin = lower, xmax = upper, color = model), width = 0.5, position = position_dodge(-0.8)) +
        geom_point(aes(color = model, shape = psig), position = position_dodge(-0.8), size = 2.5) +
        scale_color_manual(values = c(pal_bmj()(4)[4], "grey60")) +
        scale_shape_manual(values = c(21, 19)) +
        scale_x_continuous(n.breaks = 6, limits = c(0.50, 1.50)) +
        labs(title = "Dyslipidemia", y = "", x = "OR per CLR-increase in ASV", shape = "padj < 0.05") +
        theme_Publication() +
        theme(axis.text.y = element_blank()))

(dlplot <- bar + plot_spacer() + pldl + 
    plot_layout(widths = c(0.1, -0.075, 0.7), axis_titles = "collect_y"))
  
(bar <- ggplot(ht, aes(y = Tax, x = 0, fill = Family)) +
        scale_fill_manual(values = coltab, drop = FALSE) +
        geom_tile(show.legend = TRUE) +
        labs(x = "", y = "", fill = "") +
        guides(fill = guide_legend(ncol = 1)) +
        theme_bar())

(plht <- ggplot(ht, aes(x = OR, y = Tax)) +
        geom_vline(aes(xintercept = 1), color = "darkgrey") +
        geom_errorbarh(aes(xmin = lower, xmax = upper, color = model), width = 0.5, position = position_dodge(-0.8)) +
        geom_point(aes(color = model, shape = psig), position = position_dodge(-0.8), size = 2.5) +
        scale_color_manual(values = c(pal_bmj()(3)[3], "grey60")) +
        scale_shape_manual(values = c(21, 19)) +
        scale_x_continuous(n.breaks = 6, limits = c(0.50, 1.50)) +
        labs(title = "Hypertension", y = "", x = "OR per CLR-increase in ASV", shape = "padj < 0.05") +
        theme_Publication()+
        theme(axis.text.y = element_blank()))

(htplot <- bar + plot_spacer() + plht + 
    plot_layout(widths = c(0.1, -0.075, 0.7), axis_titles = "collect_y"))

dmplot / dlplot / htplot +
    plot_layout(guides = "collect", nrow = 3, heights = c(0.75, 0.25, 1.1)) +
    plot_annotation(tag_levels = list(c("A", "", "B", "", "C", ""))) &
    theme(plot.tag = element_text(face = "bold"),
          legend.key.size= unit(0.4, "cm"),
          legend.text = element_text(size = rel(1.0)))

ggsave("results/microbiome_diagnoses/diagnosispred_long_adjusted.pdf", width = 14, height = 40)
