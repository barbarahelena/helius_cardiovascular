# Table 1

# libraries
library(tableone)
library(dplyr)
library(ggplot2)
library(ggsci)
library(stringr)

theme_Publication <- function(baselinese_size=14, baselinese_family="sans") {
    library(grid)
    library(ggthemes)
    library(stringr)
    (theme_foundation(baselinese_size=baselinese_size, baselinese_family=baselinese_family)
        + theme(plot.title = element_text(face = "bold",
                                          size = rel(1.0), hjust = 0.5),
                text = element_text(),
                panel.baselineckground = element_rect(colour = NA),
                plot.baselineckground = element_rect(colour = NA),
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
                strip.baselineckground=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                strip.text = element_text(face="bold"),
                plot.caption = element_text(size = rel(0.5), face = "italic")
        ))
    
} 


##### Table 1 #####
helius <- readRDS('data/clinicaldata_wide.RDS')

table1 <- helius %>%
    dplyr::select(Age_baseline, Sex, EthnicityTot, MigrGen, FUtime, BMI_baseline, Smoking_baseline, AlcCons_baseline, 
           DM_baseline, SBP_baseline, DBP_baseline, HT_BPMed_baseline, MetSyn_baseline, 
           PPI_baseline, Metformin_baseline, Statins_baseline, TC_baseline, LDL_baseline, Trig_baseline, 
           HbA1c_baseline, SCORECVDmortNL_baseline) %>% 
    CreateTableOne(data=., test = FALSE) %>% 
    print(nonnormal=c("SCORECVDmortNL_baseline", "Trig"), noSpaces = TRUE, pDigits = 3, contDigits = 1) %>% 
    as.data.frame(.)

table2 <- helius %>% 
    dplyr::select(Age_baseline, Sex, EthnicityTot, FUtime, 
                  BMI_baseline, `BMI_follow-up`,
                  Smoking_baseline, `Smoking_follow-up`, AlcCons_baseline, `AlcCons_follow-up`,
                  DM_baseline, `DM_follow-up`, HbA1c_baseline, `HbA1c_follow-up`, 
                  SBP_baseline, `SBP_follow-up`, DBP_baseline, `DBP_follow-up`,
                  HT_BPMed_baseline, `HT_BPMed_follow-up`,
                  MetSyn_baseline, `MetSyn_follow-up`,
                  PPI_baseline, Metformin_baseline, Statins_baseline, TC_baseline, LDL_baseline, 
                  Trig_baseline, SCORECVDmortNL_baseline) %>% 
    CreateTableOne(data=., strata = 'EthnicityTot', test = FALSE, addOverall = TRUE) %>% 
    print(nonnormal=c("SCORECVDmortNL_baseline", "Trig"), noSpaces = TRUE, pDigits = 3, contDigits = 1) %>% 
    as.data.frame(.)

table(helius$DM_new, helius$EthnicityTot)
table(helius$HT_new, helius$EthnicityTot)
table(helius$MetSyn_new, helius$EthnicityTot)

helius %>% filter(!is.na(DM_new)) %>% group_by(EthnicityTot) %>% 
    summarise(rate = length(DM_new[which(DM_new == "Yes")]) / length(DM_new[which(DM_new == "No")]) * 100,
              count = length(DM_new[which(DM_new == "Yes")]))
helius %>% filter(!is.na(MetSyn_new)) %>% group_by(EthnicityTot) %>% 
    summarise(rate = length(MetSyn_new[which(MetSyn_new == "Yes")]) / length(MetSyn_new[which(MetSyn_new == "No")]) * 100,
              count = length(MetSyn_new[which(HT_new == "Yes")]))
helius %>% filter(!is.na(HT_new)) %>% group_by(EthnicityTot) %>% 
    summarise(rate = length(HT_new[which(HT_new == "Yes")]) / length(HT_new[which(HT_new == "No")]) * 100,
              count = length(HT_new[which(HT_new == "Yes")]))

helius %>% filter(!is.na(DM_new)) %>% 
    summarise(rate = length(DM_new[which(DM_new == "Yes")]) / 
                  (length(DM_new[which(DM_new == "No")])+length(DM_new[which(DM_new == "Yes")])) * 100,
              count = length(HT_new[which(DM_new == "Yes")]),
              total = nrow(.))
helius %>% filter(!is.na(MetSyn_new)) %>% 
    summarise(rate = length(MetSyn_new[which(MetSyn_new == "Yes")]) / 
                  (length(MetSyn_new[which(MetSyn_new == "No")])+length(DM_new[which(DM_new == "Yes")])) * 100,
              count = length(MetSyn_new[which(HT_new == "Yes")]),
              total = nrow(.))
helius %>% filter(!is.na(HT_new)) %>% 
    summarise(rate = length(HT_new[which(HT_new == "Yes")]) / 
                  (length(HT_new[which(HT_new == "No")])+length(DM_new[which(DM_new == "Yes")])) * 100,
              count = length(HT_new[which(HT_new == "Yes")]),
              total = nrow(.))

write.csv2(as.data.frame(table1), 'results/tables/table1.csv')
write.csv2(as.data.frame(table2), 'results/tables/table2.csv')

## Table 2
helius <- readRDS('data/clinicaldata_long.RDS')
table1 <- helius %>% 
    dplyr::select(Age, Sex, EthnicityTot, MigrGen, FUtime, BMI, Smoking, AlcCons, 
           DM, SBP, DBP, HT_BPMed, MetSyn, 
           PPI, Metformin, Statins, TC, LDL, Trig, HbA1c, SCORECVDmortNL,
           timepoint) %>% 
    CreateTableOne(data=., strata = 'timepoint', test = FALSE) %>% 
    print(nonnormal=c("SCORECVDmortNL", "Trig"), noSpaces = TRUE, pDigits = 3, contDigits = 1) %>% 
    as.data.frame(.)

write.csv2(as.data.frame(table1), 'results/tables/followuptable.csv')
