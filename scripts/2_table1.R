# Table 1

# libraries
library(tableone)
library(dplyr)

##### Supplementary Table 3B #####
helius <- readRDS('data/clinicaldata_wide.RDS')

table2 <- helius %>% 
    dplyr::select(Age_baseline, Sex, EthnicityTot, FUtime, 
                  BMI_baseline, `BMI_follow-up`,
                  Smoking_baseline, `Smoking_follow-up`, AlcCons_baseline, `AlcCons_follow-up`,
                  DM_baseline, `DM_follow-up`, HbA1c_baseline, `HbA1c_follow-up`, 
                  SBP_baseline, `SBP_follow-up`, DBP_baseline, `DBP_follow-up`,
                  HT_BPMed_baseline, `HT_BPMed_follow-up`,
                  Dyslipidemia_baseline, `Dyslipidemia_follow-up`,
                  PPI_baseline, Metformin_baseline, Statins_baseline, TC_baseline, LDL_baseline, 
                  Trig_baseline, SCORECVDmortNL_baseline) %>% 
    CreateTableOne(data=., strata = 'EthnicityTot', test = FALSE, addOverall = TRUE) %>% 
    print(nonnormal=c("SCORECVDmortNL_baseline", "Trig_baseline"), noSpaces = TRUE, 
          pDigits = 3, contDigits = 1, missing = TRUE) %>% 
    as.data.frame(.)

table(helius$DM_new, helius$EthnicityTot)
table(helius$HT_new, helius$EthnicityTot)
table(helius$Dyslip_new, helius$EthnicityTot)

helius %>% filter(!is.na(DM_new)) %>% group_by(EthnicityTot) %>% 
    summarise(rate = length(DM_new[which(DM_new == "Yes")]) / length(DM_new[which(DM_new == "No")]) * 100,
              count = length(DM_new[which(DM_new == "Yes")]))
helius %>% filter(!is.na(Dyslip_new)) %>% group_by(EthnicityTot) %>% 
    summarise(rate = length(Dyslip_new[which(Dyslip_new == "Yes")]) / length(Dyslip_new[which(Dyslip_new == "No")]) * 100,
              count = length(Dyslip_new[which(Dyslip_new == "Yes")]))
helius %>% filter(!is.na(HT_new)) %>% group_by(EthnicityTot) %>% 
    summarise(rate = length(HT_new[which(HT_new == "Yes")]) / length(HT_new[which(HT_new == "No")]) * 100,
              count = length(HT_new[which(HT_new == "Yes")]))

helius %>% filter(!is.na(DM_new)) %>% 
    summarise(rate = length(DM_new[which(DM_new == "Yes")]) / 
                  (length(DM_new[which(DM_new == "No")])+length(DM_new[which(DM_new == "Yes")])) * 100,
              count = length(HT_new[which(DM_new == "Yes")]),
              total = nrow(.))
helius %>% filter(!is.na(Dyslip_new)) %>% 
    summarise(rate = length(Dyslip_new[which(Dyslip_new == "Yes")]) / 
                  (length(Dyslip_new[which(Dyslip_new == "No")])+length(DM_new[which(DM_new == "Yes")])) * 100,
              count = length(Dyslip_new[which(Dyslip_new == "Yes")]),
              total = nrow(.))
helius %>% filter(!is.na(HT_new)) %>% 
    summarise(rate = length(HT_new[which(HT_new == "Yes")]) / 
                  (length(HT_new[which(HT_new == "No")])+length(DM_new[which(DM_new == "Yes")])) * 100,
              count = length(HT_new[which(HT_new == "Yes")]),
              total = nrow(.))

write.csv2(as.data.frame(table1), 'results/tables/table1.csv')
write.csv2(as.data.frame(table2), 'results/tables/table2.csv')

## Table 1 ##
helius <- readRDS('data/clinicaldata_long.RDS')
table1 <- helius %>% 
    dplyr::select(Age, Sex, EthnicityTot, MigrGen, FUtime, BMI, Smoking, AlcCons, 
           DM, SBP, DBP, HT_BPMed, Dyslipidemia, MetSyn,
           PPI, Metformin, AB, Statins, TC, LDL, Trig, HbA1c, SCORECVDmortNL,
           timepoint) %>% 
    CreateTableOne(data=., strata = c('timepoint', 'Sex'), test = TRUE, addOverall = TRUE) %>% 
    print(nonnormal=c("SCORECVDmortNL", "Trig"), noSpaces = TRUE, 
          missing = TRUE, pDigits = 3, contDigits = 1) %>% 
    as.data.frame(.)

write.csv2(as.data.frame(table1), 'results/tables/followuptable.csv')

## Supplementary Table 3 ##
met <- readRDS("data/metabolomics/metabolomics_paired.RDS")
metsel <- rownames(met[str_detect(rownames(met), "HELIBA"),])
metsel 
helius_met <- helius |> filter(sampleID %in% metsel) |> droplevels()
table3 <- helius_met %>% 
    dplyr::select(Age, Sex, EthnicityTot, MigrGen, FUtime, BMI, Smoking, AlcCons, 
           DM, SBP, DBP, HT_BPMed, Dyslipidemia, MetSyn,
           PPI, Metformin, AB, Statins, TC, LDL, Trig, HbA1c, SCORECVDmortNL,
           timepoint) %>% 
    CreateTableOne(data=., strata = 'EthnicityTot', test = TRUE, addOverall = TRUE) %>% 
    print(nonnormal=c("SCORECVDmortNL", "Trig"), noSpaces = TRUE, 
          missing = TRUE, pDigits = 3, contDigits = 1) %>% 
    as.data.frame(.)
write.csv2(as.data.frame(table3), 'results/tables/metabolomics.csv')
