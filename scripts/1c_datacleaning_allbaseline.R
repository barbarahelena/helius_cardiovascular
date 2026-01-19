# Data cleaning of clinical data, 16S, shotgun
## Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

## Libraries
library(tidyverse)
library(rio) 
library(haven)
library(phyloseq)
library(tableone)

#### Clinical data ####
## Open HELIUS clinical data
df <- haven::read_sav("data/240411_2_HELIUS data Barbara Verhaar -  Descriptives.sav")
names(df)
str(df$H1_FecesMicrobiome)
summary(as_factor(df$H1_FecesMicrobiome))

# Change type of variable 
yesnosmall <- function(x) fct_recode(x, "No"="nee", "Yes"="ja")
yesnocaps <- function(x) fct_recode(x, "No"="Nee", "Yes"="Ja")
zero_one <- function(x) fct_recode(x, "No"="0", "Yes"="1")

## Clean HELIUS dataframe
df_new <- df %>% 
    dplyr::select(# Demographics
                  ID, FUtime=H2_fu_time,
                  Age_BA=H1_lft, Age_FU = H2_lft,
                  Sex=H1_geslacht,  EthnicityTot=H1_EtnTotaal, 
                  Ethnicity = H1_etniciteit, MigrGen = H1_MigrGeneratie, 
                  CBS = project9639, Microbiome = H1H2_16S,
                  # Intox
                  Smoking_BA = H1_Roken, 
                  AlcCons_BA = H1_AlcoholConsumption,
                  # Physical exam
                  Physical_BA = H1_LichamelijkOnderzoekJN, Physical_FU=H2_LichamelijkOnderzoekJN,
                  BMI_BA = H1_LO_BMI, 
                  SBP_BA = H1_LO_GemBPSysZit,  
                  DBP_BA = H1_LO_GemBPDiaZit, 
                  # Hypertension
                  HT_Self_BA = H1_HT_Self, HT_Self_FU = H2_HT_Self_H1combined,
                  HT_BP_BA = H1_HT_BP, HT_BP_FU = H2_HT_BP,
                  HT_SelfBP_BA = H1_HT_SelfBP, HT_SelfBP_FU = H2_HT_Self_H1combined,
                  HT_SelfBPMed_BA = H1_HT_SelfBPMed, HT_SelfBPMed_FU = H2_HT_SelfBPMed_H1combined,
                  HT_BPMed_BA = H1_HT_BPMed, HT_BPMed_FU = H2_HT_BPMed,
                  # Metabolic diagnoses
                  MetSyn_BA = H1_MetSyn_MetabolicSyndrome, MetSyn_FU = H2_MetSyn_MetabolicSyndrome,
                  MetSyn_Obesity_BA = H1_MetSyn_CentralObesity, MetSyn_Obesity_FU = H2_MetSyn_CentralObesity, 
                  MetSyn_HighTG_BA = H1_MetSyn_HighTriglyceride, MetSyn_HighTG_FU = H2_MetSyn_HighTriglyceride, 
                  MetSyn_LowHDL_BA =H1_MetSyn_LowHDL, MetSyn_LowHDL_FU =H2_MetSyn_LowHDL,
                  MetSyn_HighBP_BA = H1_MetSyn_HighBP, MetSyn_HighBP_FU = H2_MetSyn_HighBP, 
                  MetSyn_HighGluc_BA = H1_MetSyn_HighGluc, MetSyn_HighGluc_FU = H2_MetSyn_HighGluc,
                  DM_BA = H1_Diabetes_GlucMed, DM_FU = H2_DM_GlucMed,
                  # Other drugs
                  AB_BA = H1_Antibiotica,
                  # Lab
                  Trig_BA = H1_Lab_UitslagTRIG, 
                  TC_BA = H1_Lab_UitslagCHOL,  
                  LDL_BA = H1_Lab_uitslagRLDL,  
                  HbA1c_BA = H1_Lab_UitslagIH1C, 
                  # Fecal sample 
                  SamplePresent = H1_FecesMicrobiome,
                  SampleAB_BA=H1_Feces_q2, SampleDiarrhoea_BA=H1_Feces_q3
    )

df_new2 <- df_new %>% 
    mutate(across(where(is.character), ~na_if(., c("Missing", "Missing: n.v.t.", "niet ingevuld","nvt", 
                                       "No lab result", "Missing: not applicable", "missing",
                                       "Missing: not measured", "missing",
                                       "See comments lab results", "Low (<1 mmol/L)",
                                       "Low (<0,08 mmol/L)", "Low (<0,10 mmol/L)",
                                       "Smoking status unknown", "Number or duration unknown",
                                       "Rookstatus onbekend", "Rookduur en/of aantal onbekend"))),
           sampleID_BA = str_c("HELIBA_", ID),
           sampleID_FU = str_c("HELIFU_", ID),
           ID = str_c("S",as.character(ID)),
           across(c("FUtime", "Age_BA", "Age_FU",
                   "BMI_BA",  "SBP_BA",  "DBP_BA", 
                    "Trig_BA", "TC_BA",
                    "LDL_BA","HbA1c_BA", ), as.numeric), 
           across(where(haven::is.labelled), ~haven::as_factor(.x, levels = "labels")),
           across(c("HT_BP_BA", "HT_BPMed_BA", "HT_SelfBP_BA", "HT_SelfBPMed_BA", "HT_Self_BA", "AB_BA"), yesnocaps),
           across(c("MetSyn_BA", "MetSyn_Obesity_BA",
                    "MetSyn_HighTG_BA", "MetSyn_LowHDL_BA", "MetSyn_HighGluc_BA",
                    "MetSyn_HighGluc_BA", "MetSyn_HighBP_BA"), yesnosmall),
           Sex = fct_recode(Sex, "Male" = "man", "Female" = "vrouw"),
           EthnicityTot = forcats::fct_recode(EthnicityTot, "Dutch" = "NL", "South-Asian Surinamese" = "Hind",
                                              "African Surinamese" = "Creools", "South-Asian Surinamese" = "Javaans",
                                              "Ghanaian" = "Ghanees", "Turkish" = "Turks", "Moroccan" = "Marokkaans",
                                              "Other"="Anders/onbekend",
                                              "Other"="Sur anders/onbekend"),
           across(c("Smoking_BA"), ~fct_recode(.x, "Yes" = "Ja", 
                                                             "Former smoking" = "Nee, maar vroeger wel",
                                                             "Never" = "Nee, ik heb nooit gerookt",
                                                             "Never" = "Nee, nooit gerookt")),
         Ethnicity = fct_recode(Ethnicity, "Dutch" = "Nederlands", "Surinamese" = "Surinaams",
                                "Turkish" = "Turks", "Moroccan" = "Marokkaans", "Ghanaian" = "Ghanees",
                                "Unknown" = "Onbekend", "Other" = "Anders"),
         FUtime = as.numeric(dmonths(FUtime), "years"),
         SamplePresent = case_when(SamplePresent == "Yes" ~ "Yes", SamplePresent == "No (insufficient counts or not sequenced)" ~ "No",
          .default = "No"),
         Microbiome = case_when(Microbiome == "1" ~ "Yes", .default = "No"),
         CBS = case_when(CBS == "1" ~ "Yes", .default = "No"),
         Followup = case_when(Physical_FU == "Wel H2-LO" ~ "Yes", TRUE ~ "No"),
         Dyslipidemia = case_when(MetSyn_HighTG_BA == "Yes" | MetSyn_LowHDL_BA == "Yes" ~ "Yes",
                               MetSyn_HighTG_BA == "No" & MetSyn_LowHDL_BA == "No" ~ "No")
    ) %>%
    droplevels(.)

dim(df_new2)

### Table total cohort with and without microbiome data ###
table1 <- df_new2 %>%
    dplyr::select(Age_BA, Sex, EthnicityTot, , BMI_BA, Smoking_BA, AlcCons_BA, 
           DM_BA, SBP_BA, DBP_BA, HT_BPMed_BA, MetSyn_BA, 
           TC_BA, LDL_BA, Trig_BA, 
           HbA1c_BA, SamplePresent) %>% 
    CreateTableOne(data=., strata = c("SamplePresent"), test = TRUE) %>% 
    print(nonnormal=c("Trig"), noSpaces = TRUE, 
          pDigits = 3, contDigits = 1, missing = TRUE) %>% 
    as.data.frame(.)
write.csv2(as.data.frame(table1), 'results/tables/table_totalcohort_alleth.csv')

table1 <- df_new2 %>% filter(EthnicityTot == "Dutch") %>%
    dplyr::select(Age_BA, Sex, BMI_BA, Smoking_BA, AlcCons_BA, 
           DM_BA, SBP_BA, DBP_BA, HT_BPMed_BA, MetSyn_BA, 
           TC_BA, LDL_BA, Trig_BA, 
           HbA1c_BA, SamplePresent) %>% 
    CreateTableOne(data=., strata = c("SamplePresent"), test = TRUE) %>% 
    print(nonnormal=c("Trig"), noSpaces = TRUE, 
          pDigits = 3, contDigits = 1, missing = TRUE) %>% 
    as.data.frame(.)
write.csv2(as.data.frame(table1), 'results/tables/table_totalcohort_dutch.csv')

table1 <- df_new2 %>% filter(EthnicityTot == "South-Asian Surinamese") %>%
    dplyr::select(Age_BA, Sex, BMI_BA, Smoking_BA, AlcCons_BA, 
           DM_BA, SBP_BA, DBP_BA, HT_BPMed_BA, MetSyn_BA, 
           TC_BA, LDL_BA, Trig_BA, 
           HbA1c_BA, SamplePresent) %>% 
    CreateTableOne(data=., strata = c("SamplePresent"), test = TRUE) %>% 
    print(nonnormal=c("Trig"), noSpaces = TRUE, 
          pDigits = 3, contDigits = 1, missing = TRUE) %>% 
    as.data.frame(.)
write.csv2(as.data.frame(table1), 'results/tables/table_totalcohort_sas.csv')

table1 <- df_new2 %>% filter(EthnicityTot == "African Surinamese") %>%
    dplyr::select(Age_BA, Sex, BMI_BA, Smoking_BA, AlcCons_BA, 
           DM_BA, SBP_BA, DBP_BA, HT_BPMed_BA, MetSyn_BA, 
           TC_BA, LDL_BA, Trig_BA, 
           HbA1c_BA, SamplePresent) %>% 
    CreateTableOne(data=., strata = c("SamplePresent"), test = TRUE) %>% 
    print(nonnormal=c("Trig"), noSpaces = TRUE, 
          pDigits = 3, contDigits = 1, missing = TRUE) %>% 
    as.data.frame(.)
write.csv2(as.data.frame(table1), 'results/tables/table_totalcohort_afrsur.csv')

table1 <- df_new2 %>% filter(EthnicityTot == "Turkish") %>%
    dplyr::select(Age_BA, Sex, BMI_BA, Smoking_BA, AlcCons_BA, 
           DM_BA, SBP_BA, DBP_BA, HT_BPMed_BA, MetSyn_BA, 
           TC_BA, LDL_BA, Trig_BA, 
           HbA1c_BA, SamplePresent) %>% 
    CreateTableOne(data=., strata = c("SamplePresent"), test = TRUE) %>% 
    print(nonnormal=c("Trig"), noSpaces = TRUE, 
          pDigits = 3, contDigits = 1, missing = TRUE) %>% 
    as.data.frame(.)
write.csv2(as.data.frame(table1), 'results/tables/table_totalcohort_turkish.csv')

table1 <- df_new2 %>% filter(EthnicityTot == "Moroccan") %>%
    dplyr::select(Age_BA, Sex, BMI_BA, Smoking_BA, AlcCons_BA, 
           DM_BA, SBP_BA, DBP_BA, HT_BPMed_BA, MetSyn_BA, 
           TC_BA, LDL_BA, Trig_BA, 
           HbA1c_BA, SamplePresent) %>% 
    CreateTableOne(data=., strata = c("SamplePresent"), test = TRUE) %>% 
    print(nonnormal=c("Trig"), noSpaces = TRUE, 
          pDigits = 3, contDigits = 1, missing = TRUE) %>% 
    as.data.frame(.)
write.csv2(as.data.frame(table1), 'results/tables/table_totalcohort_moroccan.csv')

table1 <- df_new2 %>% filter(EthnicityTot == "Ghanaian") %>%
    dplyr::select(Age_BA, Sex, BMI_BA, Smoking_BA, AlcCons_BA, 
           DM_BA, SBP_BA, DBP_BA, HT_BPMed_BA, MetSyn_BA, 
           TC_BA, LDL_BA, Trig_BA, 
           HbA1c_BA, SamplePresent) %>% 
    CreateTableOne(data=., strata = c("SamplePresent"), test = TRUE) %>% 
    print(nonnormal=c("Trig"), noSpaces = TRUE, 
          pDigits = 3, contDigits = 1, missing = TRUE) %>% 
    as.data.frame(.)
write.csv2(as.data.frame(table1), 'results/tables/table_totalcohort_ghanaian.csv')

#### 16S ####
df_new3 <- df_new2 |> filter(SamplePresent == "Yes")

#### With and without follow-up ###
table1 <- df_new3 %>%
      dplyr::select(Age_BA, Sex, BMI_BA, Smoking_BA, AlcCons_BA, 
           DM_BA, SBP_BA, DBP_BA, HT_BPMed_BA, MetSyn_BA, Dyslipidemia,
           TC_BA, LDL_BA, Trig_BA, 
           HbA1c_BA, Microbiome) %>% 
    CreateTableOne(data=., strata = c("Microbiome"), test = TRUE) %>% 
    print(nonnormal=c("Trig"), noSpaces = TRUE, 
          pDigits = 3, contDigits = 1, missing = TRUE) %>% 
    as.data.frame(.)
write.csv2(as.data.frame(table1), 'results/tables/table_withandwithoutfollowup.csv')
