## Cleaning metabolomics data

library(tidyverse)

## Data
helius <- readRDS('data/clinicaldata_long.RDS')
infomet <- rio::import('data/metabolomics/Info_HELIUS_metabolomics.xlsx')
colnames(infomet)[1] <- "MetabolonID"
met <- rio::import('data/metabolomics/HELIUS_metabolomics_abundance.xlsx')
colnames(met)[1] <- "sampleID"
colnames(met)[2:ncol(met)] <- infomet$CHEMICAL_NAME[match(colnames(met)[2:ncol(met)],infomet$MetabolonID)]
meta <- rio::import("data/metabolomics/HELIUS_metabolomics_metadata.xlsx")
colnames(meta)[1] <- "sampleID"
head(meta$sampleID)

meta <- meta %>% mutate(
    subjectID = str_replace_all(CLIENT_IDENTIFIER, "Helius ", "HELIUS_"),
    subjectID = case_when(
        str_detect(CLIENT_IDENTIFIER, "Covid ") ~ str_replace(subjectID, "HELIUS", "HELICOV"),
        .default = str_replace(subjectID, "HELIUS", "HELIBA")
    ),
    subjectID = str_replace_all(subjectID, "Covid ", ""),
    timepoint = case_when(
        str_detect(subjectID, "HELICOV") ~ "follow-up",
        str_detect(subjectID, "HELIBA") ~ "baseline"
    ),
    timepoint = as.factor(timepoint),
    across(c("NEG", "POLAR", "POS.EARLY", "POS.LATE"), as.factor)
)

write.csv2(meta$subjectID, 'data/metabolomics/ids_metabolomics.csv')

summary(meta$NEG)
summary(meta$POLAR)
summary(meta$POS.EARLY)
summary(meta$POS.LATE)

dim(met)
met$subjectID <- meta$subjectID[match(met$sampleID, meta$sampleID)]
rownames(met) <- met$subjectID
met$subjectID <- NULL
met$sampleID <- NULL
metmatrix <- as.matrix(met)
constants <- apply(metmatrix, 2, var)
any(constants == 0)
print('In total there are this number of constants in the data:')
length(constants[which(constants == 0)])
nameconstants <- names(constants[which(constants == 0)])
write.csv(nameconstants, "results/constant_metabolites.csv")
const <- metmatrix[, constants == 0]
metmatrix <- metmatrix[,constants != 0]
dim(metmatrix) # so 75 of 1468 lost, 1393 left
missing <- apply(metmatrix, 2, is.na)
allmissing <- apply(missing, 2, all)
any(allmissing)
missing <- apply(metmatrix, 1, is.na)
allmissing <- apply(missing, 2, all)
any(allmissing)

rownames(metmatrix)
metdf <- as.data.frame(metmatrix)
metdf$Heliusnr <- str_remove(str_remove(rownames(metdf), "HELIBA_"), "HELICOV")
metdf$Timepoint <- str_extract(rownames(metdf), "HELI[A-Z]*")
metdf$Timepoint <- as.factor(case_when(metdf$Timepoint == "HELIBA" ~ "baseline", 
                             metdf$Timepoint == "HELICOV" ~ "covid"))
metdfba <- metdf[str_detect(rownames(metmatrix), "HELIBA"),]
dim(metmatrixba)
metdfba <- metdfba %>% relocate(Timepoint) %>% relocate(Heliusnr)

saveRDS(metmatrix, "data/metabolomics/metabolomics_paired.RDS")


######################## Metabolite info ########################

summary(infomet$SUPER_PATHWAY)
infomet <- infomet %>% filter(!CHEMICAL_NAME %in% nameconstants)
infomet <- infomet %>% select(MetabolonID, sup = SUPER_PATHWAY, sub = SUB_PATHWAY,
                              metname = PLOT_NAME, platform = PLATFORM) %>% 
    mutate(sup = case_when(
        is.na(sup) ~ "Unknown pathway",
        .default = sup
            )
)
saveRDS(infomet, "data/info_metabolites.RDS")




