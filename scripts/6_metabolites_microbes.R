## Associations metabolomics and microbes

library(rio)
library(haven)
library(tidyverse)
library(tableone)
library(ggsci)
library(ggrepel)
library(tidyselect)

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
selection <- rio::import("data/ASVs_CBS.xlsx", header = FALSE)
met <- readRDS("data/metabolomics/metabolomics_paired.RDS")
metsel <- met[str_detect(rownames(met), "HELIBA"),]
metsel <- as.data.frame(metsel)
metsel <- metsel %>% dplyr::select(!starts_with("X-"))
metsel <- metsel %>% mutate(across(everything(.), ~log10(.x+1)),
                            ID = rownames(.))
mb <- readRDS("data/16S/phyloseq/rarefied/phyloseq_rarefied.RDS")
sample_names(mb) <- str_remove(sample_names(mb),"_T1")
mb <- prune_samples(str_detect(sample_names(mb), "HELIBA"), mb)
mb <- prune_samples(sample_names(mb) %in% rownames(metsel), mb)
mb <- prune_taxa(selection$...1, mb)
tax <- readRDS("data/16s/phyloseq/rarefied/taxtable_rarefied.RDS")

mb <- as.data.frame(t(as(mb@otu_table, "matrix")))
head(mb)
mbsel <- mb[, selection$...1]
# colnames(mbsel) <- make.unique(tax$Tax[match(colnames(mbsel), tax$ASV)])
mbsel <- mbsel %>% mutate(across(everything(.), ~log10(.x+1)))
mbsel$ID <- rownames(mbsel)

mbmet <- left_join(mbsel, metsel, by = "ID")
dim(mbmet)

names(mbmet)[1:13]
res <- c()
mbmet$asv <- NULL
mbmet$met <- NULL
plothist <- c()
for(a in 1:13) {
    mbmet$asv <- mbmet[[a]]
    asvname <- names(mbmet)[a]
    plothist[[asvname]] <- gghistogram(mbmet$asv) + labs(x = asvname)
    for(b in 15:(ncol(mbmet) -1)) {
        mbmet$met <- mbmet[[b]]
        metname <- names(mbmet)[b]
        cor <- cor.test(mbmet$asv, y = mbmet$met, method = "spearman")
        row <- c(asvname, metname, cor$estimate[[1]], cor$p.value)
        names(row) <- c("ASV", "metabolite", "rho", "pval")
        res <- rbind(res, row)
        mbmet$met <- NULL
    }
    mbmet$asv <- NULL
}

ggarrange(plotlist = plothist, nrow = 4, ncol = 4, 
          labels = c(LETTERS[1:13]))
ggsave("results/asv_distribution.pdf", width = 12, height = 26)

res <- as.data.frame(res)
res2 <- res %>% mutate(across(c(3:4), as.numeric)) %>% 
    group_by(ASV) %>% 
    mutate(padj = p.adjust(pval, "fdr")) %>% 
    ungroup(.)

reswide_rho <- res2 %>% dplyr::select(ASV, metabolite, rho) %>% 
                    pivot_wider(., id_cols = ASV, names_from = metabolite, values_from = rho)
reswide_p <- res2 %>% dplyr::select(ASV, metabolite, pval) %>% 
                    pivot_wider(., id_cols = ASV, names_from = metabolite, values_from = pval)
reswide_padj <- res2 %>% dplyr::select(ASV, metabolite, padj) %>% 
                    pivot_wider(., id_cols = ASV, names_from = metabolite, values_from = padj)
head(reswide_p)

all(apply(reswide_padj[2:ncol(reswide_padj)], 1, function(x) sum(x < 0.05, na.rm = TRUE) > 0) == FALSE)
tk <- apply(reswide_padj[2:ncol(reswide_padj)], 2, function(x) sum(x < 0.05, na.rm = TRUE) > 0)
summary(tk)
reswide_rho <- reswide_rho[,c(TRUE,tk)]
reswide_p <- reswide_p[,c(TRUE,tk)]
reswide_padj <- reswide_padj[,c(TRUE,tk)]

plotlist <- c()
for(a in 1:nrow(reswide_rho)){
    taxname <- tax$Tax[which(tax$ASV == reswide_rho$ASV[a])]
    resrho <- reswide_rho[a,]
    resp <- reswide_p[a,]
    respadj <- reswide_padj[a,]
    res <- rbind(resrho, resp, respadj)
    res$ASV <- NULL
    rlong <- as.data.frame(t(res))
    colnames(rlong) <- c("rho", "pval", "padj")
    rlong$pval <- log10(rlong$pval) * -1
    rlong <- rlong %>% 
        mutate(met = rownames(.),
               psig = case_when(padj < 0.05 ~ "sig", .default = "not sig"),
               dir = case_when(rho < 0 ~ "negative correlation", .default = "positive correlation" ),
               sigdir = case_when(psig == "sig" ~ dir, .default = "not sig")
        )
    topmets <- rlong %>%
        filter(padj < 0.05) %>% 
        group_by(dir) %>% arrange(padj) %>% slice(1:15) %>% ungroup(.)
    rlong <- rlong %>% mutate(met = ifelse(met %in% topmets$met, met, NA))
    if(nrow(topmets) > 0){
        colvec <- c("negative correlation" = pal_aaas()(1), "not sig" = "grey", "positive correlation" = pal_aaas()(2)[2])
        plotlist[[reswide_rho$ASV[a]]] <- ggplot(data = rlong, aes(x = rho, y = pval, color = sigdir)) +
            geom_point() +
            geom_label_repel(aes(label = met), color = "black", size = 3,
                             min.segment.length = 0.2, force = 1, box.padding = 0.4) +
            scale_color_manual(values = colvec) +
            labs(title = taxname, x = "spearman's rho", y = "-log10(p-value)", color = "") +
            theme_Publication()
    }
}


ggarrange(plotlist = plotlist[1:9], nrow = 3, ncol = 3, labels = LETTERS[1:9], common.legend = TRUE, legend = "bottom")
ggsave("results/volcano_correlations.pdf", width = 20, height = 20)

infofile <- rio::import("data/metabolomics/Info_HELIUS_metabolomics.xlsx") %>% 
    dplyr::select(metabolite = CHEMICAL_NAME, superpathway = SUPER_PATHWAY, subpathway = SUB_PATHWAY)

tot <- res2 %>% filter(padj < 0.05) %>% arrange(padj) %>% 
                left_join(., tax %>% dplyr::select(ASV, Family, Genus, Tax), by = 'ASV') %>% 
                left_join(., infofile, by = 'metabolite')
nlevels(as.factor(tot$metabolite))

tot %>% filter(superpathway != "Xenobiotics") %>% dplyr::select(metabolite, Tax, superpathway, subpathway)

(asv_indole <- tot %>% filter(metabolite == "1H-indole-7-acetic acid"))
(asv_andro <- tot %>% filter(metabolite == "5alpha-androstan-3beta,17alpha-diol disulfate"))
summary(asv_tyr$ASV %in% asv_andro$ASV)
(asv_entero <- tot %>% filter(metabolite == "enterolactone sulfate"))
(asv_tyr <- tot %>% filter(metabolite == "3-methoxytyrosine"))
(asv_urso <- tot %>% filter(metabolite == "isoursodeoxycholate"))
summary(asv_tyr$ASV %in% asv_urso$ASV)
summary(asv_urso$ASV %in% asv_tyr$ASV)

tot %>% filter(metabolite == "cinnamoylglycine") %>% select(metabolite, rho, padj, Tax, subpathway)
tot %>% filter(metabolite == "3-phenylpropionate (hydrocinnamate)") %>% select(metabolite, rho, padj, Tax, subpathway)
tot %>% filter(metabolite == "enterolactone sulfate") %>% select(metabolite, rho, padj, Tax, subpathway)
tot %>% filter(metabolite == "4-hydroxycoumarin") %>% select(metabolite, rho, padj, Tax, subpathway)
tot %>% filter(metabolite == "3-methoxytyrosine") %>% select(metabolite, rho, padj, Tax, subpathway)
tot %>% filter(metabolite == "1H-indole-7-acetic acid") %>% select(metabolite, rho, padj, Tax, subpathway)
tot %>% filter(metabolite == "isoursodeoxycholate") %>% select(metabolite, rho, padj, Tax, subpathway)

totsel <- tot %>% arrange(pval) %>% slice(1:25)
plotlist3 <- c()
for(a in 1:nrow(totsel)) {
    asvname <- totsel$ASV[a]
    taxname <- totsel$Tax[a]
    metname <- totsel$metabolite[a]
    df <- mbmet
    df$asv <- mbmet[,asvname]
    df$met <- mbmet[,metname]
    plotlist3[[a]] <- ggplot(df, aes(x = asv, y = met)) +
                        geom_point(color = "royalblue2", alpha = 0.8) +
                        geom_smooth(method = "lm", color = "firebrick") +
                        stat_cor(method = "spearman") +
                        labs(x = taxname, y = metname) +
                        theme_Publication()
}
ggarrange(plotlist = plotlist3, nrow = 5, ncol = 5, labels = LETTERS[1:25])
ggsave("results/correlationplots.pdf", width = 16, height = 16)

plotlist3[[1]]
