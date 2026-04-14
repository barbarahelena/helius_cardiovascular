# Data cleaning of clinical data, 16S, shotgun
## Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

## Libraries
library(tidyverse)
library(vegan)
library(phyloseq)
library(ggplot2)

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

# Get data
ps <- readRDS("data/16s/phyloseq/complete/phyloseq.RDS")
otu <- t(as(otu_table(ps), "matrix"))

# Rarefaction curve
rarecurve(
  otu,
  step = 1000,
  sample = min(rowSums(otu)),
  label = FALSE
)

# Generate rarefaction data
rare_data <- lapply(1:nrow(otu), function(i){
  rarecurve(otu[i, , drop = FALSE], step = 1000, sample = min(rowSums(otu)), label = FALSE)
})

# Convert rarecurve output into a tidy data frame
df <- bind_rows(lapply(1:length(rare_data), function(i){
  tibble(
    sample = rownames(otu)[i],
    reads = as.numeric(str_remove(names(rare_data[[i]][[1]]), "N")),
    richness = as.numeric(rare_data[[i]][[1]])
  )
}))

df_summary <- df %>%
  filter(reads <= 40000) %>%
  filter(reads %in% c(0, seq(1001, 40001, by = 5000))) %>%
  group_by(reads) %>%
  summarise(
    mean_richness = mean(richness),
    lowerse = mean(richness) - sd(richness) / sqrt(n()),
    higher = mean(richness) + sd(richness) / sqrt(n())
  )

ggplot(df_summary, aes(x = reads, y = mean_richness)) +
  geom_line(color = "steelblue") +
  geom_smooth(method = "loess", se = FALSE, color = "darkblue", linetype = "dashed", span = 0.3) +
  geom_ribbon(aes(ymin = lowerse, ymax = higher), fill = "steelblue", alpha = 0.2) +
  # scale_x_log10() +
  scale_y_log10() +
  coord_cartesian(xlim = c(0, 40000)) + 
  theme_Publication() +
  labs(x = "Sequencing depth", y = "ASV richness")

ggplot(df, aes(x = reads, y = richness)) +
  geom_line(aes(group = sample), alpha = 0.3, color = "steelblue") +
  geom_line(data = df_summary, aes(x = reads, y = mean_richness), color = "darkblue", size = 1) +
  theme_Publication() +
  scale_y_log10() +
  coord_cartesian(xlim = c(0, 40000)) + 
  labs(x = "Sequencing depth", y = "ASV richness") + 
  geom_vline(xintercept = 15000, linetype = "dashed", color = "black")
ggsave("results/rarefactioncurve.pdf", width = 6, height = 4)
