# ijob -A berglandlab -c10 -p standard --mem=40G

#module load gcc/11.4.0 openmpi/4.1.4  R/4.3.1;R
### libraries

library(SeqArray)
library(gdsfmt)
library(SNPRelate)
library(ggplot2)
library(dplyr)
library(tidyr)





# ---- Input files ----
genofile.fn <- "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/trimmed10bp_masked_usobtusa.gds"
metadata <- read.csv("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/2022_2024seqmetadata20250811.csv", header = TRUE)
samplestats <- read.csv("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/usdobtusa_samplestats.csv")
mitotypes <- read.csv("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/mito_types.csv")  # Save as CSV

outdir <- "/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/"

# ---- Open GDS and PCA ----
gds <- seqOpen(genofile.fn)
pca_result <- snpgdsPCA(gds, num.thread = 10, autosome.only = FALSE)
seqClose(gds)

nPCs <- min(100, ncol(pca_result$eigenvect))
pc_scores <- data.frame(
  sample.id = pca_result$sample.id,
  pca_result$eigenvect[, 1:nPCs]
)
colnames(pc_scores)[-1] <- paste0("PC", 1:nPCs)

# ---- Load metadata ----

merged <- merge(pc_scores, metadata, by.x = "sample.id", by.y = "Well", all.x = TRUE)
merged2 <- merge(merged,samplestats , by.x = "sample.id", by.y = "sampleId", all.x = TRUE)
merged3 <- merge(merged2,mitotypes , by.x = "sample.id", by.y = "CloneA", all.x = TRUE)

merged2 <- merged3

# ---- Define predictors ----
# Continuous: e.g., MeanDepth, Missingness
# Categorical: e.g., Pool, Mitotype
continuous_vars <- c("meanDepth", "missingRate")
categorical_vars <- c("accuratelocation", "date", "Group")  # adjust to actual column names



# ---- Function for linear model & ANOVA ----
results <- data.frame()

pc_cols <- grep("^PC", colnames(merged2), value = TRUE)

# Initialize results list
results <- data.frame()

pca_merged <-merged2

head(pca_merged)


for(pc in pc_cols){
  
  # Continuous variables
  for(var in continuous_vars){
    formula <- as.formula(paste(pc, "~", var))
    model <- lm(formula, data = pca_merged)
    
    # Extract p-value
    pval <- summary(model)$coefficients[2,4]  
    
    # Extract correlation coefficient (Pearson r)
    r <- cor(pca_merged[[pc]], pca_merged[[var]], use = "complete.obs")
    
    results <- rbind(results, data.frame(
      PC = pc,
      Variable = var,
      Type = "Continuous",
      p.value = pval,
      correlation = r
    ))
  }
  
  # Categorical variables
  for(var in categorical_vars){
    formula <- as.formula(paste(pc, "~", var))
    model <- lm(formula, data = pca_merged)
    
    # ANOVA p-value
    anova_res <- anova(model)
    pval <- anova_res[1, "Pr(>F)"]
    
    results <- rbind(results, data.frame(
      PC = pc,
      Variable = var,
      Type = "Categorical",
      p.value = pval,
      correlation = NA  # correlation not defined for categorical variables
    ))
  }
}








# Check results
head(results)
# ---- Save results ----
write.csv(results,
          file.path(outdir, "pc_predictor_associations.csv"),
          row.names = FALSE, quote = FALSE)

# ---- Quick plotting: p-values across PCs ----
pcplotsPvalues <- ggplot(results, aes(x = as.numeric(gsub("PC","",PC)), y = -log10(p.value), color = Variable)) +
  geom_line() +
  geom_point() +
  theme_bw(base_size = 14) +
  labs(
    x = "Principal Component",
    y = "-log10(p-value)",
    title = "Association of PCs with metadata variables"
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red")

ggsave("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/pcplotsPvalues.png", plot = pcplotsPvalues, width = 7, height = 6, dpi = 300)










results$PCnum <- as.numeric(sub("PC([0-9]+).*", "\\1", results$PC))


pcplotsPvalues <- ggplot(results, aes(x = PCnum, y = -log10(p.value), color = Variable)) +
  geom_line() +
  geom_point() +
  theme_bw(base_size = 14) +
  labs(
    x = "Principal Component",
    y = "-log10(p-value)",
    title = "PCs x metadata variables"
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red")

# Save the figure
ggsave("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/pcplotsPvalues.png",
       plot = pcplotsPvalues, width = 7, height = 6, dpi = 300)


results_cont <- results[results$Type == "Continuous", ]

# Make sure PC numbers are numeric for plotting
results_cont$PCnum <- as.numeric(sub("PC([0-9]+).*", "\\1", results_cont$PC))


pcplotsCorrelation <- ggplot(results_cont, 
                             aes(x = PCnum, y = -log10(p.value), group=Variable, color=Variable)) +
  geom_point() +
  geom_line() +
  theme_bw(base_size = 14) +
  labs(
    x = "Principal Component",
    y = "-log10(p-value)",
    color = "Correlation",
    title = "PC X p-value & correlation"
  ) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black")


ggsave("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/pcplotsCorrelation.png",
       plot = pcplotsCorrelation, width = 8, height = 6, dpi = 300)
















results3 <- data.frame()

for (pc in pc_cols) {
  formula <- as.formula(paste(pc, "~ accuratelocation * date"))
  model <- lm(formula, data = pca_merged)
  
  # R² for full model
  r2 <- summary(model)$r.squared
  
  # ANOVA p-values for each term
  anova_res <- anova(model)
  
  for (term in rownames(anova_res)) {
    pval <- anova_res[term, "Pr(>F)"]
    
    results3 <- rbind(results3, data.frame(
      PC = pc,
      Term = term,
      R2 = r2,
      p.value = pval
    ))
  }
}

# View results
head(results3)


results3$log10p <- -log10(results3$p.value)

# Order PCs numerically (remove "PC" prefix if needed)
results3$PCnum <- as.numeric(gsub("PC", "", results3$PC))
results3 <- subset(results3, Term != "Residuals")

# Plot
d <- ggplot(results3, aes(x = PCnum, y = log10p, color = Term)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_line() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  theme_bw(base_size = 14) +
  labs(
    title = "PCs x pond, date, and interaction",
    x = "Principal Component",
    y = expression(-log[10](p~value))
  ) +
  facet_wrap(~ Term, ncol = 1, scales = "free_y")

  e <- ggplot(results3, aes(x = PCnum, y = R2)) +
  geom_point(size = 2, alpha = 0.8) +
  theme_bw(base_size = 14) +
  labs(
    title = "PC by R2",
    x = "PC model}",
  )

  
f <- d / e



ggsave("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/interactions_plot.png",
       plot = d, width = 8, height = 6, dpi = 300)

ggsave("/scratch/rjp5nc/UK2022_2024/daphnia_phylo/usdobtusa_indv/r2_bymodel_plot.png",
       plot = e, width = 8, height = 6, dpi = 300)