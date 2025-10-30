###############################################################
# LAB PRACTICAL: MORPHOLOGY AND CLONALITY OF NORMAL VS TUMOUR CELLS
# Author: Hasmik Chilingaryan
# Purpose: Generate Figures & Perform Statistical Analysis
# ---------------------------------------------------------------
# Figures:
# 1. Bar Graph – White Blood Cell Composition (%)
# 2. Boxplot – Lymphocyte % vs. Sample Type
# 3. Statistical Analyses for Both Graphs (ANOVA / Kruskal-Wallis)
###############################################################

install.packages(c("ggplot2", "reshape2", "dplyr"))
library(ggplot2)
library(reshape2)
library(dplyr)
install.packages("ggpubr")
library("ggpubr")

###############################################################
# Figure 1: BAR GRAPH - Overall WBC Composition (%)
###############################################################

# Load required libraries
library(ggplot2)
library(reshape2)

# Observed counts matrix
wbc_counts <- matrix(
  c(
    83,123,16,   # Neutrophil
    14,21,2,   # Eosinophil
    0,4,0,      # Basophil
    10,21,7,    # Monocyte
    118,62,167, # Lymphocyte
    3,15,3,     # Atypical Lymphocyte
    9,0,94     # Lymphoblast (ALL)
  ),
  nrow = 7,
  byrow = TRUE
)
rownames(wbc_counts) <- c("Neutrophil","Eosinophil","Basophil",
                          "Monocyte","Lymphocyte","Atypical Lymphocyte",
                          "Lymphoblast (ALL)")
colnames(wbc_counts) <- c("A","B","C")

# Run Chi-square test
chi_result <- chisq.test(wbc_counts)
print(chi_result)

# WBC composition from class total observations
wbc <- data.frame(
  CellType = c("Neutrophil", "Eosinophil", "Basophil", 
               "Monocyte", "Lymphocyte", "Atypical Lymphocyte", "Lymphoblast (ALL)"),
  A = c(35.0, 5.9, 0.0, 4.2, 49.8, 1.3, 3.8),
  B = c(50.0, 8.5, 1.6, 8.5, 25.2, 6.1, 0.0),
  C = c(5.5, 0.7, 0.0, 2.4, 57.8, 1.0, 32.5)
)

# Convert to long format for ggplot
wbc_long <- melt(wbc, id.vars = "CellType",
                 variable.name = "Sample", value.name = "Percentage")

# Bar graph
wbc_long$CellType <- factor(
  wbc_long$CellType,
  levels = c("Neutrophil", "Eosinophil", "Basophil",
             "Monocyte", "Lymphocyte", 
             "Atypical Lymphocyte", "Lymphoblast (ALL)"),
  labels = c("Neutrophil", "Eosinophil", "Basophil",
             "Monocyte", "Lymphocyte",
             "Atypical\nLymphocyte", "Lymphoblast\n(ALL)")
)

figure1 <- ggplot(wbc_long, aes(y = CellType, x = Percentage, fill = Sample)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  coord_flip() +
  labs(
    title = "Distribution of White Blood Cell Types Across Samples",
    y = "",
    x = "% of Cell Count"
    ) +
  scale_fill_manual(
    values = c("A" = "#5DADE2", "B" = "#58D68D", "C" = "#EC7063"),
    name = "Samples"
    ) +
  scale_x_continuous(breaks = seq(0, 60, by = 10)) +
  theme_minimal(base_size = 13)

figure1 <- figure1 +
  annotate("segment",
           x = 45, xend = 45,    # Adjust horizontal range of the line
           y = 6.75, yend = 7.4,      # y = 7 corresponds to "Lymphoblast (ALL)" row
           size = 1.1, colour = "black") +
  annotate("text",
           x = 45.5, y = 7.1,   # Midpoint above the line
           label = "*****",
           size = 6, colour = "black")

# Save as a high-resolution PNG
ggsave("Figure1_WBC_Distribution.png", plot = figure1, width = 9, height = 6, dpi = 300)

# Optional: Save as PDF (for vector quality)
ggsave("Figure1_WBC_Distribution.pdf", plot = figure1, width = 9, height = 6)

###############################################################
# FIGURE 2: BOX PLOT — LYMPHOCYTE % VS SAMPLE TYPE
###############################################################

lymphocyte_data <- data.frame(
  Sample = rep(c("A","B","C"), each = 9),
  Lymphocyte = c(
    # Sample A (9 students)
    36.4, 40.0, 25.0, 44.4, 48.0, 32.0, 60.0, 100.0, 60.0,
    # Sample B (9 students)
    10.0, 24.0, 36.0, 28.6, 16.0, 28.0, 40.0, 11.1, 44.4,
    # Sample C (9 students)
    90.0, 80.0, 80.0, 73.5, 48.0, 44.0, 100.0, 100.0, 8.3
  )
)

# Statistical analysis
# 1) Normality (Shapiro–Wilk per sample)
shapiro_by_sample <- by(lymphocyte_data$Lymphocyte,
                        lymphocyte_data$Sample, shapiro.test)
print(shapiro_by_sample)

# 2) Parametric: one-way ANOVA + Tukey HSD
anova_result <- aov(Lymphocyte ~ Sample, data = lymphocyte_data)
summary(anova_result)
tukey_result <- TukeyHSD(anova_result)
print(tukey_result)

# 3) Non-parametric: Kruskal–Wallis + pairwise Wilcoxon
kruskal_result <- kruskal.test(Lymphocyte ~ Sample, data = lymphocyte_data)
print(kruskal_result)

pairwise_wilcox <- pairwise.wilcox.test(
  lymphocyte_data$Lymphocyte, lymphocyte_data$Sample,
  p.adjust.method = "bonferroni"
)
print(pairwise_wilcox)

# 4) Descriptive statistics
desc_stats <- lymphocyte_data %>%
  group_by(Sample) %>%
  summarise(
    n = n(),
    mean = mean(Lymphocyte),
    sd = sd(Lymphocyte),
    median = median(Lymphocyte),
    Q1 = quantile(Lymphocyte, 0.25),
    Q3 = quantile(Lymphocyte, 0.75)
  )

print(desc_stats)

# Boxplot with jittered individual data points
figure2 <- ggplot(lymphocyte_data, aes(x = Sample, y = Lymphocyte, fill = Sample)) +
  geom_boxplot(alpha = 0.8, width = 0.6, colour = "black", size = 0.8,
               outlier.shape = 16, outlier.size = 2) +
  geom_jitter(width = 0.12, alpha = 0.8, size = 2,
              aes(color = Sample)) +
  labs(title = "Lymphocyte Percentage Across Sample Types",
       x = "Samples", y = "Lymphocyte (%)") +
  scale_fill_manual(values = c("A"="#5DADE2","B"="#58D68D","C"="#EC7063")) +
  scale_color_manual(values = c("A"="#1F4E79","B"="#2E8B57","C"="#A93226")) +
  scale_y_continuous(breaks = seq(0, 100, by = 20)) +   # ✅ Y-axis intervals every 10%
  theme_minimal(base_size = 13) +
  theme(legend.position = "none") +
  # Overall ANOVA p-value
  stat_compare_means(method = "anova", label.y = 110) +
  # Pairwise Tukey significance
  stat_compare_means(comparisons = list(c("B","C")), method = "t.test",
                     label = "p.signif", label.y = 100)

print(figure2)

# Use ANOVA p-value for the overall comparison. Add Tukey significance for pairwise C vs B

# Save boxlpot
ggsave("Figure2_Lymphocyte_boxplot.png", plot = figure2, width = 7, height = 6, dpi = 300)
ggsave("Figure2_Lymphocyte_boxplot.pdf", plot = figure2, width = 7, height = 6)

###############################################################
# FIGURE 3: BOX PLOT — LYMPHOCYTE & LYMPHOBLAST % VS SAMPLE TYPE
###############################################################

library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)

# Original data including lymphocytes and lymphoblasts
lymph_data <- data.frame(
  Sample = rep(c("A","B","C"), each = 9),
  Lymphocyte = c(
    # Sample A
    36.4, 40.0, 25.0, 44.4, 48.0, 32.0, 60.0, 100.0, 60.0,
    # Sample B
    10.0, 24.0, 36.0, 28.6, 16.0, 28.0, 40.0, 11.1, 44.4,
    # Sample C
    90.0, 80.0, 80.0, 73.5, 48.0, 44.0, 100.0, 100.0, 8.3
  ),
  Lymphoblast = c(
    # Sample A
    0.0, 20.0, 16.7, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    # Sample B
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    # Sample C
    0.0, 0.0, 0.0, 0.0, 52.0, 56.0, 100.0, 100.0, 70.8
  )
)

# Convert to long format for faceting
lymph_long <- lymph_data %>%
  pivot_longer(cols = c(Lymphocyte, Lymphoblast),
               names_to = "CellType", values_to = "Percentage")

# Optional: nicer facet labels
lymph_long$CellType <- factor(lymph_long$CellType, 
                              levels = c("Lymphocyte", "Lymphoblast"),
                              labels = c("Lymphocyte", "Lymphoblast (ALL)"))

# Boxplot with jittered points and facets
figure3 <- ggplot(lymph_long, aes(x = Sample, y = Percentage, fill = Sample)) +
  geom_boxplot(alpha = 0.8, width = 0.6, colour = "black", size = 0.8,
               outlier.shape = 16, outlier.size = 2) +
  geom_jitter(width = 0.12, alpha = 0.8, size = 2,
              aes(color = Sample)) +
  facet_wrap(~CellType, scales = "free_y") +
  labs(title = "Comparison of Lymphocytes vs Lymphoblasts Across Samples",
       x = "Samples", y = "Percentage (%)") +
  scale_fill_manual(values = c("A"="#5DADE2","B"="#58D68D","C"="#EC7063")) +
  scale_color_manual(values = c("A"="#1F4E79","B"="#2E8B57","C"="#A93226")) +
  scale_y_continuous(breaks = seq(0, 120, by = 20)) +   # ✅ consistent y-axis ticks
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 13, face = "bold")
  ) +
  # Overall ANOVA p-values per facet
  stat_compare_means(method = "anova", label.y = 110) +
  # Pairwise comparison (B vs C)
  stat_compare_means(comparisons = list(c("B","C")), method = "t.test",
                     label = "p.signif", label.y = 100)

# Print figure
print(figure3)

# Save figure
ggsave("Figure3_Lymphocyte_Lymphoblast_boxplot.png", plot = figure3, width = 8, height = 6, dpi = 300)
ggsave("Figure3_Lymphocyte_Lymphoblast_boxplot.pdf", plot = figure3, width = 8, height = 6)
