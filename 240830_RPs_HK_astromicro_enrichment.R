# Clear environment and console
rm(list = ls())
cat("\014")

# Load required libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(viridis)
library(gridExtra)
library(grid)
library(ggpubr)
library(biomaRt)

# Load count data (already CPM normalized)
data <- read.delim("raw_counts_cpm_after_norm_reorganised_220722.txt", sep="\t", check.names=FALSE)

# Use biomaRt to get gene symbols
ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
gene_info <- getBM(attributes=c("ensembl_gene_id", "external_gene_name"),
                   filters="ensembl_gene_id",
                   values=data$ensembl_id,
                   mart=ensembl)

# Merge gene symbols with the data
data <- merge(data, gene_info, by.x="ensembl_id", by.y="ensembl_gene_id", all.x=TRUE)

# Function to identify ribosomal protein genes
is_rp_gene <- function(gene_name) {
  (
    (grepl("^Rps", gene_name) & 
       !grepl("^Rps6k", gene_name) &
       !grepl("rt", gene_name, ignore.case = TRUE) &
       gene_name != "Rps19bp1") |
      grepl("^Rpl", gene_name) | 
      gene_name %in% c("Rpsa", "Rplp0", "Rplp1", "Rplp2") |
      grepl("Rps.*like", gene_name) | 
      grepl("Rpl.*like", gene_name)
  ) &
    !grepl("-ps", gene_name)
}

# Define initial housekeeping genes list
initial_housekeeping_genes <- c("Actb", "Gapdh", "B2m", "Ppia", "Hprt", "Tbp", "Ubc", "Ywhaz", "Pgk1")

# Check which housekeeping genes are present in the data
present_hk_genes <- initial_housekeeping_genes[initial_housekeeping_genes %in% data$external_gene_name]

# Print diagnostic information about housekeeping genes
cat("Housekeeping genes check:\n")
for (gene in initial_housekeeping_genes) {
  if (gene %in% data$external_gene_name) {
    cat(gene, "is present in the data.\n")
  } else {
    cat(gene, "is NOT found in the data.\n")
  }
}

# Update housekeeping genes list
housekeeping_genes <- present_hk_genes

cat("\nFinal list of housekeeping genes for analysis:\n")
print(housekeeping_genes)

# Identify RP genes
data$is_rp <- sapply(data$external_gene_name, is_rp_gene)

# Filter for RP genes and housekeeping genes separately
rp_data <- data[data$is_rp, ]
hk_data <- data[data$external_gene_name %in% housekeeping_genes, ]

# Identify count columns
count_columns <- 2:(ncol(data) - 3)  # Excluding ensembl_id, external_gene_name, and is_rp columns

# Extract normalized RP data
normalized_rp_data <- rp_data[, count_columns]
rownames(normalized_rp_data) <- rp_data$external_gene_name

# Extract normalized housekeeping gene data
normalized_hk_data <- hk_data[, count_columns]
rownames(normalized_hk_data) <- hk_data$external_gene_name

# Filter for RT data
rt_data <- normalized_rp_data %>%
  dplyr::select(contains("_RT_"))

# Calculate fractions for each RP gene
total_counts <- colSums(rt_data)
rp_fractions <- sweep(rt_data, 2, total_counts, "/")

# Reshape RP data for plotting
rp_fractions_long <- rp_fractions %>%
  rownames_to_column("gene") %>%
  pivot_longer(cols = -gene, names_to = "sample", values_to = "fraction") %>%
  mutate(
    cell_type = ifelse(grepl("^AC", sample), "Astrocyte", "Microglia"),
    region = case_when(
      grepl("_F_", sample) ~ "Front",
      grepl("_M_", sample) ~ "Middle",
      grepl("_R_", sample) ~ "Rear"
    )
  )

# Reshape housekeeping gene data for plotting
hk_long <- normalized_hk_data %>%
  dplyr::select(contains("_RT_")) %>%
  rownames_to_column("gene") %>%
  pivot_longer(cols = -gene, names_to = "sample", values_to = "expression") %>%
  mutate(
    cell_type = ifelse(grepl("^AC", sample), "Astrocyte", "Microglia"),
    region = case_when(
      grepl("_F_", sample) ~ "Front",
      grepl("_M_", sample) ~ "Middle",
      grepl("_R_", sample) ~ "Rear"
    )
  )

# Calculate mean fractions and perform statistical tests for RP genes
mean_fractions_rp <- rp_fractions_long %>%
  group_by(gene, cell_type, region) %>%
  summarize(
    mean_fraction = mean(fraction),
    se_fraction = sd(fraction) / sqrt(n()),
    .groups = "drop"
  )

stat_test_rp <- rp_fractions_long %>%
  group_by(gene, region) %>%
  summarise(
    p_value = t.test(fraction ~ cell_type)$p.value,
    .groups = "drop"
  ) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))

mean_fractions_rp <- mean_fractions_rp %>%
  left_join(stat_test_rp, by = c("gene", "region"))

# Calculate mean expression and perform statistical tests for housekeeping genes
mean_expression_hk <- hk_long %>%
  group_by(gene, cell_type, region) %>%
  summarize(
    mean_expression = mean(expression),
    se_expression = sd(expression) / sqrt(n()),
    .groups = "drop"
  )

stat_test_hk <- hk_long %>%
  group_by(gene, region) %>%
  summarise(
    p_value = t.test(expression ~ cell_type)$p.value,
    .groups = "drop"
  ) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))

mean_expression_hk <- mean_expression_hk %>%
  left_join(stat_test_hk, by = c("gene", "region"))

# Set custom colors
custom_colors <- c("Astrocyte" = "#9B7BB8", "Microglia" = "#5A9835")

# Function to create a boxplot for a single RP gene
create_rp_boxplot <- function(data, gene, region) {
  gene_data <- stat_test_rp %>%
    filter(gene == !!gene, region == !!region)
  
  p_value <- gene_data$p_value
  
  ggboxplot(data %>% filter(gene == !!gene, region == !!region), 
            x = "cell_type", y = "fraction", fill = "cell_type",
            palette = custom_colors,
            add = "jitter") +
    stat_compare_means(method = "t.test", label.y = max(data$fraction[data$gene == gene & data$region == region]) * 1.1) +
    labs(title = paste(gene, "\np-value =", format(p_value, digits = 3)),
         x = "Cell Type", y = "Fraction") +
    theme_minimal() +
    theme(legend.position = "none") +
    scale_y_continuous(limits = c(0, NA))  # This line sets the y-axis minimum to 0
}

# Function to create a boxplot for a single housekeeping gene
create_hk_boxplot <- function(data, gene, region) {
  gene_data <- stat_test_hk %>%
    filter(gene == !!gene, region == !!region)
  
  p_value <- gene_data$p_value
  
  ggboxplot(data %>% filter(gene == !!gene, region == !!region), 
            x = "cell_type", y = "expression", fill = "cell_type",
            palette = custom_colors,
            add = "jitter") +
    stat_compare_means(method = "t.test", label.y = max(data$expression[data$gene == gene & data$region == region]) * 1.1) +
    labs(title = paste(gene, "\np-value =", format(p_value, digits = 3)),
         x = "Cell Type", y = "Expression") +
    theme_minimal() +
    theme(legend.position = "none") +
    scale_y_continuous(limits = c(0, NA))  # This line sets the y-axis minimum to 0
}

# Create boxplots for top differentially expressed RP genes in each region
regions <- c("Front", "Middle", "Rear")
p_threshold <- 0.05

for (r in regions) {
  region_sig_genes <- stat_test_rp %>%
    filter(region == r, p_value < p_threshold) %>%
    top_n(9, wt = -p_value) %>%
    pull(gene)
  
  if (length(region_sig_genes) > 0) {
    region_plots <- lapply(region_sig_genes, function(g) create_rp_boxplot(rp_fractions_long, g, r))
    region_combined <- gridExtra::grid.arrange(
      grobs = region_plots, 
      ncol = 3,
      top = textGrob(paste(r, "Region - Top Differential RP Genes"), gp = gpar(fontsize = 16, fontface = "bold"))
    )
    ggsave(paste0("differential_rp_genes_boxplot_", tolower(r), ".png"), region_combined, width = 15, height = 15)
  }
}

# Create boxplots for housekeeping genes in each region
for (r in regions) {
  hk_plots <- lapply(housekeeping_genes, function(g) create_hk_boxplot(hk_long, g, r))
  hk_combined <- gridExtra::grid.arrange(
    grobs = hk_plots, 
    ncol = 3,
    top = textGrob(paste(r, "Region - Housekeeping Genes"), gp = gpar(fontsize = 16, fontface = "bold"))
  )
  ggsave(paste0("housekeeping_genes_boxplot_", tolower(r), ".png"), hk_combined, width = 15, height = 15)
}

# Create a heatmap of mean fractions for RP genes
mean_fractions_rp_wide <- mean_fractions_rp %>%
  dplyr::select(gene, cell_type, region, mean_fraction) %>%
  pivot_wider(names_from = c(cell_type, region), values_from = mean_fraction)

pheatmap(
  as.matrix(mean_fractions_rp_wide[,-1]),
  scale = "row",
  show_rownames = TRUE,
  labels_row = mean_fractions_rp_wide$gene,
  show_colnames = TRUE,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  color = viridis(100),
  main = "Mean Fractions of RP Genes",
  filename = "rp_genes_heatmap.png",
  width = 12,
  height = 16
)

# Create a heatmap of mean expression for housekeeping genes
mean_expression_hk_wide <- mean_expression_hk %>%
  dplyr::select(gene, cell_type, region, mean_expression) %>%
  pivot_wider(names_from = c(cell_type, region), values_from = mean_expression)

pheatmap(
  as.matrix(mean_expression_hk_wide[,-1]),
  scale = "row",
  show_rownames = TRUE,
  labels_row = mean_expression_hk_wide$gene,
  show_colnames = TRUE,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  color = viridis(100),
  main = "Mean Expression of Housekeeping Genes",
  filename = "housekeeping_genes_heatmap.png",
  width = 12,
  height = 8
)

# Save summary statistics and detailed results
write.csv(mean_fractions_rp, "rp_gene_fraction_glia.csv", row.names = FALSE)
write.csv(mean_expression_hk, "housekeeping_gene_expression_detailed.csv", row.names = FALSE)

# Create density plots of p-values
ggplot(stat_test_rp, aes(x = p_value)) +
  geom_density() +
  facet_wrap(~region) +
  labs(title = "Distribution of p-values for RP genes by region", x = "p-value", y = "Density") +
  theme_minimal()
ggsave("rp_p_value_distribution.png", width = 10, height = 6)

ggplot(stat_test_hk, aes(x = p_value)) +
  geom_density() +
  facet_wrap(~region) +
  labs(title = "Distribution of p-values for housekeeping genes by region", x = "p-value", y = "Density") +
  theme_minimal()
ggsave("hk_p_value_distribution.png", width = 10, height = 6)

# Print diagnostic information
cat("RP Data summary:\n")
print(summary(rp_fractions_long))

cat("\nHousekeeping Data summary:\n")
print(summary(hk_long))

cat("\nUnique regions in the data:\n")
print(unique(rp_fractions_long$region))

cat("\nNumber of RP genes per region:\n")
print(rp_fractions_long %>% group_by(region) %>% summarise(gene_count = n_distinct(gene)))

cat("\nNumber of housekeeping genes per region:\n")
print(hk_long %>% group_by(region) %>% summarise(gene_count = n_distinct(gene)))

cat("\nSample sizes for each region and cell type combination (RP genes):\n")
print(rp_fractions_long %>% group_by(region, cell_type) %>% summarise(sample_size = n_distinct(gene)))

cat("\nSample sizes for each region and cell type combination (housekeeping genes):\n")
print(hk_long %>% group_by(region, cell_type) %>% summarise(sample_size = n_distinct(gene)))

print("Analysis complete. Check the generated PNG files for visualizations and CSV files for detailed results.")



# Load required libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(grid)

# ... [previous data loading and processing steps remain unchanged]

# Modified function to create a boxplot for a single RP gene with significance asterisks
create_rp_boxplot <- function(data, gene, region) {
  gene_data <- stat_test_rp %>%
    filter(gene == !!gene, region == !!region)
  
  p_value <- gene_data$p_value
  
  significance <- case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01 ~ "**",
    p_value < 0.05 ~ "*",
    TRUE ~ "ns"
  )
  
  max_y <- max(data$fraction[data$gene == gene & data$region == region])
  y_position <- max_y * 1.1
  
  ggboxplot(data %>% filter(gene == !!gene, region == !!region), 
            x = "cell_type", y = "fraction", fill = "cell_type",
            palette = custom_colors,
            add = "jitter") +
    stat_compare_means(comparisons = list(c("Astrocyte", "Microglia")),
                       label = "p.signif",
                       method = "t.test",
                       label.y = y_position) +
    labs(title = paste0(gene, " ", significance),
         x = NULL, y = "Fraction") +
    theme_minimal() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
  axis.text.y = element_text(size = 12),  # Increased y-axis text size
  axis.title.y = element_text(size = 14)) +  # Increased y-axis title size
    scale_y_continuous(limits = c(0, y_position * 1.2))
}

# Create boxplots for top differentially expressed RP genes in each region
regions <- c("Front", "Middle", "Rear")
p_threshold <- 0.05

plot_list <- list()

for (r in regions) {
  region_sig_genes <- stat_test_rp %>%
    filter(region == r, p_value < p_threshold) %>%
    top_n(4, wt = -p_value) %>%
    pull(gene)
  
  if (length(region_sig_genes) > 0) {
    region_plots <- lapply(region_sig_genes, function(g) create_rp_boxplot(rp_fractions_long, g, r))
    plot_list[[r]] <- region_plots
  }
}

# Combine all plots
all_plots <- do.call(c, plot_list)

# Create a common legend
legend_plot <- ggplot(rp_fractions_long, aes(x = cell_type, y = fraction, fill = cell_type)) +
  geom_boxplot() +
  scale_fill_manual(values = custom_colors, name = "Cell Type") +
  theme_minimal() +
  theme(legend.position = "bottom")

legend <- get_legend(legend_plot)

# Arrange plots in a grid
arranged_plots <- gridExtra::grid.arrange(
  grobs = c(all_plots, list(legend)),
  ncol = 4,
  nrow = 4,
  layout_matrix = rbind(
    c(1,2,3,4),
    c(5,6,7,8),
    c(9,10,11,12),
    c(13,13,13,13)
  ),
  heights = c(1,1,1,0.2),
  top = textGrob("Top Differential RP Genes by Region", gp = gpar(fontsize = 20, fontface = "bold"))
)

# Save the plot
ggsave("differential_rp_genes_boxplot_combined.png", arranged_plots, width = 15, height = 15)

print("Analysis complete. Check the generated PNG file for the updated visualization.")



# Load required libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(grid)

# ... [previous data loading and processing steps remain unchanged]

# Function to create a boxplot for a single HK gene with significance asterisks
create_hk_boxplot <- function(data, gene, region) {
  gene_data <- stat_test_hk %>%
    filter(gene == !!gene, region == !!region)
  
  p_value <- gene_data$p_value
  
  significance <- case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01 ~ "**",
    p_value < 0.05 ~ "*",
    TRUE ~ "ns"
  )
  
  max_y <- max(data$expression[data$gene == gene & data$region == region])
  y_position <- max_y * 1.1
  
  ggboxplot(data %>% filter(gene == !!gene, region == !!region), 
            x = "cell_type", y = "expression", fill = "cell_type",
            palette = custom_colors,
            add = "jitter") +
    stat_compare_means(comparisons = list(c("Astrocyte", "Microglia")),
                       label = "p.signif",
                       method = "t.test",
                       label.y = y_position) +
    labs(title = paste0(gene, " ", significance),
         x = NULL, y = "Expression") +
    theme_minimal() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 10),
          axis.title.x = element_blank(),
          axis.text.x = element_blank()) +
    scale_y_continuous(limits = c(0, y_position * 1.2))
}

# Create boxplots for HK genes in each region
regions <- c("Front", "Middle", "Rear")

plot_list <- list()

for (r in regions) {
  region_plots <- lapply(housekeeping_genes, function(g) create_hk_boxplot(hk_long, g, r))
  plot_list[[r]] <- region_plots
}

# Combine all plots
all_plots <- do.call(c, plot_list)

# Create a common legend
legend_plot <- ggplot(hk_long, aes(x = cell_type, y = expression, fill = cell_type)) +
  geom_boxplot() +
  scale_fill_manual(values = custom_colors, name = "Cell Type") +
  theme_minimal() +
  theme(legend.position = "bottom")

legend <- get_legend(legend_plot)

# Determine the number of rows and columns
n_genes <- length(housekeeping_genes)
n_cols <- min(4, n_genes)
n_rows <- ceiling(n_genes / n_cols) * 3 + 1  # 3 regions + 1 row for legend

# Arrange plots in a grid
arranged_plots <- gridExtra::grid.arrange(
  grobs = c(all_plots, list(legend)),
  ncol = n_cols,
  nrow = n_rows,
  layout_matrix = rbind(
    matrix(1:(n_genes * 3), nrow = 3, byrow = TRUE),
    rep(n_genes * 3 + 1, n_cols)
  ),
  heights = c(rep(1, n_rows - 1), 0.2),
  top = textGrob("Housekeeping Genes by Region", gp = gpar(fontsize = 16, fontface = "bold"))
)

# Save the plot
ggsave("housekeeping_genes_boxplot_combined.png", arranged_plots, width = 15, height = 5 * n_rows)

print("Analysis complete. Check the generated PNG file for the updated visualization of housekeeping genes.")

