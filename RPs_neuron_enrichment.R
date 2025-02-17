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
library(writexl)
library(biomaRt)
library(ggpubr)
library(grid)

# Load neuron count data and sample information
neuron_counts <- read.csv("SB_merged_gene_counts_all.csv", row.names = 1)
sample_info <- read.csv("SB_sample_Info.csv")
head(sample_info)
# Create a mapping between Sample_ID and Sample_Name
sample_mapping <- sample_info %>%
  group_by(Sample_Name) %>%
  mutate(replicate = row_number()) %>%
  ungroup() %>%
  mutate(new_column_name = paste(Sample_Name, replicate, sep = "_rep")) %>%
  dplyr::select(Sample_ID, new_column_name)

# Rename the columns in neuron_counts
new_colnames <- sample_mapping$new_column_name[match(colnames(neuron_counts), sample_mapping$Sample_ID)]
colnames(neuron_counts) <- new_colnames

# Use biomaRt to get gene symbols
ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
gene_info <- getBM(attributes=c("ensembl_gene_id", "external_gene_name"),
                   filters="ensembl_gene_id",
                   values=rownames(neuron_counts),
                   mart=ensembl)

# Merge gene symbols with the data
neuron_data <- neuron_counts %>%
  rownames_to_column("ensembl_id") %>%
  left_join(gene_info, by = c("ensembl_id" = "ensembl_gene_id"))

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

# Identify RP genes
neuron_data$is_rp <- sapply(neuron_data$external_gene_name, is_rp_gene)

# Filter for RP genes and RT samples
rp_data <- neuron_data %>%
  filter(is_rp) %>%
  dplyr::select(external_gene_name, contains("_RT_"))

# Calculate total RP counts for each sample
total_rp_counts <- colSums(rp_data[, -1])

# Calculate fractions for each RP gene
rp_fractions <- sweep(rp_data[, -1], 2, total_rp_counts, "/")
rp_fractions$Gene <- rp_data$external_gene_name

# Reshape data for plotting
rp_fractions_long <- rp_fractions %>%
  pivot_longer(cols = -Gene, names_to = "Sample", values_to = "Fraction") %>%
  mutate(
    Cell_Type = case_when(
      grepl("vGluT2", Sample) ~ "vGluT2",
      grepl("Gad2", Sample) ~ "Gad2",
      grepl("PV", Sample) ~ "PV",
      grepl("SST", Sample) ~ "SST"
    ),
    Region = ifelse(grepl("_Cb_", Sample), "Cerebellum", "Cerebrum"),
    Disease = case_when(
      grepl("HD_", Sample) ~ "HD",
      grepl("CJD_", Sample) ~ "CJD",
      grepl("FFI_", Sample) ~ "FFI",
      grepl("wt_", Sample) ~ "Wild-type"
    )
  )
head(rp_fractions_long)
str(rp_fractions_long)

# Function to perform comparisons
perform_comparisons <- function(data, region) {
  if (region == "Cerebellum") {
    # Compare vGluT2 vs Gad2 in Cerebellum
    data <- data %>% filter(Cell_Type %in% c("vGluT2", "Gad2"))
    if (nrow(data) < 2 || length(unique(data$Cell_Type)) < 2) {
      return(data.frame(Gene = unique(data$Gene), Comparison = "vGluT2 vs Gad2", p_value = NA, mean_diff = NA))
    }
    t_test_result <- tryCatch({
      t.test(Fraction ~ Cell_Type, data = data)
    }, error = function(e) NULL)
    if (is.null(t_test_result)) {
      return(data.frame(Gene = unique(data$Gene), Comparison = "vGluT2 vs Gad2", p_value = NA, mean_diff = NA))
    }
    return(data.frame(
      Gene = unique(data$Gene),
      Comparison = "vGluT2 vs Gad2",
      p_value = t_test_result$p.value,
      mean_diff = diff(t_test_result$estimate)
    ))
  } else if (region == "Cerebrum") {
    # Compare vGluT2 vs Gad2 and PV vs SST separately in Cerebrum
    vglut2_gad2 <- data %>% filter(Cell_Type %in% c("vGluT2", "Gad2"))
    pv_sst <- data %>% filter(Cell_Type %in% c("PV", "SST"))
    
    results <- list()
    
    if (nrow(vglut2_gad2) >= 2 && length(unique(vglut2_gad2$Cell_Type)) == 2) {
      t_test_result <- tryCatch({
        t.test(Fraction ~ Cell_Type, data = vglut2_gad2)
      }, error = function(e) NULL)
      if (!is.null(t_test_result)) {
        results$vglut2_gad2 <- data.frame(
          Gene = unique(data$Gene),
          Comparison = "vGluT2 vs Gad2",
          p_value = t_test_result$p.value,
          mean_diff = diff(t_test_result$estimate)
        )
      }
    }
    
    if (nrow(pv_sst) >= 2 && length(unique(pv_sst$Cell_Type)) == 2) {
      t_test_result <- tryCatch({
        t.test(Fraction ~ Cell_Type, data = pv_sst)
      }, error = function(e) NULL)
      if (!is.null(t_test_result)) {
        results$pv_sst <- data.frame(
          Gene = unique(data$Gene),
          Comparison = "PV vs SST",
          p_value = t_test_result$p.value,
          mean_diff = diff(t_test_result$estimate)
        )
      }
    }
    
    return(do.call(rbind, results))
  }
}

# Function to process data and create plots
process_and_plot <- function(data, condition) {
  # Calculate mean fractions
  mean_fractions <- data %>%
    group_by(Gene, Cell_Type, Region) %>%
    summarize(
      Mean_Fraction = mean(Fraction),
      SD_Fraction = sd(Fraction),
      .groups = "drop"
    )
  
  # Perform comparisons
  results <- data %>%
    group_by(Gene, Region) %>%
    do(perform_comparisons(data = ., region = .$Region[1])) %>%
    ungroup()
  
  # Apply Benjamini-Hochberg correction
  results <- results %>%
    group_by(Region, Comparison) %>%
    mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
    ungroup()
  
  # Define custom colors
  custom_colors <- c("vGluT2" = "#CC79A7", "Gad2" = "#0072B2", "PV" = "#009E73", "SST" = "#F0E442")
  
  # Function to create a bar plot for a single region
  create_region_barplot <- function(data, region_name, y_axis_limit) {
    ggplot(data, 
           aes(x = reorder(Gene, Mean_Fraction), y = Mean_Fraction, fill = Cell_Type)) +
      geom_bar(stat = "identity", position = position_dodge()) +
      geom_errorbar(aes(ymin = Mean_Fraction - SD_Fraction, 
                        ymax = Mean_Fraction + SD_Fraction),
                    position = position_dodge(0.9), width = 0.25) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
            legend.position = "top") +
      labs(title = paste(region_name, "Region -", condition),
           x = "Gene", y = "Mean Fraction", fill = "Cell Type") +
      scale_fill_manual(values = custom_colors) +
      ylim(0, y_axis_limit) +
      geom_text(aes(label = ifelse(p_adj < 0.05, "*", ""),
                    y = Mean_Fraction + SD_Fraction),
                position = position_dodge(0.9),
                vjust = -0.5)
  }
  
  # Create bar plots for each region
  regions <- c("Cerebellum", "Cerebrum")
  plots <- list()
  
  for (region in regions) {
    region_data <- mean_fractions %>%
      filter(Region == region) %>%
      left_join(results, by = c("Gene", "Region"))
    
    y_axis_limit <- max(region_data$Mean_Fraction + region_data$SD_Fraction, na.rm = TRUE) * 1.1
    
    p <- create_region_barplot(region_data, region, y_axis_limit)
    plots[[region]] <- p
    ggsave(paste0("rp_gene_fractions_", tolower(region), "_", tolower(condition), ".png"), p, width = 15, height = 8)
  }
  
  # Combine plots vertically
  combined_plot <- do.call(gridExtra::grid.arrange, c(plots, ncol = 1))
  ggsave(paste0("rp_gene_fractions_combined_vertical_", tolower(condition), ".png"), combined_plot, width = 15, height = 16)
  
  # Function to create a boxplot for a single gene
  create_gene_boxplot <- function(plot_data, gene, region, comparison) {
    gene_data <- results %>%
      filter(Gene == !!gene, Region == !!region, Comparison == !!comparison)
    
    p_value <- gene_data$p_value[1]
    
    cell_types <- strsplit(comparison, " vs ")[[1]]
    plot_data_filtered <- plot_data %>% 
      filter(Gene == !!gene, Region == !!region, Cell_Type %in% cell_types)
    
    ggboxplot(plot_data_filtered, 
              x = "Cell_Type", y = "Fraction", fill = "Cell_Type",
              palette = custom_colors,
              add = "jitter") +
      stat_compare_means(label.y = max(plot_data_filtered$Fraction, na.rm = TRUE) * 1.1) +
      labs(title = paste(gene, "-", condition, "\n", comparison, "\np-value =", format(p_value, digits = 3)),
           x = "Cell Type", y = "Fraction") +
      theme_minimal() +
      theme(legend.position = "none")
  }
  
  # Create boxplots for top differentially expressed genes in each region and comparison
  p_threshold <- 0.05
  for (r in regions) {
    if (r == "Cerebellum") {
      comparisons <- c("vGluT2 vs Gad2")
    } else {
      comparisons <- c("vGluT2 vs Gad2", "PV vs SST")
    }
    
    for (comp in comparisons) {
      region_sig_genes <- results %>%
        filter(Region == r, Comparison == comp, p_adj < p_threshold) %>%
        top_n(9, wt = -p_value) %>%
        pull(Gene)
      
      if (length(region_sig_genes) > 0) {
        region_plots <- lapply(region_sig_genes, function(g) create_gene_boxplot(data, g, r, comp))
        region_combined <- gridExtra::grid.arrange(
          grobs = region_plots, 
          ncol = 3,
          top = textGrob(paste(r, "Region -", comp, "- Top Differential Genes -", condition), gp = gpar(fontsize = 16, fontface = "bold"))
        )
        ggsave(paste0("differential_rps_boxplot_", tolower(r), "_", gsub(" vs ", "_vs_", tolower(comp)), "_", tolower(condition), ".png"), region_combined, width = 15, height = 15)
      }
    }
  }
  
  # Return results for saving
  return(list(results = results, mean_fractions = mean_fractions))
}

# Process and plot wild-type data
wt_data <- rp_fractions_long %>% filter(Disease == "Wild-type")
wt_results <- process_and_plot(wt_data, "Wild-type")

# Process and plot each disease separately
diseases <- c("HD", "CJD", "FFI")
disease_results <- list()

for (disease in diseases) {
  disease_data <- rp_fractions_long %>% filter(Disease == disease)
  disease_results[[disease]] <- process_and_plot(disease_data, disease)
}

# Save results
results_list <- list(
  "Wild_Type_Results" = wt_results$results,
  "Wild_Type_Mean_Fractions" = wt_results$mean_fractions
)

for (disease in diseases) {
  results_list[[paste0(disease, "_Results")]] <- disease_results[[disease]]$results
  results_list[[paste0(disease, "_Mean_Fractions")]] <- disease_results[[disease]]$mean_fractions
}

write_xlsx(results_list, "neuron_rp_enrichment_comparison_results_all_conditions.xlsx")

print("Analysis complete. Check the generated PNG files for visualizations and the Excel file for detailed results.")



# Load required libraries
library(dplyr)
library(ggplot2)
library(ggpubr)
library(gridExtra)

# Define custom colors for cell types
custom_colors <- c("vGluT2" = "#CC79A7", "Gad2" = "#0072B2", "PV" = "#009E73", "SST" = "#F0E442")

# Modified function to create a boxplot for a single RP gene with significance asterisks
create_rp_boxplot <- function(data, gene, region, comparison) {
  plot_data <- data %>% filter(Gene == gene, Region == region, Cell_Type %in% strsplit(comparison, " vs ")[[1]])
  
  p_value <- wilcox.test(Fraction ~ Cell_Type, data = plot_data)$p.value
  
  significance <- case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01 ~ "**",
    p_value < 0.05 ~ "*",
    TRUE ~ "ns"
  )
  
  max_y <- max(plot_data$Fraction)
  y_position <- max_y * 1.1
  
  ggboxplot(plot_data, 
            x = "Cell_Type", y = "Fraction", fill = "Cell_Type",
            palette = custom_colors,
            add = "jitter") +
    stat_compare_means(comparisons = list(strsplit(comparison, " vs ")[[1]]),
                       label = "p.signif",
                       method = "wilcox.test",
                       label.y = y_position) +
    labs(title = paste0(gene, " ", significance),
         x = NULL, y = "Fraction") +
    theme_minimal() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.y = element_text(size = 14)) +
    scale_y_continuous(limits = c(0, y_position * 1.2))
}

# Define regions and comparisons
regions <- c("Cerebellum", "Cerebrum", "Cerebrum")
comparisons <- c("vGluT2 vs Gad2", "vGluT2 vs Gad2", "PV vs SST")

# Set p-value threshold
p_threshold <- 0.05

# Create boxplots for top differentially expressed RP genes in each region and comparison
plot_list <- list()

for (i in 1:length(regions)) {
  r <- regions[i]
  comp <- comparisons[i]
  
  region_sig_genes <- wt_results$results %>%
    filter(Region == r, Comparison == comp, p_adj < p_threshold) %>%
    top_n(4, wt = -p_value) %>%
    pull(Gene)
  
  if (length(region_sig_genes) > 0) {
    region_plots <- lapply(region_sig_genes, function(g) create_rp_boxplot(wt_data, g, r, comp))
    plot_list[[paste(r, comp)]] <- region_plots
  }
}

# Combine all plots
all_plots <- do.call(c, plot_list)

# Create a common legend
legend_plot <- ggplot(wt_data, aes(x = Cell_Type, y = Fraction, fill = Cell_Type)) +
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
  top = textGrob("Top Differential RP Genes by Region and Cell Type Comparison", gp = gpar(fontsize = 20, fontface = "bold"))
)

# Save the plot
ggsave("differential_rp_genes_boxplot_neurons_combined.png", arranged_plots, width = 15, height = 15)

print("Analysis complete. Check the generated PNG file for the updated visualization of neuron RP genes.")




# Load required libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(viridis)
library(gridExtra)
library(writexl)
library(biomaRt)
library(ggpubr)
library(grid)

# Load neuron count data and sample information
neuron_counts <- read.csv("SB_merged_gene_counts_all.csv", row.names = 1)
sample_info <- read.csv("SB_sample_Info.csv")

# Create a mapping between Sample_ID and Sample_Name
sample_mapping <- sample_info %>%
  group_by(Sample_Name) %>%
  mutate(replicate = row_number()) %>%
  ungroup() %>%
  mutate(new_column_name = paste(Sample_Name, replicate, sep = "_rep")) %>%
  dplyr::select(Sample_ID, new_column_name)

# Rename the columns in neuron_counts
new_colnames <- sample_mapping$new_column_name[match(colnames(neuron_counts), sample_mapping$Sample_ID)]
colnames(neuron_counts) <- new_colnames

# Use biomaRt to get gene symbols
ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
gene_info <- getBM(attributes=c("ensembl_gene_id", "external_gene_name"),
                   filters="ensembl_gene_id",
                   values=rownames(neuron_counts),
                   mart=ensembl)

# Merge gene symbols with the data
neuron_data <- neuron_counts %>%
  rownames_to_column("ensembl_id") %>%
  left_join(gene_info, by = c("ensembl_id" = "ensembl_gene_id"))

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

# Identify RP genes
neuron_data$is_rp <- sapply(neuron_data$external_gene_name, is_rp_gene)

# Filter for RP genes and RT samples
rp_data <- neuron_data %>%
  filter(is_rp) %>%
  dplyr::select(external_gene_name, contains("_RT_"))

# Calculate total RP counts for each sample
total_rp_counts <- colSums(rp_data[, -1])

# Calculate fractions for each RP gene
rp_fractions <- sweep(rp_data[, -1], 2, total_rp_counts, "/")
rp_fractions$Gene <- rp_data$external_gene_name

# Reshape data for plotting
rp_fractions_long <- rp_fractions %>%
  pivot_longer(cols = -Gene, names_to = "Sample", values_to = "Fraction") %>%
  mutate(
    Cell_Type = case_when(
      grepl("vGluT2", Sample) ~ "vGluT2",
      grepl("Gad2", Sample) ~ "Gad2",
      grepl("PV", Sample) ~ "PV",
      grepl("SST", Sample) ~ "SST"
    ),
    Region = ifelse(grepl("_Cb_", Sample), "Cerebellum", "Cerebrum"),
    Disease = case_when(
      grepl("HD_", Sample) ~ "HD",
      grepl("CJD_", Sample) ~ "CJD",
      grepl("FFI_", Sample) ~ "FFI",
      grepl("dwt_", Sample) ~ "DWT",
      grepl("wt_", Sample) ~ "WT"
    )
  )

# Function to perform comparisons
perform_comparisons <- function(data, region) {
  if (region == "Cerebellum") {
    # Compare vGluT2 vs Gad2 in Cerebellum
    data <- data %>% filter(Cell_Type %in% c("vGluT2", "Gad2"))
    if (nrow(data) < 2 || length(unique(data$Cell_Type)) < 2) {
      return(data.frame(Gene = unique(data$Gene), Comparison = "vGluT2 vs Gad2", p_value = NA, mean_diff = NA))
    }
    t_test_result <- tryCatch({
      t.test(Fraction ~ Cell_Type, data = data)
    }, error = function(e) NULL)
    if (is.null(t_test_result)) {
      return(data.frame(Gene = unique(data$Gene), Comparison = "vGluT2 vs Gad2", p_value = NA, mean_diff = NA))
    }
    return(data.frame(
      Gene = unique(data$Gene),
      Comparison = "vGluT2 vs Gad2",
      p_value = t_test_result$p.value,
      mean_diff = diff(t_test_result$estimate)
    ))
  } else if (region == "Cerebrum") {
    # Compare vGluT2 vs Gad2 and PV vs SST separately in Cerebrum
    vglut2_gad2 <- data %>% filter(Cell_Type %in% c("vGluT2", "Gad2"))
    pv_sst <- data %>% filter(Cell_Type %in% c("PV", "SST"))
    
    results <- list()
    
    if (nrow(vglut2_gad2) >= 2 && length(unique(vglut2_gad2$Cell_Type)) == 2) {
      t_test_result <- tryCatch({
        t.test(Fraction ~ Cell_Type, data = vglut2_gad2)
      }, error = function(e) NULL)
      if (!is.null(t_test_result)) {
        results$vglut2_gad2 <- data.frame(
          Gene = unique(data$Gene),
          Comparison = "vGluT2 vs Gad2",
          p_value = t_test_result$p.value,
          mean_diff = diff(t_test_result$estimate)
        )
      }
    }
    
    if (nrow(pv_sst) >= 2 && length(unique(pv_sst$Cell_Type)) == 2) {
      t_test_result <- tryCatch({
        t.test(Fraction ~ Cell_Type, data = pv_sst)
      }, error = function(e) NULL)
      if (!is.null(t_test_result)) {
        results$pv_sst <- data.frame(
          Gene = unique(data$Gene),
          Comparison = "PV vs SST",
          p_value = t_test_result$p.value,
          mean_diff = diff(t_test_result$estimate)
        )
      }
    }
    
    return(do.call(rbind, results))
  }
}

# Function to process data and create plots
process_and_plot <- function(data, condition) {
  # Calculate mean fractions
  mean_fractions <- data %>%
    group_by(Gene, Cell_Type, Region) %>%
    summarize(
      Mean_Fraction = mean(Fraction),
      SD_Fraction = sd(Fraction),
      .groups = "drop"
    )
  
  # Perform comparisons
  results <- data %>%
    group_by(Gene, Region) %>%
    do(perform_comparisons(data = ., region = .$Region[1])) %>%
    ungroup()
  
  # Apply Benjamini-Hochberg correction
  results <- results %>%
    group_by(Region, Comparison) %>%
    mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
    ungroup()
  
  # Define custom colors
  custom_colors <- c("vGluT2" = "#DE3163", "Gad2" = "#0072B2", "PV" = "#708090", "SST" = "#FFBF00")
  
  # Function to create a bar plot for a single region
  create_region_barplot <- function(data, region_name, y_axis_limit) {
    ggplot(data, 
           aes(x = reorder(Gene, Mean_Fraction), y = Mean_Fraction, fill = Cell_Type)) +
      geom_bar(stat = "identity", position = position_dodge()) +
      geom_errorbar(aes(ymin = Mean_Fraction - SD_Fraction, 
                        ymax = Mean_Fraction + SD_Fraction),
                    position = position_dodge(0.9), width = 0.25) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
            legend.position = "top") +
      labs(title = paste(region_name, "Region -", condition),
           x = "Gene", y = "Mean Fraction", fill = "Cell Type") +
      scale_fill_manual(values = custom_colors) +
      ylim(0, y_axis_limit) +
      geom_text(aes(label = ifelse(p_adj < 0.05, "*", ""),
                    y = Mean_Fraction + SD_Fraction),
                position = position_dodge(0.9),
                vjust = -0.5)
  }
  
  # Create bar plots for each region
  regions <- c("Cerebellum", "Cerebrum")
  plots <- list()
  
  for (region in regions) {
    region_data <- mean_fractions %>%
      filter(Region == region) %>%
      left_join(results, by = c("Gene", "Region"))
    
    y_axis_limit <- max(region_data$Mean_Fraction + region_data$SD_Fraction, na.rm = TRUE) * 1.1
    
    p <- create_region_barplot(region_data, region, y_axis_limit)
    plots[[region]] <- p
    ggsave(paste0("rp_gene_fractions_", tolower(region), "_", tolower(condition), ".png"), p, width = 15, height = 8)
  }
  
  # Combine plots vertically
  combined_plot <- do.call(gridExtra::grid.arrange, c(plots, ncol = 1))
  ggsave(paste0("rp_gene_fractions_combined_vertical_", tolower(condition), ".png"), combined_plot, width = 15, height = 16)
  
  # Function to create a boxplot for a single gene
  create_gene_boxplot <- function(plot_data, gene, region, comparison) {
    gene_data <- results %>%
      filter(Gene == !!gene, Region == !!region, Comparison == !!comparison)
    
    p_value <- gene_data$p_value[1]
    
    cell_types <- strsplit(comparison, " vs ")[[1]]
    plot_data_filtered <- plot_data %>% 
      filter(Gene == !!gene, Region == !!region, Cell_Type %in% cell_types)
    
    ggboxplot(plot_data_filtered, 
              x = "Cell_Type", y = "Fraction", fill = "Cell_Type",
              palette = custom_colors,
              add = "jitter") +
      stat_compare_means(label.y = max(plot_data_filtered$Fraction, na.rm = TRUE) * 1.1) +
      labs(title = paste(gene, "-", condition, "\n", comparison, "\np-value =", format(p_value, digits = 3)),
           x = "Cell Type", y = "Fraction") +
      theme_minimal() +
      theme(legend.position = "none")
  }
  
  # Create boxplots for top differentially expressed genes in each region and comparison
  p_threshold <- 0.05
  for (r in regions) {
    if (r == "Cerebellum") {
      comparisons <- c("vGluT2 vs Gad2")
    } else {
      comparisons <- c("vGluT2 vs Gad2", "PV vs SST")
    }
    
    for (comp in comparisons) {
      region_sig_genes <- results %>%
        filter(Region == r, Comparison == comp, p_adj < p_threshold) %>%
        top_n(9, wt = -p_value) %>%
        pull(Gene)
      
      if (length(region_sig_genes) > 0) {
        region_plots <- lapply(region_sig_genes, function(g) create_gene_boxplot(data, g, r, comp))
        region_combined <- gridExtra::grid.arrange(
          grobs = region_plots, 
          ncol = 3,
          top = textGrob(paste(r, "Region -", comp, "- Top Differential Genes -", condition), gp = gpar(fontsize = 16, fontface = "bold"))
        )
        ggsave(paste0("differential_rps_boxplot_", tolower(r), "_", gsub(" vs ", "_vs_", tolower(comp)), "_", tolower(condition), ".png"), region_combined, width = 15, height = 15)
      }
    }
  }
  
  # Return results for saving
  return(list(results = results, mean_fractions = mean_fractions))
}

# Process and plot WT data
wt_data <- rp_fractions_long %>% filter(Disease == "WT")
wt_results <- process_and_plot(wt_data, "WT")

# Process and plot DWT data
dwt_data <- rp_fractions_long %>% filter(Disease == "DWT")
dwt_results <- process_and_plot(dwt_data, "DWT")

# Process and plot each disease separately
diseases <- c("HD", "CJD", "FFI")
disease_results <- list()

for (disease in diseases) {
  disease_data <- rp_fractions_long %>% filter(Disease == disease)
  disease_results[[disease]] <- process_and_plot(disease_data, disease)
}

# Save results
results_list <- list(
  "WT_Results" = wt_results$results,
  "WT_Mean_Fractions" = wt_results$mean_fractions,
  "DWT_Results" = dwt_results$results,
  "DWT_Mean_Fractions" = dwt_results$mean_fractions
)

for (disease in diseases) {
  results_list[[paste0(disease, "_Results")]] <- disease_results[[disease]]$results
  results_list[[paste0(disease, "_Mean_Fractions")]] <- disease_results[[disease]]$mean_fractions
}

write_xlsx(results_list, "neuron_rp_enrichment_comparison_results_all_conditions.xlsx")

print("Analysis complete. Check the generated PNG files for visualizations and the Excel file for detailed results.")



# Load required libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(biomaRt)

# Load neuron count data and sample information
neuron_counts <- read.csv("SB_merged_gene_counts_all.csv", row.names = 1)
sample_info <- read.csv("SB_sample_Info.csv")

# Create a mapping between Sample_ID and Sample_Name
sample_mapping <- sample_info %>%
  group_by(Sample_Name) %>%
  mutate(replicate = row_number()) %>%
  ungroup() %>%
  mutate(new_column_name = paste(Sample_Name, replicate, sep = "_rep")) %>%
  dplyr::select(Sample_ID, new_column_name)

# Rename the columns in neuron_counts
new_colnames <- sample_mapping$new_column_name[match(colnames(neuron_counts), sample_mapping$Sample_ID)]
colnames(neuron_counts) <- new_colnames

# Use biomaRt to get gene symbols
ensembl <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
gene_info <- getBM(attributes=c("ensembl_gene_id", "external_gene_name"),
                   filters="ensembl_gene_id",
                   values=rownames(neuron_counts),
                   mart=ensembl)

# Merge gene symbols with the data
neuron_data <- neuron_counts %>%
  rownames_to_column("ensembl_id") %>%
  left_join(gene_info, by = c("ensembl_id" = "ensembl_gene_id"))

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

# Identify RP genes
neuron_data$is_rp <- sapply(neuron_data$external_gene_name, is_rp_gene)

# Filter for RP genes and RT samples
rp_data <- neuron_data %>%
  filter(is_rp) %>%
  dplyr::select(external_gene_name, contains("_RT_"))

# Calculate total RP counts for each sample
total_rp_counts <- colSums(rp_data[, -1])

# Calculate fractions for each RP gene
rp_fractions <- sweep(rp_data[, -1], 2, total_rp_counts, "/")
rp_fractions$Gene <- rp_data$external_gene_name


# Modify the data preparation step to use 'Disease' instead of 'WT_Type'
rp_fractions_long <- rp_fractions %>%
  pivot_longer(cols = -Gene, names_to = "Sample", values_to = "Fraction") %>%
  mutate(
    Cell_Type = case_when(
      grepl("vGluT2", Sample) ~ "vGluT2",
      grepl("Gad2", Sample) ~ "Gad2",
      grepl("PV", Sample) ~ "PV",
      grepl("SST", Sample) ~ "SST"
    ),
    Region = ifelse(grepl("_Cb_", Sample), "Cerebellum", "Cerebrum"),
    Disease = case_when(
      grepl("HD_", Sample) ~ "HD",
      grepl("CJD_", Sample) ~ "CJD",
      grepl("FFI_", Sample) ~ "FFI",
      grepl("dwt_", Sample) ~ "WT",
      grepl("wt_", Sample) ~ "WT"
    )
  )
# Load required libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(ggpubr)

# ... [Previous data loading and processing code remains unchanged] ...

# Function to perform Wilcoxon test safely
safe_wilcox <- function(x, y) {
  if (length(x) < 2 || length(y) < 2) {
    return(NA)
  }
  tryCatch(
    wilcox.test(x, y)$p.value,
    error = function(e) NA
  )
}

# Select significant RPs
significant_rps <- wt_data %>%
  group_by(Gene, Region) %>%
  summarize(
    p_value_cerebellum = safe_wilcox(Fraction[Cell_Type == "vGluT2"], Fraction[Cell_Type == "Gad2"]),
    p_value_cerebrum_vg = safe_wilcox(Fraction[Cell_Type == "vGluT2"], Fraction[Cell_Type == "Gad2"]),
    p_value_cerebrum_ps = safe_wilcox(Fraction[Cell_Type == "PV"], Fraction[Cell_Type == "SST"]),
    .groups = "drop"
  ) %>%
  filter(!is.na(p_value_cerebellum) | !is.na(p_value_cerebrum_vg) | !is.na(p_value_cerebrum_ps)) %>%
  filter(p_value_cerebellum < 0.05 | p_value_cerebrum_vg < 0.05 | p_value_cerebrum_ps < 0.05) %>%
  arrange(p_value_cerebellum, p_value_cerebrum_vg, p_value_cerebrum_ps) %>%
  slice_head(n = 4) %>%
  pull(Gene)
# ... [Previous code for selecting significant RPs remains unchanged] ...

# Function to create boxplots
create_rp_boxplot <- function(data, gene, region, cell_types) {
  plot_data <- data %>% 
    filter(Gene == gene, Region == region, Cell_Type %in% cell_types) %>%
    mutate(Cell_Type = factor(Cell_Type, levels = cell_types))
  
  p <- ggplot(plot_data, aes(x = Cell_Type, y = Fraction, fill = Cell_Type)) +
    geom_boxplot() +
    geom_jitter(width = 0.2, alpha = 0.5) +
    scale_fill_manual(values = c("vGluT2" = "#DE3163", "Gad2" = "#0072B2", "PV" = "#708090", "SST" = "#FFBF00")) +
    theme_minimal() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  # Increased size
          axis.text.y = element_text(size = 12),  # Increased size
          plot.title = element_text(size = 14, face = "bold"),  # Increased size
          axis.title.y = element_text(size = 14)) +  # Increased size
    labs(title = gene,  # Removed region from title
         x = NULL, y = "Fraction") +
    ylim(0, NA)
  
  if (length(cell_types) == 2) {
    p <- p + stat_compare_means(comparisons = list(cell_types), 
                                method = "wilcox.test",
                                label = "p.signif")
  }
  
  return(p)
}

# Create plots for each significant gene and region
plots <- list()
for (gene in significant_rps) {
  plots[[paste(gene, "Cerebellum")]] <- create_rp_boxplot(wt_data, gene, "Cerebellum", c("vGluT2", "Gad2"))
  plots[[paste(gene, "Cerebrum_VG")]] <- create_rp_boxplot(wt_data, gene, "Cerebrum", c("vGluT2", "Gad2"))
  plots[[paste(gene, "Cerebrum_PS")]] <- create_rp_boxplot(wt_data, gene, "Cerebrum", c("PV", "SST"))
}

# Arrange plots in a 3x4 grid
arranged_plots <- grid.arrange(
  plots[[paste(significant_rps[1], "Cerebellum")]], plots[[paste(significant_rps[2], "Cerebellum")]],
  plots[[paste(significant_rps[3], "Cerebellum")]], plots[[paste(significant_rps[4], "Cerebellum")]],
  plots[[paste(significant_rps[1], "Cerebrum_VG")]], plots[[paste(significant_rps[2], "Cerebrum_VG")]],
  plots[[paste(significant_rps[3], "Cerebrum_VG")]], plots[[paste(significant_rps[4], "Cerebrum_VG")]],
  plots[[paste(significant_rps[1], "Cerebrum_PS")]], plots[[paste(significant_rps[2], "Cerebrum_PS")]],
  plots[[paste(significant_rps[3], "Cerebrum_PS")]], plots[[paste(significant_rps[4], "Cerebrum_PS")]],
  ncol = 4, nrow = 3
)

# Create a key for the regions and significance levels
key_text <- c(
  "Top row: Cerebellum (vGluT2 vs Gad2)",
  "Middle row: Cerebrum (vGluT2 vs Gad2)",
  "Bottom row: Cerebrum (PV vs SST)",
  "Significance levels: **** p<0.0001, *** p<0.001, ** p<0.01, * p<0.05, ns: not significant",
  "Wilcoxon test used for all comparisons"
)
key <- textGrob(paste(key_text, collapse = "\n"), 
                gp = gpar(fontface = "bold", fontsize = 12))  # Increased size

# Combine the plots with the key
final_plot <- grid.arrange(arranged_plots, key, 
                           heights = c(20, 2.5),  # Slightly increased height for key
                           nrow = 2)

# Save the final plot
ggsave("significant_rp_boxplots.png", final_plot, width = 16, height = 15, dpi = 300)

print("Analysis complete. Check the generated 'significant_rp_boxplots.png' file for the updated boxplots of significant RPs.")



# ... [Previous code remains unchanged until the significant_rps selection] ...

# Manually select RPs to plot
selected_rps <- c("Rpl3l", "Rps4l", "Rps21", "Rps25")

# Function to create boxplots
create_rp_boxplot <- function(data, gene, region, cell_types) {
  plot_data <- data %>% 
    filter(Gene == gene, Region == region, Cell_Type %in% cell_types) %>%
    mutate(Cell_Type = factor(Cell_Type, levels = cell_types))
  
  p <- ggplot(plot_data, aes(x = Cell_Type, y = Fraction, fill = Cell_Type)) +
    geom_boxplot() +
    geom_jitter(width = 0.2, alpha = 0.5) +
    scale_fill_manual(values = c("vGluT2" = "#DE3163", "Gad2" = "#0072B2", "PV" = "#708090", "SST" = "#FFBF00")) +
    theme_minimal() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(size = 10, face = "bold")) +
    labs(title = paste(gene, "-", region),
         x = NULL, y = "Fraction") +
    ylim(0, NA)
  
  if (length(cell_types) == 2) {
    p <- p + stat_compare_means(comparisons = list(cell_types), 
                                method = "wilcox.test",
                                label = "p.signif")
  }
  
  return(p)
}

# Create plots for each selected gene and region
plots <- list()
for (gene in selected_rps) {
  plots[[paste(gene, "Cerebellum")]] <- create_rp_boxplot(wt_data, gene, "Cerebellum", c("vGluT2", "Gad2"))
  plots[[paste(gene, "Cerebrum_VG")]] <- create_rp_boxplot(wt_data, gene, "Cerebrum", c("vGluT2", "Gad2"))
  plots[[paste(gene, "Cerebrum_PS")]] <- create_rp_boxplot(wt_data, gene, "Cerebrum", c("PV", "SST"))
}

# Arrange plots in a 3x4 grid
arranged_plots <- grid.arrange(
  plots[[paste(selected_rps[1], "Cerebellum")]], plots[[paste(selected_rps[2], "Cerebellum")]],
  plots[[paste(selected_rps[3], "Cerebellum")]], plots[[paste(selected_rps[4], "Cerebellum")]],
  plots[[paste(selected_rps[1], "Cerebrum_VG")]], plots[[paste(selected_rps[2], "Cerebrum_VG")]],
  plots[[paste(selected_rps[3], "Cerebrum_VG")]], plots[[paste(selected_rps[4], "Cerebrum_VG")]],
  plots[[paste(selected_rps[1], "Cerebrum_PS")]], plots[[paste(selected_rps[2], "Cerebrum_PS")]],
  plots[[paste(selected_rps[3], "Cerebrum_PS")]], plots[[paste(selected_rps[4], "Cerebrum_PS")]],
  ncol = 4, nrow = 3
)

# Create a key for the regions and significance levels
key_text <- c(
  "Top row: Cerebellum (vGluT2 vs Gad2)",
  "Middle row: Cerebrum (vGluT2 vs Gad2)",
  "Bottom row: Cerebrum (PV vs SST)",
  "Significance levels: **** p<0.0001, *** p<0.001, ** p<0.01, * p<0.05, ns: not significant",
  "Wilcoxon test used for all comparisons"
)
key <- textGrob(paste(key_text, collapse = "\n"), 
                gp = gpar(fontface = "bold", fontsize = 10))

# Combine the plots with the key
final_plot <- grid.arrange(arranged_plots, key, 
                           heights = c(20, 2),
                           nrow = 2)

# Save the final plot
ggsave("selected_rp_boxplots.png", final_plot, width = 16, height = 15, dpi = 300)

print("Analysis complete. Check the generated 'selected_rp_boxplots.png' file for the boxplots of selected RPs.")










