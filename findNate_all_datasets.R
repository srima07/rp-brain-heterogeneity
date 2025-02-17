# Load required libraries
library(tidyverse)
library(ggplot2)
library(patchwork)

# Function to process data
process_data <- function(file_path, dataset_name) {
  data <- read.csv(file_path, stringsAsFactors = FALSE)
  
  data_normalized <- data %>%
    mutate(
      A_norm = A_ms_205_KKK / Exn4,
      B_norm = B_ms_202_KK / Exn4,
      C_norm = C_ms_210_PKE / Exn4,
      Total_norm = A_norm + B_norm + C_norm,
      A_percent = A_norm / Total_norm * 100,
      B_percent = B_norm / Total_norm * 100,
      C_percent = C_norm / Total_norm * 100
    )
  
  if (dataset_name == "neuron") {
    data_normalized <- data_normalized %>%
      mutate(
        Disease = case_when(
          grepl("^HD", Sample) ~ "HD",
          grepl("^FFI", Sample) ~ "FFI",
          grepl("^CJD", Sample) ~ "CJD",
          grepl("^Wt|^dWt", Sample) ~ "Wt",
          TRUE ~ "Unknown"
        ),
        CellType = case_when(
          grepl("vGluT2", Sample) ~ "vGluT2",
          grepl("Gad2", Sample) ~ "Gad2",
          grepl("PV", Sample) ~ "PV",
          grepl("SST", Sample) ~ "SST",
          TRUE ~ "Unknown"
        ),
        Region = case_when(
          grepl("Cer_", Sample) ~ "Cerebrum",
          grepl("Cb", Sample) ~ "Cerebellum",
          TRUE ~ "Unknown"
        ),
        RNAType = ifelse(grepl("_TR_", Sample), "TR", "RT"),
        Group = case_when(
          RNAType == "TR" & Region == "Cerebrum" ~ "TR_Cer",
          RNAType == "TR" & Region == "Cerebellum" ~ "TR_Cb",
          TRUE ~ paste(CellType, Region, RNAType, sep = "_")
        )
      ) %>%
      filter(Disease == "Wt")  # Filter to include only wild-type samples
  } else if (dataset_name == "astromicro") {
    data_normalized <- data_normalized %>%
      mutate(
        CellType = ifelse(grepl("^AC", Sample), "AC", "MG"),
        Region = case_when(
          grepl("_F_", Sample) ~ "F",
          grepl("_M_", Sample) ~ "M",
          grepl("_R_", Sample) ~ "R"
        ),
        RNAType = ifelse(grepl("_RT_", Sample), "RT", "TR"),
        Group = paste(CellType, Region, RNAType, sep = "_")
      )
  } else if (dataset_name == "kang") {
    data_normalized <- data_normalized %>%
      mutate(
        Group = Sample
      )
  }
  
  data_normalized$Dataset <- dataset_name
  return(data_normalized)
}

# Process all datasets
neuron_data <- process_data("SB_neurons_RPs/No_zero_All_SB_isoform_counts.csv", "neuron")
astromicro_data <- process_data("astromicro_isoform_counts.csv", "astromicro")
kang_data <- process_data("kang_isoform_counts.csv", "kang")

# Combine all datasets
all_data <- bind_rows(neuron_data, astromicro_data, kang_data)

# Function to create summary statistics
create_summary <- function(data) {
  data %>%
    group_by(Dataset, Group) %>%
    summarise(across(c(A_percent, B_percent, C_percent), list(mean = mean, se = ~sd(.)/sqrt(n()))), .groups = "drop")
}

# Create summary for all datasets
all_summary <- create_summary(all_data)



# Function to order groups for each dataset
order_groups <- function(data, dataset) {
  if (dataset == "astromicro") {
    order <- c(
      "AC_F_TR", "AC_M_TR", "AC_R_TR",
      "MG_F_TR", "MG_M_TR", "MG_R_TR",
      "AC_F_RT", "AC_M_RT", "AC_R_RT",
      "MG_F_RT", "MG_M_RT", "MG_R_RT"
    )
  } else if (dataset == "neuron") {
    order <- c(
      "TR_Cer", "vGluT2_Cerebrum_RT", "Gad2_Cerebrum_RT", "PV_Cerebrum_RT", "SST_Cerebrum_RT",
      "TR_Cb", "vGluT2_Cerebellum_RT", "Gad2_Cerebellum_RT"
    )
  } else if (dataset == "kang") {
    order <- c(
      "TR", "RT_3mon_F", "RT_3mon_M", "RT_12mon_F", "RT_12mon_M", "RT_24mon_F", "RT_24mon_M",
      "App_wt", "App_Tg", "RT_GFP", "RT_Tau", "RT_PBS", "RT_IC", "RT_LPS", "CS_PBS", "CS_LPS"
    )
  }
  
  data$Group <- factor(data$Group, levels = order)
  return(data)
}

# Apply ordering to each dataset
astromicro_summary <- filter(all_summary, Dataset == "astromicro") %>% order_groups("astromicro")
neuron_summary <- filter(all_summary, Dataset == "neuron") %>% order_groups("neuron")
kang_summary <- filter(all_summary, Dataset == "kang") %>% order_groups("kang")


# Function to create overlayed (stacked) plot with consistent bar width
create_overlayed_plot <- function(data, title, x_var) {
  data_long <- data %>%
    pivot_longer(cols = ends_with("_mean"), 
                 names_to = "Isoform", 
                 values_to = "Percentage") %>%
    mutate(Isoform = str_replace(Isoform, "_percent_mean", ""))
  
  # Set a fixed bar width
  fixed_bar_width <- 0.8
  
  # Calculate the number of groups
  n_groups <- length(unique(data_long[[x_var]]))
  
  # Calculate plot width based on number of groups and fixed bar width
  plot_width <- n_groups * fixed_bar_width * 1.5  # 1.5 to account for spaces and margins
  
  p <- ggplot(data_long, aes(x = !!sym(x_var), y = Percentage, fill = Isoform)) +
    geom_bar(stat = "identity", position = "stack", width = fixed_bar_width) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_text(size = 8),
      axis.title = element_text(size = 10),
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 10),
      plot.title = element_text(size = 12),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    ) +
    labs(title = title, x = x_var, y = "Isoform Percentage") +
    scale_fill_manual(
      values = c("A" = "#006400", "B" = "#E3735E", "C" = "#800080"),
      labels = c(expression(italic("Rps24a")), expression(italic("Rps24b")), expression(italic("Rps24c")))
    ) +
    coord_cartesian(ylim = c(0, 100)) +
    scale_x_discrete(expand = expansion(mult = c(0.05, 0.2)))
  
  return(list(plot = p, width = plot_width))
}


# Function to create separated plot
create_separated_plot <- function(data, title, x_var) {
  data_long <- data %>%
    pivot_longer(cols = ends_with("_mean"), 
                 names_to = "Isoform", 
                 values_to = "Percentage") %>%
    mutate(Isoform = str_replace(Isoform, "_percent_mean", ""))
  
  data_long <- data_long %>%
    mutate(SE = case_when(
      Isoform == "A" ~ A_percent_se,
      Isoform == "B" ~ B_percent_se,
      Isoform == "C" ~ C_percent_se
    ))
  
  # Calculate the number of groups
  n_groups <- length(unique(data_long[[x_var]]))
  
  # Calculate plot width based on number of groups
  plot_width <- max(12, n_groups * 0.8)
  
  p <- ggplot(data_long, aes(x = !!sym(x_var), y = Percentage, fill = Isoform)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_errorbar(aes(ymin = Percentage - SE, ymax = Percentage + SE), 
                  width = 0.35) +
    facet_wrap(~Isoform, ncol = 1, scales = "free_y") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_text(size = 8),
      axis.title = element_text(size = 10),
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 10),
      plot.title = element_text(size = 12),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      strip.text = element_text(size = 10, face = "bold")
    ) +
    labs(title = title, x = x_var, y = "Isoform Percentage") +
    scale_fill_manual(
      values = c("A" = "#006400", "B" = "#E3735E", "C" = "#800080"),
      labels = c(expression(italic("Rps24a")), expression(italic("Rps24b")), expression(italic("Rps24c")))
    ) +
    scale_x_discrete(expand = expansion(mult = c(0.05, 0.2)))
  
  return(list(plot = p, width = plot_width))
}

# Create plots for each dataset

astromicro_separated <- create_separated_plot(astromicro_summary, "Astromicro Dataset - Separated", "Group")

neuron_separated <- create_separated_plot(neuron_summary, "Neuron Dataset (Wild-type) - Separated", "Group")

kang_separated <- create_separated_plot(kang_summary, "Kang Dataset - Separated", "Group")

# Create plots for each dataset
astromicro_overlayed <- create_overlayed_plot(astromicro_summary, "Astromicro Dataset - Overlayed", "Group")
neuron_overlayed <- create_overlayed_plot(neuron_summary, "Neuron Dataset (Wild-type) - Overlayed", "Group")
kang_overlayed <- create_overlayed_plot(kang_summary, "Kang Dataset - Overlayed", "Group")


# Save plots with dynamic widths

ggsave("astromicro_separated_ordered.png", astromicro_separated$plot, width = astromicro_separated$width, height = 12, dpi = 300)
ggsave("neuron_wt_separated_ordered.png", neuron_separated$plot, width = neuron_separated$width, height = 12, dpi = 300)
ggsave("kang_separated_ordered.png", kang_separated$plot, width = kang_separated$width, height = 12, dpi = 300)


# Save overlay plots with dynamic widths but consistent bar widths
ggsave("astromicro_overlayed_ordered.png", astromicro_overlayed$plot, width = astromicro_overlayed$width, height = 6, dpi = 300, limitsize = FALSE)
ggsave("neuron_wt_overlayed_ordered.png", neuron_overlayed$plot, width = neuron_overlayed$width, height = 6, dpi = 300, limitsize = FALSE)
ggsave("kang_overlayed_ordered.png", kang_overlayed$plot, width = kang_overlayed$width, height = 6, dpi = 300, limitsize = FALSE)




# Print summary of data
cat("Summary of data:\n")
print(all_summary %>% 
        group_by(Dataset) %>% 
        summarise(n = n(), .groups = "drop"))

# Save summary statistics
write.csv(all_summary, "all_datasets_summary_statistics_ordered.csv", row.names = FALSE)



# Load required libraries
library(tidyverse)
library(ggplot2)
library(patchwork)
library(rstatix)

# ... [Previous data loading and processing code remains the same] ...

# Function to perform statistical tests
perform_stats <- function(data, group1, group2, dataset) {
  isoforms <- c("A_percent", "B_percent", "C_percent")
  results <- data.frame()
  
  for (isoform in isoforms) {
    test_result <- data %>%
      filter(Group %in% c(group1, group2)) %>%
      t_test(as.formula(paste(isoform, "~ Group")), 
             paired = FALSE, 
             var.equal = FALSE) %>%
      add_column(Dataset = dataset, Isoform = isoform, .before = 1)
    
    results <- bind_rows(results, test_result)
  }
  
  return(results)
}

# Perform statistical tests for each dataset
stats_results <- data.frame()

# Glia data comparisons
glia_comparisons <- list(
  c("AC_F_RT", "MG_F_RT"),
  c("AC_M_RT", "MG_M_RT"),
  c("AC_R_RT", "MG_R_RT"),
  c("AC_F_TR", "AC_F_RT"),
  c("AC_M_TR", "AC_M_RT"),
  c("AC_R_TR", "AC_R_RT"),
  c("MG_F_TR", "MG_F_RT"),
  c("MG_M_TR", "MG_M_RT"),
  c("MG_R_TR", "MG_R_RT")
)

for (comp in glia_comparisons) {
  stats_results <- bind_rows(stats_results, 
                             perform_stats(astromicro_data, comp[1], comp[2], "Glia"))
}


# Neuron data comparisons (WT only)
neuron_comparisons <- list(
  # Original TR vs RT comparisons
  c("TR_Cer", "vGluT2_Cerebrum_RT"),
  c("TR_Cer", "Gad2_Cerebrum_RT"),
  c("TR_Cer", "PV_Cerebrum_RT"),
  c("TR_Cer", "SST_Cerebrum_RT"),
  c("TR_Cb", "vGluT2_Cerebellum_RT"),
  c("TR_Cb", "Gad2_Cerebellum_RT"),
  
  # New RT vs RT comparisons
  c("vGluT2_Cerebrum_RT", "Gad2_Cerebrum_RT"),
  c("vGluT2_Cerebrum_RT", "PV_Cerebrum_RT"),
  c("vGluT2_Cerebrum_RT", "SST_Cerebrum_RT"),
  c("Gad2_Cerebrum_RT", "PV_Cerebrum_RT"),
  c("Gad2_Cerebrum_RT", "SST_Cerebrum_RT"),
  c("PV_Cerebrum_RT", "SST_Cerebrum_RT"),
  c("vGluT2_Cerebellum_RT", "Gad2_Cerebellum_RT")
)


for (comp in neuron_comparisons) {
  stats_results <- bind_rows(stats_results, 
                             perform_stats(neuron_data, comp[1], comp[2], "Neuron"))
}

# Kang data comparisons
kang_comparisons <- list(
  c("RT_3mon_F", "RT_12mon_F"),
  c("RT_3mon_F", "RT_24mon_F"),
  c("RT_3mon_M", "RT_12mon_M"),
  c("RT_3mon_M", "RT_24mon_M"),
  c("App_wt", "App_Tg"),
  c("RT_GFP", "RT_Tau"),
  c("RT_PBS", "RT_IC"),
  c("RT_PBS", "RT_LPS"),
  c("CS_PBS", "CS_LPS")
)

for (comp in kang_comparisons) {
  stats_results <- bind_rows(stats_results, 
                             perform_stats(kang_data, comp[1], comp[2], "Kang"))
}

# Adjust p-values for multiple comparisons
stats_results <- stats_results %>%
  group_by(Dataset, Isoform) %>%
  mutate(p.adjusted = p.adjust(p, method = "BH")) %>%
  ungroup()

# Save statistical results
write.csv(stats_results, "isoform_statistical_analysis.csv", row.names = FALSE)

# Print summary of statistical results
cat("Summary of statistical tests:\n")
print(stats_results %>%
        group_by(Dataset) %>%
        summarise(n_tests = n(),
                  n_significant = sum(p.adjusted < 0.05),
                  .groups = "drop"))

# ... [Previous plotting code remains the same] ...

# Add significance annotations to plots
add_significance <- function(plot, stats, dataset) {
  sig_data <- stats %>%
    filter(Dataset == dataset, p.adjusted < 0.05) %>%
    mutate(y_position = 105,  # Adjust this value to position the asterisks
           label = "*")
  
  plot + 
    geom_text(data = sig_data, 
              aes(x = group1, y = y_position, label = label),
              size = 5)
}

# Update plots with significance annotations
astromicro_overlayed$plot <- add_significance(astromicro_overlayed$plot, stats_results, "Glia")
neuron_overlayed$plot <- add_significance(neuron_overlayed$plot, stats_results, "Neuron")
kang_overlayed$plot <- add_significance(kang_overlayed$plot, stats_results, "Kang")

# Save updated plots
ggsave("astromicro_overlayed_ordered_with_stats.png", astromicro_overlayed$plot, width = astromicro_overlayed$width, height = 6, dpi = 300, limitsize = FALSE)
ggsave("neuron_wt_overlayed_ordered_with_stats.png", neuron_overlayed$plot, width = neuron_overlayed$width, height = 6, dpi = 300, limitsize = FALSE)
ggsave("kang_overlayed_ordered_with_stats.png", kang_overlayed$plot, width = kang_overlayed$width, height = 6, dpi = 300, limitsize = FALSE)





