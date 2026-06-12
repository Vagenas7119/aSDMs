
#June 2026

#Development pipeline of freshwater SDMs (aSDMS)

#Post-Analysis Plots

#Sector: Final pipeline

#PhD Researcher - Georgios Vagenas (georgios.vagenas@mncn.csic.es | georgvagenas@gmail.com)



#diagnostic script to check the identical structure of all the folders in my repository


# =======================================================================
# Diagnostic Script: Check species folder structure in "regional/"
# =======================================================================

# 1. Define paths and expected structure
base_dir <- "C:/Users/geo_v/Desktop/Vagenas_aSDMs/output/regional/"
# Level 2 folders (inside each species folder)
expected_l2_folders <- c("eco", "h5", "h8", "h12")

# Level 3 folders (inside eco, h5, h8, h12)
expected_l3_folders <- c("Climate", "Hydroclimatic", "Hydromorphological")

# Level 4 files (inside Climate, Hydroclimatic, Hydromorphology)
expected_files <- c(
  "ensemble_regional.tif", 
  "eval_regional_iberia_strict.csv", 
  "regional_SDM.rds", 
  "varImp_regional.rds"
)

# 2. Check if the base directory exists
if (!dir.exists(base_dir)) {
  stop("Error: The base directory 'regional/' does not exist in the current working directory.")
}

# Get all species directories inside 'regional/' (excluding files)
species_dirs <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)

if (length(species_dirs) == 0) {
  stop("Warning: 'regional/' is empty or contains no species subdirectories.")
}

# 3. Initialize a list to store diagnostic issues
issues <- list()

# 4. Iterate over the nested structure
for (sp_dir in species_dirs) {
  sp_name <- basename(sp_dir)
  
  # A. Check Level 2 (eco, h5, h8, h12)
  for (l2_folder in expected_l2_folders) {
    l2_path <- file.path(sp_dir, l2_folder)
    
    if (!dir.exists(l2_path)) {
      issues[[length(issues) + 1]] <- data.frame(
        Species = sp_name,
        L2_Folder = l2_folder,
        L3_Folder = NA,
        Missing_Item = "ENTIRE L2 FOLDER MISSING",
        Type = "Folder"
      )
      next # Skip deeper checks if the parent folder is missing
    }
    
    # B. Check Level 3 (Climate, Hydroclimatic, Hydromorphology)
    for (l3_folder in expected_l3_folders) {
      l3_path <- file.path(l2_path, l3_folder)
      
      if (!dir.exists(l3_path)) {
        issues[[length(issues) + 1]] <- data.frame(
          Species = sp_name,
          L2_Folder = l2_folder,
          L3_Folder = l3_folder,
          Missing_Item = "ENTIRE L3 FOLDER MISSING",
          Type = "Folder"
        )
        next # Skip file checks if this intermediate folder is missing
      }
      
      # C. Check Level 4 (The 4 final files)
      for (file_name in expected_files) {
        file_path <- file.path(l3_path, file_name)
        
        if (!file.exists(file_path)) {
          issues[[length(issues) + 1]] <- data.frame(
            Species = sp_name,
            L2_Folder = l2_folder,
            L3_Folder = l3_folder,
            Missing_Item = file_name,
            Type = "File"
          )
        }
      }
    }
  }
}

# 5. Compile and report results
if (length(issues) == 0) {
  message("✅ SUCCESS: All species folders have the exact nested structure and contain all expected files!")
} else {
  # Combine the list of issues into a single data frame
  issues_df <- do.call(rbind, issues)
  
  message("⚠️ ISSUES FOUND: The following nested folders/files are missing:")
  print(issues_df)
  
  # Save the diagnostic report
  write.csv(issues_df, "structural_diagnostic_report_nested.csv", row.names = FALSE)
  message("\n-> A detailed report has been saved as 'structural_diagnostic_report_nested.csv' in your working directory.")
}



#### Pre-setting :: Libraries required to perform the analysis ####
library(dplyr)
library(readr)
library(tidyr)
library(terra)
library(raster)
library(ggplot2)
library(ggpubr)    # For violin plots & auto-statistical significance
library(sf)        # For loading and cropping shapefiles
library(patchwork) # For stitching maps and plots side-by-side
library(ggh4x)
library(rnaturalearth)
library(rnaturalearthdata)
library(rstatix) # REQUIRED for wilcox_test() and add_significance()


setwd('/Users/geo_v/Desktop/')

figdir<-"/Users/geo_v/Desktop/Vagenas_aSDMs/output/figures"
tablesdir<-"/Users/geo_v/Desktop/Vagenas_aSDMs/output/tables"

if (!dir.exists(figdir)) {
  dir.create(figdir, recursive = TRUE)
}


if (!dir.exists(tablesdir)) {
  dir.create(tablesdir, recursive = TRUE)
}


# ==============================================================================
#### GLOBAL aSDMs (N=51) | COMPARISONS PERFORMANCE METRICS - SUPPLEMENTARY MATERIAL ####
# ==============================================================================


# 1. Get ALL global evaluation files recursively
all_eval_files <- list.files(
  path = "Vagenas_aSDMs/output/global/",
  pattern = "eval_global_iberia_strict.csv",
  recursive = TRUE,
  full.names = TRUE
)

# 2. Combine all files into one master dataframe
all_metrics <- lapply(all_eval_files, function(f) {
  df <- read_csv(f, show_col_types = FALSE)
  parts <- unlist(strsplit(f, "/"))
  df$Species <- parts[length(parts) - 3]
  df$Extent  <- parts[length(parts) - 2]
  df$Set     <- parts[length(parts) - 1]
  return(df)
}) %>% bind_rows()

# 3. Pivot to compare Climate vs Hydroclimatic
# We filter only for the two sets we want to compare
comparison_table <- all_metrics %>%
  dplyr::filter(Set %in% c("Climate", "Hydroclimatic","Hydromorphological")) %>%
  dplyr::select(Species, Extent, Set, e_AUC, CBI, maxTSS, uAUC) %>%
  pivot_wider(
    names_from = Set, 
    values_from = c(e_AUC, CBI, maxTSS, uAUC)
  ) %>%
  mutate(
    dAUC = e_AUC_Hydroclimatic - e_AUC_Climate,
    dCBI = CBI_Hydroclimatic - CBI_Climate,
    dTSS = maxTSS_Hydroclimatic - maxTSS_Climate,
    duAUC = uAUC_Hydroclimatic - uAUC_Climate
  )

# 4. Statistical Tests (Looping through all 4 metrics)
metrics_to_test <- c("e_AUC", "CBI", "maxTSS", "uAUC")

cat("\n--- PAIRED WILCOXON SIGNED-RANK TEST RESULTS ---\n")
for (m in metrics_to_test) {
  hydro_col <- paste0(m, "_Hydroclimatic")
  clim_col <- paste0(m, "_Climate")
  
  # Filter out NAs for this specific pair
  data_subset <- comparison_table %>% filter(!is.na(!!sym(hydro_col)) & !is.na(!!sym(clim_col)))
  
  test <- wilcox.test(data_subset[[hydro_col]], data_subset[[clim_col]], paired = TRUE)
  
  cat(sprintf("\nMetric: %s\n", m))
  cat(sprintf("  p-value: %.5f %s\n", test$p.value, ifelse(test$p.value < 0.05, "*", "")))
  cat(sprintf("  Mean Hydro: %.4f | Mean Climate: %.4f\n", mean(data_subset[[hydro_col]]), mean(data_subset[[clim_col]])))
}

# 5. Visualizing the paired improvement

boxplot(comparison_table$duAUC, main="Distribution of Delta AUC (Hydro - Climate)", 
        ylab="Improvement in AUC", col="skyblue")
abline(h=0, col="red", lty=2)



# 1. Get the list of unique extents
extents_list <- unique(comparison_table$Extent)

cat("\n--- DETAILED STATISTICAL COMPARISON BY EXTENT ---\n")

for (ext in extents_list) {
  cat(sprintf("\n============================================\n"))
  cat(sprintf("EXTENT: %s\n", ext))
  cat(sprintf("============================================\n"))
  
  # Subset data for this extent
  ext_data <- comparison_table %>% filter(Extent == ext)
  
  # Run tests for each metric
  for (m in c("e_AUC", "CBI", "maxTSS", "uAUC")) {
    hydro_col <- paste0(m, "_Hydroclimatic")
    clim_col  <- paste0(m, "_Climate")
    
    # Filter for the specific columns and remove NAs
    subset <- ext_data %>% filter(!is.na(!!sym(hydro_col)) & !is.na(!!sym(clim_col)))
    
    if(nrow(subset) > 1) { # Need at least 2 species to run a paired test
      test <- wilcox.test(subset[[hydro_col]], subset[[clim_col]], paired = TRUE)
      
      cat(sprintf("  Metric: %s | p-value: %.5f %s\n", 
                  m, test$p.value, ifelse(test$p.value < 0.05, "*", "")))
    }
  }
}

library(ggplot2)

# Reshape data for plotting
plot_data <- comparison_table %>%
  dplyr::select(Species, Extent, e_AUC_Climate, e_AUC_Hydroclimatic) %>%
  pivot_longer(cols = c(e_AUC_Climate, e_AUC_Hydroclimatic), 
               names_to = "Model", values_to = "AUC")

# Create the faceted plot
ggplot(plot_data, aes(x = Model, y = AUC, fill = Model)) +
  geom_boxplot() +
  facet_wrap(~Extent) +
  theme_minimal() +
  labs(title = "AUC Performance by Extent: Climate vs Hydroclimatic",
       y = "AUC Score")

# ==============================================================================
# 1. DEFINE YOUR SPATIAL PATHS (Update the exact .shp filenames here!)
# ==============================================================================
iberia_path <- "Vagenas_aSDMs/input/SpatialExtents/study_area_iberia.shp"

# Add the specific .shp filenames to the end of these paths
extent_paths <- list(
  "eco" = "Vagenas_aSDMs/input/SpatialExtents/ecoregions/feow_hydrosheds.shp",
  "H5"  = "Vagenas_aSDMs/input/SpatialExtents/HydroSHEDS/H5/hybas_lake_eu_lev05_v1c/hybas_lake_eu_lev05_v1c.shp",
  "H8"  = "Vagenas_aSDMs/input/SpatialExtents/HydroSHEDS/H8/hybas_lake_eu_lev08_v1c/hybas_lake_eu_lev08_v1c.shp",
  "H12" = "Vagenas_aSDMs/input/SpatialExtents/HydroSHEDS/H12/hybas_lake_eu_lev12_v1c/hybas_lake_eu_lev12_v1c.shp"
)

# Load the base Iberian map
iberia_sf <- st_read(iberia_path, quiet = TRUE)

# 2. LOAD AND PREP DATA
all_eval_files <- list.files("Vagenas_aSDMs/output/global/", 
                             pattern = "eval_global_iberia_strict.csv", recursive = TRUE, full.names = TRUE)

all_metrics <- lapply(all_eval_files, function(f) {
  df <- read_csv(f, show_col_types = FALSE)
  parts <- unlist(strsplit(f, "/"))
  df$Species <- parts[length(parts) - 3]
  df$Extent  <- parts[length(parts) - 2]
  df$Set     <- parts[length(parts) - 1]
  return(df)
}) %>% bind_rows()

# 🚨 RENAME e_AUC to AUC here
plot_data <- all_metrics %>%
  filter(Set %in% c("Climate", "Hydroclimatic","Hydromorphological")) %>%
  dplyr::select(Species, Extent, Set, e_AUC, CBI, maxTSS, uAUC) %>%
  rename(AUC = e_AUC) %>% 
  pivot_longer(cols = c(AUC, CBI, maxTSS, uAUC), names_to = "Metric", values_to = "Score") %>%
  filter(!is.na(Score))

plot_data$Set <- factor(plot_data$Set, levels = c("Climate", "Hydroclimatic"))

# ==============================================================================
# 3. MASTER VISUALIZATION LOOP
# ==============================================================================


figdirglobal<-"/Users/geo_v/Desktop/Vagenas_aSDMs/output/figures/Global/"

figdirreg<-"/Users/geo_v/Desktop/Vagenas_aSDMs/output/figures/Regional/"

if (!dir.exists(figdirglobal)) {
  dir.create(figdirglobal, recursive = TRUE)
}


if (!dir.exists(figdirreg)) {
  dir.create(figdirreg, recursive = TRUE)
}

#mask it

# Fetch and isolate the Iberian Mainland using terra
iberia_sf <- ne_countries(country = c("Spain", "Portugal"), scale = "medium", returnclass = "sf")
iberia_vect <- vect(iberia_sf)
iberia_dissolved <- aggregate(iberia_vect)
iberia_parts <- disagg(iberia_dissolved)
land_areas <- expanse(iberia_parts, unit = "km")
iberia_mainland <- iberia_parts[which.max(land_areas), ]

# THE BRIDGE: Convert the terra mainland object back to an sf object
iberia_mainland_sf <- st_as_sf(iberia_mainland)

# Clean it once here so you don't have to repeat it inside the loop
iberia_clean <- st_zm(iberia_mainland_sf, drop = TRUE, what = "ZM")


for (ext in names(extent_paths)) {
  cat(sprintf("Processing Extent: %s\n", ext))
  
  # A. MAP PLOT (Bulletproof Topology Bypass)
  sf::sf_use_s2(FALSE)
  
  ext_sf <- st_read(extent_paths[[ext]], quiet = TRUE)
  ext_sf <- st_zm(ext_sf, drop = TRUE, what = "ZM")
  
  # Align CRS and fix topology
  st_crs(ext_sf) <- st_crs(iberia_clean)
  ext_sf <- suppressWarnings(st_buffer(st_make_valid(ext_sf), 0))
  
  # Because iberia_clean is now ONLY the mainland, this intersection acts as a perfect mask!
  ext_cropped <- tryCatch(
    suppressWarnings(st_intersection(ext_sf, iberia_clean)), 
    error = function(e) suppressWarnings(st_crop(ext_sf, st_bbox(iberia_clean)))
  )
  
  map_plot <- ggplot() +
    geom_sf(data = iberia_clean, fill = "gray95", color = "black", linewidth = 0.8) +
    geom_sf(data = ext_cropped, fill = "transparent", color = "#2c3e50", linewidth = 0.1, alpha = 0.5) +
    theme_void() +
    ggtitle(bquote(bold("Training extent:") ~ H[.(gsub("H", "", ext))])) + 
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
  
  # B. VIOLIN PLOTS (No-Crop Logic)
  ext_data <- plot_data %>% filter(Extent == ext)
  strict_paired_data <- ext_data %>% group_by(Metric, Species) %>% filter(n() == 2) %>% ungroup()
  
  if(nrow(strict_paired_data) > 2) {
    sig_legend <- "Significance: ns (p > 0.05), * (p ≤ 0.05), ** (p ≤ 0.01), *** (p ≤ 0.001), **** (p ≤ 0.0001)"
    
    violin_plot <- ggviolin(strict_paired_data, x = "Set", y = "Score", fill = "Set",
                            palette = c("#E69F00", "#56B4E9"), add = "boxplot", 
                            add.params = list(width = 0.1, fill = "white"),
                            trim = FALSE) + 
      facet_wrap(~Metric, scales = "free_y", ncol = 2) +
      stat_compare_means(paired = TRUE, method = "wilcox.test", label = "p.signif", 
                         label.x = 1.35, label.y.npc = "top") +
      theme_minimal() +
      labs(title = bquote("Model Performance Metrics" ~ "(" * H[.(gsub("H", "", ext))] * ")"),
           x = "", y = "Score", caption = sig_legend) +
      theme(plot.title = element_text(face = "bold", size = 16), 
            strip.text = element_text(face = "bold", size = 14),
            plot.caption = element_text(face = "italic", size = 10, hjust = 0, color = "gray30"), 
            legend.position = "none") +
      
      # Using 'oob = scales::squish' and 'expand' to force everything into the viewing window
      ggh4x::facetted_pos_scales(y = list(
        Metric == "AUC"    ~ scale_y_continuous(breaks = seq(0, 1, 0.2), expand = expansion(mult = c(0, 0.1))),
        Metric == "uAUC"   ~ scale_y_continuous(breaks = seq(0, 1, 0.2), expand = expansion(mult = c(0, 0.1))),
        Metric == "CBI"    ~ scale_y_continuous(breaks = seq(-1, 1, 0.5), expand = expansion(mult = 0.1)),
        Metric == "maxTSS" ~ scale_y_continuous(breaks = seq(-1, 1, 0.5), expand = expansion(mult = 0.1))
      ))
  } else {
    violin_plot <- ggplot() + annotate("text", label = "Insufficient Paired Data", x=1, y=1, size=6) + theme_void()
  }
  
  # C. STITCH AND SAVE
  master_plot <- map_plot + violin_plot + plot_layout(widths = c(1, 1.5))
  ggsave(sprintf("Vagenas_aSDMs/output/figures/Global/Extent_Comparison_%s.png", ext), 
         master_plot, width = 15, height = 9, dpi = 300)
  
  sf::sf_use_s2(TRUE)
}

# ==============================================================================
# 4. FINAL STATISTICAL SUMMARY SCRIPT
# ==============================================================================
library(effsize) # Install if needed: install.packages("effsize")

# Helper function to compute stats for a subset of data
compute_stats <- function(df, extent_name) {
  metrics <- c("e_AUC", "CBI", "maxTSS", "uAUC")
  
  results <- lapply(metrics, function(m) {
    hydro_col <- paste0(m, "_Hydroclimatic")
    clim_col  <- paste0(m, "_Climate")
    
    # Filter for pairs only
    subset <- df %>% filter(!is.na(!!sym(hydro_col)) & !is.na(!!sym(clim_col)))
    
    if (nrow(subset) > 1) {
      # Wilcoxon Paired Test
      w_test <- wilcox.test(subset[[hydro_col]], subset[[clim_col]], paired = TRUE)
      # Effect size (Cohen's d)
      d_eff  <- cohen.d(subset[[hydro_col]], subset[[clim_col]], paired = TRUE)
      
      return(data.frame(
        Extent = extent_name,
        Metric = m,
        Mean_Hydro = mean(subset[[hydro_col]]),
        Mean_Clim  = mean(subset[[clim_col]]),
        p_value    = w_test$p.value,
        Significant = w_test$p.value < 0.05,
        Cohen_d    = d_eff$estimate,
        Winner     = ifelse(mean(subset[[hydro_col]]) > mean(subset[[clim_col]]), "Hydroclimatic", "Climate")
      ))
    }
  })
  return(bind_rows(results))
}

# 1. Stats for each Extent
extent_stats <- lapply(unique(comparison_table$Extent), function(e) {
  compute_stats(comparison_table %>% filter(Extent == e), e)
}) %>% bind_rows()

# 2. Stats for Pooled Data
pooled_stats <- compute_stats(comparison_table, "Pooled")

# 3. Combine and Export
final_stats_table <- bind_rows(pooled_stats, extent_stats)

# Save to CSV
write.csv(final_stats_table, 
          "Vagenas_aSDMs/output/tables/FINAL_STATISTICAL_REPORT_GLOBAL_CLIMATE_VS_HYDROCLIMATE.csv", 
          row.names = FALSE)

print(final_stats_table)


#plot Cohen's d

library(ggplot2)
library(dplyr)

# Ensure the Extents are ordered from Fine to Coarse for the X-axis
final_stats_table$Extent <- factor(final_stats_table$Extent, 
                                   levels = c("eco","H5", "H8", "H12"))

# 1. Define the plot
cohen_plot <- ggplot(final_stats_table %>% filter(Extent != "Pooled"), 
                     aes(x = Extent, y = Cohen_d, group = Metric, color = Metric)) +
  # Add horizontal reference line at 0 (No effect)
  geom_hline(yintercept = 0, color = "black", linewidth = 1) +
  # Add lines and points
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  # Add significance indicators (asterisks) based on p-value
  geom_text(aes(label = ifelse(Significant == TRUE, "*", "")), 
            vjust = -1, size = 6, show.legend = FALSE) +
  theme_minimal() +
  labs(
    title = "Effect Size of Hydroclimatic Predictors Across Scales",
    subtitle = "Asterisks (*) denote statistical significance (p < 0.05)",
    y = "Cohen's d (Effect Size)",
    x = "Training Extent (Coarse to Fine)"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    axis.title = element_text(size = 14),
    legend.position = "bottom"
  )

# 2. Save the plot
ggsave("Vagenas_aSDMs/output/figures/Global/Cohen_d_EffectSize.png", 
       cohen_plot, width = 10, height = 6, dpi = 300)

print(cohen_plot)


###compare extents between them

# Ensure Extent is an ordered factor based on your gradient
plot_data$Extent <- factor(plot_data$Extent, levels = c("eco", "H5", "H8", "H12"))

# 1. Define all the specific pairs you want to compare
my_comparisons <- list(
  c("eco", "H5"), c("eco", "H8"), c("eco", "H12"),
  c("H5", "H8"), c("H5", "H12"), c("H8", "H12")
)
adjacent_pairs <- c("eco-H5", "H5-H8", "H8-H12")

# ==============================================================================
# PART 1: OVERALL EXTENT COMPARISON (Custom Bracket & Directional Logic)
# ==============================================================================

# A. Calculate Bounds and Stats for the global plot_data
facet_bounds <- plot_data %>%
  group_by(Metric, Set) %>%
  summarise(
    min_val = min(Score, na.rm = TRUE),
    max_val = max(Score, na.rm = TRUE),
    .groups = "drop"
  )

means_df <- plot_data %>%
  group_by(Metric, Set, Extent) %>%
  summarise(Mean_Score = mean(Score, na.rm = TRUE), .groups = "drop")

stat.test <- plot_data %>%
  group_by(Metric, Set) %>%
  wilcox_test(Score ~ Extent, comparisons = my_comparisons) %>%
  add_significance() %>%
  left_join(means_df, by = c("Metric", "Set", "group1" = "Extent")) %>%
  rename(mean1 = Mean_Score) %>%
  left_join(means_df, by = c("Metric", "Set", "group2" = "Extent")) %>%
  rename(mean2 = Mean_Score) %>%
  mutate(
    direction = ifelse(mean1 > mean2, ">", "<"),
    custom_label = ifelse(p.adj.signif == "ns", "ns", paste0(p.adj.signif, " (", direction, ")")),
    pair_name = paste(group1, group2, sep = "-"),
    is_adjacent = pair_name %in% adjacent_pairs
  ) %>%
  filter(p.adj.signif != "ns") 

# B. Initialize empty dataframes for brackets
stat.test.top <- data.frame()
stat.test.bot <- data.frame()

# C. Calculate dynamic step positioning for the brackets
if (nrow(stat.test) > 0) {
  stat.test <- stat.test %>%
    left_join(facet_bounds, by = c("Metric", "Set")) %>%
    group_by(Metric, Set, is_adjacent) %>%
    mutate(step_rank = row_number()) %>% 
    ungroup() %>%
    mutate(
      range = max_val - min_val,
      step_size = ifelse(range == 0, 0.05, range * 0.08), 
      y.position = ifelse(is_adjacent,
                          min_val - (step_size * 0.5) - (step_rank * step_size),
                          max_val + (step_size * 0.5) + (step_rank * step_size))
    )
  
  stat.test.top <- stat.test %>% filter(!is_adjacent)
  stat.test.bot <- stat.test %>% filter(is_adjacent)
}

# D. Build the base plot

# Define your reversed sequential purple gradient for Extents
extent_colors <- c(
  "eco" = "#3F007D",  # Deep Royal Purple (Stronger baseline)
  "H5"  = "#6A51A3",  # Rich Purple
  "H8"  = "#807DBA",  # Medium Purple (Your old H5)
  "H12" = "#9E9AC8"   # Solid Lilac (No longer washes out to white)
)

# D. Build the base plot
extent_comparison_plot <- ggboxplot(plot_data, x = "Extent", y = "Score", 
                                    color = "Extent", 
                                    palette = extent_colors, # <--- CHANGED HERE
                                    facet.by = c("Metric", "Set"), scales = "free_y",
                                    short.panel.labs = FALSE) +
  stat_compare_means(method = "kruskal.test", label.y.npc = "bottom")

# E. Conditionally add custom top and bottom brackets
if (nrow(stat.test.top) > 0) {
  extent_comparison_plot <- extent_comparison_plot + 
    stat_pvalue_manual(stat.test.top, label = "custom_label", tip.length = 0.01)
}
if (nrow(stat.test.bot) > 0) {
  extent_comparison_plot <- extent_comparison_plot + 
    stat_pvalue_manual(stat.test.bot, label = "custom_label", tip.length = -0.01, vjust = 1.5)
}

# F. Apply formatting and y-axis expansion
extent_comparison_plot <- extent_comparison_plot + 
  scale_y_continuous(
    labels = function(x) ifelse(is.na(x), "", ifelse(x > 1.0, "", x)),
    expand = expansion(mult = c(0.18, 0.18)) # Creates padding for top and bottom brackets
  ) +
  theme_minimal() +
  labs(title = "Performance Across Spatial Extents (All native & invasive widespread species) - [Global aSDMs]",
       subtitle = "Adjacent pairs (bottom brackets) vs Long-jump pairs (top brackets). '>' indicates Left Extent scored higher.",
       x = "Spatial Training Extent", 
       y = "Metric Score",
       caption = "Significance: * (p ≤ 0.05), ** (p ≤ 0.01), *** (p ≤ 0.001), **** (p ≤ 0.0001)") +
  theme(plot.title = element_text(face = "bold", size = 16),
        plot.caption = element_text(hjust = 0, size = 10, face = "italic", color = "gray30"),
        legend.position = "none",
        panel.spacing = unit(1, "lines"))

# G. Render the final plot
plot(extent_comparison_plot)

# H. Save the output
ggsave("Vagenas_aSDMs/output/figures/Global/Extent_Alters_Performance_Global_Overall_Kruskal_and_Wilcox.png", 
       extent_comparison_plot, width = 14, height = 10, dpi = 300)





#add the comparison between predictors at a global scale



library(dplyr)
library(ggplot2)
library(ggpubr)
library(rstatix)

# ==============================================================================
# 1. LOAD AND PREP FOR 2 SETS (GLOBAL)
# ==============================================================================
plot_data <- all_metrics %>%
  filter(Set %in% c("Climate", "Hydroclimatic")) %>% 
  dplyr::select(Species, Extent, Set, e_AUC, CBI, maxTSS, uAUC) %>%
  rename(AUC = e_AUC) %>% 
  pivot_longer(cols = c(AUC, CBI, maxTSS, uAUC), names_to = "Metric", values_to = "Score") %>%
  filter(!is.na(Score))

plot_data$Set <- factor(plot_data$Set, levels = c("Climate", "Hydroclimatic"))
plot_data$Extent <- factor(plot_data$Extent, levels = c("eco", "H5", "H8", "H12"))

# Define your exact custom colors (Only need two!)
set_colors <- c(
  "Climate" = "#E69F00", 
  "Hydroclimatic" = "#56B4E9"
)

# Only ONE comparison is possible now
my_comparisons <- list(c("Climate", "Hydroclimatic"))

# ==============================================================================
# 2. OVERALL SET COMPARISON (Custom Bracket & Directional Logic)
# ==============================================================================

means_df <- plot_data %>%
  group_by(Metric, Extent, Set) %>%
  summarise(Mean_Score = mean(Score, na.rm = TRUE), .groups = "drop")

# Wilcox test across 'Set', grouped by 'Metric' and 'Extent'
stat.test <- plot_data %>%
  group_by(Metric, Extent) %>%
  wilcox_test(Score ~ Set, comparisons = my_comparisons) %>%
  add_significance() %>%
  left_join(means_df, by = c("Metric", "Extent", "group1" = "Set")) %>%
  rename(mean1 = Mean_Score) %>%
  left_join(means_df, by = c("Metric", "Extent", "group2" = "Set")) %>%
  rename(mean2 = Mean_Score) %>%
  mutate(
    direction = ifelse(mean1 > mean2, ">", "<"),
    # CHANGED: Use p.signif instead of p.adj.signif
    custom_label = ifelse(p.signif == "ns", "ns", paste0(p.signif, " (", direction, ")")) 
  ) %>%
  # CHANGED: Use p.signif instead of p.adj.signif
  filter(p.signif != "ns") 

# Add y-position for the single bracket
if (nrow(stat.test) > 0) {
  facet_bounds <- plot_data %>%
    group_by(Metric, Extent) %>%
    summarise(max_val = max(Score, na.rm = TRUE), .groups = "drop")
  
  stat.test <- stat.test %>%
    left_join(facet_bounds, by = c("Metric", "Extent")) %>%
    mutate(y.position = max_val + 0.05)
}

# ==============================================================================
# 3. BUILD THE PLOT
# ==============================================================================

set_comparison_plot <- ggboxplot(plot_data, x = "Set", y = "Score", 
                                 color = "Set", fill = "white", palette = set_colors,
                                 facet.by = c("Metric", "Extent"), scales = "free_y",
                                 short.panel.labs = FALSE) +
  stat_compare_means(method = "wilcox.test", label.y.npc = "bottom")

if (nrow(stat.test) > 0) {
  set_comparison_plot <- set_comparison_plot + 
    stat_pvalue_manual(stat.test, label = "custom_label", tip.length = 0.01)
}

set_comparison_plot <- set_comparison_plot + 
  scale_y_continuous(
    labels = function(x) ifelse(is.na(x), "", ifelse(x > 1.0, "", x)),
    expand = expansion(mult = c(0.18, 0.18)) 
  ) +
  theme_minimal() +
  labs(title = "Performance Across Predictor Sets (All native & invasive widespread species) - [Global aSDMs]",
       subtitle = "'>' indicates Left Set scored higher.",
       x = "Predictor Set", 
       y = "Metric Score",
       caption = "Significance: * (p ≤ 0.05), ** (p ≤ 0.01), *** (p ≤ 0.001), **** (p ≤ 0.0001)") +
  theme(plot.title = element_text(face = "bold", size = 16),
        plot.caption = element_text(hjust = 0, size = 10, face = "italic", color = "gray30"),
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1), 
        panel.spacing = unit(1, "lines"))

plot(set_comparison_plot)

ggsave("Vagenas_aSDMs/output/figures/Global/Set_Alters_Performance_Global_Overall.png", 
       set_comparison_plot, width = 14, height = 10, dpi = 300)





#different for invasive and widespread


library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(rstatix)

# ==============================================================================
# 1. LOAD CLASSIFICATION AND PREP DATA (GLOBAL)
# ==============================================================================

# Load species classification
species_class <- read_csv("Vagenas_aSDMs/input/aSDMs_species_classification.csv", show_col_types = FALSE) %>%
  mutate(Species = gsub(" ", "_", Sp)) %>%
  dplyr::select(Species, Category)

# Merge with all_metrics and prep the data
plot_data_full <- all_metrics %>%
  inner_join(species_class, by = "Species") %>%  
  mutate(Category = ifelse(is.na(Category), "Unknown", Category)) %>% 
  filter(Set %in% c("Climate", "Hydroclimatic")) %>% 
  dplyr::select(Species, Category, Extent, Set, e_AUC, CBI, maxTSS, uAUC) %>%
  rename(AUC = e_AUC) %>% 
  pivot_longer(cols = c(AUC, CBI, maxTSS, uAUC), names_to = "Metric", values_to = "Score") %>%
  filter(!is.na(Score))

# Lock in the factor levels
plot_data_full$Set <- factor(plot_data_full$Set, levels = c("Climate", "Hydroclimatic"))
plot_data_full$Extent <- factor(plot_data_full$Extent, levels = c("eco", "H5", "H8", "H12"))

# Define custom colors
set_colors <- c(
  "Climate" = "#E69F00", 
  "Hydroclimatic" = "#56B4E9"
)

# Only ONE comparison is possible
my_comparisons <- list(c("Climate", "Hydroclimatic"))

# Define the exact categories you want to loop through
categories_to_plot <- c("Native Widespread", "Invasive Widespread")


# ==============================================================================
# 2. RUN THE PLOTTING LOOP
# ==============================================================================

for (cat_name in categories_to_plot) {
  
  cat(sprintf("\n======================================================\n"))
  cat(sprintf("Generating plot for: %s\n", cat_name))
  cat(sprintf("======================================================\n"))
  
  # A. Filter data for the current category
  current_data <- plot_data_full %>% filter(Category == cat_name)
  
  # Skip if there is no data for this category to avoid crashes
  if(nrow(current_data) == 0) {
    cat(sprintf("Skipping %s: No data found.\n", cat_name))
    next
  }
  
  title_text <- sprintf("Performance Across Predictor Sets (%s Global aSDMs)", cat_name)
  file_name <- sprintf("Set_Alters_Performance_Global_%s.png", gsub(" ", "_", cat_name))
  
  # B. Calculate Bounds and Stats
  means_df <- current_data %>%
    group_by(Metric, Extent, Set) %>%
    summarise(Mean_Score = mean(Score, na.rm = TRUE), .groups = "drop")
  
  # Wilcox test across 'Set', grouped by 'Metric' and 'Extent'
  stat.test <- current_data %>%
    group_by(Metric, Extent) %>%
    wilcox_test(Score ~ Set, comparisons = my_comparisons) %>%
    add_significance() %>%
    left_join(means_df, by = c("Metric", "Extent", "group1" = "Set")) %>%
    rename(mean1 = Mean_Score) %>%
    left_join(means_df, by = c("Metric", "Extent", "group2" = "Set")) %>%
    rename(mean2 = Mean_Score) %>%
    mutate(
      direction = ifelse(mean1 > mean2, ">", "<"),
      custom_label = ifelse(p.signif == "ns", "ns", paste0(p.signif, " (", direction, ")")) 
    ) %>%
    filter(p.signif != "ns") 
  
  # Add y-position for the single bracket
  if (nrow(stat.test) > 0) {
    facet_bounds <- current_data %>%
      group_by(Metric, Extent) %>%
      summarise(max_val = max(Score, na.rm = TRUE), .groups = "drop")
    
    stat.test <- stat.test %>%
      left_join(facet_bounds, by = c("Metric", "Extent")) %>%
      mutate(y.position = max_val + 0.05)
  }
  
  # C. Build the base plot
  p <- ggboxplot(current_data, x = "Set", y = "Score", 
                 color = "Set", fill = "white", palette = set_colors,
                 facet.by = c("Metric", "Extent"), scales = "free_y",
                 short.panel.labs = FALSE) +
    stat_compare_means(method = "wilcox.test", label.y.npc = "bottom")
  
  # D. Conditionally add brackets
  if (nrow(stat.test) > 0) {
    p <- p + stat_pvalue_manual(stat.test, label = "custom_label", tip.length = 0.01)
  }
  
  # E. Apply formatting
  p <- p + 
    scale_y_continuous(
      labels = function(x) ifelse(is.na(x), "", ifelse(x > 1.0, "", x)),
      expand = expansion(mult = c(0.18, 0.18)) 
    ) +
    theme_minimal() +
    labs(title = title_text,
         subtitle = "'>' indicates Left Set scored higher.",
         x = "Predictor Set", 
         y = "Metric Score",
         caption = "Significance: * (p ≤ 0.05), ** (p ≤ 0.01), *** (p ≤ 0.001), **** (p ≤ 0.0001)") +
    theme(plot.title = element_text(face = "bold", size = 16),
          plot.caption = element_text(hjust = 0, size = 10, face = "italic", color = "gray30"),
          legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1), 
          panel.spacing = unit(1, "lines"))
  
  # F. Print and Save
  print(p)
  
  save_path <- sprintf("Vagenas_aSDMs/output/figures/Global/%s", file_name)
  ggsave(save_path, p, width = 14, height = 10, dpi = 300)
  
  cat(sprintf("Saved: %s\n", save_path))
}








# ==============================================================================
# PART 2: SPECIES-SPECIFIC TRAJECTORIES (The "Deltas" Gradient)
# ==============================================================================
# The bold red line is the mean trajectory. This visualizes if performance 
# consistently drops or rises as you move from eco -> H5 -> H8 -> H12.

trajectory_plot <- ggplot(plot_data, aes(x = Extent, y = Score)) +
  # Individual species lines (faint)
  geom_line(aes(group = Species), alpha = 0.15, color = "gray40") +
  # Mean trend line (bold)
  stat_summary(aes(group = 1), fun = mean, geom = "line", linewidth = 1.5, color = "#d35400") +
  stat_summary(fun = mean, geom = "point", size = 3, color = "#d35400") +
  facet_grid(Metric ~ Set, scales = "free_y") +
  theme_bw() +
  labs(title = "Species-Level Performance Trajectories Across Scales - [Global aSDMs]",
       subtitle = "Faint lines represent individual species; bold orange line represents the mean",
       x = "Spatial Training Extent Gradient", y = "Metric Score") +
  theme(plot.title = element_text(face = "bold", size = 16),
        strip.text = element_text(face = "bold", size = 12))

plot(trajectory_plot)

ggsave("Vagenas_aSDMs/output/figures/Global/Global_Species_Trajectories_Extent.png", 
       trajectory_plot, width = 10, height = 10, dpi = 300)








# ==============================================================================
#### REGIONAL aSDMs (N=98) | COMPARISONS PERFORMANCE METRICS - SUPPLEMENTARY MATERIAL ####
# ==============================================================================






# ==============================================================================
# 1. LOAD AND PREP REGIONAL DATA
# ==============================================================================
# Load evaluation files
regional_files <- list.files("Vagenas_aSDMs/output/regional/", 
                             pattern = "eval_regional_iberia_strict.csv", 
                             recursive = TRUE, full.names = TRUE)

regional_metrics <- lapply(regional_files, function(f) {
  df <- read_csv(f, show_col_types = FALSE)
  parts <- unlist(strsplit(f, "/"))
  df$Species <- parts[length(parts) - 3]
  df$Extent  <- parts[length(parts) - 2]
  df$Set     <- parts[length(parts) - 1]
  return(df)
}) %>% bind_rows()

# Load species classification
species_class <- read_csv("Vagenas_aSDMs/input/aSDMs_species_classification.csv", 
                          show_col_types = FALSE)

# Standardize the format (replace spaces with underscores in the CSV to match regional_metrics)
species_class <- species_class %>%
  mutate(Species = gsub(" ", "_", Sp)) %>%
  dplyr::select(Species, Category)

# Now execute the master join and format
plot_data_reg <- regional_metrics %>%
  inner_join(species_class, by = "Species") %>%
  filter(Set %in% c("Climate", "Hydroclimatic", "Hydromorphological")) %>%
  dplyr::select(Species, Category, Extent, Set, e_AUC, CBI, maxTSS, uAUC) %>%
  rename(AUC = e_AUC) %>% 
  pivot_longer(cols = c(AUC, CBI, maxTSS, uAUC), names_to = "Metric", values_to = "Score") %>%
  filter(!is.na(Score))

# Verify the final object is ready for plotting
str(plot_data_reg)

# Factorize to ensure logical plotting order
plot_data_reg$Extent <- factor(plot_data_reg$Extent, levels = c("eco", "H5", "H8", "H12"))
plot_data_reg$Set <- factor(plot_data_reg$Set, levels = c("Climate", "Hydroclimatic", "Hydromorphological"))
plot_data_reg$Category <- factor(plot_data_reg$Category, levels = c("Iberian Endemic", "Native Widespread", "Invasive Widespread"))

# ==============================================================================
# 2. MASTER VISUALIZATION: SET PREVALENCE BY EXTENT & SPECIES GROUP
# ==============================================================================
# This creates a massive grid: Rows = Species Category, Columns = Metric
# X-axis = Extent, Color/Fill = Predictor Set

master_boxplot <- ggboxplot(plot_data_reg, x = "Extent", y = "Score", 
                            color = "Set", fill = "Set", alpha = 0.2, 
                            palette = c("#E69F00", "#56B4E9", "#009E73"),
                            outlier.shape = NA) +
  facet_grid(Category ~ Metric, scales = "free_y") +
  theme_bw() +
  labs(title = "Predictor Set Performance Across Spatial Extents and Species Traits",
       subtitle = "Comparing Climate, Hydroclimatic, and Hydromorphological models",
       x = "Spatial Training Extent", y = "Metric Score") +
  theme(plot.title = element_text(face = "bold", size = 16),
        strip.text.y = element_text(face = "bold", size = 10, angle = 270),
        strip.text.x = element_text(face = "bold", size = 12),
        legend.position = "bottom", legend.title = element_blank())

#plot(master_boxplot)

ggsave("Vagenas_aSDMs/output/figures/Regional/Regional_Master_Comparison.png", 
       master_boxplot, width = 16, height = 10, dpi = 300)

# ==============================================================================
# 3. STATISTICAL PROOF: WHICH SET PREVAILS?
# ==============================================================================


statistical_winners <- plot_data_reg %>%
  group_by(Category, Extent, Metric) %>%
  pairwise_wilcox_test(Score ~ Set, p.adjust.method = "holm") %>%
  arrange(Category, Metric, Extent, p.adj)

print(statistical_winners)

write.csv(statistical_winners, 
          "Vagenas_aSDMs/output/tables/Regional_Pairwise_Set_Comparison.csv", 
          row.names = FALSE)

# ==============================================================================
# 4. EXTENT COMPARISON: WHICH RESOLUTION IS BEST PER CATEGORY?
# ==============================================================================
# Find the highest mean score for each extent, grouped by category and metric

extent_summary_reg <- plot_data_reg %>%
  group_by(Category, Metric, Extent) %>%
  summarise(
    Mean_Score = mean(Score, na.rm = TRUE),
    Variance = var(Score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Category, Metric, desc(Mean_Score))

extent_summary_reg

write.csv(extent_summary_reg, 
          "Vagenas_aSDMs/output/tables/Regional_Extent_Summary.csv", 
          row.names = FALSE)


# Build a model#


library(lme4)
library(lmerTest) # Gives p-values for lmer models
library(emmeans)


# Ensure factors are set
plot_data_reg$Extent <- factor(plot_data_reg$Extent, levels = c("eco", "H5", "H8", "H12"))
plot_data_reg$Set <- factor(plot_data_reg$Set, levels = c("Climate", "Hydroclimatic", "Hydromorphological"))
plot_data_reg$Category <- factor(plot_data_reg$Category, levels = c("Iberian Endemic", "Native Widespread", "Invasive Widespread"))

metrics <- c("AUC", "CBI", "maxTSS", "uAUC")
all_emmeans <- list()
all_contrasts <- list()

# ==============================================================================
# 1. RUN THE MIXED MODELS & EXTRACT ESTIMATED MEANS
# ==============================================================================
for (m in metrics) {
  cat(sprintf("Fitting Mixed Model for: %s...\n", m))
  
  # Filter data for current metric
  df_m <- plot_data_reg %>% filter(Metric == m)
  
  # Fit the Linear Mixed-Effects Model
  # We test the three-way interaction: Set * Extent * Category
  # We use (1 | Species) as a random intercept
  model <- lmer(Score ~ Set * Extent * Category + (1 | Species), data = df_m)
  
  # Extract Estimated Marginal Means (EMMs)
  # This calculates the theoretical mean score for every combination
  emm <- emmeans(model, ~ Set | Extent * Category)
  
  # Run pairwise comparisons to see which Set wins inside each Extent/Category
  pairs <- contrast(emm, method = "pairwise", adjust = "tukey")
  
  # Save to lists
  emm_df <- as.data.frame(emm)
  emm_df$Metric <- m
  all_emmeans[[m]] <- emm_df
  
  pairs_df <- as.data.frame(pairs)
  pairs_df$Metric <- m
  all_contrasts[[m]] <- pairs_df
}

final_emmeans <- bind_rows(all_emmeans)
final_contrasts <- bind_rows(all_contrasts)

final_contrasts

# Export the statistical proof
write.csv(final_contrasts, "Vagenas_aSDMs/output/tables/LMM_Pairwise_Contrasts.csv", row.names = FALSE)


# ==============================================================================
# 2. VISUALIZE THE MODEL PREDICTIONS (THE INTERACTION PLOT)
# ==============================================================================
# This plots the model's predicted means (with 95% Confidence Intervals).
# It is much cleaner than raw boxplots because it removes species-level noise.

interaction_plot <- ggplot(final_emmeans, aes(x = Extent, y = emmean, color = Set, group = Set)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2, position = position_dodge(width = 0.5), linewidth = 1) +
  geom_line(position = position_dodge(width = 0.5), alpha = 0.5, linewidth = 1) +
  facet_grid(Metric ~ Category, scales = "free_y") +
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73")) +
  theme_bw() +
  labs(title = "Model-Predicted Performance: Interactions of Scale, Traits, and Predictors",
       subtitle = "Estimated Marginal Means (± 95% CI) from Linear Mixed-Effects Models",
       x = "Spatial Training Extent",
       y = "Predicted Metric Score (EMM)") +
  theme(plot.title = element_text(face = "bold", size = 16),
        strip.text = element_text(face = "bold", size = 11),
        legend.position = "bottom",
        legend.title = element_blank())

plot(interaction_plot)


modelsdir<-"/Users/geo_v/Desktop/Vagenas_aSDMs/output/figures/Regional/Model"

if (!dir.exists(modelsdir)) {
  dir.create(modelsdir, recursive = TRUE)
}


ggsave("Vagenas_aSDMs/output/figures/Regional/Model/LMM_Interaction_Plot.png", 
       interaction_plot, width = 15, height = 10, dpi = 300)


#Interaction plot with labels and statistical comparisons


# ==============================================================================
# 1. CALCULATE THE "UNDISPUTED DOMINANT" SET
# ==============================================================================
# Rank the means to find 1st and 2nd place
plot_labels <- final_emmeans %>%
  group_by(Metric, Category, Extent) %>%
  arrange(desc(emmean)) %>%
  mutate(Rank = row_number()) %>%
  ungroup()

# Initialize the dominance flag
plot_labels$Dominant_Label <- ""

# Loop through to check if the 1st place significantly beats the 2nd place
for (i in 1:nrow(plot_labels)) {
  if (plot_labels$Rank[i] == 1) { # Only check the top performer
    
    m <- as.character(plot_labels$Metric[i])
    c <- as.character(plot_labels$Category[i])
    e <- as.character(plot_labels$Extent[i])
    top_set <- as.character(plot_labels$Set[i])
    
    # Identify the 2nd place set for this specific scenario
    second_set <- plot_labels %>%
      filter(Metric == m, Category == c, Extent == e, Rank == 2) %>%
      pull(Set) %>%
      as.character()
    
    # Create the string patterns that emmeans outputs
    c1 <- paste(top_set, "-", second_set)
    c2 <- paste(second_set, "-", top_set)
    
    # Retrieve the p-value from your final_contrasts table
    match_row <- final_contrasts %>%
      filter(Metric == m, Category == c, Extent == e, 
             (contrast == c1 | contrast == c2))
    
    # If the p-value is < 0.05, it is statistically dominant
    if (nrow(match_row) > 0 && match_row$p.value[1] < 0.05) {
      plot_labels$Dominant_Label[i] <- "★"
    } else {
      plot_labels$Dominant_Label[i] <- "ns" # Optional: remove "ns" if you only want stars
    }
  }
}

# ==============================================================================
# 2. RE-DRAW THE INTERACTION PLOT WITH LABELS
# ==============================================================================
annotated_interaction_plot <- ggplot(plot_labels, aes(x = Extent, y = emmean, color = Set, group = Set)) +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.2, position = position_dodge(width = 0.5), linewidth = 1) +
  geom_line(position = position_dodge(width = 0.5), alpha = 0.5, linewidth = 1) +
  
  # 🚨 THE NEW PART: Add the Dominance Label above the winning point
  geom_text(aes(label = Dominant_Label, y = upper.CL + 0.02), 
            position = position_dodge(width = 0.5), 
            size = 6, show.legend = FALSE, color = "black") +
  
  facet_grid(Metric ~ Category, scales = "free_y") +
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73")) +
  theme_bw() +
  labs(title = "Model-Predicted Performance: Scale, Traits, and Predictors",
       subtitle = "★ indicates the top model significantly outperformed the second-best model (p < 0.05). 'ns' indicates a statistical tie.",
       x = "Spatial Training Extent",
       y = "Predicted Metric Score (EMM)") +
  theme(plot.title = element_text(face = "bold", size = 16),
        strip.text = element_text(face = "bold", size = 11),
        legend.position = "bottom",
        legend.title = element_blank())

ggsave("Vagenas_aSDMs/output/figures/Regional/Model/LMM_Annotated_Interaction_Plot.png", 
       annotated_interaction_plot, width = 15, height = 10, dpi = 300)

print(annotated_interaction_plot)

#Interaction across extents

library(lme4)
library(lmerTest)
library(emmeans)
library(ggplot2)
library(dplyr)
library(tidyr)

# Ensure Extent is an ordered factor for consecutive testing
plot_data_reg$Extent <- factor(plot_data_reg$Extent, levels = c("eco", "H5", "H8", "H12"))

metrics <- c("AUC", "CBI", "maxTSS", "uAUC")
all_trend_contrasts <- list()

# ==============================================================================
# 1. CALCULATE CONSECUTIVE DROPS (THE TREND)
# ==============================================================================
for (m in metrics) {
  cat(sprintf("Calculating Spatial Trends for: %s...\n", m))
  
  df_m <- plot_data_reg %>% filter(Metric == m)
  
  # Same LMM as before
  model <- lmer(Score ~ Set * Extent * Category + (1 | Species), data = df_m)
  
  # Extract EMMs, but this time condition on Set and Category
  emm_extent <- emmeans(model, ~ Extent | Set * Category)
  
  # 'consec' tests the transitions: (H5 - eco), (H8 - H5), (H12 - H8)
  # Reverse = TRUE makes it (eco - H5) so positive numbers mean a drop in performance
  trend_pairs <- contrast(emm_extent, method = "consec", reverse = TRUE)
  
  trend_df <- as.data.frame(trend_pairs)
  trend_df$Metric <- m
  all_trend_contrasts[[m]] <- trend_df
}

final_trends <- bind_rows(all_trend_contrasts)

# Clean up the output for the summary table
parsed_trends <- final_trends %>%
  mutate(
    # Name the transitions clearly
    Transition = case_when(
      contrast == "eco - H5" ~ "1. eco -> H5",
      contrast == "H5 - H8"  ~ "2. H5 -> H8",
      contrast == "H8 - H12" ~ "3. H8 -> H12"
    ),
    # If estimate is positive, performance DROPPED. If negative, it INCREASED.
    Trend_Direction = ifelse(estimate > 0, "Decline", "Improvement"),
    Significant = p.value < 0.05,
    Significance_Label = ifelse(p.value < 0.001, "***", 
                                ifelse(p.value < 0.01, "**", 
                                       ifelse(p.value < 0.05, "*", "ns")))
  ) %>%
  arrange(Metric, Category, Set, Transition)

parsed_trends

write.csv(parsed_trends, 
          "Vagenas_aSDMs/output/tables/LMM_Spatial_Trends.csv", 
          row.names = FALSE)


# ==============================================================================
# 2. VISUALIZE THE DEGRADATION RATE (THE DELTA PLOT)
# ==============================================================================
# This plot shows the EXACT change in score at each transition step.
# Drops below the zero line indicate a loss of performance.

delta_plot <- ggplot(parsed_trends, aes(x = Transition, y = -estimate, color = Set, group = Set)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 1) +
  geom_point(position = position_dodge(width = 0.5), size = 4) +
  geom_line(position = position_dodge(width = 0.5), linewidth = 1.2) +
  geom_errorbar(aes(ymin = -estimate - SE*1.96, ymax = -estimate + SE*1.96), 
                width = 0.2, position = position_dodge(width = 0.5), linewidth = 1) +
  geom_text(aes(label = Significance_Label, y = -estimate + SE*1.96 + 0.015), 
            position = position_dodge(width = 0.5), show.legend = FALSE, size = 5) +
  facet_grid(Metric ~ Category, scales = "free_y") +
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73")) +
  theme_bw() +
  labs(title = "The Resolution Cliff: Step-wise Performance Changes Across Scales",
       subtitle = "Values below zero indicate a drop in performance from the previous spatial extent. Asterisks denote significant changes.",
       x = "Spatial Transition",
       y = "Change in Score (Delta EMM)") +
  theme(plot.title = element_text(face = "bold", size = 16),
        strip.text = element_text(face = "bold", size = 11),
        axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        legend.position = "bottom", legend.title = element_blank())

plot(delta_plot)

ggsave("Vagenas_aSDMs/output/figures/Regional/Model/LMM_Delta_Trend_Plot.png", 
       delta_plot, width = 15, height = 10, dpi = 300)


###Summarize transitions###

# Ensure Extent is an ordered factor
plot_data_reg$Extent <- factor(plot_data_reg$Extent, levels = c("eco", "H5", "H8", "H12"))

metrics <- c("AUC", "CBI", "maxTSS", "uAUC")
all_pairwise_trends <- list()

# ==============================================================================
# 1. CALCULATE ALL PAIRWISE SPATIAL TRANSITIONS
# ==============================================================================
for (m in metrics) {
  cat(sprintf("Calculating All Spatial Transitions for: %s...\n", m))
  
  df_m <- plot_data_reg %>% filter(Metric == m)
  
  # Fit the model
  model <- lmer(Score ~ Set * Extent * Category + (1 | Species), data = df_m)
  
  # Extract EMMs conditioning on Set and Category
  emm_extent <- emmeans(model, ~ Extent | Set * Category)
  
  # 'pairwise' tests ALL combinations: (eco - H5), (eco - H8), (eco - H12), (H5 - H12), etc.
  all_pairs <- contrast(emm_extent, method = "pairwise")
  
  pairs_df <- as.data.frame(all_pairs)
  pairs_df$Metric <- m
  all_pairwise_trends[[m]] <- pairs_df
}

final_all_trends <- bind_rows(all_pairwise_trends)

# ==============================================================================
# 2. PARSE AND CLASSIFY THE TRENDS
# ==============================================================================
parsed_all_trends <- final_all_trends %>%
  mutate(
    # The 'estimate' is calculated as Group 1 - Group 2 (e.g., "eco - H12").
    # If the estimate is positive, it means the finer extent (eco) scored higher,
    # meaning performance DECLINED as you moved to the coarser extent.
    Trend_Direction = case_when(
      p.value >= 0.05 ~ "Insignificant (Stable)",
      p.value < 0.05 & estimate > 0 ~ "Significant Decline",
      p.value < 0.05 & estimate < 0 ~ "Significant Improvement"
    ),
    Significant = p.value < 0.05,
    Significance_Label = ifelse(p.value < 0.001, "***", 
                                ifelse(p.value < 0.01, "**", 
                                       ifelse(p.value < 0.05, "*", "ns")))
  ) %>%
  # Select and reorder columns for a clean thesis table
  dplyr::select(Metric, Category, Set, contrast, estimate, SE, p.value, Significant, Trend_Direction, Significance_Label) %>%
  arrange(Metric, Category, Set, contrast)

parsed_all_trends

# Export the master summary table
write.csv(parsed_all_trends, 
          "Vagenas_aSDMs/output/tables/LMM_All_Spatial_Transitions_Summary.csv", 
          row.names = FALSE)


##compare them

library(lme4)
library(lmerTest)
library(emmeans)
library(dplyr)
library(tidyr)
library(ggplot2)

# ==============================================================================
# 1. MERGE CATEGORIES & RUN THE LMM
# ==============================================================================
# Combine Native and Invasive into a single "Widespread" category
plot_data_combined <- plot_data_reg %>%
  mutate(Category = ifelse(Category %in% c("Native Widespread", "Invasive Widespread"), 
                           "Widespread", 
                           as.character(Category)))

metrics <- c("AUC", "CBI", "maxTSS", "uAUC")
all_emmeans <- list()
all_contrasts <- list()

for (m in metrics) {
  cat(sprintf("Fitting Combined LMM for: %s...\n", m))
  df_m <- plot_data_combined %>% filter(Metric == m)
  
  model <- lmer(Score ~ Set * Extent * Category + (1 | Species), data = df_m)
  emm <- emmeans(model, ~ Set | Extent * Category)
  pairs <- contrast(emm, method = "pairwise", adjust = "tukey")
  
  emm_df <- as.data.frame(emm)
  emm_df$Metric <- m
  all_emmeans[[m]] <- emm_df
  
  pairs_df <- as.data.frame(pairs)
  pairs_df$Metric <- m
  all_contrasts[[m]] <- pairs_df
}

final_emmeans_comb <- bind_rows(all_emmeans)
final_contrasts_comb <- bind_rows(all_contrasts)

# ==============================================================================
# 2. CALCULATE THE SUPERIORITY MATRIX
# ==============================================================================
leaderboard <- final_emmeans_comb %>%
  group_by(Metric, Category, Extent) %>%
  arrange(desc(emmean)) %>%
  mutate(Rank = row_number()) %>%
  filter(Rank %in% c(1, 2)) %>%
  summarise(
    Dominant_Set = as.character(Set[Rank == 1]),
    Dominant_Mean = round(emmean[Rank == 1], 4),
    RunnerUp_Set = as.character(Set[Rank == 2]),
    .groups = "drop"
  ) %>%
  mutate(
    c1 = paste(Dominant_Set, "-", RunnerUp_Set),
    c2 = paste(RunnerUp_Set, "-", Dominant_Set)
  )

clean_contrasts <- final_contrasts_comb %>% dplyr::select(Metric, Category, Extent, contrast, estimate, p.value)

superiority_table <- leaderboard %>%
  left_join(clean_contrasts, by = c("Metric", "Category", "Extent", "c1" = "contrast")) %>%
  left_join(clean_contrasts, by = c("Metric", "Category", "Extent", "c2" = "contrast"), suffix = c("_match1", "_match2")) %>%
  mutate(
    raw_estimate = coalesce(estimate_match1, estimate_match2),
    p_value = coalesce(p.value_match1, p.value_match2),
    Absolute_Advantage = round(abs(raw_estimate), 4)
  )

# ==============================================================================
# 3. FILTER AND PLOT (FORCING THE FULL GRID)
# ==============================================================================
significant_superiority <- superiority_table %>%
  filter(p_value < 0.05) 

# CRITICAL: Define all factor levels explicitly so ggplot knows they exist
significant_superiority$Extent <- factor(significant_superiority$Extent, levels = c("eco", "H5", "H8", "H12"))
significant_superiority$Metric <- factor(significant_superiority$Metric, levels = c("AUC", "CBI", "maxTSS", "uAUC"))
significant_superiority$Category <- factor(significant_superiority$Category, levels = c("Iberian Endemic", "Widespread"))
significant_superiority$Dominant_Set <- factor(significant_superiority$Dominant_Set, levels = c("Climate", "Hydroclimatic", "Hydromorphological"))

dominance_plot_combined <- ggplot(significant_superiority, aes(x = Extent, y = Metric, fill = Dominant_Set)) +
  geom_tile(color = "white", linewidth = 1.5) +
  geom_text(aes(label = sprintf("+%.3f", Absolute_Advantage)), 
            color = "black", fontface = "bold", size = 4.5) +
  facet_grid(~ Category, drop = FALSE) + # Forces all facet boxes to stay
  scale_x_discrete(drop = FALSE) +       # Forces H8 to stay visible
  scale_y_discrete(drop = FALSE) +       # Forces CBI to stay visible
  scale_fill_manual(values = c("Climate" = "#E69F00", 
                               "Hydroclimatic" = "#56B4E9", 
                               "Hydromorphological" = "#009E73"),
                    drop = FALSE) +      # Forces the full legend
  theme_minimal() +
  labs(
    title = "Statistically Dominant Predictor Sets: Endemic vs. Widespread",
    subtitle = "Colored tiles represent the undisputed winner (p < 0.05). Numbers indicate the absolute EMM advantage.\nBlank areas indicate a statistical tie.",
    x = "Spatial Training Extent",
    y = "Performance Metric",
    fill = "Winning Predictor Set:"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    subtitle = element_text(size = 11, color = "gray30"),
    strip.text = element_text(face = "bold", size = 13),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 13, face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 11),
    panel.grid = element_blank(), 
    panel.background = element_rect(fill = "gray95", color = NA)
  )


plot(dominance_plot_combined)



#second vs third


# ==============================================================================
# 1. EXTRACT ALL THREE RANKS
# ==============================================================================
# We need the 1st, 2nd, and 3rd place sets for every scenario
leaderboard_full <- final_emmeans_comb %>%
  group_by(Metric, Category, Extent) %>%
  arrange(desc(emmean)) %>%
  mutate(Rank = row_number()) %>%
  summarise(
    Rank1_Set = as.character(Set[Rank == 1]),
    Rank2_Set = as.character(Set[Rank == 2]),
    Rank3_Set = as.character(Set[Rank == 3]),
    .groups = "drop"
  )

# ==============================================================================
# 2. CALCULATE DUAL CONTRASTS (vs 2nd and vs 3rd)
# ==============================================================================
# Pivot the table so we have two rows per scenario: one comparing 1st to 2nd, one 1st to 3rd
comp_df <- leaderboard_full %>%
  pivot_longer(cols = c(Rank2_Set, Rank3_Set),
               names_to = "Comparison_Level",
               values_to = "Compared_Set") %>%
  mutate(
    # Create clean labels for the Y-axis
    Comparison_Label = ifelse(Comparison_Level == "Rank2_Set", "vs. 2nd Place", "vs. 3rd Place"),
    # Construct the exact strings needed to match the emmeans contrasts
    c_a = paste(Rank1_Set, "-", Compared_Set),
    c_b = paste(Compared_Set, "-", Rank1_Set)
  )

# Join the p-values and estimates from your previously generated 'clean_contrasts' table
superiority_dual <- comp_df %>%
  left_join(clean_contrasts, by = c("Metric", "Category", "Extent", "c_a" = "contrast")) %>%
  left_join(clean_contrasts, by = c("Metric", "Category", "Extent", "c_b" = "contrast"), suffix = c("_a", "_b")) %>%
  mutate(
    raw_estimate = coalesce(estimate_a, estimate_b),
    p_value = coalesce(p.value_a, p.value_b),
    Absolute_Advantage = round(abs(raw_estimate), 4),
    Significant = p_value < 0.05
  ) %>%
  dplyr::select(Metric, Category, Extent, Rank1_Set, Comparison_Label, Compared_Set, Absolute_Advantage, p_value, Significant)

# Export the detailed dual-layer CSV
# write.csv(superiority_dual, 
#           "/Users/georgevagenas/Desktop/Vagenas_aSDMs/output/regional/Dual_Superiority_Matrix.csv", 
#           row.names = FALSE)

# ==============================================================================
# 3. VISUALIZE THE DUAL-LAYER DOMINANCE MATRIX
# ==============================================================================
# Filter for true significance
plot_data_dual <- superiority_dual %>%
  filter(Significant == TRUE) 

# STRICT factor ordering to maintain grid integrity
plot_data_dual$Extent <- factor(plot_data_dual$Extent, levels = c("eco", "H5", "H8", "H12"))
plot_data_dual$Metric <- factor(plot_data_dual$Metric, levels = c("AUC", "CBI", "maxTSS", "uAUC"))
plot_data_dual$Category <- factor(plot_data_dual$Category, levels = c("Iberian Endemic", "Widespread"))
plot_data_dual$Rank1_Set <- factor(plot_data_dual$Rank1_Set, levels = c("Climate", "Hydroclimatic", "Hydromorphological"))

# We reverse the order of the comparison label so "vs 2nd" sits physically above "vs 3rd" in the grid
plot_data_dual$Comparison_Label <- factor(plot_data_dual$Comparison_Label, levels = c("vs. 3rd Place", "vs. 2nd Place"))

dual_dominance_plot <- ggplot(plot_data_dual, aes(x = Extent, y = Comparison_Label, fill = Rank1_Set)) +
  geom_tile(color = "white", linewidth = 1.2) +
  geom_text(aes(label = sprintf("+%.3f", Absolute_Advantage)), 
            color = "black", fontface = "bold", size = 4) +
  # Facet by Metric vertically, and Category horizontally
  facet_grid(Metric ~ Category, drop = FALSE, scales = "free_y") +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  scale_fill_manual(values = c("Climate" = "#E69F00", 
                               "Hydroclimatic" = "#56B4E9", 
                               "Hydromorphological" = "#009E73"),
                    drop = FALSE) +
  theme_minimal() +
  labs(
    title = "The Penalty of Error: Model Advantage Over 2nd and 3rd Place Alternatives",
    subtitle = "Colored tiles represent the undisputed Champion (p < 0.05). Numbers indicate how much it outperformed the alternative.\nBlank upper tiles indicate a tie for 1st. Blank lower tiles indicate all models performed identically.",
    x = "Spatial Training Extent",
    y = "",
    fill = "Undisputed Champion:"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    subtitle = element_text(size = 10, color = "gray30"),
    strip.text.x = element_text(face = "bold", size = 12),
    strip.text.y = element_text(face = "bold", size = 12, angle = 270),
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 10, face = "italic", color = "gray20"),
    axis.title = element_text(size = 12, face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    panel.grid = element_blank(), 
    panel.spacing.y = unit(0.5, "lines"), # Add slight breathing room between metrics
    panel.background = element_rect(fill = "gray95", color = NA)
  )

plot(dual_dominance_plot)

ggsave("Vagenas_aSDMs/output/figures/Regional/Model/LMM_Dual_Dominance_Heatmap.png", 
       dual_dominance_plot, width = 13, height = 9, dpi = 300)





# ================================================================================================================
#### VARIABLE IMPORTANCE | FIGURE 3 | COMPARISONS ACROSS PREDICTOR SETS & EXTENT - SUPPLEMENTARY MATERIAL ####
# ================================================================================================================




setwd("/Users/geo_v/Desktop/")

# 1. Load classification
species_class <- read_csv("Vagenas_aSDMs/input/aSDMs_species_classification.csv", show_col_types = FALSE) %>%
  mutate(Species = gsub(" ", "_", Sp)) %>%
  dplyr::select(Species, Category)

# 2. Extract Data (Robust for all models)
vi_files <- list.files("Vagenas_aSDMs/output/regional/", 
                       pattern = "varImp_regional.rds", recursive = TRUE, full.names = TRUE)


# 1. Pre-allocate an empty list to exactly the length of your files
vi_list <- vector("list", length(vi_files))

# 2. Extract files one by one (High-Speed & Memory-Safe)
for(i in seq_along(vi_files)) {
  f <- vi_files[i]
  parts <- unlist(strsplit(f, "/"))
  sp_name <- parts[length(parts) - 3]
  ext_name <- parts[length(parts) - 2]
  set_name <- parts[length(parts) - 1] 
  
  vi_obj <- readRDS(f)
  
  # Extract Dataframe based on object structure
  vi_df <- tryCatch({
    if(inherits(vi_obj, ".varImportance")) {
      vi_obj@varImportance
    } else if (inherits(vi_obj, ".varImportanceList") || is.list(vi_obj)) {
      do.call(rbind, lapply(vi_obj, function(x) {
        if(inherits(x, ".varImportance")) x@varImportance else as.data.frame(x)
      })) %>% group_by(variables) %>% summarise(AUCtest = mean(AUCtest, na.rm = TRUE), .groups = "drop")
    } else {
      as.data.frame(vi_obj)
    }
  }, error = function(e) return(NULL))
  
  if(!is.null(vi_df)) {
    names(vi_df)[names(vi_df) == "variables"] <- "Variable"
    names(vi_df)[names(vi_df) == "AUCtest"] <- "Mean_Imp"
    
    # Drop the temporary dataframe directly into the list slot
    vi_list[[i]] <- vi_df %>%
      mutate(Importance_Prop = Mean_Imp / sum(Mean_Imp, na.rm = TRUE),
             Species = sp_name, Extent = ext_name, Set = set_name)
  }
  
  # Clear the heavy S4 object to keep memory clean, but DO NOT run gc()
  rm(vi_obj, vi_df)
}

# 3. Bind everything together EXACTLY ONCE at the end
vi_master <- bind_rows(vi_list)

# 4. Run garbage collection once at the end
gc()


#exploration of patterns - shifts across categories - Version#1

library(sf)
library(ggtext) 
library(cowplot)

# ==============================================================================
# 0. SPATIAL PATHS & PREPARATION
# ==============================================================================
iberia_path <- "Vagenas_aSDMs/input/SpatialExtents/study_area_iberia.shp"
extent_paths <- list(
  "eco" = "Vagenas_aSDMs/input/SpatialExtents/ecoregions/feow_hydrosheds.shp",
  "H5"  = "Vagenas_aSDMs/input/SpatialExtents/HydroSHEDS/H5/hybas_lake_eu_lev05_v1c/hybas_lake_eu_lev05_v1c.shp",
  "H8"  = "Vagenas_aSDMs/input/SpatialExtents/HydroSHEDS/H8/hybas_lake_eu_lev08_v1c/hybas_lake_eu_lev08_v1c.shp",
  "H12" = "Vagenas_aSDMs/input/SpatialExtents/HydroSHEDS/H12/hybas_lake_eu_lev12_v1c/hybas_lake_eu_lev12_v1c.shp"
)

iberia_sf <- st_read(iberia_path, quiet = TRUE)
sf::sf_use_s2(FALSE)
iberia_clean <- st_zm(iberia_sf, drop = TRUE, what = "ZM")

ordered_extents <- c("eco", "H5", "H8", "H12")

# ==============================================================================
# 1. FETCH AND CREATE THE IBERIAN MAINLAND MASK
# ==============================================================================
# Fetch Spain and Portugal
iberia_sf <- ne_countries(country = c("Spain", "Portugal"), scale = "medium", returnclass = "sf")

# Convert to terra to easily isolate the mainland (drops the islands)
iberia_vect <- vect(iberia_sf)
iberia_dissolved <- aggregate(iberia_vect)
iberia_parts <- disagg(iberia_dissolved)
land_areas <- expanse(iberia_parts, unit = "km")
iberia_mainland <- iberia_parts[which.max(land_areas), ]

# Convert back to sf and strip Z/M dimensions for perfect topological masking
iberia_clean <- st_as_sf(iberia_mainland)
iberia_clean <- st_zm(iberia_clean, drop = TRUE, what = "ZM")


# ==============================================================================
# 2. GENERATE THE 4 MAPS (Dynamic Linewidth & Perfectly Cropped)
# ==============================================================================
sf::sf_use_s2(FALSE) # Turn off spherical geometry to prevent intersection errors

map_plots <- list()
extent_titles <- c("eco" = "**eco**", "H5" = "**H<sub>5</sub>**", "H8" = "**H<sub>8</sub>**", "H12" = "**H<sub>12</sub>**")

for (ext in ordered_extents) {
  cat(sprintf("Generating map for: %s\n", ext))
  
  ext_sf <- st_read(extent_paths[[ext]], quiet = TRUE)
  ext_sf <- st_zm(ext_sf, drop = TRUE, what = "ZM")
  
  # Align projection to the mask
  st_crs(ext_sf) <- st_crs(iberia_clean)
  ext_sf <- suppressWarnings(st_buffer(st_make_valid(ext_sf), 0))
  
  # THE CROP: This line slices your shapefiles using the Iberia mask!
  ext_cropped <- tryCatch(
    suppressWarnings(st_intersection(ext_sf, iberia_clean)), 
    error = function(e) suppressWarnings(st_crop(ext_sf, st_bbox(iberia_clean)))
  )
  
  # Set an ultra-thin line for the denser H8 and H12 shapefiles
  current_lwd <- ifelse(ext %in% c("H8", "H12"), 0.03, 0.08)
  
  map_plots[[ext]] <- ggplot() +
    geom_sf(data = iberia_clean, fill = "gray95", color = "black", linewidth = 0.3) +
    geom_sf(data = ext_cropped, fill = "transparent", color = "#2c3e50", linewidth = current_lwd, alpha = 0.5) +
    theme_void() +
    labs(title = extent_titles[[ext]]) +
    theme(plot.title = element_markdown(hjust = 0.5, size = 14, margin = ggplot2::margin(t=0, r=0, b=10, l=0))) 
}

sf::sf_use_s2(TRUE) # Turn spherical geometry back on

# ==============================================================================
# 1. CLEAN CLASSIFICATION & MERGE CATEGORY
# ==============================================================================
vi_master_new <- vi_master %>%
  mutate(
    Filter_Type = case_when(
      grepl("global|regional|suitability", Variable, ignore.case = TRUE) ~ "Global niche",
      grepl("hydro", Variable, ignore.case = TRUE) ~ "Hydroclimatic", 
      grepl("lka|dor|sgr|urb|for|morpho", Variable, ignore.case = TRUE) ~ "Hydromorphology",
      grepl("bio|clima|precip|temp", Variable, ignore.case = TRUE) ~ "Climate",
      TRUE ~ "Global niche" 
    )
  ) %>%
  left_join(species_class, by = "Species") %>%
  mutate(Category = ifelse(is.na(Category), "Unknown", Category))

# ==============================================================================
# 2. GENERATE THE TRAJECTORY DATA (WITH 95% CI)
# ==============================================================================
trajectory_data <- vi_master_new %>%
  mutate(Extent = factor(Extent, levels = ordered_extents)) %>%
  group_by(Category, Set, Extent, Variable, Filter_Type) %>%
  summarise(
    # 1. Count the models to calculate Standard Error properly
    N_Models = n(), 
    
    Mean_Importance = mean(Importance_Prop, na.rm = TRUE),
    SD_Importance = sd(Importance_Prop, na.rm = TRUE),
    
    # 2. Calculate the exact 95% Confidence Interval
    SE_Importance = SD_Importance / sqrt(N_Models),
    CI_95 = SE_Importance * 1.96, 
    
    .groups = "drop"
  )

trajectory_data <- trajectory_data %>%
  mutate(Category = factor(Category, levels = c("Iberian Endemic", "Native Widespread", "Invasive Widespread")))

# ==============================================================================
# 2.5 BUILD THE SHIFT PLOT (THE MISSING LINES)
# ==============================================================================
shift_plot <- ggplot(trajectory_data, aes(x = Extent, y = Mean_Importance, group = Variable, color = Filter_Type)) +
  
  # Apply the CI_95 to the error bars
  geom_errorbar(aes(ymin = pmax(0, Mean_Importance - CI_95), 
                    ymax = Mean_Importance + CI_95), 
                width = 0.15, alpha = 0.6, linewidth = 0.8) +
  
  geom_line(linewidth = 1.2, alpha = 0.8) +
  geom_point(size = 2.5) +
  facet_grid(Set ~ Category) +
  
  scale_x_discrete(
    breaks = c("eco", "H5", "H8", "H12"),
    labels = c(
      expression(bold("eco")), 
      expression(bold(H[5])), 
      expression(bold(H[8])), 
      expression(bold(H[12]))
    )
  ) + 
  
  scale_y_continuous(labels = function(x) paste0(round(x * 100, 0), "%")) +
  
  scale_color_manual(values = c("Climate" = "#E69F00", 
                                "Hydroclimatic" = "#56B4E9", 
                                "Hydromorphology" = "#009E73",
                                "Global niche" = "gray50")) +
  theme_minimal() +
  labs(x = NULL, y = "Mean Relative Importance (%)") +
  theme(
    strip.text.x = element_text(face = "bold", size = 14),
    strip.text.y = element_text(face = "bold", size = 14),
    
    axis.text.x = element_text(size = 14, color = "black"),
    axis.title.y = element_text(size = 16, face = "bold", margin = margin(r = 10)),
    panel.grid.minor = element_blank(),
    
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 14, face = "bold") 
  )

# ==============================================================================
# 3. STITCH EVERYTHING TOGETHER AND ADD THE OVERLAY BOXES
# ==============================================================================

# Create the map row with safety wrappers and increased horizontal spacing
map_row <- (
  (wrap_elements(map_plots[["eco"]] + labs(title = NULL)) + 
     ggtitle(expression(bold("eco"))) + 
     theme(plot.title = element_text(face = "bold", hjust = 0.5))) | 
    
    (wrap_elements(map_plots[["H5"]] + labs(title = NULL)) + 
       ggtitle(expression(bold(H[5]))) + 
       theme(plot.title = element_text(face = "bold", hjust = 0.5))) | 
    
    (wrap_elements(map_plots[["H8"]] + labs(title = NULL)) + 
       ggtitle(expression(bold(H[8]))) + 
       theme(plot.title = element_text(face = "bold", hjust = 0.5))) | 
    
    (wrap_elements(map_plots[["H12"]] + labs(title = NULL)) + 
       ggtitle(expression(bold(H[12]))) + 
       theme(plot.title = element_text(face = "bold", hjust = 0.5)))
) & theme(plot.margin = margin(t = 5, r = 20, b = 5, l = 20, unit = "pt"))

# Compose the main layout with patchwork
composed_plot <- shift_plot / map_row + 
  plot_layout(heights = c(3.5, 1))

# Convert to a drawing canvas and use annotate() to draw the rectangles
master_plot_with_boxes <- ggdraw(composed_plot) +
  
  # Add Climate Box (Orange)
  annotate("rect", xmin = 0.06, xmax = 0.96, ymin = 0.74, ymax = 0.95, 
           color = "black", fill = NA, linewidth = 0.5) +
  
  # Add Hydroclimatic Box (Blue)
  annotate("rect", xmin = 0.06, xmax = 0.96, ymin = 0.52, ymax = 0.73, 
           color = "black", fill = NA, linewidth = 0.5) +
  
  # Add Hydromorphology Box (Green)
  annotate("rect", xmin = 0.06, xmax = 0.96, ymin = 0.305, ymax = 0.512, 
           color = "black", fill = NA, linewidth = 0.5)


varimpdir<-"/Users/geo_v/Desktop/Vagenas_aSDMs/output/figures/Regional/VarImp"

if (!dir.exists(varimpdir)) {
  dir.create(varimpdir, recursive = TRUE)
}

# Save the final plot
ggsave("Vagenas_aSDMs/output/figures/Regional/VarImp/Master_TrajectorShifts_With_Boxes.png", 
       master_plot_with_boxes, width = 16, height = 13, dpi = 300)


#Radar plot - shifts across categories - Version#2


# 1. Clean Classification (Robust matching) & Merge Category
vi_master_new <- vi_master %>%
  mutate(
    Filter_Type = case_when(
      grepl("global|regional|suitability", Variable, ignore.case = TRUE) ~ "Global niche",
      grepl("hydro", Variable, ignore.case = TRUE) ~ "Hydroclimatic", 
      grepl("lka|dor|sgr|urb|for|morpho", Variable, ignore.case = TRUE) ~ "Hydromorphology",
      grepl("bio|clima|precip|temp", Variable, ignore.case = TRUE) ~ "Climate",
      TRUE ~ "Global niche" 
    )
  ) %>%
  left_join(species_class, by = "Species") %>%
  mutate(Category = ifelse(is.na(Category), "Unknown", Category))

# ==============================================================================
# 2. AGGREGATE DATA (The 5/6 Line Guarantee)
# ==============================================================================
radar_data_vars <- vi_master_new %>%
  mutate(
    Extent = factor(Extent, levels = c("H12", "H8", "H5", "eco")),
    # By renaming all global variables to this single string, we guarantee exactly 1 grey line
    Variable = ifelse(Filter_Type == "Global niche", "Global_Constraint", Variable)
  ) %>%
  group_by(Category, Set, Extent, Variable, Filter_Type) %>%
  summarise(Mean_Importance = mean(Importance_Prop, na.rm = TRUE), .groups = "drop")

# ==============================================================================
# 3. GENERATE THE VARIABLE-LEVEL RADAR PLOT
# ==============================================================================
ring_labels <- data.frame(
  x = 4.5, 
  y = c(0.1, 0.2, 0.3, 0.4, 0.5),
  label = c("10%", "20%", "30%", "40%", "50%")
)

radar_plot_vars <- ggplot(radar_data_vars, aes(x = Extent, y = Mean_Importance, group = Variable)) +
  
  # THE FIX: Use geom_line instead of geom_polygon so it doesn't draw lines across the middle
  geom_line(aes(color = Filter_Type), linewidth = 1.2, alpha = 0.8) +
  geom_point(aes(color = Filter_Type), size = 2) +
  
  geom_text(data = ring_labels, aes(x = x, y = y, label = label), 
            color = "gray40", size = 3.5, fontface = "bold", vjust = 0.5, 
            inherit.aes = FALSE) +
  
  coord_polar() +
  facet_grid(Set ~ Category) + 
  
  # Perfect Font Matching via HTML Markdown
  scale_x_discrete(
    breaks = c("H12", "H8", "H5", "eco"),
    labels = c(
      "H12" = "**H<sub>12</sub>**", 
      "H8"  = "**H<sub>8</sub>**", 
      "H5"  = "**H<sub>5</sub>**", 
      "eco" = "**eco**"
    )
  ) +
  
  scale_y_continuous(breaks = seq(0.1, 0.6, by = 0.1), limits = c(0, NA)) +
  
  scale_color_manual(values = c("Climate" = "#E69F00", 
                                "Hydroclimatic" = "#56B4E9", 
                                "Hydromorphology" = "#009E73",
                                "Global niche" = "gray50")) +
  theme_minimal() +
  
  labs(title = NULL, subtitle = NULL, x = NULL, y = NULL) +
  
  theme(
    strip.text.x = element_text(face = "bold", size = 14),
    strip.text.y = element_text(face = "bold", size = 14),
    
    # Allows ggtext to render the HTML labels beautifully
    axis.text.x = element_markdown(size = 14, color = "black"),
    
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    
    panel.grid.major = element_line(color = "gray80", linetype = "dashed"),
    
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 14, face = "bold"),
    
    panel.spacing.x = unit(8, "lines"), 
    panel.spacing.y = unit(3, "lines"),
    
    plot.margin = ggplot2::margin(t = 15, r = 5, b = 15, l = 5) 
  )

# ==============================================================================
# 4. STITCH EVERYTHING TOGETHER WITH PATCHWORK AND ADD BOXES
# ==============================================================================
# First, compose the layout
composed_radar_plot <- radar_plot_vars / map_row + 
  plot_layout(heights = c(3.5, 1))

# Second, convert to a drawing canvas and apply the transparent boxes
master_radar_with_boxes <- ggdraw(composed_radar_plot) +
  
  # Add Climate Box (Orange)
  annotate("rect", xmin = 0.11, xmax = 0.865, ymin = 0.74, ymax = 0.95, 
           color = "black", fill = NA, linewidth = 0.5) +
  
  # Add Hydroclimatic Box (Blue)
  annotate("rect", xmin = 0.11, xmax = 0.865, ymin = 0.52, ymax = 0.73, 
           color = "black", fill = NA, linewidth = 0.5) +
  
  # Add Hydromorphology Box (Green)
  annotate("rect", xmin = 0.11, xmax = 0.865, ymin = 0.2805, ymax = 0.512, 
           color = "black", fill = NA, linewidth = 0.5)

# print(master_radar_with_boxes)

# Save the final overlaid plot
ggsave("Vagenas_aSDMs/output/figures/Regional/VarImp/Master_Radar_With_Maps.png", 
       master_radar_with_boxes, width = 14, height = 13, dpi = 300)


#Figure 3 - Addition of Variable Names for the exploration of variables


library(forcats) # For fct_reorder

# ==============================================================================
# 1. Clean Classification & Exact Name Dictionary
# ==============================================================================
vi_master_new <- vi_master %>%
  mutate(
    Filter_Type = case_when(
      grepl("global|regional|suitability", Variable, ignore.case = TRUE) ~ "Global niche",
      grepl("hydro", Variable, ignore.case = TRUE) ~ "Hydroclimatic", 
      grepl("lka|dor|sgr|urb|for|morpho", Variable, ignore.case = TRUE) ~ "Hydromorphology",
      grepl("bio|clima|precip|temp", Variable, ignore.case = TRUE) ~ "Climate",
      TRUE ~ "Global niche" 
    )
  ) %>%
  left_join(species_class, by = "Species") %>%
  mutate(
    Category = ifelse(is.na(Category), "Unknown", Category),
    
    # Map the exact names from the ECMWF-CMIP5 / HydroSHEDS lookup table
    Var_Label = case_when(
      # Climate
      Variable == "bio1_clima" ~ "Annual Mean Air Temp.",
      Variable == "bio4_clima" ~ "Air Temp. Seasonality",
      Variable == "bio5_clima" ~ "Max Air Temp. Warmest Month",
      Variable == "bio6_clima" ~ "Min Air Temp. Coldest Month",
      Variable == "bio12_clima" ~ "Mean Annual Precip.",
      Variable == "bio15_clima" ~ "Mean Precip. Seasonality",
      Variable == "bio16_clima" ~ "Mean Precip. Wettest Qtr",
      Variable == "bio17_clima" ~ "Mean Precip. Driest Qtr",
      
      # Hydroclimate
      Variable == "bio1_hydro" ~ "Annual Mean Water Temp.",
      Variable == "bio4_hydro" ~ "Water Temp. Seasonality",
      Variable == "bio5_hydro" ~ "Max Water Temp. Warmest Month",
      Variable == "bio6_hydro" ~ "Min Water Temp. Coldest Month",
      Variable == "bio12_hydro" ~ "Mean Annual Streamflow",
      Variable == "bio15_hydro" ~ "Mean Streamflow Seasonality",
      Variable == "bio16_hydro" ~ "Mean Streamflow Wettest Qtr",
      Variable == "bio17_hydro" ~ "Mean Streamflow Driest Qtr",
      
      # Hydromorphology
      Variable == "dor_pc_pva" ~ "Degree of Regulation (%)",
      Variable == "lka_pc_use" ~ "Limnicity (%)",
      Variable == "urb_pc_use" ~ "Urban Extent (%)",
      Variable == "for_pc_use" ~ "Forest Cover (%)",
      Variable == "sgr_dk_rav" ~ "Stream Gradient",
      Variable == "pac_pc_use" ~ "Protected area Extent (%)",
      
      # Global Niche
      grepl("global", Variable, ignore.case = TRUE) ~ "Global Suitability",
      
      TRUE ~ stringr::str_to_title(gsub("_", " ", Variable)) 
    )
  )

# ==============================================================================
# 2. GENERATE THE POINT-SHIFT PLOT (HORIZONTAL, FILTERED & ORDERED)
# ==============================================================================
trajectory_data <- vi_master_new %>%
  mutate(Extent = factor(Extent, levels = ordered_extents)) %>%
  group_by(Category, Set, Extent, Variable, Var_Label) %>%
  summarise(
    N_Models = n(), 
    Mean_Importance = mean(Importance_Prop, na.rm = TRUE),
    SD_Importance = sd(Importance_Prop, na.rm = TRUE),
    SE_Importance = SD_Importance / sqrt(N_Models),
    CI_95 = SE_Importance * 1.96, 
    .groups = "drop"
  ) %>%
  # Calculate the overall average importance for each variable WITHIN its Set
  group_by(Set, Var_Label) %>%
  mutate(Set_Avg_Imp = mean(Mean_Importance, na.rm = TRUE)) %>%
  ungroup() %>%
  # Create a unique identifier so "Global Suitability" can be ranked differently in each Set
  mutate(
    Plot_Label = paste0(Var_Label, "__", Set),
    Plot_Label = fct_reorder(Plot_Label, Set_Avg_Imp) # Orders from highest to lowest on Y-axis
  )

# ==============================================================================
# 2.5 CUSTOM Y-AXIS ORDERING 
# (Global Niche first, then ranked by 'eco' Mean_Importance)
# ==============================================================================

# 1. Isolate the 'eco' extent to figure out the ranking baseline
eco_baselines <- trajectory_data %>%
  filter(Extent == "eco") %>%
  # Average across categories just in case, to get a single definitive rank per variable
  group_by(Set, Plot_Label) %>%
  summarise(Eco_Score = mean(Mean_Importance, na.rm = TRUE), .groups = "drop")

# 2. Flag the Global Niche variables so we can force them to the top
eco_baselines <- eco_baselines %>%
  mutate(Is_Global = grepl("global|niche", tolower(Plot_Label)))

# 3. Sort the dataframe to create our perfect factor levels
# We sort by Set first. 
# Then Is_Global (FALSE comes before TRUE, placing Global at the end of the factor list = TOP of the plot)
# Then Eco_Score (Ascending, so the highest score is near the end = near the TOP of the plot)
sorted_levels <- eco_baselines %>%
  arrange(Set, Is_Global, Eco_Score) %>%
  pull(Plot_Label)

# 4. Apply these new strict levels back to the main dataset
trajectory_data <- trajectory_data %>%
  mutate(Plot_Label = factor(Plot_Label, levels = sorted_levels))


#reorder trajectory data

trajectory_data <- trajectory_data %>%
  mutate(Category = factor(Category, levels = c("Iberian Endemic", "Native Widespread", "Invasive Widespread")))


# ==============================================================================

# X-axis is Importance, Y-axis is our newly ordered unique Plot_Label

# ==============================================================================

shift_plot <- ggplot(trajectory_data, aes(y = Plot_Label, x = Mean_Importance, color = Extent, group = Extent)) +
  
  
  
  geom_errorbar(aes(xmin = pmax(0, Mean_Importance - CI_95), 
                    
                    xmax = Mean_Importance + CI_95), 
                
                position = position_dodge(width = 0.7), width = 0.3, alpha = 0.6, linewidth = 0.8) +
  
  
  
  geom_point(position = position_dodge(width = 0.7), size = 2.5) +
  
  
  
  # CRITICAL FIX: scales = "free_y" and space = "free_y" drops the variables that don't belong to the Set
  
  facet_grid(Set ~ Category, scales = "free_y", space = "free_y") +
  
  
  
  scale_x_continuous(labels = function(x) paste0(round(x * 100, 0), "%")) +
  
  
  
  # CRITICAL FIX: This strips away the "__Set" tag we added for ordering, leaving just the clean name
  
  scale_y_discrete(labels = function(x) gsub("__.*", "", x)) +
  

  
  scale_color_manual(
    
    values = c( "eco" = "#3F007D",  # Deep Royal Purple (Stronger baseline)
                "H5"  = "#6A51A3",  # Rich Purple
                "H8"  = "#807DBA",  # Medium Purple (Your old H5)
                "H12" = "#9E9AC8"),
    
    labels = c(
      
      "eco" = expression(bold("eco")), 
      
      "H5" = expression(bold(H[5])), 
      
      "H8" = expression(bold(H[8])), 
      
      "H12" = expression(bold(H[12]))
      
    )
    
  ) +
  
  
  
  theme_minimal() +
  
  labs(y = NULL, x = "Mean Relative Importance (%)", color = "Spatial Extent") +
  
  theme(
    
    strip.text.x = element_text(face = "bold", size = 14),
    
    strip.text.y = element_text(face = "bold", size = 14), 
    
    axis.text.y = element_text(size = 14, color = "black", face = "bold"), 
    
    axis.text.x = element_text(size = 14, color = "black"),
    
    axis.title.x = element_text(size = 18, face = "bold"),
    
    panel.grid.minor = element_blank(),
    
    legend.position = "bottom",
    
    legend.title = element_text(face = "bold", size = 16),
    
    legend.text = element_text(size = 16) 
    
  )


# ==============================================================================

# 3. STITCH EVERYTHING TOGETHER WITH PATCHWORK AND ADD BOXES

# ==============================================================================


# Create the map row and apply margins to all maps simultaneously using '&'

map_row <- (
  
  (wrap_elements(map_plots[["eco"]] + labs(title = NULL)) + 
     
     ggtitle(expression(bold("eco"))) + 
     
     theme(plot.title = element_text(face = "bold", hjust = 0.5))) | 
    
    
    
    (wrap_elements(map_plots[["H5"]] + labs(title = NULL)) + 
       
       ggtitle(expression(bold(H[5]))) + 
       
       theme(plot.title = element_text(face = "bold", hjust = 0.5))) | 
    
    
    
    (wrap_elements(map_plots[["H8"]] + labs(title = NULL)) + 
       
       ggtitle(expression(bold(H[8]))) + 
       
       theme(plot.title = element_text(face = "bold", hjust = 0.5))) | 
    
    
    
    (wrap_elements(map_plots[["H12"]] + labs(title = NULL)) + 
       
       ggtitle(expression(bold(H[12]))) + 
       
       theme(plot.title = element_text(face = "bold", hjust = 0.5)))
  
) & theme(plot.margin = margin(t = 5, r = 20, b = 5, l = 20, unit = "pt"))


# Compose the plot with patchwork

composed_plot <- shift_plot / map_row + 
  
  plot_layout(heights = c(3.5, 1))


# Convert to a drawing canvas and apply the transparent boxes

master_plot_with_boxes <- ggdraw(composed_plot) +
  
  
  
  # Add Climate Box (Orange - top row)
  
  annotate("rect", xmin = 0.19, xmax = 0.965, ymin = 0.75, ymax = 0.96, 
           
           color = "black", fill = NA, linewidth = 0.5) +
  
  
  
  # Add Hydroclimatic Box (Blue - middle row)
  
  annotate("rect", xmin = 0.19, xmax = 0.965, ymin = 0.54, ymax = 0.74, 
           
           color = "black", fill = NA, linewidth = 0.5) +
  
  
  
  # Add Hydromorphology Box (Green - bottom row)
  
  annotate("rect", xmin = 0.19, xmax = 0.965, ymin = 0.325, ymax = 0.53, 
           
           color = "black", fill = NA, linewidth = 0.5)


# print(master_plot_with_boxes)


# Save the final plot

ggsave("Vagenas_aSDMs/output/figures/Regional/VarImp/Figure3_Variables_CI.png", 
       
       master_plot_with_boxes, width = 20, height = 12, dpi = 300) 





# ================================================================================================================
#### TRAINING EXTENT ~ PERFORMANCE METRICS | IMPORTANT SUPPLEMENTARY MATERIAL ####
# ================================================================================================================





# ==============================================================================
# 1. LOAD AND PREPARE REGIONAL DATA
# ==============================================================================
# Load all evaluation files dynamically
regional_files <- list.files("Vagenas_aSDMs/output/regional/", 
                             pattern = "eval_regional_iberia_strict.csv", 
                             recursive = TRUE, full.names = TRUE)

regional_metrics <- lapply(regional_files, function(f) {
  df <- read_csv(f, show_col_types = FALSE)
  parts <- unlist(strsplit(f, "/"))
  df$Species <- parts[length(parts) - 3]
  df$Extent  <- parts[length(parts) - 2]
  df$Set     <- parts[length(parts) - 1]
  return(df)
}) %>% bind_rows()

# Create the plot data by filtering for the 3 sets and pivoting
plot_data_reg <- regional_metrics %>%
  filter(Set %in% c("Climate", "Hydroclimatic", "Hydromorphological")) %>%
  dplyr::select(Species, Extent, Set, e_AUC, CBI, maxTSS, uAUC) %>%
  rename(AUC = e_AUC) %>% 
  pivot_longer(cols = c(AUC, CBI, maxTSS, uAUC), names_to = "Metric", values_to = "Score") %>%
  filter(!is.na(Score))

# Lock in the order of the factors so the plot is perfectly organized
plot_data_reg$Extent <- factor(plot_data_reg$Extent, levels = c("eco", "H5", "H8", "H12"))
plot_data_reg$Set <- factor(plot_data_reg$Set, levels = c("Climate", "Hydroclimatic", "Hydromorphological"))


# ==============================================================================
# 2. CALCULATE PAIRWISE STATS AND DIRECTIONAL LABELS
# ==============================================================================
# Load your species categories
species_class <- read_csv("Vagenas_aSDMs/input/aSDMs_species_classification.csv", show_col_types = FALSE) %>%
  mutate(Species = gsub(" ", "_", Sp)) %>%
  dplyr::select(Species, Category)

# Join the categories directly to plot_data_reg 
plot_data_reg <- plot_data_reg %>%
  left_join(species_class, by = "Species") %>%
  mutate(Category = ifelse(is.na(Category), "Unknown", Category))

# Define the exact comparisons and adjacent pairs for the bracket logic
my_comparisons <- list(
  c("eco", "H5"), c("eco", "H8"), c("eco", "H12"),
  c("H5", "H8"), c("H5", "H12"), c("H8", "H12")
)
adjacent_pairs <- c("eco-H5", "H5-H8", "H8-H12")

# Define the loops: First pooled, then the three specific categories
categories_to_plot <- c("Pooled", "Iberian Endemic", "Native Widespread", "Invasive Widespread")


#Extent_to_Performance_Important_Supplementary_Material

extperfdir<-"/Users/geo_v/Desktop/Vagenas_aSDMs/output/figures/Regional/Extent_Performance/"

if (!dir.exists(extperfdir)) {
  dir.create(extperfdir, recursive = TRUE)
}


# ==============================================================================
# 1. LOAD CLASSIFICATION AND FIX THE DATAFRAME
# ==============================================================================
# Load species classification and fix the spaces
species_class <- read_csv("Vagenas_aSDMs/input/aSDMs_species_classification.csv", show_col_types = FALSE) %>%
  mutate(Species = gsub(" ", "_", Sp)) %>%
  dplyr::select(Species, Category)

# Rebuild plot_data_reg with Category and fix NAs in one clean pipe
plot_data_reg <- regional_metrics %>%
  inner_join(species_class, by = "Species") %>%  
  mutate(Category = ifelse(is.na(Category), "Unknown", Category)) %>% # <-- Added right after the join
  filter(Set %in% c("Climate", "Hydroclimatic", "Hydromorphological")) %>%
  dplyr::select(Species, Category, Extent, Set, e_AUC, CBI, maxTSS, uAUC) %>%
  rename(AUC = e_AUC) %>%
  pivot_longer(cols = c(AUC, CBI, maxTSS, uAUC), names_to = "Metric", values_to = "Score") %>%
  filter(!is.na(Score))

# Lock in the order of the factors
plot_data_reg$Extent <- factor(plot_data_reg$Extent, levels = c("eco", "H5", "H8", "H12"))
plot_data_reg$Set <- factor(plot_data_reg$Set, levels = c("Climate", "Hydroclimatic", "Hydromorphological"))
plot_data_reg$Category <- factor(plot_data_reg$Category, levels = c("Iberian Endemic", "Native Widespread", "Invasive Widespread"))

# ==============================================================================
# 2. DEFINE COMPARISONS FOR THE LOOP
# ==============================================================================
my_comparisons <- list(
  c("eco", "H5"), c("eco", "H8"), c("eco", "H12"),
  c("H5", "H8"), c("H5", "H12"), c("H8", "H12")
)
adjacent_pairs <- c("eco-H5", "H5-H8", "H8-H12")

categories_to_plot <- c("Pooled", "Iberian Endemic", "Native Widespread", "Invasive Widespread")

# ==============================================================================
# 3. RUN THE PLOTTING LOOP
# ==============================================================================
for (cat_name in categories_to_plot) {
  
  cat(sprintf("\n======================================================\n"))
  cat(sprintf("Generating plot for: %s\n", cat_name))
  cat(sprintf("======================================================\n"))
  
  # A. Filter data and set titles dynamically
  if (cat_name == "Pooled") {
    current_data <- plot_data_reg
    title_text <- "Performance Across Spatial Extents (All Species Pooled) - [Regional aSDMs]"
    file_name <- "Extent_Alters_Performance_Regional_Pooled.png"
  } else {
    current_data <- plot_data_reg %>% filter(Category == cat_name)
    title_text <- sprintf("Performance Across Spatial Extents (%s) - [Regional aSDMs]", cat_name)
    file_name <- sprintf("Extent_Alters_Performance_Regional_%s.png", gsub(" ", "_", cat_name))
  }
  
  # B. Calculate Bounds and Stats for current subset
  facet_bounds <- current_data %>%
    group_by(Metric, Set) %>%
    summarise(
      min_val = min(Score, na.rm = TRUE),
      max_val = max(Score, na.rm = TRUE),
      .groups = "drop"
    )
  
  means_df <- current_data %>%
    group_by(Metric, Set, Extent) %>%
    summarise(Mean_Score = mean(Score, na.rm = TRUE), .groups = "drop")
  
  stat.test <- current_data %>%
    group_by(Metric, Set) %>%
    wilcox_test(Score ~ Extent, comparisons = my_comparisons) %>%
    add_significance() %>%
    
    left_join(means_df, by = c("Metric", "Set", "group1" = "Extent")) %>%
    rename(mean1 = Mean_Score) %>%
    left_join(means_df, by = c("Metric", "Set", "group2" = "Extent")) %>%
    rename(mean2 = Mean_Score) %>%
    
    mutate(
      direction = ifelse(mean1 > mean2, ">", "<"),
      custom_label = ifelse(p.adj.signif == "ns", "ns", paste0(p.adj.signif, " (", direction, ")")),
      pair_name = paste(group1, group2, sep = "-"),
      is_adjacent = pair_name %in% adjacent_pairs
    ) %>%
    
    filter(p.adj.signif != "ns") 
  
  # Initialize empty dataframes to prevent ggplot errors if no significance is found
  stat.test.top <- data.frame()
  stat.test.bot <- data.frame()
  
  # Only calculate bracket positioning if there are actually significant results
  if (nrow(stat.test) > 0) {
    stat.test <- stat.test %>%
      left_join(facet_bounds, by = c("Metric", "Set")) %>%
      group_by(Metric, Set, is_adjacent) %>%
      mutate(step_rank = row_number()) %>% 
      ungroup() %>%
      mutate(
        range = max_val - min_val,
        step_size = ifelse(range == 0, 0.05, range * 0.08), 
        y.position = ifelse(is_adjacent,
                            min_val - (step_size * 0.5) - (step_rank * step_size),
                            max_val + (step_size * 0.5) + (step_rank * step_size))
      )
    
    stat.test.top <- stat.test %>% filter(!is_adjacent)
    stat.test.bot <- stat.test %>% filter(is_adjacent)
  }
  
  # C. Build the base plot
  
  # Define your reversed sequential purple gradient for Extents
  extent_colors <- c( "eco" = "#3F007D",  # Deep Royal Purple (Stronger baseline)
                      "H5"  = "#6A51A3",  # Rich Purple
                      "H8"  = "#807DBA",  # Medium Purple (Your old H5)
                      "H12" = "#9E9AC8"
  )
  
  p <- ggboxplot(current_data, x = "Extent", y = "Score", 
                 color = "Extent", palette = extent_colors,
                 facet.by = c("Metric", "Set"), scales = "free_y",
                 short.panel.labs = FALSE) +
    stat_compare_means(method = "kruskal.test", label.y.npc = "bottom")
  
  # D. Conditionally add brackets
  if (nrow(stat.test.top) > 0) {
    p <- p + stat_pvalue_manual(stat.test.top, label = "custom_label", tip.length = 0.01)
  }
  if (nrow(stat.test.bot) > 0) {
    p <- p + stat_pvalue_manual(stat.test.bot, label = "custom_label", tip.length = -0.01, vjust = 1.5)
  }
  
  # E. Apply formatting
  p <- p + 
    scale_y_continuous(
      labels = function(x) ifelse(is.na(x), "", ifelse(x > 1.0, "", x)),
      expand = expansion(mult = c(0.18, 0.18)) 
    ) +
    theme_minimal() +
    labs(title = title_text,
         subtitle = "Adjacent pairs (bottom brackets) vs Long-jump pairs (top brackets). '>' indicates Left Extent scored higher.",
         x = "Spatial Training Extent", 
         y = "Metric Score",
         caption = "Significance: * (p ≤ 0.05), ** (p ≤ 0.01), *** (p ≤ 0.001), **** (p ≤ 0.0001)") +
    theme(plot.title = element_text(face = "bold", size = 16),
          plot.caption = element_text(hjust = 0, size = 10, face = "italic", color = "gray30"),
          legend.position = "none",
          panel.spacing = unit(1, "lines"))
  
  # F. Print and Save
  print(p)
  
  save_path <- sprintf("Vagenas_aSDMs/output/figures/Regional/Extent_Performance/%s", file_name)
  ggsave(save_path, p, width = 16, height = 12, dpi = 300)
  
  cat(sprintf("Saved: %s\n", save_path))
}



# ================================================================================================================
#### PREDICTOR SETS ~ PERFORMANCE METRICS | IMPORTANT SUPPLEMENTARY MATERIAL ####
# ================================================================================================================


#Critical supplementary material too to connect with Climate - Hydroclimatic - Hydromorphology 

predperfdir<-"/Users/geo_v/Desktop/Vagenas_aSDMs/output/figures/Regional/Predictor_Performance/"

if (!dir.exists(predperfdir)) {
  dir.create(predperfdir, recursive = TRUE)
}


# ==============================================================================
# 1. LOAD CLASSIFICATION AND FIX THE DATAFRAME
# ==============================================================================
# Load species classification and fix the spaces
species_class <- read_csv("Vagenas_aSDMs/input/aSDMs_species_classification.csv", show_col_types = FALSE) %>%
  mutate(Species = gsub(" ", "_", Sp)) %>%
  dplyr::select(Species, Category)

# Rebuild plot_data_reg with Category and fix NAs in one clean pipe
plot_data_reg <- regional_metrics %>%
  inner_join(species_class, by = "Species") %>%  
  mutate(Category = ifelse(is.na(Category), "Unknown", Category)) %>% 
  filter(Set %in% c("Climate", "Hydroclimatic", "Hydromorphological")) %>%
  dplyr::select(Species, Category, Extent, Set, e_AUC, CBI, maxTSS, uAUC) %>%
  rename(AUC = e_AUC) %>%
  pivot_longer(cols = c(AUC, CBI, maxTSS, uAUC), names_to = "Metric", values_to = "Score") %>%
  filter(!is.na(Score))

# Lock in the order of the factors
plot_data_reg$Extent <- factor(plot_data_reg$Extent, levels = c("eco", "H5", "H8", "H12"))
plot_data_reg$Set <- factor(plot_data_reg$Set, levels = c("Climate", "Hydroclimatic", "Hydromorphological"))
plot_data_reg$Category <- factor(plot_data_reg$Category, levels = c("Iberian Endemic", "Native Widespread", "Invasive Widespread"))

# ==============================================================================
# 2. DEFINE COMPARISONS FOR THE LOOP (MODIFIED FOR 'SET')
# ==============================================================================
# Comparing Sets instead of Extents
my_comparisons <- list(
  c("Climate", "Hydroclimatic"), 
  c("Climate", "Hydromorphological"), 
  c("Hydroclimatic", "Hydromorphological")
)

# Define which pairs are adjacent on the x-axis
adjacent_pairs <- c("Climate-Hydroclimatic", "Hydroclimatic-Hydromorphological")

categories_to_plot <- c("Pooled", "Iberian Endemic", "Native Widespread", "Invasive Widespread")

# ==============================================================================
# 3. RUN THE PLOTTING LOOP
# ==============================================================================
for (cat_name in categories_to_plot) {
  
  cat(sprintf("\n======================================================\n"))
  cat(sprintf("Generating plot for: %s\n", cat_name))
  cat(sprintf("======================================================\n"))
  
  # A. Filter data and set titles dynamically
  if (cat_name == "Pooled") {
    current_data <- plot_data_reg
    title_text <- "Performance Across Predictor Sets (All Species Pooled) - [Regional aSDMs]"
    file_name <- "Set_Alters_Performance_Regional_Pooled.png"
  } else {
    current_data <- plot_data_reg %>% filter(Category == cat_name)
    title_text <- sprintf("Performance Across Predictor Sets (%s) - [Regional aSDMs]", cat_name)
    file_name <- sprintf("Set_Alters_Performance_Regional_%s.png", gsub(" ", "_", cat_name))
  }
  
  # B. Calculate Bounds and Stats for current subset
  # Swapped grouping to Metric and Extent
  facet_bounds <- current_data %>%
    group_by(Metric, Extent) %>%
    summarise(
      min_val = min(Score, na.rm = TRUE),
      max_val = max(Score, na.rm = TRUE),
      .groups = "drop"
    )
  
  means_df <- current_data %>%
    group_by(Metric, Extent, Set) %>%
    summarise(Mean_Score = mean(Score, na.rm = TRUE), .groups = "drop")
  
  # Wilcox test across 'Set', grouped by 'Metric' and 'Extent'
  stat.test <- current_data %>%
    group_by(Metric, Extent) %>%
    wilcox_test(Score ~ Set, comparisons = my_comparisons) %>%
    add_significance() %>%
    
    # Join means to figure out the direction of the difference
    left_join(means_df, by = c("Metric", "Extent", "group1" = "Set")) %>%
    rename(mean1 = Mean_Score) %>%
    left_join(means_df, by = c("Metric", "Extent", "group2" = "Set")) %>%
    rename(mean2 = Mean_Score) %>%
    
    mutate(
      direction = ifelse(mean1 > mean2, ">", "<"),
      custom_label = ifelse(p.adj.signif == "ns", "ns", paste0(p.adj.signif, " (", direction, ")")),
      pair_name = paste(group1, group2, sep = "-"),
      is_adjacent = pair_name %in% adjacent_pairs
    ) %>%
    
    filter(p.adj.signif != "ns") 
  
  # Initialize empty dataframes to prevent ggplot errors if no significance is found
  stat.test.top <- data.frame()
  stat.test.bot <- data.frame()
  
  # Only calculate bracket positioning if there are actually significant results
  if (nrow(stat.test) > 0) {
    stat.test <- stat.test %>%
      left_join(facet_bounds, by = c("Metric", "Extent")) %>%
      group_by(Metric, Extent, is_adjacent) %>%
      mutate(step_rank = row_number()) %>% 
      ungroup() %>%
      mutate(
        range = max_val - min_val,
        step_size = ifelse(range == 0, 0.05, range * 0.08), 
        y.position = ifelse(is_adjacent,
                            min_val - (step_size * 0.5) - (step_rank * step_size),
                            max_val + (step_size * 0.5) + (step_rank * step_size))
      )
    
    stat.test.top <- stat.test %>% filter(!is_adjacent)
    stat.test.bot <- stat.test %>% filter(is_adjacent)
  }
  
  # C. Build the base plot
  # x mapped to Set, faceted by Metric and Extent
  
  # Define your exact colors for the Predictor Sets
  set_colors <- c(
    "Climate" = "#E69F00", 
    "Hydroclimatic" = "#56B4E9", 
    "Hydromorphological" = "#009E73"
  )
  
  
  p <- ggboxplot(current_data, x = "Set", y = "Score", 
                 color = "Set", palette = set_colors,
                 facet.by = c("Metric", "Extent"), scales = "free_y",
                 short.panel.labs = FALSE) +
    stat_compare_means(method = "kruskal.test", label.y.npc = "bottom")
  
  # D. Conditionally add brackets
  if (nrow(stat.test.top) > 0) {
    p <- p + stat_pvalue_manual(stat.test.top, label = "custom_label", tip.length = 0.01)
  }
  if (nrow(stat.test.bot) > 0) {
    p <- p + stat_pvalue_manual(stat.test.bot, label = "custom_label", tip.length = -0.01, vjust = 1.5)
  }
  
  # E. Apply formatting
  p <- p + 
    scale_y_continuous(
      labels = function(x) ifelse(is.na(x), "", ifelse(x > 1.0, "", x)),
      expand = expansion(mult = c(0.18, 0.18)) 
    ) +
    theme_minimal() +
    labs(title = title_text,
         subtitle = "Adjacent pairs (bottom brackets) vs Long-jump pairs (top brackets). '>' indicates Left Set scored higher.",
         x = "Predictor Set", 
         y = "Metric Score",
         caption = "Significance: * (p ≤ 0.05), ** (p ≤ 0.01), *** (p ≤ 0.001), **** (p ≤ 0.0001)") +
    theme(plot.title = element_text(face = "bold", size = 16),
          plot.caption = element_text(hjust = 0, size = 10, face = "italic", color = "gray30"),
          legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1), # Added angle for longer 'Set' labels
          panel.spacing = unit(1, "lines"))
  
  # F. Print and Save
  print(p)
  
  # Make sure you create this directory if it doesn't exist yet, or adjust the path!
  save_path <- sprintf("Vagenas_aSDMs/output/figures/Regional/Predictor_Performance/%s", file_name)
  ggsave(save_path, p, width = 16, height = 12, dpi = 300)
  
  cat(sprintf("Saved: %s\n", save_path))
}









# ================================================================================================================
#### ENSEMBLES - FINAL PRODUCTS | FIGURE 4 ####
# ================================================================================================================




setwd("/Users/geo_v/Desktop")

library(terra)

# Directories
base_dir <- "Vagenas_aSDMs/output/regional"
# New output folder to keep this multi-metric run separate from the old one
out_dir <- "Vagenas_aSDMs/output/figures/Regional/ensembles/regional_ensembles_H5_AUC" 

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

species_list <- list.dirs(base_dir, recursive = FALSE, full.names = FALSE)
sets <- c("Climate", "Hydroclimatic", "Hydromorphological")

weights_list <- list()

for (sp in species_list) {
  
  # A small dataframe to hold the 4 metrics for the 3 sets
  sp_metrics <- data.frame(
    Set = sets,
    e_AUC = numeric(3)
  )
  
  sp_rasters <- list()
  
  # --- STEP A: LOAD DATA & METRICS ---
  
  #this is when all species are fully structured with climate, hydroclimatic, hydromorphological
  
  # for (i in seq_along(sets)) {
  #   set <- sets[i]
  #   eval_file <- file.path(base_dir, sp, "H5", set, "eval_regional_iberia_strict.csv")
  #   raster_file <- file.path(base_dir, sp, "H5", set, "ensemble_regional.tif")
  #   
  #   # Extract ALL FOUR metrics but in this case only for AUC
  #   eval_data <- read_csv(eval_file, show_col_types = FALSE)
  #   sp_metrics$e_AUC[i]  <- eval_data$e_AUC[1]
  #   # sp_metrics$CBI[i]    <- eval_data$CBI[1]
  #   # sp_metrics$maxTSS[i] <- eval_data$maxTSS[1]
  #   # sp_metrics$uAUC[i]   <- eval_data$uAUC[1]
  #   
  #   sp_rasters[[set]] <- rast(raster_file)
  # }
  
  #modification for Lamperta planeri
  
  for (i in seq_along(sets)) {
    set <- sets[i]
    eval_file <- file.path(base_dir, sp, "H5", set, "eval_regional_iberia_strict.csv")
    raster_file <- file.path(base_dir, sp, "H5", set, "ensemble_regional.tif")
    
    # ---- Specific fix for Lampetra_planeri ----
    if (sp == "Lampetra_planeri" && set == "Hydromorphological") {
      sp_metrics$e_AUC[i] <- NA           # mark as missing
      # do not load the raster – leave sp_rasters[[set]] empty
      next
    }
    
    eval_data <- read_csv(eval_file, show_col_types = FALSE)
    sp_metrics$e_AUC[i]  <- eval_data$e_AUC[1]
    sp_rasters[[set]] <- rast(raster_file)
  }
  
  # --- STEP B: NORMALIZE INDIVIDUAL RASTERS (0 TO 1) ---
  for (set in names(sp_rasters)) {
    r <- sp_rasters[[set]]
    r_minmax <- minmax(r)
    r_min <- r_minmax[1, 1]
    r_max <- r_minmax[2, 1]
    
    if (r_max > r_min) {
      sp_rasters[[set]] <- (r - r_min) / (r_max - r_min)
    } else {
      sp_rasters[[set]] <- r - r_min 
    }
  }
  
  # --- STEP C: CALCULATE MULTI-METRIC WEIGHTS ---
  # 1. Calculate proportional weights (sum to 1) for each metric individually
  
  #this is when all species are fully structured with climate, hydroclimatic, hydromorphological
  
  # #In this case only for AUC
  # w_AUC    <- sp_metrics$e_AUC / sum(sp_metrics$e_AUC)
  # # w_CBI    <- sp_metrics$CBI / sum(sp_metrics$CBI)
  # # w_maxTSS <- sp_metrics$maxTSS / sum(sp_metrics$maxTSS)
  # # w_uAUC   <- sp_metrics$uAUC / sum(sp_metrics$uAUC)
  # 
  # # 2. Average the 4 proportions to enforce exactly 1/4 weight per metric
  # #In this case only for AUC
  # #final_weights <- (w_AUC + w_CBI + w_maxTSS + w_uAUC) / 4
  # final_weights<-w_AUC
  # names(final_weights) <- sets
  
  # modification for Lampetra_planeri
  if (sp == "Lampetra_planeri") {
    valid_sets <- sp_metrics$Set[!is.na(sp_metrics$e_AUC)]   # "Climate","Hydroclimatic"
    valid_auc  <- sp_metrics$e_AUC[!is.na(sp_metrics$e_AUC)]
    w_AUC <- valid_auc / sum(valid_auc)
    final_weights <- w_AUC
    names(final_weights) <- valid_sets
  } else {
    # All other species keep the original 3‑set calculation
    w_AUC    <- sp_metrics$e_AUC / sum(sp_metrics$e_AUC)
    final_weights <- w_AUC
    names(final_weights) <- sets
  }
  
  
  #this is when all species are fully structured with climate, hydroclimatic, hydromorphological
  
  # # Store percentages for the CSV
  # weights_list[[sp]] <- data.frame(
  #   Species = sp,
  #   Climate = final_weights["Climate"] * 100,
  #   Hydroclimatic = final_weights["Hydroclimatic"] * 100,
  #   Hydromorphology = final_weights["Hydromorphological"] * 100 
  # )
  
  #modification for Lamperta planeri
  
  if (sp == "Lampetra_planeri") {
    weights_list[[sp]] <- data.frame(
      Species = sp,
      Climate = final_weights["Climate"] * 100,
      Hydroclimatic = final_weights["Hydroclimatic"] * 100,
      Hydromorphology = NA   # explicitly missing
    )
  } else {
    weights_list[[sp]] <- data.frame(
      Species = sp,
      Climate = final_weights["Climate"] * 100,
      Hydroclimatic = final_weights["Hydroclimatic"] * 100,
      Hydromorphology = final_weights["Hydromorphological"] * 100 
    )
  }
  
  # --- STEP D: RASTER MAP ALGEBRA ---
  
  #this is when all species are fully structured with climate, hydroclimatic, hydromorphological
  
  # weighted_ensemble <- (sp_rasters[["Climate"]] * final_weights["Climate"]) +
  #   (sp_rasters[["Hydroclimatic"]] * final_weights["Hydroclimatic"]) +
  #   (sp_rasters[["Hydromorphological"]] * final_weights["Hydromorphological"])
  
  #modification for Lamperta planeri
  
  if (sp == "Lampetra_planeri") {
    # Only two sets exist
    weighted_ensemble <- (sp_rasters[["Climate"]] * final_weights["Climate"]) +
      (sp_rasters[["Hydroclimatic"]] * final_weights["Hydroclimatic"])
  } else {
    weighted_ensemble <- (sp_rasters[["Climate"]] * final_weights["Climate"]) +
      (sp_rasters[["Hydroclimatic"]] * final_weights["Hydroclimatic"]) +
      (sp_rasters[["Hydromorphological"]] * final_weights["Hydromorphological"])
  }
  
  # --- STEP E: FINAL NORMALIZATION 0 TO 1 ---
  ens_minmax <- minmax(weighted_ensemble)
  ens_min <- ens_minmax[1, 1]
  ens_max <- ens_minmax[2, 1]
  
  if (ens_max > ens_min) {
    normalized_ensemble <- (weighted_ensemble - ens_min) / (ens_max - ens_min)
  } else {
    normalized_ensemble <- weighted_ensemble
  }
  
  # --- STEP F: EXPORT TO THE NEW FOLDER ---
  out_raster_path <- file.path(out_dir, paste0(sp, "_H5_AUC_ensemble.tif"))
  writeRaster(normalized_ensemble, out_raster_path, overwrite = TRUE)
  
  message("Successfully processed: ", sp)
}

# --- STEP G: COMPILE AND EXPORT CSV ---
final_weights_df <- bind_rows(weights_list)

write_csv(final_weights_df, file.path(out_dir, "H5_AUC_influence_percentages.csv"))

print("Process complete! All multimetric files safely stored in 'regional_ensembles_H5_AUC'.")





#Lets begin stacking the regional setting


# ==============================================================================
# 0. SETUP AND LOAD PACKAGES
# ==============================================================================

library(tidyterra)
library(patchwork)
library(ggcorrplot)
library(rnaturalearth)

setwd("/Users/geo_v/Desktop/")

# Directories
ens_dir <- "Vagenas_aSDMs/output/figures/Regional/ensembles/regional_ensembles_H5_AUC"
out_dir <- "Vagenas_aSDMs/output/figures/Regional/ensembles/stacked_maps_AUC"

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# ==============================================================================
# 1. PREPARE DATA 
# ==============================================================================
species_class <- read_csv("Vagenas_aSDMs/input/aSDMs_species_classification.csv", show_col_types = FALSE) %>%
  mutate(Species = gsub(" ", "_", Sp)) %>%
  dplyr::select(Species, Category)

influence_df <- read_csv(file.path(ens_dir, "H5_AUC_influence_percentages.csv"), show_col_types = FALSE)

master_df <- inner_join(influence_df, species_class, by = "Species")

# ==============================================================================
# 1.5 FETCH, MERGE, AND CLEAN IBERIAN MAINLAND
# ==============================================================================
iberia_sf <- ne_countries(country = c("Spain", "Portugal"), scale = "medium", returnclass = "sf")
iberia_vect <- vect(iberia_sf)

iberia_dissolved <- aggregate(iberia_vect)
iberia_parts <- disagg(iberia_dissolved)
land_areas <- expanse(iberia_parts, unit = "km")
iberia_mainland <- iberia_parts[which.max(land_areas), ]

# ==============================================================================
# 2. CREATE STACKED RASTERS (CROPPED & NORMALIZED TO MAINLAND)
# ==============================================================================
safe_sum_rasters <- function(file_paths) {
  r_list <- lapply(file_paths, rast)
  
  master_ext <- ext(r_list[[1]])
  if (length(r_list) > 1) {
    for (i in 2:length(r_list)) {
      master_ext <- terra::union(master_ext, ext(r_list[[i]]))
    }
  }
  
  r_list_aligned <- lapply(r_list, function(r) extend(r, master_ext))
  r_stack <- rast(r_list_aligned)
  r_sum <- sum(r_stack, na.rm = TRUE)
  
  iberia_proj <- project(iberia_mainland, crs(r_sum))
  r_crop <- crop(r_sum, iberia_proj)
  r_masked <- mask(r_crop, iberia_proj)
  
  r_minmax <- minmax(r_masked)
  r_min <- r_minmax[1, 1]
  r_max <- r_minmax[2, 1]
  
  if (r_max > r_min) {
    r_normalized <- ((r_masked - r_min) / (r_max - r_min)) * 100
  } else {
    r_normalized <- r_masked * 0 
  }
  
  return(r_normalized)
}

categories <- c("Iberian Endemic", "Invasive Widespread", "Native Widespread")
stacked_rasters <- list()

for (cat in categories) {
  spp_in_cat <- master_df %>% filter(Category == cat) %>% pull(Species)
  files <- file.path(ens_dir, paste0(spp_in_cat, "_H5_AUC_ensemble.tif"))
  files <- files[file.exists(files)] 
  
  if(length(files) > 0) {
    stacked_rasters[[cat]] <- safe_sum_rasters(files)
    file_name_safe <- gsub(" ", "_", cat)
    writeRaster(stacked_rasters[[cat]], file.path(out_dir, paste0(file_name_safe, "_stacked_suitability.tif")), overwrite = TRUE)
  }
}

all_files <- list.files(ens_dir, pattern = "_H5_AUC_ensemble\\.tif$", full.names = TRUE)
stacked_rasters[["Pooled"]] <- safe_sum_rasters(all_files)
writeRaster(stacked_rasters[["Pooled"]], file.path(out_dir, "Pooled_stacked_suitability_AUC.tif"), overwrite = TRUE)


# ==============================================================================
#3. PLOT THE STACKED MAPS
# ==============================================================================
yellow_to_very_dark_orange <- colorRampPalette(c("yellow", "#993300"))(100)

plot_map <- function(r, title) {
  ggplot() +
    geom_spatraster(data = r) +
    geom_spatvector(data = iberia_mainland, fill = NA, color = "black", linewidth = 0.4) +
    scale_fill_gradientn(
      colors = yellow_to_very_dark_orange,
      name = "Stacked Suitability Index",
      limits = c(0, 100),                           
      breaks = seq(0, 100, by = 20),                
      labels = paste0(seq(0, 100, by = 20), "%"),   
      na.value = "transparent",
      guide = guide_colorbar(
        direction = "horizontal",
        label.position = "bottom",
        title.position = "top",
        title.hjust = 0.5,                          
        barwidth = unit(15, "cm"),   
        barheight = unit(1.5, "cm")    
      )
    ) +
    labs(title = title) + 
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 22),
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      plot.margin = margin(5, 5, 5, 5, unit = "pt") # Tight margins keep titles close
    )
}

n_endemic  <- nrow(master_df %>% filter(Category == "Iberian Endemic"))
n_native   <- nrow(master_df %>% filter(Category == "Native Widespread"))
n_invasive <- nrow(master_df %>% filter(Category == "Invasive Widespread"))
n_pooled   <- nrow(master_df)

p_endemic  <- plot_map(stacked_rasters[["Iberian Endemic"]], paste0("Iberian Endemic (N=", n_endemic, ")"))
p_native   <- plot_map(stacked_rasters[["Native Widespread"]], paste0("Native Widespread (N=", n_native, ")"))
p_invasive <- plot_map(stacked_rasters[["Invasive Widespread"]], paste0("Invasive Widespread (N=", n_invasive, ")"))
p_pooled   <- plot_map(stacked_rasters[["Pooled"]], paste0("All Species (N=", n_pooled, ")"))

# ==============================================================================
# 5. ARRANGE PROPORTIONAL LAYOUT WITH PATCHWORK (MATHEMATICAL GRID)
# ==============================================================================


layout_design <- c(
  patchwork::area(t = 1, l = 2, b = 1, r = 2), # p_endemic (Top Row, Middle Col)
  patchwork::area(t = 2, l = 1, b = 2, r = 1), # p_native (Middle Row, Left Col)
  patchwork::area(t = 2, l = 3, b = 2, r = 3), # p_invasive (Middle Row, Right Col)
  patchwork::area(t = 2, l = 2, b = 2, r = 2)  # p_pooled (Middle Row, Middle Col)
)

final_figure <- p_endemic + p_native + p_invasive + p_pooled +
  plot_layout(
    design = layout_design,
    widths = c(1, 2.2, 1),
    heights = c(1, 2.2),  # Only 2 rows now!
    guides = "collect"
  ) & 
  theme(
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 16),
    legend.text = element_text(size = 12),                 
    plot.background = element_rect(fill = "white", color = NA) 
  )

# Save the master figure
ggsave(file.path(out_dir, "Figure4_Groups_Ensembles_AUC.png"), 
       final_figure, width = 18, height = 14, dpi = 300, bg = "white")

#save pooled across species

writeRaster(stacked_rasters[["Pooled"]], file.path(out_dir, "Figure4_Groups_Ensembles_AUC.tif"), overwrite = TRUE)


#alternative 2.0 - Instead of grouping and ensembling per species, 
#I am grouping and ensembling per predictor category for all the species

# ==============================================================================
# 0. SETUP AND LOAD PACKAGES
# ==============================================================================

library(ggcorrplot)

# Directories
ens_dir <- "Vagenas_aSDMs/output/figures/Regional/ensembles/regional_ensembles_H5_AUC"
base_reg_dir <- "Vagenas_aSDMs/output/regional"


# ==============================================================================
# 1. PREPARE SPECIES LIST
# ==============================================================================
species_class <- read_csv("Vagenas_aSDMs/input/aSDMs_species_classification.csv", show_col_types = FALSE) %>%
  mutate(Species = gsub(" ", "_", Sp))

# Get the species that actually have valid outputs
valid_files <- list.files(ens_dir, pattern = "_H5_AUC_ensemble\\.tif$", full.names = TRUE)
valid_species <- gsub("_H5_AUC_ensemble\\.tif$", "", basename(valid_files))
species_list <- species_class %>% filter(Species %in% valid_species) %>% pull(Species)

# ==============================================================================
# 1.5 FETCH, MERGE, AND CLEAN IBERIAN MAINLAND
# ==============================================================================
iberia_sf <- ne_countries(country = c("Spain", "Portugal"), scale = "medium", returnclass = "sf")
iberia_vect <- vect(iberia_sf)

iberia_dissolved <- aggregate(iberia_vect)
iberia_parts <- disagg(iberia_dissolved)
land_areas <- expanse(iberia_parts, unit = "km")
iberia_mainland <- iberia_parts[which.max(land_areas), ]

# ==============================================================================
# 2. THE MULTIMETRIC STACKING FUNCTION (CORRECTED)
# ==============================================================================
build_multimetric_driver_stack <- function(spp_list, target_driver) {
  stack_list <- list()
  
  for (sp in spp_list) {
    file_path <- file.path(base_reg_dir, sp, "H5", target_driver, "ensemble_regional.tif")
    eval_path <- file.path(base_reg_dir, sp, "H5", target_driver, "eval_regional_iberia_strict.csv")
    
    if (file.exists(file_path) && file.exists(eval_path)) {
      
      # 1. Load the Evaluation Metrics and calculate the penalty weight
      metrics <- read_csv(eval_path, show_col_types = FALSE)
      weight <- mean(metrics$e_AUC + metrics$maxTSS + metrics$CBI + metrics$uAUC, na.rm = TRUE) / 4
      
      # 2. Load the RAW raster (Intentionally skipping the 0-1 stretch)
      r <- rast(file_path)
      
      # 3. Apply the dynamic weight directly to the raw probabilities
      r_weighted <- r * weight
      stack_list[[sp]] <- r_weighted
    }
  }
  
  if(length(stack_list) == 0) return(NULL) # Failsafe
  
  # 4. Align extents across all species
  master_ext <- ext(stack_list[[1]])
  if (length(stack_list) > 1) {
    for (i in 2:length(stack_list)) {
      master_ext <- terra::union(master_ext, ext(stack_list[[i]]))
    }
  }
  
  stack_aligned <- lapply(stack_list, function(r) extend(r, master_ext))
  r_stack <- rast(stack_aligned)
  
  # 5. Sum all weighted species together
  r_sum <- sum(r_stack, na.rm = TRUE)
  
  # 6. Crop and Mask to Iberia
  iberia_proj <- project(iberia_mainland, crs(r_sum))
  r_crop <- crop(r_sum, iberia_proj)
  r_masked <- mask(r_crop, iberia_proj)
  
  # RETURN RAW SUM (Intentionally skipping the 0-100% independent stretch)
  return(r_masked)
}

# --- GENERATE THE RAW STACKS ---
print("Building Raw Multimetric Climate Stack...")
r_climate_raw <- build_multimetric_driver_stack(species_list, "Climate")

print("Building Raw Multimetric Hydroclimatic Stack...")
r_hydrocl_raw <- build_multimetric_driver_stack(species_list, "Hydroclimatic")

print("Building Raw Multimetric Hydromorphological Stack...")
r_hydromo_raw <- build_multimetric_driver_stack(species_list, "Hydromorphological")


# --- APPLY GLOBAL STANDARDIZATION ---
print("Applying Global 0-100% Standardization across Predictor Sets...")

# Combine to find the absolute global min and max
global_stack <- c(r_climate_raw, r_hydrocl_raw, r_hydromo_raw)
global_minmax <- minmax(global_stack)

global_min <- min(global_minmax[1, ]) # The absolute lowest pixel anywhere
global_max <- max(global_minmax[2, ]) # The absolute highest pixel anywhere

# Scale each raster relative to the global bounds to preserve statistical differences
if (global_max > global_min) {
  r_climate <- ((r_climate_raw - global_min) / (global_max - global_min)) * 100
  r_hydrocl <- ((r_hydrocl_raw - global_min) / (global_max - global_min)) * 100
  r_hydromo <- ((r_hydromo_raw - global_min) / (global_max - global_min)) * 100
} else {
  r_climate <- r_climate_raw * 0
  r_hydrocl <- r_hydrocl_raw * 0
  r_hydromo <- r_hydromo_raw * 0
}

names(r_climate) <- "Climate"
names(r_hydrocl) <- "Hydroclimatic"
names(r_hydromo) <- "Hydromorphological"

# --- THE POOLED OVERALL MAP ---
# (We can independently stretch this one because it represents total cumulative suitability,
# not a direct intensity comparison against the three individual drivers)
safe_sum_overall <- function(files) {
  r_list <- lapply(files, rast)
  master_ext <- ext(r_list[[1]])
  if (length(r_list) > 1) {
    for (i in 2:length(r_list)) { master_ext <- terra::union(master_ext, ext(r_list[[i]])) }
  }
  r_list_aligned <- lapply(r_list, function(r) extend(r, master_ext))
  r_sum <- sum(rast(r_list_aligned), na.rm = TRUE)
  
  iberia_proj <- project(iberia_mainland, crs(r_sum))
  r_masked <- mask(crop(r_sum, iberia_proj), iberia_proj)
  
  rmm <- minmax(r_masked)
  return( ((r_masked - rmm[1,1]) / (rmm[2,1] - rmm[1,1])) * 100 )
}

print("Stacking Final Overall Multimetric Ensembles...")
r_pooled <- safe_sum_overall(valid_files)
names(r_pooled) <- "Overall"

# ==============================================================================
# 3. PLOT THE STACKED MAPS
# ==============================================================================
yellow_to_very_dark_orange <- colorRampPalette(c("yellow", "#993300"))(100)

plot_map <- function(r, title) {
  ggplot() +
    geom_spatraster(data = r) +
    geom_spatvector(data = iberia_mainland, fill = NA, color = "black", linewidth = 0.4) +
    scale_fill_gradientn(
      colors = yellow_to_very_dark_orange,
      name = "Stacked Suitability Index",
      limits = c(0, 100),
      breaks = seq(0, 100, by = 20),
      labels = paste0(seq(0, 100, by = 20), "%"),
      na.value = "transparent",
      guide = guide_colorbar(
        direction = "horizontal",
        label.position = "bottom",
        title.position = "top",
        title.hjust = 0.5,
        barwidth = unit(12, "cm"),
        barheight = unit(1, "cm")
      )
    ) +
    labs(title = title) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      plot.margin = margin(5, 5, 5, 5, unit = "pt")
    )
}

n_species <- length(species_list)
# 
p_pooled   <- plot_map(r_pooled, paste0("Overall Stacked Suitability (N=", n_species, ")"))
p_climate  <- plot_map(r_climate, "Climate-Driven Suitability")
p_hydrocl  <- plot_map(r_hydrocl, "Hydroclimatic-Driven Suitability")
p_hydromo  <- plot_map(r_hydromo, "Hydromorphological-Driven Suitability")

# # ==============================================================================
# # 5. ARRANGE PROPORTIONAL CROSS LAYOUT (MATHEMATICAL GRID)
# # ==============================================================================
# layout_design <- c(
#   patchwork::area(t = 1, l = 2, b = 1, r = 2), # Top Row, Middle Col
#   patchwork::area(t = 2, l = 1, b = 2, r = 1), # Middle Row, Left Col
#   patchwork::area(t = 2, l = 3, b = 2, r = 3), # Middle Row, Right Col
#   patchwork::area(t = 2, l = 2, b = 2, r = 2), # Middle Row, Middle Col
#   patchwork::area(t = 3, l = 1, b = 3, r = 1)  # Bottom Row, Left Col
# )
# 
# final_figure <- p_climate + p_hydrocl + p_hydromo + p_pooled + p_corr +
#   plot_layout(
#     design = layout_design,
#     widths = c(1, 2.2, 1),
#     heights = c(1, 1.4, 0.8),
#     guides = "collect"
#   ) &
#   theme(
#     legend.position = "bottom",
#     legend.title = element_text(face = "bold", size = 16),
#     legend.text = element_text(size = 12),
#     plot.background = element_rect(fill = "white", color = NA),
#     legend.margin = margin(t = -10, unit = "pt")
#   )
# 
# ggsave(file.path(out_dir, "Final_Multimetric_Driver_Figure_AUC.png"),
#        final_figure, width = 18, height = 14, dpi = 300, bg = "white")
# 
# print("Complete! Maps successfully built from scratch using proper global scaling.")
# 
# # Save the overall pooled raster
# terra::writeRaster(r_pooled, file.path(out_dir, "Overall_Pooled_Environmental_Suitability_AUC.tif"), overwrite = TRUE)




#alternative 3.0 - provide the transects

# ==============================================================================
# SPATIAL TRANSECT COUPLING PIPELINE
# ==============================================================================


# 1. Define Coordinates for Point A and Point B (WGS84)
# Adjust these coordinates to slice through your specific data coordinates
lon_A <- -1.56; lat_A <- 42.87  # Point A (Start, Pamplona)
lon_B <- -2.44;  lat_B <- 36.83  # Point B (End, Almeria)

# 2. Generate Ordered Sampling Points Along the Trajectory
n_samples <- 100
lon_seq <- seq(lon_A, lon_B, length.out = n_samples)
lat_seq <- seq(lat_A, lat_B, length.out = n_samples)

sampling_df <- data.frame(x = lon_seq, y = lat_seq)
sampling_vect <- vect(sampling_df, geom = c("x", "y"), crs = "+proj=longlat +datum=WGS84")

# Create a clean spatial line and endpoint markers for the map overlay
transect_line_sf <- st_sfc(st_linestring(matrix(c(lon_A, lat_A, lon_B, lat_B), ncol = 2, byrow = TRUE)), crs = 4326)
transect_line_vect <- vect(transect_line_sf)

endpoints_df <- data.frame(Label = c("A", "B"), x = c(lon_A, lon_B), y = c(lat_A, lat_B))
endpoints_vect <- vect(endpoints_df, geom = c("x", "y"), crs = "+proj=longlat +datum=WGS84")

# ==============================================================================
# 3. Extract Suitability Profiles Across Layers
# ==============================================================================
ext_clim   <- terra::extract(r_climate, sampling_vect)
ext_hydro  <- terra::extract(r_hydrocl, sampling_vect)
ext_morph  <- terra::extract(r_hydromo, sampling_vect)

# ==============================================================================
# 4. Compute Geodesic Distances along the Ellipsoid
# ==============================================================================
distances <- numeric(n_samples)
p_start <- vect(matrix(c(lon_A, lat_A), ncol = 2), crs = "+proj=longlat +datum=WGS84")

for (i in 1:n_samples) {
  p_current <- vect(matrix(c(lon_seq[i], lat_seq[i]), ncol = 2), crs = "+proj=longlat +datum=WGS84")
  distances[i] <- terra::distance(p_start, p_current, unit = "km")
}

# Consolidate and drop any edge pixels falling entirely outside the mask
profile_data <- data.frame(
  Distance_km = distances,
  Climate = ext_clim[, 2],
  Hydroclimatic = ext_hydro[, 2],
  Hydromorphological = ext_morph[, 2]
) %>% filter(!is.na(Climate) | !is.na(Hydroclimatic) | !is.na(Hydromorphological))

# Reshape for multi-line ggplot tracking
profile_long <- profile_data %>%
  pivot_longer(cols = c(Climate, Hydroclimatic, Hydromorphological),
               names_to = "Driver", values_to = "Suitability")

# Set factor levels to control the drawing order and legend order
profile_long$Driver <- factor(profile_long$Driver, 
                              levels = c("Climate", "Hydroclimatic", "Hydromorphological"))

# ==============================================================================
# 5. GRAPHICS GENERATION & COMPOSITING
# ==============================================================================
yellow_to_very_dark_orange <- colorRampPalette(c("yellow", "#993300"))(100)

# Top Panel: Pooled Suitability Base Map with Black Transparent Transect Trace
p_map_overlay <- ggplot() +
  geom_spatraster(data = r_pooled) +
  geom_spatvector(data = iberia_mainland, fill = NA, color = "black", linewidth = 0.5) +
  geom_spatvector(data = transect_line_vect, color = "black", alpha = 0.8, linewidth = 0.7) +
  geom_spatvector(data = endpoints_vect, color = "black", size = 3.5, shape = 21, fill = "white") +
  geom_spatvector_text(data = endpoints_vect, aes(label = Label), 
                       fontface = "bold", vjust = -1.2, size = 3,
                       nudge_x = -0.5, nudge_y= +0.2) +  # Adjust this value to move left (negative = left)
  scale_fill_gradientn(
    colors = yellow_to_very_dark_orange,
    name = "Pooled Suitability Index",
    limits = c(0, 100),
    na.value = "transparent"
  ) +
  labs(title = "Spatial Transect Trajectory (A to B)") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 10),
    axis.title = element_blank(), axis.text = element_blank(), panel.grid = element_blank()
  )

# Bottom Panel: Line Profile showing high spatial autocorrelation vs local variance
p_profile_lines <- ggplot(profile_long, aes(x = Distance_km, y = Suitability, color = Driver)) +
  geom_line(linewidth = 1.0, alpha = 0.9, lineend = "round") +
  scale_color_manual(
    values = c("Climate" = "#E69F00",             
               "Hydroclimatic" = "#56B4E9",       
               "Hydromorphological" = "#009E73") 
  ) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20), labels = paste0(seq(0, 100, 20), "%")) +
  labs(
    title = "Environmental Suitability Profiles along Transect Path",
    x = "Distance from Point A to B (km)",
    y = "Stacked Suitability Index",
    color = "Predictor Set"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    axis.title = element_text(face = "bold", size = 11),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

# Stack vertically using Patchwork
transect_composite_figure <- p_map_overlay / p_profile_lines + plot_layout(heights = c(1.3, 1))

print(transect_composite_figure)




### Figure 4, Combine both approaches into in grand ensemble plot


# ==============================================================================
# 0. SETUP AND LOAD PACKAGES
# # ==============================================================================
# library(dplyr)
# library(readr)
# library(terra)
# library(ggplot2)
# library(tidyterra)
# library(patchwork)
# library(ggcorrplot)
# library(rnaturalearth)
# library(grid)
# library(sf)
# library(tidyr)
# 
# setwd("/Users/geo_v/Desktop/")
# 
# # Directories
# ens_dir <- "Vagenas_aSDMs/output/figures/Regional/ensembles/regional_ensembles_H5_AUC"
# base_reg_dir <- "Vagenas_aSDMs/output/regional"
# out_dir <- "Vagenas_aSDMs/output/figures/Regional/ensembles/Figure4_stacked_unified_AUC"
# 
# if (!dir.exists(out_dir)) {
#   dir.create(out_dir, recursive = TRUE)
# }
# 
# # ==============================================================================
# # 1. PREPARE SPECIES CLASSIFICATION AND LISTS
# # ==============================================================================
# species_class <- read_csv("Vagenas_aSDMs/input/aSDMs_species_classification.csv", show_col_types = FALSE) %>%
#   mutate(Species = gsub(" ", "_", Sp))
# 
# # Get the species that actually have valid outputs
# valid_files <- list.files(ens_dir, pattern = "_H5_AUC_ensemble\\.tif$", full.names = TRUE)
# valid_species <- gsub("_H5_AUC_ensemble\\.tif$", "", basename(valid_files))
# 
# master_df <- species_class %>% filter(Species %in% valid_species)
# species_list <- master_df %>% pull(Species)
# 
# # ==============================================================================
# # 1.5 FETCH, MERGE, AND CLEAN IBERIAN MAINLAND BOUNDARY
# # ==============================================================================
# iberia_sf <- ne_countries(country = c("Spain", "Portugal"), scale = "medium", returnclass = "sf")
# iberia_vect <- vect(iberia_sf)
# 
# iberia_dissolved <- aggregate(iberia_vect)
# iberia_parts <- disagg(iberia_dissolved)
# land_areas <- expanse(iberia_parts, unit = "km")
# iberia_mainland <- iberia_parts[which.max(land_areas), ]
# 
# # ==============================================================================
# # 2. RASTER GENERATION FUNCTIONS (AUC FILTERINGED AND STANDARDIZED)
# # ==============================================================================
# 
# # FUNCTION A: Standard Sum (For Species Categories and Grand Pooled)
# safe_sum_rasters <- function(file_paths) {
#   r_list <- lapply(file_paths, rast)
#   master_ext <- ext(r_list[[1]])
#   if (length(r_list) > 1) {
#     for (i in 2:length(r_list)) { master_ext <- terra::union(master_ext, ext(r_list[[i]])) }
#   }
#   r_list_aligned <- lapply(r_list, function(r) extend(r, master_ext))
#   r_sum <- sum(rast(r_list_aligned), na.rm = TRUE)
#   
#   iberia_proj <- project(iberia_mainland, crs(r_sum))
#   r_masked <- mask(crop(r_sum, iberia_proj), iberia_proj)
#   
#   rmm <- minmax(r_masked)
#   r_norm <- if (rmm[2,1] > rmm[1,1]) ((r_masked - rmm[1,1]) / (rmm[2,1] - rmm[1,1])) * 100 else r_masked * 0
#   return(r_norm)
# }
# 
# # FUNCTION B: From-Scratch Driver Sum (Extracting RAW values weighted strictly by AUC)
# build_driver_raw_stack <- function(spp_list, target_driver) {
#   stack_list <- list()
#   for (sp in spp_list) {
#     file_path <- file.path(base_reg_dir, sp, "H5", target_driver, "ensemble_regional.tif")
#     eval_path <- file.path(base_reg_dir, sp, "H5", target_driver, "eval_regional_iberia_strict.csv")
#     
#     if (file.exists(file_path) && file.exists(eval_path)) {
#       metrics <- read_csv(eval_path, show_col_types = FALSE)
#       weight <- metrics$e_AUC[1] # Weighted strictly by AUC
#       
#       r <- rast(file_path)
#       rmm <- minmax(r)
#       r_norm <- if (rmm[2,1] > rmm[1,1]) (r - rmm[1,1]) / (rmm[2,1] - rmm[1,1]) else r * 0 
#       
#       stack_list[[sp]] <- r_norm * weight
#     }
#   }
#   
#   if(length(stack_list) == 0) return(NULL)
#   
#   master_ext <- ext(stack_list[[1]])
#   if (length(stack_list) > 1) {
#     for (i in 2:length(stack_list)) { master_ext <- terra::union(master_ext, ext(stack_list[[i]])) }
#   }
#   r_list_aligned <- lapply(stack_list, function(r) extend(r, master_ext))
#   r_sum <- sum(rast(r_list_aligned), na.rm = TRUE)
#   
#   iberia_proj <- project(iberia_mainland, crs(r_sum))
#   r_masked <- mask(crop(r_sum, iberia_proj), iberia_proj)
#   return(r_masked)
# }
# 
# # --- GENERATE SPECIES CATEGORY STACKS ---
# print("Stacking Species Category Ensembles...")
# spp_endemic  <- master_df %>% filter(Category == "Iberian Endemic") %>% pull(Species)
# r_endemic    <- safe_sum_rasters(file.path(ens_dir, paste0(spp_endemic, "_H5_AUC_ensemble.tif")))
# 
# spp_native   <- master_df %>% filter(Category == "Native Widespread") %>% pull(Species)
# r_native     <- safe_sum_rasters(file.path(ens_dir, paste0(spp_native, "_H5_AUC_ensemble.tif")))
# 
# spp_invasive <- master_df %>% filter(Category == "Invasive Widespread") %>% pull(Species)
# r_invasive   <- safe_sum_rasters(file.path(ens_dir, paste0(spp_invasive, "_H5_AUC_ensemble.tif")))
# 
# print("Stacking Master Pooled Ensemble...")
# r_pooled     <- safe_sum_rasters(valid_files)
# 
# # --- GENERATE DRIVER STACKS & APPLY PENINSULA-WIDE GLOBAL STANDARDIZATION ---
# print("Building Raw Driver Stacks...")
# r_climate_raw <- build_driver_raw_stack(species_list, "Climate")
# r_hydrocl_raw <- build_driver_raw_stack(species_list, "Hydroclimatic")
# r_hydromo_raw <- build_driver_raw_stack(species_list, "Hydromorphological")
# 
# print("Applying Cross-Predictor Global 0-100% Standardization...")
# global_stack  <- c(r_climate_raw, r_hydrocl_raw, r_hydromo_raw)
# global_minmax <- minmax(global_stack)
# global_min    <- min(global_minmax[1, ]) 
# global_max    <- max(global_minmax[2, ]) 
# 
# if (global_max > global_min) {
#   r_climate <- ((r_climate_raw - global_min) / (global_max - global_min)) * 100
#   r_hydrocl <- ((r_hydrocl_raw - global_min) / (global_max - global_min)) * 100
#   r_hydromo <- ((r_hydromo_raw - global_min) / (global_max - global_min)) * 100
# } else {
#   r_climate <- r_climate_raw * 0; r_hydrocl <- r_hydrocl_raw * 0; r_hydromo <- r_hydromo_raw * 0
# }
# 
# names(r_climate)  <- "Climate"
# names(r_hydrocl)  <- "Hydroclimatic"
# names(r_hydromo)  <- "Hydromorphological"
# names(r_pooled)   <- "Overall"
# 
# # Export TIF layers for safety
# writeRaster(r_pooled, file.path(out_dir, "Pooled_Overall_AUC_stacked.tif"), overwrite = TRUE)

# # ==============================================================================
# # 4. PLOT THE MAPS
# # ==============================================================================
# # Create map objects using the variables you already generated
# p_endemic  <- plot_map(r_endemic,  paste0("Iberian Endemic (N=", length(spp_endemic), ")"), show_legend = FALSE)
# p_native   <- plot_map(r_native,   paste0("Native Widespread (N=", length(spp_native), ")"), show_legend = FALSE)
# p_invasive <- plot_map(r_invasive, paste0("Invasive Widespread (N=", length(spp_invasive), ")"), show_legend = FALSE)
# 
# # Pooled Hub map
# # ==============================================================================
# # CORRECTED: POOLED HUB MAP WITH TRANSECT LINE
# # ==============================================================================
# p_pooled_hub <- plot_map(r_pooled, "All species (N=98)", show_legend = TRUE) + 
#   theme(panel.border = element_rect(color = "black", fill = scales::alpha("black", 0.05), linewidth = 1.5)) +
#   coord_sf(expand = FALSE, clip = "off") + 
#   
#   # 1. ADD THE TRANSECT LINE (This makes it appear on the map)
#   geom_spatvector(data = transect_line_vect, color = "black", alpha = 0.8, linewidth = 0.7) +
#   geom_spatvector(data = endpoints_vect, color = "black", size = 2.5, shape = 21, fill = "white") +
#   geom_spatvector_text(data = endpoints_vect, aes(label = Label), fontface = "bold", vjust = -1.5, size = 5) +
#   
#   # 2. Add the spokes/connector lines
#   annotation_custom(segmentsGrob(x0 = 0.5, y0 = 1.0, x1 = 0.5,  y1 = 1.25, gp = gpar(lwd = 2.0))) + 
#   annotation_custom(segmentsGrob(x0 = 0.5, y0 = 1.0, x1 = -1.2, y1 = 1.25, gp = gpar(lwd = 2.0))) + 
#   annotation_custom(segmentsGrob(x0 = 0.5, y0 = 1.0, x1 = 2.2,  y1 = 1.25, gp = gpar(lwd = 2.0))) + 
#   annotation_custom(pointsGrob(x = 0.5, y = 1.0, pch = 21, gp = gpar(fill = "black", col = "black", cex = 1.2)))
# 
# # ==============================================================================
# # 5. ASSEMBLE FINAL LAYOUT
# # ==============================================================================
# 
# # Define the layout design matrix
# layout_design <- c(
#   patchwork::area(t = 1, l = 2, b = 1, r = 2), # Endemic (Top)
#   patchwork::area(t = 2, l = 1, b = 2, r = 1), # Native (Mid Left)
#   patchwork::area(t = 2, l = 3, b = 2, r = 3), # Invasive (Mid Right)
#   patchwork::area(t = 2, l = 2, b = 2, r = 2)  # Pooled (Center)
# )
# 
# # 1. Top block (The 4 maps)
# top_block <- (p_endemic + p_native + p_invasive + p_pooled_hub) + 
#   plot_layout(design = layout_design, guides = "collect") & 
#   theme(legend.position = "bottom", legend.justification = "bottom")
# 
# # 2. Final assembly (Top block / Profile lines)
# # Ensure p_profile_lines is defined from your transect script
# final_figure <- top_block / p_profile_lines + 
#   plot_layout(heights = c(3.5, 1)) & 
#   theme(plot.background = element_rect(fill = "white", color = NA))
# 
# # 3. Save
# ggsave(file.path(out_dir, "Final_Unified_Cross_Transect_Dashboard.png"), 
#        final_figure, width = 22, height = 22, dpi = 300, bg = "white")
# 
# 
# print("Complete! Unified dashboard with top-aligned legend saved.")


# ==============================================================================
# 4. PLOT THE MAPS (WITH ENLARGED FONTS)
# ==============================================================================
# library(grid)
# 
# # 1) INCREASE FONT SIZE: Apply a global theme override for map titles
# title_theme <- theme(plot.title = element_text(size = 26, face = "bold", hjust = 0.5))
# 
# p_endemic  <- plot_map(r_endemic,  paste0("Iberian Endemic (N=", length(spp_endemic), ")"), show_legend = FALSE) + title_theme
# p_native   <- plot_map(r_native,   paste0("Native Widespread (N=", length(spp_native), ")"), show_legend = FALSE) + title_theme
# p_invasive <- plot_map(r_invasive, paste0("Invasive Widespread (N=", length(spp_invasive), ")"), show_legend = FALSE) + title_theme
# 
# 
# # ==============================================================================
# # CORRECTED: POOLED HUB MAP WITH SINGLE CONVERGENCE DOT
# # ==============================================================================
# p_pooled_hub <- plot_map(r_pooled, "All species (N=98)", show_legend = TRUE) + 
#   title_theme +
#   theme(
#     panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5),
#     # Maintain the margin so the lines have room to draw
#     plot.margin = margin(t = 60, r = 60, b = 60, l = 60, unit = "pt")
#   ) +
#   coord_sf(expand = FALSE, clip = "off") + 
#   
#   # A-B TRANSECT LINE
#   geom_spatvector(data = transect_line_vect, color = "black", alpha = 0.8, linewidth = 0.9) +
#   geom_spatvector(data = endpoints_vect, color = "black", size = 3, shape = 21, fill = "white") +
#   geom_spatvector_text(data = endpoints_vect, aes(label = Label), fontface = "bold", vjust = -1.5, size = 6) +
#   
#   # ======================================================
# # NEW CONNECTOR SPOKES (Converging at ONE dot)
# # ======================================================
# 
# # 1. The SINGLE central convergence dot at the top
# annotation_custom(pointsGrob(x = 0.5, y = 1.0, pch = 21, gp = gpar(fill = "black", col = "black", cex = 1.5))) +
#   
#   # 2. Line straight UP to Endemic
#   annotation_custom(segmentsGrob(x0 = 0.5, y0 = 1.0, x1 = 0.5, y1 = 1.25, gp = gpar(lwd = 3.0))) + 
#   
#   # 3. Line diagonally UP-LEFT to Native
#   annotation_custom(segmentsGrob(x0 = 0.5, y0 = 1.0, x1 = -0.25, y1 = 1.15, gp = gpar(lwd = 3.0))) + 
#   
#   # 4. Line diagonally UP-RIGHT to Invasive
#   annotation_custom(segmentsGrob(x0 = 0.5, y0 = 1.0, x1 = 1.25, y1 = 1.15, gp = gpar(lwd = 3.0)))
# 
# # ==============================================================================
# # 5. ASSEMBLE FINAL LAYOUT
# # ==============================================================================
# 
# # 2) CENTER AND ENLARGE BOTTOM PLOT TITLE
# p_profile_lines <- p_profile_lines + 
#   theme(
#     plot.title = element_text(hjust = 0.5, size = 26, face = "bold", margin = margin(b = 20)),
#     axis.title = element_text(size = 18, face = "bold"),
#     axis.text = element_text(size = 14),
#     legend.title = element_text(size = 20, face = "bold"),
#     legend.text = element_text(size = 16)
#   )
# 
# # Define the layout design matrix
# layout_design <- c(
#   patchwork::area(t = 1, l = 2, b = 1, r = 2), # Endemic (Top)
#   patchwork::area(t = 2, l = 1, b = 2, r = 1), # Native (Mid Left)
#   patchwork::area(t = 2, l = 3, b = 2, r = 3), # Invasive (Mid Right)
#   patchwork::area(t = 2, l = 2, b = 2, r = 2)  # Pooled (Center)
# )
# 
# # 3) MAKE CENTER MAP BIGGER: Use widths/heights inside plot_layout
# top_block <- (p_endemic + p_native + p_invasive + p_pooled_hub) + 
#   plot_layout(
#     design = layout_design, 
#     widths = c(1, 2.5, 1), # Central column is 2.5x wider
#     heights = c(1, 2.5),   # Second row is 2.5x taller
#     guides = "collect"
#   ) & 
#   theme(
#     legend.position = "bottom", 
#     legend.justification = "center",
#     legend.title = element_text(size = 20, face = "bold"),
#     legend.text = element_text(size = 16)
#   )
# 
# # Final assembly (Top block / Profile lines)
# final_figure <- top_block / p_profile_lines + 
#   plot_layout(heights = c(4, 1.2)) & 
#   theme(plot.background = element_rect(fill = "white", color = NA))
# 
# # Save
# ggsave(file.path(out_dir, "Final_Unified_Cross_Transect_Dashboard.png"), 
#        final_figure, width = 24, height = 24, dpi = 300, bg = "white")
# 
# print("Complete! Unified dashboard with visible lines, larger center map, and bolded/centered fonts saved.")
# 
# 





#test that it works
# 
# # ==============================================================================
# # 0. SETUP AND LOAD PACKAGES
# # ==============================================================================
# library(dplyr)
# library(readr)
# library(terra)
# library(ggplot2)
# library(tidyterra)
# library(patchwork)
# library(sf)
# library(tidyr)
# library(rnaturalearth)
# library(cowplot)
# library(grid)
# 
# setwd("/Users/geo_v/Desktop/")
# ens_dir <- "Vagenas_aSDMs/output/regional_ensembles_H5_AUC"
# base_reg_dir <- "Vagenas_aSDMs/output/regional"
# out_dir <- "Vagenas_aSDMs/output/stacked_maps_unified_AUC"
# 
# # ==============================================================================
# # 3. PLOT MAPS (BALANCED MARGINS TO FIX SIZING)
# # ==============================================================================
# 
# # 1) RESTORE ENDEMIC SIZE
# p_endemic  <- plot_map(r_endemic,  "Iberian Endemic (N=47)", show_legend = FALSE)
# 
# # 2) PULL SIDE MAPS UP (Without creating a huge invisible gap below them)
# p_native   <- plot_map(r_native,   "Native Widespread (N=26)", show_legend = FALSE) + 
#   theme(plot.margin = margin(t = -120, r = -20, b = 20, l = 0, unit = "pt"))
# 
# p_invasive <- plot_map(r_invasive, "Invasive Widespread (N=25)", show_legend = FALSE) + 
#   theme(plot.margin = margin(t = -120, r = 0, b = 20, l = -20, unit = "pt"))
# 
# # 3) MOVE CENTRAL HUB DOWN
# p_pooled_hub <- plot_map(r_pooled, title = NULL, show_legend = TRUE) + 
#   theme(
#     panel.border = element_rect(color = "black", fill = NA, linewidth = 1.5),
#     # The t = 80 pushes this specific map DOWN inside its grid cell
#     plot.margin = margin(t = 80, r = 20, b = 0, l = 20, unit = "pt")
#   ) +
#   coord_sf(expand = FALSE) + 
#   geom_spatvector(data = transect_line_vect, color = "black", alpha = 0.8, linewidth = 0.9) +
#   geom_spatvector(data = endpoints_vect, color = "black", size = 3, shape = 21, fill = "white") +
#   geom_spatvector_text(data = endpoints_vect, aes(label = Label), fontface = "bold", vjust = -1.5, size = 6) +
#   annotation_custom(textGrob("All species (N=98)", x = 0.5, y = 0.96, 
#                              gp = gpar(fontsize = 28, fontface = "bold", col = "black")))
# 
# 
# # ==============================================================================
# # 4. ASSEMBLE BASE PATCHWORK GRID (FIXING ENDEMIC SIZE & SPACING)
# # ==============================================================================
# 
# p_profile_lines <- p_profile_lines + 
#   theme(
#     # FIX 2: Changed t = -40 to t = 10 to add breathing room below the map block
#     plot.margin = margin(t = 5, r = 20, b = 10, l = 20, unit = "pt"),
#     plot.title = element_text(hjust = 0.5, size = 32, face = "bold", margin = margin(b = 20)),
#     axis.title = element_text(size = 24, face = "bold"), 
#     axis.text = element_text(size = 18),                 
#     legend.title = element_text(size = 26, face = "bold"), 
#     legend.text = element_text(size = 22),                 
#     legend.key.width = unit(5, "cm") 
#   )
# 
# layout_design <- "
# #A#
# BCD
# "
# 
# # FIX 1: Re-balance the grid so Endemic isn't squished
# top_block <- (p_endemic + p_native + p_pooled_hub + p_invasive) + 
#   plot_layout(
#     design = layout_design, 
#     widths = c(1.4, 2.5, 1.4), 
#     heights = c(1.7, 3), # Increased the top row height from 1.2 to 1.8 so Endemic grows larger!
#     guides = "collect"
#   ) & 
#   theme(
#     legend.position = "bottom", legend.justification = "center",
#     legend.title = element_text(size = 22, face = "bold"),
#     legend.text = element_text(size = 18),
#     legend.margin = margin(t = -10, b = 20, unit = "pt") # Added b = 20 to pad the bottom of the color band
#   )
# 
# base_figure <- top_block / p_profile_lines + 
#   plot_layout(heights = c(3.2, 1)) & 
#   theme(plot.background = element_rect(fill = "white", color = NA))
# 
# 
# # ==============================================================================
# # 5. DRAW FINAL CALIBRATED LINES WITH COWPLOT
# # ==============================================================================
# # Because we changed the grid heights to make the Endemic map bigger, the central 
# # map shifts down by a microscopic amount. I've adjusted y down to 0.68 to compensate.
# 
# final_figure <- ggdraw(base_figure) +
#   
#   # Central dot
#   draw_grob(pointsGrob(x = 0.5, y = 0.71, pch = 21, gp = gpar(fill = "black", col = "black", cex = 1.5))) +
#   
#   # Line straight UP to Endemic
#   draw_grob(segmentsGrob(x0 = 0.5, y0 = 0.71, x1 = 0.5, y1 = 0.8, gp = gpar(lwd = 3))) +
#   
#   # Line diagonally UP-LEFT to Native 
#   draw_grob(segmentsGrob(x0 = 0.5, y0 = 0.71, x1 = 0.3, y1 = 0.72, gp = gpar(lwd = 3))) +
#   
#   # Line diagonally UP-RIGHT to Invasive 
#   draw_grob(segmentsGrob(x0 = 0.5, y0 = 0.71, x1 = 0.73, y1 = 0.72, gp = gpar(lwd = 3)))
# 
# # Exporting 
# suppressWarnings({
#   ggsave(file.path(out_dir, "Final_Aligned_Dashboard_Perfected.png"), 
#          final_figure, 
#          width = 26,   
#          height = 18,  
#          dpi = 300, bg = "white")
# })
# 
# print("Success! Endemic map enlarged, spacing added below legend, and lines calibrated.")


# 
# #another final test
# 
# # ==============================================================================
# # 0. SETUP AND LOAD PACKAGES
# # ==============================================================================
# library(dplyr)
# library(readr)
# library(terra)
# library(ggplot2)
# library(tidyterra)
# library(patchwork)
# library(sf)
# library(tidyr)
# library(rnaturalearth)
# library(cowplot)
# library(grid)
# 
# setwd("/Users/geo_v/Desktop/")
# ens_dir <- "Vagenas_aSDMs/output/regional_ensembles_H5_AUC"
# base_reg_dir <- "Vagenas_aSDMs/output/regional"
# out_dir <- "Vagenas_aSDMs/output/stacked_maps_unified_AUC"
# 
# if (!dir.exists(out_dir)) {
#   dir.create(out_dir, recursive = TRUE)
# }
# 
# # ==============================================================================
# # 1. PREPARE SPECIES CLASSIFICATION AND LISTS
# # ==============================================================================
# species_class <- read_csv("Vagenas_aSDMs/input/aSDMs_species_classification.csv", show_col_types = FALSE) %>%
#   mutate(Species = gsub(" ", "_", Sp))
# 
# # Get the species that actually have valid outputs
# valid_files <- list.files(ens_dir, pattern = "_H5_AUC_ensemble\\.tif$", full.names = TRUE)
# valid_species <- gsub("_H5_AUC_ensemble\\.tif$", "", basename(valid_files))
# 
# master_df <- species_class %>% filter(Species %in% valid_species)
# species_list <- master_df %>% pull(Species)
# 
# # ==============================================================================
# # 1.5 FETCH, MERGE, AND CLEAN IBERIAN MAINLAND BOUNDARY
# # ==============================================================================
# iberia_sf <- ne_countries(country = c("Spain", "Portugal"), scale = "medium", returnclass = "sf")
# iberia_vect <- vect(iberia_sf)
# 
# iberia_dissolved <- aggregate(iberia_vect)
# iberia_parts <- disagg(iberia_dissolved)
# land_areas <- expanse(iberia_parts, unit = "km")
# iberia_mainland <- iberia_parts[which.max(land_areas), ]
# 
# # ==============================================================================
# # 2. RASTER GENERATION FUNCTIONS (AUC FILTERING AND STANDARDIZED)
# # ==============================================================================
# 
# # FUNCTION A: Standard Sum (For Species Categories and Grand Pooled)
# safe_sum_rasters <- function(file_paths) {
#   r_list <- lapply(file_paths, rast)
#   master_ext <- ext(r_list[[1]])
#   if (length(r_list) > 1) {
#     for (i in 2:length(r_list)) { master_ext <- terra::union(master_ext, ext(r_list[[i]])) }
#   }
#   r_list_aligned <- lapply(r_list, function(r) extend(r, master_ext))
#   r_sum <- sum(rast(r_list_aligned), na.rm = TRUE)
#   
#   iberia_proj <- project(iberia_mainland, crs(r_sum))
#   r_masked <- mask(crop(r_sum, iberia_proj), iberia_proj)
#   
#   rmm <- minmax(r_masked)
#   r_norm <- if (rmm[2,1] > rmm[1,1]) ((r_masked - rmm[1,1]) / (rmm[2,1] - rmm[1,1])) * 100 else r_masked * 0
#   return(r_norm)
# }
# 
# # --- GENERATE SPECIES CATEGORY STACKS ---
# print("Stacking Species Category Ensembles...")
# spp_endemic  <- master_df %>% filter(Category == "Iberian Endemic") %>% pull(Species)
# r_endemic    <- safe_sum_rasters(file.path(ens_dir, paste0(spp_endemic, "_H5_AUC_ensemble.tif")))
# 
# spp_native   <- master_df %>% filter(Category == "Native Widespread") %>% pull(Species)
# r_native     <- safe_sum_rasters(file.path(ens_dir, paste0(spp_native, "_H5_AUC_ensemble.tif")))
# 
# spp_invasive <- master_df %>% filter(Category == "Invasive Widespread") %>% pull(Species)
# r_invasive   <- safe_sum_rasters(file.path(ens_dir, paste0(spp_invasive, "_H5_AUC_ensemble.tif")))
# 
# print("Stacking Master Pooled Ensemble...")
# r_pooled     <- safe_sum_rasters(valid_files)
# 
# # ==============================================================================
# # 2.5 TRANSECT PATH PREP & PROFILE LINES GENERATION
# # ==============================================================================
# lon_A <- -8.93; lat_A <- 38.73  
# lon_B <- 2.14;  lat_B <- 41.45  
# transect_line_sf <- st_sfc(st_linestring(matrix(c(lon_A, lat_A, lon_B, lat_B), ncol = 2, byrow = TRUE)), crs = 4326)
# transect_line_vect <- vect(transect_line_sf)
# 
# endpoints_df <- data.frame(Label = c("A", "B"), x = c(lon_A, lon_B), y = c(lat_A, lat_B))
# endpoints_vect <- vect(endpoints_df, geom = c("x", "y"), crs = "+proj=longlat +datum=WGS84")
# 
# 





#STANDALONE - CORRECT FIGURE 4

# ==============================================================================
# 0. SETUP AND LOAD PACKAGES
# ==============================================================================
library(dplyr)
library(readr)
library(terra)
library(ggplot2)
library(tidyterra)
library(patchwork)
library(sf)
library(tidyr)
library(rnaturalearth)
library(cowplot)
library(grid)

yellow_to_very_dark_orange <- colorRampPalette(c("yellow", "#993300"))(100)

setwd("/Users/geo_v/Desktop/")
ens_dir <- "Vagenas_aSDMs/output/figures/Regional/ensembles/regional_ensembles_H5_AUC"
base_reg_dir <- "Vagenas_aSDMs/output/regional"
out_dir <- "Vagenas_aSDMs/output/figures/Regional/ensembles/Figure4_stacked_maps_unified_AUC"

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ==============================================================================
# 1. PREPARE DATA, SPECIES LISTS & SHAPEFILES
# ==============================================================================
species_class <- read_csv("Vagenas_aSDMs/input/aSDMs_species_classification.csv", show_col_types = FALSE) %>%
  mutate(Species = gsub(" ", "_", Sp))

valid_files <- list.files(ens_dir, pattern = "_H5_AUC_ensemble\\.tif$", full.names = TRUE)
valid_species <- gsub("_H5_AUC_ensemble\\.tif$", "", basename(valid_files))
master_df <- species_class %>% filter(Species %in% valid_species)
species_list <- master_df %>% pull(Species)

# Fetch Iberian Mainland Boundary
iberia_sf <- ne_countries(country = c("Spain", "Portugal"), scale = "medium", returnclass = "sf")
iberia_vect <- vect(iberia_sf)
iberia_dissolved <- aggregate(iberia_vect)
iberia_parts <- disagg(iberia_dissolved)
iberia_mainland <- iberia_parts[which.max(expanse(iberia_parts, unit = "km")), ]

# ==============================================================================
# 2. RASTER STACKING & PREDICTOR DRIVER EXTRACTION
# ==============================================================================
# FUNCTION A: Standard Sum (For Species Categories)
safe_sum_rasters <- function(file_paths) {
  r_list <- lapply(file_paths, rast)
  master_ext <- ext(r_list[[1]])
  if (length(r_list) > 1) for (i in 2:length(r_list)) master_ext <- terra::union(master_ext, ext(r_list[[i]]))
  
  r_list_aligned <- lapply(r_list, function(r) extend(r, master_ext))
  r_sum <- sum(rast(r_list_aligned), na.rm = TRUE)
  iberia_proj <- project(iberia_mainland, crs(r_sum))
  r_masked <- mask(crop(r_sum, iberia_proj), iberia_proj)
  
  rmm <- minmax(r_masked)
  r_norm <- if (rmm[2,1] > rmm[1,1]) ((r_masked - rmm[1,1]) / (rmm[2,1] - rmm[1,1])) * 100 else r_masked * 0
  return(r_norm)
}

# FUNCTION B: Raw Driver Extraction (For the Distance Profile)
build_driver_raw_stack <- function(spp_list, target_driver) {
  stack_list <- list()
  for (sp in spp_list) {
    file_path <- file.path(base_reg_dir, sp, "H5", target_driver, "ensemble_regional.tif")
    eval_path <- file.path(base_reg_dir, sp, "H5", target_driver, "eval_regional_iberia_strict.csv")
    
    if (file.exists(file_path) && file.exists(eval_path)) {
      metrics <- read_csv(eval_path, show_col_types = FALSE)
      weight <- metrics$e_AUC[1]
      
      r <- rast(file_path)
      rmm <- minmax(r)
      r_norm <- if (rmm[2,1] > rmm[1,1]) (r - rmm[1,1]) / (rmm[2,1] - rmm[1,1]) else r * 0 
      stack_list[[sp]] <- r_norm * weight
    }
  }
  if(length(stack_list) == 0) return(NULL)
  
  master_ext <- ext(stack_list[[1]])
  if (length(stack_list) > 1) { for (i in 2:length(stack_list)) master_ext <- terra::union(master_ext, ext(stack_list[[i]])) }
  r_list_aligned <- lapply(stack_list, function(r) extend(r, master_ext))
  r_sum <- sum(rast(r_list_aligned), na.rm = TRUE)
  
  iberia_proj <- project(iberia_mainland, crs(r_sum))
  return(mask(crop(r_sum, iberia_proj), iberia_proj))
}

print("Stacking Master Ensembles...")
spp_endemic  <- master_df %>% filter(Category == "Iberian Endemic") %>% pull(Species)
r_endemic    <- safe_sum_rasters(file.path(ens_dir, paste0(spp_endemic, "_H5_AUC_ensemble.tif")))
spp_native   <- master_df %>% filter(Category == "Native Widespread") %>% pull(Species)
r_native     <- safe_sum_rasters(file.path(ens_dir, paste0(spp_native, "_H5_AUC_ensemble.tif")))
spp_invasive <- master_df %>% filter(Category == "Invasive Widespread") %>% pull(Species)
r_invasive   <- safe_sum_rasters(file.path(ens_dir, paste0(spp_invasive, "_H5_AUC_ensemble.tif")))
r_pooled     <- safe_sum_rasters(valid_files)

print("Building Raw Driver Stacks for Transect Profile...")
r_climate_raw <- build_driver_raw_stack(species_list, "Climate")
r_hydrocl_raw <- build_driver_raw_stack(species_list, "Hydroclimatic")
r_hydromo_raw <- build_driver_raw_stack(species_list, "Hydromorphological")

print("Applying Cross-Predictor Global 0-100% Standardization...")
global_stack  <- c(r_climate_raw, r_hydrocl_raw, r_hydromo_raw)
global_minmax <- minmax(global_stack)
global_min    <- min(global_minmax[1, ]) 
global_max    <- max(global_minmax[2, ]) 

if (global_max > global_min) {
  r_climate <- ((r_climate_raw - global_min) / (global_max - global_min)) * 100
  r_hydrocl <- ((r_hydrocl_raw - global_min) / (global_max - global_min)) * 100
  r_hydromo <- ((r_hydromo_raw - global_min) / (global_max - global_min)) * 100
} else {
  r_climate <- r_climate_raw * 0; r_hydrocl <- r_hydrocl_raw * 0; r_hydromo <- r_hydromo_raw * 0
}

# ==============================================================================
# 3. SPATIAL TRANSECT EXTRACTION PIPELINE
# ==============================================================================
lon_A <- -1.56; lat_A <- 42.87  # Point A (Start, Pamplona)
lon_B <- -2.44;  lat_B <- 36.83  # Point B (End, Almeria)

n_samples <- 80
sampling_vect <- vect(data.frame(x = seq(lon_A, lon_B, length.out = n_samples), 
                                 y = seq(lat_A, lat_B, length.out = n_samples)), 
                      geom = c("x", "y"), crs = "+proj=longlat +datum=WGS84")

transect_line_vect <- vect(st_sfc(st_linestring(matrix(c(lon_A, lat_A, lon_B, lat_B), ncol = 2, byrow = TRUE)), crs = 4326))
endpoints_vect <- vect(data.frame(Label = c("A", "B"), x = c(lon_A, lon_B), y = c(lat_A, lat_B)), 
                       geom = c("x", "y"), crs = "+proj=longlat +datum=WGS84")

ext_clim  <- terra::extract(r_climate, sampling_vect)
ext_hydro <- terra::extract(r_hydrocl, sampling_vect)
ext_morph <- terra::extract(r_hydromo, sampling_vect)
ext_stacked <- terra::extract(r_pooled, sampling_vect)

distances <- numeric(n_samples)
p_start <- vect(matrix(c(lon_A, lat_A), ncol = 2), crs = "+proj=longlat +datum=WGS84")
for (i in 1:n_samples) distances[i] <- terra::distance(p_start, sampling_vect[i], unit = "km")

profile_long <- data.frame(
  Distance_km = distances, Climate = ext_clim[, 2], 
  Hydroclimatic = ext_hydro[, 2], Hydromorphological = ext_morph[, 2], Stacked = ext_stacked[,2]
) %>% 
  filter(!is.na(Climate) | !is.na(Hydroclimatic) | !is.na(Hydromorphological) | !is.na(Stacked)) %>%
  pivot_longer(cols = c(Climate, Hydroclimatic, Hydromorphological,Stacked), names_to = "Driver", values_to = "Suitability") %>%
  mutate(Driver = factor(Driver, levels = c("Climate", "Hydroclimatic", "Hydromorphological","Stacked")))

# # ==============================================================================
# # 3. PLOT MAPS (TRANSECT TWEAKS & TEXT SIZES)
# # ==============================================================================
# 
# plot_map <- function(r, title, show_legend = TRUE) {
#   p <- ggplot() +
#     geom_spatraster(data = r) +
#     geom_spatvector(data = iberia_mainland, fill = NA, color = "black", linewidth = 0.4) +
#     scale_fill_gradientn(
#       colors = yellow_to_very_dark_orange,
#       name = "Stacked Habitat Suitability",
#       limits = c(0, 100),
#       breaks = seq(0, 100, by = 20),
#       labels = paste0(seq(0, 100, by = 20), "%"),
#       na.value = "transparent",
#       guide = guide_colorbar(
#         direction = "horizontal",
#         label.position = "bottom",
#         title.position = "top",
#         title.hjust = 0.5,
#         barwidth = unit(25, "cm"),
#         barheight = unit(2, "cm"),
#         title.theme = element_text(size = 26, face = "bold"),  # Title text size
#         label.theme = element_text(size = 26, face = "bold")   # Break label text size
#       )
#     ) +
#     labs(title = title) +
#     theme_minimal() +
#     theme(
#       plot.title = element_text(hjust = 0.5, face = "bold", size = 32),
#       axis.text = element_blank(),
#       axis.title = element_blank(),
#       panel.grid = element_blank(),
#       plot.margin = margin(5, 5, 5, 5, unit = "pt"),
#       # Make legend text bold
#       legend.title = element_text(face = "bold",size=24),
#       legend.text = element_text(face = "bold",size=24)
#     )
#   
#   # Conditionally hide the legend
#   if (!show_legend) {
#     p <- p + theme(legend.position = "none")
#   }
#   
#   return(p)
# }
# 
# 
# 
# # 1) SIDE AND TOP MAPS (Margins neutralized since layout logic will handle size)
# p_endemic  <- plot_map(r_endemic,  "Iberian Endemic (N=47)", show_legend = FALSE) +
#   theme(plot.margin = margin(b = -20, unit = "pt")) # Pulls it slightly closer to the center
# 
# p_native   <- plot_map(r_native,   "Native Widespread (N=26)", show_legend = FALSE) + 
#   theme(plot.margin = margin(t = -80, r = -20, b = 20, l = 0, unit = "pt"))
# 
# p_invasive <- plot_map(r_invasive, "Invasive Widespread (N=25)", show_legend = FALSE) + 
#   theme(plot.margin = margin(t = -80, r = 0, b = 20, l = -20, unit = "pt"))
# 
# # 2) CENTRAL HUB (Thicker/transparent transect line & bigger A/B text)
# p_pooled_hub <- plot_map(r_pooled, title = NULL, show_legend = TRUE) + 
#   theme(
#     panel.border = element_rect(color = "black", fill = NA, linewidth = 2),
#     plot.margin = margin(t = 40, r = 40, b = 40, l = 40, unit = "pt")  # Increased all margins
#     ) +
#   coord_sf(expand = c(0.05, 0.05, 0.05, 0.05), clip = "off") + 
#   # Transect line: alpha 0.4 for transparency, linewidth 2.0 for thickness
#   geom_spatvector(data = transect_line_vect, color = "black", alpha = 0.5, linewidth = 2.5) +
#   geom_spatvector(data = endpoints_vect, color = "black", size = 3, shape = 21, fill = "black") +
#   # A and B Labels: size increased to 10
#   geom_spatvector_text(data = endpoints_vect, aes(label = Label), fontface = "bold", vjust = -1.2, size = 14,
#                        nudge_x = -0.5, nudge_y= -1.1) +
#   annotation_custom(textGrob("All species (N=98)", x = 0.5, y = 0.96, 
#                              gp = gpar(fontsize = 30, fontface = "bold", col = "black")))
# 
# 
# # ==============================================================================
# # 4. CREATE PROFILE PLOT, EXTRACT ITS LEGEND, AND BUILD CORRECT LAYOUT
# # ==============================================================================
# 
# # --- Build the transect profile plot (no legend yet) ---
# p_profile_lines <- ggplot(profile_long, aes(x = Distance_km, y = Suitability, color = Driver)) +
#   geom_line(linewidth = 1.5) +
#   scale_color_manual(values = c("Climate" = "#E69F00", 
#                                 "Hydroclimatic" = "#56B4E9", 
#                                 "Hydromorphological" = "#009E73","Stacked"="black")) +
#   scale_y_continuous(
#     limits = c(0, 100),           # Force y-axis from 0 to 100
#     breaks = seq(0, 100, by = 25), # Breaks every 25 units
#     labels = paste0(seq(0, 100, by = 25), "%")  # Add % to labels
#   ) +
#   labs(title = "Freshwater Stacked Habitat Suitability along transect A → B",
#        x = "Distance (km)", y = "Stacked Habitat Suitability") +
#   theme_minimal() +
#   theme(
#     plot.margin = margin(t = 5, r = 20, b = 10, l = 20, unit = "pt"),
#     plot.title = element_text(hjust = 0.5, size = 32, face = "bold", margin = margin(b = 20)),
#     axis.title = element_text(size = 27, face = "bold"),     # Already set to 27
#     axis.title.x = element_text(size = 27, face = "bold"),   # Explicitly set x-axis title
#     axis.title.y = element_text(size = 27, face = "bold"),   # Explicitly set y-axis title
#     axis.text = element_text(size = 23),                     # Already set to 23
#     axis.text.x = element_text(size = 23),                   # Explicitly set x-axis text
#     axis.text.y = element_text(size = 23),                   # Explicitly set y-axis text
#     legend.title = element_blank(),
#     legend.text = element_text(size = 27),                 
#     legend.key.width = unit(7, "cm"),
#     # Remove the legend from the plot itself — we will place it separately
#     legend.position = "none"
#   ) +
#   guides(color = guide_legend(override.aes = list(linewidth = 4)))
# 
# # Extract the legend from a copy that still has it
# p_profile_legend <- cowplot::get_legend(
#   p_profile_lines + theme(legend.position = "bottom")
# )
# 
# # Stack profile plot and its legend vertically
# profile_block <- wrap_plots(p_profile_lines, p_profile_legend, ncol = 1, heights = c(12, 1))
# 
# # --- Build the map panel (same as before) ---
# row1 <- plot_spacer() + p_endemic + plot_spacer() + 
#   plot_layout(widths = c(1.15, 1, 1.15))
# row2 <- p_native + p_pooled_hub + p_invasive + 
#   plot_layout(widths = c(1, 1.3, 1))
# top_block <- row1 / row2 + 
#   plot_layout(heights = c(1, 1.4), guides = "collect") & 
#   theme(
#     legend.position = "bottom", legend.justification = "center",
#     legend.text = element_text(size = 18),
#     legend.margin = margin(t = 15, b = 30, unit = "pt")  
#     )
# 
# # --- Combine everything with wrap_plots to avoid wrap_dims error ---
# base_figure <- wrap_plots(top_block, profile_block, ncol = 1, heights = c(3.2, 1.3)) &
#   theme(plot.background = element_rect(fill = "white", color = NA))
# 
# # ==============================================================================
# # 5. DRAW CONNECTING LINES WITH ggdraw
# # ==============================================================================
# final_figure <- ggdraw(base_figure) +
#   # Central dot
#   draw_grob(pointsGrob(x = 0.5, y = 0.705, pch = 21, gp = gpar(fill = "black", col = "black", cex = 3))) +
#   # Line straight UP to Endemic
#   draw_grob(segmentsGrob(x0 = 0.5, y0 = 0.705, x1 = 0.5, y1 = 0.77, gp = gpar(lwd = 3))) +
#   # Line diagonally UP-LEFT to Native 
#   draw_grob(segmentsGrob(x0 = 0.5, y0 = 0.705, x1 = 0.28, y1 = 0.705, gp = gpar(lwd = 3))) +
#   # Line diagonally UP-RIGHT to Invasive 
#   draw_grob(segmentsGrob(x0 = 0.5, y0 = 0.705, x1 = 0.74, y1 = 0.705, gp = gpar(lwd = 3)))
# 
# 
# # Exporting
# suppressWarnings({
#   ggsave(file.path(out_dir, "Final_Aligned_Dashboard_Perfected.png"),
#          final_figure,
#          width = 38,
#          height = 26,
#          dpi = 300, bg = "white")
# })


bbox <- sf::st_bbox(iberia_mainland)
common_xlim <- c(bbox["xmin"], bbox["xmax"])
common_ylim <- c(bbox["ymin"], bbox["ymax"])
yrange <- diff(common_ylim)

plot_map <- function(r, title, show_legend = TRUE, title_top_margin = 0,
                     xlim = NULL, ylim = NULL) {
  p <- ggplot() +
    geom_spatraster(data = r) +
    geom_spatvector(data = iberia_mainland, fill = NA, color = "black", linewidth = 0.4) +
    scale_fill_gradientn(
      colors = yellow_to_very_dark_orange,
      name = "Stacked Habitat Suitability",
      limits = c(0, 100),
      breaks = seq(0, 100, by = 20),
      labels = paste0(seq(0, 100, by = 20), "%"),
      na.value = "transparent",
      guide = guide_colorbar(
        direction      = "horizontal",
        label.position = "bottom",
        title.position = "top",
        title.hjust    = 0.5,
        barwidth       = unit(32, "cm"),
        barheight      = unit(2, "cm"),
        title.theme    = element_text(size = 50, face = "bold"),
        label.theme    = element_text(size = 45, face = "bold")
      )
    ) +
    labs(title = title) +
    coord_sf(xlim = xlim, ylim = ylim, clip = "off") +   # <-- FORCED COMMON EXTENT
    theme_minimal() +
    theme(
      plot.title   = element_text(hjust = 0.5, face = "bold", size = 64,
                                  margin = margin(t = title_top_margin)),
      axis.text    = element_blank(),
      axis.title   = element_blank(),
      panel.grid   = element_blank(),
      plot.margin  = margin(0, 0, 0, 0, "pt"),           # zero margin – map fills cell
      legend.title = element_text(face = "bold", size = 42),
      legend.text  = element_text(face = "bold", size = 42)
    )
  
  if (!show_legend) p <- p + theme(legend.position = "none")
  return(p)
}

# ----------------------------------------------------------------------
# 3. BUILD SIDE MAPS (all forced to the same extent)
# ----------------------------------------------------------------------
p_endemic  <- plot_map(r_endemic,  "Iberian Endemic (N=47)",
                       show_legend = FALSE, title_top_margin = 10,
                       xlim = common_xlim, ylim = common_ylim)

p_native   <- plot_map(r_native,   "Native Widespread (N=26)",
                       show_legend = FALSE,
                       xlim = common_xlim, ylim = common_ylim)

p_invasive <- plot_map(r_invasive, "Invasive Widespread (N=25)",
                       show_legend = FALSE,
                       xlim = common_xlim, ylim = common_ylim)

# ----------------------------------------------------------------------
# 4. CENTRAL HUB – extended y‑limit to create internal top margin
# ----------------------------------------------------------------------
central_ylim <- c(common_ylim[1], common_ylim[2] + 0.2 * yrange)  # +20% height

# Build the hub without the label first
p_pooled_hub <- plot_map(r_pooled, title = NULL, show_legend = TRUE,
                         xlim = common_xlim, ylim = central_ylim) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 2),
    plot.margin  = margin(60, 60, 60, 60, "pt")   # outer space for the border
  ) +
  geom_spatvector(data = transect_line_vect, color = "black", alpha = 0.5, linewidth = 2.5) +
  geom_spatvector(data = endpoints_vect, color = "black", size = 3, shape = 21, fill = "black") +
  geom_spatvector_text(data = endpoints_vect, aes(label = Label),
                       fontface = "bold", vjust = -1.2, size = 19,
                       nudge_x = -0.5, nudge_y = -1.1) +
  # Label at a fixed y inside the extended blank area
  geom_text(aes(x = mean(common_xlim), y = common_ylim[2] + 0.1 * yrange,
                label = "All species (N=98)"),
            size = 56 / .pt, fontface = "bold", color = "black",
            hjust = 0.5, vjust = 0.5)

# ----------------------------------------------------------------------
# 5. PROFILE PLOT (unchanged large fonts)
# ----------------------------------------------------------------------

hydro_data <- profile_long %>% filter(Driver == "Hydromorphological")
max_point  <- hydro_data[which.max(hydro_data$Suitability), ]
min_point  <- hydro_data[which.min(hydro_data$Suitability), ]

p_profile_lines <- ggplot(profile_long, aes(x = Distance_km, y = Suitability, color = Driver)) +
  geom_line(linewidth = 1.5) +
  
  # Vertical dashed line at point A (x = 0)
  annotate("segment",
           x = 0, xend = 0,
           y = 0, yend = 85,
           colour = "black", linetype = "dashed", linewidth = 0.8, alpha = 0.7) +
  # Label "A" near the top but inside the plot
  annotate("text",
           x = 0, y = 100,
           label = "A",
           size = 18, colour = "black", fontface = "bold") +
  
  # Vertical dashed line at point B (x = max distance)
  annotate("segment",
           x = max(profile_long$Distance_km), xend = max(profile_long$Distance_km),
           y = 0, yend = 85,
           colour = "black", linetype = "dashed", linewidth = 0.8, alpha = 0.7) +
  # Label "B" near the top but inside the plot
  annotate("text",
           x = max(profile_long$Distance_km), y = 100,
           label = "B",
           size = 18, colour = "black", fontface = "bold") +
  
  # Maximum – green point + black label with "Max = "
  geom_point(data = max_point, size = 5, colour = "#009E73", show.legend = FALSE) +
  geom_text(data = max_point,
            aes(label = paste0("Max = ", round(Suitability, 1), "%")),
            nudge_y = 8, 
            size = 18, colour = "black", fontface = "bold", show.legend = FALSE) +
  
  # Minimum – green point + black label with "Min = "
  geom_point(data = min_point, size = 5, colour = "#009E73", show.legend = FALSE) +
  geom_text(data = min_point,
            aes(label = paste0("Min = ", round(Suitability, 1), "%")),
            nudge_y = -8,
            size = 18, colour = "black", fontface = "bold", show.legend = FALSE) +
  
  scale_color_manual(values = c("Climate" = "#E69F00",
                                "Hydroclimatic" = "#56B4E9",
                                "Hydromorphological" = "#009E73",
                                "Stacked" = "black")) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 25),
                     labels = paste0(seq(0, 100, by = 25), "%")) +
  labs(title = "Freshwater Stacked Habitat Suitability along transect A → B",
       x = "Distance (km)", y = "Stacked Habitat Suitability") +
  theme_minimal() +
  theme(
    plot.margin      = margin(5, 20, 10, 20, "pt"),
    plot.title       = element_text(hjust = 0.5, size = 60, face = "bold", margin = margin(b = 28)),
    axis.title       = element_text(size = 46, face = "bold"),
    axis.title.x     = element_text(size = 46, face = "bold"),
    axis.title.y     = element_text(size = 46, face = "bold"),
    axis.text        = element_text(size = 40),
    axis.text.x      = element_text(size = 40),
    axis.text.y      = element_text(size = 40),
    legend.title     = element_blank(),
    legend.text      = element_text(size = 44),
    legend.key.width = unit(8, "cm"),
    legend.position  = "none"
  ) +
  guides(color = guide_legend(override.aes = list(linewidth = 5)))

# p_profile_lines <- ggplot(profile_long, aes(x = Distance_km, y = Suitability, color = Driver)) +
#   geom_line(linewidth = 1.5) +
#   scale_color_manual(values = c("Climate" = "#E69F00",
#                                 "Hydroclimatic" = "#56B4E9",
#                                 "Hydromorphological" = "#009E73",
#                                 "Stacked" = "black")) +
#   scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 25),
#                      labels = paste0(seq(0, 100, by = 25), "%")) +
#   labs(title = "Freshwater Stacked Habitat Suitability along transect A → B",
#        x = "Distance (km)", y = "Stacked Habitat Suitability") +
#   theme_minimal() +
#   theme(
#     plot.margin      = margin(5, 20, 10, 20, "pt"),
#     plot.title       = element_text(hjust = 0.5, size = 60, face = "bold", margin = margin(b = 28)),
#     axis.title       = element_text(size = 46, face = "bold"),
#     axis.title.x     = element_text(size = 46, face = "bold"),
#     axis.title.y     = element_text(size = 46, face = "bold"),
#     axis.text        = element_text(size = 40),
#     axis.text.x      = element_text(size = 40),
#     axis.text.y      = element_text(size = 40),
#     legend.title     = element_blank(),
#     legend.text      = element_text(size = 44),
#     legend.key.width = unit(8, "cm"),
#     legend.position  = "none"
#   ) +
#   guides(color = guide_legend(override.aes = list(linewidth = 5)))

p_profile_legend <- cowplot::get_legend(
  p_profile_lines + theme(legend.position = "bottom")
)
profile_block <- wrap_plots(p_profile_lines, p_profile_legend, ncol = 1, heights = c(12, 1))

# ----------------------------------------------------------------------
# 6. MAP LAYOUT (equal column widths for all side maps)
# ----------------------------------------------------------------------
side_width <- 1
hub_width  <- 1.6

row1 <- plot_spacer() + p_endemic + plot_spacer() +
  plot_layout(widths = c((side_width + hub_width)/2, side_width, (side_width + hub_width)/2))

row2 <- p_native + p_pooled_hub + p_invasive +
  plot_layout(widths = c(side_width, hub_width, side_width))

top_block <- row1 / row2 +
  plot_layout(heights = c(1.0, 1.6), guides = "collect") &
  theme(
    legend.position      = "bottom",
    legend.justification = "center",
    legend.text          = element_text(size = 34),
    legend.margin        = margin(15, 0, 30, 0, "pt")
  )

base_figure <- wrap_plots(top_block, profile_block, ncol = 1, heights = c(3.2, 1.3)) &
  theme(plot.background = element_rect(fill = "white", color = NA))

# ----------------------------------------------------------------------
# 7. CONNECTING LINES (precise centres)
# ----------------------------------------------------------------------
# Row2 total width: side_width + hub_width + side_width = 3.6
# Centre of native map:   (0.5*side_width) / 3.6 = 0.1389
# Centre of invasive map: 1 - 0.1389 = 0.8611
final_figure <- ggdraw(base_figure) +
  draw_grob(pointsGrob(x = 0.5, y = 0.714, pch = 21,
                       gp = gpar(fill = "black", col = "black", cex = 3))) +
  draw_grob(segmentsGrob(x0 = 0.5, y0 = 0.714, x1 = 0.5, y1 = 0.78,
                         gp = gpar(lwd = 3))) +
  draw_grob(segmentsGrob(x0 = 0.5, y0 = 0.714, x1 = 0.3, y1 = 0.714,
                         gp = gpar(lwd = 3))) +
  draw_grob(segmentsGrob(x0 = 0.5, y0 = 0.714, x1 = 0.74, y1 = 0.714,
                         gp = gpar(lwd = 3)))

# ----------------------------------------------------------------------
# 8. EXPORT
# ----------------------------------------------------------------------
suppressWarnings({
  ggsave(file.path(out_dir, "Final_Aligned_Dashboard_Perfected.png"),
         final_figure,
         width = 48, height = 34, dpi = 300, bg = "white")
})






###################     NOTEPAD   ##################

# 
# 
# #Compare the two stacked (species based vs environmental based pooled layers)
# 
# 
# library(terra)
# 
# # ==============================================================================
# # COMPARE POOLED RASTERS
# # ==============================================================================
# 
# # 1. Load both rasters
# r_new <- rast("/Users/geo_v/Desktop/Vagenas_aSDMs/output/stacked_maps/Pooled_stacked_suitability.tif")
# r_old <- rast("/Users/geo_v/Desktop/Vagenas_aSDMs/output/stacked_maps_alt/Overall_Pooled_Environmental_Suitability.tif")
# 
# # 2. Check if their spatial geometries (extent, resolution, CRS) match
# # This will throw an error if they don't align perfectly
# print("Checking geometry alignment...")
# compareGeom(r_new, r_old)
# 
# # 3. Calculate spatial correlation
# # Stack them and convert to a dataframe for the cor() function
# stack_compare <- c(r_new, r_old)
# names(stack_compare) <- c("New_Multimetric_Pooled", "Old_Simple_Pooled")
# 
# df_compare <- as.data.frame(stack_compare, na.rm = TRUE)
# cor_matrix <- cor(df_compare, use = "complete.obs", method = "pearson")
# 
# print("--- Pearson Correlation Matrix ---")
# print(cor_matrix)
# 
# # 4. Check for absolute identical pixel values
# # If they are exactly identical, every pixel in this difference raster will be 0
# r_diff <- r_new - r_old
# 
# print("--- Summary of Pixel Differences (New - Old) ---")
# # We use minmax to quickly see the range of differences
# print(minmax(r_diff))
# 





#### COMBINE BOTH ANALYSIS INTO ONE ####

# 
# 
# 
# # ==============================================================================
# # 0. SETUP AND LOAD PACKAGES
# # ==============================================================================
# library(dplyr)
# library(readr)
# library(terra)
# library(ggplot2)
# library(tidyterra)
# library(patchwork)
# library(ggcorrplot)
# library(rnaturalearth)
# library(grid) # Required for the connecting lines
# 
# # Directories
# ens_dir <- "Vagenas_aSDMs/output/regional_ensembles_H5_multimetric"
# base_reg_dir <- "Vagenas_aSDMs/output/regional"                      
# out_dir <- "Vagenas_aSDMs/output/stacked_maps_unified"
# 
# if (!dir.exists(out_dir)) {
#   dir.create(out_dir, recursive = TRUE)
# }
# 
# # ==============================================================================
# # 1. PREPARE DATA & SPECIES LISTS
# # ==============================================================================
# species_class <- read_csv("Vagenas_aSDMs/input/aSDMs_species_classification.csv", show_col_types = FALSE) %>%
#   mutate(Species = gsub(" ", "_", Sp))
# 
# # Get valid multimetric species
# valid_files <- list.files(ens_dir, pattern = "_H5_multimetric_ensemble\\.tif$", full.names = TRUE)
# valid_species <- gsub("_H5_multimetric_ensemble\\.tif$", "", basename(valid_files))
# master_df <- species_class %>% filter(Species %in% valid_species)
# species_list <- master_df %>% pull(Species)
# 
# # ==============================================================================
# # 1.5 FETCH, MERGE, AND CLEAN IBERIAN MAINLAND
# # ==============================================================================
# iberia_sf <- ne_countries(country = c("Spain", "Portugal"), scale = "medium", returnclass = "sf")
# iberia_vect <- vect(iberia_sf)
# 
# iberia_dissolved <- aggregate(iberia_vect)
# iberia_parts <- disagg(iberia_dissolved)
# land_areas <- expanse(iberia_parts, unit = "km")
# iberia_mainland <- iberia_parts[which.max(land_areas), ]
# 
# # ==============================================================================
# # 2. RASTER GENERATION FUNCTIONS
# # ==============================================================================
# # FUNCTION A: Standard Sum (For Species Categories and Pooled)
# safe_sum_rasters <- function(file_paths) {
#   r_list <- lapply(file_paths, rast)
#   master_ext <- ext(r_list[[1]])
#   if (length(r_list) > 1) {
#     for (i in 2:length(r_list)) { master_ext <- terra::union(master_ext, ext(r_list[[i]])) }
#   }
#   r_list_aligned <- lapply(r_list, function(r) extend(r, master_ext))
#   r_sum <- sum(rast(r_list_aligned), na.rm = TRUE)
#   
#   iberia_proj <- project(iberia_mainland, crs(r_sum))
#   r_masked <- mask(crop(r_sum, iberia_proj), iberia_proj)
#   
#   rmm <- minmax(r_masked)
#   r_norm <- if (rmm[2,1] > rmm[1,1]) ((r_masked - rmm[1,1]) / (rmm[2,1] - rmm[1,1])) * 100 else r_masked * 0
#   return(r_norm)
# }
# 
# # FUNCTION B: From-Scratch Multimetric Sum (For Environmental Predictors)
# build_multimetric_driver_stack <- function(spp_list, target_driver) {
#   stack_list <- list()
#   for (sp in spp_list) {
#     file_path <- file.path(base_reg_dir, sp, "H5", target_driver, "ensemble_regional.tif")
#     eval_path <- file.path(base_reg_dir, sp, "H5", target_driver, "eval_regional_iberia_strict.csv")
#     
#     if (file.exists(file_path) && file.exists(eval_path)) {
#       metrics <- read_csv(eval_path, show_col_types = FALSE)
#       weight <- mean(metrics$e_AUC + metrics$maxTSS + metrics$CBI + metrics$uAUC, na.rm = TRUE) / 4
#       
#       r <- rast(file_path)
#       rmm <- minmax(r)
#       r_norm <- if (rmm[2,1] > rmm[1,1]) (r - rmm[1,1]) / (rmm[2,1] - rmm[1,1]) else r * 0 
#       
#       stack_list[[sp]] <- r_norm * weight
#     }
#   }
#   
#   if(length(stack_list) == 0) return(NULL)
#   
#   master_ext <- ext(stack_list[[1]])
#   if (length(stack_list) > 1) {
#     for (i in 2:length(stack_list)) { master_ext <- terra::union(master_ext, ext(stack_list[[i]])) }
#   }
#   r_list_aligned <- lapply(stack_list, function(r) extend(r, master_ext))
#   r_sum <- sum(rast(r_list_aligned), na.rm = TRUE)
#   
#   iberia_proj <- project(iberia_mainland, crs(r_sum))
#   r_masked <- mask(crop(r_sum, iberia_proj), iberia_proj)
#   
#   final_rmm <- minmax(r_masked)
#   r_final <- if (final_rmm[2,1] > final_rmm[1,1]) ((r_masked - final_rmm[1,1]) / (final_rmm[2,1] - final_rmm[1,1])) * 100 else r_masked * 0
#   return(r_final)
# }
# 
# # --- GENERATE ALL 7 STACKS ---
# print("Building POOLED Overall Map...")
# r_pooled <- safe_sum_rasters(valid_files)
# 
# print("Building SPECIES Maps...")
# spp_endemic <- master_df %>% filter(Category == "Iberian Endemic") %>% pull(Species)
# r_endemic <- safe_sum_rasters(file.path(ens_dir, paste0(spp_endemic, "_H5_multimetric_ensemble.tif")))
# 
# spp_native <- master_df %>% filter(Category == "Native Widespread") %>% pull(Species)
# r_native <- safe_sum_rasters(file.path(ens_dir, paste0(spp_native, "_H5_multimetric_ensemble.tif")))
# 
# spp_invasive <- master_df %>% filter(Category == "Invasive Widespread") %>% pull(Species)
# r_invasive <- safe_sum_rasters(file.path(ens_dir, paste0(spp_invasive, "_H5_multimetric_ensemble.tif")))
# 
# print("Building PREDICTOR Maps (From Scratch)...")
# r_climate  <- build_multimetric_driver_stack(species_list, "Climate")
# r_hydrocl  <- build_multimetric_driver_stack(species_list, "Hydroclimatic")
# r_hydromo  <- build_multimetric_driver_stack(species_list, "Hydromorphological")
# 
# # ==============================================================================
# # 3. CALCULATE SPATIAL CORRELATIONS (DEFAULT AXIS ALIGNMENT)
# # ==============================================================================
# 
# yellow_to_very_dark_orange <- colorRampPalette(c("yellow", "#993300"))(20) 
# 
# 
# create_corrplot <- function(raster_list, name_list, title) {
#   c_ext <- ext(raster_list[[1]])
#   for (i in 2:length(raster_list)) c_ext <- terra::union(c_ext, ext(raster_list[[i]]))
#   r_list_ext <- lapply(raster_list, function(r) extend(r, c_ext))
#   
#   cor_stack <- rast(r_list_ext)
#   names(cor_stack) <- name_list
#   
#   c_mat <- cor(as.data.frame(cor_stack, na.rm = TRUE), use = "complete.obs", method = "pearson")
#   
#   # Remove the 'colors' argument from ggcorrplot and add scale_fill_gradientn
#   p <- ggcorrplot(c_mat, method = "square", type = "lower", show.diag = FALSE, 
#                   lab = TRUE, lab_size = 5, outline.color = "white", title = title) +
#     
#     # Overwrite the default ggcorrplot fill scale with your custom 100-color palette
#     scale_fill_gradientn(
#       colors = yellow_to_very_dark_orange,
#       limits = c(-1, 1), # Change to c(0, 1) if your correlations are strictly positive!
#       na.value = "transparent"
#     ) +
#     
#     guides(fill = "none") + 
#     theme(
#       plot.title = element_text(hjust = 0.5, face = "bold", size = 22, vjust = -15),
#       axis.text.x.bottom = element_text(angle = 45, hjust = 1, vjust = 1, size = 16, margin = margin(t = 5)),
#       
#       # DEFAULT Y-AXIS SETTING: Let ggplot handle spacing natively
#       axis.text.y = element_text(size = 16), 
#       
#       panel.grid = element_blank(),
#       legend.position = "none",
#       plot.margin = margin(t = 20, r = 5, b = 5, l = 5, unit = "pt")
#     )
#   
#   p$layers[[2]]$aes_params$fontface <- "bold"
#   p$layers[[2]]$data$label <- paste0(round(as.numeric(p$layers[[2]]$data$label) * 100), "%")
#   
#   # Set alpha to something like 0.6 for 40% transparency
#   p$layers[[1]]$aes_params$alpha <- 0.7
#   
#   return(p)
# }
# 
# # 1. Species Corrplot (Left Side)
# # No custom theme overrides needed anymore!
# p_corr_sp <- create_corrplot(
#   list(r_pooled, r_endemic, r_native, r_invasive), 
#   c("Overall", "Endemic", "Native Wid.", "Invasive Wid."), 
#   "Species SI Spatial Correlation (ρ)"
# )
# 
# # 2. Predictor Corrplot (Right Side - Flipped for Symmetry!)
# p_corr_pr <- create_corrplot(
#   list(r_pooled, r_climate, r_hydrocl, r_hydromo), 
#   c("Overall", "Climate", "Hydroclimatic", "Hydromorph."), 
#   "Predictor SI Spatial Correlation (ρ)"
# ) +
#   scale_x_discrete(limits = rev) +             
#   scale_y_discrete(position = "right")
# 
# 
# # ==============================================================================
# # 4. PLOT THE 7 MAPS
# # ==============================================================================
# 
# yellow_to_very_dark_orange <- colorRampPalette(c("yellow", "#993300"))(100) 
# 
# plot_map <- function(r, title) {
#   p <- ggplot() +
#     geom_spatraster(data = r) +
#     geom_spatvector(data = iberia_mainland, fill = NA, color = "black", linewidth = 0.4) +
#     scale_fill_gradientn(
#       colors = yellow_to_very_dark_orange, name = "Stacked Suitability Index (SI, %)",
#       limits = c(0, 100), breaks = seq(0, 100, by = 20), labels = paste0(seq(0, 100, by = 20), "%"), na.value = "transparent",
#       guide = guide_colorbar(
#         direction = "horizontal", label.position = "bottom", title.position = "top", 
#         title.hjust = 0.5, barwidth = unit(20, "cm"), barheight = unit(1.5, "cm")
#       )
#     ) +
#     coord_sf(expand = FALSE) +
#     theme_minimal() +
#     theme(
#       # Increased title size to 24 and added a margin to prevent clipping
#       plot.title = element_text(hjust = 0.5, face = "bold", size = 26, margin = margin(b = 15)), 
#       axis.text = element_blank(), axis.title = element_blank(), panel.grid = element_blank(),
#       plot.margin = margin(5, 5, 5, 5, unit = "pt") 
#     )
#   
#   if (!is.null(title)) {
#     p <- p + labs(title = title)
#   }
#   return(p)
# }
# 
# n_sp <- length(species_list)
# 
# p_endemic  <- plot_map(r_endemic, paste0("Iberian Endemic (N=", length(spp_endemic), ")"))
# p_native   <- plot_map(r_native, paste0("Native Widespread (N=", length(spp_native), ")"))
# p_invasive <- plot_map(r_invasive, paste0("Invasive Widespread (N=", length(spp_invasive), ")"))
# # Predictor titles now explicitly include the N count!
# p_climate  <- plot_map(r_climate, paste0("Climate-Driven Suitability (N=", n_sp, ")"))
# p_hydrocl  <- plot_map(r_hydrocl, paste0("Hydroclimatic-Driven Suitability (N=", n_sp, ")"))
# p_hydromo  <- plot_map(r_hydromo, paste0("Hydromorphological-Driven (N=", n_sp, ")"))
# 
# # ------------------------------------------------------------------------------
# # THE CENTRAL "POOLED" MAP: Adding the Transparent Box & Connecting Lines
# # ------------------------------------------------------------------------------
# p_pooled <- plot_map(r_pooled, NULL) +
#   theme(
#     panel.border = element_rect(color = "black", fill = scales::alpha("black", 0.05), linewidth = 1.2)
#   ) +
#   coord_sf(expand = FALSE, clip = "off") + # <--- UPDATE THIS LINE
#   
#   # Top Connections: ALL start at the center dot (x0=0.5, y0=1.0)
#   # Middle points straight up
#   annotation_custom(segmentsGrob(x0 = 0.5, y0 = 1.0, x1 = 0.5,  y1 = 1.15, gp = gpar(lwd = 1.5))) + 
#   # Left points far left, outside the central map (negative x1)
#   annotation_custom(segmentsGrob(x0 = 0.5, y0 = 1.0, x1 = -0.3, y1 = 1.2, gp = gpar(lwd = 1.5))) + 
#   # Right points far right, outside the central map (x1 > 1.0)
#   annotation_custom(segmentsGrob(x0 = 0.5, y0 = 1.0, x1 = 1.15,  y1 = 1.2, gp = gpar(lwd = 1.5))) +
#   
#   # Central Dot Top
#   annotation_custom(pointsGrob(x = 0.5, y = 1.0, pch = 21, gp = gpar(fill = "black", col = "black", cex = 1.5))) +
#   
#   # Bottom Connections: ALL start at the center dot (x0=0.5, y0=0.0)
#   # Middle points straight down
#   annotation_custom(segmentsGrob(x0 = 0.5, y0 = 0.0, x1 = 0.5,  y1 = -0.15, gp = gpar(lwd = 1.5))) + 
#   # Left points far left
#   annotation_custom(segmentsGrob(x0 = 0.5, y0 = 0.0, x1 = -0.3, y1 = -0.15, gp = gpar(lwd = 1.5))) + 
#   # Right points far right
#   annotation_custom(segmentsGrob(x0 = 0.5, y0 = 0.0, x1 = 1.15,  y1 = -0.15, gp = gpar(lwd = 1.5))) +
#   
#   # Central Dot Bottom
#   annotation_custom(pointsGrob(x = 0.5, y = 0.0, pch = 21, gp = gpar(fill = "black", col = "black", cex = 1.5)))
#   # # Top Connections
#   # annotation_custom(segmentsGrob(x0 = 0.5, y0 = 1.0, x1 = 0.5,  y1 = 1.15, gp = gpar(lwd = 1.5))) +
#   # annotation_custom(segmentsGrob(x0 = 0.1, y0 = 1.0, x1 = -0.2, y1 = 1.17, gp = gpar(lwd = 1.5))) +
#   # annotation_custom(segmentsGrob(x0 = 0.9, y0 = 1.0, x1 = 1.2,  y1 = 1.17, gp = gpar(lwd = 1.5))) +
#   # 
#   # # Bottom Connections (Downwards from y = 0.0 to -0.17)
#   # annotation_custom(segmentsGrob(x0 = 0.5, y0 = 0.0, x1 = 0.5,  y1 = -0.17, gp = gpar(lwd = 1.5))) +
#   # annotation_custom(segmentsGrob(x0 = 0.1, y0 = 0.0, x1 = -0.2, y1 = -0.17, gp = gpar(lwd = 1.5))) +
#   # annotation_custom(segmentsGrob(x0 = 0.9, y0 = 0.0, x1 = 1.2,  y1 = -0.17, gp = gpar(lwd = 1.5)))
#   # 
# # ==============================================================================
# # 5. ARRANGE MASTER 8x4 GRID WITH PATCHWORK
# # ==============================================================================
# layout_design <- c(
#   patchwork::area(t = 1, l = 1, b = 1, r = 2), # Top Left: Endemic
#   patchwork::area(t = 1, l = 4, b = 1, r = 5), # Top Center: Native
#   patchwork::area(t = 1, l = 7, b = 1, r = 8), # Top Right: Invasive
#   
#   patchwork::area(t = 2, l = 1, b = 3, r = 2), # Middle Left: Species Corrplot
#   patchwork::area(t = 2, l = 3, b = 3, r = 6), # THE CENTER HUB: POOLED MAP 
#   patchwork::area(t = 2, l = 7, b = 3, r = 8), # Middle Right: Predictor Corrplot (Symmetrical!)
#   
#   patchwork::area(t = 4, l = 1, b = 4, r = 2), # Bottom Left: Climate
#   patchwork::area(t = 4, l = 4, b = 4, r = 5), # Bottom Center: Hydroclimatic
#   patchwork::area(t = 4, l = 7, b = 4, r = 8)  # Bottom Right: Hydromorphological
# )
# 
# final_figure <- p_endemic + p_native + p_invasive + 
#   p_corr_sp + p_pooled + p_corr_pr + 
#   p_climate + p_hydrocl + p_hydromo +
#   plot_layout(
#     design = layout_design,
#     widths = c(1, 1, 1, 1, 1, 1, 1, 1),  
#     heights = c(1, 1, 1, 1),             # <--- EQUALIZED ROW HEIGHTS
#     guides = "collect"
#   ) & 
#   theme(
#     legend.position = "bottom",
#     legend.title = element_text(face = "bold", size = 26), 
#     legend.text = element_text(size = 18),                 
#     plot.background = element_rect(fill = "white", color = NA),
#     legend.margin = margin(t = 20, unit = "pt") 
#   )
# 
# # Changed dimensions to a perfect square to allow maps to scale horizontally
# ggsave(file.path(out_dir, "Final_Grand_Unified_Suitability.png"), 
#        final_figure, width = 24, height = 19, dpi = 300, bg = "white") # <--- HEIGHT SET TO 22
# 
# print("Complete! Grand Unified Plot successfully generated with full symmetry and box layout.")
# 
# 
# 
# 
# # CORRECT COMBINATION
# 
# 
