<div align="center">

# aSDMs :: aquatic Species Distribution Models

</div>

<img width="2344" height="976" alt="GitHub_Logo" src="https://github.com/user-attachments/assets/be772a43-b132-47c3-8ba9-983a66a731e0" />

# Repository Overview

This repository contains the complete set of analysis scripts, and final outputs (rasters and figures) from our study. Additionally, we begin with a user-friendly tutorial - a conceptual exercise that demonstrates our methodology step-by-step, making it accessible for readers to replicate and explore. The input files are included in this repository either as raw data either stored in public repositories with a link included due large sizes. The six sequential scripts are provided in the form of transferable R files. The main pipeline is included as the fifth .R object (5_pipeline_aSDMs.R).

# Project structure
```
├── 📁 data/          # 🌡️ Local data (bioclimatic vars, hydrological layers etc.)
├── 📁 Scripts        # 💻 Main analysis scripts (see structure below)
├── 📁 outputs/       # 🖨️ Maps and auxillary outputs generated through the present analysis pipeline
└── 📁 figures/       # 🖼️ Final figures included in the manuscript
```
# 🔗 Data
The baseline layers 🌐 required for the analysis can be downloaded from the Input_Layers and the Input_dataset folder here: https://saco.csic.es/s/SYTM8qZrnY2HG5q

# 💻 Scripts
The structure of the scripts for the primary analysis set is structured as: 
```
├── 🌦️ 1_Predictors_aSDMs.R                              # Environmental predictors preparation
├── 🐟 2_Biodata_GBIF_aSDMs.R                            # Species occurrence data extraction
├── 🧬 3_Taxize_Biodata_GBIF_aSDMs.R                     # Taxonomic harmonization and cleaning
├── 🗺️ 4_Figure1_aSDMs.R                                 # Study area and baseline visualizations
├── ⚙️ 5_pipeline_aSDMs.R                                # Main brute-force hierarchical SDM modeling framework
└── 📈 6_PostAnalysis_Fig3_4_plus_Supplementary_aSDMs.R  # Visualization, statistical modeling, ensembling, stacking and final figures 
```

# 📈 Outputs
The repository represents a stand-alone analysis package and contains the full set of initial data and the required script to generate the figures of the study which can be downloaded here: https://saco.csic.es/s/SYTM8qZrnY2HG5q

This repository is produced by the 6_PostAnalysis_Fig3_4_plus_Supplementary_aSDMs.R object.
```
📁 outputs/
├── 📁 figures/                         # 🗺️ All the outputs produced from this study
│   ├── 🌍 Global/                      # 🗺️ All the figure outputs produced for the native & invasive widespread species by using global aSDMs
│   └── 🏝️ Regional/                    # 🗺️ All the figure outputs produced for the native & invasive widespread species by using either global-to-regional aSDMs (native & invasive widespread; look Figure 2 below) or strictly regional aSDMs (endemic species of the study area)
│       ├── ⚙️ Model/                   # 🗺️ LLMs produced to overcome complexity of the outputs in a predictive fashion, not included in the paper
│       ├── 📊 Extent_Performance/      # 🗺️ Analysis to link extent of training, performance metrics and aSDM outputs
│       ├── 🗺️ ensembles/               # 🗺️ All the aSDM outputs for ensembled by using AUC at H5 training extent, ensembles across categories | Includes Figure 4
│       ├── 📈 Predictor_Performance/   # 🗺️ Analysis to link the predictor sets, performance metrics and aSDM outputs
│       └── 📉 VarImp/                  # 🗺️ Results to relate Variable Importance with the predictive outputs | Includes Figure 3
└── 📁 tables/                          # 🔄 Statistics produced through the aforementioned paper
```
---

# 🌊aSDM | Pipeline 🌊
📋 Description

You are an aquatic ecologist/ecohydrological engineer and you want to predict the distribution of multiple endemic species in a hydrographic network. Your task is to build an aquatic Species Distribution Model (aSDM) that predicts suitability scores or species richness. In this demo of the main pipeline (e.g., folder "Scripts" 

<div align="center">
<table>
  <tr>
    <th>Model Architecture</th>
    <th>Symbol</th>
    <th>Predictor Sets</th>
    <th>Short ID</th>
  </tr>
  <tr>
    <td>Climate aSDMs<br><em>(Widespread species)</em></td>
    <td>🌍🌡️</td>
    <td>Global Climate variables</td>
    <td><code>glob_clima_aSDM</code></td>
  </tr>
  <tr>
    <td>Hydroclimatic aSDMs<br><em>(Widespread species)</em></td>
    <td>🌍💧</td>
    <td>Global Hydroclimatic variables</td>
    <td><code>glob_hydroclima_aSDM</code></td>
  </tr>
    <tr>
    <td>Hierarchical Climate aSDMs<br><em>(Widespread & Regional species)</em></td>
    <td>🌡️</td>
    <td>Regional Climate variables<br><em>(+ global climate niche for widespread)</em></td>
    <td><code>reg_clima_aSDM</code></td>
  </tr>
     <tr>
    <td>Hierarchical Hydroclimatic aSDMs<br><em>(Widespread & Regional species)</em></td>
    <td>💧</td>
    <td>Regional Hydroclimatic variables<br><em>(+ global hydroclimatic niche for widespread species)</em></td>
    <td><code>reg_hydroclima_aSDM</code></td>
  </tr>
  <tr>
    <td>Hierarchical Hydromorphological aSDM<br><em>(Widespread & Regional species)</em></td>
    <td>🏞️</td>
    <td>Regional Hydromorphological variables<br><em>(+ global hydroclimatic niche for widespread species)</em></td>
    <td><code>reg_hydromorpho_aSDM</code></td>
  </tr>
</table>
</div>
```

### Start of Tutorial
aSDMs :: Freshwater SDMs – Complete Demo

This tutorial demonstrates the full **aSDMs workflow**:

- **Phase 1** – Global models for widespread species.
- **Phase 2** – Regional models for all species; for widespread ones the global prediction is used as an extra covariate.

The methodology is computational‑light (only 7 random species, few replicates) but exactly follows the production pipeline.  
All required files can be downloaded from: `https://saco.csic.es/s/Co8WNBa323ft3Qi`.


## Step 0 :: Required Packages
```{r}
# List of required packages

pkgs <- c("sf", "sdm", "terra", "dplyr", "tidyr", "mapview", "geodata", "raster",
          "duckdb", "arrow", "vandalico")
if (length(setdiff(pkgs, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(pkgs, rownames(installed.packages())))
}
# Load all packages silently
invisible(lapply(pkgs, library, character.only = TRUE))
```

## Step 1 :: Data Input & Species Filtering
```{r}


# set your working directory (adjust path)
setwd("path/to/your/project")

# 1.1 Occurrence data
iberia_occ <- vect("enriched_iberia_vect.shp")
global_occ <- vect("global_gbif_widespread.shp")

# 1.2 Read the species dataset, a set of points regional points is provided in the data folder (treat them as Endemics) but for full scale application global (GBIF) and regional (local databases) are needed.
species_list <- read.csv("enriched_species_grid_counts.csv")

# Keep species with ≥10 occurrences
adequate <- species_list %>% filter(Total_Grids >= 10)
robust_sp <- adequate$Sp

# subset vectors
iberia_occ <- iberia_occ[iberia_occ$Sp %in% robust_sp, ]
global_occ <- global_occ[global_occ$Sp %in% robust_sp, ]

# separate widespread vs. endemic
widespread_sp <- intersect(unique(iberia_occ$Sp), unique(global_occ$Sp))
endemic_sp    <- setdiff(unique(iberia_occ$Sp), widespread_sp)

cat("Widespread species:", length(widespread_sp), "\n")
cat("Endemic species:", length(endemic_sp), "\n")

# For demo, keep only a few species (optional)
set.seed(123)
demo_widespread <- sample(widespread_sp, 2)  # 2 widespread
demo_endemic    <- sample(endemic_sp, 5)      # 5 endemic
iberia_occ <- iberia_occ[iberia_occ$Sp %in% c(demo_widespread, demo_endemic), ]
global_occ <- global_occ[global_occ$Sp %in% demo_widespread, ]
widespread_sp <- demo_widespread
endemic_sp    <- demo_endemic
```

## Step 2 :: Predictors & Spatial Extents
```{r}

# 2.1 Predictor rasters
predictors_global   <- rast("predictors_finalized_global.tiff")
predictors_regional <- rast("predictors_finalized_regional.tiff")

# 2.2 Define predictor sets (must match raster layer names)
# Note: Regional Climate uses bio6_clima instead of bio17_clima (varies by region)
sets_global <- list(
  "Climate"       = c("bio5_clima", "bio16_clima", "bio17_clima", "bio15_clima", "bio4_clima"),
  "Hydroclimatic" = c("bio4_hydro", "bio1_hydro", "bio16_hydro", "bio17_hydro", "bio15_hydro")
)

sets_regional <- list(
  "Climate"            = c("bio5_clima", "bio16_clima", "bio6_clima", "bio15_clima", "bio4_clima"),
  "Hydroclimatic"      = c("bio4_hydro", "bio1_hydro", "bio16_hydro", "bio17_hydro", "bio15_hydro"),
  "Hydromorphological" = c("lka_pc_use", "dor_pc_pva", "sgr_dk_rav", "urb_pc_use", "for_pc_use")
)

# 2.3 Spatial extents, all the files are provided
ecoregions    <- vect("feow_hydrosheds.shp")
hydrosheds_H5 <- vect("HydroSHEDS_H5_merged.shp")  # all H5 basins
hydrosheds_H8 <- vect("HydroSHEDS_H8_merged.shp")
hydrosheds_H12<- vect("HydroSHEDS_H12_merged.shp")
study_area    <- vect("study_area_iberia.shp")

# 2.4 Training extent generation function
generate_extent <- function(occ, sp_name, eco, hydro, crop_area = NULL) {
  sp_occ <- occ[occ$Sp == sp_name, ]
  pts_eco <- terra::intersect(sp_occ, eco)
  eco_ids <- unique(pts_eco$FEOW_ID)
  sp_eco  <- eco[eco$FEOW_ID %in% eco_ids, ]
  pts_h   <- terra::intersect(sp_occ, hydro)
  h_ids   <- unique(pts_h$HYBAS_ID)
  sp_h    <- hydro[hydro$HYBAS_ID %in% h_ids, ]
  sp_poly <- crop(sp_h, sp_eco)
  sp_poly <- mask(sp_poly, sp_eco)
  sp_poly <- aggregate(sp_poly)
  if (!is.null(crop_area)) sp_poly <- crop(sp_poly, crop_area)
  return(sp_poly)
}


```

## Step 3 :: Build the Parquet Database (faster extraction, does not overloads memory)
```{r}
# 3.1 Global Parquet (for Phase 1)
global_parquet_file <- "aSDMs_Global_Parquet.parquet"

if (!file.exists(global_parquet_file)) {
  env_df <- as.data.frame(predictors_global, cells = TRUE, xy = TRUE, na.rm = TRUE)
  names(env_df) <- tolower(names(env_df))
  names(env_df) <- gsub("\\.", "_", names(env_df))
  write_parquet(env_df, global_parquet_file)
  cat("Global Parquet database created.\n")
} else {
  cat("Using existing Global Parquet file.\n")
}

# 3.2 Regional Parquet (for Phase 2) – contains hydromorphological variables
regional_parquet_file <- "aSDMs_Regional_Parquet.parquet"

if (!file.exists(regional_parquet_file)) {
  env_df_reg <- as.data.frame(predictors_regional, cells = TRUE, xy = TRUE, na.rm = TRUE)
  # ensure the column names match the raster layer names exactly
  layer_names <- names(predictors_regional)
  if (ncol(env_df_reg) == length(layer_names) + 2) {
    names(env_df_reg)[3:ncol(env_df_reg)] <- layer_names
  }
  names(env_df_reg) <- tolower(names(env_df_reg))
  names(env_df_reg) <- gsub("\\.", "_", names(env_df_reg))
  if (!"cell" %in% names(env_df_reg)) env_df_reg$cell <- 1:nrow(env_df_reg)
  write_parquet(env_df_reg, regional_parquet_file)
  cat("Regional Parquet database created.\n")
} else {
  cat("Using existing Regional Parquet file.\n")
}
```

## Step 4 :: Helper Function (Spatial Query + Background)
```{r}

# Generic function to fetch training data from a given Parquet file
get_train_data <- function(sp_name, sp_ext_poly, occ_df, t_vars, parquet_path, n_bg = NULL) {
  # 4.1 Get bounding box of the polygon
  sp_ext_coords <- as.vector(ext(sp_ext_poly))
  xmin <- sp_ext_coords['xmin']; xmax <- sp_ext_coords['xmax']
  ymin <- sp_ext_coords['ymin']; ymax <- sp_ext_coords['ymax']
  
  # 4.2 Query Parquet inside the bounding box
  con <- dbConnect(duckdb::duckdb())
  cols_sql <- paste(paste0("e.", t_vars), collapse = ", ")
  q <- sprintf("SELECT e.cell, e.x, e.y, %s FROM read_parquet('%s') e WHERE e.x BETWEEN %f AND %f AND e.y BETWEEN %f AND %f",
               cols_sql, parquet_path, xmin, xmax, ymin, ymax)
  box_data <- dbGetQuery(con, q)
  dbDisconnect(con)
  names(box_data) <- tolower(names(box_data))
  
  # 4.3 Mask to the exact polygon
  box_pts <- vect(box_data, geom = c("x", "y"), crs = crs(sp_ext_poly), keepgeom = TRUE)
  inside_pts <- box_pts[sp_ext_poly, ]
  poly_data <- as.data.frame(inside_pts)
  
  # 4.4 Presence points
  pres_ids <- occ_df[occ_df$Sp == sp_name, "cell_id"]
  valid_pres <- intersect(pres_ids, poly_data$cell)
  df_pres <- poly_data[poly_data$cell %in% valid_pres, ]
  df_pres$presence <- 1
  
  # 4.5 Background sampling (strictly from inside polygon)
  if (is.null(n_bg)) n_bg <- max(5, round(nrow(poly_data) * 0.05))
  set.seed(123)
  bg_cells <- sample(setdiff(poly_data$cell, valid_pres), n_bg, replace = FALSE)
  df_bg <- poly_data[poly_data$cell %in% bg_cells, ]
  df_bg$presence <- 0
  
  train_ready <- rbind(df_pres, df_bg)
  return(train_ready)
}


```


## Step 5 :: Prepare Output Structures
```{r}

global_output_dir <- "phase1_global"
regional_output_dir <- "phase2_regional"
dir.create(global_output_dir, recursive = TRUE)
dir.create(regional_output_dir, recursive = TRUE)

# Results table
results_all <- data.frame(
  species   = character(),
  phase     = character(),
  extent    = character(),
  set       = character(),
  AUC       = numeric(),
  COR       = numeric(),
  maxTSS    = numeric(),
  maxKappa  = numeric(),
  CBI       = numeric(),
  stringsAsFactors = FALSE
)

algo_list <- c('glm', 'brt', 'rf')
n_reps_demo <- 2   # low for speed (production uses 5)
ref_raster <- predictors_regional[[1]]

```


## Step 6 :: Phase 1 – Global Models (widespread species)
```{r}

extents_global <- list(eco = ecoregions, H5 = hydrosheds_H5, H8 = hydrosheds_H8, H12 = hydrosheds_H12)

for (sp in widespread_sp) {
  clean_sp <- gsub(" ", "_", sp)
  sp_global_occ <- as.data.frame(global_occ[global_occ$Sp == sp, ])
  # compute cell IDs using the global reference raster
  sp_global_occ$cell_id <- cellFromXY(predictors_global[[1]], crds(global_occ[global_occ$Sp == sp, ]))
  
  for (ext_name in names(extents_global)) {
    ext_poly <- generate_extent(global_occ, sp, ecoregions, extents_global[[ext_name]])
    if (is.null(ext_poly) || nrow(ext_poly) == 0) next
    
    for (set_name in names(sets_global)) {
      cat(sprintf("GLOBAL | %s | %s | %s\n", sp, ext_name, set_name))
      out_path <- file.path(global_output_dir, clean_sp, ext_name, set_name)
      dir.create(out_path, recursive = TRUE, showWarnings = FALSE)
      
      t_vars <- sets_global[[set_name]]
      # Use global Parquet for Phase 1
      train_ready <- get_train_data(sp, ext_poly, sp_global_occ, t_vars, 
                                    parquet_path = global_parquet_file)
      
      d <- sdmData(as.formula(paste0("presence ~ ", paste(t_vars, collapse = "+"), " + coords(x+y)")),
                   train = train_ready[train_ready$presence == 1, ],
                   bg    = train_ready[train_ready$presence == 0, ])
      
      m <- tryCatch({
        sdm(presence ~ ., d, methods = algo_list, replication = 'boot', n = n_reps_demo)
      }, error = function(e) { NULL })
      if (is.null(m)) next
      
      p_ens <- ensemble(m, train_ready[, c("x", "y", t_vars)],
                        setting = list(method = 'weighted', stat = 'AUC'))
      res_df <- data.frame(x = train_ready$x, y = train_ready$y, val = as.numeric(p_ens[[1]]))
      global_raster <- rast(res_df, type = "xyz", crs = crs(ext_poly))
      
      # crop to Iberia and save
      iberia_masked <- mask(crop(global_raster, study_area), study_area)
      writeRaster(iberia_masked, file.path(out_path, "ensemble_global_iberia.tif"), overwrite = TRUE)
      
      # evaluation
      e <- evaluates(d, p_ens)
      b <- sdm:::.boyce(e@observed, e@predicted)
      metrics <- data.frame(
        species = sp, phase = "global", extent = ext_name, set = set_name,
        AUC = as.numeric(e@statistics$AUC[1]),
        COR = as.numeric(e@statistics$COR[1]),
        maxTSS = max(e@threshold_based$TSS),
        maxKappa = max(e@threshold_based$Kappa),
        CBI = ifelse(!is.null(b$CBI), as.numeric(b$CBI[1]), NA)
      )
      results_all <- rbind(results_all, metrics)
      write.csv(metrics, file.path(out_path, "eval.csv"), row.names = FALSE)
      
      # variable importance
      vi <- getVarImp(m, id = "ensemble")
      write.csv(vi, file.path(out_path, "varImp.csv"))
      
      rm(d, m, global_raster, iberia_masked, p_ens); gc()
    }
  }
}


```


## Step 7 :: Phase 2 – Regional Models (all regional species)
```{r}
extents_regional <- list(eco = ecoregions, H5 = hydrosheds_H5, H8 = hydrosheds_H8, H12 = hydrosheds_H12)

# Pre‑compute cell IDs for Iberian occurrences (using regional reference raster)
iberia_occ_df <- as.data.frame(iberia_occ)
iberia_occ_df$cell_id <- cellFromXY(ref_raster, crds(iberia_occ))

for (sp in c(widespread_sp, endemic_sp)) {
  clean_sp <- gsub(" ", "_", sp)
  sp_occ <- iberia_occ_df[iberia_occ_df$Sp == sp, ]
  
  for (ext_name in names(extents_regional)) {
    ext_poly <- generate_extent(iberia_occ, sp, ecoregions, extents_regional[[ext_name]], study_area)
    if (is.null(ext_poly) || nrow(ext_poly) == 0) next
    
    for (set_name in names(sets_regional)) {
      cat(sprintf("REGIONAL | %s | %s | %s\n", sp, ext_name, set_name))
      out_path <- file.path(regional_output_dir, clean_sp, ext_name, set_name)
      dir.create(out_path, recursive = TRUE, showWarnings = FALSE)
      
      t_vars <- sets_regional[[set_name]]
      
      # Use regional Parquet for Phase 2
      train_ready <- get_train_data(sp, ext_poly, sp_occ, t_vars, 
                                    parquet_path = regional_parquet_file)
      
      # inject global suitability if species is widespread
      if (sp %in% widespread_sp) {
        global_ras <- rast(file.path(global_output_dir, clean_sp, ext_name,
                                     ifelse(set_name %in% c("Climate","Hydroclimatic"),
                                            set_name, "Hydroclimatic"),
                                     "ensemble_global_iberia.tif"))
        if (file.exists(global_ras)) {
          global_vals <- terra::extract(global_ras, cbind(train_ready$x, train_ready$y))[, 2]
          train_ready$global_suit <- ifelse(is.na(global_vals), 0, global_vals)
          t_vars <- c(t_vars, "global_suit")
        }
      }
      
      d <- sdmData(as.formula(paste0("presence ~ ", paste(t_vars, collapse = "+"), " + coords(x+y)")),
                   train = train_ready[train_ready$presence == 1, ],
                   bg    = train_ready[train_ready$presence == 0, ])
      
      m <- tryCatch({
        sdm(presence ~ ., d, methods = algo_list, replication = 'boot', n = n_reps_demo)
      }, error = function(e) { NULL })
      if (is.null(m)) next
      
      p_ens <- ensemble(m, train_ready[, c("x", "y", t_vars)],
                        setting = list(method = 'weighted', stat = 'AUC'))
      reg_raster <- rast(data.frame(x = train_ready$x, y = train_ready$y, val = as.numeric(p_ens[[1]])),
                         type = "xyz", crs = crs(ext_poly))
      writeRaster(reg_raster, file.path(out_path, "ensemble_regional.tif"), overwrite = TRUE)
      
      e <- evaluates(d, p_ens)
      b <- sdm:::.boyce(e@observed, e@predicted)
      metrics <- data.frame(
        species = sp, phase = "regional", extent = ext_name, set = set_name,
        AUC = as.numeric(e@statistics$AUC[1]),
        COR = as.numeric(e@statistics$COR[1]),
        maxTSS = max(e@threshold_based$TSS),
        maxKappa = max(e@threshold_based$Kappa),
        CBI = ifelse(!is.null(b$CBI), as.numeric(b$CBI[1]), NA)
      )
      results_all <- rbind(results_all, metrics)
      write.csv(metrics, file.path(out_path, "eval.csv"), row.names = FALSE)
      
      vi <- getVarImp(m, id = "ensemble")
      write.csv(vi, file.path(out_path, "varImp.csv"))
      
      rm(d, m, reg_raster, p_ens); gc()
    }
  }
}


```


## Step 8 :: Final Summary
```{r}
write.csv(results_all, "full_evaluation_summary.csv", row.names = FALSE)
cat("Pipeline completed. Results saved to:", global_output_dir, "and", regional_output_dir, "\n")
cat("Evaluation summary written to full_evaluation_summary.csv\n")

```

### End of Tutorial

# Manuscript Outline

## Abstract:
Aim: Species Distribution Models (SDMs) have traditionally been developed in a terrestrial context, and their application to aquatic ecosystems presents unique challenges. These include predicting species distributions across spatially constrained environments, and incorporating specific environmental drivers such as hydromorphological features. We address these challenges by exploring various spatially-explicit model training strategies, novel hierarchical model structures and different flexible predictor sets. Our goal is to establish a framework for modelling aquatic species distributions that accounts for the distinct biogeography of freshwater systems.

Innovation: We demonstrate the application of a novel analytical framework to model aquatic species distributions by investigating three critical dimensions: (i) the optimal spatial training extent for SDMs in freshwater ecosystems; (ii) the effect of multiple predictor combinations—from single variables to hierarchical sets integrating climatic, hydrological, and interactive factors—on predictive performance; and (iii) the interplay between climate and hydrology in predicting species distributions from global to regional scales. This systematic cross-comparison forms the core of our methodological framework for assessing the predictive capacity across all model configurations.

Main conclusions: We used the Iberian freshwater fish as a case study. Our results show that pre-constrained models to a species' watershed of occurrence—meaning they are trained within its boundaries—deliver higher predictive accuracy, effectively validating a spatially-constrained modelling strategy as the best practice. Furthermore, our results demonstrate that climate-based predictors consistently outperform purely hydrological ones, a finding with broad relevance for understanding freshwater species' sensitivity to large-scale environmental change. Our framework provides a basis for standardising SDMs in freshwater systems, with the rigorous protocol established serving as a foundational model for global applications.

### Keywords: 
SDMs, freshwaters, fish, hydrology, climate, watersheds, hierarchical, aquatic species

#### Citation (APA):
Vagenas, G., Matias, M., Araujo M.B. (2026). A hydroclimatic framework for the hierarchical modelling of aquatic species distributions. (Under Revision) 

#### DOI:  
[Pending]

# Figures

<img width="3561" height="1965" alt="Figure1" src="https://github.com/user-attachments/assets/f86acc11-8a29-40fd-bfc3-ea5c4ed98de2" />

**Figure 1**. Spatial distribution of species richness of the dataset used for the development of the aSDMs through a (A) global (50 arc-minute grid) to (B) regional (10 arc-minute grid) approach for the 98 freshwater fish species of the study area. Colors represent gradients of species richness (low = yellow; high = red). The finer resolution map (B) highlights richness patterns in the Iberian Peninsula.

<img width="4000" height="2250" alt="Figure2_JPG_final" src="https://github.com/user-attachments/assets/49addead-41f6-4943-86ea-eb7e25e808ae" />

**Figure 2.** Flowchart illustrating the implementation and evaluation workflow for aquatic Species Distribution Models (aSDMs), comprising nine sequential stages (i.e., I-IX), from input data preparation and modelling through to performance evaluation and the generation of stacked suitability maps.

<img width="4000" height="2250" alt="Figure3_JPG_final" src="https://github.com/user-attachments/assets/d5c72037-30e0-430c-9a73-287eaa26f638" />

**Figure 3.** Variable performarnce across different training extents and predictor settings for the freshwater fish species of the Iberian peninsula.

<img width="3602" height="2250" alt="Figure4_JPG_final" src="https://github.com/user-attachments/assets/3fc8342f-6218-4679-b135-becbc0cc8ad3" />

**Figure 4.** Stacked ensembled aSDMs for the freshwater fish species of the Iberian Peninsula. The maps represent stacked outputs derived through aSDMs using the pre-constrained h5 spatial strategy, by ensembling all the three predictor sets (i.e., climate, hydroclimatic, hydromorphology). The bottom distance-suitability trajectory chart indicates the variation of predicted suitability values across a vertical transect of the study area, indicating the baseline patterns for the thermal (orange), the hydrological (blue) and the locally influenced (green) niche for the freshwater species.

# Author: Georgios Vagenas

Name: PhD Researcher - Georgios Vagenas (georgios.vagenas@mncn.csic.es | georgvagenas@gmail.com)

Affiliation: Biogeography and Global Change Department, National Museum of Natural Sciences, CSIC, C/ Jose Gutierrez Abascal, 2, Madrid 28006, Spain

**Last modified: 15/6/2026**


