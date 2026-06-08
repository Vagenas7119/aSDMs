#June 2026
#Development pipeline of freshwater SDMs (aSDMS)

#Sector: Pipeline execution

#PhD Researcher - Georgios Vagenas (georgios.vagenas@mncn.csic.es | georgvagenas@gmail.com)

#### Pre-setting :: Libraries required to perform the analysis ####

library(terra)
library(readr)
library(tidyr)
library(dplyr)
library(sf)
library(raster)
library(ggplot2)
library(dplyr)
library(readxl)

global_gbif_widespread<-vect("/Users/georgevagenas/Desktop/Vagenas_aSDMs/input/Dataset/Final_enriched_GBIF_MITECO_SNIPAD_tax_correct_objects/global_gbif_widespread.shp")
unique(global_gbif_widespread$Sp)
unique(global_gbif_widespread$Source_DB)
global_gbif_widespread

enriched_iberia_vect<-vect("/Users/georgevagenas/Desktop/Vagenas_aSDMs/input/Dataset/Final_enriched_GBIF_MITECO_SNIPAD_tax_correct_objects/enriched_iberia_vect.shp")
unique(enriched_iberia_vect$Sp)
unique(enriched_iberia_vect$Category)
enriched_iberia_vect

species_list<-read.csv("/Users/georgevagenas/Desktop/Vagenas_aSDMs/input/Dataset/Final_enriched_GBIF_MITECO_SNIPAD_tax_correct_objects/enriched_species_grid_counts.csv")
str(species_list)

predictors_regional<-rast("/Users/georgevagenas/Desktop/Vagenas_aSDMs/input/Predictors/predictors_finalized_regional/predictors_finalized_regional.tiff")
predictors_global<-rast("/Users/georgevagenas/Desktop/Vagenas_aSDMs/input/Predictors/predictors_finalized_global/predictors_finalized_global.tiff")



# ==============================================================================
# 1. IDENTIFY SPECIES WITH < 10 GRIDS & CREATE TABLES
# ==============================================================================

# Find species with more than 10 grids
adequate_species <- species_list %>% 
  dplyr::filter(Total_Grids >=10)

# TABLE 1: The exact list of species, their category, and their grid count
cat("\n--- SPECIES WITH > 10 GRIDS ---\n")
print(adequate_species %>% dplyr::select(Sp, Category, Total_Grids))

# TABLE 2: Summary of how many species fall into each category
cat("\n--- SUMMARY BY CATEGORY ---\n")
rare_summary <- adequate_species %>% 
  dplyr::group_by(Category) %>% 
  dplyr::summarise(Number_of_Species = n())
print(rare_summary)

# ==============================================================================
# 2. SUBSET IBERIAN VECTOR (Keep species with >= 10 grids)
# ==============================================================================

# Get the list of robust species to KEEP
robust_species <- species_list %>% 
  dplyr::filter(Total_Grids >= 10) %>% 
  dplyr::pull(Sp)
unique(robust_species)
# Subset the SpatVector
thres_iberian_vect <- terra::subset(enriched_iberia_vect, enriched_iberia_vect$Sp %in% robust_species)

# ==============================================================================
# 3. SUBSET GLOBAL VECTOR (Keep widespread robust species)
# ==============================================================================

# Identify which of the robust species are actually in the global dataset
robust_global_species <- robust_species[robust_species %in% unique(global_gbif_widespread$Sp)]
unique(robust_global_species)
# Subset the global SpatVector
thres_global_gbif <- global_gbif_widespread[global_gbif_widespread$Sp %in% robust_global_species, ]
unique(thres_global_gbif$Sp)


unique(thres_iberian_vect$Sp)
unique(thres_global_gbif$Sp)

# ==============================================================================
# 4. Overlaps
# ==============================================================================

library(dplyr)

# 1. Get the unique species names currently in your thresholded global vector
global_sp_present <- unique(thres_global_gbif$Sp)

# 2. Map categories from species_list and evaluate inclusion status
widespread_audit <- species_list %>%
  # NEW: Only check species that passed the 10-grid threshold
  dplyr::filter(Total_Grids >= 10) %>%
  # Keep only the widespread categories
  dplyr::filter(Category %in% c("Native Widespread", "Invasive Widespread")) %>%
  # Check if the species made it into the thresholded global vector
  dplyr::mutate(In_Global = Sp %in% global_sp_present) %>%
  dplyr::arrange(Category, desc(In_Global), Sp)

# ==============================================================================
# PRINT AUDIT RESULTS
# ==============================================================================

cat("\n======================================================================\n")
cat("          WIDESPREAD SPECIES INCLUSION AUDIT (IBERIA vs GLOBAL)         ")
cat("\n======================================================================\n")

# Summary of Counts
summary_counts <- widespread_audit %>%
  dplyr::group_by(Category, In_Global) %>%
  dplyr::summarise(Count = n(), .groups = 'drop') %>%
  dplyr::mutate(Status = ifelse(In_Global, "Included in thres_global_gbif", "MISSING from thres_global_gbif"))

print(summary_counts %>% dplyr::select(Category, Status, Count))

cat("\n----------------------------------------------------------------------\n")
cat(" Detail: List of Widespread Species NOT INCLUDED in Global Dataset")
cat("\n----------------------------------------------------------------------\n")

missing_species <- widespread_audit %>% dplyr::filter(!In_Global)
if(nrow(missing_species) > 0) {
  print(missing_species %>% dplyr::select(Category, Sp, Total_Grids))
} else {
  cat("Success! All thresholded widespread species are perfectly represented globally.\n")
}


#all the 51 are of the enriched are included in the global


########### PIPELINE ########### 

#[INPUT 1/4]

#species occurences, break the vectors to dataframes

#Regional dataset
# 1. Extract the raw X/Y coordinates from the Iberian SpatVector

thres_iberian_df<-as.data.frame(thres_iberian_vect)
iberian_coords <- terra::crds(thres_iberian_vect)

# 2. OVERWRITE the cell IDs using the GLOBAL reference raster
# This guarantees perfect alignment with the Parquet database
ref_raster <- predictors_regional[[1]]
thres_iberian_df$cell_id <- terra::cellFromXY(ref_raster, iberian_coords)

# 3. Clean up any points that fall outside the raster (NA cells)
thres_iberian_df <- thres_iberian_df[!is.na(thres_iberian_df$cell_id), ]

# Verify it looks correct
str(thres_iberian_df)

#Global dataset
# 1. Extract the raw X/Y coordinates from the global SpatVector
thres_global_df<-as.data.frame(thres_global_gbif)
global_coords <- terra::crds(thres_global_gbif)

# 2. Calculate the cell IDs using your reference raster
ref_raster <- predictors_global[[1]]
thres_global_df$cell_id <- terra::cellFromXY(ref_raster, global_coords)

# 3. Clean up any points that fall in the ocean/outside the raster (NA cells)
thres_global_df <- thres_global_df[!is.na(thres_global_df$cell_id), ]

#Check to make sure it matches the Iberian structure now
str(thres_global_df)


#predictors
predictors_regional<-rast("/Users/georgevagenas/Desktop/Vagenas_aSDMs/input/Predictors/predictors_finalized_regional/predictors_finalized_regional.tiff")
predictors_global<-rast("/Users/georgevagenas/Desktop/Vagenas_aSDMs/input/Predictors/predictors_finalized_global/predictors_finalized_global.tiff")

#Spatial extent preparation (Ecoregion-H5-H8-H12)

library(terra)
library(dplyr)

# 1. Ecoregions
ecoregions <- vect("/Users/georgevagenas/Desktop/Vagenas_Global_aSDMs/input/SpatialExtents/ecoregions/feow_hydrosheds.shp")

# 2. Hydrosheds H5 (from your previous code)
folder_h5 <- "/Users/georgevagenas/Desktop/Vagenas_aSDMs/input/SpatialExtents/HydroSHEDS/H5/"
h5_files <- list.files(folder_h5, pattern = "\\.shp$", full.names = TRUE, recursive = TRUE)
hydrosheds_H5 <- do.call(rbind, lapply(h5_files, vect))

#dimensions  : 5393, 15  (geometries, attributes)

# 3. Hydrosheds H8
folder_h8 <- "/Users/georgevagenas/Desktop/Vagenas_aSDMs/input/SpatialExtents/HydroSHEDS/H8/"
h8_files <- list.files(folder_h8, pattern = "\\.shp$", full.names = TRUE, recursive = TRUE)
hydrosheds_H8 <- do.call(rbind, lapply(h8_files, vect))
#dimensions  : 221699, 15  (geometries, attributes)

# 4. Hydrosheds H12
folder_h12 <- "/Users/georgevagenas/Desktop/Vagenas_aSDMs/input/SpatialExtents/HydroSHEDS/H12/"
h12_files <- list.files(folder_h12, pattern = "\\.shp$", full.names = TRUE, recursive = TRUE)
hydrosheds_H12 <- do.call(rbind, lapply(h12_files, vect))
#dimensions  : 1242561, 15  (geometries, attributes)

generate_training_extent <- function(occ_points, eco_layer, hydro_layer = NULL) {
  
  # 1. Get Ecoregion intersections for all points at once
  message("Intersecting points with Ecoregions...")
  pts_eco <- terra::intersect(occ_points, eco_layer)
  
  # 2. If a hydro layer is provided, get Hydro intersections
  if (!is.null(hydro_layer)) {
    message("Intersecting points with Hydrosheds...")
    pts_hydro <- terra::intersect(occ_points, hydro_layer)
  }
  
  # 3. Loop through each unique species
  sp_list <- unique(occ_points$Sp)
  final_polys <- list()
  
  for (sp in sp_list) {
    message("Processing: ", sp)
    
    # Get unique Ecoregions for this species
    f_ids <- unique(pts_eco$FEOW_ID[pts_eco$Sp == sp])
    sp_eco_poly <- eco_layer[eco_layer$FEOW_ID %in% f_ids, ]
    
    # SCENARIO A: Ecoregion Only
    if (is.null(hydro_layer)) {
      sp_poly_diss <- aggregate(sp_eco_poly)
      sp_poly_diss$Sp <- sp
      final_polys[[sp]] <- sp_poly_diss
      
      # SCENARIO B: Hydrosheds + Ecoregion Mask  
    } else {
      h_ids <- unique(pts_hydro$HYBAS_ID[pts_hydro$Sp == sp])
      sp_hydro_poly <- hydro_layer[hydro_layer$HYBAS_ID %in% h_ids, ]
      
      # Crop and mask hydrosheds strictly to the ecoregion boundaries
      hydro_in_eco <- crop(sp_hydro_poly, sp_eco_poly)
      hydro_in_eco <- mask(hydro_in_eco, sp_eco_poly)
      
      # Dissolve into a single polygon
      sp_poly_diss <- aggregate(hydro_in_eco)
      sp_poly_diss$Sp <- sp
      final_polys[[sp]] <- sp_poly_diss
    }
  }
  
  # Bind all species polygons together and return
  message("Binding final SpatVector...")
  return(vect(final_polys))
}


out_dir <- "/Users/georgevagenas/Desktop/Vagenas_aSDMs/output/training_ext/"

# ==============================================================================
# 1. REGIONAL DATASET (thres_iberian_vect)
# ==============================================================================

# Regional Ecoregion Only
ext_regional_eco <- generate_training_extent(thres_iberian_vect, ecoregions, hydro_layer = NULL)
study_area_no_islands<-vect("/Users/georgevagenas/Desktop/Vagenas_aSDMs/input/SpatialExtents/study_area_iberia.shp")
#plot(study_area_no_islands)
ext_regional_eco<-terra::crop(ext_regional_eco,study_area_no_islands)
terra::crs(ext_regional_eco)<-"EPSG:4326"
#writeVector(ext_regional_eco, file.path(out_dir, "ext_of_training_regional_ecoregion.gpkg"), overwrite = TRUE)
ext_regional_eco<-vect("/Users/georgevagenas/Desktop/Vagenas_aSDMs/output/training_ext/ext_of_training_regional_ecoregion.gpkg")


# example
library(mapview)
unique(ext_regional_eco$Sp)
#mapview(ext_regional_eco[14])

# Regional H5
#ext_regional_H5 <- generate_training_extent(thres_iberian_vect, ecoregions, hydrosheds_H5)
ext_regional_H5<-terra::crop(ext_regional_H5,study_area_no_islands)
terra::crs(ext_regional_H5)<-"EPSG:4326"
#writeVector(ext_regional_H5, file.path(out_dir, "ext_of_training_regional_H5.gpkg"), overwrite = TRUE)
ext_regional_H5<-vect("/Users/georgevagenas/Desktop/Vagenas_aSDMs/output/training_ext/ext_of_training_regional_H5.gpkg")


#example
unique(ext_regional_H5$Sp)
#mapview(ext_regional_H5[14])


# Regional H8
#ext_regional_H8 <- generate_training_extent(thres_iberian_vect, ecoregions, hydrosheds_H8)
ext_regional_H8<-terra::crop(ext_regional_H8,study_area_no_islands)
terra::crs(ext_regional_H8)<-"EPSG:4326"
#writeVector(ext_regional_H8, file.path(out_dir, "ext_of_training_regional_H8.gpkg"), overwrite = TRUE)
ext_regional_H8<-vect("/Users/georgevagenas/Desktop/Vagenas_aSDMs/output/training_ext/ext_of_training_regional_H8.gpkg")

#example
unique(ext_regional_H8$Sp)
#mapview(ext_regional_H8[14])

# Regional H12
#ext_regional_H12 <- generate_training_extent(thres_iberian_vect, ecoregions, hydrosheds_H12)
ext_regional_H12<-terra::crop(ext_regional_H12,study_area_no_islands)
terra::crs(ext_regional_H12)<-"EPSG:4326"
#writeVector(ext_regional_H12, file.path(out_dir, "ext_of_training_regional_H12.gpkg"), overwrite = TRUE)
ext_regional_H12<-vect("/Users/georgevagenas/Desktop/Vagenas_aSDMs/output/training_ext/ext_of_training_regional_H12.gpkg")


#example
unique(ext_regional_H12$Sp)
#mapview(ext_regional_H12[14])

# ==============================================================================
# 2. GLOBAL DATASET (thres_global_gbif)
# ==============================================================================

# Global Ecoregion Only
#ext_global_eco <- generate_training_extent(thres_global_gbif, ecoregions, hydro_layer = NULL)
terra::crs(ext_global_eco)<-"EPSG:4326"
#writeVector(ext_global_eco, file.path(out_dir, "ext_of_training_global_ecoregion.gpkg"), overwrite = TRUE)
ext_global_eco<-vect("/Users/georgevagenas/Desktop/Vagenas_aSDMs/output/training_ext/ext_of_training_global_ecoregion.gpkg")


#example
unique(ext_global_eco$Sp)
#mapview(ext_global_eco[43])

# Global H5
#ext_global_H5 <- generate_training_extent(thres_global_gbif, ecoregions, hydrosheds_H5)
terra::crs(ext_global_H5)<-"EPSG:4326"
#writeVector(ext_global_H5, file.path(out_dir, "ext_of_training_global_H5.gpkg"), overwrite = TRUE)
ext_global_H5<-vect("/Users/georgevagenas/Desktop/Vagenas_aSDMs/output/training_ext/ext_of_training_global_H5.gpkg")


#example#
#mapview(ext_global_H5[43])

# Global H8
#ext_global_H8 <- generate_training_extent(thres_global_gbif, ecoregions, hydrosheds_H8)
terra::crs(ext_global_H8)<-"EPSG:4326"
#writeVector(ext_global_H8, file.path(out_dir, "ext_of_training_global_H8.gpkg"), overwrite = TRUE)
ext_global_H8<-vect("/Users/georgevagenas/Desktop/Vagenas_aSDMs/output/training_ext/ext_of_training_global_H8.gpkg")


#example
#mapview(ext_global_H8[43])

# Global H12
#ext_global_H12 <- generate_training_extent(thres_global_gbif, ecoregions, hydrosheds_H12)
terra::crs(ext_global_H12)<-"EPSG:4326"
#writeVector(ext_global_H12, file.path(out_dir, "ext_of_training_global_H12.gpkg"), overwrite = TRUE)
ext_global_H12<-vect("/Users/georgevagenas/Desktop/Vagenas_aSDMs/output/training_ext/ext_of_training_global_H12.gpkg")

#example
#mapview(ext_global_H12[43])

#cat("✅ All 8 extent files generated successfully!\n")0



#Initialize objects for the models

# 1. Extract the attributes
thres_iberian_vect
thres_global_gbif
str(thres_iberian_df)
str(thres_global_df)

#prepare the structures
# Initialize an empty dataframe to store the results
results_df <- data.frame(
  id=integer(),
  species_name= integer(),
  e_AUC = numeric(),
  e_COR = numeric(),
  t_maxSSS = numeric(),
  t_maxkappa = numeric(),
  t_prevalence = numeric(),
  CBI=numeric(),
  uAUC=numeric(), #Calculate its importance
  maxKappa=numeric(),
  maxTSS=numeric(),
  obs_prevalence=numeric(),
  stringsAsFactors = FALSE
)

#create an object
raster_list<-list()

# Initialize an empty SpatRaster object
combined_rasters <- rast()

#minimum_background_points=5

sdm_d<-list()
vi_list<-list() #the list for the variable importance
raster_list<-list() #the list for the rasters



# ==============================================================================
# 0. SETUP & LIBRARIES (STRICTLY SEQUENTIAL)
# ==============================================================================

try(parallel::stopCluster(cl), silent = TRUE) 
closeAllConnections() 

library(terra)
library(sdm)
library(dplyr)
library(duckdb)
library(arrow)
#install.packages("vandalico")
library(vandalico)

# ==============================================================================
# 1. SETTINGS & DATA LOADING
# ==============================================================================
start_time <- Sys.time()

base_dir <- "/Users/georgevagenas/Desktop/Vagenas_aSDMs/output/"

global_dir <- file.path(base_dir, "global")
regional_dir <- file.path(base_dir, "regional")
if(!dir.exists(global_dir)) dir.create(global_dir)
if(!dir.exists(regional_dir)) dir.create(regional_dir)


#Parquet generation

## 1. Define the ABSOLUTE path where you want the Parquet file saved
## (Adjust this folder path if you want it saved somewhere else)
# parquet_file <- "/Users/georgevagenas/Desktop/Vagenas_aSDMs/output/aSDMs_Full_Parquet.parquet"
# 
# # 2. Build the optimized Parquet database from your global raster stack
# cat("[SYSTEM] Building Optimized Database...\n")
# 
# env_df <- as.data.frame(predictors_global, cells=TRUE, xy=TRUE, na.rm=TRUE)
# 
# # Clean up column names (lowercase and replace dots with underscores)
# names(env_df) <- tolower(names(env_df))
# if(!"cell" %in% names(env_df)) env_df$cell <- 1:nrow(env_df)
# safe_names <- gsub("\\.", "_", names(env_df)) 
# names(env_df) <- safe_names
# 
# # Write it to the hard drive
# write_parquet(env_df, parquet_file)
# 
# cat("[SYSTEM] Parquet database built successfully at:\n", parquet_file, "\n")

parquet_file <- "/Users/georgevagenas/Desktop/Vagenas_aSDMs/output/aSDMs_Full_Parquet.parquet"

algo_list <- c('glmp', 'brt', 'maxent', 'rf', 'gam')
n_reps <- 5

study_area_no_islands <- vect("/Users/georgevagenas/Desktop/Vagenas_aSDMs/input/SpatialExtents/study_area_iberia.shp")

# ==============================================================================
# CORRECTED PREDICTOR SETS (Hardcoded to match Parquet exactly)
# ==============================================================================

sets_global <- list(
  "Climate"       = c("bio5_clima", "bio16_clima", "bio17_clima", "bio15_clima", "bio4_clima"),
  "Hydroclimatic" = c("bio4_hydro", "bio1_hydro", "bio16_hydro", "bio17_hydro", "bio15_hydro")
)

sets_regional <- list(
  "Climate"            = c("bio5_clima", "bio16_clima", "bio17_clima", "bio15_clima", "bio4_clima"),
  "Hydroclimatic"      = c("bio4_hydro", "bio1_hydro", "bio16_hydro", "bio17_hydro", "bio15_hydro"),
  "Hydromorphological" = c("lka_pc_use", "dor_pc_pva", "sgr_dk_rav", "urb_pc_use", "for_pc_use")
)

extents_global <- list(eco = ext_global_eco, H5 = ext_global_H5, H8 = ext_global_H8, H12 = ext_global_H12)
extents_regional <- list(eco = ext_regional_eco, H5 = ext_regional_H5, H8 = ext_regional_H8, H12 = ext_regional_H12)

global_species_list <- unique(thres_global_df$Sp)
regional_species_list <- unique(thres_iberian_df$Sp)

ref_raster <- predictors_regional[[1]]
iberia_all_cells <- cells(ref_raster, study_area_no_islands, touches=TRUE)[, "cell"]

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

# Database Fetcher (Pulls X and Y directly from Parquet!)
get_env_data <- function(cell_ids, vars) {
  con <- dbConnect(duckdb::duckdb())
  on.exit(dbDisconnect(con, shutdown = TRUE), add = TRUE) # Auto-closes connection
  
  cols_sql <- paste(paste0("e.", vars), collapse = ", ")
  ids_sql  <- paste(cell_ids, collapse=",")
  
  # ADDED e.x and e.y here!
  q_vars <- sprintf("SELECT e.cell, e.x, e.y, %s FROM read_parquet('%s') e WHERE e.cell IN (%s)", cols_sql, parquet_file, ids_sql)
  
  env_data <- dbGetQuery(con, q_vars)
  names(env_data) <- tolower(names(env_data))
  return(env_data)
}

extract_custom_metrics <- function(obs, pred, sp_name, model_id) {
  valid <- !is.na(pred) & !is.na(obs)
  obs <- obs[valid]; pred <- pred[valid]
  
  if(length(unique(obs)) < 2 || length(obs) < 10) {
    return(data.frame(id=model_id, species_name=sp_name, e_AUC=NA, e_COR=NA, t_maxSSS=NA, 
                      t_maxkappa=NA, t_prevalence=NA, CBI=NA, maxKappa=NA, maxTSS=NA, 
                      obs_prevalence=NA, uAUC=NA, stringsAsFactors=FALSE))
  }
  
  e <- evaluates(x=obs, p=pred)
  b <- try(sdm:::.boyce(obs, pred), silent=TRUE)
  
  e_AUC <- if(length(e@statistics$AUC)>0) as.numeric(e@statistics$AUC[1]) else NA
  e_COR <- if(length(e@statistics$COR)>0) as.numeric(e@statistics$COR[1]) else NA
  CBI <- if(!inherits(b, "try-error") && !is.null(b$CBI)) as.numeric(b$CBI[1]) else NA
  
  tb <- e@threshold_based
  maxTSS_idx <- which.max(tb$TSS)
  maxKappa_idx <- which.max(tb$Kappa)
  
  maxTSS <- tb$TSS[maxTSS_idx]
  maxKappa <- tb$Kappa[maxKappa_idx]
  
  t_maxSSS <- tb$threshold[which.max(tb$sensitivity + tb$specificity)]
  t_maxkappa <- tb$threshold[maxKappa_idx]
  
  obs_prevalence <- sum(obs == 1) / length(obs)
  t_prevalence <- tb$threshold[which.min(abs(tb$threshold - obs_prevalence))]
  
  # Calculate uniform AUC using vandalico package
  # 1. Force pure numeric vectors
  clean_obs <- as.numeric(obs)
  clean_pred <- as.numeric(pred)
  
  # 2. CRITICAL FIX: Bind with PREDICTIONS in col 1, OBSERVATIONS in col 2
  mat_data <- na.omit(cbind(clean_pred, clean_obs))
  
  # 3. Safety Check: We need BOTH 1s and 0s in the observation column (Column 2)
  if (length(unique(mat_data[, 2])) < 2) {
    uAUC_val <- NA
  } else {
    # 4. Run vandalico on pristine data
    uAUC_res <- tryCatch({
      vandalico::AUCuniform(mat = mat_data, by = 0.1, rep = 100)
    }, error = function(e) {
      NA
    })
    
    # 5. Extract safely
    uAUC_val <- if(is.list(uAUC_res) && "uAUC" %in% names(uAUC_res)) {
      as.numeric(uAUC_res$uAUC[1])
    } else if (is.numeric(uAUC_res)) {
      as.numeric(uAUC_res[1])
    } else {
      NA
    }
  }
  
  data.frame(
    id = model_id, species_name = sp_name, e_AUC = e_AUC, e_COR = e_COR,
    t_maxSSS = t_maxSSS, t_maxkappa = t_maxkappa, t_prevalence = t_prevalence,
    CBI = CBI, maxKappa = maxKappa, maxTSS = maxTSS, obs_prevalence = obs_prevalence,
    uAUC = uAUC_val, stringsAsFactors = FALSE
  )
}



# ==============================================================================
# 🎯 TARGETED DEBUG RUN: Syngnathus abaster only 
# ==============================================================================
# target_sp <- "Syngnathus abaster"
# target_ext <- c("H5","H8","H12")
# 
# # 1. Restrict the Species Lists
# # (intersect() ensures we don't accidentally add it to a list it doesn't belong to)
# global_species_list <- intersect(global_species_list, target_sp)
# regional_species_list <- intersect(regional_species_list, target_sp)
# 
# #2. Restrict the Extent Lists (Supports multiple extents)
# # Keep only the target extents that actually exist in your master lists
# valid_global_exts <- intersect(names(extents_global), target_ext)
# extents_global <- extents_global[valid_global_exts]
# 
# valid_regional_exts <- intersect(names(extents_regional), target_ext)
# extents_regional <- extents_regional[valid_regional_exts]


# # ==============================================================================
# # 🎯 RESUME RUN: Start from the species AFTER "Alosa fallax"
# # ==============================================================================
# last_completed_sp <- "Alosa fallax"
# 
# # 1. Resume Global Species List
# if (last_completed_sp %in% global_species_list) {
#   idx_g <- which(global_species_list == last_completed_sp)
#   
#   # Check if it wasn't the very last species in the list
#   if (idx_g < length(global_species_list)) {
#     global_species_list <- global_species_list[(idx_g + 1):length(global_species_list)]
#   } else {
#     global_species_list <- character(0) # Empty the list if it was the last one
#   }
# }


#Advanced resetting


# ==============================================================================
# 🎯 OUTSIDE SETUP: Set Exact Starting Point
# ==============================================================================
target_sp <- "Barbatula barbatula"

# 1. Chop the Species List
idx_g <- which(global_species_list == target_sp)
global_species_list <- global_species_list[idx_g:length(global_species_list)]

# 2. Set the initial Extents and Sets to start precisely at H5 / Hydroclimatic
run_extents <- names(extents_global)[which(names(extents_global) == "H5"):length(extents_global)]
run_sets <- names(sets_global)[which(names(sets_global) == "Hydroclimatic"):length(sets_global)]
# ==============================================================================


# ==============================================================================
# PHASE 1: GLOBAL MODELING (51 Widespread Species)
# ==============================================================================
cat("\n====================================================================\n")
cat("                   STARTING PHASE 1: GLOBAL MODELS                    \n")
cat("====================================================================\n")

for (sp in global_species_list) {
  clean_sp <- gsub(" ", "_", sp)
  sp_df <- thres_global_df[thres_global_df$Sp == sp, ]
  
  cat(sprintf("\n>>> GLOBAL PROCESSING: %s\n", sp))
  
  # Notice: we use 'run_extents' here instead of names(extents_global)
  for (ext_name in run_extents) {
    cat(sprintf("  -> Extent: %s\n", ext_name))
    
    # --- BULLETPROOF EXTENT EXTRACTION ---
    current_ext_obj <- extents_global[[ext_name]]
    if (inherits(current_ext_obj, "list")) {
      sp_ext_poly <- current_ext_obj[[sp]]
    } else {
      sp_ext_poly <- current_ext_obj[current_ext_obj$Sp == sp, ]
    }
    
    # !!! NEW FIX: FORCE CRS ALIGNMENT !!!
    if(crs(sp_ext_poly) != crs(ref_raster)) {
      sp_ext_poly <- terra::project(sp_ext_poly, crs(ref_raster))
    }
    
    # 1. STRICT SPATIAL MASKING
    extracted_data <- terra::extract(ref_raster, sp_ext_poly, cells=TRUE, touches=TRUE)
    all_valid_ids <- extracted_data$cell[!is.na(extracted_data[, 2])]
    
    # 2. Background sampling strictly within the extent
    presence_ids <- sp_df$cell_id
    n_bg_target <- max(5, round(length(all_valid_ids) * 0.05)) 
    
    set.seed(123)
    bg_ids <- sample(all_valid_ids, n_bg_target, replace = (length(all_valid_ids) < n_bg_target))
    master_cell_ids <- c(presence_ids, bg_ids)
    
    # 3. Build coordinate dataframe
    df_presence <- data.frame(cell = presence_ids, presence = 1)
    df_background <- data.frame(cell = bg_ids, presence = 0)
    base_train_df <- rbind(df_presence, df_background)
    
    # Notice: we use 'run_sets' here instead of names(sets_global)
    for (set_name in run_sets) {
      cat(sprintf("    -> Set: %s\n", set_name))
      
      out_path <- file.path(global_dir, clean_sp, ext_name, set_name)
      if(!dir.exists(out_path)) dir.create(out_path, recursive = TRUE)
      
      t_vars <- sets_global[[set_name]]
      
      # 1. PURE SPATIAL QUERY
      sp_ext_coords <- as.vector(ext(sp_ext_poly))
      xmin <- sp_ext_coords['xmin']; xmax <- sp_ext_coords['xmax']
      ymin <- sp_ext_coords['ymin']; ymax <- sp_ext_coords['ymax']
      
      con <- dbConnect(duckdb::duckdb())
      cols_sql <- paste(paste0("e.", t_vars), collapse = ", ")
      q_box <- sprintf("SELECT e.cell, e.x, e.y, %s FROM read_parquet('%s') e WHERE e.x BETWEEN %f AND %f AND e.y BETWEEN %f AND %f",
                       cols_sql, parquet_file, xmin, xmax, ymin, ymax)
      box_data <- dbGetQuery(con, q_box)
      dbDisconnect(con)
      names(box_data) <- tolower(names(box_data))
      
      # 2. ULTRA-FAST POLYGON MASKING
      box_pts <- vect(box_data, geom=c("x", "y"), crs=crs(sp_ext_poly), keepgeom=TRUE)
      inside_pts <- box_pts[sp_ext_poly, ] 
      poly_map_data <- as.data.frame(inside_pts)
      
      # 3. Identify true land cells and extract presences
      true_land_cells <- poly_map_data$cell
      presence_ids <- sp_df$cell_id
      valid_pres_ids <- intersect(presence_ids, true_land_cells)
      
      # 4. Sample background STRICTLY from true land cells
      n_bg_target <- max(5, round(length(true_land_cells) * 0.05)) 
      set.seed(123)
      bg_ids <- sample(true_land_cells, n_bg_target, replace = (length(true_land_cells) < n_bg_target))
      
      # 5. Build train_ready directly
      df_pres <- poly_map_data[poly_map_data$cell %in% valid_pres_ids, ]
      if(nrow(df_pres) > 0) df_pres$presence <- 1
      df_bg <- poly_map_data[poly_map_data$cell %in% bg_ids, ]
      df_bg$presence <- 0
      train_ready <- rbind(df_pres, df_bg)
      
      # 6. sdmData and Model Training
      d <- sdmData(as.formula(paste0("presence ~ ", paste(t_vars, collapse="+"), "+coords(x+y)")), 
                   train=train_ready[train_ready$presence==1,], bg=train_ready[train_ready$presence==0,])
      
      m_global <- tryCatch({
        sdm(presence ~ ., d, methods = algo_list, replication = 'boot', n = n_reps,
            parallelSetting=list(ncore=10, method="parallel"))
      }, error = function(e) {
        cat(sprintf("      [WARNING] sdm() failed internally: %s. Skipping.\n", e$message))
        return(NULL)
      })
      
      if(is.null(m_global) || !inherits(m_global,"sdmModels")) {
        cat("      [WARNING] Model invalid or empty. Skipping to next.\n")
        next
      }
      
      saveRDS(m_global, file.path(out_path, "global_SDM.rds"))      
      
      # [PREDICT, CROP, AND EVALUATE]
      global_raster <- NULL 
      
      p_ens_global <- tryCatch({
        ensemble(m_global, poly_map_data[, c("x", "y", t_vars)], setting=list(method='weighted', stat='AUC'))
      }, error = function(e) {
        cat("      [WARNING] ensemble() failed. Triggering pure dataframe fallback...\n")
        return(NULL)
      })
      
      if (!is.null(p_ens_global)) {
        res_df_glob <- data.frame(x=poly_map_data$x, y=poly_map_data$y, val=as.numeric(p_ens_global[[1]]))
        global_raster <- rast(res_df_glob, type="xyz", crs=crs(sp_ext_poly))
      } else {
        ev <- tryCatch(getEvaluation(m_global, stat="AUC"), error=function(e) NULL)
        if (!is.null(ev)) {
          good_models <- ev[!is.na(ev$AUC) & ev$AUC > 0.5, ]
          if (nrow(good_models) > 0) {
            raw_preds <- tryCatch({
              predict(m_global, newdata = poly_map_data[, t_vars], id = good_models$modelID)
            }, error = function(e) {
              cat(sprintf("      [WARNING] predict() failed: %s\n", e$message))
              return(NULL)
            })
            if (!is.null(raw_preds)) {
              weights <- good_models$AUC / sum(good_models$AUC)
              raw_matrix <- as.matrix(raw_preds)
              weighted_matrix <- sweep(raw_matrix, 2, weights, `*`)
              final_vals <- rowSums(weighted_matrix, na.rm = TRUE)
              res_df_glob <- data.frame(x = poly_map_data$x, y = poly_map_data$y, val = final_vals)
              global_raster <- rast(res_df_glob, type="xyz", crs=crs(sp_ext_poly))
            }
          }
        }
      }
      
      writeRaster(global_raster, file.path(out_path, paste0("ensemble_global_full.tif")), overwrite = TRUE)
      
      iberia_cropped <- crop(global_raster, study_area_no_islands)
      iberia_masked <- mask(iberia_cropped, study_area_no_islands)
      writeRaster(iberia_masked, file.path(out_path, paste0("ensemble_global_iberia_cropped.tif")), overwrite = TRUE)
      
      # ========================================================================
      # STRICT IBERIAN EVALUATION 
      # ========================================================================
      ib_pres_ids <- thres_iberian_df[thres_iberian_df$Sp == sp, "cell_id"]
      set.seed(123)
      ib_bg_ids <- sample(iberia_all_cells, max(15, round(length(iberia_all_cells)*0.05)), replace=FALSE)
      eval_obs <- c(rep(1, length(ib_pres_ids)), rep(0, length(ib_bg_ids)))
      eval_xy <- xyFromCell(ref_raster, c(ib_pres_ids, ib_bg_ids))
      eval_preds <- terra::extract(iberia_masked, eval_xy)
      pred_vals <- if(is.data.frame(eval_preds)) eval_preds[, ncol(eval_preds)] else as.numeric(eval_preds)
      
      valid_eval <- !is.na(pred_vals) & !is.na(eval_obs)
      metrics_df <- extract_custom_metrics(eval_obs[valid_eval], pred_vals[valid_eval], sp, paste0(ext_name, "_", set_name, "_global"))
      write.csv(metrics_df, file.path(out_path, "eval_global_iberia_strict.csv"), row.names=FALSE)
      
      # ========================================================================
      # BULLETPROOF VARIABLE IMPORTANCE
      # ========================================================================
      vi <- tryCatch({
        getVarImp(m_global, id="ensemble")
      }, error = function(e) {
        cat("      [WARNING] getVarImp 'ensemble' failed. Retrying with valid model IDs...\n")
        ev <- tryCatch(getEvaluation(m_global, stat="AUC"), error=function(err) NULL)
        if (!is.null(ev)) {
          good_ids <- ev$modelID[!is.na(ev$AUC) & ev$AUC > 0.5]
          if (length(good_ids) > 0) {
            tryCatch({ getVarImp(m_global, id=good_ids) }, error = function(err2) NULL)
          } else { return(NULL) }
        } else { return(NULL) }
      })
      
      if (!is.null(vi)) {
        saveRDS(vi, file.path(out_path, "varImp_global.rds"))
      } else {
        cat("      [WARNING] Variable Importance completely failed. Skipping RDS save.\n")
      }
      
      # ========================================================================
      # CLEANUP & CRASH PREVENTION (Runs at the end of every Set)
      # ========================================================================
      rm(m_global, d, global_raster, iberia_cropped, iberia_masked, p_ens_global, box_pts, inside_pts)
      gc(verbose=FALSE)
      terra::tmpFiles(current=FALSE, orphan=TRUE, remove=TRUE) # Prevents hard drive crash
      
    } # <-- END OF SET LOOP
    
    # 🎯 HEAL THE SETS LIST: Once H5 Hydroclimatic finishes, reset the list to full size!
    run_sets <- names(sets_global)
    
  } # <-- END OF EXTENT LOOP
  
  # 🎯 HEAL THE EXTENTS LIST: Once Barbatula is completely done, reset to full size!
  run_extents <- names(extents_global)
  
} # <-- END OF SPECIES LOOP


# ==============================================================================
# PHASE 2: REGIONAL MODELING (98 Species)
# ==============================================================================
cat("\n====================================================================\n")
cat("                  STARTING PHASE 2: REGIONAL MODELS                   \n")
cat("====================================================================\n")

# [NEW] Ensure correct [WIDESPREAD] vs [ENDEMIC] tags
full_global_species <- unique(thres_global_df$Sp)

for (sp in regional_species_list) {
  clean_sp <- gsub(" ", "_", sp)
  sp_df <- thres_iberian_df[thres_iberian_df$Sp == sp, ]
  
  # [OLD as silenced]
  # is_widespread <- (sp %in% global_species_list)
  # [NEW]
  is_widespread <- (sp %in% full_global_species)
  status_tag <- ifelse(is_widespread, "[WIDESPREAD - 5+1 Preds]", "[ENDEMIC - 5 Preds]")
  
  cat(sprintf("\n>>> REGIONAL PROCESSING: %s %s\n", sp, status_tag))
  
  for (ext_name in names(extents_regional)) {
    cat(sprintf("  -> Extent: %s\n", ext_name))
    
    # --- BULLETPROOF EXTENT EXTRACTION ---
    current_ext_obj <- extents_regional[[ext_name]]
    if (inherits(current_ext_obj, "list")) {
      sp_ext_poly <- current_ext_obj[[sp]]
    } else {
      sp_ext_poly <- current_ext_obj[current_ext_obj$Sp == sp, ]
    }
    
    if(is.null(sp_ext_poly) || nrow(sp_ext_poly) == 0) {
      cat("      [WARNING] Polygon not found for this species. Skipping.\n")
      next 
    }
    
    # [NEW] FORCE CRS ALIGNMENT
    if(crs(sp_ext_poly) != crs(ref_raster)) {
      sp_ext_poly <- terra::project(sp_ext_poly, crs(ref_raster))
    }
    
    for (set_name in names(sets_regional)) {
      cat(sprintf("    -> Set: %s\n", set_name))
      
      out_path <- file.path(regional_dir, clean_sp, ext_name, set_name)
      if(!dir.exists(out_path)) dir.create(out_path, recursive = TRUE)
      
      t_vars <- sets_regional[[set_name]]
      
      # 1. PURE SPATIAL QUERY
      sp_ext_coords <- as.vector(ext(sp_ext_poly))
      xmin <- sp_ext_coords['xmin']; xmax <- sp_ext_coords['xmax']
      ymin <- sp_ext_coords['ymin']; ymax <- sp_ext_coords['ymax']
      
      con <- dbConnect(duckdb::duckdb())
      cols_sql <- paste(paste0("e.", t_vars), collapse = ", ")
      q_box <- sprintf("SELECT e.cell, e.x, e.y, %s FROM read_parquet('%s') e WHERE e.x BETWEEN %f AND %f AND e.y BETWEEN %f AND %f",
                       cols_sql, parquet_file, xmin, xmax, ymin, ymax)
      box_data <- dbGetQuery(con, q_box)
      dbDisconnect(con)
      
      if(nrow(box_data) == 0) {
        cat("      [WARNING] No data found in bounding box. Skipping.\n")
        next
      }
      
      names(box_data) <- tolower(names(box_data))
      
      # 2. ULTRA-FAST POLYGON MASKING
      box_pts <- vect(box_data, geom=c("x", "y"), crs=crs(sp_ext_poly), keepgeom=TRUE)
      inside_pts <- box_pts[sp_ext_poly, ] 
      poly_map_data <- as.data.frame(inside_pts)
      
      if(nrow(poly_map_data) == 0) {
        cat("      [WARNING] No land cells inside regional polygon. Skipping.\n")
        next
      }
      
      # 3. DYNAMIC GLOBAL COVARIATE INJECTION 
      if (is_widespread) {
        global_set_match <- ifelse(set_name == "Climate", "Climate", "Hydroclimatic")
        global_tif_path <- file.path(global_dir, clean_sp, ext_name, global_set_match, "ensemble_global_iberia_cropped.tif")
        
        if(file.exists(global_tif_path)) {
          global_cov_raster <- rast(global_tif_path)
          cov_vals <- terra::extract(global_cov_raster, cbind(poly_map_data$x, poly_map_data$y))
          if(is.data.frame(cov_vals)) cov_vals <- cov_vals[, ncol(cov_vals)]
          cov_vals[is.na(cov_vals)] <- 0 
          
          poly_map_data$global_suitability <- as.numeric(cov_vals)
          t_vars <- c(t_vars, "global_suitability")
        } else {
          cat("      [WARNING] Global covariate TIF not found. Reverting to 5 predictors.\n")
        }
      }
      
      # 4. SPATIAL PRESENCE EXTRACTION 
      iberia_template <- crop(ref_raster, study_area_no_islands)
      pres_xy <- as.data.frame(xyFromCell(iberia_template, sp_df$cell_id))
      
      r_data <- poly_map_data[, c("x", "y", "cell", t_vars)]
      if ("global_suitability" %in% names(poly_map_data)) {
        r_data$global_suitability <- poly_map_data$global_suitability
      }
      poly_raster <- rast(r_data, type="xyz", crs=crs(sp_ext_poly))
      pres_env <- terra::extract(poly_raster, as.matrix(pres_xy))
      
      df_pres <- data.frame(cell = pres_env$cell, x = pres_xy$x, y = pres_xy$y)
      for (v in t_vars) { df_pres[[v]] <- pres_env[[v]] }
      if ("global_suitability" %in% names(poly_map_data)) {
        df_pres$global_suitability <- pres_env$global_suitability
      }
      
      df_pres$presence <- 1
      df_pres <- na.omit(df_pres) 
      
      if(nrow(df_pres) < 5) {
        cat(sprintf("      [WARNING] Only %d valid presences in this regional map. Skipping.\n", nrow(df_pres)))
        next
      }
      
      # 5. SAMPLE BACKGROUND
      true_land_cells <- poly_map_data$cell
      eval_bg_pool <- setdiff(true_land_cells, df_pres$cell) 
      
      n_bg_target <- max(15, round(length(true_land_cells) * 0.05)) 
      set.seed(123)
      bg_ids <- sample(eval_bg_pool, n_bg_target, replace = (length(eval_bg_pool) < n_bg_target))
      
      # 6. BUILD TRAIN_READY
      df_bg <- poly_map_data[poly_map_data$cell %in% bg_ids, ]
      df_bg$presence <- 0
      df_bg <- df_bg[, names(df_pres)] 
      train_ready <- rbind(df_pres, df_bg)
      
      # 7. MODEL TRAINING
      d <- sdmData(as.formula(paste0("presence ~ ", paste(t_vars, collapse="+"), "+coords(x+y)")), 
                   train=train_ready[train_ready$presence==1,], bg=train_ready[train_ready$presence==0,])
      
      m_regional <- tryCatch({
        sdm(presence ~ ., d, methods = algo_list, replication = 'boot', n = n_reps,
            parallelSetting=list(ncore=10, method="parallel"))
      }, error = function(e) {
        cat(sprintf("      [WARNING] sdm() failed internally: %s. Skipping.\n", e$message))
        return(NULL)
      })
      
      if(is.null(m_regional) || !inherits(m_regional,"sdmModels")) {
        cat("      [WARNING] Model invalid or empty. Skipping to next.\n")
        next
      }
      saveRDS(m_regional, file.path(out_path, "regional_SDM.rds"))      
      
      # [OLD as silenced]
      # p_ens <- ensemble(m_regional, poly_map_data[, c("x", "y", t_vars)], setting=list(method='weighted', stat='AUC'))
      # res_df <- data.frame(x=poly_map_data$x, y=poly_map_data$y, val=as.numeric(p_ens[[1]]))
      # en_raster <- rast(res_df, type="xyz", crs=crs(sp_ext_poly))
      # writeRaster(en_raster, file.path(out_path, paste0("ensemble_regional.tif")), overwrite = TRUE)
      
      # [NEW] BULLETPROOF REGIONAL PREDICTION
      en_raster <- NULL 
      p_ens <- tryCatch({
        ensemble(m_regional, poly_map_data[, c("x", "y", t_vars)], setting=list(method='weighted', stat='AUC'))
      }, error = function(e) {
        cat("      [WARNING] ensemble() failed. Triggering pure dataframe fallback...\n")
        return(NULL)
      })
      
      if (!is.null(p_ens)) {
        res_df <- data.frame(x=poly_map_data$x, y=poly_map_data$y, val=as.numeric(p_ens[[1]]))
        en_raster <- rast(res_df, type="xyz", crs=crs(sp_ext_poly))
      } else {
        ev <- tryCatch(getEvaluation(m_regional, stat="AUC"), error=function(e) NULL)
        if (!is.null(ev)) {
          good_models <- ev[!is.na(ev$AUC) & ev$AUC > 0.5, ]
          if (nrow(good_models) > 0) {
            raw_preds <- tryCatch({ predict(m_regional, newdata = poly_map_data[, t_vars], id = good_models$modelID) }, error=function(e) NULL)
            if (!is.null(raw_preds)) {
              weights <- good_models$AUC / sum(good_models$AUC)
              weighted_matrix <- sweep(as.matrix(raw_preds), 2, weights, `*`)
              en_raster <- rast(data.frame(x = poly_map_data$x, y = poly_map_data$y, val = rowSums(weighted_matrix, na.rm=TRUE)), type="xyz", crs=crs(sp_ext_poly))
            }
          }
        }
      }
      
      if (!is.null(en_raster)) {
        writeRaster(en_raster, file.path(out_path, paste0("ensemble_regional.tif")), overwrite = TRUE)
      }
      
      # [OLD as silenced]
      # p_train_eval <- ensemble(m_regional, train_ready[, c("x", "y", t_vars)], setting=list(method='weighted', stat='AUC'))
      # pred_vals <- as.numeric(p_train_eval[[1]])
      
      # [NEW] BULLETPROOF EVALUATION
      p_train_eval <- tryCatch({
        ensemble(m_regional, train_ready[, c("x", "y", t_vars)], setting=list(method='weighted', stat='AUC'))
      }, error = function(e) NULL)
      
      if (!is.null(p_train_eval)) {
        pred_vals <- as.numeric(p_train_eval[[1]])
        valid_eval <- !is.na(pred_vals) & !is.na(train_ready$presence)
        metrics_df <- extract_custom_metrics(
          train_ready$presence[valid_eval], pred_vals[valid_eval], sp, paste0(ext_name, "_", set_name, "_regional")
        )
        write.csv(metrics_df, file.path(out_path, "eval_regional_iberia_strict.csv"), row.names=FALSE)
      }
      
      # [OLD as silenced]
      # vi <- getVarImp(m_regional, id="ensemble")
      # saveRDS(vi, file.path(out_path, "varImp_regional.rds"))
      
      # [NEW] BULLETPROOF VARIABLE IMPORTANCE
      vi <- tryCatch({
        getVarImp(m_regional, id="ensemble")
      }, error = function(e) {
        cat("      [WARNING] getVarImp 'ensemble' failed. Retrying...\n")
        ev <- tryCatch(getEvaluation(m_regional, stat="AUC"), error=function(err) NULL)
        if (!is.null(ev)) {
          good_ids <- ev$modelID[!is.na(ev$AUC) & ev$AUC > 0.5]
          if (length(good_ids) > 0) {
            tryCatch({ getVarImp(m_regional, id=good_ids) }, error=function(err2) NULL)
          } else { return(NULL) }
        } else { return(NULL) }
      })
      
      if (!is.null(vi)) saveRDS(vi, file.path(out_path, "varImp_regional.rds"))
      
      # [OLD as silenced]
      # rm(m_regional, d, poly_map_data, en_raster, p_ens, box_pts, inside_pts); gc(verbose=FALSE)
      
      # [NEW] Garbage Collection & Hard Drive Crash Prevention
      rm(m_regional, d, poly_map_data, en_raster, p_ens, p_train_eval, box_pts, inside_pts)
      gc(verbose=FALSE)
      terra::tmpFiles(current=FALSE, orphan=TRUE, remove=TRUE)
    }
  }
}


end_time <- Sys.time()
cat(sprintf("\nAnalysis Done. Total execution time: %s\n", round(end_time - start_time, 2)))



