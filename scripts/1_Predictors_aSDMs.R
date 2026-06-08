#June 2026
#Development pipeline of freshwater SDMs (aSDMS)

#Sector: Predictors_aSDMs

#PhD Researcher - Georgios Vagenas (georgios.vagenas@mncn.csic.es | georgvagenas@gmail.com)

#### Pre-setting :: Libraries required to perform the analysis ####


library(raster)
library(terra)
library(dplyr)
library(tidyr)
library(ggplot2)
library(sdm)
library(sf)

#Let's start

#### SECTION_1_Global hydroclimatic layers CMPI5 (ECMWF-FutureStreams) & Hydromorphology ####

setwd("/Users/georgevagenas/Desktop/Vagenas_aSDMs/input/")

#hydro
hydro_global<-rast("Predictors/E2O_hydro_historical.tiff")
plot(hydro_global[[1]])
hydro_global

#clima
clima_global<-rast("Predictors/ECMWF_bioclim_historical_5min_downscaled.tiff")
plot(clima_global[[1]])
clima_global

hydro_global


#Align the rasters
  library(terra)

# 1. Subset by name (safest method)
bio_nums <- c(1, 4, 5, 6, 12, 15, 16, 17)
clima_sub <- clima_global[[paste0("BIO", bio_nums, "_clima")]]
hydro_sub <- hydro_global[[paste0("BIO", bio_nums, "_hydro")]]

# 2. FORCE alignment (Harmless: Changes metadata, NOT pixel values)
set.ext(hydro_sub, ext(clima_sub))
set.crs(hydro_sub, crs(clima_sub))

# 3. Combine
combined_stack <- c(clima_sub, hydro_sub)

# 4. Remove NAs where EITHER layer is missing (The Common Mask)
# This ensures every pixel has all 16 values
common_mask <- sum(combined_stack) 
final_stack <- mask(combined_stack, common_mask)

final_stack #verified correct allignmet

#writeRaster(final_stack,"Predictors/aligned_hydroclimatic_global.tiff",overwrite=FALSE)


final_stack<-rast("Predictors/aligned_hydroclimatic_global.tiff")
#hydromorpho predictors

library(sf)
library(dplyr)



library(terra)

# 1. Get the list of shapefiles (Ensure you use the dot in full.names!)
shape_folder <- "/Users/georgevagenas/Desktop/Vagenas_aSDMs/input/Predictors/RiverATLAS_Data_v10_shp/RiverATLAS_v10_shp/"
files <- list.files(path = shape_folder, pattern = "\\.shp$", full.names = TRUE)

# 2. Read all files into a list of SpatVector objects
# vect() is the terra function for reading vector data

# vect_list <- lapply(files, function(x) {
#   message("Reading: ", basename(x))
#   vect(x)
# })

# 3. Merge them all into a single global SpatVector
# do.call applies terra's fast rbind() to the whole list

#global_vect <- do.call(rbind, vect_list)

# # 4. Verify the merged object
# print(global_vect)
# plot(global_vect) # Warning: Plotting global RiverATLAS might take a moment!

# 5. Save as a single GeoPackage

# writeVector(global_vect, 
#             filename = "Predictors/RiverATLAS_Global.gpkg", 
#             filetype = "GPKG", 
#             overwrite = TRUE)



library(terra)

# 1. Define your target columns (Now including the lake parameter)
my_vars <- c("HYRIV_ID", 
             "slp_dg_uav", 
             "sgr_dk_rav", 
             "dor_pc_pva", 
             "pac_pc_use", 
             "urb_pc_use", 
             "for_pc_use", 
             "rdd_mk_uav",
             "lka_pc_use") # <-- Added here

# 2. Load the massive global vector
global_vect <- vect("Predictors/RiverATLAS_Global.gpkg")

# 3. Extract only the 8 columns that we want so that clima, hydro and hydromorpho have the same amount of candidate predictors
hdrmrph_vector <- global_vect[, my_vars]

# # 5. Save this as your final, lightweight Predictor file
# writeVector(hdrmrph_vector, 
#             filename = "Predictors/RiverATLAS_hdrmrph_vector.gpkg", 
#             filetype = "GPKG",
#             overwrite = TRUE)


hdrmrph_vector <- vect("Predictors/RiverATLAS_hdrmrph_vector.gpkg")

plot(hdrmrph_vector)

final_stack

#align both hydroclimatic and hydromorphological layers

library(terra)

# 1. Define the columns to rasterize (excluding the ID)
vars_to_rasterize <- c("slp_dg_uav", "sgr_dk_rav", "dor_pc_pva", 
                       "pac_pc_use", "urb_pc_use", "for_pc_use", 
                       "rdd_mk_uav", "lka_pc_use")

# 2. Loop through each variable one by one
# This avoids the "length > 1" error by feeding them individually
# raster_list <- lapply(vars_to_rasterize, function(var) {
#   message("Rasterizing: ", var) # This will print progress in the console!
#   rasterize(x = hdrmrph_vector, 
#             y = final_stack, 
#             field = var, 
#             fun = "mean")
# })

# 3. Stack the resulting individual rasters back into one multi-layer object
hdrmrph_raster <- do.call(c, raster_list)

# Ensure the names match perfectly
names(hdrmrph_raster) <- vars_to_rasterize

# 4. Verify output
print(hdrmrph_raster)

#writeRaster(hdrmrph_raster,"Predictors/RiverATLAS_hdrmrph_raster.tiff")

hdrmrph_raster<-rast("Predictors/RiverATLAS_hdrmrph_raster.tiff")


#remove those grids with negative values (NAs) & brief report on the removals

library(terra)

# 1. FIND the grids that are below zero
# This creates a mathematical map where 1 = "Below Zero" and 0 = "Valid"
negative_grids <- (hdrmrph_raster[["slp_dg_uav"]] < 0) | (hdrmrph_raster[["sgr_dk_rav"]] < 0)

# 2. Calculate the exact number of those negative grids
# sum() adds up all the 1s (the grids that are below zero)
grids_to_remove <- global(negative_grids, fun = "sum", na.rm = TRUE)$sum

# 3. Get the initial total of river grids before we do anything
initial_grids <- global(hdrmrph_raster[["slp_dg_uav"]], fun = "notNA")$notNA

# 4. REMOVE them from all layers
# ifel logic: If negative_grids is TRUE (1), set to NA. Otherwise, keep the original hdrmrph_raster values.
hdrmrph_clean <- ifel(negative_grids, NA, hdrmrph_raster)

# 5. Get the final count
final_grids <- global(hdrmrph_clean[["slp_dg_uav"]], fun = "notNA")$notNA

# 6. Print the Short Report
percent_removed <- round((grids_to_remove / initial_grids) * 100, 2)

cat("================ MASKING REPORT ================\n")
cat("Initial valid river grids : ", format(initial_grids, big.mark=","), "\n")
cat("Negative grids found      : ", format(grids_to_remove, big.mark=","), "\n")
cat("Final valid river grids   : ", format(final_grids, big.mark=","), "\n")
cat("Percentage removed        : ", percent_removed, "%\n")
cat("================================================\n")


#3.42% removal of negative values


library(terra)

# 1. Create the Binary Mask (1 = Keep, 0 = Remove)
# If both are >= 0, it gets a 1. If either is negative, it gets a 0.
binary_mask <- ifel(hdrmrph_raster[["slp_dg_uav"]] >= 0 & 
                      hdrmrph_raster[["sgr_dk_rav"]] >= 0, 
                    1, 0)

names(binary_mask) <- "Keep_Remove_Mask"

# 2. Check the Binary Mask
# You can look at the summary and plot it to verify where the 0s are
print(binary_mask)
plot(binary_mask, col = c("red", "green"), main = "Mask: 0 = Remove (Red), 1 = Keep (Green)")

# 3. Apply the Mask to create the Final Clean Raster
# 'maskvalues = 0' tells R to look at the binary mask and turn every '0' into NA across all 8 layers
hdrmrph_clean <- mask(hdrmrph_raster, binary_mask, maskvalues = 0)

# 4. Verify the Final Clean Stack
# The minimum values for slp_dg_uav and sgr_dk_rav should now be 0 (no more negatives)
plot(hdrmrph_clean$sgr_dk_rav)



#FINAL STEP:: Align hydromorphology (hdrmrph_clean) with hydroclimatic layers (final_stack)

library(terra)

# 1. Mask River Morphology by Climate
# We use the first layer of final_stack as the spatial cookie-cutter
hdrmrph_final <- mask(hdrmrph_clean, final_stack[[1]])

# 2. Mask Climate by River Morphology (Crucial for SDMs!)
# This removes climate pixels in the middle of landmasses where no rivers exist
final_stack_aligned <- mask(final_stack, hdrmrph_final[[1]])

# 3. Combine them into the Ultimate 24-Layer Stack
# (16 Climate/Hydro layers + 8 River Morphology layers)
asdms_predictors_global <- c(final_stack_aligned, hdrmrph_final)

# 4. Verify the perfect alignment
# Check if the number of valid cells (notNA) is identical across all 24 layers
global(asdms_predictors_global, fun = "notNA")

# 5. Save the final predictor stack for your SDM!
# writeRaster(asdms_predictors_global, 
#             filename = "Predictors/predictors_clima_hydro_hydromorpho_global_raw/asdms_predictors_global.tif",
#             overwrite = FALSE)

cat("Success! Your global SDM predictors are perfectly aligned and saved.\n")

predictors_raw<-rast("Predictors/predictors_clima_hydro_hydromorpho_global_raw/asdms_predictors_global.tif")

###standardize them from 0 to 1### in global applications we don't care but now it will supress local variations

library(terra)

# 1. Explicitly calculate the Mean and SD from the RAW raster
# global() handles the math across all millions of pixels safely
layer_means <- global(predictors_raw, fun = "mean", na.rm = TRUE)
layer_sds <- global(predictors_raw, fun = "sd", na.rm = TRUE)

# 2. Build the data frame
# global() returns a data frame where the values are in the first column
scaling_params <- data.frame(
  Variable = names(predictors_raw),
  Mean = layer_means$mean,
  SD = layer_sds$sd
)

# 3. Save the parameters to a CSV file for your future projections
#write.csv(scaling_params, "Predictors/predictors_clima_hydro_hydromorpho_global_standardized/Scaling_Parameters.csv", row.names = FALSE)

# 4. Now perform the actual scaling
predictors_scaled <- scale(predictors_raw)

# 5. Check the dataframe to ensure it worked!
print(head(scaling_params))

# 6. Save the final standardized current predictors
# writeRaster(predictors_scaled, 
#             filename = "Predictors/predictors_clima_hydro_hydromorpho_global_standardized/asdms_predictors_global_scaled.tif", 
#             overwrite = FALSE)

###FINAL STANDARDIZED LAYER

predictors_scaled<-rast("Predictors/predictors_clima_hydro_hydromorpho_global_standardized/asdms_predictors_global_scaled.tif")

plot(predictors_scaled$for_pc_use)





#### SECTION_2_PCA SELECTION OF THE 5 MOST DOMINANT AT OF EACH ####


### GLOBAL SETTING ###

library(terra)
library(raster)
library(dplyr)

# 1. Load your 24-layer master stack
# (Assuming asdms_predictors_global.tif contains your 24 scaled variables)
setwd("/Users/georgevagenas/Desktop/Vagenas_aSDMs/input/")
predictors <- rast("Predictors/predictors_clima_hydro_hydromorpho_global_standardized/asdms_predictors_global_scaled.tif")

library(terra)

# 1. Fast Sample
# Using 20,000 points - this should take < 5 seconds
set.seed(123)
env_sample <- spatSample(predictors, size = 20000, method = "random", 
                         na.rm = TRUE, as.df = TRUE)

# Remove ID column if present
env_sample <- env_sample[, names(predictors)]


# --- 1. CLIMATE GROUP ---
pca_clim <- prcomp(env_sample[, 1:8], scale. = TRUE)
# Get cumulative variance
cum_clim <- summary(pca_clim)$importance[3,]
# Variables driving the 80% variance axes
axes_clim <- which(cum_clim >= 0.80)[1]
load_clim <- abs(pca_clim$rotation[, 1:axes_clim])
# If more than 1 axis, sum them; if not, just take the first
imp_clim <- if(axes_clim > 1) rowSums(load_clim) else load_clim
top_clim <- names(sort(imp_clim, decreasing = TRUE)[1:5])

# --- 2. HYDROLOGY GROUP ---
pca_hydro <- prcomp(env_sample[, 9:16], scale. = TRUE)
cum_hydro <- summary(pca_hydro)$importance[3,]
axes_hydro <- which(cum_hydro >= 0.80)[1]
load_hydro <- abs(pca_hydro$rotation[, 1:axes_hydro])
imp_hydro <- if(axes_hydro > 1) rowSums(load_hydro) else load_hydro
top_hydro <- names(sort(imp_hydro, decreasing = TRUE)[1:5])

# --- 3. HYDROMORPHOLOGY GROUP ---
pca_morph <- prcomp(env_sample[, 17:24], scale. = TRUE)
cum_morph <- summary(pca_morph)$importance[3,]
axes_morph <- which(cum_morph >= 0.80)[1]
load_morph <- abs(pca_morph$rotation[, 1:axes_morph])
imp_morph <- if(axes_morph > 1) rowSums(load_morph) else load_morph
top_morph <- names(sort(imp_morph, decreasing = TRUE)[1:5])

# Combined Candidate List
candidate_list <- unique(c(top_clim, top_hydro, top_morph))
print(candidate_list)


library(corrplot)

# Check Climate Group
corrplot(cor(env_sample[, top_clim]), method="number", title="Climate Top 5", mar=c(0,0,1,0))

# Check Hydro Group
corrplot(cor(env_sample[, top_hydro]), method="number", title="Hydro Top 5", mar=c(0,0,1,0))

# Check Morph Group
corrplot(cor(env_sample[, top_morph]), method="number", title="Morph Top 5", mar=c(0,0,1,0))


###lets check the correlations between them now and if there are replace with the next most important



# 1. Ensure you have the 'importance' rankings saved for each group
# (Using the imp_clim, imp_hydro, imp_morph from your previous step)
rank_clim  <- names(sort(imp_clim, decreasing = TRUE))
rank_hydro <- names(sort(imp_hydro, decreasing = TRUE))
rank_morph <- names(sort(imp_morph, decreasing = TRUE))

# Your current "Top 5"
top_5_clim  <- rank_clim[1:5]
top_5_hydro <- rank_hydro[1:5]
top_5_morph <- rank_morph[1:5]



library(usdm)
library(corrplot)


select_best_5_strict <- function(full_rankings, data_sample, group_name) {
  # Take the top 8 as the initial pool
  current_selection <- full_rankings[1:8]
  
  message(paste("\n--- Processing Group:", group_name, "---"))
  
  while (length(current_selection) > 5) {
    # 1. Calculate Correlation Matrix
    cor_matrix <- cor(data_sample[, current_selection])
    
    # 2. Find the highest absolute correlation pair (excluding the diagonal)
    diag(cor_matrix) <- 0
    max_cor_val <- max(abs(cor_matrix))
    
    # 3. Check: If highest correlation is > 0.7, drop the lower-ranked variable
    if (max_cor_val > 0.7) {
      # Find which variables are involved
      indices <- which(abs(cor_matrix) == max_cor_val, arr.ind = TRUE)
      var1 <- current_selection[indices[1, 1]]
      var2 <- current_selection[indices[1, 2]]
      
      # Drop the one that is further down the PCA rankings (higher index = lower rank)
      rank1 <- which(full_rankings == var1)
      rank2 <- which(full_rankings == var2)
      var_to_drop <- if(rank1 > rank2) var1 else var2
      
      message(paste("Found Pearson r =", round(max_cor_val, 2), 
                    "between", var1, "and", var2, ". Dropping:", var_to_drop))
      
    } else {
      # 4. If no Pearson > 0.7, use VIF to find the most redundant variable
      v_res <- vif(data_sample[, current_selection])
      var_to_drop <- v_res$Variables[which.max(v_res$VIF)]
      message(paste("No high Pearson found. Dropping highest VIF:", 
                    var_to_drop, "(VIF:", round(max(v_res$VIF), 2), ")"))
    }
    
    # Remove the variable and repeat
    current_selection <- current_selection[current_selection != var_to_drop]
  }
  
  return(current_selection)
}

# 2. Execute per group
# (Assuming rank_clim, rank_hydro, rank_morph are your PCA-sorted variables)
final_clim_5  <- select_best_5_forced(rank_clim, env_sample, "Climate")
final_hydro_5 <- select_best_5_forced(rank_hydro, env_sample, "Hydrology")
final_morph_5 <- select_best_5_forced(rank_morph, env_sample, "Hydromorphology")

# 3. Combine into final stack of 15
final_15_vars <- c(final_clim_5, final_hydro_5, final_morph_5)


# A. Final VIF Table


# Set up a plotting area with 3 columns to see them side-by-side
# Alternatively, run them one by one if you want larger individual plots
par(mfrow=c(1,3)) 

# 1. Climate Group Corrplot
corrplot(cor(env_sample[, final_clim_5]), 
         method = "color", 
         type = "lower", 
         addCoef.col = "black", 
         number.cex = 0.8, 
         tl.cex = 0.9, 
         tl.col = "black",
         diag = FALSE,
         title = "Climate (Top 5)",
         mar = c(0,0,2,0))

# 2. Hydrology Group Corrplot
corrplot(cor(env_sample[, final_hydro_5]), 
         method = "color", 
         type = "lower", 
         addCoef.col = "black", 
         number.cex = 0.8, 
         tl.cex = 0.9, 
         tl.col = "black",
         diag = FALSE,
         title = "Hydrology (Top 5)",
         mar = c(0,0,2,0))

# 3. Hydromorphology Group Corrplot
corrplot(cor(env_sample[, final_morph_5]), 
         method = "color", 
         type = "lower", 
         addCoef.col = "black", 
         number.cex = 0.8, 
         tl.cex = 0.9, 
         tl.col = "black",
         diag = FALSE,
         title = "Hydromorphology (Top 5)",
         mar = c(0,0,2,0))

# Reset plotting layout to single plot
par(mfrow=c(1,1))



#So actual differences 




# --- 1. SETUP DATA STRUCTURES ---
# The original PCA-based candidates (Top 5 from each group before pruning)
candidate_list

# Your final 15 variables (The ones you officially selected to keep)
# Note: In your case, these matched the candidates, but the 'Report' 
# will confirm their internal correlations.
final_15_vars


# 1. Define the sets
candidate_list

final_15_vars

# 2. Define category indices (assuming standard order of 5 per group)
groups <- list(
  Climate = list(cand = candidate_list[1:5], final = final_15_vars[1:5]),
  Hydrology = list(cand = candidate_list[6:10], final = final_15_vars[6:10]),
  Morphology = list(cand = candidate_list[11:15], final = final_15_vars[11:15])
)

# 3. Generate Categorical Report
cat("--- CATEGORICAL SELECTION REPORT ---\n\n")

for (g in names(groups)) {
  dropped <- setdiff(groups[[g]]$cand, groups[[g]]$final)
  added   <- setdiff(groups[[g]]$final, groups[[g]]$cand)
  
  cat(paste0("[", toupper(g), " GROUP]\n"))
  cat("  Dropped: ", ifelse(length(dropped) > 0, paste(dropped, collapse=", "), "None"), "\n")
  cat("  Added:   ", ifelse(length(added) > 0, paste(added, collapse=", "), "None"), "\n\n")
}

#save the global subset setting

predictors_global<-predictors[[final_15_vars]]
writeRaster(predictors_global,"/Users/georgevagenas/Desktop/Vagenas_aSDMs/input/Predictors/predictors_finalized_global/predictors_finalized_global.tiff")





### REGIONAL (IBERIAN ECOREGIONAL) SETTING ###

library(terra)
library(raster)
library(dplyr)

# 1. Load your 24-layer master stack
# (Assuming asdms_predictors_global.tif contains your 24 scaled variables)
setwd("/Users/georgevagenas/Desktop/Vagenas_aSDMs/input/")
predictors <- rast("Predictors/predictors_clima_hydro_hydromorpho_global_standardized/asdms_predictors_global_scaled.tif")

library(terra)


#crop the global with the iberian ecoregional setting

ecoregions<-vect("/Users/georgevagenas/Desktop/Vagenas_aSDMs/input/SpatialExtents/ecoregions/feow_hydrosheds.shp")

library(terra)

# 1. Create a bounding box for Mainland Iberia
# This avoids selecting distant island territories
iberia_extent <- ext(-10, 4.5, 35.5, 44) 
iberia_poly <- as.polygons(iberia_extent, crs="EPSG:4326")

# 2. Ensure your ecoregions have the same CRS 
# (FEOW is typically WGS84 / EPSG:4326)
if(crs(ecoregions) == "") {
  crs(ecoregions) <- "EPSG:4326"
}

# 3. Spatial Intersect
# This selects any ecoregion polygon that touches the Iberian mainland box
ecoregions_iberia <- ecoregions[iberia_poly, ]

# 4. Refine the Selection
# Depending on the FEOW precision, you might get a few "neighbor" regions.
# You can view the IDs to confirm they match the Iberian basins (e.g., 401, 402, 403...)
print(ecoregions_iberia)

# 5. Visual Check
plot(ecoregions_iberia, col="lightblue", border="darkblue", lwd=1.5)
plot(iberia_poly, border="red", add=TRUE) # Show the bounding box used
text(ecoregions_iberia, "FEOW_ID", cex=0.8)


#IBERIAN FRESHWATER ECOREGIONS: 403, 412, 413, 414



library(terra)

# 1. Subset by specific FEOW IDs
# 403: Western Iberia, 412: Cantabric Coast, 413: Southern Iberia, 414: Eastern Iberia
target_ids <- c(403, 412, 413, 414)
eco_subset <- ecoregions[ecoregions$FEOW_ID %in% target_ids, ]

# 2. Merge (Dissolve) into a single study area polygon
eco_merged <- aggregate(eco_subset)

plot(eco_merged)

# 3. Remove Islands
# We convert the merged multi-polygon into individual parts (disaggregate)
# then filter by area to keep only the large mainland masses.
eco_parts <- disagg(eco_merged)

# Calculate area in square meters (or based on your CRS)
eco_parts$area_calc <- expanse(eco_parts)

# Keep only the largest parts (the mainland)
# Usually, the mainland is several orders of magnitude larger than islands.
# We take the top 1 or 2 parts depending on if the polygons are perfectly connected.
mainland_threshold <- max(eco_parts$area_calc) * 0.1 # Anything less than 10% of the max is an island
eco_mainland <- eco_parts[eco_parts$area_calc > mainland_threshold, ]

# 4. Final Aggregation
# Merge the remaining large parts back into one clean study area
study_area_iberia <- aggregate(eco_mainland)

# 5. Visual Verification
plot(study_area_iberia, col="lightgreen", border="darkgreen", main="Mainland Iberian Ecoregions (403, 412, 413, 414)")


#save the eco_iberian_merged_polygon

#writeVector(study_area_iberia,"/Users/georgevagenas/Desktop/Vagenas_aSDMs/input/SpatialExtents/study_area_iberia.shp")

study_area_iberia<-vect("/Users/georgevagenas/Desktop/Vagenas_aSDMs/input/SpatialExtents/study_area_iberia.shp")


#####

#crop the global layer for the iberia

predictors_regional<-crop(predictors,study_area_iberia,mask=T)
plot(predictors_regional)
#####


# 1. Fast Sample
# Using 20,000 points - this should take < 5 seconds
set.seed(123)
env_sample_regional <- spatSample(predictors_regional, size = 10000, method = "random", 
                         na.rm = TRUE, as.df = TRUE)

# Remove ID column if present
env_sample <- env_sample_regional[, names(predictors_regional)]


# --- 1. CLIMATE GROUP ---
pca_clim <- prcomp(env_sample[, 1:8], scale. = TRUE)
# Get cumulative variance
cum_clim <- summary(pca_clim)$importance[3,]
# Variables driving the 80% variance axes
axes_clim <- which(cum_clim >= 0.80)[1]
load_clim <- abs(pca_clim$rotation[, 1:axes_clim])
# If more than 1 axis, sum them; if not, just take the first
imp_clim <- if(axes_clim > 1) rowSums(load_clim) else load_clim
top_clim <- names(sort(imp_clim, decreasing = TRUE)[1:5])

# --- 2. HYDROLOGY GROUP ---
pca_hydro <- prcomp(env_sample[, 9:16], scale. = TRUE)
cum_hydro <- summary(pca_hydro)$importance[3,]
axes_hydro <- which(cum_hydro >= 0.80)[1]
load_hydro <- abs(pca_hydro$rotation[, 1:axes_hydro])
imp_hydro <- if(axes_hydro > 1) rowSums(load_hydro) else load_hydro
top_hydro <- names(sort(imp_hydro, decreasing = TRUE)[1:5])

# --- 3. HYDROMORPHOLOGY GROUP ---
pca_morph <- prcomp(env_sample[, 17:24], scale. = TRUE)
cum_morph <- summary(pca_morph)$importance[3,]
axes_morph <- which(cum_morph >= 0.80)[1]
load_morph <- abs(pca_morph$rotation[, 1:axes_morph])
imp_morph <- if(axes_morph > 1) rowSums(load_morph) else load_morph
top_morph <- names(sort(imp_morph, decreasing = TRUE)[1:5])

# Combined Candidate List
candidate_list_regional <- unique(c(top_clim, top_hydro, top_morph))
print(candidate_list_regional)


library(corrplot)

# Check Climate Group
corrplot(cor(env_sample[, top_clim]), method="number", title="Climate Top 5", mar=c(0,0,1,0))

# Check Hydro Group
corrplot(cor(env_sample[, top_hydro]), method="number", title="Hydro Top 5", mar=c(0,0,1,0))

# Check Morph Group
corrplot(cor(env_sample[, top_morph]), method="number", title="Morph Top 5", mar=c(0,0,1,0))


###lets check the correlations between them now and if there are replace with the next most important



# 1. Ensure you have the 'importance' rankings saved for each group
# (Using the imp_clim, imp_hydro, imp_morph from your previous step)
rank_clim  <- names(sort(imp_clim, decreasing = TRUE))
rank_hydro <- names(sort(imp_hydro, decreasing = TRUE))
rank_morph <- names(sort(imp_morph, decreasing = TRUE))

# Your current "Top 5"
top_5_clim  <- rank_clim[1:5]
top_5_hydro <- rank_hydro[1:5]
top_5_morph <- rank_morph[1:5]



library(usdm)
library(corrplot)


select_best_5_strict <- function(full_rankings, data_sample, group_name) {
  # Take the top 8 as the initial pool
  current_selection <- full_rankings[1:8]
  
  message(paste("\n--- Processing Group:", group_name, "---"))
  
  while (length(current_selection) > 5) {
    # 1. Calculate Correlation Matrix
    cor_matrix <- cor(data_sample[, current_selection])
    
    # 2. Find the highest absolute correlation pair (excluding the diagonal)
    diag(cor_matrix) <- 0
    max_cor_val <- max(abs(cor_matrix))
    
    # 3. Check: If highest correlation is > 0.7, drop the lower-ranked variable
    if (max_cor_val > 0.7) {
      # Find which variables are involved
      indices <- which(abs(cor_matrix) == max_cor_val, arr.ind = TRUE)
      var1 <- current_selection[indices[1, 1]]
      var2 <- current_selection[indices[1, 2]]
      
      # Drop the one that is further down the PCA rankings (higher index = lower rank)
      rank1 <- which(full_rankings == var1)
      rank2 <- which(full_rankings == var2)
      var_to_drop <- if(rank1 > rank2) var1 else var2
      
      message(paste("Found Pearson r =", round(max_cor_val, 2), 
                    "between", var1, "and", var2, ". Dropping:", var_to_drop))
      
    } else {
      # 4. If no Pearson > 0.7, use VIF to find the most redundant variable
      v_res <- vif(data_sample[, current_selection])
      var_to_drop <- v_res$Variables[which.max(v_res$VIF)]
      message(paste("No high Pearson found. Dropping highest VIF:", 
                    var_to_drop, "(VIF:", round(max(v_res$VIF), 2), ")"))
    }
    
    # Remove the variable and repeat
    current_selection <- current_selection[current_selection != var_to_drop]
  }
  
  return(current_selection)
}

# 2. Execute per group
# (Assuming rank_clim, rank_hydro, rank_morph are your PCA-sorted variables)
final_clim_5  <- select_best_5_forced(rank_clim, env_sample, "Climate")
final_hydro_5 <- select_best_5_forced(rank_hydro, env_sample, "Hydrology")
final_morph_5 <- select_best_5_forced(rank_morph, env_sample, "Hydromorphology")

# 3. Combine into final stack of 15
final_15_vars_regional <- c(final_clim_5, final_hydro_5, final_morph_5)


# A. Final VIF Table


# Set up a plotting area with 3 columns to see them side-by-side
# Alternatively, run them one by one if you want larger individual plots
par(mfrow=c(1,3)) 

# 1. Climate Group Corrplot
corrplot(cor(env_sample[, final_clim_5]), 
         method = "color", 
         type = "lower", 
         addCoef.col = "black", 
         number.cex = 0.8, 
         tl.cex = 0.9, 
         tl.col = "black",
         diag = FALSE,
         title = "Climate (Top 5)",
         mar = c(0,0,2,0))

# 2. Hydrology Group Corrplot
corrplot(cor(env_sample[, final_hydro_5]), 
         method = "color", 
         type = "lower", 
         addCoef.col = "black", 
         number.cex = 0.8, 
         tl.cex = 0.9, 
         tl.col = "black",
         diag = FALSE,
         title = "Hydrology (Top 5)",
         mar = c(0,0,2,0))

# 3. Hydromorphology Group Corrplot
corrplot(cor(env_sample[, final_morph_5]), 
         method = "color", 
         type = "lower", 
         addCoef.col = "black", 
         number.cex = 0.8, 
         tl.cex = 0.9, 
         tl.col = "black",
         diag = FALSE,
         title = "Hydromorphology (Top 5)",
         mar = c(0,0,2,0))

# Reset plotting layout to single plot
par(mfrow=c(1,1))






#So actual differences 




# --- 1. SETUP DATA STRUCTURES ---
# The original PCA-based candidates (Top 5 from each group before pruning)
candidate_list_regional

# Your final 15 variables (The ones you officially selected to keep)
# Note: In your case, these matched the candidates, but the 'Report' 
# will confirm their internal correlations.
final_15_vars_regional


# 1. Define the sets
candidate_list_regional

final_15_vars_regional

# 2. Define category indices (assuming standard order of 5 per group)
groups <- list(
  Climate = list(cand = candidate_list_regional[1:5], final = final_15_vars_regional[1:5]),
  Hydrology = list(cand = candidate_list_regional[6:10], final = final_15_vars_regional[6:10]),
  Morphology = list(cand = candidate_list_regional[11:15], final = final_15_vars_regional[11:15])
)

# 3. Generate Categorical Report
cat("--- CATEGORICAL SELECTION REPORT ---\n\n")

for (g in names(groups)) {
  dropped <- setdiff(groups[[g]]$cand, groups[[g]]$final)
  added   <- setdiff(groups[[g]]$final, groups[[g]]$cand)
  
  cat(paste0("[", toupper(g), " GROUP]\n"))
  cat("  Dropped: ", ifelse(length(dropped) > 0, paste(dropped, collapse=", "), "None"), "\n")
  cat("  Added:   ", ifelse(length(added) > 0, paste(added, collapse=", "), "None"), "\n\n")
}

#save the regional subset setting

predictors_regional<-predictors_regional[[final_15_vars_regional]]
writeRaster(predictors_regional,"/Users/georgevagenas/Desktop/Vagenas_aSDMs/input/Predictors/predictors_finalized_regional/predictors_finalized_regional.tiff")

predictors_regional<-rast("/Users/georgevagenas/Desktop/Vagenas_aSDMs/input/Predictors/predictors_finalized_regional/predictors_finalized_regional.tiff")


plot(predictors_regional)
plot(predictors_global)
