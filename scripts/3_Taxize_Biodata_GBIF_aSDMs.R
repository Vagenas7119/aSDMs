#June 2026
#Development pipeline of freshwater SDMs (aSDMS)

#Sector: Taxonomic update and sorting of Iberian dataset ~ Global GBIF [Final data preparation]

#PhD Researcher - Georgios Vagenas (georgios.vagenas@mncn.csic.es | georgvagenas@gmail.com)

#### Pre-setting :: Libraries required to perform the analysis ####


library(terra)
library(dplyr)
library(taxize)
library(terra)
library(tidyr)


global_dataset<-read.csv("/Users/georgevagenas/Desktop/Vagenas_aSDMs/input/Dataset/GBIF_final_strict_data_107_thinned_cleaned_373975/GBIF_final_strict_data_107_thinned_cleaned_373975.csv")
str(global_dataset)
regional_dataset<-vect("/Users/georgevagenas/Desktop/Vagenas_aSDMs/input/Dataset/vect_MITECO_SNIPAD_IBERIA/vect_dataset_IBERIA.shp")
regional_dataset

#First step is to compare the two datasets


# ==============================================================================
# 1. FORMAT GLOBAL DATASET TO MATCH REGIONAL SPATVECTOR
# ==============================================================================
message("Formatting the global dataset...")

# Filter out empty coordinates/species and create the matching columns
global_formatted <- global_dataset %>%
  # Ensure valid coordinates and species names
  filter(!is.na(decimalLongitude) & !is.na(decimalLatitude)) %>%
  filter(!is.na(species) & species != "") %>%
  
  # Create the exact columns found in regional_dataset
  mutate(
    Sp = species,
    Source_DB = "GBIF",      
    Category = "Global Data", 
    presence = 1             
  ) %>%
  
  # Keep only the matching attributes plus coordinates for spatial conversion
  dplyr::select(Sp, Source_DB, Category, presence, decimalLongitude, decimalLatitude)

# Convert to a SpatVector using WGS 84 (EPSG:4326) to match regional_dataset
global_vect <- vect(global_formatted, 
                    geom = c("decimalLongitude", "decimalLatitude"), 
                    crs = "EPSG:4326")

message("Global dataset successfully converted to SpatVector!")

# ==============================================================================
# 2. EXTRACT SPECIES LISTS
# ==============================================================================
# Get the unique species from both datasets
regional_sp <- unique(regional_dataset$Sp)
global_sp <- unique(global_vect$Sp)

# Combine them into one master list for the taxonomic check
all_unique_sp <- unique(c(regional_sp, global_sp))

cat("\n--- RAW DATASET COUNTS ---\n")
cat("Species in Regional (MITECO):", length(regional_sp), "\n")
cat("Species in Global (GBIF):    ", length(global_sp), "\n")
cat("Total unique names to check: ", length(all_unique_sp), "\n")

library(rgbif)
library(dplyr)

library(rgbif)

# ==============================================================================
# 3. TAXONOMIC CHECK WITH RGBIF (USING STANDARD LOOP)
# ==============================================================================
message("\nRunning taxonomic harmonization directly against the GBIF Backbone...")
message("Checking names one-by-one. This might take a minute depending on list size...")

# Create an empty vector to store our accepted names
accepted_names <- character(length(all_unique_sp))

# Loop through each unique name and check it against GBIF
for (i in seq_along(all_unique_sp)) {
  # Query GBIF for the single name
  gbif_result <- name_backbone(name = all_unique_sp[i])
  
  # Check if GBIF successfully matched it to an accepted species level name
  if ("species" %in% colnames(gbif_result) && !is.na(gbif_result$species)) {
    accepted_names[i] <- gbif_result$species
  } else {
    accepted_names[i] <- NA # Mark as NA if GBIF couldn't figure it out
  }
  
  # Print a little progress tracker every 10 species so you know it's not frozen!
  if (i %% 10 == 0) cat("Checked", i, "of", length(all_unique_sp), "names...\n")
}

message("Taxonomic check complete! Creating translation dictionary...")

# Create the translation dictionary (Raw Name -> Accepted GBIF Name)
name_dict <- setNames(accepted_names, all_unique_sp)

# Translate the regional list
clean_regional <- unname(name_dict[regional_sp])
# If GBIF couldn't find a match, keep the original name just in case
clean_regional <- ifelse(is.na(clean_regional), regional_sp, clean_regional)
clean_regional <- unique(clean_regional)

# Translate the global list
clean_global <- unname(name_dict[global_sp])
clean_global <- ifelse(is.na(clean_global), global_sp, clean_global)
clean_global <- unique(clean_global)



# ==============================================================================
# 4. OVERLAP AND COMPARISON
# ==============================================================================
# Calculate the intersections and differences using the harmonized names
overlap_sp <- intersect(clean_regional, clean_global)
only_regional_sp <- setdiff(clean_regional, clean_global)
only_global_sp <- setdiff(clean_global, clean_regional)

cat("\n=======================================================\n")
cat("               HARMONIZED OVERLAP RESULTS                \n")
cat("=======================================================\n")
cat("Harmonized Regional Species:     ", length(clean_regional), "\n")
cat("Harmonized Global Species:       ", length(clean_global), "\n")
cat("Overlap (Present in BOTH):       ", length(overlap_sp), "\n")
cat("Iberian Species missing in GBIF:  ", length(only_regional_sp), "\n")
cat("=======================================================\n\n")

# Print the missing ones so you can investigate them!
if(length(only_regional_sp) > 0) {
  cat("List of Regional species NOT found in the Global dataset:\n")
  print(only_regional_sp)
}



#which species were merged?

library(dplyr)

message("Hunting down the synonymized species...")

# 1. Create a data frame comparing the original MITECO names to their GBIF translations
synonym_check_df <- data.frame(
  Original_MITECO_Name = regional_sp,
  GBIF_Accepted_Name = unname(name_dict[regional_sp])
)

# 2. Apply the same NA fallback rule we used earlier
synonym_check_df$GBIF_Accepted_Name <- ifelse(
  is.na(synonym_check_df$GBIF_Accepted_Name), 
  synonym_check_df$Original_MITECO_Name, 
  synonym_check_df$GBIF_Accepted_Name
)

# 3. Find which GBIF accepted names appear more than once!
collapsed_species <- synonym_check_df %>%
  group_by(GBIF_Accepted_Name) %>%
  filter(n() > 1) %>%
  arrange(GBIF_Accepted_Name) %>%
  as.data.frame() # Print as a clean base dataframe

cat("\nHere are the Iberian species that GBIF collapsed into synonyms:\n")
print(collapsed_species)

# 
# > print(collapsed_species)
# Original_MITECO_Name    GBIF_Accepted_Name
# 1         Cobitis haasi      Cobitis paludica
# 2      Cobitis paludica      Cobitis paludica
# 3 Iberocypris palaciosi Iberocypris palaciosi
# 4    Squalius palaciosi Iberocypris palaciosi
#

###Save the updated taxonomically corrected regional dataset###

library(terra)

message("Applying harmonized taxonomy to the spatial dataset...")

# 1. Extract the original species column from your spatial dataset
current_species <- regional_dataset$Sp

# 2. Translate the names using the dictionary we built earlier
updated_species <- unname(name_dict[current_species])

# 3. Apply the fallback rule: if GBIF didn't have a match (NA), keep the original name
updated_species <- ifelse(is.na(updated_species), current_species, updated_species)

# 4. Replace the old column in the spatial dataset with the new, cleaned names
regional_dataset$Sp <- updated_species

# 5. Check to make sure it worked (the unique count should now be 120!)
cat("Unique species in updated dataset:", length(unique(regional_dataset$Sp)), "\n")

# ==============================================================================
# EXPORT THE NEW DATASET
# ==============================================================================
message("Saving the harmonized shapefile...")

# Define the new file path
# I added "_tax_check" to the filename as requested
output_path <- "C:/Users/geo_v/Desktop/MITECO_BIOCAMBIO/Dataset/vect_MITECO_SNIPAD_IBERIA_120taxonomically_correct/vect_dataset_IBERIA_tax_check.shp"

#fix geometries in the input file#

library(sf)
library(dplyr)

message("Cleaning the 20 corrupted 'ghost' points before exporting...")

# 1. Convert the modified terra object to an sf object
sf_regional <- st_as_sf(regional_dataset)

# 2. Force the species column to be plain text (prevents another type of crash)
sf_regional$Sp <- as.character(unlist(sf_regional$Sp))

# 3. DROP the empty geometries that are crashing GDAL!
sf_clean <- sf_regional %>%
  filter(!st_is_empty(.))

cat("Original rows: ", nrow(sf_regional), "\n")
cat("Cleaned rows ready to save: ", nrow(sf_clean), "\n")

# 4. Define the output path
output_path <- "C:/Users/geo_v/Desktop/MITECO_BIOCAMBIO/Dataset/vect_MITECO_SNIPAD_IBERIA_120taxonomically_corrected/tax_check_vect_dataset_IBERIA.shp"

message("Exporting safely...")

# 5. Save using st_write (which is much safer than writeVector)
#Locked
#st_write(
  obj = sf_clean, 
  dsn = output_path, 
  delete_dsn = TRUE, # Safely overwrites if the file already exists
  quiet = FALSE
)

message("Success! Your cleaned and harmonized dataset is saved.")


###save also the global


library(sf)
library(dplyr)

message("Preparing the global dataset for safe export...")

# 1. Convert the terra SpatVector to an sf object
sf_global <- st_as_sf(global_vect)

# 2. Force the text columns to be plain text to prevent GDAL panics
sf_global$Sp <- as.character(unlist(sf_global$Sp))
sf_global$Source_DB <- as.character(unlist(sf_global$Source_DB))
sf_global$Category <- as.character(unlist(sf_global$Category))

# 3. Drop any accidentally empty geometries (just to be completely safe!)
sf_clean_global <- sf_global %>%
  filter(!st_is_empty(.))

cat("Global rows ready to save: ", nrow(sf_clean_global), "\n")

# 4. Define the output path 
# Saving it in your MITECO folder next to the regional one!
global_output_path <- "C:/Users/geo_v/Desktop/MITECO_BIOCAMBIO/Dataset/vect_GLOBAL_GBIF_107taxonomically_corrected/vect_GLOBAL_GBIF_107taxonomically_corrected.shp"

message("Exporting safely via sf...")

# 5. Save using st_write
#Locked
#st_write(
  obj = sf_clean_global, 
  dsn = global_output_path, 
  delete_dsn = TRUE, # Safely overwrites if the file already exists
  quiet = FALSE
)

message(paste("Success! Saved to:", global_output_path))


####taxonomically corrected regional and global dataset####








####MERGE GRID OCCURENCES OF GBIF INTO THE REGIONAL TAXONOMICALLY CORRECTED DATASET####


#add the layers

library(terra)

regional_dataset<-vect("/Users/georgevagenas/Desktop/Vagenas_aSDMs/input/Dataset/vect_MITECO_SNIPAD_IBERIA_120taxonomically_corrected/tax_check_vect_dataset_IBERIA.shp")
regional_dataset
global_dataset<-vect("/Users/georgevagenas/Desktop/Vagenas_aSDMs/input/Dataset/vect_GLOBAL_GBIF_107taxonomically_corrected/vect_GLOBAL_GBIF_107taxonomically_corrected.shp")
global_dataset


#enrich

library(terra)
library(dplyr)
library(rnaturalearth)
library(sf)

message("Step 1: Preparing the 1:10m Mainland Mask...")

# Fetch high-resolution boundaries
world <- ne_countries(scale = 10, country = c("Spain", "Portugal"), returnclass = "sv")

# Dissolve and disaggregate to separate islands from mainland
peninsula_raw <- aggregate(world)
polys <- disagg(peninsula_raw)

# Keep only the largest polygon (the Mainland Peninsula)
mainland_mask <- polys[which.max(expanse(polys)), ]

# Final crop to ensure no stray distant territories are included
mainland_mask <- crop(mainland_mask, ext(-10.0, 4.5, 35.8, 44.0))

plot(mainland_mask, col = "lightgray", main = "Mainland Mask (Islands Removed)")


message("Step 2: Filtering species and masking to Mainland...")

# 1. Identify common species
common_species <- intersect(unique(regional_dataset$Sp), unique(global_dataset$Sp))

# 2. Pre-filter by species (reduces rows before spatial check)
reg_filtered  <- regional_dataset[regional_dataset$Sp %in% common_species, ]
glob_filtered <- global_dataset[global_dataset$Sp %in% common_species, ]

# 3. FAST SPATIAL FILTER: Use 'relate' or 'extract' instead of 'intersect'
# We check which points are "covered by" the mainland mask
message("Checking spatial overlap (this should be much faster)...")

# For Regional
is_inside_reg <- is.related(reg_filtered, mainland_mask, "coveredby")
reg_iberia    <- reg_filtered[is_inside_reg, ]

# For Global
is_inside_glob <- is.related(glob_filtered, mainland_mask, "coveredby")
glob_iberia    <- glob_filtered[is_inside_glob, ]

cat("MITECO points on mainland:", nrow(reg_iberia), "\n")
cat("GBIF points on mainland:  ", nrow(glob_iberia), "\n")



### ENDEMIC vs NON - ENDEMICS ###

endemicity<-read.csv("/Users/georgevagenas/Desktop/Vagenas_aSDMs/input/Dataset/Endemicity/Endemicity_Iberian.csv",sep=",")
str(endemicity)


library(dplyr)

# 1. Define the endemic tags to exclude
endemic_tags <- c("Endemic_Spain", "Endemic_Spain_France", "Endemic_Portugal", "Endemic_Iberian")

# 2. Get the names of the 120 species you actually have in your regional spatial data
spatial_species_120 <- unique(regional_dataset$Sp)

# 3. Filter the endemicity table to find which of those 120 are 'Widespread'
# We use Final_name to match your regional Sp column
widespread_metadata <- endemicity %>%
  filter(Final_name %in% spatial_species_120) %>%
  filter(!(Detailed_Status %in% endemic_tags))

widespread_list <- unique(widespread_metadata$Final_name)

# 4. Identify the Endemics among your 120 (to keep them MITECO-only)
endemic_list <- setdiff(spatial_species_120, widespread_list)

cat("Out of your 120 species:\n")
cat("- Widespread (Enrich with GBIF):", length(widespread_list), "\n")
cat("- Endemic (Keep MITECO only):   ", length(endemic_list), "\n")


# > cat("Out of your 120 species:\n")
# Out of your 120 species:
#   > cat("- Widespread (Enrich with GBIF):", length(widespread_list), "\n")
# - Widespread (Enrich with GBIF): 61 
# > cat("- Endemic (Keep MITECO only):   ", length(endemic_list), "\n")
# - Endemic (Keep MITECO only):    59 



#Assess only 61 species, the widespread#

# 1. Clean the widespread list of any hidden spaces or case issues
widespread_clean <- trimws(widespread_list)

# 2. Re-run the loop with an EXPLICIT firewall
enriched_list <- list()
all_sp_in_data <- unique(regional_dataset$Sp)

message("Applying strict firewall for 61 widespread species...")

library(rnaturalearth)
# Extract just Spain and Portugal
iberian_raw_sf <- world_sf %>%
  filter(admin %in% c("Spain", "Portugal"))


# Merge the two countries into a single continuous polygon
iberian_merged_sf <- st_union(iberian_raw_sf)

# Break the merged MULTIPOLYGON into individual, disconnected POLYGONs
iberia_parts <- st_cast(iberian_merged_sf, "POLYGON")

# Calculate the area of all parts and keep ONLY the single largest one (the mainland)
iberia_mainland_sf <- iberia_parts[which.max(st_area(iberia_parts)), ]

# Convert the clean mainland shapefile into a terra SpatVector
iberian_shp <- vect(iberia_mainland_sf)

r_template<-crop(predictors_regional$BIO6_clima,iberian_shp,mask=T)
r_template<-r_template*0
plot(r_template)

for (s in all_sp_in_data) {
  
  # A. REGIONAL BASE: Always included, always thinned
  s_reg <- regional_dataset[regional_dataset$Sp == s, ]
  s_reg <- s_reg[is.related(s_reg, mainland_mask, "coveredby"), ]
  s_reg$cell_id <- cellFromXY(r_template, crds(s_reg))
  s_reg_thinned <- as.data.frame(s_reg) %>% distinct(cell_id, .keep_all = TRUE)
  
  # B. THE FIREWALL: Only search GBIF if the species is in your 61 list
  if (trimws(s) %in% widespread_clean) {
    
    s_glob <- global_dataset[global_dataset$Sp == s, ]
    s_glob <- s_glob[is.related(s_glob, mainland_mask, "coveredby"), ]
    s_glob$cell_id <- cellFromXY(r_template, crds(s_glob))
    
    # Gap-filling: Only keep cells that MITECO/SNIPAD do not have
    s_glob_new_grids <- as.data.frame(s_glob) %>%
      filter(!(cell_id %in% s_reg_thinned$cell_id)) %>%
      distinct(cell_id, .keep_all = TRUE)
    
    # Store combined data for widespread species
    enriched_list[[s]] <- rbind(s_reg_thinned, s_glob_new_grids)
    
  } else {
    # PROTECT ENDEMICS: No GBIF search, no GBIF binding.
    # This species (like hispanica) stays purely regional.
    enriched_list[[s]] <- s_reg_thinned
  }
}

# 3. Overwrite your final vector with the corrected data
final_df_fixed <- do.call(rbind, enriched_list)
final_coords <- xyFromCell(r_template, final_df_fixed$cell_id)
enriched_iberia_vect <- vect(cbind(final_coords, final_df_fixed), 
                             geom = c("x", "y"), crs = "EPSG:4326")

cat("\n--- ENRICHMENT FIXED: 61 SPECIES ONLY ---\n")


#Check if indeed 61 were enriched#


# 1. Create the comparison table again from the FIXED vector
comparison_table <- as.data.frame(enriched_iberia_vect) %>%
  group_by(Sp) %>%
  summarise(
    Regional = sum(Source_DB %in% c("MITECO", "SNIPAD")),
    GBIF     = sum(Source_DB == "GBIF")
  ) %>%
  mutate(Gain = GBIF > 0)

# 2. The Final Totals
cat("Species with enrichment:", sum(comparison_table$Gain), "\n")

# 3. Check for 'Imposters' again
imposters_check <- comparison_table$Sp[comparison_table$Gain & !(comparison_table$Sp %in% widespread_clean)]
cat("Endemics accidentally enriched:", length(imposters_check), "\n")





#Assess additions




# 1. Create a raster of MITECO presence from the FINAL enriched dataset
# (This identifies cells where MITECO 'won' or was the only data)
final_miteco_points <- enriched_iberia_vect[enriched_iberia_vect$Source_DB %in% c("MITECO", "SNIPAD")]
r_final_miteco <- rasterize(final_miteco_points, r_template, field = "presence", fun = "max")

# 2. Create a raster of GBIF presence from the FINAL enriched dataset
# (This identifies the 'Gaps' that were filled)
final_gbif_points <- enriched_iberia_vect[enriched_iberia_vect$Source_DB == "GBIF"]
r_final_gbif <- rasterize(final_gbif_points, r_template, field = "presence", fun = "max")

# 3. Build the Source Map
# 1 = MITECO Only, 2 = GBIF Addition, 3 = Overlap
# NOTE: In your new loop, 'Overlap' cells are actually MITECO cells (GBIF was discarded).
# So we will map 1 = MITECO (Original/Shared) and 2 = GBIF (Pure Gain).
source_map <- rast(r_template, vals = 0)

has_miteco_final <- !is.na(r_final_miteco)
has_gbif_final   <- !is.na(r_final_gbif)

source_map[has_miteco_final] <- 1
source_map[has_gbif_final]   <- 2

# Clean up and Label
source_map <- mask(source_map, mainland_mask)
source_map[source_map == 0] <- NA
source_map <- as.factor(source_map)
levels(source_map) <- data.frame(ID=1:2, label=c("MITECO Base", "GBIF Gap-Fill"))
# 3. Define colors
richness_cols <- colorRampPalette(c("yellow", "#993300"))(50)
source_cols <- c("#FFD700", "#993300", "#4A90E2")


# Setup Layout (Top: Before/After, Bottom: Source Map)
layout(matrix(c(1, 2, 3, 3), nrow = 2, byrow = TRUE))


# Calculate Richness from the survivors
rich_before <- rasterize(regional_dataset[is.related(regional_dataset, mainland_mask, "coveredby")], 
                         r_template, field = "Sp", fun = "count")
rich_after  <- rasterize(enriched_iberia_vect, r_template, field = "Sp", fun = "count")


# --- PANEL 1 & 2: Richness ---
# Note: 'pax' and 'plg' ensure the legend fits in the margin
par(mar = c(3, 3, 3, 5))

plot(rich_before, 
     main = "MITECO & SNIPAD Richness", 
     col = richness_cols, 
     type = "continuous",
     plg = list(title = "Species", title.cex = 0.8))
plot(mainland_mask, add = TRUE, border = "black", lwd = 0.4)

plot(rich_after, 
     main = "GBIF Enriched Richness", 
     col = richness_cols,
     type = "continuous",
     plg = list(title = "Species", title.cex = 0.8))
plot(mainland_mask, add = TRUE, border = "black", lwd = 0.4)

# --- PANEL 3: Data Source Contribution ---
# For the categorical map, we use 'type=classes' to force the discrete legend
par(mar = c(4, 3, 3, 8)) # Extra right margin for the categorical labels

plot(source_map, 
     main = "Data Source Contribution", 
     col = source_cols,
     type = "classes", 
     plg = list(x = "right", cex = 0.8)) # Force legend to the right
plot(mainland_mask, add = TRUE, border = "black", lwd = 0.5)

# Reset layout
layout(1)
par(mar = c(5, 4, 4, 2) + 0.1)




###assess the contribution of GBIF###


library(dplyr)
library(tidyr)

# 1. FORCE UPDATE: Create the dataframe directly from the CORRECTED vector
# This ensures final_df contains only the data with 0 endemic enrichment
final_df <- as.data.frame(enriched_iberia_vect)

# 1. Clean the variable to ensure no hidden spaces
widespread_clean <- trimws(widespread_list)

# 2. Create the FINAL summary table from the vector currently in memory
# We use a brand new name 'PHD_FINAL_SUMMARY' to avoid any confusion with old tables
GBIF_additions <- as.data.frame(enriched_iberia_vect) %>%
  group_by(Sp) %>%
  summarise(
    Regional_Base = sum(Source_DB %in% c("MITECO", "SNIPAD")),
    GBIF_Gain     = sum(Source_DB == "GBIF"),
    Total_Presence_Grids = n()
  ) %>%
  mutate(
    Is_Enriched = GBIF_Gain > 0,
    Percent_Gain = round((GBIF_Gain / Regional_Base) * 100, 1),
    Percent_Gain = ifelse(is.infinite(Percent_Gain), 100, Percent_Gain)
  ) %>%
  arrange(desc(GBIF_Gain))

# 3. The Truth Check
cat("Final Count of Enriched Species:", sum(GBIF_additions$Is_Enriched), "\n")

#View(comparison_table)




library(dplyr)
library(tidyr)
library(ggplot2)

# 1. Prepare the data: Pivot to long format
plot_data_long <- comparison_table %>%
  # Filter for species that actually have data
  filter(Regional > 0 | GBIF > 0) %>%
  # Focus on the top 40 species by total count to keep the plot readable
  mutate(Total = Regional + GBIF) %>%
  slice_max(Total, n = 50) %>%
  # Transform columns into rows
  pivot_longer(cols = c(Regional, GBIF), 
               names_to = "Source", 
               values_to = "Count")

# 2. Check the object structure
# It should now have multiple rows per species: one for Regional, one for GBIF
head(plot_data_long)



ggplot(plot_data_long, aes(x = reorder(Sp, Count, sum), y = Count, fill = Source)) +
  geom_bar(stat = "identity", position = "stack", color = "white", linewidth = 0.1) +
  coord_flip() +
  # These names must EXACTLY match the names in the 'Source' column
  scale_fill_manual(values = c("Regional" = "#FFD700", "GBIF" = "#993300"),
                    labels = c("GBIF Gap-Fill", "Regional (MITECO/SNIPAD)")) +
  labs(
    title = "Contribution of GBIF Enrichment to Iberian Fish Data",
    subtitle = "Top 50 species by grid cell occupancy",
    x = "Species",
    y = "Total 10km Grid Cells",
    fill = "Data Source"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8, face = "italic"),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )




####THRESHOLD ON 10 AND THRESHOLD ON 15####




library(dplyr)

# 1. Calculate the status for each species relative to the thresholds
threshold_analysis <- comparison_table %>%
  mutate(
    # Total occupancy
    Total_Grids = Regional + GBIF,
    
    # Check if they crossed 10
    Base_10 = Regional >= 10,
    Total_10 = Total_Grids >= 10,
    Rescued_10 = (!Base_10 & Total_10),
    
    # Check if they crossed 15
    Base_15 = Regional >= 15,
    Total_15 = Total_Grids >= 15,
    Rescued_15 = (!Base_15 & Total_15)
  )

# 2. Summarize the Results
summary_stats <- data.frame(
  Threshold = c("10 Grids", "15 Grids"),
  Base_Only = c(sum(threshold_analysis$Base_10), sum(threshold_analysis$Base_15)),
  With_GBIF = c(sum(threshold_analysis$Total_10), sum(threshold_analysis$Total_15)),
  Rescued_by_GBIF = c(sum(threshold_analysis$Rescued_10), sum(threshold_analysis$Rescued_15))
)

print(summary_stats)

# 3. Identify the "Rescued" Species
rescued_15_names <- threshold_analysis %>% 
  filter(Rescued_15) %>% 
  dplyr::select(Sp, Regional, GBIF, Total_Grids)

cat("\n--- Species Rescued for SDM (Crossed 15 grid threshold) ---\n")
print(rescued_15_names)





# > print(summary_stats)
# Threshold Base_Only With_GBIF Rescued_by_GBIF
# 1  10 Grids        95        98               3
# 2  15 Grids        86        90               4
# > # 3. Identify the "Rescued" Species
#   > rescued_15_names <- threshold_analysis %>% 
#   +   filter(Rescued_15) %>% 
#   +   select(Sp, Regional, GBIF, Total_Grids)
# > cat("\n--- Species Rescued for SDM (Crossed 15 grid threshold) ---\n")
# 
# --- Species Rescued for SDM (Crossed 15 grid threshold) ---
#   > print(rescued_15_names)
# # A tibble: 4 × 4
# Sp                     Regional  GBIF Total_Grids
# <chr>                     <int> <int>       <int>
#   1 Blicca bjoerkna              12     5          17
# 2 Dicentrarchus labrax         11   132         143
# 3 Pomatoschistus microps       10    15          25
# 4 Syngnathus abaster           13    10          23


#wrap up

enriched_iberia_vect


library(dplyr)

# 1. Define the role of each species based on its source profile
final_stats <- as.data.frame(enriched_iberia_vect) %>%
  group_by(Sp) %>%
  summarise(
    Used_Regional = any(Source_DB %in% c("MITECO", "SNIPAD")),
    Used_GBIF     = any(Source_DB == "GBIF")
  ) %>%
  mutate(
    Final_Status = case_when(
      Used_Regional & Used_GBIF  ~ "Widespread (Enriched with GBIF)",
      Used_Regional & !Used_GBIF ~ "Regional Only (Endemic or No Gaps)",
      !Used_Regional & Used_GBIF ~ "GBIF Only (New Records)",
      TRUE                       ~ "Other"
    )
  )

# 2. Generate the count
final_summary_table <- final_stats %>%
  group_by(Final_Status) %>%
  summarise(Species_Count = n()) %>%
  mutate(Percentage = round((Species_Count / sum(Species_Count)) * 100, 1))

print(final_summary_table)
cat("Total unique species accounted for:", sum(final_summary_table$Species_Count), "\n")

unique(final_stats$Sp)


# Check length
length(unique(final_stats$Sp)) 

# Check for any underlying case-sensitivity issues (e.g., "Salmo trutta" vs "salmo trutta")
any(duplicated(tolower(unique(final_stats$Sp))))


#provide categories in the final dataset
library(openxlsx)
library(dplyr)
invasive_est<-read.xlsx("/Users/georgevagenas/Desktop/Vagenas_aSDMs/input/Dataset/Established_Invasive_Iberian_Soto_et_al_2025/ddi70071-sup-0002-tables3.xlsx",sep=",")
str(invasive_est)


# 2. Check the structure again
str(invasive_est)

# 3. Filter for Fishes (Actinopterygii)
# The Soto et al. list includes all animals; you only want the fish for your PhD
invasive_fishes <- invasive_est %>%
  filter(Group == "Fishes")

# 4. View the clean list
print(unique(invasive_fishes$scientificName))


# 1. Strip the authorities and years
# This keeps only the first two words (Genus and species)
invasive_names_clean <- sub("^(\\S+\\s+\\S+).*", "\\1", unique(invasive_fishes$scientificName))

# 2. Remove potential hybrids or messy entries (like 'x P. dragarum')
invasive_names_clean <- invasive_names_clean[!grepl(" x ", invasive_names_clean)]

cat("Cleaned invasive names sample:\n")
print(head(invasive_names_clean))



#Final classification

### SAVE FINAL OBJECTS ###

library(dplyr)

# 1. Clean the invasive names (Strip Authority/Year)
invasive_names_clean <- sub("^(\\S+\\s+\\S+).*", "\\1", unique(invasive_fishes$scientificName))
invasive_names_clean <- invasive_names_clean[!grepl(" x ", invasive_names_clean)]

# 2. Build the final classification dataframe
aSDMs_species_classification <- data.frame(Sp = unique(final_stats$Sp)) %>%
  mutate(
    Origin = case_when(
      Sp %in% trimws(widespread_list) ~ "Widespread",
      Sp %in% trimws(endemic_list)    ~ "Endemic",
      TRUE                            ~ "Other"
    ),
    Is_Invasive = Sp %in% invasive_names_clean,
    Category = case_when(
      Origin == "Endemic" ~ "Iberian Endemic",
      Is_Invasive == TRUE ~ "Invasive Widespread",
      TRUE                ~ "Native Widespread"
    )
  )

# 3. Save as CSV
#LOCKED
write.csv(aSDMs_species_classification, 
          "C:/Users/geo_v/Desktop/MITECO_BIOCAMBIO/Dataset/Final_enriched_GBIF_MITECO_SNIPAD_tax_correct_objects/aSDMs_species_classification.csv", 
          row.names = FALSE)


#Iberian file

library(dplyr)
library(terra)

# 1. Create the summary of grid counts from the final spatial object
occupancy_summary <- as.data.frame(enriched_iberia_vect) %>%
  group_by(Sp) %>%
  summarise(Total_Grids = n())

# 2. Merge with your classification categories
species_summary_final <- aSDMs_species_classification %>%
  left_join(occupancy_summary, by = "Sp") %>%
  # Replace NA with 0 if any species had no spatial data
  mutate(Total_Grids = ifelse(is.na(Total_Grids), 0, Total_Grids))

# 3. Save the CSV
write.csv(species_summary_final, 
          "C:/Users/geo_v/Desktop/MITECO_BIOCAMBIO/Dataset/Final_enriched_GBIF_MITECO_SNIPAD_tax_correct_objects/enriched_species_grid_counts.csv", 
          row.names = FALSE)

# 4. Save the Shapefile (as requested)

length(unique(enriched_iberia_vect$Sp)) #120 in total, the widespread where found to be 61

writeVector(enriched_iberia_vect, 
            "C:/Users/geo_v/Desktop/MITECO_BIOCAMBIO/Dataset/Final_enriched_GBIF_MITECO_SNIPAD_tax_correct_objects/enriched_iberia_vect.shp", 
            overwrite = TRUE)


#global GBIF file for the widespread (61 species)

# 1. Load the global dataset
global_dataset <- vect("/Users/georgevagenas/Desktop/Vagenas_aSDMs/input/Dataset/vect_GLOBAL_GBIF_107taxonomically_corrected/vect_GLOBAL_GBIF_107taxonomically_corrected.shp")

# 2. Identify the 61 names to filter
widespread_names <- species_summary_final %>%
  filter(Category %in% c("Invasive Widespread", "Native Widespread")) %>%
  pull(Sp)

# 3. Perform the filter
# NOTE: Ensure the column in your global shapefile is named 'Sp'
global_gbif_widespread <- global_dataset[global_dataset$Sp %in% widespread_names]

# Check if you have exactly 61 widespread species in the global file
length(unique(global_gbif_widespread$Sp)) 

# Check if you have exactly 120 species in the summary table
nrow(species_summary_final)

# Find the name that exists in your list but NOT in the global shapefile
missing_sp <- setdiff(widespread_names, unique(global_gbif_widespread$Sp))

cat("The missing species is:", missing_sp, "\n") 

#only two records it will be disregarded later on due to lack of records
# > cat("The missing species is:", missing_sp, "\n")
# The missing species is: Misgurnus mohoity 

# 4. Save the global widespread subset
writeVector(global_gbif_widespread, 
            "C:/Users/geo_v/Desktop/MITECO_BIOCAMBIO/Dataset/Final_enriched_GBIF_MITECO_SNIPAD_tax_correct_objects/global_gbif_widespread.shp", 
            overwrite = TRUE)

cat("Global widespread dataset (60 species) saved successfully.\n")
