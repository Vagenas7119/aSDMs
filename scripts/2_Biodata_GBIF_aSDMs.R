#June 2026
#Development pipeline of freshwater SDMs (aSDMS)

#Sector: Biodata_aSDMs [IBERIAN dataset with filtering and GBIF scanning]

#PhD Researcher - Georgios Vagenas (georgios.vagenas@mncn.csic.es | georgvagenas@gmail.com)

#### Pre-setting :: Libraries required to perform the analysis ####


library(terra)
#iberian MITECO-SNIPAD dataset - 122 species
vect_dataset_IBERIA<-vect("/Users/georgevagenas/Desktop/Vagenas_aSDMs/input/Dataset/vect_MITECO_SNIPAD_IBERIA/vect_dataset_IBERIA.shp")


#run taxonomic and habitat checking


# Install packages if you don't have them:
# install.packages(c("rfishbase", "worrms", "dplyr", "stringr", "purrr"))

library(rfishbase)
library(worrms)
library(dplyr)
library(stringr)
library(purrr)
library(terra)
library(remotes)
library(duckdb)

# Install rvest if you don't have it (this is R's premier web scraping tool)
if (!requireNamespace("rvest", quietly = TRUE)) install.packages("rvest")

library(rvest)
library(dplyr)
library(stringr)
library(purrr)

# 1. Get your unique species list
species_list <- unique(vect_dataset_IBERIA$Sp)

# 2. Define the Web Scraper Function
scrape_lentic_habitat <- possibly(function(sp) {
  # Format the name for the URL (e.g., "Cyprinus carpio" -> "Cyprinus-carpio")
  url_name <- str_replace_all(sp, " ", "-")
  url <- paste0("https://www.fishbase.se/summary/", url_name, ".html")
  
  # Download and read the live webpage HTML
  page_html <- read_html(url)
  
  # Extract all text from the page and convert to lowercase
  page_text <- page_html %>% html_text(trim = TRUE) %>% tolower()
  
  # Search the page text for our explicitly lentic keywords
  is_lentic <- str_detect(
    page_text, 
    "lake|pond|reservoir|lentic|standing water|stagnant"
  )
  
  data.frame(
    Species = sp,
    Is_Explicitly_Lentic = is_lentic,
    FishBase_URL = url
  )
}, otherwise = data.frame(Species = NA, Is_Explicitly_Lentic = NA, FishBase_URL = NA))

# 3. Run the Scraper (This will take ~30 to 60 seconds depending on internet speed)
#message("Scraping live FishBase website for ", length(species_list), " species. Please wait...")
#fb_lentic_data <- map_df(species_list, scrape_lentic_habitat)

#fb_lentic_data


# Save the data frame as a CSV file
#write.csv(fb_lentic_data, "C:/Users/geo_v/Desktop/MITECO_BIOCAMBIO/Dataset/Fishbase_Web_Scraping_122sp_MITECO_SNIPAD/Iberia_Fish_Lentic_Assessment.csv", row.names = FALSE)


fb_lentic_data<-read.csv("/Users/georgevagenas/Desktop/Vagenas_aSDMs/input/Dataset/Fishbase_Web_Scraping_122sp_MITECO_SNIPAD/Iberia_Fish_Lentic_Assessment.csv")

message("Success! File saved to your working directory.")

#Save it#

# 4. Clean up and view the results
fb_lentic_clean <- fb_lentic_data %>% filter(!is.na(Species))

cat("\n=======================================================\n")
cat("          FISHBASE LENTIC ASSESSMENT SUMMARY           \n")
cat("=======================================================\n")
cat("Total species assessed:", nrow(fb_lentic_clean), "\n")
cat("Number of explicitly Lentic species:", sum(fb_lentic_clean$Is_Explicitly_Lentic, na.rm = TRUE), "\n")
cat("=======================================================\n\n")

# Show just the lentic ones
lentic_species <- fb_lentic_clean %>% filter(Is_Explicitly_Lentic == TRUE)
print(lentic_species$Species)

#false alarm, there are no exclussively lentic in the Iberian peninsula


### tomorrow assess endemics in IB, and widespread and invasive from Soto et al. 2025 ###

###### GBIF download ######

# 1. Get your unique species list
species_list <- unique(vect_dataset_IBERIA$Sp)



# ==============================================================================
# 0. LOAD REQUIRED LIBRARIES
# ==============================================================================
#install.packages(c("rgbif", "CoordinateCleaner", "dplyr", "openxlsx", "purrr"))
library(rgbif)
library(CoordinateCleaner)
library(dplyr)
library(openxlsx)
library(purrr)

# ==============================================================================
# 1. GBIF AUTHENTICATION & TAXONOMY
# ==============================================================================
# CRITICAL: You must provide your GBIF credentials to request a formal download
user_gbif  <- "georgios_vagenas"
pwd_gbif   <- "GBIF@terra_X"
email_gbif <- "georgvagenas@gmail.com"

# Assuming 'species_list' exists in your environment from the previous script
message("Matching ", length(species_list), " species to GBIF backbone...")

# Get exact GBIF Taxon Keys for our species list
taxon_match <- name_backbone_checklist(name_data = data.frame(name = species_list))
taxon_keys <- taxon_match$usageKey[!is.na(taxon_match$usageKey)]

# ==============================================================================
# 2. TRIGGER THE GBIF DOWNLOAD
# ==============================================================================
# Define your environmental data temporal window (Update these years!)
start_year <- 1979 
end_year   <- 2025 

message("Triggering GBIF server download. This generates your DOI...")

#LOCKED
# gbif_download <- occ_download(
#   pred_in("taxonKey", taxon_keys),
#   pred("hasCoordinate", TRUE),
#   pred("hasGeospatialIssue", FALSE),
#   pred("occurrenceStatus", "PRESENT"),
#   pred_in("basisOfRecord", c("HUMAN_OBSERVATION", "OBSERVATION", "OCCURRENCE")),
#   pred_gte("year", start_year),
#   pred_lte("year", end_year),
#   format = "SIMPLE_CSV",
#   user = user_gbif, pwd = pwd_gbif, email = email_gbif
# )

# Wait for the GBIF servers to compile the data (can take 10 mins to several hours)
#LOCKED
#occ_download_wait(gbif_download)

# Import the completed dataset into R
raw_gbif_data <- occ_download_get(gbif_download) %>%
  occ_download_import()

# Load the readr package (install it with install.packages("readr") if you don't have it)
library(readr)

message("Saving 13 million records... This might take a few minutes!")

# Save the raw dataframe to a CSV file without row names #LOCKED
#write_csv(raw_gbif_data, "C:/Users/geo_v/Desktop/MITECO_BIOCAMBIO/Dataset/GBIF_2026_dataset/Raw_GBIF_Download_13M.csv")

raw_gbif_data<-read.csv("/Users/georgevagenas/Desktop/Vagenas_aSDMs/input/Dataset/GBIF_2026_dataset/Raw_GBIF_Download_13M.csv")


message("Success! Raw data saved to your working directory.")

# The correct way to extract the DOI from the meta-object #LOCKED
#master_doi <- occ_download_meta(gbif_download)$doi 

message("Master DOI successfully saved as: ", master_doi)

# ==============================================================================
# 3. COORDINATE CLEANING
# ==============================================================================


#LOAD LIBRARIES
library(dplyr)
library(CoordinateCleaner)
library(rgbif)
library(purrr)
library(openxlsx)

# If you don't have data.table installed, run: install.packages("data.table")
library(data.table) 

# ==============================================================================
# 1. SPATIAL PRE-DEDUPLICATION & CLEANING (ULTRA-FAST METHOD)
# ==============================================================================
message("1. Extracting unique spatial coordinates to save time...")

#LOCKED
# Extract only unique locations to minimize CoordinateCleaner's workload
#unique_coords <- raw_gbif_data %>%
  distinct(species, decimalLatitude, decimalLongitude)

cat("Reduced CoordinateCleaner workload to", nrow(unique_coords), "unique locations.\n")

message("2. Running CoordinateCleaner (Capitals test removed)...")

clean_flags <- clean_coordinates(
  
  x = unique_coords,
  
  lon = "decimalLongitude",
  
  lat = "decimalLatitude",
  
  species = "species",
  
  # CRITICAL LINE: This tells R to ONLY run the fast tests. 
  
  # Because "capitals" and "outliers" are missing from this list, R skips them!
  
  tests = c("centroids", "equal", "gbif","seas", "zeros"), 
  
  centroids_rad = 1000
  
)



#records removed
# Testing coordinate validity
# Flagged 0 records.
# Testing equal lat/lon
# Flagged 0 records.
# Testing zero coordinates
# Flagged 23 records.
# Testing country centroids
# Flagged 1675 records.
# Testing sea coordinates
# Reading ne_50m_land.zip from naturalearth...
# Flagged 73427 records.
# Testing GBIF headquarters, flagging records around Copenhagen
# Flagged 9 records.
# Flagged 75097 of 1268940 records, EQ = 0.06.
#



# Keep only the unique coordinates that passed the tests
valid_coords <- unique_coords[clean_flags$.summary == TRUE, ]

# ==============================================================================
# 2. RE-MERGE & TEMPORAL DEDUPLICATION
# ==============================================================================
message("3. Merging safe coordinates and applying temporal deduplication...")

# Inner join instantly drops all bad rows from the 13M dataset, 
# then we filter down to 1 record per species/pixel/month/year
final_data <- raw_gbif_data %>%
  inner_join(valid_coords, by = c("species", "decimalLatitude", "decimalLongitude")) %>%
  distinct(species, decimalLatitude, decimalLongitude, year, month, .keep_all = TRUE)

cat("Original records: ", nrow(raw_gbif_data), "\n")
cat("Cleaned & Deduplicated records: ", nrow(final_data), "\n")


# > cat("Original records: ", nrow(raw_gbif_data), "\n")
# Original records:  12998960 
# > cat("Cleaned & Deduplicated records: ", nrow(final_data), "\n")
# Cleaned & Deduplicated records:  2162377 

# Save final spatial dataset using data.table for maximum speed
message("4. Saving Cleaned Dataset to CSV...")
fwrite(final_data, "C:/Users/geo_v/Desktop/MITECO_BIOCAMBIO/Dataset/GBIF_2026_dataset/Global_Iberia_GBIF_Cleaned_Occurrences.csv")


final_data<-read.csv("/Users/georgevagenas/Desktop/Vagenas_aSDMs/input/Dataset/GBIF_2026_dataset/Global_Iberia_GBIF_Cleaned_Occurrences.csv")
# ==============================================================================
# 3. EXTRACT DOIs TO EXCEL
# ==============================================================================
message("5. Generating Dataset DOI report...")

# Create a unique list of which species came from which specific dataset key
dataset_species_map <- final_data %>%
  dplyr::select(datasetKey, species) %>%
  distinct()

# Get the unique dataset keys to query the API
unique_datasets <- unique(dataset_species_map$datasetKey)


unique_species<-unique(dataset_species_map$species)
unique_species


# Safely query GBIF to get the DOI for each specific dataset key
#LOCKED
#safe_get_doi <- possibly(function(key) {
  meta <- datasets(uuid = key)
  data.frame(
    datasetKey = key,
    Dataset_DOI = ifelse(!is.null(meta$data$doi), meta$data$doi, "No DOI Provided"),
    Dataset_Title = ifelse(!is.null(meta$data$title), meta$data$title, "Unknown Title")
  )
}, otherwise = NULL)

# Fetch all DOIs
dataset_metadata <- map_df(unique_datasets, safe_get_doi)

# Merge the DOIs back with the species list
doi_report <- dataset_species_map %>%
  left_join(dataset_metadata, by = "datasetKey") %>%
  select(Species = species, Dataset_DOI, Dataset_Title, GBIF_Dataset_Key = datasetKey) %>%
  arrange(Species)

# Export to XLSX
write.xlsx(doi_report, "C:/Users/geo_v/Desktop/MITECO_BIOCAMBIO/Dataset/GBIF_2026_dataset/GBIF_Source_DOIs_by_Species.xlsx")
message("Success! DOIs saved to 'GBIF_Source_DOIs_by_Species.xlsx'")


#spatial thin the GBIF dataset to the global_predictor setting#

str(final_data)

predictors_global<-rast("/Users/georgevagenas/Desktop/Vagenas_aSDMs/input/Dataset/predictors_finalized_global.tiff")
str(predictors_global$BIO5_clima)
plot(predictors_global$BIO5_clima)
globext<-predictors_global$BIO5_clima
globext





#### remove points outside of the predictors_global grid extent####


library(terra)
library(dplyr)

message("1. Grabbing coordinates...")
# We must use exactly X (Longitude) then Y (Latitude) for terra
coords <- final_data %>% dplyr::select(decimalLongitude, decimalLatitude)

message("2. Checking the actual pixels underneath each point...")
# This checks the physical grid cells. 
# If a point is in the ocean, it gets an NA. If it's on land, it gets the BIO5 value.
pixel_values <- terra::extract(globext, coords, ID = FALSE)

message("3. Deleting points that fell completely outside your valid grid cells...")
final_data_masked <- final_data %>%
  # Attach the extracted pixel values to your data
  bind_cols(pixel_values) %>%
  # THE MAGIC LINE: Drop any row where the pixel was NA (ocean/empty grid)
  filter(!is.na(BIO5_clima)) %>%
  # Drop the climate column so your dataset structure looks exactly like the original
  dplyr::select(-BIO5_clima) 

# Let's see how many oceanic/out-of-bounds points we killed!
cat("Original points:", nrow(final_data), "\n")
cat("Points kept (inside valid pixels):", nrow(final_data_masked), "\n")
cat("Points deleted (fell on NA pixels):", nrow(final_data) - nrow(final_data_masked), "\n")

# > cat("Original points:", nrow(final_data), "\n")
# Original points: 2162377 
# > cat("Points kept (inside valid pixels):", nrow(final_data_masked), "\n")
# Points kept (inside valid pixels): 2079968 
# > cat("Points deleted (fell on NA pixels):", nrow(final_data) - nrow(final_data_masked), "\n")
# Points deleted (fell on NA pixels): 82409 


#Spatially thinning species at the grid level


library(terra)
library(dplyr)

message("1. Identifying which grid cell every point belongs to...")
# Extract the coordinates as a matrix (terra prefers matrices for this function)
coords <- final_data_masked %>% dplyr::select(decimalLongitude, decimalLatitude) %>% as.matrix()

# Calculate the exact grid cell ID number for every single point
final_data_cells <- final_data_masked %>%
  mutate(grid_cell = terra::cellFromXY(globext, coords))

message("2. Spatially thinning: Keeping 1 record per species, per cell...")
# Group by the species and the grid cell ID, then keep just 1 row
thinned_data <- final_data_cells %>%
  group_by(species, grid_cell) %>%
  # Using slice_sample randomly picks 1 point rather than always taking the first one. 
  # This prevents time-based or observer-based bias!
  slice_sample(n = 1) %>% 
  ungroup()

message("3. Snapping points to the grid centroids...")
# Calculate the exact X and Y coordinates for the center of the remaining grid cells
cell_centroids <- terra::xyFromCell(globext, thinned_data$grid_cell)

# Overwrite the original coordinates with the perfect centroid coordinates
final_thinned_data <- thinned_data %>%
  mutate(
    decimalLongitude = cell_centroids[, "x"],
    decimalLatitude  = cell_centroids[, "y"]
  ) %>%
  dplyr::select(-grid_cell) # Remove the temporary cell ID column to keep your data clean

# Let's see the results!
cat("Records before thinning: ", nrow(final_data_masked), "\n")
cat("Final spatially thinned records: ", nrow(final_thinned_data), "\n")
cat("Total duplicate points removed: ", nrow(final_data_masked) - nrow(final_thinned_data), "\n")

#assess how many species, how many grids per species etc.


library(dplyr)

message("Generating Species Distribution Report...")

# 1. Create a detailed summary table
species_grid_report <- final_thinned_data %>%
  # Flag blank species names (records identified only to Genus or Family)
  mutate(species = ifelse(is.na(species) | trimws(species) == "", "[Unidentified to Species Level]", species)) %>%
  group_by(species) %>%
  summarise(Number_of_Grids = n(), .groups = "drop") %>%
  arrange(desc(Number_of_Grids)) # Sort from most widespread to least

# 2. Separate valid species from unidentified records for accurate math
valid_species <- species_grid_report %>% 
  filter(species != "[Unidentified to Species Level]")

unidentified_count <- species_grid_report %>% 
  filter(species == "[Unidentified to Species Level]") %>% 
  pull(Number_of_Grids)

# Handle cases where there are no unidentified records
if(length(unidentified_count) == 0) unidentified_count <- 0

# 3. Calculate Summary Statistics
total_species <- nrow(valid_species)
total_valid_grids <- sum(valid_species$Number_of_Grids)
max_grids <- max(valid_species$Number_of_Grids)
min_grids <- min(valid_species$Number_of_Grids)
median_grids <- median(valid_species$Number_of_Grids)

# 4. Print the Dashboard to the Console
cat("\n======================================================\n")
cat("          SPATIAL THINNING & SPECIES REPORT\n")
cat("======================================================\n")
cat("Total Valid Species:            ", total_species, "\n")
cat("Total Valid Grid Occurrences:   ", total_valid_grids, "\n")
cat("Unidentified/Genus-only Grids:  ", unidentified_count, "\n")
cat("------------------------------------------------------\n")
cat("Max Grids (Most widespread):    ", max_grids, "\n")
cat("Min Grids (Rarest species):     ", min_grids, "\n")
cat("Median Grids per Species:       ", median_grids, "\n")
cat("======================================================\n\n")

# 5. Show the top 10 most widespread species
message("Top 10 Most Widespread Species:")
print(head(species_grid_report, 10))


#### I queried for 122 and i received 172, i wanna assess why ####



message("1. Identifying the Intruders...")
# Find which species are in your GBIF data but NOT in your original list
intruder_species <- setdiff(valid_species$species, species_list)

cat("Found", length(intruder_species), "intruder species.\n")

message("2. Filtering dataset strictly to your 122 target species...")
# Keep ONLY rows where the 'species' column exactly matches your original 122 names
final_strict_data <- final_thinned_data %>%
  filter(species %in% species_list)

# Verify the final count
strict_species_count <- length(unique(final_strict_data$species))
cat("Total species remaining:", strict_species_count, "\n")
cat("Total model-ready grids:", nrow(final_strict_data), "\n")


#save the dataset check species richness map without antartica and greenland with the certain figure

library(data.table)
library(terra)
library(ggplot2)
library(sf) # Using sf makes projecting polygons in ggplot much easier!

# ==============================================================================
# 1. SAVE THE FINAL DATASET
# ==============================================================================
# Save as a lightning-fast CSV for your permanent records
fwrite(final_strict_data, "C:/Users/geo_v/Desktop/MITECO_BIOCAMBIO/Dataset/GBIF_final_strict_data_107_thinned_cleaned_373975/GBIF_final_strict_data_107_thinned_cleaned_373975.csv")
message("Saved final dataset to CSV.")