#June 2026
#Development pipeline of freshwater SDMs (aSDMS)

#Sector: Figure 1 - Global & Regional database

#PhD Researcher - Georgios Vagenas (georgios.vagenas@mncn.csic.es | georgvagenas@gmail.com)

#### Pre-setting :: Libraries required to perform the analysis ####


final_strict_data<-read.csv("/Users/georgevagenas/Desktop/Vagenas_aSDMs/input/Dataset/GBIF_final_strict_data_107_thinned_cleaned_373975/GBIF_final_strict_data_107_thinned_cleaned_373975.csv")

str(final_strict_data)

library(dplyr)
library(terra)

message("1. Creating the blank 50 arc-minute grid...")
blank_raster <- rast(ext(-180, 180, -90, 90), res = 50/60, crs = "EPSG:4326")

message("2. Assigning points to the new, coarser grid cells...")
coords <- final_strict_data %>% dplyr::select(decimalLongitude, decimalLatitude) %>% as.matrix()
cell_ids <- terra::cellFromXY(blank_raster, coords)

message("3. Calculating TRUE species richness per cell...")
richness_counts <- final_strict_data %>%
  # Attach the new cell IDs
  mutate(cell_50 = cell_ids) %>%
  # Drop any points that somehow fell off the map
  filter(!is.na(cell_50)) %>%
  # CRITICAL: Group by cell AND species, so each species only counts ONCE per giant cell
  group_by(cell_50, species) %>%
  summarise(presence = 1, .groups = "drop") %>%
  # Now, count the number of unique species per giant cell
  group_by(cell_50) %>%
  summarise(richness = n(), .groups = "drop")

message("4. Building the final accurate raster...")
# Create a copy of the blank raster and fill it with our accurate counts
richness_raster <- copy(blank_raster)
values(richness_raster) <- NA # Ensure it is completely empty first
richness_raster[richness_counts$cell_50] <- richness_counts$richness
names(richness_raster) <- "richness"

# Let's verify your catch!
summary(values(richness_raster, na.rm = TRUE))



# ==============================================================================
# 3. APPLY THE ROBINSON PROJECTION
# ==============================================================================


# Install the needed packages if you don't have them
# install.packages(c("rnaturalearth", "rnaturalearthdata", "sf", "dplyr", "terra"))

library(rnaturalearth)
library(sf)
library(dplyr)
library(terra)

message("1. Fetching Global Shapefiles...")
# Download the global country boundaries at medium resolution
world_sf <- ne_countries(scale = "medium", returnclass = "sf")

message("2. Creating 'land_no_antarctica'...")
# Filter out Antarctica and Greenland
world_filtered_sf <- world_sf %>%
  filter(!admin %in% c("Antarctica", "Greenland"))

# Convert to a terra SpatVector (which your previous code expects)
land_no_antarctica <- vect(world_filtered_sf)

message("3. Creating 'iberian_shp'...")
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

#plot(iberian_shp)
message("Success! Both shapefiles are now loaded in your environment.")



# Define the Robinson projection string
robinson_crs <- "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

# Project the raster to Robinson
richness_rob <- project(richness_raster, robinson_crs)

# Project your background polygons (Assuming land_no_antarctica & iberian_shp are loaded)
# Converting them to 'sf' objects makes ggplot render them perfectly
land_sf <- st_as_sf(project(land_no_antarctica, robinson_crs))
iberian_sf <- st_as_sf(project(iberian_shp, robinson_crs))

# Extract the projected raster values to a dataframe for ggplot
richness_df <- as.data.frame(richness_rob, xy = TRUE, na.rm = TRUE)
summary(richness_df$richness)
# ==============================================================================
# 4. PLOT THE GLOBAL MAP
# ==============================================================================
# Define your specific color palette
yellow_to_very_dark_orange <- colorRampPalette(c("yellow", "#993300"))(50)

global_map <- ggplot() +
  # 1. The grey background of the world
  geom_sf(data = land_sf, fill = "lightgrey", color = "grey", size = 0.4) +
  
  # 2. The species richness raster tiles
  geom_tile(data = richness_df, aes(x = x, y = y, fill = richness)) +
  
  # 3. The Iberian Peninsula outline on top
  geom_sf(data = iberian_sf, fill = NA, color = "black", linewidth = 0.4) +
  
  # 4. The custom gradient scale
  scale_fill_gradientn(
    colors = yellow_to_very_dark_orange,
    name = "Species richness",
    limits = c(0, 50),
    breaks = c(0, 10, 20, 30, 40, 50),
    labels = c("0", "10", "20", "30", "40", "50"),
    na.value = "transparent",
    guide = guide_colorbar(
      direction = "horizontal",
      label.position = "bottom",
      title.position = "top",
      title.hjust = 0.5,           # <--- THIS CENTERS THE TITLE
      barwidth = unit(6, "cm"),
      barheight = unit(0.5, "cm")
    )
  ) +
  
  # 5. Clean, minimal theme
  theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.grid.major = element_line(color = "gray90", linetype = "dashed"),
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  labs(x = NULL, y = NULL)

# Render the map
plot(global_map)


#### now add the 107 species found globally with the 122 found locally together in two species richness maps

library(sf)
library(dplyr)
library(terra)

# ==============================================================================
# 1. THE BULLETPROOF CLEAN-UP
# ==============================================================================
message("Deep cleaning geometries and attributes...")

vect_dataset_IBERIA<-vect("/Users/georgevagenas/Desktop/Vagenas_aSDMs/input/Dataset/vect_MITECO_SNIPAD_IBERIA/vect_dataset_IBERIA.shp")

# Convert the terra vector to an sf object so we can use dplyr on it
sf_IBERIA <- st_as_sf(vect_dataset_IBERIA)

# Filter out both spatial ghosts AND attribute ghosts
sf_clean <- sf_IBERIA %>%
  filter(!st_is_empty(.)) %>%       # Drops any row missing physical coordinates
  filter(!is.na(presence))          # Drops any row where 'presence' is NA

# Convert back to a terra SpatVector for the rasterize function
clean_vect_IBERIA <- vect(sf_clean)

cat("Cleaned records ready for mapping: ", nrow(clean_vect_IBERIA), "\n")

# ==============================================================================
# 2. CALCULATE TRUE RICHNESS
# ==============================================================================
message("Calculating richness from the cleaned MITECO dataset...")

# This will now run perfectly because geometries and attributes are 1:1 matched
richness_raster_miteco <- rasterize(
  x = clean_vect_IBERIA, 
  y = spatlayer, 
  field = "presence", 
  fun = "sum"
)

# Mask the MITECO raster to perfectly fit the coastline (dropping ocean grids)

library(patchwork)
regional_raster_masked <- mask(richness_raster_miteco, iberian_shp)

plot(regional_raster_masked)

# Convert to a dataframe for ggplot
regional_df <- as.data.frame(regional_raster_masked, xy = TRUE, na.rm = TRUE)
colnames(regional_df)[3] <- "Richness"

message("Success! Your regional_df is perfectly prepared.")


# ==============================================================================
# 3. BUILD THE REGIONAL INSET MAP (Map B) - TIGHT MARGINS
# ==============================================================================
message("Styling Map B with auto-tightened extent...")


iberia_mainland_wgs84<-st_as_sf(iberian_shp)

regional_map <- ggplot() +
  geom_sf(data = iberia_mainland_wgs84, fill = "white", color = NA) +
  geom_tile(data = regional_df, aes(x = x, y = y, fill = Richness)) +
  geom_sf(data = iberia_mainland_wgs84, fill = NA, color = "gray25", linewidth = 0.6) +
  
  scale_fill_gradientn(
    colors = yellow_to_very_dark_orange,
    limits = c(0, 50),
    breaks = c(0, 10, 20, 30, 40, 50),
    labels = c("0", "10", "20", "30", "40", "50"),
    na.value = "transparent"
  ) +
  
  # THE FIX: By removing xlim and ylim, we force ggplot to tightly hug the 
  # exact boundaries of the iberia_mainland_wgs84 shapefile with zero wasted space!
  coord_sf(expand = FALSE) +
  theme_minimal() +
  theme(
    plot.background = element_rect(fill = "#FFFFFF99", color = "black", linewidth = 0.4),
    panel.background = element_blank(),
    
    # We can keep the margins tight now that the extent is fixed
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5, unit = "pt"), 
    
    axis.title = element_blank(),
    axis.text = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", size = 10, hjust = 0.5, margin = margin(b = 4)),
    legend.position = "none" 
  )


# ==============================================================================
# 4. COMBINE GLOBAL (Map A) + REGIONAL INSET (Map B)
# ==============================================================================


# Now, rerun Step 4 to combine them!
final_composite_map <- global_map + 
  inset_element(
    p = regional_map,
    # 3. POSITIONING: Pushed hard left and lower
    left   = -0.01,  # Flush against the absolute left edge of the map panel
    bottom = 0.05,  # Lowered (if it overlaps your global legend, raise this to 0.08)
    right  = 0.33,  # Adjusted to keep the box proportions looking natural
    top    = 0.5,  
    align_to = "panel"
  )

print(final_composite_map)


##Save figure 1##

# Save as a 600 DPI PNG
ggsave(
  filename = "C:/Users/geo_v/Desktop/MITECO_BIOCAMBIO/Dataset/Figures/Iberian_Richness_Composite_Map.png", 
  plot = final_composite_map, 
  width = 12,       # Width in inches
  height = 8,       # Height in inches
  dpi = 600,        # High resolution!
  bg = "white"      # CRITICAL: Prevents the background from saving as transparent/black
)

message("High-res PNG saved successfully!")

#and then save as GLOBAL GBIF for the global and then merge with the regional
#for the reg




