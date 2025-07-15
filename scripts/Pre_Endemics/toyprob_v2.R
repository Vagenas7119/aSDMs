#July 2025
#Development pipeline of aquatic SDMs (aSDMS)

#Pre-Endemics | Toy Exercise 

#PhD Researcher - Georgios Vagenas (georgios.vagenas@mncn.csic.es | georgvagenas@gmail.com)

#### Pre-setting :: Libraries required to perform the analysis ####

library(sf)
library(sdm)
library(dismo)
library(dplyr)
library(tidyr)
library(mapview)
library(geodata)
library(raster)
library(RColorBrewer)
library(terra)
library(usdm) 
library(sdm)
library(randomForest)
library(parallel)

# List of required packages
pkgs <- c("sf", "sdm", "dismo", "dplyr", "tidyr", 
          "mapview", "geodata", "raster", "RColorBrewer",
          "terra", "usdm", "randomForest", "parallel")

# Install missing packages in one command
if (length(setdiff(pkgs, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(pkgs, rownames(installed.packages())))
}

# Load all packages silently
invisible(lapply(pkgs, library, character.only = TRUE))

#### Chapter 1 :: Species occurrences imported objects/files to perform the aSDMs ####

#set working directory
setwd("C:/Users/geo_v/Desktop/Metanalysis/toy_problem_aSDMs/")

#import the occurence data
miteco_sp<-vect("sevensp_toyprob.shp") #The collection and the rectification of MITECO data points are described in previous versions

miteco_sp #The dataset in total with all the occurrences for 10 species of the Iberian peninsula

terra_df<-miteco_sp

# Get unique species IDs
unique_species <- unique(terra_df$species_id)

# Count the number of unique species
numbers_sp <- length(unique_species)

#### Chapter 2 :: Environmental spatial layers imported objects/files to perform the aSDMs ####

#In this case-study we will proceed in the use of WorldClim (Bioclimatic) and RiverATLAS (Hydrological) data

#Information for WorldClim is provided, but RiverATLAS layers need to be modified and resampled to 1km (equal to Worldclim resolution)

#WorldClima top 5 vars (1km)
bioclim_global<-worldclim_global(var = "bio",res = 0.5,path = "data/")



# #Optional, but in our case required step
# #Resample the RiverATLAS (vector) as a more coarsed resolution to match the one with WorldClim (1km)
# 
# #Add RIVER ATLAS as a layer # PREPARATION FROM MULTIPLE SHAPEFILES TO ONE GEOTIFF
# 
# #Initial data preparation
# # Step 1: Conversions of the multishapefile to rasterbrick
# rivers_glob<-terra::vect("data/RIVER_ATLAS_GLOBAL/RIVER_ATLAS.shp")
# 
# # Step 2: Create an empty list to store individual rasters
# rasters_list <- list()
# 
# #names
# names(bioclim_global)
# 
# #how to manually select which bioclim vars will be included
# selected_vars <- c("wc2.1_30s_bio_15", "wc2.1_30s_bio_13", 
#                    "wc2.1_30s_bio_1", "wc2.1_30s_bio_4","wc2.1_30s_bio_3")
# 
# # Subset the RasterStack to only include the selected variables
# bioc_in <- bioclim_global[[selected_vars]]
# 
# # Step 3: Resample with a loop over each field (column) in the shapefile and rasterize
# 
# # We assume that the shapefile has several attribute fields (columns) to rasterize
# #fields <- names(rivers_glob)
# #however we need only 5 hydrological variabels, based on the variables selected across pca 1 to pca 4, for info check below final part (Variable Selection with PCA)
# #field_names<-names(rivers_glob)
# 
# #These are the selected variables
# fields<- c("lka_pc_cse","dor_pc_pva","sgr_dk_rav","run_mm_cyr","dis_m3_pmn")
#     
# for (field in fields) {
# #Step 3a: Create a template raster based on the extent and resolution of the shapefile
# template_raster <- terra::rast(bioc_in[[1]]) # Adjust resolution if needed
# # 
# # Step 3b: Rasterize each field (attribute) into the raster template
# raster_layer <- terra::rasterize(rivers_glob, template_raster, field = field)
# # 
# # Step 3c: Add the raster to the list
# rasters_list[[field]] <- raster_layer
# }
# # 
# # # Step 4: Combine all individual rasters into a SpatRaster
# rasters_rast<-rast(rasters_list)
# 
# #Inspect the structure 
# str(rasters_rast)
# nlyr(rasters_rast)
# 
# # # Step 5: Save the SpatRaster as .tif
# output_file<-"data/River_atlas_1km_5vars/top5_global_rivers_atlas_raster_list_30sec.tif"
# writeRaster(rasters_rast,output_file)

#WorldClim top 5 vars (1km)
selected_vars <- c("wc2.1_30s_bio_15", "wc2.1_30s_bio_13", 
                   "wc2.1_30s_bio_1", "wc2.1_30s_bio_4","wc2.1_30s_bio_3")

#WorldClim top 5 vars (1km)
bioclim_global<-bioclim_global[[selected_vars]]

#RiverATLAS top 5 vars (1km)
riveratl_global<-rast("data/River_atlas_1km_5vars/top5_global_rivers_atlas_raster_list_30sec.tif") #within the raster there are -999 values in stream gradient but only for greenland, where we do not expect to received data points

#This step is to crop the bioclim layer to the Iberian Peninsula extent - That part has been already done
# #Bring it to the spatial extent of the Iberian
# extent_df <- read.csv("data/extent_iberian_coordinates.csv") #UPSCALE ADJUSTMENTS p1] CHANGE THE SPATIAL EXTENT
# extent <- ext(extent_df$as.vector.extent_ib.[1], extent_df$as.vector.extent_ib.[2], extent_df$as.vector.extent_ib.[3], extent_df$as.vector.extent_ib.[4])
# 
# 
# 
# #Important, step :: obtain values for climate only for the grids of the river network
# 
# bioclim_global_rn<-crop(bioclim_global,riveratl_global,mask=TRUE) #rn=river network
# 
# #requires significant amount of team therefore I am storing the raster file
# output_file<-"data/BIOCLIM_rn_1km_5vars/bioclim_global_rn_5vars_30sec.tif"
# writeRaster(bioclim_global_rn,output_file)

bioclim_global_rn<-rast("data/BIOCLIM_rn_1km_5vars/bioclim_global_rn_5vars_30sec.tif")



#### Chapter 3 :: Initialize the objects and provide the training extent to perform the aSDMs ####

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
  maxKappa=numeric(),
  maxTSS=numeric(),
  obs_prevalence=numeric(),
  stringsAsFactors = FALSE
)

#create an object
raster_list<-list()

# Initialize an empty SpatRaster object
combined_rasters <- rast()

#minimum_background_points=2

sdm_d<-list()

#### Chapter 4 :: Allocate species occurrences to each watershed in the Iberian Peninsula training extent of modelling in the loop ####

#watersheds
watershed<-vect("data/H5_Iberian/H5_Iberian.shp") #Based on the level of analysis Should change#freshwater watesheds/ecoregions of the globe

# Check if species_id is a factor and convert it to numeric correctly
if (is.factor(terra_df$species_id)) {
  terra_df$species_id <- as.numeric(levels(terra_df$species_id))[terra_df$species_id]
} else {
  terra_df$species_id <- as.numeric(terra_df$species_id)
}

# Set CRS if not done yet (assuming WGS84)
crs(terra_df) <- "EPSG:4326"

# Convert to SpatVector
miteco_sp_vect <- terra_df

#intersect occurences with ecoregions
miteco_sp_in<-terra::intersect(miteco_sp_vect,watershed)

#provide the setting

#convert miteco to data.frame works better
data <- as.data.frame(miteco_sp_in)

# terra conversion
coords <- geom(miteco_sp_in)  # Extract coordinates

# Remove the geometry column from the data
data$geometry <- NULL

# Combine coordinates and attribute data
terra_df_o<-terra_df

terra_df_o <- cbind(coords, data)

#create an id for each species
terra_df_or <- terra_df_o %>%
  mutate(species_id = as.numeric(factor(terra_df_o$species_id)))

#create a ledger only with Species and species_id (unique)
species_name_id<-dplyr::select(terra_df_or,"species_id","HYBAS_ID")

#Remove species less than 10 obs
terra_df_filtered <- terra_df_or %>%
  group_by(species_id) %>%  # Group by species
  filter(n() >= 10) %>%     # Keep species with 500 or more observations
  ungroup() %>%              # Ungroup after filtering
  droplevels()               # Drop unused factor levels, if necessary


#check the remain species
length(levels(as.factor(terra_df_filtered$species_id)))

#select only essential vectors
terra_df<-terra_df_filtered

terra_df<-dplyr::select(terra_df,"x","y","species_id","HYBAS_ID")

#rename to preserve variables
terra_df_conv<-terra_df

#provide a species presence absence table to watersheds
presence_absence_matrix <- terra_df_conv %>%
  mutate(presence = 1) %>%
  dplyr::select(species_id, HYBAS_ID, presence) %>%
  distinct() %>%
  pivot_wider(names_from = HYBAS_ID, values_from = presence, values_fill = list(presence = 0)) %>% as.data.frame()

#order ascendingly
presence_absence_matrix <- presence_absence_matrix[order(presence_absence_matrix$species_id),]

# Convert presence-absence table to long format
presence_long <- presence_absence_matrix %>%
  pivot_longer(-species_id, names_to = "HYBAS_ID", values_to = "presence")

# Create a temporary raster for masking
masked_watershed<- list()
polygon_list<- list()
raster_list<-list()

# Initialize an empty SpatRaster object
combined_rasters_stack <- rast()

#Initialize an empty SpatRaster for the extented spatial outputs
super_iberian<-rast()

#Selection only species names log and lat
terra_df<-terra_df[,1:3]

#Provide crs
spatial_points <- SpatialPoints(terra_df[,1:2], proj4string = CRS("EPSG:4326"))

# Combine with the original data frame
spatial_points_df <- SpatialPointsDataFrame(spatial_points, data = terra_df)
terra_df <-spatial_points_df

#Provide the Iberian Shapefile from watersheds
#iberian_watersheds<-vect("data/H5_Iberian/H5_Iberian.shp") #Iberian extent, do not change
iberian_watersheds<-watershed

#Dissolve the vector
aggregated_vec<-aggregate(iberian_watersheds)

#Create an empty raster with NA values at the sizes of the Iberian so that the SDM outputs are included
resolution<-res(bioclim_global_rn)
raster_template <- rast(ext(aggregated_vec),resolution=resolution)
aggregated_raster<-rasterize(aggregated_vec,raster_template,field=1,fun="count")
aggregated_raster[!is.na(aggregated_raster)] <- NA  # Set to 0 probability
crs(aggregated_raster)<- "EPSG:4326"

combined_rasters<-resample(combined_rasters,aggregated_raster)
crs(combined_rasters)<- "EPSG:4326"

combined_rasters_stack<-resample(combined_rasters_stack,aggregated_raster)
crs(combined_rasters_stack)<- "EPSG:4326"

#### Chapter 5 :: Execution of the aSDMs | Climate PRE-CONSTRAINED [for the watersheds where the species is present] ####


#INPUT OF ENVIRONMENTAL COVARIATES

#Unlock for Clima [1]
input_cov<-bioclim_global_rn

#Unlock for Hydro [2]
#input_cov<-riveratl_global

#Unlock for Hydroclimate [3]
#input_cov<-c(bioclim_global_rn,riveratl_global)



library(sdm)
#for (i in 1:5) {
for (i in 1:numbers_sp) {
  
  #Unlock for Hierarchical-Hydro [4]
  #clima_sdm<-rast("output/test/SDMs/Pre/Pre_Clima/preh12_clima_endemics_1km.tif")  
  #input_cov<-c(riveratl_global) 
  
  #Unlock for Hierarchical-Clima [5]
  # hydro_sdm<-rast("output/test/SDMs/Pre/Pre_Hydro/preh12_hydro_endemics_1km.tif")   
  # input_cov<-c(bioclim_global_rn)  
  
  
  # Get presence data for the current species
  species_data <- presence_long %>%
    filter(species_id == unique(presence_long$species_id)[[i]] & presence == 1)  # Filter for presence
  
  #Get the  watersheds ID (HYBAS_ID)
  watershed_id <- species_data$HYBAS_ID
  
  # Get the corresponding freshwater polygon for this HYBAS_ID
  masked_watershed<- watershed[watershed$HYBAS_ID %in% watershed_id, ]
  
  
  # Connect the watersheds produced in a single polygon
  polygon_list[[i]] <- masked_watershed
  
  #Dissolve the different  watersheds to function as the training region
  dissolve<-aggregate(polygon_list[[i]],dissolve=TRUE)
  
  #plot(dissolve)
  
  # Combine the two SpatRaster objects into one
  combined_raster <- c(input_cov)
  
  ############ RUN THE SDMs #############
  
  bioc_gal_in<-combined_raster
  
  crs(dissolve)<-"EPSG:4326"
  #crop the bioclim to the extent of the watersheds where the species belongs to
  bioc_gal_in_crop<-crop(bioc_gal_in,dissolve,mask=TRUE)
  #hydro_sdm_diss<-crop(hydro_sdm[[i]],dissolve,mask=TRUE)
  #plot(bioc_gal_in_crop[[1]])
  #plot(clima_sdm_diss)
  #plot(clima_sdm[[39]])
  
  #Insert the hierarchical feature herein, when case 4 and 5 takes place #HEREIN HIERARCHICAL
  # #First crop the added covariate too
  # hydro_sdm_res<-resample(hydro_sdm_diss,bioc_gal_in_crop,method="near")
  # bioc_gal_in_crop<-c(bioc_gal_in_crop,hydro_sdm_res)
  #bioc_gal_in_crop<-merge(bioc_gal_in_crop,clima_sdm[[i]])
  
  ############ RUN THE SDMs #############
  
  
  # Define the vector of numbers, gonna be needed to name variables later on 
  numbers_sp <- levels(as.factor(terra_df$species_id))
  
  # Pre-set the dataset
  terra_df_demo<-subset(terra_df,species_id==numbers_sp[i])
  terra_df_demo$species_id[terra_df_demo$species_id == numbers_sp[i]] <- 1 #sdms operate with presence data therefore "1" as label is essential
  terra_df_demo_f <- SpatialPointsDataFrame(terra_df_demo, data = terra_df_demo@data[, "species_id", drop = FALSE])
  
  
  #generate background points equal to 5% of the training area to preserve prevalence across scales
  background_points=sum(freq(bioc_gal_in_crop[[1]]))*0.05 #this is the right formula but in this case we will reduce the number for computational efficiency
  background_points= background_points*0.5
  #plot(bioc_gal_in_crop[[1]])
  
  #mapview(bioc_gal_in_crop[[1]])
  
  #Convert SpatRaster to RasterBrick
  r_bioc_gal_in_crop<-as(bioc_gal_in_crop, "Raster")
  
  d <- sdmData(species_id~., terra_df_demo_f, predictors= r_bioc_gal_in_crop, bg = list(method='gRandom',n=round(background_points),exclude=TRUE))
  
  #plot(d,cex = 0.3)
  
  #store the background points in a vector to be evaluated at a later stage
  sdm_d[[i]]<-d
  
  #   #sdm function to fit the models / d equals to the sdm presence-background data
  #   #since there are no indepence data we use replication method through sub or bootstrap
  #   #n=1 means 1 replications --> 1 models per method
  #   #parallel computing is specified too by using 4 cores
  
  m <- sdm(species_id~., d, methods=c('glm','brt','rf'), replication=c('boot'),
           test.p=30,n=1, parallelSetting=list(ncore=4,method='parallel'))
  
  #Current prediction to the rest of the 30%
  p2 <- predict(m, r_bioc_gal_in_crop)
  
  #ensemble accumulates the various models produced through SDMs - Use of a weighted AUC threshold independent index
  #en1 <- ensemble(m, p2,setting=list(method='weighted',stat='auc',opt=2))
  en1 <- ensemble(m, p2,setting=list(method='weighted',stat='auc'))
  #plot(en1)
  
  # Evaluate the ensemble model
  e <- evaluates(d, en1)
  
  #Provide the BoyceIndex
  bc<-sdm:::.boyce(e@observed,e@predicted)
  
  # Extract the evaluation metrics
  e_AUC <- e@statistics$AUC
  e_COR <- e@statistics$COR[1]  # Extract the first COR value
  CBI<-bc$CBI
  maxTSS<-e@threshold_based$TSS[2]
  maxKappa<-e@threshold_based$Kappa[5]
  t_maxSSS <- e@threshold_based$threshold[2]  # maxSSS threshold
  t_maxkappa <- e@threshold_based$threshold[5]  # maxkappa threshold
  t_prevalence <- e@threshold_based$threshold[10]  # Prevalence threshold
  CBI<-bc$CBI
  obs_prevalence<-length(d@species$species_id@presence) / (length(d@species$species_id@background) + length(d@species$species_id@presence) )
  
  #name the raster file based on the species name
  names(en1) <- paste("Species", i)
  
  # Append the results for this species to the results dataframe
  results_df <- rbind(results_df, data.frame(
    id=i,
    species_id = names(en1),  # Add the species ID
    e_AUC = e_AUC,
    e_COR = e_COR,
    CBI=CBI,
    maxTSS=maxTSS,
    maxKappa=maxKappa,
    t_maxSSS = t_maxSSS,
    t_maxkappa = t_maxkappa,
    t_prevalence = t_prevalence,
    obs_prevalence=obs_prevalence
  ))
  
  # ppend the raster to the list
  raster_list[[i]] <- en1
  
  # #save the super_iberian_rast to be used for the hierarchical applications
  #en1_hiera<-en1
  #if(i==1){
  #  super_iberian<-c(super_iberian,en1_hiera)
  #}else{
  #  en1_hierarchical_plus<-resample(en1,super_iberian,method="near")
  #  super_iberian<-c(super_iberian,en1_hierarchical_plus)
  #}
  
  #Combine the aSDM into the Iberian extent
  crop_en1<-mask(en1,aggregated_vec)
  #crop_en1
  #plot(crop_en1)
  mask_en1<- resample(crop_en1, combined_rasters, method = "near")
  #plot(mask_en1)
  
  # Reproject raster1 to the CRS of raster2
  crs(mask_en1) = "EPSG:4326"
  
  #Keep only the extent of the Iberian Peninsula
  #mask_en_1_ib<-merge(mask_en1,aggregated_raster)
  
  #mask_en_1_ib<-extent(mask_en_1_ib,ext(aggregated_raster))
  #plot(mask_en_1_ib)
  
  #store all the species in a combon raster file
  combined_rasters_stack <- c(combined_rasters_stack,mask_en1)
  
  #plot(combined_rasters_stack)
  
  #plot(mask_en_1_ib)
  
  # Define the file name for saving and NAME IT BASE ON THE SPECIES ACCORDING TO THE INTIAL LEDGER
  file_name <- paste0("deleteplease.tif")
  
  # Save the raster file
  writeRaster(en1, filename = file_name, overwrite = TRUE)
  
  # Print a message indicating the file has been saved
  cat("Saved:", file_name, "for iteration:", i, "\n")
}

#Assess the number of layers included in the uncontsrained raster layer
nlyr(combined_rasters_stack)

#Assess the metrics produced for each species
print(results_df)

#Save the unconstrained raster layer including the ensembled aSDMs
writeRaster(combined_rasters_stack, filename ="output/preh5_clima_endemics_1km.tif",overwrite=T)

#Save the super_iberian raster layer to be included as an added covariate in the case of hierarchical models
#writeRaster(super_iberian,filename ="output/test/SDMs/unconstrained_superiberia_hydro_endemics_1km.tif")

#Save the evaluation metrics and thresholding values
write.csv(results_df, "output/preh5_clima_endemics_metrics.csv", row.names = FALSE)

#species_names endemics table
species_name_vector<-results_df$species_id

plot(combined_rasters_stack)

###############################


#### Chapter 6 :: Threshold based on evaluation metrics ####

results_df<- read.csv( "output/preh5_clima_endemics_metrics.csv")

combined_rasters_stack<-rast("output/preh5_clima_endemics_1km.tif")

#Create an empty list to store the thresholded rasters
thresholded_rasters <- list()
thresholded_rasters_comb<-rast()
# Loop over each raster and apply the corresponding threshold
for (i in 1:max(results_df$id)) {
  # Get the threshold value from the data frame for the corresponding raster
  threshold_value <- results_df$obs_prevalence[i]
  
  # Apply the threshold using app() and store the result in the list
  thresholded_rasters[[i]] <- app(combined_rasters_stack[[i]], fun = function(x) ifelse(x > threshold_value, 1, 0))
  thresholded_rasters_comb<-c(thresholded_rasters_comb,thresholded_rasters[[i]])
  names(thresholded_rasters_comb)[i]<-results_df$species_id[i]
}

plot(thresholded_rasters_comb)


#Stack the maps | I was dealing with issues due to NAs

# Mask NA values in each 
#raster to 0 (or any other value you choose)
masked_rasters <- lapply(thresholded_rasters_comb, function(r) {
  # Mask out NA values (set them to 0)
  r[is.na(r)] <- 0
  return(r)
})

# Stack the masked rasters
stacked_raster <- rast(masked_rasters)

# Optionally, plot the stacked raster
#plot(stacked_raster)

#mask with the Iberian Peninsula

final_stack<-sum(stacked_raster)

plot(final_stack)

#mask it with the Iberian Peninsula network

#get a reference
reference<-thresholded_rasters_comb[[1]]

#get only the network
sp_richness<-crop(final_stack,reference,mask=TRUE)

plot(sp_richness)
