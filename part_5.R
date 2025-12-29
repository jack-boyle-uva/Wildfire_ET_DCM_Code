#Part 5

#Adding Precipitation
#and getting evaporation/Precipitation from reference hucs
library(soilDB)
library(kgc)
library(sf)
library(dplyr)
library(stringr)
library(ggplot2)
library(leaflet)
library(leaflet.extras)
library(purrr)
library(terra)
library(dataRetrieval)
library(lubridate)
library(budyko)
library(patchwork)
library(gridExtra)
library(cowplot)
library(tidyr)
library(nhdplusTools)
library(ncdf4)
library(shiny)
library(shinythemes)
library(rnaturalearth)
library(lubridate)

################################################################################
#Step #1: Convert our gpkg file into a shp

fires<-st_read("/Users/jackboyle/Downloads/UVA/Code/projects/Wildfire_ET_project/processed_data/#6_post_pre_fire_canopy/raster_method_results/DCM_analysis/DCM_Classified_LogView.gpkg")
st_write(fires, "/Users/jackboyle/Downloads/UVA/Code/projects/Wildfire_ET_project/processed_data/#6_post_pre_fire_canopy/raster_method_results/GGE_ready_data/raster_fire.shp")


################################################################################
################################################################################
#step #2 Below is the Google Earth Engine Code for getting ET data for all fires.
#Don't run this code here... it's designed to be ran on google earth engine.


            
            #// ** SCALABLE CODE: EXPORT ONE CSV PER GEOMETRY **
            # 
            #// === 1. Load the Entire HUC Collection ===
            #var huc_polygons = ee.FeatureCollection('projects/wildfire-water-quality/assets/CanopyLoss_HUCs_Delayed');
            
            #// Get the list of features so we can iterate client-side
            #var hucList = huc_polygons.toList(huc_polygons.size());
            #var hucCount = huc_polygons.size().getInfo();
            
            #print('Total HUCs found:', hucCount);
            #print('Starting export process (This will initiate', hucCount, 'tasks)...');
            
            #// Define the date filter range once
            #var START_DATE = '2000-01-01';
            #var END_DATE = '2024-12-31';
            
            #// --- Function to Process and Export a Single HUC ---
            #  var exportTimeSeriesForHUC = function(huc) {
            
            #    // Convert the list item back to a Feature
            #    var hucFeature = ee.Feature(huc);
            #    var geom = hucFeature.geometry();
            
            #    // Get the unique ID as a server-side object
            #    var evntIdObj = hucFeature.get('Evnt_ID');
            
            #    // Get the unique ID as a client-side string for use in exports/selectors/prints
            #    var evntIdStr = evntIdObj.getInfo(); 
            
            #    // === 2. Load and filter SSEBop dataset ===
            #      var ssebop = ee.ImageCollection('OpenET/SSEBOP/CONUS/GRIDMET/MONTHLY/v2_0')
            #    .filterBounds(geom)
            #    .filterDate(START_DATE, END_DATE)
            #    .select('et');
            
            #    // === 3. Compute mean ET and set the column name to Evnt_ID ===
            #      var etSeries = ssebop.map(function(img) {
            
            #        // Calculate the mean ET (using the defensive .get() with -9999 as default)
            #        var regionMean = img.reduceRegion({
            #          reducer: ee.Reducer.mean(),
            #          geometry: geom,
            #          scale: 1000,
            #          bestEffort: true
            #        });
            
            #        // CRITICAL FIX: Use .get('et') and provide the default value as an ee.Number.
            #        // This is the most reliable way to prevent Dictionary.set crashes.
            #        var meanET = regionMean.get('et', ee.Number(-9999));
            
            #        // Create a base feature with the date property.
            #        var feature = ee.Feature(null, {
            #          'date': img.date().format('YYYY-MM')
            #        });
            
            #        // Use .set() to assign the property using the client-side string key
            #        return feature.set(evntIdStr, meanET);
            #      });
            
            #    // === 4. Export to Google Drive (Clean Output) ===
            #      Export.table.toDrive({
            #        collection: etSeries,
            #        description: 'ET_TS_' + evntIdStr,
            #        folder: 'Earth_Engine_HUC_Exports',
            #        fileFormat: 'CSV',
            #        // Selectors must include the date and the dynamic ID column name
            #        selectors: ['date', evntIdStr] 
            #      });
            #  };
            
            #// --- Execute the client-side loop ---
            #  // Iterate over the client-side list of HUC features and call the export function for each.
            #for (var i = 0; i < hucCount; i++) {
            #  var currentHuc = hucList.get(i);
            #  exportTimeSeriesForHUC(currentHuc);
            #}
            
            #// NOTE: Check the 'Tasks' tab (right side of the GEE console) to monitor the status
            #// of the 200 individual export tasks that were initiated.




#METHOD #2 ~ Raster ~ Get only the ET data from the pixels that have been burned (the AI tells me that it's possible!!)

################################################################################
#Step 2: Get Percipitation data for every fire area
        # Use Prism
        #Get monthy data from 1998-2024
                                                    # for some reason the code got AK fires???
                
                # 1. Load necessary packages (make sure stringr is installed)
                # install.packages("stringr")
                library(sf)
                library(terra)
                library(stringr) # Added for robust filename parsing
                
                # 2. Set your file paths
                wildfire_gpkg_path <- "/Users/jackboyle/Downloads/UVA/Code/projects/Wildfire_ET_project/processed_data/#6_post_pre_fire_canopy/raster_method_results/raster_method_pre_post_canopy_fire.gpkg"
                prism_data_directory <- "/Users/jackboyle/Downloads/UVA/Database/PRISM/PRISM_ppt_stable_4kmM3_198101_202411_bil"
                output_csv_directory <- "/Users/jackboyle/Downloads/UVA/Code/projects/Wildfire_ET_project/processed_data/#8_precipitation_data"
                
                # Create the output directory if it doesn't exist
                if (!dir.exists(output_csv_directory)) {
                  dir.create(output_csv_directory, recursive = TRUE)
                }
                
                # 3. Load the wildfire data
                wildfires <- st_read(wildfire_gpkg_path)
                cat("Loaded", nrow(wildfires), "wildfires from the GeoPackage.\n")
                
                # 4. Create a robust list of the PRISM precipitation files to process
                # The pattern is slightly more general to catch variations like 'stable', 'provisional', etc.
                prism_files <- list.files(path = prism_data_directory,
                                          pattern = "_bil\\.bil$", 
                                          full.names = TRUE)
                
                # --- KEY CHANGE: Use regular expressions to find the date code ---
                # str_extract will find the first sequence of 6 consecutive digits (YYYYMM) in each filename.
                # This is not dependent on file name length.
                date_codes <- str_extract(basename(prism_files), "\\d{6}")
                
                # Now, filter files based on the extracted date codes
                # This keeps only files where a valid date was found and the year is 1998 or later.
                valid_indices <- which(as.integer(substr(date_codes, 1, 4)) >= 1998 & !is.na(date_codes))
                prism_files <- prism_files[valid_indices]
                date_codes <- date_codes[valid_indices]
                
                # Sort the files chronologically based on the extracted date codes
                # This is critical to ensure the time series is in the correct order.
                order_of_files <- order(date_codes)
                prism_files <- prism_files[order_of_files]
                date_codes <- date_codes[order_of_files]
                
                cat("Found", length(prism_files), "PRISM files to process from 1998 onwards.\n\n")
                
                # 5. --- Start of the main loop to iterate through each wildfire ---
                for (i in 1:nrow(wildfires)) {
                  
                  current_wildfire <- wildfires[i, ]
                  wildfire_id <- current_wildfire$Event_ID
                  
                  cat("Processing wildfire", i, "of", nrow(wildfires), ":", wildfire_id, "\n")
                  
                  precipitation_data <- data.frame(Date = as.Date(character()),
                                                   Precipitation_mm = numeric())
                  
                  # 6. Loop through each PRISM file for the current wildfire
                  # We loop using an index to easily access the corresponding date_code
                  for (j in 1:length(prism_files)) {
                    prism_file <- prism_files[j]
                    current_date_code <- date_codes[j]
                    
                    # --- THIS IS THE ROBUST DATE CREATION ---
                    year <- as.integer(substr(current_date_code, 1, 4))
                    month <- as.integer(substr(current_date_code, 5, 6))
                    date <- as.Date(paste(year, month, "01", sep = "-"))
                    
                    prism_raster <- rast(prism_file)
                    wildfire_transformed <- st_transform(current_wildfire, crs(prism_raster))
                    mean_precipitation <- terra::extract(prism_raster, wildfire_transformed, fun = mean, na.rm = TRUE)
                    
                    precipitation_data <- rbind(precipitation_data,
                                                data.frame(Date = date,
                                                           Precipitation_mm = mean_precipitation[1, 2])) # Use [1, 2] for clarity
                  }
                  
                  # 7. Save the results for the current wildfire to a unique CSV file
                  # Using file.path() is a safer way to build file paths
                  output_filename <- file.path(output_csv_directory, paste0("precipitation_", wildfire_id, ".csv"))
                  write.csv(precipitation_data, output_filename, row.names = FALSE)
                  
                  cat("  -> Successfully created CSV file:", basename(output_filename), "\n\n")
                }
                
                cat("All wildfires have been processed.\n")
                
################################################################################
################################################################################
################################################################################
#Step 2: Filter out reference huc dataset (before we filter the fires a bit more and we can't forget to do it to the unbunred huc)
        #Also, filter so that there are max 3 hucs in every wildfire.        

                library(sf)
                library(dplyr)
                
                path_dcm_classified <- "/Users/jackboyle/Downloads/UVA/Code/projects/Wildfire_ET_project/processed_data/#6_post_pre_fire_canopy/raster_method_results/DCM_analysis/DCM_Classified_LogView.gpkg"
                path_unburned_raw <- "/Users/jackboyle/Downloads/UVA/Code/projects/Wildfire_ET_project/processed_data/#3_Unburned_huc_shape/unburned_huc_polygons.gpkg"
                output_dir <- "/Users/jackboyle/Downloads/UVA/Code/projects/Wildfire_ET_project/processed_data/#8_precipitation_data"

                if (!dir.exists(output_dir)) {
                  dir.create(output_dir, recursive = TRUE)
                }
                output_file <- file.path(output_dir, "Filtered_Unburned_HUCs_Max3.gpkg")
                dcm_data <- st_read(path_dcm_classified, quiet = TRUE)
                target_ids <- unique(dcm_data$Event_ID)
                unburned_data <- st_read(path_unburned_raw, quiet = TRUE)
                filtered_hucs <- unburned_data %>%
                  # Filter 1: Keep only rows where Event_ID exists in the DCM file
                  filter(Event_ID %in% target_ids) %>%
                  # Filter 2: Remove exact duplicates (safety check)
                  distinct(Event_ID, huc12, .keep_all = TRUE) %>%
                  # Filter 3: Keep Max 3 HUCs per Event_ID
                  group_by(Event_ID) %>%
                  slice_head(n = 3) %>%
                  ungroup()
                st_write(filtered_hucs, output_file, delete_dsn = TRUE, quiet = TRUE)
                
################################################################################
################################################################################
################################################################################
#Step 3: Precipitation data every reference huc 
                library(sf)
                library(terra)
                library(stringr)
                library(dplyr) 
                
                # ==============================================================================
                # 1. SETUP & LOAD
                # ==============================================================================
                unburned_gpkg_path <- "/Users/jackboyle/Downloads/UVA/Code/projects/Wildfire_ET_project/processed_data/#8_precipitation_data/Filtered_Unburned_HUCs_Max3.gpkg"
                prism_data_directory <- "/Users/jackboyle/Downloads/UVA/Database/PRISM/PRISM_ppt_stable_4kmM3_198101_202411_bil"
                output_csv_directory <- "/Users/jackboyle/Downloads/UVA/Code/projects/Wildfire_ET_project/processed_data/#8_precipitation_data/Unburned_huc_P_ET/P_unburned"
                
                if (!dir.exists(output_csv_directory)) dir.create(output_csv_directory, recursive = TRUE)
                
                # Load the data
                unburned_polys <- st_read(unburned_gpkg_path)
                cat("\n--- RAW DATA LOADED ---\n")
                cat("Original Row Count:", nrow(unburned_polys), "(Should be ~18,024)\n")
                
                # ==============================================================================
                # 2. THE FILTER (CRITICAL STEP)
                # ==============================================================================
                cat("Applying 'Max 3 HUCs per Fire' filter...\n")
                
                # We create a NEW variable name to be absolutely sure
                filtered_polys <- unburned_polys %>%
                  distinct(Event_ID, huc12, .keep_all = TRUE) %>% 
                  group_by(Event_ID) %>%
                  slice_head(n = 3) %>% 
                  ungroup()
                
                cat("Filtered Row Count:", nrow(filtered_polys), "\n")
                
                # ==============================================================================
                # 3. THE SAFETY CHECK (Panic Button)
                # ==============================================================================
                if (nrow(filtered_polys) > 5000) {
                  stop("⛔ STOP! The filter did not work. You still have ", nrow(filtered_polys), " rows. Do not proceed.")
                } else {
                  cat("✅ SAFETY CHECK PASSED. We successfully reduced the dataset to", nrow(filtered_polys), "rows.\n")
                }
                
                # ==============================================================================
                # 4. PREPARE PRISM FILES
                # ==============================================================================
                prism_files <- list.files(path = prism_data_directory, pattern = "_bil\\.bil$", full.names = TRUE)
                date_codes <- str_extract(basename(prism_files), "\\d{6}")
                valid_indices <- which(as.integer(substr(date_codes, 1, 4)) >= 1998 & !is.na(date_codes))
                
                prism_files <- prism_files[valid_indices][order(date_codes[valid_indices])]
                date_codes <- date_codes[valid_indices][order(date_codes[valid_indices])]
                
                # ==============================================================================
                # 5. TRANSFORM & RUN (The "Individual Loop" you like)
                # ==============================================================================
                temp_rast <- rast(prism_files[1])
                cat("\nAligning map projections (this happens ONCE)...\n")
                # Note: We use 'filtered_polys' here, NOT 'unburned_polys'
                polys_ready <- st_transform(filtered_polys, crs(temp_rast)) 
                
                cat("Starting the loop. Look at the progress below:\n\n")
                
                for (i in 1:nrow(polys_ready)) {
                  
                  # Extract single feature
                  current_poly <- polys_ready[i, ]
                  wildfire_id <- current_poly$Event_ID
                  huc_id <- current_poly$huc12
                  
                  # Print the filtered progress
                  cat("Processing", i, "of", nrow(polys_ready), "| Fire:", wildfire_id, "\n")
                  
                  # Setup storage
                  precipitation_data <- data.frame(Date = as.Date(character()), Precipitation_mm = numeric())
                  
                  # Loop through weather files
                  for (j in 1:length(prism_files)) {
                    
                    current_date_code <- date_codes[j]
                    year <- as.integer(substr(current_date_code, 1, 4))
                    month <- as.integer(substr(current_date_code, 5, 6))
                    date_obj <- as.Date(paste(year, month, "01", sep = "-"))
                    
                    prism_raster <- rast(prism_files[j])
                    
                    mean_precipitation <- terra::extract(prism_raster, current_poly, fun = mean, na.rm = TRUE)
                    val <- mean_precipitation[1, 2]
                    
                    precipitation_data <- rbind(precipitation_data, data.frame(Date = date_obj, Precipitation_mm = val))
                  }
                  
                  output_filename <- file.path(output_csv_directory, paste0("precip_", wildfire_id, "_", huc_id, ".csv"))
                  write.csv(precipitation_data, output_filename, row.names = FALSE)
                }
                
                cat("\nDone! Processed", nrow(polys_ready), "HUCs.\n")
################################################################################
################################################################################
################################################################################
#Step 4: ET data every reference huc (Use the GEE code from Part 4) 
          
    #The code below is what you would copy into Google Earth Engine to find ET:
          
                
                
         "       // =============================================================================
                  // PART 1: DIAGNOSTIC CHECK
                // =============================================================================
                  
                  // 1. Point to your uploaded asset
                // *** ACTION REQUIRED: Check this path! ***
                  // It might be "users/jackboyle/Filtered_Unburned_HUCs_gee" 
                // or "projects/wildfire-water-quality/assets/Filtered_Unburned_HUCs_gee"
                var hucs = ee.FeatureCollection('projects/wildfire-water-quality/assets/Filtered_Unburned_HUCs_gee');
                
                // 2. Print the first item to the Console
                print("--- ASSET CHECK ---");
                print("Number of HUCs:", hucs.size());
                print("First Feature Columns:", hucs.first());
                
                // INSTRUCTIONS:
                  // Look at the Console on the right. 
                // Open 'properties'. 
                // Did 'Event_ID' turn into 'Evnt_ID'? 
                  // Did 'huc12' turn into 'huc12'?
                  // If names are different, update the "Part 2" code below to match.
                // =============================================================================
                  // PART 2: UP-TO-DATE BATCH SCRIPT FOR UNBURNED HUCs (30m Resolution)
                // =============================================================================
                  
                  // 1. LOAD YOUR ASSET
                // *** ACTION REQUIRED: Ensure this path matches the one that worked in Part 1 ***
                  var hucs = ee.FeatureCollection('projects/wildfire-water-quality/assets/Filtered_Unburned_HUCs_gee');
                
                // 2. SETUP DATASETS
                // OpenET SSEBop (Evapotranspiration) - Monthly
                var et_collection = ee.ImageCollection("OpenET/SSEBOP/CONUS/GRIDMET/MONTHLY/v2_0");
                
                // BATCH SETTINGS: 50 HUCs per CSV is safe for 30m resolution
                var batchSize = 50; 
                
                // 3. THE PROCESSING LOGIC
                var processHUC = function(feature) {
                  
                  // *** ACTION REQUIRED: Check column names based on Part 1 results ***
                    // Shapefiles often shorten names (e.g., 'Event_ID' -> 'Evnt_ID')
                  var fire_id = feature.get('Event_ID'); 
                  var huc_id = feature.get('huc12');
                  var geom = feature.geometry();
                  
                  // Filter ET for this specific HUC
                  // We do not filter by date, so it grabs 2000-Present automatically
                  var et_loc = et_collection.filterBounds(geom).select('et');
                  
                  return et_loc.map(function(img) {
                    // Reduce region at 30m resolution
                    var stats = img.reduceRegion({
                      reducer: ee.Reducer.mean(),
                      geometry: geom,
                      scale: 30,         // High Resolution (Matches your Fire Analysis)
                      tileScale: 16,     // CRITICAL: Prevents "Memory Limit Exceeded" crashes
                      maxPixels: 1e13,
                      bestEffort: true
                    });
                    
                    return ee.Feature(null, {
                      'Event_ID': fire_id,
                      'huc12': huc_id,
                      'Date': img.date().format('YYYY-MM-dd'),
                      'Year': img.date().get('year'),
                      'Month': img.date().get('month'),
                      'ET_mm': stats.get('et')
                    });
                  });
                };
                
                // 4. BATCH EXECUTION LOOP (Client-Side)
                // This splits your list into chunks of 50 and creates a separate task for each.
                var hucList = hucs.toList(hucs.size());
                var totalHucs = hucs.size().getInfo();
                
                print("Total HUCs to process:", totalHucs);
                print("Creating batches of:", batchSize);
                
                for (var i = 0; i < totalHucs; i += batchSize) {
                  
                  // Create the slice for this batch (e.g., 0-50, 50-100)
                  var batch = ee.FeatureCollection(hucList.slice(i, i + batchSize));
                  
                  // Run the math and flatten the collection
                  var processed = batch.map(processHUC, true).flatten();
                  
                  // Create ONE export task for this batch
                  Export.table.toDrive({
                    collection: processed,
                    description: 'Unburned_HUC_ET_Batch_' + i,
                    folder: 'Unburned_HUC_ET_Project', // Check this folder in Drive later
                    fileFormat: 'CSV',
                    selectors: ['Event_ID', 'huc12', 'Date', 'Year', 'Month', 'ET_mm']
                  });
                }   "
                
                
################################################################################
 #The new code:               
                
              "  // =============================================================================
                  // 1. SETUP & INPUTS
                // =============================================================================
                  
                  // *** ACTION REQUIRED: Change this to your uploaded asset path ***
                  var fires = ee.FeatureCollection("users/jackboyle/Gold_Standard_Fires");
                
                // MTBS Data (Burn Severity) - Used to mask out unburned islands
                var mtbs_collection = ee.ImageCollection("USFS/GTAC/MTBS/annual_burn_severity_mosaics/v1");
                
                // OpenET SSEBop (Evapotranspiration) - 30m Native Resolution
                var et_collection = ee.ImageCollection("OpenET/SSEBOP/CONUS/GRIDMET/MONTHLY/v2_0");
                
                // Define Time Window (How many years before/after to extract?)
                // Example: 1 year pre, 5 years post.
                var years_buffer = 6; 
                
                // =============================================================================
                  // 2. THE CORE LOGIC FUNCTION
                // =============================================================================
                  var processFire = function(feature) {
                    
                    // A. GET FIRE INFO
                    var fire_id = feature.get('Event_ID'); // Ensure column name matches your upload
                    var fire_date = ee.Date(feature.get('Ig_Date')); 
                    var fire_year = fire_date.get('year');
                    var geom = feature.geometry();
                    
                    // B. GET THE BURN MASK (The "Raster Method" in Cloud)
                    // Find the MTBS mosaic for the year of the fire
                    var mtbs_year = mtbs_collection
                    .filterDate(
                      ee.Date.fromYMD(fire_year, 1, 1), 
                      ee.Date.fromYMD(fire_year, 12, 31)
                    ).first();
                    
                    // If MTBS is missing for that year, return null (Safety check)
                    // (We use ee.Algorithms.If to handle nulls safely in mapping)
                    return ee.Algorithms.If(mtbs_year, computeET(mtbs_year, fire_date, geom, fire_id), null);
                  };
                
                // --- Sub-Function to actually compute ET if MTBS exists ---
                  var computeET = function(mtbs_img, fire_date, geom, fire_id) {
                    
                    // 1. Create the Mask
                    // MTBS Classes: 1=Unburned, 2=Low, 3=Mod, 4=High
                    // We keep only 2, 3, 4 (The Burned Area)
                    // 'Severity' band usually holds the class. Check collection docs if needed.
                    var burn_mask = mtbs_img.select('Severity').gte(2).and(mtbs_img.select('Severity').lte(4));
                    
                    // Clip the mask to the fire polygon to save processing power
                    var local_mask = burn_mask.clip(geom);
                    
                    // 2. Define Time Range for this specific fire
                    // e.g. Start 1 year before, End 5 years after
                    var start_date = fire_date.advance(-1, 'year');
                    var end_date = fire_date.advance(5, 'year');
                    
                    var et_window = et_collection
                    .filterBounds(geom)
                    .filterDate(start_date, end_date)
                    .select('et');
                    
                    // 3. Map over every month of ET data
                    var time_series = et_window.map(function(et_img) {
                      
                      // *** THE GOLD STANDARD STEP ***
                        // Mask the ET image using the MTBS Burn Mask
                      // Pixels that are "Unburned Islands" become NULL and are excluded from the mean.
                      var et_masked = et_img.updateMask(local_mask);
                      
                      // Calculate Mean
                      var stats = et_masked.reduceRegion({
                        reducer: ee.Reducer.mean(),
                        geometry: geom,
                        scale: 30,  // <--- CRITICAL: Use 30m native resolution, not 1000m
                        maxPixels: 1e9
                      });
                      
                      return ee.Feature(null, {
                        'Event_ID': fire_id,
                        'Date': et_img.date().format('YYYY-MM-dd'),
                        'ET_Mean_mm': stats.get('et'),
                        'Month': et_img.date().get('month'),
                        'Year': et_img.date().get('year')
                      });
                    });
                    
                    return time_series;
                  };
                
                // =============================================================================
                  // 3. EXECUTE AND FLATTEN
                // =============================================================================
                  
                  // Run the function on all fires
                var nested_results = fires.map(processFire, true); // 'true' drops nulls
                
                // Flatten: Converts a "Collection of Collections" into one big Table
                var flat_results = nested_results.flatten();
                
                // =============================================================================
                  // 4. EXPORT (ONE TASK)
                // =============================================================================
                  
                  Export.table.toDrive({
                    collection: flat_results,
                    description: 'Gold_Standard_ET_RasterMethod',
                    folder: 'Wildfire_ET_Project', // Change to your Drive folder name
                    fileFormat: 'CSV',
                    selectors: ['Event_ID', 'Date', 'Year', 'Month', 'ET_Mean_mm']
                  });
                
                print("Export setup complete. Check Tasks tab.");"
################################################################################
################################################################################
#Step 5: Combining the data                
                # ==============================================================================
                # 0. LIBRARIES
                # ==============================================================================
                library(tidyverse)
                library(lubridate)
                library(stringr)
                library(sf)
                
                # ==============================================================================
                # 1. SETUP PATHS
                # ==============================================================================
                # Inputs
                path_fire_et_file   <- "/Users/jackboyle/Downloads/UVA/Code/projects/Wildfire_ET_project/processed_data/#7_ET_data/ET_Fire/ET_for_all_fires.csv"
                path_fire_p_folder  <- "/Users/jackboyle/Downloads/UVA/Code/projects/Wildfire_ET_project/processed_data/#8_precipitation_data/Burned_P"
                
                path_ref_et_folder  <- "/Users/jackboyle/Downloads/UVA/Code/projects/Wildfire_ET_project/processed_data/#7_ET_data/ET_unburned/New_data"
                path_ref_p_folder   <- "/Users/jackboyle/Downloads/UVA/Code/projects/Wildfire_ET_project/processed_data/#8_precipitation_data/Unburned_P/P_unburned"
                
                path_metadata_gpkg  <- "/Users/jackboyle/Downloads/UVA/Code/projects/Wildfire_ET_project/processed_data/#6_post_pre_fire_canopy/raster_method_results/DCM_analysis/DCM_Classified_LogView.gpkg"
                
                # Output
                output_dir          <- "/Users/jackboyle/Downloads/UVA/Code/projects/Wildfire_ET_project/processed_data/#9_P_ET_data"
                
                # ==============================================================================
                # 2. LOAD METADATA (Ignition Dates)
                # ==============================================================================
                message("Step 1: Loading Ignition Dates...")
                
                meta_sf <- st_read(path_metadata_gpkg, quiet = TRUE)
                metadata <- meta_sf %>%
                  st_drop_geometry() %>%
                  select(Event_ID, Ig_Date) %>%
                  mutate(Ig_Date = as.Date(Ig_Date)) %>%
                  distinct(Event_ID, .keep_all = TRUE)
                
                # ==============================================================================
                # 3. LOAD FIRE DATA (ET & P)
                # ==============================================================================
                message("Step 2: Processing FIRE Data...")
                
                # --- A. Load Fire ET (Single File) ---
                # We read assuming headers exist, then clean up
                fire_et_raw <- read_csv(path_fire_et_file, show_col_types = FALSE)
                
                # SAFETY: If row 1 looks like data (e.g., starts with AL333...), reload without headers
                if (str_starts(names(fire_et_raw)[1], "AL|TX|MT|CA|NV|WA")) {
                  fire_et_raw <- read_csv(path_fire_et_file, col_names = FALSE, show_col_types = FALSE)
                  colnames(fire_et_raw) <- c("Event_ID", "Date", "Year", "Month", "ET_Burned_mm", "ET_Whole_mm")
                } else {
                  # If headers exist, ensure consistent naming based on your screenshot
                  # Col 1=ID, Col 2=Date, Col 6=ET_Whole_mm
                  colnames(fire_et_raw)[1] <- "Event_ID"
                  colnames(fire_et_raw)[2] <- "Date"
                  colnames(fire_et_raw)[6] <- "ET_Whole_mm"
                }
                
                # CLEANING & FILTERING
                fire_et <- fire_et_raw %>%
                  select(Event_ID, Date, ET_Fire_mm = ET_Whole_mm) %>% # Select ONLY Whole_mm
                  mutate(
                    Date = ymd(Date), 
                    ET_Fire_mm = as.numeric(ET_Fire_mm)
                  ) %>%
                  filter(!is.na(ET_Fire_mm)) # <--- FILTER: Removes rows with empty data
                
                message(paste("  -> Fire ET Rows after filtering:", nrow(fire_et)))
                
                # --- B. Load Fire P (Many Files) ---
                message("  -> Reading Fire Precipitation files...")
                fire_p_files <- list.files(path_fire_p_folder, pattern = "\\.csv$", full.names = TRUE)
                
                fire_p <- map_dfr(fire_p_files, function(f) {
                  fname <- basename(f)
                  id <- str_remove(str_remove(fname, "precipitation_"), "\\.csv")
                  
                  df <- read_csv(f, show_col_types = FALSE)
                  colnames(df) <- c("Date", "P_Burned_mm")
                  df %>% mutate(Event_ID = id, Date = ymd(Date))
                })
                
                # --- C. Merge Fire ET & P (Left Join) ---
                # We KEEP fire rows even if P is missing
                fire_master <- fire_et %>%
                  left_join(fire_p, by = c("Event_ID", "Date"))
                
                # ==============================================================================
                # 4. LOAD REFERENCE DATA (ET & P)
                # ==============================================================================
                message("Step 3: Processing REFERENCE Data...")
                
                # --- A. Load Reference ET (Batches) ---
                message("  -> Reading Reference ET batches...")
                ref_et_files <- list.files(path_ref_et_folder, pattern = "\\.csv$", full.names = TRUE)
                
                ref_et <- map_dfr(ref_et_files, function(f) {
                  read_csv(f, show_col_types = FALSE) %>%
                    select(Event_ID, huc12, Date, ET_Ref_mm = ET_mm) %>%
                    mutate(Date = ymd(Date), huc12 = as.character(huc12))
                })
                
                # --- B. Load Reference P (Individual Files) ---
                message("  -> Reading Reference Precipitation files...")
                ref_p_files <- list.files(path_ref_p_folder, pattern = "\\.csv$", full.names = TRUE)
                
                ref_p <- map_dfr(ref_p_files, function(f) {
                  fname <- basename(f)
                  clean_name <- str_remove(str_remove(fname, "precip_"), "\\.csv")
                  parts <- str_split(clean_name, "_")[[1]]
                  
                  if(length(parts) < 2) return(NULL)
                  
                  fire_id <- parts[1]
                  huc_id  <- parts[2]
                  
                  df <- read_csv(f, show_col_types = FALSE)
                  colnames(df) <- c("Date", "P_Ref_mm")
                  
                  df %>% mutate(Event_ID = fire_id, huc12 = as.character(huc_id), Date = ymd(Date))
                })
                
                # --- C. Merge Reference ET & P (Inner Join) ---
                # We DROP reference rows if P is missing
                ref_master <- ref_et %>%
                  inner_join(ref_p, by = c("Event_ID", "huc12", "Date"))
                
                message(paste("  -> Combined Reference Data. Rows:", nrow(ref_master)))
                
                # ==============================================================================
                # 5. FINAL MERGE & CALCULATIONS
                # ==============================================================================
                message("Step 4: Creating Final Dataset...")
                
                # 1. Join Reference Master with Fire Master
                # This creates the Fire vs Reference comparison
                full_df <- ref_master %>%
                  inner_join(fire_master, by = c("Event_ID", "Date"))
                
                # 2. Add Ignition Dates
                full_df <- full_df %>%
                  left_join(metadata, by = "Event_ID")
                
                # 3. Calculate Timelines
                full_df <- full_df %>%
                  filter(!is.na(Ig_Date)) %>%
                  mutate(
                    Fire_Year = year(Ig_Date),
                    Current_Year = year(Date),
                    Month = month(Date),
                    
                    # Years since fire (e.g. -1, 0, 1)
                    Years_Since_Fire = Current_Year - Fire_Year,
                    
                    # Pre/Post classification
                    Period = ifelse(Date < Ig_Date, "Pre-Fire", "Post-Fire")
                  ) %>%
                  select(Event_ID, Ref_HUC12 = huc12, Date, Year=Current_Year, Month, 
                         Ig_Date, Years_Since_Fire, Period,
                         ET_Fire_mm, P_Burned_mm,  # Note: renamed to ET_Fire_mm (from Whole)
                         ET_Ref_mm, P_Ref_mm)
                
                # ==============================================================================
                # 6. EXPORT
                # ==============================================================================
                if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
                output_file <- file.path(output_dir, "ET_P_full_dataset.csv")
                
                message(paste("Step 5: Saving", nrow(full_df), "rows to CSV..."))
                write_csv(full_df, output_file)
                
                message("✅ DONE! File saved: ", output_file)
                
                
                
                
                
                
                
                
                
################################################################################                
################################################################################                
################################################################################  
#PART 5 OUTPUT: 
                #Full ET and P dataset for the wildfire zone and the reference huc zone
                
                #output<-read.csv("/Users/jackboyle/Downloads/UVA/Code/projects/Wildfire_ET_project/processed_data/#9_P_ET_data/ET_P_full_dataset.csv")
                
                #End of Part 5
                #Next up... part 6 full analysis...story: 
                  "- Show that ET post wildfire is not infulenced by P but by the burn itself"
                  "- find the magnitude of ET loss for each DCM classified fire"
                  "- Determine the difference inrecovery rates of ET for fires with High Moderate, and Immedate DCM"
                
