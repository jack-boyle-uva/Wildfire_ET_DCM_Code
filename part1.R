# PART 1
#Wildfire identification


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
library(StreamCatTools)
#Step 1: Load MTBS Perimeter
#Step 2: Filter Fires from 2000-2024
#Step 3: Filter small fires, only ones larger than 14km^2
#Step 4: Add Climate Zone for each fire
#Step 5: Add soil type details for each fire?
#Step 6: Save file as


################################################################################
#Step 1: Load MTBS
      Fire_boundary <- read_sf(
        paste0("/Users/jackboyle/Downloads/UVA/Database/WILDFIRE_WATERQUALITY",
               "/MTBS/mtbs_perimeter_data/mtbs_perims_DD.shp")
      )
################################################################################
#Step 2: Filter Fires from 2000-2024
        Fire_boundary <- Fire_boundary %>%
          mutate(Ig_Date = as.Date(Ig_Date)) %>%        
          filter(Ig_Date >= as.Date("2000-01-01"))

################################################################################
#Step 3: Filter small fires, only ones larger than 5km^2 = Cr1235.53 Acres
      Fire_boundary <- Fire_boundary %>%
        filter(
          BurnBndAc >= 1235.53                             
        )
################################################################################     
#Step 4: Filter out Prescribed Fires
      Fire_boundary<-Fire_boundary %>%
        filter(
          Incid_Type != "Prescribed Fire"
        )
################################################################################  
#Step 5: Filter out AK, HI, PR
      Fire_boundary <- Fire_boundary %>%
        filter(!grepl("^(AK|HI|PR)", Event_ID))
################################################################################   
#Step 6: Clean data
      cols_to_remove <- c("irwinID", "Map_ID", "Map_Prog", "Asmnt_Type", "Pre_ID", "Post_ID", "Perim_ID",
                          "BurnBndAc", "BurnBndLat", "BurnBndLon", "dNBR_offst",
                          "dNBR_stdDv", "NoData_T", "IncGreen_T", "Low_T",
                          "Mod_T", "High_T", "Comment")
      
      Fire_boundary_clean <- dplyr::select(Fire_boundary, -any_of(cols_to_remove))
      names(Fire_boundary_clean)
################################################################################
#Step 7: Add Climate Zone for each fire
      koppen<-rast("/Users/jackboyle/Downloads/UVA/Database/Climate_zones/koppen_geiger_tif/1991_2020/koppen_geiger_0p1.tif")
      
      fires <- Fire_boundary_clean
      fires_vect <- vect(fires)
      climate_extract <- terra::extract(
        koppen,
        fires_vect,
        fun = NULL,   
        na.rm = TRUE
      )
      climate_extract$Event_ID <- fires$Event_ID[climate_extract$ID]
      summary <- climate_extract %>%
        group_by(Event_ID) %>%
        dplyr::summarise(
          koppen_codes = paste(sort(unique(koppen_geiger_0p1)), collapse = ","),
          .groups = "drop"
        )
      Fire_boundary_clean <- fires %>%
        left_join(summary, by = "Event_ID")
      Fire_boundary_clean <- Fire_boundary_clean %>%
        mutate(
          koppen_id = sapply(
            strsplit(koppen_codes, ","),         # split "7,9" → c("7","9")
            function(x) {
              # convert to numeric and map each to its character code
              paste(sapply(as.numeric(x), kgc::getZone), collapse = ",")
            }
          )
        )
      

################################################################################
#Step 8: Adding SteamCat data


      library(FedData)
      library(terra)
      library(dplyr)
      library(sf)
      
      # ====================================================
      # 1. SETUP
      # ====================================================
      # Create directories
      dir.create("./FedData_Raw", showWarnings = FALSE)
      dir.create("./FedData_Raw/NED", showWarnings = FALSE)
      dir.create("./FedData_Raw/NLCD", showWarnings = FALSE)
      
      # Input Data (Assumes AK/HI/PR already removed)
      fires_to_process <- Fire_boundary_clean 
      
      print(paste("Fires to Process:", nrow(fires_to_process)))
      
      # Prepare storage
      results_list <- list()
      output_file <- "feddata_results_partial.csv"
      
      # ====================================================
      # 2. THE LOOP
      # ====================================================
      
      for (i in 1:nrow(fires_to_process)) {
        
        single_fire <- fires_to_process[i, ]
        fire_id <- single_fire$Event_ID
        
        # Print progress
        message(paste("Processing", i, "of", nrow(fires_to_process), "::", fire_id))
        
        tryCatch({
          
          # ----------------------------------------------------
          # 🏔️ PART A: ELEVATION (NED) - (This was already working)
          # ----------------------------------------------------
          dem <- get_ned(
            template = single_fire,
            label = fire_id,
            extraction.dir = "./FedData_Raw/NED",
            force.redo = FALSE
          )
          
          # Masking
          fire_proj <- st_transform(single_fire, crs(dem))
          dem_masked <- crop(dem, fire_proj)
          dem_masked <- mask(dem_masked, fire_proj)
          
          # Stats
          slope <- terra::terrain(dem_masked, v = "slope", unit = "degrees")
          elev_stats <- global(dem_masked, c("mean", "sd"), na.rm = TRUE)
          slope_stats <- global(slope, "mean", na.rm = TRUE)
          
          # ----------------------------------------------------
          # 🌲 PART B: LAND COVER (NLCD) - (UPDATED FIX)
          # ----------------------------------------------------
          nlcd <- get_nlcd(
            template = single_fire,
            label = fire_id,
            year = 2019, 
            extraction.dir = "./FedData_Raw/NLCD",
            force.redo = FALSE
          )
          
          # Masking
          fire_proj_nlcd <- st_transform(single_fire, crs(nlcd))
          nlcd_masked <- crop(nlcd, fire_proj_nlcd)
          nlcd_masked <- mask(nlcd_masked, fire_proj_nlcd)
          
          # Get Frequency Table (This contains Text Strings now!)
          freq_table <- terra::freq(nlcd_masked)
          colnames(freq_table) <- c("layer", "value", "count")
          
          # Calculate Total Pixels (excluding NAs)
          total_pixels <- sum(freq_table$count, na.rm = TRUE)
          
          # Helper function to match TEXT strings
          calc_pct <- function(target_names) {
            count <- freq_table %>% 
              filter(value %in% target_names) %>% 
              summarize(sum(count, na.rm=TRUE)) %>% 
              pull()
            
            if(length(count) == 0) return(0)
            if(total_pixels == 0) return(0)
            return((count / total_pixels) * 100)
          }
          
          # --- UPDATED CATEGORIES (MATCHING YOUR PRINTOUT) ---
          
          # Forest (41, 42, 43)
          pct_forest <- calc_pct(c("Deciduous Forest", "Evergreen Forest", "Mixed Forest"))
          
          # Urban (21, 22, 23, 24)
          # Note: Copied exact spelling from your diagnostic logs
          pct_urban <- calc_pct(c("Developed, Open Space", 
                                  "Developed, Low Intensity", 
                                  "Developed, Medium Intensity", 
                                  "Developed High Intensity"))
          
          # Agriculture (81, 82)
          pct_ag <- calc_pct(c("Pasture/Hay", "Cultivated Crops"))
          
          # Shrub/Grass (52, 71) - Critical for NV/WA fires
          pct_shrub <- calc_pct(c("Shrub/Scrub", "Grassland/Herbaceous"))
          
          # ----------------------------------------------------
          # 📝 PART C: SAVE TO LIST
          # ----------------------------------------------------
          results_list[[i]] <- data.frame(
            Event_ID = fire_id,
            Elev_Mean_m = elev_stats$mean,
            Elev_Roughness = elev_stats$sd,
            Slope_Mean_deg = slope_stats$mean,
            Pct_Forest_2019 = pct_forest,
            Pct_Urban_2019 = pct_urban,
            Pct_Ag_2019 = pct_ag,
            Pct_Shrub_Grass_2019 = pct_shrub
          )
          
          # Clean Memory
          rm(dem, dem_masked, slope, nlcd, nlcd_masked, fire_proj, fire_proj_nlcd)
          gc()
          
        }, error = function(e) {
          message(paste("❌ Skipped", fire_id, ":", e$message))
          results_list[[i]] <- NULL
        })
        
        # 💾 SAVE PROGRESS EVERY 50 FIRES
        if (i %% 50 == 0) {
          temp_df <- bind_rows(results_list)
          write.csv(temp_df, output_file, row.names = FALSE)
          message("--- Progress Saved to CSV ---")
        }
      }
      
      # ====================================================
      # 3. FINAL SAVE
      # ====================================================
     
      
      

################################################################################
#Step 9: Cleaning and saving metadata
      feddata_results <- bind_rows(results_list)
      Fire_boundary_Final <- left_join(fires_to_process, feddata_results, by = "Event_ID")
      
      # Save final clean CSV
      write.csv(st_drop_geometry(Fire_boundary_Final), "/Users/jackboyle/Downloads/UVA/Code/projects/Wildfire_ET_project/processed_data/#1_Original_fire_data/Metadata/Final_Fire_Landscape_Data.csv", row.names = FALSE)
      
      message("✅ PROCESSING COMPLETE! Data saved to 'Final_Fire_Landscape_Data.csv'")
      

    
      
################################################################################
#Step 10: Cleaning and saving gpkg data
      
      st_write(Fire_boundary_clean, "/Users/jackboyle/Downloads/UVA/Code/projects/Wildfire_ET_project/processed_data/#1_Original_fire_data/GPKG/Cleaned_Fire_Boundaries.gpkg", append = FALSE)
      
      message("✅ Success! 'Cleaned_Fire_Boundaries.gpkg' has been saved.")
      
      
      
      
      
      
      
      
      
      
      
      
      




