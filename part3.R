#Part 3
#giving burn seventies
#add three columns with % of high, moderate and low within the fire boundary.
#then try to find % canopy loss for every year after the wildfire.


#Step 1 = Make a Raster file for every year


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


################################################################################


#STEP #1
# ================================================================
# Extract burned-area rasters for each wildfire polygon
# ================================================================
# ---- Libraries ----
            
            library(sf)
            library(terra)
            library(stringr)
            
            # ---- Inputs from you ----
            fires_shp <- "/Users/jackboyle/Downloads/UVA/Code/projects/Wildfire_ET_project/processed_data/#3_Unburned_huc_shape/fires_hucs_with_unburned.gpkg"
            mtbs_root <- "/Users/jackboyle/Downloads/UVA/Database/WILDFIRE_WATERQUALITY/MTBS/composite_data/MTBS_BSmosaics"
            out_dir   <- "/Users/jackboyle/Downloads/UVA/Code/projects/Wildfire_ET_project/processed_data/#4_fire_burn_data_tif"
            
            # ---- Prep ----
            dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
            
            # Read fires (NAD83, MULTIPOLYGON)
            fires <- st_read(fires_shp, quiet = TRUE)
            if (!("Event_ID" %in% names(fires))) stop("Event_ID column not found in fires shapefile.")
            if (!("Ig_Date"  %in% names(fires))) stop("Ig_Date column not found in fires shapefile.")
            
            # Ensure Ig_Date is Date and derive year
            fires$Ig_Date <- as.Date(fires$Ig_Date)
            fires$year <- as.integer(format(fires$Ig_Date, "%Y"))
            if (any(is.na(fires$year))) stop("Some Ig_Date could not be parsed to a year.")
            
            # Keep only rows with a year we have a raster for (folder must exist)
            fires$mtbs_path <- file.path(mtbs_root, fires$year,
                                         paste0("mtbs_CONUS_", fires$year),
                                         paste0("mtbs_CONUS_", fires$year, ".tif"))
            
            # ---- Process by year (loads each big mosaic once) ----
            years <- sort(unique(fires$year))
            
            for (yr in years) {
              yr_rows <- fires[fires$year == yr, ]
              mtbs_tif <- unique(yr_rows$mtbs_path)
              
              if (length(mtbs_tif) != 1) {
                warning("Multiple or zero MTBS paths found for year ", yr, "; skipping this year.")
                next
              }
              if (!file.exists(mtbs_tif)) {
                warning("MTBS file missing for year ", yr, ": ", mtbs_tif)
                next
              }
              
              message("== Year ", yr, " == Loading: ", mtbs_tif)
              r_year <- rast(mtbs_tif)  # Large mosaic
              
              # Transform fires of this year to raster CRS once
              fires_yr <- st_transform(yr_rows, crs(r_year))
              v_yr <- vect(fires_yr)  # terra SpatVector
              
              # Loop fires within year
              for (i in seq_len(nrow(fires_yr))) {
                event_id <- as.character(fires_yr$Event_ID[i])
                if (is.na(event_id) || event_id == "") {
                  event_id <- paste0("fire_", yr, "_", i)
                }
                
                out_path <- file.path(out_dir, paste0(event_id, ".tif"))
                if (file.exists(out_path)) {
                  message("🟡 Exists, skipping: ", event_id)
                  next
                }
                
                message("Processing fire ", i, "/", nrow(fires_yr), " (", event_id, ")")
                
                # Geometry for this fire
                v_fire <- v_yr[i]
                
                # Crop first (fast), then mask to polygon
                r_crop <- try(crop(r_year, v_fire), silent = TRUE)
                if (inherits(r_crop, "try-error")) {
                  warning("  ❌ crop failed for ", event_id)
                  next
                }
                r_mask <- try(mask(r_crop, v_fire), silent = TRUE)
                if (inherits(r_mask, "try-error")) {
                  warning("  ❌ mask failed for ", evnet_id)
                  next
                }
                
                # Keep only classes 2, 3, 4; set everything else to NA
                # This is efficient and explicit for categorical rasters
                r_keep <- classify(
                  r_mask,
                  rcl = matrix(c(
                    -Inf, 1.9999, NA,   # <2 -> NA (covers 0 and 1)
                    2,    2,      2,
                    3,    3,      3,
                    4,    4,      4,
                    4.0001, Inf,  NA    # >4 -> NA (covers 5,6,…)
                  ), ncol = 3, byrow = TRUE),
                  include.lowest = TRUE, right = TRUE
                )
                
                # Optional quick check (silent if huge):
                # print(freq(r_keep))
                
                # Write per-fire GeoTIFF
                try({
                  writeRaster(r_keep, out_path, overwrite = TRUE)
                  message("✅ Saved: ", out_path)
                }, silent = FALSE)
              }
            }
            
            message("🔥 Finished all years.")

################################################################################
#Step #2: NLCD TCC is maped wrong... this code ajusts the pixels so its the same as MTBS
            
#Run a forloop that goes throught all . Save site data and transforms it to a beter crs
#WARNING AT THE BOTTOM (READ IT)
          
          library(terra)
          in_dir <- "/Users/jackboyle/Downloads/UVA/Database/NLC_TCC"
          files <- list.files(in_dir, pattern = "nlcd_tcc_conus_wgs84_v2023-5_\\d{8}_\\d{8}\\.tif$", 
                              recursive = TRUE, full.names = TRUE)
          out_dir <- "/Users/jackboyle/Downloads/UVA/Database/NLC_TCC/reprojected_5070"
          dir.create(out_dir, showWarnings = FALSE)
          target_crs <- "EPSG:5070"
          for (f in files) {
            r <- rast(f)
            r_5070 <- project(r, target_crs, method = "near")
            fname <- basename(f)
            out_file <- file.path(out_dir, sub(".tif", "_5070.tif", fname))
            writeRaster(r_5070, out_file, overwrite = TRUE)
            cat("✅ Reprojected:", fname, "\n")
          }
          cat("\nAll files reprojected to EPSG:5070 and saved to:\n", out_dir, "\n")
    
 #WARNING: This code will take many days to run, ask Jack for the reporjected files  
#
          #Lets see how "badly" the pixels are missalinged
################################################################################
#STEP #3: The forest cover (m^2) for every pre-fire environment

          # === Libraries ===
          library(sf)
          library(terra)
          library(dplyr)
          
          # === File paths ===
          fires_path <- "/Users/jackboyle/Downloads/UVA/Code/projects/Wildfire_ET_project/processed_data/#3_Unburned_huc_shape/fires_hucs_with_unburned.gpkg"
          burned_dir <- "/Users/jackboyle/Downloads/UVA/Code/projects/Wildfire_ET_project/processed_data/#4_fire_burn_data_tif"
          tcc_base_dir <- "/Users/jackboyle/Downloads/UVA/Database/NLC_TCC/reprojected_5070"
          out_gpkg <- "/Users/jackboyle/Downloads/UVA/Code/projects/Wildfire_ET_project/processed_data/#5_Pre_fire_canopy/prefire_canopy_snap.gpkg"
          
          # === Load shapefile ===
          fires <- st_read(fires_path, quiet = TRUE)
          
          # Ensure proper date format
          fires$Ig_Date <- as.Date(fires$Ig_Date)
          
          # Add output columns
          fires$pre_fire_canopy_m2 <- NA_real_
          fires$pre_fire_total_area_m2 <- NA_real_
          fires$pre_fire_canopy_pct <- NA_real_
          
          # === Helper: Find TCC raster for year before fire ===
          find_tcc_for_year <- function(prev_year, base_dir) {
            pattern <- paste0(prev_year, "0101_", prev_year, "1231_5070\\.tif$")
            files <- list.files(base_dir, pattern = pattern, full.names = TRUE)
            if (length(files) == 0) return(NA_character_)
            return(files[1])
          }
          
          # === Main loop ===
          for (i in seq_len(nrow(fires))) {
            ev_id <- fires$Event_ID[i]
            ig_date <- fires$Ig_Date[i]
            
            if (is.na(ig_date)) {
              message(sprintf("[%d/%d] %s: missing ignition date -> skip", i, nrow(fires), ev_id))
              next
            }
            
            fire_year <- as.integer(format(ig_date, "%Y"))
            prev_year <- fire_year - 1
            
            message(sprintf("[%d/%d] Processing %s (fire year %d, using TCC %d)", 
                            i, nrow(fires), ev_id, fire_year, prev_year))
            
            burn_file <- file.path(burned_dir, paste0(ev_id, ".tif"))
            if (!file.exists(burn_file)) {
              warning("  - Burn raster not found -> skip")
              next
            }
            
            tcc_file <- find_tcc_for_year(prev_year, tcc_base_dir)
            if (is.na(tcc_file) || !file.exists(tcc_file)) {
              warning("  - No TCC file for ", prev_year, " -> skip")
              next
            }
            
            # === Load rasters ===
            burn <- try(rast(burn_file), silent = TRUE)
            tcc  <- try(rast(tcc_file),  silent = TRUE)
            if (inherits(burn, "try-error") || inherits(tcc, "try-error")) {
              warning("  - Failed to load raster(s) -> skip")
              next
            }
            
            # === Align and snap ===
            if (!same.crs(burn, tcc)) {
              message("  - Reprojecting TCC to match burn CRS...")
              tcc <- try(project(tcc, burn, method = "near"), silent = TRUE)
              if (inherits(tcc, "try-error")) {
                warning("  - Reprojection failed -> skip")
                next
              }
            }
            
            # Safe snapping checks
            burn_origin <- try(origin(burn), silent = TRUE)
            tcc_origin  <- try(origin(tcc),  silent = TRUE)
            burn_res    <- try(res(burn),    silent = TRUE)
            tcc_res     <- try(res(tcc),     silent = TRUE)
            
            if (
              !is.numeric(burn_origin) || !is.numeric(tcc_origin) ||
              !is.numeric(burn_res) || !is.numeric(tcc_res)
            ) {
              warning("  - Invalid raster geometry -> skip")
              next
            }
            
            # Snap to burn if needed
            if (
              !isTRUE(all.equal(tcc_origin, burn_origin, tolerance = 1e-6)) ||
              !isTRUE(all.equal(tcc_res, burn_res, tolerance = 1e-6))
            ) {
              message("  - Snapping TCC to burn alignment...")
              tcc <- try(resample(tcc, burn, method = "near"), silent = TRUE)
              if (inherits(tcc, "try-error")) {
                warning("  - Snap failed -> skip fire")
                next
              }
            }
            
            # === Compute canopy area ===
            burn_mask <- burn > 0
            tcc_crop <- try(crop(tcc, burn_mask), silent = TRUE)
            if (inherits(tcc_crop, "try-error")) {
              warning("  - Crop failed -> skip")
              next
            }
            
            tcc_masked <- try(mask(tcc_crop, burn_mask), silent = TRUE)
            if (inherits(tcc_masked, "try-error")) {
              warning("  - Mask failed -> skip")
              next
            }
            
            canopy_area_raster <- tcc_masked / 100 * (30 * 30)
            total_canopy_m2 <- try(global(canopy_area_raster, "sum", na.rm = TRUE)[[1]], silent = TRUE)
            burned_area_m2  <- try(global(burn_mask, "sum", na.rm = TRUE)[[1]] * (30 * 30), silent = TRUE)
            
            if (inherits(total_canopy_m2, "try-error") || inherits(burned_area_m2, "try-error")) {
              warning("  - Failed to compute canopy -> skip")
              next
            }
            
            canopy_pct <- (total_canopy_m2 / burned_area_m2) * 100
            
            fires$pre_fire_canopy_m2[i] <- total_canopy_m2
            fires$pre_fire_total_area_m2[i] <- burned_area_m2
            fires$pre_fire_canopy_pct[i] <- canopy_pct
            
            message(sprintf("  - ✅ Canopy area: %.1f m² (%.2f%% of fire area)", 
                            total_canopy_m2, canopy_pct))
          }
          
          # === Save results ===
          st_write(fires, out_gpkg, delete_dsn = TRUE)
          message("🔥 Finished all fires — results saved to: ", out_gpkg)


#It should be noted that there are no burn rasters for ALASKA... this is something to look into in the future.

################################################################################
#STEP 4: Change the LOSS tif to 30x30 instead of 21x21
          #this will be better for our analysis
          
          library(terra)
          
          # ==============================================================================
          # 1. SETUP PATHS
          # ==============================================================================
          hansen_raw_folder  <- "/Users/jackboyle/Downloads/UVA/Database/GFC/Loss"
          mtbs_template_path <- "/Users/jackboyle/Downloads/UVA/Database/WILDFIRE_WATERQUALITY/MTBS/composite_data/MTBS_BSmosaics/2024/mtbs_CONUS_2024/mtbs_CONUS_2024.tif"
          out_file           <- "/Users/jackboyle/Downloads/UVA/Database/GFC/Loss/Hansen_US_5070_aligned.tif"
          
          # ==============================================================================
          # 2. CREATE THE VIRTUAL MOSAIC (The Dough)
          # ==============================================================================
          message("Building Virtual Mosaic from raw tiles...")
          raw_files <- list.files(hansen_raw_folder, pattern = "^Hansen_GFC.*\\.tif$", full.names = TRUE)
          if (length(raw_files) == 0) stop("No Hansen files found!")
          hansen_vrt <- vrt(raw_files)
          
          # ==============================================================================
          # 3. FIX THE TEMPLATE (The Cookie Cutter)
          # ==============================================================================
          message("Loading MTBS Template...")
          r_mtbs <- rast(mtbs_template_path)
          
          # Check the current top edge
          current_top <- ymax(r_mtbs)
          message(sprintf("Original Template Top Edge: %.0f meters", current_top))
          
          # We need to reach approx 3,170,000 meters to hit the 49th Parallel (Canada).
          # Let's extend it to 3,250,000 just to be safe.
          target_top <- 3250000
          
          if (current_top < target_top) {
            message(sprintf("⚠️ Template is too short! Extending canvas North to %.0f meters...", target_top))
            
            # Create a custom extent: Keep Left, Right, Bottom the same. Stretch Top.
            # (xmin, xmax, ymin, ymax)
            new_extent <- ext(xmin(r_mtbs), xmax(r_mtbs), ymin(r_mtbs), target_top)
            
            # 'extend' increases the raster size while keeping the grid perfectly aligned
            r_template_big <- extend(r_mtbs, new_extent)
          } else {
            message("Template is already big enough.")
            r_template_big <- r_mtbs
          }
          
          # ==============================================================================
          # 4. PROJECT AND SNAP
          # ==============================================================================
          message("Reprojecting Hansen data onto the EXPANDED template...")
          message("(This relies on Nearest Neighbor to preserve years)")
          
          project(hansen_vrt, r_template_big, 
                  method = "near", 
                  filename = out_file, 
                  overwrite = TRUE,
                  wopt = list(datatype = "INT1U"))
          
          message("✅ Done! Map saved.")
          
          # ==============================================================================
          # 5. VERIFY
          # ==============================================================================
          r_new <- rast(out_file)
          cat("\nNew Map Extent:\n")
          print(ext(r_new))
          
          if (ymax(r_new) > 3100000) {
            message("🎉 SUCCESS! The map now reaches the Canadian border.")
          } else {
            message("⚠️ Still short? Something is blocking the extent.")
          }
          
          
          #Once we have this we can compare Loss with Burn and TCC pixel by pixel (they will be fully the same)
          
          
          
################################################################################
#STEP# 4.A  A Quick Check (CRS, Resolution, Origin, Extent)
          
          library(terra)
          hansen_file <- "/Users/jackboyle/Downloads/UVA/Database/GFC/Loss/Hansen_US_5070_aligned.tif"
          mtbs_fire_file <- "/Users/jackboyle/Downloads/UVA/Code/projects/Wildfire_ET_project/processed_data/#4_fire_burn_data_tif/MT4721611299020070715.tif"
        
          hansen <- rast(hansen_file)
          mtbs_fire <- rast(mtbs_fire_file)
          
         
          cat("\nCRS identical? ", crs(hansen) == crs(mtbs_fire), "\n\n")
          cat("Resolution (Hansen): ", res(hansen), "\n")
          cat("Resolution (MTBS fire): ", res(mtbs_fire), "\n\n")
          cat("Origin (Hansen): ", origin(hansen), "\n")
          cat("Origin (MTBS fire): ", origin(mtbs_fire), "\n\n")
          print(ext(hansen))
          print(ext(mtbs_fire))
          compareGeom(hansen, mtbs_fire, stopOnError = FALSE)
          
################################################################################
#STEP #5 Using vector perimeter data, find the loss of canopy in every wildfire event for the fire_year, post_year_1, post_year_2, post_year_2 
          
          library(terra)
          library(dplyr)
          
          # --- Define paths ---
          fires_path  <- "/Users/jackboyle/Downloads/UVA/Code/projects/Wildfire_ET_project/processed_data/#5_Pre_fire_canopy/prefire_canopy_snap.gpkg"
          tif_folder  <- "/Users/jackboyle/Downloads/UVA/Code/projects/Wildfire_ET_project/processed_data/#4_fire_burn_data_tif"
          hansen_path <- "/Users/jackboyle/Downloads/UVA/Database/GFC/Loss/Hansen_US_5070_aligned.tif"
          
          output_path <- "/Users/jackboyle/Downloads/UVA/Code/projects/Wildfire_ET_project/processed_data/#6_post_pre_fire_canopy/post_pre_canopy_fire.gpkg"
          
          # --- Load data ---
          message("Loading Vectors and Rasters...")
          fires <- vect(fires_path)
          hansen_raster <- rast(hansen_path)
          tif_files <- list.files(tif_folder, pattern = "\\.tif$", full.names = TRUE)
          
          # --- Check Column Names (Safety) ---
          # We standardize to 'Event_ID' for the script to work
          if ("Evnt_ID" %in% names(fires) && !("Event_ID" %in% names(fires))) {
            names(fires)[names(fires) == "Evnt_ID"] <- "Event_ID"
          }
          if (!"Event_ID" %in% names(fires)) stop("Could not find Event_ID column!")
          
          # --- Initialize result columns (Expanded to 5 Years) ---
          fires$fire_year_loss_m2   <- NA
          fires$post_fire_year_1_m2 <- NA
          fires$post_fire_year_2_m2 <- NA
          fires$post_fire_year_3_m2 <- NA
          fires$post_fire_year_4_m2 <- NA # NEW
          fires$post_fire_year_5_m2 <- NA # NEW
          
          # --- Loop through each fire ---
          total_fires <- nrow(fires)
          
          for (i in 1:total_fires) {
            
            first_fire_id <- fires$Event_ID[i]
            fire_date     <- fires$Ig_Date[i]
            
            if (is.na(fire_date) | is.na(first_fire_id)) {
              # cat("Skipping fire", i, "- missing ID or date\n")
              next
            }
            
            fire_year <- as.numeric(format(as.Date(fire_date), "%Y"))
            
            # Optional: Print progress every 50 fires to reduce console spam
            if (i %% 50 == 0) cat("Processing fire", i, "/", total_fires, "| Year:", fire_year, "\n")
            
            # --- Find matching burn raster ---
            # fixed=TRUE ensures we don't accidentally match 'Fire_1' with 'Fire_10'
            matching_tif <- tif_files[grep(first_fire_id, tif_files, fixed=TRUE)]
            
            if (length(matching_tif) == 0) {
              # cat("  ⚠️ No matching burn raster found for:", first_fire_id, "\n")
              next
            }
            
            # --- Try loading and processing ---
            tryCatch({
              burn_raster <- rast(matching_tif[1])
              
              # Crop + mask Hansen raster (Efficient)
              hansen_crop <- crop(hansen_raster, burn_raster)
              hansen_mask <- mask(hansen_crop, burn_raster)
              
              # Calculate Pixel Area (usually 900m2 for 30m pixels)
              # We calculate it dynamically just in case resolution varies slightly
              pixel_area <- prod(res(hansen_mask))
              
              # Get pixel frequencies (Counts of each Year)
              loss_vals <- freq(hansen_mask)
              
              # Check if freq returned empty or valid data
              if (nrow(loss_vals) == 0) next
              
              loss_vals <- as.data.frame(loss_vals) %>% filter(!is.na(value))
              
              if (nrow(loss_vals) == 0) next
              
              # Convert Hansen Values (0-24) to Actual Years (2000-2024)
              # Note: Value 1 = 2001, Value 15 = 2015, etc.
              loss_summary <- loss_vals %>%
                mutate(
                  loss_year = 2000 + value,
                  lost_area_m2 = count * pixel_area
                )
              
              # --- Extract Data for Fire Year + 5 Years ---
              # We use sum(..., na.rm=TRUE) so if no pixels match, it returns 0 (Correct).
              
              fires$fire_year_loss_m2[i]   <- sum(loss_summary$lost_area_m2[loss_summary$loss_year == fire_year], na.rm = TRUE)
              
              fires$post_fire_year_1_m2[i] <- sum(loss_summary$lost_area_m2[loss_summary$loss_year == (fire_year + 1)], na.rm = TRUE)
              fires$post_fire_year_2_m2[i] <- sum(loss_summary$lost_area_m2[loss_summary$loss_year == (fire_year + 2)], na.rm = TRUE)
              fires$post_fire_year_3_m2[i] <- sum(loss_summary$lost_area_m2[loss_summary$loss_year == (fire_year + 3)], na.rm = TRUE)
              fires$post_fire_year_4_m2[i] <- sum(loss_summary$lost_area_m2[loss_summary$loss_year == (fire_year + 4)], na.rm = TRUE)
              fires$post_fire_year_5_m2[i] <- sum(loss_summary$lost_area_m2[loss_summary$loss_year == (fire_year + 5)], na.rm = TRUE)
              
            }, error = function(e) {
              cat("  ❌ Error processing", first_fire_id, ":", e$message, "\n")
            })
          }
          
          # --- Review and Save ---
          print(head(fires[, c("Event_ID", "fire_year_loss_m2", "post_fire_year_1_m2", "post_fire_year_5_m2")]))
          
          writeVector(fires, output_path, overwrite = TRUE)
          cat("\n✅ Saved updated fires with 5-year forest loss data to:\n", output_path, "\n")

         
           
          #IMPORTANT NOTE:
          #I had to run the cunk above on the RIVANA SUPERCOMPUTER because it R kept Aborting.
#################################################################################################################################
#Step 6: Filter and add the QA check:
          
          fires<-st_read("/Users/jackboyle/Downloads/UVA/Code/projects/Wildfire_ET_project/processed_data/#6_post_pre_fire_canopy/raster_method_results/post_pre_canopy_fire.gpkg")
          
          message("--- STARTING PART C: QA/QC CHECKS ---")
          
          # 1. Calculate Sum of Recorded Loss
          loss_cols <- c("fire_year_loss_m2", 
                         "post_fire_year_1_m2", "post_fire_year_2_m2", 
                         "post_fire_year_3_m2", "post_fire_year_4_m2", 
                         "post_fire_year_5_m2")
          
          # We use rowSums on the data (dropping geometry for speed)
          fires$total_recorded_loss_m2 <- rowSums(st_drop_geometry(fires)[, loss_cols], na.rm = TRUE)
          
          # 2. Logic Check
          # Impossible if: Total Loss > Total Polygon Area  OR  Total Loss > Starting Canopy
          fires$qa_status <- ifelse(
            (fires$total_recorded_loss_m2 > fires$pre_fire_total_area_m2) | 
              (fires$total_recorded_loss_m2 > fires$pre_fire_canopy_m2), 
            "Impossible", 
            "Possible"
          )
          
          # 3. Ratio Metric (Helpful for filtering later)
          fires$loss_to_canopy_ratio <- fires$total_recorded_loss_m2 / fires$pre_fire_canopy_m2
          
         
          
          final_path <- "/Users/jackboyle/Downloads/UVA/Code/projects/Wildfire_ET_project/processed_data/#6_post_pre_fire_canopy/raster_method_results/raster_method_pre_post_canopy_fire.gpkg"
          fires_clean <- fires %>%
            filter(fire_year_loss_m2 > 0) 
          st_write(fires_clean, final_path, delete_dsn=TRUE)
          
          #I would like to make a note here: many of the fires will show a QA of "Impossible"
          

###################################################################################
#############################################################################
###################################################################################
#############################################################################
#############################################################################  
############################################################################# 
#############################################################################
          #METHOD #2 ~ VECTOR ANALYSIS ~ Anti-pixel  comparison
          #In this method all we will be depending upon the entire MTBS perimeter shape
          #We will take all of the NLCD TCC trees inside the perimeter
          #We will take all of the GFC LOSS pixels inside the perimeter
          # This should be easier than the "pixel mask" method
          #Also, this method prevents the issuer of converting LOSS pixels from 21x21 to 30x30.
#############################################################################
# Step #1: Prep & Libraries
          
          library(sf)
          library(terra)
          library(dplyr)
          library(exactextractr) 
          
          # --- PATHS ---
          # Input Fire Polygons
          fires_path  <- "/Users/jackboyle/Downloads/UVA/Code/projects/Wildfire_ET_project/processed_data/#5_Pre_fire_canopy/prefire_canopy_snap.gpkg"
          
          # Big Raster Data
          hansen_path <- "/Users/jackboyle/Downloads/UVA/Database/GFC/Loss/Hansen_US_5070_aligned.tif"
          tcc_dir     <- "/Users/jackboyle/Downloads/UVA/Database/NLC_TCC/reprojected_5070" 
          
          # Outputs
          out_raw     <- "/Users/jackboyle/Downloads/UVA/Code/projects/Wildfire_ET_project/processed_data/#6_post_pre_fire_canopy/Vector_method_results/vector_method_RAW.gpkg"
          out_clean   <- "/Users/jackboyle/Downloads/UVA/Code/projects/Wildfire_ET_project/processed_data/#6_post_pre_fire_canopy/Vector_method_results/vector_method_FILTERED.gpkg"
          
          if (!dir.exists(dirname(out_raw))) dir.create(dirname(out_raw), recursive = TRUE)
          
#############################################################################  
# Step #2: Prepare the Dataframe      
          
          message("Loading Fire Polygons...")
          fires <- st_read(fires_path, quiet = TRUE)
          
          # Standardize Event_ID
          if ("Evnt_ID" %in% names(fires) && !("Event_ID" %in% names(fires))) {
            names(fires)[names(fires) == "Evnt_ID"] <- "Event_ID"
          }
          
          fires$Ig_Date <- as.Date(fires$Ig_Date)
          
          # Initialize columns
          fires$pre_fire_total_area_m2 <- NA 
          fires$pre_fire_canopy_m2     <- NA
          fires$pre_fire_canopy_pct    <- NA
          fires$fire_year_loss_m2      <- NA
          fires$post_fire_year_1_m2    <- NA
          fires$post_fire_year_2_m2    <- NA
          fires$post_fire_year_3_m2    <- NA
          fires$post_fire_year_4_m2    <- NA
          fires$post_fire_year_5_m2    <- NA
          
#############################################################################            
# Step #3: Pre-Fire Canopy (Using Vector Perimeter)
          
          message("--- STARTING PART A: PRE-FIRE CANOPY (VECTOR METHOD) ---")
          
          unique_years <- sort(unique(as.numeric(format(fires$Ig_Date, "%Y"))))
          
          for (yr in unique_years) {
            
            prev_year <- yr - 1
            message(paste("Processing fires for year:", yr, "| Looking for TCC:", prev_year))
            
            tcc_pat <- paste0(prev_year, "0101_", prev_year, "1231_5070\\.tif$")
            tcc_file <- list.files(tcc_dir, pattern = tcc_pat, full.names = TRUE)
            
            if (length(tcc_file) == 0) {
              message(paste("  ⚠️ TCC file not found for year", prev_year, "- skipping these fires."))
              next
            }
            
            r_tcc <- rast(tcc_file[1])
            
            idx <- which(as.numeric(format(fires$Ig_Date, "%Y")) == yr)
            fires_subset <- fires[idx, ]
            
            if (st_crs(fires_subset) != st_crs(r_tcc)) {
              fires_subset <- st_transform(fires_subset, st_crs(r_tcc))
            }
            
            # Calculate Mean Canopy % inside Polygon
            tcc_means <- exact_extract(r_tcc, fires_subset, "mean", progress = FALSE)
            
            # Calculate Total Area based on Polygon Geometry
            poly_areas <- as.numeric(st_area(fires_subset))
            
            fires$pre_fire_total_area_m2[idx] <- poly_areas
            fires$pre_fire_canopy_pct[idx]    <- tcc_means
            fires$pre_fire_canopy_m2[idx]     <- (tcc_means / 100) * poly_areas
            
            rm(r_tcc, fires_subset, tcc_means)
            gc()
          }
#############################################################################  
# Step #4: Post-fire Forest Loss (Using Hansen GFC + Vector)
          
          message("--- STARTING PART B: FOREST LOSS (VECTOR METHOD) ---")
          
          hansen <- rast(hansen_path)
          total <- nrow(fires)
          px_area <- 900 
          
          for (i in 1:total) {
            
            if (i %% 50 == 0) message(paste("Processing Loss for fire", i, "/", total))
            
            fire_year <- as.numeric(format(fires$Ig_Date[i], "%Y"))
            if (is.na(fire_year)) next
            
            tryCatch({
              extract_data <- exact_extract(hansen, fires[i, ], progress = FALSE)[[1]]
              
              if (nrow(extract_data) > 0) {
                loss_summary <- extract_data %>%
                  filter(value != 0) %>%
                  mutate(year = 2000 + value) %>%
                  group_by(year) %>%
                  summarise(weighted_pixels = sum(coverage_fraction)) %>%
                  mutate(area_m2 = weighted_pixels * px_area)
                
                get_loss <- function(y) {
                  val <- loss_summary$area_m2[loss_summary$year == y]
                  if (length(val) == 0) return(0) else return(val)
                }
                
                fires$fire_year_loss_m2[i]   <- get_loss(fire_year)
                fires$post_fire_year_1_m2[i] <- get_loss(fire_year + 1)
                fires$post_fire_year_2_m2[i] <- get_loss(fire_year + 2)
                fires$post_fire_year_3_m2[i] <- get_loss(fire_year + 3)
                fires$post_fire_year_4_m2[i] <- get_loss(fire_year + 4)
                fires$post_fire_year_5_m2[i] <- get_loss(fire_year + 5)
              } else {
                fires$fire_year_loss_m2[i] <- 0
              }
              
            }, error = function(e) {
              message(paste("Error on fire", i, ":", e$message))
            })
            
            if (i %% 100 == 0) gc()
          }
          
#############################################################################  
# Step #5: QA/QC 
          
          message("--- STARTING PART C: QA/QC CHECKS ---")
          
          loss_cols <- c("fire_year_loss_m2", "post_fire_year_1_m2", "post_fire_year_2_m2", 
                         "post_fire_year_3_m2", "post_fire_year_4_m2", "post_fire_year_5_m2")
          
          fires$total_recorded_loss_m2 <- rowSums(st_drop_geometry(fires)[, loss_cols], na.rm = TRUE)
          
          fires$qa_status <- ifelse(
            (fires$total_recorded_loss_m2 > fires$pre_fire_total_area_m2) | 
              (fires$total_recorded_loss_m2 > fires$pre_fire_canopy_m2), 
            "Impossible", "Possible"
          )
          
          fires$loss_to_canopy_ratio <- fires$total_recorded_loss_m2 / fires$pre_fire_canopy_m2
          
          message("QA Summary:")
          print(table(fires$qa_status))
          
#############################################################################  
# Step #6: Save RAW Data (Includes Zeros and Impossible fires)
          
          message(paste("Saving RAW results to:", out_raw))
          st_write(fires, out_raw, delete_dsn = TRUE)
          
#############################################################################  
# Step #7: Filter and Save CLEAN Data
          
          message("Filtering for valid fires...")
          
          fires_filtered <- fires %>%
            filter(fire_year_loss_m2 > 0)
          
          message(paste("Saving FILTERED results to:", out_clean))
          st_write(fires_filtered, out_clean, delete_dsn = TRUE)
          
          message("All Vector Analysis Complete.")
          
#############################################################################
#step 8: you don't really need this step
          #Compare the two methods
          
         #Just know that the raster method is actually better
