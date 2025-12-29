#PART 2
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

#Step 1 get shape files of all HUC 4's in the CONUS
      # make a list of wildfires from the 1950 onward.
      # For loop which Intersects wildfire in every huc4... 
        #and gets rid of them from Fire-Boundary so that the code can run faster
      #Then make a huc 4 shape file for each huc 4 with each fire from 1950-2024 in it.

#Step 2: run through the Fire_boundary ds in the projects folder...
        #... find huc 12s that are at lest 85% in ithe wildfire boundary
        # add the huc12 ids to a column in specific Fire_boundary row,
        # if there are multipble seperate them by a comma.

#Where will i store all of this data? I should run one analysis or one fire at a time...
#... or maybe I should first capture the fires with high delayed canopy mortaitly
#... or maybe I should look at fires with moderate severity burns first with MTBS... check hungs paper to see what classifies as  modert severity



################################################################################
#Step 1: Set the stage: Make a map of he US, all huc 4s and load all wildfires since 1950

      usa <- ne_countries(country = "United States of America", scale = "medium", returnclass = "sf")
      plot(usa["geometry"])
      hu4<-get_huc(AOI = usa, id = NULL, t_srs = NULL, buffer = 0.5, type = "huc04")
      hu4 <- hu4 %>%
        filter(!states %in% c("HI", "AK"))
      hu4 <- st_make_valid(hu4)
      
      Fire_boundary <- read_sf("/Users/jackboyle/Downloads/UVA/Database/WILDFIRE_WATERQUALITY/MTBS/mtbs_perimeter_data/mtbs_perims_DD.shp")
      
      Fire_boundary <- Fire_boundary %>%
        mutate(Ig_Date = as.Date(Ig_Date)) %>%
        filter(format(Ig_Date, "%Y") >= 1950)


#############################################################
#Step 2: Seperating all wildfires from 1950 into their huc4 zones

      Fire_boundary <- st_transform(Fire_boundary, st_crs(hu4)) %>%
        st_make_valid()
      if (!"fire_id" %in% names(Fire_boundary)) {
        Fire_boundary <- Fire_boundary %>%
          mutate(fire_id = row_number())
      }
      fires_remaining <- Fire_boundary
      outdir <- "/Users/jackboyle/Downloads/UVA/Code/projects/Wildfire_ET_project/processed_data/#2_Burned_huc_shape/Fires_in_HUC4"
      huc4_id_col <- "huc4"  # change if your column has a different name
      for (i in 1:nrow(hu4)) {
        this_huc <- hu4[i, ]
        huc_id <- this_huc[[huc4_id_col]]
        fires_in_huc <- st_filter(fires_remaining, this_huc, .predicate = st_intersects)
        if (nrow(fires_in_huc) > 0) {
          outfile <- file.path(outdir, paste0("HUC4_", huc_id, "_fires.shp"))
          st_write(fires_in_huc, outfile, delete_layer = TRUE)
          
          fires_remaining <- fires_remaining %>%
            filter(!fire_id %in% fires_in_huc$fire_id)
        }
        cat("Finished HUC4:", huc_id, " | Remaining fires:", nrow(fires_remaining), "\n")
      }

################################################################################
# Step 3: Finding Huc12s that are at least 50% within a wildfire area
      
        #Getting burned huc12s for each fire
      
      library(sf)
      library(dplyr)
      
      # --- input fire shapefile ---
      filtered_fires <- st_read(paste0(
        "/Users/jackboyle/Downloads/UVA/Code/projects/Wildfire_ET_project/processed_data/#1_Original_fire_data/GPKG/Cleaned_Fire_Boundaries.gpkg"
      ))
      
      # --- parameters ---
      huc_keep_thresh <- 0.50
      out_shp <- paste0(
        "/Users/jackboyle/Downloads/UVA/Code/projects/Wildfire_ET_project/processed_data/",
        "#2_Burned_huc_shape/filtered_fires_with_huc12.gpkg"
      )
      area_crs <- 5070   # Conus Albers equal-area
      
      # --- prepare data ---
      filtered_fires$huc12 <- NA_character_
      
      # --- main loop ---
      for (i in seq_len(nrow(filtered_fires))) {
        fire <- filtered_fires[i, , drop = FALSE]
        cat(sprintf("Processing fire %d / %d | Evnt_ID: %s\n", 
                    i, nrow(filtered_fires), fire$Evnt_ID))
        
        # --- get intersecting HUC12 polygons ---
        huc12 <- tryCatch(
          get_huc(AOI = st_transform(fire, 4326), type = "huc12"),
          error = function(e) {
            message("  get_huc() failed: ", conditionMessage(e))
            return(NULL)
          }
        )
        
        if (is.null(huc12) || nrow(huc12) == 0) {
          filtered_fires$huc12[i] <- NA_character_
          next
        }
        
        # --- intersection (disable s2 for planar operation) ---
        sf_use_s2(FALSE)
        huc12 <- st_transform(huc12, st_crs(filtered_fires))
        inter <- tryCatch(st_intersection(huc12, fire), error = function(e) NULL)
        sf_use_s2(TRUE)
        
        if (is.null(inter) || nrow(inter) == 0) {
          filtered_fires$huc12[i] <- NA_character_
          next
        }
        
        # --- area calculations ---
        inter_proj  <- st_transform(inter, area_crs)
        huc12_proj  <- st_transform(huc12, area_crs)
        
        inter_proj$overlap_area <- as.numeric(st_area(inter_proj))
        huc_area_tbl <- st_drop_geometry(huc12_proj) %>%
          mutate(huc_area = as.numeric(st_area(huc12_proj))) %>%
          dplyr::select(huc12, huc_area)
        
        inter_proj <- inter_proj %>%
          left_join(huc_area_tbl, by = "huc12") %>%
          mutate(prop_overlap = overlap_area / huc_area)
        
        # --- keep HUCs with ≥85% overlap ---
        burned <- inter_proj %>% filter(prop_overlap >= huc_keep_thresh)
        
        if (nrow(burned) == 0) {
          filtered_fires$huc12[i] <- NA_character_
          next
        }
        
        # --- store comma-separated HUC12s ---
        burned_ids <- paste(unique(burned$huc12), collapse = ",")
        filtered_fires$huc12[i] <- burned_ids
        
        cat(sprintf("  -> %d HUCs kept: %s\n", length(unique(burned$huc12)), burned_ids))
      }
      
      
      
      # --- export ---
      st_write(filtered_fires, out_shp, delete_dsn = TRUE)
      cat(" Done! Shapefile written to:", out_shp, "\n")
################################################################################
# Step 4: Getting unburned hucs12s for each fire... by  
          #Create a 1000km buffer around the wildfire
          # Identity all wildfires in the past 50 years within that buffer
          #Select HUC12s in that buffer that have not been burned in the last 50 years
          # Make sure that the Climate Zone (koppen) in the reference huc matches to the climate zone of the wildfire zone
          # Save all the data


      library(sf)
      library(stringr)
      library(dplyr)
      library(terra)
      library(nhdplusTools)
      library(kgc)
      library(rnaturalearth) # Needed for the US boundary to get HUC4s
      
      # ==============================================================================
      # 1. SETUP PATHS & PARAMETERS
      # ==============================================================================
      
      input_path <- paste0(
        "/Users/jackboyle/Downloads/UVA/Code/projects/Wildfire_ET_project/processed_data/",
        "#2_Burned_huc_shape/filtered_fires_with_huc12.gpkg"
      )
      
      out_shp <- "/Users/jackboyle/Downloads/UVA/Code/projects/Wildfire_ET_project/processed_data/#3_Unburned_huc_shape/fires_hucs_with_unburned.gpkg"
      
      huc4_shp_folder <- "/Users/jackboyle/Downloads/UVA/Code/projects/Wildfire_ET_project/processed_data/#2_Burned_huc_shape/Fires_in_HUC4"
      
      koppen_tif_path <- "/Users/jackboyle/Downloads/UVA/Database/Climate_zones/koppen_geiger_tif/1991_2020/koppen_geiger_0p1.tif"
      
      buffer_m <- 10000 
      overwrite_progress <- TRUE 
      
      # ==============================================================================
      # 2. PREPARE REFERENCE DATA
      # ==============================================================================
      
      # --- A. Load Fire Data ---
      if (file.exists(out_shp) && overwrite_progress == FALSE) {
        message(" Found existing progress file. Resuming...")
        fires_out <- st_read(out_shp, quiet = TRUE)
      } else {
        message(" Starting from scratch...")
        fires_hucs_clean <- st_read(input_path, quiet = TRUE)
        fires_out <- fires_hucs_clean
        fires_out$unburned_huc <- NA_character_
      }
      
      # --- B. Load Climate Raster ---
      message("Loading Koppen Raster...")
      if (!file.exists(koppen_tif_path)) stop("Raster file not found!")
      koppen_rast <- terra::rast(koppen_tif_path)
      
      # --- C. Prepare HUC4 Reference Map (For the "Fallback" Logic) ---
      message("Downloading/Preparing HUC4 Map of USA (for fires with missing IDs)...")
      # We check if we already have it in memory to save time, otherwise download
      if (!exists("hu4_ref")) {
        usa <- ne_countries(country = "United States of America", scale = "medium", returnclass = "sf")
        # Use tryCatch in case of download blips
        hu4_ref <- tryCatch(
          get_huc(AOI = usa, type = "huc04"),
          error = function(e) {
            stop("Failed to download HUC4 reference map. Check internet.")
          }
        )
        hu4_ref <- hu4_ref %>% 
          filter(!states %in% c("HI", "AK")) %>%
          st_make_valid() %>%
          select(huc4) # We only need the code
      }
      
      # ==============================================================================
      # 3. HELPER FUNCTION
      # ==============================================================================
      load_huc4_fires <- function(huc4_codes, folder) {
        files_to_load <- unlist(lapply(huc4_codes, function(h4) {
          pat <- paste0("^HUC4_", h4, "_fires\\.shp$")
          list.files(folder, pattern = pat, full.names = TRUE)
        }))
        files_to_load <- unique(files_to_load)
        if (length(files_to_load) == 0) return(NULL)
        
        shp_list <- lapply(files_to_load, function(f) {
          tryCatch(st_read(f, quiet = TRUE), error = function(e) NULL)
        })
        
        shp_list <- Filter(Negate(is.null), shp_list)
        if (length(shp_list) == 0) return(NULL)
        
        out <- do.call(rbind, shp_list)
        return(st_make_valid(out))
      }
      
      # ==============================================================================
      # 4. MAIN LOOP
      # ==============================================================================
      
      for (i in seq_len(nrow(fires_out))) {
        
        # RESUME CHECK
        if (overwrite_progress == FALSE && !is.na(fires_out$unburned_huc[i])) next
        
        cat(sprintf("Processing fire %d / %d | Evnt_ID: %s\n", 
                    i, nrow(fires_out), fires_out$Evnt_ID[i]))
        
        fire_row <- fires_out[i, ]
        
        # --- 1. CHECK KOPPEN ---
        fire_k_str <- as.character(fire_row[["koppen_id"]])
        if (length(fire_k_str) == 0 || is.na(fire_k_str) || fire_k_str == "NA") {
          cat(" - Fire has no Koppen data. Marking NONE.\n")
          fires_out$unburned_huc[i] <- "NONE"
          next
        }
        fire_k_vec <- str_trim(unlist(str_split(fire_k_str, ",")))
        cat(sprintf("   [Fire Climate: %s]\n", paste(fire_k_vec, collapse=", ")))
        
        fires_out$unburned_huc[i] <- "NONE" # Default
        
        # --- 2. DETERMINE HUC4 (THE NEW LOGIC) ---
        burned_cell <- fire_row[["huc12"]]
        huc4_target <- NULL
        
        # Method A: Fast string extraction
        if (length(burned_cell) > 0 && !is.na(burned_cell) && as.character(burned_cell) != "") {
          # Extract first 4 digits
          huc4_target <- substr(str_trim(as.character(burned_cell)), 1, 4)
        }
        
        # Method B: Spatial Fallback (if Method A failed)
        if (is.null(huc4_target)) {
          cat("   [!] No Burned HUC ID. Attempting spatial lookup against HUC4 map...\n")
          
          # Transform fire to WGS84 (matching hu4_ref)
          fire_wgs84_temp <- st_transform(fire_row, 4326)
          
          # Find intersection
          sf_use_s2(FALSE)
          intersecting_huc4 <- st_join(fire_wgs84_temp, hu4_ref, join = st_intersects, largest = TRUE)
          sf_use_s2(TRUE)
          
          huc4_target <- intersecting_huc4$huc4[1]
          
          if (is.na(huc4_target)) {
            cat(" - Could not locate this fire within any US HUC4. Skipping.\n")
            next
          }
          cat("   -> Spatially located in HUC4:", huc4_target, "\n")
        }
        
        # --- A. Buffer & Ring ---
        try_proj <- try({ fire_proj <- st_transform(fire_row, 5070) }, silent = TRUE)
        if (inherits(try_proj, "try-error")) next
        
        buffer_proj <- st_buffer(st_geometry(fire_proj), dist = buffer_m)
        fire_union_proj <- st_union(st_geometry(fire_proj))
        
        ring_proj <- tryCatch(st_difference(buffer_proj, fire_union_proj), error = function(e) NULL)
        if (is.null(ring_proj) || length(ring_proj) == 0) next
        
        ring_wgs84 <- st_transform(ring_proj, 4326)
        ring_sf <- st_sf(geometry = st_geometry(ring_wgs84))
        st_crs(ring_sf) <- 4326
        
        # --- B. Download HUCs ---
        huc12s <- tryCatch(get_huc(AOI = ring_sf, type = "huc12"), error = function(e) NULL)
        
        if (is.null(huc12s) || nrow(huc12s) == 0) {
          cat(" - no HUC12s returned for buffer.\n")
          next
        }
        huc12s <- st_make_valid(huc12s)
        
        # --- C. Filter 1: Remove HUCs touching fire ---
        fire_wgs84 <- st_transform(fire_proj, 4326)
        sf_use_s2(FALSE) 
        in_ring_idx <- st_intersects(huc12s, ring_sf, sparse = FALSE)[,1]
        if (!any(in_ring_idx)) { sf_use_s2(TRUE); next }
        
        huc12_in_ring <- huc12s[in_ring_idx, ]
        
        intersects_fire <- st_intersects(huc12_in_ring, fire_wgs84, sparse = FALSE)[,1]
        if (any(intersects_fire)) { huc12_in_ring <- huc12_in_ring[!intersects_fire, ] }
        sf_use_s2(TRUE)
        
        if (nrow(huc12_in_ring) == 0) next
        
        # --- D. Filter 2: KOPPEN MATCHING ---
        huc_vect <- terra::vect(huc12_in_ring)
        huc_vals_df <- terra::extract(koppen_rast, huc_vect, fun="modal", na.rm=TRUE)
        
        if (nrow(huc_vals_df) == 0) next
        
        val_col <- names(huc_vals_df)[2] 
        huc_vals_df$Koppen_Text <- sapply(huc_vals_df[[val_col]], function(x) {
          if(is.na(x)) return(NA)
          kgc::getZone(x)
        })
        
        huc12_in_ring$Koppen_Text <- huc_vals_df$Koppen_Text
        
        huc12_matched_climate <- huc12_in_ring %>% filter(Koppen_Text %in% fire_k_vec)
        
        if (nrow(huc12_matched_climate) == 0) {
          cat(" - no candidates matched Climate Zone.\n")
          next
        }
        
        # --- E. Filter 3: Historical Fires Check (USING HUC4_TARGET) ---
        # This uses the HUC4 we derived at step 2 (either string or spatial)
        huc12_final_check <- huc12_matched_climate
        
        huc4_fires_sf <- load_huc4_fires(huc4_target, huc4_shp_folder)
        
        if (is.null(huc4_fires_sf) || nrow(huc4_fires_sf) == 0) {
          unburned_huc12s_sf <- huc12_final_check
        } else {
          if (st_crs(huc4_fires_sf)$epsg != 4326) huc4_fires_sf <- st_transform(huc4_fires_sf, 4326)
          huc4_fires_sf <- st_make_valid(huc4_fires_sf)
          
          sf_use_s2(FALSE)
          intersects_any <- st_intersects(huc12_final_check, huc4_fires_sf, sparse = FALSE)
          sf_use_s2(TRUE)
          
          if (is.matrix(intersects_any)) { intersects_any_vec <- apply(intersects_any, 1, any)
          } else { intersects_any_vec <- as.logical(intersects_any) }
          
          unburned_huc12s_sf <- huc12_final_check[!intersects_any_vec, ]
        }
        
        if (nrow(unburned_huc12s_sf) == 0) {
          cat(" - all candidates historically burned.\n")
          next
        }
        
        # --- F. Store IDs ---
        possible_id_fields <- c("huc12", "HUC_12", "HUC12", "huc12ce")
        found_id_field <- intersect(possible_id_fields, names(unburned_huc12s_sf))
        if (length(found_id_field) > 0) {
          unburned_ids <- as.character(unburned_huc12s_sf[[found_id_field[1]]])
        } else { unburned_ids <- as.character(seq_len(nrow(unburned_huc12s_sf))) }
        
        fires_out$unburned_huc[i] <- paste(unburned_ids, collapse = ",")
        cat(" - found", length(unburned_ids), "unburned HUC12s.\n")
        
        if (i %% 20 == 0) st_write(fires_out, out_shp, delete_dsn = TRUE, quiet = TRUE)
      }
      
      # FINAL EXPORT
      st_write(fires_out, out_shp, delete_dsn = TRUE)
      cat("Done!\n")
      


################################################################################
#Step 5: Geting shapes for each unburned huc12
      
      
      library(sf)
      library(dplyr)
      library(tidyr)
      library(nhdplusTools)
      
      # ==============================================================================
      # 1. SETUP PATHS
      # ==============================================================================
      
      # Input: The file you just finished creating in Part 3
      input_gpkg <- "/Users/jackboyle/Downloads/UVA/Code/projects/Wildfire_ET_project/processed_data/#3_Unburned_huc_shape/fires_hucs_with_unburned.gpkg"
      
      # Output: Where the final downloaded shapes will go
      out_shp_path <- "/Users/jackboyle/Downloads/UVA/Code/projects/Wildfire_ET_project/processed_data/#3_Unburned_huc_shape/unburned_huc_polygons.gpkg"
      
      # ==============================================================================
      # 2. PREPARE THE LIST OF HUCs TO DOWNLOAD
      # ==============================================================================
      message("Reading input data...")
      fires_data <- st_read(input_gpkg, quiet = TRUE) %>% 
        st_drop_geometry() # We don't need the fire shapes anymore, just the IDs
      
      # Expand the list:
      # If row 1 has "HUC_A, HUC_B", this makes two rows: (Fire1, HUC_A) and (Fire1, HUC_B)
      huc_list <- fires_data %>%
        select(Event_ID, unburned_huc) %>%
        filter(!is.na(unburned_huc) & unburned_huc != "NONE" & unburned_huc != "") %>%
        separate_rows(unburned_huc, sep = ",") %>%
        mutate(unburned_huc = trimws(unburned_huc)) %>%
        distinct(Event_ID, unburned_huc) # Remove duplicates if any
      
      message(sprintf("Found %d unique Fire-HUC pairs to process.", nrow(huc_list)))
      
      # ==============================================================================
      # 3. RESUME LOGIC (Optional but recommended)
      # ==============================================================================
      # If the output file exists, we check what we have already done so we don't re-download.
      
      completed_hucs <- NULL
      if (file.exists(out_shp_path)) {
        message("🔄 Found existing output file. Checking for finished HUCs...")
        try({
          existing_data <- st_read(out_shp_path, quiet = TRUE)
          # Create a unique key (FireID_HUCID) to check what is done
          completed_keys <- paste0(existing_data$Event_ID, "_", existing_data$huc12)
          current_keys   <- paste0(huc_list$Event_ID, "_", huc_list$unburned_huc)
          
          # Filter list to only what is NOT done
          huc_list <- huc_list[! (current_keys %in% completed_keys), ]
          message(sprintf("Skipping already downloaded. Remaining to download: %d", nrow(huc_list)))
        }, silent = TRUE)
      }
      
      # ==============================================================================
      # 4. DOWNLOAD LOOP
      # ==============================================================================
      
      # We will collect results in a list
      results_list <- list()
      counter <- 0
      
      if (nrow(huc_list) > 0) {
        for (i in seq_len(nrow(huc_list))) {
          
          target_huc <- huc_list$unburned_huc[i]
          target_fire <- huc_list$Event_ID[i]
          
          cat(sprintf("[%d / %d] Getting HUC %s for Fire %s... ", 
                      i, nrow(huc_list), target_huc, target_fire))
          
          # Download Geometry
          huc_poly <- tryCatch(
            get_huc(id = target_huc, type = "huc12"),
            error = function(e) return(NULL)
          )
          
          if (is.null(huc_poly) || nrow(huc_poly) == 0) {
            cat("❌Failed/Empty.\n")
            next
          }
          
          # Success! Clean up the polygon
          # We rename 'huc12' column to be standard and ATTACH the Event_ID
          huc_poly_clean <- huc_poly %>%
            select(geometry) %>% # Drop extra NHD columns to keep file small
            mutate(huc12 = target_huc,
                   Event_ID = target_fire) %>%
            st_make_valid()
          
          results_list[[length(results_list) + 1]] <- huc_poly_clean
          cat("Done.\n")
          
          # --- SAVE PERIODICALLY (Every 50 HUCs) ---
          counter <- counter + 1
          if (counter %% 50 == 0) {
            if (length(results_list) > 0) {
              cat(" Saving batch to disk...\n")
              temp_sf <- do.call(rbind, results_list)
              # Append to existing file (or create new)
              st_write(temp_sf, out_shp_path, append = TRUE, quiet = TRUE)
              # Clear memory
              results_list <- list()
            }
          }
        }
        
        # --- FINAL SAVE FOR LEFTOVERS ---
        if (length(results_list) > 0) {
          cat(" Saving final batch...\n")
          temp_sf <- do.call(rbind, results_list)
          st_write(temp_sf, out_shp_path, append = TRUE, quiet = TRUE)
        }
      }
      
      message(" All downloads complete!")
      
      # ==============================================================================
      # 5. VERIFICATION
      # ==============================================================================
      # Load the final file to make sure it looks right
      final_check <- st_read(out_shp_path)
      print(head(final_check))
      message("Total polygons in file: ", nrow(final_check))
      
      
      
      
