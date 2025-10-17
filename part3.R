#Part 3
#giving burn severities
#add three columns with % of high, moderate and low within the fire boundary.
#then try to find % canopy loss for every year after the wildfire.

#Analysing Evapotrapirtation data
library(soilDB)
library(kgc)
library(sf)
library(dplyr)
library(plyr)
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


Fire_boundary <- read_sf("/Users/jackboyle/Downloads/UVA/Database/WILDFIRE_WATERQUALITY/MTBS/mtbs_perimeter_data (1)/mtbs_perims_DD.shp")

Hucs<-st_read("/Users/jackboyle/Downloads/UVA/Code/fires_hucs_with_unburned.shp")
################################################################################

        
library(terra)
library(sf)
library(dplyr)
library(lubridate)

# Base directory for MTBS mosaics
base_dir <- "/Users/jackboyle/Downloads/UVA/Database/WILDFIRE_WATERQUALITY/MTBS/composite_data/MTBS_BSmosaics"

# Empty list to hold results (we’ll later combine them into an sf object)
burn_s_list <- list()

for (i in seq_len(nrow(Hucs))) {
  event_id <- Hucs$Evnt_ID[i]
  burned_hucs <- strsplit(Hucs$brnd_hc[i], ",")[[1]]
  ig_date <- Fire_boundary %>%
    filter(Event_ID == event_id) %>%
    pull(Ig_Date)
  
  if (length(ig_date) == 0 || is.na(ig_date)) {
    message("Skipping ", event_id, ": no ignition date found.")
    next
  }
  
  year_folder <- year(ig_date)
  mtbs_path <- file.path(
    base_dir,
    paste0(year_folder),
    paste0("mtbs_CONUS_", year_folder),
    paste0("mtbs_CONUS_", year_folder, ".tif")
  )
  
  if (!file.exists(mtbs_path)) {
    message("Skipping ", event_id, ": MTBS file for ", year_folder, " not found at ", mtbs_path)
    next
  }
  
  mtbs_rast <- rast(mtbs_path)
  
  for (huc_id in burned_hucs) {
    huc_id <- trimws(huc_id)
    huc_poly <- tryCatch({
      get_huc(id = huc_id)
    }, error = function(e) {
      message("Failed to get HUC ", huc_id, " for event ", event_id)
      return(NULL)
    })
    if (is.null(huc_poly)) next
    
    # Transform to match raster CRS
    huc_poly <- st_transform(huc_poly, crs(mtbs_rast))
    
    mtbs_crop <- tryCatch({
      terra::mask(terra::crop(mtbs_rast, huc_poly), huc_poly)
    }, error = function(e) {
      message("Raster-HUC mismatch for ", huc_id, " in event ", event_id)
      return(NULL)
    })
    
    if (is.null(mtbs_crop)) next
    
    vals <- terra::values(mtbs_crop)
    vals <- vals[!is.na(vals)]
    if (length(vals) == 0) next
    
    total <- length(vals)
    low  <- sum(vals == 1) / total * 100
    mod  <- sum(vals == 2) / total * 100
    high <- sum(vals == 3) / total * 100
    unb  <- sum(vals == 4) / total * 100
    
    # Create an sf object for this HUC
    burn_entry <- st_sf(
      Event_ID = event_id,
      HUC12 = huc_id,
      Ig_Date = ig_date,
      low_severity_percent = low,
      moderate_severity_percent = mod,
      high_severity_percent = high,
      unburned_percent = unb,
      geometry = huc_poly$geometry
    )
    
    burn_s_list[[length(burn_s_list) + 1]] <- burn_entry
    message("Processed ", event_id, " - HUC ", huc_id)
  }
}

# Combine all HUC results into a single sf object
burn_s_sf <- do.call(rbind, burn_s_list)

# Define output directory
output_dir <- "/Users/jackboyle/Downloads/UVA/Code/project/shape"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Write to shapefile
output_path <- file.path(output_dir, "burn_severity_by_huc.shp")
st_write(burn_s_sf, output_path, delete_layer = TRUE)

message("✅ All done! Saved shapefile to: ", output_path)











##########################################################################
        base_dir <- "/Users/jackboyle/Downloads/UVA/Database/WILDFIRE_WATERQUALITY/MTBS/composite_data/MTBS_BSmosaics"
        burn_s <- data.frame(
          Event_ID = character(),
          HUC12 = character(),
          Ig_Date = as.Date(character()),
          low_severity_percent = numeric(),
          moderate_severity_percent = numeric(),
          high_severity_percent = numeric(),
          unburned_percent = numeric(),
          stringsAsFactors = FALSE
        )
        for (i in seq_len(nrow(Hucs))) {
          event_id <- Hucs$Evnt_ID[i]
          burned_hucs <- strsplit(Hucs$brnd_hc[i], ",")[[1]]
          ig_date <- Fire_boundary %>%
            filter(Event_ID == event_id) %>%
            pull(Ig_Date)
          
          if (length(ig_date) == 0 || is.na(ig_date)) {
            message("Skipping ", event_id, ": no ignition date found.")
            next
          }
          
          year_folder <- year(ig_date)
          mtbs_path <- file.path(
            base_dir,
            paste0(year_folder),
            paste0("mtbs_CONUS_", year_folder),
            paste0("mtbs_CONUS_", year_folder, ".tif")
          )
          
          if (!file.exists(mtbs_path)) {
            message("Skipping ", event_id, ": MTBS file for ", year_folder, " not found at ", mtbs_path)
            next
          }
          mtbs_rast <- rast(mtbs_path)
          for (huc_id in burned_hucs) {
            huc_id <- trimws(huc_id)
            huc_poly <- tryCatch({
              get_huc(id = huc_id)
            }, error = function(e) {
              message("Failed to get HUC ", huc_id, " for event ", event_id)
              return(NULL)
            })
            if (is.null(huc_poly)) next
            huc_poly <- st_transform(huc_poly, crs(mtbs_rast))
            mtbs_crop <- tryCatch({
              terra::mask(terra::crop(mtbs_rast, huc_poly), huc_poly)
            }, error = function(e) {
              message("Raster-HUC mismatch for ", huc_id, " in event ", event_id)
              return(NULL)
            })
            if (is.null(mtbs_crop)) next
            vals <- terra::values(mtbs_crop)
            vals <- vals[!is.na(vals)]
            
            if (length(vals) == 0) next
            
            total <- length(vals)
            low <- sum(vals == 1) / total * 100
            mod <- sum(vals == 2) / total * 100
            high <- sum(vals == 3) / total * 100
            unb <- sum(vals == 4) / total * 100
            
            burn_s <- rbind(
              burn_s,
              data.frame(
                Event_ID = event_id,
                HUC12 = huc_id,
                Ig_Date = ig_date,
                low_severity_percent = low,
                moderate_severity_percent = mod,
                high_severity_percent = high,
                unburned_percent = unb
              )
            )
            message("Processed ", event_id, " - HUC ", huc_id)
          }
        }
        write.csv(burn_s, "/Users/jackboyle/Downloads/UVA/Code/project/shape/burn_severity_by_huc.csv", row.names = FALSE)
        message("✅ All done! Saved to burn_severity_by_huc.csv")
        

################################################################################
# In the next step we might need google earth engine to get ET at the 30x30m landsat resolution
        
################################################################################
       
        #Ajusting the GEE dataset 
   dataset1<-read.csv("/Users/jackboyle/Downloads/SSEBop_Annual_ET_by_HUC_fullrecord_SSEbop.csv")     
        dataset1 <- dataset1 %>%
          select(-c(.geo, system.index, Ig_Date))
        dataset1 <- dataset1 %>%
          rename(Event_ID = Evnt_ID)
        dataset1 <- dataset1 %>%
          left_join(
            st_drop_geometry(Fire_boundary)[, c("Event_ID", "Ig_Date")],
            by = "Event_ID"
          )
        
################################################################################
        #Save that dataset
################################################################################  
    #Make statistical calculations and plots    
        
        
      #ok so I'll need to do ann analysis on the unburned reference hucs
        
        dataset1 <- dataset1 %>%
          rowwise() %>%
          mutate(
            dominant_severity = c("High", "Low", "Moderate", "Mostly Unburned")[which.max(c(hgh_sv_, lw_svr_, mdrt_s_, unbrnd_))]
          ) %>%
          ungroup()
        
        # --- Step 2: Compute year offset relative to fire year ---
        dataset_window <- dataset1 %>%
          mutate(
            fire_year = lubridate::year(lubridate::ymd(Ig_Date)),
            year_offset = year - fire_year
          ) %>%
          filter(year_offset >= -3 & year_offset <= 3)
        
        # --- Step 3: Label before/after ---
        dataset_window <- dataset_window %>%
          mutate(period = case_when(
            year_offset < 0 ~ "Before",
            year_offset > 0 ~ "After",
            TRUE ~ "During"
          ))
        
        # --- Step 4: Average ET for 3 years before & after per HUC and severity ---
        et_summary_3yr <- dataset_window %>%
          filter(period != "During") %>%
          group_by(HUC12, dominant_severity, period) %>%
          summarize(mean_ET = mean(mean, na.rm = TRUE), .groups = "drop") %>%
          pivot_wider(names_from = period, values_from = mean_ET) %>%
          mutate(delta_ET = After - Before)
        
        # --- Step 5: Plot (faceted by severity) ---
        ggplot(et_summary_3yr, aes(x = reorder(HUC12, delta_ET), y = delta_ET)) +
          geom_col(aes(fill = delta_ET), width = 0.8) +
          scale_fill_gradient2(
            low = "red", mid = "white", high = "blue",
            midpoint = 0,
            name = "ΔET (mm)\n(After - Before)"
          ) +
          geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
          labs(
            title = "Post-Fire Change in Annual Evapotranspiration (ET)",
            subtitle = paste(
              "Mean of 3 years before vs. 3 years after each fire",
              "Grouped by Dominant Burn Severity (HUC12 level)",
              sep = "\n"
            ),
            x = "HUC12 Watershed",
            y = "Change in ET (mm)",
            caption = paste0("Overall mean change = ",
                             round(mean(et_summary_3yr$delta_ET, na.rm = TRUE), 1),
                             " mm")
          ) +
          facet_wrap(~dominant_severity, scales = "free_x") +
          theme_minimal(base_size = 14) +
          theme(
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            plot.title = element_text(face = "bold"),
            legend.position = "right"
          )
        
        
################################################################################        
     #analysis the unburned    
        
unburned_huc12<-read_sf("/Users/jackboyle/Downloads/UVA/Code/project/data/annual_et_by_huc12._unbrn.csv") 
       
        unburned_huc12 <- unburned_huc12 %>%
          select(-any_of(c(".geo", "system:index")))
        
        Fire_boundary <- Fire_boundary %>%
          rename(Evnt_ID = Event_ID)
        unburned_huc12 <- unburned_huc12 %>%
          left_join(Fire_boundary %>% select(Evnt_ID, Ig_Date), by = "Evnt_ID")
        clean_huc <- unburned_huc12 %>%
          filter(mean != "" & !is.na(mean)) %>%         # remove empty ET values
          mutate(
            mean = as.numeric(mean),
            year = as.numeric(year),
            Ig_Date = as.Date(Ig_Date)
          ) %>%
          mutate(Fire_Year = as.numeric(format(Ig_Date, "%Y")))
        
        et_summary <- clean_huc %>%
          group_by(Evnt_ID, HUC_12, Fire_Year) %>%
          summarize(
            pre_fire_mean = mean(mean[year >= (Fire_Year - 3) & year < Fire_Year], na.rm = TRUE),
            post_fire_mean = mean(mean[year > Fire_Year & year <= (Fire_Year + 3)], na.rm = TRUE),
            .groups = "drop"
          ) %>%
          mutate(change_ET = post_fire_mean - pre_fire_mean)
        print(head(et_summary))
        summary(et_summary$change_ET)
        
       
        
        
        unburned_huc12_clean <- unburned_huc12 %>%
          mutate(mean = as.numeric(mean)) %>%       # convert to numeric
          filter(!is.na(mean))                      # drop empty mean values
        
        # Summarize by Evnt_ID
        et_summary_event <- unburned_huc12_clean %>%
          group_by(Evnt_ID, Ig_Date) %>%
          summarize(mean_ET = mean(mean, na.rm = TRUE), .groups = "drop")
        et_summary_event <- unburned_huc12_clean %>%
          mutate(
            year = as.numeric(year),
            mean = as.numeric(mean)
          ) %>%
          filter(!is.na(mean)) %>%
          group_by(Evnt_ID, Ig_Date) %>%
          summarize(
            mean_pre = mean(mean[year >= (year(Ig_Date) - 3) & year < year(Ig_Date)], na.rm = TRUE),
            mean_post = mean(mean[year > year(Ig_Date) & year <= (year(Ig_Date) + 3)], na.rm = TRUE),
            change_ET = mean_post - mean_pre,
            .groups = "drop"
          )
        
        
       
  #################The goofd plot##      
        
        
        library(dplyr)
        library(ggplot2)
        
        # Clean and summarize ET by fire event (Evnt_ID)
        clean_huc <- unburned_huc12 %>%
          mutate(
            mean = as.numeric(mean),
            year = as.numeric(year),
            Ig_Date = as.Date(Ig_Date)
          ) %>%
          filter(!is.na(mean)) %>%
          mutate(Fire_Year = as.numeric(format(Ig_Date, "%Y")))
        
        # Summarize to get 3-year pre/post fire means by fire event
        et_summary_event <- clean_huc %>%
          group_by(Evnt_ID, Fire_Year) %>%
          summarize(
            pre_fire_mean = mean(mean[year >= (Fire_Year - 3) & year < Fire_Year], na.rm = TRUE),
            post_fire_mean = mean(mean[year > Fire_Year & year <= (Fire_Year + 3)], na.rm = TRUE),
            .groups = "drop"
          ) %>%
          mutate(change_ET = post_fire_mean - pre_fire_mean) %>%
          filter(!is.na(change_ET))
        
        # Optional cleaning
        et_summary_event_clean <- et_summary_event %>%
          filter(is.finite(change_ET))
        
        # Summary stats (for caption)
        overall_mean_change <- round(mean(et_summary_event_clean$change_ET, na.rm = TRUE), 1)
        n_events <- nrow(et_summary_event_clean)
        
        # --- Plot ---
        ggplot(et_summary_event_clean, aes(x = reorder(Evnt_ID, change_ET), y = change_ET)) +
          geom_col(aes(fill = change_ET), width = 0.8) +
          scale_fill_gradient2(
            low = "red", mid = "white", high = "blue",
            midpoint = 0,
            name = "ΔET (mm)\n(Post - Pre)"
          ) +
          geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
          labs(
            title = "Change in Annual Evapotranspiration (ET) by Wildfire Event",
            subtitle = "3-Year Mean After Fire Minus 3-Year Mean Before Fire",
            x = "Wildfire Event (Evnt_ID)",
            y = "ΔET (mm)",
            caption = paste0(
              "Overall mean ΔET = ", overall_mean_change,
              " mm (", n_events, " fires)"
            )
          ) +
          theme_minimal(base_size = 14) +
          theme(
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            plot.title = element_text(face = "bold"),
            legend.position = "right"
          )
        