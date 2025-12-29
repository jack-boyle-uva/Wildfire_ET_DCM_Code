# PART 4

"Step 1: filter out wildfires that don't make any sense"
"Step 2: Find the slope of each DCM, find the K mean, classify into High, Moderate, Immedaite DCM"

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


#IMMPORTANT:
#The radio method does not work... but the ratio method should work
################################################################################
################################################################################
################################################################################
################################################################################
# PART 4: Final Data Preparation & Logic Correction
          #
          # THE LOGIC PROBLEM: "Fractional Canopy vs. Binary Loss"
          # ------------------------------------------------------------------------------
          # 1. NLCD Tree Canopy (Pre-Fire):
          #    This data is "Fractional." If a 30m pixel (900 m2) is 30% covered by trees,
          #    NLCD calculates the canopy area as: 0.30 * 900 = 270 m2.
          #
          # 2. Hansen Global Forest Change (Loss):
          #    This data is often "Binary" or "Stand-Replacing." If that same pixel burns 
          #    severely, Hansen often flags the *entire* pixel as "Loss."
          #    Loss Area = 900 m2.
          #
          # 3. The Resulting Paradox:
          #    Loss (900 m2) > Pre-Fire Canopy (270 m2).
          #    Math: 900 / 270 = 333% Loss. (Physically Impossible).
          #
          # THE SOLUTION: "Corrected Baseline"
          # ------------------------------------------------------------------------------
          # We assume that you cannot lose more vegetation than existed. 
          # If Hansen detected 900 m2 of loss, there must have been 900 m2 of vegetation 
          # there to begin with (likely including understory/shrubs that NLCD didn't count).
          #
          # Logic: Baseline = MAX(NLCD_Canopy, Hansen_Total_Loss)
          # Result: 900 / 900 = 100% Loss. (Physically Valid).

################################################################################
#Remove data from years 4 and 5
    #"The post- re tree mortality pro-
      #cess can continue up to 5 years; however, the majority of tree
      #morality in pine and larch forest occurs within 2 years after are 
      #(Vorontsov 1978, Isaev 1962) and 1–3 years after re in
      #spruce forest (Maslov 2011). Therefore, forest loss patches
      #occurring within 3 years after MODIS-detected re events
      #were considered as stand-replacement burned areas."

#Krylov, A., McCarty, J. L., Potapov, P. V., Loboda, T. V., Tyukavina, A., Turubanova, S. A., & Hansen, M. C. (2014). 
#Remote sensing estimates of stand-replacement fires in Russia, 2002–2011. Environmental Research Letters, 9(10), 
#105007. https://doi.org/10.1088/1748-9326/9/10/105007
################################################################################
################################################################################
# PART 4: DCM Analysis Suite
# Step 1: Data Cleaning & Baseline Correction
################################################################################
           
            library(scales)
      library(tidyverse)
      library(sf)
      
      input_path  <- "/Users/jackboyle/Downloads/UVA/Code/projects/Wildfire_ET_project/processed_data/#6_post_pre_fire_canopy/raster_method_results/raster_method_pre_post_canopy_fire.gpkg"
      output_path <- "/Users/jackboyle/Downloads/UVA/Code/projects/Wildfire_ET_project/processed_data/#6_post_pre_fire_canopy/raster_method_results/DCM_analysis/filterd.gpkg"
      
      Fires <- st_read(input_path, quiet = TRUE)
      
      Fire_Calculated <- Fires %>% 
        dplyr::filter(!is.na(pre_fire_canopy_m2)) %>% 
        dplyr::mutate(across(starts_with("post_fire"), ~replace_na(., 0)),
                      fire_year_loss_m2 = replace_na(fire_year_loss_m2, 0)) %>%
        dplyr::mutate(total_loss_all_years = fire_year_loss_m2 + post_fire_year_1_m2 + 
                        post_fire_year_2_m2 + post_fire_year_3_m2 + 
                        post_fire_year_4_m2 + post_fire_year_5_m2) %>%
        dplyr::mutate(
          corrected_baseline_m2 = pmax(pre_fire_canopy_m2, total_loss_all_years),
          data_mismatch_flag = ifelse(total_loss_all_years > pre_fire_canopy_m2, "Impossible (Fixed)", "Possible")
        ) %>%
        dplyr::filter(fire_year_loss_m2 > 0) %>%
        dplyr::mutate(
          pct_loss_Y0 = (fire_year_loss_m2   / corrected_baseline_m2) * 100,
          pct_loss_Y1 = (post_fire_year_1_m2 / corrected_baseline_m2) * 100,
          pct_loss_Y2 = (post_fire_year_2_m2 / corrected_baseline_m2) * 100,
          pct_loss_Y3 = (post_fire_year_3_m2 / corrected_baseline_m2) * 100,
          pct_loss_Y4 = (post_fire_year_4_m2 / corrected_baseline_m2) * 100,
          pct_loss_Y5 = (post_fire_year_5_m2 / corrected_baseline_m2) * 100
        )
      
      st_write(Fire_Calculated, output_path, delete_dsn = TRUE)

################################################################################
# Step 2: Classifying DCM by Shape (Integrated 10% Severity Filter)
################################################################################
        raw_data <- st_read(output_path, quiet = TRUE) %>% st_drop_geometry()
        
        shape_prep <- raw_data %>%
          # A. Calculate total loss in the 3-year attribution window
          dplyr::mutate(total_loss_truncated = pct_loss_Y0 + pct_loss_Y1 + pct_loss_Y2 + pct_loss_Y3) %>% 
          
          # --- THE 10% SEVERITY GATE ---
          # We filter here so the K-means algorithm only groups 'Significant' fires
          dplyr::filter(total_loss_truncated >= 10) %>% 
          
          # B. NORMALIZATION (P0-P3 sum to 1.0)
          dplyr::mutate(
            P0 = pct_loss_Y0 / total_loss_truncated,
            P1 = pct_loss_Y1 / total_loss_truncated,
            P2 = pct_loss_Y2 / total_loss_truncated,
            P3 = pct_loss_Y3 / total_loss_truncated
          ) %>%
          
          # C. SOFT CONTINUITY FILTER
          dplyr::filter(!(P1 == 0 & P2 == 0 & P3 > 0.40))
        
        # D. K-MEANS CLUSTERING
        set.seed(42)
        shape_matrix <- shape_prep %>% dplyr::select(P0, P1, P2, P3)
        k_result     <- kmeans(shape_matrix, centers = 3, nstart = 25)
        shape_prep$Cluster <- as.factor(k_result$cluster)
        
        # E. NAME THE CLASSES
        cluster_logic <- shape_prep %>%
          dplyr::group_by(Cluster) %>%
          dplyr::summarize(Avg_Delayed_Ratio = mean(P1 + P2 + P3), .groups = "drop") %>%
          dplyr::arrange(Avg_Delayed_Ratio) %>%
          dplyr::mutate(DCM_Class = factor(c("Immediate (Low Delay)", "Moderate Delayed", "High Delayed"),
                                           levels = c("Immediate (Low Delay)", "Moderate Delayed", "High Delayed")))
        
        final_data <- shape_prep %>% 
          dplyr::inner_join(cluster_logic %>% dplyr::select(Cluster, DCM_Class), by = "Cluster")
        
        # F. SAVE THE CLASSIFIED GPKG
        final_spatial <- st_read(output_path, quiet = TRUE) %>%
          dplyr::inner_join(final_data %>% dplyr::select(Event_ID, DCM_Class, P0, P1, P2, P3), by = "Event_ID")
        
        st_write(final_spatial, 
                 "/Users/jackboyle/Downloads/UVA/Code/projects/Wildfire_ET_project/processed_data/#6_post_pre_fire_canopy/raster_method_results/DCM_analysis/DCM_Classified_ShapeView.gpkg", 
                 delete_dsn = TRUE)

################################################################################
################################################################################
#End of PART 4
    #Output:  /Users/jackboyle/Downloads/UVA/Code/projects/Wildfire_ET_project/processed_data/
              #6_post_pre_fire_canopy/raster_method_results/DCM_analysis/DCM_Classified_LogView.gpkg    
          
          #Next up: ET and P for the wildfire areas, and their reference catchments
