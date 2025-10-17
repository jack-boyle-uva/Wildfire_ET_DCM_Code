# PART 1
#Wildfire identification


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
#Step 1: Load MTBS
#Step 2: Filter Fires from 2000-2024
#Step 3: Filter small fires, only ones larger than 14km^2

#Step 4: Add Climate Zone for each fire
#Step 5: Add forest type details for each fire


################################################################################
#Step 1: Load MTBS
Fire_boundary <- read_sf(
  paste0("/Users/jackboyle/Downloads/UVA/Database/WILDFIRE_WATERQUALITY",
         "/MTBS/mtbs_perimeter_data (1)/mtbs_perims_DD.shp")
)
################################################################################
#Step 2: Filter Fires from 2000-2024
Fire_boundary <- Fire_boundary %>%
  mutate(Ig_Date = as.Date(Ig_Date)) %>%        
  filter(Ig_Date >= as.Date("2000-01-01"))

################################################################################
#Step 3: Filter small fires, only ones larger than 5km^2 = 1235.53 Acres
Fire_boundary <- Fire_boundary %>%
  filter(
    BurnBndAc >= 1235.53                             
  )

################################################################################
#Step 4: Add Climate Zone for each fire
koppen<-rast("/Users/jackboyle/Downloads/UVA/Database/Climate_zones/koppen_geiger_tif/1991_2020/koppen_geiger_0p1.tif")

fires <- Fire_boundary
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
Fire_boundary <- fires %>%
  left_join(summary, by = "Event_ID")
Fire_boundary <- Fire_boundary %>%
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
#Step 6: adding soil classifications to fires





