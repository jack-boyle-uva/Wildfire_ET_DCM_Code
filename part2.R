#Step 2
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
usa <- ne_countries(country = "United States of America", scale = "medium", returnclass = "sf")
plot(usa["geometry"])
hu4<-get_huc(AOI = usa, id = NULL, t_srs = NULL, buffer = 0.5, type = "huc04")
hu4 <- hu4 %>%
  filter(!states %in% c("AK", "HI"))
hu4 <- st_make_valid(hu4)

Fire_boundary <- read_sf("/Users/jackboyle/Downloads/UVA/Database/WILDFIRE_WATERQUALITY",
         "/MTBS/mtbs_perimeter_data (1)/mtbs_perims_DD.shp")

Fire_boundary <- Fire_boundary %>%
  mutate(Ig_Date = as.Date(Ig_Date)) %>%
  filter(format(Ig_Date, "%Y") >= 1950)


#############################################################
Fire_boundary <- st_transform(Fire_boundary, st_crs(hu4)) %>%
  st_make_valid()

# Create unique fire ID if needed
if (!"fire_id" %in% names(Fire_boundary)) {
  Fire_boundary <- Fire_boundary %>%
    mutate(fire_id = row_number())
}

# Make a working copy for looping
fires_remaining <- Fire_boundary

# Output folder
outdir <- "/Users/jackboyle/Downloads/UVA/Code/project"

# HUC4 ID column name in hu4
huc4_id_col <- "huc4"  # change if your column has a different name

# Loop through HUC4s
for (i in 1:nrow(hu4)) {
  this_huc <- hu4[i, ]
  huc_id <- this_huc[[huc4_id_col]]
  
  # Select fires that intersect this HUC4
  fires_in_huc <- st_filter(fires_remaining, this_huc, .predicate = st_intersects)
  
  if (nrow(fires_in_huc) > 0) {
    # Export shapefile
    outfile <- file.path(outdir, paste0("HUC4_", huc_id, "_fires.shp"))
    st_write(fires_in_huc, outfile, delete_layer = TRUE)
    
    # Remove those fires from the working copy (not from Fire_boundary)
    fires_remaining <- fires_remaining %>%
      filter(!fire_id %in% fires_in_huc$fire_id)
  }
  
  cat("Finished HUC4:", huc_id, " | Remaining fires:", nrow(fires_remaining), "\n")
}

################################################################################
# Step 2 
        #Getting burned and unbunred huc12s for each fire

#use the filterd fires dataset
filtered_fires<-st_read("/Users/jackboyle/Downloads/UVA/Code/project/filtered_fires.shp")
results <- list()
for (i in seq_len(nrow(filtered_fires))) {
  fire <- filtered_fires[i, ]
  fire_id <- fire[["Evnt_ID"]]   # unique fire ID from your dataset
  huc12 <- tryCatch(
    get_huc(AOI = st_transform(fire, 4326), type = "huc12"),
    error = function(e) NULL
  )
  if (is.null(huc12) || nrow(huc12) == 0) {
    results[[i]] <- st_sf(
      Evnt_ID    = fire_id,
      burned_huc = NA_character_,
      geometry   = fire$geometry
    )
    next
  }
  inter <- tryCatch(st_intersection(huc12, fire), error = function(e) NULL)
  if (is.null(inter) || nrow(inter) == 0) {
    results[[i]] <- st_sf(
      Evnt_ID    = fire_id,
      burned_huc = NA_character_,
      geometry   = fire$geometry
    )
    next
  }
  inter$overlap_area <- st_area(inter)
  inter$huc_area     <- st_area(huc12)[match(inter$huc12, huc12$huc12)]
  inter$prop_overlap <- as.numeric(inter$overlap_area / inter$huc_area)
  burned <- inter[inter$prop_overlap >= 0.85, ]
  if (nrow(burned) == 0) {
    results[[i]] <- st_sf(
      Evnt_ID    = fire_id,
      burned_huc = NA_character_,
      geometry   = fire$geometry
    )
    next
  }
  burned_ids <- paste(burned$huc12, collapse = ",")
  results[[i]] <- st_sf(
    Evnt_ID    = fire_id,
    burned_huc = burned_ids,
    geometry   = fire$geometry
  )
  cat("Finished fire", i, "/", nrow(filtered_fires), "\n")
}
fires_hucs <- do.call(rbind, results)
st_write(fires_hucs, "/Users/jackboyle/Downloads/UVA/Code/projectfires_hucs.shp")
#look at how many remaining fires there are aftere getting ride odf the na values for hucs
fires_hucs_clean<-fires_hucs%>%filter(!is.na(burned_huc))

################################################################################
# Step 3 
#Getting unbunred huc12s for each fire
#1. make a 20000 meter buffer around the wildfire





library(stringr)
library(nhdplusTools)

# Parameters
buffer_m <- 20000                          # 20,000 m buffer
huc4_shp_folder <- "/Users/jackboyle/Downloads/UVA/Code/project"
out_shp <- "/Users/jackboyle/Downloads/UVA/Code/project/fires_hucs_with_unburned.shp"

# Ensure fires_hucs_clean exists
if (!exists("fires_hucs_clean")) stop("fires_hucs_clean not found in environment")

# Make a copy to modify
fires_out <- fires_hucs_clean

# Make sure geometry is valid
fires_out <- st_make_valid(fires_out)

# Ensure the known ID & burned column names exist; adjust if necessary
if (!("Evnt_ID" %in% names(fires_out))) {
  # try alternative common name
  if ("Evnt_ID" %in% names(fires_out)) {
    # nothing to do
  } else {
    # fallback to first column that looks like an ID
    warning("Evnt_ID not present. Using first column as ID. Adjust code if wrong.")
  }
}

# Function to load HUC4 fire files for a vector of HUC4 codes
load_huc4_fires <- function(huc4_codes, folder) {
  # returns a single sf object (combined) or NULL if none found
  files_to_load <- unlist(lapply(huc4_codes, function(h4) {
    # pattern: HUC4_0420_fires.shp (exact)
    pat <- paste0("^HUC4_", h4, "_fires\\.shp$")
    list.files(folder, pattern = pat, full.names = TRUE)
  }))
  files_to_load <- unique(files_to_load)
  if (length(files_to_load) == 0) return(NULL)
  # read and combine
  shp_list <- lapply(files_to_load, function(f) {
    tryCatch(st_read(f, quiet = TRUE), error = function(e) NULL)
  })
  shp_list <- Filter(Negate(is.null), shp_list)
  if (length(shp_list) == 0) return(NULL)
  out <- do.call(rbind, shp_list)
  return(st_make_valid(out))
}

# Loop through each fire row and compute unburned_huc
unburned_list <- vector("list", nrow(fires_out))

for (i in seq_len(nrow(fires_out))) {
  cat("Processing fire", i, "/", nrow(fires_out), "...\n")
  fire_row <- fires_out[i, ]
  # Unique ID for messages & output
  fire_id <- as.character(fire_row[["Evnt_ID"]])
  
  # Default: NA result (string)
  unburned_list[[i]] <- NA_character_
  
  # If burned_huc is NA or empty, just continue with NA
  burned_cell <- fire_row[["burned_huc"]]
  if (is.na(burned_cell) || is.null(burned_cell) || nchar(as.character(burned_cell)) == 0) {
    cat(" - burned_huc is NA/empty for fire", fire_id, "- skipping.\n")
    next
  }
  
  # 1) Make 20 km buffer around fire perimeter (project to EPSG:5070 for meters)
  # We'll project both the fire geometry and the buffer operations into EPSG:5070
  try({
    fire_proj  <- st_transform(fire_row, 5070)  # projected for accurate buffering
  }, silent = TRUE)
  if (!exists("fire_proj") || is.null(fire_proj)) {
    cat(" - failed to transform fire geometry for fire", fire_id, "- skipping.\n")
    next
  }
  
  # create buffer and ring (buffer minus fire)
  buffer_proj <- st_buffer(st_geometry(fire_proj), dist = buffer_m)
  # ensure single polygon for ring subtraction
  fire_union_proj <- st_union(st_geometry(fire_proj))
  ring_proj <- tryCatch(st_difference(buffer_proj, fire_union_proj),
                        error = function(e) NULL)
  if (is.null(ring_proj) || length(ring_proj) == 0) {
    cat(" - could not create buffer ring for fire", fire_id, "- skipping.\n")
    next
  }
  
  # convert ring back to WGS84 for get_huc (get_huc expects WGS84)
  ring_wgs84 <- st_transform(ring_proj, 4326)
  # turn into an sf object with a geometry column so get_huc will accept it
  ring_sf <- st_sf(geometry = st_geometry(ring_wgs84))
  st_crs(ring_sf) <- 4326
  
  # 2) Use get_huc to pull HUC12s that intersect the buffer area
  huc12s <- tryCatch(
    get_huc(AOI = ring_sf, type = "huc12"),
    error = function(e) NULL
  )
  if (is.null(huc12s) || nrow(huc12s) == 0) {
    cat(" - no HUC12s returned for buffer for fire", fire_id, "\n")
    next
  }
  huc12s <- st_make_valid(huc12s)
  
  # 3) Keep only HUC12s that actually lie in the ring (i.e., intersect ring but NOT intersect fire)
  # transform all to same CRS (use WGS84 returned by get_huc)
  fire_wgs84 <- st_transform(fire_proj, 4326)
  
  in_ring_idx <- st_intersects(huc12s, ring_sf, sparse = FALSE)[,1]
  if (!any(in_ring_idx)) {
    cat(" - no HUC12s in the buffer ring for fire", fire_id, "\n")
    next
  }
  huc12_in_ring <- huc12s[in_ring_idx, ]
  
  # remove any HUC12s that also intersect the fire polygon (we don't want HUCs inside the fire)
  intersects_fire <- st_intersects(huc12_in_ring, fire_wgs84, sparse = FALSE)[,1]
  if (any(intersects_fire)) {
    huc12_in_ring <- huc12_in_ring[!intersects_fire, ]
  }
  if (nrow(huc12_in_ring) == 0) {
    cat(" - no HUC12s between buffer and fire (after removing those touching the fire) for", fire_id, "\n")
    next
  }
  
  # 4) From burned_huc cell, extract HUC4 codes (first 4 characters of each HUC12 stored there)
  # burned_huc may contain multiple HUC12s separated by commas
  burned_huc_str <- as.character(burned_cell)
  burned_huc_vec <- str_split(burned_huc_str, pattern = ",", simplify = FALSE)[[1]]
  burned_huc_vec <- str_trim(burned_huc_vec)
  # Extract first 4 digits for each; ensure they are 4 characters (pad if needed)
  huc4_from_burned <- unique(substr(burned_huc_vec, 1, 4))
  huc4_from_burned <- huc4_from_burned[huc4_from_burned != "" & !is.na(huc4_from_burned)]
  if (length(huc4_from_burned) == 0) {
    cat(" - no valid HUC4 codes extracted from burned_huc for", fire_id, "\n")
    next
  }
  
  # 5) Load the HUC4 fire-history shapefile(s) for those HUC4s (from folder)
  huc4_fires_sf <- load_huc4_fires(huc4_from_burned, huc4_shp_folder)
  if (is.null(huc4_fires_sf) || nrow(huc4_fires_sf) == 0) {
    cat(" - no HUC4 fire-history shapefiles found for HUC4(s):", paste(huc4_from_burned, collapse = ","), "for fire", fire_id, "\n")
    # we can't determine reference unburned HUC12s if no HUC4 file is present — skip with NA
    next
  }
  
  # Ensure the HUC4 fires are in WGS84 (so intersections align)
  if (st_crs(huc4_fires_sf)$epsg != 4326) {
    huc4_fires_sf <- st_transform(huc4_fires_sf, 4326)
  }
  huc4_fires_sf <- st_make_valid(huc4_fires_sf)
  
  # 6) For each HUC12 in the ring, check whether it intersects any wildfire polygon in the combined HUC4 fires layer
  # If a HUC12 does NOT intersect any wildfire polygon within its HUC4 fire-history, then it's an "unburned reference"
  # We'll test intersection with the relevant HUC4 fires (loaded combined)
  # To speed up, we use st_intersects matrix
  intersects_any <- st_intersects(huc12_in_ring, huc4_fires_sf, sparse = FALSE)
  # intersects_any will be a matrix; reduce to logical vector indicating whether each huc12 intersects ANY huc4 fire poly
  if (is.matrix(intersects_any)) {
    intersects_any_vec <- apply(intersects_any, 1, any)
  } else {
    # if single column, ensure logical
    intersects_any_vec <- as.logical(intersects_any)
  }
  
  # Keep HUC12s that did NOT intersect any huc4 fire history polygons
  unburned_huc12s_sf <- huc12_in_ring[!intersects_any_vec, ]
  
  if (nrow(unburned_huc12s_sf) == 0) {
    cat(" - no unburned reference HUC12s found for fire", fire_id, "\n")
    unburned_list[[i]] <- NA_character_
    next
  }
  
  # Extract their HUC12 IDs (assume the field is named 'huc12' or 'HUC12' or similar)
  possible_id_fields <- c("huc12", "HUC_12", "HUC12", "huc12ce")  # try common variants
  found_id_field <- intersect(possible_id_fields, names(unburned_huc12s_sf))
  if (length(found_id_field) == 0) {
    # try to find a character column that looks like a huc12 (length 12)
    char_cols <- names(unburned_huc12s_sf)[sapply(unburned_huc12s_sf, is.character)]
    huc_guess <- NULL
    for (cc in char_cols) {
      vals <- unburned_huc12s_sf[[cc]]
      if (any(!is.na(vals) & nchar(vals) == 12)) { huc_guess <- cc; break }
    }
    if (is.null(huc_guess)) {
      warning("Couldn't find HUC12 ID column in fetched HUC12s; using row numbers. Check field names.")
      unburned_ids <- as.character(seq_len(nrow(unburned_huc12s_sf)))
    } else {
      unburned_ids <- as.character(unburned_huc12s_sf[[huc_guess]])
    }
  } else {
    unburned_ids <- as.character(unburned_huc12s_sf[[found_id_field[1]]])
  }
  
  # collapse to comma-separated string and assign to output list
  unburned_list[[i]] <- paste(unburned_ids, collapse = ",")
  cat(" - found", length(unburned_ids), "unburned HUC12s for fire", fire_id, "\n")
}

# Attach the new column to fires_out (preserving order)
fires_out$unburned_huc <- unlist(unburned_list)

# Write results to disk (GeoPackage is often nicer, but using shapefile if you prefer)
tryCatch({
  st_write(fires_out, out_shp, delete_dsn = TRUE, quiet = TRUE)
  cat("Wrote output to:", out_shp, "\n")
}, error = function(e) {
  warning("Failed to write shapefile: ", e$message)
})

# Also return fires_out to environment
fires_hucs_with_unburned <- fires_out
cat("Done. 'fires_hucs_with_unburned' created in environment with column 'unburned_huc'.\n")


 
fires_hucs_with_unburned<-fires_hucs_with_unburned%>%filter(!is.na(unburned_huc))
st_write(fires_hucs_with_unburned, "/Users/jackboyle/Downloads/UVA/Code/fires_hucs_with_unburned.shp")

################################################################################
fires_hucs_with_unburned<-st_read("/Users/jackboyle/Downloads/UVA/Code/fires_hucs_with_unburned.shp")


################################################################################
#Now we need shape files of every unburned huc, its huc12 id, and its corespoining Event_ID
huc12_shapes<-read_sf("/Users/jackboyle/Downloads/UVA/Code/project/fires_hucs_with_unburned.shp")



# Assuming you already have `unburned_huc_shapes`
# Step 1: Expand the unbrnd_ column into individual HUC12 rows
expanded_hucs <- huc12_shapes %>%
  st_drop_geometry() %>%
  separate_rows(unbrnd_, sep = ",") %>%
  rename(HUC_12 = unbrnd_) %>%
  distinct(Evnt_ID, HUC_12)

# Step 2: Create an empty list to store results
huc_results <- list()

# Step 3: Loop through each unique HUC12 and download its geometry
for (i in seq_len(nrow(expanded_hucs))) {
  huc_id <- expanded_hucs$HUC_12[i]
  evnt_id <- expanded_hucs$Evnt_ID[i]
  
  cat("Getting HUC12:", huc_id, "for Event:", evnt_id, "\n")
  
  # Try to get the polygon; skip if fails
  huc_poly <- tryCatch(
    nhdplusTools::get_huc(id = huc_id, type = "huc12"),
    error = function(e) NULL
  )
  
  if (!is.null(huc_poly) && nrow(huc_poly) > 0) {
    huc_poly <- huc_poly %>%
      mutate(Evnt_ID = evnt_id, HUC_12 = huc_id) %>%
      select(Evnt_ID, HUC_12, geometry)
    huc_results[[length(huc_results) + 1]] <- huc_poly
  }
}

# Step 4: Combine all polygons into one sf object
unburned_huc_polygons <- do.call(rbind, huc_results)

write_sf(unburned_huc_polygons,"/Users/jackboyle/Downloads/UVA/Code/project/shape/For unburned areas/HUC12_shapefiles_unburn/unburned_huc_polygons.shp")
print(unburned_huc_polygons)