# ============================================================
# Wintermoth location, VIIRS LSP, and dominant CLC per pixel
# ============================================================

library(terra)
library(sf)
library(dplyr)
install.packages("exactextractr")
library(exactextractr)

# ------------------------------------------------------------
# 1. Load & deduplicate locations
# ------------------------------------------------------------
wm_location <- read.csv("../../data/wgs_europe_metadata.csv") %>%
  distinct(pop, latitude, longitude, .keep_all = TRUE)

sin_crs <- "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"

# ------------------------------------------------------------
# 2. Build 5 x 5 km square buffers
# ------------------------------------------------------------
location_sf  <- st_as_sf(wm_location, coords = c("longitude", "latitude"), crs = 4326)
location_sin <- st_transform(location_sf, crs = sin_crs)
coords       <- st_coordinates(location_sin)
half         <- 2500

squares <- mapply(function(x, y) {
  st_polygon(list(matrix(c(
    x - half, y - half,
    x + half, y - half,
    x + half, y + half,
    x - half, y + half,
    x - half, y - half
  ), ncol = 2, byrow = TRUE)))
}, coords[, 1], coords[, 2], SIMPLIFY = FALSE)

square_buffers <- st_sf(
  st_drop_geometry(location_sin),
  geometry = st_sfc(squares, crs = sin_crs)
)

buf_vect <- vect(square_buffers)

# ------------------------------------------------------------
# 3. Load CLC GeoTIFF once
#    Replace with your own file path
# ------------------------------------------------------------
clc_rast <- rast("../../data/CLC/data/U2018_CLC2018_V2020_20u1.tif")

# If CLC is not already in sinusoidal CRS, do not project the whole thing.
# We will crop first, then project each tile-sized subset.
clc_crs <- crs(clc_rast)

# Optional: if your TIFF has multiple layers, keep the first one
if (nlyr(clc_rast) > 1) {
  clc_rast <- clc_rast[[1]]
}

terra::levels(clc_rast)
terra::cats(clc_rast)
terra::freq(clc_rast, bylayer = FALSE)



# ------------------------------------------------------------
# 4. Load VIIRS tiles for target year
# ------------------------------------------------------------
files <- list.files("../../data/VIIRS", pattern = "VNP22Q2.*\\.h5$", full.names = TRUE)
years <- as.integer(substr(basename(files), 10, 13))

target_year <- 2013
year_files  <- files[years == target_year]

extract_vnp <- function(filepath, year) {
  tile_size <- 1111950.5196
  tile_str  <- regmatches(basename(filepath), regexpr("h\\d{2}v\\d{2}", basename(filepath)))
  tile_h    <- as.integer(substr(tile_str, 2, 3))
  tile_v    <- as.integer(substr(tile_str, 5, 6))
  
  xmin <- (tile_h - 18) * tile_size
  xmax <- xmin + tile_size
  ymax <- (9 - tile_v) * tile_size
  ymin <- ymax - tile_size
  tile_ext <- ext(xmin, xmax, ymin, ymax)
  
  onset_path <- paste0('HDF5:"', filepath, '"://HDFEOS/GRIDS/Cycle_1/Data_Fields/Onset_Greenness_Increase_1')
  qc_path    <- paste0('HDF5:"', filepath, '"://HDFEOS/GRIDS/Cycle_1/Data_Fields/GLSP_QC_1')
  pgq_path   <- paste0('HDF5:"', filepath, '"://HDFEOS/GRIDS/Cycle_1/Data_Fields/PGQ_Onset_Greenness_Increase_1')
  
  onset <- flip(rast(onset_path), direction = "vertical")
  qc    <- flip(rast(qc_path),    direction = "vertical")
  pgq   <- flip(rast(pgq_path),   direction = "vertical")
  
  ext(onset) <- tile_ext; crs(onset) <- sin_crs
  ext(qc)    <- tile_ext; crs(qc)    <- sin_crs
  ext(pgq)   <- tile_ext; crs(pgq)   <- sin_crs
  
  offset          <- (year - 2000) * 366
  onset_corrected <- app(onset, function(x) ifelse(x == 32767, NA, x - offset))
  qc_clean        <- app(qc,    function(x) ifelse(x == 255,   NA, x))
  pgq_clean       <- app(pgq,   function(x) ifelse(x == 255,   NA, x))
  
  ext(onset_corrected) <- tile_ext; crs(onset_corrected) <- sin_crs
  ext(qc_clean)        <- tile_ext; crs(qc_clean)        <- sin_crs
  ext(pgq_clean)       <- tile_ext; crs(pgq_clean)       <- sin_crs
  
  names(onset_corrected) <- "Onset_Greenness_Increase"
  names(qc_clean)        <- "GLSP_QC"
  names(pgq_clean)       <- "PGQ_Onset_Greenness_Increase"
  
  list(raster = c(onset_corrected, qc_clean, pgq_clean), tile_id = tile_str)
}

tile_results <- lapply(year_files, extract_vnp, year = target_year)

# ------------------------------------------------------------
# 5. Helpers for dominant CLC class
# ------------------------------------------------------------
weighted_mode <- function(values, coverage_fraction) {
  ok <- !is.na(values) & !is.na(coverage_fraction)
  values <- values[ok]
  coverage_fraction <- coverage_fraction[ok]
  
  if (length(values) == 0) return(NA_integer_)
  
  w <- tapply(coverage_fraction, values, sum)
  as.integer(names(w)[which.max(w)])
}

make_square <- function(x, y, half_x, half_y) {
  st_polygon(list(matrix(c(
    x - half_x, y - half_y,
    x + half_x, y - half_y,
    x + half_x, y + half_y,
    x - half_x, y + half_y,
    x - half_x, y - half_y
  ), ncol = 2, byrow = TRUE)))
}

# ------------------------------------------------------------
# 6. Extract VIIRS pixels and dominant CLC per pixel
# ------------------------------------------------------------
all_extracts <- lapply(tile_results, function(tile) {
  
  r       <- tile$raster
  tile_id <- tile$tile_id
  
  # only buffers intersecting this VIIRS tile
  tile_poly     <- as.polygons(ext(r), crs = sin_crs)
  buf_intersect <- relate(buf_vect, tile_poly, relation = "intersects")
  buf_subset    <- buf_vect[buf_intersect, ]
  
  if (nrow(buf_subset) == 0) return(NULL)
  
  # crop VIIRS to intersecting buffers
  r_crop <- crop(r, buf_subset)
  
  # crop CLC first in its native CRS, then project only the tile-sized subset
  tile_poly_clc <- project(tile_poly, clc_crs)
  clc_crop_nat  <- crop(clc_rast, tile_poly_clc)
  
  # reproject only the cropped subset
  clc_crop_sin  <- project(clc_crop_nat, sin_crs, method = "near")
  
  # existing VIIRS extraction
  extracted <- extract(r_crop, buf_subset, fun = NULL, ID = TRUE, cells = TRUE, xy = TRUE)
  
  if (is.null(extracted) || nrow(extracted) == 0) return(NULL)
  
  buf_df <- as.data.frame(buf_subset)[, c("pop")]
  extracted$pop     <- buf_df[extracted$ID]
  extracted$tile_id <- tile_id
  extracted$year    <- target_year
  names(extracted)[names(extracted) == "cell"] <- "pixel_id"
  
  # one polygon per unique VIIRS pixel
  pixel_tbl <- extracted %>%
    distinct(pixel_id, x, y)
  
  hx <- res(r_crop)[1] / 2
  hy <- res(r_crop)[2] / 2
  
  pixel_sfc <- st_sfc(
    lapply(seq_len(nrow(pixel_tbl)), function(i) {
      make_square(pixel_tbl$x[i], pixel_tbl$y[i], hx, hy)
    }),
    crs = sin_crs
  )
  
  pixel_sf <- st_sf(pixel_tbl, geometry = pixel_sfc)
  
  # dominant CLC code per VIIRS pixel
  pixel_sf$clc_code <- exact_extract(
    clc_crop_sin,
    pixel_sf,
    weighted_mode
  )
  
  extracted <- extracted %>%
    left_join(
      st_drop_geometry(pixel_sf) %>% select(pixel_id, clc_code),
      by = "pixel_id"
    ) %>%
    select(-ID)
  
  extracted
})

# ------------------------------------------------------------
# 7. Combine
# ------------------------------------------------------------
final_long <- bind_rows(all_extracts) %>%
  select(
    pop, tile_id, year, pixel_id, x, y,
    Onset_Greenness_Increase, GLSP_QC, PGQ_Onset_Greenness_Increase,
    clc_code
  ) %>%
  arrange(pop, tile_id, pixel_id)

head(final_long)

# ------------------------------------------------------------
# 8. QC bit-unpacking function based on VNP22 User Guide
# ------------------------------------------------------------
decode_glsp_qc <- function(qc) {
  
  mandatory_raw <- bitwAnd(qc, 0x03)
  mandatory <- dplyr::case_when(
    mandatory_raw == 0 ~ "processed, good quality",
    mandatory_raw == 1 ~ "processed, other quality",
    mandatory_raw == 2 ~ "processed, backup algorithm",
    mandatory_raw == 3 ~ "not processed",
    TRUE ~ NA_character_
  )
  
  clim_raw <- bitwShiftR(bitwAnd(qc, 0x1C), 2)
  climatology <- dplyr::case_when(
    clim_raw == 0 ~ "climatology fill 0%",
    clim_raw == 1 ~ "climatology fill 1-5%",
    clim_raw == 2 ~ "climatology fill 6-10%",
    clim_raw == 3 ~ "climatology fill 11-15%",
    clim_raw == 4 ~ "climatology fill 16-20%",
    clim_raw == 5 ~ "climatology fill 21-30%",
    clim_raw == 6 ~ "climatology fill >31%",
    clim_raw == 7 ~ "climatology based detection",
    TRUE ~ NA_character_
  )
  
  lw_raw <- bitwShiftR(bitwAnd(qc, 0xE0), 5)
  land_water <- dplyr::case_when(
    lw_raw == 0 ~ "shallow ocean",
    lw_raw == 1 ~ "land",
    lw_raw == 2 ~ "ocean/lake coastline",
    lw_raw == 3 ~ "shallow inland water",
    lw_raw == 4 ~ "ephemeral water",
    lw_raw == 5 ~ "deep inland water",
    lw_raw == 6 ~ "moderate or continental ocean",
    lw_raw == 7 ~ "deep ocean",
    TRUE ~ NA_character_
  )
  
  data.frame(
    QC_mandatory_quality = mandatory,
    QC_climatology       = climatology,
    QC_land_water        = land_water
  )
}

qc_decoded <- decode_glsp_qc(final_long$GLSP_QC)

final_long <- bind_cols(final_long, qc_decoded)

head(final_long)

# ------------------------------------------------------------
# CLC lookup table (Level 3 full class names)
# ------------------------------------------------------------
clc_lookup <- terra::levels(clc_rast)[[1]] %>%
  dplyr::as_tibble() %>%
  dplyr::rename(clc_code = Value,
                clc_name = LABEL3) %>%
  dplyr::filter(clc_code != 48)   # drop NODATA

final_long <- final_long %>%
  left_join(clc_lookup, by = "clc_code")


write.csv(final_long, "wm_clc.csv", row.names = FALSE)




library(terra)
library(sf)
library(dplyr)
library(ggplot2)

# ------------------------------------------------------------
# One site only
# ------------------------------------------------------------
site_to_run <- "Wytham"

site_buf_sin <- square_buffers %>% filter(pop == site_to_run)
site_buf_wgs <- st_transform(site_buf_sin, wgs_crs)
site_buf_clc <- st_transform(site_buf_sin, clc_crs)

# Keep only this site's extracted VIIRS pixels
site_lsp <- final_long %>%
  filter(pop == site_to_run) %>%
  distinct(x, y, Onset_Greenness_Increase, clc_code, clc_name)

# If you want to create pixel borders explicitly, make a tile for each VIIRS pixel
# using the VIIRS grid spacing from your extracted coordinates.
# This gives visible pixel edges in the LSP plot.
viirs_dx <- min(diff(sort(unique(site_lsp$x))))
viirs_dy <- min(diff(sort(unique(site_lsp$y))))

# ------------------------------------------------------------
# 1) LSP figure: pixel-level borders, no tile borders
# ------------------------------------------------------------
p_lsp <- ggplot(site_lsp) +
  geom_tile(
    aes(x = x, y = y, fill = Onset_Greenness_Increase),
    width = viirs_dx,
    height = viirs_dy,
    colour = "grey35",
    linewidth = 0.15
  ) +
  geom_sf(
    data = site_buf_sin,
    fill = NA, colour = "red", linewidth = 0.8
  ) +
  scale_fill_gradientn(
    colours  = hcl.colors(100, "RdYlGn"),
    name     = "Day of year",
    na.value = "grey80"
  ) +
  coord_sf(crs = sin_crs) +
  labs(
    title = paste0(site_to_run, " — LSP Onset of Greenness (", target_year, ")"),
    subtitle = "Grey lines = VIIRS pixel borders"
  ) +
  theme_minimal(base_size = 10)

ggsave(
  filename = file.path(plot_dir, paste0(site_to_run, "_", target_year, "_LSP.png")),
  plot = p_lsp,
  width = 7, height = 6, dpi = 150
)

# ------------------------------------------------------------
# 2) CLC native-resolution figure
# ------------------------------------------------------------
clc_nat <- crop(clc_rast, vect(site_buf_clc))
clc_df_nat <- as.data.frame(clc_nat, xy = TRUE, na.rm = TRUE)

if (nrow(clc_df_nat) > 0) {
  names(clc_df_nat)[3] <- "clc_code"
  clc_df_nat$clc_code <- as.integer(as.character(clc_df_nat$clc_code))
  clc_df_nat <- left_join(clc_df_nat, clc_lookup, by = "clc_code")
}

p_clc_nat <- ggplot(clc_df_nat) +
  geom_raster(aes(x = x, y = y, fill = clc_name)) +
  geom_sf(
    data = site_buf_clc,
    fill = NA, colour = "red", linewidth = 0.8
  ) +
  scale_fill_viridis_d(name = "CLC class", option = "turbo", na.value = "grey80") +
  labs(title = "CLC — native resolution") +
  theme_minimal(base_size = 10) +
  theme(
    legend.text = element_text(size = 6),
    legend.key.size = unit(0.4, "cm")
  )

ggsave(
  filename = file.path(plot_dir, paste0(site_to_run, "_", target_year, "_CLC_native.png")),
  plot = p_clc_nat,
  width = 7, height = 6, dpi = 150
)

# ------------------------------------------------------------
# 3) CLC resampled to sinusoidal VIIRS grid
# ------------------------------------------------------------
clc_sin <- project(
  crop(clc_rast, vect(site_buf_clc)),
  sin_crs,
  method = "near"
)

clc_df_sin <- as.data.frame(clc_sin, xy = TRUE, na.rm = TRUE)

if (nrow(clc_df_sin) > 0) {
  names(clc_df_sin)[3] <- "clc_code"
  clc_df_sin$clc_code <- as.integer(as.character(clc_df_sin$clc_code))
  clc_df_sin <- left_join(clc_df_sin, clc_lookup, by = "clc_code")
}

p_clc_sin <- ggplot(clc_df_sin) +
  geom_raster(aes(x = x, y = y, fill = clc_name)) +
  geom_sf(
    data = site_buf_sin,
    fill = NA, colour = "red", linewidth = 0.8
  ) +
  scale_fill_viridis_d(name = "CLC class", option = "turbo", na.value = "grey80") +
  coord_sf(crs = sin_crs) +
  labs(title = "CLC — resampled to sinusoidal") +
  theme_minimal(base_size = 10) +
  theme(
    legend.text = element_text(size = 6),
    legend.key.size = unit(0.4, "cm")
  )

ggsave(
  filename = file.path(plot_dir, paste0(site_to_run, "_", target_year, "_CLC_sinu.png")),
  plot = p_clc_sin,
  width = 7, height = 6, dpi = 150
)

# ------------------------------------------------------------
# 4) Satellite basemap
# ------------------------------------------------------------
sat <- fetch_satellite_basemap(site_buf_wgs)

if (!is.null(sat)) {
  sat_df <- as.data.frame(sat, xy = TRUE, na.rm = TRUE)
  band_cols <- setdiff(names(sat_df), c("x", "y"))
  
  if (length(band_cols) >= 3) {
    names(sat_df)[match(band_cols[1:3], names(sat_df))] <- c("r", "g", "b")
    
    p_sat <- ggplot() +
      geom_raster(
        data = sat_df,
        aes(x = x, y = y, fill = rgb(r/255, g/255, b/255))
      ) +
      scale_fill_identity() +
      geom_sf(data = site_buf_wgs, fill = NA, colour = "red", linewidth = 0.8) +
      labs(title = "Satellite image (ESRI World Imagery)") +
      theme_minimal(base_size = 10)
  } else {
    p_sat <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "Satellite image unavailable or not RGB") +
      theme_void()
  }
} else {
  p_sat <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = "Satellite image unavailable") +
    theme_void()
}

ggsave(
  filename = file.path(plot_dir, paste0(site_to_run, "_", target_year, "_satellite.png")),
  plot = p_sat,
  width = 7, height = 6, dpi = 150
)
