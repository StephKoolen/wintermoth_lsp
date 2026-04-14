library(terra)
library(sf)
library(dplyr)
library(exactextractr)
library(stringr)

# ── CRS & paths ───────────────────────────────────────────────────────────────
sin_crs  <- "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"
wgs_crs  <- "EPSG:4326"

viirs_dir  <- "../../wintermoth_data/VIIRS"
clc_path   <- "../../wintermoth_data/CLC/data/U2018_CLC2018_V2020_20u1.tif"
meta_path  <- "../../wintermoth_data/wgs_europe_metadata.csv"
out_dir    <- "../output"
dir.create(out_dir, showWarnings = FALSE)

# ── Helper ────────────────────────────────────────────────────────────────────
make_square <- function(x, y, hx, hy) {
  st_polygon(list(matrix(c(
    x - hx, y - hy,
    x + hx, y - hy,
    x + hx, y + hy,
    x - hx, y + hy,
    x - hx, y - hy
  ), ncol = 2, byrow = TRUE)))
}

# ── Site buffers ──────────────────────────────────────────────────────────────
wm_location <- read.csv(meta_path) %>%
  distinct(pop, latitude, longitude, .keep_all = TRUE)

location_sf  <- st_as_sf(wm_location, coords = c("longitude", "latitude"), crs = 4326)
location_sin <- st_transform(location_sf, crs = sin_crs)
coords       <- st_coordinates(location_sin)
half         <- 2500

squares <- mapply(
  make_square,
  coords[, 1], coords[, 2],
  MoreArgs = list(hx = half, hy = half),
  SIMPLIFY = FALSE
)

square_buffers <- st_sf(
  st_drop_geometry(location_sin),
  geometry = st_sfc(squares, crs = sin_crs)
)

buf_vect <- vect(square_buffers)

# ── CLC — extract once per site ───────────────────────────────────────────────
message("Loading CLC raster...")
clc_rast <- rast(clc_path)
clc_crs  <- crs(clc_rast)
if (nlyr(clc_rast) > 1) clc_rast <- clc_rast[[1]]

message("Extracting CLC per site...")
buf_vect_clc <- project(buf_vect, clc_crs)

clc_extracted <- exact_extract(
  clc_rast,
  st_as_sf(buf_vect_clc),
  fun       = "mode",
  append_cols = FALSE
)

clc_per_site <- square_buffers %>%
  st_drop_geometry() %>%
  select(pop) %>%
  mutate(clc_code = as.integer(clc_extracted))

message("CLC extraction done.")

# ── VIIRS — extract per year ───────────────────────────────────────────────────
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

extract_vnp <- function(filepath, year) {
  tile_size <- 1111950.5196
  tile_str  <- regmatches(basename(filepath), regexpr("h\\d{2}v\\d{2}", basename(filepath)))
  tile_h    <- as.integer(substr(tile_str, 2, 3))
  tile_v    <- as.integer(substr(tile_str, 5, 6))
  
  xmin <- (tile_h - 18) * tile_size
  xmax <- xmin + tile_size
  ymax <- (9 - tile_v)  * tile_size
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

# ── Loop over all years found in the VIIRS directory ─────────────────────────
all_files <- list.files(viirs_dir, pattern = "VNP22Q2.*\\.h5$", full.names = TRUE)

if (length(all_files) == 0) stop("No VIIRS .h5 files found in: ", viirs_dir)

all_years <- as.integer(regmatches(basename(all_files), regexpr("\\d{4}", basename(all_files))))
message("Years found: ", paste(sort(unique(all_years)), collapse = ", "))

results_by_year <- list()

for (yr in sort(unique(all_years))) {
  message("\nProcessing year: ", yr)
  year_files   <- all_files[all_years == yr]
  tile_results <- lapply(year_files, extract_vnp, year = yr)
  
  year_extracts <- lapply(tile_results, function(tile) {
    r       <- tile$raster
    tile_id <- tile$tile_id
    
    tile_poly     <- as.polygons(ext(r), crs = sin_crs)
    buf_intersect <- relate(buf_vect, tile_poly, relation = "intersects")
    buf_subset    <- buf_vect[apply(buf_intersect, 1, any), ]
    if (nrow(buf_subset) == 0) return(NULL)
    
    r_crop    <- crop(r, buf_subset)
    extracted <- extract(r_crop, buf_subset, fun = "mean", na.rm = TRUE, ID = TRUE)
    if (is.null(extracted) || nrow(extracted) == 0) return(NULL)
    
    buf_df <- as.data.frame(buf_subset) %>%
      mutate(.row_id = seq_len(n()))
    
    extracted %>%
      left_join(buf_df %>% select(.row_id, pop), by = c("ID" = ".row_id")) %>%
      mutate(tile_id = tile_id, year = yr) %>%
      select(-ID)
  })
  
  year_long <- bind_rows(Filter(Negate(is.null), year_extracts))
  
  if (nrow(year_long) == 0) {
    message("  No data extracted for ", yr, " — skipping")
    next
  }
  
  results_by_year[[as.character(yr)]] <- year_long
  message("  Extracted ", nrow(year_long), " site-year rows for ", yr)
}

# ── Combine, decode QC, join CLC ──────────────────────────────────────────────
message("\nCombining all years...")

final_long <- bind_rows(results_by_year) %>%
  filter(!is.na(pop)) %>%
  left_join(clc_per_site, by = "pop") %>%
  bind_cols(decode_glsp_qc(.$GLSP_QC)) %>%
  select(
    pop, year, tile_id,
    Onset_Greenness_Increase, GLSP_QC, PGQ_Onset_Greenness_Increase,
    QC_mandatory_quality, QC_climatology, QC_land_water,
    clc_code
  ) %>%
  arrange(pop, year)

# ── Save ──────────────────────────────────────────────────────────────────────
out_path <- file.path(out_dir, "wm_lsp_clc_all_years.csv")
write.csv(final_long, out_path, row.names = FALSE)
message("Saved to: ", out_path)
message("Done. ", nrow(final_long), " rows, ", n_distinct(final_long$pop), " sites, ",
        n_distinct(final_long$year), " years.")

