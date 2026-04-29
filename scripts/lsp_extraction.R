library(terra)
library(sf)
library(dplyr)
library(tidyr)
library(exactextractr)
library(stringr)

# ── CRS & paths ───────────────────────────────────────────────────────────────
sin_crs  <- "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"
wgs_crs  <- "EPSG:4326"
out_dir  <- "../output"
dir.create(out_dir, showWarnings = FALSE)

# ── Helper functions ───────────────────────────────────────────────────────────
make_square <- function(x, y, hx, hy) {
  st_polygon(list(matrix(c(
    x - hx, y - hy,
    x + hx, y - hy,
    x + hx, y + hy,
    x - hx, y + hy,
    x - hx, y - hy
  ), ncol = 2, byrow = TRUE)))
}

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
  ymax <- (9  - tile_v) * tile_size
  ymin <- ymax - tile_size
  tile_ext <- ext(xmin, xmax, ymin, ymax)
  
  onset_path <- paste0('HDF5:"', filepath, '"://HDFEOS/GRIDS/Cycle_1/Data_Fields/Onset_Greenness_Increase_1')
  qc_path    <- paste0('HDF5:"', filepath, '"://HDFEOS/GRIDS/Cycle_1/Data_Fields/GLSP_QC_1')
  pgq_path   <- paste0('HDF5:"', filepath, '"://HDFEOS/GRIDS/Cycle_1/Data_Fields/PGQ_Onset_Greenness_Increase_1')
  
  onset <- flip(rast(onset_path), direction = "vertical")
  qc    <- flip(rast(qc_path),    direction = "vertical")
  pgq   <- flip(rast(pgq_path),   direction = "vertical")
  
  for (r in list(onset, qc, pgq)) { ext(r) <- tile_ext; crs(r) <- sin_crs }
  
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

# ── Site buffers ───────────────────────────────────────────────────────────────
square_buffers <- st_read("../../wintermoth_data/10kmbuffers.shp") %>%
  st_transform(crs = sin_crs) %>% 
  unite(pop_id, pop, ID, sep = "_")

buf_vect <- vect(square_buffers)

# ── CLC (extracted once; labels joined later) ──────────────────────────────────
clc_rast <- rast("../../wintermoth_data/CLC/data/U2018_CLC2018_V2020_20u1.tif")
clc_crs  <- crs(clc_rast)
if (nlyr(clc_rast) > 1) clc_rast <- clc_rast[[1]]

clc_lookup <- terra::levels(clc_rast)[[1]] %>%
  as_tibble() %>%
  rename(clc_code = Value, clc_name = LABEL3) %>%
  mutate(clc_code = as.integer(clc_code)) %>%
  filter(clc_code != 48)

# ── Discover all years in the VIIRS folder ─────────────────────────────────────
all_files <- list.files("../../wintermoth_data/VIIRS",
                        pattern = "VNP22Q2.*\\.h5$", full.names = TRUE)
all_years <- sort(unique(
  as.integer(regmatches(basename(all_files),
                        regexpr("\\d{4}", basename(all_files))))
))
message("Years found: ", paste(all_years, collapse = ", "))

# ── Main extraction loop (one year at a time) ──────────────────────────────────
all_years_long <- vector("list", length(all_years))

for (yi in seq_along(all_years)) {
  yr <- all_years[yi]
  message("\n── Year ", yr, " ──────────────────────────────────")
  
  year_files   <- all_files[grepl(as.character(yr), basename(all_files))]
  tile_results <- lapply(year_files, extract_vnp, year = yr)
  
  year_extracts <- lapply(tile_results, function(tile) {
    r       <- tile$raster
    tile_id <- tile$tile_id
    
    tile_poly     <- as.polygons(ext(r), crs = sin_crs)
    buf_intersect <- relate(buf_vect, tile_poly, relation = "intersects")
    buf_subset    <- buf_vect[buf_intersect, ]
    if (nrow(buf_subset) == 0) return(NULL)
    
    r_crop <- crop(r, buf_subset)
    
    # CLC: project tile footprint → native CLC crs, crop, reproject to sinusoidal
    tile_poly_clc <- project(tile_poly, clc_crs)
    clc_crop_sin  <- project(crop(clc_rast, tile_poly_clc), sin_crs, method = "near")
    
    extracted <- terra::extract(r_crop, buf_subset,
                         fun = NULL, ID = TRUE, cells = TRUE, xy = TRUE)
    if (is.null(extracted) || nrow(extracted) == 0) return(NULL)
    
    buf_df <- as.data.frame(buf_subset) %>% mutate(.row_id = seq_len(n()))
    
    extracted <- extracted %>%
      left_join(buf_df %>% select(.row_id, pop_id), by = c("ID" = ".row_id")) %>%
      rename(pixel_id = cell)
    
    extracted$tile_id <- tile_id
    extracted$year    <- yr
    
    # Build pixel polygons → extract modal CLC code
    pixel_tbl <- extracted %>% distinct(pixel_id, x, y)
    hx <- res(r_crop)[1] / 2
    hy <- res(r_crop)[2] / 2
    
    pixel_sf <- st_sf(
      pixel_tbl,
      geometry = st_sfc(
        lapply(seq_len(nrow(pixel_tbl)), function(i)
          make_square(pixel_tbl$x[i], pixel_tbl$y[i], hx, hy)),
        crs = sin_crs
      )
    )
    pixel_sf$clc_code <- as.integer(exact_extract(clc_crop_sin, pixel_sf, "mode"))
    
    extracted %>%
      left_join(st_drop_geometry(pixel_sf) %>% select(pixel_id, clc_code),
                by = "pixel_id") %>%
      select(-ID)
  })
  
  year_long <- bind_rows(Filter(Negate(is.null), year_extracts)) %>%
    select(pop_id, tile_id, year, pixel_id, x, y,
           Onset_Greenness_Increase, GLSP_QC, PGQ_Onset_Greenness_Increase,
           clc_code) %>%
    arrange(pop_id, tile_id, pixel_id)
  
  # Decode QC flags & attach CLC label
  year_long <- bind_cols(year_long, decode_glsp_qc(year_long$GLSP_QC)) %>%
    left_join(clc_lookup, by = "clc_code")
  
  all_years_long[[yi]] <- year_long
  message("  Rows extracted: ", nrow(year_long))
}

# ── Combine & write long format ────────────────────────────────────────────────
final_long <- bind_rows(all_years_long) %>%
  arrange(pop_id, year, tile_id, pixel_id)

write.csv(final_long,
          file.path(out_dir, "wm_lsp_clc_long.csv"),
          row.names = FALSE)
message("Long format written: wm_lsp_clc_long.csv  (", nrow(final_long), " rows)")

# ── Wide format ────────────────────────────────────────────────────────────────
# One row per pixel (pop_id × pixel_id), one column per year for each LSP variable.
# CLC code/name are stable across years so we take the modal value.

stable_cols <- final_long %>%
  group_by(pop_id, pixel_id) %>%
  summarise(
    x        = first(x),
    y        = first(y),
    tile_id  = first(tile_id),
    clc_code = as.integer(names(sort(table(clc_code), decreasing = TRUE))[1]),
    clc_name = as.character(names(sort(table(clc_name), decreasing = TRUE))[1]),
    .groups  = "drop"
  )

lsp_wide <- final_long %>%
  select(pop_id, pixel_id, year,
         Onset_Greenness_Increase, GLSP_QC, PGQ_Onset_Greenness_Increase,
         QC_mandatory_quality, QC_climatology, QC_land_water) %>%
  tidyr::pivot_wider(
    names_from  = year,
    values_from = c(Onset_Greenness_Increase, GLSP_QC, PGQ_Onset_Greenness_Increase,
                    QC_mandatory_quality, QC_climatology, QC_land_water),
    names_glue  = "{.value}_{year}"
  )

final_wide <- stable_cols %>%
  left_join(lsp_wide, by = c("pop_id", "pixel_id")) %>%
  arrange(pop_id, pixel_id)

write.csv(final_wide,
          file.path(out_dir, "wm_lsp_clc_wide.csv"),
          row.names = FALSE)
message("Wide format written: wm_lsp_clc_wide.csv  (", nrow(final_wide), " rows, ",
        ncol(final_wide), " columns)")

message("\nDone.")
