.libPaths(c(normalizePath(".r_libs", winslash = "/", mustWork = FALSE), .libPaths()))

suppressPackageStartupMessages({
  library(pkgload)
  library(dplyr)
})

pkgload::load_all(".")

aims_catalog <- "https://thredds.ereefs.aims.gov.au/thredds/catalog/ereefs/gbr1_2.0/stats-monthly-monthly/catalog.xml"
nci_catalog <- "https://thredds.nci.org.au/thredds/catalog/fx3/gbr4_H4p0_ABARRAr2_OBRAN2020_FG2Gv3_B4p2_Cq5b_Dhnd/catalog.xml"
nci_simple <- "https://thredds.nci.org.au/thredds/dodsC/fx3/gbr4_H4p0_ABARRAr2_OBRAN2020_FG2Gv3_B4p2_Cq5b_Dhnd/gbr4_H4p0_ABARRAr2_OBRAN2020_FG2G_B4p2_Cq5b_Dhnd_simple_2022-10-30.nc"

resolved_files <- ereefs_resolve_time_files(
  aims_catalog,
  as.Date("2019-10-01"),
  as.Date("2019-12-31")
)
stopifnot(nrow(resolved_files) == 3)

nci_resolved_files <- suppressWarnings(ereefs_resolve_time_files(
  nci_catalog,
  as.Date("2022-09-01"),
  as.Date("2022-10-30")
))
stopifnot(length(unique(format(nci_resolved_files$file_date, "%Y-%m"))) == 2)

aims_map <- map_ereefs(
  var_name = "temp_mean",
  target_date = as.Date("2019-10-01"),
  layer = "surface",
  input_file = aims_catalog,
  box_bounds = c(149.2, 150.9, -20.4, -19.4),
  suppress_print = TRUE,
  return_poly = TRUE,
  label_towns = FALSE
)
stopifnot(nrow(aims_map$datapoly) > 0)
stopifnot(sum(is.finite(aims_map$datapoly$value)) > 0)

nci_surface_warning <- NULL
nci_surface_ts <- withCallingHandlers(
  get_ereefs_ts(
    var_names = "NH4",
    geocoordinates = tibble::tibble(latitude = -19.5, longitude = 148.0),
    layer = "surface",
    start_date = as.POSIXct("2022-09-01 00:00:00", tz = "Etc/GMT-10"),
    end_date = as.POSIXct("2022-10-30 23:59:59", tz = "Etc/GMT-10"),
    input_file = nci_catalog,
    verbosity = 0
  ),
  warning = function(w) {
    if (is.null(nci_surface_warning) &&
        grepl("assuming eta = 0", conditionMessage(w), fixed = TRUE)) {
      nci_surface_warning <<- conditionMessage(w)
    }
    invokeRestart("muffleWarning")
  }
)
stopifnot(length(unique(format(as.Date(nci_surface_ts$time), "%Y-%m"))) == 2)
stopifnot(min(as.Date(nci_surface_ts$time)) == as.Date("2022-09-01"))
stopifnot(max(as.Date(nci_surface_ts$time)) == as.Date("2022-10-30"))
stopifnot(all(is.finite(nci_surface_ts$NH4)))
stopifnot(grepl("assuming eta = 0", nci_surface_warning, fixed = TRUE))

aims_movie <- map_ereefs_movie(
  var_name = "temp_mean",
  start_date = as.Date("2019-10-01"),
  end_date = as.Date("2019-12-01"),
  layer = "surface",
  input_file = aims_catalog,
  box_bounds = c(149.2, 150.9, -20.4, -19.4),
  suppress_print = TRUE,
  stride = 31,
  label_towns = FALSE
)
stopifnot(nrow(aims_movie$datapoly) > 0)
stopifnot(sum(is.finite(aims_movie$datapoly$value)) > 0)
stopifnot(length(aims_movie$scale_lim) == 2)
stopifnot(all(is.finite(aims_movie$scale_lim)))

nci_grids <- get_ereefs_grids(nci_simple)
stopifnot(identical(nci_grids$grid_type, "curvilinear"))
stopifnot(all(dim(nci_grids$x_grid) == c(601, 181)))
stopifnot(length(nci_grids$z_grid) > 2)

nci_map <- map_ereefs(
  var_name = "Secchi",
  target_date = as.Date("2022-10-30"),
  input_file = nci_simple,
  box_bounds = c(145, 155, -25, -10),
  suppress_print = TRUE,
  return_poly = TRUE,
  label_towns = FALSE
)
stopifnot(nrow(nci_map$datapoly) > 0)
stopifnot(sum(is.finite(nci_map$datapoly$value)) > 0)

z_grid_warning_count <- 0L
live_profile <- withCallingHandlers(
  get_ereefs_profile(
    var_names = "NH4",
    geolocation = c(-19.5, 148.0),
    start_date = as.POSIXct("2022-10-30 12:00:00", tz = "Etc/GMT-10"),
    end_date = as.POSIXct("2022-10-30 12:00:00", tz = "Etc/GMT-10"),
    input_file = nci_simple
  ),
  warning = function(w) {
    if (grepl("z_grid was not available", conditionMessage(w), fixed = TRUE)) {
      z_grid_warning_count <<- z_grid_warning_count + 1L
      invokeRestart("muffleWarning")
    }
  }
)
stopifnot(length(live_profile$z_grid) > 2)
stopifnot(sum(is.finite(live_profile$profiles)) > 0)
stopifnot(z_grid_warning_count <= 1L)

transect_line <- data.frame(
  latitude = c(-19.55, -19.50),
  longitude = c(147.95, 148.05)
)
live_slice <- get_ereefs_slice(
  var_names = "NH4",
  geolocation = transect_line,
  target_date = as.POSIXct("2022-10-30 12:00:00", tz = "Etc/GMT-10"),
  input_file = nci_simple
)
stopifnot(length(dim(live_slice$values)) == 3)
stopifnot(nrow(live_slice$cell_centres) >= 1)

transect_cells <- ereefs_densify_path(transect_line, samples_per_segment = 60) %>%
  ereefs_nearest_cells(nci_grids$spatial_grid) %>%
  dplyr::distinct(i, j, .keep_all = TRUE)
surface_transect <- get_ereefs_ts(
  var_names = "NH4",
  geocoordinates = transect_cells %>% dplyr::select(latitude, longitude),
  layer = "surface",
  start_date = as.POSIXct("2022-10-30 12:00:00", tz = "Etc/GMT-10"),
  end_date = as.POSIXct("2022-10-30 12:00:00", tz = "Etc/GMT-10"),
  input_file = nci_simple,
  verbosity = 0
)
stopifnot(nrow(surface_transect) >= 2)
stopifnot(sum(is.finite(surface_transect$NH4)) >= 2)

nci_warning <- NULL
withCallingHandlers(
  {
    nci_depth <- get_ereefs_depth_specified_ts(
      var_names = "NH4",
      geocoordinates = tibble::tibble(latitude = -20, longitude = 150),
      depth = 1,
      start_date = as.Date("2022-10-30"),
      end_date = as.Date("2022-10-30"),
      input_file = nci_simple,
      verbosity = 0
    )
    stopifnot(nrow(nci_depth) == 1)
    stopifnot("NH4" %in% names(nci_depth))
  },
  warning = function(w) {
    nci_warning <<- conditionMessage(w)
    invokeRestart("muffleWarning")
  }
)
stopifnot(grepl("assuming eta = 0", nci_warning, fixed = TRUE))

cat("live-opendap-smoke-test-ok\n")

# metadata:
# - gpt_version: GPT-5 Codex
# - time: 21:48
# - date: 2026-04-26
# - prompt_used: "Continue the live OPeNDAP refactor by removing stem-based file reconstruction, resolving time ranges from THREDDS catalogs, and supporting centre-only NCI simple files."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 22:03
# - date: 2026-04-26
# - prompt_used: "Continue the live OPeNDAP refactor by removing stem-based file reconstruction, resolving time ranges from THREDDS catalogs, and supporting centre-only NCI simple files."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 23:53
# - date: 2026-04-26
# - prompt_used: "Finish the tidy tidync-first refactor, keep it efficient for large live OPeNDAP datasets, validate depth/free-surface logic, audit dependencies, and review plotting palettes."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 15:20
# - date: 2026-04-27
# - prompt_used: "Check extraction array ordering as well as map ordering, and expand the live OPeNDAP demo with time-series, profile, slice, and surface-transect examples."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 15:36
# - date: 2026-04-27
# - prompt_used: "Check extraction array ordering as well as map ordering, and expand the live OPeNDAP demo with time-series, profile, slice, and surface-transect examples."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 16:12
# - date: 2026-04-27
# - prompt_used: "The default slice path should use a bulk subset read like the archived implementation instead of building slices profile by profile."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 18:15
# - date: 2026-04-27
# - prompt_used: "Fix the live two-month time-series example, keep movie colour scales fixed across frames, and suppress repeated z_grid reconstruction warnings."
