.libPaths(c(normalizePath(".r_libs", winslash = "/", mustWork = FALSE), .libPaths()))

suppressPackageStartupMessages({
  library(pkgload)
  library(dplyr)
})

pkgload::load_all(".")
source(file.path("notebooks", "create_demo_datasets.R"))

demo_paths <- create_demo_datasets(file.path("notebooks", "demo_data"))

regular_grids <- get_ereefs_grids(demo_paths$regular)
stopifnot(identical(regular_grids$grid_type, "regular"))
stopifnot(all(dim(regular_grids$x_grid) == c(5, 4)))

curvilinear_grids <- get_ereefs_grids(demo_paths$curvilinear)
stopifnot(identical(curvilinear_grids$grid_type, "curvilinear"))
stopifnot(all(dim(curvilinear_grids$x_grid) == c(4, 5)))

regular_ts <- get_ereefs_ts(
  var_names = c("temp", "salt"),
  geocoordinates = tibble::tibble(latitude = -19.75, longitude = 146.75),
  layer = "surface",
  start_date = as.POSIXct("2020-01-01 00:00:00", tz = "Etc/GMT-10"),
  end_date = as.POSIXct("2020-01-03 00:00:00", tz = "Etc/GMT-10"),
  input_file = demo_paths$regular,
  verbosity = 0
)
stopifnot(nrow(regular_ts) >= 2)
stopifnot(all(c("temp", "salt", "time") %in% names(regular_ts)))
stopifnot(isTRUE(all.equal(regular_ts$temp, c(21.9, 22.4, 22.9), tolerance = 1e-6)))
stopifnot(isTRUE(all.equal(regular_ts$salt, c(34.77, 34.74, 34.71), tolerance = 1e-6)))

surface_ts <- get_ereefs_ts(
  var_names = "temp",
  geocoordinates = tibble::tibble(latitude = -19.75, longitude = 146.75),
  layer = "surface",
  start_date = as.POSIXct("2020-01-01 00:00:00", tz = "Etc/GMT-10"),
  end_date = as.POSIXct("2020-01-03 00:00:00", tz = "Etc/GMT-10"),
  input_file = demo_paths$regular,
  verbosity = 0
)
bottom_ts <- get_ereefs_ts(
  var_names = "temp",
  geocoordinates = tibble::tibble(latitude = -19.75, longitude = 146.75),
  layer = "bottom",
  start_date = as.POSIXct("2020-01-01 00:00:00", tz = "Etc/GMT-10"),
  end_date = as.POSIXct("2020-01-03 00:00:00", tz = "Etc/GMT-10"),
  input_file = demo_paths$regular,
  verbosity = 0
)
stopifnot(all(surface_ts$temp < bottom_ts$temp))

depth_integrated <- get_ereefs_depth_integrated_ts(
  var_names = "temp",
  geocoordinates = tibble::tibble(latitude = -19.75, longitude = 146.75),
  start_date = as.POSIXct("2020-01-01 00:00:00", tz = "Etc/GMT-10"),
  end_date = as.POSIXct("2020-01-03 00:00:00", tz = "Etc/GMT-10"),
  input_file = demo_paths$regular,
  verbosity = 0
)
stopifnot(all(depth_integrated$temp >= surface_ts$temp))
stopifnot(all(depth_integrated$temp <= bottom_ts$temp))

depth_surface <- get_ereefs_depth_specified_ts(
  var_names = "temp",
  geocoordinates = tibble::tibble(latitude = -19.75, longitude = 146.75),
  depth = 1,
  start_date = as.POSIXct("2020-01-01 00:00:00", tz = "Etc/GMT-10"),
  end_date = as.POSIXct("2020-01-03 00:00:00", tz = "Etc/GMT-10"),
  input_file = demo_paths$regular,
  verbosity = 0
)
stopifnot(isTRUE(all.equal(depth_surface$temp, surface_ts$temp, tolerance = 1e-6)))

noeta_warn <- NULL
withCallingHandlers(
  {
    depth_noeta <- get_ereefs_depth_specified_ts(
      var_names = "temp",
      geocoordinates = tibble::tibble(latitude = -19.75, longitude = 146.75),
      depth = 1,
      start_date = as.POSIXct("2020-01-01 00:00:00", tz = "Etc/GMT-10"),
      end_date = as.POSIXct("2020-01-03 00:00:00", tz = "Etc/GMT-10"),
      input_file = demo_paths$regular_noeta,
      verbosity = 0
    )
    stopifnot(nrow(depth_noeta) == 3)
  },
  warning = function(w) {
    noeta_warn <<- conditionMessage(w)
    invokeRestart("muffleWarning")
  }
)
stopifnot(grepl("assuming eta = 0", noeta_warn, fixed = TRUE))

regular_map <- map_ereefs(
  var_name = "temp",
  target_date = as.Date("2020-01-02"),
  layer = 1,
  input_file = demo_paths$regular,
  box_bounds = c(146.2, 147.3, -20.4, -19.2),
  suppress_print = TRUE,
  return_poly = TRUE,
  label_towns = FALSE
)
stopifnot(nrow(regular_map$datapoly) > 0)

movie_dir <- file.path(tempdir(), "ereefs_movie_smoke")
movie_out <- map_ereefs_movie(
  var_name = "temp",
  start_date = as.Date("2020-01-01"),
  end_date = as.Date("2020-01-03"),
  layer = "surface",
  input_file = demo_paths$regular,
  output_dir = movie_dir,
  save_frames = TRUE,
  animation_format = "gif",
  stride = "daily",
  suppress_print = TRUE,
  label_towns = FALSE
)
stopifnot(length(movie_out$frame_files) == 3)
stopifnot(all(file.exists(movie_out$frame_files)))
stopifnot(file.exists(movie_out$animation_file))
stopifnot(length(movie_out$scale_lim) == 2)
stopifnot(all(is.finite(movie_out$scale_lim)))

profile <- get_ereefs_profile(
  var_names = "temp",
  geolocation = c(2L, 2L),
  start_date = as.POSIXct("2020-01-01 00:00:00", tz = "Etc/GMT-10"),
  end_date = as.POSIXct("2020-01-02 00:00:00", tz = "Etc/GMT-10"),
  input_file = demo_paths$curvilinear
)
stopifnot(length(profile$dates) >= 1)
stopifnot(length(profile$z_grid) == 3)
stopifnot(isTRUE(all.equal(
  profile$profiles,
  matrix(c(23.7, 22.8, 24.0, 23.1), nrow = 2, ncol = 2),
  tolerance = 1e-6
)))

slice <- get_ereefs_slice(
  var_names = "temp",
  geolocation = data.frame(latitude = c(-20.05, -18.7), longitude = c(147.05, 148.35)),
  target_date = as.POSIXct("2020-01-01 00:00:00", tz = "Etc/GMT-10"),
  input_file = demo_paths$curvilinear
)
stopifnot(length(dim(slice$values)) == 3)
stopifnot(nrow(slice$cell_centres) >= 1)
slice_plot <- plot_ereefs_slice(slice, var_name = "temp", scale_col = "viridis")
stopifnot(inherits(slice_plot, "ggplot"))

cat("smoke-test-ok\n")

# metadata:
# - gpt_version: GPT-5 Codex
# - time: 19:17
# - date: 2026-04-26
# - prompt_used: "Install dependencies, verify the refactored toolkit, improve efficiency for large THREDDS-served files, and build a working Jupyter demo notebook."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 19:49
# - date: 2026-04-26
# - prompt_used: "Install dependencies, verify the refactored toolkit, improve efficiency for large THREDDS-served files, and build a working Jupyter demo notebook."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 23:46
# - date: 2026-04-26
# - prompt_used: "Finish the tidy tidync-first refactor, keep it efficient for large live OPeNDAP datasets, validate depth/free-surface logic, audit dependencies, and review plotting palettes."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 15:20
# - date: 2026-04-27
# - prompt_used: "Check extraction array ordering as well as map ordering, and expand the live OPeNDAP demo with time-series, profile, slice, and surface-transect examples."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 16:58
# - date: 2026-04-27
# - prompt_used: "Make notebook plots save as well as display, and have map_ereefs_movie save frames then build an animation file."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 18:15
# - date: 2026-04-27
# - prompt_used: "Fix the live two-month time-series example, keep movie colour scales fixed across frames, and suppress repeated z_grid reconstruction warnings."
