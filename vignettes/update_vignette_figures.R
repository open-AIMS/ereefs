.libPaths(c(normalizePath(".r_libs", winslash = "/", mustWork = FALSE), .libPaths()))

suppressPackageStartupMessages({
  library(pkgload)
  library(ggplot2)
})

pkgload::load_all(".")
source(file.path("notebooks", "create_demo_datasets.R"))

demo_paths <- create_demo_datasets(file.path("notebooks", "demo_data"))
vignette_dir <- "vignettes"
nci_simple <- "https://thredds.nci.org.au/thredds/dodsC/fx3/gbr4_H4p0_ABARRAr2_OBRAN2020_FG2Gv3_B4p2_Cq5b_Dhnd/gbr4_H4p0_ABARRAr2_OBRAN2020_FG2G_B4p2_Cq5b_Dhnd_simple_2022-10-30.nc"
nci_hyd_simple <- "https://thredds.nci.org.au/thredds/dodsC/fx3/gbr4_H4p0_ABARRAr2_OBRAN2020_FG2Gv3_Dhnd/gbr4_simple_2022-10-30.nc"
aims_catalog <- "https://thredds.ereefs.aims.gov.au/thredds/catalog/ereefs/gbr1_2.0/stats-monthly-monthly/catalog.xml"

save_vignette_plot <- function(plot_obj, filename, width = 9, height = 6) {
  ggplot2::ggsave(
    filename = file.path(vignette_dir, filename),
    plot = plot_obj,
    width = width,
    height = height,
    dpi = 120
  )
}

exmp1 <- map_ereefs(
  var_name = "Secchi",
  target_date = as.Date("2022-10-30"),
  input_file = nci_simple,
  suppress_print = TRUE
)
save_vignette_plot(exmp1, "vignette-fig-exmp1-1.png")

exmp1b <- map_ereefs(
  var_name = "Secchi",
  target_date = as.Date("2022-10-30"),
  input_file = nci_simple,
  plot_style = "smooth",
  smooth_pixels = 900,
  suppress_print = TRUE
)
save_vignette_plot(exmp1b, "vignette-fig-exmp1b-1.png")

exmp2 <- map_ereefs(
  var_name = "temp_mean",
  target_date = as.Date("2019-10-01"),
  input_file = aims_catalog,
  box_bounds = c(149.2, 150.9, -20.4, -19.4),
  suppress_print = TRUE
)
save_vignette_plot(exmp2, "vignette-fig-exmp2-1.png")

exmp3 <- map_ereefs(
  var_name = "temp_mean",
  target_date = as.Date("2019-10-01"),
  input_file = aims_catalog,
  box_bounds = c(149.2, 150.9, -20.4, -19.4),
  plot_style = "smooth",
  smooth_pixels = 700,
  suppress_print = TRUE
)
save_vignette_plot(exmp3, "vignette-fig-exmp3-1.png")

exmp4 <- map_ereefs(
  var_name = "temp_mean",
  target_date = as.Date("2019-10-01"),
  input_file = aims_catalog,
  box_bounds = c(149.2, 150.9, -20.4, -19.4),
  plot_style = "smooth",
  smooth_pixels = 700,
  scale_lim = c(24, 27),
  suppress_print = TRUE
)
save_vignette_plot(exmp4, "vignette-fig-exmp4-1.png")

exmp5 <- map_ereefs(
  var_name = "NH4",
  target_date = as.Date("2022-10-30"),
  layer = -5,
  input_file = nci_simple,
  Land_map = TRUE,
  box_bounds = c(145, 150, -22, -18),
  scale_lim = c(0, 10),
  suppress_print = TRUE
)
save_vignette_plot(exmp5, "vignette-fig-exmp5-1.png")

temp_slice <- get_ereefs_slice(
  var_names = "NH4",
  target_date = as.POSIXct("2022-10-30 12:00:00", tz = "Etc/GMT-10"),
  geolocation = data.frame(latitude = c(-19.26639219, -19.26639219), longitude = c(146.805701, 147.38)),
  input_file = nci_simple
)
exmp8 <- plot_ereefs_slice(temp_slice, var_name = "NH4", var_units = "mg N m-3")
save_vignette_plot(exmp8, "vignette-fig-exmp8-1.png")

arc_transect <- data.frame(
  latitude = c(-19.26639219, -19.18, -19.02, -18.84, -18.70),
  longitude = c(146.805701, 146.95, 147.12, 147.33, 147.56)
)
arc_slice <- get_ereefs_slice(
  var_names = "temp",
  target_date = as.POSIXct("2022-10-30 12:00:00", tz = "Etc/GMT-10"),
  geolocation = arc_transect,
  input_file = nci_hyd_simple
)
arc_map <- map_ereefs(
  var_name = "temp",
  target_date = as.POSIXct("2022-10-30 12:00:00", tz = "Etc/GMT-10"),
  layer = "surface",
  input_file = nci_hyd_simple,
  box_bounds = c(146.6, 147.8, -19.5, -18.5),
  suppress_print = TRUE
) +
  geom_path(
    data = arc_transect,
    aes(x = longitude, y = latitude),
    inherit.aes = FALSE,
    colour = "black",
    linewidth = 1.6,
    lineend = "round"
  )
save_vignette_plot(arc_map, "vignette-fig-exmp8arc-map-1.png")
arc_slice_plot <- plot_ereefs_slice(
  arc_slice,
  var_name = "temp",
  scale_col = "viridis",
  var_units = "degrees C"
)
save_vignette_plot(arc_slice_plot, "vignette-fig-exmp8arc-slice-1.png")

profile_data <- get_ereefs_profile(
  var_names = "NH4",
  start_date = as.POSIXct("2022-10-30 12:00:00", tz = "Etc/GMT-10"),
  end_date = as.POSIXct("2022-10-30 12:00:00", tz = "Etc/GMT-10"),
  geolocation = c(-19.5, 148.0),
  input_file = nci_simple
)
exmp10 <- plot_ereefs_profile(
  profile_data,
  var_name = "NH4",
  target_date = as.POSIXct("2022-10-30 12:00:00", tz = "Etc/GMT-10")
)
save_vignette_plot(exmp10, "vignette-fig-exmp10-1.png")

tsdata <- get_ereefs_ts(
  var_names = c("temp", "salt"),
  geocoordinates = data.frame(latitude = -19.75, longitude = 146.75),
  layer = "surface",
  start_date = as.POSIXct("2020-01-01 00:00:00", tz = "Etc/GMT-10"),
  end_date = as.POSIXct("2020-01-03 00:00:00", tz = "Etc/GMT-10"),
  input_file = demo_paths$regular,
  verbosity = 0
)
exmp11 <- ggplot(tsdata, aes(x = time)) +
  geom_line(aes(y = temp, colour = "temp"), linewidth = 0.8) +
  geom_point(aes(y = temp, colour = "temp"), size = 1.7) +
  geom_line(aes(y = salt, colour = "salt"), linewidth = 0.8) +
  geom_point(aes(y = salt, colour = "salt"), size = 1.7) +
  labs(x = NULL, y = "value", colour = "variable") +
  theme_minimal(base_size = 12)
save_vignette_plot(exmp11, "vignette-fig-exmp11-1.png", width = 9, height = 5)

deeperdata <- get_ereefs_ts(
  var_names = "temp",
  geocoordinates = data.frame(latitude = -19.75, longitude = 146.75),
  layer = -5,
  start_date = as.POSIXct("2020-01-01 00:00:00", tz = "Etc/GMT-10"),
  end_date = as.POSIXct("2020-01-03 00:00:00", tz = "Etc/GMT-10"),
  input_file = demo_paths$regular,
  verbosity = 0
)
intdata <- get_ereefs_depth_integrated_ts(
  var_names = "temp",
  geocoordinates = data.frame(latitude = -19.75, longitude = 146.75),
  start_date = as.POSIXct("2020-01-01 00:00:00", tz = "Etc/GMT-10"),
  end_date = as.POSIXct("2020-01-03 00:00:00", tz = "Etc/GMT-10"),
  input_file = demo_paths$regular,
  verbosity = 0
)
tidalintdata <- get_ereefs_depth_specified_ts(
  var_names = "temp",
  geocoordinates = data.frame(latitude = -19.75, longitude = 146.75),
  depth = 2.0,
  start_date = as.POSIXct("2020-01-01 00:00:00", tz = "Etc/GMT-10"),
  end_date = as.POSIXct("2020-01-03 00:00:00", tz = "Etc/GMT-10"),
  input_file = demo_paths$regular,
  verbosity = 0
)
depth_plot_data <- rbind(
  data.frame(time = deeperdata$time, value = deeperdata$temp, series = "5 m below MSL"),
  data.frame(time = intdata$time, value = intdata$temp, series = "depth average"),
  data.frame(time = tidalintdata$time, value = tidalintdata$temp, series = "2 m below free surface")
)
exmp12 <- ggplot(depth_plot_data, aes(x = time, y = value, colour = series)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.7) +
  labs(x = NULL, y = "temperature") +
  theme_minimal(base_size = 12)
save_vignette_plot(exmp12, "vignette-fig-exmp12-1.png", width = 9, height = 5)

cat("vignette-figures-updated\n")

# metadata:
# - gpt_version: GPT-5 Codex
# - time: 08:41
# - date: 2026-04-28
# - prompt_used: "Update vignette PNGs to match the current code and make viridis the default colour palette."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 09:03
# - date: 2026-04-28
# - prompt_used: "Fix the vignette profile and time-series examples, avoid clipped first-map labels, and add explanatory markdown blocks throughout the notebooks."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 09:22
# - date: 2026-04-28
# - prompt_used: "Fix the remaining howto vignette map clipping, add the missing smooth-map figure, and switch the profile and slice examples to live OPeNDAP data."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 17:04
# - date: 2026-04-28
# - prompt_used: "Document multi-point curved transects in the vignette and add a live hydrodynamic arc-shaped slice example with both a map and slice figure."
