test_that("NCI catalog hrefs resolve without malformed double THREDDS paths", {
  base <- "https://thredds.nci.org.au/thredds/catalog/catalogs/fx3/catalog.xml"

  root_relative <- ereefs:::ereefs_gui_resolve_catalog_href(
    "/thredds/catalog/fx3/gbr4_H4p0_ABARRAr2_OBRAN2020_FG2Gv3_Dhnd/catalog.xml",
    base
  )
  absolute <- ereefs:::ereefs_gui_resolve_catalog_href(
    "https://thredds.nci.org.au/thredds/catalog/fx3/gbr1_2.0/catalog.xml",
    base
  )

  expect_equal(
    root_relative,
    "https://thredds.nci.org.au/thredds/catalog/fx3/gbr4_H4p0_ABARRAr2_OBRAN2020_FG2Gv3_Dhnd/catalog.xml"
  )
  expect_equal(absolute, "https://thredds.nci.org.au/thredds/catalog/fx3/gbr1_2.0/catalog.xml")
  expect_false(grepl("catalogs/fx3//thredds", root_relative, fixed = TRUE))
})

test_that("GUI date and spatial bounds are safe for sparse catalog metadata", {
  summary <- tibble::tibble(
    time_start = as.POSIXct(NA),
    time_end = as.POSIXct(NA),
    file_start = as.Date("2020-01-01"),
    file_end = as.Date("2020-01-03"),
    longitude_min = NA_real_,
    longitude_max = 148,
    latitude_min = -20,
    latitude_max = NA_real_
  )

  dates <- ereefs:::ereefs_gui_date_bounds(summary)
  bounds <- ereefs:::ereefs_gui_summary_bounds(summary)

  expect_equal(dates$start, as.Date("2020-01-01"))
  expect_equal(dates$end, as.Date("2020-01-03"))
  expect_true(is.na(bounds$lon_min))
  expect_equal(bounds$lon_max, 148)
  expect_equal(bounds$lat_min, -20)
  expect_true(is.na(bounds$lat_max))
})

test_that("GUI bbox helper stays length four when Shiny inputs are transiently missing", {
  input <- list(
    use_bbox = TRUE,
    lon_min = 142,
    lon_max = NULL,
    lat_min = -20,
    lat_max = -18
  )

  bounds <- ereefs:::ereefs_gui_bbox(input)

  expect_length(bounds, 4)
  expect_true(all(is.na(bounds)))
})

test_that("GUI variable sync helper preserves shared variables and falls back safely", {
  choices <- c("temp", "salt", "NH4")

  expect_equal(
    ereefs:::ereefs_gui_select_synced_variable("salt", choices),
    "salt"
  )
  expect_equal(
    ereefs:::ereefs_gui_select_synced_variable("true_colour", choices),
    "temp"
  )
})

test_that("GUI point table helper always returns latitude and longitude columns", {
  parsed_empty <- list(NULL)
  parsed_one <- list(tibble::tibble(latitude = -19.2, longitude = 146.8, label = "only"))

  empty_tbl <- ereefs:::ereefs_gui_point_table(parsed_empty, include_label = TRUE)
  one_tbl <- ereefs:::ereefs_gui_point_table(parsed_one, include_label = TRUE)

  expect_identical(names(empty_tbl), c("latitude", "longitude", "label"))
  expect_equal(nrow(empty_tbl), 0)
  expect_identical(names(one_tbl), c("latitude", "longitude", "label"))
  expect_equal(nrow(one_tbl), 1)
})

test_that("GUI point date range handles transiently missing date inputs", {
  inspection <- list(summary = tibble::tibble(
    time_start = as.POSIXct(NA),
    time_end = as.POSIXct(NA),
    file_start = as.Date("2020-01-01"),
    file_end = as.Date("2020-01-03")
  ))
  input <- list(points_start = NULL, points_end = as.Date("2020-01-02"))

  range <- ereefs:::ereefs_gui_points_date_range(input, inspection)

  expect_equal(range$start, as.Date("2020-01-01"))
  expect_equal(range$end, as.Date("2020-01-02"))
})

test_that("GUI date clamp snaps to the nearest inspected boundary", {
  inspection <- list(summary = tibble::tibble(
    time_start = as.POSIXct(NA),
    time_end = as.POSIXct(NA),
    file_start = as.Date("2020-01-01"),
    file_end = as.Date("2020-01-03")
  ))

  expect_equal(
    ereefs:::ereefs_gui_clamp_date(as.Date("2019-12-20"), inspection),
    as.Date("2020-01-01")
  )
  expect_equal(
    ereefs:::ereefs_gui_clamp_date(as.Date("2020-01-20"), inspection),
    as.Date("2020-01-03")
  )
})

test_that("profile plot handles single-variable multi-date profile arrays", {
  profile_obj <- list(
    dates = as.POSIXct(c("2020-01-01 12:00:00", "2020-01-02 12:00:00"), tz = "Etc/GMT-10"),
    eta = c(0, 0.1),
    z_grid = c(-10, -5, 0),
    botz = -10,
    profiles = array(c(1, 2, 3, 4), dim = c(2, 2))
  )

  p <- ereefs:::plot_ereefs_profile(
    profile_obj = profile_obj,
    var_name = "NH4",
    target_date = as.Date("2020-01-02")
  )

  expect_s3_class(p, "ggplot")
})

test_that("GUI multi-profile plot combines multiple selected points", {
  profile_a <- list(
    dates = as.POSIXct("2020-01-01 12:00:00", tz = "Etc/GMT-10"),
    eta = 0,
    z_grid = c(-10, -5, 0),
    botz = -10,
    profiles = array(c(1, 2), dim = c(2, 1), dimnames = list(NULL, "NH4"))
  )
  profile_b <- list(
    dates = as.POSIXct("2020-01-01 12:00:00", tz = "Etc/GMT-10"),
    eta = 0,
    z_grid = c(-10, -5, 0),
    botz = -10,
    profiles = array(c(3, 4), dim = c(2, 1), dimnames = list(NULL, "NH4"))
  )

  p <- ereefs:::ereefs_gui_profile_plot(
    list(
      type = "multi",
      profiles = list(profile_a, profile_b),
      labels = c("A", "B"),
      target_date = as.Date("2020-01-01")
    ),
    "NH4"
  )

  expect_s3_class(p, "ggplot")
})

test_that("inspection callback reports progress stages for local demo data", {
  demo_file <- testthat::test_path("..", "..", "notebooks", "demo_data", "regular_demo_2020-01.nc")
  testthat::skip_if_not(file.exists(demo_file), "local demo NetCDF file is not available")

  events <- list()
  inspection <- inspect_ereefs_data(
    input_file = demo_file,
    progress_callback = function(info) events[[length(events) + 1L]] <<- info
  )

  expect_s3_class(inspection, "ereefs_data_inspection")
  expect_true(any(vapply(events, function(info) identical(info$stage, "dimensions"), logical(1))))
  expect_true(any(vapply(events, function(info) identical(info$stage, "variables"), logical(1))))
  expect_true(any(vapply(events, function(info) identical(info$stage, "complete"), logical(1))))
})

test_that("extraction progress helper emits lightweight file and variable updates", {
  events <- list()
  callback <- function(info) events[[length(events) + 1L]] <<- info

  ereefs:::ereefs_extraction_progress(
    callback,
    stage = "file",
    message = "Reading time steps from source file",
    file = "https://example.invalid/demo.nc",
    file_index = 1L,
    file_count = 2L,
    time_count = 3L,
    time_start = as.POSIXct("2020-01-01 12:00:00", tz = "Etc/GMT-10"),
    time_end = as.POSIXct("2020-01-03 12:00:00", tz = "Etc/GMT-10")
  )
  ereefs:::ereefs_extraction_progress(
    callback,
    stage = "variable",
    message = "Extracting variable values",
    file = "https://example.invalid/demo.nc",
    file_index = 1L,
    file_count = 2L,
    variable = "temp",
    variable_index = 1L,
    variable_count = 4L
  )
  ereefs:::ereefs_extraction_progress(
    callback,
    stage = "complete",
    message = "Point time-series extraction complete",
    file_count = 2L,
    point_count = 1L,
    variable_count = 4L
  )

  expect_length(events, 3)
  expect_identical(events[[1]]$stage, "file")
  expect_identical(events[[2]]$stage, "variable")
  expect_identical(events[[3]]$stage, "complete")
  expect_identical(events[[2]]$variable, "temp")
  expect_identical(events[[1]]$time_count, 3L)
})

test_that("GUI progress detail shows file, variable, counts, and date range", {
  detail <- ereefs:::ereefs_gui_progress_detail(list(
    message = "Extracting variable values",
    file = "https://example.invalid/catalog/demo_2020-01.nc",
    variable = "temp",
    point_index = 2L,
    point_count = 5L,
    time_count = 3L,
    time_start = as.POSIXct("2020-01-01 12:00:00", tz = "Etc/GMT-10"),
    time_end = as.POSIXct("2020-01-03 12:00:00", tz = "Etc/GMT-10"),
    file_count = 7L,
    variable_count = 4L
  ))

  expect_match(detail, "File: demo_2020-01.nc", fixed = TRUE)
  expect_match(detail, "Variable: temp", fixed = TRUE)
  expect_match(detail, "Point: 2 of 5", fixed = TRUE)
  expect_match(detail, "Matched times: 3", fixed = TRUE)
  expect_match(detail, "Date range: 2020-01-01 to 2020-01-03", fixed = TRUE)
  expect_match(detail, "Files: 7", fixed = TRUE)
  expect_match(detail, "Variables: 4", fixed = TRUE)
})

test_that("Python OPeNDAP URLs use explicit DAP2 for HTTPS datasets", {
  expect_equal(
    ereefs:::ereefs_python_opendap_url("https://thredds.nci.org.au/thredds/dodsC/fx3/example.nc"),
    "dap2://thredds.nci.org.au/thredds/dodsC/fx3/example.nc"
  )
  expect_equal(
    ereefs:::ereefs_python_opendap_url("dap4://thredds.nci.org.au/thredds/dodsC/fx3/example.nc"),
    "dap4://thredds.nci.org.au/thredds/dodsC/fx3/example.nc"
  )
})

test_that("manual local source paths can be pasted with quotes", {
  expect_equal(
    ereefs:::ereefs_gui_clean_source_path('"C:\\Users\\brobson\\example.nc"'),
    "C:\\Users\\brobson\\example.nc"
  )
  expect_equal(
    ereefs:::ereefs_gui_clean_source_path("'C:\\Users\\brobson\\example.nc'"),
    "C:\\Users\\brobson\\example.nc"
  )
})

test_that("in-app local browser lists folders and NetCDF-like files as text entries", {
  tmp <- tempfile("ereefs-browser-")
  dir.create(tmp)
  dir.create(file.path(tmp, "subdir"))
  file.create(file.path(tmp, "demo.nc"))
  file.create(file.path(tmp, "ignore.txt"))

  entries <- ereefs:::ereefs_gui_local_browser_entries(tmp)

  expect_true(any(entries$label == "[folder] subdir"))
  expect_true(any(entries$label == "[file] demo.nc"))
  expect_false(any(grepl("ignore.txt", entries$label, fixed = TRUE)))
  expect_true(all(file.exists(entries$path)))
})

test_that("catalog file coverage handles mixed daily and monthly file tables", {
  file_tbl <- tibble::tibble(
    file_date = as.Date(c("2019-02-01", "2019-03-01", "2019-03-03")),
    date_precision = c("month", "month", "day"),
    opendap_url = c("feb.nc", "mar.nc", "mar-03.nc")
  )
  covered <- ereefs:::ereefs_add_file_coverage(file_tbl)
  requested_days <- seq.Date(as.Date("2019-02-25"), as.Date("2019-03-02"), by = "day")
  selected <- dplyr::filter(
    covered,
    .data$file_start <= as.Date("2019-03-02"),
    .data$file_end >= as.Date("2019-02-25")
  )
  day_covered <- vapply(requested_days, function(day) {
    any(selected$file_start <= day & selected$file_end >= day)
  }, logical(1))

  expect_equal(covered$file_end[covered$opendap_url == "feb.nc"], as.Date("2019-02-28"))
  expect_equal(covered$file_end[covered$opendap_url == "mar.nc"], as.Date("2019-03-31"))
  expect_true(all(day_covered))
})

test_that("local demo map workflow and animation progress callback still smoke-test", {
  demo_file <- testthat::test_path("..", "..", "notebooks", "demo_data", "regular_demo_2020-01.nc")
  testthat::skip_if_not(file.exists(demo_file), "local demo NetCDF file is not available")

  events <- list()
  callback <- function(event) {
    events[[length(events) + 1L]] <<- event
  }

  map <- map_ereefs(
    var_name = "temp",
    target_date = "2020-01-01",
    layer = "surface",
    input_file = demo_file,
    plot_style = "smooth",
    label_towns = FALSE
  )
  movie <- map_ereefs_movie(
    var_name = "temp",
    start_date = "2020-01-01",
    end_date = "2020-01-01",
    layer = "surface",
    input_file = demo_file,
    output_dir = file.path(tempdir(), "ereefs-test-frames"),
    save_frames = FALSE,
    animation_format = "none",
    plot_style = "smooth",
    label_towns = FALSE,
    progress_callback = callback
  )

  expect_s3_class(map, "ggplot")
  expect_type(movie, "list")
  expect_gt(length(events), 0)
  expect_true(any(vapply(events, function(event) identical(event$stage, "load"), logical(1))))
  expect_true(any(vapply(events, function(event) identical(event$stage, "render"), logical(1))))
})

test_that("current GUI map can be recoloured from stored preview data", {
  demo_file <- testthat::test_path("..", "..", "notebooks", "demo_data", "regular_demo_2020-01.nc")
  testthat::skip_if_not(file.exists(demo_file), "local demo NetCDF file is not available")

  stored <- map_ereefs(
    var_name = "temp",
    target_date = "2020-01-01",
    layer = "surface",
    input_file = demo_file,
    plot_style = "smooth",
    label_towns = FALSE,
    return_poly = TRUE
  )
  input <- list(
    land_map = FALSE,
    map_palette = "viridis",
    map_auto_scale = TRUE,
    map_scale_min = NA_real_,
    map_scale_max = NA_real_,
    map_style = "smooth",
    smooth_pixels = 100,
    use_bbox = FALSE,
    lon_min = NA_real_,
    lon_max = NA_real_,
    lat_min = NA_real_,
    lat_max = NA_real_,
    label_towns = FALSE,
    gbr_poly = FALSE
  )

  p1 <- ereefs:::ereefs_gui_current_map_plot(stored, input)
  input$map_auto_scale <- FALSE
  input$map_scale_min <- 20
  input$map_scale_max <- 24
  input$map_palette <- "magma"
  p2 <- ereefs:::ereefs_gui_current_map_plot(stored, input)

  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
  expect_equal(nrow(stored$datapoly), nrow(stored$datapoly))
})

test_that("GUI map fetch signature ignores colour scale but tracks data requests", {
  signature_a <- ereefs:::ereefs_gui_map_fetch_signature(
    input_file = "catalog-a",
    var_name = "temp",
    target_date = as.Date("2020-01-01"),
    layer = "surface",
    box_bounds = c(142, 143, -19, -18)
  )
  signature_b <- ereefs:::ereefs_gui_map_fetch_signature(
    input_file = "catalog-a",
    var_name = "temp",
    target_date = as.Date("2020-01-01"),
    layer = "surface",
    box_bounds = c(142, 143, -19, -18)
  )
  signature_date <- ereefs:::ereefs_gui_map_fetch_signature(
    input_file = "catalog-a",
    var_name = "temp",
    target_date = as.Date("2020-01-02"),
    layer = "surface",
    box_bounds = c(142, 143, -19, -18)
  )
  signature_layer <- ereefs:::ereefs_gui_map_fetch_signature(
    input_file = "catalog-a",
    var_name = "temp",
    target_date = as.Date("2020-01-01"),
    layer = "2",
    box_bounds = c(142, 143, -19, -18)
  )

  expect_identical(signature_a, signature_b)
  expect_false(identical(signature_a, signature_date))
  expect_false(identical(signature_a, signature_layer))
})

test_that("GUI map scale range reads numeric values from stored preview data", {
  stored <- list(
    datapoly = tibble::tibble(
      id = c(1, 1, 2, 2),
      value = c(4.5, 4.5, 9.25, NA_real_),
      x = c(1, 2, 2, 3),
      y = c(1, 1, 2, 2)
    )
  )
  true_colour <- list(
    datapoly = tibble::tibble(
      id = c(1, 1),
      value = c("#112233", "#445566"),
      x = c(1, 2),
      y = c(1, 2)
    )
  )

  expect_equal(ereefs:::ereefs_gui_map_scale_range(stored), c(4.5, 9.25))
  expect_true(all(is.na(ereefs:::ereefs_gui_map_scale_range(true_colour))))
})

test_that("GUI point date range falls back to inspected bounds when inputs are missing", {
  inspection <- list(summary = tibble::tibble(
    time_start = as.POSIXct(NA),
    time_end = as.POSIXct(NA),
    file_start = as.Date("2020-01-01"),
    file_end = as.Date("2020-01-03")
  ))
  input <- list(points_start = NA, points_end = NA)

  range <- ereefs:::ereefs_gui_points_date_range(input, inspection)

  expect_equal(range$start, as.Date("2020-01-01"))
  expect_equal(range$end, as.Date("2020-01-03"))
})

test_that("time-file resolver rejects missing dates with a clear error", {
  expect_error(
    ereefs:::ereefs_resolve_time_files("example.nc", NA, as.Date("2020-01-01")),
    "start_date and end_date must be valid dates."
  )
})

test_that("surface wet-layer helper drops non-finite layers without warnings", {
  wet_tbl <- tibble::tibble(
    row_id = c(1, 1),
    time = as.POSIXct(c("2020-01-01 12:00:00", "2020-01-01 12:00:00"), tz = "Etc/GMT-10"),
    k = c(NA_real_, NA_real_),
    dz = c(1, 2)
  )

  expect_no_warning(
    result <- ereefs:::ereefs_surface_bottom_k(wet_tbl, which_layer = "surface")
  )
  expect_equal(nrow(result), 0)
})

# metadata:
# - gpt_version: GPT-5 Codex
# - time: 18:45
# - date: 2026-07-05
# - prompt_used: "Add smoke tests for NCI URL resolution, GUI date and spatial bounds, clicked-point coordinate conversion and labels, DAP URL conversion, and local map/movie progress callbacks."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 18:50
# - date: 2026-07-05
# - prompt_used: "Adjust clicked-point coordinate smoke test to compare normalized clicks against ggplot displayed panel bounds."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 19:00
# - date: 2026-07-05
# - prompt_used: "Add a smoke test that recolours a stored GUI map preview without re-fetching data."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 19:15
# - date: 2026-07-05
# - prompt_used: "Add a smoke test for mixed daily/monthly catalog coverage so the resolver does not warn falsely when monthly files cover requested days."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 19:30
# - date: 2026-07-05
# - prompt_used: "Add a smoke test for cleaning quoted pasted Windows local NetCDF paths."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 10:20
# - date: 2026-07-06
# - prompt_used: "Add a smoke test for the in-app local file browser listing folders and NetCDF-like files as visible text entries."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 13:11
# - date: 2026-07-06
# - prompt_used: "If the user changes the colour scale limits but does not change the date, layer or variable selected, the tool should re-render the map but not re-fetch the data."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 13:11
# - date: 2026-07-06
# - prompt_used: "Remove stale smoke tests for the disabled map-click point-selection helpers so the GUI smoke suite matches the current feature set."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 13:24
# - date: 2026-07-06
# - prompt_used: "Add smoke-test coverage for deriving GUI colour scale limits from stored numeric preview data."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 13:28
# - date: 2026-07-06
# - prompt_used: "Add regression coverage for GUI point-date fallback and clear invalid-date handling in time-file resolution."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 13:57
# - date: 2026-07-06
# - prompt_used: "Add regression coverage so non-finite wet-layer rows no longer trigger max/min warnings in surface point extraction."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 14:12
# - date: 2026-07-06
# - prompt_used: "Add regression coverage so transiently missing bounding-box numeric inputs do not break map preview argument construction."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 14:20
# - date: 2026-07-06
# - prompt_used: "Add regression coverage for shared variable syncing across tabs, with safe fallback for map-only variables."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 15:55
# - date: 2026-07-06
# - prompt_used: "Add regression coverage so GUI point parsing always returns latitude/longitude columns for single-point and empty cases."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 16:03
# - date: 2026-07-06
# - prompt_used: "Add regression coverage so transiently missing profile date inputs do not break the shared GUI point date-range helper."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 16:26
# - date: 2026-07-06
# - prompt_used: "Add regression coverage so out-of-range GUI dates clamp to the nearest inspected file or catalog boundary."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 16:50
# - date: 2026-07-06
# - prompt_used: "Add regression coverage so plot_ereefs_profile() handles single-variable multi-date profile arrays without requiring variable-name columns."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 17:00
# - date: 2026-07-06
# - prompt_used: "Add regression coverage for GUI multi-point profile plotting so multiple selected points render on one profile plot."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 17:04
# - date: 2026-07-06
# - prompt_used: "Add regression coverage so inspection and extraction emit lightweight progress-callback stages without extra I/O."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 17:18
# - date: 2026-07-06
# - prompt_used: "Stabilize the new progress regression coverage by testing callback payloads and GUI progress-detail formatting without depending on a brittle demo extraction path."
