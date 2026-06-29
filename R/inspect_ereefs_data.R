#' Inspect available variables, dimensions, dates, and spatial coverage
#'
#' `inspect_ereefs_data()` gives a lightweight overview of an eReefs or
#' CSIRO-EMS NetCDF dataset before extracting large arrays. It is intended for
#' interactive exploration, especially when a user is deciding which variables,
#' dates, or grid product to use.
#'
#' The function reads NetCDF metadata, coordinate variables, and time
#' coordinates, but it does not read the full data arrays for model variables.
#' When `input_file` is a THREDDS catalog, the catalog is used to summarize the
#' available files and the first matching NetCDF file is used as a representative
#' source for variable and dimension metadata.
#'
#' @param input_file NetCDF file path, OPeNDAP URL, THREDDS catalog URI, or a
#'   legacy shortcut such as `"catalog"`, `"menu"`, or `"nci"`.
#' @param recurse_catalog Logical; whether to recurse into child catalogs when
#'   `input_file` is a THREDDS catalog. Defaults to `FALSE` for speed. If a
#'   catalog has no direct NetCDF entries, the function automatically tries a
#'   recursive lookup.
#' @param include_files Logical; whether to include the catalog/file table in
#'   the returned object when available.
#' @param max_files Maximum number of catalog/file rows to include in the
#'   returned `files` tibble. The summary still reports the full count.
#' @return An object of class `"ereefs_data_inspection"`, which is a named list
#'   containing:
#'   \describe{
#'     \item{summary}{A one-row tibble with source, grid, spatial, temporal,
#'       and variable-count information.}
#'     \item{variables}{A tibble of NetCDF variables with names, units,
#'       long names, standard names, dimension names, dimension roles, and an
#'       `is_data_variable` flag.}
#'     \item{dimensions}{A tibble of NetCDF dimensions and inferred roles.}
#'     \item{files}{A tibble of catalog/file entries when available and
#'       `include_files = TRUE`.}
#'   }
#' @export
#' @examples
#' inspect_ereefs_data("notebooks/demo_data/regular_demo_2020-01.nc")
#' \dontrun{
#' inspect_ereefs_data(
#'   "https://thredds.nci.org.au/thredds/catalog/fx3/gbr4_H4p0_ABARRAr2_OBRAN2020_FG2Gv3_Dhnd/catalog.xml"
#' )
#' }
inspect_ereefs_data <- function(input_file = "catalog",
                                recurse_catalog = FALSE,
                                include_files = TRUE,
                                max_files = 50) {
  requested_input <- ereefs_prepare_input_reference(input_file)
  files <- ereefs_inspect_files(requested_input, recurse_catalog = recurse_catalog)
  source_file <- if (nrow(files) > 0 && "opendap_url" %in% names(files)) {
    files$opendap_url[[1]]
  } else {
    requested_input[[1]]
  }

  dimensions <- ereefs_inspect_dimensions(source_file)
  variable_names <- ereefs_var_names(source_file)
  variables <- ereefs_inspect_variables(source_file, variable_names, dimensions)
  spatial <- ereefs_inspect_spatial(source_file)
  temporal <- ereefs_inspect_temporal(source_file, files)

  summary <- dplyr::tibble(
    requested_input = paste(as.character(input_file), collapse = ", "),
    inspected_source = source_file,
    source_type = if (length(requested_input) == 1L && ereefs_is_catalog_request(requested_input)) "catalog" else "netcdf",
    file_count = nrow(files),
    file_start = temporal$file_start,
    file_end = temporal$file_end,
    time_start = temporal$time_start,
    time_end = temporal$time_end,
    time_step = temporal$time_step,
    longitude_min = spatial$longitude_min,
    longitude_max = spatial$longitude_max,
    latitude_min = spatial$latitude_min,
    latitude_max = spatial$latitude_max,
    grid_type = spatial$grid_type,
    i_count = spatial$i_count,
    j_count = spatial$j_count,
    k_count = spatial$k_count,
    variable_count = nrow(variables),
    data_variable_count = sum(variables$is_data_variable, na.rm = TRUE)
  )

  out <- list(
    summary = summary,
    variables = variables,
    dimensions = dimensions,
    files = if (isTRUE(include_files)) utils::head(files, max_files) else dplyr::tibble()
  )
  class(out) <- c("ereefs_data_inspection", class(out))
  out
}

#' @export
print.ereefs_data_inspection <- function(x, ...) {
  cat("eReefs data inspection\n")
  print(x$summary, ...)
  cat("\nVariables:\n")
  print(utils::head(x$variables, 12), ...)
  if (nrow(x$variables) > 12) {
    cat("... ", nrow(x$variables) - 12, " more variables\n", sep = "")
  }
  if (!is.null(x$files) && nrow(x$files) > 0) {
    cat("\nFiles:\n")
    print(utils::head(x$files, 8), ...)
    if (nrow(x$files) > 8) {
      cat("... ", nrow(x$files) - 8, " more listed files\n", sep = "")
    }
  }
  invisible(x)
}

ereefs_inspect_files <- function(input_file, recurse_catalog = FALSE) {
  if (length(input_file) == 1L && ereefs_is_catalog_request(input_file)) {
    entries <- tryCatch(
      ereefs_catalog_entries(input_file, recurse = recurse_catalog),
      error = function(e) dplyr::tibble()
    )
    if (!nrow(entries) && !isTRUE(recurse_catalog)) {
      entries <- tryCatch(
        ereefs_catalog_entries(input_file, recurse = TRUE),
        error = function(e) dplyr::tibble()
      )
    }
    return(entries %>% dplyr::arrange(.data$file_date, .data$name))
  }

  dplyr::tibble(
    name = basename(as.character(input_file)),
    url_path = NA_character_,
    opendap_url = as.character(input_file),
    file_date = as.Date(unlist(lapply(as.character(input_file), ereefs_extract_file_date)), origin = "1970-01-01"),
    date_precision = unlist(lapply(as.character(input_file), ereefs_extract_file_date_precision))
  )
}

ereefs_inspect_dimensions <- function(input_file) {
  if (ereefs_is_remote_file(input_file) && ereefs_python_ready()) {
    ds <- ereefs_python_dataset(input_file)
    sizes <- reticulate::py_to_r(ds$sizes)
    names <- names(sizes)
    return(dplyr::tibble(
      name = names,
      length = as.integer(unlist(sizes, use.names = FALSE)),
      role = vapply(names, ereefs_dim_role, character(1))
    ))
  }

  dims <- tryCatch(ncmeta::nc_dims(input_file), error = function(e) dplyr::tibble())
  if (!nrow(dims)) {
    return(dplyr::tibble(name = character(), length = integer(), role = character()))
  }
  dplyr::tibble(
    name = dims$name,
    length = as.integer(dims$length),
    role = vapply(dims$name, ereefs_dim_role, character(1))
  )
}

ereefs_inspect_variables <- function(input_file, variable_names, dimensions) {
  empty_variables <- dplyr::tibble(
    variable = character(),
    units = character(),
    long_name = character(),
    standard_name = character(),
    dimensions = character(),
    dimension_roles = character(),
    n_dimensions = integer(),
    is_data_variable = logical()
  )
  if (!length(variable_names)) {
    return(empty_variables)
  }

  coordinate_names <- unique(c(
    dimensions$name,
    "latitude", "longitude", "lat", "lon",
    "x_centre", "y_centre", "x_grid", "y_grid", "z_grid", "zc"
  ))

  rows <- lapply(variable_names, function(var_name) {
    dims <- tryCatch(ereefs_var_dims(input_file, var_name), error = function(e) dplyr::tibble())
    dim_names <- if (nrow(dims)) dims$name else character()
    dim_roles <- if (nrow(dims)) dims$role else character()
    dplyr::tibble(
      variable = var_name,
      units = ereefs_safe_var_attr(input_file, var_name, "units"),
      long_name = ereefs_safe_var_attr(input_file, var_name, "long_name"),
      standard_name = ereefs_safe_var_attr(input_file, var_name, "standard_name"),
      dimensions = paste(dim_names, collapse = ","),
      dimension_roles = paste(stats::na.omit(dim_roles), collapse = ","),
      n_dimensions = length(dim_names),
      is_data_variable = !(var_name %in% coordinate_names) && length(dim_names) > 0
    )
  })

  dplyr::bind_rows(c(list(empty_variables), rows)) %>%
    dplyr::arrange(!.data$is_data_variable, .data$variable)
}

ereefs_safe_var_attr <- function(input_file, var_name, attribute) {
  out <- tryCatch(ereefs_var_attr(input_file, var_name, attribute), error = function(e) NA_character_)
  if (length(out) != 1 || is.na(out) || !nzchar(out)) {
    return(NA_character_)
  }
  as.character(out)
}

ereefs_inspect_spatial <- function(input_file) {
  spatial_grid <- tryCatch(ereefs_spatial_grid(input_file), error = function(e) dplyr::tibble())
  if (!nrow(spatial_grid)) {
    return(list(
      longitude_min = NA_real_,
      longitude_max = NA_real_,
      latitude_min = NA_real_,
      latitude_max = NA_real_,
      grid_type = NA_character_,
      i_count = NA_integer_,
      j_count = NA_integer_,
      k_count = NA_integer_
    ))
  }

  dims <- ereefs_inspect_dimensions(input_file)
  k_count <- dims$length[match("k", dims$role)]
  list(
    longitude_min = min(spatial_grid$longitude, na.rm = TRUE),
    longitude_max = max(spatial_grid$longitude, na.rm = TRUE),
    latitude_min = min(spatial_grid$latitude, na.rm = TRUE),
    latitude_max = max(spatial_grid$latitude, na.rm = TRUE),
    grid_type = if (ereefs_is_regular_centre_grid(spatial_grid)) "regular" else "curvilinear",
    i_count = length(unique(spatial_grid$i)),
    j_count = length(unique(spatial_grid$j)),
    k_count = if (length(k_count) && !is.na(k_count)) as.integer(k_count) else NA_integer_
  )
}

ereefs_inspect_temporal <- function(input_file, files) {
  timing <- tryCatch(get_origin_and_times(input_file), error = function(e) NULL)
  times <- if (is.null(timing)) as.POSIXct(character(), tz = "Etc/GMT-10") else timing[[2]]
  file_dates <- files$file_date[!is.na(files$file_date)]

  list(
    file_start = if (length(file_dates)) min(file_dates) else as.Date(NA),
    file_end = if (length(file_dates)) max(file_dates) else as.Date(NA),
    time_start = if (length(times)) min(times) else as.POSIXct(NA, tz = "Etc/GMT-10"),
    time_end = if (length(times)) max(times) else as.POSIXct(NA, tz = "Etc/GMT-10"),
    time_step = ereefs_describe_time_step(times)
  )
}

ereefs_describe_time_step <- function(times) {
  if (length(times) < 2) {
    return(NA_character_)
  }
  step_hours <- stats::median(as.numeric(diff(sort(times)), units = "hours"), na.rm = TRUE)
  if (!is.finite(step_hours)) {
    return(NA_character_)
  }
  if (step_hours < 48) {
    return(paste0(signif(step_hours, 4), " hours"))
  }
  paste0(signif(step_hours / 24, 4), " days")
}

# metadata:
# - gpt_version: GPT-5 Codex
# - time: 12:00
# - date: 2026-06-29
# - prompt_used: "Address GitHub issue #26 by adding inspect_ereefs_data() and update all documentation including the vignette."
