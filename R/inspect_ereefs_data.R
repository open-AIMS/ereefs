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
#' source for variable and dimension metadata. When `read_remote_metadata` is
#' `FALSE`, the function still reads variable and dimension metadata from that
#' representative file, but skips the heavier spatial and time-coordinate
#' metadata pass.
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
#' @param read_remote_metadata Logical; whether catalogue inspection should also
#'   read full remote spatial and time metadata from a representative
#'   NetCDF/OPeNDAP file. Set to `FALSE` to keep catalog inspection lighter
#'   while still listing variables and dimensions from the representative file.
#' @param progress_callback Optional function called with lightweight progress
#'   details such as the current stage, representative file, or catalog file
#'   count while inspection is running.
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
                                max_files = 50,
                                read_remote_metadata = TRUE,
                                progress_callback = NULL) {
  ereefs_inspect_progress(
    progress_callback,
    stage = "start",
    message = "Resolving selected source",
    requested_input = input_file
  )
  requested_input <- ereefs_prepare_input_reference(input_file)
  ereefs_inspect_progress(
    progress_callback,
    stage = "files",
    message = "Reading catalog or file listing",
    requested_input = requested_input
  )
  files <- ereefs_inspect_files(requested_input, recurse_catalog = recurse_catalog)
  is_catalog <- length(requested_input) == 1L && ereefs_is_catalog_request(requested_input)
  source_file <- if (nrow(files) > 0 && "opendap_url" %in% names(files)) {
    files$opendap_url[[1]]
  } else {
    requested_input[[1]]
  }
  ereefs_inspect_progress(
    progress_callback,
    stage = "representative",
    message = "Using representative file for variable metadata",
    requested_input = requested_input,
    source_file = source_file,
    file_count = nrow(files)
  )

  if (isTRUE(is_catalog) && !isTRUE(read_remote_metadata)) {
    ereefs_inspect_progress(
      progress_callback,
      stage = "dimensions",
      message = "Reading dimensions",
      source_file = source_file,
      file_count = nrow(files)
    )
    dimensions <- ereefs_inspect_dimensions(source_file)
    variable_names <- ereefs_var_names(source_file)
    ereefs_inspect_progress(
      progress_callback,
      stage = "variables",
      message = "Reading variable list",
      source_file = source_file,
      variable_count = length(variable_names)
    )
    variables <- ereefs_inspect_variables(source_file, variable_names, dimensions)
    spatial <- ereefs_empty_spatial_summary()
    ereefs_inspect_progress(
      progress_callback,
      stage = "time",
      message = "Summarising time coverage from file listing",
      file_count = nrow(files)
    )
    temporal <- ereefs_inspect_temporal_from_files(files)
  } else {
    ereefs_inspect_progress(
      progress_callback,
      stage = "dimensions",
      message = "Reading dimensions",
      source_file = source_file,
      file_count = nrow(files)
    )
    dimensions <- ereefs_inspect_dimensions(source_file)
    variable_names <- ereefs_var_names(source_file)
    ereefs_inspect_progress(
      progress_callback,
      stage = "variables",
      message = "Reading variable list",
      source_file = source_file,
      variable_count = length(variable_names)
    )
    variables <- ereefs_inspect_variables(source_file, variable_names, dimensions)
    ereefs_inspect_progress(
      progress_callback,
      stage = "spatial",
      message = "Reading spatial and grid metadata",
      source_file = source_file
    )
    spatial <- ereefs_inspect_spatial(source_file)
    ereefs_inspect_progress(
      progress_callback,
      stage = "time",
      message = "Reading temporal coverage",
      source_file = source_file,
      file_count = nrow(files)
    )
    temporal <- ereefs_inspect_temporal(source_file, files)
  }

  summary <- dplyr::tibble(
    requested_input = paste(as.character(input_file), collapse = ", "),
    inspected_source = source_file,
    source_type = if (is_catalog) "catalog" else "netcdf",
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
  ereefs_inspect_progress(
    progress_callback,
    stage = "complete",
    message = "Inspection complete",
    source_file = source_file,
    file_count = nrow(files),
    variable_count = nrow(variables)
  )
  out
}

ereefs_inspect_progress <- function(progress_callback, ...) {
  if (is.null(progress_callback) || !is.function(progress_callback)) {
    return(invisible(NULL))
  }
  tryCatch(progress_callback(list(...)), error = function(e) NULL)
  invisible(NULL)
}

ereefs_empty_dimension_table <- function() {
  dplyr::tibble(name = character(), length = integer(), role = character())
}

ereefs_empty_variable_table <- function() {
  dplyr::tibble(
    variable = character(),
    units = character(),
    long_name = character(),
    standard_name = character(),
    dimensions = character(),
    dimension_roles = character(),
    n_dimensions = integer(),
    is_data_variable = logical()
  )
}

ereefs_empty_spatial_summary <- function() {
  list(
    longitude_min = NA_real_,
    longitude_max = NA_real_,
    latitude_min = NA_real_,
    latitude_max = NA_real_,
    grid_type = NA_character_,
    i_count = NA_integer_,
    j_count = NA_integer_,
    k_count = NA_integer_
  )
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
    entries <- dplyr::tibble()
    if (isTRUE(recurse_catalog)) {
      entries <- tryCatch(
        ereefs_catalog_entries(input_file, recurse = FALSE),
        error = function(e) dplyr::tibble()
      )
      if (!nrow(entries)) {
        entries <- tryCatch(
          ereefs_catalog_entries(input_file, recurse = TRUE),
          error = function(e) dplyr::tibble()
        )
      }
    } else {
      entries <- tryCatch(
        ereefs_catalog_entries(input_file, recurse = FALSE),
        error = function(e) dplyr::tibble()
      )
      if (!nrow(entries)) {
        entries <- tryCatch(
          ereefs_catalog_entries(input_file, recurse = TRUE),
          error = function(e) dplyr::tibble()
        )
      }
    }
    return(ereefs_inspect_normalise_file_table(entries))
  }

  ereefs_inspect_normalise_file_table(dplyr::tibble(
    name = basename(as.character(input_file)),
    url_path = NA_character_,
    opendap_url = as.character(input_file),
    file_date = as.Date(unlist(lapply(as.character(input_file), ereefs_extract_file_date)), origin = "1970-01-01"),
    date_precision = unlist(lapply(as.character(input_file), ereefs_extract_file_date_precision))
  ))
}

ereefs_inspect_normalise_file_table <- function(entries) {
  if (is.null(entries) || !nrow(entries)) {
    return(dplyr::tibble(
      name = character(),
      url_path = character(),
      opendap_url = character(),
      file_date = as.Date(character()),
      date_precision = character()
    ))
  }
  if (!("name" %in% names(entries))) {
    entries$name <- if ("opendap_url" %in% names(entries)) basename(entries$opendap_url) else NA_character_
  }
  if (!("url_path" %in% names(entries))) {
    entries$url_path <- NA_character_
  }
  if (!("opendap_url" %in% names(entries))) {
    entries$opendap_url <- NA_character_
  }
  if (!("file_date" %in% names(entries))) {
    entries$file_date <- as.Date(unlist(lapply(entries$name, ereefs_extract_file_date)), origin = "1970-01-01")
  }
  if (!("date_precision" %in% names(entries))) {
    entries$date_precision <- unlist(lapply(entries$name, ereefs_extract_file_date_precision))
  }
  entries %>%
    dplyr::mutate(file_date = as.Date(.data$file_date, origin = "1970-01-01")) %>%
    dplyr::arrange(.data$file_date, .data$name)
}

ereefs_inspect_dimensions <- function(input_file) {
  if (ereefs_is_remote_file(input_file) && ereefs_python_ready()) {
    ds <- ereefs_python_dataset(input_file)
    sizes <- ereefs_python_mapping_to_named_integer(ds$sizes)
    names <- names(sizes)
    return(dplyr::tibble(
      name = names,
      length = as.integer(sizes),
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

ereefs_python_mapping_to_named_integer <- function(mapping) {
  converted <- tryCatch(reticulate::py_to_r(mapping), error = function(e) NULL)
  if (!is.null(converted) && !is.null(names(converted))) {
    values <- tryCatch(as.integer(unlist(converted, use.names = FALSE)), error = function(e) integer())
    if (length(values) == length(names(converted))) {
      return(stats::setNames(values, names(converted)))
    }
  }

  items <- tryCatch(
    reticulate::iterate(mapping$items(), function(item) reticulate::py_to_r(item)),
    error = function(e) list()
  )
  if (!length(items)) {
    return(stats::setNames(integer(), character()))
  }
  names <- vapply(items, function(item) as.character(item[[1]]), character(1))
  values <- vapply(items, function(item) as.integer(item[[2]]), integer(1))
  stats::setNames(values, names)
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

  if (ereefs_is_remote_file(input_file) && ereefs_python_ready()) {
    return(ereefs_inspect_variables_remote(input_file, variable_names, dimensions, empty_variables))
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

ereefs_inspect_variables_remote <- function(input_file, variable_names, dimensions, empty_variables) {
  ds <- ereefs_python_dataset(input_file)
  coordinate_names <- unique(c(
    dimensions$name,
    "latitude", "longitude", "lat", "lon",
    "x_centre", "y_centre", "x_grid", "y_grid", "z_grid", "zc"
  ))

  rows <- lapply(variable_names, function(var_name) {
    da <- ds[[var_name]]
    dim_names <- unlist(reticulate::py_to_r(da$dims), use.names = FALSE)
    dim_roles <- vapply(dim_names, ereefs_dim_role, character(1))
    attrs <- tryCatch(reticulate::py_to_r(da$attrs), error = function(e) list())
    dplyr::tibble(
      variable = var_name,
      units = ereefs_remote_attr_value(attrs, "units"),
      long_name = ereefs_remote_attr_value(attrs, "long_name"),
      standard_name = ereefs_remote_attr_value(attrs, "standard_name"),
      dimensions = paste(dim_names, collapse = ","),
      dimension_roles = paste(stats::na.omit(dim_roles), collapse = ","),
      n_dimensions = length(dim_names),
      is_data_variable = !(var_name %in% coordinate_names) && length(dim_names) > 0
    )
  })

  dplyr::bind_rows(c(list(empty_variables), rows)) %>%
    dplyr::arrange(!.data$is_data_variable, .data$variable)
}

ereefs_remote_attr_value <- function(attrs, attribute) {
  value <- attrs[[attribute]]
  if (is.null(value) || length(value) != 1L || is.na(value) || !nzchar(as.character(value))) {
    return(NA_character_)
  }
  as.character(value)
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
  file_temporal <- ereefs_inspect_temporal_from_files(files)

  list(
    file_start = file_temporal$file_start,
    file_end = file_temporal$file_end,
    time_start = if (length(times)) min(times) else as.POSIXct(NA, tz = "Etc/GMT-10"),
    time_end = if (length(times)) max(times) else as.POSIXct(NA, tz = "Etc/GMT-10"),
    time_step = ereefs_describe_time_step(times)
  )
}

ereefs_inspect_temporal_from_files <- function(files) {
  file_dates <- files$file_date[!is.na(files$file_date)]

  list(
    file_start = if (length(file_dates)) min(file_dates) else as.Date(NA),
    file_end = if (length(file_dates)) max(file_dates) else as.Date(NA),
    time_start = as.POSIXct(NA, tz = "Etc/GMT-10"),
    time_end = as.POSIXct(NA, tz = "Etc/GMT-10"),
    time_step = NA_character_
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
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 16:18
# - date: 2026-07-05
# - prompt_used: "Fix GUI-triggered NCI catalog inspection when catalog file tables lack file_date or related columns."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 16:13
# - date: 2026-07-05
# - prompt_used: "Harden remote inspection for xarray mapping-style dimension sizes returned by PyDAP-backed datasets."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 16:13
# - date: 2026-07-05
# - prompt_used: "Only trust Python mapping conversion when dimension sizes unlist to integers; otherwise fall back to mapping items."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 16:26
# - date: 2026-07-05
# - prompt_used: "Finalize verified live NCI catalog inspection fixes for missing file_date columns and PyDAP-backed xarray dimension mappings."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 16:40
# - date: 2026-07-05
# - prompt_used: "Prevent GUI hangs on NCI catalog inspection by allowing fast catalog-only inspection without opening a representative remote OPeNDAP file."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 12:28
# - date: 2026-07-06
# - prompt_used: "Always load representative-file variable and dimension metadata during catalog inspection so NCI catalog variables are available on the inspection page and in GUI variable menus, while keeping full remote grid/time metadata optional."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 12:37
# - date: 2026-07-06
# - prompt_used: "Speed up NCI catalog inspection by trying direct catalog file entries first and only recursing into child catalogs when the direct catalog does not expose files."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 12:46
# - date: 2026-07-06
# - prompt_used: "Speed up representative remote-file variable inspection by collecting dimensions and key attributes in one dataset-backed pass instead of re-querying each remote variable attribute separately."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 17:04
# - date: 2026-07-06
# - prompt_used: "Show richer progress detail during inspections without extra I/O, including current stage, representative file, and catalog file counts."
