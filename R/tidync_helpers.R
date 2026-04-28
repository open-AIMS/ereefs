# Internal helpers for tidync-first dataset access and geometry handling.

ereefs_dim_aliases <- list(
  i = c("i", "i_centre", "i_grid", "longitude", "lon"),
  j = c("j", "j_centre", "j_grid", "latitude", "lat"),
  k = c("k", "k_centre", "k_grid"),
  time = c("time", "t", "record", "index")
)

ereefs_dim_role <- function(name) {
  for (role in names(ereefs_dim_aliases)) {
    if (name %in% ereefs_dim_aliases[[role]]) {
      return(role)
    }
  }
  NA_character_
}

ereefs_normalise_dim_table <- function(dim_table) {
  dim_table$role <- vapply(dim_table$name, ereefs_dim_role, character(1))
  dim_table
}

ereefs_normalise_columns <- function(tbl) {
  rename_map <- c(
    latitude = "latitude",
    longitude = "longitude",
    x_centre = "longitude",
    y_centre = "latitude",
    i = "i",
    i_centre = "i",
    i_grid = "i",
    j = "j",
    j_centre = "j",
    j_grid = "j",
    k = "k",
    k_centre = "k",
    k_grid = "k",
    time = "time",
    t = "time",
    record = "time",
    index = "time"
  )

  for (old_name in names(rename_map)) {
    new_name <- rename_map[[old_name]]
    if (old_name %in% names(tbl) && !(new_name %in% names(tbl))) {
      names(tbl)[names(tbl) == old_name] <- new_name
    }
  }
  tbl
}

ereefs_activate_var <- function(input_file, var_name, select_var = var_name) {
  nc <- tidync::tidync(input_file) %>% tidync::activate(var_name)
  tidync::activate(nc, tidync::active(nc), select_var = select_var)
}

ereefs_apply_filters <- function(nc, filters) {
  if (length(filters) == 0) {
    return(nc)
  }

  for (filter_name in names(filters)) {
    filter_value <- filters[[filter_name]]
    if (length(filter_value) == 1) {
      expr <- rlang::expr(
        dplyr::between(!!rlang::sym(filter_name), !!filter_value, !!filter_value)
      )
    } else if (length(filter_value) == 0) {
      next
    } else {
      expr <- rlang::expr(
        dplyr::between(!!rlang::sym(filter_name), !!min(filter_value), !!max(filter_value))
      )
    }
    nc <- rlang::eval_bare(
      rlang::call2(tidync::hyper_filter, nc, !!!stats::setNames(list(expr), filter_name))
    )
  }
  nc
}

ereefs_read_var_array <- function(input_file, var_name, filters = list(), select_var = var_name) {
  if (ereefs_is_remote_file(input_file) && ereefs_python_ready()) {
    return(reticulate::py_to_r(ereefs_python_data_array(input_file, select_var, filters = filters)$values))
  }
  nc <- ereefs_activate_var(input_file, var_name, select_var = select_var)
  nc <- ereefs_apply_filters(nc, filters)
  tidync::hyper_array(nc, select_var = select_var)[[1]]
}

ereefs_read_var_tibble <- function(input_file, var_name, filters = list(), select_var = var_name) {
  if (ereefs_is_remote_file(input_file) && ereefs_python_ready()) {
    da <- ereefs_python_data_array(input_file, select_var, filters = filters)
    tbl <- ereefs_normalise_columns(dplyr::as_tibble(reticulate::py_to_r(da$to_dataframe(name = select_var)$reset_index())))
  } else {
    nc <- ereefs_activate_var(input_file, var_name, select_var = select_var)
    nc <- ereefs_apply_filters(nc, filters)
    tbl <- ereefs_normalise_columns(tidync::hyper_tibble(nc, select_var = select_var))
  }

  for (filter_name in names(filters)) {
    role <- ereefs_dim_role(filter_name)
    out_name <- if (is.na(role)) filter_name else role
    if (!(out_name %in% names(tbl)) && length(filters[[filter_name]]) == 1L) {
      tbl[[out_name]] <- filters[[filter_name]][[1]]
    }
  }

  tbl
}

ereefs_var_dims <- function(input_file, var_name) {
  if (ereefs_is_remote_file(input_file) && ereefs_python_ready()) {
    da <- ereefs_python_dataset(input_file)[[var_name]]
    dim_names <- unlist(reticulate::py_to_r(da$dims), use.names = FALSE)
    dim_sizes <- reticulate::py_to_r(da$sizes)
    out <- dplyr::tibble(
      name = dim_names,
      length = unname(vapply(dim_names, function(x) as.integer(dim_sizes[[x]]), integer(1)))
    )
    return(ereefs_normalise_dim_table(out))
  }
  nc <- ereefs_activate_var(input_file, var_name)
  ereefs_normalise_dim_table(tidync::hyper_dims(nc))
}

ereefs_var_attr <- function(input_file, var_name, attribute) {
  if (ereefs_is_remote_file(input_file) && ereefs_python_ready()) {
    attrs <- reticulate::py_to_r(ereefs_python_dataset(input_file)[[var_name]]$attrs)
    return(if (!is.null(attrs[[attribute]])) as.character(attrs[[attribute]]) else NA_character_)
  }
  att <- tryCatch(
    ncmeta::nc_att(input_file, var_name, attribute),
    error = function(e) NULL
  )
  if (is.null(att)) {
    return(NA_character_)
  }
  if (nrow(att) == 0) {
    return(NA_character_)
  }
  att$value[[1]]
}

ereefs_mask_array_sentinels <- function(arr, input_file, var_name) {
  if (!is.numeric(arr)) {
    return(arr)
  }

  fill_value <- suppressWarnings(as.numeric(ereefs_var_attr(input_file, var_name, "_FillValue")))
  missing_value <- suppressWarnings(as.numeric(ereefs_var_attr(input_file, var_name, "missing_value")))
  sentinel_values <- unique(c(fill_value[is.finite(fill_value)], missing_value[is.finite(missing_value)]))
  if (length(sentinel_values)) {
    for (sentinel in sentinel_values) {
      arr[arr == sentinel] <- NA_real_
    }
  }

  arr[is.finite(arr) & abs(arr) >= 1e30] <- NA_real_
  arr
}

ereefs_var_names <- function(input_file) {
  if (ereefs_is_remote_file(input_file) && ereefs_python_ready()) {
    ds <- ereefs_python_dataset(input_file)
    builtins <- reticulate::import_builtins()
    var_keys <- reticulate::py_to_r(builtins$list(ds$variables$keys()))
    coord_keys <- reticulate::py_to_r(builtins$list(ds$coords$keys()))
    return(unique(c(as.character(var_keys), as.character(coord_keys))))
  }
  vars <- tryCatch(ncmeta::nc_vars(input_file)$name, error = function(e) character())
  dims <- tryCatch(ncmeta::nc_dims(input_file)$name, error = function(e) character())
  unique(c(vars, dims))
}

ereefs_catalog_url <- function() {
  "https://thredds.nci.org.au/thredds/catalog/catalogs/fx3/catalog.xml"
}

ereefs_canonicalise_url <- function(input_file) {
  if (!is.character(input_file) || length(input_file) != 1 || is.na(input_file)) {
    return(input_file)
  }

  input_file <- sub("^http://", "https://", input_file)
  input_file <- gsub("https://dapds00\\.nci\\.org\\.au", "https://thredds.nci.org.au", input_file)
  input_file <- gsub("/thredds/catalogs/", "/thredds/catalog/catalogs/", input_file)
  input_file
}

ereefs_is_remote_file <- function(input_file) {
  is.character(input_file) && length(input_file) == 1 && grepl("^https?://", input_file)
}

.ereefs_py_cache <- new.env(parent = emptyenv())
.ereefs_grid_cache <- new.env(parent = emptyenv())
.ereefs_warning_cache <- new.env(parent = emptyenv())

ereefs_python_configure <- function() {
  Sys.setenv(RETICULATE_USE_MANAGED_VENV = "false")
  reticulate::use_python(Sys.which("python"), required = TRUE)
}

ereefs_python_ready <- function() {
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    return(FALSE)
  }

  tryCatch({
    ereefs_python_configure()
    reticulate::import("xarray", delay_load = TRUE)
    reticulate::import("pydap", delay_load = TRUE)
    TRUE
  }, error = function(e) FALSE)
}

ereefs_python_dataset <- function(input_file) {
  if (!ereefs_python_ready()) {
    stop("Python remote backend requires the reticulate, xarray, and pydap packages.")
  }

  cache_key <- gsub("[^A-Za-z0-9]+", "_", input_file)
  if (exists(cache_key, envir = .ereefs_py_cache, inherits = FALSE)) {
    return(get(cache_key, envir = .ereefs_py_cache, inherits = FALSE))
  }

  ereefs_python_configure()
  xr <- reticulate::import("xarray", delay_load = TRUE)
  ds <- xr$open_dataset(input_file, engine = "pydap", decode_times = FALSE)
  assign(cache_key, ds, envir = .ereefs_py_cache)
  ds
}

ereefs_python_indexer <- function(values) {
  if (length(values) == 1) {
    start <- as.integer(values[[1]] - 1L)
    return(reticulate::import_builtins()$slice(start, start + 1L))
  }

  start <- as.integer(min(values) - 1L)
  stop <- as.integer(max(values))
  reticulate::import_builtins()$slice(start, stop)
}

ereefs_python_data_array <- function(input_file, var_name, filters = list()) {
  da <- ereefs_python_dataset(input_file)[[var_name]]
  if (!length(filters)) {
    return(da)
  }

  indexers <- lapply(filters, ereefs_python_indexer)
  do.call(da$isel, indexers)
}

ereefs_coord_var_names <- function(input_file) {
  var_names <- ereefs_var_names(input_file)

  if (all(c("longitude", "latitude") %in% var_names)) {
    return(list(longitude = "longitude", latitude = "latitude"))
  }

  if (all(c("x_centre", "y_centre") %in% var_names)) {
    return(list(longitude = "x_centre", latitude = "y_centre"))
  }

  stop("Could not identify horizontal coordinate variables in this dataset.")
}

ereefs_coord_table <- function(input_file) {
  coord_vars <- ereefs_coord_var_names(input_file)
  longitude <- ereefs_read_var_array(input_file, coord_vars$longitude)
  latitude <- ereefs_read_var_array(input_file, coord_vars$latitude)

  if (is.null(dim(longitude)) && is.null(dim(latitude))) {
    stop("Coordinate variables must be at least one-dimensional.")
  }

  if (length(dim(longitude)) <= 1 && length(dim(latitude)) <= 1) {
    longitude <- as.numeric(longitude)
    latitude <- as.numeric(latitude)
    out <- expand.grid(
      i = seq_along(longitude),
      j = seq_along(latitude)
    )
    out$longitude <- longitude[out$i]
    out$latitude <- latitude[out$j]
    return(dplyr::as_tibble(out))
  }

  if (length(dim(longitude)) == 2 && length(dim(latitude)) == 2) {
    dims <- dim(longitude)
    dim_roles <- unname(ereefs_var_dims(input_file, coord_vars$longitude)$role[1:2])
    if (identical(dim_roles, c("j", "i"))) {
      out <- expand.grid(
        j = seq_len(dims[1]),
        i = seq_len(dims[2])
      )
    } else {
      out <- expand.grid(
        i = seq_len(dims[1]),
        j = seq_len(dims[2])
      )
    }
    out$longitude <- as.vector(longitude)
    out$latitude <- as.vector(latitude)
    return(dplyr::as_tibble(out))
  }

  stop("Unsupported coordinate layout. Expected either 1D longitude/latitude axes or 2D centre coordinates.")
}

ereefs_matrix_from_coord_table <- function(coord_table, value_col) {
  i_vals <- sort(unique(coord_table$i))
  j_vals <- sort(unique(coord_table$j))
  out <- matrix(NA_real_, nrow = length(i_vals), ncol = length(j_vals))
  idx_i <- match(coord_table$i, i_vals)
  idx_j <- match(coord_table$j, j_vals)
  out[cbind(idx_i, idx_j)] <- coord_table[[value_col]]
  out
}

ereefs_role_array <- function(arr,
                              dims,
                              target_roles,
                              filters = list(),
                              drop_singleton = TRUE) {
  if (is.null(dim(arr))) {
    return(arr)
  }

  role_names <- ifelse(is.na(dims$role), dims$name, dims$role)
  eff_lengths <- dims$length
  if (length(filters)) {
    for (idx in seq_len(nrow(dims))) {
      dim_name <- dims$name[[idx]]
      if (dim_name %in% names(filters)) {
        eff_lengths[[idx]] <- length(filters[[dim_name]])
      }
    }
  }

  present_idx <- seq_along(role_names)
  if (length(dim(arr)) != length(present_idx)) {
    present_idx <- present_idx[eff_lengths > 1L]
  }
  if (length(dim(arr)) != length(present_idx)) {
    present_idx <- tail(present_idx, length(dim(arr)))
  }
  if (length(dim(arr)) != length(present_idx)) {
    stop("Could not reconcile array dimensions with declared variable roles.")
  }

  present_roles <- role_names[present_idx]
  keep_local <- match(target_roles, present_roles)
  keep_local <- keep_local[!is.na(keep_local)]
  if (!length(keep_local)) {
    return(arr)
  }

  perm <- c(keep_local, setdiff(seq_along(present_roles), keep_local))
  arr <- aperm(arr, perm)

  if (drop_singleton) {
    keep_dims <- dim(arr) > 1L
    if (!any(keep_dims)) {
      return(as.vector(arr))
    }
    arr <- array(arr, dim = dim(arr)[keep_dims])
  }

  arr
}

ereefs_regular_edges <- function(values) {
  values <- as.numeric(values)
  if (length(values) == 1) {
    return(c(values - 0.5, values + 0.5))
  }
  midpoints <- values[-length(values)] + diff(values) / 2
  c(values[1] - diff(values)[1] / 2, midpoints, values[length(values)] + diff(values)[length(values) - 1] / 2)
}

ereefs_centre_to_corner_matrix <- function(centre_matrix) {
  nr <- nrow(centre_matrix)
  nc <- ncol(centre_matrix)
  padded <- matrix(NA_real_, nrow = nr + 2L, ncol = nc + 2L)
  padded[2:(nr + 1L), 2:(nc + 1L)] <- centre_matrix

  padded[1, 2:(nc + 1L)] <- if (nr > 1L) 2 * centre_matrix[1, ] - centre_matrix[2, ] else centre_matrix[1, ]
  padded[nr + 2L, 2:(nc + 1L)] <- if (nr > 1L) 2 * centre_matrix[nr, ] - centre_matrix[nr - 1L, ] else centre_matrix[nr, ]
  padded[2:(nr + 1L), 1] <- if (nc > 1L) 2 * centre_matrix[, 1] - centre_matrix[, 2] else centre_matrix[, 1]
  padded[2:(nr + 1L), nc + 2L] <- if (nc > 1L) 2 * centre_matrix[, nc] - centre_matrix[, nc - 1L] else centre_matrix[, nc]

  padded[1, 1] <- padded[1, 2] + padded[2, 1] - padded[2, 2]
  padded[1, nc + 2L] <- padded[1, nc + 1L] + padded[2, nc + 2L] - padded[2, nc + 1L]
  padded[nr + 2L, 1] <- padded[nr + 1L, 1] + padded[nr + 2L, 2] - padded[nr + 1L, 2]
  padded[nr + 2L, nc + 2L] <- padded[nr + 1L, nc + 2L] + padded[nr + 2L, nc + 1L] - padded[nr + 1L, nc + 1L]

  (
    padded[1:(nr + 1L), 1:(nc + 1L)] +
      padded[2:(nr + 2L), 1:(nc + 1L)] +
      padded[1:(nr + 1L), 2:(nc + 2L)] +
      padded[2:(nr + 2L), 2:(nc + 2L)]
  ) / 4
}

ereefs_is_regular_centre_grid <- function(coord_table, tolerance = 1e-8) {
  lon_mat <- ereefs_matrix_from_coord_table(coord_table, "longitude")
  lat_mat <- ereefs_matrix_from_coord_table(coord_table, "latitude")

  lon_ref <- apply(lon_mat, 1, stats::median, na.rm = TRUE)
  lat_ref <- apply(lat_mat, 2, stats::median, na.rm = TRUE)

  lon_ok <- max(abs(lon_mat - lon_ref[row(lon_mat)]), na.rm = TRUE) < tolerance
  lat_ok <- max(abs(lat_mat - lat_ref[col(lat_mat)]), na.rm = TRUE) < tolerance
  lon_ok && lat_ok
}

ereefs_build_regular_corner_grids <- function(coord_table) {
  if (!ereefs_is_regular_centre_grid(coord_table)) {
    stop("This dataset only provides cell centres, but they do not lie on a regular axis-aligned grid.")
  }

  lon_mat <- ereefs_matrix_from_coord_table(coord_table, "longitude")
  lat_mat <- ereefs_matrix_from_coord_table(coord_table, "latitude")

  lon_centres <- apply(lon_mat, 1, stats::median, na.rm = TRUE)
  lat_centres <- apply(lat_mat, 2, stats::median, na.rm = TRUE)
  lon_edges <- ereefs_regular_edges(lon_centres)
  lat_edges <- ereefs_regular_edges(lat_centres)

  x_grid <- outer(lon_edges, lat_edges, function(x, y) x)
  y_grid <- outer(lon_edges, lat_edges, function(x, y) y)

  list(x_grid = x_grid, y_grid = y_grid)
}

ereefs_build_corner_grids_from_centres <- function(coord_table) {
  lon_mat <- ereefs_matrix_from_coord_table(coord_table, "longitude")
  lat_mat <- ereefs_matrix_from_coord_table(coord_table, "latitude")

  if (ereefs_is_regular_centre_grid(coord_table)) {
    return(ereefs_build_regular_corner_grids(coord_table))
  }

  list(
    x_grid = ereefs_centre_to_corner_matrix(lon_mat),
    y_grid = ereefs_centre_to_corner_matrix(lat_mat)
  )
}

ereefs_spatial_grid <- function(input_file) {
  coord_table <- ereefs_coord_table(input_file)
  coord_table %>%
    dplyr::arrange(i, j) %>%
    dplyr::select(i, j, latitude, longitude)
}

ereefs_is_catalog_request <- function(input_file) {
  if (!is.character(input_file) || length(input_file) != 1 || is.na(input_file)) {
    return(FALSE)
  }
  lowered <- tolower(input_file)
  lowered %in% c("catalog", "menu", "nci") ||
    grepl("catalog\\.(xml|html)$", lowered)
}

ereefs_prepare_input_reference <- function(input_file) {
  if (!is.character(input_file) || !length(input_file) || any(is.na(input_file))) {
    return(input_file)
  }
  if (length(input_file) > 1L) {
    return(vapply(input_file, ereefs_canonicalise_url, character(1)))
  }

  input_file <- ereefs_canonicalise_url(input_file)
  if (grepl("catalog\\.html$", input_file, ignore.case = TRUE)) {
    return(sub("catalog\\.html$", "catalog.xml", input_file, ignore.case = TRUE))
  }
  if (grepl("^https?://", input_file) && grepl("catalog\\.xml$", input_file, ignore.case = TRUE)) {
    return(input_file)
  }
  substitute_filename(input_file)
}

ereefs_grid_cache_key <- function(input_file) {
  input_file <- ereefs_prepare_input_reference(input_file)
  if (is.character(input_file) && length(input_file) == 1L && grepl("^https?://", input_file)) {
    return(paste0("grid_", gsub("[^A-Za-z0-9]+", "_", ereefs_canonicalise_url(input_file))))
  }
  if (is.character(input_file) && length(input_file) == 1L) {
    resolved <- tryCatch(normalizePath(input_file, winslash = "/", mustWork = FALSE), error = function(e) input_file)
    return(paste0("grid_", gsub("[^A-Za-z0-9]+", "_", resolved)))
  }
  paste0("grid_", gsub("[^A-Za-z0-9]+", "_", paste(input_file, collapse = "_")))
}

ereefs_grid_family_cache_key <- function(input_file) {
  input_file <- ereefs_prepare_input_reference(input_file)
  if (!is.character(input_file) || length(input_file) != 1L) {
    return(NA_character_)
  }

  normalized <- if (grepl("^https?://", input_file)) {
    ereefs_canonicalise_url(input_file)
  } else {
    tryCatch(normalizePath(input_file, winslash = "/", mustWork = FALSE), error = function(e) input_file)
  }

  family_ref <- sub("([_-])[0-9]{4}-[0-9]{2}-[0-9]{2}(\\.nc)$", "\\1DATE\\2", normalized)
  family_ref <- sub("([_-])[0-9]{4}-[0-9]{2}(\\.nc)$", "\\1DATE\\2", family_ref)
  if (identical(family_ref, normalized)) {
    return(NA_character_)
  }
  paste0("gridfam_", gsub("[^A-Za-z0-9]+", "_", family_ref))
}

ereefs_python_requests_get <- function(url) {
  ereefs_python_configure()
  requests <- reticulate::import("requests", delay_load = TRUE)
  response <- requests$get(url, timeout = 60)
  response$raise_for_status()
  reticulate::py_to_r(response$text)
}

ereefs_catalog_xml <- function(catalog_url) {
  cache_key <- paste0("catalog_", gsub("[^A-Za-z0-9]+", "_", catalog_url))
  if (exists(cache_key, envir = .ereefs_py_cache, inherits = FALSE)) {
    return(get(cache_key, envir = .ereefs_py_cache, inherits = FALSE))
  }

  xml_text <- ereefs_python_requests_get(catalog_url)
  xml <- xml2::read_xml(xml_text)
  assign(cache_key, xml, envir = .ereefs_py_cache)
  xml
}

ereefs_catalog_base_url <- function(catalog_url) {
  sub("catalog\\.(xml|html)$", "", ereefs_canonicalise_url(catalog_url))
}

ereefs_catalog_dods_base <- function(catalog_url) {
  base_url <- ereefs_catalog_base_url(catalog_url)
  sub("/catalog/", "/dodsC/", base_url, fixed = TRUE)
}

ereefs_extract_file_date <- function(dataset_name) {
  matched <- stringr::str_extract(dataset_name, "[0-9]{4}-[0-9]{2}(-[0-9]{2})?")
  if (is.na(matched)) {
    return(as.Date(NA))
  }
  if (nchar(matched) == 7L) {
    matched <- paste0(matched, "-01")
  }
  as.Date(matched)
}

ereefs_extract_file_date_precision <- function(dataset_name) {
  matched <- stringr::str_extract(dataset_name, "[0-9]{4}-[0-9]{2}(-[0-9]{2})?")
  if (is.na(matched)) {
    return(NA_character_)
  }
  if (nchar(matched) == 10L) {
    return("day")
  }
  if (nchar(matched) == 7L) {
    return("month")
  }
  NA_character_
}

ereefs_catalog_entries <- function(catalog_url, recurse = TRUE) {
  xml <- ereefs_catalog_xml(catalog_url)
  ns <- xml2::xml_ns(xml)

  datasets <- xml2::xml_find_all(xml, ".//d1:dataset[@urlPath]", ns = ns)
  direct_entries <- dplyr::tibble(
    name = xml2::xml_attr(datasets, "name"),
    url_path = xml2::xml_attr(datasets, "urlPath")
  ) %>%
    dplyr::filter(!is.na(url_path), stringr::str_ends(url_path, ".nc")) %>%
    dplyr::mutate(
      opendap_url = paste0(sub("/$", "", ereefs_catalog_dods_base(catalog_url)), "/", basename(url_path)),
      file_date = as.Date(unlist(lapply(name, ereefs_extract_file_date)), origin = "1970-01-01"),
      date_precision = unlist(lapply(name, ereefs_extract_file_date_precision))
    )

  if (!recurse) {
    return(direct_entries)
  }

  refs <- xml2::xml_find_all(xml, ".//d1:catalogRef", ns = ns)
  ref_hrefs <- xml2::xml_attr(refs, "href")
  ref_hrefs <- ref_hrefs[!is.na(ref_hrefs)]
  if (!length(ref_hrefs)) {
    return(direct_entries)
  }

  child_catalogs <- unique(xml2::url_absolute(ref_hrefs, ereefs_catalog_base_url(catalog_url)))
  child_entries <- lapply(child_catalogs, function(child_url) {
    tryCatch(
      ereefs_catalog_entries(child_url, recurse = FALSE),
      error = function(e) NULL
    )
  })
  child_entries <- child_entries[!vapply(child_entries, is.null, logical(1))]

  dplyr::bind_rows(c(list(direct_entries), child_entries)) %>%
    dplyr::distinct(opendap_url, .keep_all = TRUE)
}

ereefs_resolve_time_files <- function(input_file, start_date, end_date) {
  start_date <- as.Date(start_date)
  end_date <- as.Date(end_date)

  if (length(input_file) > 1L) {
    file_tbl <- dplyr::tibble(opendap_url = as.character(input_file)) %>%
      dplyr::mutate(file_date = as.Date(unlist(lapply(basename(opendap_url), ereefs_extract_file_date)), origin = "1970-01-01")) %>%
      dplyr::filter(!is.na(file_date), file_date >= start_date, file_date <= end_date)
    if (!nrow(file_tbl)) {
      stop("None of the supplied files overlap the requested timeframe.")
    }
    return(file_tbl %>% dplyr::select(file_date, opendap_url))
  }

  input_file <- ereefs_prepare_input_reference(input_file)

  if (!ereefs_is_catalog_request(input_file) && !grepl("\\.xml$", input_file, ignore.case = TRUE)) {
    return(dplyr::tibble(
      file_date = start_date,
      opendap_url = input_file,
      source = "explicit"
    ))
  }

  entries <- ereefs_catalog_entries(input_file, recurse = TRUE) %>%
    dplyr::filter(!is.na(file_date))

  if (!nrow(entries)) {
    stop("No dated NetCDF datasets were found in the supplied catalog.")
  }

  month_floor <- function(x) as.Date(format(as.Date(x), "%Y-%m-01"))
  requested_months <- seq.Date(month_floor(start_date), month_floor(end_date), by = "month")
  daily_available <- any(entries$date_precision == "day", na.rm = TRUE)
  monthly_available <- any(entries$date_precision == "month", na.rm = TRUE)

  selected <- if (daily_available) {
    entries %>%
      dplyr::filter(.data$file_date >= start_date, .data$file_date <= end_date)
  } else if (monthly_available) {
    entries %>%
      dplyr::filter(.data$file_date %in% requested_months)
  } else {
    entries %>%
      dplyr::filter(.data$file_date >= month_floor(start_date), .data$file_date <= end_date)
  }

  if (daily_available && nrow(selected)) {
    requested_days <- seq.Date(start_date, end_date, by = "day")
    missing_days <- setdiff(requested_days, selected$file_date)
    if (length(missing_days)) {
      warning(
        paste0(
          "The supplied catalog does not fully cover the requested daily timeframe. ",
          "Returning the available files between ",
          min(selected$file_date),
          " and ",
          max(selected$file_date),
          "."
        )
      )
    }
  }

  if (!nrow(selected)) {
    stop("The supplied catalog did not contain files covering the requested timeframe.")
  }

  selected %>%
    dplyr::arrange(file_date) %>%
    dplyr::select(file_date, opendap_url) %>%
    dplyr::distinct()
}

ereefs_pick_file_for_date <- function(file_table, target_date) {
  target_date <- as.Date(target_date)
  if (!nrow(file_table)) {
    stop("No files were available for the requested date.")
  }

  distances <- abs(file_table$file_date - target_date)
  file_table$opendap_url[[which.min(distances)]]
}

# metadata:
# - gpt_version: GPT-5 Codex
# - time: 21:26
# - date: 2026-04-26
# - prompt_used: "Continue the live OPeNDAP refactor by removing stem-based file reconstruction, resolving time ranges from THREDDS catalogs, and supporting centre-only NCI simple files."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 21:34
# - date: 2026-04-26
# - prompt_used: "Continue the live OPeNDAP refactor by removing stem-based file reconstruction, resolving time ranges from THREDDS catalogs, and supporting centre-only NCI simple files."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 21:37
# - date: 2026-04-26
# - prompt_used: "Continue the live OPeNDAP refactor by removing stem-based file reconstruction, resolving time ranges from THREDDS catalogs, and supporting centre-only NCI simple files."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 21:39
# - date: 2026-04-26
# - prompt_used: "Continue the live OPeNDAP refactor by removing stem-based file reconstruction, resolving time ranges from THREDDS catalogs, and supporting centre-only NCI simple files."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 21:43
# - date: 2026-04-26
# - prompt_used: "Continue the live OPeNDAP refactor by removing stem-based file reconstruction, resolving time ranges from THREDDS catalogs, and supporting centre-only NCI simple files."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 21:46
# - date: 2026-04-26
# - prompt_used: "Continue the live OPeNDAP refactor by removing stem-based file reconstruction, resolving time ranges from THREDDS catalogs, and supporting centre-only NCI simple files."

ereefs_time_var_name <- function(input_file) {
  var_names <- ereefs_var_names(input_file)
  for (candidate in ereefs_dim_aliases$time) {
    if (candidate %in% var_names) {
      return(candidate)
    }
  }
  stop("Could not identify the time variable in this dataset.")
}

ereefs_z_sign <- function(input_file, override_positive = FALSE) {
  positive <- ereefs_var_attr(input_file, "botz", "positive")
  sign_value <- if (identical(positive, "down")) -1 else 1
  if (override_positive) {
    sign_value <- -sign_value
  }
  sign_value
}

ereefs_var_dim_name <- function(input_file, var_name, role) {
  dims <- ereefs_var_dims(input_file, var_name)
  match_id <- match(role, dims$role)
  if (is.na(match_id)) {
    return(NA_character_)
  }
  dims$name[[match_id]]
}

ereefs_map_reference_var <- function(input_file, var_name) {
  if (var_name %in% c("true_color", "true_colour")) {
    return("R_645")
  }
  if (var_name == "plume") {
    return("R_412")
  }
  if (var_name == "ZooT") {
    return("ZooL_N")
  }
  if (var_name == "speed") {
    vars <- ereefs_var_names(input_file)
    if ("u1" %in% vars) {
      return("u1")
    }
    return("u")
  }
  var_name
}

ereefs_bbox_cells <- function(spatial_grid, box_bounds = c(NA, NA, NA, NA)) {
  if (all(is.na(box_bounds))) {
    return(spatial_grid %>% dplyr::arrange(i, j))
  }

  spatial_grid %>%
    dplyr::filter(
      dplyr::between(longitude, box_bounds[1], box_bounds[2]),
      dplyr::between(latitude, box_bounds[3], box_bounds[4])
    ) %>%
    dplyr::arrange(i, j)
}

ereefs_spatial_subset_spec <- function(input_file,
                                       var_name,
                                       spatial_grid,
                                       box_bounds = c(NA, NA, NA, NA),
                                       padding = 0L) {
  dims <- ereefs_var_dims(input_file, var_name)
  i_name <- dims$name[match("i", dims$role)]
  j_name <- dims$name[match("j", dims$role)]
  all_i <- sort(unique(spatial_grid$i))
  all_j <- sort(unique(spatial_grid$j))
  selected_cells <- ereefs_bbox_cells(spatial_grid, box_bounds)

  if (nrow(selected_cells) == 0) {
    stop("box_bounds does not intersect any model cells in this dataset.")
  }

  if (all(is.na(box_bounds))) {
    i_vals <- all_i
    j_vals <- all_j
  } else {
    i_range <- range(selected_cells$i)
    j_range <- range(selected_cells$j)
    i_min <- max(min(all_i), i_range[[1]] - padding)
    i_max <- min(max(all_i), i_range[[2]] + padding)
    j_min <- max(min(all_j), j_range[[1]] - padding)
    j_max <- min(max(all_j), j_range[[2]] + padding)
    i_vals <- seq.int(i_min, i_max)
    j_vals <- seq.int(j_min, j_max)
  }

  filters <- list()
  if (!is.na(i_name)) {
    filters[[i_name]] <- i_vals
  }
  if (!is.na(j_name)) {
    filters[[j_name]] <- j_vals
  }

  list(
    filters = filters,
    i_vals = i_vals,
    j_vals = j_vals,
    selected_cells = selected_cells
  )
}

ereefs_array_index_grid <- function(i_vals, j_vals) {
  dplyr::as_tibble(expand.grid(i = i_vals, j = j_vals))
}

ereefs_geocoordinates_tbl <- function(geocoordinates) {
  if (is.character(geocoordinates) && length(geocoordinates) == 1L && identical(geocoordinates, "mmp")) {
    return(dplyr::as_tibble(mmp_sites))
  }

  if (is.numeric(geocoordinates) && length(geocoordinates) == 2L && is.null(dim(geocoordinates))) {
    return(dplyr::tibble(latitude = geocoordinates[[1]], longitude = geocoordinates[[2]]))
  }

  tbl <- dplyr::as_tibble(geocoordinates)
  if (!all(c("latitude", "longitude") %in% names(tbl))) {
    stop("geocoordinates must provide latitude and longitude columns, a length-2 numeric vector, or 'mmp'.")
  }
  tbl
}

ereefs_match_geocoordinates <- function(geocoordinates, spatial_grid) {
  requested <- ereefs_geocoordinates_tbl(geocoordinates)
  idx <- vapply(
    seq_len(nrow(requested)),
    function(row_id) {
      which.min(
        earth.dist(
          requested$longitude[[row_id]],
          requested$latitude[[row_id]],
          spatial_grid$longitude,
          spatial_grid$latitude
        )
      )
    },
    integer(1)
  )

  dplyr::bind_cols(
    requested,
    spatial_grid[idx, , drop = FALSE] %>%
      dplyr::rename(cell_latitude = latitude, cell_longitude = longitude)
  )
}

ereefs_layer_geometry <- function(z_grid) {
  n_layers <- length(z_grid) - 1L
  if (n_layers < 1L) {
    stop("z_grid must contain at least two interfaces.")
  }

  layer_ids <- seq_len(n_layers)
  # EMS layers are ordered from deepest k = 1 to shallowest k = max(k).
  # z_grid is the matching sequence of layer interfaces from deepest to
  # shallowest, so each layer k spans z_grid[k] to z_grid[k + 1].
  bottom_static <- z_grid[-length(z_grid)]
  top_static <- z_grid[-1L]
  top_cap <- top_static
  top_cap[[n_layers]] <- Inf

  dplyr::tibble(
    k = layer_ids,
    top_static = top_static,
    bottom_static = bottom_static,
    top_cap = top_cap
  )
}

ereefs_reconstruct_z_grid_from_zc <- function(zc,
                                              top_padding = 1e20,
                                              input_file = NA_character_) {
  zc <- sort(as.numeric(zc))
  zc <- zc[is.finite(zc)]
  if (!length(zc)) {
    return(NULL)
  }
  if (length(zc) == 1L) {
    gap <- 1
    z_grid <- c(zc - gap, max(zc + gap, top_padding))
  } else {
    internal_interfaces <- (zc[-1L] + zc[-length(zc)]) / 2
    bottom_interface <- 2 * zc[[1]] - internal_interfaces[[1]]
    top_interface <- 2 * zc[[length(zc)]] - internal_interfaces[[length(internal_interfaces)]]
    top_interface <- max(top_interface, top_padding)
    z_grid <- c(bottom_interface, internal_interfaces, top_interface)
  }

  if (!all(diff(z_grid) > 0)) {
    return(NULL)
  }

  warn_key <- paste0(
    "reconstruct_z_grid_",
    if (!is.na(input_file) && nzchar(input_file)) {
      gsub("[^A-Za-z0-9]+", "_", input_file)
    } else {
      "generic"
    }
  )
  if (!exists(warn_key, envir = .ereefs_warning_cache, inherits = FALSE)) {
    warning(
      paste0(
        "z_grid was not available",
        if (!is.na(input_file)) paste0(" in ", input_file) else "",
        ". Reconstructing layer interfaces from zc, assuming interior interfaces lie midway between zc values and setting the upper interface to 1e20 to match the standard EMS convention."
      )
    )
    assign(warn_key, TRUE, envir = .ereefs_warning_cache)
  }
  z_grid
}

ereefs_layer_index_from_depth_msl <- function(depth, z_grid) {
  geometry <- ereefs_layer_geometry(z_grid)
  if (!is.numeric(depth) || length(depth) != 1L || is.na(depth)) {
    stop("Depth below mean sea level must be a single numeric value.")
  }
  target_level <- if (depth > 0) -depth else depth
  matches <- geometry$k[
    target_level <= geometry$top_static &
      target_level >= geometry$bottom_static
  ]
  if (length(matches)) {
    return(max(matches))
  }
  if (target_level > max(geometry$top_static)) {
    return(max(geometry$k))
  }
  min(geometry$k)
}

ereefs_time_indices <- function(ds,
                                start_date,
                                end_date,
                                nearest_if_empty = FALSE,
                                warn_context = "requested time") {
  idx <- which(ds >= start_date & ds <= end_date)
  if (length(idx) || !nearest_if_empty || !length(ds)) {
    return(idx)
  }

  target_time <- start_date + (end_date - start_date) / 2
  nearest_idx <- which.min(abs(as.numeric(ds - target_time)))
  warning(sprintf(
    "No model output matched %s (%s). Using nearest available time %s.",
    warn_context,
    format(target_time, tz = "UTC", usetz = TRUE),
    format(ds[[nearest_idx]], tz = "UTC", usetz = TRUE)
  ))
  nearest_idx
}

ereefs_time_filter_values <- function(input_file,
                                      time_index,
                                      raw_time,
                                      time_dim_name = "time",
                                      time_coord_dim = TRUE) {
  if (ereefs_is_remote_file(input_file) && ereefs_python_ready()) {
    return(time_index)
  }
  if (!isTRUE(time_coord_dim) || !identical(time_dim_name, "time")) {
    return(time_index)
  }
  raw_time[time_index]
}

ereefs_eta_file <- function(main_file, eta_stem, start_date, end_date) {
  vars <- ereefs_var_names(main_file)
  if ("eta" %in% vars) {
    return(list(file = main_file, assumed_zero = FALSE))
  }

  if (length(eta_stem) == 1L && !is.na(eta_stem)) {
    eta_ref <- ereefs_prepare_input_reference(eta_stem)
    eta_files <- ereefs_resolve_time_files(eta_ref, start_date, end_date)
    return(list(
      file = ereefs_pick_file_for_date(eta_files, start_date),
      assumed_zero = FALSE
    ))
  }

  list(file = NA_character_, assumed_zero = TRUE)
}

ereefs_eta_table <- function(main_file,
                             matched_points,
                             ds,
                             time_index,
                             start_date,
                             end_date,
                             eta_stem = NA) {
  eta_source <- ereefs_eta_file(main_file, eta_stem, start_date, end_date)
  base_tbl <- expand.grid(
    row_id = matched_points$row_id,
    time = ds[time_index]
  ) %>%
    dplyr::as_tibble() %>%
    dplyr::left_join(
      matched_points %>% dplyr::select(row_id, i, j),
      by = "row_id"
    )

  if (isTRUE(eta_source$assumed_zero)) {
    warning("eta is not available in this dataset, so depth-below-free-surface calculations are assuming eta = 0.")
    return(base_tbl %>% dplyr::mutate(eta = 0))
  }

  eta_file <- eta_source$file
  eta_timing <- get_origin_and_times(eta_file)
  eta_index <- ereefs_time_indices(
    eta_timing[[2]],
    start_date,
    end_date,
    nearest_if_empty = isTRUE(all.equal(start_date, end_date)),
    warn_context = "requested eta time"
  )
  if (!length(eta_index)) {
    return(base_tbl %>% dplyr::mutate(eta = 0))
  }

  eta_dims <- ereefs_var_dims(eta_file, "eta")
  i_vals <- sort(unique(matched_points$i))
  j_vals <- sort(unique(matched_points$j))
  filters <- list()
  filters[[eta_dims$name[match("i", eta_dims$role)]]] <- i_vals
  filters[[eta_dims$name[match("j", eta_dims$role)]]] <- j_vals
  filters[[eta_dims$name[match("time", eta_dims$role)]]] <- ereefs_time_filter_values(
    eta_file,
    eta_index,
    eta_timing[[3]],
    time_dim_name = eta_dims$name[match("time", eta_dims$role)],
    time_coord_dim = eta_dims$coord_dim[match("time", eta_dims$role)]
  )
  eta_arr <- ereefs_read_var_array(eta_file, "eta", filters = filters)
  eta_arr <- ereefs_role_array(
    eta_arr,
    dims = eta_dims,
    target_roles = c("i", "j", "time"),
    drop_singleton = FALSE
  )
  eta_arr <- array(eta_arr, dim = c(length(i_vals), length(j_vals), length(eta_index)))
  eta_tbl <- dplyr::as_tibble(expand.grid(
    i = i_vals,
    j = j_vals,
    time = eta_timing[[2]][eta_index]
  )) %>%
    dplyr::mutate(eta = as.vector(eta_arr))

  if (!nrow(eta_tbl)) {
    return(base_tbl %>% dplyr::mutate(eta = 0))
  }

  base_tbl %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      eta = {
        matches <- which(eta_tbl$i == i & eta_tbl$j == j)
        if (!length(matches)) {
          0
        } else {
          time_matches <- matches[which.min(abs(as.numeric(eta_tbl$time[matches] - time)))]
          eta_tbl$eta[[time_matches]]
        }
      }
    ) %>%
    dplyr::ungroup()
}

ereefs_wet_layer_table <- function(botz_tbl, eta_tbl, z_grid) {
  geometry <- ereefs_layer_geometry(z_grid)

  expand.grid(
    row_id = unique(botz_tbl$row_id),
    time = unique(eta_tbl$time),
    k = geometry$k
  ) %>%
    dplyr::as_tibble() %>%
    dplyr::left_join(botz_tbl, by = "row_id") %>%
    dplyr::left_join(eta_tbl, by = c("row_id", "i", "j", "time")) %>%
    dplyr::left_join(geometry, by = "k") %>%
    dplyr::mutate(
      wet_top = pmin(top_cap, eta),
      wet_bottom = pmax(bottom_static, botz),
      dz = pmax(wet_top - wet_bottom, 0),
      watercol_depth = pmax(eta - botz, 0)
    )
}

ereefs_extract_target_k <- function(wet_tbl, depth) {
  wet_tbl %>%
    dplyr::mutate(
      depth_top = pmax(0, eta - wet_top),
      depth_bottom = pmax(0, eta - wet_bottom)
    ) %>%
    dplyr::group_by(row_id, time) %>%
    dplyr::summarise(
      k = {
        wet_rows <- dplyr::cur_data_all() %>% dplyr::filter(dz > 0)
        if (!nrow(wet_rows)) {
          NA_integer_
        } else {
          matches <- wet_rows$k[
            depth >= wet_rows$depth_top & depth <= wet_rows$depth_bottom
          ]
          if (length(matches)) {
            max(matches)
          } else if (depth > max(wet_rows$depth_bottom)) {
            min(wet_rows$k)
          } else {
            max(wet_rows$k)
          }
        }
      },
      .groups = "drop"
    )
}

ereefs_surface_bottom_k <- function(wet_tbl, which_layer = c("surface", "bottom")) {
  which_layer <- match.arg(which_layer)
  wet_tbl %>%
    dplyr::filter(dz > 0) %>%
    dplyr::group_by(row_id, time) %>%
    dplyr::summarise(
      k = if (which_layer == "surface") max(k) else min(k),
      .groups = "drop"
    )
}

ereefs_surface_target_k <- function(wet_tbl, var_tbl = NULL, value_name = NULL) {
  wet_surface <- ereefs_surface_bottom_k(wet_tbl, which_layer = "surface")
  if (is.null(var_tbl) || is.null(value_name)) {
    return(wet_surface)
  }

  candidate_surface <- var_tbl %>%
    dplyr::filter(is.finite(.data[[value_name]])) %>%
    dplyr::inner_join(
      wet_surface %>% dplyr::rename(k_wet = k),
      by = c("row_id", "time")
    ) %>%
    dplyr::filter(.data$k >= .data$k_wet) %>%
    dplyr::group_by(.data$row_id, .data$time) %>%
    dplyr::summarise(k = max(.data$k), .groups = "drop")

  wet_surface %>%
    dplyr::rename(k_wet = k) %>%
    dplyr::left_join(candidate_surface, by = c("row_id", "time")) %>%
    dplyr::transmute(
      row_id = .data$row_id,
      time = .data$time,
      k = dplyr::coalesce(.data$k, .data$k_wet)
    )
}

ereefs_promote_surface_fill_vector <- function(values_by_k, wet_k) {
  if (!is.numeric(wet_k) || length(wet_k) != 1L || is.na(wet_k)) {
    return(values_by_k)
  }

  k_index <- seq_along(values_by_k)
  source_k <- k_index[is.finite(values_by_k) & k_index >= wet_k]
  if (!length(source_k)) {
    return(values_by_k)
  }

  source_k <- max(source_k)
  values_by_k[[wet_k]] <- values_by_k[[source_k]]
  if (wet_k < length(values_by_k)) {
    values_by_k[(wet_k + 1L):length(values_by_k)] <- NA_real_
  }
  values_by_k
}

ereefs_resolve_input_file <- function(input_stem, ereefs_case, year, month, day = 1L) {
  if (ereefs_case[2] == "4km") {
    return(paste0(input_stem, format(as.Date(sprintf("%04d-%02d-01", year, month)), "%Y-%m"), ".nc"))
  }
  if (ereefs_case[2] == "1km") {
    return(paste0(input_stem, format(as.Date(sprintf("%04d-%02d-%02d", year, month, day)), "%Y-%m-%d"), ".nc"))
  }
  paste0(input_stem, ".nc")
}

ereefs_dates_for_stride <- function(start_date, end_date, stride) {
  if (identical(stride, "daily")) {
    return(seq(as.Date(start_date), as.Date(end_date), by = "1 day"))
  }

  stride_num <- suppressWarnings(as.numeric(stride))
  if (is.na(stride_num) || stride_num <= 0) {
    stop("stride must be 'daily' or a positive number of days.")
  }
  seq(as.Date(start_date), as.Date(end_date), by = stride_num)
}

ereefs_towns <- function() {
  dplyr::tibble(
    latitude = c(-15.47027987, -16.0899, -16.4840, -16.92303816, -19.26639219, -20.0136699, -20.07670986,
                 -20.40109791, -21.15345122, -22.82406858, -23.38031858, -23.84761069, -24.8662122, -25.54073075,
                 -26.18916037),
    longitude = c(145.2498605, 145.4622, 145.4623, 145.7710, 146.805701, 148.2475387, 146.2635394, 148.5802016,
                  149.1655418, 147.6363616, 150.5059485, 151.256349, 152.3478987, 152.7049316, 152.6581893),
    town = c("Cooktown", "Cape Tribulation", "Port Douglas", "Cairns", "Townsville", "Bowen", "Charters Towers",
             "Prosperine", "Mackay", "Clermont", "Rockhampton", "Gladstone", "Bundaberg", "Maryborough", "Gympie")
  )
}

# metadata:
# - gpt_version: GPT-5 Codex
# - time: 15:37
# - date: 2026-04-26
# - prompt_used: "Refactor the eReefs R toolkit away from ncdf4 toward tidync, add regular-grid support alongside curvilinear grids, archive existing R files, and refresh docs/examples."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 19:05
# - date: 2026-04-26
# - prompt_used: "Install dependencies, verify the refactored toolkit, improve efficiency for large THREDDS-served files, and build a working Jupyter demo notebook."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 19:49
# - date: 2026-04-26
# - prompt_used: "Install dependencies, verify the refactored toolkit, improve efficiency for large THREDDS-served files, and build a working Jupyter demo notebook."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 22:41
# - date: 2026-04-26
# - prompt_used: "Finish the tidy tidync-first refactor, keep it efficient for large live OPeNDAP datasets, validate depth/free-surface logic, audit dependencies, and review plotting palettes."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 00:04
# - date: 2026-04-27
# - prompt_used: "Finish the tidy tidync-first refactor, keep it efficient for large live OPeNDAP datasets, validate depth/free-surface logic, audit dependencies, and review plotting palettes."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 01:28
# - date: 2026-04-27
# - prompt_used: "Fix local EMS time filtering so coordinate-based time dimensions use raw values while record-style time dimensions use positional indices."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 01:40
# - date: 2026-04-27
# - prompt_used: "Add a zc-based z_grid reconstruction fallback for simple EMS files, with an explicit warning about the assumptions."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 14:09
# - date: 2026-04-27
# - prompt_used: "After comparing the reconstructed simple-file z_grid to the standard-file z_grid, reset the top interface to 1e20 for macrotidal safety."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 14:24
# - date: 2026-04-27
# - prompt_used: "Fix map/grid orientation bugs by reordering arrays explicitly by i/j/k/time roles before flattening them onto polygons."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 14:30
# - date: 2026-04-27
# - prompt_used: "Make the role-based array reordering filter-aware so singleton time and k dimensions dropped by the backend do not break map orientation fixes."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 14:33
# - date: 2026-04-27
# - prompt_used: "Finish the orientation fix by applying a full permutation before dropping singleton dimensions."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 14:39
# - date: 2026-04-27
# - prompt_used: "Mask NetCDF-style sentinel values like 1e35 before plotting so notebook maps do not show false bright edge bands."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 14:47
# - date: 2026-04-27
# - prompt_used: "Fix regular-grid corner construction so centre-only regular products build proper rectangular cells instead of vertical stripe slivers."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 15:43
# - date: 2026-04-27
# - prompt_used: "Single-time requests should fall back to the nearest available model time with a warning if there is no exact timestamp match."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 18:15
# - date: 2026-04-27
# - prompt_used: "Fix the live two-month time-series example, keep movie colour scales fixed across frames, and suppress repeated z_grid reconstruction warnings."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 12:05
# - date: 2026-04-28
# - prompt_used: "Tighten the shared vertical geometry logic so EMS k layers map to the correct z_grid interfaces and live simple-file slice plots do not place the surface layer above sea level."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 13:14
# - date: 2026-04-28
# - prompt_used: "Generalize the copied-up surface-layer convention so shared wet-layer helpers, surface extraction, and vertical arrays all treat simple-file surface values consistently."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 13:35
# - date: 2026-04-28
# - prompt_used: "Make z_grid and grid metadata reuse more efficient by caching get_ereefs_grids() results so repeated layer/profile calls do not recalculate grids for the same source file."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 13:43
# - date: 2026-04-28
# - prompt_used: "Extend grid caching across dated files in the same catalog family so z_grid and spatial grids are reused across daily/monthly siblings as well as exact file repeats."
