.libPaths(c(normalizePath(file.path("..", ".r_libs"), winslash = "/", mustWork = FALSE), .libPaths()))

create_demo_datasets <- function(output_dir = "demo_data") {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  regular_file <- file.path(output_dir, "regular_demo_2020-01.nc")
  regular_noeta_file <- file.path(output_dir, "regular_demo_noeta_2020-01.nc")
  curvilinear_file <- file.path(output_dir, "curvilinear_demo_2020-01.nc")

  make_regular_demo <- function(path, include_eta = TRUE) {
    i_vals <- c(146.25, 146.75, 147.25, 147.75)
    j_vals <- c(-20.25, -19.75, -19.25)
    k_vals <- c(-2, -8)
    z_grid <- c(0, -5, -12)
    time_vals <- c(0, 1, 2)

    dim_i <- ncdf4::ncdim_def("i", "", seq_along(i_vals), create_dimvar = TRUE)
    dim_j <- ncdf4::ncdim_def("j", "", seq_along(j_vals), create_dimvar = TRUE)
    dim_k <- ncdf4::ncdim_def("k", "", seq_along(k_vals), create_dimvar = TRUE)
    dim_k_grid <- ncdf4::ncdim_def("k_grid", "", seq_along(z_grid), create_dimvar = TRUE)
    dim_time <- ncdf4::ncdim_def("time", "days since 2020-01-01 00:00:00", time_vals, unlim = FALSE)

    var_defs <- list(
      ncdf4::ncvar_def("longitude", "degrees_east", list(dim_i), missval = NA_real_, prec = "double"),
      ncdf4::ncvar_def("latitude", "degrees_north", list(dim_j), missval = NA_real_, prec = "double"),
      ncdf4::ncvar_def("z_grid", "m", list(dim_k_grid), missval = NA_real_, prec = "double"),
      ncdf4::ncvar_def("botz", "m", list(dim_i, dim_j), missval = NA_real_, prec = "double"),
      ncdf4::ncvar_def("temp", "degC", list(dim_i, dim_j, dim_k, dim_time), missval = NA_real_, prec = "float"),
      ncdf4::ncvar_def("salt", "PSU", list(dim_i, dim_j, dim_k, dim_time), missval = NA_real_, prec = "float"),
      ncdf4::ncvar_def("u", "m s-1", list(dim_i, dim_j, dim_k, dim_time), missval = NA_real_, prec = "float"),
      ncdf4::ncvar_def("v", "m s-1", list(dim_i, dim_j, dim_k, dim_time), missval = NA_real_, prec = "float")
    )
    if (include_eta) {
      var_defs <- append(
        var_defs,
        list(ncdf4::ncvar_def("eta", "m", list(dim_i, dim_j, dim_time), missval = NA_real_, prec = "float")),
        after = 4
      )
    }

    nc <- ncdf4::nc_create(path, var_defs, force_v4 = TRUE)
    on.exit(ncdf4::nc_close(nc), add = TRUE)

    botz <- outer(seq_along(i_vals), seq_along(j_vals), function(i, j) -(8 + 0.7 * i + 0.3 * j))
    eta <- array(0, dim = c(length(i_vals), length(j_vals), length(time_vals)))
    temp <- array(0, dim = c(length(i_vals), length(j_vals), length(k_vals), length(time_vals)))
    salt <- array(0, dim = dim(temp))
    u <- array(0, dim = dim(temp))
    v <- array(0, dim = dim(temp))

    for (ti in seq_along(time_vals)) {
      for (ki in seq_along(k_vals)) {
        for (ii in seq_along(i_vals)) {
          for (jj in seq_along(j_vals)) {
            eta[ii, jj, ti] <- 0.05 * ti + 0.02 * ii - 0.01 * jj
            temp[ii, jj, ki, ti] <- 22 + 0.4 * ii - 0.3 * jj - 0.8 * (ki - 1) + 0.5 * ti
            salt[ii, jj, ki, ti] <- 34.5 + 0.1 * jj + 0.05 * ki - 0.03 * ti
            u[ii, jj, ki, ti] <- 0.08 * ii + 0.02 * ti
            v[ii, jj, ki, ti] <- -0.05 * jj + 0.01 * ti
          }
        }
      }
    }

    ncdf4::ncvar_put(nc, "longitude", i_vals)
    ncdf4::ncvar_put(nc, "latitude", j_vals)
    ncdf4::ncvar_put(nc, "z_grid", z_grid)
    ncdf4::ncvar_put(nc, "botz", botz)
    if (include_eta) ncdf4::ncvar_put(nc, "eta", eta)
    ncdf4::ncvar_put(nc, "temp", temp)
    ncdf4::ncvar_put(nc, "salt", salt)
    ncdf4::ncvar_put(nc, "u", u)
    ncdf4::ncvar_put(nc, "v", v)
    ncdf4::ncatt_put(nc, "botz", "positive", "up")
    invisible(path)
  }

  make_curvilinear_demo <- function(path) {
    ni <- 3
    nj <- 4
    nk <- 2
    time_vals <- c(0, 1)
    z_grid <- c(0, -4, -11)

    lon_centres <- matrix(rep(c(147.1, 147.7, 148.3), each = nj), nrow = ni, ncol = nj)
    lat_centres <- matrix(rep(c(-20.1, -19.6, -19.1, -18.6), times = ni), nrow = ni, ncol = nj)
    lon_centres <- lon_centres + matrix(c(0, 0.03, 0.02, 0, 0.02, 0.04, 0.03, 0.01, 0.03, 0.05, 0.04, 0.02), nrow = ni)
    lat_centres <- lat_centres + matrix(c(0, 0.02, 0.03, 0.01, 0.01, 0.03, 0.04, 0.02, 0.02, 0.04, 0.05, 0.03), nrow = ni)

    lon_edges <- matrix(NA_real_, nrow = ni + 1, ncol = nj + 1)
    lat_edges <- matrix(NA_real_, nrow = ni + 1, ncol = nj + 1)
    for (ii in seq_len(ni + 1)) {
      for (jj in seq_len(nj + 1)) {
        lon_edges[ii, jj] <- 146.8 + 0.6 * (ii - 1) + 0.03 * (jj - 1)
        lat_edges[ii, jj] <- -20.4 + 0.5 * (jj - 1) + 0.02 * (ii - 1)
      }
    }

    dim_i <- ncdf4::ncdim_def("i", "", seq_len(ni), create_dimvar = TRUE)
    dim_j <- ncdf4::ncdim_def("j", "", seq_len(nj), create_dimvar = TRUE)
    dim_i_edge <- ncdf4::ncdim_def("i_edge", "", seq_len(ni + 1), create_dimvar = TRUE)
    dim_j_edge <- ncdf4::ncdim_def("j_edge", "", seq_len(nj + 1), create_dimvar = TRUE)
    dim_k <- ncdf4::ncdim_def("k", "", seq_len(nk), create_dimvar = TRUE)
    dim_k_grid <- ncdf4::ncdim_def("k_grid", "", seq_along(z_grid), create_dimvar = TRUE)
    dim_time <- ncdf4::ncdim_def("time", "days since 2020-01-01 00:00:00", time_vals, unlim = FALSE)

    var_defs <- list(
      ncdf4::ncvar_def("longitude", "degrees_east", list(dim_i, dim_j), missval = NA_real_, prec = "double"),
      ncdf4::ncvar_def("latitude", "degrees_north", list(dim_i, dim_j), missval = NA_real_, prec = "double"),
      ncdf4::ncvar_def("x_grid", "degrees_east", list(dim_i_edge, dim_j_edge), missval = NA_real_, prec = "double"),
      ncdf4::ncvar_def("y_grid", "degrees_north", list(dim_i_edge, dim_j_edge), missval = NA_real_, prec = "double"),
      ncdf4::ncvar_def("z_grid", "m", list(dim_k_grid), missval = NA_real_, prec = "double"),
      ncdf4::ncvar_def("botz", "m", list(dim_i, dim_j), missval = NA_real_, prec = "double"),
      ncdf4::ncvar_def("eta", "m", list(dim_i, dim_j, dim_time), missval = NA_real_, prec = "float"),
      ncdf4::ncvar_def("temp", "degC", list(dim_i, dim_j, dim_k, dim_time), missval = NA_real_, prec = "float"),
      ncdf4::ncvar_def("u", "m s-1", list(dim_i, dim_j, dim_k, dim_time), missval = NA_real_, prec = "float"),
      ncdf4::ncvar_def("v", "m s-1", list(dim_i, dim_j, dim_k, dim_time), missval = NA_real_, prec = "float")
    )

    nc <- ncdf4::nc_create(path, var_defs, force_v4 = TRUE)
    on.exit(ncdf4::nc_close(nc), add = TRUE)

    botz <- outer(seq_len(ni), seq_len(nj), function(i, j) -(6 + i + 0.5 * j))
    eta <- array(0, dim = c(ni, nj, length(time_vals)))
    temp <- array(0, dim = c(ni, nj, nk, length(time_vals)))
    u <- array(0, dim = dim(temp))
    v <- array(0, dim = dim(temp))

    for (ti in seq_along(time_vals)) {
      for (ki in seq_len(nk)) {
        for (ii in seq_len(ni)) {
          for (jj in seq_len(nj)) {
            eta[ii, jj, ti] <- 0.03 * ti + 0.01 * ii + 0.02 * jj
            temp[ii, jj, ki, ti] <- 23 + 0.6 * ii - 0.4 * jj - 0.9 * (ki - 1) + 0.3 * ti
            u[ii, jj, ki, ti] <- 0.04 * ii + 0.02 * ki
            v[ii, jj, ki, ti] <- 0.03 * jj - 0.01 * ti
          }
        }
      }
    }

    ncdf4::ncvar_put(nc, "longitude", lon_centres)
    ncdf4::ncvar_put(nc, "latitude", lat_centres)
    ncdf4::ncvar_put(nc, "x_grid", lon_edges)
    ncdf4::ncvar_put(nc, "y_grid", lat_edges)
    ncdf4::ncvar_put(nc, "z_grid", z_grid)
    ncdf4::ncvar_put(nc, "botz", botz)
    ncdf4::ncvar_put(nc, "eta", eta)
    ncdf4::ncvar_put(nc, "temp", temp)
    ncdf4::ncvar_put(nc, "u", u)
    ncdf4::ncvar_put(nc, "v", v)
    ncdf4::ncatt_put(nc, "botz", "positive", "up")
    invisible(path)
  }

  make_regular_demo(regular_file, include_eta = TRUE)
  make_regular_demo(regular_noeta_file, include_eta = FALSE)
  make_curvilinear_demo(curvilinear_file)

  list(
    regular = normalizePath(regular_file, winslash = "/", mustWork = TRUE),
    regular_noeta = normalizePath(regular_noeta_file, winslash = "/", mustWork = TRUE),
    curvilinear = normalizePath(curvilinear_file, winslash = "/", mustWork = TRUE)
  )
}

# metadata:
# - gpt_version: GPT-5 Codex
# - time: 19:17
# - date: 2026-04-26
# - prompt_used: "Install dependencies, verify the refactored toolkit, improve efficiency for large THREDDS-served files, and build a working Jupyter demo notebook."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 23:43
# - date: 2026-04-26
# - prompt_used: "Finish the tidy tidync-first refactor, keep it efficient for large live OPeNDAP datasets, validate depth/free-surface logic, audit dependencies, and review plotting palettes."
