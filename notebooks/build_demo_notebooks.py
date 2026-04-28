import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
NOTEBOOKS_DIR = REPO_ROOT / "notebooks"


def lines(text: str):
    text = text.strip("\n")
    return [line + "\n" for line in text.splitlines()]


def markdown_cell(text: str, cell_id: str):
    return {
        "cell_type": "markdown",
        "id": cell_id,
        "metadata": {},
        "source": lines(text),
    }


def code_cell(text: str, cell_id: str):
    return {
        "cell_type": "code",
        "id": cell_id,
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": lines(text),
    }


def notebook_metadata(previous_metadata, prompt_used: str, time_text: str):
    previous_entry = previous_metadata.get("codex_metadata")
    history = list(previous_metadata.get("codex_metadata_history", []))
    if previous_entry is not None:
        history.append(previous_entry)
    return {
        "kernelspec": {
            "display_name": "R 4.3 (ereefs local)",
            "language": "R",
            "name": "ir43-ereefs",
        },
        "language_info": {
            "name": "R",
        },
        "codex_metadata": {
            "gpt_version": "GPT-5 Codex",
            "time": time_text,
            "date": "2026-04-27",
            "prompt_used": prompt_used,
        },
        "codex_metadata_history": history,
    }


def build_local_demo(previous_metadata):
    prompt_used = "Fix the smooth map plotting path and expand the demo notebooks so they cover the important documented ereefs workflows."
    cell_index = 0

    def cid(prefix: str):
        nonlocal cell_index
        cell_index += 1
        return f"{prefix}-{cell_index}"

    cells = [
        markdown_cell(
            """
            # eReefs demo notebook

            This notebook demonstrates the refactored `ereefs` toolkit on small synthetic EMS-style NetCDF files. It covers the main documented local workflows:

            - regular and curvilinear grid discovery
            - standard and smoother display maps
            - surface, bottom, and depth-specific time series
            - depth-averaged and depth-below-free-surface extraction
            - vertical profiles
            - transect slices where the requested line does not exactly pass through cell centres
            """,
            cid("md"),
        ),
        code_cell(
            """
            repo_root <- if (file.exists("DESCRIPTION")) getwd() else normalizePath("..", winslash = "/", mustWork = TRUE)
            .libPaths(c(file.path(repo_root, ".r_libs"), .libPaths()))

            suppressPackageStartupMessages({
              library(pkgload)
              library(dplyr)
              library(ggplot2)
              library(tibble)
            })

            pkgload::load_all(repo_root)
            source(file.path(repo_root, "notebooks", "create_demo_datasets.R"))
            demo_paths <- create_demo_datasets(file.path(repo_root, "notebooks", "demo_data"))
            demo_paths
            """,
            cid("code"),
        ),
        code_cell(
            """
            notebook_output_dir <- file.path(repo_root, "notebooks", "output", "local_demo")
            dir.create(notebook_output_dir, recursive = TRUE, showWarnings = FALSE)

            save_plot_display <- function(plot_obj, filename, width = 8, height = 6, dpi = 120) {
              out_path <- file.path(notebook_output_dir, filename)
              ggplot2::ggsave(out_path, plot = plot_obj, width = width, height = height, dpi = dpi)
              print(plot_obj)
              out_path
            }
            """,
            cid("code"),
        ),
        markdown_cell(
            """
            This block inspects the horizontal grid metadata for both demo files so we can confirm that the package is distinguishing regular and curvilinear geometries correctly.
            """,
            cid("md"),
        ),
        code_cell(
            """
            regular_grids <- get_ereefs_grids(demo_paths$regular)
            curvilinear_grids <- get_ereefs_grids(demo_paths$curvilinear)

            list(
              regular = list(
                grid_type = regular_grids$grid_type,
                x_grid_dim = dim(regular_grids$x_grid),
                y_grid_dim = dim(regular_grids$y_grid),
                n_cells = nrow(regular_grids$spatial_grid)
              ),
              curvilinear = list(
                grid_type = curvilinear_grids$grid_type,
                x_grid_dim = dim(curvilinear_grids$x_grid),
                y_grid_dim = dim(curvilinear_grids$y_grid),
                n_cells = nrow(curvilinear_grids$spatial_grid)
              )
            )
            """,
            cid("code"),
        ),
        markdown_cell(
            """
            These examples draw the same regular-grid field twice: first as polygons and then with the optional smoother rasterised display mode.
            """,
            cid("md"),
        ),
        code_cell(
            """
            box <- c(146.2, 147.3, -20.4, -19.2)

            regular_map_polygon <- map_ereefs(
              var_name = "temp",
              target_date = as.Date("2020-01-02"),
              layer = "surface",
              input_file = demo_paths$regular,
              box_bounds = box,
              scale_col = "viridis",
              suppress_print = TRUE,
              return_poly = TRUE,
              label_towns = FALSE
            )

            regular_map_plot <- plot_map(
              regular_map_polygon,
              box_bounds = box,
              scale_col = "viridis",
              suppress_print = TRUE,
              label_towns = FALSE
            )
            save_plot_display(regular_map_plot, "regular_map_polygon.png")

            regular_map_smooth <- plot_map(
              regular_map_polygon,
              box_bounds = box,
              scale_col = "viridis",
              plot_style = "smooth",
              smooth_pixels = 400,
              suppress_print = TRUE,
              label_towns = FALSE,
              var_longname = "Same data, smoother display"
            )
            save_plot_display(regular_map_smooth, "regular_map_smooth.png")
            """,
            cid("code"),
        ),
        markdown_cell(
            """
            Here we extract time series from the surface, the deepest wet layer, and a fixed depth below mean sea level at one location.
            """,
            cid("md"),
        ),
        code_cell(
            """
            ts_surface <- get_ereefs_ts(
              var_names = c("temp", "salt"),
              geocoordinates = tibble(latitude = -19.75, longitude = 146.75),
              layer = "surface",
              start_date = as.POSIXct("2020-01-01 00:00:00", tz = "Etc/GMT-10"),
              end_date = as.POSIXct("2020-01-03 00:00:00", tz = "Etc/GMT-10"),
              input_file = demo_paths$regular,
              verbosity = 0
            )

            ts_bottom <- get_ereefs_ts(
              var_names = "temp",
              geocoordinates = tibble(latitude = -19.75, longitude = 146.75),
              layer = "bottom",
              start_date = as.POSIXct("2020-01-01 00:00:00", tz = "Etc/GMT-10"),
              end_date = as.POSIXct("2020-01-03 00:00:00", tz = "Etc/GMT-10"),
              input_file = demo_paths$regular,
              verbosity = 0
            )

            ts_msl <- get_ereefs_ts(
              var_names = "temp",
              geocoordinates = tibble(latitude = -19.75, longitude = 146.75),
              layer = -5,
              start_date = as.POSIXct("2020-01-01 00:00:00", tz = "Etc/GMT-10"),
              end_date = as.POSIXct("2020-01-03 00:00:00", tz = "Etc/GMT-10"),
              input_file = demo_paths$regular,
              verbosity = 0
            )

            list(
              surface_head = head(ts_surface),
              bottom_head = head(ts_bottom),
              depth_below_msl_head = head(ts_msl)
            )
            """,
            cid("code"),
        ),
        markdown_cell(
            """
            This block demonstrates depth-averaged extraction and the fallback behaviour when a simple-format file does not provide `eta`.
            """,
            cid("md"),
        ),
        code_cell(
            """
            depth_avg <- get_ereefs_depth_integrated_ts(
              var_names = "temp",
              geocoordinates = tibble(latitude = -19.75, longitude = 146.75),
              start_date = as.POSIXct("2020-01-01 00:00:00", tz = "Etc/GMT-10"),
              end_date = as.POSIXct("2020-01-03 00:00:00", tz = "Etc/GMT-10"),
              input_file = demo_paths$regular,
              verbosity = 0
            )

            eta_warning <- NULL
            depth_free_surface <- withCallingHandlers(
              get_ereefs_depth_specified_ts(
                var_names = "temp",
                geocoordinates = tibble(latitude = -19.75, longitude = 146.75),
                depth = 2,
                start_date = as.POSIXct("2020-01-01 00:00:00", tz = "Etc/GMT-10"),
                end_date = as.POSIXct("2020-01-03 00:00:00", tz = "Etc/GMT-10"),
                input_file = demo_paths$regular_noeta,
                verbosity = 0
              ),
              warning = function(w) {
                eta_warning <<- conditionMessage(w)
                invokeRestart("muffleWarning")
              }
            )

            list(
              depth_average_head = head(depth_avg),
              depth_below_free_surface_head = head(depth_free_surface),
              no_eta_warning = eta_warning
            )
            """,
            cid("code"),
        ),
        markdown_cell(
            """
            The next example extracts a vertical profile from the curvilinear demo file and plots it at the requested time.
            """,
            cid("md"),
        ),
        code_cell(
            """
            profile <- get_ereefs_profile(
              var_names = "temp",
              geolocation = c(-19.4, 147.75),
              start_date = as.POSIXct("2020-01-01 00:00:00", tz = "Etc/GMT-10"),
              end_date = as.POSIXct("2020-01-02 00:00:00", tz = "Etc/GMT-10"),
              input_file = demo_paths$curvilinear
            )

            profile_plot <- plot_ereefs_profile(profile, var_name = "temp", target_date = as.Date("2020-01-01"))
            save_plot_display(profile_plot, "curvilinear_profile.png")
            """,
            cid("code"),
        ),
        markdown_cell(
            """
            Finally, this demo extracts a vertical slice along a transect that does not align exactly with cell centres and plots the result.
            """,
            cid("md"),
        ),
        code_cell(
            """
            slice <- get_ereefs_slice(
              var_names = "temp",
              geolocation = data.frame(
                latitude = c(-20.05, -18.7),
                longitude = c(147.05, 148.35)
              ),
              target_date = as.POSIXct("2020-01-01 00:00:00", tz = "Etc/GMT-10"),
              input_file = demo_paths$curvilinear
            )

            slice_plot <- plot_ereefs_slice(slice, var_name = "temp", scale_col = "viridis")
            save_plot_display(slice_plot, "curvilinear_slice.png", width = 9, height = 5)
            """,
            cid("code"),
        ),
    ]
    return {
        "cells": cells,
        "metadata": notebook_metadata(previous_metadata, prompt_used, "00:47"),
        "nbformat": 4,
        "nbformat_minor": 5,
    }


def build_live_demo(previous_metadata):
    prompt_used = "Fix the smooth map plotting path and expand the demo notebooks so they cover the important documented ereefs workflows."
    cell_index = 0

    def cid(prefix: str):
        nonlocal cell_index
        cell_index += 1
        return f"{prefix}-{cell_index}"

    cells = [
        markdown_cell(
            """
            # eReefs live OPeNDAP demo

            This notebook demonstrates remote access without downloading the full NetCDF files locally. It focuses on the main live workflows:

            - NCI `simple` file access through OPeNDAP
            - grid discovery for centre-only products
            - standard and smoother map rendering
            - time-aware catalog resolution across AIMS monthly files
            - surface time series across multiple live monthly files
            - vertical profiles, vertical slices, and surface transects from live data
            - depth-aware extraction when `eta` is absent
            - a short whole-of-GBR `map_ereefs_movie()` example
            """,
            cid("md"),
        ),
        code_cell(
            """
            repo_root <- normalizePath(if (basename(getwd()) == "notebooks") ".." else ".", winslash = "/")
            setwd(repo_root)
            .libPaths(c(normalizePath(".r_libs", winslash = "/", mustWork = FALSE), .libPaths()))

            suppressPackageStartupMessages({
              library(pkgload)
              library(dplyr)
              library(ggplot2)
              library(tibble)
            })

            pkgload::load_all(".")

            aims_catalog <- "https://thredds.ereefs.aims.gov.au/thredds/catalog/ereefs/gbr1_2.0/stats-monthly-monthly/catalog.xml"
            nci_catalog <- "https://thredds.nci.org.au/thredds/catalog/fx3/gbr4_H4p0_ABARRAr2_OBRAN2020_FG2Gv3_B4p2_Cq5b_Dhnd/catalog.xml"
            nci_simple <- "https://thredds.nci.org.au/thredds/dodsC/fx3/gbr4_H4p0_ABARRAr2_OBRAN2020_FG2Gv3_B4p2_Cq5b_Dhnd/gbr4_H4p0_ABARRAr2_OBRAN2020_FG2G_B4p2_Cq5b_Dhnd_simple_2022-10-30.nc"
            """,
            cid("code"),
        ),
        code_cell(
            """
            notebook_output_dir <- file.path(repo_root, "notebooks", "output", "live_demo")
            dir.create(notebook_output_dir, recursive = TRUE, showWarnings = FALSE)

            save_plot_display <- function(plot_obj, filename, width = 8, height = 6, dpi = 120) {
              out_path <- file.path(notebook_output_dir, filename)
              ggplot2::ggsave(out_path, plot = plot_obj, width = width, height = height, dpi = dpi)
              print(plot_obj)
              out_path
            }

            display_animation_file <- function(path) {
              if (is.na(path) || !file.exists(path)) {
                return(invisible(path))
              }
              if (grepl("\\\\.gif$", path, ignore.case = TRUE) && requireNamespace("magick", quietly = TRUE)) {
                print(magick::image_read(path))
              } else {
                cat(path, "\\n")
              }
              invisible(path)
            }
            """,
            cid("code"),
        ),
        markdown_cell(
            """
            This block inspects the NCI simple file grid metadata, including the reconstructed plotting geometry used for centre-based live products.
            """,
            cid("md"),
        ),
        code_cell(
            """
            nci_grids <- get_ereefs_grids(nci_simple)
            tibble(
              grid_type = nci_grids$grid_type,
              n_cells = nrow(nci_grids$spatial_grid),
              x_grid_rows = nrow(nci_grids$x_grid),
              x_grid_cols = ncol(nci_grids$x_grid)
            )
            """,
            cid("code"),
        ),
        markdown_cell(
            """
            These two maps show the same NCI Secchi field with polygon rendering and the optional smoother display mode.
            """,
            cid("md"),
        ),
        code_cell(
            """
            nci_box <- c(145, 155, -25, -10)

            nci_map <- map_ereefs(
              var_name = "Secchi",
              target_date = as.Date("2022-10-30"),
              input_file = nci_simple,
              box_bounds = nci_box,
              scale_col = "viridis",
              suppress_print = TRUE,
              return_poly = TRUE,
              label_towns = FALSE
            )

            nci_map_plot <- plot_map(
              nci_map,
              box_bounds = nci_box,
              scale_col = "viridis",
              suppress_print = TRUE,
              label_towns = FALSE
            )
            save_plot_display(nci_map_plot, "nci_secchi_polygon.png")

            nci_map_smooth <- plot_map(
              nci_map,
              box_bounds = nci_box,
              scale_col = "viridis",
              plot_style = "smooth",
              smooth_pixels = 500,
              suppress_print = TRUE,
              label_towns = FALSE,
              var_longname = "NCI Secchi depth with smoother display"
            )
            save_plot_display(nci_map_smooth, "nci_secchi_smooth.png")
            """,
            cid("code"),
        ),
        markdown_cell(
            """
            Here we compare a direct surface extraction from the NCI simple file with a depth-below-free-surface extraction that has to assume `eta = 0`.
            """,
            cid("md"),
        ),
        code_cell(
            """
            live_surface <- get_ereefs_ts(
              var_names = "Secchi",
              geocoordinates = tibble(latitude = -19.5, longitude = 148.0),
              start_date = as.POSIXct("2022-10-30 00:00:00", tz = "Etc/GMT-10"),
              end_date = as.POSIXct("2022-10-30 23:59:59", tz = "Etc/GMT-10"),
              input_file = nci_simple,
              verbosity = 0
            )

            eta_warning <- NULL
            live_depth <- withCallingHandlers(
              get_ereefs_depth_specified_ts(
                var_names = "NH4",
                geocoordinates = tibble(latitude = -19.5, longitude = 148.0),
                depth = 2,
                start_date = as.POSIXct("2022-10-30 00:00:00", tz = "Etc/GMT-10"),
                end_date = as.POSIXct("2022-10-30 23:59:59", tz = "Etc/GMT-10"),
                input_file = nci_simple,
                verbosity = 0
              ),
              warning = function(w) {
                eta_warning <<- conditionMessage(w)
                invokeRestart("muffleWarning")
              }
            )

            list(
              surface_head = head(live_surface),
              depth_head = head(live_depth),
              no_eta_warning = eta_warning
            )
            """,
            cid("code"),
        ),
        markdown_cell(
            """
            This example resolves a two-month time series across the live NCI catalog rather than a single file and then plots the extracted surface ammonium series.
            """,
            cid("md"),
        ),
        code_cell(
            """
            nci_surface_warning <- NULL
            nci_surface_ts <- withCallingHandlers(
              get_ereefs_ts(
                var_names = "NH4",
                geocoordinates = tibble(latitude = -19.5, longitude = 148.0),
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

            attr(nci_surface_ts, "warning_text") <- nci_surface_warning

            nci_surface_ts_plot <- ggplot(nci_surface_ts, aes(x = time, y = NH4)) +
              geom_line(colour = "#00798c", linewidth = 0.8) +
              geom_point(colour = "#d1495b", size = 0.9) +
              labs(
                title = "Two-month live surface time series",
                x = NULL,
                y = "Surface NH4 (mg N m-3)"
              ) +
              theme_minimal(base_size = 12)
            save_plot_display(nci_surface_ts_plot, "nci_surface_timeseries.png", width = 9, height = 5)
            """,
            cid("code"),
        ),
        markdown_cell(
            """
            This block asks the AIMS monthly catalog which files are needed to span the requested period.
            """,
            cid("md"),
        ),
        code_cell(
            """
            aims_files <- ereefs_resolve_time_files(
              aims_catalog,
              as.Date("2019-10-01"),
              as.Date("2019-12-31")
            )
            aims_files
            """,
            cid("code"),
        ),
        markdown_cell(
            """
            We then map one AIMS monthly field to show the regular-grid live workflow over a large spatial domain.
            """,
            cid("md"),
        ),
        code_cell(
            """
            aims_box <- c(142.0, 155.4, -28.7, -7.3)
            aims_map <- map_ereefs(
              var_name = "temp_mean",
              target_date = as.Date("2019-10-01"),
              layer = "surface",
              input_file = aims_catalog,
              box_bounds = aims_box,
              scale_col = "viridis",
              suppress_print = TRUE,
              return_poly = TRUE,
              label_towns = FALSE
            )

            aims_map_plot <- plot_map(
              aims_map,
              box_bounds = aims_box,
              scale_col = "viridis",
              plot_style = "smooth",
              smooth_pixels = 500,
              suppress_print = TRUE,
              label_towns = FALSE
            )
            save_plot_display(aims_map_plot, "aims_surface_map_smooth.png")
            """,
            cid("code"),
        ),
        markdown_cell(
            """
            The next block extracts a vertical profile from the live NCI simple file at one geographic location and plots the depth structure.
            """,
            cid("md"),
        ),
        code_cell(
            """
            live_profile <- get_ereefs_profile(
              var_names = "NH4",
              geolocation = c(-19.5, 148.0),
              start_date = as.POSIXct("2022-10-30 12:00:00", tz = "Etc/GMT-10"),
              end_date = as.POSIXct("2022-10-30 12:00:00", tz = "Etc/GMT-10"),
              input_file = nci_simple
            )

            live_profile_plot <- plot_ereefs_profile(
              live_profile,
              var_name = "NH4",
              target_date = as.Date("2022-10-30")
            )
            save_plot_display(live_profile_plot, "nci_vertical_profile.png")
            """,
            cid("code"),
        ),
        markdown_cell(
            """
            This block extracts a vertical slice along a short transect and plots the resulting section.
            """,
            cid("md"),
        ),
        code_cell(
            """
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

            live_slice_plot <- plot_ereefs_slice(live_slice, var_name = "NH4", scale_col = "viridis", var_units = "mg N m-3")
            save_plot_display(live_slice_plot, "nci_vertical_slice.png", width = 9, height = 5)
            """,
            cid("code"),
        ),
        markdown_cell(
            """
            Here we densify the transect, match it to nearby cells, and plot the resulting surface transect as a one-dimensional series along distance.
            """,
            cid("md"),
        ),
        code_cell(
            """
            transect_cells <- ereefs_densify_path(transect_line, samples_per_segment = 60) %>%
              ereefs_nearest_cells(get_ereefs_grids(nci_simple)$spatial_grid) %>%
              dplyr::distinct(i, j, .keep_all = TRUE)

            surface_transect <- get_ereefs_ts(
              var_names = "NH4",
              geocoordinates = transect_cells %>% dplyr::select(latitude, longitude),
              layer = "surface",
              start_date = as.POSIXct("2022-10-30 12:00:00", tz = "Etc/GMT-10"),
              end_date = as.POSIXct("2022-10-30 12:00:00", tz = "Etc/GMT-10"),
              input_file = nci_simple,
              verbosity = 0
            ) %>%
              dplyr::mutate(
                distance_km = c(
                  0,
                  cumsum(vapply(
                    seq_len(dplyr::n() - 1),
                    function(idx) earth.dist(
                      longitude[[idx]],
                      latitude[[idx]],
                      longitude[[idx + 1]],
                      latitude[[idx + 1]]
                    ),
                    numeric(1)
                  ))
                )
              )

            surface_transect_plot <- ggplot(surface_transect, aes(x = distance_km, y = NH4)) +
              geom_line(colour = "#2a6f97", linewidth = 0.8) +
              geom_point(colour = "#f4a261", size = 1.7) +
              labs(
                title = "Surface transect from live OPeNDAP data",
                x = "Distance along transect (km)",
                y = "Surface NH4 (mg N m-3)"
              ) +
              theme_minimal(base_size = 12)
            save_plot_display(surface_transect_plot, "nci_surface_transect.png", width = 9, height = 5)
            """,
            cid("code"),
        ),
        markdown_cell(
            """
            The closing example builds a short whole-of-GBR movie from AIMS monthly files, saves each frame, and shows the summary output plus animation artifact.
            """,
            cid("md"),
        ),
        code_cell(
            """
            movie_output_dir <- file.path(notebook_output_dir, "movie_frames")
            movie_demo <- map_ereefs_movie(
              var_name = "temp_mean",
              start_date = as.Date("2019-10-01"),
              end_date = as.Date("2019-11-01"),
              layer = "surface",
              output_dir = movie_output_dir,
              save_frames = TRUE,
              animation_format = "gif",
              input_file = aims_catalog,
              box_bounds = aims_box,
              plot_style = "smooth",
              smooth_pixels = 220,
              stride = 31,
              suppress_print = TRUE,
              fps = 1,
              label_towns = FALSE
            )

            movie_summary_plot <- plot_map(
              movie_demo,
              box_bounds = aims_box,
              scale_col = "viridis",
              plot_style = "smooth",
              smooth_pixels = 220,
              suppress_print = TRUE,
              label_towns = FALSE,
              var_longname = "Average surface temperature across the selected frames",
              var_units = "degC"
            )
            save_plot_display(movie_summary_plot, "aims_movie_summary.png")
            display_animation_file(movie_demo$animation_file)
            list(
              animation_file = movie_demo$animation_file,
              n_frames = length(movie_demo$frame_files),
              frame_dir = movie_demo$output_dir
            )
            """,
            cid("code"),
        ),
    ]
    return {
        "cells": cells,
        "metadata": notebook_metadata(previous_metadata, prompt_used, "00:47"),
        "nbformat": 4,
        "nbformat_minor": 5,
    }


def main():
    notebook_specs = {
        "ereefs_demo.ipynb": build_local_demo,
        "ereefs_live_opendap_demo.ipynb": build_live_demo,
    }
    for filename, builder in notebook_specs.items():
        path = NOTEBOOKS_DIR / filename
        existing = json.loads(path.read_text(encoding="utf-8")) if path.exists() else {"metadata": {}}
        notebook = builder(existing.get("metadata", {}))
        path.write_text(json.dumps(notebook, indent=2) + "\n", encoding="utf-8")
        print(f"wrote {path}")


if __name__ == "__main__":
    main()

# metadata:
# - gpt_version: GPT-5 Codex
# - time: 00:47
# - date: 2026-04-27
# - prompt_used: "Fix the smooth map plotting path and expand the demo notebooks so they cover the important documented ereefs workflows."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 01:00
# - date: 2026-04-27
# - prompt_used: "Tighten the demo notebook builder by fixing invalid example variables and adding stable notebook cell ids for execution."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 01:08
# - date: 2026-04-27
# - prompt_used: "Fix the live OPeNDAP notebook examples to use variables that actually exist in the selected NCI simple file."
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
# - time: 15:52
# - date: 2026-04-27
# - prompt_used: "The default slice path should use a bulk subset read like the archived implementation instead of building slices profile by profile."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 16:12
# - date: 2026-04-27
# - prompt_used: "The default slice path should use a bulk subset read like the archived implementation instead of building slices profile by profile."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 16:45
# - date: 2026-04-27
# - prompt_used: "The default slice path should use a bulk subset read like the archived implementation instead of building slices profile by profile."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 16:58
# - date: 2026-04-27
# - prompt_used: "Make notebook plots save as well as display, and have map_ereefs_movie save frames then build an animation file."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 17:09
# - date: 2026-04-27
# - prompt_used: "Make notebook plots save as well as display, and have map_ereefs_movie save frames then build an animation file."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 18:15
# - date: 2026-04-27
# - prompt_used: "Fix the live two-month time-series example, keep movie colour scales fixed across frames, and suppress repeated z_grid reconstruction warnings."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 09:03
# - date: 2026-04-28
# - prompt_used: "Fix the vignette profile and time-series examples, avoid clipped first-map labels, and add explanatory markdown blocks throughout the notebooks."
