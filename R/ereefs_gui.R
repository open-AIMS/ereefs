#' Launch the experimental eReefs GUI
#'
#' Opens a Shiny application for interactive eReefs data exploration. The GUI is
#' designed as a staged workflow: inspect a local NetCDF file, OPeNDAP URL, or
#' THREDDS catalog; preview a map; optionally create an animation; and extract
#' point time series or profiles with explicit slow-request confirmation.
#'
#' @param input_file Optional initial NetCDF file path, OPeNDAP URL, THREDDS
#'   catalog URI, or legacy shortcut to pre-fill in the GUI.
#' @param launch.browser Passed to [shiny::runApp()]. Defaults to `TRUE`.
#' @param port,host Optional Shiny port and host values passed to
#'   [shiny::runApp()]. Leave as `NULL` to use Shiny defaults.
#' @param ... Additional arguments passed to [shiny::runApp()].
#'
#' @return Invisibly returns the result of [shiny::runApp()].
#' @export
#'
#' @examples
#' \dontrun{
#' launch_ereefs_gui()
#' launch_ereefs_gui("notebooks/demo_data/regular_demo_2020-01.nc")
#' }
launch_ereefs_gui <- function(input_file = "",
                              launch.browser = TRUE,
                              port = NULL,
                              host = NULL,
                              ...) {
  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop("The experimental GUI requires the optional shiny package. Install it with install.packages('shiny').", call. = FALSE)
  }

  app <- ereefs_gui_app(input_file = input_file)
  run_args <- list(appDir = app, launch.browser = launch.browser)
  if (!is.null(port)) {
    run_args$port <- port
  }
  if (!is.null(host)) {
    run_args$host <- host
  }
  run_args <- c(run_args, list(...))
  invisible(do.call(shiny::runApp, run_args))
}

ereefs_gui_app <- function(input_file = "") {
  ereefs_gui_add_resource_path()
  ui <- shiny::fluidPage(
    title = "eReefs R GUI",
    shiny::tags$head(shiny::tags$style(ereefs_gui_css())),
    shiny::tags$section(
      class = "ereefs-hero",
      shiny::tags$div(
        class = "ereefs-hero-copy",
        shiny::tags$img(
          class = "ereefs-logo",
          src = "ereefs-gui/ereefs_logo.png",
          alt = "eReefs logo"
        ),
        shiny::tags$span(class = "eyebrow", "Experimental R GUI"),
        shiny::tags$h1("eReefs R GUI"),
        shiny::tags$p(
          "A guided, scriptable workspace for inspecting eReefs datasets, previewing maps, ",
          "making animations, and extracting point data without accidentally pulling huge remote arrays."
        )
      ),
      shiny::tags$div(
        class = "ereefs-hero-badges",
        shiny::tags$span("Metadata first"),
        shiny::tags$span("OPeNDAP aware"),
        shiny::tags$span("R-script export")
      )
    ),
    shiny::tags$section(
      class = "source-ribbon",
      shiny::tags$strong("Current source"),
      shiny::verbatimTextOutput("current_source", placeholder = TRUE)
    ),
    shiny::tabsetPanel(
      id = "workflow_stage",
      type = "tabs",
      shiny::tabPanel("1. Source", ereefs_gui_source_tab(input_file)),
      shiny::tabPanel("2. Map", ereefs_gui_map_tab()),
      shiny::tabPanel("3. Animation", ereefs_gui_animation_tab()),
      shiny::tabPanel("4. Points", ereefs_gui_points_tab()),
      shiny::tabPanel("5. Transect", value = "5. Transect", ereefs_gui_transect_tab()),
      shiny::tabPanel("Notes", value = "Notes", ereefs_gui_notes_tab())
    )
  )

  server <- function(input, output, session) {
    state <- shiny::reactiveValues(
      nci_catalogs = NULL,
      inspection = NULL,
      inspected_source = NULL,
      has_multiple_layers = TRUE,
      transect_tab_visible = TRUE,
      map_result = NULL,
      map_fetch_signature = NULL,
      map_script = NULL,
      map_status = "No map preview has been requested yet.",
      animation_result = NULL,
      animation_script = NULL,
      animation_progress = "No animation has been requested yet.",
      animation_frame_file = NULL,
      animation_frame_src = NULL,
      animation_media_src = NULL,
      animation_media_type = NULL,
      points_result = NULL,
      points_script = NULL,
      profile_result = NULL,
      profile_script = NULL,
      slice_result = NULL,
      slice_script = NULL,
      local_file_path = "",
      local_browser_dir = ereefs_gui_default_local_browser_dir()
    )
    suppress_nci_warning <- shiny::reactiveVal(FALSE)
    syncing_map_scale_inputs <- shiny::reactiveVal(FALSE)
    syncing_bbox_inputs <- shiny::reactiveVal(FALSE)
    syncing_variable_inputs <- shiny::reactiveVal(FALSE)

    output$has_multiple_layers <- shiny::renderText({
      tolower(as.character(isTRUE(state$has_multiple_layers)))
    })
    shiny::outputOptions(output, "has_multiple_layers", suspendWhenHidden = FALSE)

    output$map_layer_hint <- shiny::renderText({
      ereefs_gui_layer_hint(state$inspection, input$map_variable)
    })
    shiny::outputOptions(output, "map_layer_hint", suspendWhenHidden = FALSE)

    output$points_layer_hint <- shiny::renderText({
      ereefs_gui_layer_hint(state$inspection, input$points_variable)
    })
    shiny::outputOptions(output, "points_layer_hint", suspendWhenHidden = FALSE)

    output$points_start_range <- shiny::renderText({
      ereefs_gui_date_range_label(state$inspection)
    })
    shiny::outputOptions(output, "points_start_range", suspendWhenHidden = FALSE)

    output$points_end_range <- shiny::renderText({
      ereefs_gui_date_range_label(state$inspection)
    })
    shiny::outputOptions(output, "points_end_range", suspendWhenHidden = FALSE)

    output$slice_date_range <- shiny::renderText({
      ereefs_gui_date_range_label(state$inspection)
    })
    shiny::outputOptions(output, "slice_date_range", suspendWhenHidden = FALSE)

    observeEvent(input$map_variable, {
      if (isTRUE(syncing_variable_inputs())) {
        return(NULL)
      }
      ereefs_gui_sync_variable_inputs(
        session = session,
        inspection = state$inspection,
        selected_var = input$map_variable,
        syncing_flag = syncing_variable_inputs
      )
    }, ignoreInit = TRUE)

    observeEvent(input$points_variable, {
      if (isTRUE(syncing_variable_inputs())) {
        return(NULL)
      }
      ereefs_gui_sync_variable_inputs(
        session = session,
        inspection = state$inspection,
        selected_var = input$points_variable,
        syncing_flag = syncing_variable_inputs
      )
    }, ignoreInit = TRUE)

    observeEvent(input$slice_variable, {
      if (isTRUE(syncing_variable_inputs())) {
        return(NULL)
      }
      ereefs_gui_sync_variable_inputs(
        session = session,
        inspection = state$inspection,
        selected_var = input$slice_variable,
        syncing_flag = syncing_variable_inputs
      )
    }, ignoreInit = TRUE)

    observeEvent(input$points_start, {
      ereefs_gui_clamp_date_input(session, "points_start", input$points_start, state$inspection)
    }, ignoreInit = TRUE)

    observeEvent(input$points_end, {
      ereefs_gui_clamp_date_input(session, "points_end", input$points_end, state$inspection)
    }, ignoreInit = TRUE)

    observeEvent(input$slice_date, {
      ereefs_gui_clamp_date_input(session, "slice_date", input$slice_date, state$inspection)
    }, ignoreInit = TRUE)

    observe({
      should_show <- isTRUE(state$has_multiple_layers)
      if (should_show && !isTRUE(state$transect_tab_visible)) {
        shiny::insertTab(
          inputId = "workflow_stage",
          target = "Notes",
          position = "before",
          tab = shiny::tabPanel("5. Transect", value = "5. Transect", ereefs_gui_transect_tab()),
          select = FALSE,
          session = session
        )
        state$transect_tab_visible <- TRUE
      } else if (!should_show && isTRUE(state$transect_tab_visible)) {
        if (identical(input$workflow_stage %||% "", "5. Transect")) {
          shiny::updateTabsetPanel(session, "workflow_stage", selected = "4. Points")
        }
        shiny::removeTab(
          inputId = "workflow_stage",
          target = "5. Transect",
          session = session
        )
        state$transect_tab_visible <- FALSE
      }
    })

    observeEvent(input$nci_help, {
      shiny::showModal(shiny::modalDialog(
        title = "NCI eReefs catalog naming",
        easyClose = TRUE,
        footer = shiny::modalButton("Close"),
        shiny::p("NCI eReefs dataset IDs are compact descriptions of the model configuration. Read them left to right as a set of configuration tokens separated by underscores."),
        shiny::tags$ul(
          shiny::tags$li(shiny::strong("Grid: "), "`GBR4` is the regional 4 km Great Barrier Reef model; `GBR1` is the higher-resolution 1 km shelf model."),
          shiny::tags$li(shiny::strong("Hydrodynamics: "), "`H4p0` means hydrodynamic model version 4.0; older datasets often use tokens such as `H2p0`."),
          shiny::tags$li(shiny::strong("Forcing: "), "Tokens such as `ABARRAr2`, `OBRAN2020`, and `FG2Gv3` describe atmospheric, ocean-boundary, river-flow, and other forcing-data selections."),
          shiny::tags$li(shiny::strong("Biogeochemistry, tracers, or scenarios: "), "Optional tokens such as `B4p2`, `Cq5b`, `Rivers`, or pesticide/scenario codes identify extra transport, BGC, sediment, optical, catchment, or tracer configurations."),
          shiny::tags$li(shiny::strong("Run mode: "), "`Dhnd` indicates a hindcast-style dataset. Older near-real-time products may use NRT-style labels and are often superseded by current hindcasts."),
          shiny::tags$li(shiny::strong("Files: "), "Catalogs usually contain many NetCDF files, with each file representing one day or one month and including the simulated date in its filename.")
        ),
        shiny::p("Example: `GBR4_H4p0_ABARRAr2_OBRAN2020_FG2Gv3_B4p2_Cq5b_Dhnd` is a GBR4 v4.0 hydrodynamics-driven biogeochemistry/sediment hindcast with specified forcing and catchment-configuration tokens."),
        shiny::p("When in doubt, prefer current `GBR4_H4p0` datasets for new GBR4 work and use the catalog date range to confirm that the selected source covers your requested period.")
      ))
    })

    observeEvent(input$source_mode, {
      if (!identical(input$source_mode, "nci") || isTRUE(suppress_nci_warning())) {
        return(NULL)
      }
      shiny::showModal(shiny::modalDialog(
        title = "NCI catalog requests can be slow",
        easyClose = FALSE,
        shiny::p(
          "NCI THREDDS/OPeNDAP datasets are served remotely and can contain thousands of large NetCDF files. ",
          "Catalog inspection is metadata-first, but map previews, animations, profiles, and extracts may still take a while."
        ),
        shiny::p("Use smaller bounding boxes, shorter date ranges, and a single variable while exploring."),
        footer = shiny::tagList(
          shiny::actionButton("nci_warning_cancel", "Cancel", class = "btn-secondary"),
          shiny::actionButton("nci_warning_ok", "OK", class = "btn-primary"),
          shiny::actionButton("nci_warning_dont_show", "Don't show again this session", class = "btn-outline-secondary")
        )
      ))
    }, ignoreInit = TRUE)

    observeEvent(input$nci_warning_cancel, {
      shiny::removeModal()
      shiny::updateRadioButtons(session, "source_mode", selected = "manual")
    }, ignoreInit = TRUE)

    observeEvent(input$nci_warning_ok, {
      shiny::removeModal()
    }, ignoreInit = TRUE)

    observeEvent(input$nci_warning_dont_show, {
      suppress_nci_warning(TRUE)
      shiny::removeModal()
    }, ignoreInit = TRUE)

    observeEvent(input$source_mode, {
      if (!identical(input$source_mode, "nci")) {
        return(NULL)
      }
      shiny::updateSelectInput(
        session,
        "nci_catalog",
        choices = stats::setNames("", "Loading NCI catalog menu...")
      )
      state$nci_catalogs <- tryCatch(
        {
          shiny::withProgress(message = "Reading the NCI catalog menu", value = 0.5, {
            ereefs_gui_nci_catalogs()
          })
        },
        error = function(e) {
          shiny::showNotification(paste("Could not read the NCI catalog menu:", conditionMessage(e)), type = "error")
          dplyr::tibble(title = character(), url = character())
        }
      )
      shiny::updateSelectInput(
        session,
        "nci_catalog",
        choices = stats::setNames(state$nci_catalogs$url, state$nci_catalogs$title)
      )
    }, ignoreInit = TRUE)

    observeEvent(input$local_browser_home, {
      state$local_browser_dir <- ereefs_gui_default_local_browser_dir()
      shiny::updateTextInput(session, "local_browser_dir_input", value = state$local_browser_dir)
    }, ignoreInit = TRUE)

    observeEvent(input$local_browser_up, {
      parent <- dirname(state$local_browser_dir)
      if (dir.exists(parent)) {
        state$local_browser_dir <- normalizePath(parent, winslash = "\\", mustWork = FALSE)
        shiny::updateTextInput(session, "local_browser_dir_input", value = state$local_browser_dir)
      }
    }, ignoreInit = TRUE)

    observeEvent(input$local_browser_go, {
      requested <- ereefs_gui_clean_source_path(input$local_browser_dir_input)
      if (dir.exists(requested)) {
        state$local_browser_dir <- normalizePath(requested, winslash = "\\", mustWork = FALSE)
      } else {
        shiny::showNotification("That folder could not be found.", type = "warning")
      }
      shiny::updateTextInput(session, "local_browser_dir_input", value = state$local_browser_dir)
    }, ignoreInit = TRUE)

    observeEvent(input$local_browser_open, {
      selected <- input$local_browser_entry %||% ""
      if (!nzchar(selected)) {
        return(NULL)
      }
      if (dir.exists(selected)) {
        state$local_browser_dir <- normalizePath(selected, winslash = "\\", mustWork = FALSE)
        shiny::updateTextInput(session, "local_browser_dir_input", value = state$local_browser_dir)
      } else if (file.exists(selected)) {
        state$local_file_path <- normalizePath(selected, winslash = "\\", mustWork = FALSE)
        shiny::updateTextInput(session, "local_file_path_display", value = state$local_file_path)
      }
    }, ignoreInit = TRUE)

    observeEvent(input$local_browser_select, {
      selected <- input$local_browser_entry %||% ""
      if (nzchar(selected) && file.exists(selected) && !dir.exists(selected)) {
        state$local_file_path <- normalizePath(selected, winslash = "\\", mustWork = FALSE)
        shiny::updateTextInput(session, "local_file_path_display", value = state$local_file_path)
      } else {
        shiny::showNotification("Select a NetCDF or mnc file first.", type = "warning")
      }
    }, ignoreInit = TRUE)

    output$local_browser_entries <- shiny::renderUI({
      entries <- ereefs_gui_local_browser_entries(state$local_browser_dir)
      if (!nrow(entries)) {
        return(shiny::p(class = "help-text", "No readable folders or NetCDF/mnc files were found in this directory."))
      }
      shiny::selectInput(
        "local_browser_entry",
        "Folders and NetCDF/mnc files",
        choices = stats::setNames(entries$path, entries$label),
        selected = entries$path[[1]],
        size = min(14L, max(4L, nrow(entries))),
        selectize = FALSE
      )
    })

    selected_source <- shiny::reactive({
      source_mode <- input$source_mode %||% "manual"
      if (identical(source_mode, "browse") && nzchar(state$local_file_path %||% "")) {
        return(state$local_file_path)
      }
      if (identical(source_mode, "nci")) {
        shiny::req(nzchar(input$nci_catalog %||% ""))
        return(input$nci_catalog)
      }
      ereefs_gui_clean_source_path(input$manual_source %||% input_file)
    })

    output$current_source <- shiny::renderText({
      source <- selected_source()
      if (!nzchar(source %||% "")) {
        return("No NetCDF file, OPeNDAP URL, or THREDDS catalog selected yet.")
      }
      paste(
        paste("Mode:", input$source_mode %||% "manual"),
        paste("Address:", source),
        paste("Inspection:", if (isTRUE(input$read_remote_metadata)) "catalog plus full remote metadata" else "catalog plus variable list from first remote file"),
        sep = "\n"
      )
    })

    observeEvent(selected_source(), {
      source <- selected_source()
      if (identical(source, state$inspected_source)) {
        return(NULL)
      }
      state$inspection <- NULL
      state$inspected_source <- NULL
      state$has_multiple_layers <- TRUE
      state$map_result <- NULL
      state$map_fetch_signature <- NULL
      state$map_script <- NULL
      state$animation_result <- NULL
      state$animation_frame_file <- NULL
      state$animation_frame_src <- NULL
      state$animation_script <- NULL
      state$points_result <- NULL
      state$points_script <- NULL
      state$profile_result <- NULL
      state$profile_script <- NULL
      state$slice_result <- NULL
      state$slice_script <- NULL
      ereefs_gui_reset_controls(session)
    }, ignoreInit = TRUE)

    observeEvent(input$inspect_btn, {
      shiny::req(nzchar(selected_source()))
      state$inspection <- tryCatch(
        {
          shiny::withProgress(message = "Inspecting metadata", value = 0.15, {
            progress_callback <- function(info) {
              shiny::setProgress(
                value = ereefs_gui_inspection_progress_fraction(info),
                detail = ereefs_gui_progress_detail(info)
              )
              ereefs_gui_flush(session)
            }
            result <- inspect_ereefs_data(
              input_file = selected_source(),
              recurse_catalog = TRUE,
              include_files = TRUE,
              max_files = input$max_files,
              read_remote_metadata = isTRUE(input$read_remote_metadata),
              progress_callback = progress_callback
            )
            result
          })
        },
        error = function(e) {
          shiny::showNotification(paste("Inspection failed:", conditionMessage(e)), type = "error", duration = 10)
          NULL
        }
      )
      if (!is.null(state$inspection)) {
        state$inspected_source <- selected_source()
        state$has_multiple_layers <- ereefs_gui_has_multiple_layers(state$inspection)
        state$map_result <- NULL
        state$map_fetch_signature <- NULL
        state$map_script <- NULL
        state$animation_result <- NULL
        state$animation_frame_file <- NULL
        state$animation_frame_src <- NULL
        state$animation_script <- NULL
        state$points_result <- NULL
        state$points_script <- NULL
        state$profile_result <- NULL
        state$profile_script <- NULL
        state$slice_result <- NULL
        state$slice_script <- NULL
        ereefs_gui_update_from_inspection(session, state$inspection)
        shiny::updateTabsetPanel(session, "workflow_stage", selected = "2. Map")
      }
    }, ignoreInit = TRUE)

    output$source_status <- shiny::renderText({
      if (is.null(state$inspection)) {
        return("Choose a source and click Inspect dataset. If you change the source, inspect it again to refresh the variable menus and date bounds.")
      }
      summary <- state$inspection$summary
      paste(
        "Inspection complete",
        paste("Source type:", summary$source_type[[1]]),
        paste("Files:", summary$file_count[[1]]),
        paste("Variables:", summary$data_variable_count[[1]], "data variables"),
        paste("Time:", ereefs_gui_range_text(summary$time_start[[1]], summary$time_end[[1]])),
        paste("Longitude:", ereefs_gui_range_text(summary$longitude_min[[1]], summary$longitude_max[[1]])),
        paste("Latitude:", ereefs_gui_range_text(summary$latitude_min[[1]], summary$latitude_max[[1]])),
        sep = "\n"
      )
    })

    output$summary_tbl <- shiny::renderTable({
      shiny::req(state$inspection)
      state$inspection$summary
    }, striped = TRUE, bordered = TRUE, spacing = "s")

    output$variables_tbl <- shiny::renderTable({
      shiny::req(state$inspection)
      state$inspection$variables
    }, striped = TRUE, bordered = TRUE, spacing = "s")

    output$dimensions_tbl <- shiny::renderTable({
      shiny::req(state$inspection)
      state$inspection$dimensions
    }, striped = TRUE, bordered = TRUE, spacing = "s")

    output$files_tbl <- shiny::renderTable({
      shiny::req(state$inspection)
      state$inspection$files
    }, striped = TRUE, bordered = TRUE, spacing = "s")

    output$selected_call <- shiny::renderText({
      paste0(
        "inspect_ereefs_data(input_file = ",
        ereefs_gui_r_value(selected_source()),
        ", read_remote_metadata = ",
        ereefs_gui_r_value(isTRUE(input$read_remote_metadata)),
        ")"
      )
    })

    map_args <- shiny::reactive({
      shiny::req(nzchar(input$map_variable %||% ""), nzchar(selected_source()))
      ereefs_gui_map_args(
        input_file = selected_source(),
        var_name = input$map_variable,
        target_date = input$map_date,
        layer = input$map_layer %||% "surface",
        box_bounds = ereefs_gui_bbox(input),
        scale_col = input$map_palette,
        scale_lim = ereefs_gui_scale_lim(input$map_auto_scale, input$map_scale_min, input$map_scale_max),
        plot_style = input$map_style,
        smooth_pixels = input$smooth_pixels,
        label_towns = isTRUE(input$label_towns),
        land_map = isTRUE(input$land_map),
        gbr_poly = isTRUE(input$gbr_poly)
      )
    })

    map_fetch_signature <- shiny::reactive({
      shiny::req(nzchar(input$map_variable %||% ""), nzchar(selected_source()))
      ereefs_gui_map_fetch_signature(
        input_file = selected_source(),
        var_name = input$map_variable,
        target_date = input$map_date,
        layer = input$map_layer %||% "surface",
        box_bounds = ereefs_gui_bbox(input)
      )
    })

    observeEvent(list(input$map_scale_min, input$map_scale_max), {
      if (isTRUE(syncing_map_scale_inputs())) {
        return(NULL)
      }
      if (isTRUE(input$map_auto_scale)) {
        shiny::updateCheckboxInput(session, "map_auto_scale", value = FALSE)
      }
    }, ignoreInit = TRUE)

    observeEvent(list(input$lon_min, input$lon_max, input$lat_min, input$lat_max), {
      if (isTRUE(syncing_bbox_inputs())) {
        return(NULL)
      }
      if (!isTRUE(input$use_bbox)) {
        shiny::updateCheckboxInput(session, "use_bbox", value = TRUE)
      }
    }, ignoreInit = TRUE)

    observeEvent(input$use_bbox, {
      if (isTRUE(input$use_bbox) || is.null(state$inspection)) {
        return(NULL)
      }
      ereefs_gui_sync_bbox_inputs(session, state$inspection, syncing_bbox_inputs)
    }, ignoreInit = TRUE)

    observeEvent(input$map_auto_scale, {
      if (!isTRUE(input$map_auto_scale) || is.null(state$map_result)) {
        return(NULL)
      }
      ereefs_gui_sync_map_scale_inputs(session, state$map_result, syncing_map_scale_inputs)
    }, ignoreInit = TRUE)

    observeEvent(input$preview_map, {
      args <- map_args()
      fetch_signature <- map_fetch_signature()
      state$map_status <- paste(
        "Preparing map preview",
        paste("Variable:", args$var_name),
        paste("Date:", args$target_date),
        paste("Source:", args$input_file),
        sep = "\n"
      )
      state$map_script <- ereefs_gui_script("map", "map_ereefs", args)
      if (!is.null(state$map_result) && identical(state$map_fetch_signature, fetch_signature)) {
        if (isTRUE(input$map_auto_scale)) {
          ereefs_gui_sync_map_scale_inputs(session, state$map_result, syncing_map_scale_inputs)
        }
        state$map_status <- paste(
          "Map preview refreshed from stored data",
          paste("Variable:", args$var_name),
          paste("Date:", args$target_date),
          "No new NetCDF subset was fetched.",
          sep = "\n"
        )
        shiny::showNotification("Map preview refreshed from stored data.", type = "message", duration = 5)
        return(NULL)
      }
      state$map_result <- tryCatch(
        {
          shiny::withProgress(message = "Rendering map preview", value = 0.05, {
            shiny::setProgress(value = 0.15, detail = "Resolving source and matching the requested time")
            state$map_status <- "Resolving source file and matching the requested model time..."
            ereefs_gui_flush(session)
            shiny::setProgress(value = 0.45, detail = "Fetching the requested subset")
            state$map_status <- "Fetching the requested variable subset. Remote OPeNDAP sources can take a little while..."
            ereefs_gui_flush(session)
            fetch_args <- args
            fetch_args$return_poly <- TRUE
            result <- do.call(map_ereefs, fetch_args)
            shiny::setProgress(value = 0.9, detail = "Drawing map")
            state$map_status <- "Rendering the map preview in the browser..."
            ereefs_gui_flush(session)
            result
          })
        },
        error = function(e) {
          state$map_fetch_signature <- NULL
          state$map_status <- paste("Map preview failed:", conditionMessage(e))
          shiny::showNotification(paste("Map preview failed:", conditionMessage(e)), type = "error", duration = 10)
          NULL
        }
      )
      if (!is.null(state$map_result)) {
        state$map_fetch_signature <- fetch_signature
        if (isTRUE(input$map_auto_scale)) {
          ereefs_gui_sync_map_scale_inputs(session, state$map_result, syncing_map_scale_inputs)
        }
        state$map_status <- paste("Map preview ready", paste("Variable:", args$var_name), paste("Date:", args$target_date), sep = "\n")
        shiny::showNotification("Map preview ready. You can now save it or use it as the animation starting point.", type = "message")
      }
    }, ignoreInit = TRUE)

    output$map_plot <- shiny::renderPlot({
      shiny::req(state$map_result)
      print(ereefs_gui_current_map_plot(state$map_result, input))
    }, res = 110)

    output$map_call <- shiny::renderText({
      state$map_script %||% "Preview a map to generate the matching R script."
    })

    output$map_status <- shiny::renderText({
      state$map_status
    })

    output$download_map_png <- shiny::downloadHandler(
      filename = function() paste0("ereefs-map-", Sys.Date(), ".png"),
      content = function(file) {
        shiny::req(state$map_result)
        ggplot2::ggsave(file, plot = ereefs_gui_current_map_plot(state$map_result, input), width = 9, height = 6, dpi = 150)
      }
    )

    output$download_map_script <- shiny::downloadHandler(
      filename = function() paste0("ereefs-map-", Sys.Date(), ".R"),
      content = function(file) writeLines(state$map_script %||% "", file)
    )

    animation_args <- shiny::reactive({
      base <- map_args()
      base$target_date <- NULL
      base$start_date <- input$anim_start
      base$end_date <- input$anim_end
      base$output_dir <- file.path(tempdir(), "ereefs_gui_animation")
      base$save_frames <- TRUE
      base$animation_format <- input$anim_format
      base$animation_file <- NA
      base$fps <- input$anim_fps
      base$keep_frames <- isTRUE(input$anim_keep_frames)
      base$stride <- input$anim_stride
      base
    })

    output$animation_warning <- shiny::renderText({
      args <- animation_args()
      estimate <- ereefs_gui_estimate_frames(args$start_date, args$end_date, args$stride)
      if (estimate > 20 || ereefs_gui_is_remote(args$input_file)) {
        return(paste("Estimated frames:", estimate, "- This may be slow for remote or large requests."))
      }
      paste("Estimated frames:", estimate)
    })

    observeEvent(input$render_animation, {
      args <- animation_args()
      state$animation_result <- NULL
      state$animation_frame_file <- NULL
      state$animation_frame_src <- NULL
      state$animation_media_src <- NULL
      state$animation_media_type <- NULL
      state$animation_progress <- "Preparing animation..."
      args$progress_callback <- function(info) {
        state$animation_progress <- ereefs_gui_animation_progress_text(info)
        if (identical(info$stage, "saved") && !is.na(info$file) && file.exists(info$file)) {
          state$animation_frame_file <- info$file
          ereefs_gui_register_animation_frame(session, state, info$file)
        }
        shiny::setProgress(
          value = ereefs_gui_animation_progress_fraction(info),
          detail = state$animation_progress
        )
        ereefs_gui_flush(session)
      }
      state$animation_result <- tryCatch(
        {
          shiny::withProgress(message = "Creating animation", value = 0.1, {
            do.call(map_ereefs_movie, args)
          })
        },
        error = function(e) {
          state$animation_progress <- paste("Animation failed:", conditionMessage(e))
          shiny::showNotification(paste("Animation failed:", conditionMessage(e)), type = "error", duration = 10)
          NULL
        }
      )
      script_args <- args
      script_args$progress_callback <- NULL
      state$animation_script <- ereefs_gui_script("animation", "map_ereefs_movie", script_args)
      if (!is.null(state$animation_result)) {
        state$animation_progress <- "Animation complete."
        ereefs_gui_register_animation_media(session, state, state$animation_result$animation_file)
      }
    }, ignoreInit = TRUE)

    output$animation_status <- shiny::renderText({
      if (is.null(state$animation_result)) {
        return(state$animation_progress %||% "Render an animation after checking the frame estimate.")
      }
      animation_file <- state$animation_result$animation_file %||% NA_character_
      frame_count <- length(state$animation_result$frame_files %||% character())
      paste(
        state$animation_progress,
        "Animation complete.",
        paste("Animation file:", ifelse(is.na(animation_file), "not assembled", animation_file)),
        paste("Frame files retained:", frame_count),
        sep = "\n"
      )
    })

    output$animation_preview <- shiny::renderUI({
      if (identical(state$animation_media_type, "gif")) {
        return(shiny::tags$img(src = state$animation_media_src, class = "animation-media", alt = "Generated animation GIF"))
      }
      if (identical(state$animation_media_type, "mp4")) {
        return(shiny::tags$video(src = state$animation_media_src, class = "animation-media", controls = NA))
      }
      if (!is.null(state$animation_frame_src)) {
        return(shiny::tags$img(src = state$animation_frame_src, class = "animation-media", alt = "Latest generated animation frame"))
      }
      shiny::tags$p(class = "help-text", "No generated frame or final animation is available to display yet.")
    })

    output$download_animation_script <- shiny::downloadHandler(
      filename = function() paste0("ereefs-animation-", Sys.Date(), ".R"),
      content = function(file) writeLines(state$animation_script %||% "", file)
    )

    points_tbl <- shiny::reactive({
      ereefs_gui_points(input)
    })

    output$points_preview <- shiny::renderTable({
      points_tbl()
    }, striped = TRUE, bordered = TRUE, spacing = "s")

    output$points_preview_points <- shiny::renderTable({
      points_tbl()
    }, striped = TRUE, bordered = TRUE, spacing = "s")

    output$points_warning <- shiny::renderText({
      points <- points_tbl()
      if (!nrow(points)) {
        return("Add at least one latitude/longitude pair or upload a CSV.")
      }
      date_range <- ereefs_gui_points_date_range(input, state$inspection)
      days <- abs(as.numeric(date_range$end - date_range$start)) + 1
      estimated <- nrow(points) * days
      if (estimated > 60 || ereefs_gui_is_remote(selected_source())) {
        return(paste("Estimated point-days:", estimated, "- This may be slow for remote or large requests."))
      }
      paste("Estimated point-days:", estimated)
    })

    observeEvent(input$extract_points, {
      points <- points_tbl()
      shiny::req(nrow(points) > 0)
      date_range <- ereefs_gui_points_date_range(input, state$inspection)
      days <- abs(as.numeric(date_range$end - date_range$start)) + 1
      estimated <- nrow(points) * days
      args <- list(
        var_names = input$points_variable,
        geocoordinates = points[, c("latitude", "longitude"), drop = FALSE],
        layer = input$points_layer %||% "surface",
        start_date = date_range$start,
        end_date = date_range$end,
        input_file = selected_source()
      )
      state$points_result <- tryCatch(
        {
          shiny::withProgress(message = "Extracting point time series", value = 0.2, {
            args$progress_callback <- function(info) {
              shiny::setProgress(
                value = ereefs_gui_extraction_progress_fraction(info, start = 0.2, span = 0.75),
                detail = ereefs_gui_progress_detail(info)
              )
              ereefs_gui_flush(session)
            }
            result <- do.call(get_ereefs_ts, args)
            result_tbl <- dplyr::as_tibble(result)
            if ("label" %in% names(points) && all(c("latitude", "longitude") %in% names(result_tbl))) {
              result_tbl <- result_tbl %>%
                dplyr::left_join(
                  points %>% dplyr::select(latitude, longitude, label) %>% dplyr::distinct(),
                  by = c("latitude", "longitude")
                )
            }
            result_tbl
          })
        },
        error = function(e) {
          shiny::showNotification(paste("Point extraction failed:", conditionMessage(e)), type = "error", duration = 10)
          NULL
        }
      )
      state$points_script <- ereefs_gui_script("point time series", "get_ereefs_ts", args)
    }, ignoreInit = TRUE)

    output$points_plot <- shiny::renderPlot({
      shiny::req(state$points_result)
      ereefs_gui_ts_plot(state$points_result, input$points_variable)
    }, res = 110)

    output$points_tbl <- shiny::renderTable({
      shiny::req(state$points_result)
      utils::head(as.data.frame(state$points_result), 50)
    }, striped = TRUE, bordered = TRUE, spacing = "s")

    output$points_call <- shiny::renderText({
      state$points_script %||% "Extract points to generate the matching R script."
    })

    output$download_points_csv <- shiny::downloadHandler(
      filename = function() paste0("ereefs-point-timeseries-", Sys.Date(), ".csv"),
      content = function(file) utils::write.csv(state$points_result, file, row.names = FALSE)
    )

    output$download_points_png <- shiny::downloadHandler(
      filename = function() paste0("ereefs-point-timeseries-", Sys.Date(), ".png"),
      content = function(file) {
        grDevices::png(file, width = 1200, height = 800, res = 150)
        on.exit(grDevices::dev.off(), add = TRUE)
        print(ereefs_gui_ts_plot(state$points_result, input$points_variable))
      }
    )

    output$download_points_script <- shiny::downloadHandler(
      filename = function() paste0("ereefs-point-timeseries-", Sys.Date(), ".R"),
      content = function(file) writeLines(state$points_script %||% "", file)
    )

    observeEvent(input$extract_profile, {
      points <- points_tbl()
      shiny::req(nrow(points) > 0)
      date_range <- ereefs_gui_points_date_range(input, state$inspection)
      base_args <- list(
        var_names = input$points_variable,
        start_date = date_range$start,
        end_date = date_range$end,
        input_file = selected_source()
      )
      state$profile_result <- tryCatch(
        {
          shiny::withProgress(message = "Extracting vertical profile", value = 0.2, {
            if (nrow(points) == 1) {
              return(do.call(
                get_ereefs_profile,
                c(
                  base_args,
                  list(
                    geolocation = c(points$latitude[[1]], points$longitude[[1]]),
                    progress_callback = function(info) {
                      shiny::setProgress(
                        value = ereefs_gui_extraction_progress_fraction(info, start = 0.2, span = 0.75),
                        detail = ereefs_gui_progress_detail(info)
                      )
                      ereefs_gui_flush(session)
                    }
                  )
                )
              ))
            }
            profiles <- lapply(seq_len(nrow(points)), function(row_id) {
              shiny::setProgress(
                value = min(0.95, 0.2 + 0.75 * row_id / nrow(points)),
                detail = paste(
                  "Point", row_id, "of", nrow(points),
                  "-",
                  ereefs_gui_point_labels(points)[[row_id]]
                )
              )
              progress_callback <- function(info) {
                point_prefix <- paste(
                  "Point", row_id, "of", nrow(points),
                  "-",
                  ereefs_gui_point_labels(points)[[row_id]]
                )
                shiny::setProgress(
                  value = min(
                    0.95,
                    0.2 + 0.75 * ((row_id - 1) + ereefs_gui_extraction_progress_fraction(info, start = 0, span = 1)) / nrow(points)
                  ),
                  detail = paste(point_prefix, "\n", ereefs_gui_progress_detail(info))
                )
                ereefs_gui_flush(session)
              }
              do.call(
                get_ereefs_profile,
                c(
                  base_args,
                  list(
                    geolocation = c(points$latitude[[row_id]], points$longitude[[row_id]]),
                    progress_callback = progress_callback
                  )
                )
              )
            })
            list(
              type = "multi",
              profiles = profiles,
              labels = ereefs_gui_point_labels(points),
              target_date = date_range$start
            )
          })
        },
        error = function(e) {
          shiny::showNotification(paste("Profile extraction failed:", conditionMessage(e)), type = "error", duration = 10)
          NULL
        }
      )
      state$profile_script <- ereefs_gui_profile_script(points, base_args)
    }, ignoreInit = TRUE)

    output$profile_plot <- shiny::renderPlot({
      shiny::req(state$profile_result)
      print(ereefs_gui_profile_plot(state$profile_result, input$points_variable))
    }, res = 110)

    observeEvent(input$run_slice, {
      coords <- ereefs_gui_transect_points(input$transect_points)
      shiny::req(nrow(coords) >= 2)
      args <- list(
        var_names = input$slice_variable,
        geolocation = coords,
        target_date = input$slice_date,
        input_file = selected_source(),
        robust = FALSE
      )
      state$slice_result <- tryCatch(
        {
          shiny::withProgress(message = "Extracting vertical slice", value = 0.2, {
            args$progress_callback <- function(info) {
              shiny::setProgress(
                value = ereefs_gui_extraction_progress_fraction(info, start = 0.2, span = 0.75),
                detail = ereefs_gui_progress_detail(info)
              )
              ereefs_gui_flush(session)
            }
            do.call(get_ereefs_slice, args)
          })
        },
        error = function(e) {
          shiny::showNotification(paste("Slice extraction failed:", conditionMessage(e)), type = "error", duration = 10)
          NULL
        }
      )
      state$slice_script <- ereefs_gui_script("vertical slice", "get_ereefs_slice", args)
    }, ignoreInit = TRUE)

    output$slice_plot <- shiny::renderPlot({
      shiny::req(state$slice_result)
      print(plot_ereefs_slice(state$slice_result, var_name = input$slice_variable, scale_col = input$slice_palette))
    }, res = 110)

    output$slice_call <- shiny::renderText({
      state$slice_script %||% "Run a slice to generate the matching R script."
    })

    output$download_slice_png <- shiny::downloadHandler(
      filename = function() paste0("ereefs-vertical-slice-", Sys.Date(), ".png"),
      content = function(file) {
        shiny::req(state$slice_result)
        p <- plot_ereefs_slice(state$slice_result, var_name = input$slice_variable, scale_col = input$slice_palette)
        ggplot2::ggsave(file, plot = p, width = 9, height = 5, dpi = 150)
      }
    )

    output$download_slice_script <- shiny::downloadHandler(
      filename = function() paste0("ereefs-vertical-slice-", Sys.Date(), ".R"),
      content = function(file) writeLines(state$slice_script %||% "", file)
    )
  }

  shiny::shinyApp(ui = ui, server = server)
}

ereefs_gui_card <- function(title, ..., extra_class = "") {
  shiny::div(
    class = paste("ereefs-card", extra_class),
    shiny::div(class = "ereefs-card-header", title),
    shiny::div(class = "ereefs-card-body", ...)
  )
}

ereefs_gui_source_tab <- function(input_file = "") {
  shiny::fluidRow(
    shiny::column(
      5,
      ereefs_gui_card(
        "Choose a dataset",
      shiny::radioButtons(
        "source_mode",
        "Source type",
        choices = c(
          "Type file path, URL, catalog, or shortcut" = "manual",
          "Browse local NetCDF or mnc file" = "browse",
          "Browse NCI catalogs" = "nci"
        ),
        selected = "manual"
      ),
      shiny::conditionalPanel(
        "input.source_mode == 'manual'",
        shiny::textAreaInput(
          "manual_source",
          "Local path, NetCDF/OPeNDAP/THREDDS URL, or shortcut",
          value = input_file,
          rows = 4,
          placeholder = "catalog, nci, C:/path/to/file.nc, C:/path/to/list.mnc, an OPeNDAP URL, or a THREDDS catalog URL"
        )
      ),
      shiny::conditionalPanel(
        "input.source_mode == 'browse'",
        shiny::p(class = "help-text", "Large-file safe browser: filenames are shown here in the app and only the selected path is passed to R."),
        shiny::textInput("local_browser_dir_input", "Current folder", value = ereefs_gui_default_local_browser_dir()),
        shiny::div(
          class = "button-row",
          shiny::actionButton("local_browser_go", "Go to folder", class = "btn-secondary"),
          shiny::actionButton("local_browser_up", "Up one folder", class = "btn-outline-secondary"),
          shiny::actionButton("local_browser_home", "Home", class = "btn-outline-secondary")
        ),
        shiny::uiOutput("local_browser_entries"),
        shiny::div(
          class = "button-row",
          shiny::actionButton("local_browser_open", "Open folder / preview file", class = "btn-secondary"),
          shiny::actionButton("local_browser_select", "Use selected file", class = "btn-primary")
        ),
        shiny::textInput("local_file_path_display", "Selected local file path", value = "", placeholder = "No local file selected yet"),
        shiny::p(class = "help-text", "You can also paste a full file path in the manual source box. Paths are used directly and are not uploaded.")
      ),
      shiny::conditionalPanel(
        "input.source_mode == 'nci'",
        shiny::div(
          class = "button-row",
          shiny::actionButton("nci_help", "Naming help", class = "btn-outline-secondary")
        ),
        shiny::selectInput("nci_catalog", "Available NCI eReefs catalogs", choices = character())
      ),
      shiny::checkboxInput(
        "read_remote_metadata",
        "Also read full remote grid/time metadata (slower)",
        value = FALSE
      ),
      shiny::numericInput("max_files", "Maximum catalog files to list on inspection page", value = 25, min = 1, max = 500, step = 1),
      shiny::actionButton("inspect_btn", "Choose dataset", class = "btn-primary btn-lg"),
      extra_class = "source-card"
    ),
    ),
    shiny::column(
      7,
      ereefs_gui_card(
        "Metadata preview",
      shiny::verbatimTextOutput("source_status"),
      shiny::tags$details(
        shiny::tags$summary("Show R code"),
        shiny::verbatimTextOutput("selected_call")
      ),
      shiny::hr(),
      shiny::h5("Summary"),
      shiny::tableOutput("summary_tbl"),
      shiny::h5("Variables"),
      shiny::tableOutput("variables_tbl"),
      shiny::h5("Dimensions"),
      shiny::tableOutput("dimensions_tbl"),
      shiny::h5("Catalog files"),
      shiny::tableOutput("files_tbl")
    )
    )
  )
}

ereefs_gui_map_tab <- function() {
  shiny::fluidRow(
    shiny::column(
      4,
      ereefs_gui_card(
        "Map controls",
      shiny::selectInput("map_variable", "Variable", choices = character()),
      shiny::dateInput("map_date", "Time step", value = Sys.Date()),
      shiny::conditionalPanel(
        "output.has_multiple_layers == 'true'",
        shiny::div(
          class = "layer-row",
          shiny::div(
            class = "layer-input-wrap",
            shiny::textInput("map_layer", "Layer", value = "surface")
          ),
          shiny::div(
            class = "layer-hint-wrap",
            shiny::tags$label(class = "layer-hint-label", "Selected variable"),
            shiny::textOutput("map_layer_hint")
          )
        )
      ),
      shiny::checkboxInput("use_bbox", "Use bounding box", value = FALSE),
      shiny::fluidRow(
        shiny::column(
          6,
        shiny::numericInput("lon_min", "Lon min", value = NA),
        shiny::numericInput("lon_max", "Lon max", value = NA)
        ),
        shiny::column(
          6,
        shiny::numericInput("lat_min", "Lat min", value = NA),
        shiny::numericInput("lat_max", "Lat max", value = NA)
        )
      ),
      shiny::selectInput("map_palette", "Palette", choices = c("viridis", "magma", "inferno", "plasma", "cividis", "turbo", "spectral"), selected = "viridis"),
      shiny::checkboxInput("map_auto_scale", "Auto colour range", value = TRUE),
      shiny::fluidRow(
        shiny::column(
          6,
        shiny::numericInput("map_scale_min", "Colour min", value = NA),
        ),
        shiny::column(
          6,
        shiny::numericInput("map_scale_max", "Colour max", value = NA)
        )
      ),
      shiny::selectInput("map_style", "Map style", choices = c("smooth", "polygon"), selected = "polygon"),
      shiny::numericInput("smooth_pixels", "Smooth pixels", value = 600, min = 100, max = 2000, step = 100),
      shiny::checkboxInput("label_towns", "Label towns", value = TRUE),
      shiny::checkboxInput("land_map", "Land underlay", value = FALSE),
      shiny::checkboxInput("gbr_poly", "GBR polygon overlay", value = FALSE),
      shiny::actionButton("preview_map", "Preview map", class = "btn-primary")
    ),
    ),
    shiny::column(
      8,
      ereefs_gui_card(
        "Map preview",
      shiny::p(class = "help-text", "Use the Points tab to enter coordinates manually or upload a CSV of point locations."),
      shiny::verbatimTextOutput("map_status"),
      shiny::plotOutput("map_plot", height = "680px"),
      shiny::div(
        class = "button-row",
        shiny::downloadButton("download_map_png", "Save PNG"),
        shiny::downloadButton("download_map_script", "Save R script")
      ),
      shiny::tags$details(
        shiny::tags$summary("Show R code"),
        shiny::verbatimTextOutput("map_call")
      ),
      extra_class = "map-card"
    )
    )
  )
}

ereefs_gui_animation_tab <- function() {
  shiny::fluidRow(
    shiny::column(
      4,
      ereefs_gui_card(
        "Animation controls",
      shiny::p("Animation starts from the current map settings, then extends them through time."),
      shiny::dateInput("anim_start", "Start date", value = Sys.Date()),
      shiny::dateInput("anim_end", "End date", value = Sys.Date() + 2),
      shiny::selectInput("anim_stride", "Time stride", choices = c("daily", "weekly", "monthly"), selected = "daily"),
      shiny::selectInput("anim_format", "Format", choices = c("gif", "mp4", "none"), selected = "gif"),
      shiny::numericInput("anim_fps", "Frames per second", value = 2, min = 1, max = 30, step = 1),
      shiny::checkboxInput("anim_keep_frames", "Keep individual frame images", value = FALSE),
      shiny::verbatimTextOutput("animation_warning"),
      shiny::actionButton("render_animation", "Render animation", class = "btn-primary")
    ),
    ),
    shiny::column(
      8,
      ereefs_gui_card(
        "Animation output",
      shiny::verbatimTextOutput("animation_status"),
      shiny::h5("Animation preview"),
      shiny::uiOutput("animation_preview"),
      shiny::downloadButton("download_animation_script", "Save R script")
    )
    )
  )
}

ereefs_gui_points_tab <- function() {
  shiny::fluidRow(
    shiny::column(
      4,
      ereefs_gui_card(
        "Point selection",
      shiny::selectInput("points_variable", "Variable", choices = character()),
      shiny::conditionalPanel(
        "output.has_multiple_layers == 'true'",
        shiny::div(
          class = "layer-row",
          shiny::div(
            class = "layer-input-wrap",
            shiny::textInput("points_layer", "Layer", value = "surface")
          ),
          shiny::div(
            class = "layer-hint-wrap",
            shiny::tags$label(class = "layer-hint-label", "Selected variable"),
            shiny::textOutput("points_layer_hint")
          )
        )
      ),
      shiny::fluidRow(
        shiny::column(
          7,
          shiny::dateInput("points_start", "Start date", value = Sys.Date())
        ),
        shiny::column(
          5,
          shiny::tags$label(class = "layer-hint-label", "Selected file/catalog"),
          shiny::textOutput("points_start_range")
        )
      ),
      shiny::fluidRow(
        shiny::column(
          7,
          shiny::dateInput("points_end", "End date", value = Sys.Date() + 2)
        ),
        shiny::column(
          5,
          shiny::tags$label(class = "layer-hint-label", "Selected file/catalog"),
          shiny::textOutput("points_end_range")
        )
      ),
      shiny::textAreaInput(
        "points_text",
        "Manual points as latitude,longitude,label",
        value = "-19.26,146.82,Townsville shelf",
        rows = 4
      ),
      shiny::fileInput("points_csv", "Or upload CSV with latitude and longitude columns", accept = ".csv"),
      shiny::p(
        class = "help-text",
        "Tip: point selection is currently manual or CSV-based while map-click selection is disabled."
      ),
      shiny::h5("Extraction points"),
      shiny::tableOutput("points_preview_points"),
      shiny::verbatimTextOutput("points_warning"),
      shiny::div(
        class = "button-row",
        shiny::actionButton("extract_points", "Extract time series", class = "btn-primary"),
        shiny::conditionalPanel(
          "output.has_multiple_layers == 'true'",
          shiny::actionButton("extract_profile", "Extract profile(s)", class = "btn-secondary")
        )
      )
    ),
    ),
    shiny::column(
      8,
      ereefs_gui_card(
        "Point outputs",
      shiny::plotOutput("points_plot", height = "360px"),
      shiny::plotOutput("profile_plot", height = "360px"),
      shiny::div(
        class = "button-row",
        shiny::downloadButton("download_points_csv", "Save CSV"),
        shiny::downloadButton("download_points_png", "Save plot PNG"),
        shiny::downloadButton("download_points_script", "Save R script")
      ),
      shiny::tags$details(
        shiny::tags$summary("Time-series table"),
        shiny::tableOutput("points_tbl")
      ),
      shiny::tags$details(
        shiny::tags$summary("Show R code"),
        shiny::verbatimTextOutput("points_call")
      )
    )
    )
  )
}

ereefs_gui_transect_tab <- function() {
  shiny::fluidRow(
    shiny::column(
      4,
      ereefs_gui_card(
        "Transect or vertical slice",
      shiny::selectInput("slice_variable", "Variable", choices = character()),
      shiny::fluidRow(
        shiny::column(
          7,
          shiny::dateInput("slice_date", "Date", value = Sys.Date())
        ),
        shiny::column(
          5,
          shiny::tags$label(class = "layer-hint-label", "Selected file/catalog"),
          shiny::textOutput("slice_date_range")
        )
      ),
      shiny::textAreaInput(
        "transect_points",
        "Path as latitude,longitude per line",
        value = "-19.25,146.82\n-18.85,147.35",
        rows = 5
      ),
      shiny::selectInput("slice_palette", "Palette", choices = c("viridis", "magma", "inferno", "plasma", "cividis", "turbo", "spectral"), selected = "viridis"),
      shiny::actionButton("run_slice", "Run vertical slice", class = "btn-primary")
    ),
    ),
    shiny::column(
      8,
      ereefs_gui_card(
        "Slice output",
      shiny::plotOutput("slice_plot", height = "560px"),
      shiny::div(
        class = "button-row",
        shiny::downloadButton("download_slice_png", "Save PNG"),
        shiny::downloadButton("download_slice_script", "Save R script")
      ),
      shiny::tags$details(
        shiny::tags$summary("Show R code"),
        shiny::verbatimTextOutput("slice_call")
      )
    )
    )
  )
}

ereefs_gui_notes_tab <- function() {
  ereefs_gui_card(
    "Implementation notes",
    shiny::tags$p("This GUI is intentionally metadata-first. Expensive OPeNDAP reads only happen after clicking a plot, animation, or extraction button."),
    shiny::tags$p("Generated R scripts are designed to reproduce GUI actions from the command line without storing credentials or Shiny session state."),
    shiny::tags$p("Progress notes for development are saved in inst/GUI_PROGRESS_NOTES.md.")
  )
}

ereefs_gui_css <- function() {
  "
  body {
    font-family: Aptos, Candara, Calibri, Segoe UI, sans-serif;
    background:
      radial-gradient(circle at 10% 10%, rgba(216, 139, 56, 0.18), transparent 30rem),
      radial-gradient(circle at 90% 0%, rgba(11, 111, 120, 0.20), transparent 34rem),
      linear-gradient(135deg, #f7fbf8 0%, #eef7f5 45%, #fff8ed 100%);
  }
  .ereefs-hero {
    display: flex;
    justify-content: space-between;
    gap: 2rem;
    margin: 1.2rem 0 1.4rem;
    padding: 2rem;
    border-radius: 2rem;
    color: #f8fffb;
    background:
      linear-gradient(135deg, rgba(5, 55, 64, 0.96), rgba(9, 113, 123, 0.92)),
      repeating-linear-gradient(120deg, rgba(255,255,255,0.08) 0 1px, transparent 1px 18px);
    box-shadow: 0 26px 70px rgba(5, 55, 64, 0.22);
  }
  .ereefs-hero h1 {
    font-family: Georgia, Cambria, Times New Roman, serif;
    font-size: clamp(2.4rem, 6vw, 5rem);
    line-height: 0.92;
    margin: 0.2rem 0 0.8rem;
  }
  .ereefs-hero p {
    max-width: 55rem;
    font-size: 1.1rem;
    opacity: 0.95;
  }
  .ereefs-logo {
    display: block;
    width: min(210px, 55vw);
    margin-bottom: 1rem;
    padding: 0.7rem 0.9rem;
    border-radius: 1.1rem;
    background: rgba(255,255,255,0.92);
    box-shadow: 0 14px 30px rgba(0,0,0,0.14);
  }
  .eyebrow {
    letter-spacing: 0.16em;
    text-transform: uppercase;
    font-weight: 700;
    color: #ffd99a;
  }
  .ereefs-hero-badges {
    display: flex;
    flex-direction: column;
    gap: 0.7rem;
    align-self: center;
  }
  .ereefs-hero-badges span {
    border: 1px solid rgba(255,255,255,0.35);
    border-radius: 999px;
    padding: 0.45rem 0.85rem;
    background: rgba(255,255,255,0.10);
    white-space: nowrap;
  }
  .ereefs-card {
    border: 0;
    border-radius: 1.4rem;
    box-shadow: 0 18px 45px rgba(38, 82, 84, 0.12);
    overflow: hidden;
    background: rgba(255,255,255,0.92);
    margin-bottom: 1.25rem;
  }
  .ereefs-card-header {
    font-weight: 800;
    background: linear-gradient(90deg, rgba(11,111,120,0.10), rgba(216,139,56,0.10));
    padding: 0.95rem 1.15rem;
    border-bottom: 1px solid rgba(11,111,120,0.10);
  }
  .ereefs-card-body {
    padding: 1.15rem;
  }
  .nav-tabs {
    border: 0;
    gap: 0.35rem;
    margin-bottom: 1rem;
  }
  .nav-tabs .nav-link {
    border: 0;
    border-radius: 999px;
    color: #145a5f;
    background: rgba(255,255,255,0.72);
    font-weight: 700;
  }
  .nav-tabs .nav-link.active {
    color: white;
    background: linear-gradient(135deg, #0b6f78, #d88b38);
  }
  .button-row {
    display: flex;
    flex-wrap: wrap;
    gap: 0.65rem;
    align-items: center;
    margin: 0.6rem 0;
  }
  pre {
    border-radius: 1rem;
    background: #092f35;
    color: #e4fff5;
    padding: 1rem;
  }
  table {
    font-size: 0.92rem;
  }
  .help-text {
    color: #486368;
    font-size: 0.92rem;
  }
  .layer-row {
    display: flex;
    gap: 0.9rem;
    align-items: flex-end;
  }
  .layer-input-wrap {
    flex: 1 1 auto;
    min-width: 0;
  }
  .layer-hint-wrap {
    flex: 0 0 11rem;
    padding-bottom: 0.75rem;
    color: #486368;
    font-size: 0.9rem;
  }
  .layer-hint-label {
    display: block;
    margin-bottom: 0.2rem;
    font-weight: 700;
    color: #0b4d54;
  }
  .layer-hint-wrap .shiny-text-output {
    line-height: 1.25;
  }
  .source-ribbon {
    display: grid;
    grid-template-columns: auto 1fr;
    align-items: start;
    gap: 0.8rem;
    margin: -0.35rem 0 1rem;
    padding: 0.9rem 1.1rem;
    border: 1px solid rgba(11, 111, 120, 0.20);
    border-radius: 1rem;
    background: rgba(255, 255, 255, 0.72);
    box-shadow: 0 12px 30px rgba(5, 55, 64, 0.08);
  }
  .source-ribbon pre {
    margin: 0;
    padding: 0;
    border: 0;
    background: transparent;
    white-space: pre-wrap;
    word-break: break-word;
    color: #244448;
  }
  details summary {
    display: inline-flex;
    align-items: center;
    gap: 0.35rem;
    margin: 0.35rem 0 0.5rem;
    padding: 0.45rem 0.75rem;
    border-radius: 999px;
    border: 1px solid rgba(11, 111, 120, 0.28);
    background: linear-gradient(135deg, rgba(255,255,255,0.95), rgba(223,242,238,0.92));
    color: #0b4d54;
    cursor: pointer;
    font-weight: 700;
    box-shadow: 0 6px 16px rgba(5, 55, 64, 0.08);
    user-select: none;
  }
  details summary::after {
    content: 'click to expand';
    font-size: 0.78rem;
    font-weight: 600;
    color: #6f7d7f;
  }
  details[open] summary::after {
    content: 'click to hide';
  }
  details summary:hover {
    transform: translateY(-1px);
    border-color: rgba(216, 139, 56, 0.65);
    box-shadow: 0 9px 20px rgba(5, 55, 64, 0.12);
  }
  "
}

ereefs_gui_update_from_inspection <- function(session, inspection) {
  choices <- ereefs_gui_variable_choices(inspection)
  data_vars <- choices$data_vars
  map_vars <- choices$map_vars
  summary <- inspection$summary
  date_bounds <- ereefs_gui_date_bounds(summary)
  first_date <- date_bounds$start
  end_date <- date_bounds$end
  default_end <- min(first_date + 2, end_date)
  shiny::updateSelectInput(session, "map_variable", choices = map_vars, selected = map_vars[[1]])
  shiny::updateSelectInput(session, "points_variable", choices = data_vars, selected = data_vars[[1]])
  shiny::updateSelectInput(session, "slice_variable", choices = data_vars, selected = data_vars[[1]])
  shiny::updateDateInput(session, "map_date", value = first_date, min = first_date, max = end_date)
  shiny::updateDateInput(session, "anim_start", value = first_date, min = first_date, max = end_date)
  shiny::updateDateInput(session, "anim_end", value = default_end, min = first_date, max = end_date)
  shiny::updateDateInput(session, "points_start", value = first_date, min = first_date, max = end_date)
  shiny::updateDateInput(session, "points_end", value = default_end, min = first_date, max = end_date)
  shiny::updateDateInput(session, "slice_date", value = first_date, min = first_date, max = end_date)
  shiny::updateCheckboxInput(session, "use_bbox", value = FALSE)
  ereefs_gui_sync_bbox_inputs(session, inspection)
}

ereefs_gui_reset_controls <- function(session) {
  default_vars <- ereefs_gui_common_variable_choices()
  today <- Sys.Date()
  shiny::updateSelectInput(session, "map_variable", choices = default_vars, selected = default_vars[[1]])
  shiny::updateSelectInput(session, "points_variable", choices = default_vars, selected = default_vars[[1]])
  shiny::updateSelectInput(session, "slice_variable", choices = default_vars, selected = default_vars[[1]])
  shiny::updateDateInput(session, "map_date", value = today)
  shiny::updateDateInput(session, "anim_start", value = today)
  shiny::updateDateInput(session, "anim_end", value = today + 2)
  shiny::updateDateInput(session, "points_start", value = today)
  shiny::updateDateInput(session, "points_end", value = today + 2)
  shiny::updateDateInput(session, "slice_date", value = today)
}

ereefs_gui_variable_choices <- function(inspection) {
  data_vars <- inspection$variables$variable[inspection$variables$is_data_variable]
  if (!length(data_vars)) {
    data_vars <- inspection$variables$variable
  }
  if (!length(data_vars)) {
    data_vars <- ereefs_gui_common_variable_choices()
  }
  data_vars <- unique(data_vars)
  list(
    data_vars = data_vars,
    map_vars = unique(c(data_vars, ereefs_gui_available_special_map_variables(inspection)))
  )
}

ereefs_gui_sync_variable_inputs <- function(session, inspection, selected_var, syncing_flag = NULL) {
  if (!nzchar(selected_var %||% "")) {
    return(invisible(NULL))
  }
  if (is.null(inspection)) {
    return(invisible(NULL))
  }
  choices <- ereefs_gui_variable_choices(inspection)
  if (is.function(syncing_flag)) {
    syncing_flag(TRUE)
    on.exit(syncing_flag(FALSE), add = TRUE)
  }
  shiny::updateSelectInput(
    session,
    "map_variable",
    choices = choices$map_vars,
    selected = ereefs_gui_select_synced_variable(selected_var, choices$map_vars)
  )
  shiny::updateSelectInput(
    session,
    "points_variable",
    choices = choices$data_vars,
    selected = ereefs_gui_select_synced_variable(selected_var, choices$data_vars)
  )
  shiny::updateSelectInput(
    session,
    "slice_variable",
    choices = choices$data_vars,
    selected = ereefs_gui_select_synced_variable(selected_var, choices$data_vars)
  )
  invisible(selected_var)
}

ereefs_gui_select_synced_variable <- function(selected_var, choices) {
  if (selected_var %in% choices) {
    return(selected_var)
  }
  choices[[1]]
}

ereefs_gui_date_bounds <- function(summary) {
  start <- as.Date(summary$time_start[[1]])
  end <- as.Date(summary$time_end[[1]])
  if (is.na(start)) {
    start <- as.Date(summary$file_start[[1]])
  }
  if (is.na(end)) {
    end <- as.Date(summary$file_end[[1]])
  }
  if (is.na(start)) {
    start <- Sys.Date()
  }
  if (is.na(end)) {
    end <- start
  }
  if (end < start) {
    end <- start
  }
  list(start = start, end = end)
}

ereefs_gui_date_range_label <- function(inspection = NULL) {
  if (is.null(inspection) || is.null(inspection$summary)) {
    return("Range unavailable")
  }
  bounds <- ereefs_gui_date_bounds(inspection$summary)
  paste("Valid:", format(bounds$start), "to", format(bounds$end))
}

ereefs_gui_points_date_range <- function(input, inspection = NULL) {
  start <- ereefs_gui_date_scalar(input$points_start)
  end <- ereefs_gui_date_scalar(input$points_end)
  if ((is.na(start) || is.na(end)) && !is.null(inspection)) {
    fallback <- ereefs_gui_date_bounds(inspection$summary)
    if (is.na(start)) {
      start <- fallback$start
    }
    if (is.na(end)) {
      end <- fallback$end
    }
  }
  if (is.na(start) || is.na(end)) {
    stop("Choose valid start and end dates before extracting point time series.")
  }
  if (isTRUE(end < start)) {
    end <- start
  }
  list(start = start, end = end)
}

ereefs_gui_date_scalar <- function(x) {
  if (is.null(x) || length(x) == 0) {
    return(as.Date(NA))
  }
  as.Date(x[[1]])
}

ereefs_gui_clamp_date <- function(value, inspection = NULL) {
  date_value <- ereefs_gui_date_scalar(value)
  if (is.na(date_value) || is.null(inspection) || is.null(inspection$summary)) {
    return(date_value)
  }
  bounds <- ereefs_gui_date_bounds(inspection$summary)
  min(max(date_value, bounds$start), bounds$end)
}

ereefs_gui_clamp_date_input <- function(session, input_id, value, inspection = NULL) {
  clamped <- ereefs_gui_clamp_date(value, inspection)
  if (is.na(clamped)) {
    return(invisible(clamped))
  }
  current <- ereefs_gui_date_scalar(value)
  if (is.na(current) || !identical(current, clamped)) {
    shiny::updateDateInput(session, input_id, value = clamped)
  }
  invisible(clamped)
}

ereefs_gui_summary_bounds <- function(summary) {
  bound <- function(name) {
    value <- suppressWarnings(as.numeric(summary[[name]][[1]]))
    if (is.finite(value)) value else NA_real_
  }
  list(
    lon_min = bound("longitude_min"),
    lon_max = bound("longitude_max"),
    lat_min = bound("latitude_min"),
    lat_max = bound("latitude_max")
  )
}

ereefs_gui_common_variable_choices <- function() {
  c(
    "temp", "salt", "eta", "u", "v", "wspeed_u", "wspeed_v",
    "dens", "NH4", "NO3", "DIP", "Chl_a_sum", "Kd_490", "TSSM"
  )
}

ereefs_gui_available_special_map_variables <- function(inspection) {
  vars <- inspection$variables$variable
  special <- character()
  if (all(c("R_470", "R_555", "R_645") %in% vars)) {
    special <- c(special, "true_colour")
  }
  if (all(c("R_412", "R_443", "R_488", "R_531", "R_547", "R_667", "R_678") %in% vars)) {
    special <- c(special, "plume")
  }
  special
}

ereefs_gui_nci_catalogs <- function() {
  menu_url <- "https://thredds.nci.org.au/thredds/catalog/catalogs/fx3/catalog.xml"
  doc <- tryCatch(xml2::read_xml(menu_url), error = function(e) NULL)
  if (is.null(doc)) {
    return(ereefs_gui_fallback_nci_catalogs())
  }
  refs <- xml2::xml_find_all(doc, ".//*[local-name()='catalogRef']")
  if (!length(refs)) {
    return(ereefs_gui_fallback_nci_catalogs())
  }
  rows <- lapply(refs, function(ref) {
    attrs <- xml2::xml_attrs(ref)
    href <- ereefs_gui_xml_attr(attrs, c("xlink:href", "href"))
    title <- ereefs_gui_xml_attr(attrs, c("xlink:title", "title"))
    if (is.na(href) || !nzchar(href)) {
      return(dplyr::tibble(title = character(), url = character()))
    }
    if (is.na(title) || !nzchar(title)) {
      title <- href
    }
    url <- ereefs_gui_resolve_catalog_href(href, menu_url)
    dplyr::tibble(title = title, url = url)
  })
  dplyr::bind_rows(rows) %>%
    dplyr::filter(!is.na(.data$url), nzchar(.data$url)) %>%
    dplyr::arrange(.data$title)
}

ereefs_gui_resolve_catalog_href <- function(href, base_url) {
  href <- trimws(href %||% "")
  if (!nzchar(href)) {
    return(NA_character_)
  }
  if (grepl("^https?://", href, ignore.case = TRUE)) {
    return(href)
  }
  if (startsWith(href, "/")) {
    return(paste0("https://thredds.nci.org.au", href))
  }
  xml2::url_absolute(href, base_url)
}

ereefs_gui_xml_attr <- function(attrs, names) {
  for (name in names) {
    value <- attrs[name]
    if (length(value) && !is.na(value[[1]]) && nzchar(value[[1]])) {
      return(unname(value[[1]]))
    }
  }
  NA_character_
}

ereefs_gui_fallback_nci_catalogs <- function() {
  dplyr::tibble(
    title = c(
      "GBR4 hydrodynamics v4.0 hindcast",
      "GBR4 biogeochemistry v4.2 baseline catchment hindcast",
      "GBR4 biogeochemistry v4.2 baseline catchment with rivers hindcast",
      "GBR1 hydrodynamics v2.0",
      "GBR1 biogeochemistry and sediments v3.2"
    ),
    url = c(
      "https://thredds.nci.org.au/thredds/catalog/fx3/gbr4_H4p0_ABARRAr2_OBRAN2020_FG2Gv3_Dhnd/catalog.xml",
      "https://thredds.nci.org.au/thredds/catalog/fx3/gbr4_H4p0_ABARRAr2_OBRAN2020_FG2Gv3_B4p2_Cq5b_Dhnd/catalog.xml",
      "https://thredds.nci.org.au/thredds/catalog/fx3/gbr4_H4p0_ABARRAr2_OBRAN2020_FG2Gv3_B4p2_Cq5b_Dhnd_RrHdnd/catalog.xml",
      "https://thredds.nci.org.au/thredds/catalog/fx3/gbr1_2.0/catalog.xml",
      "https://thredds.nci.org.au/thredds/catalog/fx3/gbr1_bgc_GBR1_H2p0_B3p2_Cq3b_Dhnd/catalog.xml"
    )
  )
}

ereefs_gui_add_resource_path <- function() {
  resource_dir <- system.file("app/www", package = "ereefs")
  if (!nzchar(resource_dir)) {
    resource_dir <- file.path(getwd(), "inst", "app", "www")
  }
  if (dir.exists(resource_dir)) {
    shiny::addResourcePath("ereefs-gui", resource_dir)
  }
  invisible(resource_dir)
}

ereefs_gui_map_args <- function(input_file,
                                var_name,
                                target_date,
                                layer,
                                box_bounds,
                                scale_col,
                                scale_lim,
                                plot_style,
                                smooth_pixels,
                                label_towns,
                                land_map,
                                gbr_poly) {
  list(
    var_name = var_name,
    target_date = as.character(target_date),
    layer = ereefs_gui_layer_value(layer),
    input_file = input_file,
    scale_col = scale_col,
    scale_lim = scale_lim,
    plot_style = plot_style,
    smooth_pixels = smooth_pixels,
    box_bounds = box_bounds,
    label_towns = label_towns,
    land_map = land_map,
    gbr_poly = gbr_poly
  )
}

ereefs_gui_map_fetch_signature <- function(input_file,
                                           var_name,
                                           target_date,
                                           layer,
                                           box_bounds) {
  list(
    input_file = input_file,
    var_name = var_name,
    target_date = as.character(target_date),
    layer = ereefs_gui_layer_value(layer),
    box_bounds = box_bounds
  )
}

ereefs_gui_layer_value <- function(layer) {
  layer <- trimws(as.character(layer))
  numeric_layer <- suppressWarnings(as.numeric(layer))
  if (!is.na(numeric_layer) && !layer %in% c("surface", "bottom", "integrated")) {
    return(numeric_layer)
  }
  layer
}

ereefs_gui_bbox <- function(input) {
  if (!isTRUE(input$use_bbox)) {
    return(c(NA_real_, NA_real_, NA_real_, NA_real_))
  }
  bounds <- c(
    ereefs_gui_numeric_scalar(input$lon_min),
    ereefs_gui_numeric_scalar(input$lon_max),
    ereefs_gui_numeric_scalar(input$lat_min),
    ereefs_gui_numeric_scalar(input$lat_max)
  )
  if (any(!is.finite(bounds))) {
    return(c(NA_real_, NA_real_, NA_real_, NA_real_))
  }
  bounds
}

ereefs_gui_scale_lim <- function(auto, lower, upper) {
  lower <- ereefs_gui_numeric_scalar(lower)
  upper <- ereefs_gui_numeric_scalar(upper)
  if (isTRUE(auto) || !is.finite(lower) || !is.finite(upper)) {
    return(c(NA_real_, NA_real_))
  }
  c(lower, upper)
}

ereefs_gui_numeric_scalar <- function(x) {
  if (is.null(x) || length(x) == 0) {
    return(NA_real_)
  }
  value <- suppressWarnings(as.numeric(x[[1]]))
  if (!length(value) || is.na(value)) {
    return(NA_real_)
  }
  value
}

ereefs_gui_sync_bbox_inputs <- function(session, inspection, syncing_flag = NULL) {
  if (is.null(inspection) || is.null(inspection$summary)) {
    return(invisible(NULL))
  }
  bounds <- ereefs_gui_summary_bounds(inspection$summary)
  if (is.function(syncing_flag)) {
    syncing_flag(TRUE)
    on.exit(syncing_flag(FALSE), add = TRUE)
  }
  shiny::updateNumericInput(session, "lon_min", value = bounds$lon_min)
  shiny::updateNumericInput(session, "lon_max", value = bounds$lon_max)
  shiny::updateNumericInput(session, "lat_min", value = bounds$lat_min)
  shiny::updateNumericInput(session, "lat_max", value = bounds$lat_max)
  invisible(bounds)
}

ereefs_gui_map_scale_range <- function(map_result) {
  datapoly <- if (is.list(map_result) && "datapoly" %in% names(map_result)) {
    map_result$datapoly
  } else {
    map_result
  }
  if (!is.data.frame(datapoly) || !"value" %in% names(datapoly)) {
    return(c(NA_real_, NA_real_))
  }
  if (is.character(datapoly$value) || is.factor(datapoly$value)) {
    return(c(NA_real_, NA_real_))
  }
  values <- suppressWarnings(as.numeric(datapoly$value))
  values <- values[is.finite(values)]
  if (!length(values)) {
    return(c(NA_real_, NA_real_))
  }
  c(min(values), max(values))
}

ereefs_gui_sync_map_scale_inputs <- function(session, map_result, syncing_flag = NULL) {
  scale_range <- ereefs_gui_map_scale_range(map_result)
  if (!all(is.finite(scale_range))) {
    return(invisible(scale_range))
  }
  if (is.function(syncing_flag)) {
    syncing_flag(TRUE)
    on.exit(syncing_flag(FALSE), add = TRUE)
  }
  shiny::updateNumericInput(session, "map_scale_min", value = scale_range[[1]])
  shiny::updateNumericInput(session, "map_scale_max", value = scale_range[[2]])
  invisible(scale_range)
}

ereefs_gui_has_multiple_layers <- function(inspection) {
  k_count <- suppressWarnings(as.integer(inspection$summary$k_count[[1]]))
  if (is.finite(k_count)) {
    return(k_count > 1L)
  }
  dims <- inspection$dimensions
  if (!is.null(dims) && nrow(dims)) {
    k_dim <- dims$length[dims$role %in% "k"]
    k_dim <- suppressWarnings(as.integer(k_dim[[1]]))
    if (is.finite(k_dim)) {
      return(k_dim > 1L)
    }
  }
  TRUE
}

ereefs_gui_layer_hint <- function(inspection, variable) {
  if (is.null(inspection) || !nzchar(variable %||% "")) {
    return("")
  }
  vars <- inspection$variables
  hit <- which(vars$variable %in% variable)
  if (!length(hit)) {
    return("")
  }
  roles <- strsplit(vars$dimension_roles[[hit[[1]]]] %||% "", ",", fixed = TRUE)[[1]]
  roles <- trimws(stats::na.omit(roles))
  if (!("k" %in% roles)) {
    return("2D variable")
  }
  k_count <- suppressWarnings(as.integer(inspection$summary$k_count[[1]]))
  if (!is.finite(k_count)) {
    dims <- inspection$dimensions
    if (!is.null(dims) && nrow(dims)) {
      k_dim <- dims$length[dims$role %in% "k"]
      k_count <- suppressWarnings(as.integer(k_dim[[1]]))
    }
  }
  if (!is.finite(k_count) || k_count < 1L) {
    return("Depth resolved")
  }
  paste(k_count, if (k_count == 1L) "layer" else "layers")
}

ereefs_gui_points <- function(input) {
  if (!is.null(input$points_csv)) {
    csv <- tryCatch(utils::read.csv(input$points_csv$datapath, stringsAsFactors = FALSE), error = function(e) NULL)
    if (!is.null(csv) && all(c("latitude", "longitude") %in% names(csv))) {
      out <- dplyr::as_tibble(csv)
      if (!("label" %in% names(out))) {
        out$label <- paste0("point_", seq_len(nrow(out)))
      }
      return(out[, c("latitude", "longitude", "label"), drop = FALSE])
    }
  }
  text <- trimws(input$points_text %||% "")
  if (!nzchar(text)) {
    return(dplyr::tibble(latitude = numeric(), longitude = numeric(), label = character()))
  }
  rows <- strsplit(text, "\n", fixed = TRUE)[[1]]
  parsed <- lapply(seq_along(rows), function(i) {
    bits <- trimws(strsplit(rows[[i]], ",", fixed = TRUE)[[1]])
    if (length(bits) < 2) {
      return(NULL)
    }
    dplyr::tibble(
      latitude = suppressWarnings(as.numeric(bits[[1]])),
      longitude = suppressWarnings(as.numeric(bits[[2]])),
      label = if (length(bits) >= 3 && nzchar(bits[[3]])) bits[[3]] else paste0("point_", i)
    )
  })
  ereefs_gui_point_table(parsed, include_label = TRUE) %>%
    dplyr::filter(is.finite(.data$latitude), is.finite(.data$longitude))
}

ereefs_gui_transect_points <- function(text) {
  rows <- strsplit(trimws(text %||% ""), "\n", fixed = TRUE)[[1]]
  parsed <- lapply(rows, function(row) {
    bits <- trimws(strsplit(row, ",", fixed = TRUE)[[1]])
    if (length(bits) < 2) {
      return(NULL)
    }
    dplyr::tibble(
      latitude = suppressWarnings(as.numeric(bits[[1]])),
      longitude = suppressWarnings(as.numeric(bits[[2]]))
    )
  })
  ereefs_gui_point_table(parsed, include_label = FALSE) %>%
    dplyr::filter(is.finite(.data$latitude), is.finite(.data$longitude))
}

ereefs_gui_point_table <- function(parsed, include_label = TRUE) {
  out <- dplyr::bind_rows(parsed)
  required_cols <- c("latitude", "longitude", if (isTRUE(include_label)) "label")
  if (!all(required_cols %in% names(out))) {
    empty_cols <- lapply(required_cols, function(col) {
      vector(mode = if (identical(col, "label")) "character" else "numeric", length = 0)
    })
    return(dplyr::as_tibble(stats::setNames(empty_cols, required_cols)))
  }
  out[, required_cols, drop = FALSE]
}

ereefs_gui_ts_plot <- function(ts, var_name) {
  data <- dplyr::as_tibble(ts)
  metadata_cols <- c(
    "time", "date", "datetime", "latitude", "longitude", "cell_latitude",
    "cell_longitude", "watercol_depth", "i", "j", "k", "grid_index", "row_id",
    "label"
  )
  value_cols <- setdiff(
    names(data)[vapply(data, is.numeric, logical(1))],
    metadata_cols
  )
  value_col <- if (var_name %in% value_cols) {
    var_name
  } else if (length(value_cols)) {
    value_cols[[1]]
  } else {
    stop("Could not identify a numeric variable column to plot.", call. = FALSE)
  }

  if ("label" %in% names(data) && any(nzchar(data$label %||% ""))) {
    data$point_group <- data$label
  } else if (all(c("latitude", "longitude") %in% names(data))) {
    data$point_group <- sprintf("(%.5f, %.5f)", data$latitude, data$longitude)
  } else if ("row_id" %in% names(data)) {
    data$point_group <- as.character(data$row_id)
  } else {
    data$point_group <- "Selected point"
  }

  p <- ggplot2::ggplot(
    data,
    ggplot2::aes(
      x = .data$time,
      y = .data[[value_col]],
      group = .data$point_group,
      colour = .data$point_group
    )
  ) +
    ggplot2::geom_line() +
    ggplot2::labs(
      x = "Time",
      y = value_col,
      colour = if (length(unique(data$point_group)) > 1) "Point" else NULL
    ) +
    ggplot2::theme_minimal(base_size = 13)

  if (length(unique(data$point_group)) == 1) {
    p <- p + ggplot2::guides(colour = "none")
  }
  p
}

ereefs_gui_profile_values <- function(profile_obj, var_name, target_date) {
  target_date <- get_date_time(target_date)
  if (names(profile_obj)[5] == "profiles") {
    profile_dims <- dim(profile_obj$profiles)
    if (length(profile_dims) > 2) {
      colnum <- which(colnames(profile_obj$profiles) == var_name)
      if (!length(colnum)) {
        return(NULL)
      }
      dind <- which.min(abs(profile_obj$dates - target_date))
      values <- array(profile_obj$profiles[, colnum, dind], length(profile_obj$z_grid) - 1)
      eta <- profile_obj$eta[dind]
    } else if (!is.null(colnames(profile_obj$profiles)) && length(colnames(profile_obj$profiles))) {
      colnum <- which(colnames(profile_obj$profiles) == var_name)
      if (!length(colnum)) {
        return(NULL)
      }
      values <- array(profile_obj$profiles[, colnum])
      eta <- profile_obj$eta
    } else {
      dind <- which.min(abs(profile_obj$dates - target_date))
      values <- array(profile_obj$profiles[, dind], length(profile_obj$z_grid) - 1)
      eta <- profile_obj$eta[dind]
    }
  } else {
    if (is.null(dim(profile_obj$profiles))) {
      values <- profile_obj$profiles
      eta <- profile_obj$eta
    } else {
      dind <- which.min(abs(profile_obj$dates - target_date))
      values <- array(profile_obj$profiles[, dind], length(profile_obj$z_grid) - 1)
      eta <- profile_obj$eta[dind]
    }
  }
  wet <- which(!is.na(values))
  if (!length(wet)) {
    return(NULL)
  }
  dplyr::tibble(
    z = c(profile_obj$botz, profile_obj$z_grid[head(wet, -1L) + 1L], eta),
    value = c(values[wet], values[max(wet)])
  )
}

ereefs_gui_profile_plot <- function(profile_result, var_name) {
  if (is.list(profile_result) && identical(profile_result$type %||% "", "multi")) {
    target_date <- profile_result$target_date
    labels <- profile_result$labels %||% paste0("point_", seq_along(profile_result$profiles))
    data <- lapply(seq_along(profile_result$profiles), function(i) {
      tbl <- ereefs_gui_profile_values(profile_result$profiles[[i]], var_name, target_date)
      if (is.null(tbl)) {
        return(NULL)
      }
      tbl$point_group <- labels[[i]]
      tbl
    })
    data <- dplyr::bind_rows(Filter(Negate(is.null), data))
    if (!nrow(data)) {
      stop("No finite profile values were available to plot for the requested variable and date.", call. = FALSE)
    }
    p <- ggplot2::ggplot(
      data,
      ggplot2::aes(x = .data$value, y = .data$z, colour = .data$point_group, group = .data$point_group)
    ) +
      ggplot2::geom_path() +
      ggplot2::xlab(var_name) +
      ggplot2::ylab("metres above msl") +
      ggplot2::theme_minimal(base_size = 13) +
      ggplot2::labs(colour = if (length(unique(data$point_group)) > 1) "Point" else NULL)
    if (length(unique(data$point_group)) == 1) {
      p <- p + ggplot2::guides(colour = "none")
    }
    return(p)
  }
  plot_ereefs_profile(
    profile_result,
    var_name = var_name,
    target_date = profile_result$dates[[1]]
  )
}

ereefs_gui_point_labels <- function(points) {
  if ("label" %in% names(points) && any(nzchar(points$label %||% ""))) {
    return(points$label)
  }
  sprintf("(%.5f, %.5f)", points$latitude, points$longitude)
}

ereefs_gui_profile_script <- function(points, base_args) {
  if (nrow(points) == 1) {
    return(ereefs_gui_script(
      "vertical profile",
      "get_ereefs_profile",
      c(base_args, list(geolocation = c(points$latitude[[1]], points$longitude[[1]])))
    ))
  }
  paste(
    "points <- data.frame(",
    "  latitude = c(",
    paste(vapply(points$latitude, ereefs_gui_r_value, character(1)), collapse = ", "),
    "),",
    "  longitude = c(",
    paste(vapply(points$longitude, ereefs_gui_r_value, character(1)), collapse = ", "),
    "),",
    "  label = c(",
    paste(vapply(ereefs_gui_point_labels(points), ereefs_gui_r_value, character(1)), collapse = ", "),
    ")",
    ")",
    "",
    "profiles <- lapply(seq_len(nrow(points)), function(i) {",
    paste0(
      "  get_ereefs_profile(",
      "\n    var_names = ", ereefs_gui_r_value(base_args$var_names), ",",
      "\n    geolocation = c(points$latitude[[i]], points$longitude[[i]]),",
      "\n    start_date = ", ereefs_gui_r_value(base_args$start_date), ",",
      "\n    end_date = ", ereefs_gui_r_value(base_args$end_date), ",",
      "\n    input_file = ", ereefs_gui_r_value(base_args$input_file),
      "\n  )"
    ),
    "})",
    sep = "\n"
  )
}

ereefs_gui_script <- function(label, fun, args) {
  arg_lines <- vapply(names(args), function(name) {
    paste0("  ", name, " = ", ereefs_gui_r_value(args[[name]]))
  }, character(1))
  c(
    "# Reproduce eReefs GUI result from the command line.",
    paste0("# Workflow: ", label),
    "# Generated scripts do not include credentials or Shiny session state.",
    "library(ereefs)",
    "",
    paste0("result <- ", fun, "("),
    paste(arg_lines, collapse = ",\n"),
    ")",
    "",
    "result"
  ) %>%
    paste(collapse = "\n")
}

ereefs_gui_r_value <- function(x) {
  paste(capture.output(dput(x)), collapse = "")
}

ereefs_gui_range_text <- function(a, b) {
  paste(a, "to", b)
}

ereefs_gui_is_remote <- function(x) {
  grepl("^https?://|thredds|opendap|dodsC", as.character(x), ignore.case = TRUE)
}

ereefs_gui_clean_source_path <- function(x) {
  x <- trimws(as.character(x %||% ""))
  if (!nzchar(x)) {
    return(x)
  }
  x <- sub('^["\']', "", x)
  x <- sub('["\']$', "", x)
  trimws(x)
}

ereefs_gui_default_local_browser_dir <- function() {
  candidates <- c(
    getwd(),
    path.expand("~"),
    "C:/Users/brobson"
  )
  existing <- candidates[dir.exists(candidates)]
  normalizePath(existing[[1]], winslash = "\\", mustWork = FALSE)
}

ereefs_gui_local_browser_entries <- function(path) {
  path <- ereefs_gui_clean_source_path(path)
  if (!dir.exists(path)) {
    return(dplyr::tibble(label = character(), path = character(), type = character()))
  }
  entries <- tryCatch(
    list.files(path, full.names = TRUE, all.files = FALSE, no.. = TRUE),
    error = function(e) character()
  )
  if (!length(entries)) {
    return(dplyr::tibble(label = character(), path = character(), type = character()))
  }
  is_dir <- dir.exists(entries)
  is_file <- file.exists(entries) & !is_dir
  wanted_file <- is_file & grepl("\\.(nc|nc4|mnc|ncml)$", entries, ignore.case = TRUE)
  keep <- is_dir | wanted_file
  entries <- entries[keep]
  is_dir <- is_dir[keep]
  if (!length(entries)) {
    return(dplyr::tibble(label = character(), path = character(), type = character()))
  }
  dplyr::tibble(
    label = paste0(ifelse(is_dir, "[folder] ", "[file] "), basename(entries)),
    path = normalizePath(entries, winslash = "\\", mustWork = FALSE),
    type = ifelse(is_dir, "folder", "file")
  ) %>%
    dplyr::arrange(.data$type == "file", tolower(.data$label))
}

ereefs_gui_estimate_frames <- function(start_date, end_date, stride) {
  start_date <- as.Date(start_date)
  end_date <- as.Date(end_date)
  if (is.na(start_date) || is.na(end_date)) {
    return(NA_integer_)
  }
  by <- switch(stride, weekly = "1 week", monthly = "1 month", "1 day")
  length(seq(start_date, end_date, by = by))
}

ereefs_gui_current_map_plot <- function(map_result, input) {
  if (inherits(map_result, "ggplot")) {
    return(map_result)
  }
  plot_map(
    map_result,
    land_map = isTRUE(input$land_map),
    scale_col = input$map_palette,
    scale_lim = ereefs_gui_scale_lim(input$map_auto_scale, input$map_scale_min, input$map_scale_max),
    plot_style = input$map_style,
    smooth_pixels = input$smooth_pixels,
    box_bounds = ereefs_gui_display_bbox(input, map_result),
    label_towns = isTRUE(input$label_towns),
    gbr_poly = isTRUE(input$gbr_poly)
  )
}

ereefs_gui_display_bbox <- function(input, map_result) {
  bounds <- ereefs_gui_bbox(input)
  if (!all(is.na(bounds))) {
    return(bounds)
  }
  if (is.list(map_result) && "datapoly" %in% names(map_result)) {
    return(c(
      min(map_result$datapoly$x, na.rm = TRUE),
      max(map_result$datapoly$x, na.rm = TRUE),
      min(map_result$datapoly$y, na.rm = TRUE),
      max(map_result$datapoly$y, na.rm = TRUE)
    ))
  }
  bounds
}

ereefs_gui_animation_progress_text <- function(info) {
  frame_text <- if (is.finite(info$frame_index) && is.finite(info$frame_count)) {
    paste0(" (frame ", info$frame_index, " of ", info$frame_count, ")")
  } else {
    ""
  }
  date_text <- if (!is.null(info$date) && !all(is.na(info$date))) {
    paste0(" for ", as.character(info$date))
  } else {
    ""
  }
  paste0(info$message %||% info$stage, frame_text, date_text)
}

ereefs_gui_animation_progress_fraction <- function(info) {
  frame_index <- suppressWarnings(as.numeric(info$frame_index))
  frame_count <- suppressWarnings(as.numeric(info$frame_count))
  if (!is.finite(frame_index) || !is.finite(frame_count) || frame_count <= 0) {
    return(switch(info$stage, start = 0.05, scale = 0.58, assemble = 0.92, complete = 1, 0.1))
  }
  if (info$stage %in% c("load", "loaded")) {
    return(min(0.55, 0.05 + 0.45 * frame_index / frame_count))
  }
  if (info$stage %in% c("render", "saved")) {
    return(min(0.9, 0.58 + 0.3 * frame_index / frame_count))
  }
  switch(info$stage, scale = 0.57, assemble = 0.93, complete = 1, 0.1)
}

ereefs_gui_inspection_progress_fraction <- function(info) {
  switch(
    info$stage %||% "",
    start = 0.12,
    files = 0.2,
    representative = 0.3,
    dimensions = 0.45,
    variables = 0.65,
    spatial = 0.82,
    time = 0.9,
    complete = 1,
    0.15
  )
}

ereefs_gui_extraction_progress_fraction <- function(info, start = 0.2, span = 0.75) {
  stage <- info$stage %||% ""
  if (identical(stage, "complete")) {
    return(start + span)
  }
  if (identical(stage, "start")) {
    return(start + 0.05 * span)
  }
  file_index <- suppressWarnings(as.numeric(info$file_index))
  file_count <- suppressWarnings(as.numeric(info$file_count))
  variable_index <- suppressWarnings(as.numeric(info$variable_index))
  variable_count <- suppressWarnings(as.numeric(info$variable_count))
  point_index <- suppressWarnings(as.numeric(info$point_index))
  point_count <- suppressWarnings(as.numeric(info$point_count))
  if (is.finite(point_index) && is.finite(point_count) && point_count > 0) {
    return(start + span * min(0.98, point_index / point_count))
  }
  if (is.finite(file_index) && is.finite(file_count) && file_count > 0) {
    if (is.finite(variable_index) && is.finite(variable_count) && variable_count > 0) {
      progress <- ((file_index - 1) + variable_index / variable_count) / file_count
      return(start + span * min(0.98, progress))
    }
    return(start + span * min(0.98, file_index / file_count))
  }
  start + 0.1 * span
}

ereefs_gui_progress_detail <- function(info) {
  stage_message <- info$message %||% "Working"
  parts <- c(stage_message)
  if (!is.null(info$file) && !is.na(info$file) && nzchar(info$file)) {
    parts <- c(parts, paste("File:", basename(info$file)))
  }
  if (!is.null(info$variable) && !is.na(info$variable) && nzchar(info$variable)) {
    parts <- c(parts, paste("Variable:", info$variable))
  }
  if (is.finite(suppressWarnings(as.numeric(info$point_index))) && is.finite(suppressWarnings(as.numeric(info$point_count)))) {
    parts <- c(parts, paste("Point:", info$point_index, "of", info$point_count))
  }
  if (is.finite(suppressWarnings(as.numeric(info$time_count)))) {
    parts <- c(parts, paste("Matched times:", info$time_count))
  }
  if (!is.null(info$time_start) && !all(is.na(info$time_start))) {
    if (!is.null(info$time_end) && !all(is.na(info$time_end)) && !identical(info$time_start, info$time_end)) {
      parts <- c(parts, paste("Date range:", as.character(as.Date(info$time_start)), "to", as.character(as.Date(info$time_end))))
    } else {
      parts <- c(parts, paste("Date:", as.character(as.Date(info$time_start))))
    }
  }
  if (is.finite(suppressWarnings(as.numeric(info$file_count)))) {
    parts <- c(parts, paste("Files:", info$file_count))
  }
  if (is.finite(suppressWarnings(as.numeric(info$variable_count)))) {
    parts <- c(parts, paste("Variables:", info$variable_count))
  }
  paste(parts, collapse = "\n")
}

ereefs_gui_register_animation_media <- function(session, state, animation_file) {
  if (is.null(animation_file) || length(animation_file) != 1L || is.na(animation_file) || !file.exists(animation_file)) {
    state$animation_media_src <- NULL
    state$animation_media_type <- NULL
    return(invisible(NULL))
  }
  prefix <- paste0("ereefs-animation-", as.integer(Sys.time()))
  shiny::addResourcePath(prefix, normalizePath(dirname(animation_file), winslash = "/", mustWork = TRUE))
  state$animation_media_src <- paste0(prefix, "/", basename(animation_file))
  state$animation_media_type <- tolower(tools::file_ext(animation_file))
  invisible(NULL)
}

ereefs_gui_register_animation_frame <- function(session, state, frame_file) {
  if (is.null(frame_file) || length(frame_file) != 1L || is.na(frame_file) || !file.exists(frame_file)) {
    state$animation_frame_src <- NULL
    return(invisible(NULL))
  }
  prefix <- paste0("ereefs-frame-", as.integer(Sys.time()))
  shiny::addResourcePath(prefix, normalizePath(dirname(frame_file), winslash = "/", mustWork = TRUE))
  state$animation_frame_src <- paste0(prefix, "/", basename(frame_file), "?t=", as.integer(file.info(frame_file)$mtime))
  invisible(NULL)
}

ereefs_gui_flush <- function(session) {
  tryCatch({
    if (!is.null(session$flushReact) && is.function(session$flushReact)) {
      session$flushReact()
    } else {
      get("flushReact", envir = asNamespace("shiny"))()
    }
  }, error = function(e) NULL)
  invisible(NULL)
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

# metadata:
# - gpt_version: GPT-5 Codex
# - time: 12:05
# - date: 2026-07-05
# - prompt_used: "Implement the staged eReefs GUI plan with an attractive Shiny interface, progress-preserving notes, source selection, map preview, animation, point extraction, transect hooks, exports, and conservative slow-request safeguards."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 12:11
# - date: 2026-07-05
# - prompt_used: "Polish the staged GUI implementation so the launch-time source is visible in the UI and saved point time-series PNGs explicitly print ggplot output."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 12:29
# - date: 2026-07-05
# - prompt_used: "Make inspected real data variables the default GUI selections while keeping special map variables available."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 11:58
# - date: 2026-07-05
# - prompt_used: "Add the eReefs logo to the GUI launch header and make NCI catalog browsing fall back to a curated list when live HTTPS catalog loading fails."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 11:59
# - date: 2026-07-05
# - prompt_used: "Remove live Google-font handling from the GUI theme to avoid bslib cache and network failures while keeping an attractive local-font design."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 12:03
# - date: 2026-07-05
# - prompt_used: "Remove dynamic bslib theme compilation from the GUI so the page renders in restricted Windows cache environments using stable custom CSS instead."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 12:11
# - date: 2026-07-05
# - prompt_used: "Replace bslib cards and layouts with pure Shiny layout helpers to avoid bslib runtime Sass dependencies while preserving the GUI visual design."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 12:33
# - date: 2026-07-05
# - prompt_used: "Combine manual URL and local filepath source entry, add a local file browser, fix NCI catalog attribute parsing, hard-code naming help text, and rename the GUI to eReefs R GUI."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 16:13
# - date: 2026-07-05
# - prompt_used: "Fix GUI catalog inspection, local-file animation rendering, first-point profile display, remove slow-confirmation checkboxes, and add mouse point selection."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 16:13
# - date: 2026-07-05
# - prompt_used: "Match GUI slow-operation warning text to the requested wording."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 16:26
# - date: 2026-07-05
# - prompt_used: "Finalize GUI bug fixes after verification: catalog inspect, local animation, first-point profile, warning text, and mouse-selected points."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 16:40
# - date: 2026-07-05
# - prompt_used: "Make GUI NCI catalog inspection fast by default, with optional slower representative remote metadata and fallback common variable choices."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 16:43
# - date: 2026-07-05
# - prompt_used: "Show the selected NetCDF file, OPeNDAP URL, or THREDDS catalog address on every GUI workflow tab and clarify metadata-only inspection."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 16:50
# - date: 2026-07-05
# - prompt_used: "Fix malformed NCI catalog URLs by resolving absolute, root-relative, and relative THREDDS catalog hrefs correctly."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 17:05
# - date: 2026-07-05
# - prompt_used: "Avoid defaulting NCI hydrodynamic catalog maps to true_colour by separating map-only special variables from extractable data variables."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 17:20
# - date: 2026-07-05
# - prompt_used: "Move mouse point selection from a small schematic picker onto the full-sized map preview tab and overlay selected points on the preview map."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 17:30
# - date: 2026-07-05
# - prompt_used: "Rename code disclosure labels to Show R code, visually indicate they are clickable, and make Use bounding box unselected by default."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 17:40
# - date: 2026-07-05
# - prompt_used: "Show a session-scoped warning modal when users select the NCI catalog source option, with Cancel, OK, and Don't show again this session actions."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 18:05
# - date: 2026-07-05
# - prompt_used: "Recalculate/reset map spatial domain when a new source is inspected and show map-clicked points in the Points tab extraction list."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 18:25
# - date: 2026-07-05
# - prompt_used: "Allow users to edit labels for clicked map points and keep edited labels synced with the map overlay and Points tab extraction list."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 18:35
# - date: 2026-07-05
# - prompt_used: "Constrain default GUI map time step and start/end date controls to the inspected dataset date bounds."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 18:45
# - date: 2026-07-05
# - prompt_used: "Fix map-click point coordinates by converting normalized panel clicks back to displayed longitude and latitude bounds."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 19:00
# - date: 2026-07-05
# - prompt_used: "Allow users to change the current map colour scale and display style from stored preview data without re-downloading the NetCDF subset."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 19:30
# - date: 2026-07-05
# - prompt_used: "Allow large local NetCDF files to be selected by path without Shiny upload, and clean quoted pasted Windows paths."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 19:45
# - date: 2026-07-05
# - prompt_used: "Use the Tcl/Tk local file chooser for large NetCDF path selection because the native Windows choose.files dialog displayed icons without filename text."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 19:45
# - date: 2026-07-05
# - prompt_used: "Clarify GUI wording for the local file-path chooser after switching away from the native upload-limited dialog."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 10:20
# - date: 2026-07-06
# - prompt_used: "Replace the problematic native/Tcl local file dialog with an in-app Shiny local file browser that displays directory and NetCDF filenames as HTML text."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 10:25
# - date: 2026-07-06
# - prompt_used: "Remove the unused native/Tcl local file dialog helper so only the in-app text-rendered local browser remains."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 10:30
# - date: 2026-07-06
# - prompt_used: "Fix Shiny local file browser selectInput by setting selectize = FALSE when using the size argument."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 11:24
# - date: 2026-07-06
# - prompt_used: "Make polygon the default GUI map style."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 11:28
# - date: 2026-07-06
# - prompt_used: "Automatically load the NCI catalog menu when the user switches to browse NCI catalogs, and remove the manual load button."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 11:37
# - date: 2026-07-06
# - prompt_used: "Remove the recurse-into-child-catalogs checkbox and make catalog recursion always on."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 11:40
# - date: 2026-07-06
# - prompt_used: "Hide layer controls, first-point profile, and the Transect tab when inspection shows the selected dataset has only one layer."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 11:40
# - date: 2026-07-06
# - prompt_used: "Keep hidden single-layer GUI controls safe by defaulting missing layer inputs to surface in map and point extraction calls."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 11:47
# - date: 2026-07-06
# - prompt_used: "Clear stale inspection state and reset downstream GUI controls when the selected source changes so variable menus always match the currently selected catalog after re-inspection."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 11:51
# - date: 2026-07-06
# - prompt_used: "Clarify that the catalog-file limit only affects the inspection-page listing, remove the variable display limit, and rename the Inspect dataset button to Choose dataset."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 12:05
# - date: 2026-07-06
# - prompt_used: "Show the selected variable's layer status next to the Layer field in the Map and Points tabs."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 12:28
# - date: 2026-07-06
# - prompt_used: "Always load real variable metadata for catalog inspection so NCI catalog variables appear on the inspection page and in downstream variable menus, while keeping the slower full remote grid/time metadata as an optional extra."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 12:30
# - date: 2026-07-06
# - prompt_used: "Replace the dynamic uiOutput tab placeholder with proper insert/remove tab handling so the conditional Transect tab does not trigger navigation-container warnings."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 12:54
# - date: 2026-07-06
# - prompt_used: "Fix GUI point time-series plotting so it plots the selected variable rather than latitude and draws separate lines for multiple selected points."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 12:59
# - date: 2026-07-06
# - prompt_used: "Disable unreliable map-click point selection for now and keep point entry manual or CSV-based until click coordinates can be trusted."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 13:11
# - date: 2026-07-06
# - prompt_used: "If the user changes the colour scale limits but does not change the date, layer or variable selected, the tool should re-render the map but not re-fetch the data."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 13:24
# - date: 2026-07-06
# - prompt_used: "If the user enters a new value in colour min or colour max, untick auto colour range; if the user ticks auto colour range, show the current auto-derived colour min and colour max after re-rendering the map."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 13:28
# - date: 2026-07-06
# - prompt_used: "Point time-series extraction is failing with the error: 'from' must be a finite number'"
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 13:55
# - date: 2026-07-06
# - prompt_used: "Apply the same logic to the bounding box as we applied to the colour scale. If the user manually edits the Lat and Lon mins and maxes, 'Apply bounding box' should be selected. If 'Apply bounding box' is unselected, the mins and maxes should be populated with the spatial domain limits."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 14:07
# - date: 2026-07-06
# - prompt_used: "Remove the 'Alternative browser upload for small local NetCDF or mnc files' option. The files we want will always be relatively large."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 14:12
# - date: 2026-07-06
# - prompt_used: "Now I am getting 'Map preview failed: missing value where TRUE/FALSE needed'"
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 15:49
# - date: 2026-07-06
# - prompt_used: "If the variable is changed in one tab, it should be changed to match in other tabs"
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 15:55
# - date: 2026-07-06
# - prompt_used: "Fix the one-point time-series path so point parsing always returns latitude/longitude columns and does not crash the warning text when only one point is supplied."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 15:58
# - date: 2026-07-06
# - prompt_used: "Fix the empty GUI point-table fallback so it builds zero-length numeric and character columns correctly."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 16:26
# - date: 2026-07-06
# - prompt_used: "Fix the shared GUI point date-range helper and profile extraction path so transiently missing dates do not break profile runs."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 16:26
# - date: 2026-07-06
# - prompt_used: "Next to the start date and end date fields and the date field on the Transect tab, show the selected file/catalog date range, and snap out-of-range entered dates to the nearest valid date."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 16:40
# - date: 2026-07-06
# - prompt_used: "Show animation frames as they are generated, and replace the last frame with the final animation in the same preview area once assembly is complete."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 16:55
# - date: 2026-07-06
# - prompt_used: "Add the capability to plot profiles for multiple selected points in the GUI rather than only the first point."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 17:18
# - date: 2026-07-06
# - prompt_used: "During extractions and inspections, display richer progress detail such as stage, current file, variable, date span, file count, variable count, and point index without adding extra I/O."
