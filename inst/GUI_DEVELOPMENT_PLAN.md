# eReefs GUI development plan

This document describes a staged, vertical-slice approach for adding an in-R GUI
to the `ereefs` package. The goal is to make common exploration and extraction
workflows approachable without hiding the underlying R functions or encouraging
large accidental OPeNDAP downloads.

## Technology choice

Use Shiny for the first GUI because it is the standard R-native way to build
interactive apps directly from R sessions. The implemented vertical slice uses
standard Shiny layout components plus custom CSS rather than `bslib`, because
`bslib`'s runtime Sass/cache path was brittle in the restricted Windows test
environment. Follow-up testing under normal Windows `R 4.4.0` suggests the
problem is specific to the restricted Codex sandbox rather than to `bslib`
itself, but the current GUI still keeps the pure-Shiny approach because it is
working reliably. Keep Shiny in `Suggests`, not `Imports`, so command-line use
of `ereefs` stays lightweight.

Do not introduce `golem` yet. It is a good candidate if the GUI becomes a large,
production-style app with modules, configuration, deployment targets, and formal
JavaScript/CSS assets, but it would be premature for the first vertical slices.

## Development principles

- Build one working end-to-end workflow at a time.
- Prefer metadata-only or spatially bounded actions by default.
- Avoid reading large remote arrays unless the user has explicitly requested an
  extraction or plot.
- Keep every GUI action traceable to an exported package function.
- Preserve scriptability: the GUI should help users discover arguments and then
  show the equivalent R call where practical.
- Request approval before starting each new stage.

## Cross-cutting outputs and reproducibility

Every workflow that produces a result should make it easy to save or reproduce
that result without trapping the user inside the GUI.

- Add save buttons for generated maps, vertical slices, transects, profiles,
  animations, and other figures. Use PNG as the default image format, with room
  to add PDF/SVG options later where the underlying plotting code supports them.
- Add CSV download buttons for tabular outputs, especially point-extraction time
  series and profile results.
- Add an "R script" download button for each successful plot or extraction. The
  script should recreate the same result from the command line using exported
  `ereefs` functions and the selected source, variable, date range, bounding
  box, colour settings, and output path.
- Do not embed credentials, local secrets, or temporary Shiny session state in
  generated scripts. If an authenticated or private source is ever supported,
  the script should include a clear placeholder and comment rather than a token.

## Stage 1: Dataset inspection shell

Status: implemented.

User story: A user can run `launch_ereefs_gui()`, paste a NetCDF path, OPeNDAP
URL, THREDDS catalog, or legacy shortcut, and inspect available files,
variables, dimensions, and basic spatial/temporal metadata.

Scope:

- Add `launch_ereefs_gui()`.
- Keep `shiny` optional.
- Use `inspect_ereefs_data()` as the first backend.
- Show metadata tables only; do not extract model arrays.

Approval gate:

- Confirm that the launch workflow and first screen are useful enough before
  adding plotting or extraction controls.

## Stage 1b: Dataset source chooser

Status: implemented.

User story: A user can choose where the dataset comes from without needing to
remember a long THREDDS URL. They can either select a local NetCDF file, paste a
NetCDF/OPeNDAP/THREDDS URL, or browse the current NCI eReefs catalog menu at
<https://thredds.nci.org.au/thredds/catalog/catalogs/fx3/catalog.html> and
choose one of the available eReefs catalogs.

Scope:

- Add a source-mode selector with options for typed path/URL/shortcut, local
  file browser, and NCI catalog browser.
- For local files, support both typing a path in the same field used for URLs
  and browsing for a `.nc`, `.nc4`, `.mnc`, or `.ncml` file.
- For typed sources, retain the current text entry path and validate whether the
  value looks like a local file, NetCDF/OPeNDAP endpoint, THREDDS catalog, or
  package shortcut before running metadata inspection.
- For NCI catalogs, fetch and parse the `fx3` catalog menu, list catalog names
  and URLs, and pass the selected catalog URL into `inspect_ereefs_data()`.
- Mark superseded/deprecated catalogs visibly when that wording is available in
  the catalog label.
- Add a "Naming help" button or modal beside the NCI catalog chooser. This
  should link to the eReefs Marine Model Results page and summarize that dataset
  IDs encode model grid, hydrodynamic version, forcing data, biogeochemical or
  tracer configuration, catchment scenario, and hindcast/near-real-time mode.
  The help text should point users to
  <https://www.ereefs.org.au/outputs/open-access-datasets/marine-model-results>
  for the current naming and dataset guidance.
- Cache the catalog-menu lookup for the active Shiny session so the app does not
  repeatedly query THREDDS while users adjust controls.
- If live NCI catalog menu loading fails, present a curated fallback list of
  current high-value eReefs catalogs rather than leaving the selector empty.

Approval gate:

- Confirm that the source chooser makes NCI catalog selection clearer before
  adding plot previews.

## Stage 2: Map preview

Status: implemented.

User story: A user can choose a variable/date/layer from an inspected source and
produce a `map_ereefs()` preview before deciding whether to extract more data or
make an animation.

Scope:

- Populate variable/date choices from the Stage 1 inspection object.
- Default to the first available time step in the selected source or filtered
  date range.
- Use the full model spatial domain as the default map extent. Provide optional
  bounding-box controls so users can zoom to a smaller region before plotting.
- Let users edit plotting controls such as colour palette, colour scale limits,
  transformation or clipping behaviour where supported, smooth/raster versus
  polygon rendering, labels, and map title.
- Surface the equivalent `map_ereefs()` call and provide buttons to save the map
  image and a matching command-line R script.
- Include progress and friendly error handling for unavailable variables or
  dates.

## Stage 3: Animation from map preview

Status: implemented.

User story: After seeing a map preview, a user can reuse those settings to create
a short animation over a selected time window.

Scope:

- Start from the current map preview settings, including source, variable,
  bounding box, palette, colour range, labels, and plot style.
- Let users choose animation start/end times, temporal step, output format,
  frame rate or delay, and output directory.
- Show the estimated frame count and warn before large spatial or temporal
  requests, especially for remote OPeNDAP sources.
- Preserve a fixed colour scale across all frames by default, with the preview
  range used as the initial range and editable before rendering.
- Use `map_ereefs_movie()` and provide buttons to save the animation, retain or
  clean up generated frame images as selected, and download the equivalent R
  script.

Approval gate:

- Confirm that preview-to-animation feels safe and understandable before adding
  point extraction.

## Stage 4: Point extraction

Status: implemented for multi-point time series and multi-point profiles, with
map-click point selection currently disabled pending a reliable coordinate fix.

User story: A user can click or type a point and extract a surface time series
or vertical profile. They can also extract the same variable at multiple points.

Scope:

- Support manual entry of one or more latitude/longitude points, map-based point
  selection where practical, and CSV upload of point locations.
- For uploaded CSV files, validate required coordinate columns and accept an
  optional label/name column for plotting and exported results.
- Show nearest-cell feedback for each selected point before extraction so users
  can catch coordinate or grid-matching mistakes early.
- Add `get_ereefs_ts()` and `get_ereefs_profile()` workflows.
- Preview time series and profile plots, including sensible grouping or faceting
  when multiple points are selected.
- Keep date ranges explicit. If the number of points, date range, depth range, or
  estimated remote reads is large, display a clear warning near the action
  button; extraction still only runs after the user explicitly clicks the
  relevant action.
- Provide buttons to save generated figures, export time series/profile tables
  as CSV, and download the equivalent command-line R script.

## Stage 5: Transects and vertical slices

Status: implemented as a typed-coordinate workflow, with the Transect tab
suppressed when the selected source exposes only one layer.

User story: A user can define a two-point or multi-point transect and generate a
surface transect or vertical slice.

Scope:

- Support typed coordinate tables first; map-click editing can come later.
- Use `get_ereefs_slice()` and `plot_ereefs_slice()`.
- Show matched cell count before extraction.
- Keep extraction efficient by reusing the bulk-read code paths.
- Provide figure, CSV where appropriate, and R-script download buttons.

## Stage 6: Packaging, tests, and polish

Status: partly implemented.

User story: The GUI is documented, tested, and safe to ship with the package.

Scope:

- Add a vignette or README section for `launch_ereefs_gui()`.
- Add Shiny smoke tests where practical.
- Consider modularising with Shiny modules or `golem` if the app has grown
  enough to justify that structure.
- Add accessibility and usability checks.
- Keep the lightweight progress reporting model: show stage, file, variable,
  point index, and matched date span without introducing extra remote reads.

<!--
metadata:
- gpt_version: GPT-5 Codex
- time: 14:45
- date: 2026-07-05
- prompt_used: "Plan a staged in-R GUI for the eReefs package, implement stage 1, and wait for approval before continuing."
-->
<!--
metadata:
- gpt_version: GPT-5 Codex
- time: 11:20
- date: 2026-07-05
- prompt_used: "Add a planned GUI source chooser for local NetCDF files, manual URLs, NCI THREDDS catalog browsing, and naming-convention help."
-->
<!--
metadata:
- gpt_version: GPT-5 Codex
- time: 11:29
- date: 2026-07-05
- prompt_used: "Update the GUI plan so map previews default to full extent and first timestep, animations follow map preview, point extraction supports multiple/CSV points with slow-run confirmation, and outputs can be saved as images, CSV, and reproducible R scripts."
-->
<!--
metadata:
- gpt_version: GPT-5 Codex
- time: 12:36
- date: 2026-07-05
- prompt_used: "Update the GUI plan to reflect the pure-Shiny implementation, combined typed path/URL source entry, local file browser, and NCI catalog fallback."
-->
<!--
metadata:
- gpt_version: GPT-5 Codex
- time: 16:13
- date: 2026-07-05
- prompt_used: "Update GUI plan for warning-only slow-operation messaging and mouse point selection behaviour."
-->
<!--
metadata:
- gpt_version: GPT-5 Codex
- time: 16:26
- date: 2026-07-05
- prompt_used: "Finalize GUI plan metadata after implementing and verifying warning-only slow-operation messaging and mouse-selected points."
-->
<!--
metadata:
- gpt_version: GPT-5 Codex
- time: 09:41
- date: 2026-07-07
- prompt_used: "Update the GUI development plan to reflect the current implementation state, the refined bslib interpretation, multi-point profiles, disabled map-click selection, conditional transect visibility, and lightweight progress reporting."
-->
