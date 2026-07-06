# eReefs GUI progress notes

These notes track implementation progress for the staged in-R GUI so work can be
continued in another chat without having to reconstruct context from memory.

## 2026-07-05

### Security and robustness check

- No new package dependency was added for this implementation pass.
- The GUI remains metadata-first: inspecting a dataset calls
  `inspect_ereefs_data()` and does not intentionally read full model arrays.
- Potentially slow remote reads are attached to explicit action buttons for maps,
  animations, point extraction, profiles, and slices.
- Generated R scripts use selected function arguments but do not embed
  credentials, tokens, or Shiny session state.
- Remote source strings are treated as data-source URLs for package functions,
  not as code to execute.

### Stage 1 and Stage 1b: Source selection and inspection

Status: implemented in `R/ereefs_gui.R`.

- Added a visually richer Shiny application shell with workflow tabs.
- Added source modes for typed local path/URL/catalog/shortcut, local NetCDF or
  mnc file browsing, and NCI catalog browsing.
- Added an NCI catalog loader that reads the current `fx3` THREDDS catalog menu.
- Added a naming-help modal that points to the eReefs Marine Model Results page.
- Retained metadata-only inspection using `inspect_ereefs_data()`.

### Stage 2: Map preview

Status: implemented as a first vertical slice.

- Added variable, date, layer, bounding-box, palette, colour-range, map-style,
  town-label, land-map, and GBR-polygon controls.
- Defaults are updated from the inspection summary, including first available
  date and full model spatial extent.
- Added map preview with `map_ereefs()`.
- Added buttons to save map PNGs and command-line R scripts.

### Stage 3: Animation after map preview

Status: implemented as a first vertical slice.

- Animation controls reuse the current map settings.
- Added start/end date, stride, output format, frame rate, and keep-frame
  controls.
- Added frame-count estimates and confirmation for remote or larger requests.
- Added `map_ereefs_movie()` orchestration and an R-script export.

### Stage 4: Point extraction

Status: implemented for multi-point time series and conservative first-point
profiles.

- Added manual multi-point entry in `latitude,longitude,label` format.
- Added CSV upload for files with `latitude` and `longitude` columns and an
  optional `label` column.
- Added estimated point-day warnings and confirmation before large/remote
  extraction.
- Added `get_ereefs_ts()` extraction for multiple points, with plot, CSV, PNG,
  and R-script exports.
- Added `get_ereefs_profile()` for the first selected point. This is deliberately
  conservative because the backend profile function currently handles one
  location per call.

### Stage 5: Transects and vertical slices

Status: implemented as a typed-coordinate hook.

- Added a coordinate text box for multi-point paths.
- Added `get_ereefs_slice()` and `plot_ereefs_slice()` calls.
- Added PNG and R-script exports.

### Remaining checks

- Run syntax and package-load checks after this implementation pass.
- Smoke-test the app object creation with local demo NetCDF data.
- If time permits, manually launch the GUI in a browser and test the first map
  preview against `notebooks/demo_data/regular_demo_2020-01.nc`.
- Consider adding formal Shiny tests once the user has approved the interaction
  design.

### 2026-07-05 follow-up

- R 4.4.0 can create the Shiny app object and inspect
  `notebooks/demo_data/regular_demo_2020-01.nc`.
- R 4.3.2 did not have `bslib`/`shiny` available in the active library path.
- R and PowerShell HTTPS requests to NCI/eReefs URLs failed on Windows with
  schannel credential/TLS errors during this session. A WSL download succeeded
  for the public eReefs logo asset, and subsequent work should use PowerShell
  where possible.
- Added `inst/app/www/ereefs_logo.png` and wired the GUI hero header to serve it
  locally through `shiny::addResourcePath()`.
- Added a curated fallback list of major NCI eReefs catalogs so the catalog
  chooser remains useful if live THREDDS menu loading fails.
- Replaced the runtime `bslib` UI with pure Shiny layout/card helpers plus
  custom CSS because bslib's Sass/theme dependency path failed in this
  restricted Windows cache environment.
- Live Shiny response check succeeded on `http://127.0.0.1:8766` with HTTP 200.
  The rendered HTML included the eReefs hero, local logo reference, staged tabs,
  and save buttons.
- Headless Edge/Chrome screenshot attempts from PowerShell returned without
  saving an image despite the Shiny HTTP 200 response. Visual inspection can be
  done by launching the app from PowerShell and opening the URL in a normal
  browser.
- Renamed the GUI-facing title from "eReefs Explorer" to "eReefs R GUI" to avoid
  confusion with the existing eReefs Explorer product.
- Combined the manual URL and local filepath modes into one typed source field.
  The separate browser mode now accepts `.nc`, `.nc4`, `.mnc`, and `.ncml`
  files.
- Replaced the NCI naming-help link with hard-coded guidance in the modal, based
  on the eReefs Marine Model Results naming information.
- Fixed the NCI catalog parser so missing `xlink:*` attributes do not throw
  "subscript out of bounds"; the helper now safely falls back to plain
  attributes or the curated catalog list.
- Fixed catalog inspection when catalog file tables lack parsed `file_date`
  columns by normalising file metadata before temporal summaries.
- Hardened live NCI catalogue inspection against xarray/PyDAP dimension mappings
  that do not convert to ordinary R named vectors.
- Fixed GUI animation rendering from local files by removing the map-preview
  `target_date` argument before calling `map_ereefs_movie()`.
- Removed slow-operation acknowledgement checkboxes and kept visible slow-run
  warnings beside the animation and point-extraction actions.
- Added mouse point selection in the Points stage. Clicked points are combined
  with typed or CSV-uploaded points and can be cleared before extraction.
- Updated profile plotting to explicitly print the generated profile plot and
  use the first returned profile date, so clicking "Extract first-point profile"
  displays the plot when extraction succeeds.
- Changed GUI NCI inspection to use fast catalogue-only metadata by default.
  Opening the first remote OPeNDAP file for full variable/grid metadata is now
  optional because it can stall on large NCI catalogues.
- Added a persistent Current source ribbon above the workflow tabs so the active
  local NetCDF path, OPeNDAP URL, or THREDDS catalogue address remains visible
  while mapping, animating, extracting points, or building transects.
- Confirmed the fast live NCI catalogue inspection completed from the THREDDS
  catalogue listing without triggering PyDAP or reading a whole NetCDF file.
- Fixed NCI catalogue menu URL resolution so root-relative THREDDS links such
  as `/thredds/catalog/fx3/.../catalog.xml` are not incorrectly pasted under
  `/thredds/catalog/catalogs/fx3/`. The BGC catalogue URL that previously
  produced a 404 now inspects successfully through the fast catalogue path.
- Split GUI variable choices so map-only optical products (`true_colour` and
  `plume`) are only offered when inspected metadata shows the required
  reflectance bands. Fast catalogue-only NCI hydrodynamic inspections now
  default to `temp` rather than `true_colour`, avoiding missing `R_645` errors.
- Verified a small-bounding-box `temp` map from the NCI hydrodynamic catalogue
  renders successfully; PyDAP and reconstructed `z_grid` warnings are expected
  for remote simple-file subset reads.
- Moved mouse point selection from the small schematic picker to the full-sized
  Map preview. Selected points are overlaid on the preview map, included in
  saved map PNGs, and reviewed/extracted from the Points tab.

### 2026-07-06 follow-up

- Re-tested the earlier `bslib` concern explicitly under `R 4.4.0`.
- The current Windows environment can load `shiny`, `bslib`, `sass`, and
  `httpuv` successfully under `R 4.4.0`.
- A minimal `bslib::page_fluid()` plus `bslib::card()` smoke test rendered and
  saved successfully when run unsandboxed with
  `C:/Program Files/R/R-4.4.0/bin/Rscript.exe`.
- The same smoke test failed inside the Codex restricted filesystem sandbox
  with `EPERM` while `htmltools::renderTags()` was trying to stat
  `C:/Users/brobson`.
- Current interpretation: the older `bslib` issue appears to be tied to the
  restricted/sandboxed Windows execution context rather than to a general
  `bslib` incompatibility in normal `R 4.4.0` sessions.
- Keep the pure-Shiny GUI for now unless there is a specific reason to restore
  runtime `bslib`; if revisiting that decision, test both a normal Windows R
  session and the restricted Codex environment because they behave differently.

### 2026-07-06 implementation sweep

- The NCI catalog browser now auto-loads the selected catalog instead of
  requiring a separate "Load catalog" button.
- Child-catalog recursion is now always on for NCI browsing, and the old
  recurse toggle was removed from the GUI.
- The "maximum catalogs to show" control now only limits what is shown on the
  inspection page; later map/extraction tabs can still access the full catalog
  set.
- The "maximum variables to show" limit was removed so the GUI can show all
  inspected variables.
- The inspection action label was changed from "Inspect dataset" to
  "choose dataset".
- Variable selections are now synchronized across tabs so changing a variable in
  one workflow updates the others where that variable is available.
- For sources with only one layer, tabs 2-5 now hide the layer field and label,
  remove the old first-point-profile option, and suppress the Transect tab.
- The layer field now shows the number of layers available for the selected
  variable.
- The map default was changed from smooth rendering to polygon rendering.
- Polygon maps now use matching line and fill colours for each polygon.
- Map recolouring and manual colour-range edits now re-render from cached data
  when the source, date, layer, and variable have not changed.
- Manual colour-min and colour-max edits automatically disable auto-range; when
  auto-range is re-enabled, the GUI writes the derived min/max values back into
  the visible controls after rendering.
- Bounding-box controls now follow the same manual-versus-auto pattern as the
  colour controls: manual edits enable the bounding box, and disabling the box
  repopulates the full inspected spatial extent.
- The old alternative browser upload option for small local files was removed.
- Mouse-based point selection from the map preview was found to be unreliable
  and has been disabled for now rather than leaving misleading clicked points in
  the workflow.
- Point, profile, slice, map, and inspection workflows were hardened against
  transient `NULL`, missing, or out-of-range GUI inputs that had previously led
  to "missing value where TRUE/FALSE needed" and related errors.
- Start/end/date controls now show the selected source date range on the
  Points and Transect workflows, and out-of-range entered dates snap to the
  nearest valid boundary.
- Animation preview now shows frames as they are generated and then replaces the
  latest frame with the final assembled animation in the same preview area.
- Metadata inspection for catalog-backed datasets is faster because it now
  assumes all NetCDF files in a given NCI catalog share the same variables,
  grid, and layer structure, and it uses a representative file accordingly.
- Point time-series extraction was fixed so single-point and multi-point runs
  both return the requested variable values rather than plotting latitude by
  mistake.
- Profile extraction was fixed for missing-date, missing-layer, and
  no-finite-value edge cases, and the GUI can now plot profiles for multiple
  selected points instead of only the first one.
- Inspection and extraction progress reporting is now richer without adding
  extra I/O. The GUI shows the current stage plus lightweight context such as
  source file name, variable, point index, matched date span, file count, and
  variable count.

<!--
metadata:
- gpt_version: GPT-5 Codex
- time: 12:00
- date: 2026-07-05
- prompt_used: "Record GUI validation notes, Windows HTTPS limitations, local eReefs logo asset handling, and NCI catalog fallback behaviour."
-->
<!--
metadata:
- gpt_version: GPT-5 Codex
- time: 12:22
- date: 2026-07-05
- prompt_used: "Record the pure-Shiny GUI pivot, successful live HTML response check, and unresolved headless screenshot limitation."
-->
<!--
metadata:
- gpt_version: GPT-5 Codex
- time: 12:37
- date: 2026-07-05
- prompt_used: "Record source selector simplification, local file browser support, NCI catalog parser fix, hard-coded naming help, and eReefs R GUI rename."
-->
<!--
metadata:
- gpt_version: GPT-5 Codex
- time: 16:13
- date: 2026-07-05
- prompt_used: "Record fixes for GUI catalog inspection, animation argument handling, first-point profile display, warning-only slow-operation messaging, and mouse point selection."
-->
<!--
metadata:
- gpt_version: GPT-5 Codex
- time: 16:13
- date: 2026-07-05
- prompt_used: "Record live NCI inspection hardening for PyDAP-backed xarray dimension mappings."
-->
<!--
metadata:
- gpt_version: GPT-5 Codex
- time: 16:26
- date: 2026-07-05
- prompt_used: "Finalize progress notes after verified GUI bug fixes and live NCI inspection smoke test."
-->
<!--
metadata:
- gpt_version: GPT-5 Codex
- time: 16:44
- date: 2026-07-05
- prompt_used: "Record fast NCI catalogue-only inspection, persistent current-source display, and no whole-file loading behaviour."
-->
<!--
metadata:
- gpt_version: GPT-5 Codex
- time: 16:55
- date: 2026-07-05
- prompt_used: "Record fix for malformed NCI catalogue menu URLs that caused 404s during map preview."
-->
<!--
metadata:
- gpt_version: GPT-5 Codex
- time: 17:05
- date: 2026-07-05
- prompt_used: "Record variable-choice fix so NCI hydrodynamic catalogues default to temp rather than true_colour and avoid missing R_645 errors."
-->
<!--
metadata:
- gpt_version: GPT-5 Codex
- time: 17:20
- date: 2026-07-05
- prompt_used: "Record moving point selection onto the full-sized map preview and keeping Points tab for review and extraction."
-->
<!--
metadata:
- gpt_version: GPT-5 Codex
- time: 09:40
- date: 2026-07-07
- prompt_used: "Update progress notes to capture the July 6 GUI implementation sweep, including catalog browsing changes, variable/layer UX, animation preview updates, extraction fixes, and richer progress reporting."
-->

<!--
metadata:
- gpt_version: GPT-5 Codex
- time: 12:12
- date: 2026-07-05
- prompt_used: "Save progress notes while implementing the staged eReefs GUI so the work can continue in a new chat if needed."
-->
<!--
metadata:
- gpt_version: GPT-5 Codex
- time: 11:09
- date: 2026-07-06
- prompt_used: "Update the GUI progress notes to record that the Windows bslib smoke test now succeeds under unsandboxed R 4.4.0 but still fails in the restricted Codex sandbox."
-->
