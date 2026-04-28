# ereefs development changes

This file documents user-visible changes in the major refactor currently on
`tidy_dev` relative to the older `master` workflow. It is focused on
upgrade-impacting changes rather than an exhaustive commit history.

## Breaking and near-breaking changes

### Data access now follows a `tidync`/`ncmeta`-first path

- The main extraction and plotting workflows now use `tidync` and `ncmeta`
  rather than the older `ncdf4`-centric OPeNDAP path.
- This was done to keep the package working against updated THREDDS/OPeNDAP
  services, especially on NCI.
- Local NetCDF files still work, but some low-level behaviour differs because
  array handling is now more explicit about dimension roles.

### Time-series helpers now return tibbles

- `get_ereefs_ts()`
- `get_ereefs_bottom_ts()`
- `get_ereefs_depth_integrated_ts()`
- `get_ereefs_depth_specified_ts()`

These functions now return tibbles instead of base data frames in normal use.
Code that depends on base `data.frame` behaviour, row names, or exact class
matching may need to be updated.

Recommended migration:

```r
ts <- get_ereefs_ts(...)
class(ts)
# tibble/data.frame
```

If older code expects a base data frame, convert explicitly:

```r
ts_df <- as.data.frame(ts)
```

### `layer = "surface"` and `layer = "bottom"` now mean wet surface and wet bottom

For depth-resolved extraction:

- `layer = "surface"` now means the shallowest wet model layer at that time
- `layer = "bottom"` now means the deepest wet model layer at that time

This is a behavioural change from workflows that implicitly assumed the largest
or smallest `k` index always corresponded to the water surface or seabed,
regardless of tide and wetting/drying.

This is more physically robust for shallow and tidally varying cells, but it can
change extracted values relative to older scripts.

### Missing `eta` now warns and falls back to `eta = 0`

For simple-format datasets that do not include free-surface elevation:

- depth-below-free-surface helpers now warn and assume `eta = 0`

This keeps simple files usable, but it changes failure behaviour from "stop"
toward "warn and continue".

### Missing `z_grid` can now be reconstructed from `zc`

When `z_grid` is absent but `zc` is available, the package reconstructs layer
interfaces from `zc` using midpoint assumptions and resets the top interface to
`1e20` to match standard EMS practice.

This improves support for simple files and some RECOM outputs, but users should
be aware that the reconstructed vertical grid is an inferred geometry rather
than a directly stored one.

### THREDDS catalogs are now preferred over stem-based filename reconstruction

Older workflows often reconstructed monthly filenames from file stems.

Older workflows also relied heavily on interactive or menu-style dataset
selection.

The refactored code now prefers:

- a specific NetCDF/OPeNDAP URL, or
- a THREDDS catalog URL

When a catalog is provided, the package now works out which files are needed to
cover the requested period instead of reconstructing filenames from stems.

Migration guidance:

- Prefer catalog URLs for live data requests
- Prefer explicit file or catalog URIs over menu-style selection
- Treat stem-based file reconstruction as legacy behaviour
- Treat `"menu"`, `"nci"`, `"catalog"`, `"old_menu"`, and numeric menu
  selections as backward-compatibility paths rather than the primary workflow

### Single-time requests can now snap to the nearest available model time

If a requested timestamp does not exactly match an output time, the package now:

- finds the closest available time
- emits a warning
- proceeds with that time

Older code may have assumed an exact match or may have failed earlier.

### Date and time handling is broader and more explicit

The refactored code consistently accepts a wider mix of date inputs across the
main plotting and extraction functions, including:

- `Date`
- `POSIXct`
- character date or date-time strings
- numeric date stamps such as `YYYYMMDD`
- vectors such as `c(year, month, day)`

Important behaviour details:

- date-only inputs are interpreted in `Etc/GMT-10`
- date-only inputs default to midday when converted to `POSIXct`
- single-time requests may use the nearest available model time with a warning

Older scripts that assumed midnight defaults, exact timestamp matching, or only
vector-style dates should be checked.

### Grid handling now supports both curvilinear and centre-only regular grids

The package now supports:

- curvilinear EMS grids with cell corners
- regridded regular products with centre coordinates only

For centre-only grids, plotting geometry is reconstructed from the centre
coordinates instead of requiring stored corner files by default.

This changes some `get_ereefs_grids()` outputs and plotting internals:

- `get_ereefs_grids()` now returns a richer list including `spatial_grid` and
  `grid_type`
- regular-grid plotting may use synthesised corners rather than archived corner
  files

### Plot defaults have changed

- Scalar plotting defaults now use `viridis`
- `plot_ereefs_slice()` now defaults to `viridis`
- smooth maps are now an explicit display option in `map_ereefs()`

Scripts that relied on the old default palette should now set `scale_col`
explicitly if they need legacy colours.

### Plot objects are now designed for `ggplot2` workflows

The active map and slice plotting paths now consistently return `ggplot2`
objects. This makes them easier to extend, but code that relied on older
side-effect-only plotting patterns may need small changes.

### `map_ereefs_movie()` now saves frames and can assemble animations

`map_ereefs_movie()` now supports saving frame PNGs and assembling them into a
GIF or MP4-style movie workflow, rather than only generating plots ad hoc.

Notable behaviour changes:

- frame images can be saved to a folder
- one shared colour scale is kept across frames
- animation output is now part of the supported workflow

### Vertical slices and transects now accept multi-point paths

`get_ereefs_slice()` can now be used for:

- straight transects
- curved paths
- zig-zag paths

by passing a multi-row latitude/longitude object to `geolocation`.

This is more flexible than the older emphasis on simple two-point slices.

### Returned profile and slice objects are more structured around tidy workflows

The profile and slice helpers still return lists, but several components are now
used in a more tidy/tibble-oriented way internally, and some downstream code may
see differences in:

- column ordering
- object classes
- the way cell matches and cross-references are represented

If you have custom downstream scripts, check them against the current output of:

- `str(get_ereefs_profile(...))`
- `str(get_ereefs_slice(...))`
- `str(get_ereefs_grids(...))`

## Recommended migration checks

If you have older scripts, check these first:

1. Convert strict `data.frame` assumptions to accept tibbles.
2. Check any code that assumes a fixed `k` index is always the surface or
   bottom.
3. Prefer THREDDS catalog URLs over stem-based filename logic.
4. Set `scale_col` explicitly if you need a legacy palette.
5. Re-check any code that expected exact timestamp matching.
6. Re-check custom plotting or GIS code that depends on corner grids or exact
   object structure.

## Documentation pointers

- Main README: `README.md`
- Usage vignette: `vignettes/howto.Rmd`
- Comparison/background vignette: `vignettes/about.Rmd`
- Live and local demo notebooks: `notebooks/`
