#' ereefs: Some useful tools for exctracting and visualising eReefs and EMS netcdf data
#'
#' The ereefs package provides several functions:
#'    map_ereefs() Produces surface maps of model variables downloaded directly from NCI servers
#'                 or other EMS netcdf output files.
#'    map_ereefs_movie() Produces images for an animation of surface maps similar to those produced by map_ereefs.
#'    plume_class() Calculates plume optical classes. Mostly intended as a helper function for the above
#'                 two functions, but can be accessed directly.
#'    get_ereefs_ts() Extracts a time-series at any arbitrary location from an eReefs or other EMS netcdf
#'                 output file
#'    get_ereefs_depth_integrated_ts() Extracts a vertically integrated time-series at any arbitrary location
#'                 from an eReefs or other EMS netcdf output file.
#'    get_ereefs_depth_specified_ts() Extracts a time-series from a specified depth below the surface at
#'                 any arbitrary location within the EMS model domain.
#'    get_ereefs_profile() Extracts vertical profiles of specified variables from a specified geographic location,
#'                 over a specified period of time.
#'    get_ereefs_slice() Extracts a vertical slice or depth-resolved transect from an eReefs or other EMS netcdf 
#'                 output file.
#'    plot_ereefs_profile() Plots a single vertical profile obtained using get_ereefs_profile().
#'    plot_ereefs_slice() Plots a vertically-resolved transect obtained using get_ereefs_slice().
#'    plot_ereefs_zvt() Creates a contour plot showing the value of a variable over depth and time, using output
#'                 from get_ereefs_profile().
#'    poly2sp() Convert a polygon dataFrame in the format used by ggplot2::geom_plot to a spatialPolygonDataFrame,
#'                 as used in the package 'sp'.
#'    sp2raster() Wrapper to raster::rasterize to convert a spatialPolygonDataFrame (e.g. as produced by poly2sp())
#'                 to a raster that can be saved e.g. to a GeoTIF file.
#' @docType package
#' @name ereefs
NULL
