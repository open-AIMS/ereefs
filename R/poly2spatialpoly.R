#' Convert a polygon dataframe in the format used by ggplot2::geom_polygon to a SpatVector
#'
#' Assumes that the dataframe only contains rombuses.
#'
#' @param polydf A polygon in the format required by geom_polygon: a data.frame with values, x and y positions and id labels
#' @return A SpatVector as used in the package 'terra'
#' @export

poly2sv <- function(polydf) {

   polydf <- polydf %>% dplyr::filter(!is.na(value))
   polytb <- dplyr::tibble(object = polydf$id,
                    part = 1,
                    x = polydf$x,
                    y = polydf$y,
                    hole = 0,
                    value = polydf$value) %>% dplyr::arrange(object)
   sv <- terra::vect(data.matrix(polytb %>% dplyr::select(-value)), "polygons")
   terra::values(sv) <- polytb %>% dplyr::group_by(object) %>%
                            dplyr::summarize(value = value[1]) %>% 
                            dplyr::ungroup()
   return(sv)
}

#' Convert a SpatVector to a raster. (Wrapper for terra::rasterize() with some reasonable defaults for the ereefs package)
#'
#' @param sv SpatVector object, e.g. as output by poly2sv()
#' @param xmn Minimum x value to use for the grid. Defaults to 142.45
#' @param ymn Minimum y value for the grid. Defaults to -27.5
#' @param resolution Grid resolution in degrees. Defaults to 0.01
#' @param xmx Maximum x value for the grid. Defaults to the maximum x extent of
#'   `sv`.
#' @param ymx Maximum y value for the grid. Defaults to the maximum y extent of
#'   `sv`.
#' @param r A raster object with the correct grid set up (optional). Default=NULL
#' @return A raster object, as used in the package 'terra'
#' @export

sv2raster <- function(sv, xmn=142.45, ymn=-27.5, resolution=0.01, xmx=NA, ymx=NA, r=NULL) { 
   # Default settings line up with Dieter's grid but encompass full extent of eReefs domain
   bbox <- as.vector(terra::ext(sv))
   if (is.na(xmx)) xmx<-bbox[2]
   if (is.na(ymx)) ymx<-bbox[4]
   ncols <- as.integer((xmx-xmn)/resolution)
   nrows <- as.integer((ymx-ymn)/resolution)
   xmx <- xmn + resolution * ncols
   ymx <- ymn + resolution * nrows
   if (is.null(r)) r <- terra::rast(sv, ncols=ncols, nrows=nrows, xmin=xmn, xmax=xmx, ymin=ymn, ymax=ymx)
   r <- terra::rasterize(sv, r, field='value')
}

#' Backward-compatible wrapper for older function name.
#'
#' @param polydf Polygon data frame in the format used by `geom_polygon()`.
#' @return A `terra::SpatVector`.
#' @export
poly2sp <- function(polydf) {
  poly2sv(polydf)
}

#' Backward-compatible wrapper for older function name.
#'
#' @inheritParams sv2raster
#' @return A `terra` raster.
#' @export
sp2raster <- function(sv, xmn=142.45, ymn=-27.5, resolution=0.01, xmx=NA, ymx=NA, r=NULL) {
  sv2raster(sv = sv, xmn = xmn, ymn = ymn, resolution = resolution, xmx = xmx, ymx = ymx, r = r)
}

# metadata:
# - gpt_version: GPT-5 Codex
# - time: 17:05
# - date: 2026-04-26
# - prompt_used: "Install the required packages, test the refactored eReefs package against THREDDS-served data, and set up a working Jupyter demo notebook."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 14:24
# - date: 2026-06-29
# - prompt_used: "Resolve the remaining raster-helper documentation mismatches reported by R CMD check for issues #24 and #25."
