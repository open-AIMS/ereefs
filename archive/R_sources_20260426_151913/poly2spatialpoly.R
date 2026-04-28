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
#' @param xmax Maximum x value for the grid. Defaults to max value from the SpatVector extent
#' @param ymax Maximum y value for the grid. Defaults to max value from the SpatVector extent
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

