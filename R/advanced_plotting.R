#' Create a surface map of eReefs model output using a different colour scale to highlight areas where a 
#' second variable is within a specified range. Scale is not squished. Google map underlay suppressed by default.
#'
#' Creates a colour map showing concentrations of a specified eReefs model output variable at a specified
#' model layer (by default, the surface layer). The map is optionally (and by default) overlain on a Google 
#' Earth map of the region. The function should work under Linux or MacOS. The R netcdf4 library doesn't 
#' handle opendap served files under Windows, unfortunately, but can be used under Windows for locally-stored
#' model output files. By Barbara Robson (AIMS).
#'
#' References:
#'
#' Baird, M.E., Cherukuru, N., Jones, E., Margvelashvili, N., Mongin, M., Oubelkheir, K., 
#'   Ralph, P.J., Rizwi, F., Robson, B.J., Schroeder, T. and Skerratt, J., 2016. Remote-sensing 
#'   reflectance and true colour produced by a coupled hydrodynamic, optical, sediment, 
#'   biogeochemical model of the Great Barrier Reef, Australia: comparison with satellite data. 
#'   Environmental Modelling & Software, 78, pp.79-96.
#'
#' Baird, M.E., Andrewartha, J., Herzfeld, M., Jones, E.M., Margvelashvili, N., Mongin, M., Rizwi, 
#'   F.,  Skerratt, J., Soja-Wozniak, M., Wild-Allen, K., Schroe der, T., Robson, B.J., da Silva, E. 
#'   and Devlin, M., 2017. River plumes of the Great Barrier Reef: freshwater , sediment and optical 
#'   footprints quantified by the eReefs modelling system  In Syme, G., Hatton MacDonald, D., Fulton, 
#'   B. and Piantadosi, J. (eds) MODSIM2017, 22nd International Congress on Modelling and Simulation. 
#'   Modelling and Simulation Society of Australia and New Zealand, December 2017, pp. 894–900
#'
#' Devlin, M.J., McKinna, L.W., Álvarez-Romero, J.G., Petus, C., Abott, B., Harkness, P. and 
#'   Brodie, J., 2012. Mapping the pollutants in surface riverine flood plume waters in the Great 
#'   Barrier Reef, Australia. Marine pollution bulletin, 65(4-9), pp.224-235.
#'
#'
#' @return a ggplot object
#' @param var_name Short name of the variable to plot. This can be any variable in the
#'                 eReefs netcdf file that you are accessing (refer to eReefs model
#'                 documentation or extract variable names from the netcdf file for a full 
#'                 list). In addition, two special var_name values are supported: 'true_colour'
#'                 and 'plume'. 'true_colour' provides a simulated MODIS true colour
#'                 image (refer to Baird et al., 2016 for en explanation).'plume'
#'                 provides a map of calculated plume colour class as per Devlin et al. (2012)
#'                 and Baird et al. (2017). Defaults to true_colour.
#' @param filter_name Short name of the variable by which to filter output. Defaults to 'salt'.
#' @param filter_range. Show var_name only where filter_name is between the upper and lower value.
#'        Defaults to c(0, 34.5).
#' @param target_date Date to map. Can either be a) a vector of integers in the format 
#'	c(year, month, day); b) a date obtained e.g. from as.Date(); or c) a character string 
#'      formatted for input to as.Date(). Defaults to c(2018,1,30). Be careful to choose the right
#'      input file, as the function will plot the closest date in the file to the target date without
#'      complaint, however far fron the target that may be. Assumes that dates in the netcdf files are
#'      relative to 1990-01-01 (this is not checked).
#' @param layer Either an integer layer number or 'surface' to choose the surface layer. Defaults to 'surface'.
#' @param Google_map_underlay Set to TRUE (the default) to use ggmap to show a Google Map as
#'      an underlay for the model output plot. Requires the ggmap librray.
#' @param input_stem Path to the locally-stored or opendap-served netcdf file, including the stem
#'      of the file name (i.e. the start of the filename without the date or file extension) if the filename
#'       contains 'gbr1' or 'gbr4', or the full filename otherwise.  Defaults to 
#'      "http://dapds00.nci.org.au/thredds/dodsC/fx3/gbr4_bgc_GBR4_H2p0_B2p0_Chyd_Dcrt/gbr4_bgc_simple_".,
#' @param filter_stem Path to the locally-stored or opendap-served netcdf file, including the stem
#'      of the file name (i.e. the start of the filename without the date or file extension) for the filter variable. 
#'      if the filename contains 'gbr1' or 'gbr4', or the full filename otherwise. Must be on the same grid as 
#'      input_stem. Defaults to the value of inout_stem.
#' @param input_grid Name of the locally-stored or opendap-served netcdf file that contains the grid
#'      coordinates for the cell corners (x_grid and y_grid). If not specified, the function will first look for
#'      x_grid and y_grid can be found in the first INPUT_STEM file, and if not found, will check whether the size 
#'      of the variables in the input file corresponds to the size expected for GBR4 or GBR1, and load appropriate 
#'      x and y grids from data files stored in this package. Alternatively, you can provide the location of a full 
#'      (not simple-format) ereefs netcdf output file such as 
#'      "http://dapds00.nci.org.au/thredds/dodsC/fx3/gbr4_hydro_all/gbr4_all_2016-09.nc"
#' @param scale_col Vector of colours to use for low and high values in the colour scale. This can be a colour 
#'      from the ggplot colour palette or a RGB hash code. Ignored for true_colour plots. 
#'      Defaults to c('ivory', 'hotpink').
#' @param filter_col Colours to use for high values in the colour scale in areas matching
#'      the filter criteria.. This can be a colour from the ggplot colour palette or a RGB hash code.
#'      Defaults to 'blue'.
#' @param scale_lim Upper and lower bounds for colour scale. Defaults to full range of data.
#'      Ignored for true_colour plots.
#' @param zoom Value of zoom passed to ggmap(). Set to 5 if you want to show the entire extent 
#'      of eReefs models. Defaults to 6. Higher values will zoom in further.
#' @param box_bounds Minimum and maximum latitude and longitude coordinates to map. Defaults to the
#'        entire extent of the model output (though modified by the value of zoom). 
#'        Format: c(longitude_min, longitude_max, latitude_min, latitude_max).
#' @export
#' @examples
#' map_ereefs_filter()
map_ereefs_filter <- function(var_name = "true_colour", 
		       filter_name = "salt",
		       filter_range = c(0, 34.5),
		       target_date = c(2018,1,30), 
                       layer = 'surface',
		       Google_map_underlay = FALSE,
		       input_stem = "http://dapds00.nci.org.au/thredds/dodsC/fx3/gbr4_bgc_GBR4_H2p0_B2p0_Chyd_Dnrt/gbr4_bgc_simple_",
		       filter_stem = NA,
                       input_grid = NA,
                       scale_col = c('ivory', 'hotpink'),
                       filter_col = 'blue',
		       scale_lim = c(NA, NA),
                       zoom = 6,
		       box_bounds = c(NA, NA, NA, NA))
{

# Check whether this is a GBR1 or GBR4 ereefs file, or something else
ereefs_case <- 0

if (stringi::stri_detect_fixed(input_stem, 'gbr4')) {
	ereefs_case <- 4
} else if (stringi::stri_detect_fixed(input_stem ,'gbr1')) {
	ereefs_case <- 1
}

# Date to map:
if (is.vector(target_date)) {
	target_date <- as.Date(paste(target_date[1], target_date[2], target_date[3], sep='-'))
} else if (is.character(target_date)) {
	target_date <- as.Date(target_date)
}
if (ereefs_case == 4) { 
	day <- as.integer(format(target_date, '%d'))
	filename <- paste0(input_stem, format(target_date, '%Y-%m'), '.nc')
	filterfile <- paste0(filter_stem, format(target_date, '%Y-%m'), '.nc')
} else if (ereefs_case == 1) {
	day <- 1
	filename <- paste0(input_stem, format(target_date, '%Y-%m-%d'), '.nc')
	filterfile <- paste0(filter_stem, format(target_date, '%Y-%m-%d'), '.nc')
} else {
	filename <- paste0(input_stem, '.nc')
	filterfile <- paste0(filter_stem, '.nc')
	nc <- ncdf4::nc_open(filename)
        ds <- try(as.Date(ncdf4::ncvar_get(nc, "time"), origin = as.Date("1990-01-01")), TRUE)
	if (class(ds)=="try-error") {
	  # We may be dealing with an old EMS file
          print('time not found in netcdf file. Looking for t instead.')
	  ds <- as.Date(ncdf4::ncvar_get(nc, "t"), origin = as.Date("1990-01-01"))
          print('t found in netcdf file. Ignore previous error message.')
	}
	day <- which.min(abs(target_date - ds))
	ncdf4::nc_close(nc)
}

# Allow for US English:
if (var_name == "true_color") {
	var_name <- "true_colour"
}

# Check for ggmap()
if ((Google_map_underlay)&&(!requireNamespace("ggmap", quietly = TRUE))) {
  warning('Package ggmap is required to show a Google map underlay. Preparing plot with no underlay.')
  print('To avoid this message, either install ggmap or set Google_map_underlay = FALSE.')
  Google_map_underlay = FALSE
}

# Get cell grid corners
if (!is.na(input_grid)) {
  nc2 <- ncdf4::nc_open(input_grid)
  x_grid <- ncdf4::ncvar_get(nc2, "x_grid")
  y_grid <- ncdf4::ncvar_get(nc2, "y_grid")
  ncdf4::nc_close(nc2)
} else if (ereefs_case == 4) { 
     # This looks like GBR4 
     x_grid <- gbr4_x_grid
     y_grid <- gbr4_y_grid
} else if (ereefs_case == 1) { 
     # This looks like GBR1 
     x_grid <- gbr1_x_grid
     y_grid <- gbr1_y_grid
} else { 
     # Look for x_grid and y_grid in the main input file
     nc2 <- ncdf4::nc_open(filename)
     x_grid <- try(ncdf4::ncvar_get(nc2, "x_grid"), TRUE)
     if (class(x_grid) == "try-error") {
       ncdf4::nc_close(nc2)
       stop(paste("Unfamiliar netcdf file variable dimensions:", filename, 
		  "doesn't look like a GBR1 or GBR4 grid and x_grid and y_grid were not found in this file.",
		  "Please specify a filename for input_grid."))
     } 
     y_grid <- ncdf4::ncvar_get(nc2, "y_grid")
     ncdf4::nc_close(nc2)
}

dims <- dim(x_grid) - 1

outOfBox <- array(FALSE, dim=dim(x_grid))
if (!is.na(box_bounds[1])) {
  outOfBox <- apply(x_grid,2,function(x){ (x<box_bounds[1]|is.na(x)) } )
}
if (!is.na(box_bounds[2])) {
  outOfBox <- outOfBox | apply(x_grid,2,function(x){(x>box_bounds[2]|is.na(x))})
}
if (!is.na(box_bounds[3])) {
  outOfBox <- outOfBox | apply(y_grid,2,function(x){(x<box_bounds[3]|is.na(x))})
}
if (!is.na(box_bounds[4])) {
  outOfBox <- outOfBox | apply(y_grid,2,function(x){(x>box_bounds[4]|is.na(x))})
}

if (is.na(box_bounds[1])) { 
  xmin <- 1
} else {
  xmin <- which(apply(!outOfBox, 1, any))[1]
  if (length(xmin)==0) xmin <- 1
}
if (is.na(box_bounds[2])) {
  xmax <- dims[1]
} else {
  xmax <- which(apply(!outOfBox, 1, any))
  xmax <- xmax[length(xmax)]
  if ((length(xmax)==0)|(xmax > dims[1])) xmax <- dims[1]
}
if (is.na(box_bounds[3])) { 
  ymin <- 1
} else {
  ymin <- which(apply(!outOfBox, 2, any))[1]
  if (length(ymin)==0) ymin <- 1
}
if (is.na(box_bounds[4])) {
  ymax <- dims[2]
} else {
  ymax <- which(apply(!outOfBox, 2, any))
  ymax <- ymax[length(ymax)]
  if ((length(ymax)==0)|(ymax > dims[2])) ymax <- dims[2]
}

x_grid <- x_grid[xmin:(xmax+1), ymin:(ymax+1)]
y_grid <- y_grid[xmin:(xmax+1), ymin:(ymax+1)]


if (var_name=="plume") {
    #inputfile <- paste0(filename, '?R_412,R_443,R_488,R_531,R_547,R_667,R_678')
    inputfile <- filename
    nc <- ncdf4::nc_open(inputfile)
    R_412 <- ncdf4::ncvar_get(nc, "R_412", start=c(xmin,ymin,day), count=c(xmax-xmin,ymax-ymin,1))
    R_443 <- ncdf4::ncvar_get(nc, "R_443", start=c(xmin,ymin,day), count=c(xmax-xmin,ymax-ymin,1))
    R_488 <- ncdf4::ncvar_get(nc, "R_488", start=c(xmin,ymin,day), count=c(xmax-xmin,ymax-ymin,1))
    R_531 <- ncdf4::ncvar_get(nc, "R_531", start=c(xmin,ymin,day), count=c(xmax-xmin,ymax-ymin,1))
    R_547 <- ncdf4::ncvar_get(nc, "R_547", start=c(xmin,ymin,day), count=c(xmax-xmin,ymax-ymin,1))
    R_667 <- ncdf4::ncvar_get(nc, "R_667", start=c(xmin,ymin,day), count=c(xmax-xmin,ymax-ymin,1))
    R_678 <- ncdf4::ncvar_get(nc, "R_678", start=c(xmin,ymin,day), count=c(xmax-xmin,ymax-ymin,1))
    rsr <- list(R_412, R_443, R_488, R_531, R_547, R_667, R_678)
    ems_var <- plume_class(rsr)
    dims <- dim(ems_var)

} else if (var_name=="true_colour") {
    #inputfile <- paste0(filename, '?R_470,R_555,R_645,eta')
    inputfile <- filename
    nc <- ncdf4::nc_open(inputfile)
    TCbright <- 10
    R_470 <- ncdf4::ncvar_get(nc, "R_470", start=c(xmin,ymin,day), count=c(xmax-xmin,ymax-ymin,1)) * TCbright
    R_555 <- ncdf4::ncvar_get(nc, "R_555", start=c(xmin,ymin,day), count=c(xmax-xmin,ymax-ymin,1)) * TCbright
    R_645 <- ncdf4::ncvar_get(nc, "R_645", start=c(xmin,ymin,day), count=c(xmax-xmin,ymax-ymin,1)) * TCbright
    R_470[R_470>1] <- 1
    R_555[R_555>1] <- 1
    R_645[R_645>1] <- 1

    unscaledR = c(0, 30, 60, 120, 190, 255)/255;
    scaledR = c(1, 110, 160, 210, 240, 255)/255;
    scalefun <- approxfun(x=unscaledR, y=scaledR, yleft=1, yright=255)
    red <- scalefun(R_645)
    green <-scalefun(R_555)
    blue <- scalefun(R_470)
    red[is.na(red)] <- 0
    green[is.na(green)] <- 0
    blue[is.na(blue)] <- 0

    ems_var <- rgb(red, green, blue)
    ems_var[ems_var=="#000000"] <- NA
    ems_var <-array(as.character(ems_var), dim=dim(R_645))
    dims <- dim(ems_var)
} else { 
    #inputfile <- paste0(filename, '?', var_name)
    inputfile <- filename
    nc <- ncdf4::nc_open(inputfile)
      # We don't yet know the dimensions of the variable, so let's get them
      dims <- nc$var[[var_name]][['size']]
      if (is.null(dims)) stop(paste(var_name, ' not found in netcdf file.')) 
      ndims <- length(dims)
      if ((ndims > 3) && (layer == 'surface')) layer <- dims[3]
}

if (is.na(filter_stem)) {
   ncf <- nc
} else {
   ncf <- ncdf4::nc_open(filterfile)
}
dims <- ncf$var[[filter_name]][['size']]
if (is.null(dims)) stop(paste(filter_name, ' not found in netcdf file.')) 
fdims <- length(dims)
if (ndims == 4) {
   filter_var <- ncdf4::ncvar_get(ncf, filter_name, start=c(xmin,ymin,layer,day), count=c(xmax-xmin,ymax-ymin,1,1))
} else {
   filter_var_var <- ncdf4::ncvar_get(ncf, filter_name, start=c(xmin,ymin,day), count=c(xmax-xmin,ymax-ymin,1))
}


if (!((var_name == 'true_colour') || (var_name == 'plume'))) {
    if (ndims == 4) {
       ems_var <- ncdf4::ncvar_get(nc, var_name, start=c(xmin,ymin,layer,day), count=c(xmax-xmin,ymax-ymin,1,1))
    } else {
       ems_var <- ncdf4::ncvar_get(nc, var_name, start=c(xmin,ymin,day), count=c(xmax-xmin,ymax-ymin,1))
    }
}


ncdf4::nc_close(nc)
if (!is.na(filter_stem)) ncdf4::nc_close(ncf)

a <- dim(ems_var)[1]
b <- dim(ems_var)[2]

# Set up the polygon corners. 4 per polygon.
gx <- c(x_grid[1:a, 1:b], x_grid[2:(a+1), 1:b], x_grid[2:(a+1), 2:(b+1)], x_grid[1:a, 2:(b+1)])
gy <- c(y_grid[1:a, 1:b], y_grid[2:(a+1), 1:b], y_grid[2:(a+1), 2:(b+1)], y_grid[1:a, 2:(b+1)])
gx <- array(gx, dim=c(a*b,4))
gy <- array(gy, dim=c(a*b,4))

# Find and exclude points where not all corners are defined
gx_ok <- !apply(is.na(gx),1, any)
gy_ok <- !apply(is.na(gy),1, any)

range_ok <- c((filter_var >= filter_range[1]) & (filter_var <= filter_range[2]))

# Values associated with each polygon at chosen timestep
#n <- c(ems_var[,,tstep])[gx_ok&gy_ok]
if (var_name=="true_colour") { 
    n <- c(ems_var)[gx_ok&gy_ok]
    # (gx_ok and gy_ok should be identical, but let's be certain)
    gx <- c(t(gx[gx_ok&gy_ok,]))
    gy <- c(t(gy[gx_ok&gy_ok,]))
} else {
    ems_ok <- c(ems_var!="#000000")
    n <- c(ems_var)[gx_ok&gy_ok&ems_ok]
    n2 <- c(ems_var)[gx_ok&gy_ok&ems_ok&range_ok]
    # (gx_ok and gy_ok should be identical, but let's be certain)
    gx2 <- c(t(gx[gx_ok&gy_ok&ems_ok&range_ok,]))
    gy2 <- c(t(gy[gx_ok&gy_ok&ems_ok&range_ok,]))
    gx <- c(t(gx[gx_ok&gy_ok&ems_ok,]))
    gy <- c(t(gy[gx_ok&gy_ok&ems_ok,]))
}

n[n<scale_lim[1]] <- scale_lim[1]
n[n>scale_lim[2]] <- scale_lim[2]
n <- n-scale_lim[1]
n2[n2<scale_lim[1]] <- scale_lim[1]
n2[n2>scale_lim[2]] <- scale_lim[2]
n2 <- n2-scale_lim[1]
n2 <- -n2

# Unique ID for each polygon
id <- 1:length(n)

id <- as.factor(id)
values <- data.frame(id = id, value = n)
positions <- data.frame(id=rep(id, each=4), x = gx, y = gy)
datapoly <- merge(values, positions, by = c("id"))

id2 <- 1:length(n2)
id2 <- as.factor(id2)
values2 <- data.frame(id = id2, value = n2)
positions2 <- data.frame(id=rep(id2, each=4), x = gx2, y = gy2)
datapoly2 <- merge(values2, positions2, by = c("id"))

if ((var_name!="true_colour")&&(is.na(scale_lim[1]))) { 
	scale_lim <- c(min(n, na.rm=TRUE), max(n, na.rm=TRUE))
}

if (Google_map_underlay) {
  MapLocation<-c(min(gx, na.rm=TRUE)-0.5, 
 		min(gy, na.rm=TRUE)-0.5, 
 		max(gx, na.rm=TRUE)+0.5, 
 		max(gy, na.rm=TRUE)+0.5)
  myMap<-suppressWarnings(ggmap::get_map(location=MapLocation, source="google", maptype="hybrid", zoom=zoom, crop=TRUE))
  p <- ggmap::ggmap(myMap)
} else {
  p <- ggplot()
}
if (var_name=="true_colour") {
    p <- p +
        ggplot2::geom_polygon(ggplot2::aes(x=x, y=y, fill=value, group=id), data = datapoly) +
	ggplot2::scale_fill_identity()
} else {
    p <- p +
        ggplot2::geom_polygon(ggplot2::aes(x=x, y=y, fill=value, group=id), data = datapoly) +
        ggplot2::geom_polygon(ggplot2::aes(x=x, y=y, fill=value, group=id), data = datapoly2) +
        ggplot2::scale_fill_gradient2(low=filter_col,
				     mid=scale_col[1], 
				     high=scale_col[2], 
				     na.value="transparent", 
				     guide="none") +
        ggplot2::xlab('latitude') +
	ggplot2::ylab('longitude')
}
#p <- datapoly %>% ggvis(x = ~x, y = ~y, fill = ~value) %>% layer_paths() %>% group_by(id)
print(p)
return(p)
}

