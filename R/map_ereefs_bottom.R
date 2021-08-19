#' Create a benthic layer map of eReefs model output.
#'
#' Creates a colour map showing concentrations of a specified eReefs model output variable in the
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
#' @param target_date Date to map. Can either be a) a vector of integers in the format 
#'	c(year, month, day); b) a date obtained e.g. from as.Date(); or c) a character string 
#'      formatted for input to as.Date(). Defaults to c(2018,1,30). Be careful to choose the right
#'      input file, as the function will plot the closest date in the file to the target date without
#'      complaint, however far fron the target that may be. Assumes that dates in the netcdf files are
#'      relative to 1990-01-01 (this is not checked).
#' @param layer Either an integer layer number or 'surface' to choose the surface layer. Defaults to 'surface'.
#' @param Google_map_underlay Set to TRUE (the default) to use ggmap to show a Google Map as
#'      an underlay for the model output plot. Requires the ggmap librray.
#' @param input_file is the URI or file location of any of the EMS output files, 
#'        Defaults to "http://dapds00.nci.org.au/thredds/dodsC/fx3/gbr4_bgc_GBR4_H2p0_B2p0_Chyd_Dnrt/gbr4_bgc_simple_2018-03.nc"
#'        If using Windows, you will need to set this to a local inputfile stem.
#' @param input_grid Name of the locally-stored or opendap-served netcdf file that contains the grid
#'      coordinates for the cell corners (x_grid and y_grid). If not specified, the function will first look for
#'      x_grid and y_grid can be found in the first INPUT_STEM file, and if not found, will check whether the size 
#'      of the variables in the input file corresponds to the size expected for GBR4 or GBR1, and load appropriate 
#'      x and y grids from data files stored in this package. Alternatively, you can provide the location of a full 
#'      (not simple-format) ereefs netcdf output file such as 
#'      "http://dapds00.nci.org.au/thredds/dodsC/fx3/gbr4_hydro_all/gbr4_all_2016-09.nc"
#' @param scale_col Vector of colours to use for low and high values in the colour scale. This can be a colour 
#'      from the ggplot colour palette or a RGB hash code. Ignored for true_colour plots. 
#'      Defaults to c('ivory', 'coral4').
#' @param scale_lim Upper and lower bounds for colour scale. Defaults to full range of data.
#'      Ignored for true_colour plots.
#' @param zoom Value of zoom passed to ggmap(). Set to 5 if you want to show the entire extent 
#'      of eReefs models. Defaults to 6. Higher values will zoom in further.
#' @param box_bounds Minimum and maximum latitude and longitude coordinates to map. Defaults to the
#'        entire extent of the model output (though modified by the value of zoom). 
#'        Format: c(longitude_min, longitude_max, latitude_min, latitude_max).
#' @param p Handle for an existing figure if you want to add a layer instead of creating a new figure.
#'        If p is provided, Google_map_underlay is over-ridden and set to FALSE.
#' @param return_poly Instead of only the figure handle, return a list containing the figure handle and the dataframe used by geom_plot(). Default FALSE.
#' @export
#' @examples
#' map_ereefs()
#' map_ereefs('TN')
#' map_ereefs('plume', target_date=c(2011, 1, 30))
#' map_ereefs('Chl_a_sum', target_date='2016-07-15', scale_col=c('ivory', 'green4'))
#' map_ereefs('salt', box_bounds=c(145,150,-20,-15), zoom=7, scale_lim=c(32,35))
map_ereefs <- function(var_name = "true_colour", 
                       target_date = c(2018,1,30), 
                       layer = 'surface', 
                       Google_map_underlay = TRUE,
                       input_file = "http://dapds00.nci.org.au/thredds/dodsC/fx3/gbr4_bgc_GBR4_H2p0_B2p0_Chyd_Dnrt/gbr4_bgc_simple_2018-03.nc",
                       input_grid = NA,
                       scale_col = c('ivory', 'coral4'), 
                       scale_lim = c(NA, NA),
                       zoom = 6, 
                       box_bounds = c(NA, NA, NA, NA), 
                       p = NA, 
                       suppress_print = FALSE, 
                       return_poly = FALSE)
{

if (length(p)!=1) Google_map_underlay <- FALSE
if (suppress_print) Google_map_underlay <- FALSE

# Check whether this is a GBR1 or GBR4 ereefs file, or something else
ereefs_case <- get_ereefs_case(input_file)
input_stem <- get_file_stem(input_file)
check_platform_ok(input_stem)
grids <- get_ereefs_grids(input_file, input_grid)
x_grid <- grids[['x_grid']]
y_grid <- grids[['y_grid']]

# Date to map:
if (is.vector(target_date)) {
	target_date <- as.Date(paste(target_date[1], target_date[2], target_date[3], sep='-'))
} else if (is.character(target_date)) {
	target_date <- as.Date(target_date)
}
if (ereefs_case[2] == '4km') { 
	filename <- paste0(input_stem, format(target_date, '%Y-%m'), '.nc')
	nc <- ncdf4::nc_open(filename)
	if (!is.null(nc$var[['t']])) { 
	    ds <- as.Date(ncdf4::ncvar_get(nc, "t"), origin = as.Date("1990-01-01"))
        } else {
	    ds <- as.Date(ncdf4::ncvar_get(nc, "time"), origin = as.Date("1990-01-01"))
	}
	day <- which.min(abs(target_date - ds))
	ncdf4::nc_close(nc)
} else if (ereefs_case[2] == '1km') {
	day <- 1
	ds <- target_date
	filename <- paste0(input_stem, format(target_date, '%Y-%m-%d'), '.nc')
} else {
	filename <- paste0(input_stem, '.nc')
	nc <- ncdf4::nc_open(filename)
	if (!is.null(nc$var[['t']])) { 
	    ds <- as.Date(ncdf4::ncvar_get(nc, "t"), origin = as.Date("1990-01-01"))
        } else {
	    ds <- as.Date(ncdf4::ncvar_get(nc, "time"), origin = as.Date("1990-01-01"))
	}
	day <- which.min(abs(target_date - ds))
	ncdf4::nc_close(nc)
}

# Allow for US English:
if (var_name == "true_color") {
	var_name <- "true_colour"
}

# Check for ggmap()
if ((Google_map_underlay)&(!requireNamespace("ggmap", quietly = TRUE))) {
  warning('Package ggmap is required to show a Google map underlay. Preparing plot with no underlay.')
  print('To avoid this message, either install ggmap or set Google_map_underlay = FALSE.')
  Google_map_underlay = FALSE
}

# Get cell grid corners
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
    inputfile <- paste0(filename, '?R_412,R_443,R_488,R_531,R_547,R_667,R_678,longitude,latitude')
    #inputfile <- filename
    nc <- ncdf4::nc_open(inputfile)
    longitude <- ncdf4::ncvar_get(nc, 'longitude')
    latitude <-ncdf4:: ncvar_get(nc, 'latitude')
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
    var_units <- ''
    var_longname <- 'Plume colour class'

} else if (var_name=="true_colour") {
    inputfile <- paste0(filename, '?R_470,R_555,R_645,eta,longitude,latitude')
    #inputfile <- filename
    nc <- ncdf4::nc_open(inputfile)
    longitude <-ncdf4:: ncvar_get(nc, 'longitude')
    latitude <-ncdf4:: ncvar_get(nc, 'latitude')
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
    var_longname <- "Simulated true colour"
    var_units <- ""
} else if (var_name == 'ZooT') {
    inputfile <- paste0(filename, '?ZooL_N,ZooS_N,longitude,latitude')
    #inputfile <- filename
    nc <- ncdf4::nc_open(inputfile)
    longitude <-ncdf4:: ncvar_get(nc, 'longitude')
    latitude <-ncdf4:: ncvar_get(nc, 'latitude')
    # We don't yet know the dimensions of the variable, so let's get them
    dims <- nc$var[['ZooL_N']][['size']]
    if (is.null(dims)) stop(paste('ZooL_N', ' not found in netcdf file.')) 
    ndims <- length(dims)
    if ((ndims > 3) && (layer == 'surface')) layer <- dims[3]
    var_longname <- "Total zooplankton nitrogen"
    var_units <- "mg N m-3"
} else if (var_name == 'speed') {
    inputfile <- paste0(filename, '?u,v,longitude,latitude')
    #inputfile <- filename
    nc <- ncdf4::nc_open(inputfile)
    longitude <-ncdf4:: ncvar_get(nc, 'longitude')
    latitude <-ncdf4:: ncvar_get(nc, 'latitude') 
    # We don't yet know the dimensions of the variable, so let's get them 
    dims <- nc$var[['u']][['size']] 
    if (is.null(dims)) stop(paste('u', ' not found in netcdf file.')) 
    ndims <- length(dims) 
    if ((ndims > 3) && (layer == 'surface')) layer <- dims[3]
    var_longname <- "Current speed"
    var_units <- "m s-1"
} else { 
    inputfile <- paste0(filename, '?', var_name, ',longitude,latitude')
    #inputfile <- filename
    nc <- ncdf4::nc_open(inputfile)
    longitude <-ncdf4:: ncvar_get(nc, 'longitude')
    latitude <-ncdf4:: ncvar_get(nc, 'latitude') 
    # We don't yet know the dimensions of the variable, so let's get them 
    dims <- nc$var[[var_name]][['size']] 
    if (is.null(dims)) stop(paste(var_name, ' not found in netcdf file.')) 
    ndims <- length(dims) 
    if ((ndims > 3) && (layer == 'surface')) layer <- dims[3]
}

if (var_name == 'ZooT') {
    var_longname <- 'Total Zooplankton Nitrogen'
    var_units <- 'mg N m3'
    if (ndims == 4) {
       ems_var <- ncdf4::ncvar_get(nc, 'ZooL_N', start=c(xmin,ymin,layer,day), count=c(xmax-xmin,ymax-ymin,1,1))
       ems_var <- ems_var + ncdf4::ncvar_get(nc, 'ZooS_N', start=c(xmin,ymin,layer,day), count=c(xmax-xmin,ymax-ymin,1,1))
    } else {
       ems_var <- ncdf4::ncvar_get(nc, 'ZooL_N', start=c(xmin,ymin,day), count=c(xmax-xmin,ymax-ymin,1))
       ems_var <- ems_var + ncdf4::ncvar_get(nc, 'ZooS_N', start=c(xmin,ymin,day), count=c(xmax-xmin,ymax-ymin,1))
    }
} else if (var_name == 'speed') {
    var_longname <- 'Current speed'
    var_units <- 'm s-1'
    if (ndims == 4) {
       ems_var <- ncdf4::ncvar_get(nc, 'u', start=c(xmin,ymin,layer,day), count=c(xmax-xmin,ymax-ymin,1,1))
       ems_var <- sqrt(ems_var^2 + ncdf4::ncvar_get(nc, 'v', start=c(xmin,ymin,layer,day), count=c(xmax-xmin,ymax-ymin,1,1))^2)
    } else {
       ems_var <- ncdf4::ncvar_get(nc, 'ZooL_N', start=c(xmin,ymin,day), count=c(xmax-xmin,ymax-ymin,1))
       ems_var <- ems_var + ncdf4::ncvar_get(nc, 'ZooS_N', start=c(xmin,ymin,day), count=c(xmax-xmin,ymax-ymin,1))
    }
} else if (!((var_name == 'true_colour') || (var_name == 'plume'))) {
    vat <- ncdf4::ncatt_get(nc, var_name)
    var_longname <- vat$long_name
    var_units <- vat$units
    if (ndims == 4) {
       ems_var <- ncdf4::ncvar_get(nc, var_name, start=c(xmin,ymin,layer,day), count=c(xmax-xmin,ymax-ymin,1,1))
    } else {
       ems_var <- ncdf4::ncvar_get(nc, var_name, start=c(xmin,ymin,day), count=c(xmax-xmin,ymax-ymin,1))
    }
}

ncdf4::nc_close(nc)

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

# Values associated with each polygon at chosen timestep
#n <- c(ems_var[,,tstep])[gx_ok&gy_ok]
n <- c(ems_var)[gx_ok&gy_ok]
# (gx_ok and gy_ok should be identical, but let's be certain)
gx <- c(t(gx[gx_ok&gy_ok,]))
gy <- c(t(gy[gx_ok&gy_ok,]))
longitude <- c(longitude)[gx_ok&gy_ok]
latitude <- c(latitude)[gx_ok&gy_ok]

# Unique ID for each polygon
id <- 1:length(n)

id <- as.factor(id)
values <- data.frame(id = id, value = n)
positions <- data.frame(id=rep(id, each=4), x = gx, y = gy)
datapoly <- merge(values, positions, by = c("id"))

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
} else if (length(p)==1) {
  p <- ggplot2::ggplot()
}
if (var_name=="true_colour") {
  p <- p +
        ggplot2::geom_polygon(ggplot2::aes(x=x, y=y, fill=value, group=id), data = datapoly) +
	ggplot2::scale_fill_identity() 
} else {
  p <- p +
        ggplot2::geom_polygon(ggplot2::aes(x=x, y=y, fill=value, group=id), data = datapoly) +
        ggplot2::scale_fill_gradient(low=scale_col[1], 
				     high=scale_col[2], 
				     na.value="transparent", 
				     guide="colourbar",
				     limits=scale_lim,
				     name=var_units,
				     oob=scales::squish)
}
 p <- p + ggplot2::ggtitle(paste(var_longname, ds[day]))
print(p)
if (return_poly) {
  return(list(p=p, datapoly=datapoly, longitude=longitude, latitude=latitude))
} else {
  return(p)
}
}
