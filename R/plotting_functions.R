#' Calculate the number of days in the month of the specified date
#'
#' With thanks to G.Grothendieck on StackOverflow.
#'
#' @param d A date in the format given by as.Date
#' @return An integer
#' @export

daysIn <- function (d) {
	fom <- as.Date(cut(d, "month"))
	fom <- as.Date(cut(d, "month"))
	fomp1 <- as.Date(cut(fom+32, "month"))
	return(as.integer(fomp1 - fom))
}

#' Calculate and return the plume optical class from eReefs model output..
#'
#' Uses reflectance at multiple wavelengths from the model output to calculate plume colour
#' classes as defined by Devlin et al. (2012). Adapted from the Matlab function plume_detect.m 
#' by Mark Baird (CSIRO). This version by Barbara Robson (AIMS).
#'
#'
#' Devlin, M.J., McKinna, L.W., Álvarez-Romero, J.G., Petus, C., Abott, B., Harkness, P. and 
#'   Brodie, J., 2012. Mapping the pollutants in surface riverine flood plume waters in the Great 
#'   Barrier Reef, Australia. Marine pollution bulletin, 65(4-9), pp.224-235.
#'
#' @param rsr A list of 2D vectors containing the reflectances at various wavelengths
#'            extracted from an EMS netcdf file:
#'	      rsr <- list(R_412, R_443, R_488, R_531, R_547, R_667, R_678)
#' @return an array of the same size as R_412 containing values between 1 and 7 correspodning
#'         to optical plume classes.
#' @export
#' @examples
#' \dontrun{
#' plume_class(rsr)
#' }

plume_class <- function(rsr) {
    xdim <- nrow(rsr[[1]])
    ydim <- ncol(rsr[[1]])
    cl <- list(c(0.0064, 0.0093, 0.0147, 0.0242, 0.0286, 0.0245, 0.0240, 0.0050),
          c(0.0032, 0.0046, 0.0097, 0.0160, 0.0188, 0.0113, 0.0113, 0.0027),
          c(0.0031, 0.0043, 0.0092, 0.0151, 0.0178, 0.0105, 0.0105, 0.0025),
          c(0.0027, 0.0037, 0.0076, 0.0119, 0.0137, 0.0064, 0.0063, 0.0014),
          c(0.0040, 0.0045, 0.0064, 0.0065, 0.0062, 0.0013, 0.0012, 0.0002),
          c(0.0054, 0.0057, 0.0071, 0.0055, 0.0045, 0.0005, 0.0005, 0.0001),
          c(0.0130, 0.0110,  0.0070, 0.0040, 0.0030, 0.0010, 0.0010, 0.0010))

    rms <- array(0, dim=c(xdim, ydim, 7, 7))

    # Inefficient looped code. Vectorise this if it is too slow.
    for (i in 1:7) {
	# Can probably replace the j loop below with something along the lines of:
	# rms[,,i,] <- cl[[i]] - rsr # (not quite right)
	for (j in 1:7) {
	    rms[,,i,j] <- cl[[i]][j] - rsr[[j]]
	}
    }
    rms <- rms^2
    rmse <- rowSums(rms, dim=3)
    rmse[is.na(rmse)] <- 999999
    plume_class <- apply(rmse, c(1,2), which.min)
    plume_class[is.na(rsr[[1]])] <- NA
    #return(as.factor(plume_class))
    return(plume_class)
}

#' Create a surface map of eReefs model output.
#'
#' Creates a colour map showing concentrations of a specified eReefs model output variable at a specified
#' model layer (by default, the surface layer). The map is optionally (and by default) overlain on a Google 
#' Earth map of the region.
#' By Barbara Robson (AIMS).
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
#' @param layer Either a (positive) integer layer number, a negative number indicating depth below MSL (not depth below the free surface) 
#'        or 'surface' to choose the surface layer. Defaults to 'surface'.
#' @param Google_map_underlay Set to TRUE to use ggmap to show a Google Map as
#'      an underlay for the model output plot. Requires the ggmap library and an activated Google API key.
#'      Default now FALSE.
#' @param input_file is the URI or file location of any of the EMS output files, 
#'        Defaults to a menu selection. Set to "choices" to see some other pre-defined options that
#'        can be used (codenames as used in https://research.csiro.au/ereefs/models/model-outputs/access-to-raw-model-output/ )
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
#' @param label_towns Add labels for town locations to the figure. Default TRUE
#' @param strict_bounds TRUE to strictly enforce the box_bounds; FALSE for prettier edges conforming to x and y grids. Default FALSE.
#'        [dev note: An alternative would be to provide a version that plots from cell centres rather than using polygons, which should be fine other
#'         than near complex coastlines]
#' @export
#' @examples
#' \dontrun{
#' map_ereefs()
#' map_ereefs('TN')
#' map_ereefs('plume', target_date=c(2011, 1, 30))
#' map_ereefs('Chl_a_sum', target_date='2016-07-15', scale_col=c('ivory', 'green4'))
#' map_ereefs('salt', box_bounds=c(145,150,-20,-15), zoom=7, scale_lim=c(32,35))
#'}
map_ereefs <- function(var_name = "true_colour", 
                       target_date = c(2018,1,30), 
                       layer = 'surface', 
                       Google_map_underlay = FALSE,
                       input_file = "menu",
                       input_grid = NA,
                       scale_col = c('ivory', 'coral4'), 
                       scale_lim = c(NA, NA),
                       zoom = 6, 
                       box_bounds = c(NA, NA, NA, NA), 
                       p = NA, 
                       suppress_print = FALSE, 
                       return_poly = FALSE,
                       label_towns = TRUE,
                       strict_bounds = FALSE)
{

input_file <- substitute_filename(input_file)

# Check whether this is a locally-stored netcdf file or a web-served file
if (substr(input_file, 1, 4)=="http") {
  local_file = FALSE
} else local_file = TRUE

if (length(p)!=1) Google_map_underlay <- FALSE
if (suppress_print) Google_map_underlay <- FALSE

towns <- data.frame(latitude = c(-15.47027987, -16.92303816, -19.26639219, -20.0136699, -20.07670986, -20.40109791, -21.15345122, -22.82406858, -23.38031858, -23.84761069, -24.8662122, -25.54073075, -26.18916037),
                    longitude = c(145.2498605, 145.7662482, 146.805701, 148.2475387, 146.2635394, 148.5802016, 149.1655418, 147.6363616, 150.5059485, 151.256349, 152.3478987, 152.7049316, 152.6581893),
                    town = c('Cooktown', 'Cairns', 'Townsville', 'Bowen', 'Charters Towers', 'Prosperine', 'Mackay', 'Clermont', 'Rockhampton', 'Gladstone', 'Bundaberg', 'Maryborough', 'Gympie'))

# Check whether this is a GBR1 or GBR4 ereefs file, or something else
ereefs_case <- get_ereefs_case(input_file)
input_stem <- get_file_stem(input_file)
check_platform_ok(input_stem)
grids <- get_ereefs_grids(input_file, input_grid)
x_grid <- grids[['x_grid']]
y_grid <- grids[['y_grid']]

# Allow user to specify a depth below MSL by setting layer to a negative value
if (layer<=0) {
   z_grid <- grids[['z_grid']]
   layer <- which.min(z_grid<layer)
}


# Date to map:
if (is.vector(target_date)) {
	target_date <- as.Date(paste(target_date[1], target_date[2], target_date[3], sep='-'))
} else if (is.character(target_date)) {
	target_date <- as.Date(target_date)
}
if (ereefs_case == 4) { 
	filename <- paste0(input_stem, format(target_date, '%Y-%m'), '.nc')
	nc <- ncdf4::nc_open(filename)
	if (!is.null(nc$var[['t']])) { 
	    ds <- as.Date(safe_ncvar_get(nc, "t"), origin = as.Date("1990-01-01"))
        } else {
	    ds <- as.Date(safe_ncvar_get(nc, "time"), origin = as.Date("1990-01-01"))
	}
	day <- which.min(abs(target_date - ds))
	ncdf4::nc_close(nc)
} else if (ereefs_case == 1) {
	day <- 1
	ds <- target_date
	filename <- paste0(input_stem, format(target_date, '%Y-%m-%d'), '.nc')
} else {
	filename <- paste0(input_stem, '.nc')
	nc <- ncdf4::nc_open(filename)
	if (!is.null(nc$var[['t']])) { 
	    ds <- as.Date(safe_ncvar_get(nc, "t"), origin = as.Date("1990-01-01"))
        } else {
	    ds <- as.Date(safe_ncvar_get(nc, "time"), origin = as.Date("1990-01-01"))
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
    if (!local_file) {
      inputfile <- paste0(filename, '?R_412,R_443,R_488,R_531,R_547,R_667,R_678')
    } else inputfile <- filename
    nc <- ncdf4::nc_open(inputfile)
    R_412 <- safe_ncvar_get(nc, "R_412", start=c(xmin,ymin,day), count=c(xmax-xmin,ymax-ymin,1))
    R_443 <- safe_ncvar_get(nc, "R_443", start=c(xmin,ymin,day), count=c(xmax-xmin,ymax-ymin,1))
    R_488 <- safe_ncvar_get(nc, "R_488", start=c(xmin,ymin,day), count=c(xmax-xmin,ymax-ymin,1))
    R_531 <- safe_ncvar_get(nc, "R_531", start=c(xmin,ymin,day), count=c(xmax-xmin,ymax-ymin,1))
    R_547 <- safe_ncvar_get(nc, "R_547", start=c(xmin,ymin,day), count=c(xmax-xmin,ymax-ymin,1))
    R_667 <- safe_ncvar_get(nc, "R_667", start=c(xmin,ymin,day), count=c(xmax-xmin,ymax-ymin,1))
    R_678 <- safe_ncvar_get(nc, "R_678", start=c(xmin,ymin,day), count=c(xmax-xmin,ymax-ymin,1))
    rsr <- list(R_412, R_443, R_488, R_531, R_547, R_667, R_678)
    ems_var <- plume_class(rsr)
    dims <- dim(ems_var)
    var_units <- ''
    var_longname <- 'Plume colour class'

} else if (var_name=="true_colour") {
    if (!local_file) {
      inputfile <- paste0(filename, '?R_470,R_555,R_645')
    } else inputfile <- filename
    nc <- ncdf4::nc_open(inputfile)
    TCbright <- 10
    R_470 <- safe_ncvar_get(nc, "R_470", start=c(xmin,ymin,day), count=c(xmax-xmin,ymax-ymin,1)) * TCbright
    R_555 <- safe_ncvar_get(nc, "R_555", start=c(xmin,ymin,day), count=c(xmax-xmin,ymax-ymin,1)) * TCbright
    R_645 <- safe_ncvar_get(nc, "R_645", start=c(xmin,ymin,day), count=c(xmax-xmin,ymax-ymin,1)) * TCbright
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
    if (!local_file) {
      inputfile <- paste0(filename, '?ZooL_N,ZooS_N')
    } else inputfile <- filename
    nc <- ncdf4::nc_open(inputfile)
    # We don't yet know the dimensions of the variable, so let's get them
    dims <- nc$var[['ZooL_N']][['size']]
    if (is.null(dims)) stop(paste('ZooL_N', ' not found in netcdf file.')) 
    ndims <- length(dims)
    if ((ndims > 3) && (layer == 'surface')) layer <- dims[3]
    var_longname <- "Total zooplankton nitrogen"
    var_units <- "mg N m-3"
} else if (var_name == 'speed') {
    if (!local_file ) {
      inputfile <- paste0(filename, '?u,v')
    } else inputfile <- filename
    nc <- ncdf4::nc_open(inputfile)
    # We don't yet know the dimensions of the variable, so let's get them 
    dims <- nc$var[['u']][['size']] 
    if (is.null(dims)) stop(paste('u', ' not found in netcdf file.')) 
    ndims <- length(dims) 
    if ((ndims > 3) && (layer == 'surface')) layer <- dims[3]
    var_longname <- "Current speed"
    var_units <- "m s-1"
} else { 
    if (!local_file) {
      inputfile <- paste0(filename, '?', var_name)
    } else inputfile <- filename
    nc <- ncdf4::nc_open(inputfile)
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
       ems_var <- safe_ncvar_get(nc, 'ZooL_N', start=c(xmin,ymin,layer,day), count=c(xmax-xmin,ymax-ymin,1,1))
       ems_var <- ems_var + safe_ncvar_get(nc, 'ZooS_N', start=c(xmin,ymin,layer,day), count=c(xmax-xmin,ymax-ymin,1,1))
    } else {
       ems_var <- safe_ncvar_get(nc, 'ZooL_N', start=c(xmin,ymin,day), count=c(xmax-xmin,ymax-ymin,1))
       ems_var <- ems_var + safe_ncvar_get(nc, 'ZooS_N', start=c(xmin,ymin,day), count=c(xmax-xmin,ymax-ymin,1))
    }
} else if (var_name == 'speed') {
    var_longname <- 'Current speed'
    var_units <- 'm s-1'
    if (ndims == 4) {
       ems_var <- safe_ncvar_get(nc, 'u', start=c(xmin,ymin,layer,day), count=c(xmax-xmin,ymax-ymin,1,1))
       ems_var <- sqrt(ems_var^2 + safe_ncvar_get(nc, 'v', start=c(xmin,ymin,layer,day), count=c(xmax-xmin,ymax-ymin,1,1))^2)
    } else {
       ems_var <- safe_ncvar_get(nc, 'u', start=c(xmin,ymin,day), count=c(xmax-xmin,ymax-ymin,1))
       ems_var <- ems_var + safe_ncvar_get(nc, 'v', start=c(xmin,ymin,day), count=c(xmax-xmin,ymax-ymin,1))
    }
} else if (!((var_name == 'true_colour') || (var_name == 'plume'))) {
    vat <- ncdf4::ncatt_get(nc, var_name)
    var_longname <- vat$long_name
    var_units <- vat$units
    if (ndims == 4) {
       ems_var <- safe_ncvar_get(nc, var_name, start=c(xmin,ymin,layer,day), count=c(xmax-xmin,ymax-ymin,1,1))
    } else {
       ems_var <- safe_ncvar_get(nc, var_name, start=c(xmin,ymin,day), count=c(xmax-xmin,ymax-ymin,1))
    }
}

ncdf4::nc_close(nc)

nc <- ncdf4::nc_open(filename)
if (!is.null(nc$var[['botz']])) {
   botz <- safe_ncvar_get(nc, 'botz', start=c(xmin,ymin), count=c(xmax-xmin,ymax-ymin))
} else {
   botz <- array(NA, dim(ems_var))
}
if (is.null(nc$var[['latitude']])) {
# Standard EMS output file
  latitude <- safe_ncvar_get(nc, "y_centre")
  longitude <- safe_ncvar_get(nc, "x_centre")
} else { 
  # Simple format netcdf file
  latitude <- safe_ncvar_get(nc, "latitude")
  longitude <- safe_ncvar_get(nc, "longitude")
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
botz <- c(botz)[gx_ok&gy_ok]

# Unique ID for each polygon
id <- 1:length(n)

id <- as.factor(id)
values <- data.frame(id = id, value = n, depth=botz)
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
  myMap<-suppressWarnings(ggmap::get_map(location=MapLocation, source="google", maptype="satellite", zoom=zoom, crop=TRUE, scale=2))
  p <- ggmap::ggmap(myMap)
} else if (length(p)==1) {
  p <- ggplot2::ggplot()
}
if (var_name=="true_colour") {
  p <- p +
        ggplot2::geom_polygon(ggplot2::aes(x=x, y=y, fill=value, group=id), data = datapoly) +
	ggplot2::scale_fill_identity() 
} else {
   if (length(scale_col)==1) scale_col <- c('ivory', scale_col)
   p <- p + ggplot2::geom_polygon(ggplot2::aes(x=x, y=y, fill=value, group=id), data = datapoly)
   if (length(scale_col)==2) { 
      p <- p +
           ggplot2::scale_fill_gradient(low=scale_col[1], 
                                        high=scale_col[2], 
				                            na.value="transparent", 
				                            guide="colourbar",
				                            limits=scale_lim,
				                            #name=expression(paste('cells m'^'-3')),
				                            name=var_units,
				                            oob=scales::squish)
        } else { 
           p <- p +
           ggplot2::scale_fill_gradient2(low=scale_col[1], 
                                         mid=scale_col[2], 
                                         high=scale_col[3], 
                                         na.value="transparent", 
                                         guide="colourbar",
				                             limits=scale_lim,
				                             #name=expression(paste('cells m'^'-3')),
				                             name=var_units,
				                             oob=scales::squish)
        }
}

if (label_towns) {
   towns <- towns[towns$latitude>=min(gy, na.rm=TRUE),]
   towns <- towns[towns$latitude<=max(gy, na.rm=TRUE),]
   towns <- towns[towns$longitude>=min(gx, na.rm=TRUE),]
   towns <- towns[towns$longitude<=max(gx, na.rm=TRUE),]
   if (dim(towns)[1]>0) p <- p + ggplot2::geom_label(data=towns, ggplot2::aes(x=longitude, y=latitude, label=town, hjust="right"))
}

p <- p + ggplot2::ggtitle(paste(var_longname, ds[day])) +
    ggplot2::xlab('longitude') + ggplot2::ylab('latitude')
if (strict_bounds) {
  if (is.na(box_bounds[1])) box_bounds[1] <- min(positions$x)
  if (is.na(box_bounds[2])) box_bounds[2] <- max(positions$x)
  if (is.na(box_bounds[3])) box_bounds[3] <- min(positions$y)
  if (is.na(box_bounds[4])) box_bounds[4] <- max(positions$y)
  p <- p + ggplot2::xlim(box_bounds[1],box_bounds[2])+ggplot2::ylim(box_bounds[3],box_bounds[4])
}
if (!suppress_print) print(p)
if (return_poly) {
  return(list(p=p, datapoly=datapoly, longitude=longitude, latitude=latitude))
} else {
  return(p)
}
}

#' Create a series of map image files for an animation of eReefs model output AND calculate temporal
#' mean values.
#'
#' Creates and saves to disk a sequential series of colour map images showing concentrations of a specified 
#' eReefs model output variable at a specified model layer (by default, the surface layer). 
#' Also calculates the temporal mean value of each cell over the specified time (visualisation of maps can
#' be suppressed by setting suppress_print to TRUE if this is the primary desired output).
#' Maps produced are optionally (and by default) overlain on a Google Earth map of the region. 
#' Can be more efficient than calling map_ereefs multiple times if you 
#' want to produce an animation because it loads a month at a time for GBR4 runs. 
#' If output files contain multiple outputs per day, chooses the step closest to midday and uses only daily output.
#' To stitch together the images into an animation, you will need other software such as ImageMagick (recommended)
#' or imageJ.  Barbara Robson (AIMS).
#'
#' @param var_name Short name of the variable to plot. This can be any variable in the
#'                 eReefs netcdf file that you are accessing (refer to eReefs model
#'                 documentation or extract variable names from the netcdf file for a full 
#'                 list). In addition, two special var_name values are supported: 'true_colour'
#'                 and 'plume'. 'true_colour' provides a simulated MODIS true colour
#'                 image (refer to Baird et al., 2016 for en explanation).'plume'
#'                 provides a map of calculated plume colour class as per Devlin et al. (2012)
#'                 and Baird et al. (2017). Defaults to true_colour.
#' @param start_date date  for animation. Can either be a) a vector of integers in the format 
#'	c(year, month, day); b) a date obtained e.g. from as.Date(); or c) a character string 
#'      formatted for input to as.Date(). Defaults to c(2015,2,1).
#' @param end date for animation. Can either be a) a vector of integers in the format 
#'	c(year, month, day); b) a date obtained e.g. from as.Date(); or c) a character string 
#'      formatted for input to as.Date(). Defaults to c(2016,3,31).
#' @param layer Either an integer layer number or 'surface' to choose the surface layer. Defaults to 'surface'.
#' @param output_dir Path to directory in which to store images for animation. Created if necessary. Defaults
#'      to 'ToAnimate'. Images are created in this directory with filenames beginning with var_name, 
#'      followed by an underscore and then sequential numbers beginning with 100001.
#' @param Google_map_underlay Set to TRUE to use ggmap to show a Google Map as
#'      an underlay for the model output plot. Requires the ggmap librray and an activated Google API key.
#'      Default now FALSE.
#' @param input_file is the URI or file location of any of the EMS output files, 
#'        Defaults to a menu selection. Set to "choices" to see some other pre-defined options that
#'        can be used (codenames as used in https://research.csiro.au/ereefs/models/model-outputs/access-to-raw-model-output/ )
#' @param input_grid Name of the locally-stored or opendap-served netcdf file that contains the grid
#'      coordinates for the cell corners (x_grid and y_grid). If not specified, the function will first look for
#'      x_grid and y_grid can be found in the first INPUT_STEM file, and if not found, will check whether the size 
#'      of the variables in the input file corresponds to the size expected for GBR4 or GBR1, and load appropriate 
#'      x and y grids from data files stored in this package. Alternatively, you can provide the location of a full 
#'      (not simple-format) ereefs netcdf output file such as 
#'      "http://dapds00.nci.org.au/thredds/dodsC/fx3/gbr4_hydro_all/gbr4_all_2016-09.nc"
#' @param scale_col Vector of colours to use for the colour scale. This can be colours 
#'      from the ggplot colour palette or a RGB hash code. Ignored for true_colour plots. 
#'      If one value is given, low colour is set to ivory and high colour to the value given.
#'      If two values are given, these are used as low and high limit colours.
#'      If three values are given, the middle value is used to set the mid-point of the scale.
#'      Defaults to c('ivory', 'coral4').
#' @param zoom Value of zoom passed to ggmap::ggmap(). Set to 5 if you want to show the entire extent 
#'      of eReefs models. Defaults to 6. Higher values will zoom in further.
#' @param box_bounds Minimum and maximum latitude and longitude coordinates to map. Defaults to the
#'        entire extent of the model output (though modified by the value of zoom). 
#'        Format: c(longitude_min, longitude_max, latitude_min, latitude_max). It is recommended to
#'        also specify an appropriate value for zoom if specifying box_bounds.
#' @param suppress_print Set to TRUE if you don't want the plots generatedand saved. Defaults to FALSE.
#' @param stride Default 'daily', but can otherwise be set to a numeric interval indicating how many time-steps to step forward for each frame.
#' @param verbosity Set 0 for just a waitbar, 1 for more updates, 2 for debugging information. Default 0.
#' @param label_towns Add labels for town locations to the figure. Default TRUE
#' @param strict_bounds TRUE to strictly enforce the box_bounds; FALSE for prettier edges conforming to x and y grids. Default FALSE.
#' @return a data.frame formatted for use in ggplot2::geom_polygon, containing a map of the temporally averaged
#'       value of the variable specified in VAR_NAME over the selected interval.
#' @export
#' @examples
#' \dontrun{
#' map_ereefs_movie(start_date=c(2016,2,1),end_date=c(2016,2,15))
#'}
map_ereefs_movie <- function(var_name = "true_colour", 
                             start_date = c(2015,12,1), 
                             end_date = c(2016,3,31), 
                             layer = 'surface', 
                             output_dir = 'ToAnimate', 
                             Google_map_underlay = FALSE, 
                             input_file = "menu",
                             input_grid = NA, 
                             scale_col = c('ivory', 'coral4'), 
                             scale_lim = c(NA, NA), 
                             zoom = 6, 
                             box_bounds = c(NA, NA, NA, NA), 
                             suppress_print=FALSE, 
                             stride = 'daily',
                             verbosity=0, 
                             label_towns = TRUE,
                             strict_bounds = FALSE)
{
  input_file <- substitute_filename(input_file)

  # Check whether this is a locally-stored netcdf file or a web-served file
  if (substr(input_file, 1, 4)=="http") {
    local_file = FALSE
  } else local_file = TRUE

  # Check whether this is a GBR1 or GBR4 ereefs file, or something else
  ereefs_case <- get_ereefs_case(input_file) 
  if (ereefs_case==1) warning('Assuming that only one timestep is output per day/file') # find matching commented warning to fix this
  input_stem <- get_file_stem(input_file)
  check_platform_ok(input_stem)

  towns <- data.frame(latitude = c(-15.47027987, -16.92303816, -19.26639219, -20.0136699, -20.07670986, -20.40109791, -21.15345122, -22.82406858, -23.38031858, -23.84761069, -24.8662122, -25.54073075, -26.18916037),
                      longitude = c(145.2498605, 145.7662482, 146.805701, 148.2475387, 146.2635394, 148.5802016, 149.1655418, 147.6363616, 150.5059485, 151.256349, 152.3478987, 152.7049316, 152.6581893),
                      town = c('Cooktown', 'Cairns', 'Townsville', 'Bowen', 'Charters Towers', 'Prosperine', 'Mackay', 'Clermont', 'Rockhampton', 'Gladstone', 'Bundaberg', 'Maryborough', 'Gympie'))
  grids <- get_ereefs_grids(input_file, input_grid)
  x_grid <- grids[['x_grid']]
  y_grid <- grids[['y_grid']]
# Allow user to specify a depth below MSL by setting layer to a negative value
if (layer<=0) {
   z_grid <- grids[['z_grid']]
   layer <- which.min(z_grid<layer)
}


  # Dates to map:
  if (is.vector(start_date)) {
	  start_date <- as.Date(paste(start_date[1], start_date[2], start_date[3], sep='-'))
  } else if (is.character(start_date)) {
	  start_date <- as.Date(start_date)
  }
  start_day <- as.integer(format(start_date, '%d'))
  start_month <- as.integer(format(start_date, '%m'))
  start_year <- as.integer(format(start_date, '%Y'))
  
  if (is.vector(end_date)) {
	end_date <- as.Date(paste(end_date[1], end_date[2], end_date[3], sep='-'))
  } else if (is.character(end_date)) {
	end_date <- as.Date(end_date)
  }
  end_day <- as.integer(format(end_date, '%d'))
  end_month <- as.integer(format(end_date, '%m'))
  end_year <- as.integer(format(end_date, '%Y'))

  if (start_date > end_date) {
    stop('start_date must preceed end_date')
  }

  # Check for ggmap::ggmap()
  if ((Google_map_underlay)&(!requireNamespace("ggmap", quietly = TRUE))) {
    warning('Package ggmap::ggmap is required to show a Google map underlay. Preparing plot with no underlay.')
    print('To avoid this message, either install ggmap or set Google_map_underlay = FALSE.')
    Google_map_underlay <- FALSE
  }
  if (suppress_print) Google_map_underlay <- FALSE
  
  if (start_year==end_year) {
      mths <- start_month:end_month
      years <- rep(start_year, length(mths))
  } else if ((start_year + 1) == end_year) {
      mths <- c(start_month:12, 1:end_month)
      years <- c(rep(start_year, 12 - start_month + 1), rep(end_year, end_month))
  } else {
      mths <- c(start_month:12, rep(1:12, end_year - start_year - 1), 1:end_month)
      years <- c(rep(start_year, 12 - start_month + 1), 
                 rep((start_year + 1) : (end_year - 1), each=12),
                 rep(end_year, end_month))
  }


  # Allow for US English:
  if (var_name == "true_color") {
	var_name <- "true_colour"
  }
  
  dims <- dim(x_grid) - 1

  # Work out which parts of the grid are within box_bounds and which are outside
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
         
  # Find the subset of x_grid and y_grid that is inside the box and crop the grids
  # to the box_bounds
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

  # Set up the polygon corners. 4 per polygon.
  a <- xmax - xmin + 1
  b <- ymax - ymin + 1

  gx <- c(x_grid[1:a, 1:b], x_grid[2:(a+1), 1:b], x_grid[2:(a+1), 2:(b+1)], x_grid[1:a, 2:(b+1)])
  gy <- c(y_grid[1:a, 1:b], y_grid[2:(a+1), 1:b], y_grid[2:(a+1), 2:(b+1)], y_grid[1:a, 2:(b+1)])
  gx <- array(gx, dim=c(a*b,4))
  gy <- array(gy, dim=c(a*b,4))

  # Find and exclude points where not all corners are defined
  gx_ok <- !apply(is.na(gx),1, any)
  gy_ok <- !apply(is.na(gy),1, any)
  gx <- c(t(gx[gx_ok&gy_ok,]))
  gy <- c(t(gy[gx_ok&gy_ok,]))
  if (Google_map_underlay) { 
	  MapLocation<-c(min(x_grid, na.rm=TRUE)-0.5, 
                    min(y_grid, na.rm=TRUE)-0.5, 
                    max(x_grid, na.rm=TRUE)+0.5, 
                    max(y_grid, na.rm=TRUE)+0.5) 
     myMap<-suppressWarnings(ggmap::get_map(location=MapLocation, source="google", maptype="satellite", crop=TRUE, zoom=zoom, scale=2))
     #myMap<-suppressWarnings(ggmap::get_map(location=MapLocation, source="stamen", maptype="watercolor"))
  }

  # Main routine
  ndims <- 0
  icount <- 0
  mcount <- 0
  pb <- txtProgressBar(min = 0, max = 1, style = 3)
  for (month in mths) {
    mcount <- mcount + 1
    year <- years[mcount]

    if (mcount == 1) {
       from_day <- start_day
       if (ereefs_case == 0) {
	      filename <- paste0(input_stem, '.nc')
	      nc <- ncdf4::nc_open(filename)
	      if (!is.null(nc$var[['t']])) { 
	          ds <- as.Date(safe_ncvar_get(nc, "t"), origin = as.Date("1990-01-01"))
              } else {
	          ds <- as.Date(safe_ncvar_get(nc, "time"), origin = as.Date("1990-01-01"))
	      }
	      ncdf4::nc_close(nc)
       }
    } else {
       from_day <- 1
    }
    if (mcount == (length(mths))) {
       day_count <- end_day
    } else if (mcount == 1) {
       day_count <- daysIn(as.Date(paste(year, month, 1, sep='-'))) - start_day + 1
    } else if ((start_year==end_year)&&(start_month==end_month)) {
       day_count <- end_day - start_day
    } else {
       day_count <- daysIn(as.Date(paste(year, month, 1, sep='-')))
    }

    if (ereefs_case == 4) { 
	    filename <- paste0(input_stem, format(as.Date(paste(year, month, 1, sep="-")), '%Y-%m'), '.nc')
	    nc <- ncdf4::nc_open(filename)
	    if (!is.null(nc$var[['t']])) { 
	        ds <- as.Date(safe_ncvar_get(nc, "t"), origin = as.Date("1990-01-01"))
            } else {
	        ds <- as.Date(safe_ncvar_get(nc, "time"), origin = as.Date("1990-01-01"))
	    }
	    ncdf4::nc_close(nc)
       if ((ds[length(ds)] - ds[1]) > 31.5) warning('Filename looks like a monthly output file (i.e. contains two dashes) but file contains more than a month of data.')
       if ((ds[length(ds)] - ds[1]) < 27) warning('Filename looks like a monthly output file (i.e. contains two dashes) but file contains less than a month of data.')
       if(ds[2]==ds[1]) stop(paste('Error reading time from', filename, '(t[2]==t[1])'))
       tstep <- as.numeric(ds[2]-ds[1])
       dum1 <- as.integer((from_day - 0.4999999)/tstep + 1)
       dum2 <- as.integer((day_count - 1) / tstep) +1
       ds <- ds[seq(from=dum1, by=as.integer(1/tstep), to=(dum1+dum2))]
       start_array <- c(xmin, ymin, dum1)
	    count_array <- c(xmax-xmin, ymax-ymin, dum2)
	    fileslist <- 1
    } else if (ereefs_case == 1) { 
       #warning('Assuming that only one timestep is output per day/file')
	    fileslist <- from_day:(from_day+day_count-1)
       start_array <- c(xmin, ymin, 1) 
       count_array <- c(xmax-xmin, ymax-ymin, 1)
       tstep <- 1
    } else {
	    # Everything is in one file but we are only going to read a month at a time
	    # Output may be more than daily, or possibly less
	    # filename <- paste0(input_stem, '.nc') # filename has been set previously 
       tstep <- as.numeric(ds[2]-ds[1])
       from_day <- as.integer((as.Date(paste(year, month, from_day, sep="-")) - ds[1])/tstep) + 1
	    if (from_day<1) from_day <-1
	    start_array <- c(xmin, ymin, from_day)
	    count_array <- c(xmax-xmin, ymax-ymin, as.integer(day_count/tstep))
	    fileslist <- 1
    }
    if (stride == 'daily') {
       stride <- 1/tstep
    } 
    stride <- as.integer(stride)

    for (i in fileslist) {
      if (ereefs_case == 1) { 
         filename <- paste0(input_stem, format(as.Date(paste(year, month, i, sep="-")), '%Y-%m-%d'), '.nc') 
         ds <- as.Date(paste(year, month, i, sep="-", '%Y-%m-%d'))
      }
      if (verbosity>0) print(filename)
      if (var_name=="plume") {
        if (!local_file) {
          slice <- paste0('[', start_array[3]-1, ':', stride, ':', start_array[3] + count_array[3] - 2, ']', # time
                          '[', start_array[2]-1, ':', start_array[2] + count_array[2] - 1, ']', # y
                          '[', start_array[1]-1, ':', start_array[1] + count_array[1] - 1, ']') # x
          inputfile <- paste0(filename, '?R_412', slice, ',R_443', slice, ',R_488', slice, ',R_531', slice, ',R_547', slice, ',R_667', slice, ',R_678', slice)
        } else inputfile <- filename
        nc <- ncdf4::nc_open(inputfile)
        R_412 <- safe_ncvar_get(nc, "R_412")
        R_443 <- safe_ncvar_get(nc, "R_443")
        R_488 <- safe_ncvar_get(nc, "R_488")
        R_531 <- safe_ncvar_get(nc, "R_531")
        R_547 <- safe_ncvar_get(nc, "R_547")
        R_667 <- safe_ncvar_get(nc, "R_667")
        R_678 <- safe_ncvar_get(nc, "R_678")
        if (local_file) {
          R_412 <- R_412[(start_array[1] - 1) : (start_array[1] + count_array[1] - 1),
                         (start_array[2] - 1) : (start_array[2] + count_array[2] - 1),
                         seq(from = start_array[3] - 1, to = start_array[3] + count_array[3] - 2, by = stride)]
        }
        ems_var <- NA*R_678
        if (ereefs_case ==4) {
            for (day in 1:dim(R_412)[3]) {
              rsr <- list(R_412[,,day], R_443[,,day], R_488[,,day], R_531[,,day], R_547[,,day], R_667[,,day], R_678[,,day])
              ems_var[,,day] <- plume_class(rsr)
            }
        } else {
	     rsr <- list(R_412, R_433, R_488, R_532, R_547, R_668, R_678)
	     ems_var <- plume_class(rsr)
        }
        dims <- dim(ems_var)
	     var_longname <- "Plume optical class"
	     var_units <- ""
      } else if (var_name=="true_colour") {
        slice <- paste0('[', start_array[3]-1, ':', stride, ':', start_array[3] + count_array[3] - 2, ']', # time
                        '[', start_array[2]-1, ':', start_array[2] + count_array[2] - 1, ']', # y
                        '[', start_array[1]-1, ':', start_array[1] + count_array[1] - 1, ']') # x
        if (!local_file) {
          inputfile <- paste0(filename, '?R_470', slice, ',R_555', slice, ',R_645', slice)
        } else inputfile <- filename
        nc <- ncdf4::nc_open(inputfile)
        TCbright <- 10
        R_470 <- safe_ncvar_get(nc, "R_470") * TCbright
        R_555 <- safe_ncvar_get(nc, "R_555") * TCbright
        R_645 <- safe_ncvar_get(nc, "R_645") * TCbright
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
    
      } else { 

        if (ndims == 0) {
          nc <- ncdf4::nc_open(filename)
          # We don't yet know the dimensions of the variable, so let's get them
          if (var_name == "speed") {
             dims <- nc$var[['u']][['size']]
          } else if (var_name == "ZooT") {
             dims <- nc$var[['ZooL_N']][['size']]
          } else { 
             dims <- nc$var[[var_name]][['size']]
          }
          if (is.null(dims)) stop(paste(var_name, ' not found in netcdf file.')) 
          ndims <- length(dims)
          if ((ndims > 3) && (layer == 'surface')) layer <- dims[3]
          ncdf4::nc_close(nc)
        }

        if (verbosity>1) {
           print(paste('variable has ', ndims, 'dimensions'))
           print('start_array = ')
           print(start_array)
           print('count_array = ')
           print(count_array)
        }

        if (ndims == 4) {
          slice <- paste0('[', start_array[3]-1, ':', stride, ':', start_array[3] + count_array[3] - 2, ']', # time
                          '[', layer-1, ']',                                                    # layer
                          '[', start_array[2]-1, ':', start_array[2] + count_array[2] - 1, ']', # y
                          '[', start_array[1]-1, ':', start_array[1] + count_array[1] - 1, ']') # x
        } else {
          slice <- paste0('[', start_array[3]-1, ':', stride, ':', start_array[3] + count_array[3] - 2, ']', # time
                          '[', start_array[2]-1, ':', start_array[2] + count_array[2] - 1, ']', # y
                          '[', start_array[1]-1, ':', start_array[1] + count_array[1] - 1, ']') # x
        }
        if (verbosity>1) print(paste('slice = ', slice))
        if (var_name == "speed") { 
           if (!local_file) {
             inputfile <- paste0(filename, '?u', slice, ',v', slice)
           } else inputfile <- filename
           nc <- ncdf4::nc_open(inputfile)
           ems_var <- sqrt(safe_ncvar_get(nc, 'u')^2 + ncdf4::ncvar_get(nc, 'v')^2)
           vat <- ncdf4::ncatt_get(nc, 'u')
           var_longname <- 'Current speed'
           var_units <- vat$units
        } else if (var_name == "ZooT") {
           if (!local_file) {
             inputfile <- paste0(filename, '?ZooL_N', slice, ',ZooS_N', slice)
           } else inputfile <- filename
           nc <- ncdf4::nc_open(inputfile)
           ems_var <- safe_ncvar_get(nc, 'ZooL_N') + ncdf4::ncvar_get(nc, 'ZooS_N')
           vat <- ncdf4::ncatt_get(nc, 'ZooL_N')
           var_longname <- 'Total Zooplankton'
           var_units <- vat$units
        } else {
           if (!local_file) {
             inputfile <- paste0(filename, '?',var_name, slice)
           } else inputfile <- filename
           nc <- ncdf4::nc_open(inputfile)
           ems_var <- safe_ncvar_get(nc, var_name)
           vat <- ncdf4::ncatt_get(nc, var_name)
           var_longname <- vat$long_name
           var_units <- vat$units
        }
        if (local_file) { 
          if (ndims == 4) {
            ems_var <- ems_var[start_array[1] : (start_array[1] + count_array[1]),
                               start_array[2] : (start_array[2] + count_array[2]),
                               seq(from = start_array[3], to = start_array[3] + count_array[3] - 1, by = stride)] 
          }  else {
            ems_var <- ems_var[start_array[1] : (start_array[1] + count_array[1]),
                               start_array[2] : (start_array[2] + count_array[2]),
                               layer,
                               seq(from = start_array[3], to = start_array[3] + count_array[3] - 1, by = stride)] 
          }
      }
      #ds <- as.Date(safe_ncvar_get(nc, "time"), origin = as.Date("1990-01-01"))
    
      dum1 <- length(dim(ems_var))
      if (dum1==2) {
        ems_var <- array(ems_var, dim=c(dim(ems_var), 1))
      }
    
      ncdf4::nc_close(nc)
    
      d <- dim(ems_var)[3]
      for (jcount in 1:d) {
        ems_var2d <- ems_var[, , jcount]
        # Values associated with each polygon at chosen timestep
        n <- c(ems_var2d)[gx_ok&gy_ok]
        if (icount==0) {
           if (var_name=='true_colour') {
              temporal_sum <- as.integer(stringi::stri_sub(n,2,6))
           } else { 
              temporal_sum <- n
           }
        } else {
           if (var_name=='true_colour') {
              temporal_sum <- as.integer(stringi::stri_sub(n,2,6)) + temporal_sum
           } else { 
              temporal_sum <- temporal_sum + n
           }
        }
    
        # Unique ID for each polygon
        id <- as.factor(1:length(n))
        values <- data.frame(id = id, value = n)
        positions <- data.frame(id=rep(id, each=4), x = gx, y = gy)
        datapoly <- merge(values, positions, by = c("id"))
        #print('debug 1'); print(var_name)
    
        if (!suppress_print) {
            if ((var_name!="true_colour")&&(is.na(scale_lim[1]))) { 
	            scale_lim <- c(min(n, na.rm=TRUE), max(n, na.rm=TRUE))
            }
  
            if (Google_map_underlay) {
               p <- ggmap::ggmap(myMap)
	         } else {
	            p <- ggplot2::ggplot()
            }
            if (var_name=="true_colour") {
        #print('debug 2');
	            p <- p +
               ggplot2::geom_polygon(ggplot2::aes(x=x, y=y, fill=value, group=id), data = datapoly) +
	            ggplot2::scale_fill_identity()
            } else {
               p <- p +
               ggplot2::geom_polygon(ggplot2::aes(x=x, y=y, fill=value, group=id), data = datapoly) +
               if (length(scale_col)==1) scale_col <- c('ivory', scale_col)
               if (length(scale_col)==2) { 
                  p <- p +
                  ggplot2::scale_fill_gradient(low=scale_col[1],
					                                high=scale_col[2],
					                                na.value="transparent", 
					                                guide="colourbar",
					                                limits=scale_lim,
					                                name=var_units,
					                                oob=scales::squish)
               } else {
                  p <- p +
                  ggplot2::scale_fill_gradient2(low=scale_col[1],
                                                mid=scale_col[2],
					                                high=scale_col[3],
					                                na.value="transparent", 
					                                guide="colourbar",
					                                limits=scale_lim,
					                                name=var_units,
					                                oob=scales::squish)
               }
            }

            if (label_towns) {
              towns <- towns[towns$latitude>=min(gy, na.rm=TRUE),]
              towns <- towns[towns$latitude<=max(gy, na.rm=TRUE),]
              towns <- towns[towns$longitude>=min(gx, na.rm=TRUE),]
              towns <- towns[towns$longitude<=max(gx, na.rm=TRUE),]
              if (dim(towns)[1]>0) p <- p + ggplot2::geom_label(data=towns, ggplot2::aes(x=longitude, y=latitude, label=town, hjust="right"))
            }
            p <- p + ggplot2::ggtitle(paste(var_longname, ds[jcount]))
            p <- p + ggplot2::xlab("longitude") + ggplot2::ylab("latitude") 
            if (strict_bounds) {
              if (is.na(box_bounds[1])) box_bounds[1] <- min(positions$x)
              if (is.na(box_bounds[2])) box_bounds[2] <- max(positions$x)
              if (is.na(box_bounds[3])) box_bounds[3] <- min(positions$y)
              if (is.na(box_bounds[4])) box_bounds[4] <- max(positions$y)
              p <- p + ggplot2::xlim(box_bounds[1],box_bounds[2])+ggplot2::ylim(box_bounds[3],box_bounds[4])
            }
            icount <- icount + 1
            #print(p)
            if (!file.exists(output_dir)) {
               dir.create(output_dir)
            }
            fname <- paste0(output_dir, '/', var_name, '_', 100000 + icount, '.png', collapse='')
            ggplot2::ggsave(fname, p)
            rm('p')
         }  else {
            icount <- icount + 1
         }
         setTxtProgressBar(pb,icount/as.integer(end_date-start_date)/tstep*stride)
      }
    }
  }
  close(pb)
  values <- data.frame(id = id, value = temporal_sum/icount)
  datapoly <- merge(values, positions, by = c("id"))
  return(list(datapoly=datapoly, value=values$value))
}

#' Create a map using a dataframe in the format required by ggplot2::geom_plot, for instance from map_ereefs() or map_ereefs_movie()
#'
#' Plots a map figure in the same format as would be given by map_ereefs(), but using a pre-generated dataframe, instead of
#' processing data directly from ereefs netcdf files.
#' 
#' @param datapoly A dataframe in the format required by geom_plot(), as provided by map_ereefs() or map_ereefs_movie().
#' @param var_longname Character vector to use for the figure title.
#' @param var_units Units to include in the figure labelling.
#' @param Google_map_underlay Set to TRUE to use ggmap to show a Google Map as
#'      an underlay for the model output plot. Requires the ggmap library and an activated Google API key.
#'      Default now FALSE.
#' @param scale_col Vector of colours to use for the colour scale. This can be colours 
#'      from the ggplot colour palette or a RGB hash code. Ignored for true_colour plots. 
#'      If one value is given, low colour is set to ivory and high colour to the value given.
#'      If two values are given, these are used as low and high limit colours.
#'      If three values are given, the middle value is used to set the mid-point of the scale.
#'      Defaults to c('ivory', 'coral4').
#' @param scale_lim Upper and lower bounds for colour scale. Defaults to full range of data.
#'      Ignored for true_colour plots.
#' @param zoom Value of zoom passed to ggmap(). Set to 5 if you want to show the entire extent 
#'      of eReefs models. Defaults to 6. Higher values will zoom in further.
#' @param suppress_print Default FALSE. If true, don't prdocue the map image.
#' @param p Handle for an existing figure if you want to add a layer instead of creating a new figure.
#'        If p is provided, Google_map_underlay is over-ridden and set to FALSE.
#' @return p Handle for the figure generated.
#' @export
#' @examples
#' \dontrun{
#' a <- map_ereefs('TN', return_poly=TRUE)
#' plot_map(a[[2]])
#'}

plot_map <- function(datapoly,
             var_longname = '',
             var_units = '',
		       Google_map_underlay = FALSE,
             scale_col = c('ivory', 'coral4'),
		       scale_lim = c(NA, NA),
             zoom = 6,
		       p = NA,
             suppress_print = FALSE)
{
  if (class(datapoly$value)=="factor") {
     var_name <- "true_colour"
  } else {
     var_name <- "something else"
  }
  if ((var_name!="true_colour")&&(is.na(scale_lim[1]))) { 
	  scale_lim <- c(min(datapoly$value, na.rm=TRUE), max(datapoly$value, na.rm=TRUE))
  }

  if (suppress_print) Google_map_underlay <- FALSE
  if (length(p)!=1) Google_map_underlay <- FALSE
  if (Google_map_underlay) {
    MapLocation<-c(min(datapoly$x, na.rm=TRUE)-0.5, 
 		min(datapoly$y, na.rm=TRUE)-0.5, 
 		max(datapoly$x, na.rm=TRUE)+0.5, 
 		max(datapoly$y, na.rm=TRUE)+0.5)
    myMap<-suppressWarnings(ggmap::get_map(location=MapLocation, source="google", maptype="satellite", zoom=zoom, crop=TRUE, scale=2))
    p <- ggmap::ggmap(myMap)
  } else if (length(p)==1) {
    p <- ggplot2::ggplot()
  }
  if (!suppress_print) {
  if (var_name=="true_colour") {
    p <- p +
        ggplot2::geom_polygon(ggplot2::aes(x=x, y=y, fill=value, group=id), data = datapoly) +
	     ggplot2::scale_fill_identity() 
  } else {
    p <- p +
        ggplot2::geom_polygon(ggplot2::aes(x=x, y=y, fill=value, group=id), data = datapoly)
        if (length(scale_col)==1) scale_col <- c('ivory', scale_col)
        if (length(scale_col)==2) { 
           p <- p +
           ggplot2::scale_fill_gradient(low=scale_col[1], 
				                             high=scale_col[2], 
				                             na.value="transparent", 
				                             guide="colourbar",
				                             limits=scale_lim,
				                             name=var_units,
				                             oob=scales::squish)
        } else {
           p <- p +
           ggplot2::scale_fill_gradient2(low=scale_col[1], 
                                         mid=scale_col[2],
				                             high=scale_col[3], 
				                             na.value="transparent", 
				                             guide="colourbar",
				                             limits=scale_lim,
				                             name=var_units,
				                             oob=scales::squish)
        }
  }
  p <- p + ggplot2::ggtitle(var_longname)
  print(p)
  }
  return(p)
}
