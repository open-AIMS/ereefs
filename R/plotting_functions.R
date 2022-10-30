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
#' model layer (by default, the surface layer). The map is optionally (and by default) overlain on a map
#' of Queensland (if off the Queensland coast).
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
#' @param Land_map Set to TRUE to show a map of Queensland as an underlay for the model output plot. No longer requires the fgmap library 
#'      and an activated Google API key but also doesn't show maps for other locations
#'      Default now FALSE.
#' @param input_file is the URL or file location of any of the EMS output files or a THREDDS catalog URI. 
#'        Defaults to a menu selection based on current NCI catalogs. Can also be set to "nci", "menu" or "catalog" for the same behaviour.
#'        Set to "old_menu" to provide old menu options instead of menu options from the NCI catalog.
#'        Numeric values are interpreted as references to selections available from the old menu.
#'        Short codes can be used for some options (codenames as used in https://research.csiro.au/ereefs/models/model-outputs/access-to-raw-model-output/ )
#' @param input_grid Name of the locally-stored or opendap-served netcdf file that contains the grid
#'      coordinates for the cell corners (x_grid and y_grid). If not specified, the function will first look for
#'      x_grid and y_grid can be found in the first INPUT_STEM file, and if not found, will check whether the size 
#'      of the variables in the input file corresponds to the size expected for GBR4 or GBR1, and load appropriate 
#'      x and y grids from data files stored in this package. Alternatively, you can provide the location of a full 
#'      (not simple-format) ereefs netcdf output file such as 
#'      "http://dapds00.nci.org.au/thredds/dodsC/fx3/gbr4_hydro_all/gbr4_all_2016-09.nc"
#' @param scale_col Vector of colours to use for low and high values in the colour scale. This can be a colour 
#'      from the ggplot colour palette or a RGB hash code, or "spectral". Ignored for true_colour plots. 
#'      If one value is given (other than "spectral"), low colour is set to ivory and high colour to the value given.
#'      If three values are given, uses scale_fill_gradient2 (spectrum from low to high through middle value).
#'      Defaults to c('ivory', 'coral4').
#' @param scale_lim Upper and lower bounds for colour scale. Defaults to full range of data.
#'      Ignored for true_colour plots.
#' @param box_bounds Minimum and maximum latitude and longitude coordinates to map. Defaults to the
#'        entire extent of the model output (though modified by the value of zoom). 
#'        Format: c(longitude_min, longitude_max, latitude_min, latitude_max).
#' @param p Handle for an existing figure if you want to add a layer instead of creating a new figure.
#'        If p is provided, Land_map is over-ridden and set to FALSE.
#' @param return_poly Instead of only the figure handle, return a list containing the figure handle and the dataframe used by geom_plot(). Default FALSE.
#' @param label_towns Add labels for town locations to the figure. Default TRUE
#' @param strict_bounds Obsolescent: ignored
#' @param mark_points Data frame containing longitude and latitude of geolocations to mark with crosses (or a vector containing one location). Default NULL.
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
                       Land_map = FALSE,
                       input_file = "catalog",
                       input_grid = NA,
                       scale_col = c('ivory', 'coral4'), 
                       scale_lim =c(NA, NA),
                       zoom = 6, 
                       box_bounds = c(NA, NA, NA, NA), 
                       p = NA, 
                       suppress_print = FALSE, 
                       return_poly = FALSE,
                       label_towns = TRUE,
                       strict_bounds = FALSE,
                       mark_points = NULL,
                       gbr_poly = FALSE)
{

input_file <- substitute_filename(input_file)

# Check whether this is a locally-stored netcdf file or a web-served file
#if (substr(input_file, 1, 4)=="http") {
#  local_file = FALSE
#} else local_file = TRUE
# @Diego:
# I'm now setting local_file to TRUE regardless. Using this used to save a lot of memory and speed things up, but doesn't seem to 
# do so in the current ncdf4 version and causes problems with the latest ncdf4 library because it can't access variables 
# not included (such as latitude, x_centre).
# I'm leaving the code that uses local_file==FALSE in just in case for now because it was a pain to work out, but it
# should probably be deleted -- can always recover from an older version if needed.
local_file <- TRUE

if (length(p)!=1) Land_map <- FALSE
if (suppress_print) Land_map <- FALSE

towns <- data.frame(latitude = c(-15.47027987, -16.0899, -16.4840, -16.92303816, -19.26639219, -20.0136699, -20.07670986, -20.40109791, -21.15345122, -22.82406858, -23.38031858, -23.84761069, -24.8662122, -25.54073075, -26.18916037),
                    longitude = c(145.2498605, 145.4622, 145.4623, 145.7710, 146.805701, 148.2475387, 146.2635394, 148.5802016, 149.1655418, 147.6363616, 150.5059485, 151.256349, 152.3478987, 152.7049316, 152.6581893),
                    town = c('Cooktown', 'Cape Tribulation', 'Port Douglas', 'Cairns', 'Townsville', 'Bowen', 'Charters Towers', 'Prosperine', 'Mackay', 'Clermont', 'Rockhampton', 'Gladstone', 'Bundaberg', 'Maryborough', 'Gympie'))

# Check whether this is a GBR1 or GBR4 ereefs file, or something else
ereefs_case <- get_ereefs_case(input_file)
input_stem <- get_file_stem(input_file)
#check_platform_ok(input_stem)
grids <- get_ereefs_grids(input_file, input_grid)
x_grid <- grids[['x_grid']]
y_grid <- grids[['y_grid']]

# Allow user to specify a depth below MSL by setting layer to a negative value
if (layer<=0) {
   z_grid <- grids[['z_grid']]
   layer <- max(which(z_grid<layer))
}


# Date to map:
if (is.vector(target_date)) {
	target_date <- as.Date(paste(target_date[1], target_date[2], target_date[3], sep='-'))
} else if (is.character(target_date)) {
	target_date <- as.Date(target_date)
}
# @Diego: to explain why this is so convoluted:
# We have to do different things here depending on what type of eReefs output file we have. Earlier ereefs runs did not
# provide ncml files or OpeNDAP catalog information, and this is still true of some runs on servers other than the NCI server.
# In this case, we need to make some assumptions about the data structures. 4km resolution ereefs model output was served in
# monthly blocks, so we need to find the right month and then step through it. 1km resolution ereefs model output was served
# in daily netcdf files. RECOM files have all the data in one file. ncml files, where provided, can be treated as if all the
# data were in one file (by opening the shell ncml file instead of the individual nc files). Where the server provides a
# catalog and the user has chosen this approach, we can use the metadata in the catalog.
if ((!is.na(ereefs_case[1]))&&(ereefs_case[2]!="unknown")) {
  if (ereefs_case[2] == '4km') { 
	  input_file <- paste0(input_stem, format(target_date, '%Y-%m'), '.nc')
	  nc <- safe_nc_open(input_file)
	  if (!is.null(nc$var[['t']])) { 
	      ds <- as.Date(safe_ncvar_get(nc, "t"), origin = as.Date("1990-01-01"))
          } else {
	      ds <- as.Date(safe_ncvar_get(nc, "time"), origin = as.Date("1990-01-01"))
	  }
	  day <- which.min(abs(target_date - ds))
	  ncdf4::nc_close(nc)
  } else if (ereefs_case[2] == '1km') {
	  day <- 1
	  ds <- target_date
	  input_file <- paste0(input_stem, format(target_date, '%Y-%m-%d'), '.nc')
  }
} else { #recom or other netcdf or ncml file
  #input_file <- input_file
	nc <- safe_nc_open(input_file)
	if (!is.null(nc$var[['t']])) { 
	    ds <- as.Date(safe_ncvar_get(nc, "t"), origin = as.Date("1990-01-01"))
        } else {
	    ds <- as.Date(safe_ncvar_get(nc, "time"), origin = as.Date("1990-01-01"))
	}
  dum1 <- abs(target_date - ds)
  if (min(dum1)>1) warning(paste("Target date", target_date, "is", min(dum1), "days from closest available date in", input_file, ds[min(dum1)]))
	day <- which.min(abs(target_date - ds))
	ncdf4::nc_close(nc)
}
#day <- day + 6

# Allow for US English:
if (var_name == "true_color") {
	var_name <- "true_colour"
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


# There are a few options for var_name. In most cases, we just plot the variable with the given name from the netcdf file.
# There are two special cases using optical output data: "plume" calculates and plots the optical colour class of the water,
# while "true_colour" produces something that looks like a satellite image from the model output
if (var_name=="plume") {
    if (!local_file) {
      input_file <- paste0(input_file, '?R_412,R_443,R_488,R_531,R_547,R_667,R_678')
    } else input_file <- input_file
    nc <- safe_nc_open(input_file)
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
      input_file <- paste0(input_file, '?R_470,R_555,R_645')
    } else input_file <- input_file
    nc <- safe_nc_open(input_file)
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
      input_file <- paste0(input_file, '?ZooL_N,ZooS_N')
    } else input_file <- input_file
    nc <- safe_nc_open(input_file)
    # We don't yet know the dimensions of the variable, so let's get them
    dims <- nc$var[['ZooL_N']][['size']]
    if (is.null(dims)) stop(paste('ZooL_N', ' not found in netcdf file.')) 
    ndims <- length(dims)
    if ((ndims > 3) && (layer == 'surface')) layer <- dims[3]
    var_longname <- "Total zooplankton nitrogen"
    var_units <- "mg N m-3"
} else if (var_name == 'speed') {
    if (!local_file ) {
      input_file <- paste0(input_file, '?u,v')
    } else input_file <- input_file
    nc <- safe_nc_open(input_file)
    # We don't yet know the dimensions of the variable, so let's get them 
    dims <- nc$var[['u']][['size']] 
    if (is.null(dims)) stop(paste('u', ' not found in netcdf file.')) 
    ndims <- length(dims) 
    if ((ndims > 3) && (layer == 'surface')) layer <- dims[3]
    var_longname <- "Current speed"
    var_units <- "m s-1"
} else { 
    if (!local_file) {
      input_file <- paste0(input_file, '?', var_name)
    } else input_file <- input_file
    nc <- safe_nc_open(input_file)
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
       ems_var <- safe_ncvar_get(nc, 'u1', start=c(xmin,ymin,layer,day), count=c(xmax-xmin,ymax-ymin,1,1))
       ems_var <- sqrt(ems_var^2 + safe_ncvar_get(nc, 'u2', start=c(xmin,ymin,layer,day), count=c(xmax-xmin,ymax-ymin,1,1))^2)
    } else {
       ems_var <- safe_ncvar_get(nc, 'u1', start=c(xmin,ymin,day), count=c(xmax-xmin,ymax-ymin,1))
       ems_var <- ems_var + safe_ncvar_get(nc, 'u2', start=c(xmin,ymin,day), count=c(xmax-xmin,ymax-ymin,1))
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

nc <- safe_nc_open(input_file)
if (!is.null(nc$var[['botz']])) {
   botz <- safe_ncvar_get(nc, 'botz', start=c(xmin,ymin), count=c(xmax-xmin,ymax-ymin))
} else {
   botz <- array(NA, dim(ems_var))
}
if (is.null(nc$var[['latitude']])) {
# Standard EMS output file
  latitude <- safe_ncvar_get(nc, "y_centre", start=c(xmin,ymin), count=c(xmax-xmin, ymax-ymin))
  longitude <- safe_ncvar_get(nc, "x_centre", start=c(xmin,ymin), count=c(xmax-xmin, ymax-ymin))
} else { 
  # Simple format netcdf file
  latitude <- safe_ncvar_get(nc, "latitude", start=c(xmin,ymin), count=c(xmax-xmin, ymax-ymin))
  longitude <- safe_ncvar_get(nc, "longitude", start=c(xmin,ymin), count=c(xmax-xmin, ymax-ymin))
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
values <- data.frame(id = id, value = n, depth=botz, x_centre=longitude, y_centre=latitude)
positions <- data.frame(id=rep(id, each=4), x = gx, y = gy)
datapoly <- merge(values, positions, by = c("id"))

if ((var_name!="true_colour")&&(is.na(scale_lim[1]))) { 
	scale_lim <- c(min(n, na.rm=TRUE), max(n, na.rm=TRUE))
}

if (Land_map) {
  MapLocation<-c(min(gx, na.rm=TRUE)-0.5, 
 		min(gy, na.rm=TRUE)-0.5, 
 		max(gx, na.rm=TRUE)+0.5, 
 		max(gy, na.rm=TRUE)+0.5)
  p <- ggplot2::ggplot() +
           ggplot2::geom_polygon(data = map.df, colour = "black", fill="lightgrey", size=0.5, aes(x = long, y=lat, group=group))
} else if (length(p)==1) {
  p <- ggplot2::ggplot()
}

p <-  p + ggplot2::geom_polygon(ggplot2::aes(x=x, y=y, fill=value, group=id), data = datapoly) 

if (var_name=="true_colour") {
  p <- p + ggplot2::scale_fill_identity() 
} else if (scale_col[1] == "spectral") { 
  p <- p + ggplot2::scale_fill_distiller(palette = 'Spectral',
                                         na.value="transparent", 
                                         guide="colourbar", 
                                         limits=scale_lim, 
                                         name=var_units, 
                                         oob=scales::squish) 
} else if (length(scale_col)<3) { 
  if (length(scale_col)==1) scale_col <- c('ivory', scale_col) 
  p <- p + ggplot2::scale_fill_gradient(low=scale_col[1], 
                                        high=scale_col[2], 
				                                na.value="transparent", 
				                                guide="colourbar",
				                                limits=scale_lim,
				                                #name=expression(paste('cells m'^'-3')),
				                                name=var_units,
				                                oob=scales::squish) 
} else { 
  p <- p + ggplot2::scale_fill_gradient2(low=scale_col[1], 
                                         mid=scale_col[2], 
                                         high=scale_col[3], 
                                         na.value="transparent", 
                                         guide="colourbar",
     			                               limits=scale_lim,
                                         midpoint=(scale_lim[2] - scale_lim[1])/2,
                                         space="Lab",
      		                               #name=expression(paste('cells m'^'-3')),
  	    	                               name=var_units,
  		                                   oob=scales::squish) 
}

if (label_towns) {
   towns <- towns[towns$latitude>=min(gy, na.rm=TRUE),]
   towns <- towns[towns$latitude<=max(gy, na.rm=TRUE),]
   towns <- towns[towns$longitude>=min(gx, na.rm=TRUE),]
   towns <- towns[towns$longitude<=max(gx, na.rm=TRUE),]
   if (dim(towns)[1]>0) p <- p + ggplot2::geom_text(data=towns, ggplot2::aes(x=longitude, y=latitude, label=town, hjust="right"), nudge_x=-0.1) +
                                 ggplot2::geom_point(data=towns, ggplot2::aes(x=longitude, y=latitude))
}

p <- p + ggplot2::ggtitle(paste(var_longname, format(chron::chron(as.numeric(ds[day])+0.000001), "%Y-%m-%d %H:%M"))) +
    ggplot2::xlab('longitude') + ggplot2::ylab('latitude')
if (!is.null(mark_points)) {
  if (is.null(dim(mark_points))) mark_points <- data.frame(latitude = mark_points[1], longitude = mark_points[2])
  p <- p + ggplot2::geom_point(data=mark_points, ggplot2::aes(x=longitude, y=latitude), shape=4)
}
if (gbr_poly) {
  p <- p + ggplot2::geom_path(data=sdf.gbr, ggplot2::aes(y=lat, x=long, group=group))
}
if (all(is.na(box_bounds))) { 
  p <- p + ggplot2::coord_map(xlim = c(min(gx, na.rm=TRUE), max(gx, na.rm=TRUE)), ylim = c(min(gy, na.rm=TRUE), max(gy, na.rm=TRUE)))
} else { 
  p <- p + ggplot2::coord_map(xlim = box_bounds[1:2], ylim = box_bounds[3:4]) + 
    ggplot2::theme(panel.border = ggplot2::element_rect(linetype = "solid", colour="grey", fill=NA))
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
#' ALSO CALCULATES THE TEMPORAL MEAN value of each cell over the specified time (visualisation of maps can
#' be suppressed by setting suppress_print to TRUE if this is the primary desired output).
#' Maps produced are optionally overlain on a map of Queensland.
#' Can be more efficient than calling map_ereefs multiple times if you 
#' want to produce an animation because it loads a month at a time for GBR4 runs (unless selected via catalog). 
#' If output files contain multiple outputs per day, chooses the step closest to midday and uses only daily output.
#' To stitch together the images into an animation, you will need other software such as ImageMagick (recommended)
#' or imageJ.  Barbara Robson (AIMS).
#'
#' TODO: Update to use chron dates as has been done for functions in data_extraction_functions.R
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
#'	c(year, month, day, hour), c(year, month, day) or (year, month); b) a date obtained e.g. 
#'  from as.Date(); or c) a character string formatted for input to as.Date(). Defaults to c(2015,2,1).
#' @param end date for animation. Can either be a) a vector of integers in the format 
#'	c(year, month, day); b) a date obtained e.g. from as.Date(); or c) a character string 
#'      formatted for input to as.Date(). Defaults to c(2016,3,31).
#' @param layer Either an integer layer number or 'surface' to choose the surface layer. Defaults to 'surface'.
#' @param output_dir Path to directory in which to store images for animation. Created if necessary. Defaults
#'      to 'ToAnimate'. Images are created in this directory with input_files beginning with var_name, 
#'      followed by an underscore and then sequential numbers beginning with 100001.
#' @param Land_map Set to TRUE to show a land map of Queensland. Default now FALSE.
#' @param input_file is the URL or file location of any of the EMS output files or a THREDDS catalog URI. 
#'        Defaults to a menu selection based on current NCI catalogs. Can also be set to "nci", "menu" or "catalog" for the same behaviour.
#'        Set to "old_menu" to provide old menu options instead of menu options from the NCI catalog.
#'        Numeric values are interpreted as references to selections available from the old menu.
#'        Short codes can be used for some options (codenames as used in https://research.csiro.au/ereefs/models/model-outputs/access-to-raw-model-output/ )
#' @param input_grid Name of the locally-stored or opendap-served netcdf file that contains the grid
#'      coordinates for the cell corners (x_grid and y_grid). If not specified, the function will first look for
#'      x_grid and y_grid can be found in the first INPUT_STEM file, and if not found, will check whether the size 
#'      of the variables in the input file corresponds to the size expected for GBR4 or GBR1, and load appropriate 
#'      x and y grids from data files stored in this package. Alternatively, you can provide the location of a full 
#'      (not simple-format) ereefs netcdf output file such as 
#'      "http://dapds00.nci.org.au/thredds/dodsC/fx3/gbr4_hydro_all/gbr4_all_2016-09.nc"
#' @param scale_col Vector of colours to use for the colour scale. This can be colours 
#'      from the ggplot colour palette or a RGB hash code. Ignored for true_colour plots. 
#'      If set to "spectral", uses a colour spectrum from bluish to red (similar to jet but less vivid). Otherwise:
#'      If one value is given, low colour is set to ivory and high colour to the value given.
#'      If two values are given, these are used as low and high limit colours.
#'      If three values are given, the middle value is used to set the mid-point of the scale.
#'      Defaults to c('ivory', 'coral4').
#' @param box_bounds Minimum and maximum latitude and longitude coordinates to map. Defaults to the
#'        entire extent of the model output (though modified by the value of zoom). 
#'        Format: c(longitude_min, longitude_max, latitude_min, latitude_max). It is recommended to
#'        also specify an appropriate value for zoom if specifying box_bounds.
#' @param suppress_print Set to TRUE if you don't want the plots generatedand saved. Defaults to FALSE.
#' @param stride Default 'daily', but can otherwise be set to a numeric interval indicating how many time-steps to step forward for each frame.
#' @param verbosity Set 0 for just a waitbar, 1 for more updates, 2 for debugging information. Default 0.
#' @param label_towns Add labels for town locations to the figure. Default TRUE
#' @param strict_bounds Obsolescent: ignored
#' @param mark_points Data frame containing longitude and latitude of geolocations to mark with crosses (or a vector containing one location). Default NULL.
#' @param gbr_poly TRUE to show contours of approximate reef areas. Default FALSE.
#' @param add_arrows TRUE to show arrows indicating magnitude and direction of flow. Default FALSE.
#' @param max_u Velocity at which to show maximum arrow length. Default NA, in which case it will use the maximum observed velocity.
#' @param scale_arrows Factor by which to scale arrows. Values >1 result in longer arrows. Values <1 result in shorter arrows. Default 1.
#' @param show_bathy TRUE to show contours based on the bathymetry as represented in the model. Default FALSE. Requires model file to contain botz (this 
#'        requirement may be dropped in future versions for GBR1 and GBR4 runs).
#' @param contour_breaks Depths in metres to show with show_bathy. Default c(5, 10, 20).
#' @return a list that includes data.frame formatted for use in ggplot2::geom_polygon, containing a map of the temporally averaged
#'         value of the variable specified in VAR_NAME over the selected interval, plus the actual values and cell centre latitudes and longitudes.
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
                             Land_map = FALSE, 
                             input_file = "catalog",
                             input_grid = NA, 
                             scale_col = c('ivory', 'coral4'), 
                             scale_lim = c(NA, NA), 
                             zoom = 6, 
                             box_bounds = c(NA, NA, NA, NA), 
                             suppress_print=FALSE, 
                             stride = 'daily',
                             verbosity=0, 
                             label_towns = TRUE,
                             strict_bounds = FALSE,
                             mark_points = NULL,
                             gbr_poly = FALSE,
                             add_arrows = FALSE,
                             max_u = NA,
                             scale_arrows = NA,
                             show_bathy=FALSE,
                             contour_breaks=c(5,10,20))
{
  plot_eta <- FALSE
  input_file <- substitute_filename(input_file)
  if (verbosity > 1) print(paste('After substitute_filename() input_file = ', input_file))

  # Check whether this is a locally-stored netcdf file or a web-served file
  if (substr(input_file, 1, 4)=="http") {
    local_file = FALSE
  } else local_file = TRUE
  #temporary hack:
  local_file <- TRUE

  # Check whether this is a GBR1 or GBR4 ereefs file, a THREDDS catalog or something else
  ereefs_case <- get_ereefs_case(input_file) 
  if (verbosity > 1) print(paste('ereefs_case = ', ereefs_case[1], ereefs_case[2]))
  if (ereefs_case[2]=='1km') warning('Assuming that only one timestep is output per day/file') # find matching commented warning to fix this
  input_stem <- get_file_stem(input_file)
  if (verbosity > 1) print(paste('input_stem = ', input_stem))
  #check_platform_ok(input_stem)

  towns <- data.frame(latitude = c(-15.47027987, -16.0899, -16.4840, -16.92303816, -19.26639219, -20.0136699, -20.07670986, -20.40109791, -21.15345122, -22.82406858, -23.38031858, -23.84761069, -24.8662122, -25.54073075, -26.18916037),
                    longitude = c(145.2498605, 145.4622, 145.4623, 145.7710, 146.805701, 148.2475387, 146.2635394, 148.5802016, 149.1655418, 147.6363616, 150.5059485, 151.256349, 152.3478987, 152.7049316, 152.6581893),
                    town = c('Cooktown', 'Cape Tribulation', 'Port Douglas', 'Cairns', 'Townsville', 'Bowen', 'Charters Towers', 'Prosperine', 'Mackay', 'Clermont', 'Rockhampton', 'Gladstone', 'Bundaberg', 'Maryborough', 'Gympie'))

  # Dates to map. We offer quite a few date format options for flexibility.
  if (is.vector(start_date)) {
    if ((length(start_date==2)) && is.character(start_date[1])) { 
          start_date <- chron::chron(start_date[1], start_date[2], format=c('d-m-y', 'h:m:s'), 
                                     origin=c(year=1990, month=1, day=1))
    } else if (length(start_date==3)) { 
      # Set time to midday
      if (verbosity > 1) {
        print('Setting time to midday. start_date = ')
        print(start_date)
      }
      start_date <- chron::chron(paste(start_date[3], start_date[2], start_date[1], sep = '-'), 
                                 "12:00:00", format=c('d-m-y', 'h:m:s'), 
                                 origin=c(year=1990, month=1, day=1))
    } else if (length(start_date==4)) {
       if (!is.character(start_date[4])) start_date[4] <- paste0(start_date[4], ':00')
       start_date <- chron::chron(paste(start_date[3], start_date[2], start_date[1], sep = '-'), 
                                  start_date[4], format=c('d-m-y', 'h:m:s'), 
                                  origin=c(year=1990, month=1, day=1)) 
    } else {
      stop("start_date format not recognised")
    }
  } else if (is.character(start_date)) {
    start_date <- chron::chron(start_date, "12:00:00", format=c('d-m-y', 'h:m:s'), 
                                 origin=c(year=1990, month=1, day=1))
  } else if (class(start_date)[1] == "Date") {
    start_date <- chron::as.chron(start_date)
  }
  start_day <- as.integer(chron::days(start_date))
  start_tod <- as.numeric(start_date) - as.integer(start_date)
  start_month <- as.integer(months(start_date))
  start_year <- as.integer(as.character(chron::years(start_date)))

  if (is.vector(end_date)) {
    if (length(end_date==2) && is.character(end_date[1])) { 
          end_date <- chron::chron(end_date[1], end_date[2], format=c('d-m-y', 'h:m:s'), 
                                     origin=c(year=1990, month=1, day=1))
    } else if (length(end_date==3)) { 
      # Set time to midday
      end_date <- chron::chron(paste(end_date[3], end_date[2], end_date[1], sep = '-'), "12:00:00", format=c('d-m-y', 'h:m:s'), 
                                 origin=c(year=1990, month=1, day=1))
    } else if (length(end_date==4)) {
       if (!is.character(end_date[4])) end_date[4] <- paste0(end_date[4], ':00')
       end_date <- chron::chron(paste(end_date[3], end_date[2], end_date[1], sep = '-'), end_date[4], format=c('d-m-y', 'h:m:s'), 
                                  origin=c(year=1990, month=1, day=1)) 
    } else {
      stop("end_date format not recognised")
    }
  } else if (is.character(end_date)) {
    end_date <- chron::chron(end_date, "12:00:00", format=c('d-m-y', 'h:m:s'), 
                                 origin=c(year=1990, month=1, day=1))
  } else if (class(end_date)[1] == "Date") {
    end_date <- chron::as.chron(end_date)
  }
  end_day <- as.integer(chron::days(end_date))
  end_month <- as.integer(months(end_date))
  end_year <- as.integer(as.character(chron::years(end_date)))
  
  if (start_date > end_date) {
    stop('start_date must preceed end_date')
  }

  # Points of interest provided by the user to mark on the map, and possibly plot a surface elevation timeseries for:
  if (!is.null(mark_points)) {
   # If mark_points is a vector, change it into a data frame
   if (is.null(dim(mark_points))) {
     mark_points <- data.frame(latitude = mark_points[1], longitude = mark_points[2])
   }
   eta_data <- get_ereefs_ts(var_name='eta', input_file = input_file, start_date=start_date, end_date = end_date, location_latlon = mark_points)
   names(eta_data) <- c('date', 'eta')
   eta_plot <- ggplot2::ggplot(eta_data, ggplot2::aes(x=date, y=eta)) + ggplot2::geom_line() + ggplot2::ylab('surface elevation (m)')
   plot_eta <- TRUE
  }

  if (suppress_print) Land_map <- FALSE
  
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
       if (ereefs_case[2] == "4km") { 
         input_file <- paste0(input_stem, format(as.Date(paste(year, month, 1, sep="-")), '%Y-%m'), '.nc')
       } else if (ereefs_case[2] == "1km") {
         input_file <- paste0(input_stem, format(as.Date(paste(year, month, from_day, sep="-")), '%Y-%m-%d'), '.nc')
       } else { #recom or other netcdf or ncml file
         #input_file <- input_file
         ds<- get_origin_and_times(input_file)[[2]]
         ds_original <- ds
       }
       grids <- get_ereefs_grids(input_file, input_grid)
       x_grid <- grids[['x_grid']]
       y_grid <- grids[['y_grid']]

       # Allow user to specify a depth below MSL by setting layer to a negative value
       if (layer<=0) {
          z_grid <- grids[['z_grid']]
          layer <- max(which(z_grid<layer))
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

       nc <- safe_nc_open(input_file)
       if (is.null(nc$var[['latitude']])) {
       # Standard EMS output file
         latitude <- safe_ncvar_get(nc, "y_centre")
         longitude <- safe_ncvar_get(nc, "x_centre")
         botz <- safe_ncvar_get(nc, 'botz')
       } else { 
         # Simple format netcdf file
         latitude <- safe_ncvar_get(nc, "latitude")
         longitude <- safe_ncvar_get(nc, "longitude")
         botz <- NULL
         if (show_bathy) warning('Can not show bathymetry: simple format netcdf file does not contain botz')
       }
       ncdf4::nc_close(nc)

       if (add_arrows) {
         idim <- dim(latitude)[1]
         jdim <- dim(latitude)[2]
         max_arrow <- max(max(abs(longitude[idim, jdim] - longitude[1,1])/idim), max(abs(latitude[idim, jdim] - latitude[1,1])/jdim))
       }

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
       longitude <- c(longitude)[gx_ok&gy_ok]
       latitude <- c(latitude)[gx_ok&gy_ok]

    } else {
       from_day <- 1
       start_tod <- 0
    }

    if ((start_year==end_year)&&(start_month==end_month)) {
       day_count <- end_day - start_day + 1
    } else if (mcount == 1) {
       day_count <- daysIn(as.Date(paste(year, month, 1, sep='-'))) - start_day + 1
    } else if (mcount == (length(mths))) {
       day_count <- end_day
    } else {
       day_count <- daysIn(as.Date(paste(year, month, 1, sep='-')))
    }

    if (ereefs_case[2] == '4km') { 
      fileslist <- 1
	    input_file <- paste0(input_stem, format(as.Date(paste(year, month, 1, sep="-")), '%Y-%m'), '.nc')
      ds <- get_origin_and_times(input_file)[[2]]
      if ((ds[length(as.numeric(ds))] - as.numeric(ds)[1]) > 31.5) warning('Filename looks like a monthly output file (i.e. contains two dashes) but file contains more than a month of data.')
      if ((ds[length(as.numeric(ds))] - as.numeric(ds)[1]) < 27) warning('Filename looks like a monthly output file (i.e. contains two dashes) but file contains less than a month of data.')
      if(ds[2]==ds[1]) stop(paste('Error reading time from', input_file, '(t[2]==t[1])'))
      tstep <- as.numeric(ds[2]-ds[1])
      day_count <- day_count / tstep
      if (day_count > length(as.numeric(ds))) {
        warning(paste('end_date', end_date, 'is beyond available data. Ending at', max(ds)))
        day_count <- length(as.numeric(ds))
      }
      #dum1 <- as.integer((from_day - 0.4999999)/tstep + 1)
      #dum2 <- as.integer((day_count - 1) / tstep) +1
      #ds <- ds[seq(from=dum1, by=as.integer(1/tstep), to=(dum1+dum2))]
      #start_array <- c(xmin, ymin, dum1)
	    #count_array <- c(xmax-xmin, ymax-ymin, dum2)
      ds <- ds[from_day:day_count]
      start_array <- c(xmin, ymin, from_day)
	    count_array <- c(xmax-xmin, ymax-ymin, day_count)
	    fileslist <- 1
    } else if (ereefs_case[2] == '1km') { 
	    fileslist <- from_day:(from_day+day_count-1)
	    from_day <- 1
	    day_count <- 1

	    input_file <- paste0(input_stem, format(as.Date(paste(year, month, from_day, sep="-")), '%Y-%m-%d'), '.nc')
      ds <- get_origin_and_times(input_file, as_chron="FALSE")[[2]]
      dum1 <- which.min(abs(as.numeric(ds - (from_day + 0.4999999))))
      start_array <- c(xmin, ymin, dum1) 
      count_array <- c(xmax-xmin, ymax-ymin, dum1)
      tstep <- 1
    } else if ((ereefs_case[2] == "recom")|(ereefs_case[1] == "ncml")) {
	    # Everything is in one file but we are only going to read a month at a time
	    # Output may be more than daily, or possibly less
	    # input_file <- paste0(input_stem, '.nc') # input_file has been set previously 
      tstep <- as.numeric(median((ds[2:length(ds)] - ds[1:(length(ds) - 1)]), na.rm=TRUE))
      day_count <- day_count / tstep
      if (day_count > length(as.numeric(ds_original))) {
        warning(paste('end_date', end_date, 'is beyond available data. Ending at', max(ds_original)))
        day_count <- length(as.numeric(ds_original))
      }
      from_day <- (as.numeric(chron::chron(paste(year, month, from_day, sep = '-'), format=c('y-m-d'),
                                  origin=c(year=1990, month=1, day=1)) - as.numeric(ds_original)[1]) + 
                              start_tod) / tstep + 1 
      #print(paste(year, month, from_day, ds[1], start_tod, tstep))
	    if (from_day<1) from_day <-1
	    start_array <- c(xmin, ymin, floor(from_day))
	    count_array <- c(xmax-xmin, ymax-ymin, if(tstep>1) {as.integer(day_count/tstep)} else  day_count)
      #print(start_array)
      #print(count_array)
	    fileslist <- 1
    } else stop("Shouldn't happen: ereefs_case not recognised")

    if (stride == 'daily') {
      if (tstep <= 1.0) {
       stride <- 1/tstep
      } else {
        stop("Minimum timestep in netcdf file is greater than daily.")
      }
    } 
    stride <- as.integer(stride)

    for (dcount in 1:length(fileslist)) {
      if (ereefs_case[2] == '1km') { 
         input_file <- paste0(input_stem, format(as.Date(paste(year, month, fileslist[dcount], sep="-")), '%Y-%m-%d'), '.nc') 
         #ds <- as.Date(paste(year, month, fileslist[dcount], sep="-", '%Y-%m-%d'))
      }
      if (verbosity>0) print(input_file)
      # There are a few options for var_name. In most cases, we just plot the variable with the given name from the netcdf file.
      # There are two special cases using optical output data: "plume" calculates and plots the optical colour class of the water,
      # while "true_colour" produces something that looks like a satellite image from the model output
      if (var_name=="plume") {
        if (!local_file) {
          slice <- paste0('[', start_array[3]-1, ':', stride, ':', start_array[3] + count_array[3] - 2, ']', # time
                          '[', start_array[2]-1, ':', start_array[2] + count_array[2] - 1, ']', # y
                          '[', start_array[1]-1, ':', start_array[1] + count_array[1] - 1, ']') # x
          input_file <- paste0(input_file, '?R_412', slice, ',R_443', slice, ',R_488', slice, ',R_531', slice, ',R_547', slice, ',R_667', slice, ',R_678', slice)
          count_array <- c(count_array[1:2]+1, count_array[3])
        } else input_file <- input_file
        nc <- safe_nc_open(input_file)
        R_412 <- safe_ncvar_get(nc, "R_412", start=start_array, count=count_array)
        R_443 <- safe_ncvar_get(nc, "R_443", start=start_array, count=count_array)
        R_488 <- safe_ncvar_get(nc, "R_488", start=start_array, count=count_array)
        R_531 <- safe_ncvar_get(nc, "R_531", start=start_array, count=count_array)
        R_547 <- safe_ncvar_get(nc, "R_547", start=start_array, count=count_array)
        R_667 <- safe_ncvar_get(nc, "R_667", start=start_array, count=count_array)
        R_678 <- safe_ncvar_get(nc, "R_678", start=start_array, count=count_array)
        #if (local_file) {
        #  R_412 <- R_412[(start_array[1] - 1) : (start_array[1] + count_array[1] - 1),
        #                 (start_array[2] - 1) : (start_array[2] + count_array[2] - 1),
        #                 seq(from = start_array[3] - 1, to = start_array[3] + count_array[3] - 2, by = stride)]
        #}
        ems_var <- NA*R_678
        if (ereefs_case[2] == '4km') {
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
        count_array <- c(count_array[1:2]+1, count_array[3])
        if (!local_file) {
          input_file <- paste0(input_file, '?R_470', slice, ',R_555', slice, ',R_645', slice)
        } else input_file <- input_file
        nc <- safe_nc_open(input_file)
        TCbright <- 10
        R_470 <- safe_ncvar_get(nc, "R_470", start=start_array, count=count_array) * TCbright
        R_555 <- safe_ncvar_get(nc, "R_555", start=start_array, count=count_array) * TCbright
        R_645 <- safe_ncvar_get(nc, "R_645", start=start_array, count=count_array) * TCbright
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
          nc <- safe_nc_open(input_file)
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
          # If there's only one layer (e.g. a surf.nc file) then we want to reduce ndims accordingly, but not if there is only one time-step
          if ((length(dims[dims!=1])!=ndims)&&(dims[length(dims)]!=1)) ndims <- ndims - 1
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
          start_array <- c(start_array[1:2], layer-1, start_array[3])
          count_array <- c(count_array[1:2]+1, 1, count_array[3])
        } else {
          slice <- paste0('[', start_array[3]-1, ':', stride, ':', start_array[3] + count_array[3] - 2, ']', # time
                          '[', start_array[2]-1, ':', start_array[2] + count_array[2] - 1, ']', # y
                          '[', start_array[1]-1, ':', start_array[1] + count_array[1] - 1, ']') # x
          count_array <- c(count_array[1:2]+1, count_array[3])
        }
        if (verbosity>1) print(paste('slice = ', slice))
        if (var_name == "speed") { 
           if (!local_file) {
             input_file <- paste0(input_file, '?u', slice, ',v', slice)
           } else input_file <- input_file
           nc <- safe_nc_open(input_file)
           ems_var <- sqrt(safe_ncvar_get(nc, 'u', start=start_array, count=count_array)^2 + safe_ncvar_get(nc, 'v', start=start_array, count=count_array)^2)
           vat <- ncdf4::ncatt_get(nc, 'u')
           var_longname <- 'Current speed'
           var_units <- vat$units
        } else if (var_name == "ZooT") {
           if (!local_file) {
             input_file <- paste0(input_file, '?ZooL_N', slice, ',ZooS_N', slice)
           } else input_file <- input_file
           nc <- safe_nc_open(input_file)
           ems_var <- safe_ncvar_get(nc, 'ZooL_N', start=start_array, count=count_array) + safe_ncvar_get(nc, 'ZooS_N', start=start_array, count=count_array)
           vat <- ncdf4::ncatt_get(nc, 'ZooL_N')
           var_longname <- 'Total Zooplankton'
           var_units <- vat$units
        } else {
           if (!local_file) {
             if (add_arrows) {
               input_file <- paste0(input_file, '?',var_name, 'u', 'v', slice)
             } else {
               input_file <- paste0(input_file, '?',var_name, slice)
             }
           } else input_file <- input_file
           if (verbosity>1) print(paste("Before nc_open, input_file = ", input_file))
           nc <- safe_nc_open(input_file)
           browser()
           ems_var <- safe_ncvar_get(nc, var_name, start=start_array, count = count_array)
           if (add_arrows) {
             current_u <- ncvar_get(nc, 'u1', start=start_array, count = count_array)
             current_v <- ncvar_get(nc, 'u2', start=start_array, count = count_array)
             if ((dcount==1)&(is.na(max_u))) {
               max_u <- max(max(abs(c(current_u)), na.rm = TRUE), max(abs(c(current_v)), na.rm=TRUE)) 
               if (!is.na(scale_arrows)) max_u <- max_u / scale_arrows
             }
           }
           vat <- ncdf4::ncatt_get(nc, var_name)
           var_longname <- vat$long_name
           var_units <- vat$units
        }
      }
        #if (local_file) { 
        #  if (ndims == 3) {
        #    ems_var <- ems_var[start_array[1] : (start_array[1] + count_array[1]),
        #                       start_array[2] : (start_array[2] + count_array[2]),
        #                       seq(from = start_array[3], to = start_array[3] + count_array[3] - 1, by = stride)] 
        #  }  else {
        #    ems_var <- ems_var[start_array[1] : (start_array[1] + count_array[1]),
        #                       start_array[2] : (start_array[2] + count_array[2]),
        #                       layer,
        #                       seq(from = start_array[3], to = start_array[3] + count_array[3] - 1, by = stride)] 
        #  }
        #  ds <- ds_original[seq(from = start_array[3], to = start_array[3] + count_array[3] - 1, by = stride)] 
        #  if (add_arrows) {
        #    current_u <- current_u[start_array[1] : (start_array[1] + count_array[1]),
        #                           start_array[2] : (start_array[2] + count_array[2]),
        #                           seq(from = start_array[3], to = start_array[3] + count_array[3] - 1, by = stride)] 
        #    current_v <- current_v[start_array[1] : (start_array[1] + count_array[1]),
        #                 start_array[2] : (start_array[2] + count_array[2]),
        #                 seq(from = start_array[3], to = start_array[3] + count_array[3] - 1, by = stride)] 
        #  }
        #}
      #ds <- as.Date(safe_ncvar_get(nc, "time"), origin = as.Date("1990-01-01"))
    
      dum1 <- length(dim(ems_var))
      if (dum1==2) {
        ems_var <- array(ems_var, dim=c(dim(ems_var), 1))
        if (add_arrows) {
          current_u <- array(current_u, dim=c(dim(current_u), 1))
          current_v <- array(current_v, dim=c(dim(current_v), 1))
        }
      }
      if (add_arrows) {
        del_u <- current_u / max_u * max_arrow
        del_v <- current_v / max_u * max_arrow
      }
    
      ncdf4::nc_close(nc)
    
      d <- dim(ems_var)[3]
      for (jcount in 1:d) {
        ems_var2d <- ems_var[, , jcount]
        if (add_arrows) {
          del_u2d <- del_u[, , jcount]
          del_v2d <- del_v[, , jcount]
        }
        # Values associated with each polygon at chosen timestep
        n <- c(ems_var2d)[gx_ok&gy_ok]
        if (icount==0) {
           if (var_name=='true_colour') {
              temporal_sum <- col2rgb(n)^2
           } else { 
              temporal_sum <- n
           }
        } else {
           if (var_name=='true_colour') {
              temporal_sum <- col2rgb(n)^2 + temporal_sum
           } else { 
              temporal_sum <- temporal_sum + n
           }
        }
    
        # Unique ID for each polygon
        id <- as.factor(1:length(n))
        values <- data.frame(id = id, value = n)
        positions <- data.frame(id=rep(id, each=4), x = gx, y = gy)
        datapoly <- merge(values, positions, by = c("id"))
    
        if (!suppress_print) {
            if ((var_name!="true_colour")&&(is.na(scale_lim[1]))) { 
	            scale_lim <- c(min(n, na.rm=TRUE), max(n, na.rm=TRUE))
            }
  
            if (Land_map) {
               p <- ggplot2::ggplot() + geom_polygon(data = map.df, colour = "black", fill="lightgrey", size=0.5, aes(x = long, y=lat, group=group))
	          } else {
	             p <- ggplot2::ggplot()
            }
            p <- p + ggplot2::geom_polygon(ggplot2::aes(x=x, y=y, fill=value, group=id), data = datapoly)
            if (var_name=="true_colour") { 
              p <- p + ggplot2::scale_fill_identity()
            } else if (scale_col[1] == 'spectral') { 
              p <- p + ggplot2::scale_fill_distiller(palette = 'Spectral',
                                                     na.value="transparent", 
                                                     guide="colourbar", 
                                                     limits=scale_lim, 
                                                     name=var_units, 
                                                     oob=scales::squish) 
            } else if (length(scale_col)<3) { 
              if (length(scale_col) == 1) scale_col <- c('ivory', scale_col)
              p <- p + ggplot2::scale_fill_gradient(low=scale_col[1],
					                                          high=scale_col[2],
					                                          na.value="transparent", 
					                                          guide="colourbar",
					                                          limits=scale_lim,
					                                          name=var_units,
					                                          oob=scales::squish)
            } else {
                  p <- p + ggplot2::scale_fill_gradient2(low=scale_col[1],
                                                         mid=scale_col[2],
					                                               high=scale_col[3],
					                                               na.value="transparent", 
                                                         midpoint=(scale_lim[2] - scale_lim[1])/2,
                                                         space="Lab",
					                                               guide="colourbar",
					                                               limits=scale_lim,
					                                               name=var_units,
					                                               oob=scales::squish)
            }

            if (add_arrows) {
              arrow_df <- data.frame(latitude = c(latitude), longitude = c(longitude), uend = c(del_u2d) + c(longitude), vend = c(del_v2d) + c(latitude))
              p <- p + ggplot2::geom_segment(arrow_df, mapping = ggplot2::aes(x = longitude, y = latitude, xend = uend, yend = vend), arrow = ggplot2::arrow(length = ggplot2::unit(0.1, "cm")))
            }
            if ((show_bathy)&!is.null(botz)) {
              bathy_df <- data.frame(latitude = c(latitude), longitude = c(longitude), depth = c(-botz))
              reg <- marmap::griddify(bathy_df, as.integer(idim/2), as.integer(jdim/2))
              bathy <- marmap::as.bathy(reg)
              bathy_df <- marmap::as.xyz(bathy)
              names(bathy_df) <- c('latitude', 'longitude', 'depth')
              p <- p + ggplot2::geom_contour(bathy_df, mapping = ggplot2::aes(x = longitude, y = latitude, z = depth), colour='black', breaks=contour_breaks)
            }

            if (label_towns) {
              towns <- towns[towns$latitude>=min(gy, na.rm=TRUE),]
              towns <- towns[towns$latitude<=max(gy, na.rm=TRUE),]
              towns <- towns[towns$longitude>=min(gx, na.rm=TRUE),]
              towns <- towns[towns$longitude<=max(gx, na.rm=TRUE),]
              if (dim(towns)[1]>0) p <- p + ggplot2::geom_text(data=towns, ggplot2::aes(x=longitude, y=latitude, label=town, hjust="right"), nudge_x=-0.1) +
                                 ggplot2::geom_point(data=towns, ggplot2::aes(x=longitude, y=latitude))
            }
            #p <- p + ggplot2::ggtitle(paste(var_longname, format(chron::chron(as.double(ds[jcount])+0.000001), "%Y-%m-%d %H:%M")))
            p <- p + ggplot2::ggtitle(paste(var_longname, format(ds[jcount], format=c('d-m-yyyy', 'h:m'))))
            p <- p + ggplot2::xlab("longitude") + ggplot2::ylab("latitude")
            if (!is.null(mark_points)) {
              p <- p + ggplot2::geom_point(data=mark_points, ggplot2::aes(x=longitude, y=latitude), shape=4)
              if (plot_eta) {
                dind <- which(ds[jcount]==eta_data$date)
                p2 <- eta_plot + ggplot2::geom_point(data = eta_data[dind,], ggplot2::aes(x=date, y=eta), size=2, color='red')
              }
            }
            if (gbr_poly) {
              p <- p + ggplot2::geom_path(data=sdf.gbr, ggplot2::aes(y=lat, x=long, group=group)) 
            }
            if (all(is.na(box_bounds))) { 
              p <- p + ggplot2::coord_map(xlim = c(min(gx, na.rm=TRUE), max(gx, na.rm=TRUE)), ylim = c(min(gy, na.rm=TRUE), max(gy, na.rm=TRUE)))
            } else {
              p <- p + ggplot2::coord_map(xlim = box_bounds[1:2], ylim = box_bounds[3:4]) + 
                ggplot2::theme(panel.border = ggplot2::element_rect(linetype = "solid", colour="grey", fill=NA))
            }

            icount <- icount + 1

            if (plot_eta) {
              p <- cowplot::plot_grid(p, p2, ncol=1, rel_heights=c(4,1), axis='lr', align='v')
            }

            if (!file.exists(output_dir)) {
               dir.create(output_dir)
            }
            fname <- paste0(output_dir, '/', var_name, '_', 100000 + icount, '.png', collapse='')
            if (verbosity>0) print(paste(var_longname, format(ds[jcount], format=c('d-m-yyyy', 'h:m'))))
            ggplot2::ggsave(fname, p, dpi=100)
            #rm('p')
         }  else {
            icount <- icount + 1
            p <- NULL
         }
         setTxtProgressBar(pb,icount/as.integer(end_date-start_date)/tstep*stride)
      } # end jcount loop
  } # end fileslist loop
  } # end month loop
  close(pb)
  if (var_name=='true_colour') {
    temporal_sum <- rgb(sqrt(temporal_sum)/icount, maxColorValue=255)
    values <- data.frame(id = id, value = temporal_sum)
  } else {
    values <- data.frame(id = id, value = temporal_sum/icount)
  }
  datapoly <- merge(values, positions, by = c("id"))
  return(list(p=p, datapoly=datapoly, longitude=longitude, latitude=latitude))
}

#' Create a map using a dataframe in the format required by ggplot2::geom_plot, for instance from map_ereefs() or map_ereefs_movie()
#'
#' Plots a map figure in the same format as would be given by map_ereefs(), but using a pre-generated dataframe, instead of
#' processing data directly from ereefs netcdf files. Doesn't work for true_color maps.
#' 
#' @param datapoly A dataframe in the format required by geom_plot(), as provided by map_ereefs() or map_ereefs_movie().
#' @param var_longname Character vector to use for the figure title.
#' @param var_units Units to include in the figure labelling.
#' @param Land_map Set to TRUE to show a land mapof Queensland.  Default now FALSE.
#' @param scale_col Vector of colours to use for the colour scale. This can be colours 
#'      from the ggplot colour palette or a RGB hash code, or "spectral". Ignored for true_colour plots. 
#'      If set to "spectral", uses a colour spectrum from bluish to red (similar to jet but less vivid). Otherwise:
#'      If one value is given (other than "spectral"), low colour is set to ivory and high colour to the value given.
#'      If two values are given, these are used as low and high limit colours.
#'      If three values are given, the middle value is used to set the mid-point of the scale.
#'      Defaults to c('ivory', 'coral4').
#' @param scale_lim Upper and lower bounds for colour scale. Defaults to full range of data.
#'      Ignored for true_colour plots.
#' @param suppress_print Default FALSE. If true, don't prdocue the map image.
#' @param p Handle for an existing figure if you want to add a layer instead of creating a new figure.
#'        If p is provided, Land_map is over-ridden and set to FALSE.
#' @return p Handle for the figure generated.
#' @export
#' @examples
#' \dontrun{
#' a <- plot_map(p)
#' plot_map(a[[2]])
#'}

plot_map <- function(datapoly,
             var_longname = '',
             var_units = '',
  		       Land_map = FALSE,
             scale_col = c('ivory', 'coral4'),
  		       scale_lim = c(NA, NA),
             box_bounds = c(NA, NA, NA, NA), 
             label_towns = TRUE,
             zoom = 6,
  		       p = NA,
             suppress_print = FALSE,
             gbr_poly = FALSE)
{
  if ("datapoly" %in% names(datapoly)) datapoly <- datapoly$datapoly
  if (class(datapoly$value)=="factor") {
     var_name <- "true_colour"
  } else {
     var_name <- "something else"
  }
  if ((var_name!="true_colour")&&(is.na(scale_lim[1]))) { 
	  scale_lim <- c(min(datapoly$value, na.rm=TRUE), max(datapoly$value, na.rm=TRUE))
  }

  if (suppress_print) Land_map <- FALSE
  if (length(p)!=1) Land_map <- FALSE
  if (Land_map) {
    p <- ggplot2::gplot() + geom_polygon(data = map.df, colour = "black", fill="lightgrey", size=0.5, aes(x = long, y=lat, group=group))
  } else if (length(p)==1) {
    p <- ggplot2::ggplot()
  }
  if (!suppress_print) {
    p <- p + ggplot2::geom_polygon(ggplot2::aes(x=x, y=y, fill=value, group=id), data = datapoly)
    if (var_name=="true_colour") {
      p <- p + ggplot2::scale_fill_identity() 
    } else {
        if (scale_col[1] == 'spectral') { 
          p <- p + ggplot2::scale_fill_distiller(palette = 'Spectral',
                                                 na.value="transparent", 
                                                 guide="colourbar", 
                                                 limits=scale_lim, 
                                                 name=var_units, 
                                                 oob=scales::squish) 
        } else {
          if (length(scale_col)==1) scale_col <- c('ivory', scale_col)
          if (length(scale_col)<3) { 
            p <- p + ggplot2::scale_fill_gradient(low=scale_col[1], 
				                                          high=scale_col[2], 
				                                          na.value="transparent", 
				                                          guide="colourbar",
				                                          limits=scale_lim,
				                                          name=var_units,
				                                          oob=scales::squish) 
          } else { 
            p <- p + ggplot2::scale_fill_gradient2(low=scale_col[1], 
                                                   mid=scale_col[2],
				                                           high=scale_col[3], 
				                                           na.value="transparent", 
				                                           guide="colourbar",
				                                           limits=scale_lim,
                                                   midpoint=(scale_lim[2] - scale_lim[1])/2,
                                                   space="Lab",
				                                           name=var_units,
				                                           oob=scales::squish) 
          }
        }
    }
    p <- p + ggplot2::ggtitle(var_longname) + ggplot2::xlab('degrees East') + ggplot2::ylab('degrees North')
    if (label_towns) {
       towns <- data.frame(latitude = c(-15.47027987, -16.0899, -16.4840, -16.92303816, -19.26639219, -20.0136699, -20.07670986, -20.40109791, -21.15345122, -22.82406858, -23.38031858, -23.84761069, -24.8662122, -25.54073075, -26.18916037),
                    longitude = c(145.2498605, 145.4622, 145.4623, 145.7710, 146.805701, 148.2475387, 146.2635394, 148.5802016, 149.1655418, 147.6363616, 150.5059485, 151.256349, 152.3478987, 152.7049316, 152.6581893),
                    town = c('Cooktown', 'Cape Tribulation', 'Port Douglas', 'Cairns', 'Townsville', 'Bowen', 'Charters Towers', 'Prosperine', 'Mackay', 'Clermont', 'Rockhampton', 'Gladstone', 'Bundaberg', 'Maryborough', 'Gympie'))
       if (dim(towns)[1]>0) {
         p <- p + ggplot2::geom_text(data=towns, ggplot2::aes(x=longitude, y=latitude, label=town, hjust="right"), nudge_x=-0.1) +
                                 ggplot2::geom_point(data=towns, ggplot2::aes(x=longitude, y=latitude))
       }
    }
    if (all(is.na(box_bounds))) { 
      p <- p + ggplot2::coord_map(xlim = c(min(datapoly$x, na.rm=TRUE), max(datapoly$x, na.rm=TRUE)), ylim = c(min(datapoly$y, na.rm=TRUE), max(datapoly$y, na.rm=TRUE)))
    } else {
      p <- p + ggplot2::coord_map(xlim = box_bounds[1:2], ylim=box_bounds[3:4]) +
        ggplot2::theme(panel.border = ggplot2::element_rect(linetype = "solid", colour="grey", fill=NA))
    }
    if (gbr_poly) {
      p <- p + ggplot2::geom_path(data=sdf.gbr, ggplot2::aes(y=lat, x=long, group=group)) 
    }

    print(p)
  }
  return(p)
}
