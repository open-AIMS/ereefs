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
#' plume_class(rsr)

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
#'      formatted for input to as.Date(). Defaults to c(2011,1,30).
#' @param layer Either an integer layer number or 'surface' to choose the surface layer. Defaults to 'surface'.
#' @param Google_map_underlay Set to TRUE (the default) to use ggmap to show a Google Map as
#'      an underlay for the model output plot. Requires the ggmap librray.
#' @param input_stem Path to the locally-stored or opendap-served netcdf file, including the stem
#'      of the file name (i.e. the start of the filename without the date or file extension).
#'      Defaults to "http://dapds00.nci.org.au/thredds/dodsC/fx3/gbr4_bgc_GBR4_H2p0_B2p0_Chyd_Dcrt/gbr4_bgc_simple_".,
#' @param input_grid Name of the locally-stored or opendap-served netcdf file that contains the grid
#'      coordinates for the cell corners (x_grid and y_grid). This is usually the eReefs hydrodynamic model
#'      output file corresponding to the biogeochemical model output file being used. 
#'      Defaults to "http://dapds00.nci.org.au/thredds/dodsC/fx3/gbr4_hydro_all/gbr4_all_2016-09.nc"
#' @param scale_col Vector of colours to use for low and high values in the colour scale. This can be a colour 
#'      from the ggplot colour palette or a RGB hash code. Ignored for true_colour plots. 
#'      Defaults to c('ivory', 'coral4').
#' @param zoom Value of zoom passed to ggmap(). Set to 5 if you want to show the entire extent 
#'      of eReefs models. Defaults to 6. Higher values will zoom in further.
#' @export
#' @examples
#' map_ereefs()
#' map_ereefs('TN')
#' map_ereefs('plume', target_date=c(2011, 1, 30))
#' map_ereefs('Chl_a_sum', target_date='2016-07-15',, scale_col=c('ivory', 'green4'))
map_ereefs <- function(var_name = "true_colour", 
		       target_date = c(2011,1,30), 
                       layer = 'surface',
		       Google_map_underlay = TRUE,
                       input_stem = "http://dapds00.nci.org.au/thredds/dodsC/fx3/gbr4_bgc_GBR4_H2p0_B2p0_Chyd_Dcrt/gbr4_bgc_simple_",
                       input_grid = "http://dapds00.nci.org.au/thredds/dodsC/fx3/gbr4_hydro_all/gbr4_all_2016-09.nc",
                       scale_col = c('ivory', 'coral4'),
                       zoom = 6)
{

# Date to map:
if (is.vector(target_date)) {
	target_date <- as.Date(paste(target_date[1], target_date[2], target_date[3], sep='-'))
} else if (is.character(target_date)) {
	target_date <- as.Date(target_date)
}
day <- as.integer(format(target_date, '%d'))

# Allow for US English:
if (var_name == "true_color") {
	var_name <- "true_colour"
}

# Check for ggmap()
if ((Google_map_underlay)&&(!requireNamespace("ggmap", quietly = TRUE))) {
  print('Package ggmap is required to show a Google map underlay. Preparing plot with no underlay.')
  print('To avoid this message, either install ggmap or set Google_map_underlay = FALSE.')
  Google_map_underlay = FALSE
}

# We take x_grid and y_grid from a hydrodynamic model netcdf file as they aren't present in the latest bgc output
#input_grid <- 'http://cmar-02-mel.it.csiro.au:8080/opendap/scratch/rob575/gbr_924/tricho_924_2011-06.nc'
nc <- ncdf4::nc_open(input_grid)
x_grid <- ncdf4::ncvar_get(nc, "x_grid")
y_grid <- ncdf4::ncvar_get(nc, "y_grid")
ncdf4::nc_close(nc)

ndims <- 0

if (var_name=="plume") {
    inputfile <- paste0(input_stem, format(target_date, '%Y-%m'), 
		            '.nc?R_412,R_443,R_488,R_531,R_547,R_667,R_678,time')
    nc <- ncdf4::nc_open(inputfile)
    R_412 <- ncdf4::ncvar_get(nc, "R_412", start=c(1,1,day), count=c(-1,-1,1))
    R_443 <- ncdf4::ncvar_get(nc, "R_443", start=c(1,1,day), count=c(-1,-1,1))
    R_488 <- ncdf4::ncvar_get(nc, "R_488", start=c(1,1,day), count=c(-1,-1,1))
    R_531 <- ncdf4::ncvar_get(nc, "R_531", start=c(1,1,day), count=c(-1,-1,1))
    R_547 <- ncdf4::ncvar_get(nc, "R_547", start=c(1,1,day), count=c(-1,-1,1))
    R_667 <- ncdf4::ncvar_get(nc, "R_667", start=c(1,1,day), count=c(-1,-1,1))
    R_678 <- ncdf4::ncvar_get(nc, "R_678", start=c(1,1,day), count=c(-1,-1,1))
    rsr <- list(R_412, R_443, R_488, R_531, R_547, R_667, R_678)
    ems_var <- plume_class(rsr)

} else if (var_name=="true_colour") {
    inputfile <- paste0(input_stem, format(target_date, '%Y-%m'), 
		            '.nc?R_470,R_555,R_645,time,eta')
    nc <- ncdf4::nc_open(inputfile)
    TCbright <- 10
    R_470 <- ncdf4::ncvar_get(nc, "R_470", start=c(1,1,day), count=c(-1,-1,1)) * TCbright
    R_555 <- ncdf4::ncvar_get(nc, "R_555", start=c(1,1,day), count=c(-1,-1,1)) * TCbright
    R_645 <- ncdf4::ncvar_get(nc, "R_645", start=c(1,1,day), count=c(-1,-1,1)) * TCbright
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
} else { 
    inputfile <- paste0(input_stem, format(target_date, '%Y-%m'), 
		            '.nc?time,', var_name)
    nc <- ncdf4::nc_open(inputfile)
    if (ndims==0) {
      # We don't yet know the dimensions of the variable, so let's get them
      dims <- nc$var[[var_name]][['size']]
      ndims <- length(dims)
      if (layer == 'surface') layer <- dims[3]
    }
    if (ndims == 4) {
       ems_var <- ncdf4::ncvar_get(nc, var_name, start=c(1,1,layer,day), count=c(-1,-1,1,1))
    } else {
       ems_var <- ncdf4::ncvar_get(nc, var_name, start=c(1,1,day), count=c(-1,-1,1))
    }
}
#ds <- as.Date(ncdf4::ncvar_get(nc, "time"), origin = as.Date("1990-01-01"))

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
if (var_name=="true_colour") { 
    n <- c(ems_var)[gx_ok&gy_ok]
    # (gx_ok and gy_ok should be identical, but let's be certain)
    gx <- c(t(gx[gx_ok&gy_ok,]))
    gy <- c(t(gy[gx_ok&gy_ok,]))
} else {
    ems_ok <- (ems_var!="#000000")
    n <- c(ems_var)[gx_ok&gy_ok&ems_ok]
    # (gx_ok and gy_ok should be identical, but let's be certain)
    gx <- c(t(gx[gx_ok&gy_ok&ems_ok,]))
    gy <- c(t(gy[gx_ok&gy_ok&ems_ok,]))
}

# Unique ID for each polygon
id <- 1:length(n)

id <- as.factor(id)
values <- data.frame(id = id, value = n)
positions <- data.frame(id=rep(id, each=4), x = gx, y = gy)
datapoly <- merge(values, positions, by = c("id"))

if (Google_map_underlay) {
  MapLocation<-c(min(x_grid, na.rm=TRUE)-0.5, 
 		min(y_grid, na.rm=TRUE)-0.5, 
 		max(x_grid, na.rm=TRUE)+0.5, 
 		max(y_grid, na.rm=TRUE)+0.5)
  myMap<-ggmap::get_map(location=MapLocation, source="google", maptype="hybrid", zoom=zoom, crop=TRUE)
  if (var_name=="true_colour") {
    p <- ggmap::ggmap(myMap) +
        ggplot2::geom_polygon(ggplot2::aes(x=x, y=y, fill=value, group=id), data = datapoly) +
	ggplot2::scale_fill_identity()
  } else {
    p <- ggmap::ggmap(myMap) +
        ggplot2::geom_polygon(ggplot2::aes(x=x, y=y, fill=value, group=id), data = datapoly) +
        ggplot2::scale_fill_gradient(low=scale_col[1], high=scale_col[2], na.value="transparent", guide="colourbar")
  }
	#, limits=c(0,.01))
} else {
  p <- ggplot() +
       ggplot2::geom_polygon(ggplot2::aes(x=x, y=y, fill=value, group=id), data = datapoly) +
        ggplot2::scale_fill_gradient(low=scale_col[1], high=scale_col[2], na.value="transparent", guide="colourbar")
#			    limits=c(30,36))
}
print(p)
}

#' Create a series of map image files for an animation of eReefs model output.
#'
#' Creates and saves to disk a sequential series of colour map images showing concentrations of a specified 
#' eReefs model output variable at a specified model layer (by default, the surface layer). The map is optionally 
#' (and by default) overlain on a Google Earth map of the region. The function should work under Linux or MacOS. 
#' The R netcdf4 library doesn't handle opendap served files under Windows, unfortunately, but can be used under 
#' Windows for locally-stored model output files. More efficient than calling map_ereefs multiple times if you 
#' want to produce an animation because it loads a month at a time. By Barbara Robson (AIMS).
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
#'   Barrier Reef, Australia. Marine Pollution Bulletin, 65(4-9), pp.224-235.
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
#'      formatted for input to as.Date(). Defaults to c(2011,1,31).
#' @param end date for animation. Can either be a) a vector of integers in the format 
#'	c(year, month, day); b) a date obtained e.g. from as.Date(); or c) a character string 
#'      formatted for input to as.Date(). Defaults to c(2011,3,30).
#' @param layer Either an integer layer number or 'surface' to choose the surface layer. Defaults to 'surface'.
#' @param output_dir Path to directory in which to store images for animation. Created if necessary. Defaults
#'      to 'ToAnimate'. Images are created in this directory with filenames beginning with var_name, 
#'      followed by an underscore and then sequential numbers beginning with 100001.
#' @param Google_map_underlay Set to TRUE (the default) to use ggmap::ggmap to show a Google Map as
#'      an underlay for the model output plot. Requires the ggmap::ggmap librray.
#' @param input_stem Path to the locally-stored or opendap-served netcdf file, including the stem
#'      of the file name (i.e. the start of the filename without the date or file extension).
#'      Defaults to "http://dapds00.nci.org.au/thredds/dodsC/fx3/gbr4_bgc_GBR4_H2p0_B2p0_Chyd_Dcrt/gbr4_bgc_simple_".,
#' @param input_grid Name of the locally-stored or opendap-served netcdf file that contains the grid
#'      coordinates for the cell corners (x_grid and y_grid). This is usually the eReefs hydrodynamic model
#'      output file corresponding to the biogeochemical model output file being used. 
#'      Defaults to "http://dapds00.nci.org.au/thredds/dodsC/fx3/gbr4_hydro_all/gbr4_all_2016-09.nc"
#' @param ggplot2::scale_col Vector of colours to use for low and high values in the colour scale. This can be a colour 
#'      from the ggplot2::ggplot colour palette or a RGB hash code. Ignored for true_colour plots. 
#'      Defaults to c('ivory', 'coral4').
#' @param zoom Value of zoom passed to ggmap::ggmap(). Set to 5 if you want to show the entire extent 
#'      of eReefs models. Defaults to 6. Higher values will zoom in further.
#' @export
#' @examples
#' map_ereefs_movie()
map_ereefs_movie <- function(var_name = "true_colour", 
		       start_date = c(2011,1,31), 
		       end_date = c(2011,3,30), 
                       layer = 'surface',
                       output_dir = 'ToAnimate',
		       Google_map_underlay = TRUE,
                       input_stem = "http://dapds00.nci.org.au/thredds/dodsC/fx3/gbr4_bgc_GBR4_H2p0_B2p0_Chyd_Dcrt/gbr4_bgc_simple_",
                       input_grid = "http://dapds00.nci.org.au/thredds/dodsC/fx3/gbr4_hydro_all/gbr4_all_2016-09.nc",
                       scale_col = c('ivory', 'coral4'),
                       zoom = 6)
		       
                      
{

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
  stop('start_date must precede end_date')
}

# Check for ggmap::ggmap()
if ((Google_map_underlay)&&(!requireNamespace("ggmap", quietly = TRUE))) {
  print('Package ggmap::ggmap is required to show a Google map underlay. Preparing plot with no underlay.')
  print('To avoid this message, either install ggmap or set Google_map_underlay = FALSE.')
  Google_map_underlay = FALSE
}

if (start_year==end_year) {
    mths <- start_month:end_month
    years <- rep(start_year, length(mths))
} else if ((start_year + 1) == end_year) {
    mths <- c(start_month:12, 1:end_month)
    years <- c(rep(start_year, 12 - start_month + 1), rep(end_year, end_month))
} else {
    mths <- c(start_month:12, rep(1:12, end_year - start_year - 1), 1:end_month)
    years <- c(rep(start_year, 12 - start_month + 1), 
               rep(start_year + 1 : end_year - 1, each=12),
               rep(end_year, end_month))
}

# Allow for US English:
if (var_name == "true_color") {
	var_name <- "true_colour"
}

nc <- ncdf4::nc_open(input_grid)
x_grid <- ncdf4::ncvar_get(nc, "x_grid")
y_grid <- ncdf4::ncvar_get(nc, "y_grid")
ncdf4::nc_close(nc)

if (Google_map_underlay) {
  MapLocation<-c(min(x_grid, na.rm=TRUE)-0.5, 
 		min(y_grid, na.rm=TRUE)-0.5, 
 		max(x_grid, na.rm=TRUE)+0.5, 
 		max(y_grid, na.rm=TRUE)+0.5)
  myMap<-ggmap::get_map(location=MapLocation, source="google", maptype="hybrid", crop=TRUE, zoom=zoom)
}

ndims <- 0
icount <- 0
mcount <- 0
for (month in mths) {
  mcount <- mcount + 1
  year <- years[mcount]

  if (mcount == 1) {
     from_day <- start_day
  } else {
     from_day <- 1
  }
  if (mcount == (length(mths))) {
     day_count <- end_day
  } else if ((start_year==end_year)&&(start_month==end_month)) {
     day_count <- end_day - start_day
  } else {
     day_count <- -1
  }

  if (var_name=="plume_rms") {
    inputfile <- paste0(input_stem, format(as.Date(paste(year, month, 1, sep='-')), '%Y-%m'), 
		            '.nc?R_412,R_443,R_488,R_531,R_547,R_667,R_678')
    nc <- ncdf4::nc_open(inputfile)
    R_412 <- ncdf4::ncvar_get(nc, "R_412", start=c(1,1,from_day), count=c(-1,-1,day_count))
    R_443 <- ncdf4::ncvar_get(nc, "R_443", start=c(1,1,from_day), count=c(-1,-1,day_count))
    R_488 <- ncdf4::ncvar_get(nc, "R_488", start=c(1,1,from_day), count=c(-1,-1,day_count))
    R_531 <- ncdf4::ncvar_get(nc, "R_531", start=c(1,1,from_day), count=c(-1,-1,day_count))
    R_547 <- ncdf4::ncvar_get(nc, "R_547", start=c(1,1,from_day), count=c(-1,-1,day_count))
    R_667 <- ncdf4::ncvar_get(nc, "R_667", start=c(1,1,from_day), count=c(-1,-1,day_count))
    R_678 <- ncdf4::ncvar_get(nc, "R_678", start=c(1,1,from_day), count=c(-1,-1,day_count))
    ems_var <- NA*R_678
    for (day in 1:dim(R_412)[3]) {
       rsr <- list(R_412[,,day], R_443[,,day], R_488[,,day], R_531[,,day], R_547[,,day], R_667[,,day], R_678[,,day])
       ems_var[,,day] <- plume_class(rsr)
    }
  } else if (var_name=="true_colour") {
    inputfile <- paste0(input_stem, format(as.Date(paste(year, month, 1, sep='-')), '%Y-%m'), 
		            '.nc?R_470,R_555,R_645,time')
    nc <- ncdf4::nc_open(inputfile)
    TCbright <- 10
    R_470 <- ncdf4::ncvar_get(nc, "R_470", start=c(1,1,from_day), count=c(-1,-1,day_count)) * TCbright
    R_555 <- ncdf4::ncvar_get(nc, "R_555", start=c(1,1,from_day), count=c(-1,-1,day_count)) * TCbright
    R_645 <- ncdf4::ncvar_get(nc, "R_645", start=c(1,1,from_day), count=c(-1,-1,day_count)) * TCbright
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

  } else { 
    inputfile <- paste0(input_stem, format(as.Date(paste(year, month, 1, sep='-')), '%Y-%m'), 
		            '.nc?time,', var_name)
    nc <- ncdf4::nc_open(inputfile)
    if (ndims == 0) {
      # We don't yet know the dimensions of the variable, so let's get them
      dims <- nc$var[[var_name]][['size']]
      ndims <- length(dims)
      if (layer == 'surface') layer <- dims[3]
    }

    if (ndims == 4) {
       ems_var <- ncdf4::ncvar_get(nc, var_name, start=c(1,1,layer,from_day), count=c(-1,-1,1,day_count))
    } else {
       ems_var <- ncdf4::ncvar_get(nc, var_name, start=c(1,1,from_day), count=c(-1,-1,day_count))
    }
  }
  #ds <- as.Date(ncdf4::ncvar_get(nc, "time"), origin = as.Date("1990-01-01"))

  ncdf4::nc_close(nc)

  dum1 <- length(dim(ems_var))
  if (dum1==2) {
    ems_var <- array(ems_var, dim=c(dim(ems_var), 1))
  }
  a <- dim(ems_var)[1]
  b <- dim(ems_var)[2]
  d <- dim(ems_var)[3]

  if (icount==0) {
    # Set up the polygon corners. 4 per polygon.
    gx <- c(x_grid[1:a, 1:b], x_grid[2:(a+1), 1:b], x_grid[2:(a+1), 2:(b+1)], x_grid[1:a, 2:(b+1)])
    gy <- c(y_grid[1:a, 1:b], y_grid[2:(a+1), 1:b], y_grid[2:(a+1), 2:(b+1)], y_grid[1:a, 2:(b+1)])
    gx <- array(gx, dim=c(a*b,4))
    gy <- array(gy, dim=c(a*b,4))

    # Find and exclude points where not all corners are defined
    gx_ok <- !apply(is.na(gx),1, any)
    gy_ok <- !apply(is.na(gy),1, any)
    gx <- c(t(gx[gx_ok&gy_ok,]))
    gy <- c(t(gy[gx_ok&gy_ok,]))
  }
  
  for (jcount in 1:d) {
    ems_var2d <- ems_var[,,jcount]
    # Values associated with each polygon at chosen timestep
    #n <- c(ems_var[,,tstep])[gx_ok&gy_ok]
      n <- c(ems_var2d)[gx_ok&gy_ok]

    # Unique ID for each polygon
    id <- 1:length(n)

    id <- as.factor(id)
    values <- data.frame(id = id, value = n)
    positions <- data.frame(id=rep(id, each=4), x = gx, y = gy)
    datapoly <- merge(values, positions, by = c("id"))

    if (Google_map_underlay) {
      if (var_name=="true_colour") {
      p <- ggmap::ggmap(myMap) +
          ggplot2::geom_polygon(ggplot2::aes(x=x, y=y, fill=value, group=id), data = datapoly) +
	  ggplot2::scale_fill_identity()
      } else {
      p <- ggmap::ggmap(myMap) +
          ggplot2::geom_polygon(ggplot2::aes(x=x, y=y, fill=value, group=id), data = datapoly) +
          ggplot2::scale_fill_gradient(low="ivory", high="coral4", na.value="transparent", guide="colourbar")
      }
	  #, limits=c(0,.01))
    } else {
      p <- ggplot2::ggplot() +
         ggplot2::geom_polygon(ggplot2::aes(x=x, y=y, fill=value, group=id), data = datapoly) +
          ggplot2::scale_fill_gradient(low="tan4", high="green2", na.value="transparent", guide="colourbar")
    #			    limits=c(30,36))
    }
    icount <- icount + 1
    if (!file.exists(output_dir)) {
      dir.create(output_dir)
    }
    fname <- paste0(output_dir, '/', var_name, '_', 100000 + icount, '.png', collapse='')
    ggplot2::ggsave(fname, p)
  } 
}
}
