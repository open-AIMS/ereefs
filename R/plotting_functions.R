#" Calculate the number of days in the month of the specified date
#"
#" With thanks to G.Grothendieck on StackOverflow.
#"
#" @param d A date in the format given by as.Date
#" @return An integer
#" @export

daysIn <- function(d) {
 	fom <- as.Date(cut(d, "month"))
 	fomp1 <- as.Date(cut(fom + 32, "month"))
  return(as.integer(fomp1 - fom))
}

#" Calculate and return the plume optical class from eReefs model output..
#"
#" Uses reflectance at multiple wavelengths from the model output to calculate plume colour
#" classes as defined by Devlin et al. (2012). Adapted from the Matlab function plume_detect.m 
#" by Mark Baird (CSIRO). This version by Barbara Robson (AIMS).
#"
#"
#" Devlin, M.J., McKinna, L.W., Álvarez-Romero, J.G., Petus, C., Abott, B., Harkness, P. and 
#"   Brodie, J., 2012. Mapping the pollutants in surface riverine flood plume waters in the Great 
#"   Barrier Reef, Australia. Marine pollution bulletin, 65(4-9), pp.224-235.
#"
#" @param rsr A list of 2D vectors containing the reflectances at various wavelengths
#"            extracted from an EMS netcdf file:
#"	      rsr <- list(R_412, R_443, R_488, R_531, R_547, R_667, R_678)
#" @return an array of the same size as R_412 containing values between 1 and 7 correspodning
#"         to optical plume classes.
#" @export
#" @examples
#" \dontrun{
#" plume_class(rsr)
#" }

#' Classify optical plume type from modelled reflectances
#'
#' Calculates plume optical classes from eReefs/EMS reflectance outputs using
#' the Devlin et al. (2012) class definitions as adapted for this package.
#'
#' @param rsr List of 2D reflectance arrays in wavelength order:
#'   `R_412`, `R_443`, `R_488`, `R_531`, `R_547`, `R_667`, `R_678`.
#'
#' @return An integer matrix of the same horizontal shape as the reflectance
#'   inputs, with plume class values between 1 and 7.
#' @export
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

  rms <- array(0, dim = c(xdim, ydim, 7, 7))

  # Inefficient looped code. Vectorise this if it is too slow.
  for (i in 1:7) {
	# Can probably replace the j loop below with something along the lines of:
	# rms[, ,i, ] <- cl[[i]] - rsr # (not quite right)
	for (j in 1:7) {
	    rms[, , i, j] <- cl[[i]][j] - rsr[[j]]
	}
    }
    rms <- rms^2
    rmse <- rowSums(rms, dim = 3)
    rmse[is.na(rmse)] <- 999999
    plume_class <- apply(rmse, c(1, 2), which.min)
    plume_class[is.na(rsr[[1]])] <- NA
    #return(as.factor(plume_class))
    return(plume_class)
}

#" Create a surface map of eReefs model output.
#"
#" Creates a colour map showing concentrations of a specified eReefs model output variable at a specified
#" model layer (by default, the surface layer). The map is optionally (and by default) overlain on a map
#" of Queensland (if off the Queensland coast).
#" By Barbara Robson (AIMS).
#"
#" References:
#"
#" Baird, M.E., Cherukuru, N., Jones, E., Margvelashvili, N., Mongin, M., Oubelkheir, K., 
#"   Ralph, P.J., Rizwi, F., Robson, B.J., Schroeder, T. and Skerratt, J., 2016. Remote-sensing 
#"   reflectance and true colour produced by a coupled hydrodynamic, optical, sediment, 
#"   biogeochemical model of the Great Barrier Reef, Australia: comparison with satellite data. 
#"   Environmental Modelling & Software, 78, pp.79-96.
#"
#" Baird, M.E., Andrewartha, J., Herzfeld, M., Jones, E.M., Margvelashvili, N., Mongin, M., Rizwi, 
#"   F.,  Skerratt, J., Soja-Wozniak, M., Wild-Allen, K., Schroe der, T., Robson, B.J., da Silva, E. 
#"   and Devlin, M., 2017. River plumes of the Great Barrier Reef: freshwater , sediment and optical 
#"   footprints quantified by the eReefs modelling system  In Syme, G., Hatton MacDonald, D., Fulton, 
#"   B. and Piantadosi, J. (eds) MODSIM2017, 22nd International Congress on Modelling and Simulation. 
#"   Modelling and Simulation Society of Australia and New Zealand, December 2017, pp. 894–900
#"
#" Devlin, M.J., McKinna, L.W., Álvarez-Romero, J.G., Petus, C., Abott, B., Harkness, P. and 
#"   Brodie, J., 2012. Mapping the pollutants in surface riverine flood plume waters in the Great 
#"   Barrier Reef, Australia. Marine pollution bulletin, 65(4-9), pp.224-235.
#"
#"
#" @return a ggplot object
#" @param var_name Short name of the variable to plot. This can be any variable in the
#"                 eReefs netcdf file that you are accessing (refer to eReefs model
#"                 documentation or extract variable names from the netcdf file for a full 
#"                 list). In addition, two special var_name values are supported: "true_colour"
#"                 and "plume". "true_colour" provides a simulated MODIS true colour
#"                 image (refer to Baird et al., 2016 for en explanation)."plume"
#"                 provides a map of calculated plume colour class as per Devlin et al. (2012)
#"                 and Baird et al. (2017). Defaults to true_colour.
#" @param target_date Date to map. Can either be a) a vector of integers in the format 
#"	c(year, month, day); b) a date obtained e.g. from as.Date(); or c) a character string 
#"      formatted for input to as.Date(). Defaults to c(2018, 1, 30). Be careful to choose the right
#"      input file, as the function will plot the closest date in the file to the target date without
#"      complaint, however far fron the target that may be. Assumes that dates in the netcdf files are
#"      relative to 1990-01-01 (this is not checked).
#" @param layer Either a (positive) integer layer number, a negative number indicating depth below MSL (not depth below the free surface) 
#"        or "surface" to choose the surface layer. Defaults to "surface".
#" @param Land_map Set to TRUE to show a map of Queensland as an underlay for the model output plot. No longer requires the fgmap library 
#"      and an activated Google API key but also doesn't show maps for other locations
#"      Default now FALSE.
#" @param input_file is the URL or file location of any of the EMS output files or a THREDDS catalog URI. 
#"        Defaults to a menu selection based on current NCI catalogs. Can also be set to "nci", "menu" or "catalog" for the same behaviour.
#"        Set to "old_menu" to provide old menu options instead of menu options from the NCI catalog.
#"        Numeric values are interpreted as references to selections available from the old menu.
#"        Short codes can be used for some options (codenames as used in https://research.csiro.au/ereefs/models/model-outputs/access-to-raw-model-output/ )
#" @param input_grid Name of the locally-stored or opendap-served netcdf file that contains the grid
#"      coordinates for the cell corners (x_grid and y_grid). If not specified, the function will first look for
#"      x_grid and y_grid can be found in the first INPUT_STEM file, and if not found, will check whether the size 
#"      of the variables in the input file corresponds to the size expected for GBR4 or GBR1, and load appropriate 
#"      x and y grids from data files stored in this package. Alternatively, you can provide the location of a full 
#"      (not simple-format) ereefs netcdf output file such as 
#"      "http://dapds00.nci.org.au/thredds/dodsC/fx3/gbr4_hydro_all/gbr4_all_2016-09.nc"
#" @param scale_col Vector of colours to use for low and high values in the colour scale. This can be a colour 
#"      from the ggplot colour palette or a RGB hash code, or "spectral". Ignored for true_colour plots. 
#"      If one value is given (other than "spectral"), low colour is set to ivory and high colour to the value given.
#"      If three values are given, uses scale_fill_gradient2 (spectrum from low to high through middle value).
#"      Defaults to c("ivory", "coral4").
#" @param scale_lim Upper and lower bounds for colour scale. Defaults to full range of data.
#"      Ignored for true_colour plots.
#" @param box_bounds Minimum and maximum latitude and longitude coordinates to map. Defaults to the
#"        entire extent of the model output (though modified by the value of zoom). 
#"        Format: c(longitude_min, longitude_max, latitude_min, latitude_max).
#" @param p Handle for an existing figure if you want to add a layer instead of creating a new figure.
#"        If p is provided, Land_map is over-ridden and set to FALSE.
#" @param suppress_print Set to TRUE if you don't want the plots generated and saved. Defaults to TRUE.
#" @param return_poly Instead of only the figure handle, return a list containing the figure handle and the dataframe used by geom_plot(). Default FALSE.
#" @param label_towns Add labels for town locations to the figure. Default TRUE
#" @param strict_bounds Obsolescent: ignored
#" @param mark_points Data frame containing longitude and latitude of geolocations to mark with crosses (or a vector containing one location). Default NULL.
#" @export
#" @examples
#" \dontrun{
#" map_ereefs()
#" map_ereefs("TN")
#" map_ereefs("plume", target_date = c(2011, 1, 30))
#" map_ereefs("Chl_a_sum", target_date = "2016-07-15", scale_col = c("ivory", "green4"))
#" map_ereefs("salt", box_bounds = c(145, 150, -20, -15), zoom = 7, scale_lim = c(32, 35))
#"}
map_ereefs <- function(var_name = "true_colour", 
                       target_date = c(2018, 1, 30), 
                       layer = "surface", 
                       Land_map = FALSE,
                       input_file = "catalog",
                       input_grid = NA,
                       scale_col = c("ivory", "coral4"), 
                       scale_lim =c(NA, NA),
                       zoom = 6, 
                       box_bounds = c(NA, NA, NA, NA), 
                       p = NA, 
                       suppress_print = TRUE, 
                       return_poly = FALSE,
                       label_towns = TRUE,
                       strict_bounds = FALSE,
                       mark_points = NULL,
                       gbr_poly = FALSE) {
  assignList(get_params(target_date, target_date, input_file, var_name))

if (length(p)!= 1) Land_map <- FALSE # Don't add a land map if we are adding to an existing plot
if (suppress_print) Land_map <- FALSE # No need for a land map if we just want numbers, not an image

towns <- data.frame(latitude = c(-15.47027987, -16.0899, -16.4840, -16.92303816, -19.26639219, -20.0136699, -20.07670986, -20.40109791, -21.15345122, -22.82406858, -23.38031858, -23.84761069, -24.8662122, -25.54073075, -26.18916037),
                    longitude = c(145.2498605, 145.4622, 145.4623, 145.7710, 146.805701, 148.2475387, 146.2635394, 148.5802016, 149.1655418, 147.6363616, 150.5059485, 151.256349, 152.3478987, 152.7049316, 152.6581893),
                    town = c("Cooktown", "Cape Tribulation", "Port Douglas", "Cairns", "Townsville", "Bowen", "Charters Towers", "Prosperine", "Mackay", "Clermont", "Rockhampton", "Gladstone", "Bundaberg", "Maryborough", "Gympie"))

#check_platform_ok(input_stem)
grids <- get_ereefs_grids(input_file, input_grid)
x_grid <- grids[["x_grid"]]
y_grid <- grids[["y_grid"]]

# Allow user to specify a depth below MSL by setting layer to a negative value
if (layer <= 0) {
  z_grid <- grids[["z_grid"]]
  layer <- max(which(z_grid<layer))
}

day <- which.min(abs(start_date - ds))
if (min(abs(start_date-ds))>1) warning(paste("Target date", start_date, "is", min(abs(start_date-ds)), 
    "days from closest available date in", input_file, ds[min(abs(start_dat-ds))]))

# Allow for US English:
if (var_name == "true_color") {
	var_name <- "true_colour"
}

# Get cell grid corners
dims <- dim(x_grid) - 1

outOfBox <- array(FALSE, dim = dim(x_grid))
if (!is.na(box_bounds[1])) {
  outOfBox <- apply(x_grid, 2, function(x){ (x<box_bounds[1]|is.na(x)) } )
}
if (!is.na(box_bounds[2])) {
  outOfBox <- outOfBox | apply(x_grid, 2, function(x){(x>box_bounds[2]|is.na(x))})
}
if (!is.na(box_bounds[3])) {
  outOfBox <- outOfBox | apply(y_grid, 2, function(x){(x<box_bounds[3]|is.na(x))})
}
if (!is.na(box_bounds[4])) {
  outOfBox <- outOfBox | apply(y_grid, 2, function(x){(x>box_bounds[4]|is.na(x))})
}

if (is.na(box_bounds[1])) { 
  xmin <- 1
} else {
  xmin <- which(apply(!outOfBox, 1, any))[1]
  if (length(xmin) == 0) xmin <- 1
}
if (is.na(box_bounds[2])) {
  xmax <- dims[1]
} else {
  xmax <- which(apply(!outOfBox, 1, any))
  xmax <- xmax[length(xmax)]
  if ((length(xmax) == 0)|(xmax > dims[1])) xmax <- dims[1]
}
if (is.na(box_bounds[3])) { 
  ymin <- 1
} else {
  ymin <- which(apply(!outOfBox, 2, any))[1]
  if (length(ymin) == 0) ymin <- 1
}
if (is.na(box_bounds[4])) {
  ymax <- dims[2]
} else {
  ymax <- which(apply(!outOfBox, 2, any))
  ymax <- ymax[length(ymax)]
  if ((length(ymax) == 0)|(ymax > dims[2])) ymax <- dims[2]
}

x_grid <- x_grid[xmin:(xmax+1), ymin:(ymax+1)]
y_grid <- y_grid[xmin:(xmax+1), ymin:(ymax+1)]


# There are a few options for var_name. In most cases, we just plot the variable with the given name from the netcdf file.
# There are two special cases using optical output data: "plume" calculates and plots the optical colour class of the water,
# while "true_colour" produces something that looks like a satellite image from the model output
if (var_name == "plume") {
    nc <- safe_nc_open(input_file)
    R_412 <- safe_ncvar_get(nc, "R_412", start = c(xmin, ymin, day), count = c(xmax-xmin, ymax-ymin, 1))
    R_443 <- safe_ncvar_get(nc, "R_443", start = c(xmin, ymin, day), count = c(xmax-xmin, ymax-ymin, 1))
    R_488 <- safe_ncvar_get(nc, "R_488", start = c(xmin, ymin, day), count = c(xmax-xmin, ymax-ymin, 1))
    R_531 <- safe_ncvar_get(nc, "R_531", start = c(xmin, ymin, day), count = c(xmax-xmin, ymax-ymin, 1))
    R_547 <- safe_ncvar_get(nc, "R_547", start = c(xmin, ymin, day), count = c(xmax-xmin, ymax-ymin, 1))
    R_667 <- safe_ncvar_get(nc, "R_667", start = c(xmin, ymin, day), count = c(xmax-xmin, ymax-ymin, 1))
    R_678 <- safe_ncvar_get(nc, "R_678", start = c(xmin, ymin, day), count = c(xmax-xmin, ymax-ymin, 1))
    rsr <- list(R_412, R_443, R_488, R_531, R_547, R_667, R_678)
    ems_var <- plume_class(rsr)
    dims <- dim(ems_var)
    var_units <- ""
    var_longname <- "Plume colour class"

} else if (var_name == "true_colour") {
    nc <- safe_nc_open(input_file)
    TCbright <- 10
    R_470 <- safe_ncvar_get(nc, "R_470", start = c(xmin, ymin, day), count = c(xmax-xmin, ymax-ymin, 1)) * TCbright
    R_555 <- safe_ncvar_get(nc, "R_555", start = c(xmin, ymin, day), count = c(xmax-xmin, ymax-ymin, 1)) * TCbright
    R_645 <- safe_ncvar_get(nc, "R_645", start = c(xmin, ymin, day), count = c(xmax-xmin, ymax-ymin, 1)) * TCbright
    R_470[R_470>1] <- 1
    R_555[R_555>1] <- 1
    R_645[R_645>1] <- 1

    unscaledR = c(0, 30, 60, 120, 190, 255)/255;
    scaledR = c(1, 110, 160, 210, 240, 255)/255;
    scalefun <- approxfun(x = unscaledR, y = scaledR, yleft = 1, yright = 255)
    red <- scalefun(R_645)
    green <-scalefun(R_555)
    blue <- scalefun(R_470)
    red[is.na(red)] <- 0
    green[is.na(green)] <- 0
    blue[is.na(blue)] <- 0

    ems_var <- rgb(red, green, blue)
    ems_var[ems_var == "#000000"] <- NA
    ems_var <-array(as.character(ems_var), dim = dim(R_645))
    dims <- dim(ems_var)
    var_longname <- "Simulated true colour"
    var_units <- ""
} else if (var_name == "ZooT") {
    nc <- safe_nc_open(input_file)
    # We don't yet know the dimensions of the variable, so let"s get them
    dims <- nc$var[["ZooL_N"]][["size"]]
    if (is.null(dims)) stop(paste("ZooL_N", " not found in netcdf file.")) 
    ndims <- length(dims)
    if ((ndims > 3) && (layer == "surface")) layer <- dims[3]
    var_longname <- "Total zooplankton nitrogen"
    var_units <- "mg N m-3"
} else if (var_name == "speed") {
    nc <- safe_nc_open(input_file)
    # We don't yet know the dimensions of the variable, so let"s get them 
    dims <- nc$var[["u"]][["size"]] 
    if (is.null(dims)) stop(paste("u", " not found in netcdf file.")) 
    ndims <- length(dims) 
    if ((ndims > 3) && (layer == "surface")) layer <- dims[3]
    var_longname <- "Current speed"
    var_units <- "m s-1"
} else { 
    # This is a bit roundabout but we first need to activate by variable name, then check which grid is active and activate the whole grid
    nc <- tidync::tidync(input_file) %>% tidync::activate(var_names[1]) 
    nc <- nc %>% tidync::activate(tidync::active(nc), select_var = var_names)
    # We don't yet know the dimensions of the variable, so let"s get them 
    ndims <- dim(nc %>% tidync::hyper_dims())[1]
    if ((ndims > 3) && (layer == "surface")) layer <- dims[3]
}

if (var_name == "ZooT") {
    var_longname <- "Total Zooplankton Nitrogen"
    var_units <- "mg N m3"
    if (ndims == 4) {
       ems_var <- safe_ncvar_get(nc, "ZooL_N", start = c(xmin, ymin, layer, day), count = c(xmax-xmin, ymax-ymin, 1, 1))
       ems_var <- ems_var + safe_ncvar_get(nc, "ZooS_N", start = c(xmin, ymin, layer, day), count = c(xmax-xmin, ymax-ymin, 1, 1))
    } else {
       ems_var <- safe_ncvar_get(nc, "ZooL_N", start = c(xmin, ymin, day), count = c(xmax-xmin, ymax-ymin, 1))
       ems_var <- ems_var + safe_ncvar_get(nc, "ZooS_N", start = c(xmin, ymin, day), count = c(xmax-xmin, ymax-ymin, 1))
    }
} else if (var_name == "speed") {
    var_longname <- "Current speed"
    var_units <- "m s-1"
    if (ndims == 4) {
       ems_var <- safe_ncvar_get(nc, "u1", start = c(xmin, ymin, layer, day), count = c(xmax-xmin, ymax-ymin, 1, 1))
       ems_var <- sqrt(ems_var^2 + safe_ncvar_get(nc, "u2", start = c(xmin, ymin, layer, day), count = c(xmax-xmin, ymax-ymin, 1, 1))^2)
    } else {
       ems_var <- safe_ncvar_get(nc, "u1", start = c(xmin, ymin, day), count = c(xmax-xmin, ymax-ymin, 1))
       ems_var <- ems_var + safe_ncvar_get(nc, "u2", start = c(xmin, ymin, day), count = c(xmax-xmin, ymax-ymin, 1))
    }
} else if (!((var_name == "true_colour") || (var_name == "plume"))) {
    vat <- ncdf4::ncatt_get(nc, var_name)
    var_longname <- vat$long_name
    var_units <- vat$units
    if (ndims == 4) {
       ems_var <- safe_ncvar_get(nc, var_name, start = c(xmin, ymin, layer, day), count = c(xmax-xmin, ymax-ymin, 1, 1))
    } else {
       ems_var <- safe_ncvar_get(nc, var_name, start = c(xmin, ymin, day), count = c(xmax-xmin, ymax-ymin, 1))
    }
}

ncdf4::nc_close(nc)

nc <- safe_nc_open(input_file)
if (!is.null(nc$var[["botz"]])) {
  botz <- safe_ncvar_get(nc, "botz", start = c(xmin, ymin), count = c(xmax-xmin, ymax-ymin))
} else {
  botz <- array(NA, dim(ems_var))
}
if (is.null(nc$var[["latitude"]])) {
# Standard EMS output file
  latitude <- safe_ncvar_get(nc, "y_centre", start = c(xmin, ymin), count = c(xmax-xmin, ymax-ymin))
  longitude <- safe_ncvar_get(nc, "x_centre", start = c(xmin, ymin), count = c(xmax-xmin, ymax-ymin))
} else { 
  # Simple format netcdf file
  latitude <- safe_ncvar_get(nc, "latitude", start = c(xmin, ymin), count = c(xmax-xmin, ymax-ymin))
  longitude <- safe_ncvar_get(nc, "longitude", start = c(xmin, ymin), count = c(xmax-xmin, ymax-ymin))
}
ncdf4::nc_close(nc)

a <- dim(ems_var)[1]
b <- dim(ems_var)[2]

# Set up the polygon corners. 4 per polygon.
gx <- c(x_grid[1:a, 1:b], x_grid[2:(a+1), 1:b], x_grid[2:(a+1), 2:(b+1)], x_grid[1:a, 2:(b+1)])
gy <- c(y_grid[1:a, 1:b], y_grid[2:(a+1), 1:b], y_grid[2:(a+1), 2:(b+1)], y_grid[1:a, 2:(b+1)])
gx <- array(gx, dim = c(a*b, 4))
gy <- array(gy, dim = c(a*b, 4))

# Find and exclude points where not all corners are defined
gx_ok <- !apply(is.na(gx), 1, any)
gy_ok <- !apply(is.na(gy), 1, any)

# Values associated with each polygon at chosen timestep
#n <- c(ems_var[, ,tstep])[gx_ok&gy_ok]
n <- c(ems_var)[gx_ok&gy_ok]
# (gx_ok and gy_ok should be identical, but let"s be certain)
gx <- c(t(gx[gx_ok&gy_ok, ]))
gy <- c(t(gy[gx_ok&gy_ok, ]))
longitude <- c(longitude)[gx_ok&gy_ok]
latitude <- c(latitude)[gx_ok&gy_ok]
botz <- c(botz)[gx_ok&gy_ok]

# Unique ID for each polygon
id <- 1:length(n)

id <- as.factor(id)
values <- data.frame(id = id, value = n, depth = botz, x_centre = longitude, y_centre = latitude)
positions <- data.frame(id = rep(id, each = 4), x = gx, y = gy)
datapoly <- merge(values, positions, by = c("id"))

if ((var_name!= "true_colour")&&(is.na(scale_lim[1]))) { 
	scale_lim <- c(min(n, na.rm = TRUE), max(n, na.rm = TRUE))
}

if (Land_map) {
  MapLocation<-c(min(gx, na.rm = TRUE)-0.5, 
 		min(gy, na.rm = TRUE)-0.5, 
 		max(gx, na.rm = TRUE)+0.5, 
 		max(gy, na.rm = TRUE)+0.5)
  p <- ggplot2::ggplot() +
           ggplot2::geom_polygon(data = map.df, colour = "black", fill = "lightgrey", size = 0.5, ggplot2::aes(x = long, y = lat, group = group))
} else if (length(p) == 1) {
  p <- ggplot2::ggplot()
}

p <-  p + ggplot2::geom_polygon(ggplot2::aes(x = x, y = y, fill = value, group = id), data = datapoly) 

if (var_name == "true_colour") {
  p <- p + ggplot2::scale_fill_identity() 
} else if (scale_col[1] == "spectral") { 
  p <- p + ggplot2::scale_fill_distiller(palette = "Spectral",
                                         na.value = "transparent", 
                                         guide = "colourbar", 
                                         limits = scale_lim, 
                                         name = var_units, 
                                         oob = scales::squish) 
} else if (length(scale_col)<3) { 
  if (length(scale_col) == 1) scale_col <- c("ivory", scale_col) 
  p <- p + ggplot2::scale_fill_gradient(low = scale_col[1], 
                                        high = scale_col[2], 
				                                na.value = "transparent", 
				                                guide = "colourbar",
				                                limits = scale_lim,
				                                #name = expression(paste("cells m"^"-3")),
				                                name = var_units,
				                                oob = scales::squish) 
} else { 
  p <- p + ggplot2::scale_fill_gradient2(low = scale_col[1], 
                                         mid = scale_col[2], 
                                         high = scale_col[3], 
                                         na.value = "transparent", 
                                         guide = "colourbar",
     			                               limits = scale_lim,
                                         midpoint = (scale_lim[2] - scale_lim[1])/2,
                                         space = "Lab",
      		                               #name = expression(paste("cells m"^"-3")),
  	    	                               name = var_units,
  		                                   oob = scales::squish) 
}

if (label_towns) {
  towns <- towns[towns$latitude >= min(gy, na.rm = TRUE), ]
  towns <- towns[towns$latitude <= max(gy, na.rm = TRUE), ]
  towns <- towns[towns$longitude >= min(gx, na.rm = TRUE), ]
  towns <- towns[towns$longitude <= max(gx, na.rm = TRUE), ]
  if (dim(towns)[1]>0) p <- p + ggplot2::geom_text(data = towns, ggplot2::aes(x = longitude, y = latitude, label = town, hjust = "right"), nudge_x = -0.1) +
                                 ggplot2::geom_point(data = towns, ggplot2::aes(x = longitude, y = latitude))
}

p <- p + ggplot2::ggtitle(paste(var_longname, format(chron::chron(as.numeric(ds[day])+0.000001), "%Y-%m-%d %H:%M"))) +
    ggplot2::xlab("longitude") + ggplot2::ylab("latitude")
if (!is.null(mark_points)) {
  if (is.null(dim(mark_points))) mark_points <- data.frame(latitude = mark_points[1], longitude = mark_points[2])
  p <- p + ggplot2::geom_point(data = mark_points, ggplot2::aes(x = longitude, y = latitude), shape = 4)
}
if (gbr_poly) {
  p <- p + ggplot2::geom_path(data = sdf.gbr, ggplot2::aes(y = lat, x = long, group = group))
}
if (all(is.na(box_bounds))) { 
  p <- p + ggplot2::coord_map(xlim = c(min(gx, na.rm = TRUE), max(gx, na.rm = TRUE)), ylim = c(min(gy, na.rm = TRUE), max(gy, na.rm = TRUE)))
} else { 
  p <- p + ggplot2::coord_map(xlim = box_bounds[1:2], ylim = box_bounds[3:4]) + 
    ggplot2::theme(panel.border = ggplot2::element_rect(linetype = "solid", colour = "grey", fill = NA))
} 

if (!suppress_print) print(p)
if (return_poly) {
  return(list(p = p, datapoly = datapoly, longitude = longitude, latitude = latitude))
} else {
  return(p)
}
}

#" Create a series of map image files for an animation of eReefs model output AND calculate temporal
#" mean values.
#"
#" Creates and saves to disk a sequential series of colour map images showing concentrations of a specified 
#" eReefs model output variable at a specified model layer (by default, the surface layer). 
#" ALSO CALCULATES THE TEMPORAL MEAN value of each cell over the specified time (visualisation of maps can
#" be suppressed by setting suppress_print to TRUE if this is the primary desired output).
#" Maps produced are optionally overlain on a map of Queensland.
#" Can be more efficient than calling map_ereefs multiple times if you 
#" want to produce an animation because it loads a month at a time for GBR4 runs (unless selected via catalog). 
#" If output files contain multiple outputs per day, chooses the step closest to midday and uses only daily output.
#" To stitch together the images into an animation, you will need other software such as ImageMagick (recommended)
#" or imageJ.  Barbara Robson (AIMS).
#"
#" TODO: Update to use chron dates as has been done for functions in data_extraction_functions.R
#"
#" @param var_name Short name of the variable to plot. This can be any variable in the
#"                 eReefs netcdf file that you are accessing (refer to eReefs model
#"                 documentation or extract variable names from the netcdf file for a full 
#"                 list). In addition, two special var_name values are supported: "true_colour"
#"                 and "plume". "true_colour" provides a simulated MODIS true colour
#"                 image (refer to Baird et al., 2016 for en explanation)."plume"
#"                 provides a map of calculated plume colour class as per Devlin et al. (2012)
#"                 and Baird et al. (2017). Defaults to true_colour.
#" @param start_date date  for animation. Can either be a) a vector of integers in the format 
#"	c(year, month, day, hour), c(year, month, day) or (year, month); b) a date obtained e.g. 
#"  from as.Date(); or c) a character string formatted for input to as.Date(). Defaults to c(2015, 2, 1).
#" @param end date for animation. Can either be a) a vector of integers in the format 
#"	c(year, month, day); b) a date obtained e.g. from as.Date(); or c) a character string 
#"      formatted for input to as.Date(). Defaults to c(2016, 3, 31).
#" @param layer Either an integer layer number or "surface" to choose the surface layer. Defaults to "surface".
#" @param output_dir Path to directory in which to store images for animation. Created if necessary. Defaults
#"      to "ToAnimate". Images are created in this directory with input_files beginning with var_name, 
#"      followed by an underscore and then sequential numbers beginning with 100001.
#" @param Land_map Set to TRUE to show a land map of Queensland. Default now FALSE.
#" @param input_file is the URL or file location of any of the EMS output files or a THREDDS catalog URI. 
#"        Defaults to a menu selection based on current NCI catalogs. Can also be set to "nci", "menu" or "catalog" for the same behaviour.
#"        Set to "old_menu" to provide old menu options instead of menu options from the NCI catalog.
#"        Numeric values are interpreted as references to selections available from the old menu.
#"        Short codes can be used for some options (codenames as used in https://research.csiro.au/ereefs/models/model-outputs/access-to-raw-model-output/ )
#" @param input_grid Name of the locally-stored or opendap-served netcdf file that contains the grid
#"      coordinates for the cell corners (x_grid and y_grid). If not specified, the function will first look for
#"      x_grid and y_grid can be found in the first INPUT_STEM file, and if not found, will check whether the size 
#"      of the variables in the input file corresponds to the size expected for GBR4 or GBR1, and load appropriate 
#"      x and y grids from data files stored in this package. Alternatively, you can provide the location of a full 
#"      (not simple-format) ereefs netcdf output file such as 
#"      "http://dapds00.nci.org.au/thredds/dodsC/fx3/gbr4_hydro_all/gbr4_all_2016-09.nc"
#" @param scale_col Vector of colours to use for the colour scale. This can be colours 
#"      from the ggplot colour palette or a RGB hash code. Ignored for true_colour plots. 
#"      If set to "spectral", uses a colour spectrum from bluish to red (similar to jet but less vivid). Otherwise:
#"      If one value is given, low colour is set to ivory and high colour to the value given.
#"      If two values are given, these are used as low and high limit colours.
#"      If three values are given, the middle value is used to set the mid-point of the scale.
#"      Defaults to c("ivory", "coral4").
#" @param box_bounds Minimum and maximum latitude and longitude coordinates to map. Defaults to the
#"        entire extent of the model output (though modified by the value of zoom). 
#"        Format: c(longitude_min, longitude_max, latitude_min, latitude_max). It is recommended to
#"        also specify an appropriate value for zoom if specifying box_bounds.
#" @param suppress_print Set to TRUE if you don't want the plots generated and saved. Defaults to TRUE.
#" @param stride Default "daily", but can otherwise be set to a numeric interval indicating how many time-steps to step forward for each frame.
#" @param verbosity Set 0 for just a waitbar, 1 for more updates, 2 for debugging information. Default 0.
#" @param label_towns Add labels for town locations to the figure. Default TRUE
#" @param strict_bounds Obsolescent: ignored
#" @param mark_points Data frame containing longitude and latitude of geolocations to mark with crosses (or a vector containing one location). Default NULL.
#" @param gbr_poly TRUE to show contours of approximate reef areas. Default FALSE.
#" @param add_arrows TRUE to show arrows indicating magnitude and direction of flow. Default FALSE.
#" @param max_u Velocity at which to show maximum arrow length. Default NA, in which case it will use the maximum observed velocity.
#" @param scale_arrows Factor by which to scale arrows. Values >1 result in longer arrows. Values <1 result in shorter arrows. Default 1.
#" @param show_bathy TRUE to show contours based on the bathymetry as represented in the model. Default FALSE. Requires model file to contain botz (this 
#"        requirement may be dropped in future versions for GBR1 and GBR4 runs).
#" @param contour_breaks Depths in metres to show with show_bathy. Default c(5, 10, 20).
#" @return a list that includes data.frame formatted for use in ggplot2::geom_polygon, containing a map of the temporally averaged
#"         value of the variable specified in VAR_NAME over the selected interval, plus the actual values and cell centre latitudes and longitudes.
#" @export
#" @examples
#" \dontrun{
#" map_ereefs_movie(start_date = c(2016, 2, 1), end_date = c(2016, 2, 15))
#"}
map_ereefs_movie <- function(var_name = "true_colour", 
                             start_date = c(2015, 12, 1), 
                             end_date = c(2016, 3, 31), 
                             layer = "surface", 
                             output_dir = "ToAnimate", 
                             Land_map = FALSE, 
                             input_file = "catalog",
                             input_grid = NA, 
                             scale_col = c("ivory", "coral4"), 
                             scale_lim = c(NA, NA), 
                             zoom = 6, 
                             box_bounds = c(NA, NA, NA, NA), 
                             suppress_print = TRUE,
                             stride = "daily",
                             verbosity = 0, 
                             label_towns = TRUE,
                             strict_bounds = FALSE,
                             mark_points = NULL,
                             gbr_poly = FALSE,
                             add_arrows = FALSE,
                             max_u = NA,
                             scale_arrows = NA,
                             show_bathy = FALSE,
                             contour_breaks = c(5, 10, 20)) {
  plot_eta <- FALSE
  # Get parameter values and assign results from returned list to relevant variable names
  # This assigns input_file, ereefs_case, input_stem, start_date, end_date, start_tod, start_month, start_year,
  # end_date, end_day, end_month, end_year, mths, years, var_list, ereefs_origin and blank_length
  assignList(get_params(start_date, end_date, input_file, var_name))

  if (verbosity > 1) print(paste("After substitute_filename() input_file = ", input_file))

  if (ereefs_case[2] == "1km") warning("Assuming that only one timestep is output per day/file") # find matching commented warning to fix this
  #check_platform_ok(input_stem)

  towns <- data.frame(latitude = c(-15.47027987, -16.0899, -16.4840, -16.92303816, -19.26639219, -20.0136699, -20.07670986, -20.40109791, 
                                   -21.15345122, -22.82406858, -23.38031858, -23.84761069, -24.8662122, -25.54073075, -26.18916037),
                    longitude = c(145.2498605, 145.4622, 145.4623, 145.7710, 146.805701, 148.2475387, 146.2635394, 148.5802016, 149.1655418, 
                                  147.6363616, 150.5059485, 151.256349, 152.3478987, 152.7049316, 152.6581893),
                    town = c("Cooktown", "Cape Tribulation", "Port Douglas", "Cairns", "Townsville", "Bowen", "Charters Towers", 
                            "Prosperine", "Mackay", "Clermont", "Rockhampton", "Gladstone", "Bundaberg", "Maryborough", "Gympie"))

  # Points of interest provided by the user to mark on the map, and possibly plot a surface elevation timeseries for:
  if (!is.null(mark_points)) {
  # If mark_points is a vector, change it into a data frame
  if (is.null(dim(mark_points))) {
     mark_points <- data.frame(latitude = mark_points[1], longitude = mark_points[2])
  }
  eta_data <- get_ereefs_ts(var_name = "eta", input_file = input_file, start_date = start_date, end_date = end_date, location_latlon = mark_points)
  names(eta_data) <- c("date", "eta")
  eta_plot <- ggplot2::ggplot(eta_data, ggplot2::aes(x = date, y = eta)) + ggplot2::geom_line() + ggplot2::ylab("surface elevation (m)")
  plot_eta <- TRUE
  }

  if (suppress_print) Land_map <- FALSE
  
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
         input_file <- paste0(input_stem, format(as.Date(paste(year, month, 1, sep = "-")), "%Y-%m"), ".nc")
       } else if (ereefs_case[2] == "1km") {
         input_file <- paste0(input_stem, format(as.Date(paste(year, month, from_day, sep = "-")), "%Y-%m-%d"), ".nc")
       } else { #recom or other netcdf or ncml file
         #input_file <- input_file
         ds<- get_origin_and_times(input_file)[[2]]
         ds_original <- ds
       }
       grids <- get_ereefs_grids(input_file, input_grid)
       x_grid <- grids[["x_grid"]]
       y_grid <- grids[["y_grid"]]

       # Allow user to specify a depth below MSL by setting layer to a negative value
       if (layer <= 0) {
          z_grid <- grids[["z_grid"]]
          layer <- max(which(z_grid<layer))
       }

       dims <- dim(x_grid) - 1

       # Work out which parts of the grid are within box_bounds and which are outside
       outOfBox <- array(FALSE, dim = dim(x_grid))
       if (!is.na(box_bounds[1])) {
         outOfBox <- apply(x_grid, 2, function(x){ (x<box_bounds[1]|is.na(x)) } )
       }
       if (!is.na(box_bounds[2])) {
         outOfBox <- outOfBox | apply(x_grid, 2, function(x){(x>box_bounds[2]|is.na(x))})
       }
       if (!is.na(box_bounds[3])) {
         outOfBox <- outOfBox | apply(y_grid, 2, function(x){(x<box_bounds[3]|is.na(x))})
       }
       if (!is.na(box_bounds[4])) {
         outOfBox <- outOfBox | apply(y_grid, 2, function(x){(x>box_bounds[4]|is.na(x))})
       }
         
       # Find the subset of x_grid and y_grid that is inside the box and crop the grids
       # to the box_bounds
       if (is.na(box_bounds[1])) { 
        xmin <- 1
       } else {
        xmin <- which(apply(!outOfBox, 1, any))[1]
        if (length(xmin) == 0) xmin <- 1
       }
       if (is.na(box_bounds[2])) {
        xmax <- dims[1]
       } else {
         xmax <- which(apply(!outOfBox, 1, any))
         xmax <- xmax[length(xmax)]
         if ((length(xmax) == 0)|(xmax > dims[1])) xmax <- dims[1]
       }
       if (is.na(box_bounds[3])) { 
         ymin <- 1
       } else {
         ymin <- which(apply(!outOfBox, 2, any))[1]
         if (length(ymin) == 0) ymin <- 1
       }
       if (is.na(box_bounds[4])) {
         ymax <- dims[2]
       } else {
         ymax <- which(apply(!outOfBox, 2, any))
         ymax <- ymax[length(ymax)]
         if ((length(ymax) == 0)|(ymax > dims[2])) ymax <- dims[2]
       }

       x_grid <- x_grid[xmin:(xmax+1), ymin:(ymax+1)]
       y_grid <- y_grid[xmin:(xmax+1), ymin:(ymax+1)]

       nc <- safe_nc_open(input_file)
       if (is.null(nc$var[["latitude"]])) {
       # Standard EMS output file
         latitude <- safe_ncvar_get(nc, "y_centre")
         longitude <- safe_ncvar_get(nc, "x_centre")
         botz <- safe_ncvar_get(nc, "botz")
       } else { 
         # Simple format netcdf file
         latitude <- safe_ncvar_get(nc, "latitude")
         longitude <- safe_ncvar_get(nc, "longitude")
         botz <- NULL
         if (show_bathy) warning("Can not show bathymetry: simple format netcdf file does not contain botz")
       }
       ncdf4::nc_close(nc)

       if (add_arrows) {
         idim <- dim(latitude)[1]
         jdim <- dim(latitude)[2]
         max_arrow <- max(max(abs(longitude[idim, jdim] - longitude[1, 1])/idim), max(abs(latitude[idim, jdim] - latitude[1, 1])/jdim))
       }

       # Set up the polygon corners. 4 per polygon.
       a <- xmax - xmin + 1
       b <- ymax - ymin + 1
     
       gx <- c(x_grid[1:a, 1:b], x_grid[2:(a+1), 1:b], x_grid[2:(a+1), 2:(b+1)], x_grid[1:a, 2:(b+1)])
       gy <- c(y_grid[1:a, 1:b], y_grid[2:(a+1), 1:b], y_grid[2:(a+1), 2:(b+1)], y_grid[1:a, 2:(b+1)])
       gx <- array(gx, dim = c(a*b, 4))
       gy <- array(gy, dim = c(a*b, 4))

       # Find and exclude points where not all corners are defined
       gx_ok <- !apply(is.na(gx), 1, any)
       gy_ok <- !apply(is.na(gy), 1, any)
       gx <- c(t(gx[gx_ok & gy_ok, ]))
       gy <- c(t(gy[gx_ok & gy_ok, ]))
       longitude <- c(longitude)[gx_ok & gy_ok]
       latitude <- c(latitude)[gx_ok & gy_ok]

    } else {
       from_day <- 1
       start_tod <- 0
    }

    if ((start_year == end_year)&&(start_month == end_month)) {
       day_count <- end_day - start_day + 1
    } else if (mcount == 1) {
       day_count <- daysIn(as.Date(paste(year, month, 1, sep = "-"))) - start_day + 1
    } else if (mcount == (length(mths))) {
       day_count <- end_day
    } else {
       day_count <- daysIn(as.Date(paste(year, month, 1, sep = "-")))
    }

    if (ereefs_case[2] == "4km") { 
      fileslist <- 1
	    input_file <- paste0(input_stem, format(as.Date(paste(year, month, 1, sep = "-")), "%Y-%m"), ".nc")
      ds <- get_origin_and_times(input_file)[[2]]
      if ((ds[length(as.numeric(ds))] - as.numeric(ds)[1]) > 31.5) {
         warning("Filename looks like a monthly output file (i.e. contains two dashes) but file contains more than a month of data.")
      }
      if ((ds[length(as.numeric(ds))] - as.numeric(ds)[1]) < 27) {
        warning("Filename looks like a monthly output file (i.e. contains two dashes) but file contains less than a month of data.")
      }
      if (ds[2] == ds[1]) stop(paste("Error reading time from", input_file, "(t[2] == t[1])"))
      tstep <- as.numeric(ds[2] - ds[1])
      day_count <- day_count / tstep
      if (day_count > length(as.numeric(ds))) {
        warning(paste("end_date", end_date, "is beyond available data. Ending at", max(ds)))
        day_count <- length(as.numeric(ds))
      }
      #dum1 <- as.integer((from_day - 0.4999999)/tstep + 1)
      #dum2 <- as.integer((day_count - 1) / tstep) +1
      #ds <- ds[seq(from = dum1, by = as.integer(1/tstep), to = (dum1+dum2))]
      #start_array <- c(xmin, ymin, dum1)
	    #count_array <- c(xmax-xmin, ymax-ymin, dum2)
      ds <- ds[from_day:day_count]
      start_array <- c(xmin, ymin, from_day)
	    count_array <- c(xmax - xmin, ymax - ymin, day_count)
	    fileslist <- 1
    } else if (ereefs_case[2] == "1km") { 
	    fileslist <- from_day:(from_day+day_count-1)
	    from_day <- 1
	    day_count <- 1

	    input_file <- paste0(input_stem, format(as.Date(paste(year, month, from_day, sep = "-")), "%Y-%m-%d"), ".nc")
      ds <- get_origin_and_times(input_file, as_chron = "FALSE")[[2]]
      dum1 <- which.min(abs(as.numeric(ds - (from_day + 0.4999999))))
      start_array <- c(xmin, ymin, dum1) 
      count_array <- c(xmax-xmin, ymax-ymin, dum1)
      tstep <- 1
    } else if ((ereefs_case[2] == "recom")|(ereefs_case[1] == "ncml")) {
	    # Everything is in one file but we are only going to read a month at a time
	    # Output may be more than daily, or possibly less
	    # input_file <- paste0(input_stem, ".nc") # input_file has been set previously 
      tstep <- as.numeric(median((ds[2:length(ds)] - ds[1:(length(ds) - 1)]), na.rm = TRUE))
      day_count <- day_count / tstep
      if (day_count > length(as.numeric(ds_original))) {
        warning(paste("end_date", end_date, "is beyond available data. Ending at", max(ds_original)))
        day_count <- length(as.numeric(ds_original))
      }
      from_day <- (as.numeric(chron::chron(paste(year, month, from_day, sep = "-"), format = c("y-m-d"),
                                  origin = c(year = 1990, month = 1, day = 1)) - as.numeric(ds_original)[1]) + 
                              start_tod) / tstep + 1 
      #print(paste(year, month, from_day, ds[1], start_tod, tstep))
	    if (from_day<1) from_day <-1
	    start_array <- c(xmin, ymin, floor(from_day))
	    count_array <- c(xmax-xmin, ymax-ymin, if(tstep>1) {as.integer(day_count/tstep)} else  day_count)
      #print(start_array)
      #print(count_array)
	    fileslist <- 1
    } else stop("Shouldn't happen: ereefs_case not recognised")

    if (stride == "daily") {
      if (tstep <= 1.0) {
       stride <- 1/tstep
      } else {
        stop("Minimum timestep in netcdf file is greater than daily.")
      }
    } 
    stride <- as.integer(stride)

    for (dcount in 1:length(fileslist)) {
      if (ereefs_case[2] == "1km") { 
         input_file <- paste0(input_stem, format(as.Date(paste(year, month, fileslist[dcount], sep = "-")), "%Y-%m-%d"), ".nc") 
         #ds <- as.Date(paste(year, month, fileslist[dcount], sep = "-", "%Y-%m-%d"))
      }
      if (verbosity>0) print(input_file)
      # There are a few options for var_name. In most cases, we just plot the variable with the given name from the netcdf file.
      # There are two special cases using optical output data: "plume" calculates and plots the optical colour class of the water,
      # while "true_colour" produces something that looks like a satellite image from the model output
      if (var_name == "plume") {
        nc <- safe_nc_open(input_file)
        R_412 <- safe_ncvar_get(nc, "R_412", start = start_array, count = count_array)
        R_443 <- safe_ncvar_get(nc, "R_443", start = start_array, count = count_array)
        R_488 <- safe_ncvar_get(nc, "R_488", start = start_array, count = count_array)
        R_531 <- safe_ncvar_get(nc, "R_531", start = start_array, count = count_array)
        R_547 <- safe_ncvar_get(nc, "R_547", start = start_array, count = count_array)
        R_667 <- safe_ncvar_get(nc, "R_667", start = start_array, count = count_array)
        R_678 <- safe_ncvar_get(nc, "R_678", start = start_array, count = count_array)
        ems_var <- NA*R_678
        if (ereefs_case[2] == "4km") {
           for (day in 1:dim(R_412)[3]) {
             rsr <- list(R_412[, ,day], R_443[, ,day], R_488[, ,day], R_531[, ,day], R_547[, ,day], R_667[, ,day], R_678[, ,day])
             ems_var[, ,day] <- plume_class(rsr)
           }
        } else {
	         rsr <- list(R_412, R_433, R_488, R_532, R_547, R_668, R_678)
	         ems_var <- plume_class(rsr)
        }
        dims <- dim(ems_var)
	     var_longname <- "Plume optical class"
	     var_units <- ""
      } else if (var_name == "true_colour") {
        slice <- paste0("[", start_array[3]-1, ":", stride, ":", start_array[3] + count_array[3] - 2, "]", # time
                        "[", start_array[2]-1, ":", start_array[2] + count_array[2] - 1, "]", # y
                        "[", start_array[1]-1, ":", start_array[1] + count_array[1] - 1, "]") # x
        count_array <- c(count_array[1:2]+1, count_array[3])
        nc <- safe_nc_open(input_file)
        TCbright <- 10
        R_470 <- safe_ncvar_get(nc, "R_470", start = start_array, count = count_array) * TCbright
        R_555 <- safe_ncvar_get(nc, "R_555", start = start_array, count = count_array) * TCbright
        R_645 <- safe_ncvar_get(nc, "R_645", start = start_array, count = count_array) * TCbright
        R_470[R_470>1] <- 1
        R_555[R_555>1] <- 1
        R_645[R_645>1] <- 1
    
        unscaledR = c(0, 30, 60, 120, 190, 255)/255;
        scaledR = c(1, 110, 160, 210, 240, 255)/255;
        scalefun <- approxfun(x = unscaledR, y = scaledR, yleft = 1, yright = 255)
        red <- scalefun(R_645)
        green <-scalefun(R_555)
        blue <- scalefun(R_470)
        red[is.na(red)] <- 0
        green[is.na(green)] <- 0
        blue[is.na(blue)] <- 0
  
        ems_var <- rgb(red, green, blue)
        ems_var[ems_var == "#000000"] <- NA
        ems_var <-array(as.character(ems_var), dim = dim(R_645))
        dims <- dim(ems_var)

	     var_longname <- "Simulated true colour"
	     var_units <- ""
    
      } else { 

        if (ndims == 0) {
          nc <- safe_nc_open(input_file)
          # We don't yet know the dimensions of the variable, so let"s get them
          if (var_name == "speed") {
             dims <- nc$var[["u"]][["size"]]
             ndims <- nc$var[["u"]][["ndims"]]
          } else if (var_name == "ZooT") {
             dims <- nc$var[["ZooL_N"]][["size"]]
             ndims <- nc$var[["ZooL_N"]][["ndims"]]
          } else { 
             dims <- nc$var[[var_name]][["size"]]
             ndims <- nc$var[[var_name]][["ndims"]]
          }
          if (is.null(dims)) stop(paste(var_name, " not found in netcdf file.")) 
          #ndims <- length(dims)
          # If there"s only one layer (e.g. a surf.nc file) then we want to reduce ndims accordingly, but not if there is only one time-step
          #if ((length(dims[dims!= 1])!= ndims)&&(dims[length(dims)]!= 1)) ndims <- ndims - 1
          if ((ndims > 3) && (layer == "surface")) layer <- dims[3] + 1
          ncdf4::nc_close(nc)
        }

        if (verbosity>1) {
           print(paste("variable has ", ndims, "dimensions"))
           print("start_array = ")
           print(start_array)
           print("count_array = ")
           print(count_array)
        }

        if (ndims == 4) {
          slice <- paste0("[", start_array[3]-1, ":", stride, ":", floor(start_array[3] + count_array[3] - 2), "]", # time
                          "[", layer-1, "]",                                                    # layer
                          "[", start_array[2]-1, ":", floor(start_array[2] + count_array[2] - 1), "]", # y
                          "[", start_array[1]-1, ":", floor(start_array[1] + count_array[1] - 1), "]") # x
          start_array <- c(start_array[1:2], layer-1, start_array[3])
          count_array <- c(count_array[1:2]+1, 1, count_array[3])
        } else {
          slice <- paste0("[", start_array[3]-1, ":", stride, ":", floor(start_array[3] + count_array[3] - 2), "]", # time
                          "[", start_array[2]-1, ":", start_array[2] + count_array[2] - 1, "]", # y
                          "[", start_array[1]-1, ":", start_array[1] + count_array[1] - 1, "]") # x
          count_array <- c(count_array[1:2]+1, count_array[3])
        }
        if (verbosity>1) print(paste("slice = ", slice))
        if (var_name == "speed") { 
           nc <- safe_nc_open(input_file)
           ems_var <- sqrt(safe_ncvar_get(nc, "u", start = start_array, count = count_array)^2 + safe_ncvar_get(nc, "v", start = start_array, count = count_array)^2)
           vat <- ncdf4::ncatt_get(nc, "u")
           var_longname <- "Current speed"
           var_units <- vat$units
        } else if (var_name == "ZooT") {
           nc <- safe_nc_open(input_file)
           ems_var <- safe_ncvar_get(nc, "ZooL_N", start = start_array, count = count_array) + safe_ncvar_get(nc, "ZooS_N", start = start_array, count = count_array)
           vat <- ncdf4::ncatt_get(nc, "ZooL_N")
           var_longname <- "Total Zooplankton"
           var_units <- vat$units
        } else {
           if (verbosity>1) print(paste("Before nc_open, input_file = ", input_file))
           nc <- safe_nc_open(input_file)
           ems_var <- safe_ncvar_get(nc, var_name, start = start_array, count = count_array)
           if (add_arrows) {
             u_count_array <- c(count_array[1:2], 1, count_array[length(start_array)])
             u_start_array <- c(start_array[1:2], layer-1, start_array[length(start_array)])
             current_u <- safe_ncvar_get(nc, "u1", start = start_array, count = u_count_array)
             current_v <- safe_ncvar_get(nc, "u2", start = start_array, count = u_count_array)
             if ((dcount == 1)&(is.na(max_u))) {
               max_u <- max(max(abs(c(current_u)), na.rm = TRUE), max(abs(c(current_v)), na.rm = TRUE)) 
               if (!is.na(scale_arrows)) max_u <- max_u / scale_arrows
             }
           }
           vat <- ncdf4::ncatt_get(nc, var_name)
           var_longname <- vat$long_name
           var_units <- vat$units
        }
      }
    
      dum1 <- length(dim(ems_var))
      if (dum1 == 2) {
        ems_var <- array(ems_var, dim = c(dim(ems_var), 1))
        if (add_arrows) {
          current_u <- array(current_u, dim = c(dim(current_u), 1))
          current_v <- array(current_v, dim = c(dim(current_v), 1))
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
        if (icount == 0) {
           if (var_name == "true_colour") {
              temporal_sum <- col2rgb(n)^2
           } else { 
              temporal_sum <- n
           }
        } else {
           if (var_name == "true_colour") {
              temporal_sum <- col2rgb(n)^2 + temporal_sum
           } else { 
              temporal_sum <- temporal_sum + n
           }
        }
    
        # Unique ID for each polygon
        id <- as.factor(1:length(n))
        values <- data.frame(id = id, value = n)
        positions <- data.frame(id = rep(id, each = 4), x = gx, y = gy)
        datapoly <- merge(values, positions, by = c("id"))
    
        if (!suppress_print) {
            if ((var_name!= "true_colour")&&(is.na(scale_lim[1]))) { 
	            scale_lim <- c(min(n, na.rm = TRUE), max(n, na.rm = TRUE))
            }
  
            if (Land_map) {
               p <- ggplot2::ggplot() + ggplot2::geom_polygon(data = map.df, colour = "black", fill = "lightgrey", size = 0.5, ggplot2::aes(x = long, y = lat, group = group))
	          } else {
	             p <- ggplot2::ggplot()
            }
            p <- p + ggplot2::geom_polygon(ggplot2::aes(x = x, y = y, fill = value, group = id), data = datapoly)
            if (var_name == "true_colour") { 
              p <- p + ggplot2::scale_fill_identity()
            } else if (scale_col[1] == "spectral") { 
              p <- p + ggplot2::scale_fill_distiller(palette = "Spectral",
                                                     na.value = "transparent", 
                                                     guide = "colourbar", 
                                                     limits = scale_lim, 
                                                     name = var_units, 
                                                     oob = scales::squish) 
            } else if (length(scale_col)<3) { 
              if (length(scale_col) == 1) scale_col <- c("ivory", scale_col)
              p <- p + ggplot2::scale_fill_gradient(low = scale_col[1],
					                                          high = scale_col[2],
					                                          na.value = "transparent", 
					                                          guide = "colourbar",
					                                          limits = scale_lim,
					                                          name = var_units,
					                                          oob = scales::squish)
            } else {
                  p <- p + ggplot2::scale_fill_gradient2(low = scale_col[1],
                                                         mid = scale_col[2],
					                                               high = scale_col[3],
					                                               na.value = "transparent", 
                                                         midpoint = (scale_lim[2] - scale_lim[1])/2,
                                                         space = "Lab",
					                                               guide = "colourbar",
					                                               limits = scale_lim,
					                                               name = var_units,
					                                               oob = scales::squish)
            }

            if (add_arrows) {
              arrow_df <- data.frame(latitude = c(latitude), longitude = c(longitude), uend = c(del_u2d) + c(longitude), vend = c(del_v2d) + c(latitude))
              p <- p + ggplot2::geom_segment(arrow_df, mapping = ggplot2::aes(x = longitude, y = latitude, xend = uend, yend = vend), arrow = ggplot2::arrow(length = ggplot2::unit(0.1, "cm")))
            }
            if ((show_bathy)&!is.null(botz)) {
              # Not yet tested:
              bathy_df <- data.frame(latitude = c(latitude), longitude = c(longitude), depth = c(-botz))
              reg <- terra:rasterize(x = cbind(t(c(longitude)), t(c(latitude))), values = t(c(-botz)))
              p <- p + tidyterra::geom_spatratser_contour(reg)
            }

            if (label_towns) {
              towns <- towns[towns$latitude >= min(gy, na.rm = TRUE), ]
              towns <- towns[towns$latitude <= max(gy, na.rm = TRUE), ]
              towns <- towns[towns$longitude >= min(gx, na.rm = TRUE), ]
              towns <- towns[towns$longitude <= max(gx, na.rm = TRUE), ]
              if (dim(towns)[1]>0) p <- p + ggplot2::geom_text(data = towns, ggplot2::aes(x = longitude, y = latitude, label = town, hjust = "right"), nudge_x = -0.1) +
                                 ggplot2::geom_point(data = towns, ggplot2::aes(x = longitude, y = latitude))
            }
            #p <- p + ggplot2::ggtitle(paste(var_longname, format(chron::chron(as.double(ds[jcount])+0.000001), "%Y-%m-%d %H:%M")))
            p <- p + ggplot2::ggtitle(paste(var_longname, format(ds[jcount], format = c("d-m-yyyy", "h:m"))))
            p <- p + ggplot2::xlab("longitude") + ggplot2::ylab("latitude")
            if (!is.null(mark_points)) {
              p <- p + ggplot2::geom_point(data = mark_points, ggplot2::aes(x = longitude, y = latitude), shape = 4)
              if (plot_eta) {
                dind <- which(ds[jcount] == eta_data$date)
                p2 <- eta_plot + ggplot2::geom_point(data = eta_data[dind, ], ggplot2::aes(x = date, y = eta), size = 2, color = "red")
              }
            }
            if (gbr_poly) {
              p <- p + ggplot2::geom_path(data = sdf.gbr, ggplot2::aes(y = lat, x = long, group = group)) 
            }
            if (all(is.na(box_bounds))) { 
              p <- p + ggplot2::coord_map(xlim = c(min(gx, na.rm = TRUE), max(gx, na.rm = TRUE)), ylim = c(min(gy, na.rm = TRUE), max(gy, na.rm = TRUE)))
            } else {
              p <- p + ggplot2::coord_map(xlim = box_bounds[1:2], ylim = box_bounds[3:4]) + 
                ggplot2::theme(panel.border = ggplot2::element_rect(linetype = "solid", colour = "grey", fill = NA))
            }

            icount <- icount + 1

            if (plot_eta) {
              p <- cowplot::plot_grid(p, p2, ncol = 1, rel_heights = c(4, 1), axis = "lr", align = "v")
            }

            if (!file.exists(output_dir)) {
               dir.create(output_dir)
            }
            fname <- paste0(output_dir, "/", var_name, "_", 100000 + icount, ".png", collapse = "")
            if (verbosity>0) print(paste(var_longname, format(ds[jcount], format = c("d-m-yyyy", "h:m"))))
            ggplot2::ggsave(fname, p, dpi = 100)
            #rm("p")
         }  else {
            icount <- icount + 1
            p <- NULL
         }
         setTxtProgressBar(pb, icount/as.integer(end_date-start_date)/tstep*stride)
      } # end jcount loop
    } # end fileslist loop
  } # end month loop
  close(pb)
  if (var_name == "true_colour") {
    temporal_sum <- rgb(sqrt(temporal_sum)/icount, maxColorValue = 255)
    values <- data.frame(id = id, value = temporal_sum)
  } else {
    values <- data.frame(id = id, value = temporal_sum/icount)
  }
  datapoly <- merge(values, positions, by = c("id"))
  return(list(p = p, datapoly = datapoly, longitude = longitude, latitude = latitude))
}

#" Create a map using a dataframe in the format required by ggplot2::geom_plot, for instance from map_ereefs() or map_ereefs_movie()
#"
#" Plots a map figure in the same format as would be given by map_ereefs(), but using a pre-generated dataframe, instead of
#" processing data directly from ereefs netcdf files. Doesn't work for true_color maps.
#" 
#" @param datapoly A dataframe in the format required by geom_plot(), as provided by map_ereefs() or map_ereefs_movie().
#" @param var_longname Character vector to use for the figure title.
#" @param var_units Units to include in the figure labelling.
#" @param Land_map Set to TRUE to show a land mapof Queensland.  Default now FALSE.
#" @param scale_col Vector of colours to use for the colour scale. This can be colours 
#"      from the ggplot colour palette or a RGB hash code, or "spectral". Ignored for true_colour plots. 
#"      If set to "spectral", uses a colour spectrum from bluish to red (similar to jet but less vivid). Otherwise:
#"      If one value is given (other than "spectral"), low colour is set to ivory and high colour to the value given.
#"      If two values are given, these are used as low and high limit colours.
#"      If three values are given, the middle value is used to set the mid-point of the scale.
#"      Defaults to c("ivory", "coral4").
#" @param scale_lim Upper and lower bounds for colour scale. Defaults to full range of data.
#"      Ignored for true_colour plots.
#" @param suppress_print Set to TRUE if you don't want the plots generated and saved. Defaults to TRUE.
#" @param p Handle for an existing figure if you want to add a layer instead of creating a new figure.
#"        If p is provided, Land_map is over-ridden and set to FALSE.
#" @return p Handle for the figure generated.
#" @export
#" @examples
#" \dontrun{
#" a <- plot_map(p)
#" plot_map(a[[2]])
#"}

plot_map <- function(datapoly,
             var_longname = "",
             var_units = "",
  		       Land_map = FALSE,
             scale_col = c("ivory", "coral4"),
  		       scale_lim = c(NA, NA),
             box_bounds = c(NA, NA, NA, NA), 
             label_towns = TRUE,
             zoom = 6,
  		       p = NA,
             suppress_print = TRUE,
             gbr_poly = FALSE)
{
  if ("datapoly" %in% names(datapoly)) datapoly <- datapoly$datapoly
  if (class(datapoly$value) == "factor") {
     var_name <- "true_colour"
  } else {
     var_name <- "something else"
  }
  if ((var_name!= "true_colour")&&(is.na(scale_lim[1]))) { 
	  scale_lim <- c(min(datapoly$value, na.rm = TRUE), max(datapoly$value, na.rm = TRUE))
  }

  if (suppress_print) Land_map <- FALSE
  if (length(p)!= 1) Land_map <- FALSE
  if (Land_map) {
    p <- ggplot2::gplot() + geom_polygon(data = map.df, colour = "black", fill = "lightgrey", size = 0.5, aes(x = long, y = lat, group = group))
  } else if (length(p) == 1) {
    p <- ggplot2::ggplot()
  }
  if (!suppress_print) {
    p <- p + ggplot2::geom_polygon(ggplot2::aes(x = x, y = y, fill = value, group = id), data = datapoly)
    if (var_name == "true_colour") {
      p <- p + ggplot2::scale_fill_identity() 
    } else {
        if (scale_col[1] == "spectral") { 
          p <- p + ggplot2::scale_fill_distiller(palette = "Spectral",
                                                 na.value = "transparent", 
                                                 guide = "colourbar", 
                                                 limits = scale_lim, 
                                                 name = var_units, 
                                                 oob = scales::squish) 
        } else {
          if (length(scale_col) == 1) scale_col <- c("ivory", scale_col)
          if (length(scale_col)<3) { 
            p <- p + ggplot2::scale_fill_gradient(low = scale_col[1], 
				                                          high = scale_col[2], 
				                                          na.value = "transparent", 
				                                          guide = "colourbar",
				                                          limits = scale_lim,
				                                          name = var_units,
				                                          oob = scales::squish) 
          } else { 
            p <- p + ggplot2::scale_fill_gradient2(low = scale_col[1], 
                                                   mid = scale_col[2],
				                                           high = scale_col[3], 
				                                           na.value = "transparent", 
				                                           guide = "colourbar",
				                                           limits = scale_lim,
                                                   midpoint = (scale_lim[2] - scale_lim[1])/2,
                                                   space = "Lab",
				                                           name = var_units,
				                                           oob = scales::squish) 
          }
        }
    }
    p <- p + ggplot2::ggtitle(var_longname) + ggplot2::xlab("degrees East") + ggplot2::ylab("degrees North")
    if (label_towns) {
       towns <- data.frame(latitude = c(-15.47027987, -16.0899, -16.4840, -16.92303816, -19.26639219, -20.0136699, -20.07670986, -20.40109791, -21.15345122, -22.82406858, -23.38031858, -23.84761069, -24.8662122, -25.54073075, -26.18916037),
                    longitude = c(145.2498605, 145.4622, 145.4623, 145.7710, 146.805701, 148.2475387, 146.2635394, 148.5802016, 149.1655418, 147.6363616, 150.5059485, 151.256349, 152.3478987, 152.7049316, 152.6581893),
                    town = c("Cooktown", "Cape Tribulation", "Port Douglas", "Cairns", "Townsville", "Bowen", "Charters Towers", "Prosperine", "Mackay", "Clermont", "Rockhampton", "Gladstone", "Bundaberg", "Maryborough", "Gympie"))
       if (dim(towns)[1]>0) {
         p <- p + ggplot2::geom_text(data = towns, ggplot2::aes(x = longitude, y = latitude, label = town, hjust = "right"), nudge_x = -0.1) +
                                 ggplot2::geom_point(data = towns, ggplot2::aes(x = longitude, y = latitude))
       }
    }
    if (all(is.na(box_bounds))) { 
      p <- p + ggplot2::coord_map(xlim = c(min(datapoly$x, na.rm = TRUE), max(datapoly$x, na.rm = TRUE)), ylim = c(min(datapoly$y, na.rm = TRUE), max(datapoly$y, na.rm = TRUE)))
    } else {
      p <- p + ggplot2::coord_map(xlim = box_bounds[1:2], ylim = box_bounds[3:4]) +
        ggplot2::theme(panel.border = ggplot2::element_rect(linetype = "solid", colour = "grey", fill = NA))
    }
    if (gbr_poly) {
      p <- p + ggplot2::geom_path(data = sdf.gbr, ggplot2::aes(y = lat, x = long, group = group)) 
    }

    print(p)
  }
  return(p)
}

ereefs_surface_layer <- function(input_file, var_name) {
  dims <- ereefs_var_dims(input_file, var_name)
  k_rows <- dims[dims$role == "k", , drop = FALSE]
  if (nrow(k_rows) == 0) {
    return(NULL)
  }
  k_rows$length[[1]]
}

ereefs_extract_map_matrix <- function(input_file,
                                     var_name,
                                     target_index,
                                     layer = "surface",
                                     spatial_filters = list()) {
  if (identical(var_name, "true_color")) {
    var_name <- "true_colour"
  }

  pick_layer <- layer
  if (identical(layer, "surface")) {
    base_var <- dplyr::case_when(
      var_name == "speed" ~ if ("u1" %in% ereefs_var_names(input_file)) "u1" else "u",
      var_name == "ZooT" ~ "ZooL_N",
      TRUE ~ var_name
    )
    pick_layer <- ereefs_surface_layer(input_file, base_var)
  }

  timing <- get_origin_and_times(input_file)

  read_filters <- function(var_name_inner) {
    dims <- ereefs_var_dims(input_file, var_name_inner)
    filters <- spatial_filters
    time_row <- dims[dims$role == "time", , drop = FALSE]
    if (nrow(time_row) > 0) {
      filters[[time_row$name[[1]]]] <- ereefs_time_filter_values(
        input_file = input_file,
        time_index = target_index,
        raw_time = timing[[3]],
        time_dim_name = time_row$name[[1]],
        time_coord_dim = time_row$coord_dim[[1]]
      )
    }
    k_row <- dims[dims$role == "k", , drop = FALSE]
    if (nrow(k_row) > 0 && !is.null(pick_layer) && !is.character(pick_layer)) {
      filters[[k_row$name[[1]]]] <- as.integer(pick_layer)
    }
    filters
  }

  if (var_name == "plume") {
    spectral_vars <- c("R_412", "R_443", "R_488", "R_531", "R_547", "R_667", "R_678")
    rsr <- lapply(spectral_vars, function(x) {
      filters_x <- read_filters(x)
      dims_x <- ereefs_var_dims(input_file, x)
      ereefs_role_array(
        ereefs_read_var_array(input_file, x, filters = filters_x),
        dims = dims_x,
        target_roles = c("i", "j"),
        filters = filters_x
      )
    })
    matrix_value <- plume_class(rsr)
    return(list(values = matrix_value, units = "", long_name = "Plume colour class"))
  }

  if (var_name == "true_colour") {
    dims_470 <- ereefs_var_dims(input_file, "R_470")
    dims_555 <- ereefs_var_dims(input_file, "R_555")
    dims_645 <- ereefs_var_dims(input_file, "R_645")
    filters_470 <- read_filters("R_470")
    filters_555 <- read_filters("R_555")
    filters_645 <- read_filters("R_645")
    r470 <- ereefs_role_array(ereefs_read_var_array(input_file, "R_470", filters = filters_470), dims_470, c("i", "j"), filters = filters_470) * 10
    r555 <- ereefs_role_array(ereefs_read_var_array(input_file, "R_555", filters = filters_555), dims_555, c("i", "j"), filters = filters_555) * 10
    r645 <- ereefs_role_array(ereefs_read_var_array(input_file, "R_645", filters = filters_645), dims_645, c("i", "j"), filters = filters_645) * 10
    r470 <- ereefs_mask_array_sentinels(r470, input_file, "R_470")
    r555 <- ereefs_mask_array_sentinels(r555, input_file, "R_555")
    r645 <- ereefs_mask_array_sentinels(r645, input_file, "R_645")
    r470[r470 > 1] <- 1
    r555[r555 > 1] <- 1
    r645[r645 > 1] <- 1

    unscaled <- c(0, 30, 60, 120, 190, 255) / 255
    scaled <- c(1, 110, 160, 210, 240, 255) / 255
    scale_fun <- stats::approxfun(x = unscaled, y = scaled, yleft = 1, yright = 255)
    red <- scale_fun(r645)
    green <- scale_fun(r555)
    blue <- scale_fun(r470)
    red[is.na(red)] <- 0
    green[is.na(green)] <- 0
    blue[is.na(blue)] <- 0
    colour_value <- rgb(red, green, blue)
    colour_value[colour_value == "#000000"] <- NA
    return(list(values = array(as.character(colour_value), dim = dim(r645)), units = "", long_name = "Simulated true colour"))
  }

  if (var_name == "ZooT") {
    dims_z1 <- ereefs_var_dims(input_file, "ZooL_N")
    dims_z2 <- ereefs_var_dims(input_file, "ZooS_N")
    filters_z1 <- read_filters("ZooL_N")
    filters_z2 <- read_filters("ZooS_N")
    z1 <- ereefs_role_array(ereefs_read_var_array(input_file, "ZooL_N", filters = filters_z1), dims_z1, c("i", "j"), filters = filters_z1)
    z2 <- ereefs_role_array(ereefs_read_var_array(input_file, "ZooS_N", filters = filters_z2), dims_z2, c("i", "j"), filters = filters_z2)
    z1 <- ereefs_mask_array_sentinels(z1, input_file, "ZooL_N")
    z2 <- ereefs_mask_array_sentinels(z2, input_file, "ZooS_N")
    return(list(
      values = z1 + z2,
      units = ereefs_var_attr(input_file, "ZooL_N", "units"),
      long_name = "Total Zooplankton"
    ))
  }

  if (var_name == "speed") {
    u_name <- if ("u1" %in% ereefs_var_names(input_file)) "u1" else "u"
    v_name <- if ("u2" %in% ereefs_var_names(input_file)) "u2" else "v"
    dims_u <- ereefs_var_dims(input_file, u_name)
    dims_v <- ereefs_var_dims(input_file, v_name)
    filters_u <- read_filters(u_name)
    filters_v <- read_filters(v_name)
    u <- ereefs_role_array(ereefs_read_var_array(input_file, u_name, filters = filters_u), dims_u, c("i", "j"), filters = filters_u)
    v <- ereefs_role_array(ereefs_read_var_array(input_file, v_name, filters = filters_v), dims_v, c("i", "j"), filters = filters_v)
    u <- ereefs_mask_array_sentinels(u, input_file, u_name)
    v <- ereefs_mask_array_sentinels(v, input_file, v_name)
    return(list(
      values = sqrt(u^2 + v^2),
      units = ereefs_var_attr(input_file, u_name, "units"),
      long_name = "Current speed"
    ))
  }

  dims_out <- ereefs_var_dims(input_file, var_name)
  filters_out <- read_filters(var_name)
  list(
    values = ereefs_role_array(
      ereefs_read_var_array(input_file, var_name, filters = filters_out),
      dims = dims_out,
      target_roles = c("i", "j"),
      filters = filters_out
    ) %>% ereefs_mask_array_sentinels(input_file = input_file, var_name = var_name),
    units = ereefs_var_attr(input_file, var_name, "units"),
    long_name = ereefs_var_attr(input_file, var_name, "long_name")
  )
}

ereefs_build_map_geometry <- function(grids, box_bounds = c(NA, NA, NA, NA)) {
  coord_table <- grids$spatial_grid %>%
    dplyr::arrange(i, j)

  a <- nrow(grids$x_grid) - 1
  b <- ncol(grids$y_grid) - 1
  lookup <- expand.grid(i = seq_len(a), j = seq_len(b))
  lookup$id <- seq_len(nrow(lookup))
  cells <- dplyr::left_join(dplyr::as_tibble(lookup), coord_table, by = c("i", "j"))

  gx <- c(
    grids$x_grid[1:a, 1:b],
    grids$x_grid[2:(a + 1), 1:b],
    grids$x_grid[2:(a + 1), 2:(b + 1)],
    grids$x_grid[1:a, 2:(b + 1)]
  )
  gy <- c(
    grids$y_grid[1:a, 1:b],
    grids$y_grid[2:(a + 1), 1:b],
    grids$y_grid[2:(a + 1), 2:(b + 1)],
    grids$y_grid[1:a, 2:(b + 1)]
  )
  gx <- matrix(gx, nrow = a * b, ncol = 4)
  gy <- matrix(gy, nrow = a * b, ncol = 4)
  valid <- stats::complete.cases(gx) & stats::complete.cases(gy)
  cells$valid_polygon <- valid

  if (!all(is.na(box_bounds))) {
    cells <- cells %>%
      dplyr::filter(
        dplyr::between(longitude, box_bounds[1], box_bounds[2]),
        dplyr::between(latitude, box_bounds[3], box_bounds[4])
      )
  }

  keep_ids <- cells$id[cells$valid_polygon]
  pos_x <- as.vector(t(gx[keep_ids, , drop = FALSE]))
  pos_y <- as.vector(t(gy[keep_ids, , drop = FALSE]))
  positions <- dplyr::tibble(
    id = rep(keep_ids, each = 4),
    x = pos_x,
    y = pos_y
  )

  list(cells = cells, positions = positions)
}

ereefs_plot_datapoly <- function(datapoly,
                                 var_name,
                                 var_longname = "",
                                 var_units = "",
                                 Land_map = FALSE,
                                 scale_col = "viridis",
                                 scale_lim = c(NA, NA),
                                 plot_style = c("polygon", "smooth"),
                                 smooth_pixels = 600,
                                 box_bounds = c(NA, NA, NA, NA),
                                 label_towns = TRUE,
                                 p = NA,
                                 suppress_print = TRUE,
                                 gbr_poly = FALSE,
                                 mark_points = NULL) {
  plot_style <- match.arg(plot_style)
  if (!inherits(p, "ggplot")) {
    p <- ggplot2::ggplot()
  }
  if (Land_map && exists("map.df")) {
    p <- p + ggplot2::geom_polygon(
      data = map.df,
      colour = "black",
      fill = "lightgrey",
      linewidth = 0.3,
      ggplot2::aes(x = long, y = lat, group = group)
    )
  }

  if (identical(plot_style, "smooth")) {
    bbox <- if (all(is.na(box_bounds))) {
      c(
        min(datapoly$x, na.rm = TRUE),
        max(datapoly$x, na.rm = TRUE),
        min(datapoly$y, na.rm = TRUE),
        max(datapoly$y, na.rm = TRUE)
      )
    } else {
      box_bounds
    }
    width <- bbox[2] - bbox[1]
    height <- bbox[4] - bbox[3]
    if (!is.finite(width) || !is.finite(height) || width <= 0 || height <= 0) {
      plot_style <- "polygon"
    } else {
      ncols <- max(50L, as.integer(smooth_pixels))
      nrows <- max(50L, as.integer(round(smooth_pixels * height / width)))
      sv <- poly2sv(datapoly)
      r <- terra::rast(
        xmin = bbox[1],
        xmax = bbox[2],
        ymin = bbox[3],
        ymax = bbox[4],
        ncols = ncols,
        nrows = nrows
      )
      r <- terra::rasterize(sv, r, field = "value")
      raster_df <- terra::as.data.frame(r, xy = TRUE, na.rm = FALSE)
      names(raster_df)[3] <- "value"
      p <- p + ggplot2::geom_raster(
        data = raster_df,
        ggplot2::aes(x = x, y = y, fill = value),
        interpolate = TRUE
      )
    }
  }

  if (identical(plot_style, "polygon")) {
    p <- p + ggplot2::geom_polygon(
      data = datapoly,
      colour = NA,
      linewidth = 0,
      ggplot2::aes(x = x, y = y, fill = value, group = id)
    )
  }

  if (var_name == "true_colour") {
    p <- p + ggplot2::scale_fill_identity()
  } else if (scale_col[1] == "spectral") {
    p <- p + ggplot2::scale_fill_distiller(
      palette = "Spectral",
      na.value = "transparent",
      guide = "colourbar",
      limits = scale_lim,
      name = var_units,
      oob = scales::squish
    )
  } else if (scale_col[1] %in% c("viridis", "magma", "inferno", "plasma", "cividis", "turbo")) {
    palette_name <- switch(
      scale_col[1],
      viridis = "viridis",
      magma = "magma",
      inferno = "inferno",
      plasma = "plasma",
      cividis = "cividis",
      turbo = "turbo"
    )
    p <- p + ggplot2::scale_fill_gradientn(
      colours = grDevices::hcl.colors(256, palette_name),
      na.value = "transparent",
      limits = scale_lim,
      name = var_units,
      oob = scales::squish
    )
  } else if (length(scale_col) < 3) {
    if (length(scale_col) == 1) scale_col <- c("ivory", scale_col)
    p <- p + ggplot2::scale_fill_gradient(
      low = scale_col[1],
      high = scale_col[2],
      na.value = "transparent",
      guide = "colourbar",
      limits = scale_lim,
      name = var_units,
      oob = scales::squish
    )
  } else {
    p <- p + ggplot2::scale_fill_gradient2(
      low = scale_col[1],
      mid = scale_col[2],
      high = scale_col[3],
      midpoint = mean(scale_lim, na.rm = TRUE),
      na.value = "transparent",
      guide = "colourbar",
      limits = scale_lim,
      name = var_units,
      oob = scales::squish
    )
  }

  if (label_towns) {
    towns <- ereefs_towns()
    if (!all(is.na(box_bounds))) {
      towns <- towns %>%
        dplyr::filter(
          dplyr::between(longitude, box_bounds[1], box_bounds[2]),
          dplyr::between(latitude, box_bounds[3], box_bounds[4])
        )
    }
    if (nrow(towns) > 0) {
      p <- p +
        ggplot2::geom_text(data = towns, ggplot2::aes(x = longitude, y = latitude, label = town), hjust = 1.1) +
        ggplot2::geom_point(data = towns, ggplot2::aes(x = longitude, y = latitude), size = 0.8)
    }
  }

  if (!is.null(mark_points)) {
    if (is.null(dim(mark_points))) {
      mark_points <- data.frame(latitude = mark_points[1], longitude = mark_points[2])
    }
    p <- p + ggplot2::geom_point(data = mark_points, ggplot2::aes(x = longitude, y = latitude), shape = 4)
  }

  if (gbr_poly && exists("sdf.gbr")) {
    p <- p + ggplot2::geom_path(data = sdf.gbr, ggplot2::aes(y = lat, x = long, group = group))
  }

  if (all(is.na(box_bounds))) {
    xrange <- range(datapoly$x, na.rm = TRUE)
    yrange <- range(datapoly$y, na.rm = TRUE)
    xpad <- max(0.02 * diff(xrange), 0.05)
    ypad <- max(0.02 * diff(yrange), 0.05)
    p <- p + ggplot2::coord_quickmap(
      xlim = c(xrange[[1]] - xpad, xrange[[2]] + xpad),
      ylim = c(yrange[[1]] - ypad, yrange[[2]] + ypad)
    )
  } else {
    p <- p + ggplot2::coord_quickmap(xlim = box_bounds[1:2], ylim = box_bounds[3:4])
  }

  p +
    ggplot2::ggtitle(var_longname) +
    ggplot2::xlab("longitude") +
    ggplot2::ylab("latitude")
}

#' Create a map of eReefs or other EMS model output
#'
#' Extracts a single time slice from a local NetCDF file, OPeNDAP dataset, or
#' THREDDS catalog-backed workflow and returns a `ggplot2` map. The function
#' supports both curvilinear EMS grids with cell corners and regular regridded
#' products that provide cell centres only.
#'
#' @param var_name Variable to plot. Special values `"true_colour"` and
#'   `"plume"` are also supported.
#' @param target_date Date or date-time to plot. If an exact match is not
#'   available, the nearest model output time is used.
#' @param layer Layer selector. Use a positive layer index, a negative depth
#'   below mean sea level, or `"surface"`/`"bottom"`.
#' @param Land_map Logical; add a simple land underlay where available.
#' @param input_file NetCDF file path, OPeNDAP URL, or THREDDS catalog URL.
#' @param input_grid Optional alternative source for grid metadata.
#' @param scale_col Colour scale specification. Defaults to `"viridis"` for
#'   scalar maps.
#' @param scale_lim Numeric colour limits. If left as `NA`, limits are inferred
#'   from the extracted data.
#' @param plot_style Either `"polygon"` or `"smooth"`.
#' @param smooth_pixels Resolution used for `"smooth"` display maps.
#' @param zoom Deprecated legacy argument retained for backward compatibility.
#' @param box_bounds Optional bounds in the form
#'   `c(longitude_min, longitude_max, latitude_min, latitude_max)`.
#' @param p Optional existing plot object to add to.
#' @param suppress_print Logical; if `TRUE`, suppresses automatic printing.
#' @param return_poly Logical; if `TRUE`, return a list containing the plot and
#'   plotting polygons instead of only the plot object.
#' @param label_towns Logical; add town labels where available.
#' @param strict_bounds Deprecated legacy argument retained for backward
#'   compatibility.
#' @param mark_points Optional locations to mark on the map.
#' @param gbr_poly Logical; add GBR polygon overlay when available.
#'
#' @return A `ggplot2` object, or if `return_poly = TRUE`, a list containing the
#'   plot, plotting polygons, cell values, and associated metadata.
#' @export
map_ereefs <- function(var_name = "true_colour",
                       target_date = c(2018, 1, 30),
                       layer = "surface",
                       Land_map = FALSE,
                       input_file = "catalog",
                       input_grid = NA,
                       scale_col = "viridis",
                       scale_lim = c(NA, NA),
                       plot_style = c("polygon", "smooth"),
                       smooth_pixels = 600,
                       zoom = 6,
                       box_bounds = c(NA, NA, NA, NA),
                       p = NA,
                       suppress_print = TRUE,
                       return_poly = FALSE,
                       label_towns = TRUE,
                       strict_bounds = FALSE,
                       mark_points = NULL,
                       gbr_poly = FALSE) {
  plot_style <- match.arg(plot_style)
  rm(zoom, strict_bounds)
  assignList(get_params(target_date, target_date, input_file, var_name))

  map_date <- as.Date(start_date)
  file_to_use <- ereefs_pick_file_for_date(resolved_files, map_date)
  timing <- get_origin_and_times(file_to_use)
  day_index <- which.min(abs(timing[[2]] - start_date))
  if (min(abs(timing[[2]] - start_date)) > lubridate::days(1)) {
    warning("Target date is more than one day from the closest available output in the selected file.")
  }

  grids <- get_ereefs_grids(file_to_use, input_grid)
  subset_spec <- ereefs_spatial_subset_spec(
    input_file = file_to_use,
    var_name = ereefs_map_reference_var(file_to_use, var_name),
    spatial_grid = grids$spatial_grid,
    box_bounds = box_bounds
  )
  extracted <- ereefs_extract_map_matrix(
    file_to_use,
    var_name,
    day_index,
    layer = layer,
    spatial_filters = subset_spec$filters
  )
  geometry <- ereefs_build_map_geometry(grids, box_bounds = box_bounds)
  cells <- geometry$cells
  values <- as.vector(extracted$values)
  expected_cells <- length(subset_spec$i_vals) * length(subset_spec$j_vals)
  if (length(values) != expected_cells) {
    stop("Mapped values do not align with the requested horizontal subset.")
  }

  cell_values <- ereefs_array_index_grid(subset_spec$i_vals, subset_spec$j_vals) %>%
    dplyr::mutate(value = values) %>%
    dplyr::left_join(grids$spatial_grid %>% dplyr::select(i, j, longitude, latitude), by = c("i", "j")) %>%
    dplyr::inner_join(cells %>% dplyr::select(id, i, j, valid_polygon), by = c("i", "j")) %>%
    dplyr::filter(valid_polygon)

  datapoly <- dplyr::left_join(cell_values %>% dplyr::select(id, value), geometry$positions, by = "id")

  if (var_name != "true_colour" && any(is.na(scale_lim))) {
    scale_lim <- c(min(cell_values$value, na.rm = TRUE), max(cell_values$value, na.rm = TRUE))
  }

  plot_title <- paste(extracted$long_name, format(timing[[2]][day_index], "%Y-%m-%d %H:%M"))
  p <- ereefs_plot_datapoly(
    datapoly = datapoly,
    var_name = if (var_name == "true_color") "true_colour" else var_name,
    var_longname = plot_title,
    var_units = extracted$units,
    Land_map = Land_map,
    scale_col = scale_col,
    scale_lim = scale_lim,
    plot_style = plot_style,
    smooth_pixels = smooth_pixels,
    box_bounds = box_bounds,
    label_towns = label_towns,
    p = p,
    suppress_print = suppress_print,
    gbr_poly = gbr_poly,
    mark_points = mark_points
  )

  result <- list(
    p = p,
    datapoly = datapoly,
    longitude = cell_values$longitude,
    latitude = cell_values$latitude,
    cell_values = cell_values,
    var_longname = extracted$long_name,
    var_units = extracted$units,
    date = timing[[2]][day_index]
  )
  if (return_poly) {
    return(result)
  }
  p
}

#' Create a sequence of eReefs maps and optionally assemble an animation
#'
#' Repeatedly calls [map_ereefs()] over a date range, saving frame images when
#' requested and optionally assembling them into a GIF or MP4 animation.
#'
#' @inheritParams map_ereefs
#' @param start_date Start of the animation period.
#' @param end_date End of the animation period.
#' @param output_dir Directory in which to save frame images and animation
#'   outputs.
#' @param save_frames Logical; save individual frame PNG files.
#' @param animation_format One of `"none"`, `"gif"`, `"mp4"`, or legacy
#'   `"mp3"` (which is treated as `"mp4"` with a warning).
#' @param animation_file Optional output filename for the assembled animation.
#' @param fps Frames per second for assembled animations.
#' @param stride Temporal stride used to choose frames.
#' @param verbosity Verbosity level for progress messages.
#' @param add_arrows Deprecated legacy argument retained for backward
#'   compatibility.
#' @param max_u Deprecated legacy argument retained for backward compatibility.
#' @param scale_arrows Deprecated legacy argument retained for backward
#'   compatibility.
#' @param show_bathy Deprecated legacy argument retained for backward
#'   compatibility.
#' @param contour_breaks Deprecated legacy argument retained for backward
#'   compatibility.
#'
#' @return A list containing the final plot, averaged plotting polygons, saved
#'   frame filenames, animation filename, and the colour scale limits used
#'   across frames.
#' @export
map_ereefs_movie <- function(var_name = "true_colour",
                             start_date = c(2015, 12, 1),
                             end_date = c(2016, 3, 31),
                             layer = "surface",
                             output_dir = "ToAnimate",
                             save_frames = TRUE,
                             animation_format = c("none", "gif", "mp4", "mp3"),
                             animation_file = NA,
                             fps = 2,
                             Land_map = FALSE,
                             input_file = "catalog",
                             input_grid = NA,
                             scale_col = "viridis",
                             scale_lim = c(NA, NA),
                             plot_style = c("polygon", "smooth"),
                             smooth_pixels = 600,
                             zoom = 6,
                             box_bounds = c(NA, NA, NA, NA),
                             suppress_print = TRUE,
                             stride = "daily",
                             verbosity = 0,
                             label_towns = TRUE,
                             strict_bounds = FALSE,
                             mark_points = NULL,
                             gbr_poly = FALSE,
                             add_arrows = FALSE,
                             max_u = NA,
                             scale_arrows = NA,
                             show_bathy = FALSE,
                             contour_breaks = c(5, 10, 20)) {
  plot_style <- match.arg(plot_style)
  animation_format <- match.arg(animation_format)
  if (identical(animation_format, "mp3")) {
    warning("animation_format = 'mp3' is not valid for image sequences. Using 'mp4' instead.")
    animation_format <- "mp4"
  }
  rm(zoom, strict_bounds, add_arrows, max_u, scale_arrows, show_bathy, contour_breaks)
  start_date <- get_date_time(start_date)
  end_date <- get_date_time(end_date)
  date_sequence <- ereefs_dates_for_stride(start_date, end_date, stride)
  file_table <- ereefs_resolve_time_files(input_file, start_date, end_date)

  need_output_dir <- isTRUE(save_frames) || !identical(animation_format, "none")
  if (need_output_dir && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  frame_files <- character(0)

  if (is.na(animation_file)) {
    animation_file <- if (identical(animation_format, "gif")) {
      file.path(output_dir, sprintf("%s_animation.gif", var_name))
    } else if (identical(animation_format, "mp4")) {
      file.path(output_dir, sprintf("%s_animation.mp4", var_name))
    } else {
      NA_character_
    }
  }

  accumulated <- NULL
  last_frame <- NULL
  stored_frames <- vector("list", length(date_sequence))
  effective_scale_lim <- scale_lim
  for (frame_index in seq_along(date_sequence)) {
    if (verbosity > 0) {
      message("Loading frame ", frame_index, " of ", length(date_sequence), " for ", date_sequence[[frame_index]])
    }

    frame <- map_ereefs(
      var_name = var_name,
      target_date = as.Date(date_sequence[[frame_index]]),
      layer = layer,
      Land_map = Land_map,
      input_file = ereefs_pick_file_for_date(file_table, date_sequence[[frame_index]]),
      input_grid = input_grid,
      scale_col = scale_col,
      scale_lim = scale_lim,
      plot_style = plot_style,
      smooth_pixels = smooth_pixels,
      box_bounds = box_bounds,
      suppress_print = TRUE,
      return_poly = TRUE,
      label_towns = label_towns,
      mark_points = mark_points,
      gbr_poly = gbr_poly
    )
    stored_frames[[frame_index]] <- frame

    frame_values <- frame$cell_values %>% dplyr::select(id, value)
    if (is.null(accumulated)) {
      accumulated <- frame_values
    } else if (var_name == "true_colour") {
      accumulated$value <- mapply(
        function(a, b) {
          if (is.na(a)) return(b)
          if (is.na(b)) return(a)
          rgb((col2rgb(a) + col2rgb(b)) / 2, maxColorValue = 255)
        },
        accumulated$value,
        frame_values$value
      )
    } else {
      accumulated$value <- accumulated$value + frame_values$value
    }
  }

  if (var_name != "true_colour" && any(is.na(effective_scale_lim))) {
    frame_ranges <- lapply(stored_frames, function(frame) {
      vals <- frame$cell_values$value
      vals <- vals[is.finite(vals)]
      if (!length(vals)) {
        return(c(NA_real_, NA_real_))
      }
      c(min(vals), max(vals))
    })
    frame_mins <- vapply(frame_ranges, `[[`, numeric(1), 1)
    frame_maxs <- vapply(frame_ranges, `[[`, numeric(1), 2)
    if (is.na(effective_scale_lim[[1]])) {
      effective_scale_lim[[1]] <- min(frame_mins, na.rm = TRUE)
    }
    if (is.na(effective_scale_lim[[2]])) {
      effective_scale_lim[[2]] <- max(frame_maxs, na.rm = TRUE)
    }
  }

  for (frame_index in seq_along(stored_frames)) {
    frame <- stored_frames[[frame_index]]
    if (verbosity > 0) {
      message("Rendering frame ", frame_index, " of ", length(stored_frames), " for ", frame$date)
    }
    frame_values <- frame$cell_values %>% dplyr::select(id, value)
    frame_datapoly <- dplyr::left_join(frame_values, frame$datapoly %>% dplyr::select(id, x, y), by = "id")
    plot_title <- paste(frame$var_longname, format(frame$date, "%Y-%m-%d %H:%M"))
    p <- ereefs_plot_datapoly(
      datapoly = frame_datapoly,
      var_name = if (var_name == "true_color") "true_colour" else var_name,
      var_longname = plot_title,
      var_units = frame$var_units,
      Land_map = Land_map,
      scale_col = scale_col,
      scale_lim = effective_scale_lim,
      plot_style = plot_style,
      smooth_pixels = smooth_pixels,
      box_bounds = box_bounds,
      label_towns = label_towns,
      suppress_print = suppress_print,
      gbr_poly = gbr_poly,
      mark_points = mark_points
    )
    if (isTRUE(save_frames)) {
      fname <- file.path(output_dir, sprintf("%s_%05d.png", var_name, frame_index))
      ggplot2::ggsave(fname, p, dpi = 100)
      frame_files <- c(frame_files, fname)
    }
    last_frame <- p
  }

  if (!is.null(accumulated) && var_name != "true_colour") {
    accumulated$value <- accumulated$value / length(date_sequence)
  }
  frame <- stored_frames[[length(stored_frames)]]
  positions <- frame$datapoly %>% dplyr::select(id, x, y)
  datapoly <- dplyr::left_join(accumulated, positions, by = "id")
  if (!identical(animation_format, "none")) {
    if (!length(frame_files)) {
      stop("animation_format requested but no frame files were saved.")
    }
    if (!requireNamespace("magick", quietly = TRUE)) {
      stop("animation_format requested but the 'magick' package is not installed.")
    }
    animation <- magick::image_read(frame_files)
    if (identical(animation_format, "gif")) {
      animation <- magick::image_animate(animation, fps = fps)
      magick::image_write(animation, path = animation_file)
    } else if (identical(animation_format, "mp4")) {
      magick::image_write_video(animation, path = animation_file, fps = fps)
    }
  }
  list(
    p = last_frame,
    datapoly = datapoly,
    longitude = frame$longitude,
    latitude = frame$latitude,
    frame_files = frame_files,
    animation_file = animation_file,
    output_dir = output_dir,
    scale_lim = effective_scale_lim
  )
}

#' Plot an already extracted eReefs polygon map object
#'
#' A lower-level plotting helper that turns a polygon/value table, or the list
#' returned by `map_ereefs(..., return_poly = TRUE)`, into a `ggplot2` map.
#'
#' @param datapoly Polygon/value table in the format used internally by the
#'   package, or a list returned by `map_ereefs(..., return_poly = TRUE)`.
#' @param var_longname Optional long name for the plot title.
#' @param var_units Optional units label.
#' @param Land_map Logical; add a simple land underlay where available.
#' @param scale_col Colour scale specification. Defaults to `"viridis"` for
#'   scalar maps.
#' @param scale_lim Numeric colour limits. If left as `NA`, limits are inferred
#'   from the data.
#' @param plot_style Either `"polygon"` or `"smooth"`.
#' @param smooth_pixels Resolution used for `"smooth"` display maps.
#' @param box_bounds Optional bounds in the form
#'   `c(longitude_min, longitude_max, latitude_min, latitude_max)`.
#' @param label_towns Logical; add town labels where available.
#' @param zoom Deprecated legacy argument retained for backward compatibility.
#' @param p Optional existing plot object to add to.
#' @param suppress_print Logical; if `TRUE`, suppresses automatic printing.
#' @param gbr_poly Logical; add GBR polygon overlay when available.
#'
#' @return A `ggplot2` object.
#' @export
plot_map <- function(datapoly,
                     var_longname = "",
                     var_units = "",
                     Land_map = FALSE,
                     scale_col = "viridis",
                     scale_lim = c(NA, NA),
                     plot_style = c("polygon", "smooth"),
                     smooth_pixels = 600,
                     box_bounds = c(NA, NA, NA, NA),
                     label_towns = TRUE,
                     zoom = 6,
                     p = NA,
                     suppress_print = TRUE,
                     gbr_poly = FALSE) {
  plot_style <- match.arg(plot_style)
  rm(zoom)
  if ("datapoly" %in% names(datapoly)) {
    var_longname <- if (nzchar(var_longname)) var_longname else datapoly$var_longname
    var_units <- if (nzchar(var_units)) var_units else datapoly$var_units
    datapoly <- datapoly$datapoly
  }
  var_name <- if (is.character(datapoly$value)) "true_colour" else "standard"
  if (var_name != "true_colour" && any(is.na(scale_lim))) {
    scale_lim <- c(min(datapoly$value, na.rm = TRUE), max(datapoly$value, na.rm = TRUE))
  }
  ereefs_plot_datapoly(
    datapoly = datapoly,
    var_name = var_name,
    var_longname = var_longname,
    var_units = var_units,
    Land_map = Land_map,
    scale_col = scale_col,
    scale_lim = scale_lim,
    plot_style = plot_style,
    smooth_pixels = smooth_pixels,
    box_bounds = box_bounds,
    label_towns = label_towns,
    p = p,
    suppress_print = suppress_print,
    gbr_poly = gbr_poly
  )
}

# metadata:
# - gpt_version: GPT-5 Codex
# - time: 16:03
# - date: 2026-04-26
# - prompt_used: "Refactor the eReefs R toolkit away from ncdf4 toward tidync, add regular-grid support alongside curvilinear grids, archive existing R files, and refresh docs/examples."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 16:21
# - date: 2026-04-26
# - prompt_used: "Refactor the eReefs R toolkit away from ncdf4 toward tidync, add regular-grid support alongside curvilinear grids, archive existing R files, and refresh docs/examples."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 15:33
# - date: 2026-04-28
# - prompt_used: "Make sure the roxygen comments and generated help pages are up to date and match the active refactored code."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 19:10
# - date: 2026-04-26
# - prompt_used: "Install dependencies, verify the refactored toolkit, improve efficiency for large THREDDS-served files, and build a working Jupyter demo notebook."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 21:29
# - date: 2026-04-26
# - prompt_used: "Continue the live OPeNDAP refactor by removing stem-based file reconstruction, resolving time ranges from THREDDS catalogs, and supporting centre-only NCI simple files."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 00:04
# - date: 2026-04-27
# - prompt_used: "Finish the tidy tidync-first refactor, keep it efficient for large live OPeNDAP datasets, validate depth/free-surface logic, audit dependencies, and review plotting palettes."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 00:12
# - date: 2026-04-27
# - prompt_used: "Add an optional smoother map rendering mode so mapped cells do not show visible grid seams in vignette figures."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 00:31
# - date: 2026-04-27
# - prompt_used: "Fix the smooth map plotting return path and add a vignette example for the opt-in smoother display mode."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 00:39
# - date: 2026-04-27
# - prompt_used: "Fix the map plotting helper so opt-in smooth maps return a real ggplot object and can be documented safely."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 01:20
# - date: 2026-04-27
# - prompt_used: "Fix local EMS map extraction so time filters use real coordinate values instead of one-based indices."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 01:28
# - date: 2026-04-27
# - prompt_used: "Fix local EMS time filtering so coordinate-based time dimensions use raw values while record-style time dimensions use positional indices."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 14:24
# - date: 2026-04-27
# - prompt_used: "Fix map/grid orientation bugs by reordering arrays explicitly by i/j/k/time roles before flattening them onto polygons."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 14:30
# - date: 2026-04-27
# - prompt_used: "Make the role-based array reordering filter-aware so singleton time and k dimensions dropped by the backend do not break map orientation fixes."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 14:39
# - date: 2026-04-27
# - prompt_used: "Mask NetCDF-style sentinel values like 1e35 before plotting so notebook maps do not show false bright edge bands."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 16:58
# - date: 2026-04-27
# - prompt_used: "Make notebook plots save as well as display, and have map_ereefs_movie save frames then build an animation file."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 18:15
# - date: 2026-04-27
# - prompt_used: "Fix the live two-month time-series example, keep movie colour scales fixed across frames, and suppress repeated z_grid reconstruction warnings."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 08:41
# - date: 2026-04-28
# - prompt_used: "Update vignette PNGs to match the current code and make viridis the default colour palette."
# metadata:
# - gpt_version: GPT-5 Codex
# - time: 09:22
# - date: 2026-04-28
# - prompt_used: "Fix the remaining howto vignette map clipping, add the missing smooth-map figure, and switch the profile and slice examples to live OPeNDAP data."
