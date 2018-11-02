#' Determines whether the file provided looks like an example of daily (e.g. GBR1), monthly (e.g. GBR4) or other (e.g. RECOM) output
#'
#' If there are two dashes ('-') in the filename, it is assumed to be a daily output file. If there is one, it is assumed to
#' be monthly. If there are no dashes, it's something else, such as a RECOM file.
#'
#' @param filename Name of the file to examine
#' @return 1 for daily, 4 for monthly, 0 for other
#' @export
get_ereefs_case <- function(filename) {
  ext <- stringi::stri_locate_last(filename, regex='.nc')[1]
  lastfew <- substr(filename, start=ext-10, stop=ext-1)
  dashcount <- stringi::stri_count_fixed(lastfew, '-')
  if (dashcount==2) {
	  ereefs_case <- 1
  } else if (dashcount==1) {
	  ereefs_case <- 4
  } else {
	  ereefs_case <- 0
  }
  return(ereefs_case)
}

#' Returns an appropriate x_grid, y_grid and z_grid
#'
#' @param filename Name of the file to examine
#' @param input_grid Optional alternative file from which to load x_grid, y_grid and z_grid
#' @return a list containing x_grid, y_grid and z_grid
#' @export
get_ereefs_grids <- function(filename, input_grid=NA) {
	if (!is.na(input_grid)) {
		nc <- ncdf4::nc_open(input_grid)
		x_grid <- ncdf4::ncvar_get(nc, 'x_grid')
		y_grid <- ncdf4::ncvar_get(nc, 'y_grid')
		z_grid <- ncdf4::ncvar_get(nc, 'z_grid')
		ncdf4::nc_close(nc)
	} else {
		nc <- ncdf4::nc_open(filename)
		if (!is.null(nc$var[['x_grid']])) { 
		  x_grid <- ncdf4::ncvar_get(nc, 'x_grid')
		  y_grid <- ncdf4::ncvar_get(nc, 'y_grid')
		  z_grid <- ncdf4::ncvar_get(nc, 'z_grid')
                } else {
		  if (!is.null(nc$dim[['i']])) {
		    ilen <- nc$dim[['i']][['len']]
		    jlen <- nc$dim[['j']][['len']]
		    klen <- nc$dim[['k']][['len']]
                  } else {
		    ilen <- nc$dim[['i_grid']][['len']]
		    jlen <- nc$dim[['j_grid']][['len']]
		    klen <- nc$dim[['k_grid']][['len']]
                  }
		  if ((ilen==510)&(jlen==2389)) {
		    x_grid <- gbr1_x_grid
		    y_grid <- gbr1_y_grid
		    z_grid <- gbr1_z_grid
		  } else if ((ilen==600)&(jlen==180)) {
		    x_grid <- gbr4_x_grid
		    y_grid <- gbr4_y_grid
		    z_grid <- gbr4_z_grid
		  } else {
		    stop('x_grid, y_grid and z_grid not found and grid size not recognised as GBR1 or GBR4. Please specify input_grid.')
		  }
                }
		ncdf4::nc_close(nc)
	} 
	return(list(x_grid=x_grid, y_grid=y_grid, z_grid=z_grid))
}

#' Returns the 'stem' of a filename (i.e. the filename with the date and extension cut off the end)
#'
#' Utility function for ereefs package
#'
#' @param filename A netcdf filename from an EMS model run
#' @return a character string
#' @export
get_file_stem <- function(filename) {
	case <- get_ereefs_case(filename)
	if (case==0) {
		file_stem <- substr(filename, start=1, stop=stringi::stri_locate_last(filename, regex='.nc')[1]-1)
	} else if (case==1) {
		file_stem <- substr(filename, start=1, stop=stringi::stri_locate_last(filename, regex='.nc')[1]-11)
	} else {
		file_stem <- substr(filename, start=1, stop=stringi::stri_locate_last(filename, regex='.nc')[1]-8)
	}
	return(file_stem)
}

#' Check whether the platform is Windows and the filename contains "http": if so, return a warning
#'
#' Utility function for the eReefs package
#' @param input_stem A netcdf filename or file stem
#' @return TRUE or stop error
#' @export
check_platform_ok <- function(input_stem)
{
  ok <- TRUE
  webserved <- stringi::stri_startswith_fixed(input_stem, "http:")
  if ((.Platform$OS.type=="windows")&&(webserved)) {
     ok <- FALSE
     warning("Unfortunately, under Windows this function will only work with locally-stored netcdf files, not web-served files UNLESS you use a specially compiled version of the ncdf4 package -- contact b.robson@aims.gov.au for details.")
  }
  return(ok)
}

#' Check whether the given filename is a shortcut and if so, set the full filename
#'
#' Utility function for the eReefs package
#' @param input_file A netcdf filename or file stem
#' @return input_file
#' @export
substitute_filename <- function(input_file) {
  choices  <- c("GBR4-v2.0",
                "GBR4_BGC-v2.0 Chyb Dcrt",
                "GBR4_BGC-v2.0 Chyb Dnrt",
                "GBR4_BGC-v2.0 Cpre Dcrt",
                "GBR4 rivers_2.0  Dnrt   GBR4 passive",
                "GBR1-v2.0",
                "GBR1 rivers_2.0  Dnrt",
                "GBR4-v1.85",
                "GBR4_BGC-v926",
                "GBR4_BGC-v924",
                "GBR1-v1.71",
                "GBR1_BGC-v924",
                "menu")
  if (is.numeric(input_file)) {
     input_file <- choices[input_file]
  } else if (input_file == "choices") {
     print(choices)
     stop()
  } else if ((input_file == "menu")||(input_file == length(choices))) {
     selection <- utils::menu(c("Latest release 4km grid hydrodynamic model (Sept 2010-pres.)", 
                                "Latest release 4km biogeochemical model hindcast (Sept 2010 - Oct 2016)",
                                "Latest release 4km biogeochemical model near real time (Oct 2016 - pres.)",
                                "Pre-industrial catchment scanerio 4km BGC (Sept 2010 - Oct 2016)",
                                "Latest release passive river tracers (Sept 2010 - pres.)",
                                "Latest release 1km grid hydrodynamic model (Dec 2014 - pres.)",
                                "Latest release 1km grid passive river tracers (Dec 2014 - pres.)",
                                "Archived 4km hydro (v 1.85, Sept 2010-pres.)",
                                "Archived 4km bgc (v926, Sept 2010 - Dec 2014)",
                                "Archived 4km bgc (v924, Sept 2010 - Sept 2017)",
                                "Archived 1km hydro (v 1.71, Dec 2014 – Apr 2016)",
                                "Archived 1km bgc (v924, Dec 2014 – pres.)"))
     input_file <- choices[selection]
  }
  input_file <- dplyr::case_when (
# official run labels
    input_file == "GBR4-v2.0" ~ "http://dapds00.nci.org.au/thredds/dodsC/fx3/gbr4_v2/gbr4_simple_2018-10.nc",
    input_file == "GBR4_BGC-v2.0 Chyb Dcrt" ~ "http://dapds00.nci.org.au/thredds/dodsC/fx3/gbr4_bgc_GBR4_H2p0_B2p0_Chyd_Dcrt/gbr4_bgc_simple_2016-07.nc",
    input_file == "GBR4_BGC-v2.0 Chyb Dnrt" ~ "http://dapds00.nci.org.au/thredds/dodsC/fx3/gbr4_bgc_GBR4_H2p0_B2p0_Chyd_Dnrt/gbr4_bgc_simple_2017-11.nc",
    input_file == "GBR4_BGC-v2.0 Cpre Dcrt" ~ "http://dapds00.nci.org.au/thredds/dodsC/fx3/gbr4_bgc_GBR4_H2p0_B2p0_Cpre_Dcrt/gbr4_bgc_simple_2016-06.nc",
    input_file == "GBR4 rivers_2.0  Dnrt   GBR4 passive" ~ "http://dapds00.nci.org.au/thredds/dodsC/fx3/gbr4_2.0_rivers/gbr4_rivers_simple_2018-07.nc",
    input_file == "GBR1-v2.0" ~ "http://dapds00.nci.org.au/thredds/dodsC/fx3/gbr1_2.0/gbr1_simple_2018-09-23.nc",
    input_file == "GBR1 rivers_2.0  Dnrt" ~ "http://dapds00.nci.org.au/thredds/dodsC/fx3/gbr1_2.0_rivers/gbr1_rivers_simple_2018-03.nc",
    input_file == "GBR4-v1.85" ~ "http://dapds00.nci.org.au/thredds/dodsC/fx3/gbr4/gbr4_simple_2016-03.nc.html",
    input_file == "GBR4_BGC-v926" ~ "http://dapds00.nci.org.au/thredds/dodsC/fx3/gbr4_bgc_926/gbr4_bgc_simple_2014-11.nc",
    input_file == "GBR4_BGC-v924" ~ "http://dapds00.nci.org.au/thredds/dodsC/fx3/gbr4_bgc_924/gbr4_bgc_simple_2016-09.nc",
    input_file == "GBR1-v1.71" ~ "http://dapds00.nci.org.au/thredds/dodsC/fx3/gbr1/gbr1_simple_2016-03-25.nc",
    input_file == "GBR1_BGC-v924" ~ "http://dapds00.nci.org.au/thredds/dodsC/fx3/gbr1_bgc_924/gbr1_bgc_simple_2018-08-21.nc",
# additional shortcuts
    input_file == "GBR4HD" ~ "http://dapds00.nci.org.au/thredds/dodsC/fx3/gbr4_v2/gbr4_simple_2018-10.nc",
    input_file == "hd" ~ "http://dapds00.nci.org.au/thredds/dodsC/fx3/gbr4_v2/gbr4_simple_2018-10.nc",
    input_file == "GBR4BGC" ~ "http://dapds00.nci.org.au/thredds/dodsC/fx3/gbr4_bgc_GBR4_H2p0_B2p0_Chyd_Dcrt/gbr4_bgc_simple_2016-07.nc",
    input_file == "bgc" ~ "http://dapds00.nci.org.au/thredds/dodsC/fx3/gbr4_bgc_GBR4_H2p0_B2p0_Chyd_Dcrt/gbr4_bgc_simple_2016-07.nc",
    input_file == "GBR4NRT" ~ "http://dapds00.nci.org.au/thredds/dodsC/fx3/gbr4_bgc_GBR4_H2p0_B2p0_Chyd_Dnrt/gbr4_bgc_simple_2017-11.nc",
    input_file == "nrt" ~ "http://dapds00.nci.org.au/thredds/dodsC/fx3/gbr4_bgc_GBR4_H2p0_B2p0_Chyd_Dnrt/gbr4_bgc_simple_2017-11.nc",
    input_file == "rivers" ~ "http://dapds00.nci.org.au/thredds/dodsC/fx3/gbr4_2.0_rivers/gbr4_rivers_simple_2018-07.nc",
    input_file == "GBR1HD" ~ "http://dapds00.nci.org.au/thredds/dodsC/fx3/gbr1_2.0/gbr1_simple_2018-09-23.nc",
    input_file == "GBR1BGC" ~ "http://dapds00.nci.org.au/thredds/dodsC/fx3/gbr1_bgc_924/gbr1_bgc_simple_2018-08-21.nc",
# otherwise, use input file name as entered
    TRUE ~ input_file
  )
    
  return (input_file)
}

#' Extracts time series at specified locations from eReefs model output files
#'
#' Create a time-series of values of one or more selected model output variables in a specified layer of the
#' water column or sediment store (by default, the surface layer), at a specified geographic location or 
#' supplied data frame of geocoordinates within the model domain.  See also get_ereefs_depth_integrated_ts() to 
#' extract depth-integrated values and get_ereefs_depth_specified_ts() to extract values at a specified depth 
#' below the free (tidally moving) surface. Barbara Robson (AIMS).
#'
#' If you run into memory constraints, consider grouping points to be extracted within regions, and calling this once
#' for each region.
#'
#' @return a data frame containing the dates and values of extracted variables (for a single geolocation) or a
#'        list of such data frames (one list item per location) if there are multiple locations.
#' @param var_names either a single character value or a vector specifying the short names for variables that you 
#'        want from the netcdf file. Defaults to c('Chl_a_sum', 'TN').
#' @param location_latlon is a data frame containing the decimal latitude and longitude of a single desired location, or a vector containing
#'        a single latitude and longitude location. If you want to specify an x-y grid coordinate instead of a latitude and longitude, you 
#'        can: to do this, is.integer(location_latlon) must be TRUE. Defaults to c(-23.39189, 150.88852).
#' @param layer is the vertical grid layer to extract, or 'surface' to get the surface value, 'bottom' to get the
#'        value in the cell at the bottom of the water column, or 'integrated' to get a depth-integrated (mean) value.
#'        Defaults to 'surface'. Use get_ereefs_depth_specified_ts() instead if you want to specify a depth 
#'        below the free surface instead of a layer number.
#' @param start_date date  for animation. Can either be a) a vector of integers in the format 
#'	       c(year, month, day); b) a date obtained e.g. from as.Date(); or c) a character string 
#'        formatted for input to as.Date(). Defaults to c(2010,12,31).
#' @param end date for animation. Can either be a) a vector of integers in the format 
#'	       c(year, month, day); b) a date obtained e.g. from as.Date(); or c) a character string 
#'        formatted for input to as.Date(). Defaults to c(2016,10,31).
#' @param input_file is the URI or file location of any of the EMS output files, 
#'        Defaults to a menu selection. Set to "choices" to see some other pre-defined options that
#'        can be used (codenames as used in https://research.csiro.au/ereefs/models/model-outputs/access-to-raw-model-output/ )
#' @param input_grid Name of the locally-stored or opendap-served netcdf file that contains the grid
#'      coordinates for the top and bottom of each layer (z_grid). If needed (i.e. for a depth-integrated value or bottom layer)
#'      but not specified, the function will first look for z_grid in the first INPUT_STEM file, and if not found, 
#'      will check whether the size of the variables in the input file corresponds to the size expected for GBR4 or GBR1, and 
#'      load an appropriate z grid from data files stored in this package. Alternatively, you can provide the location of a full 
#'      (not simple-format) ereefs netcdf output file such as 
#'      "http://dapds00.nci.org.au/thredds/dodsC/fx3/gbr4_hydro_all/gbr4_all_2016-09.nc"
#' @param eta_stem The URI or file location of the model output files that contains the surface elevation (eta), minus the
#'       date components of the filename in the case of GBR1 or GBR4 files, and ommitting the file extension, ".nc". Needed
#'       only if eta is not in the files indicated by input_stem (e.g. some GBR1 bgc files) and you are asking for a depth-integrated
#'       time-series or a depth-specified (relative to the surface) time-series
#' @param override_positive Reverse the value of the "positive" attribute of botz for BGC files, assuming that it is
#'       incorrect. Default FALSE. Not normally needed.
#' @export
#' @examples
#' get_ereefs_ts('Chl_a_sum', location_latlon=data.frame(latitide=c(-23.39189,-18), longitude=c(150.88852, 147.5)), layer='surface', start_date=c(2010,12,31),end_date=c(2011,1,5), input_file=2)
#' get_ereefs_ts(var_names=c('Tricho_N', 'DIP', 'TP'), location_latlon=data.frame(latitide=c(-23.39189,-18), longitude=c(150.88852, 147.5)), layer='bottom', start_date="2012-07-03",end_date="2012-07-05", input_file='GBR4_BGC-v2.0 Chyb Dcrt')
#' get_ereefs_ts(var_names=c('ZooL_N', 'ZooS_N'), location_latlon=data.frame(latitide=c(-23.39189,-18), longitude=c(150.88852, 147.5)), layer=45, start_date=c(2010,12,31),end_date=c(2011,1,5), input_file="http://dapds00.nci.org.au/thredds/dodsC/fx3/gbr4_bgc_GBR4_H2p0_B2p0_Cpre_Dcrt/gbr4_bgc_simple_2016-06.nc")

get_ereefs_ts <- function(var_names=c('Chl_a_sum', 'TN'), 
                          location_latlon=c(-23.39189, 150.88852), 
                          layer='surface', 
		                    start_date = c(2010,12,31), 
		                    end_date = c(2016,10,31), 
                          input_file = "menu",
                          #input_file = "http://dapds00.nci.org.au/thredds/dodsC/fx3/gbr4_bgc_GBR4_H2p0_B2p0_Chyd_Dcrt/gbr4_bgc_simple_2010-01.nc",
                          input_grid = NA,
                          eta_stem = NA,
                          override_positive=FALSE,
                          verbosity = 1)
{
  input_file <- substitute_filename(input_file)
  if (layer=='integrated') return(get_ereefs_depth_integrated_ts(var_names, location_latlon, start_date, end_date, input_file, input_grid, eta_stem, override_positive))
  if (layer=='bottom') return(get_ereefs_bottom_ts(var_names, location_latlon, start_date, end_date, input_file, input_grid, eta_stem, override_positive))

  # Check whether netcdf output files are daily (case 1), monthly (case 4) or something else (case 0)
  ereefs_case <- get_ereefs_case(input_file)
  input_stem <- get_file_stem(input_file)
  check_platform_ok(input_stem)

  # Dates to plot
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

  var_list <- paste(var_names, collapse=",")

  if (ereefs_case == 4) {
      input_file <- paste0(input_stem, format(as.Date(paste(start_year, start_month, 1, sep='-')), '%Y-%m'), 
			  '.nc')
	nc <- ncdf4::nc_open(input_file)
	if (!is.null(nc$var[['t']])) { 
	    ds <- as.Date(ncdf4::ncvar_get(nc, "t"), origin = as.Date("1990-01-01"))
        } else {
	    ds <- as.Date(ncdf4::ncvar_get(nc, "time"), origin = as.Date("1990-01-01"))
	}
	ncdf4::nc_close(nc)
        blank_length <- as.numeric(end_date - start_date + 1) / as.numeric(ds[2] - ds[1])
			  # '.nc?latitude,longitude')
  } else if (ereefs_case == 1) {
      input_file <- paste0(input_stem, format(as.Date(paste(start_year, start_month, start_day, sep='-')), '%Y-%m-%d'), 
			  '.nc')
      blank_length <- end_date - start_date + 1
			  # '.nc?latitude,longitude')
  } else {
      input_file <- input_file
      nc <- ncdf4::nc_open(input_file)
      if (!is.null(nc$var[['t']])) {
        ds <- as.Date(ncdf4::ncvar_get(nc, "t"), origin = as.Date("1990-01-01"))
      } else {
        ds <- as.Date(ncdf4::ncvar_get(nc, "time"), origin = as.Date("1990-01-01"))
      }
      blank_length <- as.numeric(end_date - start_date + 1) / as.numeric(ds[2] - ds[1])
      ncdf4::nc_close(nc)
  }

  if (is.null(dim(location_latlon))) {
     location_latlon <- array(location_latlon, c(1,2))
  }
  if (is.integer(location_latlon)) {
     # We have specified grid coordinates rather than geocoordinates
     location_grid <- location_latlon
  } else { 
    # We have geocoordinates. Find the nearest grid-points to the sampling location

    # First, get the model grid
    nc <- ncdf4::nc_open(input_file)
    if (is.null(nc$var[['latitude']])) {
      # Not a simple format netcdf file, so assume it's a full EMS netcdf file.
      latitude <- ncdf4::ncvar_get(nc, "y_centre")
      longitude <- ncdf4::ncvar_get(nc, "x_centre")
    } else { 
      # Simple format netcdf file
      latitude <- ncdf4::ncvar_get(nc, "latitude")
      longitude <- ncdf4::ncvar_get(nc, "longitude")
    }
    ncdf4::nc_close(nc)

    if (is.null(dim(location_latlon))) {
       # Just one location
       grid_index <- (latitude - location_latlon[1])^2 + (longitude - location_latlon[2])^2 
       grid_index <- which.min(grid_index) 
    } else { 
       # Multiple locations
       if (class(location_latlon) != "data.frame") {
          # location_latlon has been provided as an array/matrix. Coerce it into a data frame for consistency.
          location_latlon <- data.frame(latitude = location_latlon[,1], longitude = location_latlon[,2])
       }
       grid_index <- apply(location_latlon,1, function(ll) which.min((latitude - ll[1])^2 + (longitude - ll[2])^2)) 
    }
    location_grid <- cbind(floor((grid_index + dim(latitude)[1]-1)/dim(latitude)[1]), 
                       (grid_index+dim(latitude)[1]-1)%%dim(latitude)[1] + 1)
  }
  numpoints <- dim(location_grid)[1]
  # Find the outer grid coordinates of the area that we need to extract from netcdf files to encompass all
  # provided geocoordinate points
  startv <- c(min(location_grid[,2]), min(location_grid[, 1]))
  countv <- c(max(location_grid[,2]), max(location_grid[, 1])) - startv + 1

  # Adjust grid locations so that they are relative to the region to be extracted instead of the whole model domain
  location_grid <- t(t(location_grid) - c(startv[2], startv[1])) + 1
  location_grid <- cbind(location_grid[,2], location_grid[,1])

  # Update grid_index so that it is also relative
  grid_index <- (location_grid[,2] - 1) * countv[1] + location_grid[,1]

  # check whether all points are within a single model grid row or column, and adjust indices accordingly
  if (countv[2] == 1) { 
   location_grid <- location_grid[,1]
  } else if (countv[1] == 1) {
   location_grid <- location_grid[,2]
  }

  # Initialise
  blanks <- rep(NA, blank_length)
  ts_frame <- array(NA, c(blank_length, length(var_names)+1, numpoints))
  colnames(ts_frame) <- c("date", var_names)

  # Loop through monthly eReefs files to extract the data
  ndims <- rep(NA, length(var_names))
  layer_actual <- rep(NA, length(var_names))
  i <- 0
  mcount <- 0
  if (verbosity>0) pb <- txtProgressBar(min = 0, max = length(mths), style = 3)
  for (month in mths) {
    mcount <- mcount + 1
    year <- years[mcount]
    if (mcount == 1) {
       from_day <- start_day
    } else {
       from_day <- 1
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
    if (ereefs_case == 4) { 
       fileslist <- 1
       input_file <- paste0(input_stem, format(as.Date(paste(year, month, 1, sep="-")), '%Y-%m'), '.nc')
	    day_count <- day_count / as.numeric(ds[2]-ds[1])
    } else if (ereefs_case == 1) {
	    fileslist <- from_day:(from_day+day_count-1)
	    from_day <- 1
	    day_count <- 1
    } else { 
       from_day <- as.integer((as.Date(paste(year, month, from_day, sep="-")) - ds[1])/as.numeric(ds[2]-ds[1])) + 1
	    if (from_day<1) from_day <-1
	    day_count <- day_count / as.numeric(ds[2]-ds[1])
	    fileslist <- 1
    }
    for (dcount in fileslist) {
      if (ereefs_case == 1) {
	    input_file <- paste0(input_stem, format(as.Date(paste(year, month, dcount, sep="-")), '%Y-%m-%d'), '.nc')
      }
      #input_file <- paste0(input_file, '?', var_list, ',time')
      nc <- ncdf4::nc_open(input_file)
      if (ereefs_case > 0) {
          if (!is.null(nc$var[['t']])) {
            d <- as.Date(ncdf4::ncvar_get(nc, "t"), origin = as.Date("1990-01-01"))[from_day:(from_day+day_count-1)]
          } else {
            d <- as.Date(ncdf4::ncvar_get(nc, "time"), origin = as.Date("1990-01-01"))[from_day:(from_day+day_count-1)]
          }
      } else { 
         d <- ds[from_day:(from_day + day_count - 1)]
      }
      im1 = i+1
      i <- i + length(d)

      ts_frame[im1:i,"date",] <- d
      for (j in 1:length(var_names)) {
          if (is.na(ndims[j])) {
             # We don't yet know the dimensions of the variable, so let's get them
              dims <- nc$var[[var_names[j]]][['size']]
              if (is.null(dims)) stop(paste(var_names[j], ' not found in netcdf file.')) 
              ndims[j] <- length(dims)
              if (layer == 'surface') {
               layer_actual[j] <- dims[3]
            } else {
               layer_actual[j] <- layer
            }
           }
          if (ndims[j] == 4) {
             wc <- ncdf4::ncvar_get(nc, var_names[j], start=c(startv,layer_actual[j],from_day), count=c(countv,1,day_count))
             #ts_frame[im1:i, j+1] <- ncdf4::ncvar_get(nc, var_names[j], start=c(location_grid[2],location_grid[1],layer_actual[j],from_day), count=c(1,1,1,day_count))
          } else {
             wc <- ncdf4::ncvar_get(nc, var_names[j], start=c(startv,from_day), count=c(countv,day_count))
             #ts_frame[im1:i, j+1] <- ncdf4::ncvar_get(nc, var_names[j], start=c(location_grid[2],location_grid[1],from_day), count=c(1,1,day_count))
          }
          wc <- array(wc, c(countv[1]*countv[2], day_count))
          ts_frame[im1:i, j+1,] <- t(wc[grid_index,])
      }
      ncdf4::nc_close(nc)
      if (verbosity>0) setTxtProgressBar(pb,mcount)
    }
  }
  if (verbosity>0) close(pb)
  ts_frame <- lapply(seq(dim(ts_frame)[3]), function(x) data.frame(date=as.Date(as.vector(ts_frame[ ,1, 1]), origin="1970-01-01"), ts_frame[ ,2:dim(ts_frame)[2] , x])) 
  if (numpoints == 1) ts_frame <- ts_frame[[1]]
  return(ts_frame)
}

#' Extracts time-series of selected variables from the bottom water-column cell at specified locations from eReefs output files
#'
#' See also get_ereefs_ts() to extract from a specified layer instead of a depth-integrated value.
#' By Barbara Robson (AIMS).
#'
#' @return a data frame containing the dates and values of extracted variables
#' @param var_names either a single character value or a vector specifying the short names for variables that you 
#'        want from the netcdf file. Defaults to c('Chl_a_sum', 'TN').
#' @param location_latlon is a vector containing the decimal latitude and longitude of the desired location. If 
#'        you want to specify an x-y grid coordinate instead of a latitude and longitude, you can: to do this, 
#'        is.integer(location_latlon) must be TRUE. Defaults to c(-23.39189, 150.88852).
#' @param start_date date  for animation. Can either be a) a vector of integers in the format 
#'	c(year, month, day); b) a date obtained e.g. from as.Date(); or c) a character string 
#'      formatted for input to as.Date(). Defaults to c(2010,12,31).
#' @param end date for animation. Can either be a) a vector of integers in the format 
#'	c(year, month, day); b) a date obtained e.g. from as.Date(); or c) a character string 
#'      formatted for input to as.Date(). Defaults to c(2016,12,31).
#' @param input_file is the URI or file location of any of the EMS output files, 
#'        Defaults to a menu selection. Set to "choices" to see some other pre-defined options that
#'        can be used (codenames as used in https://research.csiro.au/ereefs/models/model-outputs/access-to-raw-model-output/ )
#' @param input_grid Name of the locally-stored or opendap-served netcdf file that contains the grid
#'      coordinates for the top and bottom of each layer (z_grid). If not specified, the function will first look for
#'      z_grid can be found in the first INPUT_STEM file, and if not found, will check whether the size 
#'      of the variables in the input file corresponds to the size expected for GBR4 or GBR1, and load an appropriate 
#'      z grid from data files stored in this package. Alternatively, you can provide the location of a full 
#'      (not simple-format) ereefs netcdf output file such as 
#'      "http://dapds00.nci.org.au/thredds/dodsC/fx3/gbr4_hydro_all/gbr4_all_2016-09.nc"
#'      "http://dapds00.nci.org.au/thredds/dodsC/fx3/gbr4_hydro_all/gbr4_all_2016-09.nc".
#' @param eta_stem The URI or file location of the model output files that contains the surface elevation (eta), minus the
#'       date components of the filename in the case of GBR1 or GBR4 files, and ommitting the file extension, ".nc". Needed
#'       only if eta is not in the files indicated by input_stem (e.g. some GBR1 bgc files).
#' @param override_positive Reverse the value of the "positive" attribute of botz for BGC files, assuming that it is
#'       incorrect. Default FALSE
#' @export
#' @examples
#' get_ereefs_bottom_ts(c('Chl_a_sum', 'NH4'), location_latlon=data.frame(latitide=-23.39189, longitude=150.88852), layer='surface', start_date=c(2010,12,31),end_date=c(2011,1,5), input_file=2)
get_ereefs_bottom_ts <- function(var_names=c('Chl_a_sum', 'TN'), 
                          location_latlon=c(-23.39189, 150.88852), 
		                    start_date = c(2010,12,31), 
		                    end_date = c(2016,12,31), 
                          input_file = "menu",
                          #input_file = "http://dapds00.nci.org.au/thredds/dodsC/fx3/gbr4_bgc_GBR4_H2p0_B2p0_Chyd_Dcrt/gbr4_bgc_simple_2010-01.nc",
			                 input_grid = NA,
			                 eta_stem = NA,
			                 override_positive=FALSE, 
                          verbosity = 1)
{
  input_file <- substitute_filename(input_file)

  # Check whether this is a GBR1 or GBR4 ereefs file, or something else
  ereefs_case <- get_ereefs_case(input_file)
  input_stem <- get_file_stem(input_file)
  check_platform_ok(input_stem)
  z_grid <- get_ereefs_grids(input_file, input_grid)[['z_grid']]

  # Dates to plot
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

  var_list <- paste(var_names, collapse=",")

  if (ereefs_case == 4) {
      input_file <- paste0(input_stem, format(as.Date(paste(start_year, start_month, 1, sep='-')), '%Y-%m'), 
			  '.nc')
      if (!is.na(eta_stem)) etafile  <- paste0(eta_stem, format(as.Date(paste(start_year, start_month, 1, sep='-')), '%Y-%m'), 
			  '.nc')
	nc <- ncdf4::nc_open(input_file)
	if (!is.null(nc$var[['t']])) { 
	    ds <- as.Date(ncdf4::ncvar_get(nc, "t"), origin = as.Date("1990-01-01"))
        } else {
	    ds <- as.Date(ncdf4::ncvar_get(nc, "time"), origin = as.Date("1990-01-01"))
	}
	ncdf4::nc_close(nc) 
   blank_length <- as.numeric(end_date - start_date + 1) / as.numeric(ds[2] - ds[1])
  } else if (ereefs_case == 1) {
      input_file <- paste0(input_stem, format(as.Date(paste(start_year, start_month, start_day, sep='-')), '%Y-%m-%d'), 
			  '.nc')
      if (!is.na(eta_stem)) etafile <- paste0(eta_stem, format(as.Date(paste(start_year, start_month, start_day, sep='-')), '%Y-%m-%d'), 
			  '.nc') 
      blank_length <- end_date - start_date + 1
  } else {
      input_file <- input_file
      if (!is.na(eta_stem)) etafile <- paste0(eta_stem, '.nc')
      nc <- ncdf4::nc_open(input_file)
      if (!is.null(nc$var[['t']])) {
        ds <- as.Date(ncdf4::ncvar_get(nc, "t"), origin = as.Date("1990-01-01"))
      } else {
        ds <- as.Date(ncdf4::ncvar_get(nc, "time"), origin = as.Date("1990-01-01"))
      }
      blank_length <- as.numeric(end_date - start_date + 1) / as.numeric(ds[2] - ds[1])
      ncdf4::nc_close(nc)
  }

  # Initialise
  blanks <- rep(NA, blank_length)
  ts_frame <- data.frame(as.Date(blanks), array(blanks, dim=c(length(blanks), length(var_names))))
  names(ts_frame) <- c("date", var_names)

  nc <- ncdf4::nc_open(input_file)
  if (!is.na(eta_stem)) nc3 <- ncdf4::nc_open(etafile)

  if (is.integer(location_latlon)) {
     location_grid <- location_latlon
  } else { 
    # Find the nearest grid-points to the sampling location
    if (is.null(nc$var[['latitude']])) {
      latitude <- ncdf4::ncvar_get(nc, "y_centre")
      longitude <- ncdf4::ncvar_get(nc, "x_centre")
    } else { 
      latitude <- ncdf4::ncvar_get(nc, "latitude")
      longitude <- ncdf4::ncvar_get(nc, "longitude")
    }
    tmp <- (latitude - location_latlon[1])^2 + (longitude - location_latlon[2])^2 
    tmp <- which.min(tmp) 
    location_grid <- c(floor((tmp+dim(latitude)[1]-1)/dim(latitude)[1]),
		       (tmp+dim(latitude)[1]-1)%%dim(latitude)[1] + 1)
  }
  zat <- ncdf4::ncatt_get(nc, "botz")
  if (!is.null(zat$positive)) {
    if (zat$positive=="down") zsign <- -1 else zsign <- 1
    if (override_positive) zsign <- 1
  } else {
   zsign <-1
    if (override_positive) zsign <- -1
  }
  botz <- zsign * as.numeric(ncdf4::ncvar_get(nc, "botz", start=c(location_grid[2], location_grid[1]), count=c(1,1)))
  ncdf4::nc_close(nc)

  # Loop through monthly eReefs files to extract the data
  i <- 0
  mcount <- 0
  if (verbosity>0) pb <- txtProgressBar(min = 0, max = length(mths), style = 3)
  for (month in mths) {
     mcount <- mcount + 1
     year <- years[mcount]
     if (mcount == 1) {
        from_day <- start_day
     } else {
        from_day <- 1
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
     if (ereefs_case == 4) { 
        fileslist <- 1
        input_file <- paste0(input_stem, format(as.Date(paste(year, month, 1, sep="-")), '%Y-%m'), '.nc')
	     day_count <- day_count / as.numeric(ds[2]-ds[1])
        if (!is.na(eta_stem)) etafile <- paste0(eta_stem, format(as.Date(paste(year, month, 1, sep="-")), '%Y-%m'), '.nc')
     } else if (ereefs_case == 1) {
        fileslist <- from_day:(from_day+day_count-1)
        from_day <- 1
        day_count <- 1
     } else {
       from_day <- as.integer((as.Date(paste(year, month, from_day, sep="-")) - ds[1])/as.numeric(ds[2]-ds[1])) + 1
	    if (from_day<1) from_day <-1
	    day_count <- day_count / as.numeric(ds[2]-ds[1])
	    fileslist <- 1
     }
     for (dcount in fileslist) {
        if (ereefs_case == 1) {
	      input_file <- paste0(input_stem, format(as.Date(paste(year, month, dcount, sep="-")), '%Y-%m-%d'), '.nc')
         if (!is.na(eta_stem)) etafile <- paste0(eta_stem, format(as.Date(paste(year, month, dcount, sep="-")), '%Y-%m-%d'), '.nc')
        }
        #input_file <- paste0(input_file, '?', var_list, ',time,eta')
        nc <- ncdf4::nc_open(input_file)
        if (!is.na(eta_stem)) nc3 <- ncdf4::nc_open(etafile)
        # Get dates
        if (ereefs_case > 0 ) { 
           if (!is.null(nc$var[['t']])) { 
              ds <- as.Date(ncdf4::ncvar_get(nc, "t"), origin = as.Date("1990-01-01")) 
           } else { 
              ds <- as.Date(ncdf4::ncvar_get(nc, "time"), origin = as.Date("1990-01-01")) 
           }
	        d <- ds[from_day:(from_day + day_count - 1)]
        } else {
	        d <- ds[from_day:(from_day + day_count - 1)]
        }
        if (!is.null(nc$var[['eta']])) { 
           eta <- ncdf4::ncvar_get(nc, 'eta', start=c(location_grid[2], location_grid[1],from_day), count=c(1,1,day_count)) 
        } else { 
           eta <- ncdf4::ncvar_get(nc3, 'eta', start=c(location_grid[2], location_grid[1],from_day), count=c(1,1,day_count)) 
        }
        im1 = i+1
        i <- i + length(d)
        ts_frame$date[im1:i] <- d
        z <- array(z_grid[2:length(z_grid)], dim=c(length(z_grid)-1, length(eta)))
        zm1 <- array(z_grid[1:(length(z_grid)-1)], dim=c(length(z_grid)-1, length(eta)))
        eta2 <- t(array(eta, dim=c(length(eta), length(z_grid)-1)))
        dz <- 0 * z
        wet <- (eta2 > zm1) & (z > botz)           # There is water in this layer
        bottom <- (wet & (zm1<=botz))              # The bottom intersects this layer
        dz[bottom] <- 1

        for (j in 1:length(var_names)) {
          wc <- ncdf4::ncvar_get(nc, var_names[j], start=c(location_grid[2],location_grid[1],1,from_day), count=c(1,1,-1,day_count))
	       if (dim(dz)[2] == 1) wc <- array(wc, dim=dim(dz))
          # take the depth-integrated average over the water column
          ts_frame[im1:i, j+1] <- colSums(dz * wc, na.rm=TRUE)
        }
        ncdf4::nc_close(nc)
        if (!is.na(eta_stem)) ncdf4::nc_close(nc3)
        if (verbosity>0) setTxtProgressBar(pb,mcount)
    }
  }
  if (verbosity>0) close(pb)
  return(ts_frame)
}

#' Extracts depth-integrated time-series of selected variables at specified locations from eReefs output files
#'
#' See also get_ereefs_ts() to extract from a specified layer instead of a depth-integrated value.
#' By Barbara Robson (AIMS).
#'
#' @return a data frame containing the dates and values of extracted variables
#' @param var_names either a single character value or a vector specifying the short names for variables that you 
#'        want from the netcdf file. Defaults to c('Chl_a_sum', 'TN').
#' @param location_latlon is a vector containing the decimal latitude and longitude of the desired location. If 
#'        you want to specify an x-y grid coordinate instead of a latitude and longitude, you can: to do this, 
#'        is.integer(location_latlon) must be TRUE. Defaults to c(-23.39189, 150.88852).
#' @param start_date date  for animation. Can either be a) a vector of integers in the format 
#'	c(year, month, day); b) a date obtained e.g. from as.Date(); or c) a character string 
#'      formatted for input to as.Date(). Defaults to c(2010,12,31).
#' @param end date for animation. Can either be a) a vector of integers in the format 
#'	c(year, month, day); b) a date obtained e.g. from as.Date(); or c) a character string 
#'      formatted for input to as.Date(). Defaults to c(2016,12,31).
#' @param input_file is the URI or file location of any of the EMS output files, 
#'        Defaults to a menu selection. Set to "choices" to see some other pre-defined options that
#'        can be used (codenames as used in https://research.csiro.au/ereefs/models/model-outputs/access-to-raw-model-output/ )
#' @param input_grid Name of the locally-stored or opendap-served netcdf file that contains the grid
#'      coordinates for the top and bottom of each layer (z_grid). If not specified, the function will first look for
#'      z_grid can be found in the first INPUT_STEM file, and if not found, will check whether the size 
#'      of the variables in the input file corresponds to the size expected for GBR4 or GBR1, and load an appropriate 
#'      z grid from data files stored in this package. Alternatively, you can provide the location of a full 
#'      (not simple-format) ereefs netcdf output file such as 
#'      "http://dapds00.nci.org.au/thredds/dodsC/fx3/gbr4_hydro_all/gbr4_all_2016-09.nc"
#'      "http://dapds00.nci.org.au/thredds/dodsC/fx3/gbr4_hydro_all/gbr4_all_2016-09.nc".
#' @param eta_stem The URI or file location of the model output files that contains the surface elevation (eta), minus the
#'       date components of the filename in the case of GBR1 or GBR4 files, and ommitting the file extension, ".nc". Needed
#'       only if eta is not in the files indicated by input_stem (e.g. some GBR1 bgc files).
#' @param override_positive Reverse the value of the "positive" attribute of botz for BGC files, assuming that it is
#'       incorrect. Default FALSE
#' @export
#' @examples
#' get_ereefs_depth_integrated_ts(c('Chl_a_sum', 'NH4'), location_latlon=data.frame(latitide=-23.39189, longitude=150.88852), layer='surface', start_date=c(2010,12,31),end_date=c(2011,1,5), input_file=2)
get_ereefs_depth_integrated_ts <- function(var_names=c('Chl_a_sum', 'TN'), 
                          location_latlon=c(-23.39189, 150.88852), 
		                    start_date = c(2010,12,31), 
		                    end_date = c(2016,12,31), 
                          input_file = "menu",
			                 input_grid = NA,
			                 eta_stem = NA,
			                 override_positive=FALSE,
                          verbosity = 1)
{

  input_file <- substitute_filename(input_file)
  # Check whether this is a GBR1 or GBR4 ereefs file, or something else
  ereefs_case <- get_ereefs_case(input_file)
  input_stem <- get_file_stem(input_file)
  check_platform_ok(input_stem)
  z_grid <- get_ereefs_grids(input_file, input_grid)[['z_grid']]

  # Dates to plot
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

  var_list <- paste(var_names, collapse=",")

  if (ereefs_case == 4) {
      input_file <- paste0(input_stem, format(as.Date(paste(start_year, start_month, 1, sep='-')), '%Y-%m'), 
			  '.nc')
      if (verbosity>1) print(input_file)
      if (!is.na(eta_stem)) etafile  <- paste0(eta_stem, format(as.Date(paste(start_year, start_month, 1, sep='-')), '%Y-%m'), 
			  '.nc')
	nc <- ncdf4::nc_open(input_file)
	if (!is.null(nc$var[['t']])) { 
	    ds <- as.Date(ncdf4::ncvar_get(nc, "t"), origin = as.Date("1990-01-01"))
        } else {
	    ds <- as.Date(ncdf4::ncvar_get(nc, "time"), origin = as.Date("1990-01-01"))
	}
	ncdf4::nc_close(nc)
        blank_length <- as.numeric(end_date - start_date + 1) / as.numeric(ds[2] - ds[1])
			  # '.nc?latitude,longitude')
  } else if (ereefs_case == 1) {
      input_file <- paste0(input_stem, format(as.Date(paste(start_year, start_month, start_day, sep='-')), '%Y-%m-%d'), 
			  '.nc')
      if (verbosity>1) print(input_file)
      if (!is.na(eta_stem)) etafile <- paste0(eta_stem, format(as.Date(paste(start_year, start_month, start_day, sep='-')), '%Y-%m-%d'), 
			  '.nc')
      blank_length <- end_date - start_date + 1
			  # '.nc?latitude,longitude')
  } else {
      input_file <- input_file
      if (verbosity>1) print(input_file)
      if (!is.na(eta_stem)) etafile <- paste0(eta_stem, '.nc')
      nc <- ncdf4::nc_open(input_file)
      if (!is.null(nc$var[['t']])) {
        ds <- as.Date(ncdf4::ncvar_get(nc, "t"), origin = as.Date("1990-01-01"))
      } else {
        ds <- as.Date(ncdf4::ncvar_get(nc, "time"), origin = as.Date("1990-01-01"))
      }
      blank_length <- as.numeric(end_date - start_date + 1) / as.numeric(ds[2] - ds[1])
      ncdf4::nc_close(nc)
  }

  if (dim(location_latlon)[1] > 1) stop('Currently, get_ereefs_depth_integrated_ts() only supports a single location. This is on my to-do list to fix in future. Let me know if you would like this feature. b.robson@aims.gov.au')

  nc <- ncdf4::nc_open(input_file)
  if (!is.na(eta_stem)) nc3 <- ncdf4::nc_open(etafile)

  if (is.null(dim(location_latlon))) {
     location_latlon <- array(location_latlon, c(1,2))
  }
  if (is.integer(location_latlon)) {
     # We have specified grid coordinates rather than geocoordinates
     location_grid <- location_latlon
  } else { 
    # We have geocoordinates. Find the nearest grid-points to the sampling location

    # First, get the model grid
    nc <- ncdf4::nc_open(input_file)
    if (is.null(nc$var[['latitude']])) {
      # Not a simple format netcdf file, so assume it's a full EMS netcdf file.
      latitude <- ncdf4::ncvar_get(nc, "y_centre")
      longitude <- ncdf4::ncvar_get(nc, "x_centre")
    } else { 
      # Simple format netcdf file
      latitude <- ncdf4::ncvar_get(nc, "latitude")
      longitude <- ncdf4::ncvar_get(nc, "longitude")
    }
    ncdf4::nc_close(nc)

    if (is.null(dim(location_latlon))) {
       # Just one location
       grid_index <- (latitude - location_latlon[1])^2 + (longitude - location_latlon[2])^2 
       grid_index <- which.min(grid_index) 
    } else { 
       # Multiple locations
       if (class(location_latlon) != "data.frame") {
          # location_latlon has been provided as an array/matrix. Coerce it into a data frame for consistency.
          location_latlon <- data.frame(latitude = location_latlon[,1], longitude = location_latlon[,2])
       }
       grid_index <- apply(location_latlon,1, function(ll) which.min((latitude - ll[1])^2 + (longitude - ll[2])^2)) 
    }
    location_grid <- cbind(floor((grid_index + dim(latitude)[1]-1)/dim(latitude)[1]), 
                       (grid_index+dim(latitude)[1]-1)%%dim(latitude)[1] + 1)
  }
  numpoints <- dim(location_grid)[1]
  # Find the outer grid coordinates of the area that we need to extract from netcdf files to encompass all
  # provided geocoordinate points
  startv <- c(min(location_grid[,2]), min(location_grid[, 1]))
  countv <- c(max(location_grid[,2]), max(location_grid[, 1])) - startv + 1

  # Adjust grid locations so that they are relative to the region to be extracted instead of the whole model domain
  location_grid <- t(t(location_grid) - c(startv[2], startv[1])) + 1
  location_grid <- cbind(location_grid[,2], location_grid[,1])

  # Update grid_index so that it is also relative
  grid_index <- (location_grid[,2] - 1) * countv[1] + location_grid[,1]

  # check whether all points are within a single model grid row or column, and adjust indices accordingly
  if (countv[2] == 1) { 
   location_grid <- location_grid[,1]
  } else if (countv[1] == 1) {
   location_grid <- location_grid[,2]
  }

  # Initialise
  blanks <- rep(NA, blank_length)
  ts_frame <- data.frame(as.Date(blanks), array(blanks, dim=c(length(blanks), length(var_names))))
  names(ts_frame) <- c("date", var_names)

  zat <- ncdf4::ncatt_get(nc, "botz")
  if (!is.null(zat$positive)) {
    if (zat$positive=="down") zsign <- -1 else zsign <- 1
    if (override_positive) zsign <- 1
  } else {
   zsign <-1
    if (override_positive) zsign <- -1
  }
  botz <- zsign * as.numeric(ncdf4::ncvar_get(nc, "botz", start=c(location_grid[2], location_grid[1]), count=c(1,1)))
  ncdf4::nc_close(nc)

  # Loop through monthly eReefs files to extract the data
  i <- 0
  mcount <- 0
  if (verbosity>0) pb <- txtProgressBar(min = 0, max = length(mths), style = 3)
  for (month in mths) {
     mcount <- mcount + 1
     year <- years[mcount]
     if (mcount == 1) {
        from_day <- start_day
     } else {
        from_day <- 1
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
     if (ereefs_case == 4) { 
        fileslist <- 1
        input_file <- paste0(input_stem, format(as.Date(paste(year, month, 1, sep="-")), '%Y-%m'), '.nc') 
        if (verbosity) print(input_file)
        day_count <- day_count / as.numeric(ds[2]-ds[1])
        if (!is.na(eta_stem)) etafile <- paste0(eta_stem, format(as.Date(paste(year, month, 1, sep="-")), '%Y-%m'), '.nc')
     } else if (ereefs_case == 1) {
        fileslist <- from_day:(from_day+day_count-1)
        from_day <- 1
        day_count <- 1
     } else {
       from_day <- as.integer((as.Date(paste(year, month, from_day, sep="-")) - ds[1])/as.numeric(ds[2]-ds[1])) + 1
	    if (from_day<1) from_day <-1
	    day_count <- day_count / as.numeric(ds[2]-ds[1])
	    fileslist <- 1
     }
     for (dcount in fileslist) {
        if (ereefs_case == 1) {
	      input_file <- paste0(input_stem, format(as.Date(paste(year, month, dcount, sep="-")), '%Y-%m-%d'), '.nc')
         if (!is.na(eta_stem)) etafile <- paste0(eta_stem, format(as.Date(paste(year, month, dcount, sep="-")), '%Y-%m-%d'), '.nc')
         if (verbosity>1) print(input_file)
        }
        #input_file <- paste0(input_file, '?', var_list, ',time,eta')
        nc <- ncdf4::nc_open(input_file)
        if (!is.na(eta_stem)) nc3 <- ncdf4::nc_open(etafile)
        if (ereefs_case > 0) {
          if (!is.null(nc$var[['t']])) {
            d <- as.Date(ncdf4::ncvar_get(nc, "t"), origin = as.Date("1990-01-01"))[from_day:(from_day+day_count-1)]
          } else {
            d <- as.Date(ncdf4::ncvar_get(nc, "time"), origin = as.Date("1990-01-01"))[from_day:(from_day+day_count-1)]
          }
        } else { 
           d <- ds[from_day:(from_day + day_count - 1)]
        }
        if (!is.null(nc$var[['eta']])) { 
          eta <- ncdf4::ncvar_get(nc, 'eta', start=c(location_grid[2], location_grid[1],from_day), count=c(1,1,day_count)) 
        } else {
          eta <- ncdf4::ncvar_get(nc3, 'eta', start=c(location_grid[2], location_grid[1],from_day), count=c(1,1,day_count)) 
        }
        im1 = i+1
        i <- i + length(d)
        ts_frame$date[im1:i] <- d
        z <- array(z_grid[2:length(z_grid)], dim=c(length(z_grid)-1, length(eta)))
        zm1 <- array(z_grid[1:(length(z_grid)-1)], dim=c(length(z_grid)-1, length(eta)))
        eta2 <- t(array(eta, dim=c(length(eta), length(z_grid)-1)))
        dz <- 0 * z
        wet <- (eta2 > zm1) & (z > botz)           # There is water in this layer
        bottom <- (wet & (zm1<botz))               # The bottom intersects this layer
        subsurface <- (wet & (eta2 > z))           # This is a wet layer that is not the surface layer
        surface <- (wet & !subsurface)             # This is the surface layer
        singlelayer <- (surface & bottom)          # This is both the surface and bottom layer
        surface <- (surface & !bottom)             # This is the surface layer but not the bottom layer
        fullywet <- (subsurface & !bottom)         # This is a wet layer that isn't the surface or bottom layer
        submergedbottom <- (bottom & !singlelayer) # This is the bottom layer but not the surface layer
  
        # Thickness of each layer
        dz[surface] <- eta2[surface] - zm1[surface]
        dz[fullywet] <- z[fullywet] - zm1[fullywet]
        dz[submergedbottom] <- z[submergedbottom] - botz
        dz[singlelayer] <- eta2[singlelayer] - botz
  

        for (j in 1:length(var_names)) {
          wc <- ncdf4::ncvar_get(nc, var_names[j], start=c(location_grid[2],location_grid[1],1,from_day), count=c(1,1,-1,day_count))
	       if (dim(dz)[2] == 1) wc <- array(wc, dim=dim(dz))
          # take the depth-integrated average over the water column
          ts_frame[im1:i, j+1] <- colSums(dz * wc, na.rm=TRUE) / colSums(dz)
        }
        ncdf4::nc_close(nc)
        if (!is.na(eta_stem)) ncdf4::nc_close(nc3)
        if (verbosity>0) setTxtProgressBar(pb,mcount)
    }
  }
  if (verbosity>0) close(pb)
  return(ts_frame)
}

#' Extracts time-series of selected variables from eReefs output files at a specified location and depth below the
#' surface.
#'
#' See also get_ereefs_ts() to extract from a specified layer instead of a depth-integrated value and 
#' get_ereefs_depth_integrated_ts() to calculate depth-integrated values. Barbara Robson (AIMS).
#'
#' @return a data frame containing the dates and values of extracted variables
#' @param var_names either a single character value or a vector specifying the short names for variables that you 
#'        want from the netcdf file. Defaults to c('Chl_a_sum', 'TN').
#' @param location_latlon A vector containing the decimal latitude and longitude of the desired location. If 
#'        you want to specify an x-y grid coordinate instead of a latitude and longitude, you can: to do this, 
#'        is.integer(location_latlon) must be TRUE. Defaults to c(-23.39189, 150.88852).
#' @param depth  Depth in metres below the surface. Default 1.0. If the bottom of the water is shallower than the specified depth,
#'               return values from the bottom of the water column.
#' @param start_date date  for animation. Can either be a) a vector of integers in the format 
#'	c(year, month, day); b) a date obtained e.g. from as.Date(); or c) a character string 
#'      formatted for input to as.Date(). Defaults to c(2010,12,31).
#' @param end date for animation. Can either be a) a vector of integers in the format 
#'	c(year, month, day); b) a date obtained e.g. from as.Date(); or c) a character string 
#'      formatted for input to as.Date(). Defaults to c(2016,12,31).
#' @param input_file is the URI or file location of any of the EMS output files, 
#'        Defaults to a menu selection. Set to "choices" to see some other pre-defined options that
#'        can be used (codenames as used in https://research.csiro.au/ereefs/models/model-outputs/access-to-raw-model-output/ )
#' @param input_grid Name of the locally-stored or opendap-served netcdf file that contains the grid
#'      coordinates for the top and bottom of each layer (z_grid). If not specified, the function will first look for
#'      z_grid can be found in the first INPUT_STEM file, and if not found, will check whether the size 
#'      of the variables in the input file corresponds to the size expected for GBR4 or GBR1, and load an appropriate 
#'      z grid from data files stored in this package. Alternatively, you can provide the location of a full 
#'      (not simple-format) ereefs netcdf output file such as 
#'      "http://dapds00.nci.org.au/thredds/dodsC/fx3/gbr4_hydro_all/gbr4_all_2016-09.nc"
#' @param eta_stem The URI or file location of the model output files that contains the surface elevation (eta), minus the
#'       date components of the filename in the case of GBR1 or GBR4 files, and ommitting the file extension, ".nc". Needed
#'       only if eta is not in the files indicated by input_stem (e.g. some GBR1 bgc files).
#' @export
#' @examples
#' get_ereefs_depth_specified_ts(c('Chl_a_sum', 'NH4'), depth=5.0, location_latlon=data.frame(latitide=-23.39189, longitude=150.88852), layer='surface', start_date=c(2010,12,31),end_date=c(2011,1,5), input_file=2)
get_ereefs_depth_specified_ts <- function(var_names=c('Chl_a_sum', 'TN'), 
                                          location_latlon=c(-23.39189, 150.88852), 
                                          depth = 1.0, 
                                          start_date = c(2010,12,31), 
                                          end_date = c(2016,12,31), 
                                          input_file = "menu", 
                                          input_grid = NA, 
                                          eta_stem = NA, 
                                          verbosity = 1)
{

  input_file <- substitute_filename(input_file)
  # Check whether this is a GBR1 or GBR4 ereefs file, or something else
  ereefs_case <- get_ereefs_case(input_file)
  input_stem <- get_file_stem(input_file)
  check_platform_ok(input_stem)
  z_grid <- get_ereefs_grids(input_file, input_grid)[['z_grid']]

  # Dates to plot
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

  var_list <- paste(var_names, collapse=",")

  if (ereefs_case == 4) {
      input_file <- paste0(input_stem, format(as.Date(paste(start_year, start_month, 1, sep='-')), '%Y-%m'), 
			  '.nc')
      if (!is.na(eta_stem)) etafile  <- paste0(eta_stem, format(as.Date(paste(start_year, start_month, 1, sep='-')), '%Y-%m'), 
			  '.nc')
	nc <- ncdf4::nc_open(input_file)
	if (!is.null(nc$var[['t']])) { 
	    ds <- as.Date(ncdf4::ncvar_get(nc, "t"), origin = as.Date("1990-01-01"))
        } else {
	    ds <- as.Date(ncdf4::ncvar_get(nc, "time"), origin = as.Date("1990-01-01"))
	}
	ncdf4::nc_close(nc)
        blank_length <- as.numeric(end_date - start_date + 1) / as.numeric(ds[2] - ds[1])
			  # '.nc?latitude,longitude')
  } else if (ereefs_case == 1) {
      input_file <- paste0(input_stem, format(as.Date(paste(start_year, start_month, start_day, sep='-')), '%Y-%m-%d'), 
			  '.nc')
      if (!is.na(eta_stem)) etafile <- paste0(eta_stem, format(as.Date(paste(start_year, start_month, start_day, sep='-')), '%Y-%m-%d'), 
			  '.nc')
      blank_length <- end_date - start_date + 1
			  # '.nc?latitude,longitude')
  } else {
      input_file <- input_file
      if (!is.na(eta_stem)) etafile <- paste0(eta_stem, '.nc')
      nc <- ncdf4::nc_open(input_file)
      if (!is.null(nc$var[['t']])) {
        ds <- as.Date(ncdf4::ncvar_get(nc, "t"), origin = as.Date("1990-01-01"))
      } else {
        ds <- as.Date(ncdf4::ncvar_get(nc, "time"), origin = as.Date("1990-01-01"))
      }
      blank_length <- as.numeric(end_date - start_date + 1) / as.numeric(ds[2] - ds[1])
      ncdf4::nc_close(nc)
  }

  # Initialise
  blanks <- rep(NA, blank_length)
  ts_frame <- data.frame(as.Date(blanks), array(blanks, dim=c(length(blanks), length(var_names))))
  names(ts_frame) <- c("date", var_names)

  nc <- ncdf4::nc_open(input_file)

  if (is.integer(location_latlon)) {
     location_grid <- location_latlon
  } else { 
    # Find the nearest grid-points to the sampling location
    if (is.null(nc$var[['latitude']])) {
      latitude <- ncdf4::ncvar_get(nc, "y_centre")
      longitude <- ncdf4::ncvar_get(nc, "x_centre")
    } else { 
      latitude <- ncdf4::ncvar_get(nc, "latitude")
      longitude <- ncdf4::ncvar_get(nc, "longitude")
    }
    tmp <- (latitude - location_latlon[1])^2 + (longitude - location_latlon[2])^2 
    tmp <- which.min(tmp) 
    location_grid <- c(floor((tmp+dim(latitude)[1]-1)/dim(latitude)[1]),
		       (tmp+dim(latitude)[1]-1)%%dim(latitude)[1] + 1)
  }
  zat <- ncdf4::ncatt_get(nc, "botz")
  if (!is.null(zat$positive)) {
	if (zat$positive=="down") zsign <- -1 else zsign <- 1
  } else {
	zsign <-1
  }
  botz <- zsign * as.numeric(ncdf4::ncvar_get(nc, "botz", start=c(location_grid[2], location_grid[1]), count=c(1,1)))
  ncdf4::nc_close(nc)

  # Loop through monthly eReefs files to extract the data
  i <- 0
  mcount <- 0
  if (verbosity>0) pb <- txtProgressBar(min = 0, max = length(mths), style = 3)
  for (month in mths) {
     mcount <- mcount + 1
     year <- years[mcount]
     if (mcount == 1) {
        from_day <- start_day
     } else {
        from_day <- 1
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
     if (ereefs_case == 4) { 
        fileslist <- 1
        input_file <- paste0(input_stem, format(as.Date(paste(year, month, 1, sep="-")), '%Y-%m'), '.nc')
	day_count <- day_count / as.numeric(ds[2]-ds[1])
     } else if (ereefs_case == 1) {
        fileslist <- from_day:(from_day+day_count-1)
        from_day <- 1
        day_count <- 1
     } else {
            from_day <- as.integer((as.Date(paste(year, month, from_day, sep="-")) - ds[1])/as.numeric(ds[2]-ds[1])) + 1
	    if (from_day<1) from_day <-1
	    day_count <- day_count / as.numeric(ds[2]-ds[1])
	    fileslist <- 1
     }
     for (dcount in fileslist) {
        if (ereefs_case == 1) {
	      input_file <- paste0(input_stem, format(as.Date(paste(year, month, dcount, sep="-")), '%Y-%m-%d'), '.nc')
        }
        #input_file <- paste0(input_file, '?', var_list, ',time,eta')
        nc <- ncdf4::nc_open(input_file)
	if (!is.na(eta_stem)) nc3 <- ncdf4::nc_open(etafile)
        if (ereefs_case > 0) {
          if (!is.null(nc$var[['t']])) {
            d <- as.Date(ncdf4::ncvar_get(nc, "t"), origin = as.Date("1990-01-01"))[from_day:(from_day+day_count-1)]
          } else {
            d <- as.Date(ncdf4::ncvar_get(nc, "time"), origin = as.Date("1990-01-01"))[from_day:(from_day+day_count-1)]
          }
        } else { 
           d <- ds[from_day:(from_day + day_count - 1)]
        }
        if (!is.null(nc$var[['eta']])) { 
          eta <- ncdf4::ncvar_get(nc, 'eta', start=c(location_grid[2], location_grid[1],from_day), count=c(1,1,day_count))
	} else {
          eta <- ncdf4::ncvar_get(nc3, 'eta', start=c(location_grid[2], location_grid[1],from_day), count=c(1,1,day_count))
	}
        im1 = i+1
        i <- i + length(d)
        ts_frame$date[im1:i] <- d
        z <- array(z_grid[2:length(z_grid)], dim=c(length(z_grid)-1, length(eta)))
        zm1 <- array(z_grid[1:(length(z_grid)-1)], dim=c(length(z_grid)-1, length(eta)))
        eta2 <- t(array(eta, dim=c(length(eta), length(z_grid)-1)))

        # Depth relative to surface
        zr <- z - eta2
        zm1r <- zm1 - eta2
        dz <- 0 * z
        wet <- (eta2 > zm1) & (z > botz)           # There is water in this layer
        bottom <- (wet & (zm1<botz))               # The bottom intersects this layer

        shallower_top <- (-zr < depth)               # The top of this layer is above the target depth
        deeper_bottom <- (-zm1r > depth) | (bottom)  # The bottom of this layer is below the target depth or this is the bottom layer
        target <- (shallower_top & deeper_bottom)
   
        for (j in 1:length(var_names)) {
          wc <- ncdf4::ncvar_get(nc, var_names[j], start=c(location_grid[2],location_grid[1],1,from_day), count=c(1,1,-1,day_count))

	  if (dim(target)[2] == 1) wc <- array(wc, dim=dim(target))
          ts_frame[im1:i, j+1] <- colSums(target * wc, na.rm=TRUE)
        }
        ncdf4::nc_close(nc)
        if (verbosity>0) setTxtProgressBar(pb,mcount)
      }
    }
    if (verbosity>0) close(pb)
    return(ts_frame)
}
