#' Determines whether the file provided looks like an example of daily (e.g. GBR1), monthly (e.g. GBR4) or other (e.g. RECOM) output
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
		    stop('x_grid, y_grid and z_grid not found and grid size not recognised as GBR1 or GBR4. Please specifiy input_grid.')
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

#' Extracts time series at a specified location from eReefs model output files
#'
#' Create a time-series of values of one or more selected model output variables in a specified layer of the
#' water column or sediment store (by default, the surface layer), aat a specified geographic location within 
#' the model domain.  See also get_ereefs_depth_integrated_ts() to extract depth-integrated values and 
#' get_ereefs_depth_specified_ts() to extract values at a specified depth below the free (tidally moving) 
#' surface. Note that this function can use an OpenDAP URI if you are running it under Linux or MacOS, but 
#' not (as at February 2018) under Windows, because of issues with the windows version of the netcdf4 library. 
#' Under Windows, you can still use it for locally-saved files. Barbara Robson (AIMS).
#'
#'
#' @return a data frame containing the dates and values of extracted variables.
#' @param var_names either a single character value or a vector specifying the short names for variables that you 
#'        want from the netcdf file. Defaults to c('Chl_a_sum', 'TN').
#' @param location_latlon is a vector containing the decimal latitude and longitude of the desired location. If 
#'        you want to specify an x-y grid coordinate instead of a latitude and longitude, you can: to do this, 
#'        is.integer(location_latlon) must be TRUE. Defaults to c(-23.39189, 150.88852).
#' @param layer is the vertical grid layer to extract, or 'surface' to get the surface value. Defaults to 'surface'.
#' @param start_date date  for animation. Can either be a) a vector of integers in the format 
#'	c(year, month, day); b) a date obtained e.g. from as.Date(); or c) a character string 
#'      formatted for input to as.Date(). Defaults to c(2010,12,31).
#' @param end date for animation. Can either be a) a vector of integers in the format 
#'	c(year, month, day); b) a date obtained e.g. from as.Date(); or c) a character string 
#'      formatted for input to as.Date(). Defaults to c(2016,10,31).
#' @param input_file is the URI or file location of any of the EMS output files, 
#'        Defaults to "http://dapds00.nci.org.au/thredds/dodsC/fx3/gbr4_bgc_GBR4_H2p0_B2p0_Chyd_Dcrt/gbr4_bgc_simple_2010-01.nc". 
#'        If using Windows, you will need to set this to a local inputfile stem.
#' @export
#' @examples
#' get_ereefs_ts('Chl_a_sum', c(-23.39189, 150.88852), layer='surface')

get_ereefs_ts <- function(var_names=c('Chl_a_sum', 'TN'), 
                          location_latlon=c(-23.39189, 150.88852), 
                          layer='surface', 
		          start_date = c(2010,12,31), 
		          end_date = c(2016,10,31), 
                          input_file = "http://dapds00.nci.org.au/thredds/dodsC/fx3/gbr4_bgc_GBR4_H2p0_B2p0_Chyd_Dcrt/gbr4_bgc_simple_2010-01.nc")
{

  # Check whether output is daily (case 1), monthly (case 4) or something else (case 0)
  ereefs_case <- get_ereefs_case(input_file)
  input_stem <- get_file_stem(input_file)

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
      inputfile <- paste0(input_stem, format(as.Date(paste(start_year, start_month, 1, sep='-')), '%Y-%m'), 
			  '.nc')
	nc <- ncdf4::nc_open(inputfile)
	if (!is.null(nc$var[['t']])) { 
	    ds <- as.Date(ncdf4::ncvar_get(nc, "t"), origin = as.Date("1990-01-01"))
        } else {
	    ds <- as.Date(ncdf4::ncvar_get(nc, "time"), origin = as.Date("1990-01-01"))
	}
	ncdf4::nc_close(nc)
        blank_length <- as.numeric(end_date - start_date + 1) / as.numeric(ds[2] - ds[1])
			  # '.nc?latitude,longitude')
  } else if (ereefs_case == 1) {
      inputfile <- paste0(input_stem, format(as.Date(paste(start_year, start_month, start_day, sep='-')), '%Y-%m-%d'), 
			  '.nc')
      blank_length <- end_date - start_date + 1
			  # '.nc?latitude,longitude')
  } else {
      inputfile <- input_file
      nc <- ncdf4::nc_open(inputfile)
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


  if (is.integer(location_latlon)) {
     location_grid <- location_latlon
  } else { 
    # Find the nearest grid-points to the sampling location
    nc <- ncdf4::nc_open(inputfile)
    if (is.null(nc$var[['latitude']])) {
      latitude <- ncdf4::ncvar_get(nc, "y_centre")
      longitude <- ncdf4::ncvar_get(nc, "x_centre")
    } else { 
      latitude <- ncdf4::ncvar_get(nc, "latitude")
      longitude <- ncdf4::ncvar_get(nc, "longitude")
    }
    ncdf4::nc_close(nc)
    tmp <- (latitude - location_latlon[1])^2 + (longitude - location_latlon[2])^2 
    tmp <- which.min(tmp) 
    location_grid <- c(floor((tmp+dim(latitude)[1]-1)/dim(latitude)[1]),
		       (tmp+dim(latitude)[1]-1)%%dim(latitude)[1] + 1)
  }

  # Loop through monthly eReefs files to extract the data
  ndims <- rep(NA, length(var_names))
  layer_actual <- rep(NA, length(var_names))
  i <- 0
  mcount <- 0
  pb <- txtProgressBar(min = 0, max = length(mths), style = 3)
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
            inputfile <- paste0(input_stem, format(as.Date(paste(year, month, 1, sep="-")), '%Y-%m'), '.nc')
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
	    inputfile <- paste0(input_stem, format(as.Date(paste(year, month, dcount, sep="-")), '%Y-%m-%d'), '.nc')
      }
      #inputfile <- paste0(inputfile, '?', var_list, ',time')
      nc <- ncdf4::nc_open(inputfile)
      if (ereefs_case == 1) {
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
      ts_frame$date[im1:i] <- d
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
             ts_frame[im1:i, j+1] <- ncdf4::ncvar_get(nc, var_names[j], start=c(location_grid[2],location_grid[1],layer_actual[j],from_day), count=c(1,1,1,day_count))
          } else {
             ts_frame[im1:i, j+1] <- ncdf4::ncvar_get(nc, var_names[j], start=c(location_grid[2],location_grid[1],from_day), count=c(1,1,day_count))
          }
      }
      ncdf4::nc_close(nc)
      setTxtProgressBar(pb,mcount)
    }
  }
  close(pb)
  return(ts_frame)
}

#' Extracts depth-integrated time-series of selected variables at a specified location from eReefs output files
#'
#' Note that this function can use an OpenDAP URI if you are running it under Linux or MacOS, but 
#' not (as at February 2018) under Windows, because of issues with the windows version of the netcdf4 library. 
#' Under Windows, you can still use it for locally-saved files.
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
#'        Defaults to "http://dapds00.nci.org.au/thredds/dodsC/fx3/gbr4_bgc_GBR4_H2p0_B2p0_Chyd_Dcrt/gbr4_bgc_simple_2010-01.nc". 
#'        If using Windows, you will need to set this to a local inputfile stem.
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
#' get_ereefs_depth_integrated_ts('Chl_a_sum', c(-23.39189, 150.88852), layer='surface')
get_ereefs_depth_integrated_ts <- function(var_names=c('Chl_a_sum', 'TN'), 
                          location_latlon=c(-23.39189, 150.88852), 
		          start_date = c(2010,12,31), 
		          end_date = c(2016,12,31), 
                          input_file = "http://dapds00.nci.org.au/thredds/dodsC/fx3/gbr4_bgc_GBR4_H2p0_B2p0_Chyd_Dcrt/gbr4_bgc_simple_2010-01.nc",
			  input_grid = NA,
			  eta_stem = NA,
			  override_positive=FALSE)
{

  # Check whether this is a GBR1 or GBR4 ereefs file, or something else
  ereefs_case <- get_ereefs_case(input_file)
  input_stem <- get_file_stem(input_file)
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
      inputfile <- paste0(input_stem, format(as.Date(paste(start_year, start_month, 1, sep='-')), '%Y-%m'), 
			  '.nc')
      if (!is.na(eta_stem)) etafile  <- paste0(eta_stem, format(as.Date(paste(start_year, start_month, 1, sep='-')), '%Y-%m'), 
			  '.nc')
	nc <- ncdf4::nc_open(inputfile)
	if (!is.null(nc$var[['t']])) { 
	    ds <- as.Date(ncdf4::ncvar_get(nc, "t"), origin = as.Date("1990-01-01"))
        } else {
	    ds <- as.Date(ncdf4::ncvar_get(nc, "time"), origin = as.Date("1990-01-01"))
	}
	ncdf4::nc_close(nc)
        blank_length <- as.numeric(end_date - start_date + 1) / as.numeric(ds[2] - ds[1])
			  # '.nc?latitude,longitude')
  } else if (ereefs_case == 1) {
      inputfile <- paste0(input_stem, format(as.Date(paste(start_year, start_month, start_day, sep='-')), '%Y-%m-%d'), 
			  '.nc')
      if (!is.na(eta_stem)) etafile <- paste0(eta_stem, format(as.Date(paste(start_year, start_month, start_day, sep='-')), '%Y-%m-%d'), 
			  '.nc')
      blank_length <- end_date - start_date + 1
			  # '.nc?latitude,longitude')
  } else {
      inputfile <- input_file
      if (!is.na(eta_stem)) etafile <- paste0(eta_stem, '.nc')
      nc <- ncdf4::nc_open(inputfile)
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

  nc <- ncdf4::nc_open(inputfile)
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
  pb <- txtProgressBar(min = 0, max = length(mths), style = 3)
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
        inputfile <- paste0(input_stem, format(as.Date(paste(year, month, 1, sep="-")), '%Y-%m'), '.nc')
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
	      inputfile <- paste0(input_stem, format(as.Date(paste(year, month, dcount, sep="-")), '%Y-%m-%d'), '.nc')
              if (!is.na(eta_stem)) etafile <- paste0(eta_stem, format(as.Date(paste(year, month, dcount, sep="-")), '%Y-%m-%d'), '.nc')
        }
        #inputfile <- paste0(inputfile, '?', var_list, ',time,eta')
        nc <- ncdf4::nc_open(inputfile)
        if (!is.na(eta_stem)) nc3 <- ncdf4::nc_open(etafile)
        if (ereefs_case == 0) {
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
        setTxtProgressBar(pb,mcount)
    }
  }
  close(pb)
  return(ts_frame)
}

#' Extracts time-series of selected variables from eReefs output files at a specified location and depth below the
#' surface.
#'
#' Note that this function can use an OpenDAP URI if you are running it under Linux or MacOS, but 
#' not (as at February 2018) under Windows, because of issues with the windows version of the netcdf4 library. 
#' Under Windows, you can still use it for locally-saved files. See also get_ereefs_ts() to extract from a specified 
#' layer instead of a depth-integrated value and get_ereefs_depth_integrated_ts() to calculate depth-integrated values.
#' By Barbara Robson (AIMS).
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
#'        Defaults to "http://dapds00.nci.org.au/thredds/dodsC/fx3/gbr4_bgc_GBR4_H2p0_B2p0_Chyd_Dcrt/gbr4_bgc_simple_2010-01.nc". 
#'        If using Windows, you will need to set this to a local inputfile stem.
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
#' get_ereefs_depth_specified_ts('Chl_a_sum', c(-23.39189, 150.88852), depth=2.5)
get_ereefs_depth_specified_ts <- function(var_names=c('Chl_a_sum', 'TN'), 
                          location_latlon=c(-23.39189, 150.88852), 
                          depth = 1.0,
		          start_date = c(2010,12,31), 
		          end_date = c(2016,12,31), 
                          input_file = "http://dapds00.nci.org.au/thredds/dodsC/fx3/gbr4_bgc_GBR4_H2p0_B2p0_Chyd_Dcrt/gbr4_bgc_simple_2010-01.nc",
		          input_grid = NA,
			  eta_stem = NA)
{

  # Check whether this is a GBR1 or GBR4 ereefs file, or something else
  ereefs_case <- get_ereefs_case(input_file)
  input_stem <- get_file_stem(input_file)
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
      inputfile <- paste0(input_stem, format(as.Date(paste(start_year, start_month, 1, sep='-')), '%Y-%m'), 
			  '.nc')
      if (!is.na(eta_stem)) etafile  <- paste0(eta_stem, format(as.Date(paste(start_year, start_month, 1, sep='-')), '%Y-%m'), 
			  '.nc')
	nc <- ncdf4::nc_open(inputfile)
	if (!is.null(nc$var[['t']])) { 
	    ds <- as.Date(ncdf4::ncvar_get(nc, "t"), origin = as.Date("1990-01-01"))
        } else {
	    ds <- as.Date(ncdf4::ncvar_get(nc, "time"), origin = as.Date("1990-01-01"))
	}
	ncdf4::nc_close(nc)
        blank_length <- as.numeric(end_date - start_date + 1) / as.numeric(ds[2] - ds[1])
			  # '.nc?latitude,longitude')
  } else if (ereefs_case == 1) {
      inputfile <- paste0(input_stem, format(as.Date(paste(start_year, start_month, start_day, sep='-')), '%Y-%m-%d'), 
			  '.nc')
      if (!is.na(eta_stem)) etafile <- paste0(eta_stem, format(as.Date(paste(start_year, start_month, start_day, sep='-')), '%Y-%m-%d'), 
			  '.nc')
      blank_length <- end_date - start_date + 1
			  # '.nc?latitude,longitude')
  } else {
      inputfile <- input_file
      if (!is.na(eta_stem)) etafile <- paste0(eta_stem, '.nc')
      nc <- ncdf4::nc_open(inputfile)
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

  nc <- ncdf4::nc_open(inputfile)

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
  pb <- txtProgressBar(min = 0, max = length(mths), style = 3)
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
        inputfile <- paste0(input_stem, format(as.Date(paste(year, month, 1, sep="-")), '%Y-%m'), '.nc')
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
	      inputfile <- paste0(input_stem, format(as.Date(paste(year, month, dcount, sep="-")), '%Y-%m-%d'), '.nc')
        }
        #inputfile <- paste0(inputfile, '?', var_list, ',time,eta')
        nc <- ncdf4::nc_open(inputfile)
	if (!is.na(eta_stem)) nc3 <- ncdf4::nc_open(etafile)
        if (ereefs_case == 0) {
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
        setTxtProgressBar(pb,mcount)
      }
    }
    close(pb)
    return(ts_frame)
}
