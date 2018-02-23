#' Extracts time series at a specified location from eReefs model output files
#'
#' Create a time-series of values of one or more selected model output variables in a specified layer of the
#' water column or sediment store (by default, the surface layer), aat a specified geographic location within 
#' the model domain.  See also get_ereefs_z_integrated_ts() to extract depth-integrated values and 
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
#' @param input_stem is the URI or file location of the eReefs output files, ommitting the year and month parts of the filenames and ommitting ".nc" . Defaults to "http://dapds00.nci.org.au/thredds/dodsC/fx3/gbr4_bgc_GBR4_H2p0_B2p0_Chyd_Dcrt/gbr4_bgc_simple_". If using Windows, you will need to set this to a local filename stem.
#' @export
#' @examples
#' get_ereefs_ts('Chl_a_sum', c(-23.39189, 150.88852), layer='surface')

get_ereefs_ts <- function(var_names=c('Chl_a_sum', 'TN'), 
                          location_latlon=c(-23.39189, 150.88852), 
                          layer='surface', 
		          start_date = c(2010,12,31), 
		          end_date = c(2016,10,31), 
                          input_stem = "http://dapds00.nci.org.au/thredds/dodsC/fx3/gbr4_bgc_GBR4_H2p0_B2p0_Chyd_Dcrt/gbr4_bgc_simple_")
{

  # Check whether this is a GBR1 or GBR4 ereefs file, or something else
  ereefs_case <- 0

  if (stringi::stri_detect_fixed(input_stem, 'gbr4')) {
	ereefs_case <- 4
  } else if (stringi::stri_detect_fixed(input_stem ,'gbr1')) {
	ereefs_case <- 1
  } else {
	stop("Not yet implemented for other than GBR4 GBR1 and files")
  }

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

  # Initialise the data frame with the right number of NAs
  blanks <- rep(NA, end_date - start_date + 1)
  ts_frame <- data.frame(as.Date(blanks), array(blanks, dim=c(length(blanks), length(var_names))))
  names(ts_frame) <- c("date", var_names)
  if (ereefs_case == 4) {
      inputfile <- paste0(input_stem, format(as.Date(paste(start_year, start_month, 1, sep='-')), '%Y-%m'), 
			  '.nc?latitude,longitude')
  } else {
      inputfile <- paste0(input_stem, format(as.Date(paste(start_year, start_month, start_day, sep='-')), '%Y-%m-%d'), 
			  '.nc?latitude,longitude')
  }

  # Find the nearest grid-points to the sampling location
  nc <- ncdf4::nc_open(inputfile)
  latitude <- ncdf4::ncvar_get(nc, "latitude")
  longitude <- ncdf4::ncvar_get(nc, "longitude")
  ncdf4::nc_close(nc)

  if (is.integer(location_latlon)) {
     location_grid <- location_latlon
  } else { 
    tmp <- (latitude - location_latlon[1])^2 + (longitude - location_latlon[2])^2 
    tmp <- which.min(tmp) 
    location_grid <- c( floor(tmp / dim(latitude)[1]) + 1, tmp%%dim(latitude)[1])
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
    if (mcount == (length(mths))) {
       day_count <- end_day
    } else if ((start_year==end_year)&&(start_month==end_month)) {
       day_count <- end_day - start_day
    } else {
       day_count <- daysIn(as.Date(paste(year, month, 1, sep='-')))
    }
    if (ereefs_case == 4) { 
       fileslist <- 1
       filename <- paste0(input_stem, format(as.Date(paste(year, month, 1, sep="-")), '%Y-%m'), '.nc')
    } else {
	    fileslist <- from_day:(from_day+day_count-1)
	    from_day <- 1
	    day_count <- 1
    }
    for (dcount in fileslist) {
      if (ereefs_case == 1) {
	    filename <- paste0(input_stem, format(as.Date(paste(year, month, dcount, sep="-")), '%Y-%m-%d'), '.nc')
      }
      inputfile <- paste0(filename, '?', var_list, ',time')
      nc <- ncdf4::nc_open(inputfile)
      d <- ncdf4::ncvar_get(nc, "time", start=from_day, count=day_count)
      im1 = i+1
      i <- i + length(d)
      for (j in 1:length(var_names)) {
          if (is.na(ndims[j])) {
             # We don't yet know the dimensions of the variable, so let's get them
              dims <- nc$var[[var_names[j]]][['size']]
              ndims[j] <- length(dims)
              if (layer == 'surface') {
               layer_actual[j] <- dims[3]
            } else {
               layer_actual[j] <- layer
            }
           }
          ts_frame$date[im1:i] <- as.Date(d, origin = as.Date("1990-01-01"))
          if (ndims[j] == 4) {
             fred <- ncdf4::ncvar_get(nc, var_names[j], start=c(location_grid[2],location_grid[1],layer_actual[j],from_day), count=c(1,1,1,day_count))
             ts_frame[im1:i, j+1] <- ncdf4::ncvar_get(nc, var_names[j], start=c(location_grid[2],location_grid[1],layer_actual[j],from_day), count=c(1,1,1,day_count))
          } else {
             fred <- ncdf4::ncvar_get(nc, var_names[j], start=c(location_grid[2],location_grid[1],layer_actual[j],from_day), count=c(1,1,1,day_count))
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
#' @param input_stem is the URI or file location of the eReefs output files, ommitting the year and month parts of the filenames and ommitting ".nc" . Defaults to "http://dapds00.nci.org.au/thredds/dodsC/fx3/gbr4_bgc_GBR4_H2p0_B2p0_Chyd_Dcrt/gbr4_bgc_simple_". If using Windows, you will need to set this to a local filename stem.
#' @param grid_file is the location of the netcdf file that contains the vertical grid information -- likely an output file
#'    from the hydrodynamic model used to drive the BGC model. Defaults to"http://dapds00.nci.org.au/thredds/dodsC/fx3/gbr4_hydro_all/gbr4_all_2016-09.nc"
#' @export
#' @examples
#' get_ereefs_z_integrated_ts('Chl_a_sum', c(-23.39189, 150.88852), layer='surface')
get_ereefs_z_integrated_ts <- function(var_names=c('Chl_a_sum', 'TN'), 
                          location_latlon=c(-23.39189, 150.88852), 
		          start_date = c(2010,12,31), 
		          end_date = c(2016,12,31), 
                          input_stem = "http://dapds00.nci.org.au/thredds/dodsC/fx3/gbr4_bgc_GBR4_H2p0_B2p0_Chyd_Dcrt/gbr4_bgc_simple_",
		          grid_file = "http://dapds00.nci.org.au/thredds/dodsC/fx3/gbr4_hydro_all/gbr4_all_2016-09.nc") 
{

# Get the model grid
inputfile <- paste0(input_stem, '2011-10.nc?latitude,longitude')
nc <- ncdf4::nc_open(inputfile)
latitude <- ncdf4::ncvar_get(nc, "latitude")
longitude <- ncdf4::ncvar_get(nc, "longitude")
ncdf4::nc_close(nc)

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
    stop('start_date must precede end_date')
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


# Find the nearest grid-points to the sampling location
if (is.integer(location_latlon)) {
    location_grid <- location_latlon
} else { 
    tmp <- (latitude - location_latlon[1])^2 + (longitude - location_latlon[2])^2 
    tmp <- which.min(tmp) 
    location_grid <- c( floor(tmp / dim(latitude)[1]) + 1, tmp%%dim(latitude)[1])
}

# Get the Z grid information
inputfile <- paste0(grid_file, '?z_grid,botz')
nc <- ncdf4::nc_open(inputfile)
zgrid <- ncdf4::ncvar_get(nc, "z_grid")
botz <- as.numeric(ncdf4::ncvar_get(nc, "botz", start=c(location_grid[2], location_grid[1]), count=c(1,1)))
ncdf4::nc_close(nc)

var_list <- paste(var_names, collapse=",")

# Initialise the data frame with the right number of NAs
ndims <- rep(NA, length(var_names))
blanks <- rep(NA, end_date - start_date + 1)
ts_frame <- data.frame(as.Date(blanks), array(blanks, dim=c(length(blanks), length(var_names))))
names(ts_frame) <- c("date", var_names)

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
   if (mcount == (length(mths))) {
      day_count <- end_day
   } else if ((start_year==end_year)&&(start_month==end_month)) {
      day_count <- end_day - start_day
   } else {
      day_count <- -1
   }
   inputfile <- paste0(input_stem, format(as.Date(paste(year, month, 1, sep='-')), '%Y-%m'), 
				'.nc?time,eta,', paste0(var_names, collapse=',')) 
   nc <- ncdf4::nc_open(inputfile)
   d <- ncdf4::ncvar_get(nc, "time", start=from_day, count=day_count)
   eta <- ncdf4::ncvar_get(nc, 'eta', start=c(location_grid[2], location_grid[1],from_day), count=c(1,1,day_count))
   im1 = i+1
   i <- i + length(d)
   ts_frame$date[im1:i] <- as.Date(d, origin = as.Date("1990-01-01"))
   z <- array(zgrid[2:length(zgrid)], dim=c(length(zgrid)-1, length(eta)))
   zm1 <- array(zgrid[1:(length(zgrid)-1)], dim=c(length(zgrid)-1, length(eta)))
   eta2 <- t(array(eta, dim=c(length(eta), length(zgrid)-1)))
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
     # take the depth-integrated average over the water column
     ts_frame[im1:i, j+1] <- colSums(dz * wc, na.rm=TRUE) / colSums(dz)
   }
   ncdf4::nc_close(nc)
   setTxtProgressBar(pb,mcount)
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
#' layer instead of a depth-integrated value and get_ereefs_z_integrated_ts() to calculate depth-integrated values.
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
#' @param input_stem is the URI or file location of the eReefs output files, ommitting the year and month parts of the filenames and ommitting ".nc" . Defaults to "http://dapds00.nci.org.au/thredds/dodsC/fx3/gbr4_bgc_GBR4_H2p0_B2p0_Chyd_Dcrt/gbr4_bgc_simple_". If using Windows, you will need to set this to a local filename stem.
#' @param grid_file is the location of the netcdf file that contains the vertical grid information -- likely an output file
#'    from the hydrodynamic model used to drive the BGC model. Defaults to"http://dapds00.nci.org.au/thredds/dodsC/fx3/gbr4_hydro_all/gbr4_all_2016-09.nc"
#' @export
#' @examples
#' get_ereefs_z_integrated_ts('Chl_a_sum', c(-23.39189, 150.88852), layer='surface')
get_ereefs_depth_specified_ts <- function(var_names=c('Chl_a_sum', 'TN'), 
                          location_latlon=c(-23.39189, 150.88852), 
                          depth = 1.0,
		          start_date = c(2010,12,31), 
		          end_date = c(2016,12,31), 
                          input_stem = "http://dapds00.nci.org.au/thredds/dodsC/fx3/gbr4_bgc_GBR4_H2p0_B2p0_Chyd_Dcrt/gbr4_bgc_simple_",
		          grid_file = "http://dapds00.nci.org.au/thredds/dodsC/fx3/gbr4_hydro_all/gbr4_all_2016-09.nc") 
{

# Get the model grid
inputfile <- paste0(input_stem, '2011-10.nc?latitude,longitude')
nc <- ncdf4::nc_open(inputfile)
latitude <- ncdf4::ncvar_get(nc, "latitude")
longitude <- ncdf4::ncvar_get(nc, "longitude")
ncdf4::nc_close(nc)

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
    stop('start_date must precede end_date')
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


# Find the nearest grid-points to the sampling location
if (is.integer(location_latlon)) {
    location_grid <- location_latlon
} else { 
    tmp <- (latitude - location_latlon[1])^2 + (longitude - location_latlon[2])^2 
    tmp <- which.min(tmp) 
    location_grid <- c( floor(tmp / dim(latitude)[1]) + 1, tmp%%dim(latitude)[1])
}

# Get the Z grid information
inputfile <- paste0(grid_file, '?z_grid,botz')
nc <- ncdf4::nc_open(inputfile)
zgrid <- ncdf4::ncvar_get(nc, "z_grid")
botz <- as.numeric(ncdf4::ncvar_get(nc, "botz", start=c(location_grid[2], location_grid[1]), count=c(1,1)))
ncdf4::nc_close(nc)

var_list <- paste(var_names, collapse=",")

# Initialise the data frame with the right number of NAs
ndims <- rep(NA, length(var_names))
blanks <- rep(NA, end_date - start_date + 1)
ts_frame <- data.frame(as.Date(blanks), array(blanks, dim=c(length(blanks), length(var_names))))
names(ts_frame) <- c("date", var_names)

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
   if (mcount == (length(mths))) {
      day_count <- end_day
   } else if ((start_year==end_year)&&(start_month==end_month)) {
      day_count <- end_day - start_day
   } else {
      day_count <- -1
   }
   inputfile <- paste0(input_stem, format(as.Date(paste(year, month, 1, sep='-')), '%Y-%m'), 
				'.nc?time,eta,', paste0(var_names, collapse=',')) 
   nc <- ncdf4::nc_open(inputfile)
   d <- ncdf4::ncvar_get(nc, "time", start=from_day, count=day_count)
   eta <- ncdf4::ncvar_get(nc, 'eta', start=c(location_grid[2], location_grid[1],from_day), count=c(1,1,day_count))
   im1 = i+1
   i <- i + length(d)
   ts_frame$date[im1:i] <- as.Date(d, origin = as.Date("1990-01-01"))

   # Absolute depth of layer boundaries
   z <- array(zgrid[2:length(zgrid)], dim=c(length(zgrid)-1, length(eta)))
   zm1 <- array(zgrid[1:(length(zgrid)-1)], dim=c(length(zgrid)-1, length(eta)))
   eta2 <- t(array(eta, dim=c(length(eta), length(zgrid)-1)))

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
     ts_frame[im1:i, j+1] <- colSums(target * wc, na.rm=TRUE)
   }

   ncdf4::nc_close(nc)
   setTxtProgressBar(pb,mcount)
}
close(pb)
return(ts_frame)
}
