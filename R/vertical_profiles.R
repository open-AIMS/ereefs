#' Extract vertical profiles from an array of specified latitudes and longitudes over a specified time-period from an eReefs or other EMS netcdf file.
#'
#' @return a list containing a vector of dates, an array of surface elevations (eta), the vertical grid (z_grid) and a data frame of values.
#' @param var_name A vector of EMS variable names. Defailts to c('Chl_a_sum', 'TN'))
#' @param location_latlon A data frame of latitudes and longitudes.  Defaults to data.frame(latitude=c(-20, -20), longitude=c(148.5, 149)).
#'                        If length(location_lat_lon)==2, extract every grid cell along a straight line between the two points specified.
#'                        Otherwise, extract only the locations corresponding to the cells nearest the specified points.
#' @param target_date Target date to extract profile. Can be a date, or text formatted for as.Date(), or a (year, month, day) vector.
#'                   Defaults to c(2016, 02, 04).
#' @param input_file is the URI or file location of any of the EMS output files, 
#'        Defaults to "http://dapds00.nci.org.au/thredds/dodsC/fx3/gbr4_bgc_GBR4_H2p0_B2p0_Chyd_Dcrt/gbr4_bgc_simple_2016-01.nc". 
#'        If using Windows, you will need to set this to a local input_file stem.
#' @param input_grid Name of the locally-stored or opendap-served netcdf file that contains the grid
#'      coordinates for the top and bottom of each layer (z_grid). If not specified, the function will first look for
#'      z_grid can be found in the first INPUT_STEM file, and if not found, will check whether the size 
#'      of the variables in the input file corresponds to the size expected for GBR4 or GBR1, and load an appropriate 
#'      z grid from data files stored in this package. Alternatively, you can provide the location of a full 
#'      (not simple-format) ereefs netcdf output file such as 
#'      "http://dapds00.nci.org.au/thredds/dodsC/fx3/gbr4_hydro_all/gbr4_all_2016-09.nc"
#' @param eta_stem The URI or file location of the model output files that contains the surface elevation (eta), or the stem of
#'       that filename minus the date components of the filename in the case of GBR1 or GBR4 files, and ommitting the file extension, ".nc". Needed
#'       only if eta is not in the file indicated by input_file (e.g. some GBR1 bgc files).
#' @param robust If TRUE, extract one profile at a time to avoid running out of memory. Robust but slow. Default FALSE.
#' @param override_positive Reverse the value of the "positive" attribute of botz for BGC files, assuming that it is
#'       incorrect. Default FALSE
#' @export
get_ereefs_slice <- function(var_names=c('Chl_a_sum', 'TN'),
			 location_latlon=data.frame(latitude=c(-20, -20), longitude=c(148.5, 149)),
			 target_date = c(2016, 02, 04),
                         input_file = "http://dapds00.nci.org.au/thredds/dodsC/fx3/gbr4_bgc_GBR4_H2p0_B2p0_Chyd_Dcrt/gbr4_bgc_simple_2016-01.nc",
			 input_grid = NA,
			 eta_stem = NA,
			 robust = FALSE,
			 override_positive = FALSE)
{
  ereefs_case <- get_ereefs_case(input_file)
  input_stem <- get_file_stem(input_file)
  if (!is.na(eta_stem)) {
	      if (stringi::stri_detect(eta_stem, fixed='.nc')) eta_stem <- get_file_stem(eta_stem) 
  }
  grids <- get_ereefs_grids(input_file, input_grid)
  x_grid <- grids[['x_grid']]
  y_grid <- grids[['y_grid']]
  z_grid <- grids[['z_grid']]
  nc <- ncdf4::nc_open(input_file)
  if (!is.null(nc$var[['latitude']])) {
    latitude <- ncdf4::ncvar_get(nc, 'latitude')
    longitude <- ncdf4::ncvar_get(nc, 'longitude')
  } else {
    latitude <- ncdf4::ncvar_get(nc, 'x_centre')
    longitude <- ncdf4::ncvar_get(nc, 'y_centre')
  }
  ncdf4::nc_close(nc)

  if (length(location_latlon)==2) {
	a <- (dim(x_grid) - 1)[1]
	b <- (dim(x_grid) - 1)[2]
	intersected <- array(FALSE, dim=c(a,b))
        gx <- c(x_grid[1:a, 1:b], x_grid[2:(a+1), 1:b], x_grid[2:(a+1), 2:(b+1)], x_grid[1:a, 2:(b+1)])
        gy <- c(y_grid[1:a, 1:b], y_grid[2:(a+1), 1:b], y_grid[2:(a+1), 2:(b+1)], y_grid[1:a, 2:(b+1)])
        gx <- array(gx, dim=c(a*b,4))
        gy <- array(gy, dim=c(a*b,4))
	# Find the grid cells intersected by the line between the two points.

	# First, define the line:
	lon1 <- location_latlon[1,'longitude']
	lon2 <- location_latlon[2,'longitude']
	lat1 <- location_latlon[1,'latitude']
	lat2 <- location_latlon[2,'latitude']

	# Define a line through the two points
	#Y = A*x+b ; A*x +b - Y = 0
	A <- (lat2 - lat1) / (lon2 - lon1)
	b = lat1 - (A * lon1)
	# For each edge, the line intersects the edge if the value of Ax+b-Y is +ve for one vertex and -ve for the other
	# Ignore exact hits on vertices.

	# Find grid-cells along this line
	c1 <- (A * gx[,1] + b - gy[,1]) > 0
	c2 <- (A * gx[,2] + b - gy[,2]) > 0
	c3 <- (A * gx[,3] + b - gy[,3]) > 0
	c4 <- (A * gx[,4] + b - gy[,4]) > 0

	intersected[c1!=c2] <- TRUE
	intersected[c2!=c3] <- TRUE
	intersected[c3!=c4] <- TRUE
	intersected[c4!=c1] <- TRUE

	# Exclude pesky NA values
	intersected[is.na(latitude)] <- FALSE
	intersected[is.na(c1+c2+c3+c4)] <- FALSE
	intersected[is.na(rowSums(gx)+rowSums(gy))] <- FALSE

	# Exclude grid-cells beyond the ends of the line segment
	minlat <- min(lat1, lat2)
	maxlat <- max(lat1, lat2)
	minlon <- min(lon1, lon2)
	maxlon <- max(lon1, lon2)
	intersected[((gx[,1]<minlon)&(gx[,2]<minlon)&(gx[,3]<minlon)&(gx[,4]<minlon))] <- FALSE
	intersected[((gx[,1]>maxlon)&(gx[,2]>maxlon)&(gx[,3]>maxlon)&(gx[,4]>maxlon))] <- FALSE
	intersected[((gy[,1]<minlat)&(gy[,2]<minlat)&(gy[,3]<minlat)&(gy[,4]<minlat))] <- FALSE
	intersected[((gy[,1]>maxlat)&(gy[,2]>maxlat)&(gy[,3]>maxlat)&(gy[,4]>maxlat))] <- FALSE

	llind <- which(intersected)
	location_latlon <- data.frame(latitude=latitude[llind], longitude=longitude[llind])
  }
  if (dim(location_latlon)[1]==0) stop("The line segment given does not intersect any model cells on this grid.")
  i <- dim(location_latlon)[1]
  location_latlon[seq(i,1,-1),] <- location_latlon

  eta <- rep(NA, length(which(intersected)))
  botz <- rep(NA, length(which(intersected)))
  values <- array(NA, c(length(z_grid)-1, length(eta), length(var_names)))

  if (robust) {
    mydata <- get_ereefs_profile(var_names=var_names, location_latlon=as.numeric(location_latlon[1,c('latitude','longitude')]),
		     start_date = target_date, end_date = target_date, 
		     input_file = input_file, input_grid = input_grid, eta_stem = eta_stem)
    values[,1,] <- mydata$profiles
    eta[1] <- as.numeric(mydata$eta)
    grid_list <- mydata$grid_list

    for (i in 2:length(eta)) {
      print(paste('Extracting profile', i, 'of', length(eta)))
      mydata <- get_ereefs_profile(var_names=var_names, location_latlon=as.numeric(location_latlon[i,c('latitude','longitude')]),
		     start_date = target_date, end_date = target_date, 
  		     input_file = input_file, input_grid = input_grid, eta_stem = eta_stem)
      values[,i,] <- mydata$profiles
      eta[i] <- as.numeric(mydata$eta)
      botz[i] <- as.numeric(mydata$botz)
    }
  } else {
    # Date to plot
    if (is.vector(target_date)) {
      target_date <- as.Date(paste(target_date[1], target_date[2], target_date[3], sep='-'))
    } else if (is.character(target_date)) {
      target_date <- as.Date(target_date)
    }
    target_day <- as.integer(format(target_date, '%d'))
    target_month <- as.integer(format(target_date, '%m'))
    target_year <- as.integer(format(target_date, '%Y'))

    #var_list <- paste(var_names, collapse=",")

    location_grid <- cbind(floor((llind+dim(latitude)[1]-1)/dim(latitude)[1]),
		       (llind+dim(latitude)[1]-1)%%dim(latitude)[1] + 1)

    if (ereefs_case == 4) {
        input_file <- paste0(input_stem, format(as.Date(paste(target_year, target_month, 1, sep='-')), '%Y-%m'), 
			  '.nc')
        if (!is.na(eta_stem)) etafile  <- paste0(eta_stem, format(as.Date(paste(target_year, target_month, 1, sep='-')), '%Y-%m'), 
			  '.nc')
        nc <- ncdf4::nc_open(input_file)
        if (!is.null(nc$var[['t']])) {
          ds <- as.Date(ncdf4::ncvar_get(nc, "t"), origin = as.Date("1990-01-01"))
        } else {
          ds <- as.Date(ncdf4::ncvar_get(nc, "time"), origin = as.Date("1990-01-01"))
        }
	di <- which.min(abs(ds - target_date))
        ncdf4::nc_close(nc)
    } else if (ereefs_case == 1) {
        input_file <- paste0(input_stem, format(as.Date(paste(target_year, target_month, target_day, sep='-')), '%Y-%m-%d'), 
			  '.nc')
        if (!is.na(eta_stem)) etafile <- paste0(eta_stem, format(as.Date(paste(target_year, target_month, target_day, sep='-')), '%Y-%m-%d'), 
			  '.nc')
	di <- 1
    } else {
        input_file <- input_file
        if (!is.na(eta_stem)) etafile <- paste0(eta_stem, '.nc')
        nc <- ncdf4::nc_open(input_file)
        if (!is.null(nc$var[['t']])) {
          ds <- as.Date(ncdf4::ncvar_get(nc, "t"), origin = as.Date("1990-01-01"))
        } else {
          ds <- as.Date(ncdf4::ncvar_get(nc, "time"), origin = as.Date("1990-01-01"))
        }
	di <- which.min(abs(ds - target_date))
        ncdf4::nc_close(nc)
    }
    nc <- ncdf4::nc_open(input_file)
    if (!is.na(eta_stem)) {
	    nc3 <- ncdf4::nc_open(etafile)
    } else {
	    nc3 <- nc
    }

    startv <- c(min(location_grid[,2]), min(location_grid[, 1]))
    countv <- c(max(location_grid[,2]), max(location_grid[, 1])) - startv + 1
    if ((countv[1] == 1)&(countv[2] == 1)) stop('Slice matches a single grid-cell; use get_ereefs_profile() instead.')
    location_grid <- t(t(location_grid) - c(startv[2], startv[1])) + 1
    location_grid <- cbind(location_grid[,2], location_grid[,1])
    numlines <- dim(location_grid)[1]
    numcols <- 3
    if (countv[2] == 1) {
	    location_grid <- location_grid[,1]
	    numlines <- length(location_grid)
	    numcols <- 2
    } else if (countv[1] == 1) {
	    location_grid <- location_grid[,2]
	    numlines <- length(location_grid)
	    numcols <- 2
    }
    ind3d <- rep(cbind(location_grid, NA), each=(length(z_grid)-1))
    ind3d <- array(ind3d, c(numlines*(length(z_grid)-1), numcols))
    ind3d[, numcols]  <- rep(1:(length(z_grid)-1), numlines)

    values <- array(NA, dim=c(length(z_grid)-1, dim(location_latlon)[1], length(var_names)))

    zat <- ncdf4::ncatt_get(nc, "botz")
    if (!is.null(zat$positive)) {
	  if (zat$positive=="down") zsign <- -1 else zsign <- 1
          if (override_positive) zsign <- -zsign
    } else {
	  zsign <-1
    }
    botz <- zsign * as.numeric(ncdf4::ncvar_get(nc, "botz", start=startv, count=countv)[location_grid])
    eta <- as.numeric(ncdf4::ncvar_get(nc3, "eta", start=c(startv, di), count=c(countv, 1))[location_grid])

    z <- array(z_grid[2:length(z_grid)], dim=c(length(z_grid)-1, length(eta)))
    zm1 <- array(z_grid[1:(length(z_grid)-1)], dim=c(length(z_grid)-1, length(eta)))
    eta2 <- t(array(eta, dim=c(length(eta), length(z_grid)-1)))
    botz2 <- t(array(botz, dim=c(length(botz), length(z_grid)-1)))
    wet <- (eta2 > zm1) & (z > botz2)          # There is water in this layer
    dry <- !wet

    for (j in 1:length(var_names)) { 
      wc <- ncdf4::ncvar_get(nc, var_names[j], start=c(startv, 1, di), count=c(countv, -1, 1))
      wc <- ncdf4::ncvar_get(nc, var_names[j], start=c(startv, 1, di), count=c(countv, -1, 1))[ind3d]
      wc[dry] <- NA
      if (dim(z)[2] == 1) wc <- array(wc, dim=dim(z))
      values[1:(length(z_grid)-1), , j] <- wc
    }
    ncdf4::nc_close(nc)
    if (length(eta_stem)>1) ncdf4::nc_close(nc3)

  }

  dimnames(values)[[1]] <- 1:(length(z_grid)-1)
  dimnames(values)[[3]] <- var_names

  return_list <- list(eta=eta, botz=botz, z_grid=z_grid, values=values, latlon=location_latlon)
}

#' Extract vertical profiles of specified variables  from a specified latitude and longitude over a specified time-period from an eReefs or other EMS netcdf file.
#'
#' See also plot_ereefs_profile(), which relies on output from this function.
#'
#' Note that this function assumes consistent frequency of model output, even if the time-series extends across multiple output files (e.g.
#' multiple months of eReefs output). If the data in eta_stem is output on a different interval from the data in input_file, the function
#' will do its best, but surface elevation estimates may not exactly match the time-stamps in the main input file.
#'
#' @return a list containing a vector of dates, an array of surface elevations (eta), the vertical grid (z_grid) and a data frame of values.
#' @param var_name A vector of EMS variable names. Defailts to c('Chl_a_sum', 'TN'))
#' @param location_latlon Latitude and longitude of location to extract.  Defaults to c(-23.39189, 150.88852)
#' @param start_date Date from which to start extraction. Can be a date, or text formatted for as.Date(), or a (year, month, day) vector.
#'                   Defaults to c(2016, 02, 04).
#' @param end_date Date on which to end extraction, specified as for start_date. Defaults to c(2016, 03, 02).
#' @param input_file is the URI or file location of any of the EMS output files, 
#'        Defaults to "http://dapds00.nci.org.au/thredds/dodsC/fx3/gbr4_bgc_GBR4_H2p0_B2p0_Chyd_Dcrt/gbr4_bgc_simple_2016-01.nc". 
#'        If using Windows, you will need to set this to a local input_file stem.
#' @param input_grid Either a list containing the coordinates of the cell corners (x_grid, y_grid and z_grid) or the name of the                                                        
#'      locally-stored or opendap-served netcdf file that contains these. If not specified, the function will first look for                                                            
#'      z_grid can be found in the first INPUT_STEM file, and if not found, will check whether the size                                                                                 
#'      of the variables in the input file corresponds to the size expected for GBR4 or GBR1, and load an appropriate                                                                   
#'      z grid from data files stored in this package. Alternatively, you can provide the location of a full                                                                            
#'      (not simple-format) ereefs netcdf output file such as                                                                                                                           
#'      "http://dapds00.nci.org.au/thredds/dodsC/fx3/gbr4_hydro_all/gbr4_all_2016-09.nc". 
#' @param eta_stem The URI or file location of the model output files that contains the surface elevation (eta), or the stem of that
#'       filename minus the
#'       date components of the filename in the case of GBR1 or GBR4 files, and ommitting the file extension, ".nc". Needed
#'       only if eta is not in the files indicated by input_stem (e.g. some GBR1 bgc files).
#' @param squeeze Whether to reduce the number of dimensions in the output profiles array if there is only one variable and/or 
#'       only one time-step. Default TRUE.
#' @export
get_ereefs_profile <- function(var_names=c('Chl_a_sum', 'TN'),
			 location_latlon=c(-23.39189, 150.88852),
			 start_date = c(2016, 02, 04),
			 end_date = c(2016, 03, 02),
                         input_file = "http://dapds00.nci.org.au/thredds/dodsC/fx3/gbr4_bgc_GBR4_H2p0_B2p0_Chyd_Dcrt/gbr4_bgc_simple_2016-01.nc",
			 input_grid = NA,
			 eta_stem = NA,
			 squeeze = TRUE,
			 override_positive=FALSE)
{
  ereefs_case <- get_ereefs_case(input_file)
  input_stem <- get_file_stem(input_file)
  if (!is.na(eta_stem)) {
	      if (stringi::stri_detect(eta_stem, fixed='.nc')) eta_stem <- get_file_stem(eta_stem) 
  }
  z_grid <- get_ereefs_grids(input_file, input_grid)[['z_grid']]
  nc <- ncdf4::nc_open(input_file)
  if (!is.null(nc$var[['latitude']])) {
    latitude <- ncdf4::ncvar_get(nc, 'latitude')
    longitude <- ncdf4::ncvar_get(nc, 'longitude')
  } else {
    latitude <- ncdf4::ncvar_get(nc, 'x_centre')
    longitude <- ncdf4::ncvar_get(nc, 'y_centre')
  }
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

  #var_list <- paste(var_names, collapse=",")

  if (ereefs_case == 4) {
      input_file <- paste0(input_stem, format(as.Date(paste(start_year, start_month, 1, sep='-')), '%Y-%m'), 
			  '.nc')
      if (!is.na(eta_stem)) etafile  <- paste0(eta_stem, format(as.Date(paste(start_year, start_month, 1, sep='-')), '%Y-%m'), 
			  '.nc')
  } else if (ereefs_case == 1) {
      input_file <- paste0(input_stem, format(as.Date(paste(start_year, start_month, start_day, sep='-')), '%Y-%m-%d'), 
			  '.nc')
      if (!is.na(eta_stem)) etafile <- paste0(eta_stem, format(as.Date(paste(start_year, start_month, start_day, sep='-')), '%Y-%m-%d'), 
			  '.nc')
  } else {
      input_file <- input_file
      if (!is.na(eta_stem)) etafile <- paste0(eta_stem, '.nc')
  }
  nc <- ncdf4::nc_open(input_file)
  if (!is.null(nc$var[['t']])) { 
    ds <- as.Date(ncdf4::ncvar_get(nc, "t"), origin = as.Date("1990-01-01"))
  } else {
    ds <- as.Date(ncdf4::ncvar_get(nc, "time"), origin = as.Date("1990-01-01"))
  }
  if (!is.na(eta_stem)) {
    nc3 <- ncdf4::nc_open(etafile)
    if (!is.null(nc3$var[['t']])) { 
      eta_ds <- as.Date(ncdf4::ncvar_get(nc3, "t"), origin = as.Date("1990-01-01"))
    } else {
      eta_ds <- as.Date(ncdf4::ncvar_get(nc3, "time"), origin = as.Date("1990-01-01"))
    }
  } else {
    eta_ds <- ds
  }
  if (length(ds)==1) {
    blank_length <- as.numeric(end_date - start_date + 1)
  } else {
    blank_length <- as.numeric(end_date - start_date + 1) / as.numeric(ds[2] - ds[1])
  }
  if (!is.na(eta_stem)) nc3 <- ncdf4::nc_open(etafile)
  # Initialise the data frame with the right number of NAs
  blanks <- rep(NA, blank_length)
  dates <- as.Date(blanks)
  values <- array(blanks, dim=c(length(z_grid)-1, length(var_names), length(blanks)))
  colnames(values) <- var_names
  eta_record <- blanks

  if (is.integer(location_latlon)) {
     location_grid <- location_latlon
  } else { 
    # Find the nearest grid-points to the sampling location
    grid_ind <- (latitude - location_latlon[1])^2 + (longitude - location_latlon[2])^2 
    grid_ind <- which.min(grid_ind) 
    location_grid <- c(floor((grid_ind+dim(latitude)[1]-1)/dim(latitude)[1]),
		       (grid_ind+dim(latitude)[1]-1)%%dim(latitude)[1] + 1)
  }
		     
  zat <- ncdf4::ncatt_get(nc, "botz")
  if (!is.null(zat$positive)) {
	  if (zat$positive=="down") zsign <- -1 else zsign <- 1
          if (override_positive) zsign <- -zsign
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
       from_record <- which.min(ds - (start_date+0.499)) # Choose a record as close to midday as we can
       eta_from_record <- which.min(eta_ds - (start_date+0.499)) # Choose a record as close to midday as we can
     } else {
       from_record <- 1
       eta_from_record <- 1
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
        if (!is.na(eta_stem)) etafile <- paste0(eta_stem, format(as.Date(paste(year, month, 1, sep="-")), '%Y-%m'), '.nc')
     } else if (ereefs_case == 1) {
        fileslist <- 1:day_count
       day_count <- 1
     } else {
	fileslist <- 1
     }
     if (start_date == end_date) { # Assume we only want a single profile
       day_count <- 1
     } else if (length(ds)>1) {
       day_count <- day_count / as.numeric(ds[2]-ds[1])
     }
     eta_day_count <- day_count / as.numeric(eta_ds[2]-eta_ds[1])

     for (dcount in fileslist) {
        if (ereefs_case == 1) {
	      input_file <- paste0(input_stem, format(as.Date(paste(year, month, dcount, sep="-")), '%Y-%m-%d'), '.nc')
              if (!is.na(eta_stem)) etafile <- paste0(eta_stem, format(as.Date(paste(year, month, dcount, sep="-")), '%Y-%m-%d'), '.nc')
        }
        #input_file <- paste0(input_file, '?', var_list, ',time,eta')
        nc <- ncdf4::nc_open(input_file)
        if (!is.na(eta_stem)) nc3 <- ncdf4::nc_open(etafile)
        if (ereefs_case == 0) {
          d <- ncdf4::ncvar_get(nc, "time", start=from_record, count=day_count)
          d <- as.Date(d, origin = as.Date("1990-01-01"))
        } else {
	  d <- ds[from_record:(from_record + day_count - 1)]
        }
        if (!is.null(nc$var[['eta']])) { 
          eta <- ncdf4::ncvar_get(nc, 'eta', start=c(location_grid[2], location_grid[1],from_record), count=c(1,1,day_count))
	} else {
	  if (is.na(eta_stem)) stop('eta not found in netcdf file. Please specify eta_stem.')
          eta <- ncdf4::ncvar_get(nc3, 'eta', start=c(location_grid[2], location_grid[1],eta_from_record), count=c(1,1,eta_day_count))
	}
        im1 = i+1
        i <- i + length(d)
	if (length(eta_ds)==length(ds)) {
	  eta_record[im1:i] <- eta
	} else {
  	  if (length(eta_ds)<length(ds)) {
	    if (dcount==1) warning(paste('Surface elevation (eta) in', etafile, 'is output less frequently than', var_names[1], 'in', input_file,
			  '. Assuming eta always==0, though this is unlikely'))
	    eta_record[im1:i] <- 0*c(im:i)
	  } else {
	    if (length(ds)==1) {
	      ind <- which.min(abs(eta_ds-ds[1]))
	      eta_record[im1:i] <- eta[ind]
	      if ((dcount==1)&((eta_ds[ind]-ds[1])>(1/48))) warning(paste('Surface elevation (eta) in', etafile, 'is output more frequently than', var_names[1], 'in', input_file,
			  '. Using eta from closest times to output of', var_names[1], ', which is more than 30 mins away from desired time.'))
	    } else {
	      interval <- as.numeric(ds[2] - ds[1])/as.numeric(eta_ds[2] - eta_ds[1])
	      ind <- seq(from=which.min(abs(eta_ds-ds[1])), to=length(eta), by=interval)
	      if (dcount==1) {
		      tgap <- abs(eta_ds[ind] - ds)
		      if (max(tgap)>(1/48)) warning(paste('Surface elevation (eta) in', etafile, 'is output more frequently than', var_names[1], 'in', input_file,
			  '. Using eta from closest times to output of', var_names[1], ', which is sometimes more than 30 mins away from desired time.'))
	      }
	      eta_record[im1:i] <- eta[ind]
	    }
	  }
	}
        dates[im1:i] <- d
        z <- array(z_grid[2:length(z_grid)], dim=c(length(z_grid)-1, length(eta_record)))
        zm1 <- array(z_grid[1:(length(z_grid)-1)], dim=c(length(z_grid)-1, length(eta_record)))
        eta2 <- t(array(eta, dim=c(length(eta_record), length(z_grid)-1)))
        wet <- (eta2 > zm1) & (z > botz)           # There is water in this layer
	dry <- !wet                                # There is no water in this layer

        for (j in 1:length(var_names)) {
          wc <- ncdf4::ncvar_get(nc, var_names[j], start=c(location_grid[2],location_grid[1],1,from_record), count=c(1,1,-1,day_count))
	  wc[dry] <- NA
	  if (dim(z)[2] == 1) wc <- array(wc, dim=dim(z))
          values[1:(length(z_grid)-1), j, im1:i] <- wc
        }
        ncdf4::nc_close(nc)
        if (length(eta_stem)>1) ncdf4::nc_close(nc3)
        setTxtProgressBar(pb,mcount)
    }
  }
  close(pb)

  if (squeeze&(dim(values)[3] == 1)) {                          # Only one time-step
	  values <- array(values, dim=dim(values)[c(1,2)])
	  colnames(values) <- var_names
  }
  if (squeeze&(dim(values)[2] == 1)&(length(dim(values))==3)) { # Only one variable, but multiple time-steps
	  values <- array(values, dim=dim(values)[c(1,3)])
  }
  if (all(is.na(values))) warning('No wet cells in this profile. Either this is a land cell or the positive attribute of botz is incorrect (use overrid_positive=TRUE) if this is the case)')
  return_list <- list(dates=dates, eta=eta_record, z_grid=z_grid, botz=botz, profiles=values)
  return(return_list)
}

#' Plots a single vertical profile using output from get_ereefs_profile()
#'
#' Relies on output from get_ereefs_profile().
#'
#' @param profileObj A list object as output by get_ereefs_profiles(), containing dates, eta, z_grid, botz and profiles
#' @param var_name The name of the variable to plot (must be a colname in profile$profiles). Default 'Chl_a_sum'.
#' @param target_date The target date (plot the profile closest in time to this).
#' @param p The handle of an existing figure, if you don't want to create a new figure
#' @return the handle of a figure containing the vertical profile plot
#' @examples plot_ereefs_profile(get_ereefs_profile('TN'))
#' @export
plot_ereefs_profile <- function(profileObj, var_name='Chl_a_sum', target_date=c(2016,01,01), p=NA, colour='blue') {
  # Date to plot
  if (is.vector(target_date)) {
	  target_date <- as.Date(paste(target_date[1], target_date[2], target_date[3], sep='-'))
  } else if (is.character(target_date)) {
	  target_date <- as.Date(target_date)
  }
  day <- which.min(abs(target_date-profileObj$dates))
  colnum <- which(colnames(profileObj$profiles)==var_name)
  if (length(dim(profileObj$profiles))>2) {
	  dind <- which.min(abs(profileObj$dates-target_date))
          values <- array(profileObj$profiles[, colnum, dind], length(profileObj$z_grid)-1)
	  eta <- profileObj$eta[dind]
  } else {
	  values <- array(profileObj$profiles[, colnum])
	  eta <- profileObj$eta
  }
  wet <- which(!is.na(values))
  values <- c(values[wet], values[max(wet)])
  z <- c(profileObj$botz, profileObj$z_grid[wet[1:length(wet)-1]+1], eta)
  mydata <- data.frame(z=z, values=values)
  if (length(p)==1) p <- ggplot2::ggplot(mydata)
  p <- p + ggplot2::geom_line(data=mydata, ggplot2::aes(x=values, y=z), colour=colour) + xlab(var_name) + ylab('metres above msl')
  print(p)
  return(p)
}

#' Produces a coloured tile plot of a vertical slice already fetched from an eReefs or other EMS netcdf file.
#'
#' Relies on output from get_ereefs_slice().
#'
#' @param slice A list object as output by get_ereefs_slice(), containing dates, eta, z_grid, botz,
#'              a data frame of values and a data frame of latitudes and longitudes
#' @param scale_col Colours to use for low and high values. Default c("ivory", "hotpink").
#' @param scale_lim values for low and high limits of colourscale. Defaults to full range.
#' @return p handle for the generated figure
#' @export
plot_ereefs_slice <- function(slice, var_name='Chl_a_sum', scale_col=c("ivory", "hotpink"), scale_lim=NA) {
	numprofiles <- dim(slice$values)[2]
	layers <- length(slice$z_grid) - 1
	zmin <- array(slice$z_grid[1:layers], c(layers, numprofiles))
	zmax <- array(slice$z_grid[2:(layers+1)], c(layers, numprofiles))
	d <- rep(NA, numprofiles)
	for (i in 1:numprofiles) {
		zmin[zmin[,i]<slice$botz[i],i] <- slice$botz[i]
		zmin[zmin[,i]>slice$eta[i],i] <- slice$eta[i]
		zmax[zmax[,i]<slice$botz[i],i] <- slice$botz[i]
		zmax[zmax[,i]>slice$eta[i],i] <- slice$eta[i]
		if (i==1) {
			d[i] <- 0
		} else {
	          d[i] <- earth.dist(slice$latlon[i-1,'longitude'],slice$latlon[i-1,'latitude'], slice$latlon[i,'longitude'], slice$latlon[i,'latitude']) + d[i-1]
		}
	}
	dmin <- c(-d[2]/2, d[1:(length(d)-1)])
	dmin <- t(array(dmin, c(numprofiles, layers)))
	dmax <- c(d[2:length(d)], d[length(d)] + (d[length(d)] - d[length(d)-1])/2)
	dmax <- t(array(dmax, c(numprofiles, layers)))

	ind <- which(!is.na(slice$values[, , var_name]))
	values <- slice$values[,, var_name]
	if (length(scale_lim)==1) {
		scale_lim[1] <- min(c(values[ind]))
		scale_lim[2] <- max(c(values[ind]))
	}

	mydata <- data.frame(xmin=dmin[ind], xmax=dmax[ind], ymin=zmin[ind], ymax=zmax[ind], z=values[ind])
	p <- ggplot2::ggplot(data=mydata,ggplot2:: aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax, fill=z)) + 
		ggplot2::geom_rect() +
		ggplot2::scale_fill_gradient(name=var_name, low=scale_col[1], high=scale_col[2], limits=scale_lim, oob=scales::squish) +
		ggplot2::ylab('metres above msl') +
		ggplot2::xlab('kilometres from start of transect')
	plot(p)
	return(p)
}

#' Calculate rough distance in kilometers between two points
#'
#' Not exported. This is very approximate - a package is available if a more accurate distance is needed.
earth.dist <- function (long1, lat1, long2, lat2)
{
rad <- pi/180
a1 <- lat1 * rad
a2 <- long1 * rad
b1 <- lat2 * rad
b2 <- long2 * rad
dlon <- b2 - a2
dlat <- b1 - a1
a <- (sin(dlat/2))^2 + cos(a1) * cos(b1) * (sin(dlon/2))^2
c <- 2 * atan2(sqrt(a), sqrt(1 - a))
R <- 6378.145
d <- R * c
return(d)
}

#' Produces a coloured rect plot of a vertical profile over time
#'
#' Relies on output from get_ereefs_profile().
#'
#' @param profileObj A list object as output by get_ereefs_profile(), containing dates, eta, z_grid, botz,
#'              and a data frame of values.
#' @param scale_col Colours to use for low and high values. Default c("ivory", "hotpink").
#' @param scale_lim values for low and high limits of colourscale. Defaults to full range.
#' @return p handle for the generated figure
#' @export
plot_ereefs_zvt <- function(slice, var_name='Chl_a_sum', scale_col=c("ivory", "hotpink"), scale_lim=NA) {
	numprofiles <- dim(slice$profiles)[3]
	layers <- length(slice$z_grid) - 1
	zmin <- array(slice$z_grid[1:layers], c(layers, numprofiles))
	zmax <- array(slice$z_grid[2:(layers+1)], c(layers, numprofiles))
	for (i in 1:numprofiles) {
		zmin[zmin[,i]<slice$botz,i] <- slice$botz
		zmin[zmin[,i]>slice$eta[i],i] <- slice$eta[i]
		zmax[zmax[,i]<slice$botz,i] <- slice$botz
		zmax[zmax[,i]>slice$eta[i],i] <- slice$eta[i]
	}
	d <- slice$dates
	dmin <- c(d[1]-(d[2]-d[1])/2, d[1:(length(d)-1)])
	dmin <- t(array(dmin, c(numprofiles, layers)))
	dmax <- c(d[2:length(d)], d[length(d)-1] + (d[length(d)] - d[length(d)-1])/2)
	dmax <- t(array(dmax, c(numprofiles, layers)))

	ind <- which(!is.na(slice$profiles[, var_name, ]))
	profiles <- slice$profiles[, var_name, ]
	if (length(scale_lim)==1) {
		scale_lim[1] <- min(c(profiles[ind]))
		scale_lim[2] <- max(c(profiles[ind]))
	}

	mydata <- data.frame(xmin=as.Date(dmin[ind], origin='1990-01-01'), xmax=as.Date(dmax[ind], origin='1990-01-01'), 
			     ymin=zmin[ind], ymax=zmax[ind], 
			     z=profiles[ind])
	p <- ggplot2::ggplot(data=mydata, aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax, fill=z)) + 
		ggplot2::geom_rect() +
		ggplot2::scale_x_date() +
		ggplot2::ylab('metres above msl') +
		ggplot2::scale_fill_gradient(name=var_name, low=scale_col[1], high=scale_col[2], limits=scale_lim, oob=scales::squish)
	plot(p)
	return(p)
}

#' Calculate rough distance in kilometers between two points
#'
#' Not exported. This is very approximate - a package is available if a more accurate distance is needed.
earth.dist <- function (long1, lat1, long2, lat2)
{
rad <- pi/180
a1 <- lat1 * rad
a2 <- long1 * rad
b1 <- lat2 * rad
b2 <- long2 * rad
dlon <- b2 - a2
dlat <- b1 - a1
a <- (sin(dlat/2))^2 + cos(a1) * cos(b1) * (sin(dlon/2))^2
c <- 2 * atan2(sqrt(a), sqrt(1 - a))
R <- 6378.145
d <- R * c
return(d)
}
