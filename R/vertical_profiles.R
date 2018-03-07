#' Extract vertical profiles from an array of specified latitudes and longitudes over a specified time-period from an eReefs or other EMS netcdf file.
#'
#' This is slow because we extract one profile at a time to avoid memory issues with large netcdf files.
#'
#' @return a list containing a vector of dates, an array of surface elevations (eta), the vertical grid (z_grid) and a data frame of values.
#' @param var_name A vector of EMS variable names. Defailts to c('Chl_a_sum', 'TN'))
#' @param location_latlon A data frame of latitudes and longitudes.  Defaults to data.frame(latitude=c(-20, -16), longitude=c(148.5, 152.0)).
#'                        If length(location_lat_lon)==2, extract every grid cell along a straight line between the two points specified.
#'                        Otherwise, extract only the locations corresponding to the cells nearest the specified points.
#' @param target_date Target date to extract profile. Can be a date, or text formatted for as.Date(), or a (year, month, day) vector.
#'                   Defaults to c(2016, 02, 04).
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
#'       only if eta is not in the file indicated by input_file (e.g. some GBR1 bgc files).
#' @export
get_ereefs_vertical_slice <- function(var_names=c('Chl_a_sum', 'TN'),
			 location_latlon=data.frame(latitude=c(-20, -16), longitude=c(148.5, 152.0)),
			 target_date = c(2016, 02, 04),
                         input_file = "http://dapds00.nci.org.au/thredds/dodsC/fx3/gbr4_bgc_GBR4_H2p0_B2p0_Chyd_Dcrt/gbr4_bgc_simple_2010-01.nc",
			 input_grid = NA,
			 eta_stem = NA)
{
  ereefs_case <- get_ereefs_case(input_file)
  input_stem <- get_file_stem(input_file)
  grids <- get_ereefs_grids(input_file, input_grid)
  x_grid <- grids[['x_grid']]
  y_grid <- grids[['y_grid']]
  z_grid <- grids[['z_grid']]
  nc <- ncdf4::nc_open(input_file)
  if (!is.null(nc$var[['latitude']])) {
    latitude <- ncvar_get(nc, 'latitude')
    longitude <- ncvar_get(nc, 'longitude')
  } else {
    latitude <- ncvar_get(nc, 'x_centre')
    longitude <- ncvar_get(nc, 'y_centre')
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
	lon1 <- location_latlon[1,2]
	lon2 <- location_latlon[2,2]
	lat1 <- location_latlon[1,1]
	lat2 <- location_latlon[2,1]

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


	location_latlon <- data.frame(latitude=latitude[which(intersected)], longitude=longitude[which(intersected)])
  }
  if (dim(location_latlon)[1]==0) stop("The line segment given does not intersect any model cells on this grid.")

  eta <- rep(NA, length(which(intersected)))
  values <- array(NA, c(length(z_grid)-1, length(eta), length(var_names)))
  mydata <- get_ereefs_profile(var_names=var_names, location_latlon=as.numeric(location_latlon[1,1:2]),
		     start_date = target_date, end_date = target_date, 
		     input_file = input_file, input_grid = input_grid, eta_stem = eta_stem)
  values[,1,] <- mydata$profiles
  eta[1] <- as.numeric(mydata$eta)
  grid_list <- mydata$grid_list

  for (i in 2:length(eta)) {
    print(paste('Extracting profile', i, 'of', length(eta)))
    mydata <- get_ereefs_profile(var_names=var_names, location_latlon=as.numeric(location_latlon[i,1:2]),
		     start_date = target_date, end_date = target_date, 
  		     input_file = input_file, input_grid = input_grid, eta_stem = eta_stem)
    values[,i,] <- mydata$profiles
    eta[i] <- as.numeric(mydata$eta)
  }

  return_list <- list(eta=eta, z_grid=z_grid, values=values)
}

#' Extract vertical profiles of specified variables  from a specified latitude and longitude over a specified time-period from an eReefs or other EMS netcdf file.
#'
#' @return a list containing a vector of dates, an array of surface elevations (eta), the vertical grid (z_grid) and a data frame of values.
#' @param var_name A vector of EMS variable names. Defailts to c('Chl_a_sum', 'TN'))
#' @param location_latlon Latitude and longitude of location to extract.  Defaults to c(-23.39189, 150.88852)
#' @param start_date Date from which to start extraction. Can be a date, or text formatted for as.Date(), or a (year, month, day) vector.
#'                   Defaults to c(2016, 02, 04).
#' @param end_date Date on which to end extraction, specified as for start_date. Defaults to c(2016, 03, 02).
#' @param input_file is the URI or file location of any of the EMS output files, 
#'        Defaults to "http://dapds00.nci.org.au/thredds/dodsC/fx3/gbr4_bgc_GBR4_H2p0_B2p0_Chyd_Dcrt/gbr4_bgc_simple_2010-01.nc". 
#'        If using Windows, you will need to set this to a local inputfile stem.
#' @param input_grid Either a list containing the coordinates of the cell corners (x_grid, y_grid and z_grid) or the name of the                                                        
#'      locally-stored or opendap-served netcdf file that contains these. If not specified, the function will first look for                                                            
#'      z_grid can be found in the first INPUT_STEM file, and if not found, will check whether the size                                                                                 
#'      of the variables in the input file corresponds to the size expected for GBR4 or GBR1, and load an appropriate                                                                   
#'      z grid from data files stored in this package. Alternatively, you can provide the location of a full                                                                            
#'      (not simple-format) ereefs netcdf output file such as                                                                                                                           
#'      "http://dapds00.nci.org.au/thredds/dodsC/fx3/gbr4_hydro_all/gbr4_all_2016-09.nc". 
#' @param eta_stem The URI or file location of the model output files that contains the surface elevation (eta), minus the
#'       date components of the filename in the case of GBR1 or GBR4 files, and ommitting the file extension, ".nc". Needed
#'       only if eta is not in the files indicated by input_stem (e.g. some GBR1 bgc files).
#' @param squeeze Whether to reduce the number of dimensions in the output profiles array if there is only one variable and/or 
#'       only one time-step. Default TRUE.
#' @export
get_ereefs_profile <- function(var_names=c('Chl_a_sum', 'TN'),
			 location_latlon=c(-23.39189, 150.88852),
			 start_date = c(2016, 02, 04),
			 end_date = c(2016, 03, 02),
                         input_file = "http://dapds00.nci.org.au/thredds/dodsC/fx3/gbr4_bgc_GBR4_H2p0_B2p0_Chyd_Dcrt/gbr4_bgc_simple_2010-01.nc",
			 input_grid = NA,
			 eta_stem = NA,
			 squeeze = TRUE)
{
  ereefs_case <- get_ereefs_case(input_file)
  input_stem <- get_file_stem(input_file)
  z_grid <- get_ereefs_grids(input_file, input_grid)[['z_grid']]
  nc <- ncdf4::nc_open(input_file)
  if (!is.null(nc$var[['latitude']])) {
    latitude <- ncvar_get(nc, 'latitude')
    longitude <- ncvar_get(nc, 'longitude')
  } else {
    latitude <- ncvar_get(nc, 'x_centre')
    longitude <- ncvar_get(nc, 'y_centre')
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

  var_list <- paste(var_names, collapse=",")

  if (ereefs_case == 4) {
      inputfile <- paste0(input_stem, format(as.Date(paste(start_year, start_month, 1, sep='-')), '%Y-%m'), 
			  '.nc')
      if (!is.na(eta_stem)) etafile  <- paste0(eta_stem, format(as.Date(paste(start_year, start_month, 1, sep='-')), '%Y-%m'), 
			  '.nc')
      blank_length <- end_date - start_date + 1
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
  nc <- ncdf4::nc_open(inputfile)
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
    tmp <- (latitude - location_latlon[1])^2 + (longitude - location_latlon[2])^2 
    tmp <- which.min(tmp) 
    location_grid <- c(floor((tmp+dim(latitude)[1]-1)/dim(latitude)[1]),
		       (tmp+dim(latitude)[1]-1)%%dim(latitude)[1] + 1)
  }
		     
  botz <- as.numeric(ncdf4::ncvar_get(nc, "botz", start=c(location_grid[2], location_grid[1]), count=c(1,1)))
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
        if (!is.na(eta_stem)) etafile <- paste0(eta_stem, format(as.Date(paste(year, month, 1, sep="-")), '%Y-%m'), '.nc')
     } else if (ereefs_case == 1) {
        fileslist <- from_day:(from_day+day_count-1)
        from_day <- 1
        day_count <- 1
     } else {
            from_day <- as.integer((as.Date(paste(year, month, from_day, sep="-")) - ds[1])/as.numeric(ds[2]-ds[1])) + 1
	    if (from_day<1) from_day <-1
	    fileslist <- 1
	    day_count <- day_count / as.numeric(ds[2]-ds[1])
     }

     for (dcount in fileslist) {
        if (ereefs_case == 1) {
	      inputfile <- paste0(input_stem, format(as.Date(paste(year, month, dcount, sep="-")), '%Y-%m-%d'), '.nc')
              if (!is.na(eta_stem)) etafile <- paste0(eta_stem, format(as.Date(paste(year, month, dcount, sep="-")), '%Y-%m-%d'), '.nc')
        }
        #inputfile <- paste0(inputfile, '?', var_list, ',time,eta')
        nc <- ncdf4::nc_open(inputfile)
        if (!is.na(eta_stem)) nc3 <- ncdf4::nc_open(etafile)
        if (ereefs_case > 0) {
          d <- ncdf4::ncvar_get(nc, "time", start=from_day, count=day_count)
          d <- as.Date(d, origin = as.Date("1990-01-01"))
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
	eta_record[im1:i] <- eta
        dates[im1:i] <- d
        z <- array(z_grid[2:length(z_grid)], dim=c(length(z_grid)-1, length(eta)))
        zm1 <- array(z_grid[1:(length(z_grid)-1)], dim=c(length(z_grid)-1, length(eta)))
        eta2 <- t(array(eta, dim=c(length(eta), length(z_grid)-1)))
        dz <- 0 * z
        wet <- (eta2 > zm1) & (z > botz)           # There is water in this layer
	dry <- !wet                                # There is no water in this layer

        for (j in 1:length(var_names)) {
          wc <- ncdf4::ncvar_get(nc, var_names[j], start=c(location_grid[2],location_grid[1],1,from_day), count=c(1,1,-1,day_count))
	  wc[dry] <- NA
	  if (dim(dz)[2] == 1) wc <- array(wc, dim=dim(dz))
          # take the depth-integrated average over the water column
          values[1:(length(z_grid)-1), j, im1:i] <- wc
        }
        ncdf4::nc_close(nc)
        if (!is.na(eta_stem)) ncdf4::nc_close(nc3)
        setTxtProgressBar(pb,mcount)
    }
  }
  close(pb)

  if (squeeze&(dim(values)[3] == 1)) values <- array(values, dim=dim(values)[c(1,2)])
  if (squeeze&(dim(values)[2] == 1)&(length(dim(values))==3)) values <- array(values, dim=dim(values)[c(1,3)])
  if (squeeze&(dim(values)[2] == 1)) values <- array(values, dim=dim(values)[1])
  return_list <- list(dates=dates, eta=eta_record, z_grid=z_grid, botz=botz, profiles=values)
  return(return_list)
}
