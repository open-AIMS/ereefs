#' Extract a vertical slice or depth-resolved transect from an eReefs or other EMS netcdf output file.
#'
#' Extracts either a series of vertical profiles at points corresponding to a dataFrame of specified latitudes and longitudes (e.g. a cruise
#' transect) or a vertical slice along a line between two specified points, at a specified point in time.
#' Instead of interpolating, this function takes the values of cells intersected by the line segment, which can cause some striping when lines
#' transect the grid diagonally.
#'
#' @param var_name A vector of EMS variable names. Defailts to c('Chl_a_sum', 'TN'))
#' @param location_latlon A data frame of latitudes and longitudes.  Defaults to data.frame(latitude=c(-20, -20), longitude=c(148.5, 149)).
#'                        If length(location_lat_lon)==2, extract every grid cell along a straight line between the two points specified.
#'                        Otherwise, extract only the locations corresponding to the cells nearest the specified points.
#' @param target_date Target date to extract profile. Can be a date, or text formatted for as.Date(), or a (year, month, day) vector.
#'                   Defaults to c(2016, 02, 04). If target_date is a vector, 0.499 is added to the calculated date to bring it as
#'                   close to midday as possible.
#' @param input_file is the URL or file location of any of the EMS output files or a THREDDS catalog URI. 
#'        Defaults to a menu selection based on current NCI catalogs. Can also be set to "nci", "menu" or "catalog" for the same behaviour.
#'        Set to "old_menu" to provide old menu options instead of menu options from the NCI catalog.
#'        Numeric values are interpreted as references to selections available from the old menu.
#'        Short codes can be used for some options (codenames as used in https://research.csiro.au/ereefs/models/model-outputs/access-to-raw-model-output/ )
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
#' @param robust If TRUE, extract one profile at a time to avoid running out of memory. Relatively robust but slow. Default FALSE.
#' @param override_positive Reverse the value of the "positive" attribute of botz for BGC files, assuming that it is
#'       incorrect. Default FALSE
#' @return a list containing an array of surface elevations (eta), bottom depths (botz), the locations of the cell centres of the intersected cells
#'       (cell_centres), the locations of the cell edge points where the line intersects the grid edges (cell_intersections), and the vertical grid 
#'       (z_grid) and a data frame of tracer values from the model (values) and 'crossref', which indicates the starting cell in the returned values
#'       array associated with each input location.
#' @export
get_ereefs_slice <- function(var_names=c('Chl_a_sum', 'TN'),
			 location_latlon=data.frame(latitude=c(-20, -20), longitude=c(148.5, 149)),
			 target_date = c(2016, 02, 04),
                         input_file = "menu",
			 input_grid = NA,
			 eta_stem = NA,
			 robust = FALSE,
			 override_positive = FALSE)
{
  input_file <- substitute_filename(input_file)
  ereefs_case <- get_ereefs_case(input_file)
  input_stem <- get_file_stem(input_file)
  #check_platform_ok(input_stem)

  # Remember which direction our slice should be facing
  if ((location_latlon[1,1]<location_latlon[2,1]) |
      ((location_latlon[1,1]==location_latlon[2,1])&(location_latlon[1,2]<location_latlon[2,2]))) {
     latlon_dir <- 1
  } else {
     latlon_dir <- -1
  }

  if (!is.na(eta_stem)) {
	      if (stringi::stri_endswith(eta_stem, fixed='.nc')) eta_stem <- get_file_stem(eta_stem) 
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

  location_ll <- data.frame(latitude=NULL, longitude=NULL)
  location_edges <- data.frame(latitude=NULL, longitude=NULL)
  intersected <- NULL
  llind <- NULL
  crossref <- NULL
  #browser()
  for (i in 1:(dim(location_latlon)[1]-1)) {
     print(paste('transect section', i, 'of', dim(location_latlon)[1]-1))
     first_point <- (i==1)
     li <- find_intersections(location_latlon[i:(i+1),], x_grid, y_grid, latitude, longitude, first_point)
     #print(li[[1]])
     #print(li[[2]])
	   if (first_point) {
       location_ll <- li[[1]]
	     location_edges <- li[[2]]
	     intersected <- li[[3]]
	     llind <- li[[4]]
     } else if (!is.null(li[[2]])) {
	     location_edges <- rbind(location_edges, li[[2]])
	     intersected <- intersected | li[[3]] # This is boolean per cell so we want to combine them
       if (li[[4]][length(li[[4]])] == llind[length(llind)]) {
         li[[4]] <- li[[4]][1:(length(li[[4]])-1)]
         li[[1]] <- li[[1]][1:(length(li[[1]])-1), ]
       } else if (li[[4]][1] == llind[length(llind)]) {
         li[[4]] <- li[[4]][2:length(li[[4]])]
         li[[1]] <- li[[1]][2:length(li[[1]]), ]
       }
       location_ll <- rbind(location_ll, li[[1]])
	     llind <- c(llind, li[[4]])
      }
      crossref[i] <- length(llind)
  }
  if (llind[1] == llind[2]) {
    # This will happen when first_point returns a single value but subsequent points also return values
    llind <- llind[2:length(llind)]
    crossref <- crossref[2:length(crossref)]
    location_ll <- location_ll[2:dim(location_ll)[1],]
  }
  #print(llind)
  #print(location_ll)
  #browser()
  #dum1 <- location_ll; dum1$key <- 1:length(location_ll$latitude)
  #dum1 <- distinct(dum1, latitude, longitude, .keep_all=TRUE)
  #location_ll <- location_ll[dum1$key,]
  #llind <- llind[dum1$key]
  #dum1 <- location_edges; dum1$key <- 1:length(location_edges$latitude)
  #dum1 <- distinct(dum1, latitude, longitude, .keep_all=TRUE)
  #location_edges <- location_edges[dum1$key,]
  crossref <- c(1, crossref)
  location_latlon <- location_ll

  if (dim(location_latlon)[1]==0) stop("The line segment given does not intersect any model cells on this grid.")
  i <- dim(location_latlon)[1]

  eta <- rep(NA, dim(location_ll)[1])
  botz <- rep(NA, dim(location_ll)[1])
  values <- array(NA, c(length(z_grid)-1, length(eta), length(var_names)))

  if (robust) {
    #print(target_date)
    #browser()
    mydata <- get_ereefs_profile(var_names=var_names, location_latlon=as.numeric(location_latlon[1,c('latitude','longitude')]),
		     start_date = target_date, end_date = target_date, 
		     input_file = input_file, input_grid = input_grid, eta_stem = eta_stem, override_positive=override_positive)
    values[,1,] <- mydata$profiles
    eta[1] <- as.numeric(mydata$eta)
    botz[1] <- as.numeric(mydata$botz)
    grid_list <- mydata$grid_list

    for (i in 2:length(eta)) {
      print(paste('Extracting profile', i, 'of', length(eta)))
      mydata <- get_ereefs_profile(var_names=var_names, location_latlon=as.numeric(location_latlon[i,c('latitude','longitude')]),
		     start_date = target_date, end_date = target_date, 
  		     input_file = input_file, input_grid = input_grid, eta_stem = eta_stem, override_positive=override_positive)
      values[,i,] <- mydata$profiles
      eta[i] <- as.numeric(mydata$eta)
      botz[i] <- as.numeric(mydata$botz)
    }
  } else {
    if (is.vector(target_date)) {
      target_date <- as.Date(paste(target_date[1], target_date[2], target_date[3], sep='-')) + 0.499
    } else if (is.character(target_date)) {
      target_date <- as.Date(target_date)
    }
    target_day <- as.integer(format(target_date, '%d'))
    target_month <- as.integer(format(target_date, '%m'))
    target_year <- as.integer(format(target_date, '%Y'))

    #var_list <- paste(var_names, collapse=",")

    #location_grid <- arrayInd(llind, dim(latitude))
    location_grid <- cbind(floor((llind+dim(latitude)[1]-1)/dim(latitude)[1]),
		       (llind+dim(latitude)[1]-1)%%dim(latitude)[1] + 1)

    if (ereefs_case[2] == '4km') {
        input_file <- paste0(input_stem, format(as.Date(paste(target_year, target_month, 1, sep='-')), '%Y-%m'), 
			  '.nc')
        if (!is.na(eta_stem)) etafile  <- paste0(eta_stem, format(as.Date(paste(target_year, target_month, 1, sep='-')), '%Y-%m'), 
			  '.nc')
    } else if (ereefs_case[2] == '1km') {
        input_file <- paste0(input_stem, format(as.Date(paste(target_year, target_month, target_day, sep='-')), '%Y-%m-%d'), 
			  '.nc')
        if (!is.na(eta_stem)) etafile <- paste0(eta_stem, format(as.Date(paste(target_year, target_month, target_day, sep='-')), '%Y-%m-%d'), 
			  '.nc')
    } else {
        input_file <- input_file
        if (!is.na(eta_stem)) etafile <- paste0(eta_stem, '.nc')

    } 
    nc <- ncdf4::nc_open(input_file)
	  if (!is.null(nc$var[['t']])) { 
      ds <- as.Date(safe_ncvar_get(nc, "t"), origin = as.Date("1990-01-01"))
    } else { 
      ds <- as.Date(safe_ncvar_get(nc, "time"), origin = as.Date("1990-01-01"))
    }
    if (target_date < ds[1]) {
      warning(paste('Target date (', target_date, ') is before start of data. Setting target_date to', ds[1]))
      target_date <- ds[1]
    }
    if (target_date > ds[length(ds)]) {
      warning(paste('Target date (', target_date, ') is after end of data. Setting target_date to', ds[length(ds)]))
      target_date <- ds[length(ds)]
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
    from_record <- which.min(abs(ds - (target_date)))
    eta_from_record <- which.min(abs(eta_ds - (target_date)))
    ########

    # Find the outer grid coordinates of the area that we need to extract from netcdf files to encompass the slice
    startv <- c(min(location_grid[,2]), min(location_grid[, 1]))
    countv <- c(max(location_grid[,2]), max(location_grid[, 1])) - startv + 1
    if ((countv[1] == 1)&(countv[2] == 1)) stop('Slice matches a single grid-cell; use get_ereefs_profile() instead.')

    # Adjust grid locations so that they are relative to the region to be extracted instead of the whole model domain
    location_grid <- t(t(location_grid) - c(startv[2], startv[1])) + 1
    location_grid <- cbind(location_grid[,2], location_grid[,1])
    numlines <- dim(location_grid)[1]
    numcols <- 3

    # check whether all points are within a single model grid row or column, and adjust indices accordingly
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
    if (!is.null(nc$var[['eta']])) { 
       eta <- ncdf4::ncvar_get(nc, 'eta', start=c(startv,from_record), count=c(countv,1))[location_grid]
    } else { 
       if (is.na(eta_stem)) stop('eta not found in netcdf file. Please specify eta_stem.')
       eta <- as.numeric(ncdf4::ncvar_get(nc3, "eta", start=c(startv, eta_from_record), count=c(countv, 1)))
       eta <- as.numeric(ncdf4::ncvar_get(nc3, "eta", start=c(startv, eta_from_record), count=c(countv, 1))[location_grid])
    }

    z <- array(z_grid[2:length(z_grid)], dim=c(length(z_grid)-1, length(eta)))
    zm1 <- array(z_grid[1:(length(z_grid)-1)], dim=c(length(z_grid)-1, length(eta)))
    eta2 <- t(array(eta, dim=c(length(eta), length(z_grid)-1)))
    botz2 <- t(array(botz, dim=c(length(botz), length(z_grid)-1)))
    wet <- (eta2 > zm1) & (z > botz2)          # There is water in this layer
    dry <- !wet

    for (j in 1:length(var_names)) { 
      #wc <- ncdf4::ncvar_get(nc, var_names[j], start=c(startv, 1, from_record), count=c(countv, -1, 1))
      wc <- ncdf4::ncvar_get(nc, var_names[j], start=c(startv, 1, from_record), count=c(countv, -1, 1))[ind3d]
      wc[dry] <- NA
      if (dim(z)[2] == 1) wc <- array(wc, dim=dim(z))
      values[1:(length(z_grid)-1), , j] <- wc
    }
    ncdf4::nc_close(nc)
    if (length(eta_stem)>1) ncdf4::nc_close(nc3)

  }

  dimnames(values)[[1]] <- 1:(length(z_grid)-1)
  dimnames(values)[[3]] <- var_names

  if ((location_latlon[1,1]<location_latlon[length(location_latlon$latitude),1]) |
      ((location_latlon[1,1]==location_latlon[length(location_latlon$latitude),1])&(location_latlon[1,2]<location_latlon[length(location_latlon$latitude),2]))) {
     directn <- 1
  } else {
     directn <- -1
  }
  if (latlon_dir == directn) {
     ind <- seq(length(eta), 1, -1)
     eta <- eta[ind]
     botz <- botz[ind]
     location_latlon <- location_latlon[seq(dim(location_latlon)[1], 1, -1),]
     location_edges <- location_edges[seq(dim(location_edges)[1], 1, -1),]
     values[,1:length(eta),] <- values[,ind,]
  }

  return_list <- list(eta=eta, botz=botz, z_grid=z_grid, values=values, cell_centres=location_latlon, cell_intersections=location_edges, crossref = crossref)
}

#' Extract vertical profiles of specified variables  from a specified latitude and longitude over a specified time-period from an eReefs or other EMS netcdf file.
#'
#' See also plot_ereefs_profile(), which relies on output from this function.
#'
#' Note that this function assumes consistent frequency of model output, even if the time-series extends across multiple output files (e.g.
#' multiple months of eReefs output). If the data in eta_stem is output on a different interval from the data in input_file, the function
#' will do its best, but surface elevation estimates may not exactly match the time-stamps in the main input file.
#'
#' Known bugs: currently only works for one time-step at a time for 4km hydrodynamic model output. As a work-around, either call
#' once for each time-step or exctract hydrodynamic variables from the bgc model output.
#'
#' @return a list containing a vector of dates, an array of surface elevations (eta), the vertical grid (z_grid) and a data frame of values.
#' @param var_name A vector of EMS variable names. Defailts to c('Chl_a_sum', 'TN'))
#' @param location_latlon Latitude and longitude of location to extract.  Defaults to c(-23.39189, 150.88852)
#' @param start_date Date from which to start extraction. Can be a date, or text formatted for as.Date(), or a (year, month, day) vector.
#'                   Defaults to c(2016, 02, 04). If start_date is a vector, 0.499 is added to the calculated date to bring the start 
#'                   as close to midday as possible.
#' @param end_date Date on which to end extraction, specified as for start_date. Defaults to c(2016, 03, 02).
#' @param input_file is the URL or file location of any of the EMS output files or a THREDDS catalog URI. 
#'        Defaults to a menu selection based on current NCI catalogs. Can also be set to "nci", "menu" or "catalog" for the same behaviour.
#'        Set to "old_menu" to provide old menu options instead of menu options from the NCI catalog.
#'        Numeric values are interpreted as references to selections available from the old menu.
#'        Short codes can be used for some options (codenames as used in https://research.csiro.au/ereefs/models/model-outputs/access-to-raw-model-output/ )
#' @param input_grid Either a list containing the coordinates of the cell corners (x_grid, y_grid and z_grid) or the name of the                                                        
#'      locally-stored or opendap-served netcdf file that contains these. If not specified, the function will first look for                                                            
#'      z_grid can be found in the first INPUT_STEM file, and if not found, will check whether the size                                                                                 
#'      of the variables in the input file corresponds to the size expected for GBR4 or GBR1, and load an appropriate                                                                   
#'      z grid from data files stored in this package. Alternatively, you can provide the location of a full                                                                            
#'      (not simple-format) ereefs netcdf output file such as                                                                                                                           
#'      "http://dapds00.nci.org.au/thredds/dodsC/fx3/gbr4_hydro_all/gbr4_all_2016-09.nc". 
#' @param eta_stem The URI or file location of the model output files that contains the surface elevation (eta), or the stem of that
#'       filename minus the date components of the filename in the case of GBR1 or GBR4 files, and ommitting the file extension, ".nc". 
#'       Needed only if eta is not in the files indicated by input_stem (e.g. some GBR1 bgc files).
#'       Assumes that the eta files contain data on the same time-step as the input files.
#' @param squeeze Whether to reduce the number of dimensions in the output profiles array if there is only one variable and/or 
#'       only one time-step. Default TRUE.
#' @export
get_ereefs_profile <- function(var_names=c('Chl_a_sum', 'TN'),
			 location_latlon=c(-23.39189, 150.88852),
			 start_date = c(2016, 02, 04),
			 end_date = c(2016, 03, 02),
       input_file = "catalog",
			 input_grid = NA,
			 eta_stem = NA,
			 squeeze = TRUE,
			 override_positive=FALSE)
{

  # Get parameter values and assign results from returned list to relevant variable names
  # This assigns input_file, ereefs_case, input_stem, start_date, end_date, start_tod, start_month, start_year,
  # end_date, end_day, end_month, end_year, mths, years, var_list, ereefs_origin and blank_length
  assignList(get_params(start_date, end_date, input_file, var_names))

  if (length(dim(location_latlon)) > 0) stop('get_ereefs_profile() only handles one location per call. Use get_ereefs_slice() instead if appropriate.')
  #check_platform_ok(input_stem)
  if (!is.na(eta_stem)) {
	      if (stringi::stri_endswith(eta_stem, fixed='.nc')) eta_stem <- get_file_stem(eta_stem) 
  }
  z_grid <- get_ereefs_grids(input_file, input_grid)[['z_grid']]

  if (!is.na(eta_stem)) {
    if (ereefs_case[2] == '4km') {
      etafile  <- paste0(eta_stem, format(as.Date(paste(start_year, start_month, 1, sep='-')), '%Y-%m'), 
			  '.nc')
    } else if (ereefs_case[2] == '1km') {
      etafile <- paste0(eta_stem, format(as.Date(paste(start_year, start_month, start_day, sep='-')), '%Y-%m-%d'), 
			  '.nc')
    } else {
      etafile <- paste0(eta_stem, '.nc')
    }
  }

  nc <- ncdf4::nc_open(input_file)
  if (!is.na(eta_stem)) {
    nc3 <- ncdf4::nc_open(etafile)
    eta_ds <- get_origin_and_times(etafile)
  } else {
    eta_ds <- ds
  }
  if ((length(ds)!=1) & ((ds[2]-ds[1])<(10/86400))) { 
    # if it's less than 10 seconds, assume there is a duplicate output as in input_file=12
    # In this case, blank_length has already been set by get_params()
    blank_length <- as.numeric(end_date - start_date + 1) / as.numeric(ds[2] - ds[1])
  }
  # Initialise the data frame with the right number of NAs
  blanks <- rep(NA, blank_length)
  dates <- rep(get_origin_and_times(input_file)[[2]][1], blank_length)
  values <- array(blanks, dim=c(length(z_grid)-1, length(var_names), length(blanks)))
  colnames(values) <- var_names
  eta_record <- blanks

  if (is.integer(location_latlon)) {
     location_grid <- location_latlon
  } else { 
    # Find the nearest grid-points to the sampling location
    grid_ind <- (latitude - location_latlon[1])^2 + (longitude - location_latlon[2])^2 
    grid_ind <- which.min(grid_ind) 
    #location_grid <- arrayInd(grid_ind, dim(latitude))
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
       from_record <- which.min(abs(ds - (start_date)))
       eta_from_record <- which.min(abs(eta_ds - (start_date)))
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
     if (ereefs_case[2] == '4km') { 
        fileslist <- 1
        input_file <- paste0(input_stem, format(as.Date(paste(year, month, 1, sep="-")), '%Y-%m'), '.nc')
        if (!is.na(eta_stem)) etafile <- paste0(eta_stem, format(as.Date(paste(year, month, 1, sep="-")), '%Y-%m'), '.nc')
     } else if (ereefs_case[2] == '1km') {
        fileslist <- 1:day_count
       day_count <- 1
     } else {
	    fileslist <- 1
     }

     if ((start_date == end_date)|(length(ds)==1)) { # Assume we only want a single profile
       day_count <- 1
       eta_day_count <- 1
     } else if (length(ds)>1) {
       day_count <- day_count / as.numeric(ds[2]-ds[1])
       eta_day_count <- day_count / as.numeric(eta_ds[2]-eta_ds[1])
     }

     for (dcount in fileslist) {
        if (ereefs_case[2] == '1km') {
	      input_file <- paste0(input_stem, format(start_date+dcount-1, sep="-", '%Y-%m-%d'), '.nc')
              if (!is.na(eta_stem)) etafile <- paste0(eta_stem, format(start_date+dcount-1, sep="-", '%Y-%m-%d'), '.nc')
        }
        #input_file <- paste0(input_file, '?', var_list, ',time,eta')
        d <- get_origin_and_times(input_file)[[2]][from_record:(from_record+day_count-1)]
        nc <- ncdf4::nc_open(input_file)
        if (!is.na(eta_stem)) nc3 <- ncdf4::nc_open(etafile)
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
             if (start_date==end_date) {
                eta_record[im1:i] <- eta
                ind <- 1
             } else { 
                ind <- which.min(abs(eta_ds-ds[1])) 
             }
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
      z <- array(z_grid[2:length(z_grid)], dim=c(length(z_grid)-1, length(d)))
      zm1 <- array(z_grid[1:(length(z_grid)-1)], dim=c(length(z_grid)-1, length(d)))
      eta2 <- t(array(eta, dim=c(length(d), length(z_grid)-1)))
      wet <- (eta2 > zm1) & (z > botz)           # There is water in this layer
	    dry <- !wet                                # There is no water in this layer

      for (j in 1:length(var_names)) {
        wc <- ncdf4::ncvar_get(nc, var_names[j], start=c(location_grid[2],location_grid[1],1,from_record), count=c(1,1,-1,day_count))
        wc[dry] <- NA
        if (dim(z)[2] == 1) wc <- array(wc, dim=dim(z))
         values[1:(length(z_grid)-1), j, im1:(im1+day_count-1)] <- wc 
       }
       ncdf4::nc_close(nc)
       if (length(eta_stem)>1) ncdf4::nc_close(nc3)
          setTxtProgressBar(pb,mcount)
   } # end for dcount
  } #end for month
  close(pb)

  if (squeeze&(dim(values)[3] == 1)) {                          # Only one time-step
	  values <- array(values, dim=dim(values)[c(1,2)])
	  colnames(values) <- var_names
  }
  if (squeeze&(dim(values)[2] == 1)&(length(dim(values))==3)) { # Only one variable, but multiple time-steps
	  values <- array(values, dim=dim(values)[c(1,3)])
     if (all(is.na(values))) warning('No wet cells in this profile. Either the location is a land cell or the positive attribute of botz is incorrect (use override_positive=TRUE) if this is the case)')
     return_list <- list(dates=dates, eta=eta_record, z_grid=z_grid, botz=botz, profiles=values)
     names(return_list) <- c(names(return_list)[1:4], var_names)
  } else { 
     if (all(is.na(values))) warning('No wet cells in this profile. Either the location is a land cell or the positive attribute of botz is incorrect (use override_positive=TRUE) if this is the case)')
     return_list <- list(dates=dates, eta=eta_record, z_grid=z_grid, botz=botz, profiles=values)
  }
  return(return_list)
}

#' Plots a single vertical profile using output from get_ereefs_profile()
#'
#' Relies on output from get_ereefs_profile().
#'
#' @param profileObj A list object as output by get_ereefs_profiles(), containing dates, eta, z_grid, botz and profiles
#' @param var_name The name of the variable to plot (must be a colname in profile$profiles). Default 'Chl_a_sum'.
#'        If profileObj contains only one variable, var_name is ignored and the content of the profile is shown.
#' @param target_date The target date (plot the profile closest in time to this).
#' @param p The handle of an existing figure, if you don't want to create a new figure
#' @return the handle of a figure containing the vertical profile plot
#' @examples 
#' \dontrun{
#' plot_ereefs_profile(get_ereefs_profile('TN'))
#' }
#' @export
plot_ereefs_profile <- function(profileObj, var_name='Chl_a_sum', target_date=c(2016,01,01), p=NA, colour='blue') {
  # Date to plot
  # Dates to plot
  target_date <- get_chron_date(target_date)
  target_day <- as.integer(chron::days(target_date))
  target_tod <- as.numeric(target_date) - as.integer(target_date)
  target_month <- as.integer(chron:::months.default(target_date))
  target_year <- as.integer(as.character(chron::years(target_date)))

  day <- which.min(abs(target_date-profileObj$dates))
  if (names(profileObj)[5]=="profiles") { 
     colnum <- which(colnames(profileObj$profiles)==var_name)
     if (length(dim(profileObj$profiles))>2) {
	     dind <- which.min(abs(profileObj$dates-target_date)) 
        values <- array(profileObj$profiles[, colnum, dind], length(profileObj$z_grid)-1)
	     eta <- profileObj$eta[dind]
     } else {
	     values <- array(profileObj$profiles[, colnum])
	     eta <- profileObj$eta
     }
  } else { 
     # Only one variable
     var_name <- names(profileObj)[5]
     names(profileObj)[5] <- "profiles"
     if (is.null(dim(profileObj$profiles))) {
        # Only one time-step
        values <- profileObj$profiles
	     eta <- profileObj$eta
     } else {
	     dind <- which.min(abs(profileObj$dates-target_date)) 
        values <- array(profileObj$profiles[, dind], length(profileObj$z_grid)-1)
	     eta <- profileObj$eta[dind]
     }
  }
  wet <- which(!is.na(values))
  values <- c(values[wet], values[max(wet)])
  z <- c(profileObj$botz, profileObj$z_grid[wet[1:length(wet)-1]+1], eta)
  mydata <- data.frame(z=z, values=values)
  if (length(p)==1) p <- ggplot2::ggplot(mydata)
  p <- p + ggplot2::geom_path(data=mydata, ggplot2::aes(x=values, y=z), colour=colour) + ggplot2::xlab(var_name) + ggplot2::ylab('metres above msl')
  #print(p)
  return(p)
}

#' Produces a coloured tile plot of a vertical slice already fetched from an eReefs or other EMS netcdf file.
#'
#' Relies on output from get_ereefs_slice().
#'
#' @param slice A list object as output by get_ereefs_slice(), containing dates, eta, z_grid, botz,
#'              a data frame of values and a data frame of latitudes and longitudes
#' @param scale_col Vector of colours to use for low and high values in the colour scale. This can be a colour 
#'      from the ggplot colour palette or a RGB hash code, or "spectral". Ignored for true_colour plots. 
#'      Example: c('ivory', 'coral4').
#'      If one value is given (other than "spectral"), low colour is set to ivory and high colour to the value given.
#'      If three values are given, uses scale_fill_gradient2 (spectrum from low to high through middle value).
#'      Defaults to "spectral"
#' @param scale_lim values for low and high limits of colourscale. Defaults to full range.
#' @return p handle for the generated figure
#' @export
plot_ereefs_slice <- function(slice, var_name='Chl_a_sum', scale_col="spectral", scale_lim=NA, var_units="") {
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
	   d[i] <- earth.dist(slice$cell_intersections[i,'longitude'],slice$cell_intersections[i,'latitude'], slice$cell_intersections[i+1,'longitude'], slice$cell_intersections[i+1,'latitude']) 
	}

  d <- cumsum(d)
	dmin <- c(0, d[1:(numprofiles-1)])
	dmax <- d[1:numprofiles]
	dmin <- t(array(dmin, c(numprofiles, layers)))
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
		ggplot2::ylab('metres above msl') +
		ggplot2::xlab('kilometres from start of transect')
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

#' An internal utility function to find the grid line intersections of a line segment
#
#' importFrom(magrittr,"%>%")
find_intersections <- function(location_latlon, x_grid, y_grid, latitude, longitude, first_point = FALSE) {
	a <- (dim(x_grid) - 1)[1]
	b <- (dim(x_grid) - 1)[2]
	intersected <- array(FALSE, dim=c(a,b))
  # Rearrange to get an array with one row for each cell and one column for each corner of that cell (top left/bottom left/bottom right/top right)
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

	# Define a line through the two points. This probably ought to be corrected for curvature.
	#Y = A*x+b ; A*x +b - Y = 0
  if (lon2!=lon1) {
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
  } else {
    # Special case of a latitudinal line (lon1==lon2)
	  b = lon1
	  # For each edge, the line intersects the edge if the value of Ay+b-x is +ve for one vertex and -ve for the other
	  # Ignore exact hits on vertices.
  
	  # Find grid-cells along this line
	  c1 <- (b - gx[,1]) > 0
	  c2 <- (b - gx[,2]) > 0
	  c3 <- (b - gx[,3]) > 0
	  c4 <- (b - gx[,4]) > 0

	  intersected[c1!=c2] <- TRUE
	  intersected[c2!=c3] <- TRUE
	  intersected[c3!=c4] <- TRUE
	  intersected[c4!=c1] <- TRUE
  }

	# Exclude pesky NA values
  intersected[is.na(c1+c2+c3+c4)] <- FALSE
	intersected[is.na(latitude)] <- FALSE
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
   gx <- gx[intersected, ]
   gy <- gy[intersected, ]
   if (length(llind)==1) {
      gx <- t(gx)
      gy <- t(gy)
   }

   # Find the locations where the edge intersections occur (could have done this all in one go, but this is easier to follow)
   # We do this by running a line (y = cx + k) through each edge of each intersected grid cell and then finding
   # where that line intersects the previously-defined y = Ax+b line between our two transect segment end-points.
   # If lon1==lon2, y is not a function of x (which is always ==lon1)

  if (lon2!=lon1) {
    # Intersects where Ax + b == cx + k
    # y1 = c*x1 + k
    # y2 = c*x2 + k where grid-cell edges run from (x1,y1) to (x2,y2).
    # y1 - y2 = c*(x1 - x2)
    # c = (y1 - y2) / (x1 - x2)
    # We have four edges. For the first, x1 = gx[,2], x2 = gx[,1], y1 = gy[,2], y2 = gy[,1]
    c1 <- (gy[,2] - gy[,1]) / (gx[,2] - gx[,1])
    # k  = y1 - c*x1
    # k <- gy[,2] - c1 * gx[,2]
    # Intersects where Ax + b == cx + k
    # (k - b)  = Ax - cx
    # x = (k - b) / (A - c) 
    #   = (gy[,2] - c1 * gx[,2] - b) / (A - c1)
    xi <- (gy[,2] - b -gx[,2]*c1) / (A -c1)
    # Check that the intersection occurs within the relevant segment (i.e. the cell edge)
    d <- abs((xi<gx[,1]) + (xi<gx[,2]))
    xi[d>1] <- NA                     
    yi <- A * xi + b
    d <- (xi<minlon) + (xi>maxlon) + (yi<minlat) + (yi>maxlat)
    xi[d>0] <- NA
    yi[d>0] <- NA

    c1 <- (gy[,3] - gy[,2]) / (gx[,3] - gx[,2])
    x <- (gy[,3] - b -gx[,3]*c1) / (A -c1)
    d <- abs((x<gx[,3]) + (x<gx[,2]))
    x[d>1] <- NA
    y <- A * x + b
    d <- (x<minlon) + (x>maxlon) + (y<minlat) + (y>maxlat)
    x[d>0] <- NA
    y[d>0] <- NA
    xi <- c(xi, x)
    yi <- c(yi, y)

    c1 <- (gy[,4] - gy[,3]) / (gx[,4] - gx[,3])
    x <- (gy[,4] - b -gx[,4]*c1) / (A -c1)
    d <- abs((x<gx[,4]) + (x<gx[,3]))
    x[d>1] <- NA
    y <- A * x + b
    d <- (x<minlon) + (x>maxlon) + (y<minlat) + (y>maxlat)
    x[d>0] <- NA
    y[d>0] <- NA
    xi <- c(xi, x)
    yi <- c(yi, y)

    c1 <- (gy[,1] - gy[,4]) / (gx[,1] - gx[,4])
    x <- (gy[,1] - b -gx[,1]*c1) / (A -c1)
    d <- abs((x<gx[,1]) + (x<gx[,4]))
    x[d>1] <- NA
    y <- A * x + b
    d <- (x<minlon) + (x>maxlon) + (y<minlat) + (y>maxlat)
    x[d>0] <- NA
    y[d>0] <- NA
    xi <- c(xi, x)
    yi <- c(yi, y)
    yi <- yi[!is.na(xi)]
    xi <- xi[!is.na(xi)]
  } else {
    # Special case of a latitudinal line (lon1==lon2)
    #  x = lon1
    #  On our transect, y is not a function of x but on our grid cell edge, it might be, unless the grid-cell edge is also a latitudinal (North-South) line.

    c1 <- (gy[,2] - gy[,1]) / (gx[,2] - gx[,1])
    # if !is.inf(c1)
    # y1 = c*x1 + k and y2 = c*x2 + k
    # k = y1 - c * x1
    # yi = c * xi + y1 - c * x1
    #    = c * (xi - x1) + y1
    xi <- rep(lon1, length(llind))
    yi <- c1 * (xi - gx[,2]) + gy[,2]
    d <- abs((yi<gy[,2]) + (yi<gx[,1]))
    yi[d>1] <- NA
    d <- (xi<minlon) + (xi>maxlon) + (yi<minlat) + (yi>maxlat)
    xi[d>0] <- NA
    yi[d>0] <- NA

    c1 <- (gy[,3] - gy[,2]) / (gx[,3] - gx[,2])
    x <- rep(lon1, length(llind))
    y <- (gx[,3] - b -gy[,3]*c1) / (A -c1)
    d <- abs((y<gy[,3]) + (y<gx[,2]))
    y[d>1] <- NA
    d <- (x<minlon) + (x>maxlon) + (y<minlat) + (y>maxlat)
    x[d>0] <- NA
    y[d>0] <- NA
    yi <- c(yi, y)
    xi <- c(xi, x)

    c1 <- (gy[,4] - gy[,3]) / (gx[,4] - gx[,3])
    x <- rep(lon1, length(llind))
    y <- (gx[,3] - b -gy[,4]*c1) / (A -c1)
    d <- abs((y<gy[,4]) + (y<gx[,3]))
    y[d>1] <- NA
    d <- (x<minlon) + (x>maxlon) + (y<minlat) + (y>maxlat)
    x[d>0] <- NA
    y[d>0] <- NA
    yi <- c(yi, y)
    xi <- c(xi, x)

    c1 <- (gy[,1] - gy[,4]) / (gx[,1] - gx[,4])
    x <- rep(lon1, length(llind))
    y <- (gx[,4] - b -gy[,1]*c1) / (A -c1)
    d <- abs((y<gy[,1]) + (y<gx[,4]))
    y[d>1] <- NA
    d <- (x<minlon) + (x>maxlon) + (y<minlat) + (y>maxlat)
    x[d>0] <- NA
    y[d>0] <- NA
    yi <- c(yi, y)
    xi <- c(xi, x)

    d <- abs((yi<lat1) + (yi<lat2))
    yi[d>1] <- NA                     
    d <- abs((yi>lat1) + (yi>lat2))
    yi[d>1] <- NA                     
    xi <- xi[!is.na(yi)&!is.infinite(yi)]
    yi <- yi[!is.na(yi)&!is.infinite(yi)]
  }

   if (length(xi)>1) { 
     # determine the order in which the lists are sorted from the direction of the slice
     if (maxlon>minlon) {
        d <- sort(xi, index.return=TRUE)
        yi <- yi[d$ix]
        xi <- xi[d$ix]
        #yi <- c(yi[d$ix][seq(1,length(xi),3)], yi[d$ix[length(xi)]])
        #xi <- c(d$x[seq(1,length(xi),3)], d$x[length(xi)])
        if (longitude[llind[length(llind)]] < longitude[llind[1]] ) {
           yi <- yi[seq(length(xi), 1, -1)]
           xi <- xi[seq(length(xi), 1, -1)]
           #llind <- llind[seq(length(llind), 1, -1)]
        }
     } else {
        d <- sort(yi, index.return=TRUE)
        yi <- yi[d$ix]
        xi <- xi[d$ix]
        #xi <- c(xi[d$ix][seq(1,length(xi),3)], xi[d$ix[length(xi)]])
        #yi <- c(d$x[seq(1,length(xi),3)], d$x[length(xi)])
        if (latitude[llind[length(llind)]] < latitude[llind[1]] ) {
           yi <- yi[seq(length(xi), 1, -1)]
           xi <- xi[seq(length(xi), 1, -1)]
           #llind <- llind[seq(length(llind), 1, -1)]
        }
     }
     # Filter out edges that are beyond the line segment
     location_edges <- data.frame(latitude=yi, longitude=xi) %>% 
       filter(((latitude>=lat1)-(latitude<=lat2))==0, 
              ((longitude>=lon1)-(longitude<=lon2))==0)
     # Each grid-cell edge is included twice (once as the trailing edge of one cell and one as the leading edge of the next cell)
     # So we need to filter out every second edge.
     location_edges <- location_edges[seq(1,dim(location_edges)[1], by=2),]
   } else {
     llind <- NULL
     if (first_point) {
       # no cell edges are intersected, so return the cell-centre of the nearest grid cell
       # I think this is over-complicated. Is there every more than one lat1 and lon1?
       llind <- (data.frame(lat=lat1, lon=lon1) %>%
         mutate(ind = which.min((latitude - lat1)^2 + (longitude - lon1)^2)))$ind
       #llind <- apply(data.frame(latitude=lat1, longitude=lon1) ,1, function(ll) which.min((latitude - ll[1])^2 + (longitude - ll[2])^2))
     }
     location_edges <- NULL
   }
  #browser()
   location_latlon <- data.frame(latitude=latitude[llind], longitude=longitude[llind]) 
   return(list(location_latlon, location_edges, intersected, llind))
}
