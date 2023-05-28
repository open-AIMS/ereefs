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

#' Determines whether the file provided looks like a THREDDS catalog, or an example of daily (e.g. GBR1), 
#' monthly (e.g. GBR4) or other (e.g. RECOM) netcdf file output
#'
#' If the filename contains "catalog.html", it's a THREDDS catalog. 
#' If not, if it ends with ".nc", it's a netcdf file. Any trailing ".html" after the ".nc" part is clipped off.
#' If there are two dashes ('-') in the last 13 characters of the netcdf filename, it is assumed to be a daily output file. 
#' If there is one, it is assumed to be monthly. If there are no dashes, it's something else, such as a RECOM file.
#' If the output is not consistently in monthly or daily netcdf files but switches between the two, it is best to provide 
#' a THREDDS catalog URI.
#'
#' @param filename Name of the file to examine
#' @return a vector of strings. 1st value is either "nc" (for an individual netcdf file), "mnc" for a meta-netcdf file (list of .nc files), "xlm" or "ncml" (for catalogs). 
#' If the first value is "nc", the second value is "1km" for GBR1, "4km" for GBR4 or "recom" (anything else). Otherwise, second value is "unknown".
#' @export
get_ereefs_case <- function(filename) {
  # If the filename ends in ".nc.html", truncate to remove the ".html" part. These two lines should have already been taken care of
  # by substitute_filename():
  if (stringr::str_ends(filename, stringr::fixed(".nc.html"))) filename <- stringr::str_sub(filename, 1, -6)
  if (stringr::str_ends(filename, stringr::fixed("catalog.html"))) filename <- paste0(stringr::str_sub(filename, 1, -5), 'xml')

  if (stringr::str_ends(filename, stringr::fixed(".nc"))) {
    ereefs_case <- c("nc")
  } else if (stringr::str_ends(filename, stringr::fixed(".mnc"))) {
    ereefs_case <- c("mnc", "unknown")
    stop(".mnc files not yet implemented.")
  } else if (stringr::str_ends(filename, ".xml")) {
    ereefs_case <- c("xml", "unknown")
    stop("Should not get to ereefs_case xml. Should have already been converted by substitute_filename() to .ncml")
  } else if (stringr::str_ends(filename, ".ncml")) {
    ereefs_case <- c("ncml", "unknown")
  } else {
    stop(paste("Filename format (", filename, ") is not recognised as a netcdf or THREDDS catalog file"))
  }

  if (ereefs_case[1] == "nc") {
    lastfew <- stringr::str_sub(filename, start=-13, end=-3)
    dashcount <- stringi::stri_count_fixed(lastfew, '-')
    if (dashcount==2) {
	    ereefs_case <- c(ereefs_case, '1km')
    } else if (dashcount==1) {
	    ereefs_case <- c(ereefs_case, '4km')
    } else {
	    ereefs_case <- c(ereefs_case, "recom")
    }
  }
  return(ereefs_case)
}

#' Utility function copied directly from the package 'tis' by Jeff Hallman (https://github.com/cran/tis/blob/master/R/assignList.R)
#' Assigns the values in a list to variables in an environment. The variable names are taken from the names of the list, so all of the elements of the list must have non-blank names.

assignList <- function(aList, pos = -1, envir = as.environment(pos), inherits = FALSE){
  if(is.null(nms <- names(aList))) stop("names(aList) is NULL")
  if(any(nms == "")) stop("blank name")
  if(missing(envir) && pos < 0)
    envir <- parent.frame(-pos)
  for(nm in nms)
    assign(x = nm, value = aList[[nm]], envir = envir, inherits = inherits)
}

#' Set up various parameters needed by many of the data extraction and plotting functions
#'
#' @param start_date The date from which to start data extraction, in any format accepted by get_date_time(). Can be any of:
#'              c(yr, month, day)
#'              c(yr, month, day, hour) (can also add minutes and seconds to the vector)
#'              POSIX date (e.g., as.Date('1970-01-01', origin='1970-01-01'))
#'              YYYYMMDD, YYYYMMDDhh, YYYYMMDDhhmm or YYYYMMDDhhmmss
#'              character format, e.g. '1970-01-01'
#' @param end_date The end date of the data extraction request, as above.
#' @param input_file The URI or filename from which to extract data, in any format accepted by substitute_filename(), or "menu" or "catalog".
#' @param var_names A list of eReefs variable names for which to extract data.
#' @return list of parameter values that includes:
#'              input_file, the input filename, after transformation by substitute_filename()
#'              ereefs_case, information about the format of the input file, as output by get)ereefs_case()
#'              input_stem, the 'stem' of the filename in case the parent function needs to dynamically calculate filenames from an example netcdf filename.
#'              start_date, start_date transformed to a chron date
#'              start_day, the day of the month in start_date
#'              start_tod, the time of day in start_date (0.5 if not given)
#'              start_month, the integer month of the yr in start_date
#'              start_yr, the integer yr in start_date
#'              end_date, end_date transformed to a chron date
#'              end_day, the day of the month in end_date
#'              end_tod, the time of day in end_date (0.5 if not given)
#'              end_month, the integer month of the yr in end_date
#'              end_yr, the integer yr in end_date
#'              mths, a vector of months for which to extract data
#'              yrs, a vector of yrs for which to extract data
#'              var_list, a list of variable names in var_names
#'              ereefs_origin, the origin date specified for Julian dates in the netcdf file (assumed to be 1990-01-01 if not specified)
#'              ds, a date-time vector for extracted data
#'              julian_date, the time vector in the raw format provided in the netcdf file
#'              spatial_grid, a tibble containing the grid latitudes and longitudes
get_params <- function(start_date, end_date, input_file, var_names) {

  input_file <- substitute_filename(input_file)
  # Check whether netcdf output files are daily (case 1), monthly (case 4) or something else (case 0)
  ereefs_case <- get_ereefs_case(input_file)
  input_stem <- get_file_stem(input_file)

  # Dates to plot
  start_date <- get_date_time(start_date)
  start_day <- lubridate::day(start_date)
  start_tod <- as.numeric(start_date - floor_date(start_date, "day"))
  start_month <- lubridate::month(start_date)
  start_yr <- lubridate::year(start_date)

  end_date <- get_date_time(end_date)
  end_day <- lubridate::day(end_date)
  end_month <- lubridate::month(end_date)
  end_yr <- lubridate::year(end_date)

  
  if (start_date > end_date) {
    stop('start_date must preceed end_date')
  }

  if (start_yr==end_yr) {
      mths <- start_month:end_month
      yrs <- rep(start_yr, length(mths))
  } else if ((start_yr + 1) == end_yr) {
      mths <- c(start_month:12, 1:end_month)
      yrs <- c(rep(start_yr, 12 - start_month + 1), rep(end_yr, end_month))
  } else {
      mths <- c(start_month:12, rep(1:12, end_yr - start_yr - 1), 1:end_month)
      yrs <- c(rep(start_yr, 12 - start_month + 1), 
                 rep((start_yr + 1) : (end_yr - 1), each=12),
                 rep(end_yr, end_month))
  }

  var_list <- paste(var_names, collapse=",")

  if (ereefs_case[2] == '4km') {
      input_file <- paste0(input_stem, format(as.Date(paste(start_yr, start_month, 1, sep='-')), '%Y-%m'), 
			  '.nc')
    dum1 <- get_origin_and_times(input_file)
    ereefs_origin <- dum1[[1]]
    ds <- dum1[[2]]
    julian_date <- dum1[[3]]
  } else if (ereefs_case[2] == '1km') {
    input_file <- paste0(input_stem, format(as.Date(paste(start_yr, start_month, start_day, sep='-')), '%Y-%m-%d'), 
		  '.nc')
     dum1 <- get_origin_and_times(input_file)
     ereefs_origin <- dum1[[1]]
     ds <- NULL
     julian_date <- dum1[[3]]
     # We don't need ds in this case because we are taking just one output per day
  } else {
     # We are looking at a single netcdf file such as a RECOM output file or a ncml file
     input_file <- input_file
     dum1 <- get_origin_and_times(input_file)
     ereefs_origin <- dum1[[1]]
     ds <- dum1[[2]]
     julian_date <- dum1[[3]]
  }

  # Activate the grid that has two dimensions and extract the latitude and longitude variables:
  nc <- tidync::tidync(input_file) 
  nc <- nc %>% tidync::activate(as.character(dplyr::inner_join(tidync::hyper_grids(nc), data.frame('ndims' = 2), by = 'ndims') %>% dplyr::select('grid')))
  if ('latitude' %in% (tidync::hyper_vars(nc)$name)) {
    spatial_grid <- nc %>% tidync::hyper_tibble(c('latitude', 'longitude'))
  } else {
    spatial_grid <- nc %>% tidync::hyper_tibble(c('x_centre', 'y_centre')) %>% 
                           dplyr::rename(latitude = x_centre, longitude = y_centre)
  }

  params <- list(input_file = input_file,
                 ereefs_case = ereefs_case,
                 input_stem = input_stem,
                 start_date = start_date,
                 end_date = end_date,
                 start_day = start_day,
                 start_tod = start_tod,
                 start_month = start_month,
                 start_yr = start_yr,
                 end_date = end_date,
                 end_day = end_day,
                 end_month = end_month,
                 end_yr = end_yr,
                 mths = mths,
                 yrs = yrs,
                 var_list = var_list,
                 ereefs_origin = ereefs_origin,
                 ds = ds, 
                 julian_date = julian_date,
                 spatial_grid = spatial_grid)
}

#' Opens the netcdf file, input_file, and extracts the origin (reference date/time) and time dimension.
#'
#' @param input_file The name of the netcdf file (in standard or simple EMS netcdf format) from which
#' to extract the data
#' @return A list containing ereefs_origin (the reference data/time used in the netcdf file) and the time series, ds
get_origin_and_times <- function(input_file) {
    nc <- tidync::tidync(input_file) %>% tidync::activate('time')
    if ('t' %in% (tidync::hyper_vars(nc)$name)) {
      posix_origin <- stringi::stri_datetime_parse(ncmeta::nc_att(input_file ,'t', 'units')$value, "'days since 'yyyy-MM-dd HH:mm:ss", tz = "Etc/GMT-10")[1]
      ds <- (nc %>% tidync::hyper_array('t'))[[1]]
    } else { 
      posix_origin <- stringi::stri_datetime_parse(ncmeta::nc_att(input_file ,'time', 'units')$value, "'days since 'yyyy-MM-dd HH:mm:ss", tz = "Etc/GMT-10")[1]
      ds <- (nc %>% tidync::hyper_array('time'))[[1]]
    }
    julian_date <- ds
    ereefs_origin <- posix_origin

    ds <- ereefs_origin + seconds(ds*86400)
    return(list(ereefs_origin, ds, julian_date))
}

#' Returns an appropriate x_grid, y_grid and z_grid
#'
#' @param filename Name of the file to examine
#' @param input_grid Optional alternative file from which to load x_grid, y_grid and z_grid
#' @return a list containing x_grid, y_grid and z_grid
#' @export
get_ereefs_grids <- function(filename, input_grid=NA) {
	if (!is.na(input_grid)) {
    # The user has provided the name of a netcdf file from which to read the grids
    nc <- tidync::tidync(input_grid) 
    nc <- nc %>% tidync::activate(as.character(dplyr::inner_join(tidync::hyper_grids(nc), data.frame('ndims' = 2), by = 'ndims') %>% dplyr::select('grid')))
    x_grid <- (nc %>% tidync::hyper_array('x_grid'))[[1]]
    y_grid <- (nc %>% tidync::hyper_array('y_grid'))[[1]]
    z_grid <- (nc %>% tidync::hyper_array('z_grid'))[[1]]
	} else {
    nc <- tidync::tidync(filename) 
    nc <- nc %>% tidync::activate(as.character(dplyr::inner_join(tidync::hyper_grids(nc), data.frame('ndims' = 2), by = 'ndims') %>% dplyr::select('grid')))
    if ('x_grid' %in% (tidync::hyper_vars(nc)$name)) {
      # filename is in standard EMS netcdf output format. x_grid, y_grid and z_grid are provided.
      x_grid <- (nc %>% tidync::hyper_array('x_grid'))[[1]]
      y_grid <- (nc %>% tidync::hyper_array('y_grid'))[[1]]
      z_grid <- (nc %>% tidync::hyper_array('z_grid'))[[1]]
    } else {
      # Simple format netcdf file. grids are not provided. Check the dimensions to work out 
      # whether it looks like GBR4 or GBR1 and if so, load x-grid and y_grid from saved data. 
      # Otherwise, the user must provide another filename (input_grid) that contains x_grid, 
      # y_grid and z_grid. 
      dims <- nc %>% tidync::hyper_dims()
      ilen <- dims %>% dplyr::inner_join(data.frame(name=c('i', 'i_grid')), by = 'name') %>% dplyr::select('length') # Find either 'i' or 'i_grid' in the dimension names and return the length of that dimension
      jlen <- dims %>% dplyr::inner_join(data.frame(name=c('j', 'j_grid')), by = 'name') %>% dplyr::select('length')
		  if ((ilen==510)&(jlen==2389)) {
		    x_grid <- gbr1_x_grid
		    y_grid <- gbr1_y_grid
		    z_grid <- gbr1_z_grid
		  } else if ((ilen==600)&(jlen==180)) {
		    x_grid <- gbr4_x_grid
		    y_grid <- gbr4_y_grid
		    z_grid <- gbr4_z_grid
		  } else {
		    stop(paste('x_grid, y_grid and z_grid not found and grid size not recognised as GBR1 or GBR4.',
                   'Plotting using cell centres (latitude and longitude) has not yet been implemented.',
                   'Please specify input_grid (a standard format EMS file or other netcdf file that ',
                   'contains x_grid, y_grid and z_grid).'))

		  }
                }
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
  if (case[1]!='nc') {
    return(NA)
  }
  if (is.na(case[2]) | case[2] == "recom" | case[2] == "unknown") {
     # A netcdf filename, but not obviously 1km or 4km eReefs. Might be RECOM or another application
	  file_stem <- substr(filename, start=1, stop=stringi::stri_locate_last(filename, regex='.nc')[1]-1)
  } else if (case[2]=='1km') {
	  file_stem <- substr(filename, start=1, stop=stringi::stri_locate_last(filename, regex='.nc')[1]-11)
  } else {
     # must be 4km
	  file_stem <- substr(filename, start=1, stop=stringi::stri_locate_last(filename, regex='.nc')[1]-8)
  }
 return(file_stem)
}

#' Check whether the given filename is a shortcut and if so, set the full filename
#'
#' Utility function for the eReefs package
#' @param input_file is the URL or file location of any of the EMS output files or a THREDDS catalog URI. 
#'        Defaults to a menu selection based on current NCI catalogs. Can also be set to "nci", "menu" or "catalog" for the same behaviour.
#'        Set to "old_menu" to provide old menu options instead of menu options from the NCI catalog.
#'        Numeric values are interpreted as references to selections available from the old menu.
#'        Short codes can be used for some options (codenames as used in https://research.csiro.au/ereefs/models/model-outputs/access-to-raw-model-output/ )
#' @return input_file
#' @export
substitute_filename <- function(input_file = "catalog") {
  # Try to work out what the user wants if it looks a bit like a THREDDS catalog but isn't clear
  if (stringr::str_ends(input_file, stringr::fixed("catalog.html"))) input_file <- paste0(stringr::str_sub(input_file, 1, -5), 'xml')
  if (((stringr::str_starts(input_file, "http")&(stringr::str_ends(input_file, "/"))))) input_file <- paste0(input_file, "catalog.xml")

  choices <- NULL
  if ((input_file == "old_menu")|(input_file == "choices") | is.numeric(input_file)) {
    # Set up the menu of choices to display for interactive input, or to select from if the user has specified an option number
    short_names  <- c("GBR4-v2.0",
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
                      "GBR4_BGC-v3.0 Dcrt",
                      #"GBR4_BGC-v3.0 Dnrt",
                      "GBR4_BGC-v3.1",
                      "GBR4_BGC-v3.0",
                      "GBR1_BGC-v3p2surf",
                      "GBR1_BGC-v3p2",
                      "GBR4_BGC-v3p2nrtsurf",
                      "GBR4_BGC-v3p2nrt",
                      "menu") 
    if (input_file == "choices") { 
      # The user just wants a list of options
      print(short_names)
      stop()
    } 
  } else if ((input_file == "catalog") | (input_file == "nci") | (input_file == "menu")) { 
    # Get the list of eReefs data sevices from the NCI server:
    catalog <- thredds::get_catalog("https://dapds00.nci.org.au/thredds/catalogs/fx3/catalog.xml")
    choices <- catalog$get_dataset_names()
    print("Refer to https://research.csiro.au/ereefs/models/models-about/biogeochemical-simulation-naming-protocol/ for naming conventions.")
  } else if (stringr::str_ends(input_file, "xml")) {
    catalog <- thredds::get_catalog(input_file)
    choices <- catalog$get_dataset_names()
  } else if (input_file=="old_menu") {
    choices <- c("Latest release 4km grid hydrodynamic model (Sept 2010-pres.)", 
                 "Archived 4km biogeochemical model hindcast v2.0 (Sept 2010 - Oct 2016)",
                 "Archived 4km biogeochemical model near real time v2.0 (Oct 2016 - Nov 2019)",
                 "Archived Pre-industrial catchment scenario 4km BGC (Sept 2010 - Oct 2016)",
                 "Latest release passive river tracers (Sept 2010 - pres.)",
                 "Latest release 1km grid hydrodynamic model (Dec 2014 - pres.)",
                 "Latest release 1km grid passive river tracers (Dec 2014 - pres.)",
                 "Archived 4km hydro (v 1.85, Sept 2010-pres.)",
                 "Archived 4km bgc (v926, Sept 2010 - Dec 2014)",
                 "Archived 4km bgc (v924, Sept 2010 - Sept 2017)",
                 "Archived 1km hydro (v 1.71, Dec 2014 – Apr 2016)",
                 "Archived 1km bgc (v924, Dec 2014 – Nov 2019)",
                 "Archived 4km biogeochemical model hindcast v3.0 (Dec 2010 - Oct 2018)",
                 "Latest release 4km biogeochemical model v3.1 (Dec 2010 - Apr 2019)",
                 "CSIRO login required: GBR1 NRT BGC 3p0 3D (2018-09-02 to 2019-01-30)",
                 "CSIRO login required: Latest release GBR1 NRT BGC 3p2 surface (16 Oct 2019 - pres., 3x/day)",
                 "CSIRO login required: Latest release GBR1 NRT BGC 3p2 3D (16 Oct 2019 - pres., daily)",
                 "CSIRO login required: Latest release GBR4 NRT BGC surface (Oct 2019 - May 2020, 4x/day)",
                 "CSIRO login required: Latest release GBR4 NRT BGC 3D (Oct 2019 - May 2020, daily)")
  }
  if (is.numeric(input_file)) {
     input_file <- short_names[input_file]
  } else if (!is.null(choices)) {
    selection <- utils::menu(choices)
    if (input_file == "old_menu") {
      input_file <- short_names[selection]
    } else {
      input_file <- paste0(stringr::str_extract(catalog$url, "https://.[^/]*"), catalog$list_services()$ncdods[3], catalog$get_datasets()[[selection]]$get_url())
    }
  }
  # Perhaps the user has manually specified (or selected from old_menu) one of the (obsolescent) official run labels from ereefs.info
  input_file <- dplyr::case_when(
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
      input_file == "GBR4_BGC-v3.0 Dcrt" ~ "http://dapds00.nci.org.au/thredds/dodsC/fx3/gbr4_bgc_GBR4_H2p0_B3p0_Chyd_Dcrt/gbr4_bgc_all_simple_2018-10.nc",
      input_file == "GBR4_BGC-v3.1" ~ "https://regional-models.ereefs.info/thredds/dodsC/GBR4_H2p0_B3p1_Cq3b_Dhnd/all/gbr4_bgc_all_simple_2012-10.nc",
      input_file == "GBR4_BGC-v3.0" ~ "http://oa-62-cdc.it.csiro.au:8087/opendap/cache/gbr1/bgc/3.0/all/gbr1_bgc_all_2018-09-02.nc",
      input_file == "GBR1_BGC-v3p2surf" ~ "http://oa-62-cdc.it.csiro.au:8087/opendap/cache/gbr1/bgc/nrt/surf/gbr1_bgc_surf_2019-10-16.nc",
      input_file == "GBR1_BGC-v3p2" ~ "http://oa-62-cdc.it.csiro.au:8087/opendap/cache/gbr1/bgc/nrt/all/gbr1_bgc_all_2019-10-16.nc",
      input_file == "GBR4_BGC-v3p2nrtsurf" ~ "http://oa-62-cdc.it.csiro.au:8087/opendap/cache/gbr4/bgc/nrt/gbr4_bgc_surf_2019-10.nc",
      input_file == "GBR4_BGC-v3p2nrt" ~ "http://oa-62-cdc.it.csiro.au:8087/opendap/cache/gbr4/bgc/nrt/gbr4_bgc_all_2019-10.nc",
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
      TRUE ~ input_file)

  # If necessary, remove trailing ".html" from an input file name that is not a catalog.
  if (stringr::str_ends(input_file, ".nc.html")) input_file <- stringr::str_sub(input_file, 1, -6)
    
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
#' All variables to be extracted must have the same number of dimensions (e.g., either all 2D or all 3D variables).
#'
#' If you run into memory constraints, consider grouping points to be extracted within regions, and calling this once
#' for each region.
#'
#' @return a data tibble containing extracted variables, including dates and locations
#' @param var_names either a single character value or a vector specifying the short names for variables that you 
#'        want from the netcdf file. Defaults to c('Chl_a_sum', 'TN').
#' @param geocoordinates is a data frame containing the decimal latitude and longitude of a single desired location, or a vector containing
#'        a single latitude and longitude location. geocoordinates can also be set to "mmp" to extract time-series at all 
#'        Marine Monitoring Program sites. Defaults to c(-23.39189, 150.88852). Can also include other variables, e.g. for the name of each site.
#' @param layer is the vertical grid layer to extract, or 'surface' to get the surface value, 'bottom' to get the
#'        value in the cell at the bottom of the water column, or 'integrated' to get a depth-integrated (mean) value.
#'        Specify a negative value to indicate a specified depth (in metres) below MSL.
#'        Use get_ereefs_depth_specified_ts() instead if you want to specify a depth 
#'        below the free surface instead of a layer number. Defaults to 'surface'.
#' @param start_date date  for animation. Can either be a) a vector of integers in the format 
#'	       c(yr, month, day); b) a date obtained e.g. from as.Date(); or c) a character string 
#'        formatted for input to as.Date(). Defaults to c(2010,12,31).
#' @param end date for animation. Can either be a) a vector of integers in the format 
#'	       c(yr, month, day); b) a date obtained e.g. from as.Date(); or c) a character string 
#'        formatted for input to as.Date(). Defaults to c(2016,10,31).
#' @param input_file is the URL or file location of any of the EMS output files or a THREDDS catalog URI. 
#'        Defaults to a menu selection based on current NCI catalogs. Can also be set to "nci", "menu" or "catalog" for the same behaviour.
#'        Set to "old_menu" to provide old menu options instead of menu options from the NCI catalog.
#'        Numeric values are interpreted as references to selections available from the old menu.
#'        Short codes can be used for some options (codenames as used in https://research.csiro.au/ereefs/models/model-outputs/access-to-raw-model-output/ )
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
#' @param verbosity How much information to display along the way (0 to 2. Default is 1).
#' @export
#' @examples
#' \dontrun{
#' get_ereefs_ts('Chl_a_sum', geocoordinates=data.frame(latitude=c(-23.39189,-18), longitude=c(150.88852, 147.5)), layer='surface', start_date=c(2010,12,31),end_date=c(2011,1,5), input_file=2)
#' get_ereefs_ts(var_names=c('Tricho_N', 'DIP', 'TP'), geocoordinates=data.frame(latitude=c(-23.39189,-18), longitude=c(150.88852, 147.5)), layer='bottom', start_date="2012-07-03",end_date="2012-07-05", input_file='GBR4_BGC-v2.0 Chyb Dcrt')
#' get_ereefs_ts(var_names=c('ZooL_N', 'ZooS_N'), geocoordinates=data.frame(latitude=c(-23.39189,-18), longitude=c(150.88852, 147.5)), layer=45, start_date=c(2010,12,31),end_date=c(2011,1,5), input_file="http://dapds00.nci.org.au/thredds/dodsC/fx3/gbr4_bgc_GBR4_H2p0_B2p0_Cpre_Dcrt/gbr4_bgc_simple_2016-06.nc")
#' }

get_ereefs_ts <- function(var_names=c('Chl_a_sum', 'TN'), 
                          geocoordinates=c(-23.39189, 150.88852), 
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
  library(dplyr)
  if (is.null(dim(geocoordinates))) {
     geocoordinates <- array(geocoordinates, c(1,2))
  }

  # Get parameter values and assign results from returned list to relevant variable names
  # This assigns input_file, ereefs_case, input_stem, start_date, end_date, start_tod, start_month, start_yr,
  # end_date, end_day, end_month, end_yr, mths, yrs, var_list, ereefs_origin, ds, and spatial_grid
  assignList(get_params(start_date, end_date, input_file, var_names))
  if (end_date==start_date) end_date <- end_date + 1

  if (layer=='integrated') return(get_ereefs_depth_integrated_ts(var_names, geocoordinates, start_date, end_date, input_file, input_grid, eta_stem, override_positive))
  if ((layer=='bottom')&&(length(geocoordinates[,1])>1)) stop('Only one location can be given if layer==bottom')
#  if (layer=='bottom') return(get_ereefs_bottom_ts(var_names, geocoordinates, start_date, end_date, input_file, input_grid, eta_stem, override_positive))
  if (geocoordinates[[1]][1]=="mmp") {
     geocoordinates <- data.frame(latitude=mmp_sites$latitude, longitude=mmp_sites$longitude)
     mmp <- TRUE
  } else {
     mmp <- FALSE
  }

  # Make sure the geocoordinates are in the desired data.frame format. Warn if column names are absent or unclear.
  if (!("data.frame" %in% class(geocoordinates))) {
    # geocoordinates has been provided as an array/matrix. Coerce it into a data frame for consistency.
    geocoordinates <- data.frame(latitude = geocoordinates[,1], longitude = geocoordinates[,2])
  }
  ilat <- which((names(geocoordinates) == "latitude")|(names(geocoordinates) == "lat"))
  ilon <- which((names(geocoordinates) == "longitude")|(names(geocoordinates) == "lon"))
  if (length(ilat)) {
    geocoordinates <- data.frame(latitude = geocoordinates[, ilat], longitude = geocoordinates[, ilon])
  } else {
    warning("Assuming that the first column of geocoordinates is latitude and the second column is longitude")
  }
  names(geocoordinates)[1:2] <- c("latitude", "longitude")

  geocoordinates <- geocoordinates %>% dplyr::rowwise() %>% 
                                       dplyr::mutate(grid_index = which.min(earth.dist(latitude, longitude, spatial_grid$latitude, spatial_grid$longitude)))
  geocoordinates[, c("i", "j", "cell_latitude", "cell_longitude")] = spatial_grid[geocoordinates$grid_index, c("i", "j", "latitude", "longitude")]

  # Initialise
  results <- NULL

  # Loop through relevant daily or monthly eReefs files to extract the data
  ndims <- NA
  layer_actual <- rep(NA, length(var_names))
  if ((ereefs_case[2] == "1km")|(ereefs_case[2] == "4km")) {
    i <- 0
    mcount <- 0
    if (verbosity>0) pb <- txtProgressBar(min = 0, max = length(mths), style = 3)
    for (month in mths) {
      mcount <- mcount + 1
      yr <- yrs[mcount]
      if (mcount == 1) {
         from_day <- start_day
      } else {
         from_day <- 1
         start_tod <- 0
      }
      if ((start_yr==end_yr)&&(start_month==end_month)) {
         day_count <- end_day - start_day + 1
      } else if (mcount == 1) {
         day_count <- daysIn(as.Date(paste(yr, month, 1, sep='-'))) - start_day + 1
      } else if (mcount == (length(mths))) {
         day_count <- end_day
      } else {
         day_count <- daysIn(as.Date(paste(yr, month, 1, sep='-')))
      }
      if (ereefs_case[2] == '4km') { 
        fileslist <- 1
        input_file <- paste0(input_stem, format(as.Date(paste(yr, month, 1, sep="-")), '%Y-%m'), '.nc')
	      day_count <- day_count / as.numeric(ds[2]-ds[1])
        if ((ds[length(ds)] - ds[1]) > 31.5) warning('Filename looks like a monthly output file (i.e. contains two dashes) but file contains more than a month of data.')
        if ((ds[length(ds)] - ds[1]) < 27) warning('Filename looks like a monthly output file (i.e. contains two dashes) but file contains less than a month of data.')
        if(ds[2]==ds[1]) stop(paste('Error reading time from', input_file, '(t[2]==t[1])'))
        if (day_count > length(ds)) {
          warning(paste('end_date', end_date, 'is beyond available data. Ending at', ds[length(ds)]))
          day_count <- length(ds)
        }
      } else if (ereefs_case[2] == '1km') {
	      fileslist <- from_day:(from_day+day_count-1)
	      from_day <- 1
	      day_count <- 1
      } else stop("Shouldn't happen: ereefs_case not recognised")

      for (dcount in 1:length(fileslist)) {
        if (ereefs_case[2] == '1km') {
	        input_file <- paste0(input_stem, format(as.Date(paste(yr, month, dcount, sep="-")), '%Y-%m-%d'), '.nc')
        }
        # This is a bit roundabout but we first need to activate by variable name, then check which grid is active and activate the whole grid
        nc <- tidync::tidync(input_file) %>% tidync::activate(var_names[1]) %>% tidync::activate(tidync::active(nc), select_var = var_names)
        d <- get_origin_and_times(input_file)[[2]][from_day:(from_day+day_count-1)]
        im1 = i+1
        i <- i + length(d)

        if (is.na(ndims)) {
           # We don't yet know the dimensions of the variable, so let's get them
           ndims <- dim(nc %>% hyper_dims())[1]
         }
         if (layer == 'surface') { 
            layer_actual <- as.integer(nc %>% tidync::hyper_dims() %>% dplyr::filter(name=="k") %>% dplyr::select("length"))
         } else if (layer == 'bottom') { 
            # find the bottom layer -- assume it is the first layer that is not NA
            if (ndims == 4) { 
               dum1 <- nc %>% tidync::hyper_filter(i = dplyr::between(i, min(geocoordinates$i), max(geocoordinates$i)),
                                                   j = dplyr::between(j, min(geocoordinates$j), max(geocoordinates$j)),
                                                   time = index == from_day) %>% 
                              tidync::hyper_array(select_var = var_names)
               if (length(which(!is.na(dum1)))==0) stop('Location given is a dry cell at start time. It is probably on land.')
               layer <- min(which(!is.na(dum1)))
               layer_actual <- layer
               if (verbosity>0) print(paste('bottom layer = ', layer))
            } 
         } else if (layer < 0) { 
           z_grid <- get_ereefs_grids(input_file, input_grid)[['z_grid']] 
           layer <- max(which(z_grid<layer)) 
           layer_actual <- layer 
         } else {
           layer_actual <- layer 
         }

        if (ndims == 4) {
           wc <- nc %>% tidync::hyper_filter(i = dplyr::between(i, min(geocoordinates$i), max(geocoordinates$i)),
                                             j = dplyr::between(j, min(geocoordinates$j), max(geocoordinates$j)),
                                             k = k==layer_actual, 
                                             time = dplyr::between(index, from_day, from_day+day_count-1)) %>%
                        tidync::hyper_tibble(select_var = var_names)

        } else {
           wc <- nc %>% tidync::hyper_filter(i = dplyr::between(i, min(geocoordinates$i), max(geocoordinates$i)),
                                             j = dplyr::between(j, min(geocoordinates$j), max(geocoordinates$j)),
                                             time = dplyr::between(index, from_day, from_day+day_count-1)) %>%
                        tidync::hyper_tibble(select_var = var_names)
        }
        results <- rbind(results, dplyr::left_join(geocoordinates, wc, by = c("i", "j")))
        if (verbosity>0) setTxtProgressBar(pb,mcount)
      }
    }
    if (verbosity>0) close(pb)
  } else { 
    # Single nc or ncml file. No need to loop through months, but might need to check for monotonicity of timeseries 
    print(paste('No waitbar shown in this case. If you have asked for a large dataset, please be patient.',
                  'If you have asked for a long time-series from a THREDDS server,',
                  'you might encounter transfer errors, but it should self-correct.',
                  'If not, consider breaking it up into shorter sections or pointing to an individual netcdf file.'))
    if (end_date > max(ds, na.rm=TRUE)) {
      warning(paste('end_date', end_date, 'is beyond available data. Ending at', max(ds, na.rm=TRUE)))
      end_date <- max(ds, na.rm=TRUE)
    }
    if (start_date < min(ds, na.rm=TRUE)) {
      warning(paste('start_date', start_date, 'is beyond available data. Starting at', min(ds, na.rm=TRUE)))
      start_date <- min(ds, na.rm=TRUE)
    }
    nc <- tidync::tidync(input_file) %>% tidync::activate(var_names[1]) %>% tidync::activate(tidync::active(nc), select_var = var_names)
    if (is.na(ndims)) {
      # We don't yet know the dimensions of the variables, so let's get them
      dims <- nc %>% hyper_dims()
      ndims <- dim(dims)[1]
    }
    if (layer == 'surface') { 
        layer_actual <- as.integer(dims %>% inner_join(data.frame(name = c('k', 'k_grid')), by="name") %>% select('length'))
    } else if (layer == 'bottom') { 
        # find the bottom layer -- assume it is the first layer that is not NA
        if (ndims == 4) { 
          dum1 <- nc %>% tidync::hyper_filter(i = dplyr::between(i, min(geocoordinates$i), max(geocoordinates$i)),
                                              j = dplyr::between(j, min(geocoordinates$j), max(geocoordinates$j)),
                                              time = index == 1) %>%
                         tidync::hyper_array(select_var = var_names[1])
          if (length(which(!is.na(dum1[[1]])))==0) stop('Location given is a dry cell at first time step. It is probably on land.')
          layer <- min(which(!is.na(dum1[[1]])))
          layer_actual <- layer
          if (verbosity>0) print(paste('bottom layer = ', layer))
        } 
      } else if (layer < 0) { 
        z_grid <- get_ereefs_grids(input_file, input_grid)[['z_grid']] 
        layer <- max(which(z_grid<layer)) 
        layer_actual <- layer 
      } else {
        layer_actual <- layer 
      }

      if (ndims == 4) {
         wc <- nc %>% tidync::hyper_filter(i = dplyr::between(i, min(geocoordinates$i), max(geocoordinates$i)),
                                           j = dplyr::between(j, min(geocoordinates$j), max(geocoordinates$j)),
                                           k = k==layer_actual, 
                                           time = dplyr::between(time, as.numeric(start_date - ereefs_origin), as.numeric(end_date - ereefs_origin))) %>%
                      tidync::hyper_tibble(select_var = var_names)
      } else {
         wc <- nc %>% tidync::hyper_filter(i = dplyr::between(i, min(geocoordinates$i), max(geocoordinates$i)),
                                           j = dplyr::between(j, min(geocoordinates$j), max(geocoordinates$j)),
                                           time = dplyr::between(time, as.numeric(start_date - ereefs_origin), as.numeric(end_date - ereefs_origin))) %>%
                      tidync::hyper_tibble(select_var = var_names)
      }
      results <- dplyr::left_join(geocoordinates, wc, by = c("i", "j"))
    }
    results <- results %>% mutate(time = ereefs_origin + seconds(time*86400))
    results <- results %>% select(-grid_index)

  if (mmp) results <- left_join(results, mmp_sites, by = c('latitude', 'longitude'))
  return(results)
}

#' Extracts time-series of selected variables from the bottom water-column cell at specified locations from eReefs output files
#'
#' Now just a wrapper to get_ereefs_ts().
#' By Barbara Robson (AIMS).
#'
#' @return a data frame containing the dates and values of extracted variables
#' @param var_names either a single character value or a vector specifying the short names for variables that you 
#'        want from the netcdf file. Defaults to c('Chl_a_sum', 'TN').
#' @param geocoordinates is a vector containing the decimal latitude and longitude of the desired location. If 
#'        you want to specify an x-y grid coordinate instead of a latitude and longitude, you can: to do this, 
#'        is.integer(geocoordinates) must be TRUE. Defaults to c(-23.39189, 150.88852).
#' @param start_date date  for animation. Can either be a) a vector of integers in the format 
#'	c(yr, month, day); b) a date obtained e.g. from as.Date(); or c) a character string 
#'      formatted for input to as.Date(). Defaults to c(2010,12,31).
#' @param end date for animation. Can either be a) a vector of integers in the format 
#'	c(yr, month, day); b) a date obtained e.g. from as.Date(); or c) a character string 
#'      formatted for input to as.Date(). Defaults to c(2016,12,31).
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
#'      "http://dapds00.nci.org.au/thredds/dodsC/fx3/gbr4_hydro_all/gbr4_all_2016-09.nc".
#' @param eta_stem The URI or file location of the model output files that contains the surface elevation (eta), minus the
#'       date components of the filename in the case of GBR1 or GBR4 files, and ommitting the file extension, ".nc". Needed
#'       only if eta is not in the files indicated by input_stem (e.g. some GBR1 bgc files).
#' @param override_positive Reverse the value of the "positive" attribute of botz for BGC files, assuming that it is
#'       incorrect. Default FALSE
#' @export
#' @examples
#' \dontrun{
#' get_ereefs_bottom_ts(c('Chl_a_sum', 'NH4'), geocoordinates=data.frame(latitude=-23.39189, longitude=150.88852), layer='surface', start_date=c(2010,12,31),end_date=c(2011,1,5), input_file=2)
#'}
get_ereefs_bottom_ts <- function(var_names=c('Chl_a_sum', 'TN'), 
                          geocoordinates=c(-23.39189, 150.88852), 
                          layer='bottom', 
		                      start_date = c(2010,12,31), 
		                      end_date = c(2016,10,31), 
                          input_file = "menu",
                          #input_file = "http://dapds00.nci.org.au/thredds/dodsC/fx3/gbr4_bgc_GBR4_H2p0_B2p0_Chyd_Dcrt/gbr4_bgc_simple_2010-01.nc",
                          input_grid = NA,
                          eta_stem = NA,
                          override_positive=FALSE,
                          verbosity = 1)
{
  get_ereefs_ts(var_names=var_names, geocoordinates=geocoodinates, layer=layer, 
                start_date=start_date, end_date=end_date, input_file=input_file,
                input_grid=input_grid, eta_stem = eta_stem, overrride_positive=override_positive,
                verbosity=1)
}

#' Extracts depth-integrated time-series of selected variables at specified locations from eReefs output files
#'
#' See also get_ereefs_ts() to extract from a specified layer instead of a depth-integrated value.
#' By Barbara Robson (AIMS).
#'
#' @return a data frame containing the dates and values of extracted variables
#' @param var_names either a single character value or a vector specifying the short names for variables that you 
#'        want from the netcdf file. Defaults to c('Chl_a_sum', 'TN').
#' @param geocoordinates is a vector containing the decimal latitude and longitude of the desired location. If 
#'        you want to specify an x-y grid coordinate instead of a latitude and longitude, you can: to do this, 
#'        is.integer(geocoordinates) must be TRUE. Defaults to c(-23.39189, 150.88852).
#' @param start_date date  for animation. Can either be a) a vector of integers in the format 
#'	c(yr, month, day); b) a date obtained e.g. from as.Date(); or c) a character string 
#'      formatted for input to as.Date(). Defaults to c(2010,12,31).
#' @param end date for animation. Can either be a) a vector of integers in the format 
#'	c(yr, month, day); b) a date obtained e.g. from as.Date(); or c) a character string 
#'      formatted for input to as.Date(). Defaults to c(2016,12,31).
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
#'      "http://dapds00.nci.org.au/thredds/dodsC/fx3/gbr4_hydro_all/gbr4_all_2016-09.nc".
#' @param eta_stem The URI or file location of the model output files that contains the surface elevation (eta), minus the
#'       date components of the filename in the case of GBR1 or GBR4 files, and ommitting the file extension, ".nc". Needed
#'       only if eta is not in the files indicated by input_stem (e.g. some GBR1 bgc files).
#' @param override_positive Reverse the value of the "positive" attribute of botz for BGC files, assuming that it is
#'       incorrect. Default FALSE
#' @param verbosity (Defailt 1) how much do you want to know about progress?
#' @param mass (Default FALSE) Set to true if you want the mass per square metre rather than the mean concentration over depth returned
#' @param date_format (Default "date"). Set to "chron" if you'd like the date returned in chron format.
#' @export
#' @examples
#' \dontrun{
#' get_ereefs_depth_integrated_ts(c('Chl_a_sum', 'NH4'), geocoordinates=data.frame(latitude=-23.39189, longitude=150.88852), start_date=c(2010,12,31),end_date=c(2011,1,5), input_file=2)
#' get_ereefs_depth_integrated_ts(c('salt', 'temp'), geocoordinates=data.frame(latitude=c(-23.39189, -23.3), longitude=c(150.88852, 150.8)), start_date=c(2011,1,2),end_date=c(2011,1,5))
#'}
get_ereefs_depth_integrated_ts <- function(var_names=c('Chl_a_sum', 'TN'), 
                        geocoordinates=c(-23.39189, 150.88852), 
		                    start_date = c(2010,12,31), 
		                    end_date = c(2016,12,31), 
                        input_file = "menu",
			                  input_grid = NA,
			                  eta_stem = NA,
			                  override_positive=FALSE,
                        verbosity = 1,
                        mass = FALSE,
                        date_format = "date")
{
  # Get parameter values and assign results from returned list to relevant variable names
  # This assigns input_file, ereefs_case, input_stem, start_date, end_date, start_tod, start_month, start_yr,
  # end_date, end_day, end_month, end_yr, mths, yrs, var_list, ereefs_origin and spatial_grid
  assignList(get_params(start_date, end_date, input_file, var_names))


  z_grid <- get_ereefs_grids(input_file, input_grid)[['z_grid']]

  # Make sure the geocoordinates are in the desired data.frame format. Warn if column names are absent or unclear.
  if (!("data.frame" %in% class(geocoordinates))) {
    # geocoordinates has been provided as an array/matrix. Coerce it into a data frame for consistency.
    geocoordinates <- data.frame(latitude = geocoordinates[,1], longitude = geocoordinates[,2])
  }
  ilat <- which((names(geocoordinates) == "latitude")|(names(geocoordinates) == "lat"))
  ilon <- which((names(geocoordinates) == "longitude")|(names(geocoordinates) == "lon"))
  if (length(ilat)) {
    geocoordinates <- data.frame(latitude = geocoordinates[, ilat], longitude = geocoordinates[, ilon])
  } else {
    warning("Assuming that the first column of geocoordinates is latitude and the second column is longitude")
  }
  names(geocoordinates)[1:2] <- c("latitude", "longitude")

  geocoordinates <- geocoordinates %>% dplyr::rowwise() %>% 
                                       dplyr::mutate(grid_index = which.min(earth.dist(latitude, longitude, spatial_grid$latitude, spatial_grid$longitude)))
  geocoordinates[, c("i", "j", "cell_latitude", "cell_longitude")] = spatial_grid[geocoordinates$grid_index, c("i", "j", "latitude", "longitude")]

  # Initialise
  results <- NULL

  nc <- tidync::tidync(input_file)
  nc <- nc %>% tidync::activate(var_names[1]) %>% tidync::activate(tidync::active(nc), select_var = var_names)
  ndims <- dim(nc %>% hyper_dims())[1]
  if (ndims < 4) stop(paste(var_name[1], 'is not 3D in this file. Either this file contains only surface data\n',
                            'or this is a 2D variable. Cannot extract values at a specified depth.'))

  zat <- ncmeta::nc_att(input_file, 'botz', 'positive')$value
  if (!is.null(zat$positive)) {
	  if (zat=="down") zsign <- -1 else zsign <- 1
  } else {
	  zsign <-1
  }
  botz <- nc %>% tidync::activate('botz') %>% 
                 tidync::hyper_filter(i = dplyr::between(i, min(geocoordinates$i), max(geocoordinates$i)),
                                      j = dplyr::between(j, min(geocoordinates$j), max(geocoordinates$j))) %>%
                 tidync::hyper_tibble('botz') %>% 
                 dplyr::mutate(botz = botz * zsign)

  # Loop through relevant daily or monthly eReefs files to extract the data
  if ((ereefs_case[2] == "1km")|(ereefs_case[2] == "4km")) {
    i <- 0
    mcount <- 0
    if (verbosity>0) pb <- txtProgressBar(min = 0, max = length(mths), style = 3)
    for (month in mths) {
      mcount <- mcount + 1
      yr <- yrs[mcount]
      if (mcount == 1) {
         from_day <- start_day
      } else {
         from_day <- 1
         start_tod <- 0
      }
      if ((start_yr==end_yr)&&(start_month==end_month)) {
         day_count <- end_day - start_day + 1
      } else if (mcount == 1) {
         day_count <- daysIn(as.Date(paste(yr, month, 1, sep='-'))) - start_day + 1
      } else if (mcount == (length(mths))) {
         day_count <- end_day
      } else {
         day_count <- daysIn(as.Date(paste(yr, month, 1, sep='-')))
      }
      if (ereefs_case[2] == '4km') { 
        fileslist <- 1
        input_file <- paste0(input_stem, format(as.Date(paste(yr, month, 1, sep="-")), '%Y-%m'), '.nc')
	      day_count <- day_count / as.numeric(ds[2]-ds[1])
        if ((ds[length(ds)] - ds[1]) > 31.5) warning('Filename looks like a monthly output file (i.e. contains two dashes) but file contains more than a month of data.')
        if ((ds[length(ds)] - ds[1]) < 27) warning('Filename looks like a monthly output file (i.e. contains two dashes) but file contains less than a month of data.')
        if(ds[2]==ds[1]) stop(paste('Error reading time from', input_file, '(t[2]==t[1])'))
        if (day_count > length(ds)) {
          warning(paste('end_date', end_date, 'is beyond available data. Ending at', ds[length(ds)]))
          day_count <- length(ds)
        }
      } else if (ereefs_case[2] == '1km') {
	      fileslist <- from_day:(from_day+day_count-1)
	      from_day <- 1
	      day_count <- 1
      } else stop("Shouldn't happen: ereefs_case not recognised")

      for (dcount in fileslist) {
        if (ereefs_case[2] == '1km') {
	        input_file <- paste0(input_stem, format(as.Date(paste(yr, month, dcount, sep="-")), '%Y-%m-%d'), '.nc')
        }
        nc <- tidync::tidync(input_file) %>% tidync::activate(var_names[1]) 
        nc <- nc %>% tidync::activate(tidync::active(nc), select_var = var_names)

        nc_eta <- tryCatch(nc %>% tidync::activate('eta'), finally = NULL)
        if (is.null(nc_eta)) {
          if (!is.na(eta_stem)) {
            if (ereefs_case[2] == '4km') {
              etafile  <- paste0(eta_stem, format(as.Date(paste(start_yr, start_month, 1, sep='-')), '%Y-%m'), 
			          '.nc') 
            } else if (ereefs_case[2] == '1km') {
              etafile <- paste0(eta_stem, format(as.Date(paste(start_yr, start_month, start_day, sep='-')), '%Y-%m-%d'), 
			          '.nc') 
            }
          } else {
            stop(paste('eta (surface elevation) does not exist in', input_file, '\n',
                       'eta is required to calculate depth below the free surface.\n',
                       'It may exist in a standard format or hydrodynamic model output that mirrors this eReefs run.\n',
                       'with results output on the same time step. If so, specify a value for eta_stem. Otherwise, \n',
                       'you can use get_ereefs_ts() to extract data at a depth below MSL instead.'))
          }
          nc_eta <- tidync::tidync(etafile) %>% tidync::activate('eta')
        } else {
          if (!is.na(eta_stem)) warning(paste('eta exists in', input_file, ' so eta_stem is ignored.'))
        }
        eta <- nc_eta %>% tidync::hyper_filter(i = dplyr::between(i, min(geocoordinates$i), max(geocoordinates$i)),
                                            j = dplyr::between(j, min(geocoordinates$j), max(geocoordinates$j)),
                                            time = dplyr::between(index, from_day, from_day+day_count-1)) %>%
                          tidync::hyper_tibble('eta')
        botz <- right_join(botz, eta, by=c('i','j'))
        botz$botz[is.na(botz$botz)] <- 1e22

        z <- array(z_grid[2:length(z_grid)], dim=c(length(z_grid)-1, length(eta$eta)))
        zm1 <- array(z_grid[1:(length(z_grid)-1)], dim=c(length(z_grid)-1, length(eta$eta)))
        eta2 <- t(array(eta$eta, dim=c(length(eta$eta), length(z_grid)-1)))
        botz2 <- t(array(botz$botz, dim=c(length(botz$botz), length(z_grid)-1)))

        dz <- 0 * z
        wet <- (eta2 > zm1) & (z > botz2)      # There is water in this layer
        bottom <- (wet & (zm1<botz2))          # The bottom intersects this layer
        subsurface <- (wet & (eta2 > z))           # This is a wet layer that is not the surface layer
        surface <- (wet & !subsurface)             # This is the surface layer
        singlelayer <- (surface & bottom)          # This is both the surface and bottom layer
        surface <- (surface & !bottom)             # This is the surface layer but not the bottom layer
        fullywet <- (subsurface & !bottom)         # This is a wet layer that isn't the surface or bottom layer
        submergedbottom <- (bottom & !singlelayer) # This is the bottom layer but not the surface layer
  
        # Thickness of each layer
        dz[surface] <- eta2[c(surface)] - zm1[c(surface)]
        dz[fullywet] <- z[c(fullywet)] - zm1[c(fullywet)]
        dz[submergedbottom] <- z[c(submergedbottom)] - botz2[c(submergedbottom)]
        dz[singlelayer] <- eta2[c(singlelayer)] - botz2[c(singlelayer)]

        row.names(dz) <- paste0('k', 1:dim(dz)[1])
        dz <- cbind(dplyr::as_tibble(t(dz)), botz) %>% dplyr::mutate(watercol_depth = -(botz + eta)) %>% dplyr::select(-botz, -eta)
        dz$watercol_depth[dz$watercol_depth>1e21] = 0
        dz <- dz %>% tidyr::pivot_longer(cols = starts_with("k"), names_to="k", names_prefix="k", values_to="dz") %>%
          dplyr::mutate(k = as.integer(k))
  
        wc <- nc %>% tidync::hyper_filter(i = dplyr::between(i, min(geocoordinates$i), max(geocoordinates$i)),
                                          j = dplyr::between(j, min(geocoordinates$j), max(geocoordinates$j)),
                                          time = dplyr::between(index,from_day, from_day+day_count-1)) %>%
                     tidync::hyper_tibble(select_var = var_names)
        wc <- dplyr::left_join(wc, dz, by=c('i','j','time','k')) 
        results <- rbind(results, wc)
      }
    }
  } else {
        nc_eta <- tryCatch(nc %>% tidync::activate('eta'), finally = NULL)
        if (is.null(nc_eta)) {
          if (!is.na(eta_stem)) {
            etafile <- paste0(eta_stem, '.nc')
            nc_eta <- tidync::tidync(etafile) %>% tidync::activate('eta')
          } else {
            stop(paste('eta (surface elevation) does not exist in', input_file, '\n',
                       'eta is required to calculate depth below the free surface.\n',
                       'It may exist in a standard format or hydrodynamic model output that mirrors this eReefs run.\n',
                       'with results output on the same time step. If so, specify a value for eta_stem. Otherwise, \n',
                       'you can use get_ereefs_ts() to extract data at a depth below MSL instead.'))
          }
        } else {
          if (!is.na(eta_stem)) warning(paste('eta exists in', input_file, ' so eta_stem is ignored.'))
        }
        eta <- nc_eta %>% tidync::hyper_filter(i = dplyr::between(i, min(geocoordinates$i), max(geocoordinates$i)),
                                            j = dplyr::between(j, min(geocoordinates$j), max(geocoordinates$j)),
                                            time = dplyr::between(time, as.numeric(start_date - ereefs_origin), as.numeric(end_date - ereefs_origin))) %>%
                          tidync::hyper_tibble('eta')
        botz <- right_join(botz, eta, by=c('i','j'))
        botz$botz[is.na(botz$botz)] <- 1e22

        z <- array(z_grid[2:length(z_grid)], dim=c(length(z_grid)-1, length(eta$eta)))
        zm1 <- array(z_grid[1:(length(z_grid)-1)], dim=c(length(z_grid)-1, length(eta$eta)))
        eta2 <- t(array(eta$eta, dim=c(length(eta$eta), length(z_grid)-1)))
        botz2 <- t(array(botz$botz, dim=c(length(botz$botz), length(z_grid)-1)))

        dz <- 0 * z
        wet <- (eta2 > zm1) & (z > botz2)      # There is water in this layer
        bottom <- (wet & (zm1<botz2))          # The bottom intersects this layer
        subsurface <- (wet & (eta2 > z))           # This is a wet layer that is not the surface layer
        surface <- (wet & !subsurface)             # This is the surface layer
        singlelayer <- (surface & bottom)          # This is both the surface and bottom layer
        surface <- (surface & !bottom)             # This is the surface layer but not the bottom layer
        fullywet <- (subsurface & !bottom)         # This is a wet layer that isn't the surface or bottom layer
        submergedbottom <- (bottom & !singlelayer) # This is the bottom layer but not the surface layer
  
        # Thickness of each layer
        dz[surface] <- eta2[c(surface)] - zm1[c(surface)]
        dz[fullywet] <- z[c(fullywet)] - zm1[c(fullywet)]
        dz[submergedbottom] <- z[c(submergedbottom)] - botz2[c(submergedbottom)]
        dz[singlelayer] <- eta2[c(singlelayer)] - botz2[c(singlelayer)]

        row.names(dz) <- paste0('k', 1:dim(dz)[1])
        dz <- cbind(dplyr::as_tibble(t(dz)), botz) %>% dplyr::mutate(watercol_depth = -(botz + eta)) %>% dplyr::select(-botz, -eta)
        dz$watercol_depth[dz$watercol_depth>1e21] = 0
        dz <- dz %>% tidyr::pivot_longer(cols = starts_with("k"), names_to="k", names_prefix="k", values_to="dz") %>%
          dplyr::mutate(k = as.integer(k))
  
        wc <- nc %>% tidync::hyper_filter(i = dplyr::between(i, min(geocoordinates$i), max(geocoordinates$i)),
                                          j = dplyr::between(j, min(geocoordinates$j), max(geocoordinates$j)),
                                          time = dplyr::between(time, as.numeric(start_date - ereefs_origin), as.numeric(end_date - ereefs_origin))) %>%
                     tidync::hyper_tibble(select_var = var_names)
        wc <- dplyr::left_join(wc, dz, by=c('i','j','time','k')) 
        results <- rbind(results, wc)
  }

  var_summary <- function(data, var) {
    data %>% dplyr::group_by(i,j,time) %>% 
      dplyr::summarise(i =i, j=j, watercol_depth = watercol_depth, "{var}" := sum(.data[[var]] * dz / watercol_depth))
  }

  results <- results %>% var_summary(var_names[j])
  for (j in 2:length(var_names)) { 
    results <- dplyr::left_join(wc %>% var_summary(var_names[j]), results, by=c('i','j','time','watercol_depth'))
  }
  results <- dplyr::left_join(results, geocoordinates, by=c('i','j'))
  results <- results %>% mutate(time = ereefs_origin + seconds(time*86400)) %>% select(-grid_index)
  return(results)
}

#' Extracts time-series of selected variables from eReefs output files at a specified location and depth below the
#' tidal free surface. Use get_ereefs_ts() instead for a specified depth below MSL: this is faster and given the 
#' thickness of the layers, usually almost the same.
#'
#' See also get_ereefs_ts() to extract from a specified layer instead of a depth-integrated value and 
#' get_ereefs_depth_integrated_ts() to calculate depth-integrated values. Barbara Robson (AIMS).
#'
#' @return a data frame containing the dates and values of extracted variables
#' @param var_names either a single character value or a vector specifying the short names for variables that you 
#'        want from the netcdf file. Defaults to c('Chl_a_sum', 'TN').
#' @param geocoordinates A vector containing the decimal latitude and longitude of the desired location. If 
#'        you want to specify an x-y grid coordinate instead of a latitude and longitude, you can: to do this, 
#'        is.integer(geocoordinates) must be TRUE. Defaults to c(-23.39189, 150.88852).
#' @param depth  Depth in metres below the surface. Default 1.0. If the bottom of the water is shallower than the specified depth,
#'               return values from the bottom of the water column.
#' @param start_date date  for animation. Can either be a) a vector of integers in the format 
#'	c(yr, month, day); b) a date obtained e.g. from as.Date(); or c) a character string 
#'      formatted for input to as.Date(). Defaults to c(2010,12,31).
#' @param end date for animation. Can either be a) a vector of integers in the format 
#'	c(yr, month, day); b) a date obtained e.g. from as.Date(); or c) a character string 
#'      formatted for input to as.Date(). Defaults to c(2016,12,31).
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
#' @param eta_stem The URI or file location of the model output files that contains the surface elevation (eta), minus the
#'       date components of the filename in the case of GBR1 or GBR4 files, and ommitting the file extension, ".nc". Needed
#'       only if eta is not in the files indicated by input_stem (e.g. some GBR1 bgc files).
#' @export
#' @examples
#' \dontrun{
#' get_ereefs_depth_specified_ts(c('Chl_a_sum', 'NH4'), depth=5.0, geocoordinates=data.frame(latitude=-23.39189, longitude=150.88852), 
#'   start_date=c(2010,12,31),end_date=c(2011,1,5), input_file=2)
#' get_ereefs_depth_specified_ts('salt', depth=1.0, geocoordinates=data.frame(latitude=-23.39189, longitude=150.88852), 
#'   start_date=c(2010,12,31),end_date=c(2011,1,5))
#'}
get_ereefs_depth_specified_ts <- function(var_names=c('Chl_a_sum', 'TN'), 
                                          geocoordinates=c(-23.39189, 150.88852), 
                                          depth = 1.0, 
                                          start_date = c(2010,12,31), 
                                          end_date = c(2016,12,31), 
                                          input_file = "menu", 
                                          input_grid = NA, 
                                          eta_stem = NA, 
                                          verbosity = 1,
                                          date_format = "date")
{
  # Get parameter values and assign results from returned list to relevant variable names
  # This assigns input_file, ereefs_case, input_stem, start_date, end_date, start_tod, start_month, start_yr,
  # end_date, end_day, end_month, end_yr, mths, yrs, var_list, ereefs_origin and spatial_grid
  assignList(get_params(start_date, end_date, input_file, var_names))


  z_grid <- get_ereefs_grids(input_file, input_grid)[['z_grid']]

  # Make sure the geocoordinates are in the desired data.frame format. Warn if column names are absent or unclear.
  if (!("data.frame" %in% class(geocoordinates))) {
    # geocoordinates has been provided as an array/matrix. Coerce it into a data frame for consistency.
    geocoordinates <- data.frame(latitude = geocoordinates[,1], longitude = geocoordinates[,2])
  }
  ilat <- which((names(geocoordinates) == "latitude")|(names(geocoordinates) == "lat"))
  ilon <- which((names(geocoordinates) == "longitude")|(names(geocoordinates) == "lon"))
  if (length(ilat)) {
    geocoordinates <- data.frame(latitude = geocoordinates[, ilat], longitude = geocoordinates[, ilon])
  } else {
    warning("Assuming that the first column of geocoordinates is latitude and the second column is longitude")
  }
  names(geocoordinates)[1:2] <- c("latitude", "longitude")

  geocoordinates <- geocoordinates %>% dplyr::rowwise() %>% 
                                       dplyr::mutate(grid_index = which.min(earth.dist(latitude, longitude, spatial_grid$latitude, spatial_grid$longitude)))
  geocoordinates[, c("i", "j", "cell_latitude", "cell_longitude")] = spatial_grid[geocoordinates$grid_index, c("i", "j", "latitude", "longitude")]

  # Initialise
  results <- NULL

  nc <- tidync::tidync(input_file)
  nc <- nc %>% tidync::activate(var_names[1]) %>% tidync::activate(tidync::active(nc), select_var = var_names)
  ndims <- dim(nc %>% hyper_dims())[1]
  if (ndims < 4) stop(paste(var_name[1], 'is not 3D in this file. Either this file contains only surface data\n',
                            'or this is a 2D variable. Cannot extract values at a specified depth.'))

  zat <- ncmeta::nc_att(input_file, 'botz', 'positive')$value
  if (!is.null(zat$positive)) {
	  if (zat=="down") zsign <- -1 else zsign <- 1
  } else {
	  zsign <-1
  }
  botz <- nc %>% tidync::activate('botz') %>% 
                 tidync::hyper_filter(i = dplyr::between(i, min(geocoordinates$i), max(geocoordinates$i)),
                                      j = dplyr::between(j, min(geocoordinates$j), max(geocoordinates$j))) %>%
                 tidync::hyper_tibble('botz') %>% 
                 dplyr::mutate(botz = botz * zsign)

  # Loop through relevant daily or monthly eReefs files to extract the data
  if ((ereefs_case[2] == "1km")|(ereefs_case[2] == "4km")) {
    i <- 0
    mcount <- 0
    if (verbosity>0) pb <- txtProgressBar(min = 0, max = length(mths), style = 3)
    for (month in mths) {
      mcount <- mcount + 1
      yr <- yrs[mcount]
      if (mcount == 1) {
         from_day <- start_day
      } else {
         from_day <- 1
         start_tod <- 0
      }
      if ((start_yr==end_yr)&&(start_month==end_month)) {
         day_count <- end_day - start_day + 1
      } else if (mcount == 1) {
         day_count <- daysIn(as.Date(paste(yr, month, 1, sep='-'))) - start_day + 1
      } else if (mcount == (length(mths))) {
         day_count <- end_day
      } else {
         day_count <- daysIn(as.Date(paste(yr, month, 1, sep='-')))
      }
      if (ereefs_case[2] == '4km') { 
        fileslist <- 1
        input_file <- paste0(input_stem, format(as.Date(paste(yr, month, 1, sep="-")), '%Y-%m'), '.nc')
	      day_count <- day_count / as.numeric(ds[2]-ds[1])
        if ((ds[length(ds)] - ds[1]) > 31.5) warning('Filename looks like a monthly output file (i.e. contains two dashes) but file contains more than a month of data.')
        if ((ds[length(ds)] - ds[1]) < 27) warning('Filename looks like a monthly output file (i.e. contains two dashes) but file contains less than a month of data.')
        if(ds[2]==ds[1]) stop(paste('Error reading time from', input_file, '(t[2]==t[1])'))
        if (day_count > length(ds)) {
          warning(paste('end_date', end_date, 'is beyond available data. Ending at', ds[length(ds)]))
          day_count <- length(ds)
        }
      } else if (ereefs_case[2] == '1km') {
	      fileslist <- from_day:(from_day+day_count-1)
	      from_day <- 1
	      day_count <- 1
      } else stop("Shouldn't happen: ereefs_case not recognised")

      for (dcount in fileslist) {
        if (ereefs_case[2] == '1km') {
	        input_file <- paste0(input_stem, format(as.Date(paste(yr, month, dcount, sep="-")), '%Y-%m-%d'), '.nc')
        }
        nc <- tidync::tidync(input_file) %>% tidync::activate(var_names[1]) 
        nc <- nc %>% tidync::activate(tidync::active(nc), select_var = var_names)

        nc_eta <- tryCatch(nc %>% tidync::activate('eta'), finally = NULL)
        if (is.null(nc_eta)) {
          if (!is.na(eta_stem)) {
            if (ereefs_case[2] == '4km') {
              etafile  <- paste0(eta_stem, format(as.Date(paste(start_yr, start_month, 1, sep='-')), '%Y-%m'), 
			          '.nc') 
            } else if (ereefs_case[2] == '1km') {
              etafile <- paste0(eta_stem, format(as.Date(paste(start_yr, start_month, start_day, sep='-')), '%Y-%m-%d'), 
			          '.nc') 
            } else {
              etafile <- paste0(eta_stem, '.nc')
            }
          } else {
            stop(paste('eta (surface elevation) does not exist in', input_file, '\n',
                       'eta is required to calculate depth below the free surface.\n',
                       'It may exist in a standard format or hydrodynamic model output that mirrors this eReefs run.\n',
                       'with results output on the same time step. If so, specify a value for eta_stem. Otherwise, \n',
                       'you can use get_ereefs_ts() to extract data at a depth below MSL instead.'))
          }
          nc_eta <- tidync::tidync(etafile) %>% tidync::activate('eta')
        } else {
          if (!is.na(eta_stem)) warning(paste('eta exists in', input_file, ' so eta_stem is ignored.'))
        }
        eta <- nc_eta %>% tidync::hyper_filter(i = dplyr::between(i, min(geocoordinates$i), max(geocoordinates$i)),
                                            j = dplyr::between(j, min(geocoordinates$j), max(geocoordinates$j)),
                                            time = dplyr::between(index, from_day, from_day+day_count-1)) %>%
                          tidync::hyper_tibble('eta')

        z <- array(z_grid[2:length(z_grid)], dim=c(length(z_grid)-1, length(eta$eta)))
        zm1 <- array(z_grid[1:(length(z_grid)-1)], dim=c(length(z_grid)-1, length(eta$eta)))
        eta2 <- t(array(eta$eta, dim=c(length(eta$eta), length(z_grid)-1)))

        # Depth relative to surface
        zr <- z - eta2
        zm1r <- zm1 - eta2
        dz <- 0 * z
        wet <- (eta2 > zm1) & (z > botz$botz)        # There is water in this layer
        bottom <- (wet & (zm1<botz$botz))            # The bottom intersects this layer

        shallower_top <- (-zr < depth)               # The top of this layer is above the target depth
        deeper_bottom <- (-zm1r > depth) | (bottom)  # The bottom of this layer is below the target depth or this is the bottom layer
        target <- (shallower_top & deeper_bottom)    
   
        count <- 0
        for (ind in from_day:(from_day + day_count - 1)) {
          count <- count + 1
          wc <- nc %>% tidync::hyper_filter(i = dplyr::between(i, min(geocoordinates$i), max(geocoordinates$i)),
                                            j = dplyr::between(j, min(geocoordinates$j), max(geocoordinates$j)),
                                            time = index==ind,
                                            k = target[,count]) %>%
                     tidync::hyper_tibble(select_var = var_names)
          wc <- dplyr::left_join(geocoordinates, wc, by = c("i", "j"))
          results <- rbind(results, wc)
        }
        if (verbosity>0) setTxtProgressBar(pb,mcount)
      }
    }
    if (verbosity>0) close(pb)
  } else {
    # Single nc or ncml file. No need to loop through months, but might need to check for monotonicity of timeseries 
    print(paste('No waitbar shown in this case. If you have asked for a large dataset, please be patient.',
                  'If you have asked for a long time-series from a THREDDS server,',
                  'you might encounter transfer errors, but it should self-correct.',
                  'If not, consider breaking it up into shorter sections or pointing to an individual netcdf file.'))
    if (end_date > max(ds, na.rm=TRUE)) {
      warning(paste('end_date', end_date, 'is beyond available data. Ending at', max(ds, na.rm=TRUE)))
      end_date <- max(ds, na.rm=TRUE)
    }
    if (start_date < min(ds, na.rm=TRUE)) {
      warning(paste('start_date', start_date, 'is beyond available data. Starting at', min(ds, na.rm=TRUE)))
      start_date <- min(ds, na.rm=TRUE)
    }

    nc_eta <- tryCatch(nc %>% tidync::activate('eta'), finally = NULL)
    if (is.null(nc_eta)) {
      if (!is.na(eta_stem)) {
        etafile <- paste0(eta_stem, '.nc')
      } else {
        stop(paste('eta (surface elevation) does not exist in', input_file, '\n',
                   'eta is required to calculate depth below the free surface.\n',
                   'It may exist in a standard format or hydrodynamic model output that mirrors this eReefs run.\n',
                   'with results output on the same time step. If so, specify a value for eta_stem. Otherwise, \n',
                   'you can use get_ereefs_ts() to extract data at a depth below MSL instead.'))
      }
      nc_eta <- tidync::tidync(etafile) %>% tidync::activate('eta')
    } else {
      if (!is.na(eta_stem)) warning(paste('eta exists in', input_file, ' so eta_stem is ignored.'))
    }
    eta <- nc_eta %>% tidync::hyper_filter(i = dplyr::between(i, min(geocoordinates$i), max(geocoordinates$i)),
                                        j = dplyr::between(j, min(geocoordinates$j), max(geocoordinates$j)),
                                        time = dplyr::between(time, as.numeric(start_date - ereefs_origin), as.numeric(end_date - ereefs_origin)))  %>%
                      tidync::hyper_tibble('eta')

    z <- array(z_grid[2:length(z_grid)], dim=c(length(z_grid)-1, length(eta$eta)))
    zm1 <- array(z_grid[1:(length(z_grid)-1)], dim=c(length(z_grid)-1, length(eta$eta)))
    eta2 <- t(array(eta$eta, dim=c(length(eta$eta), length(z_grid)-1)))

    # Depth relative to surface
    zr <- z - eta2
    zm1r <- zm1 - eta2
    dz <- 0 * z
    wet <- (eta2 > zm1) & (z > botz$botz)           # There is water in this layer
    bottom <- (wet & (zm1<botz$botz))               # The bottom intersects this layer

    shallower_top <- (-zr < depth)               # The top of this layer is above the target depth
    deeper_bottom <- (-zm1r > depth) | (bottom)  # The bottom of this layer is below the target depth or this is the bottom layer
    target <- (shallower_top & deeper_bottom)    
   
    nc <- tidync::tidync(input_file) %>% tidync::activate(var_names[1]) %>% tidync::activate(tidync::active(nc), select_var = var_names)

    
    count <- 0
    for (ind in which(dplyr::between(ds, start_date, end_date))) { 
      count <- count + 1
      wc <- nc %>% tidync::hyper_filter(i = dplyr::between(i, min(geocoordinates$i), max(geocoordinates$i)),
                                        j = dplyr::between(j, min(geocoordinates$j), max(geocoordinates$j)),
                                        k = target[, count], 
                                        time = index==ind) %>%
                 tidync::hyper_tibble(select_var = var_names)
      wc <- dplyr::left_join(geocoordinates, wc, by = c("i", "j"))
      results <- rbind(results, wc)
    }
  }
  results <- results %>% mutate(time = ereefs_origin + seconds(time*86400))
  results <- results %>% select(-grid_index)
  return(results)
}

#' A wrapper to ncdf4::ncvar_get() that will pause and try again several times (defaulting to 12)
#' if it is a web-served netcdf file and at first it fails, to overcome temporary net access errors or DAP errors.
#'
#' Parameters before 'tries' are passed unlist(month.day.yr(start_date))
#'
#' @param tries number of times to retry (increasing pause length by one second each time. Default 4 
#' @return variable extracted using ncvar_get()
#' @export
safe_ncvar_get <- function(nc,varid=NA, start=NA, count=NA, verbose=FALSE,
 signedbyte=TRUE, collapse_degen=TRUE, raw_datavals=FALSE, tries=20) {
   if (substr(nc$filename, 1, 4)!="http") {
     myvar <- ncdf4::ncvar_get(nc, varid, start, count, verbose, signedbyte, collapse_degen, raw_datavals)
   } else {
     myvar <- try(ncdf4::ncvar_get(nc, varid, start, count, verbose, signedbyte, collapse_degen, raw_datavals))
     trywait = 1
     while ((class(myvar)[[1]]=='try-error')&(trywait<=(tries*2))) { 
        print(paste('retrying in ', trywait, 'second(s)')) 
        Sys.sleep(trywait) 
        trywait <- trywait+1 
        myvar <- try(ncdf4::ncvar_get(nc, varid, start, count, verbose, signedbyte, collapse_degen, raw_datavals))
     }
     if (trywait>(tries*2)) stop(paste('Cannot access netcdf file', nc$filename))
   }
   return(myvar)
}

#' A wrapper to ncdf4::nc_open() that will pause and try again up to 119 times
#' if at first it fails, to overcome temporary net access errors or DAP errors.
#'
#' Parameters before 'tries' are passed through to ncvar_get
#'
#' @param tries number of times to retry (increasing pause length by one second each time. Default 4 
#' @return variable extracted using ncvar_get()
#' @export
safe_nc_open <- function(filename, tries=4) {
   
   if (substr(filename, 1, 4)!="http") {
    nc <- ncdf4::nc_open(filename)
   } else {
    nc <- try(ncdf4::nc_open(filename))
    trywait = 1
    while ((class(nc)[[1]]=='try-error')&(trywait<=(tries+1))) { 
       warning(paste('Trouble opening', filename,
              'This probably means the netcdf file name is incorrect or target date out of range for this eReefs run,\n',
              'but could be a network connection issue...\n', 
              '  Retrying in ', trywait, 'second(s). This will be attempt', trywait+1, 'of', tries+1)) 
       Sys.sleep(trywait) 
       trywait <- trywait+1 
       nc <- try(ncdf4::nc_open(filename))
    }
    if (trywait>(tries+1)) stop(paste('Cannot open netcdf file', filename))
   }
   return(nc)
}

#' A simple function to convert a date provided in any of several formats to a date-time object
#' @param d The date of interest. Can be any of:
#'              c(yr, month, day)
#'              c(yr, month, day, hour) 
#'              c(yr, month, day, hour, minute)
#'              c(yr, month, day, hour, minute, second) 
#'              YYYYMMDD
#'              YYYYMMDDHHMM
#'              YYYYMMDDHHMMSS
#'              Date format date (e.g., as.Date('1970-01-01', origin='1970-01-01'))
#'              character format, e.g. '1970-01-01'
#' @return date in date-time format (for lubridate)
#' @export
get_date_time <- function(d) {
  if (is.vector(d)) {
    if (length(d)==1) {
      d <- switch(nchar(d)/2 - 2, 
                  lubridate::ymd(d, tz = "Etc/GMT-10"),
                  lubridate::ymd_h(d, tz = "Etc/GMT-10"),
                  lubridate::ymd_hm(d, tz = "Etc/GMT-10"),
                  lubridate::ymd_hms(d, tz = "Etc/GMT-10"))
    } else {
      if (length(d)<4) d[4] <- 12 # midday
      if (length(d)<5) d[5] <- 00
      if (length(d)<6) d[6] <- 00
      d <- lubridate::ymd_hms(paste0(d[1], '-', d[2], '-', d[3], ' ', d[4], ":", d[5], ":", d[6]), tz="Etc/GMT-10")
    }
  } else if (is.character(d)) {
    d <- lubridate::ymd_h(paste(d, '12'), tz="Etc/GMT-10")
  } else if (class(d)[1] == "Date") {
    d <- lubridate::ymd_h(paste(format(d, '%Y-%m-%d'), '12'), tz="Etc/GMT-10")
  }
}
