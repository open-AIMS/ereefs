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
#' @return 1 for daily, 4 for monthly, 0 for other
#' @export
get_ereefs_case <- function(filename) {
  # If the filename ends in ".nc.html", truncate to remove the ".html" part
  if (stringr::str_ends(filename, stringr::fixed(".nc.html"))) filename <- stringr::str_sub(filename, 1, -6)

  if (stringr::str_ends(filename, stringr::fixed(".nc"))) {
    ereefs_case <- "nc"
  } else if (stringr::str_ends(filename, stringr::fixed(".mnc"))) {
    ereefs_case <- "mnc"
    warn(".mnc files not yet implemented.")
  } else if (stringr::str_ends(filename, "catalog.html")) {
    ereefs_case <- c("thredds_catalog", "nci")
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

#' Opens the netcdf file, input_file, and extracts the origin (reference date/time) and time dimension.
#'
#' @param input_file The name of the netcdf file (in standard or simple EMS netcdf format) from which
#' to extract the data
#' @return A list containing ereefs_origin (the reference data/time used in the netcdf file) and the time series, ds
get_origin_and_times <- function(input_file, as_chron="TRUE") {
	  nc <- safe_nc_open(input_file)
    if (!is.null(nc$var[['t']])) { 
      posix_origin <- stringi::stri_datetime_parse(ncdf4::ncatt_get(nc ,'t'), "'days since 'yyyy-MM-dd HH:mm:ss")[1]
      ereefs_origin <- as.numeric(as.Date('1990-01-01') - as.Date(posix_origin), origin=c(year=1990, month=1, day=1))-1 
      ds <- safe_ncvar_get(nc, "t") 
    } else { 
      posix_origin <- stringi::stri_datetime_parse(ncdf4::ncatt_get(nc ,'time'), "'days since 'yyyy-MM-dd HH:mm:ss")[1]
      ereefs_origin <- as.numeric(as.Date('1990-01-01') - as.Date(posix_origin), origin=c(year=1990, month=1, day=1))-1 
      ds <- safe_ncvar_get(nc, "time") 
    }
    ncdf4::nc_close(nc) 
    if (as_chron) {
      ereefs_origin <- ereefs_origin + chron::chron('1990-01-01', origin=c(year=1990, month=1, day=1), format='y-m-d')
      ds <- ds + ereefs_origin
    } else {
      ds <- as.Date(ds, origin = as.Date("1990-01-01"))
    }
    return(list(ereefs_origin, ds))
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
		nc <- safe_nc_open(input_grid)
		x_grid <- safe_ncvar_get(nc, 'x_grid')
		y_grid <- safe_ncvar_get(nc, 'y_grid')
		z_grid <- safe_ncvar_get(nc, 'z_grid')
		ncdf4::nc_close(nc)
	} else {
		nc <- safe_nc_open(filename)
		if (!is.null(nc$var[['x_grid']])) { 
      # filename is in standard EMS netcdf output format. x_grid, y_grid and z_grid are provided.
		  x_grid <- safe_ncvar_get(nc, 'x_grid')
		  y_grid <- safe_ncvar_get(nc, 'y_grid')
		  z_grid <- safe_ncvar_get(nc, 'z_grid')
    } else {
      # simple format netcdf file. grids are not provided. Check the dimensions to work out 
      # whether it looks like GBR4 or GBR1 and if so, load x-grid and y_grid from saved data. 
      # Otherwise, the user must provide another filename (input_grid) that contains x_grid, 
      # y_grid and z_grid. 
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
		    stop(paste('x_grid, y_grid and z_grid not found and grid size not recognised as GBR1 or GBR4.',
                   'Plotting using cell centres (latitude and longitude) has not yet been implemented.',
                   'Please specify input_grid (a standard format EMS file or other netcdf file that ',
                   'contains x_grid, y_grid and z_grid).'))

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
  if (case[1]!='nc') {
    return(NA)
  }
  if (is.na(case[2]) | case[2] == "recom") {
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

#' Check whether the platform is Windows and the filename contains "http": if so, return a warning
#' Now always returns TRUE as the CRAN version of ncdf4 for Windows can now handle THREDDS services
#'
#' Utility function for the eReefs package
#' @param input_stem A netcdf filename or file stem
#' @return TRUE or stop error
#' @export
check_platform_ok <- function(input_stem)
{
  ok <- TRUE
  #webserved <- stringi::stri_startswith_fixed(input_stem, "http:")
  #if ((.Platform$OS.type=="windows")&&(webserved)) {
     #ok <- FALSE
     #if (packageVersion("ncdf4")!="1.16.9001") {
        #warning("Please check that you have installed the mdsumner branch of the ncdf4 package. This package will not work reliably under Windows using the CRAN version of ncdf4.")
        #warning("To install the mdsumner version, use: install('devtools'); devtools::install_github('mdsumner/ncdf4')")
     #}
  #}
  return(ok)
}

#' Check whether the given filename is a shortcut and if so, set the full filename
#'
#' Utility function for the eReefs package
#' @param input_file A netcdf filename, file stem or THREDDS catalog URI, or "menu" or "old_menu", or a run label from ereefs.info such as "GBR4-v2.0"
#'                   Default is "menu".
#' @return input_file
#' @export
substitute_filename <- function(input_file = "menu") {
  choices <- NA
  if ((input_file == "menu")|(input_file == "old_menu")|(input_file == "choices") | is.numeric(input_file)) {
    # Set up the menu of choices to display for interactive input, or to select from if the user has specified an option number
    if (input_file == "old_menu") {
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
                    "GBR4_BGC-v3.0 Dcrt",
                    #"GBR4_BGC-v3.0 Dnrt",
                    "GBR4_BGC-v3.1",
                    "GBR4_BGC-v3.0",
                    "GBR1_BGC-v3p2surf",
                    "GBR1_BGC-v3p2",
                    "GBR4_BGC-v3p2nrtsurf",
                    "GBR4_BGC-v3p2nrt",
                    "menu")
    } else {
      # Get the list of eReefs data sevices from the NCI server:
      services <- thredds::tds_list_datasets("https://dapds00.nci.org.au/thredds/catalogs/fx3/catalog.html")
      # Trim the list to only show the catalogues
      services <- services[which(services$type=="catalog"), ]
      choices <- services$dataset
      #choices[length(choices) + 1] <- "menu"
      # I'm probably missing something, but the following returns paths that will work:
      paths <- stringr::str_replace(services$path, "catalogs/fx3//thredds/", "") 
    }
  }
  if ((stringr::str_ends(input_file, "catalog.html"))|((stringr::str_starts(input_file, "http")&(stringr::str_ends(input_file, "/"))))) {
      services <- thredds::tds_list_datasets(input_file)
      # Trim the list to only show the catalogues
      services <- services[which(services$type=="catalog"), ]
      choices <- services$dataset
      #choices[length(choices) + 1] <- "menu"
      # I'm probably missing something, but the following returns paths that will work:
      paths <- stringr::str_replace(services$path, "catalogs/fx3//thredds/", "") 
  } else if (input_file == "choices") {
  # The user just wants a list of options
   print(choices)
   stop()
  } else if (input_file=="old_menu") {
      selection <- utils::menu(c("Latest release 4km grid hydrodynamic model (Sept 2010-pres.)", 
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
                               "CSIRO login required: Latest release GBR4 NRT BGC 3D (Oct 2019 - May 2020, daily)"
                              ))
    input_file <- choices[selection]
  } else if ((input_file == "menu")|(input_file == length(choices))) {
    # We are using the NCI catalog to provide options
    print("Refer to https://research.csiro.au/ereefs/models/models-about/biogeochemical-simulation-naming-protocol/ for naming conventions.")
    selection <- utils::menu(choices)
    input_file  <- paths[selection]
  }
  if (is.numeric(input_file)) {
     input_file <- choices[input_file]
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
      TRUE ~ input_file
    )
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
#' If you run into memory constraints, consider grouping points to be extracted within regions, and calling this once
#' for each region.
#'
#' @return a data frame containing the dates and values of extracted variables and location(s) or a
#'        list of such data frames (one list item per location) if return_list=TRUE (for backward compatibility).
#' @param var_names either a single character value or a vector specifying the short names for variables that you 
#'        want from the netcdf file. Defaults to c('Chl_a_sum', 'TN').
#' @param location_latlon is a data frame containing the decimal latitude and longitude of a single desired location, or a vector containing
#'        a single latitude and longitude location. If you want to specify an x-y grid coordinate instead of a latitude and longitude, you 
#'        can: to do this, is.integer(location_latlon) must be TRUE. location_latlon can also be set to "mmp" to extract time-series at all 
#'        Marine Monitoring Program sites. Defaults to c(-23.39189, 150.88852).
#' @param layer is the vertical grid layer to extract, or 'surface' to get the surface value, 'bottom' to get the
#'        value in the cell at the bottom of the water column, or 'integrated' to get a depth-integrated (mean) value.
#'        Specify a negative value to indicate a specified depth (in metres) below MSL.
#'        Use get_ereefs_depth_specified_ts() instead if you want to specify a depth 
#'        below the free surface instead of a layer number. Defaults to 'surface'.
#' @param start_date date  for animation. Can either be a) a vector of integers in the format 
#'	       c(year, month, day); b) a date obtained e.g. from as.Date(); or c) a character string 
#'        formatted for input to as.Date(). Defaults to c(2010,12,31).
#' @param end date for animation. Can either be a) a vector of integers in the format 
#'	       c(year, month, day); b) a date obtained e.g. from as.Date(); or c) a character string 
#'        formatted for input to as.Date(). Defaults to c(2016,10,31).
#' @param input_file is the URI or file location of any of the EMS output files or a THREDDS catalog URI. 
#'        Defaults to a menu selection. Set to "choices" to see some other pre-defined options that
#'        can be used (codenames as used in https://research.csiro.au/ereefs/models/model-outputs/access-to-raw-model-output/ )
#'        Set to "old_menu" to provide old menu options instead of menu options from the NCI catalog.
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
#' @param return_list Default FALSE. Set to true if you want a list of dataframes returned (one df per geolocation). Included for
#'                    backward compatibility.
#' @param date_format Default "date". Set to "chron" to have the date returned in chron format.
#' @export
#' @examples
#' \dontrun{
#' get_ereefs_ts('Chl_a_sum', location_latlon=data.frame(latitide=c(-23.39189,-18), longitude=c(150.88852, 147.5)), layer='surface', start_date=c(2010,12,31),end_date=c(2011,1,5), input_file=2)
#' get_ereefs_ts(var_names=c('Tricho_N', 'DIP', 'TP'), location_latlon=data.frame(latitide=c(-23.39189,-18), longitude=c(150.88852, 147.5)), layer='bottom', start_date="2012-07-03",end_date="2012-07-05", input_file='GBR4_BGC-v2.0 Chyb Dcrt')
#' get_ereefs_ts(var_names=c('ZooL_N', 'ZooS_N'), location_latlon=data.frame(latitide=c(-23.39189,-18), longitude=c(150.88852, 147.5)), layer=45, start_date=c(2010,12,31),end_date=c(2011,1,5), input_file="http://dapds00.nci.org.au/thredds/dodsC/fx3/gbr4_bgc_GBR4_H2p0_B2p0_Cpre_Dcrt/gbr4_bgc_simple_2016-06.nc")
#' }

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
                          verbosity = 1,
                          return_list = FALSE,
                          date_format = "date")
{
  if (is.null(dim(location_latlon))) {
     location_latlon <- array(location_latlon, c(1,2))
  }
  input_file <- substitute_filename(input_file)
  if (layer=='integrated') return(get_ereefs_depth_integrated_ts(var_names, location_latlon, start_date, end_date, input_file, input_grid, eta_stem, override_positive))
  if ((layer=='bottom')&&(length(location_latlon[,1])>1)) stop('Only one location can be given if layer==bottom')
#  if (layer=='bottom') return(get_ereefs_bottom_ts(var_names, location_latlon, start_date, end_date, input_file, input_grid, eta_stem, override_positive))
  if (layer < 0) { 
     z_grid <- get_ereefs_grids(input_file, input_grid)[['z_grid']]
     layer <- which.min(z_grid<layer)
  }
  if (is.character(location_latlon)&&(location_latlon=="mmp")) {
     location_latlon <- data.frame(latitude=mmp_sites$latitude, longitude=mmp_sites$longitude)
     mmp <- TRUE
  } else {
     mmp <- FALSE
  }

  # Check whether netcdf output files are daily (case 1), monthly (case 4) or something else (case 0)
  ereefs_case <- get_ereefs_case(input_file)
  input_stem <- get_file_stem(input_file)

  # Dates to plot
  if (is.vector(start_date)) {
    if ((length(start_date==2)) && is.character(start_date[1])) { 
          start_date <- chron::chron(start_date[1], start_date[2], format=c('d-m-y', 'h:m:s'), 
                                     origin=c(year=1990, month=1, day=1))
    } else if (length(start_date==3)) { 
      # Set time to midday
      start_date <- chron::chron(paste(start_date[3], start_date[2], start_date[1], sep = '-'), "12:00:00", format=c('d-m-y', 'h:m:s'), 
                                 origin=c(year=1990, month=1, day=1))
    } else if (length(start_date==4)) {
       if (!is.character(start_date[4])) start_date[4] <- paste0(start_date[4], ':00')
       start_date <- chron::chron(paste(start_date[3], start_date[2], start_date[1], sep = '-'), start_date[4], format=c('d-m-y', 'h:m:s'), 
                                  origin=c(year=1990, month=1, day=1)) 
    } else {
      stop("start_date format not recognised")
    }
  } else if (is.character(start_date)) {
    start_date <- chron::chron(start_date, "12:00:00", format=c('d-m-y', 'h:m:s'), 
                                 origin=c(year=1990, month=1, day=1))
  } else if (class(start_date)[1] == "Date") {
    start_date <- chron::chron(as.character(start_date, format="%d-%m-%y"), "12:00:00", format=c('d-m-y', 'h:m:s'), 
                                 origin=c(year=1990, month=1, day=1))
    #start_date <- chron::as.chron(start_date) # I don't know why this does not work.
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
    end_date <- chron::chron(as.character(end_date, format="%d-%m-%y"), "12:00:00", format=c('d-m-y', 'h:m:s'), 
                                 origin=c(year=1990, month=1, day=1))
    #end_date <- chron::as.chron(end_date)
  }
  end_day <- as.integer(chron::days(end_date))
  end_month <- as.integer(months(end_date))
  end_year <- as.integer(as.character(chron::years(end_date)))
  
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

  # Should make this a function as this block of code is also used by other functions
  if (ereefs_case[2] == '4km') {
      input_file <- paste0(input_stem, format(as.Date(paste(start_year, start_month, 1, sep='-')), '%Y-%m'), 
			  '.nc')
    dum1 <- get_origin_and_times(input_file)
    ereefs_origin <- dum1[[1]]
    ds <- dum1[[2]]
    blank_length <- as.numeric(end_date - start_date + 1) / as.numeric(ds[2] - ds[1])
  } else if (ereefs_case[2] == '1km') {
      input_file <- paste0(input_stem, format(as.Date(paste(start_year, start_month, start_day, sep='-')), '%Y-%m-%d'), 
			  '.nc')
      if (verbosity>1) print(input_file)
      if (!is.na(eta_stem)) etafile <- paste0(eta_stem, format(as.Date(paste(start_year, start_month, start_day, sep='-')), '%Y-%m-%d'), 
			  '.nc')
      blank_length <- end_date - start_date + 1
      ereefs_origin <- get_origin_and_times(input_file)[[1]]
      # We don't need ds in this case because we are taking just one output per day
  } else if (ereefs_case[1] == "thredds_catalog") {
    # (We are currently in get_ereefs_ts())
      if (!is.na(eta_stem)) etafile <- paste0(eta_stem, '.nc')
      catalog_list <- thredds::tds_list_datasets(input_file)
      catalog_list <- catalog_list[stringr::str_ends(catalog_list$path, "nc"), ]$path
      catalog_list <- stringr::str_replace(catalog_list, "catalog/.*dataset=fx3-", stringr::fixed("dodsC/fx3/"))
      # Special case to remove unwanted additional files from certain catalogues
      in_range <- stringr::str_which(catalog_list, "recom_wc", negate=TRUE)
      catalog_list <- catalog_list[in_range]
      catalog_startdates <- chron::chron(rep(1, length(catalog_list)), origin=c(year=1990,month=1,day=1), format="y-m-d")
      catalog_enddates <- catalog_startdates

      # Pare down the catalog_list to include only those that cover the date range of interest, and set up a list of
      # the times of outputs in each of these relevant files. This is slow.
      # We could save this time by saving the results to sysdata.rda so we can look them up for known catalogs and only 
      # doing this if we encounter an unknown catalog. Alternatively, we could save some time by assuming that the list is
      # in temporal order and only looking for the first and last file needed.
      if (verbosity>0) print("Finding out which files in the catalog are relevant to the time range...")
      catalog_times <- vector("list", length(catalog_list))
      in_range <- rep(FALSE, length(catalog_list))
      blank_length <- 0
      for (i in 1:length(catalog_list)) {
         if (verbosity >1) print(paste("File", i, catalog_list[i]))
         catalog_times[[i]] <- get_origin_and_times(catalog_list[i])[[2]]
         catalog_startdates[i] <- catalog_times[[i]][1]
         catalog_enddates[i] <- catalog_times[[i]][length(catalog_times[[i]])]
         dum1 <- length(which((catalog_times[[i]] >= start_date)&(catalog_times[[i]] <= end_date)))
         if (dum1 >0) {
           blank_length <- blank_length + dum1
           in_range[i] <- TRUE
         }
      }

      # Let's make sure these are in temporal order. Note that sort.chron() ignores index.return so we need to convert to numeric.
      ix <- sort(as.numeric(catalog_startdates[in_range]), index.return=TRUE)$ix
      ix <- which(in_range)[ix]
      catalog_list <- catalog_list[ix]
      catalog_times <- catalog_times[ix]
      catalog_startdates <- catalog_startdates[ix]
      catalog_enddates <- catalog_enddates[ix]
      icatalog <- 1
      input_file <- catalog_list[icatalog]
      ds <- catalog_times[[icatalog]]
      ereefs_origin <- get_origin_and_times(catalog_list[icatalog])[[1]]
      blank_length <- as.numeric(end_date - start_date + 1) / as.numeric(ds[2] - ds[1]) #+ 0.5/(as.numeric(ds[2] - ds[1]))
  } else {
      # We are looking at a single netcdf file such as a RECOM output file
      if (!is.na(eta_stem)) etafile <- paste0(eta_stem, '.nc')
      input_file <- input_file
      dum1 <- get_origin_and_times(input_file)
      ereefs_origin <- dum1[[1]]
      ds <- dum1[[2]]
      blank_length <- as.numeric(end_date - start_date + 1) / as.numeric(ds[2] - ds[1]) #+ 0.5/(as.numeric(ds[2] - ds[1]))
  }

  if (is.integer(location_latlon)) {
     # We have specified grid coordinates rather than geocoordinates
     location_grid <- location_latlon
  } else { 
    # We have geocoordinates. Find the nearest grid-points to the sampling location

    # First, get the model grid
    nc <- safe_nc_open(input_file)
    if (is.null(nc$var[['latitude']])) {
      # Not a simple format netcdf file, so assume it's a full EMS netcdf file.
      latitude <- safe_ncvar_get(nc, "y_centre")
      longitude <- safe_ncvar_get(nc, "x_centre")
    } else { 
      # Simple format netcdf file
      latitude <- safe_ncvar_get(nc, "latitude")
      longitude <- safe_ncvar_get(nc, "longitude")
    }
    ncdf4::nc_close(nc)

    if (is.null(dim(location_latlon))) {
       # Just one location
       grid_index <- (latitude - location_latlon[1])^2 + (longitude - location_latlon[2])^2 
       grid_index <- which.min(grid_index) 
    } else { 
       # Multiple locations
       if (class(location_latlon)[1] != "data.frame") {
          # location_latlon has been provided as an array/matrix. Coerce it into a data frame for consistency.
          location_latlon <- data.frame(latitude = location_latlon[,1], longitude = location_latlon[,2])
       }
       grid_index <- apply(location_latlon,1, function(ll) which.min((latitude - ll[1])^2 + (longitude - ll[2])^2)) 
    }
    #location_grid <- arrayInd(grid_index, dim(latitude))
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
  #ts_frame <- data.frame(blanks, array(blanks, dim=c(length(blanks), length(var_names))))
  colnames(ts_frame) <- c("date", var_names)

  # Loop through relevant daily or monthly eReefs files to extract the data
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
    } else if ((ereefs_case[2] == 'recom')|(ereefs_case[1] == "ncml")) { 
      day_count <- day_count / as.numeric(ds[2]-ds[1])
      if (day_count > length(ds)) {
        warning(paste('end_date', end_date, 'is beyond available data. Ending at', ds[length(ds)]))
        day_count <- length(ds)
      }
      from_day <- (as.numeric(chron::chron(paste(year, month, from_day, sep = '-'), format=c('y-m-d'),
                                  origin=c(year=1990, month=1, day=1)) - ds[1]) + 
                              start_tod) / as.numeric(ds[2] - ds[1]) + 1 
	    if (from_day<1) from_day <-1
	    fileslist <- 1
    } else if (ereefs_case[1] == "thredds_catalog") { 
      month_startdate <- chron::chron(paste(year, month, 1, sep = '-'), format = 'y-m-d',
                                      origin=c(year=1990, month=1, day=1))
      month_enddate <- chron::chron(paste(year, month, daysIn(as.Date(paste(year, month, 1, sep='-'))) , sep = '-'), format = 'y-m-d',
                                      origin=c(year=1990, month=1, day=1))
      day_count <- day_count / as.numeric(ds[2]-ds[1])
      ix <- which((catalog_startdates >= month_startdate) & (catalog_enddates <= (month_enddate + 0.999)))
      icatalog <- ix[1]
      fileslist <- catalog_list[ix]
      if (end_date > catalog_enddates[length(catalog_enddates)]) {
        warning(paste('end_date', end_date, 'is beyond available data. Ending at', catalog_enddates[length(catalog_enddates)]))
        day_count <- day_count - (end_date - catalog_enddates[length(catalog_enddates)])
      }
      from_day <- (as.numeric(chron::chron(paste(year, month, from_day, sep = '-'), format='y-m-d',
                               origin=c(year=1990, month=1, day=1)) - catalog_startdates[ix[1]]) 
                   + start_tod) / as.numeric(ds[2] - ds[1]) + 1 
      if (from_day<1) from_day <- 1
    } else stop("Shouldn't happen: ereefs_case not recognised")

    for (dcount in 1:length(fileslist)) {
      if (ereefs_case[2] == '1km') {
	      input_file <- paste0(input_stem, format(as.Date(paste(year, month, fileslist[dcount], sep="-")), '%Y-%m-%d'), '.nc')
      } else if (ereefs_case[1] == "thredds_catalog") {
        ## Are we past the end date of the current netcdf file from the catalog? If so, move on to the next one.
        #if (chron::chron(paste(year, month, dcount, sep = '-'), format = c('y-m-d'), origin = c(year=1990, month=1, day=1)) > 
        #    catalog_times[[icatalog]][length(catalog_times[[icatalog]])]) {
          icatalog <- icatalog[dcount]
          input_file <- catalog_list[icatalog]
          ds <- catalog_times[[icatalog]]
        #}
      }
      #input_file <- paste0(input_file, '?', var_list, ',time')
      nc <- safe_nc_open(input_file)
      if ((ereefs_case[2] == "1km")||(ereefs_case[2] == "4km")) {
          if (!is.null(nc$var[['t']])) {
            d <- as.Date(safe_ncvar_get(nc, "t"), origin = as.Date("1990-01-01"))[from_day:(from_day+day_count-1)]
          } else {
            d <- as.Date(safe_ncvar_get(nc, "time"), origin = as.Date("1990-01-01"))[from_day:(from_day+day_count-1)]
          }
      } else if ((ereefs_case[2] == "recom")|(ereefs_case[1] == "ncml")) { 
         d <- ds[from_day:(from_day + day_count - 1)]
      } else if (ereefs_case[1] == "thredds_catalog") {
         d <- ds[from_day:(from_day + day_count - 1)]
      } else stop("Shouldn't happen: ereefs_case not recognised")
      im1 = i+1
      i <- i + length(d)

      ts_frame[im1:i,"date",] <- d
      for (j in 1:length(var_names)) {
          if (is.na(ndims[j])) {
             # We don't yet know the dimensions of the variable, so let's get them
              dims <- nc$var[[var_names[j]]][['size']]
              if (is.null(dims)) stop(paste(var_names[j], ' not found in netcdf file.')) 
              ndims[j] <- length(dims)
           }
           if (layer == 'surface') { 
              dims <- nc$var[[var_names[j]]][['size']]
              layer_actual <- dims[3] 
           } else if (layer == 'bottom') { 
              # find the bottom layer -- assume it is the first layer that is not NA
              if (ndims[j] == 4) { 
                 dum1 <- safe_ncvar_get(nc, var_names[j], start=c(startv,1,from_day), count=c(countv,-1,1)) 
                 if (length(which(!is.na(dum1)))==0) stop('Location given is a dry cell at start time. It is probably on land.')
                 layer <- min(which(!is.na(dum1)))
                 layer_actual <- layer
                 if (verbosity>0) print(paste('bottom layer = ', layer))
              } 
           } else { 
              layer_actual <- layer 
           }

          if (ndims[j] == 4) {
             wc <- safe_ncvar_get(nc, var_names[j], start=c(startv,layer_actual,from_day), count=c(countv,1,day_count))
             #ts_frame[im1:i, j+1] <- safe_ncvar_get(nc, var_names[j], start=c(location_grid[2],location_grid[1],layer_actual[j],from_day), count=c(1,1,1,day_count))
          } else {
             #print(c(startv,from_day))
             #print(c(countv,day_count))
             wc <- safe_ncvar_get(nc, var_names[j], start=c(startv,from_day), count=c(countv,day_count))
             #ts_frame[im1:i, j+1] <- safe_ncvar_get(nc, var_names[j], start=c(location_grid[2],location_grid[1],from_day), count=c(1,1,day_count))
          }
          #print(wc)
          wc <- array(wc, c(countv[1]*countv[2], day_count))
          ts_frame[im1:i, j+1,] <- t(wc[grid_index,])
      }
      ncdf4::nc_close(nc)
      if (verbosity>0) setTxtProgressBar(pb,mcount)
    }
  }
  if (verbosity>0) close(pb)

  ts_frame = lapply(seq(dim(ts_frame)[3]),
                   function(x) data.frame(date = chron::chron(as.vector(ts_frame[, 1, x]),
                                    format=c('y-m-d', 'h:m:s')),
                                 ts_frame[, 2:dim(ts_frame)[2], x]))
  if (date_format == "date") ts_frame$date <- lapply(seq(dim(ts_frame)[3]),
                   function(x) data.frame(date = as.Date(ts_frame$date) + (as.numeric(ts_frame$date) - floor(as.numeric(ts_frame$date))),
                                 ts_frame[, 2:dim(ts_frame)[2], x]))
  if (numpoints == 1) ts_frame <- data.frame(ts_frame[[1]])

  #ts_frame <- lapply(seq(dim(ts_frame)[3]), function(x) data.frame(date=as.Date(as.vector(ts_frame[ ,1, 1]), origin="1970-01-01"), ts_frame[ ,2:dim(ts_frame)[2] , x])) 
    #ts_frame$date <- (chron::chron(chron::as.chron(ts_frame$date, origin=c(year=1990, month=1, day=1)), origin=c(1,1,1990), format=c('y-m-d', 'h:m:s')))
    #ts_frame$date <- ts_frame$date + as.numeric(ereefs_origin)
  if (mmp) names(ts_frame) <- mmp_sites$Name
  if (length(var_names)==1) names(ts_frame)[2] <- var_names
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
#' \dontrun{
#' get_ereefs_bottom_ts(c('Chl_a_sum', 'NH4'), location_latlon=data.frame(latitide=-23.39189, longitude=150.88852), layer='surface', start_date=c(2010,12,31),end_date=c(2011,1,5), input_file=2)
#'}
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
  #check_platform_ok(input_stem)
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

  if (ereefs_case[2] == '4km') {
      input_file <- paste0(input_stem, format(as.Date(paste(start_year, start_month, 1, sep='-')), '%Y-%m'), 
			  '.nc')
      if (!is.na(eta_stem)) etafile  <- paste0(eta_stem, format(as.Date(paste(start_year, start_month, 1, sep='-')), '%Y-%m'), 
			  '.nc')
      ds <- get_origin_and_times(input_file)[[2]]
      blank_length <- as.numeric(end_date - start_date + 1) / as.numeric(ds[2] - ds[1])
  } else if (ereefs_case[2] == '1km') {
      input_file <- paste0(input_stem, format(as.Date(paste(start_year, start_month, start_day, sep='-')), '%Y-%m-%d'), 
			  '.nc')
      if (!is.na(eta_stem)) etafile <- paste0(eta_stem, format(as.Date(paste(start_year, start_month, start_day, sep='-')), '%Y-%m-%d'), 
			  '.nc') 
      blank_length <- end_date - start_date + 1
  } else if (ereefs_case[1] == "thredds_catalog") {
    # (We are currently in get_ereefs_bottom_ts())
      if (!is.na(eta_stem)) etafile <- paste0(eta_stem, '.nc')
      catalog_list <- thredds::tds_list_datasets(input_file)
      catalog_list <- catalog_list[stringr::str_ends(catalog_list$path, "nc"), ]$path
      catalog_list <- stringr::str_replace(catalog_list, "catalog/.*dataset=fx3-", stringr::fixed("dodsC/fx3/"))
      # Special case to remove unwanted additional files from certain catalogues
      in_range <- stringr::str_which(catalog_list, "recom_wc", negate=TRUE)
      catalog_list <- catalog_list[in_range]
      catalog_startdates <- chron::chron(rep(1, length(catalog_list)), origin=c(year=1990,month=1,day=1), format="y-m-d")
      catalog_enddates <- catalog_startdates

      # Pare down the catalog_list to include only those that cover the date range of interest, and set up a list of
      # the times of outputs in each of these relevant files. This is slow.
      # We could save this time by saving the results to sysdata.rda so we can look them up for known catalogs and only 
      # doing this if we encounter an unknown catalog. Alternatively, we could save some time by assuming that the list is
      # in temporal order and only looking for the first and last file needed.
      if (verbosity>0) print("Finding out which files in the catalog are relevant to the time range...")
      catalog_times <- vector("list", length(catalog_list))
      in_range <- rep(FALSE, length(catalog_list))
      blank_length <- 0
      for (i in 1:length(catalog_list)) {
         if (verbosity >1) print(paste("File", i, catalog_list[i]))
         catalog_times[[i]] <- get_origin_and_times(catalog_list[i])[[2]]
         catalog_startdates[i] <- catalog_times[[i]][1]
         catalog_enddates[i] <- catalog_times[[i]][length(catalog_times[[i]])]
         dum1 <- length(which((catalog_times[[i]] >= start_date)&(catalog_times[[i]] <= end_date)))
         if (dum1 >0) {
           blank_length <- blank_length + dum1
           in_range[i] <- TRUE
         }
      }

      # Let's make sure these are in temporal order. Note that sort.chron() ignores index.return so we need to convert to numeric.
      ix <- sort(as.numeric(catalog_startdates[in_range]), index.return=TRUE)$ix
      ix <- which(in_range)[ix]
      catalog_list <- catalog_list[ix]
      catalog_times <- catalog_times[ix]
      catalog_startdates <- catalog_startdates[ix]
      catalog_enddates <- catalog_enddates[ix]
      icatalog <- 1
      input_file <- catalog_list[icatalog]
      ds <- catalog_times[[icatalog]]
      blank_length <- as.numeric(end_date - start_date + 1) / as.numeric(ds[2] - ds[1]) #+ 0.5/(as.numeric(ds[2] - ds[1]))
  } else {
      # We are looking at a single netcdf file such as a RECOM output file
      if (!is.na(eta_stem)) etafile <- paste0(eta_stem, '.nc')
      input_file <- input_file
      dum1 <- get_origin_and_times(input_file)
      ereefs_origin <- dum1[[1]]
      ds <- dum1[[2]]
      blank_length <- as.numeric(end_date - start_date + 1) / as.numeric(ds[2] - ds[1]) #+ 0.5/(as.numeric(ds[2] - ds[1]))
  }

  # Initialise
  blanks <- rep(NA, blank_length)
  ts_frame <- data.frame(as.Date(blanks), array(blanks, dim=c(length(blanks), length(var_names))))
  names(ts_frame) <- c("date", var_names)

  nc <- safe_nc_open(input_file)
  if (!is.na(eta_stem)) nc3 <- safe_nc_open(etafile)

  if (is.integer(location_latlon)) {
     location_grid <- location_latlon
  } else { 
    # Find the nearest grid-points to the sampling location
    if (is.null(nc$var[['latitude']])) {
      latitude <- safe_ncvar_get(nc, "y_centre")
      longitude <- safe_ncvar_get(nc, "x_centre")
    } else { 
      latitude <- safe_ncvar_get(nc, "latitude")
      longitude <- safe_ncvar_get(nc, "longitude")
    }
    tmp <- (latitude - location_latlon[1])^2 + (longitude - location_latlon[2])^2 
    tmp <- which.min(tmp) 
    #location_grid <- arrayInd(tmp, dim(latitude))
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
  botz <- zsign * as.numeric(safe_ncvar_get(nc, "botz", start=c(location_grid[2], location_grid[1]), count=c(1,1)))
  ncdf4::nc_close(nc)

  # Loop through relevant eReefs files to extract the data
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
     if (ereefs_case[2] == '4km') { 
        fileslist <- 1
        input_file <- paste0(input_stem, format(as.Date(paste(year, month, 1, sep="-")), '%Y-%m'), '.nc')
	     day_count <- day_count / as.numeric(ds[2]-ds[1])
        if (!is.na(eta_stem)) etafile <- paste0(eta_stem, format(as.Date(paste(year, month, 1, sep="-")), '%Y-%m'), '.nc')
     } else if (ereefs_case[2] == '1km') {
        fileslist <- from_day:(from_day+day_count-1)
        from_day <- 1
        day_count <- 1
    } else if ((ereefs_case[2] == 'recom')|(ereefs_case[1] == "ncml")) { 
      day_count <- day_count / as.numeric(ds[2]-ds[1])
      if (day_count > length(ds)) {
        warning(paste('end_date', end_date, 'is beyond available data. Ending at', ds[length(ds)]))
        day_count <- length(ds)
      }
      from_day <- (as.numeric(chron::chron(paste(year, month, from_day, sep = '-'), format=c('y-m-d'),
                                  origin=c(year=1990, month=1, day=1)) - ds[1]) + 
                              start_tod) / as.numeric(ds[2] - ds[1]) + 1 
	    if (from_day<1) from_day <-1
	    fileslist <- 1
    } else if (ereefs_case[1] == "thredds_catalog") { 
      month_startdate <- chron::chron(paste(year, month, 1, sep = '-'), format = 'y-m-d',
                                      origin=c(year=1990, month=1, day=1))
      month_enddate <- chron::chron(paste(year, month, daysIn(as.Date(paste(year, month, 1, sep='-'))) , sep = '-'), format = 'y-m-d',
                                      origin=c(year=1990, month=1, day=1))
      day_count <- day_count / as.numeric(ds[2]-ds[1])
      ix <- which((catalog_startdates >= month_startdate) & (catalog_enddates <= (month_enddate + 0.999)))
      icatalog <- ix[1]
      fileslist <- catalog_list[ix]
      if (end_date > catalog_enddates[length(catalog_enddates)]) {
        warning(paste('end_date', end_date, 'is beyond available data. Ending at', catalog_enddates[length(catalog_enddates)]))
        day_count <- day_count - (end_date - catalog_enddates[length(catalog_enddates)])
      }
      from_day <- (as.numeric(chron::chron(paste(year, month, from_day, sep = '-'), format='y-m-d',
                               origin=c(year=1990, month=1, day=1)) - catalog_startdates[ix[1]]) 
                   + start_tod) / as.numeric(ds[2] - ds[1]) + 1 
      if (from_day<1) from_day <- 1
    } else stop("Shouldn't happen: ereefs_case not recognised")

     for (dcount in 1:length(fileslist)) {
        if (ereefs_case[2] == '1km') {
	        input_file <- paste0(input_stem, format(as.Date(paste(year, month, fileslist[dcount], sep="-")), '%Y-%m-%d'), '.nc')
          if (!is.na(eta_stem)) etafile <- paste0(eta_stem, format(as.Date(paste(year, month, fileslist[dcount], sep="-")), '%Y-%m-%d'), '.nc')
        }
        #input_file <- paste0(input_file, '?', var_list, ',time,eta')
        nc <- safe_nc_open(input_file)
        if (!is.na(eta_stem)) nc3 <- safe_nc_open(etafile)
        # Get dates
        if ((ereefs_case[2] == "1km")||(ereefs_case[2] == "4km")) {
           if (!is.null(nc$var[['t']])) { 
              ds <- as.Date(safe_ncvar_get(nc, "t"), origin = as.Date("1990-01-01")) 
           } else { 
              ds <- as.Date(safe_ncvar_get(nc, "time"), origin = as.Date("1990-01-01")) 
           }
	        d <- ds[from_day:(from_day + day_count - 1)]
        } else if ((ereefs_case[2] == "recom")|(ereefs_case[1] == "ncml")) { 
           d <- ds[from_day:(from_day + day_count - 1)]
        } else if (ereefs_case[1] == "thredds_catalog") {
           d <- ds[from_day:(from_day + day_count - 1)]
        } else stop("Shouldn't happen: ereefs_case not recognised")

        if (!is.null(nc$var[['eta']])) { 
           eta <- safe_ncvar_get(nc, 'eta', start=c(location_grid[2], location_grid[1],from_day), count=c(1,1,day_count)) 
        } else { 
           eta <- safe_ncvar_get(nc3, 'eta', start=c(location_grid[2], location_grid[1],from_day), count=c(1,1,day_count)) 
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
          wc <- safe_ncvar_get(nc, var_names[j], start=c(location_grid[2],location_grid[1],1,from_day), count=c(1,1,-1,day_count))
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
#' @param verbosity (Defailt 1) how much do you want to know about progress?
#' @param mass (Default FALSE) Set to true if you want the mass per square metre rather than the mean concentration over depth returned
#' @param date_format (Default "date"). Set to "chron" if you'd like the date returned in chron format.
#' @export
#' @examples
#' \dontrun{
#' get_ereefs_depth_integrated_ts(c('Chl_a_sum', 'NH4'), location_latlon=data.frame(latitide=-23.39189, longitude=150.88852), layer='surface', start_date=c(2010,12,31),end_date=c(2011,1,5), input_file=2)
#'}
get_ereefs_depth_integrated_ts <- function(var_names=c('Chl_a_sum', 'TN'), 
                        location_latlon=c(-23.39189, 150.88852), 
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
  input_file <- substitute_filename(input_file)
  # Check whether this is a GBR1 or GBR4 ereefs file, or something else
  ereefs_case <- get_ereefs_case(input_file)
  input_stem <- get_file_stem(input_file)
  #check_platform_ok(input_stem)
  z_grid <- get_ereefs_grids(input_file, input_grid)[['z_grid']]

  # Dates to plot
  if (is.vector(start_date)) {
    if ((length(start_date==2)) && is.character(start_date[1])) { 
          start_date <- chron::chron(start_date[1], start_date[2], format=c('d-m-y', 'h:m:s'), 
                                     origin=c(year=1990, month=1, day=1))
    } else if (length(start_date==3)) { 
      # Set time to midday
      start_date <- chron::chron(paste(start_date[3], start_date[2], start_date[1], sep = '-'), "12:00:00", format=c('d-m-y', 'h:m:s'), 
                                 origin=c(year=1990, month=1, day=1))
    } else if (length(start_date==4)) {
       if (!is.character(start_date[4])) start_date[4] <- paste0(start_date[4], ':00')
       start_date <- chron::chron(paste(start_date[3], start_date[2], start_date[1], sep = '-'), start_date[4], format=c('d-m-y', 'h:m:s'), 
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
 # Note that origin format must be m/d/y
  if (ereefs_case[2] == '4km') {
      input_file <- paste0(input_stem, format(as.Date(paste(start_year, start_month, 1, sep='-')), '%Y-%m'), 
			  '.nc')
      if (verbosity>1) print(input_file)
      if (!is.na(eta_stem)) etafile  <- paste0(eta_stem, format(as.Date(paste(start_year, start_month, 1, sep='-')), '%Y-%m'), 
			  '.nc') 
      dum1 <- get_origin_and_times(input_file)
      ereefs_origin <- dum1[[1]]
      ds <- dum1[[2]]
      blank_length <- as.numeric(end_date - start_date + 1) / as.numeric(ds[2] - ds[1])
  } else if (ereefs_case[2] == '1km') {
      input_file <- paste0(input_stem, format(as.Date(paste(start_year, start_month, start_day, sep='-')), '%Y-%m-%d'), 
			  '.nc')
      if (verbosity>1) print(input_file)
      if (!is.na(eta_stem)) etafile <- paste0(eta_stem, format(as.Date(paste(start_year, start_month, start_day, sep='-')), '%Y-%m-%d'), 
			  '.nc')
      blank_length <- end_date - start_date + 1
      dum1 <- get_origin_and_times(input_file)
      ereefs_origin <- dum1[[1]]
      ds <- dum1[[2]]
  } else if (ereefs_case[1] == "thredds_catalog") {
    # (We are currently in get_ereefs_depth_integrated_ts())
      if (!is.na(eta_stem)) etafile <- paste0(eta_stem, '.nc')
      catalog_list <- thredds::tds_list_datasets(input_file)
      catalog_list <- catalog_list[stringr::str_ends(catalog_list$path, "nc"), ]$path
      catalog_list <- stringr::str_replace(catalog_list, "catalog/.*dataset=fx3-", stringr::fixed("dodsC/fx3/"))
      # Special case to remove unwanted additional files from certain catalogues
      in_range <- stringr::str_which(catalog_list, "recom_wc", negate=TRUE)
      catalog_list <- catalog_list[in_range]
      catalog_startdates <- chron::chron(rep(1, length(catalog_list)), origin=c(year=1990,month=1,day=1), format="y-m-d")
      catalog_enddates <- catalog_startdates

      # Pare down the catalog_list to include only those that cover the date range of interest, and set up a list of
      # the times of outputs in each of these relevant files. This is slow.
      # We could save this time by saving the results to sysdata.rda so we can look them up for known catalogs and only 
      # doing this if we encounter an unknown catalog. Alternatively, we could save some time by assuming that the list is
      # in temporal order and only looking for the first and last file needed.
      if (verbosity>0) print("Finding out which files in the catalog are relevant to the time range...")
      catalog_times <- vector("list", length(catalog_list))
      in_range <- rep(FALSE, length(catalog_list))
      blank_length <- 0
      for (i in 1:length(catalog_list)) {
         if (verbosity >1) print(paste("File", i, catalog_list[i]))
         catalog_times[[i]] <- get_origin_and_times(catalog_list[i])[[2]]
         catalog_startdates[i] <- catalog_times[[i]][1]
         catalog_enddates[i] <- catalog_times[[i]][length(catalog_times[[i]])]
         dum1 <- length(which((catalog_times[[i]] >= start_date)&(catalog_times[[i]] <= end_date)))
         if (dum1 >0) {
           blank_length <- blank_length + dum1
           in_range[i] <- TRUE
         }
      }

      # Let's make sure these are in temporal order. Note that sort.chron() ignores index.return so we need to convert to numeric.
      ix <- sort(as.numeric(catalog_startdates[in_range]), index.return=TRUE)$ix
      ix <- which(in_range)[ix]
      catalog_list <- catalog_list[ix]
      catalog_times <- catalog_times[ix]
      catalog_startdates <- catalog_startdates[ix]
      catalog_enddates <- catalog_enddates[ix]
      icatalog <- 1
      input_file <- catalog_list[icatalog]
      ds <- catalog_times[[icatalog]]
      blank_length <- as.numeric(end_date - start_date + 1) / as.numeric(ds[2] - ds[1]) #+ 0.5/(as.numeric(ds[2] - ds[1]))
  } else {
      # We are looking at a single netcdf file such as a RECOM output file
      if (!is.na(eta_stem)) etafile <- paste0(eta_stem, '.nc')
      input_file <- input_file
      dum1 <- get_origin_and_times(input_file)
      ereefs_origin <- dum1[[1]]
      ds <- dum1[[2]]
      blank_length <- as.numeric(end_date - start_date + 1) / as.numeric(ds[2] - ds[1]) #+ 0.5/(as.numeric(ds[2] - ds[1]))
  }

  if (!is.null(dim(location_latlon))) { 
    if (dim(location_latlon)[1] > 1) stop('Currently, get_ereefs_depth_integrated_ts() only supports a single location. This is on my to-do list to fix in future. Let me know if you would like this feature. b.robson@aims.gov.au')
  }

  nc <- safe_nc_open(input_file)
  if (!is.na(eta_stem)) nc3 <- safe_nc_open(etafile)

  if (is.null(dim(location_latlon))) {
     location_latlon <- array(location_latlon, c(1,2))
  }
  if (is.integer(location_latlon)) {
     # We have specified grid coordinates rather than geocoordinates
     location_grid <- location_latlon
  } else { 
    # We have geocoordinates. Find the nearest grid-points to the sampling location

    # First, get the model grid
    #nc <- safe_nc_open(input_file)
    if (is.null(nc$var[['latitude']])) {
      # Not a simple format netcdf file, so assume it's a full EMS netcdf file.
      latitude <- safe_ncvar_get(nc, "y_centre")
      longitude <- safe_ncvar_get(nc, "x_centre")
    } else { 
      # Simple format netcdf file
      latitude <- safe_ncvar_get(nc, "latitude")
      longitude <- safe_ncvar_get(nc, "longitude")
    }

    if (is.null(dim(location_latlon))) {
       # Just one location
       grid_index <- (latitude - location_latlon[1])^2 + (longitude - location_latlon[2])^2 
       grid_index <- which.min(grid_index) 
    } else { 
       # Multiple locations
       if (class(location_latlon)[1] != "data.frame") {
          # location_latlon has been provided as an array/matrix. Coerce it into a data frame for consistency.
          location_latlon <- data.frame(latitude = location_latlon[,1], longitude = location_latlon[,2])
       }
       grid_index <- apply(location_latlon,1, function(ll) which.min((latitude - ll[1])^2 + (longitude - ll[2])^2)) 
    }
    #location_grid <- arrayInd(grid_index, dim(latitude))
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
  if ((countv[2] == 1)&&(countv[1] != 1)) { 
   location_grid <- location_grid[,1]
  } else if ((countv[1] == 1)&&(countv[2] != 1)) {
   location_grid <- location_grid[,2]
  }

  # Initialise
  blanks <- rep(NA, blank_length)
  ts_frame <- data.frame(blanks, array(blanks, dim=c(length(blanks), length(var_names))))
#  if (!is.list(ts_frame))  
    names(ts_frame) <- c("date", var_names)

  zat <- ncdf4::ncatt_get(nc, "botz")
  if (!is.null(zat$positive)) {
    if (zat$positive=="down") zsign <- -1 else zsign <- 1
    if (override_positive) zsign <- 1
  } else {
   zsign <-1
    if (override_positive) zsign <- -1
  }
  botz <- zsign * as.numeric(safe_ncvar_get(nc, "botz", start=c(startv), count=c(countv)))
  #botz <- zsign * as.numeric(safe_ncvar_get(nc, "botz", start=c(location_grid[2], location_grid[1]), count=c(1,1)))
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
        if (verbosity) print(input_file)
        day_count <- day_count / as.numeric(ds[2]-ds[1])
        if (day_count > length(ds)) {
          warning(paste('end_date', end_date, 'is beyond available data. Ending at', ds[length(ds)]))
          day_count <- length(ds)
        }
        if (!is.na(eta_stem)) etafile <- paste0(eta_stem, format(as.Date(paste(year, month, 1, sep="-")), '%Y-%m'), '.nc')
     } else if (ereefs_case[2] == '1km') {
        fileslist <- from_day:(from_day+day_count-1)
        from_day <- 1
        day_count <- 1
     } else if ((ereefs_case[2] == 'recom')|(ereefs_case[1] == "ncml")) { 
       day_count <- day_count / as.numeric(ds[2]-ds[1])
       if (day_count > length(ds)) {
         warning(paste('end_date', end_date, 'is beyond available data. Ending at', ds[length(ds)]))
         day_count <- length(ds)
       }
       from_day <- (as.numeric(chron::chron(paste(year, month, from_day, sep = '-'), format=c('y-m-d'),
                                   origin=c(year=1990, month=1, day=1)) - ds[1]) + 
                               start_tod) / as.numeric(ds[2] - ds[1]) + 1 
	     if (from_day<1) from_day <-1
	     fileslist <- 1
     } else if (ereefs_case[1] == "thredds_catalog") { 
       month_startdate <- chron::chron(paste(year, month, 1, sep = '-'), format = 'y-m-d',
                                       origin=c(year=1990, month=1, day=1))
       month_enddate <- chron::chron(paste(year, month, daysIn(as.Date(paste(year, month, 1, sep='-'))) , sep = '-'), format = 'y-m-d',
                                       origin=c(year=1990, month=1, day=1))
       day_count <- day_count / as.numeric(ds[2]-ds[1])
       ix <- which((catalog_startdates >= month_startdate) & (catalog_enddates <= (month_enddate + 0.999)))
       icatalog <- ix[1]
       fileslist <- catalog_list[ix]
       if (end_date > catalog_enddates[length(catalog_enddates)]) {
         warning(paste('end_date', end_date, 'is beyond available data. Ending at', catalog_enddates[length(catalog_enddates)]))
         day_count <- day_count - (end_date - catalog_enddates[length(catalog_enddates)])
       }
       from_day <- (as.numeric(chron::chron(paste(year, month, from_day, sep = '-'), format='y-m-d',
                                origin=c(year=1990, month=1, day=1)) - catalog_startdates[ix[1]]) 
                    + start_tod) / as.numeric(ds[2] - ds[1]) + 1 
       if (from_day<1) from_day <- 1
     } else stop("Shouldn't happen: ereefs_case not recognised")

     for (dcount in 1:length(fileslist)) {
        if (ereefs_case[2] == '1km') {
	        input_file <- paste0(input_stem, format(as.Date(paste(year, month, fileslist[dcount], sep="-")), '%Y-%m-%d'), '.nc')
          if (!is.na(eta_stem)) etafile <- paste0(eta_stem, format(as.Date(paste(year, month, fileslist[dcount], sep="-")), '%Y-%m-%d'), '.nc')
          if (verbosity>1) print(input_file)
        } else if (ereefs_case[1] == "thredds_catalog") {
          icatalog <- icatalog[dcount]
          input_file <- catalog_list[icatalog]
          ds <- catalog_times[[icatalog]]
        }
        #input_file <- paste0(input_file, '?', var_list, ',time,eta')
        nc <- safe_nc_open(input_file)
        if (!is.na(eta_stem)) {
          nc3 <- safe_nc_open(etafile)
        } else if (is.null(nc$var[['eta']])) { 
          stop("Simple format files do not include surface elevation (needed for depth integration or depth below surface). Please either use a standard format file or provide another filename as eta_stem that contains matching eta data (e.g. from a hydrodynamic run).")
        }
        if ((ereefs_case[2] == "1km")||(ereefs_case[2] == "4km")) {
          if (!is.null(nc$var[['t']])) {
            d <- (safe_ncvar_get(nc, "t") + ereefs_origin)[from_day:(from_day+day_count-1)]
          } else {
            d <- (safe_ncvar_get(nc, "time") + ereefs_origin)[from_day:(from_day+day_count-1)]
          }
        } else if ((ereefs_case[2] == "recom")|(ereefs_case[1] == "ncml")) { 
           d <- ds[from_day:(from_day + day_count - 1)]
        } else if (ereefs_case[1] == "thredds_catalog") {
           d <- ds[from_day:(from_day + day_count - 1)]
        } else stop("Shouldn't happen: ereefs_case not recognised")

        if (!is.null(nc$var[['eta']])) { 
          eta <- safe_ncvar_get(nc, 'eta', start=c(startv,from_day), count=c(countv,day_count)) 
          #eta <- safe_ncvar_get(nc, 'eta', start=c(location_grid[2], location_grid[1],from_day), count=c(1,1,day_count)) 
        } else {
          if (!is.na(eta_stem)) {
            nc3 <- safe_nc_open(etafile)
          } else if (is.null(nc$var[['eta']])) { 
            stop("Simple format files do not include surface elevation (needed for depth integration or depth below surface). Please either use a standard format file or provide another filename as eta_stem that contains matching eta data (e.g. from a hydrodynamic run).")
          }
          eta <- safe_ncvar_get(nc3, 'eta', start=c(startv,from_day), count=c(countv,day_count)) 
          #eta <- safe_ncvar_get(nc3, 'eta', start=c(location_grid[2], location_grid[1],from_day), count=c(1,1,day_count)) 
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
           wc <- safe_ncvar_get(nc, var_names[j], start=c(startv,1,from_day), count=c(countv,-1,day_count))
          #wc <- safe_ncvar_get(nc, var_names[j], start=c(location_grid[2],location_grid[1],1,from_day), count=c(1,1,-1,day_count))
	        if (dim(dz)[2] == 1) wc <- array(wc, dim=dim(dz))
          if (mass) {
            ts_frame[im1:i, j+1] <- colSums(dz * wc, na.rm=TRUE)
          } else { 
            # take the depth-integrated average over the water column
            ts_frame[im1:i, j+1] <- colSums(dz * wc, na.rm=TRUE) / colSums(dz) 
          }
        }
        ncdf4::nc_close(nc)
        if (!is.na(eta_stem)) ncdf4::nc_close(nc3)
        if (verbosity>0) setTxtProgressBar(pb,mcount)
    }
  }
  if (verbosity>0) close(pb)
  ts_frame$date <- (chron::chron(chron::as.chron(ts_frame$date, origin=c(year=1990, month=1, day=1)), origin=c(1,1,1990), format=c('y-m-d', 'h:m:s')))
  if (date_format == "date") ts_frame$date <- as.Date(ts_frame$date) + (as.numeric(ts_frame$date) - floor(as.numeric(ts_frame$date)))
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
#' \dontrun{
#' get_ereefs_depth_specified_ts(c('Chl_a_sum', 'NH4'), depth=5.0, location_latlon=data.frame(latitide=-23.39189, longitude=150.88852), layer='surface', start_date=c(2010,12,31),end_date=c(2011,1,5), input_file=2)
#'}
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
  #check_platform_ok(input_stem)
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

  if (ereefs_case[2] == '4km') {
      input_file <- paste0(input_stem, format(as.Date(paste(start_year, start_month, 1, sep='-')), '%Y-%m'), 
			  '.nc')
      if (!is.na(eta_stem)) etafile  <- paste0(eta_stem, format(as.Date(paste(start_year, start_month, 1, sep='-')), '%Y-%m'), 
			  '.nc')
	    nc <- safe_nc_open(input_file)
	    if (!is.null(nc$var[['t']])) { 
	        ds <- as.Date(safe_ncvar_get(nc, "t"), origin = as.Date("1990-01-01"))
            } else {
	        ds <- as.Date(safe_ncvar_get(nc, "time"), origin = as.Date("1990-01-01"))
	    }
	    ncdf4::nc_close(nc)
      blank_length <- as.numeric(end_date - start_date + 1) / as.numeric(ds[2] - ds[1])
	    # '.nc?latitude,longitude')
  } else if (ereefs_case[2] == '1km') {
      input_file <- paste0(input_stem, format(as.Date(paste(start_year, start_month, start_day, sep='-')), '%Y-%m-%d'), 
			  '.nc')
      if (!is.na(eta_stem)) etafile <- paste0(eta_stem, format(as.Date(paste(start_year, start_month, start_day, sep='-')), '%Y-%m-%d'), 
			  '.nc')
      blank_length <- end_date - start_date + 1
			  # '.nc?latitude,longitude')
  } else if (ereefs_case[1] == "thredds_catalog") {
    # (We are currently in get_ereefs_depth_specified_ts())
      if (!is.na(eta_stem)) etafile <- paste0(eta_stem, '.nc')
      catalog_list <- thredds::tds_list_datasets(input_file)
      catalog_list <- catalog_list[stringr::str_ends(catalog_list$path, "nc"), ]$path
      catalog_list <- stringr::str_replace(catalog_list, "catalog/.*dataset=fx3-", stringr::fixed("dodsC/fx3/"))
      # Special case to remove unwanted additional files from certain catalogues
      in_range <- stringr::str_which(catalog_list, "recom_wc", negate=TRUE)
      catalog_list <- catalog_list[in_range]
      catalog_startdates <- chron::chron(rep(1, length(catalog_list)), origin=c(year=1990,month=1,day=1), format="y-m-d")
      catalog_enddates <- catalog_startdates

      # Pare down the catalog_list to include only those that cover the date range of interest, and set up a list of
      # the times of outputs in each of these relevant files. This is slow.
      # We could save this time by saving the results to sysdata.rda so we can look them up for known catalogs and only 
      # doing this if we encounter an unknown catalog. Alternatively, we could save some time by assuming that the list is
      # in temporal order and only looking for the first and last file needed.
      if (verbosity>0) print("Finding out which files in the catalog are relevant to the time range...")
      catalog_times <- vector("list", length(catalog_list))
      in_range <- rep(FALSE, length(catalog_list))
      blank_length <- 0
      for (i in 1:length(catalog_list)) {
         if (verbosity >1) print(paste("File", i, catalog_list[i]))
         catalog_times[[i]] <- get_origin_and_times(catalog_list[i])[[2]]
         catalog_startdates[i] <- catalog_times[[i]][1]
         catalog_enddates[i] <- catalog_times[[i]][length(catalog_times[[i]])]
         dum1 <- length(which((catalog_times[[i]] >= start_date)&(catalog_times[[i]] <= end_date)))
         if (dum1 >0) {
           blank_length <- blank_length + dum1
           in_range[i] <- TRUE
         }
      }

      # Let's make sure these are in temporal order. Note that sort.chron() ignores index.return so we need to convert to numeric.
      ix <- sort(as.numeric(catalog_startdates[in_range]), index.return=TRUE)$ix
      ix <- which(in_range)[ix]
      catalog_list <- catalog_list[ix]
      catalog_times <- catalog_times[ix]
      catalog_startdates <- catalog_startdates[ix]
      catalog_enddates <- catalog_enddates[ix]
      icatalog <- 1
      input_file <- catalog_list[icatalog]
      ds <- catalog_times[[icatalog]]
      blank_length <- as.numeric(end_date - start_date + 1) / as.numeric(ds[2] - ds[1]) #+ 0.5/(as.numeric(ds[2] - ds[1]))
  } else {
      # We are looking at a single netcdf file such as a RECOM output file
      if (!is.na(eta_stem)) etafile <- paste0(eta_stem, '.nc')
      input_file <- input_file
      dum1 <- get_origin_and_times(input_file)
      ereefs_origin <- dum1[[1]]
      ds <- dum1[[2]]
      blank_length <- as.numeric(end_date - start_date + 1) / as.numeric(ds[2] - ds[1]) #+ 0.5/(as.numeric(ds[2] - ds[1]))
  }

  # Initialise
  blanks <- rep(NA, blank_length)
  ts_frame <- data.frame(as.Date(blanks), array(blanks, dim=c(length(blanks), length(var_names))))
  names(ts_frame) <- c("date", var_names)

  nc <- safe_nc_open(input_file)

  if (is.integer(location_latlon)) {
     location_grid <- location_latlon
  } else { 
    # Find the nearest grid-points to the sampling location
    if (is.null(nc$var[['latitude']])) {
      latitude <- safe_ncvar_get(nc, "y_centre")
      longitude <- safe_ncvar_get(nc, "x_centre")
    } else { 
      latitude <- safe_ncvar_get(nc, "latitude")
      longitude <- safe_ncvar_get(nc, "longitude")
    }
    tmp <- (latitude - location_latlon[1])^2 + (longitude - location_latlon[2])^2 
    tmp <- which.min(tmp) 
    #location_grid <- arrayInd(tmp, dim(latitude))
    location_grid <- c(floor((tmp+dim(latitude)[1]-1)/dim(latitude)[1]),
		       (tmp+dim(latitude)[1]-1)%%dim(latitude)[1] + 1)
  }
  zat <- ncdf4::ncatt_get(nc, "botz")
  if (!is.null(zat$positive)) {
	if (zat$positive=="down") zsign <- -1 else zsign <- 1
  } else {
	zsign <-1
  }
  botz <- zsign * as.numeric(safe_ncvar_get(nc, "botz", start=c(location_grid[2], location_grid[1]), count=c(1,1)))
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
     if (ereefs_case[2] == '4km') { 
        fileslist <- 1
        input_file <- paste0(input_stem, format(as.Date(paste(year, month, 1, sep="-")), '%Y-%m'), '.nc') 
        day_count <- day_count / as.numeric(ds[2]-ds[1])
     } else if (ereefs_case[2] == '1km') {
        fileslist <- from_day:(from_day+day_count-1)
        from_day <- 1
        day_count <- 1
     } else if ((ereefs_case[2] == 'recom')|(ereefs_case[1] == "ncml")) { 
       day_count <- day_count / as.numeric(ds[2]-ds[1])
       if (day_count > length(ds)) {
         warning(paste('end_date', end_date, 'is beyond available data. Ending at', ds[length(ds)]))
         day_count <- length(ds)
       }
       from_day <- (as.numeric(chron::chron(paste(year, month, from_day, sep = '-'), format=c('y-m-d'),
                                   origin=c(year=1990, month=1, day=1)) - ds[1]) + 
                               start_tod) / as.numeric(ds[2] - ds[1]) + 1 
	     if (from_day<1) from_day <-1
	     fileslist <- 1
     } else if (ereefs_case[1] == "thredds_catalog") { 
       month_startdate <- chron::chron(paste(year, month, 1, sep = '-'), format = 'y-m-d',
                                       origin=c(year=1990, month=1, day=1))
       month_enddate <- chron::chron(paste(year, month, daysIn(as.Date(paste(year, month, 1, sep='-'))) , sep = '-'), format = 'y-m-d',
                                       origin=c(year=1990, month=1, day=1))
       day_count <- day_count / as.numeric(ds[2]-ds[1])
       ix <- which((catalog_startdates >= month_startdate) & (catalog_enddates <= (month_enddate + 0.999)))
       icatalog <- ix[1]
       fileslist <- catalog_list[ix]
       if (end_date > catalog_enddates[length(catalog_enddates)]) {
         warning(paste('end_date', end_date, 'is beyond available data. Ending at', catalog_enddates[length(catalog_enddates)]))
         day_count <- day_count - (end_date - catalog_enddates[length(catalog_enddates)])
       }
       from_day <- (as.numeric(chron::chron(paste(year, month, from_day, sep = '-'), format='y-m-d',
                               origin=c(year=1990, month=1, day=1)) - catalog_startdates[ix[1]]) 
                    + start_tod) / as.numeric(ds[2] - ds[1]) + 1 
       if (from_day<1) from_day <- 1
     } else stop("Shouldn't happen: ereefs_case not recognised")

     for (dcount in fileslist) {
        if (ereefs_case[2] == '1km') {
	        input_file <- paste0(input_stem, format(as.Date(paste(year, month, fileslist[dcount], sep="-")), '%Y-%m-%d'), '.nc')
        } else if (ereefs_case[1] == "thredds_catalog") {
          icatalog <- icatalog[dcount]
          input_file <- catalog_list[icatalog]
          ds <- catalog_times[[icatalog]]
        }
        #input_file <- paste0(input_file, '?', var_list, ',time,eta')
        nc <- safe_nc_open(input_file)
	      if (!is.na(eta_stem)) nc3 <- safe_nc_open(etafile)
        if ((ereefs_case[2] == "1km")||(ereefs_case[2] == "4km")) {
          if (!is.null(nc$var[['t']])) {
            d <- as.Date(safe_ncvar_get(nc, "t"), origin = as.Date("1990-01-01"))[from_day:(from_day+day_count-1)]
          } else {
            d <- as.Date(safe_ncvar_get(nc, "time"), origin = as.Date("1990-01-01"))[from_day:(from_day+day_count-1)]
          }
        } else { 
           d <- ds[from_day:(from_day + day_count - 1)]
        }
        if (!is.null(nc$var[['eta']])) { 
          eta <- safe_ncvar_get(nc, 'eta', start=c(location_grid[2], location_grid[1],from_day), count=c(1,1,day_count)) 
        } else {
          eta <- safe_ncvar_get(nc3, 'eta', start=c(location_grid[2], location_grid[1],from_day), count=c(1,1,day_count)) 
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
           wc <- safe_ncvar_get(nc, var_names[j], start=c(location_grid[2],location_grid[1],1,from_day), count=c(1,1,-1,day_count)) 
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

#' A wrapper to ncdf4::ncvar_get() that will pause and try again several times (defaulting to 12)
#' if it is a web-served netcdf file and at first it fails, to overcome temporary net access errors or DAP errors.
#'
#' Parameters before 'tries' are passed through to ncvar_get
#'
#' @param tries number of times to retry (increasing pause length by one second each time. Default 4 
#' @return variable extracted using ncvar_get()
#' @export
safe_ncvar_get <- function(nc,varid=NA, start=NA, count=NA, verbose=FALSE,
 signedbyte=TRUE, collapse_degen=TRUE, raw_datavals=FALSE, tries=11) {
   if (substr(nc$filename, 1, 4)!="http") {
     myvar <- ncdf4::ncvar_get(nc, varid, start, count, verbose, signedbyte, collapse_degen, raw_datavals)
   } else {
     myvar <- try(ncdf4::ncvar_get(nc, varid, start, count, verbose, signedbyte, collapse_degen, raw_datavals))
     trywait = 1
     while ((class(myvar)=='try-error')&&(trywait<=(tries+1))) { 
        print(paste('retrying in ', trywait, 'second(s)')) 
        Sys.sleep(trywait) 
        trywait <- trywait+1 
        myvar <- try(ncdf4::ncvar_get(nc, varid, start, count, verbose, signedbyte, collapse_degen, raw_datavals))
     }
     if (trywait>(tries+1)) stop(paste('Cannot access netcdf file', nc$filename))
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
    while ((class(nc)=='try-error')&&(trywait<=(tries+1))) { 
       warning(paste('This probably means the netcdf file name is incorrect or target date out of range for this eReefs run,\n',
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
