library(ereefs)

# 1. Making maps and animations from eReefs output
#
# Create a true colour image map from optical model output of the GBR on a specific date (select menu option 9 for high-resolution or 11 for lower resolution):
p <- map_ereefs(target_date = c(2019, 2, 28))

# Map optical colour classes (note the default colour scheme isn't great for this, so you might want to add a different colour scale after producing the figure handle):
p <- map_ereefs(var_name = "plume", target_date = c(2019, 2, 28))

# Map the extent of the influence of the Burdekin River (This time, choose MENU OPTION 8 or 7):
# THE FOLLOWING OUGHT TO WORK BUT IS CURRENTLY BROKEN. This will be an issue with the parsing of the filenames provided by the catalog. Need to fix.
p <- map_ereefs(var_name = "bur", target_date = c(2022, 2, 28))

# Meanwhile, manually specifying a single netcdf inut file name does still work:
p <- map_ereefs(var_name = "bur", target_date = c(2022, 2, 28), input_file="https://dapds00.nci.org.au/thredds/dodsC/fx3/gbr1_2.0_rivers/gbr1_rivers_simple_2022-08-28.nc.html")

# Change the limits of the colour scale so that the max colour intensity is achieved when 1% of water is river water, and change the colours:
p <- map_ereefs(var_name = "bur", target_date = c(2022, 2, 28), scale_lim = c(0, 0.01), input_file="https://dapds00.nci.org.au/thredds/dodsC/fx3/gbr1_2.0_rivers/gbr1_rivers_simple_2022-08-28.nc.html")
#p <- map_ereefs(var_name = "burdekin", target_date = c(2019, 2, 28), scale_lim = c(0, 0.1), scale_col = c('ivory', 'dodgerblue'))

# Map total chlorophyll a using a spectral colour scale (MENU OPTION 11):
p <- map_ereefs(var_name = "Chl_a_sum", target_date = c(2019, 2, 28), scale_col = "spectral")

# Map total chlorophyll 5m below MSL instead of at the surface:
p <- map_ereefs(var_name = "Chl_a_sum", target_date = c(2019, 2, 28), scale_col = "spectral", layer = -5)

# Plot surface chlorophyll again, but this time add a land map underlay:
p <- map_ereefs(var_name = "Chl_a_sum", target_date = c(2019, 2, 28), scale_col = "spectral", Land_map = TRUE)

# Zoom in to the area near Townsville (use MENU OPTION 9 for higher resolution):
p <- map_ereefs(var_name = "Chl_a_sum", target_date = c(2019, 2, 28), scale_col = "spectral", Land_map = TRUE, box_bounds = c(145, 150, -22, -18))

# Create the image files needed for a true colour animation of surface salinity in this area (then use animate.sh or similar to generate the animation itself):
# (salt_list is a list that includes the plot handle as well as the temporal mean salinity for each point in the map over the period of the animation, plus the cell centre geolocations)
# By default, images are saved to the directory ToAnimate. (MENU OPTION 9)
salt_list <- map_ereefs_movie(var_name = "salt", start_date = c(2019, 2, 15), end_date = c(2019, 3, 10), Land_map = TRUE, box_bounds = c(145, 150, -22, -18), scale_col = "spectral", scale_lim = c(30, 35))

# 2. Plotting vertical profiles and slices through the data
#
# Extract the data from a vertical slice of temperature along a line segment at latitude 20 degrees South (select MENU OPTION 9):
# Note that a slice can also be defined along a long, curvy line, for example a boat track, by adding more points to location_latlon. You can
# also get data from multiple variables by making var_names a vector (e.g., var_names = c("temp", "salt")):
temp_slice <- get_ereefs_slice(var_names = "temp", target_date = c(2022, 8, 1), location_latlon = data.frame(latitude = c(-20, -20), longitude=c(145, 150))) 

# Visualise the results:
# Todo: add "spectral" colour scale option to this function
p <- plot_ereefs_slice(temp_slice, var_name="temp")

# Extract a vertical profile of chlorophyll a and nitrate at a single location and time:
# THIS DOESN"T CURRENTLY WORK BECAUSE I HAVEN'T YET UPDATED THE DATE FORMAT TO CHRON IN THIS FUNCTION
  profile_data <- get_ereefs_profile(var_names = c("Chl_a_sum", "NO3"), start_date = c(2019, 2, 10), end_date = c(2019, 2, 10), location_latlon = c(-23.39, 150.89))

  # Visualise the results for NO3:
  p <- plot_ereefs_profile(profile_data, var_name = "NO3", target_date = c(2019, 2, 10))
  
  # Exctract a vertical profile of TN at the same location for each date in September 2022:
  profile_data <- get_ereefs_profile(var_names = "TN", start_date = c(2022, 9, 1), end_date = c(2022, 9, 30), location_latlon = c(-23.39, 150.89))
  
  # Visualise the results:
  p <- plot_ereefs_zvt(profile_data, var_name = "TN")

# Extract a time-series of data

