# Run these two lines the first time if you don't already have the ereefs and ncdf4 libraries:
# install.packages('devtools', 'ncdf4')
# devtools::install_github('BarbaraRobson/ereefs')
library('ereefs')
library('ncdf4')
#load('sysdata.rda')
input_file <- substitute_filename('GBR4-v2.0')
nc <- nc_open(input_file)
botz <- ncvar_get(nc, 'botz')
latitude <- ncvar_get(nc, 'latitude')
longitude <- ncvar_get(nc, 'longitude')
nc_close(nc)

times_avail <- as.Date(as.Date('2011-01-01'):as.Date('2018-10-15'), origin='1970-01-01')

# dimension is 600x180. Avoid the bottom 7 rows and the right-hand 7 rows to move away from the open boundary.
max_x <- dim(botz)[2] - 7
max_y <- dim(botz)[1] - 7

n <- 1000 # target number of samples
n2 <- n * 1.5 # sample slightly more than the target so we can remove duplicates and dry cells

sample_points <- data.frame(x = sample(1:max_x, n2, replace = TRUE), 
                            y = sample(1:max_y, n2, replace = TRUE), 
                            d = sort(sample(times_avail, n2, replace = TRUE)))
# remove duplicates
sample_points <- sample_points[!duplicated(sample_points),]
# remove dry cells
depth = c(botz)[sample_points$y + (sample_points$x - 1) * dim(botz)[1]]
sample_points <- sample_points[!is.na(depth), ][1:n,]
depth <- depth[!is.na(depth)][1:n]

ereefs_data <- data.frame(date = sample_points$d,
                          depth = depth,
                          latitude = c(latitude)[sample_points$y + (sample_points$x - 1) * dim(latitude)[1]],
                          longitude = c(longitude)[sample_points$y + (sample_points$x - 1) * dim(longitude)[1]],
                          Tricho_Chl = NA*depth,
                          Chl_a_sum = NA*depth,
                          temperature = NA*depth,
                          TN = NA*depth,
                          TP = NA*depth,
                          PAR = NA*depth)
var_names <-c('Tricho_Chl', 'Chl_a_sum', 'temp', 'TN', 'TP', 'PAR')

if (any(ereefs_data$date > as.Date('2016-10-31'))) { 
   if (any(ereefs_data$date <= as.Date('2016-10-31'))) {
     arch <- 1:max(which(ereefs_data$date <= as.Date('2016-10-31')))
     nrt <- (arch[length(arch)] + 1):n
   } else {
     arch <- NULL
     nrt <- 1:n
   }
} else {
   arch <- 1:n
   nrt <- NULL
}


pb <- txtProgressBar(min = 0, max = n, style = 3)
#for (i in 11571:length(arch)) {
for (i in arch) {
   ereefs_data[i, 5:dim(ereefs_data)[2]] <-  get_ereefs_ts(var_names, 
                                                          location_latlon = c(sample_points$x[i], sample_points$y[i]),
                                                          start_date = ereefs_data$date[i], 
                                                          end_date = ereefs_data$date[i], 
                                                          input_file = 2,
                                                          verbosity = 0)[,2]
        setTxtProgressBar(pb,i)
}
for (i in nrt) {
#for (i in 12845:nrt[length(nrt)]) {
   #ereefs_data[i, 5:dim(ereefs_data)[2]] <-  get_ereefs_depth_integrated_ts(var_names, 
   ereefs_data[i, 5:dim(ereefs_data)[2]] <-  get_ereefs_ts(var_names, 
                                                          location_latlon = c(sample_points$x[i], sample_points$y[i]),
                                                          start_date = ereefs_data$date[i], 
                                                          end_date = ereefs_data$date[i], 
                                                          input_file = 3,
                                                          verbosity = 0)[,2]
        setTxtProgressBar(pb,i)
}
close(pb)
