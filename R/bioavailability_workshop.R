library(ggplot2)
source('../../multiplot.R')
load('sysdata.rda')
load('filenames.rda')
d <- c(2017,4,5) # Post Cyclone Nathan
box_bounds <- c(144,152,-25,-17) # North of Rockhampton
transect <- data.frame(latitude=c(-21.05, -21.5), longitude=c(149.8, 149.3))
salt <- map_ereefs(var_name='salt', input_file=gbr4hd, target_date=d, box_bounds=bb, zoom=7, scale_col=c("ivory", "royalblue"), scale_lim=c(32,35))
salt <- map_ereefs(var_name='salt', input_file=gbr1hd, p=salt, target_date=d, box_bounds=bb, zoom=7, scale_col=c("ivory", "royalblue"), scale_lim=c(32,35))
salts <- get_ereefs_slice(var_name='salt', input_file=gbr1hd, target_date=d, location_latlon=transect)
saltslice <- plot_ereefs_slice(salts, var_name='salt', scale_col=c("ivory", "royalblue"), scale_lim=c(32,35)) + xlim(0,40) + ylim(-25,1)
TN <- map_ereefs(var_name='TN', input_file=gbr4bgc, target_date=d, box_bounds=bb, zoom=7, scale_col=c("ivory", "royalblue"), scale_lim=c(0,200))
TN <- map_ereefs(var_name='TN', input_file=gbr1bgc, p=TN, target_date=d, box_bounds=bb, zoom=7, scale_col=c("ivory", "royalblue"), scale_lim=c(0,200))

s <- get_ereefs_slice(var_name=c('PhyS_N', 'PhyL_N', "Tricho_N", "ZooL_N", "ZooS_N", 'TN', 'DIN', 'NO3', 'NH4', 'DOR_N', 'DetPL_N', 'DetBL_N', 'DetR_N'), input_file=gbr1bgc, target_date=d, location_latlon=transect, eta_stem=gbr1hd, override_positive=TRUE)

        s2 <- s
	dum1 <- s$values
	dum1[,,1] <- dum1[,,1] + dum1[,,2] + dum1[,,3]
	dum1[,,4] <- dum1[,,4] + dum1[,,5]
	s2$values <- dum1
	dimnames(s2$values)[[3]][1] <- "Phy_N"
	dimnames(s2$values)[[3]][4] <- "Zoo_N"

tnslice <- plot_ereefs_slice(s, var_name='TN', scale_col=c("ivory", "coral4"), scale_lim=c(0,200)) + xlim(0,40) + ylim(-25,1)
Physlice <- plot_ereefs_slice(s2, var_name='Phy_N', scale_col=c("ivory", "forestgreen"), scale_lim=c(2,15)) + xlim(0,40) + ylim(-25,1)
Zooslice <- plot_ereefs_slice(s2, var_name='Zoo_N', scale_col=c("ivory", "purple4"), scale_lim=c(2,15)) + xlim(0,40) + ylim(-25,1)
DINslice <- plot_ereefs_slice(s, var_name='DIN', scale_col=c("ivory", "coral4"), scale_lim=c(0,20)) + xlim(0,40) + ylim(-25,1)
DONslice <- plot_ereefs_slice(s, var_name='DOR_N', scale_col=c("ivory", "chocolate4"), scale_lim=c(0,200)) + xlim(0,40) + ylim(-25,1)
DetPLslice <- plot_ereefs_slice(s, var_name='DetPL_N', scale_col=c("ivory", "coral4"), scale_lim=c(0,200)) + xlim(0,40) + ylim(-25,1)
NH4slice <- plot_ereefs_slice(s, var_name='NH4', scale_col=c("ivory", "coral4"), scale_lim=c(0,200)) + xlim(0,40) + ylim(-25,1)
NO3slice <- plot_ereefs_slice(s, var_name='NO3', scale_col=c("ivory", "coral4"), scale_lim=c(0,200)) + xlim(0,40) + ylim(-25,1)
PhyLslice <- plot_ereefs_slice(s, var_name='PhyL_N', scale_col=c("ivory", "coral4"), scale_lim=c(0,200)) + xlim(0,40) + ylim(-25,1)
png('slices.png', width=1000, height=1200)
multiplot(saltslice, DINslice, Physlice, Zooslice, DONslice)
dev.off()

DOR_N <- map_ereefs(var_name='DOR_N', input_file=gbr4bgc, target_date=d, box_bounds=bb, zoom=7, scale_col=c("ivory", "chocolate4"), scale_lim=c(0,100))
DOR_N <- map_ereefs(var_name='DOR_N', input_file=gbr1bgc, p=DOR_N, target_date=d, box_bounds=bb, zoom=7, scale_col=c("ivory", "chocolate4"), scale_lim=c(0,200))
NH4 <- map_ereefs(var_name='NH4', input_file=gbr4bgc, target_date=d, box_bounds=bb, zoom=7, scale_col=c("ivory", "coral4"), scale_lim=c(0,20))
NH4 <- map_ereefs(var_name='NH4', input_file=gbr1bgc, p=NH4, target_date=d, box_bounds=bb, zoom=7, scale_col=c("ivory", "coral4"), scale_lim=c(0,20))
NO3 <- map_ereefs(var_name='NO3', input_file=gbr4bgc, target_date=d, box_bounds=bb, zoom=7, scale_col=c("ivory", "coral4"), scale_lim=c(0,20))
NO3 <- map_ereefs(var_name='NO3', input_file=gbr1bgc, p=NO3, target_date=d, box_bounds=bb, zoom=7, scale_col=c("ivory", "coral4"), scale_lim=c(0,20))
DIN <- map_ereefs(var_name='DIN', input_file=gbr4bgc, target_date=d, box_bounds=bb, zoom=7, scale_col=c("ivory", "coral4"), scale_lim=c(0,40))
DIN <- map_ereefs(var_name='DIN', input_file=gbr1bgc, p=DIN, target_date=d, box_bounds=bb, zoom=7, scale_col=c("ivory", "coral4"), scale_lim=c(0,40))
Phy <- map_ereefs(var_name='PhyL_N', input_file=gbr4bgc, target_date=d, box_bounds=bb, zoom=7, scale_col=c("ivory", "forestgreen"))
Phy <- map_ereefs(var_name='PhyL_N', input_file=gbr1bgc, p=Phy, target_date=d, box_bounds=bb, zoom=7, scale_col=c("ivory", "forestgreen"))
Zoo <- map_ereefs(var_name='ZooL_N', input_file=gbr4bgc, target_date=d, box_bounds=bb, zoom=7, scale_col=c("ivory", "purple4"))
Zoo <- map_ereefs(var_name='ZooL_N', input_file=gbr1bgc, p=Zoo, target_date=d, box_bounds=bb, zoom=7, scale_col=c("ivory", "purple4"))
true <- map_ereefs(input_file=gbr4bgc, target_date=d, box_bounds=bb, zoom=7)
true <- map_ereefs(input_file=gbr1bgc, p=true, target_date=d, box_bounds=bb, zoom=7)
plume <- map_ereefs(var_name='plume', input_file=gbr4bgc, target_date=d, box_bounds=bb, zoom=7, scale_col=c("purple4", "ivory"), scale_lim=c(0,6))
plume <- map_ereefs(var_name='plume', input_file=gbr1bgc, p=plume, target_date=d, box_bounds=bb, zoom=7, scale_col=c("purple4", "ivory"), scale_lim=c(0,6))

pNng('maps.png', width=1000, height=1200)
multiplot(salt, DIN, Phy, Zoo, DOR_N, cols=2)
dev.off()

map_ereefs_movie(var_name='DIN', input_file=gbr1bgc, start_date=c(2017,3,15), end_date=c(2017,4,15), scale_lim=c(0,40), box_bounds=bb, scale_col=c("ivory","coral4"), zoom=7)
map_ereefs_movie(var_name='PhyL_N', input_file=gbr1bgc, start_date=c(2017,3,15), end_date=c(2017,4,15), scale_lim=c(0,20), box_bounds=bb, scale_col=c("ivory","forestgreen"), zoom=7)
map_ereefs_movie(input_file=gbr1bgc, start_date=c(2017,3,15), end_date=c(2017,4,15), box_bounds=bb, zoom=7)
