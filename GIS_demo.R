#######################################################################
##################### GIS in R - R Seminar Fall 2015 ##################
#######################################################################
#install.packages('ggplot2')
#install.packages('sp')
#install.pacakges('rgdal')
#install.packages('raster')
#install.packages('rgeos')
#install.packages('RColorBrewer')
library(raster);library(ggplot2);library(RColorBrewer);library(rgeos);library(data.table)

################################ Shapefiles, rasters, and simple GIS data summaries in R ##################################
#R is great for *some* GIS tasks. If you're exploring data and want to quickly demo a visalization or click around a map, try
#qGIS or arcMap. If you're working with a bunch of files, want to automate repetitive tasks, or just love ggplot2, R is
#a good choice.

#Reading in shapefiles (range maps of the two most common hummingbirds in Pacific NW - Anna's and Rufous)
anhu <- shapefile("./ranges/Calypte_anna_22688199.shp")
ruhu <- shapefile("./ranges/Selasphorus_rufus_22688296.shp")

#check out the data
summary(anhu)
head(anhu@data)

#plot with base graphics
plot(anhu)

#Looks weird. The map has separate polygons for breeding, postbreeding, and wintering ranges in the same shapefile. 
#Filter polygons by the associated data with the "[ ]" function.  
anhu.b <- anhu[anhu@data$SEASONAL == 1,]
ruhu.b <- ruhu[ruhu@data$SEASONAL == 2,]
plot(anhu.b)

#standard vector functions are all available through the rgeos package
plot(union(anhu.b,ruhu.b))
plot(intersect(anhu.b,ruhu.b))

#change projections with spTransform. 3395 is a mercator projection, 4087 is equidistant cylindrical. 
anhu.proj <- spTransform(anhu.b, CRS("+init=epsg:3395"))
ruhu.proj <- spTransform(ruhu.b,CRS("+init=epsg:3395"))

#Buffering. If the spatial object is projected in a lat/long system, the unit of width is meters. Otherwise, map units. 
plot(buffer(anhu.proj,width=1e4))

#check out the area of some shapefiles
gArea(buffer(anhu.proj,width=1e4))
gArea(anhu.proj)/gArea(ruhu.proj)

#######Plotting shapefiles with ggplot2: 
#Convert sp objects to dataframes with fortify()
anhu.df <- fortify(anhu.b)
ruhu.df <- fortify(ruhu.b)
#load in a country outlines file and crop to the plot extent 
map <- shapefile("cntry06/cntry06.shp")
admin <- crop(map,c(-150,-105,25,65))
admin.df <- fortify(admin)
#plot with geom_polygon()
ggplot()+coord_map(xlim=c(-150,-105),ylim=c(25,65))+theme_bw()+theme(panel.grid=element_blank())+
  geom_polygon(data=anhu.df,aes(x=long,y=lat,group=group),fill="violet")+
  geom_polygon(data=ruhu.df,aes(x=long,y=lat,group=group),fill="orangered",alpha=0.7)+
  geom_polygon(data=admin.df,aes(x=long,y=lat,group=group),fill=NA,col="black")
#pretty! 

#################################################################################
##################### Loop 1: Regional species diversity heatmap ################
#################################################################################
#let's build a raster where the cell value equals the number of species within a clade that occur in a given area. 
##build an empty raster on a lat/long grid at 10min resolution
r <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, res=1/6, vals=0)
west <- crop(r,c(-150,-30,-42,55))
#test out rasterizing with just the anhu range.
anhu.r <- rasterize(anhu.b,west)
plot(anhu.r)

#summarizing multiple shapefiles: loop over all the shapefiles in the rangeMaps folder, rasterize, add the rasters.  
files <- list.files("./ranges/")
files <- files[grep(".shp",files)]
n.species <- crop(r,c(-150,-30,-42,55))

for (i in files){
  range <- shapefile(paste("./ranges/",i,sep=""))
  range.r <- rasterize(range,west,1,background=0)
  n.species <- range.r+n.species
}

plot(n.species)
#writeRaster(n.species,"~/Desktop/beeHum_div_10min.tif")


##############################################################################################################################
############## loop #2: visualizing hummingbird migration ####################################################################
### See Supp 2015 (http://www.esajournals.org/doi/abs/10.1890/ES14-00290.1) for a better version of this #####################
##############################################################################################################################

#read in rufous hummingbird eBird reports and a country outlines map
ruhu <- read.delim("ebd_rufhum_relNov-2014/ebd_rufhum_relNov-2014.txt")
ruhu$date <- as.Date(ruhu$OBSERVATION.DATE,"%m/%d/%Y")
ruhu$month <- as.numeric(substr(ruhu$date,6,7))
ruhu$monthName <- months(ruhu$date)
ruhuW <- subset(ruhu,LONGITUDE < -96)
ext <- extent(c(-145,-60,10,62))
map <- shapefile("cntry06/cntry06.shp")

#read in an effort raster (total # reports per grid cell at 1 deg resolution)
r <- raster(xmn=-180, xmx=180, ymn=-90, ymx=90, res=1, vals=0)
effort <- crop(resample(raster("effort.tif"),r),ext)
r.ruhu <- crop(r,ext)

##loop over months: subset occurrence records, rasterize, divide by effort, write to file
#as a static figure
old.par <- par()
pdf(file="./ruhu_freq/ruhu_monthly.pdf",width=6.5,height=5)
par(mfrow=c(3,4),mar=c(2.5,1,1,1),mgp=c(1.5, .5, 0),oma=c(0,0,3,0))
for(i in c(1:12)){
  a <- subset(ruhu,month == i)
  locs <- SpatialPoints(data.frame(a$LONGITUDE,a$LATITUDE))
  frequency <- rasterize(locs,r.ruhu,fun="count")/effort
  plot(frequency,col=brewer.pal(n=7,name="YlOrRd"),legend=F,axes=F,breaks=c(1e-2,5e-2,1e-1,5e-1,1,5,10),
       main=paste(a$monthName[1]),xaxs="i", yaxs="i")+
          plot(map,col=NA,add=TRUE)+
          mtext("Rufous Hummingbird Report Frequency",outer = TRUE,side = 3,cex = 1.2,line = 1)
}
dev.off()
par(old.par)

#or output new pdf's for each month, then gif it up in photoshop
for(i in c(1:12)){
  a <- subset(ruhu,month == i)
  locs <- SpatialPoints(data.frame(a$LONGITUDE,a$LATITUDE))
  frequency <- rasterize(locs,r.ruhu,fun="count")/effort
  pdf(paste("./ruhu_freq/freq_",i,".pdf",sep=""))
  plot(frequency,col=brewer.pal(n=7,name="YlOrRd"),legend=F,axes=F,breaks=c(1e-2,5e-2,1e-1,5e-1,1,5,10),
       main=paste(a$monthName[1]),xaxs="i", yaxs="i")+
    plot(map,col=NA,add=TRUE)+
    mtext("Rufous Hummingbird Report Frequency",outer = TRUE,side = 3,cex = 1.2,line = 1)
  dev.off()
}






