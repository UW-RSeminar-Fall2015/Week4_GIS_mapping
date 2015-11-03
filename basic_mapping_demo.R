###########################
#### FUN WITH MAPS!!!! ####
###########################

### what you'll need

getwd()
install.packages('ggplot2')
install.packages('ggmap')
install.packages('mapdata')
install.packages('maps')
library(ggplot2)
library(ggmap)
library(mapdata)
library(maps)

################################
#### making basic maps in R ####
################################

map("worldHires", "Mexico") #pick your basemap, define your country
map("worldHires", "Mexico", col="grey90", fill=TRUE) #example visual tweak
map("worldHires", xlim=c(-130, -53), ylim=c(15,35)) #define by lat / long instead

localities <- read.csv("PABU.csv") #format is a column "lat" and a column "long" with data in decimal degrees
colnames(localities) <-c("num","lat","long","species")
map("worldHires", "Mexico", col="grey90", fill=TRUE)
points(localities$long, localities$lat, pch=19, col="red", cex=0.5) #plot localities data, choose aesthetic parameters

##################################################
#### more of the same but better with ggplot2 ####
##################################################

map <- map_data("world", "Mexico") #pick basemap -- higher res options available 
ggplot() + theme_bw() + geom_path(data=map, aes(x=long, y=lat, group=group)) + coord_map() #look at the basemap
ggplot() + geom_point(data=localities, aes(x=long, y=lat)) #look at the points in space

#### and together now!
ggplot() +
  geom_path(data=map, aes(x=long, y=lat, group=group)) +
  geom_point(data=localities, aes(x=long, y=lat)) 

ggplot() + coord_map()+ #hold ratios / project constant
  geom_path(data=map, aes(x=long, y=lat, group=group)) +
  geom_point(data=localities, aes(x=long, y=lat)) 

#### add more graphical parameters
ggplot() + coord_map()+
  geom_path(data=map, aes(x=long, y=lat, group=group)) +
  geom_point(data=localities, aes(x=long, y=lat, col=species, size=num)) #color code points by species, scale by size

ggplot() + coord_map()+
  geom_path(data=map, aes(x=long, y=lat, group=group)) +
  geom_point(data=localities, aes(x=long, y=lat, col=species, size=num)) +
  scale_size_continuous(range = c(3,13)) #tweak acceptable range of point sizes

ggplot() + coord_map()+
  geom_path(data=map, aes(x=long, y=lat, group=group)) +
  geom_point(data=localities, aes(x=long, y=lat, col=species, size=num)) +
  scale_size_continuous(range = c(3,13)) + 
  theme_bw() # remove greyscale background

ggplot() + coord_map()+
  geom_path(data=map, aes(x=long, y=lat, group=group)) +
  geom_point(data=localities, aes(x=long, y=lat, col=species, size=num)) +
  scale_size_continuous(range = c(3,13)) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #remove background grid

ggplot() + coord_map()+
  geom_path(data=map, aes(x=long, y=lat, group=group)) +
  geom_point(data=localities, aes(x=long, y=lat, col=species, size=num)) +
  scale_size_continuous(range = c(3,13)) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.title = element_blank(),legend.text = element_text(face = "italic"), axis.title.x = element_blank(), axis.title.y = element_blank()) #remove labels

#### a few other tricks

ggplot() + coord_map()+
  geom_path(data=map, aes(x=long, y=lat, group=group)) +
  #geom_point(data=localities, aes(x=long, y=lat, col=species, size=num)) +
  #scale_size_continuous(range = c(3,13)) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.title = element_blank(),legend.text = element_text(face = "italic"), axis.title.x = element_blank(), axis.title.y = element_blank()) + 
  geom_bin2d(data=localities,aes(x=long,y=lat)) #"rasterize" your data (record density)

ggplot() + coord_map()+
  geom_path(data=map, aes(x=long, y=lat, group=group)) +
  #geom_point(data=localities, aes(x=long, y=lat, col=species, size=num)) +
  #scale_size_continuous(range = c(3,13)) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.title = element_blank(),legend.text = element_text(face = "italic"), axis.title.x = element_blank(), axis.title.y = element_blank()) + 
  stat_summary2d(data=localities,aes(x=long,y=lat,z=num,fun="mean")) # visualize a summary statistic of it 

ggplot() + coord_map()+
  geom_path(data=map, aes(x=long, y=lat, group=group)) +
  #geom_point(data=localities, aes(x=long, y=lat, col=species, size=num)) +
  #scale_size_continuous(range = c(3,13)) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.title = element_blank(),legend.text = element_text(face = "italic"), axis.title.x = element_blank(), axis.title.y = element_blank()) + 
  stat_density2d(data=localities,aes(x=long,y=lat)) ### cool topographic overlay of density 





