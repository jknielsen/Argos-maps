#Makes a map of most likely location of PSAT stranded on shore
#Works by rasterizing individual error ellipses from Argos and
#assigning a precision estimate to each (higher values to better estimates).
#Then all rasters are summed so that grid cells with numerous good estimates will have a higher 
#value than grid cells with poor estimates.
#Outputs a raster that can be exported to GIS for overlay on imagery, making maps, etc.

#Julie Nielsen
#julie.nielsen@gmail.com
#Kingfisher Marine Research, LLC
#updated June 2018


#Libraries required
library(sp)
library(spatstat)
library(plotrix)
library(raster)
library(rgdal) 

#Export data from Argos, but first:

## 1. *** Make sure you check the box for "Diagnostic data" before you download! **** ##
## 2. Check the box for "Most significant or best on satellite pass"
## 3. Make sure records are sorted by increasing date before you export
## 4. Save as .csv

#Example data are provided here for a tag washed ashore in a Norwegian fjord

#Import data:
argos<- read.csv("ArgosData_Norway.csv")
str(argos)

#Select the useful columns
argos<- argos[,c(1,4,5,6,9,10,11,15,17,18,19,20,28:32)]

#By column name, in case order is changed in the future
#Column names may differ, though, so by column # is probably better
# argos<- argos[,c("Platform.ID.No.","Pass.dur...s.","Msg.Date","Sat.","Loc..date","Longitude","Latitude",
#                  "Loc..quality","Long..1","Lat..sol..1","Long..2","Lat..sol..2",
#                  "Error.radius","Semi.major.axis","Semi.minor.axis",
#                  "Ellipse.orientation","GDOP")]


dimnames(argos)[[2]]<- c("ID","duration","msg_date","sat","loc_date","Long","Lat",
                         "qual","long1","lat1","long2","lat2",
                        "error_radius","major_axis","minor_axis","orientation","GDOP")

#Pick out individual tag if desired, not necessary if only one tag was downloaded from Argos
#argos<- argos[argos$ID==129847,]

#get rid of duplicates
# This step is not necessary if the "most significant or best on satellite pass" box was checked for export
#argos<- argos[!duplicated(argos$GDOP),]#risk of throwing out entries with the same GDOP; check to make sure!

#get rid of NAs
argos<- argos[!is.na(argos$GDOP),]

#Assign record number
argos$Record<- 1:nrow(argos)

#Add UTM locs to argos records
#For Norway tag:
proj<- "+proj=tmerc +lat_0=58 +lon_0=6.166666666666667000 +k=1 +x_0=0 +y_0=0 +a=6377492.018 +b=6356173.508712696 +units=m +no_defs " 

argos[,c("X","Y")]<- project(cbind(argos$Long,argos$Lat),proj=proj)

#Change from semi for ellipse function
argos$new_maj<- argos$major_axis*2
argos$new_min<- argos$minor_axis*2

#remove records before tag was beached
rec.no<- 60 #Tag was floating for locations up to 60

plot(argos$X[argos$Record>rec.no],argos$Y[argos$Record>rec.no])

plot(argos$X,argos$Y)

#reselect so that only beached records are included in the following procedures
argos<- argos[argos$Record>rec.no,]

#Assign a score to each record (e.g., weighting by quality)
#You want the lowest GDOP value (interpreted roughly in meters)
#You also want the lowest error radius (also interpreted in meters)
#I tend to trust the GDOP more because the error estimates are ellipses, not radii
#A good estimate should have a very low GDOP and a very low error radius

#This is somewhat of an art now,but I look at the histograms and see if there are natural groups of
#GDOP values that are better than others
hist(argos$GDOP)
hist(argos$GDOP[argos$GDOP<1000])#zoom in
hist(argos$error_radius)
hist(argos$error_radius[argos$error_radius<1000]) 

argos$score[argos$GDOP>1000]<- 1
argos$score[argos$GDOP<1000]<- 10
argos$score[argos$GDOP<600]<- 20		
argos$score[argos$GDOP<300]<- 30 
argos$score[argos$error_radius<300]<- 30

hist(argos$score)

#Look at records with score to make sure everything makes sense
argos[,c("Lat","Long","qual","loc_date","error_radius","GDOP","score","Record")]

#Use geometric mean of best estimates to assign grid dimensions and find center of grid
table(argos$qual)

#loc.best<- argos[argos$qual=="3",] #If there are a lot of 3's
loc.best<- argos[argos$qual=="3"|argos$qual=="2",]
#loc.best<- argos[argos$qual=="1",]#Typically few 3's with Desert Star

x.gm<- sqrt(sum(loc.best$X^2/nrow(loc.best)))
y.gm<- sqrt(sum(loc.best$Y^2/nrow(loc.best)))

#Check it on the locations plot
plot(argos$X,argos$Y)
points(x.gm,y.gm,col="red",pch=19)

#Specify grid that encompasses potential location of stranded tag
dist<- 1000
from.x<-   x.gm-dist 
to.x<-     x.gm+dist  
from.y<-   y.gm-dist
to.y<-     y.gm+dist
grid.spacing<- 10
area.grid <- expand.grid(seq(from.x,to.x,grid.spacing), seq(from.y,to.y,grid.spacing))
point.x<- area.grid[,1]
point.y<- area.grid[,2]
str(area.grid)

#Make polygons for each record
# ** Note: run function for ellipseSP (below) first
x<- argos$X
y<- argos$Y
a<- argos$minor_axis
b<- argos$major_axis
phi<- 180-argos$orientation
ellipses<- lapply(1:nrow(argos), function(i) ellipseSP(x[i], y[i], a[i], 
            b[i], phi = phi[i])) 

#Check out ellipses with points - if they don't all overlap, some may still be drifting toward shore
#And additional points should be removed by rec.no
plot(argos$X,argos$Y)
for(i in 1:nrow(argos)){
  plot(ellipses[[i]],add=T)
}

#Make array
grid.array= array(0, dim=c(length(seq(from.x,to.x,grid.spacing)),length(seq(from.y,to.y,grid.spacing)),nrow(argos)))

#For each polygon, find grid cells that intersect
for (i in 1:length(ellipses)){
   poly.coords<- ellipses[[i]]@polygons[[1]]@Polygons[[1]]@coords
   poly.x<- poly.coords[,1]
   poly.y<- poly.coords[,2]
   pp<- point.in.polygon(point.x,point.y,poly.x,poly.y)  
   pp[pp>0]<- argos[i,"score"]
   grid.array[,,i]<- pp
}

#Sum all grids in array
grid.sum<- apply(grid.array,1:2,sum)

#Quick view
image(grid.spacing* (1:nrow(grid.sum)), grid.spacing* (1:ncol(grid.sum)), grid.sum, col = topo.colors(12))

#Export as raster - can be imported into ArcGIS
rstr=list()
rstr$x=seq(from.x,to.x,grid.spacing)
rstr$y=seq(from.y,to.y,grid.spacing)
rstr$z=grid.sum
str(rstr)
image(rstr,asp=1)
psat.rstr=raster(rstr, crs = proj);plot(psat.rstr)
writeRaster(psat.rstr, 'Norway.asc')
showWKT(proj4string(psat.rstr), file = "Norway.prj")

#Once in Arc, I draw polygons around the area with the highest probability and find the location
#At the center of the highest probability area.

#There are ways to overlay grid on satellite images in R too, so code may be expanded to include this later.



##############################################
#### Diagnostics
# just look at one record
i<- 5
poly.coords<- ellipses[[i]]@polygons[[1]]@Polygons[[1]]@coords
   poly.x<- poly.coords[,1]
   poly.y<- poly.coords[,2]
   pp<- point.in.polygon(point.x,point.y,poly.x,poly.y)
   pp[pp>0]<- argos[i,"score"]  
pp.m<- matrix(pp, nrow=length(seq(from.x,to.x,grid.spacing)),ncol=length(seq(from.y,to.y,grid.spacing)))
image(20* (1:nrow(grid.sum)), 20* (1:ncol(grid.sum)), pp.m, col = topo.colors(12))



# Look at all ellipses
plot(argos$X, argos$Y, main="test draw.ellipse", xlim=c(from.x,to.x),
   ylim=c(from.y,to.y))
draw.ellipse(argos$X, argos$Y, argos$minor_axis,argos$major_axis, angle=180-argos$orientation)
library(plotrix)


########################################
##### **** Ellipse function ***** ######
########################################

## modified code from Rolf Turner and Tony Fischbach, found on line

require(sp)
require(spatstat)
ellipseSP <- function (x0, y0, a, b, phi=pi/2, npts = 1000, D=T,...) {
## Returns SpatialPolygons ellipse
## phi defaults to 90 degrees
## default = degrees (D = T)
## whereby a lies in the N-S and
## b in the E-W directions
  if (D==T) phi = (phi/180)*pi
  theta <- seq(0, 2 * pi, length = npts)
  xh <- a * cos(theta)
  yh <- b * sin(theta)
  co <- cos(phi)
  si <- sin(phi)
  x<- x0 + co*xh - si*yh
  x[npts]<-x[1] ## force closure of polygon
  y<- y0 + si*xh + co*yh #This was x0 in Tony's code....
  y[npts]<-y[1] ## force closure of polygon
  return(SpatialPolygons(list(Polygons(list(with(list(x = x, y),Polygon(cbind(x,y)))),1))))
}

## Example plot
eSP <- ellipseSP(x0=0, y0=0, a=2, b=7,phi=360-93,D=T)
plot(eSP)
str(eSP)




