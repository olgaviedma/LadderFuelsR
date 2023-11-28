![](https://github.com/olgaviedma/LadderFuelsR/blob/master/readme/cover.png)<br/>

[![CRAN](https://www.r-pkg.org/badges/version/LadderFuelsR)](https://cran.r-project.org/package=LadderFuelsR)
![Github](https://img.shields.io/badge/Github-0.0.1-green.svg)
![licence](https://img.shields.io/badge/Licence-GPL--3-blue.svg) 
![Downloads](https://cranlogs.r-pkg.org/badges/grand-total/LadderFuelsR)
[![Build Status](https://travis-ci.com/olgaviedma/LadderFuelsR.svg?token=Jqizwyc6gBxNafNccTdU&branch=master)](https://travis-ci.com/olgaviedma/LadderFuelsR)


**LadderFuelsR: An R Package  for vertical fuel continuity analysis using LiDAR data.**

Authors: Olga Viedma  

The LadderFuelsR package is an automated tool for vertical fuel continuity analysis using LiDAR data that can be applied on multiple tree species and for large-scale studies. We include a suite of tools integrated into a standard workflow for deriving outputs meaningful to vertical fuel continuity analysis and canopy base height detection. The workflow consisted of: 1) calculating the vertical height profiles (VHP) of each segmented tree; 2) identifying gaps and fuel layers; 3) estimating the distance between fuel layers; and 4) retrieving the fuel layers base height (FBH) and depth. Additionally, other functions recalculate previous metrics after considering distances > 1 m and calculate the canopy base height (CBH) as the FBH located at the largest- and at the last-distance. Moreover, the package calculates: i) the percentage of LAD comprised in each fuel layer, ii) remove fuel layers with LAD percentage < 25 %, iii) recalculate the distances among the reminder ones, and iv) identify the CBH as the FBH with the highest LAD percentage. On the other hand, when the vertical height profiles (VHP) showed only one fuel layer, it identifies a possible CBH performing a segmented linear regression (breaking points) on the cumulative sum of LAD as a function of height. Finally, a collection of plotting functions is developed to represent: i) the VHP with the initial gaps and fuel layers; ii) the FBHs, depths and gaps with distances > 1 m and, iii) the FBHs and depths after applying the breaking point method over trees with only one fuel layer..

# Getting Started

## Installation
```r
#The CRAN version:
install.packages("LadderFuelsR")

# The development version:
#install.packages("remotes")
library(remotes)
install_github("https://github.com/olgaviedma/LadderFuelsR", dependencies = TRUE)

# loading rGEDI package
library(LadderFuelsR)

```    


## Rquired libraries

```{r pressure, echo=FALSE}
if (!require("pacman")) install.packages("pacman")
pacman::p_load(plyr,
               dplyr,
               tidyr,
               stringr,
               stringi,
               purrr,
               rlang,
               tidyverse,
               magrittr,
               sp,
               sf,
               raster,
               data.table,
               rgdal,
               viridis,
               lidR,
               leafR,
               segmented,
               lidRplugins,
               ggplot2,
               gt,
               gridExtra,
               patchwork,
               SSBtools,
               tibble,
               rgl,
               rglwidget,
               LadderFuelsR,
               magrittr,
               gdata)

```

## 1. Computing Canopy hieght model (CHM) using lidR package

```{r CHM pitfree 0.5 m, echo=TRUE, message=FALSE, warning=FALSE}

LIDAR_dir <- file.path(system.file("extdata", package = "LadderFuelsR"))
LIDAR_out<-file.path(system.file("extdata", package = "LadderFuelsR"))

las_list <-list.files (LIDAR_dir, pattern = "*.las", full.names = TRUE)
las_list1 <-list.files (LIDAR_dir, pattern = "*.las", full.names = FALSE)
short_name1<-stri_sub(las_list1, 1,-5)

lidar_maps<-lapply(las_list, function(X) lidR::readLAS(X, filter = "-drop_z_below 0"))

#plot(lidar_maps[[1]])

chm_pitfree_list<-list()
outputnames2<-NULL

for (i in seq_along(lidar_maps)) {
  
  lidar_maps1<-lidar_maps[[i]]
  
  chm_pitfree<- grid_canopy(lidar_maps1, res=0.5,pitfree( c(0,2,5,10,15,20,25,30,35,40), c(0,1.5), subcircle=0.15))
  chm_pitfree[chm_pitfree > 40] <- NA
  chm_pitfree[chm_pitfree < 0] <- 0
  chm_pitfree1 <- projectRaster(chm_pitfree, crs=26916)
  
  chm_pitfree_list[[i]]<-chm_pitfree1
  
  outputnames2<- paste0(LIDAR_out ,"/", short_name1[i] ,"_CHM05m.tif",collapse = NULL,sep="")
  writeRaster(chm_pitfree_list[[i]],outputnames2, overwrite=T)
  
}
col <- height.colors(25)
plot(chm_pitfree_list[[1]],col=col)

```

## 2.Detecting individual tree top from the lidar-derived CHM

```{r Tree tops detection, echo=TRUE, message=FALSE, warning=FALSE}

ttop_dir<- file.path(system.file("extdata", package = "LadderFuelsR"))

las_dir <- file.path(system.file("extdata", package = "LadderFuelsR"))
las_list <-list.files (las_dir, pattern = "*.las", full.names = TRUE)
las_list1 <-list.files (las_dir, pattern = "*.las", full.names = FALSE)
las_files<-lapply(las_list, function(X) lidR::readLAS(X, filter = "-drop_z_below 0"))


names1<-list()
filename_post2<-list()
chm_ttop_list2<-list()

for (j in seq_along(las_list)){
  
  short_name<-stri_sub(las_list1[j], 1, -5)
  names1<-rbind(names1,short_name)
  filename_post2 <- paste0(names1,"_multi2_tree_tops" )
  
  las_file <- lidR::readLAS(las_list[j],"-drop_z_below 0")
  
  # parameters
  ws= 2.5
  hmin = 2
  res=0.5
  ttops_multichm = find_trees(las_file, multichm(res = res, dist_2d = 2,ws= ws, layer_thickness = 0.3,dist_3d = 1, hmin = hmin))
  #ttops_multichm = find_trees(las_file, multichm(res = 2, dist_2d = 2,ws= 5, layer_thickness = 0.3,dist_3d = 1, hmin = 4))
  proj4string(ttops_multichm) <- CRS('+init=EPSG:26916')
  
  chm_ttop_list2[[j]]<-ttops_multichm
  
  writeOGR(chm_ttop_list2[[j]], dsn=ttop_dir ,layer = filename_post2[[j]],driver = 'ESRI Shapefile',
           overwrite_layer = T)
}

las1<-las_files[[1]]
chm1<-chm_pitfree_list[[1]]
ttops1<-chm_ttop_list2[[1]]

#plot(chm1, col = height.colors(50))
#plot(ttops1, add = TRUE, pch= 16, col = "black", main = "Tree tops over the CHM")

# Create an rgl point cloud
x<-add_treetops3d(plot(las_files[[1]], bg = "white", size = 4), ttops1)

# Customize the plot orientation
rgl.viewpoint(theta = 0, phi = 0, fov = 60, zoom = 0.75)

# Convert the rgl scene to an HTML widget
rglwidget(elementId = "x", width = 800, height = 600)

```

## 3. Individual tree crown deliniation (Silva et al. 2016) 

```{r Crowns Silva, echo=TRUE, message=FALSE, warning=FALSE}

crowns_dir <- file.path(system.file("extdata", package = "LadderFuelsR"))

ttop_dir <- file.path(system.file("extdata", package = "LadderFuelsR"))
ttop_list <-list.files (ttop_dir, pattern = "*_multi2_tree_tops.shp", full.names = TRUE)
ttops_multichm<-lapply(ttop_list, function (x) st_read(x))

chm_dir <- file.path(system.file("extdata", package = "LadderFuelsR"))
chm_list <-list.files (chm_dir, pattern = "_CHM05m.tif", full.names = TRUE)
las_dir <- file.path(system.file("extdata", package = "LadderFuelsR"))
las_list <-list.files (las_dir, pattern = "*.las", full.names = TRUE)


names1<-list()
filename_post1<-list()
crowns_silva_list<-list()

for (j in seq_along(las_list)){
  
  short_name<-gsub(".*/","",las_list[j])
  short_name1<-stri_sub(short_name, 1, -11)
  names1<-rbind(names1,short_name1)
  filename_post1 <- paste0(names1,"crown2m_silva.las" )

  chm05m <- raster(chm_list[j])
  ttops<-ttops_multichm[[j]]
  las_file <- lidR::readLAS(las_list[j])
  
  algo_silva1 <-silva2016(chm05m, ttops, max_cr_factor = 0.6, exclusion = 0.3, ID = "treeID")
  crowns_silva_las1 <-segment_trees(las_file, algo_silva1, attribute = "treeID", uniqueness = "incremental")
  crowns_silva_las2<-filter_poi(crowns_silva_las1, !is.na(treeID))
  
  crowns_silva_list[[j]]<-crowns_silva_las2

  lidR::writeLAS(crowns_silva_list[[j]], file = paste0(crowns_dir, "/", filename_post1[[j]], sep=""))
}

my_palette <- colorRampPalette(col)
x1<-plot(crowns_silva_list[[1]], color = "treeID", pal = my_palette, bg = "white")

# Customize the plot orientation
rgl.viewpoint(theta = 0, phi = 0, fov = 10, zoom = 0.75)

# Convert the rgl scene to an HTML widget
rglwidget(elementId = "x1", width = 800, height = 600)

```

## 4. Definiding function for computing crown-level metrics

```{r tree metrics function, echo=TRUE}

custom_crown_metrics <- function(z, i) { # user-defined function
  metrics <- list(
     dz = 1, 
     th = 1,
     z_max = max(z),   # max height
     z_min = min(z),   # min height
     z_mean = mean(z),   # mean height
     z_sd = sd(z), # vertical variability of points
     z_q1=quantile(z, probs = 0.01),
     z_q5=quantile(z, probs = 0.05),
    z_q25=quantile(z, probs = 0.25),
    z_q50=quantile(z, probs = 0.50),
    z_q75=quantile(z, probs = 0.75),
    z_q95=quantile(z, probs = 0.95),
     crr=(mean(z)-min(z))/(max(z)-min(z))
   )
   return(metrics) # output
}
ccm = ~custom_crown_metrics(z = Z, i = Intensity)

```

## 5.Computing crown level standard metrics within all trees detected

```{r tree and crown standard and own metrics, echo=TRUE, message=FALSE, warning=FALSE}

stats_out<-file.path(system.file("extdata", package = "LadderFuelsR"))
las_dir <- file.path(system.file("extdata", package = "LadderFuelsR"))
las_list <-list.files (las_dir, pattern = "*_crown2m_silva.las", full.names = TRUE)

names2<-list()
filename_post4<-list()
crowns_metrics_list<-list()

for (j in seq_along(las_list)){
  
  short_name<-stri_sub(las_list[j], 1, -5)
  short_name1<-gsub(".*/","",short_name)
  names2<-rbind(names2,short_name1)
  filename_post4 <- paste0(names2,"_metrics_convex.shp" )
  
  las_file1 <- lidR::readLAS(las_list[j])
  las_file2<-filter_poi(las_file1, !is.na(treeID))
  las_file3<-filter_poi(las_file2, Z >= 1)
  
  metrics1 = crown_metrics(las_file3,func = .stdtreemetrics, geom = "convex")
  crown_diam<-data.frame(sqrt(metrics1$convhull_area/ pi) * 2)
  names(crown_diam)<-"crown_diam"
  metrics2 = crown_metrics(las_file3,func = ccm, geom = "convex") #concave
  metrics_all <- dplyr::bind_cols(list(metrics1,crown_diam,metrics2))
  metrics_all1 <- metrics_all[,c(1:4,6,10:21)]
  names(metrics_all1)<-c("treeID", "Z", "npoints", "convhull_area", "crown_diam", "z_max", "z_min", "z_mean","z_sd", "z_q1","z_q5", "z_q25","z_q50","z_q75", "z_q95", "crr", "geometry" )
  
  metrics_all2<- st_as_sf(metrics_all1)
  
  crowns_metrics_list[[j]]<-metrics_all2
  
  st_write(crowns_metrics_list[[j]], paste0(dsn=stats_out, "/", filename_post4[[j]],sep=""), driver = 'ESRI Shapefile',append=F)
  
}

ttops1<-st_as_sf(chm_ttop_list2[[1]])
crowns1<-st_as_sf(crowns_metrics_list[[1]])
ttops_within_crowns <- st_intersection(ttops1, crowns1)

# Set the size of the plotting device
par(mfrow = c(1, 1), mar = c(1, 1, 1, 1), pin = c(5, 4)) 
plot(st_geometry(crowns1), pch = 16, col = "green")
plot(ttops_within_crowns, add = TRUE, pch= 16, col = "darkblue", main = "Tree tops over the crowns")


```
![](https://github.com/carlos-alberto-silva/rGEDI/blob/master/readme/fig3.png)

### 6. Checking for no overlapping crown poygons

```{r no overlapping crown polygons, echo=TRUE, message=FALSE, warning=FALSE}

dir_out<-file.path(system.file("extdata", package = "LadderFuelsR"))
tree_dir <- file.path(system.file("extdata", package = "LadderFuelsR"))
tree_list <-list.files (tree_dir, pattern = "*_crown2m_silva_metrics_convex.shp", full.names = TRUE)
tree_files <-lapply(tree_list, function (X) st_read(X))


polys1<-tree_files[[1]]

for(i in 1:nrow(polys1)) {
  other_polys <- st_union(polys1[-i, ]$geometry)
  polys1[i, ]$geometry <- st_difference(polys1[i, ]$geometry, other_polys)
}
st_write(polys1,paste0(tree_dir, "/", "crown2m_silva_convex_no_overlap.shp"), append=FALSE)

# Set the size of the plotting device
#par(mfrow = c(1, 1), mar = c(1, 1, 1, 1), pin = c(5, 4)) 
#plot(st_geometry(polys1), pch = 16, col = "orange")
#plot(ttops_within_crowns, add = TRUE, pch= 16, col = "blue", main = "Tree tops over the crowns")


```
![](https://github.com/carlos-alberto-silva/rGEDI/blob/master/readme/fig3.png)


## 7.Extracting individual trees within the crown polygons

```{r cropLAS files with no overlapping crowns, echo=TRUE, message=FALSE, warning=FALSE}

dir_out<-file.path(system.file("extdata", package = "LadderFuelsR"))

las_dir <- file.path(system.file("extdata", package = "LadderFuelsR"))
las_list <-list.files (las_dir, pattern = "*_crown2m_silva.las", full.names = TRUE)
las_files <-lapply(las_list, function (X) lidR::readLAS(X))

tree_dir <- file.path(system.file("extdata", package = "LadderFuelsR"))
tree_list <-list.files (tree_dir, pattern = "*_no_overlap.shp", full.names = TRUE)
tree_files <-lapply(tree_list, function (X) st_read(X))

###################################

names1 <- list()
filename_post3 <- list()

for (j in seq_along(las_list)) {
  short_name <- stri_sub(las_list[j])
  short_name1 <- gsub(".*/", "", short_name)
  short_name2 <- stri_sub(short_name1, 1, 11)
  
  print(short_name2)
  
  if (is.na(short_name2)) {
    # Handle NA values or skip them
    next
  }
  
  names1[[j]] <- short_name2
  filename_post3[[j]] <- paste0(short_name2, "_CROWN.las")
  
  las1 <- las_files[[j]] 
  trees1 <- tree_files[[j]] 
  
  trees_ID <- trees1 %>% dplyr::select(treeID)
  n <- nrow(trees_ID)
  
  crown_cort <- vector("list", length=n)
  
  for (i in 1:n) {
    kk <- trees_ID[i,]
    crown_cort[[i]] = clip_roi(las1, kk)
  
    lidR::writeLAS(crown_cort[[i]], paste0(dir_out,"/", i, "_", filename_post3[j]))
  }
}

my_palette <- colorRampPalette(col)
x2<-plot(crown_cort[[1]], color = "Z", pal = my_palette, bg = "black", size = 2.5)

# Customize the plot orientation
rgl.viewpoint(theta = 0, phi = 0, fov = 60, zoom = 0.75)

# Convert the rgl scene to an HTML widget
rglwidget(elementId = "x2", width = 400, height = 600)

```
![](https://github.com/carlos-alberto-silva/rGEDI/blob/master/readme/fig3.png)


## 8. Computing LAI-LAD derived metrics within all trees extracted

```{r LAI and LAD tree metrics, echo=TRUE, message=FALSE, warning=FALSE}

stats_out<-file.path(system.file("extdata", package = "LadderFuelsR"))

las_dir1 <- file.path(system.file("extdata", package = "leafR"))
las_list1 <- list.files(las_dir1, pattern = "*_CROWN.las", full.names = TRUE, ignore.case = TRUE)

# create a vector to hold the file names of .las files with more than 3 points
files_with_more_than_3_points <- c()

# loop through each file
for (file in las_list1) {
  las_data <- lidR::readLAS(file)
  las_data1<-filter_poi(las_data, Z >= 1)
  
  # skip to next file if there was a problem reading
  if (is.null(las_data1)) next
  
  # check if it contains more than three points
  if (las_data1@header$`Number of point records` > 3) {
    files_with_more_than_3_points <- c(files_with_more_than_3_points, file)
  }
}

# Creates a data frame of the 3D voxels information (xyz) with Leaf Area Density values
short_name1<-NULL
profile_list<-NULL
lidar_lai_list<-NULL
understory_lai_list<-NULL
LAHV_metric_list<-NULL

for (j in seq_along(files_with_more_than_3_points)){
  
  short_name<-stri_sub(files_with_more_than_3_points[j], 1, -5)
  short_name1<-gsub(".*/","",short_name)
  
  normlas_file<-files_with_more_than_3_points[[j]]

  VOXELS_LAD = lad.voxels(normlas_file, grain.size = 2)
  
  lad_profile = lad.profile(VOXELS_LAD, relative = F)
  lai_tot = lai(lad_profile)
  understory_lai <- lai(lad_profile, min = 0.3, max = 2.5)
  LAHV_metric<- LAHV(lad_profile, LAI.weighting = FALSE, height.weighting = FALSE) # Leaf Area height volume (wood volume, AGB)
  
  lad_profile1 = data.frame(lad_profile, treeID = short_name1)
  lai_tot1 = data.frame(lai_tot, treeID = short_name1)
  understory_lai1 = data.frame(understory_lai, treeID = short_name1)
  LAHV_metric1 = data.frame(LAHV_metric, treeID = short_name1)
  
  profile_list<-rbind(profile_list, lad_profile1)
  lidar_lai_list<-rbind(lidar_lai_list,lai_tot1)
  understory_lai_list <-rbind(understory_lai_list,understory_lai1)
  LAHV_metric_list<-rbind(LAHV_metric_list,LAHV_metric1)
}

write.table(profile_list,file =paste(stats_out, "/","alltrees_LAD_profile_voxels2m",".txt", sep=""), sep="\t", row.names = FALSE)
write.table(lidar_lai_list,file =paste(stats_out, "/","alltrees_LAI_profile_voxels2m_",".txt", sep=""), sep="\t", row.names = FALSE)
write.table(understory_lai_list,file =paste(stats_out, "/","alltrees_understory_profile_voxels2m",".txt", sep=""), sep="\t", row.names = FALSE)
write.table(LAHV_metric_list,file =paste(stats_out, "/","alltrees_LAHV_profile_voxels2m",".txt", sep=""), sep="\t", row.names = FALSE)
#head(profile_list)
```
##           shot_number latitude_bin0 latitude_lastbin longitude_bin0 longitude_lastbin elevation_bin0
##  1: 19640002800109382     -13.75903        -13.75901      -44.17219         -44.17219       784.8348
##  2: 19640003000109383     -13.75862        -13.75859      -44.17188         -44.17188       799.0491
##  3: 19640003200109384     -13.75821        -13.75818      -44.17156         -44.17156       814.4647
##  4: 19640003400109385     -13.75780        -13.75777      -44.17124         -44.17124       820.1437
##  5: 19640003600109386     -13.75738        -13.75736      -44.17093         -44.17093       821.7012
##  6: 19640003800109387     -13.75697        -13.75695      -44.17061         -44.17061       823.2526


## 9. Depurating tree lab profiles (\> 3 height values)

```{r depurating LAD databases, echo=TRUE, message=FALSE, warning=FALSE}

LAD_path<- file.path(system.file("extdata", package = "LadderFuelsR"),"alltrees_LAD_profile_voxels2m.txt")
stats_tot<- read.table(LAD_path,sep="\t", header=T)
cols <- c('treeID')
stats_tot[cols] <- lapply(stats_tot[cols], function (x) as.factor(x))
stats_tot$lad<-round(stats_tot$lad,digits = 4)

cases <- data.frame(table(stats_tot$treeID))
cases1 <-cases[cases$Freq > 3, ]
names(cases1)<-c("treeID", "Freq")

stats_tot1 <- stats_tot[stats_tot$treeID %in% cases1$treeID, ]
stats_tot2 <- data.frame(stats_tot1 %>% replace(is.na(.), 0.01))

output_file <- file.path(system.file("extdata", package = "LadderFuelsR"), "LAD_profiles.txt")
write.table(stats_tot2, file = output_file, sep = "\t", row.names = FALSE)

```

## 10. Computing gaps and fuel layers base height (FBH)

```{r Gaps and Fuel layers Base Height (fbh), echo=TRUE, message=FALSE, warning=FALSE}

# LAD profiles derived from normalized ALS data after applying [lad.profile()] function
data_path <- file.path(system.file("extdata", package = "LadderFuelsR"), "LAD_profiles.txt")
LAD_profiles <- read.table(data_path, sep = "\t", header = TRUE)
LAD_profiles$treeID <- factor(LAD_profiles$treeID)

trees_name1 <- as.character(LAD_profiles$treeID)
trees_name2 <- factor(unique(trees_name1))

metrics_precentile_list<-list()
for (i in levels(trees_name2)) {
  tree2 <- LAD_profiles |> dplyr::filter(treeID == i)
  metrics_precentil <- get_gaps_fbhs(tree2)
  metrics_precentile_list[[i]] <- metrics_precentil
}

metrics_all_percentil <- dplyr::bind_rows(metrics_precentile_list)
metrics_all_percentil$treeID <- factor(metrics_all_percentil$treeID)

# Remove the row with all NA values from the original data frame
# First remove "treeID" and "treeID1" columns
tree_metrics_no_treeID <- metrics_all_percentil[, -which(names(metrics_all_percentil) == c("treeID","treeID1"))]

# Check if any row has all NA values
rows_with_all_NA_or_zero <- apply(tree_metrics_no_treeID, 1, function(row) all(is.na(row) | row == 0))

# Get the row index with all NA values
row_index <- which(rows_with_all_NA_or_zero)

# Remove the row with all NA values from the original data frame
if (length(row_index) > 0) {
  tree_metrics_filtered <- metrics_all_percentil[-row_index, ]
} else {
  tree_metrics_filtered <- metrics_all_percentil
}
write.table(tree_metrics_filtered, file= file.path(system.file("extdata", package = "LadderFuelsR"), "1_gaps_fbhs_metrics.txt"),
sep = "\t",row.names = FALSE)

#head(tree_metrics_filtered)

```
##           shot_number latitude_bin0 latitude_lastbin longitude_bin0 longitude_lastbin elevation_bin0
##  1: 19640002800109382     -13.75903        -13.75901      -44.17219         -44.17219       784.8348
##  2: 19640003000109383     -13.75862        -13.75859      -44.17188         -44.17188       799.0491
##  3: 19640003200109384     -13.75821        -13.75818      -44.17156         -44.17156       814.4647
##  4: 19640003400109385     -13.75780        -13.75777      -44.17124         -44.17124       820.1437
##  5: 19640003600109386     -13.75738        -13.75736      -44.17093         -44.17093       821.7012
##  6: 19640003800109387     -13.75697        -13.75695      -44.17061         -44.17061       823.2526


## 11. Computing distance between fuel layers

```{r Distances (and their heights) between fuel layers, echo=TRUE, message=FALSE, warning=FALSE}

# Tree metrics derived from get_gaps_fbhs() function
gaps_fbhs_metrics_path <- system.file("extdata", "1_gaps_fbhs_metrics.txt", package = "LadderFuelsR")
gaps_fbhs_metrics <- read.table(gaps_fbhs_metrics_path, sep="\t", header=TRUE)

gaps_fbhs_metrics$treeID <- factor(gaps_fbhs_metrics$treeID)
trees_name1 <- as.character(gaps_fbhs_metrics$treeID)
trees_name2 <- factor(unique(trees_name1))

metrics_distance_list <- list()

for (i in levels(trees_name2)) {

  # Filter data for each tree
  tree2 <- gaps_fbhs_metrics |> dplyr::filter(treeID == i)

  # Get distance metrics for each tree
  metrics_distance <- get_distance(tree2)
  metrics_distance_list[[i]] <- metrics_distance
}

# Combine the individual data frames
metrics_all_distance <- dplyr::bind_rows(metrics_distance_list)

distance_path <- file.path(system.file("extdata", package = "LadderFuelsR"), "2_distance_metrics.txt")
write.table(metrics_all_distance, file=distance_path, sep="\t", row.names=FALSE)

metrics_all_distance <- metrics_all_distance[, order(names(metrics_all_distance))]
#head(metrics_all_distance)

```
##           shot_number latitude_bin0 latitude_lastbin longitude_bin0 longitude_lastbin elevation_bin0
##  1: 19640002800109382     -13.75903        -13.75901      -44.17219         -44.17219       784.8348
##  2: 19640003000109383     -13.75862        -13.75859      -44.17188         -44.17188       799.0491
##  3: 19640003200109384     -13.75821        -13.75818      -44.17156         -44.17156       814.4647
##  4: 19640003400109385     -13.75780        -13.75777      -44.17124         -44.17124       820.1437
##  5: 19640003600109386     -13.75738        -13.75736      -44.17093         -44.17093       821.7012
##  6: 19640003800109387     -13.75697        -13.75695      -44.17061         -44.17061       823.2526


## 12. Computing fuel layer depth

```{r Distane between fuel layers, echo=TRUE, message=FALSE, warning=FALSE}

library(dplyr)
library(magrittr)

# LAD profiles derived from normalized ALS data after applying [lad.profile()] function
data_path <- file.path(system.file("extdata", package = "LadderFuelsR"), "LAD_profiles.txt")
LAD_profiles <- read.table(data_path, sep = "\t", header = TRUE)
LAD_profiles$treeID <- factor(LAD_profiles$treeID)

# Tree metrics derived from get_distance() function
distance_path <- file.path(system.file("extdata", package = "LadderFuelsR"), "2_distance_metrics.txt")
distance_metrics <- read.table(distance_path, sep = "\t", header = TRUE)
distance_metrics$treeID <- factor(distance_metrics$treeID)

metrics_depth_list <- list()

for (i in levels(LAD_profiles$treeID)){

  tree1 <- LAD_profiles |> dplyr::filter(treeID == i)
  tree2 <- distance_metrics |> dplyr::filter(treeID == i)

  # Get depths for each tree
  metrics_depth <- get_depths(tree1, tree2)
  metrics_depth_list[[i]] <- metrics_depth
}

# Combine the individual data frames
metrics_all_depth <- dplyr::bind_rows(metrics_depth_list)

depth_path <- file.path(system.file("extdata", package = "LadderFuelsR"), "3_depth_metrics.txt")
write.table(metrics_all_depth, file = depth_path, sep = "\t", row.names = FALSE)

metrics_all_depth <- metrics_all_depth[, order(names(metrics_all_depth))]
#head(metrics_all_depth)

```
##           shot_number latitude_bin0 latitude_lastbin longitude_bin0 longitude_lastbin elevation_bin0
##  1: 19640002800109382     -13.75903        -13.75901      -44.17219         -44.17219       784.8348
##  2: 19640003000109383     -13.75862        -13.75859      -44.17188         -44.17188       799.0491
##  3: 19640003200109384     -13.75821        -13.75818      -44.17156         -44.17156       814.4647
##  4: 19640003400109385     -13.75780        -13.75777      -44.17124         -44.17124       820.1437
##  5: 19640003600109386     -13.75738        -13.75736      -44.17093         -44.17093       821.7012
##  6: 19640003800109387     -13.75697        -13.75695      -44.17061         -44.17061       823.2526


## 13. Plotting gaps and fuel layer base height (FBH)

```{r Plots Gaps and Fuel layers Base Height (fbh), echo=TRUE, message=FALSE, warning=FALSE}

library(LadderFuelsR)
library(ggplot2)
library(lattice)

# LAD profiles derived from normalized ALS data after applying [lad.profile()] function
data_path <- file.path(system.file("extdata", package = "LadderFuelsR"), "LAD_profiles.txt")
LAD_profiles <- read.table(data_path, sep = "\t", header = TRUE)
LAD_profiles$treeID <- factor(LAD_profiles$treeID)

# Tree metrics derived from get_depths() function
depth_path <- file.path(system.file("extdata", package = "LadderFuelsR"), "3_depth_metrics.txt")
depth_metrics <- read.table(depth_path, sep = "\t", header = TRUE)
depth_metrics$treeID <- factor(depth_metrics$treeID)

# Generate plots for gaps and fbhs
plots_gaps_fbhs <- get_plots_gap_fbh(LAD_profiles, depth_metrics)

# Save plots for each tree
for (name in names(plots_gaps_fbhs)) {
  #print(plots_gaps_fbhs[[name]])  # print the plot
  ggsave(file.path(system.file("extdata", package = "LadderFuelsR"),  paste0( name, "_gap_fbh", ".tiff")), plot = plots_gaps_fbhs[[name]])
}

par(mfrow = c(1, 1))
# Plot the first panel
plot(plots_gaps_fbhs[[1]],width = 5, height = 5)

```

## 14.Computig fuel layer base height (FBH) after reomving distances =1

```{r Fuels base height after removing distances equal 1 m, echo=TRUE, message=FALSE, warning=FALSE}

library(SSBtools)
library(dplyr)
library(magrittr)

# Tree metrics derived from get_depths() function
depth_path <- file.path(system.file("extdata", package = "LadderFuelsR"), "3_depth_metrics.txt")
depth_metrics <- read.table(depth_path, sep = "\t", header = TRUE)
depth_metrics$treeID <- factor(depth_metrics$treeID)

trees_name1 <- as.character(depth_metrics$treeID)
trees_name2 <- factor(unique(trees_name1))

fbh_corr_list <- list()

for (i in levels(trees_name2)){

  # Filter data for each tree
  tree3 <- depth_metrics |> dplyr::filter(treeID == i)

  # Get real fbh for each tree
  fbh_corr <- get_real_fbh(tree3)

  # Store fbh values in a list
  fbh_corr_list[[i]] <- fbh_corr
}

# Combine fbh values for all trees
fbh_corr_all <- dplyr::bind_rows(fbh_corr_list)
fbh_corr_all$treeID <- factor(fbh_corr_all$treeID)

# Reorder columns
# Get original column names
original_column_names <- colnames(fbh_corr_all)

# Specify prefixes
prefixes <- c("treeID", "Hdist", "Hcbh", "Hdepth", "dist", "depth", "max_height")

# Initialize vector to store new order
new_order <- c()

# Loop over prefixes
for (prefix in prefixes) {
  # Find column names matching the current prefix
  matching_columns <- grep(paste0("^", prefix), original_column_names, value = TRUE)
  # Append to the new order
  new_order <- c(new_order, matching_columns)
}

# Reorder values
fbh_corr_all <- fbh_corr_all[, new_order]

# Save the reordered data
fbh_corr_path <- file.path(system.file("extdata", package = "LadderFuelsR"), "4_fbh_metrics_corr.txt")
write.table(fbh_corr_all, file = fbh_corr_path, sep = "\t", row.names = FALSE)
#head(fbh_corr_all)


```
##           shot_number latitude_bin0 latitude_lastbin longitude_bin0 longitude_lastbin elevation_bin0
##  1: 19640002800109382     -13.75903        -13.75901      -44.17219         -44.17219       784.8348
##  2: 19640003000109383     -13.75862        -13.75859      -44.17188         -44.17188       799.0491
##  3: 19640003200109384     -13.75821        -13.75818      -44.17156         -44.17156       814.4647
##  4: 19640003400109385     -13.75780        -13.75777      -44.17124         -44.17124       820.1437
##  5: 19640003600109386     -13.75738        -13.75736      -44.17093         -44.17093       821.7012
##  6: 19640003800109387     -13.75697        -13.75695      -44.17061         -44.17061       823.2526



## 15.Computing fuel lauer depth after removing distance =1

```{r Fuel layers depth after removinG distances equal 1 m, echo=TRUE, message=FALSE, warning=FALSE}

library(dplyr)
library(magrittr)
library(tidyr)
# Tree metrics derived from get_real_fbh() function
fbhcor_path <- file.path(system.file("extdata", package = "LadderFuelsR"), "4_fbh_metrics_corr.txt")
fbh_metrics_corr <- read.table(fbhcor_path, sep = "\t", header = TRUE)
fbh_metrics_corr$treeID <- factor(fbh_metrics_corr$treeID)

trees_name1 <- as.character(fbh_metrics_corr$treeID)
trees_name2 <- factor(unique(trees_name1))

depth_metrics_corr <- lapply(levels(trees_name2), function(i) {
  # Filter data for each tree
  tree2 <- fbh_metrics_corr |> dplyr::filter(treeID == i)
  # Get real depths for each tree
  get_real_depths(tree2)
})

# Combine depth values for all trees
depth_corr_all <- dplyr::bind_rows(depth_metrics_corr)
depth_corr_path <- file.path(system.file("extdata", package = "LadderFuelsR"), "5_tree_depth_metrics_corr.txt")
write.table(depth_corr_all, file = depth_corr_path, sep = "\t", row.names = FALSE)

#head(depth_corr_all)

```
##           shot_number latitude_bin0 latitude_lastbin longitude_bin0 longitude_lastbin elevation_bin0
##  1: 19640002800109382     -13.75903        -13.75901      -44.17219         -44.17219       784.8348
##  2: 19640003000109383     -13.75862        -13.75859      -44.17188         -44.17188       799.0491
##  3: 19640003200109384     -13.75821        -13.75818      -44.17156         -44.17156       814.4647
##  4: 19640003400109385     -13.75780        -13.75777      -44.17124         -44.17124       820.1437
##  5: 19640003600109386     -13.75738        -13.75736      -44.17093         -44.17093       821.7012
##  6: 19640003800109387     -13.75697        -13.75695      -44.17061         -44.17061       823.2526



## 16.Computing fuel layer distance (\> 1 M) and CBH based on maximum and last distance

```{r Fuel layers distances after removing distances equal 1 m, echo=TRUE, message=FALSE, warning=FALSE}

library(dplyr)
library(magrittr)
library(stringr)

# Tree metrics derived from get_real_depths() function
depth_corr_path <- file.path(system.file("extdata", package = "LadderFuelsR"), "5_tree_depth_metrics_corr.txt")
depth_metrics_corr <- read.table(depth_corr_path, sep = "\t", header = TRUE)
depth_metrics_corr$treeID <- factor(depth_metrics_corr$treeID)

trees_name1 <- as.character(depth_metrics_corr$treeID)
trees_name2 <- factor(unique(trees_name1))

distance_metrics_corr <- lapply(levels(trees_name2), function(i) {
  # Filter data for each tree
  tree2 <- depth_metrics_corr |> dplyr::filter(treeID == i)
  # Get effective gap for each tree
  get_effective_gap(tree2)
})

# Combine the individual data frames
distances_corr_all <- dplyr::bind_rows(distance_metrics_corr)

# =======================================================================#
# REORDER COLUMNS:
# =======================================================================#
# Get original column names
original_column_names <- colnames(distances_corr_all)

# Specify prefixes
prefixes <- c("treeID", "Hcbh", "dptf", "Hdptf", "effdist", "dist", "Hdist", "max_Hcbh", "max_dptf", "max_Hdptf", "last_Hcbh", "last_dptf", "last_Hdptf", "max_height")

# Initialize vector to store new order
new_order <- c()

# Loop over prefixes
for (prefix in prefixes) {
  # Find column names matching the current prefix
  matching_columns <- grep(paste0("^", prefix), original_column_names, value = TRUE)

  # Extract numeric suffixes and order the columns based on these suffixes
  numeric_suffixes <- as.numeric(gsub(paste0("^", prefix), "", matching_columns))
  matching_columns <- matching_columns[order(numeric_suffixes)]

  # Append to new order
  new_order <- c(new_order, matching_columns)
}

# Reorder values
distances_corr_all1 <- distances_corr_all[, new_order]
distance_corr_path <- file.path(system.file("extdata", package = "LadderFuelsR"), "6_tree_distances_metrics_corr.txt")
write.table(distances_corr_all1, file = distance_corr_path, sep = "\t", row.names = FALSE)

#head(distances_corr_all1)

```
##           shot_number latitude_bin0 latitude_lastbin longitude_bin0 longitude_lastbin elevation_bin0
##  1: 19640002800109382     -13.75903        -13.75901      -44.17219         -44.17219       784.8348
##  2: 19640003000109383     -13.75862        -13.75859      -44.17188         -44.17188       799.0491
##  3: 19640003200109384     -13.75821        -13.75818      -44.17156         -44.17156       814.4647
##  4: 19640003400109385     -13.75780        -13.75777      -44.17124         -44.17124       820.1437
##  5: 19640003600109386     -13.75738        -13.75736      -44.17093         -44.17093       821.7012
##  6: 19640003800109387     -13.75697        -13.75695      -44.17061         -44.17061       823.2526


## 17.Computing fuel layer LAD in percentage (\> 25 %) and CBH based on maximum LAD percentage

```{r Fuels LAD percentage and canopy base height (CBH) based on maximum LAD percentage (distances greater than 1 m), echo=TRUE, message=FALSE, warning=FALSE}

library(dplyr)
library(magrittr)

# LAD profiles derived from normalized ALS data after applying [lad.profile()] function
data_path <- system.file("extdata", "LAD_profiles.txt", package = "LadderFuelsR")
LAD_profiles <- read.table(data_path, sep = "\t", header = TRUE)
LAD_profiles$treeID <- factor(LAD_profiles$treeID)

# Tree metrics derived from get_effective_gap() function
distcor_path <- file.path(system.file("extdata", package = "LadderFuelsR"), "6_tree_distances_metrics_corr.txt")
distances_metrics_corr <- read.table(distcor_path, sep = "\t", header = TRUE)
distances_metrics_corr$treeID <- factor(distances_metrics_corr$treeID)

trees_name1 <- as.character(distances_metrics_corr$treeID)
trees_name2 <- factor(unique(trees_name1))

LAD_metrics1 <- list()
LAD_metrics2 <- list()

for (i in levels(trees_name2)) {
  # Filter data for each tree
  tree1 <- LAD_profiles |> dplyr::filter(treeID == i)
  tree2 <- distances_metrics_corr |> dplyr::filter(treeID == i)

  # Get LAD metrics for each tree
  LAD_metrics <- get_layers_lad(tree1, tree2)
  LAD_metrics1[[i]] <- LAD_metrics$df1
  LAD_metrics2[[i]] <- LAD_metrics$df2
}

LAD_metrics_all1 <- dplyr::bind_rows(LAD_metrics1)
LAD_metrics_all2 <- dplyr::bind_rows(LAD_metrics2)

# List of data frames
LAD_metrics_list <- list(LAD_metrics_all1, LAD_metrics_all2)

# Initialize an empty list to store reordered data frames
reordered_list <- list()

# Specify prefixes (adjust accordingly)
prefixes <- c("treeID", "Hdist", "Hcbh", "effdist", "dptf", "Hdptf", "max", "last")

# Loop over each data frame
for (i in seq_along(LAD_metrics_list)) {

  LAD_metrics_all <- LAD_metrics_list[[i]]

  # Get original column names
  original_column_names <- colnames(LAD_metrics_all)

  # Initialize vector to store new order
  new_order <- c()

  # Loop over prefixes
  for (prefix in prefixes) {
    # Find column names matching the current prefix
    matching_columns <- grep(paste0("^", prefix), original_column_names, value = TRUE)

    # Extract numeric suffixes and order the columns based on these suffixes
    numeric_suffixes <- as.numeric(gsub(paste0("^", prefix), "", matching_columns))
    
       # Order the columns based on numeric suffixes
    matching_columns <- matching_columns[order(numeric_suffixes)]

    # Append to new order
    new_order <- c(new_order, matching_columns)
  }
  # Reorder columns
  LAD_metrics_all <- LAD_metrics_all[, new_order]
  # Store the reordered data frame in the list
  reordered_list[[i]] <- LAD_metrics_all
 }
  # Write the reordered data frame to a file
 fuels_lad_all_path <- file.path(system.file("extdata", package = "LadderFuelsR"), "7_fuels_lad_all.txt")
 write.table(reordered_list[[1]], file = fuels_lad_all_path, sep = "\t", row.names = FALSE)
 fuels_lad_gt25_path <- file.path(system.file("extdata", package = "LadderFuelsR"), "7_fuels_lad_gt25perc.txt")
 write.table(reordered_list[[2]], file = fuels_lad_gt25_path, sep = "\t", row.names = FALSE)
 
 #head(reordered_list[[2]])
 
```
##           shot_number latitude_bin0 latitude_lastbin longitude_bin0 longitude_lastbin elevation_bin0
##  1: 19640002800109382     -13.75903        -13.75901      -44.17219         -44.17219       784.8348
##  2: 19640003000109383     -13.75862        -13.75859      -44.17188         -44.17188       799.0491
##  3: 19640003200109384     -13.75821        -13.75818      -44.17156         -44.17156       814.4647
##  4: 19640003400109385     -13.75780        -13.75777      -44.17124         -44.17124       820.1437
##  5: 19640003600109386     -13.75738        -13.75736      -44.17093         -44.17093       821.7012
##  6: 19640003800109387     -13.75697        -13.75695      -44.17061         -44.17061       823.2526




## 18. Plotting fuel layers with LAD \> 25 % and CBH based on maximum LAD percentage

```{r Plots of fuel layers with LAD percentage greater than 25 and the canopy base height (CBH) based on the maximum LAD percentage, echo=TRUE, message=FALSE, warning=FALSE}

library(ggplot2)

# LAD profiles derived from normalized ALS data after applying [lad.profile()] function
data_path <- file.path(system.file("extdata", package = "LadderFuelsR"), "LAD_profiles.txt")
LAD_profiles <- read.table(data_path, sep = "\t", header = TRUE)
LAD_profiles$treeID <- factor(LAD_profiles$treeID)

# Tree metrics derived from get_layers_lad() function
LAD_gt25p_path <- file.path(system.file("extdata", package = "LadderFuelsR"), "7_fuels_lad_gt25perc.txt")
fuels_LAD_metrics <- read.table(LAD_gt25p_path, sep = "\t", header = TRUE)
fuels_LAD_metrics$treeID <- factor(fuels_LAD_metrics$treeID)

trees_name1 <- as.character(fuels_LAD_metrics$treeID)
trees_name2 <- factor(unique(trees_name1))

# Generate plots for fuels LAD metrics
plots_trees_LAD <- get_plots_cbh_LAD(LAD_profiles, fuels_LAD_metrics)

# Save plots for each tree
for (name in names(plots_trees_LAD)) {
  plots <- plots_trees_LAD[[name]]

  if (!is.null(plots)) {
    #print(paste("Saving plot for tree:", name))
    ggsave(file.path(system.file("extdata", package = "LadderFuelsR"),  paste0( name, "_LAD_25perc", ".tiff")), plot = plots,
           width = 6, height = 6, units = "in")
  }}

par(mfrow = c(1, 1))
# Plot the first panel
plot(plots_trees_LAD[[1]],width = 6, height = 6)

```
![](https://github.com/carlos-alberto-silva/rGEDI/blob/master/readme/fig3.png)


## 19. Compute CBH based on tree breaking point method and LAD percentage

```{r CBH and the LAD percentage below and above the CBH using the breaking point method, echo=TRUE, message=FALSE, warning=FALSE}

library(dplyr)
library(magrittr)

# LAD profiles derived from normalized ALS data after applying [lad.profile()] function
data_path <- file.path(system.file("extdata", package = "LadderFuelsR"), "LAD_profiles.txt")
LAD_profiles <- read.table(data_path, sep = "\t", header = TRUE)
LAD_profiles$treeID <- factor(LAD_profiles$treeID)

# Tree metrics derived from get_effective_gap() function
distcor_path <- file.path(system.file("extdata", package = "LadderFuelsR"), "6_tree_distances_metrics_corr.txt")
distances_metrics_corr <- read.table(distcor_path, sep = "\t", header = TRUE)
distances_metrics_corr$treeID <- factor(distances_metrics_corr$treeID)

trees_name1 <- as.character(distances_metrics_corr$treeID)
trees_name2 <- factor(unique(trees_name1))

cum_LAD_metrics_list <- list()

for (i in levels(trees_name2)) {
  # Filter data for each tree
  tree1 <- LAD_profiles |> dplyr::filter(treeID == i)
  tree2 <- distances_metrics_corr |> dplyr::filter(treeID == i)

  # Get cumulative LAD metrics for each tree
  cum_LAD_metrics <- get_cum_break(tree1, tree2)
  cum_LAD_metrics_list[[i]] <- cum_LAD_metrics
}

# Combine the individual data frames
cum_LAD_metrics_all <- dplyr::bind_rows(cum_LAD_metrics_list) 

# =======================================================================#
# REORDER COLUMNS
# =======================================================================#

# Get original column names
original_column_names <- colnames(cum_LAD_metrics_all)

# Specify prefixes (adjust accordingly)
prefixes <- c("treeID", "Hcbh", "below", "above", "depth", "Hdepth", "dptf", "Hdptf", "Hdist", "effdist", "max", "last", "cumlad")

# Initialize vector to store new order
new_order <- c()

# Loop over prefixes
for (prefix in prefixes) {
  # Find column names matching the current prefix
  matching_columns <- grep(paste0("^", prefix), original_column_names, value = TRUE)

  # Extract numeric suffixes and order the columns based on these suffixes
  numeric_suffixes <- as.numeric(gsub(paste0("^", prefix), "", matching_columns))
  matching_columns <- matching_columns[order(numeric_suffixes)]

  # Append to new order
  new_order <- c(new_order, matching_columns)
}

# Reorder columns
cum_LAD_metrics_order <- cum_LAD_metrics_all[, new_order]
cum_LAD_metrics_path <- file.path(system.file("extdata", package = "LadderFuelsR"), "8_cbh_breaking_point_lad.txt")
write.table(cum_LAD_metrics_order, file = cum_LAD_metrics_path, sep = "\t", row.names = FALSE)

#head(cum_LAD_metrics_order)

```
##           shot_number latitude_bin0 latitude_lastbin longitude_bin0 longitude_lastbin elevation_bin0
##  1: 19640002800109382     -13.75903        -13.75901      -44.17219         -44.17219       784.8348
##  2: 19640003000109383     -13.75862        -13.75859      -44.17188         -44.17188       799.0491
##  3: 19640003200109384     -13.75821        -13.75818      -44.17156         -44.17156       814.4647
##  4: 19640003400109385     -13.75780        -13.75777      -44.17124         -44.17124       820.1437
##  5: 19640003600109386     -13.75738        -13.75736      -44.17093         -44.17093       821.7012
##  6: 19640003800109387     -13.75697        -13.75695      -44.17061         -44.17061       823.2526


## 20. Plotting CBH based on the breaking point method and LAD percentage

```{r Plots of the CBH and the LAD percentage below and above the CBH using the breaking point method, echo=TRUE, message=FALSE, warning=FALSE}

library(ggplot2)
# LAD profiles derived from normalized ALS data after applying [lad.profile()] function
data_path <- file.path(system.file("extdata", package = "LadderFuelsR"), "LAD_profiles.txt")
LAD_profiles <- read.table(data_path, sep = "\t", header = TRUE)
LAD_profiles$treeID <- factor(LAD_profiles$treeID)

# Tree metrics derived from get_cum_break() function
cum_LAD_metrics_path <- file.path(system.file("extdata", package = "LadderFuelsR"), "8_cbh_breaking_point_lad.txt")
cum_LAD_metrics <- read.table(cum_LAD_metrics_path, sep = "\t", header = TRUE)
cum_LAD_metrics$treeID <- factor(cum_LAD_metrics$treeID)

# Generate cumulative LAD plots
plots_trees_cumlad <- get_plots_cumm(LAD_profiles, cum_LAD_metrics)

# Save plots for each tree
for (name in names(plots_trees_cumlad)) {
  plots <- plots_trees_cumlad[[name]]

  if (!is.null(plots)) {
    #print(paste("Saving plot for tree:", name))
    ggsave(file.path(system.file("extdata", package = "LadderFuelsR"),  paste0( name, "_cumm_Hcbh", ".tiff")), plot = plots)
  }}

par(mfrow = c(1, 1))
# Plot the first panel
plot(plots_trees_cumlad[[1]],width = 6, height = 6)

```
![](https://github.com/carlos-alberto-silva/rGEDI/blob/master/readme/fig3.png)

## 21. Joining fuel ladder properties with crown polygons

```{r Joining crown polygons and ladder fuels metrics, echo=TRUE, message=FALSE, warning=FALSE}

# Tree metrics derived from get_layers_lad() function
LAD_gt25p_path <- file.path(system.file("extdata", package = "LadderFuelsR"), "7_fuels_lad_gt25perc.txt")
fuels_LAD_metrics <- read.table(LAD_gt25p_path, sep = "\t", header = TRUE)
fuels_LAD_metrics$treeID <- factor(fuels_LAD_metrics$treeID)

# No overlapping crown polygins (output from step 6)
trees_path <- file.path(system.file("extdata", package = "LadderFuelsR"), "crown2m_silva_convex_no_overlap.shp")
tree_file <- st_read(trees_path)
tree_file$treeID1 <- factor(tree_file$treeID)

crowns_properties<-merge (tree_file,fuels_LAD_metrics, by="treeID1")


crowns_properties$maxlad_Hcbh_factor <- cut(crowns_properties$maxlad_Hcbh, breaks = 5)

# Plotting with a discrete legend

red_to_green_palette <- colorRampPalette(c("red", "green"))

ggplot() +
  geom_sf(data = crowns_properties, aes(fill = maxlad_Hcbh_factor)) +
  scale_fill_manual(values = red_to_green_palette(5)) +
  theme_minimal() +
  labs(title = "Tree Crowns", fill = "maxlad_Hcbh")


```
![](https://github.com/carlos-alberto-silva/rGEDI/blob/master/readme/fig3.png)

# Acknowledgements
We gratefully acknowledge funding from XXXX, Carlos Silva was supported by the Carbon Monitoring System funding (CMS, grant 22-CMS22-0015). 

# Reporting Issues 
Please report any issue regarding the LadderFuelsR package to Dr. Olga Viedma (olga.viedma@uclm.es)

# Citing rICESat2Veg
Viedma,O.M;C.LadderFuelsR: LadderFuelsR: An R Package  for vertical fuel continuity analysis using LiDAR data.version 0.0.1, accessed on November. 22 2023, available at: <https://CRAN.R-project.org/package=LadderFuelsR>

# Disclaimer
**LadderFuelsR package comes with no guarantee, expressed or implied, and the authors hold no responsibility for its use or reliability of its outputs.**
