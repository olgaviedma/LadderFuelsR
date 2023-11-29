![](https://github.com/olgaviedma/LadderfuelsR/blob/master/Readme/LadderFuels_image.png)<br/>
![](https://github.com/olgaviedma/LadderfuelsR/blob/master/Readme/CHM pitfree 0.5 m-1.png)<br/>

[![CRAN](https://www.r-pkg.org/badges/version/LadderFuelsR)](https://cran.r-project.org/package=LadderFuelsR) ![Github](https://img.shields.io/badge/Github-0.0.1-green.svg) ![licence](https://img.shields.io/badge/Licence-GPL--3-blue.svg) ![Downloads](https://cranlogs.r-pkg.org/badges/grand-total/LadderFuelsR) [![Build Status](https://travis-ci.com/olgaviedma/LadderFuelsR.svg?token=Jqizwyc6gBxNafNccTdU&branch=master)](https://travis-ci.com/olgaviedma/LadderFuelsR)

**LadderFuelsR: An R Package for vertical fuel continuity analysis using LiDAR data.**

Authors: Olga Viedma, Carlos Silva and JM Moreno

The LadderFuelsR package is an automated tool for vertical fuel continuity analysis using LiDAR data that can be applied on multiple tree species and for large-scale studies. We include a suite of tools integrated into a standard workflow for deriving outputs meaningful to vertical fuel continuity analysis and canopy base height detection. The workflow consisted of: 1) calculating the vertical height profiles (VHP) of each segmented tree; 2) identifying gaps and fuel layers; 3) estimating the distance between fuel layers; and 4) retrieving the fuel layers base height (FBH) and depth. Additionally, other functions recalculate previous metrics after considering distances \> 1 m and calculate the canopy base height (CBH) as the FBH located at the largest- and at the last-distance. Moreover, the package calculates: i) the percentage of LAD comprised in each fuel layer, ii) remove fuel layers with LAD percentage \< 25 %, iii) recalculate the distances among the reminder ones, and iv) identify the CBH as the FBH with the highest LAD percentage. On the other hand, when the vertical height profiles (VHP) showed only one fuel layer, it identifies a possible CBH performing a segmented linear regression (breaking points) on the cumulative sum of LAD as a function of height. Finally, a collection of plotting functions is developed to represent: i) the VHP with the initial gaps and fuel layers; ii) the FBHs, depths and gaps with distances \> 1 m and, iii) the FBHs and depths after applying the breaking point method over trees with only one fuel layer.

# Getting Started

## Installation
``` r
#The CRAN version:
install.packages("LadderFuelsR")

# The development version:
#install.packages("remotes")
library(remotes)
install_github("https://github.com/olgaviedma/LadderFuelsR", dependencies = TRUE)

# loading LadderFuelsR package
library(LadderFuelsR)
```

## Required libraries

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
               sf,
               terra,
               data.table,
               rgdal,
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

## 1. Computing Canopy height model (CHM) using lidR package
```{r CHM pitfree 0.5 m, echo=TRUE, message=FALSE, warning=FALSE}

 LIDAR_dir <- file.path(system.file("extdata", package = "LadderFuelsR"))
 lidar_file<- lidR::readLAS(file.path(LIDAR_dir, "Eglin_zone1_clipped_000000.las"), filter = "-drop_z_below 0")

  chm_pitfree<- grid_canopy(lidar_file, res=0.5,pitfree( c(0,2,5,10,15,20,25,30,35,40), c(0,1.5), subcircle=0.15))
  chm_pitfree[chm_pitfree > 40] <- NA
  chm_pitfree[chm_pitfree < 0] <- 0
  chm_pitfree1 <- projectRaster(chm_pitfree, crs=26916)
  
col <- height.colors(25)
plot(chm_pitfree1,col=col)

```
<img align="right" src="https://github.com/olgaviedma/LadderFuelsR/blob/master/Readme/CHM%20pitfree%200.5%20m-1.png>

## 2.Detecting individual tree top from the lidar-derived CHM
```{r Tree tops detection, echo=TRUE, message=FALSE, warning=FALSE}

  # parameters
  ws= 2.5
  hmin = 2
  res=0.5
  ttops_multichm = find_trees(lidar_file, multichm(res = res, dist_2d = 2,ws= ws, layer_thickness = 0.3,dist_3d = 1, hmin = hmin))
  proj4string(ttops_multichm) <- CRS('+init=EPSG:26916')

# Create an rgl point cloud
x<-add_treetops3d(plot(lidar_file, bg = "white", size = 4), ttops_multichm)
# Customize the plot orientation
rgl.viewpoint(theta = 0, phi = 0, fov = 60, zoom = 0.75)
# Convert the rgl scene to an HTML widget
rglwidget(elementId = "x", width = 800, height = 600)

```
## 3. Individual tree crown deliniation (Silva et al. 2016)
```{r Crowns Silva, echo=TRUE, message=FALSE, warning=FALSE}

  algo_silva1 <-silva2016(chm_pitfree1, ttops_multichm, max_cr_factor = 0.6, exclusion = 0.3, ID = "treeID")
  crowns_silva_las1 <-segment_trees(lidar_file, algo_silva1, attribute = "treeID", uniqueness = "incremental")
  crowns_silva_las2<-filter_poi(crowns_silva_las1, !is.na(treeID))
  
my_palette <- colorRampPalette(col)
x1<-plot(crowns_silva_las2, color = "treeID", pal = my_palette, bg = "white")

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
     z_max = max(z),# max height
     z_min = min(z),# min height
     z_mean = mean(z),# mean height
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
 crowns_silva_filter<-filter_poi(crowns_silva_las2, Z >= 1)
  
  metrics1 = crown_metrics(crowns_silva_filter,func = .stdtreemetrics, geom = "convex")
  crown_diam<-data.frame(sqrt(metrics1$convhull_area/ pi) * 2)
  names(crown_diam)<-"crown_diam"
  metrics2 = crown_metrics(crowns_silva_filter,func = ccm, geom = "convex") #concave
  metrics_all <- dplyr::bind_cols(list(metrics1,crown_diam,metrics2))
  metrics_all1 <- metrics_all[,c(1:4,6,10:21)]
  names(metrics_all1)<-c("treeID", "Z", "npoints", "convhull_area", "crown_diam", "z_max", "z_min", "z_mean","z_sd", "z_q1","z_q5", "z_q25","z_q50","z_q75", "z_q95", "crr", "geometry" )
  
  tree_crowns <- st_as_sf(metrics_all1)
  
ttops1<-st_as_sf(ttops_multichm)
crowns1<-st_as_sf(tree_crowns)
ttops_within_crowns <- st_intersection(ttops1, crowns1)

# Set the size of the plotting device
par(mfrow = c(1, 1), mar = c(1, 1, 1, 1), pin = c(5, 4))
plot(st_geometry(crowns1), pch = 16, col = "green")
plot(ttops_within_crowns, add = TRUE, pch= 16, col = "darkblue", main = "Tree tops over the crowns")
```
## 6.Crop las files with crown polygons
```{r cropLAS files with no overlapping crowns, echo=TRUE, message=FALSE, warning=FALSE}

dir_out<-file.path(system.file("extdata", package = "leafR"))

  trees_ID <- tree_crowns %>% dplyr::select(treeID)
  n <- nrow(trees_ID)
  
  crown_cort <- vector("list", length=n)
  
  for (i in 1:n) {
    kk <- trees_ID[i,]
    crown_cort[[i]] = clip_roi(crowns_silva_las2, kk)
  
    lidR::writeLAS(crown_cort[[i]], paste0(dir_out,"/", i, "_CROWN.las"))
  }

my_palette <- colorRampPalette(col)
x2<-plot(crown_cort[[1]], color = "Z", pal = my_palette, bg = "black", size = 2.5)

# Customize the plot orientation
rgl.viewpoint(theta = 0, phi = 0, fov = 60, zoom = 0.75)

# Convert the rgl scene to an HTML widget
rglwidget(elementId = "x2", width = 400, height = 600)
```
## 7.LAI-LAD metrics by Trees
```{r LAI and LAD tree metrics, echo=TRUE, message=FALSE, warning=FALSE}

las_dir1 <- file.path(system.file("extdata", package = "leafR"))
las_list1 <- list.files(las_dir1, pattern = "*_CROWN.las", full.names = TRUE, ignore.case = TRUE)

# create a vector to hold the file names of .las files with more than 10 points
files_with_more_than_10_points <- c()

# loop through each file
for (file in las_list1) {
  las_data <- lidR::readLAS(file)
  las_data1<-filter_poi(las_data, Z >= 1)
  
  # skip to next file if there was a problem reading
  if (is.null(las_data1)) next
  
  # check if it contains more than three points
  if (las_data1@header$`Number of point records` > 10) {
  files_with_more_than_10_points <- c(files_with_more_than_10_points, file)
  }
}

# Creates a data frame of the 3D voxels information (xyz) with Leaf Area Density values
short_name1<-NULL
profile_list<-NULL
lidar_lai_list<-NULL
understory_lai_list<-NULL
LAHV_metric_list<-NULL

for (j in seq_along(files_with_more_than_10_points)){
  
  short_name<-stri_sub(files_with_more_than_10_points[j], 1, -5)
  short_name1<-gsub(".*/","",short_name)
  
  normlas_file<-files_with_more_than_10_points[[j]]

  VOXELS_LAD = lad.voxels(normlas_file, grain.size = 2)
  
  lad_profile = lad.profile(VOXELS_LAD, relative = F)
  lai_tot = lai(lad_profile)
  understory_lai <- lai(lad_profile, min = 0.3, max = 2.5)
  LAHV_metric<- LAHV(lad_profile, LAI.weighting = FALSE, height.weighting = FALSE)
  
  lad_profile1 = data.frame(lad_profile, treeID = short_name1)
  lai_tot1 = data.frame(lai_tot, treeID = short_name1)
  understory_lai1 = data.frame(understory_lai, treeID = short_name1)
  LAHV_metric1 = data.frame(LAHV_metric, treeID = short_name1)
  
  profile_list<-rbind(profile_list, lad_profile1)
  lidar_lai_list<-rbind(lidar_lai_list,lai_tot1)
  understory_lai_list <-rbind(understory_lai_list,understory_lai1)
  LAHV_metric_list<-rbind(LAHV_metric_list,LAHV_metric1)
}

head(profile_list,10)
```
## 8.Depurating Tree LAD profiles
```{r depurating LAD databases, echo=TRUE, message=FALSE, warning=FALSE}

cols <- c('treeID')
profile_list[cols] <- lapply(profile_list[cols], function (x) as.factor(x))
profile_list$lad<-round(profile_list$lad,digits = 4)

cases <- data.frame(table(profile_list$treeID))
cases1 <-cases[cases$Freq > 5, ]
names(cases1)<-c("treeID", "Freq")

profile_list1 <- profile_list[profile_list$treeID %in% cases1$treeID, ]
profile_list2 <- data.frame(profile_list1 %>% replace(is.na(.), 0.01))

```
## 9.Gaps and Fuel Layers Base Height (FBH)
```{r Gaps and Fuel layers Base Height (fbh), echo=TRUE, message=FALSE, warning=FALSE}

# LAD profiles derived from normalized ALS data after applying [lad.profile()] function from leafR package
profile_list2$treeID <- factor(profile_list2$treeID)

trees_name1 <- as.character(profile_list2$treeID)
trees_name2 <- factor(unique(trees_name1))

gaps_fbhs_list<-list()
for (i in levels(trees_name2)) {
  tree2 <- profile_list2 |> dplyr::filter(treeID == i)
  gaps_fbhs <- get_gaps_fbhs(tree2)
  gaps_fbhs_list[[i]] <- gaps_fbhs
}

gaps_fbhs_list1 <- dplyr::bind_rows(gaps_fbhs_list)
gaps_fbhs_list1$treeID <- factor(gaps_fbhs_list1$treeID)

# Remove the row with all NA values from the original data frame
# First remove "treeID" and "treeID1" columns
gaps_fbhs_list1_no_treeID <- gaps_fbhs_list1[, -which(names(gaps_fbhs_list1) == c("treeID","treeID1"))]
# Check if any row has all NA values
rows_with_all_NA_or_zero <- apply(gaps_fbhs_list1_no_treeID, 1, function(row) all(is.na(row) | row == 0))
# Get the row index with all NA values
row_index <- which(rows_with_all_NA_or_zero)

# Remove the row with all NA values from the original data frame
if (length(row_index) > 0) {
  gaps_fbhs_metrics <- gaps_fbhs_list1[-row_index, ]
} else {
  gaps_fbhs_metrics <- gaps_fbhs_list1
}
rownames(gaps_fbhs_metrics) <- NULL
head(gaps_fbhs_metrics)
```
## 10.Distance between Fuel Layers
```{r Distances (and their heights) between fuel layers, echo=TRUE, message=FALSE, warning=FALSE}

# Tree metrics derived from get_gaps_fbhs() function

numeric_vars <- setdiff(names(gaps_fbhs_metrics), c("treeID", "treeID1"))
gaps_fbhs_metrics[numeric_vars] <- lapply(gaps_fbhs_metrics[numeric_vars], function(x) as.numeric(ifelse(x == "NA", NA, x)))

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
distance_metrics <- dplyr::bind_rows(metrics_distance_list)
distance_metrics <- distance_metrics[, order(names(distance_metrics))]
rownames(distance_metrics) <- NULL
head(distance_metrics)
```
## 11.Fuel Layers Depth
```{r Distane between fuel layers, echo=TRUE, message=FALSE, warning=FALSE}

library(dplyr)
library(magrittr)

# LAD profiles derived from normalized ALS data after applying [lad.profile()] function
profile_list2$treeID <- factor(profile_list2$treeID)
# Tree metrics derived from get_distance() function
distance_metrics$treeID <- factor(distance_metrics$treeID)

metrics_depth_list <- list()

for (i in levels(profile_list2$treeID)){

  tree1 <- profile_list2 |> dplyr::filter(treeID == i)
  tree2 <- distance_metrics |> dplyr::filter(treeID == i)

  # Get depths for each tree
  metrics_depth <- get_depths(tree1, tree2)
  metrics_depth_list[[i]] <- metrics_depth
}

# Combine the individual data frames
depth_metrics <- dplyr::bind_rows(metrics_depth_list)

depth_metrics <- depth_metrics[, order(names(depth_metrics))]
rownames(depth_metrics) <- NULL
head(depth_metrics)
```
## 12.Plot Gaps and Fuel Layers Base Height (FBH)
```{r Plots Gaps and Fuel layers Base Height (fbh), echo=TRUE, message=FALSE, warning=FALSE}

library(LadderFuelsR)
library(ggplot2)
library(lattice)

# LAD profiles derived from normalized ALS data after applying [lad.profile()] function
profile_list2$treeID <- factor(profile_list2$treeID)
# Tree metrics derived from get_depths() function
depth_metrics$treeID <- factor(depth_metrics$treeID)

# Generate plots for gaps and fbhs
plots_gaps_fbhs <- get_plots_gap_fbh(profile_list2, depth_metrics)

par(mfrow = c(2, 2))
# Plot in RED are the GAPS and in GREEN the FBHs
plot(plots_gaps_fbhs[[1]])
plot(plots_gaps_fbhs[[2]])
plot(plots_gaps_fbhs[[3]])
plot(plots_gaps_fbhs[[4]])
```
## 13.Fuel Layers Base Height (FBH) after removing distances = 1
```{r Fuels base height after removing distances equal 1 m, echo=TRUE, message=FALSE, warning=FALSE}

library(SSBtools)
library(dplyr)
library(magrittr)

# Tree metrics derived from get_depths() function
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
fbh_metrics_corr <- dplyr::bind_rows(fbh_corr_list)
fbh_metrics_corr$treeID <- factor(fbh_metrics_corr$treeID)

# Reorder columns
# Get original column names
original_column_names <- colnames(fbh_metrics_corr)

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
fbh_metrics_corr <- fbh_metrics_corr[, new_order]
rownames(fbh_metrics_corr) <- NULL
head(fbh_metrics_corr)
```
## 14.Fuel Layers Depth after removing distances = 1
```{r Fuel layers depth after removinG distances equal 1 m, echo=TRUE, message=FALSE, warning=FALSE}

library(dplyr)
library(magrittr)
library(tidyr)

# Tree metrics derived from get_real_fbh() function
fbh_metrics_corr$treeID <- factor(fbh_metrics_corr$treeID)

trees_name1 <- as.character(fbh_metrics_corr$treeID)
trees_name2 <- factor(unique(trees_name1))

depth_metrics_corr_list <- lapply(levels(trees_name2), function(i) {
  # Filter data for each tree
  tree2 <- fbh_metrics_corr |> dplyr::filter(treeID == i)
  # Get real depths for each tree
  get_real_depths(tree2)
})

# Combine depth values for all trees
depth_metrics_corr <- dplyr::bind_rows(depth_metrics_corr_list)
rownames(depth_metrics_corr) <- NULL
head(depth_metrics_corr)
```
## 15.Distance between Fuel Layers after removing distances = 1 and CBH based ON MAX- AND LAST- Distance
```{r Fuel layers distances after removing distances equal 1 m, echo=TRUE, message=FALSE, warning=FALSE}

library(dplyr)
library(magrittr)
library(stringr)

# Tree metrics derived from get_real_depths() function
depth_metrics_corr$treeID <- factor(depth_metrics_corr$treeID)

trees_name1 <- as.character(depth_metrics_corr$treeID)
trees_name2 <- factor(unique(trees_name1))

distance_metrics_corr_list <- lapply(levels(trees_name2), function(i) {
  # Filter data for each tree
  tree2 <- depth_metrics_corr |> dplyr::filter(treeID == i)
  # Get effective gap for each tree
  get_effective_gap(tree2)
})

# Combine the individual data frames
distances_metrics_corr <- dplyr::bind_rows(distance_metrics_corr_list)

# =======================================================================#
# REORDER COLUMNS:
# =======================================================================#
# Get original column names
original_column_names <- colnames(distances_metrics_corr)

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
distances_metrics_corr1 <- distances_metrics_corr[, new_order]
# Unlist the data frame
distances_metrics_corr2 <- as.data.frame(lapply(distances_metrics_corr1, function(x) unlist(x)))
rownames(distances_metrics_corr2) <- NULL
head(distances_metrics_corr2)
```
## 16.Fuels LAD percentage (greater than 25 percent) and CBH based on MAX LAD percentage
```{r Fuels LAD percentage and canopy base height (CBH) based on maximum LAD percentage (distances greater than 1 m), echo=TRUE, message=FALSE, warning=FALSE}

library(dplyr)
library(magrittr)

# LAD profiles derived from normalized ALS data after applying [lad.profile()] function
profile_list2$treeID <- factor(profile_list2$treeID)

# Tree metrics derived from get_effective_gap() function
distances_metrics_corr2$treeID <- factor(distances_metrics_corr2$treeID)

trees_name1 <- as.character(distances_metrics_corr2$treeID)
trees_name2 <- factor(unique(trees_name1))

LAD_metrics1 <- list()
LAD_metrics2 <- list()

for (i in levels(trees_name2)) {
  # Filter data for each tree
  tree1 <- profile_list2 |> dplyr::filter(treeID == i)
  tree2 <- distances_metrics_corr2 |> dplyr::filter(treeID == i)

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
fuels_LAD_metrics <- list()

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
  fuels_LAD_metrics[[i]] <- LAD_metrics_all
 }
rownames(fuels_LAD_metrics[[1]]) <- NULL
rownames(fuels_LAD_metrics[[2]]) <- NULL

head(fuels_LAD_metrics[[1]])
head(fuels_LAD_metrics[[2]])
```
## 17.Plot Fuel Layers with LAD percentage greater than 25 percent and CBH based on MAX LAD percentage
```{r Plots of fuel layers with LAD percentage greater than 25 and the canopy base height (CBH) based on the maximum LAD percentage, echo=TRUE, message=FALSE, warning=FALSE}

library(ggplot2)

# LAD profiles derived from normalized ALS data after applying [lad.profile()] function
profile_list2$treeID <- factor(profile_list2$treeID)
# Tree metrics derived from get_layers_lad() function
LAD_gt25p <- fuels_LAD_metrics[[2]]

trees_name1 <- as.character(LAD_gt25p$treeID)
trees_name2 <- factor(unique(trees_name1))

# Generate plots for fuels LAD metrics
plots_trees_LAD <- get_plots_cbh_LAD(profile_list2, LAD_gt25p)

par(mfrow = c(2, 2))
plot(plots_trees_LAD[[1]])
plot(plots_trees_LAD[[2]])
plot(plots_trees_LAD[[3]])
plot(plots_trees_LAD[[4]])
```
## 18.CBH based on the Breaking Point method and LAD percentage
```{r CBH and the LAD percentage below and above the CBH using the breaking point method, echo=TRUE, message=FALSE, warning=FALSE}

library(dplyr)
library(magrittr)

# LAD profiles derived from normalized ALS data after applying [lad.profile()] function
profile_list2$treeID <- factor(profile_list2$treeID)

# Tree metrics derived from get_effective_gap() function
distances_metrics_corr2$treeID <- factor(distances_metrics_corr2$treeID)

trees_name1 <- as.character(distances_metrics_corr2$treeID)
trees_name2 <- factor(unique(trees_name1))

cum_LAD_metrics_list <- list()

for (i in levels(trees_name2)) {
  # Filter data for each tree
  tree1 <- profile_list2 |> dplyr::filter(treeID == i)
  tree2 <- distances_metrics_corr2 |> dplyr::filter(treeID == i)

  # Get cumulative LAD metrics for each tree
  cum_LAD_metrics_all <- get_cum_break(tree1, tree2)
  cum_LAD_metrics_list[[i]] <- cum_LAD_metrics_all
}

# Combine the individual data frames
cum_LAD_metrics <- dplyr::bind_rows(cum_LAD_metrics_list)

# =======================================================================#
# REORDER COLUMNS
# =======================================================================#

# Get original column names
original_column_names <- colnames(cum_LAD_metrics)

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
cum_LAD_metrics <- cum_LAD_metrics[, new_order]

## when % LAD is < 75 % below or above the BP (Breaking Point), Hcbh1 is not calculated and only it is given the Hcbh_bp (Hcbh from the BP)
rownames(cum_LAD_metrics) <- NULL
head(cum_LAD_metrics)
```
## 19.Plot CBH based on the Breaking Point method and LAD percentage
```{r Plots of the CBH and the LAD percentage below and above the CBH using the breaking point method, echo=TRUE, message=FALSE, warning=FALSE}

library(ggplot2)

# LAD profiles derived from normalized ALS data after applying [lad.profile()] function
profile_list2$treeID <- factor(profile_list2$treeID)

# Tree metrics derived from get_cum_break() function
cum_LAD_metrics$treeID <- factor(cum_LAD_metrics$treeID)

# Generate cumulative LAD plots
plots_trees_cumlad <- get_plots_cumm(profile_list2, cum_LAD_metrics)

par(mfrow = c(2, 2))
plot(plots_trees_cumlad[[1]])
plot(plots_trees_cumlad[[2]])
plot(plots_trees_cumlad[[3]])
plot(plots_trees_cumlad[[4]])

```
## 20. Joinning Fuel ladder properties with Crown polygons
```{r Joining crown polygons and ladder fuels metrics, echo=TRUE, message=FALSE, warning=FALSE}

# Tree metrics derived from get_layers_lad() function
LAD_gt25p <- fuels_LAD_metrics[[2]]
LAD_gt25p$treeID1 <- factor(LAD_gt25p$treeID1)

# crown polygons (output from step 4)
tree_crowns$treeID1 <- factor(tree_crowns$treeID)

crowns_properties<-merge (tree_crowns,LAD_gt25p, by="treeID1")

crowns_properties$maxlad_Hcbh_factor <- cut(crowns_properties$maxlad_Hcbh, breaks = 5)

# Plotting with a discrete legend

palette <- colorRampPalette(c("orange", "dark green"))

ggplot() +
  geom_sf(data = crowns_properties, aes(fill = maxlad_Hcbh_factor)) +
  scale_fill_manual(values = palette(5)) +
  theme_minimal() +
  labs(title = "Tree Crowns", fill = "maxlad_Hcbh")

```

# Acknowledgements

We gratefully acknowledge funding from project INFORICAM (PID2020-119402RB-I00), funded by the Spanish MCIN/AEI/ 10.13039/501100011033 and by the “European Union NextGenerationEU/PRTR”. Carlos Silva was supported by the NASA's Carbon Monitoring System funding (CMS, grant 22-CMS22-0015).

# Reporting Issues

Please report any issue regarding the LadderFuelsR package to Dr. Olga Viedma ([olga.viedma\@uclm.es](mailto:olga.viedma@uclm.es){.email})

# Citing LadderFuelsR

Viedma,O.;C.LadderFuelsR: LadderFuelsR: An R Package for vertical fuel continuity analysis using LiDAR data.version 0.0.1, accessed on November. 22 2023, available at: <https://CRAN.R-project.org/package=LadderFuelsR>

# Disclaimer

**LadderFuelsR package comes with no guarantee, expressed or implied, and the authors hold no responsibility for its use or reliability of its outputs.**
