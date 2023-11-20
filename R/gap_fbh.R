#' Gaps and Fuel layers Base Height (fbh)
#' @description This function calculates gaps and fuel layers base height (fbh) as the difference in percentiles between consecutive LAD values along the vertical tree profile (VTP).
#' Negative differences are linked to gaps and positive differences to fuel base height.
#' @usage get_gaps_fbhs (LAD_profiles)
#'
#' @param LAD_profiles original tree Leaf Area Index (LAD) profile (output of [lad.profile()] function in the \emph{leafR} package.
#' An object of the class text
#' @return A data frame giving the height of gaps and fuel layers bases in meters.
#' @author Olga Viedma.
#'
#' @details
#' # List of tree metrics:
#' \itemize{
#' \item treeID: tree ID with strings and numeric values
#' \item treeID1: tree ID with only numeric values
#' \item cbh - Height of the fuel layer base height (m)
#' \item gap - Height of gap between fuel layers (m)
#' \item gap_lad: LAD value in the gap height
#' \item gap_perc - Percentage of LAD in the gap height
#' \item cbh_lad - LAD value in the fuel base height
#' \item cbh_perc - Percentage of LAD in the fuel base height
#' \item max_height - Maximum height of the tree profile
#' }
#' @examples
#' ## Not run:
#' library(dplyr)
#' library(magrittr)
#'
#' # LAD profiles derived from normalized ALS data after applying [lad.profile()] function
#' data_path <- file.path(system.file("extdata", package = "LadderFuelsR"), "LAD_profiles.txt")
#' LAD_profiles <- read.table(data_path, sep = "\t", header = TRUE)
#' LAD_profiles$treeID <- factor(LAD_profiles$treeID)
#'
#' trees_name1 <- as.character(LAD_profiles$treeID)
#' trees_name2 <- factor(unique(trees_name1))
#'
#' metrics_precentile_list<-list()
#' for (i in levels(trees_name2)) {
#'   tree2 <- LAD_profiles |> dplyr::filter(treeID == i)
#'   metrics_precentil <- get_gaps_fbhs(tree2)
#'   metrics_precentile_list[[i]] <- metrics_precentil
#' }
#'
#' metrics_all_percentil <- dplyr::bind_rows(metrics_precentile_list)
#' metrics_all_percentil$treeID <- factor(metrics_all_percentil$treeID)
#'
#' # Remove the row with all NA values from the original data frame
#' # First remove "treeID" and "treeID1" columns
#' tree_metrics_no_treeID <- metrics_all_percentil[, -which(names(metrics_all_percentil) == c("treeID","treeID1"))]
#'
#' # Check if any row has all NA values
#' rows_with_all_NA_or_zero <- apply(tree_metrics_no_treeID, 1, function(row) all(is.na(row) | row == 0))
#'
#' # Get the row index with all NA values
#' row_index <- which(rows_with_all_NA_or_zero)
#'
#' # Remove the row with all NA values from the original data frame
#' if (length(row_index) > 0) {
#'   tree_metrics_filtered <- metrics_all_percentil[-row_index, ]
#' } else {
#'   tree_metrics_filtered <- metrics_all_percentil
#' }
#' write.table(tree_metrics_filtered, file= file.path(system.file("extdata", package = "LadderFuelsR"), "1_gaps_fbhs_metrics.txt"),
#' sep = "\t",row.names = FALSE)
#' ## End(Not run)
#' @export get_gaps_fbhs
#' @importFrom dplyr group_by summarise mutate arrange
#' @importFrom magrittr %>%
get_gaps_fbhs<- function (LAD_profiles) {

  df<- LAD_profiles

  treeID<-"treeID"
  #print(paste("treeID:", df[[treeID]][1]))  # Debugging line

    df$height<-as.numeric(df$height)
     df$treeID<-factor(df$treeID)
trees_name1a<- as.character(df$treeID)
trees_name3<- factor(unique(trees_name1a))

df_ord<-df[with(df, order(lad)), ]

all_equal <- length(unique(df_ord$lad)) == 1

if(all_equal) {

  treeID<-unique(factor(df$treeID))

  crown_height_data <- data.frame(NA)
  names(crown_height_data) <- "cbh0"
  gaps_height_data <- data.frame(NA)
  names(gaps_height_data) <- "gap0"
  distance_data <- data.frame(NA)
  names(distance_data) <- "dist0"
  depth0 <- data.frame(NA)
  names(depth0) <- "depth0"
  depth_data <- data.frame(NA)
  names(depth_data) <- "depth01"

  metrics_tree<-cbind.data.frame(crown_height_data, gaps_height_data, distance_data, depth0,depth_data, treeID)

} else {

  PERCENTIL_Z <- df %>%
    group_by(treeID) %>%
    summarise(
      P5 = quantile(lad, probs = 0.05, na.rm = TRUE),
      P25 = quantile(lad, probs = 0.25, na.rm = TRUE),
      P50 = quantile(lad, probs = 0.50, na.rm = TRUE),
      P75 = quantile(lad, probs = 0.75, na.rm = TRUE),
      P90 = quantile(lad, probs = 0.90, na.rm = TRUE),
      P95 = quantile(lad, probs = 0.95, na.rm = TRUE),
      P99 = quantile(lad, probs = 0.99, na.rm = TRUE)
    )

  x1<-df$height
  y1<-df$lad

  # Identify missing and infinite values in x and y
missing_x <- is.na(x1)
missing_y <- is.na(y1) | is.infinite(y1)

# Remove missing and infinite values from x and y
x <- x1[!missing_x & !missing_y]
y <- y1[!missing_x & !missing_y]


fit <- smooth.spline(x, y) # Fit a smoothing spline
y_second_deriv <- predict(fit, fit$x, deriv = 2) # Calculate the second derivative

#plot(y_second_deriv$x, y_second_deriv$y, type = "l", lwd = 2, col = "red", xlab = "X", ylab = "f''(x)", main = "Plot of second derivative of f(x)")

base_2drivative<-data.frame(do.call(cbind , y_second_deriv)) # convert the list into a dataframe
base_2drivative$y<-round(base_2drivative$y, digits=10)

critical_points<-base_2drivative[,2]# Extract the values of the second derivative
base_2drivative2<-cbind.data.frame(df[,c(1:3)],critical_points)

gaps_perc<- with(base_2drivative2,
          ifelse(lad <= PERCENTIL_Z$P5 , "5",
          ifelse(lad > PERCENTIL_Z$P5 & lad <= PERCENTIL_Z$P25, "25",
         ifelse(lad > PERCENTIL_Z$P25 & lad <= PERCENTIL_Z$P50, "50",
         ifelse(lad > PERCENTIL_Z$P50 & lad <= PERCENTIL_Z$P75, "75",
         ifelse(lad > PERCENTIL_Z$P75 & lad <= PERCENTIL_Z$P90, "90",
          ifelse(lad > PERCENTIL_Z$P90 & lad <= PERCENTIL_Z$P95, "95",
          ifelse(lad > PERCENTIL_Z$P95 & lad <= PERCENTIL_Z$P99, "99",
          ifelse(lad >  PERCENTIL_Z$P99, "100",NA)) )))))))

gaps_perc1 <- data.frame(percentil = as.numeric(gaps_perc))
gaps_perc2<-cbind.data.frame(base_2drivative2,gaps_perc1)

percentil<-gaps_perc2$percentil

#######################################
##GAPS
#######################################

diffs <- diff(percentil)  # calculate differences between adjacent elements
ind <- which(diffs < 0) + 1  # get indices of elements with a difference less than 0, add 1 to get indices of the corresponding consecutive elements

group <- cumsum(c(1, diff(ind) != 1)) # Use cumsum() to create a grouping variable for consecutive values
split_x <- split(ind, group)# split the vector into a list of vectors based on the grouping variable

gaps1 <- lapply(split_x, function (x) gaps_perc2[x,]) ## extract values from original df
gaps2<-do.call(rbind,gaps1)

#####
varname <- "lad"# define the variable to use for subsetting: minimum LAD
minvals <- lapply(gaps1, function(df) min(df[[varname]])) #  the minimum value of the variable in each data frame

subsetlist <- lapply(gaps1, function(df) {
  minval <- min(df[[varname]])
  subset(df, df[[varname]] == minval)
}) # subset the list with the min LAD

gaps5<-do.call(rbind,subsetlist)
gaps5a<-data.frame(t(gaps5))

# Subset the dataframe based on the condition

gaps5b<-dplyr::filter(gaps5, percentil <= 25)

# filter only percentil == 5
  gaps_perc_5 <- gaps_perc2 %>%
  filter(percentil == 5)

# Create an empty list to store the subsetted data frames
subset_list <- list()

# Identify consecutive groups based on the "height" column
consecutive_groups <- cumsum(c(TRUE, diff(gaps_perc_5$height) != 1))

# Iterate over each consecutive group
for (group in unique(consecutive_groups)) {
  # Subset the current group
  group_subset <- gaps_perc_5[consecutive_groups == group, ]

  # Add the first and last rows of the subset to the list
  subset_list[[group]] <- rbind(group_subset[1, ], group_subset[nrow(group_subset), ])
}

# Combine the subsetted data frames into a single data frame
result <- do.call(rbind, subset_list)
result_unique <- unique(result)


  ############################################

# filter only percentil <= 25
  gaps_perc_25 <- gaps_perc2 %>%
  filter(percentil <= 25)

# Create an empty list to store the subsetted data frames
subset_list1 <- list()

# Identify consecutive groups based on the "height" column
consecutive_groups1 <- cumsum(c(TRUE, diff(gaps_perc_25$height) != 1))

# Iterate over each consecutive group
for (group in unique(consecutive_groups1)) {
  # Subset the current group
  group_subset1 <- gaps_perc_25[consecutive_groups1 == group, ]

  # Add the first and last rows of the subset to the list
  subset_list1[[group]] <- rbind(group_subset1[1, ], group_subset1[nrow(group_subset1), ])
}

# Combine the subsetted data frames into a single data frame
result1 <- do.call(rbind, subset_list1)
result_unique1 <- unique(result1)

############################################

    gaps5i<-rbind.data.frame(gaps5b,result_unique,result_unique1)

     gaps5j <- gaps5i[!duplicated(gaps5i), ]
gaps5k<-gaps5j[with(gaps5j, order(height)), ]

   gaps6<-data.frame(t(gaps5k))


#######################################
##CROWNS-BASE
#######################################

diffs <- diff(percentil)  # calculate differences between adjacent elements
diff_vec <- c(percentil[1], diff(percentil))

ind1 <- which(diff_vec > 0)  # get indices of elements with a difference greater than 0, add 1 to get indices of the corresponding consecutive elements

group1 <- cumsum(c(1, diff(ind1) != 1)) # Use cumsum() to create a grouping variable for consecutive values
split_y <- split(ind1, group1)# split the vector into a list of vectors based on the grouping variable

crown1 <- lapply(split_y, function (x) df[x,]) ## extract values from original df
crown2<-do.call(rbind,crown1)

crown2a <- lapply(split_y, function (x) gaps_perc2[x,]) ## extract values from original df
crown2b<-do.call(rbind,crown2a)

#######################################

subset_crown2a <- lapply(crown2a, function(x) x[x$percentil > 5, ])

crown3<-do.call(rbind,subset_crown2a)
crown3a<-data.frame(t(crown3))

varname <- "lad"# define the variable to use for subsetting: minimum LAD
minvals1 <- lapply(subset_crown2a, function(df) min(df[[varname]])) #  the minimum value of the variable in each data frame

subsetlist1 <- lapply(subset_crown2a, function(df) {
  minval1 <- min(df[[varname]])
  subset(df, df[[varname]] == minval1)
}) # subset the list with the minimum LAD

crown3b<-do.call(rbind,subsetlist1)
crown3c<-data.frame(t(crown3b))


# filter only percentil > 5
perc_5 <- gaps_perc2 %>%
  filter(percentil > 5)


if (any(diff(perc_5$height) == 1)) {

# find the first and last consecutive rows based on the "height" column
consecutive_rows1 <- which(diff(perc_5$height) == 1)
first_consecutive_row1 <- min(consecutive_rows1)
last_consecutive_row1 <- max(consecutive_rows1)

# subset the original data frame using the first and last consecutive rows
crown_first_last <- perc_5[c(first_consecutive_row1, last_consecutive_row1), ]

crown3c<-rbind.data.frame(crown3b,crown_first_last)
crown3d <- crown3c[!duplicated(crown3c), ]
crown3e<-crown3d[with(crown3d, order(height)), ]

crown4<-data.frame(t(crown3e))

} else {

  crown4<-data.frame(t(crown3b))
}

############ ADAPT THE GAPS TO THE CBH (PREFEReNCE THE CBHs because there are some problems with distance calculus)

if (exists("crown3e") && !is.null(crown3e) && length(crown3e) > 0) {
# Select heights in gaps5k that are not in crown3e
gaps5l <- setdiff(gaps5k$height, crown3e$height)
gaps5m <- df[df$height %in% gaps5l, ]
 gaps6<-data.frame(t(gaps5m))

} else {
  gaps6<-data.frame(t(gaps5k))
}

#######################################
## EXTRACT LAD AND PERCENTIL FROM GAPS AND CBHS ##############33
#######################################
gaps6a<-data.frame(t(gaps6))
gaps6a$height <- as.numeric(gaps6a$height)
gaps_perc2$height <- as.numeric(gaps_perc2$height)

merged_gaps <- merge(gaps6a, gaps_perc2, by=c("height", "treeID"))
merged_gaps1<-merged_gaps[,-4]
merged_gaps1$type<-rep(c("gap"), nrow(merged_gaps1))

# Add a numeric suffix to each 'gap' value
merged_gaps1 <- merged_gaps1 %>%
  mutate(type = ifelse(type == "gap", paste0("gap", seq_along(type[type == "gap"])), type))
merged_gaps2<-data.frame(t(merged_gaps1))

gaps_lad <- data.frame(merged_gaps2[3,])
  gaps_lad1 <- as.data.frame(lapply(gaps_lad, as.numeric))
   colnames(gaps_lad1) <- paste0("gap_lad", seq_along(gaps_lad1))

   gaps_perc <- data.frame(merged_gaps2[5,])
  gaps_perc1 <- as.data.frame(lapply(gaps_perc, as.numeric))
   colnames(gaps_perc1) <- paste0("gap_perc", seq_along(gaps_perc1))

#######################################
crown4a<-data.frame(t(crown4))
crown4a$type<-rep(c("cbh"), nrow(crown4a))

# Add a numeric suffix to each 'cbh' value
merged_crown1 <- crown4a %>%
  mutate(type = ifelse(type == "cbh", paste0("cbh", seq_along(type[type == "cbh"])), type))
merged_crown2<-data.frame(t(merged_crown1))

crown_lad <- data.frame(merged_crown2[2,])
  crown_lad1 <- as.data.frame(lapply(crown_lad, as.numeric))
   colnames(crown_lad1) <- paste0("cbh_lad", seq_along(crown_lad1))

   crown_perc <- data.frame(merged_crown2[5,])
  crown_perc1 <- as.data.frame(lapply(crown_perc, as.numeric))
   colnames(crown_perc1) <- paste0("cbh_perc", seq_along(crown_perc1))

   #######################################
   #######################################

   # Check if gaps6 exists and has data
   if (!is.null(gaps6) && nrow(gaps6) > 0) {
     gaps_height_t <- gaps6[1,] %>%
       as.data.frame() %>%
       dplyr::mutate_all(as.numeric) %>%
       t() %>%
       as.data.frame() %>%
       dplyr::mutate(type = "gap")
   } else {
     gaps_height_t <- tibble(V1 = NA, type = "gap")
   }


   # Check if crown4 exists and contains data
   if (!is.null(crown4) && nrow(crown4) > 0) {
     crown_height_t <- crown4[1,] %>%
       as.data.frame() %>%
       dplyr::mutate_all(as.numeric) %>%
       t() %>%
       as.data.frame() %>%
       dplyr::mutate(type = "cbh")
   } else {
     crown_height_t <- tibble(V1 = NA, type = "cbh")
   }

   # Rename the column in gaps_height_t to match the column name in crown_height_t
   colnames(gaps_height_t)[colnames(gaps_height_t) == "V1"] <- "height"
   colnames(crown_height_t)[colnames(crown_height_t) == "V1"] <- "height"

   # Now you can combine the data frames
   combined_df <- rbind(gaps_height_t, crown_height_t)
   combined_df_ord <- combined_df[order(combined_df$height), ]

   combined_df_ord1<-na.omit(combined_df_ord)

   ##########################################################

   # Transpose the dataframe
   kk <- as.data.frame(t(combined_df_ord1))

   # Set the column names based on the second row and then remove that row
   colnames(kk) <- kk[2, ]
   kk <- kk[-2, ]

   # Rename columns to add indices to duplicate names
   names(kk) <- paste0(names(kk), ave(seq_along(names(kk)), names(kk), FUN = seq_along))

   # Convert the first row to numeric
   kk[1, ] <- as.data.frame(lapply(kk[1, ], as.numeric))

   # Make a copy of kk without the second row (which doesn't exist anymore)
   kk_copy <- kk

treeID<-unique(factor(df$treeID))

max_height<-data.frame(max(df$height))
names(max_height)="max_height"

    gap_cbh_metrics <- cbind.data.frame(treeID, kk_copy, gaps_lad1,crown_perc1, crown_lad1,max_height)

    gap_cbh_metrics$treeID1 <- sub("_.*", "", gap_cbh_metrics$treeID)
    gap_cbh_metrics$treeID1 <- as.numeric(gap_cbh_metrics$treeID1)
    gap_cbh_metrics <- arrange(gap_cbh_metrics, treeID1)
   }

return(gap_cbh_metrics)
}

