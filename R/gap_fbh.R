#' Gaps and Fuel layers Base Height (FBH)
#' @description This function calculates gaps and fuel layers base height (FBH) as the difference in percentiles between consecutive LAD values along the vertical tree profile (VTP).
#' Negative differences are linked to gaps and positive differences to fuel base height.
#' @usage get_gaps_fbhs (LAD_profiles, step=1, min_height=1.5,
#' perc_gap= 25, perc_base= 25, verbose=TRUE)
#' @param LAD_profiles original tree Leaf Area Density (LAD) profile (output of [lad.profile()] function in the \emph{leafR} package.
#' An object of the class text.
#' @param step Numeric value for the actual height bin step (in meters).
#' @param min_height Numeric value for the actual minimum base height (in meters).
#' @param perc_gap Numeric value of the percentile threshold used to identify gaps (default percentile 25th).
#' @param perc_base Numeric value of the percentile threshold used to identify fuels layers base height (default percentile 25th).
#' @param verbose Logical, indicating whether to display informational messages (default is TRUE).
#' @return A data frame giving the height of gaps and fuel layers bases in meters.
#' @author Olga Viedma, Carlos Silva, JM Moreno and A.T. Hudak
#'
#'@details
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
#'
#' @examples
#' library(magrittr)
#' library(dplyr)
#'
#' # LAD profiles derived from normalized ALS data after applying [lad.profile()] function
#' LAD_profiles <- read.table(system.file("extdata", "LAD_profiles.txt", package = "LadderFuelsR"),
#' header = TRUE)
#' LAD_profiles$treeID <- factor(LAD_profiles$treeID)
#'
#' trees_name1 <- as.character(LAD_profiles$treeID)
#' trees_name2 <- factor(unique(trees_name1))
#'
#' metrics_precentile_list1<-list()
#'
#' for (i in levels(trees_name2)) {
#' tree1 <- LAD_profiles |> dplyr::filter(treeID == i)
#' metrics_precentil <- get_gaps_fbhs(tree1, step=1,
#' min_height=1.5,
#' perc_gap= 25,perc_base= 25,
#' verbose=TRUE)
#' metrics_precentile_list1[[i]] <- metrics_precentil
#' }
#'
#' metrics_all_percentil <- dplyr::bind_rows(metrics_precentile_list1)
#' metrics_all_percentil$treeID <- factor(metrics_all_percentil$treeID)
#'
#' # Remove the row with all NA values from the original data frame
#' # First remove "treeID" and "treeID1" columns
#' no_treeID <- metrics_all_percentil[, -which(names(metrics_all_percentil) == c("treeID","treeID1"))]
#'
#' # Check if any row has all NA values
#' NA_or_zero <- apply(no_treeID, 1, function(row) all(is.na(row) | row == 0))
#'
#' # Get the row index with all NA values
#' row_index <- which(NA_or_zero)
#'
#' # Remove the row with all NA values from the original data frame
#' if (length(row_index) > 0) {
#' gap_cbh_metrics <- metrics_all_percentil[-row_index, ]
#' } else {
#' gap_cbh_metrics <- metrics_all_percentil
#' }
#' @importFrom dplyr select_if group_by summarise summarize mutate arrange rename rename_with filter slice slice_tail ungroup distinct
#' across matches row_number all_of vars bind_cols case_when left_join mutate if_else lag n_distinct
#' @importFrom segmented segmented seg.control
#' @importFrom magrittr %>%
#' @importFrom stats ave dist lm na.omit predict quantile setNames smooth.spline
#' @importFrom utils tail
#' @importFrom tidyselect starts_with everything one_of
#' @importFrom stringr str_extract str_match str_detect str_remove_all
#' @importFrom tibble tibble
#' @importFrom tidyr pivot_longer fill pivot_wider replace_na
#' @importFrom gdata startsWith
#' @importFrom ggplot2 aes geom_line geom_path geom_point geom_polygon geom_text geom_vline ggtitle coord_flip theme_bw
#' theme element_text xlab ylab ggplot xlim
#' @export
get_gaps_fbhs<- function (LAD_profiles, step=1,
                          min_height=1.5,
                          perc_gap= 25,perc_base= 25,
                          verbose=TRUE) {

  df <- LAD_profiles

  # check required columns exist in data and stop if not
  nms_df <- sort(names(df))
  nms_needed <- sort(c("treeID", "height", "lad"))
  nms_miss <- nms_needed[!nms_needed %in% nms_df]
  if(length(nms_miss)>0){
    stop(paste(
      "Columns not found in data. Supply missing columns:"
      , paste(nms_miss, collapse = ", ")
    ))
  }

  # reorder the columns so that it doesn't matter how the input data is structured
  df <- df %>% dplyr::relocate(height, lad, treeID)


  if(min_height==0){
    min_height <-0.5

    # Ensure the column starts with a negative value
    if (df$height[1] < min_height) {
      # Calculate the shift value
      shift_value <- abs(df$height[1])

      # Adjust the column to start from 0
      df$height <- df$height + shift_value
    }


    # Ensure the column starts with a negative value
    if (df$height[1] > min_height) {
      # Calculate the shift value
      shift_value1 <- abs(df$height[1])

      # Adjust the column to start from 0
      df$height <- df$height - shift_value1
    }
  }



    treeID <- "treeID"

    if (verbose) {
      message("Unique treeIDs:", paste(unique(df$treeID), collapse = ", "))
    }

    df$height <- as.numeric(df$height)
    df$treeID <- factor(df$treeID)
    trees_name1a <- as.character(df$treeID)
    trees_name3 <- factor(unique(trees_name1a))

    lad<-df$lad
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

  # Calculate percentiles for each treeID
  PERCENTIL_Z <- df %>%
    dplyr::group_by(treeID) %>%
    dplyr::summarise(across(lad, list(
      P5 = ~ quantile(.x, probs = 0.05, na.rm = TRUE),
      P10 = ~ quantile(.x, probs = 0.10, na.rm = TRUE),
      P15 = ~ quantile(.x, probs = 0.15, na.rm = TRUE),
      P20 = ~ quantile(.x, probs = 0.20, na.rm = TRUE),
      P25 = ~ quantile(.x, probs = 0.25, na.rm = TRUE),
      P30 = ~ quantile(.x, probs = 0.30, na.rm = TRUE),
      P35 = ~ quantile(.x, probs = 0.35, na.rm = TRUE),
      P40 = ~ quantile(.x, probs = 0.40, na.rm = TRUE),
      P45 = ~ quantile(.x, probs = 0.45, na.rm = TRUE),
      P50 = ~ quantile(.x, probs = 0.50, na.rm = TRUE),
      P55 = ~ quantile(.x, probs = 0.55, na.rm = TRUE),
      P60 = ~ quantile(.x, probs = 0.60, na.rm = TRUE),
      P65 = ~ quantile(.x, probs = 0.65, na.rm = TRUE),
      P70 = ~ quantile(.x, probs = 0.70, na.rm = TRUE),
      P75 = ~ quantile(.x, probs = 0.75, na.rm = TRUE),
      P80 = ~ quantile(.x, probs = 0.80, na.rm = TRUE),
      P85 = ~ quantile(.x, probs = 0.85, na.rm = TRUE),
      P90 = ~ quantile(.x, probs = 0.90, na.rm = TRUE),
      P95 = ~ quantile(.x, probs = 0.95, na.rm = TRUE),
      P99 = ~ quantile(.x, probs = 0.99, na.rm = TRUE)
    ))) %>%
    dplyr::ungroup()

  #print(PERCENTIL_Z)

  x1 <- df$height
  y1 <- df$lad

    # Identify missing and infinite values in x and y
    missing_x <- is.na(x1)
    missing_y <- is.na(y1) | is.infinite(y1)

    # Remove missing and infinite values from x and y
    x <- x1[!missing_x & !missing_y]
    y <- y1[!missing_x & !missing_y]

    # Fit a smoothing spline to the data
    fit <- smooth.spline(x, y)

    # Calculate the second derivative for all original x values
    y_second_deriv <- predict(fit, x, deriv = 2)

    # Convert the second derivative results into a data frame
    base_2drivative <- data.frame(x = y_second_deriv$x, y = round(y_second_deriv$y, digits = 10))

    # Merge the second derivative with the original data frame
    critical_points <- base_2drivative[, 2]  # Extract the values of the second derivative
    base_2drivative2 <- cbind.data.frame(df[, c(1:3)], critical_points)

    # Calculate gaps_perc using dplyr::case_when
    gaps_perc <- df %>%
      mutate(
        gaps_perc = case_when(
          lad <= PERCENTIL_Z$lad_P5 ~ "5",
          lad > PERCENTIL_Z$lad_P5 & lad <= PERCENTIL_Z$lad_P10 ~ "10",
          lad > PERCENTIL_Z$lad_P10 & lad <= PERCENTIL_Z$lad_P15 ~ "15",
          lad > PERCENTIL_Z$lad_P15 & lad <= PERCENTIL_Z$lad_P20 ~ "20",
          lad > PERCENTIL_Z$lad_P20 & lad <= PERCENTIL_Z$lad_P25 ~ "25",
          lad > PERCENTIL_Z$lad_P25 & lad <= PERCENTIL_Z$lad_P30 ~ "30",
          lad > PERCENTIL_Z$lad_P30 & lad <= PERCENTIL_Z$lad_P35 ~ "35",
          lad > PERCENTIL_Z$lad_P35 & lad <= PERCENTIL_Z$lad_P40 ~ "40",
          lad > PERCENTIL_Z$lad_P40 & lad <= PERCENTIL_Z$lad_P45 ~ "45",
          lad > PERCENTIL_Z$lad_P45 & lad <= PERCENTIL_Z$lad_P50 ~ "50",
          lad > PERCENTIL_Z$lad_P50 & lad <= PERCENTIL_Z$lad_P55 ~ "55",
          lad > PERCENTIL_Z$lad_P55 & lad <= PERCENTIL_Z$lad_P60 ~ "60",
          lad > PERCENTIL_Z$lad_P60 & lad <= PERCENTIL_Z$lad_P65 ~ "65",
          lad > PERCENTIL_Z$lad_P65 & lad <= PERCENTIL_Z$lad_P70 ~ "70",
          lad > PERCENTIL_Z$lad_P70 & lad <= PERCENTIL_Z$lad_P75 ~ "75",
          lad > PERCENTIL_Z$lad_P75 & lad <= PERCENTIL_Z$lad_P80 ~ "80",
          lad > PERCENTIL_Z$lad_P80 & lad <= PERCENTIL_Z$lad_P85 ~ "85",
          lad > PERCENTIL_Z$lad_P85 & lad <= PERCENTIL_Z$lad_P90 ~ "90",
          lad > PERCENTIL_Z$lad_P90 & lad <= PERCENTIL_Z$lad_P95 ~ "95",
          lad > PERCENTIL_Z$lad_P95 & lad <= PERCENTIL_Z$lad_P99 ~ "99",
          lad > PERCENTIL_Z$lad_P99 ~ "100",
          TRUE ~ NA_character_
        )
      )

    # Check the distribution of gaps_perc
   #print(table(gaps_perc$gaps_perc))


    # Combine results into final data frame
    gaps_perc1 <- gaps_perc %>% mutate(percentil = as.numeric(gaps_perc))
    gaps_perc2 <- bind_cols(base_2drivative2, percentil = gaps_perc1$percentil)

    percentil <- gaps_perc2$percentil

    #######################################
    ##GAPS
    #######################################

    diffs <- diff(percentil)  # calculate differences between adjacent elements
    ind <- which(diffs < 0)+1 # get indices of elements with a difference less than 0, add 1 to get indices of the corresponding consecutive elements

    group <- cumsum(c(1, diff(ind) != 1)) # Use cumsum() to create a grouping variable for consecutive values
    split_x <- split(ind, group)# split the vector into a list of vectors based on the grouping variable

    gaps1 <- lapply(split_x, function (x) gaps_perc2[x,]) ## extract values from original df
    gaps2<-do.call(rbind,gaps1)

    subsetlist <- lapply(gaps1, function(df) {
      minval <- min(df$lad[is.finite(df$lad)])
      subset(df, df$lad == minval)
    }) # subset the list with the min LAD

    gaps5<-do.call(rbind,subsetlist)
    gaps5a<-data.frame(t(gaps5))

    # Subset the dataframe based on the condition

    perc_gap_value <-perc_gap ## default 25

    gaps5b<-dplyr::filter(gaps5, percentil <= perc_gap_value)

        ############################################

    # filter only percentil == 5
    gaps_perc_5 <- gaps_perc2 %>%
      dplyr::filter(percentil == 5)

    # Create an empty list to store the subsetted data frames
    subset_list <- list()

    # Identify consecutive groups based on the "height" column
    consecutive_groups <- cumsum(c(TRUE, diff(gaps_perc_5$height) != step))

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
      dplyr::filter(percentil <= perc_gap_value)

    # Create an empty list to store the subsetted data frames
    subset_list1 <- list()

    # Identify consecutive groups based on the "height" column
    consecutive_groups1 <- cumsum(c(TRUE, diff(gaps_perc_25$height) != step))

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

    gaps5i<-rbind.data.frame(gaps5b,result_unique, result_unique1)

    gaps5j <- gaps5i[!duplicated(gaps5i), ]
    gaps5k<-gaps5j[with(gaps5j, order(height)), ]

    gaps6<-data.frame(t(gaps5k))


    #######################################
    ##CROWNS-BASE
    #######################################

    diffs <- diff(percentil)  # calculate differences between adjacent elements
    #diff_vec <- c(percentil[1], diff(percentil))
    #ind1 <- which(diff_vec > 5)  # (before was 0) get indices of elements with a difference greater than 0

    ind1 <- which(diffs > 5) + 1

    group1 <- cumsum(c(1, diff(ind1) != 1)) # Use cumsum() to create a grouping variable for consecutive values
    split_y <- split(ind1, group1)# split the vector into a list of vectors based on the grouping variable

    crown1 <- lapply(split_y, function (x) df[x,]) ## extract values from original df
    crown2<-do.call(rbind,crown1)

    crown2a <- lapply(split_y, function (x) gaps_perc2[x,]) ## extract values from original df
    crown2b<-do.call(rbind,crown2a)

    subsetlist1 <- lapply(crown2a, function(df) {
      if (any(is.finite(df$lad) & !is.na(df$lad))) {
        minval1 <- min(df$lad[is.finite(df$lad)])
        subset(df, df$lad == minval1)
      } else {
        NULL  # or any other value that indicates no valid subset
      }
    }) # subset the list with the minimum LAD

    crown2c<-do.call(rbind,subsetlist)

    ################################

    perc_base_value <-perc_base ## default 5

    subset_crown2a <- lapply(subsetlist1, function(x) x[x$percentil > perc_base_value, ])

    crown3<-do.call(rbind,subset_crown2a)


  ####################################################################
        # filter only percentil > 5

    perc_base_value <-perc_base ## default 5

    perc_5 <- gaps_perc2 %>%
      dplyr::filter(percentil > perc_base_value)

      # Create an empty list to store the subsetted data frames
    subset_list2 <- list()

    # Identify consecutive groups based on the "height" column
    consecutive_groups2 <- cumsum(c(TRUE, diff(perc_5$height) != step))

    # Iterate over each consecutive group
    for (group in unique(consecutive_groups2)) {
      # Subset the current group
      group_subset2 <- perc_5[consecutive_groups2 == group, ]

      # Add the first and last rows of the subset to the list
      subset_list2[[group]] <- rbind(group_subset2[1, ], group_subset2[nrow(group_subset2), ])
    }

    # Combine the subsetted data frames into a single data frame
    result2 <- do.call(rbind, subset_list2)
    result_unique2 <- unique(result2)

    ################################################3

    crown3b<-rbind.data.frame(crown3,result_unique2)

    crown3b <- crown3b[!duplicated(crown3b), ]
    crown3b<-crown3b[with(crown3b, order(height)), ]

       crown4<-data.frame(t(crown3b))


  ############ ADAPT THE GAPS TO THE CBH (PREFEReNCE THE CBHs because there are some problems with distance calculus)

if (exists("crown3b") && !is.null(crown3b) && length(crown3b) > 0) {
# Select heights in gaps5k that are not in crown3e
gaps5l <- setdiff(gaps5k$height, crown3b$height)
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
  dplyr::mutate(type = ifelse(type == "gap", paste0("gap", seq_along(type[type == "gap"])), type))
merged_gaps2<-data.frame(t(merged_gaps1))
type<-merged_gaps1$type

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
  dplyr::mutate(type = ifelse(type == "cbh", paste0("cbh", seq_along(type[type == "cbh"])), type))
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
   if (!is.null(gaps5m) && nrow(gaps5m) > 0) {
     gaps_height_t <- gaps5m %>%
       dplyr::select(height) %>%
       dplyr::mutate(type = "gap") %>%
       dplyr::rename(V1=height)
   } else {
     gaps_height_t <- tibble(V1 = NA, type = "gap")
   }


   # Check if crown4 exists and contains data
   if (!is.null(crown3b) && nrow(crown3b) > 0) {
     crown_height_t <- crown3b %>%
       dplyr::select(height) %>%
       dplyr::mutate(type = "cbh") %>%
       dplyr::rename(V1=height)
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
    gap_cbh_metrics <- dplyr::arrange(gap_cbh_metrics, treeID1)

    # make copy of treeID
    df$treeID1 <- df$treeID
    treeID1<-unique(factor(df$treeID1))

    # Apply the code
    round_columns <- grep("^(cbh|gap)", names(gap_cbh_metrics), value = TRUE)
    numeric_columns <- sapply(gap_cbh_metrics[round_columns], is.numeric)
    numeric_round_columns <- round_columns[numeric_columns]

    gap_cbh_metrics[numeric_round_columns] <- lapply(gap_cbh_metrics[numeric_round_columns], floor)


}

return(gap_cbh_metrics)
}
