#' Distances between fuel layers
#' @description This function calculates distances (and their heights) between fuel layers as the difference between consecutive gaps and fuel bases
#' (the gap height always must be lower than the fuel base height).
#' @usage get_distance (gap_cbh_metrics,gaps_perc,step=1,min_height=1.5,verbose=TRUE)
#' @param gap_cbh_metrics data frame with gaps (distances) and fuel base heights (output of [get_gaps_fbhs()] function).
#' An object of the class text.
#' @param gaps_perc data frame with Leaf Area Density (LAD) percentiles for each height values (output of [calculate_gaps_perc()] function).
#' An object of the class text.
#' @param step Numeric value for the actual height bin step (in meters).
#' @param min_height Numeric value for the actual minimum base height (in meters).
#' @param verbose Logical, indicating whether to display informational messages (default is TRUE).
#' @return A data frame giving distances (and their heights) between fuel layers in meters.
#' @author Olga Viedma, Carlos Silva, JM Moreno and A.T. Hudak
#'
#'@details
#' # List of tree metrics:
#' \itemize{
#' \item treeID: tree ID with strings and numeric values
#' \item treeID1: tree ID with only numeric values
#' \item cbh - Height of the fuel layers base height (m)
#' \item gap - Height of gaps between consecutive fuel layers (m)
#' \item dist: Distance between consecutive fuel layers (m)
#' \item Hdist - Height of the distance between consecutive fuel layers (m)
#' \item max_height - Maximum height of the tree profile
#' }
#'
#' @examples
#' library(magrittr)
#' library(gdata)
#' library(dplyr)
#'
#' # Before running this example, make sure to run get_gaps_fbhs().
#' if (interactive()) {
#' gap_cbh_metrics <- get_gaps_fbhs()
#' LadderFuelsR::gap_cbh_metrics$treeID <- factor(LadderFuelsR::gap_cbh_metrics$treeID)
#'
#' # Before running this example, make sure to run calculate_gaps_perc().
#' LadderFuelsR::gaps_perc$treeID <- factor(LadderFuelsR::gaps_perc$treeID)
#'
#' trees_name1 <- as.character(gaps_perc$treeID)
#' trees_name2 <- factor(unique(trees_name1))
#'
#' metrics_distance_list <- list()
#'
#' for (i in levels(trees_name2)) {
#'
#' # Filter data for each tree
#' tree1 <- gap_cbh_metrics |> dplyr::filter(treeID == i)
#' tree2 <- gaps_perc |> dplyr::filter(treeID == i)
#' # Get distance metrics for each tree
#' metrics_distance <- get_distance(tree1, tree2, step=1, min_height=1.5)
#' metrics_distance_list[[i]] <- metrics_distance
#' }
#' # Combine the individual data frames
#' distance_metrics <- dplyr::bind_rows(metrics_distance_list)
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
#' @seealso \code{\link{get_gaps_fbhs}}
#' @seealso \code{\link{calculate_gaps_perc}}
#' @export
get_distance <- function (gap_cbh_metrics,gaps_perc, step=1, min_height=1.5, verbose = TRUE) {

  if(min_height==0){
    min_height <-0.5}

    gaps_perc2<-gaps_perc
    df <- gap_cbh_metrics

    gaps_perc2$treeID <- factor(gaps_perc2$treeID)
    df$treeID <- factor(df$treeID)

    if (verbose) {
      message("Unique treeIDs:", paste(unique(gaps_perc2$treeID), collapse = ", "))
    }

    treeID<-factor(df$treeID)
    treeID1<-factor(df$treeID1)

    df <- df %>%
    dplyr::mutate_at(
      vars(-treeID),  # Exclude the 'treeID' column
      as.numeric
    )

  df1 <- df[, !colSums(is.na(df)) > 0]
  # Select only numeric columns
  df1_numeric <- df1 %>% dplyr::select_if(is.numeric)


  # Assuming that columns starting with "gap" or "cbh" are the ones you want to keep
  columns_to_keep <- names(df1_numeric)[grepl("^gap\\d+$|^cbh\\d+$", names(df1_numeric))]
    # Subset the data frame
  kk_copy <- df1_numeric[, columns_to_keep]

  # Sort the column names based on their values in the first row
  sorted_columns <- names(kk_copy)[order(unlist(kk_copy))]

  # Reorder the columns in ascending order
  kk_copy <- kk_copy[, sorted_columns]

  cbh_cols <- which(grepl("cbh", colnames(kk_copy)))
  gap_cols <- which(grepl("gap", colnames(kk_copy)))


  ### difference between consecutive gaps and any cbh starting from the first gap:
  # calculate the difference between the first consecutive 'gap' column in a series and the next 'cbh' column,
  #while the 'Hdist' variable should hold the value of the last consecutive 'gap' column

  if (is.null(kk_copy[, grep("^gap", colnames(kk_copy))]) || length(gap_cols) == 0 || length(cbh_cols) == 0) {
    distance_data <- data.frame(c(NA))
    names(distance_data) <- "dist_0"

  } else if (length(cbh_cols) > 1) {
    distance_data <- data.frame()  # initialize distance_data before the loop
    i <- 1
    first_gap_col <- NULL
    last_gap_col <- NULL
    cbh_col <- NULL
    diff_col_names <- list()

    while (i <= ncol(kk_copy)) {
      col_name <- names(kk_copy)[i]

      # Skip past consecutive 'gap' columns and store only the last one
      if (startsWith(col_name, "gap")) {
        if (is.null(first_gap_col)) {
          first_gap_col <- col_name
        }
        while (i <= ncol(kk_copy) && startsWith(names(kk_copy)[i], "gap")) {
          last_gap_col <- names(kk_copy)[i]
          i <- i + 1
        }
      }

      # Now we're at the first 'cbh' after a series of 'gap's or at the end of the dataframe
      if (i <= ncol(kk_copy) && startsWith(names(kk_copy)[i], "cbh")) {
        cbh_col <- names(kk_copy)[i]

        if (!is.null(first_gap_col) && !is.null(last_gap_col)) {
          diff_col_name <- paste0(first_gap_col, "_", cbh_col)
          diff_col_names[[diff_col_name]] <- abs(kk_copy[[cbh_col]] - kk_copy[[first_gap_col]])

          Hdist <- kk_copy[[last_gap_col]]
          Hdist_name <- paste0("Hdist_", last_gap_col)
          diff_col_names[[Hdist_name]] <- Hdist  # add Hdist to diff_col_names
        }

        first_gap_col <- NULL
        last_gap_col <- NULL
        cbh_col <- NULL
      }

      # Break out of the loop if there are no more columns
      if (i > ncol(kk_copy)) {
        break
      }

      i <- i + 1
    }

    if (length(diff_col_names) == 0) {
      distance_data <- data.frame(Hdist = kk_copy[[last_gap_col]])
    } else {
      distance_data <- as.data.frame(diff_col_names)  # construct data frame after the loop
    }
  }

 #####################################

####################################

  if(length(gap_cols) >= 1 & length(cbh_cols) == 0 && (!exists("distance_data") || is.null(distance_data) || ncol(distance_data) == 0 && nrow(distance_data) == 0 ||
                                                       all(is.na(distance_data))))  {

    # Identify the first and last indices of the consecutive "gap" columns
    first_gap_index <- min(gap_cols)
    last_gap_index <- max(gap_cols)

    # Subset the consecutive "gap" columns
    gap_subset <- kk_copy[, gap_cols, drop = FALSE]

    # Calculate the difference between the last and the previous consecutive "gap" columns
    last_gap_value <- gap_subset[, last_gap_index]
    previous_gap_value <- gap_subset[, last_gap_index - 1]
    gap_difference <- last_gap_value - previous_gap_value

    distance_data <- as.data.frame(gap_difference)
    Hdist1<-last_gap_value
    distance_data<-cbind.data.frame(distance_data,Hdist1)

  }

  #### IF THERE ARE GAPS BEFORE and AFTER A CBH VALUE: CALCULATE THE GAP BELOW THE CBH #################


  if (length(gap_cols) > 1 && length(cbh_cols) >= 1 ) {

    split_by_prefix_consecutive <- function(df) {
      prefixes <- sapply(strsplit(names(df), "[0-9]+"), `[`, 1) # Split names by numeric part

      # We'll use rle() to find where the prefixes change:
      r <- rle(prefixes)

      # Then, we'll use the lengths from r to split the data frame's columns:
      df_list <- split.default(df, rep(seq_along(r$lengths), r$lengths))

      # Name the list elements based on the prefix in each data frame:
      names(df_list) <- r$values

      return(df_list)
    }

    df_list <- split_by_prefix_consecutive(kk_copy)

    max_cbh_col <- max(kk_copy[, cbh_cols])
    min_cbh_col <- min(kk_copy[, cbh_cols])

    # Iterate over each dataframe in df_list
    for (df_name in names(df_list)) {

      # Check if the dataframe has "gap" columns
      if (startsWith(df_name, "gap")) {

        gap_df <- df_list[[df_name]]

        # Check if there are at least 2 columns in gap_df
        if (ncol(gap_df) >= 2) {

          # Check if the maximum cbh column value is greater than the last column of gap_df
          if (max_cbh_col > gap_df[,ncol(gap_df)] && min_cbh_col > gap_df[,ncol(gap_df)]) {

            gap_difference1 <- min_cbh_col - gap_df[,1]
            distance_data1 <- as.data.frame(gap_difference1)
            Hdist1 <- gap_df[,2]
            distance_data1 <- cbind.data.frame(distance_data1,Hdist1)
          }
        }
      }
    }
  }

  ##################################

  if ((exists("distance_data") && !is.null(distance_data) && ncol(distance_data) != 0 && nrow(distance_data) != 0) &&
      (exists("distance_data1") && !is.null(distance_data1) && ncol(distance_data1) != 0 && nrow(distance_data1) != 0)) {

    # Get column names with prefixes in both dataframes
    gap_cols_distance_data <- names(distance_data)[startsWith(names(distance_data), "gap")]
    Hdist_cols_distance_data <- names(distance_data)[startsWith(names(distance_data), "Hdist")]

    gap_cols_distance_data1 <- names(distance_data1)[startsWith(names(distance_data1), "gap")]
    Hdist_cols_distance_data1 <- names(distance_data1)[startsWith(names(distance_data1), "Hdist")]

    # Initialize the equality check results
    gap_equal <- TRUE
    Hdist_equal <- TRUE

    # Check equality for each column
    for (gap_col in seq_along(gap_cols_distance_data)) {
      if(gap_col <= length(gap_cols_distance_data1)){
        gap_equal <- gap_equal & (distance_data[, gap_cols_distance_data[gap_col]] == distance_data1[, gap_cols_distance_data1[gap_col]])
      }
    }

    for (Hdist_col in seq_along(Hdist_cols_distance_data)) {
      if(Hdist_col <= length(Hdist_cols_distance_data1)){
        Hdist_equal <- Hdist_equal & (distance_data[, Hdist_cols_distance_data[Hdist_col]] == distance_data1[, Hdist_cols_distance_data1[Hdist_col]])
      }
    }

    # If any of the columns are not equal, bind the dataframes
    if (!gap_equal | !Hdist_equal) {
      distance_data <- cbind(distance_data, distance_data1)
    }
  }

  #### IF THERE ARE GAPS BEFORE and AFTER A CBH VALUE: CALCULATE THE GAP ABOVE THE CBH #################

  # Find the last set of consecutive "gap" columns

  # Find columns that start with "gap" or "cbh"
  gap_cols <- grep("^gap", colnames(kk_copy))
  cbh_cols <- grep("^cbh", colnames(kk_copy))

  # Check the number of gap columns is greater than 1 and number of cbh columns equals to 1
  if(length(gap_cols) > 1 & length(cbh_cols) == 1 ) {

    # Check the number of gap columns is greater than 1
    if (length(gap_cols) > 1) {
      # Get the column names
      col_names <- colnames(kk_copy)

      # Identify columns starting with "gap" and "cbh"
      gap_cols <- grep("^gap", col_names)
      cbh_cols <- grep("^cbh", col_names)

      # Determine the "gap" column sets based on the presence of a "cbh" column
      set_breaks <- c(0, cbh_cols, length(col_names) + 1)
      gap_sets <- list()

      for (i in 1:(length(set_breaks) - 1)) {
        set_start <- set_breaks[i] + 1
        set_end <- set_breaks[i + 1] - 1
        if (set_end >= set_start) {
          gap_sets[[i]] <- col_names[set_start:set_end]
        }
      }

      # If there's more than one set of consecutive "gap" columns
      if (length(gap_sets) > 1) {
        # Get the last set of "gap" columns
        last_set <- gap_sets[[length(gap_sets)]]

        # Calculate the difference between the last and first column in this set
        last_cons_gap <- kk_copy[, last_set, drop = FALSE]
        gap_difference <- last_cons_gap[, length(last_set)] - last_cons_gap[, 1]
        distance_data2 <- data.frame(dist = gap_difference, Hdist = last_cons_gap[, length(last_set)])
      } else {
        # If there's only one set of consecutive "gap" columns, assign NA values
        distance_data2 <- data.frame(dist = NA, Hdist = NA)
      }
    }

  } else if (length(gap_cols) > 1) {  # Other case where gap_cols > 1

    # Subset the first and last value of the "gap" columns
    subset_values <- kk_copy[, gap_cols]
    first_value <- subset_values[1]
    last_value <- subset_values[length(subset_values)]

    gap_difference1 <- last_value - first_value
    distance_data3 <- as.data.frame(gap_difference1)
    Hdist4<-last_value
  }

    # Add Hdist4 column to the data frame if it exists and is not empty
    if (exists("Hdist4") && exists("distance_data3") && !is.null(Hdist4) && ncol(Hdist4) != 0 && nrow(Hdist4) != 0 && any(!is.na(Hdist4))) {
      distance_data3 <- cbind(distance_data3,Hdist4)
      colnames(distance_data3)<-c("dist", "Hdist")
    }


  ###############################


  if (!exists("distance_data")) {
    if (exists("distance_data1")) {
      if(!is.null(distance_data1) || ncol(distance_data1) != 0 || nrow(distance_data1) != 0 || any(!is.na(distance_data1))) {
        if(!exists("distance_data2") || is.null(distance_data2) || any(is.na(distance_data2))) {
          if (!exists("distance_data3") || is.null(distance_data3) || any(is.na(distance_data3))) {

            distance_data <-data.frame(distance_data1)
          }}}}}


  if (!exists("distance_data")) {
    if (exists("distance_data2")) {
      if(!is.null(distance_data2) || ncol(distance_data2) != 0 || nrow(distance_data2) != 0 || any(!is.na(distance_data2))) {
        if(!exists("distance_data1") || is.null(distance_data1) || any(is.na(distance_data1))) {
          if (!exists("distance_data3") || is.null(distance_data3) || any(is.na(distance_data3))) {


            distance_data <-data.frame(distance_data2)
          }}}}}


  if (!exists("distance_data")) {
    if (exists("distance_data3")) {
      if(!is.null(distance_data3) || ncol(distance_data3) != 0 || nrow(distance_data3) != 0 || any(!is.na(distance_data3))) {
        if(!exists("distance_data1") || is.null(distance_data1) || any(is.na(distance_data1))) {
          if (!exists("distance_data2") || is.null(distance_data2) || any(is.na(distance_data2))) {

            distance_data <-data.frame(distance_data3)
          }}}}}


  if (!exists("distance_data")) {
    if (exists("distance_data1")) {
      if(!is.null(distance_data1) || ncol(distance_data1) != 0 || nrow(distance_data1) != 0 || any(!is.na(distance_data1))) {
        if(exists("distance_data2")) {
          if(!is.null(distance_data2) || any(!is.na(distance_data2))) {
            if (!exists("distance_data3") || is.null(distance_data3) || any(is.na(distance_data3))) {

              distance_data <-data.frame(distance_data1, distance_data2)
            }}}}}}


  if (!exists("distance_data")) {
    if (exists("distance_data1")) {
      if(!is.null(distance_data1) || ncol(distance_data1) != 0 || nrow(distance_data1) != 0 || any(!is.na(distance_data1))) {
        if(exists("distance_data3")) {
          if(!is.null(distance_data3) || any(!is.na(distance_data3))) {
            if (!exists("distance_data2") || is.null(distance_data2) || any(is.na(distance_data2))) {

              distance_data <-data.frame(distance_data1, distance_data3)
            }}}}}}


  #################################
  # Initialize `Hdist_equal` as TRUE
  Hdist_equal <- TRUE

  # Check if `distance_data` exists, and if it does not, initialize it
  if (!exists("distance_data")) {
    distance_data <- data.frame()
  }

    condition1 <- is.null(distance_data) || ncol(distance_data) == 0 || nrow(distance_data) == 0 || any(is.na(distance_data))
    condition2 <- exists("distance_data1") && !is.null(distance_data1) && ncol(distance_data1) != 0 && nrow(distance_data1) != 0 && any(!is.na(distance_data1))
    condition3 <- exists("distance_data2") && !is.null(distance_data2) && ncol(distance_data2) != 0 && nrow(distance_data2) != 0 && any(!is.na(distance_data2))

    if (condition1 && condition2 && condition3) {

    # Get the column names with 'Hdist_' prefix from both data frames
    Hdist_cols_data1 <- grep("^Hdist", names(distance_data1), value = TRUE)
    Hdist_cols_data2 <- grep("^Hdist", names(distance_data2), value = TRUE)

    # Loop through each 'Hdist_' column from the first data frame and compare it with each 'Hdist_' column from the second data frame
    for (Hdist_col_data1 in Hdist_cols_data1) {
      for (Hdist_col_data2 in Hdist_cols_data2) {
        # Check if the unique values in the current pair of columns are equal
        Hdist_equal <- Hdist_equal & all(distance_data1[, Hdist_col_data1] == distance_data2[, Hdist_col_data2])
      }
    }

    # If any of the pairs of columns are not equal, bind the data frames
    # Initialize `Hdist_equal` as TRUE
    Hdist_equal <- TRUE

    # Check if `distance_data` exists, and if it does not, initialize it
    if (!exists("distance_data")) {
      distance_data <- data.frame()
    }

    condition1 <- is.null(distance_data) || ncol(distance_data) == 0 || nrow(distance_data) == 0 || any(is.na(distance_data))
    condition2 <- exists("distance_data1") && !is.null(distance_data1) && ncol(distance_data1) != 0 && nrow(distance_data1) != 0 && any(!is.na(distance_data1))
    condition3 <- exists("distance_data2") && !is.null(distance_data2) && ncol(distance_data2) != 0 && nrow(distance_data2) != 0 && any(!is.na(distance_data2))

    if (condition1 && condition2 && condition3) {


      # Get the column names with 'Hdist_' prefix from both data frames
      Hdist_cols_data1 <- grep("^Hdist", names(distance_data1), value = TRUE)
      Hdist_cols_data2 <- grep("^Hdist", names(distance_data2), value = TRUE)

      # Loop through each 'Hdist_' column from the first data frame and compare it with each 'Hdist_' column from the second data frame
      for (Hdist_col_data1 in Hdist_cols_data1) {
        for (Hdist_col_data2 in Hdist_cols_data2) {
          # Check if the unique values in the current pair of columns are equal
          Hdist_equal <- Hdist_equal & all(distance_data1[, Hdist_col_data1] == distance_data2[, Hdist_col_data2])

          # If Hdist values are equal, remove the corresponding Hdist and its immediately previous gap column from distance_data1
          if (Hdist_equal) {
            # Get the index of the current Hdist column in distance_data1
            Hdist_index <- which(names(distance_data1) == Hdist_col_data1)

            # If the index of the current Hdist column is greater than 1
            if (Hdist_index > 1) {
              # Get the name of the immediately previous column
              prev_gap_col <- names(distance_data1)[Hdist_index - 1]

              # If the name of the immediately previous column starts with 'gap'
              if (startsWith(prev_gap_col, "gap")) {
                # Remove the current Hdist and its immediately previous gap column
                distance_data1 <- distance_data1[, !(names(distance_data1) %in% c(prev_gap_col, Hdist_col_data1))]
              }
            }

            # Reset Hdist_equal to TRUE for the next comparison
            Hdist_equal <- TRUE
          }
        }
      }

      # Bind the data frames
      distance_data <- cbind.data.frame(distance_data1, distance_data2)
    }
  }

  #############################

  # Initialize `Hdist_equal` as TRUE
  Hdist_equal <- TRUE

  # Check if `distance_data` exists, and if it does not, initialize it
  if (!exists("distance_data")) {
    distance_data <- data.frame()
  }

    condition1 <- is.null(distance_data) || ncol(distance_data) == 0 || nrow(distance_data) == 0 || any(is.na(distance_data))
    condition2 <- exists("distance_data1") && !is.null(distance_data1) && ncol(distance_data1) != 0 && nrow(distance_data1) != 0 && any(!is.na(distance_data1))
    condition3 <- exists("distance_data3") && !is.null(distance_data3) && ncol(distance_data3) != 0 && nrow(distance_data3) != 0 && any(!is.na(distance_data3))

    if (condition1 && condition2 && condition3) {


    # Get the column names with 'Hdist_' prefix from both data frames
    Hdist_cols_data1 <- grep("^Hdist", names(distance_data1), value = TRUE)
    Hdist_cols_data2 <- grep("^Hdist", names(distance_data3), value = TRUE)

    # Loop through each 'Hdist_' column from the first data frame and compare it with each 'Hdist_' column from the second data frame
    for (Hdist_col_data1 in Hdist_cols_data1) {
      for (Hdist_col_data2 in Hdist_cols_data2) {
        # Check if the unique values in the current pair of columns are equal
        Hdist_equal <- Hdist_equal & all(distance_data1[, Hdist_col_data1] == distance_data3[, Hdist_col_data2])
      }
    }

    # If any of the pairs of columns are not equal, bind the data frames
    # Initialize `Hdist_equal` as TRUE
    Hdist_equal <- TRUE

    # Check if `distance_data` exists, and if it does not, initialize it
    if (!exists("distance_data")) {
      distance_data <- data.frame()
    }


      condition1 <- is.null(distance_data) || ncol(distance_data) == 0 || nrow(distance_data) == 0 || any(is.na(distance_data))
      condition2 <- exists("distance_data1") && !is.null(distance_data1) && ncol(distance_data1) != 0 && nrow(distance_data1) != 0 && any(!is.na(distance_data1))
      condition3 <- exists("distance_data3") && !is.null(distance_data3) && ncol(distance_data3) != 0 && nrow(distance_data3) != 0 && any(!is.na(distance_data3))

      if (condition1 && condition2 && condition3) {

      # Get the column names with 'Hdist_' prefix from both data frames
      Hdist_cols_data1 <- grep("^Hdist", names(distance_data1), value = TRUE)
      Hdist_cols_data2 <- grep("^Hdist", names(distance_data3), value = TRUE)

      # Loop through each 'Hdist_' column from the first data frame and compare it with each 'Hdist_' column from the second data frame
      for (Hdist_col_data1 in Hdist_cols_data1) {
        for (Hdist_col_data2 in Hdist_cols_data2) {
          # Check if the unique values in the current pair of columns are equal
          Hdist_equal <- Hdist_equal & all(distance_data1[, Hdist_col_data1] == distance_data3[, Hdist_col_data2])

          # If Hdist values are equal, remove the corresponding Hdist and its immediately previous gap column from distance_data1
          if (Hdist_equal) {
            # Get the index of the current Hdist column in distance_data1
            Hdist_index <- which(names(distance_data1) == Hdist_col_data1)

            # If the index of the current Hdist column is greater than 1
            if (Hdist_index > 1) {
              # Get the name of the immediately previous column
              prev_gap_col <- names(distance_data1)[Hdist_index - 1]

              # If the name of the immediately previous column starts with 'gap'
              if (startsWith(prev_gap_col, "gap")) {
                # Remove the current Hdist and its immediately previous gap column
                distance_data1 <- distance_data1[, !(names(distance_data1) %in% c(prev_gap_col, Hdist_col_data1))]
              }
            }

            # Reset Hdist_equal to TRUE for the next comparison
            Hdist_equal <- TRUE
          }
        }
      }

      # Bind the data frames
      distance_data <- cbind.data.frame(distance_data1, distance_data3)
    }
  }

  #################################

  if (exists("distance_data") && nrow(distance_data) != 0 && any(!is.na(distance_data))) {

    if (exists("distance_data1")) {
      if (length(distance_data1) > 0 && any(!is.na(distance_data1))) {

        if (exists("distance_data2")) {
          if (length(distance_data2) > 0 && any(!is.na(distance_data2))) {

            if (exists("distance_data3")) {
              if (length(distance_data3) > 0 && any(!is.na(distance_data3))) {

                if(ncol(distance_data)==1) {

                  if(distance_data2[,2] == distance_data3[,2] && distance_data1[,2] == distance_data2[,2] && distance_data1[,2] == distance_data3[,2]) {

                    distance_data <- distance_data1

                  } else {
                    distance_data <- distance_data
                  }
                }
              }}}}}}}


  if ((exists("distance_data") && nrow(distance_data) != 0 && any(!is.na(distance_data)))) {

    if (!exists("distance_data1")) {

      if (exists("distance_data2")) {
        if (length(distance_data2) > 0 && any(!is.na(distance_data2))) {

          if (exists("distance_data3")) {
            if (length(distance_data3) > 0 && any(!is.na(distance_data3))) {

              if(ncol(distance_data)==1) {

                if(distance_data2[,2] == distance_data3[,2]) {

                  distance_data <- distance_data2

                } else {
                  distance_data <- distance_data
                }
              }
            }}}}}}


  if ((exists("distance_data") && nrow(distance_data) != 0 && any(!is.na(distance_data)))) {

    if (!exists("distance_data1")) {

      if (!exists("distance_data2")) {

        if (exists("distance_data3") && length(distance_data3) > 0 && any(!is.na(distance_data3))) {

          if(ncol(distance_data)==1) {

            distance_data <- distance_data3

          } else {
            distance_data <- distance_data
          }
        }
      }}}


  if ((exists("distance_data") && nrow(distance_data) != 0 && any(!is.na(distance_data)))) {

    if  (exists("distance_data1") && length(distance_data1) != 0 && any(!is.na(distance_data1))) {
      if (exists("distance_data2") && length(distance_data2) != 0 && any(!is.na(distance_data2))) {
        if (!exists("distance_data3")) {

          if(ncol(distance_data)==1) {

            if(distance_data1[,2] != distance_data2[,2]) {

              distance_data <- distance_data1

            } else {
              distance_data <- distance_data
            }
          }
        }}}}


  if ((exists("distance_data") && nrow(distance_data) != 0 && any(!is.na(distance_data)))) {

    if (exists("distance_data1") && length(distance_data1) == 0 && any(is.na(distance_data1))) {
      if (!exists("distance_data2")) {
        if (exists("distance_data3") && length(distance_data3) != 0 && any(!is.na(distance_data3))) {

          if(ncol(distance_data)==1) {

            if(distance_data1[,2] != distance_data3[,2]) {

              distance_data <- distance_data1

            } else {
              distance_data <- distance_data
            }
          }
        }}}
    }


  if (length(gap_cols) > 1 && length(cbh_cols) > 1 && ((exists("distance_data") || !is.null(distance_data) || (ncol(distance_data) != 0 &&
  nrow(distance_data) != 0)) || any(!is.na(distance_data)))) {

    height<-gaps_perc2$height
    percent2a <- gaps_perc2 %>% dplyr::filter(height < min(kk_copy[, cbh_cols]))

    # Check the conditions and remove the first column if they are met
    if (nrow(percent2a) != 0 && any(!is.na(percent2a)) && any(percent2a$percentil > 25)) {
      distance_data <- distance_data %>%
        dplyr::select(-1)  # Remove the first column
    }
  }


  #### IF THERE ARE CONSECUTIVE GAPS AFTER CBH COLUMNS (ABOVE DISTANCES) #################

  gap_cols <- grep("^gap", colnames(kk_copy))
  cbh_cols <- grep("^cbh", colnames(kk_copy))

  # remove missing values in gap_cols or cbh_cols
  gap_cols <- gap_cols[!is.na(kk_copy[,gap_cols])]
  cbh_cols <- cbh_cols[!is.na(kk_copy[,cbh_cols])]

  ###################################

  if (length(gap_cols) > 1 && length(cbh_cols) > 1 &&
      ((exists("distance_data") || !is.null(distance_data) ||
        (ncol(distance_data) != 0 && nrow(distance_data) != 0)) ||
       any(!is.na(distance_data)))) {

    if (length(gap_cols) > 1 && any(diff(gap_cols) > 1)) {

      # Get the positions of 'gap' and 'cbh' columns
      gap_positions <- grep("^gap", names(kk_copy))
      cbh_positions <- grep("^cbh", names(kk_copy))

      # Find the last position of 'cbh' column
      last_cbh_position <- max(cbh_positions)

      # Check if there's any 'gap' column after the last 'cbh' column
      if (any(gap_positions > last_cbh_position)) {

        # Get the first 'gap' column that appears after the `last_cbh_position`
        first_gap_after_cbh <- min(gap_positions[gap_positions > last_cbh_position])

        # Get the set of consecutive 'gap' columns from `first_gap_after_cbh` to the last column
        last_set_cols <- names(kk_copy)[first_gap_after_cbh:length(kk_copy)]

        # Subset the last set of consecutive "gap" columns
        last_cons_gap <- kk_copy[, last_set_cols, drop = FALSE]

        # Retrieve the first and last gap column of the subsetted data frame
        first_gap_column <- last_cons_gap[, 1, drop = FALSE]
        last_gap_column <- last_cons_gap[, length(last_set_cols), drop = FALSE]

        max_cbh_col<-max(kk_copy[, cbh_cols])

        if (max_cbh_col < first_gap_column) {
          # Calculate the difference between the last and first gap columns
          gap_difference2 <- last_gap_column - first_gap_column

          distance_data4 <- as.data.frame(gap_difference2)

          Hdist5<-last_gap_column


          # Get the column index of the column that starts with "gap"
          gap_col1<- grep("^gap", colnames(Hdist5))

          # Rename the column
          if (length(gap_col1) > 0) {
            for (i in gap_col1) {
              colnames(Hdist5)[i] <- paste0("Hdist_", colnames(Hdist5)[i])

            }
          }

          distance_data4<-cbind.data.frame(distance_data4,Hdist5)
        }
      }
    }
  }



    condition1 <- exists("distance_data") && !is.null(distance_data) && ncol(distance_data) != 0 && nrow(distance_data) != 0 && any(!is.na(distance_data))
    condition2 <- exists("distance_data4") && !is.null(distance_data4) && ncol(distance_data4) != 0 && nrow(distance_data4) != 0 && any(!is.na(distance_data4))

    if (condition1 && condition2) {


    # Get the column names with 'Hdist_' prefix from both data frames
    Hdist_cols_data <- grep("^Hdist", names(distance_data), value = TRUE)
    Hdist_cols_data4 <- grep("^Hdist", names(distance_data4), value = TRUE)

    # Initialize the equality check results
    Hdist_equal <- TRUE

    # Loop through each 'Hdist_' column from the first data frame and compare it with each 'Hdist_' column from the second data frame
    for (Hdist_col_data in Hdist_cols_data) {
      for (Hdist_col_data4 in Hdist_cols_data4) {
        # Check if the unique values in the current pair of columns are equal
        Hdist_equal <- Hdist_equal & all(distance_data[, Hdist_col_data] == distance_data4[, Hdist_col_data4])
      }
    }

    # If any of the pairs of columns are not equal, bind the data frames
    if (!Hdist_equal) {
      distance_data <- cbind(distance_data, distance_data4)
    }
  }


  if(exists("distance_data") && !is.null(distance_data) && ncol(distance_data) != 0 && nrow(distance_data) != 0 &&  any(!is.na(distance_data))
     && all(distance_data==0)) {

    distance_data <-data.frame(c(NA))
    names(distance_data)= "dist_0"
  }

  if (length(gap_cols) > 1 && length(cbh_cols) > 1 && ((!exists("distance_data") || is.null(distance_data) || (ncol(distance_data) == 0 && nrow(distance_data) == 0)) ||
                                                       all(is.na(distance_data)))) {

    diff2 <- max(kk_copy[, gap_cols]) - min(kk_copy[, gap_cols])
    Hdist6<-max(kk_copy[, gap_cols])
    distance_data <- cbind.data.frame(diff2,Hdist6)
  }



  if(exists("distance_data") && exists("distance_data1") && all(is.na(distance_data)) && all(!is.na(distance_data1)) && (!exists("distance_data3") && !exists("distance_data4"))){
    Hdist_cols <- grep("^(Hdist)", names(distance_data1), value = TRUE)
    distance_data <- distance_data1[Hdist_cols]
  }

  if (length(distance_data) == 0) {
    distance_data <- data.frame(distance_data = 0)
    names(distance_data)="dist0"
  }

 ##################################
    hdist_columns <- grep("^Hdist", names(distance_data), value = TRUE)
    max_value <- max(unlist(distance_data[hdist_columns]), na.rm = TRUE)

    cbh_cols <- grepl("^cbh", names(kk_copy))
    max_cbh_value <- max(kk_copy[, cbh_cols])

    # Check if all values in Hdist columns are greater than max_cbh_value
    all_greater <- all(distance_data[, hdist_columns] > max_cbh_value)

    if (all_greater) {
      distance_data<-NA
    }

  #######################################################

    # Check if distance_data is not empty
    if (length(distance_data) != 0 && any(!is.na(distance_data))) {

     # Function to remove columns if Hdist > last cbh column value
    remove_hdist_if_greater_than_cbh <- function(distance_data, kk_copy) {

      # Find the value of the last cbh column in kk_copy
      cbh_cols <- grep("^cbh", names(kk_copy), value = TRUE)
      last_cbh_value <- kk_copy[[cbh_cols[length(cbh_cols)]]]

      # Get column names in distance_data
      col_names <- colnames(distance_data)

      # Find columns to remove
      cols_to_remove <- c()

      for (col in col_names) {
        if (startsWith(col, "Hdist")) {
          # Get the numeric suffix
          suffix <- gsub("\\D", "", col)
          dist_col <- paste0("gap", suffix)

          if (any(distance_data[[col]] > last_cbh_value)) {
            cols_to_remove <- c(cols_to_remove, col, dist_col)
          }
        }
      }

      # Remove columns
      distance_data <- distance_data[, !colnames(distance_data) %in% cols_to_remove, drop = FALSE]

      return(distance_data)
    }

    # Apply the function
    distance_data_prueba <- remove_hdist_if_greater_than_cbh(distance_data, kk_copy)

    distance_data<-distance_data_prueba

  ##############################################################

    pair_cbh_Hdist <- function(df) {
      # Identify cbh and Hdist columns
      cbh_cols <- grep("^cbh", names(df), value = TRUE)
      Hdist_cols <- grep("^Hdist", names(df), value = TRUE)

      # Debug ##print initial column lists
      # cat("Initial cbh columns:", cbh_cols, "\n")
      # cat("Initial Hdist columns:", Hdist_cols, "\n")

      # Initialize lists to store paired columns
      paired_cbh_cols <- list()
      paired_hdist_cols <- list()

      # Handle Hdist columns with values greater than 0
      df_with_values <- df
      for (hdist_col in Hdist_cols) {
        hdist_values <- df_with_values[[hdist_col]]
        if (any(hdist_values > 0)) {
          #  cat("Processing", hdist_col, "with values:", hdist_values[hdist_values > 0], "\n")

          matched <- FALSE
          for (cbh_col in cbh_cols) {
            cbh_values <- df_with_values[[cbh_col]]
            for (i in which(hdist_values > 0)) {
              hdist_value <- hdist_values[i]
              if (cbh_values[i] >= hdist_value) {
                # Rename Hdist column with the suffix of matched cbh column
                suffix <- gsub("^cbh", "", cbh_col)
                new_hdist_col <- paste0("Hdist", suffix)

                # Check if the new Hdist column name already exists
                if (new_hdist_col %in% names(df_with_values)) {
                  # Rename the current Hdist column to HdistX (incremental number)
                  new_hdist_col <- make.unique(new_hdist_col, sep = "_")
                }

                # cat("Renaming", hdist_col, "to", new_hdist_col, "\n")

                names(df_with_values)[names(df_with_values) == hdist_col] <- new_hdist_col

                # Store paired columns
                paired_cbh_cols[[new_hdist_col]] <- cbh_col
                paired_hdist_cols[[new_hdist_col]] <- new_hdist_col

                # Remove from remaining Hdist columns
                Hdist_cols <- setdiff(Hdist_cols, hdist_col)

                # Mark as matched
                matched <- TRUE
                break
              }
            }
            if (matched) break
          }

          if (!matched) {
            # If no match found, rename with the original suffix
            suffix <- gsub("^Hdist", "", hdist_col)
            new_hdist_col <- paste0("Hdist", suffix)

            # Check if the new Hdist column name already exists
            if (new_hdist_col %in% names(df_with_values)) {
              # Rename the current Hdist column to HdistX (incremental number)
              new_hdist_col <- make.unique(new_hdist_col, sep = "_")
            }

            #  cat("Renaming", hdist_col, "to", new_hdist_col, "\n")

            names(df_with_values)[names(df_with_values) == hdist_col] <- new_hdist_col

            # Store paired columns
            paired_hdist_cols[[new_hdist_col]] <- new_hdist_col

            # Remove from remaining Hdist columns
            Hdist_cols <- setdiff(Hdist_cols, hdist_col)
          }
        }
      }

      # Debug ##print paired columns
      #cat("Paired cbh columns:", unlist(paired_cbh_cols), "\n")
      #cat("Paired Hdist columns:", unlist(paired_hdist_cols), "\n")

      # Create a new dataframe for remaining Hdist columns with values equal to 0
      remaining_cbh_cols <- setdiff(cbh_cols, unlist(paired_cbh_cols))
      df_remaining <- df_with_values[, c(remaining_cbh_cols, Hdist_cols)]

      if(length(df_remaining > 0)) {

        # Rename remaining Hdist columns sequentially with the remaining cbh columns
        for (i in seq_along(Hdist_cols)) {
          hdist_col <- Hdist_cols[i]
          if (i <= length(remaining_cbh_cols)) {
            cbh_col <- remaining_cbh_cols[i]
            # Rename Hdist column with the suffix of remaining cbh column
            suffix <- gsub("^cbh", "", cbh_col)
            new_hdist_col <- paste0("Hdist", suffix)

            # Check if the new Hdist column name already exists
            if (new_hdist_col %in% names(df_remaining)) {
              # Rename the current Hdist column to HdistX (incremental number)
              new_hdist_col <- make.unique(new_hdist_col, sep = "_")
            }

            names(df_remaining)[names(df_remaining) == hdist_col] <- new_hdist_col

            #cat("Renaming", hdist_col, "to", new_hdist_col, "\n")
          }
        }

        # ##print remaining columns for debugging
        #cat("Remaining Hdist columns with values 0:", Hdist_cols, "\n")
        #cat("Remaining cbh columns:", remaining_cbh_cols, "\n")

        # Return the modified dataframe with values and the cleaned dataframe
        return(list(df_with_values, df_remaining))

      } else {

        return(list(df_with_values))
      }
        # Return the modified dataframe with values and the cleaned dataframe
      return(list(df_with_values, df_remaining))
    }


    treeID<-unique(factor(df$treeID))
    treeID1<-unique(factor(df$treeID1))

    max_height<-data.frame(df$max_height)
    names(max_height)<-"max_height"


    ##################################################3

    distance_data1<-distance_data

    # Extract the names of the columns
    col_names <- names(distance_data1)

    # Identify columns starting with "gap" and containing "cbh"
    gap_cbh_cols <- grep("^gap.*cbh", colnames(distance_data1))

    # Check if any columns match the pattern
    if (length(gap_cbh_cols) > 0) {
      # Extract the numeric suffix from the "cbh" part and create new names
      new_col_names <- col_names
      for (col in gap_cbh_cols) {
        num_suffix <- gsub("^gap.*_cbh(\\d+)$", "\\1", col_names[col])
        new_col_names[col] <- paste0("dist", num_suffix)
      }

      # Rename the columns in the dataframe
      names(distance_data1) <- new_col_names
    }
    distance_data2 <-distance_data1


    # Identify columns named "gap" with only a numeric suffix
    gap_numeric_cols <- grep("^gap\\d+$", colnames(distance_data2))
    col_names <- names(distance_data2)

    # Check if any columns match the pattern
    if (length(gap_numeric_cols) > 0) {
      # Extract the numeric suffix from the "gap" part and create new names
      new_col_names <- col_names
      for (col in gap_numeric_cols) {
        num_suffix <- gsub("^gap(\\d+)$", "\\1", col_names[col])
        new_col_names[col] <- paste0("dist", num_suffix)
      }

      # Rename the columns in the dataframe
      names(distance_data2) <- new_col_names
    }

      distance_data1 <- distance_data2

    ###############################################

    # Identify columns starting with "gap_"
    gap_cols <- grep("^gap_", colnames(distance_data1))

    # Check if any columns match the pattern
    if (length(gap_cols) > 0) {
      # Extract the numeric suffix from the "gap_" part and create new names
      new_col_names <- colnames(distance_data1)
      for (col in gap_cols) {
        num_suffix <- gsub("^gap_(\\D*)(\\d+)$", "\\2", colnames(distance_data1)[col])
        new_col_names[col] <- paste0("dist", num_suffix)
      }

      # Rename the columns in the dataframe
      names(distance_data1) <- new_col_names

      # Remove original "gap_" columns
      distance_data1 <- distance_data1[, !grepl("^gap_", colnames(distance_data1))]
    }

    ###############################################

      # Check if there are any Hdist_gap columns
      Hdist_gap_cols <- grep("^Hdist_gap", colnames(distance_data1))

      # If Hdist_gap columns are found, proceed with renaming
      if (length(Hdist_gap_cols) > 0) {

        # Function to rename Hdist_gap columns
        rename_hdist_gap_columns <- function(df) {
          # Get column names
          col_names <- colnames(df)

          # Initialize new names vector
          new_col_names <- col_names

          # Iterate over columns
          for (i in seq_along(col_names)) {
            col <- col_names[i]
            if (startsWith(col, "Hdist_gap")) {
              # Find the preceding dist column
              prev_dist_col_index <- which(col_names == col) - 1

              # Check if the preceding column is a dist column
              if (prev_dist_col_index > 0 && startsWith(col_names[prev_dist_col_index], "dist")) {
                # Extract numeric suffix from the dist column
                suffix <- gsub("\\D", "", col_names[prev_dist_col_index])

                # Rename Hdist_gap column to Hdist with the same suffix as the preceding dist column
                new_col_name <- paste0("Hdist", suffix)
                new_col_names[i] <- new_col_name
              } else {
               # warning(paste0("No corresponding dist column found for ", col, ". Skipping renaming."))
              }
            }
          }

          # Rename columns in the dataframe
          colnames(df) <- new_col_names

          return(df)
        }

        # Apply the function
        distance_data_kk <- rename_hdist_gap_columns(distance_data1)
      }

      if (exists("distance_data_kk")) {

      distance_data_kk1<-distance_data_kk

      Hdist_gap_cols1 <- grep("^Hdist_gap", colnames(distance_data_kk1), value =T)

      # Check if there is exactly one Hdist_gap column, no "dist" column exists, and gap1 > cbh1

      if ((length(Hdist_gap_cols1) == 1) && (!any(grepl("^dist1", colnames(distance_data_kk1)))) && (min(kk_copy$gap1) > min(kk_copy$cbh1))) {

        # Extract the value from Hdist column in distance_data_kk1
        hdist_value <- distance_data_kk1[, Hdist_gap_cols1]

        # Find the cbh columns in kk_copy
        cbh_cols <- grep("^cbh", colnames(kk_copy), value = TRUE)

        # Find the closest cbh value that is greater than the hdist_value
        closest_cbh_col <- cbh_cols[which.min(ifelse(kk_copy[, cbh_cols] > hdist_value, kk_copy[, cbh_cols], Inf))]

        # Extract the suffix from the closest cbh column
        suffix <- gsub("\\D", "", closest_cbh_col)

        # Rename the Hdist_gap column in distance_data_kk1 with the new name
        new_hdist_col <- paste0("Hdist", suffix)
        colnames(distance_data_kk1)[colnames(distance_data_kk1) == Hdist_gap_cols1] <- new_hdist_col


        # Find the cbh columns in kk_copy
        cbh_cols <- grep("^cbh", colnames(kk_copy), value = TRUE)

        # Find the closest cbh value that is less than the hdist_value
        closest_cbh_col <- cbh_cols[which.max(ifelse(kk_copy[, cbh_cols] < hdist_value, kk_copy[, cbh_cols], -Inf))]
        closest_cbh_col_value <-kk_copy[, closest_cbh_col]

        # Check if the corresponding dist column does not exist
        dist_col_name <- paste0("dist", suffix)
        if (!dist_col_name %in% colnames(distance_data_kk1)) {
          # Create the corresponding dist column with the same suffix and value equal to floor of Hdist column
          distance_data_kk1[[dist_col_name]] <- floor(hdist_value - closest_cbh_col_value)
        }

        # #print the updated distance_data_kk
        #print(distance_data_kk1)
      }
      }

      if (exists("distance_data_kk1")){
      # Function to order columns
      order_columns <- function(df) {
        col_names <- colnames(df)

        # Extract Hdist and dist columns
        hdist_cols <- grep("^Hdist", col_names, value = TRUE)
        dist_cols <- grep("^dist", col_names, value = TRUE)

        # Extract numeric suffixes
        hdist_suffixes <- as.numeric(gsub("\\D", "", hdist_cols))
        dist_suffixes <- as.numeric(gsub("\\D", "", dist_cols))

        # Create a named vector for sorting
        all_suffixes <- sort(unique(c(hdist_suffixes, dist_suffixes)))

        # Initialize an empty vector for ordered column names
        ordered_cols <- character(0)

        # Loop through each suffix and add corresponding Hdist and dist columns
        for (suffix in all_suffixes) {
          hdist_col <- paste0("Hdist", suffix)
          dist_col <- paste0("dist", suffix)

          if (hdist_col %in% col_names) {
            ordered_cols <- c(ordered_cols, hdist_col)
          }
          if (dist_col %in% col_names) {
            ordered_cols <- c(ordered_cols, dist_col)
          }
        }

        # Reorder the dataframe columns
        df <- df[, ordered_cols, drop = FALSE]

        return(df)
      }

      # Apply the function
      distance_data_kk1 <- order_columns(distance_data_kk1)
      }

      if (!exists("distance_data_kk") && !exists("distance_data_kk1")) {
      # If Hdist_gap columns are found, proceed with renaming
      if (length(Hdist_gap_cols) == 1 && ("dist" %in% colnames(distance_data1))) {

        rename_hdist_gap_columns1 <- function(df) {
          # Get column names
          col_names <- colnames(df)

          # Initialize new names vector
          new_col_names <- col_names

          # Iterate over columns
          for (i in seq_along(col_names)) {
            col <- col_names[i]
            if (startsWith(col, "Hdist_gap")) {
              # Find the corresponding dist column
              dist_col_index <- which(startsWith(col_names, "dist"))

              if (length(dist_col_index) > 0) {
                # Extract numeric suffix from the dist column
                dist_col <- col_names[dist_col_index[1]]
                suffix <- gsub("\\D", "", dist_col)

                # Rename Hdist_gap column to Hdist with the same suffix as the dist column
                new_col_name <- paste0("Hdist", suffix)
                new_col_names[i] <- new_col_name
              } else {
               # warning(paste0("No corresponding dist column found for ", col, ". Skipping renaming."))
              }
            }
          }

          # Rename columns in the dataframe
          colnames(df) <- new_col_names

          return(df)
        }

        # Apply the function
        distance_data_kk2 <- rename_hdist_gap_columns1(distance_data1)
      }
      }

      if (!exists("distance_data_kk") && !exists("distance_data_kk1") && !exists("distance_data_kk2")) {
      if (length(Hdist_gap_cols) == 1 && (!"dist" %in% colnames(distance_data1))) {

         # Function to rename Hdist_gap column and create corresponding dist column
      rename_and_create_dist_column <- function(df) {
        col_names <- colnames(df)

        # Find the Hdist_gap column
        hdist_gap_col <- grep("^Hdist_gap", col_names, value = TRUE)

        if (length(hdist_gap_col) == 1) {
          # Extract numeric suffix from the Hdist_gap column
          suffix <- gsub("\\D", "", hdist_gap_col)

          # Rename Hdist_gap column to Hdist with the same suffix
          new_col_name <- paste0("Hdist", suffix)
          colnames(df)[colnames(df) == hdist_gap_col] <- new_col_name

          # Create the corresponding dist column with the same suffix and value equal to floor of Hdist column
          dist_col_name <- paste0("dist", suffix)
          df[[dist_col_name]] <- floor(df[[new_col_name]])
        } else {
         # warning("There is more than one Hdist_gap column or none found.")
        }

        return(df)
      }
      # Apply the function
      distance_data_kk3 <- rename_and_create_dist_column(distance_data1)

      }
      }


      if (!exists("distance_data_kk") && !exists("distance_data_kk1") && !exists("distance_data_kk2") && !exists("distance_data_kk3")) {

         Hdist_cols <- grep("^Hdist", colnames(distance_data1))

        if (length(Hdist_cols) == 1 && (!"dist" %in% colnames(distance_data1))) {

          # Function to rename Hdist_gap column and create corresponding dist column
          rename_and_create_dist_column <- function(df) {
            col_names <- colnames(df)

            # Find the Hdist_gap column
            hdist_col <- grep("^Hdist", col_names, value = TRUE)

            if (length(hdist_col) == 1) {
              # Extract numeric suffix from the Hdist_gap column
              suffix <- gsub("\\D", "", hdist_col)

              # Rename Hdist_gap column to Hdist with the same suffix
              new_col_name <- paste0("Hdist", suffix)
              colnames(df)[colnames(df) == hdist_col] <- new_col_name

              # Create the corresponding dist column with the same suffix and value equal to floor of Hdist column
              dist_col_name <- paste0("dist", suffix)
              df[[dist_col_name]] <- floor(df[[new_col_name]])
            } else {
              # warning("There is more than one Hdist_gap column or none found.")
            }

            return(df)
          }
          # Apply the function
          distance_data_kk4 <- rename_and_create_dist_column(distance_data1)

        }
      }


      if (exists("distance_data_kk") && !exists("distance_data_kk1")) {
        distance_data <-distance_data_kk
      }

      if (exists("distance_data_kk") && exists("distance_data_kk1")) {
        distance_data <-distance_data_kk1
      }

      if (exists("distance_data_kk2")) {
        distance_data <-distance_data_kk2
      }

      if (exists("distance_data_kk3")) {
        distance_data <-distance_data_kk3
      }

      if (exists("distance_data_kk4")) {
        distance_data <-distance_data_kk4
      }

      if (!exists("distance_data_kk") && !exists("distance_data_kk1") && !exists("distance_data_kk2") && !exists("distance_data_kk3")  && !exists("distance_data_kk4")) {
        distance_data <-distance_data1
      }

   ########################################


      if (min(kk_copy$cbh1) <= min_height) {

        distance_data$Hdist1 <-kk_copy$cbh1
        distance_data$dist1 <- 0
      }

      # Find the first cbh value in kk_copy
      first_cbh <- kk_copy$cbh1

      # Check if there is any Hdist column with a value less than the first cbh
      hdist_cols <- grep("^Hdist", colnames(distance_data), value = TRUE)
      any_hdist_less_than_first_cbh <- any(sapply(distance_data[, hdist_cols, drop = FALSE], function(x) x < first_cbh))

      # If no Hdist column has a value less than the first cbh, create new columns Hdist1 and dist1
      if (!any_hdist_less_than_first_cbh && first_cbh > min_height) {
        distance_data$Hdist1 <- first_cbh -step
        distance_data$dist1 <- floor(distance_data$Hdist1)
      }

      if((distance_data$Hdist1 <= min_height) && (distance_data$dist1 > distance_data$Hdist1)) {
        distance_data$Hdist1<- (distance_data$Hdist1 +distance_data$dist1) -step
      }

      if((distance_data$dist1 > floor(distance_data$Hdist1))) {
        distance_data$dist1<- distance_data$dist1 -step
      }


    # Filter out "Hdist" columns with values greater than max_cbh_value

    cbh_cols <- grepl("^cbh", names(kk_copy))

    # Extract the maximum value from the "cbh" columns
    max_cbh_value <- max(kk_copy[, cbh_cols])

    # Identify columns starting with "Hdist_gap"
    hdist_cols <- grepl("^Hdist", names(distance_data))
    dist_cols <- grepl("^dist", names(distance_data))

    # Filter out "Hdist_gap" columns with values greater than max_cbh_value
    cols_to_keep <- !(hdist_cols & distance_data[1, ] > max_cbh_value)

    # Subset the dataframe to keep only the required columns
    distance_data <- distance_data[, cols_to_keep]

   ##################################
    # Identify columns starting with "dist" and "Hdist"
    dist_cols <- grep("^dist", names(distance_data), value = TRUE)
    Hdist_cols <- grep("^Hdist", names(distance_data), value = TRUE)

    # Extract numeric suffixes from dist and Hdist columns
    dist_suffixes <- gsub("\\D", "", dist_cols)
    Hdist_suffixes <- gsub("\\D", "", Hdist_cols)

    # Find dist columns where suffixes do not match with Hdist columns
    dist_to_remove <- dist_cols[!(dist_suffixes %in% Hdist_suffixes)]

    # Remove mismatched dist columns
    distance_data_df <- distance_data[, !(names(distance_data) %in% dist_to_remove)]

    ##################################
     ########################################################

    # Extract column names
    dist_columns <- names(distance_data_df)[grepl("^dist", names(distance_data_df))]
    Hdist_columns <- names(distance_data_df)[grepl("^Hdist", names(distance_data_df))]

    # Extract numeric suffixes
    dist_suffixes <- sub("^dist", "", dist_columns)
    Hdist_suffixes <- sub("^Hdist", "", Hdist_columns)

    # Create a mapping for Hdist suffixes to match dist suffixes
    suffix_mapping <- match(dist_suffixes, Hdist_suffixes)
    new_Hdist_columns <- paste0("Hdist", dist_suffixes)

    # Rename Hdist columns
    names(distance_data_df)[grepl("^Hdist", names(distance_data_df))] <- new_Hdist_columns

    ########################################

distance_metrics <- cbind.data.frame(treeID1,treeID, kk_copy, distance_data_df,max_height)

# List of columns to be rounded
dist_columns <- grep("^dist", names(distance_metrics), value = TRUE)

# Apply round function to the selected columns
distance_metrics[dist_columns] <- lapply(distance_metrics[dist_columns], floor)

########################################
distance_metrics1<-distance_metrics

# Extract column names
cbh_cols <- grep("^cbh", names(distance_metrics1), value = TRUE)

# Function to determine if any Hdist <= cbh value
check_hdist_le_cbh <- function(hdist_values, cbh_value) {
  any(!is.na(hdist_values) & hdist_values <= cbh_value)
}

# Iterate over cbh columns
for (cbh_col in cbh_cols) {
  hdist_col_suffix <- sub("^cbh", "", cbh_col)
  hdist_col <- paste0("Hdist", hdist_col_suffix)
  dist_col_suffix <- sub("^cbh", "", cbh_col)
  dist_col <- paste0("dist", hdist_col_suffix)

  # Check if hdist_col exists
  if (hdist_col %in% names(distance_metrics1)) {
    cbh_values <- distance_metrics1[[cbh_col]]
    hdist_values <- distance_metrics1[[hdist_col]]
    dist_values <- distance_metrics1[[dist_col]]

    # Check if any Hdist value is <= cbh value
    if (any(hdist_values <= cbh_values)) {
      # Update Hdist column
      distance_metrics1[[hdist_col]] <- hdist_values
    }
  } else {
    # Create new Hdist column with 0 if hdist_col doesn't exist
    distance_metrics1[[hdist_col]] <- 0
  }

  # Check if dist_col exists
  if (dist_col %in% names(distance_metrics1)) {
    dist_values <- distance_metrics1[[dist_col]]

    # Check if any dist value is >= 0
    if (any(dist_values >= 0)) {
      # Update dist column
      distance_metrics1[[dist_col]] <- dist_values
    }
  } else {
    # Create new dist column with 0 if dist_col doesn't exist
    distance_metrics1[[dist_col]] <- 0
  }
}


##########################################################
#df<- distance_metrics1

## remove Hdist values > max(cbh) and the last gap and cbh
# Identify columns starting with "cbh" and calculate their maximum
cbh_cols <- grep("^cbh", names(distance_metrics1), value = TRUE)
max_cbh <- max(unlist(distance_metrics1[cbh_cols]))

# Identify columns starting with "Hdist"
remove_cols <- grep("^Hdist", names(distance_metrics1), value = TRUE)

# Keep the columns that should not be removed
remove_cols_to_keep <- character()

for (col in remove_cols) {
  # Check if any value in the column is less than or equal to max_cbh
  if (any(distance_metrics1[[col]] <= max_cbh)) {
    remove_cols_to_keep <- c(remove_cols_to_keep, col)
  }
}

# Columns to remove
cols_to_remove <- setdiff(remove_cols, remove_cols_to_keep)

if (length(cols_to_remove) > 0) {
# Remove the identified columns
distance_metrics2 <- distance_metrics1[, !(names(distance_metrics1) %in% cols_to_remove)]
} else {
  distance_metrics2 <- distance_metrics1
}

# Extract the number and names of cbh columns to match the length of Hdist cols

cbh_columns <- grep("^cbh", names(distance_metrics2), value = TRUE)
hdist_columns <- grep("^Hdist", names(distance_metrics2), value = TRUE)

if (length(cbh_columns) > length(hdist_columns)) {

# Determine the number of cbh columns
num_cbh_cols <- length(cbh_columns)

# Generate a new Hdist column name
new_Hdist_col <- paste0("Hdist", num_cbh_cols + 1)

# Add Hdist column with initial value 0
distance_metrics2[[new_Hdist_col]] <- rep(0, nrow(distance_metrics2))

}

###################################################
#Let's refine the approach to ensure that each Hdist column with a value greater than 0 is correctly paired with its corresponding cbh column,
#and then rename the remaining Hdist columns (with a value of 0) with the suffixes of the unpaired cbh columns.

result <- pair_cbh_Hdist(distance_metrics2)

# Extract the modified dataframe with values
distance_metrics_modified <- result[[1]]

# Function to retain Hdist columns where at least one value is non-zero

retain_non_zero_hdist_columns <- function(df) {
  hdist_cols <- grepl("^Hdist", names(df))  # Identify columns starting with 'Hdist'

  # Subset dataframe to retain only Hdist columns with at least one non-zero value
  df <- df[, c(names(df)[!hdist_cols], names(df)[hdist_cols][colSums(df[hdist_cols]) != 0])]

  return(df)
}

# Apply the function to your dataframe
distance_metrics_modified1 <- retain_non_zero_hdist_columns(distance_metrics_modified)

#############################

# Check how many elements in the result list are data frames
num_dfs <- sum(sapply(result, is.data.frame))

# Perform actions based on the number of data frames
if (num_dfs > 1) {

# Extract the cleaned dataframe with remaining columns
distance_metrics_remaining <- result[[2]]

distance_metrics_remaining1 <- distance_metrics_remaining %>% dplyr::select(matches("^Hdist"))

distance_metrics_def<-data.frame(distance_metrics_modified1,distance_metrics_remaining1)


} else {

  distance_metrics_def<-data.frame(distance_metrics_modified1)
  }

  } else {

      if (min(kk_copy$cbh1) > min_height) {

        kk_copy$Hdist1 <- min(kk_copy$cbh1) -step
        kk_copy$dist1 <- floor(kk_copy$Hdist1)
      }

      if (min(kk_copy$cbh1) <= min_height) {

        kk_copy$Hdist1 <-kk_copy$cbh1
        kk_copy$dist1 <- 0
      }


    # Extract column names
    cbh_cols <- grep("^cbh", names(kk_copy), value = TRUE)

    # Function to determine if any Hdist <= cbh value
    check_hdist_le_cbh <- function(hdist_values, cbh_value) {
      any(!is.na(hdist_values) & hdist_values <= cbh_value)
    }

    # Iterate over cbh columns
    for (cbh_col in cbh_cols) {
      hdist_col_suffix <- sub("^cbh", "", cbh_col)
      hdist_col <- paste0("Hdist", hdist_col_suffix)
      dist_col_suffix <- sub("^cbh", "", cbh_col)
      dist_col <- paste0("dist", hdist_col_suffix)

      # Check if hdist_col exists
      if (hdist_col %in% names(kk_copy)) {
        cbh_values <- kk_copy[[cbh_col]]
        hdist_values <- kk_copy[[hdist_col]]
        dist_values <- kk_copy[[dist_col]]

        # Check if any Hdist value is <= cbh value
        if (any(hdist_values <= cbh_values)) {
          # Update Hdist column
          kk_copy[[hdist_col]] <- hdist_values
        }
      } else {
        # Create new Hdist column with 0 if hdist_col doesn't exist
        kk_copy[[hdist_col]] <- 0
      }

      # Check if dist_col exists
      if (dist_col %in% names(kk_copy)) {
        dist_values <- kk_copy[[dist_col]]

        # Check if any dist value is >= 0
        if (any(dist_values >= 0)) {
          # Update dist column
          kk_copy[[dist_col]] <- dist_values
        }
      } else {
        # Create new dist column with 0 if dist_col doesn't exist
        kk_copy[[dist_col]] <- 0
      }
    }


    treeID<-unique(factor(df$treeID))
    treeID1<-unique(factor(df$treeID1))

    max_height<-data.frame(df$max_height)
    names(max_height)<-"max_height"

    distance_metrics_def <- cbind.data.frame(treeID, treeID1, kk_copy, max_height)

  }

 ##########################################################

    # Identify the columns for cbh and gap values
    cbh_cols <- grep("^cbh", colnames(distance_metrics_def), value = TRUE)
    gap_cols <- grep("^gap", colnames(distance_metrics_def), value = TRUE)

    # Initialize a data frame to store the results
    results <- data.frame(gap = numeric(), cbh = numeric())

    # Loop through each gap value
    for (gap_col in gap_cols) {
      gap_value <- distance_metrics_def[, gap_col]

      # Find all cbh values greater than the current gap value
      valid_cbhs <- distance_metrics_def[, cbh_cols][distance_metrics_def[, cbh_cols] > gap_value]

      if (length(valid_cbhs) > 0) {
        # Find the smallest cbh value that is greater than the gap value
        closest_cbh <- min(valid_cbhs)

        # Add the result to the data frame
        results <- rbind(results, data.frame(gap = gap_value, cbh = closest_cbh))
      } else {
        # In case no cbh value is greater than the gap value
        results <- rbind(results, data.frame(gap = gap_value, cbh = NA))
      }
    }

    # Compute the difference and store it in a new column called dist
    # Calculate the new 'dist' variable and keep only one row per duplicated cbh
    results1 <- results %>%
      dplyr::group_by(cbh) %>%
      dplyr::filter(gap == min(gap)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(dist = cbh - gap)

    results2 <- results %>%
      dplyr::group_by(cbh) %>%
      dplyr::filter(gap == max(gap)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(Hdist = gap)

    # Extract the desired columns
    cbh <- results1$cbh
    dist <- results1$dist
    Hdist <- results2$Hdist

    # Combine the columns into a new data frame
    resultsfi <- cbind.data.frame(cbh, dist, Hdist)



    # Loop through each result
    for (i in 1:nrow(resultsfi)) {
      Hdist_value <- resultsfi$Hdist[i]
      cbh_value <- resultsfi$cbh[i]
      dist_value <- resultsfi$dist[i]

      if (!is.na(cbh_value)) {
        # Find the corresponding Hdist column name based on cbh_value
        hdist_col <- paste0("Hdist", which(distance_metrics_def[1, grep("^cbh", names(distance_metrics_def))] == cbh_value))
        dist_col <- paste0("dist", sub("Hdist", "", hdist_col))

        # Check if the gap_value is already in the Hdist column
        if (!(Hdist_value %in% distance_metrics_def[, hdist_col])) {
          # Add the gap value to the Hdist column
          distance_metrics_def[, hdist_col] <- Hdist_value
          distance_metrics_def[, dist_col] <- dist_value
        }
      }
    }


#####################################################################

    # Find the cbh and Hdist columns: Update Hdist values
    cbh_cols <- grep("^cbh", colnames(distance_metrics_def), value = TRUE)
    hdist_cols <- grep("^Hdist", colnames(distance_metrics_def), value = TRUE)

    # Iterate over cbh columns
    for (cbh_col in cbh_cols) {
      # Extract the numeric suffix
      suffix <- sub("^cbh", "", cbh_col)

      # Find the corresponding Hdist column
      hdist_col <- paste0("Hdist", suffix)

      # Check if the Hdist column exists
      if (hdist_col %in% colnames(distance_metrics_def) ) {
        if(any(distance_metrics_def[[hdist_col]] > min_height)) {
        # Update Hdist values
        distance_metrics_def[[hdist_col]] <- ifelse(distance_metrics_def[[hdist_col]] != distance_metrics_def[[cbh_col]] - step,
                                                    distance_metrics_def[[cbh_col]] - step,
                                                    distance_metrics_def[[hdist_col]])
        }
      }
    }

    if(distance_metrics_def$cbh1 > min_height) {
      distance_metrics_def$Hdist1 <-distance_metrics_def$cbh1 -step
      distance_metrics_def$dist1 <-floor(distance_metrics_def$Hdist1)
    }

    if(distance_metrics_def$cbh1 <= min_height) {
      distance_metrics_def$Hdist1 <-min_height
      distance_metrics_def$dist1 <-0
    }


return(distance_metrics_def)
}

