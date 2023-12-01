#'
#' Distances between fuel layers
#' @description This function calculates distances (and their heights) between fuel layers as the difference between consecutive gaps and fuel bases
#' (the gap height always must be lower than the fuel base height).
#' @usage get_distance (gap_cbh_metrics)
#'
#' @param gap_cbh_metrics data frame with gaps (distances) and fuel base heights (output of [get_gaps_fbhs()] function).
#' An object of the class text
#' @return A data frame giving distances (and their heights) between fuel layers in meters.
#' @author Olga Viedma, Carlos Silva and JM Moreno
#'
#' @details
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
#' @examples
#' ## Not run:
#' library(magrittr)
#' library(gdata)
#' library(dplyr)
#'
#' # Load the effective_distances object
#' if (interactive()) {
#'   gap_cbh_metrics <- get_gaps_fbhs()
#'   LadderFuelsR::gap_cbh_metrics$treeID <- factor(LadderFuelsR::gap_cbh_metrics$treeID)
#'
#' metrics_distance_list <- list()
#'
#' for (i in levels(trees_name2)) {
#'
#'   # Filter data for each tree
#'   tree2 <- gap_cbh_metrics |> dplyr::filter(treeID == i)
#'
#'   # Get distance metrics for each tree
#'   metrics_distance <- get_distance(tree2)
#'   metrics_distance_list[[i]] <- metrics_distance
#' }
#'
#' # Combine the individual data frames
#' distance_metrics <- dplyr::bind_rows(metrics_distance_list)
#'}
#' ## End(Not run)
#' @export get_distance
#' @importFrom dplyr select_if group_by summarise summarize mutate arrange rename rename_with filter slice ungroup distinct
#' @importFrom magrittr %>%
#' @importFrom gdata startsWith
#' @include gap_fbh.R
get_distance <- function (gap_cbh_metrics) {

  df <- gap_cbh_metrics %>%
    dplyr::mutate_at(
      vars(-treeID),  # Exclude the 'treeID' column
      as.numeric
    )

  df1 <- df[, !colSums(is.na(df)) > 0]
  # Select only numeric columns
  df1_numeric <- df1 %>% dplyr::select_if(is.numeric)

  #print(paste("treeID:", df1_numeric$treeID))

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

####################################

  if(length(gap_cols) >= 1 & length(cbh_cols) == 0 && (!exists("distance_data") || is.null(distance_data) || ncol(distance_data) == 0 && nrow(distance_data) == 0 || is.na(distance_data)))  {

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
  gap_cols <- grep("^Hgap", colnames(kk_copy))
  cbh_cols <- grep("^Hcbh", colnames(kk_copy))

  # Check the number of gap columns is greater than 1 and number of cbh columns equals to 1
  if(length(gap_cols) > 1 & length(cbh_cols) == 1 ) {

    # Check the number of gap columns is greater than 1
    if (length(gap_cols) > 1) {
      # Get the column names
      col_names <- colnames(kk_copy)

      # Identify columns starting with "gap" and "cbh"
      gap_cols <- grep("^Hgap", col_names)
      cbh_cols <- grep("^Hcbh", col_names)

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

    # Add Hdist4 column to the data frame if it exists and is not empty
    if (exists("Hdist4") && !is.null(Hdist4) && ncol(Hdist4) != 0 && nrow(Hdist4) != 0 && !is.na(Hdist4)) {
      distance_data3 <- cbind(distance_data3,Hdist4)
      colnames(distance_data3)<-c("dist", "Hdist")
    }
  }

  ###############################


  if (!exists("distance_data")) {
    if (exists("distance_data1")) {
      if(!is.null(distance_data1) || ncol(distance_data1) != 0 || nrow(distance_data1) != 0 || !is.na(distance_data1)) {
        if(!exists("distance_data2") || is.null(distance_data2) || any(is.na(distance_data2))) {
          if (!exists("distance_data3") || is.null(distance_data3) || any(is.na(distance_data3))) {

            distance_data <-data.frame(distance_data1)
          }}}}}



  if (!exists("distance_data")) {
    if (exists("distance_data2")) {
      if(!is.null(distance_data2) || ncol(distance_data2) != 0 || nrow(distance_data2) != 0 || !is.na(distance_data2)) {
        if(!exists("distance_data1") || is.null(distance_data1) || any(is.na(distance_data1))) {
          if (!exists("distance_data3") || is.null(distance_data3) || any(is.na(distance_data3))) {


            distance_data <-data.frame(distance_data2)
          }}}}}


  if (!exists("distance_data")) {
    if (exists("distance_data3")) {
      if(!is.null(distance_data3) || ncol(distance_data3) != 0 || nrow(distance_data3) != 0 || !is.na(distance_data3)) {
        if(!exists("distance_data1") || is.null(distance_data1) || any(is.na(distance_data1))) {
          if (!exists("distance_data2") || is.null(distance_data2) || any(is.na(distance_data2))) {

            distance_data <-data.frame(distance_data3)
          }}}}}


  if (!exists("distance_data")) {
    if (exists("distance_data1")) {
      if(!is.null(distance_data1) || ncol(distance_data1) != 0 || nrow(distance_data1) != 0 || !is.na(distance_data1)) {
        if(exists("distance_data2")) {
          if(!is.null(distance_data2) || any(!is.na(distance_data2))) {
            if (!exists("distance_data3") || is.null(distance_data3) || any(is.na(distance_data3))) {

              distance_data <-data.frame(distance_data1, distance_data2)
            }}}}}}


  if (!exists("distance_data")) {
    if (exists("distance_data1")) {
      if(!is.null(distance_data1) || ncol(distance_data1) != 0 || nrow(distance_data1) != 0 || !is.na(distance_data1)) {
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

  if ((is.null(distance_data) || ncol(distance_data) == 0 || nrow(distance_data) == 0 || is.na(distance_data)) &&
      (exists("distance_data1") && !is.null(distance_data1) && ncol(distance_data1) != 0 && nrow(distance_data1) != 0 && !is.na(distance_data1)) &&
      (exists("distance_data2") && !is.null(distance_data2) && ncol(distance_data2) != 0 && nrow(distance_data2) != 0 && !is.na(distance_data2))) {


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

    if ((is.null(distance_data) || ncol(distance_data) == 0 || nrow(distance_data) == 0 || is.na(distance_data)) &&
        (exists("distance_data1") && !is.null(distance_data1) && ncol(distance_data1) != 0 && nrow(distance_data1) != 0 && !is.na(distance_data1)) &&
        (exists("distance_data2") && !is.null(distance_data2) && ncol(distance_data2) != 0 && nrow(distance_data2) != 0 && !is.na(distance_data2))) {

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

  if ((is.null(distance_data) || ncol(distance_data) == 0 || nrow(distance_data) == 0 || is.na(distance_data)) &&
      (exists("distance_data1") && !is.null(distance_data1) && ncol(distance_data1) != 0 && nrow(distance_data1) != 0 && !is.na(distance_data1)) &&
      (exists("distance_data3") && !is.null(distance_data3) && ncol(distance_data3) != 0 && nrow(distance_data3) != 0 && !is.na(distance_data3))) {


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

    if ((is.null(distance_data) || ncol(distance_data) == 0 || nrow(distance_data) == 0 || is.na(distance_data)) &&
        (exists("distance_data1") && !is.null(distance_data1) && ncol(distance_data1) != 0 && nrow(distance_data1) != 0 && !is.na(distance_data1)) &&
        (exists("distance_data3") && !is.null(distance_data3) && ncol(distance_data3) != 0 && nrow(distance_data3) != 0 && !is.na(distance_data3))) {

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
      if (!is.na(distance_data1) && length(distance_data1) > 0 && any(!is.na(distance_data1))) {

        if (exists("distance_data2")) {
          if (!is.na(distance_data2) && length(distance_data2) > 0 && any(!is.na(distance_data2))) {

            if (exists("distance_data3")) {
              if (!is.na(distance_data3) && length(distance_data3) > 0 && any(!is.na(distance_data3))) {

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
        if (!is.na(distance_data2) && length(distance_data2) > 0 && any(!is.na(distance_data2))) {

          if (exists("distance_data3")) {
            if (!is.na(distance_data3) && length(distance_data3) > 0 && any(!is.na(distance_data3))) {

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

        if (exists("distance_data3") && !is.na(distance_data3) && length(distance_data3) > 0 && any(!is.na(distance_data3))) {

          if(ncol(distance_data)==1) {

            distance_data <- distance_data3

          } else {
            distance_data <- distance_data
          }
        }
      }}}


  if ((exists("distance_data") && nrow(distance_data) != 0 && any(!is.na(distance_data)))) {

    if  (exists("distance_data1") && !is.na(distance_data1) && length(distance_data1) != 0 && any(!is.na(distance_data1))) {
      if (exists("distance_data2") && !is.na(distance_data2) && length(distance_data2) != 0 && any(!is.na(distance_data2))) {
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

    if (exists("distance_data1") && is.na(distance_data1) && length(distance_data1) == 0 && any(is.na(distance_data1))) {
      if (!exists("distance_data2")) {
        if (exists("distance_data3") && !is.na(distance_data3) && length(distance_data3) != 0 && any(!is.na(distance_data3))) {

          if(ncol(distance_data)==1) {

            if(distance_data1[,2] != distance_data3[,2]) {

              distance_data <- distance_data1

            } else {
              distance_data <- distance_data
            }
          }
        }}}}


  if (length(gap_cols) > 1 && length(cbh_cols) > 1 && ((exists("distance_data") || !is.null(distance_data) || (ncol(distance_data) != 0 && nrow(distance_data) != 0)) || !is.na(distance_data))) {

    percent2a <- gaps_perc2 %>% dplyr::filter(height < min(kk_copy[, cbh_cols]))

    if (nrow(percent2a) != 0 && any(!is.na(percent2a)) && any(percent2a$percentil > 25)) {
      distance_data <- as.data.frame(distance_data[,-1])
    }
  }


  #### IF THERE ARE CONSECUTIVE GAPS AFTER CBH COLUMNS (ABOVE DISTANCES) #################

  gap_cols <- grep("^Hgap", colnames(kk_copy))
  cbh_cols <- grep("^Hcbh", colnames(kk_copy))

  # remove missing values in gap_cols or cbh_cols
  gap_cols <- gap_cols[!is.na(kk_copy[,gap_cols])]
  cbh_cols <- cbh_cols[!is.na(kk_copy[,cbh_cols])]

  ###################################

  if (length(gap_cols) > 1 && length(cbh_cols) > 1 &&
      ((exists("distance_data") || !is.null(distance_data) ||
        (ncol(distance_data) != 0 && nrow(distance_data) != 0)) ||
       !is.na(distance_data))) {

    if (length(gap_cols) > 1 && any(diff(gap_cols) > 1)) {

      # Get the positions of 'gap' and 'cbh' columns
      gap_positions <- grep("^Hgap", names(kk_copy))
      cbh_positions <- grep("^Hcbh", names(kk_copy))

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
          gap_col1<- grep("^Hgap", colnames(Hdist5))

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



  if ((exists("distance_data") && !is.null(distance_data) && ncol(distance_data) != 0 && nrow(distance_data) != 0 && !is.na(distance_data)) &&
      (exists("distance_data4") && !is.null(distance_data4) && ncol(distance_data4) != 0 && nrow(distance_data4) != 0 && !is.na(distance_data4))) {

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


  if(exists("distance_data") && !is.null(distance_data) && ncol(distance_data) != 0 && nrow(distance_data) != 0 && !is.na(distance_data)
     && all(distance_data==0)) {

    distance_data <-data.frame(c(NA))
    names(distance_data)= "dist_0"
  }

  if (length(gap_cols) > 1 && length(cbh_cols) > 1 && ((!exists("distance_data") || is.null(distance_data) || (ncol(distance_data) == 0 && nrow(distance_data) == 0)) || is.na(distance_data))) {

    diff2 <- max(kk_copy[, gap_cols]) - min(kk_copy[, gap_cols])
    Hdist6<-max(kk_copy[, gap_cols])
    distance_data <- cbind.data.frame(diff2,Hdist6)
  }



  if(all(!is.na(distance_data))){

    Hdist_cols <- grep("^(Hdist)", names(distance_data), value = TRUE)
    Hdist_vals <- distance_data[, Hdist_cols, drop = FALSE]

    # Convert the dataframe into a vector
    vec <- unlist(Hdist_vals)

    # Identify duplicated values
    dups <- unique(vec[duplicated(vec)])

    # For each duplicated value, remove the last occurrence
    for (dup in dups) {
      # Get the indices of the duplicated values
      indices <- which(vec == dup)
      # Remove the last occurrence
      vec[tail(indices, n=1)] <- NA
    }

    # Convert the vector back to a dataframe
    Hdist_vals_clean <- as.data.frame(t(vec))
    Hdist_vals_clean <- Hdist_vals_clean[, !apply(is.na(Hdist_vals_clean), 2, all)]

    distance_data_p <- distance_data[,-which(names(distance_data) %in% Hdist_cols)]
    distance_data<-cbind(distance_data_p,Hdist_vals_clean)

  }


  if(all(is.na(distance_data)) && all(!is.na(distance_data1)) && (!exists("distance_data3") && !exists("distance_data4"))){
    Hdist_cols <- grep("^(Hdist)", names(distance_data1), value = TRUE)
    distance_data <- distance_data1[Hdist_cols]
  }

  if (length(distance_data) == 0) {
    distance_data <- data.frame(distance_data = 0)
    names(distance_data)="dist0"
  }


  # Check if distance_data is not empty
  if (length(distance_data) != 0) {

    distance_data <- data.frame(distance_data)

    old_names_distance_data <- colnames(distance_data)
    new_names_distance_data <- vector(mode="character", length=length(old_names_distance_data))  # initialize new_names_distance_data

    # Initialize counters for 'dist' and 'Hdist' column indices
    dist_index <- 1
    Hdist_index <- 1

    for (i in seq_along(old_names_distance_data)) {
      if (startsWith(old_names_distance_data[i], "gap")) {
        new_names_distance_data[i] <- paste("dist", dist_index, sep="")
        dist_index <- dist_index + 1
      } else if (startsWith(old_names_distance_data[i], "Hdist")) {
        new_names_distance_data[i] <- paste("Hdist", Hdist_index, sep="")
        Hdist_index <- Hdist_index + 1
      } else if (startsWith(old_names_distance_data[i], "dist")) {
        new_names_distance_data[i] <- paste("dist", dist_index, sep="")
        dist_index <- dist_index + 1
      }

      colnames(distance_data) <- new_names_distance_data
    }
  }


  treeID<-unique(factor(df$treeID))
  treeID1<-unique(factor(df$treeID1))

  max_height<-data.frame(df$max_height)
  names(max_height)<-"max_height"

distance_metrics <- cbind.data.frame(treeID1,treeID, kk_copy, distance_data,max_height)


return(distance_metrics)
}

