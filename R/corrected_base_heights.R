# Declare global variables to avoid R CMD check warnings
utils::globalVariables(c("treeID2", "max1"))
#' Fuels base height recalculation after after considering distances greater than any number of height bin steps
#' @description
#' This function reshapes fuel layers after removing distances equal to any number of height bin steps,
#' keeping the first "base height" from those consecutive ones separated by such distance.
#' @usage
#' get_real_fbh(depth_metrics, step= 1, number_steps = 1, min_height=1.5, verbose=TRUE)
#' @param depth_metrics
#' Tree metrics with distances, fuel base heights, and depths
#' (output of [get_depths()] function). An object of the class text.
#' @param step Numeric value for the actual height bin step (in meters).
#' @param number_steps Numeric value for the number of height bin steps that can be jumped to reshape fuels layers.
#' @param min_height Numeric value for the actual minimum base height (in meters).
#' @param verbose Logical, indicating whether to display informational messages (default is TRUE).
#' @return
#' A data frame giving the first "base height" from those consecutive ones separated by the number of height bin steps indicated in the function.
#' The function returns new fuel layers separated by distances greater than the indicated number of steps.
#' @author Olga Viedma, Carlos Silva, JM Moreno and A.T. Hudak
#'
#' @details
#' # List of tree metrics:
#' \itemize{
#'   \item treeID: tree ID with strings and numeric values
#'   \item treeID1: tree ID with only numeric values
#'   \item dist: Distance between consecutive fuel layers (m)
#'   \item Hdist - Height of the distance between consecutive fuel layers (m)
#'   \item Hcbh - Base height of each fuel separated by a distance greater than the certain number of steps
#'   \item depth - Depth of fuel layers (m)
#'   \item Hdepth - Height of the depth of fuel layers (m)
#'   \item max_height - Maximum height of the tree profile
#' }
#'
#' @examples
#' library(magrittr)
#' library(dplyr)
#' #Before running this example, make sure to run get_depths()
#' if (interactive()) {
#' depth_metrics <- get_depths()
#' LadderFuelsR::depth_metrics$treeID <- factor(LadderFuelsR::depth_metrics$treeID)
#'
#' trees_name1 <- as.character(depth_metrics$treeID)
#' trees_name2 <- factor(unique(trees_name1))
#'
#' fbh_corr_list <- list()
#'
#' for (i in levels(trees_name2)){
#' # Filter data for each tree
#' tree3 <- depth_metrics |> dplyr::filter(treeID == i)
#' # Get real fbh for each tree
#' fbh_corr <- get_real_fbh(tree3, step= 1, number_steps = 1, min_height=1.5, verbose=TRUE)
#' # Store fbh values in a list
#' fbh_corr_list[[i]] <- fbh_corr
#' }
#'
#' # Combine fbh values for all trees
#' effective_fbh <- dplyr::bind_rows(fbh_corr_list)
#' effective_fbh$treeID <- factor(effective_fbh$treeID)
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
#' @seealso \code{\link{get_depths}}
#' @export
get_real_fbh <- function (depth_metrics, step= 1, number_steps = 1, min_height=1.5,verbose = TRUE) {

  if(min_height==0){
    min_height <-0.5
  }

    df<- depth_metrics

    remove_cols <- function(df) {
      cols_to_remove <- colSums(is.na(df)) > 0
      df <- df[, !cols_to_remove, drop = FALSE]
      return(df)
    }

    df<-remove_cols(df)

if (verbose) {
  message("Unique treeIDs:", paste(unique(df$treeID), collapse = ", "))
}

###################################################

  for (i in 1:nrow(df)) {

    effective_fbh <- NULL  # Initialize to NULL for each row

    current_row <- df[1, ]

    # Extract numeric suffixes from column names
    numeric_suffixes <- as.numeric(gsub("\\D", "", colnames(current_row)))

    # Identify non-numeric columns and replace their numeric_suffixes with Inf
    non_numeric_cols <- is.na(numeric_suffixes)
    numeric_suffixes[non_numeric_cols] <- Inf

    # Order columns by numeric suffix, placing non-numeric columns at the end
    ordered_cols <- order(numeric_suffixes)

    # Reorder the dataframe columns
    current_row <- current_row[, ordered_cols]

    # Make column names unique for the rest of the columns
    colnames(current_row) <- make.names(colnames(current_row))

    hHdist_cols <- grep("^Hdist\\d+$", colnames(current_row))
    hcbh_cols <- grep("^cbh\\d+$", colnames(current_row))

    hHdist_na <- is.na(current_row[hHdist_cols])
    hcbh_na <- is.na(current_row[hcbh_cols])

    # If all values in current row for HHdist and Hcbh columns are NA, skip this iteration
    if (all(hHdist_na) && all(hcbh_na)) {
      next
    }

    df2 <- current_row

  ################################################3

    # Get the column names that start with "Hcbh" followed by a number
    hcbh_cols <- grep("^cbh[0-9]+$", names(df2), value = TRUE)

    # Initialize an empty vector to store the column indices
    hcbh_cols_numeric <- numeric()

    # Iterate through each of the "Hcbh" columns
    for(col in hcbh_cols){
      # Check if the column contains only numeric values
      if(all(sapply(df2[[col]], is.numeric))){
        # If it does, add its index to the vector
        hcbh_cols_numeric <- c(hcbh_cols_numeric, which(names(df2) == col))
      }
    }

    ###################################3

    # Remove columns with only NA values
    df2 <- df2[, colSums(!is.na(df2)) > 0]

    dist_cols <- grep("^dist", names(df2), value = TRUE)
    if (length(dist_cols)==0) {
      df2$dist0<- 0
    }

    # Find columns starting with "Hdist"
    hdist_cols <- grep("^Hdist", names(df2), value = TRUE)
    hdist_vals <- df2 %>% dplyr::select(all_of(hdist_cols))

    # Find columns starting with "cbh"
    hcbh_cols <- grep("^cbh", names(df2), value = TRUE)
    hcbh_vals <- df2 %>% dplyr::select(all_of(hcbh_cols))

    # Check if the first cbh value is 0 and replace it with min_height
    if (hcbh_vals[[1]] <= min_height) {
      hcbh_vals[[1]] <- min_height
      df2$cbh1 <- min_height
    }
    hcbh_vals <- df2 %>% dplyr::select(all_of(hcbh_cols))


    #####################################################
    hdist_cols <- (grep("^Hdist", names(df2), value = TRUE))
    hdist_vals <- df2[, hdist_cols]

    hcbh_cols <- grep("^cbh", names(df2), value = TRUE)
    hcbh_vals <- df2[, hcbh_cols]

    dist_cols <- grep("^dist", names(df2), value = TRUE)
    dist_vals <- df2[, dist_cols]

    hdepth_cols <- grep("^Hdepth", names(df2), value = TRUE)
    hdepth_vals <- df2[, hdepth_cols]

    depth_cols <- grep("^depth", names(df2), value = TRUE)
    depth_vals <- df2[, depth_cols]

    # Convert any dist values of 0 to 1
   # dist_vals[dist_vals == 0] <- number_steps


    # Function to extract numeric part and order column names correctly
    order_columns <- function(cols) {
      # Extract numeric part from the column names
      num_part <- as.numeric(sub("^[^0-9]+", "", cols))
      # Order the columns based on the numeric part
      ordered_cols <- cols[order(num_part)]
      return(ordered_cols)
    }


    update_hcbh_values <- function(df, number_steps) {
      hcbh_cols <- grep("^cbh", names(df), value = TRUE)
      dist_cols <- grep("^dist", names(df), value = TRUE)
      hdist_cols <- grep("^Hdist", names(df), value = TRUE)

      # #print column names for debugging
      #print(paste("hcbh_cols:", paste(hcbh_cols, collapse = ", ")))
      #print(paste("dist_cols:", paste(dist_cols, collapse = ", ")))
      #print(paste("hdist_cols:", paste(hdist_cols, collapse = ", ")))

      # Check if the first Hdist column value is less than the first hcbh column value
      first_Hdist_col <- hdist_cols[1]
      first_hcbh_col <- hcbh_cols[1]
      min_height <- 0  # Define min_height if it is not defined elsewhere

      if (df[, first_Hdist_col] <= df[, first_hcbh_col] && df[, first_Hdist_col] > min_height) {
        start_index <- 1  # Start from index 1
      } else {
        start_index <- 2  # Start from Hdist2
      }

      # Define new_Hcbh using cbh columns
      new_Hcbh <- df[, hcbh_cols, drop = FALSE]

      # Ensure hdist_cols and dist_cols have the same length as hcbh_cols
      if (length(hdist_cols) > length(hcbh_cols)) {
        hdist_cols <- hdist_cols[1:length(hcbh_cols)]
        dist_cols <- dist_cols[1:length(hcbh_cols)]
      }

      for (i in start_index:length(hdist_cols)) {
        current_hdist_col <- hdist_cols[i]

        # Debug: #print the current hdist column being accessed
        #print(paste("Accessing column:", current_hdist_col))

        current_hdist_value <- df[, current_hdist_col]

        # Find hcbh columns that are greater than or equal to current_hdist_value
        relevant_hcbh_cols <- hcbh_cols[sapply(df[, hcbh_cols, drop = FALSE], function(x) any(x >= current_hdist_value))]

        if (length(relevant_hcbh_cols) > 0) {
          corresponding_dist_col <- dist_cols[i]
          current_dist <- df[, corresponding_dist_col]

          # Update the corresponding Hcbh column if dist is greater than number_steps
          match_found <- FALSE
          for (hcbh_col in relevant_hcbh_cols) {
            if (df[hcbh_col] >= current_hdist_value && current_dist > number_steps) {
              new_Hcbh[[hcbh_col]] <- df[[hcbh_col]]  # Update the corresponding Hcbh column
              match_found <- TRUE
              break
            }
          }

          # If no match found, use the last updated hcbh column
          if (!match_found && i > 1) {
            new_Hcbh[[hcbh_cols[i]]] <- new_Hcbh[[hcbh_cols[i - 1]]]  # Use the last updated Hcbh column
          }
        } else if (i > 1) {
          new_Hcbh[[hcbh_cols[i]]] <- new_Hcbh[[hcbh_cols[i - 1]]]  # Use the last updated Hcbh column
        }
      }

      # Update the original dataframe with the new Hcbh values
      for (col in hcbh_cols) {
        if (col %in% names(df)) {
          df[[col]] <- new_Hcbh[[col]]
        }
      }

      return(df)
    }

    # Apply the function to the dataframe
    updated_df2 <- update_hcbh_values(df2, number_steps = number_steps)



    #############################
    # Identify cbh columns with suffixes greater than the last Hdist suffix

    # Get the names of columns starting with "Hdist"
    hdist_cols <- grep("^Hdist", names(updated_df2), value = TRUE)

    # Extract the numeric suffix from the last Hdist column
    last_hdist_suffix <- as.numeric(gsub("Hdist", "", tail(hdist_cols, 1)))

    # Get the names of columns starting with "cbh" with numeric suffixes
    cbh_cols <- grep("^cbh\\d+$", names(updated_df2), value = TRUE)

    # Extract the numeric suffixes from the cbh columns
    cbh_suffixes <- as.numeric(gsub("cbh", "", cbh_cols))

    # Identify cbh columns with suffixes greater than the last Hdist suffix
    cbh_cols_remove <- cbh_cols[cbh_suffixes > last_hdist_suffix]

    # Remove the identified cbh columns from the dataframe
    updated_df2 <- updated_df2[, !names(updated_df2) %in% cbh_cols_remove]

    ###########################  CORRECT HDIST ####################

    # Find the number of Hdist columns
    num_columns <- grep("^Hdist", names(updated_df2), value = TRUE)

    if (length(num_columns) > 2){
    # Iterate through each Hdist column except the first one
    for (i in 2:length(num_columns)) {
      Hdist_col <- num_columns[i]
      Hcbh_col <- sub("Hdist", "cbh", Hdist_col)

      # Replace Hdisti with previous Hdist value if Hdisti > Hcbhi
      updated_df2[[Hdist_col]] <- ifelse(
        updated_df2[[Hdist_col]] > updated_df2[[Hcbh_col]],
        ifelse(
          is.na(updated_df2[[num_columns[i - 1]]]),
          updated_df2[[Hdist_col]],
          updated_df2[[num_columns[i - 1]]]
        ),
        updated_df2[[Hdist_col]]
      )
    }
    }

       #################################################3

    if (exists("updated_df2")) {
      effective_fbh <- updated_df2
    } else {
      effective_fbh <- df2
    }

    effective_fbh <- effective_fbh[, colSums(!is.na(effective_fbh)) > 0]

    rename_cbh_columns <- function(df) {
      cbh_cols <- grep("^cbh", names(df), value = TRUE)
      new_names <- sub("^cbh", "Hcbh", cbh_cols)
      names(df)[names(df) %in% cbh_cols] <- new_names
      return(df)
    }

    # Apply the renaming function to the dataframe
    effective_fbh <- rename_cbh_columns(effective_fbh)

    # Extract columns that start with 'Hcbh'
    Hcbh_cols <- grep("^Hcbh", names(effective_fbh), value = TRUE)

    # Check the conditions
    first_value <- effective_fbh[1, Hcbh_cols[1]]
    subsequent_values <- effective_fbh[1, Hcbh_cols[-1]]

    condition1 <- any(first_value == subsequent_values + number_steps)
    condition2 <- length(unique(as.numeric(subsequent_values))) == 1

    if(condition1 && condition2) {
      # Set all Hcbh column values to the first Hcbh column value
      effective_fbh[1, Hcbh_cols] <- first_value
    }


    ##################  rename columns

      # Extract unique prefixes
    prefixes <- unique(gsub("([a-zA-Z]+).*", "\\1", names(effective_fbh)))

    # Rename the columns based on the extracted prefixes
    for (prefix in prefixes) {
      # Identify columns with the current prefix
      cols <- grep(paste0("^", prefix), names(effective_fbh))

      # Generate new column names with consecutive suffixes
      new_names <- paste0(prefix, 1:length(cols))

      # Assign new names to the columns
      names(effective_fbh)[cols] <- new_names
    }

    effective_fbh$treeID <- effective_fbh$treeID2
    effective_fbh$treeID1 <- sub("_.*", "", effective_fbh$treeID)
    effective_fbh$treeID1 <- as.numeric(effective_fbh$treeID1)
    effective_fbh$max_height <- effective_fbh$max1

    effective_fbh <- effective_fbh %>%
      dplyr::select(-treeID2, -max1)


      if(!exists("effective_fbh")) {
      next  # skip the current iteration if effective_fbh is not generated
    }
  }

  return(effective_fbh)
}
