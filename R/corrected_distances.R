#' Effective Distances between fuel layers
#' @description
#' This function recalculates the distance between fuel layers after considering distances greater than any number of height bin steps.
#' @usage
#' get_effective_gap(effective_depth, number_steps = 1, min_height= 1.5, verbose=TRUE)
#' @param effective_depth
#' Tree metrics with the recalculated depth values after considering distances greater than the actual height bin step (output of [get_real_depths()]
#' function). An object of the class data frame.
#' @param number_steps Numeric value for the number of height bin steps that can be jumped to reshape fuels layers.
#' @param min_height Numeric value for the actual minimum base height (in meters).
#' @param verbose Logical, indicating whether to display informational messages (default is TRUE).
#' @return
#' A data frame giving the effective distances (> any number of steps) between consecutive fuel layers.
#' @author Olga Viedma, Carlos Silva, JM Moreno and A.T. Hudak
#'
#' @details
#' List of tree metrics:
#' \itemize{
#'   \item treeID: tree ID with strings and numeric values
#'   \item treeID1: tree ID with only numeric values
#'   \item dist: Distance between consecutive fuel layers (m)
#'   \item dptf: Depth of fuel layers (m) after considering distances greater than the actual height bin step
#'   \item effdist: Effective distance between consecutive fuel layers (m) after considering distances greater than any number of steps
#'   \item Hcbh: Base height of each fuel separated by a distance greater than the certain number of steps
#'   \item Hdist: Height of the distance (> any number of steps) between consecutive fuel layers (m)
#'   \item Hdptf: Height of the depth of fuel layers (m) after considering distances greater than the actual step
#'   \item max_height: Maximum height of the tree
#' }
#'
#' @examples
#' library(magrittr)
#' library(stringr)
#' library(dplyr)
#'
#' # Before running this example, make sure to run get_real_depths().
#' if (interactive()) {
#' effective_depth <- get_real_depths()
#' LadderFuelsR::effective_depth$treeID <- factor(LadderFuelsR::effective_depth$treeID)
#'
#' trees_name1 <- as.character(effective_depth$treeID)
#' trees_name2 <- factor(unique(trees_name1))
#'
#' corr_distance_metrics_list <- list()
#'
#' for (i in levels(trees_name2)) {
#' tree1 <- effective_depth |> dplyr::filter(treeID == i)
#' corr_distance_metrics <- get_effective_gap(tree1, number_steps = 1, min_height= 1.5, verbose=TRUE)
#' corr_distance_metrics_list[[i]] <- corr_distance_metrics
#' }
#'
#' # Combine the individual data frames
#' effective_distances <- dplyr::bind_rows(corr_distance_metrics_list)
#'
#' # Get original column names
#' original_column_names <- colnames(effective_distances)
#'
#' # Specify prefixes
#' desired_order <- c("treeID", "Hcbh", "dptf","effdist","dist", "Hdist", "Hdptf", "max_")
#'
#'# Identify unique prefixes
#' prefixes <- unique(sub("^([a-zA-Z]+).*", "\\1", original_column_names))
#' # Initialize vector to store new order
#' new_order <- c()
#'
#' # Loop over desired order of prefixes
#' for (prefix in desired_order) {
#'  # Find column names matching the current prefix
#' matching_columns <- grep(paste0("^", prefix), original_column_names, value = TRUE)
#' # Append to the new order
#' new_order <- c(new_order, matching_columns)
#' }
#' effective_distances <- effective_distances[, new_order]
#' }
#' @importFrom dplyr select_if group_by summarise summarize mutate arrange rename rename_with filter slice slice_tail ungroup distinct
#' across matches row_number all_of vars last
#' @importFrom segmented segmented seg.control
#' @importFrom magrittr %>%
#' @importFrom stats ave dist lm na.omit predict quantile setNames smooth.spline
#' @importFrom utils tail
#' @importFrom tidyselect starts_with everything one_of
#' @importFrom stringr str_extract str_match str_detect str_remove_all
#' @importFrom tibble tibble
#' @importFrom tidyr pivot_longer fill
#' @importFrom gdata startsWith
#' @importFrom ggplot2 aes geom_line geom_path geom_point geom_polygon geom_text geom_vline ggtitle coord_flip theme_bw
#' theme element_text xlab ylab ggplot
#' @seealso \code{\link{get_real_depths}}
#' @export
get_effective_gap <- function(effective_depth, number_steps = 1, min_height= 1.5, verbose = TRUE) {

  if(min_height==0){
    min_height <-0.5
  }

  df <- effective_depth

if (verbose) {
  message("Unique treeIDs:", paste(unique(df$treeID), collapse = ", "))
}

  df5b <- df[, colSums(!is.na(df)) > 0]

  # Identify columns starting with 'Hcbh' and remove duplicates
  hcbh_cols <- grep("^Hcbh", colnames(df5b), value = TRUE)
  # Identify which columns have duplicated values
  duplicated_values <- duplicated(as.numeric(df5b[1, hcbh_cols]))
  # Set all but the first occurrence to NA
  df5b[1, hcbh_cols[duplicated_values]] <- NA

  df5b <- df5b[, colSums(!is.na(df5b)) > 0]

  # Extract unique prefixes
  prefixes <- unique(gsub("([a-zA-Z]+).*", "\\1", names(df5b)))
  # Rename the columns based on the extracted prefixes
  for (prefix in prefixes) {
    # Identify columns with the current prefix
    cols <- grep(paste0("^", prefix), names(df5b))
    # Generate new column names with consecutive suffixes
    new_names <- paste0(prefix, 1:length(cols))
    # Assign new names to the columns
    names(df5b)[cols] <- new_names
  }

  hcbh_cols <- sort(grep("^Hcbh", names(df5b), value = TRUE))
  hcbh_vals <- df5b[1, hcbh_cols, drop = FALSE]

  hdist_cols <- sort(grep("^Hdist", names(df5b), value = TRUE))
  hdist_vals <- df5b[1, hdist_cols, drop = FALSE]

  dist_cols <- sort(grep("^dist", names(df5b), value = TRUE))
  dist_vals <- df5b[1, dist_cols, drop = FALSE]


  #############################################################33
  # Check if there's only one dist column
  dist_cols <- grep("^dist", names(df5b))
  if (length(dist_cols) == 1) {
    # If there's only one dist column, set effective_dist1 to dist1
    effective_dist1 <- data.frame(effdist = df5b[1, dist_cols])
    colnames(effective_dist1) <- "effdist1"
  } else {
    # Initialize effdist with zeros
    effdist <- numeric(ncol(df5b[grep("dist", names(df5b))]))

    hcbh_vals <- df5b[grep("^Hcbh", names(df5b))]
    dist_vals <- df5b[grep("^dist", names(df5b))]

    if (ncol(hcbh_vals) > 1) {
      # Adjust the iteration limit based on the number of columns in hcbh_vals
      end_col <- ncol(hcbh_vals) - 1

      # Start iterating from the first column
      for (i in 1:end_col) {
        # Check if neither of the values is NA before comparing
        if (!is.na(hcbh_vals[1, i]) && !is.na(hcbh_vals[1, i + 1])) {
          # Compare consecutive Hcbh columns
          if (hcbh_vals[1, i] != hcbh_vals[1, i + 1]) {
            if (i == 1 && dist_vals[1, i] <= number_steps) {
              effdist[i] <- 0  # Set to 0 if the first dist value <= number_steps
              effdist[i + 1] <- dist_vals[1, i + 1]  # Proceed with the next dist value
            } else if (dist_vals[1, i] > number_steps && dist_vals[1, i + 1] > number_steps) {
              effdist[i] <- dist_vals[1, i]
              effdist[i + 1] <- dist_vals[1, i + 1]
            }
          }
        }
      }

      # If only one dist and Hcbh column, assign the dist value to effdist if dist > 1
      if (ncol(dist_vals) == 1 && dist_vals[1, 1] > number_steps) {
        effdist[1] <- dist_vals[1, 1]
      }
    }

    # Remove trailing zeros and keep only the relevant distances
    effdist <- effdist[effdist > 0]

    # Get the numeric suffixes of the original dist columns analyzed
    dist_suffixes <- which(dist_vals > number_steps)
    dist_suffixes <- gsub("dist", "", names(dist_vals)[dist_suffixes])

    # Create a data frame for the effective distances
    effective_dist1 <- data.frame(t(effdist))

    # Name the columns of effective_dist1 with "effdist" and the numeric suffixes
    colnames(effective_dist1) <- paste0("effdist", dist_suffixes)
  }

    ################################################

  if (exists("effective_dist1") && ncol(effective_dist1) > 0) {
    df6 <- cbind.data.frame(df5b, effective_dist1)
  }

  if (length(effective_dist1) == 0 && ncol(hcbh_vals) ==1 && ncol(dist_vals) ==1) {

    effdist1 <-data.frame(df5b$dist1)
    names(effdist1)<-"effdist1"
    effective_dist1 <- effdist1
    df6 <- cbind.data.frame(df5b, effective_dist1)
  }

  if(!exists ("df6")){
    df6<-df5b
  }

  # Check if dist1 <= number_steps, if so, set effdist1 to 0
  if (df6$dist1 < number_steps) {
    df6$effdist1 <- 0
  }

  if(!"effdist1" %in% colnames(df6) && "dist1" %in% colnames(df6)){
    df6$effdist1 <- ceiling(df6$Hcbh1 - min_height)

  }

  ######################################################

  # Remove columns with only NA values
  df6a <- df6[, colSums(!is.na(df6)) > 0]
  df6a <- df6a[, order(colnames(df6a))]

  # Extract unique prefixes
  prefixes <- unique(gsub("([a-zA-Z]+).*", "\\1", names(df6a)))

  # Rename the columns based on the extracted prefixes
  for (prefix in prefixes) {
    # Identify columns with the current prefix
    cols <- grep(paste0("^", prefix), names(df6a))

    # Generate new column names with consecutive suffixes
    new_names <- paste0(prefix, 1:length(cols))

    # Assign new names to the columns
    names(df6a)[cols] <- new_names
  }

  # Adjust the regex to only match effdist columns with a numeric suffix
  if (sum(grepl("Hdptf[0-9]+$", colnames(df6a))) > sum(grepl("Hcbh[0-9]+$", colnames(df6a)))) {
    # Find the last 'Hdptf' column with numeric suffix
    last_Hdptf_col <- tail(grep("Hdptf[0-9]+$", colnames(df6a)), 1)
    # Drop the column
    df6a <- df6a[,-last_Hdptf_col]
  }


  # Check if Hdist column exists and create it if it doesn't
  if (!"Hdist1" %in% colnames(df6a)) {
    df6a <- df6a %>%
      mutate(Hdist1 = min_height)
  }

  # Extracting column names
  col_names <- names(df6a)

  # Finding indices of relevant columns
  hdist_indices <- grep("^Hdist[0-9]*$", col_names)
  hdptf_indices <- grep("^Hdptf[0-9]*$", col_names)
  hcbh_indices <- grep("^Hcbh[0-9]*$", col_names)
  effdist_indices <- grep("^effdist[0-9]*$", col_names)
  dist_indices <- grep("^dist[0-9]*$", col_names)

  # Convert effdist columns to numeric
  df6a[, effdist_indices] <- lapply(df6a[, effdist_indices], as.numeric)


  # Check condition and modify columns if TRUE
  if (df6a$Hcbh1 == min_height && (df6a$Hdist1 > df6a$Hcbh1)) {
    # Create a new data frame with updated column names
    df6a_updated <- df6a

    # Create a list to store new column names
    new_col_names <- col_names

    # Rename Hdist columns with incremented suffixes
    for (i in hdist_indices) {
      old_name <- col_names[i]
      # Extract the suffix, handle NA case
      old_suffix <- as.numeric(sub("Hdist", "", old_name))
      new_suffix <-  old_suffix + 1
      new_name <- paste0("Hdist", new_suffix)

      # Rename column if the new name is valid (not NA)
      if (!is.na(new_suffix)) {
        new_col_names[i] <- new_name
      }
    }

    # Remove "Hdist" if it exists in the new column names list
    if ("Hdist" %in% new_col_names) {
      new_col_names <- new_col_names[new_col_names != "Hdist"]
    }

    colnames(df6a_updated) <- new_col_names

    new_col_names1 <- colnames(df6a_updated)

    # Rename Hdist columns with incremented suffixes
    for (i in dist_indices) {
      old_name <- colnames(df6a_updated)[i]
      # Extract the suffix, handle NA case
      old_suffix <- as.numeric(sub("dist", "", old_name))
      new_suffix <- ifelse(is.na(old_suffix), NA, old_suffix + 1)
      new_name <- paste0("dist", new_suffix)

      # Rename column if the new name is valid (not NA)
      if (!is.na(new_suffix)) {
        new_col_names1[i] <- new_name
      }
    }

    colnames(df6a_updated) <- new_col_names1

     new_col_names2 <-  colnames(df6a_updated)

    # Rename Hdist columns with incremented suffixes
    for (i in effdist_indices) {
      old_name <-  colnames(df6a_updated)[i]
      # Extract the suffix, handle NA case
      old_suffix <- as.numeric(sub("effdist", "", old_name))
      new_suffix <- ifelse(is.na(old_suffix), NA, old_suffix + 1)
      new_name <- paste0("effdist", new_suffix)

      # Rename column if the new name is valid (not NA)
      if (!is.na(new_suffix)) {
        new_col_names2[i] <- new_name
      }
    }

    colnames(df6a_updated) <- new_col_names2


    # Check if 'dist' exists in the new column names list
    if (!("effdist1" %in%  colnames(df6a_updated) )) {
      if(df6a_updated$Hcbh1 == min_height) {
      df6a <- df6a_updated %>%
        mutate(effdist1 =0)
    }}

    # Check if 'dist' exists in the new column names list
    if (!("Hdist1" %in%  colnames(df6a))) {
      if(df6a_updated$Hcbh1 == min_height) {
      df6a <- df6a %>%
        mutate(Hdist1 = min_height)
    }}

    # Check if 'dist' exists in the new column names list
    if (!("dist1" %in% colnames(df6a))) {
      if(df6a_updated$Hcbh1 == min_height) {
      df6a <- df6a %>%
        mutate(dist1 = 0)
    }}
}


  if (length(effdist_indices) > length(hdist_indices)) {
    cols_to_remove <- c(tail(effdist_indices, 1))
    df6a <- df6a[, -cols_to_remove, drop = FALSE]
  }


 ######################################################3

  max_height<-data.frame(df5b$max1)
  names(max_height)<-"max_height"

  if((!"max_height" %in% colnames(df6a))) {
    df6f<-data.frame(df6a, max_height)
  }

  treeID <-df6f$treeID
  treeID1 <-df6f$treeID1

  cols_to_exclude <- grep("treeID|treeID1)", names(df6f), value = TRUE)
  trees<-df6f[ , (names(df6f) %in% cols_to_exclude)]
  df6f <- df6f[ , !(names(df6f) %in% cols_to_exclude)]

  df6f <- data.frame (trees,df6f)

  if("treeID2" %in% colnames (df6f)){
    treeID2 <-df6f$treeID2
    df6f <- df6f %>%
      dplyr::rename(
        treeID1= treeID2,
        treeID = treeID1)
  }

  cols_to_exclude1 <- grep("max1", names(df6f), value = TRUE)
  df6f <- df6f[ , !(names(df6f) %in% cols_to_exclude1)]

  effective_distances<-df6f

  # Remove list attributes from columns
  effective_distances[] <- lapply(effective_distances, function(x) {
    if (is.list(attributes(x)))
      attributes(x) <- attributes(x)[!names(attributes(x)) %in% c("dimnames")]
    return(x)
  })

  # Loop through each column
  for (col in names(effective_distances)) {
    # Convert each column to a vector
    effective_distances[[col]] <- unlist(effective_distances[[col]])
  }

  # Identify columns with matrix-like structure
  matrix_columns <- sapply(effective_distances, function(x) is.matrix(x) && nrow(x) > 1)

  # Extract the numeric values from matrix columns
  effective_distances[matrix_columns] <- lapply(effective_distances[matrix_columns], function(x) as.numeric(x[, 1]))

  # Convert data frame to remove list attributes
  effective_distances <- data.frame(effective_distances)

  # Get the columns that start with "treeID"
  treeID_columns <- grep("^treeID", names(effective_distances), value = TRUE)

  # Convert all variables to numeric except "treeID" columns
  effective_distances[, !names(effective_distances) %in% treeID_columns] <-
    lapply(effective_distances[, !names(effective_distances) %in% treeID_columns], as.numeric)



  return(effective_distances)
}

