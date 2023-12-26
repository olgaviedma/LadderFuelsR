#' Effective Distances and CBHs based on maximum and last distance
#' @description
#' This function recalculates the distance between fuel layers after removing distances = 1 m and determines the CBH based on the fuel layers
#' with the highest and the last distance.
#' @usage
#' get_effective_gap(effective_depth, verbose=TRUE)
#' @param effective_depth
#' Tree metrics with the recalculated depth values considering distances > 1 m (output of [get_real_depths()] function).
#' An object of the class data frame.
#' @param verbose Logical, indicating whether to display informational messages (default is TRUE).
#' @return
#' A data frame giving the effective distances (> 1 m) between consecutive fuel layers, identifying the Canopy Base Height (CBH) of the fuel layer
#' at the maximum distance and at the last distance.
#' @author
#' Olga Viedma, Carlos Silva and JM Moreno
#'
#' @details
#' List of tree metrics:
#' \itemize{
#'   \item treeID: tree ID with strings and numeric values
#'   \item treeID1: tree ID with only numeric values
#'   \item dist: Distance between consecutive fuel layers (m)
#'   \item Hdist: Height of the distance between consecutive fuel layers (m)
#'   \item Hcbh: Height of the base of each fuel layer (m)
#'   \item effdist: Effective distance between consecutive fuel layers (m) (> 1 m)
#'   \item dptf: Depth of fuel layers (m) after removing distances equal to 1 m
#'   \item Hdptf: Height of the depth of fuel layers (m) after removing distances equal to 1 m
#'   \item max_Hcbh: Canopy base height of the segmented tree based on the maximum distance found in its profile
#'   \item last_Hcbh: Canopy base height of the segmented tree based on the last distance found in its profile
#'   \item max_height: Maximum height of the tree profile
#'   \item max_: Values of distance and fuel depth and their corresponding heights at the maximum distance found in the tree profile
#'   \item last_: Values of distance and fuel depth and their corresponding heights at the last distance found in its profile
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
#' corr_distance_metrics <- get_effective_gap(tree1, verbose=TRUE)
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
#' desired_order <- c("treeID", "Hcbh", "dptf","effdist","dist", "Hdist", "Hdptf", "max_","last_")
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
#' @importFrom stringr str_extract str_match str_detect
#' @importFrom tibble tibble
#' @importFrom tidyr pivot_longer fill
#' @importFrom gdata startsWith
#' @importFrom ggplot2 aes geom_line geom_path geom_point geom_polygon geom_text geom_vline ggtitle coord_flip theme_bw
#' theme element_text xlab ylab ggplot
#' @seealso \code{\link{get_real_depths}}
#' @export
get_effective_gap <- function(effective_depth, verbose = TRUE) {

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

  if (length(dist_cols) == 0) {
    df5b$dist1 <- 0
  }


  # Initialize effdist with zeros
  effdist <- numeric(ncol(df5b[grep("Hcbh", names(df5b))]))

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
        if (hcbh_vals[1, i] != hcbh_vals[1, i + 1] && dist_vals[1, i] > 1) {
          effdist[i] <- dist_vals[1, i]
        }
      }
    }

    # For the last column, take the dist value if there are no more paired Hcbh values
    if (dist_vals[1, ncol(dist_vals)] > 1) {
      effdist[ncol(hcbh_vals)] <- dist_vals[1, ncol(dist_vals)]
    }
  } else {
    # If only one dist and Hcbh column, assign the dist value to effdist if dist > 1
    if (dist_vals[1, 1] > 1) {
      effdist[1] <- dist_vals[1, 1]
    }
  }

  effective_dist1 <- data.frame(t(effdist))
  # Assign column names
  colnames(effective_dist1) <- paste0("effdist", 1:ncol(effective_dist1))

  df6 <- cbind.data.frame(df5b, effective_dist1)


  # 1. Remove `effdist` values that are 0
  effdist_cols <- grep("^effdist", names(df6), value = TRUE)
  df6 <- df6[, !(names(df6) %in% effdist_cols[df6[effdist_cols] == 0])]

  Hcbh_cols <- grep("^Hcbh", names(df6), value = TRUE)
  hcbh_data <- df6[Hcbh_cols]

  # Transpose, remove duplicates, and transpose back
  hcbh_data_cleaned <- t(unique(t(hcbh_data)))

  original_colnames <- colnames(hcbh_data)
  cleaned_colnames <- colnames(hcbh_data_cleaned)

  # Determine which columns to replace
  cols_to_replace <- original_colnames[original_colnames %in% cleaned_colnames]

  # Replace in the original df
  df6[cols_to_replace] <- hcbh_data_cleaned

  df6 <- df6[, !names(df6) %in% original_colnames[!original_colnames %in% cleaned_colnames]]

  # Select "Hdepth" columns
  Hdepth_cols <- names(df6)[str_detect(names(df6), "^Hdptf")]

  # Replace 0s with NAs in "Hdepth" columns
  df6[Hdepth_cols] <- lapply(df6[Hdepth_cols], function(x) ifelse(x == 0, NA, x))

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


  # Extracting column names
  col_names <- names(df6a)

  # Finding indices of Hdist, Hdptf, and effdist columns
  hdist_indices <- grep("^Hdist[0-9]*$", col_names)
  hdptf_indices <- grep("^Hdptf[0-9]*$", col_names)
  hcbh_indices <- grep("^Hcbh[0-9]*$", col_names)
  effdist_indices <- grep("^effdist[0-9]*$", col_names)
  dist_indices <- grep("^dist[0-9]*$", col_names)


  # Remove the columns
  if (df6a$Hcbh1 == 1.5 && length(hdist_indices) == length(hcbh_indices)) {
    cols_to_remove <- c(tail(effdist_indices, 1), tail(hdist_indices, 1), tail(dist_indices, 1))
    df6a <- df6a[, -cols_to_remove, drop = FALSE]
  }

  if (df6a$Hcbh1 > 1.5 && length(hdist_indices) > length(hcbh_indices)) {
    cols_to_remove <- c(tail(effdist_indices, 1), tail(hdist_indices, 1), tail(dist_indices, 1))
    df6a <- df6a[, -cols_to_remove, drop = FALSE]
  }

  if (length(effdist_indices) > length(hdist_indices)) {
    cols_to_remove <- c(tail(effdist_indices, 1))
    df6a <- df6a[, -cols_to_remove, drop = FALSE]
  }

  #############################################3333
  effdist_cols <- grep("^effdist", names(df6a), value = TRUE)

  if(length(effdist_cols) == 0) {

    df6f<- df6a

    hcbh_cols <- sort(names(df6f)[grep("Hcbh", names(df6f))])
    dptf_cols <- sort(names(df6f)[grep("^dptf", names(df6f))])
    hdptf_cols <- sort(names(df6f)[grep("^Hdptf", names(df6f))])


    if (hcbh_cols %in% colnames(df6f) & dptf_cols %in% colnames(df6f)  & hdptf_cols %in% colnames(df6f)) {
      last_df1 <- data.frame(
        last_Hcbh = df6f[[hcbh_cols]],
        last_dptf = df6f[[dptf_cols]],
        last_Hdptf = df6f[[hdptf_cols]])

      max_df1 <- data.frame(
        max_Hcbh = df6f[[hcbh_cols]],
        max_dptf = df6f[[dptf_cols]],
        max_Hdptf = df6f[[hdptf_cols]])

      names(last_df1) <-  c("last_Hcbh", "last_dptf","last_Hdptf")  # And update here
      names(max_df1) <-  c("max_Hcbh", "max_dptf","max_Hdptf")  # And update here

      df6f<-cbind (df6f,max_df1,last_df1)

      treeID1<-df6f$treeID1

      df6f <- df6f %>%
        dplyr::rename(
          treeID= treeID2,
          treeID1 = treeID1)

      max_height<-data.frame(df5b$max1)
      names(max_height)<-"max_height"

      if(!("max_height" %in% colnames(df6f))) {

        df6f<-data.frame(df6f, max_height)
      }
    }
  }


  ######################### THE HCBH AND HDIST CORRESPONDING TO THE MAXIMUM EFFECTIVE DISTANCE BY POSITION ################################

  if(length(effdist_cols) > 0) {

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

    #########################################
    effdist_cols <- sort(grep("^effdist", names(df6a), value = TRUE))
    effdist_vals <- df6a[1, effdist_cols, drop = FALSE]


    if(length(effdist_cols) > 0 && any(df6a[, effdist_cols] > 1)) {

      # First, select the 'effdist' columns
      effdist_cols <- names(df6a)[str_detect(names(df6a), "^effdist")]
      # Extract effdist values from the first row
      effdist_values <- df6a[1, grep("effdist", names(df6a))]
      # Find the index of the column with max effdist
      max_effdist_col_index <- tail(which(effdist_values == max(effdist_values)), n=1)

      # Define the column names for Hcbh, Hdist, dptf, and Hdptf
      hcbh_cols <- sort(names(df6a)[grep("Hcbh", names(df6a))])
      hdist_cols <- sort(names(df6a)[grep("Hdist", names(df6a))])
      dptf_cols <- sort(names(df6a)[grep("^dptf", names(df6a))])
      hdptf_cols <- sort(names(df6a)[grep("^Hdptf", names(df6a))])

      if (df6a$effdist1 == 1 && all(sapply(df6a[effdist_cols], function(x) x[1] == df6a$effdist1[1]))) {

        max_df <- df6a %>%
          dplyr::rename_with(~ ifelse(. %in% c("Hcbh1", "dptf1", "Hdist1", "Hdptf1", "effdist1"), paste0("max_", str_remove(., "\\d+$")), .))

      }

      if (df6a$effdist1[1] == 1 & all(df6a[effdist_cols][-1] != df6a$effdist1[1])) {
        df6a$effdist1 <-NULL
      }

      #### rename ##############3
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

      effdist_cols <- sort(grep("^effdist", names(df6a), value = TRUE))
      effdist_vals <- df6a[1, effdist_cols, drop = FALSE]
      hcbh_cols <- sort(names(df6a)[grep("Hcbh", names(df6a))])
      hdist_cols <- sort(names(df6a)[grep("Hdist", names(df6a))])
      dptf_cols <- sort(names(df6a)[grep("^dptf", names(df6a))])
      hdptf_cols <- sort(names(df6a)[grep("^Hdptf", names(df6a))])
      #####################################################3

      if (length(effdist_cols) == 1 & df6a[[hcbh_cols[1]]][1] == 1.5 & length(hdist_cols) > 1) {
        # Create a new dataframe with the required columns
        df6ab <- df6a[,c( hcbh_cols[2], hdist_cols[2], dptf_cols[2], hdptf_cols[2], effdist_cols[1])]
        # Rename columns with prefix "max_" and remove suffix
        max_df <- df6ab %>% dplyr::rename_with(.fn = ~ paste0("max_", str_remove(., "\\d+$")))
      }

      if (length(effdist_cols) > 1 & df6a[[hcbh_cols[1]]][1] == 1.5 && length(hdist_cols) >= 1 && !exists("max_df")) {

        suffix <- max_effdist_col_index
        suffix1 <- as.character(as.numeric(suffix) +1)
        #suffix2 <- as.character(as.numeric(suffix) +2)

        # Get the suffixes of the columns
        hcbh_suffixes <- str_extract(hcbh_cols, "\\d+$")
        hdist_suffixes <- str_extract(hdist_cols, "\\d+$")
        dptf_suffixes <- str_extract(dptf_cols, "\\d+$")
        hdptf_suffixes <- str_extract(hdptf_cols, "\\d+$")
        effdist_suffixes <- str_extract(effdist_cols, "\\d+$")

        # Get the columns whose suffixes are greater or equal to the max effdist suffix
        hcbh_col <- hcbh_cols[which.max(as.numeric(hcbh_suffixes) >= as.numeric(suffix1))]
        hdist_col <- hdist_cols[which.max(as.numeric(hdist_suffixes) >= as.numeric(suffix))]
        dptf_col <- dptf_cols[which.max(as.numeric(dptf_suffixes) >= as.numeric(suffix1))]
        hdptf_col <- hdptf_cols[which.max(as.numeric(hdptf_suffixes) >= as.numeric(suffix1))]
        effdist_col <- effdist_cols[which.max(as.numeric(effdist_suffixes) >= as.numeric(suffix))]

        # Create a new dataframe with the required columns
        df6ab <- df6a[,c(hcbh_col, hdist_col, dptf_col, hdptf_col,effdist_col)]

        # Rename columns with prefix "max_" and remove suffix
        max_df <- df6ab %>% dplyr::rename_with(.fn = ~ paste0("max_", str_remove(., "\\d+$")))

      }

      if (length(effdist_cols) > 1 & df6a[[hcbh_cols[1]]][1] > 1.5 && length(hdist_cols) >= 1 && df6a[[dist_cols[1]]][1] > 1 && !exists("max_df")) {

        suffix <- max_effdist_col_index
        #suffix2 <- as.character(as.numeric(suffix) +2)

        # Get the suffixes of the columns
        hcbh_suffixes <- str_extract(hcbh_cols, "\\d+$")
        hdist_suffixes <- str_extract(hdist_cols, "\\d+$")
        dptf_suffixes <- str_extract(dptf_cols, "\\d+$")
        hdptf_suffixes <- str_extract(hdptf_cols, "\\d+$")
        effdist_suffixes <- str_extract(effdist_cols, "\\d+$")

        # Get the columns whose suffixes are greater or equal to the max effdist suffix
        hcbh_col <- hcbh_cols[which.max(as.numeric(hcbh_suffixes) >= as.numeric(suffix))]
        hdist_col <- hdist_cols[which.max(as.numeric(hdist_suffixes) >= as.numeric(suffix))]
        dptf_col <- dptf_cols[which.max(as.numeric(dptf_suffixes) >= as.numeric(suffix))]
        hdptf_col <- hdptf_cols[which.max(as.numeric(hdptf_suffixes) >= as.numeric(suffix))]
        effdist_col <- effdist_cols[which.max(as.numeric(effdist_suffixes) >= as.numeric(suffix))]

        # Create a new dataframe with the required columns
        df6ab <- df6a[,c(hcbh_col, hdist_col, dptf_col, hdptf_col,effdist_col)]

        # Rename columns with prefix "max_" and remove suffix
        max_df <- df6ab %>% dplyr::rename_with(.fn = ~ paste0("max_", str_remove(., "\\d+$")))
      }

      if (length(effdist_cols) > 1 & df6a[[hcbh_cols[1]]][1] > 1.5 && length(hdist_cols) >= 1 && df6a[[dist_cols[1]]][1] == 1 && !exists("max_df")) {

        suffix <- max_effdist_col_index
        suffix2 <- as.character(as.numeric(suffix) +1)

        # Get the suffixes of the columns
        hcbh_suffixes <- str_extract(hcbh_cols, "\\d+$")
        hdist_suffixes <- str_extract(hdist_cols, "\\d+$")
        dptf_suffixes <- str_extract(dptf_cols, "\\d+$")
        hdptf_suffixes <- str_extract(hdptf_cols, "\\d+$")
        effdist_suffixes <- str_extract(effdist_cols, "\\d+$")

        # Get the columns whose suffixes are greater or equal to the max effdist suffix
        hcbh_col <- hcbh_cols[which.max(as.numeric(hcbh_suffixes) >= as.numeric(suffix2))]
        hdist_col <- hdist_cols[which.max(as.numeric(hdist_suffixes) >= as.numeric(suffix2))]
        dptf_col <- dptf_cols[which.max(as.numeric(dptf_suffixes) >= as.numeric(suffix2))]
        hdptf_col <- hdptf_cols[which.max(as.numeric(hdptf_suffixes) >= as.numeric(suffix2))]
        effdist_col <- effdist_cols[which.max(as.numeric(effdist_suffixes) >= as.numeric(suffix))]

        # Create a new dataframe with the required columns
        df6ab <- df6a[,c(hcbh_col, hdist_col, dptf_col, hdptf_col,effdist_col)]

        # Rename columns with prefix "max_" and remove suffix
        max_df <- df6ab %>% dplyr::rename_with(.fn = ~ paste0("max_", str_remove(., "\\d+$")))
      }


    }

    effdist_cols <- sort(grep("^effdist", names(df6a), value = TRUE))
    effdist_vals <- df6a[1, effdist_cols, drop = FALSE]

    if(length(effdist_cols) == 1 | all(df6a[, effdist_cols] %in% c(0, 1)) && !exists("max_df")) { ### max-last Hdist

      # Define the column names for Hcbh, Hdist, dptf, and Hdptf

      hcbh_cols <- grep("^Hcbh", names(df6a), value=TRUE)
      last_hcbh_col <- hcbh_cols[length(hcbh_cols)]
      hdist_cols <- grep("^Hdist", names(df6a), value=TRUE)
      last_hdist_col <- hdist_cols[length(hdist_cols)]
      dptf_cols <- grep("^dptf", names(df6a), value=TRUE)
      last_dptf_col <- dptf_cols[length(dptf_cols)]
      hdptf_cols <- grep("^Hdptf", names(df6a), value=TRUE)
      last_hdptf_col <- hdptf_cols[length(hdptf_cols)]
      effdist_cols <- grep("^effdist", names(df6a), value=TRUE)
      last_effdist_col <- effdist_cols[length(effdist_cols)]

      max_df <- df6a %>%
        dplyr::select(max_Hcbh = {{last_hcbh_col}},
                      max_Hdptf = {{last_hdptf_col}},
                      max_dptf = {{last_dptf_col}},
                      max_Hdist = {{last_hdist_col}},
                      max_effdist = {{last_effdist_col}})
    }

    ##############################################################
    ####### the last Hcbh with numeric values ####################
    ######################################

    # Identify columns which only have zeros
    #cols_to_remove <- sapply(df6a, function(col) all(col == 0))
    # Replace NA values with FALSE in cols_to_remove
    #cols_to_remove[is.na(cols_to_remove)] <- FALSE
    # Remove those columns
    #df6a <- df6a[, !cols_to_remove]

    effdist_cols <- names(df6a)[str_detect(names(df6a), "^effdist")]
    effdist_values <- df6a[1, grep("effdist", names(df6a))]


    if(length(effdist_cols) > 0 && any(df6a[, effdist_cols] > 1)) {

      # Get the column names that start with "Hcbh" followed by a number
      hcbh_cols <- grep("^Hcbh[0-9]+$", names(df6a), value = TRUE)

      # Initialize an empty vector to store the column indices
      hcbh_cols_numeric <- numeric()

      # Iterate through each of the "Hcbh" columns
      for(col in hcbh_cols){
        # Check if the column contains only numeric values
        if(all(sapply(df6a[[col]], is.numeric))){
          # If it does, add its index to the vector
          hcbh_cols_numeric <- c(hcbh_cols_numeric, which(names(df6a) == col))
        }
      }

      # The last "Hcbh" column with numeric values will be:
      last_Hcbh_numeric_col <- names(df6a)[max(hcbh_cols_numeric)]
      last_Hcbh_numeric_val<-df6a[[last_Hcbh_numeric_col]]

      # Extract the suffix from the last_Hcbh_numeric_col
      suffix <- gsub("Hcbh", "", last_Hcbh_numeric_col)

      # Form the corresponding "Hdist" and "Hdepth" column names
      last_Hdist_col <- paste0("Hdist", suffix)
      last_Hdepth_col <- paste0("Hdptf", suffix)
      last_dptf_col <- paste0("dptf", suffix)
      last_effdist_col <- paste0("effdist", suffix)

      suffix1 <- as.numeric(str_extract(last_effdist_col, "\\d+$"))

      # Assuming 'suffix' and 'last_effdist_col' have been defined before this code block
      last_effdist_col1 <- NULL

      # If last_effdist_col doesn't exist in the dataframe
      if (!(last_effdist_col %in% names(df6a)) && !(last_Hdist_col %in% names(df6a)) && suffix > 0) {
        suffix2 <- suffix1 - 1
        last_effdist_col1 <- paste0("effdist", suffix2)
        last_Hdist_col1 <- paste0("Hdist", suffix2)
        #print(paste("Trying column:", last_effdist_col1)) # Diagnostic print
      }

      if (!is.null(last_effdist_col1)) {
        previous_effdist <- last_effdist_col1
      }

      # Check the first condition

      effdist_columns <- grep("^effdist", names(df6a), value = TRUE)

      if (length(effdist_columns) == 1) {

        last_effdist_col1 <- effdist_columns[1]
        last_effdist1 <- df6a[[last_effdist_col1]]

        hcbh_cols <- sort(names(df6a)[grep("Hcbh", names(df6a))])
        hdist_cols <- sort(names(df6a)[grep("Hdist", names(df6a))])
        dptf_cols <- sort(names(df6a)[grep("^dptf", names(df6a))])
        hdptf_cols <- sort(names(df6a)[grep("^Hdptf", names(df6a))])

        # Extract the suffix from the last_Hcbh_numeric_col
        suffix <- gsub("Hcbh", "", last(hcbh_cols))

        # Form the corresponding "Hdist" and "Hdepth" column names
        last_Hdist_col <- paste0("Hdist", suffix)
        last_Hdepth_col <- paste0("Hdptf", suffix)
        last_dptf_col <- paste0("dptf", suffix)
        last_effdist_col <- paste0("effdist", suffix)
        last_Hcbh_col <- paste0("Hcbh", suffix)


        if (last_Hdepth_col %in% colnames(df6a) & last_Hdist_col %in% colnames(df6a)) {
          last_df <- data.frame(
            last_Hcbh = df6a[[last_Hcbh_col]],
            last_Hdist = df6a[[last_Hdist_col]],
            last_dptf = df6a[[last_dptf_col]],
            last_Hdptf = df6a[[last_Hdepth_col]],  # Again change to last_Hdptf
            last_effdist = df6a[[last_effdist_col1]])

          names(last_df) <-  c("last_Hcbh", "last_Hdist","last_dptf","last_Hdptf","last_effdist")  # And update here
        }

        if (!last_Hdist_col %in% colnames(df6a)) {
          last_df <- data.frame(
            last_Hcbh = df6a[[last_Hcbh_col]],
            last_Hdist = df6a$Hdist1,
            last_dptf = df6a[[last_dptf_col]],
            last_Hdptf = df6a[[last_Hdepth_col]],  # Again change to last_Hdptf
            last_effdist = df6a[[last_effdist_col1]])

          names(last_df) <-  c("last_Hcbh", "last_Hdist","last_dptf","last_Hdptf","last_effdist")  # And update here
        }

        Hdptf1<- "Hdptf1"
        if (!last_Hdepth_col %in% colnames (df6a) && "Hdptf1" %in% colnames (df6a)) {
          last_df <- data.frame(
            last_Hcbh = df6a[[last_Hcbh_numeric_col]],
            last_Hdist = df6a[[last_Hdist_col]],
            last_Hdptf = df6a[[Hdptf1]],  # Change this to last_Hdptf to match the actual column name
            last_effdist = last_effdist1
          )

          names(last_df) <- c("last_Hcbh", "last_Hdist","last_Hdptf", "last_effdist")  # And update here accordingly
        }
      }


      if(length(effdist) > 0 && !exists("last_df")) {

        # First, select the 'effdist' columns
        effdist_cols <- names(df6a)[str_detect(names(df6a), "^effdist")]
        # Extract effdist values from the first row
        effdist_values <- df6a[1, grep("effdist", names(df6a))]
        # Find the index of the column with last effdist
        last_effdist_col_index <- tail(which(effdist_values == last(effdist_values)), n=1)

        # Define the column names for Hcbh, Hdist, dptf, and Hdptf
        hcbh_cols <- sort(names(df6a)[grep("Hcbh", names(df6a))])
        hdist_cols <- sort(names(df6a)[grep("Hdist", names(df6a))])
        dptf_cols <- sort(names(df6a)[grep("^dptf", names(df6a))])
        hdptf_cols <- sort(names(df6a)[grep("^Hdptf", names(df6a))])

        if (df6a$effdist1 == 1 && all(sapply(df6a[effdist_cols], function(x) x[1] == df6a$effdist1[1]))) {

          last_df <- df6a %>%
            dplyr::rename_with(~ ifelse(. %in% c("Hcbh", "dptf", "Hdist", "Hdptf", "effdist"), paste0("last_", str_remove(., "\\d+$")), .))
        }

        if (length(effdist_cols) == 1 & df6a[[hcbh_cols[1]]][1] == 1.5 && length(hdist_cols) >= 1) {
          # Create a new dataframe with the required columns
          last1 <- df6a[,c( hcbh_cols[2], hdist_cols[2], dptf_cols[2], hdptf_cols[2], effdist_cols [1])]
          # Rename columns with prefix "last_" and remove suffix
          last_df <- last1 %>% dplyr::rename_with(.fn = ~ paste0("last_", str_remove(., "\\d+$")))

        }

        if (length(effdist_cols) > 1 & df6a[[hcbh_cols[1]]][1] == 1.5 && length(hdist_cols) >= 1) {

          suffix <- last_effdist_col_index
          suffix1 <- as.character(as.numeric(suffix) +1)
          #suffix2 <- as.character(as.numeric(suffix) +2)

          # Get the suffixes of the columns
          hcbh_suffixes <- as.numeric(gsub("Hcbh", "", hcbh_cols))
          hdist_suffixes <- as.numeric(gsub("Hdist", "", hdist_cols))
          dptf_suffixes <- as.numeric(gsub("dptf", "", dptf_cols))
          hdptf_suffixes <- as.numeric(gsub("Hdptf", "", hdptf_cols))
          effdist_suffixes <- as.numeric(gsub("effdist", "", effdist_cols))

          # Get the columns whose suffixes are greater or equal to the last effdist suffix
          hcbh_col_p <- hcbh_suffixes >= as.numeric(suffix1)
          hdist_col_p <- hdist_suffixes >= as.numeric(suffix)
          dptf_col_p <- dptf_suffixes >= as.numeric(suffix1)
          hdptf_col_p <- hdptf_suffixes >= as.numeric(suffix1)
          effdist_col_p <- effdist_suffixes >= as.numeric(suffix)

          hcbh_col <- hcbh_cols[hcbh_col_p]
          hdist_col<- last(hdist_cols[hdist_col_p])
          dptf_col<- dptf_cols[dptf_col_p]
          hdptf_col<- hdptf_cols[hdptf_col_p]
          effdist_col<- effdist_cols[effdist_col_p]

          # Check if hcbh_col, dptf_col, hdptf_col exist in df6a
          if (length(hcbh_col) == 0 | length(dptf_col) == 0 | length(hdptf_col) == 0) {
            hcbh_col <- last(hcbh_cols)
            dptf_col <- last(dptf_cols)
            hdptf_col <- last(hdptf_cols)
            hdist_col <- last(hdist_col)
          }

          # Create a new dataframe with the required columns
          last1 <- df6a[, c(hcbh_col, hdist_col, dptf_col, hdptf_col, effdist_col)]
          # Rename columns with prefix "last_" and remove suffix
          last_df <- last1 %>% dplyr::rename_with(.fn = ~ paste0("last_", str_remove(., "\\d+$")))
        }

        if (length(effdist_cols) > 1 & df6a[[hcbh_cols[1]]][1] > 1.5 && length(hdist_cols) >= 1) {

          # Define the column names for Hcbh, Hdist, dptf, and Hdptf
          hcbh_cols <- sort(names(df6a)[grep("Hcbh", names(df6a))])
          hdist_cols <- sort(names(df6a)[grep("Hdist", names(df6a))])
          dptf_cols <- sort(names(df6a)[grep("^dptf", names(df6a))])
          hdptf_cols <- sort(names(df6a)[grep("^Hdptf", names(df6a))])

          effdist_cols <- names(df6a)[str_detect(names(df6a), "^effdist")]
          effdist_values <- df6a[1, grep("effdist", names(df6a))]
          last_effdist_col_index <- tail(which(effdist_values == last(effdist_values)), n=1)

          suffix <- last_effdist_col_index
          #suffix2 <- as.character(as.numeric(suffix) +2)

          # Get the suffixes of the columns
          hcbh_suffixes <- as.numeric(gsub("Hcbh", "", hcbh_cols))
          hdist_suffixes <- as.numeric(gsub("Hdist", "", hdist_cols))
          dptf_suffixes <- as.numeric(gsub("dptf", "", dptf_cols))
          hdptf_suffixes <- as.numeric(gsub("Hdptf", "", hdptf_cols))
          effdist_suffixes <- as.numeric(gsub("effdist", "", effdist_cols))

          # Get the columns whose suffixes are greater or equal to the last effdist suffix
          hcbh_col_p <- hcbh_suffixes >= as.numeric(suffix)
          dptf_col_p <- dptf_suffixes >= as.numeric(suffix)
          hdptf_col_p <- hdptf_suffixes >= as.numeric(suffix)
          hdist_col_p <- hdist_suffixes >= as.numeric(suffix)
          effdist_col_p <- effdist_suffixes >= as.numeric(suffix)

          hcbh_col <- last(hcbh_cols[hcbh_col_p])
          hdist_col<- last(hdist_cols[hdist_col_p])
          dptf_col<- last(dptf_cols[dptf_col_p])
          hdptf_col<- last(hdptf_cols[hdptf_col_p])
          effdist_col<- effdist_cols[effdist_col_p]

          # Create a new dataframe with the required columns
          last1 <- df6a[,c(hcbh_col, hdist_col, dptf_col, hdptf_col,effdist_col)]
          # Rename columns with prefix "last_" and remove suffix
          last_df <- last1 %>% dplyr::rename_with(.fn = ~ paste0("last_", str_remove(., "\\d+$")))


          if(!hcbh_col %in% colnames(df6a) && !dptf_col %in% colnames(df6a) && !hdptf_col %in% colnames(df6a) ) {

            hcbh_col <- last (hcbh_cols)
            dptf_col<- last (dptf_cols)
            hdptf_col<- last (hdptf_cols)

            # Create a new dataframe with the required columns
            last1 <- df6a[,c(hcbh_col, dptf_col, hdptf_col)]
            # Rename columns with prefix "last_" and remove suffix
            last_df <- last1 %>% dplyr::rename_with(.fn = ~ paste0("last_", str_remove(., "\\d+$")))
          }

        } }

      effdist_cols <- names(df6a)[str_detect(names(df6a), "^effdist")]

      if(length(effdist_cols) == 0 | all(df6a[, effdist_cols] %in% c(0, 1)) && !exists("last_df") ) { ### max-last Hdist
        # Define the column names for Hcbh, Hdist, dptf, and Hdptf

        hcbh_cols <- grep("^Hcbh", names(df6a), value=TRUE)
        last_hcbh_col <- hcbh_cols[length(hcbh_cols)]
        hdist_cols <- grep("^Hdist", names(df6a), value=TRUE)
        last_hdist_col <- hdist_cols[length(hdist_cols)]
        dptf_cols <- grep("^dptf", names(df6a), value=TRUE)
        last_dptf_col <- dptf_cols[length(dptf_cols)]
        hdptf_cols <- grep("^Hdptf", names(df6a), value=TRUE)
        last_hdptf_col <- hdptf_cols[length(hdptf_cols)]
        effdist_cols <- grep("^effdist", names(df6a), value=TRUE)
        last_effdist_col <- effdist_cols[length(effdist_cols)]

        last_df <- df6a %>%
          dplyr::select(last_Hcbh = {{last_hcbh_col}},
                        last_Hdptf = {{last_hdptf_col}},
                        last_dptf = {{last_dptf_col}},
                        last_Hdist = {{last_hdist_col}},
                        last_effdist = {{last_effdist_col}}
          )
      }

    }

    df6f<-data.frame(df6a, max_df, last_df)


    ###################  rename columns

    # Exclude columns with prefixes: last_, max_, treeID
    cols_to_exclude <- grep("^(last_|max_|treeID)", names(df6f), value = TRUE)
    excluded_data <- df6f[, cols_to_exclude] # Store the excluded columns
    df6f <- df6f[ , !(names(df6f) %in% cols_to_exclude)]

    # Extract unique prefixes
    prefixes <- unique(gsub("([a-zA-Z]+).*", "\\1", names(df6f)))

    # Rename the columns based on the extracted prefixes
    for (prefix in prefixes) {
      # Identify columns with the current prefix
      cols <- grep(paste0("^", prefix), names(df6f))

      # Generate new column names with consecutive suffixes
      new_names <- paste0(prefix, 1:length(cols))

      # Assign new names to the columns
      names(df6f)[cols] <- new_names
    }

    df6f <- cbind(df6f, excluded_data)

    ###################

    # Find duplicated column names
    dup_columns <- names(df6f)[duplicated(names(df6f)) | duplicated(names(df6f), fromLast = TRUE)]
    # Keep one and remove the duplicates
    for (col in dup_columns) {
      columns_to_check <- which(names(df6f) == col)
      # Assuming columns with same name have the same values, then
      # keep the first column and remove the rest.
      df6f <- df6f[-columns_to_check[-1]]
    }

    # Adjust the regex to only match effdist columns with a numeric suffix
    if (df6f$Hcbh1[1] == 1.5 && sum(grepl("effdist[0-9]+$", colnames(df6f))) >= sum(grepl("Hcbh[0-9]+$", colnames(df6f)))) {
      # Find the last 'effdist' column with numeric suffix
      last_effdist_col <- tail(grep("effdist[0-9]+$", colnames(df6f)), 1)
      # Drop the column
      df6f <- df6f[,-last_effdist_col]
    }

    if(df6f$effdist1 == 1 || df6f$effdist1 == 0 && df6f$Hdptf1 == 0.5) {
      df6f$Hdepth1<-NULL
    }
  }
  #########################################

  max_height<-data.frame(df5b$max1)
  names(max_height)<-"max_height"

  if(!("max_height" %in% colnames(df6f))) {

    df6f<-data.frame(df6f, max_height)
  }

  cols_to_exclude <- grep("treeID|treeID1)", names(df6f), value = TRUE)
  trees<-df6f[ , (names(df6f) %in% cols_to_exclude)]
  df6f <- df6f[ , !(names(df6f) %in% cols_to_exclude)]

  df6f <- data.frame (trees,df6f)

  if("treeID2" %in% colnames (df6f)){
    treeID2 <-df6f$treeID2
    df6f <- df6f %>%
      dplyr::rename(
        treeID= treeID2,
        treeID1 = treeID1)
  }

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

