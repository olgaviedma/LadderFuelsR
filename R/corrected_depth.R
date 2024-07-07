# Declare global variables to avoid R CMD check warnings
utils::globalVariables(c("Hcbh", "Hdepth", "Hdepth1", "Hdist", "Hdist1",
                         "cbh", "cbh1", "depth", "depth1", "dist1",
                         "gap", "index1", "last_cbh", "last_nonzero_Hdist",
                         "max_dist", "temp_Hdepth", "temp_depth",
                         "total_depth", "total_dist"))
#' Effective fuel layers depth
#' @description This function recalculates fuel layers depth after considering distances greater than the actual height bin step.
#' @usage get_real_depths (effective_fbh, step=1, min_height=1.5, verbose=TRUE)
#' @param effective_fbh tree metrics with the recalculated base height of fuel layers after considering distances greater than any number of height bin steps
#' (output of [get_real_fbh()] function).An object of the class text.
#' @param step Numeric value for the actual height bin step (in meters).
#' @param min_height Numeric value for the actual minimum base height (in meters).
#' @param verbose Logical, indicating whether to display informational messages (default is TRUE).
#' @return A data frame giving new fuel layers depth after considering distances greater than the actual height bin step.
#' @author Olga Viedma, Carlos Silva, JM Moreno and A.T. Hudak
#'
#' @details
#' # List of tree metrics:
#' \itemize{
#' \item treeID: tree ID with strings and numeric values
#' \item treeID1: tree ID with only numeric values
#' \item dist: Distance between consecutive fuel layers (m)
#' \item Hdist - Height of the distance between consecutive fuel layers (m)
#' \item Hcbh - Base height of each fuel separated by a distance greater than the certain number of steps
#' \item dptf - Depth of fuel layers (m) after considering distances greater than the actual height bin step
#' \item Hdptf - Height of the depth of fuel layers (m) after considering distances greater than the actual height bin step
#' \item max_height - Maximum height of the tree profile
#' }
#' @examples
#' library(magrittr)
#' library(tidyr)
#' library(dplyr)
#'
#' # Before running this example, make sure to run get_real_fbh().
#' if (interactive()) {
#' effective_fbh <- get_real_fbh()
#' LadderFuelsR::effective_fbh$treeID <- factor(LadderFuelsR::effective_fbh$treeID)
#'
#' trees_name1 <- as.character(effective_fbh$treeID)
#' trees_name2 <- factor(unique(trees_name1))
#'
#' depth_metrics_corr_list <- list()
#'
#' for (i in levels(trees_name2)){
#' # Filter data for each tree
#' tree3 <- effective_fbh |> dplyr::filter(treeID == i)
#' # Get real depths for each tree
#' depth_metrics_corr <- get_real_depths(tree3, step=1, min_height=1.5,verbose=TRUE)
#' depth_metrics_corr_list[[i]] <- depth_metrics_corr
#' }
#'
#' # Combine depth values for all trees
#' effective_depth <- dplyr::bind_rows(depth_metrics_corr_list)
#'
#' # Reorder columns
#' original_column_names <- colnames(effective_depth)
#'
#' # Specify prefixes
#' desired_order <- c("treeID", "Hcbh", "dptf", "dist", "Hdist", "Hdptf", "max_height")
#'
#' # Identify unique prefixes
#' prefixes <- unique(sub("^([a-zA-Z]+).*", "\\1", original_column_names))
#' # Initialize vector to store new order
#' new_order <- c()
#'
#' # Loop over desired order of prefixes
#' for (prefix in desired_order) {
#'   # Find column names matching the current prefix
#'   matching_columns <- grep(paste0("^", prefix), original_column_names, value = TRUE)
#'   # Append to the new order
#'   new_order <- c(new_order, matching_columns)
#' }
#' effective_depth <- effective_depth[, new_order]
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
#' @seealso \code{\link{get_renamed0_df}}
#' @seealso \code{\link{get_real_fbh}}
#' @export
get_real_depths <- function (effective_fbh,  step=1, min_height=1.5, verbose=TRUE) {

  if(min_height==0){
    min_height <-0.5
  }

  #remove the columns from the dataframe df2a which contain only NA values.
  df<- effective_fbh
  df2a <- df[, colSums(!is.na(df)) > 0]

  if (verbose) {
    message("Unique treeIDs:", paste(unique(df2a$treeID), collapse = ", "))
  }

  # Rename cbh columns to Hcbh with their numeric suffixes
  names(df2a) <- gsub("^cbh", "Hcbh", names(df2a))

  ##check if the first "dist" column has a value greater than 1. If this is true, then it will remove depth0 in computation

  # Identify the "dist" and "depth" column names
  dist_cols <- grep("^dist", colnames(df2a), value = TRUE)

  if (length(dist_cols) == 0) {
    df2a$dist1 <- 0
  }

  if (length(dist_cols) == 1  && all(df2a$Hcbh1 <= min_height)) {
    df2a$dist1 <- 0
    df2a$Hdist1 <- min_height
  }


  if (length(dist_cols) == 1 && all(df2a$Hcbh1 > min_height) && all(df2a$Hcbh1 > step)) {
    df2a$Hdist1 <- abs(df2a$Hcbh1 - step)
    df2a$dist1 <-  floor(df2a$Hcbh1 - step)
  }

  if (length(dist_cols) == 1 && all(df2a$Hcbh1 > min_height) && all(df2a$Hcbh1 == step)) {
    df2a$Hdist1 <- df2a$Hcbh1
    df2a$dist1 <-  df2a$Hdist1
  }


  if (all(df2a$Hcbh1 > min_height) && all(df2a$Hcbh1 > step)) {
    df2a$Hdist1 <- abs(df2a$Hcbh1 - step)
    df2a$dist1 <-  floor(df2a$Hcbh1 - step)
  }

  if (all(df2a$Hcbh1 > min_height) && all(df2a$Hcbh1 == step)) {
    df2a$Hdist1 <- df2a$Hcbh1
    df2a$dist1 <-  df2a$Hdist1
  }

  if (all(df2a$Hcbh1 <= min_height)) {
    df2a$dist1 <- 0
    df2a$Hdist1 <- min_height
  }

  if (all(df2a$Hcbh1 <= min_height) && all(df2a$Hdepth1 == 0)) {
    df2a$dist1 <- 0
    df2a$Hdist1 <- min_height
    df2a$Hdepth1 <- min_height
  }


  df4 <- df2a

  # Function to extract prefix and suffix
  extract_prefix_suffix <- function(colname) {
    match <- regexpr("\\d+$", colname)
    if (match == -1) {
      prefix <- colname
      suffix <- NA
    } else {
      prefix <- substr(colname, 1, match - 1)
      suffix <- as.integer(substr(colname, match, nchar(colname)))
    }
    list(prefix = prefix, suffix = suffix)
  }

  # Apply the function to the column names
  col_info <- sapply(names(df2a), extract_prefix_suffix, simplify = FALSE)

  # Convert to a dataframe for better readability
  col_info_df <- do.call(rbind, lapply(col_info, as.data.frame))
  col_info_df <- cbind(col_name = names(df2a), col_info_df)
  rownames(col_info_df) <- NULL


  col_info_df$colname <- names(df2a)

  # Ordenar por prefijo y luego por sufijo numérico
  col_info_df <- col_info_df[order(col_info_df$prefix, col_info_df$suffix), ]

  # Reordenar las columnas del dataframe original según el orden calculado
  df4 <- df4[, col_info_df$colname]

    #########################################################3
  # Identify Hcbh and Hdepth columns
  Hcbh_cols <- grep("^Hcbh", names(df4), value = TRUE)
  Hdepth_cols <- grep("^Hdepth", names(df4), value = TRUE)
  depth_cols <- grep("^depth", names(df4), value = TRUE)

  # Check if more Hcbh columns than Hdepth columns
  if (length(Hcbh_cols) > length(Hdepth_cols)) {
    # Calculate how many new Hdepth columns to create
    num_new_Hdepth <- length(Hcbh_cols) - length(Hdepth_cols)

    # Get the values from the last existing Hdepth column
    last_Hdepth_col <- df4[[tail(Hdepth_cols, 1)]]

    # Create new Hdepth columns with the values from the last existing column
    for (i in seq_len(num_new_Hdepth)) {
      new_col_name <- paste0("Hdepth", length(Hdepth_cols) + i)
      df4[[new_col_name]] <- last_Hdepth_col
    }
  }

  # Check if more Hcbh columns than depth columns
  if (length(Hcbh_cols) > length(depth_cols)) {
    # Calculate how many new depth columns to create
    num_new_depth <- length(Hcbh_cols) - length(depth_cols)

    # Get the values from the last existing depth column
    last_depth_col <- df4[[tail(depth_cols, 1)]]

    # Create new depth columns with the values from the last existing column
    for (i in seq_len(num_new_depth)) {
      new_col_name <- paste0("depth", length(depth_cols) + i)
      df4[[new_col_name]] <- last_depth_col
    }
  }

    #######################################3

  # Identify the "dist" and "depth" column names
  dist_cols <- grep("^dist", colnames(df4), value = TRUE)
  depth_cols <- grep("^depth", colnames(df4), value = TRUE)
  Hdepth_cols <- grep("^Hdepth", colnames(df4), value = TRUE)
  Hcbh_cols <- grep("^Hcbh", colnames(df4), value = TRUE)

  cols_to_transpose <- colnames(df4)
  # Exclude both treeID and treeID1 columns from transposition
  cols_to_transpose <- setdiff(cols_to_transpose, c("treeID", "treeID1","treeID2", "max1", "max_height"))

  # Transpose the dataframe excluding treeID and treeID1 columns
  df_transposed0 <- df4 %>%
    dplyr::select(all_of(cols_to_transpose)) %>%
    pivot_longer(cols = everything(), names_to = c(".value", "index"), names_pattern = "(\\D+)(\\d+)?")
  index<-df_transposed0$index

  df_transposed0 <- df_transposed0 %>%
    dplyr::arrange(as.numeric(index))

  Hdepth<- df_transposed0$Hdepth

  ###############################################

  df_transposed <- df_transposed0

  # Define a function to update depth based on conditions
  update_depth_if_needed <- function(df) {
    for (i in seq_len(nrow(df))) {
      if (is.na(df$depth[i]) && !is.na(df$Hdepth[i])) {
        df$depth[i] <- abs(df$Hdepth[i] - df$Hcbh[i])
      }
    }
    return(df)
  }

  df_transposed <- update_depth_if_needed(df_transposed)

 #####################################
  # Function to move Hdepth values based on Hdist condition

  move_Hdepth <- function(df) {
    # Identify the Hdepth values to move and their original indices
    Hdepth_values <- df$Hdepth[df$Hdepth > 0 & !is.na(df$Hdepth)]
    original_indices <- which(df$Hdepth > 0 & !is.na(df$Hdepth))

    # Identify the indices where Hdist > 0 and are not NA
    Hdist_indices <- which(df$Hdist > 0 & !is.na(df$Hdist))

    # Ensure there are enough indices to move Hdepth values to
    num_to_move <- min(length(Hdepth_values), length(Hdist_indices))

    # Create a new vector to store the moved Hdepth values
    new_Hdepth <- df$Hdepth

    # Move each Hdepth value to the corresponding row where Hdist > 0
    for (i in seq_len(num_to_move)) {
      new_Hdepth[Hdist_indices[i]] <- Hdepth_values[i]
    }

    # Assign the updated Hdepth values back to the dataframe
    df$Hdepth <- new_Hdepth

    # Convert Hdepth values to NA where Hdist is 0 or NA
    df$Hdepth <- if_else(df$Hdist == 0 | is.na(df$Hdist), NA_real_, df$Hdepth)

    return(df)
  }

  # Apply the function
  df_transposed_updated <- move_Hdepth(df_transposed)

  ########################################

  move_depth <- function(df) {
    # Identify the depth values to move and their original indices
    depth_values <- df$depth[df$depth > 0 & !is.na(df$depth)]
    original_indices <- which(df$depth > 0 & !is.na(df$depth))

    # Identify the indices where Hdist > 0 and are not NA
    Hdist_indices <- which(df$Hdist > 0 & !is.na(df$Hdist))

    # Ensure there are enough indices to move depth values to
    num_to_move <- min(length(depth_values), length(Hdist_indices))

    # Create a new vector to store the moved depth values
    new_depth <- df$depth

    # Move each depth value to the corresponding row where Hdist > 0
    for (i in seq_len(num_to_move)) {
      new_depth[Hdist_indices[i]] <- depth_values[i]
    }

    # Assign the updated depth values back to the dataframe
    df$depth <- new_depth

    # Convert depth values to NA where Hdist is 0 or NA
    df$depth <- if_else(df$Hdist == 0 | is.na(df$Hdist), NA_real_, df$depth)

    return(df)
  }

  # Apply the function
  df_transposed_updated1 <- move_depth(df_transposed_updated)

  ########################################

  # Convert NA values in Hdist to 0
  df_transposed_updated1 <- df_transposed_updated1 %>%
    dplyr::mutate(Hdist = replace_na(Hdist, 0))

  # Propagate last non-NA Hcbh value into NA rows
  df_transposed_updated1 <- df_transposed_updated1 %>%
    dplyr::mutate(Hcbh = if_else(is.na(Hcbh), lag(Hcbh, default = Hcbh[1]), Hcbh))

  ########################################

  # Propagate the last depth value > 0 into NA or 0 values
  df_transposed_updated1 <- df_transposed_updated1 %>%
    group_by(Hcbh) %>%
    mutate(temp_depth = if_else(depth > 0, depth, NA_real_)) %>%
    fill(temp_depth, .direction = "down") %>%
    mutate(depth = if_else(is.na(depth) | depth == 0, temp_depth, depth)) %>%
    dplyr::select(-temp_depth) %>%
    ungroup()

  # Propagate the last Hdist value > 0
  df_transposed_updated1 <- df_transposed_updated1 %>%
    dplyr::group_by(Hcbh) %>%
    dplyr::mutate(last_nonzero_Hdist = last(na.omit(Hdist[Hdist > 0]), default = dplyr::first(na.omit(Hdist[Hdist > 0]))),
           Hdist = if_else(is.na(Hdist) | Hdist == 0, last_nonzero_Hdist, Hdist)) %>%
    dplyr::select(-last_nonzero_Hdist) %>%
    dplyr::ungroup()

  # Propagate the last Hdepth value > 0
  df_transposed_updated1 <- df_transposed_updated1 %>%
    dplyr::group_by(Hcbh) %>%
    dplyr::mutate(temp_Hdepth = if_else(Hdepth > 0, Hdepth, lag(Hdepth, default = Hdepth[1]))) %>%
    fill(temp_Hdepth, .direction = "down") %>%
    dplyr::mutate(Hdepth = if_else(is.na(Hdepth) | Hdepth == 0, temp_Hdepth, Hdepth)) %>%
    dplyr::select(-temp_Hdepth) %>%
    dplyr::ungroup()


  # Keep only the last row of duplicated Hcbh values and retain maximum dist
  # Assign the maximum dist value to all rows in groups before filtering

  # Compute the maximum dist for each Hcbh group
  max_dist_per_Hcbh <- df_transposed_updated1 %>%
    dplyr::group_by(Hcbh) %>%
    dplyr::summarize(max_dist = max(dist), .groups = 'drop')

  # Join the max dist back to the original dataframe and update dist
  df_transposed_updated1 <- df_transposed_updated1 %>%
    left_join(max_dist_per_Hcbh, by = "Hcbh") %>%
    dplyr::mutate(dist = if_else(Hcbh %in% Hcbh[duplicated(Hcbh)], max_dist, dist)) %>%
    dplyr::select(-max_dist)


  #########################################

  # Filter to keep rows with max Hdepth for each unique Hcbh, avoiding NA values in Hdepth
  df_transposed_updated1 <- df_transposed_updated1 %>%
    group_by(Hcbh) %>%
    dplyr::filter(
      !is.na(Hdepth) & Hdepth == max(Hdepth, na.rm = TRUE))  %>%
    ungroup()

  #########################################

  # Check if all Hdist values are non-NA and greater than 0
  if (all(!is.na(df_transposed_updated1$Hdist)) && all(df_transposed_updated1$Hdist > 0)) {

    # Summarize dist and depth for duplicated Hcbh values
    df_transposed_summarized <- df_transposed_updated1 %>%
      dplyr::group_by(Hcbh) %>%
      dplyr::summarize(
        total_depth = sum(depth) + sum(dist[-dplyr::n()]), # sum of all depth + sum of dist excluding last
        total_dist = sum(dist),
        Hdepth = max(Hdepth),
        Hdist = max(Hdist),
        .groups = 'drop'
      ) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(index = row_number()) %>%
      dplyr::select(index, total_depth, total_dist, Hcbh, Hdepth, Hdist) %>%
      dplyr::rename(depth = total_depth, dist = total_dist)

    # #print the summarized dataframe
    #print(df_transposed_summarized)
  } else {
    #print("Not all Hdist values are non-NA and greater than 0.")
  }

  if (exists ("df_transposed_summarized")) {

    df_transposed_updated1 <-df_transposed_summarized
  } else {

    df_transposed_updated1 <-df_transposed_updated1
  }

  #########################################

  # Filter to keep rows with max Hdepth for each unique Hcbh, avoiding NA values in Hdepth
  df_transposed_updated1 <- df_transposed_updated1 %>%
    group_by(Hcbh) %>%
    dplyr::filter(
      !is.na(Hdepth) & Hdepth == max(Hdepth, na.rm = TRUE))  %>%
    ungroup()

  # Keep only the first row with duplicated values in all variables
  df_transposed_updated1 <- df_transposed_updated1 %>%
    distinct(depth, Hcbh, Hdepth, Hdist, .keep_all = TRUE)


  update_hdist <- function(df) {
    # Check if dataframe has exactly one row
    if (nrow(df) == 1) {
      # Check if Hdist > Hcbh
      if (df$Hdist > df$Hcbh) {
        # Update Hdist to be equal to Hcbh
        df$Hdist <- df$Hcbh -step
      }
    }
    return(df)
  }
  df_transposed_updated1 <- update_hdist(df_transposed_updated1)


  all_equal<- n_distinct(df_transposed_updated1$Hcbh) == 1

  if (nrow(df_transposed_updated1) ==2  && all_equal) {
  # Loop to Hdist > Hcbh with the previous non-NA value
  for (i in 1:nrow(df_transposed_updated1)) {
    df_transposed_updated1$depth[i] <- round(df_transposed_updated1$dist[i + 1] +  df_transposed_updated1$dist [i]  + df_transposed_updated1$depth [i + 1] + df_transposed_updated1$depth [i], 0)
  }


  df_transposed_updated1 <- df_transposed_updated1 %>%
    group_by(Hdist) %>%
    summarise(
      dist = max(dist),
      Hcbh = max(Hcbh),
      depth = max(depth, na.rm =T),
      Hdepth = max(Hdepth)
    )
  }


  ##########################################33

  cols_to_transpose <- grep("^Hdepth", colnames(df_transposed_updated1), value = TRUE)

  df_transposed_updated1$index <- 1:nrow(df_transposed_updated1)

  df_wide <- df_transposed_updated1 %>%
    pivot_wider(names_from = index, values_from = c(all_of(cols_to_transpose), dist, depth, Hcbh,Hdist))

  df_single1 <- df_wide %>%
    dplyr::summarize(across(everything(), ~na.omit(.x)[1]))

  df_single1  <- df_single1 %>%
    dplyr::select_if(~!any(is.na(.)))

  df5b<-cbind.data.frame(df_single1)

  df5b <- df5b %>%
    dplyr::rename_with(~gsub("^Hdepth", "Hdptf", .), starts_with("Hdepth"))

  df5b <- df5b %>%
    dplyr::rename_with(~gsub("^depth", "dptf", .), starts_with("depth"))

  # Remove underscores from column names
  df5b <- df5b %>%
    rename_with(~ str_remove_all(., "_"), everything())


 #######################################

  update_hdepth <- function(df) {
    # Get column names of Hdepth and Hcbh
    hdepth_cols <- grep("^Hdptf\\d+$", names(df), value = TRUE)
    hcbh_cols <- grep("^Hcbh\\d+$", names(df), value = TRUE)

    # Iterate over each pair of columns
    for (i in seq_along(hdepth_cols)) {
      hdepth_col <- hdepth_cols[i]
      hcbh_col <- hcbh_cols[i]

      # Update Hdepth column where Hdepth < Hcbh
      df[[hdepth_col]] <- pmax(df[[hdepth_col]], df[[hcbh_col]])
    }

    return(df)
  }
  df5b <- update_hdepth(df5b)

     ###################  rename columns
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


  if (df5b$Hcbh1 == min_height) {
    df5b$Hdist1 <- min_height
    df5b$dist1 <- 0
  }


  ######### CORRECT DPTF  ##################################################

  correct_dptf <- function(df) {
    # Get column names related to Hcbh and Hdist
    hcbh_cols <- grep("^Hcbh", names(df), value = TRUE)
    hdptf_cols <- grep("^Hdptf", names(df), value = TRUE)

    # Initialize a vector to store new dptf column names
    new_dptf_cols <- character(0)

    # Iterate over each pair of Hcbh and Hdist columns
    for (i in seq_along(hcbh_cols)) {
      hcbh_col <- hcbh_cols[i]
      hdptf_col <- hdptf_cols[i]

      # Update dptf values based on the formula: dptf = hdptf - hcbh
      index <- sub("^Hcbh_", "", hcbh_col)
      dptf_col <- paste0("dptf", index)
      df[[dptf_col]] <- df[[hdptf_col]] - df[[hcbh_col]]

      # Store the new dptf column name
      new_dptf_cols <- c(new_dptf_cols, dptf_col)
    }

    # Remove old dptf columns that have corresponding Hdptf columns
    for (dptf_col in new_dptf_cols) {
      hdptf_col <- gsub("^dptf", "Hdptf", dptf_col)
      if (hdptf_col %in% colnames(df)) {
        df[[dptf_col]] <- df[[dptf_col]]
        df[[dptf_col]] <- NULL
      }
    }

    # Update dptf columns with the values of dptfHcbh columns
    dptfHcbh_cols <- grep("^dptfHcbh", names(df), value = TRUE)
    for (dptfHcbh_col in dptfHcbh_cols) {
      dptf_col <- gsub("^dptfHcbh", "dptf", dptfHcbh_col)
      df[[dptf_col]] <- df[[dptfHcbh_col]]
    }

    # Remove dptfHcbh columns
    df <- df[, !grepl("^dptfHcbh", names(df))]

    # #print the updated dataframe
    ##print(df)
  }

  # Apply the function to your dataframe (replace df5b with your dataframe name)
  df5b <- correct_dptf(df5b)

  # Replace 0 values with 1 in dptf columns
  df5b <- df5b %>%
    dplyr::mutate(across(starts_with("dptf"), ~ ifelse(. == 0, 1, .)))

 ############################################

   treeID<-unique(factor(df$treeID))
  if(!"treeID" %in% colnames(df5b)) {
    df5b <- cbind(df[c("treeID","treeID1")],df5b )
  }

  max_height<-data.frame(df$max_height)
  names(max_height)<-"max_height"

  if(!"max_height" %in% colnames(df5b)) {
    df5b <- cbind(df5b,df[c("max_height")])
  }

  effective_depth<-df5b

  dptf_cols <- grep("^dptf", colnames(effective_depth), value = TRUE)

  # Apply floor to all identified columns
  effective_depth <- effective_depth %>%
    mutate(across(all_of(dptf_cols), ceiling))

  return (effective_depth)

}
