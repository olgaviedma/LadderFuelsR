#' Leaf Area Density (LAD) percentage comprised in each effective fuel layer
#' @description This function calculates the percentage of Leaf Area Density (LAD) within each fuel layer (first output)
#' and removes those fuel layers with LAD percentage less than a specified threshold (default 10%), recalculating the distances and
#' the depth of the remaining ones (second output).
#' @usage get_layers_lad(LAD_profiles, effective_distances,
#' threshold=10, step = 1, min_height= 1.5, verbose=TRUE)
#' @param LAD_profiles Original tree Leaf Area Density (LAD) profile (output of [lad.profile()] function in the \emph{leafR} package).
#' An object of the class text.
#' @param effective_distances Tree metrics of fuel layers giving the effective distances (> any number of steps) between consecutive fuel layers
#' (output of [get_effective_gap()] function). An object of the class text.
#' @param threshold Numeric value for the minimum required LAD percentage in a fuel layer. The default threshold is 10.
#' @param step Numeric value for the actual height bin step (in meters).
#' @param min_height Numeric value for the actual minimum base height (in meters).
#' @param verbose Logical, indicating whether to display informational messages (default is TRUE).
#' @return A data frame identifying the fuel layers with their corresponding LAD percentage.
#' @author Olga Viedma, Carlos Silva, JM Moreno and A.T. Hudak
#'
#' @details
#'\itemize{
#'   \item treeID: tree ID with strings and numeric values
#'   \item treeID1: tree ID with only numeric values
#'   \item dptf: Depth of fuel layers (m) after considering distances greater than the actual height bin step
#'   \item effdist: Effective distance between consecutive fuel layers (m) after considering distances greater than any number of steps
#'   \item Hcbh: Base height of each fuel separated by a distance greater than the certain number of steps
#'   \item Hdptf: Height of the depth of fuel layers (m) after considering distances greater than the actual step
#'   \item Hdist: Height of the distance (> any number of steps) between consecutive fuel layers (m)
#'   \item Hcbh_Hdptf - Percentage of LAD values comprised in each effective fuel layer
#'   \item max_height - Maximum height of the tree profile
#'   \item nlayers - Number of effective fuel layers
#' }
#' @examples
#' library(magrittr)
#' library(gdata)
#' library(dplyr)
#' library(stringr)
#'
#' # LAD profiles derived from normalized ALS data after applying [lad.profile()] function
#' LAD_profiles <- read.table(system.file("extdata", "LAD_profiles.txt", package = "LadderFuelsR"),
#' header = TRUE)
#' LAD_profiles$treeID <- factor(LAD_profiles$treeID)
#'
#' # Before running this example, make sure to run get_effective_gap().
#' if (interactive()) {
#' effective_distances <- get_effective_gap()
#' LadderFuelsR::effective_distances$treeID <- factor(LadderFuelsR::effective_distances$treeID)
#'
#' trees_name1 <- as.character(effective_distances$treeID)
#' trees_name2 <- factor(unique(trees_name1))
#'
#' LAD_metrics1 <- list()
#' LAD_metrics2 <- list()
#'
#' for (i in levels(trees_name2)) {
#' # Filter data for each tree
#' tree1 <- LAD_profiles |> dplyr::filter(treeID == i)
#' tree2 <- effective_distances |> dplyr::filter(treeID == i)
#'
#' # Get LAD metrics for each tree
#' LAD_metrics <- get_layers_lad(tree1, tree2,
#' threshold=10,
#' step = 1,min_height= 1.5,
#' verbose=TRUE)
#'
#' LAD_metrics1[[i]] <- LAD_metrics$df1
#' LAD_metrics2[[i]] <- LAD_metrics$df2
#' }
#'
#' all_LAD <- dplyr::bind_rows(LAD_metrics1)
#' effective_LAD <- dplyr::bind_rows(LAD_metrics2)
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
#' @seealso \code{\link{get_renamed_df}}
#' @seealso \code{\link{get_effective_gap}}
#' @export
get_layers_lad <- function(LAD_profiles,
                           effective_distances,
                           threshold=10,
                           step = 1,min_height= 1.5,
                           verbose=TRUE) {

  df_orig <- LAD_profiles
  effectiv_gaps<- effective_distances


  if(min_height==0){
    min_height <-0.5

    # Ensure the column starts with a negative value
    if (df_orig$height[1] < min_height) {
      # Calculate the shift value
      shift_value <- abs(df_orig$height[1])

      # Adjust the column to start from 0
      df_orig$height <- df_orig$height + shift_value
    }


    # Ensure the column starts with a negative value
    if (df_orig$height[1] > min_height) {
      # Calculate the shift value
      shift_value1 <- abs(df_orig$height[1])

      # Adjust the column to start from 0
      df_orig$height <- df_orig$height - shift_value1
    }
  }


  effectiv_gaps <- effectiv_gaps[, !apply(effectiv_gaps, 2, function(x) all(is.na(x)))]

  if(all(is.na(effectiv_gaps$Hdist1))){
    effectiv_gaps$Hdist1<-0
  }

  if (!("effdist1" %in% colnames(effectiv_gaps))) {
    effectiv_gaps$effdist1 <- 0
  }

  if (("effdist1" %in% colnames(effectiv_gaps)) && is.na(effectiv_gaps$effdist1)) {
    effectiv_gaps$effdist1 <- 0
  }

  df_effective1 <- effectiv_gaps

  treeID<-"treeID"
  treeID1<-"treeID1"

  ##print(paste(unique(df_effective1$treeID), collapse = ", "))

 if (verbose) { message("Unique treeIDs:", paste(unique(df_effective1$treeID), collapse = ", "))}

  ######################################
  Hcbh_cols <- grep("^Hcbh\\d+$", names(df_effective1), value = TRUE)
  Hdist_cols <- grep("^Hdist\\d+$", names(df_effective1), value = TRUE)
  Hdptf_cols <- grep("^Hdptf\\d+$", names(df_effective1), value = TRUE)
  dptf_cols <- grep("^dptf\\d+$", names(df_effective1), value = TRUE)
  effdist_cols <- grep("^effdist\\d+$", names(df_effective1), value = TRUE)

  # Select columns with prefixes: last_, max_, treeID
  cols_to_append <- grep("^(ltreeID)", names(df_effective1), value = TRUE)
  append_df <- df_effective1[ , cols_to_append]

  # Exclude columns with prefixes: last_, max_, treeID
  cols_to_exclude <- grep("^(treeID)", names(df_effective1), value = TRUE)
  df_effective1 <- df_effective1[ , !(names(df_effective1) %in% cols_to_exclude)]

  # Extract unique prefixes
  prefixes <- unique(gsub("([a-zA-Z]+).*", "\\1", names(df_effective1)))

  # Rename the columns based on the extracted prefixes
  for (prefix in prefixes) {
    # Identify columns with the current prefix
    cols <- grep(paste0("^", prefix), names(df_effective1))

    # Generate new column names with consecutive suffixes
    new_names <- paste0(prefix, 1:length(cols))

    # Assign new names to the columns
    names(df_effective1)[cols] <- new_names
  }

  ###################################3

  # Calculate the total "lad" in the series
  total_lad <- sum(df_orig$lad)

  # Empty list to store percentages for each height range
  percentage1 <- list()
  # Empty list to store column names
  col_names1 <- list()

  # Iterate over each pair of "Hcbh" and "Hdepth" columns in df_effective1
  # Determine the unique suffixes based on 'Hcbh' columns
  suffixes <- gsub("Hcbh", "", grep("^Hcbh\\d+$", names(df_effective1), value = TRUE))


  # Iterate over each unique suffix
  for (i in suffixes) {

    # Get the column names for Hcbh and Hdepth
    hcbh_colname <- paste0("Hcbh", i)
    hdepth_colname <- paste0("Hdptf", i)

    # Get the height range values and extract the numeric values
    start_height <- as.numeric(df_effective1[1, hcbh_colname])
    end_height <- as.numeric(df_effective1[1, hdepth_colname]) +1

    # Subset the target data frame based on the height range
    subset_df <- df_orig[df_orig$height >= start_height & df_orig$height <= end_height, ]

    # Calculate the percentage of total "lad" values within the height range
    lad_values <- subset_df$lad
    range_percentage <- sum(lad_values) / total_lad * 100

    # Store the percentage value
    percentage1[[i]] <- range_percentage

    # Create the column name by concatenating Hcbh and Hdepth column names
    col_name <- paste0(hcbh_colname, "_", hdepth_colname)
    col_names1[[i]] <- col_name
  }


  # Filter out elements equal to 0 from percentage1 and corresponding col_names1
  nonzero_indices <- sapply(percentage1, function(x) x != 0)
  percentage1_filtered <- percentage1[nonzero_indices]
  col_names1_filtered <- col_names1[nonzero_indices]

  # Create the dataframe
  output_df <- data.frame(Percentage = unlist(percentage1_filtered))
  output_df <- data.frame(t(output_df))
  colnames(output_df) <- unlist(col_names1_filtered)

  output_df1<-data.frame(df_effective1,output_df)

   ##############################################################

  merged_df <- output_df1

  #####################################################################
  # 1. Identify Hcbh_Hdptf columns with values < 25
  #####################################################################

  merged_df1 <- merged_df

  # Ensure last_lad_suffix is extracted correctly
  lad_columns2 <- grep("^Hcbh\\d+_Hdptf\\d+$", names(merged_df1))
  lad_columns3 <- names(merged_df1)[lad_columns2]
  lad_suffixes <- as.numeric(str_extract(lad_columns3, "\\d+$"))

  # Get the last suffix
  last_lad_suffix <- max(lad_suffixes)

  # Find the columns with value less than 25
  threshold_value <- threshold

  # Use the threshold dynamically in your code
  cols_to_remove <- sapply(merged_df1[, lad_columns3], function(col) any(col < threshold_value))


  all_effdist_cols <- grep("^effdist\\d+$", names(merged_df1), value = TRUE)
  all_effdist_suffixes <- as.numeric(stringr::str_extract(all_effdist_cols, "\\d+$"))
  suffixes_to_remove <- sort(as.numeric(stringr::str_extract(names(cols_to_remove[cols_to_remove]), "\\d+$")))

  cols_no_remove <- sapply(merged_df1[, lad_columns3], function(col) any(col > threshold_value))

  if (any(cols_no_remove)) {

     # Make a copy of the original dataframe to avoid overwriting it
    merged_df1 <- merged_df

    ##########################################################################
    # IF ANY CHANGE
    ##########################################################################

  if (any(cols_to_remove)) {

        # Function to update effdist columns before removing specified columns
    update_effdist_and_remove <- function(merged_df1, threshold_value) {
      # Identify Hcbh_Hdptf columns
      lad_columns2 <- grep("^Hcbh\\d+_Hdptf\\d+$", names(merged_df1))
      lad_columns3 <- names(merged_df1)[lad_columns2]

      # Find columns to remove based on threshold value
      cols_to_remove <- sapply(merged_df1[, lad_columns3], function(col) any(col < threshold_value))
      suffixes_to_remove <- sort(as.numeric(str_extract(names(cols_to_remove[cols_to_remove]), "\\d+$")))

      # Prepare columns to remove including their corresponding dist and dptf columns
      columns_to_remove <- c()
      next_dist_columns_to_remove <- c()
      for (suffix in suffixes_to_remove) {
        dist_col <- paste0("dist", suffix)
        effdist_col <- paste0("effdist", suffix)
        dptf_col <- paste0("dptf", suffix)
        hdptf_col <- paste0("Hdptf", suffix)
        hcbh_col <- paste0("Hcbh", suffix)
        hdist_col <- paste0("Hdist", suffix)
        hcbh_hdptf_col <- paste0("Hcbh", suffix, "_Hdptf", suffix)
        next_dist_col <- paste0("dist", suffix + 1)
        next_effdist_col <- paste0("effdist", suffix + 1)  # Next effdist column

        columns_to_remove <- c(columns_to_remove, dist_col, hdist_col, dptf_col, hdptf_col, hcbh_col, hcbh_hdptf_col)

        if (next_dist_col %in% names(merged_df1)) {
          next_dist_columns_to_remove <- c(next_dist_columns_to_remove, next_dist_col)
        }

        # Also consider removing the next effdist column if the next dist column exists
        if (next_effdist_col %in% names(merged_df1)) {
          next_dist_columns_to_remove <- c(next_dist_columns_to_remove, next_effdist_col)
        }
      }

      # Ensure only existing columns are included in the removal list
      columns_to_remove <- columns_to_remove[columns_to_remove %in% names(merged_df1)]
      next_dist_columns_to_remove <- next_dist_columns_to_remove[next_dist_columns_to_remove %in% names(merged_df1)]

      # Debug: #print columns to be removed
      #print("Columns to remove:")
      #print(columns_to_remove)
      #print("Next dist columns to remove:")
      #print(next_dist_columns_to_remove)

      # Update effdist columns before removing
      for (suffix in suffixes_to_remove) {
        dist_col <- paste0("dist", suffix)
        dptf_col <- paste0("dptf", suffix)
        hdptf_col <- paste0("Hdptf", suffix)
        effdist_col <- paste0("effdist", suffix)
        next_dist_col <- paste0("dist", suffix + 1)  # Next dist column

        if (dist_col %in% names(merged_df1) && dptf_col %in% names(merged_df1) && effdist_col %in% names(merged_df1)) {
          merged_df1[, effdist_col] <- merged_df1[, dist_col] + merged_df1[, dptf_col]
          if (next_dist_col %in% names(merged_df1)) {
            merged_df1[, effdist_col] <- merged_df1[, effdist_col] + merged_df1[, next_dist_col]
          }
        }
      }


      # Remove the identified columns
      merged_df1 <- merged_df1 %>% dplyr::select(-one_of(columns_to_remove))

      # Debug: #print columns before removing next dist columns
      #print("Columns before removing next dist columns:")
      #print(names(merged_df1))

      # Remove the next dist columns if they exist
      merged_df1 <- merged_df1 %>% dplyr::select(-one_of(next_dist_columns_to_remove))

      # Debug: #print columns after removing next dist columns
      #print("Columns after removing next dist columns:")
      #print(names(merged_df1))

      # Remove dist and Hdist columns if they exist
      dist_cols <- grep("^dist\\d+$", names(merged_df1), value = TRUE)
      hdist_cols <- grep("^Hdist", names(merged_df1), value = TRUE)

      if (all(dist_cols %in% colnames(merged_df1)) && all(hdist_cols %in% colnames(merged_df1))) {
        merged_df1 <- merged_df1 %>%
          dplyr::select(-starts_with("dist"))
      }

      return(merged_df1)
    }


    merged_df1 <- update_effdist_and_remove(merged_df1, threshold)
    #merged_df1<- get_renamed_df (merged_df1)

    col_names<- colnames(merged_df1)
    # Extract numeric suffixes from Hcbh and effdist columns
    hcbh_suffixes <- as.numeric(gsub("Hcbh", "", col_names[grep("^Hcbh", col_names)]))
    hcbh_suffixes <- hcbh_suffixes[!is.na(hcbh_suffixes)]  # Remove NAs

    effdist_suffixes <- as.numeric(gsub("effdist", "", col_names[grep("^effdist", col_names)]))
    effdist_suffixes <- effdist_suffixes[!is.na(effdist_suffixes)]  # Remove NAs

    # Identify the maximum suffix in the Hcbh columns
    max_hcbh_suffix <- max(hcbh_suffixes, na.rm = TRUE)

    # Identify effdist columns with a suffix greater than the maximum Hcbh suffix
    invalid_effdist_cols <- paste0("effdist", effdist_suffixes[effdist_suffixes > max_hcbh_suffix])

    # Remove invalid effdist columns
    merged_df1 <- merged_df1[, !colnames(merged_df1) %in% invalid_effdist_cols]

  }


  ###############################################

  if (identical(merged_df, merged_df1)) {


    dist_cols <- grep("^dist\\d+$", names(merged_df1), value = TRUE)
    hdist_cols <- grep("^Hdist", names(merged_df1), value = TRUE)

    if (all(dist_cols %in% colnames (merged_df1)) && all(hdist_cols %in% colnames (merged_df1))){
      merged_df1 <- merged_df1 %>%
        dplyr::select(-starts_with("dist"))

    }
  }

##########################################################################

  tree_columns <- grep("^tree", colnames(effectiv_gaps), value = TRUE)

  for (tree_col in tree_columns) {
    if (!(tree_col %in% colnames(merged_df1))) {
      merged_df1 <- cbind(effectiv_gaps[[tree_col]], merged_df1)
      colnames(merged_df1)[1] <- tree_col
    }
  }

  for (tree_col in tree_columns) {
    if (!(tree_col %in% colnames(merged_df))) {
      merged_df <- cbind(effectiv_gaps[[tree_col]], merged_df)
      colnames(merged_df)[1] <- tree_col
    }
  }

  max_height<-data.frame(effectiv_gaps$max_height)
  names(max_height)<-"max_height"

  if(!"max_height" %in% colnames(merged_df1)) {
    merged_df1 <- cbind(merged_df1, effectiv_gaps[c("max_height")])
  }
  if(!"max_height" %in% colnames(merged_df)) {
    merged_df <- cbind(merged_df, effectiv_gaps[c("max_height")])
  }

  all_LAD <- merged_df
  effective_LAD <- merged_df1

  effective_LAD<- get_renamed_df (effective_LAD)

  cbh_columns1 <- grep("^Hcbh\\d+$", names(all_LAD))  # Extracts numeric cbh columns
  nlayers1 <- sum(sapply(all_LAD[, cbh_columns1], is.numeric))  # Count numeric values

  cbh_columns2 <- grep("^Hcbh\\d+$", names(effective_LAD))  # Extracts numeric cbh columns
  nlayers2 <- sum(sapply(effective_LAD[, cbh_columns2], is.numeric))  # Count numeric values

  # Add a new column to tree1 with the calculated nlayers1 value
  all_LAD$nlayers <- nlayers1
  effective_LAD$nlayers <- nlayers2

  #######################################
  # Check if Hcbh1 is equal to min_height
  # Ensure effective_LAD1 and min_height are defined
  effective_LAD1 <- effective_LAD

  # Extract column indices for Hcbh, Hdist, Hdptf, and effdist
  hcbh_cols <- grep("^Hcbh\\d+$", colnames(effective_LAD1))
  num_hcbh <- length(hcbh_cols)

  # Check if the first Hcbh value is equal to min_height and there is more than one Hcbh column
  if (effective_LAD1$Hcbh1[1] == min_height && num_hcbh > 1) {
    # Get the column names
    col_names <- names(effective_LAD1)

    # Find the indices of columns matching the pattern
    hdist_indices <- grep("^Hdist\\d+$", col_names)
    hdptf_indices <- grep("^Hdptf\\d+$", col_names)
    effdist_indices <- grep("^effdist\\d+$", col_names)

    # Sort indices to ensure correct order
    hdist_indices <- sort(hdist_indices)
    hdptf_indices <- sort(hdptf_indices)
    effdist_indices <- sort(effdist_indices)


    # Check if we have enough effdist columns to update
    if (length(effdist_indices) >= 2) {
      # Calculate and assign effdist values
      for (i in 2:length(effdist_indices)) {
        current_hdptf_col <- col_names[hdptf_indices[i - 1]]
        next_hdist_col <- col_names[hdist_indices[i]]
        effdist_col <- col_names[effdist_indices[i]]

        # Calculate effdist values for columns starting from effdist2
        effdist_value <- effective_LAD1[[next_hdist_col]] - effective_LAD1[[current_hdptf_col]]

        # Assign the calculated value explicitly
        effective_LAD1[[effdist_col]] <- effdist_value
      }
    }
  }



  # Check if Hcbh1 is equal to min_height
  if (effective_LAD1$Hcbh1[1] > min_height && num_hcbh == 1) {
    # Get the column names
    col_names <- names(effective_LAD1)

     effdist_indices <- grep("^effdist\\d+$", col_names)
      effdist_indices <- sort(effdist_indices)

    # Set the first effdist value to Hcbh1
      effective_LAD1[[col_names[effdist_indices[1]]]] <- floor(effective_LAD1$Hdist1)
  }


  # Make a copy of effective_LAD for manipulation
  effective_LAD2 <- effective_LAD1

  # Find Hcbh and effdist columns
  effdist_cols <- grep("^effdist\\d+$", colnames(effective_LAD2))
  hcbh_cols <- grep("^Hcbh\\d+$", colnames(effective_LAD2))

  # Count the number of Hcbh and effdist columns
  num_effdist <- length(effdist_cols)
  num_hcbh <- length(hcbh_cols)

  # Check if there are more effdist columns than Hcbh columns
  if (num_effdist > num_hcbh) {
    # Sort effdist indices in ascending order to remove from the first extra column
    effdist_indices <- sort(effdist_cols, decreasing = FALSE)

    # Remove extra effdist columns starting from the last one
    cols_to_remove <- effdist_indices[(num_effdist - num_hcbh + 1):num_effdist]
    effective_LAD2 <- effective_LAD2[, -cols_to_remove]
  }


  # Return them in a list
  return(list(df1 = all_LAD, df2 = effective_LAD2))

}

}
