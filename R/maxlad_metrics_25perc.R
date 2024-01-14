#' Leaf Area Density (LAD) percentage comprised in each effective fuel layer
#' @description This function calculates the percentage of Leaf Area Density (LAD) within each fuel layer (first output)
#' and removes those fuel layers with LAD percentage less than a specified threshold, recalculating the distances and
#' the depth of the remaining ones (second output).
#' @usage get_layers_lad(LAD_profiles, effective_distances, threshold=25, verbose=TRUE)
#' @param LAD_profiles Original tree Leaf Area Density (LAD) profile (output of [lad.profile()] function in the \emph{leafR} package).
#' An object of the class text.
#' @param effective_distances Tree metrics of fuel layers separated by distances greater than 1 m (output of [get_effective_gap()] function).
#' An object of the class text.
#' @param threshold Numeric value for the minimum required LAD percentage in a fuel layer. The default threshold is 25.
#' @param verbose Logical, indicating whether to display informational messages (default is TRUE).
#' @return A data frame identifying the fuel layers with their corresponding LAD percentage.
#' @author Olga Viedma, Carlos Silva and JM Moreno
#' @details
#'\itemize{
#' \item treeID: tree ID with strings and numeric values
#' \item treeID1: tree ID with only numeric values
#' \item Hcbh - Height of the base of each effective fuel layer (m)
#' \item Hdist - Height of the distance between consecutive fuel layers (m)
#' \item effdist - Distance between consecutive fuel layers (m)
#' \item dptf - Depth of the effective fuel layers (m) at distances greater than 1 m
#' \item Hdptf - Height of the depth of fuel layers (m) at distances greater than 1 m
#' \item Hcbh_Hdptf - Percentage of LAD values comprised in each effective fuel layer
#' \item max_height - Maximum height of the tree profile
#' \item nlayers - Number of effective fuel layers
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
#' LAD_metrics <- get_layers_lad(tree1, tree2, threshold=25, verbose=TRUE)
#' LAD_metrics1[[i]] <- LAD_metrics$df1
#' LAD_metrics2[[i]] <- LAD_metrics$df2
#' }
#'
#' all_LAD <- dplyr::bind_rows(LAD_metrics1)
#' effective_LAD <- dplyr::bind_rows(LAD_metrics2)
#' }
#' @importFrom dplyr select_if group_by summarise summarize mutate arrange rename rename_with filter slice slice_tail ungroup distinct
#' across matches row_number all_of vars first
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
#' @seealso \code{\link{get_renamed_df}}
#' @seealso \code{\link{get_effective_gap}}
#' @seealso \code{\link{remove_no_flayer_noconsec}}
#' @export
get_layers_lad <- function(LAD_profiles, effective_distances,threshold=25, verbose=TRUE) {

  df_orig <- LAD_profiles
  effectiv_gaps<- effective_distances

  if(is.na(effectiv_gaps$Hdist1)){
    effectiv_gaps$Hdist1<-0.5
  }

  if (!("effdist1" %in% colnames(effectiv_gaps))) {
    effectiv_gaps$effdist1 <- 0
  }

  if (("effdist1" %in% colnames(effectiv_gaps)) && is.na(effectiv_gaps$effdist1)) {
    effectiv_gaps$effdist1 <- 0
  }

  df_effective1 <- effectiv_gaps[, !apply(effectiv_gaps, 2, function(x) all(is.na(x)))]
  treeID<-"treeID"
  treeID1<-"treeID1"

  #print(paste(unique(df_effective1$treeID), collapse = ", "))

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

  ##########################################################################
  # IF ANY CHANGE
  ##########################################################################

  if (any(cols_to_remove)) {

    # Extract the numeric suffixes and sort them
    suffixes_to_remove <- sort(as.numeric(stringr::str_extract(names(cols_to_remove[cols_to_remove]), "\\d+$")))

    #########################################################
    ### CONSECUTIVE SUFFIXES TO REMOVE
    #########################################################
    # Make a copy of the original dataframe to avoid overwriting it
    merged_df1 <- merged_df

    lad_columns2 <- grep("^Hcbh\\d+_Hdptf\\d+$", names(merged_df1),value=TRUE)
    lad_suffixes <- as.numeric(str_extract(lad_columns2, "\\d+$"))

    lad_no_remove <- sapply(merged_df1[, lad_columns2], function(col) any(col > threshold_value))
    suffixes_no_remove <- sort(as.numeric(stringr::str_extract(names(lad_no_remove[lad_no_remove]), "\\d+$")))

    colnames_to_extract <- names(lad_no_remove)[lad_no_remove]
    # Extract the values from the desired columns
    lad_values_noremove <- merged_df1[, colnames_to_extract, drop=FALSE]

    # Collapsing multiple patterns into a single pattern
    pattern <- paste(suffixes_no_remove, "$", sep = "", collapse = "|")

    # Use sapply to get matches for each pattern
    matching_cols_list <- sapply(pattern, function(pat) grep(pat, names(merged_df1), value=TRUE))

    # Unlist and make unique the result
    matching_cols <- unique(unlist(matching_cols_list))


    # Find the first Hcbh column from matching_cols
    first_cbh_gt5_col <- first(grep("^Hcbh\\d+$", matching_cols, value = TRUE))

    # extract the values of that column from merged_df1
    if (!is.null(first_cbh_gt5_col) && first_cbh_gt5_col %in% names(merged_df1)) {
      first_cbh_gt5_value <- merged_df1[, first_cbh_gt5_col]
    } else {
      first_cbh_gt5_value <- NA
    }


    last_cbh_gt5_col <- last(grep("^Hcbh\\d+$", matching_cols, value = TRUE))

    if (!is.null(last_cbh_gt5_col) && last_cbh_gt5_col %in% colnames(merged_df1)) {
      last_cbh_gt5_value <- merged_df1[, last_cbh_gt5_col]
    } else {
      last_cbh_gt5_value <- NA  # or some other default value
    }

    ###################################33

    threshold_value <- threshold

    # Find the columns with value less than 5
    cols_to_remove <- sapply(merged_df1[, lad_columns2], function(col) any(col < threshold_value))
    suffixes_to_remove <- sort(as.numeric(stringr::str_extract(names(cols_to_remove[cols_to_remove]), "\\d+$")))

    colnames_to_extract <- names(cols_to_remove)[cols_to_remove]
    # Extract the values from the desired columns
    lad_values_remove <- merged_df1[, colnames_to_extract, drop=FALSE]

    pattern1 <- paste0(suffixes_to_remove, "$")
    # Get column names matching the pattern
    matching_cols1 <- unique(unlist(sapply(pattern1, function(pat) {
      grep(pat, names(merged_df1), value=TRUE)
    })))

    # extract the matching columns
    data_remove <- merged_df1[, matching_cols1]
    # Find indices of columns with names starting with "Hdist"
    hdist_indices <- grep("^Hdist", colnames(data_remove))

    if (length(hdist_indices)> 0) {
      # Get the maximum index (i.e., the last "Hdist" column)
      last_hdist_index <- max(hdist_indices)
      # Remove the column corresponding to that index
      data_remove <- data_remove[, -last_hdist_index, drop = FALSE]
    }

    #####################################

    # Identify blocks of consecutive columns
    consec_blocks <- split(suffixes_to_remove, cumsum(c(1, diff(suffixes_to_remove) != 1)))
    has_single_element_block <- any(sapply(consec_blocks, length) == 1)
    has_some_element_block <- any(sapply(consec_blocks, length) > 1)

    # Identify number of element in blocks
    some_element_blocks <- consec_blocks[sapply(consec_blocks, length) > 1]
    single_element_blocks <- consec_blocks[sapply(consec_blocks, length) == 1]

    # Suffixes from some_element_blocks
    suffixes <-lapply(some_element_blocks, unlist)
    suffixes <- unlist(suffixes)

    pattern1 <- paste0(suffixes, "$")
    # Get column names matching the pattern
    matching_cols1 <- unique(unlist(sapply(pattern1, function(pat) {
      grep(pat, names(merged_df1), value=TRUE)
    })))

    # Subset the dataframe using matching columns
    block_values <- merged_df1[, matching_cols1]
    # Extracting column names to be removed
    cols_to_remove <- grep("^(effdist|Hdist)", names(block_values), value = TRUE)
    # Removing columns from the dataframe
    block_values <- block_values[, !names(block_values) %in% cols_to_remove]

    first_hdist_remove_col <- paste0("Hdist", first(suffixes))

    if (first_hdist_remove_col %in% colnames(merged_df1)) {
      first_hdist_remove_value <- merged_df1[, first_hdist_remove_col]
    } else {
      first_hdist_remove_value <- NA  # or some other default value
    }


    previous_hdptf_noremove_col <- paste0("Hdptf", first(suffixes)-1)

    if (previous_hdptf_noremove_col %in% colnames(merged_df1)) {
      previous_hdptf_noremove_value <- merged_df1[, previous_hdptf_noremove_col]
    } else {
      previous_hdptf_noremove_value <- NA  # or some other default value
    }

    next_hdptf_noremove_col <- paste0("Hdptf", last(suffixes) + 1)

    if (next_hdptf_noremove_col %in% colnames(merged_df1)) {
      next_hdptf_noremove_value <- merged_df1[, next_hdptf_noremove_col]
    } else {
      next_hdptf_noremove_value <- NA  # or some other default value
    }

    next_cbh_gt5_col <- paste0("Hcbh", last(suffixes) + 1)

    if (!is.null(next_cbh_gt5_col) && next_cbh_gt5_col %in% colnames(merged_df1)) {
      next_cbh_gt5_value <- merged_df1[, next_cbh_gt5_col]
    } else {
      next_cbh_gt5_value <- NA  # or some other default value
    }

    next_cbh_noremove_col <- paste0("Hcbh", last(suffixes) + 1)

    if (next_cbh_noremove_col %in% colnames(merged_df1)) {
      next_cbh_noremove_value <- merged_df1[, next_cbh_noremove_col]
    } else {
      next_cbh_noremove_value <- NA  # or some other default value
    }

    ##########################################################

    if (has_some_element_block == TRUE) {

      block <- some_element_blocks[[1]]

      first_consec_suffix <- min(block)
      last_consec_suffix <- max(block)
      next_column_suffix <- last_consec_suffix + 1
      next2_column_suffix <- last_consec_suffix + 2
      previous_column_suffix <- first_consec_suffix - 1
      first_col_name <- paste0("effdist", first_consec_suffix)
      last_col_name <- paste0("effdist", last_consec_suffix)

      if (last_col_name %in% colnames(merged_df1)) {
        last_col_value <- merged_df1[, last_col_name]
      } else {
        last_col_value <- NA  # or some other default value
      }

      valid_cols <- sapply(block, function(suf) paste0("effdist", suf) %in% colnames(merged_df1))

      sufix_block1 <- paste0("effdist", first(block))
      sufix_block2 <- paste0("effdist", last(block))


      effdist_values <- sapply(block[valid_cols], function(suf) merged_df1[[paste0("effdist", suf)]])
      dptf_values <- sapply(block[valid_cols], function(suf) {
        col_name <- paste0("dptf", suf)
        if (col_name %in% names(merged_df1)) {
          return(merged_df1[[col_name]])
        } else {
          return(NA) # or some other default value
        }
      })

      Hdist_remove_cols<- paste0("Hdist", block)
      first_Hdist_remove_cols <-first(Hdist_remove_cols)
      first_Hdist_remove_values <-merged_df1[[first_Hdist_remove_cols]]

      effdist_sum <- sum(effdist_values)
      dptf_sum <- sum(dptf_values, na.rm=TRUE)

      next_effdist_col <- paste0("effdist", next_column_suffix)
      next2_effdist_col <- paste0("effdist", next_column_suffix + 1)
      previous_effdist_col <- paste0("effdist", previous_column_suffix)


      if (next_effdist_col %in% colnames(merged_df1)) {
        next_effdist_value <- merged_df1[, next_effdist_col]
      } else {
        next_effdist_value <- NA  # or some other default value
      }

      if (previous_effdist_col %in% colnames(merged_df1)) {
        previous_effdist_value <- merged_df1[, previous_effdist_col]
      } else {
        previous_effdist_value <- NA  # or some other default value
      }

      Hdist_block1 <-paste0("Hdist", first(block))
      Hdist_block2 <-paste0("Hdist", last(block))

      if (Hdist_block1 %in% colnames(merged_df1)) {
        Hdist_block1_value <- merged_df1[, Hdist_block1]
      } else {
        Hdist_block1_value <- NA  # or some other default value
      }

      if (Hdist_block2 %in% colnames(merged_df1)) {
        Hdist_block2_value <- merged_df1[, Hdist_block2]
      } else {
        Hdist_block2_value <- NA  # or some other default value
      }


      Hdist_next_colname <- paste0("Hdist", next_column_suffix)
      dist_next_colname <- paste0("dist", next_column_suffix)
      Hdptf_next_colname <- paste0("Hdptf", next_column_suffix)

      if (Hdist_next_colname %in% colnames(merged_df1)) {
        Hdist_next_value <- merged_df1[, Hdist_next_colname]
      } else {
        Hdist_next_value <- NA  # or some other default value
      }

      Hdist_next2_colname <- paste0("Hdist", next2_column_suffix)

      if (Hdist_next2_colname %in% colnames(merged_df1)) {
        Hdist_next2_value <- merged_df1[, Hdist_next2_colname]
      } else {
        Hdist_next2_value <- NA  # or some other default value
      }

      if (Hdptf_next_colname %in% colnames(merged_df1)) {
        Hdptf_next_value <- merged_df1[, Hdptf_next_colname]
      } else {
        Hdptf_next_value <- NA  # or some other default value
      }

      if (dist_next_colname %in% colnames(merged_df1)) {
        dist_next_value <- merged_df1[, dist_next_colname]
      } else {
        dist_next_value <- NA  # or some other default value
      }

      Hdist_previous_colname <- paste0("Hdist", previous_column_suffix)
      dist_previous_colname <- paste0("dist", previous_column_suffix)
      Hdptf_previous_colname <- paste0("Hdptf", previous_column_suffix)

      if (any(Hdist_previous_colname %in% colnames(merged_df1))) {
        Hdist_previous_value <- merged_df1[, Hdist_previous_colname]
      } else {
        Hdist_previous_value <- NA  # or some other default value
      }


      if (any(dist_previous_colname %in% colnames(merged_df1))) {
        dist_previous_value <- merged_df1[, dist_previous_colname]
      } else {
        dist_previous_value <- NA  # or some other default value
      }

      if (Hdptf_previous_colname %in% colnames(merged_df1)) {
        Hdptf_previous_value <- merged_df1[, Hdptf_previous_colname]
      } else {
        Hdptf_previous_value <- NA  # or some other default value
      }

      Hdptf_block1 <-paste0("Hdptf", first(block))
      Hdptf_block2 <-paste0("Hdptf", last(block))

      if (Hdptf_block1 %in% colnames(merged_df1)) {
        Hdptf_block1_value <- merged_df1[, Hdptf_block1]
      } else {
        Hdptf_block1_value <- NA  # or some other default value
      }

      if (Hdptf_block2 %in% colnames(merged_df1)) {
        Hdptf_block2_value <- merged_df1[, Hdptf_block2]
      } else {
        Hdptf_block2_value <- NA  # or some other default value
      }

      previous_dptf_col <- paste0("dptf", previous_column_suffix)
      if (all(previous_dptf_col %in% colnames(merged_df1))) {
        previous_dptf_value <- merged_df1[, first(previous_dptf_col)]
      } else {
        previous_dptf_value <- NA  # or some other default value
      }


      Hcbh_cols <- grep("^Hcbh\\d+$", names(merged_df1), value = TRUE)

      if (all(Hcbh_cols %in% colnames(merged_df1))) {
        Hcbh1 <- merged_df1[, first(Hcbh_cols)]
      } else {
        Hcbh1 <- NA  # or some other default value
      }

      first_Hcbh_col <- grep("^Hcbh\\d+$", names(merged_df1), value = TRUE)[1]

      first_Hcbh_value <- merged_df1[, first_Hcbh_col]


      all_effdist_cols <- grep("^effdist\\d+$", names(merged_df1), value = TRUE)
      all_effdist_values <-merged_df1[,all_effdist_cols]

      #########################

      cols_to_check <- c(Hdist_next_colname, first_cbh_gt5_col,first(Hcbh_cols),Hdist_block2,all_effdist_cols)
      if (all(cols_to_check %in% colnames(merged_df1))) {

        if (first_consec_suffix == first(lad_suffixes) &&
            !is.na(Hdist_next_value) &&
            !is.na(first_cbh_gt5_value) &&
            first_cbh_gt5_value  > 1.5  &&
            Hcbh1 == 1.5  &&
            length(suffixes_no_remove)== 1 &&
            (Hdist_block2_value  < first_cbh_gt5_value) && all (sapply(all_effdist_values, function(x) x != "1"))) {

          merged_df1[first_col_name] <- dptf_sum + effdist_sum
          merged_df1 <- merged_df1[, !names(merged_df1) %in% data_remove]

        }}


      cols_to_check <- c(Hdist_next_colname, first_cbh_gt5_col,first(Hcbh_cols))
      if (all(cols_to_check %in% colnames(merged_df1))) {

        if (first_consec_suffix == first(lad_suffixes) &&
            !is.na(Hdist_next_value) &&
            !is.na(first_cbh_gt5_value) &&
            first_cbh_gt5_value  > 1.5  &&
            Hcbh1 == 1.5  &&
            (Hdist_next_value > first_cbh_gt5_value) ) { #&& length(suffixes_no_remove==1)

          merged_df1[first_col_name] <- dptf_sum + effdist_sum

        }}

      cols_to_check <- c(first_cbh_gt5_col,first(Hcbh_cols))
      if (all(cols_to_check %in% colnames(merged_df1))) {

        if (first_consec_suffix == first(lad_suffixes) &&
            is.na(Hdist_next_value) &&
            !is.na(first_cbh_gt5_value) &&
            first_cbh_gt5_value  > 1.5  &&
            Hcbh1 == 1.5) {

          merged_df1[first_col_name] <- dptf_sum + effdist_sum

        }}


      cols_to_check <- c(Hdist_next_colname,  first_cbh_gt5_col, next_hdptf_noremove_col,Hdist_next2_colname,first(Hcbh_cols))
      if (all(cols_to_check %in% colnames(merged_df1))) {

        if (first_consec_suffix == first(lad_suffixes) &&
            !is.na(Hdist_next2_value) &&
            !is.na(next_hdptf_noremove_value) &&
            Hcbh1 == 1.5  &&
            (Hdist_next_value < first_cbh_gt5_value) &&
            length(suffixes_no_remove) > 1  &&
            (Hdist_next2_value > next_hdptf_noremove_value)) {  # corresponds to next effdist col

          merged_df1[first_col_name] <- dptf_sum + effdist_sum
        }}


      cols_to_check <- c(Hdist_next_colname, first_cbh_gt5_col,first(Hcbh_cols),Hdist_block2,all_effdist_cols,next_effdist_col,Hdist_next_colname)
      if (all(cols_to_check %in% colnames(merged_df1))) {

        if (first_consec_suffix == first(lad_suffixes) &&
            !is.na(Hdist_next_value) &&
            !is.na(first_cbh_gt5_value) &&
            first_cbh_gt5_value  > 1.5  &&
            Hcbh1 == 1.5  &&
            length(suffixes_no_remove)== 1 &&
            (Hdist_next_value < first_cbh_gt5_value) &&
            (Hdist_block2_value  < first_cbh_gt5_value) && any (sapply(all_effdist_values, function(x) x == "1"))) {

          merged_df1[first_col_name] <- dptf_sum + effdist_sum + merged_df1[[next_effdist_col]]
          merged_df1[[next_effdist_col]] <- NULL
        }}


      cols_to_check <- c(next_effdist_col,Hdist_next_colname, first_cbh_gt5_col,first(Hcbh_cols))
      if (all(cols_to_check %in% colnames(merged_df1))) {

        if (first_consec_suffix == first(lad_suffixes) &&
            !is.na(Hdist_next_value) &&
            !is.na(first_cbh_gt5_value) &&
            first_cbh_gt5_value  > 1.5  &&
            Hcbh1 > 1.5  &&
            (Hdist_next_value < first_cbh_gt5_value) &&
            all(sapply(all_effdist_values, function(x) x != "1")) &&
            length(suffixes_no_remove) == 1) {

          merged_df1[first_col_name] <- dptf_sum + effdist_sum + merged_df1[[next_effdist_col]]
          merged_df1[[next_effdist_col]] <- NULL

        }}


      cols_to_check <- c(next_effdist_col, Hdist_next_colname,  first_cbh_gt5_col,first(Hcbh_cols))
      if (all(cols_to_check %in% colnames(merged_df1))) {

        if (first_consec_suffix == first(lad_suffixes) &&
            (!next2_effdist_col %in%colnames(merged_df1)) &&
            is.na(Hdist_next2_value) &&
            !is.na(next_hdptf_noremove_value) &&
            first_cbh_gt5_value  > 1.5  &&
            Hcbh1 > 1.5  &&
            (Hdist_next_value < first_cbh_gt5_value) &&
            any (sapply(all_effdist_values, function(x) x == "1"))) {

          merged_df1[first_col_name] <- dptf_sum + effdist_sum + merged_df1[[next_effdist_col]]
          merged_df1[[next_effdist_col]] <- NULL

        }}


      cols_to_check <- c(next_effdist_col, Hdist_next_colname,  first_cbh_gt5_col, next_hdptf_noremove_col,Hdist_next2_colname,first(Hcbh_cols))
      if (all(cols_to_check %in% colnames(merged_df1))) {

        if (first_consec_suffix == first(lad_suffixes) &&
            !is.na(Hdist_next2_value) &&
            !is.na(next_hdptf_noremove_value) &&
            Hcbh1 > 1.5  &&
            (Hdist_next_value < first_cbh_gt5_value) &&
            length(suffixes_no_remove) > 1  &&
            (Hdist_next2_value > next_hdptf_noremove_value)) {  # corresponds to next effdist col

          merged_df1[first_col_name] <- dptf_sum + effdist_sum + merged_df1[[next_effdist_col]]
          merged_df1[[next_effdist_col]] <- NULL

        }}


      cols_to_check <- c( Hdist_next_colname,  first_cbh_gt5_col, next_hdptf_noremove_col,Hdist_next2_colname)
      if (all(cols_to_check %in% colnames(merged_df1))) {

        if (first_consec_suffix == first(lad_suffixes) &&
            !is.na(Hdist_next2_value) &&
            !is.na(next_hdptf_noremove_value) &&
            first_cbh_gt5_value  > 1.5  &&
            (Hdist_next_value > first_cbh_gt5_value) &&
            length(suffixes_no_remove) > 1  &&
            (Hdist_next2_value < next_hdptf_noremove_value)) {

          merged_df1[first_col_name] <- dptf_sum + effdist_sum + merged_df1[[next_effdist_col]]
          merged_df1[[next_effdist_col]] <- NULL

        }}


      cols_to_check <- c(next_effdist_col, Hdist_next_colname,  first_cbh_gt5_col, next_hdptf_noremove_col,first(Hcbh_cols))
      if (all(cols_to_check %in% colnames(merged_df1))) {

        if (first_consec_suffix == first(lad_suffixes) &&
            !is.na(next_hdptf_noremove_value) &&
            Hcbh1 > 1.5  &&
            (Hdist_next_value < next_cbh_noremove_value) &&
            length(suffixes_no_remove) > 1) {  # corresponds to next effdist col

          merged_df1[first_col_name] <- dptf_sum + effdist_sum + merged_df1[[next_effdist_col]]
          merged_df1[[next_effdist_col]] <- NULL
        }}


      cols_to_check <- c(next_effdist_col, next2_effdist_col, Hdist_next_colname,  first_cbh_gt5_col, next_hdptf_noremove_col,Hdist_next2_colname,first(Hcbh_cols))
      if (all(cols_to_check %in% colnames(merged_df1))) {

        if (first_consec_suffix == first(lad_suffixes) &&
            !is.na(Hdist_next2_value) &&
            !is.na(next_hdptf_noremove_value) &&
            first_cbh_gt5_value  > 1.5  &&
            Hcbh1 > 1.5  &&
            (Hdist_next_value < first_cbh_gt5_value) &&
            length(suffixes_no_remove) > 1  &&
            (Hdist_next2_value < next_hdptf_noremove_value)) {

          merged_df1[first_col_name] <- dptf_sum + effdist_sum + merged_df1[[next_effdist_col]]  + merged_df1[[next2_effdist_col]]
          merged_df1[[next_effdist_col]] <- NULL
          merged_df1[[next2_effdist_col]] <- NULL

        }}


      cols_to_check <- c(next_effdist_col, next2_effdist_col, Hdist_next_colname,  first_cbh_gt5_col,first(Hcbh_cols))
      if (all(cols_to_check %in% colnames(merged_df1))) {

        if (first_consec_suffix == first(lad_suffixes) &&
            is.na(Hdist_next2_value) &&
            !is.na(next_hdptf_noremove_value) &&
            first_cbh_gt5_value  > 1.5  &&
            Hcbh1 > 1.5  &&
            (Hdist_next_value < first_cbh_gt5_value) &&
            any (sapply(all_effdist_values, function(x) x == "1"))) {

          merged_df1[first_col_name] <- dptf_sum + effdist_sum + merged_df1[[next_effdist_col]]  + merged_df1[[next2_effdist_col]]
          merged_df1[[next_effdist_col]] <- NULL
          merged_df1[[next2_effdist_col]] <- NULL

        }}


      cols_to_check <- c(first_hdist_remove_col,  last_cbh_gt5_col, previous_hdptf_noremove_col)
      if (all(cols_to_check %in% colnames(merged_df1))) {

        if (first_consec_suffix == first(lad_suffixes) &&
            !is.na(first_hdist_remove_value) && !is.na(previous_hdptf_noremove_value) &&
            (first_hdist_remove_value > previous_hdptf_noremove_value) &&
            length(suffixes_no_remove) ==1 &&
            last_cbh_gt5_value > 1.5) {

          merged_df1 <- merged_df1[, !(names(merged_df1) %in% matching_cols1)]
          merged_df1[[previous_effdist_col]] <- NULL
        }
        }

      ###############################################3
      ###############################################3

      last_col_name <- paste0("effdist", last(block))
      hdptf_cols <- grep("^Hdptf\\d+$", names(merged_df1), value = TRUE)

      # if (last_col_name %in% colnames(merged_df1)) {

      cols_to_check <- c(first_Hdist_remove_cols , previous_hdptf_noremove_col,Hdist_previous_colname)
      if (all(cols_to_check %in% colnames(merged_df1))) {

        if (isTRUE(first_consec_suffix != first(lad_suffixes) && (first_Hdist_remove_values  > previous_hdptf_noremove_value)
                   && length(suffixes_no_remove) ==1) && (!last_col_name %in% colnames(merged_df1))) {

          merged_df1 <- merged_df1[, !(names(merged_df1) %in% matching_cols1)]
          merged_df1[[previous_effdist_col]] <- NULL

        }}


      cols_to_check <- c(first_Hdist_remove_cols , previous_hdptf_noremove_col,Hdist_previous_colname)
      if (all(cols_to_check %in% colnames(merged_df1))) {

        if (isTRUE(first_consec_suffix != first(lad_suffixes) && (first_Hdist_remove_values  > previous_hdptf_noremove_value)
                   && length(suffixes_no_remove) > 1) && (!last_col_name %in% colnames(merged_df1))) {

          merged_df1 <- merged_df1[, !(names(merged_df1) %in% matching_cols1)]
          merged_df1[[previous_effdist_col]] <- NULL

        }}


      cols_to_check <- c(first_Hdist_remove_cols , previous_hdptf_noremove_col)
      if (all(cols_to_check %in% colnames(merged_df1))) {

        if (isTRUE(first_consec_suffix != first(lad_suffixes) && (first_Hdist_remove_values  > previous_hdptf_noremove_value)
                   && length(suffixes_no_remove) ==1) && (last_col_name %in% colnames(merged_df1))) {

          merged_df1 <- merged_df1[, !(names(merged_df1) %in% matching_cols1)]

        }}


      cols_to_check <- c(Hdist_previous_colname, previous_hdptf_noremove_col, Hdist_block1,last_cbh_gt5_col,next_hdptf_noremove_col)
      if (all(cols_to_check %in% colnames(merged_df1))) {

        if (isTRUE(first_consec_suffix != first(lad_suffixes) &&
                   (Hdist_previous_value < previous_hdptf_noremove_value) &&
                   (first_Hdist_remove_values  <  next_hdptf_noremove_value) &&
                   (Hdist_block1_value  > last_cbh_gt5_value))) {

          merged_df1[sufix_block1] <- dptf_sum + effdist_sum
          #merged_df1 <- merged_df1[, !(names(merged_df1) %in% matching_cols1)]

        }}


      cols_to_check <- c(previous_effdist_col, Hdist_previous_colname,previous_hdptf_noremove_col,Hdptf_next_colname,first_Hcbh_col,previous_dptf_col)
      if (all(cols_to_check %in% colnames(merged_df1))) {

        if (isTRUE(first_consec_suffix != first(lad_suffixes) &&
                   (Hdist_previous_value > previous_hdptf_noremove_value) &&
                   length (suffixes_no_remove) > 1 &&
                   first_Hcbh_value == 1.5) &&
            merged_df1[[Hdist_block2]] < Hdptf_next_value ) {


          merged_df1[sufix_block1] <- dptf_sum + effdist_sum + merged_df1[[previous_effdist_col]]
          merged_df1[[previous_effdist_col]] <- NULL
        }}


      cols_to_check <- c(previous_effdist_col, Hdist_previous_colname , Hdist_next_colname, previous_hdptf_noremove_col, next_hdptf_noremove_col)
      if (all(cols_to_check %in% colnames(merged_df1))) {

        if (isTRUE(first_consec_suffix != first(lad_suffixes) && (Hdist_previous_value == 1.5) &&
                   (Hdist_next_value < next_hdptf_noremove_value) &&
                   is.null(merged_df1[[next_effdist_col]]) &&
                   Hdist_block1_value  > previous_hdptf_noremove_value)) {

          merged_df1[sufix_block1] <- dptf_sum + effdist_sum + merged_df1[[previous_effdist_col]]
          merged_df1[[previous_effdist_col]] <- NULL
        }}


      cols_to_check <- c(previous_effdist_col, Hdist_previous_colname, Hdist_next_colname, previous_hdptf_noremove_col, next_hdptf_noremove_col)
      if (all(cols_to_check %in% colnames(merged_df1))) {

        if (isTRUE(first_consec_suffix != first(lad_suffixes) &&
                   (Hdist_previous_value > previous_hdptf_noremove_value) &&
                   (Hdist_next_value > next_hdptf_noremove_value))) {

          merged_df1[sufix_block1] <- dptf_sum + effdist_sum +  merged_df1[[previous_effdist_col]]
          merged_df1[[previous_effdist_col]] <- NULL
        }}


      cols_to_check <- c(previous_effdist_col, Hdist_previous_colname, previous_hdptf_noremove_col)
      if (all(cols_to_check %in% colnames(merged_df1))) {

        if (isTRUE(first_consec_suffix != first(lad_suffixes) &&
                   (Hdist_previous_value > previous_hdptf_noremove_value) &&
                   (!Hdist_next_colname %in% colnames(merged_df1)))) {

          merged_df1[sufix_block1] <- dptf_sum + effdist_sum +  merged_df1[[previous_effdist_col]]
          merged_df1[[previous_effdist_col]] <- NULL
        }}


      cols_to_check <- c(next_effdist_col, Hdist_previous_colname, Hdist_next_colname, previous_hdptf_noremove_col,next_cbh_gt5_col)
      if (all(cols_to_check %in% colnames(merged_df1))) {

        if (isTRUE(first_consec_suffix != first(lad_suffixes) &&
                   Hdist_previous_value < previous_hdptf_noremove_value &&
                   Hdist_next_value < next_cbh_gt5_value)) {

          merged_df1[sufix_block1] <- dptf_sum + effdist_sum +  merged_df1[[next_effdist_col]]
          merged_df1[[next_effdist_col]] <- NULL
        }}


      cols_to_check <- c(next_effdist_col, previous_effdist_col, Hdist_previous_colname , Hdist_next_colname, previous_hdptf_noremove_col, next_hdptf_noremove_col)
      if (all(cols_to_check %in% colnames(merged_df1))) {

        if (isTRUE(first_consec_suffix != first(lad_suffixes) && (Hdist_previous_value  > previous_hdptf_noremove_value) && (Hdist_next_value < next_hdptf_noremove_value) && first_Hcbh_value > 1.5)) {

          merged_df1[sufix_block1] <- dptf_sum + effdist_sum + merged_df1[[previous_effdist_col]] + merged_df1[[next_effdist_col]]
          merged_df1[[previous_effdist_col]] <- NULL
          merged_df1[[next_effdist_col]] <- NULL

        }}

      #####################3
      # Extract column names from block_values
      cols_to_remove <- names(block_values)

      # Remove these columns from merged_df1
      merged_df1 <- merged_df1[, !names(merged_df1) %in% cols_to_remove]
      merged_df1<- get_renamed_df (merged_df1)

      ######################################################
      ########## MATCHING THE HDIST COLUMNS WITH HCBH COLUMNS  ##############################################

      hcbh_cols <- grep("^Hcbh[0-9]+$", names(merged_df1), value=T)
      hcbh_vals <-merged_df1[,hcbh_cols]

      hdist_cols <- grep("^Hdist[0-9]+$", names(merged_df1), value=T)
      hdist_vals <-merged_df1[,hdist_cols]

      # Include additional columns you want to keep (modify this according to your needs)
      other_cols <- merged_df1[, -which(names(merged_df1) %in% c(hdist_cols, hcbh_cols))]

      # Determine the number of "Hdist" columns to keep
      hdist_to_keep <- min(length(hdist_cols), length(hcbh_cols))

      # Select the columns to keep
      merged_df2 <- merged_df1[, c(hcbh_cols, hdist_cols[1:hdist_to_keep])]
      merged_df1 <- data.frame(merged_df2,other_cols)


      ########################################################
      # Replace Hdist columns where the value is greater than the corresponding Hcbh column minus one
      ########################################################

      # Identify Hcbh and Hdist columns
      hcbh_cols <- grep("^Hcbh[0-9]+$", names(merged_df1), value = TRUE)
      hdist_cols <- grep("^Hdist[0-9]+$", names(merged_df1), value = TRUE)

      # Determine the number of Hdist columns
      num_hdist_cols <- length(hdist_cols)

      # Determine the number of Hcbh columns
      num_hcbh_cols <- length(hcbh_cols)

      # If there are fewer Hcbh columns than Hdist columns, add the last Hcbh value as a new column
      if (num_hcbh_cols < num_hdist_cols) {
        last_hcbh_value <- tail(merged_df1[[tail(hcbh_cols, 1)]], 1)
        new_hcbh_col <- rep(last_hcbh_value, each = nrow(merged_df1))
        col_name <- paste0("Hcbh", num_hcbh_cols + 1)
        merged_df1[[col_name]] <- new_hcbh_col
      }

      # Replace Hdist columns where the value is different from the corresponding Hcbh column minus one
      for (i in 1:num_hdist_cols) {
        condition <- merged_df1[, hdist_cols[i]] != (merged_df1[[paste0("Hcbh", i)]] - 1)
        merged_df1[, hdist_cols[i]] <- ifelse(condition, merged_df1[[paste0("Hcbh", i)]] - 1, merged_df1[, hdist_cols[i]])
      }

      merged_df1<- get_renamed_df (merged_df1)


      }

    #########################################################
    ### NON-CONSECUTIVE SUFFIXES TO REMOVE
    #########################################################


    if (identical(merged_df1,merged_df) || length(suffixes_to_remove) > 0 ) {

      merged_df1<-remove_no_flayer_noconsec (merged_df1, threshold=10)
      merged_df1<-remove_no_flayer_noconsec (merged_df1, threshold=10)
      merged_df1<-remove_no_flayer_noconsec (merged_df1, threshold=10)
    }

  }

    ########################################################################

    if ((!"Hdist1" %in% names(merged_df1)  && (merged_df1$Hcbh1 <= 2.5))) {
      merged_df1$Hdist1 <-(merged_df1$Hcbh1- 1)
    }

    if ("Hdist1" %in% names(merged_df1) && (merged_df1$Hcbh1 == 1.5) && (merged_df1$Hdist1 > merged_df1$Hcbh1)) {
      merged_df1$Hdist0 <-0.5
    }
    if ("Hdist1" %in% names(merged_df1) && (merged_df1$Hcbh1 == 2.5) && (merged_df1$Hdist1 > merged_df1$Hcbh1)) {
      merged_df1$Hdist0 <-1.5
    }

    Hdist_cols <- grep("^Hdist[0-9]+$", names(merged_df1))

    # Sort effdist columns numerically
    Hdist_cols <- Hdist_cols[order(as.numeric(gsub("Hdist", "", names(merged_df1)[Hdist_cols])))]

    # Determine the new numeric suffix
    new_suffix <- seq_along(Hdist_cols)

    # Generate the new column names
    new_col_name <- paste0("Hdist", new_suffix)

    # Rename effdist columns directly using indexing
    names(merged_df1)[Hdist_cols] <- new_col_name

    all_cols <- names(merged_df1)
    numeric_suffix <- as.numeric(gsub("[^0-9]", "", all_cols))
    ordered_cols <- all_cols[order(numeric_suffix)]

    # Reorder the columns in df6a
    merged_df1 <- merged_df1[, ordered_cols]

    ######################################################
    ########## MATCHING THE HDIST COLUMNS WITH HCBH COLUMNS  ##############################################

    hcbh_cols <- grep("^Hcbh[0-9]+$", names(merged_df1), value=T)
    hcbh_vals <-merged_df1[,hcbh_cols]

    hdist_cols <- grep("^Hdist[0-9]+$", names(merged_df1), value=T)
    hdist_vals <-merged_df1[,hdist_cols]

    # Include additional columns you want to keep (modify this according to your needs)
    other_cols <- merged_df1[, -which(names(merged_df1) %in% c(hdist_cols, hcbh_cols))]

    # Determine the number of "Hdist" columns to keep
    hdist_to_keep <- min(length(hdist_cols), length(hcbh_cols))

    # Select the columns to keep
    merged_df2 <- merged_df1[, c(hcbh_cols, hdist_cols[1:hdist_to_keep])]
    merged_df1 <- data.frame(merged_df2,other_cols)


    ########################################################
    # Replace Hdist columns where the value is greater than the corresponding Hcbh column minus one
    ########################################################

    # Identify Hcbh and Hdist columns
    hcbh_cols <- grep("^Hcbh[0-9]+$", names(merged_df1), value = TRUE)
    hdist_cols <- grep("^Hdist[0-9]+$", names(merged_df1), value = TRUE)

    # Determine the number of Hdist columns
    num_hdist_cols <- length(hdist_cols)

    # Determine the number of Hcbh columns
    num_hcbh_cols <- length(hcbh_cols)

    # If there are fewer Hcbh columns than Hdist columns, add the last Hcbh value as a new column
    if (num_hcbh_cols < num_hdist_cols) {
      last_hcbh_value <- tail(merged_df1[[tail(hcbh_cols, 1)]], 1)
      new_hcbh_col <- rep(last_hcbh_value, each = nrow(merged_df1))
      col_name <- paste0("Hcbh", num_hcbh_cols + 1)
      merged_df1[[col_name]] <- new_hcbh_col
    }

    # Replace Hdist columns where the value is different from the corresponding Hcbh column minus one
    for (i in 1:num_hdist_cols) {
      condition <- merged_df1[, hdist_cols[i]] != (merged_df1[[paste0("Hcbh", i)]] - 1)
      merged_df1[, hdist_cols[i]] <- ifelse(condition, merged_df1[[paste0("Hcbh", i)]] - 1, merged_df1[, hdist_cols[i]])
    }

    merged_df1<- get_renamed_df (merged_df1)

    if ((!"Hdist1" %in% names(merged_df1)  && (merged_df1$Hcbh1 <= 2.5))) {
      merged_df1$Hdist1 <-(merged_df1$Hcbh1- 1)
    }

    if ("Hdist1" %in% names(merged_df1)) {
      merged_df1$Hdist1 <-(merged_df1$Hcbh1- 1)
    }

    ########################################################

    effdist_cols <- grep("^effdist[0-9]+$", names(merged_df1), value=T)

    if ((length(hdist_cols) > length(effdist_cols)) && (merged_df1$Hcbh1 == 1.5 && merged_df1$Hdist1 > merged_df1$Hcbh1)) {
      merged_df1$Hdist1 <- 0.5
    }

    if(merged_df1$Hcbh1 == 1.5 && merged_df1$Hdist1 > merged_df1$Hcbh1) {
      merged_df1$Hdist0 <- 0.5

      # Identify effdist columns
      Hdist_cols <- grep("^Hdist[0-9]+$", names(merged_df1))

      # Sort effdist columns numerically
      Hdist_cols <- Hdist_cols[order(as.numeric(gsub("Hdist", "", names(merged_df1)[Hdist_cols])))]

      # Determine the new numeric suffix
      new_suffix <- seq_along(Hdist_cols)

      # Generate the new column names
      new_col_name <- paste0("Hdist", new_suffix)

      # Rename effdist columns directly using indexing
      names(merged_df1)[Hdist_cols] <- new_col_name

      all_cols <- names(merged_df1)
      numeric_suffix <- as.numeric(gsub("[^0-9]", "", all_cols))
      ordered_cols <- all_cols[order(numeric_suffix)]

      # Reorder the columns in df6a
      merged_df1 <- merged_df1[, ordered_cols]
    }

    dist_cols <- grep("^dist[0-9]+$", names(merged_df1), value=T)
    merged_df1 <- merged_df1[, !names(merged_df1) %in% dist_cols]
    merged_df1<- get_renamed_df (merged_df1)

    ##########################################################################

    if ("effdist1" %in% colnames(merged_df1) && merged_df1$Hcbh1 <= 2.5 && merged_df1$effdist1 > 0) {
      merged_df1$effdist0 <- 0}

    if (!"effdist1" %in% colnames(merged_df1) && merged_df1$Hcbh1 <= 2.5) {
      merged_df1$effdist0 <- 0}

    # Identify effdist columns
    effdist_cols <- grep("^effdist[0-9]+$", names(merged_df1))

    # Sort effdist columns numerically
    effdist_cols <- effdist_cols[order(as.numeric(gsub("effdist", "", names(merged_df1)[effdist_cols])))]

    # Determine the new numeric suffix
    new_suffix <- seq_along(effdist_cols)

    # Generate the new column names
    new_col_name <- paste0("effdist", new_suffix)

    # Rename effdist columns directly using indexing
    names(merged_df1)[effdist_cols] <- new_col_name

    all_cols <- names(merged_df1)
    numeric_suffix <- as.numeric(gsub("[^0-9]", "", all_cols))
    ordered_cols <- all_cols[order(numeric_suffix)]

    # Reorder the columns in df6a
    merged_df1 <- merged_df1[, ordered_cols]


 } else {

    merged_df1 <- merged_df

    }

  ##########################################################################
  # IF NOT ANY CHANGE
  ##########################################################################

  if (identical(merged_df, merged_df1)) {


    dist_cols <- grep("^dist\\d+$", names(merged_df1), value = TRUE)
    hdist_cols <- grep("^Hdist", names(merged_df1), value = TRUE)

    if (all(dist_cols %in% colnames (merged_df1)) && all(hdist_cols %in% colnames (merged_df1))){
      merged_df1 <- merged_df1 %>%
        dplyr::select(-starts_with("dist"))

    }

    merged_df1<- get_renamed_df (merged_df1)

  }

  ####################################

  #if(!all(cols_to_append  %in% colnames(merged_df1))){
  # merged_df1<- cbind( merged_df1, append_df) }

  #if(!all(cols_to_append  %in% colnames(merged_df))){
  #  merged_df<- cbind( merged_df, append_df)}

  if (!("effdist1" %in% colnames(merged_df1))) {

    merged_df1$effdist1 <- 0

    # Identify effdist columns
    effdist_cols <- grep("^effdist[0-9]+$", names(merged_df1))

    # Sort effdist columns numerically
    effdist_cols <- effdist_cols[order(as.numeric(gsub("effdist", "", names(merged_df1)[effdist_cols])))]

    # Determine the new numeric suffix
    new_suffix <- seq_along(effdist_cols)

    # Generate the new column names
    new_col_name <- paste0("effdist", new_suffix)

    # Rename effdist columns directly using indexing
    names(merged_df1)[effdist_cols] <- new_col_name
  }

  all_effdist_cols <- grep("^effdist\\d+$", names(merged_df1), value = TRUE)
  dist_cols <- grep("^dist\\d+$", names(merged_df1), value = TRUE)
  hdist_cols <- grep("^Hdist\\d+$", names(merged_df1), value = TRUE)
  hcbh_cols <- grep("^Hcbh[0-9]+$", names(merged_df1), value = TRUE)

  if (length(all_effdist_cols) > 0 && length(Hdist_cols) > 0) {

    if (length(dist_cols) > 0) {

      merged_df1 <- merged_df1[, -which(names(merged_df1) %in% dist_cols)]

      merged_df1<- get_renamed_df (merged_df1)
    }

    #################################################
    if(("Hdist1" %in% names(merged_df1)) && sum(grepl("^Hcbh[0-9]+$", hcbh_cols)) == 1 && merged_df1$Hdist1 > merged_df1$Hcbh1) {
      merged_df1$Hdist1 <- merged_df1$Hcbh1 -1
    }

    if(("Hdist1" %in% names(merged_df1)) && sum(grepl("^Hcbh[0-9]+$", hcbh_cols)) > 1 && merged_df1$Hdist1 > merged_df1$Hcbh1) {
      merged_df1$Hdist0 <- merged_df1$Hcbh1 -1
    }

    #################################################

    Hdist_cols <- grep("^Hdist[0-9]+$", names(merged_df1))

    # Sort effdist columns numerically
    Hdist_cols <- Hdist_cols[order(as.numeric(gsub("Hdist", "", names(merged_df1)[Hdist_cols])))]

    # Determine the new numeric suffix
    new_suffix <- seq_along(Hdist_cols)

    # Generate the new column names
    new_col_name <- paste0("Hdist", new_suffix)

    # Rename effdist columns directly using indexing
    names(merged_df1)[Hdist_cols] <- new_col_name

    all_cols <- names(merged_df1)
    numeric_suffix <- as.numeric(gsub("[^0-9]", "", all_cols))
    ordered_cols <- all_cols[order(numeric_suffix)]

    # Reorder the columns in df6a
    merged_df1 <- merged_df1[, ordered_cols]

    #################################################

    if (("effdist1" %in% names(merged_df1) && !is.na(merged_df1$effdist1))  && sum(grepl("^Hcbh[0-9]+$", hcbh_cols)) == 1 &&
        ("Hdist1" %in% names(merged_df1))) {
      if(merged_df1$effdist1 != (merged_df1$Hdist1-0.5)) {
        merged_df1$effdist1 <-(merged_df1$Hdist1-0.5)
      }
    }

    if (("effdist1" %in% names(merged_df1) && !is.na(merged_df1$effdist1)) && (!"Hdist1" %in% names(merged_df1))) {
      merged_df1$Hdist1 <- merged_df1$Hcbh1 -1
      merged_df1$effdist1 <-(merged_df1$Hdist1-0.5)
    }

    if (("effdist1" %in% names(merged_df1) && !is.na(merged_df1$effdist1))  && sum(grepl("^Hcbh[0-9]+$", hcbh_cols)) > 1 &&
        (merged_df1$Hcbh1 <= 2.5) && merged_df1$effdist1 != 0 ) {
      merged_df1$effdist0 <-(merged_df1$Hdist1-0.5)

    # Identify effdist columns
    effdist_cols <- grep("^effdist[0-9]+$", names(merged_df1))

    # Sort effdist columns numerically
    effdist_cols <- effdist_cols[order(as.numeric(gsub("effdist", "", names(merged_df1)[effdist_cols])))]

    # Determine the new numeric suffix
    new_suffix <- seq_along(effdist_cols)

    # Generate the new column names
    new_col_name <- paste0("effdist", new_suffix)

    # Rename effdist columns directly using indexing
    names(merged_df1)[effdist_cols] <- new_col_name

    all_cols <- names(merged_df1)
    numeric_suffix <- as.numeric(gsub("[^0-9]", "", all_cols))
    ordered_cols <- all_cols[order(numeric_suffix)]

    # Reorder the columns in df6a
    merged_df1 <- merged_df1[, ordered_cols]

    }}

  if ((!"effdist1" %in% names(merged_df1)  && (merged_df1$Hcbh1 <= 2.5))) {
    merged_df1$effdist1 <-(merged_df1$Hdist1-0.5)
}

  ##########################################################################

  if(merged_df1$Hcbh1 == 1.5 && merged_df1$Hdist1 > merged_df1$Hcbh1) {
    merged_df1$Hdist0 <- 0.5

    # Identify effdist columns
    Hdist_cols <- grep("^Hdist[0-9]+$", names(merged_df1))

    # Sort effdist columns numerically
    Hdist_cols <- Hdist_cols[order(as.numeric(gsub("Hdist", "", names(merged_df1)[Hdist_cols])))]

    # Determine the new numeric suffix
    new_suffix <- seq_along(Hdist_cols)

    # Generate the new column names
    new_col_name <- paste0("Hdist", new_suffix)

    # Rename effdist columns directly using indexing
    names(merged_df1)[Hdist_cols] <- new_col_name

    all_cols <- names(merged_df1)
    numeric_suffix <- as.numeric(gsub("[^0-9]", "", all_cols))
    ordered_cols <- all_cols[order(numeric_suffix)]

    # Reorder the columns in df6a
    merged_df1 <- merged_df1[, ordered_cols]
  }

  ##########################################################################

  if(merged_df1$Hcbh1 == 1.5 && merged_df1$effdist1 > 0) {
    merged_df1$effdist0 <- 0


    # Identify effdist columns
    effdist_cols <- grep("^effdist[0-9]+$", names(merged_df1))

    # Sort effdist columns numerically
    effdist_cols <- effdist_cols[order(as.numeric(gsub("effdist", "", names(merged_df1)[effdist_cols])))]

    # Determine the new numeric suffix
    new_suffix <- seq_along(effdist_cols)

    # Generate the new column names
    new_col_name <- paste0("effdist", new_suffix)

    # Rename effdist columns directly using indexing
    names(merged_df1)[effdist_cols] <- new_col_name

    all_cols <- names(merged_df1)
    numeric_suffix <- as.numeric(gsub("[^0-9]", "", all_cols))
    ordered_cols <- all_cols[order(numeric_suffix)]

    # Reorder the columns in df6a
    merged_df1 <- merged_df1[, ordered_cols]
  }

  ##########################################################################

  if(merged_df$Hcbh1 == 1.5 && merged_df$Hdist1 > merged_df$Hcbh1) {
    merged_df$Hdist0 <- 0.5

    # Identify effdist columns
    Hdist_cols <- grep("^Hdist[0-9]+$", names(merged_df))

    # Sort effdist columns numerically
    Hdist_cols <- Hdist_cols[order(as.numeric(gsub("Hdist", "", names(merged_df)[Hdist_cols])))]

    # Determine the new numeric suffix
    new_suffix <- seq_along(Hdist_cols)

    # Generate the new column names
    new_col_name <- paste0("Hdist", new_suffix)

    # Rename effdist columns directly using indexing
    names(merged_df)[Hdist_cols] <- new_col_name

    all_cols <- names(merged_df)
    numeric_suffix <- as.numeric(gsub("[^0-9]", "", all_cols))
    ordered_cols <- all_cols[order(numeric_suffix)]

    # Reorder the columns in df6a
    merged_df <- merged_df[, ordered_cols]

  }

  ##########################################################################

  if(merged_df$Hcbh1 == 1.5 && merged_df$effdist1 > 0) {
    merged_df$effdist0 <- 0

    # Identify effdist columns
    effdist_cols <- grep("^effdist[0-9]+$", names(merged_df))

    # Sort effdist columns numerically
    effdist_cols <- effdist_cols[order(as.numeric(gsub("effdist", "", names(merged_df)[effdist_cols])))]

    # Determine the new numeric suffix
    new_suffix <- seq_along(effdist_cols)

    # Generate the new column names
    new_col_name <- paste0("effdist", new_suffix)

    # Rename effdist columns directly using indexing
    names(merged_df)[effdist_cols] <- new_col_name

    all_cols <- names(merged_df)
    numeric_suffix <- as.numeric(gsub("[^0-9]", "", all_cols))
    ordered_cols <- all_cols[order(numeric_suffix)]

    # Reorder the columns in df6a
    merged_df <- merged_df[, ordered_cols]
  }


  ########  MATCHING EFFDIST COLUMNS TO HCBH COLUMNS###################

  effdist_cols <- grep("^effdist\\d+$", names(merged_df1), value = TRUE)

  if(length(effdist_cols) > 0) {

  hcbh_cols <- grep("^Hcbh[0-9]+$", names(merged_df1), value = TRUE)
  hcbh_vals <- merged_df1[, hcbh_cols]

  hdist_cols <- grep("^Hdist[0-9]+$", names(merged_df1), value = TRUE)
  hdist_vals <- merged_df1[, hdist_cols]

  # Include additional columns
  other_cols <- merged_df1[, -which(names(merged_df1) %in% effdist_cols)]

  # Combine all the columns to keep
  cols_to_keep1 <- c(effdist_cols[1:min(length(effdist_cols), length(hcbh_cols))])
  cols_to_keep_df <- merged_df1 %>% dplyr::select(all_of(cols_to_keep1))


  if (all(!is.na(cols_to_keep1)) && exists("cols_to_keep_df")) {
    # Select the columns to keep
    merged_df1 <- data.frame(other_cols, cols_to_keep_df)
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

  cbh_columns1 <- grep("^Hcbh\\d+$", names(all_LAD))  # Extracts numeric cbh columns
  nlayers1 <- sum(sapply(all_LAD[, cbh_columns1], is.numeric))  # Count numeric values

  cbh_columns2 <- grep("^Hcbh\\d+$", names(effective_LAD))  # Extracts numeric cbh columns
  nlayers2 <- sum(sapply(effective_LAD[, cbh_columns2], is.numeric))  # Count numeric values

  # Add a new column to tree1 with the calculated nlayers1 value
  all_LAD$nlayers <- nlayers1
  effective_LAD$nlayers <- nlayers2


  # Return them in a list
  return(list(df1 = all_LAD, df2 = effective_LAD))

}


