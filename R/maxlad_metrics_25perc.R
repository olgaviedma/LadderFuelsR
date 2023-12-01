#' Fuels LAD percentage and canopy base height (CBH) based on maximum LAD percentage (distances > 1 m)
#'
#' @description
#' This function calculates the percentage of leaf area density (LAD) within each fuel layer (first output),
#' and removes those fuel layers with LAD percentage less than 25, recalculating the distances of the remaining ones.
#' It determines the canopy base height (CBH) as the fuel layer with the highest LAD percentage (second output).
#'
#' @usage
#' get_layers_lad(LAD_profiles, effective_distances)
#'
#' @param LAD_profiles
#' Original tree Leaf Area Index (LAD) profile (output of [lad.profile()] function in the \emph{leafR} package).
#' An object of the class text.
#'
#' @param effective_distances
#' Tree metrics of fuel layers separated by distances greater than 1 m (output of [get_effective_gap()] function).
#' An object of the class text.
#'
#' @return
#' A data frame identifying the canopy base height (CBH) of the fuel layer with maximum LAD percentage and other fuel layers with their corresponding LAD percentage.
#'
#' @author
#' Olga Viedma, Carlos Silva and JM Moreno
#'
#' @details
#' # List of tree metrics:
#' \itemize{
#'   \item treeID: tree ID with strings and numeric values
#'   \item treeID1: tree ID with only numeric values
#'   \item Hdist - Height of the distance between consecutive fuel layers (m)
#'   \item Hcbh - Height of the base of each fuel layer (m)
#'   \item effdist - Distance between consecutive fuel layers (m)
#'   \item dptf - Depth of fuel layers (m) after removing distances equal to 1 m
#'   \item Hdptf - Height of the depth of fuel layers (m) after removing distances equal to 1 m
#'   \item Hcbh_Hdptf - Percentage of LAD values comprised in each fuel layer
#'   \item maxlad_Hcbh - Height of the CBH of the segmented tree based on the maximum LAD percentage
#'   \item max_Hcbh - Height of the CBH of the segmented tree based on the maximum distance found in its profile
#'   \item last_Hcbh - Height of the CBH of the segmented tree based on the last distance found in its profile
#'   \item maxlad_ - Values of distance and fuel depth and their corresponding heights at the maximum LAD percentage
#'   \item max_ - Values of distance and fuel depth and their corresponding heights at the maximum distance found in the tree profile
#'   \item last_ - Values of distance and fuel depth and their corresponding heights at the last distance found in its profile
#'   \item max_height - Maximum height of the tree profile
#' }
#'
#' @examples
#' ## Not run:
#' library(dplyr)
#' library(magrittr)
#' library(gdata)
#'
#' # LAD profiles derived from normalized ALS data after applying [lad.profile()] function
#' LAD_profiles$treeID <- factor(LAD_profiles$treeID)
#'
#' ## Not run:
#' # Load or create the effective_distances object
#' if (interactive()) {
#'   effective_distances <- get_effective_gap()
#'   LadderFuelsR::effective_distances$treeID <- factor(LadderFuelsR::effective_distances$treeID)
#'
#' trees_name1 <- as.character(effective_distances$treeID)
#' trees_name2 <- factor(unique(trees_name1))
#'
#' LAD_metrics1 <- list()
#' LAD_metrics2 <- list()
#'
#' for (i in levels(trees_name2)) {
#'   # Filter data for each tree
#'   tree1 <- LAD_profiles |> dplyr::filter(treeID == i)
#'   tree2 <- effective_distances |> dplyr::filter(treeID == i)
#'
#'   # Get LAD metrics for each tree
#'   LAD_metrics <- get_layers_lad(tree1, tree2)
#'   LAD_metrics1[[i]] <- LAD_metrics$df1
#'   LAD_metrics2[[i]] <- LAD_metrics$df2
#' }
#'
#' all_LAD <- dplyr::bind_rows(LAD_metrics1)
#' effective_LAD <- dplyr::bind_rows(LAD_metrics2)
#'}
#' ## End(Not run)
#'
#' @export get_layers_lad
#' @importFrom dplyr select_if group_by summarise mutate arrange rename rename_with filter slice ungroup
#' @importFrom magrittr %>%
#' @importFrom SSBtools RbindAll
#' @importFrom gdata startsWith
#' @include gap_fbh.R
#' @include distances_calculation.R
#' @include depths_calculation.R
#' @include corrected_base_heights.R
#' @include corrected_depth.R
#' @include corrected_distances.R
#' @seealso \code{\link{get_renamed_df}}
get_layers_lad <- function(LAD_profiles, effective_distances) {

  df_orig <- LAD_profiles
  effectiv_gaps<- effective_distances

  df_effective1 <- effectiv_gaps[, !apply(effectiv_gaps, 2, function(x) all(is.na(x)))]
  treeID<-"treeID"
  treeID1<-"treeID1"
  print(paste("treeID:", df_effective1[[treeID]]))  # Debugging line

  ######################################
  Hcbh_cols <- grep("^Hcbh\\d+$", names(df_effective1), value = TRUE)
  Hdist_cols <- grep("^Hdist\\d+$", names(df_effective1), value = TRUE)
  Hdptf_cols <- grep("^Hdptf\\d+$", names(df_effective1), value = TRUE)
  dptf_cols <- grep("^dptf\\d+$", names(df_effective1), value = TRUE)
  effdist_cols <- grep("^effdist\\d+$", names(df_effective1), value = TRUE)

  # Select columns with prefixes: last_, max_, treeID
  cols_to_append <- grep("^(last_|max_|treeID)", names(df_effective1), value = TRUE)
  append_df <- df_effective1[ , cols_to_append]

  # Exclude columns with prefixes: last_, max_, treeID
  cols_to_exclude <- grep("^(last_|max_|treeID)", names(df_effective1), value = TRUE)
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
    end_height <- as.numeric(df_effective1[1, hdepth_colname])

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
  cols_to_remove <- sapply(merged_df1[, lad_columns3], function(col) any(col < 25))

  all_effdist_cols <- grep("^effdist\\d+$", names(merged_df1), value = TRUE)
  all_effdist_suffixes <- as.numeric(stringr::str_extract(all_effdist_cols, "\\d+$"))
  suffixes_to_remove <- sort(as.numeric(stringr::str_extract(names(cols_to_remove[cols_to_remove]), "\\d+$")))

  #print(paste("Last LAD Suffix:", last_lad_suffix))  # Debugging line

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

    lad_no_remove <- sapply(merged_df1[, lad_columns2], function(col) any(col > 25))
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

    # Find the columns with value less than 5
    cols_to_remove <- sapply(merged_df1[, lad_columns2], function(col) any(col < 25))
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


      # && all (sapply(all_effdist_values, function(x) x != "1")

      #########################

      cols_to_check <- c(Hdist_next_colname, first_cbh_gt5_col,first(Hcbh_cols),Hdist_block2,all_effdist_cols)
      if (all(cols_to_check %in% colnames(merged_df1))) {

        if (first_consec_suffix == first(lad_suffixes) &&
            !is.na(Hdist_next_value) &&
            !is.na(first_cbh_gt5_value) &&
            first_cbh_gt5_value  > 1.5  &&
            Hcbh1 == 1.5  &&
            length(suffixes_no_remove)== 1 &&
            (Hdist_block2_value  < first_cbh_gt5_value) && all (sapply(all_effdist_values, function(x) x != "1"))) { #&& length(suffixes_no_remove==1)

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
            Hcbh1 == 1.5) { #&& length(suffixes_no_remove==1)

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
            (Hdist_block2_value  < first_cbh_gt5_value) && any (sapply(all_effdist_values, function(x) x == "1"))) { #&& length(suffixes_no_remove==1)

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
        }}

      #######################

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


      ####################################

      all_effdist_cols <- grep("^effdist\\d+$", names(merged_df1), value = TRUE)

      if (length(all_effdist_cols) > 0) {
        for (col in all_effdist_cols) {
          if (any(merged_df1[, col] == 1)) {
            merged_df1[, col] <- NULL
          }}
      }


      dist_cols <- grep("^dist\\d+$", names(merged_df1), value = TRUE)
      hdist_cols <- grep("^Hdist", names(merged_df1), value = TRUE)

      if (all(dist_cols %in% colnames (merged_df1)) && all(hdist_cols %in% colnames (merged_df1))){
        merged_df1 <- merged_df1 %>%
          dplyr::select(-starts_with("dist"))

      }
      merged_df1<- get_renamed_df (merged_df1)


      if (length(hdist_cols) > length(all_effdist_cols) && merged_df1$Hcbh1 == 1.5) {
        merged_df1 <- merged_df1[, !names(merged_df1) %in% last(hdist_cols)]
      }

      merged_df1<- get_renamed_df (merged_df1)



      ####################################

      lad_columns2 <- grep("^Hcbh\\d+_Hdptf\\d+$", names(merged_df1),value=TRUE)
      lad_suffixes <- as.numeric(str_extract(lad_columns2, "\\d+$"))

      lad_no_remove <- sapply(merged_df1[, lad_columns2], function(col) any(col > 25))
      suffixes_no_remove <- sort(as.numeric(stringr::str_extract(names(lad_no_remove[lad_no_remove]), "\\d+$")))

      colnames_to_extract <- names(lad_no_remove)[lad_no_remove]
      # Extract the values from the desired columns
      lad_values_noremove <- merged_df1[, colnames_to_extract, drop=FALSE]

      pattern <- paste0(suffixes_no_remove, "$")

      # Use sapply to get matches for each pattern
      matching_cols_list <- sapply(pattern, function(pat) grep(pat, names(merged_df1), value=TRUE))

      # Unlist and make unique the result
      matching_cols <- unique(unlist(matching_cols_list))

      # Extract column name
      first_cbh_gt5_col <- first(grep("^Hcbh\\d+$", matching_cols, value = TRUE))


      if (first_cbh_gt5_col %in% colnames(merged_df1)) {
        first_cbh_gt5_value <- merged_df1[, first_cbh_gt5_col]
      } else {
        first_cbh_gt5_value <- NA  # or some other default value
      }

      last_cbh_gt5_col <- last(grep("^Hcbh\\d+$", matching_cols, value = TRUE))

      if (last_cbh_gt5_col %in% colnames(merged_df1)) {
        last_cbh_gt5_value <- merged_df1[, last_cbh_gt5_col]
      } else {
        last_cbh_gt5_value <- NA  # or some other default value
      }

      #########################################
      # Find the columns with value less than 5
      cols_to_remove <- sapply(merged_df1[, lad_columns2], function(col) any(col < 25))
      suffixes_to_remove <- sort(as.numeric(stringr::str_extract(names(cols_to_remove[cols_to_remove]), "\\d+$")))

      pattern2 <- paste0(suffixes_to_remove, "$")
      # Get column names matching the pattern
      matching_cols2 <- unique(unlist(sapply(pattern2, function(pat) {
        grep(pat, names(merged_df1), value=TRUE)
      })))

      cols_to_remove_noeffdist <- matching_cols2[!grepl("^(effdist|treeID)", matching_cols2)]

      # Subset the dataframe using matching columns
      remove_noeffdist_values <- merged_df1[, cols_to_remove_noeffdist]

      if (length(suffixes_to_remove) > 0) {

        first_no_consenc_suffix <- min(suffixes_to_remove)
        last_no_consenc_suffix <- max(suffixes_to_remove)
        next_column_suffix <- last_no_consenc_suffix + 1
        previous_column_suffix <- first_no_consenc_suffix - 1
        first_col_name <- paste0("effdist", first_no_consenc_suffix)
        last_col_name <- paste0("effdist", last_no_consenc_suffix)
        last_lad_suffix<-last(lad_suffixes)


        if (last_col_name %in% colnames(merged_df1)) {
          last_col_value <- merged_df1[, last_col_name]
        } else {
          last_col_value <- NA  # or some other default value
        }


        next_effdist_col <- paste0("effdist", next_column_suffix)
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


        Hdist_next_colname <- paste0("Hdist", next_column_suffix)
        Hdptf_next_colname <- paste0("Hdptf", next_column_suffix)

        if (Hdist_next_colname %in% colnames(merged_df1)) {
          Hdist_next_value <- merged_df1[, Hdist_next_colname]
        } else {
          Hdist_next_value <- NA  # or some other default value
        }

        if (Hdptf_next_colname %in% colnames(merged_df1)) {
          Hdptf_next_value <- merged_df1[, Hdptf_next_colname]
        } else {
          Hdptf_next_value <- NA  # or some other default value
        }

        Hdist_previous_colname <- paste0("Hdist", previous_column_suffix)
        Hdptf_previous_colname <- paste0("Hdptf", previous_column_suffix)

        if (any(Hdist_previous_colname %in% colnames(merged_df1))) {
          Hdist_previous_value <- merged_df1[, Hdist_previous_colname]
        } else {
          Hdist_previous_value <- NA  # or some other default value
        }

        if (Hdptf_previous_colname %in% colnames(merged_df1)) {
          Hdptf_previous_value <- merged_df1[, Hdptf_previous_colname]
        } else {
          Hdptf_previous_value <- NA  # or some other default value
        }

        dist_col1 <- paste0("dist", first_no_consenc_suffix)

        if (dist_col1 %in% colnames(merged_df1)) {
          dist1_value <- merged_df1[, dist_col1]
        } else {
          dist1_value <- NA  # or some other default value
        }

        lad_columns2 <- grep("^Hcbh\\d+_Hdptf\\d+$", names(merged_df1),value=TRUE)
        lad_suffixes <- as.numeric(str_extract(lad_columns2, "\\d+$"))


        ####################################################

        if (any(suffixes_to_remove == 1)) {


          first_suffix_eq1 <- first(suffixes_to_remove[suffixes_to_remove == 1])
          suffixes_to_remove1 <- first_suffix_eq1

          pattern2 <- paste0(suffixes_to_remove1, "$")
          # Get column names matching the pattern
          matching_cols3 <- unique(unlist(sapply(pattern2, function(pat) {
            grep(pat, names(merged_df1), value=TRUE)
          })))
          cols_to_remove_noeffdist1 <- matching_cols3[!grepl("^(effdist|Hdist|treeID)", matching_cols3)]

          next_column_suffix <- suffixes_to_remove1 + 1
          next_col_name <- paste0("effdist", next_column_suffix)

          effdist_cols <- paste0("effdist", suffixes_to_remove1)

          if (!is.null(effdist_cols) && length(effdist_cols) > 0 && effdist_cols %in% colnames(merged_df1)) {

            effdist_values <- sapply(suffixes_to_remove1, function(suf) merged_df1[[paste0("effdist", suf)]])
            dptf_values <- sapply(suffixes_to_remove1, function(suf) merged_df1[[paste0("dptf", suf)]])

            effdist_values <- sapply(effdist_values, function(x) if(is.null(x)) 0 else x)
            dptf_values <- sapply(dptf_values, function(x) if(is.null(x)) 0 else x)

            effdist_sum <- sum(effdist_values)
            dptf_sum <- sum(dptf_values)

            next_effdist_col <- paste0("effdist", next_column_suffix)

            if (next_effdist_col %in% colnames(merged_df1)) {
              next_effdist_value <- merged_df1[, next_effdist_col]
            } else {
              next_effdist_value <- NA  # or some other default value
            }

            next2_effdist_col <- paste0("effdist", next_column_suffix +1)

            if (next2_effdist_col %in% colnames(merged_df1)) {
              next2_effdist_value <- merged_df1[, next2_effdist_col]
            } else {
              next2_effdist_value <- NA  # or some other default value
            }

            Hdist_next_colname <- paste0("Hdist", next_column_suffix)

            if (Hdist_next_colname %in% colnames(merged_df1)) {
              Hdist_next_value <- merged_df1[, Hdist_next_colname]
            } else {
              Hdist_next_value <- NA  # or some other default value
            }

            Hdist_next2_colname <- paste0("Hdist", next_column_suffix +1)

            if (Hdist_next2_colname %in% colnames(merged_df1)) {
              Hdist_next2_value <- merged_df1[, Hdist_next2_colname]
            } else {
              Hdist_next2_value <- NA  # or some other default value
            }


            next_hdptf_noremove_col <- paste0("Hdptf", suffixes_to_remove1 + 1)

            if (next_hdptf_noremove_col %in% colnames(merged_df1)) {
              next_hdptf_noremove_value <- merged_df1[, next_hdptf_noremove_col]
            } else {
              next_hdptf_noremove_value <- NA  # or some other default value
            }

            all_effdist_cols <- grep("^effdist\\d+$", names(merged_df1), value = TRUE)
            all_effdist_values <-merged_df1[,all_effdist_cols]

            Hcbh_cols <- grep("^Hcbh\\d+$", names(merged_df1), value = TRUE)

            if (all(Hcbh_cols %in% colnames(merged_df1))) {
              Hcbh1 <- merged_df1[, first(Hcbh_cols)]
            } else {
              Hcbh1 <- NA  # or some other default value
            }

            lad_columns2 <- grep("^Hcbh\\d+_Hdptf\\d+$", names(merged_df1),value=TRUE)
            lad_suffixes <- as.numeric(str_extract(lad_columns2, "\\d+$"))

            ###############################################################################3


            cols_to_check <- c(first_cbh_gt5_col,dist_col1)
            if (all(cols_to_check %in% colnames(merged_df1))) {

              if (isTRUE(first_suffix_eq1 == first(lad_suffixes)) && first_cbh_gt5_value  > 1.5 &&
                  !next_effdist_col %in% colnames(merged_df1) && dist1_value > 1) {

                merged_df1[first_col_name] <- dptf_sum + effdist_sum
                merged_df1 <- merged_df1[, !(names(merged_df1) %in% cols_to_remove_noeffdist1)]
              }}



            cols_to_check <- c(first_cbh_gt5_col, first(Hcbh_cols),Hdist_next_colname)
            if (all(cols_to_check %in% colnames(merged_df1))) {

              if (isTRUE(first_suffix_eq1 == first(lad_suffixes)) &&
                  first_cbh_gt5_value  > 1.5 && Hcbh1 ==1.5 &&
                  !Hdist_next2_colname %in% colnames(merged_df1) &&
                  Hdist_next_value < first_cbh_gt5_value && all (sapply(all_effdist_values, function(x) x != "1"))) {

                merged_df1[first_col_name] <- dptf_sum + effdist_sum
                merged_df1 <- merged_df1[, !(names(merged_df1) %in% cols_to_remove_noeffdist1)]
              }}


            cols_to_check <- c(first_cbh_gt5_col,Hdist_next_colname)
            if (all(cols_to_check %in% colnames(merged_df1))) {


              if (isTRUE(first_suffix_eq1 == first(lad_suffixes) &&
                         first_cbh_gt5_value  > 1.5 &&
                         (Hdist_next_value > first_cbh_gt5_value) &&
                         !Hdist_next2_colname %in% colnames(merged_df1) &&
                         all(sapply(all_effdist_values, function(x) x != "1")) &&
                         length(lad_values_noremove) > 1) ) {

                merged_df1[first_col_name] <- dptf_sum + effdist_sum
                merged_df1 <- merged_df1[, !(names(merged_df1) %in% cols_to_remove_noeffdist1)]
              }}



            cols_to_check <- c(first_cbh_gt5_col, first(Hcbh_cols),next_effdist_col,Hdist_next2_colname)
            if (all(cols_to_check %in% colnames(merged_df1))) {

              if (isTRUE(first_suffix_eq1 == first(lad_suffixes)) &&
                  first_cbh_gt5_value > 1.5 && Hcbh1== 1.5 &&  ## first distance == 0.5
                  next_effdist_value==1 &&
                  Hdist_next2_value > first_cbh_gt5_value && any(sapply(all_effdist_values, function(x) x == "1"))) {

                merged_df1[first_col_name] <- dptf_sum + effdist_sum
                merged_df1 <- merged_df1[, !(names(merged_df1) %in% cols_to_remove_noeffdist1)]
              }}



            cols_to_check <- c(first_cbh_gt5_col, first(Hcbh_cols),next_effdist_col,Hdist_next2_colname)
            if (all(cols_to_check %in% colnames(merged_df1))) {

              if (isTRUE(first_suffix_eq1 == first(lad_suffixes)) &&
                  first_cbh_gt5_value > 1.5 && Hcbh1== 1.5 &&  ## first distance == 0.5
                  next_effdist_value > 1 &&
                  Hdist_next2_value > first_cbh_gt5_value && all (sapply(all_effdist_values, function(x) x != "1"))) {

                merged_df1[first_col_name] <- dptf_sum + effdist_sum
                merged_df1 <- merged_df1[, !(names(merged_df1) %in% cols_to_remove_noeffdist1)]
              }}


            cols_to_check <- c(first_cbh_gt5_col, first(Hcbh_cols))
            if (all(cols_to_check %in% colnames(merged_df1))) {

              if (isTRUE(first_suffix_eq1 == first(lad_suffixes)) &&
                  first_cbh_gt5_value > 1.5 && Hcbh1== 1.5 &&  ## first distance == 0.5
                  !next_effdist_col %in% colnames(merged_df1) &&
                  !Hdist_next2_colname %in% colnames(merged_df1) && all (sapply(all_effdist_values, function(x) x != "1"))) {

                merged_df1[first_col_name] <- dptf_sum + effdist_sum
                merged_df1 <- merged_df1[, !(names(merged_df1) %in% cols_to_remove_noeffdist1)]
              }}

            #################################

            cols_to_check <- c(first_cbh_gt5_col, first(Hcbh_cols),next_effdist_col,Hdist_next2_colname,next_hdptf_noremove_col)
            if (all(cols_to_check %in% colnames(merged_df1))) {

              if (isTRUE(first_suffix_eq1 == first(lad_suffixes)) &&
                  first_cbh_gt5_value  > 1.5 && Hcbh1 ==1.5 &&  ## first distance == 0.5
                  Hdist_next2_value < next_hdptf_noremove_value &&
                  length(lad_values_noremove) == 1) {

                merged_df1[first_col_name] <- dptf_sum + effdist_sum + merged_df1[[next_effdist_col]]
                merged_df1 <- merged_df1[, !(names(merged_df1) %in% cols_to_remove_noeffdist1)]
                merged_df1[[next_effdist_col]] <- NULL
              }}

            #################################

            cols_to_check <- c(first_cbh_gt5_col, first(Hcbh_cols),next_effdist_col,Hdist_next_colname, Hdist_next2_colname,next_hdptf_noremove_col)
            if (all(cols_to_check %in% colnames(merged_df1))) {

              if (isTRUE(first_suffix_eq1 == first(lad_suffixes) &&
                         first_cbh_gt5_value  > 1.5 && Hcbh1 ==1.5 &&  ## first distance == 0.5
                         (Hdist_next_value < first_cbh_gt5_value) &&
                         Hdist_next2_value > next_hdptf_noremove_value &&
                         any(sapply(all_effdist_values, function(x) x == "1")) &&
                         length(lad_values_noremove) > 1)) {

                merged_df1[first_col_name] <- dptf_sum + effdist_sum + merged_df1[[next_effdist_col]]
                merged_df1 <- merged_df1[, !(names(merged_df1) %in% cols_to_remove_noeffdist1)]
                merged_df1[[next_effdist_col]] <- NULL
              }}

            #################################

            cols_to_check <- c(first_cbh_gt5_col,next_effdist_col,Hdist_next_colname,Hdist_next2_colname,next_hdptf_noremove_col)
            if (all(cols_to_check %in% colnames(merged_df1))) {


              if (isTRUE(first_suffix_eq1 == first(lad_suffixes) &&
                         first_cbh_gt5_value  > 1.5 && Hcbh1 > 1.5 &&
                         (Hdist_next_value < first_cbh_gt5_value) &&
                         Hdist_next2_value > next_hdptf_noremove_value )) {

                merged_df1[first_col_name] <- dptf_sum + effdist_sum  + merged_df1[[next_effdist_col]]
                merged_df1 <- merged_df1[, !(names(merged_df1) %in% cols_to_remove_noeffdist1)]
                merged_df1[[next_effdist_col]] <- NULL
              }}


            ########################

            cols_to_check <- c(first_cbh_gt5_col,next_effdist_col,Hdist_next_colname,Hdist_next2_colname,next_hdptf_noremove_col)
            if (all(cols_to_check %in% colnames(merged_df1))) {


              if (isTRUE(first_suffix_eq1 == first(lad_suffixes) &&
                         first_cbh_gt5_value  > 1.5 &&
                         (Hdist_next_value < first_cbh_gt5_value) &&
                         Hdist_next2_value < next_hdptf_noremove_value )) {

                merged_df1[first_col_name] <- dptf_sum + effdist_sum  + merged_df1[[next_effdist_col]]
                merged_df1 <- merged_df1[, !(names(merged_df1) %in% cols_to_remove_noeffdist1)]
                merged_df1[[next_effdist_col]] <- NULL
              }}

            ########################

            cols_to_check <- c(first_cbh_gt5_col,next_effdist_col,Hdist_next_colname,next_effdist_col,Hcbh_cols)
            if (all(cols_to_check %in% colnames(merged_df1))) {


              if (isTRUE(first_suffix_eq1 == first(lad_suffixes) &&
                         first_cbh_gt5_value  > 1.5 && Hcbh1 == 1.5 &&
                         (Hdist_next_value < first_cbh_gt5_value) &&
                         !Hdist_next2_colname %in% colnames(merged_df1) &&
                         length(lad_values_noremove) == 1)) {

                merged_df1[first_col_name] <- dptf_sum + effdist_sum  + merged_df1[[next_effdist_col]]
                merged_df1 <- merged_df1[, !(names(merged_df1) %in% cols_to_remove_noeffdist1)]
                merged_df1[[next_effdist_col]] <- NULL

              }}


            ########################

            cols_to_check <- c(first_cbh_gt5_col,next_effdist_col,Hdist_next_colname,all_effdist_cols)
            if (all(cols_to_check %in% colnames(merged_df1))) {


              if (isTRUE(first_suffix_eq1 == first(lad_suffixes)) &&
                  first_cbh_gt5_value  > 1.5 &&
                  (Hdist_next_value < first_cbh_gt5_value) &&
                  (!Hdist_next2_colname %in% colnames(merged_df1)) &&
                  all (sapply(all_effdist_values, function(x) x != "1"))) {

                merged_df1[first_col_name] <- dptf_sum + effdist_sum  + merged_df1[[next_effdist_col]]
                merged_df1 <- merged_df1[, !(names(merged_df1) %in% cols_to_remove_noeffdist1)]
                merged_df1[[next_effdist_col]] <- NULL

              }}


            ########################
            cols_to_check <- c(first_cbh_gt5_col,next_effdist_col,Hdist_next_colname,next_effdist_col,next2_effdist_col)
            if (all(cols_to_check %in% colnames(merged_df1))) {


              if (isTRUE(first_suffix_eq1 == first(lad_suffixes) &&
                         first_cbh_gt5_value  > 1.5 &&
                         (Hdist_next_value < first_cbh_gt5_value) &&
                         !Hdist_next2_colname %in% colnames(merged_df1) &&
                         any(sapply(all_effdist_values, function(x) x == "1"))) ) {

                merged_df1[first_col_name] <- dptf_sum + effdist_sum  + merged_df1[[next_effdist_col]] + merged_df1[[next2_effdist_col]]
                merged_df1 <- merged_df1[, !(names(merged_df1) %in% cols_to_remove_noeffdist1)]
                merged_df1[[next_effdist_col]] <- NULL
                merged_df1[[next2_effdist_col]] <- NULL
              }}


            ########################
            cols_to_check <- c(first_cbh_gt5_col,next_effdist_col,Hdist_next_colname,next_effdist_col)
            if (all(cols_to_check %in% colnames(merged_df1))) {


              if (isTRUE(first_suffix_eq1 == first(lad_suffixes) &&
                         first_cbh_gt5_value  > 1.5 &&
                         (Hdist_next_value < first_cbh_gt5_value) &&
                         !Hdist_next2_colname %in% colnames(merged_df1) &&
                         any(sapply(all_effdist_values, function(x) x == "1")) &&
                         length(lad_values_noremove)== 1) ) {

                merged_df1[first_col_name] <- dptf_sum + effdist_sum  + merged_df1[[next_effdist_col]]
                merged_df1 <- merged_df1[, !(names(merged_df1) %in% cols_to_remove_noeffdist1)]
                merged_df1[[next_effdist_col]] <- NULL
              }}


            ########################
            cols_to_check <- c(first_cbh_gt5_col,Hdist_next_colname, dist_col1)
            if (all(cols_to_check %in% colnames(merged_df1))) {


              if (isTRUE(first_suffix_eq1 == first(lad_suffixes) &&
                         first_cbh_gt5_value  > 1.5 &&
                         (Hdist_next_value < first_cbh_gt5_value) &&
                         !Hdist_next2_colname %in% colnames(merged_df1) &&
                         all(sapply(all_effdist_values, function(x) x != "1")) &&
                         length(lad_values_noremove)== 1 && dist1_value == 1) ) {

                merged_df1[first_col_name] <- dist1_value + dptf_sum + effdist_sum
                merged_df1 <- merged_df1[, !(names(merged_df1) %in% cols_to_remove_noeffdist1)]
                merged_df1[[next_effdist_col]] <- NULL
              }}

            ########################

            merged_df1<- get_renamed_df (merged_df1)

          }
        }

        #############################
        #############################

        if (any(suffixes_to_remove > 1)) {

          last_col_name <- paste0("effdist", last_no_consenc_suffix)

          if (last_col_name %in% colnames(merged_df1)) {

            lad_columns2 <- grep("^Hcbh\\d+_Hdptf\\d+$", names(merged_df1),value=TRUE)
            lad_suffixes <- as.numeric(str_extract(lad_columns2, "\\d+$"))

            cols_to_remove1 <- sapply(merged_df1[, lad_columns2], function(col) any(col < 25))
            suffixes_to_remove2 <- sort(as.numeric(stringr::str_extract(names(cols_to_remove1[cols_to_remove1]), "\\d+$")))

            lad_no_remove <- sapply(merged_df1[, lad_columns2], function(col) any(col > 25))
            suffixes_no_remove <- sort(as.numeric(stringr::str_extract(names(lad_no_remove[lad_no_remove]), "\\d+$")))

            colnames_to_extract <- names(lad_no_remove)[lad_no_remove]
            # Extract the values from the desired columns
            lad_values_noremove <- merged_df1[, colnames_to_extract, drop=FALSE]

            next_column_suffix<-  suffixes_to_remove2 +1
            previous_column_suffix<-  suffixes_to_remove2 -1

            effdist_cols <- paste0("effdist", suffixes_to_remove2)

            if (!is.null(effdist_cols) && length(effdist_cols) > 0 && all(effdist_cols %in% colnames(merged_df1))) {

              next_effdist_col <- paste0("effdist", next_column_suffix)
              previous_effdist_col <- paste0("effdist", previous_column_suffix)

            }

            for(suffix in suffixes_to_remove2) {

              pattern2 <- paste0(suffix, "$")
              matching_cols3 <- grep(pattern2, names(merged_df1), value=TRUE)
              cols_to_remove_noeffdist <- matching_cols3[!grepl("^(effdist|treeID)", matching_cols3)]

              # Compute the next and previous suffixes
              next_column_suffix <- suffix + 1
              previous_column_suffix <- suffix - 1

              # Get the values from columns corresponding to the current suffix
              effdist_value <- merged_df1[[paste0("effdist", suffix)]]
              dptf_value <- merged_df1[[paste0("dptf", suffix)]]

              # Use is.null to handle absent columns
              if(is.null(effdist_value)) effdist_value <- 0
              if(is.null(dptf_value)) dptf_value <- 0

              # Compute the sum values for the current suffix
              effdist_sum <- effdist_value
              dptf_sum <- dptf_value


              previous_effdist_col <- paste0("effdist", previous_column_suffix)

              Hdist_previous_colname <- paste0("Hdist", previous_column_suffix)
              Hdptf_previous_colname <- paste0("Hdptf", previous_column_suffix)


              if (any(Hdist_previous_colname %in% colnames(merged_df1))) {
                Hdist_previous_value <- merged_df1[, Hdist_previous_colname]
              } else {
                Hdist_previous_value <- NA  # or some other default value
              }

              if (Hdptf_previous_colname %in% colnames(merged_df1)) {
                Hdptf_previous_value <- merged_df1[, Hdptf_previous_colname]
              } else {
                Hdptf_previous_value <- NA  # or some other default value
              }


              next_effdist_col <- paste0("effdist", next_column_suffix)

              if (next_effdist_col %in% colnames(merged_df1)) {
                next_effdist_value <- merged_df1[, next_effdist_col]
              } else {
                next_effdist_value <- NA  # or some other default value
              }

              Hdist_next_colname <- paste0("Hdist", next_column_suffix)

              if (Hdist_next_colname %in% colnames(merged_df1)) {
                Hdist_next_value <- merged_df1[, Hdist_next_colname]
              } else {
                Hdist_next_value <- NA  # or some other default value
              }

              Hdist_col <- paste0("Hdist", suffix)

              if (Hdist_col %in% colnames(merged_df1)) {
                Hdist_values <- merged_df1[, Hdist_col]
              } else {
                Hdist_values <- NA  # or some other default value
              }

              previous_hdptf_noremove_col <- paste0("Hdptf", first(suffix)-1)

              if (previous_hdptf_noremove_col %in% colnames(merged_df1)) {
                previous_hdptf_noremove_value <- merged_df1[, previous_hdptf_noremove_col]
              } else {
                previous_hdptf_noremove_value <- NA  # or some other default value
              }

              next_hdptf_noremove_col <- paste0("Hdptf", first(suffix) + 1)

              if (next_hdptf_noremove_col %in% colnames(merged_df1)) {
                next_hdptf_noremove_value <- merged_df1[, next_hdptf_noremove_col]
              } else {
                next_hdptf_noremove_value <- NA  # or some other default value
              }


              next_suffix_toremove<-suffix +1
              next_cbh_noremove_col <- paste0("Hcbh", next_suffix_toremove)

              if (next_cbh_noremove_col %in% colnames(merged_df1)) {
                next_cbh_noremove_value <- merged_df1[, next_cbh_noremove_col]
              } else {
                next_cbh_noremove_value <- NA  # or some other default value
              }


              first_Hcbh_col <- grep("^Hcbh\\d+$", names(merged_df1), value = TRUE)[1]

              first_Hcbh_value <- merged_df1[, first_Hcbh_col]

              sufix_col_name <- paste0("effdist", suffix)

              #####################################################

              cols_to_check <- c(previous_hdptf_noremove_col)
              if (all(cols_to_check %in% colnames(merged_df1))) {

                if (isTRUE(suffix != first(lad_suffixes) && (Hdist_values > previous_hdptf_noremove_value)) &&
                    suffix >= last(lad_suffixes) ) {

                  merged_df1 <- merged_df1[, !(names(merged_df1) %in% cols_to_remove_noeffdist)]
                  cols_to_remove <- match(effdist_cols, names(merged_df1))
                  merged_df1 <- merged_df1[,-cols_to_remove]

                }}

              ##############################


              cols_to_check <- c(next_effdist_col, previous_effdist_col, Hdist_previous_colname, Hdist_next_colname, previous_hdptf_noremove_col,     next_cbh_noremove_col,first_Hcbh_col)
              if (all(cols_to_check %in% colnames(merged_df1))) {

                if (isTRUE(suffix != first(lad_suffixes) &&
                           (Hdist_previous_value > previous_hdptf_noremove_value) &&
                           (Hdist_next_value < next_cbh_noremove_value)) &&
                    length (suffixes_no_remove) > 1 &&
                    first_Hcbh_value == 1.5) {


                  merged_df1[sufix_col_name] <- dptf_sum + effdist_sum + merged_df1[[previous_effdist_col]] + merged_df1[[next_effdist_col]]
                  merged_df1[[previous_effdist_col]] <- NULL
                  merged_df1[[next_effdist_col]] <- NULL

                  merged_df1 <- merged_df1[, !(names(merged_df1) %in% cols_to_remove_noeffdist)]


                } else if (previous_effdist_col %in% colnames(merged_df1) && !next_effdist_col %in% colnames(merged_df1)) {
                  merged_df1[sufix_col_name] <- dptf_sum + effdist_sum + merged_df1[[previous_effdist_col]]
                  merged_df1[[previous_effdist_col]] <- NULL

                  merged_df1 <- merged_df1[, !(names(merged_df1) %in% cols_to_remove_noeffdist)]
                }}

              ###############################


              cols_to_check <- c(next_effdist_col, previous_effdist_col, Hdist_col, Hdist_previous_colname, Hdist_next_colname, previous_hdptf_noremove_col, next_cbh_noremove_col,first_Hcbh_col)
              if (all(cols_to_check %in% colnames(merged_df1))) {

                if (isTRUE(suffix != first(lad_suffixes) &&
                           (Hdist_values > previous_hdptf_noremove_value) &&
                           (Hdist_next_value < next_cbh_noremove_value)) &&
                    length (suffixes_no_remove) > 1  &&
                    first_Hcbh_value == 1.5 &&
                    Hdist_previous_value > previous_hdptf_noremove_value) {


                  merged_df1[sufix_col_name] <- dptf_sum + effdist_sum + merged_df1[[previous_effdist_col]] + merged_df1[[next_effdist_col]]
                  merged_df1[[previous_effdist_col]] <- NULL
                  merged_df1[[next_effdist_col]] <- NULL

                  merged_df1 <- merged_df1[, !(names(merged_df1) %in% cols_to_remove_noeffdist)]


                }  else if (previous_effdist_col %in% colnames(merged_df1) && ! next_effdist_col %in% colnames(merged_df1)) {
                  merged_df1[sufix_col_name] <- dptf_sum + effdist_sum + merged_df1[[previous_effdist_col]]
                  merged_df1[[previous_effdist_col]] <- NULL

                  merged_df1 <- merged_df1[, !(names(merged_df1) %in% cols_to_remove_noeffdist)]

                }}


              ###############################

              cols_to_check <- c(next_effdist_col, Hdist_previous_colname, Hdist_next_colname, previous_hdptf_noremove_col,  next_cbh_noremove_col,first_Hcbh_col)
              if (all(cols_to_check %in% colnames(merged_df1))) {

                if (isTRUE(suffix != first(lad_suffixes) &&
                           (Hdist_previous_value < previous_hdptf_noremove_value) &&
                           (Hdist_next_value < next_cbh_noremove_value)) &&
                    length (suffixes_no_remove) >= 2 && first_Hcbh_value > 1.5) {

                  merged_df1[sufix_col_name] <- dptf_sum + effdist_sum + merged_df1[[next_effdist_col]]
                  merged_df1[[next_effdist_col]] <- NULL

                  merged_df1 <- merged_df1[, !(names(merged_df1) %in% cols_to_remove_noeffdist)]
                }}



              ##############################

              cols_to_check <- c(next_effdist_col, previous_effdist_col, Hdist_col, Hdist_previous_colname, Hdist_next_colname, previous_hdptf_noremove_col, next_cbh_noremove_col)
              if (all(cols_to_check %in% colnames(merged_df1))) {

                if (isTRUE(suffix != first(lad_suffixes) &&
                           (Hdist_values > previous_hdptf_noremove_value) &&
                           (Hdist_next_value < next_cbh_noremove_value)) &&
                    length (suffixes_no_remove) > 1  &&
                    first_Hcbh_value == 1.5 &&
                    Hdist_previous_value < previous_hdptf_noremove_value) {


                  merged_df1[sufix_col_name] <- dptf_sum + effdist_sum +  merged_df1[[next_effdist_col]]
                  merged_df1[[next_effdist_col]] <- NULL

                  merged_df1 <- merged_df1[, !(names(merged_df1) %in% cols_to_remove_noeffdist)]


                }  else if (previous_effdist_col %in% colnames(merged_df1) && ! next_effdist_col %in% colnames(merged_df1)) {
                  merged_df1[last_col_name] <- dptf_sum + effdist_sum + merged_df1[[previous_effdist_col]]
                  merged_df1[[previous_effdist_col]] <- NULL

                  merged_df1 <- merged_df1[, !(names(merged_df1) %in% cols_to_remove_noeffdist)]

                }}

              ##############################

              cols_to_check <- c(previous_effdist_col, Hdist_col, Hdist_previous_colname, Hdist_next_colname, previous_hdptf_noremove_col, next_cbh_noremove_col)
              if (all(cols_to_check %in% colnames(merged_df1))) {

                if (isTRUE(suffix != first(lad_suffixes) &&
                           (Hdist_values > previous_hdptf_noremove_value) &&
                           (Hdist_next_value < next_cbh_noremove_value)) &&
                    !next_effdist_col %in% colnames(merged_df1) &&
                    length (suffixes_no_remove) > 1  &&
                    first_Hcbh_value == 1.5) {

                  merged_df1[sufix_col_name] <- dptf_sum + effdist_sum + merged_df1[[previous_effdist_col]]
                  merged_df1[[previous_effdist_col]] <- NULL

                  merged_df1 <- merged_df1[, !(names(merged_df1) %in% cols_to_remove_noeffdist)]

                }}


              ##############################

              cols_to_check <- c( previous_effdist_col, Hdist_previous_colname, Hdist_next_colname, previous_hdptf_noremove_col, next_hdptf_noremove_col)
              if (all(cols_to_check %in% colnames(merged_df1))) {

                if (isTRUE(suffix != first(lad_suffixes) &&
                           (Hdist_previous_value > previous_hdptf_noremove_value) && (
                             Hdist_next_value > next_hdptf_noremove_value))) {

                  merged_df1[sufix_col_name] <- dptf_sum + effdist_sum +  merged_df1[[previous_effdist_col]]
                  merged_df1[[previous_effdist_col]] <- NULL

                  merged_df1 <- merged_df1[, !(names(merged_df1) %in% cols_to_remove_noeffdist)]
                }}


              ##############################

              cols_to_check <- c( previous_effdist_col, Hdist_previous_colname, previous_hdptf_noremove_col)
              if (all(cols_to_check %in% colnames(merged_df1))) {

                if (isTRUE(suffix != first(lad_suffixes) &&
                           (Hdist_previous_value > previous_hdptf_noremove_value) &&
                           (!Hdist_next_colname %in% colnames(merged_df1)))) {

                  merged_df1[sufix_col_name] <- dptf_sum + effdist_sum +  merged_df1[[previous_effdist_col]]
                  merged_df1[[previous_effdist_col]] <- NULL

                  merged_df1 <- merged_df1[, !(names(merged_df1) %in% cols_to_remove_noeffdist)]

                }}


              ##############################


              cols_to_check <- c( previous_effdist_col, Hdist_col, previous_hdptf_noremove_col)
              if (all(cols_to_check %in% colnames(merged_df1))) {

                if (isTRUE(suffix != first(lad_suffixes) &&
                           (Hdist_values > previous_hdptf_noremove_value) &&
                           (!next_effdist_col %in% colnames(merged_df1))&&
                           first_Hcbh_value >= 1.5)) {


                  merged_df1[sufix_col_name] <- dptf_sum + effdist_sum +  merged_df1[[previous_effdist_col]]
                  merged_df1[[previous_effdist_col]] <- NULL

                  merged_df1 <- merged_df1[, !(names(merged_df1) %in% cols_to_remove_noeffdist)]
                }}

              merged_df1<- get_renamed_df (merged_df1)

            }

          } else {

            pattern <- paste0(".*", last_no_consenc_suffix, "$")  # Create a pattern that matches the suffix at the end of a string
            selected_cols <- grep(pattern, names(merged_df1), value=TRUE)
            merged_df1 <- merged_df1[, !(names(merged_df1) %in% selected_cols)]

            second_last_no_consenc_suffix <-last_no_consenc_suffix -1
            effdist_to_remove <- paste0("effdist", second_last_no_consenc_suffix)
            merged_df1 <- merged_df1[, !(names(merged_df1) %in% effdist_to_remove)]

            merged_df1<- get_renamed_df (merged_df1)

          }
        }

        all_effdist_cols <- grep("^effdist\\d+$", names(merged_df1), value = TRUE)

        if (length(all_effdist_cols) > 0) {
          for (col in all_effdist_cols) {
            if (any(merged_df1[, col] == 1)) {
              merged_df1[, col] <- NULL
            }}
        }


        dist_cols <- grep("^dist\\d+$", names(merged_df1), value = TRUE)
        hdist_cols <- grep("^Hdist", names(merged_df1), value = TRUE)

        if (all(dist_cols %in% colnames (merged_df1)) && all(hdist_cols %in% colnames (merged_df1))){
          merged_df1 <- merged_df1 %>%
            dplyr::select(-starts_with("dist"))

        }
        merged_df1<- get_renamed_df (merged_df1)


        if (length(hdist_cols) > length(all_effdist_cols) && merged_df1$Hcbh1 == 1.5) {
          merged_df1 <- merged_df1[, !names(merged_df1) %in% last(hdist_cols)]
        }

        merged_df1<- get_renamed_df (merged_df1)


      }
    }

    #########################################################
    ### NON-CONSECUTIVE SUFFIXES TO REMOVE
    #########################################################

    if (identical(merged_df1,merged_df)) {

      merged_df1 <-merged_df

      lad_columns2 <- grep("^Hcbh\\d+_Hdptf\\d+$", names(merged_df1),value=TRUE)
      lad_suffixes <- as.numeric(str_extract(lad_columns2, "\\d+$"))

      lad_no_remove <- sapply(merged_df1[, lad_columns2], function(col) any(col > 25))
      suffixes_no_remove <- sort(as.numeric(stringr::str_extract(names(lad_no_remove[lad_no_remove]), "\\d+$")))

      colnames_to_extract <- names(lad_no_remove)[lad_no_remove]
      # Extract the values from the desired columns
      lad_values_noremove <- merged_df1[, colnames_to_extract, drop=FALSE]

      pattern <- paste0(suffixes_no_remove, "$")

      # Use sapply to get matches for each pattern
      matching_cols_list <- sapply(pattern, function(pat) grep(pat, names(merged_df1), value=TRUE))

      # Unlist and make unique the result
      matching_cols <- unique(unlist(matching_cols_list))

      # Extract column name
      first_cbh_gt5_col <- first(grep("^Hcbh\\d+$", matching_cols, value = TRUE))


      if (first_cbh_gt5_col %in% colnames(merged_df1)) {
        first_cbh_gt5_value <- merged_df1[, first_cbh_gt5_col]
      } else {
        first_cbh_gt5_value <- NA  # or some other default value
      }

      last_cbh_gt5_col <- last(grep("^Hcbh\\d+$", matching_cols, value = TRUE))

      if (last_cbh_gt5_col %in% colnames(merged_df1)) {
        last_cbh_gt5_value <- merged_df1[, last_cbh_gt5_col]
      } else {
        last_cbh_gt5_value <- NA  # or some other default value
      }


      #########################################
      # Find the columns with value less than 5
      cols_to_remove <- sapply(merged_df1[, lad_columns2], function(col) any(col < 25))
      suffixes_to_remove <- sort(as.numeric(stringr::str_extract(names(cols_to_remove[cols_to_remove]), "\\d+$")))

      pattern2 <- paste0(suffixes_to_remove, "$")
      # Get column names matching the pattern
      matching_cols2 <- unique(unlist(sapply(pattern2, function(pat) {
        grep(pat, names(merged_df1), value=TRUE)
      })))

      cols_to_remove_noeffdist <- matching_cols2[!grepl("^(effdist|treeID)", matching_cols2)]

      # Subset the dataframe using matching columns
      remove_noeffdist_values <- merged_df1[, cols_to_remove_noeffdist]


      first_no_consenc_suffix <- min(suffixes_to_remove)
      last_no_consenc_suffix <- max(suffixes_to_remove)
      next_column_suffix <- last_no_consenc_suffix + 1
      previous_column_suffix <- first_no_consenc_suffix - 1
      first_col_name <- paste0("effdist", first_no_consenc_suffix)
      last_col_name <- paste0("effdist", last_no_consenc_suffix)
      last_lad_suffix<-last(lad_suffixes)


      if (last_col_name %in% colnames(merged_df1)) {
        last_col_value <- merged_df1[, last_col_name]
      } else {
        last_col_value <- NA  # or some other default value
      }

      if(length(suffixes_to_remove) > 0) {


        next_effdist_col <- paste0("effdist", next_column_suffix)
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


        Hdist_next_colname <- paste0("Hdist", next_column_suffix)
        Hdptf_next_colname <- paste0("Hdptf", next_column_suffix)

        if (Hdist_next_colname %in% colnames(merged_df1)) {
          Hdist_next_value <- merged_df1[, Hdist_next_colname]
        } else {
          Hdist_next_value <- NA  # or some other default value
        }

        if (Hdptf_next_colname %in% colnames(merged_df1)) {
          Hdptf_next_value <- merged_df1[, Hdptf_next_colname]
        } else {
          Hdptf_next_value <- NA  # or some other default value
        }

        Hdist_previous_colname <- paste0("Hdist", previous_column_suffix)
        Hdptf_previous_colname <- paste0("Hdptf", previous_column_suffix)

        if (any(Hdist_previous_colname %in% colnames(merged_df1))) {
          Hdist_previous_value <- merged_df1[, Hdist_previous_colname]
        } else {
          Hdist_previous_value <- NA  # or some other default value
        }

        if (Hdptf_previous_colname %in% colnames(merged_df1)) {
          Hdptf_previous_value <- merged_df1[, Hdptf_previous_colname]
        } else {
          Hdptf_previous_value <- NA  # or some other default value
        }

        dist_col1 <- paste0("dist", first_no_consenc_suffix)

        if (dist_col1 %in% colnames(merged_df1)) {
          dist1_value <- merged_df1[, dist_col1]
        } else {
          dist1_value <- NA  # or some other default value
        }

        lad_columns2 <- grep("^Hcbh\\d+_Hdptf\\d+$", names(merged_df1),value=TRUE)
        lad_suffixes <- as.numeric(str_extract(lad_columns2, "\\d+$"))


        ####################################################

        if (any(suffixes_to_remove == 1)) {


          first_suffix_eq1 <- first(suffixes_to_remove[suffixes_to_remove == 1])
          suffixes_to_remove1 <- first_suffix_eq1

          pattern2 <- paste0(suffixes_to_remove1, "$")
          # Get column names matching the pattern
          matching_cols3 <- unique(unlist(sapply(pattern2, function(pat) {
            grep(pat, names(merged_df1), value=TRUE)
          })))
          cols_to_remove_noeffdist1 <- matching_cols3[!grepl("^(effdist|Hdist|treeID)", matching_cols3)]

          next_column_suffix <- suffixes_to_remove1 + 1
          next_col_name <- paste0("effdist", next_column_suffix)

          effdist_cols <- paste0("effdist", suffixes_to_remove1)

          if (!is.null(effdist_cols) && length(effdist_cols) > 0 && effdist_cols %in% colnames(merged_df1)) {

            effdist_values <- sapply(suffixes_to_remove1, function(suf) merged_df1[[paste0("effdist", suf)]])
            dptf_values <- sapply(suffixes_to_remove1, function(suf) merged_df1[[paste0("dptf", suf)]])

            effdist_values <- sapply(effdist_values, function(x) if(is.null(x)) 0 else x)
            dptf_values <- sapply(dptf_values, function(x) if(is.null(x)) 0 else x)

            effdist_sum <- sum(effdist_values)
            dptf_sum <- sum(dptf_values)

            next_effdist_col <- paste0("effdist", next_column_suffix)

            if (next_effdist_col %in% colnames(merged_df1)) {
              next_effdist_value <- merged_df1[, next_effdist_col]
            } else {
              next_effdist_value <- NA  # or some other default value
            }

            next2_effdist_col <- paste0("effdist", next_column_suffix +1)

            if (next2_effdist_col %in% colnames(merged_df1)) {
              next2_effdist_value <- merged_df1[, next2_effdist_col]
            } else {
              next2_effdist_value <- NA  # or some other default value
            }

            Hdist_next_colname <- paste0("Hdist", next_column_suffix)

            if (Hdist_next_colname %in% colnames(merged_df1)) {
              Hdist_next_value <- merged_df1[, Hdist_next_colname]
            } else {
              Hdist_next_value <- NA  # or some other default value
            }

            Hdist_next2_colname <- paste0("Hdist", next_column_suffix +1)

            if (Hdist_next2_colname %in% colnames(merged_df1)) {
              Hdist_next2_value <- merged_df1[, Hdist_next2_colname]
            } else {
              Hdist_next2_value <- NA  # or some other default value
            }


            next_hdptf_noremove_col <- paste0("Hdptf", suffixes_to_remove1 + 1)

            if (next_hdptf_noremove_col %in% colnames(merged_df1)) {
              next_hdptf_noremove_value <- merged_df1[, next_hdptf_noremove_col]
            } else {
              next_hdptf_noremove_value <- NA  # or some other default value
            }

            all_effdist_cols <- grep("^effdist\\d+$", names(merged_df1), value = TRUE)
            all_effdist_values <-merged_df1[,all_effdist_cols]

            Hcbh_cols <- grep("^Hcbh\\d+$", names(merged_df1), value = TRUE)

            if (all(Hcbh_cols %in% colnames(merged_df1))) {
              Hcbh1 <- merged_df1[, first(Hcbh_cols)]
            } else {
              Hcbh1 <- NA  # or some other default value
            }

            lad_columns2 <- grep("^Hcbh\\d+_Hdptf\\d+$", names(merged_df1),value=TRUE)
            lad_suffixes <- as.numeric(str_extract(lad_columns2, "\\d+$"))

            ###############################################################################3


            cols_to_check <- c(first_cbh_gt5_col,dist_col1)
            if (all(cols_to_check %in% colnames(merged_df1))) {

              if (isTRUE(first_suffix_eq1 == first(lad_suffixes)) && first_cbh_gt5_value  > 1.5 &&
                  !next_effdist_col %in% colnames(merged_df1) && dist1_value > 1) {

                merged_df1[first_col_name] <- dptf_sum + effdist_sum
                merged_df1 <- merged_df1[, !(names(merged_df1) %in% cols_to_remove_noeffdist1)]
              }}



            cols_to_check <- c(first_cbh_gt5_col, first(Hcbh_cols),Hdist_next_colname)
            if (all(cols_to_check %in% colnames(merged_df1))) {

              if (isTRUE(first_suffix_eq1 == first(lad_suffixes)) &&
                  first_cbh_gt5_value  > 1.5 && Hcbh1 ==1.5 &&
                  !Hdist_next2_colname %in% colnames(merged_df1) &&
                  Hdist_next_value < first_cbh_gt5_value && all (sapply(all_effdist_values, function(x) x != "1"))) {

                merged_df1[first_col_name] <- dptf_sum + effdist_sum
                merged_df1 <- merged_df1[, !(names(merged_df1) %in% cols_to_remove_noeffdist1)]
              }}


            cols_to_check <- c(first_cbh_gt5_col,Hdist_next_colname)
            if (all(cols_to_check %in% colnames(merged_df1))) {


              if (isTRUE(first_suffix_eq1 == first(lad_suffixes) &&
                         first_cbh_gt5_value  > 1.5 &&
                         (Hdist_next_value > first_cbh_gt5_value) &&
                         !Hdist_next2_colname %in% colnames(merged_df1) &&
                         all(sapply(all_effdist_values, function(x) x != "1")) &&
                         length(lad_values_noremove) > 1) ) {

                merged_df1[first_col_name] <- dptf_sum + effdist_sum
                merged_df1 <- merged_df1[, !(names(merged_df1) %in% cols_to_remove_noeffdist1)]
              }}



            cols_to_check <- c(first_cbh_gt5_col, first(Hcbh_cols),next_effdist_col,Hdist_next2_colname)
            if (all(cols_to_check %in% colnames(merged_df1))) {

              if (isTRUE(first_suffix_eq1 == first(lad_suffixes)) &&
                  first_cbh_gt5_value > 1.5 && Hcbh1== 1.5 &&  ## first distance == 0.5
                  next_effdist_value==1 &&
                  Hdist_next2_value > first_cbh_gt5_value && any(sapply(all_effdist_values, function(x) x == "1"))) {

                merged_df1[first_col_name] <- dptf_sum + effdist_sum
                merged_df1 <- merged_df1[, !(names(merged_df1) %in% cols_to_remove_noeffdist1)]
              }}



            cols_to_check <- c(first_cbh_gt5_col, first(Hcbh_cols),next_effdist_col,Hdist_next2_colname)
            if (all(cols_to_check %in% colnames(merged_df1))) {

              if (isTRUE(first_suffix_eq1 == first(lad_suffixes)) &&
                  first_cbh_gt5_value > 1.5 && Hcbh1== 1.5 &&  ## first distance == 0.5
                  next_effdist_value > 1 &&
                  Hdist_next2_value > first_cbh_gt5_value && all (sapply(all_effdist_values, function(x) x != "1"))) {

                merged_df1[first_col_name] <- dptf_sum + effdist_sum
                merged_df1 <- merged_df1[, !(names(merged_df1) %in% cols_to_remove_noeffdist1)]
              }}




            cols_to_check <- c(first_cbh_gt5_col, first(Hcbh_cols))
            if (all(cols_to_check %in% colnames(merged_df1))) {

              if (isTRUE(first_suffix_eq1 == first(lad_suffixes)) &&
                  first_cbh_gt5_value > 1.5 && Hcbh1== 1.5 &&  ## first distance == 0.5
                  !next_effdist_col %in% colnames(merged_df1) &&
                  !Hdist_next2_colname %in% colnames(merged_df1) && all (sapply(all_effdist_values, function(x) x != "1"))) {

                merged_df1[first_col_name] <- dptf_sum + effdist_sum
                merged_df1 <- merged_df1[, !(names(merged_df1) %in% cols_to_remove_noeffdist1)]
              }}

            #################################

            cols_to_check <- c(first_cbh_gt5_col, first(Hcbh_cols),next_effdist_col,Hdist_next2_colname,next_hdptf_noremove_col)
            if (all(cols_to_check %in% colnames(merged_df1))) {

              if (isTRUE(first_suffix_eq1 == first(lad_suffixes)) &&
                  first_cbh_gt5_value  > 1.5 && Hcbh1 ==1.5 &&  ## first distance == 0.5
                  Hdist_next2_value < next_hdptf_noremove_value &&
                  length(lad_values_noremove) == 1) {

                merged_df1[first_col_name] <- dptf_sum + effdist_sum + merged_df1[[next_effdist_col]]
                merged_df1 <- merged_df1[, !(names(merged_df1) %in% cols_to_remove_noeffdist1)]
                merged_df1[[next_effdist_col]] <- NULL
              }}

            #################################

            cols_to_check <- c(first_cbh_gt5_col, first(Hcbh_cols),next_effdist_col,Hdist_next_colname, Hdist_next2_colname,next_hdptf_noremove_col)
            if (all(cols_to_check %in% colnames(merged_df1))) {

              if (isTRUE(first_suffix_eq1 == first(lad_suffixes) &&
                         first_cbh_gt5_value  > 1.5 && Hcbh1 ==1.5 &&  ## first distance == 0.5
                         (Hdist_next_value < first_cbh_gt5_value) &&
                         Hdist_next2_value > next_hdptf_noremove_value &&
                         any(sapply(all_effdist_values, function(x) x == "1")) &&
                         length(lad_values_noremove) > 1)) {

                merged_df1[first_col_name] <- dptf_sum + effdist_sum + merged_df1[[next_effdist_col]]
                merged_df1 <- merged_df1[, !(names(merged_df1) %in% cols_to_remove_noeffdist1)]
                merged_df1[[next_effdist_col]] <- NULL
              }}

            #################################

            cols_to_check <- c(first_cbh_gt5_col,next_effdist_col,Hdist_next_colname,Hdist_next2_colname,next_hdptf_noremove_col)
            if (all(cols_to_check %in% colnames(merged_df1))) {


              if (isTRUE(first_suffix_eq1 == first(lad_suffixes) &&
                         first_cbh_gt5_value  > 1.5 && Hcbh1 > 1.5 &&
                         (Hdist_next_value < first_cbh_gt5_value) &&
                         Hdist_next2_value > next_hdptf_noremove_value )) {

                merged_df1[first_col_name] <- dptf_sum + effdist_sum  + merged_df1[[next_effdist_col]]
                merged_df1 <- merged_df1[, !(names(merged_df1) %in% cols_to_remove_noeffdist1)]
                merged_df1[[next_effdist_col]] <- NULL
              }}


            ########################

            cols_to_check <- c(first_cbh_gt5_col,next_effdist_col,Hdist_next_colname,Hdist_next2_colname,next_hdptf_noremove_col)
            if (all(cols_to_check %in% colnames(merged_df1))) {


              if (isTRUE(first_suffix_eq1 == first(lad_suffixes) &&
                         first_cbh_gt5_value  > 1.5 &&
                         (Hdist_next_value < first_cbh_gt5_value) &&
                         Hdist_next2_value < next_hdptf_noremove_value )) {

                merged_df1[first_col_name] <- dptf_sum + effdist_sum  + merged_df1[[next_effdist_col]]
                merged_df1 <- merged_df1[, !(names(merged_df1) %in% cols_to_remove_noeffdist1)]
                merged_df1[[next_effdist_col]] <- NULL
              }}

            ########################

            cols_to_check <- c(first_cbh_gt5_col,next_effdist_col,Hdist_next_colname,next_effdist_col,Hcbh_cols)
            if (all(cols_to_check %in% colnames(merged_df1))) {


              if (isTRUE(first_suffix_eq1 == first(lad_suffixes) &&
                         first_cbh_gt5_value  > 1.5 && Hcbh1 == 1.5 &&
                         (Hdist_next_value < first_cbh_gt5_value) &&
                         !Hdist_next2_colname %in% colnames(merged_df1) &&
                         length(lad_values_noremove) == 1)) {

                merged_df1[first_col_name] <- dptf_sum + effdist_sum  + merged_df1[[next_effdist_col]]
                merged_df1 <- merged_df1[, !(names(merged_df1) %in% cols_to_remove_noeffdist1)]
                merged_df1[[next_effdist_col]] <- NULL

              }}


            ########################

            cols_to_check <- c(first_cbh_gt5_col,next_effdist_col,Hdist_next_colname,all_effdist_cols)
            if (all(cols_to_check %in% colnames(merged_df1))) {


              if (isTRUE(first_suffix_eq1 == first(lad_suffixes)) &&
                  first_cbh_gt5_value  > 1.5 &&
                  (Hdist_next_value < first_cbh_gt5_value) &&
                  (!Hdist_next2_colname %in% colnames(merged_df1)) &&
                  all (sapply(all_effdist_values, function(x) x != "1"))) {

                merged_df1[first_col_name] <- dptf_sum + effdist_sum  + merged_df1[[next_effdist_col]]
                merged_df1 <- merged_df1[, !(names(merged_df1) %in% cols_to_remove_noeffdist1)]
                merged_df1[[next_effdist_col]] <- NULL

              }}


            ########################
            cols_to_check <- c(first_cbh_gt5_col,next_effdist_col,Hdist_next_colname,next_effdist_col,next2_effdist_col)
            if (all(cols_to_check %in% colnames(merged_df1))) {


              if (isTRUE(first_suffix_eq1 == first(lad_suffixes) &&
                         first_cbh_gt5_value  > 1.5 &&
                         (Hdist_next_value < first_cbh_gt5_value) &&
                         !Hdist_next2_colname %in% colnames(merged_df1) &&
                         any(sapply(all_effdist_values, function(x) x == "1"))) ) {

                merged_df1[first_col_name] <- dptf_sum + effdist_sum  + merged_df1[[next_effdist_col]] + merged_df1[[next2_effdist_col]]
                merged_df1 <- merged_df1[, !(names(merged_df1) %in% cols_to_remove_noeffdist1)]
                merged_df1[[next_effdist_col]] <- NULL
                merged_df1[[next2_effdist_col]] <- NULL
              }}


            ########################
            cols_to_check <- c(first_cbh_gt5_col,next_effdist_col,Hdist_next_colname,next_effdist_col)
            if (all(cols_to_check %in% colnames(merged_df1))) {


              if (isTRUE(first_suffix_eq1 == first(lad_suffixes) &&
                         first_cbh_gt5_value  > 1.5 &&
                         (Hdist_next_value < first_cbh_gt5_value) &&
                         !Hdist_next2_colname %in% colnames(merged_df1) &&
                         any(sapply(all_effdist_values, function(x) x == "1")) &&
                         length(lad_values_noremove)== 1) ) {

                merged_df1[first_col_name] <- dptf_sum + effdist_sum  + merged_df1[[next_effdist_col]]
                merged_df1 <- merged_df1[, !(names(merged_df1) %in% cols_to_remove_noeffdist1)]
                merged_df1[[next_effdist_col]] <- NULL
              }}


            ########################
            cols_to_check <- c(first_cbh_gt5_col,Hdist_next_colname, dist_col1)
            if (all(cols_to_check %in% colnames(merged_df1))) {


              if (isTRUE(first_suffix_eq1 == first(lad_suffixes) &&
                         first_cbh_gt5_value  > 1.5 &&
                         (Hdist_next_value < first_cbh_gt5_value) &&
                         !Hdist_next2_colname %in% colnames(merged_df1) &&
                         all(sapply(all_effdist_values, function(x) x != "1")) &&
                         length(lad_values_noremove)== 1 && dist1_value == 1) ) {

                merged_df1[first_col_name] <- dist1_value + dptf_sum + effdist_sum
                merged_df1 <- merged_df1[, !(names(merged_df1) %in% cols_to_remove_noeffdist1)]
                merged_df1[[next_effdist_col]] <- NULL
              }}

            ########################

            merged_df1<- get_renamed_df (merged_df1)

          }
        }

        #############################
        #############################

        if (any(suffixes_to_remove > 1)) {

          last_col_name <- paste0("effdist", last_no_consenc_suffix)

          if (last_col_name %in% colnames(merged_df1)) {

            lad_columns2 <- grep("^Hcbh\\d+_Hdptf\\d+$", names(merged_df1),value=TRUE)
            lad_suffixes <- as.numeric(str_extract(lad_columns2, "\\d+$"))

            cols_to_remove1 <- sapply(merged_df1[, lad_columns2], function(col) any(col < 25))
            suffixes_to_remove2 <- sort(as.numeric(stringr::str_extract(names(cols_to_remove1[cols_to_remove1]), "\\d+$")))

            lad_no_remove <- sapply(merged_df1[, lad_columns2], function(col) any(col > 25))
            suffixes_no_remove <- sort(as.numeric(stringr::str_extract(names(lad_no_remove[lad_no_remove]), "\\d+$")))

            colnames_to_extract <- names(lad_no_remove)[lad_no_remove]
            # Extract the values from the desired columns
            lad_values_noremove <- merged_df1[, colnames_to_extract, drop=FALSE]

            next_column_suffix<-  suffixes_to_remove2 +1
            previous_column_suffix<-  suffixes_to_remove2 -1

            effdist_cols <- paste0("effdist", suffixes_to_remove2)

            if (!is.null(effdist_cols) && length(effdist_cols) > 0 && all(effdist_cols %in% colnames(merged_df1))) {

              next_effdist_col <- paste0("effdist", next_column_suffix)
              previous_effdist_col <- paste0("effdist", previous_column_suffix)

            }

            for(suffix in suffixes_to_remove2) {

              pattern2 <- paste0(suffix, "$")
              matching_cols3 <- grep(pattern2, names(merged_df1), value=TRUE)
              cols_to_remove_noeffdist <- matching_cols3[!grepl("^(effdist|treeID)", matching_cols3)]

              # Compute the next and previous suffixes
              next_column_suffix <- suffix + 1
              previous_column_suffix <- suffix - 1

              # Get the values from columns corresponding to the current suffix
              effdist_value <- merged_df1[[paste0("effdist", suffix)]]
              dptf_value <- merged_df1[[paste0("dptf", suffix)]]

              # Use is.null to handle absent columns
              if(is.null(effdist_value)) effdist_value <- 0
              if(is.null(dptf_value)) dptf_value <- 0

              # Compute the sum values for the current suffix
              effdist_sum <- effdist_value
              dptf_sum <- dptf_value


              previous_effdist_col <- paste0("effdist", previous_column_suffix)

              Hdist_previous_colname <- paste0("Hdist", previous_column_suffix)
              Hdptf_previous_colname <- paste0("Hdptf", previous_column_suffix)


              if (any(Hdist_previous_colname %in% colnames(merged_df1))) {
                Hdist_previous_value <- merged_df1[, Hdist_previous_colname]
              } else {
                Hdist_previous_value <- NA  # or some other default value
              }

              if (Hdptf_previous_colname %in% colnames(merged_df1)) {
                Hdptf_previous_value <- merged_df1[, Hdptf_previous_colname]
              } else {
                Hdptf_previous_value <- NA  # or some other default value
              }


              next_effdist_col <- paste0("effdist", next_column_suffix)

              if (next_effdist_col %in% colnames(merged_df1)) {
                next_effdist_value <- merged_df1[, next_effdist_col]
              } else {
                next_effdist_value <- NA  # or some other default value
              }

              Hdist_next_colname <- paste0("Hdist", next_column_suffix)

              if (Hdist_next_colname %in% colnames(merged_df1)) {
                Hdist_next_value <- merged_df1[, Hdist_next_colname]
              } else {
                Hdist_next_value <- NA  # or some other default value
              }

              Hdist_col <- paste0("Hdist", suffix)

              if (Hdist_col %in% colnames(merged_df1)) {
                Hdist_values <- merged_df1[, Hdist_col]
              } else {
                Hdist_values <- NA  # or some other default value
              }

              previous_hdptf_noremove_col <- paste0("Hdptf", first(suffix)-1)

              if (previous_hdptf_noremove_col %in% colnames(merged_df1)) {
                previous_hdptf_noremove_value <- merged_df1[, previous_hdptf_noremove_col]
              } else {
                previous_hdptf_noremove_value <- NA  # or some other default value
              }

              next_hdptf_noremove_col <- paste0("Hdptf", first(suffix) + 1)

              if (next_hdptf_noremove_col %in% colnames(merged_df1)) {
                next_hdptf_noremove_value <- merged_df1[, next_hdptf_noremove_col]
              } else {
                next_hdptf_noremove_value <- NA  # or some other default value
              }


              next_suffix_toremove<-suffix +1
              next_cbh_noremove_col <- paste0("Hcbh", next_suffix_toremove)

              if (next_cbh_noremove_col %in% colnames(merged_df1)) {
                next_cbh_noremove_value <- merged_df1[, next_cbh_noremove_col]
              } else {
                next_cbh_noremove_value <- NA  # or some other default value
              }


              first_Hcbh_col <- grep("^Hcbh\\d+$", names(merged_df1), value = TRUE)[1]

              first_Hcbh_value <- merged_df1[, first_Hcbh_col]

              sufix_col_name <- paste0("effdist", suffix)

              #####################################################

              cols_to_check <- c(previous_hdptf_noremove_col)
              if (all(cols_to_check %in% colnames(merged_df1))) {

                if (isTRUE(suffix != first(lad_suffixes) && (Hdist_values > previous_hdptf_noremove_value)) &&
                    suffix >= last(lad_suffixes) ) {

                  merged_df1 <- merged_df1[, !(names(merged_df1) %in% cols_to_remove_noeffdist)]
                  cols_to_remove <- match(effdist_cols, names(merged_df1))
                  merged_df1 <- merged_df1[,-cols_to_remove]

                }}

              ##############################


              cols_to_check <- c(next_effdist_col, previous_effdist_col, Hdist_previous_colname, Hdist_next_colname, previous_hdptf_noremove_col,     next_cbh_noremove_col,first_Hcbh_col)
              if (all(cols_to_check %in% colnames(merged_df1))) {

                if (isTRUE(suffix != first(lad_suffixes) &&
                           (Hdist_previous_value > previous_hdptf_noremove_value) &&
                           (Hdist_next_value < next_cbh_noremove_value)) &&
                    length (suffixes_no_remove) > 1 &&
                    first_Hcbh_value == 1.5) {


                  merged_df1[sufix_col_name] <- dptf_sum + effdist_sum + merged_df1[[previous_effdist_col]] + merged_df1[[next_effdist_col]]
                  merged_df1[[previous_effdist_col]] <- NULL
                  merged_df1[[next_effdist_col]] <- NULL

                  merged_df1 <- merged_df1[, !(names(merged_df1) %in% cols_to_remove_noeffdist)]


                } else if (previous_effdist_col %in% colnames(merged_df1) && !next_effdist_col %in% colnames(merged_df1)) {
                  merged_df1[sufix_col_name] <- dptf_sum + effdist_sum + merged_df1[[previous_effdist_col]]
                  merged_df1[[previous_effdist_col]] <- NULL

                  merged_df1 <- merged_df1[, !(names(merged_df1) %in% cols_to_remove_noeffdist)]
                }}

              ###############################


              cols_to_check <- c(next_effdist_col, previous_effdist_col, Hdist_col, Hdist_previous_colname, Hdist_next_colname, previous_hdptf_noremove_col, next_cbh_noremove_col,first_Hcbh_col)
              if (all(cols_to_check %in% colnames(merged_df1))) {

                if (isTRUE(suffix != first(lad_suffixes) &&
                           (Hdist_values > previous_hdptf_noremove_value) &&
                           (Hdist_next_value < next_cbh_noremove_value)) &&
                    length (suffixes_no_remove) > 1  &&
                    first_Hcbh_value == 1.5 &&
                    Hdist_previous_value > previous_hdptf_noremove_value) {


                  merged_df1[sufix_col_name] <- dptf_sum + effdist_sum + merged_df1[[previous_effdist_col]] + merged_df1[[next_effdist_col]]
                  merged_df1[[previous_effdist_col]] <- NULL
                  merged_df1[[next_effdist_col]] <- NULL

                  merged_df1 <- merged_df1[, !(names(merged_df1) %in% cols_to_remove_noeffdist)]


                }  else if (previous_effdist_col %in% colnames(merged_df1) && ! next_effdist_col %in% colnames(merged_df1)) {
                  merged_df1[sufix_col_name] <- dptf_sum + effdist_sum + merged_df1[[previous_effdist_col]]
                  merged_df1[[previous_effdist_col]] <- NULL

                  merged_df1 <- merged_df1[, !(names(merged_df1) %in% cols_to_remove_noeffdist)]

                }}


              ###############################

              cols_to_check <- c(next_effdist_col, Hdist_previous_colname, Hdist_next_colname, previous_hdptf_noremove_col,  next_cbh_noremove_col,first_Hcbh_col)
              if (all(cols_to_check %in% colnames(merged_df1))) {

                if (isTRUE(suffix != first(lad_suffixes) &&
                           (Hdist_previous_value < previous_hdptf_noremove_value) &&
                           (Hdist_next_value < next_cbh_noremove_value)) &&
                    length (suffixes_no_remove) >= 2 && first_Hcbh_value > 1.5) {

                  merged_df1[sufix_col_name] <- dptf_sum + effdist_sum + merged_df1[[next_effdist_col]]
                  merged_df1[[next_effdist_col]] <- NULL

                  merged_df1 <- merged_df1[, !(names(merged_df1) %in% cols_to_remove_noeffdist)]
                }}



              ##############################

              cols_to_check <- c(next_effdist_col, previous_effdist_col, Hdist_col, Hdist_previous_colname, Hdist_next_colname, previous_hdptf_noremove_col, next_cbh_noremove_col)
              if (all(cols_to_check %in% colnames(merged_df1))) {

                if (isTRUE(suffix != first(lad_suffixes) &&
                           (Hdist_values > previous_hdptf_noremove_value) &&
                           (Hdist_next_value < next_cbh_noremove_value)) &&
                    length (suffixes_no_remove) > 1  &&
                    first_Hcbh_value == 1.5 &&
                    Hdist_previous_value < previous_hdptf_noremove_value) {


                  merged_df1[sufix_col_name] <- dptf_sum + effdist_sum +  merged_df1[[next_effdist_col]]
                  merged_df1[[next_effdist_col]] <- NULL

                  merged_df1 <- merged_df1[, !(names(merged_df1) %in% cols_to_remove_noeffdist)]


                }  else if (previous_effdist_col %in% colnames(merged_df1) && ! next_effdist_col %in% colnames(merged_df1)) {
                  merged_df1[last_col_name] <- dptf_sum + effdist_sum + merged_df1[[previous_effdist_col]]
                  merged_df1[[previous_effdist_col]] <- NULL

                  merged_df1 <- merged_df1[, !(names(merged_df1) %in% cols_to_remove_noeffdist)]

                }}

              ##############################

              cols_to_check <- c(previous_effdist_col, Hdist_col, Hdist_previous_colname, Hdist_next_colname, previous_hdptf_noremove_col, next_cbh_noremove_col)
              if (all(cols_to_check %in% colnames(merged_df1))) {

                if (isTRUE(suffix != first(lad_suffixes) &&
                           (Hdist_values > previous_hdptf_noremove_value) &&
                           (Hdist_next_value < next_cbh_noremove_value)) &&
                    !next_effdist_col %in% colnames(merged_df1) &&
                    length (suffixes_no_remove) > 1  &&
                    first_Hcbh_value == 1.5) {

                  merged_df1[sufix_col_name] <- dptf_sum + effdist_sum + merged_df1[[previous_effdist_col]]
                  merged_df1[[previous_effdist_col]] <- NULL

                  merged_df1 <- merged_df1[, !(names(merged_df1) %in% cols_to_remove_noeffdist)]

                }}


              ##############################

              cols_to_check <- c( previous_effdist_col, Hdist_previous_colname, Hdist_next_colname, previous_hdptf_noremove_col, next_hdptf_noremove_col)
              if (all(cols_to_check %in% colnames(merged_df1))) {

                if (isTRUE(suffix != first(lad_suffixes) &&
                           (Hdist_previous_value > previous_hdptf_noremove_value) && (
                             Hdist_next_value > next_hdptf_noremove_value))) {

                  merged_df1[sufix_col_name] <- dptf_sum + effdist_sum +  merged_df1[[previous_effdist_col]]
                  merged_df1[[previous_effdist_col]] <- NULL

                  merged_df1 <- merged_df1[, !(names(merged_df1) %in% cols_to_remove_noeffdist)]
                }}


              ##############################

              cols_to_check <- c( previous_effdist_col, Hdist_previous_colname, previous_hdptf_noremove_col)
              if (all(cols_to_check %in% colnames(merged_df1))) {

                if (isTRUE(suffix != first(lad_suffixes) &&
                           (Hdist_previous_value > previous_hdptf_noremove_value) &&
                           (!Hdist_next_colname %in% colnames(merged_df1)))) {

                  merged_df1[sufix_col_name] <- dptf_sum + effdist_sum +  merged_df1[[previous_effdist_col]]
                  merged_df1[[previous_effdist_col]] <- NULL

                  merged_df1 <- merged_df1[, !(names(merged_df1) %in% cols_to_remove_noeffdist)]

                }}


              ##############################


              cols_to_check <- c( previous_effdist_col, Hdist_col, previous_hdptf_noremove_col)
              if (all(cols_to_check %in% colnames(merged_df1))) {

                if (isTRUE(suffix != first(lad_suffixes) &&
                           (Hdist_values > previous_hdptf_noremove_value) &&
                           (!next_effdist_col %in% colnames(merged_df1))&&
                           first_Hcbh_value >= 1.5)) {


                  merged_df1[sufix_col_name] <- dptf_sum + effdist_sum +  merged_df1[[previous_effdist_col]]
                  merged_df1[[previous_effdist_col]] <- NULL

                  merged_df1 <- merged_df1[, !(names(merged_df1) %in% cols_to_remove_noeffdist)]
                }}

              merged_df1<- get_renamed_df (merged_df1)

            }

          } else {

            pattern <- paste0(".*", last_no_consenc_suffix, "$")  # Create a pattern that matches the suffix at the end of a string
            selected_cols <- grep(pattern, names(merged_df1), value=TRUE)
            merged_df1 <- merged_df1[, !(names(merged_df1) %in% selected_cols)]

            second_last_no_consenc_suffix <-last_no_consenc_suffix -1
            effdist_to_remove <- paste0("effdist", second_last_no_consenc_suffix)
            merged_df1 <- merged_df1[, !(names(merged_df1) %in% effdist_to_remove)]

            merged_df1<- get_renamed_df (merged_df1)

          }
        }

        all_effdist_cols <- grep("^effdist\\d+$", names(merged_df1), value = TRUE)

        if (length(all_effdist_cols) > 0) {
          for (col in all_effdist_cols) {
            if (any(merged_df1[, col] == 1)) {
              merged_df1[, col] <- NULL
            }}
        }


        dist_cols <- grep("^dist\\d+$", names(merged_df1), value = TRUE)
        hdist_cols <- grep("^Hdist", names(merged_df1), value = TRUE)

        if (all(dist_cols %in% colnames (merged_df1)) && all(hdist_cols %in% colnames (merged_df1))){
          merged_df1 <- merged_df1 %>%
            dplyr::select(-starts_with("dist"))

        }
        merged_df1<- get_renamed_df (merged_df1)


        if (length(hdist_cols) > length(all_effdist_cols) && merged_df1$Hcbh1 == 1.5) {
          merged_df1 <- merged_df1[, !names(merged_df1) %in% last(hdist_cols)]
        }

        merged_df1<- get_renamed_df (merged_df1)



      }
    }
  }

  ##########################################################################
  # IF NOT ANY CHANGE
  ##########################################################################

  if (identical(merged_df, merged_df1)) {

    all_effdist_cols <- grep("^effdist\\d+$", names(merged_df1), value = TRUE)

    if (length(all_effdist_cols) > 0) {
      for (col in all_effdist_cols) {
        if (any(merged_df1[, col] == 1)) {
          merged_df1[, col] <- NULL
        }}
    }

    dist_cols <- grep("^dist\\d+$", names(merged_df1), value = TRUE)
    hdist_cols <- grep("^Hdist", names(merged_df1), value = TRUE)

    if (all(dist_cols %in% colnames (merged_df1)) && all(hdist_cols %in% colnames (merged_df1))){
      merged_df1 <- merged_df1 %>%
        dplyr::select(-starts_with("dist"))

    }

    merged_df1<- get_renamed_df (merged_df1)

  }


  ####################################

  if(!all(cols_to_append  %in% colnames(merged_df1))){
    merged_df1<- cbind( merged_df1, append_df)
  }


  if(!all(cols_to_append  %in% colnames(merged_df))){
    merged_df<- cbind( merged_df, append_df)
  }


  all_effdist_cols <- grep("^effdist\\d+$", names(merged_df1), value = TRUE)
  dist_cols <- grep("^dist\\d+$", names(merged_df1), value = TRUE)
  hdist_cols <- grep("^Hdist\\d+$", names(merged_df1), value = TRUE)

  if (length(all_effdist_cols) > 0 && length(hdist_cols) > 0) {

    for (col in all_effdist_cols) {
      if (any(merged_df1[, col] == 1)) {
        # Extract the numeric suffix
        suffix <- sub("^effdist", "", col)

        # Identify corresponding Hdist and dist columns by the suffix
        corresponding_hdist <- paste0("Hdist", suffix)
        corresponding_dist <- paste0("dist", suffix)

        # Set them to NULL
        merged_df1[, col] <- NULL
        if(corresponding_hdist %in% names(merged_df1)) {
          merged_df1[, corresponding_hdist] <- NULL
        }
        if(corresponding_dist %in% names(merged_df1)) {
          merged_df1[, corresponding_dist] <- NULL
        }
      }
    }


    merged_df1<- get_renamed_df (merged_df1)

    ###########################################################
    dist_cols <- grep("^dist\\d+$", names(merged_df1), value = TRUE)

    if (length(dist_cols) > 0) {
      for (col in dist_cols) {
        if (any(merged_df1[, col] == 1)) {
          merged_df1[, col] <- NULL
        }}
    }

    merged_df1<- get_renamed_df (merged_df1)

    ##########################################################################

    # For each column starting with "Hdist"
    hcbh_cols <- grep("^Hcbh\\d+$", colnames(merged_df1), value = TRUE)
    hdist_cols <- grep("^Hdist\\d+$", colnames(merged_df1), value = TRUE)

    for (j in seq_along(hdist_cols)) {
      if (j < length(hcbh_cols) && merged_df1$Hcbh1 == 1.5) { # Ensure we're not exceeding the bounds of hcbh columns
        hcbh_col <- hcbh_cols[j + 1] # Accessing the next Hcbh column
        hdist_col <- hdist_cols[j]

        # Adjust value if the condition is met
        if (merged_df1[1, hdist_col] != merged_df1[1, hcbh_col] - 1) {
          merged_df1[1, hdist_col] <- merged_df1[1, hcbh_col] - 1
        }
      }
    }

    # Iterate through all the Hdist and Hcbh columns
    for (j in seq_along(hdist_cols)) {
      # Ensure we're not exceeding the bounds of hcbh columns
      if (j <= length(hcbh_cols)) {
        hcbh_col <- hcbh_cols[j]
        hdist_col <- hdist_cols[j]

        # Check the condition for Hcbh1 (this assumes that the condition for Hcbh1 holds for all other Hcbh columns)
        if (merged_df1$Hcbh1 > 1.5) {
          # Set Hdist value to Hcbh - 1
          merged_df1[1, hdist_col] <- merged_df1[1, hcbh_col] - 1
        }
      }
    }

    ########################################
    hdist_cols <- grep("^Hdist\\d+$", colnames(merged_df1), value = TRUE)

    # Only proceed if there's more than one Hdist column
    if (length(hdist_cols) > 1) {

      cols_to_remove <- c()  # Initialize an empty vector

      # Loop through Hdist columns starting from the second one
      for (j in 2:length(hdist_cols)) {
        current_value <- merged_df1[1, hdist_cols[j]]
        previous_value <- merged_df1[1, hdist_cols[j-1]]

        # Check if neither of the values is NA and if they are equal
        if (!is.na(current_value) && !is.na(previous_value) && current_value == previous_value) {
          cols_to_remove <- c(cols_to_remove, hdist_cols[j])
        }
      }

      # Remove the columns outside the loop
      merged_df1 <- merged_df1[, !names(merged_df1) %in% cols_to_remove, drop = FALSE]
    }


    merged_df1<- get_renamed_df (merged_df1)

  }

  ######## Hcbh with max % LAD    ###########################

  # Check if df_sub is not empty and contains "Hdist" in column names
  col_names <- names(merged_df1)
  hcbh_cols <- grep("^Hcbh\\d+$", col_names, value = TRUE)
  hdptf_cols <- grep("^Hdptf\\d+$", col_names, value = TRUE)
  hdist_cols <- grep("^Hdist\\d+$", col_names, value = TRUE)
  effdist_cols <- grep("^effdist\\d+$", col_names, value = TRUE)


  # Create a vector of all Hdist values
  hcbh_values <- unlist(merged_df1[1, hcbh_cols])
  hdist_values <- unlist(merged_df1[1, hdist_cols])
  hdepth_values <- unlist(merged_df1[1, hdptf_cols])
  effdist_values <- unlist(merged_df1[1, effdist_cols])

  # Identify lad columns based on naming pattern
  lad_columns <- grep("^Hcbh\\d+_Hdptf\\d+$", names(merged_df1), value = TRUE)

  # Flatten all lad values into a vector
  lad_values <- as.vector(unlist(merged_df1[, lad_columns]))

  # Get the maximum value across all lad columns
  max_lad_value <- max(lad_values, na.rm = TRUE)

  # Get the index (column number) of the last occurrence of the max value
  index_max_lad <- max(which(lad_values == max_lad_value))

  # Now you can use index_max_lad to get the corresponding column
  max_lad_column <- lad_columns[index_max_lad]

  # Get the max lad value from the first row of the dataframe
  max_lad_values <- unlist(merged_df1[1, max_lad_column])

  # Extract suffix from max_lad_column
  suffix <- as.numeric(sub(".*Hcbh(\\d+)_.*", "\\1", max_lad_column))
  suffix1<-suffix-1

  # Use the suffix to get the corresponding Hcbh, Hdist and Hdepth columns
  max_Hcbh_column <- paste0("Hcbh", suffix)
  max_Hdist_column <- paste0("Hdist", suffix)
  max_Hdist_column_minusone <- paste0("Hdist", suffix1)
  max_effdist_column_minusone <- paste0("effdist", suffix1)
  max_Hdepth_column <- paste0("Hdptf", suffix)
  max_effdist_column <- paste0("effdist", suffix)
  max_dptf_column <- paste0("dptf", suffix)

  if (max_Hdist_column %in% colnames(merged_df1) && max_effdist_column %in% colnames(merged_df1)) {
    max_df <- data.frame(
      maxlad_Hcbh = merged_df1[[max_Hcbh_column]],
      maxlad_Hdist = merged_df1[[max_Hdist_column]],
      maxlad_Hdptf = merged_df1[[max_Hdepth_column]],
      maxlad_dptf = merged_df1[[max_dptf_column]],
      maxlad_effdist = merged_df1[[max_effdist_column]]
    )
    names(max_df) <- c("maxlad_Hcbh", "maxlad_Hdist", "maxlad_Hdptf","maxlad_dptf", "maxlad_effdist")

  }


  if (max_Hdist_column %in% colnames(merged_df1) && max_effdist_column %in% colnames(merged_df1) && merged_df1$Hcbh1 == 1.5) {
    max_df <- data.frame(
      maxlad_Hcbh = merged_df1[[max_Hcbh_column]],
      maxlad_Hdist = merged_df1[[max_Hdist_column]],
      maxlad_Hdptf = merged_df1[[max_Hdepth_column]],
      maxlad_dptf = merged_df1[[max_dptf_column]],
      maxlad_effdist = merged_df1[[max_effdist_column]]
    )
    names(max_df) <- c("maxlad_Hcbh", "maxlad_Hdist", "maxlad_Hdptf","maxlad_dptf", "maxlad_effdist")

  }

  if (!max_Hdist_column %in% colnames(merged_df1) && all(max_Hdist_column_minusone %in% colnames(merged_df1)) && all(max_effdist_column_minusone %in% colnames(merged_df1))) {
    # Check if the columns exist; if not, use default values
    max_Hdist_val <- ifelse(max_Hdist_column_minusone %in% colnames(merged_df1), merged_df1[[max_Hdist_column_minusone]], NA)
    max_effdist_val <- ifelse(max_effdist_column_minusone %in% colnames(merged_df1), merged_df1[[max_effdist_column_minusone]], NA)

    max_df <- data.frame(
      maxlad_Hcbh = merged_df1[[max_Hcbh_column]],
      maxlad_Hdist = merged_df1[[max_Hdist_column_minusone]],
      maxlad_Hdptf = merged_df1[[max_Hdepth_column]],
      maxlad_dptf = merged_df1[[max_dptf_column]],
      maxlad_effdist = merged_df1[[max_effdist_column_minusone]]
    )
    names(max_df) <- c("maxlad_Hcbh", "maxlad_Hdist", "maxlad_Hdptf", "maxlad_dptf","maxlad_effdist")

  }

  if (!max_Hdist_column %in% colnames(merged_df1) && !max_effdist_column %in% colnames(merged_df1) && !exists ("max_df")) {
    max_df <- data.frame(
      maxlad_Hcbh = merged_df1[[max_Hcbh_column]],
      maxlad_dptf = merged_df1[[max_dptf_column]],
      maxlad_Hdptf = merged_df1[[max_Hdepth_column]],
      maxlad_effdist = 0
    )
    names(max_df) <- c("maxlad_Hcbh","maxlad_dptf", "maxlad_Hdptf","maxlad_effdist")
  }

  Hdist1<-"Hdist1"

  if (max_Hdist_column %in% colnames(merged_df1) && !max_effdist_column %in% colnames(merged_df1) && Hdist1 %in% colnames(merged_df1)) {
    if (merged_df1$Hdist1 > 1.5) {
      max_df <- data.frame(
        maxlad_Hcbh = merged_df1[[max_Hcbh_column]],
        maxlad_Hdist = merged_df1[[max_Hdist_column]],
        maxlad_Hdptf = merged_df1[[max_Hdepth_column]],
        maxlad_dptf = merged_df1[[max_dptf_column]],
        maxlad_effdist = 0
      )
      names(max_df) <- c("maxlad_Hcbh","maxlad_Hdist", "maxlad_Hdptf","maxlad_dptf","maxlad_effdist")
    }}


  if (max_Hdist_column %in% colnames(merged_df1) && !max_effdist_column %in% colnames(merged_df1) && Hdist1 %in% colnames(merged_df1) && max_effdist_column_minusone %in% colnames(merged_df1)) {
    if (merged_df1$Hdist1 == 1.5) {
      max_df <- data.frame(
        maxlad_Hcbh = merged_df1[[max_Hcbh_column]],
        maxlad_Hdist = merged_df1[[max_Hdist_column]],
        maxlad_Hdptf = merged_df1[[max_Hdepth_column]],
        maxlad_dptf = merged_df1[[max_dptf_column]],
        maxlad_effdist = merged_df1[[max_effdist_column_minusone]]
      )
      names(max_df) <- c("maxlad_Hcbh","maxlad_Hdist", "maxlad_Hdptf","maxlad_dptf","maxlad_effdist")
    }}

  if (max_Hdist_column %in% colnames(merged_df1) && !max_effdist_column %in% colnames(merged_df1) && Hdist1 %in% colnames(merged_df1) && !max_effdist_column_minusone %in% colnames(merged_df1)) {
    if (merged_df1$Hdist1 == 1.5) {
      max_df <- data.frame(
        maxlad_Hcbh = merged_df1[[max_Hcbh_column]],
        maxlad_Hdist = merged_df1[[max_Hdist_column]],
        maxlad_Hdptf = merged_df1[[max_Hdepth_column]],
        maxlad_dptf = merged_df1[[max_dptf_column]],
        maxlad_effdist = 0
      )
      names(max_df) <- c("maxlad_Hcbh","maxlad_Hdist", "maxlad_Hdptf","maxlad_dptf","maxlad_effdist")
    }}



  maxlad_Hdist<-"maxlad_Hdist"

  if (exists("max_df") && maxlad_Hdist %in% colnames(max_df)) {
    if(max_Hdist_column %in% colnames(merged_df1) && !is.na(max_df[[maxlad_Hdist]])) {
      if (max_df$maxlad_Hdist > max_df$maxlad_Hdptf) {

        max_df$maxlad_Hdist<- (max_df$maxlad_Hdptf - max_df$maxlad_dptf)-1
      } }

    if ("maxlad_Hdist"%in% colnames(max_df)) {
      if(max_df$maxlad_Hdist < 1 ) {
        max_df$maxlad_Hdist<- 0
      }}

    if (suffix == 1 && merged_df1$Hcbh1 == 1.5){
      max_df$maxlad_effdist<- 0
    }

    maxlad_Hcbh<-"maxlad_Hcbh"
    maxlad_effdist<-"maxlad_effdist"
    maxlad_Hdist<-"maxlad_Hdist"

    if (maxlad_Hcbh %in% colnames(max_df) && !is.na(max_df[[maxlad_effdist]])) {
      if (max_df$maxlad_Hcbh[1] == 1.5 || (maxlad_Hdist %in% colnames(max_df) && max_df$maxlad_Hdist[1] == 1.5)) {
        max_df$maxlad_effdist[1] <- 0
      }
    }

    if (maxlad_Hcbh %in% colnames(max_df) && !is.na(max_df[[maxlad_Hdist]])) {
      max_df$maxlad_Hdist[1] <- max_df$maxlad_Hcbh[1] -1
    }
  }

  merged_df1<-cbind(merged_df1, max_df)

  #######################################

  if (length(all_effdist_cols) > 0 && length(hdist_cols) > 0) {

    # Extract relevant columns
    hcbh_cols <- grep("^Hcbh\\d+$", colnames(merged_df1), value = TRUE)
    effdist_cols <- grep("^effdist\\d+$", colnames(merged_df1), value = TRUE)
    hdist_cols <- grep("^Hdist\\d+$", colnames(merged_df1), value = TRUE)

    # Determine number of effdist columns that should be present
    if (merged_df1$Hcbh1 > 1.5 && merged_df1$Hdist1 > 1.5) {
      effdist_required <- length(hcbh_cols)
    } else if (merged_df1$Hcbh1 == 1.5 || (merged_df1$Hcbh1 > 1.5 && merged_df1$Hdist1 == 1.5)) {
      effdist_required <- length(hcbh_cols) - 1
    }

    # Determine the effdist columns to remove
    effdist_to_remove <- tail(effdist_cols, -effdist_required)

    # Repeat for Hdist columns
    if (merged_df1$Hcbh1 > 1.5) {
      hdist_required <- length(hcbh_cols)
    } else if (merged_df1$Hcbh1 == 1.5) {
      hdist_required <- length(hcbh_cols) - 1
    }

    hdist_to_remove <- tail(hdist_cols, -hdist_required)

    # Combine the columns to remove
    cols_to_removefi <- c(effdist_to_remove, hdist_to_remove)

    # Remove those columns
    merged_df1 <- merged_df1[ , !(names(merged_df1) %in% cols_to_removefi)]

    ######################################

    # Extracting column names related to Hcbh, Hdist, Hdptf, and effdist
    hcbh_cols <- grep("^Hcbh", colnames(merged_df1), value = TRUE)
    hdist_cols <- grep("^Hdist", colnames(merged_df1), value = TRUE)
    hdptf_cols <- grep("^Hdptf", colnames(merged_df1), value = TRUE)
    effdist_cols <- grep("^effdist", colnames(merged_df1), value = TRUE)

    # Identify the effdist column that matches maxlad_effdist before calculations
    matching_index <- which(merged_df1[1, effdist_cols] == merged_df1$maxlad_effdist)[1]

    if (merged_df1$Hcbh1 > 2.5) {
      for (j in seq_along(effdist_cols)) {
        if (j == 1) {
          merged_df1[1, effdist_cols[j]] <- merged_df1[[hcbh_cols[j]]] - 1.5
        } else {
          merged_df1[1, effdist_cols[j]] <- merged_df1[[hdist_cols[j]]] - merged_df1[[hdptf_cols[j - 1]]]
        }
      }
    }

    if (merged_df1$Hcbh1 <= 2.5 && length(hdist_cols) > 0) {
      for (j in seq_along(effdist_cols)) {
        if (merged_df1$Hdist1 == 1.5) {
          merged_df1[1, effdist_cols[j]] <- merged_df1[[hdist_cols[j+1]]] - merged_df1[[hdptf_cols[j]]]
        }
      }
    }

    # If there's only one effdist column, update maxlad_effdist with its value
    if (length(effdist_cols) == 1) {
      merged_df1$maxlad_effdist <- merged_df1[1, effdist_cols]
    } else if (!is.na(matching_index)) {
      merged_df1$maxlad_effdist <- merged_df1[1, effdist_cols[matching_index]]
    }

    if (length(effdist_cols) == 1 && length(hdist_cols) == 1 && merged_df1$Hcbh1 <= 2.5) {
      merged_df1[1, effdist_cols[1]] <- merged_df1[[hdist_cols[1]]] - merged_df1[[hdptf_cols[1]]]
      merged_df1$maxlad_effdist <- merged_df1[[hdist_cols[1]]] - merged_df1[[hdptf_cols[1]]]
    }

    # Remove sub-columns if they exist
    cols_to_remove <- grep("^maxlad_effdist\\.", colnames(merged_df1), value = TRUE)
    merged_df1 <- merged_df1[, !names(merged_df1) %in% cols_to_remove]

  }

  # Remove sub-columns if they exist
  cols_to_remove1 <- grep("^max1", colnames(merged_df1), value = TRUE)
  merged_df1 <- merged_df1[, !names(merged_df1) %in% cols_to_remove1]

  ##########################################################################
  all_LAD<-merged_df
  effective_LAD<-merged_df1
  # Return them in a list
  return(list(df1 = all_LAD, df2 = effective_LAD))

}




