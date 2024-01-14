#' Remove non-consecutive fuel layers with a Leaf Area Density (LAD) less than a specified threshold
#' @description This function removes non-consecutive fuel layers with a Leaf Area Density (LAD) less than a specified threshold and
#' recalculates the distances and the depth of remaining fuel layers. Don't run it. It is an internal function.
#' @usage remove_no_flayer_noconsec(df, threshold = 25)
#' @param df Internal data frame derived from [get_layers_lad()] function.
#' @param threshold Numeric value for the minimum required LAD percentage in a fuel layer. The default threshold is 25.
#' @return No return value. The function is called for side effects.
#' @author Olga Viedma, Carlos Silva and JM Moreno
#' @examples
#' # remove_no_flayer_noconsec() removes fuel layers with a LAD percentage less than a threshold
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
#' @export
remove_no_flayer_noconsec <-function(df, threshold=25) {

  merged_df1<-df

  lad_columns2 <- grep("^Hcbh\\d+_Hdptf\\d+$", names(merged_df1),value=TRUE)
  lad_suffixes <- as.numeric(str_extract(lad_columns2, "\\d+$"))

  threshold_value <- threshold

  # Find the columns with value less than 5
  cols_to_remove <- sapply(merged_df1[, lad_columns2], function(col) any(col < threshold_value))
  suffixes_to_remove <- sort(as.numeric(stringr::str_extract(names(cols_to_remove[cols_to_remove]), "\\d+$")))

  lad_no_remove <- sapply(merged_df1[, lad_columns2], function(col) any(col > threshold_value))
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

      cols_to_remove_noeffdist2 <- matching_cols3[!grepl("^(effdist|Hdist|treeID)", matching_cols3)]
      cols_to_remove_noeffdist1 <- matching_cols3[!grepl("^(effdist|treeID)", matching_cols3)]

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
              !next_effdist_col %in% colnames(merged_df1) && dist1_value > 1 && merged_df1$Hcbh1 <=2.5) {

            merged_df1[first_col_name] <- dptf_sum + effdist_sum
            merged_df1 <- merged_df1[, !(names(merged_df1) %in% cols_to_remove_noeffdist2)]
          }}

        #################################

        cols_to_check <- c(first_cbh_gt5_col,dist_col1)
        if (all(cols_to_check %in% colnames(merged_df1))) {

          if (isTRUE(first_suffix_eq1 == first(lad_suffixes)) && first_cbh_gt5_value  > 1.5 &&
              !next_effdist_col %in% colnames(merged_df1) && dist1_value > 1 && merged_df1$Hcbh1 > 2.5) {

            merged_df1[first_col_name] <- dptf_sum + effdist_sum
            merged_df1 <- merged_df1[, !(names(merged_df1) %in% cols_to_remove_noeffdist1)]
          }}

        #################################

        cols_to_check <- c(first_cbh_gt5_col, first(Hcbh_cols),Hdist_next_colname)
        if (all(cols_to_check %in% colnames(merged_df1))) {

          if (isTRUE(first_suffix_eq1 == first(lad_suffixes)) &&
              first_cbh_gt5_value  > 1.5 && Hcbh1<= 2.5 &&
              !Hdist_next2_colname %in% colnames(merged_df1) &&
              Hdist_next_value < first_cbh_gt5_value && all (sapply(all_effdist_values, function(x) x != "1"))) {

            merged_df1[first_col_name] <- dptf_sum + effdist_sum
            merged_df1 <- merged_df1[, !(names(merged_df1) %in% cols_to_remove_noeffdist2)]
          }}

        #################################

        cols_to_check <- c(first_cbh_gt5_col, first(Hcbh_cols))
        if (all(cols_to_check %in% colnames(merged_df1))) {

          if (isTRUE(first_suffix_eq1 == first(lad_suffixes)) &&
              first_cbh_gt5_value  > 1.5 && Hcbh1<= 2.5 &&
              !Hdist_next2_colname %in% colnames(merged_df1) &&
              all (sapply(all_effdist_values, function(x) x != "1"))) {

            merged_df1[first_col_name] <- dptf_sum + effdist_sum
            merged_df1 <- merged_df1[, !(names(merged_df1) %in% cols_to_remove_noeffdist2)]
          }}

        #################################

        cols_to_check <- c(first_cbh_gt5_col,Hdist_next_colname)
        if (all(cols_to_check %in% colnames(merged_df1))) {


          if (isTRUE(first_suffix_eq1 == first(lad_suffixes) &&
                     first_cbh_gt5_value  > 1.5 &&
                     (Hdist_next_value > first_cbh_gt5_value) &&
                     !Hdist_next2_colname %in% colnames(merged_df1) &&
                     all(sapply(all_effdist_values, function(x) x != "1")) &&
                     length(lad_values_noremove) > 1) ) {

            merged_df1[first_col_name] <- dptf_sum + effdist_sum
            merged_df1 <- merged_df1[, !(names(merged_df1) %in% cols_to_remove_noeffdist2)]
          }}

        #################################

        cols_to_check <- c(first_cbh_gt5_col, first(Hcbh_cols),next_effdist_col,Hdist_next2_colname)
        if (all(cols_to_check %in% colnames(merged_df1))) {

          if (isTRUE(first_suffix_eq1 == first(lad_suffixes)) &&
              first_cbh_gt5_value > 1.5 && Hcbh1<= 2.5 &&  ## first distance == 0.5
              next_effdist_value==1 &&
              Hdist_next2_value > first_cbh_gt5_value && any(sapply(all_effdist_values, function(x) x == "1"))) {

            merged_df1[first_col_name] <- dptf_sum + effdist_sum
            merged_df1 <- merged_df1[, !(names(merged_df1) %in% cols_to_remove_noeffdist2)]
          }}

        #################################

        cols_to_check <- c(first_cbh_gt5_col, first(Hcbh_cols),next_effdist_col,Hdist_next2_colname)
        if (all(cols_to_check %in% colnames(merged_df1))) {

          if (isTRUE(first_suffix_eq1 == first(lad_suffixes)) &&
              first_cbh_gt5_value > 1.5 && Hcbh1<= 2.5 &&  ## first distance == 0.5
              next_effdist_value > 1 &&
              Hdist_next2_value > first_cbh_gt5_value && all (sapply(all_effdist_values, function(x) x != "1"))) {

            merged_df1[first_col_name] <- dptf_sum + effdist_sum
            merged_df1 <- merged_df1[, !(names(merged_df1) %in% cols_to_remove_noeffdist2)]
          }}

        #################################

        cols_to_check <- c(first_cbh_gt5_col, first(Hcbh_cols))
        if (all(cols_to_check %in% colnames(merged_df1))) {

          if (isTRUE(first_suffix_eq1 == first(lad_suffixes)) &&
              first_cbh_gt5_value > 1.5 && Hcbh1<= 2.5 &&  ## first distance == 0.5
              !next_effdist_col %in% colnames(merged_df1) &&
              !Hdist_next2_colname %in% colnames(merged_df1) && all (sapply(all_effdist_values, function(x) x != "1"))) {

            merged_df1[first_col_name] <- dptf_sum + effdist_sum
            merged_df1 <- merged_df1[, !(names(merged_df1) %in% cols_to_remove_noeffdist2)]
          }}

        #################################

        cols_to_check <- c(first_cbh_gt5_col, first(Hcbh_cols),next_effdist_col,Hdist_next2_colname,next_hdptf_noremove_col)
        if (all(cols_to_check %in% colnames(merged_df1))) {

          if (isTRUE(first_suffix_eq1 == first(lad_suffixes)) &&
              first_cbh_gt5_value  > 1.5 && Hcbh1<= 2.5 &&  ## first distance == 0.5
              Hdist_next2_value < next_hdptf_noremove_value &&
              length(lad_values_noremove) == 1) {

            merged_df1[first_col_name] <- dptf_sum + effdist_sum + merged_df1[[next_effdist_col]]
            merged_df1 <- merged_df1[, !(names(merged_df1) %in% cols_to_remove_noeffdist2)]
            merged_df1[[next_effdist_col]] <- NULL
          }}

        #################################

        cols_to_check <- c(first_cbh_gt5_col, first(Hcbh_cols),next_effdist_col,Hdist_next_colname, Hdist_next2_colname,next_hdptf_noremove_col)
        if (all(cols_to_check %in% colnames(merged_df1))) {

          if (isTRUE(first_suffix_eq1 == first(lad_suffixes) &&
                     first_cbh_gt5_value  > 1.5 && Hcbh1<= 2.5 &&  ## first distance == 0.5
                     (Hdist_next_value < first_cbh_gt5_value) &&
                     Hdist_next2_value > next_hdptf_noremove_value &&
                     any(sapply(all_effdist_values, function(x) x == "1")) &&
                     length(lad_values_noremove) > 1)) {

            merged_df1[first_col_name] <- dptf_sum + effdist_sum + merged_df1[[next_effdist_col]]
            merged_df1 <- merged_df1[, !(names(merged_df1) %in% cols_to_remove_noeffdist2)]
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
            #merged_df1[[next_effdist_col]] <- NULL
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
                     first_cbh_gt5_value  > 1.5 && Hcbh1<= 2.5 &&
                     (Hdist_next_value < first_cbh_gt5_value) &&
                     !Hdist_next2_colname %in% colnames(merged_df1) &&
                     length(lad_values_noremove) == 1)) {

            merged_df1[first_col_name] <- dptf_sum + effdist_sum  + merged_df1[[next_effdist_col]]
            merged_df1 <- merged_df1[, !(names(merged_df1) %in% cols_to_remove_noeffdist2)]
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

    lad_columns2 <- grep("^Hcbh\\d+_Hdptf\\d+$", names(merged_df1),value=TRUE)
    lad_suffixes <- as.numeric(str_extract(lad_columns2, "\\d+$"))

    threshold_value <- threshold
    lad_no_remove <- sapply(merged_df1[, lad_columns2], function(col) any(col > threshold_value))
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

    #################################################

    # Find the columns with value less than  athreshold
    threshold_value <- threshold
    cols_to_remove <- sapply(merged_df1[, lad_columns2], function(col) any(col < threshold_value))
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


    if (any(suffixes_to_remove > 1)) {

      last_col_name <- paste0("effdist", last_no_consenc_suffix)

      if (last_col_name %in% colnames(merged_df1)) {

        lad_columns2 <- grep("^Hcbh\\d+_Hdptf\\d+$", names(merged_df1),value=TRUE)
        lad_suffixes <- as.numeric(str_extract(lad_columns2, "\\d+$"))

        threshold_value <- threshold
        cols_to_remove1 <- sapply(merged_df1[, lad_columns2], function(col) any(col < threshold_value))
        suffixes_to_remove2 <- sort(as.numeric(stringr::str_extract(names(cols_to_remove1[cols_to_remove1]), "\\d+$")))

        lad_no_remove <- sapply(merged_df1[, lad_columns2], function(col) any(col > threshold_value))
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

        suffix<-suffixes_to_remove2[1]

          pattern2 <- paste0(suffix, "$")
          matching_cols3 <- grep(pattern2, names(merged_df1), value=TRUE)
          cols_to_remove_noeffdist <- matching_cols3[!grepl("^(effdist|treeID)", matching_cols3)]
          cols_to_remove_noeffdist1 <- matching_cols3[!grepl("^(effdist|Hdist|treeID)", matching_cols3)]## fro Hcbh <=2.5

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

          cols_to_check <- c(previous_hdptf_noremove_col,effdist_cols)
          if (all(cols_to_check %in% colnames(merged_df1))) {

            if (isTRUE(suffix != first(lad_suffixes) && (Hdist_values > previous_hdptf_noremove_value)) &&
                suffix >= last(lad_suffixes) ) {

              merged_df1 <- merged_df1[, !(names(merged_df1) %in% cols_to_remove_noeffdist)]
              cols_to_remove <- match(effdist_cols, names(merged_df1))
              merged_df1<- merged_df1[,-cols_to_remove]

            }}

          ##############################

          cols_to_check <- c(next_effdist_col, previous_effdist_col, Hdist_previous_colname, Hdist_next_colname, previous_hdptf_noremove_col,
                             next_cbh_noremove_col,first_Hcbh_col)
          if (all(cols_to_check %in% colnames(merged_df1))) {

            if (isTRUE(suffix != first(lad_suffixes) &&
                       (Hdist_previous_value > previous_hdptf_noremove_value) &&
                       (Hdist_next_value < next_cbh_noremove_value)) &&
                length (suffixes_no_remove) > 1 &&
                first_Hcbh_value <= 2.5) {


              merged_df1[sufix_col_name] <- dptf_sum + effdist_sum + merged_df1[[previous_effdist_col]] + merged_df1[[next_effdist_col]]
              merged_df1[[previous_effdist_col]] <- NULL
              merged_df1[[next_effdist_col]] <- NULL

              merged_df1<- merged_df1[, !(names(merged_df1) %in% cols_to_remove_noeffdist1)]


            } else if (previous_effdist_col %in% colnames(merged_df1) && !next_effdist_col %in% colnames(merged_df1)) {
              merged_df1[sufix_col_name] <- dptf_sum + effdist_sum + merged_df1[[previous_effdist_col]]
              merged_df1[[previous_effdist_col]] <- NULL

              merged_df1<- merged_df1[, !(names(merged_df1) %in% cols_to_remove_noeffdist1)]
            }}

          ###############################

          cols_to_check <- c(next_effdist_col, previous_effdist_col, Hdist_col, Hdist_previous_colname, Hdist_next_colname, previous_hdptf_noremove_col, next_cbh_noremove_col,first_Hcbh_col)
          if (all(cols_to_check %in% colnames(merged_df1))) {

            if (isTRUE(suffix != first(lad_suffixes) &&
                       (Hdist_values > previous_hdptf_noremove_value) &&
                       (Hdist_next_value < next_cbh_noremove_value)) &&
                length (suffixes_no_remove) > 1  &&
                first_Hcbh_value <= 2.5 &&
                Hdist_previous_value > previous_hdptf_noremove_value) {


              merged_df1[sufix_col_name] <- dptf_sum + effdist_sum + merged_df1[[previous_effdist_col]] + merged_df1[[next_effdist_col]]
              merged_df1[[previous_effdist_col]] <- NULL
              merged_df1[[next_effdist_col]] <- NULL

              merged_df1<- merged_df1[, !(names(merged_df1) %in% cols_to_remove_noeffdist1)]


            }  else if (previous_effdist_col %in% colnames(merged_df1) && ! next_effdist_col %in% colnames(merged_df1)) {
              merged_df1[sufix_col_name] <- dptf_sum + effdist_sum + merged_df1[[previous_effdist_col]]
              merged_df1[[previous_effdist_col]] <- NULL

              merged_df1<- merged_df1[, !(names(merged_df1) %in% cols_to_remove_noeffdist1)]

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

              merged_df1<- merged_df1[, !(names(merged_df1) %in% cols_to_remove_noeffdist)]
            }}

          ##############################

          cols_to_check <- c(next_effdist_col, previous_effdist_col, Hdist_col, Hdist_previous_colname, Hdist_next_colname, previous_hdptf_noremove_col, next_cbh_noremove_col)
          if (all(cols_to_check %in% colnames(merged_df1))) {

            if (isTRUE(suffix != first(lad_suffixes) &&
                       (Hdist_values > previous_hdptf_noremove_value) &&
                       (Hdist_next_value < next_cbh_noremove_value)) &&
                length (suffixes_no_remove) > 1  &&
                first_Hcbh_value <= 2.5 &&
                Hdist_previous_value < previous_hdptf_noremove_value) {


              merged_df1[sufix_col_name] <- dptf_sum + effdist_sum +  merged_df1[[next_effdist_col]]
              merged_df1[[next_effdist_col]] <- NULL

              merged_df1<- merged_df1[, !(names(merged_df1) %in% cols_to_remove_noeffdist1)]


            }  else if (previous_effdist_col %in% colnames(merged_df1) && ! next_effdist_col %in% colnames(merged_df1)) {
              merged_df1[last_col_name] <- dptf_sum + effdist_sum + merged_df1[[previous_effdist_col]]
              merged_df1[[previous_effdist_col]] <- NULL

              merged_df1<- merged_df1[, !(names(merged_df1) %in% cols_to_remove_noeffdist1)]

            }}

          ##############################

          cols_to_check <- c(previous_effdist_col, Hdist_col, Hdist_previous_colname, Hdist_next_colname, previous_hdptf_noremove_col, next_cbh_noremove_col)
          if (all(cols_to_check %in% colnames(merged_df1))) {

            if (isTRUE(suffix != first(lad_suffixes) &&
                       (Hdist_values > previous_hdptf_noremove_value) &&
                       (Hdist_next_value < next_cbh_noremove_value)) &&
                !next_effdist_col %in% colnames(merged_df1) &&
                length (suffixes_no_remove) > 1  &&
                first_Hcbh_value <= 2.5) {

              merged_df1[sufix_col_name] <- dptf_sum + effdist_sum + merged_df1[[previous_effdist_col]]
              merged_df1[[previous_effdist_col]] <- NULL

              merged_df1<- merged_df1[, !(names(merged_df1) %in% cols_to_remove_noeffdist1)]

            }}

          ##############################

          cols_to_check <- c( previous_effdist_col, Hdist_previous_colname, Hdist_next_colname, previous_hdptf_noremove_col, next_hdptf_noremove_col)
          if (all(cols_to_check %in% colnames(merged_df1))) {

            if (isTRUE(suffix != first(lad_suffixes) &&
                       (Hdist_previous_value > previous_hdptf_noremove_value) && (
                         Hdist_next_value > next_hdptf_noremove_value))) {

              merged_df1[sufix_col_name] <- dptf_sum + effdist_sum +  merged_df1[[previous_effdist_col]]
              merged_df1[[previous_effdist_col]] <- NULL

              merged_df1<- merged_df1[, !(names(merged_df1) %in% cols_to_remove_noeffdist1)]
            }}

          ##############################

          cols_to_check <- c( previous_effdist_col, Hdist_previous_colname, previous_hdptf_noremove_col)
          if (all(cols_to_check %in% colnames(merged_df1))) {

            if (isTRUE(suffix != first(lad_suffixes) &&
                       (Hdist_previous_value > previous_hdptf_noremove_value) &&
                       (!Hdist_next_colname %in% colnames(merged_df1)))) {

              merged_df1[sufix_col_name] <- dptf_sum + effdist_sum +  merged_df1[[previous_effdist_col]]
              merged_df1[[previous_effdist_col]] <- NULL

              merged_df1<- merged_df1[, !(names(merged_df1) %in% cols_to_remove_noeffdist1)]

            }}

          ##############################

          cols_to_check <- c( previous_effdist_col, Hdist_col, previous_hdptf_noremove_col)
          if (all(cols_to_check %in% colnames(merged_df1))) {

            if (isTRUE(suffix != first(lad_suffixes) &&
                       (Hdist_values > previous_hdptf_noremove_value) &&
                       (!next_effdist_col %in% colnames(merged_df1))&&
                       first_Hcbh_value <= 2.5)) {


              merged_df1[sufix_col_name] <- dptf_sum + effdist_sum +  merged_df1[[previous_effdist_col]]
              merged_df1[[previous_effdist_col]] <- NULL

              merged_df1<- merged_df1[, !(names(merged_df1) %in% cols_to_remove_noeffdist1)]
            }}

          merged_df1<- get_renamed_df (merged_df1)


      } else {

        pattern <- paste0(".*", last_no_consenc_suffix, "$")  # Create a pattern that matches the suffix at the end of a string
        selected_cols <- grep(pattern, names(merged_df1), value=TRUE)
        merged_df1 <- merged_df1[, !(names(merged_df1) %in% selected_cols)]

        second_last_no_consenc_suffix <-last_no_consenc_suffix -1
        effdist_to_remove <- paste0("effdist", second_last_no_consenc_suffix)
        merged_df1 <- merged_df1[, !(names(merged_df1) %in% effdist_to_remove)]
      }
    }


    merged_df1<- get_renamed_df (merged_df1)

  }

    ################################################


  return(merged_df1)
}
