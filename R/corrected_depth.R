#'
#' Effective fuel layers depth
#'
#' @description This function recalculates fuel layers depth after removing distances = 1 m
#' @usage get_real_depths (effective_fbh)
#' @param effective_fbh tree metrics with the recalculated base height of fuel layers after considering distances > 1 m  (output of [get_real_fbh()] function).
#' An object of the class text
#' @return A data frame giving the fuel layers depth after removing distances = 1 m.
#' @author Olga Viedma, Carlos Silva and JM Moreno
#'
#' @details
#' # List of tree metrics:
#' \itemize{
#' \item treeID: tree ID with strings and numeric values
#' \item treeID1: tree ID with only numeric values
#' \item dist: Distance between consecutive fuel layers (m)
#' \item Hdist - Height of the distance between consecutive fuel layers (m)
#' \item Hcbh - Height of the base of each fuel layer (m)
#' \item dptf - Depth of fuel layers (m) after removing distances equal 1 m
#' \item Hdptf - Height of the depth of fuel layers (m) after removing distances equal 1 m
#' \item max_height - Maximum height of the tree profile
#' }
#' @examples
#' ## Not run:
#' library(magrittr)
#' library(tidyr)
#' library(dplyr)
#'
#' # Load the effective_fbh object
#' if (interactive()) {
#'   effective_fbh <- get_real_fbh()
#'   LadderFuelsR::effective_fbh$treeID <- factor(LadderFuelsR::effective_fbh$treeID)
#'
#' # Tree metrics derived from get_real_fbh() function
#' effective_fbh$treeID <- factor(effective_fbh$treeID)
#'
#' trees_name1 <- as.character(effective_fbh$treeID)
#' trees_name2 <- factor(unique(trees_name1))
#'
#' depth_metrics_corr_list <- list()
#'
#' for (i in levels(trees_name2)){
#'   # Filter data for each tree
#'   tree3 <- effective_fbh |> dplyr::filter(treeID == i)
#'   # Get real depths for each tree
#'   depth_metrics_corr <- get_real_depths(tree3)
#'   depth_metrics_corr_list[[i]] <- depth_metrics_corr
#' }
#'
#' # Combine depth values for all trees
#' effective_depth <- dplyr::bind_rows(depth_metrics_corr_list)
#' }
#' ## End(Not run)
#'
#' @export get_real_depths
#' @importFrom dplyr select_if group_by summarise summarize mutate arrange rename rename_with filter slice slice_tail ungroup distinct
#' across matches row_number all_of vars n
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
#' @include gap_fbh.R
#' @include distances_calculation.R
#' @include depths_calculation.R
#' @include corrected_base_heights.R
#' @keywords internal
#' @seealso \code{\link{get_renamed0_df}}
get_real_depths <- function (effective_fbh) {

  #remove the columns from the dataframe df2a which contain only NA values.
  df<- effective_fbh
  df2a <- df[, colSums(!is.na(df)) > 0]
  df2a <- df2a[, !(names(df2a) %in% "depth01")]

  print(paste("Unique treeIDs:", paste(unique(df2a$treeID), collapse = ", ")))

  ##check if the first "dist" column has a value greater than 1. If this is true, then it will remove depth0 in computation

  # Identify the "dist" and "depth" column names
  dist_cols <- grep("^dist", colnames(df2a), value = TRUE)

  if (length(dist_cols) == 0) {
    df2a$dist1 <- 1
  }

  if (length(dist_cols) == 1 && all(df2a$dist1==0)) {
    df2a$dist1 <- 1
  }

  for (col in names(df2a)) {
    # Check if column name starts with "dist"
    if (startsWith(col, "dist")) {
      # Replace 0 values with NA
      df2a[df2a[, col] == 0, col] <- 1}}


  # First, create a vector of the Hdepth column names
  hdepth_cols <- grep("^Hdepth", names(df2a), value = TRUE)
  zero_cols <- character(0)
  # If there's more than one Hdepth column
  if(length(hdepth_cols) > 1) {
    # Determine which Hdepth columns have a 0
    zero_cols <- hdepth_cols[apply(df2a[hdepth_cols], 2, function(x) all(x == 0))]
  }

  # Drop those columns
  df4 <- df2a[!names(df2a) %in% zero_cols]


  # RENAME: Extract unique prefixes
  prefixes <- unique(gsub("([a-zA-Z]+).*", "\\1", names(df4)))

  # Rename the columns based on the extracted prefixes
  for (prefix in prefixes) {
    # Identify columns with the current prefix
    cols <- grep(paste0("^", prefix), names(df4))

    # Generate new column names with consecutive suffixes
    new_names <- paste0(prefix, 1:length(cols))

    # Assign new names to the columns
    names(df4)[cols] <- new_names
  }

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
  df_transposed0 <- df_transposed0 %>%
    dplyr::filter(!is.na(Hdepth))

  ############################################
  df_transposed <- df_transposed0

  # Special case for the last row if its dist is NA
  if (is.na(df_transposed$dist[nrow(df_transposed)])) {
    df_transposed$dist[nrow(df_transposed)] <- 1
  }

  if (nrow(df_transposed) > 1) {
    dptf<-df_transposed$dptf
    Hdepth1<-df_transposed$Hdepth1
    df_transposed$dptf <- NA
    df_transposed$Hdepth1 <- NA

    if (any(df_transposed$dist > 1)) {

      # Find the position of the first occurrence of dist > 1
      first_dist_gt_one <- which(df_transposed$dist > 1)[1]
      # If there's no dist value > 1, set the position to beyond the end of the dataframe
      if (is.na(first_dist_gt_one)) {
        first_dist_gt_one <- nrow(df_transposed) + 1
      }
      # Check if all values before this position are 1
      all_ones_before <- all(df_transposed$dist[1:(first_dist_gt_one - 1)] == 1)
      # Check if all values from this position onwards are > 1
      all_gt_one_after <- all(df_transposed$dist[first_dist_gt_one:length(df_transposed$dist)] > 1)
      # Combine both conditions
      conditions_met <- all_ones_before & all_gt_one_after

      if (conditions_met && length(all_ones_before) > 1) {

        i <- 1

        while (i <= nrow(df_transposed)) {
          depth_sum <- 0
          dist_sum <- 0
          j = i

          while (j <= nrow(df_transposed) && df_transposed$dist[j] == 1) {
            depth_sum <- depth_sum + df_transposed$depth[j]
            dist_sum <- dist_sum + df_transposed$dist[j]
            j <- j + 1
          }

          df_transposed$dptf[i:(j-1)] <- dist_sum + depth_sum
          df_transposed$Hdepth1[i:(j-1)] <- df_transposed$Hdepth[j-1]

          if (j <= nrow(df_transposed) && df_transposed$dist[j] > 1) {
            df_transposed$dptf[j] <- df_transposed$depth[j]
            df_transposed$Hdepth1[j] <- df_transposed$Hdepth[j]
            j <- j + 1
          }

          i <- j
        }
      } else if (all(df_transposed$dist == 1)) {

        if (df_transposed$Hcbh[1] == 1.5) {
          df_transposed$dptf <- sum(df_transposed$depth) + sum(df_transposed$dist)
          df_transposed$Hdepth1 <- tail(df_transposed$Hdepth, n = 1)
        } else if (df_transposed$Hcbh[1] > 1.5) {
          df_transposed$dptf <- sum(df_transposed$depth[1:nrow(df_transposed)]) + sum(df_transposed$dist[2:nrow(df_transposed)])
          df_transposed$Hdepth1 <- tail(df_transposed$Hdepth, n = 1)
        }
      } }



    if (any(is.na(df_transposed$dptf)) || any(is.na(df_transposed$Hdepth1))) {

      start_row <- 1

      if (nrow(df_transposed) == 2 && df_transposed$Hcbh[1] > 1.5 && !is.na(df_transposed$dist[2]) && df_transposed$dist[2] > 1) {
        df_transposed$dptf[1] <- df_transposed$depth[1]
        df_transposed$Hdepth1[1] <- df_transposed$Hdepth[1]
        df_transposed$dptf[2] <- df_transposed$depth[2]
        df_transposed$Hdepth1[2] <- df_transposed$Hdepth[2]
      }
    }


    if (any(is.na(df_transposed$dptf)) || any(is.na(df_transposed$Hdepth1))) {

      df_transposed$dptf <- NA
      df_transposed$Hdepth1 <- NA

      if (df_transposed$Hcbh[1] > 1.5) {
        # Find rows where dist = 1 starting from the second row
        indices_dist_1 <- which(df_transposed$dist[-1] == 1) + 1

        if (length(indices_dist_1) > 0 && indices_dist_1[1] == 2 && all(diff(indices_dist_1) == 1)) {
          depth_sum <- df_transposed$depth[1] + sum(df_transposed$depth[indices_dist_1]) + sum(df_transposed$dist[indices_dist_1])
          df_transposed$dptf[c(1, indices_dist_1)] <- depth_sum
          df_transposed$Hdepth1[c(1, indices_dist_1)] <- df_transposed$Hdepth[max(indices_dist_1)]
        } else {
          df_transposed$dptf[1] <- df_transposed$depth[1]
          df_transposed$Hdepth1[1] <- df_transposed$Hdepth[1]
          df_transposed$dptf[indices_dist_1] <- df_transposed$depth[indices_dist_1]
          df_transposed$Hdepth1[indices_dist_1] <- df_transposed$Hdepth[indices_dist_1]
        }

        # Update rows where dist > 1
        indices_dist_gt_1 <- which(df_transposed$dist > 1 & seq_len(nrow(df_transposed)) != 1)
        df_transposed$dptf[indices_dist_gt_1] <- df_transposed$depth[indices_dist_gt_1]
        df_transposed$Hdepth1[indices_dist_gt_1] <- df_transposed$Hdepth[indices_dist_gt_1]
      }}



    if (any(is.na(df_transposed$dptf)) || any(is.na(df_transposed$Hdepth1))) {


      if (nrow(df_transposed) > 2 && df_transposed$Hcbh[1] == 1.5 && df_transposed$dist[1] == 1 &&
          !is.na(df_transposed$dist[2]) && df_transposed$dist[2] > 1) {
        df_transposed$dptf[1:2] <- df_transposed$depth[1] + df_transposed$dist[1] + df_transposed$depth[2]
        df_transposed$Hdepth1[1:2] <- df_transposed$Hdepth[2]
        start_row <- 3

      }

      i <- start_row

      while (i <= nrow(df_transposed)) {
        depth_sum <- 0
        dist_sum <- 0
        j = i

        # While dist is 1, sum the depths and dists
        while (j <= nrow(df_transposed) && df_transposed$dist[j] == 1) {
          depth_sum <- depth_sum + df_transposed$depth[j]
          dist_sum <- dist_sum + df_transposed$dist[j]
          j <- j + 1
        }

        # If we've reached beyond the dataframe or encountered a dist > 1
        if (j > nrow(df_transposed) || (j <= nrow(df_transposed) && df_transposed$dist[j] > 1)) {
          # If the next row has dist > 1 and is not the last row, add its depth to the sum
          if (j <= nrow(df_transposed) && df_transposed$dist[j] > 1 && j != nrow(df_transposed)) {
            depth_sum <- depth_sum + df_transposed$depth[j]
          }

          # Assign accumulated values to the sequence
          df_transposed$dptf[i:min(j, nrow(df_transposed))] <- dist_sum + depth_sum
          df_transposed$Hdepth1[i:min(j, nrow(df_transposed))] <- df_transposed$Hdepth[min(j, nrow(df_transposed))]

          # If we're on the last row and its dist value is 1
          if (j == nrow(df_transposed) && df_transposed$dist[j] == 1) {
            if (df_transposed$dist[j-1] > 1) {
              df_transposed$dptf[j] <- df_transposed$depth[j]
              df_transposed$Hdepth1[j] <- df_transposed$Hdepth[j]
            } else {
              preceding_dist_1_rows = which(df_transposed$dist[1:j] == 1)
              total_depth_for_last_rows = sum(df_transposed$depth[preceding_dist_1_rows])
              total_dist_for_last_rows = sum(df_transposed$dist[preceding_dist_1_rows])

              df_transposed$dptf[j] <- total_depth_for_last_rows + total_dist_for_last_rows
              df_transposed$Hdepth1[j] <- df_transposed$Hdepth[j]
            }
            break
          } else if (j == nrow(df_transposed) && df_transposed$dist[j] > 1) {
            df_transposed$dptf[j] <- df_transposed$depth[j]
            df_transposed$Hdepth1[j] <- df_transposed$Hdepth[j]
            break
          }
        }

        # Increment for next iteration
        i <- j + 1
      }}

  }

  ##########################################3

  df_transposed1 <- df_transposed

  if (nrow(df_transposed1) > 1) {

    df_transposed1 <- df_transposed1 %>%
      dplyr::mutate(dist = ifelse(is.na(dist), 1, dist)) %>%
      dplyr::filter(!is.na(dptf), Hdepth != 0.5)

    #df_transposed1 <- df_transposed1 %>% dplyr::select(-Hdepth)
    df_transposed1 <- df_transposed1 %>%
      dplyr::rename(Hdepth2 = Hdepth)%>%
      dplyr::rename(Hdepth = Hdepth1)


    ##########################################3

    #We'll begin from the second row. If dist is greater than 1, we'll keep its depth value aside. For the next rows where dist is 1, we'll sum their depth and dist #values. This summation continues until we either reach the end of the dataframe or encounter a row where dist is greater than 1 again. In the latter case, we add the previously kept aside depth value to our sum.we also update the Hdepth values for the rows that were affected by the computation using the last Hdepth2 value from each sequence.

    if (df_transposed1$Hcbh[1] > 1.5) {
      i <- 2

      while (i <= nrow(df_transposed1)) {
        # If current row's dist > 1, look ahead
        if (df_transposed1$dist[i] > 1) {
          depth_sum <- df_transposed1$depth[i]
          j <- i + 1
          last_Hdepth2 <- df_transposed1$Hdepth2[i]

          # Check the next rows as long as dist == 1
          while (j <= nrow(df_transposed1) && df_transposed1$dist[j] == 1) {
            depth_sum <- depth_sum + df_transposed1$depth[j] + df_transposed1$dist[j]
            last_Hdepth2 <- df_transposed1$Hdepth2[j]  # update the last Hdepth2 value
            j <- j + 1
          }

          # Update the dptf column
          df_transposed1$dptf[i:(j-1)] <- depth_sum
          # Update the Hdepth value
          df_transposed1$Hdepth[i:(j-1)] <- last_Hdepth2

          i <- j
        } else {
          i <- i + 1
        }
      }}


    #################################################

    df_transposed1 <- df_transposed1 %>%
      dplyr::group_by(Hdepth, dptf) %>%
      dplyr::filter(dist == max(dist)) %>%
      dplyr::slice(n()) %>%  # This will select the last row if there are multiple rows with the max dist.
      dplyr::ungroup()


    ###########################################################3
    Hcbh<-df_transposed1$Hcbh
    df_transposed2 <- df_transposed1 %>% dplyr::select(-Hcbh)

    # Extract Hcbh columns from df4
    Hcbh_cols <- grep("^Hcbh", colnames(df4), value = TRUE)
    hcbh_data <- df4[, Hcbh_cols, drop = FALSE]
    hcbh_data_unique <- unique(as.vector(t(hcbh_data)))

    # If all dptf values are the same OR dptf values are duplicated but Hdepth values differ
    if (length(unique(df_transposed2$dptf)) == 1 || (length(unique(df_transposed2$dptf)) < nrow(df_transposed2) && length(unique(df_transposed2$Hdepth)) != 1)) {
      # If there are more unique Hcbh values than rows, truncate the Hcbh values
      if (length(hcbh_data_unique) >= nrow(df_transposed2)) {
        df_transposed2$Hcbh <- hcbh_data_unique[1:nrow(df_transposed2)]
      } else {
        # Use all available Hcbh values and extend the last one for remaining rows
        num_remaining_rows <- nrow(df_transposed2) - length(hcbh_data_unique)
        extended_values <- rep(hcbh_data_unique[length(hcbh_data_unique)], num_remaining_rows)
        df_transposed2$Hcbh <- c(hcbh_data_unique, extended_values)
      }
    } else {
      # Previous logic for different dptf values

      unique_dptf <- unique(df_transposed2$dptf)
      dptf_counts <- table(df_transposed2$dptf)
      dptf_counts_df <- as.data.frame(dptf_counts)

      extended_hcbh_values <- vector("double", length = 0)

      for (i in seq_along(hcbh_data_unique)) {
        if (unique_dptf[i] %in% dptf_counts_df$Var1) {
          count <- as.numeric(dptf_counts_df[dptf_counts_df$Var1 == unique_dptf[i], 2])
        } else {
          count <- 0
        }

        hcbh_value <- hcbh_data_unique[i]
        extended_hcbh_values <- c(extended_hcbh_values, rep(hcbh_value, count))
      }

      if (length(unique_dptf) > length(hcbh_data_unique)) {
        for (i in (length(hcbh_data_unique) + 1):length(unique_dptf)) {
          if (unique_dptf[i] %in% dptf_counts_df$Var1) {
            count <- as.numeric(dptf_counts_df[dptf_counts_df$Var1 == unique_dptf[i], 2])
            hcbh_value <- tail(hcbh_data_unique, 1)
            extended_hcbh_values <- c(extended_hcbh_values, rep(hcbh_value, count))
          }
        }
      }

      df_transposed2$Hcbh <- extended_hcbh_values
    }


    # Reorder columns to match the original dataframe
    df_transposed2 <- df_transposed2[, names(df_transposed1)]
    Hdist<-df_transposed2$Hdist
    df_transposed2 <- df_transposed2 %>%
      fill(Hdist, .direction = "down")
  }

  #####################################################

  if(!exists("df_transposed2")) {

    df_transposed2 <-df_transposed1
  }

  if(!exists("df_transposed3")) {

    df_transposed3 <-df_transposed2
  }

  ########################################

  #We want to check if the first Hcbh value > 1.5 and the next dist value is 1.
  #If the condition is met, we should sum all depth values until we find a dist row value > 1.
  #This sum should update the depth of the first row.
  #We also need to update the Hdptf of the first row with the Hdepth value corresponding to the last dist value of 1.
  #After updating, remove all the rows with consecutive dist values of 1 after the first row.

  if (nrow(df_transposed3) > 1) {

    if(df_transposed3$Hcbh[1] > 1.5 && df_transposed3$dist[2] == 1 && any(df_transposed3$dist > 1) && nrow(df_transposed3) > 2) {

      sum_val <- df_transposed3$depth[1]  # Begin sum from the first row
      last_Hdepth <- df_transposed3$Hdepth[1]
      rows_to_remove <- NULL  # Start with an empty vector

      for(i in 2:nrow(df_transposed3)) {
        if(i <= nrow(df_transposed3) && df_transposed3$dist[i] == 1 && nrow(df_transposed3) > 2) {
          sum_val <- sum_val + df_transposed3$depth[i] + df_transposed3$dist[i]
          last_Hdepth <- df_transposed3$Hdepth2[i]
          rows_to_remove <- c(rows_to_remove, i)
        } else {
          break
        }
      }

      # Update the depth and Hdptf values of the first row after the sequence of dist=1
      df_transposed3$dptf[1] <- sum_val
      df_transposed3$Hdepth[1] <- last_Hdepth

      # Remove rows
      df_transposed3 <- df_transposed3[-rows_to_remove, ]
    }

    ################################################

    if(all(df_transposed3$dist == 1)) {
      # Preserve the first Hdist
      first_Hdist <- df_transposed3$Hdist[1]
      # Keep only the last row
      df_transposed3 <- tail(df_transposed3, 1)
      # Replace Hdist and depth with required values
      df_transposed3$Hdist <- first_Hdist
    }

  }
  ################################################

  if (nrow(df_transposed3) > 1) {

    if(df_transposed3$Hcbh[1] > 1.5 && df_transposed3$dist[2] == 1 && nrow(df_transposed3) == 2) {

      df_transposed3$Hdepth[1:2] <- df_transposed3$Hdepth[2]
      df_transposed3$dptf[1:2] <- df_transposed3$Hdepth[2] - df_transposed3$Hcbh[1]
    }


    if(df_transposed3$Hcbh[1] > 1.5 && df_transposed3$dist[2] == 1 && nrow(df_transposed3) == 2 && df_transposed3$Hdist[2] < df_transposed3$Hdepth[2] &&
       df_transposed3$dptf[1] != df_transposed3$dptf[2]) {

      sum_val <- 0  # Initialize sum_val to 0
      for(i in 1:nrow(df_transposed3)) {
        sum_val <- sum_val + df_transposed3$depth[i]
      }

      df_transposed3$dptf[2] <- sum_val + df_transposed3$dist[2]
      df_transposed3$Hdist[2] <- df_transposed3$Hdist[1]
      df_transposed3$dist[2] <- df_transposed3$dist[1]
      df_transposed3 <- df_transposed3[2, ]
    }

    if(df_transposed3$Hcbh[1] > 1.5 && df_transposed3$dist[2] > 1  && df_transposed3$Hdepth[1] >= df_transposed3$Hcbh[2]) {
      df_transposed3$Hdepth[1]<-df_transposed3$Hdist[2] - df_transposed3$dist[2]
    }

    ###############################################
    df_transposed4 <- df_transposed3 %>% dplyr::filter(!is.na(dptf))

    df_transposed4 <- df_transposed4 %>%
      dplyr::mutate(original_order = row_number()) %>%
      dplyr::group_by(dist,dptf, Hdepth) %>%
      dplyr::slice_tail(n = 1) %>%
      dplyr::ungroup()

    original_order<-df_transposed4$original_order

    df_transposed4 <- df_transposed4 %>%
      dplyr::arrange(original_order) %>%
      dplyr::select(-original_order)

    ################################################

    if(df_transposed4$Hcbh[1] == 1.5 && df_transposed4$dist[2] == 1 && df_transposed4$Hdepth[1] == df_transposed4$Hdepth[2] &&
       df_transposed4$dptf[2] > df_transposed4$dptf[1]){
      df_transposed4$dptf[1]<-df_transposed4$dptf[1] <-df_transposed4$dptf[2]
    }
    ##Check if there are still any duplicated Hdepth values

    df_transposed4 <- df_transposed4 %>%
      dplyr::group_by(Hdepth, dptf) %>%
      dplyr::filter(dist == max(dist)) %>%
      dplyr::slice(n()) %>%  # This will select the last row if there are multiple rows with the max dist.
      dplyr::ungroup()

    if(df_transposed4$Hcbh[1] > 1.5){
      df_transposed4$Hdist[1]<-df_transposed4$Hcbh [1] - 1
    }

    ######################################################

    df_transposed5 <- df_transposed4 %>% dplyr::select(-Hcbh)

    # Extract Hcbh columns from df4
    Hcbh_cols <- grep("^Hcbh", colnames(df4), value = TRUE)
    hcbh_data <- df4[, Hcbh_cols, drop = FALSE]
    hcbh_data_unique <- unique(as.vector(t(hcbh_data)))

    # If all dptf values are the same OR dptf values are duplicated but Hdepth values differ
    if (length(unique(df_transposed5$dptf)) == 1 || (length(unique(df_transposed5$dptf)) < nrow(df_transposed5) && length(unique(df_transposed5$Hdepth)) != 1)) {
      # If there are more unique Hcbh values than rows, truncate the Hcbh values
      if (length(hcbh_data_unique) >= nrow(df_transposed5)) {
        df_transposed5$Hcbh <- hcbh_data_unique[1:nrow(df_transposed5)]
      } else {
        # Use all available Hcbh values and extend the last one for remaining rows
        num_remaining_rows <- nrow(df_transposed5) - length(hcbh_data_unique)
        extended_values <- rep(hcbh_data_unique[length(hcbh_data_unique)], num_remaining_rows)
        df_transposed5$Hcbh <- c(hcbh_data_unique, extended_values)
      }
    } else {
      # Previous logic for different dptf values

      unique_dptf <- unique(df_transposed5$dptf)
      dptf_counts <- table(df_transposed5$dptf)
      dptf_counts_df <- as.data.frame(dptf_counts)

      extended_hcbh_values <- vector("double", length = 0)

      for (i in seq_along(hcbh_data_unique)) {
        if (unique_dptf[i] %in% dptf_counts_df$Var1) {
          count <- as.numeric(dptf_counts_df[dptf_counts_df$Var1 == unique_dptf[i], 2])
        } else {
          count <- 0
        }

        hcbh_value <- hcbh_data_unique[i]
        extended_hcbh_values <- c(extended_hcbh_values, rep(hcbh_value, count))
      }

      if (length(unique_dptf) > length(hcbh_data_unique)) {
        for (i in (length(hcbh_data_unique) + 1):length(unique_dptf)) {
          if (unique_dptf[i] %in% dptf_counts_df$Var1) {
            count <- as.numeric(dptf_counts_df[dptf_counts_df$Var1 == unique_dptf[i], 2])
            hcbh_value <- tail(hcbh_data_unique, 1)
            extended_hcbh_values <- c(extended_hcbh_values, rep(hcbh_value, count))
          }
        }
      }

      df_transposed5$Hcbh <- extended_hcbh_values
    }

    # Reorder columns to match the original dataframe
    df_transposed5 <- df_transposed5[, names(df_transposed1)]

    ###############################################
    #For each duplicated Hcbh value, you want to update the dist value with the maximum dist of the first occurrence of that duplicated Hcbh value.
    #For rows where Hcbh is duplicated:
    #You want to update the dist value to the maximum dist for that group.
    #You want to update the dptf value to the maximum dptf for that group.
    #You want to keep only the last row of these duplicated Hcbh values.
    #For rows where Hcbh is unique:


    df_transposed5 <- df_transposed5 %>%
      dplyr::group_by(Hcbh) %>%
      dplyr::mutate(is_duplicated = duplicated(Hcbh) | duplicated(Hcbh, fromLast = TRUE))

    is_duplicated<-df_transposed5$is_duplicated

    df_transposed5 <- df_transposed5 %>%
      # For duplicated Hcbh values, compute the max dptf and max dist
      dplyr::mutate(dptf = ifelse(is_duplicated, max(dptf, na.rm = TRUE), dptf),
             dist = ifelse(is_duplicated, max(dist, na.rm = TRUE), dist),
             Hdist = ifelse(is_duplicated,
                            min(Hdist[is_duplicated], na.rm = TRUE),
                            Hdist)) %>%
      # If Hcbh is duplicated, keep only the last row. Otherwise, keep all rows.
      dplyr::filter(ifelse(is_duplicated, row_number() == n(), TRUE)) %>%
      dplyr::ungroup() %>%
      dplyr::select(-is_duplicated)  # Drop the helper column

    Hdepth2<-df_transposed5$Hdepth2
    df_transposed5 <- df_transposed5 %>% dplyr::select(-Hdepth2)

  } else {
    df_transposed5 <-df_transposed3

    if("Hdepth2" %in%colnames(df_transposed5)){
      df_transposed5 <- df_transposed5 %>% dplyr::select(-Hdepth2)
    }}


  ###############################################

  if (nrow(df_transposed5) > 1) {
    if (df_transposed5$Hcbh[1] == 1.5 && df_transposed5$dist[1] == 1 && df_transposed5$dist[2] > 1) {
      if (df_transposed5$dptf[1] > df_transposed5$dptf[2]) {
        if (df_transposed5$Hdepth[1] < df_transposed5$Hdepth[2]) {
          df_transposed5$dptf[2] <- df_transposed5$dptf[1]
          df_transposed5$Hdepth[1] <- df_transposed5$Hdepth[2]
        }}}
  }

  ###############################################

  if (nrow(df_transposed5) == 1 && !"dptf" %in% colnames (df_transposed5)){
    df_transposed5$dptf <- df_transposed5$depth
  }
  depth<-df_transposed5$depth

  ###############################################
  df_transposed5 <- df_transposed5 %>%
    dplyr::group_by(Hdepth) %>%
    dplyr::mutate(is_duplicated = duplicated(Hdepth) | duplicated(Hdepth, fromLast = TRUE)) %>%
    # For duplicated Hcbh values, compute the max dptf and max dist
    dplyr::mutate(dptf = ifelse(is_duplicated, max(dptf, na.rm = TRUE), dptf),
           dist = ifelse(is_duplicated, max(dist, na.rm = TRUE), dist),
           Hdist = ifelse(is_duplicated,max(Hdist, na.rm = TRUE), Hdist) ) %>%
    # If Hcbh is duplicated, keep only the last row. Otherwise, keep all rows.
    dplyr::filter(ifelse(is_duplicated, row_number() == n(), TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::select(-is_duplicated)  # Drop the helper column


  ###############################################

  cols_to_transpose <- colnames(df_transposed5)
  cols_to_transpose <- setdiff(cols_to_transpose, c("Hcbh"))

  df_wide <- df_transposed5  %>%
    dplyr::select(all_of(cols_to_transpose)) %>%
    tidyr::pivot_wider(names_from = index, values_from = c(Hdepth, Hdist, depth, dist, dptf))

  df_single1 <- df_wide %>%
    dplyr::summarize(across(everything(), ~na.omit(.x)[1]))

  df_single1 <- df_single1 %>%
    dplyr::rename_with(~ gsub("_", "", .))

  df_single1 <- df_single1 %>%
    dplyr::select_if(~!any(is.na(.)))

  Hcbh_cols <- grep("^Hcbh", colnames(df4), value = TRUE)
  Hcbh_vals <- df4[, Hcbh_cols, drop = FALSE]

  all_equal <- apply(Hcbh_vals, 1, function(x) length(unique(x)) == 1)

  # If all rows have the same values across columns
  if (all(all_equal)) {
    Hcbh_vals <- Hcbh_vals[,1, drop=FALSE]
    colnames(Hcbh_vals) <- "Hcbh1"
  }


  df5b<-cbind.data.frame(df_single1,Hcbh_vals)


  ###################  rename columns

  # Exclude columns with prefixes: treeID
  cols_to_exclude <- grep("^(max_|treeID)", names(df5b), value = TRUE)
  df5b <- df5b[ , !(names(df5b) %in% cols_to_exclude)]

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

  # Rename columns Hdepth

  df5b <- df5b %>%
    dplyr::rename_with(~gsub("^Hdepth", "Hdptf", .), starts_with("Hdepth"))

  # Remove original depth columns from df5b
  depth_data <- df5b[grep("^depth", names(df5b))]
  df5b <- df5b[ , !(names(df5b) %in% names(depth_data))]

  ###############   UPDATE DPTF VALUES: HDPTF - HCBH == DPTF (REMOVING HDIST OVER THE LAST CBH) #######3

  Hcbh_col_names <- grep("^Hcbh", colnames(df5b), value = TRUE)
  dptf_col_names <- grep("^dptf", colnames(df5b), value = TRUE)
  Hdptf_col_names <- grep("^Hdptf", colnames(df5b), value = TRUE)


  # Get the positions (indices) of the dptf columns
  dptf_positions <- grep("^dptf", names(df5b))

  # Loop over dptf column positions
  for (pos in dptf_positions) {
    # Extract column names based on positions
    dptf_col_name <- names(df5b)[pos]

    # Construct the corresponding Hdptf and Hcbh column names
    suffix <- gsub("dptf", "", dptf_col_name)
    Hdptf_col_name <- paste0("Hdptf", suffix)
    Hcbh_col_name <- paste0("Hcbh", suffix)

    # Ensure the Hcbh column is unlisted
    df5b[[Hcbh_col_name]] <- unlist(df5b[[Hcbh_col_name]])

    # Check if Hdptf and Hcbh columns with the current suffix exist in the dataframe
    if (Hdptf_col_name %in% names(df5b) && Hcbh_col_name %in% names(df5b)) {
      # Check the condition and update the dptf value if necessary
      if ((df5b[[dptf_col_name]] == df5b[[Hdptf_col_name]] - df5b[[Hcbh_col_name]]) && (df5b[[Hdptf_col_name]] - df5b[[Hcbh_col_name]] == 0)) {
        df5b[[dptf_col_name]] <- 1
      }
      if (df5b[[dptf_col_name]] > df5b[[Hdptf_col_name]] - df5b[[Hcbh_col_name]]) {
        df5b[[dptf_col_name]] <- df5b[[Hdptf_col_name]] - df5b[[Hcbh_col_name]]
      }
    }
  }

  # Check if dptf == 0
  # Extract the suffixes from dptf columns

  dptf_suffixes <- gsub("dptf", "", grep("^dptf", names(df5b), value = TRUE))

  # Loop over the effdist suffixes
  for (suffix in dptf_suffixes) {
    # Construct dptf column name using the current suffix
    dptf_col_name <- paste0("dptf", suffix)

    # Check if the dptf column with the current suffix exists in the dataframe
    if (dptf_col_name %in% names(df5b)) {
      # If any value in that column is 0, update the value to 1
      df5b[df5b[[dptf_col_name]] == 0, dptf_col_name] <- 1
    }
  }


  ##############################################

  hcbh_cols <- grep("^Hcbh", colnames(df5b))
  duplicated_cols <- which(duplicated(as.list(df5b[hcbh_cols])))
  # Remove those columns
  df5b <- df5b %>% dplyr::select(-one_of(colnames(df5b)[hcbh_cols[duplicated_cols]]))

  hdist_cols <- grep("^Hdist", colnames(df5b))
  duplicated_cols1 <- which(duplicated(as.list(df5b[hdist_cols])))
  # Remove those columns
  df5b <- df5b %>% dplyr::select(-one_of(colnames(df5b)[hdist_cols[duplicated_cols1]]))

  ##########################################

  hcbh_cols <- grep("^Hcbh", colnames(df5b), value=TRUE)
  dist_cols <- grep("^dist", colnames(df5b), value=TRUE)
  dptf_cols <- grep("^dptf", colnames(df5b), value=TRUE)
  hdist_cols <- grep("^Hdist", colnames(df5b), value=TRUE)

  if(df5b$Hcbh1 == 1.5 && length(dist_cols) == length(hcbh_cols)) {
    # Get the last column name from dist_cols vector
    last_dist_column_name <- tail(dist_cols, 1)
    # Drop the column
    df5b <- df5b %>% dplyr::select(-all_of(c(last_dist_column_name)))
  }

  if(df5b$Hcbh1 > 1.5 && length(dist_cols) > length(hcbh_cols)) {
    # Get the last column name from dist_cols vector
    last_dist_column_name <- tail(dist_cols, 1)
    # Drop the column
    df5b <- df5b %>% dplyr::select(-all_of(c(last_dist_column_name)))
  }

  ##############################################
  if(df5b$Hcbh1 == 1.5 && length(hdist_cols) == length(hcbh_cols)) {
    # Get the last column name from dist_cols vector
    last_hdist_column_name <- tail(hdist_cols, 1)
    # Drop the column
    df5b <- df5b %>% dplyr::select(-all_of(c(last_hdist_column_name)))
  }

  if(df5b$Hcbh1 > 1.5 && length(hdist_cols) > length(hcbh_cols)) {
    # Get the last column name from dist_cols vector
    last_hdist_column_name <- tail(hdist_cols, 1)
    # Drop the column
    df5b <- df5b %>% dplyr::select(-all_of(c(last_hdist_column_name)))
  }
  ##############################################

  if(all(hdist_cols %in% colnames(df5b)) && all(dist_cols %in% colnames(df5b))) {
    if (length(hdist_cols) == 1 && (df5b$Hdist1 > df5b$Hdptf1) && (df5b$Hcbh1 > 1.5)) {
      df5b$Hdist1 <- NULL
      df5b$dist1 <- NULL
    }
  }

  if(all(hdist_cols %in% colnames(df5b)) && all(dist_cols %in% colnames(df5b))) {
    if (length(hdist_cols) == 1 && (df5b$dist1 == 1) && (df5b$Hcbh1 > 1.5)) {
      df5b$Hdist1 <- 1.5
    }
  }

  if(all(hdist_cols %in% colnames(df5b)) && length(hdist_cols) == 1 && all(!dist_cols %in% colnames(df5b))) {
    df5b$Hdist1 <- NULL
  }

  df5b<-get_renamed0_df(df5b)

  for (i in 1:nrow(df5b)) {
    # For each column starting with "dptf"
    dptf_cols <- grep("^dptf", colnames(df5b), value = TRUE)
    hcbh_cols <- grep("^Hcbh", colnames(df5b), value = TRUE)
    hdptf_cols <- grep("^Hdptf", colnames(df5b), value = TRUE)

    for (j in seq_along(dptf_cols)) {
      dptf_col <- dptf_cols[j]
      hcbh_col <- hcbh_cols[j]
      hdptf_col <- hdptf_cols[j]

      # Recalculate dptf values
      if (df5b[i, hdptf_col] - df5b[i, hcbh_col] != df5b[i, dptf_col]) {
        df5b[i, dptf_col] <- df5b[i, hdptf_col] - df5b[i, hcbh_col]
      }
    }
  }

  dptf_cols <- grep("^dptf", colnames(df5b), value = TRUE)
  for (j in seq_along(dptf_cols)) {
    dptf_col <- dptf_cols[j]
    if(df5b[i, dptf_col] == 0) {
      df5b[i, dptf_col]<- 1
    }}

  ####################################################
  # For each column starting with "Hdist"
  hcbh_cols <- grep("^Hcbh\\d+$", colnames(df5b), value = TRUE)
  hdist_cols <- grep("^Hdist\\d+$", colnames(df5b), value = TRUE)

  for (j in seq_along(hdist_cols)) {
    if (j < length(hcbh_cols) && df5b$Hcbh1 == 1.5) { # Ensure we're not exceeding the bounds of hcbh columns
      hdist_col <- hdist_cols[j]
      hcbh_col <- hcbh_cols[j + 1]
      df5b[1, hdist_col] <- df5b[1, hcbh_col] -  1
    }
  }

  for (j in seq_along(hdist_cols)) {
    if (j < length(hcbh_cols) && df5b$Hcbh1 > 1.5) { # Ensure we're not exceeding the bounds of hcbh columns
      hdist_col <- hdist_cols[j]
      hcbh_col <- hcbh_cols[j]
      df5b[1, hdist_col] <- df5b[1, hcbh_col] -  1
    }
  }


  treeID<-unique(factor(df2a$treeID))
  if(!"treeID" %in% colnames(df5b)) {
    df5b <- cbind(df2a[c("treeID","treeID1")],df5b )
  }

  max_height<-data.frame(df2a$max_height)
  names(max_height)<-"max_height"

  if(!"max_height" %in% colnames(df5b)) {
    df5b <- cbind(df5b,df2a[c("max_height")])
  }

  # Get original column names
  original_column_names <- colnames(df5b)

  # Specify prefixes (adjust accordingly)
  prefixes <- c("treeID", "Hcbh", "dptf","Hdptf","dist","Hdist","max_height")

  # Initialize vector to store new order
  new_order <- c()

  # Loop over prefixes
  for (prefix in prefixes) {
    # Find column names matching current prefix
    matching_columns <- grep(paste0("^", prefix), original_column_names, value = TRUE)

    # Extract numeric suffixes and order the columns based on these suffixes
    numeric_suffixes <- as.numeric(gsub(paste0("^", prefix), "", matching_columns))
    matching_columns <- matching_columns[order(numeric_suffixes)]

    # Append to new order
    new_order <- c(new_order, matching_columns)
  }

  # Reorder values by treeID
  df5b <- df5b[, new_order]
  effective_depth<-df5b

  return (effective_depth)

}
