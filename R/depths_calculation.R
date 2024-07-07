# Declare global variables to avoid R CMD check warnings
utils::globalVariables(c("Hcbh", "Hdepth", "Hdepth1", "Hdist", "Hdist1",
                         "cbh", "cbh1", "depth", "depth1", "dist1",
                         "gap", "index1", "last_cbh", "last_nonzero_Hdist",
                         "max_dist", "temp_Hdepth", "temp_depth",
                         "total_depth", "total_dist"))
#' Fuels depth in meters
#' @description This function calculates fuels depth as the difference between gaps interleaved between fuel layers minus one step if the fuel depths are greater than one step.
#' @usage get_depths (LAD_profiles, distance_metrics,step= 1,min_height= 1.5, verbose=TRUE)
#' @param LAD_profiles original tree Leaf Area Density (LAD) profile (output of [lad.profile()] function in the \emph{leafR} package.
#' An object of the class text.
#' @param distance_metrics tree metrics with gaps (distances) and fuel base heights (output of [get_distance()] function).
#' An object of the class text.
#' @param step Numeric value for the actual height bin step (in meters).
#' @param min_height Numeric value for the actual minimum base height (in meters).
#' @param verbose Logical, indicating whether to display informational messages (default is TRUE).
#' @return A data frame giving fuel layers depth and the height of the depths in meters.
#' @author Olga Viedma, Carlos Silva, JM Moreno and A.T. Hudak
#'
#'@details
#' # List of tree metrics:
#' \itemize{
#' \item treeID: tree ID with strings and numeric values
#' \item treeID1: tree ID with only numeric values
#' \item cbh - Height of the fuel layers base height (m)
#' \item gap - Height of gaps between consecutive fuel layers (m)
#' \item dist: Distance between consecutive fuel layers (m)
#' \item Hdist - Height of the distance between consecutive fuel layers (m)
#' \item depth - Depth of fuel layers (m)
#' \item Hdepth - Height of the depth of fuel layers (m)
#' \item max_height - Maximum height of the tree profile
#' }
#'
#' @examples
#' library(magrittr)
#' library(dplyr)
#'
#' # Before running this example, make sure to run get_distance().
#' if (interactive()) {
#' distance_metrics <- get_distance()
#' LadderFuelsR::LAD_profiles$treeID <- factor(LadderFuelsR::LAD_profiles$treeID)
#' LadderFuelsR::distance_metrics$treeID <- factor(LadderFuelsR::distance_metrics$treeID)
#'
#' metrics_depth_list <- list()
#'
#' for (i in levels(LAD_profiles$treeID)){
#'
#' tree1 <- LAD_profiles |> dplyr::filter(treeID == i)
#' tree2 <- distance_metrics |> dplyr::filter(treeID == i)
#'
#' # Get depths for each tree
#' metrics_depth <- get_depths(tree1, tree2,step= 1,min_height= 1.5, verbose=TRUE)
#' metrics_depth_list[[i]] <- metrics_depth
#' }
#'
#' # Combine the individual data frames
#' depth_metrics <- dplyr::bind_rows(metrics_depth_list)
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
#' @seealso \code{\link{get_distance}}
#' @export
get_depths <- function (LAD_profiles,distance_metrics,step= 1,min_height= 1.5, verbose=TRUE) {

    df1 <- LAD_profiles

  df1$height<-as.numeric(df1$height)
  df1$treeID<-factor(df1$treeID)
  trees_name1a<- as.character(df1$treeID)
  trees_name3<- factor(unique(trees_name1a))


  if(min_height==0){
    min_height <-0.5

    # Ensure the column starts with a negative value
    if (df1$height[1] < min_height) {
      # Calculate the shift value
      shift_value <- abs(df1$height[1])

      # Adjust the column to start from 0
      df1$height <- df1$height + shift_value
    }


    # Ensure the column starts with a negative value
    if (df1$height[1] > min_height) {
      # Calculate the shift value
      shift_value1 <- abs(df1$height[1])

      # Adjust the column to start from 0
      df1$height <- df1$height - shift_value1
    }
  }



  if (verbose) {
    message("Unique treeIDs:", paste(unique(df1$treeID), collapse = ", "))
  }

  lad<-df1$lad

  df1_ord<-df1[with(df1, order(lad)), ]

  # Calculate percentiles for each treeID
  PERCENTIL_Z <- df1 %>%
    dplyr::group_by(treeID) %>%
    dplyr::summarise(across(lad, list(
      P5 = ~ quantile(.x, probs = 0.05, na.rm = TRUE),
      P10 = ~ quantile(.x, probs = 0.10, na.rm = TRUE),
      P15 = ~ quantile(.x, probs = 0.15, na.rm = TRUE),
      P20 = ~ quantile(.x, probs = 0.20, na.rm = TRUE),
      P25 = ~ quantile(.x, probs = 0.25, na.rm = TRUE),
      P30 = ~ quantile(.x, probs = 0.30, na.rm = TRUE),
      P35 = ~ quantile(.x, probs = 0.35, na.rm = TRUE),
      P40 = ~ quantile(.x, probs = 0.40, na.rm = TRUE),
      P45 = ~ quantile(.x, probs = 0.45, na.rm = TRUE),
      P50 = ~ quantile(.x, probs = 0.50, na.rm = TRUE),
      P55 = ~ quantile(.x, probs = 0.55, na.rm = TRUE),
      P60 = ~ quantile(.x, probs = 0.60, na.rm = TRUE),
      P65 = ~ quantile(.x, probs = 0.65, na.rm = TRUE),
      P70 = ~ quantile(.x, probs = 0.70, na.rm = TRUE),
      P75 = ~ quantile(.x, probs = 0.75, na.rm = TRUE),
      P80 = ~ quantile(.x, probs = 0.80, na.rm = TRUE),
      P85 = ~ quantile(.x, probs = 0.85, na.rm = TRUE),
      P90 = ~ quantile(.x, probs = 0.90, na.rm = TRUE),
      P95 = ~ quantile(.x, probs = 0.95, na.rm = TRUE),
      P99 = ~ quantile(.x, probs = 0.99, na.rm = TRUE)
    ))) %>%
    dplyr::ungroup()

  #print(PERCENTIL_Z)

  x1 <- df1$height
  y1 <- df1$lad

  # Identify missing and infinite values in x and y
  missing_x <- is.na(x1)
  missing_y <- is.na(y1) | is.infinite(y1)

  # Remove missing and infinite values from x and y
  x <- x1[!missing_x & !missing_y]
  y <- y1[!missing_x & !missing_y]

  # Fit a smoothing spline to the data
  fit <- smooth.spline(x, y)

  # Calculate the second derivative for all original x values
  y_second_deriv <- predict(fit, x, deriv = 2)

  # Convert the second derivative results into a data frame
  base_2drivative <- data.frame(x = y_second_deriv$x, y = round(y_second_deriv$y, digits = 10))

  # Merge the second derivative with the original data frame
  critical_points <- base_2drivative[, 2]  # Extract the values of the second derivative
  base_2drivative2 <- cbind.data.frame(df1[, c(1:3)], critical_points)

  # Calculate gaps_perc using dplyr::case_when
  gaps_perc <- df1 %>%
    dplyr::mutate(
      gaps_perc = dplyr::case_when(
        lad <= PERCENTIL_Z$lad_P5 ~ "5",
        lad > PERCENTIL_Z$lad_P5 & lad <= PERCENTIL_Z$lad_P10 ~ "10",
        lad > PERCENTIL_Z$lad_P10 & lad <= PERCENTIL_Z$lad_P15 ~ "15",
        lad > PERCENTIL_Z$lad_P15 & lad <= PERCENTIL_Z$lad_P20 ~ "20",
        lad > PERCENTIL_Z$lad_P20 & lad <= PERCENTIL_Z$lad_P25 ~ "25",
        lad > PERCENTIL_Z$lad_P25 & lad <= PERCENTIL_Z$lad_P30 ~ "30",
        lad > PERCENTIL_Z$lad_P30 & lad <= PERCENTIL_Z$lad_P35 ~ "35",
        lad > PERCENTIL_Z$lad_P35 & lad <= PERCENTIL_Z$lad_P40 ~ "40",
        lad > PERCENTIL_Z$lad_P40 & lad <= PERCENTIL_Z$lad_P45 ~ "45",
        lad > PERCENTIL_Z$lad_P45 & lad <= PERCENTIL_Z$lad_P50 ~ "50",
        lad > PERCENTIL_Z$lad_P50 & lad <= PERCENTIL_Z$lad_P55 ~ "55",
        lad > PERCENTIL_Z$lad_P55 & lad <= PERCENTIL_Z$lad_P60 ~ "60",
        lad > PERCENTIL_Z$lad_P60 & lad <= PERCENTIL_Z$lad_P65 ~ "65",
        lad > PERCENTIL_Z$lad_P65 & lad <= PERCENTIL_Z$lad_P70 ~ "70",
        lad > PERCENTIL_Z$lad_P70 & lad <= PERCENTIL_Z$lad_P75 ~ "75",
        lad > PERCENTIL_Z$lad_P75 & lad <= PERCENTIL_Z$lad_P80 ~ "80",
        lad > PERCENTIL_Z$lad_P80 & lad <= PERCENTIL_Z$lad_P85 ~ "85",
        lad > PERCENTIL_Z$lad_P85 & lad <= PERCENTIL_Z$lad_P90 ~ "90",
        lad > PERCENTIL_Z$lad_P90 & lad <= PERCENTIL_Z$lad_P95 ~ "95",
        lad > PERCENTIL_Z$lad_P95 & lad <= PERCENTIL_Z$lad_P99 ~ "99",
        lad > PERCENTIL_Z$lad_P99 ~ "100",
        TRUE ~ NA_character_
      )
    )

  # Check the distribution of gaps_perc
  #print(table(gaps_perc$gaps_perc))


  # Combine results into final data frame
  gaps_perc1 <- gaps_perc %>% mutate(percentil = as.numeric(gaps_perc))
  gaps_perc2 <- bind_cols(base_2drivative2, percentil = gaps_perc1$percentil)

  percentil <- gaps_perc2$percentil

  df <- distance_metrics

  df <- df[, !colSums(is.na(df)) > 0]
  # Select only numeric columns
  df1_numeric <- df %>% dplyr::select_if(is.numeric)


  # Assuming that columns starting with "Hdist" or "cbh" are the ones you want to keep
  columns_to_keep <- names(df1_numeric)[grepl("^dist\\d+$|^Hdist\\d+$|^cbh\\d+$", names(df1_numeric))]
  # Subset the data frame
  kk_copy <- df1_numeric[, columns_to_keep]

  # Sort the column names based on their values in the first row
  sorted_columns <- names(kk_copy)[order(unlist(kk_copy))]

  # Reorder the columns in ascending order
  kk_copy <- kk_copy[, sorted_columns]

  cbh_cols <- which(grepl("cbh", colnames(kk_copy)))
  gap_cols <- which(grepl("Hdist", colnames(kk_copy)))

  if(kk_copy$cbh1 <=min_height){
    kk_copy$cbh1<-min_height
    kk_copy$dist1<-min_height
    kk_copy$Hdist1<-min_height
  }

  # Extract the numeric suffixes from the column names
  extract_suffix <- function(name) {
    as.numeric(sub("^[^0-9]*", "", name))
  }

  # Order the columns by their numeric suffixes
  ordered_columns <- names(kk_copy)[order(sapply(names(kk_copy), extract_suffix))]

  # Reorder the dataframe columns based on the ordered suffixes
  kk_copy <- kk_copy[, ordered_columns]

  if (kk_copy$cbh1 == min_height) {
    kk_copy$dist1 <- min_height
  }

  if (kk_copy$Hdist1 == 0) {
    kk_copy$Hdist1 <-min_height
  }

  if (kk_copy$Hdist1 == min_height) {
    kk_copy$dist1 <-min_height
  }

  cols_to_transpose <- colnames(kk_copy)

  df_transposed0 <- kk_copy %>%
    dplyr::select(all_of(cols_to_transpose)) %>%
    pivot_longer(cols = everything(), names_to = c(".value", "index"), names_pattern = "(\\D+)(\\d+)?")
  index<-df_transposed0$index

  df_transposed0 <- df_transposed0 %>%
    dplyr::arrange(as.numeric(index))


  # Initialize new columns
  df_transposed0 <- df_transposed0 %>%
    mutate(depth = NA, Hdepth = NA)

#################################################333
  # Fill down the Hdist values
  df_transposed1 <- df_transposed0 %>%
    dplyr::mutate(Hdist = ifelse(Hdist == 0, NA, Hdist)) %>%  # Replace 0s with NA for fill
    fill(Hdist, .direction = "down")                  # Fill NA values downwards

  # #print the updated dataframe
  #print(df_transposed1)

    ####################################################3

  # Update the first row of each set of duplicate Hdist rows with the last cbh value
  df_updated <- df_transposed1 %>%
    dplyr::group_by(Hdist) %>%
    dplyr::mutate(last_cbh = last(cbh)) %>%
    dplyr::mutate(depth = ifelse(row_number() == 1, last_cbh, cbh)) %>%
    dplyr::select(-last_cbh)


  remove_last_row_if_equal_Hdist <- function(df) {
    # Check if dataframe has exactly 2 rows
    if (nrow(df) == 2) {
      # Check if all Hdist values are equal
      if (all(df$Hdist[-nrow(df)] == df$Hdist[nrow(df)])) {
        # Remove the last row
        df <- df[-nrow(df), ]
      }
    }
    return(df)
  }
  df_updated <- remove_last_row_if_equal_Hdist(df_updated)


  update_Hdepth_if_conditions_met <- function(df, min_height) {
    # Check if dataframe has exactly 1 row
    if (nrow(df) == 1) {
      # Check if all Hcbh values are equal to min_height
      if (df$cbh == min_height) {
        # Set Hdepth to depth
        df$Hdepth <- df$depth
      }
    }
    return(df)
  }
  df_updated <- update_Hdepth_if_conditions_met(df_updated, min_height)

  ####################################################

  if ((nrow(df_updated) == 1) && (any(!is.na(df_updated$depth)) || any(is.na(df_updated$Hdepth)))) {
    df_updated$Hdepth <- df_updated$depth
    df_updated$depth <- floor(df_updated$depth)

     }

  ####################################################

  if ((nrow(df_updated) > 1) && (all(!is.na(df_updated$depth)) || all(!is.na(df_updated$Hdepth)))) {

    df_updated <- df_updated %>%
      dplyr::group_by(Hdist) %>%
      dplyr::summarise(
        dist = max(dist, na.rm = TRUE),
        cbh = max(cbh, na.rm = TRUE),
        depth = max(depth, na.rm = TRUE),
        Hdepth = if (all(is.na(Hdepth))) NA_real_ else max(Hdepth, na.rm = TRUE)
      )


    # Loop to Hdist > Hcbh with the previous non-NA value
    for (i in 1:nrow(df_updated)) {
      if (!is.na(df_updated$Hdist[i]) && !is.na(df_updated$cbh[i]) &&  !is.na(df_updated$depth[i]) && df_updated$depth[i] > df_updated$cbh[i]) {
      df_updated$depth[i] <- round(df_updated$cbh [i] - df_updated$Hdist[i], 0)
      df_updated$Hdepth[i] <-df_updated$cbh [i] + df_updated$depth[i]
      }}

      for (i in 1:nrow(df_updated)) {
        if (!is.na(df_updated$Hdist[i]) && !is.na(df_updated$cbh[i]) &&  !is.na(df_updated$depth[i]) && df_updated$depth[i] <= df_updated$cbh[i]) {
          df_updated$depth[i] <- round(df_updated$cbh [i] - df_updated$Hdist[i], 0)
          df_updated$Hdepth[i] <-df_updated$cbh [i]
        }}

    if (nrow(df_updated) > 1) {
    # Loop to Hdist > Hcbh with the previous non-NA value
    for (i in 1:(nrow(df_updated) - 1)) {
      if (!is.na(df_updated$Hdist[i]) && !is.na(df_updated$cbh[i]) && df_updated$Hdist[i] > df_updated$cbh[i]) {
        df_updated$depth[i] <- round(df_updated$Hdist[i + 1] - df_updated$cbh[i] - df_updated$dist[i + 1], 0)
        df_updated$Hdepth[i] <- df_updated$Hdist[i + 1] - df_updated$dist[i + 1]
      }
    }

    # Replace last Hdist values with max_height if they are NA or 0
    last_row <- nrow(df_updated)
    if (is.na(df_updated$depth[last_row]) || df_updated$depth[last_row] == 0) {
      df_updated$depth[last_row] <- round((df_updated$cbh[last_row] - df_updated$Hdist[last_row]) - step, 0)
      df_updated$Hdepth[last_row] <- df_updated$cbh[last_row]
    }
    }


    # Select the first row of each set of duplicated Hdist rows where dist > 0
    first_rows <- df_transposed0 %>%
      dplyr::filter(dist > 0) %>%         # Keep rows where dist > 0
      dplyr::group_by(Hdist) %>%          # Group by Hdist
      slice(1) %>%                 # Select the first row of each group
      dplyr::ungroup()                    # Remove grouping


    cbh_subset <- first_rows %>%
      dplyr::select(cbh)


    df_updated_nocbh <- df_updated %>%
      dplyr::select(-cbh)

    # Check length of df_transposed2_nocbh and cbh_subset
    if (nrow(df_updated_nocbh) > nrow(cbh_subset)) {
      # Trim df_transposed2_nocbh from the last row
      df_updated_nocbh <- df_updated_nocbh %>%
        slice(1:nrow(cbh_subset))
    }

    df_updated<-cbind(df_updated_nocbh,cbh_subset)

    # Loop to Hdist > Hcbh with the previous non-NA value
    for (i in 1:nrow(df_updated)) {
      if (df_updated$depth [i] == 0) {
        df_updated$depth[i] <- 1
      }}

    df_transposed3<- df_updated

    if (nrow(df_transposed3) > 1) {

      if ((df_transposed3$depth[1] != round(df_transposed3$Hdist[1 + 1] -  df_transposed3$cbh [1] -df_transposed3$dist [1 + 1], 0))) {
    # Loop to Hdist > Hcbh with the previous non-NA value{
      df_transposed3$depth[1] <- round(df_transposed3$Hdist[1 + 1] -  df_transposed3$cbh [1] -df_transposed3$dist [1 + 1], 0)
      }
    }


    # Apply the condition to the first row
    if (df_transposed3$cbh[1] == min_height && df_transposed3$Hdepth[1] < df_transposed3$depth[1]) {
      df_transposed3$Hdepth[1] <- df_transposed3$cbh[1] + df_transposed3$depth[1]
    }

   ###################

    if (nrow(df_transposed3) > 1) {

      # Iterate over rows to apply the condition and update
      for (i in 1:(nrow(df_transposed3) - 1)) {
        if (!is.na(df_transposed3$Hdist[i]) && !is.na(df_transposed3$cbh[i]) && df_transposed3$Hdist[i] > df_transposed3$cbh[i]) {
        # Calculate the expected depth and Hdepth values
        expected_depth <- round(df_transposed3$Hdist[i + 1] - df_transposed3$cbh[i] - df_transposed3$dist[i + 1], 0)
        expected_Hdepth <- df_transposed3$Hdist[i + 1] - df_transposed3$dist[i + 1]

        # Check and update the depth and Hdepth values if they do not match the expected values
        if (df_transposed3$depth[i] != expected_depth) {
          df_transposed3$depth[i] <- expected_depth
        }
        if (df_transposed3$Hdepth[i] != expected_Hdepth) {
          df_transposed3$Hdepth[i] <- expected_Hdepth
        }
      }
      } }

    # Update depth based on the condition
    df_transposed3 <- df_transposed3 %>%
      mutate(depth = ifelse(dist > 0, Hdepth - cbh, depth))


    # Loop to Hdist > Hcbh with the previous non-NA value
    for (i in 1:nrow(df_transposed3)) {
      if (df_transposed3$depth [i] == 0) {
        df_transposed3$depth[i] <- 1
      }}


       ######################################

    cols_to_transpose <- grep("^Hdepth", colnames(df_transposed3), value = TRUE)

    df_transposed3$index <- 1:nrow(df_transposed3)

    # Pivot wider
    df_wide <- df_transposed3 %>%
      pivot_wider(names_from = index, values_from = c(all_of(cols_to_transpose), dist, Hdist,depth, cbh))


    df_single1 <- df_wide %>%
      dplyr::summarize(across(everything(), ~na.omit(.x)[1]))

    df_single1  <- df_single1 %>%
      dplyr::select_if(~!any(is.na(.)))

    depth_data<-cbind.data.frame(df_single1)

    # Remove underscores from column names
    depth_data <- depth_data %>%
      rename_with(~ str_remove_all(., "_"), everything())

    treeID<-unique(factor(df$treeID))
    treeID1<-unique(factor(df$treeID1))

    max_height<-data.frame(df$max_height)
    names(max_height)="max_height"

    depth_metrics <- cbind.data.frame(treeID, treeID1, depth_data,max_height)

  }

  ##########################################################3

  if (!exists("depth_metrics")){

  if (any(is.na(df_updated$depth)) || any(is.na(df_updated$Hdepth))) {

    # Group by Hdist and summarize to get max dist and max cbh
  df_updated <- df_updated %>%
    group_by(Hdist) %>%
    summarise(
      dist = max(dist),
      cbh = max(cbh)
    )

  # Update cbh with Hdist where Hdist = 1.5
  df_transposed2 <- df_updated %>%
    mutate(cbh = ifelse(Hdist == min_height, Hdist, cbh))

  # Initialize new columns
  df_transposed2 <- df_transposed2 %>%
    mutate(depth = NA, Hdepth = NA)

  ####################################################3

  # Loop to Hdist > Hcbh with the previous non-NA value
  for (i in 1:nrow(df_transposed2)) {
    df_transposed2$depth[i] <- round(df_transposed2$Hdist[i + 1] -  df_transposed2$cbh [i] -df_transposed2$dist [i + 1], 0)
    df_transposed2$Hdepth[i] <- df_transposed2$Hdist[i + 1] - df_transposed2$dist [i + 1]
    }


  # Replace last Hdist values with max_height if they are NA or 0
  last_row <- nrow(df_transposed2)
  if (is.na(df_transposed2$depth [last_row]) || df_transposed2$depth[last_row] == 0) {
    df_transposed2$depth[last_row] <-  df_transposed2$cbh[last_row] - df_transposed2$Hdist[last_row]
    df_transposed2$Hdepth[last_row] <-  df_transposed2$cbh[last_row]
  }

  # Loop to Hdist > Hcbh with the previous non-NA value
  for (i in 1:nrow(df_transposed2)) {
    if (df_transposed2$depth [i] == 0) {
    df_transposed2$depth[i] <- 1
  }}

  # Select the first row of each set of duplicated Hdist rows where dist > 0
  first_rows <- df_transposed0 %>%
    dplyr::filter(dist > 0) %>%         # Keep rows where dist > 0
    group_by(Hdist) %>%          # Group by Hdist
    slice(1) %>%                 # Select the first row of each group
    ungroup()                    # Remove grouping


  cbh_subset <- first_rows %>%
    dplyr::select(cbh)


  df_transposed2_nocbh <- df_transposed2 %>%
    dplyr::select(-cbh)

  # Check length of df_transposed2_nocbh and cbh_subset
  if (nrow(df_transposed2_nocbh) > nrow(cbh_subset)) {
    # Trim df_transposed2_nocbh from the last row
    df_transposed2_nocbh <- df_transposed2_nocbh %>%
      slice(1:nrow(cbh_subset))
  }

  df_transposed3<-cbind(df_transposed2_nocbh,cbh_subset)

  # Update depth based on the condition
  df_transposed3 <- df_transposed3 %>%
    mutate(depth = ifelse(dist > 0, Hdepth - cbh, depth))

  # Update depth = 0
  df_transposed3 <- df_transposed3 %>%
    mutate(depth = ifelse(depth == 0, 1, depth))

  ###############################################

  cols_to_transpose <- grep("^Hdepth", colnames(df_transposed3), value = TRUE)

  df_transposed3$index <- 1:nrow(df_transposed3)

  # Pivot wider
  df_wide <- df_transposed3 %>%
    pivot_wider(names_from = index, values_from = c(all_of(cols_to_transpose), dist, Hdist,depth, cbh))


  df_single1 <- df_wide %>%
    dplyr::summarize(across(everything(), ~na.omit(.x)[1]))

  df_single1  <- df_single1 %>%
    dplyr::select_if(~!any(is.na(.)))

  depth_data<-cbind.data.frame(df_single1)

  # Remove underscores from column names
  depth_data <- depth_data %>%
    rename_with(~ str_remove_all(., "_"), everything())

  treeID<-unique(factor(df$treeID))
  treeID1<-unique(factor(df$treeID1))

  max_height<-data.frame(df$max_height)
  names(max_height)="max_height"

  depth_metrics <- cbind.data.frame(treeID, treeID1, depth_data,max_height)


  } else {

    add_suffix_and_convert <- function(df, suffix = 1) {
      # Add suffix to each column name
      colnames(df) <- paste0(colnames(df), suffix)

      # Convert tibble to dataframe
      df <- as.data.frame(df)

      return(df)
    }

    df_updated <- add_suffix_and_convert(df_updated)

    df_updated$depth1 <- floor(df_updated$depth1)

    df_updated <- df_updated %>%
      dplyr::rename(index = index1)

     # Update depth = 0
     df_updated <- df_updated %>%
       mutate(depth1 = ifelse(depth1 == 0, 1, depth1))

    df_transposed3<-df_updated


  # Define columns to transpose based on the existing columns
  cols_to_transpose <- grep("^Hdepth", colnames(df_transposed3), value = TRUE)

  # Ensure the index column is present and correctly named
  df_transposed3$index <- 1:nrow(df_transposed3)

  # Correctly pivot the dataframe
  df_wide <- df_transposed3 %>%
    pivot_wider(names_from = index, values_from = c(dist1, cbh1, Hdist1, depth1, Hdepth1))


  df_single1 <- df_wide %>%
    dplyr::summarize(across(everything(), ~na.omit(.x)[1]))

  df_single1  <- df_single1 %>%
    dplyr::select_if(~!any(is.na(.)))

  depth_data<-cbind.data.frame(df_single1)

  # Remove underscores from column names
  depth_data <- depth_data %>%
    rename_with(~ str_remove_all(., "_1"), everything())

  treeID<-unique(factor(df$treeID))
  treeID1<-unique(factor(df$treeID1))

  max_height<-data.frame(df$max_height)
  names(max_height)="max_height"

  depth_metrics <- cbind.data.frame(treeID, treeID1, depth_data,max_height)
  }
  }

  return(depth_metrics)
}
