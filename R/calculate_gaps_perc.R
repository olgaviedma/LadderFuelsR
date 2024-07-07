#' Compute the percentile value of each height
#' @description This function calculates the percentile value of each height
#' @usage calculate_gaps_perc (LAD_profiles,min_height=1.5)
#' @param LAD_profiles original tree Leaf Area Density (LAD) profile (output of [lad.profile()] function in the \emph{leafR} package.
#' An object of the class text
#' @param min_height Numeric value for the actual minimum base height (in meters).
#' @return A data frame giving the percentile value of each height.
#' @author Olga Viedma, Carlos Silva, JM Moreno and A.T. Hudak
#'
#' @examples
#' library(magrittr)
#' library(dplyr)
#'
#' # LAD profiles derived from normalized ALS data after applying [lad.profile()] function
#' LAD_profiles <- read.table(system.file("extdata", "LAD_profiles.txt", package = "LadderFuelsR"),
#' header = TRUE)
#' LAD_profiles$treeID <- factor(LAD_profiles$treeID)
#'
#' trees_name1 <- as.character(LAD_profiles$treeID)
#' trees_name2 <- factor(unique(trees_name1))
#'
#' percentile_list1<-list()
#'
#' for (i in levels(trees_name2)) {
#' tree1 <- LAD_profiles |> dplyr::filter(treeID == i)
#' percentiles <- calculate_gaps_perc (tree1,min_height=1.5)
#' percentile_list1[[i]] <- percentiles
#' }
#' gaps_perc <- dplyr::bind_rows(percentile_list1)
#' gaps_perc$treeID <- factor(gaps_perc$treeID)
#'
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
#' @export
calculate_gaps_perc <- function(LAD_profiles, min_height=1.5) {

  df <- LAD_profiles

  if(min_height==0){
    min_height <-0.5

    # Ensure the column starts with a negative value
    if (df$height[1] < min_height) {
      # Calculate the shift value
      shift_value <- abs(df$height[1])

      # Adjust the column to start from 0
      df$height <- df$height + shift_value
    }


    # Ensure the column starts with a negative value
    if (df$height[1] > min_height) {
      # Calculate the shift value
      shift_value1 <- abs(df$height[1])

      # Adjust the column to start from 0
      df$height <- df$height - shift_value1
    }
  }


  treeID <- "treeID"
  print(paste("treeID:", df[[treeID]][1]))  # Debugging line
  df$height <- as.numeric(df$height)
  df$treeID <- factor(df$treeID)
  trees_name1a <- as.character(df$treeID)
  trees_name3 <- factor(unique(trees_name1a))

  lad<-df$lad
  df_ord<-df[with(df, order(lad)), ]

  all_equal <- length(unique(df_ord$lad)) == 1

  if(all_equal) {

       gaps_perc2<-df_ord

  } else {


    # Calculate percentiles for each treeID
    PERCENTIL_Z <- df %>%
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

    x1 <- df$height
    y1 <- df$lad

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
    base_2drivative2 <- cbind.data.frame(df[, c(1:3)], critical_points)

    # Calculate gaps_perc using dplyr::case_when
    gaps_perc <- df %>%
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
    gaps_perc1 <- gaps_perc %>% dplyr::mutate(percentil = as.numeric(gaps_perc))
    gaps_perc2 <- bind_cols(base_2drivative2, percentil = gaps_perc1$percentil)

    }
  return(gaps_perc2)
}
