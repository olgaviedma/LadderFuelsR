#' Compute the percentile value of each height
#' @description This function calculates the percentile value of each height
#' @usage calculate_gaps_perc2 (LAD_profiles)
#' @param LAD_profiles original tree Leaf Area Density (LAD) profile (output of [lad.profile()] function in the \emph{leafR} package.
#' An object of the class text
#' @return A data frame giving the percentile value of each height.
#' @author Olga Viedma, Carlos Silva and JM Moreno
#'
#' @examples
#' ## Not run:
#' library(magrittr)
#' library(dplyr)
#'
#' # LAD profiles derived from normalized ALS data after applying [lad.profile()] function
#' data_path <- file.path("D:/OLGA/R_library/LadderFuelsR/extdata/LAD_profiles.txt")
#' LAD_profiles <- read.table(data_path, sep = "\t", header = TRUE)
#' LAD_profiles$treeID <- factor(LAD_profiles$treeID)
#'
#' trees_name1 <- as.character(LAD_profiles$treeID)
#' trees_name2 <- factor(unique(trees_name1))
#'
#' percentile_list1<-list()
#'
#' for (i in levels(trees_name2)) {
#'   tree1 <- LAD_profiles |> dplyr::filter(treeID == i)
#'   percentiles <- calculate_gaps_perc2(tree1)
#'   percentile_list1[[i]] <- percentiles
#'}
#' gaps_perc2 <- dplyr::bind_rows(percentile_list1)
#' gaps_perc2$treeID <- factor(gaps_perc2$treeID)
#'
#' ## End(Not run)
#'
#' @export calculate_gaps_perc2
#' @importFrom dplyr select_if group_by summarise summarize mutate arrange rename rename_with filter slice slice_tail ungroup distinct
#' across matches row_number all_of vars
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
calculate_gaps_perc2 <- function(LAD_profiles) {

  df <- LAD_profiles

  treeID <- "treeID"
  print(paste("treeID:", df[[treeID]][1]))  # Debugging line

  df$height <- as.numeric(df$height)
  df$treeID <- factor(df$treeID)
  trees_name1a <- as.character(df$treeID)
  trees_name3 <- factor(unique(trees_name1a))

  lad<-df$lad

  df_ord <- df[with(df, order(lad)), ]

  all_equal <- length(unique(df_ord$lad)) == 1

  if (all_equal) {
    return(NULL)  # Skip the file, return NULL
  }

  PERCENTIL_Z <- df %>%
    dplyr::group_by(treeID) %>%
    dplyr::summarise(
      P5 = quantile(lad, probs = 0.05, na.rm = TRUE),
      P25 = quantile(lad, probs = 0.25, na.rm = TRUE),
      P50 = quantile(lad, probs = 0.50, na.rm = TRUE),
      P75 = quantile(lad, probs = 0.75, na.rm = TRUE),
      P90 = quantile(lad, probs = 0.90, na.rm = TRUE),
      P95 = quantile(lad, probs = 0.95, na.rm = TRUE),
      P99 = quantile(lad, probs = 0.99, na.rm = TRUE)
    )

  x1 <- df$height
  y1 <- df$lad

  # Identify missing and infinite values in x and y
  missing_x <- is.na(x1)
  missing_y <- is.na(y1) | is.infinite(y1)

  # Remove missing and infinite values from x and y
  x <- x1[!missing_x & !missing_y]
  y <- y1[!missing_x & !missing_y]

  fit <- smooth.spline(x, y)  # Fit a smoothing spline
  y_second_deriv <- predict(fit, fit$x, deriv = 2)  # Calculate the second derivative

  base_2drivative <- data.frame(do.call(cbind, y_second_deriv))  # convert the list into a dataframe
  base_2drivative$y <- round(base_2drivative$y, digits = 10)

  critical_points <- base_2drivative[, 2]  # Extract the values of the second derivative
  base_2drivative2 <- cbind.data.frame(df[, c(1:3)], critical_points)

  gaps_perc <- with(base_2drivative2,
                    ifelse(lad <= PERCENTIL_Z$P5 , "5",
                           ifelse(lad > PERCENTIL_Z$P5 & lad <= PERCENTIL_Z$P25, "25",
                                  ifelse(lad > PERCENTIL_Z$P25 & lad <= PERCENTIL_Z$P50, "50",
                                         ifelse(lad > PERCENTIL_Z$P50 & lad <= PERCENTIL_Z$P75, "75",
                                                ifelse(lad > PERCENTIL_Z$P75 & lad <= PERCENTIL_Z$P90, "90",
                                                       ifelse(lad > PERCENTIL_Z$P90 & lad <= PERCENTIL_Z$P95, "95",
                                                              ifelse(lad > PERCENTIL_Z$P95 & lad <= PERCENTIL_Z$P99, "99",
                                                                     ifelse(lad >  PERCENTIL_Z$P99, "100", NA)) )))))))

  gaps_perc1 <- data.frame(percentil = as.numeric(gaps_perc))
  gaps_perc2 <- cbind.data.frame(base_2drivative2, gaps_perc1)

  return(gaps_perc2)
}
