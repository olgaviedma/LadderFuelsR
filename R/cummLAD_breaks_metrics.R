#' CBH estimation using the breaking point method and the LAD percentage below and above the CBH
#' @description This function calculates the canopy base height (CBH) of the vertical tree profile (VTP) using a segmented regression model fitted
#' to the cumulative LAD values as a function of height.The function also calculates the percentage of LAD values below and above the identified
#' CBH or breaking point.
#' @usage get_cum_break(LAD_profiles, cbh_metrics, threshold=75, verbose=TRUE)
#' @param LAD_profiles Original tree Leaf Area Density (LAD) profile (output of [lad.profile()] function in the \emph{leafR} package).
#' An object of the class data frame.
#' @param cbh_metrics CBH metrics based on three criteria: maximum LAD percentage, maximum distance and last distance (output of [get_cbh_metrics()] function).
#' An object of the class data frame.
#' @param threshold Numeric value of the LAD percentage below or above the breaking point to set the CBH (default 75).
#' @param verbose Logical, indicating whether to display informational messages (default is TRUE).
#' @return A data frame identifying the CBH of the vertical tree profile (VTP) based on the breaking point method and
#' the percentage of LAD values below and above the identified CBH or breaking point.
#' @author Olga Viedma, Carlos Silva and JM Moreno
#' @details
#' # List of tree metrics:
#' \itemize{
#' \item treeID: tree ID with strings and numeric values
#' \item treeID1: tree ID with only numeric values
#' \item Hdist: Height of the distance between the ground and the CBH or breaking point (m)
#' \item Hcbh_brpt: Height of the CBH based on the breaking point method (m)
#' \item below_hcbhbp: Percentage of LAD values below the CBH or breaking point
#' \item above_hcbhbp: Percentage of LAD values above the CBH or breaking point
#' \item bp_hcbh: Height of the CBH based on the breaking point method or on the maximum LAD criterium if there is not breaking point (m)
#' \item bp_Hdptf: Height of the canopy layer depth using the breaking point method or the maximum LAD criterium (m)
#' \item bp_dptf: Depth of the CBH using the breaking point method or the maximum LAD criterium (m)
#' \item bp_Hdist: Height of the distance between the CBH and the ground using the breaking point method or the maximum LAD criterium (m)
#' \item bp_effdist: Distance between the CBH and the ground using the breaking point method or the maximum LAD criterium (m)
#' \item bp_lad: Percentage of LAD comprised by the canopy layer
#' \item cumlad: Cumulative LAD values at the CBH or breaking point
#' \item max_height: Maximum height of the tree profile
#' }
#' @examples
#' library(magrittr)
#' library(segmented)
#' library(gdata)
#' library(dplyr)
#'
#' # LAD profiles derived from normalized ALS data after applying [lad.profile()] function
#' LAD_profiles <- read.table(system.file("extdata", "LAD_profiles.txt", package = "LadderFuelsR"),
#' header = TRUE)
#' LAD_profiles$treeID <- factor(LAD_profiles$treeID)
#'
#' # Before running this example, make sure to run get_cbh_metrics().
#' if (interactive()) {
#' cbh_metrics <- get_cbh_dist()
#' LadderFuelsR::cbh_metrics$treeID <- factor(LadderFuelsR::cbh_metrics$treeID)
#'
#' trees_name1 <- as.character(cbh_metrics$treeID)
#' trees_name2 <- factor(unique(trees_name1))
#'
#' cum_LAD_metrics_list <- list()
#'
#' for (i in levels(trees_name2)) {
#' # Filter data for each tree
#' tree1 <- LAD_profiles |> dplyr::filter(treeID == i)
#' tree2 <- cbh_metrics |> dplyr::filter(treeID == i)
#'
#' # Get cumulative LAD metrics for each tree
#' cum_LAD_metrics <- get_cum_break(tree1, tree2, threshold=75, verbose=TRUE)
#' cum_LAD_metrics_list[[i]] <- cum_LAD_metrics
#' }
#'
#' # Combine the individual data frames
#' cummulative_LAD <- dplyr::bind_rows(cum_LAD_metrics_list)
#' }
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
#' @seealso \code{\link{get_cbh_metrics}}
#' @export
get_cum_break <- function(LAD_profiles, cbh_metrics, threshold=75, verbose=TRUE) {

  # Initialize the result data frame
  closest_row <- data.frame()

  df_effective<- cbh_metrics
  df_orig<- LAD_profiles

  #print(paste(unique(df_effective$treeID), collapse = ", "))
  if (verbose) {message("Unique treeIDs:", paste(unique(df_effective$treeID), collapse = ", "))}

  tryCatch({
    # Validate input data frames
    stopifnot(is.data.frame(df_orig), is.data.frame(df_effective))
    stopifnot("height" %in% colnames(df_orig), "lad" %in% colnames(df_orig))

    # Remove columns with only NAs from df_effective
    df_effective <- df_effective[, colSums(!is.na(df_effective)) > 0]

    treeID<-"treeID"
    treeID1<-"treeID1"

    # Add more noise to x-values

    df_orig$height <- jitter(df_orig$height, amount = 1)

    # Calculate the cumulative sum
    df_orig$cumulative_value <- cumsum(df_orig$lad)

    # Fit a segmented regression with automatic breakpoint estimation
    seg.mod <- segmented(lm(cumulative_value ~ height, data = df_orig),
                         seg.Z = ~height,
                         control = seg.control(fix.npsi = TRUE))

    # Check if breakpoint estimation was successful
    if (is.null(seg.mod$psi) || any(is.na(seg.mod$coef))) {
      # If breakpoint estimation failed, try alternative modeling approach (e.g., polynomial regression)
      poly.mod <- lm(cumulative_value ~ poly(height, degree = 2), data = df_orig)
      breakpoint <- poly.mod$coefficients[3]  # Adjusted to get the appropriate coefficient
    } else {
      # Get the breakpoint value
      breakpoint <- seg.mod$psi["psi1.height", "Est."]
          }

    # Get the index of the value in df$height that is closest to the breakpoint
    index <- which.min(abs(df_orig$height - breakpoint))

    # Handle breakpoints at the boundary
    if (index <= 1 || index >= nrow(df_orig)) {
      warning("Breakpoint at the boundary detected. Consider an alternative modeling approach or data transformation.")
      closest_row <- data.frame(Hcbh_bp = NA, cumlad = NA)
    } else {
      # Subset the row from the data frame that is closest to the breakpoint
      closest_row <- df_orig[index, c("height", "cumulative_value")]
      names(closest_row) <- c("Hcbh_brpt", "cumlad")
    }


  }, error = function(e) {
    # Error handling code block
    error_row <- df_orig[is.na(df_orig$cumulative_value), ]
    message("An error occurred in row ", which(is.na(df_orig$cumulative_value)),
            " with treeID '", error_row$treeID, "': ", conditionMessage(e))
    closest_row <- df_effective
  })


  if(length(closest_row) > 0 && all(!is.na(closest_row))) {

    total_lad <- sum(df_orig$lad)

    Hcbh_brpt <- closest_row[1, grep("^Hcbh_brpt$", names(closest_row))]

    below_df1 <- df_orig[df_orig$height <= Hcbh_brpt, ]
    above_df1 <- df_orig[df_orig$height >= Hcbh_brpt, ]

    lad_below <- below_df1$lad
    lad_above <- above_df1$lad

    below_percentage <- sum(lad_below) / total_lad * 100
    above_percentage <- sum(lad_above) / total_lad * 100

    # Store the percentages below and above Hcbh_bp
    percentage2 <- data.frame(below_percentage, above_percentage)

    # Create the column name by concatenating Hcbh and Hdepth column names
    names(percentage2) <- c("below_hcbhbp", "above_hcbhbp")

    output_df2<-data.frame(closest_row,percentage2)

    lad_columns <- grep("^Hcbh_brpt$", names(output_df2))
    lad_columns1 <- names(output_df2)[lad_columns]

    # Check how many lad_columns in the current row have numeric values
    numeric_count <- sum(sapply(output_df2[1, lad_columns1], function(x) is.numeric(x) & !is.na(x)))
  }

    threshold_value <- threshold

    if (exists("output_df2") && length(output_df2) > 0) {

    # Check the specified conditions
    if(!is.na(output_df2$above_hcbhbp[1]) && (output_df2$above_hcbhbp[1] >= threshold_value) && numeric_count == 1) {
      output_df2$bp_Hcbh[1] <- round(output_df2$Hcbh_brpt[1], 1)
      output_df2$bp_Hdptf[1] <- round(df_effective$max_height[1], 1)
      output_df2$bp_Hdist[1] <- round(output_df2$Hcbh_brpt[1] - 1, 1)
      output_df2$bp_effdist[1] <- round(output_df2$Hcbh_brpt[1] - 1, 0)
      output_df2$bp_dptf[1] <- round(df_effective$max_height[1] - output_df2$Hcbh_brpt[1], 0)
      output_df2$bp_lad[1] <- output_df2$above_hcbhbp[1]

    }

    # Check the specified conditions and modify columns if necessary
    if (!is.na(output_df2$below_hcbhbp[1]) && output_df2$below_hcbhbp[1] >= threshold_value && numeric_count == 1) {
      output_df2$bp_Hcbh[1] <- 1.5
      output_df2$bp_Hdptf[1] <- round(output_df2$Hcbh_brpt[1], 1)
      output_df2$bp_Hdist[1] <- 0.5
      output_df2$bp_effdist[1] <- 0
      output_df2$bp_dptf[1] <- round(output_df2$Hcbh_brpt[1], 0)
      output_df2$bp_lad[1] <- round(output_df2$below_hcbhbp[1], 1)
    }


    if(!is.na(output_df2$below_hcbhbp[1]) && !is.na(output_df2$above_hcbhbp[1]) && (output_df2$below_hcbhbp[1]< threshold_value) &&
       (output_df2$above_hcbhbp[1]< threshold_value) && numeric_count == 1) {

       # Check if df_effective$maxlad_Hcbh has a valid value
      if (!is.na(df_effective$maxlad_Hcbh[1])) {
        # Modify values in output_df2 based on conditions
        output_df2$bp_Hcbh[1] <- round(df_effective$maxlad_Hcbh[1], 1)
        output_df2$bp_Hdptf[1] <- round(df_effective$maxlad_Hdptf[1], 1)
        output_df2$bp_Hdist[1] <- round(df_effective$maxlad_Hdist[1], 1)
        output_df2$bp_effdist[1] <- round(df_effective$maxlad_effdist[1], 0)
        output_df2$bp_dptf[1] <- round(df_effective$maxlad_dptf[1], 0)
        output_df2$bp_lad[1] <- round(df_effective$maxlad_lad[1], 1)
      }
    }


  } else {

      # Initialize output_df2 with at least one row and necessary columns
      output_df2 <- data.frame(
        bp_Hcbh = numeric(1),
        bp_Hdptf = numeric(1),
        bp_Hdist = numeric(1),
        bp_effdist = numeric(1),
        bp_dptf = numeric(1),
        bp_lad = numeric(1),
        stringsAsFactors = FALSE
      )

    output_df2$bp_Hcbh[1] <-  round(df_effective$maxlad_Hcbh[1], 1)
    output_df2$bp_Hdptf[1] <- round(df_effective$maxlad_Hdptf[1], 1)
    output_df2$bp_Hdist[1] <- round(df_effective$maxlad_Hdist[1], 1)
    output_df2$bp_effdist[1] <- round(df_effective$maxlad_effdist[1], 0)
    output_df2$bp_dptf[1] <- round(df_effective$maxlad_dptf[1], 0)
    output_df2$bp_lad[1] <- round(df_effective$maxlad_lad[1], 1)
  }

  treeID<-unique(factor(df_effective$treeID))
  if(!"treeID" %in% colnames(output_df2)) {
    output_df2 <- cbind(df_effective[c("treeID","treeID1")],output_df2 )
  }

  max_height<-data.frame(df_effective$max_height)
  names(max_height)<-"max_height"

  if(!"max_height" %in% colnames(output_df2)) {
    output_df2 <- cbind(output_df2, df_effective[c("max_height")])
  }

  cummulative_LAD <-output_df2


  return(cummulative_LAD)
}
