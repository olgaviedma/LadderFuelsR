#'
#' CBH and the LAD percentage below and above the CBH using the breaking point method
#'
#' @description
#' This function calculates the canopy base height (CBH) of the vertical tree profile (VTP) using a segmented regression model fitted to
#' the cumulative LAD values as a function of height.The function also calculates the percentage of LAD values below and above the identified CBH or breaking point.
#'
#' @usage
#' get_cum_break(LAD_profiles, effective_distances)
#'
#' @param LAD_profiles
#' Original tree Leaf Area Index (LAD) profile (output of [lad.profile()] function in the \emph{leafR} package).
#' An object of the class data frame.
#'
#' @param effective_distances
#' Tree metrics of fuel layers separated by distances greater than 1 m (output of [get_effective_gap()] function).
#' An object of the class data frame.
#'
#' @return
#' A data frame identifying the CBH of the vertical tree profile (VTP) based on the breaking point identified by the segmented regression model,
#' and the percentage of LAD values below and above the identified CBH or breaking point.
#'
#' @author
#' Olga Viedma, Carlos Silva and JM Moreno
#'
#' @details
#' \itemize{
#'   \item treeID: tree ID with strings and numeric values
#'   \item treeID1: tree ID with only numeric values
#'   \item Hdist: Height of the distance between the ground and the CBH or breaking point (m)
#'   \item Hcbh: Height of the CBH or breaking point (m)
#'   \item Hcbh_bp: Height of the CBH or breaking point (m)
#'   \item effdist: Distance between the ground and the CBH or breaking point (m)
#'   \item dptf: Depth of the CBH or breaking point (m)
#'   \item Hdptf: Height of the depth of the CBH or breaking point (m)
#'   \item below_hcbhbp: Percentage of LAD values below the CBH or breaking point
#'   \item above_hcbhbp: Percentage of LAD values above the CBH or breaking point
#'   \item maxlad_Hcbh: Height of the CBH or breaking point with maximum LAD percentage
#'   \item max_Hcbh: Height of the CBH or breaking point at maximum height distance
#'   \item last_Hcbh: Height of the CBH or breaking point at the last distance
#'   \item maxlad_: Values of distance and fuel depth and their corresponding heights for the CBH or breaking point with maximum LAD percentage
#'   \item max_: Values of distance and fuel depth and their corresponding heights for the CBH or breaking point at maximum height distance
#'   \item last_: Values of distance and fuel depth and their corresponding heights for the CBH or breaking point at the last distance
#'   \item cumlad: Cumulative LAD values at the CBH or breaking point
#'   \item max_height: Maximum height of the tree profile
#' }
#'
#' @examples
#' ## Not run:
#' library(magrittr)
#' library(segmented)
#' library(gdata)
#' library(dplyr)
#'
#' # LAD profiles derived from normalized ALS data after applying [lad.profile()] function
#' LAD_profiles$treeID <- factor(LAD_profiles$treeID)
#'
#' # Load the effective_distances object
#' if (interactive()) {
#'   effective_distances <- get_effective_gap()
#'   LadderFuelsR::effective_distances$treeID <- factor(LadderFuelsR::effective_distances$treeID)
#'
#'   trees_name1 <- as.character(effective_distances$treeID)
#'   trees_name2 <- factor(unique(trees_name1))
#'
#'   cum_LAD_metrics_list <- list()
#'
#'   for (i in levels(trees_name2)) {
#'     # Filter data for each tree
#'     tree1 <- LAD_profiles |> dplyr::filter(treeID == i)
#'     tree2 <- effective_distances |> dplyr::filter(treeID == i)
#'
#'     # Get cumulative LAD metrics for each tree
#'     cum_LAD_metrics <- get_cum_break(tree1, tree2)
#'     cum_LAD_metrics_list[[i]] <- cum_LAD_metrics
#'   }
#'
#'   # Combine the individual data frames
#'   cummulative_LAD <- dplyr::bind_rows(cum_LAD_metrics_list)
#' }
#' ## End(Not run)
#'
#' @export get_cum_break
#' @importFrom dplyr select_if group_by summarise summarize mutate arrange rename rename_with filter slice ungroup distinct
#' @importFrom magrittr %>%
#' @importFrom segmented segmented seg.control
#' @importFrom SSBtools RbindAll
#' @importFrom gdata startsWith
#' @include gap_fbh.R
#' @include distances_calculation.R
#' @include depths_calculation.R
#' @include corrected_base_heights.R
#' @include corrected_depth.R
#' @include corrected_distances.R
get_cum_break <- function(LAD_profiles, effective_distances) {

  # Initialize the result data frame
  closest_row <- data.frame()

  df_effective<- effective_distances
  df_orig<- LAD_profiles

  tryCatch({
    # Validate input data frames
    stopifnot(is.data.frame(df_orig), is.data.frame(df_effective))
    stopifnot("height" %in% colnames(df_orig), "lad" %in% colnames(df_orig))

    # Remove columns with only NAs from df_effective
    df_effective <- df_effective[, colSums(!is.na(df_effective)) > 0]

    treeID<-"treeID"
    treeID1<-"treeID1"
    print(paste("treeID:", df_effective[[treeID]]))  # Debugging line

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
      names(closest_row) <- c("Hcbh_bp", "cumlad")
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

    hcbh_bp <- closest_row[1, grep("^Hcbh_bp$", names(closest_row))]

    below_df1 <- df_orig[df_orig$height <= hcbh_bp, ]
    above_df1 <- df_orig[df_orig$height >= hcbh_bp, ]

    lad_below <- below_df1$lad
    lad_above <- above_df1$lad

    below_percentage <- sum(lad_below) / total_lad * 100
    above_percentage <- sum(lad_above) / total_lad * 100

    # Store the percentages below and above Hcbh_bp
    percentage2 <- data.frame(below_percentage, above_percentage)

    # Create the column name by concatenating Hcbh and Hdepth column names
    names(percentage2) <- c("below_hcbhbp", "above_hcbhbp")

    output_df2<-data.frame(closest_row,percentage2)

    lad_columns <- grep("^Hcbh_bp$", names(output_df2))
    lad_columns1 <- names(output_df2)[lad_columns]

    # Check how many lad_columns in the current row have numeric values
    numeric_count <- sum(sapply(output_df2[1, lad_columns1], function(x) is.numeric(x) & !is.na(x)))

    # Check the specified conditions
    if(!is.na(output_df2$above_hcbhbp[1]) && output_df2$above_hcbhbp[1] > 75 && numeric_count == 1) {
      output_df2$maxlad_Hcbh[1] <- round(output_df2$Hcbh_bp[1], 1)
      output_df2$maxlad_Hdptf[1] <- round(df_effective$max_height[1], 1)
      output_df2$maxlad_Hdist[1] <- round(output_df2$Hcbh_bp[1] - 1, 1)
      output_df2$maxlad_effdist[1] <- round(output_df2$Hcbh_bp[1] - 1, 0)

      output_df2$Hcbh1_Hdptf1[1] <- output_df2$above_hcbhbp[1]
      output_df2$Hcbh1[1] <- output_df2$Hcbh_bp[1]
      output_df2$Hdptf1[1] <- df_effective$max_height[1]
      output_df2$dptf1[1] <- round(df_effective$max_height[1] - output_df2$Hcbh_bp[1], 0)
      output_df2$Hdist1[1] <- round(output_df2$Hcbh_bp[1] - 1, 1)
      output_df2$effdist1[1] <- round(output_df2$Hcbh_bp[1] - 1, 0)

      output_df2$last_effdist[1] <- round(output_df2$Hcbh_bp[1] - 1, 0)
      output_df2$last_Hcbh[1] <- round(output_df2$Hcbh_bp[1], 1)
      output_df2$last_Hdptf[1] <- round(df_effective$max_height[1], 1)
      output_df2$last_Hdist[1] <- round(output_df2$Hcbh_bp[1] - 1, 1)

      output_df2$max_dptf[1] <- round(df_effective$max_height[1] - output_df2$Hcbh_bp[1], 0)
      output_df2$max_Hcbh[1] <- round(output_df2$Hcbh_bp[1], 1)
      output_df2$max_Hdptf[1] <- round(df_effective$max_height[1], 1)
      output_df2$max_Hdist[1] <- round(output_df2$Hcbh_bp[1] - 1, 1)
    }



    # Check the specified conditions
    if(!is.na(output_df2$below_hcbhbp[1]) && output_df2$below_hcbhbp[1] > 75 && numeric_count == 1) {
      output_df2$maxlad_Hcbh[1] <- 1.5
      output_df2$maxlad_Hdptf[1] <- round(output_df2$Hcbh_bp[1], 1)
      output_df2$maxlad_Hdist[1] <- 0.5
      output_df2$maxlad_effdist[1] <- NA
      output_df2$Hcbh1_Hdptf1[1] <- round(output_df2$below_hcbhbp[1], 1)

      output_df2$Hcbh1_Hdptf1[1] <- round(output_df2$below_hcbhbp[1], 1)
      output_df2$Hcbh1[1] <- 1.5
      output_df2$Hdptf1[1] <- round(output_df2$Hcbh_bp[1], 1)
      output_df2$dptf1[1] <- round(output_df2$Hcbh_bp[1], 0)
      output_df2$Hdist1[1] <- 0.5
      output_df2$effdist1[1] <- NA

      output_df2$last_effdist[1] <- NA
      output_df2$last_Hcbh[1] <- 1.5
      output_df2$last_Hdptf[1] <- round(output_df2$Hcbh_bp[1], 1)
      output_df2$last_Hdist[1] <- 0.5

      output_df2$max_dptf[1] <- round(output_df2$Hcbh_bp[1], 0)
      output_df2$max_Hcbh[1] <- 1.5
      output_df2$max_Hdptf[1] <- round(output_df2$Hcbh_bp[1], 1)
      output_df2$max_Hdist[1] <- 0.5
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

  } else {

    output_df2 <- data.frame(matrix(ncol = 0, nrow = nrow(df_effective)))

    # Add the "treeID" and "treeID1" columns to output_df2
    if (!("treeID" %in% colnames(output_df2))) {
      output_df2 <- cbind(df_effective[c("treeID", "treeID1")], output_df2)
    }

    max_height<-data.frame(df_effective$max_height)
    names(max_height)<-"max_height"

    if(!"max_height" %in% colnames(output_df2)) {
      output_df2 <- cbind(output_df2, df_effective[c("max_height")])
    }
  }
  cummulative_LAD <-output_df2

  return(cummulative_LAD)
}
