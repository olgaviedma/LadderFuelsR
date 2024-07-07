#' Methods to estimated the Crown Base Height of a tree: maximum LAD percentage, maximum distance and the last distance
#' @description
#' This function determines the CBH of a segmented tree using three criteria: maximum LAD percentage, maximum distance and the last distance. In
#' the case of maximum LAD, if the first fuel layer meets the condition, then the function will be select the next fuel layer to avoid confusing
#' with the understory vegetation.
#' @usage
#' get_cbh_metrics(effective_LAD, min_height= 1.5, verbose=TRUE)
#' @param effective_LAD
#' Tree metrics with gaps (distances), fuel base heights, and depths of fuel layers with LAD percentage greater than a threshold (10%).
#' (output of [get_layers_lad()] function).
#' An object of the class text.
#' @param min_height Numeric value for the actual minimum base height (in meters).
#' @param verbose Logical, indicating whether to display informational messages (default is TRUE).
#' @return
#' A data frame giving the Crown Base Height (CBH) of a tree using three criteria: maximum LAD percentage, maximum distance and the last distance.
#' @author Olga Viedma, Carlos Silva, JM Moreno and A.T. Hudak
#'
#' @details
#' List of tree metrics:
#' \itemize{
#'   \item treeID: tree ID with strings and numeric values
#'   \item treeID1: tree ID with only numeric values
#'   \item dptf: Depth of fuel layers (m) after considering distances greater than the actual height bin step
#'   \item effdist: Effective distance between consecutive fuel layers (m) after considering distances greater than any number of steps
#'   \item Hcbh: Base height of each fuel separated by a distance greater than the certain number of steps
#'   \item Hdptf: Height of the depth of fuel layers (m) after considering distances greater than the actual step
#'   \item Hdist: Height of the distance (> any number of steps) between consecutive fuel layers (m)
#'   \item Hcbh_Hdptf - Percentage of LAD values comprised in each effective fuel layer
#'   \item maxlad_Hcbh - Height of the CBH of the segmented tree based on the maximum LAD percentage
#'   \item max_Hcbh - Height of the CBH of the segmented tree based on the maximum distance found in its profile
#'   \item last_Hcbh - Height of the CBH of the segmented tree based on the last distance found in its profile
#'   \item maxlad_ - Values of distance and fuel depth and their corresponding heights at the maximum LAD percentage
#'   \item max_ - Values of distance and fuel depth and their corresponding heights at the maximum distance
#'   \item last_ - Values of distance and fuel depth and their corresponding heights at the last distance
#'   \item nlayers - Number of effective fuel layers
#'   \item max_height - Maximum height of the tree profile
#' }
#'
#' @examples
#' library(magrittr)
#' library(stringr)
#' library(dplyr)
#'
#' # Before running this example, make sure to run get_real_depths().
#' if (interactive()) {
#' effective_LAD <- get_layers_lad()
#' LadderFuelsR::effective_LAD$treeID <- factor(LadderFuelsR::effective_LAD$treeID)
#'
#' trees_name1 <- as.character(effective_LAD$treeID)
#' trees_name2 <- factor(unique(trees_name1))
#'
#' cbh_dist_list <- list()
#'
#' for (i in levels(trees_name2)) {
#' tree1 <- effective_LAD |> dplyr::filter(treeID == i)
#' cbh_dist_metrics <- get_cbh_metrics(tree1, min_height= 1.5, verbose=TRUE)
#' cbh_dist_list[[i]] <- cbh_dist_metrics
#' }
#'
#' # Combine the individual data frames
#' cbh_metrics <- dplyr::bind_rows(cbh_dist_list)
#'
#' # Get original column names
#' original_column_names <- colnames(cbh_metrics)
#'
#' # Specify prefixes
#' desired_order <- c("treeID", "Hcbh", "dptf","effdist","dist", "Hdist", "Hdptf", "max_","last_",
#' "maxlad_", "nlayers")
#'
#'# Identify unique prefixes
#' prefixes <- unique(sub("^([a-zA-Z]+).*", "\\1", original_column_names))
#' # Initialize vector to store new order
#' new_order <- c()
#'
#' # Loop over desired order of prefixes
#' for (prefix in desired_order) {
#'  # Find column names matching the current prefix
#' matching_columns <- grep(paste0("^", prefix), original_column_names, value = TRUE)
#' # Append to the new order
#' new_order <- c(new_order, matching_columns)
#' }
#' # Reorder values
#' cbh_metrics <- cbh_metrics[, new_order]
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
#' @seealso \code{\link{get_layers_lad}}
#' @export
get_cbh_metrics <- function(effective_LAD, min_height= 1.5, verbose=TRUE) {

  if(min_height==0){
    min_height <-0.5
  }

df6a<- effective_LAD

#print(paste(unique(df6a$treeID), collapse = ", "))

if (verbose) {message("Unique treeIDs:", paste(unique(df6a$treeID), collapse = ", "))}

#########################################################
######## Hcbh with max % LAD    ###########################
#########################################################

all_cols <- names(df6a)
numeric_suffix <- as.numeric(gsub("[^0-9]", "", all_cols))
ordered_cols <- all_cols[order(numeric_suffix)]

# Reorder the columns in df6a
df6a <- df6a[, ordered_cols]
if (!("Hdist1" %in% colnames(df6a))) {
  df6a$Hdist1 <- 0.5
}

if(all(is.na(df6a$effdist1)) && df6a$Hcbh1 == min_height) {
  df6a$effdist1<-0
}

df6a <- df6a[, colSums(!is.na(df6a)) > 0]

# Check if df_sub is not empty and contains "Hdist" in column names
col_names <- names(df6a)
hcbh_cols <- grep("^Hcbh\\d+$", col_names, value = TRUE)
hdptf_cols <- grep("^Hdptf\\d+$", col_names, value = TRUE)
hdist_cols <- grep("^Hdist\\d+$", col_names, value = TRUE)
effdist_cols <- grep("^effdist\\d+$", col_names, value = TRUE)

# Create vectors of all Hcbh, Hdist, Hdepth, and effdist values
hcbh_values <- unlist(df6a[1, hcbh_cols])
hdist_values <- unlist(df6a[1, hdist_cols])
hdepth_values <- unlist(df6a[1, hdptf_cols])
effdist_values <- unlist(df6a[1, effdist_cols])

# Identify lad columns based on naming pattern
lad_columns <- grep("^Hcbh\\d+_Hdptf\\d+$", names(df6a), value = TRUE)

# Flatten all lad values into a vector
lad_values <- as.vector(unlist(df6a[, lad_columns]))

# Get the maximum value across all lad columns
max_lad_value <- max(lad_values, na.rm = TRUE)

# Get the index (column number) of the last occurrence of the max value
index_max_lad <- max(which(lad_values == max_lad_value))

# Now you can use index_max_lad to get the corresponding column
max_lad_column <- lad_columns[index_max_lad]

# Get the max lad value from the first row of the dataframe
max_lad_values <- unlist(df6a[1, max_lad_column])

# Extract suffix from max_lad_column
suffix <- as.numeric(sub(".*Hcbh(\\d+)_Hdptf\\d+$", "\\1", max_lad_column))

# Check for any NA values in suffix or df6a$nlayers
if (any(!is.na(suffix)) || any(!is.na(df6a$nlayers))) {

# Check if the maximum lad column coincides with Hcbh1 and nlayers > 1
if (suffix == 1 && df6a$nlayers > 1) {
  suffix <- 2
}

# Update the max_lad_column based on the adjusted suffix
max_lad_column <- paste0("Hcbh", suffix, "_Hdptf", suffix)

# Use the suffix to get the corresponding Hcbh, Hdist and Hdepth columns
max_Hcbh_column <- paste0("Hcbh", suffix)
max_Hdist_column <- paste0("Hdist", suffix)
max_Hdepth_column <- paste0("Hdptf", suffix)
max_effdist_column <- paste0("effdist", suffix)
max_dptf_column <- paste0("dptf", suffix)

# Create the data frame with the selected columns
maxlad_df <- data.frame(
  maxlad_Hcbh = df6a[[max_Hcbh_column]],
  maxlad_Hdist = df6a[[max_Hdist_column]],
  maxlad_Hdptf = df6a[[max_Hdepth_column]],
  maxlad_dptf = df6a[[max_dptf_column]],
  maxlad_effdist = df6a[[max_effdist_column]],
  maxlad_lad = df6a[[max_lad_column]]
)

# Rename the columns for clarity
names(maxlad_df) <- c("maxlad_Hcbh", "maxlad_Hdist", "maxlad_Hdptf", "maxlad_dptf", "maxlad_effdist", "maxlad_lad")
}

#########################################################
######## Hcbh with max height   ###########################
#########################################################

maxlad_cols <- sort(grep("^maxlad_", names(df6a), value = TRUE))
df6a <- df6a[ , !(names(df6a) %in% maxlad_cols)]
df6a <- df6a[, colSums(is.na(df6a)) == 0]

    effdist_cols <- sort(grep("^effdist", names(df6a), value = TRUE))
    effdist_vals <- df6a[1, effdist_cols, drop = FALSE]

    dist_cols <- sort(grep("^dist", names(df6a), value = TRUE))
    dist_vals <- df6a[1, dist_cols, drop = FALSE]

  if(length(effdist_cols) > 0 ) {

    # First, select the 'effdist' columns
    effdist_cols <- names(df6a)[str_detect(names(df6a), "^effdist")]
    # Extract effdist values from the first row
    effdist_values <- df6a[1, grep("effdist", names(df6a))]
    # Find the index of the column with max effdist
    max_effdist_col_index <- tail(which(effdist_values == max(effdist_values, na.rm = TRUE)), n = 1)

    # Define the column names for Hcbh, Hdist, dptf, and Hdptf
    hcbh_cols <- grep("^Hcbh[0-9]+$", names(df6a), value = TRUE)
    hdist_cols <- sort(names(df6a)[grep("Hdist", names(df6a))])
    hdptf_cols <- sort(names(df6a)[grep("^Hdptf", names(df6a))])
    dptf_cols <- sort(names(df6a)[grep("^dptf", names(df6a))])
    lad_cols <- sort(names(df6a)[grep("^Hcbh\\d+_Hdptf\\d+$", names(df6a))])

      suffix <- max_effdist_col_index

      # Get the suffixes of the columns
      hcbh_suffixes <- str_extract(hcbh_cols, "\\d+$")
      hdist_suffixes <- str_extract(hdist_cols, "\\d+$")
      hdptf_suffixes <- str_extract(hdptf_cols, "\\d+$")
      dptf_suffixes <- str_extract(dptf_cols, "\\d+$")
      effdist_suffixes <- str_extract(effdist_cols, "\\d+$")
      lad_suffixes <- str_extract(lad_columns, "\\d+$")

      # Get the columns whose suffixes are greater or equal to the max effdist suffix
      hcbh_col <- hcbh_cols[which.max(as.numeric(hcbh_suffixes) >= as.numeric(suffix))]
      hdist_col <- hdist_cols[which.max(as.numeric(hdist_suffixes) >= as.numeric(suffix))]
      hdptf_col <- hdptf_cols[which.max(as.numeric(hdptf_suffixes) >= as.numeric(suffix))]
      dptf_col <- dptf_cols[which.max(as.numeric(dptf_suffixes) >= as.numeric(suffix))]
      effdist_col <- effdist_cols[which.max(as.numeric(effdist_suffixes) >= as.numeric(suffix))]
      lad_col <- lad_cols[which.max(as.numeric(lad_suffixes) >= as.numeric(suffix))]

      # Create a new dataframe with the required columns
      df6ab <- df6a[,c(hcbh_col, hdist_col, hdptf_col, dptf_col, effdist_col,lad_col)]
      max_df <- df6ab
      # Rename columns with prefix "max_" and remove suffix
      #max_df <- df6ab %>% dplyr::rename_with(.fn = ~ paste0("max_", str_remove(., "\\d+$")))
      names(max_df) <- c("max_Hcbh", "max_Hdist", "max_Hdptf","max_dptf", "max_effdist","max_lad")
    }


  ##############################################################
    ####### the last Hcbh with numeric values ####################
    ######################################

    if(length(effdist_cols) > 0 ) {

    hcbh_cols <- grep("^Hcbh[0-9]+$", names(df6a), value = TRUE)
    last_hcbh_col <- hcbh_cols[length(hcbh_cols)]
    hdist_cols <- grep("^Hdist", names(df6a), value=TRUE)
    last_hdist_col <- hdist_cols[length(hdist_cols)]
    dptf_cols <- grep("^dptf", names(df6a), value=TRUE)
    last_dptf_col <- dptf_cols[length(dptf_cols)]
    hdptf_cols <- grep("^Hdptf", names(df6a), value=TRUE)
    last_hdptf_col <- hdptf_cols[length(hdptf_cols)]
    last_effdist_col <- effdist_cols[length(effdist_cols)]
    lad_col <- grep("^Hcbh\\d+_Hdptf\\d+$", names(df6a), value = TRUE)
    last_lad_col <- lad_col[length(lad_col)]

      # Define the column names for Hcbh, Hdist, dptf, and Hdptf
      last_df <- df6a %>%
        dplyr::select(last_Hcbh = {{last_hcbh_col}},
                      last_Hdist = {{last_hdist_col}},
                      last_Hdptf = {{last_hdptf_col}},
                      last_dptf = {{last_dptf_col}},
                      last_effdist = {{last_effdist_col}},
                      last_lad = {{last_lad_col}})
    }


  df6f<-data.frame(df6a, maxlad_df,max_df, last_df)

 cbh_metrics <- data.frame(df6f)

  return(cbh_metrics)
 }

