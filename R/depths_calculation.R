#' Fuels depth in meters
#' @description This function calculates fuels depth as the difference between gaps interleaved between fuel layers minus 1 if the fuel depths are greater than 1.
#' @usage get_depths (LAD_profiles, distance_metrics)
#' @param LAD_profiles original tree Leaf Area Density (LAD) profile (output of [lad.profile()] function in the \emph{leafR} package.
#' An object of the class text
#' @param distance_metrics tree metrics with gaps (distances) and fuel base heights (output of [get_distance()] function).
#' An object of the class text
#' @return A data frame giving fuel layers depth and the height of the depths in meters.
#' @author Olga Viedma, Carlos Silva and JM Moreno
#'
#' @details
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
#' @examples
#' \dontrun{
#' library(magrittr)
#' library(dplyr)
#'
#' # Load the effective_distances object
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
#' metrics_depth <- get_depths(tree1, tree2)
#' metrics_depth_list[[i]] <- metrics_depth
#' }
#'
#' # Combine the individual data frames
#' depth_metrics <- dplyr::bind_rows(metrics_depth_list)
#' }
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
#' @export
get_depths <- function (LAD_profiles,distance_metrics) {

   df1 <- LAD_profiles

  df1$height<-as.numeric(df1$height)
  df1$treeID<-factor(df1$treeID)
  trees_name1a<- as.character(df1$treeID)
  trees_name3<- factor(unique(trees_name1a))

  lad<-df1$lad

  df1_ord<-df1[with(df1, order(lad)), ]

  PERCENTIL_Z <- df1_ord %>%
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

    x1<-df1$height
    y1<-df1$lad

    # Identify missing and infinite values in x and y
    missing_x <- is.na(x1)
    missing_y <- is.na(y1) | is.infinite(y1)

    # Remove missing and infinite values from x and y
    x <- x1[!missing_x & !missing_y]
    y <- y1[!missing_x & !missing_y]

    gaps_perc<- with(df1,
                     ifelse(lad <= PERCENTIL_Z$P5 , "5",
                            ifelse(lad > PERCENTIL_Z$P5 & lad <= PERCENTIL_Z$P25, "25",
                                   ifelse(lad > PERCENTIL_Z$P25 & lad <= PERCENTIL_Z$P50, "50",
                                          ifelse(lad > PERCENTIL_Z$P50 & lad <= PERCENTIL_Z$P75, "75",
                                                 ifelse(lad > PERCENTIL_Z$P75 & lad <= PERCENTIL_Z$P90, "90",
                                                        ifelse(lad > PERCENTIL_Z$P90 & lad <= PERCENTIL_Z$P95, "95",
                                                               ifelse(lad > PERCENTIL_Z$P95 & lad <= PERCENTIL_Z$P99, "99",
                                                                      ifelse(lad >  PERCENTIL_Z$P99, "100",NA)) )))))))

    gaps_perc1 <- data.frame(percentil = as.numeric(gaps_perc))
    gaps_perc2<-cbind.data.frame(df1,gaps_perc1)

  df <- distance_metrics

  df <- df[, !colSums(is.na(df)) > 0]
  # Select only numeric columns
  df1_numeric <- df %>% dplyr::select_if(is.numeric)

  #print(paste("treeID1:", df1_numeric$treeID1))

  # Assuming that columns starting with "gap" or "cbh" are the ones you want to keep
  columns_to_keep <- names(df1_numeric)[grepl("^gap\\d+$|^cbh\\d+$", names(df1_numeric))]
  # Subset the data frame
  kk_copy <- df1_numeric[, columns_to_keep]

  # Sort the column names based on their values in the first row
  sorted_columns <- names(kk_copy)[order(unlist(kk_copy))]

  # Reorder the columns in ascending order
  kk_copy <- kk_copy[, sorted_columns]

  cbh_cols <- which(grepl("cbh", colnames(kk_copy)))
  gap_cols <- which(grepl("gap", colnames(kk_copy)))


depth_file <- data.frame()  # create an empty data frame

if(length(gap_cols) <= 1 & length(cbh_cols) == 0) {
  # No cbh columns, set depth_data to NA
  depth_file <- data.frame(NA)
  names(depth_file) <- "depth_crown_0"
}

if(length(gap_cols) <= 1 & length(cbh_cols) != 0 & max(kk_copy[, gap_cols]) < max(kk_copy[, cbh_cols])) {
  depth_file <-data.frame(max(kk_copy[, cbh_cols]) - min(kk_copy[, gap_cols]))
  colnames(depth_file) <- "depth"
  hdepth<-data.frame(max(kk_copy[, cbh_cols]))
  colnames(hdepth) <- "Hdepth"
  depth_file<-cbind.data.frame(hdepth,depth_file)

}

if(length(gap_cols) <= 1 & length(cbh_cols) != 0 & min(kk_copy[, gap_cols]) < min(kk_copy[, cbh_cols]) && nrow(depth_file) == 0 && ncol(depth_file) == 0) {
  depth_file <-data.frame(abs(min(kk_copy[, gap_cols]) - max(kk_copy[, cbh_cols])))
  names(depth_file) <- "depth"
  hdepth<-data.frame(max(kk_copy[, cbh_cols]))
  names(hdepth) <- "Hdepth"
  depth_file<-cbind.data.frame(hdepth,depth_file)

}

if(length(gap_cols) <= 1 & length(cbh_cols) != 0 & max(kk_copy[, gap_cols]) > max(kk_copy[, cbh_cols])&& nrow(depth_file) == 0 && ncol(depth_file) == 0) {
  depth_file <-data.frame(max(kk_copy[, gap_cols]) - min(kk_copy[, cbh_cols]))
  names(depth_file) <- "depth"
  hdepth<-data.frame(max(kk_copy[, gap_cols])-1)
  names(hdepth) <- "Hdepth"
  depth_file<-cbind.data.frame(hdepth,depth_file)

}


new_depth_file <- depth_file  # Create a new dataframe to store the modified results

if (length(gap_cols) > 1 && length(cbh_cols) != 0 && nrow(new_depth_file) == 0 && ncol(new_depth_file) == 0) {
  gap_cols <- gap_cols[!is.na(kk_copy[, gap_cols])]
  cbh_cols <- cbh_cols[!is.na(kk_copy[, cbh_cols])]

  gap_indices <- sort(gap_cols)
  gap_indices <- gap_indices[!is.na(gap_indices)]
  depth_gap_values <- kk_copy[, gap_indices]

  gap_diffs <- numeric(length(gap_indices) - 1)
  cbh_index_values1 <- numeric(0)  # Create an empty vector to store cbh_index values

  found_cbh_between_gaps <- FALSE  # This flag will be set to TRUE if we find a cbh index between any two gap indices

  for (i in 1:(length(gap_indices) - 1)) {
    gap_index1 <- gap_indices[i]
    gap_index2 <- gap_indices[i + 1]
    diffs <- kk_copy[, gap_index2] - kk_copy[, gap_index1]

    if (any(!is.na(cbh_cols))) {
      diffs1 <- numeric(length(gap_indices) - 1)

      for (j in 1:length(cbh_cols)) {
        cbh_index <- cbh_cols[j]

        if (cbh_index > gap_index1 && cbh_index < gap_index2) {
          diffs1 <- diffs[!is.na(kk_copy[, cbh_index])]
          cbh_index_values1 <- c(cbh_index_values1, cbh_index)  # Store cbh_index value
          found_cbh_between_gaps <- TRUE
        }
      }

      gap_diffs[i] <- abs(diffs1)
    }

    Hdepth <- kk_copy[1, gap_index2] - 1
    new_depth_file <- rbind(new_depth_file, c(Hdepth, t(gap_diffs)))
  }

  if (!found_cbh_between_gaps) {
    new_depth_file <- data.frame()
  } else {
    colnames(new_depth_file) <- c("Hdepth", paste0("X", 1:(length(gap_indices) - 1)))
  }
}

#new_depth_file <- new_depth_file[!c(FALSE, diff(new_depth_file$Hdepth) == 1), ]


# Depurating output format

if (exists("new_depth_file") && !is.null(new_depth_file) && nrow(new_depth_file) > 1 && ncol(new_depth_file) != 0) {

  # Check if all X columns are equal to 0 or NA
  all_zeros <- apply(new_depth_file[, -1], 1, function(x) all(is.na(x) | x == 0))

  # Remove rows where all X columns are equal to 0 or NA
  new_depth_file <- new_depth_file[!all_zeros, ]

  # Replace zeros with NA (excluding 'Hdepth' column)
  new_depth_file[, -1][new_depth_file[, -1] == 0] <- NA

  # Generate 'depth' column that contains the last non-NA value from each row
  new_depth_file$depth <- apply(new_depth_file[, -1, drop = FALSE], 1, function(x) {
    non_na_vals <- tail(na.omit(x), n = 1)  # exclude NA values and get the last non-NA value
    if (length(non_na_vals) == 0) NA else non_na_vals
  })

  # Create a 'group' column by concatenating the values in X columns in each row
  new_depth_file$group <- apply(new_depth_file[,-c(1, ncol(new_depth_file))], 1, paste, collapse = "_")

  # Generate a 'group_id' column to identify unique groups
  new_depth_file$group_id <- match(new_depth_file$group, unique(new_depth_file$group))

  # For each unique 'group_id', keep only the first row
  new_depth_file <- new_depth_file[!duplicated(new_depth_file$group_id),]

  # Remove unnecessary columns
  new_depth_file <- new_depth_file[,-c(ncol(new_depth_file)-1, ncol(new_depth_file))]

  # Define the final data frame
  group_list <- list()

  # For each unique 'Hdepth' value, append the Hdepth/depth pair to 'group_list'
  for(i in seq_along(unique(new_depth_file$Hdepth))){
    group_df <- new_depth_file[new_depth_file$Hdepth == unique(new_depth_file$Hdepth)[i], c("Hdepth", "depth")]

    # Set column names
    colnames(group_df) <- paste0(c("Hdepth_", "depth_"), i)

    # Add the data frame to 'group_list'
    group_list[[length(group_list)+1]] <- group_df
  }

  # Combine all data frames in 'group_list' into 'depth_file'
  depth_file <- do.call(cbind, group_list)

} else if (nrow(new_depth_file) == 1) {

  if (startsWith(names(new_depth_file)[1], "Hdepth")) {
    depth_file <- new_depth_file[, c(1, 2)]
    colnames(depth_file) <- c("Hdepth_1", "depth_1")
  } else {
    depth_file <- new_depth_file[, c(2, 1)]
    colnames(depth_file) <- c("Hdepth_1", "depth_1")
  }
}



if(exists("cbh_index_values1")){

  depth_cbh_index_values<-kk_copy[,cbh_index_values1]

} else {

  depth_cbh_index_values<-data.frame(NA)
  names(depth_cbh_index_values)<-"depths_cbh"
}


if (!exists("depth_file")) {
  depth_file <- data.frame(NA)
  names(depth_file) <- "depth_crown_0"
}

#########################
cbh_cols <- grep("^cbh", names(kk_copy))
gap_cols <- grep("^gap", names(kk_copy))

if (any(diff(cbh_cols) == 1)) {

  consec_cbh_cols <- which(diff(cbh_cols) == 1)
  first_consec_cbh_col <- min(consec_cbh_cols)
  last_consec_cbh_col <- max(consec_cbh_cols) + 1 # add 1 to include last column

  # Extract only the consecutive cbh columns
  consec_cbh_vals <- kk_copy[, cbh_cols[first_consec_cbh_col:last_consec_cbh_col]]

} else {
  consec_cbh_vals <- data.frame(NA)
  names(consec_cbh_vals) <- "depth_crown_0"
}

# Check if there are any consecutive "cbh" columns
if (!is.na(consec_cbh_vals) && length(gap_cols) > 1) {

  # Get the starting and ending values of the consecutive "cbh" columns
  consec_cbh_start <- min(consec_cbh_vals)
  consec_cbh_tail <- max(consec_cbh_vals)
  gaps_vals <- kk_copy[, gap_cols]

  # Find the position of the last gap column before the last consecutive cbh column
  last_gap_col <- max(gap_cols[gap_cols < cbh_cols[last_consec_cbh_col]])

  # Subset the last gap column

  if (!is.finite(last_gap_col) || any(is.na(consec_cbh_vals))) {

    consec_cbh_vals <- data.frame(NA)

  }else {

    last_gap_vals <- kk_copy[, last_gap_col]
  }

  if (exists("last_gap_vals") && any(!is.na(consec_cbh_vals))) {

    if (last_gap_vals < consec_cbh_start && max(gaps_vals) < consec_cbh_tail && length(gap_cols) > 1) {
      # Get the difference between the last "gap" column and the ending value of the consecutive "cbh" columns
      depth_file11 <- abs(consec_cbh_tail - last_gap_vals)
      hdepth11<-data.frame(consec_cbh_tail)
      names(hdepth11)<-"Hdepth11"
      depth_file11<-cbind.data.frame(hdepth11,depth_file11)
    } else {
      depth_file11 <- data.frame(NA)
      names(depth_file11) <- "depth_crown_0"
    }
  }
}

##################################################

cbh_cols <- grep("^cbh", names(kk_copy))
gap_cols <- grep("^gap", names(kk_copy))

if (max(kk_copy[,gap_cols]) < max(kk_copy[,cbh_cols])) {

  depth_file12 <- abs(max(kk_copy[,cbh_cols]) - max(kk_copy[,gap_cols]))
  hdepth12<-data.frame(max(kk_copy[,cbh_cols]))
  names(hdepth12)<-"Hdepth12"
  depth_file12<-cbind.data.frame(hdepth12, depth_file12)

} else {
  depth_file12 <- data.frame(NA)
  names(depth_file12) <- "depth_crown_0"
}

if (!exists("depth_file11")) {
  depth_file11 <- data.frame(NA)
  names(depth_file11) <- "depth_crown_0"
}
if (!exists("depth_file12")) {
  depth_file12 <- data.frame(NA)
  names(depth_file12) <- "depth_crown_0"
}

##########################################

last_Hdepth_col_name <- tail(grep("Hdepth", names(depth_file), value = TRUE), 1)
last_depth <- depth_file[, last_Hdepth_col_name]

if(exists("last_depth") || !is.na(last_depth)) {

  last_depth <- depth_file[, last_Hdepth_col_name]

  # All elements in depth_file are not equal to 0 or missing values
  if ((nrow(depth_file) != 0 && any(!is.na(depth_file))) &&
      !is.na(depth_file11) && length(depth_file11) > 0 && any(!is.na(depth_file11)) &&
      !is.na(depth_file12) && length(depth_file12) > 0 && exists("last_gap_vals")) {

    if(last_gap_vals != max(kk_copy[, gap_cols]) &&
       depth_file11[,1] != depth_file12[,1] &&
       (last_depth != depth_file11[,1] && last_depth != depth_file12[,1])) {

      depth_file <- cbind.data.frame(depth_file, depth_file11, depth_file12)

    } else {
      depth_file <- depth_file
    }
  }

  if (length(depth_file11) > 0 && any(!is.na(depth_file11)) && !is.na(last_depth) &&
      nrow(depth_file) != 0 &&
      (is.na(depth_file12) || length(depth_file12) == 0)) {

    if(last_depth != depth_file11[,1]){

      depth_file <- cbind.data.frame(depth_file, depth_file11)

    } else {
      depth_file <- depth_file
    }
  }

  if ((length(depth_file11) > 0 || any(!is.na(depth_file11))) &&
      (!is.na(depth_file12) || length(depth_file12) != 0)  &&
      nrow(depth_file) != 0 &&
      (!is.na(last_depth) && !is.na(depth_file11) && last_depth !=  depth_file11[,1]) &&
      (!is.na(last_depth) && !is.na(depth_file12) && last_depth !=  depth_file12[,1])) {

    if(depth_file11[,1] == depth_file12[,1]){

      depth_file <- cbind.data.frame(depth_file, depth_file11)

    } else {
      depth_file <- depth_file
    }
  }


  if ((nrow(depth_file) != 0 ) &&
      (!is.na(depth_file11) && length(depth_file11) > 0 && any(!is.na(depth_file11))) &&
      (!is.na(depth_file12) && length(depth_file12) > 0)) {

    if(depth_file11[,1] != depth_file12[,1]) {

      depth_file <- cbind.data.frame(depth_file, depth_file11,depth_file12)

    } else {
      depth_file <- depth_file
    }
  }


  if ((nrow(depth_file) != 0 ) &&
      (!is.na(depth_file11) || length(depth_file11) != 0 || any(!is.na(depth_file11))) &&
      (is.na(depth_file12) || length(depth_file12) == 0) &&
      any(!is.na(last_depth)) && !is.na(depth_file11)) {

    if(last_depth !=  depth_file11[,1] && !("depth_file11" %in% colnames(depth_file))){

      depth_file <- cbind.data.frame(depth_file, depth_file11)

    } else {
      depth_file <- depth_file
    }
  }


  if (nrow(depth_file) != 0 &&
      (!is.na(depth_file11) || length(depth_file11) != 0 || any(!is.na(depth_file11))) &&
      any(!is.na(last_depth)) && !is.na(depth_file11)) {

    if ((last_gap_vals < consec_cbh_start && max(gaps_vals) < consec_cbh_tail) && !("depth_file11" %in% colnames(depth_file))) {
      depth_file <- cbind.data.frame(depth_file, depth_file11)
    } else {
      depth_file <- depth_file
    }
  }


  if (nrow(depth_file) != 0 &&
      (is.na(depth_file11) || length(depth_file11) == 0 || any(is.na(depth_file11))) &&
      (!is.na(depth_file12) || length(depth_file12) != 0) &&
      any(!is.na(last_depth)) && !is.na(depth_file12)) {

    if (last_depth != depth_file12[,1]) {
      depth_file <- cbind.data.frame(depth_file, depth_file12)
    }
  }


  if (nrow(depth_file) != 0 &&
      (!is.na(depth_file11) || length(depth_file11) == 0 || any(!is.na(depth_file11))) &&
      (!is.na(depth_file12) || length(depth_file12) == 0) &&
      any(!is.na(last_depth)) && !is.na(depth_file12) &&
      (!is.na(last_depth) && !is.na(depth_file12))) {

    if (any(!is.na(depth_file11)) && depth_file11[,1] != depth_file12[,1] &&
        (!is.na(last_depth) && last_depth != depth_file12[,1])) {

      depth_file <- cbind.data.frame(depth_file, depth_file12)
    } else {
      depth_file <- depth_file
    }
  }

  if (nrow(depth_file) != 0 && length(gap_cols) > 1 &&
      (is.na(depth_file11) || length(depth_file11) == 0 || any(is.na(depth_file11))) &&
      (!is.na(depth_file12) || length(depth_file12) > 0 || any(!is.na(depth_file12)))) {

    if (!is.na(depth_file12) && last_depth == depth_file12[,1] && (max(kk_copy[,gap_cols]) < max(kk_copy[,cbh_cols]))) {
      depth_file <- cbind.data.frame(depth_file,depth_file12)
    } else {
      depth_file <- depth_file
    }
  }
}


if ((is.na (depth_file)|| nrow(depth_file) == 0) &&
    (!is.na(depth_file11) || length(depth_file11) != 0 || any(!is.na(depth_file11))) &&
    (!is.na(depth_file12) || length(depth_file12) > 0 || any(!is.na(depth_file12)))) {

  if (!is.na(depth_file11) && !is.na(depth_file12) && depth_file11[,1] != depth_file12[,1]) {
    depth_file <- cbind.data.frame(depth_file11, depth_file12)
  } else {
    depth_file <- depth_file11
  }
}

if ((is.na (depth_file)||nrow(depth_file) == 0) &&
    (is.na(depth_file11) || length(depth_file11) == 0 || any(is.na(depth_file11))) &&
    (all(!is.na(depth_file12)) || length(depth_file12) > 0 || any(!is.na(depth_file12)))) {

  if (any(!is.na(depth_file12))) {
    depth_file <- cbind.data.frame(depth_file,depth_file12)
  } else {
    depth_file <- depth_file
  }
}

if (!exists("depth_file") || nrow(depth_file) == 0) {
  depth_file <- data.frame(NA)
  names(depth_file) <- "depth_crown_0"
}

# Remove columns with zeros or NAs
depth_file <- data.frame(depth_file[, colSums(depth_file != 0, na.rm = TRUE) > 0])

if (length(depth_file)==0) {
  depth_file <- data.frame(NA)
  names(depth_file) <- "depth_crown_0"
}


### DIFFERENCE BETWEEN DEPTH AND DISTANCE (1) TO GET THE REAL "LAYER DEPTH"  #########################3

# Check if "depth_file" exists and if it has more than 0 rows
if (exists("depth_file") && nrow(depth_file) > 0) {

  # Identify the column names to subtract 1 from
  columns_to_subtract <- grep("^depth", names(depth_file), value = TRUE)

  # Subtract 1 from the specified columns if they are greater than 1
  for (column in columns_to_subtract) {
    depth_file[[column]] <- ifelse(depth_file[[column]] > 1, depth_file[[column]] - 1, depth_file[[column]])
  }
}

if (exists("depth_file") && !is.null(depth_file) && is.data.frame(depth_file) && nrow(depth_file) != 0) {

  depth_file_transposed <- data.frame(t(depth_file))

  if (!is.na(depth_file_transposed[1,1]) && depth_file_transposed[1,1] == 0) {

    # Get the relevant columns from depth_file_transposed, skipping the rows where t.depth_file. is 0
    depth_file_col <- depth_file_transposed[depth_file_transposed$t.depth_file. != 0, "t.depth_file."]

    # Calculate the differences
    result <- depth_file_col - 1

    # Convert the result back to a data frame
    depth_file <- as.data.frame(t(result))
  }
}


if (!exists("depth_file") &&  is.null(depth_file) && !is.data.frame(depth_file) && nrow(depth_file) == 0) {
  depth_file <- data.frame(NA)
  names(depth_file) <- "depth_crown_0"
}


if ((exists("depth_file") || is.data.frame(depth_file) && any(!is.na(depth_file)) && length(depth_file) != 0 && !is.null(depth_file))) {
  depth_data <- depth_file
}

######################################
#### DEPTH TOT: DIFFERENCE BETWEEN THE FIRST CONSECUTIVE CBH COLUMN AND THE FIRST GAP ENCOUNTERES
###########################################################33

cbh_cols <- which(grepl("cbh", colnames(kk_copy)))
gap_cols <- which(grepl("gap", colnames(kk_copy)))

gap_cols <- gap_cols[!is.na(kk_copy[, gap_cols])]
cbh_cols <- cbh_cols[!is.na(kk_copy[, cbh_cols])]


if (nrow(kk_copy) !=0 && length(gap_cols)!= 0 && length(cbh_cols)== 0 ) {

  if (min(kk_copy[,gap_cols]) > 1.5) {
    depth1 <-  min(kk_copy[,gap_cols])-0.5
    depth0 <- data.frame(depth1)
    names(depth0) <- "depth0"
    Hdepth0<-data.frame(min(kk_copy[,gap_cols])-1)
    names(Hdepth0)<-"Hdepth0"
    depth0<-cbind.data.frame(Hdepth0,depth0)
  }

  if (min(kk_copy[,gap_cols]) == 1.5) {

    min_gap <- min(kk_copy[, gap_cols])
    next_gap_col <- gap_cols[which(kk_copy[, gap_cols] > min_gap)[1]]
    next_gap <- kk_copy[[next_gap_col]]

    depth1 <- next_gap - min_gap
    depth0 <- data.frame(depth1)
    names(depth0) <- "depth0"
    Hdepth0<-data.frame(next_gap-1)
    names(Hdepth0)<-"Hdepth0"
    depth0<-cbind.data.frame(Hdepth0,depth0)
  }
}


if (nrow(kk_copy) != 0 && length(gap_cols) == 0 && length(cbh_cols) != 0 && (!exists("depth0") || nrow(depth0) == 0)) {

  if (any(diff(cbh_cols) == 1)) {
    consec_cbh_cols <- diff(cbh_cols) == 1
    depth1 <- max(df$height) - min(kk_copy[, cbh_cols])
    depth0 <- data.frame(depth1)
    names(depth0) <- "depth0"
    Hdepth0<-data.frame(max(df$height))
    names(Hdepth0)<-"Hdepth0"
    depth0<-cbind.data.frame(Hdepth0,depth0)


  } else if (length(cbh_cols) == 1 && min(kk_copy[, cbh_cols]) == max(df$height)){

    depth1 <-  min(kk_copy[, cbh_cols])-0.5
    depth0 <- data.frame(depth1)
    names(depth0) <- "depth0"
    Hdepth0<-data.frame(min(kk_copy[, cbh_cols]))
    names(Hdepth0)<-"Hdepth0"
    depth0<-cbind.data.frame(Hdepth0,depth0)

  }else if (length(cbh_cols) == 1 && min(kk_copy[, cbh_cols]) < max(df$height)){

    depth1 <-  max(df$height) - min(kk_copy[, cbh_cols])
    depth0 <- data.frame(depth1)
    names(depth0) <- "depth0"
    Hdepth0<-data.frame(max(df$height))
    names(Hdepth0)<-"Hdepth0"
    depth0<-cbind.data.frame(Hdepth0,depth0)

  }
}

if (any(diff(cbh_cols) == 1) && min(kk_copy[, cbh_cols]) >= 1.5 && min(kk_copy[, gap_cols]) > min(kk_copy[, cbh_cols]) && (!exists("depth0") ||  nrow(depth0)==0)) {
  consec_cbh_cols <- diff(cbh_cols) == 1
  if (any(consec_cbh_cols)) {
    consec_cbh_start <- min(kk_copy[1, consec_cbh_cols])
    if (length(consec_cbh_start) > 0 && any(min(kk_copy[, gap_cols]) > consec_cbh_start)) {
      depth1 <- (min(kk_copy[, gap_cols]) - consec_cbh_start)
      depth0 <- data.frame(depth1)
      names(depth0) <- "depth0"
      Hdepth0<-data.frame(min(kk_copy[, gap_cols])-1)
      names(Hdepth0)<-"Hdepth0"
      depth0<-cbind.data.frame(Hdepth0,depth0)
    }

  }else {

    depth0 <- depth0

  }
}

##################################################33

if (nrow(kk_copy) != 0 && length(gap_cols) != 0 && length(cbh_cols) != 0 && (!exists("depth0") || nrow(depth0) == 0)) {

  if (!is.na(depth_file[, 1]) && depth_file[, 1] != 0 && !is.na(min(kk_copy[, gap_cols])) && !is.na(min(kk_copy[, cbh_cols])) &&
      min(kk_copy[, gap_cols]) > min(kk_copy[, cbh_cols])) {

    if (any(min(kk_copy[, cbh_cols]) > 1.5) && any(min(kk_copy[, gap_cols]) > min(kk_copy[, cbh_cols]))) {

       height<-gaps_perc2$height
      # percent1 <- gaps_perc2 %>% dplyr::filter(height < min(kk_copy[, cbh_cols]))
      percent2 <- gaps_perc2 %>% dplyr::filter(height < min(kk_copy[, gap_cols]))

    } else {

      percent2 <- data.frame(NA)
      names(percent2) <- "percentil"
    }
    if(nrow(percent2) != 0 && any(!is.na(percent2)) && any(percent2$percentil <= 5) && (!exists("depth0") || nrow(depth0)==0)) {

      percentil<-percent2$percentil
      height_percent1<-percent2|> dplyr::filter(percentil <= 5)
      runs <- rle(height_percent1$percentil)

      if (any(runs$lengths > 1)) {
        first_percent_col <- max (height_percent1$height)
        depth1 <-  min(kk_copy[,gap_cols]) - first_percent_col
        depth0 <- data.frame(depth1)
        names(depth0) <- "depth0"
        Hdepth0<-data.frame(min(kk_copy[, gap_cols])-1)
        names(Hdepth0)<-"Hdepth0"
        depth0<-cbind.data.frame(Hdepth0,depth0)

      } else {
        depth1 <-   min(kk_copy[,gap_cols]) - max (height_percent1$height)
        depth0 <- data.frame(depth1)
        names(depth0) <- "depth0"
        names(depth0) <- "depth0"
        Hdepth0<-data.frame(min(kk_copy[, gap_cols])-1)
        names(Hdepth0)<-"Hdepth0"
        depth0<-cbind.data.frame(Hdepth0,depth0)

      }
    } else if(all(nrow(percent2) != 0 && any(!is.na(percent2)) && any(percent2$percentil > 25))){

      height_percent2<-percent2|> dplyr::filter(percentil > 25)

      depth1 <-  (min(kk_copy[,gap_cols]) -min(height_percent2$height))-1.5
      depth0 <- data.frame(depth1)
      names(depth0) <- "depth0"
      Hdepth0<-data.frame(min(kk_copy[, gap_cols])-1)
      names(Hdepth0)<-"Hdepth0"
      depth0<-cbind.data.frame(Hdepth0,depth0)

    }
  }
}

if (min(kk_copy[, gap_cols]) < min(kk_copy[, cbh_cols]) && (!exists("depth0") || nrow(depth0) == 0)) {

  if (any(diff(gap_cols) == 1)) {

    # Find the first set of consecutive "gap" columns
    first_set <- rle(grepl("^gap", names(kk_copy)))
    first_set_length <- first_set$lengths[1]
    first_set_cols <- gap_cols[1:first_set_length]

    # Subset the first set of consecutive "gap" columns
    first_cons_gap <- kk_copy[, first_set_cols, drop = FALSE]

    # Retrieve the first column of the subsetted data frame
    first_gap_column <- first_cons_gap[, 1]

    # Retrieve the last value of the first set of consecutive "gap" columns (gap2)
    last_gap_value <- first_cons_gap[, first_set_length]

    # Filter the rows where the height is less than the last gap value
    percent3 <- gaps_perc2[gaps_perc2$height < last_gap_value, ]

    # Filter the remaining rows where the percentile is greater than 25
    height_percent3 <- percent3[percent3$percentil > 25, ]

    depth0 <- height_percent3$height - first_gap_column
    depth0 <- data.frame(depth0)
    names(depth0) <- "depth0"
    Hdepth0<-data.frame(height_percent3$height)
    names(Hdepth0)<-"Hdepth0"
    depth0<-cbind.data.frame(Hdepth0,depth0)

  } else {

    depth1 <- c(min(kk_copy[, gap_cols]) - 1.5)
    depth0 <- data.frame(depth1)
    names(depth0) <- "depth0"
    Hdepth0<-data.frame(min(kk_copy[, gap_cols])-1)
    names(Hdepth0)<-"Hdepth0"
    depth0<-cbind.data.frame(Hdepth0,depth0)

  }
}


if (any(length(cbh_cols) == 1 && min(kk_copy[, cbh_cols]) == 1.5) && (!exists("depth0")|| nrow(depth0)==0)) {
  # if there is only one cbh column, calculate difference between it and first gap column
  cbh_col <- cbh_cols[1]
  first_gap_col <- min(gap_cols[gap_cols > cbh_col])
  if (!is.na(cbh_col) && !is.na(first_gap_col)) {
    depth1 <- kk_copy[1, first_gap_col] - kk_copy[1, cbh_col]
    depth0 <- data.frame(depth1)
    names(depth0) <- "depth0"
    Hdepth0<-data.frame(kk_copy[1, first_gap_col]-1)
    names(Hdepth0)<-"Hdepth0"
    depth0<-cbind.data.frame(Hdepth0,depth0)


  }
}


if (any(min(kk_copy[,cbh_cols]) > 1.5) &&  any(min(kk_copy[,gap_cols]) <  min(kk_copy[,cbh_cols])) && (!exists("depth0") || nrow(depth0)==0)) {
  depth1 <-  min(kk_copy[,gap_cols])-1.5
  depth0 <- data.frame(depth1)
  names(depth0) <- "depth0"
  Hdepth0<-data.frame(min(kk_copy[,gap_cols])-1)
  names(Hdepth0)<-"Hdepth0"
  depth0<-cbind.data.frame(Hdepth0,depth0)
}


if (any(min(kk_copy[,cbh_cols]) == 1.5) &&  any(min(kk_copy[,gap_cols]) >  min(kk_copy[,cbh_cols])) && (!exists("depth0") || nrow(depth0)==0)) {
  depth1 <-  min(kk_copy[,gap_cols])-min(kk_copy[,cbh_cols])
  depth0 <- data.frame(depth1)
  names(depth0) <- "depth0"
  names(depth0) <- "depth0"
  Hdepth0<-data.frame(min(kk_copy[,gap_cols])-1)
  names(Hdepth0)<-"Hdepth0"
  depth0<-cbind.data.frame(Hdepth0,depth0)
}


if ( !exists("depth0")){
  depth0 <- data.frame(NA)
  names(depth0) <- "depth0"
}

if(any(is.na(depth_data)) && any(is.na(depth0)) || length(depth_data) == 0 ) {

  depth_data <- data.frame(NA)
  names(depth_data) <- "depth_crown_0"

}

#####################################################
#### ADJUST DEPTH DATA ACCORDING TO DEPTH0 ############

if(exists("depth_cbh_index_values") && any(!is.na(depth_cbh_index_values))) {

  depth_cbh_index_values<-data.frame(depth_cbh_index_values)
  min_cbh_depth<-depth_cbh_index_values[,1]

} else {

  depth_cbh_index_values<-data.frame(NA)
  names(depth_cbh_index_values)<-"depths_cbh"
}


min_cbh<-data.frame(min(kk_copy[, cbh_cols]))
max_cbh<-data.frame(max(kk_copy[, cbh_cols]))


if (exists("depth_data") && any(!is.na(depth_cbh_index_values)) && any(!is.na(min_cbh))) {
  if (!is.na(min_cbh) && !is.na(depth_cbh_index_values[, 1])) {
    if (min_cbh >= depth_cbh_index_values[, 1]) {
      depth_data <- data.frame(depth_data)
    }
  }
}


if (exists("depth_data") && length(cbh_cols) > 1 && any(!is.na(depth_cbh_index_values)) && any(!is.na(min_cbh)) && (exists("depth0") || nrow(depth0) != 0)) {
  if (!is.na(min_cbh) && any(!is.na(depth_cbh_index_values[, 1])) && any(!is.na(depth0[, 1])) && all(depth0[, 1] == depth_data[, 1])) {
    if (exists("last_gap_vals") && last_gap_vals < consec_cbh_start && max(gaps_vals) < consec_cbh_tail) {
      depth_data <- depth_data
    }
  }
}

if (exists("depth_data") && length(cbh_cols) > 1 && any(!is.na(depth_cbh_index_values)) && any(!is.na(min_cbh)) && (exists("depth0") || nrow(depth0) != 0)) {
  if (!is.na(min_cbh) && any(!is.na(depth_cbh_index_values[, 1])) && any(!is.na(depth0[, 1])) && all(depth0[, 1] == depth_data[, 1])) {
    if (min_cbh == depth_cbh_index_values[, 1] && last_gap_vals < consec_cbh_start && max(gaps_vals) > consec_cbh_tail) {
      depth_data <- data.frame(depth_data[, -c(1:2)])
    }
  }
} else {
  depth_data <- depth_data
}


if (exists("depth_data") && any(is.na(depth_cbh_index_values)) && !is.na(min_cbh)) {

  depth_data <- data.frame(depth_data)

  if (any(!is.na(depth0)) && any(!is.na(depth_data)) &&
      (any(depth0[, 1] == depth_data[, 1])) && max(kk_copy[,gap_cols]) > max(kk_copy[,cbh_cols])) {

    depth_data <- data.frame(depth_data[,- c(1:2)])

  } else if (max(kk_copy[,gap_cols]) < max(kk_copy[,cbh_cols])) {
    depth_data <- depth_data
  }
}

if (any(!is.na(depth0)) && any(!is.na(depth_data)) &&
    (any(depth0[, 1] == depth_data[, 1] +1)) && max(kk_copy[,gap_cols]) > max(kk_copy[,cbh_cols]) && length(gap_cols) == 1) {

  depth_data <- data.frame(depth_data[,- c(1:2)])

} else {
  depth_data <- depth_data
}


if (exists("depth_data")) {
  if (any(!is.na(depth_data)) && any(depth_data[, 1] == 0 )) {

    depth_data <- depth_data[,- c(1:2)]
  }
}

# Remove columns with zeros or NAs
depth_data <- data.frame(depth_data[, colSums(depth_data != 0, na.rm = TRUE) > 0])


#######################################

  if( nrow(kk_copy) ==0) {
    depth0 <- data.frame(max(df$height) - min(df$height))
    names(depth0) <- "depth0"
  }

  if( nrow(kk_copy) ==0) {
    depth_data <- data.frame(NA)
    names(depth_data) <- "depth_crown_0"
  }

  #######################################

  # Check if depth is not empty
  if (length(depth_data) != 0) {

    depth_data <- data.frame(depth_data)

    old_names_depth_data <- colnames(depth_data)
    new_names_depth_data <- vector(mode="character", length=length(old_names_depth_data))  # initialize new_names_depth_data

    # Initialize counter for column indices
    column_index <- 1

    for (i in seq_along(old_names_depth_data)) {
      # For odd indices, assign 'depth' prefix
      if (i %% 2 != 0) {
        new_names_depth_data[i] <- paste("Hdepth", column_index, sep="")
      }
      # For even indices, assign 'Hdepth' prefix
      else {
        new_names_depth_data[i] <- paste("depth", column_index, sep="")
        # Increase column index every time an 'Hdepth' is assigned (i.e., every complete pair)
        column_index <- column_index + 1
      }
    }

    colnames(depth_data) <- new_names_depth_data
  }


treeID<-unique(factor(df$treeID))
treeID1<-unique(factor(df$treeID1))

max_height<-data.frame(df$max_height)
names(max_height)="max_height"

    depth_metrics <- cbind.data.frame(df, depth0, depth_data)

return(depth_metrics)
}







