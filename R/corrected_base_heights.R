#'
#' Fuels base height after removing distances = 1 m
#'
#' @description
#' This function recalculates fuels base height after removing distances = 1 m,
#' keeping the first "base height" from those consecutive ones separated by a distance = 1.
#'
#' @usage
#' get_real_fbh(depth_metrics)
#'
#' @param depth_metrics
#' Tree metrics with gaps (distances), fuel base heights, and depths
#' (output of [get_depths()] function). An object of the class text.
#'
#' @return
#' A data frame giving the first "base height" from those consecutive ones separated by a distance = 1.
#' This value replaces the values of the next base heights if they are separated by a distance = 1.
#'
#' @author
#' Olga Viedma, Carlos Silva and JM Moreno
#'
#' @details
#' # List of tree metrics:
#' \itemize{
#'   \item treeID: tree ID with strings and numeric values
#'   \item treeID1: tree ID with only numeric values
#'   \item dist: Distance between consecutive fuel layers (m)
#'   \item Hdist - Height of the distance between consecutive fuel layers (m)
#'   \item Hcbh - Height of the first base height from those ones separated by a distance = 1.
#'   \item depth - Depth of fuel layers (m)
#'   \item Hdepth - Height of the depth of fuel layers (m)
#'   \item max_height - Maximum height of the tree profile
#' }
#'
#' @examples
#' ## Not run:
#' library(SSBtools)
#' library(dplyr)
#' library(magrittr)
#'
#' # Tree metrics derived from get_depths() function
#' if (interactive()) {
#'   depth_metrics$treeID <- factor(depth_metrics$treeID)
#' }
#'
#' # Load or create the effective_depth object
#' if (interactive()) {
#'   depth_metrics <- get_depths()
#'   LadderFuelsR::depth_metrics$treeID <- factor(LadderFuelsR::depth_metrics$treeID)
#'
#' trees_name1 <- as.character(depth_metrics$treeID)
#' trees_name2 <- factor(unique(trees_name1))
#'
#' fbh_corr_list <- list()
#'
#' for (i in levels(trees_name2)){
#'   # Filter data for each tree
#'   tree3 <- depth_metrics |> dplyr::filter(treeID == i)
#'   # Get real fbh for each tree
#'   fbh_corr <- get_real_fbh(tree3)
#'   # Store fbh values in a list
#'   fbh_corr_list[[i]] <- fbh_corr
#' }
#'
#' # Combine fbh values for all trees
#' effective_fbh <- dplyr::bind_rows(fbh_corr_list)
#' effective_fbh$treeID <- factor(effective_fbh$treeID)
#' }
#' ## End(Not run)
#'
#' @export get_real_fbh
#' @importFrom dplyr group_by summarise mutate arrange
#' @importFrom magrittr %>%
#' @importFrom SSBtools RbindAll
#' @importFrom gdata startsWith
#' @include gap_fbh.R
#' @include distances_calculation.R
#' @include depths_calculation.R
get_real_fbh <- function (depth_metrics) {

  df<- depth_metrics
  df <- df[, !colSums(is.na(df)) > 0]

       for (i in 1:nrow(df)) {

    effective_fbh <- NULL  # Initialize to NULL for each row
    current_row <- df[i, ]

    hgap_cols <- grep("^gap\\d+$", colnames(current_row))
    hcbh_cols <- grep("^cbh\\d+$", colnames(current_row))

    hgap_na <- is.na(current_row[hgap_cols])
    hcbh_na <- is.na(current_row[hcbh_cols])

    # If all values in current row for Hgap and Hcbh columns are NA, skip this iteration
    if (all(hgap_na) && all(hcbh_na)) {
      next
    }

    df2 <- df

    print(paste("Unique treeIDs:", paste(unique(df2$treeID), collapse = ", ")))

    ###################  rename columns

    # Exclude columns with prefixes: last_, max_, treeID
    cols_to_exclude <- grep("^(max_|treeID)", names(df2), value = TRUE)

    df2 <- df2[ , !(names(df2) %in% cols_to_exclude)]

    # Extract unique prefixes
    prefixes <- unique(gsub("([a-zA-Z]+).*", "\\1", names(df2)))

    # Rename the columns based on the extracted prefixes
    for (prefix in prefixes) {
      # Identify columns with the current prefix
      cols <- grep(paste0("^", prefix), names(df2))

      # Generate new column names with consecutive suffixes
      new_names <- paste0(prefix, 1:length(cols))

      # Assign new names to the columns
      names(df2)[cols] <- new_names
    }


    # Get the column names that start with "Hcbh" followed by a number
    hcbh_cols <- grep("^cbh[0-9]+$", names(df2), value = TRUE)

    # Initialize an empty vector to store the column indices
    hcbh_cols_numeric <- numeric()

    # Iterate through each of the "Hcbh" columns
    for(col in hcbh_cols){
      # Check if the column contains only numeric values
      if(all(sapply(df2[[col]], is.numeric))){
        # If it does, add its index to the vector
        hcbh_cols_numeric <- c(hcbh_cols_numeric, which(names(df2) == col))
      }
    }

    ###################################3

    # Remove columns with only NA values
    df2 <- df2[, colSums(!is.na(df2)) > 0]
    #df2[is.na(df2)] <- 0

    dist_cols <- grep("^dist", names(df2), value = TRUE)
    if (length(dist_cols)==0) {
      df2$dist1<- 0
    }

    # Identify the columns named with 'Hdepth' prefix and values of 0
    cols_to_remove <- grep("^Hdepth", names(df2))[df2[1, grep("^Hdepth", names(df2))] == 0]

    # Check if any columns are identified for removal
    if(length(cols_to_remove) > 0) {
      df2 <- df2[, -cols_to_remove, drop = FALSE]
    }

    # Find columns starting with "Hdist"
    hdist_cols <- (grep("^Hdist", names(df2), value = TRUE))
    hdist_vals <- df2[, hdist_cols]

    hcbh_cols <- grep("^cbh", names(df2), value = TRUE)
    hcbh_vals <- df2[, hcbh_cols]

    dist_cols <- grep("^dist", names(df2), value = TRUE)
    dist_vals <- df2[, dist_cols]

    hdepth_cols <- grep("^Hdepth", names(df2), value = TRUE)
    hdepth_vals <- df2[, hdepth_cols]

    depth_cols <- grep("^depth", names(df2), value = TRUE)
    depth_vals <- df2[, depth_cols]


    ###################  rename columns
    # Extract unique prefixes
    prefixes <- unique(gsub("([a-zA-Z]+).*", "\\1", names(df2)))

    # Rename the columns based on the extracted prefixes
    for (prefix in prefixes) {
      # Identify columns with the current prefix
      cols <- grep(paste0("^", prefix), names(df2))

      # Generate new column names with consecutive suffixes
      new_names <- paste0(prefix, 1:length(cols))

      # Assign new names to the columns
      names(df2)[cols] <- new_names
    }

    hdist_cols <- (grep("^Hdist", names(df2), value = TRUE))
    hdist_vals <- df2[, hdist_cols]

    hcbh_cols <- grep("^cbh", names(df2), value = TRUE)
    hcbh_vals <- df2[, hcbh_cols]

    dist_cols <- grep("^dist", names(df2), value = TRUE)
    dist_vals <- df2[, dist_cols]

    hdepth_cols <- grep("^Hdepth", names(df2), value = TRUE)
    hdepth_vals <- df2[, hdepth_cols]

    depth_cols <- grep("^depth", names(df2), value = TRUE)
    depth_vals <- df2[, depth_cols]

    if (hcbh_vals[[1]] > 1.5 && hdepth_vals [[1]] == 0.5) {
      hdepth_vals [[1]]  <- NULL
      depth_vals [[1]]  <- NULL
    }


    #####################################################
    # Convert vector to one-row data frame if necessary
    if (is.vector(hcbh_vals)) hcbh_vals <- as.data.frame(t(hcbh_vals))
    if (is.vector(hdist_vals)) hdist_vals <- as.data.frame(t(hdist_vals))
    if (is.vector(dist_vals)) dist_vals <- as.data.frame(t(dist_vals))

    # Convert any dist values of 0 to 1
    dist_vals[dist_vals == 0] <- 1

    new_Hcbh <- rep(hcbh_vals[1, 1], ncol(hcbh_vals))

    # Check if all dist values are 1
    if(!all(dist_vals == 1)) {
      for (i in 1:ncol(hdist_vals)) {
        current_hdist <- hdist_vals[1, i]

        # Find the closest Hcbh value that's greater than the current Hdist
        idx <- which(hcbh_vals[1, ] > current_hdist)[1]

        #print(paste("Checking Hdist", i, "=", current_hdist, "Found Hcbh idx:", idx))

        # Check if idx is valid
        if (!is.null(idx) && !is.na(idx) && idx > 0) {

          # If we're at the first position, use the first hcbh value
          if (idx == 1) {
            prev_value <- hcbh_vals[1, 1]
          } else {
            prev_value <- new_Hcbh[idx - 1]
          }

          # If dist value is 1, propagate the last recorded Hcbh value
          if (dist_vals[1, i] == 1) {
            new_Hcbh[idx] <- prev_value
            #print(paste("Dist value is 1. Propagating Hcbh value:", new_Hcbh[idx]))
          }
          # If dist value is >1 and hdist is < the next Hcbh value, update new_Hcbh and propagate
          else if (hdist_vals[1, i] < hcbh_vals[1, idx]) {
            new_Hcbh[idx:ncol(hcbh_vals)] <- hcbh_vals[1, idx]
            #print(paste("Updating Hcbh to:", hcbh_vals[1, idx], "and propagating"))
          }
        }
      }}

    # Transpose the dataframe
    transposed_df <- as.data.frame(new_Hcbh)
    # Remove duplicated rows
    #distinct_df <- transposed_df %>% distinct()
    # Transpose back
    new_transposed_df <- as.data.frame(t(transposed_df))

    colnames(new_transposed_df) <- paste0("Hcbh", seq_len(ncol(new_transposed_df)))
    new_Hcbh_df<-new_transposed_df


    #################################################3

    prefixes <- c("treeID", "Hdist", "Hdepth", "dist", "depth", "max_height")
    df_subset <- df %>%
      dplyr::select(matches(paste0("^", paste(prefixes, collapse = "|"))))

    if (exists("new_Hcbh_df")) {
      effective_fbh <- cbind.data.frame(df_subset, new_Hcbh_df)
    } else {
      effective_fbh <- df_subset
    }

    effective_fbh <- effective_fbh[, colSums(!is.na(effective_fbh)) > 0]

    # Extract columns that start with 'Hcbh'
    Hcbh_cols <- grep("^Hcbh", names(effective_fbh), value = TRUE)

    # Check the conditions
    first_value <- effective_fbh[1, Hcbh_cols[1]]
    subsequent_values <- effective_fbh[1, Hcbh_cols[-1]]

    condition1 <- any(first_value == subsequent_values + 1)
    condition2 <- length(unique(as.numeric(subsequent_values))) == 1

    if(condition1 && condition2) {
      # Set all Hcbh column values to the first Hcbh column value
      effective_fbh[1, Hcbh_cols] <- first_value
    }

    if ("Hdepth0" %in% colnames (effective_fbh) && "depth0" %in% colnames (effective_fbh)) {
      if(any(effective_fbh$Hdepth0 == 0.5) && any(effective_fbh$depth0 == 0)) {
        effective_fbh$Hdepth0 <-NULL
        effective_fbh$depth0 <-NULL
      }}

    effective_fbh <- effective_fbh[, apply(effective_fbh, 2, function(x) x != 0)]

    ##################  rename columns

    # Exclude columns with prefixes: last_, max_, treeID
    cols_to_exclude <- grep("^(sum_dist_|treeID)", names(effective_fbh), value = TRUE)

    effective_fbh <- effective_fbh[ , !(names(effective_fbh) %in% cols_to_exclude)]

    # Extract unique prefixes
    prefixes <- unique(gsub("([a-zA-Z]+).*", "\\1", names(effective_fbh)))

    # Rename the columns based on the extracted prefixes
    for (prefix in prefixes) {
      # Identify columns with the current prefix
      cols <- grep(paste0("^", prefix), names(effective_fbh))

      # Generate new column names with consecutive suffixes
      new_names <- paste0(prefix, 1:length(cols))

      # Assign new names to the columns
      names(effective_fbh)[cols] <- new_names
    }


    treeID<-unique(factor(df$treeID))
    if(!"treeID" %in% colnames(effective_fbh)) {
      effective_fbh <- cbind(df[c("treeID","treeID1")],effective_fbh )
    }

    max_height<-data.frame(df$max_height)
    names(max_height)<-"max_height"

    if(!"max_height" %in% colnames(effective_fbh)) {
      effective_fbh <- cbind(effective_fbh, df[c("max_height")])
    }


    hdist_cols <- (grep("^Hdist", names(effective_fbh), value = TRUE))
    hcbh_cols <- grep("^Hcbh", names(effective_fbh), value = TRUE)
    dist_cols <- grep("^dist", names(effective_fbh), value = TRUE)
    hdepth_cols <- grep("^Hdepth", names(effective_fbh), value = TRUE)
    depth_cols <- grep("^depth", names(effective_fbh), value = TRUE)


    # Convert all relevant columns into numeric vectors
    effective_fbh[hcbh_cols] <- lapply(effective_fbh[hcbh_cols], as.numeric)
    effective_fbh[hdist_cols] <- lapply(effective_fbh[hdist_cols], as.numeric)
    effective_fbh[hdepth_cols] <- lapply(effective_fbh[hdepth_cols], as.numeric)
    effective_fbh[depth_cols] <- lapply(effective_fbh[depth_cols], as.numeric)

    # Iterate over Hcbh columns
    for (hcbh_col in hcbh_cols) {
      # Iterate over Hdist columns
      for (hdist_col in hdist_cols) {
        # Check if the Hcbh value is equal to any Hdist value
        if (effective_fbh[[hcbh_col]] == effective_fbh[[hdist_col]]) {
          effective_fbh[[hcbh_col]] <- effective_fbh[[hcbh_col]] + 1
        }
      }
    }

    if(!exists("effective_fbh")) {
      next  # skip the current iteration if effective_fbh is not generated
    }
  }

  return(effective_fbh)
}
