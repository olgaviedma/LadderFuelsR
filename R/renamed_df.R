#'
#' Rename and reorder columns (II)
#'
#' @description This function deals with concatenated column names, reorders columns and appends numeric suffixes.
#' Don´t run it. It is an internal function.
#'
#' @usage get_renamed_df (df)
#'
#' @param df internal data frame derived from [get_layers_lad()] function
#'
#' @examples
#' ## Not run:
#' # get_renamed_df function deals with concatenated column names, reorders columns and appends numeric suffixes
#' ## End(Not run)
#' @export get_renamed_df
#' @importFrom dplyr select_if
#' @importFrom stringr str_match
get_renamed_df <- function(df) {

  df1 <- df %>%
    dplyr::select_if(~!all(is.na(.) | . == 0))

  # Get indices of columns containing "hcbh_hdptf"
  Hcbh_Hdptf_indices <- grep("Hcbh[0-9]_Hdptf[0-9]", names(df1), value = FALSE)

  # Create two dataframes
  df1_Hcbh_Hdptf <- df1[, Hcbh_Hdptf_indices, drop = FALSE]
  df1_no_Hcbh_Hdptf <- df1[, -Hcbh_Hdptf_indices, drop = FALSE]

  # Get column names
  col_names <- names(df1_no_Hcbh_Hdptf)

  # Extract prefixes and suffixes
  col_parts <- str_match(col_names, "(\\D+)(\\d*)")
  prefixes <- col_parts[, 2]
  suffixes <- col_parts[, 3]

  # Find which column names have numeric suffix
  has_numeric_suffix <- suffixes != ""

  # Define a function to rename the columns with numeric suffixes
  rename_cols <- function(prefix) {
    # Find column names with the prefix
    cols <- col_names[prefixes == prefix & has_numeric_suffix]

    # Extract the numeric suffixes and convert to integers
    nums <- as.integer(str_extract(cols, "\\d+$"))

    # Generate the new names
    new_names <- paste0(prefix, seq_along(nums))

    # Return a named vector of new names
    setNames(new_names, cols)
  }

  # Apply the function to each prefix and unlist the results
  new_names <- unlist(lapply(unique(prefixes[has_numeric_suffix]), rename_cols))

  # Rename the columns in the dataframe
  names(df1_no_Hcbh_Hdptf) <- ifelse(names(df1_no_Hcbh_Hdptf) %in% names(new_names), new_names[names(df1_no_Hcbh_Hdptf)], names(df1_no_Hcbh_Hdptf))

  #########################################3
  # Special cases with double prefix and suffix
  # Initialize count to 1
  count <- 1

  # Create a character vector to hold new column names
  new_col_names <- character(length(names(df1_Hcbh_Hdptf)))

  # Iterate over all column names
  for (i in seq_along(names(df1_Hcbh_Hdptf))) {

    # Get the current column name
    col <- names(df1_Hcbh_Hdptf)[i]

    # Split by underscore
    parts <- unlist(strsplit(col, "_"))

    # Get prefixes (strip numeric suffix)
    prefix1 <- sub("[0-9]+$", "", parts[1])
    prefix2 <- sub("[0-9]+$", "", parts[2])

    # Construct new column name using the counter
    new_col_name <- paste0(prefix1, count, "_", prefix2, count)

    # Store the new column name
    new_col_names[i] <- new_col_name

    # Increment counter
    count <- count + 1
  }

  # Update the column names of the dataframe
  names(df1_Hcbh_Hdptf) <- new_col_names

  merged_df<-data.frame(df1_no_Hcbh_Hdptf,df1_Hcbh_Hdptf)

  #################################################

  return(merged_df)
}
