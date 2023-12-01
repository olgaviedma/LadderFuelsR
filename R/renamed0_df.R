#'
#' Rename and reorder columns (I)
#'
#' @description This function reorders columns and appends numeric suffixes. Don´t run it. It is an internal function.
#'
#' @usage get_renamed0_df (df)
#'
#' @param df internal data frame derived from [get_real_depths()] function
#'
#' @examples
#' ## Not run:
#' library(dplyr)
#' # get_renamed0_df function reorders columns and appends numeric suffixes
#' ## End(Not run)
#' @export get_renamed0_df
#' @importFrom dplyr select_if group_by summarise summarize mutate arrange rename rename_with filter slice slice_tail ungroup distinct
#' @importFrom stringr str_match
get_renamed0_df <- function(df) {

  df1 <- df %>%
    dplyr::select_if(~!all(is.na(.) | . == 0))

  # Exclude columns with prefixes: last_, max_, treeID
  cols_to_exclude <- grep("^(last_|max_|treeID)", names(df1), value = TRUE)
  df1 <- df1[ , !(names(df1) %in% cols_to_exclude)]

  # Extract unique prefixes
  prefixes <- unique(gsub("([a-zA-Z]+).*", "\\1", names(df1)))

  # Rename the columns based on the extracted prefixes
  for (prefix in prefixes) {
    # Identify columns with the current prefix
    cols <- grep(paste0("^", prefix), names(df1))

    # Generate new column names with consecutive suffixes
    new_names <- paste0(prefix, 1:length(cols))

    # Assign new names to the columns
    names(df1)[cols] <- new_names
  }

  # Select columns with prefixes: last_, max_, treeID
  cols_to_append <- grep("^(last_|max_|treeID)", names(df1), value = TRUE)
  append_df <- df1[ , cols_to_append]

  # Append the columns to df1
  df1 <- cbind(df1, append_df)

  return(df1)
}
