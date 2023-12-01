#' LAD Profiles Dataset
#'
#' Description of the LAD Profiles dataset.
#'
#' @format A data frame with columns: height, lad, treeID.
#'
#' @examples
#' data(LAD_profiles)
#'
#' @export
LAD_profiles <- read.table(system.file("extdata", "LAD_profiles.txt", package = "LadderFuelsR"), header = TRUE)

