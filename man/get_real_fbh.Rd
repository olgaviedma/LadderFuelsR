% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/corrected_base_heights.R
\name{get_real_fbh}
\alias{get_real_fbh}
\title{Fuels base height after removing distances = 1 m}
\usage{
get_real_fbh(depth_metrics,verbose=TRUE)
}
\arguments{
\item{depth_metrics}{Tree metrics with gaps (distances), fuel base heights, and depths
(output of [get_depths()] function). An object of the class text.}

\item{verbose}{Logical, indicating whether to display informational messages (default is TRUE).}
}
\value{
A data frame giving the first "base height" from those consecutive ones separated by a distance = 1.
This value replaces the values of the next base heights if they are separated by a distance = 1.
}
\description{
This function recalculates fuels base height after removing distances = 1 m,
keeping the first "base height" from those consecutive ones separated by a distance = 1.
}
\details{
# List of tree metrics:
\itemize{
  \item treeID: tree ID with strings and numeric values
  \item treeID1: tree ID with only numeric values
  \item dist: Distance between consecutive fuel layers (m)
  \item Hdist - Height of the distance between consecutive fuel layers (m)
  \item Hcbh - Height of the first base height from those ones separated by a distance = 1.
  \item depth - Depth of fuel layers (m)
  \item Hdepth - Height of the depth of fuel layers (m)
  \item max_height - Maximum height of the tree profile
}
}
\examples{
library(magrittr)
library(dplyr)
#Before running this example, make sure to run get_depths()
if (interactive()) {
depth_metrics <- get_depths()
LadderFuelsR::depth_metrics$treeID <- factor(LadderFuelsR::depth_metrics$treeID)

trees_name1 <- as.character(depth_metrics$treeID)
trees_name2 <- factor(unique(trees_name1))

fbh_corr_list <- list()

for (i in levels(trees_name2)){
# Filter data for each tree
tree3 <- depth_metrics |> dplyr::filter(treeID == i)
# Get real fbh for each tree
fbh_corr <- get_real_fbh(tree3, verbose=TRUE)
# Store fbh values in a list
fbh_corr_list[[i]] <- fbh_corr
}

# Combine fbh values for all trees
effective_fbh <- dplyr::bind_rows(fbh_corr_list)
effective_fbh$treeID <- factor(effective_fbh$treeID)
}
}
\seealso{
\code{\link{get_depths}}
}
\author{
Olga Viedma, Carlos Silva and JM Moreno
}
