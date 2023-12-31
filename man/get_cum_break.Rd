% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cummLAD_breaks_metrics.R
\name{get_cum_break}
\alias{get_cum_break}
\title{CBH and the LAD percentage below and above the CBH using the breaking point method}
\usage{
get_cum_break(LAD_profiles, effective_distances, verbose=TRUE)
}
\arguments{
\item{LAD_profiles}{Original tree Leaf Area Density (LAD) profile (output of [lad.profile()] function in the \emph{leafR} package).
An object of the class data frame.}

\item{effective_distances}{Tree metrics of fuel layers separated by distances greater than 1 m (output of [get_effective_gap()] function).
An object of the class data frame.}

\item{verbose}{Logical, indicating whether to display informational messages (default is TRUE).}
}
\value{
A data frame identifying the Canopy Base Height (CBH) of the vertical tree profile (VTP) based on the breaking point identified by the segmented regression model,
and the percentage of LAD values below and above the identified CBH or breaking point.
}
\description{
This function calculates the canopy base height (CBH) of the vertical tree profile (VTP) using a segmented regression model fitted to
the cumulative LAD values as a function of height.The function also calculates the percentage of LAD values below and above the identified CBH or breaking point.
}
\details{
\itemize{
  \item treeID: tree ID with strings and numeric values
  \item treeID1: tree ID with only numeric values
  \item Hdist: Height of the distance between the ground and the CBH or breaking point (m)
  \item Hcbh: Height of the CBH or breaking point (m)
  \item Hcbh_bp: Height of the CBH or breaking point (m)
  \item effdist: Distance between the ground and the CBH or breaking point (m)
  \item dptf: Depth of the CBH or breaking point (m)
  \item Hdptf: Height of the depth of the CBH or breaking point (m)
  \item below_hcbhbp: Percentage of LAD values below the CBH or breaking point
  \item above_hcbhbp: Percentage of LAD values above the CBH or breaking point
  \item maxlad_Hcbh: Height of the CBH or breaking point with maximum LAD percentage
  \item max_Hcbh: Height of the CBH or breaking point at maximum height distance
  \item last_Hcbh: Height of the CBH or breaking point at the last distance
  \item maxlad_: Values of distance and fuel depth and their corresponding heights for the CBH or breaking point with maximum LAD percentage
  \item max_: Values of distance and fuel depth and their corresponding heights for the CBH or breaking point at maximum height distance
  \item last_: Values of distance and fuel depth and their corresponding heights for the CBH or breaking point at the last distance
  \item cumlad: Cumulative LAD values at the CBH or breaking point
  \item max_height: Maximum height of the tree profile
}
}
\examples{
library(magrittr)
library(segmented)
library(gdata)
library(dplyr)

# LAD profiles derived from normalized ALS data after applying [lad.profile()] function
LAD_profiles <- read.table(system.file("extdata", "LAD_profiles.txt", package = "LadderFuelsR"),
header = TRUE)
LAD_profiles$treeID <- factor(LAD_profiles$treeID)

# Before running this example, make sure to run get_effective_gap().
if (interactive()) {
effective_distances <- get_effective_gap()
LadderFuelsR::effective_distances$treeID <- factor(LadderFuelsR::effective_distances$treeID)

trees_name1 <- as.character(effective_distances$treeID)
trees_name2 <- factor(unique(trees_name1))

cum_LAD_metrics_list <- list()

for (i in levels(trees_name2)) {
# Filter data for each tree
tree1 <- LAD_profiles |> dplyr::filter(treeID == i)
tree2 <- effective_distances |> dplyr::filter(treeID == i)

# Get cumulative LAD metrics for each tree
cum_LAD_metrics <- get_cum_break(tree1, tree2)
cum_LAD_metrics_list[[i]] <- cum_LAD_metrics
}

# Combine the individual data frames
cummulative_LAD <- dplyr::bind_rows(cum_LAD_metrics_list)
}
}
\seealso{
\code{\link{get_effective_gap}}
}
\author{
Olga Viedma, Carlos Silva and JM Moreno
}
