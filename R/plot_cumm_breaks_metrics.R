#' Plot of tree profiles with the canopy base height (CBH) based on the breaking point method
#' @description This function plots the canopy base height (CBH) based on breaking point over the cummulative LAD values and gives the LAD percentage below and above that breaking point
#' @usage get_plots_cumm(LAD_profiles, cummulative_LAD)
#' @param LAD_profiles original tree Leaf Area Density (LAD) profile (output of [lad.profile()] function from leafR package).
#' An object of the class text
#' @param cummulative_LAD tree metrics derived from using breaking points on cummulative LAD (output of [get_cum_break()] function).
#' An object of the class text
#' @return A plot of the CBH and LAD percentage below and above the CBH
#' @author Olga Viedma, Carlos Silva and JM Moreno
#'
#' @examples
#' \dontrun{
#' library(ggplot2)
#' library(dplyr)
#'
#' # LAD profiles derived from normalized ALS data after applying [lad.profile()] function
#' LAD_profiles <- read.table(system.file("extdata", "LAD_profiles.txt", package = "LadderFuelsR"),
#' header = TRUE)
#' LAD_profiles$treeID <- factor(LAD_profiles$treeID)
#'
#' # Load the cummulative_LAD object
#' if (interactive()) {
#' cummulative_LAD <- get_cum_break()
#' LadderFuelsR::cummulative_LAD$treeID <- factor(LadderFuelsR::cummulative_LAD$treeID)
#'
#' # Tree metrics derived from get_cum_break() function
#' cummulative_LAD$treeID <- factor(cummulative_LAD$treeID)
#'
#' # Generate cumulative LAD plots
#' plots_trees_cumlad <- get_plots_cumm(LAD_profiles, cummulative_LAD)
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
get_plots_cumm <- function(LAD_profiles, cummulative_LAD) {

  df_orig<- LAD_profiles
  df_effective1 <- cummulative_LAD

  #  Ensure treeID columns are factors
  df_orig$treeID <- factor(df_orig$treeID)
  df_effective1$treeID <- factor(df_effective1$treeID)
  treeID<-factor(df_orig$treeID)

  #Remove duplicates and columns with all NA values
  df_effective1 <- df_effective1 %>% dplyr::distinct(treeID, .keep_all = TRUE)

  df_effective1 <- df_effective1[, !apply(is.na(df_effective1), 2, all)]

  trees_name1<- as.character(df_effective1$treeID)
  trees_name2<- factor(unique(trees_name1))

  #Print the number of unique trees
  #print(paste("Number of unique trees: ", length(trees_name2)))

  #Initialize the list for the plots with annotations
  plot_with_annotations_list <- list()

  for (i in levels(trees_name2)) {

    #  Filter data for each tree
    tree_data <- df_orig[df_orig$treeID == i, ]
    df_effective1_tree <- df_effective1[df_effective1$treeID == i, ]

    #  Print the tree being processed
    #print(paste("Processing tree: ", i))

    height <- tree_data$height
    lad <- tree_data$lad

    min_y <- min(tree_data$lad, na.rm = TRUE)
    max_y <- max(tree_data$lad, na.rm = TRUE)


    if (nrow(df_effective1_tree) > 0 && any(!is.na(df_effective1_tree$Hcbh_bp))) {

      Hcbh_bp <- as.numeric(as.character(df_effective1_tree$Hcbh_bp))
      Hcbh_bp <- round(Hcbh_bp, 1)

      hcbh_cols <- df_effective1_tree[, grepl("^Hcbh_bp", names(df_effective1_tree))]
      # columns_to_check <- names(hcbh_cols)[!(names(hcbh_cols) %in% c("Hcbh1", "Hcbh_bp"))]
      is_only_Hcbh <- sum(!is.na(as.numeric(unlist(hcbh_cols)))) >=1


      #print(paste("is_only_Hcbh for tree ", i, ": ", is_only_Hcbh))

      if (is_only_Hcbh) {
        bp2a <- ggplot()  #  Initialize bp2 inside the loop

        tryCatch({
          x = height
          y = lad
          bp2a <- bp2a +
            geom_line(data = tree_data, aes(x = height, y = lad), color = "black", size = 0.5) +
            geom_point(data = tree_data, aes(x = height, y = lad), color = "black", size = 1.5)

          line_data_bp <- data.frame(x = c(Hcbh_bp, Hcbh_bp), y = c(min_y, max_y))
          bp2a <- bp2a +
            geom_path(data = line_data_bp, aes(x = x, y = y), color = "dark green", linewidth = 1, linetype = "solid", na.rm = TRUE)


          bp2a <- bp2a +
            theme_bw() +
            theme(
              axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 1, color = "black", size = 14, family = "sans"),
              axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, color = "black", size = 14, family = "sans"),
              axis.title.x = element_text(size = 14, family = "sans", color = "black", face = "bold"),
              axis.title.y = element_text(size = 14, family = "sans", color = "black", face = "bold")
            )

          bp2a <- bp2a + xlab("Height") +
            ylab("LAD") +
            ggtitle(paste0("tree_", i)) +
            coord_flip()

          below_hcbhbp <- as.numeric(as.character(df_effective1_tree$below_hcbhbp))
          above_hcbhbp <- as.numeric(as.character(df_effective1_tree$above_hcbhbp))
          label_below <- round(below_hcbhbp, 1)
          label_above <- round(above_hcbhbp, 1)

          below_hcbhbp1 <- as.character(label_below)
          above_hcbhbp1 <- as.character(label_above)

          Hcbh_bp1 <- as.numeric(as.character(df_effective1_tree$Hcbh_bp))
          label_cbh_bp <- round(Hcbh_bp1, 1)
          label_cbh_bp1<- as.character(label_cbh_bp)



          bp2_annotations1 <- bp2a +
            geom_text(data = data.frame(Hcbh_bp = Hcbh_bp, min_y = min_y, below_hcbhbp1 = below_hcbhbp1),
                      aes(x = Hcbh_bp, y = min_y, label = paste0(below_hcbhbp1,"","%")),
                      color = "black", hjust = -1, vjust = 2, size = 6) +
            geom_text(data = data.frame(Hcbh_bp = Hcbh_bp, max_y = max_y, above_hcbhbp1 = above_hcbhbp1),
                      aes(x = Hcbh_bp, y = max_y, label = paste0(above_hcbhbp1,"","%")),
                      color = "black", hjust = 1, vjust = -1, size = 6) +
            geom_text(data = data.frame(Hcbh_bp = Hcbh_bp, max_y = max_y, label_cbh_bp1 = label_cbh_bp1),
                      aes(x = Hcbh_bp, y = max_y, label =  paste0("CBH_bp ="," ",label_cbh_bp1," m")),
                      color = "black", hjust = 1, vjust = 2, size = 6)


          plot_with_annotations_list[[i]] <- bp2_annotations1  #  Store plot with annotations separately
          #print(paste("Plot for tree ", i, " created successfully"))

        }, error = function(e) {
          #print(paste("Error occurred for tree:", i))
          #print(e)
        })
      }
    }
  }

  return(plot_with_annotations_list)
}


