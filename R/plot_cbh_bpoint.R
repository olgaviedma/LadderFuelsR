#' Plot of the Canopy Base Height (CBH) based on the breaking point method
#' @description This function plots the canopy base height (CBH) based on breaking point over the cummulative LAD values and gives the LAD percentage
#' of the canopy layer
#' @usage get_plots_cbh_bp(LAD_profiles, cummulative_LAD)
#' @param LAD_profiles original tree Leaf Area Density (LAD) profile (output of [lad.profile()] function from leafR package).
#' An object of the class text.
#' @param cummulative_LAD tree metrics derived from using breaking points on cummulative LAD (output of [get_cum_break()] function).
#' An object of the class text.
#' @return A plot of the Canopy Base Height (CBH) based on the breaking point method and Leaf Area Density (LAD) percentage of the canopy layer.
#' @author Olga Viedma, Carlos Silva and JM Moreno
#'
#' @examples
#' library(ggplot2)
#' library(dplyr)
#'
#' # LAD profiles derived from normalized ALS data after applying [lad.profile()] function
#' LAD_profiles <- read.table(system.file("extdata", "LAD_profiles.txt", package = "LadderFuelsR"),
#' header = TRUE)
#' LAD_profiles$treeID <- factor(LAD_profiles$treeID)
#'
#' # Before running this example, make sure to run get_cum_break().
#' if (interactive()) {
#' cummulative_LAD <- get_cum_break()
#' LadderFuelsR::cummulative_LAD$treeID <- factor(LadderFuelsR::cummulative_LAD$treeID)
#'
#' # Generate cumulative LAD plots
#' plots_cbh_bp <- get_plots_cbh_bp(LAD_profiles, cummulative_LAD)
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
#' @seealso \code{\link{get_cum_break}}
#' @export
get_plots_cbh_bp <- function(LAD_profiles, cummulative_LAD) {

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


  plot_with_annotations_list <- list()

  for (i in levels(trees_name2)) {

    tree_data <- df_orig %>%
      dplyr::filter(treeID == i) %>%
      dplyr::mutate(lad = as.numeric(lad)) %>%
      dplyr::filter(!is.na(lad))

    height <- tree_data$height
    lad <- tree_data$lad

    df_effective1 <- cummulative_LAD %>% dplyr::filter(treeID == i)

    bp_CBH <- round(as.numeric(as.character(df_effective1$bp_Hcbh)), 1)
    bp_Hdepth <- as.numeric(as.character(df_effective1$bp_Hdptf))

    min_y <- min(tree_data$lad, na.rm = TRUE)
    max_y <- max(tree_data$lad, na.rm = TRUE)

    x<-tree_data$height
    y<-tree_data$lad

    tryCatch({
      bp2 <- ggplot(tree_data, aes(x = height)) +
        geom_line(aes(y = lad), color = "black", linewidth = 0.5) +
        geom_point(data = tree_data, aes(x = height, y = lad), color = "black", size = 1.5)

      if (!is.na(min_y) && !is.na(max_y)) {

        tryCatch({

          if (!any(is.na(bp_CBH)) && !any(is.na(bp_Hdepth))) {
            if (bp_CBH != bp_Hdepth) {
              polygon_data_1 <- data.frame(x = c(bp_CBH, bp_CBH, bp_Hdepth, bp_Hdepth),
                                           y = c(min_y, max_y, max_y, min_y))
              bp2 <- bp2 +
                geom_polygon(data = polygon_data_1,
                             aes(x = x, y = y), fill = "dark green", alpha = 0.3)
            } else {
              line_data_1 <- data.frame(x = c(bp_CBH, bp_Hdepth),
                                        y = c(min_y, max_y))
              bp2 <- bp2 +
                geom_path(data = line_data_1,
                          aes(x = x, y = y), color = "dark green", size = 1, linetype = "solid")
            } }
        }, error = function(e) {})


        bp2 <- bp2 +
          theme_bw() +
          theme(
            axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 1, color = "black", size = 14, family = "sans"),
            axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, color = "black", size = 14, family = "sans"),
            axis.title.x = element_text(size = 14, family = "sans", color = "black", face = "bold"),
            axis.title.y = element_text(size = 14, family = "sans", color = "black", face = "bold")) +
          xlab("Height") +
          ylab("LAD") +
          ggtitle(paste0("tree_", i)) +
          coord_flip()


        Hcbh1_Hdptf1 <- as.numeric(as.character(df_effective1$bp_lad))
        label_Hcbh1_Hdptf1 <- round(Hcbh1_Hdptf1, 1)
        Hcbh1_Hdptf1a <- paste0(as.character(label_Hcbh1_Hdptf1),"","%")


        CBH1_label<- paste0("CBH ="," ",bp_CBH,"m")
        Depth1_label<- paste0("Depth ="," ",bp_Hdepth,"m")


        bp2_annotations <- bp2


        if (any(!is.na(bp_CBH)) && any(!is.na(Hcbh1_Hdptf1a))) {

          y_1 = min_y
          bp2_annotations <- bp2_annotations + geom_text(data = data.frame(bp_CBH = bp_CBH, y_1 = min_y , Hcbh1_Hdptf1a = Hcbh1_Hdptf1a),
                                                         aes(x = bp_CBH,y = y_1, label = Hcbh1_Hdptf1a),
                                                         color = "black", hjust = -2.5, vjust = 0, size = 5)
          y_1 = max_y
          bp2_annotations <- bp2_annotations + geom_text(data = data.frame(bp_CBH = bp_CBH, y_1 = max_y , CBH1_label = CBH1_label),
                                                         aes(x = bp_CBH,y = y_1, label = CBH1_label),
                                                         color = "black", hjust = 1, vjust = 0, size = 5)
          y_1 = max_y
          bp2_annotations <- bp2_annotations + geom_text(data = data.frame(bp_Hdepth = bp_Hdepth, y_1 = max_y , Depth1_label = Depth1_label),
                                                         aes(x = bp_Hdepth,y = y_1, label = Depth1_label),
                                                         color = "black", hjust = 2, vjust = 1, size = 5)

        }



        plot_with_annotations_list[[i]] <- bp2_annotations  # Store plot with annotations separately
        #print(paste("Plot for tree ", i, " created successfully"))
      }

    }, error = function(e) {
      #print(paste("Error occurred for tree:", i))
      #print(e)
    })
  }

  return(plot_with_annotations_list)  # Changed from plot_with_annotations_list to plot_list

}

