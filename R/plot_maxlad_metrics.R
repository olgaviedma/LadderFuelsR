#' Plots of fuel layers with LAD percentage > 5 and the canopy base height (CBH) based on the maximum LAD percentage
#' @description
#' This function plots the CBH of the fuel layer with the maximum LAD percentage and other fuel layers with LAD percentage greater than 5.
#' @usage
#' get_plots_cbh_LAD(LAD_profiles, effective_LAD)
#' @param LAD_profiles
#' Original tree Leaf Area Density (LAD) profile (output of [lad.profile()] function in the \emph{leafR} package).
#' An object of the class text.
#' @param effective_LAD
#' Tree metrics with gaps (distances), fuel base heights, and depths of fuel layers with LAD percentage greater than 5
#' (output of [get_layers_lad()] function).
#' An object of the class text.
#' @return
#' A plot drawing the CBH of the fuel layer with the maximum LAD percentage and other fuel layers with LAD percentage greater than 5.
#' @author
#' Olga Viedma, Carlos Silva and JM Moreno
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
#' # Load the effective_LAD object
#' if (interactive()) {
#' effective_LAD <- get_layers_lad()
#' LadderFuelsR::effective_LAD$treeID <- factor(LadderFuelsR::effective_LAD$treeID)
#'
#' trees_name1 <- as.character(effective_LAD$treeID)
#' trees_name2 <- factor(unique(trees_name1))
#'
#' # Generate plots for fuels LAD metrics
#' plots_trees_LAD <- get_plots_cbh_LAD(LAD_profiles, effective_LAD)
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
get_plots_cbh_LAD <- function (LAD_profiles, effective_LAD) {

  df_orig <- LAD_profiles

  df_orig$treeID <- factor(df_orig$treeID)
  trees_name1a <- as.character(df_orig$treeID)
  trees_name3 <- factor(unique(trees_name1a))
  treeID<-df_orig$treeID

  plot_with_annotations_list <- list()

  for (i in levels(trees_name3)) {

    tree_data <- df_orig %>%
      dplyr::filter(treeID == i) %>%
      dplyr::mutate(lad = as.numeric(lad)) %>%
      dplyr::filter(!is.na(lad))

    height <- tree_data$height
    lad <- tree_data$lad

    df_effective1 <- effective_LAD %>% dplyr::filter(treeID == i)

    CBH_1 <- round(as.numeric(as.character(df_effective1$Hcbh1)), 1)
    CBH_2 <- round(as.numeric(as.character(df_effective1$Hcbh2)), 1)
    CBH_3 <- round(as.numeric(as.character(df_effective1$Hcbh3)), 1)
    CBH_4 <- round(as.numeric(as.character(df_effective1$Hcbh4)), 1)
    CBH_5 <- round(as.numeric(as.character(df_effective1$Hcbh5)), 1)
    CBH_6 <- round(as.numeric(as.character(df_effective1$Hcbh6)), 1)

    depth_1 <- as.numeric(as.character(df_effective1$Hdptf1))
    depth_2 <- as.numeric(as.character(df_effective1$Hdptf2))
    depth_3 <- as.numeric(as.character(df_effective1$Hdptf3))
    depth_4 <- as.numeric(as.character(df_effective1$Hdptf4))
    depth_5 <- as.numeric(as.character(df_effective1$Hdptf5))
    depth_6 <- as.numeric(as.character(df_effective1$Hdptf6))


    min_y <- min(tree_data$lad, na.rm = TRUE)
    max_y <- max(tree_data$lad, na.rm = TRUE)

    #print(CBH_0)
    #print(depth_0)
    #print(min_y)
    #print(max_y)

    x<-tree_data$height
    y<-tree_data$lad

    tryCatch({
      bp2 <- ggplot(tree_data, aes(x = height)) +
        geom_line(aes(y = lad), color = "black", linewidth = 0.5) +
        geom_point(data = tree_data, aes(x = height, y = lad), color = "black", size = 1.5)

      if (!is.na(min_y) && !is.na(max_y)) {

        tryCatch({

          if (!any(is.na(CBH_1)) && !any(is.na(depth_1))) {
            if (CBH_1 != depth_1) {
              polygon_data_1 <- data.frame(x = c(CBH_1, CBH_1, depth_1, depth_1),
                                           y = c(min_y, max_y, max_y, min_y))
              bp2 <- bp2 +
                geom_polygon(data = polygon_data_1,
                             aes(x = x, y = y), fill = "dark green", alpha = 0.3)
            } else {
              line_data_1 <- data.frame(x = c(CBH_1, depth_1),
                                        y = c(min_y, max_y))
              bp2 <- bp2 +
                geom_path(data = line_data_1,
                          aes(x = x, y = y), color = "dark green", size = 1, linetype = "solid")
            } }
        }, error = function(e) {})


        tryCatch({

          if (!any(is.na(CBH_2)) && !any(is.na(depth_2))) {
            if (CBH_2 != depth_2) {

              polygon_data_2 <- data.frame(x = c(CBH_2, CBH_2, depth_2, depth_2),
                                           y = c(min_y, max_y, max_y, min_y))
              bp2 <- bp2 +
                geom_polygon(data = polygon_data_2,
                             aes(x = x, y = y), fill = "dark green", alpha = 0.3)
            } else {
              line_data_2 <- data.frame(x = c(CBH_2, depth_2),
                                        y = c(min_y, max_y))
              bp2 <- bp2 +
                geom_path(data = line_data_2,
                          aes(x = x, y = y), color = "dark green", size = 1, linetype = "solid")
            }}

        }, error = function(e) {})


        tryCatch({

          if (!any(is.na(CBH_3)) && !any(is.na(depth_3))) {
            if (CBH_3 != depth_3) {

              polygon_data_3 <- data.frame(x = c(CBH_3, CBH_3, depth_3, depth_3),
                                           y = c(min_y, max_y, max_y, min_y))
              bp2 <- bp2 +
                geom_polygon(data = polygon_data_3,
                             aes(x = x, y = y), fill = "dark green", alpha = 0.3)
            } else {
              line_data_3 <- data.frame(x = c(CBH_3, depth_3),
                                        y = c(min_y, max_y))
              bp2 <- bp2 +
                geom_path(data = line_data_3,
                          aes(x = x, y = y), color = "dark green", size = 1, linetype = "solid")
            }}
        }, error = function(e) {})


        tryCatch({

          if (!any(is.na(CBH_4)) && !any(is.na(depth_4))) {
            if (CBH_4 != depth_4) {

              polygon_data_4 <- data.frame(x = c(CBH_4, CBH_4, depth_4, depth_4),
                                           y = c(min_y, max_y, max_y, min_y))
              bp2 <- bp2 +
                geom_polygon(data = polygon_data_4,
                             aes(x = x, y = y), fill = "dark green", alpha = 0.3)
            } else {
              line_data_4 <- data.frame(x = c(CBH_4, depth_4),
                                        y = c(min_y, max_y))
              bp2 <- bp2 +
                geom_path(data = line_data_4,
                          aes(x = x, y = y), color = "dark green", size = 1, linetype = "solid")
            }}
        }, error = function(e) {})


        tryCatch({

          if (!any(is.na(CBH_5)) && !any(is.na(depth_5))) {
            if (CBH_5 != depth_5) {

              polygon_data_5 <- data.frame(x = c(CBH_5, CBH_5, depth_5, depth_5),
                                           y = c(min_y, max_y, max_y, min_y))
              bp2 <- bp2 +
                geom_polygon(data = polygon_data_5,
                             aes(x = x, y = y), fill = "dark green", alpha = 0.3)
            } else {
              line_data_5<- data.frame(x = c(CBH_5, depth_5),
                                       y = c(min_y, max_y))
              bp2 <- bp2 +
                geom_path(data = line_data_5,
                          aes(x = x, y = y), color = "dark green", size = 1, linetype = "solid")
            }}

        }, error = function(e) {})


        tryCatch({

          if (!any(is.na(CBH_6)) && !any(is.na(depth_6))) {
            if (CBH_6 != depth_6) {

              polygon_data_6 <- data.frame(x = c(CBH_6, CBH_6, depth_6, depth_6),
                                           y = c(min_y, max_y, max_y, min_y))
              bp2 <- bp2 +
                geom_polygon(data = polygon_data_6,
                             aes(x = x, y = y), fill = "dark green", alpha = 0.3)
            } else {
              line_data_6 <- data.frame(x = c(CBH_6, depth_6),
                                        y = c(min_y, max_y))
              bp2 <- bp2 +
                geom_path(data = line_data_6,
                          aes(x = x, y = y), color = "dark green", size = 1, linetype = "solid")
            }}
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


        Hcbh1_Hdptf1 <- as.numeric(as.character(df_effective1$Hcbh1_Hdptf1))
        label_Hcbh1_Hdptf1 <- round(Hcbh1_Hdptf1, 1)
        Hcbh1_Hdptf1a <- paste0(as.character(label_Hcbh1_Hdptf1),"","%")

        Hcbh2_Hdptf2 <- as.numeric(as.character(df_effective1$Hcbh2_Hdptf2))
        label_Hcbh2_Hdptf2 <- round(Hcbh2_Hdptf2, 1)
        Hcbh2_Hdptf2a <- paste0(as.character(label_Hcbh2_Hdptf2),"","%")

        Hcbh3_Hdptf3 <- as.numeric(as.character(df_effective1$Hcbh3_Hdptf3))
        label_Hcbh3_Hdptf3 <- round(Hcbh3_Hdptf3, 1)
        Hcbh3_Hdptf3a <-  paste0(as.character(label_Hcbh3_Hdptf3),"","%")

        Hcbh4_Hdptf4 <- as.numeric(as.character(df_effective1$Hcbh4_Hdptf4))
        label_Hcbh4_Hdptf4 <- round(Hcbh4_Hdptf4, 1)
        Hcbh4_Hdptf4a <- paste0(as.character(label_Hcbh4_Hdptf4),"","%")

        Hcbh5_Hdptf5 <- as.numeric(as.character(df_effective1$Hcbh5_Hdptf5))
        label_Hcbh5_Hdptf5 <- round(Hcbh5_Hdptf5, 1)
        Hcbh5_Hdptf5a <-  paste0(as.character(label_Hcbh5_Hdptf5),"","%")

        Hcbh6_Hdptf6 <- as.numeric(as.character(df_effective1$Hcbh6_Hdptf6))
        label_Hcbh6_Hdptf6 <- round(Hcbh6_Hdptf6, 1)
        Hcbh6_Hdptf6a <- paste0(as.character(label_Hcbh6_Hdptf6),"","%")

        CBH1_label<- paste0("CBH ="," ",CBH_1,"m")
        Depth1_label<- paste0("Depth ="," ",depth_1,"m")
        CBH2_label<- paste0("CBH ="," ",CBH_2,"m")
        Depth2_label<- paste0("Depth ="," ",depth_2,"m")
        CBH3_label<- paste0("CBH ="," ",CBH_3,"m")
        Depth3_label<- paste0("Depth ="," ",depth_3,"m")
        CBH4_label<- paste0("CBH ="," ",CBH_4,"m")
        Depth4_label<- paste0("Depth ="," ",depth_4,"m")
        CBH5_label<- paste0("CBH ="," ",CBH_5,"m")
        Depth5_label<- paste0("Depth ="," ",depth_5,"m")
        CBH6_label<- paste0("CBH ="," ",CBH_6,"m")
        Depth6_label<- paste0("Depth ="," ",depth_6,"m")


        bp2_annotations <- bp2


        if (any(!is.na(CBH_1)) && any(!is.na(Hcbh1_Hdptf1a))) {

          y_1 = min_y
          bp2_annotations <- bp2_annotations + geom_text(data = data.frame(CBH_1 = CBH_1, y_1 = min_y , Hcbh1_Hdptf1a = Hcbh1_Hdptf1a),
                                                         aes(x = CBH_1,y = y_1, label = Hcbh1_Hdptf1a),
                                                         color = "black", hjust = -2.5, vjust = 0, size = 5)
          y_1 = max_y
          bp2_annotations <- bp2_annotations + geom_text(data = data.frame(CBH_1 = CBH_1, y_1 = max_y , CBH1_label = CBH1_label),
                                                         aes(x = CBH_1,y = y_1, label = CBH1_label),
                                                         color = "black", hjust = 1, vjust = 0, size = 5)
          y_1 = max_y
          bp2_annotations <- bp2_annotations + geom_text(data = data.frame(depth_1 = depth_1, y_1 = max_y , Depth1_label = Depth1_label),
                                                         aes(x = depth_1,y = y_1, label = Depth1_label),
                                                         color = "black", hjust = 2, vjust = 1, size = 5)

        }
        if (any(!is.na(CBH_2)) && any(!is.na(Hcbh2_Hdptf2a))){

          y_2 = min_y

          bp2_annotations <- bp2_annotations + geom_text(data = data.frame(CBH_2 = CBH_2, y_2 = min_y , Hcbh2_Hdptf2a = Hcbh2_Hdptf2a),
                                                         aes(x = CBH_2,y = y_2, label = Hcbh2_Hdptf2a),
                                                         color = "black", hjust = -2.5, vjust = 0, size = 5)
          y_2 = max_y
          bp2_annotations <- bp2_annotations + geom_text(data = data.frame(CBH_2 = CBH_2, y_2 = max_y , CBH2_label = CBH2_label),
                                                         aes(x = CBH_2,y = y_2, label = CBH2_label),
                                                         color = "black", hjust = 1, vjust = 0, size = 5)
          y_2 = max_y
          bp2_annotations <- bp2_annotations + geom_text(data = data.frame(depth_2 = depth_2, y_2 = max_y , Depth2_label = Depth2_label),
                                                         aes(x = depth_2,y = y_2, label = Depth2_label),
                                                         color = "black", hjust = 2, vjust = 1, size = 5)
        }

        if (any(!is.na(CBH_3)) && any(!is.na(Hcbh3_Hdptf3a))) {

          y_3 = min_y
          bp2_annotations <- bp2_annotations + geom_text(data = data.frame(CBH_3 = CBH_3, y_3 = min_y , Hcbh3_Hdptf3a = Hcbh3_Hdptf3a),
                                                         aes(x = CBH_3,y = y_3, label = Hcbh3_Hdptf3a),
                                                         color = "black", hjust = -2.5, vjust = 0, size = 5)
          y_3 = max_y
          bp2_annotations <- bp2_annotations + geom_text(data = data.frame(CBH_3 = CBH_3, y_3 = max_y , CBH3_label = CBH3_label),
                                                         aes(x = CBH_3,y = y_3, label = CBH3_label),
                                                         color = "black", hjust = 1, vjust = 0, size = 5)
          y_3 = max_y
          bp2_annotations <- bp2_annotations + geom_text(data = data.frame(depth_3 = depth_3, y_3 = max_y , Depth3_label = Depth3_label),
                                                         aes(x = depth_3,y = y_3, label = Depth3_label),
                                                         color = "black", hjust =2, vjust = 1, size = 5)
        }


        if (any(!is.na(CBH_4)) && any(!is.na(Hcbh4_Hdptf4a))) {

          y_4 = min_y
          bp2_annotations <- bp2_annotations +geom_text(data = data.frame(CBH_4 = CBH_4, y_4 = min_y , Hcbh4_Hdptf4a = Hcbh4_Hdptf4a),
                                                        aes(x = CBH_4,y = y_4, label = Hcbh4_Hdptf4a),
                                                        color = "black", hjust =-2.5, vjust = 0, size = 5)
          y_4 = max_y
          bp2_annotations <- bp2_annotations + geom_text(data = data.frame(CBH_4 = CBH_4, y_4 = max_y , CBH4_label = CBH4_label),
                                                         aes(x = CBH_4,y = y_4, label = CBH4_label),
                                                         color = "black", hjust = 1, vjust = 0, size = 5)
          y_4 = max_y
          bp2_annotations <- bp2_annotations + geom_text(data = data.frame(depth_4 = depth_4, y_4 = max_y , Depth4_label = Depth4_label),
                                                         aes(x = depth_4,y = y_4, label = Depth4_label),
                                                         color = "black", hjust = 2, vjust = 1, size = 5)
        }


        if (any(!is.na(CBH_5)) && any(!is.na(Hcbh5_Hdptf5a))) {

          y_5 = min_y
          bp2_annotations <- bp2_annotations + geom_text(data = data.frame(CBH_5 = CBH_5, y_5 = min_y, Hcbh5_Hdptf5a = Hcbh5_Hdptf5a),
                                                         aes(x = CBH_5,y = y_5, label = Hcbh5_Hdptf5a),
                                                         color = "black", hjust = -2.5, vjust = 0, size = 5)
          y_5 = max_y
          bp2_annotations <- bp2_annotations + geom_text(data = data.frame(CBH_5 = CBH_5, y_5 = max_y , CBH5_label = CBH5_label),
                                                         aes(x = CBH_5,y = y_5, label = CBH5_label),
                                                         color = "black", hjust = 1, vjust = 0, size = 5)
          y_5 = max_y
          bp2_annotations <- bp2_annotations + geom_text(data = data.frame(depth_5 = depth_5, y_5 = max_y , Depth5_label = Depth5_label),
                                                         aes(x = depth_5,y = y_5, label = Depth5_label),
                                                         color = "black", hjust = 2, vjust = 1, size = 5)
        }

        if (any(!is.na(CBH_6)) && any(!is.na(Hcbh6_Hdptf6a))) {

          y_6 = min_y
          bp2_annotations <- bp2_annotations +geom_text(data = data.frame(CBH_6 = CBH_6, y_6 = min_y, Hcbh6_Hdptf6a = Hcbh6_Hdptf6a),
                                                        aes(x = CBH_6,y = y_6, label = Hcbh6_Hdptf6a),
                                                        color = "black", hjust = -2, vjust = 0, size = 5)
          y_6 = max_y
          bp2_annotations <- bp2_annotations + geom_text(data = data.frame(CBH_6 = CBH_6, y_6 = max_y , CBH6_label = CBH6_label),
                                                         aes(x = CBH_6,y = y_6, label = CBH6_label),
                                                         color = "black", hjust = 1, vjust = 0, size = 5)
          y_6 = max_y
          bp2_annotations <- bp2_annotations + geom_text(data = data.frame(depth_6 = depth_6, y_6 = max_y , Depth6_label = Depth6_label),
                                                         aes(x = depth_6,y = y_6, label = Depth6_label),
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

