#' Plots the Crown Base Height (CBH) based on the last distance criterium
#' @description
#' This function plots the CBH of a segmented tree based on the fuel layer located at the last distance.
#' @usage
#' get_plots_cbh_lastdist(LAD_profiles, cbh_metrics, min_height=1.5)
#' @param LAD_profiles
#' Original tree Leaf Area Density (LAD) profile (output of [lad.profile()] function in the \emph{leafR} package).
#' An object of the class text.
#' @param cbh_metrics
#' CBH metrics based on three criteria: maximum LAD percentage, maximum distance and last distance.
#' (output of [get_cbh_metrics()] function).
#' An object of the class text.
#' @param min_height Numeric value for the actual minimum base height (in meters).
#' @return
#' A plot drawing the Crown Base Height (CBH) of the fuel layer located at the last distance.
#' @author Olga Viedma, Carlos Silva, JM Moreno and A.T. Hudak
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
#' # Before running this example, make sure to run get_cbh_metrics().
#' if (interactive()) {
#' cbh_metrics <- get_cbh_metrics()
#' LadderFuelsR::cbh_metrics$treeID <- factor(LadderFuelsR::cbh_metrics$treeID)
#'
#' trees_name1 <- as.character(cbh_metrics$treeID)
#' trees_name2 <- factor(unique(trees_name1))
#'
#' # Generate plots for CBH based on the fuel layer at the last distance
#' plots_cbh_lastdist <- get_plots_cbh_lastdist(LAD_profiles, cbh_metrics, min_height=1.5)
#' }
#' @importFrom dplyr select_if group_by summarise summarize mutate arrange rename rename_with filter slice slice_tail ungroup distinct
#' across matches row_number all_of vars bind_cols case_when left_join mutate if_else lag n_distinct
#' @importFrom segmented segmented seg.control
#' @importFrom magrittr %>%
#' @importFrom stats ave dist lm na.omit predict quantile setNames smooth.spline
#' @importFrom utils tail
#' @importFrom tidyselect starts_with everything one_of
#' @importFrom stringr str_extract str_match str_detect str_remove_all
#' @importFrom tibble tibble
#' @importFrom tidyr pivot_longer fill pivot_wider replace_na
#' @importFrom gdata startsWith
#' @importFrom ggplot2 aes geom_line geom_path geom_point geom_polygon geom_text geom_vline ggtitle coord_flip theme_bw
#' theme element_text xlab ylab ggplot xlim
#' @seealso \code{\link{get_cbh_metrics}}
#' @export
get_plots_cbh_lastdist <- function (LAD_profiles, cbh_metrics,min_height=1.5) {

  df_orig <- LAD_profiles


  if(min_height==0){
    min_height <-0.5

    # Ensure the column starts with a negative value
    if (df_orig$height[1] < min_height) {
      # Calculate the shift value
      shift_value <- abs(df_orig$height[1])

      # Adjust the column to start from 0
      df_orig$height <- df_orig$height + shift_value
    }


    # Ensure the column starts with a negative value
    if (df_orig$height[1] > min_height) {
      # Calculate the shift value
      shift_value1 <- abs(df_orig$height[1])

      # Adjust the column to start from 0
      df_orig$height <- df_orig$height - shift_value1
    }
  }


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

    df_effective1 <- cbh_metrics %>% dplyr::filter(treeID == i)

    last_CBH <- round(as.numeric(as.character(df_effective1$last_Hcbh)), 1)
    last_Hdepth <- as.numeric(as.character(df_effective1$last_Hdptf))

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

          if (!any(is.na(last_CBH)) && !any(is.na(last_Hdepth))) {
            if (last_CBH != last_Hdepth) {
              polygon_data_1 <- data.frame(x = c(last_CBH, last_CBH, last_Hdepth, last_Hdepth),
                                           y = c(min_y, max_y, max_y, min_y))
              bp2 <- bp2 +
                geom_polygon(data = polygon_data_1,
                             aes(x = x, y = y), fill = "dark green", alpha = 0.3)
            } else {
              line_data_1 <- data.frame(x = c(last_CBH, last_Hdepth),
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


        Hcbh1_Hdptf1 <- as.numeric(as.character(df_effective1$last_lad))
        label_Hcbh1_Hdptf1 <- round(Hcbh1_Hdptf1, 1)
        Hcbh1_Hdptf1a <- paste0(as.character(label_Hcbh1_Hdptf1),"","%")


        CBH1_label<- paste0("CBH ="," ",last_CBH,"m")
        Depth1_label<- paste0("Depth ="," ",last_Hdepth,"m")


        bp2_annotations <- bp2


        if (any(!is.na(last_CBH)) && any(!is.na(Hcbh1_Hdptf1a))) {

          y_1 = min_y
          bp2_annotations <- bp2_annotations + geom_text(data = data.frame(last_CBH = last_CBH, y_1 = min_y , Hcbh1_Hdptf1a = Hcbh1_Hdptf1a),
                                                         aes(x = last_CBH,y = y_1, label = Hcbh1_Hdptf1a),
                                                         color = "black", hjust = -2.5, vjust = 0, size = 5)
          y_1 = max_y
          bp2_annotations <- bp2_annotations + geom_text(data = data.frame(last_CBH = last_CBH, y_1 = max_y , CBH1_label = CBH1_label),
                                                         aes(x = last_CBH,y = y_1, label = CBH1_label),
                                                         color = "black", hjust = 1, vjust = 0, size = 5)
          y_1 = max_y
          bp2_annotations <- bp2_annotations + geom_text(data = data.frame(last_Hdepth = last_Hdepth, y_1 = max_y , Depth1_label = Depth1_label),
                                                         aes(x = last_Hdepth,y = y_1, label = Depth1_label),
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

