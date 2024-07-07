#' Plots of tree profiles with gaps and fuel layers base height (fbh)
#' @description This function plots gaps and fuel layers base height (fbh) in the vertical tree profile (VTP).
#' @usage get_plots_gap_fbh (LAD_profiles,gap_cbh_metrics,min_height=1.5)
#' @param LAD_profiles original tree Leaf Area Density (LAD) profile (output of [lad.profile()] function in the \emph{leafR} package.
#' An object of the class text
#' @param gap_cbh_metrics data frame with gaps (distances) and fuel base heights (output of [get_gaps_fbhs()] function).
#' An object of the class text.
#' @param min_height Numeric value for the actual minimum base height (in meters).
#' @return A plot drawing by lines the height of gaps and fuel layers bases in tiff format.
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
#' # Before running this example, make sure to run get_gaps_fbhs().
#' if (interactive()) {
#' gap_cbh_metrics <- get_gaps_fbhs()
#' LadderFuelsR::gap_cbh_metrics$treeID <- factor(LadderFuelsR::gap_cbh_metrics$treeID)
#'
#' # Generate plots for gaps and fbhs
#' plots_gaps_fbhs <- get_plots_gap_fbh(LAD_profiles, gap_cbh_metrics, min_height=1.5)
#' }
#' @importFrom dplyr select_if group_by summarise summarize mutate arrange rename rename_with filter slice slice_tail ungroup distinct
#' across matches row_number all_of vars
#' @importFrom segmented segmented seg.control
#' @importFrom magrittr %>%
#' @importFrom stats ave dist lm na.omit predict quantile setNames smooth.spline
#' @importFrom utils tail
#' @importFrom tidyselect starts_with everything one_of
#' @importFrom stringr str_extract str_match str_detect str_remove_all
#' @importFrom tibble tibble
#' @importFrom tidyr pivot_longer fill
#' @importFrom gdata startsWith
#' @importFrom ggplot2 aes geom_line geom_path geom_point geom_polygon geom_text geom_vline ggtitle coord_flip theme_bw
#' theme element_text xlab ylab ggplot
#' @seealso \code{\link{get_gaps_fbhs}}
#' @export
get_plots_gap_fbh <- function (LAD_profiles,gap_cbh_metrics,min_height=1.5) {

  df_orig<-LAD_profiles


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


  df_orig$treeID<-factor(df_orig$treeID)
  treeID<-factor(df_orig$treeID)

  trees_name1a<- as.character(df_orig$treeID)
  trees_name3<- factor(unique(trees_name1a))

  df_metrics<- gap_cbh_metrics


  plot_list <- list()

  for (i in levels(trees_name3)){

    tree_data <- df_orig %>%
      dplyr::filter(treeID == i ) %>%
      dplyr::mutate(lad = as.numeric(lad)) %>%
      dplyr::filter(!is.na(lad))

    height<-df_orig$height
    lad<-df_orig$lad

    df_metrics11<-df_metrics %>% dplyr::filter(treeID == i )

    CBH_1<-as.numeric(as.character(df_metrics11$cbh1))
    CBH_2<-as.numeric(as.character(df_metrics11$cbh2))
    CBH_3<-as.numeric(as.character(df_metrics11$cbh3))
    CBH_4<-as.numeric(as.character(df_metrics11$cbh4))
    CBH_5<-as.numeric(as.character(df_metrics11$cbh5))
    CBH_6<-as.numeric(as.character(df_metrics11$cbh6))
    CBH_7<-as.numeric(as.character(df_metrics11$cbh7))
    CBH_8<-as.numeric(as.character(df_metrics11$cbh8))
    CBH_9<-as.numeric(as.character(df_metrics11$cbh9))

    GAP_1<-as.numeric(as.character(df_metrics11$gap1))
    GAP_2<-as.numeric(as.character(df_metrics11$gap2))
    GAP_3<-as.numeric(as.character(df_metrics11$gap3))
    GAP_4<-as.numeric(as.character(df_metrics11$gap4))
    GAP_5<-as.numeric(as.character(df_metrics11$gap5))
    GAP_6<-as.numeric(as.character(df_metrics11$gap6))
    GAP_7<-as.numeric(as.character(df_metrics11$gap7))
    GAP_8<-as.numeric(as.character(df_metrics11$gap8))
    GAP_9<-as.numeric(as.character(df_metrics11$gap9))



    bp <- ggplot(tree_data, aes(x = height)) +
      geom_line(aes(y = lad), color = "black", linewidth = 1.1) +
      geom_point(aes(y = lad), color = "black", size = 2.5) +

      geom_vline(xintercept = GAP_1 , linetype="dotted", color = "red", linewidth=1) +
      geom_vline(xintercept = GAP_2, linetype="dotted", color = "red", linewidth=1) +
      geom_vline(xintercept = GAP_3, linetype="dotted", color = "red", linewidth=1) +
      geom_vline(xintercept = GAP_4, linetype="dotted", color = "red", linewidth=1) +
      geom_vline(xintercept = GAP_5, linetype="dotted", color = "red", linewidth=1) +
      geom_vline(xintercept = GAP_6, linetype="dotted", color = "red", linewidth=1) +
      geom_vline(xintercept = GAP_7, linetype="dotted", color = "red", linewidth=1) +
      geom_vline(xintercept = GAP_8, linetype="dotted", color = "red", linewidth=1) +
      geom_vline(xintercept = GAP_9, linetype="dotted", color = "red", linewidth=1) +


      geom_vline(xintercept = CBH_1, color="dark green", linetype="twodash", linewidth=1) +
      geom_vline(xintercept = CBH_2, color="dark green", linetype="twodash", linewidth=1) +
      geom_vline(xintercept = CBH_3, color="dark green", linetype="twodash", linewidth=1) +
      geom_vline(xintercept = CBH_4, color="dark green", linetype="twodash", linewidth=1) +
      geom_vline(xintercept = CBH_5, color="dark green", linetype="twodash", linewidth=1) +
      geom_vline(xintercept = CBH_6, color="dark green", linetype="twodash", linewidth=1) +
      geom_vline(xintercept = CBH_7, color="dark green", linetype="twodash", linewidth=1) +
      geom_vline(xintercept = CBH_8, color="dark green", linetype="twodash", linewidth=1) +
      geom_vline(xintercept = CBH_9, color="dark green", linetype="twodash", linewidth=1) +

      ggtitle(paste0("tree_", i)) +

      theme_bw() +
      coord_flip() +
      theme(axis.title.x=element_text(size = 14, family = "sans", color = "black", face = "bold"),    # Adjust font size of x-axis label here
            axis.title.y=element_text(size = 14, family = "sans", color = "black", face = "bold"),    # Adjust font size of y-axis label here
            axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 1, color = "black", size = 14, family = "sans"),   # Size for x axis text
            axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 1, color = "black", size = 14, family = "sans"))

    plot_list[[i]] <- bp
    #print(paste("Plot for tree ", i, " created successfully"))
  }
  return(plot_list)
}



