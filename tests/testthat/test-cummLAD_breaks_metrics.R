
test_data <- read.table(system.file("extdata", "LAD_profiles.txt", package = "LadderFuelsR"), header = TRUE)
test_data$treeID <- factor(test_data$treeID)
trees_name1 <- as.character(test_data$treeID)
trees_name2 <- factor(unique(trees_name1))
test_data <- test_data |> dplyr::filter(treeID == "1_Eglin_zone1_CROWN")

test_data1 <- data.frame(
  treeID1 = 1,
  treeID = "1_CROWN",
  Hcbh1 = 5.5,
  Hcbh1_Hdptf1 = 100,
  dptf1 = 9,
  effdist1 = 4,
  Hdist1 = 4.5,
  Hdptf1 = 14.5,
  maxlad_Hcbh = 5.5,
  maxlad_Hdist = 4.5,
  maxlad_Hdptf = 14.5,
  maxlad_dptf = 9,
  maxlad_effdist = 4,
  maxlad_lad = 100,
  max_height = 14.5,
  max_Hcbh = 5.5,
  max_Hdist = 4.5,
  max_Hdptf = 14.5,
  max_dptf = 9,
  max_effdist = 4,
  max_lad = 100,
  last_Hcbh = 5.5,
  last_Hdist = 4.5,
  last_Hdptf = 14.5,
  last_dptf = 9,
  last_effdist = 4,
  last_lad = 100,
  nlayers = 1
)

set.seed(125)
output <- get_cum_break(test_data, test_data1, threshold=75, min_height= 1.5)

expected_output <- data.frame(
  treeID = "1_CROWN",
  treeID1 = 1,
  Hcbh_brpt = 6.43039,
  cumlad = 0.0172,
  below_hcbhbp = 1.347012,
  above_hcbhbp = 98.65299,
  bp_Hcbh = 6.4,
  bp_Hdptf = 14.5,
  bp_Hdist = 5.4,
  bp_effdist = 5,
  bp_dptf = 8,
  bp_lad = 98.65299,
  max_height = 14.5
)
rownames(output) <- seq_len(nrow(expected_output))

# Write the test
test_that("get_cbh_metrics returns expected output", {
  expect_equal(output, expected_output, tolerance = 1e-6)
})

