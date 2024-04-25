
test_data <- read.table(system.file("extdata", "LAD_profiles.txt", package = "LadderFuelsR"), header = TRUE)
test_data$treeID <- factor(test_data$treeID)
trees_name1 <- as.character(test_data$treeID)
trees_name2 <- factor(unique(trees_name1))
test_data <- test_data |> dplyr::filter(treeID == "1_Eglin_zone1_CROWN")

test_data1 <- data.frame(
  treeID1 = factor("1"),
  Hcbh1 = 6.5,
  Hdist1 = 5.5,
  dptf1 = 8,
  Hdptf1 = 14.5,
  effdist1 = 5,
  Hcbh1_Hdptf1 = 98.65299,
  treeID = factor("1_Eglin_zone1_CROWN"),
  max_height = 15.5,
  nlayers = 1,
  maxlad_Hcbh = 6.5,
  maxlad_Hdist = 5.5,
  maxlad_Hdptf = 14.5,
  maxlad_dptf = 8,
  maxlad_effdist = 5,
  maxlad_lad = 98.65299,
  max_Hcbh = 6.5,
  max_Hdist = 5.5,
  max_Hdptf = 14.5,
  max_dptf = 8,
  max_effdist = 5,
  max_lad = 98.65299,
  last_Hcbh = 6.5,
  last_Hdist = 5.5,
  last_Hdptf = 14.5,
  last_dptf = 8,
  last_effdist = 5,
  last_lad = 98.65299
)

set.seed(125)
output <- get_cum_break(test_data, test_data1, threshold=75)

expected_output <- data.frame(
  treeID = factor("1_Eglin_zone1_CROWN"),
  treeID1 = factor("1"),
  Hcbh_brpt = 6.43039,
  cumlad = 0.0172,
  below_hcbhbp = 1.347012,
  above_hcbhbp = 98.65299,
  bp_Hcbh = 6.4,
  bp_Hdptf = 15.5,
  bp_Hdist = 5.4,
  bp_effdist = 5,
  bp_dptf = 9,
  bp_lad = 98.65299,
  max_height = 15.5
)

rownames(output) <- seq_len(nrow(expected_output))

# Write the test
test_that("get_cbh_metrics returns expected output", {
  expect_equal(output, expected_output, tolerance = 1e-6)
})

