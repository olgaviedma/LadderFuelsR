
test_data <- read.table(system.file("extdata", "LAD_profiles.txt", package = "LadderFuelsR"), header = TRUE)
test_data$treeID <- factor(test_data$treeID)
trees_name1 <- as.character(test_data$treeID)
trees_name2 <- factor(unique(trees_name1))
test_data <- test_data |> dplyr::filter(treeID == "1_Eglin_zone1_CROWN")

test_data1 <- data.frame(
  treeID1 = 1,
  treeID = "1_CROWN",
  Hcbh1 = 5.5,
  dptf1 = 9,
  Hdptf1 = 14.5,
  effdist1 = 4,
  dist1 = 4,
  Hdist1 = 4.5,
  max_height = 14.5
)

output <- suppressWarnings(get_layers_lad(test_data, test_data1, threshold=10, step = 1, min_height= 1.5))


# Define the expected output for df1
expected_df1 <- data.frame(
  treeID = "1_CROWN",
  treeID1 = 1,
  Hcbh1 = 5.5,
  dptf1 = 9,
  Hdptf1 = 14.5,
  effdist1 = 4,
  dist1 = 4,
  Hdist1 = 4.5,
  max1 = 14.5,
  Hcbh1_Hdptf1 = 98.65299,
  max_height = 14.5,
  nlayers = 1
)

# Define the expected output for df2
expected_df2 <- data.frame(
  treeID = "1_CROWN",
  treeID1 = 1,
  Hcbh1 = 5.5,
  dptf1 = 9,
  Hdptf1 = 14.5,
  effdist1 = 4,
  Hdist1 = 4.5,
  max1 = 14.5,
  max_height = 14.5,
  Hcbh1_Hdptf1 = 98.65299,
  nlayers = 1
)

# Reset row names of actual data frame
rownames(output$df1) <- seq_len(nrow(expected_df1))
# Write the test for df1
test_that("get_layers_lad returns expected output with test data for df1", {
  expect_equal(output$df1, expected_df1, tolerance = 1e-6)
})

rownames(output$df2) <- seq_len(nrow(expected_df2))
test_that("get_layers_lad returns expected output with test data for df2", {
  expect_equal(output$df2, expected_df2, tolerance = 1e-6) # Adjust tolerance as needed
})
