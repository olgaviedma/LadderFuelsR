
test_data <- read.table(system.file("extdata", "LAD_profiles.txt", package = "LadderFuelsR"), header = TRUE)
test_data$treeID <- factor(test_data$treeID)
trees_name1 <- as.character(test_data$treeID)
trees_name2 <- factor(unique(trees_name1))
test_data <- test_data |> dplyr::filter(treeID == "1_Eglin_zone1_CROWN")

test_data1 <- data.frame(
  treeID1 = factor("1", levels = "1"),
  treeID = factor("1_Eglin_zone1_CROWN"),
  cbh1 = 1.5,
  gap1 = 2.5,
  gap2 = 5.5,
  cbh2 = 6.5,
  cbh3 = 9.5,
  cbh4 = 14.5,
  dist1 = 4,
  Hdist1 = 5.5,
  max_height = 15.5
)

# Call the function with the test data
output <- get_depths(test_data, test_data1)

# Define expected output with correct data types and levels
expected_output <- data.frame(
  treeID1 = factor("1", levels = "1"),
  treeID = factor("1_Eglin_zone1_CROWN"),
  cbh1 = 1.5,
  gap1 = 2.5,
  gap2 = 5.5,
  cbh2 = 6.5,
  cbh3 = 9.5,
  cbh4 = 14.5,
  dist1 = 4,
  Hdist1 = 5.5,
  max_height = 15.5,
  Hdepth0 = 1.5,  # Adding Hdepth0, depth0, Hdepth1, and depth1 columns to expected_output
  depth0 = 1,
  Hdepth1 = 14.5,
  depth1 = 8
)

# Write the corrected test
test_that("get_depth returns expected output with test data", {
  expect_equal(as.data.frame(output), expected_output)
})
