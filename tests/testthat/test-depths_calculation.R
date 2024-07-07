# Load test data
test_data <- read.table(system.file("extdata", "LAD_profiles.txt", package = "LadderFuelsR"), header = TRUE)
test_data$treeID <- factor(test_data$treeID)  # Ensure treeID is a factor

# Filter test data to specific treeID
test_data <- test_data %>%
  filter(treeID == "1_Eglin_zone1_CROWN")

# Define test data 1 (with expected input format)
test_data1 <- data.frame(
  cbh1 = 5.5,
  cbh2 = 7.5,
  cbh3 = 10.5,
  cbh4 = 12.5,
  cbh5 = 14.5,
  dist1 = 4,
  dist2 = 0,
  dist3 = 0,
  dist4 = 0,
  dist5 = 0,
  gap1 = 1.5,
  gap2 = 4.5,
  Hdist1 = 4.5,
  Hdist2 = 0,
  Hdist3 = 0,
  Hdist4 = 0,
  Hdist5 = 0,
  max_height = 14.5,
  treeID = "1_CROWN",
  treeID1 = 1
)

# Call the function with the test data
output <- get_depths(test_data, test_data1, step = 1, min_height = 1.5)
cat("Unique treeIDs:", unique(test_data$treeID), "\n")

# Define expected output with correct data types and levels
expected_output <- data.frame(
  treeID = "1_CROWN",
  treeID1 = factor(1, levels = levels(output$treeID1)),
  Hdepth1 = 14.5,
  dist1 = 4,
  Hdist1 = 4.5,
  depth1 = 9,
  cbh1 = 5.5,
  max_height = 14.5
)

# Convert treeID in output to character
output$treeID <- as.character(output$treeID)

# Write the corrected test
test_that("get_depth returns expected output with test data", {
  expect_equal(output, expected_output)
})
