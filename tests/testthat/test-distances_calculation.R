# Define the test data
test_data <- data.frame(
  treeID = "1_CROWN",
  gap1 = 1.5,
  gap2 = 4.5,
  cbh1 = 5.5,
  cbh2 = 7.5,
  cbh3 = 10.5,
  cbh4 = 12.5,
  cbh5 = 14.5,
  gap_lad1 = 0,
  gap_lad2 = 0,
  cbh_perc1 = 50,
  cbh_perc2 = 95,
  cbh_perc3 = 70,
  cbh_perc4 = 80,
  cbh_perc5 = 100,
  cbh_lad1 = 0,
  cbh_lad2 = 0,
  cbh_lad3 = 0,
  cbh_lad4 = 0,
  cbh_lad5 = 0,
  max_height = 14.5,
  treeID1 = 1
)

perc_data <- data.frame(
  height = c(1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 13.5, 14.5),
  lad = c(0.0000, 0.0000, 0.0000, 0.0000, 0.0294, 0.0558, 0.5521, 0.0185, 0.0172, 0.1221, 0.0985, 0.2506, 0.3198, 0.6785),
  treeID = rep("1_CROWN", 14),
  critical_points = c(-0.0000060210, 0.0057978945, 0.0138412870, 0.0184321963, 0.0092044330, -0.0194855834, -0.0766810373, 0.0018473451, 0.0385413648, 0.0442142395, 0.0549649985, 0.0461459480, 0.0366307724, -0.0000511847),
  percentil = c(5, 5, 5, 5, 50, 55, 95, 40, 35, 70, 65, 80, 85, 100)
)

# Call the function with the test data
output <- get_distance(test_data, perc_data, step = 1, min_height = 1.5)
cat("Unique treeIDs:", unique(test_data$treeID), "\n")

# Define expected output with correct data types and levels
expected_output <-  data.frame(
  treeID1 = factor(1),
  treeID = "1_CROWN",
  gap1 = 1.5,
  gap2 = 4.5,
  cbh1 = 5.5,
  cbh2 = 7.5,
  cbh3 = 10.5,
  cbh4 = 12.5,
  cbh5 = 14.5,
  dist1 = 4,
  max_height = 14.5,
  dist2 = 0,
  dist3 = 0,
  dist4 = 0,
  dist5 = 0,
  Hdist1 = 4.5,
  Hdist2 = 0,
  Hdist3 = 0,
  Hdist4 = 0,
  Hdist5 = 0
)

# Write the corrected test
test_that("get_distance returns expected output with test data", {
  # Convert treeID in output to character
  output$treeID <- as.character(output$treeID)

  # Ensure treeID1 has correct levels based on output$treeID1
  expected_output$treeID1 <- factor(expected_output$treeID1, levels = levels(output$treeID1))

  expect_equal(output, expected_output)
})

