# Define the test data
test_data <- data.frame(
  treeID = "1_Eglin_zone1_CROWN",
  cbh1 = 1.5,
  gap1 = 2.5,
  gap2 = 5.5,
  cbh2 = 6.5,
  cbh3 = 9.5,
  cbh4 = 14.5,
  dist1 = 4,
  Hdist1 = 5.5,
  max_height = 15.5,
  treeID1 = 1
)

perc_data <- data.frame(
  height = c(1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 13.5, 14.5, 15.5),
  lad = c(0.0172, 0.0000, 0.0000, 0.0000, 0.0000, 0.1014, 0.2985, 0.0161, 0.1933, 0.2202, 0.2082, 0.0686, 0.0670, 0.0333, 0.0531),
  treeID = rep("1_Eglin_zone1_CROWN", 15),
  critical_points = c(-0.0000017595, 0.0015369406, 0.0029446957, 0.0032695284, 0.0014307565, -0.0037842900, -0.0090418045, -0.0064759674, -0.0096333229, -0.0108601690, -0.0089205902, -0.0038542723, -0.0011138678, 0.0003046800, -0.0000008436),
  percentil = c(50, 5, 5, 5, 5, 75, 100, 50, 90, 95, 90, 75, 75, 50, 50)
)

# Call the function with the test data
output <- get_distance(test_data, perc_data)

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
  max_height = 15.5
)

# Write the corrected test
test_that("get_distance returns expected output with test data", {
  expect_equal(as.data.frame(output), expected_output)
})

