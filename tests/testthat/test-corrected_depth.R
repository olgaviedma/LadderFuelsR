
test_data <- data.frame(
  treeID = factor("1_Eglin_zone1_CROWN"),
  treeID1 = factor("1"),
  Hdepth1 = 1.5,
  depth1 = 1,
  dist1 = 4,
  Hdist1 = 5.5,
  Hdepth2 = 14.5,
  depth2 = 8,
  max1 = 15.5,
  Hcbh1 = 1.5,
  Hcbh2 = 6.5,
  Hcbh3 = 6.5,
  Hcbh4 = 6.5,
  max_height=15.5
)


output <- suppressWarnings(get_real_depths(test_data))

expected_output <- data.frame(
  treeID1 = factor("1", levels = "1"),
  treeID = factor("1_Eglin_zone1_CROWN"),
  Hcbh1 = 1.5,
  Hcbh2 = 6.5,
  dptf1 = 1,
  dptf2 = 8,
  Hdptf1 = 1.5,
  Hdptf2 = 14.5,
  dist1 = 4,
  Hdist1 = 5.5,
  max_height = 15.5
)

# Write the test
test_that("get_real_depth returns expected output with test data", {
  expect_equal(output, expected_output)
})
