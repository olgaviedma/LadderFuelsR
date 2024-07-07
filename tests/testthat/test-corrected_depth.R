
test_data <- data.frame(
  treeID1 = factor("1", levels = "1"),
  treeID2 = factor("1_CROWN"),
  treeID = factor("1_CROWN"),
  Hdist1 = 4.5,
  Hcbh1 = 5.5,
  Hdepth1 = 14.5,
  dist1 = 4,
  depth1 = 9,
  max_height = 14.5
)


output <- suppressWarnings(get_real_depths(test_data, step=1, min_height=1.5))

expected_output <- data.frame(
  treeID = factor("1_CROWN"),
  treeID1 = factor("1", levels = "1"),
  Hdptf1 = 14.5,
  dist1 = 4,
  dptf1 = 9,
  Hcbh1 = 5.5,
  Hdist1 = 4.5,
  max_height = 14.5
)

# Write the test
test_that("get_real_depth returns expected output with test data", {
  expect_equal(output, expected_output)
})
