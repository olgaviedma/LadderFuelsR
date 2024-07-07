
test_data <- data.frame(
  cbh1 = 5.5,
  depth1 = 9,
  dist1 = 4,
  Hdepth1 = 14.5,
  Hdist1 = 4.5,
  max_height = 14.5,
  treeID = "1_CROWN",
  treeID1 = 1
)

# Call the function with the test data
output <- get_real_fbh(test_data, step= 1, number_steps = 1, min_height=1.5)

expected_output <- data.frame(
  Hcbh1 = 5.5,
  depth1 = 9,
  dist1 = 4,
  Hdepth1 = 14.5,
  Hdist1 = 4.5,
  treeID1 = 1,
  treeID = "1_CROWN",
  max_height = 14.5
)


# Write the corrected test, ignoring row names during comparison
test_that("get_real_fbh returns expected output with test data", {
  expect_equal(output, expected_output, ignore_attr = T)
})


