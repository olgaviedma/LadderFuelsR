
test_data <- data.frame(
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

# Call the function with the test data
output <- get_real_fbh(test_data)

expected_output <- data.frame(
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


# Write the corrected test, ignoring row names during comparison
test_that("get_real_fbh returns expected output with test data", {
  expect_equal(output, expected_output, ignore_attr = T)
})


