
test_data <- data.frame(
  treeID = factor("1_Eglin_zone1_CROWN"),
  treeID1 = factor("1"),
  Hcbh1 = 6.5,
  Hdist1 = 5.5,
  dptf1 = 8,
  Hdptf1 = 14.5,
  Hcbh1_Hdptf1 = 98.65299, # Adjusted to match precision
  effdist1 = 5,
  max_height = 15.5,
  nlayers = 1
)

output <- get_cbh_metrics(test_data)

expected_output <- data.frame(
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

# Write the test
test_that("get_cbh_metrics returns expected output", {
  expect_equal(output, expected_output, tolerance = 1e-6)
})
