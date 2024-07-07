
test_data <- data.frame(
  treeID1 = 1,
  treeID = "1_CROWN",
  Hdist1 = 4.5,
  Hcbh1 = 5.5,
  Hcbh1_Hdptf1 = 100,
  effdist1 = 4,
  dptf1 = 9,
  Hdptf1 = 14.5,
  max1 = 14.5,
  max_height = 14.5,
  nlayers = 1
)


output <- get_cbh_metrics(test_data,min_height= 1.5)

expected_output <- data.frame(
  treeID1 = 1,
  Hdist1 = 4.5,
  Hcbh1 = 5.5,
  effdist1 = 4,
  dptf1 = 9,
  Hdptf1 = 14.5,
  max1 = 14.5,
  Hcbh1_Hdptf1 = 100,
  treeID = "1_CROWN",
  max_height = 14.5,
  nlayers = 1,
  maxlad_Hcbh = 5.5,
  maxlad_Hdist = 4.5,
  maxlad_Hdptf = 14.5,
  maxlad_dptf = 9,
  maxlad_effdist = 4,
  maxlad_lad = 100,
  max_Hcbh = 5.5,
  max_Hdist = 4.5,
  max_Hdptf = 14.5,
  max_dptf = 9,
  max_effdist = 4,
  max_lad = 100,
  last_Hcbh = 5.5,
  last_Hdist = 4.5,
  last_Hdptf = 14.5,
  last_dptf = 9,
  last_effdist = 4,
  last_lad = 100
)

# Write the test
test_that("get_cbh_metrics returns expected output", {
  expect_equal(output, expected_output, tolerance = 1e-6)
})
