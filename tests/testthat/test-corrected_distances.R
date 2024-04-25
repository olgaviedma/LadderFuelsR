
test_data <- data.frame(
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

output <- get_effective_gap(test_data)

original_column_names <- colnames(output)

# Specify prefixes
desired_order <- c("treeID", "Hcbh", "dptf","effdist","dist", "Hdist", "Hdptf", "max_")

# Identify unique prefixes
prefixes <- unique(sub("^([a-zA-Z]+).*", "\\\\1", original_column_names))
# Initialize vector to store new order
new_order <- c()

# Loop over desired order of prefixes
for (prefix in desired_order) {
  # Find column names matching the current prefix
  matching_columns <- grep(paste0("^", prefix), original_column_names, value = TRUE)
  # Append to the new order
  new_order <- c(new_order, matching_columns)
}
output <- output[, new_order]


expected_output <- data.frame(
  treeID1 = factor("1", levels = "1"),
  treeID = factor("1_Eglin_zone1_CROWN"),
  Hcbh1 = 1.5,
  Hcbh2 = 6.5,
  dptf1 = 1,
  dptf2 = 8,
  effdist1 = 4,
  dist1 = 4,
  Hdist1 = 5.5,
  Hdptf1 = 1.5,
  Hdptf2 = 14.5,
  max_height = 15.5
)

# Write the test
test_that("get_effective_gap returns expected output with test data", {
  expect_equal(output, expected_output)
})
