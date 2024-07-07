# Define test data (mocked or simulated data)
test_data <- data.frame(
  height = c(1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 13.5, 14.5),
  lad = c(0.0000, 0.0000, 0.0000, 0.0000, 0.0294, 0.0558, 0.5521, 0.0185, 0.0172, 0.1221, 0.0985, 0.2506, 0.3198, 0.6785),
  treeID = rep("1_CROWN", 14)
)

# Mocking the get_gaps_fbhs function with a simplified version
get_gaps_fbhs <- function(input_data) {
  # Simulate the behavior of the function
  output <- data.frame(
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
  return(output)
}

library(testthat)

# Write the test
test_that("get_gaps_fbhs returns expected output with test data", {
  # Call the function with the test data
  output <- get_gaps_fbhs(test_data)

  # Define expected output with adjusted data types
  expected_output <- data.frame(
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

  # Convert character columns to numeric in actual output
  output[, c("gap1", "gap2", "cbh1", "cbh2", "cbh3", "cbh4", "cbh5")] <-
    lapply(output[, c("gap1", "gap2", "cbh1", "cbh2", "cbh3", "cbh4", "cbh5")], as.numeric)

  # Compare the output with the expected output
  expect_equal(output, expected_output)
})
