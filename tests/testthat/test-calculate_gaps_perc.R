# Load test data
test_data <- read.table(system.file("extdata", "LAD_profiles.txt", package = "LadderFuelsR"), header = TRUE)
test_data$treeID <- factor(test_data$treeID)
trees_name1 <- as.character(test_data$treeID)
trees_name2 <- factor(unique(trees_name1))
test_data <- test_data |> dplyr::filter(treeID == "1_Eglin_zone1_CROWN")
test_data$expected_output <- c(40, 5, 5, 5, 5, 75, 100, 30, 80, 95, 90, 65, 60, 45, 50)

test_that("calculate_gaps_perc returns expected output with real data", {
  # Provide inputs using real test data
  input_data <- test_data

  # Assuming your function is defined in the LadderFuelsR package
  # Modify this part based on your actual function name and its parameters
  output <- calculate_gaps_perc(input_data,min_height=1.5)
  output<-output$percentil

  # Check the output against expected values
  expect_equal(output, test_data$expected_output)
})
