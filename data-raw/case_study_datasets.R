library(usethis)

case_study1 <- read.csv("data-raw/case_study1_dataset.csv", fileEncoding = "UTF-8-BOM")
case_study2 <- read.csv("data-raw/case_study2_dataset.csv", fileEncoding = "UTF-8-BOM")
case_study3 <- read.csv("data-raw/case_study3_dataset.csv", fileEncoding = "UTF-8-BOM")
coral_data <- read.csv("data-raw/extra_examples_coral_data.csv", fileEncoding = "UTF-8-BOM")


usethis::use_data(case_study1, overwrite = TRUE)
usethis::use_data(case_study2, overwrite = TRUE)
usethis::use_data(case_study3, overwrite = TRUE)
usethis::use_data(coral_data, overwrite = TRUE)