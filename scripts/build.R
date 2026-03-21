# Build helper script for OptiSparseMET

source("data-raw/generate_example_data.R")
devtools::document()
devtools::test()
devtools::check()
devtools::install()

unloadNamespace("OptiSparseMET")
pkgdown::build_site()

devtools::build()


