# Build helper script for OptiSparseMET

source("data-raw/generate_example_data.R")
devtools::document()
devtools::test()
devtools::check()
devtools::install()

unloadNamespace("OptiSparseMET")
pkgdown::build_site()

pkgdown::build_site(override = list(template = list(favicon = FALSE)))

devtools::build()


# Create the correct folder
dir.create("pkgdown/favicon", showWarnings = FALSE)

# Copy all favicon files from assets/ to favicon/
file.copy(
  from      = list.files("pkgdown/assets", full.names = TRUE),
  to        = "pkgdown/favicon/",
  overwrite = TRUE
)

list.files("pkgdown/favicon")
