# Build helper script for OptiSparseMET
remove.packages("OptiSparseMET")

unlink(file.path(.libPaths()[1], "00LOCK-OptiSparseMET"), recursive = TRUE)

# 2. Delete the compiled DLL from the source tree (src/*.so / src/*.dll).
#    Do this BEFORE the restart so there is nothing to unload conflicts with.
#devtools::clean_dll()

# 3. Restart R.  All loaded DLLs are released, file locks are cleared.
.rs.restartR()

# -- after restart, continue in SESSION B --------------------------------------
devtools::document()

source("data-raw/generate_example_data.R")
devtools::document()

source("scripts/diagnose_pev.R")


devtools::install()

devtools::test()

Sys.setenv(RUN_LIVE_API_TESTS = "true"); devtools::test(filter = "network-fetch"); Sys.unsetenv("RUN_LIVE_API_TESTS")


devtools::check()


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
