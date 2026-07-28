#' Catalogue of enviromic variable names, descriptions and units
#'
#' @description
#' The weather (NASA POWER) and soil (SoilGrids) fetchers return columns with
#' fixed but not always self-explanatory names (e.g. `phh2o`, `bdod`,
#' `total_precip`). A breeder choosing the *limiting variables* for
#' [build_environment_relationship()] needs to know which names are available
#' and what they mean. `enviromic_variable_catalog()` returns that reference
#' table, so the exact strings to pass to the `variables` argument can be read
#' off directly. The assembled matrix's own `colnames()` always lists what is
#' actually present for a given fetch; this catalogue explains each code.
#'
#' @param source Which variables to list: `"all"` (default), `"weather"`, or
#'   `"soil"`.
#'
#' @return A data frame with columns `variable` (the exact column name to use),
#'   `source`, `description`, and `units` (conventional units after the
#'   package's SoilGrids conversion factors are applied).
#'
#' @seealso [build_enviromic_covariates()], [fetch_soilgrids()],
#'   [build_environment_relationship()].
#' @examples
#' enviromic_variable_catalog("soil")
#' # pick limiting variables by their catalogued names:
#' # build_environment_relationship(X, source = "enviromic",
#' #   variables = c("mean_temp", "total_precip", "clay", "phh2o"))
#' @export
enviromic_variable_catalog <- function(source = c("all", "weather", "soil")) {
  source <- match.arg(source)

  # The first six are the default fetched columns (friendly names). The rest are
  # NASA POWER parameter codes you can request via `weather_pars` in
  # build_enviromic_covariates(); the fetched column keeps the POWER code.
  weather <- data.frame(
    variable = c("mean_temp", "max_temp", "min_temp", "total_precip",
                 "mean_radiation", "mean_humidity",
                 "T2MDEW", "T2MWET", "WS2M", "WS10M", "PS", "QV2M", "TS",
                 "GWETTOP", "GWETROOT", "GWETPROF", "EVPTRNS", "EVLAND",
                 "ALLSKY_SFC_PAR_TOT", "CLRSKY_SFC_SW_DWN"),
    source = "weather",
    description = c("Mean 2 m air temperature over the window (POWER T2M)",
                   "Mean daily maximum 2 m air temperature (T2M_MAX)",
                   "Mean daily minimum 2 m air temperature (T2M_MIN)",
                   "Total precipitation over the window (PRECTOTCORR)",
                   "Mean all-sky downward shortwave radiation (ALLSKY_SFC_SW_DWN)",
                   "Mean 2 m relative humidity (RH2M)",
                   "2 m dew/frost point temperature",
                   "2 m wet-bulb temperature",
                   "2 m wind speed",
                   "10 m wind speed",
                   "Surface pressure",
                   "2 m specific humidity",
                   "Earth skin temperature",
                   "Top-layer (0-5 cm) soil wetness (fraction)",
                   "Root-zone soil wetness (fraction)",
                   "Profile soil moisture (fraction)",
                   "Evapotranspiration",
                   "Land evaporation",
                   "All-sky photosynthetically active radiation",
                   "Clear-sky downward shortwave radiation"),
    units = c("degC", "degC", "degC", "mm", "W/m^2", "%",
              "degC", "degC", "m/s", "m/s", "kPa", "kg/kg", "degC",
              "fraction", "fraction", "fraction", "mm/day", "mm/day",
              "W/m^2", "W/m^2"),
    stringsAsFactors = FALSE)

  soil <- data.frame(
    variable = c("bdod", "cec", "cfvo", "clay", "nitrogen", "ocd", "ocs",
                 "phh2o", "sand", "silt", "soc",
                 "wv0010", "wv0033", "wv1500"),
    source = "soil",
    description = c("Bulk density of the fine earth fraction",
                   "Cation exchange capacity (at pH 7)",
                   "Volumetric fraction of coarse fragments (> 2 mm)",
                   "Clay (< 0.002 mm) mass fraction",
                   "Total nitrogen",
                   "Organic carbon density",
                   "Organic carbon stock",
                   "Soil pH in water",
                   "Sand (0.05-2 mm) mass fraction",
                   "Silt (0.002-0.05 mm) mass fraction",
                   "Soil organic carbon content",
                   "Volumetric water content at 10 kPa",
                   "Volumetric water content at 33 kPa",
                   "Volumetric water content at 1500 kPa"),
    units = c("kg/dm^3", "cmol(c)/kg", "vol%", "%", "g/kg", "kg/m^3",
              "kg/m^2", "-", "%", "%", "g/kg", "vol%", "vol%", "vol%"),
    stringsAsFactors = FALSE)

  out <- switch(source,
                weather = weather,
                soil = soil,
                all = rbind(weather, soil))
  rownames(out) <- NULL
  out
}
