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
#' @param names Whether `variable` should contain output-column names
#'   (`"output"`, the backward-compatible default) or source API codes
#'   (`"api"`). The explicit `api_code` and `output_name` columns are always
#'   returned.
#'
#' @return A data frame with separate `api_code` and `output_name` columns,
#'   advisory `default_aggregation` metadata for weather, and property-specific
#'   depth, quantile, conversion, profile, and stock metadata for soil.
#'
#' @seealso [build_enviromic_covariates()], [fetch_soilgrids()],
#'   [build_environment_relationship()].
#' @examples
#' enviromic_variable_catalog("soil")
#' # pick limiting variables by their catalogued names:
#' # build_environment_relationship(X, source = "enviromic",
#' #   variables = c("mean_temp", "total_precip", "clay", "phh2o"))
#' @export
enviromic_variable_catalog <- function(source = c("all", "weather", "soil"),
                                       names = c("output", "api")) {
  source <- match.arg(source)
  name_view <- match.arg(names)

  api_code <- c("T2M", "T2M_MAX", "T2M_MIN", "PRECTOTCORR",
                "ALLSKY_SFC_SW_DWN", "RH2M", "T2MDEW", "T2MWET", "WS2M",
                "WS10M", "PS", "QV2M", "TS", "GWETTOP", "GWETROOT",
                "GWETPROF", "EVPTRNS", "EVLAND", "ALLSKY_SFC_PAR_TOT",
                "CLRSKY_SFC_SW_DWN")
  output_name <- vapply(api_code, .power_colname, character(1))
  weather <- data.frame(
    variable = if (name_view == "api") api_code else output_name,
    api_code = api_code,
    output_name = output_name,
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
    default_aggregation = c("mean", "mean", "mean", "sum", "mean", "mean",
                            "mean", "mean", "mean", "mean", "mean", "mean",
                            "mean", "mean", "mean", "mean", "sum", "sum",
                            "sum", "mean"),
    supported_depths = NA_character_, supported_quantiles = NA_character_,
    conversion_factor = NA_real_, default_depth = NA_character_,
    profile_variable = NA, stock_variable = NA,
    stringsAsFactors = FALSE)

  soil_meta <- .soilgrids_property_metadata()
  soil <- data.frame(
    variable = soil_meta$property,
    api_code = NA_character_,
    output_name = soil_meta$property,
    source = "soil",
    description = soil_meta$description,
    units = soil_meta$units,
    default_aggregation = NA_character_,
    supported_depths = vapply(soil_meta$supported_depths, paste,
                              collapse = ";", character(1)),
    supported_quantiles = vapply(soil_meta$supported_quantiles, paste,
                                 collapse = ";", character(1)),
    conversion_factor = soil_meta$conversion_factor,
    default_depth = soil_meta$default_depth,
    profile_variable = soil_meta$profile_variable,
    stock_variable = soil_meta$stock_variable,
    stringsAsFactors = FALSE)

  out <- switch(source,
                weather = weather,
                soil = soil,
                all = rbind(weather, soil))
  rownames(out) <- NULL
  out
}


#' List supported NASA POWER weather parameters
#'
#' Returns every supported POWER API code, named by the corresponding output
#' column. The result can be passed directly (after `unname()`) to
#' `weather_pars` or `pars` arguments.
#'
#' @return A named character vector of NASA POWER API codes.
#' @seealso [enviromic_variable_catalog()], [fetch_weather_series()].
#' @export
available_weather_parameters <- function() {
  catalog <- enviromic_variable_catalog("weather", names = "api")
  stats::setNames(catalog$api_code, catalog$output_name)
}
