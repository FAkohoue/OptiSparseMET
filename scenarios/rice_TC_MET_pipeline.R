###############################################################################
## OptiSparseMET -- end-to-end pipeline for a real planning scenario
##
## 300 testcross (TC) hybrids = 150 restorers (12 families) x 2 females,
## all genotyped with 800 SNPs. Evaluate across 7 sites:
##   - 4 sites in Ecuador     : NO field-capacity limit, 10 g seed/plot, direct
##   - Saldana (Colombia)     : 200 plots (incl. checks), 16 g/plot, direct
##   - Palmira (CIAT)         : 200 plots (incl. checks), 20 g/plot, transplant
##   - Peru                   : 182 plots (incl. checks), 16 g/plot, transplant
## Seed availability per TC: 30-300 g. Only management info = planting mode.
## Design is p-rep or alpha-lattice ONLY (never augmented), subject to constraints.
##
## We are PLANNING: the trials have not run, so weather is historical DAILY data
## and Sigma_E is calibrated from valid historical MET values when available.
## Without such values Sigma_E is identity; enviromic kernels remain descriptive
## and are evaluated separately as sensitivity structures. Run section by section;
## anything that errors, paste back.
###############################################################################

library(OptiSparseMET)
set.seed(2026)

## This script never guesses whether a production file is authoritative. Set
## `use_simulated_inputs <- FALSE` only after placing all four input files below.
## `offline_validation_mode` avoids network calls and keeps optimisation short;
## it is intended for an installation smoke test, not for scientific inference.
use_simulated_inputs    <- TRUE
use_dominance_kernel    <- TRUE   # appropriate for these hybrids; FALSE = additive only
offline_validation_mode <- identical(Sys.getenv("OPTISPARSEMET_OFFLINE"), "true")
run_design_optimisation <- !offline_validation_mode
outdir <- "scenarios/outputs"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

required_input_files <- list(
  parental_markers = "scenarios/inputs/rice_parental_markers.csv",
  crosses          = "scenarios/inputs/rice_testcrosses.csv",
  seed_inventory   = "scenarios/inputs/rice_seed_inventory.csv",
  sites            = "scenarios/inputs/rice_sites.csv"
)
optional_input_files <- list(
  management       = "scenarios/inputs/rice_management.csv",
  phenology        = "scenarios/inputs/rice_phenology.csv",
  station_weather  = "scenarios/inputs/rice_station_weather.csv",
  historical_met   = "scenarios/inputs/rice_historical_met_values.csv",
  weather_features = "scenarios/inputs/rice_weather_features.csv",
  soil_features    = "scenarios/inputs/rice_soil_features.csv",
  environmental_interactions =
    "scenarios/inputs/rice_environmental_interactions.csv"
)
## Optional schemas:
## management: environment + planned management variables (one row/site).
## phenology: environment, stage, and start_date/end_date (or start_day/end_day).
## station_weather: environment, year, date (preferred) or day, weather columns.
## historical_met: genotype, environment, value, and optionally year; `value`
## must be an adjusted genetic value/BLUE, not an unadjusted plot observation.
## weather_features / soil_features: environment + audited precomputed
## covariates; these take precedence over network retrieval.
## environmental_interactions: interaction, parent, modality, variable.
## This optional internal file defines a short, pre-specified set of dedicated
## variable interactions. It is not an exhaustive interaction search.

## Replace these priors with posterior means from historical METs. They are kept
## explicit so additive, dominance, and residual scales cannot be conflated.
variance_prior <- list(sigma2_A = 2, sigma2_D = 1, sigma2_e = 1)
simulation_replicates <- if (offline_validation_mode) 8L else 100L

## ---- PLANNING CONFIG (edit these) ------------------------------------------
## How many PAST seasons of DAILY weather to characterise each site with, the
## intended growing window (month-day), which NASA POWER variables to pull, and
## the fetch retry count. All exposed here -- nothing hard-coded downstream.
weather_years  <- 2015:2024                       # <- number of past years
## SITE-SPECIFIC growing windows (month-day). Peru's season differs from
## Colombia/Ecuador and crosses the calendar year (Nov -> Feb).
weather_window <- data.frame(
  environment = c("Ecu_Babahoyo", "Ecu_Daule", "Ecu_Samborondon", "Ecu_Quevedo",
                  "Saldana", "Palmira", "Peru_Chiclayo"),
  start_md = c(rep("06-01", 6), "11-01"),
  end_md   = c(rep("09-30", 6), "02-28"),          # Peru end in the NEXT year
  stringsAsFactors = FALSE)
weather_pars   <- c("T2M", "T2M_MAX", "T2M_MIN",              # temperature
                    "PRECTOTCORR", "RH2M", "T2MDEW",         # rain, humidity, dew
                    "ALLSKY_SFC_SW_DWN", "WS2M",             # radiation, wind
                    "GWETTOP", "GWETROOT", "GWETPROF", "EVPTRNS")  # soil water + ET
weather_max_tries <- 4L
weather_cache  <- file.path(outdir, "power_cache")  # cache -> fast re-runs
weather_workers <- 1L                               # >1 needs future.apply (parallel)
weather_min_coverage <- 0.80
weather_bandwidth_multipliers <- c(0.5, 1, 2)
soil_depths <- c("0-5cm", "5-15cm", "15-30cm", "30-60cm",
                 "60-100cm")
soil_quantiles <- c("Q0.05", "Q0.5", "Q0.95")
soil_properties <- c("bdod", "cec", "cfvo", "clay", "nitrogen", "phh2o",
                     "sand", "silt", "soc", "wv0033", "wv1500")
soil_root_depth_cm <- 100
soil_spatial_summary <- "mean_sd"

## DESIGN CONFIG (edit these)
## n_common: give a number to fix it, or leave NULL to optimize its size,
## identities, and treatment-by-site replication jointly. "Common" requires
## presence at every site; it does NOT require equal replication across sites.
n_common <- NULL
## For unrestricted Ecuador sites, evaluate a nested interval and stop at the
## first diminishing-return knee. Change the interval and threshold in section 6
## during sensitivity analysis.
ecuador_cap_interval <- c(70, 220)
ecuador_cap_step     <- 25

## ===========================================================================
## 1. GERMPLASM DATA
## ===========================================================================
if (use_simulated_inputs) {
  n_restorers <- 150L; n_fam <- 12L; n_snp <- 800L
  restorers <- sprintf("R%03d", seq_len(n_restorers))
  females <- c("FA", "FB")
  fam_of_restorer <- setNames(
    paste0("Fam", sprintf("%02d",
                           sample(seq_len(n_fam), n_restorers, replace = TRUE))),
    restorers)
  maf <- runif(n_snp, 0.05, 0.5)
  parents <- c(restorers, females)
  Mpar <- vapply(maf, function(p) 2L * rbinom(length(parents), 1, p),
                 numeric(length(parents)))
  rownames(Mpar) <- parents
  colnames(Mpar) <- sprintf("SNP%04d", seq_len(n_snp))

  tc_grid <- expand.grid(Restorer = restorers, Female = females,
                         stringsAsFactors = FALSE)
  tc_ids <- sprintf("TC_%s_%s", tc_grid$Restorer, tc_grid$Female)
  tc_grid$Treatment <- tc_ids
  tc_grid$Family <- unname(fam_of_restorer[tc_grid$Restorer])
  seed_info <- data.frame(
    Treatment = tc_ids,
    SeedAvailable = round(runif(length(tc_ids), 30, 300)),
    stringsAsFactors = FALSE)
} else {
  missing_inputs <- names(required_input_files)[
    !file.exists(unlist(required_input_files))]
  if (length(missing_inputs))
    stop("Missing production input file(s): ",
         paste(missing_inputs, collapse = ", "), ".")

  marker_df <- read.csv(required_input_files$parental_markers, check.names = FALSE,
                        stringsAsFactors = FALSE)
  tc_grid <- read.csv(required_input_files$crosses, stringsAsFactors = FALSE)
  seed_info <- read.csv(required_input_files$seed_inventory, stringsAsFactors = FALSE)
  required_cross <- c("Treatment", "Restorer", "Female", "Family")
  if (!"Parent" %in% names(marker_df) ||
      !all(required_cross %in% names(tc_grid)) ||
      !all(c("Treatment", "SeedAvailable") %in% names(seed_info)))
    stop("Production schemas must contain Parent + SNP columns; ",
         "Treatment/Restorer/Female/Family; and Treatment/SeedAvailable.")
  Mpar <- data.matrix(marker_df[, setdiff(names(marker_df), "Parent"),
                               drop = FALSE])
  rownames(Mpar) <- as.character(marker_df$Parent)
  tc_ids <- as.character(tc_grid$Treatment)
}

if (anyDuplicated(rownames(Mpar)) || anyDuplicated(tc_ids) ||
    any(!is.finite(Mpar)) || any(!Mpar %in% 0:2))
  stop("Parental marker dosages must be finite 0/1/2 with unique parent IDs.")
if (!all(c(tc_grid$Restorer, tc_grid$Female) %in% rownames(Mpar)))
  stop("Every Restorer and Female must occur in the parental marker matrix.")
if (anyDuplicated(seed_info$Treatment) ||
    !all(tc_ids %in% seed_info$Treatment) ||
    any(!is.finite(seed_info$SeedAvailable)) ||
    any(seed_info$SeedAvailable < 0))
  stop("Seed inventory must uniquely cover every testcross with non-negative grams.")
seed_info <- seed_info[match(tc_ids, seed_info$Treatment), , drop = FALSE]
treatment_info <- data.frame(
  Treatment = tc_ids, Family = as.character(tc_grid$Family),
  stringsAsFactors = FALSE)

## Hybrid dosages are used only for the optional hybrid dominance kernel.
Mtc <- 0.5 * (Mpar[tc_grid$Restorer, , drop = FALSE] +
                Mpar[tc_grid$Female, , drop = FALSE])
rownames(Mtc) <- tc_ids
write.csv(treatment_info, file.path(outdir, "treatment_info.csv"), row.names = FALSE)
write.csv(seed_info,      file.path(outdir, "seed_info.csv"),      row.names = FALSE)
cat("Germplasm:", length(tc_ids), "testcrosses;",
    length(unique(tc_grid$Restorer)), "restorers;",
    length(unique(tc_grid$Family)), "families;",
    length(unique(tc_grid$Female)), "females\n")

## ===========================================================================
## 2. SITES (GPS, capacity, seed/plot, planting mode)
## ===========================================================================
if (use_simulated_inputs) {
  sites <- data.frame(
    environment = c("Ecu_Babahoyo", "Ecu_Daule", "Ecu_Samborondon", "Ecu_Quevedo",
                    "Saldana", "Palmira", "Peru_Chiclayo"),
    country = c("Ecuador","Ecuador","Ecuador","Ecuador",
                "Colombia","Colombia","Peru"),
    latitude = c(-1.80, -1.86, -1.96, -1.03,  3.93,  3.50, -6.77),
    longitude = c(-79.53,-79.98,-79.72,-79.46,-75.02,-76.35,-79.84),
    total_plots = c(NA, NA, NA, NA, 200, 200, 182),
    seed_per_plot = c(10,10,10,10, 16, 20, 16),
    planting_mode = c("direct","direct","direct","direct",
                      "direct","transplanting","transplanting"),
    stringsAsFactors = FALSE)
} else {
  sites <- read.csv(required_input_files$sites, stringsAsFactors = FALSE)
}
required_site <- c("environment", "country", "latitude", "longitude",
                   "total_plots", "seed_per_plot", "planting_mode")
if (!all(required_site %in% names(sites)) ||
    anyDuplicated(sites$environment) ||
    any(!is.finite(sites$latitude)) || any(!is.finite(sites$longitude)) ||
    any(!is.finite(sites$seed_per_plot)) || any(sites$seed_per_plot <= 0))
  stop("Site data must contain unique valid environments, coordinates, ",
       "plot limits, positive seed rates, and planting modes.")
rownames(sites) <- sites$environment
if (!all(c("Saldana", "Palmira", "Peru_Chiclayo") %in% sites$environment) ||
    sum(is.na(sites$total_plots)) != 4L)
  stop("This rice scenario requires the three named constrained partners and ",
       "exactly four unconstrained sites.")
if (!setequal(weather_window$environment, sites$environment))
  stop("`weather_window` must cover exactly the environments in `sites`.")
write.csv(sites, file.path(outdir, "sites.csv"), row.names = FALSE)
print(sites[, c("environment","country","total_plots","seed_per_plot","planting_mode")])

## ===========================================================================
## 3. GENOMIC RELATIONSHIP among the testcrosses
##    Additive relationship is derived from the parental GRM on the parental
##    allele-frequency scale. Dominance is an explicit, optional hybrid/SCA
##    component; ordinary inbred-line trials should leave it disabled.
## ===========================================================================
G_parent_A <- build_relationship_matrix(
  markers = Mpar, type = "genomic", method = "AGHmatrix",
  relationship = "additive", ploidy = 2, tuneup = TRUE)
G_A <- build_hybrid_relationship(
  G_parent_A, tc_grid[, c("Restorer", "Female")], hybrid_ids = tc_ids)
G_A <- tune_relationship_matrix(G_A, method = "auto")

if (use_dominance_kernel) {
  if (any(!Mtc %in% 0:2))
    stop("Dominance requires integer 0/1/2 hybrid dosages; disable ",
         "`use_dominance_kernel` or supply valid inbred-parent crosses.")
  G_D <- build_relationship_matrix(
    markers = Mtc, type = "genomic", method = "AGHmatrix",
    relationship = "dominance", ploidy = 2, tuneup = TRUE)
  G <- combine_relationship_matrices(
    list(additive = G_A, dominance = G_D),
    weights = c(variance_prior$sigma2_A, variance_prior$sigma2_D))
  sigma_g2 <- variance_prior$sigma2_A + variance_prior$sigma2_D
  relationship_label <- "additive + optional hybrid dominance"
} else {
  G <- G_A
  sigma_g2 <- variance_prior$sigma2_A
  relationship_label <- "additive"
}
G <- tune_relationship_matrix(G, method = "auto")
cat("Relationship:", relationship_label, "|", nrow(G), "hybrids\n")
## population-structure PCs (ASRgenomics::kinship.pca, eigen fallback)
pcs <- kinship_pca(G, ncp = 13)

## ===========================================================================
## 4. ENVIRONMENTAL COVARIATES + RELATIONSHIP (planning: historical DAILY weather)
## ===========================================================================
mgmt <- if (!use_simulated_inputs &&
            file.exists(optional_input_files$management)) {
  read.csv(optional_input_files$management, check.names = FALSE,
           stringsAsFactors = FALSE)
} else {
  data.frame(environment = sites$environment,
             planting_mode = sites$planting_mode,
             stringsAsFactors = FALSE)
}
if (!"environment" %in% names(mgmt) || anyDuplicated(mgmt$environment) ||
    !setequal(as.character(mgmt$environment), sites$environment))
  stop("Management data must uniquely cover all environments. Prefer site-year ",
       "records of irrigation, N/P/K dose and timing, establishment, density, ",
       "water regime, and crop protection when these are available.")
mgmt <- mgmt[match(sites$environment, mgmt$environment), , drop = FALSE]
site_ll <- sites[, c("environment","latitude","longitude")]
station_daily <- if (!use_simulated_inputs &&
                     file.exists(optional_input_files$station_weather))
  read.csv(optional_input_files$station_weather, check.names = FALSE,
           stringsAsFactors = FALSE) else NULL
phenology <- if (!use_simulated_inputs &&
                 file.exists(optional_input_files$phenology))
  read.csv(optional_input_files$phenology, stringsAsFactors = FALSE) else NULL

## Planning-stage weather MUST be DAILY, to keep within-season day-to-day
## variation and its timing. NASA POWER "climatology" is only MONTHLY/annual,
## so instead we characterise each site from several PAST seasons of DAILY data
## (historical_envirotype -> fetch_weather_series at daily resolution). This
## returns the TYPICAL daily-derived profile (across years), tails/trends, and
## interannual VARIABILITY. Each source is fetched independently so a soil
## outage cannot discard valid weather, management, or geographic information.
## NOTE on soil: SoilGrids is a single STATIC model (no annual layers), so soil
## texture/pH/bulk-density are fetched once. DYNAMIC soil-water across years is
## captured by the daily NASA POWER GWET* / EVPTRNS variables above. Supply your
## own multi-year soil via `soil = <data.frame>` if you have field measurements.
env_hist <- NULL
weather_block <- NULL
weather_error <- NULL
if (!use_simulated_inputs &&
    file.exists(optional_input_files$weather_features)) {
  wf <- read.csv(optional_input_files$weather_features, check.names = FALSE,
                 stringsAsFactors = FALSE)
  if (!"environment" %in% names(wf) || anyDuplicated(wf$environment))
    stop("Precomputed weather features need unique environment rows.")
  rownames(wf) <- wf$environment
  weather_block <- data.matrix(wf[, setdiff(names(wf), "environment"),
                                  drop = FALSE])
} else if (!offline_validation_mode) {
  env_hist <- tryCatch(
    historical_envirotype(
    sites = site_ll, years = weather_years, window = weather_window,
    pars = weather_pars, max_tries = weather_max_tries,
    cache_dir = weather_cache, workers = weather_workers,
    station_daily = station_daily,
    station_correction = if (is.null(station_daily)) "none" else "mean_bias",
    feature_function = "rice",
    min_daily_coverage = weather_min_coverage,
    daily_missing_action = "impute",
    min_years = max(3L, ceiling(length(weather_years) * 0.60)),
    summary_components = c("mean", "sd", "q10", "q50", "q90", "trend"),
    envirotype = list(
      windows = list(vegetative = c(0, 40), reproductive = c(40, 80),
                     ripening = c(80, 120)),
      phenology = phenology,
      time_scale = "calendar",
      stats = c("mean", "sd", "min", "max", "sum"))),
    error = function(e) {
      weather_error <<- conditionMessage(e)
      message("Historical weather unavailable: ", weather_error)
      NULL
    })
  if (!is.null(env_hist)) weather_block <- env_hist$combined
}

soil_raw <- NULL
soil_block <- NULL
soil_error <- NULL
if (!use_simulated_inputs && file.exists(optional_input_files$soil_features)) {
  sf <- read.csv(optional_input_files$soil_features, check.names = FALSE,
                 stringsAsFactors = FALSE)
  if (!"environment" %in% names(sf) || anyDuplicated(sf$environment))
    stop("Precomputed soil features need unique environment rows.")
  rownames(sf) <- sf$environment
  soil_block <- data.matrix(sf[, setdiff(names(sf), "environment"),
                               drop = FALSE])
} else if (!offline_validation_mode) {
  soil_raw <- tryCatch(
    fetch_soilgrids(
      site_ll, backend = "wcs", properties = soil_properties,
      depth = soil_depths, quantile = soil_quantiles,
      buffer_m = 1000, spatial_summary = soil_spatial_summary
    ),
    error = function(e) {
      soil_error <<- conditionMessage(e)
      message("Soil profile unavailable: ", soil_error)
      NULL
    })
  if (!is.null(soil_raw)) {
    soil_block <- tryCatch(
      soil_profile_features(soil_raw, root_depth_cm = soil_root_depth_cm),
      error = function(e) {
        soil_error <<- conditionMessage(e)
        message("Soil profile summarisation unavailable: ", soil_error)
        NULL
      })
  }
}

## Historical MET values are optional. They calibrate central covariance
## weights, but their absence must not suppress environmental interaction
## hypotheses.
historical_met <- if (!use_simulated_inputs &&
                      file.exists(optional_input_files$historical_met))
  read.csv(optional_input_files$historical_met, stringsAsFactors = FALSE) else NULL
dedicated_interactions <- if (
  !use_simulated_inputs &&
  file.exists(optional_input_files$environmental_interactions)
) {
  read.csv(
    optional_input_files$environmental_interactions,
    stringsAsFactors = FALSE, check.names = FALSE
  )
} else NULL

kernel_build <- build_environment_kernels(
  weather = weather_block, soil = soil_block, management = mgmt,
  geography = site_ll,
  environments = sites$environment,
  bandwidth_multipliers = weather_bandwidth_multipliers,
  min_coverage = weather_min_coverage,
  missing_action = "impute", impute = "median",
  include_interactions = FALSE
)
environment_kernels <- kernel_build$kernels

stab <- NULL
if (!is.null(env_hist) && "weather" %in% names(environment_kernels)) {
  cat("Weather characterised from", env_hist$n_years,
      "past seasons of DAILY data (stages, tails, trends and variability)\n")
  cat("Interannual variability (mean across features) per site:\n")
  print(round(rowMeans(env_hist$variability, na.rm = TRUE), 3))
  stab <- tryCatch(
    assess_envirotype_stability(
      env_hist, kernel = "gaussian",
      bandwidth_multipliers = weather_bandwidth_multipliers,
      n_boot = 500L, seed = 2026
    ),
    error = function(e) {
      message("Weather stability assessment unavailable: ",
              conditionMessage(e))
      NULL
    })
  if (!is.null(stab)) {
    Ktyp <- environment_kernels$weather
    Kyear <- stab$consensus_D[sites$environment, sites$environment,
                             drop = FALSE]
    environment_kernels$weather <- consensus_relationship(
      list(historical_summary = Ktyp, yearly_consensus = Kyear),
      method = "rv_weighted", use_off_diagonal = TRUE
    )
    cat("Weather stability:", round(stab$mean_stability, 3),
        "bootstrap", paste(round(stab$mean_stability_ci, 3), collapse = " to "),
        "\n")
  }
}

## Always construct interaction structures. Historical adjusted genetic
## values/BLUEs, when supplied, determine only whether an interaction receives
## weight in the central genetic covariance. Without them, Sigma_E is identity
## and every main or interaction kernel remains an unweighted structural
## uncertainty scenario. Environmental similarity is never converted into an
## assumed genetic covariance by arbitrary prior weights.
main_names <- intersect(
  c("weather", "soil", "management", "geography"),
  names(environment_kernels)
)
environment_kernels <- add_environment_kernel_interactions(
  environment_kernels[main_names], mode = "anova"
)
variable_kernel_build <- NULL
if (!is.null(dedicated_interactions)) {
  variable_kernel_build <- build_variable_interaction_kernels(
    covariates = kernel_build$covariates,
    interactions = dedicated_interactions,
    hierarchy = "strong",
    orthogonalise = TRUE,
    max_tensor_columns = 100L
  )
  environment_kernels <- c(
    environment_kernels, variable_kernel_build$kernels
  )
}
present <- names(environment_kernels)
primary_kernel_names <- intersect(
  c("weather", "soil", "management", "geography"), present
)
primary_environment_kernels <- environment_kernels[primary_kernel_names]
variable_interaction_assessment <- NULL
eligible_interactions <- NULL
if (!is.null(variable_kernel_build)) {
  variable_interaction_assessment <- assess_variable_interactions(
    environment_kernels,
    historical = historical_met,
    year_col = if (!is.null(historical_met) &&
                   "year" %in% names(historical_met)) "year" else NULL,
    p_adjust = "holm",
    n_boot = if (is.null(historical_met)) 0L else 500L,
    seed = 2026
  )
  all_interaction_names <- names(environment_kernels)[vapply(
    names(environment_kernels),
    function(nm) {
      identical(
        attr(environment_kernels[[nm]], "kernel_role"), "interaction"
      ) || grepl("_x_", nm, fixed = TRUE)
    },
    logical(1)
  )]
  broad_interaction_names <- setdiff(
    all_interaction_names,
    names(variable_kernel_build$hierarchy)
  )
  eligible_interactions <- unique(c(
    broad_interaction_names,
    variable_interaction_assessment$supported
  ))
}
environment_covariance <- calibrate_environment_covariance(
  environment_kernels, historical = historical_met,
  year_col = if (!is.null(historical_met) &&
                 "year" %in% names(historical_met)) "year" else NULL,
  ## With historical MET values, ridge = 0 makes the reported weights entirely
  ## data-calibrated rather than shrunk toward user-defined modality weights.
  ridge = 0, eligible_interactions = eligible_interactions,
  interaction_policy = "evidence",
  n_boot = if (is.null(historical_met)) 0L else 500L,
  seed = 2026
)
Sigma_E <- environment_covariance$Sigma_E
if (!is.null(stab) && length(stab$per_year_D)) {
  weather_year_candidates <- list()
  weather_interaction_names <- intersect(
    c("weather_x_soil", "weather_x_management"),
    names(environment_kernels)
  )
  weather_interaction_year_candidates <- stats::setNames(
    vector("list", length(weather_interaction_names)),
    weather_interaction_names
  )
  for (yy in names(stab$per_year_D)) {
    Ky <- environment_kernels[intersect(
      c("weather", "soil", "management", "geography"),
      names(environment_kernels)
    )]
    Ky$weather <- stab$per_year_D[[yy]][sites$environment,
                                       sites$environment, drop = FALSE]
    ## Year-specific weather may change both its main-effect structure and its
    ## joint structure with soil and management, irrespective of whether the
    ## central interaction model passed the historical evidence gate.
    Ky <- add_environment_kernel_interactions(Ky, mode = "anova")
    for (interaction_name in weather_interaction_names) {
      weather_interaction_year_candidates[[interaction_name]][[yy]] <-
        Ky[[interaction_name]]
    }
    if (identical(environment_covariance$status, "historically_calibrated")) {
      feature_weights <- environment_covariance$weights[names(Ky)]
      identity_weight_realised <- environment_covariance$weights["identity"]
      if (sum(feature_weights) > .Machine$double.eps &&
          identity_weight_realised < 1 - .Machine$double.eps) {
        weather_year_candidates[[paste0("weather_year_", yy)]] <-
          combine_environment_kernels(
            Ky, weights = feature_weights / sum(feature_weights),
            identity_weight = identity_weight_realised
          )
      } else {
        weather_year_candidates[[paste0("weather_year_", yy)]] <- Ky$weather
      }
    } else {
      ## Without MET evidence, retain a weather-year kernel as a separate
      ## sensitivity structure; do not mix it with other modalities.
      weather_year_candidates[[paste0("weather_year_", yy)]] <- Ky$weather
    }
  }
  if (length(weather_year_candidates)) {
    dist_from_central <- vapply(
      weather_year_candidates,
      function(S) sqrt(sum((S - Sigma_E)^2)), numeric(1)
    )
    environment_covariance$candidates <- c(
      environment_covariance$candidates, weather_year_candidates
    )
    environment_covariance$candidates$weather_anomalous_year <-
      weather_year_candidates[[which.max(dist_from_central)]]
  }
  ## Retain the most atypical observed-year structure for every interaction
  ## involving weather. This captures across-year changes in joint similarity
  ## even when that interaction has zero central weight.
  for (interaction_name in names(weather_interaction_year_candidates)) {
    year_set <- weather_interaction_year_candidates[[interaction_name]]
    if (!length(year_set)) next
    central_interaction <- environment_kernels[[interaction_name]]
    distance_from_interaction <- vapply(
      year_set,
      function(S) sqrt(sum((S - central_interaction)^2)),
      numeric(1)
    )
    environment_covariance$candidates[[
      paste0(interaction_name, "_anomalous_year")
    ]] <- year_set[[which.max(distance_from_interaction)]]
  }
}
candidate_priority <- c(
  "calibrated", "independent", "bootstrap_q10", "bootstrap_q90",
  paste0("kernel_", c("weather", "soil", "management", "geography",
                      "weather_x_soil", "weather_x_management",
                      "soil_x_management")),
  "weather_anomalous_year", "weather_x_soil_anomalous_year",
  "weather_x_management_anomalous_year"
)
candidate_keep <- intersect(candidate_priority,
                            names(environment_covariance$candidates))
Sigma_E_sensitivity <- environment_covariance$candidates[candidate_keep]
if (!length(Sigma_E_sensitivity))
  Sigma_E_sensitivity <- list(central = Sigma_E)
design_uncertainty_scenarios <- robust_scenarios(
  sigma_g2 = sigma_g2,
  sigma_e2 = variance_prior$sigma2_e * c(0.5, 1, 2),
  sigmaE_shrink = c(0, 0.25, 0.5),
  Sigma_E_candidates = Sigma_E_sensitivity
)
capacity_uncertainty_scenarios <- robust_scenarios(
  sigma_g2 = sigma_g2,
  sigma_e2 = variance_prior$sigma2_e,
  sigmaE_shrink = 0,
  Sigma_E_candidates = Sigma_E_sensitivity
)
## Structural alternatives do not have defensible probabilities. Maximin makes
## them operational without manufacturing weights from the presence or absence
## of historical MET responses.
covariance_aggregate <- "min"

## D is a descriptive multimodal environmental relationship for clustering and
## representativeness. Sigma_E is the separately labelled calibrated covariance,
## or identity without historical MET evidence, used by MET calculations.
if (identical(environment_covariance$status, "historically_calibrated")) {
  d_weights <- environment_covariance$weights[present]
  D <- if (sum(d_weights) > 0) {
    combine_environment_kernels(
      environment_kernels, weights = d_weights / sum(d_weights),
      identity_weight = 0
    )
  } else {
    consensus_environment_kernels(primary_environment_kernels)
  }
} else {
  ## Use each primary modality exactly once. Interaction kernels remain separate
  ## sensitivity structures so duplicated weather/soil information cannot act as
  ## an implicit weight in descriptive clustering.
  D <- consensus_environment_kernels(primary_environment_kernels)
}
dimnames(D) <- list(sites$environment, sites$environment)
dimnames(Sigma_E) <- dimnames(D)
cat("Environmental covariance status:", environment_covariance$status, "\n")
if (is.null(environment_covariance$weights)) {
  cat("Environmental kernel weights: not estimated (no historical MET).\n")
} else {
  print(round(environment_covariance$weights, 3))
}
cat("D / Sigma_E:", nrow(D), "x", ncol(D), "\n"); print(round(D, 2))

write.csv(kernel_build$audit,
          file.path(outdir, "environment_data_audit.csv"), row.names = FALSE)
write.csv(kernel_build$provenance,
          file.path(outdir, "environment_data_provenance.csv"), row.names = FALSE)
write.csv(kernel_build$imputation,
          file.path(outdir, "environment_data_imputation.csv"), row.names = FALSE)
write.csv(kernel_build$block_diagnostics,
          file.path(outdir, "environment_block_diagnostics.csv"),
          row.names = FALSE)
write.csv(kernel_build$kernel_agreement,
          file.path(outdir, "environment_kernel_agreement.csv"),
          row.names = FALSE)
if (!is.null(variable_kernel_build)) {
  write.csv(
    variable_kernel_build$audit,
    file.path(outdir, "environment_variable_interaction_audit.csv"),
    row.names = FALSE
  )
  write.csv(
    variable_kernel_build$variable_ledger,
    file.path(outdir, "environment_variable_ledger.csv"),
    row.names = FALSE
  )
}
if (!is.null(variable_interaction_assessment)) {
  write.csv(
    variable_interaction_assessment$evidence,
    file.path(outdir, "environment_variable_interaction_evidence.csv"),
    row.names = FALSE
  )
}
bandwidth_ledger <- do.call(rbind, lapply(names(kernel_build$bandwidths),
                                         function(nm) data.frame(
  kernel = nm,
  bandwidth = as.numeric(kernel_build$bandwidths[[nm]]),
  stringsAsFactors = FALSE
)))
write.csv(bandwidth_ledger,
          file.path(outdir, "environment_kernel_bandwidths.csv"),
          row.names = FALSE)
weather_provenance <- if (!is.null(env_hist)) {
  wp <- env_hist$weather_provenance
  wp <- wp[!vapply(wp, is.null, logical(1))]
  if (length(wp)) do.call(rbind, wp) else NULL
} else NULL
if (!is.null(weather_provenance))
  write.csv(weather_provenance,
            file.path(outdir, "weather_source_provenance.csv"),
            row.names = FALSE)
if (!is.null(env_hist) && !is.null(env_hist$daily_audit))
  write.csv(env_hist$daily_audit,
            file.path(outdir, "weather_daily_audit.csv"), row.names = FALSE)
if (!is.null(env_hist) && !is.null(env_hist$daily_imputation))
  write.csv(env_hist$daily_imputation,
            file.path(outdir, "weather_daily_imputation.csv"),
            row.names = FALSE)
soil_provenance <- attr(soil_raw, "provenance")
if (!is.null(soil_provenance))
  write.csv(soil_provenance,
            file.path(outdir, "soil_source_provenance.csv"),
            row.names = FALSE)
soil_feature_audit <- attr(soil_block, "soil_audit")
if (!is.null(soil_feature_audit))
  write.csv(soil_feature_audit,
            file.path(outdir, "soil_feature_audit.csv"), row.names = FALSE)
source_manifest <- data.frame(
  source = c("weather", "soil", "management", "geography"),
  available = c(!is.null(weather_block), !is.null(soil_block), TRUE, TRUE),
  input_file = c(
    if (file.exists(optional_input_files$weather_features))
      optional_input_files$weather_features else "",
    if (file.exists(optional_input_files$soil_features))
      optional_input_files$soil_features else "",
    if (file.exists(optional_input_files$management))
      optional_input_files$management else "",
    if (!use_simulated_inputs) required_input_files$sites else ""
  ),
  error = c(if (is.null(weather_error)) "" else weather_error,
            if (is.null(soil_error)) "" else soil_error, "", ""),
  stringsAsFactors = FALSE)
source_manifest$md5 <- vapply(source_manifest$input_file, function(f)
  if (nzchar(f) && file.exists(f)) unname(tools::md5sum(f)) else "",
  character(1))
write.csv(source_manifest,
  file.path(outdir, "environment_source_manifest.csv"), row.names = FALSE)
write.csv(environment_covariance$diagnostics,
          file.path(outdir, "environment_covariance_diagnostics.csv"),
          row.names = FALSE)
write.csv(environment_covariance$interaction_evidence,
          file.path(outdir, "environment_interaction_evidence.csv"),
          row.names = FALSE)
cat("Environmental interaction status:",
    environment_covariance$interaction_status, "\n")
weight_ledger <- if (is.null(environment_covariance$weights)) {
  data.frame(
    kernel = c(names(environment_kernels), "identity"),
    weight = NA_real_,
    status = "not_estimated_without_historical_MET",
    stringsAsFactors = FALSE
  )
} else {
  data.frame(
    kernel = names(environment_covariance$weights),
    weight = as.numeric(environment_covariance$weights),
    status = "historically_calibrated",
    stringsAsFactors = FALSE
  )
}
write.csv(weight_ledger,
          file.path(outdir, "environment_kernel_weights.csv"), row.names = FALSE)
if (!is.null(environment_covariance$cross_validation))
  write.csv(environment_covariance$cross_validation,
            file.path(outdir, "environment_covariance_cross_validation.csv"),
            row.names = FALSE)
if (!is.null(stab))
  write.csv(data.frame(
    mean_stability = stab$mean_stability,
    ci_lower = stab$mean_stability_ci["lower"],
    ci_upper = stab$mean_stability_ci["upper"]),
    file.path(outdir, "weather_kernel_stability.csv"), row.names = FALSE)
for (nm in names(environment_kernels))
  write.csv(data.frame(environment = rownames(environment_kernels[[nm]]),
                       environment_kernels[[nm]], check.names = FALSE),
            file.path(outdir, paste0("kernel_", nm, ".csv")),
            row.names = FALSE)
for (nm in names(environment_covariance$candidates))
  write.csv(data.frame(
    environment = rownames(environment_covariance$candidates[[nm]]),
    environment_covariance$candidates[[nm]], check.names = FALSE),
    file.path(outdir, paste0("Sigma_E_candidate_", nm, ".csv")),
    row.names = FALSE)

## Mega-environments are inferred rather than fixed by the analyst. Each
## modality is one evidence block; repeated weather years remain within the
## weather block, so a long weather archive cannot outvote soil, management,
## geography, or historical MET evidence. The selector compares hierarchical
## clustering with kernel-PCA k-means over every defensible k, then requires
## minimum cluster size, silhouette separation, and adjusted-Rand stability.
## If no candidate is stable, all sites remain explicitly unpartitioned.
grouping_relationships <- list()
grouping_blocks <- character(0)
for (nm in names(primary_environment_kernels)) {
  key <- paste0("modality_", nm)
  grouping_relationships[[key]] <- primary_environment_kernels[[nm]]
  grouping_blocks[key] <- nm
}
if (!is.null(stab) && length(stab$per_year_D)) {
  for (yy in names(stab$per_year_D)) {
    key <- paste0("weather_year_", yy)
    grouping_relationships[[key]] <-
      stab$per_year_D[[yy]][sites$environment, sites$environment,
                            drop = FALSE]
    grouping_blocks[key] <- "weather"
  }
}
if (identical(environment_covariance$status, "historically_calibrated") &&
    is.matrix(environment_covariance$target) &&
    all(is.finite(environment_covariance$target))) {
  grouping_relationships$historical_MET <- environment_covariance$target
  grouping_blocks["historical_MET"] <- "historical_MET"
}
mega <- infer_mega_environments(
  D,
  relationships = grouping_relationships,
  relationship_groups = grouping_blocks,
  Sigma_E = Sigma_E,
  n_boot = 500L,
  seed = 2026L
)
cat("Mega-environment inference:", mega$status, "-", mega$reason, "\n")
print(mega$membership)
write.csv(mega$diagnostics,
          file.path(outdir, "mega_environment_diagnostics.csv"),
          row.names = FALSE)
genetic_mega_supported <- mega$hard_groups &&
  identical(environment_covariance$status, "historically_calibrated") &&
  "historical_MET" %in% names(grouping_relationships)
write.csv(data.frame(
  Status = mega$status,
  StableEnvironmentalGrouping = mega$hard_groups,
  GeneticMegaEnvironmentSupported = genetic_mega_supported,
  InferredGroups = mega$n_clusters,
  SelectedMethod = if (is.na(mega$method)) "" else mega$method,
  RelationshipMatrices = mega$resampling$n_relationships,
  IndependentEvidenceBlocks = mega$resampling$n_relationship_blocks,
  BootstrapDraws = mega$resampling$n_boot_used,
  Reason = mega$reason,
  stringsAsFactors = FALSE
), file.path(outdir, "mega_environment_status.csv"), row.names = FALSE)
mega_label <- if (identical(mega$status, "unstable")) {
  rep("UNPARTITIONED", nrow(D))
} else if (identical(mega$status, "provisional")) {
  paste0("PROVISIONAL_ME", unname(mega$membership[rownames(D)]))
} else if (!genetic_mega_supported) {
  paste0("ENVIRONMENTAL_STRATUM_", unname(mega$membership[rownames(D)]))
} else {
  paste0("ME", unname(mega$membership[rownames(D)]))
}

## All seven contracted sites remain in the MET. This table quantifies which
## sites best represent the full TPE, without dropping capacity-constrained
## partners merely because another site is highly correlated with them.
environment_priority <- data.frame(
  Environment = rownames(D),
  TPERepresentativeness = rowMeans(D),
  MegaEnvironment = mega_label,
  MegaEnvironmentStatus = mega$status,
  StableEnvironmentalGrouping = mega$hard_groups,
  GeneticMegaEnvironmentSupported = genetic_mega_supported,
  stringsAsFactors = FALSE)
environment_priority <- environment_priority[
  order(environment_priority$TPERepresentativeness, decreasing = TRUE), ]
write.csv(environment_priority,
          file.path(outdir, "environment_representativeness.csv"),
          row.names = FALSE)

## ===========================================================================
## 5. FIELD-CAPACITY RECONCILIATION (plots -> test entries)
## ===========================================================================
n_checks <- 3L; n_blocks <- 4L
prep_check_plots <- n_checks * n_blocks
alpha_reps <- 2L
alpha_check_plots <- alpha_reps * n_checks * n_blocks
prep_frac <- 0.3                                  # 30% of entries replicated (p-rep)

## Designs are p-rep or alpha-lattice ONLY (never augmented). Capacity-limited
## sites use PARTIALLY REPLICATED designs, so distinct entries fit the budget as
## (budget - checks) / (1 + prep_frac); `prep` entries are then replicated twice.
lim_sites <- c("Saldana", "Palmira", "Peru_Chiclayo")
lim_info <- lapply(lim_sites, function(s) {
  B <- sites[s, "total_plots"]
  E <- test_entry_capacity(B, n_checks = n_checks, n_blocks = n_blocks,
                           avg_reps_per_entry = 1 + prep_frac)
  list(entries = E, prep = (B - prep_check_plots) - E, budget = B)
})
names(lim_info) <- lim_sites
cap_limited <- vapply(lim_info, function(x) as.integer(x$entries), integer(1))

ecuador_sites <- sites$environment[is.na(sites$total_plots)]
capacity_for_sweep <- stats::setNames(integer(nrow(sites)), sites$environment)
capacity_for_sweep[ecuador_sites] <- ecuador_cap_interval[1]
capacity_for_sweep[lim_sites] <- cap_limited
plot_multiplier <- stats::setNames(
  ifelse(sites$environment %in% ecuador_sites, alpha_reps, 1 + prep_frac),
  sites$environment)
check_overhead <- stats::setNames(
  ifelse(sites$environment %in% ecuador_sites,
         alpha_check_plots, prep_check_plots),
  sites$environment)

## Only Ecuador varies in this sweep. Partner commitments remain fixed, and the
## same nested entry order and Monte Carlo stream are used at every candidate.
sc <- suggest_site_capacity(
  G, Sigma_E,
  candidate_plots = seq(ecuador_cap_interval[1], ecuador_cap_interval[2],
                        by = ecuador_cap_step),
  scope = "subset", focal_envs = ecuador_sites,
  site_capacities = capacity_for_sweep,
  plots_per_entry = plot_multiplier,
  check_plots_per_site = check_overhead,
  sigma_g2 = sigma_g2, sigma_e2 = variance_prior$sigma2_e,
  n_sim = simulation_replicates, select = "diminishing",
  min_gain = 0.0005,
  ## Capacity is maximin across structural covariance candidates at the nominal
  ## variance ratio. The fuller variance/shrinkage grid is used below for the
  ## common set, allocation refinement, and final evaluation.
  robust = capacity_uncertainty_scenarios,
  robust_aggregate = covariance_aggregate,
  seed = 1)
cap_ecuador <- min(sc$recommended_plots, length(tc_ids))
write.csv(sc$table, file.path(outdir, "ecuador_capacity_sweep.csv"),
          row.names = FALSE)
cat("Ecuador entries (nested alpha-capacity sweep):", cap_ecuador,
    "| limited-site entries (p-rep):\n"); print(cap_limited)

capacity <- capacity_for_sweep
capacity[ecuador_sites] <- cap_ecuador
print(capacity)

## ===========================================================================
## 6. ROBUST GLOBAL COMMON SET (BINARY PRESENCE)
## ===========================================================================
seed_required_per_plot <- data.frame(
  Environment = sites$environment, SeedRequiredPerPlot = sites$seed_per_plot,
  stringsAsFactors = FALSE)
minimum_seed_reserve <- 0
cat("One inventory will fund the entire network; reserve:",
    minimum_seed_reserve, "g per testcross\n")

## A common TC occurs at all seven sites and counts once for every environment
## pair, irrespective of its later local plot count. The common-set optimizer
## chooses size and identities only. Ecuador alpha designs subsequently require
## two plots for every allocated entry; partner p-rep sites start with one plot
## and assign additional plots from the remaining network seed. Variance and
## Sigma_E uncertainty use maximin aggregation. Even with historical MET
## calibration, unsupported or weakly
## estimable interactions remain physically plausible structural scenarios and
## therefore receive neither arbitrary probabilities nor silent exclusion.
## Genetic effective sample size remains a separate diversity criterion.
mandatory_local_plots <- stats::setNames(
  ifelse(sites$environment %in% ecuador_sites, alpha_reps, 1L),
  sites$environment)
common_seed_cost <- stats::setNames(
  sites$seed_per_plot * mandatory_local_plots,
  sites$environment)
common_scenarios <- design_uncertainty_scenarios

common_opt <- optimize_common_treatments(
  G = G, Sigma_E = Sigma_E,
  treatment_info = treatment_info, seed_info = seed_info,
  ## Effective mandatory seed debit for one presence at each environment.
  seed_required_per_plot = common_seed_cost,
  entry_capacities = capacity,
  n_common = n_common,
  minimum_seed_buffer = minimum_seed_reserve,
  scenarios = common_scenarios, aggregate = covariance_aggregate,
  cvar_alpha = 0.25,
  target_se = 0.15, min_per_family = 1,
  objective_weights = c(reliability = 0.40, connectivity = 0.30,
                        genetic_diversity = 0.15,
                        testing_breadth = 0.15),
  pair_aggregate = "maximin", seed = 11)
common_treatments <- common_opt$selected
n_common <- common_opt$n_common
if (any(capacity < n_common))
  stop("At least one site cannot hold the optimized common treatments.")
write.csv(common_opt$comparison,
          file.path(outdir, "common_count_optimization.csv"),
          row.names = FALSE)
write.csv(common_opt$common_presence,
          file.path(outdir, "common_presence_matrix.csv"))
write.csv(common_opt$seed_ledger,
          file.path(outdir, "common_seed_ledger.csv"), row.names = FALSE)
write.csv(common_opt$pairwise_connectivity,
          file.path(outdir, "common_pairwise_connectivity.csv"),
          row.names = FALSE)
write.csv(common_opt$selection_diagnostics,
          file.path(outdir, "common_selection_diagnostics.csv"),
          row.names = FALSE)
cat("Optimized common TCs:", n_common,
    if (isTRUE(common_opt$rationale$selected_count_was_user_fixed))
      "(count fixed by user)" else "(count optimized)", "\n")
cat("Common presence is binary and fixed at every environment:\n")
print(colSums(common_opt$common_presence))

## ===========================================================================
## 7. SPARSE ALLOCATION + WITHIN-ENVIRONMENT DESIGNS
## ===========================================================================
## Design is p-rep or alpha-lattice ONLY (never augmented). Capacity-limited
## sites -> PARTIALLY REPLICATED (p-rep, `prep` entries replicated twice, fits the
## fixed budget). Unlimited Ecuador sites -> ALPHA-LATTICE (2 reps). The field
## grid is the smallest exact rectangle for the planned number of occupied plots.
exact_grid <- function(n) {
  divisors <- which(n %% seq_len(floor(sqrt(n))) == 0L)
  nr <- max(divisors)
  c(n_rows = nr, n_cols = n %/% nr)
}

env_specs <- lapply(sites$environment, function(s) {
  base <- list(check_treatments = paste0("CHK", 1:n_checks),
               check_families = rep("CHECK", n_checks),
               cluster_source = "Family", use_dispersion = FALSE,
               order = "row", serpentine = TRUE)
  if (is.na(sites[s, "total_plots"])) {
    used <- alpha_reps * (capacity[s] + n_blocks * n_checks)
    grid <- exact_grid(used)
    c(base, list(design = "met_alpha_rc_stream", n_reps = alpha_reps,
                 n_blocks_per_rep = n_blocks,
                 n_rows = unname(grid["n_rows"]),
                 n_cols = unname(grid["n_cols"]), verbose = FALSE))
  } else {
    grid <- exact_grid(sites[s, "total_plots"])
    c(base, list(design = "met_prep_famoptg", replication_mode = "p_rep",
                 desired_replications = 2L, max_prep = lim_info[[s]]$prep,
                 max_extra_replication_plots = lim_info[[s]]$prep,
                 shortage_action = "downgrade", n_blocks = n_blocks,
                 n_rows = unname(grid["n_rows"]),
                 n_cols = unname(grid["n_cols"])))
  }
})
names(env_specs) <- sites$environment

plan <- plan_sparse_met_design(
  treatments = tc_ids, environments = sites$environment,
  allocation_method = "random_balanced",
  n_test_entries_per_environment = capacity,
  common_treatments = common_treatments,
  balance = "env_pair", balance_iter = 3000L, balance_seed = 71,
  Sigma_E = Sigma_E_sensitivity, pair_aggregate = "maximin",
  env_design_specs = env_specs,
  treatment_info = treatment_info, seed_info = seed_info,
  seed_required_per_plot = seed_required_per_plot,
  minimum_seed_buffer = minimum_seed_reserve,
  allocation_group_source = "Family", seed = 7)
M <- plan$sparse_allocation$allocation_matrix
cat("Allocation matrix:", nrow(M), "testcrosses x", ncol(M),
    "sites | site sizes:\n")
print(colSums(M))
write.csv(M, file.path(outdir, "allocation_matrix.csv"))

## ---- FIELD BOOKS + FIELD MAPS per site (and combined) ----------------------
draw_field_map <- function(fb, title) {
  if (is.null(fb) || !all(c("Row", "Column") %in% names(fb))) {
    plot.new(); title(title); return(invisible())
  }
  nr <- max(fb$Row, na.rm = TRUE); nc <- max(fb$Column, na.rm = TRUE)
  chk_col <- intersect(c("Check", "IsCheck", "is_check"), names(fb))
  is_check <- if (length(chk_col)) as.logical(fb[[chk_col[1]]]) else
    grepl("CHK|CHECK", as.character(fb$Treatment), ignore.case = TRUE)
  m <- matrix(NA_real_, nr, nc)
  for (i in seq_len(nrow(fb)))
    if (is.finite(fb$Row[i]) && is.finite(fb$Column[i]))
      m[fb$Row[i], fb$Column[i]] <- if (isTRUE(is_check[i])) 2 else 1
  image(seq_len(nr), seq_len(nc), m, col = c("#3AA0FF", "#D9534F"),
        axes = FALSE, xlab = "", ylab = "", main = title, zlim = c(1, 2))
  graphics::box()
}

if (!is.null(plan)) {
  ed <- plan$environment_designs
  ## individual field books per site -> CSV
  for (s in names(ed)) {
    fb <- ed[[s]]$field_book
    if (!is.null(fb))
      write.csv(fb, file.path(outdir, paste0("fieldbook_", s, ".csv")),
                row.names = FALSE)
  }
  ## combined MET field book -> CSV
  if (!is.null(plan$combined_field_book))
    write.csv(plan$combined_field_book,
              file.path(outdir, "fieldbook_combined.csv"), row.names = FALSE)
  cat("Field books written:", length(ed), "per-site CSVs + fieldbook_combined.csv\n")
  ## a peek at the first site's field book
  cat("Field book head (", names(ed)[1], "):\n", sep = "")
  print(utils::head(ed[[1]]$field_book, 6))
  ## FIELD MAPS (spatial layout, blue = entry, red = check) -> PDF
  pdf(file.path(outdir, "field_maps.pdf"), width = 11, height = 6)
  op <- par(mfrow = c(2, 4), mar = c(2, 2, 2.5, 1))
  for (s in names(ed)) draw_field_map(ed[[s]]$field_book, s)
  par(op); dev.off()
  cat("Field maps -> field_maps.pdf\n")
}

## The statistical model must use delivered fieldbooks, not intended average
## replication. Count each candidate's actual plots in each environment.
reps_matrix <- matrix(
  0, nrow = length(tc_ids), ncol = nrow(sites),
  dimnames = list(tc_ids, sites$environment))
physical_plots <- stats::setNames(integer(nrow(sites)), sites$environment)
for (s in sites$environment) {
  fb <- plan$environment_designs[[s]]$field_book
  trt <- as.character(fb$Treatment)
  physical_plots[s] <- sum(!is.na(trt) & nzchar(trt))
  tt <- table(trt[trt %in% tc_ids])
  reps_matrix[names(tt), s] <- as.numeric(tt)
}
M_actual <- 1L * (reps_matrix > 0)
if (!identical(unname(M_actual), unname(M)))
  stop("Delivered fieldbooks do not reproduce the planned incidence matrix.")

expected_physical <- stats::setNames(
  ifelse(sites$environment %in% ecuador_sites,
         alpha_reps * (capacity + n_blocks * n_checks),
         sites$total_plots),
  sites$environment)
if (any(physical_plots != expected_physical))
  stop("Seed-feasible local replication did not fill the planned field size: ",
       paste(names(physical_plots)[physical_plots != expected_physical],
             collapse = ", "), ". Reduce capacity/replication or increase seed.")
if (is.null(plan$seed_ledger) || any(!plan$seed_ledger$Feasible) ||
    min(plan$seed_ledger$SeedRemaining) < minimum_seed_reserve - 1e-8)
  stop("Final fieldbooks violate the network-wide seed reserve.")
write.csv(plan$seed_ledger, file.path(outdir, "seed_ledger.csv"),
          row.names = FALSE)
cat("Physical plots per site:\n"); print(physical_plots)
cat("Total candidate seed allocated:", sum(plan$seed_ledger$SeedAllocated),
    "g | minimum reserve:", min(plan$seed_ledger$SeedRemaining), "g\n")

## ===========================================================================
## 8. MET-LEVEL INFORMATION, SIMULATION (accuracy / gain / reliability)
## ===========================================================================
efficiency <- stats::setNames(rep(1, nrow(sites)), sites$environment)
if ("eff_A" %in% names(plan$environment_summary)) {
  ee <- plan$environment_summary$eff_A
  names(ee) <- plan$environment_summary$Environment
  ok <- is.finite(ee) & ee > 0 & ee <= 1
  efficiency[names(ee)[ok]] <- ee[ok]
}
info <- met_information(
  M_actual, G = G, Sigma_E = Sigma_E, reps = reps_matrix,
  env_efficiency = efficiency, sigma_g2 = sigma_g2,
  sigma_e2 = variance_prior$sigma2_e)
sim <- simulate_met(
  M_actual, G = G, Sigma_E = Sigma_E, reps = reps_matrix,
  env_efficiency = efficiency, sigma_g2 = sigma_g2,
  sigma_e2 = variance_prior$sigma2_e,
  n_sim = simulation_replicates, select_fraction = 0.1, seed = 3)
cat(sprintf(paste0("MET: mean_PEV=%.3f CDmean=%.3f | accuracy=%.3f ",
                   "(95%% CI %.3f, %.3f) gain=%.3f (95%% CI %.3f, %.3f)\n"),
            info$mean_PEV, info$CDmean, sim$accuracy_mean,
            sim$accuracy_ci95[1], sim$accuracy_ci95[2], sim$gain_mean,
            sim$gain_ci95[1], sim$gain_ci95[2]))

## A mega-environment target is evaluated only when a stable grouping also has
## historical genetic-response support. Purely environmental or provisional
## strata remain diagnostics and never alter treatment allocation, common-set
## optimisation, or breeding-value targets.
if (genetic_mega_supported) {
  mega_simulation <- do.call(rbind, lapply(
    seq_along(mega$clusters),
    function(i) {
      members <- mega$clusters[[i]]
      z <- simulate_met(
        M_actual, G = G, Sigma_E = Sigma_E, reps = reps_matrix,
        env_efficiency = efficiency, sigma_g2 = sigma_g2,
        sigma_e2 = variance_prior$sigma2_e,
        n_sim = simulation_replicates, select_fraction = 0.1,
        bv_target = "mega_environment", target_envs = members,
        seed = 3000L + i
      )
      data.frame(
        MegaEnvironment = names(mega$clusters)[i],
        NEnvironments = length(members),
        Environments = paste(members, collapse = ";"),
        AccuracyMean = z$accuracy_mean,
        AccuracySE = z$accuracy_se,
        AccuracyCI_Lower = z$accuracy_ci95[1],
        AccuracyCI_Upper = z$accuracy_ci95[2],
        GainMean = z$gain_mean,
        GainSE = z$gain_se,
        GainCI_Lower = z$gain_ci95[1],
        GainCI_Upper = z$gain_ci95[2],
        stringsAsFactors = FALSE
      )
    }
  ))
  write.csv(mega_simulation,
            file.path(outdir, "mega_environment_simulation.csv"),
            row.names = FALSE)
}

## ===========================================================================
## 9. CONSTRAINED OPTIMISATION + ROBUST PRIOR EVALUATION
## ===========================================================================
mandatory_seed_cost <- stats::setNames(
  sites$seed_per_plot *
    ifelse(sites$environment %in% ecuador_sites, alpha_reps, 1),
  sites$environment)
mandatory_seed_used <- rowSums(sweep(M, 2L, mandatory_seed_cost, `*`))
extra_replication_seed <- pmax(
  0, plan$seed_ledger$SeedAllocated[
    match(rownames(M), plan$seed_ledger$Treatment)] - mandatory_seed_used)
optimizer_seed_info <- data.frame(
  Treatment = rownames(M),
  SeedAvailable = plan$seed_ledger$SeedAvailable[
    match(rownames(M), plan$seed_ledger$Treatment)] - extra_replication_seed,
  stringsAsFactors = FALSE)

opt <- NULL
if (run_design_optimisation) {
  opt <- optimize_design(
    M, G = G, Sigma_E = Sigma_E,
    objective = list(weights = list(gain = 1, reliability = 0, cost = 0)),
    sigma_g2 = sigma_g2, sigma_e2 = variance_prior$sigma2_e,
    robust = design_uncertainty_scenarios,
    robust_aggregate = covariance_aggregate,
    preserve = "margins",                 # exact site sizes and TC replication
    seed_available = optimizer_seed_info,
    seed_required_per_environment = mandatory_seed_cost,
    minimum_seed_buffer = minimum_seed_reserve,
    environment_capacities = capacity,
    n_starts = 4, iters = 250, seed = 5)
  cat(sprintf("optimize_design: start=%.4f optimised=%.4f\n",
              opt$score_start, opt$score))
}

## Evaluate the deliverable fieldbooks, including their actual replication,
## across plausible residual variances and uncertainty in Sigma_E.
scen <- design_uncertainty_scenarios
rob <- robust_design_score(
  M_actual, G = G, Sigma_E = Sigma_E, scenarios = scen,
  aggregate = covariance_aggregate,
  reps = reps_matrix, env_efficiency = efficiency,
  weights = list(gain = 1))
cat(sprintf("robust %s score across %d variance scenarios: %.4f\n",
            toupper(covariance_aggregate), length(scen), rob$score))

## ===========================================================================
## 10. FEASIBLE PARETO FRONTIER FROM NESTED REDUCTIONS OF THE FINAL DESIGN
## ===========================================================================
## Each candidate is made only by removing Ecuador allocations from the
## seed-feasible final design. Partner sites, p-rep choices, common treatments,
## and the seed reserve remain untouched. A removal is allowed only when the TC
## still occurs elsewhere, so full population coverage is preserved.
make_reduced_design <- function(target_ecuador_entries) {
  Mr <- M_actual
  Rr <- reps_matrix
  for (s in ecuador_sites) {
    remove_n <- sum(Mr[, s]) - target_ecuador_entries
    if (remove_n < 0L) return(NULL)
    if (remove_n > 0L) {
      removable <- which(Mr[, s] == 1L & rowSums(Mr) > 1L &
                           !rownames(Mr) %in% common_treatments)
      removable <- removable[
        order(rowSums(Mr)[removable], decreasing = TRUE)]
      if (length(removable) < remove_n) return(NULL)
      drop <- removable[seq_len(remove_n)]
      Mr[drop, s] <- 0L
      Rr[drop, s] <- 0
    }
  }
  list(M = Mr, reps = Rr)
}

frontier_candidates <- sort(unique(c(
  ecuador_cap_interval[1],
  seq(ecuador_cap_interval[1], cap_ecuador, by = ecuador_cap_step),
  cap_ecuador)))
frontier_rows <- lapply(frontier_candidates, function(k) {
  dsg <- make_reduced_design(k)
  if (is.null(dsg)) return(NULL)
  obj <- design_objective(
    dsg$M, G = G, Sigma_E = Sigma_E, reps = dsg$reps,
    env_efficiency = efficiency, sigma_g2 = sigma_g2,
    sigma_e2 = variance_prior$sigma2_e,
    weights = list(gain = 1, reliability = 0, cost = 0))
  data.frame(EcuadorEntries = k,
             PhysicalPlots = sum(dsg$reps) + sum(check_overhead),
             Reliability = obj$reliability, Gain = obj$gain)
})
frontier <- do.call(rbind, frontier_rows)
is_dominated <- vapply(seq_len(nrow(frontier)), function(i) {
  any(frontier$PhysicalPlots <= frontier$PhysicalPlots[i] &
        frontier$Gain >= frontier$Gain[i] &
        (frontier$PhysicalPlots < frontier$PhysicalPlots[i] |
           frontier$Gain > frontier$Gain[i]))
}, logical(1))
frontier <- frontier[!is_dominated, , drop = FALSE]
write.csv(frontier, file.path(outdir, "feasible_pareto_frontier.csv"),
          row.names = FALSE)
print(frontier)
pdf(file.path(outdir, "pareto_frontier.pdf"), width = 7, height = 5)
plot(frontier$PhysicalPlots, frontier$Gain, type = "b", pch = 19,
     xlab = "Physical plots", ylab = "Expected genetic gain",
     main = "Seed-feasible nested frontier")
dev.off()

## ===========================================================================
## 11. VISUALISATIONS
## ===========================================================================
pdf(file.path(outdir, "pipeline_visuals.pdf"), width = 8, height = 6)
op <- par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
## (a) environment relationship heatmap
image(seq_len(nrow(D)), seq_len(ncol(D)), D[, ncol(D):1], axes = FALSE,
      xlab = "", ylab = "", main = "Multimodal environment relationship (D)")
axis(1, seq_len(nrow(D)), sites$environment, las = 2, cex.axis = 0.6)
## (b) allocation incidence (TC x site) -- subsample rows for legibility
sub <- sort(sample(nrow(M), 60))
image(seq_len(ncol(M)), seq_len(length(sub)), t(M[sub, ]), col = c("grey90","#0C4C80"),
      axes = FALSE, xlab = "site", ylab = "TC (subsample)",
      main = "Sparse allocation incidence")
axis(1, seq_len(ncol(M)), colnames(M), las = 2, cex.axis = 0.6)
## (c) per-site test-entry capacity
barplot(colSums(M), las = 2, cex.names = 0.6, col = "#3AA0FF",
        main = "Test entries per site", ylab = "entries")
## (d) mega-environment membership
plot(sites$longitude, sites$latitude,
     col = mega$membership[sites$environment], pch = 19, cex = 2,
     xlab = "longitude", ylab = "latitude",
     main = paste("Mega-environments:", mega$status))
text(sites$longitude, sites$latitude, sites$environment, pos = 3, cex = 0.55)
par(op); dev.off()

cat("\nDONE. Outputs written to", normalizePath(outdir), "\n")
cat("  - robust common-set, replication, connectivity, and seed diagnostics\n")
cat("  - allocation, fieldbooks, seed ledger, and environment summaries\n")
cat("  - weather/soil provenance, QC/imputation ledgers, modality kernels,\n")
cat("    calibration diagnostics, and Sigma_E sensitivity candidates\n")
cat("  - data-inferred mega-environment stability diagnostics and status\n")
cat("  - capacity sweep, feasible frontier, field maps, and diagnostics\n")
