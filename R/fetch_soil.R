#' Retrieve SoilGrids soil properties at site coordinates
#'
#' @description
#' Production-grade access to ISRIC SoilGrids predictions at a set of site
#' coordinates. The public REST "properties/query" endpoint is a beta service
#' that is rate-limited and intermittently returns only null values, so it is
#' fragile for a pipeline. ISRIC instead recommends the OGC **Web Coverage
#' Service (WCS)** for extracting spatial subsets, and **WebDAV** for direct
#' access to the global VRT layers. `fetch_soilgrids()` supports both of those
#' robust backends (plus the legacy REST endpoint and purely local files) behind
#' a single interface.
#'
#' @details
#' SoilGrids stores every property as scaled integers to save space; the raw
#' values must be divided by a property-specific conversion factor to obtain
#' conventional units (e.g. clay is stored in g/kg and divided by 10 to reach
#' g/100g = \%). `fetch_soilgrids()` applies these factors automatically unless
#' `apply_conversion = FALSE`. The conversion factors are:
#' \tabular{lll}{
#'   Property \tab Mapped units \tab Factor \cr
#'   bdod \tab cg/cm3 \tab 100 \cr
#'   nitrogen \tab cg/kg \tab 100 \cr
#'   cec, cfvo, clay, sand, silt, soc, ocd, ocs, phh2o \tab various \tab 10 \cr
#'   wv0010, wv0033, wv1500 \tab 10^-3 \tab 10
#' }
#'
#' Note that SoilGrids is a single *modelled snapshot*, not a time series. Some
#' soil properties change slowly (texture, bulk density) while others change with
#' time and management (organic carbon, nitrogen, pH, nutrients, salinity), and
#' soil moisture changes daily. For time-resolved soil use: the daily NASA POWER
#' `GWET*` soil-wetness variables (via [fetch_weather_series()]), satellite soil
#' moisture (SMAP, ESA-CCI, Copernicus), point profiles (ISRIC WoSIS), or -- best
#' -- your own dated field measurements, which you can pass directly as the
#' `soil` argument of [build_enviromic_covariates()] (including one set per year).
#'
#' Predictions exist for six standard depth intervals
#' (`0-5cm`, `5-15cm`, `15-30cm`, `30-60cm`, `60-100cm`, `100-200cm`) and four
#' prediction statistics (`mean`, `Q0.05`, `Q0.5`, `Q0.95`). SoilGrids native
#' data are on the Interrupted Goode Homolosine grid (EPSG:152160) at 250 m
#' resolution; the WCS and WebDAV backends reproject site coordinates onto that
#' grid with [terra::project()] and read only the cells needed via GDAL's
#' `/vsicurl` virtual file system, so no full-grid download occurs.
#'
#' @param sites A data frame with columns `environment`, `latitude`, `longitude`
#'   (WGS84 decimal degrees).
#' @param backend One of `"wcs"` (default; OGC Web Coverage Service, needs
#'   \pkg{terra}), `"webdav"` (crop the global VRT over `/vsicurl` with
#'   \pkg{gdalUtilities} then read locally; falls back to `"wcs"` if that
#'   package is absent), `"rest"` (legacy beta properties/query endpoint, needs
#'   \pkg{jsonlite}), or `"local"` (read from local raster files).
#' @param properties Character vector of SoilGrids property names to retrieve.
#' @param depth One or more of the six standard depth intervals.
#' @param quantile One or more prediction statistics: `"mean"`, `"Q0.05"`,
#'   `"Q0.5"` or `"Q0.95"`. Retrieving Q0.05/Q0.5/Q0.95 makes model uncertainty
#'   available to [soil_profile_features()].
#' @param local_paths Named list/vector mapping each property to a local raster
#'   file, required when `backend = "local"`.
#' @param buffer_m Half-width (metres) of the window requested around each site
#'   on the 250 m Homolosine grid (WCS `GetCoverage` subset, or WebDAV
#'   `gdal_translate` crop). Default 1000.
#' @param resolution_m Output cell size (metres) for the WebDAV crop. Default
#'   250 (SoilGrids native resolution).
#' @param apply_conversion Logical; divide raw values by the SoilGrids
#'   conversion factor to obtain conventional units. Default `TRUE`.
#' @param spatial_summary `"point"` extracts the cell containing the coordinate;
#'   `"mean_sd"` returns the neighbourhood mean and SD from the requested
#'   buffer, representing local spatial heterogeneity.
#'
#' @return A data frame with one row per site (row names = `environment`) and one
#'   column per property for scalar point requests. Profile/uncertainty requests
#'   use `property__depth__quantile__summary` names. Returns `NULL` if nothing
#'   usable was retrieved (with a
#'   warning naming the cause). The `"rest"`, `"wcs"`, `"webdav"` and `"local"`
#'   raster backends require the \pkg{terra} package (except `"rest"`, which
#'   requires \pkg{jsonlite}).
#'
#' @references
#' ISRIC SoilGrids FAQ and WCS/WebDAV access guides, \url{https://www.isric.org/},
#' \url{https://maps.isric.org/}, \url{https://files.isric.org/soilgrids/latest/data/}.
#'
#' @seealso [build_enviromic_covariates()], [build_environment_relationship()].
#' @examples
#' \dontrun{
#' sites <- data.frame(environment = c("E1", "E2"),
#'                     latitude = c(9.06, -1.29), longitude = c(7.49, 36.82))
#' # Robust OGC Web Coverage Service backend (default), topsoil mean values:
#' soil <- fetch_soilgrids(sites, backend = "wcs",
#'                         properties = c("clay", "sand", "phh2o", "soc"),
#'                         depth = "0-5cm", quantile = "mean")
#' # Direct VRT access over WebDAV instead:
#' soil <- fetch_soilgrids(sites, backend = "webdav", properties = "clay")
#' }
#' @export
fetch_soilgrids <- function(sites,
                            backend = c("wcs", "webdav", "rest", "local"),
                            properties = c("bdod", "cec", "cfvo", "clay",
                                           "nitrogen", "ocd", "phh2o", "sand",
                                           "silt", "soc"),
                            depth = "0-5cm", quantile = "mean",
                            local_paths = NULL, buffer_m = 1000,
                            resolution_m = 250, apply_conversion = TRUE,
                            spatial_summary = c("point", "mean_sd")) {
  backend <- match.arg(backend)
  spatial_summary <- match.arg(spatial_summary)
  .validate_soil_request(depth, quantile, properties)

  if (backend == "rest") {
    if (length(depth) != 1L || depth != "0-5cm" ||
        length(quantile) != 1L || quantile != "mean" ||
        spatial_summary != "point")
      stop("The legacy REST backend supports only depth = '0-5cm', ",
           "quantile = 'mean', and spatial_summary = 'point'. Use WCS, ",
           "WebDAV, or local rasters for soil profiles and uncertainty.")
    s <- .fetch_soil_soilgrids(sites)          # legacy beta endpoint
    return(s)
  }

  if (!requireNamespace("terra", quietly = TRUE)) {
    if (backend == "local")
      stop("Package 'terra' is required for backend = 'local'.")
    if (length(depth) > 1L || length(quantile) > 1L ||
        spatial_summary != "point") {
      warning("Package 'terra' is required for multi-depth, multi-quantile, ",
              "or neighbourhood SoilGrids retrieval; soil was not fetched.",
              call. = FALSE)
      return(NULL)
    }
    warning("Package 'terra' is required for soil backend '", backend,
            "'; falling back to the REST endpoint.", call. = FALSE)
    return(.fetch_soil_soilgrids(sites))
  }
  # The WebDAV backend crops the remote VRT with gdalUtilities; if it is not
  # installed, fall back to the (terra-only) WCS backend.
  if (backend == "webdav" &&
      !requireNamespace("gdalUtilities", quietly = TRUE)) {
    warning("Package 'gdalUtilities' is required for the WebDAV backend; ",
            "falling back to the WCS backend.", call. = FALSE)
    backend <- "wcs"
  }
  if (!all(c("latitude", "longitude") %in% names(sites))) {
    warning("`sites` needs `latitude` and `longitude` for soil fetch.",
            call. = FALSE)
    return(NULL)
  }
  combos <- expand.grid(
    property = properties, depth = depth, quantile = quantile,
    stringsAsFactors = FALSE
  )
  if (backend == "local" && is.null(local_paths))
    stop("`local_paths` is required when backend = 'local'.")
  if (backend == "local") {
    exact <- paste(combos$property, combos$depth, combos$quantile, sep = "__")
    legacy_ok <- nrow(combos) == length(properties) &&
      length(depth) == 1L && length(quantile) == 1L &&
      all(properties %in% names(local_paths))
    if (!legacy_ok && !all(exact %in% names(local_paths)))
      stop("For multi-depth/quantile local retrieval, `local_paths` must name ",
           "every property__depth__quantile combination.")
  }

  envs   <- as.character(sites$environment)
  pts_ll <- terra::vect(
    data.frame(longitude = as.numeric(sites$longitude),
               latitude  = as.numeric(sites$latitude)),
    geom = c("longitude", "latitude"), crs = "EPSG:4326")

  # For the SoilGrids-native backends (wcs, webdav) project the SITE POINTS to
  # the ISRIC Interrupted Goode Homolosine grid once. We never derive the CRS
  # from the downloaded raster: the WCS GeoTIFF response carries no CRS that
  # terra recognises, so terra::crs(raster) is invalid and projecting to it
  # fails. All geometry is anchored on the known Homolosine PROJ string instead.
  crds <- if (backend %in% c("wcs", "webdav")) .project_sites_igh(pts_ll) else NULL

  cols <- lapply(seq_len(nrow(combos)), function(k) {
    p <- combos$property[k]; dep <- combos$depth[k]; qua <- combos$quantile[k]
    path_key <- paste(p, dep, qua, sep = "__")
    local_path <- if (path_key %in% names(local_paths))
      local_paths[[path_key]] else local_paths[[p]]
    vals <- switch(backend,
      wcs    = .soil_extract_wcs(p, dep, qua, crds, buffer_m,
                                 spatial_summary),
      webdav = .soil_extract_webdav(p, dep, qua, crds,
                                    buffer_m, resolution_m, spatial_summary),
      local  = tryCatch(
        .soil_extract_local(local_path, pts_ll, buffer_m, spatial_summary),
        error = function(e) {
          warning(sprintf("Soil fetch failed for '%s' (local): %s",
                          p, conditionMessage(e)), call. = FALSE)
          matrix(NA_real_, length(envs),
                 if (spatial_summary == "point") 1L else 2L)
        }))
    if (apply_conversion) vals <- vals / .soilgrids_conversion_factor(p)
    vals <- as.matrix(vals)
    if (spatial_summary == "point") colnames(vals) <- "point"
    else colnames(vals) <- c("mean", "sd")
    vals
  })

  legacy_names <- nrow(combos) == length(properties) &&
    length(depth) == 1L && length(quantile) == 1L &&
    spatial_summary == "point"
  named_cols <- lapply(seq_along(cols), function(k) {
    M <- cols[[k]]
    prefix <- if (legacy_names) combos$property[k] else paste(
      combos$property[k], combos$depth[k], combos$quantile[k], sep = "__"
    )
    colnames(M) <- if (legacy_names) prefix else
      paste(prefix, colnames(M), sep = "__")
    M
  })
  out <- as.data.frame(do.call(cbind, named_cols), stringsAsFactors = FALSE,
                       check.names = FALSE)
  rownames(out) <- envs
  attr(out, "soilgrids_request") <- combos
  attr(out, "spatial_summary") <- spatial_summary
  attr(out, "provenance") <- data.frame(
    source = "ISRIC SoilGrids",
    backend = backend,
    model_release = "latest",
    retrieved_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    native_resolution_m = 250,
    buffer_m = buffer_m,
    requested_resolution_m = resolution_m,
    spatial_summary = spatial_summary,
    depths = paste(depth, collapse = ";"),
    quantiles = paste(quantile, collapse = ";"),
    properties = paste(properties, collapse = ";"),
    apply_conversion = apply_conversion,
    stringsAsFactors = FALSE
  )

  if (all(vapply(out, function(col) all(is.na(col)), logical(1)))) {
    warning("SoilGrids returned only null/NA values (service may be degraded ",
            "or points fall outside coverage); skipping soil.", call. = FALSE)
    return(NULL)
  }
  out
}


# ---- constants --------------------------------------------------------------

# Interrupted Goode Homolosine projection SoilGrids is stored in (EPSG:152160).
.SOILGRIDS_IGH <- "+proj=igh +lat_0=0 +lon_0=0 +datum=WGS84 +units=m +no_defs"

.SOILGRIDS_DEPTHS <- c("0-5cm", "5-15cm", "15-30cm",
                       "30-60cm", "60-100cm", "100-200cm")

.SOILGRIDS_QUANTILES <- c("mean", "Q0.05", "Q0.5", "Q0.95")

.SOILGRIDS_PROPERTIES <- c("bdod", "cec", "cfvo", "clay", "nitrogen",
                           "ocd", "ocs", "phh2o", "sand", "silt", "soc",
                           "wv0010", "wv0033", "wv1500")


# ---- helpers ----------------------------------------------------------------

# Divide raw SoilGrids integers by these to reach conventional units.
.soilgrids_conversion_factor <- function(prop) {
  cf <- c(bdod = 100, nitrogen = 100,
          cec = 10, cfvo = 10, clay = 10, sand = 10, silt = 10,
          soc = 10, ocd = 10, ocs = 10, phh2o = 10,
          wv0010 = 10, wv0033 = 10, wv1500 = 10)
  f <- cf[prop]
  if (is.na(f)) 1 else unname(f)
}

.validate_soil_request <- function(depth, quantile, properties) {
  if (!length(depth) || any(!depth %in% .SOILGRIDS_DEPTHS))
    stop("`depth` values must be among: ",
         paste(.SOILGRIDS_DEPTHS, collapse = ", "), ".")
  if (!length(quantile) || any(!quantile %in% .SOILGRIDS_QUANTILES))
    stop("`quantile` must be one of: ",
         paste(.SOILGRIDS_QUANTILES, collapse = ", "), ".")
  unknown <- setdiff(properties, .SOILGRIDS_PROPERTIES)
  if (length(unknown))
    stop("Unknown SoilGrids propert", if (length(unknown) > 1) "ies: " else "y: ",
         paste(unknown, collapse = ", "), ".")
  invisible(TRUE)
}

# Build a WCS 2.0.1 GetCoverage URL for a property over a Homolosine bbox.
.wcs_getcoverage_url <- function(prop, depth, quantile, bbox) {
  paste0(
    sprintf("https://maps.isric.org/mapserv?map=/map/%s.map", prop),
    "&SERVICE=WCS&VERSION=2.0.1&REQUEST=GetCoverage",
    sprintf("&COVERAGEID=%s_%s_%s", prop, depth, quantile),
    "&FORMAT=GEOTIFF_INT16",
    sprintf("&SUBSET=X(%f,%f)", bbox[1], bbox[2]),
    sprintf("&SUBSET=Y(%f,%f)", bbox[3], bbox[4]),
    "&SUBSETTINGCRS=http://www.opengis.net/def/crs/EPSG/0/152160",
    "&OUTPUTCRS=http://www.opengis.net/def/crs/EPSG/0/152160")
}

# GDAL /vsicurl root for the SoilGrids WebDAV store, with retries enabled and
# directory listing disabled (per the ISRIC WebDAV access guide).
.SOILGRIDS_WEBDAV_ROOT <- paste0(
  "/vsicurl?max_retry=3&retry_delay=1&list_dir=no&",
  "url=https://files.isric.org/soilgrids/latest/data/")

# Build a GDAL /vsicurl path to a SoilGrids WebDAV VRT layer.
.webdav_vrt_url <- function(prop, depth, quantile) {
  paste0(.SOILGRIDS_WEBDAV_ROOT, prop, "/", prop, "_", depth, "_",
         quantile, ".vrt")
}

# Project WGS84 site points onto the SoilGrids Homolosine grid and return their
# x/y coordinates (metres). Fails loudly if projection produces non-finite
# coordinates.
.project_sites_igh <- function(pts_ll) {
  pts_igh <- terra::project(pts_ll, .SOILGRIDS_IGH)
  crds <- terra::crds(pts_igh, df = TRUE)          # columns x, y
  if (nrow(crds) != nrow(pts_ll) || any(!is.finite(as.matrix(crds))))
    stop("Projection of site coordinates to the SoilGrids Homolosine CRS failed.",
         call. = FALSE)
  crds
}

# A single point expressed directly in the Homolosine CRS.
.soil_point_igh <- function(x, y) {
  terra::vect(data.frame(x = x, y = y), geom = c("x", "y"),
              crs = .SOILGRIDS_IGH)
}

# WCS backend: one small GetCoverage window per site, extracted at the projected
# Homolosine point. Errors are caught per site so one bad point cannot void the
# whole property.
.soil_extract_wcs <- function(prop, depth, quantile, crds, buffer_m,
                              spatial_summary = "point") {
  out <- lapply(seq_len(nrow(crds)), function(i) {
    x <- crds$x[i]; y <- crds$y[i]
    tryCatch({
      url <- .wcs_getcoverage_url(prop, depth, quantile,
                                  c(x - buffer_m, x + buffer_m,
                                    y - buffer_m, y + buffer_m))
      r <- terra::rast(url)
      if (spatial_summary == "point") {
        v <- terra::extract(r, .soil_point_igh(x, y), ID = FALSE)
        ans <- as.numeric(v[[1]])[1]
        if (is.finite(ans) && ans <= -32000) ans <- NA_real_
        c(point = ans)
      } else {
        v <- as.numeric(terra::values(r))
        v <- v[is.finite(v) & v > -32000]
        c(mean = if (length(v)) mean(v) else NA_real_,
          sd = if (length(v) > 1L) stats::sd(v) else NA_real_)
      }
    }, error = function(e) {
      warning(sprintf("Soil fetch failed for '%s' (wcs) at site %d: %s",
                      prop, i, conditionMessage(e)), call. = FALSE)
      if (spatial_summary == "point") c(point = NA_real_) else
        c(mean = NA_real_, sd = NA_real_)
    })
  })
  do.call(rbind, out)
}

# WebDAV backend (ISRIC official workflow): for each site, crop the remote VRT
# to a small local GeoTIFF with gdal_translate (projwin in native Homolosine),
# then read the local file and extract the point. Reading a single cell from the
# global remote VRT directly is unreliable; the local crop is the supported path.
.soil_extract_webdav <- function(prop, depth, quantile, crds,
                                 buffer_m, resolution_m,
                                 spatial_summary = "point") {
  vrt <- .webdav_vrt_url(prop, depth, quantile)
  out_all <- lapply(seq_len(nrow(crds)), function(i) {
    x <- crds$x[i]; y <- crds$y[i]
    tryCatch({
      out <- tempfile(pattern = "soilgrids_webdav_", fileext = ".tif")
      on.exit(unlink(out), add = TRUE)
      gdalUtilities::gdal_translate(
        src_dataset = vrt, dst_dataset = out,
        tr = c(resolution_m, resolution_m),
        # projwin order: ulx, uly, lrx, lry (metres, Homolosine)
        projwin = c(x - buffer_m, y + buffer_m, x + buffer_m, y - buffer_m),
        projwin_srs = .SOILGRIDS_IGH,
        co = c("TILED=YES", "COMPRESS=DEFLATE"))
      if (!file.exists(out) || is.na(file.info(out)$size) ||
          file.info(out)$size <= 0)
        stop("GDAL did not create a valid local SoilGrids raster.")
      r <- terra::rast(out)
      if (spatial_summary == "point") {
        v <- terra::extract(r, .soil_point_igh(x, y), ID = FALSE)
        ans <- as.numeric(v[[1]])[1]
        if (is.finite(ans) && ans <= -32000) ans <- NA_real_
        c(point = ans)
      } else {
        v <- as.numeric(terra::values(r))
        v <- v[is.finite(v) & v > -32000]
        c(mean = if (length(v)) mean(v) else NA_real_,
          sd = if (length(v) > 1L) stats::sd(v) else NA_real_)
      }
    }, error = function(e) {
      warning(sprintf("Soil fetch failed for '%s' (webdav) at site %d: %s",
                      prop, i, conditionMessage(e)), call. = FALSE)
      if (spatial_summary == "point") c(point = NA_real_) else
        c(mean = NA_real_, sd = NA_real_)
    })
  })
  do.call(rbind, out_all)
}

# Local backend: the file carries its own (valid) CRS, so here we DO project the
# WGS84 points into the raster's CRS.
.soil_extract_local <- function(path, pts_ll, buffer_m = 1000,
                                spatial_summary = "point") {
  r <- terra::rast(path)
  p <- terra::project(pts_ll, terra::crs(r))
  if (spatial_summary == "point") {
    v <- terra::extract(r, p, ID = FALSE)
    return(matrix(as.numeric(v[[1]]), ncol = 1L,
                  dimnames = list(NULL, "point")))
  }
  av <- terra::extract(r, p, buffer = buffer_m, fun = mean,
                       na.rm = TRUE, ID = FALSE)
  sdv <- terra::extract(r, p, buffer = buffer_m, fun = stats::sd,
                        na.rm = TRUE, ID = FALSE)
  cbind(mean = as.numeric(av[[1]]), sd = as.numeric(sdv[[1]]))
}


#' Summarise multi-depth SoilGrids predictions over the crop root zone
#'
#' @description
#' Converts the multi-depth output of [fetch_soilgrids()] into root-zone
#' covariates. Layer values are weighted by their thickness within
#' `root_depth_cm`. Prediction-interval widths (Q0.95 - Q0.05), neighbourhood
#' heterogeneity, plant-available water (`wv0033 - wv1500`), and isometric
#' log-ratio coordinates for the sand-silt-clay composition are included where
#' the necessary source variables exist.
#'
#' @param soil Multi-depth `fetch_soilgrids()` output.
#' @param root_depth_cm Effective rooting depth in cm.
#' @param texture_zero Small positive replacement used before the texture
#'   compositional log-ratio transform.
#'
#' @return An environment-by-feature numeric matrix with source-layer metadata
#'   and quantile/AWC/texture audit attributes.
#' @export
soil_profile_features <- function(soil, root_depth_cm = 100,
                                  texture_zero = 0.01) {
  S <- as.data.frame(soil, check.names = FALSE)
  if (is.null(rownames(S)))
    stop("`soil` needs environment row names.")
  if (!is.numeric(root_depth_cm) || length(root_depth_cm) != 1L ||
      !is.finite(root_depth_cm) || root_depth_cm <= 0 || root_depth_cm > 200)
    stop("`root_depth_cm` must be one number in (0, 200].")
  if (!is.numeric(texture_zero) || length(texture_zero) != 1L ||
      !is.finite(texture_zero) || texture_zero <= 0)
    stop("`texture_zero` must be finite and positive.")

  parsed <- lapply(names(S), .parse_soil_profile_name)
  keep <- !vapply(parsed, is.null, logical(1))
  if (!any(keep))
    stop("No multi-depth SoilGrids columns found. Fetch more than one depth or ",
         "quantile so columns follow property__depth__quantile__summary.")
  meta <- do.call(rbind, parsed[keep])
  meta$column <- names(S)[keep]
  bounds <- do.call(rbind, lapply(meta$depth, .soil_depth_bounds))
  meta$lower_cm <- bounds[, 1L]; meta$upper_cm <- bounds[, 2L]
  meta$thickness_in_root <- pmax(
    0, pmin(meta$upper_cm, root_depth_cm) - meta$lower_cm
  )
  meta <- meta[meta$thickness_in_root > 0, , drop = FALSE]
  if (!nrow(meta))
    stop("None of the supplied soil layers overlaps the requested root depth.")

  groups <- split(seq_len(nrow(meta)),
                  interaction(meta$property, meta$quantile, meta$summary,
                              drop = TRUE, lex.order = TRUE))
  out <- list()
  for (ii in groups) {
    mm <- meta[ii, , drop = FALSE]
    vals <- data.matrix(S[, mm$column, drop = FALSE])
    w <- mm$thickness_in_root
    ans <- apply(vals, 1L, function(z) {
      ok <- is.finite(z) & is.finite(w) & w > 0
      if (!any(ok)) return(NA_real_)
      stats::weighted.mean(z[ok], w[ok])
    })
    qname <- gsub("\\.", "", tolower(mm$quantile[1L]))
    sname <- mm$summary[1L]
    if (sname == "point" && mm$quantile[1L] == "mean") {
      suffix <- "root_mean"
    } else {
      suffix <- paste0("root_", qname, "_", sname)
    }
    out[[paste(mm$property[1L], suffix, sep = "_")]] <- ans
  }
  F <- as.data.frame(out, check.names = FALSE, stringsAsFactors = FALSE)
  rownames(F) <- rownames(S)
  audit_rows <- list()

  props <- unique(meta$property)
  for (p in props) {
    lo <- grep(paste0("^", p, "_root_q005_(point|mean)$"), names(F),
               value = TRUE)
    hi <- grep(paste0("^", p, "_root_q095_(point|mean)$"), names(F),
               value = TRUE)
    if (length(lo) && length(hi)) {
      width <- F[[hi[1L]]] - F[[lo[1L]]]
      crossed <- is.finite(width) & width < 0
      audit_rows[[length(audit_rows) + 1L]] <- data.frame(
        check = paste0(p, "_quantile_order"),
        n_flagged = sum(crossed),
        status = if (any(crossed)) "failed_set_missing" else "ok",
        stringsAsFactors = FALSE
      )
      if (any(crossed))
        warning("Crossing SoilGrids prediction quantiles for property '", p,
                "' were set to missing.", call. = FALSE)
      width[crossed] <- NA_real_
      F[[paste0(p, "_root_prediction_width")]] <- width
    }
  }
  value_col <- function(p) {
    candidates <- c(
      paste0(p, "_root_mean"),
      paste0(p, "_root_mean_mean"),
      paste0(p, "_root_q05_point"),
      paste0(p, "_root_q05_mean")
    )
    intersect(candidates, names(F))[1L]
  }
  fc <- value_col("wv0033"); wp <- value_col("wv1500")
  if (!is.na(fc) && !is.na(wp)) {
    awc <- F[[fc]] - F[[wp]]
    invalid_awc <- is.finite(awc) & awc < 0
    audit_rows[[length(audit_rows) + 1L]] <- data.frame(
      check = "available_water_capacity_nonnegative",
      n_flagged = sum(invalid_awc),
      status = if (any(invalid_awc)) "failed_set_missing" else "ok",
      stringsAsFactors = FALSE
    )
    if (any(invalid_awc))
      warning("Negative root-zone available-water values were set to missing.",
              call. = FALSE)
    awc[invalid_awc] <- NA_real_
    F$available_water_capacity_root <- awc
  }

  sand <- value_col("sand"); silt <- value_col("silt"); clay <- value_col("clay")
  if (!anyNA(c(sand, silt, clay))) {
    a <- pmax(texture_zero, F[[sand]])
    b <- pmax(texture_zero, F[[silt]])
    clay_part <- pmax(texture_zero, F[[clay]])
    total <- a + b + clay_part
    audit_rows[[length(audit_rows) + 1L]] <- data.frame(
      check = "texture_closure_before_ilr",
      n_flagged = sum(abs(total - 100) > 10, na.rm = TRUE),
      status = "closed_to_unit_sum",
      stringsAsFactors = FALSE
    )
    a <- a / total; b <- b / total; clay_part <- clay_part / total
    F$texture_ilr_sand_vs_silt <- sqrt(1 / 2) * log(a / b)
    F$texture_ilr_coarse_vs_clay <-
      sqrt(2 / 3) * log(sqrt(a * b) / clay_part)
  }
  M <- data.matrix(F)
  attr(M, "soil_layers") <- meta
  attr(M, "root_depth_cm") <- root_depth_cm
  attr(M, "soil_audit") <- if (length(audit_rows))
    do.call(rbind, audit_rows) else data.frame()
  M
}

.parse_soil_profile_name <- function(x) {
  parts <- strsplit(x, "__", fixed = TRUE)[[1L]]
  if (length(parts) != 4L) return(NULL)
  if (!parts[2L] %in% .SOILGRIDS_DEPTHS ||
      !parts[3L] %in% .SOILGRIDS_QUANTILES ||
      !parts[4L] %in% c("point", "mean", "sd"))
    return(NULL)
  data.frame(property = parts[1L], depth = parts[2L],
             quantile = parts[3L], summary = parts[4L],
             stringsAsFactors = FALSE)
}

.soil_depth_bounds <- function(x) {
  as.numeric(strsplit(sub("cm$", "", x), "-", fixed = TRUE)[[1L]])
}
