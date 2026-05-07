# Internal: elemental deltas for supported modifications
.mod_adjust <- list(
  "P_15.9949" = c(C = 0, H = 0, N = 0, O = 1, S = 0),
  "M_15.9949" = c(C = 0, H = 0, N = 0, O = 1, S = 0),
  "C_57.0215" = c(C = 2, H = 3, N = 1, O = 1, S = 0),
  "n_42.0106" = c(C = 2, H = 2, N = 0, O = 1, S = 0)
)

# Internal: parse a modified peptide sequence and return bare sequence + elemental delta
.parse_mods <- function(mod_seq, mod_adjust = .mod_adjust) {
  delta <- c(C = 0, H = 0, N = 0, O = 0, S = 0)

  if (grepl("^n\\[", mod_seq)) {
    delta <- delta + mod_adjust[["n_42.0106"]]
    mod_seq <- sub("^n\\[[0-9.]+\\]", "", mod_seq)
  }

  hits <- gregexpr("([A-Z])\\[([0-9.]+)\\]", mod_seq, perl = TRUE)
  matches <- regmatches(mod_seq, hits)[[1]]
  for (m in matches) {
    aa <- sub("\\[.*", "", m)
    val <- sub(".*\\[([0-9.]+)\\]", "\\1", m)
    key <- paste0(aa, "_", val)
    if (key %in% names(mod_adjust)) delta <- delta + mod_adjust[[key]]
  }

  bare <- gsub("\\[[0-9.]+\\]", "", mod_seq)
  list(bare_seq = bare, delta = delta)
}

#' Compute Ion Masses and Isotopic Distributions
#'
#' For each row of an LC-MS data frame, computes the monoisotopic m/z for each
#' requested adduct, the molecular formula, and (optionally) the isotopic
#' distribution (M+0 through M+4).
#'
#' @param lcms A data frame with a modified-peptide-sequence column.
#' @param mod_col Character. Column name for modified sequences. Auto-detected
#'   if `NULL` (default).
#' @param adducts Named numeric vector of adduct mass offsets (Da). Names become
#'   column names in the output. Default: `[M+H]+`, `[M+Na]+`, `[M+K]+`.
#' @param mod_adjust Named list of elemental deltas for each modification key.
#'   Extend the built-in list to support additional modifications.
#' @param calc_isotopes Logical. If `TRUE` (default), computes the full isotopic
#'   distribution (M+0 to M+4) via `OrgMassSpecR::IsotopicDistribution`.
#'   Set to `FALSE` to skip this — much faster for large datasets when isotope
#'   columns are not needed.
#' @param output_path Character path for an optional simplified CSV. Default
#'   `NULL` (no file written).
#'
#' @return The `lcms` data frame with added columns: one m/z column per adduct,
#'   `MH_plus` (alias for `[M+H]+`), `Formula`. When `calc_isotopes = TRUE`,
#'   also adds `Max_Peak_mz`, `Max_Peak_Label`, and isotope columns
#'   `M0_mz`--`M4_mz` / `M0_%`--`M4_%`.
#' @export
#'
#' @examples
#' \dontrun{
#' lcms <- read.csv("combined_modified_peptide.csv")
#'
#' # Full calculation (default)
#' lcms_ann <- ComputeIon(lcms)
#'
#' # Fast: adduct m/z and formula only, skip isotopes
#' lcms_ann <- ComputeIon(lcms, calc_isotopes = FALSE)
#'
#' # Add ammonium adduct
#' my_adducts <- c("[M+H]+" = 1.00728, "[M+NH4]+" = 18.034)
#' lcms_ann <- ComputeIon(lcms, adducts = my_adducts)
#' }
ComputeIon <- function(
  lcms,
  mod_col = NULL,
  adducts = c("[M+H]+" = 1.00728, "[M+Na]+" = 22.98922, "[M+K]+" = 38.96316),
  mod_adjust = .mod_adjust,
  calc_isotopes = TRUE,
  output_path = NULL
) {
  if (!requireNamespace("OrgMassSpecR", quietly = TRUE)) {
    stop("Package 'OrgMassSpecR' is required. Install it first.")
  }

  # Resolve modified-sequence column
  if (is.null(mod_col)) {
    mod_col <- grep(
      "Modified.Sequence|Modified_Sequence|ModifiedSequence",
      names(lcms),
      value = TRUE
    )[1]
    if (is.na(mod_col)) {
      stop("Could not find a Modified Sequence column. Set mod_col manually.")
    }
  }

  PROTON <- 1.00728

  # Memoise slow OrgMassSpecR calls so repeated formulas are computed only once
  if (!requireNamespace("memoise", quietly = TRUE)) {
    message("Install 'memoise' for faster caching: install.packages('memoise')")
    mono_fn <- OrgMassSpecR::MonoisotopicMass
    iso_fn <- OrgMassSpecR::IsotopicDistribution
  } else {
    mono_fn <- memoise::memoise(OrgMassSpecR::MonoisotopicMass)
    iso_fn <- memoise::memoise(OrgMassSpecR::IsotopicDistribution)
  }

  n <- nrow(lcms)
  results <- vector("list", n)
  pb <- txtProgressBar(min = 0, max = n, style = 3)

  for (i in seq_len(n)) {
    ms <- lcms[[mod_col]][i]
    parsed <- .parse_mods(ms, mod_adjust)
    base_formula <- OrgMassSpecR::ConvertPeptide(parsed$bare_seq)

    adj <- list(
      C = base_formula$C + parsed$delta["C"],
      H = base_formula$H + parsed$delta["H"],
      N = base_formula$N + parsed$delta["N"],
      O = base_formula$O + parsed$delta["O"],
      S = base_formula$S + parsed$delta["S"]
    )

    fstr <- paste0(
      "C",
      adj$C,
      "H",
      adj$H,
      "N",
      adj$N,
      "O",
      adj$O,
      if (adj$S > 0) paste0("S", adj$S) else ""
    )

    neutral_mono <- mono_fn(formula = adj, charge = 0)
    adduct_mz <- neutral_mono + adducts

    if (calc_isotopes) {
      iso <- iso_fn(formula = adj, charge = 1)
      iso$mz <- iso$mz + PROTON
      iso$peak_label <- paste0("M+", seq_len(nrow(iso)) - 1)
      iso$pct <- round(iso$intensity / max(iso$intensity) * 100, 2)
      max_idx <- which.max(iso$intensity)

      results[[i]] <- list(
        adduct_mz = adduct_mz,
        formula = fstr,
        iso = iso,
        max_peak_mz = iso$mz[max_idx],
        max_peak_label = iso$peak_label[max_idx]
      )
    } else {
      results[[i]] <- list(
        adduct_mz = adduct_mz,
        formula = fstr,
        iso = NULL,
        max_peak_mz = NA_real_,
        max_peak_label = NA_character_
      )
    }

    setTxtProgressBar(pb, i)
  }
  close(pb)

  # Attach adduct columns
  safe_names <- sub("_+$", "", gsub("[^A-Za-z0-9]", "_", names(adducts)))
  for (i in seq_along(adducts)) {
    lcms[[safe_names[i]]] <- vapply(
      results,
      function(x) x$adduct_mz[i],
      numeric(1)
    )
  }
  # Backwards-compatible MH_plus alias
  mh_idx <- which(names(adducts) == "[M+H]+")
  if (length(mh_idx) == 1L) {
    lcms$MH_plus <- lcms[[safe_names[mh_idx]]]
  }

  lcms$Formula <- vapply(results, function(x) x$formula, character(1))
  lcms$Max_Peak_mz <- vapply(results, function(x) x$max_peak_mz, numeric(1))
  lcms$Max_Peak_Label <- vapply(
    results,
    function(x) x$max_peak_label,
    character(1)
  )

  if (calc_isotopes) {
    iso_mat <- t(vapply(
      results,
      function(r) {
        df <- r$iso
        n <- min(5L, nrow(df))
        c(df$mz[seq_len(n)], df$pct[seq_len(n)])
      },
      numeric(10)
    ))
    colnames(iso_mat) <- c(paste0("M", 0:4, "_mz"), paste0("M", 0:4, "_%"))
    lcms <- cbind(lcms, as.data.frame(iso_mat))
  }

  stopifnot(
    "MH_plus is empty" = !is.null(lcms$MH_plus),
    "MH_plus has NAs" = !any(is.na(lcms$MH_plus)),
    "Max_Peak_Label is empty" = !calc_isotopes ||
      length(lcms$Max_Peak_Label) == nrow(lcms)
  )

  if (!is.null(output_path)) {
    simp <- data.frame(
      Modified.Sequence = lcms[[mod_col]],
      Name = if ("Gene" %in% names(lcms)) lcms$Gene else NA,
      Formula = gsub("\\s+", "", lcms$Formula),
      `MH+` = lcms$MH_plus,
      check.names = FALSE
    )
    if (calc_isotopes) {
      simp$`M0_%` <- lcms$`M0_%`
      simp$`M1_%` <- lcms$`M1_%`
      simp$`M2_%` <- lcms$`M2_%`
      simp$`M3_%` <- lcms$`M3_%`
      simp$`M4_%` <- lcms$`M4_%`
    }
    write.csv(simp, output_path, row.names = FALSE)
    message("Simplified table written to: ", output_path)
  }

  invisible(lcms)
}
