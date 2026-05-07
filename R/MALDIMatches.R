# Internal matching helper
.match_maldi_to_ref <- function(maldi_mz, ref_table, tol_ppm) {
  out_list <- vector("list", length(maldi_mz))
  k <- 0L

  for (i in seq_along(maldi_mz)) {
    mz_obs <- maldi_mz[i]
    ppm_err <- abs(mz_obs - ref_table$ref_mz) / ref_table$ref_mz * 1e6
    hits <- which(ppm_err <= tol_ppm)

    for (h in hits) {
      k <- k + 1L
      out_list[[k]] <- data.frame(
        mz_MALDI = mz_obs,
        closest_value_LCMS = ref_table$ref_mz[h],
        ppm_error = round(ppm_err[h], 4),
        iso_peak_matched = ref_table$iso_peak[h],
        ref_row = ref_table$row_index[h],
        stringsAsFactors = FALSE
      )
    }
  }

  if (k == 0L) {
    warning("No matches found within ", tol_ppm, " ppm.")
    return(data.frame(
      mz_MALDI = numeric(0),
      closest_value_LCMS = numeric(0),
      ppm_error = numeric(0),
      iso_peak_matched = character(0),
      ref_row = integer(0)
    ))
  }

  do.call(rbind, out_list[seq_len(k)])
}


#' Match MALDI Peaks to LC-MS Peptide Ions
#'
#' Matches a vector of MALDI observed m/z values against the monoisotopic
#' (M+0) and most-abundant isotope peaks from an annotated LC-MS data frame.
#' All hits within the specified ppm tolerance are returned.
#'
#' @param maldi_mz Numeric vector of MALDI observed m/z values.
#' @param lcms A data frame annotated by [ComputeIon()], containing at minimum
#'   `MH_plus`, `Max_Peak_mz`, and `Max_Peak_Label`.
#' @param tolerance Numeric. Mass accuracy in ppm. Default `10`.
#' @param output_path Character path for an optional CSV of matched results.
#'   Default `NULL` (no file written).
#'
#' @return A `data.table` with columns `mz_MALDI`, `closest_value_LCMS`,
#'   `ppm_error`, `iso_peak_matched`, and all annotation columns from `lcms`.
#' @export
#' @import data.table
NULL
#'
#' @examples
#' \dontrun{
#' maldi <- read.table("T-ReX_v1.csv", header = TRUE, sep = ";",
#'                     skip = 8, fill = TRUE)
#' lcms_ann <- ComputeIon(read.csv("combined_modified_peptide.csv"))
#' matches  <- MALDIMatches(maldi$m.z, lcms_ann, tolerance = 10)
#' }
MALDIMatches <- function(maldi_mz, lcms, tolerance = 10, output_path = NULL) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required. Install it first.")
  }

  required <- c("MH_plus", "Max_Peak_mz", "Max_Peak_Label")
  missing <- setdiff(required, names(lcms))
  if (length(missing)) {
    stop(
      "lcms is missing columns: ",
      paste(missing, collapse = ", "),
      ". Run ComputeIon() first."
    )
  }

  # Build expanded reference table
  ref_mono <- data.frame(
    ref_mz = lcms$MH_plus,
    iso_peak = "M+0",
    row_index = seq_len(nrow(lcms)),
    stringsAsFactors = FALSE
  )

  is_not_mono <- lcms$Max_Peak_Label != "M+0"
  message(
    "Peptides where most-abundant peak != M+0: ",
    sum(is_not_mono),
    " of ",
    nrow(lcms)
  )

  ref_all <- if (sum(is_not_mono) > 0) {
    ref_max <- data.frame(
      ref_mz = lcms$Max_Peak_mz[is_not_mono],
      iso_peak = lcms$Max_Peak_Label[is_not_mono],
      row_index = which(is_not_mono),
      stringsAsFactors = FALSE
    )
    rbind(ref_mono, ref_max)
  } else {
    ref_mono
  }

  message(
    "Expanded reference entries: ",
    nrow(ref_all),
    " (",
    nrow(ref_mono),
    " monoisotopic + ",
    nrow(ref_all) - nrow(ref_mono),
    " most-abundant)"
  )

  # Run matching
  matches <- .match_maldi_to_ref(maldi_mz, ref_all, tolerance)

  message("Total matches: ", nrow(matches))
  message("  Matched to M+0:      ", sum(matches$iso_peak_matched == "M+0"))
  message("  Matched to non-M+0:  ", sum(matches$iso_peak_matched != "M+0"))

  # Merge back to full LCMS annotation
  matches_dt <- data.table::data.table(matches)
  LCMS_dt <- data.table::data.table(lcms)
  LCMS_dt[, row_index := .I]

  merged_dt <- merge(
    matches_dt,
    LCMS_dt,
    by.x = "ref_row",
    by.y = "row_index",
    all.x = TRUE
  )
  merged_dt[, ref_row := NULL]

  key_cols <- c(
    "mz_MALDI",
    "closest_value_LCMS",
    "ppm_error",
    "iso_peak_matched"
  )
  data.table::setcolorder(
    merged_dt,
    c(key_cols, setdiff(names(merged_dt), key_cols))
  )

  if (!is.null(output_path)) {
    write.csv(merged_dt, output_path, row.names = FALSE)
    message("Matched table written to: ", output_path)
  }

  invisible(merged_dt)
}
