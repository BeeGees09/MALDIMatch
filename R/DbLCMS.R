#' Create an LC-MS Database for METASPACE
#'
#' Builds a three-column METASPACE-compatible database table from an annotated
#' LC-MS data frame produced by [ComputeIon()].
#'
#' @param lcms A data frame containing at minimum the columns `Gene`,
#'   a modified-sequence column (auto-detected), and `Formula`.
#' @param mod_col Character. Name of the modified-sequence column. Auto-detected
#'   if `NULL` (default), consistent with [ComputeIon()].
#' @param output_path Character path for an optional tab-delimited TSV.
#'   Default `NULL` (no file written).
#'
#' @return A data frame with columns `id`, `name`
#'   (`Gene|Modified.Sequence`), and `formula`.
#' @export
#'
#' @examples
#' \dontrun{
#' lcms_ann <- ComputeIon(read.csv("combined_modified_peptide.csv"))
#' db <- DbLCMS(lcms_ann, output_path = "Kidney_ECM_LCMS_DB.tsv")
#' }
DbLCMS <- function(lcms, mod_col = NULL, output_path = NULL) {

  # Resolve modified-sequence column
  if (is.null(mod_col)) {
    mod_col <- grep(
      "Modified.Sequence|Modified_Sequence|ModifiedSequence",
      names(lcms), value = TRUE)[1]
    if (is.na(mod_col))
      stop("Could not find a Modified Sequence column. Set mod_col manually.")
  }

  gene_col <- if ("Gene" %in% names(lcms)) "Gene" else {
    g <- grep("^[Gg]ene", names(lcms), value = TRUE)[1]
    if (is.na(g))
      stop("Could not find a Gene column in lcms.")
    g
  }

  if (!"Formula" %in% names(lcms))
    stop("'Formula' column not found. Run ComputeIon() first.")

  LCMS_DB <- data.frame(
    id      = seq_len(nrow(lcms)),
    name    = paste0(lcms[[gene_col]], "|", lcms[[mod_col]]),
    formula = gsub("\\s+", "", lcms$Formula),
    stringsAsFactors = FALSE
  )

  if (!is.null(output_path)) {
    write.table(LCMS_DB, output_path, sep = "\t", row.names = FALSE,
                quote = FALSE)
    message("LCMS database written to: ", output_path)
  }

  invisible(LCMS_DB)
}

