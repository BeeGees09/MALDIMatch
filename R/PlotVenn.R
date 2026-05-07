#' Plot a Pairwise Venn Diagram of MALDI vs LC-MS Matches
#'
#' Draws a Euler/Venn diagram showing the overlap between MALDI m/z values
#' and LC-MS peptide ions, using the match count from [MALDIMatches()] as
#' the intersection.
#'
#' @param maldi_mz Numeric vector of all MALDI observed m/z values (sets the
#'   MALDI circle area).
#' @param lcms A data frame annotated by [ComputeIon()], used to count total
#'   LC-MS entries via the `MH_plus` column.
#' @param matches A data frame or `data.table` returned by [MALDIMatches()].
#'   The number of rows is used as the cross-area (intersection).
#' @param output_path Character path to save a PNG. Default `NULL` (no file
#'   written).
#' @param colors Character vector of length 2 giving fill colours for the MALDI
#'   and LC-MS circles respectively. Default `c("lightblue", "lightgreen")`.
#' @param ... Additional arguments passed to
#'   [VennDiagram::draw.pairwise.venn()].
#'
#' @return The Venn diagram grob, invisibly.
#' @export
#'
#' @examples
#' \dontrun{
#' PlotVenn(maldi$m.z, lcms_ann, matches,
#'          output_path = "venn_overlap.png")
#' }
PlotVenn <- function(maldi_mz, lcms, matches,
                     output_path = NULL,
                     colors = c("lightblue", "lightgreen"),
                     ...) {

  if (!requireNamespace("VennDiagram", quietly = TRUE))
    stop("Package 'VennDiagram' is required. Install it first.")
  if (!requireNamespace("grid", quietly = TRUE))
    stop("Package 'grid' is required. Install it first.")

  total_maldi   <- length(maldi_mz)
  total_lcms    <- length(lcms$MH_plus)
  total_matches <- nrow(matches)

  venn_plot <- VennDiagram::draw.pairwise.venn(
    area1      = total_maldi,
    area2      = total_lcms,
    cross.area = total_matches,
    category   = c("MALDI m/z", "LC-MS m/z"),
    fill       = colors,
    alpha      = 0.5,
    scaled     = TRUE,
    euler.d    = TRUE,
    cat.pos    = c(0, 0),
    cat.dist   = c(0.025, 0.025),
    cex        = 1.5,
    cat.cex    = 1.5,
    cat.col    = c("blue", "darkgreen"),
    ...
  )

  grid::grid.draw(venn_plot)

  if (!is.null(output_path)) {
    grDevices::png(output_path, width = 6, height = 6,
                   units = "in", res = 300)
    grid::grid.draw(venn_plot)
    grDevices::dev.off()
    message("Venn diagram saved to: ", output_path)
  }

  invisible(venn_plot)
}

