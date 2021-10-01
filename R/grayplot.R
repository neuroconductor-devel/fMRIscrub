#' Grayplot of data matrix
#' 
#' Plot a matrix with \code{graphics::image}. For fMRI data, this is the "grayplot"
#'  coined by (Power, 2017). The \code{graphics} package is required.
#' 
#' @param x The \eqn{T x V} numeric data matrix, or a \code{"xifti"} object.
#'  In the plot, the \eqn{T} index will increase from left to right, and the 
#'  \eqn{V} will increase from bottom to top.
#' @param qcut Sets blackpoint at the \code{qcut} quantile, and the
#'  whitepoint at the \code{1-qcut} quantile. Default: \code{.1}. This is
#'  equivalent to setting the color range between the 10% and 90% quantiles.
#'  The quantiles are computed across the entire data matrix after any 
#'  centering or scaling. 
#' 
#'  Must be between 0 and .49. If \code{0} or \code{NULL} (default), do not 
#'  clamp the data values.
#' @param fname A \code{.pdf} (highly recommended) or \code{.png} file path
#'  to write the grayplot to. If \code{NULL} (default), return the plot directly
#'  instead of writing a file.
#' @param center,scale Center and scale the data? Cannot scale without centering.
#'  If \code{x} is fMRI data which has not otherwise been centered or scaled,
#'  it is recommended to center but not scale it (default).
#' @param colors \code{"gray255"} (default) will use a grayscale color ramp
#'  from black to white. Otherwise, this should be a character vector of
#'  color names to use. 
#' @param sortSub If \code{x} is a \code{"xifti"} object with subcortical data,
#'  should the voxels be sorted by structure alphabetically? Default: \code{TRUE}.
#' 
#'  Colors will be assigned from the lowest to the highest data value, after 
#'  any clamping of the data values by \code{qcut}.
#' @param ... Additional arguments to \code{pdf} or \code{png}, such as width
#'  and height.
#' 
#' @importFrom grDevices dev.off pdf png gray.colors
#' 
#' @return The image or \code{NULL, invisibly} if a file was written.
#' 
#' @section References:
#'  \itemize{
#'    \item{Power, J. D. A simple but useful way to assess fMRI scan qualities. NeuroImage 154, 150-158 (2017).}
#' }
#' 
#' @export
#' 
grayplot <- function(
  x, qcut=.1, fname=NULL, center=TRUE, scale=FALSE, colors="gray255", sortSub=TRUE, ...){

  if (!requireNamespace("graphics", quietly = TRUE)) {
    stop("Package \"graphics\" needed since `svd` failed. Please install it.", call. = FALSE)
  }

  # Get T x V matrix.
  x <- as.matrix2(x, sortSub=sortSub)

  # Center or scale.
  if (center | scale) { x <- scale(x, center=center, scale=scale) }

  # Orient correctly for `image`
  x <- x[,seq(ncol(x), 1)]

  # Clamp to quantile cutoffs.
  if (!is.null(qcut)) {
    qcut <- as.numeric(qcut)
    if (qcut > 1e-8) {
      stopifnot(qcut < .49)
      x[x < quantile(x, qcut)] <- quantile(x, qcut)
      x[x > quantile(x, 1-qcut)] <- quantile(x, 1-qcut)
    }
  }

  # Check colors.
  if (identical(colors, "gray255")) { 
    colors <- gray.colors(255, start=0, end=1)
  }

  # Check file name
  if (!is.null(fname)) {
    fname <- as.character(fname)
    ftype <- ifelse(endsWith(fname, ".png"), "png", "pdf")
    if (ftype == "pdf" && !endsWith(fname, ".pdf")) {
      fname <- paste0(fname, ".pdf")
    }

    switch(ftype, png=png, pdf=pdf)(fname, ...)
  }

  graphics::image(
    x, 
    axes=FALSE, frame.plot = FALSE, xlab="", ylab="", ann=FALSE, 
    useRaster=TRUE, col=colors
  )

  if (!is.null(fname)) { 
    dev.off()
  }

  invisible(NULL)
}