#' Is this an integer?
#' 
#' @param x The putative integer
#' @param nneg Require \code{x>=0} (non-negative) too?
#' @return Logical indicating whether \code{x} is an integer
#' 
#' @keywords internal
is_integer <- function(x, nneg=FALSE){
  y <- FALSE
  if (length(x)==1 && is.numeric(x)) {
    if (x%%1==0) {
      if (x>=0 || !nneg) { y <- TRUE }
    }
  } 
  y
}

#' Is this numeric vector constant?
#' 
#' @param x The numeric vector
#' @param TOL minimum range of \code{x} to be considered non-constant.
#'  Default: \code{1e-8}
#' 
#' @return Is \code{x} constant? 
#' 
#' @keywords internal
is_constant <- function(x, TOL=1e-8) {
  abs(max(x) - min(x)) < TOL
}

#' Check design matrix
#' 
#' @param design The design matrix
#' 
#' @return The (modified) design matrix
#' 
#' @keywords internal
check_design_matrix <- function(design, T_=nrow(design)) {
  class(design) <- "numeric"
  if (identical(design, 1)) { design <- matrix(1, nrow=T_) }
  design <- as.matrix(design)
  stopifnot(nrow(design) == T_)
  # Set constant columns (intercept regressor) to 1, and scale the other columns.
  design_const_mask <- apply(design, 2, is_constant)
  design[,design_const_mask] <- 1
  design[,!design_const_mask] <- scale(design[,!design_const_mask])
  design
}

#' Robust scaling
#' 
#' Centers and scales the columns of a matrix robustly
#'
#' Centers each column on its median, and scales each column by its median
#' absolute deviation (MAD). Constant-valued columns are set to \code{NA}
#' (or removed if \code{drop_const}) and a warning is raised. If all 
#' MADs are zero, an error is raised.
#'
#' @param mat A numerical matrix.
#' @param TOL minimum MAD to consider a column non-constant.
#'  Default: \code{1e-8}
#' @param drop_const Drop
#'
#' @return The input matrix with its columns centered and scaled.
#'
#' @importFrom robustbase rowMedians
scale_med <- function(mat, TOL=1e-8, drop_const=TRUE){
  # Transpose.
  mat <- t(mat)

  #	Center.
  mat <- mat - c(rowMedians(mat, na.rm=TRUE))

  # Scale.
  mad <- 1.4826 * rowMedians(abs(mat), na.rm=TRUE)
  mad <- as.numeric(mad)
  const_mask <- mad < TOL
  if (any(const_mask)) {
    if (all(const_mask)) {
    stop("All columns are zero-variance.\n")
    } else {
      warning(paste0(
        "Warning: ", sum(const_mask),
        " constant columns (out of ", length(const_mask),
        " ). These will be removed.\n"
      ))
    }
  }
  mad <- mad[!const_mask]
  mat[const_mask,] <- NA
  mat[!const_mask,] <- mat[!const_mask,] / mad

  # Revert transpose.
  mat <- t(mat)

  if (drop_const) { mat <- mat[!const_mask,] }

  mat
}

#' Wrapper to common functions for reading NIFTIs
#' 
#' Tries \code{RNifti::readNifti}, then \code{oro.nifti::readNIfTI}. If
#'  neither package is available an error is raised.
#' 
#' @param nifti_fname The file name of the NIFTI.
#' @return The NIFTI
#' @keywords internal
read_nifti <- function(nifti_fname){
  if (requireNamespace("RNifti", quietly = TRUE)) {
    return(RNifti::readNifti(nifti_fname))
  } else if (requireNamespace("oro.nifti", quietly = TRUE)) {
    return(oro.nifti::readNIfTI(nifti_fname, reorient=FALSE))
  } else {
    stop(
      "Package \"RNifti\" or \"oro.nifti\" needed to read", nifti_fname, 
      ". Please install at least one", call. = FALSE
    )
  }
}

#' Convert to \eqn{T} by \eqn{V} matrix
#' 
#' A \code{"xifti"} is VxT, whereas \code{scrub} and \code{grayplot} accept a
#'  TxV matrix. This function calls \code{as.matrix} and transposes the data
#'  if it is a \code{"xifti"}.
#' 
#' @param x The object to coerce to a matrix
#' @param sortSub Sort subcortex by labels? Default: \code{FALSE}
#' @return x as a matrix.
#' @keywords internal
as.matrix2 <- function(x, sortSub=FALSE) {
  if (inherits(x, "xifti")) {
    if (sortSub && !is.null(x$data$subcort)) {
      x$data$subcort <- x$data$subcort[order(x$meta$subcort$labels),]
    }
    x <- t(as.matrix(x))
  } else {
    x <- as.matrix(x)
  }

  x
}