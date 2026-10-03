# Node-matrix assembly ----------------------------------------------------

# Coerce supported inputs into a numeric time x nodes matrix.
#' @keywords internal
#' @noRd
.as_node_matrix <- function(X) {
  if (methods::is(X, "PhysioExperiment") ||
      methods::is(X, "MultiPhysioExperiment"))
    return(physioNodeMatrix(X))
  if (is.matrix(X) && is.numeric(X)) return(X)
  if (is.data.frame(X)) {
    num <- vapply(X, is.numeric, logical(1))
    if (!all(num))
      stop("data.frame columns must all be numeric.", call. = FALSE)
    return(as.matrix(X))
  }
  if (is.list(X)) return(physioNodeMatrix(X))
  if (is.numeric(X))
    stop("X must have at least two columns (nodes).", call. = FALSE)
  stop("X must be a numeric matrix, data frame, or named list of vectors.",
       call. = FALSE)
}

#' Assemble a physiological node matrix
#'
#' Builds the time-by-nodes matrix that [timeDelayStability()] consumes, where
#' each node is one physiological variable (an organ system or a derived feature
#' such as an EEG band-power, heart rate, respiration amplitude, or an EMG
#' envelope). Inputs are a named list of equal-length numeric vectors, a numeric
#' matrix, or a data frame.
#'
#' In Network Physiology, coupling is assessed between the *dynamics* of these
#' variables, so nodes are typically instantaneous features sampled on a common,
#' relatively slow grid (e.g. 1 Hz), not the raw high-rate waveforms. Assemble
#' such features with the relevant `Physio*` package and pass them here.
#'
#' @param x A named list of equal-length numeric vectors, a numeric matrix
#'   (time x nodes), a data frame of numeric columns, or a `PhysioExperiment` /
#'   `MultiPhysioExperiment` object (channels/modalities become nodes).
#' @param nodes Optional character vector of node names (defaults to the names of
#'   `x`).
#' @param standardize If `TRUE`, z-normalise each column (default `FALSE`; TDS
#'   z-normalises per window internally, so this is usually unnecessary).
#' @param assay For a `PhysioExperiment`, which assay to read (default: its
#'   default assay).
#' @param channels For a `PhysioExperiment`, an optional subset of channels.
#' @param feature Per-channel feature: `"raw"` (default) or `"envelope"` (the
#'   analytic-signal amplitude — often the right node series for coupling).
#' @param target_rate Optional rate (Hz) to down-sample each channel to, by
#'   block-averaging. Required for a `MultiPhysioExperiment` whose experiments
#'   have differing sampling rates.
#'
#' @return A numeric matrix, time (rows) by nodes (columns), with column names.
#'
#' @examples
#' nodes <- list(
#'   hr   = cumsum(rnorm(500)),
#'   resp = sin(seq(0, 20 * pi, length.out = 500)),
#'   emg  = abs(rnorm(500))
#' )
#' M <- physioNodeMatrix(nodes)
#' dim(M)
#'
#' @export
physioNodeMatrix <- function(x, nodes = NULL, standardize = FALSE,
                             assay = NULL, channels = NULL,
                             feature = c("raw", "envelope"),
                             target_rate = NULL) {
  feature <- match.arg(feature)
  if (methods::is(x, "PhysioExperiment")) {
    M <- .pe_node_matrix(x, assay, channels, feature, target_rate)
  } else if (methods::is(x, "MultiPhysioExperiment")) {
    M <- .mpe_node_matrix(x, feature, target_rate)
  } else if (is.matrix(x) || is.data.frame(x)) {
    M <- .as_node_matrix(x)
  } else if (is.list(x)) {
    if (!length(x)) stop("`x` is empty.", call. = FALSE)
    ok <- vapply(x, function(v) is.numeric(v) && is.null(dim(v)), logical(1))
    if (!all(ok))
      stop("every element of `x` must be a numeric vector.", call. = FALSE)
    lens <- vapply(x, length, integer(1))
    if (length(unique(lens)) != 1L)
      stop(sprintf("all vectors must be the same length (got %s).",
                   paste(range(lens), collapse = "-")), call. = FALSE)
    M <- do.call(cbind, x)
    if (is.null(colnames(M)))
      colnames(M) <- names(x) %||% paste0("V", seq_along(x))
  } else {
    stop("`x` must be a named list, matrix, or data frame.", call. = FALSE)
  }
  if (!is.null(nodes)) {
    if (length(nodes) != ncol(M))
      stop("`nodes` length must equal the number of columns.", call. = FALSE)
    colnames(M) <- nodes
  }
  if (standardize) M <- apply(M, 2L, .znorm)
  M
}

#' @keywords internal
#' @noRd
`%||%` <- function(a, b) if (is.null(a)) b else a
