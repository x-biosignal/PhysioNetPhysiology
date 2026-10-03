# PhysioExperiment / MultiPhysioExperiment integration --------------------
# Turn ecosystem data objects into node matrices. Channels of a PhysioExperiment
# (or one node-set per experiment in a MultiPhysioExperiment) become nodes. Kept
# in Suggests: these helpers only run when a caller passes such an object, at
# which point PhysioCore / PhysioCrossModal are necessarily loaded.

# block-average down-sampling to a target rate
#' @keywords internal
#' @noRd
.downsample_mat <- function(M, fs, target) {
  if (is.na(fs) || is.null(target) || target >= fs) return(M)
  factor <- max(1L, as.integer(round(fs / target)))
  if (factor <= 1L) return(M)
  n <- nrow(M); nb <- n %/% factor
  if (nb < 1L) return(M)
  keep <- seq_len(nb * factor)
  grp <- rep(seq_len(nb), each = factor)
  cn <- colnames(M)
  out <- vapply(seq_len(ncol(M)),
                function(j) as.numeric(tapply(M[keep, j], grp, mean)),
                numeric(nb))
  colnames(out) <- cn
  out
}

#' @keywords internal
#' @noRd
.pe_node_matrix <- function(pe, assay = NULL, channels = NULL,
                            feature = "raw", target_rate = NULL) {
  an <- assay %||% PhysioCore::defaultAssay(pe)
  M <- as.matrix(SummarizedExperiment::assay(pe, an))     # time x channels
  cn <- tryCatch(PhysioCore::channelNames(pe), error = function(e) NULL)
  if (!is.null(cn) && length(cn) == ncol(M)) colnames(M) <- cn
  if (!is.null(channels)) M <- M[, channels, drop = FALSE]
  if (feature == "envelope") {
    cc <- colnames(M)
    M <- vapply(seq_len(ncol(M)), function(j) Mod(.analytic(M[, j])),
                numeric(nrow(M)))               # amplitude envelope
    colnames(M) <- cc
  }
  fs <- tryCatch(PhysioCore::samplingRate(pe), error = function(e) NA_real_)
  M <- .downsample_mat(M, fs, target_rate)
  M
}

#' @keywords internal
#' @noRd
.mpe_node_matrix <- function(mpe, feature = "raw", target_rate = NULL) {
  exps <- PhysioCrossModal::experiments(mpe)
  if (!length(exps)) stop("MultiPhysioExperiment has no experiments.",
                          call. = FALSE)
  nm <- names(exps) %||% paste0("exp", seq_along(exps))
  rates <- vapply(exps, function(e)
    tryCatch(PhysioCore::samplingRate(e), error = function(err) NA_real_),
    numeric(1))
  if (is.null(target_rate) && length(unique(stats::na.omit(rates))) > 1L)
    stop("experiments have different sampling rates; set `target_rate`.",
         call. = FALSE)
  mats <- lapply(seq_along(exps), function(i) {
    m <- .pe_node_matrix(exps[[i]], NULL, NULL, feature, target_rate)
    base <- colnames(m) %||% paste0("V", seq_len(ncol(m)))
    colnames(m) <- paste(nm[i], base, sep = ".")
    m
  })
  L <- min(vapply(mats, nrow, integer(1)))
  do.call(cbind, lapply(mats, function(m) m[seq_len(L), , drop = FALSE]))
}
