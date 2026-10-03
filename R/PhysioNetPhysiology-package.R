#' PhysioNetPhysiology: Network Physiology for the Physio ecosystem
#'
#' Dynamic-interaction analysis of coupled physiological systems. The package
#' centres on **Time Delay Stability (TDS)** — a method that detects stable
#' temporal couplings between signals by tracking, over sliding windows, the lag
#' that maximises their cross-correlation and rewarding periods where that lag
#' stays put. Stable couplings become links in an *organ-interaction network*
#' whose edges are significance-tested against phase-randomised surrogates, and
#' the package quantifies how such networks **reconfigure across physiological
#' states** (wake / sleep stages, rest / task, ...).
#'
#' Core entry points:
#' \itemize{
#'   \item [timeDelayStability()] — the TDS engine over a multivariate series.
#'   \item [tdsSurrogateThreshold()] — surrogate-based significance threshold.
#'   \item [tdsNetwork()] — build a thresholded organ-interaction network.
#'   \item [tdsNetworkByState()] — state-resolved networks + reconfiguration.
#'   \item [physioNodeMatrix()] — assemble a node matrix of physiological systems.
#' }
#'
#' @references
#' Bashan A, Bartsch RP, Kantelhardt JW, Havlin S, Ivanov PC (2012).
#' Network physiology reveals relations between network topology and
#' physiological function. \emph{Nature Communications} 3:702.
#' \doi{10.1038/ncomms1705}
#'
#' Bartsch RP, Liu KKL, Bashan A, Ivanov PC (2015). Network physiology: how
#' organ systems dynamically interact. \emph{PLoS ONE} 10(11):e0142143.
#' \doi{10.1371/journal.pone.0142143}
#'
#' @keywords internal
#' @useDynLib PhysioNetPhysiology, .registration = TRUE
#' @importFrom Rcpp sourceCpp
"_PACKAGE"
