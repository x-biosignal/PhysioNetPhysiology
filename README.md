# PhysioNetPhysiology

**Network Physiology for the Physio ecosystem** — dynamic-interaction analysis
of coupled physiological systems (brain, heart, respiration, muscle).

Single-modality tools cannot ask how *different* organ systems couple and
reorganise together. This package does, on top of the multimodal
`PhysioExperiment` data model, via **Time Delay Stability (TDS)**
(Bashan et al. 2012): it tracks, over sliding windows, the time lag that
maximises each signal pair's cross-correlation, and rewards periods where that
lag stays put — the fingerprint of a genuine dynamic coupling. Stable couplings
become links in an *organ-interaction network*, and the package measures how
that network **reconfigures across physiological states**.

```r
library(PhysioNetPhysiology)

# Four physiological systems on a common (slow) grid. Here a shared driver
# couples brain/heart/resp at fixed lags, while muscle stays uncoupled; in
# practice these come from PhysioEEG/PhysioECG/PhysioEMG.
set.seed(42)
n <- 1200
drive <- cumsum(rnorm(n)); drive <- (drive - mean(drive)) / sd(drive)
lag <- function(v, k) c(rep(v[1], k), head(v, -k))
X <- physioNodeMatrix(list(
  brain  = drive          + rnorm(n, sd = 0.05),
  heart  = lag(drive, 2)  + rnorm(n, sd = 0.05),
  resp   = lag(drive, 4)  + rnorm(n, sd = 0.05),
  muscle = rnorm(n)
))

# 1. Time Delay Stability between every pair
tds <- timeDelayStability(X, sampling_rate = 1,
                          window_sec = 60, max_lag_sec = 8)
summary(tds)

# 2. Organ-interaction network, links significance-tested vs surrogates
net <- tdsNetwork(tds, surrogate = TRUE, X = X)
plotTDSnetwork(net)

# 3. How the network reconfigures across states (e.g. sleep stages)
sleep_stage <- rep(c("wake", "sleep"), each = n / 2)
bystate <- tdsNetworkByState(X, states = sleep_stage, sampling_rate = 1,
                             window_sec = 60, max_lag_sec = 8)
plotTDSreconfiguration(bystate)
```

## Installation

The ecosystem builds on Bioconductor, so its repositories have to be on the
list as well -- without them the install stops at `SummarizedExperiment`.

```r
install.packages("BiocManager", repos = "https://cloud.r-project.org")
install.packages(
  "PhysioNetPhysiology",
  repos = c("https://x-biosignal.r-universe.dev", BiocManager::repositories())
)
```

From GitHub instead:

```r
install.packages("remotes", repos = "https://cloud.r-project.org")
remotes::install_github("x-biosignal/PhysioNetPhysiology")
```

## Why TDS

TDS is robust to the amplitude non-stationarity of physiological signals and to
differing signal types (it operates on the *timing* of coupling, not its
strength), which is why it underlies the Network-Physiology programme. Links are
significance-tested against **phase-randomised surrogates** that preserve each
signal's own spectrum while destroying cross-signal timing.

## Scope

The core engine works on any multivariate time series; the `physioNodeMatrix()`
helper and the `PhysioExperiment` integration make organ-system assembly
convenient. Reference: Bashan et al. (2012) *Nat Commun* 3:702; Bartsch et al.
(2015) *PLoS ONE* 10:e0142143.
