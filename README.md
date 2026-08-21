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
# nodes = physiological systems sampled on a common (slow) grid
X <- physioNodeMatrix(list(
  brain = eeg_beta_power,   # e.g. from PhysioEEG
  heart = heart_rate,       # e.g. from PhysioECG
  resp  = respiration,      # e.g. from PhysioECG
  muscle = emg_envelope     # e.g. from PhysioEMG
))

# 1. Time Delay Stability between every pair
tds <- timeDelayStability(X, sampling_rate = 1,
                          window_sec = 60, max_lag_sec = 4)
summary(tds)

# 2. Organ-interaction network, links significance-tested vs surrogates
net <- tdsNetwork(tds, surrogate = TRUE, X = X)
plotTDSnetwork(net)

# 3. How the network reconfigures across states (e.g. sleep stages)
bystate <- tdsNetworkByState(X, states = sleep_stage,
                             window_sec = 60, max_lag_sec = 4)
byState_recon <- bystate$reconfiguration
plotTDSreconfiguration(byState_recon)
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
