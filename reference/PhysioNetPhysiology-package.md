# PhysioNetPhysiology: Network Physiology for the Physio ecosystem

Dynamic-interaction analysis of coupled physiological systems. The
package centres on **Time Delay Stability (TDS)** — a method that
detects stable temporal couplings between signals by tracking, over
sliding windows, the lag that maximises their cross-correlation and
rewarding periods where that lag stays put. Stable couplings become
links in an *organ-interaction network* whose edges are
significance-tested against phase-randomised surrogates, and the package
quantifies how such networks **reconfigure across physiological states**
(wake / sleep stages, rest / task, ...).

## Details

Core entry points:

- [`timeDelayStability()`](https://x-biosignal.github.io/PhysioNetPhysiology/reference/timeDelayStability.md)
  — the TDS engine over a multivariate series.

- [`tdsSurrogateThreshold()`](https://x-biosignal.github.io/PhysioNetPhysiology/reference/tdsSurrogateThreshold.md)
  — surrogate-based significance threshold.

- [`tdsNetwork()`](https://x-biosignal.github.io/PhysioNetPhysiology/reference/tdsNetwork.md)
  — build a thresholded organ-interaction network.

- [`tdsNetworkByState()`](https://x-biosignal.github.io/PhysioNetPhysiology/reference/tdsNetworkByState.md)
  — state-resolved networks + reconfiguration.

- [`physioNodeMatrix()`](https://x-biosignal.github.io/PhysioNetPhysiology/reference/physioNodeMatrix.md)
  — assemble a node matrix of physiological systems.

## References

Bashan A, Bartsch RP, Kantelhardt JW, Havlin S, Ivanov PC (2012).
Network physiology reveals relations between network topology and
physiological function. *Nature Communications* 3:702.
[doi:10.1038/ncomms1705](https://doi.org/10.1038/ncomms1705)

Bartsch RP, Liu KKL, Bashan A, Ivanov PC (2015). Network physiology: how
organ systems dynamically interact. *PLoS ONE* 10(11):e0142143.
[doi:10.1371/journal.pone.0142143](https://doi.org/10.1371/journal.pone.0142143)

## See also

Useful links:

- <https://github.com/x-biosignal/PhysioNetPhysiology>

- <https://x-biosignal.r-universe.dev/PhysioNetPhysiology>

- <https://x-biosignal.github.io/PhysioNetPhysiology>

- Report bugs at
  <https://github.com/x-biosignal/PhysioNetPhysiology/issues>

## Author

**Maintainer**: Yusuke Matsui <mail.to.matsui@gmail.com>
