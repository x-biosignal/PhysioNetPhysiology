# PhysioNetPhysiology 0.3.0

Classical graph-topology metrics and community detection for the physiological
networks from `tdsNetwork()` (previously the package had only degree/strength and
partition-consuming role metrics; `igraph` was suggested but unused).

- `networkMetrics()` computes the standard node-level (degree, strength, weighted
  clustering, betweenness, closeness, eigenvector centrality, local efficiency)
  and global (density, transitivity, characteristic path length, global/local
  efficiency, degree assortativity, modularity) measures. Edge weights are
  coupling strengths; path-based metrics use `1/weight` distances.
- `communityDetection()` partitions the network (Louvain / fast-greedy / walktrap
  / leading-eigenvector / infomap) and reports modularity Q. Its membership feeds
  the existing `moduleRoles()` / `participationCoefficient()`, closing the loop
  those role metrics needed (they required an externally supplied partition).
- `smallWorldness()` reports the Humphries-Gurney sigma and Telford omega against
  degree-preserving random and lattice references.
- These delegate the graph algorithms to `igraph` (an existing Suggests). Verified
  on known structure: two-block community recovery (Q ~ 0.49), and a
  Watts-Strogatz graph flagged small-world (sigma > 1) while a random graph is not.

# PhysioNetPhysiology 0.2.0

- The Time Delay Stability engine's windowed peak-lag cross-correlation is now
  implemented in C++ (Rcpp): about 3x faster than the pure-R loop and
  byte-identical to it (verified in the test suite). This compounds in the
  surrogate-threshold path, which runs the engine many times.

# PhysioNetPhysiology 0.1.0

* Initial release: Network Physiology for the Physio ecosystem.
* `timeDelayStability()` — the core Time Delay Stability (TDS) engine over a
  multivariate physiological series (Bashan et al. 2012).
* `tdsSurrogateThreshold()` — phase-randomised surrogate significance threshold
  for TDS links.
* `tdsNetwork()` — thresholded, weighted organ-interaction networks.
* `tdsNetworkByState()` / `tdsReconfiguration()` — state-resolved networks and
  their reconfiguration across physiological states (e.g. sleep stages).
* `physioNodeMatrix()` — assemble a node matrix of physiological systems.
* Visualisation: `plotDelayStability()`, `plotTDSnetwork()`,
  `plotTDSreconfiguration()`.
