#include <Rcpp.h>
#include <cmath>
#include <vector>
using namespace Rcpp;

// Windowed peak-lag cross-correlation for one signal pair -- the TDS engine's
// inner loop. Mirrors the pure-R .pair_tds() exactly: per window, z-normalise
// (sample sd, n-1) both segments, then find the lag maximising |cross-corr|
// (cc = sum(a*b)/win). long double accumulation matches R's sum() precision.
//
// [[Rcpp::export]]
List tds_peaklag_cpp(NumericVector x, NumericVector y,
                     IntegerVector starts, int win, int max_lag) {
  int nt = starts.size();
  IntegerVector tau0(nt);
  NumericVector peak(nt);
  std::vector<double> xw(win), yw(win);

  for (int t = 0; t < nt; ++t) {
    int s0 = starts[t] - 1;                      // R 1-based -> 0-based

    long double mx = 0.0L, my = 0.0L;
    for (int i = 0; i < win; ++i) { mx += x[s0 + i]; my += y[s0 + i]; }
    mx /= win; my /= win;

    long double vx = 0.0L, vy = 0.0L;
    for (int i = 0; i < win; ++i) {
      long double dx = x[s0 + i] - mx, dy = y[s0 + i] - my;
      vx += dx * dx; vy += dy * dy;
    }
    double sx = (double) std::sqrt(vx / (win - 1));
    double sy = (double) std::sqrt(vy / (win - 1));
    for (int i = 0; i < win; ++i) {
      xw[i] = (sx > 0.0) ? (double)((x[s0 + i] - mx) / sx) : 0.0;
      yw[i] = (sy > 0.0) ? (double)((y[s0 + i] - my) / sy) : 0.0;
    }

    double best = -1.0; int blag = 0; double bcc = 0.0;
    for (int lag = -max_lag; lag <= max_lag; ++lag) {
      long double cc = 0.0L;
      if (lag >= 0) {
        for (int i = 0; i < win - lag; ++i) cc += (long double) xw[i] * yw[i + lag];
      } else {
        for (int i = 0; i < win + lag; ++i) cc += (long double) xw[i - lag] * yw[i];
      }
      double ccd = (double)(cc / win);
      double acc = std::fabs(ccd);
      if (acc > best) { best = acc; blag = lag; bcc = ccd; }
    }
    tau0[t] = blag;
    peak[t] = bcc;
  }
  return List::create(Named("tau0") = tau0, Named("peak") = peak);
}
