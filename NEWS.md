Version 0.5.2
----------------------------------------------------------------------

- Fixed Rcpp namespace import so Rcpp symbols are properly loaded at runtime
- Updated testthat comparisons against `stats::stl()` to use explicit matched settings and platform-robust tolerances
- Added periodic-mode regression test coverage
- Corrected centering in `plot_seasonal()` so cycle-subseries are centered by the seasonal component mean
- Enabled explicit native routine registration (`useDynLib` with `.registration = TRUE`)
- Removed stale `LazyData` field from DESCRIPTION

Version 0.5
----------------------------------------------------------------------

- Updates to be CRAN-ready (0.5.1)
- Major overhaul of code to match more modern R package development standards (0.5.0)
- Update old C interface to Rcpp (0.5.0)
