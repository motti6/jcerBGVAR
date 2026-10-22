# Implementation status: preJCERBGVAR 0.1.1.9000

Updated: 2026-10-22 UTC. Historical BGVAR validation and revised synthetic validation executed.
Development package. Repository: https://github.com/motti6/jcerBGVAR. Stable release / CRAN submission not performed.
Author/maintainer: motohiro watanabe <watanabemotosan@gmail.com>. Standard GPL-3 adopted by user instruction.
Package namespace: preJCERBGVAR. Historical jcer_* functions and JCERBGVAR S3 classes retained for compatibility.

## Completed in this revision

jcer_irf_bgvar() reuses the unchanged BGVAR fitted posterior and its country-error GIRF definition.
Selected-vector propagation replaces intermediate full-impulse matrix calculations; setup and five quantiles
are included in timing. The MNIW fitter, generalized/recursive/sign IRFs and diagnostics remain available.
BGVAR eerData/W.trade0012: 4/12/20-country NG, 4-country MN with long rate, and 4-country NG with SV.
1000 saved draws + 1000 burn-in, five repetitions per case, warm-up excluded, one core.
All 25 runs compared every retained IRF draw and five quantiles against the public BGVAR IRF API.
All responses, one US.stir shock, horizons 0..12, per-draw +0.25 percentage point normalization.
Maximum historical raw-draw discrepancy: 3.108624e-15.
Stable-draw sensitivity uses exactly the same strictly stable subset for both methods.
Separate MNIW comparison aligns reduced-world generalized shocks and labels changes in priors/SV/covariance.
Historical manifests record transformations, countries, periods and excluded trade mass. Observations are not bundled.

## Synthetic improvements

Primary DGP and fitting omit contemporaneous foreign regressors (G0=I, predetermined lags).
Two persistence regimes (.4/.75), 500 independent datasets each, 250 posterior draws per dataset.
The earlier simultaneous DGP is retained separately for 200 diagnostic datasets, not a correctness certificate.
Report pointwise coverage and MCSE; aggregate MCSE clusters response/horizon cells by independent dataset.
Exclude mechanically zero foreign h0 responses from coverage. Separate coefficient coverage is also recorded.
The independent lag-only validation assembly matches the general fitter for three seeds and all 300 world draws.
Matched synthetic IRF timing includes raw-matrix factorization, impact, normalization, recurrence and quantiles.
Old 3.6x timing had unequal setup boundaries and is superseded; legacy tables remain archived.

## Checks

13 test blocks / 81 passing expectations; no test failures/errors/warnings.
R CMD check --as-cran --no-manual: Status: 1 NOTE
All four vignettes built. README example rerun in a fresh R session against the archive-installed build.
--no-manual omits LaTeX compilation. _R_CHECK_FORCE_SUGGESTS_=false; optional pkgdown unavailable.
Development version is intentionally marked .9000; incoming/CRAN notes are recorded in the check log.
Local Windows 11 / R 4.6.1 / i7-1165G7 checks are recorded separately from remote CI.
Windows/macOS/Linux GitHub Actions workflow: https://github.com/motti6/jcerBGVAR/actions/workflows/R-CMD-check.yaml.
Consult the individual workflow run for current cross-platform outcomes; this file does not infer success from configuration alone.

## Limits

Post-processing speedup preserves the fitted posterior; it is not an MCMC inference speedup.
Public BGVAR baseline also returns auxiliary median structural summaries; timing scope discloses this.
NG/MN/SV fits can have low stable-world draw fractions. Timing keeps a fixed budget and reports these rates.
Single-chain functional ESS does not establish convergence; benchmarks do not certify policy conclusions.
Market-rate generalized innovations are not an identified Fed causal policy shock.
Two matched DGPs do not establish universal calibration. Inner posterior quantiles use 250 draws.
Peak process memory unmeasured; result object sizes are not peak RAM.
Normal-Gamma SVB/ELBO, hybrid, GPU, C++ inference, full simultaneous posterior and heterogeneous variables remain unimplemented.
method='vb'/'hybrid' fails explicitly. No fabricated authors, DOI, organization approval or GitHub identity.
GitHub Desktop's bundled Git is available and the specified development origin is configured.
Original JSS full-text retrieval timed out; current manual/vignette and installed BGVAR 2.7.0 behavior were checked.
Stable-release gates remain in RELEASE_CHECKLIST.md. Historical numerical results were obtained under package name JCERBGVAR 0.1.1.

## Evidence and reproduction

- docs/HISTORICAL_VALIDATION.md: current comparisons, timings, numerical errors and limitations.
- benchmarks/results/historical_20261022/: repetitions, summaries, sensitivities, distribution comparisons, figures.
- benchmarks/results/synthetic_revised_20261022/: 1200-dataset validation, pointwise MCSE, corrected timing.
- artifacts/logs/: execution/test/build/check logs; preJCERBGVAR.Rcheck/00check.log is the clean check report.
- artifacts/preJCERBGVAR_0.1.1.9000.tar.gz and artifacts/SHA256SUMS.txt: build/checksum.

Run tools/reproduce.ps1, or the R scripts listed in docs/HISTORICAL_VALIDATION.md.
Package execution uses base/recommended R packages; optional historical validation requires BGVAR.
