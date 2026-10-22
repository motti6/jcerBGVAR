# References and implementation provenance

All production numerical code is independently written for this project.
No BGVAR code, example observations, figures or tables are redistributed.

- Baseline package and conceptual GVAR comparison: Boeck, M., Feldkircher, M.,
  and Huber, F. (2022), *BGVAR: Bayesian Global Vector Autoregressions with
  Shrinkage Priors in R*, Journal of Statistical Software 104(9), 1--28.
  [DOI](https://doi.org/10.18637/jss.v104.i09). Bibliography/abstract verified;
  publisher full-paper access timed out. Current implementation API checked
  against the [CRAN manual](https://cran.r-project.org/web/packages/BGVAR/BGVAR.pdf)
  and installed package 2.7.0.
- Research background only, IRGA not implemented: Feldkircher, M., Huber, F.,
  Koop, G., and Pfarrhofer, M. (2022), *Approximate Bayesian inference and
  forecasting in huge-dimensional multicountry VARs*, International Economic
  Review. [Publisher](https://onlinelibrary.wiley.com/doi/full/10.1111/iere.12577).
- Research background only, panel VI not implemented: ter Steege, L. (2024),
  *Variational inference for Bayesian panel VAR models*, ECB Working Paper 2991.
  [Primary paper](https://www.ecb.europa.eu/pub/pdf/scpwps/ecb.wp2991~f5f2cf91dc.en.pdf).
- Research background only, shrinkage VI not implemented: Gefang, D., Koop, G.,
  and Poon, A. (2023), *Forecasting using variational Bayesian inference in large
  vector autoregressions with hierarchical shrinkage*, International Journal of
  Forecasting 39(1), 346--363.
  [Publisher](https://www.sciencedirect.com/science/article/abs/pii/S0169207021001862),
  DOI 10.1016/j.ijforecast.2021.11.012.
- Advanced sign/rank search not implemented: Chan, J., Matthes, C., and Yu, X.
  (2026), *Large structural VARs with multiple sign and ranking restrictions*,
  Quantitative Economics 17(3), 709--740.
  [Publisher](https://onlinelibrary.wiley.com/doi/full/10.3982/QE2529).
  This software uses elementary single-column Haar rejection with an explicitly
  conditional rotation target, not that paper's joint accept/reject algorithm.

MNIW algebra, the country-to-world maps and the independent reference Gibbs
conditionals implemented here are fully stated in vignettes/methodology.Rmd.
No external speed ratios are reused. Package author and maintainer: motohiro
watanabe <watanabemotosan@gmail.com>. See inst/CITATION for the development
software citation and RELEASE_CHECKLIST.md for stable-release gates.
