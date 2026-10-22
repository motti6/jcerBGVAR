# Baseline review (2026-10-22 UTC)

Environment: Windows 11 x64, R 4.6.1, Intel i7-1165G7, eight logical CPUs.
Rtools45 exists, but this package needs no compiler. Workspace and its parents
contained no applicable AGENTS.md. No Git executable/repository was found;
no commit ID or remote is fabricated. BGVAR 2.7.0 and roxygen2 8.1.1 were
installed into the excluded project-local .library. testthat 3.3.2, knitr 1.51
and rmarkdown 2.31 were already available. Pandoc 3.8.3 comes from local RStudio.

## BGVAR baseline and boundaries

The [CRAN manual](https://cran.r-project.org/web/packages/BGVAR/BGVAR.pdf)
lists Minnesota, SSVS and Normal-Gamma priors, optional stochastic volatility,
and compiled dependencies including Rcpp, RcppParallel and stochvol. The
estimation API is `bgvar(Data,W,plag,draws,burnin,prior,SV,...,expert)`;
post-processing is `irf(x,n.ahead,shockinfo,quantiles,expert,verbose)`.
Identification instructions are constructed with `get_shockinfo("girf"/"chol"/"sign")`
and `add_shockinfo()`. Local `args()` output is retained in artifacts/logs/bgvar_api.txt.
The installed source confirms separate country estimation and global stacking
before IRF processing. The default compiled IRF path combines rotation searches
and propagation; no independently timed sign-only public call exists. We report
the combined stage and leave sign-only time NA rather than infer a difference.
This release does not assume that MCMC is the dominant expense.

The [official vignette](https://www.stats.bris.ac.uk/R/web/packages/BGVAR/vignettes/examples.html)
describes country VARs augmented with weighted foreign variables, treated as
weakly exogenous, then combined into a global representation. It supports
heterogeneous variable sets and variable-specific weights. In contrast, this
initial implementation requires a common variable set. `data("eerData",package="BGVAR")`
loads both eerData and W.trade0012; the data cover 43 countries and 76 quarters,
1995Q1--2013Q4. `stir` is a short-term market rate. These data remain inside
the optional dependency and are never bundled here. Public expert options
include country parallelism and saved country draws. Signs and generalized
responses are downstream identification choices, not country priors.

The [JSS paper](https://www.jstatsoft.org/article/view/v104i09)
is Boeck, Feldkircher and Huber (2022), 104(9), 1--28,
DOI 10.18637/jss.v104.i09. Direct publisher HTML/PDF retrieval timed out in
this environment; its abstract and bibliographic metadata were checked against
the [author's university record](https://cris.fau.de/publications/328166831/).
The abstract describes shrinkage, stochastic volatility, compiled acceleration
and structural/forecasting tools. Full-paper-specific claims beyond the
available abstract are not made. Actual API details above were checked against
the installed 2.7.0 package and current official manual, not inferred from 2022.

## Fair comparisons

JCERBGVAR uses an MNIW prior and full within-country constant covariance.
BGVAR MN/NG/SSVS and SV choices define different models; comparing their IRFs
does not isolate numerical error. The independent Gibbs reference uses exactly
the same MNIW country likelihood/prior as direct sampling. Only this comparison
supports statements about algorithmic sampling efficiency. World assembly and
IRF recurrence have separate deterministic inverse-formula tests.

The optional benchmark times country estimation using the supported applyfun
hook. A verified local stacking function is temporarily traced for profiling,
then untraced. Only public functions perform estimation/IRFs. No BGVAR source
is incorporated into this package. Compiled sign search and propagation are
timed jointly; memory peaks and standalone sign times remain unmeasured.

## Related research and release policy

[ECB WP2991](https://www.ecb.europa.eu/pub/pdf/scpwps/ecb.wp2991~f5f2cf91dc.en.pdf)
studies approximate variational panel VAR inference. It motivates a later
same-model VI validation; it is not the MNIW algorithm implemented here.
Other verified research references and non-implemented extensions are listed
in CITATIONS.md. No speed claim from a paper is transferred to this project.
The [R extension manual](https://cran.r-project.org/doc/manuals/r-release/R-exts.html)
guides package layout, documentation, vignettes and build/check.
[CRAN policy](https://cran.r-project.org/web/packages/policies.html) and the
[submission checklist](https://cran.r-project.org/web/packages/submission_checklist.html)
require valid maintainer metadata, rights, portability and clean checks.
Author and maintainer placeholders explicitly prevent submission readiness.
