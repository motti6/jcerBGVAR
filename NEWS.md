# preJCERBGVAR 0.1.1.9000 (development)

* Rename the R package to preJCERBGVAR; retain jcer_* functions and historical S3 classes for compatibility.
* Confirm motohiro watanabe as author and maintainer, with the supplied public contact.
* Set motti6/jcerBGVAR as the development repository and use standard GPL-3.
* Add development contribution/versioning guidance and exclude local libraries, logs and raw posterior binaries from Git.

# JCERBGVAR 0.1.1 (previous development name)

* Add selected-vector IRFs reusing an unchanged BGVAR posterior and matching its country-error GIRF.
* Distinguish BGVAR country-error shocks from reduced-world generalized innovations.
* Add lag-only synthetic DGPs consistent with conditional exogeneity and configurable persistence.
* Add historical eerData validation and timing, with setup/quantiles included in both IRF paths.

# JCERBGVAR 0.1.0

* Independent conditional MNIW country posterior sampling and linked world systems.
* Known stable synthetic GVAR, generalized/recursive and single-shock sign IRFs.
* Explicit posterior target, interest-rate units, stability filtering diagnostics.
* Reproducible country workers, numerical tests and independent Gibbs validation.
* VB, hybrid, GPU and full simultaneous likelihood are research roadmap items.
