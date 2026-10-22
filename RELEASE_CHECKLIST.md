# preJCERBGVAR development publication

- [x] R package name: preJCERBGVAR; development version: 0.1.1.9000.
- [x] Author/maintainer: motohiro watanabe <watanabemotosan@gmail.com>.
- [x] Development repository: https://github.com/motti6/jcerBGVAR.
- [x] Standard GPL-3 chosen under the user's instruction; full text in COPYING.
- [x] DESCRIPTION, R citation and CITATION.cff contain the confirmed metadata.
- [x] jcer_* functions and old S3 class names retained for saved-object compatibility.
- [x] BGVAR historical-data and revised synthetic validation scripts/results provided.
- [x] Git excludes local libraries, archives, working logs and binary posterior draws.
- [x] No official institutional endorsement or stable-release status is asserted.
- [ ] Push the reviewed development source to the specified GitHub repository.
- [ ] Inspect actual Windows/macOS/Linux CI results after pushing; local Windows results are recorded separately.

## Before a stable release or CRAN submission

- [ ] Verify supported BGVAR versions and add changes to compatibility tests.
- [ ] Review empirical stability, MCMC convergence, prior sensitivity and intended use cases.
- [ ] Complete independent installation/code review and any additional organization-specific release process.
- [ ] Recheck naming availability, documentation, citation version and release notes.
- [ ] Run submission checks with the intended dependency/toolchain environment and explain remaining notes.
- [ ] Deliberately choose a stable version and tag/release it. No CRAN submission is authorized by this development setup.

Unimplemented SVB/hybrid/GPU methods remain future work and are not advertised
as available. See CONTRIBUTING.md for branches, versioning and API changes.
