# Developing preJCERBGVAR

This repository hosts a development package, not a stable release or a CRAN
submission. The repository name is `motti6/jcerBGVAR`; the R package name is
`preJCERBGVAR`. The local workspace may still be named `JCER_BGVAR`.

Author and maintainer: motohiro watanabe <watanabemotosan@gmail.com>.
License: standard GNU GPL version 3 (COPYING).

## Versioning and changes

The current version is `0.1.1.9000`. The `.9000` suffix marks ongoing development.
Use short feature branches and pull requests for subsequent changes. Keep the
main branch installable; a stable release/tag is a separate deliberate action.
Change package metadata in DESCRIPTION and update CITATION.cff when releasing.
The R citation reads the installed version automatically.

Functions retain their `jcer_*` names. S3 classes (`JCERBGVAR`,
`JCERBGVAR_irf`, etc.) intentionally retain their historical names so saved
validation objects remain readable after the package rename. Record breaking
API or class changes in NEWS.md with a migration example.

## Validation

From the repository root:

```r
source("tools/bootstrap.R")
source("tools/document-test.R")
source("tools/install-local.R")
```

The light tests run during package checking. Heavy benchmarks remain separate;
`tools/reproduce.ps1` runs the full local historical/synthetic workflow.
Record the seed, package versions, model/prior, shock definition, retained draw
count and timing boundaries for new performance claims. Keep model changes
separate from numerical optimizations of an unchanged posterior.

Preserve historical CSV results rather than replacing them with new values
without a labelled rerun. Update reports and NEWS when assumptions change.
Raw draws, local libraries and working logs are intentionally ignored by Git;
the benchmark scripts regenerate them. Package source tarballs are built locally.

## Stable release gate

Check supported BGVAR versions, run Windows/macOS/Linux CI, review empirical
stability/convergence limits, and update RELEASE_CHECKLIST.md before tagging a
stable release. Do not advertise unimplemented VB/hybrid methods as available.
