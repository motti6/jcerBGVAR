param([string]$RBin = 'C:\Program Files\R\R-4.6.1\bin')
$ErrorActionPreference = 'Stop'
$env:LC_ALL = 'Japanese_Japan.utf8'
$env:TZ = 'UTC'
$env:OMP_NUM_THREADS = '1'
$env:OPENBLAS_NUM_THREADS = '1'
Set-Location -LiteralPath (Split-Path -Parent $PSScriptRoot)
$rscript = Join-Path $RBin 'Rscript.exe'
$rexe = Join-Path $RBin 'R.exe'
# Capture process streams via Start-Process so PowerShell does not reinterpret
# R's normal stderr progress output as a failed command. No visible window.
function Invoke-RTask([string]$Name, [string[]]$RArguments, [string]$Executable = $rscript) {
  $p = Start-Process -FilePath $Executable -ArgumentList $RArguments -Wait -PassThru -WindowStyle Hidden `
    -RedirectStandardOutput "artifacts/logs/$Name.stdout.log" -RedirectStandardError "artifacts/logs/$Name.stderr.log"
  if ($p.ExitCode -ne 0) { throw "$Name failed with exit code $($p.ExitCode)" }
}
New-Item -ItemType Directory -Path 'artifacts/logs' -Force | Out-Null
Invoke-RTask 'bootstrap' @('--vanilla', 'tools/bootstrap.R')
Invoke-RTask 'unit-tests' @('--vanilla', 'tools/document-test.R')
Invoke-RTask 'install' @('--vanilla', 'tools/install-local.R')
Invoke-RTask 'parallel' @('--vanilla', 'tools/parallel-smoke.R')
# Benchmark jobs are sequential to avoid CPU contention between measurements.
Invoke-RTask 'historical-validation' @('--vanilla', 'benchmarks/run_historical_validation.R')
Invoke-RTask 'synthetic-engine-check' @('--vanilla', 'tools/check-synthetic-engine.R')
Invoke-RTask 'synthetic-revision-validation' @('--vanilla', 'benchmarks/run_synthetic_validation.R')
Invoke-RTask 'irf-matched-benchmark' @('--vanilla', 'benchmarks/run_irf_benchmark_matched.R')
Invoke-RTask 'historical-summary' @('--vanilla', 'tools/summarize-historical-validation.R')
$env:R_LIBS = (Resolve-Path '.library').Path
$pandocDir = Join-Path $env:LOCALAPPDATA 'Programs\RStudio\resources\app\bin\quarto\bin\tools'
if (Test-Path -LiteralPath (Join-Path $pandocDir 'pandoc.exe')) {
  $env:RSTUDIO_PANDOC = $pandocDir
  $env:PATH = "$pandocDir;$env:PATH"
}
$env:_R_CHECK_FORCE_SUGGESTS_ = 'false'
Invoke-RTask 'build' @('CMD', 'build', '.') $rexe
$packageVersion = ((Get-Content -LiteralPath DESCRIPTION | Where-Object { $_ -match '^Version:' }) -replace '^Version:\s*', '').Trim()
$archiveName = "preJCERBGVAR_$packageVersion.tar.gz"
Invoke-RTask 'check' @('CMD', 'check', '--as-cran', '--no-manual', $archiveName) $rexe
Copy-Item -LiteralPath $archiveName -Destination "artifacts/$archiveName" -Force
Invoke-RTask 'install-archive' @('CMD', 'INSTALL', '--library=.library', "artifacts/$archiveName") $rexe
Invoke-RTask 'readme-example' @('--vanilla', 'tools/readme-example.R')
Invoke-RTask 'final-status' @('--vanilla', 'tools/finalize-validation.R')
Invoke-RTask 'archive-audit' @('--vanilla', 'tools/archive-audit.R')
Get-FileHash -LiteralPath "artifacts/$archiveName" -Algorithm SHA256 |
  ForEach-Object { "$($_.Hash)  $archiveName" } |
  Set-Content -LiteralPath 'artifacts/SHA256SUMS.txt' -Encoding ASCII
Write-Output 'Completed. Review docs/RESULTS.md and preJCERBGVAR.Rcheck/00check.log.'
