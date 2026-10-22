Sys.setlocale("LC_CTYPE", ".UTF-8")
.libPaths(c(normalizePath(".library"), .libPaths()))
hdir <- "benchmarks/results/historical_20261022"
sdir <- "benchmarks/results/synthetic_revised_20261022"
h <- read.csv(file.path(hdir, "same_bgvar_timing_raw.csv"))
m <- read.csv(file.path(hdir, "different_model_timing_raw.csv"))
d <- read.csv(file.path(hdir, "irf_distribution_differences.csv"))
ss <- read.csv(file.path(sdir, "scenario_summary.csv"))
st <- read.csv(file.path(sdir, "matched_irf_timing.csv"))
stopifnot(nrow(h) == 25, all(h$draws == 1000), all(table(h$case) == 5),
          all(ss$independent_replications[ss$case != "simultaneous_diagnostic"] == 500),
          ss$independent_replications[ss$case == "simultaneous_diagnostic"] == 200)
numeric_cols <- c("fit_seconds", "bgvar_irf_seconds", "fast_irf_seconds", "bgvar_total_seconds", "fast_total_seconds", "stable_draw_rate")
hm <- aggregate(h[, numeric_cols], h["case"], median)
hm$irf_speedup <- hm$bgvar_irf_seconds / hm$fast_irf_seconds
hm$total_speedup <- hm$bgvar_total_seconds / hm$fast_total_seconds
hm$max_absolute_error <- vapply(hm$case, function(cs) max(h$max_draw_absolute_error[h$case == cs]), numeric(1))
hm$max_quantile_error <- vapply(hm$case, function(cs) max(h$max_quantile_absolute_error[h$case == cs]), numeric(1))
write.csv(hm, file.path(hdir, "same_bgvar_timing_summary.csv"), row.names = FALSE)
hv <- aggregate(h[, numeric_cols], h["case"], var)
write.csv(hv, file.path(hdir, "same_bgvar_timing_variances.csv"), row.names = FALSE)
mm <- aggregate(m[, c("bgvar_fit_seconds", "mniw_fit_seconds", "mniw_total_seconds", "practical_time_ratio",
  "median_irf_rmse", "interval_endpoint_mae", "median_sign_agreement", "mniw_stable_rate")], m["case"], median)
write.csv(mm, file.path(hdir, "different_model_summary.csv"), row.names = FALSE)
points <- list()
for (cs in hm$case) {
  obj <- readRDS(file.path(hdir, paste0(cs, "_irf_results.rds")))
  aa <- obj$bgvar_world_definition; bb <- obj$mniw_world_definition
  for (rr in aa$response) for (hh in aa$horizon) {
    i <- match(rr, aa$response); j <- match(rr, bb$response); hidx <- hh + 1L
    points[[length(points) + 1L]] <- data.frame(case = cs, response = rr, horizon = hh,
      bgvar_median = aa$quantiles[3, i, hidx], mniw_median = bb$quantiles[3, j, hidx],
      bgvar_q05 = aa$quantiles[1, i, hidx], bgvar_q95 = aa$quantiles[5, i, hidx],
      mniw_q05 = bb$quantiles[1, j, hidx], mniw_q95 = bb$quantiles[5, j, hidx],
      bgvar_positive_probability = mean(aa$draws[i, hidx, ] > 0),
      mniw_positive_probability = mean(bb$draws[j, hidx, ] > 0),
      note = "same world-generalized shock, different priors/covariance model; final repetition")
  }
}
write.csv(do.call(rbind, points), file.path(hdir, "different_model_pointwise.csv"), row.names = FALSE)
stm <- aggregate(st[, c("dense_seconds", "vector_seconds", "speedup", "max_absolute_error")], st["countries"], median)
fmt <- function(x) formatC(x, format = "fg", digits = 4)
table_lines <- c("| 条件 | BGVAR IRF 秒 | 選択ベクトルIRF 秒 | IRF計時比 | 従来の全体 秒 | 置換後の全体 秒 | 全体計時比 |",
  "|---|---:|---:|---:|---:|---:|---:|")
for (cs in c("NG_4", "NG_12", "NG_20", "MN_4_long", "NG_4_SV")) {
  z <- hm[hm$case == cs, ]
  table_lines <- c(table_lines, paste0("| ", cs, " | ", paste(fmt(unlist(z[, c("bgvar_irf_seconds", "fast_irf_seconds", "irf_speedup", "bgvar_total_seconds", "fast_total_seconds", "total_speedup")])), collapse = " | "), " |"))
}
report <- c("# BGVAR付属データによる高速化評価 (0.1.1 / 2026-10-22 UTC)", "",
  "主比較は、既存BGVARの推定をそのまま使い、IRF計算のみを選択ベクトル再帰へ置き換える比較。",
  "MNIWへの推定法変更は別表に分け、同一事後分布の高速化と混同しない。", "",
  "## データと共通条件", "",
  "インストール済みBGVAR 2.7.0のeerDataとW.trade0012を実行時に読み込み。観測値の複製配布はしない。",
  "元標本1995Q1〜2013Q4の76期。yは100*diff(log GDP)、Dpは100倍、短期/長期金利は",
  "入力単位をdecimalと明示して100倍しpercent pointsへ。変換後75期、p=1、切片あり・トレンドなし。",
  "国別データは共通標本と共通変数へ揃える。対象国を除外した貿易ウェイトは行ごとに明示的に再正規化。",
  "4国はUS/JP/EA/KR。12/20国は同じ先頭4国に、付属データ順で必要変数の揃う国を追加する。",
  "GDP・物価・短期金利の3変数が主比較。MN_4_longだけ長期金利を加える。各国・期間・残存貿易比率は*_manifest.csv。",
  "地域集計EAと一部構成国を同時に含む拡張例は計算規模の検証であり、国際政策モデルの推奨仕様ではない。", "",
  "stirは短期市場金利で、政策決定系列そのものではない。金利単位は分析設定として明示し、値から自動判定しない。", "",
  "## 同じBGVAR事後分布での速度と数値的一致", "",
  "NGはNormal-Gamma、MNはMinnesota。主比較SV=FALSE、NG_4_SVはSV=TRUE。各1000保持ドロー、",
  "1000 burn-in、thin=1、1コア、ウォームアップ除外、独立seedで5反復。表は各所要時間の中央値の比。",
  "US.stirの初期応答をドローごとに+0.25pp(25bp)へ規格化。1ショック・全応答・h=0..12・",
  "全IRFドローと同じ5分位点を計算。推定は共有し、全体時間=同じ推定計時+それぞれのIRF計時。", "", table_lines, "",
  paste0("全25回の最大ドロー絶対差=", format(max(h$max_draw_absolute_error), scientific = TRUE),
    "、最大分位点絶対差=", format(max(h$max_quantile_absolute_error), scientific = TRUE), "。"),
  "この比較では事前分布・MCMC・共分散ドロー・国際ウェイト・安定性選択を変更しない。",
  "同じ既存BGVAR事後分布を対象としているので、ここでの差は数値演算誤差として評価できる。",
  "IRF計時には入力抽出、初期ショック計算、規格化、再帰、分位点計算を含む。実行順は交互にする。",
  "BGVAR公開IRFは補助的な構造中央値等も返すため、その付随コストも含む公開API間の比較である。",
  "低水準の再帰ループだけの速度比や、MCMC自体の高速化を表す値ではない。", "",
  "## 安定性・精度上の限界", "",
  "主計時では同数の保持ドローを処理するためeigen=FALSEとし、非定常ドローも保持した。",
  "厳密に安定なドローの比率は次の通り。これが低い条件のIRFを政策分析として解釈してはならない。",
  paste(vapply(seq_len(nrow(hm)), function(i) paste0("- ", hm$case[i], ": ", fmt(100 * hm$stable_draw_rate[i]), "%"), character(1)), collapse = "\n"),
  "各条件の初回では、両手法を同じ安定ドロー部分集合に限定してIRF一致も検査した(*_stable_sensitivity.csv)。",
  "5反復の計時だけではMCMC収束を保証しない。selected_draw_diagnostics.csvのESSは選択されたIRF関数の",
  "単一チェーン記述統計であり、包括的な収束診断ではない。ピークRAMは未測定で、結果オブジェクト容量と区別。", "",
  "## 従来手法との違い", "",
  "| 項目 | BGVAR公開IRF → jcer_irf_bgvar | BGVAR推定 → MNIW推定 |",
  "|---|---|---|",
  "| 事前分布・事後ドロー | 同一 | NG/MNとMNIWで異なる |",
  "| 分散モデル | 同じ保存済み共分散、SVケースも維持 | MNIWは定数分散。SVモデルとは異なる |",
  "| ショック定義 | BGVARのcountry-error GIRFを維持 | 世界縮約形一般化ショックへ揃えて比較 |",
  "| 計算の変更 | 全インパルス行列の中間計算を避け、1ショックのベクトルを伝播 | MCMCから共役直接ドローへ変更 |",
  "| 分布差の解釈 | 浮動小数点誤差 | 事前・分散モデルの違いを含む |", "",
  "BGVAR 2.7.0のGIRFはG0^{-1} Sigma_e e_j、世界縮約形GIRFはSigma_u e_jに比例する。",
  "両者は一般に同じではない。definition='bgvar'を主速度比較に固定する。これはローカル実装を",
  "実際に確認した上でのショック定義の区別で、既存BGVARの誤りを意味しない。",
  "同じ事後ドローでショック定義だけを変えた差もirf_distribution_differences.csvへ保存。", "",
  "### MNIWへの変更を含む別モデル比較", "",
  "世界縮約形の一般化ショックと25bpを共通化。下表のIRF差を計算誤差と呼ばない。",
  "MNIWのprior_controlは既定値であり、BGVARのNG/MN階層事前分布には一致させていない。",
  "| 条件 | BGVAR推定 秒 | MNIW推定 秒 | 世界IRF中央値RMSE | 区間端点平均絶対差 |",
  "|---|---:|---:|---:|---:|")
for (i in seq_len(nrow(mm))) report <- c(report, paste0("| ", mm$case[i], " | ", paste(fmt(unlist(mm[i, c("bgvar_fit_seconds", "mniw_fit_seconds", "median_irf_rmse", "interval_endpoint_mae")])), collapse = " | "), " |"))
report <- c(report, "", "## 合成データ検証の改善", "",
  "主検証は同時点外国係数B_i0=0・推定も同時点外国変数なしとし、尤度とDGPを整合させた。",
  "持続性0.4と0.75を各500独立データセットで検証。従来の同時点DGPは200回の仮定違反診断として分ける。",
  "各データセット250事後ドロー、3国2変数,T=120。機械的にゼロとなる海外h=0は被覆率に含めない。",
  "被覆率とMonte Carlo標準誤差は点別に保存し、全体集計では同一データ内の応答・期間の相関を考慮する。",
  "| シナリオ | 独立反復数 | 68%被覆率 | MC標準誤差 | 90%被覆率 | MC標準誤差 |",
  "|---|---:|---:|---:|---:|---:|")
for (i in seq_len(nrow(ss))) report <- c(report, paste0("| ", ss$case[i], " | ", paste(c(ss$independent_replications[i], fmt(unlist(ss[i, c("coverage68", "coverage68_mcse", "coverage90", "coverage90_mcse")]))), collapse = " | "), " |"))
report <- c(report, "", "これらは選択したDGP・既定事前分布・250事後ドローに対する結果。名目被覆率の一般保証ではない。",
  "従来の20反復・同時点DGPの集約被覆率を主な精度根拠として使わない。",
  "旧来の約3.6倍IRF比較も、分解と事前準備の計時範囲が異なっていたため主張を差し替える。",
  "修正後は両経路とも生のG0/G_l/Sigma_eから始め、分解・ショック・規格化・再帰・分位点を含む。",
  "| 国数 | 密行列参照 秒 | ベクトル再帰 秒 | 計時比(各反復の中央値) | 最大絶対差 |",
  "|---:|---:|---:|---:|---:|")
for (i in seq_len(nrow(stm))) report <- c(report, paste0("| ", paste(c(stm$countries[i], fmt(unlist(stm[i, -1]))), collapse = " | "), " |"))
report <- c(report, "", "## 再実行・成果物", "", "```powershell",
  "Rscript benchmarks/run_historical_validation.R", "Rscript tools/check-synthetic-engine.R",
  "Rscript benchmarks/run_synthetic_validation.R", "Rscript benchmarks/run_irf_benchmark_matched.R",
  "Rscript tools/summarize-historical-validation.R", "```", "",
  "各スクリプトはJCER_BGVAR直下で実行。タイミング用ジョブは重ねず実行する。",
  "CSV・IRF分布RDS・PDF/PNG: benchmarks/results/historical_20261022/。",
  "合成検証: benchmarks/results/synthetic_revised_20261022/。生ログ: artifacts/logs/。",
  "[BGVAR公式マニュアル](https://cran.r-project.org/web/packages/BGVAR/BGVAR.pdf)と",
  "[公式vignette](https://www.stats.bris.ac.uk/R/web/packages/BGVAR/vignettes/examples.html)でデータ/APIを確認した。",
  "一般化ショックへの応答であり、識別済み米国金融政策の因果効果という主張は行わない。")
if (!file.exists("docs/RESULTS_legacy_0.1.0.md")) file.copy("docs/RESULTS.md", "docs/RESULTS_legacy_0.1.0.md")
writeLines(report, "docs/HISTORICAL_VALIDATION.md", useBytes = TRUE)
writeLines(report, "docs/RESULTS.md", useBytes = TRUE)
readme <- readLines("README.md", encoding = "UTF-8")
mark <- match("<!-- MEASURED_RESULTS -->", readme)
if (!is.na(mark)) readme <- readme[seq_len(mark - 1L)]
writeLines(c(readme, "", "<!-- MEASURED_RESULTS -->", "## 付属データでの実測", "",
  "同じBGVAR事後分布・ショック定義・1000保持ドロー、各5反復の中央値。全体時間は同じ推定時間を含みます。", "",
  table_lines, "", "事前分布をMNIWへ変更した比較、安定性、合成DGPと計時範囲の修正は",
  "[付属データ検証報告](docs/HISTORICAL_VALIDATION.md)を参照。旧3.6倍という合成比較は主張を差し替えています。"), "README.md", useBytes = TRUE)

for (ext in c("pdf", "png")) {
  path <- file.path(hdir, paste0("timing_comparison.", ext))
  if (ext == "pdf") pdf(path, 10, 6) else png(path, 1400, 850)
  par(mfrow = c(1, 2), mar = c(7, 4, 4, 1))
  barplot(t(as.matrix(hm[, c("bgvar_irf_seconds", "fast_irf_seconds")])), beside = TRUE, names.arg = hm$case, las = 2,
    col = c("grey50", "steelblue"), ylab = "Median seconds", main = "HISTORICAL eerData | unchanged BGVAR draws\nIRF including setup and five quantiles")
  legend("topleft", c("BGVAR public IRF", "Selected vector"), fill = c("grey50", "steelblue"), cex = .7)
  barplot(t(as.matrix(hm[, c("bgvar_total_seconds", "fast_total_seconds")])), beside = TRUE, names.arg = hm$case, las = 2,
    col = c("grey50", "steelblue"), ylab = "Median seconds", main = "HISTORICAL eerData | shared inference\nFit plus IRF")
  dev.off()
  for (cs in c("NG_4", "MN_4_long")) {
    obj <- readRDS(file.path(hdir, paste0(cs, "_irf_results.rds")))
    a <- obj$bgvar_world_definition; b <- obj$mniw_world_definition
    resp <- intersect(c("JP.y", "JP.Dp", "JP.stir", "JP.ltir"), a$response)
    path <- file.path(hdir, paste0(cs, "_model_irf_comparison.", ext))
    if (ext == "pdf") pdf(path, 9, 7) else png(path, 1200, 950)
    par(mfrow = c(2, 2), mar = c(3, 4, 4, 1))
    for (rr in resp) {
      i <- match(rr, a$response); j <- match(rr, b$response)
      qa <- matrix(a$quantiles[, i, ], 5); qb <- matrix(b$quantiles[, j, ], 5)
      plot(a$horizon, qa[3, ], type = "n", ylim = range(qa[c(1, 5), ], qb[c(1, 5), ]),
        xlab = "Horizon", ylab = rr, main = paste("HISTORICAL eerData |", rr))
      polygon(c(a$horizon, rev(a$horizon)), c(qa[1, ], rev(qa[5, ])), col = adjustcolor("grey50", .25), border = NA)
      polygon(c(b$horizon, rev(b$horizon)), c(qb[1, ], rev(qb[5, ])), col = adjustcolor("steelblue", .2), border = NA)
      lines(a$horizon, qa[3, ], col = "grey30", lwd = 2); lines(b$horizon, qb[3, ], col = "steelblue", lwd = 2)
      abline(h = 0, lty = 3)
      mtext("World generalized +25bp | DIFFERENT priors/models | last repetition", side = 3, cex = .55)
      legend("topright", c("BGVAR", "MNIW"), col = c("grey30", "steelblue"), lty = 1, cex = .7)
    }
    dev.off()
  }
}
cat(paste(table_lines, collapse = "\n"), "\n")
