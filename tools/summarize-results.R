Sys.setlocale("LC_CTYPE", ".UTF-8")
.libPaths(c(normalizePath(".library"), .libPaths()))
library(preJCERBGVAR)
out <- "benchmarks/results"
b <- read.csv(file.path(out, "benchmark_results.csv"))
bg <- read.csv(file.path(out, "bgvar_comparison.csv"))
g <- read.csv(file.path(out, "gibbs_comparison.csv"))
gi <- read.csv(file.path(out, "gibbs_irf_comparison.csv"))
a <- read.csv(file.path(out, "accuracy_results.csv"))
v <- read.csv(file.path(out, "irf_algorithm_comparison.csv"))
c <- read.csv(file.path(out, "core_scaling.csv"))
p <- read.csv(file.path(out, "predictive_log_scores.csv"))
med <- aggregate(b[, c("preprocessing", "country_estimation", "world_construction", "irf", "total", "fit_bytes")],
                 b[, c("n_countries", "n_obs", "lags", "cores")], median)
variance <- aggregate(b[, c("country_estimation", "world_construction", "irf", "total")],
                 b[, c("n_countries", "n_obs", "lags", "cores")], var)
write.csv(med, file.path(out, "benchmark_medians.csv"), row.names = FALSE)
write.csv(variance, file.path(out, "benchmark_variances.csv"), row.names = FALSE)
vm <- aggregate(v[, c("dense_reference_seconds", "selected_vector_seconds", "speedup", "max_absolute_error")], v["n_countries"], median)
cm <- aggregate(c[, c("country_estimation", "world_construction", "irf", "total")], c["cores"], median)
write.csv(cm, file.path(out, "core_scaling_medians.csv"), row.names = FALSE)
small <- subset(med, n_obs == 120 & lags == 1)
fmt <- function(z) formatC(z, digits = 4, format = "fg")
tab <- c("| 国数 | 国別推定 秒 | 世界構築 秒 | IRF 秒 | 全体 秒 |",
         "|---:|---:|---:|---:|---:|")
for (i in seq_len(nrow(small))) tab <- c(tab, paste0("| ", paste(c(small$n_countries[i], fmt(unlist(small[i, c("country_estimation", "world_construction", "irf", "total")]))), collapse = " | "), " |"))
summary <- c("# 実測結果 (2026-10-22 UTC)", "",
  "Windows 11 / R 4.6.1 / i7-1165G7、100ドロー、各条件5反復の中央値。すべて合成データ。",
  "計時は他のベンチマークやRチェックと重ねず実施。パッケージロード時間・コンパイルを除外。",
  "3/6/12/20国 × T=80/120/200 × p=1/2 = 24条件。以下はT=120,p=1,cores=1。", "", tab, "",
  paste0("同一MNIW国別モデル、直接2500ドロー対Gibbs2500ドロー(+burn-in1000)の比較: ",
    "Gibbs/直接の計時比 ", fmt(g$speedup_draw_budget), " 倍。直接 ", fmt(g$direct_seconds),
    "秒、Gibbs1チェーン中央値 ", fmt(g$gibbs_seconds_per_chain), "秒。"),
  "4本のGibbsチェーンは同じモデル。実効標本数と独立ドロー数は区別する。",
  paste0("選択係数Rhat=", fmt(g$selected_coefficient_Rhat), ", Gibbs ESS/秒=", fmt(g$gibbs_ESS_per_second),
    ", 直接ESS/秒=", fmt(g$direct_ESS_per_second), "。"),
  paste0("係数事後平均の平均絶対差=", fmt(g$mean_abs_posterior_mean_difference),
    "、分散の平均絶対差=", fmt(g$mean_abs_variance_difference), "。"),
  paste0("日本2変数・h=0..12・25bp IRF: Gibbs対直接の中央値RMSE=", fmt(gi$median_rmse),
    "、90%区間端点平均絶対差=", fmt(gi$interval_endpoint_mae), "。"),
  paste0("条件付き1期先ホールドアウト10点の平均対数予測スコア差(直接-Gibbs)=",
    fmt(mean(p$direct_log_score - p$gibbs_log_score)), "。外国変数の実現経路を所与とする評価。"), "",
  "## 同じ世界ドローに対する決定論的IRFアルゴリズム比較", "",
  "比較参照はG0を一度だけ逆行列化し全インパルス行列を再帰計算。新方式は選択ショックのベクトル再帰。",
  "モデル・事前・ドロー・ショックは同一。計算結果の差は浮動小数点誤差。",
  "| 国数 | 密行列参照 秒 | ベクトル再帰 秒 | 計時比 | 最大絶対差 |",
  "|---:|---:|---:|---:|---:|")
for (i in seq_len(nrow(vm))) summary <- c(summary, paste0("| ", paste(c(vm$n_countries[i], fmt(unlist(vm[i, -1]))), collapse = " | "), " |"))
summary <- c(summary, "", "## 被覆率と解釈", "",
  paste0("20回の既知同時GVARのMonte Carlo: IRF RMSE平均=", fmt(mean(a$rmse)),
    "、68%点別区間被覆率=", fmt(mean(a$coverage68)), "、90%被覆率=", fmt(mean(a$coverage90)),
    "、ゼロでない真のIRFの正符号確率平均=", fmt(mean(a$sign_probability_truth)), "。"),
  "20反復は探索的。条件付きモデルの厳密なサンプリングは、外国変数に内生性があり得る同時DGPでの",
  "名目被覆率を保証しない。識別済み政策因果効果の検証と混同しない。", "",
  "## BGVARとの実務比較", "",
  paste0("BGVAR 2.7.0 MN, SV=FALSE, 100保持ドロー+100burn-in: 推定中央値 ", fmt(median(bg$fit_total)),
    "秒、国別 ", fmt(median(bg$country_estimation)), "秒、世界構築 ", fmt(median(bg$world_stacking)),
    "秒、一般化IRF ", fmt(median(bg$generalized_irf)), "秒、符号識別+IRF ",
    fmt(median(bg$sign_identification_plus_irf)), "秒。"),
  "MNIWとは事前分布・共分散モデルが異なる。単純な時間比やIRF差を同一モデルの高速化・近似誤差と呼ばない。",
  "BGVARの公開APIでは符号探索とIRFが融合。符号探索だけの時間はNAのまま保存。", "",
  "## コア数・未測定事項", "",
  "12国,T=120,p=1の国別ワーカー比較。起動時間を含むため、小さい処理で並列化が遅くなることがある。",
  paste(capture.output(print(cm, row.names = FALSE)), collapse = "\n"), "",
  "ピークRAM未測定。fit_bytesはRオブジェクト容量。コア数の全条件との直積や40国は未実行。",
  "Structured VB/Normal-Gamma MCMC/hybrid/GPU/C++は未実装。ELBOやSVB精度の結果は存在しない。",
  "実データの因果IRFは未実行。BGVARの歴史データは読み込み確認のみ。",
  "PDF/PNG・CSV・生ログはbenchmarks/resultsおよびartifacts/logsに保存。")
writeLines(summary, "docs/RESULTS.md", useBytes = TRUE)
readme <- readLines("README.md", encoding = "UTF-8")
mark <- match("<!-- MEASURED_RESULTS -->", readme)
if (!is.na(mark)) readme <- readme[seq_len(mark - 1L)]
writeLines(c(readme, "", "<!-- MEASURED_RESULTS -->", "## 実測値", "", tab, "",
  paste0("条件: T=120,p=1,100ドロー,5反復中央値。国別MNIWのGibbs/直接計時比は ", fmt(g$speedup_draw_budget),
    " 倍、同モデルの日本IRF中央値RMSEは ", fmt(gi$median_rmse), "。"),
  "全24条件、コア数比較、同じ世界ドローのIRFアルゴリズム比較、BGVARとの仕様差を含む比較は",
  "[実測詳細](docs/RESULTS.md)を参照。", "",
  paste0("20反復Monte Carloの90%点別被覆率は ", fmt(mean(a$coverage90)),
    "。条件付きモデルの仮定と反復数の限界を伴います。")), "README.md", useBytes = TRUE)
# Source-only result plots; every figure carries the synthetic-data label.
for (ext in c("pdf", "png")) {
  path <- file.path(out, paste0("model_comparison.", ext))
  if (ext == "pdf") pdf(path, 8, 5) else png(path, 1000, 650)
  barplot(c(median(bg$fit_total), subset(small, n_countries == 3)$total),
    names.arg = c("BGVAR MN fit", "MNIW fit+IRF+plot"), col = c("orange", "steelblue"),
    ylab = "Median seconds", main = "SYNTHETIC DATA | practical comparison: DIFFERENT priors/models")
  dev.off()
}
cat(summary, sep = "\n")
