# preJCERBGVAR — development version 0.1.1.9000

条件付き Bayesian Global VAR の独立したR実装です。正式なRパッケージ名は
`preJCERBGVAR`。共役MNIW事前分布からバーンイン不要の独立事後ドローを生成し、
世界モデル、一般化IRF、明示的な順序に基づく再帰IRF、単一ショックの符号制約を提供します。
既存の `BGVAR` 推定結果には `jcer_irf_bgvar()` を使い、推定法・事前分布・事後ドローを変えずにIRF計算を置き換えられます。
組織の公式承認や安定版のリリースを意味しません。

開発リポジトリ: https://github.com/motti6/jcerBGVAR
安定版のリリース前で、APIや推定手法は変更される可能性があります。変更方針は[CONTRIBUTING.md](CONTRIBUTING.md)に記載しています。

GitHubの開発版は次のコマンドで導入できます。

```r
install.packages("remotes")
remotes::install_github("motti6/jcerBGVAR")
library(preJCERBGVAR)
```

プロジェクトフォルダで、初回は `Rscript tools/install-local.R` を実行してください。
依存開発ツールの任意導入は `Rscript tools/bootstrap.R`。本体の実行依存はR標準パッケージのみです。
以降のRセッションでは `.libPaths(c(normalizePath(".library"), .libPaths()))` を最初に実行します。
ビルドしたソースは `install.packages("artifacts/preJCERBGVAR_0.1.1.9000.tar.gz", repos=NULL, type="source")`
でも導入できます。C++コンパイラは不要です。

```r
library(preJCERBGVAR)
sim <- jcer_simulate_gvar(n_countries = 3, n_obs = 120, seed = 42, contemporaneous = FALSE)
fit <- jcer_bgvar(sim$data, sim$W, method = "conjugate", prior = "mniw",
                  lags = 1, draws = 1000, seed = 42, contemporaneous = FALSE)
print(summary(fit))
ir <- jcer_irf(fit, impulse = "US.r_short",
               response = c("JP.r_short", "JP.y"),
               identification = "generalized", horizon = 12)
plot(ir)
```

すべて合成データです。25bpへの規格化は `jcer_irf(..., shock_bp=25)`。
金利単位 `rate_unit="percent"` なら0.25、`"decimal"`なら0.0025です。
単位はデータ値から推測せず、推定時に明示します。
一般化IRFを純粋なFed政策ショックの因果効果とは解釈しません。

| 機能 | 状態 |
|---|---|
| 国別MNIW直接事後ドロー・世界モデル | 実装、数値テストあり |
| 既存BGVARの推定結果のIRF高速化 | `jcer_irf_bgvar()`。同じ事後ドローとBGVARのショック定義を維持 |
| 一般化IRF・明示的順序の再帰IRF | 実装、独立再帰計算との照合あり |
| 単一ショック符号制約 | Haar回転・受理数記録・不足時エラー |
| 複数国ワーカー | 国別シードで再現可能。起動時間込みで評価 |
| BGVAR比較 | 任意依存。事前・分散構造が異なる実務比較 |
| Normal-Gamma Structured VB / hybrid | 未実装。`method="vb"` はエラー |

条件付き国別事後分布の積を推定します。完全な世界同時尤度の事後分布とは異なり、
外国変数の弱外生性を仮定します。v0.1では国ごとの変数集合は共通です。
同一モデルでの推定法比較と、異なる事前分布のBGVARとの比較を区別します。
非定常ドローは既定で保持し、`stability="condition"`は棄却数と対象の変更を記録します。
符号制約の分位点はパラメータと回転の不確実性を混合し、識別集合の上下限とは異なります。

計測結果は `benchmarks/results/`、検証状況は `IMPLEMENTATION_STATUS.md`。
再実行は `Rscript benchmarks/run_benchmarks.R` と
`Rscript benchmarks/run_simulations.R`。`--full`で大規模グリッドを指定できます。
実測要約は検証実行後に本ファイル末尾へ記載します。

数学的仕様と制約は4本のvignette、`docs/MODEL_ASSUMPTIONS.md`、
`docs/BASELINE_REVIEW.md`を参照してください。
GPL-3。著者・保守担当者: motohiro watanabe (watanabemotosan@gmail.com)。
引用は `citation("preJCERBGVAR")`。参照論文と実装の対応は `CITATIONS.md`。

保存済み結果との互換性のため、S3クラス名 `JCERBGVAR` / `JCERBGVAR_irf` と関数名 `jcer_*` は維持しています。
下記の性能数値は名称変更前の0.1.1で取得した検証記録です。今回の開発版では名称・配布情報を変更し、計算方法は維持しています。


<!-- MEASURED_RESULTS -->
## 付属データでの実測

同じBGVAR事後分布・ショック定義・1000保持ドロー、各5反復の中央値。全体時間は同じ推定時間を含みます。

| 条件 | BGVAR IRF 秒 | 選択ベクトルIRF 秒 | IRF計時比 | 従来の全体 秒 | 置換後の全体 秒 | 全体計時比 |
|---|---:|---:|---:|---:|---:|---:|
| NG_4 |  4.41 |  0.27 | 16.33 |  8.95 |  4.61 | 1.941 |
| NG_12 | 11.72 |  0.59 | 19.86 | 26.15 | 14.91 | 1.754 |
| NG_20 | 26.59 |  1.25 | 21.27 | 53.91 | 28.57 | 1.887 |
| MN_4_long |   5.8 |  0.33 | 17.58 | 10.55 |  5.08 | 2.077 |
| NG_4_SV |  5.23 |  0.25 | 20.92 | 16.25 | 11.56 | 1.406 |

事前分布をMNIWへ変更した比較、安定性、合成DGPと計時範囲の修正は
[付属データ検証報告](docs/HISTORICAL_VALIDATION.md)を参照。旧3.6倍という合成比較は主張を差し替えています。
