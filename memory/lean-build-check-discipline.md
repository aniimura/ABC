---
name: lean-build-check-discipline
description: 合否判定は必ず `cd /d/Math_ABC3/lean && lake build`。`lake env lean FILE` が通っても build は落ちうる。
metadata:
  type: feedback
---

Lean の合否は **`cd /d/Math_ABC3/lean && lake build`** で判定する。ゲートは
**`cd /d/Math_ABC3 && node tools/check.mjs`** を単独で実行し、終了コードと `PASS` を目視する。

**Why:** 2026-08-17 に 2 回同じ事故を起こした。
1. `cd` を付けずに `lake build` を走らせ、`D:\Math_ABC3\lakefile.lean が無い` で
   失敗したのを見落としたまま commit した(壊れたビルドを push)。
   ★**並行セッションの子が全体ビルドで検出**し、原因まで特定してくれた。
2. `lake env lean FILE` は通るのに `lake build` が落ちた
   (セクション変数 `F` が仮定の型に現れて自動で含まれる、という差が出た)。

**How to apply:**
- 各コマンドの先頭に **絶対パスの `cd` を必ず書く**(`cd` の効果は次の呼び出しに残るので、
  前の呼び出しの cwd に依存しない)。
- 単一ファイルの反復には `lake env lean` を使ってよいが、
  ★**commit の前に必ず `lake build`(全体)を通す。**
- ゲートの出力は `tail` などにパイプしない —— パイプすると終了コードが `tail` のものになる。
  [[frdi-verbatim-ascii-only]] と同じく、**器具の出力を自分の目で読む**のが規律。
