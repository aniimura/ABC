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

★★★**3 例目(2026-08-17): auto-bound universe。**
`structure Foo where Point : Type u` と書いたとき、
`lake env lean` は **`u` を自動で束縛して通す**が、`lake build` は
`unknown universe level u` で落ちた。★**`universe u v` を明示する**こと。

★★これで「`lake env lean` と `lake build` が食い違う」実例は **3 種**になった:
セクション変数の自動包含 / auto-bound universe / `cd` の欠落。
**共通点は「`lake env lean` の方が緩い」**である。

## ★★★[[parallel-session-sweeps-my-files]] との調停

あちらは「`lake env lean` が通ったら**即コミット**」と言う(共有ワークツリーで
巻き込まれないため)。★★両立させる順番は:

1. `lake env lean FILE` が通る → ★**すぐコミットする**
2. その後 `lake build` を回す
3. 落ちたら **次のコミットで直す**

★★★**「build まで確かめてから 1 回で綺麗にコミット」を狙ってはならない**——
その間に並行セッションに巻き込まれる(2026-08-17 に実際に起きた)。

**How to apply:**
- 各コマンドの先頭に **絶対パスの `cd` を必ず書く**(`cd` の効果は次の呼び出しに残るので、
  前の呼び出しの cwd に依存しない)。
- 単一ファイルの反復には `lake env lean` を使ってよいが、
  ★**commit の前に必ず `lake build`(全体)を通す。**
- ゲートの出力は `tail` などにパイプしない —— パイプすると終了コードが `tail` のものになる。
  [[frdi-verbatim-ascii-only]] と同じく、**器具の出力を自分の目で読む**のが規律。
