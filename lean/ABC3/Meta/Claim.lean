/-!
# 台帳の型 — 出典・実装待ち・飛躍

本プロジェクトは台帳を JSON ではなく **Lean の宣言そのもの**として持つ。
ノードとファイルが別々に生まれないので、被覆の乖離が構造的に起こりえない。

`tools/check.mjs` は、ここで定義した型の**存在**を機械的に検査する。
**型が付くかどうかは Lean が検査する**——この分業を守ること。
check.mjs は「`Foo.nonvacuous` という宣言があるか」しか見ておらず、
それが本当に `Nonempty Foo` を証明しているかは `lake build` が保証する。
-/

namespace ABC3.Meta

/-- 原典の一つの主張への参照。`1_Structured` の構造単位と 1:1 に対応する。

`sectionId` は `1_Structured/<論文>/section-N.html` の `<section id="...">`。
これにより Lean の宣言から、PDF 目視確認済みの逐語まで一本で辿れる。 -/
structure Source where
  /-- `ResearchPaper/papers.json` のタグ(例: `"pGC"`) -/
  paper : String
  /-- **物理**ページ番号(`pdftoppm -f N` が指すページ) -/
  pdfPage : Nat
  /-- 項目名(例: `"Proposition 1.1"`) -/
  item : String
  /-- `1_Structured` の `<section id="...">` -/
  sectionId : String
  deriving Repr, DecidableEq

/-- **load-bearing** の印——主定理へ実際に消費される系列であること。

`Skeleton` の宣言 `foo` に `foo.loadBearing` を付けると、`tools/check.mjs` は
その宣言に **G3(負の対照)** `foo.negControl` を要求する。

荷重配分の規律: G1〜G5 全部を全 skeleton に課すのは有害(薄く均一に儀式を撒くと、
効かせるべき箇所の検査が相対的に軽くなる)。**この印を付けたものにだけ全ゲートを課す。** -/
structure LoadBearing where
  /-- どの結果へ向けて消費されるか(消費者を名指しできないなら load-bearing ではない) -/
  consumer : String
  deriving Repr

/-- **負の対照**(G3)の記録。

主張している性質を **1つだけ** 落とした対照が、実際に破れることを確認した記録。
破れないなら、その性質は statement に効いていない——飾りである。

`witness` には、対照が破れることを示す Lean の宣言名を書く。 -/
structure NegControl where
  /-- 落とした性質(1つだけ) -/
  dropped : String
  /-- 対照が破れることを示す宣言名 -/
  witness : String
  deriving Repr

/-- `Interface` の構造体が、まだ非空虚 witness を持てない理由。

**空欄で済ませないための型**——「まだ作れない」は許されるが、
「何を待っているか言えない」は許されない。 -/
structure WaitingFor where
  /-- 何の実装を待っているか(例: `"局所類体論の相互法則"`) -/
  what : String
  /-- Track B のどの作業に対応するか -/
  trackB : String
  deriving Repr

/-- 飛躍の分類。**既定は `modelError`**——`sourceGap` は最後の手段。 -/
inductive GapClass
  /-- ① こちらのモデル化の誤り -/
  | modelError
  /-- ② 必要な数学が未構築 -/
  | missingMath
  /-- ③ 原典側の飛躍。`falsifier` を書けないなら、これを名乗ってはならない -/
  | sourceGap
  deriving Repr, DecidableEq

/-- 原典が確立していないと判明した段の記録。

`Gap/` の各構造体はこの記録を伴う。**`falsifier` は必須**——
「何が起きればこれが①②に落ちるか」を書けないうちは、まだ③ではない。 -/
structure GapRecord where
  source : Source
  classification : GapClass
  /-- 何が起きればこの分類が覆るか -/
  falsifier : String
  deriving Repr

end ABC3.Meta
