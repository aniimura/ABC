/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Meta.Claim
import Mathlib.Data.Finset.Basic
import Mathlib.Data.Real.Basic
import Mathlib.Data.Fintype.Basic
import Mathlib.Algebra.BigOperators.Ring.Finset

/-!
# [CorrHyp] §1–§6 の語彙を受ける `Interface`

原典: S. Mochizuki, *Correspondences on Hyperbolic Curves* [CorrHyp]、物理 p.3–p.18。
`ResearchPaper/corrhyp-goal.md` の §1 0/5 〜 §6 0/1(合計 0/24)に対応するスケルトンの土台。

★フィールドはすべて原文に実際に現れる語と対応する。§0(Theorem A/B/C)は
`corrhyp-goal.md` のとおり後続定理の再掲なので対象外。

## ★★逐語との対応関係についての読み替え(逸脱、記録)

原文は双曲曲線 `X` と、その付随 Riemann 面 `𝒳`・普遍被覆 `𝒳̃`・解析的商 `𝒴`・
代数的商 `Y`・モジュライスタック `M_{g,r}` を書き分けるが、`ResearchPaper/papers.json`
(2026-09-04、200dpi 目視 p.3-12・p.17)のとおり pdftotext ではこれらがすべて同じ
ASCII 文字に潰れる。★原文自身も `Definition 2.1`-`2.3` で「X, 𝒳, or Γ」と三通りの
呼び方を並記しており、どれで述べても同じ性質だと明言している。ゆえに本 `Interface`
ではこれらを単一の `Space`(k 上の 1 次元(スタッキーもありうる)空間)として扱う。

## ★★Skeleton 段階での抽象化(逸脱、記録)

`Definition 2.2`(Margulis-arithmetic)・`Definition 2.3`(Shimura-arithmetic)は
代数群・Galois コホモロジー・四元数環を使う重い定義だが、本論文の以降の議論は
これらを**専ら「arithmetic かどうか」という述語としてしか使わない**
(実際に代数群 `G` を分解するのは `Proposition 2.4` 自身の証明の中だけであり、
その証明は本 Skeleton では `sorry` の対象——`Proposition 2.4` という 1 個の
`theorem` に閉じ込められる)。ゆえに `MargulisArithmetic` / `ShimuraArithmetic` を
抽象述語として posit する。Track B で中身を組み立てる時が来たら、この 2 つの
定義だけを差し替えれば済む設計にした。 -/

namespace ABC3.Interface.CorrHyp

open ABC3.Meta

/-- `Definition 5.2` の `Notation` 段落が導入する、スタック `Y` の「型」
`(g_Y; r_Y; {i_σ}_{σ∈Σ_Y})`。

原文 (CorrHyp p.12):
> Let us write `g_Y` for the genus of the compactification of `Y^c`, and `r_Y`
> for the number of points that need to be added to `Y^c` to compactify it. Let
> us write `Σ_Y` for the set of points of `Y^c` over which `Y → Y^c` is not
> étale.

★`Σ_Y` の有限性は原文に明記されないが、`Lemma 5.4`/`5.5` が `e_Y` の和を
有限和として扱っているので implicit に仮定されている。ここでは `Fintype` として
明示的に posit する(逸脱、記録)。 -/
structure StackType where
  /-- `g_Y`。 -/
  g : ℕ
  /-- `r_Y`。 -/
  r : ℕ
  /-- `Σ_Y`(分岐点の集合)。 -/
  Sigma : Type
  fintypeSigma : Fintype Sigma
  /-- `i_σ`(`σ ∈ Σ_Y` での `Y → Y^c` の分岐指数、常に `≥ 2`)。 -/
  i : Sigma → ℕ

attribute [instance] StackType.fintypeSigma

/-- `e_Y ≝ 2g_Y − 2 + r_Y + Σ_{σ∈Σ_Y} j_σ`、`j_σ ≝ (i_σ−1)/i_σ`(`Y` の Euler 標数)。

原文 (CorrHyp p.12):
> `e_Y ≝ 2g_Y − 2 + r_Y + Σ_{σ∈Σ_Y} j_σ` … Thus, one may think of `e_Y` as the
> Euler characteristic of `Y`. -/
noncomputable def StackType.e (t : StackType) : ℚ :=
  2 * (t.g : ℚ) - 2 + (t.r : ℚ) +
    ∑ σ ∈ (Finset.univ : Finset t.Sigma), ((t.i σ : ℚ) - 1) / (t.i σ : ℚ)

/-- [CorrHyp] §1–§6 の語彙を受ける `Interface`。フィールドの由来は各フィールドの
docstring、および `ABC3.Skeleton.CorrHyp.*` の各節点の docstring を見よ。 -/
structure HyperbolicCurveData where
  /-- k 上の双曲曲線・その商であるスタック(`Y`)・モジュライスタック(`M_{g,r}`)を
  まとめて表す台。原文の `X` / `𝒳` / `Y` / `𝒴` / `M_{g,r}` はすべてここに落ちる
  (ファイル冒頭の読み替えを見よ)。 -/
  Space : Type
  /-- 有限 étale 射(原文 `finite, étale morphism`)。 -/
  FEt : Space → Space → Type
  /-- 有限 étale 射の合成。 -/
  comp : ∀ {A B C : Space}, FEt A B → FEt B C → FEt A C
  /-- 有限 étale 射に沿ったファイバー積(`Definition 1.4` の `C₁ ×_Y C₂`)。 -/
  pullback : ∀ {A B C : Space}, FEt A C → FEt B C → Space
  pbFst : ∀ {A B C : Space} (f : FEt A C) (g : FEt B C), FEt (pullback f g) A
  pbSnd : ∀ {A B C : Space} (f : FEt A C) (g : FEt B C), FEt (pullback f g) B
  /-- 同型(原文の "up to isomorphism")。 -/
  Iso : Space → Space → Prop
  /-- 双曲曲線の型 `(g, r)`(原文 §1 冒頭)。 -/
  type : Space → ℕ × ℕ
  /-- `PSL₂(ℝ)⁰` の部分群 `Γ`(原文 §2、Fuchsian 群)。 -/
  Fuchsian : Type
  /-- `X ↦ Γ_X`(原文 §2、`Π ↪ Aut(H) = PSL₂(ℝ)⁰` の像)。 -/
  Gamma : Space → Fuchsian
  /-- `Comm(Γ)`(原文 §2 の可換化群)。 -/
  Comm : Fuchsian → Fuchsian
  /-- `Γ₁ が Γ₂ の中で有限指数`。 -/
  FiniteIndexIn : Fuchsian → Fuchsian → Prop
  /-- `Γ₁ ⊆ Γ₂`(部分群、`Proposition 3.2` で使う)。 -/
  Sub : Fuchsian → Fuchsian → Prop
  /-- `Definition 2.2` の Margulis-arithmetic(ファイル冒頭の抽象化を見よ)。 -/
  MargulisArithmetic : Fuchsian → Prop
  /-- `Definition 2.3` の Shimura-arithmetic(同上)。 -/
  ShimuraArithmetic : Fuchsian → Prop
  /-- `C` が連結(`Proposition 3.2` の前で「for simplicity」仮定される)。 -/
  IsConnected : Space → Prop
  /-- `hyperbolic core`(`Definition 3.1` / `Lemma 5.1` / `Definition 5.2`)。
  ★構成そのものは `Lemma 5.1` の内容なので、ここでは対象と自然な射だけを posit する。 -/
  core : Space → Space
  coreMap : (X : Space) → FEt X (core X)
  /-- `Space` の自己同型(`§4` の `Aut(X)`、`§6` の `M_{g,r}` の自己同型)。 -/
  Aut : Space → Type
  /-- 自明な自己同型(恒等射)。 -/
  idAut : (M : Space) → Aut M
  /-- `k` から `K` への係数拡大 `(−)_K`(`§4`)。 -/
  Ext : Space → Space
  extFEt : ∀ {A B : Space}, FEt A B → FEt (Ext A) (Ext B)
  /-- `Y` の `(g_Y; r_Y; {i_σ})`(`Definition 5.2` の `Notation`)。 -/
  stackType : Space → StackType
  /-- モジュライスタック `M_{g,r}`(`§6`)。 -/
  ModuliStack : ℕ → ℕ → Space
  /-- 「概スキームである」(`Theorem 6.1` の結論の前半)。 -/
  IsGenericallyScheme : Space → Prop
  /-- 有限 étale 射の次数(`§5` の `d`)。 -/
  deg : ∀ {A B : Space}, FEt A B → ℕ
  /-- `(M_{g,r})_k` の中で「開かつ稠密な部分スタック」に属する
  (`Theorem 5.3` の `U ⊆ (M_{g,r})_k`)。 -/
  IsOpenDenseIn : Set Space → ℕ → ℕ → Prop
  /-- `(M_{g,r})_k` の中で「構成可能」(`Lemma 5.6`)。 -/
  IsConstructibleIn : Set Space → ℕ → ℕ → Prop

end ABC3.Interface.CorrHyp
