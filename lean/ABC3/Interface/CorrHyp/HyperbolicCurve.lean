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
> morphism Y →Y c. Let us write gY for the genus of the compactification of Y c, and rY
> for the number of points that need to be added to Y c to compactify it. Let us write ΣY
> for the set of points of Y c over which Y →Y c is not ´etale.

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

原文 (CorrHyp p.12、`eY` の定義式に続く一文——式そのものは表組みで複数行に
分かれて抽出されるため逐語から外し、この一文だけを引く):
> Thus, one may think of eY as the Euler characteristic of Y . -/
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
  /-- 恒等射(`Found/CorrHyp/FuchsianGroup.lean` の `isFiniteIndexIn_refl` が
  具体化する——`Definition 1.5` 直後の「isogeny は同値関係」の反射性に要る)。 -/
  idFEt : ∀ A : Space, FEt A A
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
  /-- `PSL₂(ℝ)⁰` の部分群 `Γ`(原文 §2)。★離散性は要求しない——原文
  (p.5)自身が `Comm(Γ) := {γ ∈ PSL₂(ℝ)⁰ | (γ·Γ·γ⁻¹) ∼ Γ}` を「`PSL₂(ℝ)⁰`
  の部分群」として定義しており、離散性を課していない。`Γ` が伝統的な
  意味での Fuchsian 群(離散)であることは `IsDiscrete` で別途 posit する
  ——`Comm(Γ)` は non-arithmetic な `Γ` に対してのみ離散になる
  (Margulis の二分法、`corrhyp-goal.md` 2026-09-04 続報)。 -/
  Fuchsian : Type
  /-- `X ↦ Γ_X`(原文 §2、`Π ↪ Aut(H) = PSL₂(ℝ)⁰` の像)。 -/
  Gamma : Space → Fuchsian
  /-- `Γ` が(伝統的な意味で)離散であること。`Gamma X` は常にこれを満たす
  (`Gamma_isDiscrete`)。`Comm(Γ)` はこれを満たすとは限らない——満たす
  ことが `¬ MargulisArithmetic Γ ∧ ¬ ShimuraArithmetic Γ` と同値、という
  のが Margulis の二分法([Marg])の主張(`Theorem 2.5` 相当、まだ
  Skeleton には現れない補助事実)。 -/
  IsDiscrete : Fuchsian → Prop
  Gamma_isDiscrete : ∀ X : Space, IsDiscrete (Gamma X)
  /-- `Comm(Γ)`(原文 §2 の commensurator 群、p.5 の定義そのまま——
  `PSL₂(ℝ)⁰` の部分群として常に well-defined、離散とは限らない)。 -/
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

/-- ★**Track B は何を作らねばならないか。**

`StackType` は `g` と `r` を固定して `Sigma := Empty` とすれば非空虚(`e = 2g-2+r`)
なので、実は非空虚 witness 自体は作れる——ここでは代わりに、それが
「原文の `Y` の型と一致する」という**意味のある**非空虚性ではないことを明記する。 -/
def StackType.waiting : WaitingFor :=
  { what := "StackType.nonvacuousを構成すること自体は容易(Sigma := Emptyでe = 2g-2+r)だが、それは『Yの型として実際に出現するデータ』であることを何も保証しない。意味のある非空虚性はHyperbolicCurveDataのcore/stackTypeが具体的なCorrHyp/GenEllのような対象で実装されて初めて言える"
    trackB := "未着手。GenEllのFound/GenEll/相当のCorrHyp向け実装が要る" }

/-- ★**Track B は何を作らねばならないか。**

本構造体は双曲曲線・Fuchsian 群・モジュライスタック等、原文 §1-§6 のほぼ全語彙を
一度に posit している。個々のフィールドが「意味のある」対象(実際の双曲曲線の圏、
実際の `PSL₂(ℝ)` の部分群の圏、等)で実装されて初めて非空虚性に意味が出る。 -/
def HyperbolicCurveData.waiting : WaitingFor :=
  { what := "双曲曲線の圏(Space・FEt・pullback)、Fuchsian群の圏(PSL2(R)の部分群、Comm)、Margulis/Shimura-arithmeticの実装、モジュライスタックM_{g,r}の構成。2026-09-04時点でmathlibにhyperbolic curve・PSL2(R)・moduli stack of curvesの直接対応物は未実測"
    trackB := "未着手。Skeleton化(本ファイル群)が2026-09-04に先行——GenEllと同じ順序(Track A→Track B)" }

end ABC3.Interface.CorrHyp
