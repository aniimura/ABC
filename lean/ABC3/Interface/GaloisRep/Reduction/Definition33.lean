/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Interface.GaloisRep.Representation
import Mathlib.NumberTheory.NumberField.Basic
import Mathlib.AlgebraicGeometry.EllipticCurve.Reduction
import Mathlib.Analysis.SpecialFunctions.Log.Basic
import Mathlib.GroupTheory.QuotientGroup.Basic
import Mathlib.NumberTheory.NumberField.InfinitePlace.Embeddings
import Mathlib.NumberTheory.ModularForms.Discriminant
import Mathlib.NumberTheory.ModularForms.EisensteinSeries.Basic

/-!
# Reduction —— `[GenEll] Definition 3.3` の分

☆もとの 1 枚を**条なしの項目ごと**に割ったものである（第 1458、案 a）。
★「1 ファイル = 1 ノード」を回復するための分割で、中身は動かしていない。
-/

namespace ABC3.Interface.GaloisRep

open ABC3.Meta WeierstrassCurve NumberField

/-! ## ★★★G6 —— Tate 曲線と局所高さ -/

/-- **(G6)** 乗法還元をもつ楕円曲線の **Tate 母数 `q_E`** と**局所高さ** `v_K(q_E)`。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

★★★**原文の `Definition 3.3` そのものである。**
★`Lemma 3.2`(局所の階数 1 部分群)がこの上で述べられる。

★★★**2026-08-17: 退化検査で作り直した。**以前は `LocalField : Type` /
`Curve : LocalField → Type` と**世界ごと posit** していたので、
`PUnit` と定数 `1` で埋まってしまった(実測済)。
★**mathlib の `WeierstrassCurve` と `Valuation` に接地し**、
**Tate 一意化そのもの**を要求する形に改めた。 -/
structure TateCurveData where
  /-- ★★**Tate 母数** `q_E ∈ Kˣ`。

  ★★★**述語は posit しない**——mathlib の
  `WeierstrassCurve.HasSplitMultiplicativeReduction` に接地する。
  ★これで「`HasSplitMultRed := False`」という**空虚 witness も不可能**になる。 -/
  tateParam : {R : Type} → [CommRing R] → [IsDomain R] → [IsDiscreteValuationRing R] →
    [IsAdicComplete (IsLocalRing.maximalIdeal R) R] →
    {K : Type} → [Field K] → [Algebra R K] → [IsFractionRing R K] →
    (W : WeierstrassCurve K) → [W.IsElliptic] → [W.IsMinimal R] →
    W.HasSplitMultiplicativeReduction R → Kˣ
  /-- ★★★**Tate 一意化** `E(K) ≅ Kˣ / q^ℤ`。

  ★★★**これが内容の本体である。**`q` が実際に一意化していることを要求する。 -/
  uniformization : ∀ {R : Type} [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
    [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
    {K : Type} [Field K] [DecidableEq K] [Algebra R K] [IsFractionRing R K]
    (W : WeierstrassCurve K) [W.IsElliptic] [W.IsMinimal R]
    (h : W.HasSplitMultiplicativeReduction R),
    Nonempty (W.toAffine.Point ≃+
      Additive (Kˣ ⧸ Subgroup.zpowers (tateParam W h : Kˣ)))
  /-- ★★★**局所高さ** `v_K(q_E) ∈ ℤ_{>0}`(原文 `Definition 3.3`)。

  ★★★★**2026-08-26 の訂正(逸脱の記録)**——以前は付値
  `v : Kˣ →* Multiplicative ℤ` を**任意に受けていた**。しかしそれだと
  `v := 1`(自明な準同型、`toAdd (1 q) = 0`)を入れたとき
  `localHeight_eq` が `localHeight = 0` を、`localHeight_pos` が `0 < localHeight` を
  要求して**衝突する**——すなわち**構造そのものが充足不能**だった。
  ★離散付値環の**正規化付値**(mathlib の `IsDiscreteValuationRing.maximalIdeal` の
  adic 付値)に固定して直す。★★これは**弱めるのではなく強める**訂正である
  (`localHeight` が `q` から一意に決まる)。

  ★★★★**2026-08-26 の第 2 の訂正**——`Δ ≠ 0`(`IsElliptic`)を要求する。
  mathlib の `HasMultiplicativeReduction` は `v(Δ) < 1` としか言わず、
  **`Δ = 0` を排除しない**(`v(0) = 0 < 1`)。実際 `y² + xy = x³` は
  `Δ = 0`・`c₄ = 1`・接線の 2 次式 `X² + X` が分解、で**すべての仮定を満たす**が、
  これは**結節三次曲線**であり滑らかな点の群は `Kˣ` そのもの——
  Tate 母数は `q = 1` にあたるので**局所高さが `0`** になり `localHeight_pos` と衝突する。 -/
  localHeight : {R : Type} → [CommRing R] → [IsDomain R] → [IsDiscreteValuationRing R] →
    [IsAdicComplete (IsLocalRing.maximalIdeal R) R] →
    {K : Type} → [Field K] → [Algebra R K] → [IsFractionRing R K] →
    (W : WeierstrassCurve K) → [W.IsElliptic] → [W.IsMinimal R] →
    W.HasSplitMultiplicativeReduction R → ℕ
  /-- ★★★**局所高さは Tate 母数の付値そのもの**——定義を型に出す。

  ★正規化付値では `v(q) = ofAdd (−localHeight)`(乗法的な書き方)である。
  ★★これがあるので `localHeight := 1` のような定数 witness は通らない。 -/
  localHeight_eq : ∀ {R : Type} [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
    [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
    {K : Type} [Field K] [Algebra R K] [IsFractionRing R K]
    (W : WeierstrassCurve K) [W.IsElliptic] [W.IsMinimal R]
    (h : W.HasSplitMultiplicativeReduction R),
    IsDedekindDomain.HeightOneSpectrum.valuation K (IsDiscreteValuationRing.maximalIdeal R)
        ((tateParam W h : Kˣ) : K)
      = (Multiplicative.ofAdd (-(localHeight W h : ℤ)) : Multiplicative ℤ)
  /-- ★局所高さは正(原文が `∈ Z_{>0}` と書いている)。 -/
  localHeight_pos : ∀ {R : Type} [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
    [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
    {K : Type} [Field K] [Algebra R K] [IsFractionRing R K]
    (W : WeierstrassCurve K) [W.IsElliptic] [W.IsMinimal R]
    (h : W.HasSplitMultiplicativeReduction R), 0 < localHeight W h

def TateCurveData.waiting : WaitingFor :=
  { what := "(G6) Tate 曲線 E(K̄) ≅ K̄^x/q^Z、Tate 母数 q_E、および局所高さ v_K(q_E)(= [GenEll] Definition 3.3)"
    trackB := "Found/GaloisRep — ★mathlib に Tate 曲線は無い(2026-08-16 実測)。★★**FLT にも無い**——`FLT/TateCurve/TateCurve.lean` は 20 行の入口宣言だけ(2026-08-17、clone して計数)。★★★原文の典拠は [FC](Faltings-Chai, *Degenerations of Abelian Varieties*)Chapter III, Corollary 7.3 で、これも mathlib に無い" }

/-! ## ★★G7 —— 半安定還元とモデル -/

/-- **半安定**であること —— 離散付値環 `R` の上で**良還元または乗法還元**。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

★★★**posit ではなく mathlib の還元型による定義である。**
`Reduction.lean` の三分律(良・乗法・加法)のうち加法還元を除いたものが半安定。 -/
def IsSemiStableAt (R : Type) [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
    {K : Type} [Field K] [Algebra R K] [IsFractionRing R K]
    (W : WeierstrassCurve K) [W.IsMinimal R] : Prop :=
  W.HasGoodReduction R ∨ W.HasMultiplicativeReduction R

/-- ★半安定とは「加法還元でない」ことである(三分律から)。 -/
theorem isSemiStableAt_iff_not_additive (R : Type) [CommRing R] [IsDomain R]
    [IsDiscreteValuationRing R] {K : Type} [Field K] [Algebra R K] [IsFractionRing R K]
    (W : WeierstrassCurve K) [W.IsMinimal R] :
    IsSemiStableAt R W ↔ ¬ W.HasAdditiveReduction R := by
  constructor
  · rintro (h | h)
    · exact h.not_hasAdditiveReduction
    · exact h.not_hasAdditiveReduction
  · intro h
    rcases WeierstrassCurve.hasGoodReduction_or_hasMultiplicativeReduction_or_hasAdditiveReduction
      R (W := W) with hg | hm | ha
    · exact Or.inl hg
    · exact Or.inr hm
    · exact absurd ha h

/-- **(G7)** 半安定還元と、`𝓞_L` 上への延長が与える **Néron 微分** `ω_E`。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

★★原文は「`E` は `L` のすべての有限素点で半安定還元をもつ」を繰り返し仮定し、
そこから `𝓞_L` 上の 1 次元半アーベルスキームへの延長を得る。

★★★**その延長が使われるのは `ω_E` を作るためだけ**である
(`Proposition 3.4` の `ht^Falt = deg(ω_E)`)。ゆえに `ω_E` を型に出す。

★★★**`ω_E` は階数 1** —— これが退化(`PUnit` は階数 0)を殺す。 -/
structure SemistableModelData where
  /-- 台となる Tate 曲線。 -/
  toTateCurveData : TateCurveData
  /-- ★★**Néron 微分が定める分数イデアル** `ω_E ⊆ L`。

  ★★★★**2026-08-26 の訂正(4 つ目の穴を塞ぐ)**——以前は
  `omega : … → Type` と**階数 1 だけ**を課していた。
  `Check/GaloisRep/OmegaNondegenerate.lean` が示したとおり、それは
  **`ω_E := 𝓞_L`(曲線を無視する定数)で満たせてしまう**。
  ★ここでは `ω_E` を **`L` の分数イデアル**として持ち、
  **変数変換での変わり方**を課す——これが「その欄と入力データの両方を言及する条件」である。

  ★★不変微分 `ω = dx/(2y + a₁x + a₃)` は変数変換で `ω' = u·ω` と変わるので、
  Néron 微分を `ω` で割った係数は `u⁻¹` 倍され、分数イデアルも `u⁻¹` 倍される。 -/
  omegaFrac : (L : Type) → [Field L] → [NumberField L] → WeierstrassCurve L →
    Submodule (𝓞 L) L
  /-- ★分数イデアルであること(有限生成)——`L` 全体は f.g. でないので、これが効く。

  ★★★★**2026-08-26 の訂正**——`[W.IsElliptic]` を付けた。`Δ = 0` の曲線に
  Néron 微分は無い(第 304 で `TateCurveData` に入れたのと同じ理由)。
  ★欄そのものは全曲線で定義しておき、**性質だけを楕円曲線に限る**。 -/
  omegaFrac_fg : ∀ (L : Type) [Field L] [NumberField L] (W : WeierstrassCurve L) [W.IsElliptic],
    (omegaFrac L W).FG
  /-- ★零でない(階数 1 の代わり)。 -/
  omegaFrac_ne_bot : ∀ (L : Type) [Field L] [NumberField L] (W : WeierstrassCurve L)
    [W.IsElliptic], omegaFrac L W ≠ ⊥
  /-- ★★★★**変数変換で `u⁻¹` 倍される**——これが `ω_E` を曲線に結びつける条件。 -/
  omegaFrac_variableChange : ∀ (L : Type) [Field L] [NumberField L] (W : WeierstrassCurve L)
    [W.IsElliptic] (C : WeierstrassCurve.VariableChange L) (x : L),
    x ∈ omegaFrac L (C • W) ↔ ((C.u : L) * x) ∈ omegaFrac L W
  /-- 大域の半安定性(原文が繰り返し仮定するもの)。 -/
  SemiStable : (L : Type) → [Field L] → [NumberField L] → WeierstrassCurve L → Prop
  /-- ★★★**その意味を固定する**——各離散付値環の上で半安定であること。

  ★これがあるので `SemiStable := True` のような自由な posit にはならない。 -/
  semiStable_iff : ∀ (L : Type) [Field L] [NumberField L] (W : WeierstrassCurve L),
    SemiStable L W ↔
      ∀ (R : Type) [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
        [Algebra R L] [IsFractionRing R L] [W.IsMinimal R], IsSemiStableAt R W

def SemistableModelData.waiting : WaitingFor :=
  { what := "(G7) 半安定還元(★mathlib の還元型で定義済——posit ではない)と、𝓞_L 上への延長が与える Néron 微分の分数イデアル omega_E(★2026-08-26: 変数変換での変わり方を課して曲線に結びつけた)"
    trackB := "Found/GaloisRep — ★mathlib は `AlgebraicGeometry/EllipticCurve/Reduction.lean` を持つ(2026-08-17 実測)が、**Néron モデル・半アーベルスキームは無い**。★★原文は延長の存在を [FC] Chapter I, Proposition 2.7 に帰している" }

/-! ## ★出典の紐付け(`.src`) -/

def TateCurveData.src : Source :=
  { paper := "GenEll", pdfPage := 15, item := "Definition 3.3(Tate 母数と局所高さの定式化のみ)",
    sectionId := "genell-def-3-3" }

def SemistableModelData.src : Source :=
  { paper := "GenEll", pdfPage := 15, item := "Definition 3.3(半安定還元と 𝓞_L 上のモデルのみ)",
    sectionId := "genell-def-3-3" }

end ABC3.Interface.GaloisRep
