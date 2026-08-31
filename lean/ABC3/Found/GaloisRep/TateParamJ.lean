/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.TateJInv
import ABC3.Found.GaloisRep.TateCurveWitness
import ABC3.Found.GaloisRep.LocalHeightDelta

/-!
# Galois (G6) 第 882 ブロック —— **★★★★★★★★★★曲線の水準で「Tate 母数は `j` で決まる」**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> parameter qE of E satisfies the relation qE = qEl ; in particular, we have

## ★★★★★★★★★★これは何か

第 878 の `tateParam_injective` は**Tate 曲線同士**で述べていた。
しかし `Lemma 3.5` が実際に持っているのは、**`K` の上の曲線 `W`・`W′`**と
その Tate 母数 `tateParamR W h` である。

★本ブロックはその隔たりを埋める:

    `W.j = W′.j`  ⟹  `tateParamR W h = tateParamR W′ h′`

☆道は 3 段である。

| 段 | 中身 | 出どころ |
|---|---|---|
| 1 | `j` の一致を分母を払った形に直す | `c4_cube_mul_Delta_of_j_eq`（第 879） |
| 2 | `K` の等式を整モデルに戻す | `integralModel_c₄_eq`・`integralModel_Δ_eq`（mathlib） |
| 3 | 変数変換で Tate モデルに移す | `variableChange_c₄`・`variableChange_Δ`（mathlib） |

★★**変数変換の `u` は消える**——`Δ ↦ u⁻¹¹²Δ`、`c₄ ↦ u⁻¹⁴c₄` なので
`Δ·c₄′³` と `Δ′·c₄³` は**同じ重み `u⁻¹¹²u′⁻¹¹²`** を拾う。
☆だから `j` を経由せずに直接**積の形**で払える（`IsElliptic` のインスタンスが
`rw` の motive を壊す問題を避けられる）。

| 定理 | 内容 |
|---|---|
| `tateModel_c₄` / `tateModel_Δ` | Tate モデルの `c₄`・`Δ` を整モデルで書く |
| `integralModel_c4_cube_mul_Delta` | ★分母を払った等式を `R` に戻す |
| `tateParamR_eq_of_j` | ★★★★★★★★★★**`j` が同じなら Tate 母数も同じ** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve IsDedekindDomain

variable {R : Type} [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
  [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
  {K : Type} [Field K] [CharZero K] [Algebra R K] [IsFractionRing R K]

/-! ## ★整モデルに戻す -/

/-- ★**`j` の一致を整モデルの等式に戻す**。

☆`algebraMap R K` は単射（`IsFractionRing.injective`）なので、
`K` で成り立つ多項式的な等式は `R` に降りる。 -/
theorem integralModel_c4_cube_mul_Delta (W W' : WeierstrassCurve K)
    [W.IsElliptic] [W'.IsElliptic] [IsIntegral R W] [IsIntegral R W']
    (hj : W.j = W'.j) :
    (integralModel R W).c₄ ^ 3 * (integralModel R W').Δ
      = (integralModel R W').c₄ ^ 3 * (integralModel R W).Δ := by
  have hK := ABC3.Found.GenEll.c4_cube_mul_Delta_of_j_eq W W' hj
  refine IsFractionRing.injective R K ?_
  rw [map_mul, map_mul, map_pow, map_pow, integralModel_c₄_eq, integralModel_c₄_eq,
    integralModel_Δ_eq, integralModel_Δ_eq]
  exact hK

/-! ## ★★★★★★★★★★Tate 母数は `j` で決まる -/

/-- ★★★★★★★★★★**`j` が同じなら Tate 母数も同じ**（曲線の水準）。

原文 (GenEll p.15):
> parameter qE of E satisfies the relation qE = qEl ; in particular, we have

★これが `Lemma 3.5` の葉 1「`E′` の Tate 母数は `q^l`」を
**`j` の一致 1 本**に落とす道具である。 -/
theorem tateParamR_eq_of_j (W W' : WeierstrassCurve K)
    [W.IsElliptic] [W'.IsElliptic] [W.IsMinimal R] [W'.IsMinimal R]
    (h : W.HasSplitMultiplicativeReduction R) (h' : W'.HasSplitMultiplicativeReduction R)
    (hj : W.j = W'.j) :
    tateParamR W h = tateParamR W' h' := by
  obtain ⟨hq, C, hne, hCE⟩ := tateParamR_spec W h
  obtain ⟨hq', C', hne', hCE'⟩ := tateParamR_spec W' h'
  have hM := integralModel_c4_cube_mul_Delta (R := R) W W' hj
  -- ★Tate モデルの `c₄`・`Δ` を整モデルで書く
  have hc4 : (tateCurveAt (tateParamR W h) hq).c₄
      = ((C.u⁻¹ : Rˣ) : R) ^ 4 * (integralModel R W).c₄ := by
    rw [← hCE, WeierstrassCurve.variableChange_c₄]
  have hc4' : (tateCurveAt (tateParamR W' h') hq').c₄
      = ((C'.u⁻¹ : Rˣ) : R) ^ 4 * (integralModel R W').c₄ := by
    rw [← hCE', WeierstrassCurve.variableChange_c₄]
  have hd : (tateCurveAt (tateParamR W h) hq).Δ
      = ((C.u⁻¹ : Rˣ) : R) ^ 12 * (integralModel R W).Δ := by
    rw [← hCE, WeierstrassCurve.variableChange_Δ]
  have hd' : (tateCurveAt (tateParamR W' h') hq').Δ
      = ((C'.u⁻¹ : Rˣ) : R) ^ 12 * (integralModel R W').Δ := by
    rw [← hCE', WeierstrassCurve.variableChange_Δ]
  refine tateParam_injective (I := IsLocalRing.maximalIdeal R) hq hq' ?_
  rw [hc4, hc4', hd, hd']
  linear_combination
    (-(((C.u⁻¹ : Rˣ) : R) ^ 12) * ((C'.u⁻¹ : Rˣ) : R) ^ 12) * hM

/-- ★★★★★★★★★★★★**`W` の `j` が `E_{q₀}` の `j` なら Tate 母数は `q₀`**。

原文 (GenEll p.15):
> parameter qE of E satisfies the relation qE = qEl ; in particular, we have

★これが `Lemma 3.5` の葉 1 が実際に使う形である——
`E′` が Vélu の商と同じ `j` を持つことだけから `q_{E′} = q^l` が出る。

☆道は `tateParamR_eq_of_j` と同じだが、相手が**直接 Tate 曲線**なので
向こう側の変数変換が要らない。 -/
theorem tateParamR_eq_of_j_tateCurveAt (W : WeierstrassCurve K)
    [W.IsElliptic] [W.IsMinimal R] (h : W.HasSplitMultiplicativeReduction R)
    (q₀ : R) (hq₀ : q₀ ∈ IsLocalRing.maximalIdeal R)
    [((tateCurveAt q₀ hq₀).map (algebraMap R K)).IsElliptic]
    (hj : W.j = ((tateCurveAt q₀ hq₀).map (algebraMap R K)).j) :
    tateParamR W h = q₀ := by
  obtain ⟨hq, C, hne, hCE⟩ := tateParamR_spec W h
  -- ★段 1: `j` の一致を分母を払った形にし、`R` に降ろす
  have hK := ABC3.Found.GenEll.c4_cube_mul_Delta_of_j_eq W
    ((tateCurveAt q₀ hq₀).map (algebraMap R K)) hj
  rw [WeierstrassCurve.map_c₄, WeierstrassCurve.map_Δ] at hK
  have hR : (integralModel R W).c₄ ^ 3 * (tateCurveAt q₀ hq₀).Δ
      = (tateCurveAt q₀ hq₀).c₄ ^ 3 * (integralModel R W).Δ := by
    refine IsFractionRing.injective R K ?_
    rw [map_mul, map_mul, map_pow, map_pow, integralModel_c₄_eq, integralModel_Δ_eq]
    exact hK
  -- ★段 2: Tate モデルの `c₄`・`Δ` を整モデルで書く
  have hc4 : (tateCurveAt (tateParamR W h) hq).c₄
      = ((C.u⁻¹ : Rˣ) : R) ^ 4 * (integralModel R W).c₄ := by
    rw [← hCE, WeierstrassCurve.variableChange_c₄]
  have hd : (tateCurveAt (tateParamR W h) hq).Δ
      = ((C.u⁻¹ : Rˣ) : R) ^ 12 * (integralModel R W).Δ := by
    rw [← hCE, WeierstrassCurve.variableChange_Δ]
  -- ★段 3: 第 878 の単射性
  refine tateParam_injective (I := IsLocalRing.maximalIdeal R) hq hq₀ ?_
  rw [hc4, hd]
  linear_combination (-(((C.u⁻¹ : Rˣ) : R) ^ 12)) * hR

/-! ## ★★★★★★★★`v(1/j) = v(q)` -/

/-- ★★★★★★★★**`1/j` の付値は `q` の付値に等しい**。

原文 (GenEll p.15):
> parameter qE of E satisfies the relation qE = qEl ; in particular, we have

☆第 930 の `evalAdic tateJinvSeries q hq = q·u`（`u` は `R` の単元）と、
単元の付値が `0` であること（`tateDvrVal_eq_zero_of_isUnit`、証明済み）から即座である。

★これが `v_p(j) = −v_p(q)` の中身であり、
`Lemma 3.5` の連鎖の結論（`j` の一致）を `jExp` に繋ぐ橋である。 -/
theorem vAdd_evalAdic_tateJinvSeries (q : R) (hq : q ∈ IsLocalRing.maximalIdeal R)
    (hq0 : q ≠ 0)
    (hne : algebraMap R K (evalAdic tateJinvSeries q hq) ≠ 0)
    (hqne : algebraMap R K q ≠ 0) :
    vAdd (tateDvrVal R K)
        (Units.mk0 (algebraMap R K (evalAdic tateJinvSeries q hq)) hne)
      = vAdd (tateDvrVal R K) (Units.mk0 (algebraMap R K q) hqne) := by
  obtain ⟨u, hu, hqu⟩ := evalAdic_tateJinvSeries_eq_mul_unit (I := IsLocalRing.maximalIdeal R)
    q hq
  have hune : algebraMap R K u ≠ 0 := by
    intro h0
    apply hne
    rw [hqu, map_mul, h0, mul_zero]
  have hsplit : Units.mk0 (algebraMap R K (evalAdic tateJinvSeries q hq)) hne
      = Units.mk0 (algebraMap R K q) hqne * Units.mk0 (algebraMap R K u) hune := by
    refine Units.ext ?_
    show algebraMap R K (evalAdic tateJinvSeries q hq)
      = algebraMap R K q * algebraMap R K u
    rw [hqu, map_mul]
  rw [hsplit, vAdd_mul, tateDvrVal_eq_zero_of_isUnit (K := K) u hu hune, add_zero]

/-! ## ★出典の紐付け(`.src`) -/

def vAdd_evalAdic_tateJinvSeries.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(1/j の付値は q の付値に等しい。★無条件)",
    sectionId := "genell-lemma-3-2" }

def integralModel_c4_cube_mul_Delta.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(j の一致を整モデルの等式に戻す。★無条件)",
    sectionId := "genell-lemma-3-2" }

def tateParamR_eq_of_j_tateCurveAt.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(j が E_{q₀} の j なら Tate 母数は q₀。★無条件)",
    sectionId := "genell-lemma-3-2" }

def tateParamR_eq_of_j.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(j が同じなら Tate 母数も同じ——曲線の水準。★無条件)",
    sectionId := "genell-lemma-3-2" }

end ABC3.Found.GaloisRep
