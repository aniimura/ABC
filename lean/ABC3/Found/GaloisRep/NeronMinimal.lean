import ABC3.Found.GaloisRep.NeronExp

/-!
# Galois (G7) 第 313 ブロック —— **★★★★★★`v(Δ) < 12` なら極小**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★到達点

> 整モデルで **`v(Δ) < 12`** なら、その曲線は**すでに極小**である
> (`isMinimal_of_vAdd_Delta_lt`)。したがって `minimalExp = 0`。

★★★これが「**ほとんどの素点で Néron 指数が `0`**」の芯である——
`v(Δ) = 0` の素点(有限個を除く全部)では、この判定がそのまま効く。

## ★★★★★なぜ `12` か——`Δ` は `u⁻¹²` で動く

変数変換で `Δ ↦ u⁻¹²Δ`、すなわち加法付値では

    v(Δ') = v(Δ) − 12·v(u)

★整モデルなら `v(Δ') ≥ 0` だから `12·v(u) ≤ v(Δ) < 12`、よって **`v(u) ≤ 0`**。
★★一方「`Δ` の付値を下げられる」なら `v(u) ≥ 0`。合わせて **`v(u) = 0`**——
つまり判別式の付値は下がらない。

## ★★★★極小性は「最大性」なので、比較で示す

mathlib の `IsMinimal` は `MaximalFor`(整モデルの中で `v(Δ)` が最大)。
★示すべきは「どの整モデル `C • W` も `v(Δ)` を真に下げられない」ことで、
上の不等式 2 本を `omega` に渡すだけで済む。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `valuation_le_iff_vAdd_le` | ★★乗法的付値の順序と `vAdd` の順序 |
| `integral_Delta_vAdd_nonneg` | ★★整モデルの `v(Δ) ≥ 0` |
| `vAdd_Delta_variableChange` | ★★★`v(Δ') = v(Δ) − 12v(u)` |
| `isMinimal_of_vAdd_Delta_lt` | ★★★★★★**`v(Δ) < 12` なら極小** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve

variable {R : Type} [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
  {K : Type} [Field K] [Algebra R K] [IsFractionRing R K]

/-- ★★乗法的付値の順序は `vAdd` の逆順序。 -/
theorem valuation_le_iff_vAdd_le (x y : Kˣ) :
    (IsDiscreteValuationRing.maximalIdeal R).valuation K (x : K)
        ≤ (IsDiscreteValuationRing.maximalIdeal R).valuation K (y : K)
      ↔ vAdd (tateDvrVal R K) y ≤ vAdd (tateDvrVal R K) x := by
  rw [valuation_eq_ofAdd_neg_vAdd, valuation_eq_ofAdd_neg_vAdd, WithZero.coe_le_coe,
    Multiplicative.ofAdd_le]
  omega

/-- ★★整モデルの判別式の付値は非負。 -/
theorem integral_Delta_vAdd_nonneg (W : WeierstrassCurve K) [IsIntegral R W] (hΔ : W.Δ ≠ 0) :
    0 ≤ vAdd (tateDvrVal R K) (Units.mk0 W.Δ hΔ) := by
  refine tateDvrVal_nonneg_of_mem _ ⟨(integralModel R W).Δ, ?_⟩
  exact WeierstrassCurve.integralModel_Δ_eq R W

/-- ★★★変数変換での判別式の付値の動き。 -/
theorem vAdd_Delta_variableChange (W : WeierstrassCurve K) (hΔ : W.Δ ≠ 0) (C : VariableChange K) :
    vAdd (tateDvrVal R K) (Units.mk0 ((C • W).Δ) (variableChange_Delta_ne_zero W hΔ C))
      = vAdd (tateDvrVal R K) (Units.mk0 W.Δ hΔ) - 12 * vAdd (tateDvrVal R K) C.u := by
  have hu : (Units.mk0 ((C • W).Δ) (variableChange_Delta_ne_zero W hΔ C))
      = C.u⁻¹ ^ 12 * Units.mk0 W.Δ hΔ := by
    refine Units.ext ?_
    show (C • W).Δ = _
    rw [WeierstrassCurve.variableChange_Δ]
    push_cast
    simp
  rw [hu, vAdd_mul, vAdd_pow, vAdd_inv]
  omega

set_option maxHeartbeats 1600000 in
/-- ★★★★★★**`v(Δ) < 12` の整モデルは極小**。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem isMinimal_of_vAdd_Delta_lt (W : WeierstrassCurve K) [hint : IsIntegral R W]
    (hΔ : W.Δ ≠ 0) (h : vAdd (tateDvrVal R K) (Units.mk0 W.Δ hΔ) < 12) : IsMinimal R W := by
  rw [WeierstrassCurve.isMinimal_iff]
  constructor
  · show IsIntegral R ((1 : VariableChange K) • W)
    rw [one_smul]
    exact hint
  · intro C hCint hle
    have hCint' : IsIntegral R (C • W) := hCint
    have hCΔ : (C • W).Δ ≠ 0 := variableChange_Delta_ne_zero W hΔ C
    have h0 : WeierstrassCurve.valuation_Δ_aux R ((1 : VariableChange K) • W)
        ≤ WeierstrassCurve.valuation_Δ_aux R (C • W) := hle
    rw [one_smul] at h0
    rw [WeierstrassCurve.valuation_Δ_aux, dif_pos hint, WeierstrassCurve.valuation_Δ_aux,
      dif_pos hCint'] at h0
    have hle2 : (IsDiscreteValuationRing.maximalIdeal R).valuation K W.Δ
        ≤ (IsDiscreteValuationRing.maximalIdeal R).valuation K ((C • W).Δ) := h0
    have hbridge : vAdd (tateDvrVal R K) (Units.mk0 ((C • W).Δ) hCΔ)
        ≤ vAdd (tateDvrVal R K) (Units.mk0 W.Δ hΔ) := by
      rw [← valuation_le_iff_vAdd_le]
      exact hle2
    have hnn : 0 ≤ vAdd (tateDvrVal R K) (Units.mk0 ((C • W).Δ) hCΔ) :=
      integral_Delta_vAdd_nonneg (R := R) (C • W) hCΔ
    have hchg := vAdd_Delta_variableChange (R := R) W hΔ C
    show WeierstrassCurve.valuation_Δ_aux R (C • W)
      ≤ WeierstrassCurve.valuation_Δ_aux R ((1 : VariableChange K) • W)
    rw [one_smul, WeierstrassCurve.valuation_Δ_aux, dif_pos hCint',
      WeierstrassCurve.valuation_Δ_aux, dif_pos hint]
    show (IsDiscreteValuationRing.maximalIdeal R).valuation K ((C • W).Δ)
      ≤ (IsDiscreteValuationRing.maximalIdeal R).valuation K W.Δ
    exact (valuation_le_iff_vAdd_le (R := R) (Units.mk0 ((C • W).Δ) hCΔ)
      (Units.mk0 W.Δ hΔ)).2 (by omega)

/-- ★★★★★★**`v(Δ) < 12` なら Néron 指数は `0`**。 -/
theorem minimalExp_eq_zero_of_vAdd_Delta_lt (W : WeierstrassCurve K) [IsIntegral R W]
    (hΔ : W.Δ ≠ 0) (h : vAdd (tateDvrVal R K) (Units.mk0 W.Δ hΔ) < 12) :
    minimalExp R W = 0 := by
  haveI := isMinimal_of_vAdd_Delta_lt W hΔ h
  exact minimalExp_of_isMinimal W hΔ

/-! ## ★出典の紐付け(`.src`) -/

def isMinimal_of_vAdd_Delta_lt.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(v(Δ) < 12 なら極小)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
