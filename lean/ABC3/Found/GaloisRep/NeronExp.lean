import ABC3.Found.GaloisRep.TateCurveWitness

/-!
# Galois (G7) 第 312 ブロック —— **★★★★★★局所の Néron 指数**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★到達点

> 極小モデルへ移す変数変換の **`u` の付値は一意**であり(`minimal_u_vAdd_eq`)、
> それを `minimalExp R W` と書く。変数変換で **`v(u₀)` だけずれる**
> (`minimalExp_variableChange`)

★★★これが G7 の `ω_E`(Néron 微分の分数イデアル)の**各素点での指数**である。

## ★★★★★なぜ一意か——判別式が押さえている

`IsMinimal` は「整モデルの中で `v(Δ)` を最大にする」という条件(mathlib の定義)。
★2 つの極小モデル `C • W`、`C' • W` があると、`C' * C⁻¹` を `C • W` に当てて
**互いに `≤` を言い合える**ので `v(Δ)` は一致する。
★★`Δ(C • W) = u⁻¹²·Δ(W)` だから `12·v(u)` が決まり、**`v(u)` が一意**になる。

## ★★★★変数変換でのずれ方

`C • (C₀ • W)` が極小なら `(C * C₀) • W` も極小で、`u` は掛け算になる:

    minimalExp R (C₀ • W) = minimalExp R W − v(C₀.u)

★★★これが界面(第 311)の `omegaFrac_variableChange` の**局所版**である。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `vAdd_eq_of_valuation_eq` | ★★乗法的付値と `vAdd` の橋 |
| `minimal_val_Δ_eq` | ★★★★極小モデルの `v(Δ)` は一致 |
| `minimal_u_vAdd_eq` | ★★★★★**`v(u)` は一意** |
| `minimalExp`・`minimalExp_variableChange` | ★★★★★★**局所の Néron 指数** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve

variable {R : Type} [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
  {K : Type} [Field K] [Algebra R K] [IsFractionRing R K]

/-! ## ★★付値と `vAdd` の橋 -/

theorem vAdd_pow (v : Kˣ →* Multiplicative ℤ) (u : Kˣ) (n : ℕ) :
    vAdd v (u ^ n) = n * vAdd v u := by
  rw [← zpow_natCast, vAdd_zpow]

/-- ★★乗法的付値が等しければ `vAdd` も等しい。 -/
theorem vAdd_eq_of_valuation_eq (x y : Kˣ)
    (h : (IsDiscreteValuationRing.maximalIdeal R).valuation K (x : K)
       = (IsDiscreteValuationRing.maximalIdeal R).valuation K (y : K)) :
    vAdd (tateDvrVal R K) x = vAdd (tateDvrVal R K) y := by
  rw [valuation_eq_ofAdd_neg_vAdd, valuation_eq_ofAdd_neg_vAdd] at h
  have h2 : (Multiplicative.ofAdd (-(vAdd (tateDvrVal R K) x)) : Multiplicative ℤ)
      = (Multiplicative.ofAdd (-(vAdd (tateDvrVal R K) y)) : Multiplicative ℤ) :=
    WithZero.coe_injective h
  have h3 := congrArg Multiplicative.toAdd h2
  simpa using h3

/-! ## ★★★★極小モデルの `v(Δ)` は一致する -/

set_option maxHeartbeats 1600000 in
/-- ★★★★**2 つの極小モデルの `v(Δ)` は一致する**——互いに `≤` を言い合えるから。 -/
theorem minimal_val_Δ_eq (W : WeierstrassCurve K) (C C' : VariableChange K)
    (h : IsMinimal R (C • W)) (h' : IsMinimal R (C' • W)) :
    WeierstrassCurve.valuation_Δ_aux R (C • W) = WeierstrassCurve.valuation_Δ_aux R (C' • W) := by
  have hm := h.val_Δ_maximal
  have hm' := h'.val_Δ_maximal
  have hint : IsIntegral R (C' • W) := inferInstance
  have hint' : IsIntegral R (C • W) := inferInstance
  have hsm : (C' * C⁻¹) • (C • W) = C' • W := by
    rw [mul_smul, inv_smul_smul]
  have hsm' : (C * C'⁻¹) • (C' • W) = C • W := by
    rw [mul_smul, inv_smul_smul]
  rcases le_total (WeierstrassCurve.valuation_Δ_aux R (C • W))
    (WeierstrassCurve.valuation_Δ_aux R (C' • W)) with hle | hle
  · refine le_antisymm hle ?_
    have hstep := hm.2 (j := C' * C⁻¹)
      (show IsIntegral R ((C' * C⁻¹) • (C • W)) by rw [hsm]; exact hint)
      (show WeierstrassCurve.valuation_Δ_aux R ((1 : VariableChange K) • (C • W))
          ≤ WeierstrassCurve.valuation_Δ_aux R ((C' * C⁻¹) • (C • W)) by
        rw [one_smul, hsm]; exact hle)
    have hstep' : WeierstrassCurve.valuation_Δ_aux R ((C' * C⁻¹) • (C • W))
        ≤ WeierstrassCurve.valuation_Δ_aux R ((1 : VariableChange K) • (C • W)) := hstep
    rw [hsm, one_smul] at hstep'
    exact hstep'
  · refine le_antisymm ?_ hle
    have hstep := hm'.2 (j := C * C'⁻¹)
      (show IsIntegral R ((C * C'⁻¹) • (C' • W)) by rw [hsm']; exact hint')
      (show WeierstrassCurve.valuation_Δ_aux R ((1 : VariableChange K) • (C' • W))
          ≤ WeierstrassCurve.valuation_Δ_aux R ((C * C'⁻¹) • (C' • W)) by
        rw [one_smul, hsm']; exact hle)
    have hstep' : WeierstrassCurve.valuation_Δ_aux R ((C * C'⁻¹) • (C' • W))
        ≤ WeierstrassCurve.valuation_Δ_aux R ((1 : VariableChange K) • (C' • W)) := hstep
    rw [hsm', one_smul] at hstep'
    exact hstep'

set_option maxHeartbeats 1600000 in
/-- ★★★★★**極小モデルへの尺度 `u` の付値は一意**。 -/
theorem minimal_u_vAdd_eq (W : WeierstrassCurve K) (hΔ : W.Δ ≠ 0) (C C' : VariableChange K)
    (h : IsMinimal R (C • W)) (h' : IsMinimal R (C' • W)) :
    vAdd (tateDvrVal R K) C.u = vAdd (tateDvrVal R K) C'.u := by
  have hval := minimal_val_Δ_eq W C C' h h'
  have hi : IsIntegral R (C • W) := inferInstance
  have hi' : IsIntegral R (C' • W) := inferInstance
  rw [WeierstrassCurve.valuation_Δ_aux, dif_pos hi, WeierstrassCurve.valuation_Δ_aux,
    dif_pos hi'] at hval
  have hv := congrArg Subtype.val hval
  simp only at hv
  set dW : Kˣ := Units.mk0 W.Δ hΔ with hdW
  have hu1 : ((C.u⁻¹ ^ 12 * dW : Kˣ) : K) = (C • W).Δ := by
    rw [WeierstrassCurve.variableChange_Δ]
    push_cast
    rw [hdW]
    simp
  have hu2 : ((C'.u⁻¹ ^ 12 * dW : Kˣ) : K) = (C' • W).Δ := by
    rw [WeierstrassCurve.variableChange_Δ]
    push_cast
    rw [hdW]
    simp
  have hkey : vAdd (tateDvrVal R K) (C.u⁻¹ ^ 12 * dW)
      = vAdd (tateDvrVal R K) (C'.u⁻¹ ^ 12 * dW) := by
    refine vAdd_eq_of_valuation_eq _ _ ?_
    rw [hu1, hu2]
    exact hv
  rw [vAdd_mul, vAdd_mul, vAdd_pow, vAdd_pow, vAdd_inv, vAdd_inv] at hkey
  omega

/-! ## ★★★★★★局所の Néron 指数 -/

/-- ★★★★★★**局所の Néron 指数**——極小モデルへ移す尺度 `u` の付値。 -/
noncomputable def minimalExp (R : Type) [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
    {K : Type} [Field K] [Algebra R K] [IsFractionRing R K] (W : WeierstrassCurve K) : ℤ :=
  vAdd (tateDvrVal R K) (WeierstrassCurve.exists_isMinimal R W).choose.u

/-- ★★★どの極小化変換で測っても同じ。 -/
theorem minimalExp_eq (W : WeierstrassCurve K) (hΔ : W.Δ ≠ 0) (C : VariableChange K)
    (h : IsMinimal R (C • W)) : minimalExp R W = vAdd (tateDvrVal R K) C.u :=
  minimal_u_vAdd_eq W hΔ _ C (WeierstrassCurve.exists_isMinimal R W).choose_spec h

theorem variableChange_Delta_ne_zero (W : WeierstrassCurve K) (hΔ : W.Δ ≠ 0)
    (C : VariableChange K) : (C • W).Δ ≠ 0 := by
  rw [WeierstrassCurve.variableChange_Δ]
  exact mul_ne_zero (pow_ne_zero _ (Units.ne_zero _)) hΔ

/-- ★★極小な曲線の指数は `0`。 -/
theorem minimalExp_of_isMinimal (W : WeierstrassCurve K) (hΔ : W.Δ ≠ 0) [h : IsMinimal R W] :
    minimalExp R W = 0 := by
  have h1 : IsMinimal R ((1 : VariableChange K) • W) := by
    rw [one_smul]
    exact h
  rw [minimalExp_eq W hΔ 1 h1]
  show vAdd (tateDvrVal R K) 1 = 0
  rw [vAdd, map_one]
  rfl

/-- ★★★★★★**変数変換でのずれ方**——界面 `omegaFrac_variableChange` の局所版。 -/
theorem minimalExp_variableChange (W : WeierstrassCurve K) (hΔ : W.Δ ≠ 0)
    (C₀ : VariableChange K) :
    minimalExp R (C₀ • W) = minimalExp R W - vAdd (tateDvrVal R K) C₀.u := by
  obtain ⟨C, hC⟩ := WeierstrassCurve.exists_isMinimal R (C₀ • W)
  have h1 : minimalExp R (C₀ • W) = vAdd (tateDvrVal R K) C.u :=
    minimalExp_eq _ (variableChange_Delta_ne_zero W hΔ C₀) C hC
  have h2 : IsMinimal R ((C * C₀) • W) := by
    rw [mul_smul]
    exact hC
  have h3 : minimalExp R W = vAdd (tateDvrVal R K) (C * C₀).u := minimalExp_eq W hΔ _ h2
  rw [h1, h3, show (C * C₀).u = C.u * C₀.u from rfl, vAdd_mul]
  omega

/-! ## ★出典の紐付け(`.src`) -/

def minimal_u_vAdd_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(極小モデルへの尺度の付値の一意性)",
    sectionId := "genell-def-3-3" }

def minimalExp_variableChange.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(局所の Néron 指数の変数変換則)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
