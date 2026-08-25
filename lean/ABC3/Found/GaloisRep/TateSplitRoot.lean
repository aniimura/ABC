import ABC3.Found.GaloisRep.TateMultRed
import Mathlib.RingTheory.Henselian

/-!
# Galois (G6) 第 305 ブロック —— **★★★★★★★★分裂条件から接線の傾きが持ち上がる**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★★★到達点

> 分裂乗法還元なら、接線の 2 次式の根が **`R` の中に取れる**(`exists_tangent_root`)

★★★mathlib の `HasSplitMultiplicativeReduction` は**剰余体で分解する**としか言わない。
それを **Hensel で `R` まで持ち上げる**のが本ブロックである。
★G6 の第 2 段(Tate 標準形への正規化)の最初の一歩にあたる。

## ★★★★★★根が単純であることが要——判別式は `−c₄c₆`

Hensel は**単純根**を要求する。接線の 2 次式

    c₄X² + a₁c₄X − (54b₆ − 3b₂b₄ + a₂c₄)

の判別式は、`b` と `c` の定義を展開すると

    (a₁c₄)² + 4c₄(54b₆ − 3b₂b₄ + a₂c₄) = −c₄c₆        (`ring` で出る)

★★★★**乗法還元では `c₄` も `c₆` も単元**なので、判別式は剰余体で `0` でない——
すなわち**根は単純**である。
★`c₆` が単元なのは `1728Δ = c₄³ − c₆²` から:`c₄³` は単元、`1728Δ ∈ 𝔪`。

## ★★★★モニックに直してから Hensel

mathlib の `HenselianRing.is_henselian` は**モニック**多項式にしか使えない。
★`c₄` が単元なので `c₄⁻¹` を掛けて `X² + a₁X − c₄⁻¹Q` に直し、
得られた根に `c₄` を掛けて戻す。★★`IsAdicComplete.henselianRing` がそのまま使えた。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `integralModel_c4_isUnit`・`integralModel_c6_isUnit` | ★★★整モデルの `c₄, c₆` は単元 |
| `integralModel_Delta_mem` | ★★★整モデルの `Δ` は `𝔪` の元 |
| `tangent_disc_eq` | ★★★★★判別式は `−c₄c₆` |
| `exists_tangent_root` | ★★★★★★★★**接線の傾きは `R` に持ち上がる** |
-/

namespace ABC3.Found.GaloisRep

open IsDedekindDomain WeierstrassCurve Polynomial

section Units

variable {R : Type} [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
  {K : Type} [Field K] [Algebra R K] [IsFractionRing R K]

/-- ★★★乗法還元なら整モデルの `c₄` は単元。 -/
theorem integralModel_c4_isUnit (W : WeierstrassCurve K)
    [hmul : W.HasMultiplicativeReduction R] :
    IsUnit (integralModel R W).c₄ := by
  obtain ⟨r, hru, hr⟩ := exists_isUnit_of_valuation_eq_one W.c₄ hmul.multiplicativeReduction
  have hEq : r = (integralModel R W).c₄ := by
    apply IsFractionRing.injective R K
    rw [hr, WeierstrassCurve.integralModel_c₄_eq]
  rwa [hEq] at hru

/-- ★★★乗法還元なら整モデルの `Δ` は `𝔪` の元。 -/
theorem integralModel_Delta_mem (W : WeierstrassCurve K)
    [hmul : W.HasMultiplicativeReduction R] (hΔ0 : W.Δ ≠ 0) :
    (integralModel R W).Δ ∈ IsLocalRing.maximalIdeal R := by
  obtain ⟨d, hdm, hd⟩ := exists_mem_of_valuation_lt_one W.Δ hΔ0 hmul.badReduction
  have hEq : d = (integralModel R W).Δ := by
    apply IsFractionRing.injective R K
    rw [hd, WeierstrassCurve.integralModel_Δ_eq]
  rwa [hEq] at hdm

/-- ★★★乗法還元なら整モデルの `c₆` も単元(`1728Δ = c₄³ − c₆²` から)。 -/
theorem integralModel_c6_isUnit (W : WeierstrassCurve K)
    [hmul : W.HasMultiplicativeReduction R] (hΔ0 : W.Δ ≠ 0) :
    IsUnit (integralModel R W).c₆ := by
  have h4 : IsUnit (integralModel R W).c₄ := integralModel_c4_isUnit W
  have hΔ : (integralModel R W).Δ ∈ IsLocalRing.maximalIdeal R := integralModel_Delta_mem W hΔ0
  have hrel : (integralModel R W).c₆ ^ 2
      = (integralModel R W).c₄ ^ 3 - 1728 * (integralModel R W).Δ := by
    have h := (integralModel R W).c_relation
    linear_combination h
  have hsq : IsUnit ((integralModel R W).c₆ ^ 2) := by
    rw [hrel, ← IsLocalRing.notMem_maximalIdeal]
    intro hc
    have h3 : (integralModel R W).c₄ ^ 3 ∈ IsLocalRing.maximalIdeal R := by
      have heq : (integralModel R W).c₄ ^ 3
          = ((integralModel R W).c₄ ^ 3 - 1728 * (integralModel R W).Δ)
            + 1728 * (integralModel R W).Δ := by ring
      rw [heq]
      exact Ideal.add_mem _ hc (Ideal.mul_mem_left _ _ hΔ)
    exact (IsLocalRing.notMem_maximalIdeal.2 (h4.pow 3)) h3
  exact (isUnit_pow_iff (by norm_num)).1 hsq

end Units

/-! ## ★★★★★判別式は `−c₄c₆` -/

/-- ★★★★★**接線の 2 次式の判別式は `−c₄c₆`**——`b`・`c` を展開すれば `ring` で出る。 -/
theorem tangent_disc_eq {R : Type} [CommRing R] (E : WeierstrassCurve R) :
    (E.a₁ * E.c₄) ^ 2 + 4 * E.c₄ * (54 * E.b₆ - 3 * E.b₂ * E.b₄ + E.a₂ * E.c₄)
      = -(E.c₄ * E.c₆) := by
  simp only [WeierstrassCurve.c₄, WeierstrassCurve.c₆, WeierstrassCurve.b₂, WeierstrassCurve.b₄,
    WeierstrassCurve.b₆]
  ring

/-! ## ★★★★★★★★接線の傾きは `R` に持ち上がる -/

section Lift

variable {R : Type} [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
  [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
  {K : Type} [Field K] [Algebra R K] [IsFractionRing R K]

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★★**分裂乗法還元なら接線の 2 次式は `R` に根をもつ**。

★剰余体での分解(mathlib の定義)を Hensel で `R` まで持ち上げる。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem exists_tangent_root (W : WeierstrassCurve K) [hmul : W.HasMultiplicativeReduction R]
    (hsplit : W.HasSplitMultiplicativeReduction R) (hΔ0 : W.Δ ≠ 0) :
    ∃ α : R, (integralModel R W).c₄ * α ^ 2
        + (integralModel R W).a₁ * (integralModel R W).c₄ * α
        - (54 * (integralModel R W).b₆ - 3 * (integralModel R W).b₂ * (integralModel R W).b₄
          + (integralModel R W).a₂ * (integralModel R W).c₄) = 0 := by
  have h4 : IsUnit (integralModel R W).c₄ := integralModel_c4_isUnit W
  have h6 : IsUnit (integralModel R W).c₆ := integralModel_c6_isUnit W hΔ0
  have hsp := hsplit.splitMultiplicativeReduction
  set E := integralModel R W with hE
  set Q : R := 54 * E.b₆ - 3 * E.b₂ * E.b₄ + E.a₂ * E.c₄ with hQ
  set φ := algebraMap R (IsLocalRing.ResidueField R) with hφ
  have hmapped : (Polynomial.map φ
      (Polynomial.C E.c₄ * Polynomial.X ^ 2 + Polynomial.C (E.a₁ * E.c₄) * Polynomial.X
        - Polynomial.C Q))
      = Polynomial.C (φ E.c₄) * Polynomial.X ^ 2 + Polynomial.C (φ (E.a₁ * E.c₄)) * Polynomial.X
        + Polynomial.C (φ (-Q)) := by
    simp only [Polynomial.map_add, Polynomial.map_sub, Polynomial.map_mul, Polynomial.map_pow,
      Polynomial.map_C, Polynomial.map_X, map_neg]
    ring
  have hc4ne : φ E.c₄ ≠ 0 := (h4.map φ).ne_zero
  have hc6ne : φ E.c₆ ≠ 0 := (h6.map φ).ne_zero
  have hdeg : (Polynomial.map φ
      (Polynomial.C E.c₄ * Polynomial.X ^ 2 + Polynomial.C (E.a₁ * E.c₄) * Polynomial.X
        - Polynomial.C Q)).degree ≠ 0 := by
    rw [hmapped, Polynomial.degree_quadratic hc4ne]
    norm_num
  obtain ⟨b, hb⟩ := Polynomial.Splits.exists_eval_eq_zero hsp hdeg
  rw [hmapped] at hb
  simp only [Polynomial.eval_add, Polynomial.eval_mul, Polynomial.eval_pow, Polynomial.eval_C,
    Polynomial.eval_X, Polynomial.eval_neg, map_neg, map_mul] at hb
  have hident : (E.a₁ * E.c₄) ^ 2 + 4 * E.c₄ * Q = -(E.c₄ * E.c₆) := tangent_disc_eq E
  have hidentk : (φ E.a₁ * φ E.c₄) ^ 2 + 4 * φ E.c₄ * φ Q = -(φ E.c₄ * φ E.c₆) := by
    have h := congrArg φ hident
    simp only [map_add, map_mul, map_pow, map_neg, map_ofNat] at h
    exact h
  have hder : 2 * b + φ E.a₁ ≠ 0 := by
    intro hc
    have hkey : (φ E.c₄ * (2 * b + φ E.a₁)) ^ 2 = -(φ E.c₄ * φ E.c₆) := by
      linear_combination 4 * φ E.c₄ * hb + hidentk
    rw [hc, mul_zero] at hkey
    have hz : φ E.c₄ * φ E.c₆ = 0 := by linear_combination hkey
    rcases mul_eq_zero.1 hz with h | h
    · exact hc4ne h
    · exact hc6ne h
  set v : R := ((h4.unit⁻¹ : Rˣ) : R) with hv
  have hvc : v * E.c₄ = 1 := by
    have h2 : ((h4.unit⁻¹ : Rˣ) : R) * ((h4.unit : Rˣ) : R) = 1 := by
      rw [← Units.val_mul, inv_mul_cancel]
      rfl
    rw [h4.unit_spec] at h2
    exact h2
  have hvck : φ v * φ E.c₄ = 1 := by
    rw [← map_mul, hvc, map_one]
  set f : Polynomial R := Polynomial.X ^ 2 + Polynomial.C E.a₁ * Polynomial.X
    - Polynomial.C (v * Q) with hf
  have hmonic : f.Monic := by
    rw [hf]
    monicity!
  obtain ⟨a₀, ha₀⟩ := IsLocalRing.residue_surjective (R := R) b
  have ha₀' : φ a₀ = b := ha₀
  have heval : Polynomial.eval a₀ f ∈ IsLocalRing.maximalIdeal R := by
    rw [← Ideal.Quotient.eq_zero_iff_mem]
    show φ (Polynomial.eval a₀ f) = 0
    rw [hf]
    simp only [Polynomial.eval_add, Polynomial.eval_sub, Polynomial.eval_mul, Polynomial.eval_pow,
      Polynomial.eval_C, Polynomial.eval_X, map_add, map_sub, map_mul, map_pow, ha₀']
    linear_combination φ v * hb - (b ^ 2 + φ E.a₁ * b) * hvck
  have hderiv : Polynomial.derivative f = 2 * Polynomial.X + Polynomial.C E.a₁ := by
    rw [hf]
    simp [Polynomial.derivative_sub, Polynomial.derivative_add, Polynomial.derivative_mul]
    ring
  have hunit : IsUnit ((Ideal.Quotient.mk (IsLocalRing.maximalIdeal R))
      (Polynomial.eval a₀ (Polynomial.derivative f))) := by
    rw [hderiv]
    simp only [Polynomial.eval_add, Polynomial.eval_mul, Polynomial.eval_C, Polynomial.eval_X,
      Polynomial.eval_ofNat]
    show IsUnit (φ (2 * a₀ + E.a₁))
    rw [map_add, map_mul, ha₀']
    refine Ne.isUnit ?_
    simpa [map_ofNat] using hder
  obtain ⟨α, hα, _⟩ := (IsAdicComplete.henselianRing R (IsLocalRing.maximalIdeal R)).is_henselian
    f hmonic a₀ heval hunit
  refine ⟨α, ?_⟩
  have hαr : α ^ 2 + E.a₁ * α - v * Q = 0 := by
    have h := hα
    rw [Polynomial.IsRoot, hf] at h
    simpa using h
  linear_combination E.c₄ * hαr + Q * hvc

end Lift

/-! ## ★出典の紐付け(`.src`) -/

def tangent_disc_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(接線の 2 次式の判別式)",
    sectionId := "genell-def-3-3" }

def exists_tangent_root.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(分裂条件から接線の傾きを持ち上げる)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
