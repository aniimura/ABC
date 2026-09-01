/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.MulOrder
import ABC3.Found.GaloisRep.NonDegen
import ABC3.Found.GaloisRep.PointValuation
import Mathlib.AlgebraicGeometry.EllipticCurve.DivisionPolynomial.Degree

/-!
# 第 1072 ブロック —— **素数位数の点は `preΨ_l` の根**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★★★★★これは何か —— 形式群は要らなかった

第 1067-1071 で「良い素点での捻れ点の整性には形式群が要る」と測定したが、
**在庫を引くと橋はプロジェクトに既にあった**:

    exists_divisor_root（`MulOrder.lean:127`、第 53 ブロック）
      n • P = 0  ⟹  ∃ k ∣ n, 1 ≤ k ∧ ΨSq_k(x_P) = 0

☆これは Theorem 3.8（Weil 対）のために積まれたもので、体上で完全に一般である。

★★★★これで良い素点側は**形式群を経由せずに**閉じる:

| 段 | 内容 | 出どころ |
|---|---|---|
| 1 | `l • Q = 0` ⟹ `∃ k ∣ l`, `ΨSq_k(x) = 0` | `exists_divisor_root`（在庫） |
| 2 | `k = 1` は `ΨSq_1 = 1` で潰れるので `k = l` | `ΨSq_one` |
| 3 | `l` は奇なので `ΨSq_l = preΨ_l²`、よって `preΨ_l(x) = 0` | `PSq_eval`（在庫） |
| 4 | `preΨ_l` の主係数は `l`、次数は `(l²−1)/2` | mathlib `leadingCoeff_preΨ` |
| 5 | `p ∤ l` なら `l` は単元、よって `x` は `R` 上整 | 整閉性 |
| 6 | `y` も方程式のモニック 2 次式の根なので `R` の元 | 第 1070 の `ValAtLeast` |

★本ブロックは段 1-3 である。
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve Polynomial WeierstrassCurve.Affine

universe u

variable {F : Type u} [Field F] [DecidableEq F] (W : WeierstrassCurve F)

/-- ★★★★★★★★★★★★★★★★**素数位数の点の `x` 座標は `preΨ_l` の根**（第 1072）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆`exists_divisor_root`（第 53）で `k ∣ l` に絞り、`k = 1` を `ΨSq_1 = 1` で潰し、
`l` が奇であることから `ΨSq_l = preΨ_l²` に落とす。 -/
theorem preΨ_eval_eq_zero_of_addOrderOf_prime {x y : F} (h : W.toAffine.Nonsingular x y)
    {l : ℕ} (hl : l.Prime) (hodd : l ≠ 2)
    (hQ : addOrderOf (Point.some x y h) = l) :
    (W.preΨ (l : ℤ)).eval x = 0 := by
  have hn : l • (Point.some x y h) = 0 := by
    rw [← hQ]; exact addOrderOf_nsmul_eq_zero _
  obtain ⟨k, hk1, hkdvd, hkroot⟩ := exists_divisor_root W h l hl.one_lt.le hn
  have hkl : k = l := by
    rcases (Nat.Prime.eq_one_or_self_of_dvd hl k hkdvd) with rfl | rfl
    · exfalso
      rw [show ((1 : ℕ) : ℤ) = 1 from rfl, WeierstrassCurve.ΨSq_one] at hkroot
      simp at hkroot
    · rfl
  subst hkl
  rw [PSq_eval] at hkroot
  have hoddl : Odd (k : ℤ) := by
    have : Odd k := hl.odd_of_ne_two hodd
    exact_mod_cast this
  rw [if_neg (Int.not_even_iff_odd.2 hoddl), mul_one] at hkroot
  exact pow_eq_zero_iff (two_ne_zero).elim |>.mp hkroot


/-! ## ★★★★★★★★★★★★第 1073 —— 根から整性へ -/

section Integral

open IsDedekindDomain NumberField

/-- ☆**主係数が単元なら根は整**。 -/
theorem isIntegral_of_isUnit_leadingCoeff {R S : Type*} [CommRing R] [CommRing S] [Algebra R S]
    {x : S} {q : Polynomial R} (hq : Polynomial.aeval x q = 0)
    (hu : IsUnit q.leadingCoeff) : IsIntegral R x := by
  simpa [smul_smul] using (isIntegral_leadingCoeff_smul _ _ hq).smul ((hu.unit⁻¹ : Rˣ) : R)

variable {L : Type} [Field L] [NumberField L] [DecidableEq L]

/-- ★★★★★★★★★★★★★★★★★★★★**良い素点では捕れ点の `x` 座標は整**（第 1073）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆`preΨ_l` の主係数は `l` で、`p ∤ l` ならそれは単元。
★よって `x` は `R` 上整で、`R`（付値環）は整閉なので `x ∈ R`。 -/
theorem mem_primeSubring_x_of_addOrderOf_prime (p : HeightOneSpectrum (𝓞 L))
    (W : WeierstrassCurve L) [WeierstrassCurve.IsIntegral (primeSubring p) W]
    {l : ℕ} (hl : l.Prime) (hodd : l ≠ 2) (hlu : IsUnit ((l : primeSubring p)))
    {x y : L} (h : W.toAffine.Nonsingular x y)
    (hQ : addOrderOf (WeierstrassCurve.Affine.Point.some x y h) = l) :
    x ∈ primeSubring p := by
  have hroot : (W.preΨ (l : ℤ)).eval x = 0 :=
    preΨ_eval_eq_zero_of_addOrderOf_prime W h hl hodd hQ
  have hbc : (WeierstrassCurve.integralModel (primeSubring p) W).baseChange L = W :=
    WeierstrassCurve.baseChange_integralModel_eq (primeSubring p) W
  have hmap : W.preΨ (l : ℤ)
      = ((WeierstrassCurve.integralModel (primeSubring p) W).preΨ (l : ℤ)).map
          (algebraMap (primeSubring p) L) := by
    conv_lhs => rw [← hbc]
    exact (WeierstrassCurve.integralModel (primeSubring p) W).map_preΨ
      (algebraMap (primeSubring p) L) (l : ℤ)
  have hnt : Nontrivial (primeSubring p) := inferInstance
  have hne : (((l : ℤ)) : primeSubring p) ≠ 0 := by
    have : (((l : ℤ)) : primeSubring p) = ((l : ℕ) : primeSubring p) := by push_cast; ring
    rw [this]
    exact hlu.ne_zero
  have hoddl : Odd ((l : ℤ)) := by
    have : Odd l := hl.odd_of_ne_two hodd
    exact_mod_cast this
  have hlc : ((WeierstrassCurve.integralModel (primeSubring p) W).preΨ (l : ℤ)).leadingCoeff
      = (((l : ℤ)) : primeSubring p) := by
    rw [(WeierstrassCurve.integralModel (primeSubring p) W).leadingCoeff_preΨ hne,
      if_neg (Int.not_even_iff_odd.2 hoddl)]
  have hu : IsUnit
      ((WeierstrassCurve.integralModel (primeSubring p) W).preΨ (l : ℤ)).leadingCoeff := by
    rw [hlc]
    have : (((l : ℤ)) : primeSubring p) = ((l : ℕ) : primeSubring p) := by push_cast; ring
    rw [this]; exact hlu
  have haeval : Polynomial.aeval x
      ((WeierstrassCurve.integralModel (primeSubring p) W).preΨ (l : ℤ)) = 0 := by
    rw [Polynomial.aeval_def, ← Polynomial.eval_map, ← hmap]
    exact hroot
  obtain ⟨z, hz⟩ := IsIntegrallyClosed.isIntegral_iff.1
    (isIntegral_of_isUnit_leadingCoeff haeval hu)
  exact hz ▸ z.2

/-- ★★★★★★★★★★★★**`x` が整なら `y` も整**（第 1073）。

☆`y` が深ければ左辺の付値は `2 v(y) < 0`、右辺は `≥ 0` で矛盾する。 -/
theorem mem_primeSubring_y_of_mem_x (p : HeightOneSpectrum (𝓞 L))
    (W : WeierstrassCurve L) [WeierstrassCurve.IsIntegral (primeSubring p) W]
    {x y : L} (heq : W.toAffine.Equation x y) (hx : x ∈ primeSubring p) :
    y ∈ primeSubring p := by
  by_contra hcon
  have hy : y ≠ 0 := fun h0 => hcon (by rw [h0]; exact zero_mem _)
  have hb : valAdd p (Units.mk0 y hy) < 0 := by
    by_contra hge
    rw [not_lt] at hge
    exact hcon ((mem_primeSubring_iff p y).2 ((valAdd_nonneg_iff p _).1 hge))
  obtain ⟨h1', h2', h3', h4', h6'⟩ := mem_primeSubring_of_isIntegral p W
  have h1 : ValAtLeast p 0 W.a₁ := valAtLeast_of_mem h1'
  have h2 : ValAtLeast p 0 W.a₂ := valAtLeast_of_mem h2'
  have h3 : ValAtLeast p 0 W.a₃ := valAtLeast_of_mem h3'
  have h4 : ValAtLeast p 0 W.a₄ := valAtLeast_of_mem h4'
  have h6 : ValAtLeast p 0 W.a₆ := valAtLeast_of_mem h6'
  set b : ℤ := valAdd p (Units.mk0 y hy) with hbdef
  have hxA : ValAtLeast p 0 x := valAtLeast_of_mem hx
  have hyA : ValAtLeast p b y := valAtLeast_unit p (Units.mk0 y hy)
  have hy2ne : y ^ 2 ≠ 0 := pow_ne_zero _ hy
  have hy2val : valAdd p (Units.mk0 (y ^ 2) hy2ne) = 2 * b := by
    have hmul : valAdd p (Units.mk0 y hy ^ 2) = ((2 : ℕ) : ℤ) * b := valAdd_pow p _ 2
    have hc : valAdd p (Units.mk0 (y ^ 2) hy2ne) = valAdd p (Units.mk0 y hy ^ 2) :=
      valAdd_eq_of_valuation_eq p _ _ (by simp)
    rw [hc, hmul]; norm_num
  have hrest : ValAtLeast p b (W.a₁ * x * y + W.a₃ * y) := by
    refine valAtLeast_add ?_ ?_
    · refine valAtLeast_mono ?_ (valAtLeast_mul (valAtLeast_mul h1 hxA) hyA); omega
    · refine valAtLeast_mono ?_ (valAtLeast_mul h3 hyA); omega
  have hlt : valAdd p (Units.mk0 (y ^ 2) hy2ne) < b := by omega
  have hLne : y ^ 2 + (W.a₁ * x * y + W.a₃ * y) ≠ 0 :=
    add_ne_zero_of_valAdd_lt hy2ne hrest hlt
  have hLval : valAdd p (Units.mk0 (y ^ 2 + (W.a₁ * x * y + W.a₃ * y)) hLne) = 2 * b := by
    rw [valAdd_add_eq_of_lt hy2ne hLne hrest hlt, hy2val]
  rw [WeierstrassCurve.Affine.equation_iff] at heq
  have heq' : y ^ 2 + (W.a₁ * x * y + W.a₃ * y)
      = x ^ 3 + W.a₂ * x ^ 2 + W.a₄ * x + W.a₆ := by linear_combination heq
  have hRne : x ^ 3 + W.a₂ * x ^ 2 + W.a₄ * x + W.a₆ ≠ 0 := by rw [← heq']; exact hLne
  have hx2 : ValAtLeast p 0 (x ^ 2) := by
    have hh := valAtLeast_mul hxA hxA
    have hsq : x * x = x ^ 2 := by ring
    rw [hsq] at hh
    refine valAtLeast_mono ?_ hh; omega
  have hx3 : ValAtLeast p 0 (x ^ 3) := by
    have hh := valAtLeast_mul hxA hx2
    have hcu : x * x ^ 2 = x ^ 3 := by ring
    rw [hcu] at hh
    refine valAtLeast_mono ?_ hh; omega
  have hRhs : ValAtLeast p 0 (x ^ 3 + W.a₂ * x ^ 2 + W.a₄ * x + W.a₆) := by
    refine valAtLeast_add (valAtLeast_add (valAtLeast_add ?_ ?_) ?_) ?_
    · exact hx3
    · refine valAtLeast_mono ?_ (valAtLeast_mul h2 hx2); omega
    · refine valAtLeast_mono ?_ (valAtLeast_mul h4 hxA); omega
    · exact h6
  have h0 := hRhs hRne
  rw [← valAdd_congr p hLne hRne heq', hLval] at h0
  omega

end Integral

/-! ## ★出典の紐付け(`.src`) -/

def preΨ_eval_eq_zero_of_addOrderOf_prime.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(素数位数の点の x 座標は preΨ_l の根。★無条件)",
    sectionId := "genell-lemma-3-5" }

end ABC3.Found.GaloisRep
