/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.MulOrder
import ABC3.Found.GaloisRep.NonDegen
import ABC3.Found.GaloisRep.PointValuation
import ABC3.Found.GenEll.VeluPointSet
import ABC3.Found.GenEll.SymmSum
import ABC3.Found.GaloisRep.TateVeluPoints
import ABC3.Found.GenEll.VeluImage
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

open IsDedekindDomain NumberField ABC3.Found.GenEll
open scoped Classical

/-- ☆**主係数が単元なら根は整**。 -/
theorem isIntegral_of_isUnit_leadingCoeff {R S : Type*} [CommRing R] [CommRing S] [Algebra R S]
    {x : S} {q : Polynomial R} (hq : Polynomial.aeval x q = 0)
    (hu : IsUnit q.leadingCoeff) : IsIntegral R x := by
  simpa [smul_smul] using (isIntegral_leadingCoeff_smul _ _ hq).smul ((hu.unit⁻¹ : Rˣ) : R)

variable {L : Type} [Field L] [NumberField L]

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


/-- ☆**`v_p(n) = 0` なら `n` は `primeSubring p` で単元**。 -/
theorem isUnit_natCast_primeSubring (p : HeightOneSpectrum (𝓞 L)) {n : ℕ}
    (hn : ((n : L)) ≠ 0) (hval : valAdd p (Units.mk0 ((n : L)) hn) = 0) :
    IsUnit ((n : primeSubring p)) := by
  have hinv : valAdd p (Units.mk0 ((n : L)) hn)⁻¹ = 0 := by
    rw [valAdd_inv, hval, neg_zero]
  have hmem : ((n : L))⁻¹ ∈ primeSubring p := by
    rw [mem_primeSubring_iff]
    have h := (valAdd_nonneg_iff p (Units.mk0 ((n : L)) hn)⁻¹).1 (le_of_eq hinv.symm)
    simpa using h
  refine ⟨⟨(n : primeSubring p), ⟨((n : L))⁻¹, hmem⟩, ?_, ?_⟩, rfl⟩ <;>
    · apply Subtype.ext
      push_cast
      field_simp

/-- ★★★★★★★★★★★★★★★★**`p ∣ l` でも `l·x` は整**（第 1076）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆`preΨ_l` の主係数は `l` なので、単元でなくても
`isIntegral_leadingCoeff_smul` が `l • x` の整性を与える。
★これが `p ∣ l` での局所評価 `v_p(x) ≥ −v_p(l)` の正体である。 -/
theorem natCast_mul_x_mem_of_addOrderOf_prime (p : HeightOneSpectrum (𝓞 L))
    (E : WeierstrassCurve L) [WeierstrassCurve.IsIntegral (primeSubring p) E]
    {l : ℕ} (hl : l.Prime) (hodd : l ≠ 2)
    {x y : L} (h : E.toAffine.Nonsingular x y)
    (hQ : addOrderOf (WeierstrassCurve.Affine.Point.some x y h) = l) :
    ((l : L) * x) ∈ primeSubring p := by
  have hroot : (E.preΨ (l : ℤ)).eval x = 0 :=
    preΨ_eval_eq_zero_of_addOrderOf_prime E h hl hodd hQ
  have hbc : (WeierstrassCurve.integralModel (primeSubring p) E).baseChange L = E :=
    WeierstrassCurve.baseChange_integralModel_eq (primeSubring p) E
  have hmap : E.preΨ (l : ℤ)
      = ((WeierstrassCurve.integralModel (primeSubring p) E).preΨ (l : ℤ)).map
          (algebraMap (primeSubring p) L) := by
    conv_lhs => rw [← hbc]
    exact (WeierstrassCurve.integralModel (primeSubring p) E).map_preΨ
      (algebraMap (primeSubring p) L) (l : ℤ)
  have haeval : Polynomial.aeval x
      ((WeierstrassCurve.integralModel (primeSubring p) E).preΨ (l : ℤ)) = 0 := by
    rw [Polynomial.aeval_def, ← Polynomial.eval_map, ← hmap]
    exact hroot
  have hoddl : Odd ((l : ℤ)) := by
    have : Odd l := hl.odd_of_ne_two hodd
    exact_mod_cast this
  have hne : (((l : ℤ)) : primeSubring p) ≠ 0 := by
    have hc : (((l : ℤ)) : primeSubring p) = ((l : ℕ) : primeSubring p) := by push_cast; ring
    rw [hc]
    simpa using (Nat.cast_ne_zero (R := primeSubring p)).2 hl.ne_zero
  have hlc : ((WeierstrassCurve.integralModel (primeSubring p) E).preΨ (l : ℤ)).leadingCoeff
      = (((l : ℤ)) : primeSubring p) := by
    rw [(WeierstrassCurve.integralModel (primeSubring p) E).leadingCoeff_preΨ hne,
      if_neg (Int.not_even_iff_odd.2 hoddl)]
  have hint : IsIntegral (primeSubring p)
      (((WeierstrassCurve.integralModel (primeSubring p) E).preΨ (l : ℤ)).leadingCoeff • x) :=
    isIntegral_leadingCoeff_smul _ _ haeval
  obtain ⟨z, hz⟩ := IsIntegrallyClosed.isIntegral_iff.1 hint
  have hsm : ((WeierstrassCurve.integralModel (primeSubring p) E).preΨ (l : ℤ)).leadingCoeff • x
      = (l : L) * x := by
    rw [hlc, Algebra.smul_def]
    congr 1
    push_cast
    ring
  rw [hsm] at hz
  exact hz ▸ z.2

/-! ## ★★★★★★★★★★★★★★★★第 1074 —— すべての倍点の座標が整 -/

/-- ★★★★★★★★★★★★★★★★**位数 `l` の点のすべての倍点の座標は整**（第 1074）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆`l` が素数なので `0 < k < l` なら `k • Q` も位数 `l` である。 -/
theorem pointCoords_mem_of_addOrderOf_prime (p : HeightOneSpectrum (𝓞 L))
    (E : WeierstrassCurve L) [WeierstrassCurve.IsIntegral (primeSubring p) E]
    {l : ℕ} (hl : l.Prime) (hodd : l ≠ 2) (hlu : IsUnit ((l : primeSubring p)))
    (Q : E.toAffine.Point) (hQ : addOrderOf Q = l)
    {k : ℕ} (hk0 : k ≠ 0) (hkl : k < l) :
    (pointCoords (k • Q)).1 ∈ primeSubring p ∧ (pointCoords (k • Q)).2 ∈ primeSubring p := by
  have hkne : k • Q ≠ 0 := nsmul_ne_zero_of_lt_addOrderOf hQ hk0 hkl
  have hlz : l • Q = 0 := by rw [← hQ]; exact addOrderOf_nsmul_eq_zero Q
  have hdvd : addOrderOf (k • Q) ∣ l := by
    refine addOrderOf_dvd_of_nsmul_eq_zero ?_
    rw [smul_comm, hlz, smul_zero]
  have hord : addOrderOf (k • Q) = l := by
    rcases (Nat.Prime.eq_one_or_self_of_dvd hl _ hdvd) with h1 | h1
    · exact absurd (AddMonoid.addOrderOf_eq_one_iff.1 h1) hkne
    · exact h1
  rcases hkQ : k • Q with _ | ⟨x, y, h⟩
  · exact absurd hkQ hkne
  · have hord' : addOrderOf (WeierstrassCurve.Affine.Point.some x y h) = l := hkQ ▸ hord
    have hx := mem_primeSubring_x_of_addOrderOf_prime p E hl hodd hlu h hord'
    refine ⟨?_, ?_⟩
    · simpa only [hkQ, pointCoords_some] using hx
    · simpa only [hkQ, pointCoords_some] using mem_primeSubring_y_of_mem_x p E h.1 hx


open Finset in
/-- ★★★★★★★★★★★★★★★★★★★★**[GenEll] 良い素点では Vélu の商は整**（第 1074）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆座標が `R` の元（第 1073）、`v` は和そのもの、`w` は `exists_veluW_of_inv`（第 960）。
★形式群は一度も使っていない。 -/
theorem isIntegral_veluQuotientFull_of_addOrderOf_prime (p : HeightOneSpectrum (𝓞 L))
    (E : WeierstrassCurve L) [hE : WeierstrassCurve.IsIntegral (primeSubring p) E]
    {l : ℕ} (hl : l.Prime) (hodd : l ≠ 2) (hlu : IsUnit ((l : primeSubring p)))
    (Q : E.toAffine.Point) (hQ : addOrderOf Q = l) :
    (veluQuotientFull E (((range l).erase 0).image
      (fun k : ℕ => pointCoords (k • Q)))).IsIntegral (primeSubring p) := by
  classical
  obtain ⟨Wi, hWi⟩ := hE.integral
  obtain ⟨m, rfl⟩ : ∃ m, l = 2 * m + 1 := hl.odd_of_ne_two hodd
  have hlz : (2 * m + 1) • Q = 0 := by rw [← hQ]; exact addOrderOf_nsmul_eq_zero Q
  have hmem : ∀ k ∈ (range (2 * m + 1)).erase 0,
      (pointCoords (k • Q)).1 ∈ primeSubring p ∧
        (pointCoords (k • Q)).2 ∈ primeSubring p := by
    intro k hk
    rw [mem_erase, mem_range] at hk
    exact pointCoords_mem_of_addOrderOf_prime p E hl hodd hlu Q hQ hk.1 hk.2
  set X : ℕ → primeSubring p := fun i =>
    if h : (pointCoords (i • Q)).1 ∈ primeSubring p then ⟨(pointCoords (i • Q)).1, h⟩
    else 0 with hXdef
  set Y : ℕ → primeSubring p := fun i =>
    if h : (pointCoords (i • Q)).2 ∈ primeSubring p then ⟨(pointCoords (i • Q)).2, h⟩
    else 0 with hYdef
  have hXc : ∀ i ∈ (range (2 * m + 1)).erase 0,
      algebraMap (primeSubring p) L (X i) = (pointCoords (i • Q)).1 := by
    intro i hi
    simp only [hXdef, dif_pos (hmem i hi).1]
    rfl
  have hYc : ∀ i ∈ (range (2 * m + 1)).erase 0,
      algebraMap (primeSubring p) L (Y i) = (pointCoords (i • Q)).2 := by
    intro i hi
    simp only [hYdef, dif_pos (hmem i hi).2]
    rfl
  have hP : ∀ i ∈ (range (2 * m + 1)).erase 0, pointCoords (i • Q)
      = ((algebraMap (primeSubring p) L (X i),
          algebraMap (primeSubring p) L (Y i)) : L × L) := by
    intro i hi
    rw [hXc i hi, hYc i hi]
  -- ☆添字の反転は点の反転
  have hsub : ∀ i ∈ Icc 1 m, (2 * m + 1 - i) ∈ (range (2 * m + 1)).erase 0 := by
    intro i hi
    rw [mem_Icc] at hi
    rw [mem_erase, mem_range]
    omega
  have hin : ∀ i ∈ Icc 1 m, i ∈ (range (2 * m + 1)).erase 0 := by
    intro i hi
    rw [mem_Icc] at hi
    rw [mem_erase, mem_range]
    omega
  have hneg : ∀ i ∈ Icc 1 m,
      pointCoords ((2 * m + 1 - i) • Q)
        = ((pointCoords (i • Q)).1,
           (Wi.map (algebraMap (primeSubring p) L)).toAffine.negY
             (pointCoords (i • Q)).1 (pointCoords (i • Q)).2) := by
    intro i hi
    rw [mem_Icc] at hi
    have hkne : i • Q ≠ 0 :=
      nsmul_ne_zero_of_lt_addOrderOf hQ (by omega) (by omega)
    have := nsmul_eq_neg_nsmul_of_addOrderOf hlz (by omega : i ≤ 2 * m + 1)
    rw [this]
    have hEq : E = Wi.map (algebraMap (primeSubring p) L) := hWi
    subst hEq
    exact pointCoords_neg hkne
  have hXinv : ∀ i ∈ Icc 1 m, X (2 * m + 1 - i) = X i := by
    intro i hi
    have h1 := hXc _ (hsub i hi)
    have h2 := hXc i (hin i hi)
    apply Subtype.ext
    have : algebraMap (primeSubring p) L (X (2 * m + 1 - i))
        = algebraMap (primeSubring p) L (X i) := by
      rw [h1, h2, hneg i hi]
    exact this
  have hYinv : ∀ i ∈ Icc 1 m, Y (2 * m + 1 - i) = Wi.toAffine.negY (X i) (Y i) := by
    intro i hi
    apply Subtype.ext
    have hL : algebraMap (primeSubring p) L (Y (2 * m + 1 - i))
        = algebraMap (primeSubring p) L (Wi.toAffine.negY (X i) (Y i)) := by
      rw [hYc _ (hsub i hi), hneg i hi]
      rw [← hXc i (hin i hi), ← hYc i (hin i hi)]
      exact (WeierstrassCurve.Affine.map_negY (W' := Wi)
        (algebraMap (primeSubring p) L) (X i) (Y i)).symm
    exact hL
  obtain ⟨w, hw⟩ := ABC3.Found.GenEll.exists_veluW_of_inv Wi m X Y hXinv hYinv
  -- ★単射性
  have hinj : ∀ i ∈ (range (2 * m + 1)).erase 0, ∀ j ∈ (range (2 * m + 1)).erase 0,
      ((algebraMap (primeSubring p) L (X i), algebraMap (primeSubring p) L (Y i)) : L × L)
        = ((algebraMap (primeSubring p) L (X j), algebraMap (primeSubring p) L (Y j)) : L × L)
      → i = j := by
    intro i hi j hj hij
    rw [mem_erase, mem_range] at hi hj
    have hne_i : i • Q ≠ 0 := nsmul_ne_zero_of_lt_addOrderOf hQ hi.1 hi.2
    have hne_j : j • Q ≠ 0 := nsmul_ne_zero_of_lt_addOrderOf hQ hj.1 hj.2
    have hEq : i • Q = j • Q := by
      refine pointCoords_injective hne_i hne_j ?_
      rw [hP i (by rw [mem_erase, mem_range]; exact hi),
        hP j (by rw [mem_erase, mem_range]; exact hj)]
      exact hij
    rcases le_total i j with hle | hle
    · have h1 : (j - i) • Q + i • Q = j • Q := by
        rw [← add_nsmul, Nat.sub_add_cancel hle]
      rw [← hEq] at h1
      have h2 : (j - i) • Q = 0 :=
        add_right_cancel (b := i • Q) (by rw [h1, zero_add])
      have h3 : addOrderOf Q ∣ (j - i) := addOrderOf_dvd_of_nsmul_eq_zero h2
      rw [hQ] at h3
      have := Nat.eq_zero_of_dvd_of_lt h3
      omega
    · have h1 : (i - j) • Q + j • Q = i • Q := by
        rw [← add_nsmul, Nat.sub_add_cancel hle]
      rw [hEq] at h1
      have h2 : (i - j) • Q = 0 :=
        add_right_cancel (b := j • Q) (by rw [h1, zero_add])
      have h3 : addOrderOf Q ∣ (i - j) := addOrderOf_dvd_of_nsmul_eq_zero h2
      rw [hQ] at h3
      have := Nat.eq_zero_of_dvd_of_lt h3
      omega
  -- ★組み立て
  have hS : ((range (2 * m + 1)).erase 0).image (fun k : ℕ => pointCoords (k • Q))
      = ((range (2 * m + 1)).erase 0).image
          (fun i : ℕ => ((algebraMap (primeSubring p) L (X i),
                          algebraMap (primeSubring p) L (Y i)) : L × L)) :=
    Finset.image_congr (fun i hi => hP i hi)
  refine ⟨ABC3.Found.GenEll.veluCurve Wi
    (∑ i ∈ (range (2 * m + 1)).erase 0, ABC3.Found.GenEll.veluV2 Wi (X i) (Y i)) w, ?_⟩
  rw [hS, hWi]
  exact ABC3.Found.GenEll.veluQuotientFull_image_eq Wi ((range (2 * m + 1)).erase 0) X Y hinj
    _ w (two_ne_zero) rfl hw

end Integral

/-! ## ★出典の紐付け(`.src`) -/

def preΨ_eval_eq_zero_of_addOrderOf_prime.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(素数位数の点の x 座標は preΨ_l の根。★無条件)",
    sectionId := "genell-lemma-3-5" }

end ABC3.Found.GaloisRep
