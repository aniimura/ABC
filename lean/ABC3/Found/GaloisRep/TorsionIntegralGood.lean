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

/-! ## ★★★★★★★★★★★★第 1077 —— `p ∣ l` での深さの上限 -/

/-- ★★★★★★★★★★★★★★★★**捕れ点の深さは `v_p(l)/2` 以下**（第 1077）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆第 1076 で `l·x` が整だから `v(x) ≥ −v_p(l)`、
第 1070 で `v(x) = −2m`・`v(y) = −3m` なので `2m ≤ v_p(l)`。
★よって `m ≤ M` で `v(x) ≥ −2M`・`v(y) ≥ −3M`。 -/
theorem valAtLeast_pointCoords_of_le (p : HeightOneSpectrum (𝓞 L))
    (E : WeierstrassCurve L) [WeierstrassCurve.IsIntegral (primeSubring p) E]
    {l : ℕ} (hl : l.Prime) (hodd : l ≠ 2)
    (Q : E.toAffine.Point) (hQ : addOrderOf Q = l)
    (M : ℤ) (hM0 : 0 ≤ M)
    (hlne : ((l : L)) ≠ 0) (hM : valAdd p (Units.mk0 ((l : L)) hlne) ≤ 2 * M)
    {k : ℕ} (hk0 : k ≠ 0) (hkl : k < l) :
    ValAtLeast p (-2 * M) (pointCoords (k • Q)).1 ∧
      ValAtLeast p (-3 * M) (pointCoords (k • Q)).2 := by
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
    -- ☆`x` が整なら両方整
    by_cases hxint : x ∈ primeSubring p
    · have hyint := mem_primeSubring_y_of_mem_x p E h.1 hxint
      refine ⟨?_, ?_⟩ <;> simp only [hkQ, pointCoords_some]
      · exact valAtLeast_mono (by omega) (valAtLeast_of_mem hxint)
      · exact valAtLeast_mono (by omega) (valAtLeast_of_mem hyint)
    -- ★深い場合
    · have hx0 : x ≠ 0 := fun hc => hxint (by rw [hc]; exact zero_mem _)
      have hneg : valAdd p (Units.mk0 x hx0) < 0 := by
        by_contra hge
        rw [not_lt] at hge
        exact hxint ((mem_primeSubring_iff p x).2 ((valAdd_nonneg_iff p _).1 hge))
      have hy0 : y ≠ 0 := y_ne_zero_of_valAdd_x_neg p E hx0 h.1 hneg
      obtain ⟨m, hm0, hmx, hmy⟩ := exists_depth_of_valAdd_x_neg p E hx0 hy0 h.1 hneg
      -- ☆`l·x` は整
      have hlx := natCast_mul_x_mem_of_addOrderOf_prime p E hl hodd h hord'
      have hlxne : ((l : L)) * x ≠ 0 := mul_ne_zero hlne hx0
      have hlxv : 0 ≤ valAdd p (Units.mk0 (((l : L)) * x) hlxne) :=
        (valAdd_nonneg_iff p _).2 ((mem_primeSubring_iff p _).1 hlx)
      have hsplit : valAdd p (Units.mk0 (((l : L)) * x) hlxne)
          = valAdd p (Units.mk0 ((l : L)) hlne) + valAdd p (Units.mk0 x hx0) := by
        rw [← valAdd_mul p (Units.mk0 ((l : L)) hlne) (Units.mk0 x hx0)]
        exact valAdd_eq_of_valuation_eq p _ _ (by simp)
      rw [hsplit] at hlxv
      have hmle : (m : ℤ) ≤ M := by omega
      refine ⟨?_, ?_⟩ <;> simp only [hkQ, pointCoords_some]
      · refine valAtLeast_mono ?_ (valAtLeast_unit p (Units.mk0 x hx0)); omega
      · refine valAtLeast_mono ?_ (valAtLeast_unit p (Units.mk0 y hy0)); omega

/-! ## ★★★★★★★★★★★★第 1078 —— Vélu の係数の付値評価 -/

/-- ☆有限和も `ValAtLeast` を保つ。 -/
theorem valAtLeast_sum {ι : Type*} (p : HeightOneSpectrum (𝓞 L)) {s : Finset ι}
    {f : ι → L} {n : ℤ} (h : ∀ i ∈ s, ValAtLeast p n (f i)) :
    ValAtLeast p n (∑ i ∈ s, f i) :=
  Finset.sum_induction f (ValAtLeast p n) (fun _ _ ha hb => valAtLeast_add ha hb)
    (valAtLeast_zero p n) h

/-- ★★★★★★★★★★★★★★★★**Vélu の `v`・`w` の付値評価**（第 1078）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆`v(x) ≥ −2M`・`v(y) ≥ −3M` なら `v(v) ≥ −4M`・`v(w) ≥ −6M`。
★`w` の `/2` は `v_p(2) = 0` のとき無害である（`l` が奇なら `p ∣ l` でそう）。 -/
theorem valAtLeast_veluVFull_veluWFull (p : HeightOneSpectrum (𝓞 L))
    (E : WeierstrassCurve L) [WeierstrassCurve.IsIntegral (primeSubring p) E]
    (S : Finset (L × L)) (M : ℤ) (hM0 : 0 ≤ M)
    (hS : ∀ z ∈ S, ValAtLeast p (-2 * M) z.1 ∧ ValAtLeast p (-3 * M) z.2)
    (h2 : ValAtLeast p 0 ((2 : L)⁻¹)) :
    ValAtLeast p (-4 * M) (veluVFull E S) ∧ ValAtLeast p (-6 * M) (veluWFull E S) := by
  obtain ⟨ha1, ha2, ha3, ha4, _⟩ := mem_primeSubring_of_isIntegral p E
  have h1 : ValAtLeast p 0 E.a₁ := valAtLeast_of_mem ha1
  have h2' : ValAtLeast p 0 E.a₂ := valAtLeast_of_mem ha2
  have h3 : ValAtLeast p 0 E.a₃ := valAtLeast_of_mem ha3
  have h4 : ValAtLeast p 0 E.a₄ := valAtLeast_of_mem ha4
  have hc3 : ValAtLeast p 0 ((3 : L)) := by
    refine valAtLeast_of_mem ?_
    simpa using (3 : primeSubring p).2
  have hc2 : ValAtLeast p 0 ((2 : L)) := by
    refine valAtLeast_of_mem ?_
    simpa using (2 : primeSubring p).2
  -- ☆各点での `veluV2`・`veluGy`・`veluU`
  have hV2 : ∀ z ∈ S, ValAtLeast p (-4 * M) (veluV2 E z.1 z.2) := by
    intro z hz
    obtain ⟨hx, hy⟩ := hS z hz
    have hx2 : ValAtLeast p (-4 * M) (z.1 ^ 2) := by
      have hh := valAtLeast_mul hx hx
      have hsq : z.1 * z.1 = z.1 ^ 2 := by ring
      rw [hsq] at hh
      refine valAtLeast_mono ?_ hh; omega
    show ValAtLeast p (-4 * M) (veluGx E z.1 z.2)
    rw [veluGx, sub_eq_add_neg]
    refine valAtLeast_add (valAtLeast_add (valAtLeast_add ?_ ?_) ?_) ?_
    · refine valAtLeast_mono ?_ (valAtLeast_mul hc3 hx2); omega
    · refine valAtLeast_mono ?_ (valAtLeast_mul (valAtLeast_mul hc2 h2') hx); omega
    · refine valAtLeast_mono ?_ h4; omega
    · refine valAtLeast_mono ?_ (valAtLeast_neg (valAtLeast_mul h1 hy)); omega
  have hU : ∀ z ∈ S, ValAtLeast p (-6 * M) (veluU E z.1 z.2) := by
    intro z hz
    obtain ⟨hx, hy⟩ := hS z hz
    have hgy : ValAtLeast p (-3 * M) (veluGy E z.1 z.2) := by
      rw [veluGy, sub_eq_add_neg, sub_eq_add_neg]
      refine valAtLeast_add (valAtLeast_add ?_ ?_) ?_
      · refine valAtLeast_mono ?_ (valAtLeast_mul (valAtLeast_neg hc2) hy); omega
      · refine valAtLeast_mono ?_ (valAtLeast_neg (valAtLeast_mul h1 hx)); omega
      · refine valAtLeast_mono ?_ (valAtLeast_neg h3); omega
    rw [veluU]
    have hh := valAtLeast_mul hgy hgy
    have hsq : veluGy E z.1 z.2 * veluGy E z.1 z.2 = veluGy E z.1 z.2 ^ 2 := by ring
    rw [hsq] at hh
    refine valAtLeast_mono ?_ hh; omega
  refine ⟨?_, ?_⟩
  · rw [veluVFull]
    refine valAtLeast_sum p (fun z hz => ?_)
    refine valAtLeast_mono ?_ (hV2 z hz); omega
  · rw [veluWFull]
    refine valAtLeast_sum p (fun z hz => ?_)
    refine valAtLeast_add ?_ ?_
    · rw [div_eq_mul_inv]
      refine valAtLeast_mono ?_ (valAtLeast_mul (hU z hz) h2); omega
    · refine valAtLeast_mono ?_ (valAtLeast_mul (hV2 z hz) (hS z hz).1); omega

/-! ## ★★★★★★★★★★★★第 1079 —— 付値評価から `neronExp` の下限へ -/

/-- ☆`ValAtLeast p 0` なら付値環の元。 -/
theorem mem_primeSubring_of_valAtLeast {p : HeightOneSpectrum (𝓞 L)} {z : L}
    (h : ValAtLeast p 0 z) : z ∈ primeSubring p := by
  rcases eq_or_ne z 0 with rfl | hz
  · exact zero_mem _
  · exact (mem_primeSubring_iff p z).2 ((valAdd_nonneg_iff p _).1 (h hz))

/-- ☆整モデルにする変数変換があれば `neronExp ≥ valAdd u`。 -/
theorem valAdd_u_le_neronExp (p : HeightOneSpectrum (𝓞 L)) (W : WeierstrassCurve L)
    (hΔ : W.Δ ≠ 0) (C : WeierstrassCurve.VariableChange L)
    (hint : (C • W).IsIntegral (primeSubring p)) : valAdd p C.u ≤ neronExp p W := by
  have h := neronExp_nonneg p (C • W) (variableChange_Delta_ne_zero W hΔ C) hint
  rw [neronExp_variableChange p W hΔ C] at h
  omega

/-- ★★★★★★★★★★★★★★★★**係数の付値から `neronExp` の下限**（第 1079）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆`v(a_i) ≥ −i·M` で `v(c) = M` なる `c` があれば、
`u = c⁻¹` の変数変換で整になるので `neronExp ≥ −M`。 -/
theorem neronExp_ge_of_valAtLeast (p : HeightOneSpectrum (𝓞 L)) (W : WeierstrassCurve L)
    (hΔ : W.Δ ≠ 0) {c : L} (hc : c ≠ 0) (M : ℤ)
    (hcM : valAdd p (Units.mk0 c hc) = M)
    (h1 : ValAtLeast p (-M) W.a₁) (h2 : ValAtLeast p (-2 * M) W.a₂)
    (h3 : ValAtLeast p (-3 * M) W.a₃) (h4 : ValAtLeast p (-4 * M) W.a₄)
    (h6 : ValAtLeast p (-6 * M) W.a₆) :
    -M ≤ neronExp p W := by
  set C : WeierstrassCurve.VariableChange L :=
    { u := (Units.mk0 c hc)⁻¹, r := 0, s := 0, t := 0 } with hCdef
  have hCu : valAdd p C.u = -M := by
    rw [hCdef]
    show valAdd p (Units.mk0 c hc)⁻¹ = -M
    rw [valAdd_inv, hcM]
  have hstep : ∀ (i : ℕ) (a : L), ValAtLeast p (-(i : ℤ) * M) a →
      (c ^ i * a) ∈ primeSubring p := by
    intro i a ha
    refine mem_primeSubring_of_valAtLeast ?_
    have hpow : ValAtLeast p ((i : ℤ) * M) (c ^ i) := by
      have h0 := valAtLeast_unit p (Units.mk0 c hc ^ i)
      rw [valAdd_pow, hcM] at h0
      simpa using h0
    exact valAtLeast_mono (le_of_eq (by ring)) (valAtLeast_mul hpow ha)
  have hint : (C • W).IsIntegral (primeSubring p) := by
    refine isIntegral_of_mem _ _ ?_ ?_ ?_ ?_ ?_
    · rw [WeierstrassCurve.variableChange_a₁, hCdef]
      simpa using hstep 1 W.a₁ (valAtLeast_mono (by omega) h1)
    · rw [WeierstrassCurve.variableChange_a₂, hCdef]
      simpa using hstep 2 W.a₂ (valAtLeast_mono (by omega) h2)
    · rw [WeierstrassCurve.variableChange_a₃, hCdef]
      simpa using hstep 3 W.a₃ (valAtLeast_mono (by omega) h3)
    · rw [WeierstrassCurve.variableChange_a₄, hCdef]
      simpa using hstep 4 W.a₄ (valAtLeast_mono (by omega) h4)
    · rw [WeierstrassCurve.variableChange_a₆, hCdef]
      simpa using hstep 6 W.a₆ (valAtLeast_mono (by omega) h6)
  have := valAdd_u_le_neronExp p W hΔ C hint
  omega

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


/-! ## ★★★★★★★★第 1081 —— `p ∣ l`（`l` 奇）なら `v_p(2) = 0` -/

/-- ☆整数は付値環の元。 -/
theorem valAtLeast_intCast (p : HeightOneSpectrum (𝓞 L)) (a : ℤ) :
    ValAtLeast p 0 ((a : L)) := by
  refine valAtLeast_of_mem ?_
  simpa using ((a : primeSubring p)).2

/-- ★★★★★★★★**`l` が奇素数で `p ∣ l` なら `v_p(2) = 0`**（第 1081）。

☆Bézout: `2a + lb = 1`。両方の付値が正なら `1` の付値も正になり矛盾する。 -/
theorem valAtLeast_two_inv_of_dvd (p : HeightOneSpectrum (𝓞 L))
    {l : ℕ} (hl : l.Prime) (hodd : l ≠ 2)
    (hlne : ((l : L)) ≠ 0) (hpos : 0 < valAdd p (Units.mk0 ((l : L)) hlne)) :
    ValAtLeast p 0 ((2 : L)⁻¹) := by
  have h2ne : ((2 : L)) ≠ 0 := two_ne_zero
  have h2nn : 0 ≤ valAdd p (Units.mk0 ((2 : L)) h2ne) := by
    rw [valAdd_nonneg_iff]
    refine (mem_primeSubring_iff p _).1 ?_
    simpa using ((2 : primeSubring p)).2
  have h2z : valAdd p (Units.mk0 ((2 : L)) h2ne) = 0 := by
    by_contra hne
    have hge : 1 ≤ valAdd p (Units.mk0 ((2 : L)) h2ne) := by omega
    -- ☆Bézout
    obtain ⟨a, b, hab⟩ : ∃ a b : ℤ, a * 2 + b * (l : ℤ) = 1 := by
      have hnc : Nat.Coprime 2 l :=
        (Nat.coprime_primes Nat.prime_two hl).2 (fun hc => hodd hc.symm)
      have hcop : IsCoprime (2 : ℤ) (l : ℤ) := by
        rw [Int.isCoprime_iff_gcd_eq_one]
        show Int.gcd 2 (l : ℤ) = 1
        unfold Int.gcd
        simpa using hnc
      obtain ⟨a, b, h⟩ := hcop
      exact ⟨a, b, h⟩
    have hone : ((a : L)) * ((2 : L)) + ((b : L)) * ((l : L)) = 1 := by
      have := congrArg (fun z : ℤ => ((z : L))) hab
      push_cast at this
      exact this
    have hA : ValAtLeast p 1 (((a : L)) * ((2 : L))) := by
      refine valAtLeast_mono ?_ (valAtLeast_mul (valAtLeast_intCast p a)
        (valAtLeast_mono hge (valAtLeast_unit p (Units.mk0 ((2 : L)) h2ne))))
      omega
    have hB : ValAtLeast p 1 (((b : L)) * ((l : L))) := by
      refine valAtLeast_mono ?_ (valAtLeast_mul (valAtLeast_intCast p b)
        (valAtLeast_mono hpos (valAtLeast_unit p (Units.mk0 ((l : L)) hlne))))
      omega
    have hsum := valAtLeast_add hA hB
    rw [hone] at hsum
    have hfin := hsum one_ne_zero
    have h1z : valAdd p (Units.mk0 ((1 : L)) one_ne_zero) = 0 := by
      rw [← valAdd_one p]
      exact valAdd_eq_of_valuation_eq p _ _ (by simp)
    omega
  intro hz
  have hinv : valAdd p (Units.mk0 ((2 : L)) h2ne)⁻¹ = 0 := by
    rw [valAdd_inv, h2z, neg_zero]
  refine le_of_eq ?_
  rw [← hinv]
  exact (valAdd_eq_of_valuation_eq p _ _ (by simp)).symm

open Finset in
/-- ★★★★★★★★★★★★★★★★★★★★**[GenEll] `p ∣ l` でも `neronExp p E′ ≥ −v_p(l)`**（第 1080）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆第 1077（座標の深さ）→ 第 1078（`v`・`w` の付値）→ 第 1079（変数変換）。
★`h2` は `v_p(2) = 0`。`l` が奇素数で `p ∣ l` なら自動的に成り立つ。 -/
theorem neronExp_veluQuotientFull_ge (p : HeightOneSpectrum (𝓞 L))
    (E : WeierstrassCurve L) [WeierstrassCurve.IsIntegral (primeSubring p) E]
    {l : ℕ} (hl : l.Prime) (hodd : l ≠ 2)
    (Q : E.toAffine.Point) (hQ : addOrderOf Q = l)
    (hlne : ((l : L)) ≠ 0) (h2 : ValAtLeast p 0 ((2 : L)⁻¹))
    (hΔ : (veluQuotientFull E (((range l).erase 0).image
      (fun k : ℕ => pointCoords (k • Q)))).Δ ≠ 0) :
    -(valAdd p (Units.mk0 ((l : L)) hlne))
      ≤ neronExp p (veluQuotientFull E (((range l).erase 0).image
          (fun k : ℕ => pointCoords (k • Q)))) := by
  classical
  set M : ℤ := valAdd p (Units.mk0 ((l : L)) hlne) with hMdef
  have hM0 : 0 ≤ M := by
    rw [hMdef, valAdd_nonneg_iff]
    refine (mem_primeSubring_iff p _).1 ?_
    simpa using ((l : primeSubring p)).2
  have hS : ∀ z ∈ ((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)),
      ValAtLeast p (-2 * M) z.1 ∧ ValAtLeast p (-3 * M) z.2 := by
    intro z hz
    obtain ⟨k, hk, rfl⟩ := Finset.mem_image.1 hz
    rw [mem_erase, mem_range] at hk
    exact valAtLeast_pointCoords_of_le p E hl hodd Q hQ M hM0 hlne (by omega) hk.1 hk.2
  obtain ⟨hv, hw⟩ := valAtLeast_veluVFull_veluWFull p E _ M hM0 hS h2
  obtain ⟨ha1, ha2, ha3, ha4, ha6⟩ := mem_primeSubring_of_isIntegral p E
  have k1 : ValAtLeast p 0 E.a₁ := valAtLeast_of_mem ha1
  have k2 : ValAtLeast p 0 E.a₂ := valAtLeast_of_mem ha2
  have k3 : ValAtLeast p 0 E.a₃ := valAtLeast_of_mem ha3
  have k4 : ValAtLeast p 0 E.a₄ := valAtLeast_of_mem ha4
  have k6 : ValAtLeast p 0 E.a₆ := valAtLeast_of_mem ha6
  have kb2 : ValAtLeast p 0 E.b₂ := by
    rw [WeierstrassCurve.b₂]
    refine valAtLeast_add ?_ ?_
    · have hh := valAtLeast_mul k1 k1
      have hsq : E.a₁ * E.a₁ = E.a₁ ^ 2 := by ring
      rw [hsq] at hh
      exact valAtLeast_mono (by omega) hh
    · have hc4 : ValAtLeast p 0 ((4 : L)) := by
        refine valAtLeast_of_mem ?_
        simpa using (4 : primeSubring p).2
      exact valAtLeast_mono (by omega) (valAtLeast_mul hc4 k2)
  have hc5 : ValAtLeast p 0 ((5 : L)) := by
    refine valAtLeast_of_mem ?_
    simpa using (5 : primeSubring p).2
  have hc7 : ValAtLeast p 0 ((7 : L)) := by
    refine valAtLeast_of_mem ?_
    simpa using (7 : primeSubring p).2
  refine neronExp_ge_of_valAtLeast p _ hΔ hlne M hMdef.symm ?_ ?_ ?_ ?_ ?_
  · show ValAtLeast p (-M) (veluCurve E _ _).a₁
    rw [veluCurve_a₁]; exact valAtLeast_mono (by omega) k1
  · show ValAtLeast p (-2 * M) (veluCurve E _ _).a₂
    rw [veluCurve_a₂]; exact valAtLeast_mono (by omega) k2
  · show ValAtLeast p (-3 * M) (veluCurve E _ _).a₃
    rw [veluCurve_a₃]; exact valAtLeast_mono (by omega) k3
  · show ValAtLeast p (-4 * M) (veluCurve E _ _).a₄
    rw [veluCurve_a₄, sub_eq_add_neg]
    refine valAtLeast_add ?_ ?_
    · exact valAtLeast_mono (by omega) k4
    · exact valAtLeast_mono (by omega) (valAtLeast_neg (valAtLeast_mul hc5 hv))
  · show ValAtLeast p (-6 * M) (veluCurve E _ _).a₆
    rw [veluCurve_a₆, sub_eq_add_neg, sub_eq_add_neg]
    refine valAtLeast_add (valAtLeast_add ?_ ?_) ?_
    · exact valAtLeast_mono (by omega) k6
    · exact valAtLeast_mono (by omega) (valAtLeast_neg (valAtLeast_mul kb2 hv))
    · exact valAtLeast_mono (by omega) (valAtLeast_neg (valAtLeast_mul hc7 hw))


open Finset in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★**[GenEll] `neronExp p E − neronExp p E′ ≤ v_p(l)`**（第 1082）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★★**大域極小モデルは不要である**——差は変数変換で不変なので、
**各素点ごとに別々に** `p` で極小なモデルを取ってよい。
☆`p ∤ l` なら第 1074（`E′` が整）、`p ∣ l` なら第 1080。 -/
theorem neronExp_sub_le_valAdd_natCast (p : HeightOneSpectrum (𝓞 L))
    (E E' : WeierstrassCurve L) [E.IsElliptic] [E'.IsElliptic]
    {l : ℕ} (hl : l.Prime) (hodd : l ≠ 2)
    (Q : E.toAffine.Point) (hQ : addOrderOf Q = l)
    (hE' : E' = veluQuotientFull E (((range l).erase 0).image
      (fun k : ℕ => pointCoords (k • Q))))
    (hlne : ((l : L)) ≠ 0) :
    neronExp p E - neronExp p E' ≤ valAdd p (Units.mk0 ((l : L)) hlne) := by
  classical
  have hΔE : E.Δ ≠ 0 := E.isUnit_Δ.ne_zero
  have hΔE' : E'.Δ ≠ 0 := E'.isUnit_Δ.ne_zero
  obtain ⟨C, hC⟩ := WeierstrassCurve.exists_isMinimal (primeSubring p) E
  haveI hCE : (C • E).IsElliptic := by
    rw [WeierstrassCurve.isElliptic_iff, WeierstrassCurve.variableChange_Δ]
    exact ((C.u⁻¹).isUnit.pow 12).mul E.isUnit_Δ
  haveI hCE' : (C • E').IsElliptic := by
    rw [WeierstrassCurve.isElliptic_iff, WeierstrassCurve.variableChange_Δ]
    exact ((C.u⁻¹).isUnit.pow 12).mul E'.isUnit_Δ
  haveI := hC
  haveI hCint : (C • E).IsIntegral (primeSubring p) := inferInstance
  -- ☆差は変数変換で不変（第 1053）
  have hinv : neronExp p (C • E) - neronExp p (C • E') = neronExp p E - neronExp p E' := by
    rw [neronExp_variableChange p E hΔE C, neronExp_variableChange p E' hΔE' C]
    ring
  -- ☆`C • E` は極小なので `neronExp = 0`
  have hzero : neronExp p (C • E) = 0 := by
    have h := neronExp_eq p (C • E) (variableChange_Delta_ne_zero E hΔE C) 1 (by
      rw [one_smul]; exact hC)
    rw [h]
    show valAdd p (1 : Lˣ) = 0
    exact valAdd_one p
  -- ☆Vélu の商を変数変換で運ぶ（第 969）
  have hQ' : addOrderOf (ABC3.Found.GenEll.vcPoint C E Q) = l := by
    rw [ABC3.Found.GenEll.addOrderOf_vcPoint C E Q]; exact hQ
  have hEq : C • E' = veluQuotientFull (C • E) (((range l).erase 0).image
      (fun k : ℕ => pointCoords (k • ABC3.Found.GenEll.vcPoint C E Q))) :=
    veluQuotientFull_vcPoint_eq C E E' hQ (two_ne_zero) hE'
  have hΔv : (veluQuotientFull (C • E) (((range l).erase 0).image
      (fun k : ℕ => pointCoords (k • ABC3.Found.GenEll.vcPoint C E Q)))).Δ ≠ 0 := by
    rw [← hEq]; exact (C • E').isUnit_Δ.ne_zero
  have hM0 : 0 ≤ valAdd p (Units.mk0 ((l : L)) hlne) := by
    rw [valAdd_nonneg_iff]
    refine (mem_primeSubring_iff p _).1 ?_
    simpa using ((l : primeSubring p)).2
  -- ★場合分け
  have hge : -(valAdd p (Units.mk0 ((l : L)) hlne)) ≤ neronExp p (C • E') := by
    rcases eq_or_lt_of_le hM0 with hz | hpos
    · -- ☆`p ∤ l`: 第 1074 で `E'` は整
      haveI hint : (C • E').IsIntegral (primeSubring p) := by
        rw [hEq]
        exact isIntegral_veluQuotientFull_of_addOrderOf_prime p (C • E) hl hodd
          (isUnit_natCast_primeSubring p hlne hz.symm) (ABC3.Found.GenEll.vcPoint C E Q) hQ'
      have := neronExp_nonneg p (C • E') (variableChange_Delta_ne_zero E' hΔE' C) hint
      omega
    · -- ★`p ∣ l`: 第 1080
      have h2 := valAtLeast_two_inv_of_dvd p hl hodd hlne hpos
      have := neronExp_veluQuotientFull_ge p (C • E) hl hodd (ABC3.Found.GenEll.vcPoint C E Q) hQ' hlne h2 hΔv
      rw [hEq]
      exact this
  omega

end Integral

/-! ## ★出典の紐付け(`.src`) -/

def preΨ_eval_eq_zero_of_addOrderOf_prime.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(素数位数の点の x 座標は preΨ_l の根。★無条件)",
    sectionId := "genell-lemma-3-5" }

end ABC3.Found.GaloisRep
