import ABC3.Found.PGC.LubinTateEndoStepAssembly

/-!
# `𝒪_K` 作用への拡張(1変数版)の次数ごとの再帰構成の出発点(`sorry` 無し)

`Found/PGC/LubinTateEndoStepAssembly.lean::exists_next_step_endo` は `n ≠ 0`
を要求するため、`φ_0 := 0` から `n=0` で踏み出すことはできない
(2変数版 `LubinTateBaseCase.lean` と同じ事情)。代わりに、線形関数
`φ_1 := a•X`(`a` は構成したい `[a]_f` の目標とする線形係数)が実際に
障害を次数 `≤1` まで消していることを直接示し、次数ごとの再帰の
**出発点**とする。
-/

namespace ABC3.Found.PGC

/-- `PowerSeries.subst (0 : PowerSeries A) f` は(`f` の定数項が0のとき)
恒等的に0——`PowerSeries.coeff_subst` の有限和が `d=0` の項だけに潰れ、
それが `hf0` で消える。2変数版 `coeff_subst_zero_eq_zero` の1変数版。 -/
theorem coeff_subst_zero_eq_zero_1var {A : Type*} [CommRing A] (f : PowerSeries A)
    (hf0 : PowerSeries.coeff 0 f = 0) (k : ℕ) :
    PowerSeries.coeff k (PowerSeries.subst (0 : PowerSeries A) f) = 0 := by
  have h0HS : PowerSeries.HasSubst (0 : PowerSeries A) := by
    show IsNilpotent (PowerSeries.constantCoeff (0 : PowerSeries A))
    rw [map_zero]; exact IsNilpotent.zero
  rw [show PowerSeries.coeff k (PowerSeries.subst (0 : PowerSeries A) f) =
    MvPowerSeries.coeff (Finsupp.single () k) (PowerSeries.subst (0 : PowerSeries A) f) from rfl]
  rw [PowerSeries.coeff_subst h0HS]
  rw [finsum_eq_sum_of_support_subset _ (s := ({0} : Finset ℕ)) (fun d hd => by
    simp only [Function.mem_support] at hd
    simp only [Finset.coe_singleton, Set.mem_singleton_iff]
    by_contra hcon
    exact hd (by rw [zero_pow hcon]; simp))]
  rw [Finset.sum_singleton, pow_zero, hf0, zero_smul]

/-- ★★**出発点 `a•X` は次数ごとの再帰の出発点として使える**: `g≡πX(mod
deg2)`・`f≡πX(mod deg2)` のとき、`Obstruction(a•X) := f.subst(a•X) −
g.subst(a•X)` は次数 `≤1` の範囲で消える。f側は `coeff_subst_linearize_1var`
(`Φ:=0`)+ `coeff_subst_zero_eq_zero_1var`、g側は `PowerSeries.subst g
(a•X) = a•g` という完全な等式(近似不要、`a•X` が「厳密に線形」だから)
から得る。 -/
theorem base_case_invariant_endo {A : Type*} [CommRing A] {π a : A}
    (g : PowerSeries A) (hg0 : PowerSeries.constantCoeff g = 0)
    (hg1 : PowerSeries.coeff 1 g = π)
    (f : PowerSeries A) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (k : ℕ) (hk : k ≤ 1) :
    PowerSeries.coeff k
      (PowerSeries.subst (a • (PowerSeries.X : PowerSeries A)) f -
        PowerSeries.subst g (a • (PowerSeries.X : PowerSeries A))) = 0 := by
  set φ₁ : PowerSeries A := a • (PowerSeries.X : PowerSeries A) with hφ₁_def
  have hφ₁cc : PowerSeries.constantCoeff φ₁ = 0 := by
    rw [hφ₁_def]
    show MvPowerSeries.constantCoeff (a • (PowerSeries.X : PowerSeries A)) = 0
    rw [MvPowerSeries.constantCoeff_smul]
    show a • MvPowerSeries.constantCoeff (MvPowerSeries.X () : MvPowerSeries Unit A) = 0
    rw [MvPowerSeries.constantCoeff_X, smul_zero]
  have hφ₁order : (1 : ℕ∞) ≤ MvPowerSeries.order (φ₁ : PowerSeries A) :=
    MvPowerSeries.one_le_order_iff_constCoeff_eq_zero.mpr hφ₁cc
  have h0order : (1 : ℕ∞) ≤ MvPowerSeries.order (0 : PowerSeries A) := by simp [MvPowerSeries.order_zero]
  have h00 : PowerSeries.constantCoeff (0 : PowerSeries A) = 0 := map_zero _
  have hFside := coeff_subst_linearize_1var h00 hφ₁cc h0order hφ₁order le_rfl f π hf0 hf1 k
    (by exact_mod_cast hk)
  rw [zero_add] at hFside
  have hf0zero : PowerSeries.coeff k (PowerSeries.subst (0 : PowerSeries A) f) = 0 :=
    coeff_subst_zero_eq_zero_1var f hf0 k
  rw [hf0zero, sub_zero] at hFside
  have hHSg : PowerSeries.HasSubst g := by
    show IsNilpotent (PowerSeries.constantCoeff g); rw [hg0]; exact IsNilpotent.zero
  have hGeq : PowerSeries.subst g φ₁ = a • g := by
    rw [hφ₁_def, PowerSeries.subst_smul hHSg, PowerSeries.subst_X hHSg]
  rw [map_sub, hFside, hGeq, PowerSeries.coeff_smul, PowerSeries.coeff_X, PowerSeries.coeff_smul]
  split_ifs with h
  · subst h; rw [hg1, smul_eq_mul, smul_eq_mul, mul_one]; ring
  · have hk0 : k = 0 := by omega
    subst hk0
    rw [smul_zero, mul_zero, PowerSeries.coeff_zero_eq_constantCoeff_apply, hg0, smul_zero, sub_zero]

end ABC3.Found.PGC
