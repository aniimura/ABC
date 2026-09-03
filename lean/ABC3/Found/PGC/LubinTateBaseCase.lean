import ABC3.Found.PGC.LubinTateStepAssembly

/-!
# Lubin-Tate 形式群の存在補題: 次数ごとの再帰構成の出発点(`sorry` 無し)

`Found/PGC/LubinTateStepAssembly.lean::exists_next_step` は `n ≠ 0` を要求
するため(`exists_step_solution` が `1−π^n` の単数性に `n≠0` を使う)、
`Φ_0 := 0` から `n=0` で踏み出すことはできない。代わりに、加法的形式群
`Φ_1 := X_0+X_1`(いちばん自然な「次数1までの解」)が実際に障害を次数 `≤1`
まで消していることを直接示し、次数ごとの再帰の**出発点**とする。
-/

namespace ABC3.Found.PGC

/-- `PowerSeries.subst (0 : MvPowerSeries (Fin 2) A) f` は(`f` の定数項が0のとき)
恒等的に0——`coeff_subst` の有限和が `d=0` の項だけに潰れ、それが `hf0` で消える。 -/
theorem coeff_subst_zero_eq_zero {A : Type*} [CommRing A] (f : PowerSeries A)
    (hf0 : PowerSeries.coeff 0 f = 0) (e : Fin 2 →₀ ℕ) :
    MvPowerSeries.coeff e (PowerSeries.subst (0 : MvPowerSeries (Fin 2) A) f) = 0 := by
  have h0HS : PowerSeries.HasSubst (0 : MvPowerSeries (Fin 2) A) := by
    show IsNilpotent (MvPowerSeries.constantCoeff (0 : MvPowerSeries (Fin 2) A))
    rw [map_zero]; exact IsNilpotent.zero
  rw [PowerSeries.coeff_subst h0HS]
  rw [finsum_eq_sum_of_support_subset _ (s := ({0} : Finset ℕ)) (fun d hd => by
    simp only [Function.mem_support] at hd
    simp only [Finset.coe_singleton, Set.mem_singleton_iff]
    by_contra hcon
    exact hd (by rw [zero_pow hcon]; simp))]
  by_cases he : e = 0
  · simp [he, hf0]
  · simp [MvPowerSeries.coeff_one, he]

/-- 加法的形式群 `X_0+X_1` は次数1の斉次式。 -/
theorem isHomogeneous_add_X0_X1 {A : Type*} [CommRing A] :
    ∀ d : Fin 2 →₀ ℕ, Finsupp.degree d ≠ 1 →
      MvPowerSeries.coeff d ((MvPowerSeries.X 0 : MvPowerSeries (Fin 2) A) + MvPowerSeries.X 1) = 0 := by
  intro d hd
  rw [map_add, MvPowerSeries.coeff_X, MvPowerSeries.coeff_X]
  have hd0 : d ≠ Finsupp.single (0 : Fin 2) 1 := by
    intro heq; apply hd; rw [heq, finsupp_degree_fin2]; simp
  have hd1 : d ≠ Finsupp.single (1 : Fin 2) 1 := by
    intro heq; apply hd; rw [heq, finsupp_degree_fin2]; simp
  simp [hd0, hd1]

/-- ★★**加法的形式群 `X_0+X_1` は次数ごとの再帰の出発点として使える**:
`g≡πX(mod deg2)`・`f≡πX(mod deg2)` のとき、`Obstruction(X_0+X_1) :=
(X_0+X_1).subst(g,g) − f.subst(X_0+X_1)` は次数 `≤1` の範囲で消える。
g側は `coeff_subst_g_linearize`(`n:=0`、`X_0+X_1` 自身が次数1の斉次式)、
f側は `coeff_subst_linearize`(`Φ:=0`)+ `coeff_subst_zero_eq_zero` の
組み合わせで得る。 -/
theorem base_case_invariant {A : Type*} [CommRing A] {π : A}
    (g : PowerSeries A) (hg0 : PowerSeries.constantCoeff g = 0)
    (hg1 : PowerSeries.coeff 1 g = π)
    (f : PowerSeries A) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (e : Fin 2 →₀ ℕ) (he : Finsupp.degree e ≤ 1) :
    MvPowerSeries.coeff e
      (MvPowerSeries.subst (fun i => PowerSeries.subst (MvPowerSeries.X i) g)
          ((MvPowerSeries.X 0 : MvPowerSeries (Fin 2) A) + MvPowerSeries.X 1) -
        PowerSeries.subst ((MvPowerSeries.X 0 : MvPowerSeries (Fin 2) A) + MvPowerSeries.X 1) f) = 0 := by
  set Φ₁ : MvPowerSeries (Fin 2) A := (MvPowerSeries.X 0 : MvPowerSeries (Fin 2) A) + MvPowerSeries.X 1
    with hΦ₁_def
  have hGside : MvPowerSeries.coeff e
      (MvPowerSeries.subst (fun i => PowerSeries.subst (MvPowerSeries.X i) g) Φ₁) =
      π ^ (0 + 1) * MvPowerSeries.coeff e Φ₁ := by
    apply coeff_subst_g_linearize isHomogeneous_add_X0_X1 π hg0 hg1 e
    simpa using he
  have hΦ0order : (1 : ℕ∞) ≤ (0 : MvPowerSeries (Fin 2) A).order := by
    simp [MvPowerSeries.order_zero]
  have hΦ10 : MvPowerSeries.constantCoeff (0 : MvPowerSeries (Fin 2) A) = 0 := map_zero _
  have hΦ1order : (1 : ℕ∞) ≤ Φ₁.order := by
    apply MvPowerSeries.one_le_order_iff_constCoeff_eq_zero.mpr
    rw [hΦ₁_def, map_add, MvPowerSeries.constantCoeff_X, MvPowerSeries.constantCoeff_X]
    ring
  have hFside : MvPowerSeries.coeff e (PowerSeries.subst ((0 : MvPowerSeries (Fin 2) A) + Φ₁) f) -
      MvPowerSeries.coeff e (PowerSeries.subst (0 : MvPowerSeries (Fin 2) A) f) =
      π * MvPowerSeries.coeff e Φ₁ :=
    coeff_subst_linearize hΦ10 (by rw [hΦ₁_def, map_add, MvPowerSeries.constantCoeff_X,
      MvPowerSeries.constantCoeff_X]; ring) hΦ0order hΦ1order le_rfl f π hf0 hf1 e
      (by exact_mod_cast he)
  rw [zero_add] at hFside
  have hzero : MvPowerSeries.coeff e (PowerSeries.subst (0 : MvPowerSeries (Fin 2) A) f) = 0 :=
    coeff_subst_zero_eq_zero f hf0 e
  rw [hzero, sub_zero] at hFside
  rw [map_sub, hGside, hFside]
  ring

end ABC3.Found.PGC
