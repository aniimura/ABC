import ABC3.Found.PGC.LubinTateSequence

/-!
# Lubin-Tate 形式群の存在補題: 近似列の極限(進行中)

`Found/PGC/LubinTateSequence.lean::ΦSeq` の係数の**安定性**(`k` が十分
大きいとき、与えられた次数の係数がそれ以上変化しない)を確立し、極限
`F` を係数関数として直接定義し、`F` が関数等式を exact に満たすことを
示す。
-/

namespace ABC3.Found.PGC

variable {A : Type*} [CommRing A] [IsLocalRing A] [IsDomain A]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField A) pp] [Fintype (IsLocalRing.ResidueField A)]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField A) = pp ^ ff)
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π})
    (g : PowerSeries A) (hg0 : PowerSeries.constantCoeff g = 0)
    (hg1 : PowerSeries.coeff 1 g = π)
    (hgres : PowerSeries.map (IsLocalRing.residue A) g = PowerSeries.X ^ (pp ^ ff))
    (f : PowerSeries A) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hfres : PowerSeries.map (IsLocalRing.residue A) f = PowerSeries.X ^ (pp ^ ff))

include hq hπmax hg0 hg1 hgres hf0 hf1 hfres in
/-- `ΦSeq` の逐次関係を露出させる: `ΦSeq (k+1)` は `ΦSeq k` に、次数 `k+2`
の斉次式 `φ` を足したもの。 -/
theorem ΦSeq_succ_eq (k : ℕ) :
    ∃ φ : MvPowerSeries (Fin 2) A,
      (∀ d : Fin 2 →₀ ℕ, Finsupp.degree d ≠ k + 2 → MvPowerSeries.coeff d φ = 0) ∧
      (ΦSeq hq hπmax g hg0 hg1 hgres f hf0 hf1 hfres (k + 1)).1 =
        (ΦSeq hq hπmax g hg0 hg1 hgres f hf0 hf1 hfres k).1 + φ := by
  set prev := ΦSeq hq hπmax g hg0 hg1 hgres f hf0 hf1 hfres k with hprev_def
  set φex := obstructionVanishesUpTo_step hq hπmax g hg0 hg1 hgres f hf0 hf1 hfres
    prev.2.2 (Nat.succ_ne_zero k) prev.2.1 with hφex_def
  refine ⟨φex.choose, φex.choose_spec.1, ?_⟩
  show (ΦSeq hq hπmax g hg0 hg1 hgres f hf0 hf1 hfres (k + 1)).1 = prev.1 + φex.choose
  rfl

include hq hπmax hg0 hg1 hgres hf0 hf1 hfres in
/-- `ΦSeq` の差の次数は、離れているほど高い: `k≤m` のとき
`(ΦSeq m).1 − (ΦSeq k).1` の次数は `≥k+2`。 -/
theorem order_diff_ΦSeq_ge (k : ℕ) : ∀ m, k ≤ m →
    ((k + 2 : ℕ) : ℕ∞) ≤
      ((ΦSeq hq hπmax g hg0 hg1 hgres f hf0 hf1 hfres m).1 -
        (ΦSeq hq hπmax g hg0 hg1 hgres f hf0 hf1 hfres k).1).order := by
  intro m hkm
  induction m, hkm using Nat.le_induction with
  | base => simp [MvPowerSeries.order_zero]
  | succ n hn ih =>
    obtain ⟨φ, hφhomog, hφeq⟩ := ΦSeq_succ_eq hq hπmax g hg0 hg1 hgres f hf0 hf1 hfres n
    have hstep : (ΦSeq hq hπmax g hg0 hg1 hgres f hf0 hf1 hfres (n + 1)).1 -
        (ΦSeq hq hπmax g hg0 hg1 hgres f hf0 hf1 hfres k).1 =
        ((ΦSeq hq hπmax g hg0 hg1 hgres f hf0 hf1 hfres n).1 -
          (ΦSeq hq hπmax g hg0 hg1 hgres f hf0 hf1 hfres k).1) + φ := by
      rw [hφeq]; ring
    rw [hstep]
    have hφorder : ((n + 2 : ℕ) : ℕ∞) ≤ φ.order := by
      apply MvPowerSeries.nat_le_order
      intro d hd
      exact hφhomog d (by omega)
    have hnk : ((k + 2 : ℕ) : ℕ∞) ≤ ((n + 2 : ℕ) : ℕ∞) := by exact_mod_cast (by omega : k + 2 ≤ n + 2)
    have hmin : ((k + 2 : ℕ) : ℕ∞) ≤
        min ((ΦSeq hq hπmax g hg0 hg1 hgres f hf0 hf1 hfres n).1 -
              (ΦSeq hq hπmax g hg0 hg1 hgres f hf0 hf1 hfres k).1).order φ.order := by
      rw [le_min_iff]; exact ⟨ih, le_trans hnk hφorder⟩
    exact le_trans hmin MvPowerSeries.min_order_le_add

include hq hπmax hg0 hg1 hgres hf0 hf1 hfres in
/-- **係数の安定性**: `degree e ≤ k+1` のとき、`ΦSeq m`(`m≥k`)の `e` 係数は
`ΦSeq k` のそれと一致する——`k` を超えて先の項を足しても、次数 `e` の係数は
もう変化しない。 -/
theorem coeff_ΦSeq_stable (e : Fin 2 →₀ ℕ) (k m : ℕ) (hkm : k ≤ m)
    (he : Finsupp.degree e ≤ k + 1) :
    MvPowerSeries.coeff e (ΦSeq hq hπmax g hg0 hg1 hgres f hf0 hf1 hfres m).1 =
      MvPowerSeries.coeff e (ΦSeq hq hπmax g hg0 hg1 hgres f hf0 hf1 hfres k).1 := by
  have hord := order_diff_ΦSeq_ge hq hπmax g hg0 hg1 hgres f hf0 hf1 hfres k m hkm
  have hlt : ((Finsupp.degree e : ℕ) : ℕ∞) <
      ((ΦSeq hq hπmax g hg0 hg1 hgres f hf0 hf1 hfres m).1 -
        (ΦSeq hq hπmax g hg0 hg1 hgres f hf0 hf1 hfres k).1).order := by
    calc ((Finsupp.degree e : ℕ) : ℕ∞) ≤ ((k + 1 : ℕ) : ℕ∞) := by exact_mod_cast he
      _ < ((k + 2 : ℕ) : ℕ∞) := by exact_mod_cast (by omega : k + 1 < k + 2)
      _ ≤ _ := hord
  have hz := MvPowerSeries.coeff_of_lt_order hlt
  rw [map_sub] at hz
  exact sub_eq_zero.mp hz

include hq hπmax hg0 hg1 hgres hf0 hf1 hfres in
/-- ★★★**Lubin-Tate 形式群**。`ΦSeq` の極限——`e` の係数は、`e` の次数まで
進んだ近似 `ΦSeq (degree e)` から直接読み取る(`coeff_ΦSeq_stable` により、
それより先の近似でも同じ値になる)。 -/
noncomputable def LubinTateF : MvPowerSeries (Fin 2) A :=
  fun e => MvPowerSeries.coeff e
    (ΦSeq hq hπmax g hg0 hg1 hgres f hf0 hf1 hfres (Finsupp.degree e)).1

include hq hπmax hg0 hg1 hgres hf0 hf1 hfres in
theorem coeff_LubinTateF (e : Fin 2 →₀ ℕ) :
    MvPowerSeries.coeff e (LubinTateF hq hπmax g hg0 hg1 hgres f hf0 hf1 hfres) =
      MvPowerSeries.coeff e (ΦSeq hq hπmax g hg0 hg1 hgres f hf0 hf1 hfres (Finsupp.degree e)).1 :=
  rfl

include hq hπmax hg0 hg1 hgres hf0 hf1 hfres in
theorem coeff_LubinTateF_eq_ΦSeq (e : Fin 2 →₀ ℕ) (m : ℕ) (he : Finsupp.degree e ≤ m) :
    MvPowerSeries.coeff e (LubinTateF hq hπmax g hg0 hg1 hgres f hf0 hf1 hfres) =
      MvPowerSeries.coeff e (ΦSeq hq hπmax g hg0 hg1 hgres f hf0 hf1 hfres m).1 := by
  rw [coeff_LubinTateF]
  exact (coeff_ΦSeq_stable hq hπmax g hg0 hg1 hgres f hf0 hf1 hfres e (Finsupp.degree e) m he
    (by omega)).symm

include hq hπmax hg0 hg1 hgres hf0 hf1 hfres in
/-- `F − ΦSeq m` の次数は `≥m+1`。 -/
theorem order_diff_LubinTateF_ΦSeq_ge (m : ℕ) :
    ((m + 1 : ℕ) : ℕ∞) ≤
      (LubinTateF hq hπmax g hg0 hg1 hgres f hf0 hf1 hfres -
        (ΦSeq hq hπmax g hg0 hg1 hgres f hf0 hf1 hfres m).1).order := by
  apply MvPowerSeries.nat_le_order
  intro d hd
  have hdm : Finsupp.degree d ≤ m := by omega
  have heq := coeff_LubinTateF_eq_ΦSeq hq hπmax g hg0 hg1 hgres f hf0 hf1 hfres d m hdm
  rw [map_sub, heq, sub_self]

include hg0 in
omit [IsLocalRing A] [IsDomain A] [Fintype (IsLocalRing.ResidueField A)] in
/-- `a := fun i => g.subst(X_i)` の各成分は次数 `≥1`。 -/
theorem hai_order_g (i : Fin 2) :
    (1 : ℕ∞) ≤ (PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries (Fin 2) A) g).order := by
  apply MvPowerSeries.one_le_order_iff_constCoeff_eq_zero.mpr
  exact PowerSeries.constantCoeff_subst_eq_zero (MvPowerSeries.constantCoeff_X i) g hg0

include hq hπmax hg0 hg1 hgres hf0 hf1 hfres in
/-- **g側の安定性**: `F` を `g` で代入した結果は、`e` の次数以上まで進んだ
`ΦSeq m` を代入した結果と、次数 `e` で一致する。 -/
theorem coeff_gsubst_LubinTateF_eq_ΦSeq (e : Fin 2 →₀ ℕ) (m : ℕ) (he : Finsupp.degree e ≤ m) :
    MvPowerSeries.coeff e
        (MvPowerSeries.subst (fun i => PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries (Fin 2) A) g)
          (LubinTateF hq hπmax g hg0 hg1 hgres f hf0 hf1 hfres)) =
      MvPowerSeries.coeff e
        (MvPowerSeries.subst (fun i => PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries (Fin 2) A) g)
          (ΦSeq hq hπmax g hg0 hg1 hgres f hf0 hf1 hfres m).1) := by
  have haHS := hasSubst_g_subst_X g hg0
  have hsub : MvPowerSeries.subst (fun i => PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries (Fin 2) A) g)
        (LubinTateF hq hπmax g hg0 hg1 hgres f hf0 hf1 hfres) -
      MvPowerSeries.subst (fun i => PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries (Fin 2) A) g)
        (ΦSeq hq hπmax g hg0 hg1 hgres f hf0 hf1 hfres m).1 =
      MvPowerSeries.subst (fun i => PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries (Fin 2) A) g)
        (LubinTateF hq hπmax g hg0 hg1 hgres f hf0 hf1 hfres -
          (ΦSeq hq hπmax g hg0 hg1 hgres f hf0 hf1 hfres m).1) :=
    (MvPowerSeries.subst_sub haHS _ _).symm
  have horder : ((m + 1 : ℕ) : ℕ∞) ≤
      (MvPowerSeries.subst (fun i => PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries (Fin 2) A) g)
        (LubinTateF hq hπmax g hg0 hg1 hgres f hf0 hf1 hfres -
          (ΦSeq hq hπmax g hg0 hg1 hgres f hf0 hf1 hfres m).1)).order := by
    have hle := MvPowerSeries.le_order_subst haHS
      (LubinTateF hq hπmax g hg0 hg1 hgres f hf0 hf1 hfres -
        (ΦSeq hq hπmax g hg0 hg1 hgres f hf0 hf1 hfres m).1)
    have hinf : (1 : ℕ∞) ≤ ⨅ i, (PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries (Fin 2) A) g).order :=
      le_iInf (hai_order_g g hg0)
    have hdiffge := order_diff_LubinTateF_ΦSeq_ge hq hπmax g hg0 hg1 hgres f hf0 hf1 hfres m
    calc ((m + 1 : ℕ) : ℕ∞) = 1 * ((m + 1 : ℕ) : ℕ∞) := by ring
      _ ≤ (⨅ i, (PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries (Fin 2) A) g).order) *
            ((m + 1 : ℕ) : ℕ∞) := by gcongr
      _ ≤ (⨅ i, (PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries (Fin 2) A) g).order) *
            (LubinTateF hq hπmax g hg0 hg1 hgres f hf0 hf1 hfres -
              (ΦSeq hq hπmax g hg0 hg1 hgres f hf0 hf1 hfres m).1).order := by gcongr
      _ ≤ _ := hle
  have hz : MvPowerSeries.coeff e
      (MvPowerSeries.subst (fun i => PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries (Fin 2) A) g)
        (LubinTateF hq hπmax g hg0 hg1 hgres f hf0 hf1 hfres -
          (ΦSeq hq hπmax g hg0 hg1 hgres f hf0 hf1 hfres m).1)) = 0 := by
    apply MvPowerSeries.coeff_of_lt_order
    calc ((Finsupp.degree e : ℕ) : ℕ∞) ≤ ((m : ℕ) : ℕ∞) := by exact_mod_cast he
      _ < ((m + 1 : ℕ) : ℕ∞) := by exact_mod_cast (by omega : m < m + 1)
      _ ≤ _ := horder
  rw [← hsub, map_sub] at hz
  exact sub_eq_zero.mp hz

include hq hπmax hg0 hg1 hgres hf0 hf1 hfres in
/-- **f側の安定性**: `f` に `F` を代入した結果は、`e` の次数以上まで進んだ
`ΦSeq m` を代入した結果と、次数 `e` で一致する。`coeff_subst_linearize`
(f側の線形化)を、`degree e < φ.order` のとき `π·coeff e φ=0` になることで
再利用する。 -/
theorem coeff_fsubst_LubinTateF_eq_ΦSeq (e : Fin 2 →₀ ℕ) (m : ℕ) (he : Finsupp.degree e ≤ m) :
    MvPowerSeries.coeff e (PowerSeries.subst (LubinTateF hq hπmax g hg0 hg1 hgres f hf0 hf1 hfres) f) =
      MvPowerSeries.coeff e
        (PowerSeries.subst (ΦSeq hq hπmax g hg0 hg1 hgres f hf0 hf1 hfres m).1 f) := by
  set Φm := (ΦSeq hq hπmax g hg0 hg1 hgres f hf0 hf1 hfres m).1 with hΦm_def
  set φ := LubinTateF hq hπmax g hg0 hg1 hgres f hf0 hf1 hfres - Φm with hφ_def
  have hΦmadd : Φm + φ = LubinTateF hq hπmax g hg0 hg1 hgres f hf0 hf1 hfres := by rw [hφ_def]; ring
  have hΦm0 : MvPowerSeries.constantCoeff Φm = 0 :=
    (ΦSeq hq hπmax g hg0 hg1 hgres f hf0 hf1 hfres m).2.2
  have hΦmorder : (1 : ℕ∞) ≤ Φm.order :=
    MvPowerSeries.one_le_order_iff_constCoeff_eq_zero.mpr hΦm0
  have hφorder : ((m + 1 : ℕ) : ℕ∞) ≤ φ.order :=
    order_diff_LubinTateF_ΦSeq_ge hq hπmax g hg0 hg1 hgres f hf0 hf1 hfres m
  have hφ0 : MvPowerSeries.constantCoeff φ = 0 := by
    rw [← MvPowerSeries.coeff_zero_eq_constantCoeff_apply]
    apply MvPowerSeries.coeff_of_lt_order
    calc ((Finsupp.degree (0 : Fin 2 →₀ ℕ) : ℕ) : ℕ∞) = (0 : ℕ∞) := by simp
      _ < ((m + 1 : ℕ) : ℕ∞) := by exact_mod_cast (by omega : 0 < m + 1)
      _ ≤ _ := hφorder
  have hlin := coeff_subst_linearize hΦm0 hφ0 hΦmorder hφorder (by omega : 1 ≤ m + 1) f π hf0 hf1 e
    (by exact_mod_cast (by omega : Finsupp.degree e ≤ m + 1))
  rw [hΦmadd] at hlin
  have hφe0 : MvPowerSeries.coeff e φ = 0 := by
    apply MvPowerSeries.coeff_of_lt_order
    calc ((Finsupp.degree e : ℕ) : ℕ∞) ≤ ((m : ℕ) : ℕ∞) := by exact_mod_cast he
      _ < ((m + 1 : ℕ) : ℕ∞) := by exact_mod_cast (by omega : m < m + 1)
      _ ≤ _ := hφorder
  rw [hφe0, mul_zero] at hlin
  exact sub_eq_zero.mp hlin

include hq hπmax hg0 hg1 hgres hf0 hf1 hfres in
/-- ★★★★★★★★**Lubin-Tate 形式群の存在補題**。`F := LubinTateF` は関数等式
`F.subst(g,g) = f.subst(F)` を(係数ごとに)exact に満たす。`e` の次数 `m`
まで進んだ `ΦSeq m` に対する g側・f側の安定性(`F` と `ΦSeq m` は次数 `e`
で同じ値を与える)と、`ΦSeq m` 自身が満たす不変量(次数 `m+1` まで障害が
消えている)を組み合わせるだけで出る。 -/
theorem coeff_LubinTateF_functional_equation (e : Fin 2 →₀ ℕ) :
    MvPowerSeries.coeff e
        (MvPowerSeries.subst (fun i => PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries (Fin 2) A) g)
          (LubinTateF hq hπmax g hg0 hg1 hgres f hf0 hf1 hfres)) =
      MvPowerSeries.coeff e
        (PowerSeries.subst (LubinTateF hq hπmax g hg0 hg1 hgres f hf0 hf1 hfres) f) := by
  set m := Finsupp.degree e with hm_def
  have h1 := coeff_gsubst_LubinTateF_eq_ΦSeq hq hπmax g hg0 hg1 hgres f hf0 hf1 hfres e m le_rfl
  have h2 := coeff_fsubst_LubinTateF_eq_ΦSeq hq hπmax g hg0 hg1 hgres f hf0 hf1 hfres e m le_rfl
  have h3 := (ΦSeq hq hπmax g hg0 hg1 hgres f hf0 hf1 hfres m).2.1 e (by omega)
  rw [map_sub] at h3
  rw [h1, h2]
  exact sub_eq_zero.mp h3

include hq hπmax hg0 hg1 hgres hf0 hf1 hfres in
/-- ★★★★★★★★★**Lubin-Tate 形式群の存在補題(まとめ)**。`F` は関数等式
`F.subst(g,g) = f.subst(F)` を power series の等式として満たす——
`coeff_LubinTateF_functional_equation` を `MvPowerSeries.ext` で束ねるだけ。 -/
theorem LubinTateF_functional_equation :
    MvPowerSeries.subst (fun i => PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries (Fin 2) A) g)
        (LubinTateF hq hπmax g hg0 hg1 hgres f hf0 hf1 hfres) =
      PowerSeries.subst (LubinTateF hq hπmax g hg0 hg1 hgres f hf0 hf1 hfres) f :=
  MvPowerSeries.ext (coeff_LubinTateF_functional_equation hq hπmax g hg0 hg1 hgres f hf0 hf1 hfres)

end ABC3.Found.PGC
