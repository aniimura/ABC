import ABC3.Found.PGC.LubinTateActionEndomorphism

/-!
# `𝒪_K` 作用への拡張: `a↦[a]_f` の乗法との両立性(`sorry` 無し)

節目(2)(`Found/PGC/LubinTateActionEndomorphism.lean::formalGroupLaw_
endomorphism_of_action`)に続き、`a↦[a]_f` が環準同型であることの
**乗法側**——`[ab]_f = [a]_f ∘ [b]_f`(合成として)——を確立する。

## 証明の筋

`α:=[ab]_f`・`β:=[a]_f∘[b]_f` を構成し、両者が同じ次数1の係数(いずれも
`ab`)・同じ関数等式(`f` と可換)を満たすことを示し、**1変数の**一意性
補題(`Found/PGC/LubinTateUniqueness.lean::powerSeries_uniqueness`、
2変数版と違い `[a]_f` 自身は1変数の冪級数なのでこちらを使う)で
`α=β` を結論する。

`β` の関数等式は `[a]_f`・`[b]_f` それぞれが `f` と可換であること
(`LubinTateAction_functional_equation`)を2回使うだけで出る——
`([a]_f∘[b]_f)∘f = [a]_f∘([b]_f∘f) = [a]_f∘(f∘[b]_f) = ([a]_f∘f)∘[b]_f
= (f∘[a]_f)∘[b]_f = f∘([a]_f∘[b]_f)`、という素朴な計算だが、Lean では
1変数どうしの代入の結合律(`subst_comp_subst_single_gen`、節目(2)で
確立済み)を4回組み合わせる形になる。次数1の係数は新規の補助補題
`coeff_one_subst_1var`(代入の次数1係数は係数の積になる、という
1変数版の事実)で計算する。

加法側(`[a+b]_f = F_f([a]_f,[b]_f)`)は別途の課題として残る——2変数の
一意性補題を要するが、関数等式の導出に「代入の族との合成」という
もう1段一般的な補題が要ると見立てている。
-/

namespace ABC3.Found.PGC

/-- **代入の次数1係数は係数の積になる(1変数版)**: `v` の定数項が0のとき
`coeff 1 (v.subst(p)) = coeff 1 p * coeff 1 v`——次数≥2の項は代入後の
次数を押し上げるので効かない、という標準的な次数勘定。 -/
theorem coeff_one_subst_1var {A : Type*} [CommRing A] {v p : PowerSeries A}
    (hv0 : PowerSeries.constantCoeff v = 0) :
    PowerSeries.coeff 1 (PowerSeries.subst v p) = PowerSeries.coeff 1 p * PowerSeries.coeff 1 v := by
  have hvHS : PowerSeries.HasSubst v := by
    show IsNilpotent (PowerSeries.constantCoeff v); rw [hv0]; exact IsNilpotent.zero
  show MvPowerSeries.coeff (Finsupp.single () 1) (PowerSeries.subst v p) = _
  rw [PowerSeries.coeff_subst hvHS]
  rw [finsum_eq_sum_of_support_subset _ (s := ({1} : Finset ℕ)) (fun d hd => by
    simp only [Function.mem_support] at hd
    simp only [Finset.coe_singleton, Set.mem_singleton_iff]
    by_contra hcon
    apply hd
    have hvorder : (1 : ℕ∞) ≤ MvPowerSeries.order (v : PowerSeries A) :=
      MvPowerSeries.one_le_order_iff_constCoeff_eq_zero.mpr hv0
    rcases Nat.eq_zero_or_pos d with hd0 | hdpos
    · subst hd0
      have h01 : MvPowerSeries.coeff (Finsupp.single () (1 : ℕ)) (v ^ 0 : PowerSeries A) = 0 := by
        show PowerSeries.coeff 1 (1 : PowerSeries A) = 0
        simp
      rw [pow_zero] at h01 ⊢
      rw [h01, smul_zero]
    · have hd2 : 2 ≤ d := by omega
      have hvdorder : ((d : ℕ) : ℕ∞) ≤ MvPowerSeries.order (v ^ d : PowerSeries A) := by
        calc ((d : ℕ) : ℕ∞) = d • (1 : ℕ∞) := by simp
          _ ≤ d • MvPowerSeries.order (v : PowerSeries A) := by gcongr
          _ ≤ _ := MvPowerSeries.le_order_pow d
      have hvdorder' : ((d : ℕ) : ℕ∞) ≤ PowerSeries.order (v ^ d : PowerSeries A) := by
        rw [PowerSeries.order_eq_order]; exact hvdorder
      have h1lt : ((1 : ℕ) : ℕ∞) < ((d : ℕ) : ℕ∞) := by exact_mod_cast hd2
      have hz : MvPowerSeries.coeff (Finsupp.single () (1 : ℕ)) (v ^ d : PowerSeries A) = 0 :=
        PowerSeries.coeff_of_lt_order 1 (lt_of_lt_of_le h1lt hvdorder')
      rw [hz, smul_zero])]
  rw [Finset.sum_singleton, pow_one, smul_eq_mul]
  show PowerSeries.coeff 1 p * PowerSeries.coeff 1 v = PowerSeries.coeff 1 p * PowerSeries.coeff 1 v
  rfl

/-- ★★★★★★★**`a↦[a]_f` の乗法との両立性**: `[ab]_f = [a]_f ∘ [b]_f`
(冪級数の合成として)。1変数の一意性補題(`powerSeries_uniqueness`)を、
`α:=[ab]_f`・`β:=[a]_f∘[b]_f` が同じ次数1係数(`ab`)・同じ関数等式
(`f` と可換)を満たすことに適用して得る。 -/
theorem LubinTateAction_comp {A : Type*} [CommRing A] [IsLocalRing A] [IsDomain A]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField A) pp] [Fintype (IsLocalRing.ResidueField A)]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField A) = pp ^ ff)
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π}) (hπne0 : π ≠ 0)
    (f : PowerSeries A) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue A) f = PowerSeries.X ^ (pp ^ ff)) (a b : A) :
    LubinTateAction hq hπmax f hf0 hf1 hf (a * b) =
      PowerSeries.subst (LubinTateAction hq hπmax f hf0 hf1 hf b) (LubinTateAction hq hπmax f hf0 hf1 hf a) := by
  have hf0c : PowerSeries.constantCoeff f = 0 := by
    rw [← PowerSeries.coeff_zero_eq_constantCoeff_apply]; exact hf0
  have ha0 : PowerSeries.constantCoeff (LubinTateAction hq hπmax f hf0 hf1 hf a) = 0 :=
    constantCoeff_LubinTateAction hq hπmax f hf0 hf1 hf a
  have hb0 : PowerSeries.constantCoeff (LubinTateAction hq hπmax f hf0 hf1 hf b) = 0 :=
    constantCoeff_LubinTateAction hq hπmax f hf0 hf1 hf b
  have hab0 : PowerSeries.constantCoeff (LubinTateAction hq hπmax f hf0 hf1 hf (a * b)) = 0 :=
    constantCoeff_LubinTateAction hq hπmax f hf0 hf1 hf (a * b)
  have hβ0 : PowerSeries.constantCoeff
      (PowerSeries.subst (LubinTateAction hq hπmax f hf0 hf1 hf b) (LubinTateAction hq hπmax f hf0 hf1 hf a)) = 0 :=
    PowerSeries.constantCoeff_subst_eq_zero hb0 (LubinTateAction hq hπmax f hf0 hf1 hf a) ha0
  apply powerSeries_uniqueness hπmax hπne0 hf0c hf1 hab0 hβ0
  · rw [coeff_one_LubinTateAction hq hπmax f hf0 hf1 hf (a * b)]
    rw [coeff_one_subst_1var hb0, coeff_one_LubinTateAction hq hπmax f hf0 hf1 hf a,
      coeff_one_LubinTateAction hq hπmax f hf0 hf1 hf b]
  · exact (LubinTateAction_functional_equation hq hπmax f hf0 hf1 hf (a * b)).symm
  · calc PowerSeries.subst f
          (PowerSeries.subst (LubinTateAction hq hπmax f hf0 hf1 hf b) (LubinTateAction hq hπmax f hf0 hf1 hf a))
        = PowerSeries.subst (PowerSeries.subst f (LubinTateAction hq hπmax f hf0 hf1 hf b))
            (LubinTateAction hq hπmax f hf0 hf1 hf a) :=
          subst_comp_subst_single_gen f (LubinTateAction hq hπmax f hf0 hf1 hf b)
            (LubinTateAction hq hπmax f hf0 hf1 hf a) hf0c hb0
      _ = PowerSeries.subst (PowerSeries.subst (LubinTateAction hq hπmax f hf0 hf1 hf b) f)
            (LubinTateAction hq hπmax f hf0 hf1 hf a) := by
          rw [← LubinTateAction_functional_equation hq hπmax f hf0 hf1 hf b]
      _ = PowerSeries.subst (LubinTateAction hq hπmax f hf0 hf1 hf b)
            (PowerSeries.subst f (LubinTateAction hq hπmax f hf0 hf1 hf a)) :=
          (subst_comp_subst_single_gen (LubinTateAction hq hπmax f hf0 hf1 hf b) f
            (LubinTateAction hq hπmax f hf0 hf1 hf a) hb0 hf0c).symm
      _ = PowerSeries.subst (LubinTateAction hq hπmax f hf0 hf1 hf b)
            (PowerSeries.subst (LubinTateAction hq hπmax f hf0 hf1 hf a) f) := by
          rw [← LubinTateAction_functional_equation hq hπmax f hf0 hf1 hf a]
      _ = PowerSeries.subst
            (PowerSeries.subst (LubinTateAction hq hπmax f hf0 hf1 hf b) (LubinTateAction hq hπmax f hf0 hf1 hf a)) f :=
          subst_comp_subst_single_gen (LubinTateAction hq hπmax f hf0 hf1 hf b)
            (LubinTateAction hq hπmax f hf0 hf1 hf a) f hb0 ha0

end ABC3.Found.PGC
