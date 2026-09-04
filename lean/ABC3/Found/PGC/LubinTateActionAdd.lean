import ABC3.Found.PGC.LubinTateActionMul

/-!
# `𝒪_K` 作用への拡張: `a↦[a]_f` の加法との両立性(`sorry` 無し)

`Found/PGC/LubinTateActionMul.lean::LubinTateAction_comp`(乗法側)に続き、
`a↦[a]_f` が環準同型であることの**加法側**——
`[a+b]_f = F_f([a]_f,[b]_f)`——を確立する。これで節目(2b)
(`a↦[a]_f : 𝒪_K → End(F_f)` が環準同型であること)が全体として完成する。

## 証明の筋

`α:=[a+b]_f`・`δ:=F_f([a]_f,[b]_f)`(`family:=fun i=>if i=0 then [a]_f
else [b]_f` を `F_f` へ代入したもの、`Unit` 添字の値なので `δ` 自身
1変数の冪級数)を構成し、両者が同じ次数1の係数(いずれも `a+b`)・
同じ関数等式(`f` と可換)を満たすことを示し、**1変数の**一意性補題
(`powerSeries_uniqueness`、乗法側と同じもの——`δ` が2変数の値ではなく
`Unit` 添字、すなわち正真正銘1変数の冪級数だから)で `α=δ` を結論する。

鍵になったのは新規の一般補題 `subst_value_comp_family_general`
(「1つの値の代入」と「代入の族」の合成則のもう半分——
`LubinTateActionAssociativity.lean::subst_family_comp_value_general` が
`c(v(p))=p(c(v))` の形なのに対し、こちらは `v(c(Φ))=Φ(fun i=>v(c_i))`
の形、`MvPowerSeries.subst_comp_subst_apply` を `PowerSeries.subst_def`
で両側から橋渡しするだけで出る)。これと `subst_preserves_functional_
equation`(結合律・節目(2)で確立済み)・各成分が `f` と可換であること
(`LubinTateAction_functional_equation`)を組み合わせて `δ` の関数等式
を導く。次数1の係数は `coeff_single_subst_degree_one`(結合律で確立済み)
を `j:=()` で適用して計算する。
-/

namespace ABC3.Found.PGC

/-! ### 部品0: 「1つの値の代入」と「代入の族」の合成則(もう半分) -/

/-- **`subst_family_comp_value_general` のもう半分**: `v` を1つの値として
代入する操作は、`family` 越しの代入と可換に振る舞う——
`v.subst(family.subst(Φ)) = Φ.subst(fun i=>v.subst(family i))`。
`PowerSeries.subst_def` で両側から `MvPowerSeries.subst_comp_subst_apply`
に橋渡しする(`subst_comp_subst_single_gen` と同じ手法)。 -/
theorem subst_value_comp_family_general {A σ : Type*} [CommRing A] [Finite σ]
    {family : σ → PowerSeries A} (hfamilyHS : MvPowerSeries.HasSubst family)
    {v : PowerSeries A} (hv0 : PowerSeries.constantCoeff v = 0)
    (Φ : MvPowerSeries σ A) :
    PowerSeries.subst v (MvPowerSeries.subst family Φ) =
      MvPowerSeries.subst (fun i => PowerSeries.subst v (family i)) Φ := by
  have hvHS : MvPowerSeries.HasSubst (fun _ : Unit => v) := hasSubst_const_general hv0
  rw [PowerSeries.subst_def v (MvPowerSeries.subst family Φ),
    MvPowerSeries.subst_comp_subst_apply hfamilyHS hvHS]
  rfl

/-! ### 部品1: `family:=fun i=>if i=0 then [a]_f else [b]_f` の代入可能性 -/

/-- `family` は代入可能(`HasSubst`)——各成分の定数項が0であることから。 -/
theorem hasSubst_actionFam2 {A : Type*} [CommRing A] [IsLocalRing A] [IsDomain A]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField A) pp] [Fintype (IsLocalRing.ResidueField A)]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField A) = pp ^ ff)
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π})
    (f : PowerSeries A) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue A) f = PowerSeries.X ^ (pp ^ ff)) (a b : A) :
    MvPowerSeries.HasSubst (fun i : Fin 2 => if i = 0 then LubinTateAction hq hπmax f hf0 hf1 hf a
      else LubinTateAction hq hπmax f hf0 hf1 hf b) := by
  constructor
  · intro i
    show IsNilpotent (MvPowerSeries.constantCoeff (if i = 0 then LubinTateAction hq hπmax f hf0 hf1 hf a
      else LubinTateAction hq hπmax f hf0 hf1 hf b))
    split
    · show IsNilpotent (PowerSeries.constantCoeff (LubinTateAction hq hπmax f hf0 hf1 hf a))
      rw [constantCoeff_LubinTateAction]; exact IsNilpotent.zero
    · show IsNilpotent (PowerSeries.constantCoeff (LubinTateAction hq hπmax f hf0 hf1 hf b))
      rw [constantCoeff_LubinTateAction]; exact IsNilpotent.zero
  · intro d; exact Set.toFinite _

/-! ### 部品2: `δ:=F_f([a]_f,[b]_f)` の定数項・次数1係数 -/

theorem constantCoeff_actionAdd_delta {A : Type*} [CommRing A] [IsLocalRing A] [IsDomain A]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField A) pp] [Fintype (IsLocalRing.ResidueField A)]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField A) = pp ^ ff)
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π})
    (f : PowerSeries A) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue A) f = PowerSeries.X ^ (pp ^ ff)) (a b : A) :
    PowerSeries.constantCoeff
      (MvPowerSeries.subst (fun i : Fin 2 => if i = 0 then LubinTateAction hq hπmax f hf0 hf1 hf a
        else LubinTateAction hq hπmax f hf0 hf1 hf b) (formalGroupLaw hq hπmax f hf0 hf1 hf)) = 0 := by
  apply MvPowerSeries.constantCoeff_subst_eq_zero (hasSubst_actionFam2 hq hπmax f hf0 hf1 hf a b)
  · intro i
    show MvPowerSeries.constantCoeff (if i = 0 then LubinTateAction hq hπmax f hf0 hf1 hf a
      else LubinTateAction hq hπmax f hf0 hf1 hf b) = 0
    split
    · exact constantCoeff_LubinTateAction hq hπmax f hf0 hf1 hf a
    · exact constantCoeff_LubinTateAction hq hπmax f hf0 hf1 hf b
  · exact constantCoeff_formalGroupLaw hq hπmax f hf0 hf1 hf

/-- `δ` の次数1の係数は `a+b`——`coeff_single_subst_degree_one` を `j:=()`
で適用し、`F_f` 自身の次数1係数(いずれも `1`)と `family` の次数1係数
(`a`・`b`)の双線型な和として計算する。 -/
theorem coeff_one_actionAdd_delta {A : Type*} [CommRing A] [IsLocalRing A] [IsDomain A]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField A) pp] [Fintype (IsLocalRing.ResidueField A)]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField A) = pp ^ ff)
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π})
    (f : PowerSeries A) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue A) f = PowerSeries.X ^ (pp ^ ff)) (a b : A) :
    PowerSeries.coeff 1
      (MvPowerSeries.subst (fun i : Fin 2 => if i = 0 then LubinTateAction hq hπmax f hf0 hf1 hf a
        else LubinTateAction hq hπmax f hf0 hf1 hf b) (formalGroupLaw hq hπmax f hf0 hf1 hf)) = a + b := by
  set family : Fin 2 → PowerSeries A := fun i => if i = 0 then LubinTateAction hq hπmax f hf0 hf1 hf a
    else LubinTateAction hq hπmax f hf0 hf1 hf b with hfamily_def
  have hfamilyHS : MvPowerSeries.HasSubst family := hasSubst_actionFam2 hq hπmax f hf0 hf1 hf a b
  have ha_order : ∀ i : Fin 2, (1 : ℕ∞) ≤ MvPowerSeries.order (family i) := by
    intro i
    apply MvPowerSeries.one_le_order_iff_constCoeff_eq_zero.mpr
    rw [hfamily_def]
    show PowerSeries.constantCoeff (if i = 0 then LubinTateAction hq hπmax f hf0 hf1 hf a
      else LubinTateAction hq hπmax f hf0 hf1 hf b) = 0
    split
    · exact constantCoeff_LubinTateAction hq hπmax f hf0 hf1 hf a
    · exact constantCoeff_LubinTateAction hq hπmax f hf0 hf1 hf b
  have hΦ0 : MvPowerSeries.constantCoeff (formalGroupLaw hq hπmax f hf0 hf1 hf) = 0 :=
    constantCoeff_formalGroupLaw hq hπmax f hf0 hf1 hf
  show MvPowerSeries.coeff (Finsupp.single () 1)
    (MvPowerSeries.subst family (formalGroupLaw hq hπmax f hf0 hf1 hf)) = a + b
  rw [coeff_single_subst_degree_one hΦ0 ha_order hfamilyHS ()]
  rw [Fin.sum_univ_two]
  rw [coeff_single01_formalGroupLaw hq hπmax f hf0 hf1 hf, coeff_single10_formalGroupLaw hq hπmax f hf0 hf1 hf]
  show (1 : A) * MvPowerSeries.coeff (Finsupp.single () 1) (family 0) +
      1 * MvPowerSeries.coeff (Finsupp.single () 1) (family 1) = a + b
  rw [hfamily_def]
  show (1 : A) * PowerSeries.coeff 1 (LubinTateAction hq hπmax f hf0 hf1 hf a)
    + 1 * PowerSeries.coeff 1 (LubinTateAction hq hπmax f hf0 hf1 hf b) = a + b
  rw [coeff_one_LubinTateAction hq hπmax f hf0 hf1 hf a, coeff_one_LubinTateAction hq hπmax f hf0 hf1 hf b]
  ring

/-! ### 部品3: `δ` の関数等式(`f` と可換) -/

/-- `δ:=F_f([a]_f,[b]_f)` は `f` と可換——`subst_preserves_functional_
equation`(`F_f` 自身の関数等式)と `subst_value_comp_family_general`
(部品0)・各成分の可換性(`LubinTateAction_functional_equation`)を
組み合わせる。 -/
theorem actionAdd_delta_functional_equation {A : Type*} [CommRing A] [IsLocalRing A] [IsDomain A]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField A) pp] [Fintype (IsLocalRing.ResidueField A)]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField A) = pp ^ ff)
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π})
    (f : PowerSeries A) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue A) f = PowerSeries.X ^ (pp ^ ff)) (a b : A) :
    PowerSeries.subst f
      (MvPowerSeries.subst (fun i : Fin 2 => if i = 0 then LubinTateAction hq hπmax f hf0 hf1 hf a
        else LubinTateAction hq hπmax f hf0 hf1 hf b) (formalGroupLaw hq hπmax f hf0 hf1 hf)) =
    PowerSeries.subst
      (MvPowerSeries.subst (fun i : Fin 2 => if i = 0 then LubinTateAction hq hπmax f hf0 hf1 hf a
        else LubinTateAction hq hπmax f hf0 hf1 hf b) (formalGroupLaw hq hπmax f hf0 hf1 hf)) f := by
  set family : Fin 2 → PowerSeries A := fun i => if i = 0 then LubinTateAction hq hπmax f hf0 hf1 hf a
    else LubinTateAction hq hπmax f hf0 hf1 hf b with hfamily_def
  have hf0c : PowerSeries.constantCoeff f = 0 := by
    rw [← PowerSeries.coeff_zero_eq_constantCoeff_apply]; exact hf0
  have hfamilyHS : MvPowerSeries.HasSubst family := hasSubst_actionFam2 hq hπmax f hf0 hf1 hf a b
  have hΦ0 : MvPowerSeries.constantCoeff (formalGroupLaw hq hπmax f hf0 hf1 hf) = 0 :=
    constantCoeff_formalGroupLaw hq hπmax f hf0 hf1 hf
  have hkey := subst_preserves_functional_equation hf0 hΦ0
    (formalGroupLaw_f_isEndomorphism hq hπmax f hf0 hf1 hf) hfamilyHS
  calc PowerSeries.subst f (MvPowerSeries.subst family (formalGroupLaw hq hπmax f hf0 hf1 hf))
      = MvPowerSeries.subst (fun i => PowerSeries.subst f (family i)) (formalGroupLaw hq hπmax f hf0 hf1 hf) :=
        subst_value_comp_family_general hfamilyHS hf0c (formalGroupLaw hq hπmax f hf0 hf1 hf)
    _ = MvPowerSeries.subst (fun i => PowerSeries.subst (family i) f) (formalGroupLaw hq hπmax f hf0 hf1 hf) := by
        congr 1
        funext i
        by_cases hi : i = 0
        · subst hi
          show PowerSeries.subst f (LubinTateAction hq hπmax f hf0 hf1 hf a) =
            PowerSeries.subst (LubinTateAction hq hπmax f hf0 hf1 hf a) f
          exact (LubinTateAction_functional_equation hq hπmax f hf0 hf1 hf a).symm
        · have hi1 : i = 1 := by omega
          subst hi1
          show PowerSeries.subst f (LubinTateAction hq hπmax f hf0 hf1 hf b) =
            PowerSeries.subst (LubinTateAction hq hπmax f hf0 hf1 hf b) f
          exact (LubinTateAction_functional_equation hq hπmax f hf0 hf1 hf b).symm
    _ = PowerSeries.subst (MvPowerSeries.subst family (formalGroupLaw hq hπmax f hf0 hf1 hf)) f := hkey

/-! ### 部品4: 結論 -/

/-- ★★★★★★★**`a↦[a]_f` の加法との両立性**: `[a+b]_f = F_f([a]_f,[b]_f)`。
1変数の一意性補題(`powerSeries_uniqueness`——`δ` は `Unit` 添字、
すなわち正真正銘1変数の冪級数なので、乗法側と同じくこちらを使う)を、
`α:=[a+b]_f`・`δ:=F_f([a]_f,[b]_f)` が同じ次数1係数(`a+b`)・同じ
関数等式(`f` と可換)を満たすことに適用して得る。これで `a↦[a]_f`
が環準同型であること(節目2b)が乗法側(`LubinTateAction_comp`)と
合わせて全体として完成する。 -/
theorem LubinTateAction_add {A : Type*} [CommRing A] [IsLocalRing A] [IsDomain A]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField A) pp] [Fintype (IsLocalRing.ResidueField A)]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField A) = pp ^ ff)
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π}) (hπne0 : π ≠ 0)
    (f : PowerSeries A) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue A) f = PowerSeries.X ^ (pp ^ ff)) (a b : A) :
    LubinTateAction hq hπmax f hf0 hf1 hf (a + b) =
      MvPowerSeries.subst (fun i : Fin 2 => if i = 0 then LubinTateAction hq hπmax f hf0 hf1 hf a
        else LubinTateAction hq hπmax f hf0 hf1 hf b) (formalGroupLaw hq hπmax f hf0 hf1 hf) := by
  have hf0c : PowerSeries.constantCoeff f = 0 := by
    rw [← PowerSeries.coeff_zero_eq_constantCoeff_apply]; exact hf0
  apply powerSeries_uniqueness hπmax hπne0 hf0c hf1
    (constantCoeff_LubinTateAction hq hπmax f hf0 hf1 hf (a + b))
    (constantCoeff_actionAdd_delta hq hπmax f hf0 hf1 hf a b)
  · rw [coeff_one_LubinTateAction hq hπmax f hf0 hf1 hf (a + b),
      coeff_one_actionAdd_delta hq hπmax f hf0 hf1 hf a b]
  · exact (LubinTateAction_functional_equation hq hπmax f hf0 hf1 hf (a + b)).symm
  · exact actionAdd_delta_functional_equation hq hπmax f hf0 hf1 hf a b

end ABC3.Found.PGC
