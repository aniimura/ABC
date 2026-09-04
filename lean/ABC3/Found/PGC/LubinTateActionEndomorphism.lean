import ABC3.Found.PGC.LubinTateEndoLimit
import ABC3.Found.PGC.LubinTateAssociativity

/-!
# `𝒪_K` 作用への拡張・節目(2): `[a]_f` が `F_f` の自己準同型であること(`sorry` 無し)

`ResearchPaper/blocked-leaves.json` の `progress2026_09_04q`〜`s`・§9-1346 で
分解した節目(2)——`F_f([a]X,[a]Y)=[a]F_f(X,Y)`(`[a]_f` が形式群 `F_f` の
自己準同型であること)——を確立する。§9-1340 で見立てた通り、これは
`Found/PGC/LubinTateEndoLimit.lean::LubinTateEndo_functional_equation`
(`g:=f` に特殊化したものを `LubinTateAction` と呼ぶ)と、可換律で確立した
2変数の一意性補題 `mvPowerSeries_uniqueness` から**既存の道具だけ**で
形式的に出る——新しい次数ごとの構成は一切不要だった。

## 証明の筋

`Ψ_1 := F_f([a]X,[a]Y)`・`Ψ_2 := [a]F_f(X,Y)` を構成し、両者が
(a) 同じ次数1の係数(いずれも全て `a`)、(b) 同じ関数等式
(`F_f` 自身の関数等式)を満たすことを示し、一意性補題で `Ψ_1=Ψ_2` を
結論する。鍵になったのは:

- `subst_comp_subst_single_gen`: 1変数どうしの代入の合成規則
  `c.subst(b.subst(f)) = (c.subst(b)).subst(f)`(`MvPowerSeries.subst_
  comp_subst_apply` を `PowerSeries.subst_def` で橋渡しして得る、
  `Unit` 添字での特殊化)。
- これと `LubinTateAction` 自身の関数等式(`f∘[a]=[a]∘f`)を組み合わせて、
  「`X_i` の点で評価した可換性」`key_commute_at_Xi` を導き、
  `Ψ_1`・`Ψ_2` それぞれの関数等式を `subst_preserves_functional_equation`
  (可換律・結合律で確立済み)と組み合わせて得る。
- 次数1の係数は `coeff_single_subst_degree_one`(結合律で確立済み)を
  2段(`aFam2 i` 自身の次数1係数、`F_f` の次数1係数)に適用して計算する。
-/

namespace ABC3.Found.PGC

/-! ### 部品0: 1変数どうしの代入の合成規則 -/

/-- `c.subst(b.subst(f)) = (c.subst(b)).subst(f)`——1変数の`PowerSeries.subst`
どうしの合成規則(`c` の行き先は一般の `MvPowerSeries τ R` でよい)。
`PowerSeries.subst_def` で `MvPowerSeries.subst_comp_subst_apply` に橋渡しする。 -/
theorem subst_comp_subst_single_gen {R τ : Type*} [CommRing R] [Finite τ]
    (c : MvPowerSeries τ R) (b f : PowerSeries R)
    (hc0 : MvPowerSeries.constantCoeff c = 0) (hb0 : PowerSeries.constantCoeff b = 0) :
    PowerSeries.subst c (PowerSeries.subst b f) = PowerSeries.subst (PowerSeries.subst c b) f := by
  have hcHS : MvPowerSeries.HasSubst (fun _ : Unit => c) := hasSubst_const_general hc0
  have hbHS : MvPowerSeries.HasSubst (fun _ : Unit => b) := hasSubst_const_general hb0
  rw [PowerSeries.subst_def c, PowerSeries.subst_def b, MvPowerSeries.subst_comp_subst_apply hbHS hcHS]
  rw [PowerSeries.subst_def]
  congr 1

/-! ### 部品1: `[a]_f`(`LubinTateEndo` の `g:=f` 特殊化) -/

/-- `[a]_f := LubinTateEndo`(`g:=f`)。`f` 自身と合成として可換
(`f∘[a]_f=[a]_f∘f`)な、次数1の係数が `a` の1変数冪級数。 -/
noncomputable def LubinTateAction {A : Type*} [CommRing A] [IsLocalRing A] [IsDomain A]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField A) pp] [Fintype (IsLocalRing.ResidueField A)]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField A) = pp ^ ff)
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π})
    (f : PowerSeries A) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue A) f = PowerSeries.X ^ (pp ^ ff)) (a : A) :
    PowerSeries A :=
  LubinTateEndo hq hπmax f (PowerSeries.coeff_zero_eq_constantCoeff_apply f ▸ hf0) hf1 hf f hf0 hf1 hf a

theorem LubinTateAction_functional_equation {A : Type*} [CommRing A] [IsLocalRing A] [IsDomain A]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField A) pp] [Fintype (IsLocalRing.ResidueField A)]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField A) = pp ^ ff)
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π})
    (f : PowerSeries A) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue A) f = PowerSeries.X ^ (pp ^ ff)) (a : A) :
    PowerSeries.subst (LubinTateAction hq hπmax f hf0 hf1 hf a) f =
      PowerSeries.subst f (LubinTateAction hq hπmax f hf0 hf1 hf a) :=
  LubinTateEndo_functional_equation hq hπmax f (PowerSeries.coeff_zero_eq_constantCoeff_apply f ▸ hf0) hf1 hf
    f hf0 hf1 hf a

theorem constantCoeff_LubinTateAction {A : Type*} [CommRing A] [IsLocalRing A] [IsDomain A]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField A) pp] [Fintype (IsLocalRing.ResidueField A)]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField A) = pp ^ ff)
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π})
    (f : PowerSeries A) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue A) f = PowerSeries.X ^ (pp ^ ff)) (a : A) :
    PowerSeries.constantCoeff (LubinTateAction hq hπmax f hf0 hf1 hf a) = 0 :=
  constantCoeff_LubinTateEndo hq hπmax f (PowerSeries.coeff_zero_eq_constantCoeff_apply f ▸ hf0) hf1 hf
    f hf0 hf1 hf a

/-- `LubinTateAction` の次数1の係数は目標どおり `a`。 -/
theorem coeff_one_LubinTateAction {A : Type*} [CommRing A] [IsLocalRing A] [IsDomain A]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField A) pp] [Fintype (IsLocalRing.ResidueField A)]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField A) = pp ^ ff)
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π})
    (f : PowerSeries A) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue A) f = PowerSeries.X ^ (pp ^ ff)) (a : A) :
    PowerSeries.coeff 1 (LubinTateAction hq hπmax f hf0 hf1 hf a) = a :=
  coeff_one_LubinTateEndo hq hπmax f (PowerSeries.coeff_zero_eq_constantCoeff_apply f ▸ hf0) hf1 hf f hf0 hf1 hf a

/-! ### 部品2: `X_i` の点で評価した可換性 -/

/-- `f∘[a]_f=[a]_f∘f`(`LubinTateAction_functional_equation`)を `X_i` の点
(`MvPowerSeries (Fin 2) A` の値)で評価したもの: `f([a](X_i))=[a](f(X_i))`。 -/
theorem key_commute_at_Xi {A : Type*} [CommRing A] [IsLocalRing A] [IsDomain A]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField A) pp] [Fintype (IsLocalRing.ResidueField A)]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField A) = pp ^ ff)
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π})
    (f : PowerSeries A) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue A) f = PowerSeries.X ^ (pp ^ ff)) (a : A) (i : Fin 2) :
    PowerSeries.subst (PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries (Fin 2) A)
      (LubinTateAction hq hπmax f hf0 hf1 hf a)) f =
      PowerSeries.subst (PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries (Fin 2) A) f)
        (LubinTateAction hq hπmax f hf0 hf1 hf a) := by
  have heq := LubinTateAction_functional_equation hq hπmax f hf0 hf1 hf a
  have hc0 : MvPowerSeries.constantCoeff (MvPowerSeries.X i : MvPowerSeries (Fin 2) A) = 0 :=
    MvPowerSeries.constantCoeff_X i
  have h1 := subst_comp_subst_single_gen (MvPowerSeries.X i : MvPowerSeries (Fin 2) A)
    (LubinTateAction hq hπmax f hf0 hf1 hf a) f hc0 (constantCoeff_LubinTateAction hq hπmax f hf0 hf1 hf a)
  have hf0c : PowerSeries.constantCoeff f = 0 := by
    rw [← PowerSeries.coeff_zero_eq_constantCoeff_apply]; exact hf0
  have h2 := subst_comp_subst_single_gen (MvPowerSeries.X i : MvPowerSeries (Fin 2) A)
    f (LubinTateAction hq hπmax f hf0 hf1 hf a) hc0 hf0c
  rw [← h1, ← h2, heq]

/-- `subst(fFam2)(aFam2 i) = subst(aFam2 i) f`——チェインルール
(`subst_family_comp_value_general`)と `key_commute_at_Xi` を合わせたもの。 -/
theorem subst_fFam2_aFam2_eq {A : Type*} [CommRing A] [IsLocalRing A] [IsDomain A]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField A) pp] [Fintype (IsLocalRing.ResidueField A)]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField A) = pp ^ ff)
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π})
    (f : PowerSeries A) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue A) f = PowerSeries.X ^ (pp ^ ff)) (a : A) (i : Fin 2) :
    MvPowerSeries.subst (fun j : Fin 2 => PowerSeries.subst (MvPowerSeries.X j : MvPowerSeries (Fin 2) A) f)
        (PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries (Fin 2) A) (LubinTateAction hq hπmax f hf0 hf1 hf a)) =
      PowerSeries.subst (PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries (Fin 2) A)
        (LubinTateAction hq hπmax f hf0 hf1 hf a)) f := by
  have hf0c : PowerSeries.constantCoeff f = 0 := by
    rw [← PowerSeries.coeff_zero_eq_constantCoeff_apply]; exact hf0
  have hfFam2HS : MvPowerSeries.HasSubst
      (fun j : Fin 2 => PowerSeries.subst (MvPowerSeries.X j : MvPowerSeries (Fin 2) A) f) :=
    hasSubst_g_subst_X f hf0c
  rw [subst_family_comp_value_general hfFam2HS (MvPowerSeries.constantCoeff_X i)
    (LubinTateAction hq hπmax f hf0 hf1 hf a)]
  have hXi : MvPowerSeries.subst (fun j : Fin 2 => PowerSeries.subst (MvPowerSeries.X j : MvPowerSeries (Fin 2) A) f)
      (MvPowerSeries.X i : MvPowerSeries (Fin 2) A) =
      PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries (Fin 2) A) f := MvPowerSeries.subst_X hfFam2HS i
  rw [hXi]
  exact (key_commute_at_Xi hq hπmax f hf0 hf1 hf a i).symm

/-! ### 部品3: `Ψ_1:=F_f([a]X,[a]Y)`・`Ψ_2:=[a]F_f(X,Y)` の関数等式 -/

/-- `Ψ_1 := F_f([a]X,[a]Y)` は `F_f` 自身と同じ関数等式を満たす。 -/
theorem Psi1_functional_equation {A : Type*} [CommRing A] [IsLocalRing A] [IsDomain A]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField A) pp] [Fintype (IsLocalRing.ResidueField A)]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField A) = pp ^ ff)
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π})
    (f : PowerSeries A) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue A) f = PowerSeries.X ^ (pp ^ ff)) (a : A) :
    MvPowerSeries.subst (fun i : Fin 2 => PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries (Fin 2) A) f)
        (MvPowerSeries.subst (fun i : Fin 2 => PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries (Fin 2) A)
          (LubinTateAction hq hπmax f hf0 hf1 hf a)) (formalGroupLaw hq hπmax f hf0 hf1 hf)) =
      PowerSeries.subst (MvPowerSeries.subst (fun i : Fin 2 => PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries (Fin 2) A)
          (LubinTateAction hq hπmax f hf0 hf1 hf a)) (formalGroupLaw hq hπmax f hf0 hf1 hf)) f := by
  set aFam2 : Fin 2 → MvPowerSeries (Fin 2) A :=
    fun i => PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries (Fin 2) A) (LubinTateAction hq hπmax f hf0 hf1 hf a)
    with haFam2_def
  have hf0c : PowerSeries.constantCoeff f = 0 := by
    rw [← PowerSeries.coeff_zero_eq_constantCoeff_apply]; exact hf0
  have haFam2HS : MvPowerSeries.HasSubst aFam2 := hasSubst_g_subst_X (LubinTateAction hq hπmax f hf0 hf1 hf a)
    (constantCoeff_LubinTateAction hq hπmax f hf0 hf1 hf a)
  have hfFam2HS : MvPowerSeries.HasSubst
      (fun i : Fin 2 => PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries (Fin 2) A) f) :=
    hasSubst_g_subst_X f hf0c
  have hkey := subst_preserves_functional_equation hf0
    (constantCoeff_formalGroupLaw hq hπmax f hf0 hf1 hf)
    (formalGroupLaw_f_isEndomorphism hq hπmax f hf0 hf1 hf) haFam2HS
  rw [MvPowerSeries.subst_comp_subst_apply haFam2HS hfFam2HS]
  rw [show (fun s : Fin 2 => MvPowerSeries.subst
        (fun i : Fin 2 => PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries (Fin 2) A) f) (aFam2 s)) =
      (fun i : Fin 2 => PowerSeries.subst (aFam2 i) f) from
    funext (fun i => subst_fFam2_aFam2_eq hq hπmax f hf0 hf1 hf a i)]
  exact hkey

/-- `Ψ_2 := [a]F_f(X,Y)` も `F_f` 自身と同じ関数等式を満たす。 -/
theorem Psi2_functional_equation {A : Type*} [CommRing A] [IsLocalRing A] [IsDomain A]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField A) pp] [Fintype (IsLocalRing.ResidueField A)]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField A) = pp ^ ff)
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π})
    (f : PowerSeries A) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue A) f = PowerSeries.X ^ (pp ^ ff)) (a : A) :
    MvPowerSeries.subst (fun i : Fin 2 => PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries (Fin 2) A) f)
        (PowerSeries.subst (formalGroupLaw hq hπmax f hf0 hf1 hf) (LubinTateAction hq hπmax f hf0 hf1 hf a)) =
      PowerSeries.subst (PowerSeries.subst (formalGroupLaw hq hπmax f hf0 hf1 hf)
        (LubinTateAction hq hπmax f hf0 hf1 hf a)) f := by
  have hf0c : PowerSeries.constantCoeff f = 0 := by
    rw [← PowerSeries.coeff_zero_eq_constantCoeff_apply]; exact hf0
  have hfFam2HS : MvPowerSeries.HasSubst
      (fun i : Fin 2 => PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries (Fin 2) A) f) :=
    hasSubst_g_subst_X f hf0c
  have hFf0 : MvPowerSeries.constantCoeff (formalGroupLaw hq hπmax f hf0 hf1 hf) = 0 :=
    constantCoeff_formalGroupLaw hq hπmax f hf0 hf1 hf
  have ha0 : PowerSeries.constantCoeff (LubinTateAction hq hπmax f hf0 hf1 hf a) = 0 :=
    constantCoeff_LubinTateAction hq hπmax f hf0 hf1 hf a
  calc MvPowerSeries.subst (fun i : Fin 2 => PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries (Fin 2) A) f)
        (PowerSeries.subst (formalGroupLaw hq hπmax f hf0 hf1 hf) (LubinTateAction hq hπmax f hf0 hf1 hf a))
      = PowerSeries.subst (MvPowerSeries.subst
          (fun i : Fin 2 => PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries (Fin 2) A) f)
          (formalGroupLaw hq hπmax f hf0 hf1 hf)) (LubinTateAction hq hπmax f hf0 hf1 hf a) :=
        subst_family_comp_value_general hfFam2HS hFf0 (LubinTateAction hq hπmax f hf0 hf1 hf a)
    _ = PowerSeries.subst (PowerSeries.subst (formalGroupLaw hq hπmax f hf0 hf1 hf) f)
          (LubinTateAction hq hπmax f hf0 hf1 hf a) := by
        rw [formalGroupLaw_f_isEndomorphism hq hπmax f hf0 hf1 hf]
    _ = PowerSeries.subst (formalGroupLaw hq hπmax f hf0 hf1 hf)
          (PowerSeries.subst f (LubinTateAction hq hπmax f hf0 hf1 hf a)) :=
        (subst_comp_subst_single_gen (formalGroupLaw hq hπmax f hf0 hf1 hf) f
          (LubinTateAction hq hπmax f hf0 hf1 hf a) hFf0 hf0c).symm
    _ = PowerSeries.subst (formalGroupLaw hq hπmax f hf0 hf1 hf)
          (PowerSeries.subst (LubinTateAction hq hπmax f hf0 hf1 hf a) f) := by
        rw [LubinTateAction_functional_equation hq hπmax f hf0 hf1 hf a]
    _ = PowerSeries.subst (PowerSeries.subst (formalGroupLaw hq hπmax f hf0 hf1 hf)
          (LubinTateAction hq hπmax f hf0 hf1 hf a)) f :=
        subst_comp_subst_single_gen (formalGroupLaw hq hπmax f hf0 hf1 hf)
          (LubinTateAction hq hπmax f hf0 hf1 hf a) f hFf0 ha0

/-! ### 部品4: `Ψ_1`・`Ψ_2` の次数1の係数(いずれも `a`) -/

theorem coeff_single_Psi1_eq_a {A : Type*} [CommRing A] [IsLocalRing A] [IsDomain A]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField A) pp] [Fintype (IsLocalRing.ResidueField A)]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField A) = pp ^ ff)
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π})
    (f : PowerSeries A) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue A) f = PowerSeries.X ^ (pp ^ ff)) (a : A) (j : Fin 2) :
    MvPowerSeries.coeff (Finsupp.single j 1)
      (MvPowerSeries.subst (fun i : Fin 2 => PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries (Fin 2) A)
        (LubinTateAction hq hπmax f hf0 hf1 hf a)) (formalGroupLaw hq hπmax f hf0 hf1 hf)) = a := by
  set aFam2 : Fin 2 → MvPowerSeries (Fin 2) A :=
    fun i => PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries (Fin 2) A) (LubinTateAction hq hπmax f hf0 hf1 hf a)
    with haFam2_def
  have ha0 : PowerSeries.constantCoeff (LubinTateAction hq hπmax f hf0 hf1 hf a) = 0 :=
    constantCoeff_LubinTateAction hq hπmax f hf0 hf1 hf a
  have haiorder : ∀ i : Fin 2, (1 : ℕ∞) ≤ MvPowerSeries.order (aFam2 i) := by
    intro i
    apply MvPowerSeries.one_le_order_iff_constCoeff_eq_zero.mpr
    show MvPowerSeries.constantCoeff (PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries (Fin 2) A)
      (LubinTateAction hq hπmax f hf0 hf1 hf a)) = 0
    exact PowerSeries.constantCoeff_subst_eq_zero (MvPowerSeries.constantCoeff_X i)
      (LubinTateAction hq hπmax f hf0 hf1 hf a) ha0
  have haHS : MvPowerSeries.HasSubst aFam2 := hasSubst_g_subst_X (LubinTateAction hq hπmax f hf0 hf1 hf a) ha0
  rw [coeff_single_subst_degree_one (constantCoeff_formalGroupLaw hq hπmax f hf0 hf1 hf) haiorder haHS j]
  rw [Fin.sum_univ_two, coeff_single01_formalGroupLaw hq hπmax f hf0 hf1 hf,
    coeff_single10_formalGroupLaw hq hπmax f hf0 hf1 hf, one_mul, one_mul]
  have hcoeff_aFam2 : ∀ i : Fin 2, MvPowerSeries.coeff (Finsupp.single j 1) (aFam2 i) =
      if j = i then a else 0 := by
    intro i
    show MvPowerSeries.coeff (Finsupp.single j 1)
      (PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries (Fin 2) A) (LubinTateAction hq hπmax f hf0 hf1 hf a)) = _
    have hXi_order : (1 : ℕ∞) ≤ MvPowerSeries.order (MvPowerSeries.X i : MvPowerSeries (Fin 2) A) :=
      MvPowerSeries.one_le_order_iff_constCoeff_eq_zero.mpr (MvPowerSeries.constantCoeff_X i)
    have hthis := coeff_single_subst_degree_one (σ := Fin 2) (τ := Unit)
      (Φ := LubinTateAction hq hπmax f hf0 hf1 hf a) ha0
      (a := fun _ : Unit => (MvPowerSeries.X i : MvPowerSeries (Fin 2) A))
      (fun _ => hXi_order) (hasSubst_const_general (MvPowerSeries.constantCoeff_X i)) j
    rw [Fintype.sum_unique] at hthis
    show MvPowerSeries.coeff (Finsupp.single j 1)
      (MvPowerSeries.subst (fun _ : Unit => (MvPowerSeries.X i : MvPowerSeries (Fin 2) A))
        (LubinTateAction hq hπmax f hf0 hf1 hf a)) = _
    rw [hthis]
    have h1 : MvPowerSeries.coeff (Finsupp.single (default : Unit) 1) (LubinTateAction hq hπmax f hf0 hf1 hf a) = a := by
      show PowerSeries.coeff 1 (LubinTateAction hq hπmax f hf0 hf1 hf a) = a
      exact coeff_one_LubinTateAction hq hπmax f hf0 hf1 hf a
    rw [h1, MvPowerSeries.coeff_X]
    simp only [Finsupp.single_left_inj (one_ne_zero (α := ℕ))]
    split_ifs <;> ring
  rw [hcoeff_aFam2 0, hcoeff_aFam2 1]
  fin_cases j <;> simp

theorem coeff_single_Psi2_eq_a {A : Type*} [CommRing A] [IsLocalRing A] [IsDomain A]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField A) pp] [Fintype (IsLocalRing.ResidueField A)]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField A) = pp ^ ff)
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π})
    (f : PowerSeries A) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue A) f = PowerSeries.X ^ (pp ^ ff)) (a : A) (j : Fin 2) :
    MvPowerSeries.coeff (Finsupp.single j 1)
      (PowerSeries.subst (formalGroupLaw hq hπmax f hf0 hf1 hf) (LubinTateAction hq hπmax f hf0 hf1 hf a)) = a := by
  have ha0 : PowerSeries.constantCoeff (LubinTateAction hq hπmax f hf0 hf1 hf a) = 0 :=
    constantCoeff_LubinTateAction hq hπmax f hf0 hf1 hf a
  have hFf0 : MvPowerSeries.constantCoeff (formalGroupLaw hq hπmax f hf0 hf1 hf) = 0 :=
    constantCoeff_formalGroupLaw hq hπmax f hf0 hf1 hf
  have hFforder : (1 : ℕ∞) ≤ MvPowerSeries.order (formalGroupLaw hq hπmax f hf0 hf1 hf : MvPowerSeries (Fin 2) A) :=
    MvPowerSeries.one_le_order_iff_constCoeff_eq_zero.mpr hFf0
  have hHSFf : MvPowerSeries.HasSubst (fun _ : Unit => (formalGroupLaw hq hπmax f hf0 hf1 hf : MvPowerSeries (Fin 2) A)) :=
    hasSubst_const_general hFf0
  have hthis := coeff_single_subst_degree_one (σ := Fin 2) (τ := Unit)
    (Φ := LubinTateAction hq hπmax f hf0 hf1 hf a) ha0
    (a := fun _ : Unit => (formalGroupLaw hq hπmax f hf0 hf1 hf : MvPowerSeries (Fin 2) A))
    (fun _ => hFforder) hHSFf j
  rw [Fintype.sum_unique] at hthis
  show MvPowerSeries.coeff (Finsupp.single j 1)
    (MvPowerSeries.subst (fun _ : Unit => (formalGroupLaw hq hπmax f hf0 hf1 hf : MvPowerSeries (Fin 2) A))
      (LubinTateAction hq hπmax f hf0 hf1 hf a)) = a
  rw [hthis]
  have h1 : MvPowerSeries.coeff (Finsupp.single (default : Unit) 1) (LubinTateAction hq hπmax f hf0 hf1 hf a) = a :=
    coeff_one_LubinTateAction hq hπmax f hf0 hf1 hf a
  rw [h1]
  by_cases hj : j = 0
  · rw [hj, coeff_single01_formalGroupLaw hq hπmax f hf0 hf1 hf, mul_one]
  · have hj1 : j = 1 := by omega
    rw [hj1, coeff_single10_formalGroupLaw hq hπmax f hf0 hf1 hf, mul_one]

/-! ### ★★★★★★★★本題: 節目(2)——`[a]_f` は `F_f` の自己準同型 -/

/-- ★★★★★★★★**節目(2)完成**: `F_f([a]X,[a]Y) = [a]F_f(X,Y)`——
`LubinTateAction`(`[a]_f`)は `formalGroupLaw`(`F_f`)の自己準同型。
`Ψ_1:=F_f([a]X,[a]Y)`・`Ψ_2:=[a]F_f(X,Y)` が同じ次数1係数(全て `a`)・
同じ関数等式を満たすことから、2変数の一意性補題(可換律で確立済みの
`mvPowerSeries_uniqueness`)で結論する。これで `a↦[a]_f` が
`𝒪_K→End(F_f)` の環準同型であることの第一段(準同型であることの
必要条件)が確立した。 -/
theorem formalGroupLaw_endomorphism_of_action {A : Type*} [CommRing A] [IsLocalRing A] [IsDomain A]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField A) pp] [Fintype (IsLocalRing.ResidueField A)]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField A) = pp ^ ff)
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π}) (hπne0 : π ≠ 0)
    (f : PowerSeries A) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue A) f = PowerSeries.X ^ (pp ^ ff)) (a : A) :
    MvPowerSeries.subst (fun i : Fin 2 => PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries (Fin 2) A)
        (LubinTateAction hq hπmax f hf0 hf1 hf a)) (formalGroupLaw hq hπmax f hf0 hf1 hf) =
      PowerSeries.subst (formalGroupLaw hq hπmax f hf0 hf1 hf) (LubinTateAction hq hπmax f hf0 hf1 hf a) := by
  have ha0 : PowerSeries.constantCoeff (LubinTateAction hq hπmax f hf0 hf1 hf a) = 0 :=
    constantCoeff_LubinTateAction hq hπmax f hf0 hf1 hf a
  have hFf0 : MvPowerSeries.constantCoeff (formalGroupLaw hq hπmax f hf0 hf1 hf) = 0 :=
    constantCoeff_formalGroupLaw hq hπmax f hf0 hf1 hf
  have haHS : MvPowerSeries.HasSubst
      (fun i : Fin 2 => PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries (Fin 2) A)
        (LubinTateAction hq hπmax f hf0 hf1 hf a)) :=
    hasSubst_g_subst_X (LubinTateAction hq hπmax f hf0 hf1 hf a) ha0
  have hΨ1_0 : MvPowerSeries.constantCoeff
      (MvPowerSeries.subst (fun i : Fin 2 => PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries (Fin 2) A)
        (LubinTateAction hq hπmax f hf0 hf1 hf a)) (formalGroupLaw hq hπmax f hf0 hf1 hf)) = 0 := by
    apply MvPowerSeries.constantCoeff_subst_eq_zero haHS
    · intro i
      exact PowerSeries.constantCoeff_subst_eq_zero (MvPowerSeries.constantCoeff_X i)
        (LubinTateAction hq hπmax f hf0 hf1 hf a) ha0
    · exact hFf0
  have hΨ2_0 : MvPowerSeries.constantCoeff
      (PowerSeries.subst (formalGroupLaw hq hπmax f hf0 hf1 hf) (LubinTateAction hq hπmax f hf0 hf1 hf a)) = 0 :=
    PowerSeries.constantCoeff_subst_eq_zero hFf0 (LubinTateAction hq hπmax f hf0 hf1 hf a) ha0
  apply mvPowerSeries_uniqueness hπmax hπne0 hf0 hf1 hΨ1_0 hΨ2_0
  · intro i
    by_cases hi : i = 0
    · rw [hi, coeff_single_Psi1_eq_a hq hπmax f hf0 hf1 hf a 0, coeff_single_Psi2_eq_a hq hπmax f hf0 hf1 hf a 0]
    · have hi1 : i = 1 := by omega
      rw [hi1, coeff_single_Psi1_eq_a hq hπmax f hf0 hf1 hf a 1, coeff_single_Psi2_eq_a hq hπmax f hf0 hf1 hf a 1]
  · exact Psi1_functional_equation hq hπmax f hf0 hf1 hf a
  · exact Psi2_functional_equation hq hπmax f hf0 hf1 hf a

end ABC3.Found.PGC
