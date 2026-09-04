import ABC3.Found.PGC.LubinTateGeneralUniqueness

/-!
# Lubin-Tate 形式群法則: 結合律 `F_f(F_f(X,Y),Z)=F_f(X,F_f(Y,Z))`

`Found/PGC/LubinTateGeneralUniqueness.lean::mvPowerSeries_uniqueness_general`
(任意の有限添字型、`Fin 3` を含む)を実際に3変数へ適用して、古典的な
Lubin-Tate 形式群法則の最後の性質——結合律——を確立する。

## 方針

`G := F_f(F_f(X,Y),Z)`・`H := F_f(X,F_f(Y,Z))` はどちらも `MvPowerSeries
(Fin 3) A`。この2つが同じ次数1の係数(いずれも全て `1`——`F_f` の次数1
部分が `X+Y` と対称だから)・同じ関数等式(`F_f` 自身の関数等式を
2回合成して得る)を満たすことを示し、`mvPowerSeries_uniqueness_general`
で `G=H` を結論する。

鍵になった一般補題は `subst_preserves_functional_equation`——「`Φ` が
自身の関数等式 `Φ.subst(aFam)=f.subst(Φ)`(`aFam:=`各変数に `f` を
適用する族)を満たすならば、任意の代入族 `b` を経由しても同じ形の
関数等式 `Φ.subst(b).subst(...) = ...` を保つ」という
`swap_preserves_functional_equation`(可換律)の一般化。これと
`coeff_single_subst_degree_one`(次数1の係数は代入前後で線型に振る舞う
という一般公式)の2つの新しい一般補題だけで、`G`・`H` それぞれの
関数等式・次数1係数がほぼ機械的に出た。
-/

namespace ABC3.Found.PGC

/-! ### 部品0: 代入の次数1係数についての一般公式 -/

/-- `d.prod` の次数下界: `a` の各成分の次数が `≥1` ならば、`d.prod(a^·)`
の次数は `Finsupp.degree d` 以上——`order_prod_pow_sub_prod_pow_ge_finset`
の「差」を取らない単純版。 -/
theorem order_finsuppProd_ge_degree {A σ τ : Type*} [CommRing A] {a : τ → MvPowerSeries σ A}
    (ha_order : ∀ i, 1 ≤ (a i).order) (d : τ →₀ ℕ) :
    ((Finsupp.degree d : ℕ) : ℕ∞) ≤ (d.prod fun s m => (a s) ^ m).order := by
  have hprod_eq : d.prod (fun s m => (a s) ^ m) = ∏ i ∈ d.support, (a i) ^ (d i) := rfl
  rw [hprod_eq, Finsupp.degree_apply]
  rw [Nat.cast_sum]
  calc ∑ i ∈ d.support, (d i : ℕ∞) ≤ ∑ i ∈ d.support, (a i ^ d i).order := by
        apply Finset.sum_le_sum
        intro i _
        calc (d i : ℕ∞) = d i • (1 : ℕ∞) := by simp
          _ ≤ d i • (a i).order := by gcongr; exact ha_order i
          _ ≤ (a i ^ d i).order := MvPowerSeries.le_order_pow (d i)
    _ ≤ (∏ i ∈ d.support, (a i) ^ (d i)).order := MvPowerSeries.le_order_prod _ d.support

/-- ★★★★★★**代入の次数1係数の一般公式**: `Φ` の定数項が0、代入族 `a` の
各成分の次数が `≥1` のとき、`subst a Φ` の次数1(添字 `j`)の係数は、
`Φ` の次数1の係数と `a` の次数1の係数の「双線型な」和 `∑_i Φ_i • a_i(j)`
に一致する——次数 `≥2` の項は代入後の次数を押し上げるので効かない、
という標準的な次数勘定。 -/
theorem coeff_single_subst_degree_one {A σ τ : Type*} [CommRing A] [Fintype τ] [DecidableEq τ]
    {Φ : MvPowerSeries τ A} (hΦ0 : MvPowerSeries.constantCoeff Φ = 0)
    {a : τ → MvPowerSeries σ A} (ha_order : ∀ i, 1 ≤ (a i).order)
    (haHS : MvPowerSeries.HasSubst a) (j : σ) :
    MvPowerSeries.coeff (Finsupp.single j 1) (MvPowerSeries.subst a Φ) =
      ∑ i : τ, MvPowerSeries.coeff (Finsupp.single i 1) Φ * MvPowerSeries.coeff (Finsupp.single j 1) (a i) := by
  rw [MvPowerSeries.coeff_subst haHS]
  rw [finsum_eq_sum_of_support_subset _ (s := Finset.univ.image (fun i : τ => Finsupp.single i (1 : ℕ)))
    (fun d hd => by
      simp only [Function.mem_support] at hd
      simp only [Finset.coe_image, Finset.coe_univ, Set.image_univ, Set.mem_range]
      by_contra hcon
      apply hd
      by_cases hd0 : d = 0
      · rw [hd0]
        show MvPowerSeries.coeff 0 Φ • MvPowerSeries.coeff (Finsupp.single j 1) (1 : MvPowerSeries σ A) = 0
        rw [← MvPowerSeries.coeff_zero_eq_constantCoeff_apply] at hΦ0
        rw [hΦ0]; simp
      · have hdeg_not_le1 : ¬ Finsupp.degree d ≤ 1 := by
          intro hle
          rcases finsupp_degree_le_one_cases d hle with h0 | ⟨i, hi⟩
          · exact hd0 h0
          · exact hcon ⟨i, hi.symm⟩
        have hdeg2 : 2 ≤ Finsupp.degree d := by omega
        have horder := order_finsuppProd_ge_degree ha_order d
        have hlt : ((Finsupp.degree (Finsupp.single j (1 : ℕ)) : ℕ) : ℕ∞) <
            (d.prod fun s m => (a s) ^ m).order := by
          calc ((Finsupp.degree (Finsupp.single j (1 : ℕ)) : ℕ) : ℕ∞) = 1 := by
                simp [Finsupp.degree_single]
            _ < ((2 : ℕ) : ℕ∞) := by exact_mod_cast (by norm_num : (1 : ℕ) < 2)
            _ ≤ ((Finsupp.degree d : ℕ) : ℕ∞) := by exact_mod_cast hdeg2
            _ ≤ _ := horder
        have hz := MvPowerSeries.coeff_of_lt_order hlt
        rw [hz]; simp)]
  rw [Finset.sum_image (fun i _ i' _ hii' => Finsupp.single_left_injective one_ne_zero hii')]
  apply Finset.sum_congr rfl
  intro i _
  show MvPowerSeries.coeff (Finsupp.single i 1) Φ •
      MvPowerSeries.coeff (Finsupp.single j 1) ((Finsupp.single i (1 : ℕ)).prod fun s m => (a s) ^ m) = _
  rw [Finsupp.prod_single_index (by simp)]
  simp

/-! ### 部品1: 代入越しに関数等式を保つ(可換律の `swap_preserves_functional_equation` の一般化) -/

/-- `hasSubst_const`(`Fin 2` 固定)の一般化。 -/
theorem hasSubst_const_general {A : Type*} [CommRing A] {σ τ : Type*} [Finite σ] {v : MvPowerSeries τ A}
    (hv0 : MvPowerSeries.constantCoeff v = 0) :
    MvPowerSeries.HasSubst (fun _ : σ => v) := by
  constructor
  · intro _; show IsNilpotent (MvPowerSeries.constantCoeff v); rw [hv0]; exact IsNilpotent.zero
  · intro d; exact Set.toFinite _

/-- `subst_family_comp_value`(`Fin 2` 固定)の一般化。 -/
theorem subst_family_comp_value_general {A τ σ : Type*} [CommRing A] [Finite σ]
    {c : σ → MvPowerSeries τ A} (hc : MvPowerSeries.HasSubst c) {v : MvPowerSeries σ A}
    (hv0 : MvPowerSeries.constantCoeff v = 0) (p : PowerSeries A) :
    MvPowerSeries.subst c (PowerSeries.subst v p) = PowerSeries.subst (MvPowerSeries.subst c v) p := by
  rw [PowerSeries.subst_def v p, MvPowerSeries.subst_comp_subst_apply (hasSubst_const_general hv0) hc,
    ← PowerSeries.subst_def]

/-- ★★★★★★★★★★**関数等式は代入越しに保たれる**(`swap_preserves_functional_
equation` の一般化、`swap` という特別な代入だけでなく**任意の**代入族
`b` について成り立つ): `Φ` が自身の関数等式 `Φ.subst(aFam)=f.subst(Φ)`
(`aFam:=`各変数に `f` を適用する族)を満たすならば、`Φ.subst(fun i=>
f.subst(b i)) = f.subst(Φ.subst(b))`。3段の `subst_comp_subst_apply`/
`subst_family_comp_value_general` の合成で閉じる——`swap_preserves_
functional_equation` の証明と全く同じ構造。 -/
theorem subst_preserves_functional_equation {A τ σ : Type*} [CommRing A] [Finite τ] [Finite σ]
    {f : PowerSeries A} (hf0 : PowerSeries.coeff 0 f = 0)
    {Φ : MvPowerSeries τ A} (hΦ0 : MvPowerSeries.constantCoeff Φ = 0)
    (heqΦ : MvPowerSeries.subst (fun i => PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries τ A) f) Φ =
      PowerSeries.subst Φ f)
    {b : τ → MvPowerSeries σ A} (hbHS : MvPowerSeries.HasSubst b) :
    MvPowerSeries.subst (fun i => PowerSeries.subst (b i) f) Φ =
      PowerSeries.subst (MvPowerSeries.subst b Φ) f := by
  have hf0c : PowerSeries.constantCoeff f = 0 := by
    rw [← PowerSeries.coeff_zero_eq_constantCoeff_apply]; exact hf0
  have hHSf : PowerSeries.HasSubst f := by
    show IsNilpotent (PowerSeries.constantCoeff f); rw [hf0c]; exact IsNilpotent.zero
  have hHSa : MvPowerSeries.HasSubst (fun i => PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries τ A) f) :=
    hasSubst_g_subst_X f hf0c
  set a : τ → MvPowerSeries τ A := fun i => PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries τ A) f with ha_def
  have hstep1 : MvPowerSeries.subst (fun i => PowerSeries.subst (b i) f) Φ =
      MvPowerSeries.subst b (MvPowerSeries.subst a Φ) := by
    rw [MvPowerSeries.subst_comp_subst_apply hHSa hbHS]
    congr 1
    funext i
    rw [ha_def]
    show PowerSeries.subst (b i) f = MvPowerSeries.subst b (PowerSeries.subst (MvPowerSeries.X i) f)
    rw [subst_family_comp_value_general hbHS (MvPowerSeries.constantCoeff_X i) f, MvPowerSeries.subst_X hbHS i]
  have hstep2 : MvPowerSeries.subst b (PowerSeries.subst Φ f) = PowerSeries.subst (MvPowerSeries.subst b Φ) f :=
    subst_family_comp_value_general hbHS hΦ0 f
  rw [hstep1, heqΦ, hstep2]

/-! ### 部品2: `F_f(X,Y)`・`F_f(Y,Z)` を3変数の世界へ埋め込む -/

/-- `X_0,X_1` を使う代入族——`F_f(X,Y)` を3変数(`X,Y,Z`)の世界に
埋め込むためのもの。 -/
noncomputable def substXY {A : Type*} [CommRing A] : Fin 2 → MvPowerSeries (Fin 3) A :=
  fun i => MvPowerSeries.X (if i = 0 then (0 : Fin 3) else 1)

/-- `X_1,X_2` を使う代入族——`F_f(Y,Z)` を3変数の世界に埋め込むためのもの。 -/
noncomputable def substYZ {A : Type*} [CommRing A] : Fin 2 → MvPowerSeries (Fin 3) A :=
  fun i => MvPowerSeries.X (if i = 0 then (1 : Fin 3) else 2)

theorem hasSubst_substXY {A : Type*} [CommRing A] : MvPowerSeries.HasSubst (substXY (A := A)) := by
  constructor
  · intro i
    show IsNilpotent (MvPowerSeries.constantCoeff
      (MvPowerSeries.X (if i = 0 then (0 : Fin 3) else 1) : MvPowerSeries (Fin 3) A))
    rw [MvPowerSeries.constantCoeff_X]; exact IsNilpotent.zero
  · intro d; exact Set.toFinite _

theorem hasSubst_substYZ {A : Type*} [CommRing A] : MvPowerSeries.HasSubst (substYZ (A := A)) := by
  constructor
  · intro i
    show IsNilpotent (MvPowerSeries.constantCoeff
      (MvPowerSeries.X (if i = 0 then (1 : Fin 3) else 2) : MvPowerSeries (Fin 3) A))
    rw [MvPowerSeries.constantCoeff_X]; exact IsNilpotent.zero
  · intro d; exact Set.toFinite _

variable {A : Type*} [CommRing A] [IsLocalRing A] [IsDomain A]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField A) pp] [Fintype (IsLocalRing.ResidueField A)]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField A) = pp ^ ff)
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π})
    (f : PowerSeries A) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hfres : PowerSeries.map (IsLocalRing.residue A) f = PowerSeries.X ^ (pp ^ ff))

include hq hπmax hf0 hf1 hfres in
/-- `F_f(X,Y)` を3変数の世界に埋め込んだもの。 -/
noncomputable def Ffxy : MvPowerSeries (Fin 3) A :=
  MvPowerSeries.subst substXY (formalGroupLaw hq hπmax f hf0 hf1 hfres)

include hq hπmax hf0 hf1 hfres in
/-- `F_f(Y,Z)` を3変数の世界に埋め込んだもの。 -/
noncomputable def Ffyz : MvPowerSeries (Fin 3) A :=
  MvPowerSeries.subst substYZ (formalGroupLaw hq hπmax f hf0 hf1 hfres)

include hq hπmax hf0 hf1 hfres in
theorem constantCoeff_Ffxy : MvPowerSeries.constantCoeff (Ffxy hq hπmax f hf0 hf1 hfres) = 0 := by
  apply MvPowerSeries.constantCoeff_subst_eq_zero hasSubst_substXY
  · intro i
    exact MvPowerSeries.constantCoeff_X _
  · exact constantCoeff_formalGroupLaw hq hπmax f hf0 hf1 hfres

include hq hπmax hf0 hf1 hfres in
theorem constantCoeff_Ffyz : MvPowerSeries.constantCoeff (Ffyz hq hπmax f hf0 hf1 hfres) = 0 := by
  apply MvPowerSeries.constantCoeff_subst_eq_zero hasSubst_substYZ
  · intro i
    exact MvPowerSeries.constantCoeff_X _
  · exact constantCoeff_formalGroupLaw hq hπmax f hf0 hf1 hfres

include hq hπmax hf0 hf1 hfres in
/-- `Ffxy` も(3変数の`f`族に対して)`F_f` と同じ形の関数等式を満たす
——`F_f` の自己準同型性を `substXY` 越しに引き継ぐだけ
(`subst_preserves_functional_equation`)。 -/
theorem Ffxy_functional_equation :
    MvPowerSeries.subst (fun i => PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries (Fin 3) A) f)
        (Ffxy hq hπmax f hf0 hf1 hfres) =
      PowerSeries.subst (Ffxy hq hπmax f hf0 hf1 hfres) f := by
  have hf0c : PowerSeries.constantCoeff f = 0 := by
    rw [← PowerSeries.coeff_zero_eq_constantCoeff_apply]; exact hf0
  have hHSfFam3 : MvPowerSeries.HasSubst
      (fun i => PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries (Fin 3) A) f) := hasSubst_g_subst_X f hf0c
  have hkey := subst_preserves_functional_equation hf0
    (constantCoeff_formalGroupLaw hq hπmax f hf0 hf1 hfres)
    (formalGroupLaw_f_isEndomorphism hq hπmax f hf0 hf1 hfres) hasSubst_substXY
  show MvPowerSeries.subst (fun i => PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries (Fin 3) A) f)
      (MvPowerSeries.subst substXY (formalGroupLaw hq hπmax f hf0 hf1 hfres)) = _
  rw [MvPowerSeries.subst_comp_subst_apply hasSubst_substXY hHSfFam3]
  have hfam_eq : (fun s : Fin 2 => MvPowerSeries.subst
      (fun i => PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries (Fin 3) A) f) ((substXY (A := A)) s)) =
      (fun i : Fin 2 => PowerSeries.subst ((substXY (A := A)) i) f) := by
    funext i
    show MvPowerSeries.subst (fun j => PowerSeries.subst (MvPowerSeries.X j : MvPowerSeries (Fin 3) A) f)
      (MvPowerSeries.X (if i = 0 then (0 : Fin 3) else 1) : MvPowerSeries (Fin 3) A)
      = PowerSeries.subst (MvPowerSeries.X (if i = 0 then (0 : Fin 3) else 1) : MvPowerSeries (Fin 3) A) f
    rw [MvPowerSeries.subst_X hHSfFam3]
  rw [hfam_eq, hkey]
  rfl

include hq hπmax hf0 hf1 hfres in
/-- `Ffyz` も同じ関数等式を満たす(`Ffxy_functional_equation` の対)。 -/
theorem Ffyz_functional_equation :
    MvPowerSeries.subst (fun i => PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries (Fin 3) A) f)
        (Ffyz hq hπmax f hf0 hf1 hfres) =
      PowerSeries.subst (Ffyz hq hπmax f hf0 hf1 hfres) f := by
  have hf0c : PowerSeries.constantCoeff f = 0 := by
    rw [← PowerSeries.coeff_zero_eq_constantCoeff_apply]; exact hf0
  have hHSfFam3 : MvPowerSeries.HasSubst
      (fun i => PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries (Fin 3) A) f) := hasSubst_g_subst_X f hf0c
  have hkey := subst_preserves_functional_equation hf0
    (constantCoeff_formalGroupLaw hq hπmax f hf0 hf1 hfres)
    (formalGroupLaw_f_isEndomorphism hq hπmax f hf0 hf1 hfres) hasSubst_substYZ
  show MvPowerSeries.subst (fun i => PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries (Fin 3) A) f)
      (MvPowerSeries.subst substYZ (formalGroupLaw hq hπmax f hf0 hf1 hfres)) = _
  rw [MvPowerSeries.subst_comp_subst_apply hasSubst_substYZ hHSfFam3]
  have hfam_eq : (fun s : Fin 2 => MvPowerSeries.subst
      (fun i => PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries (Fin 3) A) f) ((substYZ (A := A)) s)) =
      (fun i : Fin 2 => PowerSeries.subst ((substYZ (A := A)) i) f) := by
    funext i
    show MvPowerSeries.subst (fun j => PowerSeries.subst (MvPowerSeries.X j : MvPowerSeries (Fin 3) A) f)
      (MvPowerSeries.X (if i = 0 then (1 : Fin 3) else 2) : MvPowerSeries (Fin 3) A)
      = PowerSeries.subst (MvPowerSeries.X (if i = 0 then (1 : Fin 3) else 2) : MvPowerSeries (Fin 3) A) f
    rw [MvPowerSeries.subst_X hHSfFam3]
  rw [hfam_eq, hkey]
  rfl

/-! ### 部品3: `G:=F_f(F_f(X,Y),Z)`・`H:=F_f(X,F_f(Y,Z))` -/

include hq hπmax hf0 hf1 hfres in
/-- `G := F_f(F_f(X,Y),Z)`。 -/
noncomputable def assocG : MvPowerSeries (Fin 3) A :=
  MvPowerSeries.subst (fun i : Fin 2 => if i = 0 then Ffxy hq hπmax f hf0 hf1 hfres
      else (MvPowerSeries.X 2 : MvPowerSeries (Fin 3) A))
    (formalGroupLaw hq hπmax f hf0 hf1 hfres)

include hq hπmax hf0 hf1 hfres in
/-- `H := F_f(X,F_f(Y,Z))`。 -/
noncomputable def assocH : MvPowerSeries (Fin 3) A :=
  MvPowerSeries.subst (fun i : Fin 2 => if i = 0 then (MvPowerSeries.X 0 : MvPowerSeries (Fin 3) A)
      else Ffyz hq hπmax f hf0 hf1 hfres)
    (formalGroupLaw hq hπmax f hf0 hf1 hfres)

include hq hπmax hf0 hf1 hfres in
theorem hasSubst_gFam2 : MvPowerSeries.HasSubst (fun i : Fin 2 => if i = 0 then Ffxy hq hπmax f hf0 hf1 hfres
      else (MvPowerSeries.X 2 : MvPowerSeries (Fin 3) A)) := by
  constructor
  · intro i
    show IsNilpotent (MvPowerSeries.constantCoeff (if i = 0 then Ffxy hq hπmax f hf0 hf1 hfres
      else (MvPowerSeries.X 2 : MvPowerSeries (Fin 3) A)))
    split
    · rw [constantCoeff_Ffxy hq hπmax f hf0 hf1 hfres]; exact IsNilpotent.zero
    · rw [MvPowerSeries.constantCoeff_X]; exact IsNilpotent.zero
  · intro d; exact Set.toFinite _

include hq hπmax hf0 hf1 hfres in
theorem hasSubst_hFam2 : MvPowerSeries.HasSubst (fun i : Fin 2 => if i = 0 then (MvPowerSeries.X 0 : MvPowerSeries (Fin 3) A)
      else Ffyz hq hπmax f hf0 hf1 hfres) := by
  constructor
  · intro i
    show IsNilpotent (MvPowerSeries.constantCoeff (if i = 0 then (MvPowerSeries.X 0 : MvPowerSeries (Fin 3) A)
      else Ffyz hq hπmax f hf0 hf1 hfres))
    split
    · rw [MvPowerSeries.constantCoeff_X]; exact IsNilpotent.zero
    · rw [constantCoeff_Ffyz hq hπmax f hf0 hf1 hfres]; exact IsNilpotent.zero
  · intro d; exact Set.toFinite _

include hq hπmax hf0 hf1 hfres in
theorem constantCoeff_assocG : MvPowerSeries.constantCoeff (assocG hq hπmax f hf0 hf1 hfres) = 0 := by
  apply MvPowerSeries.constantCoeff_subst_eq_zero (hasSubst_gFam2 hq hπmax f hf0 hf1 hfres)
  · intro i
    show MvPowerSeries.constantCoeff (if i = 0 then Ffxy hq hπmax f hf0 hf1 hfres
      else (MvPowerSeries.X 2 : MvPowerSeries (Fin 3) A)) = 0
    split
    · exact constantCoeff_Ffxy hq hπmax f hf0 hf1 hfres
    · exact MvPowerSeries.constantCoeff_X _
  · exact constantCoeff_formalGroupLaw hq hπmax f hf0 hf1 hfres

include hq hπmax hf0 hf1 hfres in
theorem constantCoeff_assocH : MvPowerSeries.constantCoeff (assocH hq hπmax f hf0 hf1 hfres) = 0 := by
  apply MvPowerSeries.constantCoeff_subst_eq_zero (hasSubst_hFam2 hq hπmax f hf0 hf1 hfres)
  · intro i
    show MvPowerSeries.constantCoeff (if i = 0 then (MvPowerSeries.X 0 : MvPowerSeries (Fin 3) A)
      else Ffyz hq hπmax f hf0 hf1 hfres) = 0
    split
    · exact MvPowerSeries.constantCoeff_X _
    · exact constantCoeff_Ffyz hq hπmax f hf0 hf1 hfres
  · exact constantCoeff_formalGroupLaw hq hπmax f hf0 hf1 hfres

include hq hπmax hf0 hf1 hfres in
/-- `G` は `Ffxy` の関数等式(`Ffxy_functional_equation`)を経由して
`F_f` 自身の関数等式を2回合成したものを満たす。 -/
theorem assocG_functional_equation :
    MvPowerSeries.subst (fun i => PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries (Fin 3) A) f)
        (assocG hq hπmax f hf0 hf1 hfres) =
      PowerSeries.subst (assocG hq hπmax f hf0 hf1 hfres) f := by
  have hf0c : PowerSeries.constantCoeff f = 0 := by
    rw [← PowerSeries.coeff_zero_eq_constantCoeff_apply]; exact hf0
  have hHSfFam3 : MvPowerSeries.HasSubst
      (fun i => PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries (Fin 3) A) f) := hasSubst_g_subst_X f hf0c
  have hkey := subst_preserves_functional_equation hf0
    (constantCoeff_formalGroupLaw hq hπmax f hf0 hf1 hfres)
    (formalGroupLaw_f_isEndomorphism hq hπmax f hf0 hf1 hfres) (hasSubst_gFam2 hq hπmax f hf0 hf1 hfres)
  show MvPowerSeries.subst (fun i => PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries (Fin 3) A) f)
      (MvPowerSeries.subst (fun i : Fin 2 => if i = 0 then Ffxy hq hπmax f hf0 hf1 hfres
        else (MvPowerSeries.X 2 : MvPowerSeries (Fin 3) A)) (formalGroupLaw hq hπmax f hf0 hf1 hfres)) = _
  rw [MvPowerSeries.subst_comp_subst_apply (hasSubst_gFam2 hq hπmax f hf0 hf1 hfres) hHSfFam3]
  have hfam_eq : (fun s : Fin 2 => MvPowerSeries.subst
      (fun i => PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries (Fin 3) A) f)
      (if s = 0 then Ffxy hq hπmax f hf0 hf1 hfres else (MvPowerSeries.X 2 : MvPowerSeries (Fin 3) A))) =
      (fun i : Fin 2 => PowerSeries.subst
        (if i = 0 then Ffxy hq hπmax f hf0 hf1 hfres else (MvPowerSeries.X 2 : MvPowerSeries (Fin 3) A)) f) := by
    funext i
    split
    · exact Ffxy_functional_equation hq hπmax f hf0 hf1 hfres
    · exact MvPowerSeries.subst_X hHSfFam3 2
  rw [hfam_eq, hkey]
  rfl

include hq hπmax hf0 hf1 hfres in
/-- `H` も同じ関数等式を満たす(`assocG_functional_equation` の対)。 -/
theorem assocH_functional_equation :
    MvPowerSeries.subst (fun i => PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries (Fin 3) A) f)
        (assocH hq hπmax f hf0 hf1 hfres) =
      PowerSeries.subst (assocH hq hπmax f hf0 hf1 hfres) f := by
  have hf0c : PowerSeries.constantCoeff f = 0 := by
    rw [← PowerSeries.coeff_zero_eq_constantCoeff_apply]; exact hf0
  have hHSfFam3 : MvPowerSeries.HasSubst
      (fun i => PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries (Fin 3) A) f) := hasSubst_g_subst_X f hf0c
  have hkey := subst_preserves_functional_equation hf0
    (constantCoeff_formalGroupLaw hq hπmax f hf0 hf1 hfres)
    (formalGroupLaw_f_isEndomorphism hq hπmax f hf0 hf1 hfres) (hasSubst_hFam2 hq hπmax f hf0 hf1 hfres)
  show MvPowerSeries.subst (fun i => PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries (Fin 3) A) f)
      (MvPowerSeries.subst (fun i : Fin 2 => if i = 0 then (MvPowerSeries.X 0 : MvPowerSeries (Fin 3) A)
        else Ffyz hq hπmax f hf0 hf1 hfres) (formalGroupLaw hq hπmax f hf0 hf1 hfres)) = _
  rw [MvPowerSeries.subst_comp_subst_apply (hasSubst_hFam2 hq hπmax f hf0 hf1 hfres) hHSfFam3]
  have hfam_eq : (fun s : Fin 2 => MvPowerSeries.subst
      (fun i => PowerSeries.subst (MvPowerSeries.X i : MvPowerSeries (Fin 3) A) f)
      (if s = 0 then (MvPowerSeries.X 0 : MvPowerSeries (Fin 3) A) else Ffyz hq hπmax f hf0 hf1 hfres)) =
      (fun i : Fin 2 => PowerSeries.subst
        (if i = 0 then (MvPowerSeries.X 0 : MvPowerSeries (Fin 3) A) else Ffyz hq hπmax f hf0 hf1 hfres) f) := by
    funext i
    split
    · exact MvPowerSeries.subst_X hHSfFam3 0
    · exact Ffyz_functional_equation hq hπmax f hf0 hf1 hfres
  rw [hfam_eq, hkey]
  rfl

/-! ### 部品4: `G`・`H` の次数1の係数(いずれも全て `1`) -/

include hq hπmax hf0 hf1 hfres in
theorem coeff_single_Ffxy (j : Fin 3) :
    MvPowerSeries.coeff (Finsupp.single j 1) (Ffxy hq hπmax f hf0 hf1 hfres) =
      MvPowerSeries.coeff (Finsupp.single j 1) (substXY (A := A) 0) +
        MvPowerSeries.coeff (Finsupp.single j 1) (substXY (A := A) 1) := by
  have hord : ∀ i : Fin 2, (1 : ℕ∞) ≤ ((substXY (A := A)) i).order := by
    intro i
    apply MvPowerSeries.one_le_order_iff_constCoeff_eq_zero.mpr
    exact MvPowerSeries.constantCoeff_X _
  show MvPowerSeries.coeff (Finsupp.single j 1)
    (MvPowerSeries.subst substXY (formalGroupLaw hq hπmax f hf0 hf1 hfres)) = _
  rw [coeff_single_subst_degree_one (constantCoeff_formalGroupLaw hq hπmax f hf0 hf1 hfres) hord
    hasSubst_substXY j]
  rw [Fin.sum_univ_two]
  rw [coeff_single01_formalGroupLaw hq hπmax f hf0 hf1 hfres,
    coeff_single10_formalGroupLaw hq hπmax f hf0 hf1 hfres]
  ring

include hq hπmax hf0 hf1 hfres in
theorem coeff_single_Ffyz (j : Fin 3) :
    MvPowerSeries.coeff (Finsupp.single j 1) (Ffyz hq hπmax f hf0 hf1 hfres) =
      MvPowerSeries.coeff (Finsupp.single j 1) (substYZ (A := A) 0) +
        MvPowerSeries.coeff (Finsupp.single j 1) (substYZ (A := A) 1) := by
  have hord : ∀ i : Fin 2, (1 : ℕ∞) ≤ ((substYZ (A := A)) i).order := by
    intro i
    apply MvPowerSeries.one_le_order_iff_constCoeff_eq_zero.mpr
    exact MvPowerSeries.constantCoeff_X _
  show MvPowerSeries.coeff (Finsupp.single j 1)
    (MvPowerSeries.subst substYZ (formalGroupLaw hq hπmax f hf0 hf1 hfres)) = _
  rw [coeff_single_subst_degree_one (constantCoeff_formalGroupLaw hq hπmax f hf0 hf1 hfres) hord
    hasSubst_substYZ j]
  rw [Fin.sum_univ_two]
  rw [coeff_single01_formalGroupLaw hq hπmax f hf0 hf1 hfres,
    coeff_single10_formalGroupLaw hq hπmax f hf0 hf1 hfres]
  ring

include hq hπmax hf0 hf1 hfres in
theorem coeff_single_assocG_eq_one (j : Fin 3) :
    MvPowerSeries.coeff (Finsupp.single j 1) (assocG hq hπmax f hf0 hf1 hfres) = 1 := by
  have hord : ∀ i : Fin 2, (1 : ℕ∞) ≤ ((fun i : Fin 2 => if i = 0 then Ffxy hq hπmax f hf0 hf1 hfres
      else (MvPowerSeries.X 2 : MvPowerSeries (Fin 3) A)) i).order := by
    intro i
    apply MvPowerSeries.one_le_order_iff_constCoeff_eq_zero.mpr
    show MvPowerSeries.constantCoeff (if i = 0 then Ffxy hq hπmax f hf0 hf1 hfres
      else (MvPowerSeries.X 2 : MvPowerSeries (Fin 3) A)) = 0
    split
    · exact constantCoeff_Ffxy hq hπmax f hf0 hf1 hfres
    · exact MvPowerSeries.constantCoeff_X _
  show MvPowerSeries.coeff (Finsupp.single j 1)
    (MvPowerSeries.subst (fun i : Fin 2 => if i = 0 then Ffxy hq hπmax f hf0 hf1 hfres
      else (MvPowerSeries.X 2 : MvPowerSeries (Fin 3) A)) (formalGroupLaw hq hπmax f hf0 hf1 hfres)) = 1
  rw [coeff_single_subst_degree_one (constantCoeff_formalGroupLaw hq hπmax f hf0 hf1 hfres) hord
    (hasSubst_gFam2 hq hπmax f hf0 hf1 hfres) j]
  rw [Fin.sum_univ_two, if_pos rfl, if_neg (by decide : (1 : Fin 2) ≠ 0)]
  rw [coeff_single01_formalGroupLaw hq hπmax f hf0 hf1 hfres,
    coeff_single10_formalGroupLaw hq hπmax f hf0 hf1 hfres, one_mul, one_mul]
  rw [coeff_single_Ffxy hq hπmax f hf0 hf1 hfres j]
  show (MvPowerSeries.coeff (Finsupp.single j 1)) (MvPowerSeries.X (0 : Fin 3)) +
      (MvPowerSeries.coeff (Finsupp.single j 1)) (MvPowerSeries.X (1 : Fin 3)) +
      (MvPowerSeries.coeff (Finsupp.single j 1)) (MvPowerSeries.X (2 : Fin 3)) = 1
  simp only [MvPowerSeries.coeff_X, Finsupp.single_left_inj (one_ne_zero (α := ℕ))]
  fin_cases j <;> simp

include hq hπmax hf0 hf1 hfres in
theorem coeff_single_assocH_eq_one (j : Fin 3) :
    MvPowerSeries.coeff (Finsupp.single j 1) (assocH hq hπmax f hf0 hf1 hfres) = 1 := by
  have hord : ∀ i : Fin 2, (1 : ℕ∞) ≤ ((fun i : Fin 2 => if i = 0 then (MvPowerSeries.X 0 : MvPowerSeries (Fin 3) A)
      else Ffyz hq hπmax f hf0 hf1 hfres) i).order := by
    intro i
    apply MvPowerSeries.one_le_order_iff_constCoeff_eq_zero.mpr
    show MvPowerSeries.constantCoeff (if i = 0 then (MvPowerSeries.X 0 : MvPowerSeries (Fin 3) A)
      else Ffyz hq hπmax f hf0 hf1 hfres) = 0
    split
    · exact MvPowerSeries.constantCoeff_X _
    · exact constantCoeff_Ffyz hq hπmax f hf0 hf1 hfres
  show MvPowerSeries.coeff (Finsupp.single j 1)
    (MvPowerSeries.subst (fun i : Fin 2 => if i = 0 then (MvPowerSeries.X 0 : MvPowerSeries (Fin 3) A)
      else Ffyz hq hπmax f hf0 hf1 hfres) (formalGroupLaw hq hπmax f hf0 hf1 hfres)) = 1
  rw [coeff_single_subst_degree_one (constantCoeff_formalGroupLaw hq hπmax f hf0 hf1 hfres) hord
    (hasSubst_hFam2 hq hπmax f hf0 hf1 hfres) j]
  rw [Fin.sum_univ_two, if_pos rfl, if_neg (by decide : (1 : Fin 2) ≠ 0)]
  rw [coeff_single01_formalGroupLaw hq hπmax f hf0 hf1 hfres,
    coeff_single10_formalGroupLaw hq hπmax f hf0 hf1 hfres, one_mul, one_mul]
  rw [coeff_single_Ffyz hq hπmax f hf0 hf1 hfres j]
  show (MvPowerSeries.coeff (Finsupp.single j 1)) (MvPowerSeries.X (0 : Fin 3)) +
      ((MvPowerSeries.coeff (Finsupp.single j 1)) (MvPowerSeries.X (1 : Fin 3)) +
      (MvPowerSeries.coeff (Finsupp.single j 1)) (MvPowerSeries.X (2 : Fin 3))) = 1
  simp only [MvPowerSeries.coeff_X, Finsupp.single_left_inj (one_ne_zero (α := ℕ))]
  fin_cases j <;> simp

include hq hπmax hf0 hf1 hfres in
/-- ★★★★★★★★★★★★★★★★**Lubin-Tate 形式群法則の結合律**:
`F_f(F_f(X,Y),Z) = F_f(X,F_f(Y,Z))`。`G:=F_f(F_f(X,Y),Z)`・
`H:=F_f(X,F_f(Y,Z))` が同じ次数1係数(いずれも全て `1`)・同じ関数等式
(`assocG_functional_equation`・`assocH_functional_equation`)を満たす
ので、3変数の一意性補題(`mvPowerSeries_uniqueness_general`)で
`G=H` が出る。これで古典的な Lubin-Tate 形式群法則の3性質
(単位元則2つ・可換律・結合律)が全て sorry 無しで揃った。 -/
theorem formalGroupLaw_associative (hπne0 : π ≠ 0) :
    assocG hq hπmax f hf0 hf1 hfres = assocH hq hπmax f hf0 hf1 hfres := by
  apply mvPowerSeries_uniqueness_general hπmax hπne0 hf0 hf1
    (constantCoeff_assocG hq hπmax f hf0 hf1 hfres) (constantCoeff_assocH hq hπmax f hf0 hf1 hfres)
  · intro j
    rw [coeff_single_assocG_eq_one hq hπmax f hf0 hf1 hfres j,
      coeff_single_assocH_eq_one hq hπmax f hf0 hf1 hfres j]
  · exact assocG_functional_equation hq hπmax f hf0 hf1 hfres
  · exact assocH_functional_equation hq hπmax f hf0 hf1 hfres

end ABC3.Found.PGC
