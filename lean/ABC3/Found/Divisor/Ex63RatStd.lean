/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Divisor.Ex63Types
import ABC3.Found.Divisor.ArithSuppSplit
import ABC3.Found.FrdI.Thm64RatStd
import ABC3.Found.FrdI.Prop55PfRatStd
import ABC3.Found.FrdI.Prop55UnTrRatStd

/-!
# [FrdI] Theorem 6.4, (i) —— 算術 Frobenioid は **rationally standard 型**

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.114。

原文 (FrdI p.115):
> C is of [strictly] rational type, and that every object of (Cun-tr)birat is Frobenius-

## ★★原文が 1 行で畳んだ 2 つ

原文 `Theorem 6.4, (i)` の証明は

> Also, it is immediate from the definition of B that C is of [strictly] rational type,
> and that every object of (C^un-tr)^birat is Frobenius-compact.

と書く。★この `immediate` の中身が本ファイルである。

### ★(a) strictly rational —— `hsp`

`Definition 4.5, (ii)` は「どの素点 `p` にも `a, b ∈ Φ(L)` と `y ∈ B(L)` があって
`toGp a − toGp b = ±Div_B(y)`、`p ∈ Supp a`、`p ∉ Supp b`」を要求する。

★中身は **「`v` で `|y|_v ≠ 1` となる `y ∈ L^×` を取る」**だけである:

| 素点 | 取り方 |
|---|---|
| 有限素点 `𝔭` | `exists_ordFin_eq_one`(一様化元、在庫) |
| 無限素点 `w` | ★**`y = 2`** —— `w(2) = 2 ≠ 1` |

★あとは `exists_split_suppElt_R`(`ArithSuppSplit.lean`)で正負に割るだけ。

### ★★(b) `(𝒞^un-tr)^birat` の Frobenius-compact 対象

★★★**幾何 (`Theorem 6.2, (iii)`) の入口は算術では閉じている** ——
`Φ^birat(L) = Div_B(L^×)`(主因子)は積公式で次数 0、`Φ(L)` は有効なので

    Φ(L) ∩ Φ^birat(L) = {0}

であり、在庫の `birat_isFrobeniusCompact_of_primary` が要求する
「`Φ` の非零元が `Φ^birat` に入る」は**満たしようがない**。

★逃げ道は `Thm64RatStd.lean` に置いた:
**底の自己射が恒等しかない対象**(ここでは `⊥` = `F` 自身)を選べば
`σ = Φ(Base θ⁻¹) = id` が無料で出るので、準素元の仮定は要らない。
★`Φ^birat` の無限位数の元は `Div_B(2)` を 2 本の pre-step の差として作る。

## ★本ファイルで閉じること

| 定理 | 中身 |
|---|---|
| `infinitePlace_two` | ★`w(2) = 2` |
| `exists_arithDiv_ne_zero_at` | ★どの素点にも `arithDiv y` が非零な `y` がある |
| `toGp_sub_eq_divBGalois` | ★`Div_B` と有効因子の差の橋 |
| `ex63_hsp` | ★★`Definition 4.5, (ii)` の入力 |
| `ex63_not_isOfGroupLikeType` | ★`Φ ≠ 0` なので group-like でない |
| `botSub_end_eq_id` / `ex63_phi_map_bot_eq_id` | ★★`⊥` の自己射は恒等だけ |
| `divBGalois_nsmul_ne_zero` | ★`Div_B(y)` は無限位数 |
| `ex63_ratStd` | ★★★★★★**算術 Frobenioid は rationally standard 型** |
-/

namespace ABC3.Found.Divisor

open CategoryTheory NumberField ABC3.Found.FrdI ABC3.Meta

/-! ## ★1. どの素点にも「そこで非自明な有理関数」がある -/

/-- ★**無限素点での `|2| = 2`** —— `w = place φ` に戻して `‖φ 2‖ = ‖(2:ℂ)‖` を見るだけ。 -/
theorem infinitePlace_two (L : Type*) [Field L] [NumberField L] (w : InfinitePlace L) :
    w (2 : L) = 2 := by
  obtain ⟨φ, hφ⟩ := w.2
  show w.1 (2 : L) = 2
  rw [← hφ]
  show ‖φ (2 : L)‖ = 2
  rw [map_ofNat]
  norm_num

/-- ★★**どの素点 `v` にも `arithDiv y v ≠ 0` となる `y ∈ L^×` がある**。

★原文 `Theorem 6.4, (i)` の「immediate from the definition of B」の中身の半分。
★有限素点は一様化元、無限素点は `y = 2` でよい。 -/
theorem exists_arithDiv_ne_zero_at (L : Type*) [Field L] [NumberField L] (v : ArithPlace L) :
    ∃ y : Lˣ, arithDiv y v ≠ 0 := by
  rcases v with w | w
  · obtain ⟨π, hπ0, hπ1⟩ := exists_ordFin_eq_one (L := L) (FinitePlace.maximalIdeal w)
    refine ⟨Units.mk0 π hπ0, ?_⟩
    rw [arithDiv_apply]
    have hw : w = FinitePlace.mk (FinitePlace.maximalIdeal w) :=
      (FinitePlace.mk_maximalIdeal w).symm
    rw [show ((Units.mk0 π hπ0 : Lˣ) : L) = π from rfl, hw,
      arithPlaceLog_finite _ π hπ0, hπ1]
    have hpos := log_absNorm_pos (L := L) (FinitePlace.maximalIdeal w)
    push_cast
    intro hc
    rw [one_mul] at hc
    exact absurd hc (ne_of_gt hpos)
  · refine ⟨Units.mk0 (2 : L) two_ne_zero, ?_⟩
    rw [arithDiv_apply]
    show -(w.mult : ℝ) * Real.log (w ((Units.mk0 (2 : L) two_ne_zero : Lˣ) : L)) ≠ 0
    rw [show ((Units.mk0 (2 : L) two_ne_zero : Lˣ) : L) = (2 : L) from rfl,
      infinitePlace_two L w]
    have h1 : (0 : ℝ) < Real.log 2 := Real.log_pos (by norm_num)
    have h2 : 0 < w.mult := w.mult_pos
    have h3 : (0 : ℝ) < (w.mult : ℝ) := by exact_mod_cast h2
    intro hc
    nlinarith [hc]

/-! ## ★2. `Div_B` と有効因子の差の橋 -/

section Ex63

variable {F Kbar : Type} [Field F] [NumberField F] [Field Kbar] [Algebra F Kbar]
  [IsGalois F Kbar]

/-- ★★**有効因子の差が `arithDiv y` なら、その差は `Div_B(y)` である**。

★`Φ(L)^gp ≃ Γ`(`effRGpEquiv`)の単射性 1 本。 -/
theorem toGp_sub_eq_divBGalois (A : (FinSub F Kbar)ᵒᵖ)
    (y : (bmonGalois F Kbar).val A)
    (a b : effR (arithDivGroup A.unop.toIF))
    (hab : ((a : ArithPlace A.unop.toIF →₀ ℝ)) - ((b : ArithPlace A.unop.toIF →₀ ℝ))
      = arithDiv (Additive.toMul y)) :
    toGp _ a - toGp _ b = divBGalois A y := by
  refine effRGpHom_injective (arithDivGroup A.unop.toIF) ?_
  have h1 : effRGpHom (arithDivGroup A.unop.toIF) (toGp _ a - toGp _ b)
      = ⟨(a : ArithPlace A.unop.toIF →₀ ℝ) - (b : ArithPlace A.unop.toIF →₀ ℝ),
          by rw [hab]; exact arithDiv_mem_arithDivGroup _⟩ := by
    rw [effRGpHom, map_sub, gpLift_toGp, gpLift_toGp]
    rfl
  rw [h1]
  show _ = (effRGpEquiv (arithDivGroup A.unop.toIF) isGenSubgroupR_arithDivGroup)
      ((effRGpEquiv (arithDivGroup A.unop.toIF) isGenSubgroupR_arithDivGroup).symm
        (arithDivGroupHom A.unop.toIF y))
  rw [AddEquiv.apply_symm_apply]
  exact Subtype.ext hab

/-- ★★**`Div_B(y)` は `arithDiv y ≠ 0` なら無限位数**(`Φ^gp` は実係数なので捩れなし)。 -/
theorem divBGalois_nsmul_ne_zero (A : (FinSub F Kbar)ᵒᵖ)
    (y : ((A.unop.toIF))ˣ) (v : ArithPlace A.unop.toIF) (hy : arithDiv y v ≠ 0)
    (n : ℕ) (hn : 0 < n) :
    n • (divBGalois A (Additive.ofMul y)) ≠ 0 := by
  intro hc
  have h1 := congrArg (phiGpHom (arithDatumGalois F Kbar) finSubOp_isOfFSMType) hc
  rw [map_nsmul, map_zero] at h1
  have h2 : Subtype.val (phiGpHom (arithDatumGalois F Kbar) finSubOp_isOfFSMType
      (divBGalois A (Additive.ofMul y))) = arithDiv y :=
    phiGpHom_divBGalois A (Additive.ofMul y)
  have h4 := congrArg
    (fun w : ((arithDatumGalois F Kbar).grp A) => (Subtype.val w) v) h1
  simp only [AddSubgroup.coe_nsmul, Finsupp.smul_apply, ZeroMemClass.coe_zero,
    Finsupp.coe_zero, nsmul_eq_mul, h2] at h4
  rcases mul_eq_zero.mp h4 with hz | hz
  · exact absurd (Nat.cast_eq_zero.mp hz) hn.ne'
  · exact hy hz

/-! ## ★3. `⊥`(= `F` 自身)の自己射は恒等だけ

★★`botSub_end_eq_id` は `Sec6GaloisCat.lean` へ移した(`Example 6.1` と共用するため)。 -/

end Ex63

section Ex63Std

variable (F Kbar : Type) [Field F] [NumberField F] [Field Kbar] [Algebra F Kbar]
  [IsGalois F Kbar]

/-- ★★★**`Φ` は `⊥` の上では恒等しか誘導しない** —— これが `σ = id` の中身。 -/
theorem ex63_phi_map_bot_eq_id
    (e : (Opposite.op (botSub F Kbar)) ⟶ (Opposite.op (botSub F Kbar))) :
    (ex63ModelData (F := F) (Kbar := Kbar)).phi.map e
      = AddMonoidHom.id ((ex63ModelData (F := F) (Kbar := Kbar)).phi.val
          (Opposite.op (botSub F Kbar))) := by
  have he : e = 𝟙 _ := Quiver.Hom.unop_inj (botSub_end_eq_id e.unop)
  rw [he]
  ext x
  exact MonoidOn.map_id _ _ x

/-- ★★**因子分解写像 `ι`** —— `Definition 2.4, (i)` の (c)(d)。実係数の `iotaR` そのもの。 -/
noncomputable def ex63Iota :
    ∀ Y : (FinSub F Kbar)ᵒᵖ, Prime ((ex63ModelData (F := F) (Kbar := Kbar)).phi.val Y)
      → Pf ((ex63ModelData (F := F) (Kbar := Kbar)).phi.val Y) → NNReal :=
  fun Y => iotaR (Γ := arithDivGroup Y.unop.toIF) isGenSubgroupR_arithDivGroup

/-- ★★★★**[FrdI] Definition 4.5, (ii) の入力**(`strictly rational` の中身)。

原文 (FrdI p.115):
> Also, it is immediate from the definition of B that C is of [strictly] rational type, -/
theorem ex63_hsp :
    ∀ (A : ModelData.Obj (ex63ModelData (F := F) (Kbar := Kbar)))
      (p : Prime ((ex63ModelData (F := F) (Kbar := Kbar)).phi.val
        ((ModelData.modelPre (ex63Hyp F Kbar)).toElem.obj A).base)),
      ∃ (a b : (ex63ModelData (F := F) (Kbar := Kbar)).phi.val
          ((ModelData.modelPre (ex63Hyp F Kbar)).toElem.obj A).base)
        (y : (ModelData.modelRatFnData (ex63Hyp F Kbar)).bmon.val
          ((ModelData.modelPre (ex63Hyp F Kbar)).toElem.obj A).base),
        (toGp _ a - toGp _ b = (ModelData.modelRatFnData (ex63Hyp F Kbar)).divB _ y ∨
          toGp _ a - toGp _ b = -((ModelData.modelRatFnData (ex63Hyp F Kbar)).divB _ y)) ∧
        p ∈ SuppElt (ex63Iota F Kbar _) a ∧ p ∉ SuppElt (ex63Iota F Kbar _) b := by
  classical
  intro A p
  obtain ⟨y, hy⟩ := exists_arithDiv_ne_zero_at A.base.unop.toIF
    (effRPrimeEquiv (Γ := arithDivGroup A.base.unop.toIF) isGenSubgroupR_arithDivGroup p)
  obtain ⟨a, b, z, hzmem, hzeq, hab, hap, hbp⟩ :=
    exists_split_suppElt_R (Γ := arithDivGroup A.base.unop.toIF)
      isCoordwiseR_arithDivGroup isGenSubgroupR_arithDivGroup
      (arithDiv_mem_arithDivGroup y) p hy
  refine ⟨a, b, Additive.ofMul y, ?_, hap, hbp⟩
  have hsub : ((a : effR (arithDivGroup A.base.unop.toIF)) : ArithPlace _ →₀ ℝ)
      - ((b : effR (arithDivGroup A.base.unop.toIF)) : ArithPlace _ →₀ ℝ) = z := by
    rw [hab]; abel
  rcases hzeq with rfl | rfl
  · exact Or.inl (toGp_sub_eq_divBGalois A.base (Additive.ofMul y) a b hsub)
  · have h2 := toGp_sub_eq_divBGalois A.base (Additive.ofMul y) b a
      (by rw [← neg_sub, hsub]; abel)
    have h3 : toGp (effR (arithDivGroup A.base.unop.toIF)) a
        - toGp (effR (arithDivGroup A.base.unop.toIF)) b
        = -(divBGalois A.base (Additive.ofMul y)) := by
      rw [← h2]; abel
    exact Or.inr h3

/-- ★★★**[FrdI] Theorem 6.4, (i)** —— 算術 Frobenioid は **group-like 型でない**。

原文 (FrdI p.114):
> and rationally standard type, but not of group-like type; D is Frobenius-slim

★中身は「`Φ(L) ≠ 0` かつ sharp」だけ(`Φ(L)` は数体なら無限素点を必ず持つ)。 -/
theorem ex63_not_isOfGroupLikeType :
    ¬ IsOfGroupLikeType (ModelData.modelPre (ex63Hyp F Kbar)) := by
  intro hgl
  have hcc := hgl ⟨Opposite.op (botSub F Kbar), 0⟩
  obtain ⟨w⟩ := (inferInstance : Nonempty (InfinitePlace ((botSub F Kbar).toIF)))
  obtain ⟨xx, hxx⟩ := exists_ne_zero_arithEff (L := ((botSub F Kbar).toIF)) w
  exact hxx (eq_zero_of_isGroupLike_of_isSharp hcc (isSharp_effR _) xx)

/-- ★★★★**[FrdI] Theorem 6.4, (i)** —— `𝒞`・`𝒞^pf`・`𝒞^un-tr`・`(𝒞^pf)^un-tr` が
**group-like 型でない**。

原文 (FrdI p.114):
> and rationally standard type, but not of group-like type; D is Frobenius-slim

★伝播は `Proposition 5.5, (iii)`(在庫)。
★★`𝒞^rlf` の側は `scModel_not_isOfGroupLikeType` に係数拡大の入力
(`hcharInj` / `hint` / `hdiv` / `hperf` …)が要るので、`p55iii-ratstd` 節点に回した。 -/
theorem ex63_not_groupLike_family :
    ¬ IsOfGroupLikeType (ModelData.modelPre (ex63Hyp F Kbar))
      ∧ ¬ IsOfGroupLikeType (pfRootPre (ModelData.modelPre (ex63Hyp F Kbar))
          (ModelData.model_frobenioid (ex63Hyp F Kbar)).core)
      ∧ ¬ IsOfGroupLikeType (unTrPre (ModelData.modelPre (ex63Hyp F Kbar))
          (ModelData.model_frobenioid (ex63Hyp F Kbar)).core)
      ∧ ¬ IsOfGroupLikeType (unTrPre (pfRootPre (ModelData.modelPre (ex63Hyp F Kbar))
            (ModelData.model_frobenioid (ex63Hyp F Kbar)).core)
          (pfRootCore (F := (ModelData.model_frobenioid (ex63Hyp F Kbar)).core)
            (isOfFrobeniusIsotropicType_of_isotropic
              (ModelData.model_isotropicType (ex63Hyp F Kbar))))) := by
  haveI := (ex63Hyp F Kbar).connectedD
  have h0 := ex63_not_isOfGroupLikeType F Kbar
  exact ⟨h0, pfRoot_not_isOfGroupLikeType h0, unTr_not_isOfGroupLikeType _ h0,
    unTr_not_isOfGroupLikeType _ (pfRoot_not_isOfGroupLikeType h0)⟩

/-- ★★★★★★★**[FrdI] Theorem 6.4, (i)** ——
**算術 Frobenioid `C_{F̄/F}` は rationally standard 型**。

原文 (FrdI p.115):
> C is of [strictly] rational type, and that every object of (Cun-tr)birat is Frobenius-

## ★★★手順

1. `hsp`(strictly rational)—— `ex63_hsp`。
2. `hgl`(group-like なら compact 対象)—— 前件が偽(`ex63_not_isOfGroupLikeType`)。
3. `(𝒞^un-tr)^birat` の Frobenius-compact 対象 —— **底 `⊥` の上の対象**を選ぶ。
   * `σ = id` は `ex63_phi_map_bot_eq_id`。
   * `Φ^birat` の無限位数の元は `Div_B(2)`。
     ★これを **2 本の pre-step の差**として `Φ^birat` に入れるのが
     `toGp_div_sub_mem_phiBiratAt_of_preStep` である
     (`X = (⊥, 0)`、`A' = (⊥, Div_B(2) の負部分)`、
      `f` の除数が正部分・`g` の除数が負部分)。 -/
theorem ex63_ratStd :
    IsOfRationallyStandardType (ModelData.modelPre (ex63Hyp F Kbar))
      (ModelData.model_frobenioid (ex63Hyp F Kbar)) (ex63Iota F Kbar) := by
  classical
  set L := (botSub F Kbar).toIF with hL
  set Y : (FinSub F Kbar)ᵒᵖ := Opposite.op (botSub F Kbar) with hY
  set M := ex63ModelData (F := F) (Kbar := Kbar) with hM
  set h := ex63Hyp F Kbar with hh
  set P := ModelData.modelPre h with hP
  set G' := ModelData.model_frobenioid h with hG'
  set Fc := G'.core with hFc
  set v0 : ArithPlace L := Sum.inr (Classical.arbitrary (InfinitePlace L)) with hv0
  obtain ⟨y, hy0⟩ := exists_arithDiv_ne_zero_at L v0
  obtain ⟨a0, b0, ha0, hb0, hab0⟩ := exists_effR_sub (arithDivGroup L)
    isGenSubgroupR_arithDivGroup (arithDiv y) (arithDiv_mem_arithDivGroup y)
  set a : M.phi.val Y := ⟨a0, ha0⟩ with hadef
  set b : M.phi.val Y := ⟨b0, hb0⟩ with hbdef
  have hkey : toGp (M.phi.val Y) a - toGp (M.phi.val Y) b = M.divB Y (Additive.ofMul y) :=
    toGp_sub_eq_divBGalois Y (Additive.ofMul y) a b (by
      show a0 - b0 = arithDiv y
      rw [hab0])
  set X : ModelData.Obj M := ⟨Y, 0⟩ with hXdef
  set A' : ModelData.Obj M := ⟨Y, toGp (M.phi.val Y) b⟩ with hA'def
  set gm : X ⟶ A' :=
    { base := 𝟙 Y, div := b, deg := 1, u := 0,
      cond := by
        show ((1 : ℕ+) : ℕ) • (0 : Gp (M.phi.val Y)) + toGp (M.phi.val Y) b
          = M.phi.gpMapOn (𝟙 Y) (toGp (M.phi.val Y) b) + M.divB Y 0
        rw [M.phi.gpMapOn_id, map_zero, smul_zero, zero_add, add_zero] } with hgmdef
  set fm : X ⟶ A' :=
    { base := 𝟙 Y, div := a, deg := 1, u := Additive.ofMul y,
      cond := by
        show ((1 : ℕ+) : ℕ) • (0 : Gp (M.phi.val Y)) + toGp (M.phi.val Y) a
          = M.phi.gpMapOn (𝟙 Y) (toGp (M.phi.val Y) b) + M.divB Y (Additive.ofMul y)
        rw [M.phi.gpMapOn_id, smul_zero, zero_add, ← hkey]
        abel } with hfmdef
  set Xi : Istr P := ⟨X, ModelData.model_isotropicType h X⟩ with hXi
  set A'i : Istr P := ⟨A', ModelData.model_isotropicType h A'⟩ with hA'i
  set Q := unTrPre P Fc with hQ
  set GQ := unTr_frobenioid P Fc G' with hGQ
  set fu : (show UnTr P from Xi) ⟶ (show UnTr P from A'i) :=
    (istrToUnTr P).map (ObjectProperty.homMk fm) with hfu
  set gu : (show UnTr P from Xi) ⟶ (show UnTr P from A'i) :=
    (istrToUnTr P).map (ObjectProperty.homMk gm) with hgu
  haveI : IsIso (Q.Base fu) := by
    show IsIso (𝟙 Y)
    infer_instance
  have hmem := toGp_div_sub_mem_phiBiratAt_of_preStep Q GQ
    (fun Z => (unTr_isOfModelType Fc G').2 Z)
    (X := show UnTr P from Xi) (A := show UnTr P from A'i)
    (unTr_isotropic P Fc _) fu gu rfl rfl rfl
  refine ModelData.model_isOfRationallyStandardType_of_baseId h (ex63Iota F Kbar)
    (ModelData.modelRatFnData h) (ex63_hsp F Kbar)
    (isOfFSMFFType_of_isOfFSMType finSubOp_isOfFSMType)
    arithDatumGalois_isNonDilatingOn
    (fun hgl => absurd hgl (ex63_not_isOfGroupLikeType F Kbar))
    (show BiratCat Q GQ from (show UnTr P from Xi))
    (fun θ => ex63_phi_map_bot_eq_id F Kbar (biratBase θ.inv))
    (M.divB Y (Additive.ofMul y)) ?_ ?_
  · rw [← hkey]
    exact hmem
  · intro n hn
    exact divBGalois_nsmul_ne_zero Y y v0 hy0 n hn

/-- ★★★★★★★**[FrdI] Theorem 6.4, (i)(型の 3 条)** ——
算術 Frobenioid `C_{F̄/F}` と派生圏について

1. `C, C^pf, C^istr, C^un-tr, (C^pf)^un-tr` は **isotropic 型**
2. `C` は **rationally standard 型**
3. `C, C^pf, C^un-tr, (C^pf)^un-tr` は **group-like 型でない**

原文 (FrdI p.114):
> and rationally standard type, but not of group-like type; D is Frobenius-slim

★★2 の `C^pf` / `C^rlf` / `C^un-tr` / `(C^pf)^un-tr` への伝播は原文が
`Proposition 5.5, (iii)` に委ねている段であり、依存グラフでは
**`p55iii-ratstd`** 節点(`prop55` 鎖)に分けてある。 -/
theorem ex63_thm64_i_types :
    (IsOfIsotropicType (ModelData.modelPre (ex63Hyp F Kbar))
      ∧ IsOfFrobeniusIsotropicType (ModelData.modelPre (ex63Hyp F Kbar))
      ∧ IsOfIsotropicType (pfRootPre (ModelData.modelPre (ex63Hyp F Kbar))
          (ModelData.model_frobenioid (ex63Hyp F Kbar)).core)
      ∧ IsOfIsotropicType (istrPre (ModelData.modelPre (ex63Hyp F Kbar))
          (ModelData.model_frobenioid (ex63Hyp F Kbar)).core)
      ∧ IsOfFrobeniusIsotropicType (unTrPre (ModelData.modelPre (ex63Hyp F Kbar))
          (ModelData.model_frobenioid (ex63Hyp F Kbar)).core)
      ∧ IsOfIsotropicType (unTrPre (ModelData.modelPre (ex63Hyp F Kbar))
          (ModelData.model_frobenioid (ex63Hyp F Kbar)).core)
      ∧ IsOfIsotropicType
          (unTrPre (pfRootPre (ModelData.modelPre (ex63Hyp F Kbar))
              (ModelData.model_frobenioid (ex63Hyp F Kbar)).core)
            (pfRootCore (F := (ModelData.model_frobenioid (ex63Hyp F Kbar)).core)
              (isOfFrobeniusIsotropicType_of_isotropic
                (ModelData.model_isotropicType (ex63Hyp F Kbar))))))
    ∧ IsOfRationallyStandardType (ModelData.modelPre (ex63Hyp F Kbar))
        (ModelData.model_frobenioid (ex63Hyp F Kbar)) (ex63Iota F Kbar)
    ∧ (¬ IsOfGroupLikeType (ModelData.modelPre (ex63Hyp F Kbar))
      ∧ ¬ IsOfGroupLikeType (pfRootPre (ModelData.modelPre (ex63Hyp F Kbar))
          (ModelData.model_frobenioid (ex63Hyp F Kbar)).core)
      ∧ ¬ IsOfGroupLikeType (unTrPre (ModelData.modelPre (ex63Hyp F Kbar))
          (ModelData.model_frobenioid (ex63Hyp F Kbar)).core)
      ∧ ¬ IsOfGroupLikeType (unTrPre (pfRootPre (ModelData.modelPre (ex63Hyp F Kbar))
            (ModelData.model_frobenioid (ex63Hyp F Kbar)).core)
          (pfRootCore (F := (ModelData.model_frobenioid (ex63Hyp F Kbar)).core)
            (isOfFrobeniusIsotropicType_of_isotropic
              (ModelData.model_isotropicType (ex63Hyp F Kbar)))))) :=
  ⟨ex63_isotropic_family F Kbar, ex63_ratStd F Kbar, ex63_not_groupLike_family F Kbar⟩

/-- ★★★★**[FrdI] Theorem 6.4, (i)** —— 算術 Frobenioid は **standard 型**。 -/
theorem ex63_standardType :
    IsOfStandardType ((FinSub F Kbar)ᵒᵖ) _
      (ModelData.modelPre (ex63Hyp F Kbar))
      (ModelData.model_frobenioid (ex63Hyp F Kbar)).core :=
  ModelData.model_isOfStandardType (ex63Hyp F Kbar)
    (ModelData.model_frobenioid (ex63Hyp F Kbar)).core
    (isOfFSMFFType_of_isOfFSMType finSubOp_isOfFSMType)
    arithDatumGalois_isNonDilatingOn
    (fun hgl => absurd hgl (ex63_not_isOfGroupLikeType F Kbar))

/-- ★★★★★★**[FrdI] Proposition 5.5, (iii)** —— 算術の **`𝒞^un-tr` も rationally standard**。

★一般形は `Thm64RatStd.lean` の `unTr_isOfRationallyStandardType`。
★★`Theorem 6.4, (i)` が並べる 5 圏のうち **`𝒞^un-tr` の 1 つ**がこれで閉じる。 -/
theorem ex63_unTr_ratStd :
    IsOfRationallyStandardType
      (unTrPre (ModelData.modelPre (ex63Hyp F Kbar))
        (ModelData.model_frobenioid (ex63Hyp F Kbar)).core)
      (unTr_frobenioid (ModelData.modelPre (ex63Hyp F Kbar))
        (ModelData.model_frobenioid (ex63Hyp F Kbar)).core
        (ModelData.model_frobenioid (ex63Hyp F Kbar)))
      (ex63Iota F Kbar) := by
  classical
  haveI := (ex63Hyp F Kbar).connectedD
  set L := (botSub F Kbar).toIF with hL
  set Y : (FinSub F Kbar)ᵒᵖ := Opposite.op (botSub F Kbar) with hY
  set M := ex63ModelData (F := F) (Kbar := Kbar) with hM
  set h := ex63Hyp F Kbar with hh
  set X0 : ModelData.Obj M := ⟨Y, 0⟩ with hX0
  set Xi : Istr (ModelData.modelPre h) := ⟨X0, ModelData.model_isotropicType h X0⟩ with hXi
  set Xi2 : Istr (unTrPre (ModelData.modelPre h)
      (ModelData.model_frobenioid h).core) :=
    ⟨show UnTr (ModelData.modelPre h) from Xi,
      unTr_isotropic (ModelData.modelPre h) (ModelData.model_frobenioid h).core _⟩ with hXi2
  set v0 : ArithPlace L := Sum.inr (Classical.arbitrary (InfinitePlace L)) with hv0
  obtain ⟨y, hy⟩ := exists_arithDiv_ne_zero_at L v0
  exact unTr_isOfRationallyStandardType h (ex63Iota F Kbar) (ex63_hsp F Kbar)
    (fun A => (isDivisorial_arithEff (L := A.unop.toIF)).1.1)
    finSubOp_isOfFSMType
    (ex63_not_isOfGroupLikeType F Kbar)
    (ex63_standardType F Kbar)
    (show UnTr (unTrPre (ModelData.modelPre h)
      (ModelData.model_frobenioid h).core) from Xi2)
    (fun e => ex63_phi_map_bot_eq_id F Kbar e)
    (Additive.ofMul y)
    (fun n hn => divBGalois_nsmul_ne_zero Y y v0 hy n hn)

set_option maxHeartbeats 1000000 in
/-- ★★★★★★**[FrdI] Proposition 5.5, (iii)** —— 算術の **`𝒞^pf` も rationally standard**。

★一般形は `Prop55PfRatStd.lean` の `pfRoot_isOfRationallyStandardType`。
★★`Theorem 6.4, (i)` が並べる 5 圏のうち **`𝒞^pf` の 1 つ**がこれで閉じる
(残るは `𝒞^rlf` と `(𝒞^pf)^un-tr`)。

★(b) の Frobenius-compact 対象は `𝒞^un-tr` のときと同じく **底 `⊥`(= `F` 自身)の対象**で取る。
`Div_B(2)` が `(Φ^pf)^gp` でも無限位数であることは `nsmul_gpMap_pfOf_ne_zero` による。 -/
theorem ex63_pf_ratStd :
    IsOfRationallyStandardType
      (pfRootPre (ModelData.modelPre (ex63Hyp F Kbar))
        (ModelData.model_frobenioid (ex63Hyp F Kbar)).core)
      (pfRoot_frobenioid (F := (ModelData.model_frobenioid (ex63Hyp F Kbar)).core)
        (isOfFrobeniusIsotropicType_of_isotropic
          (ModelData.model_isotropicType (ex63Hyp F Kbar)))
        (ModelData.model_frobenioid (ex63Hyp F Kbar)))
      (fun Y => iotaPf ((ModelData.modelPre (ex63Hyp F Kbar)).divisorial Y)
        (ex63Iota F Kbar Y)) := by
  classical
  haveI := (ex63Hyp F Kbar).connectedD
  set L := (botSub F Kbar).toIF with hL
  set Y : (FinSub F Kbar)ᵒᵖ := Opposite.op (botSub F Kbar) with hY
  set h := ex63Hyp F Kbar with hh
  set v0 : ArithPlace L := Sum.inr (Classical.arbitrary (InfinitePlace L)) with hv0
  obtain ⟨y, hy⟩ := exists_arithDiv_ne_zero_at L v0
  exact pfRoot_isOfRationallyStandardType
    (isOfFrobeniusIsotropicType_of_isotropic (ModelData.model_isotropicType h))
    (ModelData.model_isotropicType h)
    (ex63_isotropic_family F Kbar).2.2.1
    _
    (ModelData.model_birat_frobenioidCore h)
    (ex63Iota F Kbar)
    (ModelData.model_isOfBiratFrobNormalizedType h)
    (fun A => isStrictlyRational_of_divB (ModelData.modelRatFnData h) (ex63Iota F Kbar) A
      (ex63_hsp F Kbar A))
    (ex63_not_isOfGroupLikeType F Kbar)
    (ex63_standardType F Kbar)
    (unTr_pf_biratCompact_of_baseId h ⟨Y, 0⟩ (Additive.ofMul y)
      (fun e => ex63_phi_map_bot_eq_id F Kbar e)
      (fun n hn => nsmul_gpMap_pfOf_ne_zero ((ModelData.modelPre h).divisorial Y) n
        (divBGalois_nsmul_ne_zero Y y v0 hy n hn)))

set_option maxHeartbeats 1000000 in
/-- ★★★★★★**[FrdI] Proposition 5.5, (iii)** ——
算術の **`(𝒞^pf)^un-tr` も rationally standard**。

★一般形は `Prop55UnTrRatStd.lean` の `pfRoot_unTr_isOfRationallyStandardType`。
★★`Theorem 6.4, (i)` が並べる 5 圏のうち **4 つ目**がこれで閉じる
(`𝒞`・`𝒞^un-tr`・`𝒞^pf`・`(𝒞^pf)^un-tr`)。残るは `𝒞^rlf` である。 -/
theorem ex63_pf_unTr_ratStd :
    IsOfRationallyStandardType
      (unTrPre (pfRootPre (ModelData.modelPre (ex63Hyp F Kbar))
          (ModelData.model_frobenioid (ex63Hyp F Kbar)).core)
        (pfRoot_frobenioid (F := (ModelData.model_frobenioid (ex63Hyp F Kbar)).core)
          (isOfFrobeniusIsotropicType_of_isotropic
            (ModelData.model_isotropicType (ex63Hyp F Kbar)))
          (ModelData.model_frobenioid (ex63Hyp F Kbar))).core)
      (unTr_frobenioid (pfRootPre (ModelData.modelPre (ex63Hyp F Kbar))
          (ModelData.model_frobenioid (ex63Hyp F Kbar)).core)
        (pfRoot_frobenioid (F := (ModelData.model_frobenioid (ex63Hyp F Kbar)).core)
          (isOfFrobeniusIsotropicType_of_isotropic
            (ModelData.model_isotropicType (ex63Hyp F Kbar)))
          (ModelData.model_frobenioid (ex63Hyp F Kbar))).core
        (pfRoot_frobenioid (F := (ModelData.model_frobenioid (ex63Hyp F Kbar)).core)
          (isOfFrobeniusIsotropicType_of_isotropic
            (ModelData.model_isotropicType (ex63Hyp F Kbar)))
          (ModelData.model_frobenioid (ex63Hyp F Kbar))))
      (fun Y => iotaPf ((ModelData.modelPre (ex63Hyp F Kbar)).divisorial Y)
        (ex63Iota F Kbar Y)) := by
  classical
  haveI := (ex63Hyp F Kbar).connectedD
  set L := (botSub F Kbar).toIF with hL
  set Y : (FinSub F Kbar)ᵒᵖ := Opposite.op (botSub F Kbar) with hY
  set v0 : ArithPlace L := Sum.inr (Classical.arbitrary (InfinitePlace L)) with hv0
  obtain ⟨y, hy⟩ := exists_arithDiv_ne_zero_at L v0
  exact pfRoot_unTr_isOfRationallyStandardType
    (isOfFrobeniusIsotropicType_of_isotropic
      (ModelData.model_isotropicType (ex63Hyp F Kbar)))
    (ex63_isotropic_family F Kbar).2.2.1
    (pfRoot_frobenioid (F := (ModelData.model_frobenioid (ex63Hyp F Kbar)).core)
      (isOfFrobeniusIsotropicType_of_isotropic
        (ModelData.model_isotropicType (ex63Hyp F Kbar)))
      (ModelData.model_frobenioid (ex63Hyp F Kbar)))
    (ModelData.model_birat_frobenioidCore (ex63Hyp F Kbar))
    (ex63Iota F Kbar)
    (fun A => isStrictlyRational_of_divB (ModelData.modelRatFnData (ex63Hyp F Kbar))
      (ex63Iota F Kbar) A (ex63_hsp F Kbar A))
    finSubOp_isOfFSMType
    (ex63_not_isOfGroupLikeType F Kbar)
    (ex63_standardType F Kbar)
    ⟨Y, 0⟩
    (fun e => ex63_phi_map_bot_eq_id F Kbar e)
    ((ModelData.modelRatFnData (ex63Hyp F Kbar)).divB _ (Additive.ofMul y))
    (divB_mem_phiBiratAt (ModelData.modelRatFnData (ex63Hyp F Kbar))
      (A := (⟨Y, 0⟩ : ModelData.Obj _)) (Additive.ofMul y))
    (fun n hn => nsmul_gpMap_pfOf_ne_zero
      ((ModelData.modelPre (ex63Hyp F Kbar)).divisorial Y) n
      (divBGalois_nsmul_ne_zero Y y v0 hy n hn))

end Ex63Std

/-! ### ★出典の紐付け(`.src`)と、証明が要求するもの(`.needs`) -/

def ex63_standardType.src : Source :=
  { paper := "FrdI", pdfPage := 114,
    item := "Theorem 6.4, (i) — 算術 Frobenioid は standard 型",
    sectionId := "frdi-thm-6-4" }

def ex63_unTr_ratStd.src : Source :=
  { paper := "FrdI", pdfPage := 105,
    item := "Proposition 5.5, (iii) — 算術の 𝒞^un-tr も rationally standard",
    sectionId := "frdi-prop-5-5" }

def ex63_pf_ratStd.src : Source :=
  { paper := "FrdI", pdfPage := 105,
    item := "Proposition 5.5, (iii) — 算術の 𝒞^pf も rationally standard",
    sectionId := "frdi-prop-5-5" }

def ex63_pf_unTr_ratStd.src : Source :=
  { paper := "FrdI", pdfPage := 105,
    item := "Proposition 5.5, (iii) — 算術の (𝒞^pf)^un-tr も rationally standard",
    sectionId := "frdi-prop-5-5" }

def ex63_pf_ratStd.needs : List ProofObligation :=
  [ .citation "[ABC3]" "pfRoot_isOfRationallyStandardType(一般形。4 条を集める)"
      (.inProject "ABC3" "ABC3.Found.FrdI.pfRoot_isOfRationallyStandardType") 105,
    .citation "[ABC3]" "unTr_pf_biratCompact_of_baseId((b) の Frobenius-compact 対象)"
      (.inProject "ABC3" "ABC3.Found.FrdI.unTr_pf_biratCompact_of_baseId") 114,
    .citation "[ABC3]" "ModelData.model_birat_frobenioidCore(𝒞^birat の 21 条)"
      (.inProject "ABC3" "ABC3.Found.FrdI.ModelData.model_birat_frobenioidCore") 85,
    .citation "[ABC3]" "nsmul_gpMap_pfOf_ne_zero(Div_B(2) は (Φ^pf)^gp でも無限位数)"
      (.inProject "ABC3" "ABC3.Found.FrdI.nsmul_gpMap_pfOf_ne_zero") 105 ]

def exists_arithDiv_ne_zero_at.src : Source :=
  { paper := "FrdI", pdfPage := 114,
    item := "Theorem 6.4, (i) — どの素点にもそこで非自明な有理関数がある",
    sectionId := "frdi-thm-6-4" }

def ex63_hsp.src : Source :=
  { paper := "FrdI", pdfPage := 114,
    item := "Theorem 6.4, (i) — 算術 Frobenioid は strictly rational 型",
    sectionId := "frdi-thm-6-4" }

def ex63_hsp.needs : List ProofObligation :=
  [ .citation "[ABC3]" "exists_ordFin_eq_one(有限素点の一様化元)"
      (.inProject "ABC3" "ABC3.Found.Divisor.exists_ordFin_eq_one") 114,
    .citation "[ABC3]" "exists_split_suppElt_R(正負への分割、実係数)"
      (.inProject "ABC3" "ABC3.Found.Divisor.exists_split_suppElt_R") 114,
    .implicitStep
      "★原文は「immediate from the definition of B」。無限素点で y = 2 を取る段が畳まれている" 114 ]

def ex63_not_isOfGroupLikeType.src : Source :=
  { paper := "FrdI", pdfPage := 114,
    item := "Theorem 6.4, (i) — 算術 Frobenioid は group-like 型でない",
    sectionId := "frdi-thm-6-4" }

def ex63_thm64_i_types.src : Source :=
  { paper := "FrdI", pdfPage := 114,
    item := "Theorem 6.4, (i) — 5 つの圏の型(isotropic / rationally standard / not group-like)",
    sectionId := "frdi-thm-6-4" }

def ex63_thm64_i_types.needs : List ProofObligation :=
  [ .citation "[ABC3]" "ex63_isotropic_family(isotropic の側)"
      (.inProject "ABC3" "ABC3.Found.Divisor.ex63_isotropic_family") 114,
    .citation "[ABC3]" "ex63_ratStd(rationally standard の側)"
      (.inProject "ABC3" "ABC3.Found.Divisor.ex63_ratStd") 114,
    .citation "[FrdI]" "Proposition 5.5, (iii) — rationally standard 性の 𝒞^rlf への伝播"
      (.absent ("2026-08-25 実測(第 2 版)。5 圏のうち 𝒞 / 𝒞^un-tr / 𝒞^pf / (𝒞^pf)^un-tr は" ++
        "ex63_ratStd / ex63_unTr_ratStd / ex63_pf_ratStd / ex63_pf_unTr_ratStd で閉じた。" ++
        "残るは 𝒞^rlf(実化)のみ。依存グラフの p55iii-ratstd 節点に分けてある")) 114,
    .citation "[ABC3]" "not group-like の伝播(𝒞^pf / 𝒞^un-tr / (𝒞^pf)^un-tr)"
      (.inProject "ABC3" "ABC3.Found.FrdI.pfRoot_not_isOfGroupLikeType") 114 ]

def ex63_ratStd.src : Source :=
  { paper := "FrdI", pdfPage := 114,
    item := "Theorem 6.4, (i) — 算術 Frobenioid は rationally standard 型",
    sectionId := "frdi-thm-6-4" }

def ex63_ratStd.needs : List ProofObligation :=
  [ .citation "[ABC3]" "model_isOfRationallyStandardType_of_baseId"
      (.inProject "ABC3"
        "ABC3.Found.FrdI.ModelData.model_isOfRationallyStandardType_of_baseId") 114,
    .citation "[ABC3]" "toGp_div_sub_mem_phiBiratAt_of_preStep(Φ^birat の元を作る)"
      (.inProject "ABC3" "ABC3.Found.FrdI.toGp_div_sub_mem_phiBiratAt_of_preStep") 114,
    .citation "[ABC3]" "arithDatumGalois_isNonDilatingOn"
      (.inProject "ABC3" "ABC3.Found.Divisor.arithDatumGalois_isNonDilatingOn") 114,
    .implicitStep
      ("★逸脱の記録(2026-08-25)—— 原文は「every object of (C^un-tr)^birat is Frobenius-compact」" ++
       "と**すべての対象**について言うが、Definition 4.5, (iii), (b) は**存在**しか要求しない。" ++
       "本形式化は底が ⊥(= F 自身)の対象 1 つを取り、そこで σ = id を無料で得る。" ++
       "★在庫の birat_isFrobeniusCompact_of_primary は「Φ の非零元が Φ^birat に入る」を要求するが、" ++
       "算術では積公式より Φ(L) ∩ Φ^birat(L) = {0} なので満たしようがない。") 114 ]

end ABC3.Found.Divisor
