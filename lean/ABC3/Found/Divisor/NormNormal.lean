/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Divisor.SchemeNormalCriterion
import Mathlib.RingTheory.DedekindDomain.IntegralClosure

/-!
# `V[L]` は正規である(鎖 `normalize` の `normalization-proper-normal`)

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.109。

原文 (FrdI p.109):
> a proper normal variety], then let us write DL for the set of prime divisors of V [L]

## ★★三段で出る

| 段 | mathlib の在庫 |
|---|---|
| アフィン局所の切断 = 整閉包 | `Scheme.Hom.normalizationObjIso` |
| 整閉包は整閉 | `integralClosure.isIntegrallyClosedOfFiniteExtension` |
| アフィン局所の整閉性 ⟹ 正規 | `isNormalScheme_of_exists_affine`(在庫) |

## ★★★繋ぎの 1 段が本ファイルの中身

`normalizationObjIso` が与えるのは `integralClosure Γ(V,U) Γ(Spec L, f ⁻¹ᵁ U)` であって、
`integralClosure Γ(V,U) L` ではない。★両者を繋ぐのは

  `Γ(Spec L, f ⁻¹ᵁ U) = Γ(Spec L, ⊤) ≅ L`   (`f ⁻¹ᵁ U = ⊤` は `specToV_preimage_eq_top`)

であり、★★**この同型が `Γ(V,U)`-代数同型であること**が要る。それは

  `f.appLE U ⊤ ≫ ΓSpecIso.hom = germ ≫ (K → L)`

という 1 本の等式で、`f = Spec.map φ ≫ V.fromSpecStalk ξ` を
`Scheme.fromSpecStalk_app` と `ΓSpecIso_naturality` で開くだけで出る。

## ★本ファイルで閉じること

| 定理 | 中身 |
|---|---|
| `integralClosureCongr` | 整閉包は代数同型で移る |
| `specToV_appLE_ΓSpecIso` | ★★★**繋ぎの 1 本**(`f.appLE ≫ ΓSpecIso = germ ≫ 拡大`) |
| `sectionsAlgEquiv` | ★`Γ(Spec L, f ⁻¹ᵁ U) ≃ₐ[Γ(V,U)] L` |
| `isIntegrallyClosed_normSections` | ★★`Γ(V[L], D ⁻¹ᵁ U)` は整閉 |
| `normObj_isNormalScheme` | ★★★★★★**`V[L]` は正規** |
-/

namespace ABC3.Found.Divisor

open AlgebraicGeometry CategoryTheory ABC3.Found.FrdI

universe u

/-! ## ★1. 整閉包は代数同型で移る -/

section IntClosure

variable {R S S' : Type*} [CommRing R] [CommRing S] [CommRing S']
  [Algebra R S] [Algebra R S']

/-- ★★**整閉包は `R`-代数同型で移る**。 -/
def integralClosureCongr (e : S ≃ₐ[R] S') :
    integralClosure R S ≃ₐ[R] integralClosure R S' where
  toFun x := ⟨e x, IsIntegral.map (e : S →ₐ[R] S') x.2⟩
  invFun y := ⟨e.symm y, IsIntegral.map (e.symm : S' →ₐ[R] S) y.2⟩
  left_inv x := Subtype.ext (e.symm_apply_apply _)
  right_inv y := Subtype.ext (e.apply_symm_apply _)
  map_mul' x y := Subtype.ext (map_mul e _ _)
  map_add' x y := Subtype.ext (map_add e _ _)
  commutes' r := Subtype.ext (e.commutes r)

end IntClosure

/-! ## ★2. 繋ぎの 1 本 -/

variable (V : Scheme.{u}) [IsIntegral V] {Kbar : Type u} [Field Kbar]
  [Algebra V.functionField Kbar]

/-- ★生成点の茎からの射による、空でない開集合の逆像も全体。 -/
theorem fromSpecStalk_preimage_eq_top {U : V.Opens} (hU : (U : Set V).Nonempty) :
    (V.fromSpecStalk (genericPoint (V : Type u))) ⁻¹ᵁ U = ⊤ := by
  refine TopologicalSpace.Opens.ext (Set.eq_univ_iff_forall.mpr fun y => ?_)
  show (V.fromSpecStalk (genericPoint (V : Type u))).base y ∈ U
  have hmem : (V.fromSpecStalk (genericPoint (V : Type u))).base y
      ∈ Set.range (V.fromSpecStalk (genericPoint (V : Type u))).base := ⟨y, rfl⟩
  rw [Scheme.range_fromSpecStalk] at hmem
  exact hmem.mem_open U.2 (genericPoint_mem_of_nonempty V hU)

set_option backward.isDefEq.respectTransparency false in
/-- ★★**茎からの射の切断への作用は `germ` そのもの**(`⊤` へ制限した形)。 -/
theorem fromSpecStalk_appLE_ΓSpecIso {X : Scheme.{u}} {x : X} {U : X.Opens} (hxU : x ∈ U)
    (h : (X.fromSpecStalk x) ⁻¹ᵁ U = ⊤) :
    (X.fromSpecStalk x).appLE U ⊤ (le_of_eq h.symm)
        ≫ (Scheme.ΓSpecIso (X.presheaf.stalk x)).hom
      = X.presheaf.germ U x hxU := by
  rw [Scheme.Hom.appLE, Scheme.fromSpecStalk_app hxU]
  simp only [Category.assoc]
  rw [← Functor.map_comp_assoc, Subsingleton.elim
    ((homOfLE (le_top : X.fromSpecStalk x ⁻¹ᵁ U ≤ ⊤)).op ≫
      (homOfLE (le_of_eq h.symm)).op)
    (𝟙 (Opposite.op (⊤ : (Spec (X.presheaf.stalk x)).Opens))),
    CategoryTheory.Functor.map_id, Category.id_comp, Iso.inv_hom_id, Category.comp_id]

/-- ★★**前に `Spec.map φ` を合成すると、切断への作用に `φ` が付く**。 -/
theorem comp_appLE_ΓSpecIso {Y : Scheme.{u}} {R S : CommRingCat.{u}} (φ : R ⟶ S)
    (g : Spec R ⟶ Y) {U : Y.Opens} (hg : g ⁻¹ᵁ U = ⊤)
    (h : (Spec.map φ ≫ g) ⁻¹ᵁ U = ⊤) :
    (Spec.map φ ≫ g).appLE U ⊤ (le_of_eq h.symm) ≫ (Scheme.ΓSpecIso S).hom
      = (g.appLE U ⊤ (le_of_eq hg.symm) ≫ (Scheme.ΓSpecIso R).hom) ≫ φ := by
  have hmap : (Spec.map φ).appTop = (Spec.map φ).appLE ⊤ ⊤ (le_of_eq (by simp)) := by
    rw [Scheme.Hom.appLE, Scheme.Hom.appTop]
    simp
  simp only [Category.assoc]
  rw [← Scheme.ΓSpecIso_naturality φ, ← Category.assoc, hmap,
    Scheme.Hom.appLE_comp_appLE]

set_option backward.isDefEq.respectTransparency false in
set_option maxHeartbeats 1000000 in
/-- ★★★★**繋ぎの 1 本** —— `Spec L → V` の切断への作用は、
生成点への `germ` に体の拡大を合成したものである。 -/
theorem specToV_appLE_ΓSpecIso (L : FinSub V.functionField Kbar) {U : V.Opens}
    (hU : (U : Set V).Nonempty) :
    (specToV V L).appLE U ⊤ (le_of_eq (specToV_preimage_eq_top V L hU).symm)
        ≫ (Scheme.ΓSpecIso (CommRingCat.of L.toIF)).hom
      = V.presheaf.germ U (genericPoint (V : Type u)) (genericPoint_mem_of_nonempty V hU)
          ≫ CommRingCat.ofHom (algebraMap V.functionField L.toIF) := by
  have hxU : genericPoint (V : Type u) ∈ U := genericPoint_mem_of_nonempty V hU
  have hg : (V.fromSpecStalk (genericPoint (V : Type u))) ⁻¹ᵁ U = ⊤ :=
    fromSpecStalk_preimage_eq_top V hU
  show (Spec.map (CommRingCat.ofHom (algebraMap V.functionField L.toIF))
      ≫ V.fromSpecStalk (genericPoint (V : Type u))).appLE U ⊤
        (le_of_eq (specToV_preimage_eq_top V L hU).symm)
      ≫ (Scheme.ΓSpecIso (CommRingCat.of L.toIF)).hom = _
  rw [comp_appLE_ΓSpecIso _ _ hg (specToV_preimage_eq_top V L hU),
    fromSpecStalk_appLE_ΓSpecIso hxU hg]

/-! ## ★3. `Γ(Spec L, f ⁻¹ᵁ U) ≅ L` -/

/-- ★`f ⁻¹ᵁ U = ⊤` から来る同型 `Γ(Spec L, f ⁻¹ᵁ U) ≅ L`。 -/
noncomputable def sectionsIso (L : FinSub V.functionField Kbar) {U : V.Opens}
    (hU : (U : Set V).Nonempty) :
    Γ(Spec (CommRingCat.of L.toIF), (specToV V L) ⁻¹ᵁ U) ≅ CommRingCat.of L.toIF :=
  (Spec (CommRingCat.of L.toIF)).presheaf.mapIso
      (eqToIso (congrArg Opposite.op (specToV_preimage_eq_top V L hU)))
    ≪≫ Scheme.ΓSpecIso _

/-- ★★**`sectionsIso` は `Γ(V,U)`-代数の射である** —— 繋ぎの 1 本の言い換え。 -/
theorem app_comp_sectionsIso (L : FinSub V.functionField Kbar) {U : V.Opens}
    (hU : (U : Set V).Nonempty) :
    (specToV V L).app U ≫ (sectionsIso V L hU).hom
      = V.presheaf.germ U (genericPoint (V : Type u)) (genericPoint_mem_of_nonempty V hU)
          ≫ CommRingCat.ofHom (algebraMap V.functionField L.toIF) := by
  rw [← specToV_appLE_ΓSpecIso V L hU, Scheme.Hom.appLE, sectionsIso]
  simp only [Iso.trans_hom, Functor.mapIso_hom, Category.assoc]
  congr 2

/-! ## ★4. アフィン開の切断は整閉 -/

set_option maxHeartbeats 1000000 in
/-- ★★★★**`Γ(V[L], D ⁻¹ᵁ U)` は整閉**(`U` は空でないアフィン開)。

★★中身は 3 本:
* `normalizationObjIso` —— 切断は整閉包である
* `sectionsIso` / `app_comp_sectionsIso` —— `Γ(Spec L, f ⁻¹ᵁ U) ≅ L` は `Γ(V,U)`-代数同型
* `integralClosure.isIntegrallyClosedOfFiniteExtension` —— 整閉包は整閉 -/
theorem isIntegrallyClosed_normSections (L : FinSub V.functionField Kbar) {U : V.Opens}
    (hUaff : IsAffineOpen U) (hU : (U : Set V).Nonempty) :
    IsIntegrallyClosed Γ(normObj V L, normDown V L ⁻¹ᵁ U) := by
  haveI : Nonempty U := ⟨⟨hU.some, hU.some_mem⟩⟩
  haveI := L.fin
  letI : Algebra Γ(V, U) L.toIF :=
    ((algebraMap V.functionField L.toIF).comp (V.germToFunctionField U).hom).toAlgebra
  haveI : IsScalarTower Γ(V, U) V.functionField L.toIF :=
    IsScalarTower.of_algebraMap_eq' (by
      show ((algebraMap V.functionField L.toIF).comp (V.germToFunctionField U).hom) = _
      rfl)
  haveI := functionField_isFractionRing_of_isAffineOpen (X := V) U hUaff
  haveI : IsIntegrallyClosed (integralClosure Γ(V, U) L.toIF) :=
    integralClosure.isIntegrallyClosedOfFiniteExtension V.functionField
  letI := ((specToV V L).app U).hom.toAlgebra
  let e : Γ(Spec (CommRingCat.of L.toIF), (specToV V L) ⁻¹ᵁ U) ≃ₐ[Γ(V, U)] L.toIF :=
    { (sectionsIso V L hU).commRingCatIsoToRingEquiv with
      commutes' := fun r => by
        show (sectionsIso V L hU).hom.hom (((specToV V L).app U).hom r) = _
        exact congrArg (fun t : Γ(V, U) ⟶ CommRingCat.of L.toIF => CommRingCat.Hom.hom t r)
          (app_comp_sectionsIso V L hU) }
  haveI : IsIntegrallyClosed (integralClosure Γ(V, U)
      Γ(Spec (CommRingCat.of L.toIF), (specToV V L) ⁻¹ᵁ U)) :=
    IsIntegrallyClosed.of_equiv (integralClosureCongr e).symm.toRingEquiv
  exact IsIntegrallyClosed.of_equiv
    (Scheme.Hom.normalizationObjIso (specToV V L) hUaff).symm.commRingCatIsoToRingEquiv

/-! ## ★5. `V[L]` は正規 -/

/-- ★★★★★★**[FrdI] Example 6.1 —— `V[L]` は正規である**。

原文 (FrdI p.109):
> a proper normal variety], then let us write DL for the set of prime divisors of V [L]

★★各点 `x` について、`x` の像を含む `V` のアフィン開 `U` を取れば
`D ⁻¹ᵁ U` は `x` を含むアフィン開(`D` は整射、したがってアフィン射)であり、
その切断は整閉である(`isIntegrallyClosed_normSections`)。 -/
theorem normObj_isNormalScheme (L : FinSub V.functionField Kbar) :
    IsNormalScheme (normObj V L) := by
  refine isNormalScheme_of_exists_affine _ (fun x => ?_)
  obtain ⟨_, ⟨U, hUaff, rfl⟩, hxU, -⟩ :=
    V.isBasis_affineOpens.exists_subset_of_mem_open
      (Set.mem_univ ((normDown V L).base x)) isOpen_univ
  refine ⟨(normDown V L) ⁻¹ᵁ U, hUaff.preimage (normDown V L), hxU, ?_⟩
  exact isIntegrallyClosed_normSections V L hUaff ⟨(normDown V L).base x, hxU⟩

/-! ### ★出典の紐付け -/

/-- ★★★★★★locator —— `Example 6.1` の「`V[L]` は正規多様体である」。 -/
def normObj_isNormalScheme.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 109,
    item := "Example 6.1 — V[L] は正規である",
    sectionId := "frdi-example-6-1" }

end ABC3.Found.Divisor
