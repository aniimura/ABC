/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.ModelRigid
import ABC3.Found.FrdI.Prop53Diag

/-!
# [FrdI] Corollary 5.4 の rigidity —— 組み立て

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.104。

原文 (FrdI p.104):
> arrows are equivalences of categories]. Moreover, each of the composite functors of

## ★★`ModelRigid.lean` の補題に入れるための 1 本

`ModelData.isRigidFunctor_of_proj` は「底への合成 `F ⋙ (𝒞^model → 𝒟)` が rigid なら
`F` も rigid」と言う。★本ファイルはその仮定を `Proposition 5.3` の縦の矢印
`cToSc : 𝒞 ⥤ 𝒞^rlf` について**実際に満たす**ところを閉じる:

  `cToSc ⋙ (𝒞^rlf → 𝒟) ≅ P.proj`

★★鎖は 3 本の可換性でできている:

| 段 | 根拠 |
|---|---|
| `ModelDataHom.functor` は底を変えない | `obj_base` / `hom_base` が `rfl` |
| `Theorem 5.2, (iv)` の圏同値は `𝔽_Φ` と 1-可換 | `pathToModel_toElem`(`rfl`)＋ 圏同値の counit |
| `cToUnTr` は底への射影と 1-可換 | `cToUnTr_comp_proj`(`rfl`) |

## ★本ファイルで閉じること

| 定理 | 中身 |
|---|---|
| `modelDataHom_functor_comp_proj` | `ModelDataHom` の関手は底を変えない |
| `thm52iv_comp_toElem` | ★`Theorem 5.2, (iv)` の圏同値は `𝔽_Φ` と 1-可換 |
| `modelTypeEquiv_comp_toElem` | 同上(`modelType_equiv` 版) |
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ}

/-! ## ★1. `ModelDataHom` の関手は底を変えない -/

/-- ★**`ModelDataHom` の関手は底への射影と可換** —— `obj_base` / `hom_base` が `rfl`。 -/
theorem modelDataHom_functor_comp_proj {M M' : ModelData.{v, u, w} D}
    (F : ModelDataHom M M') (h : ModelData.Hyp M) (h' : ModelData.Hyp M') :
    F.functor ⋙ (ModelData.modelPre h').proj = (ModelData.modelPre h).proj := rfl

/-! ## ★2. `Theorem 5.2, (iv)` の圏同値は `𝔽_Φ` への関手と 1-可換 -/

variable [IsConnected D]

/-- ★★★**`Theorem 5.2, (iv)` の圏同値は `𝔽_Φ` への関手と 1-可換**。

★`pathToModel_toElem`(`rfl`)を圏同値の counit で `𝒞` へ降ろすだけ。 -/
noncomputable def thm52iv_comp_toElem {G : Frobenioid P} (R : RatFnData P G)
    (hiso : ∀ Y : C, IsIsotropic P Y)
    (hfn : ∀ Z : BiratCat P G, IsFrobeniusNormalized (biratPre P G) Z)
    (Sb : BaseSection P) (Fs : ℕ+ →* SectionEnd Sb) (hFs : IsFrobeniusSection Sb Fs) :
    (thm_5_2_iv R hiso hfn Sb Fs hFs).functor ⋙ ModelData.toElem ≅ P.toElem := by
  haveI := pathForget_isEquivalence G Sb
  haveI := pathToModel_isEquivalence R hiso hfn Sb Fs hFs
  exact Functor.associator _ _ _ ≪≫
    Functor.isoWhiskerLeft (pathForget Sb).asEquivalence.inverse
      (eqToIso (pathToModel_toElem R hiso hfn Sb Fs hFs)) ≪≫
    (Functor.associator _ _ _).symm ≪≫
    Functor.isoWhiskerRight (pathForget Sb).asEquivalence.counitIso P.toElem ≪≫
    Functor.leftUnitor _

/-- ★★**`modelType_equiv` 版**。 -/
noncomputable def modelTypeEquiv_comp_toElem {G : Frobenioid P} (R : RatFnData P G)
    (hiso : ∀ Y : C, IsIsotropic P Y) (hm : IsOfModelType C P G) :
    (modelType_equiv R hiso hm).functor ⋙ ModelData.toElem ≅ P.toElem :=
  thm52iv_comp_toElem R hiso (fun Z => hm.2 Z) (Classical.choice hm.1).sec
    (Classical.choice hm.1).frob (Classical.choice hm.1).isFrobSection

/-! ## ★3. `Proposition 5.3` の縦の矢印は底への射影と 1-可換 -/

variable {S : Type} [CommSemiring S]

variable (S) in
/-- ★★★★★**`cToSc ⋙ (𝒞^rlf → 𝒟) ≅ P.proj`** ——
`ModelData.isRigidFunctor_of_proj` の仮定を満たす 1 本。

★3 段のうち 2 段は `rfl`(`modelDataHom_functor_comp_proj` と `cToUnTr_comp_proj`)、
残る 1 段が `Theorem 5.2, (iv)` の圏同値の counit である。 -/
noncomputable def cToSc_comp_proj (Fc : FrobenioidCore P) (G : Frobenioid P)
    (hiso : ∀ X : C, IsIsotropic P X)
    (hint : ∀ A : D, IsIntegralMonoid (Φ.val A))
    (hcharInj : ∀ {A B : D} (α : B ⟶ A),
      IsCharacteristicallyInjective (scMap (S := S) (Φ.map α)))
    (hintS : ∀ A : D, IsIntegralMonoid (ScT S (Φ.val A)))
    (hfsmD : IsOfFSMType D)
    (hM' : ModelData.Hyp (scModel S (unTr_frobenioid P Fc G) (unTr_isotropic P Fc)
      (fun Z => (unTr_isOfModelType Fc G).2 Z) hcharInj hintS hfsmD)) :
    cToSc (S := S) Fc G hiso hint hcharInj hintS hfsmD ⋙ (ModelData.modelPre hM').proj
      ≅ P.proj :=
  Functor.isoWhiskerLeft (cToUnTr hiso)
      (Functor.isoWhiskerRight (modelTypeEquiv_comp_toElem
        (unTr_ratFnData Fc G hint hfsmD) (unTr_isotropic P Fc) (unTr_isOfModelType Fc G))
        ElemFrobCat.proj)
    ≪≫ eqToIso (cToUnTr_comp_proj hiso)

variable (S) in
/-- ★★★★★★**[FrdI] Corollary 5.4 の rigidity** ——
`Proposition 5.3` の縦の矢印 `𝒞 ⥤ 𝒞^rlf` は **rigid**。

★★`Corollary 4.11, (i)` の道(充満性に依る)は使えないが、
model Frobenioid の 4 成分をつぶす道(`ModelData.isRigidFunctor_of_proj`)なら通る。
★`Div_B` の単射性は `scModel` では**包含そのもの**なので無料。

原文 (FrdI p.104):
> arrows are equivalences of categories]. Moreover, each of the composite functors of -/
theorem isRigidFunctor_cToSc (Fc : FrobenioidCore P) (G : Frobenioid P)
    (hiso : ∀ X : C, IsIsotropic P X)
    (hint : ∀ A : D, IsIntegralMonoid (Φ.val A))
    (hcharInj : ∀ {A B : D} (α : B ⟶ A),
      IsCharacteristicallyInjective (scMap (S := S) (Φ.map α)))
    (hintS : ∀ A : D, IsIntegralMonoid (ScT S (Φ.val A)))
    (hfsmD : IsOfFSMType D)
    (hM' : ModelData.Hyp (scModel S (unTr_frobenioid P Fc G) (unTr_isotropic P Fc)
      (fun Z => (unTr_isOfModelType Fc G).2 Z) hcharInj hintS hfsmD))
    (hrig : IsRigidFunctor P.proj) :
    IsRigidFunctor (cToSc (S := S) Fc G hiso hint hcharInj hintS hfsmD) :=
  ModelData.isRigidFunctor_of_proj hM' (fun _ => Subtype.coe_injective) _
    (isRigidFunctor_of_iso
      (cToSc_comp_proj S Fc G hiso hint hcharInj hintS hfsmD hM') hrig)

/-- ★★★★★★**合成 `𝒞₁ ⥤ 𝒞₂ ⥤ 𝒞₂^rlf` も rigid**(`Ψ` が圏同値のとき)。 -/
theorem isRigidFunctor_comp_cToSc {C₁ : Type u2} [Category.{v2} C₁] (Ψ : C₁ ⥤ C)
    [Ψ.IsEquivalence] (Fc : FrobenioidCore P) (G : Frobenioid P)
    (hiso : ∀ X : C, IsIsotropic P X)
    (hint : ∀ A : D, IsIntegralMonoid (Φ.val A))
    (hcharInj : ∀ {A B : D} (α : B ⟶ A),
      IsCharacteristicallyInjective (scMap (S := S) (Φ.map α)))
    (hintS : ∀ A : D, IsIntegralMonoid (ScT S (Φ.val A)))
    (hfsmD : IsOfFSMType D)
    (hM' : ModelData.Hyp (scModel S (unTr_frobenioid P Fc G) (unTr_isotropic P Fc)
      (fun Z => (unTr_isOfModelType Fc G).2 Z) hcharInj hintS hfsmD))
    (hrig : IsRigidFunctor P.proj) :
    IsRigidFunctor (Ψ ⋙ cToSc (S := S) Fc G hiso hint hcharInj hintS hfsmD) :=
  isRigidFunctor_comp_of_isEquivalence Ψ _
    (isRigidFunctor_cToSc S Fc G hiso hint hcharInj hintS hfsmD hM' hrig)

/-! ### ★出典の紐付け -/

/-- ★★★locator —— `Corollary 5.4` の rigidity の組み立て(★**条つき**)。 -/
def thm52iv_comp_toElem.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 104,
    item := "Corollary 5.4 — 縦の矢印は 𝔽_Φ への関手と 1-可換(rigidity の組み立て)",
    sectionId := "frdi-cor-5-4" }

/-! ## ★4. 1-一意性の後半 —— **同型の一意性**

原文 (FrdI p.104):
> arrows are equivalences of categories]. Moreover, each of the composite functors of

★★Mochizuki の「1-unique」は「**一意な**同型を除いて一意」の意味である。
★その**後半(同型の一意性)は rigidity からただちに出る** ——
`IsRigidFunctor F` は「`F ≅ F` は恒等だけ」なので、
2 つの同型 `α β : G ≅ F` に対し `α.symm ≪≫ β : F ≅ F` が恒等になり `α = β`。 -/

universe uu1 vv1 uu2 vv2

/-- ★★★**rigid な関手への同型は一意**。 -/
theorem iso_unique_of_rigid {E₁ : Type uu1} [Category.{vv1} E₁] {E₂ : Type uu2}
    [Category.{vv2} E₂] {G H : E₁ ⥤ E₂} (hH : IsRigidFunctor H) (α β : G ≅ H) : α = β := by
  have h : α.symm ≪≫ β = Iso.refl H := hH _
  have hh : α.inv ≫ β.hom = 𝟙 H := congrArg Iso.hom h
  refine Iso.ext ?_
  have h2 := congrArg (fun t => α.hom ≫ t) hh
  simp only [← Category.assoc, Iso.hom_inv_id, Category.id_comp, Category.comp_id] at h2
  exact h2.symm

/-- ★★同じく、rigid な関手**から**の同型も一意。 -/
theorem iso_unique_of_rigid' {E₁ : Type uu1} [Category.{vv1} E₁] {E₂ : Type uu2}
    [Category.{vv2} E₂] {G H : E₁ ⥤ E₂} (hG : IsRigidFunctor G) (α β : G ≅ H) : α = β := by
  have h := iso_unique_of_rigid hG α.symm β.symm
  have := congrArg Iso.symm h
  simpa using this

variable (S) in
/-- ★★★★★★**[FrdI] Corollary 5.4 の 1-可換性の同型は一意**。

★合成 `𝒞₁ ⥤ 𝒞₂ ⥤ 𝒞₂^rlf` が rigid(`isRigidFunctor_comp_cToSc`)なので、
そこへの同型は一意である。

原文 (FrdI p.104):
> arrows are equivalences of categories]. Moreover, each of the composite functors of -/
theorem cor54_comm_iso_unique {C₁ : Type u2} [Category.{v2} C₁] (Ψ : C₁ ⥤ C)
    [Ψ.IsEquivalence] (Fc : FrobenioidCore P) (G : Frobenioid P)
    (hiso : ∀ X : C, IsIsotropic P X)
    (hint : ∀ A : D, IsIntegralMonoid (Φ.val A))
    (hcharInj : ∀ {A B : D} (α : B ⟶ A),
      IsCharacteristicallyInjective (scMap (S := S) (Φ.map α)))
    (hintS : ∀ A : D, IsIntegralMonoid (ScT S (Φ.val A)))
    (hfsmD : IsOfFSMType D)
    (hM' : ModelData.Hyp (scModel S (unTr_frobenioid P Fc G) (unTr_isotropic P Fc)
      (fun Z => (unTr_isOfModelType Fc G).2 Z) hcharInj hintS hfsmD))
    (hrig : IsRigidFunctor P.proj)
    {H : C₁ ⥤ ScModelObj S (unTr_frobenioid P Fc G) (unTr_isotropic P Fc)
      (fun Z => (unTr_isOfModelType Fc G).2 Z) hcharInj hintS hfsmD}
    (α β : H ≅ Ψ ⋙ cToSc (S := S) Fc G hiso hint hcharInj hintS hfsmD) : α = β :=
  iso_unique_of_rigid
    (isRigidFunctor_comp_cToSc (S := S) Ψ Fc G hiso hint hcharInj hintS hfsmD hM' hrig) α β

variable (S) in
/-- ★★★★★**縦の矢印への同型も一意**(`𝒞 ⥤ 𝒞^rlf` が rigid だから)。 -/
theorem cor54_vert_iso_unique (Fc : FrobenioidCore P) (G : Frobenioid P)
    (hiso : ∀ X : C, IsIsotropic P X)
    (hint : ∀ A : D, IsIntegralMonoid (Φ.val A))
    (hcharInj : ∀ {A B : D} (α : B ⟶ A),
      IsCharacteristicallyInjective (scMap (S := S) (Φ.map α)))
    (hintS : ∀ A : D, IsIntegralMonoid (ScT S (Φ.val A)))
    (hfsmD : IsOfFSMType D)
    (hM' : ModelData.Hyp (scModel S (unTr_frobenioid P Fc G) (unTr_isotropic P Fc)
      (fun Z => (unTr_isOfModelType Fc G).2 Z) hcharInj hintS hfsmD))
    (hrig : IsRigidFunctor P.proj)
    {H : C ⥤ ScModelObj S (unTr_frobenioid P Fc G) (unTr_isotropic P Fc)
      (fun Z => (unTr_isOfModelType Fc G).2 Z) hcharInj hintS hfsmD}
    (α β : H ≅ cToSc (S := S) Fc G hiso hint hcharInj hintS hfsmD) : α = β :=
  iso_unique_of_rigid
    (isRigidFunctor_cToSc (S := S) Fc G hiso hint hcharInj hintS hfsmD hM' hrig) α β

/-- ★★★★★locator —— `Corollary 5.4` の 1-一意性(同型の一意性)。 -/
def cor54_comm_iso_unique.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 104,
    item := "Corollary 5.4 — 1-可換性の同型は一意(rigidity から)",
    sectionId := "frdi-cor-5-4" }

end ABC3.Found.FrdI
