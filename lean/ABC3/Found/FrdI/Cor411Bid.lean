/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Cor411Otimes

/-!
# [FrdI] Corollary 4.11, (ii) —— base-identity 自己射の輸送(一般形)

原文 (FrdI p.93):
> such that every endomorphism [of an object of Cibirat] induced by

## ★★★★★測って分かった一般形(2026-08-19)

`Cor411Otimes.lean` の `otimes_map_of_divSlim` は `𝒪^×` について書いたが、
★**同じ骨が「base-identity 自己射」一般について通る**。
しかも **Div-slim の試験に使う単系 `Φ₀` は、pre-Frobenioid の単系と別でよい** ——
原文 (ii) が `𝒞^birat`(単系は `0_𝒟`)の base-identity 自己射を
**もとの `Φᵢ`** で試験するのはこの形である。

★手筋(原文どおり):
1. `Base α = 𝟙` を pull-back に沿って一意に持ち上げる(`Proposition 1.11, (iii)`)。
2. その族はスライスの**自己同型**を定める(自然性は pull-back の単射性)。
3. `Φ₀` が恒等へ送る ＋ **Div-slim** ⟹ 恒等 ⟹ 終対象で `Base (Ψ α) = 𝟙`。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

section BaseIdPull

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D}

/-- ★★**base-identity 自己射の pull-back に沿った持ち上げ**。 -/
noncomputable def baseIdPull (Q : PreFrobenioid C Φ) {A B : C} (z : B ⟶ A)
    (hz : IsPullBack Q z) (α : A ⟶ A) (hα : IsBaseIdentity Q α) : B ⟶ B :=
  ((prop_1_11_iii Q z hz (α : End A) (𝟙 _) (by
    rw [show Q.Base (((α : End A)) : A ⟶ A) = Q.Base (𝟙 A) from hα, Q.Base_id,
      Category.comp_id, Category.id_comp])).choose : B ⟶ B)

theorem baseIdPull_base (Q : PreFrobenioid C Φ) {A B : C} (z : B ⟶ A)
    (hz : IsPullBack Q z) (α : A ⟶ A) (hα : IsBaseIdentity Q α) :
    Q.Base (baseIdPull Q z hz α hα) = 𝟙 _ :=
  ((prop_1_11_iii Q z hz (α : End A) (𝟙 _) (by
    rw [show Q.Base (((α : End A)) : A ⟶ A) = Q.Base (𝟙 A) from hα, Q.Base_id,
      Category.comp_id, Category.id_comp])).choose_spec.1).1

theorem baseIdPull_isBaseIdentity (Q : PreFrobenioid C Φ) {A B : C} (z : B ⟶ A)
    (hz : IsPullBack Q z) (α : A ⟶ A) (hα : IsBaseIdentity Q α) :
    IsBaseIdentity Q (baseIdPull Q z hz α hα) := by
  show Q.Base (baseIdPull Q z hz α hα) = Q.Base (𝟙 B)
  rw [Q.Base_id]
  exact baseIdPull_base Q z hz α hα

theorem baseIdPull_spec (Q : PreFrobenioid C Φ) {A B : C} (z : B ⟶ A)
    (hz : IsPullBack Q z) (α : A ⟶ A) (hα : IsBaseIdentity Q α) :
    z ≫ α = baseIdPull Q z hz α hα ≫ z :=
  ((prop_1_11_iii Q z hz (α : End A) (𝟙 _) (by
    rw [show Q.Base (((α : End A)) : A ⟶ A) = Q.Base (𝟙 A) from hα, Q.Base_id,
      Category.comp_id, Category.id_comp])).choose_spec.1).2

/-- ★★★**自然性** —— pull-back の単射性から。 -/
theorem baseIdPull_natural (Q : PreFrobenioid C Φ) {A : C} (α : A ⟶ A)
    (hα : IsBaseIdentity Q α) {Z W : C} (f : Z ⟶ W) (w : W ⟶ A) (hw : IsPullBack Q w)
    (hz : IsPullBack Q (f ≫ w)) :
    baseIdPull Q (f ≫ w) hz α hα ≫ f = f ≫ baseIdPull Q w hw α hα := by
  have hsZ := baseIdPull_spec Q (f ≫ w) hz α hα
  have hsW := baseIdPull_spec Q w hw α hα
  have hcomp : (baseIdPull Q (f ≫ w) hz α hα ≫ f) ≫ w
      = (f ≫ baseIdPull Q w hw α hα) ≫ w := by
    calc (baseIdPull Q (f ≫ w) hz α hα ≫ f) ≫ w
        = baseIdPull Q (f ≫ w) hz α hα ≫ (f ≫ w) := by rw [Category.assoc]
      _ = (f ≫ w) ≫ α := hsZ.symm
      _ = f ≫ (w ≫ α) := by rw [Category.assoc]
      _ = f ≫ (baseIdPull Q w hw α hα ≫ w) := by rw [hsW]
      _ = (f ≫ baseIdPull Q w hw α hα) ≫ w := by rw [Category.assoc]
  have hbase : Q.Base (baseIdPull Q (f ≫ w) hz α hα ≫ f)
      = Q.Base (f ≫ baseIdPull Q w hw α hα) := by
    rw [Q.Base_comp, Q.Base_comp, baseIdPull_base, baseIdPull_base,
      Category.id_comp, Category.comp_id]
  exact (hw Z).1 (Subtype.ext (Prod.ext hcomp hbase))

theorem baseIdPull_congr (Q : PreFrobenioid C Φ) {A B : C} {z z' : B ⟶ A} (h : z = z')
    (hz : IsPullBack Q z) (hz' : IsPullBack Q z') (α : A ⟶ A) (hα : IsBaseIdentity Q α) :
    baseIdPull Q z hz α hα = baseIdPull Q z' hz' α hα := by
  subst h; rfl

def baseIdPull.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 93,
    item := "Corollary 4.11, (ii) — base-identity 自己射の pull-back に沿った持ち上げ",
    sectionId := "frdi-cor-4-11" }

end BaseIdPull

/-! ## ★2. スライスの自己同型と Div-slim -/

section BidTransport

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D}
  {D₂ : Type u} [Category.{v} D₂] {C₂ : Type u2} [Category.{v2} C₂]
  {Φ' : MonoidOn.{v, u, w} D₂}

/-- ★スライスの各対象での成分。 -/
noncomputable def baseIdSliceComp (Q : PreFrobenioid C Φ) (Q₂ : PreFrobenioid C₂ Φ')
    (Ψ : C ⥤ C₂) {A : C} (α : A ⟶ A) (hα : IsBaseIdentity Q α)
    (Z : Over (⟨A⟩ : PlBk Q)) :
    (Q₂.toElem.obj (Ψ.obj Z.left.obj)).base ⟶ (Q₂.toElem.obj (Ψ.obj Z.left.obj)).base :=
  Q₂.Base (Ψ.map (baseIdPull Q Z.hom.hom Z.hom.property α hα))

/-- ★★★★**スライスの自己同型**。 -/
noncomputable def baseIdSliceIso (Q : PreFrobenioid C Φ) (Q₂ : PreFrobenioid C₂ Φ')
    (Ψ : C ⥤ C₂)
    (hPB : ∀ {X Y : C} (f : X ⟶ Y), IsPullBack Q f → IsPullBack Q₂ (Ψ.map f))
    (hBI : ∀ {X : C} (β : X ⟶ X), IsBaseIdentity Q β → IsIso (Q₂.Base (Ψ.map β)))
    {A : C} (α : A ⟶ A) (hα : IsBaseIdentity Q α) :
    (plBkSlicePsi Q Q₂ Ψ hPB A ⋙ plBkOverFunctor Q₂ (Ψ.obj A)
      ⋙ Over.forget ((Q₂.toElem.obj (Ψ.obj A)).base))
    ≅ (plBkSlicePsi Q Q₂ Ψ hPB A ⋙ plBkOverFunctor Q₂ (Ψ.obj A)
      ⋙ Over.forget ((Q₂.toElem.obj (Ψ.obj A)).base)) :=
  NatIso.ofComponents
    (fun Z => @asIso _ _ _ _ (baseIdSliceComp Q Q₂ Ψ α hα Z)
      (hBI _ (baseIdPull_isBaseIdentity Q Z.hom.hom Z.hom.property α hα)))
    (fun {Z W} f => by
      have hw : f.left.hom ≫ W.hom.hom = Z.hom.hom :=
        congrArg InducedWideCategory.Hom.hom (Over.w f)
      have hz : IsPullBack Q (f.left.hom ≫ W.hom.hom) := by
        rw [hw]; exact Z.hom.property
      have hcong : baseIdPull Q Z.hom.hom Z.hom.property α hα
          = baseIdPull Q (f.left.hom ≫ W.hom.hom) hz α hα :=
        baseIdPull_congr Q hw.symm _ _ α hα
      have hnat := baseIdPull_natural Q α hα f.left.hom W.hom.hom W.hom.property hz
      show Q₂.Base (Ψ.map f.left.hom)
          ≫ Q₂.Base (Ψ.map (baseIdPull Q W.hom.hom W.hom.property α hα))
        = Q₂.Base (Ψ.map (baseIdPull Q Z.hom.hom Z.hom.property α hα))
          ≫ Q₂.Base (Ψ.map f.left.hom)
      rw [hcong, ← Q₂.Base_comp, ← Q₂.Base_comp, ← Ψ.map_comp, ← Ψ.map_comp, hnat])

set_option maxHeartbeats 1000000 in
/-- ★★★★★★**base-identity 自己射は Div-slim の下で `Ψ` に保たれる**(一般形)。

★試験に使う単系 `Φ₀` は `Q₂` の単系と**別でよい** ——
原文 (ii) が `𝒞^birat` の base-identity 自己射を `Φᵢ` で試験するのはこの形。 -/
theorem baseIdentity_map_of_divSlim (Q : PreFrobenioid C Φ) (F : FrobenioidCore Q)
    (Q₂ : PreFrobenioid C₂ Φ') (F₂ : FrobenioidCore Q₂) (Ψ : C ≌ C₂)
    (hPB : ∀ {X Y : C} (f : X ⟶ Y), IsPullBack Q f → IsPullBack Q₂ (Ψ.functor.map f))
    (hPB' : ∀ {X Y : C} (f : X ⟶ Y), IsPullBack Q₂ (Ψ.functor.map f) → IsPullBack Q f)
    {Φ₀ : MonoidOn.{v, u, w} D₂} (hds : IsDivSlim Φ₀)
    (hBI : ∀ {X : C} (β : X ⟶ X), IsBaseIdentity Q β →
      IsIso (Q₂.Base (Ψ.functor.map β)))
    (hΦ : ∀ {X : C} (β : X ⟶ X), IsBaseIdentity Q β →
      ∀ x : Φ₀.val (Q₂.proj.obj (Ψ.functor.obj X)),
        Φ₀.map (Q₂.Base (Ψ.functor.map β)) x = x)
    {A : C} (α : A ⟶ A) (hα : IsBaseIdentity Q α) :
    IsBaseIdentity Q₂ (Ψ.functor.map α) := by
  haveI := F₂.plBkEquiv (Ψ.functor.obj A)
  haveI := plBkPsi_isEquivalence Q Q₂ Ψ.functor hPB hPB'
  haveI : (plBkSlicePsi Q Q₂ Ψ.functor hPB A).IsEquivalence :=
    inferInstanceAs (Functor.IsEquivalence (Over.post (X := (⟨A⟩ : PlBk Q))
      (plBkPsi Q Q₂ Ψ.functor hPB)))
  haveI : Functor.IsEquivalence (plBkSlicePsi Q Q₂ Ψ.functor hPB A
      ⋙ plBkOverFunctor Q₂ (Ψ.functor.obj A)) := inferInstance
  set e := (plBkSlicePsi Q Q₂ Ψ.functor hPB A
    ⋙ plBkOverFunctor Q₂ (Ψ.functor.obj A)).asEquivalence with he
  set Gf := Over.forget ((Q₂.toElem.obj (Ψ.functor.obj A)).base) with hGf
  set θ : e.functor ⋙ Gf ≅ e.functor ⋙ Gf :=
    baseIdSliceIso Q Q₂ Ψ.functor hPB hBI α hα with hθ
  set X : e.inverse ⋙ e.functor ⋙ Gf ≅ Gf := e.invFunIdAssoc Gf with hX
  have hcomp : ∀ (Z : Over (⟨A⟩ : PlBk Q)) (x : Φ₀.val ((e.functor ⋙ Gf).obj Z)),
      Φ₀.map (θ.hom.app Z) x = x := fun Z x =>
    hΦ _ (baseIdPull_isBaseIdentity Q Z.hom.hom Z.hom.property α hα) x
  have hconj : overPhiAut Φ₀ ((Q₂.toElem.obj (Ψ.functor.obj A)).base)
      (X.symm ≪≫ Functor.isoWhiskerLeft e.inverse θ ≪≫ X) = 1 := by
    refine overPhiAut_eq_one Φ₀ _ (fun Z x => ?_)
    show Φ₀.map (X.inv.app Z ≫ θ.hom.app (e.inverse.obj Z) ≫ X.hom.app Z) x = x
    have hkey := hcomp (e.inverse.obj Z) (Φ₀.map (X.hom.app Z) x)
    calc Φ₀.map (X.inv.app Z ≫ θ.hom.app (e.inverse.obj Z) ≫ X.hom.app Z) x
        = Φ₀.map (X.inv.app Z) (Φ₀.map (θ.hom.app (e.inverse.obj Z))
            (Φ₀.map (X.hom.app Z) x)) := by
            rw [Φ₀.map_comp, Φ₀.map_comp]
            rfl
      _ = Φ₀.map (X.inv.app Z) (Φ₀.map (X.hom.app Z) x) :=
            congrArg (fun t => Φ₀.map (X.inv.app Z) t) hkey
      _ = Φ₀.map (X.inv.app Z ≫ X.hom.app Z) x :=
            (Φ₀.map_comp (X.hom.app Z) (X.inv.app Z) x).symm
      _ = x := by rw [X.inv_hom_id_app, Φ₀.map_id]
  have h1 : (X.symm ≪≫ Functor.isoWhiskerLeft e.inverse θ ≪≫ X)
      = (1 : Aut Gf) := hds _ (hconj.trans (overPhiAut_one Φ₀ _).symm)
  have hθrefl : θ = Iso.refl _ := eq_refl_of_conj_eq_refl e θ h1
  set Z₀ : Over (⟨A⟩ : PlBk Q) := Over.mk (𝟙 (⟨A⟩ : PlBk Q)) with hZ₀
  have hZ₀c : θ.hom.app Z₀ = 𝟙 ((e.functor ⋙ Gf).obj Z₀) :=
    congrArg (fun t : (e.functor ⋙ Gf) ≅ (e.functor ⋙ Gf) => t.hom.app Z₀) hθrefl
  have hpull : baseIdPull Q Z₀.hom.hom Z₀.hom.property α hα = α := by
    have hs := baseIdPull_spec Q Z₀.hom.hom Z₀.hom.property α hα
    have hs2 : (𝟙 Z₀.left.obj) ≫ α
        = baseIdPull Q Z₀.hom.hom Z₀.hom.property α hα ≫ (𝟙 Z₀.left.obj) := hs
    rw [Category.id_comp, Category.comp_id] at hs2
    exact hs2.symm
  show Q₂.Base (Ψ.functor.map α) = Q₂.Base (𝟙 (Ψ.functor.obj A))
  rw [Q₂.Base_id, ← hpull]
  exact hZ₀c

def baseIdentity_map_of_divSlim.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 93,
    item := "Corollary 4.11, (ii) — base-identity 自己射の輸送",
    sectionId := "frdi-cor-4-11" }

end BidTransport

end ABC3.Found.FrdI
