/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Cor411Base
import ABC3.Found.FrdI.Cor411Untr
import ABC3.Found.FrdI.Prop111
import ABC3.Found.FrdI.Thm34VBase

/-!
# [FrdI] Corollary 4.11, (i) —— `𝒪^×` は pull-back に沿って持ち上がる

原文 (FrdI p.92):
> determines an automorphism [cf. the proof of Proposition

## ★★★★★測って分かった道(2026-08-19)

原文は `f ∈ 𝒪^×(A)` が `Aut((𝒞^pl-bk)_A → 𝒞)` の元を定めることから始める。
★我々の語彙では、それは **`𝒪^×` を pull-back に沿って一意に持ち上げる**ことである
(`Proposition 1.11, (iv)`、在庫の `otriPullHom`)。

★★`otriPullHom : 𝒪^▷(A) →* 𝒪^▷(B)` は**単系準同型**なので、
**単元は単元へ写る**(`IsUnit.map`)——`𝒪^×` の持ち上げはただで出る。

★★★**自然性**は pull-back の**単射性**から出る:
`z = f ≫ w` のとき `β_Z ≫ f` と `f ≫ β_W` は
`(- ≫ w, Base -)` が同じ値を取るので一致する。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

section OtimesPull

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ)

/-- ★★`𝒪^×` の元は `𝒪^▷` の中で可逆。 -/
theorem isUnit_otri_of_otimes {A : C} {δ : End A} (h : δ ∈ OTimes P A) :
    IsUnit (⟨δ, h.1⟩ : OTri P A) := by
  haveI : IsIso ((δ : A ⟶ A)) := (CategoryTheory.isUnit_iff_isIso δ).mp h.2
  have hb : P.Base ((δ : A ⟶ A)) = 𝟙 _ := by
    have h2 : P.Base ((δ : A ⟶ A)) = P.Base (𝟙 A) := h.1.1
    rw [h2, P.Base_id]
  have hinv : ((inv ((δ : A ⟶ A))) : End A) ∈ OTri P A := by
    refine ⟨?_, degFr_of_isIso P _⟩
    show P.Base (inv ((δ : A ⟶ A))) = P.Base (𝟙 A)
    have h1 : P.Base ((δ : A ⟶ A)) ≫ P.Base (inv ((δ : A ⟶ A))) = P.Base (𝟙 A) := by
      rw [← P.Base_comp, IsIso.hom_inv_id]
    rw [hb, Category.id_comp] at h1
    exact h1
  refine ⟨⟨⟨δ, h.1⟩, ⟨((inv ((δ : A ⟶ A))) : End A), hinv⟩, ?_, ?_⟩, rfl⟩
  · refine Subtype.ext ?_
    show ((inv ((δ : A ⟶ A))) : A ⟶ A) ≫ ((δ : A ⟶ A)) = 𝟙 A
    exact IsIso.inv_hom_id _
  · refine Subtype.ext ?_
    show ((δ : A ⟶ A)) ≫ ((inv ((δ : A ⟶ A))) : A ⟶ A) = 𝟙 A
    exact IsIso.hom_inv_id _

variable (F : FrobenioidCore P) (hiso : ∀ X : C, IsIsotropic P X)

include hiso in
/-- ★pull-back 射は co-angular(isotropic 型なので `Proposition 1.4, (i)`)。 -/
theorem coAngular_of_isotropic {A B : C} (φ : B ⟶ A) : IsCoAngular P φ :=
  prop_1_4_i P φ (fun Y _ => hiso Y)

include F hiso in
/-- ★★★**`𝒪^×` の pull-back に沿った持ち上げ**。 -/
noncomputable def otimesPull {A B : C} (φ : B ⟶ A) (hpb : IsPullBack P φ)
    (δ : OTimes P A) : OTri P B :=
  otriPullHom P F φ (coAngular_of_isotropic P hiso φ) (F.pullBackLB φ hpb).2
    ⟨(δ : End A), δ.2.1⟩

include F hiso in
/-- ★持ち上げの四角形。 -/
theorem otimesPull_spec {A B : C} (φ : B ⟶ A) (hpb : IsPullBack P φ) (δ : OTimes P A) :
    φ ≫ ((δ : End A) : A ⟶ A)
      = ((otimesPull P F hiso φ hpb δ : End B) : B ⟶ B) ≫ φ :=
  otriPull_spec P F φ (coAngular_of_isotropic P hiso φ) (F.pullBackLB φ hpb).2
    ⟨(δ : End A), δ.2.1⟩

include F hiso in
/-- ★★持ち上げも `𝒪^×` に入る(`otriPullHom` が単系準同型だから)。 -/
theorem otimesPull_mem {A B : C} (φ : B ⟶ A) (hpb : IsPullBack P φ) (δ : OTimes P A) :
    ((otimesPull P F hiso φ hpb δ : End B)) ∈ OTimes P B := by
  refine ⟨(otimesPull P F hiso φ hpb δ).2, ?_⟩
  have h := (isUnit_otri_of_otimes P δ.2).map (otriPullHom P F φ
    (coAngular_of_isotropic P hiso φ) (F.pullBackLB φ hpb).2)
  exact h.map (OTri P B).subtype

include F hiso in
/-- ★★★★**自然性** —— pull-back の**単射性**から出る。 -/
theorem otimesPull_natural {A : C} (δ : OTimes P A)
    {Z W : C} (f : Z ⟶ W) (w : W ⟶ A) (hw : IsPullBack P w)
    (hz : IsPullBack P (f ≫ w)) :
    ((otimesPull P F hiso (f ≫ w) hz δ : End Z) : Z ⟶ Z) ≫ f
      = f ≫ ((otimesPull P F hiso w hw δ : End W) : W ⟶ W) := by
  set βZ := ((otimesPull P F hiso (f ≫ w) hz δ : End Z) : Z ⟶ Z) with hβZ
  set βW := ((otimesPull P F hiso w hw δ : End W) : W ⟶ W) with hβW
  have hspecZ : (f ≫ w) ≫ ((δ : End A) : A ⟶ A) = βZ ≫ (f ≫ w) :=
    otimesPull_spec P F hiso (f ≫ w) hz δ
  have hspecW : w ≫ ((δ : End A) : A ⟶ A) = βW ≫ w :=
    otimesPull_spec P F hiso w hw δ
  -- ★2 つの射は `(- ≫ w, Base -)` で同じ値を取る
  have hcomp : (βZ ≫ f) ≫ w = (f ≫ βW) ≫ w := by
    calc (βZ ≫ f) ≫ w = βZ ≫ (f ≫ w) := by rw [Category.assoc]
      _ = (f ≫ w) ≫ ((δ : End A) : A ⟶ A) := hspecZ.symm
      _ = f ≫ (w ≫ ((δ : End A) : A ⟶ A)) := by rw [Category.assoc]
      _ = f ≫ (βW ≫ w) := by rw [hspecW]
      _ = (f ≫ βW) ≫ w := by rw [Category.assoc]
  have hbZ : P.Base βZ = 𝟙 _ := by
    have := (otimesPull P F hiso (f ≫ w) hz δ).2.1
    show P.Base βZ = 𝟙 _
    rw [show P.Base βZ = P.Base (𝟙 Z) from this, P.Base_id]
  have hbW : P.Base βW = 𝟙 _ := by
    have := (otimesPull P F hiso w hw δ).2.1
    show P.Base βW = 𝟙 _
    rw [show P.Base βW = P.Base (𝟙 W) from this, P.Base_id]
  have hbase : P.Base (βZ ≫ f) = P.Base (f ≫ βW) := by
    rw [P.Base_comp, P.Base_comp, hbZ, hbW, Category.id_comp, Category.comp_id]
  exact (hw Z).1 (Subtype.ext (Prod.ext hcomp hbase))

def otimesPull.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 92,
    item := "Corollary 4.11, (i) — 𝒪^× の pull-back に沿った持ち上げ",
    sectionId := "frdi-cor-4-11" }

end OtimesPull

/-! ## ★2. スライスの自然同型 -/

section PsiSlice

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ}
  {D₂ : Type u} [Category.{v} D₂] {C₂ : Type u2} [Category.{v2} C₂]
  {Φ₂ : MonoidOn.{v, u, w} D₂} {P₂ : PreFrobenioid C₂ Φ₂}

/-- ★`otimesPull` は射の等式で移る。 -/
theorem otimesPull_congr (P : PreFrobenioid C Φ) (F : FrobenioidCore P)
    (hiso : ∀ X : C, IsIsotropic P X) {A B : C} {φ φ' : B ⟶ A} (h : φ = φ')
    (hpb : IsPullBack P φ) (hpb' : IsPullBack P φ') (δ : OTimes P A) :
    otimesPull P F hiso φ hpb δ = otimesPull P F hiso φ' hpb' δ := by
  subst h; rfl

/-- ★★スライスの各対象での成分 —— `Base (Ψ (α の持ち上げ))`。 -/
noncomputable def otimesSliceComp (P : PreFrobenioid C Φ) (F : FrobenioidCore P)
    (hiso : ∀ X : C, IsIsotropic P X) (P₂ : PreFrobenioid C₂ Φ₂) (Ψ : C ⥤ C₂)
    {A : C} (δ : OTimes P A) (Z : Over (⟨A⟩ : PlBk P)) :
    (P₂.toElem.obj (Ψ.obj Z.left.obj)).base ⟶ (P₂.toElem.obj (Ψ.obj Z.left.obj)).base :=
  P₂.Base (Ψ.map (((otimesPull P F hiso Z.hom.hom Z.hom.property δ : End Z.left.obj)
    : Z.left.obj ⟶ Z.left.obj)))

/-- ★成分は同型 —— 逆射は `Base (Ψ (逆射))`。 -/
theorem otimesSliceComp_isIso (P : PreFrobenioid C Φ) (F : FrobenioidCore P)
    (hiso : ∀ X : C, IsIsotropic P X) (P₂ : PreFrobenioid C₂ Φ₂) (Ψ : C ⥤ C₂)
    {A : C} (δ : OTimes P A) (Z : Over (⟨A⟩ : PlBk P)) :
    IsIso (otimesSliceComp P F hiso P₂ Ψ δ Z) := by
  have hu : ((otimesPull P F hiso Z.hom.hom Z.hom.property δ : End Z.left.obj))
      ∈ OTimes P Z.left.obj := otimesPull_mem P F hiso Z.hom.hom Z.hom.property δ
  haveI hg : IsIso (((otimesPull P F hiso Z.hom.hom Z.hom.property δ : End Z.left.obj)
      : Z.left.obj ⟶ Z.left.obj)) := (CategoryTheory.isUnit_iff_isIso _).mp hu.2
  set g : Z.left.obj ⟶ Z.left.obj :=
    ((otimesPull P F hiso Z.hom.hom Z.hom.property δ : End Z.left.obj)
      : Z.left.obj ⟶ Z.left.obj) with hgdef
  refine ⟨P₂.Base (Ψ.map (inv g)), ?_, ?_⟩
  · show P₂.Base (Ψ.map g) ≫ P₂.Base (Ψ.map (inv g)) = 𝟙 _
    rw [← P₂.Base_comp, ← Ψ.map_comp, IsIso.hom_inv_id, CategoryTheory.Functor.map_id,
      P₂.Base_id]
  · show P₂.Base (Ψ.map (inv g)) ≫ P₂.Base (Ψ.map g) = 𝟙 _
    rw [← P₂.Base_comp, ← Ψ.map_comp, IsIso.inv_hom_id, CategoryTheory.Functor.map_id,
      P₂.Base_id]

/-- ★★★★**スライスの自然同型** ——
`(𝒞^pl-bk)_A --Ψ--> (𝒞₂^pl-bk)_{ΨA} --> 𝒟₂` の自己同型。 -/
noncomputable def otimesSliceIso (P : PreFrobenioid C Φ) (F : FrobenioidCore P)
    (hiso : ∀ X : C, IsIsotropic P X) (P₂ : PreFrobenioid C₂ Φ₂) (Ψ : C ⥤ C₂)
    (hPB : ∀ {X Y : C} (f : X ⟶ Y), IsPullBack P f → IsPullBack P₂ (Ψ.map f))
    {A : C} (δ : OTimes P A) :
    (plBkSlicePsi P P₂ Ψ hPB A ⋙ plBkOverFunctor P₂ (Ψ.obj A)
      ⋙ Over.forget ((P₂.toElem.obj (Ψ.obj A)).base))
    ≅ (plBkSlicePsi P P₂ Ψ hPB A ⋙ plBkOverFunctor P₂ (Ψ.obj A)
      ⋙ Over.forget ((P₂.toElem.obj (Ψ.obj A)).base)) :=
  NatIso.ofComponents
    (fun Z => @asIso _ _ _ _ (otimesSliceComp P F hiso P₂ Ψ δ Z)
      (otimesSliceComp_isIso P F hiso P₂ Ψ δ Z))
    (fun {Z W} f => by
      have hw : f.left.hom ≫ W.hom.hom = Z.hom.hom :=
        congrArg InducedWideCategory.Hom.hom (Over.w f)
      have hz : IsPullBack P (f.left.hom ≫ W.hom.hom) := by
        rw [hw]; exact Z.hom.property
      have hcong : otimesPull P F hiso Z.hom.hom Z.hom.property δ
          = otimesPull P F hiso (f.left.hom ≫ W.hom.hom) hz δ :=
        otimesPull_congr P F hiso hw.symm _ _ δ
      have hnat := otimesPull_natural P F hiso δ f.left.hom W.hom.hom W.hom.property hz
      show P₂.Base (Ψ.map f.left.hom)
          ≫ P₂.Base (Ψ.map (((otimesPull P F hiso W.hom.hom W.hom.property δ
              : End W.left.obj) : W.left.obj ⟶ W.left.obj)))
        = P₂.Base (Ψ.map (((otimesPull P F hiso Z.hom.hom Z.hom.property δ
              : End Z.left.obj) : Z.left.obj ⟶ Z.left.obj)))
          ≫ P₂.Base (Ψ.map f.left.hom)
      rw [hcong, ← P₂.Base_comp, ← P₂.Base_comp, ← Ψ.map_comp, ← Ψ.map_comp, hnat])

def otimesSliceIso.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 92,
    item := "Corollary 4.11, (i) — 𝒪^× が定めるスライスの自己同型",
    sectionId := "frdi-cor-4-11" }

end PsiSlice

/-! ## ★3. Div-slim を当てる —— `Ψ` は `𝒪^×` を保つ -/

section PsiOtimes

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ}
  {D₂ : Type u} [Category.{v} D₂] {C₂ : Type u2} [Category.{v2} C₂]
  {Φ₂ : MonoidOn.{v, u, w} D₂} {P₂ : PreFrobenioid C₂ Φ₂}

set_option maxHeartbeats 1000000 in
/-- ★★★★★★**[FrdI] Corollary 4.11, (i) の核** —— `Ψ` は `𝒪^×` を保つ。

原文 (FrdI p.93):
> and Div-identity automorphisms [cf. Theorem 4.2, (i); our assumption that the

★★手筋(原文どおり):
1. `δ ∈ 𝒪^×(A)` を pull-back に沿って持ち上げ、スライスの**自己同型**を作る(`otimesSliceIso`)。
2. `Ψ` で移すと各成分は **Div-identity**(`Theorem 4.2, (i)`)なので `Φ₂` が恒等へ送る。
3. `𝒟₂` が **Div-slim** なので自己同型は恒等 ⟹ 終対象で見て `Base (Ψ δ) = 𝟙`。 -/
theorem otimes_map_of_divSlim (P : PreFrobenioid C Φ) (F : FrobenioidCore P)
    (hiso : ∀ X : C, IsIsotropic P X) (P₂ : PreFrobenioid C₂ Φ₂) (F₂ : FrobenioidCore P₂)
    (Ψ : C ≌ C₂)
    (hPB : ∀ {X Y : C} (f : X ⟶ Y), IsPullBack P f → IsPullBack P₂ (Ψ.functor.map f))
    (hPB' : ∀ {X Y : C} (f : X ⟶ Y), IsPullBack P₂ (Ψ.functor.map f) → IsPullBack P f)
    (hds₂ : IsDivSlim Φ₂)
    (hDivId : ∀ {X : C} (α : X ⟶ X), IsDivIdentity P α →
      IsDivIdentity P₂ (Ψ.functor.map α))
    {A : C} (δ : End A) (h : δ ∈ OTimes P A) :
    ((Ψ.functor.map ((δ : A ⟶ A))) : End (Ψ.functor.obj A)) ∈ OTimes P₂ (Ψ.functor.obj A) := by
  haveI := F₂.plBkEquiv (Ψ.functor.obj A)
  haveI := plBkPsi_isEquivalence P P₂ Ψ.functor hPB hPB'
  haveI : (plBkSlicePsi P P₂ Ψ.functor hPB A).IsEquivalence :=
    inferInstanceAs (Functor.IsEquivalence (Over.post (X := (⟨A⟩ : PlBk P))
      (plBkPsi P P₂ Ψ.functor hPB)))
  haveI : Functor.IsEquivalence (plBkSlicePsi P P₂ Ψ.functor hPB A
      ⋙ plBkOverFunctor P₂ (Ψ.functor.obj A)) := inferInstance
  set δ' : OTimes P A := ⟨δ, h⟩ with hδ'
  set e := (plBkSlicePsi P P₂ Ψ.functor hPB A
    ⋙ plBkOverFunctor P₂ (Ψ.functor.obj A)).asEquivalence with he
  set Gf := Over.forget ((P₂.toElem.obj (Ψ.functor.obj A)).base) with hGf
  set θ : e.functor ⋙ Gf ≅ e.functor ⋙ Gf :=
    otimesSliceIso P F hiso P₂ Ψ.functor hPB δ' with hθ
  set X : e.inverse ⋙ e.functor ⋙ Gf ≅ Gf := e.invFunIdAssoc Gf with hX
  -- ★各成分は `Φ₂` が恒等へ送る(`Theorem 4.2, (i)`)
  have hcomp : ∀ (Z : Over (⟨A⟩ : PlBk P)) (x : Φ₂.val ((e.functor ⋙ Gf).obj Z)),
      Φ₂.map (θ.hom.app Z) x = x := by
    intro Z x
    have hbi : IsDivIdentity P
        (((otimesPull P F hiso Z.hom.hom Z.hom.property δ' : End Z.left.obj)
          : Z.left.obj ⟶ Z.left.obj)) := by
      show Φ.map (P.Base _) = Φ.map (P.Base (𝟙 _))
      rw [show P.Base (((otimesPull P F hiso Z.hom.hom Z.hom.property δ'
          : End Z.left.obj) : Z.left.obj ⟶ Z.left.obj)) = P.Base (𝟙 Z.left.obj) from
        (otimesPull P F hiso Z.hom.hom Z.hom.property δ').2.1]
    have hbi₂ := hDivId _ hbi
    have h2 : Φ₂.map (P₂.Base (Ψ.functor.map
        (((otimesPull P F hiso Z.hom.hom Z.hom.property δ' : End Z.left.obj)
          : Z.left.obj ⟶ Z.left.obj))))
        = Φ₂.map (P₂.Base (𝟙 (Ψ.functor.obj Z.left.obj))) := hbi₂
    have h3 := congrArg (fun t : Φ₂.val _ →+ Φ₂.val _ => t x) h2
    refine h3.trans ?_
    show Φ₂.map (P₂.Base (𝟙 (Ψ.functor.obj Z.left.obj))) x = x
    rw [P₂.Base_id, Φ₂.map_id]
  -- ★共役して Div-slim を当てる
  have hconj : overPhiAut Φ₂ ((P₂.toElem.obj (Ψ.functor.obj A)).base)
      (X.symm ≪≫ Functor.isoWhiskerLeft e.inverse θ ≪≫ X) = 1 := by
    refine overPhiAut_eq_one Φ₂ _ (fun Z x => ?_)
    show Φ₂.map (X.inv.app Z ≫ θ.hom.app (e.inverse.obj Z) ≫ X.hom.app Z) x = x
    have hkey := hcomp (e.inverse.obj Z) (Φ₂.map (X.hom.app Z) x)
    calc Φ₂.map (X.inv.app Z ≫ θ.hom.app (e.inverse.obj Z) ≫ X.hom.app Z) x
        = Φ₂.map (X.inv.app Z) (Φ₂.map (θ.hom.app (e.inverse.obj Z))
            (Φ₂.map (X.hom.app Z) x)) := by
            rw [Φ₂.map_comp, Φ₂.map_comp]
            rfl
      _ = Φ₂.map (X.inv.app Z) (Φ₂.map (X.hom.app Z) x) :=
            congrArg (fun t => Φ₂.map (X.inv.app Z) t) hkey
      _ = Φ₂.map (X.inv.app Z ≫ X.hom.app Z) x :=
            (Φ₂.map_comp (X.hom.app Z) (X.inv.app Z) x).symm
      _ = x := by rw [X.inv_hom_id_app, Φ₂.map_id]
  have h1 : (X.symm ≪≫ Functor.isoWhiskerLeft e.inverse θ ≪≫ X)
      = (1 : Aut Gf) := hds₂ _ (hconj.trans (overPhiAut_one Φ₂ _).symm)
  have hθrefl : θ = Iso.refl _ := eq_refl_of_conj_eq_refl e θ h1
  -- ★終対象で見る
  set Z₀ : Over (⟨A⟩ : PlBk P) := Over.mk (𝟙 (⟨A⟩ : PlBk P)) with hZ₀
  have hZ₀c : θ.hom.app Z₀ = 𝟙 ((e.functor ⋙ Gf).obj Z₀) :=
    congrArg (fun t : (e.functor ⋙ Gf) ≅ (e.functor ⋙ Gf) => t.hom.app Z₀) hθrefl
  have hpull : ((otimesPull P F hiso Z₀.hom.hom Z₀.hom.property δ' : End Z₀.left.obj)
      : Z₀.left.obj ⟶ Z₀.left.obj) = (δ : A ⟶ A) := by
    have hs := otimesPull_spec P F hiso Z₀.hom.hom Z₀.hom.property δ'
    have hs2 : (𝟙 Z₀.left.obj) ≫ ((δ' : End A) : A ⟶ A)
        = (((otimesPull P F hiso Z₀.hom.hom Z₀.hom.property δ' : End Z₀.left.obj)
            : Z₀.left.obj ⟶ Z₀.left.obj)) ≫ (𝟙 Z₀.left.obj) := hs
    rw [Category.id_comp, Category.comp_id] at hs2
    exact hs2.symm
  refine (mem_otimes_iff P₂ _).mpr ⟨?_, ?_⟩
  · have hu : IsIso ((δ : A ⟶ A)) := (CategoryTheory.isUnit_iff_isIso _).mp h.2
    haveI := hu
    haveI : IsIso (Ψ.functor.map ((δ : A ⟶ A))) := Ψ.functor.map_isIso _
    exact (CategoryTheory.isUnit_iff_isIso _).mpr inferInstance
  · show P₂.Base (Ψ.functor.map ((δ : A ⟶ A))) = P₂.Base (𝟙 (Ψ.functor.obj A))
    rw [P₂.Base_id, ← hpull]
    exact hZ₀c

def otimes_map_of_divSlim.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 93,
    item := "Corollary 4.11, (i) — Ψ は 𝒪^× を保つ",
    sectionId := "frdi-cor-4-11" }

end PsiOtimes

/-! ## ★4. `Ψ^un-tr` —— `Corollary 4.11, (i)` の関手 -/

section PsiUnTr

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ}
  {D₂ : Type u} [Category.{v} D₂] {C₂ : Type u2} [Category.{v2} C₂]
  {Φ₂ : MonoidOn.{v, u, w} D₂} {P₂ : PreFrobenioid C₂ Φ₂}

/-- ★★★★★★**[FrdI] Corollary 4.11, (i)** —— `Ψ^un-tr` の構成。

原文 (FrdI p.91):
> (i) There exists a 1-unique functor

★★**底の 1-可換図式を経由しない** ——`Div-slim` ＋ `Theorem 4.2, (i)` から
直に `𝒪^×` の保存が出るので、`Corollary 4.11, (ii)` に依存しない。 -/
noncomputable def psiUnTrOfDivSlim (P : PreFrobenioid C Φ) (F : FrobenioidCore P)
    (hiso : ∀ X : C, IsIsotropic P X) (P₂ : PreFrobenioid C₂ Φ₂) (F₂ : FrobenioidCore P₂)
    (Ψ : C ≌ C₂)
    (hPB : ∀ {X Y : C} (f : X ⟶ Y), IsPullBack P f → IsPullBack P₂ (Ψ.functor.map f))
    (hPB' : ∀ {X Y : C} (f : X ⟶ Y), IsPullBack P₂ (Ψ.functor.map f) → IsPullBack P f)
    (hds₂ : IsDivSlim Φ₂)
    (hDivId : ∀ {X : C} (α : X ⟶ X), IsDivIdentity P α →
      IsDivIdentity P₂ (Ψ.functor.map α))
    (h₁ : IsOfQuasiIsotropicType C P) (h₂ : IsOfQuasiIsotropicType C₂ P₂) :
    UnTr P ⥤ UnTr P₂ :=
  haveI := Ψ.isEquivalence_functor
  psiUnTr Ψ.functor h₁ h₂
    (fun α₁ α₂ hh => toElem_map_congr_of_otimes Ψ.functor F hiso
      (fun _X δ hδ => otimes_map_of_divSlim P F hiso P₂ F₂ Ψ hPB hPB' hds₂ hDivId δ hδ)
      α₁ α₂ hh)

/-- ★★★★**1-可換図式**(`rfl`)。 -/
theorem psiUnTrOfDivSlim_square (P : PreFrobenioid C Φ) (F : FrobenioidCore P)
    (hiso : ∀ X : C, IsIsotropic P X) (P₂ : PreFrobenioid C₂ Φ₂) (F₂ : FrobenioidCore P₂)
    (Ψ : C ≌ C₂)
    (hPB : ∀ {X Y : C} (f : X ⟶ Y), IsPullBack P f → IsPullBack P₂ (Ψ.functor.map f))
    (hPB' : ∀ {X Y : C} (f : X ⟶ Y), IsPullBack P₂ (Ψ.functor.map f) → IsPullBack P f)
    (hds₂ : IsDivSlim Φ₂)
    (hDivId : ∀ {X : C} (α : X ⟶ X), IsDivIdentity P α →
      IsDivIdentity P₂ (Ψ.functor.map α))
    (h₁ : IsOfQuasiIsotropicType C P) (h₂ : IsOfQuasiIsotropicType C₂ P₂) :
    haveI := Ψ.isEquivalence_functor
    psiIstr Ψ.functor P P₂ h₁ h₂ ⋙ istrToUnTr P₂
      = istrToUnTr P ⋙ psiUnTrOfDivSlim P F hiso P₂ F₂ Ψ hPB hPB' hds₂ hDivId h₁ h₂ :=
  rfl

def psiUnTrOfDivSlim.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 91,
    item := "Corollary 4.11, (i) — Ψ^un-tr（Div-slim からの構成）",
    sectionId := "frdi-cor-4-11" }

end PsiUnTr

end ABC3.Found.FrdI
