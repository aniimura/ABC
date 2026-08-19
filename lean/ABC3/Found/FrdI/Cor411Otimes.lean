/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Cor411Base
import ABC3.Found.FrdI.Cor411Untr
import ABC3.Found.FrdI.Prop111

/-!
# [FrdI] Corollary 4.11, (i) —— `𝒪^×` は pull-back に沿って持ち上がる

原文 (FrdI p.92):
> an element f ∈ O×(A) determines an automorphism [cf. the proof of Proposition

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

end ABC3.Found.FrdI
