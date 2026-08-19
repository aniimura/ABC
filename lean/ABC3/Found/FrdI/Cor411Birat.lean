/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Cor411Bid
import ABC3.Found.FrdI.Cor410Birat
import ABC3.Found.FrdI.Thm49
import ABC3.Found.FrdI.Prop44Gl

/-!
# [FrdI] Corollary 4.11, (ii) —— `Ψ^birat` は base-identity 自己射を保つ

原文 (FrdI p.93):
> preserves the base-identity endomorphisms [hence, in particular, that

## ★★★★★測って分かった還元(2026-08-19)

`Proposition 4.4, (iv)` の辞書で
`IsBaseIdentity (biratPre P G) (HomBirat.mk Z φ) ↔ Base ζ = Base φ`
(`ζ := Z.unop.hom.hom` は co-angular pre-step)である。
★`biratBase (HomBirat.mk Z φ) = inv (Base ζ) ≫ Base φ` なので、
**`Φ₂` がこれを恒等へ送ること**は
`Φ₂.map (Base (Ψ ζ)) = Φ₂.map (Base (Ψ φ))`、すなわち
★★**`Ψ` が Div-equivalence を保つこと**(`Theorem 4.9` の系 `divEquivalent_map`)
に他ならない。

★★★あとは `Cor411Bid.lean` の一般形 `baseIdentity_map_of_divSlim` を
**試験単系 `Φ₀ := Φ₂`** で当てれば、Div-slim が base-identity へ上げてくれる。
-/

namespace ABC3.Found.FrdI

open CategoryTheory CategoryTheory.Limits

universe v u w u2 v2

variable {D₁ : Type u} [Category.{v} D₁] {C₁ : Type u2} [Category.{v2} C₁]
  {Φ₁ : MonoidOn.{v, u, w} D₁} {P₁ : PreFrobenioid C₁ Φ₁} (G₁ : Frobenioid P₁)
  {D₂ : Type u} [Category.{v} D₂] {C₂ : Type u2} [Category.{v2} C₂]
  {Φ₂ : MonoidOn.{v, u, w} D₂} {P₂ : PreFrobenioid C₂ Φ₂} (G₂ : Frobenioid P₂)

variable (Ψ : C₁ ⥤ C₂)
  (hfwd : ∀ {X Y : C₁} (f : X ⟶ Y), coaPreProp P₁ f → coaPreProp P₂ (Ψ.map f))

section BiratBaseId

include hfwd in
/-- ★★★**`Ψ^birat` が base-identity 自己射に与える底の射**の明示形。 -/
theorem biratPsi_base_mk {A : C₁} (Z : IdxBirat P₁ G₁ A) (φ : Z.unop.left.obj ⟶ A)
    (hz : IsIso (P₂.Base (Ψ.map Z.unop.hom.hom))) :
    (biratPre P₂ G₂).Base
        ((psiBiratCor G₁ G₂ Ψ hfwd).map (HomBirat.mk Z φ))
      = inv (P₂.Base (Ψ.map Z.unop.hom.hom)) ≫ P₂.Base (Ψ.map φ) := by
  show biratBase (biratPsiMap G₁ G₂ Ψ hfwd A A (HomBirat.mk Z φ)) = _
  rw [biratPsiMap_mk G₁ G₂ Ψ hfwd A A Z φ, biratBase_mk]
  show (haveI := ((idxBiratPsi G₁ G₂ Ψ hfwd A).obj Z).unop.hom.property.2.2
    inv (P₂.Base (Ψ.map Z.unop.hom.hom)) ≫ P₂.Base (Ψ.map φ)) = _
  congr 1

include hfwd in
/-- ★★★★★**`Φ₂` は `Ψ^birat` が base-identity 自己射に与える底の射を恒等へ送る**。

★これが `Corollary 4.11, (ii)` の要 ——`Theorem 4.9` の系 `divEquivalent_map` を使う。 -/
theorem biratPsi_baseId_divTrivial
    (hDivEq : ∀ {X Y : C₁} (f g : X ⟶ Y), DivEquivalent P₁ f g →
      DivEquivalent P₂ (Ψ.map f) (Ψ.map g))
    (hPS : ∀ {X Y : C₁} (f : X ⟶ Y), IsPreStep P₁ f → IsPreStep P₂ (Ψ.map f))
    {A : C₁} (β : (show BiratCat P₁ G₁ from A) ⟶ (show BiratCat P₁ G₁ from A))
    (hβ : IsBaseIdentity (biratPre P₁ G₁) β)
    (x : Φ₂.val ((biratPre P₂ G₂).proj.obj ((psiBiratCor G₁ G₂ Ψ hfwd).obj A))) :
    Φ₂.map ((biratPre P₂ G₂).Base ((psiBiratCor G₁ G₂ Ψ hfwd).map β)) x = x := by
  obtain ⟨Z, φ, rfl, hbe⟩ := (birat_isBaseIdentity_iff (P := P₁) (G := G₁) β).mp hβ
  -- `ζ` は co-angular pre-step、`Ψ ζ` も pre-step なので底は同型
  have hζs : IsPreStep P₁ Z.unop.hom.hom := Z.unop.hom.property.2
  haveI hz : IsIso (P₂.Base (Ψ.map Z.unop.hom.hom)) := (hPS _ hζs).2
  -- `Base ζ = Base φ` ⟹ Div-equivalent ⟹ `Ψ` で移る
  have hde : DivEquivalent P₁ Z.unop.hom.hom φ := by
    show Φ₁.map (P₁.Base Z.unop.hom.hom) = Φ₁.map (P₁.Base φ)
    rw [show P₁.Base Z.unop.hom.hom = P₁.Base φ from hbe]
    rfl
  have hde₂ : Φ₂.map (P₂.Base (Ψ.map Z.unop.hom.hom))
      = Φ₂.map (P₂.Base (Ψ.map φ)) := hDivEq _ _ hde
  refine Eq.trans (congrArg (fun t => Φ₂.map t x)
    (biratPsi_base_mk G₁ G₂ Ψ hfwd Z φ hz)) ?_
  have step1 : Φ₂.map (inv (P₂.Base (Ψ.map Z.unop.hom.hom)) ≫ P₂.Base (Ψ.map φ)) x
      = Φ₂.map (inv (P₂.Base (Ψ.map Z.unop.hom.hom))) (Φ₂.map (P₂.Base (Ψ.map φ)) x) :=
    Φ₂.map_comp _ _ x
  have step2 : Φ₂.map (P₂.Base (Ψ.map φ)) x
      = Φ₂.map (P₂.Base (Ψ.map Z.unop.hom.hom)) x :=
    congrArg (fun t : Φ₂.val _ →+ Φ₂.val _ => t x) hde₂.symm
  have step3 : Φ₂.map (inv (P₂.Base (Ψ.map Z.unop.hom.hom)))
      (Φ₂.map (P₂.Base (Ψ.map Z.unop.hom.hom)) x) = x :=
    (Φ₂.map_comp (P₂.Base (Ψ.map Z.unop.hom.hom))
      (inv (P₂.Base (Ψ.map Z.unop.hom.hom))) x).symm.trans (by
        rw [IsIso.inv_hom_id]; exact Φ₂.map_id _ x)
  exact step1.trans
    ((congrArg (fun t => Φ₂.map (inv (P₂.Base (Ψ.map Z.unop.hom.hom))) t) step2).trans step3)

include hfwd in
/-- ★★`Ψ^birat` が base-identity 自己射に与える底の射は**同型**。 -/
theorem biratPsi_baseId_isIso
    (hPS : ∀ {X Y : C₁} (f : X ⟶ Y), IsPreStep P₁ f → IsPreStep P₂ (Ψ.map f))
    (hBIso : ∀ {X Y : C₁} (f : X ⟶ Y), IsBaseIsomorphism P₁ f →
      IsBaseIsomorphism P₂ (Ψ.map f))
    {A : C₁} (β : (show BiratCat P₁ G₁ from A) ⟶ (show BiratCat P₁ G₁ from A))
    (hβ : IsBaseIdentity (biratPre P₁ G₁) β) :
    IsIso ((biratPre P₂ G₂).Base ((psiBiratCor G₁ G₂ Ψ hfwd).map β)) := by
  obtain ⟨Z, φ, rfl, hbe⟩ := (birat_isBaseIdentity_iff (P := P₁) (G := G₁) β).mp hβ
  have hζs : IsPreStep P₁ Z.unop.hom.hom := Z.unop.hom.property.2
  haveI hz : IsIso (P₂.Base (Ψ.map Z.unop.hom.hom)) := (hPS _ hζs).2
  have hφb : IsBaseIsomorphism P₁ φ := by
    show IsIso (P₁.Base φ)
    rw [← show P₁.Base Z.unop.hom.hom = P₁.Base φ from hbe]
    exact hζs.2
  haveI : IsIso (P₂.Base (Ψ.map φ)) := hBIso _ hφb
  have hbm := biratPsi_base_mk G₁ G₂ Ψ hfwd Z φ hz
  have hii : IsIso (inv (P₂.Base (Ψ.map Z.unop.hom.hom)) ≫ P₂.Base (Ψ.map φ)) :=
    inferInstance
  exact hbm ▸ hii

def biratPsi_baseId_divTrivial.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 93,
    item := "Corollary 4.11, (ii) — Ψ^birat と base-identity 自己射",
    sectionId := "frdi-cor-4-11" }

/-! ## ★2. 組み立て —— `Ψ^birat` は base-identity 自己射を保つ -/

include hfwd in
set_option maxHeartbeats 1000000 in
/-- ★★★★★★**[FrdI] Corollary 4.11, (ii) の段 1** ——
`Ψ^birat` は base-identity 自己射を保つ。

原文 (FrdI p.93):
> preserves the base-identity endomorphisms [hence, in particular, that -/
theorem biratPsi_baseIdentity_map [Ψ.IsEquivalence]
    (Fc₂ : FrobenioidCore (biratPre P₂ G₂)) (hds₂ : IsDivSlim Φ₂)
    (hbwd : ∀ {X Y : C₁} (f : X ⟶ Y), coaPreProp P₂ (Ψ.map f) → coaPreProp P₁ f)
    (hPBb : ∀ {X Y : BiratCat P₁ G₁} (f : X ⟶ Y), IsPullBack (biratPre P₁ G₁) f →
      IsPullBack (biratPre P₂ G₂) ((psiBiratCor G₁ G₂ Ψ hfwd).map f))
    (hPBb' : ∀ {X Y : BiratCat P₁ G₁} (f : X ⟶ Y),
      IsPullBack (biratPre P₂ G₂) ((psiBiratCor G₁ G₂ Ψ hfwd).map f) →
        IsPullBack (biratPre P₁ G₁) f)
    (hDivEq : ∀ {X Y : C₁} (f g : X ⟶ Y), DivEquivalent P₁ f g →
      DivEquivalent P₂ (Ψ.map f) (Ψ.map g))
    (hPS : ∀ {X Y : C₁} (f : X ⟶ Y), IsPreStep P₁ f → IsPreStep P₂ (Ψ.map f))
    (hBIso : ∀ {X Y : C₁} (f : X ⟶ Y), IsBaseIsomorphism P₁ f →
      IsBaseIsomorphism P₂ (Ψ.map f))
    {A : C₁} (β : (show BiratCat P₁ G₁ from A) ⟶ (show BiratCat P₁ G₁ from A))
    (hβ : IsBaseIdentity (biratPre P₁ G₁) β) :
    IsBaseIdentity (biratPre P₂ G₂) ((psiBiratCor G₁ G₂ Ψ hfwd).map β) := by
  haveI := psiBiratCor_isEquivalence G₁ G₂ Ψ hfwd hbwd
  exact baseIdentity_map_of_divSlim (biratPre P₁ G₁) (biratPre P₂ G₂) Fc₂
    (psiBiratCor G₁ G₂ Ψ hfwd).asEquivalence hPBb hPBb' hds₂
    (fun {X} β' hβ' => biratPsi_baseId_isIso G₁ G₂ Ψ hfwd hPS hBIso β' hβ')
    (fun {X} β' hβ' x => biratPsi_baseId_divTrivial G₁ G₂ Ψ hfwd hDivEq hPS β' hβ' x)
    β hβ

def biratPsi_baseIdentity_map.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 93,
    item := "Corollary 4.11, (ii) — Ψ^birat は base-identity 自己射を保つ",
    sectionId := "frdi-cor-4-11" }

include hfwd in
/-- ★★★★★**[FrdI] Corollary 4.11, (ii) の段 2** —— `Ψ^birat` は `𝒪^×` を保つ。

原文 (FrdI p.93):
> preserves the base-identity endomorphisms [hence, in particular, that -/
theorem biratPsi_otimes_map [Ψ.IsEquivalence]
    (Fc₂ : FrobenioidCore (biratPre P₂ G₂)) (hds₂ : IsDivSlim Φ₂)
    (hbwd : ∀ {X Y : C₁} (f : X ⟶ Y), coaPreProp P₂ (Ψ.map f) → coaPreProp P₁ f)
    (hPBb : ∀ {X Y : BiratCat P₁ G₁} (f : X ⟶ Y), IsPullBack (biratPre P₁ G₁) f →
      IsPullBack (biratPre P₂ G₂) ((psiBiratCor G₁ G₂ Ψ hfwd).map f))
    (hPBb' : ∀ {X Y : BiratCat P₁ G₁} (f : X ⟶ Y),
      IsPullBack (biratPre P₂ G₂) ((psiBiratCor G₁ G₂ Ψ hfwd).map f) →
        IsPullBack (biratPre P₁ G₁) f)
    (hDivEq : ∀ {X Y : C₁} (f g : X ⟶ Y), DivEquivalent P₁ f g →
      DivEquivalent P₂ (Ψ.map f) (Ψ.map g))
    (hPS : ∀ {X Y : C₁} (f : X ⟶ Y), IsPreStep P₁ f → IsPreStep P₂ (Ψ.map f))
    (hBIso : ∀ {X Y : C₁} (f : X ⟶ Y), IsBaseIsomorphism P₁ f →
      IsBaseIsomorphism P₂ (Ψ.map f))
    (X : BiratCat P₁ G₁) (δ : End X) (h : δ ∈ OTimes (biratPre P₁ G₁) X) :
    (((psiBiratCor G₁ G₂ Ψ hfwd).map ((δ : X ⟶ X)))
      : End ((psiBiratCor G₁ G₂ Ψ hfwd).obj X))
      ∈ OTimes (biratPre P₂ G₂) ((psiBiratCor G₁ G₂ Ψ hfwd).obj X) :=
  (mem_otimes_iff (biratPre P₂ G₂) _).mpr
    ⟨h.2.map (Functor.mapEnd X (psiBiratCor G₁ G₂ Ψ hfwd)),
     biratPsi_baseIdentity_map G₁ G₂ Ψ hfwd Fc₂ hds₂ hbwd hPBb hPBb' hDivEq hPS hBIso
       ((δ : X ⟶ X)) h.1.1⟩

def biratPsi_otimes_map.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 93,
    item := "Corollary 4.11, (ii) — Ψ^birat は 𝒪^× を保つ",
    sectionId := "frdi-cor-4-11" }

end BiratBaseId

end ABC3.Found.FrdI
