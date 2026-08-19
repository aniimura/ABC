/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Cor411Birat
import ABC3.Found.FrdI.Prop311
import ABC3.Found.FrdI.Prop48Gl

/-!
# [FrdI] Corollary 4.11, (ii) —— `(𝒞^birat)^un-tr` の型と `Ψ_Base`

原文 (FrdI p.94):
> alences of categories]. Since, moreover, the Frobenioids (Cibirat)un-tr are of isotropic,

## ★★★★★測って分かったこと(2026-08-19)

原文は `(𝒞ᵢ^birat)^un-tr` が **isotropic ＋ unit-trivial ＋ group-like** 型であることを
使って `Proposition 3.11, (iii)` を当てる。★3 条件はすべて在庫から出る:

| 条件 | 出どころ |
|---|---|
| isotropic | `unTr_isotropic`(`Prop33UnTr.lean`)——`𝒞^un-tr` の対象はすべて isotropic |
| unit-trivial | `unTr_unitTrivial`(`UnTr.lean`) |
| group-like | ★`𝒞^birat` の単系は **`0_𝒟`**(1 元)なので `isOfGroupLikeType_of_phiTrivial` |

★★`Proposition 3.11, (iii)` を丸ごと使うと `IsSlimCat D₂` を要求されるが、
**存在・1-可換・1-一意は `psiBase` / `psiBaseCommute` / `psiBaseUniq` で
slim 無しに取れる**(rigidity だけが slim を要る ——原文どおり)。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

section UnTrBiratType

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ} (G : Frobenioid P)

/-- ★★`(𝒞^birat)^un-tr` も group-like 型。 -/
theorem untrBirat_isOfGroupLikeType (Fc : FrobenioidCore (biratPre P G)) :
    IsOfGroupLikeType (unTrPre (biratPre P G) Fc) :=
  isOfGroupLikeType_of_phiTrivial (fun _ _ => Subsingleton.elim (α := PUnit.{w + 1}) _ _)

/-- ★★`(𝒞^birat)^un-tr` は isotropic 型。 -/
theorem untrBirat_isOfIsotropicType (Fc : FrobenioidCore (biratPre P G)) :
    IsOfIsotropicType (unTrPre (biratPre P G) Fc) :=
  fun B => unTr_isotropic (biratPre P G) Fc B

/-- ★★`(𝒞^birat)^un-tr` は unit-trivial 型。 -/
theorem untrBirat_isOfUnitTrivialType (Fc : FrobenioidCore (biratPre P G)) :
    IsOfUnitTrivialType (unTrPre (biratPre P G) Fc) :=
  fun B => unTr_unitTrivial (biratPre P G) Fc B

def untrBirat_isOfGroupLikeType.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 94,
    item := "Corollary 4.11, (ii) — (𝒞^birat)^un-tr の 3 つの型",
    sectionId := "frdi-cor-4-11" }

end UnTrBiratType

/-! ## ★2. `hbi` —— `(Ψ^birat)^un-tr` は base-identity 自己射を保つ -/

section PsiBaseCor

variable {D₁ : Type u} [Category.{v} D₁] {C₁ : Type u2} [Category.{v2} C₁]
  {Φ₁ : MonoidOn.{v, u, w} D₁} {P₁ : PreFrobenioid C₁ Φ₁} (G₁ : Frobenioid P₁)
  {D₂ : Type u} [Category.{v} D₂] {C₂ : Type u2} [Category.{v2} C₂]
  {Φ₂ : MonoidOn.{v, u, w} D₂} {P₂ : PreFrobenioid C₂ Φ₂} (G₂ : Frobenioid P₂)
  (Ψ : C₁ ⥤ C₂)
  (hfwd : ∀ {X Y : C₁} (f : X ⟶ Y), coaPreProp P₁ f → coaPreProp P₂ (Ψ.map f))

include hfwd in
set_option maxHeartbeats 1000000 in
/-- ★★★★★**`(Ψ^birat)^un-tr` は base-identity 自己射を保つ**(`Proposition 3.11, (iii)` の `hbi`)。 -/
theorem psiBiratUnTr_baseIdentity [Ψ.IsEquivalence]
    (Fc₁ : FrobenioidCore (biratPre P₁ G₁)) (Fc₂ : FrobenioidCore (biratPre P₂ G₂))
    (hiso₁ : IsOfIsotropicType P₁) (hiso₂ : IsOfIsotropicType P₂)
    (hfn₁ : ∀ X : BiratCat P₁ G₁, IsFrobeniusNormalized (biratPre P₁ G₁) X)
    (hfn₂ : ∀ X : BiratCat P₂ G₂, IsFrobeniusNormalized (biratPre P₂ G₂) X)
    (hds₂ : IsDivSlim Φ₂)
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
    {A : UnTr (biratPre P₁ G₁)} (u : A ⟶ A)
    (hu : IsBaseIdentity (unTrPre (biratPre P₁ G₁) Fc₁) u) :
    IsBaseIdentity (unTrPre (biratPre P₂ G₂) Fc₂)
      ((psiBiratUnTr G₁ G₂ Ψ hfwd Fc₁ Fc₂ hiso₁ hiso₂ hfn₁ hfn₂ hds₂ hbwd hPBb hPBb'
        hDivEq hPS hBIso).map u) := by
  induction u using Quotient.inductionOn with
  | _ α =>
    have hα : IsBaseIdentity (biratPre P₁ G₁) α := hu
    exact biratPsi_baseIdentity_map G₁ G₂ Ψ hfwd Fc₂ hds₂ hbwd hPBb hPBb' hDivEq hPS hBIso
      α hα

def psiBiratUnTr_baseIdentity.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 94,
    item := "Corollary 4.11, (ii) — (Ψ^birat)^un-tr は base-identity 自己射を保つ",
    sectionId := "frdi-cor-4-11" }

/-! ## ★3. `Ψ_Base` の構成と 1-可換図式 -/

variable (Fc₁ : FrobenioidCore (biratPre P₁ G₁)) (Fc₂ : FrobenioidCore (biratPre P₂ G₂))
  (hiso₁ : IsOfIsotropicType P₁) (hiso₂ : IsOfIsotropicType P₂)
  (hfn₁ : ∀ X : BiratCat P₁ G₁, IsFrobeniusNormalized (biratPre P₁ G₁) X)
  (hfn₂ : ∀ X : BiratCat P₂ G₂, IsFrobeniusNormalized (biratPre P₂ G₂) X)
  (hds₂ : IsDivSlim Φ₂)
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

include hfn₁ in
/-- ★`(𝒞₁^birat)^un-tr ≌ 𝒟₁ × 𝒩`(`Proposition 3.11, (i)` の帰結)。 -/
theorem untrBirat_toProdCat_isEquivalence :
    (toProdCat (unTrPre (biratPre P₁ G₁) Fc₁)).IsEquivalence :=
  toElem_comp_prod_isEquivalence (unTrPre (biratPre P₁ G₁) Fc₁)
    (unTr_frobenioidCore (biratPre P₁ G₁) Fc₁)
    (unTr_frobenioid (biratPre P₁ G₁) Fc₁
      (birat_frobenioid_of_frobNormalized P₁ G₁ hfn₁))
    (untrBirat_isOfIsotropicType G₁ Fc₁) (untrBirat_isOfUnitTrivialType G₁ Fc₁)
    (untrBirat_isOfGroupLikeType G₁ Fc₁)

include hfwd hiso₁ hiso₂ hfn₂ hds₂ hbwd hPBb hPBb' hDivEq hPS hBIso in
/-- ★★★★★★**[FrdI] Corollary 4.11, (ii) の `Ψ_Base`** ——
`(𝒞ᵢ^birat)^un-tr ≌ 𝒟ᵢ × 𝒩` を経由して作る。 -/
noncomputable def psiBaseBirat [Ψ.IsEquivalence] : D₁ ⥤ D₂ :=
  haveI := untrBirat_toProdCat_isEquivalence G₁ Fc₁ hfn₁
  psiBase (psiBiratUnTr G₁ G₂ Ψ hfwd Fc₁ Fc₂ hiso₁ hiso₂ hfn₁ hfn₂ hds₂ hbwd hPBb hPBb'
      hDivEq hPS hBIso)
    (unTrPre (biratPre P₁ G₁) Fc₁) (unTrPre (biratPre P₂ G₂) Fc₂)

include hfwd hiso₁ hiso₂ hfn₂ hds₂ hbwd hPBb hPBb' hDivEq hPS hBIso in
set_option maxHeartbeats 1000000 in
/-- ★★★★★★**1-可換図式**(`(𝒞^birat)^un-tr` の層で)。 -/
noncomputable def psiBaseBirat_commute [Ψ.IsEquivalence] :
    (unTrPre (biratPre P₁ G₁) Fc₁).proj
        ⋙ psiBaseBirat G₁ G₂ Ψ hfwd Fc₁ Fc₂ hiso₁ hiso₂ hfn₁ hfn₂ hds₂ hbwd hPBb hPBb'
            hDivEq hPS hBIso
      ≅ psiBiratUnTr G₁ G₂ Ψ hfwd Fc₁ Fc₂ hiso₁ hiso₂ hfn₁ hfn₂ hds₂ hbwd hPBb hPBb'
            hDivEq hPS hBIso
          ⋙ (unTrPre (biratPre P₂ G₂) Fc₂).proj :=
  haveI := untrBirat_toProdCat_isEquivalence G₁ Fc₁ hfn₁
  psiBaseCommute _ _ _
    (fun {A} u hu => psiBiratUnTr_baseIdentity G₁ G₂ Ψ hfwd Fc₁ Fc₂ hiso₁ hiso₂ hfn₁ hfn₂
      hds₂ hbwd hPBb hPBb' hDivEq hPS hBIso u hu)

def psiBaseBirat.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 94,
    item := "Corollary 4.11, (ii) — Ψ_Base の構成と 1-可換図式",
    sectionId := "frdi-cor-4-11" }

end PsiBaseCor


/-! ## ★4. 鎖 `𝒞 → 𝒞^birat → (𝒞^birat)^un-tr` と射影の一致 -/

section Chain

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ} (G : Frobenioid P)

/-- ★`𝒞^birat` はすべての対象が isotropic なので `Istr` へ持ち上がる。 -/
noncomputable def biratToIstr (hiso : IsOfIsotropicType P) :
    BiratCat P G ⥤ Istr (biratPre P G) where
  obj A := ⟨A, birat_isOfIsotropicType hiso A⟩
  map f := ObjectProperty.homMk f
  map_id _ := rfl
  map_comp _ _ := rfl

/-- ★★**鎖** `𝒞 → 𝒞^birat → (𝒞^birat)^un-tr`。 -/
noncomputable def biratUnTrChain (hiso : IsOfIsotropicType P) :
    C ⥤ UnTr (biratPre P G) :=
  toBiratCat P G ⋙ biratToIstr G hiso ⋙ istrToUnTr (biratPre P G)

/-- ★★★**鎖と射影の合成はもとの射影**。 -/
theorem biratUnTrChain_proj (hiso : IsOfIsotropicType P)
    (Fc : FrobenioidCore (biratPre P G)) :
    biratUnTrChain G hiso ⋙ (unTrPre (biratPre P G) Fc).proj = P.proj := by
  refine CategoryTheory.Functor.ext (fun A => rfl) (fun A B f => ?_)
  show biratBase (toHomBirat (P := P) (G := G) f) = eqToHom rfl ≫ P.Base f ≫ eqToHom rfl
  simp only [eqToHom_refl, Category.id_comp, Category.comp_id]
  exact biratBase_toHomBirat f

def biratUnTrChain.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 94,
    item := "Corollary 4.11, (ii) — 𝒞 → 𝒞^birat → (𝒞^birat)^un-tr の鎖",
    sectionId := "frdi-cor-4-11" }

end Chain

/-! ## ★5. 鎖の 1-可換性と `Corollary 4.11, (ii)` の図式 -/

section ChainSquare

variable {D₁ : Type u} [Category.{v} D₁] {C₁ : Type u2} [Category.{v2} C₁]
  {Φ₁ : MonoidOn.{v, u, w} D₁} {P₁ : PreFrobenioid C₁ Φ₁} (G₁ : Frobenioid P₁)
  {D₂ : Type u} [Category.{v} D₂] {C₂ : Type u2} [Category.{v2} C₂]
  {Φ₂ : MonoidOn.{v, u, w} D₂} {P₂ : PreFrobenioid C₂ Φ₂} (G₂ : Frobenioid P₂)
  (Ψ : C₁ ⥤ C₂)
  (hfwd : ∀ {X Y : C₁} (f : X ⟶ Y), coaPreProp P₁ f → coaPreProp P₂ (Ψ.map f))
  (Fc₁ : FrobenioidCore (biratPre P₁ G₁)) (Fc₂ : FrobenioidCore (biratPre P₂ G₂))
  (hiso₁ : IsOfIsotropicType P₁) (hiso₂ : IsOfIsotropicType P₂)
  (hfn₁ : ∀ X : BiratCat P₁ G₁, IsFrobeniusNormalized (biratPre P₁ G₁) X)
  (hfn₂ : ∀ X : BiratCat P₂ G₂, IsFrobeniusNormalized (biratPre P₂ G₂) X)
  (hds₂ : IsDivSlim Φ₂)
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

include hfwd hiso₁ hiso₂ hfn₁ hfn₂ hds₂ hbwd hPBb hPBb' hDivEq hPS hBIso in
set_option maxHeartbeats 1000000 in
/-- ★★★★★**鎖の 1-可換性**(厳密な等式)。

★★`𝒞 → 𝒞^birat` の四角形(`psiBiratCor_square`)と
`𝒞^istr → 𝒞^un-tr` の四角形(`psiUnTr_square`、どちらも**等式**)を貼り合わせるだけ。
★間の `biratToIstr` は対象に isotropic の証明を付けるだけなので `rfl` で通り抜ける。 -/
theorem biratUnTrChain_square [Ψ.IsEquivalence] :
    biratUnTrChain G₁ hiso₁
        ⋙ psiBiratUnTr G₁ G₂ Ψ hfwd Fc₁ Fc₂ hiso₁ hiso₂ hfn₁ hfn₂ hds₂ hbwd hPBb hPBb'
            hDivEq hPS hBIso
      = Ψ ⋙ biratUnTrChain G₂ hiso₂ := by
  rw [show biratUnTrChain G₁ hiso₁
        ⋙ psiBiratUnTr G₁ G₂ Ψ hfwd Fc₁ Fc₂ hiso₁ hiso₂ hfn₁ hfn₂ hds₂ hbwd hPBb hPBb'
            hDivEq hPS hBIso
      = toBiratCat P₁ G₁ ⋙ psiBiratCor G₁ G₂ Ψ hfwd
          ⋙ biratToIstr G₂ hiso₂ ⋙ istrToUnTr (biratPre P₂ G₂) from rfl,
    ← Functor.assoc (toBiratCat P₁ G₁) (psiBiratCor G₁ G₂ Ψ hfwd),
    ← psiBiratCor_square G₁ G₂ Ψ hfwd]
  rfl

include hfwd hiso₁ hiso₂ hfn₁ hfn₂ hds₂ hbwd hPBb hPBb' hDivEq hPS hBIso in
set_option maxHeartbeats 1000000 in
/-- ★★★★★★**[FrdI] Corollary 4.11, (ii) の 1-可換図式** ——
`P₁.proj ⋙ Ψ_Base ≅ Ψ ⋙ P₂.proj`。 -/
noncomputable def cor_4_11_ii_square [Ψ.IsEquivalence] :
    P₁.proj ⋙ psiBaseBirat G₁ G₂ Ψ hfwd Fc₁ Fc₂ hiso₁ hiso₂ hfn₁ hfn₂ hds₂ hbwd hPBb hPBb'
        hDivEq hPS hBIso
      ≅ Ψ ⋙ P₂.proj := by
  have e1 : biratUnTrChain G₁ hiso₁
      ⋙ ((unTrPre (biratPre P₁ G₁) Fc₁).proj
        ⋙ psiBaseBirat G₁ G₂ Ψ hfwd Fc₁ Fc₂ hiso₁ hiso₂ hfn₁ hfn₂ hds₂ hbwd hPBb hPBb'
            hDivEq hPS hBIso)
      = P₁.proj ⋙ psiBaseBirat G₁ G₂ Ψ hfwd Fc₁ Fc₂ hiso₁ hiso₂ hfn₁ hfn₂ hds₂ hbwd hPBb
          hPBb' hDivEq hPS hBIso := by
    rw [← Functor.assoc, biratUnTrChain_proj G₁ hiso₁ Fc₁]
  have e3 : (biratUnTrChain G₁ hiso₁
        ⋙ psiBiratUnTr G₁ G₂ Ψ hfwd Fc₁ Fc₂ hiso₁ hiso₂ hfn₁ hfn₂ hds₂ hbwd hPBb hPBb'
            hDivEq hPS hBIso)
      ⋙ (unTrPre (biratPre P₂ G₂) Fc₂).proj = Ψ ⋙ P₂.proj := by
    rw [biratUnTrChain_square G₁ G₂ Ψ hfwd Fc₁ Fc₂ hiso₁ hiso₂ hfn₁ hfn₂ hds₂ hbwd hPBb
      hPBb' hDivEq hPS hBIso, Functor.assoc, biratUnTrChain_proj G₂ hiso₂ Fc₂]
  exact eqToIso e1.symm
    ≪≫ Functor.isoWhiskerLeft (biratUnTrChain G₁ hiso₁)
        (psiBaseBirat_commute G₁ G₂ Ψ hfwd Fc₁ Fc₂ hiso₁ hiso₂ hfn₁ hfn₂ hds₂ hbwd hPBb
          hPBb' hDivEq hPS hBIso)
    ≪≫ eqToIso e3

def cor_4_11_ii_square.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 91,
    item := "Corollary 4.11, (ii) — Ψ_Base の 1-可換図式",
    sectionId := "frdi-cor-4-11" }

end ChainSquare

/-! ## ★6. 前合成が関手の同型を反射する(一般補題) -/

section Precomp

universe v₃ u₃

/-- ★★★★**本質的全射かつ充満な関手との前合成は関手の同型を反射する**。

★★`Proposition 3.11, (iii)` の `projPrecompIso` は
`𝒞 ≌ 𝒟 × 𝒩`(group-like 型)を要求するが、
`Corollary 4.11, (ii)` の 1-一意性は**一般の `𝒞`** で要る。
★`P.proj` は `baseSurj` で本質的全射、`Proposition 1.11, (i)` で充満なので、
この一般補題で足りる。 -/
noncomputable def natIsoOfWhiskerLeft {Cc : Type u2} [Category.{v2} Cc]
    {Dd : Type u} [Category.{v} Dd] {Ee : Type u₃} [Category.{v₃} Ee]
    (F : Cc ⥤ Dd) [F.EssSurj] [F.Full] {G G' : Dd ⥤ Ee}
    (h : F ⋙ G ≅ F ⋙ G') : G ≅ G' := by
  -- ★綴りを固定した成分族(`(F ⋙ G).obj` と `G.obj (F.obj _)` の食い違いを避ける)
  let h₀ : ∀ Z : Cc, G.obj (F.obj Z) ≅ G'.obj (F.obj Z) := fun Z => h.app Z
  have hh₀ : ∀ Z, (h₀ Z).hom = h.hom.app Z := fun _ => rfl
  have hn : ∀ {Z W : Cc} (g : Z ⟶ W),
      G.map (F.map g) ≫ (h₀ W).hom = (h₀ Z).hom ≫ G'.map (F.map g) := by
    intro Z W g
    rw [hh₀ W, hh₀ Z]
    exact h.hom.naturality g
  refine NatIso.ofComponents
    (fun X => (G.mapIso (F.objObjPreimageIso X)).symm ≪≫ h₀ (F.objPreimage X)
      ≪≫ G'.mapIso (F.objObjPreimageIso X)) (fun {X Y} f => ?_)
  set pX := F.objPreimage X with hpX
  set pY := F.objPreimage Y with hpY
  set eX := F.objObjPreimageIso X with heX
  set eY := F.objObjPreimageIso Y with heY
  have hnat : G.map (eX.hom ≫ f ≫ eY.inv) ≫ (h₀ pY).hom
      = (h₀ pX).hom ≫ G'.map (eX.hom ≫ f ≫ eY.inv) := by
    have hg := hn (F.preimage (eX.hom ≫ f ≫ eY.inv))
    rwa [F.map_preimage] at hg
  rw [G.map_comp, G.map_comp, G'.map_comp, G'.map_comp] at hnat
  show G.map f ≫ (G.map eY.inv ≫ (h₀ pY).hom ≫ G'.map eY.hom)
    = (G.map eX.inv ≫ (h₀ pX).hom ≫ G'.map eX.hom) ≫ G'.map f
  have hEY : G'.map eY.inv ≫ G'.map eY.hom = 𝟙 (G'.obj Y) := by
    rw [← Functor.map_comp, eY.inv_hom_id, G'.map_id]
  have hEX : G.map eX.hom ≫ G.map eX.inv = 𝟙 (G.obj (F.obj pX)) := by
    rw [← Functor.map_comp, eX.hom_inv_id, G.map_id]
  refine (cancel_epi (G.map eX.hom)).mp ?_
  simp only [← Category.assoc] at hnat ⊢
  rw [hnat]
  simp only [Category.assoc, hEY, hEX, Category.comp_id]
  exact (Category.id_comp _).symm

def natIsoOfWhiskerLeft.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 94,
    item := "Corollary 4.11, (ii) — 前合成が関手の同型を反射する",
    sectionId := "frdi-cor-4-11" }

end Precomp


/-! ## ★7. `Corollary 4.11, (ii)` の 1-一意性 -/

section Uniq

variable {D₁ : Type u} [Category.{v} D₁] {C₁ : Type u2} [Category.{v2} C₁]
  {Φ₁ : MonoidOn.{v, u, w} D₁} {P₁ : PreFrobenioid C₁ Φ₁} (G₁ : Frobenioid P₁)
  {D₂ : Type u} [Category.{v} D₂] {C₂ : Type u2} [Category.{v2} C₂]
  {Φ₂ : MonoidOn.{v, u, w} D₂} {P₂ : PreFrobenioid C₂ Φ₂} (G₂ : Frobenioid P₂)
  (Ψ : C₁ ⥤ C₂)
  (hfwd : ∀ {X Y : C₁} (f : X ⟶ Y), coaPreProp P₁ f → coaPreProp P₂ (Ψ.map f))
  (Fc₁ : FrobenioidCore (biratPre P₁ G₁)) (Fc₂ : FrobenioidCore (biratPre P₂ G₂))
  (hiso₁ : IsOfIsotropicType P₁) (hiso₂ : IsOfIsotropicType P₂)
  (hfn₁ : ∀ X : BiratCat P₁ G₁, IsFrobeniusNormalized (biratPre P₁ G₁) X)
  (hfn₂ : ∀ X : BiratCat P₂ G₂, IsFrobeniusNormalized (biratPre P₂ G₂) X)
  (hds₂ : IsDivSlim Φ₂)
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

include hfwd hiso₁ hiso₂ hfn₁ hfn₂ hds₂ hbwd hPBb hPBb' hDivEq hPS hBIso in
/-- ★★★★★★**[FrdI] Corollary 4.11, (ii) の 1-一意性** ——
図式に入る関手は `Ψ_Base` と同型である。

★★在庫の `projPrecompIsoGen`(`Thm34VBase.lean`、`Theorem 3.4, (v)` の 1-一意性で作った
「`P.proj` との前合成は関手の同型を反射する」)がそのまま効く ——
★`base_three_factor` の 3 分解が「`Base` の像で自然 ⟹ `𝒟` 全体で自然」を担う。 -/
noncomputable def cor_4_11_ii_uniq [Ψ.IsEquivalence] (F₁ : FrobenioidCore P₁)
    (Gf : D₁ ⥤ D₂) (hG : P₁.proj ⋙ Gf ≅ Ψ ⋙ P₂.proj) :
    Gf ≅ psiBaseBirat G₁ G₂ Ψ hfwd Fc₁ Fc₂ hiso₁ hiso₂ hfn₁ hfn₂ hds₂ hbwd hPBb hPBb'
      hDivEq hPS hBIso :=
  projPrecompIsoGen P₁ F₁
    (hG ≪≫ (cor_4_11_ii_square G₁ G₂ Ψ hfwd Fc₁ Fc₂ hiso₁ hiso₂ hfn₁ hfn₂ hds₂ hbwd hPBb
      hPBb' hDivEq hPS hBIso).symm)

def cor_4_11_ii_uniq.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 91,
    item := "Corollary 4.11, (ii) — Ψ_Base の 1-一意性",
    sectionId := "frdi-cor-4-11" }

include hfwd hiso₁ hiso₂ hfn₁ hfn₂ hds₂ hbwd hPBb hPBb' hDivEq hPS hBIso in
set_option maxHeartbeats 1000000 in
/-- ★★★★★★**[FrdI] Corollary 4.11, (ii) の「両方が圏同値」** ——
`Ψ_Base : 𝒟₁ ⥤ 𝒟₂` は圏同値。

★`hUE'`・`hbi'` は擬逆側の条(原文が (ii) で `Ψ` **とその擬逆**に課すもの)。 -/
noncomputable def cor_4_11_ii_equivalence [Ψ.IsEquivalence]
    [(psiBiratUnTr G₁ G₂ Ψ hfwd Fc₁ Fc₂ hiso₁ hiso₂ hfn₁ hfn₂ hds₂ hbwd hPBb hPBb'
      hDivEq hPS hBIso).IsEquivalence]
    (hbi' : ∀ {X : UnTr (biratPre P₂ G₂)} (u : X ⟶ X),
      IsBaseIdentity (unTrPre (biratPre P₂ G₂) Fc₂) u →
        IsBaseIdentity (unTrPre (biratPre P₁ G₁) Fc₁)
          ((psiBiratUnTr G₁ G₂ Ψ hfwd Fc₁ Fc₂ hiso₁ hiso₂ hfn₁ hfn₂ hds₂ hbwd hPBb hPBb'
            hDivEq hPS hBIso).inv.map u)) :
    D₁ ≌ D₂ := by
  haveI := untrBirat_toProdCat_isEquivalence G₁ Fc₁ hfn₁
  haveI := untrBirat_toProdCat_isEquivalence G₂ Fc₂ hfn₂
  exact psiBaseEquivalence _ (unTrPre (biratPre P₁ G₁) Fc₁) (unTrPre (biratPre P₂ G₂) Fc₂)
    (fun {A} u hu => psiBiratUnTr_baseIdentity G₁ G₂ Ψ hfwd Fc₁ Fc₂ hiso₁ hiso₂ hfn₁ hfn₂
      hds₂ hbwd hPBb hPBb' hDivEq hPS hBIso u hu) hbi'

include hfwd hiso₁ hiso₂ hfn₁ hfn₂ hds₂ hbwd hPBb hPBb' hDivEq hPS hBIso in
set_option maxHeartbeats 1000000 in
/-- ★★★★★★**[FrdI] Corollary 4.11, (ii) の rigidity** ——
`𝒟₂` が slim なら図式の**両方の合成関手**が rigid。

原文 (FrdI p.92):
> [where the vertical arrows are the natural projection functors; the horizontal -/
theorem cor_4_11_ii_rigid [Ψ.IsEquivalence] (F₂ : FrobenioidCore P₂)
    (hslim : IsSlimCat D₂) :
    IsRigidFunctor (Ψ ⋙ P₂.proj) ∧
      IsRigidFunctor (P₁.proj ⋙ psiBaseBirat G₁ G₂ Ψ hfwd Fc₁ Fc₂ hiso₁ hiso₂ hfn₁ hfn₂
        hds₂ hbwd hPBb hPBb' hDivEq hPS hBIso) := by
  have h1 : IsRigidFunctor (Ψ ⋙ P₂.proj) :=
    isRigidFunctor_comp_of_isEquivalence Ψ P₂.proj (prop_1_13_i_global P₂ F₂ hslim)
  exact ⟨h1, isRigidFunctor_of_iso
    (cor_4_11_ii_square G₁ G₂ Ψ hfwd Fc₁ Fc₂ hiso₁ hiso₂ hfn₁ hfn₂ hds₂ hbwd hPBb hPBb'
      hDivEq hPS hBIso) h1⟩

def cor_4_11_ii_equivalence.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 92,
    item := "Corollary 4.11, (ii) — Ψ_Base は圏同値・合成関手は rigid",
    sectionId := "frdi-cor-4-11" }

end Uniq

end ABC3.Found.FrdI
