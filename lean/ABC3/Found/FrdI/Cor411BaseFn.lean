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

end ABC3.Found.FrdI
