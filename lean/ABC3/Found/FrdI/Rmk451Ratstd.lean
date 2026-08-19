/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Rmk451Birat
import ABC3.Found.FrdI.Prop44Core

/-!
# [FrdI] `Remark 4.5.1` の (a) —— `𝒞^istr` の双有理化は `𝒞` のそれと同じ

原文 (FrdI p.86):
> that if C is of rationally standard type (respectively, of standard type), then so is

## ★★★★★見立てを覆した(2026-08-19)

`Rmk451Birat.lean` の ★9 は
「`𝒞^istr` は `𝒞` の真の充満部分圏なので比較は圏同値ではない ——
`A` が isotropic でも `A′ → A` が pre-step のとき `A′` が isotropic とは限らない」
と書いていた。★**これは pre-step 一般については正しいが、
双有理化の添字圏が使うのは co-angular pre-step であり、そちらでは成り立つ。**

## ★測った筋(3 本の在庫を繋ぐだけ)

`φ : Z ⟶ A` が **co-angular pre-step**で `A` が isotropic なら:

| 段 | 在庫 |
|---|---|
| `A` が `𝒞^birat` で isotropic | `birat_isIsotropic_iff` |
| `φ` は `𝒞^birat` で**同型** | `birat_isIso_of_coaPre`(`Proposition 4.4, (iv)` の中心) |
| 同型で isotropy が移る | `isIsotropic_of_iso` |
| `Z` が `𝒞` で isotropic | `birat_isIsotropic_iff` を戻す |

★★★したがって **isotropic 対象への co-angular pre-step の域は isotropic**であり、
`SliceA` は `𝒞` と `𝒞^istr` で**一致する**。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ) (G : Frobenioid P)

/-! ## ★1. ★★★★★isotropic 対象への co-angular pre-step の域は isotropic -/

include G in
/-- ★★★★★**isotropic 対象への co-angular pre-step の域は isotropic**。

原文 (FrdI p.83):
> An object of C maps to an isotropic object of Cbirat if and only if

★★`𝒞^birat` を経由する —— co-angular pre-step はそこで**同型**になるので、
isotropy がそのまま逆向きに戻ってくる。
★**pre-step 一般では成り立たない**(co-angularity が要る)。 -/
theorem isIsotropic_dom_of_coaPreStep {Z A : C} (φ : Z ⟶ A)
    (hc : IsCoAngular P φ) (hs : IsPreStep P φ) (hA : IsIsotropic P A) :
    IsIsotropic P Z := by
  have hiso : IsIso ((toBiratCat P G).map φ) := birat_isIso_of_coaPre φ hc hs
  letI := hiso
  have hAb : IsIsotropic (biratPre P G) (show BiratCat P G from A) :=
    (birat_isIsotropic_iff P G A).mpr hA
  have e : (show BiratCat P G from Z) ≅ (show BiratCat P G from A) :=
    @asIso _ _ _ _ ((toBiratCat P G).map φ) hiso
  have hZb : IsIsotropic (biratPre P G) (show BiratCat P G from Z) :=
    isIsotropic_of_iso (biratPre P G) e hAb
  exact (birat_isIsotropic_iff P G Z).mp hZb

def isIsotropic_dom_of_coaPreStep.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 83,
    item := "Proposition 4.4, (iv) — isotropic 対象への co-angular pre-step の域は isotropic",
    sectionId := "frdi-prop-4-4" }

/-! ## ★2. 系 —— `𝒞^istr` の双有理化の添字圏は `𝒞` のものと一致する

★`IdxBirat P G A = (SliceA P G A)ᵒᵖ = (Over (coaPreObj P G A))ᵒᵖ` の対象は
`A` への co-angular pre-step である。★上の定理により、`A` が isotropic なら
**その域はすべて isotropic**、すなわち `𝒞^istr` の側にも同じ対象がある。 -/

include G in
/-- ★★`A` が isotropic なら、`A` への co-angular pre-step は `𝒞^istr` の中に収まる。 -/
theorem coaPreProp_dom_isotropic {Z A : C} (hA : IsIsotropic P A) {φ : Z ⟶ A}
    (h : coaPreProp P φ) : IsIsotropic P Z :=
  isIsotropic_dom_of_coaPreStep P G φ h.1 h.2 hA

end ABC3.Found.FrdI
