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

/-! ## ★3. `coaPreProp` の辞書 —— `𝒞^istr` と `𝒞` で一致する

★★`𝒞^istr` では co-angular が**自動**(`istr_coAngular`)であり、
`𝒞` の側でも **isotropic 対象から出る射は co-angular**
(`isCoAngular_of_isotropic_dom`)。★したがって両側で co-angular 性は消え、
残るのは `IsPreStep`(= linear ∧ base-isomorphism)だけで、
これは `degFr` と `Base` が一致するので `Iff.rfl` である。 -/

variable (F : FrobenioidCore P)

include F in
/-- ★`𝒞^istr` の pre-step は `𝒞` の pre-step —— `degFr` と `Base` が一致するので定義的。 -/
theorem istr_isPreStep_iff {X Y : Istr P} (g : X ⟶ Y) :
    IsPreStep (istrPre P F) g ↔ IsPreStep P g.hom := Iff.rfl

include F in
/-- ★★★**`coaPreProp` の辞書** —— 双有理化の添字圏が両側で一致する根拠。 -/
theorem istr_coaPreProp_iff {X Y : Istr P} (g : X ⟶ Y) :
    coaPreProp (istrPre P F) g ↔ coaPreProp P g.hom :=
  ⟨fun h => ⟨isCoAngular_of_isotropic_dom P F X.property g.hom, h.2⟩,
   fun h => ⟨istr_coAngular P F g, h.2⟩⟩

include F G in
/-- ★★★★★**添字圏が一致する** —— `A` が isotropic なら、`A` への co-angular pre-step は
**域まで込めて `𝒞^istr` の中で尽きる**。

★★`IdxBirat P G A.obj` の対象は `A.obj` への co-angular pre-step であり、
本定理よりその**域はすべて isotropic**。したがって
`IdxBirat (istrPre P F) _ A` と対象が一対一に対応する。
★これが `ratstd-istr` の残り 2 条(birationally Frobenius-normalized 型・rational 型)を
`Rmk451Birat.lean` の手法で取れる根拠である。 -/
theorem coaPre_into_istr {A : Istr P} {Z : C} (φ : Z ⟶ A.obj) (h : coaPreProp P φ) :
    IsIsotropic P Z :=
  isIsotropic_dom_of_coaPreStep P G φ h.1 h.2 A.property

def istr_coaPreProp_iff.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 83,
    item := "Proposition 4.4, (iv) — coaPreProp は 𝒞^istr と 𝒞 で一致する",
    sectionId := "frdi-prop-4-4" }

/-! ## ★4. 添字圏の同値 —— `(𝒞^istr)^birat` と `𝒞^birat` の `Hom` が一致する根拠

★★`SliceA` は**細い圏**(`slice_hom_ext` —— 構造射が mono だから)なので、
同値であることの中身は**本質的全射性と充満性**だけである。 -/

/-- ★★`𝒞^istr` の `𝒞^{coa-pre}` から `𝒞` のそれへの関手(包含)。 -/
noncomputable def coaPreIstr :
    CoaPre (istrPre P F) (istr_frobenioid P F G) ⥤ CoaPre P G where
  obj X := ⟨X.obj.obj⟩
  map {X Y} f := ⟨f.hom.hom, (istr_coaPreProp_iff P F f.hom).mp f.property⟩
  map_id X := WideSubcategory.hom_ext _ rfl
  map_comp {X Y Z} f g := WideSubcategory.hom_ext _ rfl

/-- ★スライスの対応。 -/
noncomputable abbrev sliceIstr (A : Istr P) :
    SliceA (istrPre P F) (istr_frobenioid P F G) A ⥤ SliceA P G A.obj :=
  Over.post (coaPreIstr P G F)

/-- ★添字圏(反対圏)の対応。 -/
noncomputable abbrev idxIstr (A : Istr P) :
    IdxBirat (istrPre P F) (istr_frobenioid P F G) A ⥤ IdxBirat P G A.obj :=
  (sliceIstr P G F A).op

include F G in
/-- ★★★★★**本質的全射** —— `A` が isotropic なので、`A` への co-angular pre-step は
**域まで込めて `𝒞^istr` の中で尽きる**(`isIsotropic_dom_of_coaPreStep`)。 -/
theorem sliceIstr_essSurj (A : Istr P) : (sliceIstr P G F A).EssSurj where
  mem_essImage Z := by
    have hiso : IsIsotropic P Z.left.obj :=
      isIsotropic_dom_of_coaPreStep P G Z.hom.hom Z.hom.property.1 Z.hom.property.2 A.property
    let Z0 : Istr P := ⟨Z.left.obj, hiso⟩
    let g : Z0 ⟶ A := InducedCategory.homMk Z.hom.hom
    have hg : coaPreProp (istrPre P F) g := (istr_coaPreProp_iff P F g).mpr Z.hom.property
    let Zc : CoaPre (istrPre P F) (istr_frobenioid P F G) := ⟨Z0⟩
    let gc : Zc ⟶ coaPreObj (istrPre P F) (istr_frobenioid P F G) A := ⟨g, hg⟩
    exact ⟨Over.mk gc, ⟨Iso.refl _⟩⟩

include F G in
/-- ★★**充満** —— `Istr P` は充満部分圏なので射がそのまま持ち上がる。 -/
theorem sliceIstr_full (A : Istr P) : (sliceIstr P G F A).Full where
  map_surjective {Z W} u := by
    let h0 : Z.left.obj ⟶ W.left.obj := InducedCategory.homMk u.left.hom
    have hh : coaPreProp (istrPre P F) h0 := (istr_coaPreProp_iff P F h0).mpr u.left.property
    have hw : u.left.hom ≫ W.hom.hom.hom = Z.hom.hom.hom :=
      congrArg (fun t : ((sliceIstr P G F A).obj Z).left ⟶ coaPreObj P G A.obj =>
        t.hom) (Over.w u)
    refine ⟨Over.homMk (⟨h0, hh⟩ : Z.left ⟶ W.left) ?_, ?_⟩
    · exact WideSubcategory.hom_ext _ (InducedCategory.hom_ext hw)
    · exact slice_hom_ext _ _

include F G in
/-- ★**忠実** —— `SliceA` は細い圏なので自明。 -/
theorem sliceIstr_faithful (A : Istr P) : (sliceIstr P G F A).Faithful where
  map_injective _ := slice_hom_ext _ _

include F G in
/-- ★★★★★**添字圏の同値** —— これで `Hom^birat` が両側で一致する。 -/
theorem sliceIstr_isEquivalence (A : Istr P) : (sliceIstr P G F A).IsEquivalence where
  faithful := sliceIstr_faithful P G F A
  full := sliceIstr_full P G F A
  essSurj := sliceIstr_essSurj P G F A

include F G in
/-- ★添字圏(反対圏)でも同値。 -/
theorem idxIstr_isEquivalence (A : Istr P) : (idxIstr P G F A).IsEquivalence := by
  haveI := sliceIstr_isEquivalence P G F A
  exact inferInstanceAs ((sliceIstr P G F A).op.IsEquivalence)

def sliceIstr_isEquivalence.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 83,
    item := "Proposition 4.4, (iv) — 双有理化の添字圏は 𝒞^istr と 𝒞 で一致する",
    sectionId := "frdi-prop-4-4" }

end ABC3.Found.FrdI
