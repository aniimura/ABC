/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Rmk451Birat
import ABC3.Found.FrdI.Prop44Core
import ABC3.Found.FrdI.Prop44Phi
import ABC3.Found.FrdI.MonoidTransport
import Mathlib.CategoryTheory.Limits.Final

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

open CategoryTheory CategoryTheory.Limits
open scoped NNReal

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

/-! ## ★5. `Hom^birat` の写像 —— `Rmk451Birat.lean` の手法をなぞる

★★`homFunctorBirat` の対象は `ULift (Z.unop.left.obj ⟶ B)` であり、
`Istr P` の hom は `C` の hom と**定義的に同じ**(`InducedCategory`)なので、
余錐の成分は `.hom` を取るだけである。 -/

/-- ★★降ろすための余錐。自然性は `HomBirat.mk_map` 1 本。 -/
noncomputable def biratIstrCocone (A B : Istr P) :
    Cocone (homFunctorBirat (istrPre P F) (istr_frobenioid P F G) A B) :=
  Cocone.mk (HomBirat P G A.obj B.obj)
    { app := fun Z => TypeCat.ofHom fun φ =>
        HomBirat.mk ((idxIstr P G F A).obj Z) φ.down.hom
      naturality := fun Z W u => by
        ext φ
        show HomBirat.mk ((idxIstr P G F A).obj W) (u.unop.left.hom ≫ φ.down).hom
          = HomBirat.mk ((idxIstr P G F A).obj Z) φ.down.hom
        exact HomBirat.mk_map ((idxIstr P G F A).map u) φ.down.hom }

/-- ★★★`Hom^birat` の写像。 -/
noncomputable def biratIstrMap (A B : Istr P) :
    HomBirat (istrPre P F) (istr_frobenioid P F G) A B → HomBirat P G A.obj B.obj :=
  fun z => colimit.desc _ (biratIstrCocone P G F A B) z

/-- ★代表元での値 —— 以後の計算の入口。 -/
theorem biratIstrMap_mk (A B : Istr P)
    (Z : IdxBirat (istrPre P F) (istr_frobenioid P F G) A) (φ : Z.unop.left.obj ⟶ B) :
    biratIstrMap P G F A B (HomBirat.mk Z φ)
      = HomBirat.mk ((idxIstr P G F A).obj Z) φ.hom := by
  show colimit.desc _ (biratIstrCocone P G F A B)
      (colimit.ι (homFunctorBirat (istrPre P F) (istr_frobenioid P F G) A B) Z
        (ULift.up φ)) = _
  rw [← types_comp_apply
    (colimit.ι (homFunctorBirat (istrPre P F) (istr_frobenioid P F G) A B) Z)
    (colimit.desc _ (biratIstrCocone P G F A B)), colimit.ι_desc]
  rfl

/-! ## ★6. 余錐が余極限であること —— 添字圏の同値から

★★`idxIstr` は同値なので **final** であり、
`colimit.cocone` を whisker したものは余極限のままである。 -/

include F G in
/-- ★`idxIstr` は final(同値だから)。 -/
theorem idxIstr_final (A : Istr P) : (idxIstr P G F A).Final := by
  haveI := idxIstr_isEquivalence P G F A
  infer_instance

/-! ## ★6. 全単射 —— 添字圏の同値から直に出す

★★`HomBirat.exists_rep` / `HomBirat.mk_map` / `HomBirat.eq_iff` の 3 本で足りる。
★余極限の一般論(`Functor.Final`)を経由するより、**代表元で書く方が短い**。 -/

include F G in
/-- ★★★**全射** —— 添字を `idxIstr` の原像へ引き戻し、`mk_map` で移す。 -/
theorem biratIstrMap_surjective (A B : Istr P) :
    Function.Surjective (biratIstrMap P G F A B) := by
  haveI := idxIstr_isEquivalence P G F A
  intro z
  obtain ⟨Z, φ, rfl⟩ := HomBirat.exists_rep z
  let Z' := (idxIstr P G F A).objPreimage Z
  let e : (idxIstr P G F A).obj Z' ≅ Z := (idxIstr P G F A).objObjPreimageIso Z
  refine ⟨HomBirat.mk Z' (InducedCategory.homMk (e.inv.unop.left.hom ≫ φ)), ?_⟩
  rw [biratIstrMap_mk]
  exact HomBirat.mk_map e.inv φ

include F G in
/-- ★★★**単射** —— 共通の下界を `idxIstr` の原像へ引き戻し、充満性で射を持ち上げる。 -/
theorem biratIstrMap_injective (A B : Istr P) :
    Function.Injective (biratIstrMap P G F A B) := by
  haveI := idxIstr_isEquivalence P G F A
  intro z₁ z₂ h
  obtain ⟨Z₁, φ₁, rfl⟩ := HomBirat.exists_rep z₁
  obtain ⟨Z₂, φ₂, rfl⟩ := HomBirat.exists_rep z₂
  rw [biratIstrMap_mk, biratIstrMap_mk] at h
  obtain ⟨V, u, v, hV⟩ := HomBirat.eq_iff.mp h
  let V' := (idxIstr P G F A).objPreimage V
  let e : (idxIstr P G F A).obj V' ≅ V := (idxIstr P G F A).objObjPreimageIso V
  obtain ⟨u', hu'⟩ := (idxIstr P G F A).map_surjective (u ≫ e.inv)
  obtain ⟨v', hv'⟩ := (idxIstr P G F A).map_surjective (v ≫ e.inv)
  refine HomBirat.sound V' u' v' (InducedCategory.hom_ext ?_)
  have hu := congrArg (fun t : (idxIstr P G F A).obj Z₁ ⟶ (idxIstr P G F A).obj V' =>
    t.unop.left.hom) hu'
  have hv := congrArg (fun t : (idxIstr P G F A).obj Z₂ ⟶ (idxIstr P G F A).obj V' =>
    t.unop.left.hom) hv'
  show u'.unop.left.hom.hom ≫ φ₁.hom = v'.unop.left.hom.hom ≫ φ₂.hom
  rw [show u'.unop.left.hom.hom = (u ≫ e.inv).unop.left.hom from hu,
    show v'.unop.left.hom.hom = (v ≫ e.inv).unop.left.hom from hv]
  show (e.inv.unop.left.hom ≫ u.unop.left.hom) ≫ φ₁.hom
    = (e.inv.unop.left.hom ≫ v.unop.left.hom) ≫ φ₂.hom
  rw [Category.assoc, Category.assoc, hV]

include F G in
/-- ★★★★★**`Hom^birat` は `𝒞^istr` と `𝒞` で一致する**。 -/
theorem biratIstrMap_bijective (A B : Istr P) :
    Function.Bijective (biratIstrMap P G F A B) :=
  ⟨biratIstrMap_injective P G F A B, biratIstrMap_surjective P G F A B⟩

def biratIstrMap_bijective.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 86,
    item := "Remark 4.5.1 — Hom^birat は 𝒞^istr と 𝒞 で一致する",
    sectionId := "frdi-remark-4-5-1" }


/-! ## ★7. 合成の保存 —— `Rmk451Birat.lean` の `biratPull_alpha_transport` と同型

★★`BiratCat` の合成は Ore 四角形(`biratPull*`)を `choose` で取るので、
**両側で選択が違いうる**。★しかし添字圏が**細い**(`idxBirat_hom_ext`)ので、
共通の上界へ送れば 2 つの選択は一致する。 -/

set_option maxHeartbeats 1000000 in
include F G in
/-- ★★★引き戻しの `α` 成分の移送 —— `γ` 成分の一致から `α` 成分の一致が出る。

★`W` の構造射が **mono**(pre-step)なので消約できる。 -/
theorem biratPullAlpha_transport {A B : Istr P}
    (Z : IdxBirat (istrPre P F) (istr_frobenioid P F G) A) (φ : Z.unop.left.obj ⟶ B)
    (W : IdxBirat (istrPre P F) (istr_frobenioid P F G) B)
    {V : IdxBirat P G A.obj}
    (c : (idxIstr P G F A).obj (biratPullIdx (istr_frobenioid P F G).core Z φ W) ⟶ V)
    (c' : biratPullIdx G.core ((idxIstr P G F A).obj Z) φ.hom
      ((idxIstr P G F B).obj W) ⟶ V)
    (hkey : c.unop.left.hom ≫ (biratPullGamma (istr_frobenioid P F G).core Z φ W).hom
        = c'.unop.left.hom ≫ biratPullGamma G.core ((idxIstr P G F A).obj Z) φ.hom
            ((idxIstr P G F B).obj W)) :
    c.unop.left.hom ≫ (biratPullAlpha (istr_frobenioid P F G).core Z φ W).hom
      = c'.unop.left.hom ≫ biratPullAlpha G.core ((idxIstr P G F A).obj Z) φ.hom
            ((idxIstr P G F B).obj W) := by
  haveI hb : Mono ((idxIstr P G F B).obj W).unop.hom.hom :=
    G.core.preStepMono _ ((idxIstr P G F B).obj W).unop.hom.property.2
  refine (cancel_mono ((idxIstr P G F B).obj W).unop.hom.hom).mp ?_
  have h1 := biratPull_sq (istr_frobenioid P F G).core Z φ W
  have h2 := biratPull_sq G.core ((idxIstr P G F A).obj Z) φ.hom ((idxIstr P G F B).obj W)
  have e1 : (biratPullAlpha (istr_frobenioid P F G).core Z φ W).hom
      ≫ ((idxIstr P G F B).obj W).unop.hom.hom
      = (biratPullGamma (istr_frobenioid P F G).core Z φ W).hom ≫ φ.hom :=
    (congrArg (fun t : biratPullObj (istr_frobenioid P F G).core Z φ W ⟶ B =>
      t.hom) h1).symm
  refine Eq.trans (Category.assoc _ _ _) ?_
  refine Eq.trans (congrArg (fun t => c.unop.left.hom ≫ t) e1) ?_
  refine Eq.trans (Category.assoc _ _ _).symm ?_
  refine Eq.trans (congrArg (fun t => t ≫ φ.hom) hkey) ?_
  refine Eq.trans (Category.assoc _ _ _) ?_
  refine Eq.trans (congrArg (fun t => c'.unop.left.hom ≫ t) h2) ?_
  exact (Category.assoc _ _ _).symm

set_option maxHeartbeats 1000000 in
include F G in
/-- ★★★代表元での合成の保存。★細さ(`idxBirat_hom_ext`)がすべてを潰す。 -/
theorem biratIstrMap_comp_mk {A B E : Istr P}
    (Z : IdxBirat (istrPre P F) (istr_frobenioid P F G) A) (φ : Z.unop.left.obj ⟶ B)
    (W : IdxBirat (istrPre P F) (istr_frobenioid P F G) B) (ψ : W.unop.left.obj ⟶ E) :
    HomBirat.mk ((idxIstr P G F A).obj
        (biratPullIdx (istr_frobenioid P F G).core Z φ W))
      ((biratPullAlpha (istr_frobenioid P F G).core Z φ W ≫ ψ).hom)
    = HomBirat.mk
        (biratPullIdx G.core ((idxIstr P G F A).obj Z) φ.hom ((idxIstr P G F B).obj W))
        (biratPullAlpha G.core ((idxIstr P G F A).obj Z) φ.hom
          ((idxIstr P G F B).obj W) ≫ ψ.hom) := by
  refine HomBirat.sound _
    (IsFiltered.leftToMax ((idxIstr P G F A).obj
      (biratPullIdx (istr_frobenioid P F G).core Z φ W))
      (biratPullIdx G.core ((idxIstr P G F A).obj Z) φ.hom ((idxIstr P G F B).obj W)))
    (IsFiltered.rightToMax ((idxIstr P G F A).obj
      (biratPullIdx (istr_frobenioid P F G).core Z φ W))
      (biratPullIdx G.core ((idxIstr P G F A).obj Z) φ.hom ((idxIstr P G F B).obj W))) ?_
  have hkey := congrArg (fun t => t.unop.left.hom)
    (idxBirat_hom_ext
      ((idxIstr P G F A).map (biratPullHom (istr_frobenioid P F G).core Z φ W)
        ≫ IsFiltered.leftToMax _ _)
      (biratPullHom G.core ((idxIstr P G F A).obj Z) φ.hom ((idxIstr P G F B).obj W)
        ≫ IsFiltered.rightToMax _ _))
  have hmid := biratPullAlpha_transport P G F Z φ W _ _ hkey
  simp only [← Category.assoc]
  exact Eq.trans (Category.assoc _ _ _).symm
    (congrArg (fun t => t ≫ ψ.hom) hmid)

set_option maxHeartbeats 1000000 in
include F G in
/-- ★★★★`biratIstrMap` は合成を保つ。 -/
theorem biratIstrMap_comp {A B E : Istr P}
    (f : HomBirat (istrPre P F) (istr_frobenioid P F G) A B)
    (g : HomBirat (istrPre P F) (istr_frobenioid P F G) B E) :
    biratIstrMap P G F A E
        (compBirat (istrPre P F) (istr_frobenioid P F G)
          (istr_frobenioid P F G).core f g)
      = compBirat P G G.core (biratIstrMap P G F A B f) (biratIstrMap P G F B E g) := by
  obtain ⟨Z, φ, rfl⟩ := HomBirat.exists_rep f
  obtain ⟨W, ψ, rfl⟩ := HomBirat.exists_rep g
  rw [compBirat_mk, biratIstrMap_mk, biratIstrMap_mk, biratIstrMap_mk, compBirat_mk]
  exact biratIstrMap_comp_mk P G F Z φ W ψ

include F G in
/-- ★`biratIstrMap` は恒等を保つ。 -/
theorem biratIstrMap_id (A : Istr P) :
    biratIstrMap P G F A A (toHomBirat (𝟙 A)) = toHomBirat (𝟙 A.obj) := by
  show biratIstrMap P G F A A
      (HomBirat.mk (idxBiratOne (istrPre P F) (istr_frobenioid P F G) A) (𝟙 A)) = _
  rw [biratIstrMap_mk]
  rfl

def biratIstrMap_comp.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 86,
    item := "Remark 4.5.1 — 双有理化の比較は合成を保つ",
    sectionId := "frdi-remark-4-5-1" }

/-! ## ★8. 比較関手 —— `(𝒞^istr)^birat ⥤ 𝒞^birat`

★★これは**充満忠実**であり(`biratIstrMap_bijective`)、
`Base`・`degFr` を保つ。★対象は isotropic なものに限られるので圏同値ではないが、
**`End` の単系同型**が立つので ∀ 条件も ∃ 条件も両向きに降りる。 -/

/-- ★★★★**比較関手** —— `(𝒞^istr)^birat ⥤ 𝒞^birat`。 -/
noncomputable def biratIstrFunctor :
    BiratCat (istrPre P F) (istr_frobenioid P F G) ⥤ BiratCat P G where
  obj X := (show Istr P from X).obj
  map {A B} f := biratIstrMap P G F A B f
  map_id A := biratIstrMap_id P G F A
  map_comp f g := biratIstrMap_comp P G F f g

instance biratIstrFunctor_full : (biratIstrFunctor P G F).Full where
  map_surjective {A B} f := (biratIstrMap_bijective P G F A B).2 f

instance biratIstrFunctor_faithful : (biratIstrFunctor P G F).Faithful where
  map_injective {A B} {_ _} h := (biratIstrMap_bijective P G F A B).1 h

include F G in
/-- ★★`biratIstrMap` は Frobenius 次数を保つ。 -/
theorem biratIstrMap_degFr (A B : Istr P)
    (z : HomBirat (istrPre P F) (istr_frobenioid P F G) A B) :
    biratDeg (biratIstrMap P G F A B z) = biratDeg z := by
  obtain ⟨Z, φ, rfl⟩ := HomBirat.exists_rep z
  rw [biratIstrMap_mk, biratDeg_mk, biratDeg_mk]
  exact (istr_compat_degFr P F φ).symm

include F G in
/-- ★★`biratIstrMap` は底を保つ。 -/
theorem biratIstrMap_base (A B : Istr P)
    (z : HomBirat (istrPre P F) (istr_frobenioid P F G) A B) :
    biratBase (biratIstrMap P G F A B z) = biratBase z := by
  obtain ⟨Z, φ, rfl⟩ := HomBirat.exists_rep z
  rw [biratIstrMap_mk, biratBase_mk, biratBase_mk, sliceBaseOf_eq, sliceBaseOf_eq]
  rfl

def biratIstrFunctor.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 86,
    item := "Remark 4.5.1 — (𝒞^istr)^birat ⥤ 𝒞^birat は充満忠実で構造を保つ",
    sectionId := "frdi-remark-4-5-1" }

/-! ## ★9. `Definition 4.5, (iii), (a)` の第 1 条 —— birationally Frobenius-normalized 型

原文 (FrdI p.86):
> that if C is of rationally standard type (respectively, of standard type), then so is

★★`IsFrobeniusNormalized` は **∀ 条件**なので、
充満忠実で `Base`・`degFr` を保つ関手があれば**そのまま降りる**。 -/

include F G in
/-- ★`biratIstrFunctor` は `(biratPre _).Base` を保つ。 -/
theorem biratIstrFunctor_Base {X Y : BiratCat (istrPre P F) (istr_frobenioid P F G)}
    (f : X ⟶ Y) :
    (biratPre P G).Base ((biratIstrFunctor P G F).map f)
      = (biratPre (istrPre P F) (istr_frobenioid P F G)).Base f :=
  biratIstrMap_base P G F X Y f

include F G in
/-- ★`biratIstrFunctor` は `(biratPre _).degFr` を保つ。 -/
theorem biratIstrFunctor_degFr {X Y : BiratCat (istrPre P F) (istr_frobenioid P F G)}
    (f : X ⟶ Y) :
    (biratPre P G).degFr ((biratIstrFunctor P G F).map f)
      = (biratPre (istrPre P F) (istr_frobenioid P F G)).degFr f :=
  biratIstrMap_degFr P G F X Y f

include F G in
/-- ★★`𝒪^▷` は `biratIstrFunctor` で移る。 -/
theorem biratIstrFunctor_mem_otri {X : BiratCat (istrPre P F) (istr_frobenioid P F G)}
    {α : End X} (h : α ∈ OTri (biratPre (istrPre P F) (istr_frobenioid P F G)) X) :
    (CategoryTheory.Functor.mapEnd X (biratIstrFunctor P G F)) α
      ∈ OTri (biratPre P G) ((biratIstrFunctor P G F).obj X) := by
  refine ⟨?_, ?_⟩
  · show (biratPre P G).Base ((biratIstrFunctor P G F).map (α : X ⟶ X))
      = (biratPre P G).Base (𝟙 _)
    rw [biratIstrFunctor_Base, (biratPre P G).Base_id, h.1,
      (biratPre (istrPre P F) (istr_frobenioid P F G)).Base_id]
    rfl
  · show (biratPre P G).degFr ((biratIstrFunctor P G F).map (α : X ⟶ X)) = 1
    rw [biratIstrFunctor_degFr]
    exact h.2

include F G in
/-- ★★★★★**(a) の第 1 条** —— birationally Frobenius-normalized 型は `𝒞^istr` へ移る。

★★`IsFrobeniusNormalized` は ∀ 条件なので、
**充満忠実性 ＋ `Base`・`degFr` の保存**だけで降りる。 -/
theorem istr_biratFrobNormalizedType
    (h : IsOfBirationallyFrobeniusNormalizedType C P G) :
    IsOfBirationallyFrobeniusNormalizedType (Istr P) (istrPre P F)
      (istr_frobenioid P F G) := by
  intro A
  set X : BiratCat (istrPre P F) (istr_frobenioid P F G) :=
    (toBiratCat (istrPre P F) (istr_frobenioid P F G)).obj A with hX
  intro φ hφ α hα
  have hΘ := h A.obj
  have hbi : IsBaseIdentity (biratPre P G)
      ((CategoryTheory.Functor.mapEnd X (biratIstrFunctor P G F)) φ) := by
    show (biratPre P G).Base ((biratIstrFunctor P G F).map (φ : X ⟶ X))
      = (biratPre P G).Base (𝟙 _)
    rw [biratIstrFunctor_Base, (biratPre P G).Base_id, hφ,
      (biratPre (istrPre P F) (istr_frobenioid P F G)).Base_id]
    rfl
  have hmain := hΘ _ hbi _ (biratIstrFunctor_mem_otri P G F hα)
  have hd : (biratPre P G).degFr
      ((CategoryTheory.Functor.mapEnd X (biratIstrFunctor P G F)) φ)
      = (biratPre (istrPre P F) (istr_frobenioid P F G)).degFr φ :=
    biratIstrFunctor_degFr P G F _
  have hpow := congrArg (fun n : ℕ+ =>
    (((CategoryTheory.Functor.mapEnd X (biratIstrFunctor P G F)) α ^ ((n : ℕ+) : ℕ)
      : End ((biratIstrFunctor P G F).obj X)) :
        (biratIstrFunctor P G F).obj X ⟶ (biratIstrFunctor P G F).obj X)) hd
  have hmain2 := (congrArg (fun t =>
    (((CategoryTheory.Functor.mapEnd X (biratIstrFunctor P G F)) φ) :
      (biratIstrFunctor P G F).obj X ⟶ (biratIstrFunctor P G F).obj X) ≫ t)
      hpow).symm.trans hmain
  refine (biratIstrFunctor P G F).map_injective ?_
  rw [CategoryTheory.Functor.map_comp, CategoryTheory.Functor.map_comp]
  refine Eq.trans ?_ hmain2
  exact congrArg (fun t =>
    (((CategoryTheory.Functor.mapEnd X (biratIstrFunctor P G F)) φ) :
      (biratIstrFunctor P G F).obj X ⟶ (biratIstrFunctor P G F).obj X) ≫ t)
    (map_pow (CategoryTheory.Functor.mapEnd X (biratIstrFunctor P G F)) α _)

def istr_biratFrobNormalizedType.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 86,
    item := "Remark 4.5.1 — birationally Frobenius-normalized 型は 𝒞^istr で保たれる",
    sectionId := "frdi-remark-4-5-1" }

/-! ## ★10. isotropic hull に沿った自己射の移送 —— 一般の pre-Frobenioid で

★★`Definition 1.2, (iv)` の isotropic hull は `∃!-lift` を持つので、
**`𝒪^▷` も `𝒪^×` もそのまま向こう側へ移る**。
★これは Frobenioid 一般の話であり、`𝒞^birat` に当てるために先に一般形で取る。 -/

section HullTransport

universe v' u' w' u2' v2'

variable {D' : Type u'} [Category.{v'} D'] {C' : Type u2'} [Category.{v2'} C']
  {Φ' : MonoidOn.{v', u', w'} D'} (Q : PreFrobenioid C' Φ')

/-- ★★★**isotropic hull に沿って自己射が一意に移る**。 -/
theorem hull_transport_end {B H : C'} {k : B ⟶ H} (hk : IsIsotropicHull Q k) (δ : End B) :
    ∃! δ' : End H, ((δ : B ⟶ B) ≫ k) = k ≫ (δ' : H ⟶ H) := by
  obtain ⟨-, -, hHiso, hlift⟩ := hk
  obtain ⟨β, hβ, huniq⟩ := hlift H hHiso ((δ : B ⟶ B) ≫ k)
  exact ⟨β, hβ, fun γ hγ => huniq γ hγ⟩

/-- ★★移った自己射は `𝒪^▷` に入る。 -/
theorem hull_transport_mem_otri {B H : C'} {k : B ⟶ H} (hk : IsIsotropicHull Q k)
    {δ : End B} (hδ : δ ∈ OTri Q B) {δ' : End H}
    (h : ((δ : B ⟶ B) ≫ k) = k ≫ (δ' : H ⟶ H)) : δ' ∈ OTri Q H := by
  haveI hki : IsIso (Q.Base k) := hk.2.1.2
  refine ⟨?_, ?_⟩
  · have hb := congrArg Q.Base h
    rw [Q.Base_comp, Q.Base_comp,
      show Q.Base (δ : B ⟶ B) = Q.Base (𝟙 B) from hδ.1, Q.Base_id, Category.id_comp] at hb
    show Q.Base (δ' : H ⟶ H) = Q.Base (𝟙 H)
    rw [Q.Base_id]
    refine (cancel_epi (Q.Base k)).mp ?_
    rw [Category.comp_id]
    exact hb.symm
  · have hd := congrArg Q.degFr h
    rw [Q.degFr_comp, Q.degFr_comp,
      show Q.degFr (δ : B ⟶ B) = 1 from hδ.2, mul_one] at hd
    show Q.degFr (δ' : H ⟶ H) = 1
    exact (mul_right_cancel (a := Q.degFr (δ' : H ⟶ H)) (b := Q.degFr k)
      (c := 1) (by rw [one_mul]; exact hd.symm))

/-- ★★★移した自己射は**単元**になる —— `∃!` の一意性から。 -/
theorem hull_transport_isUnit {B H : C'} {k : B ⟶ H} (hk : IsIsotropicHull Q k)
    {δ ε : End B} {δ' ε' : End H}
    (hde : ((δ : B ⟶ B) ≫ (ε : B ⟶ B)) = 𝟙 B)
    (hed : ((ε : B ⟶ B) ≫ (δ : B ⟶ B)) = 𝟙 B)
    (h : ((δ : B ⟶ B) ≫ k) = k ≫ (δ' : H ⟶ H))
    (h' : ((ε : B ⟶ B) ≫ k) = k ≫ (ε' : H ⟶ H)) :
    IsUnit δ' := by
  obtain ⟨-, -, hHiso, hlift⟩ := hk
  obtain ⟨β0, -, huniq⟩ := hlift H hHiso k
  have hid : (𝟙 H : H ⟶ H) = β0 := huniq _ (Category.comp_id k).symm
  have e2 : k ≫ ((ε' : H ⟶ H) ≫ (δ' : H ⟶ H)) = k := by
    rw [← Category.assoc, ← h', Category.assoc, ← h, ← Category.assoc, hed, Category.id_comp]
  have e1 : k ≫ ((δ' : H ⟶ H) ≫ (ε' : H ⟶ H)) = k := by
    rw [← Category.assoc, ← h, Category.assoc, ← h', ← Category.assoc, hde, Category.id_comp]
  have hone2 : ((ε' : H ⟶ H) ≫ (δ' : H ⟶ H)) = 𝟙 H := by
    rw [huniq _ e2.symm, ← hid]
  have hone1 : ((δ' : H ⟶ H) ≫ (ε' : H ⟶ H)) = 𝟙 H := by
    rw [huniq _ e1.symm, ← hid]
  exact ⟨⟨δ', ε', hone2, hone1⟩, rfl⟩

/-- ★★★★**移送は `Div` を底の同型で運ぶ** —— hull は等長なので余計な項が出ない。

★`Div(δ) = Φ(Base k)(Div(δ'))`。★これが `Φ^birat` の移送の実体である。 -/
theorem hull_transport_Div {B H : C'} {k : B ⟶ H} (hk : IsIsotropicHull Q k)
    {δ : End B} (hδ : δ ∈ OTri Q B) {δ' : End H} (hδ' : δ' ∈ OTri Q H)
    (h : ((δ : B ⟶ B) ≫ k) = k ≫ (δ' : H ⟶ H)) :
    Q.Div (δ : B ⟶ B) = Φ'.map (Q.Base k) (Q.Div (δ' : H ⟶ H)) := by
  have hdv := congrArg Q.Div h
  rw [Q.Div_comp, Q.Div_comp] at hdv
  rw [show Q.Base (δ : B ⟶ B) = Q.Base (𝟙 B) from hδ.1, Q.Base_id, Φ'.map_id] at hdv
  rw [show Q.degFr k = 1 from hk.2.1.1, show Q.degFr (δ' : H ⟶ H) = 1 from hδ'.2] at hdv
  rw [show Q.Div k = 0 from hk.1] at hdv
  simpa using hdv

end HullTransport

def hull_transport_end.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 22,
    item := "Definition 1.2, (iv) — isotropic hull に沿った自己射の移送",
    sectionId := "frdi-def-1-2-iv" }

/-! ## ★11. `Definition 4.5, (iii), (a)` の第 2 条 —— rational 型の引き戻しの側

原文 (FrdI p.86):
> Definition 2.4, (i), (d)]. We shall say that A is rational if there exists a pull-back

★★`IsRational A := ∃ B (φ : B ⟶ A), IsPullBack φ ∧ IsStrictlyRational B` であり、
**`B` が isotropic とは限らない**のが `𝒞^istr` へ移すときの障害だった。
★★★等長化 `B^istr` を取れば引き戻しの側は**在庫だけで立つ**。 -/

include F G in
/-- ★★★★**引き戻しの側** —— `B^istr` から `A` への pull-back が立つ。

★`isotropification_pullBack`(`Prop19.lean`)＋ `hullMap_isIso`(`A` が isotropic)
＋ `IsPullBack.comp` の 3 本。 -/
theorem istr_rational_pullBack {A : Istr P} {B : C} (φ : B ⟶ A.obj)
    (hpb : IsPullBack P φ) :
    ∃ ψ : hullIstr P F B ⟶ A, IsPullBack (istrPre P F) ψ := by
  haveI : IsIso (hullMap P F A.obj) := hullMap_isIso P F A.obj A.property
  have h1 : IsPullBack P (istrMap P F φ) := isotropification_pullBack P F φ hpb
  have h2 : IsPullBack P (inv (hullMap P F A.obj)) :=
    isPullBack_of_isIso P (inv (hullMap P F A.obj))
  refine ⟨InducedCategory.homMk (istrMap P F φ ≫ inv (hullMap P F A.obj)), ?_⟩
  exact istr_isPullBack_of P F _ (IsPullBack.comp P h1 h2)

def istr_rational_pullBack.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 86,
    item := "Remark 4.5.1 — rational 型の引き戻しの側",
    sectionId := "frdi-remark-4-5-1" }

/-! ## ★12. isotropic hull の一意性と同型による移送 —— 一般の pre-Frobenioid で

★★`𝒞^birat` の hull(`birat_isotropicHullExists`)は **`G.core` が選んだ** hull の像であり、
`F` が選んだ `hullMap P F B` とは選択が違いうる。
★★★しかし **isotropic hull は同型を除いて一意**なので、その差は同型に吸収される。 -/

section HullUniq

universe v'' u'' w'' u2'' v2''

variable {D'' : Type u''} [Category.{v''} D''] {C'' : Type u2''} [Category.{v2''} C'']
  {Φ'' : MonoidOn.{v'', u'', w''} D''} (Q : PreFrobenioid C'' Φ'')

/-! ### ★在庫にあった(見落とし 12 件目)

`isIsotropicHull_comp_iso`(`Prop19.lean:238`)——
**isotropic hull に同型を後置しても isotropic hull**。
★`isotropification_*` の名前でしか検索しておらず見落とした。
★★`Proposition 1.9` は等長化そのものなので、
**hull まわりはまず `Prop19.lean` を通しで読むべき**である。 -/

/-- ★★★**isotropic hull は同型を除いて一意** —— 2 つの hull の間に構造射を保つ同型がある。 -/
theorem isIsotropicHull_unique {A B B' : C''} {φ : A ⟶ B} {φ' : A ⟶ B'}
    (h : IsIsotropicHull Q φ) (h' : IsIsotropicHull Q φ') :
    ∃ θ : B ≅ B', φ ≫ θ.hom = φ' := by
  obtain ⟨-, -, hiso, hlift⟩ := h
  obtain ⟨-, -, hiso', hlift'⟩ := h'
  obtain ⟨u, hu, huq⟩ := hlift B' hiso' φ'
  obtain ⟨v, hv, hvq⟩ := hlift' B hiso φ
  have h1 : u ≫ v = 𝟙 B := by
    refine (hlift B hiso φ).unique ?_ (Category.comp_id φ).symm
    rw [← Category.assoc, ← hu]; exact hv
  have h2 : v ≫ u = 𝟙 B' := by
    refine (hlift' B' hiso' φ').unique ?_ (Category.comp_id φ').symm
    rw [← Category.assoc, ← hv]; exact hu
  exact ⟨⟨u, v, h1, h2⟩, hu.symm⟩

end HullUniq

def isIsotropicHull_unique.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 22,
    item := "Definition 1.2, (iv) — isotropic hull は同型を除いて一意",
    sectionId := "frdi-def-1-2-iv" }

/-! ## ★13. 等長化 hull は `𝒞^birat` でも isotropic hull

原文 (FrdI p.83):
> An object of C maps to an isotropic object of Cbirat if and only if

★★在庫の `birat_isotropicHullExists` は **`G.core` が選んだ** hull を与えるが、
引き戻しに使うのは **`F` が選んだ `hullMap P F B`** である。
★★★両者の差は `isIsotropicHull_unique` ではなく、
**`𝒞^birat` で pre-step が同型になる**こと(`birat_isIso_of_preStep_of_isotropic`)で潰れる。 -/

include F G in
/-- ★★★★★**`hullMap` の birat 化は `𝒞^birat` の isotropic hull**。

★`G.core` の hull `k₀` から降ろした `β` が
「isotropic 対象から出る pre-step」なので**同型**になり、
`isIsotropicHull_comp_iso`(在庫)で閉じる。 -/
theorem birat_map_hullMap_isIsotropicHull
    (hfn : ∀ X : BiratCat P G, IsFrobeniusNormalized (biratPre P G) X) (B : C) :
    IsIsotropicHull (biratPre P G) ((toBiratCat P G).map (hullMap P F B)) := by
  obtain ⟨H₀, k₀, hk₀⟩ := birat_isotropicHullExists (P := P) (G := G) ((toBiratCat P G).obj B)
  obtain ⟨β, hβ, -⟩ := hk₀.2.2.2 ((toBiratCat P G).obj (hullObj P F B))
    (birat_isotropic_up (P := P) (G := G) (hullMap_spec P F B).2.2.1)
    ((toBiratCat P G).map (hullMap P F B))
  have hd0 : (biratPre P G).degFr ((toBiratCat P G).map (hullMap P F B)) = 1 :=
    (birat_isLinear_iff (P := P) (G := G) (hullMap P F B)).mpr (hullMap_spec P F B).2.1.1
  have hb0 : IsIso ((biratPre P G).Base ((toBiratCat P G).map (hullMap P F B))) :=
    (birat_isBaseIsomorphism_iff (P := P) (G := G) (hullMap P F B)).mpr (hullMap_spec P F B).2.1.2
  have hdβ : (biratPre P G).degFr β = 1 := by
    have h1 := congrArg (biratPre P G).degFr hβ
    rw [hd0, (biratPre P G).degFr_comp,
      show (biratPre P G).degFr k₀ = 1 from hk₀.2.1.1, mul_one] at h1
    exact h1.symm
  have hbβ : IsBaseIsomorphism (biratPre P G) β := by
    have h2 := congrArg (biratPre P G).Base hβ
    rw [(biratPre P G).Base_comp] at h2
    haveI : IsIso ((biratPre P G).Base k₀) := hk₀.2.1.2
    haveI : IsIso ((biratPre P G).Base k₀ ≫ (biratPre P G).Base β) := by
      rw [← h2]; exact hb0
    exact IsIso.of_isIso_comp_left ((biratPre P G).Base k₀) ((biratPre P G).Base β)
  haveI : IsIso β :=
    birat_isIso_of_preStep_of_isotropic P G hfn hk₀.2.2.1 ⟨hdβ, hbβ⟩
  have hfin := isIsotropicHull_comp_iso (biratPre P G) k₀ β hk₀
  rwa [← hβ] at hfin

def birat_map_hullMap_isIsotropicHull.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 83,
    item := "Proposition 4.4, (iv) — 等長化 hull は 𝒞^birat でも isotropic hull",
    sectionId := "frdi-prop-4-4" }

/-! ## ★14. `Div^gp` の hull に沿った移送

★★★**測ったこと**: `biratPre P G` の単系は `trivialOn D`(自明)なので、
`hull_transport_Div` は `𝒞^birat` では**何も言わない**。
★`phiBiratAt` が使うのは `biratDivGp`(`Gp (Φ.val ...)` へ行く別の写像)であり、
そちらの合成則(`biratDivGp_comp'`)で書き直す必要がある。 -/

variable {P G} in
/-- ★★★★★**`Div^gp` は hull で底の同型に沿って移る**。

★`biratDivGp (f ≫ g) = Φ^gp(Base f)(Div^gp g) + deg(g) • Div^gp f` を
`δ ≫ k = k ≫ δ'` の両辺に当て、`Div^gp k` を消約する。 -/
theorem biratDivGp_hull_transport {B H : BiratCat P G} (k : B ⟶ H)
    (hkd : biratDeg k = 1)
    {δ : End B} (hδ : δ ∈ OTri (biratPre P G) B)
    {δ' : End H} (hδ' : δ' ∈ OTri (biratPre P G) H)
    (h : ((δ : B ⟶ B) ≫ k) = k ≫ (δ' : H ⟶ H)) :
    biratDivGp ((δ : B ⟶ B)) = gpMap _ (Φ.map (biratBase k)) (biratDivGp ((δ' : H ⟶ H))) := by
  have hbδ : biratBase ((δ : B ⟶ B)) = 𝟙 _ := by
    have h1 : (biratPre P G).Base ((δ : B ⟶ B)) = (biratPre P G).Base (𝟙 B) := hδ.1
    rw [(biratPre P G).Base_id] at h1
    exact h1
  have hdδ' : biratDeg ((δ' : H ⟶ H)) = 1 := hδ'.2
  have hcomp := congrArg (biratDivGp (P := P) (G := G)) h
  rw [biratDivGp_comp', biratDivGp_comp', hbδ, hkd, hdδ'] at hcomp
  have hmapid : (gpMap (Φ.val (P.toElem.obj (biratDown P G B)).base)
      (Φ.map (𝟙 (P.toElem.obj (biratDown P G B)).base))) (biratDivGp k) = biratDivGp k := by
    have hmid : Φ.map (𝟙 (P.toElem.obj (biratDown P G B)).base)
        = AddMonoidHom.id (Φ.val (P.toElem.obj (biratDown P G B)).base) :=
      AddMonoidHom.ext fun z => Φ.map_id _ z
    rw [hmid, gpMap_id]
    rfl
  simp only [one_smul, PNat.one_coe] at hcomp
  rw [hmapid] at hcomp
  -- `Div^gp k + Div^gp δ = Φ^gp(Base k)(Div^gp δ') + Div^gp k`
  rw [add_comm (biratDivGp k)] at hcomp
  exact add_right_cancel hcomp

/-! ## ★15. 比較関手は `Div^gp` を保つ —— `Φ^birat` が両側で一致する

★★`sliceDivGpOf` は `P.Div` と `P.Base` と `P.degFr` だけで書かれており、
`istr_compat_Div` / `istr_compat_Base` / `istr_compat_degFr` がすべて `rfl` なので、
**そのまま一致する**。 -/

include F G in
/-- ★★★★`biratIstrMap` は `Div^gp` を保つ。 -/
theorem biratIstrMap_divGp (A B : Istr P)
    (z : HomBirat (istrPre P F) (istr_frobenioid P F G) A B) :
    biratDivGp (biratIstrMap P G F A B z) = biratDivGp z := by
  obtain ⟨Z, φ, rfl⟩ := HomBirat.exists_rep z
  rw [biratIstrMap_mk, biratDivGp_mk, biratDivGp_mk, sliceDivGpOf_eq, sliceDivGpOf_eq]
  rfl

def biratIstrMap_divGp.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 86,
    item := "Remark 4.5.1 — Φ^birat は 𝒞^istr と 𝒞 で一致する",
    sectionId := "frdi-remark-4-5-1" }

/-- ★全単射な単系準同型は単元性を反射する。 -/
theorem isUnit_of_bijective_map {M₀ N₀ : Type*} [Monoid M₀] [Monoid N₀] (f : M₀ →* N₀)
    (hf : Function.Bijective f) {a : M₀} (h : IsUnit (f a)) : IsUnit a := by
  obtain ⟨u, hu⟩ := h
  obtain ⟨b, hb⟩ := hf.2 (↑u⁻¹ : N₀)
  refine ⟨⟨a, b, ?_, ?_⟩, rfl⟩
  · exact hf.1 (by rw [map_mul, hb, ← hu, u.mul_inv, map_one])
  · exact hf.1 (by rw [map_mul, hb, ← hu, u.inv_mul, map_one])

include F G in
/-- ★★★★★**`Φ^birat` の元は `𝒞^istr` 側へ持ち上がる** ——
比較関手が充満忠実で `Base`・`degFr`・`Div^gp` を保つから。 -/
theorem mem_phiBiratAt_istr_of (A : Istr P)
    {x : Gp (Φ.val (P.toElem.obj A.obj).base)}
    (hx : x ∈ phiBiratAt P G ((toBiratCat P G).obj A.obj)) :
    x ∈ phiBiratAt (istrPre P F) (istr_frobenioid P F G)
        ((toBiratCat (istrPre P F) (istr_frobenioid P F G)).obj A) := by
  set X : BiratCat (istrPre P F) (istr_frobenioid P F G) :=
    (toBiratCat (istrPre P F) (istr_frobenioid P F G)).obj A with hX
  obtain ⟨ε, hε, hεx⟩ := hx
  obtain ⟨δ, hδeq⟩ := (biratIstrFunctor P G F).map_surjective (ε : _ ⟶ _)
  have hmap : (CategoryTheory.Functor.mapEnd X (biratIstrFunctor P G F)) δ = ε := hδeq
  refine ⟨δ, ⟨⟨?_, ?_⟩, ?_⟩, ?_⟩
  · show (biratPre (istrPre P F) (istr_frobenioid P F G)).Base (δ : X ⟶ X)
      = (biratPre (istrPre P F) (istr_frobenioid P F G)).Base (𝟙 X)
    rw [(biratPre (istrPre P F) (istr_frobenioid P F G)).Base_id,
      ← biratIstrFunctor_Base P G F (δ : X ⟶ X), hδeq]
    have h1 : (biratPre P G).Base (ε : _ ⟶ _) = (biratPre P G).Base (𝟙 _) := hε.1.1
    rw [(biratPre P G).Base_id] at h1
    exact h1
  · show (biratPre (istrPre P F) (istr_frobenioid P F G)).degFr (δ : X ⟶ X) = 1
    rw [← biratIstrFunctor_degFr P G F (δ : X ⟶ X), hδeq]
    exact hε.1.2
  · refine isUnit_of_bijective_map (CategoryTheory.Functor.mapEnd X (biratIstrFunctor P G F))
      (biratIstrMap_bijective P G F A A) ?_
    rw [hmap]
    exact hε.2
  · have h2 : biratDivGp (biratIstrMap P G F A A (δ : X ⟶ X)) = biratDivGp (δ : X ⟶ X) :=
      biratIstrMap_divGp P G F A A (δ : X ⟶ X)
    refine h2.symm.trans ?_
    rw [show biratIstrMap P G F A A (δ : X ⟶ X) = (ε : _ ⟶ _) from hδeq]
    exact hεx

include F G in
/-- ★★★★★**`Φ^birat(B)` の元は `B^istr` 側から底の同型で来る**。

★`hull_transport_end` で自己射を移し、`hull_transport_mem_otri` /
`hull_transport_isUnit` で `𝒪^×` に入ることを見て、
`biratDivGp_hull_transport` で `Div^gp` を運ぶ。 -/
theorem phiBiratAt_hull
    (hfn : ∀ X : BiratCat P G, IsFrobeniusNormalized (biratPre P G) X) (B : C)
    {x : Gp (Φ.val (P.toElem.obj B).base)}
    (hx : x ∈ phiBiratAt P G ((toBiratCat P G).obj B)) :
    ∃ y ∈ phiBiratAt P G ((toBiratCat P G).obj (hullObj P F B)),
      x = gpMap _ (Φ.map (biratBase ((toBiratCat P G).map (hullMap P F B)))) y := by
  obtain ⟨δ, hδ, rfl⟩ := hx
  have hk := birat_map_hullMap_isIsotropicHull P G F hfn B
  obtain ⟨δ', hδ'k, -⟩ := hull_transport_end (biratPre P G) hk δ
  have hotri := hull_transport_mem_otri (biratPre P G) hk hδ.1 hδ'k
  obtain ⟨u, hu⟩ := hδ.2
  obtain ⟨ε', hε'k, -⟩ := hull_transport_end (biratPre P G) hk ((↑u⁻¹ : End _))
  have hde : ((δ : _ ⟶ _) ≫ ((↑u⁻¹ : End _) : _ ⟶ _)) = 𝟙 _ := by
    have h1 := u.inv_val
    rw [hu] at h1
    exact h1
  have hed : (((↑u⁻¹ : End _) : _ ⟶ _) ≫ (δ : _ ⟶ _)) = 𝟙 _ := by
    have h1 := u.val_inv
    rw [hu] at h1
    exact h1
  have hunit := hull_transport_isUnit (biratPre P G) hk hde hed hδ'k hε'k
  refine ⟨biratDivGp ((δ' : _ ⟶ _)), ⟨δ', ⟨hotri, hunit⟩, rfl⟩, ?_⟩
  exact biratDivGp_hull_transport _ hk.2.1.1 hδ.1 hotri hδ'k

def phiBiratAt_hull.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 86,
    item := "Remark 4.5.1 — Φ^birat は等長化 hull を通って移る",
    sectionId := "frdi-remark-4-5-1" }

/-! ## ★16. 組み上げ —— `Definition 4.5, (iii), (a)` の第 2 条

原文 (FrdI p.86):
> Definition 2.4, (i), (d)]. We shall say that A is rational if there exists a pull-back
-/

/-- ★単系同型から `Gp` の全単射。 -/
theorem gpMap_bijective_of_addEquiv {M₁ N₁ : Type w} [AddCommMonoid M₁] [AddCommMonoid N₁]
    (e : M₁ ≃+ N₁) : Function.Bijective (gpMap M₁ e.toAddMonoidHom) := by
  have h1 : (e.symm.toAddMonoidHom).comp (e.toAddMonoidHom) = AddMonoidHom.id M₁ :=
    AddMonoidHom.ext fun z => e.symm_apply_apply z
  have h2 : (e.toAddMonoidHom).comp (e.symm.toAddMonoidHom) = AddMonoidHom.id N₁ :=
    AddMonoidHom.ext fun z => e.apply_symm_apply z
  have key : (gpMap N₁ e.symm.toAddMonoidHom).comp (gpMap M₁ e.toAddMonoidHom)
      = AddMonoidHom.id (Gp M₁) := by rw [← gpMap_comp, h1, gpMap_id]
  have key2 : (gpMap M₁ e.toAddMonoidHom).comp (gpMap N₁ e.symm.toAddMonoidHom)
      = AddMonoidHom.id (Gp N₁) := by rw [← gpMap_comp, h2, gpMap_id]
  exact Function.bijective_iff_has_inverse.mpr ⟨gpMap N₁ e.symm.toAddMonoidHom,
    fun z => congrArg (fun f : Gp M₁ →+ Gp M₁ => f z) key,
    fun z => congrArg (fun f : Gp N₁ →+ Gp N₁ => f z) key2⟩

include F G in
/-- ★★★★★**(a) の第 2 条** —— rational 型は `𝒞^istr` へ移る。

★引き戻しは `istr_rational_pullBack`、`Φ^birat` は `phiBiratAt_hull` ＋
`mem_phiBiratAt_istr_of`、`Supp` は `mem_supp_map` で運ぶ。 -/
theorem istr_rationalType
    {ι : ∀ Y : D, Prime (Φ.val Y) → Pf (Φ.val Y) → ℝ≥0}
    (Hperf : ∀ Y : D, IsPerfFactorialWith (Φ.val Y) (ι Y))
    (hfn : ∀ X : BiratCat P G, IsFrobeniusNormalized (biratPre P G) X)
    (h : IsOfRationalType C P G ι) :
    IsOfRationalType (Istr P) (istrPre P F) (istr_frobenioid P F G) ι := by
  intro A
  obtain ⟨B, φ, hpb, hsr⟩ := h A.obj
  obtain ⟨ψ, hψ⟩ := istr_rational_pullBack P G F φ hpb
  refine ⟨hullIstr P F B, ψ, hψ, ?_⟩
  haveI hbi : IsIso (P.Base (hullMap P F B)) := (hullMap_spec P F B).2.1.2
  set eiso : (P.toElem.obj B).base ≅ (P.toElem.obj (hullObj P F B)).base :=
    asIso (P.Base (hullMap P F B)) with heiso
  set e : Φ.val (P.toElem.obj (hullObj P F B)).base ≃+ Φ.val (P.toElem.obj B).base :=
    AddEquiv.ofBijective (Φ.map eiso.hom) (Φ.map_bijective_of_iso eiso) with he
  intro p
  obtain ⟨a₁, b₁, hmem, hpa, hpb'⟩ := hsr (primeEquiv e p)
  refine ⟨e.symm a₁, e.symm b₁, ?_, ?_, ?_⟩
  · -- `Φ^birat` の元であること
    obtain ⟨y, hy, hxy⟩ := phiBiratAt_hull P G F hfn B hmem
    have hgb : biratBase ((toBiratCat P G).map (hullMap P F B)) = P.Base (hullMap P F B) :=
      biratBase_toHomBirat (P := P) (G := G) (hullMap P F B)
    rw [hgb] at hxy
    have hey : gpMap _ e.toAddMonoidHom y = toGp _ a₁ - toGp _ b₁ := hxy.symm
    have hinj := (gpMap_bijective_of_addEquiv e).1
    have hgoal : y = toGp _ (e.symm a₁) - toGp _ (e.symm b₁) := by
      refine hinj ?_
      rw [hey, map_sub, gpMap_toGp, gpMap_toGp]
      show toGp _ a₁ - toGp _ b₁ = toGp _ (e (e.symm a₁)) - toGp _ (e (e.symm b₁))
      rw [e.apply_symm_apply, e.apply_symm_apply]
    have hres := mem_phiBiratAt_istr_of P G F (hullIstr P F B) hy
    rw [hgoal] at hres
    exact hres
  · -- `Supp` の条件(正)
    refine (mem_supp_map (Hperf _) (Hperf _) e (Pf.mk (e.symm a₁) 1) p).mpr ?_
    show primeEquiv e p ∈ Supp (factorMap (ι _) (Pf.mk (e (e.symm a₁)) 1))
    rw [e.apply_symm_apply]
    exact hpa
  · -- `Supp` の条件(負)
    intro hcon
    refine hpb' ?_
    have := (mem_supp_map (Hperf _) (Hperf _) e (Pf.mk (e.symm b₁) 1) p).mp hcon
    show primeEquiv e p ∈ Supp (factorMap (ι _) (Pf.mk b₁ 1))
    rw [← e.apply_symm_apply b₁]
    exact this

def istr_rationalType.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 86,
    item := "Remark 4.5.1 — rational 型は 𝒞^istr で保たれる",
    sectionId := "frdi-remark-4-5-1" }

/-! ## ★17. ★★★★★★`Remark 4.5.1` —— rationally standard 型が `𝒞^istr` で保たれる

原文 (FrdI p.86):
> that if C is of rationally standard type (respectively, of standard type), then so is
-/

include G in
/-- ★★★★★★**[FrdI] Remark 4.5.1** —— `𝒞` が rationally standard 型なら `𝒞^istr` もそうである。

★4 条の内訳:

| 条 | 宣言 |
|---|---|
| (a) birationally Frobenius-normalized 型 | `istr_biratFrobNormalizedType` |
| (a) rational 型 | `istr_rationalType` |
| (a) standard 型 | `istr_standardType`(在庫) |
| (b) `(𝒞^un-tr)^birat` の Frobenius-compact 対象 | `istr_unTrBiratCompact`(在庫) |
-/
theorem istr_rationallyStandardType
    {ι : ∀ Y : D, Prime (Φ.val Y) → Pf (Φ.val Y) → ℝ≥0}
    (Hperf : ∀ Y : D, IsPerfFactorialWith (Φ.val Y) (ι Y))
    (h : IsOfRationallyStandardType P G ι) :
    IsOfRationallyStandardType (istrPre P G.core)
      (istr_frobenioid P G.core G) ι where
  biratFrobNormalized := istr_biratFrobNormalizedType P G G.core h.biratFrobNormalized
  rational := istr_rationalType P G G.core Hperf h.biratFrobNormalized h.rational
  standard := istr_standardType (F := G.core) h.standard
  unTrBiratCompact := istr_unTrBiratCompact (F := G.core) G h.unTrBiratCompact

def istr_rationallyStandardType.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 86, item := "Remark 4.5.1",
    sectionId := "frdi-remark-4-5-1" }

end ABC3.Found.FrdI
