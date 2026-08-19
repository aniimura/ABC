/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Rmk451Birat
import ABC3.Found.FrdI.Prop44Core
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

include F G in
/-- ★★`biratIstrMap` は底を保つ。 -/
theorem biratIstrMap_base (A B : Istr P)
    (z : HomBirat (istrPre P F) (istr_frobenioid P F G) A B) :
    biratBase (biratIstrMap P G F A B z) = biratBase z := by
  obtain ⟨Z, φ, rfl⟩ := HomBirat.exists_rep z
  rw [biratIstrMap_mk, biratBase_mk, biratBase_mk, sliceBaseOf_eq, sliceBaseOf_eq]

def biratIstrFunctor.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 86,
    item := "Remark 4.5.1 — (𝒞^istr)^birat ⥤ 𝒞^birat は充満忠実で構造を保つ",
    sectionId := "frdi-remark-4-5-1" }

end ABC3.Found.FrdI
