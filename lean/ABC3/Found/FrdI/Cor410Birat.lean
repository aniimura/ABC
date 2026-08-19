/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Cor410
import ABC3.Found.FrdI.Prop48Nf

/-!
# [FrdI] Corollary 4.10 —— `Ψ^birat` の構成

原文 (FrdI p.91):
> an equivalence of categories. Then there exists a 1-unique functor

★★原文は「follows immediately from the deﬁnition of `𝒞ᵢ^birat` and the fact that
`Ψ` preserves co-angular pre-steps [cf. Theorem 3.4, (ii)]」と書く。

## ★雛形は `Prop48Nf.lean`

`Proposition 4.8, (ii)` で `naiveFrob`(同じ圏の自己関手)について
`coaPreNfMap → idxBiratNfMap → biratNfMap → psiBirat` を組んである。
★★本ファイルはそれを**2 つの圏の間**へ移すだけである。

★`Ψ` が co-angular pre-step を保ち・反射することは
`Theorem 3.4, (ii)` の内容なので、**仮定として受ける**。
-/

namespace ABC3.Found.FrdI

open CategoryTheory CategoryTheory.Limits

universe v u w u2 v2

/-! ★★**両側の宇宙を揃える** —— `Cocone` の頂点は関手の行き先と同じ宇宙でなければならず、
`HomBirat P₁ G₁ A B : Type (max v u2 v2)` と
`HomBirat P₂ G₂ _ _ : Type (max v' u2' v2')` が食い違うと余錐が組めない(実測)。
★`Corollary 4.10` の用法では両側は同じ大きさなので実害は無い。 -/
variable {D₁ : Type u} [Category.{v} D₁] {C₁ : Type u2} [Category.{v2} C₁]
  {Φ₁ : MonoidOn.{v, u, w} D₁} {P₁ : PreFrobenioid C₁ Φ₁} (G₁ : Frobenioid P₁)
  {D₂ : Type u} [Category.{v} D₂] {C₂ : Type u2} [Category.{v2} C₂]
  {Φ₂ : MonoidOn.{v, u, w} D₂} {P₂ : PreFrobenioid C₂ Φ₂} (G₂ : Frobenioid P₂)

/-! ## ★1. `𝒞^{coa-pre}` の間の関手 -/

variable (Ψ : C₁ ⥤ C₂)
  (hfwd : ∀ {X Y : C₁} (f : X ⟶ Y), coaPreProp P₁ f → coaPreProp P₂ (Ψ.map f))

include hfwd in
/-- ★★**`Ψ` は `𝒞^{coa-pre}` の間の関手を誘導する**。 -/
def coaPrePsi : CoaPre P₁ G₁ ⥤ CoaPre P₂ G₂ where
  obj X := ⟨Ψ.obj X.obj⟩
  map {X Y} f := ⟨Ψ.map f.hom, hfwd f.hom f.property⟩
  map_id X := by
    refine WideSubcategory.hom_ext _ ?_
    show Ψ.map (𝟙 X.obj) = 𝟙 _
    exact CategoryTheory.Functor.map_id _ _
  map_comp {X Y Z} f g := by
    refine WideSubcategory.hom_ext _ ?_
    show Ψ.map (f.hom ≫ g.hom) = Ψ.map f.hom ≫ Ψ.map g.hom
    exact CategoryTheory.Functor.map_comp _ _ _

include hfwd in
/-- ★忠実 —— `Ψ` の忠実性そのもの。 -/
theorem coaPrePsi_faithful [Ψ.Faithful] : (coaPrePsi G₁ G₂ Ψ hfwd).Faithful where
  map_injective {_ _ f g} h :=
    WideSubcategory.hom_ext _ (Ψ.map_injective
      (congrArg (fun t : (coaPrePsi G₁ G₂ Ψ hfwd).obj _ ⟶ _ => t.hom) h))

include hfwd in
/-- ★★充満 —— `Ψ` の充満性 ＋ `coaPreProp` の**反射**。 -/
theorem coaPrePsi_full [Ψ.Full]
    (hbwd : ∀ {X Y : C₁} (f : X ⟶ Y), coaPreProp P₂ (Ψ.map f) → coaPreProp P₁ f) :
    (coaPrePsi G₁ G₂ Ψ hfwd).Full where
  map_surjective {X Y} g := by
    obtain ⟨f₀, hf₀⟩ := Ψ.map_surjective g.hom
    refine ⟨⟨f₀, hbwd f₀ ?_⟩, WideSubcategory.hom_ext _ hf₀⟩
    rw [hf₀]
    exact g.property

include hfwd in
/-- ★本質的全射 —— `Ψ` の本質的全射性を同型ごと包む。 -/
theorem coaPrePsi_essSurj [Ψ.EssSurj] : (coaPrePsi G₁ G₂ Ψ hfwd).EssSurj where
  mem_essImage Y := by
    obtain ⟨A, ⟨ε⟩⟩ := Functor.EssSurj.mem_essImage (F := Ψ) Y.obj
    exact ⟨(⟨A⟩ : CoaPre P₁ G₁), ⟨coaPreIsoOfIso P₂ G₂ ε⟩⟩

include hfwd in
/-- ★★★★**`𝒞^{coa-pre}` の間の関手は圏同値**。 -/
theorem coaPrePsi_isEquivalence [Ψ.IsEquivalence]
    (hbwd : ∀ {X Y : C₁} (f : X ⟶ Y), coaPreProp P₂ (Ψ.map f) → coaPreProp P₁ f) :
    (coaPrePsi G₁ G₂ Ψ hfwd).IsEquivalence where
  faithful := coaPrePsi_faithful G₁ G₂ Ψ hfwd
  full := coaPrePsi_full G₁ G₂ Ψ hfwd hbwd
  essSurj := coaPrePsi_essSurj G₁ G₂ Ψ hfwd

/-! ## ★2. 添字圏の対応 -/

section Idx

variable (A : C₁)

include hfwd in
/-- ★スライスの対応。 -/
def slicePsi : SliceA P₁ G₁ A ⥤ SliceA P₂ G₂ (Ψ.obj A) :=
  Over.post (coaPrePsi G₁ G₂ Ψ hfwd)

include hfwd in
/-- ★添字圏(反対圏)の対応。 -/
def idxBiratPsi : IdxBirat P₁ G₁ A ⥤ IdxBirat P₂ G₂ (Ψ.obj A) :=
  (slicePsi G₁ G₂ Ψ hfwd A).op

include hfwd in
/-- ★★★添字圏の対応も圏同値。 -/
theorem idxBiratPsi_isEquivalence [Ψ.IsEquivalence]
    (hbwd : ∀ {X Y : C₁} (f : X ⟶ Y), coaPreProp P₂ (Ψ.map f) → coaPreProp P₁ f) :
    (idxBiratPsi G₁ G₂ Ψ hfwd A).IsEquivalence := by
  haveI := coaPrePsi_isEquivalence G₁ G₂ Ψ hfwd hbwd
  exact inferInstanceAs
    (Functor.IsEquivalence (Over.post (X := coaPreObj P₁ G₁ A) (coaPrePsi G₁ G₂ Ψ hfwd)).op)

end Idx

def coaPrePsi.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 91,
    item := "Corollary 4.10 — Ψ が誘導する 𝒞^{coa-pre} の間の圏同値",
    sectionId := "frdi-cor-4-10" }

/-! ## ★3. `Hom^birat` の写像

★遷移は**前合成**なので、naturality は `Ψ` の関手性そのもの。 -/

include hfwd in
/-- ★★降ろすための余錐。 -/
noncomputable def biratPsiCocone (A B : C₁) :
    Cocone (homFunctorBirat P₁ G₁ A B) :=
  Cocone.mk (HomBirat P₂ G₂ (Ψ.obj A) (Ψ.obj B))
    { app := fun Z => TypeCat.ofHom fun φ =>
        HomBirat.mk ((idxBiratPsi G₁ G₂ Ψ hfwd A).obj Z) (Ψ.map φ.down)
      naturality := fun Z W u => by
        ext φ
        show HomBirat.mk ((idxBiratPsi G₁ G₂ Ψ hfwd A).obj W)
              (Ψ.map (u.unop.left.hom ≫ φ.down))
          = HomBirat.mk ((idxBiratPsi G₁ G₂ Ψ hfwd A).obj Z) (Ψ.map φ.down)
        rw [CategoryTheory.Functor.map_comp]
        exact HomBirat.mk_map ((idxBiratPsi G₁ G₂ Ψ hfwd A).map u) (Ψ.map φ.down) }

include hfwd in
/-- ★★★`Hom^birat` の写像。 -/
noncomputable def biratPsiMap (A B : C₁) :
    HomBirat P₁ G₁ A B → HomBirat P₂ G₂ (Ψ.obj A) (Ψ.obj B) :=
  fun z => colimit.desc _ (biratPsiCocone G₁ G₂ Ψ hfwd A B) z

include hfwd in
/-- ★代表元での値。 -/
theorem biratPsiMap_mk (A B : C₁) (Z : IdxBirat P₁ G₁ A) (φ : Z.unop.left.obj ⟶ B) :
    biratPsiMap G₁ G₂ Ψ hfwd A B (HomBirat.mk Z φ)
      = HomBirat.mk ((idxBiratPsi G₁ G₂ Ψ hfwd A).obj Z) (Ψ.map φ) := by
  show colimit.desc _ (biratPsiCocone G₁ G₂ Ψ hfwd A B)
      (colimit.ι (homFunctorBirat P₁ G₁ A B) Z (ULift.up φ)) = _
  rw [← types_comp_apply
    (colimit.ι (homFunctorBirat P₁ G₁ A B) Z)
    (colimit.desc _ (biratPsiCocone G₁ G₂ Ψ hfwd A B)), colimit.ι_desc]
  rfl

def biratPsiMap.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 91,
    item := "Corollary 4.10 — Hom^birat の写像",
    sectionId := "frdi-cor-4-10" }

/-! ## ★4. 合成の保存

★★`BiratCat` の合成は Ore 四角形(`biratPull*`)を `choose` で取るので、
両側で選択が違いうる。★添字圏が**細い**(`idxBirat_hom_ext`)ので共通の上界で潰れる。 -/

set_option maxHeartbeats 1000000 in
include hfwd in
/-- ★★★引き戻しの `α` 成分の移送。 -/
theorem biratPsi_pullAlpha_transport {A B : C₁}
    (Z : IdxBirat P₁ G₁ A) (φ : Z.unop.left.obj ⟶ B) (W : IdxBirat P₁ G₁ B)
    {V : IdxBirat P₂ G₂ (Ψ.obj A)}
    (c : (idxBiratPsi G₁ G₂ Ψ hfwd A).obj (biratPullIdx G₁.core Z φ W) ⟶ V)
    (c' : biratPullIdx G₂.core ((idxBiratPsi G₁ G₂ Ψ hfwd A).obj Z) (Ψ.map φ)
      ((idxBiratPsi G₁ G₂ Ψ hfwd B).obj W) ⟶ V)
    (hkey : c.unop.left.hom ≫ Ψ.map (biratPullGamma G₁.core Z φ W)
        = c'.unop.left.hom ≫ biratPullGamma G₂.core
            ((idxBiratPsi G₁ G₂ Ψ hfwd A).obj Z) (Ψ.map φ)
            ((idxBiratPsi G₁ G₂ Ψ hfwd B).obj W)) :
    c.unop.left.hom ≫ Ψ.map (biratPullAlpha G₁.core Z φ W)
      = c'.unop.left.hom ≫ biratPullAlpha G₂.core
          ((idxBiratPsi G₁ G₂ Ψ hfwd A).obj Z) (Ψ.map φ)
          ((idxBiratPsi G₁ G₂ Ψ hfwd B).obj W) := by
  haveI hb : Mono ((idxBiratPsi G₁ G₂ Ψ hfwd B).obj W).unop.hom.hom :=
    G₂.core.preStepMono _ ((idxBiratPsi G₁ G₂ Ψ hfwd B).obj W).unop.hom.property.2
  refine (cancel_mono ((idxBiratPsi G₁ G₂ Ψ hfwd B).obj W).unop.hom.hom).mp ?_
  have h1 := biratPull_sq G₁.core Z φ W
  have h2 := biratPull_sq G₂.core ((idxBiratPsi G₁ G₂ Ψ hfwd A).obj Z) (Ψ.map φ)
    ((idxBiratPsi G₁ G₂ Ψ hfwd B).obj W)
  have e1 : Ψ.map (biratPullAlpha G₁.core Z φ W)
      ≫ ((idxBiratPsi G₁ G₂ Ψ hfwd B).obj W).unop.hom.hom
      = Ψ.map (biratPullGamma G₁.core Z φ W) ≫ Ψ.map φ := by
    have h := congrArg Ψ.map h1
    rw [CategoryTheory.Functor.map_comp] at h
    exact (h.trans (CategoryTheory.Functor.map_comp _ _ _)).symm
  refine Eq.trans (Category.assoc _ _ _) ?_
  refine Eq.trans (congrArg (fun t => c.unop.left.hom ≫ t) e1) ?_
  refine Eq.trans (Category.assoc _ _ _).symm ?_
  refine Eq.trans (congrArg (fun t => t ≫ Ψ.map φ) hkey) ?_
  refine Eq.trans (Category.assoc _ _ _) ?_
  refine Eq.trans (congrArg (fun t => c'.unop.left.hom ≫ t) h2) ?_
  exact (Category.assoc _ _ _).symm

set_option maxHeartbeats 1000000 in
include hfwd in
/-- ★★★代表元での合成の保存。 -/
theorem biratPsiMap_comp_mk {A B E : C₁}
    (Z : IdxBirat P₁ G₁ A) (φ : Z.unop.left.obj ⟶ B)
    (W : IdxBirat P₁ G₁ B) (ψ : W.unop.left.obj ⟶ E) :
    HomBirat.mk ((idxBiratPsi G₁ G₂ Ψ hfwd A).obj (biratPullIdx G₁.core Z φ W))
        (Ψ.map (biratPullAlpha G₁.core Z φ W ≫ ψ))
      = HomBirat.mk
          (biratPullIdx G₂.core ((idxBiratPsi G₁ G₂ Ψ hfwd A).obj Z) (Ψ.map φ)
            ((idxBiratPsi G₁ G₂ Ψ hfwd B).obj W))
          (biratPullAlpha G₂.core ((idxBiratPsi G₁ G₂ Ψ hfwd A).obj Z) (Ψ.map φ)
            ((idxBiratPsi G₁ G₂ Ψ hfwd B).obj W) ≫ Ψ.map ψ) := by
  refine HomBirat.sound _
    (IsFiltered.leftToMax ((idxBiratPsi G₁ G₂ Ψ hfwd A).obj (biratPullIdx G₁.core Z φ W))
      (biratPullIdx G₂.core ((idxBiratPsi G₁ G₂ Ψ hfwd A).obj Z) (Ψ.map φ)
        ((idxBiratPsi G₁ G₂ Ψ hfwd B).obj W)))
    (IsFiltered.rightToMax ((idxBiratPsi G₁ G₂ Ψ hfwd A).obj (biratPullIdx G₁.core Z φ W))
      (biratPullIdx G₂.core ((idxBiratPsi G₁ G₂ Ψ hfwd A).obj Z) (Ψ.map φ)
        ((idxBiratPsi G₁ G₂ Ψ hfwd B).obj W))) ?_
  have hkey := congrArg (fun t => t.unop.left.hom)
    (idxBirat_hom_ext
      ((idxBiratPsi G₁ G₂ Ψ hfwd A).map (biratPullHom G₁.core Z φ W)
        ≫ IsFiltered.leftToMax _ _)
      (biratPullHom G₂.core ((idxBiratPsi G₁ G₂ Ψ hfwd A).obj Z) (Ψ.map φ)
        ((idxBiratPsi G₁ G₂ Ψ hfwd B).obj W) ≫ IsFiltered.rightToMax _ _))
  have hmid := biratPsi_pullAlpha_transport G₁ G₂ Ψ hfwd Z φ W _ _ hkey
  simp only [CategoryTheory.Functor.map_comp, ← Category.assoc]
  exact Eq.trans (Category.assoc _ _ _).symm
    (congrArg (fun t => t ≫ Ψ.map ψ) hmid)

set_option maxHeartbeats 1000000 in
include hfwd in
/-- ★★★★`biratPsiMap` は合成を保つ。 -/
theorem biratPsiMap_comp {A B E : C₁}
    (f : HomBirat P₁ G₁ A B) (g : HomBirat P₁ G₁ B E) :
    biratPsiMap G₁ G₂ Ψ hfwd A E (compBirat P₁ G₁ G₁.core f g)
      = compBirat P₂ G₂ G₂.core (biratPsiMap G₁ G₂ Ψ hfwd A B f)
          (biratPsiMap G₁ G₂ Ψ hfwd B E g) := by
  obtain ⟨Z, φ, rfl⟩ := HomBirat.exists_rep f
  obtain ⟨W, ψ, rfl⟩ := HomBirat.exists_rep g
  rw [compBirat_mk, biratPsiMap_mk, biratPsiMap_mk, biratPsiMap_mk, compBirat_mk]
  exact biratPsiMap_comp_mk G₁ G₂ Ψ hfwd Z φ W ψ

include hfwd in
/-- ★添字の「頂点」も対応する。 -/
theorem idxBiratPsi_one (A : C₁) :
    (idxBiratPsi G₁ G₂ Ψ hfwd A).obj (idxBiratOne P₁ G₁ A)
      = idxBiratOne P₂ G₂ (Ψ.obj A) :=
  congrArg Opposite.op (congrArg Over.mk ((coaPrePsi G₁ G₂ Ψ hfwd).map_id (coaPreObj P₁ G₁ A)))

include hfwd in
/-- ★★★**`𝒞 → 𝒞^birat` と可換** —— これが `map_id` と 1-可換図式を与える。

★対象の等式で `rw` すると motive が壊れるので、
**添字圏の射を 1 本作って `HomBirat.mk_map` で移す**。 -/
theorem biratPsiMap_toHomBirat {A B : C₁} (φ : A ⟶ B) :
    biratPsiMap G₁ G₂ Ψ hfwd A B (toHomBirat (P := P₁) (G := G₁) φ)
      = toHomBirat (P := P₂) (G := G₂) (Ψ.map φ) := by
  let u : idxBiratOne P₂ G₂ (Ψ.obj A)
      ⟶ (idxBiratPsi G₁ G₂ Ψ hfwd A).obj (idxBiratOne P₁ G₁ A) :=
    Quiver.Hom.op (Over.homMk (𝟙 _) (by
      refine (Category.comp_id _).trans ?_
      exact ((coaPrePsi G₁ G₂ Ψ hfwd).map_id (coaPreObj P₁ G₁ A)).symm))
  have h1 : biratPsiMap G₁ G₂ Ψ hfwd A B (HomBirat.mk (idxBiratOne P₁ G₁ A) φ)
      = HomBirat.mk ((idxBiratPsi G₁ G₂ Ψ hfwd A).obj (idxBiratOne P₁ G₁ A)) (Ψ.map φ) :=
    biratPsiMap_mk G₁ G₂ Ψ hfwd A B (idxBiratOne P₁ G₁ A) φ
  refine h1.trans (Eq.trans ?_ (HomBirat.mk_map u (Ψ.map φ)))
  exact congrArg (HomBirat.mk ((idxBiratPsi G₁ G₂ Ψ hfwd A).obj (idxBiratOne P₁ G₁ A)))
    (Category.id_comp (Ψ.map φ)).symm

include hfwd in
/-- ★`biratPsiMap` は恒等を保つ。 -/
theorem biratPsiMap_id (A : C₁) :
    biratPsiMap G₁ G₂ Ψ hfwd A A (toHomBirat (𝟙 A)) = toHomBirat (𝟙 (Ψ.obj A)) :=
  (biratPsiMap_toHomBirat G₁ G₂ Ψ hfwd (𝟙 A)).trans
    (congrArg (toHomBirat (P := P₂) (G := G₂)) (CategoryTheory.Functor.map_id Ψ A))

/-! ## ★5. ★★★★★`Ψ^birat` -/

include hfwd in
/-- ★★★★★**[FrdI] Corollary 4.10** —— `Ψ^birat : 𝒞₁^birat ⥤ 𝒞₂^birat`。 -/
noncomputable def psiBiratCor : BiratCat P₁ G₁ ⥤ BiratCat P₂ G₂ where
  obj A := Ψ.obj A
  map {A B} f := biratPsiMap G₁ G₂ Ψ hfwd A B f
  map_id A := biratPsiMap_id G₁ G₂ Ψ hfwd A
  map_comp f g := biratPsiMap_comp G₁ G₂ Ψ hfwd f g

def psiBiratCor.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 91,
    item := "Corollary 4.10 — Ψ^birat の構成",
    sectionId := "frdi-cor-4-10" }

/-! ## ★6. 1-可換図式と圏同値性 -/

include hfwd in
/-- ★★★★**1-可換図式** —— `𝒞ᵢ → 𝒞ᵢ^birat` と可換(**厳密な等式**)。 -/
theorem psiBiratCor_square :
    Ψ ⋙ toBiratCat P₂ G₂ = toBiratCat P₁ G₁ ⋙ psiBiratCor G₁ G₂ Ψ hfwd := by
  refine CategoryTheory.Functor.ext (fun A => rfl) (fun A B φ => ?_)
  show (toBiratCat P₂ G₂).map (Ψ.map φ)
    = eqToHom rfl ≫ biratPsiMap G₁ G₂ Ψ hfwd A B ((toBiratCat P₁ G₁).map φ) ≫ eqToHom rfl
  simp only [eqToHom_refl, Category.id_comp, Category.comp_id]
  exact (biratPsiMap_toHomBirat G₁ G₂ Ψ hfwd φ).symm

include hfwd in
/-- ★★★`Ψ^birat` は本質的全射。 -/
theorem psiBiratCor_essSurj [Ψ.EssSurj] : (psiBiratCor G₁ G₂ Ψ hfwd).EssSurj := by
  constructor
  intro B
  obtain ⟨A, ⟨ε⟩⟩ := Functor.EssSurj.mem_essImage (F := Ψ) (biratDown P₂ G₂ B)
  exact ⟨A, ⟨(toBiratCat P₂ G₂).mapIso ε⟩⟩

set_option maxHeartbeats 1000000 in
include hfwd in
/-- ★★★★`Ψ^birat` は忠実。★添字圏が圏同値なので共通の上界を像の中に取り直せる。 -/
theorem psiBiratCor_faithful [Ψ.IsEquivalence]
    (hbwd : ∀ {X Y : C₁} (f : X ⟶ Y), coaPreProp P₂ (Ψ.map f) → coaPreProp P₁ f) :
    (psiBiratCor G₁ G₂ Ψ hfwd).Faithful where
  map_injective {X Y} {x y} h := by
    haveI := idxBiratPsi_isEquivalence G₁ G₂ Ψ hfwd (biratDown P₁ G₁ X) hbwd
    obtain ⟨Z, φ, rfl⟩ := HomBirat.exists_rep x
    obtain ⟨W, ψ, rfl⟩ := HomBirat.exists_rep y
    have h' : HomBirat.mk ((idxBiratPsi G₁ G₂ Ψ hfwd (biratDown P₁ G₁ X)).obj Z) (Ψ.map φ)
        = HomBirat.mk ((idxBiratPsi G₁ G₂ Ψ hfwd (biratDown P₁ G₁ X)).obj W) (Ψ.map ψ) :=
      (biratPsiMap_mk G₁ G₂ Ψ hfwd _ _ Z φ).symm.trans
        (h.trans (biratPsiMap_mk G₁ G₂ Ψ hfwd _ _ W ψ))
    obtain ⟨V, u, v, hV⟩ := HomBirat.eq_iff.mp h'
    obtain ⟨V₀, ⟨θ⟩⟩ := Functor.EssSurj.mem_essImage
      (F := idxBiratPsi G₁ G₂ Ψ hfwd (biratDown P₁ G₁ X)) V
    obtain ⟨u₀, hu₀⟩ :=
      (idxBiratPsi G₁ G₂ Ψ hfwd (biratDown P₁ G₁ X)).map_surjective (u ≫ θ.inv)
    obtain ⟨v₀, hv₀⟩ :=
      (idxBiratPsi G₁ G₂ Ψ hfwd (biratDown P₁ G₁ X)).map_surjective (v ≫ θ.inv)
    have hu : Ψ.map u₀.unop.left.hom = θ.inv.unop.left.hom ≫ u.unop.left.hom :=
      congrArg (fun t : (idxBiratPsi G₁ G₂ Ψ hfwd (biratDown P₁ G₁ X)).obj Z ⟶
        (idxBiratPsi G₁ G₂ Ψ hfwd (biratDown P₁ G₁ X)).obj V₀ => t.unop.left.hom) hu₀
    have hv : Ψ.map v₀.unop.left.hom = θ.inv.unop.left.hom ≫ v.unop.left.hom :=
      congrArg (fun t : (idxBiratPsi G₁ G₂ Ψ hfwd (biratDown P₁ G₁ X)).obj W ⟶
        (idxBiratPsi G₁ G₂ Ψ hfwd (biratDown P₁ G₁ X)).obj V₀ => t.unop.left.hom) hv₀
    refine HomBirat.eq_iff.mpr ⟨V₀, u₀, v₀, Ψ.map_injective ?_⟩
    have e3 : Ψ.map (u₀.unop.left.hom ≫ φ) = Ψ.map u₀.unop.left.hom ≫ Ψ.map φ :=
      CategoryTheory.Functor.map_comp _ _ _
    have e4 : Ψ.map (v₀.unop.left.hom ≫ ψ) = Ψ.map v₀.unop.left.hom ≫ Ψ.map ψ :=
      CategoryTheory.Functor.map_comp _ _ _
    have g1 : Ψ.map u₀.unop.left.hom ≫ Ψ.map φ
        = (θ.inv.unop.left.hom ≫ u.unop.left.hom) ≫ Ψ.map φ :=
      congrArg (fun t => t ≫ Ψ.map φ) hu
    have g3 : θ.inv.unop.left.hom ≫ (u.unop.left.hom ≫ Ψ.map φ)
        = θ.inv.unop.left.hom ≫ (v.unop.left.hom ≫ Ψ.map ψ) :=
      congrArg (fun t => θ.inv.unop.left.hom ≫ t) hV
    have g5 : (θ.inv.unop.left.hom ≫ v.unop.left.hom) ≫ Ψ.map ψ
        = Ψ.map v₀.unop.left.hom ≫ Ψ.map ψ :=
      congrArg (fun t => t ≫ Ψ.map ψ) hv.symm
    exact e3.trans (g1.trans ((Category.assoc _ _ _).trans
      (g3.trans ((Category.assoc _ _ _).symm.trans (g5.trans e4.symm)))))

include hfwd in
/-- ★★★★`Ψ^birat` は充満。 -/
theorem psiBiratCor_full [Ψ.IsEquivalence]
    (hbwd : ∀ {X Y : C₁} (f : X ⟶ Y), coaPreProp P₂ (Ψ.map f) → coaPreProp P₁ f) :
    (psiBiratCor G₁ G₂ Ψ hfwd).Full where
  map_surjective {X Y} h := by
    haveI := idxBiratPsi_isEquivalence G₁ G₂ Ψ hfwd (biratDown P₁ G₁ X) hbwd
    obtain ⟨V, χ, rfl⟩ := HomBirat.exists_rep h
    obtain ⟨V₀, ⟨θ⟩⟩ := Functor.EssSurj.mem_essImage
      (F := idxBiratPsi G₁ G₂ Ψ hfwd (biratDown P₁ G₁ X)) V
    obtain ⟨χ₀, hχ₀⟩ := Ψ.map_surjective (θ.inv.unop.left.hom ≫ χ)
    refine ⟨HomBirat.mk V₀ χ₀, ?_⟩
    refine (biratPsiMap_mk G₁ G₂ Ψ hfwd _ _ V₀ χ₀).trans ?_
    exact (congrArg (HomBirat.mk
      ((idxBiratPsi G₁ G₂ Ψ hfwd (biratDown P₁ G₁ X)).obj V₀)) hχ₀).trans
      (HomBirat.mk_map θ.inv χ)

include hfwd in
/-- ★★★★★**[FrdI] Corollary 4.10** —— `Ψ^birat` は**圏同値**。 -/
theorem psiBiratCor_isEquivalence [Ψ.IsEquivalence]
    (hbwd : ∀ {X Y : C₁} (f : X ⟶ Y), coaPreProp P₂ (Ψ.map f) → coaPreProp P₁ f) :
    (psiBiratCor G₁ G₂ Ψ hfwd).IsEquivalence where
  faithful := psiBiratCor_faithful G₁ G₂ Ψ hfwd hbwd
  full := psiBiratCor_full G₁ G₂ Ψ hfwd hbwd
  essSurj := psiBiratCor_essSurj G₁ G₂ Ψ hfwd

/-- ★★★★★★**[FrdI] Corollary 4.10** —— 条なしの locator。 -/
def cor_4_10.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 90, item := "Corollary 4.10",
    sectionId := "frdi-cor-4-10" }

end ABC3.Found.FrdI
