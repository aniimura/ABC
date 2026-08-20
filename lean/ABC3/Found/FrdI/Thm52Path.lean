/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Thm51FrTr
import ABC3.Found.FrdI.Thm52Birat

/-!
# [FrdI] Theorem 5.2, (iv) の第 1 段 —— `F-𝒫-path` の圏

原文 (FrdI p.101):
> such that A ∈ Ob(P) as an F P-path for C. Write

原文 (FrdI p.101):
> for the category C whose objects are objects of C equipped with an F P-path, and

★★**なぜ path つきの圏を経由するのか**:
`𝒞 → 𝒞^model` を直に作ろうとすると、対象の類 `α` と射の単元成分 `u` に
**選択**が入って関手性が壊れる。★path を対象の一部にすれば、どちらも
**path から決まる**ので関手になる。★そして `𝒞̃ → 𝒞` は忘却で圏同値である
(path は存在し、射は `𝒞` の射そのもの)。

★本ファイルではその第 1 段:
* `FPPath` —— path の型
* `PathCat` —— path つきの対象のなす圏
* `pathForget` が圏同値であること
* path が定める **`Φ^gp(Base X)` の類**(model の `α`)と、
  `𝒞^birat` での同型 `X^birat ≅ ref^birat`(その `Div^gp` がちょうど類)
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ}

/-! ## ★1. `F-𝒫-path` -/

/-- ★★**`F-𝒫-path`** —— `Theorem 5.2, (iv)` の証明で使う「pre-step の対」。

原文 (FrdI p.101):
> (B → A, A → C)
-/
structure FPPath (S : BaseSection P) (X : C) where
  /-- `𝒫` の対象 `A` -/
  ref : C
  /-- `A ∈ Ob(𝒫)` -/
  ref_mem : S.objP ref
  /-- 頂点 `B` -/
  vertex : C
  /-- `B → A` -/
  toRef : vertex ⟶ ref
  /-- `B → X` -/
  toObj : vertex ⟶ X
  toRef_preStep : IsPreStep P toRef
  toObj_preStep : IsPreStep P toObj

/-- ★★`F-𝒫-path` つきの対象のなす圏。射は `𝒞` の射そのもの。 -/
structure PathCat (S : BaseSection P) where
  obj : C
  path : FPPath S obj

instance (S : BaseSection P) : Category.{v2} (PathCat S) where
  Hom X Y := X.obj ⟶ Y.obj
  id X := 𝟙 X.obj
  comp f g := f ≫ g

/-- ★忘却関手 `𝒞̃ → 𝒞`。 -/
def pathForget (S : BaseSection P) : PathCat S ⥤ C where
  obj X := X.obj
  map f := f
  map_id _ := rfl
  map_comp _ _ := rfl

/-- ★どの対象にも `F-𝒫-path` が付く。 -/
theorem exists_fpPath (G : Frobenioid P) (S : BaseSection P) (X : C) :
    Nonempty (FPPath S X) := by
  obtain ⟨A, hA, ⟨e⟩⟩ := S.essSurjP ((P.toElem.obj X).base)
  haveI : IsIso e.hom := e.isIso_hom
  obtain ⟨W, φ, ψ, hφ, hψ, -⟩ := G.core.preStepSpan A X e.hom e.isIso_hom
  exact ⟨⟨A, hA, W, φ, ψ, hφ, hψ⟩⟩

/-- ★★★**忘却関手 `𝒞̃ → 𝒞` は圏同値**。

原文 (FrdI p.101):
> we have a natural functor C → C [obtained by forgetting the F P-paths],
-/
theorem pathForget_isEquivalence (G : Frobenioid P) (S : BaseSection P) :
    (pathForget S).IsEquivalence := by
  refine ⟨⟨fun {X Y} {f g} h => h⟩, ⟨fun {X Y} f => ⟨f, rfl⟩⟩, ⟨fun X => ?_⟩⟩
  obtain ⟨p⟩ := exists_fpPath G S X
  exact ⟨⟨X, p⟩, ⟨Iso.refl X⟩⟩

/-! ## ★2. path が定める類 -/

/-- ★★`F-𝒫-path` が定める **`Φ^gp(Base X)` の類** —— model Frobenioid の `α`。

★★模型で検算すると: `X = (d, α)`、`ref = (d, 0)`、`vertex = (d, w)` のとき
model Frobenioid の射の条件 `deg·α_A + Div = Φ(Base)(α_B) + Div_B(u)` から
`Div(toObj) = α − w`、`Div(toRef) = −w` なので
`Φ(Base toObj)⁻¹(Div(toRef) − Div(toObj)) = −α`。
★**したがって `α` は `spanCls` の符号を反転したもの**である
(この符号は模型 `𝒟 = pt`, `Φ = ℕ`, `B = 0` で `birat_divGp_sub_mem` と
突き合わせて確定した)。 -/
noncomputable def FPPath.cls {S : BaseSection P} {X : C} (p : FPPath S X) :
    Gp (Φ.val (P.toElem.obj X).base) :=
  -spanCls p.toObj p.toObj_preStep.2 p.toRef

/-- ★`[δ₁]⁻¹ ≫ [δ₂]` は `𝒞^birat` の同型(終域が違ってもよい版)。 -/
theorem birat_mk_isIso' (G : Frobenioid P) {W A B : C}
    (δ₁ : W ⟶ A) (δ₂ : W ⟶ B) (hc₁ : IsCoAngular P δ₁) (hs₁ : IsPreStep P δ₁)
    (hc₂ : IsCoAngular P δ₂) (hs₂ : IsPreStep P δ₂) :
    IsIso (show (toBiratCat P G).obj A ⟶ (toBiratCat P G).obj B from
      HomBirat.mk (idxBiratMk P G δ₁ hc₁ hs₁) δ₂) := by
  haveI h1 : IsIso ((toBiratCat P G).map δ₁) := birat_isIso_of_coaPre δ₁ hc₁ hs₁
  haveI h2 : IsIso ((toBiratCat P G).map δ₂) := birat_isIso_of_coaPre δ₂ hc₂ hs₂
  haveI : IsIso ((toBiratCat P G).map δ₁ ≫ (show (toBiratCat P G).obj A ⟶ (toBiratCat P G).obj B
      from HomBirat.mk (idxBiratMk P G δ₁ hc₁ hs₁) δ₂)) := by
    have hcomp : ((toBiratCat P G).map δ₁ ≫ (show (toBiratCat P G).obj A ⟶ (toBiratCat P G).obj B
        from HomBirat.mk (idxBiratMk P G δ₁ hc₁ hs₁) δ₂)) = (toBiratCat P G).map δ₂ :=
      birat_toHom_comp_mk (G := G) δ₁ hc₁ hs₁ δ₂
    rw [hcomp]; exact h2
  exact IsIso.of_isIso_comp_left ((toBiratCat P G).map δ₁) _

/-- ★★`F-𝒫-path` が `𝒞^birat` で定める同型 `X^birat ≅ ref^birat`。 -/
noncomputable def FPPath.biratIso (G : Frobenioid P) (hiso : ∀ Y : C, IsIsotropic P Y)
    {S : BaseSection P} {X : C} (p : FPPath S X) : HomBirat P G X p.ref :=
  HomBirat.mk (idxBiratMk P G p.toObj (prop_1_4_i P _ (fun Y _ => hiso Y)) p.toObj_preStep)
    p.toRef

theorem FPPath.biratIso_isIso (G : Frobenioid P) (hiso : ∀ Y : C, IsIsotropic P Y)
    {S : BaseSection P} {X : C} (p : FPPath S X) :
    IsIso (show (toBiratCat P G).obj X ⟶ (toBiratCat P G).obj p.ref from p.biratIso G hiso) :=
  birat_mk_isIso' G p.toObj p.toRef (prop_1_4_i P _ (fun Y _ => hiso Y)) p.toObj_preStep
    (prop_1_4_i P _ (fun Y _ => hiso Y)) p.toRef_preStep

/-- ★★★**その `Div^gp` はちょうど類の符号反転**。

★model の `α` は `−Div^gp(π_X)` である。 -/
theorem FPPath.biratDivGp_biratIso (G : Frobenioid P) (hiso : ∀ Y : C, IsIsotropic P Y)
    {S : BaseSection P} {X : C} (p : FPPath S X) :
    biratDivGp (p.biratIso G hiso) = -p.cls := by
  rw [FPPath.biratIso, biratDivGp_mk, FPPath.cls, neg_neg,
    spanCls_eq_sliceDivGpOf _ _ _ p.toRef_preStep.1]
  exact sliceDivGpOf_congr rfl _ _ _

/-- ★その底。 -/
theorem FPPath.biratBase_biratIso (G : Frobenioid P) (hiso : ∀ Y : C, IsIsotropic P Y)
    {S : BaseSection P} {X : C} (p : FPPath S X) :
    biratBase (p.biratIso G hiso)
      = @inv _ _ _ _ (P.Base p.toObj) p.toObj_preStep.2 ≫ P.Base p.toRef := by
  rw [FPPath.biratIso, biratBase_mk]
  exact sliceBaseOf_eq _ _ _

/-- ★その次数は 1。 -/
theorem FPPath.biratDeg_biratIso (G : Frobenioid P) (hiso : ∀ Y : C, IsIsotropic P Y)
    {S : BaseSection P} {X : C} (p : FPPath S X) :
    biratDeg (p.biratIso G hiso) = 1 := by
  rw [FPPath.biratIso, biratDeg_mk]
  exact p.toRef_preStep.1

/-! ### ★出典の紐付け(`.src`) -/

/-- ★locator —— `Theorem 5.2, (iv)` の `F-𝒫-path` の圏。 -/
def pathForget_isEquivalence.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 101,
    item := "Theorem 5.2, (iv) — F-𝒫-path の圏と 𝒞 との同値",
    sectionId := "frdi-thm-5-2" }

end ABC3.Found.FrdI
