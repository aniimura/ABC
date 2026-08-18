/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Thm34EndBs
import ABC3.Found.FrdI.Def31Pf

/-!
# [FrdI] Theorem 3.4, (iii) —— `Ψ^pf` の構成

原文 (FrdI p.64):
> Definition 3.1, (iii)] that we obtain a 1-unique 1-commutative diagram as in the

## ★★測ったこと(2026-08-19)

`𝒞^pf` は `Def31Pf.lean` で**余極限**として構成されている:

- 対象は `𝒞` の対象そのもの(`PfCat _P _F := C`)
- 射は添字圏 `IdxPf P F A B = Under (biFrObj P F A B)` 上の余極限
- 添字圏の底は `BiFr = WideSubcategory (biFrPropOf P F)`、すなわち
  「同次数の Frobenius 型射の対」がなす `𝒞 × 𝒞` の広い部分圏

★★★**`iv-endbsiso` の `plBkMap` と同じ形**である ——
射のクラスを保つ関手が広い部分圏の間の関手を誘導し、
`Under.post` / `Over.post` でスライスへ降りる。
★本ファイルはその第 1 段(`biFrMap`)を取る。
-/

universe v u w u2 v2

namespace ABC3.Found.FrdI

open CategoryTheory

variable {D₁ : Type u} [Category.{v} D₁] {C₁ : Type u2} [Category.{v2} C₁]
  {Φ₁ : MonoidOn.{v, u, w} D₁} {P₁ : PreFrobenioid C₁ Φ₁} (F₁ : FrobenioidCore P₁)
  {D₂ : Type u} [Category.{v} D₂] {C₂ : Type u2} [Category.{v2} C₂]
  {Φ₂ : MonoidOn.{v, u, w} D₂} {P₂ : PreFrobenioid C₂ Φ₂} (F₂ : FrobenioidCore P₂)

/-! ## ★段 1 —— `𝒞^{bi-Fr}` の間の関手 -/

/-- ★★★★**Frobenius 型と次数の一致を保つ関手は `𝒞^{bi-Fr}` の間の関手を誘導する**。

★`biFrProp` は「両成分が Frobenius 型」かつ「次数が一致」なので、
`Ψ` について要るのはちょうどその 2 つ。
★★これは `iv-endbsiso` の `plBkMap` と**同じ形**である。 -/
def biFrMap (Ψ : C₁ ⥤ C₂)
    (hFT : ∀ {X Y : C₁} (f : X ⟶ Y), IsFrobeniusType P₁ f → IsFrobeniusType P₂ (Ψ.map f))
    (hdegEq : ∀ {X Y X' Y' : C₁} (f : X ⟶ Y) (g : X' ⟶ Y'),
      P₁.degFr f = P₁.degFr g → P₂.degFr (Ψ.map f) = P₂.degFr (Ψ.map g)) :
    BiFr P₁ F₁ ⥤ BiFr P₂ F₂ where
  obj X := ⟨(Ψ.obj X.obj.1, Ψ.obj X.obj.2)⟩
  map f := ⟨(Ψ.map f.hom.1, Ψ.map f.hom.2),
    hFT _ f.property.1, hFT _ f.property.2.1, hdegEq _ _ f.property.2.2⟩
  map_id X := WideSubcategory.hom_ext _ (Prod.ext (Ψ.map_id X.obj.1) (Ψ.map_id X.obj.2))
  map_comp f g := WideSubcategory.hom_ext _
    (Prod.ext (Ψ.map_comp f.hom.1 g.hom.1) (Ψ.map_comp f.hom.2 g.hom.2))

@[simp]
theorem biFrMap_obj (Ψ : C₁ ⥤ C₂) (hFT) (hdegEq) (X : BiFr P₁ F₁) :
    (biFrMap F₁ F₂ Ψ hFT hdegEq).obj X = ⟨(Ψ.obj X.obj.1, Ψ.obj X.obj.2)⟩ := rfl

/-- ★★**添字圏の底が対応する** —— `biFrObj` は `biFrMap` で `biFrObj` に写る。 -/
theorem biFrMap_biFrObj (Ψ : C₁ ⥤ C₂) (hFT) (hdegEq) (A B : C₁) :
    (biFrMap F₁ F₂ Ψ hFT hdegEq).obj (biFrObj P₁ F₁ A B)
      = biFrObj P₂ F₂ (Ψ.obj A) (Ψ.obj B) := rfl

/-- ★★★★**添字圏の間の関手** —— `Under.post` で降ろすだけ。

★`IdxPf P F A B = Under (biFrObj P F A B)` なので、
`biFrMap` を `Under.post` に通せばそのまま添字圏の写像になる。 -/
def idxPfMap (Ψ : C₁ ⥤ C₂)
    (hFT : ∀ {X Y : C₁} (f : X ⟶ Y), IsFrobeniusType P₁ f → IsFrobeniusType P₂ (Ψ.map f))
    (hdegEq : ∀ {X Y X' Y' : C₁} (f : X ⟶ Y) (g : X' ⟶ Y'),
      P₁.degFr f = P₁.degFr g → P₂.degFr (Ψ.map f) = P₂.degFr (Ψ.map g)) (A B : C₁) :
    IdxPf P₁ F₁ A B ⥤ IdxPf P₂ F₂ (Ψ.obj A) (Ψ.obj B) :=
  Under.post (X := biFrObj P₁ F₁ A B) (biFrMap F₁ F₂ Ψ hFT hdegEq)

def biFrMap.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 64,
    item := "Theorem 3.4, (iii) — Ψ^pf の添字圏の写像",
    sectionId := "frdi-thm-3-4" }

end ABC3.Found.FrdI
