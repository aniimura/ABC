/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Cor54
import ABC3.Found.FrdI.Prop53Base

/-!
# [FrdI] Corollary 5.4 の最後の文 —— `Ψ^rlf` の作り方は完全化と 1-可換

原文 (FrdI p.104):
> with the 1-commutative diagram of Proposition 5.3 [involving perfections, unit-

★`Corollary 5.4` の最後の文は「`Ψ` から `Ψ^rlf` を作る操作は、
`Proposition 5.3` の 1-可換図式(完全化・単位自明化などが出てくるもの)と
**1-可換である**」と言う。

★★我々の書き方では、`Proposition 5.3` の図式の下の行の右の矢印は
**係数の拡大 `σ : S ⟶ S'`**(`Prop53Base.lean` の `scBaseFunctor`)であり、
`Ψ^rlf` は**因子の単系の同型 `η` を `⊗_S` したもの**(`Cor54.lean` の `psiSc`)である。
したがって主張は次の四角形の可換性になる:

```
(𝒞₁^un-tr)^pf --Ψ^pf--> (𝒞₂^un-tr)^pf
     |                        |
   係数の拡大               係数の拡大
     v                        v
   𝒞₁^rlf     --Ψ^rlf-->    𝒞₂^rlf
```

★★★**中身は `scBase σ` と `scMap θ` が可換であること 1 本だけ**である ——
どちらもテンソル積の別々の因子に効くのだから当たり前で、
それが `scBase_scMap`(`Prop53Base.lean`)である。

## ★まだ実装していない条(記録)

`Corollary 5.4` の残りは `Cor54.lean` の冒頭に書いた通り
(1. `hbirat` を `Corollary 4.11` から導く / 2. 1-一意性 / 3. rigidity /
4. `Proposition 5.3` の**縦の矢印**との 1-可換性)。
★本ファイルが閉じるのは最後の文(**完全化・実化との 1-可換性**)である。
-/

namespace ABC3.Found.FrdI

open CategoryTheory
open scoped TensorProduct

universe v u w u2 v2

section Compat

variable {S S' : Type} [CommSemiring S] [CommSemiring S']

variable {D₁ : Type u} [Category.{v} D₁] {C₁ : Type u2} [Category.{v2} C₁]
  {Φ₁ : MonoidOn.{v, u, w} D₁} {P₁ : PreFrobenioid C₁ Φ₁}
  {D₂ : Type u} [Category.{v} D₂] {C₂ : Type u2} [Category.{v2} C₂]
  {Φ₂ : MonoidOn.{v, u, w} D₂} {P₂ : PreFrobenioid C₂ Φ₂}

/-- ★★★★★★**`Corollary 5.4` の最後の文**(`ModelData` の射の段)——
「`Ψ` から `Ψ^rlf` を作る操作」と「係数の拡大」は可換。

★中身は `scBase σ (scMap θ x) = scMap θ (scBase σ x)`(`scBase_scMap`)だけである。 -/
theorem scBaseHom_compOver_scModelHomOver (σ : S →+* S')
    (ΨB : D₁ ⥤ D₂) (η : Φ₁.functor ≅ ΨB.op ⋙ Φ₂.functor)
    (G₁ : Frobenioid P₁) (hiso₁ : ∀ Y : C₁, IsIsotropic P₁ Y)
    (hfn₁ : ∀ X : BiratCat P₁ G₁, IsFrobeniusNormalized (biratPre P₁ G₁) X)
    (hcharInj₁ : ∀ {A B : D₁} (α : B ⟶ A),
      IsCharacteristicallyInjective (scMap (S := S) (Φ₁.map α)))
    (hint₁ : ∀ A : D₁, IsIntegralMonoid (ScT S (Φ₁.val A)))
    (hcharInj₁' : ∀ {A B : D₁} (α : B ⟶ A),
      IsCharacteristicallyInjective (scMap (S := S') (Φ₁.map α)))
    (hint₁' : ∀ A : D₁, IsIntegralMonoid (ScT S' (Φ₁.val A)))
    (G₂ : Frobenioid P₂) (hiso₂ : ∀ Y : C₂, IsIsotropic P₂ Y)
    (hfn₂ : ∀ X : BiratCat P₂ G₂, IsFrobeniusNormalized (biratPre P₂ G₂) X)
    (hcharInj₂ : ∀ {A B : D₂} (α : B ⟶ A),
      IsCharacteristicallyInjective (scMap (S := S) (Φ₂.map α)))
    (hint₂ : ∀ A : D₂, IsIntegralMonoid (ScT S (Φ₂.val A)))
    (hcharInj₂' : ∀ {A B : D₂} (α : B ⟶ A),
      IsCharacteristicallyInjective (scMap (S := S') (Φ₂.map α)))
    (hint₂' : ∀ A : D₂, IsIntegralMonoid (ScT S' (Φ₂.val A)))
    (hfsmD₁ : IsOfFSMType D₁) (hfsmD₂ : IsOfFSMType D₂)
    (hbirat : ∀ (d : D₁) (y : Gp (Φ₁.val d)), y ∈ phiBiratOn G₁ d →
      gpMap _ (phiIsoApp ΨB η d) y ∈ phiBiratOn G₂ (ΨB.obj d)) :
    (scBaseHom σ G₁ hiso₁ hfn₁ hcharInj₁ hint₁ hcharInj₁' hint₁' hfsmD₁).compOver
        (scModelHomOver S' ΨB η G₁ hiso₁ hfn₁ hcharInj₁' hint₁'
          G₂ hiso₂ hfn₂ hcharInj₂' hint₂' hfsmD₁ hfsmD₂ hbirat)
      = (scModelHomOver S ΨB η G₁ hiso₁ hfn₁ hcharInj₁ hint₁
          G₂ hiso₂ hfn₂ hcharInj₂ hint₂ hfsmD₁ hfsmD₂ hbirat).compHom
        (scBaseHom σ G₂ hiso₂ hfn₂ hcharInj₂ hint₂ hcharInj₂' hint₂' hfsmD₂) := by
  refine ModelDataHomOver.ext' (fun d => ?_) (fun d => ?_)
  · exact AddMonoidHom.ext (fun x => (scBase_scMap σ (phiIsoApp ΨB η d) x).symm)
  · exact AddMonoidHom.ext (fun x =>
      Subtype.ext (gpMap_scBase_scMap σ (phiIsoApp ΨB η d) _).symm)

/-- ★★★★★★**`Corollary 5.4` の最後の文**(関手の等式として)。

```
(𝒞₁^un-tr)^pf --Ψ^pf--> (𝒞₂^un-tr)^pf --係数の拡大--> 𝒞₂^rlf
        ＝
(𝒞₁^un-tr)^pf --係数の拡大--> 𝒞₁^rlf --Ψ^rlf--> 𝒞₂^rlf
```
-/
theorem scBaseFunctor_comp_psiSc (σ : S →+* S')
    (ΨB : D₁ ⥤ D₂) (η : Φ₁.functor ≅ ΨB.op ⋙ Φ₂.functor)
    (G₁ : Frobenioid P₁) (hiso₁ : ∀ Y : C₁, IsIsotropic P₁ Y)
    (hfn₁ : ∀ X : BiratCat P₁ G₁, IsFrobeniusNormalized (biratPre P₁ G₁) X)
    (hcharInj₁ : ∀ {A B : D₁} (α : B ⟶ A),
      IsCharacteristicallyInjective (scMap (S := S) (Φ₁.map α)))
    (hint₁ : ∀ A : D₁, IsIntegralMonoid (ScT S (Φ₁.val A)))
    (hcharInj₁' : ∀ {A B : D₁} (α : B ⟶ A),
      IsCharacteristicallyInjective (scMap (S := S') (Φ₁.map α)))
    (hint₁' : ∀ A : D₁, IsIntegralMonoid (ScT S' (Φ₁.val A)))
    (G₂ : Frobenioid P₂) (hiso₂ : ∀ Y : C₂, IsIsotropic P₂ Y)
    (hfn₂ : ∀ X : BiratCat P₂ G₂, IsFrobeniusNormalized (biratPre P₂ G₂) X)
    (hcharInj₂ : ∀ {A B : D₂} (α : B ⟶ A),
      IsCharacteristicallyInjective (scMap (S := S) (Φ₂.map α)))
    (hint₂ : ∀ A : D₂, IsIntegralMonoid (ScT S (Φ₂.val A)))
    (hcharInj₂' : ∀ {A B : D₂} (α : B ⟶ A),
      IsCharacteristicallyInjective (scMap (S := S') (Φ₂.map α)))
    (hint₂' : ∀ A : D₂, IsIntegralMonoid (ScT S' (Φ₂.val A)))
    (hfsmD₁ : IsOfFSMType D₁) (hfsmD₂ : IsOfFSMType D₂)
    (hbirat : ∀ (d : D₁) (y : Gp (Φ₁.val d)), y ∈ phiBiratOn G₁ d →
      gpMap _ (phiIsoApp ΨB η d) y ∈ phiBiratOn G₂ (ΨB.obj d)) :
    scBaseFunctor σ G₁ hiso₁ hfn₁ hcharInj₁ hint₁ hcharInj₁' hint₁' hfsmD₁
        ⋙ psiSc S' ΨB η G₁ hiso₁ hfn₁ hcharInj₁' hint₁'
          G₂ hiso₂ hfn₂ hcharInj₂' hint₂' hfsmD₁ hfsmD₂ hbirat
      = psiSc S ΨB η G₁ hiso₁ hfn₁ hcharInj₁ hint₁
          G₂ hiso₂ hfn₂ hcharInj₂ hint₂ hfsmD₁ hfsmD₂ hbirat
        ⋙ scBaseFunctor σ G₂ hiso₂ hfn₂ hcharInj₂ hint₂ hcharInj₂' hint₂' hfsmD₂ := by
  unfold scBaseFunctor psiSc
  rw [← ModelDataHom.compOver_functor, ← ModelDataHomOver.compHom_functor,
    scBaseHom_compOver_scModelHomOver]

end Compat

/-! ### ★出典の紐付け -/

/-- ★★★★locator —— `Corollary 5.4` の最後の文
「`Ψ^rlf` の作り方は完全化・実化と 1-可換」。 -/
def scBaseFunctor_comp_psiSc.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 104,
    item := "Corollary 5.4 — Ψ^rlf の作り方は完全化・実化と 1-可換",
    sectionId := "frdi-cor-5-4" }

end ABC3.Found.FrdI
