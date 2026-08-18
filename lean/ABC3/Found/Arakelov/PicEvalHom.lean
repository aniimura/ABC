import ABC3.Found.Arakelov.PicEvalBil

/-!
# Arakelov (B2) 第 175 ブロック —— ★★★★★★★**評価射**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★★`F ⊗ F^∨ ⟶ 𝟙_` が書けた

    evHom : F ⊗ (dualPresheaf F) ⟶ 𝟙_
      app W = TensorProduct.lift (evBil F W)

★★自然性は **`φ` の自然性そのもの**である
——`Over.homMk (homOfLE h) : Over.mk (homOfLE h) ⟶ Over.mk (𝟙 W)` に沿って取り、
`fVal x` を代入するだけ。

## ★★★逃げ道 2 つ

| 症状 | 逃げ道 |
|---|---|
| `rw [ConcreteCategory.comp_apply]` が型検査で落ちる | ★**`erw`** |
| 残りの等式 | ★★**`congrArg … hnat` 一発**(すべて defeq) |

★★★★第 173 の両側の橋(`fVal`)のおかげで、
`(F ⊗ F^∨).map` と `φ.app` の合成が**定義から一致**した。
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace

variable {X : Scheme.{u}} (F : X.PresheafOfModules)

set_option maxHeartbeats 1000000 in
/-- ★★評価射。 -/
noncomputable def evHom :
    F ⊗ (dualPresheaf F) ⟶ 𝟙_ (X.PresheafOfModules) where
  app W := ModuleCat.ofHom (TensorProduct.lift (evBil F W))
  naturality := by
    intro W W' f
    refine ModuleCat.MonoidalCategory.tensor_ext ?_
    intro x φ
    erw [ConcreteCategory.comp_apply, ConcreteCategory.comp_apply]
    have hnat := φ.naturality
      (Over.homMk (homOfLE (leOfHom f.unop))
        : (Over.mk (homOfLE (leOfHom f.unop)) : Over W.unop) ⟶ Over.mk (𝟙 W.unop)).op
    exact congrArg (fun (m : _ ⟶ _) => (ModuleCat.Hom.hom m) (fVal F W x)) hnat


/-! ## ★出典の紐付け(`.src`) -/

def evHom.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——評価射)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
