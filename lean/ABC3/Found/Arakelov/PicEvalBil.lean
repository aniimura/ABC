import ABC3.Found.Arakelov.PicFVal

/-!
# Arakelov (B2) 第 174 ブロック —— ★★★★★★**評価の双線型写像**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★逆層の同型 `F ⊗ F^∨ ≅ 𝟙_` の第 1 歩

    evBil W : F(W) →ₗ F^∨(W) →ₗ 𝒪(W)
      x ↦ φ ↦ φ.app t x

★4 つの線型性がすべて必要である:

| 欄 | 根拠 |
|---|---|
| `φ` について加法的 | ★`rfl` |
| `φ` について線型 | ★第 172(`unitMul_app_apply`)+ `mul_comm` |
| `x` について加法的 | ★第 173(`fVal_add`) |
| `x` について線型 | ★第 173(`fVal_smul`)+ `unitVal_smul` |

★★第 173 の**両側の橋**が無ければ書けなかった。

★★★`maxHeartbeats 1000000` が要る——`LinearMap.ext` の単一化が重い。
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace

variable {X : Scheme.{u}} (F : X.PresheafOfModules)

set_option maxHeartbeats 1000000 in
/-- ★★評価の双線型写像。 -/
noncomputable def evBil (W : (X.Opens)ᵒᵖ) :
    (F.obj W : Type u) →ₗ[((X.presheaf.obj W) : Type u)]
      ((dualPresheaf F).obj W : Type u) →ₗ[((X.presheaf.obj W) : Type u)]
        ((𝟙_ (X.PresheafOfModules)).obj W : Type u) where
  toFun x :=
    { toFun := fun φ => unitVal W.unop ((φ.app (op (Over.mk (𝟙 W.unop)))).hom (fVal F W x))
      map_add' := fun a b => rfl
      map_smul' := fun c φ => by
        have h1 : (c • φ) = φ ≫ unitMul W.unop c := dual_smul_eq F W.unop c φ
        rw [h1, PresheafOfModules.comp_app, ConcreteCategory.comp_apply]
        exact (unitMul_app_apply W.unop c _).trans (mul_comm _ _) }
  map_add' := fun a b => by
    refine LinearMap.ext fun φ => ?_
    show unitVal W.unop ((φ.app (op (Over.mk (𝟙 W.unop)))).hom (fVal F W (a + b))) = _
    rw [fVal_add, map_add]
    rfl
  map_smul' := fun c x => by
    refine LinearMap.ext fun φ => ?_
    show unitVal W.unop ((φ.app (op (Over.mk (𝟙 W.unop)))).hom (fVal F W (c • x))) = _
    rw [fVal_smul, map_smul, unitVal_smul]
    rfl


/-! ## ★出典の紐付け(`.src`) -/

def evBil.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——評価の双線型写像)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
