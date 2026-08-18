import ABC3.Found.Arakelov.PicDualSmul

/-!
# Arakelov (B2) 第 168 ブロック —— ★★★★★★**制限の半線型性**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★双対を前層に組む最後の部品

    (restrictOnFunctor h).map (unitMul U c) = unitMul V (c|_V)

★これが「双対の制限が半線型である」ことの中身である。

## ★★逃げ道——**`unitMul` を展開しない**

`unitMul U c` の各成分の具体形(`x ↦ x · c|_W`)を出そうとすると
`freeYonedaEquiv.symm` の展開が要って重い。

★★**自然性を使う**と展開が要らない:

    (unitMul U c) の自然性を Over.homMk (homOfLE h) に沿って取り、
    unitOne U を代入する

と、左辺が `(unitMul U c).app (Over.mk (homOfLE h)) (unitOne V)`(`res 1 = 1`)、
右辺が `(unitEnd U (unitMul U c))|_V = c|_V`(★第 165)になる。

★★★**自然性 1 本で済んだ。**
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace

variable {X : Scheme.{u}} {U V : X.Opens}

set_option maxHeartbeats 1000000 in
theorem unitMul_res (h : V ≤ U) (c : (Γ(X, U) : Type u)) :
    (restrictOnFunctor h).map (unitMul U c)
      = unitMul V ((X.presheaf.map (homOfLE h).op).hom c) := by
  apply (unitEndEquiv V).injective
  rw [show (unitEndEquiv V) (unitMul V ((X.presheaf.map (homOfLE h).op).hom c))
      = (X.presheaf.map (homOfLE h).op).hom c from unitEnd_unitMul V _]
  have hnat := (unitMul U c).naturality
    (Over.homMk (homOfLE h) : (Over.mk (homOfLE h) : Over U) ⟶ Over.mk (𝟙 U)).op
  have hone : ((𝟙_ (PresheafModulesOn X U)).map
      (Over.homMk (homOfLE h) : (Over.mk (homOfLE h) : Over U) ⟶ Over.mk (𝟙 U)).op).hom
      (unitOne U) = unitOne V := by
    show (X.presheaf.map _).hom (1 : (Γ(X, U) : Type u)) = (1 : (Γ(X, V) : Type u))
    exact map_one _
  have happ := congrArg (fun (m : _ ⟶ _) => (ModuleCat.Hom.hom m) (unitOne U)) hnat
  simp only [ConcreteCategory.comp_apply] at happ
  rw [hone] at happ
  show unitVal V (((unitMul U c).app (op (Over.mk (homOfLE h)))).hom (unitOne V)) = _
  exact happ.trans (congrArg ((X.presheaf.map (homOfLE h).op).hom) (unitEnd_unitMul U c))


/-! ## ★出典の紐付け(`.src`) -/

def unitMul_res.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——制限の半線型性)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
