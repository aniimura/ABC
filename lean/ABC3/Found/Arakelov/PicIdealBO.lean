import ABC3.Found.Arakelov.PicIdealSheaf

/-!
# Arakelov (B2) 第 152 ブロック —— **基本開集合での切断**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★局所自明性へ向けた第 1 歩

    idealSections D (X.basicOpen f) = (D.ideal A).map (制限)

★アフィン開 `A` の基本開集合では、切断は元のイデアルの**像**である。

★★これで「`D.ideal A` が可逆 ⟹ `idealSections` は局所自由」の議論が
**第 130 ブロック(`bijective_smul_liftGen`)を再利用して**書ける
——`Γ(tilde M, ·)` のときと同じ形だからである。

## ★★逃げ道——`rw` ではなく `Eq.trans`

`rw [D.map_ideal …]` は型検査(instances 透明度)で落ちる。
★`(idealSections_affine …).trans (D.map_ideal …).symm` と**項で書けば**通る
——[[exact-term-over-rw]] の 5 例目。
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite TopologicalSpace

variable {X : Scheme.{u}} (D : X.IdealSheafData)

/-- ★基本開集合での切断は元のイデアルの像である。 -/
theorem idealSections_basicOpen {A : X.affineOpens} (f : (Γ(X, A.1) : Type u)) :
    idealSections D (X.basicOpen f)
      = (D.ideal A).map (X.presheaf.map (homOfLE (X.basicOpen_le f)).op).hom := by
  exact (idealSections_affine D (X.affineBasicOpen f)).trans
    (D.map_ideal (U := X.affineBasicOpen f) (V := A) (X.basicOpen_le f)).symm


/-! ## ★出典の紐付け(`.src`) -/

def idealSections_basicOpen.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——基本開集合での切断)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
