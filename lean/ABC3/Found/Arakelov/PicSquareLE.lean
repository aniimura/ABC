import ABC3.Found.Arakelov.PicCongrApp

/-!
# Arakelov (B2) 第 226 ブロック —— ★★★★**可換正方形を任意の座標で読む**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★§9-256 の「座標系 `⊤` を捨てる」が完成した

    (Spec.map appLE ≫ hB.fromSpec).appLE B W e₁ = (hA.fromSpec ≫ f).appLE B W e₂

★**`W` は任意**である——第 213 は `appTop`(`⊤` 固定)だったが、
本ブロックは**どの開集合でも**読める。

## ★★機構は `appLE` の**終域が固定**であること

`app U` の終域は `Γ(X, f⁻¹U)` で射に依存するが、
★`appLE U V e` の終域は **`Γ(X, V)`** で**射に依存しない**。
★★したがって等式が**素直に書ける**——依存位置が消える。

★★★摩擦 #2(始域が違う)の逃げ道「分解を挟む」と同じ発想である
——**依存しない座標を選ぶ**。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `square_appLE` | ★★★★**可換正方形を任意の座標で読む** |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite TopologicalSpace

variable {X Y : Scheme.{u}} (f : X ⟶ Y)

/-- ★★★★**可換正方形を任意の座標 `W` で読む**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★`appLE` の終域は射に依存しないので、依存位置が消える。 -/
theorem square_appLE {A : X.Opens} {B : Y.Opens}
    (hA : IsAffineOpen A) (hB : IsAffineOpen B) (i : A ≤ f ⁻¹ᵁ B)
    (W : (Spec Γ(X, A)).Opens)
    (e1 : W ≤ (Spec.map (f.appLE B A i) ≫ hB.fromSpec) ⁻¹ᵁ B)
    (e2 : W ≤ (hA.fromSpec ≫ f) ⁻¹ᵁ B) :
    (Spec.map (f.appLE B A i) ≫ hB.fromSpec).appLE B W e1
      = (hA.fromSpec ≫ f).appLE B W e2 := by
  show Scheme.Hom.app _ B ≫ _ = Scheme.Hom.app _ B ≫ _
  rw [Scheme.Hom.congr_app (IsAffineOpen.SpecMap_appLE_fromSpec f hB hA i) B,
    Category.assoc, ← Functor.map_comp, ← op_comp]
  congr 1

/-! ## ★出典の紐付け(`.src`) -/

def square_appLE.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——可換正方形を任意の座標で読む)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
