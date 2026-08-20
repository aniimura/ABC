import ABC3.Found.Arakelov.PicUnfoldLE

/-!
# Arakelov (B2) 第 234 ブロック —— **可換正方形を `Y` 側の任意の開集合で読む**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★「座標を固定しない」を **`Y` 側にも**適用する

第 226 ブロック `square_appLE` は `X` 側の座標 `W` を自由にしたが、
`Y` 側は `B` に固定していた。★ところが残っていた 1 点は

    hB.fromSpec ''ᵁ ⊤   と   B

の**綴りの違い**(摩擦 #7)であり、`Y` 側の座標が固定されているせいで
そこに `rw` を差し込めなかった。

★★逃げ道は第 218・219・226・228 と**まったく同じ**である
——**`Y` 側の開集合 `V` もパラメータにする**。証明は 1 文字も変わらない。

★★★★これで「座標・型・仮定は、要るまで固定しない」が **5 回目**の当たりになった。

| 定理 | 内容 |
|---|---|
| `square_appLE_gen` | ★★★★**両側の座標を自由にした可換正方形** |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite TopologicalSpace

variable {X Y : Scheme.{u}} (f : X ⟶ Y)

/-- ★★★★**可換正方形を両側の任意の座標で読む**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★第 226 との違いは `V` を `B` に固定していないことだけである
——`appLE` の終域は射に依存しないので、依存位置は最初から無い。 -/
theorem square_appLE_gen {A : X.Opens} {B : Y.Opens}
    (hA : IsAffineOpen A) (hB : IsAffineOpen B) (i : A ≤ f ⁻¹ᵁ B)
    (V : Y.Opens) (W : (Spec Γ(X, A)).Opens)
    (e1 : W ≤ (Spec.map (f.appLE B A i) ≫ hB.fromSpec) ⁻¹ᵁ V)
    (e2 : W ≤ (hA.fromSpec ≫ f) ⁻¹ᵁ V) :
    (Spec.map (f.appLE B A i) ≫ hB.fromSpec).appLE V W e1
      = (hA.fromSpec ≫ f).appLE V W e2 := by
  show Scheme.Hom.app _ V ≫ _ = Scheme.Hom.app _ V ≫ _
  rw [Scheme.Hom.congr_app (IsAffineOpen.SpecMap_appLE_fromSpec f hB hA i) V,
    Category.assoc, ← Functor.map_comp, ← op_comp]
  congr 1

/-! ## ★出典の紐付け(`.src`) -/

def square_appLE_gen.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——可換正方形を両側の任意の座標で読む)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
