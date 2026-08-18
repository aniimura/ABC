import ABC3.Found.Arakelov.PicImgClosure
import ABC3.Found.Arakelov.PicComapImage

/-!
# Arakelov (B2) 第 213 ブロック —— **可換正方形を切断で読む**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★§9-238 で特定した B2 最後の 1 点

    Spec.map (f.appLE B A i) ≫ hB.fromSpec = hA.fromSpec ≫ f      ★mathlib

を **`appTop` に落として両辺を分解する**:

| 辺 | 分解 |
|---|---|
| 左 | `hB.fromSpec.appTop ≫ (Spec.map appLE).appTop` |
| 右 | `f.appTop ≫ hA.fromSpec.appTop` |

★★これで「切断を `Spec` 側へ運んでから `appLE` を当てる」と
「先に引き戻してから `Spec` 側へ運ぶ」が**一致する**と言える。

## ★★`appTop` は関手的

`Scheme.Hom.comp_appTop` が両辺の分解を与える。★射の等式に `appTop` を当てる
(`congrArg`)だけで、切断レベルの両立が出る——**新しい数学は要らなかった**。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `appLE_fromSpec_sections` | ★★★可換正方形の `appTop` 版 |
| `lhs_decomp` / `rhs_decomp` | ★両辺の分解 |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite TopologicalSpace

variable {X Y : Scheme.{u}} (f : X ⟶ Y)

/-- ★★★**可換正方形を `appTop` で読む**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★これが B2 最後の 1 点である。 -/
theorem appLE_fromSpec_sections {A : X.Opens} {B : Y.Opens}
    (hA : IsAffineOpen A) (hB : IsAffineOpen B) (i : A ≤ f ⁻¹ᵁ B) :
    (Spec.map (f.appLE B A i) ≫ hB.fromSpec).appTop
      = (hA.fromSpec ≫ f).appTop :=
  congrArg (fun (m : _ ⟶ _) => Scheme.Hom.appTop m)
    (IsAffineOpen.SpecMap_appLE_fromSpec f hB hA i)

/-- ★**左辺の分解**。 -/
theorem lhs_decomp {A : X.Opens} {B : Y.Opens}
    (hB : IsAffineOpen B) (i : A ≤ f ⁻¹ᵁ B) :
    (Spec.map (f.appLE B A i) ≫ hB.fromSpec).appTop
      = hB.fromSpec.appTop ≫ (Spec.map (f.appLE B A i)).appTop :=
  Scheme.Hom.comp_appTop _ _

/-- ★**右辺の分解**。 -/
theorem rhs_decomp {A : X.Opens} (hA : IsAffineOpen A) :
    (hA.fromSpec ≫ f).appTop = f.appTop ≫ hA.fromSpec.appTop :=
  Scheme.Hom.comp_appTop _ _

/-! ## ★出典の紐付け(`.src`) -/

def appLE_fromSpec_sections.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——可換正方形を切断で読む)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
