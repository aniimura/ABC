import ABC3.Found.Arakelov.PicAppIsoInv

/-!
# Arakelov (B2) 第 225 ブロック —— ★★★**`congr_app` が依存位置を解く**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★摩擦 #4(依存位置)の**別の逃げ道**

`f.app U` の**終域が `f` に依存する**(`Γ(X, f ⁻¹ᵁ U)`)ので、
射の等式 `f = g` から `f.app U = g.app U` を `congrArg` で作ろうとすると型が合わない。

★★mathlib の **`Scheme.Hom.congr_app`** がこれを解く:

    congr_app (e : f = g) (U) : f.app U = g.app U ≫ X.presheaf.map (eqToHom …).op

★**`eqToHom` を右に押し出す**形で等式を与える——依存位置を「転送つき」に変換している。

## ★★これで `⊤` を使わない座標系が組める

第 213 は `appTop`(`⊤` 固定)で書いたが、本ブロックの形なら**任意の `U`** で書ける。
★§9-256 で「座標系 `⊤` を捨てる」と決めた筋が通る。

## ★★摩擦 #4 の逃げ道が 2 つになった

| 逃げ道 | 場面 |
|---|---|
| mathlib を grep し直す | 第 222(`app_appIso_inv` が在った) |
| ★**`congr_app` で `eqToHom` に変換** | ★本ブロック |

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `appLE_fromSpec_pre` | ★★★可換正方形を任意の `B` で読む |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite TopologicalSpace

variable {X Y : Scheme.{u}} (f : X ⟶ Y)

/-- ★★★**可換正方形を任意の `B` で読む**(`congr_app` 版)。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★`app` の終域が射に依存するので、`congrArg` ではなく `congr_app` を使う。 -/
theorem appLE_fromSpec_pre {A : X.Opens} {B : Y.Opens}
    (hA : IsAffineOpen A) (hB : IsAffineOpen B) (i : A ≤ f ⁻¹ᵁ B) :
    (Spec.map (f.appLE B A i) ≫ hB.fromSpec).app B
      = (hA.fromSpec ≫ f).app B ≫ (Spec Γ(X, A)).presheaf.map (eqToHom (by
          rw [IsAffineOpen.SpecMap_appLE_fromSpec f hB hA i])).op :=
  Scheme.Hom.congr_app (IsAffineOpen.SpecMap_appLE_fromSpec f hB hA i) B

/-! ## ★出典の紐付け(`.src`) -/

def appLE_fromSpec_pre.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——congr_app が依存位置を解く)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
