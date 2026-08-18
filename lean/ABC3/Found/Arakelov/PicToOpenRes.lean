import ABC3.Found.Arakelov.PicBasicCover

/-!
# Arakelov (B1) 第 118 ブロック —— **切断の制限は局所化と両立する**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★生成元の制限が生成元であることの土台

`M_g ≅ R_g` の生成元を `D(g)` の切断にしたものを `D(h·g)` へ制限すると、
`M_{h·g}` の生成元になる——これが局所自明性の最後の一点である。

★★その土台は **`tilde.toOpen` が制限と両立する**ことであり、
mathlib の `tilde.toOpen_res` が**`rfl`** で与える(2026-08-24 実測)。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `toOpen_res'` | ★`toOpen` は制限と両立(mathlib の在庫に名を付ける) |

## ★★★残り

★`tildeAwayEquiv`(第 85)どうしが制限で対応することを、
`IsLocalizedModule` の一意性(第 86 の `linearMap_ext` と同型の議論)で言えば、
生成元の制限が生成元であることが出る。
★★そこから第 112(局所化は全単射を保つ)+ 第 103(生成元の乗法は全単射)で
第 115 の仮定が埋まり、`IsLocallyTrivial` が出る。
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite TopologicalSpace

variable (R : CommRingCat.{u}) (M : ModuleCat.{u} (R : Type u))

/-- ★**`toOpen` は制限と両立する**——mathlib の在庫に名を付ける(`rfl`)。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★これが「生成元の制限は生成元」の土台である。 -/
theorem toOpen_res' (U W : (Spec R).Opens) (i : W ⟶ U) :
    tilde.toOpen M U ≫ (modulesSpecToSheaf.obj (tilde M)).presheaf.map i.op
      = tilde.toOpen M W := tilde.toOpen_res _ _ _ _

/-! ## ★出典の紐付け(`.src`) -/

def toOpen_res'.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——切断の制限は局所化と両立する)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
