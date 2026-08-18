import ABC3.Found.Arakelov.PicUnitOn

/-!
# Arakelov (B1) 第 57 ブロック —— **制限の推移律は `rfl`**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★旗を立てて測ったら杞憂だった

§9-61 で「制限の推移律は `Over A → Over B` の関手が要るので **`rfl` とは限らない**
——測ること」と旗を立てた。

★★**測った結果(2026-08-18)**: **3 つとも `rfl`** である。

| 問い | 結果 |
|---|---|
| `Over.map (homOfLE h) ⋙ Over.forget A = Over.forget B` | ★`rfl` |
| `(M|_A)|_B = M|_B` | ★★`rfl` |

★★★どちらも **precomposition の合成**だからである。

## ★★これで局所自明性の輸送が繋がる

`V ≤ f⁻¹W` のとき

    (f^*_pre L)|_V = ((f^*_pre L)|_{f⁻¹W})|_V ≅ (𝟙_)|_V = 𝟙_

★第 54–56 ブロックで `(f^*_pre L)|_{f⁻¹W} ≅ 𝟙_` が出ているので、
**制限するだけ**で `V` の上でも自明になる。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `overMap_forget` | ★site の合成(`rfl`) |
| `restrictOnFunctor` | ★`PresheafModulesOn X A ⥤ PresheafModulesOn X B` |
| `restrict_trans` | ★★★**制限の推移律**(`rfl`) |
| `restrictOnUnit` | ★★単位の制限は単位(`rfl`) |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace Limits

variable {X : Scheme.{u}} {A B : X.Opens} (h : B ≤ A)

/-- ★**site の合成**——`Over B ⥤ Over A ⥤ X.Opens` は `Over B ⥤ X.Opens` である(`rfl`)。 -/
theorem overMap_forget :
    Over.map (homOfLE h) ⋙ Over.forget A = Over.forget B := rfl

/-- ★**`Over A` 上の前層加群を `Over B` へ制限する関手**(`B ≤ A`)。 -/
noncomputable abbrev restrictOnFunctor :
    PresheafModulesOn X A ⥤ PresheafModulesOn X B :=
  PresheafOfModules.pushforward₀OfCommRingCat (Over.map (homOfLE h))
    ((Over.forget A).op ⋙ X.presheaf)

/-- ★★★**制限の推移律**——`(M|_A)|_B = M|_B`(`rfl`)。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★どちらも precomposition の合成だから**定義から等しい**。
★これで「`f⁻¹W` の上で自明」から「その下の `V` の上でも自明」が出る。 -/
theorem restrict_trans (M : X.PresheafOfModules) :
    (restrictOnFunctor h).obj ((restrictPresheafFunctor X A).obj M)
      = (restrictPresheafFunctor X B).obj M := rfl

/-- ★★**単位の制限は単位である**(`rfl`)。 -/
theorem restrictOnUnit :
    (restrictOnFunctor h).obj (𝟙_ (PresheafModulesOn X A))
      = 𝟙_ (PresheafModulesOn X B) := rfl

/-! ## ★出典の紐付け(`.src`) -/

def restrict_trans.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——制限の推移律)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
