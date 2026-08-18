import ABC3.Found.Arakelov.PicRestrictColimit

/-!
# Arakelov (B1) 第 44 ブロック —— **開集合へ制限した引き戻し**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★Beck–Chevalley の第 2 段の下ごしらえ

§9-53 の道:

    (f^*_pre L)|_{f⁻¹V} ≅ (f|)^*_pre (L|_V)

★右辺の **`(f|)^*_pre`**(制限した引き戻し)をまず作る。

## ★★機構

| 部品 | 何か | 在庫 |
|---|---|---|
| `Over.post (Opens.map f.base)` | `Over V ⥤ Over (f⁻¹V)` | ★mathlib |
| `restrictedC` | 環前層の射(`f.c` を `Over` に落とす) | ★本ブロック(2 行) |
| `pullbackPreOn` | 制限した引き戻し | ★★第 20 ブロックの `alphaR` を当てる |

★★★`Over.post` は「`W ≤ V` を `f⁻¹W ≤ f⁻¹V` に送る」関手である
——これが「`f` を `V` の上に制限したもの」の site 側の姿である。

## ★★本ブロックで取れるもの

| 定義 | 内容 |
|---|---|
| `overPost` | ★`Over V ⥤ Over (f⁻¹V)` |
| `restrictedC` | ★制限した構造層の射 |
| `pullbackPreOn` | ★★★**制限した引き戻し `(f|)^*_pre`** |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace Limits

variable {X Y : Scheme.{u}} (f : X ⟶ Y) (V : Y.Opens)

/-! ## ★site の側 -/

/-- ★**`W ≤ V` を `f⁻¹W ≤ f⁻¹V` に送る関手**——`f` を `V` の上に制限したものの site 側の姿。 -/
noncomputable abbrev overPost : Over V ⥤ Over ((Opens.map f.base).obj V) :=
  Over.post (Opens.map f.base)

/-- ★**制限した構造層の射**——`f.c` を `Over` に落としたものである。 -/
noncomputable def restrictedC :
    ((Over.forget V).op ⋙ Y.presheaf)
      ⟶ (overPost f V).op ⋙
        ((Over.forget ((Opens.map f.base).obj V)).op ⋙ X.presheaf) where
  app W := f.c.app (op (W.unop.left))
  naturality := by
    intro W W' φ
    exact f.c.naturality ((Over.forget V).map φ.unop).op

/-! ## ★★★制限した引き戻し -/

/-- ★★★**制限した引き戻し `(f|)^*_pre`**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★第 20 ブロックの `alphaR`(`CommRingCat` の射を `RingCat` へ送る)を
制限した site に当てたものである。

★★★これが Beck–Chevalley の右辺である。 -/
noncomputable abbrev pullbackPreOn :
    PresheafModulesOn Y V ⥤ PresheafModulesOn X ((Opens.map f.base).obj V) :=
  PresheafOfModules.pullback
    (alphaR (overPost f V)
      ((Over.forget ((Opens.map f.base).obj V)).op ⋙ X.presheaf)
      ((Over.forget V).op ⋙ Y.presheaf)
      (restrictedC f V))

/-! ## ★出典の紐付け(`.src`) -/

def restrictedC.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——制限した構造層の射)",
    sectionId := "genell-def-1-1-i" }

def pullbackPreOn.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——開集合へ制限した引き戻し)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
