import ABC3.Found.Arakelov.PicTildeInv

/-!
# Arakelov (B1) 第 93 ブロック —— **切断は単位からの射を与える**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★局所自明性に要る器具

`IsLocallyTrivial` は `(restrict V).obj P ≅ 𝟙_` を要求する。
★★その同型を作るには、まず**射**を作らねばならない。

★★★**単位からの射は `V` での切断そのものである**:

    Hom(𝟙_, P)  ≅  Hom(free (yoneda (Over.mk (𝟙 V))), P)  ≅  P(V)

| 段 | 出典 |
|---|---|
| `𝟙_ ≅ free (yoneda 終対象)` | ★第 55 ブロック(一般形) |
| `Hom(free (yoneda X), P) ≅ P(X)` | ★mathlib `freeYonedaEquiv` |

## ★★本ブロックで取れるもの

| 定義 | 内容 |
|---|---|
| `unitHomOfSection` | ★★★★**切断から `𝟙_ ⟶ P`** |

## ★★★次

★可逆加群 `M` は `D(r)` の上で `M_r ≅ R_r`(第 76 ブロック)なので、
生成元に対応する切断が取れる。★★それが与える `𝟙_ ⟶ (tilde M)|_{D r}` が
**同型**であることを示せば局所自明性が出る。
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace

variable {X : Scheme.{u}} (V : X.Opens) (P : PresheafModulesOn X V)

/-- ★★★★**切断は単位からの射を与える**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★`Over V` の終対象での切断 `s` に対し、`c ↦ c · s` を与える射である。 -/
noncomputable def unitHomOfSection (s : P.obj (op (Over.mk (𝟙 V)))) :
    𝟙_ (PresheafModulesOn X V) ⟶ P :=
  (freeYonedaTermIso (R := (Over.forget V).op ⋙ X.presheaf) (Over.mk (𝟙 V))
      (overTerminalUnique V)).inv
    ≫ PresheafOfModules.freeYonedaEquiv.symm s

/-! ## ★出典の紐付け(`.src`) -/

def unitHomOfSection.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——切断は単位からの射を与える)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
