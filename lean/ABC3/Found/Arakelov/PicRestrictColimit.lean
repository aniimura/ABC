import ABC3.Found.Arakelov.PicFreeTop

/-!
# Arakelov (B1) 第 43 ブロック —— **制限は余極限を保つ**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★Beck–Chevalley の第 1 段

§9-52 で測った通り、第 41 ブロックの仮定を示すには

    (f^*_pre L)|_{f⁻¹V} ≅ (f|)^*_pre (L|_V)

が要る(mathlib に無い)。★★同型性は**第 28・29 ブロックの器具**で持ち上がるが、
その前提として**両辺が余極限を保つ**ことが要る。

| 側 | 余極限を保つか | 根拠 |
|---|---|---|
| `f^*_pre` | ★保つ | 左随伴(第 23 ブロック) |
| 制限 | ★★保つ | **前層の余極限は切断ごと**(本ブロック) |

## ★★機構

`restrictPresheafFunctor X U` は `Over.forget U` に沿った precomposition である。
★前層加群の余極限は切断ごとに計算されるので、precomposition は**何も壊さない**。

★★mathlib の `evaluationJointlyReflectsColimits`(余極限は切断ごとに判定できる)
と `evaluation` の余極限保存を組み合わせるだけである。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `restrictPresheafFunctor_preservesColimits` | ★★★**制限は余極限を保つ** |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace Limits

variable (X : Scheme.{u}) (U : X.Opens)

/-- ★★★**制限は余極限を保つ**——前層の余極限は切断ごとだからである。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★これが Beck–Chevalley(`(f^*_pre L)|_{f⁻¹V} ≅ (f|)^*_pre (L|_V)`)を
第 29 ブロックの器具で持ち上げるための前提である。 -/
noncomputable def restrictPresheafFunctor_preservesColimits :
    PreservesColimitsOfSize.{u, u} (restrictPresheafFunctor X U) where
  preservesColimitsOfShape := ⟨⟨fun hc =>
    ⟨PresheafOfModules.evaluationJointlyReflectsColimits _ _
      (fun Z => isColimitOfPreserves
        (PresheafOfModules.evaluation _ ((Over.forget U).op.obj Z)) hc)⟩⟩⟩

/-! ## ★出典の紐付け(`.src`) -/

def restrictPresheafFunctor_preservesColimits.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——制限が余極限を保つこと)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
