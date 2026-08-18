import ABC3.Found.Arakelov.PicDualSheaf
import ABC3.Found.Arakelov.PicDualTrivial
import ABC3.Found.Arakelov.PicSheafTensor

/-!
# Arakelov (B2) 第 180 ブロック —— **双対は層加群であり局所自明**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★`InvSheaf` の 5 欄のうち 3 欄が埋まる

    structure InvSheaf (X : Scheme) where
      carrier : X.Modules       -- ← idealSheaf D(第 151)
      inv     : X.Modules       -- ← ★本ブロック(dualModules)
      isInv   : ...             -- ← 評価射の同型(次)
      trivial : ...             -- ← 第 162
      invTrivial : ...          -- ← ★本ブロック

## ★★局所自明性は「同じ篩を使い回す」だけ

`dualTrivialIso`(第 171)は `(F|_V ≅ 𝟙_) ⟹ (𝟙_ ≅ (dual F)|_V)` を与える。
★したがって `F` の自明化篩を**そのまま**使えばよい。

## ★★本ブロックで取れるもの

| 定義・定理 | 内容 |
|---|---|
| `isLocallyTrivial_dualPresheaf` | ★★★**双対も局所自明** |
| `dualModules` | ★★★**双対の層加群** |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace

variable {X : Scheme.{u}}

/-- ★★★**双対も局所自明である**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★自明化篩は `F` のものをそのまま使う——第 171 の `dualTrivialIso` が
各葉で逆向きの同型を与えるからである。 -/
theorem isLocallyTrivial_dualPresheaf (F : X.PresheafOfModules) (hF : IsLocallyTrivial X F) :
    IsLocallyTrivial X (dualPresheaf F) := by
  intro U
  obtain ⟨S, hS, hiso⟩ := hF U
  refine ⟨S, hS, fun V i hi => ?_⟩
  obtain ⟨e⟩ := hiso i hi
  exact ⟨(dualTrivialIso F V e).symm⟩

/-- ★★★**双対の層加群**——第 179 で層であることを示したので `X.Modules` に載る。 -/
noncomputable def dualModules (F : X.Modules) : X.Modules where
  val := dualPresheaf F.val
  isSheaf := isSheaf_dualPresheaf F.val

@[simp] theorem dualModules_val (F : X.Modules) : (dualModules F).val = dualPresheaf F.val := rfl

/-! ## ★出典の紐付け(`.src`) -/

def dualModules.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——双対の層加群)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
