import ABC3.Found.Arakelov.ArcGreenMetric
import ABC3.Found.GenEll.ArchConj

/-!
# Arakelov (C3) 第 293 ブロック —— **`ι_X` との両立**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★Green 表示だと「共役両立」もそのまま出る

Interface の

    IsConjCompatible X L m ↔ ∀ p, logMetric X L m (conj X p) = logMetric X L m p

は、`logMetric = green` なので **`IsConjInvariant m.green` そのもの**である
(`Found/GenEll/ArchConj.lean` に既に在る述語)。

★★★`isConjCompatible_iff` は **`Iff.rfl`** で出る。

## ★★本ブロックで取れるもの

| 定義・定理 | 内容 |
|---|---|
| `GreenMetric.IsConjCompatible` | ★共役と両立する計量 |
| `isConjCompatible_iff` | ★★`Iff.rfl` |
| `isConjCompatible_scale` | ★★定数倍で保たれる |
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite ABC3.Found.GenEll

variable {X : Scheme.{0}} {F : X.Modules}

namespace GreenMetric

/-- ★**`ι_X` と両立する計量**——Green 関数が共役不変であること。 -/
def IsConjCompatible (m : GreenMetric X F) : Prop := IsConjInvariant m.green

/-- ★★その意味——Interface の `isConjCompatible_iff` そのもの。 -/
theorem isConjCompatible_iff (m : GreenMetric X F) :
    m.IsConjCompatible ↔ ∀ p : complexPoints X, m.green (conjPoint p) = m.green p := Iff.rfl

/-- ★★共役両立性は定数倍で保たれる。 -/
theorem isConjCompatible_scale (c : ℝ) (m : GreenMetric X F) (h : m.IsConjCompatible) :
    (scale c m).IsConjCompatible := by
  intro p
  show m.green (conjPoint p) + c = m.green p + c
  rw [h p]

/-- ★定数関数 `0` の Green 計量は共役両立である。 -/
theorem isConjCompatible_of_green_const (m : GreenMetric X F) (c : ℝ)
    (h : ∀ p, m.green p = c) : m.IsConjCompatible := by
  intro p
  rw [h (conjPoint p), h p]

end GreenMetric

/-! ## ★出典の紐付け(`.src`) -/

def GreenMetric.isConjCompatible_iff.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C3——ι_X と両立する計量)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
