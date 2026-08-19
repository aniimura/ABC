import ABC3.Found.Arakelov.ArcGlueCont
import ABC3.Found.Arakelov.ArcOpenMap
import ABC3.Found.Arakelov.ArcLiftOpen

/-!
# Arakelov (C3) 第 274 ブロック —— ★★★**`V` を通る点の集合は開**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★第 273 の仮定 `ho` を落とす

    arcOpenSet V = {p | p ⁻¹ᵁ V = ⊤} = Set.range (· ≫ V.ι)

★等式は第 259 の 2 つ(`liftToOpen_fac` と `comp_preimage_eq_top`)で出る。
★★開性は **C1 の `isOpen_range_comp_ι`**(第 30 ブロック台)がそのまま効く。

★★★**C1 で作った器具が、C3 の貼り合わせでそのまま使えた**——
`X^arc` の位相を「アフィン chart の `⨆`」で定義したことの配当である。

| 定理 | 内容 |
|---|---|
| `arcOpenSet_eq_range` | ★像としての表示 |
| `isOpen_arcOpenSet` | ★★★開性 |
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite TopologicalSpace

variable {X : Scheme.{0}} (V : X.Opens)

/-- ★`V` を通る点の集合は `V.ι` の像である。 -/
theorem arcOpenSet_eq_range :
    arcOpenSet V = Set.range (fun q : Spec (CommRingCat.of ℂ) ⟶ V.toScheme => q ≫ V.ι) := by
  ext p
  constructor
  · intro h
    exact ⟨liftToOpenOfTop V p h, liftToOpen_fac V p h⟩
  · rintro ⟨q, rfl⟩
    exact comp_preimage_eq_top V q

/-- ★★★`V` を通る点の集合は開である。 -/
theorem isOpen_arcOpenSet : @IsOpen _ (arcTopology X) (arcOpenSet V) := by
  rw [arcOpenSet_eq_range]
  exact isOpen_range_comp_ι V


/-! ## ★出典の紐付け(`.src`) -/

def isOpen_arcOpenSet.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C——V を通る点の集合は開)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
