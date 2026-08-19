import ABC3.Found.Arakelov.ArcContMetric

/-!
# Arakelov (C3) 第 284 ブロック —— **被覆から 1 の分割が出る**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★`Spec ℂ` が 1 点であることが効く

`X` の開被覆 `U` があれば、`X^arc` の被覆 `arcOpenSet (U i)` が得られる。
★★理由は単純で、**`ℂ` は体なので `Spec ℂ` は 1 点**(mathlib の
`instance [Field R] : Unique (PrimeSpectrum R)`)——
したがって `p` の像が `U i` に入れば `p ⁻¹ᵁ (U i) = ⊤` である。

## ★★コンパクト Hausdorff から 1 の分割へ

| 要求 | mathlib |
|---|---|
| `NormalSpace` | ★コンパクト + T2 から instance で出る |
| `ParacompactSpace` | ★コンパクトから instance で出る |
| 1 の分割の存在 | ★★`PartitionOfUnity.exists_isSubordinate` |

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `preimage_eq_top_of_mem` | ★像が入れば逆像は全体 |
| `arcOpenSet_cover` | ★★`X` の被覆から `X^arc` の被覆 |
| `exists_pou_of_cover` | ★★★★**1 の分割の存在** |
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace

variable {X : Scheme.{0}} {ι : Type}

/-- ★**像がすべて `W` に入るなら逆像は全体である**。 -/
theorem preimage_eq_top_of_mem (p : Spec (CommRingCat.of ℂ) ⟶ X) (W : X.Opens)
    (h : ∀ y, p.base y ∈ W) : p ⁻¹ᵁ W = ⊤ := by
  ext y
  exact ⟨fun _ => trivial, fun _ => h y⟩

/-- ★★**`X` の開被覆から `X^arc` の被覆が出る**——`Spec ℂ` が 1 点だから。 -/
theorem arcOpenSet_cover (U : ι → X.Opens) (hU : ∀ x : X, ∃ i, x ∈ U i) :
    (Set.univ : Set (Spec (CommRingCat.of ℂ) ⟶ X)) ⊆ ⋃ i, arcOpenSet (U i) := by
  intro p _
  haveI : Unique (Spec (CommRingCat.of ℂ) : Type) := inferInstanceAs (Unique (PrimeSpectrum ℂ))
  obtain ⟨i, hi⟩ := hU (p.base default)
  refine Set.mem_iUnion.2 ⟨i, ?_⟩
  refine preimage_eq_top_of_mem p (U i) (fun y => ?_)
  rw [Subsingleton.elim y default]
  exact hi

/-- ★★★★**1 の分割の存在**——`X^arc` がコンパクト Hausdorff なら出る。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これが (C3) の `metric_nonempty` に要る唯一の点集合位相の入力である。 -/
theorem exists_pou_of_cover (U : ι → X.Opens) (hU : ∀ x : X, ∃ i, x ∈ U i)
    (hc : @CompactSpace (Spec (CommRingCat.of ℂ) ⟶ X) (arcTopology X))
    (ht : @T2Space (Spec (CommRingCat.of ℂ) ⟶ X) (arcTopology X)) :
    ∃ f : @PartitionOfUnity ι (Spec (CommRingCat.of ℂ) ⟶ X) (arcTopology X) Set.univ,
      @PartitionOfUnity.IsSubordinate ι (Spec (CommRingCat.of ℂ) ⟶ X) (arcTopology X)
        Set.univ f (fun i => arcOpenSet (U i)) := by
  letI := arcTopology X
  haveI := hc
  haveI := ht
  exact PartitionOfUnity.exists_isSubordinate isClosed_univ (fun i => arcOpenSet (U i))
    (fun i => isOpen_arcOpenSet (U i)) (arcOpenSet_cover U hU)

/-! ## ★出典の紐付け(`.src`) -/

def exists_pou_of_cover.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C3——被覆から 1 の分割へ)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
