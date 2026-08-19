import ABC3.Found.Arakelov.ArcLogMetric

/-!
# Arakelov (C3) 第 291 ブロック —— **1 の分割は引き戻せる**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★§9-335 の案 (b) の第 1 歩

「貼り合わせ計量の `V` への制限 = 制限したデータの貼り合わせ計量」を示すために、
まず **1 の分割を `V^arc` へ引き戻す**。

## ★★測った——開集合の引き戻しは `rfl`

| 問い | 結果 |
|---|---|
| `(q ≫ V.ι) ⁻¹ᵁ W = q ⁻¹ᵁ (V.ι ⁻¹ᵁ W)` | ★`rfl` |
| `(· ≫ V.ι) ⁻¹' arcOpenSet W = arcOpenSet (V.ι ⁻¹ᵁ W)` | ★★`rfl` |

★★★弧空間の開集合の対応が**定義から**出るので、従属性の移送が軽い。

## ★★本ブロックで取れるもの

| 定義・定理 | 内容 |
|---|---|
| `LocallyFinite.preimage_of_continuous` | ★局所有限性は連続写像で引き戻せる |
| `pouPullback` | ★★★**1 の分割の引き戻し** |
| `pouPullback_isSubordinate` | ★★★従属性も引き戻せる |
| `arcOpenSet_pullback` | ★弧開集合の引き戻し(`rfl`) |
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace

/-- ★**局所有限性は連続写像で引き戻せる**。 -/
theorem LocallyFinite.preimage_of_continuous {ι Y Z : Type} [TopologicalSpace Y]
    [TopologicalSpace Z] {s : ι → Set Z} (hs : LocallyFinite s) {g : Y → Z} (hg : Continuous g) :
    LocallyFinite (fun i => g ⁻¹' (s i)) := by
  intro y
  obtain ⟨t, ht, hfin⟩ := hs (g y)
  refine ⟨g ⁻¹' t, hg.continuousAt.preimage_mem_nhds ht, hfin.subset (fun i hi => ?_)⟩
  obtain ⟨z, hz1, hz2⟩ := hi
  exact ⟨g z, hz1, hz2⟩

/-- ★★★**1 の分割の引き戻し**。 -/
noncomputable def pouPullback {ι Y Z : Type} [TopologicalSpace Y] [TopologicalSpace Z]
    (f : PartitionOfUnity ι Z Set.univ) (g : C(Y, Z)) : PartitionOfUnity ι Y Set.univ where
  toFun := fun i => (f i).comp g
  locallyFinite' := LocallyFinite.preimage_of_continuous f.locallyFinite g.continuous
  nonneg' := fun i y => f.nonneg i (g y)
  sum_eq_one' := fun y _ => f.sum_eq_one (Set.mem_univ (g y))
  sum_le_one' := fun y => f.sum_le_one' (g y)

@[simp] theorem pouPullback_apply {ι Y Z : Type} [TopologicalSpace Y] [TopologicalSpace Z]
    (f : PartitionOfUnity ι Z Set.univ) (g : C(Y, Z)) (i : ι) (y : Y) :
    pouPullback f g i y = f i (g y) := rfl

/-- ★★★**従属性も引き戻せる**。 -/
theorem pouPullback_isSubordinate {ι Y Z : Type} [TopologicalSpace Y] [TopologicalSpace Z]
    (f : PartitionOfUnity ι Z Set.univ) (g : C(Y, Z)) (O : ι → Set Z)
    (h : f.IsSubordinate O) : (pouPullback f g).IsSubordinate (fun i => g ⁻¹' (O i)) := by
  intro i
  refine subset_trans ?_ (Set.preimage_mono (h i))
  exact subset_trans (closure_mono (subset_of_eq rfl))
    (g.continuous.closure_preimage_subset (Function.support (f i)))

variable {X : Scheme.{0}} (V : X.Opens)

/-- ★**弧開集合の引き戻しは開集合の引き戻し**(`rfl`)。 -/
theorem arcOpenSet_pullback (W : X.Opens) :
    (fun q : Spec (CommRingCat.of ℂ) ⟶ V.toScheme => q ≫ V.ι) ⁻¹' (arcOpenSet W)
      = arcOpenSet ((Opens.map V.ι.base).obj W) := rfl

/-! ## ★出典の紐付け(`.src`) -/

def pouPullback.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C3——1 の分割の引き戻し)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
