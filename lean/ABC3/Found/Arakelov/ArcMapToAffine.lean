import ABC3.Found.Arakelov.ArcMapAffine
import ABC3.Found.Arakelov.ArcContCriterion

/-!
# Arakelov (D1) 第 298 ブロック —— **アフィンへの射は弧空間で連続**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★弧空間の関手性——第 2 歩

第 297 は `Spec A ⟶ Spec B` の場合だった。本ブロックでは**始域を一般のスキームに広げる**:

    g : X ⟶ Spec B  ⟹  (· ≫ g) : X^arc ⟶ (Spec B)^arc  は連続

★★第 279 の `continuous_of_charts`(「各アフィン chart で連続なら連続」)で
アフィン chart に落とし、第 297 を当てるだけである。

★`arcTopologyOpen U = induced (· ≫ isoSpec.hom) (arcTopologyAffine Γ(X,U))` なので、
chart の上では `isoSpec` を挟んで第 297 の形になる。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `continuous_comp_toAffine` | ★★★★**アフィンへの射は弧空間で連続** |
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite TopologicalSpace

/-- ★★★★**アフィンへの射は弧空間で連続である**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★各アフィン chart で第 297 に落ちる。 -/
theorem continuous_comp_toAffine {X : Scheme.{0}} {B : CommRingCat.{0}} (g : X ⟶ Spec B) :
    @Continuous _ _ (arcTopology X) (arcTopologyAffine B)
      (fun p : Spec (CommRingCat.of ℂ) ⟶ X => p ≫ g) := by
  letI := arcTopologyAffine B
  refine continuous_of_charts (fun p => p ≫ g) (fun U => ?_)
  letI := arcTopologyOpen U
  letI := arcTopologyAffine (Γ(X, U.1))
  have he : (fun p : Spec (CommRingCat.of ℂ) ⟶ U.1.toScheme => (p ≫ U.1.ι) ≫ g)
      = (fun r : Spec (CommRingCat.of ℂ) ⟶ Spec (Γ(X, U.1)) =>
          r ≫ (U.2.isoSpec.inv ≫ U.1.ι ≫ g))
        ∘ (fun p : Spec (CommRingCat.of ℂ) ⟶ U.1.toScheme => p ≫ U.2.isoSpec.hom) := by
    funext p
    show (p ≫ U.1.ι) ≫ g = (p ≫ U.2.isoSpec.hom) ≫ (U.2.isoSpec.inv ≫ U.1.ι ≫ g)
    simp only [Category.assoc, Iso.hom_inv_id_assoc]
  rw [he]
  exact (continuous_comp_affine (U.2.isoSpec.inv ≫ U.1.ι ≫ g)).comp continuous_induced_dom

/-! ## ★出典の紐付け(`.src`) -/

def continuous_comp_toAffine.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C——アフィンへの射が弧空間で連続であること)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
