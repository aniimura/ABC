import ABC3.Found.Arakelov.ArcMapToAffine
import ABC3.Found.Arakelov.ArcEmbedding
import ABC3.Found.Arakelov.ArcCoverPou

/-!
# Arakelov (D1) 第 298 ブロック —— **★★★★★★弧空間は関手的である**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★任意の射で `X^arc ⟶ Y^arc` は連続

    f : X ⟶ Y  ⟹  (· ≫ f) : X^arc ⟶ Y^arc  は連続

★これが (D1) の `pullback` に要る——算術直線束の引き戻しは
**Green 関数も引き戻す**ので、連続性が保たれねばならない。

## ★★機構——局所性で終わる

| 段 | 使うもの |
|---|---|
| 1 | `Y` のアフィン開被覆から `V ∋ f(p)` を取る |
| 2 | `W := f ⁻¹ᵁ V` として `p ∈ arcOpenSet W`(第 284 `preimage_eq_top_of_mem`) |
| 3 | `arcOpenSet W = range (· ≫ W.ι)`(第 274)+ 第 275 `continuousOn_range_of_comp` |
| 4 | `W.ι ≫ f = (f ∣_ V) ≫ V.ι`(mathlib `morphismRestrict_ι`) |
| 5 | `(· ≫ (f∣_V))` は第 297、`(· ≫ V.ι)` は第 277 の開埋め込み版 |

★★★連続性は**点ごと**の性質なので、`ContinuousAt` に落として被覆で押し切れる。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `continuous_comp_scheme` | ★★★★★★**弧空間の関手性** |
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite TopologicalSpace

/-- ★★★★★★**弧空間は関手的である**——任意の射で `(· ≫ f)` は連続。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★(D1) の `pullback` が Green 関数を引き戻せるのは、これによる。 -/
theorem continuous_comp_scheme {X Y : Scheme.{0}} (f : X ⟶ Y) :
    @Continuous _ _ (arcTopology X) (arcTopology Y)
      (fun p : Spec (CommRingCat.of ℂ) ⟶ X => p ≫ f) := by
  letI := arcTopology X
  letI := arcTopology Y
  refine continuous_iff_continuousAt.2 (fun p => ?_)
  have hy : f.base (p.base default) ∈ (⊤ : Y.Opens) := trivial
  rw [← iSup_affineOpens_eq_top Y] at hy
  obtain ⟨V, hV⟩ := Opens.mem_iSup.1 hy
  have hpW : p ∈ arcOpenSet (f ⁻¹ᵁ V.1) := by
    refine preimage_eq_top_of_mem p (f ⁻¹ᵁ V.1) (fun z => ?_)
    rw [Subsingleton.elim z default]
    exact hV
  refine ContinuousOn.continuousAt ?_ ((isOpen_arcOpenSet (f ⁻¹ᵁ V.1)).mem_nhds hpW)
  rw [arcOpenSet_eq_range]
  refine continuousOn_range_of_comp (f ⁻¹ᵁ V.1) _ ?_
  letI := arcTopology (f ⁻¹ᵁ V.1).toScheme
  letI := arcTopology V.1.toScheme
  have he : (fun q : Spec (CommRingCat.of ℂ) ⟶ (f ⁻¹ᵁ V.1).toScheme =>
        (q ≫ (f ⁻¹ᵁ V.1).ι) ≫ f)
      = (fun r : Spec (CommRingCat.of ℂ) ⟶ V.1.toScheme => r ≫ V.1.ι)
        ∘ (fun q => q ≫ (f ∣_ V.1)) := by
    funext q
    show (q ≫ (f ⁻¹ᵁ V.1).ι) ≫ f = (q ≫ (f ∣_ V.1)) ≫ V.1.ι
    rw [Category.assoc, Category.assoc, morphismRestrict_ι]
  rw [he]
  refine (continuous_comp_openImmersion V.1.ι).comp ?_
  rw [← arcTopologyOpen_eq_arcTopology V]
  refine continuous_induced_rng.2 ?_
  letI := arcTopologyAffine Γ(Y, V.1)
  have he2 : (fun q : Spec (CommRingCat.of ℂ) ⟶ (f ⁻¹ᵁ V.1).toScheme =>
        (q ≫ (f ∣_ V.1)) ≫ V.2.isoSpec.hom)
      = fun q => q ≫ ((f ∣_ V.1) ≫ V.2.isoSpec.hom) :=
    funext (fun q => Category.assoc _ _ _)
  show Continuous (fun q => (q ≫ (f ∣_ V.1)) ≫ V.2.isoSpec.hom)
  rw [he2]
  exact continuous_comp_toAffine ((f ∣_ V.1) ≫ V.2.isoSpec.hom)

/-! ## ★出典の紐付け(`.src`) -/

def continuous_comp_scheme.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C——弧空間が関手的であること)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
