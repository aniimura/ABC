import ABC3.Found.Arakelov.ArcTrivNorm
import ABC3.Found.Arakelov.ArcTopology

/-!
# Arakelov (C3) 第 249 ブロック —— **`X^arc` 上の連続性の判定**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★連続性は chart ごとに落ちる

`arcTopology X` は **`⨆`(coinduced)**で定義してある。
★`⨆` から**出る**写像の連続性は、各成分からの連続性に落ちる:

    Continuous[⨆ tᵢ] g  ⟺  ∀ i, Continuous[tᵢ] g          (`continuous_iSup_dom`)
    Continuous[coinduced f t] g  ⟺  Continuous[t] (g ∘ f)  (`continuous_coinduced_dom`)

★★したがって「各アフィン chart で連続」を示せばよい。

★★★さらに `arcTopologyOpen U` は `induced (· ≫ isoSpec.hom) (arcTopologyAffine …)` なので、
**アフィンの上の連続関数への分解**を与えれば済む——これが `continuous_of_charts_factor` である。

★★★★これが `Interface` の `normSection_continuous`(§9-286 系で追加)を
満たすための唯一の入口になる。

| 定理 | 内容 |
|---|---|
| `continuous_of_charts` | ★★★各 chart で連続なら連続 |
| `continuous_of_charts_factor` | ★★★★★**アフィンの連続関数に分解すればよい** |
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite TopologicalSpace

variable {X : Scheme.{0}}

/-- ★★★**`X^arc` からの写像は、各アフィン chart で連続なら連続である**。 -/
theorem continuous_of_charts {Y : Type} [TopologicalSpace Y]
    (g : (Spec (CommRingCat.of ℂ) ⟶ X) → Y)
    (h : ∀ U : X.affineOpens,
      @Continuous _ _ (arcTopologyOpen U) _ (fun p => g (p ≫ U.1.ι))) :
    @Continuous _ _ (arcTopology X) _ g := by
  refine continuous_iSup_dom.2 fun U => ?_
  letI := arcTopologyOpen U
  refine continuous_coinduced_dom.2 ?_
  exact h U

/-- ★★★★★**使いやすい形**——各 chart で「アフィンの上の連続関数に分解する」ことを示せばよい。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★`arcTopologyOpen U = induced (· ≫ isoSpec.hom) (arcTopologyAffine …)` を使う。 -/
theorem continuous_of_charts_factor {Y : Type} [TopologicalSpace Y]
    (g : (Spec (CommRingCat.of ℂ) ⟶ X) → Y)
    (G : ∀ U : X.affineOpens,
      (Spec (CommRingCat.of ℂ) ⟶ Spec (X.presheaf.obj (op U.1))) → Y)
    (hG : ∀ U : X.affineOpens,
      @Continuous _ _ (arcTopologyAffine (X.presheaf.obj (op U.1))) _ (G U))
    (hfac : ∀ (U : X.affineOpens) (p : Spec (CommRingCat.of ℂ) ⟶ U.1.toScheme),
      g (p ≫ U.1.ι) = G U (p ≫ U.2.isoSpec.hom)) :
    @Continuous _ _ (arcTopology X) _ g := by
  refine continuous_of_charts g fun U => ?_
  letI := arcTopologyOpen U
  letI := arcTopologyAffine (X.presheaf.obj (op U.1))
  have he : (fun p : Spec (CommRingCat.of ℂ) ⟶ U.1.toScheme => g (p ≫ U.1.ι))
      = (G U) ∘ (fun p => p ≫ U.2.isoSpec.hom) := funext (hfac U)
  rw [he]
  exact (hG U).comp continuous_induced_dom

/-! ## ★出典の紐付け(`.src`) -/

def continuous_of_charts_factor.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C——X^arc 上の連続性の判定)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
