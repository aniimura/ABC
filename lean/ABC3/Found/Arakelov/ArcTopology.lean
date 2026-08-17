import ABC3.Found.Arakelov.ArcTopologyAffine

/-!
# Arakelov (C1) の第四段 —— **一般の `X^arc` の位相**(貼り合わせ)(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★構成

アフィンでの位相は `ArcTopologyAffine.lean` で得た。一般の `X` は貼り合わせる:

1. アフィン開 `U ⊆ X` に対し、`U ≅ Spec Γ(X,U)`(`IsAffineOpen.isoSpec`)で
   アフィンの位相を**輸送**する(`arcTopologyOpen`)。
2. 各 `U` の chart `Arc U → Arc X`(`p ↦ p ≫ U.ι`)で **coinduced** を取り、
   すべての `U` について **`⨅`(共通部分)** を取る(`arcTopology`)。

★★★**`⨅` である理由**: mathlib では `t ≤ s ⟺ t の開集合が少ない`(粗い)なので、
`⨅` は**開集合の共通部分**を取る。★「`V` が開 ⟺ どの chart でも `V∩U` が開」
という貼り合わせの条件が、ちょうどこれである。

## ★★本ファイルが取るもの

| 定理 | C1 のどの場か |
|---|---|
| `arcTopology` | `topology X`(一般) |
| `continuous_conjPoint` | ★★★**`conj_continuous`(一般の `X`)** |

★`conjPoint` が一般でも連続なのは、★★**複素共役が `Spec ℂ` の側に作用する**からである
——どの chart も保たれる(`conjPoint_comp`)。
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory ABC3.Found.GenEll TopologicalSpace

/-! ## ★★アフィン開の上の位相(輸送) -/

/-- ★**アフィン開 `U ⊆ X` の複素点の位相** —— `U ≅ Spec Γ(X,U)` で輸送する。 -/
@[reducible] noncomputable def arcTopologyOpen {X : Scheme.{0}} (U : X.affineOpens) :
    TopologicalSpace (Spec (CommRingCat.of ℂ) ⟶ U.1.toScheme) :=
  TopologicalSpace.induced (fun p => p ≫ U.2.isoSpec.hom)
    (arcTopologyAffine (X.presheaf.obj (Opposite.op U.1)))

/-! ## ★★★一般の位相 -/

/-- ★★★**一般の `X` の複素点の位相** —— アフィン開被覆から貼り合わせる。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★`V` が開 ⟺ どのアフィン chart に引き戻しても開、という条件である。 -/
@[reducible] noncomputable def arcTopology (X : Scheme.{0}) :
    TopologicalSpace (Spec (CommRingCat.of ℂ) ⟶ X) :=
  ⨅ U : X.affineOpens,
    TopologicalSpace.coinduced (fun p : Spec (CommRingCat.of ℂ) ⟶ U.1.toScheme => p ≫ U.1.ι)
      (arcTopologyOpen U)

/-- ★各 chart の coinduced は `arcTopology` より細かい(上にある)。 -/
theorem arcTopology_le {X : Scheme.{0}} (U : X.affineOpens) :
    arcTopology X ≤
      TopologicalSpace.coinduced
        (fun p : Spec (CommRingCat.of ℂ) ⟶ U.1.toScheme => p ≫ U.1.ι)
        (arcTopologyOpen U) :=
  iInf_le _ U

/-! ## ★★★複素共役は一般でも連続 -/

/-- ★**アフィン開の上でも `ι_X` は連続**。

★輸送した位相なので、アフィンの場合(`continuous_conjPoint_affine`)から出る。 -/
theorem continuous_conjPoint_open {X : Scheme.{0}} (U : X.affineOpens) :
    @Continuous _ _ (arcTopologyOpen U) (arcTopologyOpen U)
      (conjPoint (X := U.1.toScheme)) := by
  letI := arcTopologyOpen U
  letI := arcTopologyAffine (X.presheaf.obj (Opposite.op U.1))
  refine continuous_induced_rng.2 ?_
  have h : (fun p : Spec (CommRingCat.of ℂ) ⟶ U.1.toScheme =>
      conjPoint p ≫ U.2.isoSpec.hom)
      = fun p => conjPoint (p ≫ U.2.isoSpec.hom) := by
    funext p; rw [conjPoint_comp]
  simp only [Function.comp_def, h]
  exact (continuous_conjPoint_affine _).comp continuous_induced_dom

/-- ★★★**[GenEll] Definition 1.1, (i)** —— `ι_X` は連続である(一般の `X`)。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★`ArcSpaceData.conj_continuous` そのもの。
★機構は「**複素共役は `Spec ℂ` の側に作用する**ので、どのアフィン chart も保たれる」
(`conjPoint_comp`)である。 -/
theorem continuous_conjPoint {X : Scheme.{0}} :
    @Continuous _ _ (arcTopology X) (arcTopology X) (conjPoint (X := X)) := by
  refine continuous_iInf_rng.2 fun U => ?_
  refine continuous_iInf_dom (i := U) ?_
  letI := arcTopologyOpen U
  letI : TopologicalSpace (Spec (CommRingCat.of ℂ) ⟶ X) :=
    TopologicalSpace.coinduced
      (fun p : Spec (CommRingCat.of ℂ) ⟶ U.1.toScheme => p ≫ U.1.ι) (arcTopologyOpen U)
  refine continuous_coinduced_dom.2 ?_
  have h : (conjPoint (X := X)) ∘
        (fun p : Spec (CommRingCat.of ℂ) ⟶ U.1.toScheme => p ≫ U.1.ι)
      = (fun q : Spec (CommRingCat.of ℂ) ⟶ U.1.toScheme => q ≫ U.1.ι) ∘
        (conjPoint (X := U.1.toScheme)) := by
    funext p
    exact conjPoint_comp U.1.ι p
  rw [h]
  exact Continuous.comp continuous_coinduced_rng (continuous_conjPoint_open U)

/-! ## ★出典の紐付け(`.src`) -/

def arcTopology.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C——X^arc の位相の構成のみ)",
    sectionId := "genell-def-1-1-i" }

def continuous_conjPoint.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C——ι_X の連続性)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
