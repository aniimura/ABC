import ABC3.Found.Arakelov.ArcTopology
import ABC3.Found.Arakelov.ArcFunctorial

/-!
# Arakelov (C1) の第六段 —— **`topology_affine`**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★両向きが「chart の連続性」から出る

`Interface/Arakelov/ArcSpace.lean` の C1 は

    topology (Spec A) = induced (evalAffine A) Pi.topologicalSpace

を要求する。★我々の `arcTopology` は `⨆`(アフィン開被覆の coinduced)なので、
証明すべきは

    ⨆ U, coinduced (chart U) (arcTopologyOpen U) = arcTopologyAffine A

★★★**両向きとも `continuous_comp_affine` に帰着する**:

| 向き | 使うもの |
|---|---|
| `⨆ ≤ affine` | ★各 chart の連続性(`continuous_iff_coinduced_le`) |
| `affine ≤ ⨆` | ★★`U = ⊤` の chart の**逆写像**の連続性 |

★★`a ≤ b ⟺ a が細かい`(`⊥` が離散であることで確定。2026-08-17 実測)。

## ★★chart の連続性はどこから来るか

`chart U p = p ≫ U.ι` を

    p ≫ U.ι = (p ≫ isoSpec.hom) ≫ (isoSpec.inv ≫ U.ι)

と分解する。★`isoSpec.inv ≫ U.ι : Spec Γ(Spec A, U) ⟶ Spec A` は**アフィン射**なので、
`ArcFunctorial.lean` の `continuous_comp_affine` が使える。
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory TopologicalSpace

/-! ## ★★★chart の連続性(アフィンの場合) -/

/-- ★★★**アフィン `Spec A` のどの chart も連続である**。

★分解 `p ≫ U.ι = (p ≫ isoSpec.hom) ≫ (isoSpec.inv ≫ U.ι)` を使う。
★★`isoSpec.inv ≫ U.ι` はアフィン射なので `continuous_comp_affine` が効く。 -/
theorem continuous_chart_affine (A : CommRingCat.{0}) (U : (Spec A).affineOpens) :
    @Continuous _ _ (arcTopologyOpen U) (arcTopologyAffine A)
      (fun p : Spec (CommRingCat.of ℂ) ⟶ U.1.toScheme => p ≫ U.1.ι) := by
  letI := arcTopologyOpen U
  letI := arcTopologyAffine ((Spec A).presheaf.obj (Opposite.op U.1))
  letI := arcTopologyAffine A
  have hsplit : (fun p : Spec (CommRingCat.of ℂ) ⟶ U.1.toScheme => p ≫ U.1.ι)
      = (fun q : Spec (CommRingCat.of ℂ) ⟶
            Spec ((Spec A).presheaf.obj (Opposite.op U.1)) =>
          q ≫ (U.2.isoSpec.inv ≫ U.1.ι))
        ∘ (fun p : Spec (CommRingCat.of ℂ) ⟶ U.1.toScheme => p ≫ U.2.isoSpec.hom) := by
    funext p
    simp only [Function.comp_apply, Category.assoc, Iso.hom_inv_id_assoc]
  rw [hsplit]
  exact Continuous.comp (continuous_comp_affine _) continuous_induced_dom

/-! ## ★★逆向きの連続性(`⊤` の chart の逆写像) -/

/-- ★**`⊤` の chart の逆写像**——`q ↦ q ≫ topIso.inv`。 -/
noncomputable def chartTopInv (A : CommRingCat.{0})
    (q : Spec (CommRingCat.of ℂ) ⟶ Spec A) :
    Spec (CommRingCat.of ℂ) ⟶ (⊤ : (Spec A).Opens).toScheme :=
  q ≫ (Spec A).topIso.inv

/-- ★★**逆写像は連続**。

★`topIso.inv ≫ isoSpec.hom : Spec A ⟶ Spec Γ(Spec A, ⊤)` は**アフィン射**なので
`continuous_comp_affine` が効く。 -/
theorem continuous_chartTopInv (A : CommRingCat.{0}) :
    @Continuous _ _ (arcTopologyAffine A)
      (arcTopologyOpen ⟨⊤, isAffineOpen_top (Spec A)⟩) (chartTopInv A) := by
  letI := arcTopologyAffine A
  letI := arcTopologyOpen (⟨⊤, isAffineOpen_top (Spec A)⟩ : (Spec A).affineOpens)
  letI := arcTopologyAffine ((Spec A).presheaf.obj (Opposite.op (⊤ : (Spec A).Opens)))
  refine continuous_induced_rng.2 ?_
  have h : (fun q : Spec (CommRingCat.of ℂ) ⟶ Spec A =>
        chartTopInv A q ≫ (isAffineOpen_top (Spec A)).isoSpec.hom)
      = fun q => q ≫ ((Spec A).topIso.inv ≫ (isAffineOpen_top (Spec A)).isoSpec.hom) := by
    funext q; rw [chartTopInv, Category.assoc]
  simp only [Function.comp_def, h]
  exact continuous_comp_affine _

/-- ★**`chart ⊤` と `chartTopInv` は互いに逆**。 -/
theorem chart_chartTopInv (A : CommRingCat.{0})
    (q : Spec (CommRingCat.of ℂ) ⟶ Spec A) :
    chartTopInv A q ≫ (⊤ : (Spec A).Opens).ι = q := by
  rw [chartTopInv, Category.assoc, Scheme.toIso_inv_ι, Category.comp_id]

/-! ## ★★★`topology_affine` -/

/-- ★★★**[GenEll] Definition 1.1, (i)** —— アフィンでの `X^arc` の位相は各点収束である。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★`ArcSpaceData.topology_affine` そのもの:

    arcTopology (Spec A) = induced (evalAffine A) Pi.topologicalSpace

★★両向きとも `continuous_comp_affine`(アフィン射との合成の連続性)に帰着する。 -/
theorem arcTopology_spec (A : CommRingCat.{0}) :
    arcTopology (Spec A) = arcTopologyAffine A := by
  refine le_antisymm ?_ ?_
  · -- `⨆ ≤ affine`: 各 chart の連続性(`continuous_iff_coinduced_le`)
    refine iSup_le fun U => ?_
    exact continuous_iff_coinduced_le.1 (continuous_chart_affine A U)
  · -- `affine ≤ ⨆`: `U = ⊤` を経由し、逆写像の連続性を使う
    letI := arcTopologyAffine A
    letI := arcTopologyOpen (⟨⊤, isAffineOpen_top (Spec A)⟩ : (Spec A).affineOpens)
    refine le_trans ?_ (le_arcTopology ⟨⊤, isAffineOpen_top (Spec A)⟩)
    refine fun V hV => ?_
    have h2 : V = chartTopInv A ⁻¹'
        ((fun p : Spec (CommRingCat.of ℂ) ⟶ (⊤ : (Spec A).Opens).toScheme =>
          p ≫ (⊤ : (Spec A).Opens).ι) ⁻¹' V) := by
      ext q
      simp only [Set.mem_preimage, chart_chartTopInv]
    rw [h2]
    exact (continuous_chartTopInv A).isOpen_preimage _ hV

/-! ## ★出典の紐付け(`.src`) -/

def continuous_chart_affine.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C——アフィンでの chart の連続性)",
    sectionId := "genell-def-1-1-i" }

def arcTopology_spec.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C——アフィンでの X^arc の位相が各点収束であること)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
