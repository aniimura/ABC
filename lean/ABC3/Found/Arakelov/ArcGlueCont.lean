import ABC3.Found.Arakelov.ArcGlueX
import Mathlib.Topology.PartitionOfUnity

/-!
# Arakelov (C3) 第 273 ブロック —— ★★★★★**貼り合わせたノルムは連続**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★mathlib の 1 の分割がそのまま効いた

    PartitionOfUnity.IsSubordinate.continuous_finsum_smul
      (ho : ∀ i, IsOpen (U i)) (hf : f.IsSubordinate U) (hg : ∀ i, ContinuousOn (g i) (U i))
      : Continuous fun x => ∑ᶠ i, f i x • g i x

★`E = ℝ` では `•` は `*` なので、`gluedNormX` の形にそのまま合う。

★★§9-285 で「コンパクト Hausdorff ⟹ `NormalSpace` ∧ `ParacompactSpace` は無料」と
実測してあるので、1 の分割の**存在**も在庫で出る。

## ★★残る 2 つの仮定

| 仮定 | 内容 |
|---|---|
| `ho` | ★`arcOpenSet (U i)`(`U i` を通る点の集合)が `X^arc` で開 |
| `hg` | ★★局所ノルムが `U i` 上で連続(第 270 を `X` の点へ移す) |

★★★どちらも C1 の `topology_openImmersion` と第 271 の同型から出る見込みである。

| 定義・定理 | 内容 |
|---|---|
| `arcOpenSet` | ★`V` を通る複素点の集合 |
| `continuous_gluedNormX` | ★★★★★**貼り合わせたノルムの連続性** |
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite TopologicalSpace
open scoped Classical

variable {X : Scheme.{0}} {ι : Type}
/-- ★`V` を通る複素点の集合。 -/
def arcOpenSet (V : X.Opens) : Set (Spec (CommRingCat.of ℂ) ⟶ X) :=
  {p | p ⁻¹ᵁ V = ⊤}

/-- ★★★★★貼り合わせたノルムは連続である。 -/
theorem continuous_gluedNormX (F : X.Modules) (U : ι → X.Opens)
    (e : ∀ i, (restrictPresheafFunctor X (U i)).obj F.val
      ≅ 𝟙_ (PresheafModulesOn X (U i)))
    (f : @PartitionOfUnity ι (Spec (CommRingCat.of ℂ) ⟶ X) (arcTopology X) Set.univ)
    (ho : ∀ i, @IsOpen _ (arcTopology X) (arcOpenSet (U i)))
    (hsub : @PartitionOfUnity.IsSubordinate ι (Spec (CommRingCat.of ℂ) ⟶ X) (arcTopology X)
      Set.univ f (fun i => arcOpenSet (U i)))
    (φ : ∀ p : Spec (CommRingCat.of ℂ) ⟶ X, ↥(arcFiber p F))
    (hg : ∀ i, @ContinuousOn _ ℝ (arcTopology X) _
      (fun p => extNormX F U e i p (φ p)) (arcOpenSet (U i))) :
    @Continuous _ ℝ (arcTopology X) _
      (fun p => gluedNormX F U e (fun i p => f i p) p (φ p)) := by
  letI := arcTopology X
  have h := PartitionOfUnity.IsSubordinate.continuous_finsum_smul (f := f) ho hsub hg
  exact h


/-! ## ★出典の紐付け(`.src`) -/

def continuous_gluedNormX.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C——貼り合わせたノルムは連続)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
