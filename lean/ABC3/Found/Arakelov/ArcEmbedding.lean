import ABC3.Found.Arakelov.ArcOpenSetOpen
import ABC3.Found.Arakelov.ArcOpenImmersion
import Mathlib.Topology.Maps.Basic
import Mathlib.Topology.Defs.Induced

/-!
# Arakelov (C3) 第 275 ブロック —— ★★★★**開埋め込みに沿った連続性の移送**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★第 273 の仮定 `hg` を落とすための器具

`(· ≫ V.ι) : V^arc → X^arc` は**開埋め込み**である:

| 条件 | 出どころ |
|---|---|
| inducing | ★C1 `arcTopology_openImmersion` |
| injective | ★C1 `comp_openImmersion_injective` |
| 像が開 | ★第 274 `isOpen_range_comp_ι` |

★★したがって `ContinuousOn g (range (· ≫ V.ι))` は
**`Continuous (g ∘ (· ≫ V.ι))` から出る**(mathlib `IsOpenEmbedding.continuousAt_iff`)。

★★★これで第 270(`V.toScheme` 上の連続性)を `X^arc` の開部分集合上へ移せる。

## ★摩擦 —— `IsOpenEmbedding` は `Topology` 名前空間にある

`import` を 2 つ足しても `unknownIdentifier` が消えず、
★原因は **`namespace Topology` の中**だったこと。`open Topology` で通る。
★★**「import したのに見えない」ときは名前空間も疑う**——
本セッションの import 非推移(3 例)とは別の原因である。

| 定理 | 内容 |
|---|---|
| `isOpenEmbedding_comp_ι` | ★★開埋め込み |
| `continuousOn_range_of_comp` | ★★★★連続性の移送 |
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite TopologicalSpace Topology

variable {X : Scheme.{0}} (V : X.Opens)

/-- ★★`(· ≫ V.ι)` は開埋め込みである。 -/
theorem isOpenEmbedding_comp_ι :
    @IsOpenEmbedding _ _ (arcTopology V.toScheme) (arcTopology X)
      (fun q : Spec (CommRingCat.of ℂ) ⟶ V.toScheme => q ≫ V.ι) := by
  letI := arcTopology V.toScheme
  letI := arcTopology X
  refine ⟨⟨⟨?_⟩, comp_openImmersion_injective V.ι⟩, ?_⟩
  · exact arcTopology_openImmersion V.ι
  · exact isOpen_range_comp_ι V

/-- ★★★★像の上での連続性は合成の連続性から出る。 -/
theorem continuousOn_range_of_comp {C : Type} [TopologicalSpace C]
    (g : (Spec (CommRingCat.of ℂ) ⟶ X) → C)
    (h : @Continuous _ C (arcTopology V.toScheme) _
      (fun q => g (q ≫ V.ι))) :
    @ContinuousOn _ C (arcTopology X) _ g
      (Set.range (fun q : Spec (CommRingCat.of ℂ) ⟶ V.toScheme => q ≫ V.ι)) := by
  letI := arcTopology V.toScheme
  letI := arcTopology X
  rintro y ⟨x, rfl⟩
  exact (((isOpenEmbedding_comp_ι V).continuousAt_iff (g := g)).1
    h.continuousAt).continuousWithinAt


/-! ## ★出典の紐付け(`.src`) -/

def continuousOn_range_of_comp.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C——開埋め込みに沿った連続性の移送)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
