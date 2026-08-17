import ABC3.Found.Arakelov.ArcConjInvol
import Mathlib.Analysis.Complex.Basic

/-!
# Arakelov (C1) の第三段 —— **アフィンの `X^arc` の位相**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★これが `ArcSpaceData.topology_affine` である

`Interface/Arakelov/ArcSpace.lean` の C1 は

    topology (Spec A) = induced (evalAffine A) Pi.topologicalSpace

を要求する。★本ファイルはその右辺を**定義として置き**、性質を証明する。

## ★★なぜこの位相なのか

★アフィン `Spec A` の複素点は環準同型 `A → ℂ` そのもの(`ArcEval.lean` の
`evalHom_injective` / `evalHom_Spec_map` で全単射)。
★★その上の自然な位相は**各点収束**であり、それが `ℂ^A` の積位相からの誘導である。

★★★**離散位相ではない。**退化検査で「離散位相なら `conj_continuous` は自明」と
分かったので、C1 はこの条件を課している。

## ★★本ファイルが取るもの

| 定理 | C1 のどの場か |
|---|---|
| `arcTopologyAffine` | `topology (Spec A)` |
| `continuous_conjPoint_affine` | ★★`conj_continuous`(アフィンの場合) |
| `t2Space_arcAffine` | ★分離性(`ArcModel` の議論が要求する) |
| `continuous_evalAffine` | 各切断は連続関数である |
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory ABC3.Found.GenEll

/-! ## ★★アフィンの位相 -/

/-- ★★★**アフィン `Spec A` の複素点の位相** —— 各点収束。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★`ℂ^A` の積位相を `evalAffine` で引き戻したもの。 -/
@[reducible] noncomputable def arcTopologyAffine (A : CommRingCat.{0}) :
    TopologicalSpace (Spec (CommRingCat.of ℂ) ⟶ Spec A) :=
  TopologicalSpace.induced (evalAffine A) Pi.topologicalSpace

/-- ★**定義そのもの**——`ArcSpaceData.topology_affine` が要求する等式。 -/
theorem arcTopologyAffine_eq (A : CommRingCat.{0}) :
    arcTopologyAffine A
      = TopologicalSpace.induced (evalAffine A) Pi.topologicalSpace := rfl

/-! ## ★切断は連続関数である -/

/-- ★**各切断 `a : A` は `X^arc` 上の連続関数を定める**。

★★これが「`X^arc` 上の連続関数」という原文の言い回しの内実である。 -/
theorem continuous_evalAffine (A : CommRingCat.{0}) (a : A) :
    @Continuous _ _ (arcTopologyAffine A) _ (fun p => evalAffine A p a) := by
  letI := arcTopologyAffine A
  exact (continuous_apply a).comp continuous_induced_dom

/-! ## ★★★複素共役は連続 -/

/-- ★★★**[GenEll] Definition 1.1, (i)** —— `ι_X` は連続である(アフィンの場合)。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★`ArcSpaceData.conj_continuous` のアフィンの場合。
★機構は `evalAffine_conjPoint`(`ι_X` は値を複素共役にする)と
`ℂ` 上の複素共役の連続性である。 -/
theorem continuous_conjPoint_affine (A : CommRingCat.{0}) :
    @Continuous _ _ (arcTopologyAffine A) (arcTopologyAffine A)
      (conjPoint (X := Spec A)) := by
  letI := arcTopologyAffine A
  refine continuous_induced_rng.2 (continuous_pi fun a => ?_)
  simp only [Function.comp_apply, evalAffine_conjPoint]
  exact Complex.continuous_conj.comp (continuous_evalAffine A a)

/-! ## ★★分離性 -/

/-- ★★**アフィンの `X^arc` は Hausdorff**。

★機構は 2 つ:評価が点を決めること(`evalHom_injective`)と、`ℂ^A` が Hausdorff。
★★★`Example 1.3, (ii)` の「compact domain」が意味を持つには分離性が要る。 -/
theorem t2Space_arcAffine (A : CommRingCat.{0}) :
    @T2Space _ (arcTopologyAffine A) := by
  letI := arcTopologyAffine A
  have hinj : Function.Injective (evalAffine A) := fun p q h =>
    eq_of_evalAffine_eq A (fun a => congrFun h a)
  exact Topology.IsEmbedding.t2Space ⟨⟨rfl⟩, hinj⟩

/-! ## ★出典の紐付け(`.src`) -/

def arcTopologyAffine.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C——アフィンでの X^arc の位相のみ)",
    sectionId := "genell-def-1-1-i" }

def continuous_conjPoint_affine.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C——ι_X の連続性、アフィンの場合)",
    sectionId := "genell-def-1-1-i" }

def t2Space_arcAffine.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C——アフィンでの X^arc の分離性)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
