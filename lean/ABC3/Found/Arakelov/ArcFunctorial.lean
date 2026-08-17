import ABC3.Found.Arakelov.ArcTopologyAffine

/-!
# Arakelov (C1) の第五段 —— **評価の関手性と chart の連続性**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★これが `topology_affine` の核である

`Interface/Arakelov/ArcSpace.lean` の C1 は

    topology (Spec A) = induced (evalAffine A) Pi.topologicalSpace

を要求する。★`ArcTopology.lean` の `arcTopology` は `⨅`(アフィン開被覆)で
定義したので、この等式には**両向き**が要る:

| 向き | 意味 | 必要なもの |
|---|---|---|
| `⨅ ≤ induced` | `⊤` も chart なので自明 | ★`⊤` の chart が同相 |
| `induced ≤ ⨅` | **どの chart も連続** | ★★★**本ファイル** |

## ★★★機構: `Spec.preimage` の関手性

アフィン射 `f : Spec B ⟶ Spec A` に沿って

    evalHom A (p ≫ f) = Spec.preimage f ≫ evalHom B p

が成り立つ(`Spec` が反変な充満忠実関手だから)。★したがって

    a(p ≫ f) = (f^♯ a)(p)

——**切断を引き戻してから評価しても同じ**。★★これで chart の連続性が出る。
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory

/-! ## ★★★評価の関手性 -/

/-- ★★★**評価はアフィン射に沿って関手的である**。

    `evalHom A (p ≫ f) = Spec.preimage f ≫ evalHom B p`

★`Spec` は反変なので合成の順が入れ替わる。 -/
theorem evalHom_comp {A B : CommRingCat.{0}} (f : Spec B ⟶ Spec A)
    (p : Spec (CommRingCat.of ℂ) ⟶ Spec B) :
    evalHom A (p ≫ f) = Spec.preimage f ≫ evalHom B p := by
  refine Spec.map_injective ?_
  rw [Spec_map_evalHom, Spec.map_comp, Spec_map_evalHom, Spec.map_preimage]

/-- ★★**切断を引き戻してから評価しても同じ**。

    `a(p ≫ f) = (f^♯ a)(p)` -/
theorem evalAffine_comp {A B : CommRingCat.{0}} (f : Spec B ⟶ Spec A)
    (p : Spec (CommRingCat.of ℂ) ⟶ Spec B) (a : A) :
    evalAffine A (p ≫ f) a = evalAffine B p ((Spec.preimage f).hom a) := by
  rw [evalAffine, evalHom_comp]
  rfl

/-! ## ★★★合成は連続 -/

/-- ★★★**アフィン射との合成は連続写像を定める**。

★★これが chart の連続性であり、`topology_affine` の `induced ≤ ⨅` の向きを出す。 -/
theorem continuous_comp_affine {A B : CommRingCat.{0}} (f : Spec B ⟶ Spec A) :
    @Continuous _ _ (arcTopologyAffine B) (arcTopologyAffine A)
      (fun p : Spec (CommRingCat.of ℂ) ⟶ Spec B => p ≫ f) := by
  letI := arcTopologyAffine B
  letI := arcTopologyAffine A
  refine continuous_induced_rng.2 (continuous_pi fun a => ?_)
  have h : (fun p : Spec (CommRingCat.of ℂ) ⟶ Spec B => evalAffine A (p ≫ f) a)
      = fun p => evalAffine B p ((Spec.preimage f).hom a) := by
    funext p; exact evalAffine_comp f p a
  simp only [Function.comp_apply, h]
  exact continuous_evalAffine B _

/-! ## ★出典の紐付け(`.src`) -/

def evalHom_comp.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C——評価がアフィン射に沿って関手的であること)",
    sectionId := "genell-def-1-1-i" }

def continuous_comp_affine.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C——アフィン射との合成の連続性)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
