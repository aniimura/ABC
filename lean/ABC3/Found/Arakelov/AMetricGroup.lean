/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Arakelov.AMetricIso
import ABC3.Meta.Claim

/-!
# 算術直線束の**単位律**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★★★`L̄ ⊗ Ō_X ≅ L̄`（等長）

原文の `APic(X)` が群であることの欄のうち、**単位律**である。

    `IsIsometry (L * 1) L (ρ_ L.sheaf)`   （`isIsometry_mul_one`）

★機構は 2 つ:

1. ★`pullTriv (ρ_ A) V eA = tensorTriv eA (baseTriv X V)`
   ——右単位子は「基準の自明化とのテンソル」である。
   その中身は mathlib の `Functor.OplaxMonoidal.right_unitality_hom` と
   `MonoidalCategory.unitors_equal`・`rightUnitor_naturality`。
2. ★★`structLocalMetric` は基準の自明化で `h = 1`（`structLocalMetric_h_baseTriv`）。

## ★残っている段（明示）

★★**結合律**——`pullTriv (α_ A B C) V (tensorTriv eA (tensorTriv eB eC))`
`= tensorTriv (tensorTriv eA eB) eC` は **`rfl` ではない**（2026-08-28 実測）。
★機構は分かっている（`Functor.OplaxMonoidal.oplax_associativity` ＋
`associator_naturality` ＋ `triangle` ＋ `unitors_equal`）が、
`simp` は `λ_` を関手の lax 構造へ展開してしまい**逆方向へ進む**ので、
手で並べる必要がある。★★これは 1 ブロックぶんの仕事として残す。

★★★**逆元（双対計量）** と **商そのもの**もまだである。
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory MonoidalCategory Opposite
open ABC3.Found.GenEll

variable {X : Scheme.{0}}

/-- ★★★★**右単位子は「基準の自明化とのテンソル」である**。

★機構は mathlib の `right_unitality_hom` ＋ `unitors_equal` ＋ `rightUnitor_naturality`。 -/
theorem pullTriv_rightUnitor {A : X.PresheafOfModules} (V : X.Opens)
    (eA : (restrictPresheafFunctor X V).obj A ≅ 𝟙_ (PresheafModulesOn X V)) :
    pullTriv (ρ_ A) V eA = tensorTriv eA (baseTriv X V) := by
  ext1
  show ((restrictPresheafFunctor X V).mapIso (ρ_ A)).hom ≫ eA.hom
    = (tensorTriv eA (baseTriv X V)).hom
  simp only [tensorTriv, baseTriv, restrictPresheafTensor, restrictPresheafUnit,
    Iso.trans_hom, Iso.symm_hom, tensorIso_hom, Functor.Monoidal.μIso_inv,
    Category.assoc, MonoidalCategory.tensorHom_def']
  rw [MonoidalCategory.unitors_equal, MonoidalCategory.rightUnitor_naturality]
  have hru : Functor.OplaxMonoidal.δ (restrictPresheafFunctor X V) A (𝟙_ X.PresheafOfModules)
        ≫ ((restrictPresheafFunctor X V).obj A
            ◁ (Functor.Monoidal.εIso (restrictPresheafFunctor X V)).inv)
        ≫ (ρ_ ((restrictPresheafFunctor X V).obj A)).hom
        ≫ eA.hom
      = (restrictPresheafFunctor X V).map (ρ_ A).hom ≫ eA.hom :=
    Functor.OplaxMonoidal.right_unitality_hom_assoc
      (F := restrictPresheafFunctor X V) A eA.hom
  rw [hru]
  rfl

/-- ★★★★★★★★**単位律** `L̄ ⊗ Ō_X ≅ L̄`（等長）。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★`Ō_X` の基準ノルムは基準の自明化で `1` なので、テンソル積の `h` は変わらない。 -/
theorem isIsometry_mul_one (L : AMetric X) :
    IsIsometry (L * 1) L (ρ_ L.sheaf) := by
  intro V e p hp
  have hten : (L * 1).metric.h V (tensorTriv e (baseTriv X V)) p
      = L.metric.h V e p * (1 : AMetric X).metric.h V (baseTriv X V) p :=
    isTensorOf_tensor L.triv (1 : AMetric X).triv L.metric (1 : AMetric X).metric V e
      (baseTriv X V) p hp
  have hone : (1 : AMetric X).metric.h V (baseTriv X V) p = 1 :=
    structLocalMetric_h_baseTriv X V p hp
  have hpull : pullTriv (ρ_ L.sheaf) V e = tensorTriv e (baseTriv X V) :=
    pullTriv_rightUnitor V e
  have hstep : (L * 1).metric.h V (pullTriv (ρ_ L.sheaf) V e) p
      = (L * 1).metric.h V (tensorTriv e (baseTriv X V)) p :=
    congrArg (fun t => (L * 1).metric.h V t p) hpull
  exact hstep.trans (hten.trans (by rw [hone, mul_one]))

/-- ★★**したがって同値類の上で `[L̄] · 1 = [L̄]`**。 -/
theorem isometric_mul_one (L : AMetric X) : Isometric (L * 1) L :=
  ⟨ρ_ L.sheaf, isIsometry_mul_one L⟩

/-! ### ★出典の紐付け(`.src`) -/

def pullTriv_rightUnitor.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(右単位子は基準の自明化とのテンソルであること)",
    sectionId := "genell-def-1-1-i" }

def isIsometry_mul_one.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(単位律 L̄ ⊗ Ō_X ≅ L̄——等長)",
    sectionId := "genell-def-1-1-i" }

def isIsometry_mul_one.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "isTensorOf_tensor(構成した計量がテンソル積の計量であること)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.isTensorOf_tensor") 3,
    .citation "[ABC3]" "structLocalMetric_h_baseTriv(Ō_X は基準の自明化で h = 1)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.structLocalMetric_h_baseTriv") 3,
    .citation "[mathlib]" "Functor.OplaxMonoidal.right_unitality_hom"
      (.inMathlib "CategoryTheory.Functor.OplaxMonoidal.right_unitality_hom") 3,
    .implicitStep
      ("★★★★**結合律は残っている**——" ++
       "pullTriv (α_ A B C) V (tensorTriv eA (tensorTriv eB eC)) " ++
       "= tensorTriv (tensorTriv eA eB) eC は rfl ではない(2026-08-28 実測)。" ++
       "★機構は分かっている(oplax_associativity ＋ associator_naturality ＋ " ++
       "triangle ＋ unitors_equal)が、simp は λ_ を関手の lax 構造へ展開して" ++
       "逆方向へ進むので手で並べる必要がある") 3,
    .implicitStep
      ("★★★逆元(双対計量)と商そのものもまだである。" ++
       "★Definition 1.1 の項目全体には (ii) の deg_F も要る") 3 ]

end ABC3.Found.Arakelov
