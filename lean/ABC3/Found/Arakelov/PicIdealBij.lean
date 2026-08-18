import ABC3.Found.Arakelov.PicIdealScal

/-!
# Arakelov (B2) 第 158 ブロック —— **基本開集合での全単射性**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★第 131 の `idealSections` 版——**第 130 をそのまま再利用**

第 130(`bijective_smul_liftGen`)は `R`・`M` について一般に述べてあるので、
`R := Γ(X, A)`、`M := D.ideal A` を代入するだけでよい。

    Γ(X, D(t·g)) --awayRingEquivX--> R_{t·g} --(· • y)--> I_{t·g}
      --idealAwayEquivScalar--> idealSections D (D(t·g))

★★3 段の合成で、両端は全単射(環同型・線型同型)。**一発で通った。**
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite TopologicalSpace

variable {X : Scheme.{u}} (D : X.IdealSheafData) (A : X.affineOpens)

theorem bijective_smul_idealGen (g t : (Γ(X, A.1) : Type u))
    (e : LocalizedModule (Submonoid.powers g) (D.ideal A)
      ≃ₗ[Localization (Submonoid.powers g)] Localization (Submonoid.powers g)) :
    Function.Bijective (fun c : (Γ(X, X.basicOpen (t * g)) : Type u) =>
      c • (idealAwayEquivScalar D A (t * g)
        (liftAwayMapA (Γ(X, A.1) : Type u) g t (D.ideal A) (e.symm 1)))) := by
  have h1 := bijective_smul_liftGen (Γ(X, A.1) : Type u) g t (D.ideal A) e
  have hcomp : (fun c : (Γ(X, X.basicOpen (t * g)) : Type u) =>
      c • (idealAwayEquivScalar D A (t * g)
        (liftAwayMapA (Γ(X, A.1) : Type u) g t (D.ideal A) (e.symm 1))))
      = (idealAwayEquivScalar D A (t * g))
        ∘ (fun d : Localization (Submonoid.powers (t * g)) =>
            d • liftAwayMapA (Γ(X, A.1) : Type u) g t (D.ideal A) (e.symm 1))
        ∘ (awayRingEquivX A (t * g)) := by
    funext c
    show c • (idealAwayEquivScalar D A (t * g)) _ = _
    rw [← map_smul]
    rfl
  rw [hcomp]
  exact ((idealAwayEquivScalar D A (t * g)).bijective.comp h1).comp
    (awayRingEquivX A (t * g)).bijective


/-! ## ★出典の紐付け(`.src`) -/

def bijective_smul_idealGen.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——基本開集合での全単射性)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
