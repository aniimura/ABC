import ABC3.Found.Arakelov.ArcEval
import ABC3.Found.GenEll.ArchConj

/-!
# Arakelov (C1) の第二段 —— **複素共役 `ι_X` は対合であり、評価と両立する**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★`Interface/Arakelov/ArcSpace.lean` の C1 が要求するもの

| 場 | 状態 |
|---|---|
| `Arc X` | ★既存(`complexPoints`) |
| `evalAffine` | ★★済(`ArcEval.lean`) |
| `conj` | ★既存(`ArchConj.lean` の `conjPoint`) |
| **`conj_involutive`** | ★★★**本ファイル** |
| `conj_continuous` | 位相を作ってから |
| `topology_affine` / `topology_openImmersion` | 次段 |

## ★★対合性はどこから来るか

`conjPoint p = conjSpec ≫ p` であり、`conjSpec = Spec.map (starRingEnd ℂ)`。
★`starRingEnd ℂ` は対合(`Complex.conj_conj`)なので、
`Spec` の関手性から `conjSpec ≫ conjSpec = 𝟙` が出る。

★★★**`Spec` は反変**なので合成の順が入れ替わることに注意する
——ここは `Spec.map_comp` の向きを見誤りやすい。

## ★★評価との両立

★★★アフィンでは **`ι_X` は値を複素共役にする**:

    evalAffine A (conjPoint p) a = conj (evalAffine A p a)

★これが「`ι_X` と両立する計量」(C3 の `IsConjCompatible`)の意味を
アフィンで具体化したものである。
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory ABC3.Found.GenEll

/-! ## ★★★`conjSpec` は対合 -/

/-- ★★**`Spec ℂ` の複素共役は対合**。

★`starRingEnd ℂ` が対合であることから、`Spec` の関手性で出る。 -/
@[simp] theorem conjSpec_conjSpec :
    conjSpec ≫ conjSpec = 𝟙 (Spec (CommRingCat.of ℂ)) := by
  rw [conjSpec, ← Spec.map_comp]
  convert Spec.map_id (CommRingCat.of ℂ)
  ext x
  exact Complex.conj_conj x

/-! ## ★★★`ι_X` は対合 -/

/-- ★★★**[GenEll] Definition 1.1, (i)** —— 複素共役 `ι_X` は対合である。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★`Interface/Arakelov/ArcSpace.lean` の `ArcSpaceData.conj_involutive` そのもの。 -/
@[simp] theorem conjPoint_conjPoint {X : Scheme.{0}} (p : complexPoints X) :
    conjPoint (conjPoint p) = p := by
  rw [conjPoint, conjPoint, ← Category.assoc, conjSpec_conjSpec, Category.id_comp]

/-- ★**`ι_X` は全単射**(対合だから)。 -/
theorem conjPoint_involutive {X : Scheme.{0}} :
    Function.Involutive (conjPoint (X := X)) :=
  conjPoint_conjPoint

theorem conjPoint_bijective {X : Scheme.{0}} :
    Function.Bijective (conjPoint (X := X)) :=
  conjPoint_involutive.bijective

/-! ## ★★★評価との両立 -/

/-- ★**`ι_X` は射に沿って自然である**。 -/
theorem conjPoint_comp {X Y : Scheme.{0}} (f : X ⟶ Y) (p : complexPoints X) :
    conjPoint (p ≫ f) = conjPoint p ≫ f := by
  rw [conjPoint, conjPoint, Category.assoc]

/-- ★★★**アフィンでは `ι_X` は値を複素共役にする**。

    `a(ι_X p) = conj (a(p))`

★★これが「計量が `ι_X` と両立する」(C3)の意味をアフィンで具体化したものである。 -/
theorem evalAffine_conjPoint (A : CommRingCat.{0})
    (p : Spec (CommRingCat.of ℂ) ⟶ Spec A) (a : A) :
    evalAffine A (conjPoint p) a = starRingEnd ℂ (evalAffine A p a) := by
  have h : evalHom A (conjPoint p)
      = evalHom A p ≫ CommRingCat.ofHom (starRingEnd ℂ) := by
    refine Spec.map_injective ?_
    rw [Spec_map_evalHom, Spec.map_comp, Spec_map_evalHom]
    rfl
  rw [evalAffine, h]
  rfl

/-! ## ★出典の紐付け(`.src`) -/

def conjPoint_conjPoint.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C——ι_X が対合であること)",
    sectionId := "genell-def-1-1-i" }

def evalAffine_conjPoint.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C——ι_X が評価を複素共役にすること)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
