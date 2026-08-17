import ABC3.Found.Arakelov.ArcTopologyAffine

/-!
# Arakelov (C1) の部品 —— **基本開集合の複素点は開集合**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★これが `topology_openImmersion` の残り 1 向きの核である

残っているのは「**`(· ≫ f)` が開写像**」という主張である。★その証明では
`Arc W` の中で「点が `W ∩ f(X)` に落ちる」という部分集合が開であることを使う。

★★アフィン `W = Spec A` では `W ∩ f(X)` は開なので基本開集合 `D(g)` の合併であり、

    {p : Arc (Spec A) | g(p) ≠ 0}

が開であることに帰着する。★★★**それが本ファイルである。**

## ★★機構

`g(p) = evalAffine A p g` は **`p` について連続**(`continuous_evalAffine`)。
★したがって `{p | g(p) ≠ 0}` は `ℂ \ {0}`(開)の逆像である。
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory

/-! ## ★★★基本開集合 -/

/-- ★**切断 `g` が消えない複素点の集合**(基本開集合 `D(g)` の複素点)。 -/
def arcBasicOpen (A : CommRingCat.{0}) (g : A) :
    Set (Spec (CommRingCat.of ℂ) ⟶ Spec A) :=
  {p | evalAffine A p g ≠ 0}

/-- ★★★**基本開集合の複素点は開集合である**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★機構は `continuous_evalAffine`(切断は連続関数)と `ℂ \ {0}` が開であること。
★★★これが「`X^arc` の位相が正しい」ことの具体的な現れである
——離散位相ならこの主張は無内容になる。 -/
theorem isOpen_arcBasicOpen (A : CommRingCat.{0}) (g : A) :
    @IsOpen _ (arcTopologyAffine A) (arcBasicOpen A g) := by
  letI := arcTopologyAffine A
  have h : arcBasicOpen A g = (fun p => evalAffine A p g) ⁻¹' {z : ℂ | z ≠ 0} := rfl
  rw [h]
  exact (continuous_evalAffine A g).isOpen_preimage _ isOpen_ne

/-- ★**逆に、`g` が消える点の集合は閉集合**。 -/
theorem isClosed_arcZeroLocus (A : CommRingCat.{0}) (g : A) :
    @IsClosed _ (arcTopologyAffine A) {p | evalAffine A p g = 0} := by
  letI := arcTopologyAffine A
  have h : {p | evalAffine A p g = 0}
      = (fun p => evalAffine A p g) ⁻¹' {z : ℂ | z = 0} := rfl
  rw [h]
  exact IsClosed.preimage (continuous_evalAffine A g) (isClosed_singleton (x := (0:ℂ)))

/-- ★★**基本開集合の合併も開**(開被覆を扱うために要る)。 -/
theorem isOpen_iUnion_arcBasicOpen (A : CommRingCat.{0}) {ι : Type} (gs : ι → A) :
    @IsOpen _ (arcTopologyAffine A) (⋃ i, arcBasicOpen A (gs i)) := by
  letI := arcTopologyAffine A
  exact isOpen_iUnion fun i => isOpen_arcBasicOpen A (gs i)

/-! ## ★★★局所化の延長の解析的核心 -/

/-- ★★★**`φ ↦ φ(b)/φ(s)` は `{φ(s) ≠ 0}` 上で連続である**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★**これが C1 の残る 1 点(局所化の延長の連続性)の解析的核心である。**

`Arc (Spec A) = Hom_Ring(A, ℂ)` の上で、基本開集合 `D(s)` への持ち上げは
`b/sⁿ ↦ φ(b)/φ(s)ⁿ` で与えられる。★その各成分の連続性が本定理である。

★機構は `continuous_evalAffine`(切断は連続)と
**分母が消えない所での除算の連続性**(`ContinuousOn.div`)。 -/
theorem continuousOn_div_evalAffine (A : CommRingCat.{0}) (b s : A) :
    @ContinuousOn _ _ (arcTopologyAffine A) _
      (fun p => evalAffine A p b / evalAffine A p s) (arcBasicOpen A s) := by
  letI := arcTopologyAffine A
  exact ContinuousOn.div (continuous_evalAffine A b).continuousOn
    (continuous_evalAffine A s).continuousOn (fun p hp => hp)

/-- ★★**`φ(b)/φ(s)ⁿ` の形でも連続**(局所化の一般の元に対応)。 -/
theorem continuousOn_div_pow_evalAffine (A : CommRingCat.{0}) (b s : A) (n : ℕ) :
    @ContinuousOn _ _ (arcTopologyAffine A) _
      (fun p => evalAffine A p b / (evalAffine A p s) ^ n) (arcBasicOpen A s) := by
  letI := arcTopologyAffine A
  refine ContinuousOn.div (continuous_evalAffine A b).continuousOn
    ((continuous_evalAffine A s).continuousOn.pow n) (fun p hp => ?_)
  exact pow_ne_zero n hp

/-! ## ★出典の紐付け(`.src`) -/

def isOpen_arcBasicOpen.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C——基本開集合の複素点が開であること)",
    sectionId := "genell-def-1-1-i" }

def isClosed_arcZeroLocus.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C——切断の零点集合が閉であること)",
    sectionId := "genell-def-1-1-i" }

def continuousOn_div_evalAffine.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C——局所化の延長の連続性の核心)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
