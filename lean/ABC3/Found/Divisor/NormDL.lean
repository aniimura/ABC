/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Divisor.NormFunctor

/-!
# `D_L` —— `D_K` の素因子に写る `V[L]` の素因子(鎖 `normalize` の `DL-set`)

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.109。

原文 (FrdI p.109):
> a proper normal variety], then let us write DL for the set of prime divisors of V [L]

原文 (FrdI p.109):
> that map into [possibly subvarieties of codimension ≥1 of] prime divisors of DK.

## ★★原文の括弧書きが効く所

原文は `D_L` を「`D_K` の素因子**に写る**(余次元が上がってもよい)`V[L]` の素因子」と定める。
★★**「写る」は像が素因子の閉包に入ること**であり、**余次元は保たれなくてよい**
——原文が角括弧で明記しているとおりである。

★スキームの言葉では「`x` の像が `y` の特殊化である」、すなわち
`normDown x ⤳ y` が成り立たない…ではなく、**`normDown x` が `y` の閉包に入る**、
つまり `y ⤳ normDown x`(`y` が `normDown x` へ特殊化する側)である。

## ★本ファイルで閉じること

| 定義 | 中身 |
|---|---|
| `MapsIntoDivisor` | 素因子 `x` が素因子 `y` に写ること |
| `DLSet` | ★★**`D_L`** |
| `DLSet_mono` | `D_K` を大きくすれば `D_L` も大きくなる |
-/

namespace ABC3.Found.Divisor

open AlgebraicGeometry CategoryTheory ABC3.Found.FrdI

universe u

variable (V : Scheme.{u}) [IsIntegral V] {Kbar : Type u} [Field Kbar]
  [Algebra V.functionField Kbar]

/-- ★**素因子 `x`(`V[L]` 上)が素因子 `y`(`V` 上)に写る** ——
`x` の像が `y` の閉包に入ること。

★★原文の「[possibly subvarieties of codimension ≥1 of]」がここに効く ——
**余次元は保たれなくてよい**ので、`=` ではなく**閉包への所属**で書く。 -/
def MapsIntoDivisor {L : FinSub V.functionField Kbar}
    (x : PrimeDivisorPt (normObj V L)) (y : PrimeDivisorPt V) : Prop :=
  (normDown V L).base x.1 ∈ closure ({y.1} : Set V)

/-- ★★★**`D_L`** —— `D_K` の素因子に写る `V[L]` の素因子の全体。 -/
def DLSet (DK : Set (PrimeDivisorPt V)) (L : FinSub V.functionField Kbar) :
    Set (PrimeDivisorPt (normObj V L)) :=
  {x | ∃ y ∈ DK, MapsIntoDivisor V x y}

theorem mem_DLSet_iff (DK : Set (PrimeDivisorPt V)) (L : FinSub V.functionField Kbar)
    (x : PrimeDivisorPt (normObj V L)) :
    x ∈ DLSet V DK L ↔ ∃ y ∈ DK, MapsIntoDivisor V x y := Iff.rfl

/-- ★`D_K` を大きくすれば `D_L` も大きくなる。 -/
theorem DLSet_mono {DK DK' : Set (PrimeDivisorPt V)} (h : DK ⊆ DK')
    (L : FinSub V.functionField Kbar) : DLSet V DK L ⊆ DLSet V DK' L := by
  rintro x ⟨y, hy, hxy⟩
  exact ⟨y, h hy, hxy⟩

/-- ★`D_K = ∅` なら `D_L = ∅`。 -/
@[simp] theorem DLSet_empty (L : FinSub V.functionField Kbar) :
    DLSet V (∅ : Set (PrimeDivisorPt V)) L = ∅ := by
  ext x
  simp [DLSet]

/-! ### ★出典の紐付け -/

/-- ★★★locator —— `Example 6.1` の `D_L`。 -/
def DLSet.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 109,
    item := "Example 6.1 — D_L(D_K の素因子に写る V[L] の素因子)",
    sectionId := "frdi-example-6-1" }

end ABC3.Found.Divisor
