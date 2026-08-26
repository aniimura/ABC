/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Divisor.SchemeWeilOrd
import Mathlib.RingTheory.Ideal.Height
import Mathlib.RingTheory.Ideal.MinimalPrime.Basic

/-!
# `Ideal.height = 1` と `HeightOnePrime` の橋(鎖 `weil` の `div-finite` への接続)

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.109。

原文 (FrdI p.109):
> normal [geometrically integral] variety over a field k; K the function field of V ;

## ★★2 つの「高さ 1」を突き合わせる

我々は**高さ 1 の素イデアル**を 2 通りに書いてきた:

| 場所 | 綴り |
|---|---|
| `HeightOneDVR.lean`(環の層、`div(a)` を組んだ側) | `p ≠ ⊥ ∧ ∀ q 素, q ≠ ⊥ → q ≤ p → q = p` |
| `SchemeWeilOrd.lean`(スキームの層、`isCodimOnePt_spec_iff`) | `Ideal.height p = 1` |

★★**両者が同値であることを示す**。これで
「スキームの余次元 1 の点」→「`Ideal.height = 1`」→「`HeightOnePrime`」→
`divOfElem`(`HeightOneDVR.lean` の因子)がひと繋がりになる。

## ★★★中身

`Ideal.height_le_iff` は「`height p ≤ n ⟺ `p` の下の素 `q` はすべて `height q < n`」と言う。
`n = 1` では「`p` の下の素はすべて高さ `0`」、すなわち
`minimalPrimes R`(整域では `{⊥}`)の元である。
★下からの `1 ≤ height p` は `p ∉ minimalPrimes R`、すなわち `p ≠ ⊥` そのもの。

## ★本ファイルで閉じること

| 定理 | 中身 |
|---|---|
| `height_eq_one_iff` | ★★**2 つの「高さ 1」は同値** |
| `HeightOnePrime.ofHeightEqOne` | `height = 1` から `HeightOnePrime` を作る |
| `isCodimOnePt_spec_iff_heightOne` | ★`Spec R` の余次元 1 の点 ⟺ `HeightOnePrime` |
| `primeDivisorPtSpecEquiv` | `Spec R` の素因子 ≃ `HeightOnePrime R` |
| `weilDivOfFrac` / `weilDivOfFrac_mul` | ★★**`div : K^× → WeilDiv (Spec R)`**（鎖 `weil` の `weil-group` のアフィンの場合）|
-/

namespace ABC3.Found.Divisor

open AlgebraicGeometry Ideal

universe u

/-! ## ★1. 2 つの「高さ 1」の同値 -/

/-- ★★★**`Ideal.height p = 1` と「非零かつ下に非零素が無い」は同値**(整域)。

★`≤ 1` は `Ideal.height_le_iff`(下の素はすべて高さ `0` = 極小 = `⊥`)、
`1 ≤` は `Ideal.height_eq_zero_iff`(高さ `0` ⟺ 極小 ⟺ `⊥`)。 -/
theorem height_eq_one_iff {R : Type*} [CommRing R] [IsDomain R] (p : Ideal R) [p.IsPrime] :
    p.height = 1 ↔ p ≠ ⊥ ∧ ∀ q : Ideal R, q.IsPrime → q ≠ ⊥ → q ≤ p → q = p := by
  have hmin : minimalPrimes R = {(⊥ : Ideal R)} := IsDomain.minimalPrimes_eq_singleton_bot R
  constructor
  · intro hp
    refine ⟨?_, ?_⟩
    · intro hbot
      rw [hbot, Ideal.height_bot] at hp
      exact zero_ne_one hp
    · intro q hq hq0 hqp
      by_contra hne
      haveI := hq
      have hlt : q < p := lt_of_le_of_ne hqp hne
      have h1 := (Ideal.height_le_iff.mp hp.le) q hq hlt
      have h0 : q.height = 0 := by
        have h2 : q.height < 1 := by exact_mod_cast h1
        by_contra hne
        exact absurd h2 (not_lt.mpr (Order.one_le_iff_ne_zero.mpr hne))
      have : q ∈ minimalPrimes R := Ideal.height_eq_zero_iff.mp h0
      rw [hmin, Set.mem_singleton_iff] at this
      exact hq0 this
  · rintro ⟨hp0, hmn⟩
    refine le_antisymm ?_ ?_
    · refine Ideal.height_le_iff.mpr (fun q hq hlt => ?_)
      haveI := hq
      have hq0 : q = ⊥ := by
        by_contra hne
        exact absurd (hmn q hq hne hlt.le) (ne_of_lt hlt)
      rw [hq0, Ideal.height_bot]
      exact zero_lt_one
    · rw [Order.one_le_iff_ne_zero, Ne, Ideal.height_eq_zero_iff, hmin,
        Set.mem_singleton_iff]
      exact hp0

/-! ## ★2. `HeightOnePrime` を作る -/

/-- ★**`Ideal.height = 1` から `HeightOnePrime` を作る**。 -/
def HeightOnePrime.ofHeightEqOne {R : Type*} [CommRing R] [IsDomain R] (p : Ideal R)
    (hpp : p.IsPrime) (hp : p.height = 1) : HeightOnePrime R :=
  haveI := hpp
  { asIdeal := p
    isPrime := hpp
    ne_bot := ((height_eq_one_iff p).mp hp).1
    minimal := ((height_eq_one_iff p).mp hp).2 }

/-- ★逆向き —— `HeightOnePrime` の高さは `1`。 -/
theorem HeightOnePrime.height_eq_one {R : Type*} [CommRing R] [IsDomain R]
    (v : HeightOnePrime R) : v.asIdeal.height = 1 :=
  haveI := v.isPrime
  (height_eq_one_iff v.asIdeal).mpr ⟨v.ne_bot, v.minimal⟩

/-! ## ★3. スキームの側との接続 -/

/-- ★★★**`Spec R` の余次元 1 の点は、ちょうど `HeightOnePrime R` の元**。

★★これで「スキームの素因子」と `HeightOneDVR.lean` の `divOfElem` がひと繋がりになる。 -/
theorem isCodimOnePt_spec_iff_heightOne (R : CommRingCat.{u}) [IsDomain R] (x : Spec R) :
    IsCodimOnePt (Spec R) x ↔
      (x.asIdeal ≠ ⊥ ∧ ∀ q : Ideal R, q.IsPrime → q ≠ ⊥ → q ≤ x.asIdeal → q = x.asIdeal) := by
  rw [isCodimOnePt_spec_iff]
  exact height_eq_one_iff x.asIdeal

/-! ## ★4. `Spec R` の Weil 因子 -/

/-- ★★**`Spec R` の素因子は `HeightOnePrime R` と 1 対 1**。 -/
def primeDivisorPtSpecEquiv (R : CommRingCat.{u}) [IsDomain R] :
    PrimeDivisorPt (Spec R) ≃ HeightOnePrime R where
  toFun x := HeightOnePrime.ofHeightEqOne x.1.asIdeal x.1.isPrime
    ((isCodimOnePt_spec_iff R x.1).mp x.2)
  invFun v := ⟨⟨v.asIdeal, v.isPrime⟩,
    (isCodimOnePt_spec_iff R ⟨v.asIdeal, v.isPrime⟩).mpr v.height_eq_one⟩
  left_inv _ := rfl
  right_inv _ := rfl

open scoped Classical in
/-- ★★★★★**`Spec R` の有理関数の Weil 因子** ——
`HeightOneDVR.lean` の `divOfFrac` を素因子の対応で送ったもの。

★★これで鎖 `weil` の `weil-group` が【アフィンの場合に限って】閉じる。 -/
noncomputable def weilDivOfFrac (R : CommRingCat.{u}) [IsDomain R] [IsNoetherianRing R]
    [IsIntegrallyClosed R] (x : FractionRing R) (hx : x ≠ 0) : WeilDiv (Spec R) :=
  Finsupp.equivMapDomain (primeDivisorPtSpecEquiv R).symm (divOfFrac x hx)

theorem weilDivOfFrac_apply (R : CommRingCat.{u}) [IsDomain R] [IsNoetherianRing R]
    [IsIntegrallyClosed R] (x : FractionRing R) (hx : x ≠ 0) (v : PrimeDivisorPt (Spec R)) :
    weilDivOfFrac R x hx v = divOfFrac x hx (primeDivisorPtSpecEquiv R v) := rfl

/-- ★★★**`div` は乗法を加法へ移す**（`Spec R` の層）。 -/
theorem weilDivOfFrac_mul (R : CommRingCat.{u}) [IsDomain R] [IsNoetherianRing R]
    [IsIntegrallyClosed R] {x y : FractionRing R} (hx : x ≠ 0) (hy : y ≠ 0) :
    weilDivOfFrac R (x * y) (mul_ne_zero hx hy)
      = weilDivOfFrac R x hx + weilDivOfFrac R y hy := by
  refine Finsupp.ext fun v => ?_
  rw [weilDivOfFrac_apply, Finsupp.add_apply, weilDivOfFrac_apply, weilDivOfFrac_apply,
    divOfFrac_mul hx hy, Finsupp.add_apply]

/-! ### ★出典の紐付け -/

/-- ★★★locator —— `Example 6.1` の「余次元 1 の点 = 高さ 1 の素イデアル」の 2 つの綴りの一致。 -/
def height_eq_one_iff.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 109,
    item := "Example 6.1 — 高さ 1 の素イデアルの 2 つの綴りの一致",
    sectionId := "frdi-example-6-1" }

/-- ★★★locator —— `Example 6.1` の `div : K^× → WeilDiv`（アフィンの場合）。 -/
def weilDivOfFrac.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 109,
    item := "Example 6.1 — 有理関数の Weil 因子 div(f)（Spec R の層）",
    sectionId := "frdi-example-6-1" }

end ABC3.Found.Divisor
