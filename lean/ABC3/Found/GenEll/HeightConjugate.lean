import ABC3.Found.GenEll.DenominatorBound
import Mathlib.FieldTheory.IntermediateField.Adjoin.Basic

/-!
# [GenEll] Proposition 1.4, (iv) —— 高さから**共役と分母**を抑える(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.6。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

## ★★これが古典的 Northcott の最後の橋である

`Found/GenEll/NorthcottDegree.lean` は
**「次数 `≤ d`・共役 `≤ B` の代数的整数は有限個」**(体を固定しない)を取った。
残るのは「**高さ**が有界」から「**共役が有界**」「**分母が有界**」を出すことである。

★本ファイルがそれを取る:

> `H_K(x) ≤ B` なら
> - `x` の**共役はすべて `B` 以下**(`|φ(x)| ≤ H_K(x)`)
> - **`B` 以下の自然数 `m`** があって `m·x` は代数的整数

★★どちらも既に持っている道具の言い換えである:
- 共役: `NorthcottNF.lean` の `infinitePlace_le_mulHeight₁`
- 分母: `DenominatorBound.lean` の `denIdeal` と `N(𝔡_x) ≤ H(x)`

★**共役の集合が `{φ(x) | φ : K →+* ℂ}` であること**は mathlib の
`NumberField.Embeddings.range_eval_eq_rootSet_minpoly` である。
`x : K` の最小多項式と、その `ℂ` での像 `α` の最小多項式が一致することは
`minpoly.algebraMap_eq`(単射な代数写像で最小多項式は変わらない)。
-/

namespace ABC3.Found.GenEll

open NumberField Polynomial

section Conjugates

variable {K : Type*} [Field K] [NumberField K] [Algebra K ℂ]

/-- ★**共役はすべて高さ以下** —— `α = x` の `ℂ` での像とすると、
`minpoly ℚ α` の根はすべて絶対値が `H_K(x)` 以下。

★根の集合が `{φ(x)}` であることは mathlib
(`range_eval_eq_rootSet_minpoly`)、
各 `‖φ(x)‖ ≤ H_K(x)` は `NorthcottNF.lean` の `infinitePlace_le_mulHeight₁`。 -/
theorem norm_le_of_mem_rootSet_minpoly' (x : K) {B : ℝ}
    (hB : ∀ φ : K →+* ℂ, ‖φ x‖ ≤ B) :
    ∀ z ∈ (minpoly ℚ (algebraMap K ℂ x)).rootSet ℂ, ‖z‖ ≤ B := by
  intro z hz
  have hmin : minpoly ℚ (algebraMap K ℂ x) = minpoly ℚ x :=
    minpoly.algebraMap_eq (algebraMap K ℂ).injective x
  rw [hmin, ← NumberField.Embeddings.range_eval_eq_rootSet_minpoly K ℂ x] at hz
  obtain ⟨φ, rfl⟩ := hz
  exact hB φ

/-- ★**高さ版** —— 上の系。 -/
theorem norm_le_of_mem_rootSet_minpoly (x : K) {B : ℝ}
    (hB : Height.mulHeight₁ x ≤ B) :
    ∀ z ∈ (minpoly ℚ (algebraMap K ℂ x)).rootSet ℂ, ‖z‖ ≤ B :=
  norm_le_of_mem_rootSet_minpoly' x (fun φ => by
    calc ‖φ x‖ = (InfinitePlace.mk φ) x := (InfinitePlace.apply φ x).symm
      _ ≤ Height.mulHeight₁ x := infinitePlace_le_mulHeight₁ _ x
      _ ≤ B := hB)

open scoped Classical in
/-- ★★**分母は高さ以下** —— `H_K(x) ≤ B` なら `B` 以下の自然数 `m` があって
`m·x` の `ℂ` での像は**代数的整数**。

★`DenominatorBound.lean` の `N(𝔡_x) ≤ H(x)` と `Ideal.absNorm I ∈ I` そのもの。 -/
theorem exists_denom_isIntegral (x : K) (hx : x ≠ 0) {B : ℝ}
    (hB : Height.mulHeight₁ x ≤ B) :
    ∃ m : ℕ, 0 < m ∧ (m : ℝ) ≤ B ∧ IsIntegral ℤ ((m : ℂ) * algebraMap K ℂ x) := by
  classical
  set u : Kˣ := Units.mk0 x hx with hu
  set m : ℕ := Ideal.absNorm (denIdeal u) with hm
  have hmle : (m : ℝ) ≤ B := le_trans (absNorm_denIdeal_le_mulHeight₁ u) hB
  have hmne : m ≠ 0 := by
    rw [hm]
    simpa [Ideal.absNorm_eq_zero_iff] using denIdeal_ne_bot u
  have hmem : ((m : ℕ) : 𝓞 K) ∈ denIdeal u := Ideal.absNorm_mem (denIdeal u)
  have hmne' : ((m : ℕ) : 𝓞 K) ≠ 0 := Nat.cast_ne_zero.2 hmne
  obtain ⟨z, hz⟩ := mul_isIntegral_of_mem_denIdeal (x := u) hmne' hmem
  rw [map_natCast] at hz
  refine ⟨m, Nat.pos_of_ne_zero hmne, hmle, ?_⟩
  -- `algebraMap (𝓞 K) K z = (m : K) * x` を `ℂ` へ運ぶ
  have hzK : algebraMap (𝓞 K) K z = (m : K) * x := hz
  have hzC : algebraMap K ℂ (algebraMap (𝓞 K) K z) = (m : ℂ) * algebraMap K ℂ x := by
    rw [hzK, map_mul, map_natCast]
  rw [← hzC]
  exact (NumberField.RingOfIntegers.isIntegral_coe z).map
    ((algebraMap K ℂ).toIntAlgHom)

end Conjugates

/-! ## ★出典の紐付け(`.src`) -/

def norm_le_of_mem_rootSet_minpoly.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(高さ → 共役・分母の評価)",
    sectionId := "genell-prop-1-4" }

end ABC3.Found.GenEll
