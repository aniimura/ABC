import ABC3.Found.GenEll.NorthcottDegree
import ABC3.Found.GenEll.HeightConjugate
import Mathlib.FieldTheory.IntermediateField.Adjoin.Basic

/-!
# [GenEll] Proposition 1.4, (iv) —— **古典的 Northcott の定理**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.6。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

## ★★これが原文の `X(ℚ̄)^{≤d}` の形である

> **次数が `d` 以下かつ高さが `B` 以下の代数的数は有限個**

★体は固定されていない——`[F:ℚ] ≤ d` なる数体すべてを走る。
**次数 `≤ d` の数体は無限個ある**ので、固定体の結果を並べても届かない。

## ★機構 —— 3 つの部品を繋ぐだけ

1. `HeightConjugate.lean`: `H(x) ≤ B` から**共役 ≤ B** と**分母 `m ≤ B`**
2. `m·α` は**代数的整数**で、共役は `≤ B²`、次数は `≤ d`
3. `NorthcottDegree.lean`: そういう代数的整数は**有限個**
4. `α = (m·α)/m` で、`m` も有限個(`m ≤ B`)

★★**新しい数学は 1 つも無い。** 本日作った部品を繋いだだけである。

## ★何を主張していないか

原文 `Proposition 1.4, (iv)` は `X` 上の**算術直線束の高さ** `ht_L` についての主張で、
そちらには `Definition 1.1`(→ `X^arc`、複素解析空間)が要る。
★★本定理が取ったのは **`X = ℙ¹` の場合の高さ**についての Northcott 性である。
`.src` は条つきである。
-/

namespace ABC3.Found.GenEll

open NumberField Polynomial IntermediateField

/-- ★★★**古典的 Northcott の定理**。

> `ℚ(α)` の次数が `d` 以下、かつ `H_{ℚ(α)}(α) ≤ B` なる `α ∈ ℂ` は有限個。

★★体を固定していない——これが原文の `X(ℚ̄)^{≤d}` の形である。 -/
theorem finite_of_finrank_le_of_mulHeight₁_le (d : ℕ) (B : ℝ) :
    {α : ℂ | ∃ h : IsIntegral ℚ α,
      Module.finrank ℚ ℚ⟮α⟯ ≤ d ∧
      (haveI := IntermediateField.adjoin.finiteDimensional h
       haveI : NumberField ℚ⟮α⟯ := ⟨⟩
       Height.mulHeight₁ (⟨α, IntermediateField.mem_adjoin_simple_self ℚ α⟩ : ℚ⟮α⟯)
         ≤ B)}.Finite := by
  classical
  set S : Set ℂ := {β : ℂ | IsIntegral ℤ β ∧ (minpoly ℤ β).natDegree ≤ d ∧
      ∀ z ∈ ((minpoly ℤ β).map (algebraMap ℤ ℂ)).roots, ‖z‖ ≤ B * B} with hSdef
  have hSfin : S.Finite := finite_isIntegral_natDegree_le_norm_le d (B * B)
  refine Set.Finite.subset (Set.Finite.insert 0 (Set.Finite.biUnion
    (Set.finite_Icc 1 (⌊B⌋₊))
    (fun m _ => hSfin.image (fun β : ℂ => β / (m : ℂ))))) ?_
  rintro α ⟨h, hdeg, hB⟩
  rcases eq_or_ne α 0 with rfl | hα0
  · exact Set.mem_insert _ _
  refine Set.mem_insert_of_mem _ ?_
  haveI hfd := IntermediateField.adjoin.finiteDimensional h
  haveI : NumberField ℚ⟮α⟯ := ⟨⟩
  set x : ℚ⟮α⟯ := ⟨α, IntermediateField.mem_adjoin_simple_self ℚ α⟩ with hxdef
  have hxC : algebraMap ℚ⟮α⟯ ℂ x = α := rfl
  have hx0 : x ≠ 0 := by
    intro hc
    exact hα0 (by rw [← hxC, hc, map_zero])
  -- 分母を取る
  obtain ⟨m, hm0, hmB, hmint⟩ := exists_denom_isIntegral x hx0 hB
  have hB1 : (1 : ℝ) ≤ B := le_trans (by exact_mod_cast hm0) hmB
  set β : ℂ := (m : ℂ) * α with hβdef
  -- `β` を `ℚ⟮α⟯` の元として持つ
  set y : ℚ⟮α⟯ := (m : ℚ⟮α⟯) * x with hydef
  have hyC : algebraMap ℚ⟮α⟯ ℂ y = β := by
    rw [hydef, map_mul, map_natCast, hxC]
  -- ★共役の評価
  have hconj : ∀ z ∈ (minpoly ℚ β).rootSet ℂ, ‖z‖ ≤ B * B := by
    rw [← hyC]
    refine norm_le_of_mem_rootSet_minpoly' y (fun φ => ?_)
    have hφx : ‖φ x‖ ≤ B := by
      calc ‖φ x‖ = (InfinitePlace.mk φ) x := (InfinitePlace.apply φ x).symm
        _ ≤ Height.mulHeight₁ x := infinitePlace_le_mulHeight₁ _ x
        _ ≤ B := hB
    have : ‖φ y‖ = (m : ℝ) * ‖φ x‖ := by
      rw [hydef, map_mul, map_natCast, norm_mul, Complex.norm_natCast]
    rw [this]
    exact mul_le_mul hmB hφx (norm_nonneg _) (le_trans zero_le_one hB1)
  -- ★`β` は代数的整数
  have hβint : IsIntegral ℤ β := hmint
  -- ★次数の評価 —— `β` は `ℚ⟮α⟯` の元なので、最小多項式の次数は `[ℚ(α):ℚ]` 以下
  have hdegβ : (minpoly ℤ β).natDegree ≤ d := by
    have h1 : minpoly ℚ β = (minpoly ℤ β).map (algebraMap ℤ ℚ) :=
      minpoly.isIntegrallyClosed_eq_field_fractions' ℚ hβint
    have h2 : (minpoly ℤ β).natDegree = (minpoly ℚ β).natDegree := by
      rw [h1, (minpoly.monic hβint).natDegree_map]
    have h3 : minpoly ℚ β = minpoly ℚ y := by
      rw [← hyC]
      exact minpoly.algebraMap_eq (algebraMap ℚ⟮α⟯ ℂ).injective y
    have h4 : (minpoly ℚ y).natDegree ≤ Module.finrank ℚ ℚ⟮α⟯ := minpoly.natDegree_le y
    rw [h2, h3]
    exact le_trans h4 hdeg
  -- ★根の集合の言い換え
  have hroots : ∀ z ∈ ((minpoly ℤ β).map (algebraMap ℤ ℂ)).roots, ‖z‖ ≤ B * B := by
    intro z hz
    refine hconj z ?_
    have h1 : minpoly ℚ β = (minpoly ℤ β).map (algebraMap ℤ ℚ) :=
      minpoly.isIntegrallyClosed_eq_field_fractions' ℚ hβint
    have h5 : (minpoly ℚ β).map (algebraMap ℚ ℂ) = (minpoly ℤ β).map (algebraMap ℤ ℂ) := by
      rw [h1, Polynomial.map_map, ← IsScalarTower.algebraMap_eq]
    rw [Polynomial.rootSet_def, Finset.mem_coe, Multiset.mem_toFinset, Polynomial.aroots, h5]
    exact hz
  simp only [Set.mem_iUnion, Set.mem_image, exists_prop]
  refine ⟨m, ⟨hm0, Nat.le_floor hmB⟩, β, ⟨hβint, hdegβ, hroots⟩, ?_⟩
  rw [hβdef]
  have hmC : (m : ℂ) ≠ 0 := Nat.cast_ne_zero.2 hm0.ne'
  field_simp

/-! ## ★出典の紐付け(`.src`) -/

def finite_of_finrank_le_of_mulHeight₁_le.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(古典的 Northcott——`X = ℙ¹` の高さの場合)",
    sectionId := "genell-prop-1-4" }

end ABC3.Found.GenEll
