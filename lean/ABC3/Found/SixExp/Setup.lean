/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.SixExp.SiegelConst

/-!
# 背理法の入口 —— 「6 個が代数的」から数体を作る

★`ResearchPaper/frdi-decomposition.json` の `sixexp` チェーン。
最終目標は [FrdI] `Lemma 6.5, (ii)`(six exponentials theorem)。

## ★中身

six exponentials は背理法で示す。「`exp(x_j y_i)`(6 個)がすべて代数的」から出発して

1. **数体を作る** —— `sixField v = ℚ(v₀₀, …, v₁₂) ⊆ ℂ`。
   有限個の代数的元で生成するので有限次(`IntermediateField.finiteDimensional_adjoin`)、
   したがって `NumberField`(`sixField_numberField`)。
   複素数への埋め込み `σ` は包含そのもの(`sixFieldHom`)。
2. **共通分母を取る** —— `b ∈ 𝓞_K \ {0}` で `b·v j i ∈ 𝓞_K`
   (`exists_common_denom`。`IsLocalization.exist_integer_multiples`)。
   ★これで `DenomClear.lean` の `denomBase` が使える。
3. **鋭い数え上げ** —— `gap_of_base`。

## ★`gap_of_base` がなぜ要るか

`Counting.lean` の `exists_gap_index` は「`2^N > A^b`」を要求するので、
★**`b` が `N` とともに増える場合(実際そうなる)に使えない**。
`gap_of_base` は**初段の不等式だけ**を要求し、あとは

  `(m+1)³ − m³ = 3m² + 3m + 1 ≥ 3N²`

で押していく。★この形なら「`L' ≈ N^{3/2}` を取ると初段が `N³log2 > O(N^{5/2})` で成り立つ」
という**本来の数え上げ**がそのまま渡せる。
-/

namespace ABC3.Found.SixExp

open NumberField IntermediateField

/-! ## ★1. 鋭い数え上げ -/

/-- ★★★**鋭い数え上げ** —— 初段で成り立ち、`Z` が `2^{3N²}` を超えなければ、
**すべての段で成り立つ**。 -/
theorem gap_of_base {K Z : ℝ} (_hK : 0 < K) (hZ1 : 1 ≤ Z) {N : ℕ}
    (hbase : K * Z ^ N < 2 ^ (N ^ 3))
    (hZ : Z ≤ 2 ^ (3 * N ^ 2)) :
    ∀ m, N ≤ m → K * Z ^ m < 2 ^ (m ^ 3) := by
  intro m hm
  induction m, hm using Nat.le_induction with
  | base => exact hbase
  | succ n hn ih =>
    have hZ0 : (0:ℝ) < Z := lt_of_lt_of_le zero_lt_one hZ1
    have hexp : (3 : ℕ) * N ^ 2 ≤ (n + 1) ^ 3 - n ^ 3 := by
      have h1 : (n + 1) ^ 3 - n ^ 3 = 3 * n ^ 2 + 3 * n + 1 := by ring_nf; omega
      rw [h1]
      have h2 : N ^ 2 ≤ n ^ 2 := Nat.pow_le_pow_left hn 2
      omega
    have hZle : Z ≤ 2 ^ ((n + 1) ^ 3 - n ^ 3) :=
      le_trans hZ (pow_le_pow_right₀ (by norm_num) hexp)
    calc K * Z ^ (n + 1) = (K * Z ^ n) * Z := by ring
      _ < 2 ^ (n ^ 3) * Z := mul_lt_mul_of_pos_right ih hZ0
      _ ≤ 2 ^ (n ^ 3) * 2 ^ ((n + 1) ^ 3 - n ^ 3) :=
          mul_le_mul_of_nonneg_left hZle (by positivity)
      _ = 2 ^ ((n + 1) ^ 3) := by
          rw [← pow_add]
          congr 1
          have h3 : n ^ 3 ≤ (n + 1) ^ 3 := Nat.pow_le_pow_left (by omega) 3
          omega

/-! ## ★2. 6 個の代数的数が生成する数体 -/

/-- ★6 個の複素数が生成する部分体。 -/
noncomputable def sixField (v : Fin 2 → Fin 3 → ℂ) : IntermediateField ℚ ℂ :=
  IntermediateField.adjoin ℚ (Set.range (fun p : Fin 2 × Fin 3 => v p.1 p.2))

theorem mem_sixField (v : Fin 2 → Fin 3 → ℂ) (j : Fin 2) (i : Fin 3) : v j i ∈ sixField v :=
  IntermediateField.subset_adjoin ℚ _ ⟨(j, i), rfl⟩

/-- ★生成元(数体の元として)。 -/
noncomputable def sixFieldElt (v : Fin 2 → Fin 3 → ℂ) (j : Fin 2) (i : Fin 3) :
    ↥(sixField v) := ⟨v j i, mem_sixField v j i⟩

/-- ★複素数への埋め込み —— 包含そのもの。 -/
noncomputable def sixFieldHom (v : Fin 2 → Fin 3 → ℂ) : ↥(sixField v) →+* ℂ :=
  (sixField v).val.toRingHom

@[simp] theorem sixFieldHom_elt (v : Fin 2 → Fin 3 → ℂ) (j : Fin 2) (i : Fin 3) :
    sixFieldHom v (sixFieldElt v j i) = v j i := rfl

theorem sixFieldHom_injective (v : Fin 2 → Fin 3 → ℂ) :
    Function.Injective (sixFieldHom v) := (sixField v).val.injective

/-- ★★有限個の代数的元で生成するので**有限次**。 -/
theorem sixField_finiteDimensional (v : Fin 2 → Fin 3 → ℂ)
    (hv : ∀ j i, IsAlgebraic ℚ (v j i)) : FiniteDimensional ℚ ↥(sixField v) := by
  classical
  haveI : Finite ↥(Set.range (fun p : Fin 2 × Fin 3 => v p.1 p.2)) :=
    (Set.finite_range _).to_subtype
  refine IntermediateField.finiteDimensional_adjoin ?_
  rintro _ ⟨p, rfl⟩
  exact (hv p.1 p.2).isIntegral

/-- ★★したがって**数体**。 -/
theorem sixField_numberField (v : Fin 2 → Fin 3 → ℂ)
    (hv : ∀ j i, IsAlgebraic ℚ (v j i)) : NumberField ↥(sixField v) :=
  haveI := sixField_finiteDimensional v hv
  ⟨⟩

/-! ## ★3. 共通分母 -/

/-- ★★**共通分母** —— 数体の有限個の元に共通の分母 `b ∈ 𝓞_K \ {0}` が取れる。

★これで `DenomClear.lean` の `denomBase` が使える。 -/
theorem exists_common_denom (K : Type*) [Field K] [NumberField K] (v : Fin 2 → Fin 3 → K) :
    ∃ b : 𝓞 K, b ≠ 0 ∧ ∀ j i, ∃ w : 𝓞 K, (w : K) = (b : K) * v j i := by
  classical
  obtain ⟨b, hb⟩ := IsLocalization.exist_integer_multiples (nonZeroDivisors (𝓞 K))
    (Finset.univ : Finset (Fin 2 × Fin 3)) (fun p => v p.1 p.2)
  refine ⟨(b : 𝓞 K), nonZeroDivisors.coe_ne_zero b, ?_⟩
  intro j i
  obtain ⟨w, hw⟩ := hb (j, i) (Finset.mem_univ _)
  refine ⟨w, ?_⟩
  have h1 : ((w : 𝓞 K) : K) = (algebraMap (𝓞 K) K) w := rfl
  rw [h1, hw]
  simp only [Algebra.smul_def]

end ABC3.Found.SixExp
