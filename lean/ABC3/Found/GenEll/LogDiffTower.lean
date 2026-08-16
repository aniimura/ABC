import ABC3.Found.GenEll.LogDiffValue

/-!
# [GenEll] Definition 1.5, (iii) —— **なぜ minimal field of definition なのか**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.8。

原文 (GenEll p.8):
> (iii) Let x ∈ X(F ) ⊆ X(Q), where F is a minimal field of definition of x. Then the diﬀerent ideal of F determines an eﬀective arithmetic divisor

## ★★原文が "minimal" と書く理由を、定理にする

`Found/GenEll/LogDiffValue.lean` で

> `log-diff = log|disc F| / [F:ℚ]`

を取った。★しかしそれは**`F` を 1 つ選んだときの値**であって、
「なぜ最小の `F` を選ぶのか」には答えていなかった。

★★本ファイルが答える:

> **`F ⊆ K` なら `log-diff(F) ≤ log-diff(K)`**

——正規化した log-different は**体を大きくすると増える**。
だから `x` の定義体として**最小のもの**を取ると `log-diff_X(x)` が**最小**になり、
それが原文の `x` に対する不変量として意味を持つ。

★**"minimal" は恣意的な選択ではなく、最小値を取る操作である。**

## ★機構 —— 差積の tower 公式

mathlib に **`NumberField.natAbs_discr_eq_absNorm_differentIdeal_mul_natAbs_discr_pow`** が
**ある**(2026-08-17 実測):

> `|disc K| = N(𝔡_{K/F}) · |disc F|^{[K:F]}`

対数を取ると `log|disc K| = log N(𝔡_{K/F}) + [K:F]·log|disc F|`。
`[K:ℚ] = [F:ℚ]·[K:F]` で割ると

> `log-diff(K) = log N(𝔡_{K/F})/[K:ℚ] + log-diff(F)`

★**差は `log N(𝔡_{K/F})/[K:ℚ] ≥ 0` そのもの**である。
不等式ではなく**等式**が取れるので、増加分が何であるかまで分かる。
`N(𝔡) ≥ 1` は `𝔡 ≠ 0`(`differentIdeal_ne_bot`)から出る。

★★**`K/F` が不分岐なら `𝔡_{K/F} = (1)` で差は 0** ——
これも上の等式から読める(本ファイルでは取らない)。
-/

namespace ABC3.Found.GenEll

open NumberField IsDedekindDomain

section Tower

variable (F K : Type*) [Field F] [NumberField F] [Field K] [NumberField K] [Algebra F K]

/-- ★**差積の tower 公式の対数版**。

> `log|disc K| = log N(𝔡_{K/F}) + [K:F]·log|disc F|`

★mathlib の `natAbs_discr_eq_absNorm_differentIdeal_mul_natAbs_discr_pow` の対数を取っただけ。 -/
theorem log_discr_tower :
    Real.log ((NumberField.discr K).natAbs)
      = Real.log (Ideal.absNorm (differentIdeal (𝓞 F) (𝓞 K)))
        + (Module.finrank F K : ℝ) * Real.log ((NumberField.discr F).natAbs) := by
  have hd : (NumberField.discr K).natAbs
      = Ideal.absNorm (differentIdeal (𝓞 F) (𝓞 K))
        * (NumberField.discr F).natAbs ^ Module.finrank F K :=
    NumberField.natAbs_discr_eq_absNorm_differentIdeal_mul_natAbs_discr_pow F (𝓞 F) K (𝓞 K)
  have hNne : (Ideal.absNorm (differentIdeal (𝓞 F) (𝓞 K)) : ℝ) ≠ 0 := by
    have : differentIdeal (𝓞 F) (𝓞 K) ≠ 0 := differentIdeal_ne_bot
    have := (Ideal.absNorm_eq_zero_iff (I := differentIdeal (𝓞 F) (𝓞 K))).not.2 this
    exact_mod_cast this
  have hdFne : ((NumberField.discr F).natAbs : ℝ) ≠ 0 := by
    have : (NumberField.discr F) ≠ 0 := NumberField.discr_ne_zero F
    have : (NumberField.discr F).natAbs ≠ 0 := Int.natAbs_ne_zero.2 this
    exact_mod_cast this
  rw [hd]
  push_cast
  rw [Real.log_mul hNne (pow_ne_zero _ hdFne), Real.log_pow]

/-- ★**`N(𝔡_{K/F}) ≥ 1`**、したがって `log N(𝔡_{K/F}) ≥ 0`。 -/
theorem log_absNorm_differentIdeal_nonneg :
    0 ≤ Real.log (Ideal.absNorm (differentIdeal (𝓞 F) (𝓞 K))) := by
  have hne : differentIdeal (𝓞 F) (𝓞 K) ≠ 0 := differentIdeal_ne_bot
  have h0 : Ideal.absNorm (differentIdeal (𝓞 F) (𝓞 K)) ≠ 0 :=
    (Ideal.absNorm_eq_zero_iff (I := differentIdeal (𝓞 F) (𝓞 K))).not.2 hne
  have h1n : 1 ≤ Ideal.absNorm (differentIdeal (𝓞 F) (𝓞 K)) := Nat.one_le_iff_ne_zero.2 h0
  have h1 : (1 : ℝ) ≤ (Ideal.absNorm (differentIdeal (𝓞 F) (𝓞 K)) : ℝ) := by
    exact_mod_cast h1n
  exact Real.log_nonneg h1

/-- ★★**正規化 log-different の増加分は `log N(𝔡_{K/F})/[K:ℚ]` ちょうどである**。

> `log-diff(K) = log N(𝔡_{K/F})/[K:ℚ] + log-diff(F)`

★不等式ではなく**等式**を取る。増加分が何であるかまで見える。 -/
theorem logDiffOfField_tower :
    logDiffOfField K
      = Real.log (Ideal.absNorm (differentIdeal (𝓞 F) (𝓞 K)))
          / (Module.finrank ℚ K : ℝ)
        + logDiffOfField F := by
  haveI : FiniteDimensional F K := Module.Finite.of_restrictScalars_finite ℚ F K
  have htower : Module.finrank ℚ F * Module.finrank F K = Module.finrank ℚ K :=
    Module.finrank_mul_finrank ℚ F K
  have hFK : (0 : ℝ) < (Module.finrank F K : ℝ) := by
    exact_mod_cast Module.finrank_pos (R := F) (M := K)
  have hF : (0 : ℝ) < (Module.finrank ℚ F : ℝ) := by
    exact_mod_cast Module.finrank_pos (R := ℚ) (M := F)
  rw [logDiffOfField_eq, logDiffOfField_eq, log_discr_tower F K, ← htower]
  push_cast
  field_simp

/-- ★★★**正規化 log-different は体を大きくすると増える**。

> **`F ⊆ K` なら `log-diff(F) ≤ log-diff(K)`**

★★これが原文 Definition 1.5, (iii) の "minimal field of deﬁnition" の**理由**である。
`x` の定義体はいくらでも大きく取れるが、`log-diff` はそのたびに増えるので、
**最小の定義体でだけ `x` の不変量になる**。

原文 (GenEll p.8):
> (iii) Let x ∈ X(F ) ⊆ X(Q), where F is a minimal field of definition of x. Then the diﬀerent ideal of F determines an eﬀective arithmetic divisor -/
theorem logDiffOfField_le : logDiffOfField F ≤ logDiffOfField K := by
  rw [logDiffOfField_tower F K]
  have hK : (0 : ℝ) < (Module.finrank ℚ K : ℝ) := by
    exact_mod_cast Module.finrank_pos (R := ℚ) (M := K)
  have := log_absNorm_differentIdeal_nonneg F K
  have hdiv : 0 ≤ Real.log (Ideal.absNorm (differentIdeal (𝓞 F) (𝓞 K)))
      / (Module.finrank ℚ K : ℝ) := div_nonneg this hK.le
  linarith

end Tower

/-! ## ★出典の紐付け(`.src`) -/

def logDiffOfField_le.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 8, item := "Definition 1.5, (iii)",
    sectionId := "genell-def-1-5" }

end ABC3.Found.GenEll
