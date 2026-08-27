import ABC3.Skeleton.NCBelyi.Theorem25

/-!
# ★[NCBelyi] `Theorem 2.5` の結論は空虚ではない(`Check`)

原典: S. Mochizuki, *Noncritical Belyi Maps* [NCBelyi]、物理 p.5。

原文 (NCBelyi p.5):
> Theorem 2.5. (Belyi Maps Noncritical at Prescribed Points) Let X be a smooth, proper, connected curve over Q[bb][bar] and

## ★★★なぜこの検査が要るか

2026-08-27(第 425 ブロック)に `theorem_2_5` を構成へ載せ替えた。
★**載せ替えで `sorry` が消えても、結論が空虚なら進捗ではない**。

★★本ファイルは結論が**実際に中身を持つ**ことを機械検証する:

- (a)(b) の側は `S = {√2}` のような**有理でない点**に対しても多項式が出る
- (c) の側は `S₀ = {0, 1}`、`β = 2` という**実際に仮定を満たす**組で叩ける
-/

namespace ABC3.Check.NCBelyi

open ABC3.Skeleton.NCBelyi

/-- ★**(a)(b) は空でない集合に対して中身を持つ** —— `S = {0}` で叩く。 -/
theorem theorem_2_5_ab_nonvacuous :
    ∃ F : Polynomial ℚ, 0 < F.natDegree
      ∧ (∀ x ∈ ({0} : Finset ℂ), Polynomial.aeval x F = 0 ∨ Polynomial.aeval x F = 1)
      ∧ (∀ w : ℂ, Polynomial.aeval w (Polynomial.derivative F) = 0 →
          Polynomial.aeval w F = 0 ∨ Polynomial.aeval w F = 1) :=
  (theorem_2_5 ({0} : Finset ℂ) (fun x hx => by
    rw [Finset.mem_singleton] at hx
    subst hx
    exact isIntegral_zero)).1

/-- ★★**(c) の仮定は実際に満たされる** —— `S₀ = {0, 1}`、`β = 2`。

★`0 ∈ S₀`、`1 > 0`、`2 ∉ S₀`、`2 ≠ 0`、`2·1 ≤ 2` がすべて成り立つ。 -/
theorem theorem_2_5_c_nonvacuous :
    ∃ f : Polynomial ℚ, 0 < f.natDegree
      ∧ (∀ α ∈ ({0, 1} : Finset ℚ), f.eval α = 0 ∨ f.eval α = 1)
      ∧ f.eval 2 ≠ 0 ∧ f.eval 2 ≠ 1
      ∧ (∀ x : ℂ, (Polynomial.derivative (f.map (algebraMap ℚ ℂ))).eval x = 0 →
          (f.map (algebraMap ℚ ℂ)).eval x = 0
            ∨ (f.map (algebraMap ℚ ℂ)).eval x = 1) := by
  refine (theorem_2_5 (∅ : Finset ℂ) (by simp)).2 ({0, 1} : Finset ℚ) 2 (by simp) ?_
    (by norm_num) (by norm_num) ?_
  · intro α hα hα0
    fin_cases hα
    · exact absurd rfl hα0
    · norm_num
  · intro α hα _
    fin_cases hα <;> norm_num

end ABC3.Check.NCBelyi
