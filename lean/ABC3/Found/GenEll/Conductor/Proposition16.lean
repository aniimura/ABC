/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.ProductFormula
import ABC3.Found.GenEll.Conductor.Definition15

/-!
# Conductor —— `[GenEll] Proposition 1.6` の分

☆もとの 1 枚を**条なしの項目ごと**に割ったものである（第 1458、案 a）。
★「1 ファイル = 1 ノード」を回復するための分割で、中身は動かしていない。
-/

namespace ABC3.Found.GenEll

open NumberField IsDedekindDomain
variable {F : Type*} [Field F] [NumberField F]

@[simp] theorem ADivRed_fin (a : ADiv F) (v : FinitePlace F) :
    (ADivRed a).fin v = if 0 < a.fin v then (1 : ℤ) else 0 := rfl

@[simp] theorem ADivRed_arc (a : ADiv F) : (ADivRed a).arc = 0 := rfl

/-- ★`(−)_red` は**冪等**である。原文の `E = E_red` が「被約化の不動点」であることの確認。 -/
theorem adivRed_idem (a : ADiv F) : ADivRed (ADivRed a) = ADivRed a := by
  refine Prod.ext ?_ rfl
  ext v
  simp only [ADivRed, ADiv.fin, Finsupp.mapRange_apply]
  split_ifs <;> omega

/-- ★被約化は**有効**である。 -/
theorem adivRed_isEffective (a : ADiv F) : (ADivRed a).IsEffective := by
  constructor
  · intro v
    rw [ADivRed_fin]
    split_ifs <;> omega
  · intro v
    simp

/-! ## ★★被約化は次数を増やさない —— `Proposition 1.6` の非アルキメデス側 -/

/-- ★各有限素点の剰余体は 2 元以上なので `log q_v > 0`。 -/
theorem log_residueCard_pos (v : FinitePlace F) : 0 < Real.log (residueCard v) := by
  have h1 : 1 < residueCard v := by
    simpa [residueCard] using NumberField.HeightOneSpectrum.one_lt_absNorm v
  exact Real.log_pos (by exact_mod_cast h1)

/-- ★★**`deg((a)_red) ≤ deg(a)`**(`a` が有効なとき)。

原文 (GenEll p.9):
> Proposition 1.6. (Conductor Bounded by the Height) Let D ⊆ X be an effective Cartier divisor,

★**これが `Proposition 1.6` の証明の非アルキメデス側**である
(「the contributions at the nonarchimedean primes, from the definition of log-cond_D
[i.e., involving "(−)_red"]」)。

★機構は単純: 有効なら各係数 `n ≥ 0` で、被約化は `n` を `min(n,1)` に潰す。
`log q_v > 0` なので各項が減り、アルキメデス側は被約化で `0` になるが元は非負である。 -/
theorem deg_adivRed_le (a : ADiv F) (ha : a.IsEffective) : deg (ADivRed a) ≤ deg a := by
  classical
  -- 有限側: 台を `a.fin.support` に揃えて項ごとに比べる
  have hsub : (ADivRed a).fin.support ⊆ a.fin.support := by
    intro v hv
    simp only [Finsupp.mem_support_iff, ADivRed_fin] at hv ⊢
    intro hcon
    rw [hcon] at hv
    simp at hv
  have hfin1 : (ADivRed a).fin.sum (fun v n => (n : ℝ) * Real.log (residueCard v))
      = ∑ v ∈ a.fin.support,
          (((ADivRed a).fin v : ℤ) : ℝ) * Real.log (residueCard v) :=
    Finsupp.sum_of_support_subset _ hsub _ (fun v _ => by simp)
  have hfin2 : a.fin.sum (fun v n => (n : ℝ) * Real.log (residueCard v))
      = ∑ v ∈ a.fin.support, ((a.fin v : ℤ) : ℝ) * Real.log (residueCard v) := rfl
  have hle : ∑ v ∈ a.fin.support,
        (((ADivRed a).fin v : ℤ) : ℝ) * Real.log (residueCard v)
      ≤ ∑ v ∈ a.fin.support, ((a.fin v : ℤ) : ℝ) * Real.log (residueCard v) := by
    refine Finset.sum_le_sum fun v _ => ?_
    have hn : 0 ≤ a.fin v := ha.1 v
    have hcoef : (((ADivRed a).fin v : ℤ) : ℝ) ≤ ((a.fin v : ℤ) : ℝ) := by
      rw [ADivRed_fin]
      split_ifs with h
      · exact_mod_cast h
      · exact_mod_cast hn
    exact mul_le_mul_of_nonneg_right hcoef (log_residueCard_pos v).le
  -- アルキメデス側: 被約化は `0`、元は非負
  have harc : (0 : ℝ) ≤ a.arc.sum (fun _ r => r) :=
    Finset.sum_nonneg fun v _ => ha.2 v
  have harc0 : (ADivRed a).arc.sum (fun _ r => r) = 0 := by
    rw [ADivRed_arc]; simp
  rw [deg, deg, harc0, hfin1, hfin2]
  linarith

/-- ★★**条なしにしてはならない。**

`deg_adivRed_le` は `Proposition 1.6` の**非アルキメデス側だけ**である。
アルキメデス側(「コンパクト空間 `X^arc` 上の連続関数 `|s|_L` が有界」)は
複素解析空間を要求し、mathlib に 0 件である。

★★**2026-08-17 に自分で条なしを付けてしまい、`genell-progress` が 4/24 → 5/24 に
誤って動いた。** 直後の自己監査で発覚し、ここを条つきに直した。
★これは `.src` の 2 値規則が**まさに防ぐために存在する**誤りである。 -/
def deg_adivRed_le.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 9, item := "Proposition 1.6(証明の非アルキメデス側のみ)",
    sectionId := "genell-prop-1-6" }

end ABC3.Found.GenEll
