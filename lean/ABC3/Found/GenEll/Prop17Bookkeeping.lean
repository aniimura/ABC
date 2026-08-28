/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.Prop17Fundamental
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★★★★★★★★★★★★★★正規化次数の帳簿（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.9。

原文 (GenEll p.9):
> Proposition 1.7. (Conductors and Log Differents) Let

## ★★★★★★★★★★★★★★★★★★★★★★★★★★★★これは何か —— 左の `≲` が等式である理由

`§9-970` で測ったとおり、`Proposition 1.7, (i)` の**左の `≲` は各素点で等式**である。
その等式の中身は**純粋に算術**であり、本ファイルがそれを取る:

    `(∑_i (e_i − 1)·f_i·L) / ([K:F]·[F:ℚ])  =  (1 − (∑_i f_i)/[K:F]) · (L/[F:ℚ])`

（`L ≔ log N(p)`、`∑_i e_i f_i = [K:F]` は**基本等式**）。

★左辺が `log-diff(K) − log-diff(F)` の `p` での寄与、
右辺が `log-cond_E(z) − log-cond_D(y)` の `p` での寄与である。

## ★★★機構 —— 基本等式で分子が畳める

    `∑_i (e_i − 1) f_i = ∑_i e_i f_i − ∑_i f_i = [K:F] − ∑_i f_i`

★あとは `[K:ℚ] = [K:F]·[F:ℚ]` で割るだけである。

## ★これで `Proposition 1.7` に残るのは

★★★**幾何の側の紐**だけである:

* `log N(P) = f_P · log N(p)`（剰余次数の定義——`SigmaBound.lean` に近い形がある）
* `∑_{P∣p} e_P f_P = [K:F]`（mathlib の `Ideal.sum_ramification_inertia`）
* 「`D_ℚ = φ_ℚ^{-1}(E_ℚ)_red`」（原文の条件 (b)）を導手の台の対応に翻訳すること

★数値と代数の側は本日（`§9-954`〜`§9-971`）すべて閉じた。
-/

namespace ABC3.Found.GenEll

/-! ## ★★★★★★★★★★★★★★★★★★★★★★★★★★★★等式の中身 -/

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★**`log-diff` の寄与 ＝
`log-cond_E − log-cond_D` の寄与**（基本等式のもとで）。

原文 (GenEll p.9):
> Proposition 1.7. (Conductors and Log Differents) Let

    `(∑_i (e_i − 1)·f_i·L) / (n·m)  =  (1 − (∑_i f_i)/n) · (L/m)`

（`n = [K:F]`、`m = [F:ℚ]`、`L = log N(p)`、`∑_i e_i f_i = n`）。

★★これが原文の『prime-to-`Σ` 部分では **`=` と `≤` であって `≲` ではない**』の
**算術の中身**である。 -/
theorem diff_contrib_eq (ι : Type*) [Fintype ι] (ee ff : ι → ℕ) (n m : ℕ)
    (he1 : ∀ i, 1 ≤ ee i)
    (hsum : ∑ i, ee i * ff i = n) (hn : 0 < n) (hm : 0 < m) (L : ℝ) :
    (∑ i, ((ee i - 1 : ℕ) : ℝ) * (ff i : ℝ) * L) / ((n : ℝ) * (m : ℝ))
      = (1 - (∑ i, (ff i : ℝ)) / (n : ℝ)) * (L / (m : ℝ)) := by
  have hnR : (0:ℝ) < (n : ℝ) := by exact_mod_cast hn
  have hmR : (0:ℝ) < (m : ℝ) := by exact_mod_cast hm
  have hcast : ∀ i, ((ee i - 1 : ℕ) : ℝ) = (ee i : ℝ) - 1 := by
    intro i
    have h := he1 i
    push_cast [Nat.cast_sub h]
    ring
  have hnum : (∑ i, ((ee i - 1 : ℕ) : ℝ) * (ff i : ℝ) * L)
      = ((n : ℝ) - ∑ i, (ff i : ℝ)) * L := by
    have hs : (∑ i, (ee i : ℝ) * (ff i : ℝ)) = (n : ℝ) := by
      have h := congrArg (fun k : ℕ => (k : ℝ)) hsum
      push_cast at h
      exact h
    calc (∑ i, ((ee i - 1 : ℕ) : ℝ) * (ff i : ℝ) * L)
        = ∑ i, (((ee i : ℝ) - 1) * (ff i : ℝ)) * L := by
          exact Finset.sum_congr rfl (fun i _ => by rw [hcast i])
      _ = (∑ i, ((ee i : ℝ) * (ff i : ℝ) - (ff i : ℝ))) * L := by
          rw [Finset.sum_mul]
          exact Finset.sum_congr rfl (fun i _ => by ring)
      _ = ((∑ i, (ee i : ℝ) * (ff i : ℝ)) - ∑ i, (ff i : ℝ)) * L := by
          rw [Finset.sum_sub_distrib]
      _ = ((n : ℝ) - ∑ i, (ff i : ℝ)) * L := by rw [hs]
  rw [hnum]
  field_simp

/-! ## ★出典の紐付け(`.src`) -/

def diff_contrib_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 9,
    item := "Proposition 1.7(log-diff の寄与 ＝ log-cond_E − log-cond_D の寄与)",
    sectionId := "genell-prop-1-7" }

def diff_contrib_eq.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]" "Ideal.sum_ramification_inertia(基本等式 ∑ e_P f_P = [K:F])"
      (.inMathlib "Ideal.sum_ramification_inertia") 3,
    .implicitStep
      ("★★★★★★測定(2026-08-29): Proposition 1.7, (i) の左の ≲ が等式である理由は" ++
       "**純粋に算術**である——基本等式で分子が畳める:" ++
       "∑ (e_i − 1) f_i = ∑ e_i f_i − ∑ f_i = [K:F] − ∑ f_i。" ++
       "あとは [K:ℚ] = [K:F]·[F:ℚ] で割るだけ。" ++
       "★これが原文の『prime-to-Σ 部分では = と ≤ であって ≲ ではない』の中身である") 6,
    .implicitStep
      ("★★これで Proposition 1.7 に残るのは**幾何の側の紐**だけである: " ++
       "log N(P) = f_P·log N(p)、∑_{P|p} e_P f_P = [K:F](mathlib)、" ++
       "そして原文の条件 (b)(D_ℚ = φ_ℚ^{-1}(E_ℚ)_red)を導手の台の対応に翻訳すること。" ++
       "★数値と代数の側は本日(§9-954〜§9-971)すべて閉じた") 5 ]

end ABC3.Found.GenEll
