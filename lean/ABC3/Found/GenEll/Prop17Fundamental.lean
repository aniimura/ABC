/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.Prop17Sandwich
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★★★★★★★★★★★★★`Proposition 1.7` は基本等式で読める（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.9。

原文 (GenEll p.9):
> Proposition 1.7. (Conductors and Log Differents) Let

## ★★★★★★★★★★★★★★★★★★★★★★★★★★★これは何か —— 原文の読み直し

`ResearchPaper/1_Structured` の原文（p.9）を読み直したところ、
`Proposition 1.7, (i)` の設定は

* `φ : Y → Z` は生成的有限、`D ⊆ Y`・`E ⊆ Z` は有効 `ℤ`-平坦 Cartier 因子
* (b) `D_ℚ = φ_ℚ^{-1}(E_ℚ)_red`
* (c) `φ_ℚ` は `(U_Y)_ℚ → (U_Z)_ℚ` 上**有限エタール**
* (d) `D_ℚ` の各点での分岐指数は `e` を**割る**

であり、結論は

    `log-cond_E − log-cond_D  ≲  log-diff_Y − log-diff_Z  ≲  (1 − 1/e)·log-cond_E`

である。★★★★★★**`log-cond_D` は分母を大きくする項ではなく、
基本等式で `log-diff` の差と釣り合う項である**。

## ★★★★★測定 —— 左の `≲` は各素点で**等式**である

点 `y ∈ U_Y(ℚ̄)` の定義体を `K`、その像 `z` の定義体を `F` とし、
`p` を `F` の素点、`P ∣ p` を `K` の素点、`e_P`・`f_P` をその分岐指数・剰余次数とすると:

| 量 | `p` での寄与（正規化後） |
|---|---|
| `log-cond_E(z)` | `log N(p) / [F:ℚ]` |
| `log-cond_D(y)` | `(∑_{P∣p} f_P / [K:F]) · log N(p) / [F:ℚ]` |
| `log-diff(K) − log-diff(F)` | `(∑_{P∣p} (e_P − 1) f_P / [K:F]) · log N(p) / [F:ℚ]` |

★**基本等式** `∑_{P∣p} e_P f_P = [K:F]`（mathlib の `Ideal.sum_ramification_inertia`）より

    `∑ (e_P − 1) f_P / [K:F] = 1 − (∑ f_P)/[K:F]`

——すなわち **`log-diff` の差 ＝ `log-cond_E − log-cond_D`（等式）**。
★★★これが原文が『prime-to-`Σ` 部分では **`=` と `≤` であって `≲` ではない**』と
注記していることの中身である。

## ★★★右の `≲` は `e_P ≤ e` だけ

`e_P ∣ e`（条件 (d)）より `e_P ≤ e` なので `∑ e_P f_P ≤ e·∑ f_P`、
したがって `(∑ f_P)/[K:F] ≥ 1/e`、ゆえに

    `1 − (∑ f_P)/[K:F]  ≤  1 − 1/e`

★★**本ファイルはこの数値の核を証明する**（`one_sub_ratio_le`）。

## ★残っている段（明示）

★★★残るのは**正規化次数の帳簿**だけである:

* `log N(P) = f_P · log N(p)`（`SigmaBound.lean` の `log_residueCard_eq_inertiaDeg_mul` が近い）
* `[K:ℚ] = [K:F]·[F:ℚ]`
* `∑_{P∣p} e_P f_P = [K:F]`（mathlib の `Ideal.sum_ramification_inertia`）

★これらを繋げば `Proposition 1.7, (i)` は**両側とも**出る。
`§9-969` の挟み込み（`(1−1/e)` と `(e−1)`）は `log-cond_D` を使わない**粗い版**であり、
本ファイルの読みが**原文どおりの版**への道である。
-/

namespace ABC3.Found.GenEll

/-! ## ★★★★数値の核 -/

/-- ★★★★**`e_P ≤ e` なら `∑ e_P f_P ≤ e·∑ f_P`**——基本等式の右辺を押さえる形。 -/
theorem sum_ram_le_mul (ι : Type*) [Fintype ι] (ee ff : ι → ℕ) (e : ℕ)
    (hle : ∀ i, ee i ≤ e) :
    (∑ i, (ee i : ℝ) * (ff i : ℝ)) ≤ (e : ℝ) * ∑ i, (ff i : ℝ) := by
  rw [Finset.mul_sum]
  refine Finset.sum_le_sum (fun i _ => ?_)
  have h1 : ((ee i : ℝ)) ≤ (e : ℝ) := by exact_mod_cast hle i
  have h2 : (0:ℝ) ≤ (ff i : ℝ) := by positivity
  exact mul_le_mul_of_nonneg_right h1 h2

/-- ★★★★★★★★★★★★★★**`Proposition 1.7` の右の `≲` の数値の核**。

原文 (GenEll p.9):
> Proposition 1.7. (Conductors and Log Differents) Let

    `1 − (∑ f_P)/(∑ e_P f_P)  ≤  1 − 1/e`   （`e_P ≤ e` のとき）

★左辺は**基本等式** `∑ e_P f_P = [K:F]` のもとで
`log-diff` の差／`log-cond_E` の比であり、右辺が原文の係数である。
★★したがって原文の右の `≲` は **`e_P ≤ e` だけ**から出る。 -/
theorem one_sub_ratio_le (ι : Type*) [Fintype ι] (ee ff : ι → ℕ) (e : ℕ) (he : 0 < e)
    (hle : ∀ i, ee i ≤ e) (heepos : ∀ i, 0 < ee i)
    (hfpos : (0:ℝ) < ∑ i, (ff i : ℝ)) :
    1 - (∑ i, (ff i : ℝ)) / (∑ i, (ee i : ℝ) * (ff i : ℝ)) ≤ 1 - 1 / (e : ℝ) := by
  have hepos : (0:ℝ) < (e : ℝ) := by exact_mod_cast he
  have hsum := sum_ram_le_mul ι ee ff e hle
  have hden : (0:ℝ) < ∑ i, (ee i : ℝ) * (ff i : ℝ) := by
    refine lt_of_lt_of_le hfpos (Finset.sum_le_sum (fun i _ => ?_))
    have h1 : (1:ℝ) ≤ (ee i : ℝ) := by exact_mod_cast heepos i
    have h2 : (0:ℝ) ≤ (ff i : ℝ) := by positivity
    nlinarith
  have hkey : 1 / (e : ℝ) ≤ (∑ i, (ff i : ℝ)) / (∑ i, (ee i : ℝ) * (ff i : ℝ)) :=
    (div_le_div_iff₀ hepos hden).mpr (by linarith)
  linarith

/-! ## ★出典の紐付け(`.src`) -/

def sum_ram_le_mul.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 9,
    item := "Proposition 1.7(e_P ≤ e なら ∑ e_P f_P ≤ e·∑ f_P)",
    sectionId := "genell-prop-1-7" }

def one_sub_ratio_le.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 9,
    item := "Proposition 1.7(右の ≲ の数値の核——基本等式のもとで)",
    sectionId := "genell-prop-1-7" }

def one_sub_ratio_le.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]" "Ideal.sum_ramification_inertia(基本等式 ∑ e_P f_P = [K:F])"
      (.inMathlib "Ideal.sum_ramification_inertia") 3,
    .implicitStep
      ("★★★★★★測定(2026-08-29、原文 p.9 を 1_Structured で読み直し): " ++
       "Proposition 1.7, (i) の**左の ≲ は各素点で等式である**。" ++
       "log-cond_E(z) の p での寄与が log N(p)/[F:ℚ]、" ++
       "log-cond_D(y) が (∑ f_P/[K:F])·log N(p)/[F:ℚ]、" ++
       "log-diff の差が (∑ (e_P−1) f_P/[K:F])·log N(p)/[F:ℚ] であり、" ++
       "基本等式 ∑ e_P f_P = [K:F] から**両者が一致する**。" ++
       "★これが原文の『prime-to-Σ 部分では = と ≤ であって ≲ ではない』の中身である") 7,
    .implicitStep
      ("★★右の ≲ は e_P ∣ e(条件 (d))から e_P ≤ e、" ++
       "したがって (∑ f_P)/[K:F] ≥ 1/e——**本ファイルが証明する数値の核**である") 4,
    .implicitStep
      ("★残るのは正規化次数の帳簿だけである: " ++
       "log N(P) = f_P·log N(p)、[K:ℚ] = [K:F]·[F:ℚ]、∑ e_P f_P = [K:F]。" ++
       "★★§9-969 の挟み込み((1−1/e) と (e−1))は log-cond_D を使わない**粗い版**であり、" ++
       "本ファイルの読みが原文どおりの版への道である") 5 ]

end ABC3.Found.GenEll
