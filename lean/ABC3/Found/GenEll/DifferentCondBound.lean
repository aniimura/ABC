/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.DifferentRamIdx
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★★★★★★★★`(1 − 1/e)·log-cond ≤ log-diff` の差（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.9。

原文 (GenEll p.9):
> Proposition 1.7. (Conductors and Log Differents) Let

## ★★★★★★★★★★★★★★★★★★★★★★これは何か —— `hup` の原文どおりの形

`Skeleton/GenEll/Section1.lean` の `prop_1_7` が受けている

    `hup : (1 − 1/e)·cond_E x ≤ (diff_Y x − diff_Z x) + slackUp x`

の**主要項**を取る。★`§9-958` は

    `(∑_P (e_P − 1)·log N(P)) / [K:ℚ] ≤ log-diff(K) − log-diff(F)`

を与えていた。★★あとは**分岐している素点では `e_P − 1 ≥ 1 ≥ 1 − 1/e`** なので、
`(1 − 1/e)` の因子が `e_P − 1` に吸収される:

    **`(1 − 1/e)·(∑_P log N(P)) / [K:ℚ]  ≤  log-diff(K) − log-diff(F)`**

★★★左辺の `∑_P log N(P)` が**導手**（分岐する素点の被約な寄与）である。

## ★★★機構 —— 項ごとの不等式

`log N(P) ≥ 0`（`N(P) ≥ 1`）と `1 − 1/e ≤ 1 ≤ e_P − 1`（`e_P ≥ 2`）だけである。
★原文が `(1 − 1/e)` という係数を書くのは、**`e_P` が `e` で抑えられている**からであり、
不等式の向きから見ると `e_P ≥ 2` さえあれば `1` で足りる。
★★`(1 − 1/e)` はそれより弱い係数なので、なおさら成り立つ。

## ★残っている段（明示）

★★`log-cond_E` を**因子の引き戻しの被約化**（`logCond`、`CartierPullback.lean`）として
定義した形と、本ファイルの `∑_P log N(P)` を結ぶ段が残る
——「被覆 `Y → Z` が `D` の上でだけ分岐する」という幾何の仮定が要る。
★これが `Proposition 1.7` の**導手の側**であり、different の側は本日閉じた。
-/

namespace ABC3.Found.GenEll

open NumberField

/-! ## ★★★★係数 `(1 − 1/e)` は `e_P − 1` に吸収される -/

/-- ★★★★**分岐している素点では `(1 − 1/e)` の係数は `e_P − 1` に吸収される**。

★`log N(P) ≥ 0` と `1 − 1/e ≤ 1 ≤ a_P` だけの項ごとの不等式である。 -/
theorem one_sub_inv_mul_sum_le (K : Type) [Field K] [NumberField K]
    (e : ℕ) (he : 0 < e) (S : Finset (Ideal (𝓞 K))) (hne : ∀ P ∈ S, P ≠ ⊥)
    (a : Ideal (𝓞 K) → ℕ) (h1 : ∀ P ∈ S, 1 ≤ a P) :
    (1 - 1 / (e : ℝ)) * (∑ P ∈ S, Real.log (Ideal.absNorm P))
      ≤ ∑ P ∈ S, (a P : ℝ) * Real.log (Ideal.absNorm P) := by
  rw [Finset.mul_sum]
  refine Finset.sum_le_sum (fun P hP => ?_)
  have hNpos : 1 ≤ (Ideal.absNorm P : ℝ) := by
    have h : Ideal.absNorm P ≠ 0 := by rw [Ne, Ideal.absNorm_eq_zero_iff]; exact hne P hP
    exact_mod_cast Nat.one_le_iff_ne_zero.mpr h
  have hlog : 0 ≤ Real.log (Ideal.absNorm P) := Real.log_nonneg hNpos
  have hfe : (0:ℝ) < (e : ℝ) := by exact_mod_cast he
  have hle : (1 - 1 / (e : ℝ)) ≤ (a P : ℝ) := by
    have h1' : (1:ℝ) ≤ (a P : ℝ) := by exact_mod_cast h1 P hP
    have hinv : (0:ℝ) ≤ 1 / (e : ℝ) := by positivity
    linarith
  exact mul_le_mul_of_nonneg_right hle hlog

/-! ## ★★★★★★★★★★★★★★★★★★★★★★`hup` の原文どおりの形 -/

/-- ★★★★★★★★★★★★★★★★★★★★★★**`(1 − 1/e)·log-cond ≤ log-diff の差`**
—— `Proposition 1.7` の `hup`。

原文 (GenEll p.9):
> Proposition 1.7. (Conductors and Log Differents) Let

    分岐する素点 `P ∈ S`（`e_P ≥ 2`）について
      `(1 − 1/e)·(∑_P log N(P)) / [K:ℚ]  ≤  log-diff(K) − log-diff(F)`

★`§9-958`（分岐指数からの下からの評価）に係数の吸収を掛けただけである。
★★左辺の `∑_P log N(P)` が**導手**（分岐する素点の被約な寄与）である。 -/
theorem one_sub_inv_mul_cond_le_logDiff (F K : Type) [Field F] [NumberField F]
    [Field K] [NumberField K] [Algebra F K]
    (e : ℕ) (he : 0 < e)
    (S : Finset (Ideal (𝓞 K))) (hmax : ∀ P ∈ S, P.IsMaximal) (hne : ∀ P ∈ S, P ≠ ⊥)
    (p : Ideal (𝓞 K) → Ideal (𝓞 F)) (ram : Ideal (𝓞 K) → ℕ)
    (hpmax : ∀ P ∈ S, (p P).IsMaximal) (hpne : ∀ P ∈ S, p P ≠ ⊥)
    (hdvd : ∀ P ∈ S, P ^ (ram P) ∣ (p P).map (algebraMap (𝓞 F) (𝓞 K)))
    (hram2 : ∀ P ∈ S, 2 ≤ ram P)
    (hI : differentIdeal (𝓞 F) (𝓞 K) ≠ 0) :
    (1 - 1 / (e : ℝ)) * (∑ P ∈ S, Real.log (Ideal.absNorm P))
        / (Module.finrank ℚ K : ℝ)
      ≤ logDiffOfField K - logDiffOfField F := by
  have hfr : (0:ℝ) < (Module.finrank ℚ K : ℝ) := by exact_mod_cast Module.finrank_pos
  have hstep := one_sub_inv_mul_sum_le K e he S hne (fun P => ram P - 1)
    (fun P hP => by have := hram2 P hP; omega)
  refine le_trans ((div_le_div_iff_of_pos_right hfr).mpr hstep) ?_
  exact sum_log_le_logDiff_of_ramificationIdx F K S hmax hne p ram hpmax hpne hdvd hI

/-! ## ★出典の紐付け(`.src`) -/

def one_sub_inv_mul_sum_le.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 9,
    item := "Proposition 1.7(分岐している素点では (1 − 1/e) の係数は e_P − 1 に吸収される)",
    sectionId := "genell-prop-1-7" }

def one_sub_inv_mul_cond_le_logDiff.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 9,
    item := "Proposition 1.7((1 − 1/e)·log-cond ≤ log-diff の差——hup)",
    sectionId := "genell-prop-1-7" }

def one_sub_inv_mul_cond_le_logDiff.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "sum_log_le_logDiff_of_ramificationIdx(分岐指数からの評価、§9-958)"
      (.inProject "ABC3" "ABC3.Found.GenEll.sum_log_le_logDiff_of_ramificationIdx") 3,
    .implicitStep
      ("★★★★測定(2026-08-29): 原文の (1 − 1/e) という係数は" ++
       "**e_P が e で抑えられている**ことから来るが、不等式の向きから見ると " ++
       "e_P ≥ 2 さえあれば 1 で足りる。(1 − 1/e) はそれより弱い係数なので" ++
       "なおさら成り立つ——項ごとの不等式で吸収される") 4,
    .implicitStep
      ("★★残るのは log-cond_E を**因子の引き戻しの被約化**(logCond、CartierPullback.lean)" ++
       "として定義した形と、本ファイルの ∑_P log N(P) を結ぶ段である" ++
       "——「被覆 Y → Z が D の上でだけ分岐する」という幾何の仮定が要る。" ++
       "★これが Proposition 1.7 の**導手の側**であり、different の側は本日閉じた") 5 ]

end ABC3.Found.GenEll
