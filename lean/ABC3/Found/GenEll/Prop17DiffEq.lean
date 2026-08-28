/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.Prop17Norms
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★`log-diff` の差は**等式**である（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.9。

原文 (GenEll p.9):
> Proposition 1.7. (Conductors and Log Differents) Let

## ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★これは何か —— 両側から挟んで等号

本日の 2 本を合わせると**等式**になる:

| 向き | 出どころ |
|---|---|
| `≥` | `§9-958`（`pow_sub_one_dvd_differentIdeal`——`P^{e−1} ∣ 𝔡`） |
| `≤` | `§9-968`（完全分岐・馴なら `𝔭^e ∤ 𝔡`） |

★★★★したがって

    **`log-diff(K) − log-diff(F) = (∑_{P∈S} (e_P − 1)·log N(P)) / [K:ℚ]`**

——`§9-971` の帳簿が要求していた**左辺そのもの**である。

## ★これで `Proposition 1.7, (i)` の左の `≲` は

| 部品 | 状態 |
|---|---|
| `log-diff` の差 ＝ `∑ (e_P−1) log N(P) / [K:ℚ]` | ★**本ファイル** |
| `∑ (e_i−1) f_i L /(n·m) = (1 − (∑ f_i)/n)·(L/m)` | ★`§9-971` |
| `log N(P) = f_P·log N(p)` | ★`§9-972` |
| `∑_{P∣p} e_P f_P = [K:F]` | ★mathlib |
| `log-cond` ＝ 台の上の `∑ log N` | ★`§9-966` |
| 条件 (b) を導手の台の対応に翻訳 | ☆**残る 1 本** |

★★★★★★数値・代数の側は**全部揃った**。残るのは原文の条件 (b)
（`D_ℚ = φ_ℚ^{-1}(E_ℚ)_red`）を台の対応に読み替える段だけである。
-/

namespace ABC3.Found.GenEll

open NumberField

/-! ## ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★等式 -/

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★**`log-diff` の差は等式である**
（分岐がすべて完全分岐・馴のとき）。

原文 (GenEll p.9):
> Proposition 1.7. (Conductors and Log Differents) Let

    `log-diff(K) − log-diff(F) = (∑_{P∈S} (e_P − 1)·log N(P)) / [K:ℚ]`

★`≥` は `§9-958`（mathlib の `pow_sub_one_dvd_differentIdeal`）、
`≤` は `§9-968`（完全分岐・馴なら `𝔭^e ∤ 𝔡`）。
★★これが `§9-971` の帳簿が要求していた左辺そのものである。 -/
theorem logDiff_tower_eq_sum_of_totallyRamified_tame (F K : Type) [Field F] [NumberField F]
    [Field K] [NumberField K] [Algebra F K] [DecidableEq (Ideal (𝓞 K))]
    (S : Finset (Ideal (𝓞 K))) (hmax : ∀ P ∈ S, P.IsMaximal) (hne : ∀ P ∈ S, P ≠ ⊥)
    (pr : Ideal (𝓞 K) → Ideal (𝓞 F)) (ram : Ideal (𝓞 K) → ℕ)
    (hpmax : ∀ P ∈ S, (pr P).IsMaximal) (hpne : ∀ P ∈ S, pr P ≠ ⊥)
    (hram1 : ∀ P ∈ S, 1 ≤ ram P)
    (hfac : ∀ P ∈ S, P ^ (ram P) = Ideal.map (algebraMap (𝓞 F) (𝓞 K)) (pr P))
    (htame : ∀ P ∈ S, (Algebra.intTrace (𝓞 F) (𝓞 K)) 1 ∉ pr P)
    (hI : differentIdeal (𝓞 F) (𝓞 K) ≠ 0)
    (hprodne : (∏ P ∈ S, P ^ (ram P - 1)) ≠ 0)
    (hramsupp : ∀ Q : Ideal (𝓞 K), ∀ _ : Q.IsPrime,
      ¬ Algebra.IsUnramifiedAt (𝓞 F) Q → Q ∈ S) :
    logDiffOfField K - logDiffOfField F
      = (∑ P ∈ S, ((ram P - 1 : ℕ) : ℝ) * Real.log (Ideal.absNorm P))
        / (Module.finrank ℚ K : ℝ) := by
  refine le_antisymm ?_ ?_
  · exact logDiff_tower_le_sum_of_totallyRamified_tame F K S hne pr ram hram1 hfac htame hI
      hprodne hramsupp
  · refine sum_log_le_logDiff_of_ramificationIdx F K S hmax hne pr ram hpmax hpne ?_ hI
    intro P hP
    rw [hfac P hP]

/-! ## ★出典の紐付け(`.src`) -/

def logDiff_tower_eq_sum_of_totallyRamified_tame.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 9,
    item := "Proposition 1.7(log-diff の差は等式である——完全分岐・馴のとき)",
    sectionId := "genell-prop-1-7" }

def logDiff_tower_eq_sum_of_totallyRamified_tame.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "sum_log_le_logDiff_of_ramificationIdx(≥ の側、§9-958)"
      (.inProject "ABC3" "ABC3.Found.GenEll.sum_log_le_logDiff_of_ramificationIdx") 3,
    .citation "[ABC3]" "logDiff_tower_le_sum_of_totallyRamified_tame(≤ の側、§9-968)"
      (.inProject "ABC3" "ABC3.Found.GenEll.logDiff_tower_le_sum_of_totallyRamified_tame") 3,
    .implicitStep
      ("★★★★★★測定(2026-08-29): 本日の 2 本(§9-958 の ≥ と §9-968 の ≤)を合わせると" ++
       "**等式**になる——log-diff(K) − log-diff(F) = (∑ (e_P−1)·log N(P))/[K:ℚ]。" ++
       "★これが §9-971 の帳簿が要求していた左辺そのものである") 6,
    .implicitStep
      ("★★★Proposition 1.7, (i) の左の ≲ に要る部品は" ++
       "本ファイル・§9-971(算術)・§9-972(ノルム)・mathlib(基本等式)・§9-966(導手)で" ++
       "**全部揃った**。残るのは原文の条件 (b)(D_ℚ = φ_ℚ^{-1}(E_ℚ)_red)を" ++
       "導手の台の対応に読み替える段**だけ**である") 5 ]

end ABC3.Found.GenEll
