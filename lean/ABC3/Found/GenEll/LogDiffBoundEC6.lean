/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.VerticalBound
import ABC3.Found.GenEll.LogDiffTower
import ABC3.Found.GenEll.CyclotomicTame
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★段 EC6 の核 —— elementary claim は `log-diff` を一様に抑える（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.10。

原文 (GenEll p.10):
> integer n such that for any finite Galois extension L/K of finite extensions

## ★★★★★★★★★★これは何か —— 「有界である」の中身

`Proposition 1.7` の証明で、原文は Σ の上の食い違いについて

> to complete the proof of assertion (i), it suffices to show that the "portion over Σ"
> of the quantity log-diff_Y − log-diff_Z is bounded in U_Y(Q̄)

と書き、そのために elementary claim（`𝔡_{L/K} ⊇ p^n·O_L`）を出す。

★**なぜ elementary claim で「有界」が出るのか**を型で書くと:

    `p^n ∈ 𝔡` ⟹ `log-diff_K − log-diff_F ≤ n·log p`

★★右辺は **`K` にも `F` にも依らない**——これが「一様に有界」の中身である。

## ★★★機構 —— 在庫 2 本の合成

| 道具 | 役割 |
|---|---|
| `logDiffOfField_tower`（`§9-LogDiffTower`） | `log-diff_K − log-diff_F = log N(𝔡)/[K:ℚ]`（**等式**） |
| `degNormalized_idealADiv_le_log_of_natCast_mem`（`§9-VerticalBound`） | `(m : 𝓞_K) ∈ I` ⟹ `deg_normalized(I) ≤ log m` |

★★★★**測定（2026-08-28）**: 後者は既に手元にあった
——`§9-VerticalBound` は「2 つの `ℤ`-モデルの高さの差」のために作ったものだが、
**そのまま elementary claim の消費側になる**。
★CLAUDE.md の在庫の規律（木を grep する前にまず引く）がまた当たった。

## ★残っている段（明示）

★★elementary claim 本体（`p^n ∈ 𝔡` の `n` を `p` と `d` だけから作る）は
`EC1`–`EC5`（`§9-893`〜`§9-896`）の部品を塔で組めば出る。
★本ファイルはその**消費側**を先に用意したものである。
-/

namespace ABC3.Found.GenEll

open NumberField

/-! ## ★★★★★★★★★★`log-diff` の差は `log m` で抑えられる -/

/-- ★★★★★★★★★**`m ∈ 𝔡` なら `log-diff` の差は `log m` 以下**。

原文 (GenEll p.10):
> integer n such that for any finite Galois extension L/K of finite extensions

★`logDiffOfField_tower`（**等式**）に `§9-VerticalBound` の評価を差し込むだけである。
★★右辺は `K` にも `F` にも依らない。 -/
theorem logDiffOfField_sub_le_log (F K : Type) [Field F] [NumberField F] [Field K]
    [NumberField K] [Algebra F K] (m : ℕ) (hm : m ≠ 0)
    (hne : differentIdeal (𝓞 F) (𝓞 K) ≠ 0)
    (hmem : ((m : ℕ) : 𝓞 K) ∈ differentIdeal (𝓞 F) (𝓞 K)) :
    logDiffOfField K - logDiffOfField F ≤ Real.log m := by
  have hb := degNormalized_idealADiv_le_log_of_natCast_mem K
    (differentIdeal (𝓞 F) (𝓞 K)) m hm hmem
  rw [degNormalized, deg_idealADiv K _ hne] at hb
  rw [logDiffOfField_tower F K]
  linarith

/-- ★★★★★★★★★★**elementary claim の形** —— `p^n ∈ 𝔡` なら差は `n·log p` 以下。

原文 (GenEll p.10):
> integer n such that for any finite Galois extension L/K of finite extensions

★これが原文の「the portion over Σ … is bounded」の中身である
——`n` が `p` と `d` だけで決まる（`§9-894` の `exponent_le_log`）ので、
**点にも定義体にも依らない一様な上界**になる。 -/
theorem logDiffOfField_sub_le_of_pow_mem (F K : Type) [Field F] [NumberField F] [Field K]
    [NumberField K] [Algebra F K] (p n : ℕ) (hp : p ≠ 0)
    (hne : differentIdeal (𝓞 F) (𝓞 K) ≠ 0)
    (hmem : ((p ^ n : ℕ) : 𝓞 K) ∈ differentIdeal (𝓞 F) (𝓞 K)) :
    logDiffOfField K - logDiffOfField F ≤ (n : ℝ) * Real.log p := by
  have h := logDiffOfField_sub_le_log F K (p ^ n) (pow_ne_zero n hp) hne hmem
  rwa [Nat.cast_pow, Real.log_pow] at h

/-! ## ★出典の紐付け(`.src`) -/

def logDiffOfField_sub_le_log.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 10,
    item := "Proposition 1.7(m ∈ 𝔡 なら log-diff の差は log m 以下)",
    sectionId := "genell-prop-1-7" }

def logDiffOfField_sub_le_of_pow_mem.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 10,
    item := "Proposition 1.7(elementary claim の形——p^n ∈ 𝔡 なら差は n·log p 以下)",
    sectionId := "genell-prop-1-7" }

def logDiffOfField_sub_le_of_pow_mem.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "logDiffOfField_tower(log-diff の差は log N(𝔡)/[K:ℚ]、**等式**)"
      (.inProject "ABC3" "ABC3.Found.GenEll.logDiffOfField_tower") 1,
    .citation "[ABC3]" "degNormalized_idealADiv_le_log_of_natCast_mem(§9-VerticalBound)"
      (.inProject "ABC3" "ABC3.Found.GenEll.degNormalized_idealADiv_le_log_of_natCast_mem") 1,
    .implicitStep
      ("★★★★測定(2026-08-28): 後者は既に手元にあった" ++
       "——§9-VerticalBound は「2 つの ℤ-モデルの高さの差」のために作ったものだが、" ++
       "**そのまま elementary claim の消費側になる**。" ++
       "★CLAUDE.md の在庫の規律(木を grep する前にまず引く)がまた当たった") 2,
    .implicitStep
      ("★★elementary claim 本体(p^n ∈ 𝔡 の n を p と d だけから作る)は " ++
       "EC1–EC5(§9-893〜§9-896)の部品を塔で組めば出る。" ++
       "★本ファイルはその**消費側**を先に用意したものである") 3 ]

end ABC3.Found.GenEll
