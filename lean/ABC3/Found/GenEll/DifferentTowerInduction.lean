/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.PGroupChain
import ABC3.Found.GenEll.DifferentKummer
import ABC3.Found.GenEll.LogDiffBoundEC6
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★塔を組む帰納 —— elementary claim の `p`-群の側（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.10。

原文 (GenEll p.10):
> integer n such that for any finite Galois extension L/K of finite extensions

## ★★★★★★★★★★★★★これは何か —— 部品を塔で組む

`§9-892` で割った elementary claim の段のうち、部品は出揃っていた:

| 部品 | 場所 |
|---|---|
| 次数 `p` の段（`p^p ∈ 𝔡`） | `§9-DifferentKummer`（`EC2`） |
| 塔の合成 | `§9-DifferentKummer`（`EC3`） |
| Galois な中間体 | `§9-900`（`EC4`） |
| 指数の一様化 | `§9-898`（`EC6`） |

★本ファイルは**それらを帰納で組む**:

    `[L:K] = p^k` かつ Galois ⟹ `p^{k·p} ∈ 𝔡_{L/K}`

★★さらに `[L:K] ≤ d` なら `k ≤ log_p d` なので、
**`N ≔ p·log_p d` で一様に**`p^N ∈ 𝔡_{L/K}` が成り立つ。

## ★★★機構 —— 中間体で 1 段ずつ

★`§9-900` が `[L:M] = p`（Galois）かつ `[M:K] = p^k`（Galois）な中間体を与える。
★★帰納法の仮定を `M/K` に、次数 `p` の段を `L/M` に当て、
different の乗法性（mathlib の `differentIdeal_eq_differentIdeal_mul_differentIdeal`）で継ぐ。

## ★★★★測定の記録

★★★★★**`base`（`[L:K] = 1` なら `𝔡 = ⊤`）は仮定として受けている**。
mathlib に `differentIdeal … = ⊤` の判定（不分岐との同値、自明拡大の場合）が
**見当たらなかった**（`Algebra.isUnramified_iff_differentIdeal_eq_top`・
`differentIdeal_eq_top_iff`・`differentIdeal_self` いずれも無し、2026-08-28 実測）。
★中身は「次数 1 の拡大は同型なので different は単位イデアル」であり、
数学としては自明だが、mathlib の語彙に無いので**明示的に受けた**。

## ★★★★★★後日の追記（2026-08-28、`§9-905`）—— `base` は外れた

★★★★★上の測定は**「その名前で無い」でしかなかった**。
`differentIdeal` を主語にした判定は確かに無いが、**判別式を主語にすれば全部あった**:

    `[L:K] = 1` ⟹ `K ≃ₐ[ℚ] L` ⟹ `|disc L| = |disc K|`
                ⟹ `N(𝔡) = 1`（`|disc L| = N(𝔡)·|disc K|^{[L:K]}` から）⟹ `𝔡 = ⊤`

★`base` を外した形は `Found/GenEll/DifferentTrivialDegree.lean`
（`pow_mem_differentIdeal_of_pow_finrank'`・`pow_mem_uniform_of_finrank_le'`）。
★★本ファイルの形（`base` を仮定に持つ）は**一般形として残す**——
局所体版など `NumberField` 以外に当てるときに要る。
-/

namespace ABC3.Found.GenEll

open NumberField IntermediateField

/-! ## ★★★★★★★★★★★★★塔の帰納 -/

/-- ★★★★★★★★★★★★★**次数 `p^k` の Galois 拡大では `p^{k·p} ∈ 𝔡`**。

原文 (GenEll p.10):
> integer n such that for any finite Galois extension L/K of finite extensions

★`§9-900` の Galois な中間体で 1 段ずつ降り、
`§9-DifferentKummer` の塔公式で継ぐ帰納である。
★★`base` と `step` は `EC1`/`EC2` の内容であり、仮定として受けている。 -/
theorem pow_mem_differentIdeal_of_pow_finrank (p : ℕ) [Fact p.Prime]
    (base : ∀ (K L : Type) [Field K] [NumberField K] [Field L] [NumberField L]
      [Algebra K L] [IsGalois K L], Module.finrank K L = 1 →
      differentIdeal (𝓞 K) (𝓞 L) = ⊤)
    (step : ∀ (K L : Type) [Field K] [NumberField K] [Field L] [NumberField L]
      [Algebra K L] [IsGalois K L], Module.finrank K L = p →
      ((p : 𝓞 L)) ^ p ∈ differentIdeal (𝓞 K) (𝓞 L)) :
    ∀ (k : ℕ) (K L : Type) [Field K] [NumberField K] [Field L] [NumberField L]
      [Algebra K L] [IsGalois K L], Module.finrank K L = p ^ k →
      ((p : 𝓞 L)) ^ (k * p) ∈ differentIdeal (𝓞 K) (𝓞 L) := by
  intro k
  induction k with
  | zero =>
      intro K L _ _ _ _ _ _ h
      rw [base K L (by simpa using h)]
      exact Submodule.mem_top
  | succ k ih =>
      intro K L _ _ _ _ _ _ h
      obtain ⟨M, hML, hKM, hgal⟩ := exists_intermediateField_galois p K L h
      haveI := hgal
      have hAB := ih K M hKM
      have hBC := step M L hML
      have hmul := differentIdeal_eq_differentIdeal_mul_differentIdeal (𝓞 K) (𝓞 M) (𝓞 L)
      have h2 := pow_mem_differentIdeal_tower (𝓞 K) (𝓞 M) (𝓞 L) p (k * p) p hmul hAB hBC
      have heq : p + k * p = (k + 1) * p := by ring
      rwa [heq] at h2

/-! ## ★★★★★★★★★★★★★★一様な形 -/

/-- ★★★★★★★★★★★★★★**elementary claim の `p`-群の側**（一様な形）。

原文 (GenEll p.10):
> integer n such that for any finite Galois extension L/K of finite extensions

    `[L:K] = p^k ≤ d` ⟹ `p^N ∈ 𝔡_{L/K}`   （`N ≔ p·log_p d` 以上なら何でもよい）

★右辺の `N` は **`K` にも `L` にも依らない**——これが原文の
`Fix a prime number p and a positive integer d. Then there exists a positive integer n`
の「there exists」の中身である。 -/
theorem pow_mem_uniform_of_finrank_le (p d : ℕ) [Fact p.Prime]
    (base : ∀ (K L : Type) [Field K] [NumberField K] [Field L] [NumberField L]
      [Algebra K L] [IsGalois K L], Module.finrank K L = 1 →
      differentIdeal (𝓞 K) (𝓞 L) = ⊤)
    (step : ∀ (K L : Type) [Field K] [NumberField K] [Field L] [NumberField L]
      [Algebra K L] [IsGalois K L], Module.finrank K L = p →
      ((p : 𝓞 L)) ^ p ∈ differentIdeal (𝓞 K) (𝓞 L))
    (hd : d ≠ 0) (N : ℕ) (hN : Nat.log p d * p ≤ N)
    (K L : Type) [Field K] [NumberField K] [Field L] [NumberField L]
    [Algebra K L] [IsGalois K L] (k : ℕ) (hk : Module.finrank K L = p ^ k)
    (hle : Module.finrank K L ≤ d) :
    ((p : 𝓞 L)) ^ N ∈ differentIdeal (𝓞 K) (𝓞 L) := by
  have hp2 : 2 ≤ p := (Fact.out : p.Prime).two_le
  refine pow_mem_uniform _ p d k N hp2 hd (by rw [← hk]; exact hle) hN ?_
  exact pow_mem_differentIdeal_of_pow_finrank p base step k K L hk

/-! ## ★出典の紐付け(`.src`) -/

def pow_mem_differentIdeal_of_pow_finrank.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 10,
    item := "Proposition 1.7(elementary claim——次数 p^k の Galois 拡大では p^{k·p} ∈ 𝔡)",
    sectionId := "genell-prop-1-7" }

def pow_mem_uniform_of_finrank_le.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 10,
    item := "Proposition 1.7(elementary claim の p-群の側——一様な形)",
    sectionId := "genell-prop-1-7" }

def pow_mem_uniform_of_finrank_le.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "exists_intermediateField_galois(Galois な中間体、§9-900)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_intermediateField_galois") 2,
    .citation "[ABC3]" "pow_mem_differentIdeal_tower(塔で継ぐ、§9-DifferentKummer)"
      (.inProject "ABC3" "ABC3.Found.GenEll.pow_mem_differentIdeal_tower") 2,
    .citation "[mathlib]" "differentIdeal_eq_differentIdeal_mul_differentIdeal(different の乗法性)"
      (.inMathlib "differentIdeal_eq_differentIdeal_mul_differentIdeal") 1,
    .implicitStep
      ("★★★★★測定: **base([L:K] = 1 なら 𝔡 = ⊤)は仮定として受けている**。" ++
       "mathlib に differentIdeal … = ⊤ の判定が**見当たらなかった**" ++
       "(Algebra.isUnramified_iff_differentIdeal_eq_top・differentIdeal_eq_top_iff・" ++
       "differentIdeal_self いずれも無し、2026-08-28 実測)。" ++
       "★中身は「次数 1 の拡大は同型なので different は単位イデアル」で数学としては自明だが、" ++
       "mathlib の語彙に無いので明示的に受けた") 3,
    .implicitStep
      ("★★step(次数 p の Galois 拡大で p^p ∈ 𝔡)は EC1/EC2 の内容であり、" ++
       "§9-DifferentKummer の pow_mem_differentIdeal_of_kummer と " ++
       "mem_differentIdeal_of_isUnit_natCast を場合分けで繋げば出る") 3,
    .implicitStep
      ("★★★これで elementary claim の p-群の側は組み上がった。" ++
       "残るのは ζ_p の添加(§9-895・§9-896 が馴であることを言っている)を" ++
       "塔の一段として挟む配管である") 3 ]

end ABC3.Found.GenEll
