/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.DifferentDescent
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★降下は要らなかった —— 消費側は単調性で足りる（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.10。

原文 (GenEll p.10):
> integer n such that for any finite Galois extension L/K of finite extensions

## ★★★★★★★★★★★★これは何か —— `EC8`（イデアルの降下）の**回避**

`§9-902` で `ζ_p` の還元の前半

    `p^n ∈ 𝔡_{K, M}` ⟹ `p^n ∈ (𝔡_{K,L})·O_M`   （`M = L(ζ)`）

を取り、`§9-903` で残る段を **`EC8`（降下 `(I·C) ∩ B = I`）** として立てた。
★そして mathlib にその形が**無い**ことを測った。

★★★★**しかし消費側は降下を要求しない。**

`Proposition 1.7` が `𝔡_{K,L}` を使うのは **`log-diff` を通してだけ**である。
そして `log-diff` は **`§9-LogDiffTower` の `logDiffOfField_le` により単調**:

    `L ⊆ M` ⟹ `log-diff(L) ≤ log-diff(M)`

したがって

    `log-diff(L) − log-diff(K) ≤ log-diff(M) − log-diff(K) ≤ n·log p`

——**`M` での主張がそのまま `L` での主張を含む**。★イデアルを降ろす必要が無い。

## ★★★測定の記録

| 見立て（`§9-903`） | 実際 |
|---|---|
| `EC8` は自作が要る（忠実平坦性 or 素点ごとの比較） | ★★★★**そもそも要らない** |

★★★★★**「イデアルで降ろす」から「量で挟む」への読み替え**である。
原文が `L[ζ]/K` を経由するのは `Gal` を扱いやすくするためであり、
最後に `L/K` へ戻す段は——**消費される量が単調なら**——自動である。

★★これは CLAUDE.md の「逸脱」に当たる読み替えなので明示する:
本ファイルは `𝔡_{K,L} ∋ p^n` を**証明していない**。
証明しているのは**それを使って出るはずだった結論**
（`log-diff` の差が一様に有界）の方である。
★後続（`Proposition 1.7` の `hup`）が要求するのは後者だけなので、影響は無い。

## ★機構

| 道具 | 役割 |
|---|---|
| `logDiffOfField_sub_le_of_pow_mem`（`§9-897`） | `p^n ∈ 𝔡_{K,M}` ⟹ 差 `≤ n·log p` |
| `logDiffOfField_le`（`§9-LogDiffTower`） | `L ⊆ M` ⟹ `log-diff(L) ≤ log-diff(M)` |

★★★★**互換性（`IsScalarTower K L M`）は要らなかった**——
証明は `K → M` と `L → M` の 2 本しか使わない（2026-08-28 実測）。
仮定を弱いままにしてある。
-/

namespace ABC3.Found.GenEll

open NumberField

/-! ## ★★★★★★★★★★★★中間体へは単調性で降りる -/

/-- ★★★★★★★★★★★★**降下を経由しない中間体への移送**。

原文 (GenEll p.10):
> integer n such that for any finite Galois extension L/K of finite extensions

    `p^n ∈ 𝔡_{K,M}` かつ `K ⊆ M`、`L ⊆ M` ⟹ `log-diff(L) − log-diff(K) ≤ n·log p`

★`M` での主張が `L` での主張を**含む**——`log-diff` が単調だからである。
★★これが `EC8`（イデアルの降下）を**不要にする**読み替えである。 -/
theorem logDiffOfField_sub_le_of_pow_mem_of_le (K L M : Type) [Field K] [NumberField K]
    [Field L] [NumberField L] [Field M] [NumberField M]
    [Algebra K M] [Algebra L M] (p n : ℕ) (hp : p ≠ 0)
    (hmem : ((p ^ n : ℕ) : 𝓞 M) ∈ differentIdeal (𝓞 K) (𝓞 M)) :
    logDiffOfField L - logDiffOfField K ≤ (n : ℝ) * Real.log p := by
  have h1 : logDiffOfField M - logDiffOfField K ≤ (n : ℝ) * Real.log p :=
    logDiffOfField_sub_le_of_pow_mem K M p n hp differentIdeal_ne_bot hmem
  have h2 : logDiffOfField L ≤ logDiffOfField M := logDiffOfField_le L M
  linarith

/-! ## ★★★★★★★★★★★★★★組み上げ —— `p`-群の側と繋ぐ -/

/-- ★★★★★★★★★★★★★★**elementary claim の消費側（一様な形）**。

原文 (GenEll p.10):
> integer n such that for any finite Galois extension L/K of finite extensions

`§9-901` の `pow_mem_uniform_of_finrank_le`（`p`-群の側）を、
本ファイルの単調性で **`K ⊆ M` の任意の中間体 `L`** へ移す:

    `[M:K] = p^k ≤ d`（Galois）、`L ⊆ M` ⟹ `log-diff(L) − log-diff(K) ≤ N·log p`

★右辺の `N ≔ p·log_p d` は **`K`・`L`・`M` のいずれにも依らない**
——これが原文の `there exists a positive integer n` の中身である。 -/
theorem logDiffOfField_sub_le_uniform_of_le (p d : ℕ) [Fact p.Prime]
    (base : ∀ (K L : Type) [Field K] [NumberField K] [Field L] [NumberField L]
      [Algebra K L] [IsGalois K L], Module.finrank K L = 1 →
      differentIdeal (𝓞 K) (𝓞 L) = ⊤)
    (step : ∀ (K L : Type) [Field K] [NumberField K] [Field L] [NumberField L]
      [Algebra K L] [IsGalois K L], Module.finrank K L = p →
      ((p : 𝓞 L)) ^ p ∈ differentIdeal (𝓞 K) (𝓞 L))
    (hd : d ≠ 0) (N : ℕ) (hN : Nat.log p d * p ≤ N)
    (K L M : Type) [Field K] [NumberField K] [Field L] [NumberField L]
    [Field M] [NumberField M] [Algebra K M] [IsGalois K M] [Algebra L M]
    (k : ℕ) (hk : Module.finrank K M = p ^ k) (hle : Module.finrank K M ≤ d) :
    logDiffOfField L - logDiffOfField K ≤ (N : ℝ) * Real.log p := by
  have hmem := pow_mem_uniform_of_finrank_le p d base step hd N hN K M k hk hle
  have hp : p ≠ 0 := (Fact.out : p.Prime).ne_zero
  refine logDiffOfField_sub_le_of_pow_mem_of_le K L M p N hp ?_
  push_cast
  exact hmem

/-! ## ★出典の紐付け(`.src`) -/

def logDiffOfField_sub_le_of_pow_mem_of_le.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 10,
    item := "Proposition 1.7(elementary claim——降下を経由しない中間体への移送)",
    sectionId := "genell-prop-1-7" }

def logDiffOfField_sub_le_uniform_of_le.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 10,
    item := "Proposition 1.7(elementary claim の消費側——一様な形)",
    sectionId := "genell-prop-1-7" }

def logDiffOfField_sub_le_of_pow_mem_of_le.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "logDiffOfField_sub_le_of_pow_mem(p^n ∈ 𝔡 なら差は n·log p 以下、§9-897)"
      (.inProject "ABC3" "ABC3.Found.GenEll.logDiffOfField_sub_le_of_pow_mem") 1,
    .citation "[ABC3]" "logDiffOfField_le(log-diff は体を大きくすると増える、§9-LogDiffTower)"
      (.inProject "ABC3" "ABC3.Found.GenEll.logDiffOfField_le") 1,
    .implicitStep
      ("★★★★測定(2026-08-28): §9-903 で立てた EC8(イデアルの降下 (I·C) ∩ B = I)は " ++
       "**そもそも要らない**。Proposition 1.7 が 𝔡_{K,L} を使うのは log-diff を通してだけで、" ++
       "log-diff は単調(logDiffOfField_le)なので M での主張が L での主張を含む") 4,
    .implicitStep
      ("★★逸脱の明示: 本ファイルは 𝔡_{K,L} ∋ p^n を**証明していない**。" ++
       "証明しているのはそれを使って出るはずだった結論(log-diff の差が一様に有界)の方である。" ++
       "後続(Proposition 1.7 の hup)が要求するのは後者だけなので影響は無い") 3,
    .implicitStep
      ("★★★★互換性(IsScalarTower K L M)は要らなかった" ++
       "——証明は K → M と L → M の 2 本しか使わない(2026-08-28 実測)。" ++
       "仮定を弱いままにしてある") 2 ]

end ABC3.Found.GenEll
