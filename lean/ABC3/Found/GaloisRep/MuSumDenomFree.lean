/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.TateDSeries
import ABC3.Found.GenEll.CyclotomicDenomFree
import ABC3.Meta.Claim

/-!
# 第 1097 ブロック —— **μ_l の頭項の和を分母なしで書く**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★これは何か

第 1092 の `natCast_pow_mul_sum_inverse` を、Tate の**頭項**に当てる。

☆頭項はすべて `(t の多項式) * Ring.inverse (1 − t)^k` の形をしている:

| 項 | 分子 | `k` |
|---|---|---|
| `tateXterm` | `t` | 2 |
| `tateYterm` | `t²` | 3 |
| `tateDXterm` | `t(1+t)` | 3 |
| `tateDYterm` | `t²(2+t)` | 4 |
| `tateD2Xterm` | `t(1+4t+t²)` | 4 |
| `tateD3Xterm` | `t(1+11t+11t²+t³)` | 5 |

★したがって `l^k` を掛ければ `Ring.inverse` が消え、
`∑_i (分子) · ∏_{j≠i}(1 − ζ^j)^k` という**分母を含まない形**になる。

★★**これが `p ∣ l`（`1 − ζ^i` が単元でない）でも意味を持つ形である。**
☆尾（`tateXtail` 等）は `q^n u ∈ 𝔪` なので `1 − q^n u` が単元であり、はじめから問題ない。
-/

namespace ABC3.Found.GaloisRep

open Finset ABC3.Meta ABC3.Found.GenEll

variable {R : Type} [CommRing R] [IsDomain R] {l : ℕ} {ζ : R}

/-- ★★★★★★★★**`X` の頭項の和を分母なしで**（第 1097）。 -/
theorem natCast_pow_mul_sum_tateXterm (hl : 0 < l) (hζ : IsPrimitiveRoot ζ l)
    (hu : ∀ i ∈ (range l).erase 0, IsUnit (1 - ζ ^ i)) :
    (l : R) ^ 2 * ∑ i ∈ (range l).erase 0, tateXterm (ζ ^ i)
      = ∑ i ∈ (range l).erase 0,
          ζ ^ i * ∏ j ∈ ((range l).erase 0).erase i, (1 - ζ ^ j) ^ 2 :=
  natCast_pow_mul_sum_inverse hl hζ 2 hu (fun i => ζ ^ i)

/-- ★★★★★★**`Y` の頭項の和を分母なしで**（第 1097）。 -/
theorem natCast_pow_mul_sum_tateYterm (hl : 0 < l) (hζ : IsPrimitiveRoot ζ l)
    (hu : ∀ i ∈ (range l).erase 0, IsUnit (1 - ζ ^ i)) :
    (l : R) ^ 3 * ∑ i ∈ (range l).erase 0, tateYterm (ζ ^ i)
      = ∑ i ∈ (range l).erase 0,
          (ζ ^ i) ^ 2 * ∏ j ∈ ((range l).erase 0).erase i, (1 - ζ ^ j) ^ 3 :=
  natCast_pow_mul_sum_inverse hl hζ 3 hu (fun i => (ζ ^ i) ^ 2)

/-- ★★★★★★**`DX` の頭項の和を分母なしで**（第 1097）。 -/
theorem natCast_pow_mul_sum_tateDXterm (hl : 0 < l) (hζ : IsPrimitiveRoot ζ l)
    (hu : ∀ i ∈ (range l).erase 0, IsUnit (1 - ζ ^ i)) :
    (l : R) ^ 3 * ∑ i ∈ (range l).erase 0, tateDXterm (ζ ^ i)
      = ∑ i ∈ (range l).erase 0,
          (ζ ^ i * (1 + ζ ^ i)) * ∏ j ∈ ((range l).erase 0).erase i, (1 - ζ ^ j) ^ 3 :=
  natCast_pow_mul_sum_inverse hl hζ 3 hu (fun i => ζ ^ i * (1 + ζ ^ i))

/-- ★★★★★★**`DY` の頭項の和を分母なしで**（第 1097）。 -/
theorem natCast_pow_mul_sum_tateDYterm (hl : 0 < l) (hζ : IsPrimitiveRoot ζ l)
    (hu : ∀ i ∈ (range l).erase 0, IsUnit (1 - ζ ^ i)) :
    (l : R) ^ 4 * ∑ i ∈ (range l).erase 0, tateDYterm (ζ ^ i)
      = ∑ i ∈ (range l).erase 0,
          ((ζ ^ i) ^ 2 * (2 + ζ ^ i))
            * ∏ j ∈ ((range l).erase 0).erase i, (1 - ζ ^ j) ^ 4 :=
  natCast_pow_mul_sum_inverse hl hζ 4 hu (fun i => (ζ ^ i) ^ 2 * (2 + ζ ^ i))

/-- ★★★★★★★★**`D²X` の頭項の和を分母なしで**（第 1097）。 -/
theorem natCast_pow_mul_sum_tateD2Xterm (hl : 0 < l) (hζ : IsPrimitiveRoot ζ l)
    (hu : ∀ i ∈ (range l).erase 0, IsUnit (1 - ζ ^ i)) :
    (l : R) ^ 4 * ∑ i ∈ (range l).erase 0, tateD2Xterm (ζ ^ i)
      = ∑ i ∈ (range l).erase 0,
          (ζ ^ i * (1 + 4 * ζ ^ i + (ζ ^ i) ^ 2))
            * ∏ j ∈ ((range l).erase 0).erase i, (1 - ζ ^ j) ^ 4 :=
  natCast_pow_mul_sum_inverse hl hζ 4 hu (fun i => ζ ^ i * (1 + 4 * ζ ^ i + (ζ ^ i) ^ 2))

/-- ★★★★★★★★**`D³X` の頭項の和を分母なしで**（第 1097）。 -/
theorem natCast_pow_mul_sum_tateD3Xterm (hl : 0 < l) (hζ : IsPrimitiveRoot ζ l)
    (hu : ∀ i ∈ (range l).erase 0, IsUnit (1 - ζ ^ i)) :
    (l : R) ^ 5 * ∑ i ∈ (range l).erase 0, tateD3Xterm (ζ ^ i)
      = ∑ i ∈ (range l).erase 0,
          (ζ ^ i * (1 + 11 * ζ ^ i + 11 * (ζ ^ i) ^ 2 + (ζ ^ i) ^ 3))
            * ∏ j ∈ ((range l).erase 0).erase i, (1 - ζ ^ j) ^ 5 :=
  natCast_pow_mul_sum_inverse hl hζ 5 hu
    (fun i => ζ ^ i * (1 + 11 * ζ ^ i + 11 * (ζ ^ i) ^ 2 + (ζ ^ i) ^ 3))

/-! ## ★出典の紐付け(`.src`) -/

def natCast_pow_mul_sum_tateXterm.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(X の頭項の和を分母なしで——l² を掛ける)",
    sectionId := "genell-lemma-3-5" }

def natCast_pow_mul_sum_tateD2Xterm.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(D²X の頭項の和を分母なしで——l⁴ を掛ける)",
    sectionId := "genell-lemma-3-5" }

def natCast_pow_mul_sum_tateD3Xterm.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(D³X の頭項の和を分母なしで——l⁵ を掛ける)",
    sectionId := "genell-lemma-3-5" }

def natCast_pow_mul_sum_tateXterm.needs : List ProofObligation :=
  [ .citation "[ABC3]" "natCast_pow_mul_sum_inverse(第 1092、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.natCast_pow_mul_sum_inverse") 1,
    .implicitStep
      ("★★**2026-09-01（第 1097）**——これで頭項の側は `p ∣ l` でも意味を持つ形になった。" ++
       "☆残るのは `c4_velu_tate`・`c6_velu_tate` の恒等式を" ++
       "**右辺（分母なし）の言葉で証明し直す**ことである。" ++
       "★尾は `1 − q^n u` が単元なので問題ない。") 30 ]

end ABC3.Found.GaloisRep
