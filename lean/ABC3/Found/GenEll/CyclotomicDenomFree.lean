/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.CyclotomicUnits
import ABC3.Meta.Claim

/-!
# 第 1092 ブロック —— **μ_l の和の分母払い**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★これは何か

`Lemma 3.5` の悪い素点の計算は `1 − ζ^i` の可逆性に寄りかかっている。
☆`p ∣ l` ではそれが破れる（`∏(1 − ζ^i) = l` なので `1 − ζ^i` は非単元）。

★第 277（`CollDenomFree.lean`）が共線性で使った技法——**行ごとに分母を払う**——を
μ_l の和に当てるための**土台**が本ファイルである:

    `l^k · ∑_i f_i · (1 − ζ^i)^{−k}  =  ∑_i f_i · ∏_{j ≠ i} (1 − ζ^j)^k`

☆右辺は `Ring.inverse` を含まない。したがって**右辺の形で恒等式を証明しておけば**、
`1 − ζ^i` が単元でない場合にも意味を持ち、最後に商体で `l^k` を割ればよい
（商体は標数 0 なので `l ≠ 0`）。

★★これが `d + 1 < l` を外す道の第 1 歩である。
-/

namespace ABC3.Found.GenEll

open Finset ABC3.Meta

/-- ☆`∏_{i}(1 − ζ^i)^k = l^k`。 -/
theorem prod_one_sub_pow_erase_pow {R : Type} [CommRing R] [IsDomain R] {l : ℕ} (hl : 0 < l)
    {ζ : R} (hζ : IsPrimitiveRoot ζ l) (k : ℕ) :
    ∏ i ∈ (range l).erase 0, (1 - ζ ^ i) ^ k = (l : R) ^ k := by
  rw [Finset.prod_pow, prod_one_sub_pow_erase hl hζ]

/-- ★★★★★★★★★★★★★★★★**μ_l の和の分母払い**（第 1092）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆`l^k` を掛けると `Ring.inverse (1 − ζ^i)^k` が消えて、
`∏_{j ≠ i}(1 − ζ^j)^k` という**分母を含まない**重みに変わる。
★右辺は `1 − ζ^i` が単元でなくても意味を持つ。 -/
theorem natCast_pow_mul_sum_inverse {R : Type} [CommRing R] [IsDomain R] {l : ℕ} (hl : 0 < l)
    {ζ : R} (hζ : IsPrimitiveRoot ζ l) (k : ℕ)
    (hu : ∀ i ∈ (range l).erase 0, IsUnit (1 - ζ ^ i)) (f : ℕ → R) :
    (l : R) ^ k * ∑ i ∈ (range l).erase 0, f i * Ring.inverse (1 - ζ ^ i) ^ k
      = ∑ i ∈ (range l).erase 0,
          f i * ∏ j ∈ ((range l).erase 0).erase i, (1 - ζ ^ j) ^ k := by
  rw [← prod_one_sub_pow_erase_pow hl hζ k, Finset.mul_sum]
  refine Finset.sum_congr rfl (fun i hi => ?_)
  have hsplit : ∏ j ∈ (range l).erase 0, (1 - ζ ^ j) ^ k
      = (1 - ζ ^ i) ^ k * ∏ j ∈ ((range l).erase 0).erase i, (1 - ζ ^ j) ^ k :=
    (Finset.mul_prod_erase _ _ hi).symm
  rw [hsplit]
  have hinv : (1 - ζ ^ i) ^ k * Ring.inverse (1 - ζ ^ i) ^ k = 1 := by
    rw [← mul_pow, Ring.mul_inverse_cancel _ (hu i hi), one_pow]
  have hgoal : ((1 - ζ ^ i) ^ k * ∏ j ∈ ((range l).erase 0).erase i, (1 - ζ ^ j) ^ k)
        * (f i * Ring.inverse (1 - ζ ^ i) ^ k)
      = ((1 - ζ ^ i) ^ k * Ring.inverse (1 - ζ ^ i) ^ k) *
          (f i * ∏ j ∈ ((range l).erase 0).erase i, (1 - ζ ^ j) ^ k) := by ring
  rw [hgoal, hinv, one_mul]

/-! ## ★出典の紐付け(`.src`) -/

def prod_one_sub_pow_erase_pow.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(∏(1 − ζ^i)^k = l^k。★無条件)",
    sectionId := "genell-lemma-3-5" }

def natCast_pow_mul_sum_inverse.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(μ_l の和の分母払い——l^k を掛けると Ring.inverse が消える)",
    sectionId := "genell-lemma-3-5" }

def natCast_pow_mul_sum_inverse.needs : List ProofObligation :=
  [ .citation "[ABC3]" "prod_one_sub_pow_erase(∏(1 − ζ^i) = l、第 950 系、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.prod_one_sub_pow_erase") 1,
    .implicitStep
      ("★★**2026-09-01（第 1092）**——これは `d + 1 < l` を外す道の第 1 歩である。" ++
       "☆右辺（分母を含まない形）で `c4_velu_tate`・`c6_velu_tate` の恒等式を" ++
       "証明し直せば、`1 − ζ^i` が単元でない `p ∣ l` の悪い素点でも意味を持つ。" ++
       "★最後に商体で `l^k` を割る（商体は標数 0 なので `l ≠ 0`）。" ++
       "☆残りの規模は 40-80 ブロック。") 40 ]

end ABC3.Found.GenEll
