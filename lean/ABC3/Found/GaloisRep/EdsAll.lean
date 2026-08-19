import ABC3.Found.GaloisRep.EdsOddStep

/-!
# Galois (G1) 第 25 ブロック —— **`T(k)` がすべての `k` で成り立つ**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★mathlib の `normEDSRec` に基底 5 つと 2 段を渡すだけ

| 入力 | 出所 |
|---|---|
| `P 0`〜`P 4` | ★第 22 ブロック |
| 偶数段 | ★第 23 ブロック |
| 奇数段 | ★第 24 ブロック |

★★**曲線が入るのは `hd : d = b⁴ + A·c·b` ただ 1 つ**である。
-/

namespace ABC3.Found.GaloisRep

open scoped Classical

universe u

variable {R : Type u} [CommRing R] [CharP R 2] (b c d A : R)

/-- ★★★★★★**`T(k)` は自然数のすべての `k` で成り立つ**。 -/
theorem edsTarget_nat (hd : d = b ^ 4 + A * c * b) (n : ℕ) :
    EdsTarget b c d A (n : ℤ) := by
  induction n using normEDSRec with
  | zero => exact edsTarget_zero b c d A
  | one => exact_mod_cast edsTarget_one b c d A
  | two => exact_mod_cast edsTarget_two b c d A hd
  | three => exact_mod_cast edsTarget_three b c d A hd
  | four => exact_mod_cast edsTarget_four b c d A hd
  | even m _ ih2 ih3 ih4 _ =>
    push_cast at ih2 ih3 ih4 ⊢
    exact edsTarget_even_step b c d A ((m : ℤ) + 3) ih3
      (by rw [show (m : ℤ) + 3 + 1 = (m : ℤ) + 4 by ring]; exact ih4)
      (by rw [show (m : ℤ) + 3 - 1 = (m : ℤ) + 2 by ring]; exact ih2)
  | odd m ih1 ih2 ih3 _ =>
    push_cast at ih1 ih2 ih3 ⊢
    exact edsTarget_odd_step b c d A ((m : ℤ) + 2) ih2
      (by rw [show (m : ℤ) + 2 + 1 = (m : ℤ) + 3 by ring]; exact ih3)
      (by rw [show (m : ℤ) + 2 - 1 = (m : ℤ) + 1 by ring]; exact ih1)

/-! ## ★出典の紐付け(`.src`) -/

def edsTarget_nat.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(ω_n の分子——帰納法の結論)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
