import ABC3.Found.GaloisRep.Psi3Char2

/-!
# Galois (G1) 第 14 ブロック —— **標数 2 での `preΨ₄`**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★恒等式の展開に要る 2 本目

    preΨ₄ = 2X⁶ + b₂X⁵ + 5b₄X⁴ + 10b₆X³ + 10b₈X² + (b₂b₈−b₄b₆)X + (b₄b₈−b₆²)
          = b₂X⁵ + b₄X⁴ + (b₂b₈+b₄b₆)X + (b₄b₈+b₆²)      （標数 2)

★§9-277 の手計算で「左辺」として使った形である。

## ★★標数 2 の数値は 3 つ

| リテラル | 標数 2 |
|---|---|
| `2`, `10` | ★`0` |
| `5` | ★`1` |

★どれも `rw [map_ofNat]` で `C` の中に落としてから計算する。

## ★★★引き算を足し算に直す 1 行が要った

`ring` は `b₂b₈ − b₄b₆` と `b₂b₈ + b₄b₆` を**同一視しない**。
★標数 2 では `x − y = x + y` なので、その補題を `have` で作って `rw` する:

    hsub : ∀ x y : R, x - y = x + y

★★`neg_eq_of_add_eq_zero_left` に `y + y = 0`(`two_mul` + 標数 2)を渡す。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `preP4_char_two` | ★★標数 2 での `preΨ₄` |
-/

universe u

namespace ABC3.Found.GaloisRep

open WeierstrassCurve Polynomial

variable {R : Type u} [CommRing R] [CharP R 2] (W : WeierstrassCurve R)

/-- ★★**標数 2 での `preΨ₄`**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★§9-277 の手計算で「左辺」として使った形である。 -/
theorem preP4_char_two :
    W.preΨ₄ = C W.b₂ * X ^ 5 + C W.b₄ * X ^ 4
      + C (W.b₂ * W.b₈ + W.b₄ * W.b₆) * X + C (W.b₄ * W.b₈ + W.b₆ ^ 2) := by
  rw [WeierstrassCurve.preΨ₄]
  have h : (2 : R) = 0 := CharP.cast_eq_zero R 2
  have h2p : (2 : R[X]) = 0 := by rw [show (2 : R[X]) = C (2 : R) by rw [map_ofNat], h, map_zero]
  have h5 : (5 : R[X]) = 1 := by
    rw [show (5 : R[X]) = C (5 : R) by rw [map_ofNat]]
    have : (5 : R) = 2 * 2 + 1 := by norm_num
    rw [this, h, mul_zero, zero_add, map_one]
  have h10 : (10 : R[X]) = 0 := by
    rw [show (10 : R[X]) = C (10 : R) by rw [map_ofNat]]
    have : (10 : R) = 2 * 5 := by norm_num
    rw [this, h, zero_mul, map_zero]
  rw [h2p, h5, h10]
  have hsub : ∀ x y : R, x - y = x + y := fun x y => by
    have hyy : y + y = 0 := by rw [← two_mul, h, zero_mul]
    rw [sub_eq_add_neg, neg_eq_of_add_eq_zero_left hyy]
  rw [hsub, hsub]
  ring

/-! ## ★出典の紐付け(`.src`) -/

def preP4_char_two.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(G1——標数 2 での preΨ₄)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
