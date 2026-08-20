import ABC3.Found.GaloisRep.PreP4Char2

/-!
# Galois (G1) 第 15 ブロック —— **標数 2 で `ω₂ = 0`**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

## ★★★★★★手計算(§9-277)が Lean に載った

§9-277 で手で展開した等式——

    右辺 = (a₁X+a₃)⁴ − a₁(a₁X+a₃)Ψ₃
         = b₂X⁵ + b₄X⁴ + ★(b₂b₆ + b₄²)X² + (b₂b₈+b₄b₆)X + (b₄b₈+b₆²)
    左辺 = preΨ₄
         = b₂X⁵ + b₄X⁴ +                    (b₂b₈+b₄b₆)X + (b₄b₈+b₆²)

★★★差は **`(b₂b₆ + b₄²)X²` の 1 項だけ**であり、標数 2 で
`b₂b₆ = a₁²a₃² = (a₁a₃)² = b₄²`(第 12 ブロック `b2_b6_eq_b4_sq`)なので消える。

## ★★`ring` が閉じなかった理由と逃げ道(摩擦 #8)

`b₂ → a₁²` 等を **`C` を押し込む前**に書き換えないと、`ring` は
`C (C (a₁^2))` と `C (C a₁)` を**別の原子**として扱う。順序を入れ替えたうえで、
なお残る差は

    右辺 − 左辺 = 2a₁⁴X⁴ + 6a₁³a₃X³ + 8a₁²a₃²X² + 4a₁a₃³X = 2·K

★★★`ring` は環の公理しか知らないので `2 = 0` を渡す必要がある
——`linear_combination (−K) * h₂` で一発で閉じた。

| 定理 | 内容 |
|---|---|
| `psi_two` / `psi_three` | ★`ψ₂`, `ψ₃` の値 |
| `psiComp_two` | ★`ψ^c₂ = preΨ₄` |
| `omegaNum_two_char_two` | ★★★★★★**標数 2 で `ω₂ = 0`** |
-/

universe u

namespace ABC3.Found.GaloisRep

open WeierstrassCurve Polynomial Polynomial.Bivariate

section General

variable {R : Type u} [CommRing R] (W : WeierstrassCurve R)

/-- ★`ψ 2 = ψ₂`。 -/
theorem psi_two : W.ψ 2 = W.ψ₂ := by
  rw [WeierstrassCurve.ψ]; exact normEDS_two _ _ _

/-- ★`ψ 3 = C Ψ₃`。 -/
theorem psi_three : W.ψ 3 = C W.Ψ₃ := by
  rw [WeierstrassCurve.ψ]; exact normEDS_three _ _ _

/-- ★補因子の 2 での値——`ψ^c₂ = C preΨ₄`。 -/
theorem psiComp_two : psiComp W 2 = C W.preΨ₄ := by
  rw [psiComp]; exact complEDS₂_two _ _ _

end General

section CharTwo

variable {R : Type u} [CommRing R] [CharP R 2] (W : WeierstrassCurve R)

/-- ★★★★★★**標数 2 で `ω₂ = 0`**。

★★§9-277 の手計算——両辺の差は `(b₂b₆ + b₄²)X²` の 1 項だけで、
標数 2 では `b₂b₆ = b₄²` なので消える——を Lean に載せたもの。

★★★`ring` には `2 = 0` を明示的に渡す(`linear_combination`)。 -/
theorem omegaNum_two_char_two : omegaNum W 2 = 0 := by
  rw [omegaNum_eq_zero_iff', psiComp_two, psi_two]
  norm_num
  rw [psi2_char_two, preP4_char_two, Psi3_char_two]
  rw [b2_char_two, b4_char_two, b6_char_two]
  simp only [map_add, map_mul, map_pow]
  have h : (2 : R) = 0 := CharP.cast_eq_zero R 2
  have h2 : (2 : R[X][Y]) = 0 := by
    rw [show (2 : R[X][Y]) = C (C (2:R)) by rw [map_ofNat, map_ofNat], h, map_zero, map_zero]
  have hsub : ∀ x y : R[X][Y], x - y = x + y := fun x y => by
    have hyy : y + y = 0 := by rw [← two_mul, h2, zero_mul]
    rw [sub_eq_add_neg, neg_eq_of_add_eq_zero_left hyy]
  rw [hsub]
  linear_combination (-(C (C W.a₁) ^ 4 * C X ^ 4
      + 3 * C (C W.a₁) ^ 3 * C (C W.a₃) * C X ^ 3
      + 4 * C (C W.a₁) ^ 2 * C (C W.a₃) ^ 2 * C X ^ 2
      + 2 * C (C W.a₁) * C (C W.a₃) ^ 3 * C X)) * h2

end CharTwo

/-! ## ★出典の紐付け(`.src`) -/

def omegaNum_two_char_two.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(G1——ω の 2 での消滅、標数 2)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
