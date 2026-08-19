import ABC3.Found.GaloisRep.OmegaIdentity
import ABC3.Found.GaloisRep.CharTwo

/-!
# Galois (G1) 第 16 ブロック —— **標数 2 での目標式は `ψ` と `a₁` だけで書ける**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★目標を**帰納向きの形**に直す

G1 の律速は `omegaNum n = 0`(標数 2)であった。第 10 ブロックの `omegaNum_eq_zero_iff'` は

    psiComp n = ψₙ · (a₁ · (X ψₙ² − ψ₍ₙ₊₁₎ψ₍ₙ₋₁₎) + a₃ ψₙ²)

という形だが、★★標数 2 では **`ψ₂ = a₁X + a₃`** なので `X` と `a₃` が `ψ₂` に畳める:

    psiComp n = ψₙ · (ψ₂ · ψₙ² + a₁ · ψ₍ₙ₊₁₎ · ψ₍ₙ₋₁₎)

★★★右辺が **`ψ` と `a₁` だけ**になった——`normEDS` の漸化式で帰納するのに向いた形である。

## ★★これがなぜ効くか

| 前の形 | 新しい形 |
|---|---|
| `X`・`a₁`・`a₃` の 3 つが混じる | ★`ψ₂` と `a₁` の 2 つだけ |
| `φₙ` を経由する | ★`ψ` だけ |

★`complEDS₂` の定義は `preNormEDS` の多項式式なので、
**両辺とも `preNormEDS` の言葉で閉じている**。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `omegaNum_eq_zero_iff_char_two` | ★★★★**標数 2 での目標式**(`ψ` と `a₁` だけ) |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve Polynomial Polynomial.Bivariate

universe u

variable {R : Type u} [CommRing R] [CharP R 2] (W : WeierstrassCurve R)

/-- ★★★★**標数 2 での目標式**——`X` と `a₃` は `ψ₂` に畳める。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★★★これが `normEDS` の漸化式で帰納するための形である。 -/
theorem omegaNum_eq_zero_iff_char_two (n : ℤ) :
    omegaNum W n = 0
      ↔ psiComp W n
        = W.ψ n * (W.ψ₂ * W.ψ n ^ 2 + C (C W.a₁) * (W.ψ (n + 1) * W.ψ (n - 1))) := by
  rw [omegaNum_eq_zero_iff' W n, psi2_char_two]
  have hsub : ∀ x y : R[X][Y], x - y = x + y := fun x y => by
    have h2R : (2 : R) = 0 := CharP.cast_eq_zero R 2
    have h2 : (2 : R[X][Y]) = 0 := by
      rw [show (2 : R[X][Y]) = C (C (2 : R)) by rw [map_ofNat, map_ofNat],
        h2R, map_zero, map_zero]
    have hyy : y + y = 0 := by rw [← two_mul, h2, zero_mul]
    rw [sub_eq_add_neg, neg_eq_of_add_eq_zero_left hyy]
  constructor
  · intro h
    rw [h]
    rw [map_add, map_mul]
    rw [hsub]
    ring
  · intro h
    rw [h]
    rw [map_add, map_mul]
    rw [hsub]
    ring

/-! ## ★出典の紐付け(`.src`) -/

def omegaNum_eq_zero_iff_char_two.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(ω_n の分子が標数 2 で消えること——帰納向きの形)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
