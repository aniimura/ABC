import ABC3.Found.GaloisRep.ModTwo

/-!
# Galois (G1) 第 6 ブロック —— **標数 2 での計算**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★標数 2 では `2Y` が消える

    ψ₂ = 2Y + (a₁X + a₃)   ⟹   ψ₂ = a₁X + a₃      （標数 2)

★これが `omegaNum` が標数 2 で消える筋の出発点である。

## ★★初期値も消える

`omegaNum 0 = 2` なので、標数 2 では **`omegaNum 0 = 0`** である。

## ★★摩擦 —— `(2 : R[X][Y]) = C (C 2)` は `simp` で出ない

★`simp` は「no progress」と言う。★★`rw [map_ofNat, map_ofNat]` で
**`C` が `ofNat` を保つ**ことを 2 段使うと通る。

★★★数値リテラルの持ち上げは、**`simp` より `map_ofNat` のほうが確実**である
——多段の `C` では `simp` が止まる。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `psi2_char_two` | ★★標数 2 では `ψ₂ = a₁X + a₃` |
| `omegaNum_zero_char_two` | ★標数 2 では `omegaNum 0 = 0` |
-/

universe u

namespace ABC3.Found.GaloisRep

open WeierstrassCurve Polynomial Polynomial.Bivariate

variable {R : Type u} [CommRing R] [CharP R 2] (W : WeierstrassCurve R)

/-- ★★**標数 2 では `ψ₂ = a₁X + a₃`**——`2Y` が消える。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★★★これが `omegaNum` が標数 2 で消える筋の出発点である。 -/
theorem psi2_char_two : W.ψ₂ = C (C W.a₁ * X + C W.a₃) := by
  rw [psi2_val]
  have : (2 : R) = 0 := CharP.cast_eq_zero R 2
  rw [this]
  simp

/-- ★**標数 2 では `omegaNum 0 = 0`**。 -/
theorem omegaNum_zero_char_two : omegaNum W 0 = 0 := by
  rw [omegaNum_zero]
  have h : (2 : R) = 0 := CharP.cast_eq_zero R 2
  have h2 : (2 : R[X][Y]) = C (C (2 : R)) := by
    rw [map_ofNat, map_ofNat]
  rw [h2, h, map_zero, map_zero]

/-! ## ★出典の紐付け(`.src`) -/

def psi2_char_two.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(G1——標数 2 での計算)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
