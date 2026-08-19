import ABC3.Found.GaloisRep.BInvChar2

/-!
# Galois (G1) 第 13 ブロック —— **標数 2 での `Ψ₃`**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★恒等式の展開に要る 1 本目

    Ψ₃ = 3X⁴ + b₂X³ + 3b₄X² + 3b₆X + b₈
       = X⁴ + b₂X³ + b₄X² + b₆X + b₈          （標数 2、3 = 1)

★§9-277 の手計算で使った形である。

## ★★摩擦 —— `(3 : R[X]) = C (3 : R)` も `simp` で出ない

★第 6 の `(2 : R[X][Y]) = C (C 2)` と**同じ症状**である。
★★`rw [map_ofNat]` で通る——**数値リテラルの持ち上げは `map_ofNat`** が定石。

★★★この摩擦は**この区間で 2 度目**なので、
[[ring-instance-two-paths]] に「**数値リテラル**」の欄を足すべきである。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `Psi3_char_two` | ★★標数 2 での `Ψ₃` |
-/

universe u

namespace ABC3.Found.GaloisRep

open WeierstrassCurve Polynomial

variable {R : Type u} [CommRing R] [CharP R 2] (W : WeierstrassCurve R)

/-- ★★**標数 2 での `Ψ₃ = X⁴ + b₂X³ + b₄X² + b₆X + b₈`**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★§9-277 の手計算で使った形である。 -/
theorem Psi3_char_two :
    W.Ψ₃ = X ^ 4 + C W.b₂ * X ^ 3 + C W.b₄ * X ^ 2 + C W.b₆ * X + C W.b₈ := by
  rw [WeierstrassCurve.Ψ₃]
  have h : (2 : R) = 0 := CharP.cast_eq_zero R 2
  have h3 : (3 : R) = 1 := by
    have : (3 : R) = 2 + 1 := by norm_num
    rw [this, h, zero_add]
  have h3p : (3 : R[X]) = 1 := by
    rw [show (3 : R[X]) = C (3 : R) by rw [map_ofNat], h3, map_one]
  rw [h3p]
  ring

/-! ## ★出典の紐付け(`.src`) -/

def Psi3_char_two.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(G1——標数 2 での Ψ₃)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
