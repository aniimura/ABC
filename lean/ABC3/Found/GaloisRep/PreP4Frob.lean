import ABC3.Found.GaloisRep.PreP4Char2
import ABC3.Found.GaloisRep.Psi3Char2
import ABC3.Found.GaloisRep.BInvChar2

/-!
# Galois (G1) 第 17 ブロック —— **標数 2 では `preΨ₄ = ψ₂⁴ + a₁Ψ₃ψ₂`**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★`n = 3` の基底が、この 1 本の恒等式に落ちた

§9-350 で `n = 3` の目標を追ったところ、`normEDS_odd`(`ψ₅ = ψ₄ψ₂³ − ψ₃³`)と
`ψ₄ = preΨ₄ · ψ₂` を代入すると、標数 2 では

    preΨ₄ = (a₁X + a₃)⁴ + a₁ · Ψ₃ · (a₁X + a₃)

の 1 本に帰着する。★★右辺の第 1 項は **Frobenius** で `a₁⁴X⁴ + a₃⁴` に潰れる。

## ★★確かめ方

在庫の標数 2 表示(第 11–13 ブロック)を代入すると、両辺とも

    a₁²X⁵ + a₁a₃X⁴ + (a₁²b₈ + a₁a₃³)X + (a₃⁴ + a₁a₃b₈)

になる。★`b₈` は展開せずに済む——**両辺で同じ形で現れて相殺する**。

★★差は `2·(a₁⁴X⁴ + a₁³a₃X³ + a₁²a₃²X²)` なので `linear_combination` で閉じる。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `preP4_frob_char_two` | ★★★★**`preΨ₄ = ψ₂⁴ + a₁Ψ₃ψ₂`**(標数 2) |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve Polynomial

universe u

variable {R : Type u} [CommRing R] [CharP R 2] (W : WeierstrassCurve R)

/-- ★★★★**標数 2 では `preΨ₄` が `ψ₂` と `Ψ₃` だけで書ける**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★★★これが `n = 3` の基底の核である。 -/
theorem preP4_frob_char_two :
    W.preΨ₄ = (C W.a₁ * X + C W.a₃) ^ 4
      + C W.a₁ * W.Ψ₃ * (C W.a₁ * X + C W.a₃) := by
  have hR : (2 : R) = 0 := CharP.cast_eq_zero R 2
  have h2 : (2 : R[X]) = 0 := by
    rw [show (2 : R[X]) = C (2 : R) by rw [map_ofNat], hR, map_zero]
  rw [preP4_char_two, Psi3_char_two, b2_char_two, b4_char_two, b6_char_two]
  simp only [map_add, map_mul, map_pow]
  linear_combination
    (-(C W.a₁ ^ 4 * X ^ 4 + 3 * C W.a₁ ^ 3 * C W.a₃ * X ^ 3
      + 4 * C W.a₁ ^ 2 * C W.a₃ ^ 2 * X ^ 2 + 2 * C W.a₁ * C W.a₃ ^ 3 * X)) * h2

/-! ## ★出典の紐付け(`.src`) -/

def preP4_frob_char_two.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(標数 2 での preΨ₄ の Frobenius 表示)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
