import ABC3.Found.GaloisRep.OmegaCharForm
import ABC3.Found.GaloisRep.PreP4Frob
import ABC3.Found.GaloisRep.PsiRec
import ABC3.Found.GaloisRep.OmegaTwo

/-!
# Galois (G1) 第 18 ブロック —— **`omegaNum 3 = 0`(標数 2)**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★基底が `n = 0, 1, 2, 3` まで揃う

第 17 ブロックの `preΨ₄ = ψ₂⁴ + a₁Ψ₃ψ₂` を代入するだけで閉じる:

    psiComp 3 = ψ₂(preΨ₄·ψ₂⁴ + ψ₃³ + preΨ₄²)
    目標      = ψ₂(ψ₃³ + a₁·preΨ₄·ψ₂·ψ₃)

★★差は `preΨ₄·(ψ₂⁴ + preΨ₄ + a₁ψ₃ψ₂)` であり、第 17 でこれは `2ψ₂⁴` に等しい
——標数 2 で 0 である。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `psi_four` | ★`ψ₄ = preΨ₄ · ψ₂` |
| `psiComp_three` | ★★`psiComp 3` の明示形 |
| `omegaNum_three_char_two` | ★★★★**`omegaNum 3 = 0`** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve Polynomial Polynomial.Bivariate

universe u

section General

variable {R : Type u} [CommRing R] (W : WeierstrassCurve R)

/-- ★**`ψ₄ = preΨ₄ · ψ₂`**。 -/
theorem psi_four : W.ψ 4 = C W.preΨ₄ * W.ψ₂ := by
  norm_num [WeierstrassCurve.ψ, normEDS, Int.odd_iff, Int.even_iff]

/-- ★★**`psiComp 3` の明示形**。 -/
theorem psiComp_three :
    psiComp W 3 = W.ψ 5 * W.ψ₂ - C W.preΨ₄ ^ 2 * W.ψ₂ := by
  have h5 : W.ψ 5 = preNormEDS (W.ψ₂ ^ 4) (C W.Ψ₃) (C W.preΨ₄) 5 := by
    norm_num [WeierstrassCurve.ψ, normEDS, Int.odd_iff, Int.even_iff]
  rw [psiComp, complEDS₂_three, h5]

end General

section CharTwo

variable {R : Type u} [CommRing R] [CharP R 2] (W : WeierstrassCurve R)

/-- ★★★★**`omegaNum 3 = 0`(標数 2)**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★★★第 17 ブロックの Frobenius 表示を代入するだけで閉じる。 -/
theorem omegaNum_three_char_two : omegaNum W 3 = 0 := by
  have hR : (2 : R) = 0 := CharP.cast_eq_zero R 2
  have h2 : (2 : R[X][Y]) = 0 := by
    rw [show (2 : R[X][Y]) = C (C (2 : R)) by rw [map_ofNat, map_ofNat], hR, map_zero, map_zero]
  have hd : C W.preΨ₄ = W.ψ₂ ^ 4 + C (C W.a₁) * C W.Ψ₃ * W.ψ₂ := by
    rw [psi2_char_two, ← map_pow, ← map_mul, ← map_mul, ← map_add]
    exact congrArg C (preP4_frob_char_two W)
  have h5 : W.ψ 5 = C W.preΨ₄ * W.ψ₂ * W.ψ₂ ^ 3 - C W.Ψ₃ ^ 3 := by
    have h := psi_odd W 2
    norm_num at h
    exact h
  rw [omegaNum_eq_zero_iff_char_two]
  rw [show (3 : ℤ) + 1 = 4 by norm_num, show (3 : ℤ) - 1 = 2 by norm_num]
  rw [psiComp_three, h5, psi_three, psi_four, psi_two]
  linear_combination (-(C W.preΨ₄ * W.ψ₂)) * hd
    + (-(C W.Ψ₃ ^ 3 * W.ψ₂) - C (C W.a₁) * C W.preΨ₄ * W.ψ₂ ^ 2 * C W.Ψ₃) * h2

end CharTwo

/-! ## ★出典の紐付け(`.src`) -/

def omegaNum_three_char_two.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(ω_n の分子が標数 2 で消えること——n = 3)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
