import ABC3.Found.GaloisRep.UniversalF2
import ABC3.Found.GaloisRep.OmegaThree

/-!
# Galois (G1) 第 20 ブロック —— **`omegaNum 4 = 0`(普遍曲線、標数 2)**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★`normEDSRec` の基底 `P 0`〜`P 4` が揃った

mathlib の **`normEDSRec`** は `P 0, P 1, P 2, P 3, P 4` と偶奇の帰納段を要求する。
★第 13–15・18 で `n = 0, 1, 2, 3` が済んでいたので、残るは `P 4` だった。

## ★★`ψ₂` を 2 回掛けた形で示す

`complEDS₂` の 4 の値は `preNormEDS … 6` を含み、`ψ₆` に直すのに `ψ₂` の約分が要る。
★そこで **`ψ₂²` を掛けた形**(`key4`)を任意の標数 2 の環で示し、
★★**普遍曲線の上でだけ約す**(第 19 の `uCurveF2_psi2_ne_zero`)。

## ★★★手で作った証明証明書

`ψ₅ = preΨ₄·ψ₂⁴ − Ψ₃³`、`ψ₆ψ₂ = ψ₂²Ψ₃ψ₅ − Ψ₃(preΨ₄ψ₂)²` を代入し、
第 17 の `preΨ₄ = ψ₂⁴ + a₁Ψ₃ψ₂` で `preΨ₄` を消すと、差は

    2·(ψ₂¹⁰Ψ₃³ + a₁ψ₂⁷Ψ₃⁴ − ψ₂²Ψ₃⁶ − ψ₂¹⁸
        − 3a₁Ψ₃ψ₂¹⁵ − 3a₁²Ψ₃²ψ₂¹² − a₁³Ψ₃³ψ₂⁹)

になる。★**この係数を手で計算して `linear_combination` に渡した**——一発で通った。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `key4` | ★★★`ψ₂²` を掛けた形(任意の標数 2 の環) |
| `omegaNum_four_char_two_univ` | ★★★★**`omegaNum 4 = 0`**(普遍曲線) |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve Polynomial Polynomial.Bivariate

universe u

variable {R : Type u} [CommRing R] [CharP R 2] (W : WeierstrassCurve R)

/-- ★★★**`ψ₂²` を掛けた形**——任意の標数 2 の環で成り立つ。 -/
theorem key4 :
    psiComp W 4 * W.ψ₂ ^ 2
      = (W.ψ 4 * (W.ψ₂ * W.ψ 4 ^ 2
          + C (C W.a₁) * (W.ψ 5 * W.ψ 3))) * W.ψ₂ ^ 2 := by
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
  have h6 : W.ψ 6 * W.ψ₂
      = W.ψ₂ ^ 2 * C W.Ψ₃ * W.ψ 5 - C W.Ψ₃ * (C W.preΨ₄ * W.ψ₂) ^ 2 := by
    have h := psi_even W 3
    norm_num at h
    exact h
  have h4 : psiComp W 4 * W.ψ₂ = C W.Ψ₃ ^ 2 * W.ψ 6 - W.ψ₂ * W.ψ 5 ^ 2 := by
    have h := psiComp_mul_psi2 W 4
    norm_num at h
    exact h
  have poly :
      C W.Ψ₃ ^ 2 * (W.ψ₂ ^ 2 * C W.Ψ₃ * (C W.preΨ₄ * W.ψ₂ * W.ψ₂ ^ 3 - C W.Ψ₃ ^ 3)
            - C W.Ψ₃ * (C W.preΨ₄ * W.ψ₂) ^ 2)
          - W.ψ₂ ^ 2 * (C W.preΨ₄ * W.ψ₂ * W.ψ₂ ^ 3 - C W.Ψ₃ ^ 3) ^ 2
        = (C W.preΨ₄ * W.ψ₂ * (W.ψ₂ * (C W.preΨ₄ * W.ψ₂) ^ 2
            + C (C W.a₁) * ((C W.preΨ₄ * W.ψ₂ * W.ψ₂ ^ 3 - C W.Ψ₃ ^ 3) * C W.Ψ₃)))
          * W.ψ₂ ^ 2 := by
    rw [hd]
    linear_combination
      (W.ψ₂ ^ 10 * C W.Ψ₃ ^ 3 + C (C W.a₁) * W.ψ₂ ^ 7 * C W.Ψ₃ ^ 4
        - W.ψ₂ ^ 2 * C W.Ψ₃ ^ 6 - W.ψ₂ ^ 18
        - 3 * C (C W.a₁) * C W.Ψ₃ * W.ψ₂ ^ 15
        - 3 * C (C W.a₁) ^ 2 * C W.Ψ₃ ^ 2 * W.ψ₂ ^ 12
        - C (C W.a₁) ^ 3 * C W.Ψ₃ ^ 3 * W.ψ₂ ^ 9) * h2
  calc psiComp W 4 * W.ψ₂ ^ 2
      = (psiComp W 4 * W.ψ₂) * W.ψ₂ := by ring
    _ = (C W.Ψ₃ ^ 2 * W.ψ 6 - W.ψ₂ * W.ψ 5 ^ 2) * W.ψ₂ := by rw [h4]
    _ = C W.Ψ₃ ^ 2 * (W.ψ 6 * W.ψ₂) - W.ψ₂ ^ 2 * W.ψ 5 ^ 2 := by ring
    _ = _ := by rw [h6, h5, psi_four, psi_three]; exact poly


/-- ★★★★**`omegaNum 4 = 0`(普遍曲線)**——`ψ₂ ≠ 0` なので約せる。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★★★これで `normEDSRec` の基底 `P 0`〜`P 4` が揃った。 -/
theorem omegaNum_four_char_two_univ : omegaNum uCurveF2 4 = 0 := by
  haveI : Fact (Nat.Prime 2) := ⟨Nat.prime_two⟩
  haveI : IsDomain UF2Ring := inferInstance
  haveI : IsDomain (Polynomial UF2Ring) := inferInstance
  haveI : IsDomain (Polynomial (Polynomial UF2Ring)) := inferInstance
  rw [omegaNum_eq_zero_iff_char_two]
  rw [show (4 : ℤ) + 1 = 5 by norm_num, show (4 : ℤ) - 1 = 3 by norm_num]
  exact mul_right_cancel₀ (pow_ne_zero 2 uCurveF2_psi2_ne_zero) (key4 uCurveF2)

/-! ## ★出典の紐付け(`.src`) -/

def omegaNum_four_char_two_univ.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(ω_n の分子が標数 2 で消えること——n = 4)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
