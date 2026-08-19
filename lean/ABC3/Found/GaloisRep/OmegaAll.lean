import ABC3.Found.GaloisRep.EdsAll
import ABC3.Found.GaloisRep.OmegaCharForm
import ABC3.Found.GaloisRep.OmegaThree
import Mathlib.Algebra.CharP.Algebra

/-!
# Galois (G1) 第 26 ブロック —— **★★★★★★★★`omegaNum = 0`(標数 2、すべての `n`)**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★mathlib の TODO が閉じた

mathlib の `DivisionPolynomial/Basic.lean` は `ω_n` を **TODO** と書き、
FLT も `n_torsion_card` を sorry のままにしている。
★その律速が **`2 ∣ ω_n の分子`**、すなわち標数 2 で `omegaNum = 0` であった。

## ★★第 25 の `T(k)` を `omegaNum` に翻訳する

| `k` の偶奇 | `complEDS₂` の因子 | `ψ` の因子 | 変換 |
|---|---|---|---|
| 偶 | `1` | `ψ₂` | ★そのまま |
| 奇 | `ψ₂` | `1` | ★★`ψ₂` を掛ける |

★★★どちらも**割り算が要らない**——普遍曲線を経由せずに任意の標数 2 の環で出る。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `omegaNum_char_two_nat` | ★★★自然数の `n` で `omegaNum = 0` |
| `omegaNum_char_two_all` | ★★★★★★★★**すべての `n : ℤ` で `omegaNum = 0`** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve Polynomial Polynomial.Bivariate

open scoped Classical

universe u

variable {R : Type u} [CommRing R] [CharP R 2] (W : WeierstrassCurve R)

noncomputable instance : CharP (Polynomial R) 2 :=
  charP_of_injective_ringHom (Polynomial.C_injective) 2

noncomputable instance : CharP (R[X][Y]) 2 :=
  charP_of_injective_ringHom (Polynomial.C_injective) 2

/-- ★`T(k)` の曲線への当てはめ。 -/
theorem edsTarget_curve (n : ℕ) :
    EdsTarget W.ψ₂ (C W.Ψ₃) (C W.preΨ₄) (C (C W.a₁)) (n : ℤ) := by
  refine edsTarget_nat W.ψ₂ (C W.Ψ₃) (C W.preΨ₄) (C (C W.a₁)) ?_ n
  rw [psi2_char_two, ← map_pow, ← map_mul, ← map_mul, ← map_add]
  exact congrArg C (preP4_frob_char_two W)

/-- ★★★**自然数の `n` で `omegaNum = 0`**。 -/
theorem omegaNum_char_two_nat (n : ℕ) : omegaNum W (n : ℤ) = 0 := by
  have hz : (2 : R[X][Y]) = 0 := CharP.cast_eq_zero _ 2
  have hT := edsTarget_curve W n
  simp only [EdsTarget, edsEps, edsP] at hT
  rw [omegaNum_eq_zero_iff_char_two]
  show complEDS₂ W.ψ₂ (C W.Ψ₃) (C W.preΨ₄) (n : ℤ) = _
  rw [complEDS₂]
  simp only [WeierstrassCurve.ψ, normEDS]
  by_cases hn : Even (n : ℤ)
  · rw [if_pos hn] at hT ⊢
    rw [if_pos hn, if_neg (by simpa [Int.even_add_one] using hn),
      if_neg (by simpa [Int.even_sub] using hn)]
    linear_combination hT
      + (-(preNormEDS (W.ψ₂ ^ 4) (C W.Ψ₃) (C W.preΨ₄) ((n : ℤ) - 2)
          * preNormEDS (W.ψ₂ ^ 4) (C W.Ψ₃) (C W.preΨ₄) ((n : ℤ) + 1) ^ 2)) * hz
  · rw [if_neg hn] at hT ⊢
    rw [if_neg hn, if_pos (by simpa [Int.even_add_one] using hn),
      if_pos (by simpa [Int.even_sub] using hn)]
    linear_combination W.ψ₂ * hT
      + (-(W.ψ₂ * preNormEDS (W.ψ₂ ^ 4) (C W.Ψ₃) (C W.preΨ₄) ((n : ℤ) - 2)
          * preNormEDS (W.ψ₂ ^ 4) (C W.Ψ₃) (C W.preΨ₄) ((n : ℤ) + 1) ^ 2)) * hz

/-- ★`ω_n` の分子は符号反転で不変(標数 2)。 -/
theorem omegaNum_neg_char_two (n : ℤ) : omegaNum W (-n) = omegaNum W n := by
  have hz : (2 : R[X][Y]) = 0 := CharP.cast_eq_zero _ 2
  have hpsiN : ∀ k : ℤ, W.ψ (-k) = -W.ψ k := fun k => normEDS_neg _ _ _ k
  have hpsi : W.ψ (-n) = -W.ψ n := hpsiN n
  have hcomp : psiComp W (-n) = psiComp W n := complEDS₂_neg _ _ _ n
  have hphi : W.φ (-n) = W.φ n := by
    rw [phi_def, phi_def, hpsi]
    rw [show -n + 1 = -(n - 1) by ring, show -n - 1 = -(n + 1) by ring,
      hpsiN, hpsiN]
    ring
  show psiComp W (-n) - W.ψ (-n) * _ = psiComp W n - W.ψ n * _
  rw [hcomp, hpsi, hphi]
  linear_combination (W.ψ n * (C (C W.a₁) * W.φ n + C (C W.a₃) * W.ψ n ^ 2)) * hz

/-- ★★★★★★★★**すべての `n : ℤ` で `omegaNum = 0`**(標数 2)。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★★★これが mathlib の `ω_n` の TODO の律速であった。 -/
theorem omegaNum_char_two_all (n : ℤ) : omegaNum W n = 0 := by
  obtain ⟨m, hm | hm⟩ := n.eq_nat_or_neg
  · rw [hm]; exact omegaNum_char_two_nat W m
  · rw [hm, omegaNum_neg_char_two]; exact omegaNum_char_two_nat W m

/-! ## ★出典の紐付け(`.src`) -/

def edsTarget_curve.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(ω_n の分子——曲線への当てはめ)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
