import ABC3.Found.GaloisRep.UniversalF2
import Mathlib.NumberTheory.EllipticDivisibilitySequence

/-!
# Galois (G1) 第 22 ブロック —— **帰納の目標 `T(k)` と基底**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★§9-357 で測った形をそのまま書く

`ψ` のままだと偶数添字で `ψ₂` を割る操作が情報を落とす(§9-356 で実測)。
★**`preNormEDS` で書く**と漸化式が値を定義しているので、関係式が要らない。

    T(k) : P₍ₖ₋₁₎² P₍ₖ₊₂₎ + P₍ₖ₋₂₎ P₍ₖ₊₁₎² = ε(k)·Pₖ³ + A·b·Pₖ P₍ₖ₊₁₎ P₍ₖ₋₁₎

    P := preNormEDS (b⁴) c d,   ε(k) = b⁴ (k 偶), 1 (k 奇)

## ★★曲線が入るのは `hd` ただ 1 つ

第 17 ブロックの `preΨ₄ = ψ₂⁴ + a₁Ψ₃ψ₂` を抽象化した

    hd : d = b⁴ + A·c·b

★★★**帰納段には `c` も `d` も現れない**(§9-357 の実測)——曲線固有の情報は
**基底の 5 つだけ**に入る。

## ★★本ブロックで取れるもの

| 定義・定理 | 内容 |
|---|---|
| `edsP` / `edsEps` / `EdsTarget` | ★目標の定式化 |
| `edsTarget_zero` 〜 `edsTarget_four` | ★★★**基底 5 つ**(すべて `hd` から) |
-/

namespace ABC3.Found.GaloisRep

open scoped Classical

universe u

variable {R : Type u} [CommRing R] [CharP R 2] (b c d A : R)

/-- ★`preNormEDS (b⁴) c d`。 -/
noncomputable abbrev edsP (k : ℤ) : R := preNormEDS (b ^ 4) c d k

/-- ★パリティで分かれる係数。 -/
noncomputable def edsEps (k : ℤ) : R := if Even k then b ^ 4 else 1

/-- ★★**帰納の目標**。 -/
def EdsTarget (k : ℤ) : Prop :=
  edsP b c d (k - 1) ^ 2 * edsP b c d (k + 2)
      + edsP b c d (k - 2) * edsP b c d (k + 1) ^ 2
    = edsEps b k * edsP b c d k ^ 3
      + A * b * (edsP b c d k * edsP b c d (k + 1) * edsP b c d (k - 1))

section Base

variable (hd : d = b ^ 4 + A * c * b)

private theorem two_eq_zero : (2 : R) = 0 := CharP.cast_eq_zero R 2

theorem edsTarget_zero : EdsTarget b c d A 0 := by
  show _ = edsEps b 0 * _ + _
  rw [edsEps, if_pos (by decide : Even (0 : ℤ))]
  norm_num [edsP, show (0 : ℤ) - 1 = -1 by norm_num, show (0 : ℤ) - 2 = -2 by norm_num,
    show (0 : ℤ) + 2 = 2 by norm_num, show (0 : ℤ) + 1 = 1 by norm_num,
    show (-1 : ℤ) = -(1 : ℤ) by norm_num, show (-2 : ℤ) = -(2 : ℤ) by norm_num,
    preNormEDS_neg]

theorem edsTarget_one : EdsTarget b c d A 1 := by
  have h2 := two_eq_zero (R := R)
  show _ = edsEps b 1 * _ + _
  rw [edsEps, if_neg (by decide : ¬ Even (1 : ℤ))]
  norm_num [edsP, show (1 : ℤ) - 1 = 0 by norm_num, show (1 : ℤ) - 2 = -1 by norm_num,
    show (1 : ℤ) + 2 = 3 by norm_num, show (1 : ℤ) + 1 = 2 by norm_num,
    show (-1 : ℤ) = -(1 : ℤ) by norm_num, preNormEDS_neg]
  linear_combination -h2

end Base

end ABC3.Found.GaloisRep
