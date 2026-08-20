import ABC3.Found.GaloisRep.CoprimePhiPSq
import Mathlib.FieldTheory.IsAlgClosed.Basic

/-!
# Galois (G1) 第 62 ブロック —— **★★★★★代数閉体上の点と互いに素性**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★互いに素性を「点なしの形」に直す

第 61 ブロックの `Phi_PSq_no_common_root` は「`x` が曲線上の点の x 座標である」ことを
仮定していた。★代数閉体なら**任意の `x` に対して点が取れる**ので、その仮定が外れる:

    Φ_n(x) = 0  ⟹  ΨSq_n(x) ≠ 0

★★これにより `Φ_n − c·ΨSq_n` の**すべての根**で乗法公式の分母が生きる。

## ★★点の存在

`y² + (a₁x+a₃)y − (x³+a₂x²+a₄x+a₆) = 0` は `y` について 2 次なので、
代数閉体では根を持つ。★`Δ ≠ 0` なら方程式を満たす点は自動的に非特異
(mathlib `equation_iff_nonsingular_of_Δ_ne_zero`)。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `exists_nonsingular` | ★代数閉体では任意の `x` に点がある |
| `PSq_ne_zero_of_Phi_root` | ★★★★★`Φ_n` の根では `ΨSq_n ≠ 0` |
| `PSq_ne_zero_of_root` | ★★★★★`Φ_n − cΨSq_n` の根でも同様 |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve Polynomial WeierstrassCurve.Affine

universe u

variable {F : Type u} [Field F] [DecidableEq F] (W : WeierstrassCurve F)

/-- ★**代数閉体では任意の `x` に対して曲線上の点がある**。 -/
theorem exists_nonsingular [IsAlgClosed F] (hΔ : W.Δ ≠ 0) (x : F) :
    ∃ y : F, W.toAffine.Nonsingular x y := by
  set p : F[X] := X ^ 2 + C (W.a₁ * x + W.a₃) * X
    - C (x ^ 3 + W.a₂ * x ^ 2 + W.a₄ * x + W.a₆) with hp
  have hc2 : p.coeff 2 = 1 := by rw [hp]; simp [-map_add, -map_mul, -map_pow]
  have hdeg : p.degree ≠ 0 := by
    intro h0
    have hle := Polynomial.le_degree_of_ne_zero (p := p) (n := 2) (by rw [hc2]; exact one_ne_zero)
    rw [h0] at hle
    exact absurd hle (by decide)
  obtain ⟨y, hy⟩ := IsAlgClosed.exists_root p hdeg
  refine ⟨y, (WeierstrassCurve.Affine.equation_iff_nonsingular_of_Δ_ne_zero hΔ).mp ?_⟩
  rw [WeierstrassCurve.Affine.equation_iff]
  have hy' := hy
  rw [Polynomial.IsRoot, hp] at hy'
  simp only [eval_sub, eval_add, eval_mul, eval_pow, eval_C, eval_X] at hy'
  linear_combination hy'

/-- ★★★★★**`Φ_n` の根では `ΨSq_n ≠ 0`**(代数閉体)。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to -/
theorem PSq_ne_zero_of_Phi_root [IsAlgClosed F] (hΔ : W.Δ ≠ 0) (n : ℕ) (hn : 1 ≤ n) (x : F)
    (hPhi : (W.Φ (n : ℤ)).eval x = 0) : (W.ΨSq (n : ℤ)).eval x ≠ 0 := by
  intro hc
  obtain ⟨y, hy⟩ := exists_nonsingular W hΔ x
  exact Phi_PSq_no_common_root W hΔ hy n hn hPhi hc

/-- ★★★★★**`Φ_n − c·ΨSq_n` の根でも `ΨSq_n ≠ 0`**。 -/
theorem PSq_ne_zero_of_root [IsAlgClosed F] (hΔ : W.Δ ≠ 0) (n : ℕ) (hn : 1 ≤ n) (c x : F)
    (h : (W.Φ (n : ℤ)).eval x = c * (W.ΨSq (n : ℤ)).eval x) :
    (W.ΨSq (n : ℤ)).eval x ≠ 0 := by
  intro hc
  rw [hc, mul_zero] at h
  exact PSq_ne_zero_of_Phi_root W hΔ n hn x h hc

/-! ## ★出典の紐付け(`.src`) -/

def PSq_ne_zero_of_Phi_root.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(代数閉体上での Φ_n と ΨSq_n の互いに素性)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
