import ABC3.Found.GaloisRep.TorsionMul

/-!
# Galois (G1) 第 39 ブロック —— **★★★★★★`4P = O ⟹ ΨSq₄(x) = 0`**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★倍化公式の**合成**で `n = 4` が出た

    4P = O  ⟺  2(2P) = O  ⟹  Ψ₂Sq(x(2P)) = 0

★`x(2P) = Φ₂/Ψ₂Sq`(第 30)を代入し、分母 `Ψ₂Sq³` を払うと

    4Φ₂³ + b₂Φ₂²Ψ₂Sq + 2b₄Φ₂Ψ₂Sq² + b₆Ψ₂Sq³ = preΨ₄²

★★★**ぴったり `preΨ₄²` になる**(計算機で確認、b-不変量の関係
`4b₈ = b₂b₆ − b₄²` を使う)。したがって `preΨ₄(x) = 0`、
`ΨSq₄ = preΨ₄²·Ψ₂Sq` なので `ΨSq₄(x) = 0`。

## ★★★これが `2^k` へ伸びる型である

★倍化公式を**合成する**だけで `n → 2n` が出る。
★★`E[2^k]` の根はこの繰り返しで得られる。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `psi2Sq_comp_phi2` | ★★★合成した分母払いが `preΨ₄²` になること |
| `four_torsion_root` | ★★★★**`4P = O ⟹ ΨSq₄(x) = 0`** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve Polynomial WeierstrassCurve.Affine

universe u

variable {F : Type u} [Field F] [DecidableEq F] (W : WeierstrassCurve F)

/-- ★`Ψ₂Sq` を `Φ₂/Ψ₂Sq` で合成し分母を払うと `preΨ₄²` になる。 -/
theorem psi2Sq_comp_phi2 (x : F) :
    4 * ((W.Φ 2).eval x) ^ 3 + W.b₂ * ((W.Φ 2).eval x) ^ 2 * (W.Ψ₂Sq.eval x)
        + 2 * W.b₄ * ((W.Φ 2).eval x) * (W.Ψ₂Sq.eval x) ^ 2
        + W.b₆ * (W.Ψ₂Sq.eval x) ^ 3
      = (W.preΨ₄.eval x) ^ 2 := by
  rw [WeierstrassCurve.Φ_two, WeierstrassCurve.Ψ₂Sq, WeierstrassCurve.preΨ₄]
  simp only [eval_add, eval_sub, eval_mul, eval_pow, eval_C, eval_X, eval_ofNat]
  linear_combination (W.b₂ ^ 2 * x ^ 6 + 6 * W.b₂ * W.b₄ * x ^ 5 + 7 * W.b₂ * W.b₆ * x ^ 4
    + 8 * W.b₂ * x ^ 7 + 8 * W.b₄ ^ 2 * x ^ 4 + 16 * W.b₄ * W.b₆ * x ^ 3 + 28 * W.b₄ * x ^ 6
    + 7 * W.b₆ ^ 2 * x ^ 2 + 38 * W.b₆ * x ^ 5 + W.b₈ ^ 2
    + 2 * W.b₈ * x * (2 * W.b₂ * x ^ 2 + 4 * W.b₄ * x + 3 * W.b₆ + 11 * x ^ 3)
    + 13 * x ^ 8) * (-W.b_relation)

theorem four_torsion_root {x y : F} (h : W.toAffine.Nonsingular x y)
    (h4 : (4 : ℕ) • (Point.some x y h) = 0) : (W.ΨSq 4).IsRoot x := by
  rw [Polynomial.IsRoot, WeierstrassCurve.ΨSq_four, eval_mul, eval_pow]
  by_cases hy : y = W.toAffine.negY x y
  · rw [(psi2Sq_eval_eq_zero_iff W h.left).2 hy]; ring
  · have hQ : (2 : ℕ) • (Point.some x y h)
        = Point.some _ _ (nonsingular_add h h fun hxy => hy hxy.right) := by
      rw [two_nsmul]; exact Point.add_self_of_Y_ne hy
    have h2Q : (2 : ℕ) • ((2 : ℕ) • (Point.some x y h)) = 0 := by
      rw [smul_smul]; exact h4
    rw [hQ] at h2Q
    have hroot := (two_smul_eq_zero_iff W (nonsingular_add h h fun hxy => hy hxy.right)).1 h2Q
    have hdbl := doubling_x W h.left hy
    rw [WeierstrassCurve.Ψ₂Sq] at hroot
    simp only [eval_add, eval_mul, eval_pow, eval_C, eval_X, eval_ofNat] at hroot
    have key : (W.preΨ₄.eval x) ^ 2 = 0 := by
      rw [← psi2Sq_comp_phi2 W x, ← hdbl]
      linear_combination (W.Ψ₂Sq.eval x) ^ 3 * hroot
    rw [pow_eq_zero_iff two_ne_zero] at key
    rw [key]
    ring


/-! ## ★出典の紐付け(`.src`) -/

def four_torsion_root.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(4P = O ならば ΨSq₄(x) = 0)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
