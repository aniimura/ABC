import ABC3.Found.GaloisRep.XDiff

/-!
# Galois (G1) 第 45 ブロック —— **★★★★★★`y` 側の基底 `ψ(2P)·d³ = preΨ₄`**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★同時帰納の `y` 側の基底

§9-386 で導いた形

    ψ(nP) · ΨSq_n(x)² = preΨ_{2n}(x) · ψ(P)

の `n = 2` である。★`ΨSq_2 = Ψ₂Sq = d²`、`ψ(P) = d` なので

    ψ(2P) · d⁴ = preΨ₄ · d    ⟺    ψ(2P) · d³ = preΨ₄

## ★★通し方——第 30 と同じ型

★**除算を先に潰す**。`ℓ = N/d` を含む部分だけ `field_simp` で払い(`hkey`)、
曲線の方程式は `linear_combination` で入れる。

★★係数(24 項)は Python で剰余計算して得た——**一発で通った**。

## ★★★これで同時帰納の部品が全部揃った

| 部品 | ブロック |
|---|---|
| `M_x` の基底(倍化公式) | ★第 30 |
| **`M_y` の基底** | ★★**本ブロック** |
| 分母 `x − x(nP)` | ★第 44 |
| `nP = O` / `(n+1)P = O` の場合 | ★第 41–43 |
| `ω_n` | ★★第 27–28 |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve Polynomial WeierstrassCurve.Affine

universe u

variable {F : Type u} [Field F] [DecidableEq F] (W : WeierstrassCurve F)

theorem psi_double {x y : F} (h : W.toAffine.Equation x y)
    (hy : y ≠ W.toAffine.negY x y) :
    (W.toAffine.addY x x y (W.toAffine.slope x x y y)
      - W.toAffine.negY (W.toAffine.addX x x (W.toAffine.slope x x y y))
          (W.toAffine.addY x x y (W.toAffine.slope x x y y)))
      * (y - W.toAffine.negY x y) ^ 3 = W.preΨ₄.eval x := by
  have hd : y - W.toAffine.negY x y ≠ 0 := sub_ne_zero.mpr hy
  set d := y - W.toAffine.negY x y with hdd
  set N := 3 * x ^ 2 + 2 * W.a₂ * x + W.a₄ - W.a₁ * y with hN
  have hkey : (W.toAffine.addY x x y (N / d)
      - W.toAffine.negY (W.toAffine.addX x x (N / d)) (W.toAffine.addY x x y (N / d))) * d ^ 3
      = -2 * N * ((N ^ 2 + W.a₁ * N * d - (W.a₂ + 2 * x) * d ^ 2) - x * d ^ 2)
        - (2 * y + W.a₃) * d ^ 3
        - W.a₁ * (N ^ 2 + W.a₁ * N * d - (W.a₂ + 2 * x) * d ^ 2) * d := by
    rw [WeierstrassCurve.Affine.addY, WeierstrassCurve.Affine.negAddY,
      WeierstrassCurve.Affine.addX, WeierstrassCurve.Affine.negY,
      WeierstrassCurve.Affine.negY]
    field_simp
    ring
  rw [W.toAffine.slope_of_Y_ne rfl hy, ← hN, ← hdd, hkey, WeierstrassCurve.preΨ₄]
  simp only [eval_add, eval_mul, eval_pow, eval_C, eval_X, eval_ofNat]
  rw [WeierstrassCurve.b₂, WeierstrassCurve.b₄, WeierstrassCurve.b₆, WeierstrassCurve.b₈, hN, hdd,
    WeierstrassCurve.Affine.negY]
  rw [WeierstrassCurve.Affine.Equation] at h
  simp only [WeierstrassCurve.Affine.polynomial, evalEval, eval_add, eval_sub, eval_mul,
    eval_pow, eval_C, eval_X] at h
  linear_combination (W.a₁ ^ 4 * x + W.a₁ ^ 3 * W.a₃ + 8 * W.a₁ ^ 2 * W.a₂ * x
    + 2 * W.a₁ ^ 2 * W.a₄ + 10 * W.a₁ ^ 2 * x ^ 2 + 4 * W.a₁ * W.a₂ * W.a₃
    - 4 * W.a₁ * W.a₃ * x + 16 * W.a₂ ^ 2 * x + 8 * W.a₂ * W.a₄ + 56 * W.a₂ * x ^ 2
    - 8 * W.a₃ ^ 2 + 8 * W.a₄ * x - 16 * W.a₆ + 56 * x ^ 3 - 16 * y ^ 2
    - 16 * y * (W.a₁ * x + W.a₃)) * h


/-! ## ★出典の紐付け(`.src`) -/

def psi_double.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(乗法公式の同時帰納——y 側の基底)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
