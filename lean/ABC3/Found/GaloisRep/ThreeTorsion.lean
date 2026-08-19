import ABC3.Found.GaloisRep.FiniteTorsion

/-!
# Galois (G1) 第 36 ブロック —— **★★★★★★3-捩れと `E[3]` の有限性**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★倍化公式から `n = 3` が出た

    3 • P = 0  ⟹  Ψ₃(x_P) = 0

★★機構は短い:

| 段 | 内容 |
|---|---|
| 1 | `3•P = 0` かつ `2•P = 0` なら `P = 0`(矛盾)——よって `y ≠ negY x y` |
| 2 | `3•P = 0` ⟹ `2•P = −P` ⟹ **`x(2P) = x`** |
| 3 | 第 30 の倍化公式: `x(2P)·Ψ₂Sq(x) = Φ₂(x)` |
| 4 | ★★★**`Φ₂ = X·Ψ₂Sq − Ψ₃`** ⟹ `Ψ₃(x) = 0` |

★★★★段 4 の恒等式が鍵である——`Φ` と `Ψ` の定義から `ring` で出る。

## ★★これで `E[3]` も有限になる

第 35 の `finite_torsion_of_formula` に代入するだけ(`ΨSq 3 = Ψ₃²`)。

## ★★★一般の `n` への含意

★`n = 2`(第 31)と `n = 3`(本ブロック)は**別々の機構**で出た。
★★一般の `n` には乗法公式の帰納が要る——ただし
**倍化公式(第 30)がその基底**であることが本ブロックで確かめられた。
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve Polynomial WeierstrassCurve.Affine

universe u

variable {F : Type u} [Field F] [DecidableEq F] (W : WeierstrassCurve F)

/-- `Φ₂ = X·Ψ₂Sq − Ψ₃`。 -/
theorem Phi_two_eq : W.Φ 2 = X * W.Ψ₂Sq - W.Ψ₃ := by
  rw [WeierstrassCurve.Φ_two, WeierstrassCurve.Ψ₂Sq, WeierstrassCurve.Ψ₃]
  simp only [map_mul, map_ofNat]
  ring

theorem three_torsion_root (h2 : (2 : F) ≠ 0) {x y : F} (h : W.toAffine.Nonsingular x y)
    (h3 : (3 : ℕ) • (Point.some x y h) = 0) : W.Ψ₃.eval x = 0 := by
  have hP0 : (Point.some x y h) ≠ 0 := Point.some_ne_zero h
  have hy : y ≠ W.toAffine.negY x y := by
    intro hyy
    have h2P : (2 : ℕ) • (Point.some x y h) = 0 :=
      (two_smul_eq_zero_iff W h).2 ((psi2Sq_eval_eq_zero_iff W h.left).2 hyy)
    refine hP0 ?_
    have : (3 : ℕ) • (Point.some x y h) - (2 : ℕ) • (Point.some x y h)
        = Point.some x y h := by
      rw [show (3 : ℕ) = 2 + 1 from rfl, add_smul, one_smul]
      abel
    rw [h3, h2P, sub_zero] at this
    exact this.symm
  have hneg : (2 : ℕ) • (Point.some x y h) = -(Point.some x y h) := by
    rw [eq_neg_iff_add_eq_zero, ← succ_nsmul]
    exact h3
  rw [two_nsmul, Point.add_self_of_Y_ne hy, Point.neg_some h, Point.some.injEq] at hneg
  have hx : W.toAffine.addX x x (W.toAffine.slope x x y y) = x := hneg.1
  have hdbl := doubling_x W h.left hy
  rw [hx, Phi_two_eq W] at hdbl
  simp only [eval_sub, eval_mul, eval_X] at hdbl
  linear_combination hdbl


/-- ★★★★**`E[3]` は有限**。 -/
theorem finite_three_torsion (h2 : (2 : F) ≠ 0) (h3 : (3 : F) ≠ 0) :
    {P : W.toAffine.Point | (3 : ℕ) • P = 0}.Finite := by
  refine finite_torsion_of_formula W (n := (3 : ℤ)) (by exact_mod_cast h3) ?_
  intro x y h hP
  rw [Polynomial.IsRoot, WeierstrassCurve.ΨSq_three, eval_pow,
    three_torsion_root W h2 h hP]
  ring

/-! ## ★出典の紐付け(`.src`) -/

def three_torsion_root.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(3-捩れと E[3] の有限性)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
