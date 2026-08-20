import ABC3.Found.GaloisRep.NonDegen

/-!
# Galois (G1) 第 51 ブロック —— **★★★★★★帰納段(x 側・y 側)**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★体の水準に切り出す

点の水準の議論(mathlib の `Point.add`)と、多項式の水準の恒等式(第 48・49 ブロック)を
繋ぐのが本ブロックである。★**体の中の等式だけ**を扱うので、点の場合分けと分離できる。

| 仮定 | 由来 |
|---|---|
| `hxm`, `hxmm` | ★帰納の仮定(x 側) |
| `hym` | ★帰納の仮定(y 側) |
| `hnew` | ★mathlib `addX_eq_addX_negY_sub` |
| `hnewy` | ★mathlib `addY_sub_negY_addY` |
| `hpsi2` | ★曲線の式(`ψ² = Ψ₂Sq(x)`) |

## ★★★x 側の鎖

    (x − x_m)·ΨSq_m = preΨ_{m+1}preΨ_{m−1}f_m           (第 44)
    ⟹ (x − x_m)²·ΨSq_m² = ΨSq_{m+1}·ΨSq_{m−1}          (第 50)
    ⟹ ΨSq_{m+1}(Φ_{m−1} − x_new ΨSq_{m−1}) = preΨ_{2m}Ψ₂Sq
    ⟹ x_new·ΨSq_{m+1} = Φ_{m+1}                        (第 48)

## ★★★y 側の鎖

    ψ_new·(x − x_m) = ψ(x_m − x_new) − ψ_m(x − x_new)   (mathlib)
    ×ΨSq_m²ΨSq_{m+1}² して (Y)(第 49)を使う
    ⟹ ψ_new·ΨSq_{m+1}² = preΨ_{2(m+1)}·ψ

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `xstep` | ★★★★★**x 側の帰納段** |
| `ystep` | ★★★★★★**y 側の帰納段** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve Polynomial

universe u

variable {F : Type u} [Field F] (W : WeierstrassCurve F)

/-- ★★★★★**x 側の帰納段**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★`hnew` は mathlib の `addX_eq_addX_negY_sub` を分母払いした形である。 -/
theorem xstep (m : ℤ) (x xm xmm xnew psim psi : F)
    (hxm : xm * (W.ΨSq m).eval x = (W.Φ m).eval x)
    (hxmm : xmm * (W.ΨSq (m - 1)).eval x = (W.Φ (m - 1)).eval x)
    (hym : psim * (W.ΨSq m).eval x ^ 2 = (W.preΨ (2 * m)).eval x * psi)
    (hpsi2 : psi ^ 2 = W.Ψ₂Sq.eval x)
    (hnew : (x - xm) ^ 2 * (xmm - xnew) = psim * psi)
    (hAmm : (W.ΨSq (m - 1)).eval x ≠ 0) :
    xnew * (W.ΨSq (m + 1)).eval x = (W.Φ (m + 1)).eval x := by
  have hd : (x - xm) * (W.ΨSq m).eval x
      = (W.preΨ (m + 1)).eval x * (W.preΨ (m - 1)).eval x
        * (if Even m then (1 : F) else W.Ψ₂Sq.eval x) := by
    rw [sub_mul, hxm]; exact x_mul_PSq_sub_Phi W m x
  have hD2 : (x - xm) ^ 2 * (W.ΨSq m).eval x ^ 2
      = (W.ΨSq (m + 1)).eval x * (W.ΨSq (m - 1)).eval x := by
    rw [PSq_mul_PSq]; linear_combination ((x - xm) * (W.ΨSq m).eval x
      + (W.preΨ (m + 1)).eval x * (W.preΨ (m - 1)).eval x
        * (if Even m then (1 : F) else W.Ψ₂Sq.eval x)) * hd
  have hcross := congrArg (Polynomial.eval x) (Phi_succ_mul_PSq_pred W m)
  simp only [eval_mul, eval_sub] at hcross
  have step2 : (W.ΨSq (m + 1)).eval x * ((W.Φ (m - 1)).eval x - xnew * (W.ΨSq (m - 1)).eval x)
      * (W.ΨSq (m - 1)).eval x
      = (W.preΨ (2 * m)).eval x * W.Ψ₂Sq.eval x * (W.ΨSq (m - 1)).eval x := by
    linear_combination ((W.ΨSq m).eval x ^ 2 * (W.ΨSq (m - 1)).eval x) * hnew
      - ((xmm - xnew) * (W.ΨSq (m - 1)).eval x) * hD2
      + (psi * (W.ΨSq (m - 1)).eval x) * hym
      + ((W.preΨ (2 * m)).eval x * (W.ΨSq (m - 1)).eval x) * hpsi2
      - ((W.ΨSq (m + 1)).eval x * (W.ΨSq (m - 1)).eval x) * hxmm
  have step3 : (W.ΨSq (m + 1)).eval x * ((W.Φ (m - 1)).eval x - xnew * (W.ΨSq (m - 1)).eval x)
      = (W.preΨ (2 * m)).eval x * W.Ψ₂Sq.eval x := mul_right_cancel₀ hAmm step2
  refine mul_right_cancel₀ hAmm ?_
  linear_combination -step3 - hcross

/-- ★★★★★★**y 側の帰納段**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★`hnewy` は mathlib の `addY_sub_negY_addY` を分母払いした形である。
★★第 49 ブロックの (Y) がここで効く。 -/
theorem ystep (m : ℤ) (x xm xnew psi psim psinew : F)
    (hxm : xm * (W.ΨSq m).eval x = (W.Φ m).eval x)
    (hxnew : xnew * (W.ΨSq (m + 1)).eval x = (W.Φ (m + 1)).eval x)
    (hym : psim * (W.ΨSq m).eval x ^ 2 = (W.preΨ (2 * m)).eval x * psi)
    (hnewy : psinew * (x - xm) = psi * (xm - xnew) - psim * (x - xnew))
    (hAm : (W.ΨSq m).eval x ≠ 0) (hxne : x - xm ≠ 0) :
    psinew * (W.ΨSq (m + 1)).eval x ^ 2 = (W.preΨ (2 * (m + 1))).eval x * psi := by
  have hys := congrArg (Polynomial.eval x) (y_side W m)
  simp only [eval_mul, eval_sub, apply_ite (Polynomial.eval x), eval_one] at hys
  have hQ : (W.preΨ (m + 2)).eval x * (W.preΨ m).eval x
        * (if Even (m + 1) then (1 : F) else W.Ψ₂Sq.eval x)
      = (x - xnew) * (W.ΨSq (m + 1)).eval x := by
    have h := x_mul_PSq_sub_Phi W (m + 1) x
    rw [show m + 1 + 1 = m + 2 by ring, show m + 1 - 1 = m by ring] at h
    rw [← h, sub_mul, hxnew]
  have hD : (W.preΨ (m + 1)).eval x * (W.preΨ (m - 1)).eval x
        * (if Even m then (1 : F) else W.Ψ₂Sq.eval x)
      = (x - xm) * (W.ΨSq m).eval x := by
    rw [← x_mul_PSq_sub_Phi W m x, sub_mul, hxm]
  rw [hQ, hD] at hys
  refine mul_right_cancel₀ (mul_ne_zero hxne (pow_ne_zero 2 hAm)) ?_
  linear_combination ((W.ΨSq (m + 1)).eval x ^ 2 * (W.ΨSq m).eval x ^ 2) * hnewy
    - ((W.ΨSq (m + 1)).eval x ^ 2 * (x - xnew)) * hym
    + psi * hys
    + (psi * (W.ΨSq (m + 1)).eval x ^ 2 * (W.ΨSq m).eval x) * hxm
    - (psi * (W.ΨSq m).eval x ^ 2 * (W.ΨSq (m + 1)).eval x) * hxnew

/-! ## ★出典の紐付け(`.src`) -/

def xstep.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(乗法公式の帰納段——体の水準)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
