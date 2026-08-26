import ABC3.Found.GaloisRep.TateLambda

/-!
# Galois (G6) 第 264 ブロック —— **★★★★★★★★葉 (e) の同定の段**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★★★母数座標が一致すれば点が一致する

葉 (e)(全射性)は 2 段に分かれる:

| 段 | 内容 | 状態 |
|---|---|---|
| 存在 | `Λ(a) = c` を満たす `a` を縮小写像で作る | 次のブロック |
| **同定** | **その `a` の点が目標の点に等しい** | ★本ブロック |

本ブロックは**同定の段**である:

> `(x,y)` が曲線上、`x` と `x+y` が単元、`Λ(a) = y·(x+y)⁻¹` なら
> `X(a,w) = x` かつ `Y(a,w) = y`

## ★★★★★★段取り

1. `Λ(a) = Y·(X+Y)⁻¹` と `c = y·(x+y)⁻¹` の等式を、両辺に `(X+Y)(x+y)` を掛けて
   **逆元を消す**:`Y(x+y) = y(X+Y)`。
2. 展開すると **`Y x = y X`**——第 263 の `curve_point_unique_of_rel` の仮定の形。
3. `(X,Y)` は曲線上(葉 (b)、第 240 の `tate_equation`)、`X` は単元
   (`isUnit_tateXpair`——主要項 `a/(1−a)²` が単元)。
4. 第 263 の一意性で `X = x`、`Y = y`。

★★逆元を「掛けて消す」だけで済むのは、`X + Y` と `x + y` が**両方とも単元**だから
である(第 262 の `isUnit_tateXpair_add_tateYpair`)。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `isUnit_tateXpair` | ★★★★★`X` は単元 |
| `tate_point_eq_of_lambda` | ★★★★★★★★**母数座標が一致すれば点が一致する** |
-/

namespace ABC3.Found.GaloisRep

open Complex Real

variable {R : Type} [CommRing R] {I : Ideal R}

/-- ★★★★★**`X` は単元**——`a`・`1 − a` が単元、`w ∈ I` のとき。 -/
theorem isUnit_tateXpair [IsAdicComplete I R] (a w q : R) (hq : q ∈ I) (hw : w ∈ I)
    (ha : IsUnit a) (hu : IsUnit (1 - a)) : IsUnit (tateXpair a w q hq) := by
  have hf : IsUnit (tateXterm a) := by
    have hmul : IsUnit ((1 - a) ^ 2 * tateXterm a) := by
      rw [mul_tateXterm' hu]; exact ha
    exact isUnit_of_mul_isUnit_right hmul
  have hd := tateXpair_sub_mem a w q hq hw
  have hsplit : tateXpair a w q hq = tateXterm a + (tateXpair a w q hq - tateXterm a) := by ring
  rw [hsplit]
  exact isUnit_add_mem hf hd

set_option maxHeartbeats 800000 in
/-- ★★★★★★★★**母数座標が一致すれば点が一致する**——葉 (e) の同定の段。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem tate_point_eq_of_lambda [IsAdicComplete I R] (x y q : R) (hq : q ∈ I)
    (hx : IsUnit x) (hxy : IsUnit (x + y))
    (he : y ^ 2 + x * y = x ^ 3 + (tateCurveAt q hq).a₄ * x + (tateCurveAt q hq).a₆)
    (a w : R) (hw : w ∈ I) (haw : a * w = q) (ha : IsUnit a) (hu : IsUnit (1 - a))
    (hΛ : tateLambda a w q hq = y * Ring.inverse (x + y)) :
    tateXpair a w q hq = x ∧ tateYpair a w q hq = y := by
  have hS : IsUnit (tateXpair a w q hq + tateYpair a w q hq) :=
    isUnit_tateXpair_add_tateYpair a w q hq hw ha hu
  have hSinv : (tateXpair a w q hq + tateYpair a w q hq)
      * Ring.inverse (tateXpair a w q hq + tateYpair a w q hq) = 1 :=
    Ring.mul_inverse_cancel _ hS
  have hsinv : (x + y) * Ring.inverse (x + y) = 1 := Ring.mul_inverse_cancel _ hxy
  have hcross : tateYpair a w q hq * (x + y)
      = y * (tateXpair a w q hq + tateYpair a w q hq) := by
    have h1 := congrArg (fun z => z * ((tateXpair a w q hq + tateYpair a w q hq) * (x + y))) hΛ
    simp only [tateLambda] at h1
    calc tateYpair a w q hq * (x + y)
        = tateYpair a w q hq
          * Ring.inverse (tateXpair a w q hq + tateYpair a w q hq)
          * ((tateXpair a w q hq + tateYpair a w q hq) * (x + y)) := by
          rw [show tateYpair a w q hq
              * Ring.inverse (tateXpair a w q hq + tateYpair a w q hq)
              * ((tateXpair a w q hq + tateYpair a w q hq) * (x + y))
            = tateYpair a w q hq * ((tateXpair a w q hq + tateYpair a w q hq)
              * Ring.inverse (tateXpair a w q hq + tateYpair a w q hq)) * (x + y) by ring,
            hSinv, mul_one]
      _ = y * Ring.inverse (x + y) * ((tateXpair a w q hq + tateYpair a w q hq) * (x + y)) := h1
      _ = y * (tateXpair a w q hq + tateYpair a w q hq) := by
          rw [show y * Ring.inverse (x + y)
              * ((tateXpair a w q hq + tateYpair a w q hq) * (x + y))
            = y * ((x + y) * Ring.inverse (x + y))
              * (tateXpair a w q hq + tateYpair a w q hq) by ring, hsinv, mul_one]
  have hrel : tateYpair a w q hq * x = y * tateXpair a w q hq := by linear_combination hcross
  have heq := tate_equation a w q hq haw hu (isUnit_one_sub hw)
  exact curve_point_unique_of_rel (I := I) (tateCurveAt_a4_mem q hq) (tateCurveAt_a6_mem q hq)
    (by linear_combination heq) he (isUnit_tateXpair a w q hq hw ha hu) hx hrel

/-! ## ★出典の紐付け(`.src`) -/

def tate_point_eq_of_lambda.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——母数座標が一致すれば点が一致する)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
