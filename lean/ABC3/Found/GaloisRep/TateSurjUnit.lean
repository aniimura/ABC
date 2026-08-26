import ABC3.Found.GaloisRep.TateContract

/-!
# Galois (G6) 第 266 ブロック —— **★★★★★★★★★葉 (e)(単元の領域)**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★★★★到達点

> `(x, y)` が `E_q` の上の点で **`x` が単元**なら、単元 `a` と `w ∈ I` が在って
> `a·w = q`、`X(a,w) = x`、`Y(a,w) = y`(`tate_surjective`)

★これが葉 (e) の**単元の領域**での全射性である。

## ★★★★★★不動点で `a` を作る

第 263 の `Λ(a) = a + κ(a)`(`κ(a) ∈ I`)と第 265 の
`κ(a) − κ(b) ∈ I^{j+1}` を使い、`a = c + b`(`b ∈ I`)と置いて

    b = −κ(c + b)

の不動点を第 102 の `exists_fixedPoint_of_contraction` で作る(`exists_tate_lambda_eq`)。
★★`c := y·(x+y)⁻¹` が目標の母数座標である。

## ★★★仮定は `x` が単元だけでよい

`y(y+x) = x³ + a₄x + a₆` で、右辺は `x³`(単元)+ `I` の元なので**単元**。よって

    `y` も `x + y` も単元

が**曲線の式から出る**(`isUnit_of_mul_isUnit_left`/`_right`)。★`c` と `1 − c` の
単元性もそこから出る(`1 − c = x·(x+y)⁻¹`)。

## ★残る領域

`x` が単元でない点は 2 種類ある。

| 領域 | `u` | `x = X(u)` |
|---|---|---|
| ★済 | `u` が単元、`1 − u` が単元 | 単元 |
| 未 | `v(u) > 0`(環帯の内部) | `x ∈ I` |
| 未 | `u ≡ 1`(原点の近く) | `x ∉ R`(`K` の中で極) |

★第 2 の領域は `X ≈ u + q/u` なので別の縮小写像、第 3 の領域は形式群である。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `tateLambda_sub_eq` | ★★★★★`Λ(a) − a = ε(a)·inv(S(a))` |
| `tateS_diff_mem` | ★★★★★`S = X + Y` の差 |
| `tateKappa_diff_mem` | ★★★★★★★`κ` の差は 1 つ位が上がる |
| `exists_tate_lambda_eq` | ★★★★★★★★**`Λ(a) = c` の解の存在** |
| `tate_surjective` | ★★★★★★★★★**葉 (e)(単元の領域)** |
-/

namespace ABC3.Found.GaloisRep

open Complex Real

variable {R : Type} [CommRing R] {I : Ideal R}

/-! ## ★★★★★誤差項と `S` の差 -/

/-- ★★★★★`Λ(a) − a = ε(a)·inv(S(a))`。 -/
theorem tateLambda_sub_eq [IsAdicComplete I R] (a w q : R) (hq : q ∈ I) (hw : w ∈ I)
    (ha : IsUnit a) (hu : IsUnit (1 - a)) :
    tateLambda a w q hq - a
      = (tateYpair a w q hq - a * (tateXpair a w q hq + tateYpair a w q hq))
        * Ring.inverse (tateXpair a w q hq + tateYpair a w q hq) := by
  have hS : IsUnit (tateXpair a w q hq + tateYpair a w q hq) :=
    isUnit_tateXpair_add_tateYpair a w q hq hw ha hu
  rw [tateLambda]
  have h2 : a = a * ((tateXpair a w q hq + tateYpair a w q hq)
      * Ring.inverse (tateXpair a w q hq + tateYpair a w q hq)) := by
    rw [Ring.mul_inverse_cancel _ hS, mul_one]
  calc tateYpair a w q hq * Ring.inverse (tateXpair a w q hq + tateYpair a w q hq) - a
      = tateYpair a w q hq * Ring.inverse (tateXpair a w q hq + tateYpair a w q hq)
        - a * ((tateXpair a w q hq + tateYpair a w q hq)
          * Ring.inverse (tateXpair a w q hq + tateYpair a w q hq)) := by rw [← h2]
    _ = _ := by ring

set_option maxHeartbeats 1200000 in
/-- ★★★★★`S = X + Y` の差。 -/
theorem tateS_diff_mem [IsAdicComplete I R] {k : ℕ} (a b q : R) (hq : q ∈ I)
    (hau : IsUnit a) (hbu : IsUnit b) (ha1 : IsUnit (1 - a)) (hb1 : IsUnit (1 - b))
    (hab : a - b ∈ I ^ k) :
    (tateXpair a (wOf q a) q hq + tateYpair a (wOf q a) q hq)
      - (tateXpair b (wOf q b) q hq + tateYpair b (wOf q b) q hq) ∈ I ^ k := by
  have hwa : wOf q a ∈ I := wOf_mem a hq
  have hwb : wOf q b ∈ I := wOf_mem b hq
  have hwab : wOf q a - wOf q b ∈ I ^ (k + 1) := wOf_diff_mem hq hau hbu hab
  have hq1 : q ∈ I ^ 1 := by rwa [pow_one]
  have hup : ∀ (n : ℕ), ∀ x : R, x ∈ I ^ (k + n) → x ∈ I ^ k := fun n x hx =>
    Ideal.pow_le_pow_right (by omega) hx
  have dfa : tateXterm a - tateXterm b ∈ I ^ k := tateXterm_diff_mem_gen ha1 hb1 hab
  have dga : tateYterm a - tateYterm b ∈ I ^ k := tateYterm_diff_mem_gen ha1 hb1 hab
  have dTfa : tateXtail a q hq - tateXtail b q hq ∈ I ^ k :=
    hup 1 _ (tateXtail_diff_mem (m := 1) a b q hq hq1 hab)
  have dTga : tateYtail a q hq - tateYtail b q hq ∈ I ^ k :=
    hup 1 _ (tateYtail_diff_mem (m := 1) a b q hq hq1 hab)
  have dgw : tateYterm (wOf q a) - tateYterm (wOf q b) ∈ I ^ k :=
    hup 1 _ (tateYterm_diff_mem_gen (isUnit_one_sub hwa) (isUnit_one_sub hwb) hwab)
  have dTgw : tateYtail (wOf q a) q hq - tateYtail (wOf q b) q hq ∈ I ^ k :=
    hup 2 _ (tateYtail_diff_mem (m := 1) (wOf q a) (wOf q b) q hq hq1 hwab)
  have hkey : (tateXpair a (wOf q a) q hq + tateYpair a (wOf q a) q hq)
      - (tateXpair b (wOf q b) q hq + tateYpair b (wOf q b) q hq)
      = (tateXterm a - tateXterm b) + (tateYterm a - tateYterm b)
        + (tateXtail a q hq - tateXtail b q hq) + (tateYtail a q hq - tateYtail b q hq)
        - (tateYterm (wOf q a) - tateYterm (wOf q b))
        - (tateYtail (wOf q a) q hq - tateYtail (wOf q b) q hq) := by
    rw [tateXpair, tateYpair, tateXpair, tateYpair]
    ring
  rw [hkey]
  exact Ideal.sub_mem _ (Ideal.sub_mem _ (Ideal.add_mem _ (Ideal.add_mem _
    (Ideal.add_mem _ dfa dga) dTfa) dTga) dgw) dTgw

set_option maxHeartbeats 1200000 in
/-- ★★★★★★★**`κ = Λ − id` の差は 1 つ位が上がる**。 -/
theorem tateKappa_diff_mem [IsAdicComplete I R] {k : ℕ} (a b q : R) (hq : q ∈ I)
    (hau : IsUnit a) (hbu : IsUnit b) (ha1 : IsUnit (1 - a)) (hb1 : IsUnit (1 - b))
    (hab : a - b ∈ I ^ k) :
    (tateLambda a (wOf q a) q hq - a) - (tateLambda b (wOf q b) q hq - b) ∈ I ^ (k + 1) := by
  have hwa : wOf q a ∈ I := wOf_mem a hq
  have hwb : wOf q b ∈ I := wOf_mem b hq
  have hSa : IsUnit (tateXpair a (wOf q a) q hq + tateYpair a (wOf q a) q hq) :=
    isUnit_tateXpair_add_tateYpair a (wOf q a) q hq hwa hau ha1
  have hSb : IsUnit (tateXpair b (wOf q b) q hq + tateYpair b (wOf q b) q hq) :=
    isUnit_tateXpair_add_tateYpair b (wOf q b) q hq hwb hbu hb1
  have hdS : (tateXpair a (wOf q a) q hq + tateYpair a (wOf q a) q hq)
      - (tateXpair b (wOf q b) q hq + tateYpair b (wOf q b) q hq) ∈ I ^ k :=
    tateS_diff_mem a b q hq hau hbu ha1 hb1 hab
  have hdinv : Ring.inverse (tateXpair a (wOf q a) q hq + tateYpair a (wOf q a) q hq)
      - Ring.inverse (tateXpair b (wOf q b) q hq + tateYpair b (wOf q b) q hq) ∈ I ^ k :=
    inverse_diff_mem hSa hSb hdS
  have hdε := tateEps_diff_mem a b q hq hau hbu ha1 hb1 hab
  have hεb : tateYpair b (wOf q b) q hq
      - b * (tateXpair b (wOf q b) q hq + tateYpair b (wOf q b) q hq) ∈ I :=
    tateYpair_sub_mul_mem b (wOf q b) q hq hwb hb1
  rw [tateLambda_sub_eq a (wOf q a) q hq hwa hau ha1,
    tateLambda_sub_eq b (wOf q b) q hq hwb hbu hb1]
  have hkey : (tateYpair a (wOf q a) q hq
        - a * (tateXpair a (wOf q a) q hq + tateYpair a (wOf q a) q hq))
      * Ring.inverse (tateXpair a (wOf q a) q hq + tateYpair a (wOf q a) q hq)
    - (tateYpair b (wOf q b) q hq
        - b * (tateXpair b (wOf q b) q hq + tateYpair b (wOf q b) q hq))
      * Ring.inverse (tateXpair b (wOf q b) q hq + tateYpair b (wOf q b) q hq)
    = ((tateYpair a (wOf q a) q hq
          - a * (tateXpair a (wOf q a) q hq + tateYpair a (wOf q a) q hq))
        - (tateYpair b (wOf q b) q hq
          - b * (tateXpair b (wOf q b) q hq + tateYpair b (wOf q b) q hq)))
      * Ring.inverse (tateXpair a (wOf q a) q hq + tateYpair a (wOf q a) q hq)
      + (Ring.inverse (tateXpair a (wOf q a) q hq + tateYpair a (wOf q a) q hq)
          - Ring.inverse (tateXpair b (wOf q b) q hq + tateYpair b (wOf q b) q hq))
        * (tateYpair b (wOf q b) q hq
          - b * (tateXpair b (wOf q b) q hq + tateYpair b (wOf q b) q hq)) := by
    ring
  rw [hkey]
  refine Ideal.add_mem _ (Ideal.mul_mem_right _ _ hdε) ?_
  rw [pow_succ]
  exact Ideal.mul_mem_mul hdinv hεb

/-! ## ★★★★★★★★不動点 -/

set_option maxHeartbeats 1200000 in
/-- ★★★★★★★★**`Λ(a) = c` の解の存在**——第 102 の縮小写像定理。 -/
theorem exists_tate_lambda_eq [IsAdicComplete I R] (q c : R) (hq : q ∈ I)
    (hc : IsUnit c) (hc1 : IsUnit (1 - c)) :
    ∃ a : R, IsUnit a ∧ IsUnit (1 - a) ∧ a - c ∈ I ∧ tateLambda a (wOf q a) q hq = c := by
  have hunit : ∀ b : R, b ∈ I → IsUnit (c + b) := fun b hb => isUnit_add_mem hc hb
  have hunit1 : ∀ b : R, b ∈ I → IsUnit (1 - (c + b)) := fun b hb => by
    have h : 1 - (c + b) = (1 - c) + (-b) := by ring
    rw [h]
    exact isUnit_add_mem hc1 (neg_mem hb)
  set F : R → R := fun b => -(tateLambda (c + b) (wOf q (c + b)) q hq - (c + b)) with hFdef
  have hFI : ∀ b ∈ I, F b ∈ I := by
    intro b hb
    simp only [hFdef]
    exact neg_mem (tateLambda_sub_mem (c + b) (wOf q (c + b)) q hq (wOf_mem _ hq)
      (hunit b hb) (hunit1 b hb))
  have hcontr : ∀ x ∈ I, ∀ y ∈ I, ∀ k : ℕ, x - y ∈ I ^ k → F x - F y ∈ I ^ (k + 1) := by
    intro x hx y hy k hxy
    have hxy' : (c + x) - (c + y) ∈ I ^ k := by
      have h : (c + x) - (c + y) = x - y := by ring
      rw [h]; exact hxy
    have h := neg_mem (tateKappa_diff_mem (c + x) (c + y) q hq (hunit x hx) (hunit y hy)
      (hunit1 x hx) (hunit1 y hy) hxy')
    simp only [hFdef]
    convert h using 1
    ring
  obtain ⟨b, hb, hfix⟩ := exists_fixedPoint_of_contraction F hFI (hFI 0 (Ideal.zero_mem I)) hcontr
  refine ⟨c + b, hunit b hb, hunit1 b hb, ?_, ?_⟩
  · have h : (c + b) - c = b := by ring
    rw [h]; exact hb
  · simp only [hFdef] at hfix
    linear_combination -hfix

/-! ## ★★★★★★★★★葉 (e)(単元の領域) -/

set_option maxHeartbeats 1200000 in
/-- ★★★★★★★★★**葉 (e) —— Tate 座標は全射**(`x` が単元の領域)。

★仮定は `x` が単元だけでよい:`y(y+x) = x³ + a₄x + a₆` の右辺は `x³`(単元)+ `I` の元
なので単元、よって `y` も `x + y` も**曲線の式から**単元になる。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem tate_surjective [IsAdicComplete I R] (x y q : R) (hq : q ∈ I) (hx : IsUnit x)
    (he : y ^ 2 + x * y = x ^ 3 + (tateCurveAt q hq).a₄ * x + (tateCurveAt q hq).a₆) :
    ∃ a w : R, IsUnit a ∧ w ∈ I ∧ a * w = q ∧
      tateXpair a w q hq = x ∧ tateYpair a w q hq = y := by
  have h4 := tateCurveAt_a4_mem (I := I) q hq
  have h6 := tateCurveAt_a6_mem (I := I) q hq
  have hprod : IsUnit (y * (y + x)) := by
    have h : y * (y + x) = x ^ 3 + ((tateCurveAt q hq).a₄ * x + (tateCurveAt q hq).a₆) := by
      linear_combination he
    rw [h]
    exact isUnit_add_mem (hx.pow 3)
      (Ideal.add_mem _ (Ideal.mul_mem_right _ _ h4) h6)
  have hy : IsUnit y := isUnit_of_mul_isUnit_left hprod
  have hyx : IsUnit (y + x) := isUnit_of_mul_isUnit_right hprod
  have hxy : IsUnit (x + y) := by rwa [add_comm]
  have hxyinv : (x + y) * Ring.inverse (x + y) = 1 := Ring.mul_inverse_cancel _ hxy
  have hcu : IsUnit (y * Ring.inverse (x + y)) := hy.mul (isUnit_ringInverse.mpr hxy)
  have hc1 : IsUnit (1 - y * Ring.inverse (x + y)) := by
    have h : 1 - y * Ring.inverse (x + y) = x * Ring.inverse (x + y) := by
      calc 1 - y * Ring.inverse (x + y)
          = (x + y) * Ring.inverse (x + y) - y * Ring.inverse (x + y) := by rw [hxyinv]
        _ = x * Ring.inverse (x + y) := by ring
    rw [h]
    exact hx.mul (isUnit_ringInverse.mpr hxy)
  obtain ⟨a, hau, ha1, _, hΛ⟩ := exists_tate_lambda_eq q (y * Ring.inverse (x + y)) hq hcu hc1
  have hpt := tate_point_eq_of_lambda x y q hq hx hxy he a (wOf q a) (wOf_mem a hq)
    (mul_wOf q a hau) hau ha1 hΛ
  exact ⟨a, wOf q a, hau, wOf_mem a hq, mul_wOf q a hau, hpt.1, hpt.2⟩

/-! ## ★出典の紐付け(`.src`) -/

def exists_tate_lambda_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——母数座標の方程式の解の存在)",
    sectionId := "genell-def-3-3" }

def tate_surjective.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——葉 (e)、単元の領域での全射性)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
