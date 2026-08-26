import ABC3.Found.GaloisRep.TateSurjOrigin

/-!
# Galois (G6) 第 289 ブロック —— **★★★★★★★★★★全射性が全領域でそろった**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★★★★★到達点

> `E_q(K)` のどのアフィン点 `(x,y)` にも対 `(a,w)` が在って
> `X_K(a,w) = x`、`Y_K(a,w) = y`(`tate_point_surjective`)

★★★3 領域を通じた**一つの主張**である。仮定は局所環の二分法・`I` の素性・
モニック 2 次式の根が `R` に入ること(整閉性)・付値環の性質だけ。

## ★★★★★場合分けは `x` が `R` から来るかどうか

| `x` | 領域 | 使うもの |
|---|---|---|
| `R` の像 | 単元・環帯 | `tate_surjective_integral'`(第 286 の強めた版) |
| `R` の像でない | 原点近傍 | `tate_surjective_origin`(第 288) |

★★`x ∈ R` なら **`y` も自動で `R`**——`y² + x·y − (x³+a₄x+a₆) = 0` がモニックだから。
★★★整閉性を「モニック 2 次式の根」の形で受けたので、`x` からも `y` からも同じ仮定が使える
(第 287 の `γ`、ここの `y`)。

## ★★★第 266 を強める

第 266 の `tate_surjective` は `IsUnit (1 − a)` を結論に出していなかったが、
証明の中では `exists_tate_lambda_eq` からすでに得ている。`K` の水準に上げるには
`tateXK_eq`(第 275)に `IsUnit (1 − a)` が要るので、出すようにした(`tate_surjective'`)。
★**証明が持っている情報を結論に出しておく**——後から要る。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `tate_surjective'` | ★★★★★★第 266 を強めた版 |
| `tate_surjective_integral'` | ★★★★★★★★整な点(強めた版) |
| `tate_point_surjective` | ★★★★★★★★★★**全射性(全領域、`K` の水準)** |
-/

namespace ABC3.Found.GaloisRep

open Complex Real

/-! ## ★★★第 266 を強める -/

section Strong

variable {R : Type} [CommRing R] {I : Ideal R} [IsAdicComplete I R]

/-- ★★★★★★第 266 を強めた版——`1 − a` が単元であることも出す。 -/
theorem tate_surjective' (x y q : R) (hq : q ∈ I) (hx : IsUnit x)
    (he : y ^ 2 + x * y = x ^ 3 + (tateCurveAt q hq).a₄ * x + (tateCurveAt q hq).a₆) :
    ∃ a w : R, IsUnit a ∧ IsUnit (1 - a) ∧ w ∈ I ∧ a * w = q ∧
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
  have hcu : IsUnit (y * Ring.inverse (x + y)) := hy.mul (_root_.isUnit_ringInverse.mpr hxy)
  have hc1 : IsUnit (1 - y * Ring.inverse (x + y)) := by
    have h : 1 - y * Ring.inverse (x + y) = x * Ring.inverse (x + y) := by
      calc 1 - y * Ring.inverse (x + y)
          = (x + y) * Ring.inverse (x + y) - y * Ring.inverse (x + y) := by rw [hxyinv]
        _ = x * Ring.inverse (x + y) := by ring
    rw [h]
    exact hx.mul (_root_.isUnit_ringInverse.mpr hxy)
  obtain ⟨a, hau, ha1, _, hΛ⟩ := exists_tate_lambda_eq q (y * Ring.inverse (x + y)) hq hcu hc1
  have hpt := tate_point_eq_of_lambda x y q hq hx hxy he a (wOf q a) (wOf_mem a hq)
    (mul_wOf q a hau) hau ha1 hΛ
  exact ⟨a, wOf q a, hau, ha1, wOf_mem a hq, mul_wOf q a hau, hpt.1, hpt.2⟩

/-- ★★★★★★★★**整な点についての全射性(強めた版)**。 -/
theorem tate_surjective_integral' (hloc : ∀ x : R, IsUnit x ∨ x ∈ I) (hI : I.IsPrime)
    (x y q : R) (hq : q ∈ I)
    (he : y ^ 2 + x * y = x ^ 3 + (tateCurveAt q hq).a₄ * x + (tateCurveAt q hq).a₆) :
    ∃ a w : R, a * w = q ∧ w ∈ I ∧ IsUnit (1 - a) ∧
      tateXpair a w q hq = x ∧ tateYpair a w q hq = y := by
  rcases hloc x with hxu | hxm
  · obtain ⟨a, w, _, ha1, hwI, haw, hX, hY⟩ := tate_surjective' x y q hq hxu he
    exact ⟨a, w, haw, hwI, ha1, hX, hY⟩
  · obtain ⟨a, haI, w, hwI, haw, hX, hY⟩ :=
      tate_surjective_annulus x y q hq hxm (tate_y_mem_of_x_mem hI x y q hq hxm he) he
    exact ⟨a, w, haw, hwI, isUnit_one_sub haI, hX, hY⟩

end Strong

/-! ## ★★★★★★★★★★全射性(全領域、`K` の水準) -/

section All

variable {R K : Type} [CommRing R] {I : Ideal R} [IsAdicComplete I R] [Field K] [Algebra R K]

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★★★★**全射性(全領域、`K` の水準)**——3 領域を通じた一つの主張。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem tate_point_surjective (hinj : Function.Injective (algebraMap R K))
    (hloc : ∀ x : R, IsUnit x ∨ x ∈ I) (hI : I.IsPrime)
    (hquad : ∀ (t : K) (b c : R),
      t ^ 2 + algebraMap R K b * t + algebraMap R K c = 0 → ∃ r : R, algebraMap R K r = t)
    (hvalring : ∀ t : K, t ≠ 0 →
      (∃ r : R, algebraMap R K r = t) ∨ (∃ r ∈ I, algebraMap R K r = t⁻¹))
    (q : R) (hq : q ∈ I) (x y : K)
    (he : y ^ 2 + x * y = x ^ 3 + algebraMap R K ((tateCurveAt q hq).a₄) * x
      + algebraMap R K ((tateCurveAt q hq).a₆)) :
    ∃ a w : R, a * w = q ∧ w ∈ I ∧ algebraMap R K (1 - a) ≠ 0 ∧
      tateXK a w q hq = x ∧ tateYK (K := K) a w q hq = y := by
  set A := (tateCurveAt q hq).a₄ with hA
  set B := (tateCurveAt q hq).a₆ with hB
  have hint : ∀ x₀ : R, algebraMap R K x₀ = x →
      ∃ a w : R, a * w = q ∧ w ∈ I ∧ algebraMap R K (1 - a) ≠ 0 ∧
        tateXK a w q hq = x ∧ tateYK (K := K) a w q hq = y := by
    intro x₀ hx₀
    obtain ⟨y₀, hy₀⟩ : ∃ y₀ : R, algebraMap R K y₀ = y := by
      refine hquad y x₀ (-(x₀ ^ 3 + A * x₀ + B)) ?_
      rw [map_neg, map_add, map_add, map_mul, map_pow, hx₀]
      linear_combination he
    have heR : y₀ ^ 2 + x₀ * y₀ = x₀ ^ 3 + A * x₀ + B := by
      refine hinj ?_
      rw [map_add, map_add, map_add, map_mul, map_mul, map_pow, map_pow, hx₀, hy₀]
      exact he
    obtain ⟨a, w, haw, hwI, ha1, hX, hY⟩ := tate_surjective_integral' hloc hI x₀ y₀ q hq heR
    refine ⟨a, w, haw, hwI, (ha1.map (algebraMap R K)).ne_zero, ?_, ?_⟩
    · rw [tateXK_eq a w q hq ha1, hX, hx₀]
    · rw [tateYK_eq a w q hq ha1, hY, hy₀]
  by_cases hx0 : x = 0
  · exact hint 0 (by rw [map_zero, hx0])
  · rcases hvalring x hx0 with ⟨x₀, hx₀⟩ | ⟨α₀, hα, hαx⟩
    · exact hint x₀ hx₀
    · obtain ⟨u, hu, hne, hX, hY⟩ :=
        tate_surjective_origin hinj hI hquad q hq x y hx0 α₀ hα hαx he
      exact ⟨u, wOf q u, mul_wOf q u (isUnit_of_one_sub_mem hu), wOf_mem u hq, hne, hX, hY⟩

end All

/-! ## ★出典の紐付け(`.src`) -/

def tate_point_surjective.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——全射性(全領域))",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
