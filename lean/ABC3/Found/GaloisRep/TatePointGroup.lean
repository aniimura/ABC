import ABC3.Found.GaloisRep.TateDVR

/-!
# Galois (G6) 第 261 ブロック —— **★★★★★★★★★Tate 点の 3 つ組は群として 0**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★★★★到達点

> `u v w = q`、`u, v, w ∈ I` が相異なるなら
> `P(u) + P(v) + P(w) = 0` (`E_q(K)` の中で)

これが**葉 (c)(準同型)の最終形**である。第 256(共線性)・第 257(群法則)・
第 259(単射性)・第 260(`x` 座標の相異性)がここで合流する。

## ★★基底変換

`tateXpair`・`tateYpair` は完備環 `R` の値だが、mathlib の群法則は**体**を要求する。
そこで `R → K`(分数体)へ送る:

    (tateCurveAt q hq).map (algebraMap R K)

★係数は `a₁ = 1`、`a₂ = a₃ = 0`、`a₄ = f(a₄)`、`a₆ = f(a₆)` になるので、
葉 (b) の `tate_equation` をそのまま送れば `Equation` が出る(`tate_equation_map`)。
★★`Δ ≠ 0` なら `equation_iff_nonsingular_of_Δ_ne_zero` で `Nonsingular` に上がる。

## ★★★相異性は 3 回の当てはめで

第 260 の `tateXpair_ne_of_ne` は `(u, v, w)` の順に依存するので、
`(u,v,w)`・`(u,w,v)`・`(v,w,u)` の 3 通りに当てる。
★相方の積は可換なので `mul_comm` の書き換えだけで揃う。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `tate_equation_map` | ★★★★★基底変換した曲線の上にある |
| `nonsingular_tate_point` | ★★★★★`Δ ≠ 0` なら非特異 |
| `tate_points_add_eq_zero` | ★★★★★★★★★**3 つ組は群として 0** |
-/

namespace ABC3.Found.GaloisRep

open Complex Real WeierstrassCurve WeierstrassCurve.Affine

/-! ## ★★★★★基底変換 -/

section BaseChange

variable {R K : Type} [CommRing R] {I : Ideal R} [IsAdicComplete I R]
  [Field K] [Algebra R K]

/-- ★★★★★**Tate 点は基底変換した曲線の上にある**。 -/
theorem tate_equation_map (a w q : R) (hq : q ∈ I) (haw : a * w = q)
    (ha : a ∈ I) (hw : w ∈ I) :
    ((tateCurveAt q hq).map (algebraMap R K)).toAffine.Equation
      (algebraMap R K (tateXpair a w q hq)) (algebraMap R K (tateYpair a w q hq)) := by
  rw [WeierstrassCurve.Affine.equation_iff]
  have h := tate_equation a w q hq haw (isUnit_one_sub ha) (isUnit_one_sub hw)
  have h2 := congrArg (algebraMap R K) h
  simp only [map_add, map_mul, map_pow] at h2
  have e1 : ((tateCurveAt q hq).map (algebraMap R K)).toAffine.a₁ = 1 := by
    show algebraMap R K ((tateCurveAt q hq).a₁) = 1
    rw [show (tateCurveAt q hq).a₁ = 1 by
      simp [tateCurveAt, tateCurve, WeierstrassCurve.map], map_one]
  have e2 : ((tateCurveAt q hq).map (algebraMap R K)).toAffine.a₂ = 0 := by
    show algebraMap R K ((tateCurveAt q hq).a₂) = 0
    rw [show (tateCurveAt q hq).a₂ = 0 by
      simp [tateCurveAt, tateCurve, WeierstrassCurve.map], map_zero]
  have e3 : ((tateCurveAt q hq).map (algebraMap R K)).toAffine.a₃ = 0 := by
    show algebraMap R K ((tateCurveAt q hq).a₃) = 0
    rw [show (tateCurveAt q hq).a₃ = 0 by
      simp [tateCurveAt, tateCurve, WeierstrassCurve.map], map_zero]
  have e4 : ((tateCurveAt q hq).map (algebraMap R K)).toAffine.a₄
      = algebraMap R K ((tateCurveAt q hq).a₄) := rfl
  have e6 : ((tateCurveAt q hq).map (algebraMap R K)).toAffine.a₆
      = algebraMap R K ((tateCurveAt q hq).a₆) := rfl
  rw [e1, e2, e3, e4, e6]
  linear_combination h2

/-- ★★★★★**Tate 点は非特異**——`Δ ≠ 0` から。 -/
theorem nonsingular_tate_point (a w q : R) (hq : q ∈ I) (haw : a * w = q)
    (ha : a ∈ I) (hw : w ∈ I)
    (hΔ : ((tateCurveAt q hq).map (algebraMap R K)).toAffine.Δ ≠ 0) :
    ((tateCurveAt q hq).map (algebraMap R K)).toAffine.Nonsingular
      (algebraMap R K (tateXpair a w q hq)) (algebraMap R K (tateYpair a w q hq)) :=
  (WeierstrassCurve.Affine.equation_iff_nonsingular_of_Δ_ne_zero hΔ).mp
    (tate_equation_map a w q hq haw ha hw)

end BaseChange

/-! ## ★★★★★★★★★3 つ組は群として 0 -/

section Group

variable {R K : Type} [CommRing R] [IsDomain R] {I : Ideal R} [IsAdicComplete I R]
  [Field K] [DecidableEq K] [Algebra R K]

set_option maxHeartbeats 800000 in
/-- ★★★★★★★★★**Tate 点の 3 つ組は群として 0 になる**——`u v w = q` のとき。

★第 256(共線性)・第 257(群法則)・第 259(単射性)・第 260(`x` の相異性)が合流する。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem tate_points_add_eq_zero (hinj : Function.Injective (algebraMap R K))
    (u v w q : R) (hq : q ∈ I) (hu : u ∈ I) (hv : v ∈ I) (hw : w ∈ I)
    (huvw : u * v * w = q) (hu0 : u ≠ 0) (hv0 : v ≠ 0)
    (huv : u ≠ v) (huw : u ≠ w) (hvw : v ≠ w)
    (hcp : ∀ i, IsUnit (1 - collPts u v w i))
    (n₁ : ((tateCurveAt q hq).map (algebraMap R K)).toAffine.Nonsingular
      (algebraMap R K (tateXpair u (v * w) q hq)) (algebraMap R K (tateYpair u (v * w) q hq)))
    (n₂ : ((tateCurveAt q hq).map (algebraMap R K)).toAffine.Nonsingular
      (algebraMap R K (tateXpair v (u * w) q hq)) (algebraMap R K (tateYpair v (u * w) q hq)))
    (n₃ : ((tateCurveAt q hq).map (algebraMap R K)).toAffine.Nonsingular
      (algebraMap R K (tateXpair w (u * v) q hq)) (algebraMap R K (tateYpair w (u * v) q hq))) :
    Point.some _ _ n₁ + Point.some _ _ n₂ + Point.some _ _ n₃ = 0 := by
  have hd0 := TateCollinearSection.tate_collinear u v w q hq huvw hcp
  have hd := congrArg (algebraMap R K) hd0
  simp only [map_add, map_mul, map_sub, map_zero] at hd
  have h12 : algebraMap R K (tateXpair u (v * w) q hq)
      ≠ algebraMap R K (tateXpair v (u * w) q hq) :=
    fun h => tateXpair_ne_of_ne u v w q hq hu hv hw huvw hu0 huv (hinj h)
  have h13 : algebraMap R K (tateXpair u (v * w) q hq)
      ≠ algebraMap R K (tateXpair w (u * v) q hq) := by
    have hbase := tateXpair_ne_of_ne u w v q hq hu hw hv (by rw [← huvw]; ring) hu0 huw
    rw [show w * v = v * w from mul_comm w v] at hbase
    exact fun h => hbase (hinj h)
  have h23 : algebraMap R K (tateXpair v (u * w) q hq)
      ≠ algebraMap R K (tateXpair w (u * v) q hq) := by
    have hbase := tateXpair_ne_of_ne v w u q hq hv hw hu (by rw [← huvw]; ring) hv0 hvw
    rw [show w * u = u * w from mul_comm w u, show v * u = u * v from mul_comm v u] at hbase
    exact fun h => hbase (hinj h)
  exact add_add_eq_zero_of_collDet n₁ n₂ n₃ hd h12 h13 h23

end Group

/-! ## ★出典の紐付け(`.src`) -/

def tate_equation_map.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——基底変換した曲線の上にある)",
    sectionId := "genell-def-3-3" }

def tate_points_add_eq_zero.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——3 点の和は群として 0)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
