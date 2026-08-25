import ABC3.Found.GaloisRep.TateQUnique

/-!
# Galois (G6) 第 272 ブロック —— **★★★★★★★★★混合の場合の群法則**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★★★★到達点

第 261 の群法則は **`u, v, w` がすべて `I` の元**(環帯)の場合だった。
本ブロックは **`u` が単元・`v, w ∈ I`** の場合を閉じる。

`v(u) + v(v) + v(w) = v(q) > 0` なので、3 つのうち少なくとも 1 つは `I` に入る。
したがって場合は「1 つが `I`」「2 つが `I`」「3 つとも `I`」の 3 通りで、
本ブロックは **2 つが `I`** の場合にあたる。

## ★★★★★★単元は `I` の元と自動的に異なる

群法則には 3 点の `x` 座標の相異性が要る。★★これまでは単射性(第 259・第 267)を
使ってきたが、**混合の場合は相異性がただで出る**:

    X₁ = X(u, vw) は単元(第 264)、  X₂ = X(v, uw)・X₃ = X(w, uv) は I の元(第 271)

`I` は真のイデアルだから、単元は `I` に入らない(`ne_of_isUnit_of_mem`)。
★★★**「領域が違えば座標が違う」**——単射性を経由しないでよい。
残る `X₂ ≠ X₃` だけが第 260 の環帯の単射性を要る。

## ★★★★★★第 260 を緩める

第 260 の `tateXpair_ne_of_ne` は `u, v, w` の 3 つとも `I` を要求していたが、
証明が実際に使うのは**相方の所属**(`v*w ∈ I`、`u*w ∈ I`)と `w ≠ 1` だけである。
`u` が単元でも `w ∈ I` なら `w*u ∈ I` なので通る(`tateXpair_ne_of_ne'`)。
★`w ≠ 1` は `IsUnit (1 − u)` から出る(`1 − 1 = 0` は単元でない)。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `ne_of_isUnit_of_mem` | ★★★★**単元は `I` の元と異なる** |
| `tateXpair_ne_of_ne'` | ★★★★★★第 260 の緩めた版 |
| `tate_three_ne_mixed` | ★★★★★★**混合の場合の相異性** |
| `collUnits_mixed` | ★★★混合の場合の 6 点の単元性 |
| `tate_equation_map'`・`nonsingular_tate_point'` | ★★★★★第 261 の緩めた版 |
| `tate_points_add_eq_zero_mixed` | ★★★★★★★★★**混合の場合の群法則** |
-/

namespace ABC3.Found.GaloisRep

open Complex Real WeierstrassCurve WeierstrassCurve.Affine

variable {R : Type} [CommRing R] {I : Ideal R}

/-! ## ★★★★単元は `I` の元と異なる -/

/-- ★★★★**単元は `I` の元と異なる**——`I` は真のイデアルだから。

★`I ≠ ⊤` は `IsAdicComplete` から出る:`1 ∈ I` なら `1 − 1 = 0` が単元になってしまう。 -/
theorem ne_of_isUnit_of_mem [IsAdicComplete I R] [Nontrivial R] {x y : R}
    (hx : IsUnit x) (hy : y ∈ I) : x ≠ y := by
  intro h
  rw [h] at hx
  obtain ⟨v, hv⟩ := hx
  have h1 : (1 : R) ∈ I := by
    have h2 := Ideal.mul_mem_left I ((v⁻¹ : Rˣ) : R) hy
    rw [← hv, ← Units.val_mul, inv_mul_cancel, Units.val_one] at h2
    exact h2
  have h3 := isUnit_one_sub (I := I) h1
  rw [sub_self] at h3
  exact not_isUnit_zero h3

/-! ## ★★★★★★混合の場合の相異性 -/

/-- ★★★★★★第 260 の緩めた版——`u ≠ 0`・`u ≠ v` と**相方の所属**だけでよい。 -/
theorem tateXpair_ne_of_ne' [IsAdicComplete I R] [IsDomain R] (u v w q : R) (hq : q ∈ I)
    (hu : u ∈ I) (hv : v ∈ I) (hvw : v * w ∈ I) (huw : u * w ∈ I) (hw1 : w ≠ 1)
    (hu0 : u ≠ 0) (huv : u ≠ v) (huvw : u * v * w = q) :
    tateXpair u (v * w) q hq ≠ tateXpair v (u * w) q hq := by
  intro heq
  have h1 : u * (v * w) = q := by rw [← huvw]; ring
  have h2 : v * (u * w) = q := by rw [← huvw]; ring
  rcases tate_inj_X u (v * w) v (u * w) q hq hu hvw hv huw h1 h2 heq with ⟨h, _⟩ | ⟨h, _⟩
  · exact huv h
  · have hz : u * (1 - w) = 0 := by linear_combination h
    rcases mul_eq_zero.1 hz with h' | h'
    · exact hu0 h'
    · exact hw1 (by linear_combination -h')

/-- ★★★★★★**混合の場合の相異性**——単元 1 つと `I` の元 2 つ。

★`X₁` は単元・`X₂`, `X₃` は `I` の元なので、`X₁ ≠ X₂`・`X₁ ≠ X₃` はただで出る。 -/
theorem tate_three_ne_mixed [IsAdicComplete I R] [IsDomain R] (u v w q : R) (hq : q ∈ I)
    (hu1 : IsUnit u) (hu2 : IsUnit (1 - u)) (hv : v ∈ I) (hw : w ∈ I) (huvw : u * v * w = q)
    (hv0 : v ≠ 0) (hvw : v ≠ w) :
    tateXpair u (v * w) q hq ≠ tateXpair v (u * w) q hq
      ∧ tateXpair u (v * w) q hq ≠ tateXpair w (u * v) q hq
      ∧ tateXpair v (u * w) q hq ≠ tateXpair w (u * v) q hq := by
  have hvwI : v * w ∈ I := Ideal.mul_mem_right _ _ hv
  have huwI : u * w ∈ I := Ideal.mul_mem_left _ _ hw
  have huvI : u * v ∈ I := Ideal.mul_mem_left _ _ hv
  have hX1 : IsUnit (tateXpair u (v * w) q hq) :=
    isUnit_tateXpair u (v * w) q hq hvwI hu1 hu2
  have hX2 : tateXpair v (u * w) q hq ∈ I := tateXpair_mem v (u * w) q hq hv huwI
  have hX3 : tateXpair w (u * v) q hq ∈ I := tateXpair_mem w (u * v) q hq hw huvI
  have hu1' : u ≠ 1 := by
    intro h
    rw [h, sub_self] at hu2
    exact not_isUnit_zero hu2
  refine ⟨ne_of_isUnit_of_mem hX1 hX2, ne_of_isUnit_of_mem hX1 hX3, ?_⟩
  have hbase := tateXpair_ne_of_ne' v w u q hq hv hw
    (by rw [show w * u = u * w from mul_comm w u]; exact huwI)
    (by rw [show v * u = u * v from mul_comm v u]; exact huvI)
    hu1' hv0 hvw (by rw [← huvw]; ring)
  rw [show w * u = u * w from mul_comm w u, show v * u = u * v from mul_comm v u] at hbase
  exact hbase

/-! ## ★★★混合の場合の 6 点の単元性 -/

theorem collUnits_mixed [IsAdicComplete I R] {u v w : R} (hu2 : IsUnit (1 - u)) (hv : v ∈ I)
    (hw : w ∈ I) : ∀ i, IsUnit (1 - collPts u v w i) := by
  intro i
  fin_cases i
  · exact hu2
  · exact isUnit_one_sub hv
  · exact isUnit_one_sub hw
  · exact isUnit_one_sub (Ideal.mul_mem_right _ _ hv)
  · exact isUnit_one_sub (Ideal.mul_mem_left _ _ hw)
  · exact isUnit_one_sub (Ideal.mul_mem_left _ _ hv)

/-! ## ★★★★★基底変換(緩めた版) -/

section BaseChange

variable {K : Type} [IsAdicComplete I R] [Field K] [Algebra R K]

/-- ★★★★★第 261 の緩めた版——`1 − a`・`1 − w` が単元であればよい。 -/
theorem tate_equation_map' (a w q : R) (hq : q ∈ I) (haw : a * w = q)
    (ha : IsUnit (1 - a)) (hw : IsUnit (1 - w)) :
    ((tateCurveAt q hq).map (algebraMap R K)).toAffine.Equation
      (algebraMap R K (tateXpair a w q hq)) (algebraMap R K (tateYpair a w q hq)) := by
  rw [WeierstrassCurve.Affine.equation_iff]
  have h := tate_equation a w q hq haw ha hw
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

theorem nonsingular_tate_point' (a w q : R) (hq : q ∈ I) (haw : a * w = q)
    (ha : IsUnit (1 - a)) (hw : IsUnit (1 - w))
    (hΔ : ((tateCurveAt q hq).map (algebraMap R K)).toAffine.Δ ≠ 0) :
    ((tateCurveAt q hq).map (algebraMap R K)).toAffine.Nonsingular
      (algebraMap R K (tateXpair a w q hq)) (algebraMap R K (tateYpair a w q hq)) :=
  (WeierstrassCurve.Affine.equation_iff_nonsingular_of_Δ_ne_zero hΔ).mp
    (tate_equation_map' a w q hq haw ha hw)

end BaseChange

/-! ## ★★★★★★★★★混合の場合の群法則 -/

section Group

variable {K : Type} [IsDomain R] [IsAdicComplete I R] [Field K] [DecidableEq K] [Algebra R K]

set_option maxHeartbeats 800000 in
/-- ★★★★★★★★★**混合の場合の群法則**——単元 1 つと `I` の元 2 つ。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem tate_points_add_eq_zero_mixed (hinj : Function.Injective (algebraMap R K))
    (u v w q : R) (hq : q ∈ I) (hu1 : IsUnit u) (hu2 : IsUnit (1 - u))
    (hv : v ∈ I) (hw : w ∈ I) (huvw : u * v * w = q) (hv0 : v ≠ 0) (hvw : v ≠ w)
    (n₁ : ((tateCurveAt q hq).map (algebraMap R K)).toAffine.Nonsingular
      (algebraMap R K (tateXpair u (v * w) q hq)) (algebraMap R K (tateYpair u (v * w) q hq)))
    (n₂ : ((tateCurveAt q hq).map (algebraMap R K)).toAffine.Nonsingular
      (algebraMap R K (tateXpair v (u * w) q hq)) (algebraMap R K (tateYpair v (u * w) q hq)))
    (n₃ : ((tateCurveAt q hq).map (algebraMap R K)).toAffine.Nonsingular
      (algebraMap R K (tateXpair w (u * v) q hq)) (algebraMap R K (tateYpair w (u * v) q hq))) :
    Point.some _ _ n₁ + Point.some _ _ n₂ + Point.some _ _ n₃ = 0 := by
  have hd0 := TateCollinearSection.tate_collinear u v w q hq huvw (collUnits_mixed hu2 hv hw)
  have hd := congrArg (algebraMap R K) hd0
  simp only [map_add, map_mul, map_sub, map_zero] at hd
  obtain ⟨h12, h13, h23⟩ := tate_three_ne_mixed u v w q hq hu1 hu2 hv hw huvw hv0 hvw
  exact add_add_eq_zero_of_collDet n₁ n₂ n₃ hd (fun h => h12 (hinj h))
    (fun h => h13 (hinj h)) (fun h => h23 (hinj h))

end Group

/-! ## ★出典の紐付け(`.src`) -/

def ne_of_isUnit_of_mem.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——単元は I の元と異なる)",
    sectionId := "genell-def-3-3" }

def tate_points_add_eq_zero_mixed.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——混合の場合の群法則)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
