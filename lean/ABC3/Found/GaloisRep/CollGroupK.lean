import ABC3.Found.GaloisRep.CollOriginUniv

/-!
# Galois (G6) 第 279 ブロック —— **★★★★★★★★★★原点近傍を含む群法則**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★★★★★到達点

> `u·v·w = q`、相方 `vw, uw, uv` が `I` の元、母数の `1 − u, 1 − v, 1 − w` が
> `K` で `0` でない —— このとき **3 点の和は群として 0**(`tate_points_add_eq_zero_K`)

★★★★母数そのものは `I` に居なくても、単元でなくてもよい。**`≡ 1` でもよい**。
第 278 の `collDefectE_eq_zero` を `K` に落として割り算するだけである。

## ★★★★★どの対も片側は必ず単元

`p·r = q ∈ 𝔪` なので、局所環では **`p`, `r` の少なくとも一方が `𝔪` に入る**。
`p ≡ 1` なら `p` は単元だから `r = q/p ∈ 𝔪`、したがって `1 − r` は単元。
★★★**対の両側が同時に `≡ 1` になることはない**——これが分母払いを片側で済ませる根拠。

    tateXpairEE p r q = (1−p)(1−r)³ · tateXpairE p r q        (`1 − r` が単元なら)

★これで第 275 の `X_K = XE/(1−p)²` と第 277 の両側払いが接がる。

## ★★★★正規化した 3 つ組では相方は必ず `𝔪` に入る

基本領域 `0 ≤ v(p_i) < v(q)` で `Σ v(p_i) = v(q)` なら
`v(相方_i) = v(q) − v(p_i) > 0` ★★したがって相方はすべて `𝔪` の元である。
★同時に「単元は高々 1 つ」も出る(2 つが 0 なら第 3 が `v(q)` になって領域外)。

| `I` に入る母数の個数 | 群法則 |
|---|---|
| 2 個以上(正規化した 3 つ組) | **本ブロック**(原点近傍を含む) |
| 1 個(単元 2 つ、`Σ v = 0` の場合) | 第 273 |

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `tateXpairEE_eq_one_side` | ★★★★★両側払いと片側払いの接ぎ木 |
| `map_tateXpairEE`・`map_tateYpairEE` | ★★★★`K` での `XE = M·X` |
| `collDet_K_eq_zero` | ★★★★★★★★★★**`K` の水準の共線性** |
| `tate_points_add_eq_zero_K` | ★★★★★★★★★★**原点近傍を含む群法則** |
-/

namespace ABC3.Found.GaloisRep

open Complex Real WeierstrassCurve WeierstrassCurve.Affine

/-! ## ★★★★★両側払いと片側払いの接ぎ木 -/

section Bridge

variable {R : Type} [CommRing R] {I : Ideal R}

/-- ★★★★★相方が `I` の元なら、両側払いは片側払いの `(1−p)(1−r)³` 倍。 -/
theorem tateXpairEE_eq_one_side [IsAdicComplete I R] (p r q : R) (hq : q ∈ I)
    (hr : IsUnit (1 - r)) :
    tateXpairEE p r q hq = (1 - p) * (1 - r) ^ 3 * tateXpairE p r q hq := by
  rw [tateXpairEE, tateXpairE, tateM]
  linear_combination (-((1 - p) ^ 3 * (1 - r))) * mul_tateXterm' hr

theorem tateYpairEE_eq_one_side [IsAdicComplete I R] (p r q : R) (hq : q ∈ I)
    (hr : IsUnit (1 - r)) :
    tateYpairEE p r q hq = (1 - r) ^ 3 * tateYpairE p r q hq := by
  rw [tateYpairEE, tateYpairE, tateM]
  linear_combination ((1 - p) ^ 3 * (1 - r)) * mul_tateXterm' hr
    + ((1 - p) ^ 3) * mul_tateYterm' hr

end Bridge

/-! ## ★★★★`K` での `XE = M·X` -/

section Field

variable {R K : Type} [CommRing R] {I : Ideal R} [IsAdicComplete I R] [Field K] [Algebra R K]

theorem map_tateXpairEE (p r q : R) (hq : q ∈ I) (hr : IsUnit (1 - r))
    (hp : algebraMap R K (1 - p) ≠ 0) :
    algebraMap R K (tateXpairEE p r q hq)
      = algebraMap R K (tateM p r) * tateXK p r q hq := by
  rw [tateXpairEE_eq_one_side p r q hq hr, tateXK, tateM]
  simp only [map_mul, map_pow]
  field_simp

theorem map_tateYpairEE (p r q : R) (hq : q ∈ I) (hr : IsUnit (1 - r))
    (hp : algebraMap R K (1 - p) ≠ 0) :
    algebraMap R K (tateYpairEE p r q hq)
      = algebraMap R K (tateM p r) * tateYK p r q hq := by
  rw [tateYpairEE_eq_one_side p r q hq hr, tateYK, tateM]
  simp only [map_mul, map_pow]
  field_simp

theorem map_tateM_ne_zero (p r : R) (hp : algebraMap R K (1 - p) ≠ 0)
    (hr : algebraMap R K (1 - r) ≠ 0) : algebraMap R K (tateM p r) ≠ 0 := by
  rw [tateM]
  simp only [map_mul, map_pow]
  exact mul_ne_zero (pow_ne_zero _ hp) (pow_ne_zero _ hr)

/-! ## ★★★★★★★★★★`K` の水準の共線性 -/

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★★★★**`K` の水準の共線性**——母数が `≡ 1` でもよい。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem collDet_K_eq_zero (u v w q : R) (hq : q ∈ I) (huvw : u * v * w = q)
    (hvw : v * w ∈ I) (huw : u * w ∈ I) (huv : u * v ∈ I)
    (hu : algebraMap R K (1 - u) ≠ 0) (hv : algebraMap R K (1 - v) ≠ 0)
    (hw : algebraMap R K (1 - w) ≠ 0) :
    tateXK u (v * w) q hq * (tateYK (K := K) v (u * w) q hq - tateYK (K := K) w (u * v) q hq)
      + tateXK v (u * w) q hq
        * (tateYK (K := K) w (u * v) q hq - tateYK (K := K) u (v * w) q hq)
      + tateXK w (u * v) q hq
        * (tateYK (K := K) u (v * w) q hq - tateYK (K := K) v (u * w) q hq) = 0 := by
  have uvw : IsUnit (1 - v * w) := isUnit_one_sub hvw
  have uuw : IsUnit (1 - u * w) := isUnit_one_sub huw
  have uuv : IsUnit (1 - u * v) := isUnit_one_sub huv
  have nvw : algebraMap R K (1 - v * w) ≠ 0 := (uvw.map (algebraMap R K)).ne_zero
  have nuw : algebraMap R K (1 - u * w) ≠ 0 := (uuw.map (algebraMap R K)).ne_zero
  have nuv : algebraMap R K (1 - u * v) ≠ 0 := (uuv.map (algebraMap R K)).ne_zero
  have h0 := congrArg (algebraMap R K) (collDefectE_eq_zero u v w q hq huvw)
  rw [collDefectE, map_zero] at h0
  simp only [map_add, map_sub, map_mul] at h0
  rw [map_tateXpairEE u (v * w) q hq uvw hu, map_tateYpairEE u (v * w) q hq uvw hu,
    map_tateXpairEE v (u * w) q hq uuw hv, map_tateYpairEE v (u * w) q hq uuw hv,
    map_tateXpairEE w (u * v) q hq uuv hw, map_tateYpairEE w (u * v) q hq uuv hw] at h0
  have hM : (algebraMap R K (tateM u (v * w)) * algebraMap R K (tateM v (u * w))
      * algebraMap R K (tateM w (u * v))) ≠ 0 :=
    mul_ne_zero (mul_ne_zero (map_tateM_ne_zero u (v * w) hu nvw)
      (map_tateM_ne_zero v (u * w) hv nuw)) (map_tateM_ne_zero w (u * v) hw nuv)
  have hkey : (algebraMap R K (tateM u (v * w)) * algebraMap R K (tateM v (u * w))
      * algebraMap R K (tateM w (u * v)))
      * (tateXK u (v * w) q hq
          * (tateYK (K := K) v (u * w) q hq - tateYK (K := K) w (u * v) q hq)
        + tateXK v (u * w) q hq
          * (tateYK (K := K) w (u * v) q hq - tateYK (K := K) u (v * w) q hq)
        + tateXK w (u * v) q hq
          * (tateYK (K := K) u (v * w) q hq - tateYK (K := K) v (u * w) q hq)) = 0 := by
    linear_combination h0
  exact (mul_eq_zero.1 hkey).resolve_left hM

end Field

/-! ## ★★★★★★★★★★原点近傍を含む群法則 -/

section Group

variable {R K : Type} [CommRing R] {I : Ideal R} [IsAdicComplete I R] [Field K]
  [DecidableEq K] [Algebra R K]

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★★★★**原点近傍を含む群法則**——母数が `≡ 1` でもよい。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem tate_points_add_eq_zero_K (u v w q : R) (hq : q ∈ I) (huvw : u * v * w = q)
    (hvw : v * w ∈ I) (huw : u * w ∈ I) (huv : u * v ∈ I)
    (hu : algebraMap R K (1 - u) ≠ 0) (hv : algebraMap R K (1 - v) ≠ 0)
    (hw : algebraMap R K (1 - w) ≠ 0)
    (h12 : tateXK u (v * w) q hq ≠ (tateXK v (u * w) q hq : K))
    (h13 : tateXK u (v * w) q hq ≠ (tateXK w (u * v) q hq : K))
    (h23 : tateXK v (u * w) q hq ≠ (tateXK w (u * v) q hq : K))
    (n₁ : ((tateCurveAt q hq).map (algebraMap R K)).toAffine.Nonsingular
      (tateXK u (v * w) q hq) (tateYK (K := K) u (v * w) q hq))
    (n₂ : ((tateCurveAt q hq).map (algebraMap R K)).toAffine.Nonsingular
      (tateXK v (u * w) q hq) (tateYK (K := K) v (u * w) q hq))
    (n₃ : ((tateCurveAt q hq).map (algebraMap R K)).toAffine.Nonsingular
      (tateXK w (u * v) q hq) (tateYK (K := K) w (u * v) q hq)) :
    Point.some _ _ n₁ + Point.some _ _ n₂ + Point.some _ _ n₃ = 0 :=
  add_add_eq_zero_of_collDet n₁ n₂ n₃
    (collDet_K_eq_zero u v w q hq huvw hvw huw huv hu hv hw) h12 h13 h23

end Group

/-! ## ★出典の紐付け(`.src`) -/

def collDet_K_eq_zero.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——K の水準の共線性)",
    sectionId := "genell-def-3-3" }

def tate_points_add_eq_zero_K.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——原点近傍を含む群法則)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
