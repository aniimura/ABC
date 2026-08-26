import ABC3.Found.GaloisRep.TateInversion

/-!
# Galois (G6) 第 282 ブロック —— **★★★★★★3 領域を通じた一つの点の式**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★一つの式で 3 領域を書く

    tatePtPair a w q hq haw hw hne hΔ := Point.some (X_K(a,w)) (Y_K(a,w)) …

仮定は **`a·w = q`、`1 − w` が単元、`1 − a` が `K` で `0` でない**——これだけである。

| 領域 | `a` | `w` | `1 − a` | `1 − w` |
|---|---|---|---|---|
| 単元 | 単元、`≢ 1` | `𝔪` | 単元 | 単元 |
| 環帯 | `𝔪` | `𝔪` | 単元 | 単元 |
| 原点近傍 | `≡ 1` | `𝔪` | **非単元だが `≠ 0`** | 単元 |

★★★対の片側が必ず `𝔪` に入る(第 279)ので `1 − w` が単元という仮定は
**正規化した対では自動**である。

## ★★★★★★★★★★この一つの式で群法則も反転則も書ける

| 定理 | 内容 |
|---|---|
| `tatePtPair_add_add_eq_zero` | 3 点の和は群として 0(第 279) |
| `tatePtPair_swap` | `P(q/u) = −P(u)`(環帯・単元) |
| `tatePtPair_ringInverse` | `P(1/u) = −P(u)`(全領域、第 281) |
| `tatePtPair_ne_zero` | 母数から来る点は原点でない |

★★反転則が 2 つあるのは、類 `c⁻¹` の正規化代表元が
`v(u) > 0` なら `q/u`、`v(u) = 0` なら `1/u` だからである。

## ★★次に要るもの

`Φ : Kˣ/q^ℤ → E_q(K)` を作るには、この点が**類だけで決まる**ことが要る。
第 213 の `normRep`(基本領域の代表元)と `pair_eq_of_same_class` がその道具である。
-/

namespace ABC3.Found.GaloisRep

open Complex Real WeierstrassCurve WeierstrassCurve.Affine

variable {R K : Type} [CommRing R] {I : Ideal R} [IsAdicComplete I R] [Field K]
  [DecidableEq K] [Algebra R K]

/-- ★★★★★★**対 `(a,w)` から作る `E_q(K)` の点**——3 領域を通じて一つの式。 -/
noncomputable def tatePtPair (a w q : R) (hq : q ∈ I) (haw : a * w = q)
    (hw : IsUnit (1 - w)) (hne : algebraMap R K (1 - a) ≠ 0)
    (hΔ : ((tateCurveAt q hq).map (algebraMap R K)).toAffine.Δ ≠ 0) :
    ((tateCurveAt q hq).map (algebraMap R K)).toAffine.Point :=
  Point.some _ _ (nonsingular_tateK a w q hq haw hw hne hΔ)

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★★★★**3 点の和は群として 0**(対の形)。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem tatePtPair_add_add_eq_zero (u v w q : R) (hq : q ∈ I) (huvw : u * v * w = q)
    (hvw : v * w ∈ I) (huw : u * w ∈ I) (huv : u * v ∈ I)
    (hu : algebraMap R K (1 - u) ≠ 0) (hv : algebraMap R K (1 - v) ≠ 0)
    (hw : algebraMap R K (1 - w) ≠ 0)
    (h12 : tateXK u (v * w) q hq ≠ (tateXK v (u * w) q hq : K))
    (h13 : tateXK u (v * w) q hq ≠ (tateXK w (u * v) q hq : K))
    (h23 : tateXK v (u * w) q hq ≠ (tateXK w (u * v) q hq : K))
    (hΔ : ((tateCurveAt q hq).map (algebraMap R K)).toAffine.Δ ≠ 0) :
    tatePtPair u (v * w) q hq (by rw [← huvw]; ring) (isUnit_one_sub hvw) hu hΔ
      + tatePtPair v (u * w) q hq (by rw [← huvw]; ring) (isUnit_one_sub huw) hv hΔ
      + tatePtPair w (u * v) q hq (by rw [← huvw]; ring) (isUnit_one_sub huv) hw hΔ = 0 :=
  tate_points_add_eq_zero_K u v w q hq huvw hvw huw huv hu hv hw h12 h13 h23 _ _ _

set_option maxHeartbeats 1600000 in
/-- ★★★★★★**相方に入れ替えると反転する**(環帯・単元の領域)。 -/
theorem tatePtPair_swap (a w q : R) (hq : q ∈ I) (haw : a * w = q)
    (ha : IsUnit (1 - a)) (hw : IsUnit (1 - w))
    (hΔ : ((tateCurveAt q hq).map (algebraMap R K)).toAffine.Δ ≠ 0) :
    tatePtPair w a q hq (by rw [← haw]; ring) ha ((hw.map (algebraMap R K)).ne_zero) hΔ
      = -tatePtPair a w q hq haw hw ((ha.map (algebraMap R K)).ne_zero) hΔ := by
  rw [tatePtPair, tatePtPair, Point.neg_some]
  have e1 : ((tateCurveAt q hq).map (algebraMap R K)).toAffine.a₁ = 1 := by
    show algebraMap R K ((tateCurveAt q hq).a₁) = 1
    rw [show (tateCurveAt q hq).a₁ = 1 by
      simp [tateCurveAt, tateCurve, WeierstrassCurve.map], map_one]
  have e3 : ((tateCurveAt q hq).map (algebraMap R K)).toAffine.a₃ = 0 := by
    show algebraMap R K ((tateCurveAt q hq).a₃) = 0
    rw [show (tateCurveAt q hq).a₃ = 0 by
      simp [tateCurveAt, tateCurve, WeierstrassCurve.map], map_zero]
  simp only [Point.some.injEq]
  constructor
  · rw [tateXK_eq w a q hq hw, tateXK_eq a w q hq ha, tateXpair_symm]
  · rw [WeierstrassCurve.Affine.negY, e1, e3, tateYK_eq w a q hq hw, tateXK_eq a w q hq ha,
      tateYK_eq a w q hq ha, tateYpair_swap]
    simp only [map_sub, map_neg]
    ring

/-- ★★★★★★★★**`P(1/u) = −P(u)`**(対の形、全領域)。 -/
theorem tatePtPair_ringInverse (a q : R) (hq : q ∈ I) (ha : IsUnit a)
    (hne : algebraMap R K (1 - a) ≠ 0)
    (hΔ : ((tateCurveAt q hq).map (algebraMap R K)).toAffine.Δ ≠ 0) :
    tatePtPair (Ring.inverse a) (q * a) q hq
        (by rw [show Ring.inverse a * (q * a) = q * (Ring.inverse a * a) by ring,
          Ring.inverse_mul_cancel _ ha, mul_one])
        (isUnit_one_sub (Ideal.mul_mem_right _ _ hq))
        (map_one_sub_ringInverse_ne_zero ha hne) hΔ
      = -tatePtPair a (wOf q a) q hq (mul_wOf q a ha) (isUnit_one_sub (wOf_mem a hq)) hne hΔ :=
  tate_point_ringInverse_eq_neg a q hq ha hne hΔ

/-- ★★母数から来る点はけっして原点ではない。 -/
theorem tatePtPair_ne_zero (a w q : R) (hq : q ∈ I) (haw : a * w = q)
    (hw : IsUnit (1 - w)) (hne : algebraMap R K (1 - a) ≠ 0)
    (hΔ : ((tateCurveAt q hq).map (algebraMap R K)).toAffine.Δ ≠ 0) :
    tatePtPair a w q hq haw hw hne hΔ ≠ 0 := by
  rw [tatePtPair]
  simp

/-! ## ★出典の紐付け(`.src`) -/

def tatePtPair.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——3 領域を通じた一つの点の式)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
