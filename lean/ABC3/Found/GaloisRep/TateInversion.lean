import ABC3.Found.GaloisRep.TateSeparate

/-!
# Galois (G6) 第 281 ブロック —— **★★★★★★★★全領域で成り立つ反転則**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★★★到達点

> `u` が単元なら **`P(1/u) = −P(u)`**(`tate_point_ringInverse_eq_neg`)

★★★`1 − u` が単元である必要が**無い**——原点近傍でも成り立つ。
準同型性の退化した場合(`c₂ = c₁⁻¹`)と倍化の補助母数の議論に要る。

## ★★★★★★分母を払えば反転則も通る

`f(1/t) = f(t)`・`g(1/t) = −g(t) − f(t)` は `1 − t` が単元でないと書けない
(`Ring.inverse (1−t)` が意味を持たない)。★★しかし**分母を払った形なら通る**:

    XE(1/u, qu) = XE(u, q/u)·(1/u)²
    YE(1/u, qu) = (1/u)³·((1−u)·XE(u, q/u) + YE(u, q/u))

★★★`1 − 1/u = −(1−u)/u` なので `(1−1/u)² = (1−u)²/u²`——**分母の変換も `u` の冪だけ**。
これを `K` で割れば `X(1/u) = X(u)`、`Y(1/u) = −X(u) − Y(u)` になる。

## ★★★★★尾のずらしが両側をつなぐ

    T(u) = f(qu) + T(qu)        (`tateXtail_shift`)

★これで `X(1/u, qu)` の `f(qu) + T(qu)` が `X(u, q/u)` の `T(u)` に化ける。
★★`adicSum_shift`(在庫)がそのまま使えた。

## ★在庫の引き当て(2 件)

- `tateXterm_inv`・`tateYterm_inv` は **`TateXY.lean` に既にあった**(体の `⁻¹` 版)。
  本ブロックのものは一般の環の `Ring.inverse` 版なので `_ringInverse` と改名した。
- `isUnit_ringInverse` は **mathlib にある**(`IsUnit (Ring.inverse a) ↔ IsUnit a`)。
  自作せず `.mpr` を使った。

★着手前の `grep` は**補題名の単位で**行うこと(第 270 と同じ轍)。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `tateXterm_ringInverse`・`tateYterm_ringInverse` | ★★★★★項の反転則(`Ring.inverse` 版) |
| `tateXtail_shift`・`tateYtail_shift` | ★★★★★尾のずらし |
| `tateXpairE_ringInverse`・`tateYpairE_ringInverse` | ★★★★★★分母を払った反転則 |
| `tateXK_ringInverse`・`tateYK_ringInverse` | ★★★★★★★`K` の水準の反転則 |
| `tate_point_swap_eq_neg` | ★★★★★★相方への入れ替えは反転 |
| `tate_point_ringInverse_eq_neg` | ★★★★★★★★**全領域で `P(1/u) = −P(u)`** |
-/

namespace ABC3.Found.GaloisRep

open Complex Real WeierstrassCurve WeierstrassCurve.Affine

/-! ## ★★★★★項の反転則 -/

section Term

variable {R : Type} [CommRing R]

theorem one_sub_ringInverse_eq {t : R} (ht : IsUnit t) :
    (1 : R) - Ring.inverse t = -((1 - t) * Ring.inverse t) := by
  have hti : Ring.inverse t * t = 1 := Ring.inverse_mul_cancel _ ht
  linear_combination -hti

theorem isUnit_one_sub_ringInverse {t : R} (ht : IsUnit t) (h1 : IsUnit (1 - t)) :
    IsUnit (1 - Ring.inverse t) := by
  rw [one_sub_ringInverse_eq ht]
  exact (h1.mul (isUnit_ringInverse.mpr ht)).neg

/-- ★★★★★**`f(1/t) = f(t)`**(`Ring.inverse` 版)——`X` は `u ↦ 1/u` で不変。 -/
theorem tateXterm_ringInverse {t : R} (ht : IsUnit t) (h1 : IsUnit (1 - t)) :
    tateXterm (Ring.inverse t) = tateXterm t := by
  have hti : Ring.inverse t * t = 1 := Ring.inverse_mul_cancel _ ht
  have hinv := isUnit_one_sub_ringInverse ht h1
  have h2 := mul_tateXterm' h1
  have h3 : Ring.inverse t ^ 2 * t = Ring.inverse t := by
    rw [pow_succ, mul_assoc, hti, mul_one, pow_one]
  have hkey : (1 - Ring.inverse t) ^ 2 * tateXterm t = Ring.inverse t := by
    rw [one_sub_ringInverse_eq ht]
    calc (-((1 - t) * Ring.inverse t)) ^ 2 * tateXterm t
        = Ring.inverse t ^ 2 * ((1 - t) ^ 2 * tateXterm t) := by ring
      _ = Ring.inverse t ^ 2 * t := by rw [h2]
      _ = Ring.inverse t := h3
  have h4 := mul_tateXterm' hinv
  exact (hinv.pow 2).mul_left_cancel (by rw [h4, hkey])

/-- ★★★★★**`g(1/t) = −g(t) − f(t)`**(`Ring.inverse` 版)。 -/
theorem tateYterm_ringInverse {t : R} (ht : IsUnit t) (h1 : IsUnit (1 - t)) :
    tateYterm (Ring.inverse t) = -tateYterm t - tateXterm t := by
  have hti : Ring.inverse t * t = 1 := Ring.inverse_mul_cancel _ ht
  have hinv := isUnit_one_sub_ringInverse ht h1
  have hX := mul_tateXterm' h1
  have hY := mul_tateYterm' h1
  have h3 : Ring.inverse t ^ 3 * t = Ring.inverse t ^ 2 := by
    rw [pow_succ, mul_assoc, hti, mul_one]
  have hkey : (1 - Ring.inverse t) ^ 3 * (-tateYterm t - tateXterm t)
      = Ring.inverse t ^ 2 := by
    rw [one_sub_ringInverse_eq ht]
    calc (-((1 - t) * Ring.inverse t)) ^ 3 * (-tateYterm t - tateXterm t)
        = Ring.inverse t ^ 3 * ((1 - t) ^ 3 * tateYterm t)
          + Ring.inverse t ^ 3 * ((1 - t) * ((1 - t) ^ 2 * tateXterm t)) := by ring
      _ = Ring.inverse t ^ 3 * t ^ 2 + Ring.inverse t ^ 3 * ((1 - t) * t) := by rw [hY, hX]
      _ = Ring.inverse t ^ 3 * t := by ring
      _ = Ring.inverse t ^ 2 := h3
  have h4 := mul_tateYterm' hinv
  exact (hinv.pow 3).mul_left_cancel (by rw [h4, hkey])

end Term

/-! ## ★★★★★尾のずらし -/

section Shift

variable {R : Type} [CommRing R] {I : Ideal R}

/-- ★★★★★**尾のずらし**——`T(u) = f(qu) + T(qu)`。 -/
theorem tateXtail_shift [IsAdicComplete I R] (u q : R) (hq : q ∈ I) :
    tateXtail u q hq = tateXterm (q * u) + tateXtail (q * u) q hq := by
  rw [tateXtail, adicSum_shift]
  congr 1
  · norm_num
  · rw [tateXtail]
    refine adicSum_congr _ _ fun n => ?_
    congr 1
    ring

theorem tateYtail_shift [IsAdicComplete I R] (u q : R) (hq : q ∈ I) :
    tateYtail u q hq = tateYterm (q * u) + tateYtail (q * u) q hq := by
  rw [tateYtail, adicSum_shift]
  congr 1
  · norm_num
  · rw [tateYtail]
    refine adicSum_congr _ _ fun n => ?_
    congr 1
    ring

/-! ## ★★★★★★分母を払った反転則 -/

/-- ★★★★★★**`X` は `u ↦ 1/u` で不変**(分母を払った形)。 -/
theorem tateXpairE_ringInverse [IsAdicComplete I R] (a q : R) (hq : q ∈ I) (ha : IsUnit a) :
    tateXpairE (Ring.inverse a) (q * a) q hq
      = tateXpairE a (wOf q a) q hq * Ring.inverse a ^ 2 := by
  have hti : Ring.inverse a * a = 1 := Ring.inverse_mul_cancel _ ha
  rw [tateXpairE, tateXpairE, wOf, one_sub_ringInverse_eq ha,
    tateXtail_shift a q hq, tateXtail_shift (Ring.inverse a) q hq]
  linear_combination (-Ring.inverse a) * hti

/-- ★★★★★★**`Y` は `u ↦ 1/u` で反転する**(分母を払った形)。 -/
theorem tateYpairE_ringInverse [IsAdicComplete I R] (a q : R) (hq : q ∈ I) (ha : IsUnit a) :
    tateYpairE (Ring.inverse a) (q * a) q hq
      = Ring.inverse a ^ 3
        * ((1 - a) * tateXpairE a (wOf q a) q hq + tateYpairE a (wOf q a) q hq) := by
  have hti : Ring.inverse a * a = 1 := Ring.inverse_mul_cancel _ ha
  rw [tateYpairE, tateYpairE, tateXpairE, wOf, one_sub_ringInverse_eq ha,
    tateXtail_shift a q hq, tateYtail_shift a q hq, tateYtail_shift (Ring.inverse a) q hq]
  linear_combination (-Ring.inverse a ^ 2) * hti

end Shift

/-! ## ★★★★★★★`K` の水準の反転則 -/

section FieldLevel

variable {R K : Type} [CommRing R] {I : Ideal R} [IsAdicComplete I R] [Field K] [Algebra R K]

/-- ★★★★★★★**`X(1/u) = X(u)`**——`K` の水準、全領域で成り立つ。 -/
theorem tateXK_ringInverse (a q : R) (hq : q ∈ I) (ha : IsUnit a)
    (hne : algebraMap R K (1 - a) ≠ 0) :
    tateXK (Ring.inverse a) (q * a) q hq = (tateXK a (wOf q a) q hq : K) := by
  have hai : algebraMap R K (Ring.inverse a) ≠ 0 :=
    ((isUnit_ringInverse.mpr ha).map (algebraMap R K)).ne_zero
  rw [tateXK, tateXK, tateXpairE_ringInverse a q hq ha, one_sub_ringInverse_eq ha]
  simp only [map_mul, map_pow, map_neg, map_sub, map_one]
  field_simp

/-- ★★★★★★★**`Y(1/u) = −X(u) − Y(u)`**——`K` の水準、全領域で成り立つ。 -/
theorem tateYK_ringInverse (a q : R) (hq : q ∈ I) (ha : IsUnit a)
    (hne : algebraMap R K (1 - a) ≠ 0) :
    tateYK (Ring.inverse a) (q * a) q hq
      = -tateXK a (wOf q a) q hq - (tateYK a (wOf q a) q hq : K) := by
  have hai : algebraMap R K (Ring.inverse a) ≠ 0 :=
    ((isUnit_ringInverse.mpr ha).map (algebraMap R K)).ne_zero
  rw [tateYK, tateXK, tateYK, tateYpairE_ringInverse a q hq ha, one_sub_ringInverse_eq ha]
  simp only [map_mul, map_pow, map_neg, map_sub, map_one, map_add]
  field_simp
  ring

theorem map_one_sub_ringInverse_ne_zero {a : R} (ha : IsUnit a)
    (hne : algebraMap R K (1 - a) ≠ 0) : algebraMap R K (1 - Ring.inverse a) ≠ 0 := by
  rw [one_sub_ringInverse_eq ha]
  simp only [map_neg, map_mul, ne_eq, neg_eq_zero, mul_eq_zero, not_or]
  exact ⟨hne, ((isUnit_ringInverse.mpr ha).map (algebraMap R K)).ne_zero⟩

end FieldLevel

/-! ## ★★★★★★★★点の水準の反転則 -/

section Point

variable {R K : Type} [CommRing R] {I : Ideal R} [IsAdicComplete I R] [Field K]
  [DecidableEq K] [Algebra R K]

/-- ★★★★★★**相方に入れ替えると点は反転する**——`(X(q/u), Y(q/u)) = −(X(u), Y(u))`。 -/
theorem tate_point_swap_eq_neg (a w q : R) (hq : q ∈ I) (haw : a * w = q)
    (ha : IsUnit (1 - a)) (hw : IsUnit (1 - w))
    (hΔ : ((tateCurveAt q hq).map (algebraMap R K)).toAffine.Δ ≠ 0) :
    Point.some _ _ (nonsingular_tate_point' (K := K) w a q hq (by rw [← haw]; ring) hw ha hΔ)
      = -Point.some _ _ (nonsingular_tate_point' (K := K) a w q hq haw ha hw hΔ) := by
  rw [Point.neg_some]
  have e1 : ((tateCurveAt q hq).map (algebraMap R K)).toAffine.a₁ = 1 := by
    show algebraMap R K ((tateCurveAt q hq).a₁) = 1
    rw [show (tateCurveAt q hq).a₁ = 1 by
      simp [tateCurveAt, tateCurve, WeierstrassCurve.map], map_one]
  have e3 : ((tateCurveAt q hq).map (algebraMap R K)).toAffine.a₃ = 0 := by
    show algebraMap R K ((tateCurveAt q hq).a₃) = 0
    rw [show (tateCurveAt q hq).a₃ = 0 by
      simp [tateCurveAt, tateCurve, WeierstrassCurve.map], map_zero]
  have hX : algebraMap R K (tateXpair w a q hq) = algebraMap R K (tateXpair a w q hq) := by
    rw [tateXpair_symm]
  have hY : algebraMap R K (tateYpair w a q hq)
      = ((tateCurveAt q hq).map (algebraMap R K)).toAffine.negY
        (algebraMap R K (tateXpair a w q hq)) (algebraMap R K (tateYpair a w q hq)) := by
    rw [WeierstrassCurve.Affine.negY, e1, e3, tateYpair_swap]
    simp only [map_sub, map_neg]
    ring
  simp only [Point.some.injEq]
  exact ⟨hX, hY⟩

/-- ★★★★★★★★**`P(1/u) = −P(u)`**——全領域で成り立つ反転則。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem tate_point_ringInverse_eq_neg (a q : R) (hq : q ∈ I) (ha : IsUnit a)
    (hne : algebraMap R K (1 - a) ≠ 0)
    (hΔ : ((tateCurveAt q hq).map (algebraMap R K)).toAffine.Δ ≠ 0) :
    Point.some _ _ (nonsingular_tateK (K := K) (Ring.inverse a) (q * a) q hq
        (by rw [show Ring.inverse a * (q * a) = q * (Ring.inverse a * a) by ring,
          Ring.inverse_mul_cancel _ ha, mul_one])
        (isUnit_one_sub (Ideal.mul_mem_right _ _ hq))
        (map_one_sub_ringInverse_ne_zero ha hne) hΔ)
      = -Point.some _ _ (nonsingular_tateK (K := K) a (wOf q a) q hq (mul_wOf q a ha)
        (isUnit_one_sub (wOf_mem a hq)) hne hΔ) := by
  rw [Point.neg_some]
  have e1 : ((tateCurveAt q hq).map (algebraMap R K)).toAffine.a₁ = 1 := by
    show algebraMap R K ((tateCurveAt q hq).a₁) = 1
    rw [show (tateCurveAt q hq).a₁ = 1 by
      simp [tateCurveAt, tateCurve, WeierstrassCurve.map], map_one]
  have e3 : ((tateCurveAt q hq).map (algebraMap R K)).toAffine.a₃ = 0 := by
    show algebraMap R K ((tateCurveAt q hq).a₃) = 0
    rw [show (tateCurveAt q hq).a₃ = 0 by
      simp [tateCurveAt, tateCurve, WeierstrassCurve.map], map_zero]
  simp only [Point.some.injEq]
  refine ⟨tateXK_ringInverse a q hq ha hne, ?_⟩
  rw [WeierstrassCurve.Affine.negY, e1, e3, tateYK_ringInverse a q hq ha hne]
  ring

end Point

/-! ## ★出典の紐付け(`.src`) -/

def tateXtail_shift.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——尾のずらし)",
    sectionId := "genell-def-3-3" }

def tate_point_ringInverse_eq_neg.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——全領域で P(1/u) = -P(u))",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
