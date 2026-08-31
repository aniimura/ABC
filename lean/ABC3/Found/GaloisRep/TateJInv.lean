/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.TateJ
import ABC3.Found.GaloisRep.TateMultRed
import ABC3.Found.GaloisRep.TateLinearQ

/-!
# Galois (G6) 第 877 ブロック —— **★★★★★★★★`1/j = q + O(q²)` と単射性**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> parameter qE of E satisfies the relation qE = qEl ; in particular, we have

## ★★★★★★★★★★これは何か

`tateParam_quot_mu`（第 873 の葉）に残る義務のひとつが「`q ↦ j(E_q)` の単射性」である。
在庫を見ると必要なものはほぼ揃っていた:

* `tateCurve_Delta_eq` / `coeff_zero_tateDelta` / `coeff_one_tateDelta`（`Δ = q + O(q²)`）
* `tateC4_isUnit` / `constantCoeff_tateC4`（`c₄ = 1 + O(q)`、単元）
* `tateInvJ_eq_X_mul_unit`（`Δ = X·w·c₄³`、`w` は単元）
* `evalAdic_injective_of_coeff_one`（第 875）

★本ブロックはそれらを繋いで **`q ↦ (1/j)(q)` は `I` の上で単射**を出す。

| 定義・定理 | 内容 |
|---|---|
| `tateJinvSeries` | `Δ · (c₄³)⁻¹` |
| `tateJinvSeries_mul_c4` | ★`(1/j)·c₄³ = Δ` |
| `coeff_one_tateJinvSeries` | ★★★**1 次の係数は `1`** |
| `tateJinv_injective` | ★★★★★★**単射性** |
-/

namespace ABC3.Found.GaloisRep

variable {R : Type} [CommRing R] {I : Ideal R}

theorem coeff_one_mul' (f g : PowerSeries ℤ) :
    PowerSeries.coeff 1 (f * g)
      = PowerSeries.constantCoeff f * PowerSeries.coeff 1 g
        + PowerSeries.coeff 1 f * PowerSeries.constantCoeff g := by
  rw [PowerSeries.coeff_mul, Finset.Nat.sum_antidiagonal_eq_sum_range_succ_mk,
    Finset.sum_range_succ, Finset.sum_range_succ, Finset.sum_range_zero]
  simp

/-! ## ★★★★★★★★`1/j` のべき級数 -/

/-- ★`c₄³` を単元として名前を付ける（参照を安定させるため）。 -/
noncomputable def tateC4CubeUnit : (PowerSeries ℤ)ˣ := (tateC4_isUnit.pow 3).unit

theorem tateC4CubeUnit_val : (tateC4CubeUnit : PowerSeries ℤ) = tateCurve.c₄ ^ 3 :=
  (tateC4_isUnit.pow 3).unit_spec

/-- ★★★★**`1/j = Δ · (c₄³)⁻¹`**。`c₄` は単元なので逆が取れる。 -/
noncomputable def tateJinvSeries : PowerSeries ℤ :=
  tateCurve.Δ * ((tateC4CubeUnit⁻¹ : (PowerSeries ℤ)ˣ) : PowerSeries ℤ)

theorem c4inv_mul_c4 :
    ((tateC4CubeUnit⁻¹ : (PowerSeries ℤ)ˣ) : PowerSeries ℤ) * tateCurve.c₄ ^ 3 = 1 := by
  rw [← tateC4CubeUnit_val, ← Units.val_mul, inv_mul_cancel, Units.val_one]

/-- ★★★**`(1/j)·c₄³ = Δ`**。 -/
theorem tateJinvSeries_mul_c4 :
    tateJinvSeries * tateCurve.c₄ ^ 3 = tateCurve.Δ := by
  rw [tateJinvSeries, mul_assoc, c4inv_mul_c4, mul_one]

theorem constantCoeff_c4inv :
    PowerSeries.constantCoeff ((tateC4CubeUnit⁻¹ : (PowerSeries ℤ)ˣ) : PowerSeries ℤ) = 1 := by
  have h := congrArg PowerSeries.constantCoeff c4inv_mul_c4
  rw [map_mul, map_pow, constantCoeff_tateC4, map_one] at h
  simpa using h

theorem constantCoeff_tateJinvSeries :
    PowerSeries.constantCoeff tateJinvSeries = 0 := by
  rw [tateJinvSeries, map_mul, constantCoeff_c4inv,
    ← PowerSeries.coeff_zero_eq_constantCoeff_apply, coeff_zero_tateDelta]
  ring

/-- ★★★★★★**`1/j` の 1 次の係数は `1`**——`Δ = q + O(q²)` と `c₄ = 1 + O(q)` から。 -/
theorem coeff_one_tateJinvSeries : PowerSeries.coeff 1 tateJinvSeries = 1 := by
  rw [tateJinvSeries, coeff_one_mul', constantCoeff_c4inv,
    ← PowerSeries.coeff_zero_eq_constantCoeff_apply, coeff_zero_tateDelta,
    coeff_one_tateDelta]
  ring

/-! ## ★★★★★★★★単射性 -/

/-- ★★★★★★★★**`q ↦ (1/j)(q)` は `I` の上で単射**。

★これが `Lemma 3.5` の残る義務「`q ↦ j(E_q)` の単射性」の中身である。 -/
theorem tateJinv_injective [IsAdicComplete I R] {a b : R} (ha : a ∈ I) (hb : b ∈ I)
    (h : evalAdic tateJinvSeries a ha = evalAdic tateJinvSeries b hb) : a = b := by
  refine evalAdic_injective_of_coeff_one tateJinvSeries ?_ coeff_one_tateJinvSeries ha hb h
  rw [PowerSeries.coeff_zero_eq_constantCoeff_apply]
  exact constantCoeff_tateJinvSeries

/-- ★★★★**`(1/j)(q)·c₄(E_q)³ = Δ(E_q)`**——べき級数の恒等式を評価しただけ。 -/
theorem evalAdic_tateJinvSeries_mul_c4 [IsAdicComplete I R] (q : R) (hq : q ∈ I) :
    evalAdic tateJinvSeries q hq * (tateCurveAt q hq).c₄ ^ 3 = (tateCurveAt q hq).Δ := by
  rw [tateCurveAt_c4, tateCurveAt_Delta]
  show evalAdicHom q hq tateJinvSeries * _ = _
  rw [← map_pow, ← map_mul, tateJinvSeries_mul_c4]

/-- ★★★★★★★★★★**Tate 母数は `j` で決まる**——分母を払った形。

`j = c₄³/Δ` なので `j(E_q) = j(E_{q′})` は
`Δ(E_q)·c₄(E_{q′})³ = Δ(E_{q′})·c₄(E_q)³` と同じことである。

★これが `Lemma 3.5` の残る義務「`q ↦ j(E_q)` の単射性」そのものである。 -/
theorem tateParam_injective [IsAdicComplete I R] {q q' : R} (hq : q ∈ I) (hq' : q' ∈ I)
    (h : (tateCurveAt q hq).Δ * (tateCurveAt q' hq').c₄ ^ 3
       = (tateCurveAt q' hq').Δ * (tateCurveAt q hq).c₄ ^ 3) : q = q' := by
  have e1 := evalAdic_tateJinvSeries_mul_c4 q hq
  have e2 := evalAdic_tateJinvSeries_mul_c4 q' hq'
  have hu : IsUnit ((tateCurveAt q hq).c₄ ^ 3 * (tateCurveAt q' hq').c₄ ^ 3) :=
    ((tateCurveAt_c4_isUnit q hq).pow 3).mul ((tateCurveAt_c4_isUnit q' hq').pow 3)
  refine tateJinv_injective hq hq' ?_
  refine hu.mul_left_cancel ?_
  have hlhs : ((tateCurveAt q hq).c₄ ^ 3 * (tateCurveAt q' hq').c₄ ^ 3)
      * evalAdic tateJinvSeries q hq
      = (tateCurveAt q hq).Δ * (tateCurveAt q' hq').c₄ ^ 3 := by
    rw [← e1]
    ring
  have hrhs : ((tateCurveAt q hq).c₄ ^ 3 * (tateCurveAt q' hq').c₄ ^ 3)
      * evalAdic tateJinvSeries q' hq'
      = (tateCurveAt q' hq').Δ * (tateCurveAt q hq).c₄ ^ 3 := by
    rw [← e2]
    ring
  rw [hlhs, hrhs]
  exact h

/-! ## ★出典の紐付け(`.src`) -/

def tateJinvSeries.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(1/j = Δ·(c₄³)⁻¹ のべき級数。★無条件)",
    sectionId := "genell-lemma-3-2" }

def coeff_one_tateJinvSeries.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(1/j の 1 次の係数は 1。★無条件)",
    sectionId := "genell-lemma-3-2" }

def tateJinv_injective.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(q ↦ (1/j)(q) は単射。★無条件)",
    sectionId := "genell-lemma-3-2" }

end ABC3.Found.GaloisRep
