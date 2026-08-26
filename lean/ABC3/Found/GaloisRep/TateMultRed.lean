import ABC3.Found.GaloisRep.TateVarChange
import Mathlib.AlgebraicGeometry.EllipticCurve.Reduction

/-!
# Galois (G6) 第 304 ブロック —— **★★★★★★★★乗法還元から Tate 母数が取れる**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★★★到達点

> `W` が乗法還元をもち `Δ ≠ 0` なら、`q ∈ 𝔪`、`q ≠ 0` が在って
> **`j(E_q) = j(W)`**(`exists_tateParam_of_mult`)

★★★第 100 の `exists_tateParam`(`j` の反転)を、**mathlib の還元の言葉に接続**した。

## ★★★★★★mathlib の乗法還元は 2 つの不等式

    v(Δ) < 1        (`badReduction`)
    v(c₄) = 1       (`multiplicativeReduction`)

★どちらも `IsDiscreteValuationRing.maximalIdeal R` の adic 付値である。
★★第 302 で `tateDvrVal` を**その付値に揃えた**ので、第 301 の 2 条件
(`dvr_mem_of_nonneg`・`dvr_mem_max_of_pos`)がそのまま使える:

| mathlib の条件 | `R` の言葉 |
|---|---|
| `v(c₄) = 1` | `c₄` は `R` の**単元**から来る |
| `v(Δ) < 1` | `Δ` は **`𝔪` の元**から来る |

## ★★★★★★★`Δ ≠ 0` は仮定に要る——結節三次曲線の反例

★★★**mathlib の `HasMultiplicativeReduction` は `Δ = 0` を排除しない**
(`v(0) = 0 < 1` だから)。実際 `y² + xy = x³`(`⟨1,0,0,0,0⟩`)は

* `Δ = 0`、`c₄ = 1`(単元)
* 接線の 2 次式は `X² + X` で**分解する**(分裂乗法還元の条件を満たす)
* `Δ` が全モデルで `0` なので `IsMinimal` も満たす

★★★★これは**結節三次曲線**であり、滑らかな点の群は `Kˣ` そのもの——
Tate 母数は `q = 1` にあたるので**局所高さが `0`** になる。
★したがって界面 `TateCurveData` は `Δ ≠ 0`(`IsElliptic`)を要求しなければならない。
★★本ブロックでその訂正も行った(§9-617)。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `exists_isUnit_of_valuation_eq_one` | ★★★★付値 `1` は `R` の単元 |
| `exists_mem_of_valuation_lt_one` | ★★★★付値 `< 1` は `𝔪` の元 |
| `tateCurveAt_c4`・`tateCurveAt_c4_isUnit` | ★★`E_q` の `c₄` は単元 |
| `exists_tateParam_of_mult` | ★★★★★★★★**乗法還元 ⟹ Tate 母数** |
-/

namespace ABC3.Found.GaloisRep

open IsDedekindDomain WeierstrassCurve

section Val

variable {R : Type} [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
  {K : Type} [Field K] [Algebra R K] [IsFractionRing R K]

/-- ★★★★**付値が `1` なら `R` の単元から来る**。 -/
theorem exists_isUnit_of_valuation_eq_one (x : K)
    (h : (IsDiscreteValuationRing.maximalIdeal R).valuation K x = 1) :
    ∃ r : R, IsUnit r ∧ algebraMap R K r = x := by
  have hx0 : x ≠ 0 := by
    intro hc
    rw [hc, map_zero] at h
    exact zero_ne_one h
  set xu : Kˣ := Units.mk0 x hx0 with hxu
  have hval : (IsDiscreteValuationRing.maximalIdeal R).valuation K ((xu : K)) = 1 := h
  have h1 : 0 ≤ vAdd (tateDvrVal R K) xu := (tateDvrVal_nonneg_iff xu).2 (le_of_eq hval)
  have hinv : (IsDiscreteValuationRing.maximalIdeal R).valuation K ((xu⁻¹ : Kˣ) : K) = 1 := by
    rw [Units.val_inv_eq_inv_val, map_inv₀, hval, inv_one]
  have h2 : 0 ≤ vAdd (tateDvrVal R K) xu⁻¹ := (tateDvrVal_nonneg_iff xu⁻¹).2 (le_of_eq hinv)
  rw [vAdd_inv] at h2
  obtain ⟨r, hr⟩ := dvr_mem_of_nonneg xu h1
  refine ⟨r, ?_, hr⟩
  by_contra hu
  have hm : r ∈ IsLocalRing.maximalIdeal R := (IsLocalRing.mem_maximalIdeal r).2 hu
  have := tateDvrVal_pos_of_mem_max xu ⟨r, hm, hr⟩
  omega

/-- ★★★★**付値が `1` 未満なら `𝔪` の元から来る**。 -/
theorem exists_mem_of_valuation_lt_one (x : K) (hx0 : x ≠ 0)
    (h : (IsDiscreteValuationRing.maximalIdeal R).valuation K x < 1) :
    ∃ r ∈ IsLocalRing.maximalIdeal R, algebraMap R K r = x := by
  set xu : Kˣ := Units.mk0 x hx0 with hxu
  have hpos : 0 < vAdd (tateDvrVal R K) xu := (tateDvrVal_pos_iff xu).2 h
  exact dvr_mem_max_of_pos xu hpos

end Val

section C4

variable {R : Type} [CommRing R] {I : Ideal R} [IsAdicComplete I R]

/-- ★`E_q` の `c₄` は特殊化と両立する。 -/
theorem tateCurveAt_c4 (q : R) (hq : q ∈ I) :
    (tateCurveAt q hq).c₄ = evalAdicHom q hq tateCurve.c₄ := by
  rw [tateCurveAt, WeierstrassCurve.map_c₄]

/-- ★★`E_q` の `c₄` は単元(定数項が `1` だから)。 -/
theorem tateCurveAt_c4_isUnit (q : R) (hq : q ∈ I) :
    IsUnit (tateCurveAt q hq).c₄ := by
  rw [tateCurveAt_c4]
  exact tateC4_isUnit.map (evalAdicHom q hq)

end C4

section MultRed

variable {R : Type} [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
  [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
  {K : Type} [Field K] [Algebra R K] [IsFractionRing R K]

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★★**乗法還元をもつ楕円曲線には Tate 母数がある**。

結論は分母を払った形で書いた——`Δ_q · c₄(W)³ = Δ(W) · c₄(E_q)³` は
**`j(E_q) = j(W)`** に他ならない(`j = c₄³/Δ`)。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem exists_tateParam_of_mult (W : WeierstrassCurve K)
    [hmul : W.HasMultiplicativeReduction R] (hΔ0 : W.Δ ≠ 0) :
    ∃ (q : R) (hq : q ∈ IsLocalRing.maximalIdeal R), q ≠ 0 ∧
      algebraMap R K (tateCurveAt q hq).Δ * W.c₄ ^ 3
        = W.Δ * algebraMap R K (tateCurveAt q hq).c₄ ^ 3 := by
  obtain ⟨c, hcu, hc⟩ := exists_isUnit_of_valuation_eq_one W.c₄ hmul.multiplicativeReduction
  obtain ⟨d, hdm, hd⟩ := exists_mem_of_valuation_lt_one W.Δ hΔ0 hmul.badReduction
  have hd0 : d ≠ 0 := by
    intro h
    rw [h, map_zero] at hd
    exact hΔ0 hd.symm
  have hunit : algebraMap R K ((hcu.unit⁻¹ : Rˣ) : R) * W.c₄ = 1 := by
    rw [← hc, ← map_mul]
    have h2 : ((hcu.unit⁻¹ : Rˣ) : R) * ((hcu.unit : Rˣ) : R) = 1 := by
      rw [← Units.val_mul, inv_mul_cancel]
      rfl
    rw [hcu.unit_spec] at h2
    rw [h2, map_one]
  set t : R := d * ((hcu.unit⁻¹ : Rˣ) : R) ^ 3 with ht
  have htm : t ∈ IsLocalRing.maximalIdeal R := Ideal.mul_mem_right _ _ hdm
  have ht0 : t ≠ 0 := mul_ne_zero hd0 (pow_ne_zero _ (Units.ne_zero _))
  obtain ⟨q, hq, hqe⟩ := exists_tateParam (I := IsLocalRing.maximalIdeal R) htm
  have hΔq : (tateCurveAt q hq).Δ = t * (tateCurveAt q hq).c₄ ^ 3 := by
    rw [tateCurveAt_Delta, tateCurveAt_c4, ← map_pow]
    exact hqe
  refine ⟨q, hq, ?_, ?_⟩
  · intro hq0
    obtain ⟨v, hv, hΔeq⟩ := tateCurveAt_Delta_eq_mul_unit q hq
    have hzero : (tateCurveAt q hq).Δ = 0 := by
      rw [hΔeq, hq0, zero_mul]
    rw [hzero] at hΔq
    have hc4u : IsUnit ((tateCurveAt q hq).c₄ ^ 3) := (tateCurveAt_c4_isUnit q hq).pow 3
    rcases mul_eq_zero.1 hΔq.symm with h | h
    · exact ht0 h
    · rw [h] at hc4u
      exact not_isUnit_zero hc4u
  · have key := congrArg (algebraMap R K) hΔq
    rw [map_mul, map_pow, ht, map_mul, map_pow, hd] at key
    linear_combination (W.c₄ ^ 3) * key
      + W.Δ * (algebraMap R K (tateCurveAt q hq).c₄) ^ 3
        * ((algebraMap R K ((hcu.unit⁻¹ : Rˣ) : R) * W.c₄) ^ 2
          + (algebraMap R K ((hcu.unit⁻¹ : Rˣ) : R) * W.c₄) + 1) * hunit

end MultRed

/-! ## ★出典の紐付け(`.src`) -/

def exists_tateParam_of_mult.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(乗法還元から Tate 母数を取る)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
