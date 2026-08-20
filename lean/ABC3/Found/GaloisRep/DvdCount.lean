import ABC3.Found.GaloisRep.HyperInvValuation
import ABC3.Found.GaloisRep.PullbackAlg

/-!
# Galois (G5) 第 149 ブロック —— **★★★★★★★D1' が閉じた**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★`div(f_P ∘ [n])` が `n` で割れる

    ∀ v,  n ∣ count_v((μ f_P))

★2 つの場合を合流させて閉じた:

| 場合 | 判定 | 内容 |
|---|---|---|
| A | `∀ r, w(μ r) ≤ 1` | 第 143(水準イデアルへの帰納法) |
| B | `∃ r, 1 < w(μ r)` | 第 144–148(無限遠、偶奇、対合) |

★★**分岐指数も、`deg[n] = n²` も、`#E[n] = n²` も、場所の分類定理も使っていない。**

## ★★★★★場合 B の締め

第 141 の `f_P · f_{−P} = c(x − x_P)^n` に `count` を当てると

    count(μ f_P) + count(μ f_{−P}) = n · count(μ x)

★第 147 で `ι(f_P)` と `f_{−P}` は同伴、第 148 で `count(μ ι a) = count(μ a)`。
したがって `count(μ f_{−P}) = count(μ f_P)` となり

    2·count(μ f_P) = n · count(μ x)

★★第 145 で `count(μ x)` は偶数だから **`n ∣ count(μ f_P)`**。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `count_unit_image_eq_zero` | ★★単元の像の count は 0(第 128) |
| `count_genX_sub_const` | ★★場合 B で `count(μ(x − c)) = count(μ x)` |
| `count_hyperInv` | ★★★★`count(μ ι a) = count(μ a)` |
| `dvd_count_case_B` | ★★★★★★**場合 B の結論** |
| `exists_fNeg` | ★`f_{−P}` の存在(第 126 を `−P` に当てる) |
| `dvd_count_pullback` | ★★★★★★★**D1'** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine Polynomial IsDedekindDomain nonZeroDivisors

variable {F : Type} [Field F] (W : WeierstrassCurve.Affine F)
  [inst : IsDedekindDomain W.CoordinateRing]

/-! ## ★★`count` の小道具 -/

/-- ★★単元の像の `count` は 0(第 128 により単元は定数)。 -/
theorem count_unit_image_eq_zero [W.IsElliptic] (v : HeightOneSpectrum W.CoordinateRing)
    (μ : W.CoordinateRing →+* W.FunctionField)
    (hμF : ∀ c : F, μ (algebraMap F W.CoordinateRing c) = algebraMap F W.FunctionField c)
    {u : W.CoordinateRing} (hu : IsUnit u) :
    FractionalIdeal.count W.FunctionField v
      (FractionalIdeal.spanSingleton W.CoordinateRing⁰ (μ u)) = 0 := by
  obtain ⟨c, hc0, hc⟩ := isUnit_coordinateRing hu
  have hval : v.valuation W.FunctionField (μ u) = 1 := by
    rw [hc, hμF, IsScalarTower.algebraMap_apply F W.CoordinateRing W.FunctionField]
    exact valuation_algebraMap_isUnit v ((isUnit_iff_ne_zero.2 hc0).map
      (algebraMap F W.CoordinateRing))
  have hne : μ u ≠ 0 := by
    intro h0
    rw [h0, Valuation.map_zero] at hval
    exact absurd hval (by simp)
  rw [valuation_eq_exp_neg_count (K := W.FunctionField) v hne] at hval
  have hz : WithZero.exp (-FractionalIdeal.count W.FunctionField v
      (FractionalIdeal.spanSingleton W.CoordinateRing⁰ (μ u))) = WithZero.exp (0 : ℤ) := by
    rw [hval, WithZero.exp_zero]
  have h2 := WithZero.exp_inj.1 hz
  omega

/-- ★定数の像の `count` は 0。 -/
theorem count_algebraMap_const (v : HeightOneSpectrum W.CoordinateRing) {c : F} (hc : c ≠ 0) :
    FractionalIdeal.count W.FunctionField v
      (FractionalIdeal.spanSingleton W.CoordinateRing⁰ (algebraMap F W.FunctionField c)) = 0 := by
  have hval : v.valuation W.FunctionField (algebraMap F W.FunctionField c) = 1 := by
    rw [IsScalarTower.algebraMap_apply F W.CoordinateRing W.FunctionField]
    exact valuation_algebraMap_isUnit v ((isUnit_iff_ne_zero.2 hc).map
      (algebraMap F W.CoordinateRing))
  have hne : algebraMap F W.FunctionField c ≠ 0 := by
    intro h0
    rw [h0, Valuation.map_zero] at hval
    exact absurd hval (by simp)
  rw [valuation_eq_exp_neg_count (K := W.FunctionField) v hne] at hval
  have hz : WithZero.exp (-FractionalIdeal.count W.FunctionField v
      (FractionalIdeal.spanSingleton W.CoordinateRing⁰ (algebraMap F W.FunctionField c)))
      = WithZero.exp (0 : ℤ) := by rw [hval, WithZero.exp_zero]
  have := WithZero.exp_inj.1 hz
  omega

theorem count_mul' (v : HeightOneSpectrum W.CoordinateRing) {a b : W.FunctionField}
    (ha : a ≠ 0) (hb : b ≠ 0) :
    FractionalIdeal.count W.FunctionField v
        (FractionalIdeal.spanSingleton W.CoordinateRing⁰ (a * b))
      = FractionalIdeal.count W.FunctionField v
          (FractionalIdeal.spanSingleton W.CoordinateRing⁰ a)
        + FractionalIdeal.count W.FunctionField v
          (FractionalIdeal.spanSingleton W.CoordinateRing⁰ b) := by
  rw [← FractionalIdeal.spanSingleton_mul_spanSingleton,
    FractionalIdeal.count_mul W.FunctionField v
      (FractionalIdeal.spanSingleton_ne_zero_iff.2 ha)
      (FractionalIdeal.spanSingleton_ne_zero_iff.2 hb)]

theorem count_pow' (v : HeightOneSpectrum W.CoordinateRing) (a : W.FunctionField) (n : ℕ) :
    FractionalIdeal.count W.FunctionField v
        (FractionalIdeal.spanSingleton W.CoordinateRing⁰ (a ^ n))
      = n * FractionalIdeal.count W.FunctionField v
        (FractionalIdeal.spanSingleton W.CoordinateRing⁰ a) := by
  rw [← FractionalIdeal.spanSingleton_pow, FractionalIdeal.count_pow]

/-- ★★場合 B では `count(μ(x − c)) = count(μ x)`。 -/
theorem count_genX_sub_const (v : HeightOneSpectrum W.CoordinateRing)
    (μ : W.CoordinateRing →+* W.FunctionField)
    (hx : 1 < v.valuation W.FunctionField (μ (genX W))) (c : F) :
    FractionalIdeal.count W.FunctionField v (FractionalIdeal.spanSingleton W.CoordinateRing⁰
        (μ (genX W) - algebraMap F W.FunctionField c))
      = FractionalIdeal.count W.FunctionField v
        (FractionalIdeal.spanSingleton W.CoordinateRing⁰ (μ (genX W))) := by
  set w := v.valuation W.FunctionField with hw
  have hcv : w (algebraMap F W.FunctionField c) ≤ 1 := valuation_algebraMap_field W v c
  have hval : w (μ (genX W) - algebraMap F W.FunctionField c) = w (μ (genX W)) := by
    rw [sub_eq_add_neg, Valuation.map_add_of_distinct_val w (by
      rw [Valuation.map_neg]
      exact ne_of_gt (lt_of_le_of_lt hcv hx)), Valuation.map_neg]
    exact max_eq_left (le_of_lt (lt_of_le_of_lt hcv hx))
  have hx0 : μ (genX W) ≠ 0 := by
    intro h0
    rw [h0, Valuation.map_zero] at hx
    exact absurd hx (by simp)
  have hs0 : μ (genX W) - algebraMap F W.FunctionField c ≠ 0 := by
    intro h0
    rw [h0, Valuation.map_zero] at hval
    exact absurd hval.symm (ne_of_gt (lt_trans zero_lt_one hx))
  rw [valuation_eq_exp_neg_count (K := W.FunctionField) v hs0,
    valuation_eq_exp_neg_count (K := W.FunctionField) v hx0, WithZero.exp_inj] at hval
  omega

/-- ★★★★`count(μ ι a) = count(μ a)`(場合 B)。 -/
theorem count_hyperInv (h2 : IsUnit (2 : F)) (v : HeightOneSpectrum W.CoordinateRing)
    (μ : W.CoordinateRing →+* W.FunctionField) (hμinj : Function.Injective μ)
    (hμF : ∀ c : F, μ (algebraMap F W.CoordinateRing c) = algebraMap F W.FunctionField c)
    (hx : 1 < v.valuation W.FunctionField (μ (genX W))) {a : W.CoordinateRing} (ha : a ≠ 0) :
    FractionalIdeal.count W.FunctionField v
        (FractionalIdeal.spanSingleton W.CoordinateRing⁰ (μ (hyperInv W a)))
      = FractionalIdeal.count W.FunctionField v
        (FractionalIdeal.spanSingleton W.CoordinateRing⁰ (μ a)) := by
  have hval := valuation_hyperInv W h2 v μ hμF hx a
  have ha0 : μ a ≠ 0 := fun h0 => ha (hμinj (by rw [h0, map_zero]))
  have hia0 : μ (hyperInv W a) ≠ 0 := by
    intro h0
    rw [h0, Valuation.map_zero] at hval
    exact ((Valuation.zero_iff _).1 hval.symm) |> ha0
  rw [valuation_eq_exp_neg_count (K := W.FunctionField) v hia0,
    valuation_eq_exp_neg_count (K := W.FunctionField) v ha0, WithZero.exp_inj] at hval
  omega

/-! ## ★★★★★★場合 B の結論 -/

/-- ★★★★★★**場合 B —— 指数は `n` で割れる**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★第 141 の `f_P · f_{−P} = c(x − x_P)^n` に `count` を当て、
第 147(`ι(f_P) ~ f_{−P}`)と第 148(`count(μ ι a) = count(μ a)`)で
`count(μ f_{−P}) = count(μ f_P)` に潰し、第 145(`count(μ x)` は偶数)で割る。 -/
theorem dvd_count_case_B (h2 : IsUnit (2 : F)) [W.IsElliptic]
    (v : HeightOneSpectrum W.CoordinateRing)
    (μ : W.CoordinateRing →+* W.FunctionField) (hμinj : Function.Injective μ)
    (hμF : ∀ c : F, μ (algebraMap F W.CoordinateRing c) = algebraMap F W.FunctionField c)
    (hx : 1 < v.valuation W.FunctionField (μ (genX W)))
    {x y : F} (h : W.Nonsingular x y) (n : ℕ) (fP fN : W.CoordinateRing)
    (hfP : (CoordinateRing.XYIdeal W x (Polynomial.C y)) ^ n = Ideal.span {fP})
    (hfN : (CoordinateRing.XYIdeal W x (Polynomial.C (W.negY x y))) ^ n = Ideal.span {fN}) :
    (n : ℤ) ∣ FractionalIdeal.count W.FunctionField v
      (FractionalIdeal.spanSingleton W.CoordinateRing⁰ (μ fP)) := by
  obtain ⟨c, hc0, hprod⟩ := mu_fP_mul_mu_fNegP W h n fP fN hfP hfN μ hμF
  have hfP0 : fP ≠ 0 := fP_ne_zero W n fP hfP
  have hfN0 : fN ≠ 0 := fP_ne_zero W n fN hfN
  have hmP : μ fP ≠ 0 := fun h0 => hfP0 (hμinj (by rw [h0, map_zero]))
  have hmN : μ fN ≠ 0 := fun h0 => hfN0 (hμinj (by rw [h0, map_zero]))
  have hcK : algebraMap F W.FunctionField c ≠ 0 := fun h0 => hc0
    ((algebraMap F W.FunctionField).injective (by rw [h0, map_zero]))
  have hxc : μ (genX W) - algebraMap F W.FunctionField x ≠ 0 := by
    intro h0
    have hcs := count_genX_sub_const W v μ hx x
    rw [h0] at hcs
    simp only [FractionalIdeal.spanSingleton_zero, FractionalIdeal.count_zero] at hcs
    have hneg := count_genX_neg W v μ hx
    omega
  have hcount := congrArg (fun t => FractionalIdeal.count W.FunctionField v
    (FractionalIdeal.spanSingleton W.CoordinateRing⁰ t)) hprod
  rw [count_mul' W v hmP hmN, count_mul' W v hcK (pow_ne_zero n hxc), count_pow' W v,
    count_algebraMap_const W v hc0, count_genX_sub_const W v μ hx x, zero_add] at hcount
  obtain ⟨u, hu⟩ := hyperInv_fP_assoc W n fP fN hfP hfN
  have hNeq : FractionalIdeal.count W.FunctionField v
      (FractionalIdeal.spanSingleton W.CoordinateRing⁰ (μ fN))
      = FractionalIdeal.count W.FunctionField v
        (FractionalIdeal.spanSingleton W.CoordinateRing⁰ (μ fP)) := by
    rw [← hu, map_mul, count_mul' W v
      (fun h0 => hmN (by rw [← hu, map_mul, h0, zero_mul]))
      (fun h0 => by
        have hun : IsUnit (μ (u : W.CoordinateRing)) := u.isUnit.map μ
        rw [h0] at hun
        exact absurd hun not_isUnit_zero),
      count_unit_image_eq_zero W v μ hμF u.isUnit, add_zero,
      count_hyperInv W h2 v μ hμinj hμF hx hfP0]
  rw [hNeq] at hcount
  obtain ⟨k, hk⟩ := even_count_genX W h2 v μ hμF hx
  refine ⟨k, ?_⟩
  rw [hk] at hcount
  have hr : (n : ℤ) * (2 * k) = 2 * ((n : ℤ) * k) := by ring
  rw [hr] at hcount
  linarith

/-! ## ★★★★★★★D1' -/

/-- ★`f_{−P}` の存在——第 126 を `−P` に当てる。 -/
theorem exists_fNeg [W.IsElliptic] [DecidableEq F] {x y : F} (h : W.Nonsingular x y) (n : ℕ)
    (hP : n • (Point.some x y h) = 0) :
    ∃ fN : W.CoordinateRing, fN ≠ 0 ∧
      (CoordinateRing.XYIdeal W x (Polynomial.C (W.negY x y))) ^ n = Ideal.span {fN} := by
  have hns : W.Nonsingular x (W.negY x y) := (nonsingular_neg ..).mpr h
  have hnegP : n • (Point.some x (W.negY x y) hns) = 0 := by
    have hneg : Point.some x (W.negY x y) hns = -(Point.some x y h) := rfl
    rw [hneg, smul_neg, hP, neg_zero]
  exact xyIdeal_pow_isPrincipal_integral hns n hnegP

/-- ★★★★★★★**D1' —— 各素点での指数は `n` で割れる**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★2 つの場合を合流させた。★★**分岐指数も `deg[n] = n²` も `#E[n] = n²` も
場所の分類定理も使っていない。** -/
theorem dvd_count_pullback (h2 : IsUnit (2 : F)) [W.IsElliptic] [DecidableEq F]
    (μ : W.CoordinateRing →+* W.FunctionField) (hμinj : Function.Injective μ)
    (hμF : ∀ c : F, μ (algebraMap F W.CoordinateRing c) = algebraMap F W.FunctionField c)
    {x y : F} (h : W.Nonsingular x y) (n : ℕ) (hP : n • (Point.some x y h) = 0)
    (fP : W.CoordinateRing)
    (hfP : (CoordinateRing.XYIdeal W x (Polynomial.C y)) ^ n = Ideal.span {fP})
    (v : HeightOneSpectrum W.CoordinateRing) :
    (n : ℤ) ∣ FractionalIdeal.count W.FunctionField v
      (FractionalIdeal.spanSingleton W.CoordinateRing⁰ (μ fP)) := by
  by_cases hcase : ∀ r : W.CoordinateRing, v.valuation W.FunctionField (μ r) ≤ 1
  · refine dvd_count_of_span_pair_pow v μ hμinj hcase
      (CoordinateRing.XClass_ne_zero x) (CoordinateRing.YClass_ne_zero (Polynomial.C y))
      (fP_ne_zero W n fP hfP) ?_
    rw [← hfP, CoordinateRing.XYIdeal]
  · rw [not_forall] at hcase
    obtain ⟨r, hr⟩ := hcase
    have hx : 1 < v.valuation W.FunctionField (μ (genX W)) :=
      one_lt_valuation_genX W v μ hμF (not_le.1 hr)
    obtain ⟨fN, hfN0, hfN⟩ := exists_fNeg W h n hP
    exact dvd_count_case_B W h2 v μ hμinj hμF hx h n fP fN hfP hfN

/-! ## ★出典の紐付け(`.src`) -/

def dvd_count_pullback.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——div(f_P ∘ [n]) の各素点での指数が n で割れること)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
