import ABC3.Found.GaloisRep.ValuationEven

/-!
# Galois (G5) 第 146 ブロック —— **★★★★場合 B の基底計算**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★mathlib が偶奇構造を持っていた

場合 B の残りは「`a = p·1 + q·y` の count が `−m·deg N(a)` であること」である。
★**2026-08-20 実測**: mathlib の

    CoordinateRing.degree_norm_smul_basis :
      deg N(p•1 + q•y) = max (2·deg p) (2·deg q + 3)

が**まさにその偶奇構造**を与えている。★★`count(μ(p•1 + q•y))` の側も

    min(deg p · c_x, deg q · c_x + c_y)   ただし c_x = −2m, c_y = −3m

であり、偶数倍と奇数倍なので**決して一致しない**——これが
「`ν` は `F(x)` 上への制限で一意に決まる」ことの正体である。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `count_algebraMap_poly` | ★★★★**`count(μ p(x)) = deg p · count(μ x)`** |
| `valuation_genX_lt_genZ` | ★★★`w(μ x) < w(μ z)` |
| `valuation_genY_eq_genZ` | ★★★★**`w(μ y) = w(μ z)`**(`2y = z − a₁x − a₃`) |

★`valuation_genY_eq_genZ` は「`z` の項が単独で最大」から出る
(mathlib の `Valuation.map_add_of_distinct_val`)。
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine Polynomial IsDedekindDomain nonZeroDivisors

variable {F : Type} [Field F] (W : WeierstrassCurve.Affine F)
  [inst : IsDedekindDomain W.CoordinateRing]

/-- ★★★★**`count(μ p(x)) = deg p · count(μ x)`**(場合 B)。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★主係数は 0 でない定数なので付値 1、したがって最高次が単独で最大になる。 -/
theorem count_algebraMap_poly (v : HeightOneSpectrum W.CoordinateRing)
    (μ : W.CoordinateRing →+* W.FunctionField)
    (hμF : ∀ c : F, μ (algebraMap F W.CoordinateRing c) = algebraMap F W.FunctionField c)
    (hx : 1 < v.valuation W.FunctionField (μ (genX W)))
    {p : Polynomial F} (hp : p ≠ 0) :
    FractionalIdeal.count W.FunctionField v (FractionalIdeal.spanSingleton W.CoordinateRing⁰
        (μ (algebraMap (Polynomial F) W.CoordinateRing p)))
      = (p.natDegree : ℤ) * FractionalIdeal.count W.FunctionField v
        (FractionalIdeal.spanSingleton W.CoordinateRing⁰ (μ (genX W))) := by
  set w := v.valuation W.FunctionField with hw
  have hcomp : μ.comp (algebraMap F W.CoordinateRing) = algebraMap F W.FunctionField :=
    RingHom.ext hμF
  have hleadval : w (algebraMap F W.FunctionField p.leadingCoeff) = 1 := by
    rw [IsScalarTower.algebraMap_apply F W.CoordinateRing W.FunctionField]
    exact valuation_algebraMap_isUnit v ((isUnit_iff_ne_zero.2
      (Polynomial.leadingCoeff_ne_zero.2 hp)).map (algebraMap F W.CoordinateRing))
  have hev : μ (algebraMap (Polynomial F) W.CoordinateRing p)
      = Polynomial.eval₂ (algebraMap F W.FunctionField) (μ (genX W)) p := by
    rw [← eval₂_genX_eq_algebraMap, Polynomial.hom_eval₂, hcomp]
  have hval : w (μ (algebraMap (Polynomial F) W.CoordinateRing p))
      = w (μ (genX W)) ^ p.natDegree := by
    rw [hev]
    exact valuation_eval₂_of_leadingCoeff w _ (valuation_algebraMap_field W v) hx hp hleadval
  have hx0 : μ (genX W) ≠ 0 := by
    intro h0
    rw [h0, Valuation.map_zero] at hx
    exact absurd hx (by simp)
  have hp0 : μ (algebraMap (Polynomial F) W.CoordinateRing p) ≠ 0 := by
    intro h0
    rw [h0, Valuation.map_zero] at hval
    exact absurd hval.symm (ne_of_gt (pow_pos (lt_trans zero_lt_one hx) _))
  rw [valuation_eq_exp_neg_count (K := W.FunctionField) v hp0,
    valuation_eq_exp_neg_count (K := W.FunctionField) v hx0,
    ← WithZero.exp_nsmul, WithZero.exp_inj, nsmul_eq_mul, mul_neg] at hval
  exact neg_inj.1 hval

/-- ★★★場合 B では `w(μ x) < w(μ z)`。 -/
theorem valuation_genX_lt_genZ (h2 : IsUnit (2 : F)) (v : HeightOneSpectrum W.CoordinateRing)
    (μ : W.CoordinateRing →+* W.FunctionField)
    (hμF : ∀ c : F, μ (algebraMap F W.CoordinateRing c) = algebraMap F W.FunctionField c)
    (hx : 1 < v.valuation W.FunctionField (μ (genX W))) :
    v.valuation W.FunctionField (μ (genX W)) < v.valuation W.FunctionField (μ (genZ W)) := by
  set w := v.valuation W.FunctionField with hw
  have hkey := valuation_genZ_sq W h2 v μ hμF hx
  by_contra hle
  rw [not_lt] at hle
  have h3 : w (μ (genX W)) ^ 3 ≤ w (μ (genX W)) ^ 2 := by
    rw [← hkey]; exact pow_le_pow_left' hle 2
  exact absurd h3 (not_le.2 (pow_lt_pow_right₀ hx (by norm_num)))

/-- ★★★★**場合 B では `w(μ y) = w(μ z)`**——`2y = z − a₁x − a₃` で `z` の項が単独で最大。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★これにより mathlib の基底 `{1, y}` での計算が `z` での計算と一致する。 -/
theorem valuation_genY_eq_genZ (h2 : IsUnit (2 : F)) (v : HeightOneSpectrum W.CoordinateRing)
    (μ : W.CoordinateRing →+* W.FunctionField)
    (hμF : ∀ c : F, μ (algebraMap F W.CoordinateRing c) = algebraMap F W.FunctionField c)
    (hx : 1 < v.valuation W.FunctionField (μ (genX W))) :
    v.valuation W.FunctionField (μ (genY W)) = v.valuation W.FunctionField (μ (genZ W)) := by
  set w := v.valuation W.FunctionField with hw
  have hxz := valuation_genX_lt_genZ W h2 v μ hμF hx
  have hdecomp : (2 : W.FunctionField) * μ (genY W)
      = μ (genZ W) - (μ (algebraMap F W.CoordinateRing W.a₁) * μ (genX W)
        + μ (algebraMap F W.CoordinateRing W.a₃)) := by
    simp only [genZ, map_add, map_mul, map_ofNat]
    ring
  have hR : w (μ (algebraMap F W.CoordinateRing W.a₁) * μ (genX W)
      + μ (algebraMap F W.CoordinateRing W.a₃)) < w (μ (genZ W)) := by
    refine lt_of_le_of_lt (le_trans (Valuation.map_add w _ _) (max_le ?_ ?_)) hxz
    · rw [Valuation.map_mul, hμF]
      calc w (algebraMap F W.FunctionField W.a₁) * w (μ (genX W))
          ≤ 1 * w (μ (genX W)) := mul_le_mul' (valuation_algebraMap_field W v _) le_rfl
        _ = w (μ (genX W)) := one_mul _
    · rw [hμF]
      exact le_trans (valuation_algebraMap_field W v _) (le_of_lt hx)
  have h2u : w ((2 : W.FunctionField)) = 1 := by
    have h2e : (2 : W.FunctionField) = algebraMap F W.FunctionField 2 := (map_ofNat _ _).symm
    rw [h2e, IsScalarTower.algebraMap_apply F W.CoordinateRing W.FunctionField]
    exact valuation_algebraMap_isUnit v (h2.map (algebraMap F W.CoordinateRing))
  have hsub : w (μ (genZ W) - (μ (algebraMap F W.CoordinateRing W.a₁) * μ (genX W)
      + μ (algebraMap F W.CoordinateRing W.a₃))) = w (μ (genZ W)) := by
    rw [sub_eq_add_neg, Valuation.map_add_of_distinct_val w (by
      rw [Valuation.map_neg]; exact ne_of_gt hR), Valuation.map_neg]
    exact max_eq_left (le_of_lt hR)
  have hfin := congrArg w hdecomp
  rw [Valuation.map_mul, h2u, one_mul, hsub] at hfin
  exact hfin

/-! ## ★出典の紐付け(`.src`) -/

def count_algebraMap_poly.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——無限遠の上での多項式の付値)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
