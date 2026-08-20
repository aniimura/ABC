import ABC3.Found.GaloisRep.HyperInv

/-!
# Galois (G5) 第 148 ブロック —— **★★★★★★対合は付値を保つ**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★場合 B の核心

    w(μ ι(a)) = w(μ a)          (場合 B、すなわち `1 < w(μ x)`)

★これで `count(μ f_P) = count(μ f_{−P})` が出て、第 141 の
`f_P · f_{−P} = c(x − x_P)^n` と合わせれば `2·count(μ f_P) = n·count(μ x)`、
第 145(`count(μ x)` は偶数)から **`n ∣ count(μ f_P)`** となる。

## ★★★★★機構——偶奇で「決して一致しない」

mathlib の基底 `{1, y}` で `a = p·1 + q·y` と書くと:

| 項 | count |
|---|---|
| `p(x)` | `deg p · c_x`(**偶数**倍) |
| `q(x)·y` | `deg q · c_x + c_y`(**奇数**倍) |

★根拠は `2c_y = 3c_x`(第 145・146)と `c_x < 0`。
`2·deg p = 2·deg q + 3` は成り立たないので、**2 項の付値は決して一致しない**。
★★したがって和の付値は max であり、`y` を `ι(y)` に取り替えても
(`w(μ ι y) = w(μ y)` なので)**同じ max になる**。

★★★これが「`ν` は `F(x)` 上への制限で一意に決まる」ことの Lean での正体である。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `valuation_algebraMap_poly_eq` | ★★★`w(μ p(x)) = w(μ x)^{deg p}` |
| `count_genX_neg` | ★★`count(μ x) < 0` |
| `count_genY_relation` | ★★★★**`2·count(μ y) = 3·count(μ x)`** |
| `valuation_hyperInv_genY` | ★★★★`w(μ ι(y)) = w(μ y)` |
| `valuation_poly_ne_genY_term` | ★★★★★**2 項の付値は一致しない**(偶奇) |
| `valuation_hyperInv` | ★★★★★★**`w(μ ι(a)) = w(μ a)`** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine Polynomial IsDedekindDomain nonZeroDivisors

variable {F : Type} [Field F] (W : WeierstrassCurve.Affine F)
  [inst : IsDedekindDomain W.CoordinateRing]

/-! ## ★★★多項式の付値と `count` -/

/-- ★★★`w(μ p(x)) = w(μ x)^{deg p}`(場合 B)。 -/
theorem valuation_algebraMap_poly_eq (v : HeightOneSpectrum W.CoordinateRing)
    (μ : W.CoordinateRing →+* W.FunctionField)
    (hμF : ∀ c : F, μ (algebraMap F W.CoordinateRing c) = algebraMap F W.FunctionField c)
    (hx : 1 < v.valuation W.FunctionField (μ (genX W)))
    {p : Polynomial F} (hp : p ≠ 0) :
    v.valuation W.FunctionField (μ (algebraMap (Polynomial F) W.CoordinateRing p))
      = (v.valuation W.FunctionField (μ (genX W))) ^ p.natDegree := by
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
  rw [hev]
  exact valuation_eval₂_of_leadingCoeff w _ (valuation_algebraMap_field W v) hx hp hleadval

theorem valuation_algebraMap_poly_ne_zero (v : HeightOneSpectrum W.CoordinateRing)
    (μ : W.CoordinateRing →+* W.FunctionField)
    (hμF : ∀ c : F, μ (algebraMap F W.CoordinateRing c) = algebraMap F W.FunctionField c)
    (hx : 1 < v.valuation W.FunctionField (μ (genX W)))
    {p : Polynomial F} (hp : p ≠ 0) :
    v.valuation W.FunctionField (μ (algebraMap (Polynomial F) W.CoordinateRing p)) ≠ 0 := by
  rw [valuation_algebraMap_poly_eq W v μ hμF hx hp]
  exact ne_of_gt (pow_pos (lt_trans zero_lt_one hx) _)

theorem valuation_genY_ne_zero (h2 : IsUnit (2 : F)) (v : HeightOneSpectrum W.CoordinateRing)
    (μ : W.CoordinateRing →+* W.FunctionField)
    (hμF : ∀ c : F, μ (algebraMap F W.CoordinateRing c) = algebraMap F W.FunctionField c)
    (hx : 1 < v.valuation W.FunctionField (μ (genX W))) :
    v.valuation W.FunctionField (μ (genY W)) ≠ 0 := by
  rw [valuation_genY_eq_genZ W h2 v μ hμF hx]
  intro h0
  have hkey := valuation_genZ_sq W h2 v μ hμF hx
  rw [h0, zero_pow (by norm_num)] at hkey
  exact absurd hkey.symm (ne_of_gt (pow_pos (lt_trans zero_lt_one hx) 3))

/-! ## ★★★★`count` の関係式 -/

/-- ★★場合 B では `count(μ x) < 0`。 -/
theorem count_genX_neg (v : HeightOneSpectrum W.CoordinateRing)
    (μ : W.CoordinateRing →+* W.FunctionField)
    (hx : 1 < v.valuation W.FunctionField (μ (genX W))) :
    FractionalIdeal.count W.FunctionField v
      (FractionalIdeal.spanSingleton W.CoordinateRing⁰ (μ (genX W))) < 0 := by
  have hx0 : μ (genX W) ≠ 0 := by
    intro h0
    rw [h0, Valuation.map_zero] at hx
    exact absurd hx (by simp)
  have hb := valuation_eq_exp_neg_count (K := W.FunctionField) v hx0
  rw [hb] at hx
  have hpos : (0 : ℤ) < -FractionalIdeal.count W.FunctionField v
      (FractionalIdeal.spanSingleton W.CoordinateRing⁰ (μ (genX W))) := by
    have h1 : (WithZero.exp (0 : ℤ)) < WithZero.exp
        (-FractionalIdeal.count W.FunctionField v
          (FractionalIdeal.spanSingleton W.CoordinateRing⁰ (μ (genX W)))) := by
      rw [WithZero.exp_zero]; exact hx
    exact WithZero.exp_lt_exp.1 h1
  omega

/-- ★★★★**`2·count(μ y) = 3·count(μ x)`**——これが偶奇の根拠である。 -/
theorem count_genY_relation (h2 : IsUnit (2 : F)) (v : HeightOneSpectrum W.CoordinateRing)
    (μ : W.CoordinateRing →+* W.FunctionField)
    (hμF : ∀ c : F, μ (algebraMap F W.CoordinateRing c) = algebraMap F W.FunctionField c)
    (hx : 1 < v.valuation W.FunctionField (μ (genX W))) :
    2 * FractionalIdeal.count W.FunctionField v
        (FractionalIdeal.spanSingleton W.CoordinateRing⁰ (μ (genY W)))
      = 3 * FractionalIdeal.count W.FunctionField v
        (FractionalIdeal.spanSingleton W.CoordinateRing⁰ (μ (genX W))) := by
  set w := v.valuation W.FunctionField with hw
  have hkey : w (μ (genY W)) ^ 2 = w (μ (genX W)) ^ 3 := by
    rw [valuation_genY_eq_genZ W h2 v μ hμF hx]
    exact valuation_genZ_sq W h2 v μ hμF hx
  have hx0 : μ (genX W) ≠ 0 := by
    intro h0
    rw [h0, Valuation.map_zero] at hx
    exact absurd hx (by simp)
  have hy0 : μ (genY W) ≠ 0 := by
    intro h0
    rw [h0, Valuation.map_zero, zero_pow (by norm_num)] at hkey
    exact absurd hkey.symm (ne_of_gt (pow_pos (lt_trans zero_lt_one hx) 3))
  rw [valuation_eq_exp_neg_count (K := W.FunctionField) v hx0,
    valuation_eq_exp_neg_count (K := W.FunctionField) v hy0,
    ← WithZero.exp_nsmul, ← WithZero.exp_nsmul, WithZero.exp_inj,
    nsmul_eq_mul, nsmul_eq_mul] at hkey
  omega

/-! ## ★★★★対合は付値を保つ -/

/-- ★★★★場合 B では `w(μ ι(y)) = w(μ y)`——`y` の項が単独で最大。 -/
theorem valuation_hyperInv_genY (h2 : IsUnit (2 : F)) (v : HeightOneSpectrum W.CoordinateRing)
    (μ : W.CoordinateRing →+* W.FunctionField)
    (hμF : ∀ c : F, μ (algebraMap F W.CoordinateRing c) = algebraMap F W.FunctionField c)
    (hx : 1 < v.valuation W.FunctionField (μ (genX W))) :
    v.valuation W.FunctionField (μ (hyperInv W (genY W)))
      = v.valuation W.FunctionField (μ (genY W)) := by
  set w := v.valuation W.FunctionField with hw
  have hxy : w (μ (genX W)) < w (μ (genY W)) := by
    rw [valuation_genY_eq_genZ W h2 v μ hμF hx]
    exact valuation_genX_lt_genZ W h2 v μ hμF hx
  have hdec : μ (hyperInv W (genY W))
      = -μ (genY W) - (μ (algebraMap F W.CoordinateRing W.a₁) * μ (genX W)
        + μ (algebraMap F W.CoordinateRing W.a₃)) := by
    rw [hyperInv_genY_eq, map_sub, map_sub, map_mul, map_neg]
    ring
  have hR : w (μ (algebraMap F W.CoordinateRing W.a₁) * μ (genX W)
      + μ (algebraMap F W.CoordinateRing W.a₃)) < w (μ (genY W)) := by
    refine lt_of_le_of_lt (le_trans (Valuation.map_add w _ _) (max_le ?_ ?_)) hxy
    · rw [Valuation.map_mul, hμF]
      calc w (algebraMap F W.FunctionField W.a₁) * w (μ (genX W))
          ≤ 1 * w (μ (genX W)) := mul_le_mul' (valuation_algebraMap_field W v _) le_rfl
        _ = w (μ (genX W)) := one_mul _
    · rw [hμF]
      exact le_trans (valuation_algebraMap_field W v _) (le_of_lt hx)
  rw [hdec, sub_eq_add_neg,
    Valuation.map_add_of_distinct_val w (by
      rw [Valuation.map_neg, Valuation.map_neg]; exact ne_of_gt hR),
    Valuation.map_neg, Valuation.map_neg]
  exact max_eq_left (le_of_lt hR)

/-- ★★★★★**2 項の付値は決して一致しない**——偶奇による。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★`2·deg p = 2·deg q + 3` は成り立たない。 -/
theorem valuation_poly_ne_genY_term (h2 : IsUnit (2 : F))
    (v : HeightOneSpectrum W.CoordinateRing)
    (μ : W.CoordinateRing →+* W.FunctionField)
    (hμF : ∀ c : F, μ (algebraMap F W.CoordinateRing c) = algebraMap F W.FunctionField c)
    (hx : 1 < v.valuation W.FunctionField (μ (genX W)))
    {p q : Polynomial F} (hp : p ≠ 0) (hq : q ≠ 0) :
    v.valuation W.FunctionField (μ (algebraMap (Polynomial F) W.CoordinateRing p))
      ≠ v.valuation W.FunctionField
        (μ (algebraMap (Polynomial F) W.CoordinateRing q) * μ (genY W)) := by
  set w := v.valuation W.FunctionField with hw
  intro heq
  have hAne : μ (algebraMap (Polynomial F) W.CoordinateRing p) ≠ 0 := fun h0 =>
    valuation_algebraMap_poly_ne_zero W v μ hμF hx hp (by rw [h0, Valuation.map_zero])
  have hqne : μ (algebraMap (Polynomial F) W.CoordinateRing q) ≠ 0 := fun h0 =>
    valuation_algebraMap_poly_ne_zero W v μ hμF hx hq (by rw [h0, Valuation.map_zero])
  have hyne : μ (genY W) ≠ 0 := fun h0 =>
    valuation_genY_ne_zero W h2 v μ hμF hx (by rw [h0, Valuation.map_zero])
  have hBne : μ (algebraMap (Polynomial F) W.CoordinateRing q) * μ (genY W) ≠ 0 :=
    mul_ne_zero hqne hyne
  rw [valuation_eq_exp_neg_count (K := W.FunctionField) v hAne,
    valuation_eq_exp_neg_count (K := W.FunctionField) v hBne, WithZero.exp_inj] at heq
  have hsplit : FractionalIdeal.count W.FunctionField v (FractionalIdeal.spanSingleton
        W.CoordinateRing⁰ (μ (algebraMap (Polynomial F) W.CoordinateRing q) * μ (genY W)))
      = FractionalIdeal.count W.FunctionField v (FractionalIdeal.spanSingleton
          W.CoordinateRing⁰ (μ (algebraMap (Polynomial F) W.CoordinateRing q)))
        + FractionalIdeal.count W.FunctionField v (FractionalIdeal.spanSingleton
          W.CoordinateRing⁰ (μ (genY W))) := by
    rw [← FractionalIdeal.spanSingleton_mul_spanSingleton,
      FractionalIdeal.count_mul W.FunctionField v
        (FractionalIdeal.spanSingleton_ne_zero_iff.2 hqne)
        (FractionalIdeal.spanSingleton_ne_zero_iff.2 hyne)]
  rw [hsplit, count_algebraMap_poly W v μ hμF hx hp,
    count_algebraMap_poly W v μ hμF hx hq] at heq
  have hcy := count_genY_relation W h2 v μ hμF hx
  have hcx := count_genX_neg W v μ hx
  have hcancel : ((2 * p.natDegree : ℕ) : ℤ) = ((2 * q.natDegree + 3 : ℕ) : ℤ) := by
    push_cast
    refine mul_right_cancel₀ (ne_of_lt hcx) ?_
    linarith [heq, hcy]
  have hnat : 2 * p.natDegree = 2 * q.natDegree + 3 := by exact_mod_cast hcancel
  omega

/-- ★★★★★★**場合 B では `w(μ ι(a)) = w(μ a)`**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★基底 `{1, y}` で分解し、2 項の付値が決して一致しないことから和の付値は max。
★★`y` を `ι(y)` に取り替えても同じ max になる。 -/
theorem valuation_hyperInv (h2 : IsUnit (2 : F)) (v : HeightOneSpectrum W.CoordinateRing)
    (μ : W.CoordinateRing →+* W.FunctionField)
    (hμF : ∀ c : F, μ (algebraMap F W.CoordinateRing c) = algebraMap F W.FunctionField c)
    (hx : 1 < v.valuation W.FunctionField (μ (genX W))) (a : W.CoordinateRing) :
    v.valuation W.FunctionField (μ (hyperInv W a))
      = v.valuation W.FunctionField (μ a) := by
  set w := v.valuation W.FunctionField with hw
  obtain ⟨p, q, hpq⟩ := CoordinateRing.exists_smul_basis_eq a
  have ha : a = algebraMap (Polynomial F) W.CoordinateRing p
      + algebraMap (Polynomial F) W.CoordinateRing q * genY W := by
    rw [← hpq, Algebra.smul_def, Algebra.smul_def, mul_one]; rfl
  have hia : hyperInv W a = algebraMap (Polynomial F) W.CoordinateRing p
      + algebraMap (Polynomial F) W.CoordinateRing q * hyperInv W (genY W) := by
    rw [ha, map_add, map_mul, hyperInv_algebraMap_poly, hyperInv_algebraMap_poly]
  have hL : μ (hyperInv W a) = μ (algebraMap (Polynomial F) W.CoordinateRing p)
      + μ (algebraMap (Polynomial F) W.CoordinateRing q) * μ (hyperInv W (genY W)) := by
    rw [hia, map_add, map_mul]
  have hRt : μ a = μ (algebraMap (Polynomial F) W.CoordinateRing p)
      + μ (algebraMap (Polynomial F) W.CoordinateRing q) * μ (genY W) := by
    rw [ha, map_add, map_mul]
  have hB : w (μ (algebraMap (Polynomial F) W.CoordinateRing q) * μ (hyperInv W (genY W)))
      = w (μ (algebraMap (Polynomial F) W.CoordinateRing q) * μ (genY W)) := by
    rw [Valuation.map_mul, Valuation.map_mul, valuation_hyperInv_genY W h2 v μ hμF hx]
  rcases eq_or_ne (w (μ (algebraMap (Polynomial F) W.CoordinateRing p)))
      (w (μ (algebraMap (Polynomial F) W.CoordinateRing q) * μ (genY W))) with heq | hne
  · by_cases hp : p = 0
    · by_cases hq : q = 0
      · rw [hL, hRt, hp, hq]; simp
      · exfalso
        rw [hp, map_zero, map_zero, Valuation.map_zero] at heq
        have hqne : μ (algebraMap (Polynomial F) W.CoordinateRing q) ≠ 0 := fun h0 =>
          valuation_algebraMap_poly_ne_zero W v μ hμF hx hq (by rw [h0, Valuation.map_zero])
        have hyne : μ (genY W) ≠ 0 := fun h0 =>
          valuation_genY_ne_zero W h2 v μ hμF hx (by rw [h0, Valuation.map_zero])
        exact mul_ne_zero hqne hyne ((Valuation.zero_iff w).1 heq.symm)
    · by_cases hq : q = 0
      · exfalso
        rw [hq, map_zero, map_zero, zero_mul, Valuation.map_zero] at heq
        exact valuation_algebraMap_poly_ne_zero W v μ hμF hx hp heq
      · exact absurd heq (valuation_poly_ne_genY_term W h2 v μ hμF hx hp hq)
  · rw [hL, hRt, Valuation.map_add_of_distinct_val w (by rw [hB]; exact hne),
      Valuation.map_add_of_distinct_val w hne, hB]

/-! ## ★出典の紐付け(`.src`) -/

def valuation_hyperInv.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——超楕円対合が無限遠の上の付値を保つこと)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
