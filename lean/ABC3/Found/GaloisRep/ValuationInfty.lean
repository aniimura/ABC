import ABC3.Found.GaloisRep.ValuationCount

/-!
# Galois (G5) 第 144 ブロック —— **★★★★★D1' 場合 B の入口**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★場合 B は「`x` の付値が 1 より大きい」で始まる

第 143 ブロックで場合 A(`∀ r, w(μ r) ≤ 1`)を閉じた。
★残るのは**その否定**、すなわち `∃ r, 1 < w(μ r)` の場合である。

★★本ブロックは**そこから `1 < w(μ x)` が出る**ことを示す:

    r が `F[x]` 上整      (第 116)
      p モニック、`p(r) = 0`
      もし `w(μ x) ≤ 1` なら `F[x]` の元の付値はすべて 1 以下
      ⟹ `w(p(μ r)) = w(μ r)^{deg p} ≠ 0`   (最高次が単独で最大)
      ⟹ `w(0) ≠ 0` となり矛盾

★★★**これは「無限遠点の上にいる」ことの初等的な言い換え**である。
場所の分類定理も射影モデルも使わない。

## ★★★★道具立て

| 補題 | 内容 |
|---|---|
| `valuation_algebraMap_isUnit` | ★★単元の付値は 1 |
| `valuation_eval₂_le_one` | ★★★係数と `t` の付値が 1 以下なら `p(t)` も 1 以下 |
| `valuation_eval₂_of_monic` | ★★★★**モニックなら `w(p(t)) = w(t)^{deg p}`**(`w t > 1`) |
| `valuation_algebraMap_field` | ★定数の付値は 1 以下 |
| `one_lt_valuation_genX` | ★★★★★**場合 B では `1 < w(μ x)`** |

★`valuation_eval₂_of_monic` は mathlib の `Valuation.map_sum_eq_of_lt`
(和の中に単独最大があれば付値はそれに一致)で出る。
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine Polynomial IsDedekindDomain nonZeroDivisors

/-! ## ★★単元の付値 -/

variable {R K : Type} [CommRing R] [IsDedekindDomain R] [Field K] [Algebra R K]
  [IsFractionRing R K]

/-- ★★**単元の付値は 1 である**。 -/
theorem valuation_algebraMap_isUnit (v : HeightOneSpectrum R) {u : R} (hu : IsUnit u) :
    v.valuation K (algebraMap R K u) = 1 := by
  obtain ⟨uu, rfl⟩ := hu
  rw [v.valuation_of_algebraMap]
  have h1 : v.intValuation ((uu : Rˣ) : R) ≤ 1 := v.intValuation_le_one _
  have h2 : v.intValuation ((uu⁻¹ : Rˣ) : R) ≤ 1 := v.intValuation_le_one _
  have hmul : v.intValuation ((uu : Rˣ) : R) * v.intValuation ((uu⁻¹ : Rˣ) : R) = 1 := by
    rw [← Valuation.map_mul, ← Units.val_mul, mul_inv_cancel, Units.val_one, Valuation.map_one]
  by_contra hne
  have hlt : v.intValuation ((uu : Rˣ) : R) < 1 := lt_of_le_of_ne h1 hne
  have hcontra : v.intValuation ((uu : Rˣ) : R) * v.intValuation ((uu⁻¹ : Rˣ) : R)
      ≤ v.intValuation ((uu : Rˣ) : R) := by
    calc v.intValuation ((uu : Rˣ) : R) * v.intValuation ((uu⁻¹ : Rˣ) : R)
        ≤ v.intValuation ((uu : Rˣ) : R) * 1 := by exact mul_le_mul_left' h2 _
      _ = v.intValuation ((uu : Rˣ) : R) := mul_one _
  rw [hmul] at hcontra
  exact absurd hcontra (not_le.2 hlt)

/-! ## ★★★多項式の付値 -/

section Poly

variable {A L Γ₀ : Type} [CommRing A] [Field L] [LinearOrderedCommGroupWithZero Γ₀]

/-- ★★★係数と `t` の付値が 1 以下なら `p(t)` の付値も 1 以下。 -/
theorem valuation_eval₂_le_one (w : Valuation L Γ₀) (σ : A →+* L)
    (hσ : ∀ c : A, w (σ c) ≤ 1) {t : L} (ht : w t ≤ 1) (p : Polynomial A) :
    w (Polynomial.eval₂ σ t p) ≤ 1 := by
  classical
  have hsum : Polynomial.eval₂ σ t p = ∑ i ∈ p.support, σ (p.coeff i) * t ^ i := by
    rw [Polynomial.eval₂_eq_sum, Polynomial.sum]
  rw [hsum]
  refine Valuation.map_sum_le w fun i _ => ?_
  rw [Valuation.map_mul, Valuation.map_pow]
  calc w (σ (p.coeff i)) * w t ^ i ≤ 1 * 1 := mul_le_mul' (hσ _) (pow_le_one' ht i)
    _ = 1 := one_mul 1

/-- ★★★★**モニックなら `w(p(t)) = w(t)^{deg p}`**(`w t > 1`)。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★最高次の項が単独で最大になるので、mathlib の `Valuation.map_sum_eq_of_lt` が効く。 -/
theorem valuation_eval₂_of_monic [Nontrivial A] (w : Valuation L Γ₀) (σ : A →+* L)
    (hσ : ∀ c : A, w (σ c) ≤ 1) {t : L} (ht : 1 < w t)
    {p : Polynomial A} (hp : p.Monic) :
    w (Polynomial.eval₂ σ t p) = (w t) ^ p.natDegree := by
  classical
  have hp0 : p ≠ 0 := hp.ne_zero
  have hsum : Polynomial.eval₂ σ t p = ∑ i ∈ p.support, σ (p.coeff i) * t ^ i := by
    rw [Polynomial.eval₂_eq_sum, Polynomial.sum]
  have hdeg : p.natDegree ∈ p.support := Polynomial.natDegree_mem_support_of_nonzero hp0
  have htop : w (σ (p.coeff p.natDegree) * t ^ p.natDegree) = (w t) ^ p.natDegree := by
    rw [Polynomial.coeff_natDegree, hp.leadingCoeff, map_one, one_mul, Valuation.map_pow]
  rw [hsum, Valuation.map_sum_eq_of_lt w hdeg ?_, htop]
  intro i hi
  have hi' : i ∈ p.support ∧ i ≠ p.natDegree := by simpa [Finset.mem_sdiff] using hi
  rw [htop, Valuation.map_mul, Valuation.map_pow]
  calc w (σ (p.coeff i)) * w t ^ i ≤ 1 * w t ^ i := mul_le_mul' (hσ _) le_rfl
    _ = w t ^ i := one_mul _
    _ < (w t) ^ p.natDegree := pow_lt_pow_right₀ ht
        (lt_of_le_of_ne (Polynomial.le_natDegree_of_mem_supp i hi'.1) hi'.2)

end Poly

/-! ## ★★★★★場合 B の入口 -/

variable {F : Type} [Field F] (W : WeierstrassCurve.Affine F)
  [inst : IsDedekindDomain W.CoordinateRing]

/-- ★定数の付値は 1 以下。 -/
theorem valuation_algebraMap_field (v : HeightOneSpectrum W.CoordinateRing) (c : F) :
    v.valuation W.FunctionField (algebraMap F W.FunctionField c) ≤ 1 := by
  rcases eq_or_ne c 0 with rfl | hc
  · rw [map_zero, Valuation.map_zero]; exact zero_le_one
  · rw [IsScalarTower.algebraMap_apply F W.CoordinateRing W.FunctionField]
    exact le_of_eq (valuation_algebraMap_isUnit v
      ((isUnit_iff_ne_zero.2 hc).map (algebraMap F W.CoordinateRing)))

/-- ★★★★★**場合 B では `1 < w(μ x)`**——整拡大から。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★`w(μ x) ≤ 1` なら `F[x]` の元の付値はすべて 1 以下になり、
`r` のモニック方程式に `valuation_eval₂_of_monic` を当てると
`w(0) = w(μ r)^{deg} ≠ 0` となって矛盾する。
★★**これが「無限遠点の上にいる」ことの初等的な言い換え**である。 -/
theorem one_lt_valuation_genX (v : HeightOneSpectrum W.CoordinateRing)
    (μ : W.CoordinateRing →+* W.FunctionField)
    (hμF : ∀ c : F, μ (algebraMap F W.CoordinateRing c) = algebraMap F W.FunctionField c)
    {r : W.CoordinateRing} (hr : 1 < v.valuation W.FunctionField (μ r)) :
    1 < v.valuation W.FunctionField (μ (genX W)) := by
  set w := v.valuation W.FunctionField with hw
  by_contra hx
  rw [not_lt] at hx
  have hσ : ∀ q : Polynomial F,
      w (μ (algebraMap (Polynomial F) W.CoordinateRing q)) ≤ 1 := by
    intro q
    rw [← eval₂_genX_eq_algebraMap, Polynomial.hom_eval₂]
    refine valuation_eval₂_le_one w _ (fun c => ?_) hx q
    rw [RingHom.comp_apply, hμF]
    exact valuation_algebraMap_field W v c
  obtain ⟨p, hpm, hpe⟩ := Algebra.IsIntegral.isIntegral (R := Polynomial F) r
  have h0 : μ (Polynomial.eval₂ (algebraMap (Polynomial F) W.CoordinateRing) r p) = 0 := by
    rw [hpe, map_zero]
  rw [Polynomial.hom_eval₂] at h0
  have hval := valuation_eval₂_of_monic w (μ.comp (algebraMap (Polynomial F) W.CoordinateRing))
    (fun q => by rw [RingHom.comp_apply]; exact hσ q) hr hpm
  rw [h0, Valuation.map_zero] at hval
  exact absurd hval.symm (ne_of_gt (pow_pos (lt_trans zero_lt_one hr) p.natDegree))

/-! ## ★出典の紐付け(`.src`) -/

def one_lt_valuation_genX.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——引き戻しが極を持つ素点では x の付値が 1 より大きいこと)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
