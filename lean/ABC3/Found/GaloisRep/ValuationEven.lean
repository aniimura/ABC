import ABC3.Found.GaloisRep.ValuationInfty

/-!
# Galois (G5) 第 145 ブロック —— **★★★★★場合 B で `ν(x)` は偶数**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★奇数次が効く

第 144 で場合 B の入口(`1 < w(μ x)`)を出した。★本ブロックはそこから

    w(μ z)² = w(μ x)³        ⟹        `count_v((μ x))` は偶数

を出す。★★根拠は **`deg Ψ₂Sq = 3`(奇数)**である(第 142)。
`z² = Ψ₂Sq(x)`(第 129)に付値を当てると `2ν(z) = 3ν(x)` となり、
`gcd(2,3) = 1` から `2 ∣ ν(x)`。

★★★**これが「無限遠点で分岐する」ことの初等的な形**である。
射影モデルも場所の分類も使わない——`Ψ₂Sq` の次数が奇数であることだけを使う。

## ★★★★道具の一般化

第 144 の `valuation_eval₂_of_monic`(モニック)を
**主係数の付値が 1 の場合**に広げた(`valuation_eval₂_of_leadingCoeff`)。
★`Ψ₂Sq` の主係数は `4` で、標数 ≠ 2 なら単元だから付値は 1 である。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `valuation_eval₂_of_leadingCoeff` | ★★★★主係数の付値が 1 なら `w(p(t)) = w(t)^{deg p}` |
| `valuation_genZ_sq` | ★★★★★**`w(μ z)² = w(μ x)³`** |
| `even_count_genX` | ★★★★★**場合 B では `count_v((μ x))` は偶数** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine Polynomial IsDedekindDomain nonZeroDivisors

/-! ## ★★★★主係数版 -/

section Poly

variable {A L Γ₀ : Type} [CommRing A] [Field L] [LinearOrderedCommGroupWithZero Γ₀]

/-- ★★★★**主係数の付値が 1 なら `w(p(t)) = w(t)^{deg p}`**(`w t > 1`)。

★第 144 の `valuation_eval₂_of_monic` の一般化。 -/
theorem valuation_eval₂_of_leadingCoeff (w : Valuation L Γ₀) (σ : A →+* L)
    (hσ : ∀ c : A, w (σ c) ≤ 1) {t : L} (ht : 1 < w t)
    {p : Polynomial A} (hp0 : p ≠ 0) (hlead : w (σ p.leadingCoeff) = 1) :
    w (Polynomial.eval₂ σ t p) = (w t) ^ p.natDegree := by
  classical
  have hsum : Polynomial.eval₂ σ t p = ∑ i ∈ p.support, σ (p.coeff i) * t ^ i := by
    rw [Polynomial.eval₂_eq_sum, Polynomial.sum]
  have hdeg : p.natDegree ∈ p.support := Polynomial.natDegree_mem_support_of_nonzero hp0
  have htop : w (σ (p.coeff p.natDegree) * t ^ p.natDegree) = (w t) ^ p.natDegree := by
    rw [Valuation.map_mul, Valuation.map_pow, ← Polynomial.leadingCoeff, hlead, one_mul]
  rw [hsum, Valuation.map_sum_eq_of_lt w hdeg ?_, htop]
  intro i hi
  have hi' : i ∈ p.support ∧ i ≠ p.natDegree := by simpa [Finset.mem_sdiff] using hi
  rw [htop, Valuation.map_mul, Valuation.map_pow]
  calc w (σ (p.coeff i)) * w t ^ i ≤ 1 * w t ^ i := mul_le_mul' (hσ _) le_rfl
    _ = w t ^ i := one_mul _
    _ < (w t) ^ p.natDegree := pow_lt_pow_right₀ ht
        (lt_of_le_of_ne (Polynomial.le_natDegree_of_mem_supp i hi'.1) hi'.2)

end Poly

/-! ## ★★★★★`w(μ z)² = w(μ x)³` -/

variable {F : Type} [Field F] (W : WeierstrassCurve.Affine F)
  [inst : IsDedekindDomain W.CoordinateRing]

/-- ★★★★★**`w(μ z)² = w(μ x)³`**——`z² = Ψ₂Sq(x)` と `deg Ψ₂Sq = 3` から。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★`Ψ₂Sq` の主係数は `4`(単元)なので、最高次が単独で最大になる。 -/
theorem valuation_genZ_sq (h2 : IsUnit (2 : F)) (v : HeightOneSpectrum W.CoordinateRing)
    (μ : W.CoordinateRing →+* W.FunctionField)
    (hμF : ∀ c : F, μ (algebraMap F W.CoordinateRing c) = algebraMap F W.FunctionField c)
    (hx : 1 < v.valuation W.FunctionField (μ (genX W))) :
    (v.valuation W.FunctionField (μ (genZ W))) ^ 2
      = (v.valuation W.FunctionField (μ (genX W))) ^ 3 := by
  set w := v.valuation W.FunctionField with hw
  have h4 : (4 : F) ≠ 0 := by
    have h44 : (4 : F) = 2 * 2 := by norm_num
    rw [h44]; exact mul_ne_zero h2.ne_zero h2.ne_zero
  have hnd : (Ψ₂Sq W).natDegree = 3 := natDegree_Psi2Sq W h2
  have hne0 : (Ψ₂Sq W) ≠ 0 := by
    intro h0
    have hc := coeff_Ψ₂Sq W
    rw [h0, Polynomial.coeff_zero] at hc
    exact h4 hc.symm
  have hlead : (Ψ₂Sq W).leadingCoeff = 4 := by
    show (Ψ₂Sq W).coeff (Ψ₂Sq W).natDegree = 4
    rw [hnd]; exact coeff_Ψ₂Sq W
  have hleadval : w (algebraMap F W.FunctionField (Ψ₂Sq W).leadingCoeff) = 1 := by
    rw [hlead, IsScalarTower.algebraMap_apply F W.CoordinateRing W.FunctionField]
    exact valuation_algebraMap_isUnit v
      ((isUnit_iff_ne_zero.2 h4).map (algebraMap F W.CoordinateRing))
  have hcomp : μ.comp (algebraMap F W.CoordinateRing) = algebraMap F W.FunctionField :=
    RingHom.ext hμF
  have hq := congrArg μ (genZ_sq W)
  rw [map_pow, Polynomial.hom_eval₂, hcomp] at hq
  rw [← Valuation.map_pow, hq,
    valuation_eval₂_of_leadingCoeff w _ (valuation_algebraMap_field W v) hx hne0 hleadval, hnd]

/-- ★★★★★**場合 B では `count_v((μ x))` は偶数である**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★`2·count(μ z) = 3·count(μ x)` と `gcd(2,3) = 1` から。
★★これが「無限遠点で分岐する」ことの初等的な形である。 -/
theorem even_count_genX (h2 : IsUnit (2 : F)) (v : HeightOneSpectrum W.CoordinateRing)
    (μ : W.CoordinateRing →+* W.FunctionField)
    (hμF : ∀ c : F, μ (algebraMap F W.CoordinateRing c) = algebraMap F W.FunctionField c)
    (hx : 1 < v.valuation W.FunctionField (μ (genX W))) :
    (2 : ℤ) ∣ FractionalIdeal.count W.FunctionField v
      (FractionalIdeal.spanSingleton W.CoordinateRing⁰ (μ (genX W))) := by
  set w := v.valuation W.FunctionField with hw
  have hkey := valuation_genZ_sq W h2 v μ hμF hx
  have hx0 : μ (genX W) ≠ 0 := by
    intro h0
    rw [h0, Valuation.map_zero] at hx
    exact absurd hx (by simp)
  have hz0 : μ (genZ W) ≠ 0 := by
    intro h0
    rw [h0, Valuation.map_zero, zero_pow (by norm_num)] at hkey
    have hpos : (0 : WithZero (Multiplicative ℤ)) < w (μ (genX W)) ^ 3 :=
      pow_pos (lt_trans zero_lt_one hx) 3
    exact absurd hkey.symm (ne_of_gt hpos)
  rw [valuation_eq_exp_neg_count (K := W.FunctionField) v hx0,
    valuation_eq_exp_neg_count (K := W.FunctionField) v hz0,
    ← WithZero.exp_nsmul, ← WithZero.exp_nsmul, WithZero.exp_inj] at hkey
  rw [nsmul_eq_mul, nsmul_eq_mul] at hkey
  omega

/-! ## ★出典の紐付け(`.src`) -/

def even_count_genX.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——無限遠の上では x の付値が偶数であること)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
