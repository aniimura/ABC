import ABC3.Found.GaloisRep.NeronWitness

/-!
# Galois (G8) 第 319 ブロック —— **★★★★★★★無限遠因子の次数 `deg∞`**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.18。

原文 (GenEll p.18):
> First, observe that if v is any local height of EL, then d · deg∞([EL]) ≥

## ★★★★★★★到達点

> **`deg∞(E) = (1/[L:ℚ]) Σ_p v_p(Δ_min) · log N(p)`** を構成し、**非負性**を示した
> (`degInfOf`・`degInfOf_nonneg`)

★★★これが G8 の `degInf` 欄の中身である。

## ★★★★★★局所高さの和ではなく**極小判別式**で書く

原文の `deg∞` は乗法還元の素点での局所高さの和だが、**局所高さ `v(q_E)` は
極小判別式の付値 `v(Δ_min)` に等しい**(第 101 の `addVal_tateDelta`:`Δ = q·(単元)`)。
★★★★そこで**すべての素点にわたる `v_p(Δ_min)` の和**として定義する:

* 乗法還元の素点では局所高さそのもの
* 良還元の素点では `0`(`Δ_min` が単元)
* 加法還元の素点では正の寄与(原文の `deg∞` は半安定を仮定するので同じ)

★★★★★**完備化を経由せずに定義できる**のが利点である——
局所高さの定義(第 310)は完備 DVR を要求するが、`v_p(Δ_min)` は `L` の素点だけで書ける。

## ★★★★極小判別式の指数

    minDeltaExp p W := v_p(Δ_W) − 12·(Néron 指数)

★`Δ_min = u⁻¹²Δ_W` だから、これはちょうど `v_p(Δ_min)` である(`minDeltaExp_eq`)。
★★**非負**(極小モデルは整だから)で、**台は有限**(第 314-315)。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `valAdd_one`・`valAdd_inv`・`valAdd_pow` | ★★`valAdd` の代数 |
| `valAdd_nonneg_iff` | ★★★付値 `≤ 1` と `valAdd ≥ 0` |
| `minDeltaExp`・`minDeltaExp_eq` | ★★★★★**極小判別式の指数** |
| `minDeltaExp_nonneg`・`minDeltaExp_finite` | ★★★★非負性と有限台 |
| `degInfOf`・`degInfOf_nonneg` | ★★★★★★★**`deg∞` の構成** |
-/

namespace ABC3.Found.GaloisRep

open IsDedekindDomain NumberField WeierstrassCurve
open scoped nonZeroDivisors Classical

variable {L : Type} [Field L] [NumberField L]

/-! ## ★★`valAdd` の代数 -/

theorem valAdd_one (p : HeightOneSpectrum (𝓞 L)) : valAdd p (1 : Lˣ) = 0 := by
  rw [valAdd_eq_zero_iff]
  show (p.valuation L) (1 : L) = 1
  exact map_one _

theorem valAdd_inv (p : HeightOneSpectrum (𝓞 L)) (x : Lˣ) : valAdd p x⁻¹ = -valAdd p x := by
  have h : valAdd p (x * x⁻¹) = valAdd p x + valAdd p x⁻¹ := valAdd_mul p x x⁻¹
  rw [mul_inv_cancel, valAdd_one] at h
  omega

theorem valAdd_pow (p : HeightOneSpectrum (𝓞 L)) (x : Lˣ) (n : ℕ) :
    valAdd p (x ^ n) = n * valAdd p x := by
  induction n with
  | zero => simp [valAdd_one]
  | succ k ih =>
    rw [pow_succ, valAdd_mul, ih]
    push_cast
    ring

/-- ★★★付値が `1` 以下であることと `valAdd` が非負であることは同値。 -/
theorem valAdd_nonneg_iff (p : HeightOneSpectrum (𝓞 L)) (x : Lˣ) :
    0 ≤ valAdd p x ↔ (p.valuation L) (x : L) ≤ 1 := by
  rw [valAdd, neg_nonneg]
  constructor
  · intro h
    rw [← WithZero.coe_unzero (valuationP_ne_zero p x), WithZero.coe_le_one]
    exact Multiplicative.toAdd_le.1 (by simpa using h)
  · intro h
    rw [← WithZero.coe_unzero (valuationP_ne_zero p x), WithZero.coe_le_one] at h
    simpa using Multiplicative.toAdd_le.2 h

/-! ## ★★★★★極小判別式の指数 -/

/-- ★★★★★**極小判別式の指数** `v_p(Δ_min)`。 -/
noncomputable def minDeltaExp (p : HeightOneSpectrum (𝓞 L)) (W : WeierstrassCurve L) : ℤ :=
  if h : W.Δ = 0 then 0 else valAdd p (Units.mk0 W.Δ h) - 12 * neronExp p W

set_option maxHeartbeats 1600000 in
/-- ★★★★これはちょうど極小モデルの判別式の付値である。 -/
theorem minDeltaExp_eq (p : HeightOneSpectrum (𝓞 L)) (W : WeierstrassCurve L) (hΔ : W.Δ ≠ 0)
    (C : VariableChange L) (h : IsMinimal (primeSubring p) (C • W)) :
    minDeltaExp p W = valAdd p (Units.mk0 ((C • W).Δ) (variableChange_Delta_ne_zero W hΔ C)) := by
  rw [minDeltaExp, dif_neg hΔ, neronExp_eq p W hΔ C h]
  have hu : (Units.mk0 ((C • W).Δ) (variableChange_Delta_ne_zero W hΔ C))
      = C.u⁻¹ ^ 12 * Units.mk0 W.Δ hΔ := by
    refine Units.ext ?_
    show (C • W).Δ = _
    rw [WeierstrassCurve.variableChange_Δ]
    push_cast
    simp
  rw [hu, valAdd_mul, valAdd_pow, valAdd_inv]
  omega

/-- ★★★★**非負**——極小モデルは整だから。 -/
theorem minDeltaExp_nonneg (p : HeightOneSpectrum (𝓞 L)) (W : WeierstrassCurve L) :
    0 ≤ minDeltaExp p W := by
  by_cases hΔ : W.Δ = 0
  · rw [minDeltaExp, dif_pos hΔ]
  · obtain ⟨C, hC⟩ := WeierstrassCurve.exists_isMinimal (primeSubring p) W
    rw [minDeltaExp_eq p W hΔ C hC, valAdd_nonneg_iff]
    haveI : WeierstrassCurve.IsIntegral (primeSubring p) (C • W) := inferInstance
    have hmem : ((C • W).Δ) ∈ primeSubring p := by
      have h1 := WeierstrassCurve.integralModel_Δ_eq (primeSubring p) (C • W)
      rw [← h1]
      exact ((integralModel (primeSubring p) (C • W)).Δ).2
    exact (mem_primeSubring_iff p _).1 hmem

set_option maxHeartbeats 1600000 in
/-- ★★★★**台は有限**。 -/
theorem minDeltaExp_finite (W : WeierstrassCurve L) (hΔ : W.Δ ≠ 0) :
    {p : HeightOneSpectrum (𝓞 L) | minDeltaExp p W ≠ 0}.Finite := by
  refine Set.Finite.subset (((IsDedekindDomain.HeightOneSpectrum.Support.finite (𝓞 L) W.Δ).union
    (IsDedekindDomain.HeightOneSpectrum.Support.finite (𝓞 L) W.Δ⁻¹)).union
    (finite_bad_primes' W hΔ)) ?_
  intro p hp
  simp only [Set.mem_setOf_eq] at hp
  by_contra hnot
  simp only [Set.mem_union, not_or, Set.mem_setOf_eq] at hnot
  obtain ⟨⟨n1, n2⟩, n3⟩ := hnot
  apply hp
  rw [minDeltaExp, dif_neg hΔ]
  have hne : neronExp p W = 0 := by
    by_contra hc
    exact n3 hc
  rw [hne]
  have hv1 : ¬ (1 < (p.valuation L) W.Δ) := n1
  have hv2 : ¬ (1 < (p.valuation L) W.Δ⁻¹) := n2
  have hz : valAdd p (Units.mk0 W.Δ hΔ) = 0 := by
    have h1 : 0 ≤ valAdd p (Units.mk0 W.Δ hΔ) := by
      rw [valAdd_nonneg_iff]
      exact not_lt.1 hv1
    have h2 : 0 ≤ valAdd p (Units.mk0 W.Δ hΔ)⁻¹ := by
      rw [valAdd_nonneg_iff]
      have hinv : (((Units.mk0 W.Δ hΔ)⁻¹ : Lˣ) : L) = (W.Δ)⁻¹ := by
        rw [Units.val_inv_eq_inv_val]
        rfl
      rw [hinv]
      exact not_lt.1 hv2
    rw [valAdd_inv] at h2
    omega
  rw [hz]
  omega

/-! ## ★★★★★★★`deg∞` の構成 -/

theorem log_absNorm_nonneg (p : HeightOneSpectrum (𝓞 L)) :
    0 ≤ Real.log (Ideal.absNorm p.asIdeal) := by
  refine Real.log_nonneg ?_
  have h : Ideal.absNorm p.asIdeal ≠ 0 := by
    rw [Ne, Ideal.absNorm_eq_zero_iff]
    exact p.ne_bot
  have h1 : 1 ≤ Ideal.absNorm p.asIdeal := Nat.one_le_iff_ne_zero.2 h
  exact_mod_cast h1

/-- ★★★★★★★**無限遠因子の次数** `deg∞(E)`。

原文 (GenEll p.18):
> First, observe that if v is any local height of EL, then d · deg∞([EL]) ≥ -/
noncomputable def degInfOf (L : Type) [Field L] [NumberField L] (W : WeierstrassCurve L) : ℝ :=
  (∑ᶠ p : HeightOneSpectrum (𝓞 L),
      (minDeltaExp p W : ℝ) * Real.log (Ideal.absNorm p.asIdeal))
    / (Module.finrank ℚ L)

/-- ★★★★★**`deg∞ ≥ 0`**——各項が非負だから。 -/
theorem degInfOf_nonneg (W : WeierstrassCurve L) : 0 ≤ degInfOf L W := by
  rw [degInfOf]
  refine div_nonneg ?_ (by positivity)
  refine finsum_nonneg (fun p => ?_)
  exact mul_nonneg (by exact_mod_cast minDeltaExp_nonneg p W) (log_absNorm_nonneg p)

/-! ## ★出典の紐付け(`.src`) -/

def minDeltaExp.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 18,
    item := "Proposition 3.4(極小判別式の指数)",
    sectionId := "genell-prop-3-4" }

def degInfOf.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 18,
    item := "Proposition 3.4(無限遠因子の次数 deg∞)",
    sectionId := "genell-prop-3-4" }

end ABC3.Found.GaloisRep
