import ABC3.Found.GaloisRep.TateDoubling
import ABC3.Found.GaloisRep.TateNonDeg
import Mathlib.RingTheory.DedekindDomain.AdicValuation
import Mathlib.RingTheory.Valuation.ValuationRing

/-!
# Galois (G6) 第 301 ブロック —— **★★★★★★★★★★★離散付値環から `TateSetup` を作る**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★★★★★★到達点

> `R` が完備離散付値環、`K` がその分数体なら **`E_q(K) ≅ Kˣ/q^ℤ`**
> (`tate_uniformization_dvr`)

★★★第 245–300 ブロックは `TateSetup`(抽象の付値と 2 つの引き戻し条件)の上で
進めてきた。**その仮定がすべて `IsDiscreteValuationRing` から出る**ことを示す。

## ★★★★★★仮定は 5 つとも標準的な事実に落ちる

| `TateSetup`/仮定 | 出どころ |
|---|---|
| `v : Kˣ →* Multiplicative ℤ` | 高さ 1 スペクトルの adic 付値(符号を反転) |
| `hmem0`(`v ≥ 0 ⟹ 環の元`) | **DVR は付値環**(`ValuationRing.isInteger_or_isInteger`) |
| `hmem`(`v > 0 ⟹ 𝔪 の元`) | 単元なら `v = 0` だから |
| `hloc`(単元か `𝔪`) | 局所環の定義 |
| `hI`(`𝔪` は素) | 極大だから |
| `hquad`(2 次方程式の根は環の元) | **整閉**(`IsIntegrallyClosed.isIntegral_iff`) |
| `hvalring`(`t` か `t⁻¹` が環の元) | **付値環**(同上) |
| `hΔ`(`Δ ≠ 0`) | `Δ = q·(単元)`(第 101)と単射性 |

★★★★**付値そのものは 2 つの不等式の橋にしか使っていない**:
`0 ≤ vAdd v x ↔ v(x) ≤ 1` と `0 < vAdd v x ↔ v(x) < 1`。

## ★★符号について

`HeightOneSpectrum.valuation` は**環の元で `≤ 1`** になる向き(乗法的)である。
`TateSetup` は**`q` で正**になる古典的な `ord` を要求するので、
`unzero` を取って**符号を反転**した(`tateDvrVal`)。

## ★在庫との関係

`ABC3.Found.Divisor.dvrSpectrum` は同じ対象だが、そちらは Scheme 論を
引き込む鎖の上にある。★`GaloisRep` を Scheme 論から独立に保つため、
**ここでは 3 行の定義を置き直した**(`tateSpectrum`)。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `tateDvrVal` | ★★★★★DVR の分数体の正規化付値 |
| `dvr_mem_of_nonneg`・`dvr_mem_max_of_pos` | ★★★★★★引き戻しの 2 条件 |
| `dvrTateSetup` | ★★★★★★★★**`TateSetup` の構成** |
| `tate_uniformization_dvr` | ★★★★★★★★★★★**`E_q(K) ≅ Kˣ/q^ℤ`(DVR)** |
-/

namespace ABC3.Found.GaloisRep

open Complex Real WeierstrassCurve WeierstrassCurve.Affine QuotientGroup IsDedekindDomain

/-! ## ★★`WithZero (Multiplicative ℤ)` と `ℤ` の橋 -/

/-- ★★`x ≤ 1` は `toAdd ≤ 0`。 -/
theorem withZero_le_one_iff_nonpos (x : WithZero (Multiplicative ℤ)) (h : x ≠ 0) :
    x ≤ 1 ↔ Multiplicative.toAdd (WithZero.unzero h) ≤ 0 := by
  conv_lhs => rw [← WithZero.coe_unzero h]
  rw [WithZero.coe_le_one]
  constructor
  · intro hx
    simpa using Multiplicative.toAdd_le.2 hx
  · intro hx
    exact Multiplicative.toAdd_le.1 (by simpa using hx)

/-- ★★`x < 1` は `toAdd < 0`。 -/
theorem withZero_lt_one_iff_neg (x : WithZero (Multiplicative ℤ)) (h : x ≠ 0) :
    x < 1 ↔ Multiplicative.toAdd (WithZero.unzero h) < 0 := by
  conv_lhs => rw [← WithZero.coe_unzero h]
  rw [← WithZero.coe_one, WithZero.coe_lt_coe]
  constructor
  · intro hx
    simpa using Multiplicative.toAdd_lt.2 hx
  · intro hx
    exact Multiplicative.toAdd_lt.1 (by simpa using hx)

/-! ## ★★★★★DVR の分数体の正規化付値 -/

/-- ★★DVR の高さ 1 スペクトル——極大イデアルただ 1 つ。

★`ABC3.Found.Divisor.dvrSpectrum` と同じ対象だが、そちらは Scheme 論の鎖の上にある。 -/
noncomputable def tateSpectrum (R : Type) [CommRing R] [IsDomain R]
    [IsDiscreteValuationRing R] : HeightOneSpectrum R :=
  ⟨IsLocalRing.maximalIdeal R, inferInstance, IsDiscreteValuationRing.not_a_field R⟩

section Dvr

variable {R : Type} [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
  {K : Type} [Field K] [Algebra R K] [IsFractionRing R K]

/-- ★単元の付値は `0` でない。 -/
theorem tateDvrValuation_ne_zero (x : Kˣ) : ((tateSpectrum R).valuation K) (x : K) ≠ 0 :=
  (Valuation.ne_zero_iff _).2 x.ne_zero

/-- ★★★★★**DVR の分数体の正規化付値** `Kˣ →* Multiplicative ℤ`。

★`HeightOneSpectrum.valuation` は環の元で `≤ 1` なので、**符号を反転**して
古典的な `ord`(環の元で非負)に合わせる。 -/
noncomputable def tateDvrVal (R : Type) [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
    (K : Type) [Field K] [Algebra R K] [IsFractionRing R K] : Kˣ →* Multiplicative ℤ where
  toFun x := (WithZero.unzero (tateDvrValuation_ne_zero (R := R) (K := K) x))⁻¹
  map_one' := by
    have h : WithZero.unzero (tateDvrValuation_ne_zero (R := R) (K := K) 1) = 1 := by
      apply WithZero.coe_injective
      rw [WithZero.coe_unzero]
      simp
    rw [h, inv_one]
  map_mul' x y := by
    have h : WithZero.unzero (tateDvrValuation_ne_zero (R := R) (K := K) (x * y))
        = WithZero.unzero (tateDvrValuation_ne_zero (R := R) (K := K) x)
          * WithZero.unzero (tateDvrValuation_ne_zero (R := R) (K := K) y) := by
      apply WithZero.coe_injective
      rw [WithZero.coe_unzero, WithZero.coe_mul, WithZero.coe_unzero, WithZero.coe_unzero]
      simp
    rw [h, mul_inv]

theorem vAdd_tateDvrVal (x : Kˣ) :
    vAdd (tateDvrVal R K) x
      = -Multiplicative.toAdd
          (WithZero.unzero (tateDvrValuation_ne_zero (R := R) (K := K) x)) := by
  rw [vAdd, tateDvrVal]
  simp

/-- ★★★★`0 ≤ vAdd v x` は `v(x) ≤ 1`。 -/
theorem tateDvrVal_nonneg_iff (x : Kˣ) :
    0 ≤ vAdd (tateDvrVal R K) x ↔ ((tateSpectrum R).valuation K) (x : K) ≤ 1 := by
  rw [vAdd_tateDvrVal, neg_nonneg,
    ← withZero_le_one_iff_nonpos _ (tateDvrValuation_ne_zero (R := R) (K := K) x)]

/-- ★★★★`0 < vAdd v x` は `v(x) < 1`。 -/
theorem tateDvrVal_pos_iff (x : Kˣ) :
    0 < vAdd (tateDvrVal R K) x ↔ ((tateSpectrum R).valuation K) (x : K) < 1 := by
  rw [vAdd_tateDvrVal, neg_pos,
    ← withZero_lt_one_iff_neg _ (tateDvrValuation_ne_zero (R := R) (K := K) x)]

/-! ## ★★★★★★引き戻しの 2 条件 -/

/-- ★★★★環の元の付値は非負。 -/
theorem tateDvrVal_nonneg_of_mem (t : Kˣ) (h : ∃ r : R, algebraMap R K r = (t : K)) :
    0 ≤ vAdd (tateDvrVal R K) t := by
  obtain ⟨r, hr⟩ := h
  rw [tateDvrVal_nonneg_iff, ← hr]
  exact HeightOneSpectrum.valuation_le_one _ r

/-- ★★★★`𝔪` の元の付値は正。 -/
theorem tateDvrVal_pos_of_mem_max (t : Kˣ)
    (h : ∃ r ∈ IsLocalRing.maximalIdeal R, algebraMap R K r = (t : K)) :
    0 < vAdd (tateDvrVal R K) t := by
  obtain ⟨r, hrm, hr⟩ := h
  rw [tateDvrVal_pos_iff, ← hr, HeightOneSpectrum.valuation_of_algebraMap,
    HeightOneSpectrum.intValuation_lt_one_iff_dvd]
  exact (Ideal.dvd_span_singleton).2 hrm

/-- ★★★★★★**付値が非負なら環の元**——DVR は付値環である。 -/
theorem dvr_mem_of_nonneg (x : Kˣ) (h : 0 ≤ vAdd (tateDvrVal R K) x) :
    ∃ y : R, algebraMap R K y = (x : K) := by
  rcases ValuationRing.isInteger_or_isInteger R (x : K) with ⟨y, hy⟩ | ⟨y, hy⟩
  · exact ⟨y, hy⟩
  · by_cases hu : IsUnit y
    · obtain ⟨u, rfl⟩ := hu
      refine ⟨((u⁻¹ : Rˣ) : R), ?_⟩
      have h1 : algebraMap R K ((u⁻¹ : Rˣ) : R) * algebraMap R K ((u : Rˣ) : R) = 1 := by
        rw [← map_mul]
        simp
      rw [hy] at h1
      have hx0 : (x : K) ≠ 0 := x.ne_zero
      field_simp at h1
      exact h1
    · exfalso
      have hxinv : algebraMap R K y = ((x⁻¹ : Kˣ) : K) := by
        rw [hy, Units.val_inv_eq_inv_val]
      have hpos : 0 < vAdd (tateDvrVal R K) x⁻¹ :=
        tateDvrVal_pos_of_mem_max _ ⟨y, (IsLocalRing.mem_maximalIdeal y).2 hu, hxinv⟩
      rw [vAdd_inv] at hpos
      omega

/-- ★★★★★★**付値が正なら `𝔪` の元**。 -/
theorem dvr_mem_max_of_pos (x : Kˣ) (h : 0 < vAdd (tateDvrVal R K) x) :
    ∃ y ∈ IsLocalRing.maximalIdeal R, algebraMap R K y = (x : K) := by
  obtain ⟨y, hy⟩ := dvr_mem_of_nonneg x (le_of_lt h)
  refine ⟨y, ?_, hy⟩
  by_contra hmem
  have hu : IsUnit y := by
    by_contra hnu
    exact hmem ((IsLocalRing.mem_maximalIdeal y).2 hnu)
  obtain ⟨u, rfl⟩ := hu
  have h1 : algebraMap R K ((u⁻¹ : Rˣ) : R) = ((x⁻¹ : Kˣ) : K) := by
    have h2 : algebraMap R K ((u⁻¹ : Rˣ) : R) * algebraMap R K ((u : Rˣ) : R) = 1 := by
      rw [← map_mul]
      simp
    rw [hy] at h2
    rw [Units.val_inv_eq_inv_val]
    field_simp
    exact h2
  have h3 : 0 ≤ vAdd (tateDvrVal R K) x⁻¹ := tateDvrVal_nonneg_of_mem _ ⟨((u⁻¹ : Rˣ) : R), h1⟩
  rw [vAdd_inv] at h3
  omega

/-! ## ★★★★★★★★`TateSetup` の構成 -/

/-- ★★★★★★★★**離散付値環から `TateSetup` を作る**。 -/
noncomputable def dvrTateSetup (R : Type) [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
    (K : Type) [Field K] [Algebra R K] [IsFractionRing R K]
    {q : R} (hq : q ∈ IsLocalRing.maximalIdeal R) (hq0 : q ≠ 0) :
    TateSetup R (IsLocalRing.maximalIdeal R) K where
  v := tateDvrVal R K
  Q := Units.mk0 (algebraMap R K q) ((map_ne_zero_iff _ (IsFractionRing.injective R K)).2 hq0)
  hQ := tateDvrVal_pos_of_mem_max _ ⟨q, hq, rfl⟩
  q := q
  hq := hq
  hQq := rfl
  hinj := IsFractionRing.injective R K
  hmem := dvr_mem_max_of_pos
  hmem0 := dvr_mem_of_nonneg

/-! ## ★★★★残りの仮定 -/

/-- ★★局所環では単元か極大イデアルの元。 -/
theorem dvr_hloc (x : R) : IsUnit x ∨ x ∈ IsLocalRing.maximalIdeal R := by
  by_cases h : IsUnit x
  · exact Or.inl h
  · exact Or.inr ((IsLocalRing.mem_maximalIdeal x).2 h)

/-- ★★★★**2 次方程式の根は環の元**——整閉性そのもの。 -/
theorem dvr_hquad (t : K) (b c : R)
    (h : t ^ 2 + algebraMap R K b * t + algebraMap R K c = 0) :
    ∃ r : R, algebraMap R K r = t := by
  refine IsIntegrallyClosed.isIntegral_iff.1 ⟨Polynomial.X ^ 2 + Polynomial.C b * Polynomial.X
    + Polynomial.C c, ?_, ?_⟩
  · monicity!
  · simpa using h

/-- ★★★★**`t` か `t⁻¹` が環の元**——付値環そのもの。 -/
theorem dvr_hvalring (t : K) (ht : t ≠ 0) :
    (∃ r : R, algebraMap R K r = t)
      ∨ (∃ r ∈ IsLocalRing.maximalIdeal R, algebraMap R K r = t⁻¹) := by
  by_cases hint : ∃ r : R, algebraMap R K r = t
  · exact Or.inl hint
  rcases ValuationRing.isInteger_or_isInteger R t with ⟨y, hy⟩ | ⟨y, hy⟩
  · exact absurd ⟨y, hy⟩ hint
  · refine Or.inr ⟨y, ?_, hy⟩
    refine (IsLocalRing.mem_maximalIdeal y).2 ?_
    intro hu
    obtain ⟨u, rfl⟩ := hu
    refine hint ⟨((u⁻¹ : Rˣ) : R), ?_⟩
    have h2 : algebraMap R K ((u⁻¹ : Rˣ) : R) * algebraMap R K ((u : Rˣ) : R) = 1 := by
      rw [← map_mul]
      simp
    rw [hy] at h2
    field_simp at h2
    exact h2

end Dvr

/-! ## ★★★★★★★★★★★`E_q(K) ≅ Kˣ/q^ℤ`(離散付値環) -/

section Uniformization

variable {R : Type} [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
  [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
  {K : Type} [Field K] [DecidableEq K] [Algebra R K] [IsFractionRing R K]

/-- ★★★★`Δ(E_q) ≠ 0`——`K` の水準で。 -/
theorem dvr_hDelta {q : R} (hq : q ∈ IsLocalRing.maximalIdeal R) (hq0 : q ≠ 0) :
    ((tateCurveAt q hq).map (algebraMap R K)).toAffine.Δ ≠ 0 := by
  show ((tateCurveAt q hq).map (algebraMap R K)).Δ ≠ 0
  rw [WeierstrassCurve.map_Δ]
  exact (map_ne_zero_iff _ (IsFractionRing.injective R K)).2
    (tateCurveAt_Delta_ne_zero hq hq0)

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★★★★★**Tate 一意化(離散付値環)** `E_q(K) ≅ Kˣ/q^ℤ`。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem tate_uniformization_dvr {q : R} (hq : q ∈ IsLocalRing.maximalIdeal R) (hq0 : q ≠ 0) :
    Nonempty (((tateCurveAt q hq).map (algebraMap R K)).toAffine.Point ≃+
      Additive (Kˣ ⧸ Subgroup.zpowers
        (Units.mk0 (algebraMap R K q)
          ((map_ne_zero_iff _ (IsFractionRing.injective R K)).2 hq0)))) :=
  tate_uniformization_all (dvrTateSetup R K hq hq0) dvr_hloc inferInstance dvr_hquad dvr_hvalring
    (fun t h => tateDvrVal_nonneg_of_mem t h) (fun t h => tateDvrVal_pos_of_mem_max t h) hq0
    (dvr_hDelta hq hq0)

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★★★★★**Tate 一意化(底変換の向き)**。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem tate_uniformization_baseChange {q : R} (hq : q ∈ IsLocalRing.maximalIdeal R)
    (hq0 : q ≠ 0) :
    Nonempty (((tateCurveAt q hq).baseChange K).toAffine.Point ≃+
      Additive (Kˣ ⧸ Subgroup.zpowers
        (Units.mk0 (algebraMap R K q)
          ((map_ne_zero_iff _ (IsFractionRing.injective R K)).2 hq0)))) :=
  tate_uniformization_dvr hq hq0

end Uniformization

/-! ## ★出典の紐付け(`.src`) -/

def tateDvrVal.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——DVR の分数体の正規化付値)",
    sectionId := "genell-def-3-3" }

def dvrTateSetup.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化——離散付値環からの TateSetup)",
    sectionId := "genell-def-3-3" }

def tate_uniformization_dvr.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化 E_q(K) ≅ Kˣ/q^ℤ)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
