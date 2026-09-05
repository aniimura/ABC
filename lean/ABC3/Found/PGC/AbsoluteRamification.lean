import ABC3.Found.PGC.UnitsSplit

/-!
# `K/ℚ_p` の分岐: `e·f = [K:ℚ_p]` と `q = p^f`

`Found/PGC/UnramifiedExtension.lean` で拡大 `K(x)/K` について示した
`e·f = [K(x):K]` を、**底の `K/ℚ_p`** について繰り返す。材料は
`Found/PGC/CarrierIntegersFree.lean`(`𝒪_K` は `ℤ_p` 上有限)。

```
e(K/ℚ_p) · f(K/ℚ_p) = [K:ℚ_p],    q = |𝓀_K| = p^{f(K/ℚ_p)}
```

## なぜ要るか

[pGC] Proposition 1.2 は `Γ_K` の位相群としての構造から `q` と `[K:ℚ_p]`
の**両方**が決まると主張する。`Found/PGC/UnitsSplit.lean` までで

```
K^×  ≅  ℤ × 𝓀^× × (1+𝔪_K)
```

と分解し、第二因子の位数 `q-1` から `q` が、第三因子の `ℤ_p`-階数から
`[K:ℚ_p]` が読めるところまで来た。本ファイルは、その二つの量が
`e·f = [K:ℚ_p]`・`q = p^f` という古典的関係で結ばれていることを確定させる
——`q` と `[K:ℚ_p]` は独立な量ではなく、`(e,f)` の組で決まる。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open scoped NormedField Valued

variable {p : ℕ} [Fact p.Prime]

/-- `ℤ_[p] → 𝒪_K` は局所準同型(ノルムが保たれるから)。 -/
instance isLocalHom_algebraMap_padicInt (K : PAdicLocalField p) :
    IsLocalHom (algebraMap ℤ_[p] 𝒪[K.carrier]) := by
  constructor
  intro a ha
  by_contra hcon
  have hlt : ‖a‖ < 1 := by
    rcases lt_or_eq_of_le a.2 with h | h
    · exact h
    · exact absurd (PadicInt.isUnit_iff.mpr h) hcon
  refine absurd ha ((Valuation.Integer.not_isUnit_iff_valuation_lt_one
    (v := (Valued.v : Valuation K.carrier NNReal))
    (x := (algebraMap ℤ_[p] 𝒪[K.carrier] a))).mpr ?_)
  show Valued.v (algebraMap ℤ_[p] K.carrier a) < 1
  have hv : Valued.v (algebraMap ℤ_[p] K.carrier a)
      = (‖algebraMap ℤ_[p] K.carrier a‖₊ : NNReal) := NNReal.eq rfl
  rw [hv]
  have hn : ‖algebraMap ℤ_[p] K.carrier a‖ < 1 := by
    show ‖algebraMap ℚ_[p] K.carrier ((a : ℚ_[p]))‖ < 1
    rw [norm_algebraMap']
    exact hlt
  exact_mod_cast hn

/-- **絶対分岐指数** `e(K/ℚ_p)`。 -/
noncomputable def absoluteRamificationIndex (K : PAdicLocalField p) : ℕ :=
  (IsLocalRing.maximalIdeal ℤ_[p]).ramificationIdx (IsLocalRing.maximalIdeal 𝒪[K.carrier])

/-- **絶対慣性次数** `f(K/ℚ_p)`。 -/
noncomputable def absoluteInertiaDegree (K : PAdicLocalField p) : ℕ :=
  (IsLocalRing.maximalIdeal ℤ_[p]).inertiaDeg (IsLocalRing.maximalIdeal 𝒪[K.carrier])

/-- **★`e·f = [K:ℚ_p]`**。 -/
theorem absoluteRamificationIndex_mul_absoluteInertiaDegree (K : PAdicLocalField p) :
    absoluteRamificationIndex K * absoluteInertiaDegree K = Module.finrank ℚ_[p] K.carrier := by
  haveI := valuationRing_isDVR K
  haveI : IsFractionRing 𝒪[K.carrier] K.carrier :=
    ValuationRing.instIsFractionRingInteger (K := K.carrier) Valued.v
  haveI := module_finite_carrierIntegers K
  exact Ideal.ramificationIdx_mul_inertiaDeg_of_isLocalRing 𝒪[K.carrier] ℚ_[p] K.carrier
    (IsDiscreteValuationRing.not_a_field ℤ_[p])

/-- `f` は剰余体の拡大次数 `[𝓀_K : 𝔽_p]`。 -/
theorem absoluteInertiaDegree_eq_finrank (K : PAdicLocalField p) :
    absoluteInertiaDegree K = Module.finrank (IsLocalRing.ResidueField ℤ_[p]) 𝓀[K.carrier] :=
  Ideal.inertiaDeg_algebraMap _ _

/-- **★`q = p^f`**——剰余体の位数は `p` の `f` 乗。 -/
theorem residueCard_eq_pow (K : PAdicLocalField p) :
    Nat.card 𝓀[K.carrier] = p ^ (absoluteInertiaDegree K) := by
  haveI : Finite 𝓀[K.carrier] := residueField_finite K
  haveI : Fintype 𝓀[K.carrier] := Fintype.ofFinite _
  haveI : Fintype (IsLocalRing.ResidueField ℤ_[p]) :=
    Fintype.ofEquiv (ZMod p) (PadicInt.residueField (p := p)).toEquiv.symm
  have hcard : Fintype.card (IsLocalRing.ResidueField ℤ_[p]) = p := by
    rw [Fintype.card_congr (PadicInt.residueField (p := p)).toEquiv]
    haveI : NeZero p := ⟨Nat.Prime.ne_zero Fact.out⟩
    exact ZMod.card p
  rw [absoluteInertiaDegree_eq_finrank K, Nat.card_eq_fintype_card,
    Module.card_eq_pow_finrank (K := IsLocalRing.ResidueField ℤ_[p]) (V := 𝓀[K.carrier]), hcard]

end ABC3.Found.PGC
