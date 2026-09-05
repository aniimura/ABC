import ABC3.Found.PGC.PrincipalUnitsRank

/-!
# `K^× ≅ ℤ × 𝒪_K^×`——付値による分裂

素元 `ϖ`(`𝒪_K` の既約元)を一つ固定すると

```
Multiplicative ℤ × 𝒪_K^×  ≃*  K^×,   (n, u) ↦ ϖ^n · u
```

が群同型になる。全射性は mathlib の
`IsDiscreteValuationRing.exists_units_eq_smul_zpow_of_irreducible`、
単射性は `0 < ‖ϖ‖ < 1` から `‖ϖ‖^n = 1 ⟹ n = 0`。

## なぜ要るか

古典的局所類体論の相互律 `Γ_K^ab ≅ (K^×)^∧` の右辺を分解する第一歩。
`Found/Teichmuller.lean`(`𝒪_K^× ≅ 𝓀^× × (1+𝔪_K)`)、
`Found/PGC/PrincipalUnitsLog.lean`(`1+𝔪_K ≃* 加法群`)、
`Found/PGC/PrincipalUnitsRank.lean`(その `ℤ_p`-階数 = `[K:ℚ_p]`)と
合わせると

```
K^×  ≅  ℤ × 𝓀^× × (1+𝔪_K),   1+𝔪_K は ℤ_p 上階数 [K:ℚ_p]
```

——[pGC] Proposition 1.2 が `Γ_K^ab` の群構造から `q` と `[K:ℚ_p]` を
読み取るときの、右辺側の完全な分解。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open scoped NormedField Valued

variable {p : ℕ} [Fact p.Prime]

theorem algebraMap_irreducible_ne_zero (K : PAdicLocalField p) {ϖ : 𝒪[K.carrier]}
    (hϖ : Irreducible ϖ) : algebraMap 𝒪[K.carrier] K.carrier ϖ ≠ 0 := by
  intro hcon
  exact hϖ.ne_zero (Subtype.ext hcon)

/-- `𝒪_K` の単数のノルムは `1`。 -/
theorem norm_unit_carrierIntegers (K : PAdicLocalField p) (u : (𝒪[K.carrier])ˣ) :
    ‖((u : 𝒪[K.carrier]) : K.carrier)‖ = 1 := by
  have h1 : ((u : 𝒪[K.carrier]) : K.carrier) * ((u⁻¹ : (𝒪[K.carrier])ˣ) : K.carrier) = 1 :=
    congrArg (fun z : 𝒪[K.carrier] => (z : K.carrier)) u.mul_inv
  have hu : ‖((u : 𝒪[K.carrier]) : K.carrier)‖ ≤ 1 := by
    have hv : Valued.v ((u : 𝒪[K.carrier]) : K.carrier)
        = (‖((u : 𝒪[K.carrier]) : K.carrier)‖₊ : NNReal) := NNReal.eq rfl
    have h2 : Valued.v ((u : 𝒪[K.carrier]) : K.carrier) ≤ 1 := (u : 𝒪[K.carrier]).2
    rw [hv] at h2
    exact_mod_cast h2
  have hiu : ‖((u⁻¹ : (𝒪[K.carrier])ˣ) : K.carrier)‖ ≤ 1 := by
    have hv : Valued.v ((u⁻¹ : (𝒪[K.carrier])ˣ) : K.carrier)
        = (‖((u⁻¹ : (𝒪[K.carrier])ˣ) : K.carrier)‖₊ : NNReal) := NNReal.eq rfl
    have h2 : Valued.v ((u⁻¹ : (𝒪[K.carrier])ˣ) : K.carrier) ≤ 1 :=
      ((u⁻¹ : (𝒪[K.carrier])ˣ) : 𝒪[K.carrier]).2
    rw [hv] at h2
    exact_mod_cast h2
  have hprod := congrArg norm h1
  rw [norm_mul, norm_one] at hprod
  nlinarith [norm_nonneg ((u : 𝒪[K.carrier]) : K.carrier),
    norm_nonneg ((u⁻¹ : (𝒪[K.carrier])ˣ) : K.carrier)]

/-- `(n, u) ↦ ϖ^n · u`。 -/
noncomputable def unitsSplitHom (K : PAdicLocalField p) {ϖ : 𝒪[K.carrier]}
    (hϖ : Irreducible ϖ) : Multiplicative ℤ × (𝒪[K.carrier])ˣ →* (K.carrier)ˣ where
  toFun q := (Units.mk0 (algebraMap 𝒪[K.carrier] K.carrier ϖ)
      (algebraMap_irreducible_ne_zero K hϖ)) ^ (Multiplicative.toAdd q.1)
    * Units.map (algebraMap 𝒪[K.carrier] K.carrier : 𝒪[K.carrier] →+* K.carrier).toMonoidHom q.2
  map_one' := by simp
  map_mul' a b := by
    simp only [Prod.fst_mul, Prod.snd_mul, map_mul, toAdd_mul, zpow_add]
    simp only [mul_assoc, mul_comm, mul_left_comm]

theorem unitsSplitHom_val (K : PAdicLocalField p) {ϖ : 𝒪[K.carrier]} (hϖ : Irreducible ϖ)
    (n : ℤ) (u : (𝒪[K.carrier])ˣ) :
    ((unitsSplitHom K hϖ (Multiplicative.ofAdd n, u) : (K.carrier)ˣ) : K.carrier)
      = (algebraMap 𝒪[K.carrier] K.carrier ϖ) ^ n * ((u : 𝒪[K.carrier]) : K.carrier) := by
  simp only [unitsSplitHom, MonoidHom.coe_mk, OneHom.coe_mk, Units.val_mul,
    Units.val_zpow_eq_zpow_val, Units.val_mk0, toAdd_ofAdd, Units.coe_map,
    RingHom.toMonoidHom_eq_coe, MonoidHom.coe_coe]
  rfl

theorem unitsSplitHom_surjective (K : PAdicLocalField p) {ϖ : 𝒪[K.carrier]}
    (hϖ : Irreducible ϖ) : Function.Surjective (unitsSplitHom K hϖ) := by
  haveI := valuationRing_isDVR K
  haveI : IsFractionRing 𝒪[K.carrier] K.carrier :=
    ValuationRing.instIsFractionRingInteger (K := K.carrier) Valued.v
  intro x
  obtain ⟨n, u, hnu⟩ := IsDiscreteValuationRing.exists_units_eq_smul_zpow_of_irreducible
    (R := 𝒪[K.carrier]) (K := K.carrier) hϖ (x := (x : K.carrier)) x.ne_zero
  refine ⟨(Multiplicative.ofAdd n, u), ?_⟩
  apply Units.ext
  rw [unitsSplitHom_val K hϖ n u, hnu]
  rw [show (u • (algebraMap 𝒪[K.carrier] K.carrier ϖ) ^ n)
      = ((u : 𝒪[K.carrier]) : K.carrier) * (algebraMap 𝒪[K.carrier] K.carrier ϖ) ^ n by
    rw [Units.smul_def, Algebra.smul_def]; rfl]
  ring

theorem unitsSplitHom_injective (K : PAdicLocalField p) {ϖ : 𝒪[K.carrier]}
    (hϖ : Irreducible ϖ) : Function.Injective (unitsSplitHom K hϖ) := by
  haveI := valuationRing_isDVR K
  rw [injective_iff_map_eq_one]
  rintro ⟨m, u⟩ h
  have hval := congrArg (fun z : (K.carrier)ˣ => (z : K.carrier)) h
  rw [show ((m, u) : Multiplicative ℤ × (𝒪[K.carrier])ˣ)
      = (Multiplicative.ofAdd (Multiplicative.toAdd m), u) from rfl,
    unitsSplitHom_val K hϖ (Multiplicative.toAdd m) u] at hval
  simp only [Units.val_one] at hval
  have hnorm := congrArg norm hval
  rw [norm_mul, norm_zpow, norm_unit_carrierIntegers K u, mul_one, norm_one] at hnorm
  have hpos : 0 < ‖(algebraMap 𝒪[K.carrier] K.carrier ϖ)‖ :=
    Valued.integer.norm_irreducible_pos hϖ
  have hlt : ‖(algebraMap 𝒪[K.carrier] K.carrier ϖ)‖ < 1 :=
    Valued.integer.norm_irreducible_lt_one hϖ
  have hn0 : Multiplicative.toAdd m = 0 :=
    (zpow_eq_one_iff_right₀ (le_of_lt hpos) (ne_of_lt hlt)).mp hnorm
  have hm1 : m = 1 := by
    have hmm : m = Multiplicative.ofAdd (Multiplicative.toAdd m) := rfl
    rw [hmm, hn0]; rfl
  rw [hm1] at hval
  simp only [toAdd_one, zpow_zero, one_mul] at hval
  have hu1 : u = 1 := by
    apply Units.ext
    apply Subtype.ext
    exact hval
  rw [Prod.ext_iff]
  exact ⟨hm1, hu1⟩

/-- **★★`K^× ≅ ℤ × 𝒪_K^×`**——素元 `ϖ` を一つ固定した分裂。 -/
noncomputable def unitsSplitEquiv (K : PAdicLocalField p) {ϖ : 𝒪[K.carrier]}
    (hϖ : Irreducible ϖ) : Multiplicative ℤ × (𝒪[K.carrier])ˣ ≃* (K.carrier)ˣ :=
  MulEquiv.ofBijective (unitsSplitHom K hϖ)
    ⟨unitsSplitHom_injective K hϖ, unitsSplitHom_surjective K hϖ⟩

/-- 素元は存在する(`𝒪_K` は離散付値環)。 -/
theorem exists_unitsSplitEquiv (K : PAdicLocalField p) :
    ∃ ϖ : 𝒪[K.carrier], Irreducible ϖ
      ∧ Nonempty (Multiplicative ℤ × (𝒪[K.carrier])ˣ ≃* (K.carrier)ˣ) := by
  haveI := valuationRing_isDVR K
  obtain ⟨ϖ, hϖ⟩ := IsDiscreteValuationRing.exists_irreducible 𝒪[K.carrier]
  exact ⟨ϖ, hϖ, ⟨unitsSplitEquiv K hϖ⟩⟩

end ABC3.Found.PGC
