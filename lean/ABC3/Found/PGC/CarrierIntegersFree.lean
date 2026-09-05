import ABC3.Found.PGC.UnramifiedExtension

/-!
# `𝒪_K` は `ℤ_p` 上の階数 `[K:ℚ_p]` の自由加群

`Found/PGC/UnramifiedExtension.lean` で `𝒪[K.carrier]` が
`𝒪[K.carrier]`-拡大の整閉包であることをスペクトルノルムで示した筋を、
そのまま**底の `ℤ_[p] ⊆ ℚ_[p]`** に対して繰り返す:

```
IsIntegral ℤ_[p] y  ↔  ‖y‖ ≤ 1        (y : K.carrier)
⟹ IsIntegralClosure 𝒪[K.carrier] ℤ_[p] K.carrier
⟹ Module.Finite / Module.Free ℤ_[p] 𝒪[K.carrier]
⟹ Module.finrank ℤ_[p] 𝒪[K.carrier] = Module.finrank ℚ_[p] K.carrier = [K:ℚ_p]
```

## なぜ要るか

`Found/PGC/PrincipalUnitsLog.lean` で主単数群 `1+𝔪_K` が加法群
`{y | ‖y‖ ≤ 1/4}` と同型であることを示した。その加法群は `𝒪_K` の
(`π^k` 倍による)相似形なので、**`ℤ_p`-加群としての階数が `[K:ℚ_p]`**。
[pGC] Proposition 1.2 が `Γ_K^ab` の群構造から `[K:ℚ_p]` を読み取るのは
この階数——本ファイルはその「階数 = `[K:ℚ_p]`」を確定させる。

`q`(剰余体の位数)の側は `Found/Teichmuller.lean` の
`𝒪_K^× ≅ 𝓀^× × (1+𝔪_K)` の第一因子から読める。

★配管: `ℤ_[p]` は `{x : ℚ_[p] // ‖x‖ ≤ 1}`、`PadicInt.subring p` は
同じ台の `Subring ℚ_[p]`——両者は `rfl` で一致するが `Ring` インスタンスの
**経路が違う**ので `rw` は落ちる(`tools/lean-idioms.md` #55)。
`eval₂` の式は部分環の形で作ってから `exact` で渡す。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open scoped NormedField Valued

variable {p : ℕ} [Fact p.Prime]

/-- `ℚ_[p]` を経由した `ℤ_[p]`-代数構造。 -/
noncomputable instance algebraPadicIntCarrier (K : PAdicLocalField p) :
    Algebra ℤ_[p] K.carrier :=
  ((algebraMap ℚ_[p] K.carrier).comp (algebraMap ℤ_[p] ℚ_[p])).toAlgebra

instance isScalarTowerPadicIntCarrier (K : PAdicLocalField p) :
    IsScalarTower ℤ_[p] ℚ_[p] K.carrier :=
  IsScalarTower.of_algebraMap_eq (fun _ => rfl)

/-- **ノルム≤1 ⟹ `ℤ_[p]` 上整**——最小多項式の係数がすべてノルム≤1
(`spectralValue_le_one_iff`)だから `PadicInt.subring` に落ちる。 -/
theorem isIntegral_of_norm_le_one_base (K : PAdicLocalField p) (y : K.carrier) (hy : ‖y‖ ≤ 1) :
    IsIntegral ℤ_[p] y := by
  have hmonic : (minpoly ℚ_[p] y).Monic := minpoly.monic (IsIntegral.of_finite ℚ_[p] y)
  have hcoeff : ∀ n : ℕ, ‖(minpoly ℚ_[p] y).coeff n‖ ≤ 1 := by
    rw [← spectralValue_le_one_iff hmonic]
    have h1 : spectralNorm ℚ_[p] K.carrier y = spectralValue (minpoly ℚ_[p] y) := rfl
    rw [← h1]
    exact le_of_eq_of_le (NormedAlgebra.norm_eq_spectralNorm ℚ_[p] y).symm hy
  have hsub : ↑(minpoly ℚ_[p] y).coeffs ⊆ ((PadicInt.subring p) : Set ℚ_[p]) := by
    intro c hc
    obtain ⟨n, _, rfl⟩ := Polynomial.mem_coeffs_iff.mp hc
    exact hcoeff n
  have h3 : Polynomial.eval₂
      ((algebraMap ℚ_[p] K.carrier).comp (Subring.subtype (PadicInt.subring p))) y
      ((minpoly ℚ_[p] y).toSubring (PadicInt.subring p) hsub) = 0 := by
    rw [← Polynomial.eval₂_map, Polynomial.map_toSubring]
    exact minpoly.aeval ℚ_[p] y
  exact ⟨(minpoly ℚ_[p] y).toSubring (PadicInt.subring p) hsub,
    (Polynomial.monic_toSubring _ _ _).mpr hmonic, h3⟩

/-- **`ℤ_[p]` 上整 ⟹ ノルム≤1**——`norm_root_le_spectralValue`。 -/
theorem norm_le_one_of_isIntegral_base (K : PAdicLocalField p) (y : K.carrier)
    (hy : IsIntegral ℤ_[p] y) : ‖y‖ ≤ 1 := by
  obtain ⟨q, hqm, hq⟩ := hy
  have hmap : (q.map (Subring.subtype (PadicInt.subring p))).Monic := hqm.map _
  have hroot : Polynomial.aeval y (q.map (Subring.subtype (PadicInt.subring p))) = 0 := by
    rw [Polynomial.aeval_def, Polynomial.eval₂_map]
    exact hq
  have hcoeff : ∀ n : ℕ, ‖(q.map (Subring.subtype (PadicInt.subring p))).coeff n‖ ≤ 1 := by
    intro n
    rw [Polynomial.coeff_map]
    exact (q.coeff n).2
  have hsv : spectralValue (q.map (Subring.subtype (PadicInt.subring p))) ≤ 1 :=
    (spectralValue_le_one_iff hmap).mpr hcoeff
  have hle : (spectralAlgNorm ℚ_[p] K.carrier) y
      ≤ spectralValue (q.map (Subring.subtype (PadicInt.subring p))) :=
    norm_root_le_spectralValue spectralAlgNorm_isPowMul isNonarchimedean_spectralNorm hmap hroot
  rw [spectralAlgNorm_def] at hle
  calc ‖y‖ = spectralNorm ℚ_[p] K.carrier y := NormedAlgebra.norm_eq_spectralNorm ℚ_[p] y
    _ ≤ _ := hle
    _ ≤ 1 := hsv

/-- **`𝒪[K.carrier]` は `ℤ_[p]` の `K.carrier` における整閉包**。 -/
theorem isIntegralClosure_carrierIntegers (K : PAdicLocalField p) :
    IsIntegralClosure 𝒪[K.carrier] ℤ_[p] K.carrier := by
  constructor
  · exact Subtype.val_injective
  · intro z
    constructor
    · intro hz
      exact ⟨⟨z, by
        rw [Valuation.mem_integer_iff]
        have hv : Valued.v z = (‖z‖₊ : NNReal) := NNReal.eq rfl
        rw [hv]
        exact_mod_cast norm_le_one_of_isIntegral_base K z hz⟩, rfl⟩
    · rintro ⟨w, rfl⟩
      refine isIntegral_of_norm_le_one_base K _ ?_
      have hw : Valued.v (w : K.carrier) ≤ 1 := w.2
      have hv : Valued.v (w : K.carrier) = (‖(w : K.carrier)‖₊ : NNReal) := NNReal.eq rfl
      rw [hv] at hw
      exact_mod_cast hw

/-- `ℤ_[p] → 𝒪[K.carrier]`(ノルムが保たれるので像は整数環に入る)。 -/
noncomputable instance algebraPadicIntCarrierIntegers (K : PAdicLocalField p) :
    Algebra ℤ_[p] 𝒪[K.carrier] :=
  RingHom.toAlgebra
    { toFun := fun a => ⟨algebraMap ℤ_[p] K.carrier a, by
        rw [Valuation.mem_integer_iff]
        have hv : Valued.v (algebraMap ℤ_[p] K.carrier a)
            = (‖algebraMap ℤ_[p] K.carrier a‖₊ : NNReal) := NNReal.eq rfl
        rw [hv]
        have hle : ‖algebraMap ℤ_[p] K.carrier a‖ ≤ 1 := by
          show ‖algebraMap ℚ_[p] K.carrier ((a : ℚ_[p]))‖ ≤ 1
          rw [norm_algebraMap']
          exact a.2
        exact_mod_cast hle⟩
      map_one' := by apply Subtype.ext; exact map_one _
      map_mul' := fun a b => by apply Subtype.ext; exact map_mul _ a b
      map_zero' := by apply Subtype.ext; exact map_zero _
      map_add' := fun a b => by apply Subtype.ext; exact map_add _ a b }

instance isScalarTowerPadicIntCarrierIntegers (K : PAdicLocalField p) :
    IsScalarTower ℤ_[p] 𝒪[K.carrier] K.carrier :=
  IsScalarTower.of_algebraMap_eq (fun _ => rfl)

theorem injective_algebraMap_padicInt (K : PAdicLocalField p) :
    Function.Injective (algebraMap ℤ_[p] K.carrier) := by
  have h : (algebraMap ℤ_[p] K.carrier)
      = (algebraMap ℚ_[p] K.carrier).comp (algebraMap ℤ_[p] ℚ_[p]) := rfl
  rw [h, RingHom.coe_comp]
  exact Function.Injective.comp (algebraMap ℚ_[p] K.carrier).injective Subtype.val_injective

/-- **`𝒪_K` は `ℤ_[p]` 上有限**。 -/
theorem module_finite_carrierIntegers (K : PAdicLocalField p) :
    Module.Finite ℤ_[p] 𝒪[K.carrier] := by
  haveI := isIntegralClosure_carrierIntegers K
  exact IsIntegralClosure.finite ℤ_[p] ℚ_[p] K.carrier 𝒪[K.carrier]

/-- **`𝒪_K` は `ℤ_[p]` 上自由**(`ℤ_[p]` は PID)。 -/
theorem module_free_carrierIntegers (K : PAdicLocalField p) :
    Module.Free ℤ_[p] 𝒪[K.carrier] := by
  haveI := isIntegralClosure_carrierIntegers K
  haveI : Module.IsTorsionFree ℤ_[p] K.carrier :=
    Module.isTorsionFree_iff_algebraMap_injective.mpr (injective_algebraMap_padicInt K)
  exact IsIntegralClosure.module_free ℤ_[p] ℚ_[p] K.carrier 𝒪[K.carrier]

/-- **★★階数は `[K:ℚ_p]`**——`𝒪_K` は `ℤ_p` 上の階数 `[K:ℚ_p]` の自由加群。 -/
theorem finrank_carrierIntegers (K : PAdicLocalField p) :
    Module.finrank ℤ_[p] 𝒪[K.carrier] = Module.finrank ℚ_[p] K.carrier := by
  haveI := isIntegralClosure_carrierIntegers K
  haveI : Module.IsTorsionFree ℤ_[p] K.carrier :=
    Module.isTorsionFree_iff_algebraMap_injective.mpr (injective_algebraMap_padicInt K)
  exact IsIntegralClosure.rank ℤ_[p] ℚ_[p] K.carrier 𝒪[K.carrier]

end ABC3.Found.PGC
