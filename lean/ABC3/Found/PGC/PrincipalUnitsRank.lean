import ABC3.Found.PGC.CarrierIntegersFree
import ABC3.Found.PGC.PrincipalUnitsLog

/-!
# 主単数群の `ℤ_p`-階数は `[K:ℚ_p]`

`Found/PGC/PrincipalUnitsLog.lean` で

```
smallPrincipalUnits K = {u ∈ K^× | ‖u-1‖ ≤ 1/4}  ≃*  Multiplicative (smallBall K)
```

(p 進対数)を得た。本ファイルはその右辺 `smallBall K = {y | ‖y‖ ≤ 1/4}` に
`ℤ_p`-加群の構造を入れ、**自由で階数が `[K:ℚ_p]`** であることを示す。

## 筋

`smallBall` は `𝒪_K` の中に挟まれる:

```
p² · 𝒪_K  ⊆  smallBall  ⊆  𝒪_K
```

(`‖p‖ ≤ 1/2` だから `‖p²x‖ ≤ 1/4`。逆は `1/4 ≤ 1`。)両側とも単射
`ℤ_p`-線型写像なので `finrank` は挟み撃ちで
`finrank ℤ_[p] 𝒪_K = [K:ℚ_p]`(`Found/PGC/CarrierIntegersFree.lean`)に
一致する。有限性は `𝒪_K` が `ℤ_p` 上 Noether 加群であることから。

これで [pGC] Proposition 1.2 が `Γ_K^ab` の群構造から読み取る二つの量
——`q`(`Found/Teichmuller.lean` の `𝒪_K^× ≅ 𝓀^× × (1+𝔪_K)` の第一因子)と
`[K:ℚ_p]`(第二因子の `ℤ_p`-階数)——の**両方**が形式化された量として
手元に揃う。

★記録: 本ファイルを書くために `Found/PGC/PadicLogMul.lean` の
`coeff_pow_eq_zero_of_lt` を `coeff_polynomial_pow_eq_zero_of_lt` へ改名した
——`Found/PGC/AdjoinIntegers.lean` に同名(こちらは `PowerSeries` 版)が
あり、**両方を import した瞬間に `environment already contains` で落ちて
いた**。Lubin-Tate 系(AdjoinIntegers)と p 進対数系(PadicLog*)の二つの
枝は、この改名まで一つのファイルから同時に使えなかった。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open scoped NormedField Valued

variable {p : ℕ} [Fact p.Prime]

theorem norm_algebraMap_padicInt_le_one (K : PAdicLocalField p) (c : ℤ_[p]) :
    ‖algebraMap ℤ_[p] K.carrier c‖ ≤ 1 := by
  show ‖algebraMap ℚ_[p] K.carrier ((c : ℚ_[p]))‖ ≤ 1
  rw [norm_algebraMap']
  exact c.2

/-- `‖p‖ ≤ 1/2`(`p ≥ 2` だから)。 -/
theorem norm_natCast_p_le_half (K : PAdicLocalField p) :
    ‖((p : ℕ) : K.carrier)‖ ≤ 1 / 2 := by
  have hp : ((p : ℕ) : K.carrier) = algebraMap ℚ_[p] K.carrier ((p : ℕ) : ℚ_[p]) := by
    push_cast; rfl
  rw [hp, norm_algebraMap', Padic.norm_p]
  have h2 : (2 : ℝ) ≤ (p : ℝ) := by exact_mod_cast (Fact.out : p.Prime).two_le
  rw [inv_le_comm₀ (by linarith) (by norm_num)]
  linarith

/-- 半径 `1/4` の球を `ℤ_p`-部分加群として見たもの(台は `smallBall K` と同じ)。 -/
def smallBallSubmodule (K : PAdicLocalField p) : Submodule ℤ_[p] K.carrier where
  carrier := {y : K.carrier | ‖y‖ ≤ 1 / 4}
  zero_mem' := by simp
  add_mem' := by
    intro a b ha hb
    simp only [Set.mem_setOf_eq] at *
    exact le_trans (IsUltrametricDist.norm_add_le_max a b) (max_le ha hb)
  smul_mem' := by
    intro c y hy
    simp only [Set.mem_setOf_eq] at *
    rw [Algebra.smul_def, norm_mul]
    calc ‖algebraMap ℤ_[p] K.carrier c‖ * ‖y‖
        ≤ 1 * (1 / 4) :=
          mul_le_mul (norm_algebraMap_padicInt_le_one K c) hy (norm_nonneg y) (by norm_num)
      _ = 1 / 4 := by ring

/-- 台は `PrincipalUnitsLog.lean` の `smallBall K` と同じ。 -/
theorem smallBallSubmodule_carrier (K : PAdicLocalField p) :
    (smallBallSubmodule K : Set K.carrier) = (smallBall K : Set K.carrier) := rfl

/-- 包含 `smallBall ↪ 𝒪_K`。 -/
noncomputable def smallBallIncl (K : PAdicLocalField p) :
    (smallBallSubmodule K) →ₗ[ℤ_[p]] 𝒪[K.carrier] where
  toFun y := ⟨(y : K.carrier), by
    rw [Valuation.mem_integer_iff]
    have hv : Valued.v (y : K.carrier) = (‖(y : K.carrier)‖₊ : NNReal) := NNReal.eq rfl
    rw [hv]
    have hy : ‖(y : K.carrier)‖ ≤ 1 / 4 := y.2
    have hy1 : ‖(y : K.carrier)‖ ≤ 1 := le_trans hy (by norm_num)
    exact_mod_cast hy1⟩
  map_add' a b := by apply Subtype.ext; rfl
  map_smul' c y := by apply Subtype.ext; rfl

/-- `p²` 倍 `𝒪_K ↪ smallBall`。 -/
noncomputable def smallBallMul (K : PAdicLocalField p) :
    𝒪[K.carrier] →ₗ[ℤ_[p]] (smallBallSubmodule K) where
  toFun x := ⟨((p : ℕ) : K.carrier) ^ 2 * (x : K.carrier), by
    show ‖((p : ℕ) : K.carrier) ^ 2 * (x : K.carrier)‖ ≤ 1 / 4
    rw [norm_mul, norm_pow]
    have hx : ‖(x : K.carrier)‖ ≤ 1 := by
      have hv : Valued.v (x : K.carrier) = (‖(x : K.carrier)‖₊ : NNReal) := NNReal.eq rfl
      have h1 : Valued.v (x : K.carrier) ≤ 1 := x.2
      rw [hv] at h1
      exact_mod_cast h1
    have hp : ‖((p : ℕ) : K.carrier)‖ ≤ 1 / 2 := norm_natCast_p_le_half K
    have hpsq : ‖((p : ℕ) : K.carrier)‖ ^ 2 ≤ 1 / 4 := by
      nlinarith [norm_nonneg ((p : ℕ) : K.carrier)]
    calc ‖((p : ℕ) : K.carrier)‖ ^ 2 * ‖(x : K.carrier)‖
        ≤ (1 / 4) * 1 := mul_le_mul hpsq hx (norm_nonneg _) (by norm_num)
      _ = 1 / 4 := by ring⟩
  map_add' a b := by apply Subtype.ext; show _ = _; push_cast; ring
  map_smul' c x := by
    apply Subtype.ext
    show ((p : ℕ) : K.carrier) ^ 2 * ((algebraMap ℤ_[p] 𝒪[K.carrier] c * x : 𝒪[K.carrier]) : K.carrier)
      = algebraMap ℤ_[p] K.carrier c * (((p : ℕ) : K.carrier) ^ 2 * (x : K.carrier))
    push_cast
    rw [show ((algebraMap ℤ_[p] 𝒪[K.carrier] c : 𝒪[K.carrier]) : K.carrier)
      = algebraMap ℤ_[p] K.carrier c from rfl]
    ring

theorem smallBallIncl_injective (K : PAdicLocalField p) :
    Function.Injective (smallBallIncl K) := by
  intro a b hab
  apply Subtype.ext
  exact congrArg (fun z : 𝒪[K.carrier] => (z : K.carrier)) hab

theorem smallBallMul_injective (K : PAdicLocalField p) :
    Function.Injective (smallBallMul K) := by
  intro a b hab
  have h : ((p : ℕ) : K.carrier) ^ 2 * (a : K.carrier)
      = ((p : ℕ) : K.carrier) ^ 2 * (b : K.carrier) :=
    congrArg (fun z : smallBallSubmodule K => (z : K.carrier)) hab
  haveI : CharZero K.carrier :=
    charZero_of_injective_algebraMap (algebraMap ℚ_[p] K.carrier).injective
  have hp0 : ((p : ℕ) : K.carrier) ≠ 0 := by
    have hpn : (p : ℕ) ≠ 0 := Nat.Prime.ne_zero Fact.out
    exact_mod_cast (Nat.cast_ne_zero (R := K.carrier)).mpr hpn
  apply Subtype.ext
  exact mul_left_cancel₀ (pow_ne_zero 2 hp0) h

/-- `smallBall` は `ℤ_p` 上有限(`𝒪_K` が Noether 加群だから)。 -/
theorem module_finite_smallBall (K : PAdicLocalField p) :
    Module.Finite ℤ_[p] (smallBallSubmodule K) := by
  haveI := module_finite_carrierIntegers K
  haveI : IsNoetherian ℤ_[p] 𝒪[K.carrier] :=
    isNoetherian_of_isNoetherianRing_of_finite ℤ_[p] 𝒪[K.carrier]
  exact Module.Finite.of_injective (smallBallIncl K) (smallBallIncl_injective K)

/-- `smallBall` は `ℤ_p` 上自由(捩れなし+有限生成、`ℤ_p` は PID)。 -/
theorem module_free_smallBall (K : PAdicLocalField p) :
    Module.Free ℤ_[p] (smallBallSubmodule K) := by
  haveI := module_finite_smallBall K
  haveI : Module.IsTorsionFree ℤ_[p] K.carrier :=
    Module.isTorsionFree_iff_algebraMap_injective.mpr (injective_algebraMap_padicInt K)
  exact Module.free_of_finite_type_torsion_free'

/-- **★★★階数は `[K:ℚ_p]`**——`p²·𝒪_K ⊆ smallBall ⊆ 𝒪_K` の挟み撃ち。
`PrincipalUnitsLog.padicLogUnitsEquiv` と合わせると、主単数群
`1+𝔪_K`(の十分小さい層)は `ℤ_p` 上の階数 `[K:ℚ_p]` の自由加群と
(乗法群として)同型。 -/
theorem finrank_smallBall (K : PAdicLocalField p) :
    Module.finrank ℤ_[p] (smallBallSubmodule K) = Module.finrank ℚ_[p] K.carrier := by
  haveI := module_finite_carrierIntegers K
  haveI := module_finite_smallBall K
  have h1 : Module.finrank ℤ_[p] (smallBallSubmodule K) ≤ Module.finrank ℤ_[p] 𝒪[K.carrier] :=
    LinearMap.finrank_le_finrank_of_injective (smallBallIncl_injective K)
  have h2 : Module.finrank ℤ_[p] 𝒪[K.carrier] ≤ Module.finrank ℤ_[p] (smallBallSubmodule K) :=
    LinearMap.finrank_le_finrank_of_injective (smallBallMul_injective K)
  have h3 := finrank_carrierIntegers K
  omega

end ABC3.Found.PGC
