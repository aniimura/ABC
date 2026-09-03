import ABC3.Found.PGC.LubinTateDivisibility
import ABC3.Found.PGC.ValuationRingDVR
import Mathlib.RingTheory.MvPowerSeries.Order

/-!
# Lubin-Tate 形式群の存在補題: 次数 `n+1` の障害方程式の可解性(`sorry` 無し)

`Found/PGC/LubinTateDivisibility.lean::residue_divides_R` で確立した可除性
(剰余項 `R_n` が一意化元 `π` で割り切れること)を使って、実際に次数 `n+1` の
障害方程式 `(π − π^{n+1})·φ = R_n` が**解ける**ことを示す——次数ごとの
再帰構成の**1ステップ分**が完成する。

部品4(`exists_step_solution_for_R`)はこれを **`Φ` が何であっても**(次数
`n` までの整合性という帰納的仮定を課さずに)、`residue_divides_R`(可除性、
`Φ` 任意)と `MvPowerSeries.homogeneousComponent` の線型性(`map` と可換、
`map_homogeneousComponent`)を組み合わせて直接与える——次数 `n+1` の斉次成分
`homogeneousComponent (n+1) (Obstruction Φ)` に対する障害方程式が解ける。

## 続報(2026-09-04): 1ステップ全体が完成した

本ファイルの「(任意の `Φ` から出発した)1ステップの可解性」に、f側の線形化
(`Found/PGC/LubinTateLinearization.lean`)・g側の線形化(`Found/PGC/
LubinTateGLinearization.lean`)を組み合わせて、`Found/PGC/
LubinTateStepAssembly.lean::exists_next_step` として**1ステップ全体**
(「障害が次数 `n` まで消えている」→「障害が次数 `n+1` まで消えている」)を
sorry 無しで完成させた。残るのは `Nat.rec` で無限次まで実際に構成し、
最終的な `F` が関数等式を exact に満たすことを示す組み立て作業のみ。
-/

namespace ABC3.Found.PGC

/-! ### 部品1: power series 全体としての可除性(係数ごとの可除性から) -/

/-- 各係数が `ϖ` で割り切れるなら、power series 全体が `ϖ` のスカラー倍として
書ける。`Classical.choice` で係数ごとに商を選ぶ——`MvPowerSeries σ A` が
定義上ただの関数 `(σ→₀ℕ)→A` であることを使う。 -/
theorem scalar_dvd_of_forall_coeff_dvd {A σ : Type*} [CommRing A] {ϖ : A}
    (φ : MvPowerSeries σ A) (h : ∀ d, ϖ ∣ MvPowerSeries.coeff d φ) :
    ∃ ψ : MvPowerSeries σ A, φ = ϖ • ψ := by
  classical
  refine ⟨(fun d => (h d).choose : MvPowerSeries σ A), ?_⟩
  ext d
  show MvPowerSeries.coeff d φ = ϖ * (fun d => (h d).choose) d
  exact (h d).choose_spec

/-- 剰余体への還元が 0 であることから、極大イデアルの生成元 `π` による
可除性(power series 全体として)を得る——`IsLocalRing.residue_eq_zero_iff` +
`Ideal.mem_span_singleton'` + 上の部品1。 -/
theorem exists_scalar_dvd_of_map_residue_eq_zero {A σ : Type*} [CommRing A] [IsLocalRing A]
    {π : A} (hπ : IsLocalRing.maximalIdeal A = Ideal.span {π})
    {R : MvPowerSeries σ A} (hR : MvPowerSeries.map (IsLocalRing.residue A) R = 0) :
    ∃ c, R = π • c := by
  apply scalar_dvd_of_forall_coeff_dvd
  intro d
  have hcoeff : IsLocalRing.residue A (MvPowerSeries.coeff d R) = 0 := by
    have := congrArg (fun x => MvPowerSeries.coeff d x) hR
    simpa [MvPowerSeries.coeff_map] using this
  rw [IsLocalRing.residue_eq_zero_iff, hπ, Ideal.mem_span_singleton'] at hcoeff
  obtain ⟨c, hc⟩ := hcoeff
  exact ⟨c, by rw [← hc, mul_comm]⟩

/-! ### 部品2: `1 − π^n` は単数(局所環で極大イデアルの元を1から引く) -/

theorem one_sub_pow_isUnit {A : Type*} [CommRing A] [IsLocalRing A] {π : A}
    (hπ : π ∈ IsLocalRing.maximalIdeal A) (n : ℕ) (hn : n ≠ 0) :
    IsUnit (1 - π ^ n) := by
  apply IsLocalRing.isUnit_one_sub_self_of_mem_nonunits
  rw [← IsLocalRing.mem_maximalIdeal]
  exact Ideal.pow_mem_of_mem _ hπ n (Nat.pos_of_ne_zero hn)

/-! ### ★★次数 `n+1` の障害方程式は解ける -/

/-- ★★**`(π − π^{n+1})·φ = R_n` の可解性**。`π − π^{n+1} = π(1−π^n)` と分解し、
`π | R_n`(`hRn`)から `R_n = π•c` を取り、`(1−π^n)` が単数(`one_sub_pow_isUnit`)
であることからその逆数を `c` に掛けるだけで解ける——`R_n` が `π` で割り切れて
さえいれば、`(π−π^{n+1})` 自身は零因子的でも構わない(実際 `π−π^{n+1}` は
`n=0` を除き非零因子だが、この証明はそれを経由しない)。 -/
theorem exists_step_solution {A : Type*} [CommRing A] [IsDomain A] [IsLocalRing A] {π : A}
    (hπ : π ∈ IsLocalRing.maximalIdeal A) {n : ℕ} (hn : n ≠ 0)
    {Rn : MvPowerSeries (Fin 2) A} (hRn : ∃ c, Rn = π • c) :
    ∃ φ : MvPowerSeries (Fin 2) A, (π - π ^ (n + 1)) • φ = Rn := by
  obtain ⟨c, hc⟩ := hRn
  obtain ⟨u, hu⟩ := one_sub_pow_isUnit hπ n hn
  refine ⟨((u⁻¹ : Aˣ) : A) • c, ?_⟩
  have hfact : π - π ^ (n + 1) = π * (1 - π ^ n) := by ring
  rw [hfact, hc, smul_smul]
  congr 1
  rw [← hu, mul_assoc, Units.mul_inv, mul_one]

/-! ### 部品3: `MvPowerSeries.map` と `homogeneousComponent` は可換 -/

theorem map_homogeneousComponent {R S σ : Type*} [Semiring R] [Semiring S] (h : R →+* S)
    (p : ℕ) (f : MvPowerSeries σ R) :
    MvPowerSeries.map h (MvPowerSeries.homogeneousComponent p f) =
      MvPowerSeries.homogeneousComponent p (MvPowerSeries.map h f) := by
  ext d
  rw [MvPowerSeries.coeff_map, MvPowerSeries.coeff_homogeneousComponent,
    MvPowerSeries.coeff_homogeneousComponent]
  split_ifs with hd
  · rw [MvPowerSeries.coeff_map]
  · exact map_zero h

/-! ### ★★部品4: 任意の `Φ` から出発して次数 `n+1` の障害方程式が解ける -/

/-- ★★**`Φ` が何であっても(次数ごとの整合性という帰納的仮定を課さずに)、
次数 `n+1` の障害方程式が解ける**。`residue_divides_R`(可除性、`Φ` 任意)を
`homogeneousComponent (n+1)` の線型性(`map_homogeneousComponent`)と組み合わせ、
`exists_step_solution` に渡すだけ。 -/
theorem exists_step_solution_for_R {A : Type*} [CommRing A] [IsLocalRing A] [IsDomain A]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField A) pp] [Fintype (IsLocalRing.ResidueField A)]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField A) = pp ^ ff)
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π})
    (g : PowerSeries A) (hg0 : PowerSeries.constantCoeff g = 0) (f : PowerSeries A)
    (hf : PowerSeries.map (IsLocalRing.residue A) f = PowerSeries.X ^ (pp ^ ff))
    (hg : PowerSeries.map (IsLocalRing.residue A) g = PowerSeries.X ^ (pp ^ ff))
    (Φ : MvPowerSeries (Fin 2) A) (hΦ0 : MvPowerSeries.constantCoeff Φ = 0) (n : ℕ) (hn : n ≠ 0) :
    ∃ φ : MvPowerSeries (Fin 2) A,
      (π - π ^ (n + 1)) • φ =
        MvPowerSeries.homogeneousComponent (n + 1)
          (MvPowerSeries.subst (fun i => PowerSeries.subst (MvPowerSeries.X i) g) Φ -
            PowerSeries.subst Φ f) := by
  have hR0 : MvPowerSeries.map (IsLocalRing.residue A)
      (MvPowerSeries.subst (fun i => PowerSeries.subst (MvPowerSeries.X i) g) Φ -
        PowerSeries.subst Φ f) = 0 :=
    residue_divides_R hq g hg0 f hf hg Φ hΦ0
  have hRcomp0 : MvPowerSeries.map (IsLocalRing.residue A)
      (MvPowerSeries.homogeneousComponent (n + 1)
        (MvPowerSeries.subst (fun i => PowerSeries.subst (MvPowerSeries.X i) g) Φ -
          PowerSeries.subst Φ f)) = 0 := by
    rw [map_homogeneousComponent, hR0, map_zero]
  obtain ⟨c, hc⟩ := exists_scalar_dvd_of_map_residue_eq_zero hπmax hRcomp0
  exact exists_step_solution (hπmax ▸ Ideal.mem_span_singleton_self π) hn ⟨c, hc⟩

end ABC3.Found.PGC
