import ABC3.Found.PGC.LubinTateDivisibility
import ABC3.Found.PGC.ValuationRingDVR

/-!
# Lubin-Tate 形式群の存在補題: 次数 `n+1` の障害方程式の可解性(`sorry` 無し)

`Found/PGC/LubinTateDivisibility.lean::residue_divides_R` で確立した可除性
(剰余項 `R_n` が一意化元 `π` で割り切れること)を使って、実際に次数 `n+1` の
障害方程式 `(π − π^{n+1})·φ = R_n` が**解ける**ことを示す——次数ごとの
再帰構成の**1ステップ分**が完成する。

## まだ無いもの

本ファイルは「1ステップの可解性」までを確立する。実際の構成
(`Φ : ℕ → MvPowerSeries (Fin 2) A` を `Nat.rec` で組み立て、各 `Φ (n+1)` が
`Φ n` に**次数 `n+1` ちょうどの斉次成分**を足したものであることを保ち、
最終的な極限 `F` が関数等式を exact に満たすことを示す)は、まだ残る——
特に「`f(Φ_n + φ_{n+1}) − (Φ_n + φ_{n+1})(g,g)` の次数 `n+2` 以上の誤差項が
消える」という段は、本ファイルでは検証していない(紙の上での議論は
`ResearchPaper/blocked-leaves.json` の `progress2026_09_04e` に記録済みだが、
Lean での形式化は別途の作業として残る)。
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

end ABC3.Found.PGC
