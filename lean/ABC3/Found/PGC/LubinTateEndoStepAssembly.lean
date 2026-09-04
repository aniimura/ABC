import ABC3.Found.PGC.LubinTateEndoLinearization
import ABC3.Found.PGC.LubinTateStep

/-!
# `𝒪_K` 作用への拡張(1変数版)の次数ごとの1ステップの組み立て(`sorry` 無し)

節目(1)([a]_f:f∘[a]_f=[a]_f∘fを満たす1変数冪級数の存在、`ResearchPaper/
blocked-leaves.json` の `progress2026_09_04q`〜`s`)に要る3部品——可除性
(`Found/PGC/LubinTateEndoDivisibility.lean::residue_divides_R_endo`)・
f側の線形化・g側の線形化(`Found/PGC/LubinTateEndoLinearization.lean`)
——を実際に組み合わせて、「`φ` の障害 `Obstruction(φ):=f.subst(φ)−
g.subst(φ)` が次数 `≤n` の範囲で消えている」→「`φ':=φ+c•X^{n+1}`
(`c∈A` はスカラー)の障害が次数 `≤n+1` の範囲で消える」という**帰納法の
1ステップ全体**を確立する。

## 2変数版との違い: 「取り出す」操作すら要らない

2変数版(`LubinTateStepAssembly.lean::exists_next_step`)は可除性が返す解
`φ₀`(それ自体は次数 `n+1` の斉次式とは限らない)から
`homogeneousComponent (n+1) φ₀` を「取り出す」ひと手間が要った。1変数版
では次数 `n+1` の「斉次式」は単に `c•X^{n+1}`(スカラー `c` そのもの)
なので、可除性が与える `Obstruction(φ) = π•c₀`(`c₀:PowerSeries A`)から
係数 `r:=coeff(n+1)(c₀)` を直接取り出すだけで済む——`homogeneousComponent`
もその線型性の議論も一切不要。

`(π−π^{n+1})•φ = R_n` に対応する式は、ここではスカラー方程式
`c·(π−π^{n+1}) = −Obstruction(φ)(n+1)`(`π|Obstruction(φ)(n+1)`と
`1−π^n` が単数であることから解ける)になり、`exists_step_solution`
(`LubinTateStep.lean`)と全く同じ機構(`π=π(1−π^n)`への分解+単数の
逆元)がそのままスカラーレベルで使える。

## まだ無いもの

本ファイルは**1ステップ**の帰納法補題を確立する。実際に `Nat.rec` で
`φ:ℕ→PowerSeries A` の無限列を構成し、極限 `φ_∞` が関数等式
`f.subst(φ_∞)=φ_∞.subst(g)` を**exact に**満たすことを示す段は、
2変数版(`LubinTateSequence.lean`〜`LubinTateLimit.lean`)と同型だが
1変数へ作り直す必要があり、別途の課題として残る。
-/

namespace ABC3.Found.PGC

/-- ★★★**次数ごとの帰納法の1ステップ(1変数版)**。`Obstruction φ :=
f.subst(φ) − g.subst(φ)` が次数 `≤n` の範囲で消えているとき、あるスカラー
`c` が存在して、`Obstruction (φ+c•X^{n+1})` は次数 `≤n+1` の範囲で消える。
3部品(可除性・f側線形化・g側線形化)の組み合わせで、`c` を可除性が返す
商の次数 `n+1` 成分として取ることで得られる。 -/
theorem exists_next_step_endo {A : Type*} [CommRing A] [IsLocalRing A] [IsDomain A]
    {pp : ℕ} [ExpChar (IsLocalRing.ResidueField A) pp] [Fintype (IsLocalRing.ResidueField A)]
    {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField A) = pp ^ ff)
    {π : A} (hπmax : IsLocalRing.maximalIdeal A = Ideal.span {π})
    (g : PowerSeries A) (hg0 : PowerSeries.constantCoeff g = 0)
    (hg1 : PowerSeries.coeff 1 g = π)
    (hg : PowerSeries.map (IsLocalRing.residue A) g = PowerSeries.X ^ (pp ^ ff))
    (f : PowerSeries A) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (hf : PowerSeries.map (IsLocalRing.residue A) f = PowerSeries.X ^ (pp ^ ff))
    {φ : PowerSeries A} (hφ0 : PowerSeries.constantCoeff φ = 0)
    {n : ℕ} (hn : n ≠ 0)
    (hinv : ∀ k ≤ n, PowerSeries.coeff k (PowerSeries.subst φ f - PowerSeries.subst g φ) = 0) :
    ∃ c : A,
      PowerSeries.constantCoeff (φ + c • PowerSeries.X ^ (n + 1)) = 0 ∧
      ∀ k ≤ n + 1, PowerSeries.coeff k
        (PowerSeries.subst (φ + c • PowerSeries.X ^ (n + 1)) f -
          PowerSeries.subst g (φ + c • PowerSeries.X ^ (n + 1))) = 0 := by
  have hRdvd := residue_divides_R_endo hq g hg0 f hf hg φ hφ0
  obtain ⟨c₀, hc₀⟩ := exists_scalar_dvd_of_map_residue_eq_zero hπmax hRdvd
  set r := PowerSeries.coeff (n + 1) c₀ with hr_def
  have hπn_mem : π ∈ IsLocalRing.maximalIdeal A := hπmax ▸ Ideal.mem_span_singleton_self π
  obtain ⟨u, hu⟩ := one_sub_pow_isUnit hπn_mem n hn
  set c : A := -r * ((u⁻¹ : Aˣ) : A) with hc_def
  have hnew0 : PowerSeries.constantCoeff (c • (PowerSeries.X : PowerSeries A) ^ (n + 1)) = 0 := by
    rw [PowerSeries.constantCoeff_smul]
    have hX0 : PowerSeries.constantCoeff ((PowerSeries.X : PowerSeries A) ^ (n + 1)) = 0 := by
      rw [map_pow, PowerSeries.constantCoeff_X]
      exact zero_pow (Nat.succ_ne_zero n)
    rw [hX0, smul_zero]
  refine ⟨c, ?_, ?_⟩
  · rw [map_add, hφ0, zero_add]; exact hnew0
  · intro k hk
    have hforder : (1 : ℕ∞) ≤ MvPowerSeries.order (φ : PowerSeries A) :=
      MvPowerSeries.one_le_order_iff_constCoeff_eq_zero.mpr hφ0
    have hXord : MvPowerSeries.order ((PowerSeries.X : PowerSeries A) ^ (n + 1)) = ((n + 1 : ℕ) : ℕ∞) := by
      rw [← PowerSeries.order_eq_order]
      exact PowerSeries.order_eq.mpr ⟨fun i hi => by
          have hcx := PowerSeries.coeff_X_pow (R := A) i (n + 1)
          rw [if_pos (by exact_mod_cast hi)] at hcx
          rw [hcx]; exact one_ne_zero, fun i hi => by
          rw [PowerSeries.coeff_X_pow]
          rw [if_neg (by intro heq; rw [heq] at hi; exact absurd hi (lt_irrefl _))]⟩
    have hneworder : ((n + 1 : ℕ) : ℕ∞) ≤ MvPowerSeries.order (c • (PowerSeries.X : PowerSeries A) ^ (n + 1)) :=
      hXord ▸ MvPowerSeries.le_order_smul
    have hflin := coeff_subst_linearize_1var hφ0 hnew0 hforder hneworder (Nat.succ_pos n) f π hf0 hf1 k
      (by exact_mod_cast hk)
    have hfnew : PowerSeries.coeff k (c • (PowerSeries.X : PowerSeries A) ^ (n + 1)) =
        if k = n + 1 then c else 0 := by
      rw [PowerSeries.coeff_smul, PowerSeries.coeff_X_pow]
      split_ifs with h <;> simp
    have hgsplit : PowerSeries.subst g (φ + c • (PowerSeries.X : PowerSeries A) ^ (n + 1)) =
        PowerSeries.subst g φ + c • g ^ (n + 1) := by
      have hHSg : PowerSeries.HasSubst g := by
        show IsNilpotent (PowerSeries.constantCoeff g); rw [hg0]; exact IsNilpotent.zero
      rw [PowerSeries.subst_add hHSg, PowerSeries.subst_smul hHSg, PowerSeries.subst_pow hHSg,
        PowerSeries.subst_X hHSg]
    have hgnew : PowerSeries.coeff k (c • g ^ (n + 1)) = if k = n + 1 then c * π ^ (n + 1) else 0 := by
      rw [PowerSeries.coeff_smul]
      split_ifs with h
      · rw [h, coeff_pow_self_eq_pow hg0 hg1 n, smul_eq_mul]
      · rw [coeff_lt_pow_self_eq_zero hg0 n (by omega), smul_zero]
    have hflin' : PowerSeries.coeff k (PowerSeries.subst (φ + c • PowerSeries.X ^ (n + 1)) f) =
        PowerSeries.coeff k (PowerSeries.subst φ f) + π * PowerSeries.coeff k (c • PowerSeries.X ^ (n + 1)) := by
      linear_combination hflin
    rw [map_sub, hgsplit, map_add, hflin', hfnew, hgnew]
    by_cases hkn1 : k = n + 1
    · subst hkn1
      simp only [if_true]
      have hc0eq : PowerSeries.coeff (n + 1) (PowerSeries.subst φ f - PowerSeries.subst g φ)
          = PowerSeries.coeff (n + 1) (π • c₀) := by rw [hc₀]
      rw [map_sub, PowerSeries.coeff_smul] at hc0eq
      have hobs' : PowerSeries.coeff (n + 1) (PowerSeries.subst φ f) -
          PowerSeries.coeff (n + 1) (PowerSeries.subst g φ) = π * r := by
        rw [hc0eq, hr_def, smul_eq_mul]
      have hcsolve : c * (1 - π ^ n) = -r := by
        rw [hc_def, ← hu, mul_assoc, Units.inv_mul, mul_one]
      have hπeq : π - π ^ (n + 1) = π * (1 - π ^ n) := by ring
      linear_combination hobs' + c * hπeq + π * hcsolve
    · simp only [if_neg hkn1]
      have hobs0 := hinv k (by omega)
      rw [map_sub] at hobs0
      linear_combination hobs0

end ABC3.Found.PGC
