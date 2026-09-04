import ABC3.Found.PGC.LubinTateEndoDivisibility
import ABC3.Found.PGC.LubinTateGLinearization
import ABC3.Found.PGC.LubinTateLinearization

/-!
# `𝒪_K` 作用への拡張(1変数版)の両側の線形化(`sorry` 無し)

節目(1)([a]_f:f∘[a]_f=[a]_f∘fを満たす1変数冪級数の存在、`ResearchPaper/
blocked-leaves.json` の `progress2026_09_04q`〜`s`)に要る**両側の線形化**を
確立する。可除性(`Found/PGC/LubinTateEndoDivisibility.lean::residue_divides_
R_endo`)と合わせ、次数ごとの1ステップの組み立てに要る3部品のうち2つが揃う。

## f側(既存の一般補題がそのまま使える)

`Found/PGC/LubinTateLinearization.lean::coeff_subst_linearize` は既に
任意の添字型 `σ` について一般化済み(§9-1336)なので、`σ:=Unit`
(`PowerSeries A = MvPowerSeries Unit A` という `abbrev`)として**そのまま**
適用できる——`coeff_subst_linearize_1var` として1行で済んだ。

## g側(2変数版より単純になる)

2変数版の `Found/PGC/LubinTateGLinearization.lean::coeff_subst_g_linearize`
は斉次多項式の展開(`Fin.prod_univ_two`)を要したが、1変数版で本当に要るのは
「`g≡πX(mod deg2)` のとき `g^{n+1}` の次数 `n+1` の係数が `π^{n+1}`」という
**より単純な事実**(`coeff_pow_self_eq_pow`)だけ——`order_pow_sub_pow_ge'`
(§9-1336 で一般化済み)を `a:=g, a':=π•X` として直接適用するだけで出る。

## 技術的な発見: `PowerSeries.order` と `MvPowerSeries.order` の橋渡し

当初「`PowerSeries.order` と `MvPowerSeries.order` は defeq」と見積もって
いたが、実測(`#print`)したところ誤りだった——`PowerSeries.order` は
`Nat.find` を使う独立した定義であり、`MvPowerSeries.order` とは defeq でない。
mathlib に既に橋渡し補題 `PowerSeries.order_eq_order : φ.order =
MvPowerSeries.order φ` があり(`RingTheory/PowerSeries/Order.lean`)、これを
経由すれば `order_pow_sub_pow_ge'` 等の一般補題がそのまま使える
(`ResearchPaper/blocked-leaves.json` の `progress2026_09_04s` で記録した
「2択」のうち (i) が実際には既存の1補題で解決した)。

## まだ無いもの

両側の線形化(f側・g側)のみを確立する。実際に次数ごとの1ステップを組み立て
(可除性から得た `π | R_n(n+1)` を線形化と組み合わせて `φ_{n+1}:=φ_n+c•X^{n+1}`
の `c` を解く)・`Nat.rec` で近似列を構成し・極限が関数等式を exact に満たす
ことを示す段は、2変数版(`LubinTateStepAssembly.lean`〜`LubinTateLimit.lean`)
と同型だが1変数へ作り直す必要があり、別途の課題として残る。
-/

namespace ABC3.Found.PGC

/-! ### f側: 既存の一般補題を `σ:=Unit` で直接適用 -/

/-- `coeff_subst_linearize`(任意の添字型 `σ`)を `σ:=Unit` として1変数
(`PowerSeries A`)に適用したもの。 -/
theorem coeff_subst_linearize_1var {A : Type*} [CommRing A] {Φ φ : PowerSeries A}
    (hΦ0 : PowerSeries.constantCoeff Φ = 0) (hφ0 : PowerSeries.constantCoeff φ = 0)
    (hΦ : 1 ≤ MvPowerSeries.order Φ) {m : ℕ} (hφord : (m : ℕ∞) ≤ MvPowerSeries.order φ) (hm : 1 ≤ m)
    (f : PowerSeries A) (π : A) (hf0 : PowerSeries.coeff 0 f = 0) (hf1 : PowerSeries.coeff 1 f = π)
    (n : ℕ) (he : (n : ℕ∞) ≤ m) :
    PowerSeries.coeff n (PowerSeries.subst (Φ + φ) f) - PowerSeries.coeff n (PowerSeries.subst Φ f)
      = π * PowerSeries.coeff n φ :=
  coeff_subst_linearize (σ := Unit) hΦ0 hφ0 hΦ hφord hm f π hf0 hf1 (Finsupp.single () n)
    (by simpa using he)

/-! ### g側: `g^{n+1}` の最高次係数 -/

theorem order_sub_smul_X_ge_two {A : Type*} [CommRing A] {g : PowerSeries A} {π : A}
    (hg0 : PowerSeries.constantCoeff g = 0) (hg1 : PowerSeries.coeff 1 g = π) :
    (2 : ℕ∞) ≤ MvPowerSeries.order (g - π • PowerSeries.X : PowerSeries A) := by
  rw [← PowerSeries.order_eq_order]
  apply PowerSeries.nat_le_order
  intro k hk
  interval_cases k
  · rw [map_sub]
    have h0 : PowerSeries.coeff 0 g = 0 := by
      rw [PowerSeries.coeff_zero_eq_constantCoeff_apply, hg0]
    rw [h0]; simp
  · rw [map_sub, hg1]; simp [PowerSeries.coeff_one_X]

/-- ★★**g側の線形化(1変数版)**: `g≡πX(mod deg2)` のとき `g^{n+1}` の
次数 `n+1` の係数は `π^{n+1}`。2変数版の `coeff_subst_g_linearize` に
対応するが、`g^{n+1}` 自身の最高次係数という、より直接的な形。 -/
theorem coeff_pow_self_eq_pow {A : Type*} [CommRing A] {g : PowerSeries A} {π : A}
    (hg0 : PowerSeries.constantCoeff g = 0) (hg1 : PowerSeries.coeff 1 g = π) (n : ℕ) :
    PowerSeries.coeff (n + 1) (g ^ (n + 1)) = π ^ (n + 1) := by
  have hgorder : (1 : ℕ∞) ≤ MvPowerSeries.order (g : PowerSeries A) :=
    MvPowerSeries.one_le_order_iff_constCoeff_eq_zero.mpr hg0
  have hpiXorder : (1 : ℕ∞) ≤ MvPowerSeries.order (π • PowerSeries.X : PowerSeries A) := by
    apply MvPowerSeries.one_le_order_iff_constCoeff_eq_zero.mpr
    show MvPowerSeries.constantCoeff (π • PowerSeries.X : PowerSeries A) = 0
    rw [MvPowerSeries.constantCoeff_smul]
    show π • MvPowerSeries.constantCoeff (MvPowerSeries.X () : MvPowerSeries Unit A) = 0
    rw [MvPowerSeries.constantCoeff_X, smul_zero]
  have hdiff := order_sub_smul_X_ge_two hg0 hg1
  have hpow := order_pow_sub_pow_ge' (a := (g : PowerSeries A)) (a' := π • PowerSeries.X)
    hgorder hpiXorder (n + 1)
  have hbound : ((n : ℕ) : ℕ∞) + 2 ≤
      MvPowerSeries.order (g ^ (n + 1) - (π • PowerSeries.X) ^ (n + 1) : PowerSeries A) := by
    calc ((n : ℕ) : ℕ∞) + 2 = 2 + ((n + 1 - 1 : ℕ) : ℕ∞) := by push_cast; ring
      _ ≤ MvPowerSeries.order (g - π • PowerSeries.X : PowerSeries A) + ((n + 1 - 1 : ℕ) : ℕ∞) := by gcongr
      _ ≤ _ := hpow
  have hbound' : ((n : ℕ) : ℕ∞) + 2 ≤
      PowerSeries.order (g ^ (n + 1) - (π • PowerSeries.X) ^ (n + 1) : PowerSeries A) := by
    rw [PowerSeries.order_eq_order]; exact hbound
  have hlt : ((n + 1 : ℕ) : ℕ∞) < ((n : ℕ) : ℕ∞) + 2 := by
    exact_mod_cast (by omega : n + 1 < n + 2)
  have hz := PowerSeries.coeff_of_lt_order (n + 1) (lt_of_lt_of_le hlt hbound')
  rw [map_sub] at hz
  have hzz : PowerSeries.coeff (n + 1) ((π • PowerSeries.X : PowerSeries A) ^ (n + 1)) = π ^ (n + 1) := by
    rw [smul_pow]
    show PowerSeries.coeff (n + 1) (π ^ (n + 1) • (PowerSeries.X : PowerSeries A) ^ (n + 1)) = _
    rw [PowerSeries.coeff_smul, PowerSeries.X_pow_eq, PowerSeries.coeff_monomial]
    simp
  rw [hzz] at hz
  exact sub_eq_zero.mp hz

/-- `coeff_pow_self_eq_pow` の対: `g^{n+1}` の次数 `e<n+1` の係数は0
——`g` の次数が `≥1` なので `g^{n+1}` の次数は `≥n+1`、という単純な
order 勘定だけで出る(`coeff_pow_self_eq_pow` の主張より簡単)。 -/
theorem coeff_lt_pow_self_eq_zero {A : Type*} [CommRing A] {g : PowerSeries A}
    (hg0 : PowerSeries.constantCoeff g = 0) (n : ℕ) {e : ℕ} (he : e < n + 1) :
    PowerSeries.coeff e (g ^ (n + 1)) = 0 := by
  have hgorder : (1 : ℕ∞) ≤ MvPowerSeries.order (g : PowerSeries A) :=
    MvPowerSeries.one_le_order_iff_constCoeff_eq_zero.mpr hg0
  have hpoworder : ((n + 1 : ℕ) : ℕ∞) ≤ MvPowerSeries.order (g ^ (n + 1) : PowerSeries A) := by
    calc ((n + 1 : ℕ) : ℕ∞) = (n + 1) • (1 : ℕ∞) := by simp
      _ ≤ (n + 1) • MvPowerSeries.order (g : PowerSeries A) := by gcongr
      _ ≤ _ := MvPowerSeries.le_order_pow (n + 1)
  have hpoworder' : ((n + 1 : ℕ) : ℕ∞) ≤ PowerSeries.order (g ^ (n + 1) : PowerSeries A) := by
    rw [PowerSeries.order_eq_order]; exact hpoworder
  exact PowerSeries.coeff_of_lt_order e (lt_of_lt_of_le (by exact_mod_cast he) hpoworder')

end ABC3.Found.PGC
