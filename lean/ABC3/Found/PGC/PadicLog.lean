import ABC3.Found.PGC.LocalFieldNorm
import ABC3.Found.PGC.QpResidueField

/-!
# p 進対数(級数の収束のみ、`sorry` 無し)

原文が §1(Proposition 1.2)・§2(Proposition 2.1・2.2)で境界外入力として使う
「p進対数」——`U_K`(またはより一般に `𝔪_K`)上で
`log(1+x) ≝ Σ_{n≥1} (-1)^{n+1} xⁿ/n` が収束し、`(1+𝔪_K, ×) → (K, +)` の
準同型を与える、という古典的事実。

## ★到達点(2026-09-04): 級数の**収束**のみ確立した

mathlib 実測(2026-08-14・2026-09-04 再測定)で確認済みのとおり、mathlib には
`Real.log`/`Complex.log` の p 進版(`PadicInt.log` 相当)が**存在しない**。
`RingTheory/PowerSeries/Log.lean` に**形式的**冪級数としての `PowerSeries.log`
はあるが、これを「収束する評価」に繋ぐ経路(`MvPowerSeries.HasSubst`/`HasEval`)
は未接続——本ファイルはその形式的冪級数の道を使わず、**直接** `‖x‖ < 1` での
項別収束を示す初等的な道を取った(逸脱として記録: mathlib の形式的冪級数
インフラよりも、非アルキメデス的な「項が 0 に収束すれば和が存在する」という
`NonarchimedeanGroup` 系の道具のほうが短かった)。

**まだ無いもの(★これが `.needs` に残る本体)**: 準同型性
`log((1+x)(1+y)) = log(1+x) + log(1+y)`。これが無いと「対数」ではなく
「収束する冪級数」に留まる。§1 Proposition 1.2・§2 Proposition 2.1・2.2 が
本当に必要とするのはこの準同型性(と、それが単射になる範囲の特定)である。

## 鍵となった mathlib の道具(実測、2026-09-04)

- `Summable f ↔ Tendsto f cofinite (𝓝 1)`(`NonarchimedeanGroup.multipliable_iff_tendsto_cofinite_one`、
  `Topology/Algebra/InfiniteSum/Nonarchimedean.lean`)は**乗法的**版のみで、加法版は
  無かった——`summable_of_tendsto_cofinite_zero` として自分で汎用の形(任意の完備
  非アルキメデス的加法ノルム群)に仕立て直した。
- `Finset.Nonempty.norm_sum_le_sup'_norm`(超距離不等式の有限和版)。
- `Padic.norm_eq_zpow_neg_valuation` + `Padic.valuation_natCast` + `pow_padicValNat_dvd`
  から `‖(n:ℚ_[p])‖⁻¹ ≤ n` を出し、`‖xⁿ/n‖ ≤ n·‖x‖ⁿ` という評価に落とす。
  この評価さえあれば、`summable_pow_mul_geometric_of_norm_lt_one` の実解析側の
  結果(`‖x‖ⁿ` が幾何級数的に減衰する)と合わせて、ちょうど押し出し(squeeze)で
  項別収束が出る——**局所類体論も高次分岐群も使わない**、純粋な解析の議論。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open scoped ABC3.Found.PGC NormedField Valued

/-- **一般形**: 完備な非アルキメデス的加法ノルム群では、`cofinite` 上で 0 に収束する族は
和を持つ。乗法版(`NonarchimedeanGroup.multipliable_iff_tendsto_cofinite_one`)の
加法版が mathlib に無かったので、超距離不等式(`Finset.Nonempty.norm_sum_le_sup'_norm`)
から直接組んだ。 -/
theorem summable_of_tendsto_cofinite_zero {E : Type*} [SeminormedAddCommGroup E]
    [IsUltrametricDist E] [CompleteSpace E] {ι : Type*} {f : ι → E}
    (hf : Filter.Tendsto f Filter.cofinite (nhds 0)) : Summable f := by
  rw [summable_iff_cauchySeq_finset, cauchySeq_finset_iff_vanishing_norm]
  intro ε hε
  have hnorm : Filter.Tendsto (‖f ·‖) Filter.cofinite (nhds 0) := by
    have := hf.norm; simpa using this
  have hev : ∀ᶠ i in Filter.cofinite, ‖f i‖ < ε :=
    (Metric.tendsto_nhds.mp hnorm ε hε).mono (fun i hi => by simpa using hi)
  rw [Filter.eventually_cofinite] at hev
  refine ⟨hev.toFinset, fun t ht => ?_⟩
  have hmem : ∀ i ∈ t, ‖f i‖ < ε := by
    intro i hi
    by_contra hc
    exact Finset.disjoint_left.mp ht hi (hev.mem_toFinset.mpr hc)
  rcases t.eq_empty_or_nonempty with rfl | hne
  · simpa using hε
  · calc ‖∑ i ∈ t, f i‖ ≤ t.sup' hne (‖f ·‖) := Finset.Nonempty.norm_sum_le_sup'_norm hne f
      _ < ε := by rw [Finset.sup'_lt_iff]; exact hmem

variable {p : ℕ} [Fact p.Prime]

/-- log 級数の第 n 項: `n = 0` では 0(定数項)、それ以外は `(-1)^{n+1} xⁿ/n`。 -/
noncomputable def logTerm (K : PAdicLocalField p) (x : K.carrier) (n : ℕ) : K.carrier :=
  if n = 0 then 0 else (-1 : K.carrier) ^ (n + 1) * x ^ n / (n : K.carrier)

/-- 第 n 項の評価: `‖xⁿ/n‖ ≤ n·‖x‖ⁿ`。`‖(n:K)‖⁻¹ ≤ n` が核心(`ℚ_[p]` へ落として
`Padic.norm_eq_zpow_neg_valuation` と `p^{v_p(n)} ∣ n` から出す)。 -/
theorem logTerm_norm_le (K : PAdicLocalField p) {x : K.carrier} (n : ℕ) :
    ‖logTerm K x n‖ ≤ (n : ℝ) * ‖x‖ ^ n := by
  rcases eq_or_ne n 0 with rfl | hn
  · simp [logTerm]
  · unfold logTerm
    rw [if_neg hn]
    have hnum : ‖(-1 : K.carrier) ^ (n + 1) * x ^ n‖ = ‖x‖ ^ n := by simp
    have hnormK : ‖(n : K.carrier)‖ = ‖(n : ℚ_[p])‖ := by
      rw [show ((n : ℕ) : K.carrier) = algebraMap ℚ_[p] K.carrier ((n : ℕ) : ℚ_[p]) from
            (map_natCast _ _).symm, ABC3.Found.PGC.norm_algebraMap]
    have hinv : ‖(n : ℚ_[p])‖⁻¹ ≤ (n : ℝ) := by
      have h1 := Padic.norm_eq_zpow_neg_valuation (p := p) (x := (n : ℚ_[p])) (by exact_mod_cast hn)
      rw [Padic.valuation_natCast] at h1
      rw [h1, zpow_neg, inv_inv, zpow_natCast]
      have hdvd : p ^ padicValNat p n ∣ n := pow_padicValNat_dvd
      exact_mod_cast Nat.le_of_dvd (Nat.pos_of_ne_zero hn) hdvd
    rw [norm_div, hnum, hnormK]
    calc ‖x‖ ^ n * ‖(n : ℚ_[p])‖⁻¹ ≤ ‖x‖ ^ n * (n : ℝ) := mul_le_mul_of_nonneg_left hinv (by positivity)
      _ = (n : ℝ) * ‖x‖ ^ n := by ring

/-- **log 級数は `‖x‖ < 1` で収束する**。`sorry` 無し。 -/
theorem logTerm_summable (K : PAdicLocalField p) {x : K.carrier} (hx : ‖x‖ < 1) :
    Summable (logTerm K x) := by
  apply summable_of_tendsto_cofinite_zero
  rw [Nat.cofinite_eq_atTop]
  have hxabs : ‖(‖x‖ : ℝ)‖ < 1 := by rwa [Real.norm_eq_abs, abs_of_nonneg (norm_nonneg x)]
  have hs : Summable (fun n : ℕ => (n : ℝ) * ‖x‖ ^ n) := by
    have := summable_pow_mul_geometric_of_norm_lt_one (R := ℝ) 1 hxabs
    simpa using this
  have hgeom : Filter.Tendsto (fun n : ℕ => (n : ℝ) * ‖x‖ ^ n) Filter.atTop (nhds 0) :=
    hs.tendsto_atTop_zero
  have hnorm0 : Filter.Tendsto (fun n : ℕ => ‖logTerm K x n‖) Filter.atTop (nhds 0) :=
    squeeze_zero (fun n => norm_nonneg _) (logTerm_norm_le K) hgeom
  exact tendsto_zero_iff_norm_tendsto_zero.mpr hnorm0

/-- **p 進対数**(収束級数としての値)。`‖x‖ ≥ 1` では(級数が発散しうるので)`0` を返す
——原文が使うのは `‖x‖ < 1` の範囲だけなので、この場合分けは実害がない。 -/
noncomputable def padicLog (K : PAdicLocalField p) (x : K.carrier) : K.carrier :=
  ∑' n, logTerm K x n

/-- `padicLog` は実際に `logTerm` の和になっている(`‖x‖ < 1` のとき)。 -/
theorem hasSum_padicLog (K : PAdicLocalField p) {x : K.carrier} (hx : ‖x‖ < 1) :
    HasSum (logTerm K x) (padicLog K x) :=
  (logTerm_summable K hx).hasSum

end ABC3.Found.PGC
