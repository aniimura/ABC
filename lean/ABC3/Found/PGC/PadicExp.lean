import ABC3.Found.PGC.PadicLog
import Mathlib.Topology.Algebra.InfiniteSum.Nonarchimedean

/-!
# p 進指数関数(収束 + 加法性、`sorry` 無し)

`Found/PGC/PadicLog.lean` の p 進対数は収束のみ確立していた(準同型性
`log((1+x)(1+y))=log(1+x)+log(1+y)` は未解決のまま残した、`memory/
padic-log-additivity-blocked.md` 参照)。**本ファイルは指数関数の側で
同じことを最後までやり切った**——`exp(x) ≝ Σ xⁿ/n!` が十分小さい `x` で収束し、
かつ **`exp(x+y) = exp(x)·exp(y)`(準同型性)まで sorry 無しで証明した**。

## なぜ指数のほうが対数より易しかったか

対数の準同型性 `log(uv)=log(u)+log(v)` の古典的証明は微分(`d/dx log(1+x)=1/(1+x)`)
を経由する——形式冪級数の微分と「導関数が 0 なら定数」という議論を要し、
mathlib にも p進版はおろか一般の完備非アルキメデス的環でも直接使える形が無かった
(`memory/padic-log-additivity-blocked.md` に記録した2つの近道がどちらも
別の実装コストを要した)。

**指数の加法性 `exp(x+y)=exp(x)exp(y)` は微分を経由しない**——Cauchy 積
(`Summable.tsum_mul_tsum_eq_tsum_sum_antidiagonal`)と**二項定理**
(`add_pow`)だけで閉じる、純粋に組み合わせ的な議論である:
```
exp(x)·exp(y) = Σ_{(k,l)} xᵏ/k! · yˡ/l! = Σₙ Σ_{k+l=n} xᵏyˡ/(k!l!)
             = Σₙ (1/n!) Σ_{k+l=n} C(n,k) xᵏyˡ = Σₙ (x+y)ⁿ/n! = exp(x+y)
```
中央の等号(`Σ_{k+l=n} xᵏyˡ/(k!l!) = (x+y)ⁿ/n!`)が二項定理そのもの
(`antidiagonal_sum_expTerm`)。

## 収束域(`log` より狭い)

`n!` の分母があるぶん、`log` の収束域(`‖x‖<1`)より狭い——標準的には
`v_p(x) > 1/(p-1)` を要する。ここでは実数冪(`Real.rpow`)を避けて
`‖x‖^(p-1)·p < 1` という同値な(自然数冪だけで書ける)条件を使った
(両辺を `(p-1)` 乗すれば正確に一致する境界になる——逸脱ではなく、
`Real.rpow` を経由しない書き方を選んだだけ)。

## まだ無いもの(次にここへ戻るときの入口)

`exp` と `log` が互いに逆であること(`exp(log(1+x))=1+x`・`log(exp(x))=x`)
——これが分かれば、`log(uv)=log(u)+log(v)` は `exp` の加法性から**微分無しで**
出る(`u=log(1+x), v=log(1+y)` とおき、`exp(u+v)=exp(u)exp(v)=(1+x)(1+y)`
の両辺に `log` を当てて `u+v` を回収する、という筋)。この互逆性の証明は
`exp`/`log` の係数どうしの合成公式(スターリング数的な組み合わせ論)を要し、
加法性ほど単純ではない——次の一手として記録するに留める。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open scoped ABC3.Found.PGC NormedField Valued

variable {p : ℕ} [Fact p.Prime]

/-- `K.carrier` は非アルキメデス的環——`Summable.mul_of_nonarchimedean` を使うために要る。
`NonarchimedeanAddGroup`(`IsUltrametricDist` から自動)の中身をそのまま転用する。 -/
instance padicNonarchimedeanRing (K : PAdicLocalField p) : NonarchimedeanRing K.carrier where
  is_nonarchimedean := NonarchimedeanAddGroup.is_nonarchimedean

/-- `K.carrier` は標数 0(`ℚ_[p]` からの埋め込みが単射であることから)。 -/
theorem padicCharZero (K : PAdicLocalField p) : CharZero K.carrier :=
  charZero_of_injective_algebraMap (algebraMap ℚ_[p] K.carrier).injective

/-- exp 級数の第 n 項: `xⁿ/n!`。 -/
noncomputable def expTerm (K : PAdicLocalField p) (x : K.carrier) (n : ℕ) : K.carrier :=
  x ^ n / (Nat.factorial n : K.carrier)

/-- `‖(n!:ℚ_[p])‖⁻¹ ^ (p-1) ≤ p^n`——Legendre の公式の書き換え
`(p-1)·v_p(n!) = n - s_p(n) ≤ n`(`sub_one_mul_padicValNat_factorial`)から。
`(p-1)` 乗するのは実数冪を避けるため(`p^{1/(p-1)}` を経由しない)。 -/
theorem expFactorial_norm_inv_pow_le (n : ℕ) :
    (‖(Nat.factorial n : ℚ_[p])‖⁻¹) ^ (p - 1) ≤ (p : ℝ) ^ n := by
  rcases eq_or_ne n 0 with rfl | hn
  · simp
  have h1 := Padic.norm_eq_zpow_neg_valuation (p := p) (x := (Nat.factorial n : ℚ_[p]))
    (by exact_mod_cast Nat.factorial_ne_zero n)
  rw [Padic.valuation_natCast] at h1
  rw [h1, zpow_neg, inv_inv, zpow_natCast, ← pow_mul]
  have hb : padicValNat p (Nat.factorial n) * (p - 1) ≤ n := by
    rw [mul_comm, sub_one_mul_padicValNat_factorial]; omega
  exact pow_le_pow_right₀ (by exact_mod_cast (Fact.out : p.Prime).one_lt.le) hb

/-- `‖f n‖^k → 0`(`k` 固定 `≠ 0`)ならば `‖f n‖ → 0`——`t ↦ t^k` が `ℝ≥0` 上
単調な全単射であることから、ε-N で直接示す。`Real.rpow` を経由しない。 -/
theorem tendsto_zero_of_pow_tendsto_zero {k : ℕ} (hk : k ≠ 0) {a : ℕ → ℝ} (ha : ∀ n, 0 ≤ a n)
    (h : Filter.Tendsto (fun n => (a n) ^ k) Filter.atTop (nhds 0)) :
    Filter.Tendsto a Filter.atTop (nhds 0) := by
  rw [Metric.tendsto_atTop] at h ⊢
  intro ε hε
  obtain ⟨N, hN⟩ := h (ε ^ k) (by positivity)
  refine ⟨N, fun n hn => ?_⟩
  have hthis := hN n hn
  simp only [Real.dist_eq, sub_zero] at hthis ⊢
  rw [abs_of_nonneg (ha n)]
  rw [abs_of_nonneg (pow_nonneg (ha n) k)] at hthis
  exact (pow_lt_pow_iff_left₀ (ha n) hε.le hk).mp hthis

/-- 第 n 項の評価(`(p-1)` 乗した形): `‖xⁿ/n!‖^{p-1} ≤ (‖x‖^{p-1}·p)ⁿ`。 -/
theorem expTerm_norm_pow_le (K : PAdicLocalField p) {x : K.carrier} (n : ℕ) :
    ‖expTerm K x n‖ ^ (p - 1) ≤ (‖x‖ ^ (p - 1) * (p : ℝ)) ^ n := by
  unfold expTerm
  have hnormK : ‖(Nat.factorial n : K.carrier)‖ = ‖(Nat.factorial n : ℚ_[p])‖ := by
    rw [show ((Nat.factorial n : ℕ) : K.carrier)
          = algebraMap ℚ_[p] K.carrier ((Nat.factorial n : ℕ) : ℚ_[p]) from
          (map_natCast _ _).symm, ABC3.Found.PGC.norm_algebraMap]
  have hle := expFactorial_norm_inv_pow_le (p := p) n
  rw [norm_div, norm_pow, hnormK, div_pow, div_eq_mul_inv, ← inv_pow]
  calc (‖x‖ ^ n) ^ (p - 1) * (‖(Nat.factorial n : ℚ_[p])‖⁻¹) ^ (p - 1)
      ≤ (‖x‖ ^ n) ^ (p - 1) * (p : ℝ) ^ n := mul_le_mul_of_nonneg_left hle (by positivity)
    _ = (‖x‖ ^ (p - 1) * (p : ℝ)) ^ n := by
        rw [mul_pow, ← pow_mul, ← pow_mul, mul_comm n (p - 1), mul_comm (p - 1) n]

/-- **exp 級数は `‖x‖^{p-1}·p < 1` で収束する**。`sorry` 無し。 -/
theorem expTerm_summable (K : PAdicLocalField p) {x : K.carrier}
    (hx : ‖x‖ ^ (p - 1) * (p : ℝ) < 1) : Summable (expTerm K x) := by
  apply summable_of_tendsto_cofinite_zero
  rw [Nat.cofinite_eq_atTop]
  have hc0 : 0 ≤ ‖x‖ ^ (p - 1) * (p : ℝ) := by positivity
  have hgeom : Filter.Tendsto (fun n : ℕ => (‖x‖ ^ (p - 1) * (p : ℝ)) ^ n) Filter.atTop (nhds 0) :=
    tendsto_pow_atTop_nhds_zero_of_lt_one hc0 hx
  have hpow0 : Filter.Tendsto (fun n : ℕ => ‖expTerm K x n‖ ^ (p - 1)) Filter.atTop (nhds 0) :=
    squeeze_zero (fun n => by positivity) (expTerm_norm_pow_le K) hgeom
  have hp1 : p - 1 ≠ 0 := by have := (Fact.out : p.Prime).two_le; omega
  have hnorm0 : Filter.Tendsto (fun n : ℕ => ‖expTerm K x n‖) Filter.atTop (nhds 0) :=
    tendsto_zero_of_pow_tendsto_zero hp1 (fun n => norm_nonneg _) hpow0
  exact tendsto_zero_iff_norm_tendsto_zero.mpr hnorm0

/-- **二項定理そのもの**: `Σ_{k+l=n} xᵏ/k!·yˡ/l! = (x+y)ⁿ/n!`。準同型性の核心。 -/
theorem antidiagonal_sum_expTerm (K : PAdicLocalField p) (x y : K.carrier) (n : ℕ) :
    ∑ kl ∈ Finset.antidiagonal n, expTerm K x kl.1 * expTerm K y kl.2 = expTerm K (x + y) n := by
  haveI := padicCharZero K
  unfold expTerm
  rw [Finset.Nat.sum_antidiagonal_eq_sum_range_succ
    (fun k l => x ^ k / (Nat.factorial k : K.carrier) * (y ^ l / (Nat.factorial l : K.carrier)))]
  rw [add_pow, show Finset.range n.succ = Finset.range (n + 1) from rfl,
      eq_div_iff (Nat.cast_ne_zero.mpr (Nat.factorial_ne_zero n) : (Nat.factorial n : K.carrier) ≠ 0),
      Finset.sum_mul]
  apply Finset.sum_congr rfl
  intro k hk
  simp only [Finset.mem_range, Nat.lt_succ_iff] at hk
  have hfact : n.choose k * Nat.factorial k * Nat.factorial (n - k) = Nat.factorial n :=
    Nat.choose_mul_factorial_mul_factorial hk
  have hfactK : (n.choose k : K.carrier) * (Nat.factorial k : K.carrier) * (Nat.factorial (n - k) : K.carrier)
      = (Nat.factorial n : K.carrier) := by
    have h2 := congrArg (Nat.cast : ℕ → K.carrier) hfact
    push_cast at h2; exact h2
  have hkfac : (Nat.factorial k : K.carrier) ≠ 0 := Nat.cast_ne_zero.mpr (Nat.factorial_ne_zero k)
  have hlfac : (Nat.factorial (n - k) : K.carrier) ≠ 0 := Nat.cast_ne_zero.mpr (Nat.factorial_ne_zero (n - k))
  rw [div_mul_div_comm, div_mul_eq_mul_div, div_eq_iff (mul_ne_zero hkfac hlfac), ← hfactK]
  ring

/-- **p 進指数関数**。 -/
noncomputable def padicExp (K : PAdicLocalField p) (x : K.carrier) : K.carrier :=
  ∑' n, expTerm K x n

/-- **`padicExp` は加法的**(`sorry` 無し): `exp(x+y) = exp(x)·exp(y)`。

★`x+y` 自身の収束は要らない——`antidiagonal_sum_expTerm` が各 `n` ごとの
**恒等式**(収束の有無によらず成り立つ)なので、`∑'n, expTerm K (x+y) n` への
書き換えは `tsum_congr` だけで済む。`x` と `y` それぞれの収束(`hx`・`hy`)のみで
`padicExp K x * padicExp K y` の値を経由して結論が出る——超距離不等式により
`x+y` も自動的に同じ収束域に入る(`IsUltrametricDist.norm_add_le_max`)ので、
実害はない。 -/
theorem padicExp_add (K : PAdicLocalField p) {x y : K.carrier}
    (hx : ‖x‖ ^ (p - 1) * (p : ℝ) < 1) (hy : ‖y‖ ^ (p - 1) * (p : ℝ) < 1) :
    padicExp K (x + y) = padicExp K x * padicExp K y := by
  have hSx := expTerm_summable K hx
  have hSy := expTerm_summable K hy
  have hSprod : Summable (fun i : ℕ × ℕ => expTerm K x i.1 * expTerm K y i.2) :=
    hSx.mul_of_nonarchimedean hSy
  have key := hSx.tsum_mul_tsum_eq_tsum_sum_antidiagonal hSy hSprod
  rw [show (∑' n, ∑ kl ∈ Finset.antidiagonal n, expTerm K x kl.1 * expTerm K y kl.2)
        = ∑' n, expTerm K (x + y) n from tsum_congr (antidiagonal_sum_expTerm K x y)] at key
  unfold padicExp
  rw [key]

end ABC3.Found.PGC
