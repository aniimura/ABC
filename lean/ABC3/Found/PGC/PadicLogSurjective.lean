import ABC3.Found.PGC.PadicLogInjective
import Mathlib.Topology.MetricSpace.Contracting

/-!
# p進対数は半径 `1/4` の球からそれ自身への全単射(`sorry` 無し)

`Found/PGC/PadicLogInjective.lean` で確立した単射性に続き、**全射性**——
「`padicLog K` の像が `K` の開部分群に一致する」——を一般の局所体 `K` で確立する。
Proposition 2.1 の folklore 入力「p進対数が `U_K`(捩れを法として)を `K` の開
部分群に写す」は、これで**完全に**(単射・全射とも)片付いたことになる。

## 戦略(`exp∘log=id` の互逆性を経由しない)

`memory/padic-log-additivity-blocked.md` に記録したとおり、形式冪級数の
`exp∘log=id` はスターリング数的な組み合わせ論を要し重い。ここでは代わりに
**p進版の逆写像定理**(Newton 法/縮小写像)を直接使う——`padicLog K x = x + t(x)`
(`t(x) := ∑_{n≥2} logTerm K x n`、尾部)と書けるとき、`t` 自身が
**縮小写像**(Lipschitz 定数 `1/2 < 1`、`Found/PGC/PadicLogInjective.lean` の
`logTerm_tail_norm_le` の「差」版)であることを直接示し、mathlib の
`ContractingWith.efixedPoint'`(Banach の不動点定理、`Topology/MetricSpace/
Contracting.lean`)を `f(x) := w - t(x)` に適用する——不動点 `x` はちょうど
`padicLog K x = w` を満たす。

`x^n - y^n = (x-y)·∑_{i<n} xⁱyⁿ⁻¹⁻ⁱ`(mathlib `geom_sum₂_mul`)という古典的
因数分解が、超距離の下で「差」の評価(縮小写像であること)の鍵になる——
このおかげで `exp`/`log` の係数の組み合わせ論を一切経由しない。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC

variable {p : ℕ} [Fact p.Prime]

/-! ### 部品1: `xⁿ - yⁿ` の因数分解による評価(汎用、超距離のみを使う) -/

/-- `‖x‖,‖y‖ ≤ r` なら `‖xⁿ - yⁿ‖ ≤ ‖x - y‖ · r^(n-1)`。
`geom_sum₂_mul`(`x^n-y^n = (x-y)·∑ᵢ xⁱyⁿ⁻¹⁻ⁱ`)を超距離不等式で評価する。 -/
theorem norm_pow_sub_pow_le {K : PAdicLocalField p} {x y : K.carrier} {r : ℝ}
    (hr0 : 0 ≤ r) (hx : ‖x‖ ≤ r) (hy : ‖y‖ ≤ r) (n : ℕ) :
    ‖x ^ n - y ^ n‖ ≤ ‖x - y‖ * r ^ (n - 1) := by
  rcases n with _ | m
  · simp
  · have hfact : (∑ i ∈ Finset.range (m + 1), x ^ i * y ^ (m + 1 - 1 - i)) * (x - y)
        = x ^ (m + 1) - y ^ (m + 1) := geom_sum₂_mul x y (m + 1)
    have hsum : ‖∑ i ∈ Finset.range (m + 1), x ^ i * y ^ (m + 1 - 1 - i)‖ ≤ r ^ m := by
      refine IsUltrametricDist.norm_sum_le_of_forall_le_of_nonneg (by positivity) ?_
      intro i hi
      rw [Finset.mem_range] at hi
      have hile : i ≤ m := by omega
      rw [norm_mul, norm_pow, norm_pow]
      calc ‖x‖ ^ i * ‖y‖ ^ (m + 1 - 1 - i) ≤ r ^ i * r ^ (m + 1 - 1 - i) := by gcongr
        _ = r ^ (i + (m + 1 - 1 - i)) := by rw [pow_add]
        _ = r ^ m := by congr 1; omega
    calc ‖x ^ (m + 1) - y ^ (m + 1)‖
        = ‖(∑ i ∈ Finset.range (m + 1), x ^ i * y ^ (m + 1 - 1 - i)) * (x - y)‖ := by rw [hfact]
      _ = ‖∑ i ∈ Finset.range (m + 1), x ^ i * y ^ (m + 1 - 1 - i)‖ * ‖x - y‖ := norm_mul _ _
      _ ≤ r ^ m * ‖x - y‖ := mul_le_mul_of_nonneg_right hsum (norm_nonneg _)
      _ = ‖x - y‖ * r ^ m := mul_comm _ _
      _ = ‖x - y‖ * r ^ (m + 1 - 1) := by norm_num

/-- 対数級数の項の「差」の評価: `‖logTerm K x n - logTerm K y n‖ ≤ n·‖x-y‖·r^(n-1)`。 -/
theorem logTerm_diff_norm_le {K : PAdicLocalField p} {x y : K.carrier} {r : ℝ}
    (hr0 : 0 ≤ r) (hx : ‖x‖ ≤ r) (hy : ‖y‖ ≤ r) {n : ℕ} (hn : n ≠ 0) :
    ‖logTerm K x n - logTerm K y n‖ ≤ (n : ℝ) * ‖x - y‖ * r ^ (n - 1) := by
  have heq : logTerm K x n - logTerm K y n
      = (-1 : K.carrier) ^ (n + 1) * (x ^ n - y ^ n) / (n : K.carrier) := by
    unfold logTerm
    rw [if_neg hn, if_neg hn]
    field_simp
  rw [heq, norm_div, norm_mul, norm_pow, norm_neg, norm_one, one_pow, one_mul]
  have h1 : (1 : ℝ) / n ≤ ‖((n : ℕ) : K.carrier)‖ := one_div_le_norm_natCast K hn
  have hpos : 0 < ‖((n : ℕ) : K.carrier)‖ := lt_of_lt_of_le (by positivity) h1
  rw [div_le_iff₀ hpos]
  calc ‖x ^ n - y ^ n‖ ≤ ‖x - y‖ * r ^ (n - 1) := norm_pow_sub_pow_le hr0 hx hy n
    _ = ‖x - y‖ * r ^ (n - 1) * ((n : ℝ) * (1 / n)) := by
          have hn' : (n : ℝ) ≠ 0 := Nat.cast_ne_zero.mpr hn
          field_simp
    _ ≤ ‖x - y‖ * r ^ (n - 1) * ((n : ℝ) * ‖((n : ℕ) : K.carrier)‖) := by
          have hpos2 : (0 : ℝ) ≤ ‖x - y‖ * r ^ (n - 1) := by positivity
          have hn0 : (0 : ℝ) ≤ (n : ℝ) := Nat.cast_nonneg n
          gcongr
    _ = (n : ℝ) * ‖x - y‖ * r ^ (n - 1) * ‖((n : ℕ) : K.carrier)‖ := by ring

/-- 実数の不等式 `(m+2)·(1/4)^(m+1) ≤ 1/2`(`m + 2 ≤ 2·4^m` の書き換え)。
超距離の sup 評価だけで尾部のリプシッツ定数 `1/2` が**そのまま**出る鍵。 -/
private theorem key_bound (m : ℕ) : ((m : ℝ) + 2) * (1 / 4 : ℝ) ^ (m + 1) ≤ 1 / 2 := by
  have h : (m : ℝ) + 2 ≤ 2 * 4 ^ m := by
    induction m with
    | zero => norm_num
    | succ k ih =>
      have h4 : (0 : ℝ) ≤ 4 ^ k := by positivity
      have hpow : (4 : ℝ) ^ (k + 1) = 4 * 4 ^ k := by ring
      push_cast [hpow]
      linarith [ih, h4, Nat.cast_nonneg (α := ℝ) k]
  have hpow : (1 / 4 : ℝ) ^ (m + 1) = (1 / 4) ^ m * (1 / 4) := by ring
  rw [hpow]
  have heq : (1 / 4 : ℝ) ^ m = ((4 : ℝ) ^ m)⁻¹ := by rw [one_div, inv_pow]
  calc ((m : ℝ) + 2) * ((1 / 4 : ℝ) ^ m * (1 / 4)) = ((m : ℝ) + 2) * (1 / 4) * (1 / 4) ^ m := by
        ring
    _ ≤ (2 * 4 ^ m) * (1 / 4) * (1 / 4) ^ m := by
        have : (0 : ℝ) ≤ (1 / 4 : ℝ) * (1 / 4) ^ m := by positivity
        nlinarith [h, this]
    _ = 1 / 2 := by rw [heq]; field_simp; norm_num

/-! ### 部品2: 尾部は縮小写像(Lipschitz 定数 `1/2`) -/

/-- ★**`‖x‖,‖y‖ ≤ 1/4` の範囲で、対数級数の尾部(次数 `≥2`)は Lipschitz 定数
`1/2` を持つ**。`x^n-y^n` の因数分解評価 + 超距離の sup 評価(実際の和ではない)
のおかげで、リプシッツ定数がぴったり `1/2`(`n=2` の項で達成、`n≥3` はさらに小さい)
と分かる——`aux_le_pow'`(`Found/PGC/PadicLogInjective.lean`)と同型の不等式が
ここでも(添字をずらして)効く。 -/
theorem tail_diff_norm_le {K : PAdicLocalField p} {x y : K.carrier}
    (hx : ‖x‖ ≤ 1 / 4) (hy : ‖y‖ ≤ 1 / 4) :
    ‖(∑' n : ℕ, logTerm K x (n + 2)) - ∑' n : ℕ, logTerm K y (n + 2)‖ ≤ (1 / 2) * ‖x - y‖ := by
  have hx' : ‖x‖ < 1 := lt_of_le_of_lt hx (by norm_num)
  have hy' : ‖y‖ < 1 := lt_of_le_of_lt hy (by norm_num)
  have hsx : Summable (fun n : ℕ => logTerm K x (n + 2)) :=
    (summable_nat_add_iff 2).mpr (logTerm_summable K hx')
  have hsy : Summable (fun n : ℕ => logTerm K y (n + 2)) :=
    (summable_nat_add_iff 2).mpr (logTerm_summable K hy')
  rw [← hsx.tsum_sub hsy]
  refine IsUltrametricDist.norm_tsum_le_of_forall_le_of_nonneg (by positivity) ?_
  intro n
  have hbound := logTerm_diff_norm_le (r := (1 : ℝ) / 4) (by norm_num) hx hy
    (n := n + 2) (Nat.succ_ne_zero (n + 1))
  have hsimp : n + 2 - 1 = n + 1 := by omega
  rw [hsimp] at hbound
  calc ‖logTerm K x (n + 2) - logTerm K y (n + 2)‖
      ≤ ((n : ℝ) + 2) * ‖x - y‖ * (1 / 4) ^ (n + 1) := by push_cast at hbound ⊢; linarith [hbound]
    _ = (((n : ℝ) + 2) * (1 / 4) ^ (n + 1)) * ‖x - y‖ := by ring
    _ ≤ (1 / 2) * ‖x - y‖ := mul_le_mul_of_nonneg_right (key_bound n) (norm_nonneg _)

/-! ### 部品3: ★★Banach の不動点定理を当てる —— 全射性 -/

/-- ★★**`padicLog K` は半径 `1/4` の球からそれ自身への全射**。
`f(x) := w - t(x)`(`t` = 対数級数の尾部)の不動点を`ContractingWith.efixedPoint'`
(mathlib、`Topology/MetricSpace/Contracting.lean`)で取る——`f` は球を球に写し
(`t` の像が `‖x‖≤1/4` のとき `≤1/4` に収まる)、`1/2`-縮小写像(`tail_diff_norm_le`)
なので不動点が存在し、`x = f(x) = w - t(x)` すなわち `padicLog K x = x+t(x) = w`。 -/
theorem padicLog_surjOn (K : PAdicLocalField p) :
    Set.SurjOn (padicLog K) {x : K.carrier | ‖x‖ ≤ 1 / 4} {w : K.carrier | ‖w‖ ≤ 1 / 4} := by
  intro w hw
  simp only [Set.mem_setOf_eq] at hw
  set s : Set K.carrier := {x : K.carrier | ‖x‖ ≤ 1 / 4} with hs_def
  set f : K.carrier → K.carrier := fun x => w - ∑' n : ℕ, logTerm K x (n + 2) with hf_def
  have hsclosed : s = Metric.closedBall (0 : K.carrier) (1 / 4) := by
    ext x; simp [hs_def, Metric.mem_closedBall, dist_eq_norm]
  have hsc : IsComplete s := by rw [hsclosed]; exact Metric.isClosed_closedBall.isComplete
  have hMaps : Set.MapsTo f s s := by
    intro x hx
    simp only [hs_def, Set.mem_setOf_eq] at hx ⊢
    have htail : ‖∑' n : ℕ, logTerm K x (n + 2)‖ ≤ 2 * ‖x‖ ^ 2 :=
      IsUltrametricDist.norm_tsum_le_of_forall_le_of_nonneg (by positivity)
        (fun n => logTerm_tail_norm_le K hx n)
    have hxsq : ‖x‖ * ‖x‖ ≤ (1 / 4) * (1 / 4) := mul_le_mul hx hx (norm_nonneg x) (by norm_num)
    have htail' : ‖∑' n : ℕ, logTerm K x (n + 2)‖ ≤ 1 / 4 := by nlinarith [sq_nonneg ‖x‖]
    calc ‖f x‖ = ‖w - ∑' n : ℕ, logTerm K x (n + 2)‖ := rfl
      _ ≤ max ‖w‖ ‖∑' n : ℕ, logTerm K x (n + 2)‖ := by
          rw [sub_eq_add_neg]; exact (IsUltrametricDist.norm_add_le_max w _).trans_eq (by rw [norm_neg])
      _ ≤ 1 / 4 := max_le hw htail'
  have hLip : LipschitzWith (1 / 2 : NNReal) (hMaps.restrict f s s) := by
    apply LipschitzWith.of_dist_le_mul
    rintro ⟨x, hx⟩ ⟨y, hy⟩
    simp only [hs_def, Set.mem_setOf_eq] at hx hy
    show dist (f x) (f y) ≤ ((1 / 2 : NNReal) : ℝ) * dist x y
    rw [dist_eq_norm, dist_eq_norm]
    have hb := tail_diff_norm_le hx hy
    have heqfxy : f x - f y = (∑' n : ℕ, logTerm K y (n + 2)) - ∑' n : ℕ, logTerm K x (n + 2) := by
      simp only [hf_def]; ring
    rw [heqfxy, norm_sub_rev]
    push_cast
    linarith [hb]
  have hcw : ContractingWith (1 / 2 : NNReal) (hMaps.restrict f s s) := ⟨by norm_num, hLip⟩
  have hne : s.Nonempty := ⟨0, by simp [hs_def]⟩
  obtain ⟨x0, hx0s⟩ := hne
  have hx0 : edist x0 (f x0) ≠ ⊤ := edist_ne_top _ _
  set x := ContractingWith.efixedPoint' f hsc hMaps hcw x0 hx0s hx0 with hxdef
  have hxs : x ∈ s := ContractingWith.efixedPoint_mem' hsc hMaps hcw hx0s hx0
  have hxfix : Function.IsFixedPt f x := ContractingWith.efixedPoint_isFixedPt' hsc hMaps hcw hx0s hx0
  refine ⟨x, hxs, ?_⟩
  have hxeq : x = w - ∑' n : ℕ, logTerm K x (n + 2) := hxfix.symm
  simp only [hs_def, Set.mem_setOf_eq] at hxs
  have hlogx : padicLog K x = x + ∑' n : ℕ, logTerm K x (n + 2) :=
    padicLog_eq_add_tail K (lt_of_le_of_lt hxs (by norm_num))
  rw [hlogx]
  linear_combination hxeq

/-- ★★★**`padicLog K` は半径 `1/4` の球からそれ自身への全単射**——
Proposition 2.1 の folklore 入力(「p進対数が `U_K` を `K` の開部分群に写す」)を
一般の局所体 `K` について**完全に**確立した(単射: `padicLog_injOn`、
全射: `padicLog_surjOn`)。`{y : K.carrier | ‖y‖ ≤ 1/4}` は超距離不等式により
`(K,+)` の(開かつコンパクトな)部分群である。 -/
theorem padicLog_bijOn (K : PAdicLocalField p) :
    Set.BijOn (padicLog K) {x : K.carrier | ‖x‖ ≤ 1 / 4} {w : K.carrier | ‖w‖ ≤ 1 / 4} := by
  refine ⟨?_, padicLog_injOn K, padicLog_surjOn K⟩
  intro x hx
  simp only [Set.mem_setOf_eq] at hx ⊢
  have htail : ‖∑' n : ℕ, logTerm K x (n + 2)‖ ≤ 2 * ‖x‖ ^ 2 :=
    IsUltrametricDist.norm_tsum_le_of_forall_le_of_nonneg (by positivity)
      (fun n => logTerm_tail_norm_le K hx n)
  have hxsq : ‖x‖ * ‖x‖ ≤ (1 / 4) * (1 / 4) := mul_le_mul hx hx (norm_nonneg x) (by norm_num)
  have htail' : ‖∑' n : ℕ, logTerm K x (n + 2)‖ ≤ 1 / 4 := by nlinarith [sq_nonneg ‖x‖]
  have hlogx : padicLog K x = x + ∑' n : ℕ, logTerm K x (n + 2) :=
    padicLog_eq_add_tail K (lt_of_le_of_lt hx (by norm_num))
  rw [hlogx]
  exact (IsUltrametricDist.norm_add_le_max x _).trans (max_le hx htail')

end ABC3.Found.PGC
