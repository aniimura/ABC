import ABC3.Found.PGC.PadicLogMul

/-!
# p進対数は `‖x‖ ≤ 1/4` の球上で単射(`sorry` 無し)

`Found/PGC/PadicLogMul.lean` で確立した `padicLog_mul`(準同型性)に続き、
Proposition 2.1 の `.needs` の folklore 入力「p進対数が `U_K`(捩れを法として)を
`K` の開部分群に写す」の**単射性**側を、一般の局所体 `K` について確立する。

## 戦略

`ABC3.Found.IUTchIII.norm_logOneAdd_eq`(ℚ_p 限定)と同じ計算——「`‖x‖ ≤ 1/4` なら
`padicLog K x` の先頭項(次数1、= `x` 自身)が残りの項(次数 `≥2`、`≤ 2‖x‖²`)を
超距離的に支配するので `‖padicLog K x‖ = ‖x‖`」——を一般の `K` に一般化する
(`Found/PGC/PadicLog.lean` の `logTerm`/`logTerm_norm_le` を土台にすれば、
実数側の不等式(`aux_le_pow'` 等)はそのまま流用できる)。

これと `padicLog_mul` を組み合わせると、`padicLog K` の**核が自明**
(`‖x‖≤1/4` の範囲で `padicLog K x = 0 → x = 0`)になり、そこから
「乗法的な円周演算 `x * y := x + y + x*y`」に関する**単射性**が出る——
2元 `x, y`(`‖x‖,‖y‖≤1/4`)から `z := (x-y)/(1+y)` を作ると
`x = z + y + z*y` かつ `‖z‖ = ‖x-y‖ ≤ 1/4` なので、
`padicLog K x = padicLog K y` なら `padicLog_mul` で `padicLog K z = 0`、
核の自明性で `z = 0`、よって `x = y`。

## まだ無いもの

- **全射性**(`padicLog` の像が `K` の開部分群に一致すること)。`padicExp` を
  逆写像として使えば得られる見込みだが、`exp∘log=id` の互逆性は未着手
  (`memory/padic-log-additivity-blocked.md` 参照)。
- 半径 `1/4` は十分条件であって、原文が要求する「`U_K` を(捩れを法として)」の
  正確な範囲(`1+𝔪_K` 全体、`𝔪_K` の分岐指数に応じて半径が変わりうる)との
  一致は未検証——ℚ_p の場合の `1/4` がそのまま一般の `K` でも十分なのは、
  この証明が `p` にも `K` の分岐指数にも依存しない純粋な超距離不等式だから
  (`aux_le_pow'` は実数の不等式であって `K` の情報を一切使わない)。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC

variable {p : ℕ} [Fact p.Prime]

/-- 実数の不等式 `n + 2 ≤ 2·4ⁿ`(`ABC3.Found.IUTchIII` の `private` 版と同じ内容——
`K` にも `p` にも依存しない純粋な帰納法)。 -/
private theorem aux_le_pow' (n : ℕ) : (n : ℝ) + 2 ≤ 2 * 4 ^ n := by
  induction n with
  | zero => norm_num
  | succ k ih =>
    have h4 : (0 : ℝ) ≤ 4 ^ k := by positivity
    have hpow : (4 : ℝ) ^ (k + 1) = 4 * 4 ^ k := by ring
    push_cast [hpow]
    linarith [ih, h4, Nat.cast_nonneg (α := ℝ) k]

/-- 次数 `≥2` の項の一様評価: `‖x‖ ≤ 1/4` なら `logTerm K x (n+2) ≤ 2‖x‖²`。 -/
theorem logTerm_tail_norm_le (K : PAdicLocalField p) {x : K.carrier}
    (hx : ‖x‖ ≤ 1 / 4) (n : ℕ) : ‖logTerm K x (n + 2)‖ ≤ 2 * ‖x‖ ^ 2 := by
  have h0 : (0 : ℝ) ≤ ‖x‖ := norm_nonneg x
  refine (logTerm_norm_le K (n + 2)).trans ?_
  have hstep : ‖x‖ ^ n ≤ (1 / 4 : ℝ) ^ n := by gcongr
  have hkey : ((n : ℝ) + 2) * ‖x‖ ^ n ≤ 2 := by
    calc ((n : ℝ) + 2) * ‖x‖ ^ n
        ≤ ((n : ℝ) + 2) * (1 / 4 : ℝ) ^ n := by gcongr
      _ = ((n : ℝ) + 2) / 4 ^ n := by rw [div_pow, one_pow]; ring
      _ ≤ 2 := by
          rw [div_le_iff₀ (by positivity)]
          linarith [aux_le_pow' n]
  calc ((n + 2 : ℕ) : ℝ) * ‖x‖ ^ (n + 2)
      = (((n : ℝ) + 2) * ‖x‖ ^ n) * ‖x‖ ^ 2 := by push_cast; ring
    _ ≤ 2 * ‖x‖ ^ 2 := by gcongr

/-- `padicLog K x` を先頭2項(次数0=0・次数1=x)と、次数 `≥2` の尾部に分ける。 -/
theorem padicLog_eq_add_tail (K : PAdicLocalField p) {x : K.carrier} (hx : ‖x‖ < 1) :
    padicLog K x = x + ∑' n : ℕ, logTerm K x (n + 2) := by
  have hs : Summable (logTerm K x) := logTerm_summable K hx
  have hs1 : Summable (fun n => logTerm K x (n + 1)) := (summable_nat_add_iff 1).mpr hs
  rw [padicLog, hs.tsum_eq_zero_add]
  have h0 : logTerm K x 0 = 0 := by simp [logTerm]
  rw [h0, zero_add, hs1.tsum_eq_zero_add]
  have h1 : logTerm K x (0 + 1) = x := by simp [logTerm]
  rw [h1]

/-- ★**`‖x‖ ≤ 1/4` なら `‖padicLog K x‖ = ‖x‖`**(一般の `K`)。
`ABC3.Found.IUTchIII.norm_logOneAdd_eq`(ℚ_p 限定)の一般化。 -/
theorem norm_padicLog_eq (K : PAdicLocalField p) {x : K.carrier}
    (hx : ‖x‖ ≤ 1 / 4) (hx0 : x ≠ 0) : ‖padicLog K x‖ = ‖x‖ := by
  have hlt1 : ‖x‖ < 1 := lt_of_le_of_lt hx (by norm_num)
  have hpos : 0 < ‖x‖ := norm_pos_iff.mpr hx0
  have htail : ‖∑' n : ℕ, logTerm K x (n + 2)‖ ≤ 2 * ‖x‖ ^ 2 :=
    IsUltrametricDist.norm_tsum_le_of_forall_le_of_nonneg (by positivity)
      (fun n => logTerm_tail_norm_le K hx n)
  have hsmall : 2 * ‖x‖ ^ 2 < ‖x‖ := by nlinarith
  have hne : ‖x‖ ≠ ‖∑' n : ℕ, logTerm K x (n + 2)‖ := (ne_of_lt (lt_of_le_of_lt htail hsmall)).symm
  rw [padicLog_eq_add_tail K hlt1, IsUltrametricDist.norm_add_eq_max_of_norm_ne_norm hne]
  exact max_eq_left (le_of_lt (lt_of_le_of_lt htail hsmall))

/-- ★★**`padicLog K` は `‖x‖ ≤ 1/4` の球上で単射**(一般の `K`)。

Proposition 2.1 の `.needs` の folklore 入力(「p進対数が `U_K` を `K` の開部分群に
写す」)の単射性側。`padicLog_mul`(乗法的な円周演算との両立性)+
`norm_padicLog_eq`(核が自明)から、2元の差を円周演算で「割った」商 `z` を作れば
出る。全射性(像が実際に開部分群になること)は依然として未着手
(module docstring 参照)。 -/
theorem padicLog_injOn (K : PAdicLocalField p) :
    Set.InjOn (padicLog K) {x : K.carrier | ‖x‖ ≤ 1 / 4} := by
  intro x hx y hy hxy
  simp only [Set.mem_setOf_eq] at hx hy
  by_contra hne
  have hnormy1 : ‖(1 : K.carrier) + y‖ = 1 := by
    rw [IsUltrametricDist.norm_add_eq_max_of_norm_ne_norm (S := K.carrier)]
    · rw [norm_one]; exact max_eq_left (le_trans hy (by norm_num))
    · rw [norm_one]; exact ne_of_gt (lt_of_le_of_lt hy (by norm_num))
  have hy1 : (1 : K.carrier) + y ≠ 0 := by
    intro h
    rw [h, norm_zero] at hnormy1
    exact one_ne_zero hnormy1.symm
  obtain ⟨z, hclear⟩ : ∃ z : K.carrier, z * (1 + y) = x - y :=
    ⟨(x - y) / (1 + y), div_mul_cancel₀ (x - y) hy1⟩
  have hxyz : x = z + y + z * y := by linear_combination -hclear
  have hzy : ‖z‖ = ‖x - y‖ := by
    have heq : z = (x - y) / (1 + y) := (eq_div_iff hy1).mpr hclear
    rw [heq, norm_div, hnormy1, div_one]
  have hxysub : ‖x - y‖ ≤ 1 / 4 := by
    calc ‖x - y‖ = ‖x + (-y)‖ := by rw [sub_eq_add_neg]
      _ ≤ max ‖x‖ ‖(-y)‖ := IsUltrametricDist.norm_add_le_max x (-y)
      _ = max ‖x‖ ‖y‖ := by rw [norm_neg]
      _ ≤ 1 / 4 := max_le hx hy
  have hzle : ‖z‖ ≤ 1 / 4 := hzy ▸ hxysub
  have hz0 : z ≠ 0 := by
    intro h0
    apply hne
    have hxy0 : x - y = 0 := by rw [← hclear, h0, zero_mul]
    exact sub_eq_zero.mp hxy0
  have hlog : padicLog K x = padicLog K z + padicLog K y := by
    rw [hxyz]
    exact padicLog_mul (lt_of_le_of_lt hzle (by norm_num)) (lt_of_le_of_lt hy (by norm_num))
  rw [hxy] at hlog
  have hzlog : padicLog K z = 0 := by
    have h2 : padicLog K z + padicLog K y = 0 + padicLog K y := by rw [zero_add]; exact hlog.symm
    exact add_right_cancel h2
  have hfinal : ‖padicLog K z‖ = ‖z‖ := norm_padicLog_eq K hzle hz0
  rw [hzlog, norm_zero] at hfinal
  exact hz0 (norm_eq_zero.mp hfinal.symm)

end ABC3.Found.PGC
