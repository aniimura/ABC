import ABC3.Found.PGC.PadicLog
import ABC3.Found.PGC.PadicExp
import ABC3.Found.IUTchIII.PadicLog
import ABC3.Found.IUTchIII.PowerSeriesLog
import Mathlib.RingTheory.PowerSeries.Substitution

/-!
# 一般の局所体 `K` 上の p進対数の乗法性(`sorry` 無し)

`Found/PGC/PadicLog.lean` は `K : PAdicLocalField p`(ℚ_p の任意の有限次拡大)上の
p進対数 `padicLog K : K.carrier → K.carrier` の**収束**のみを確立し、
「まだ無いもの」として**準同型性** `log((1+x)(1+y)) = log(1+x)+log(1+y)` を
`.needs`(`Skeleton/PGC/Section1.lean` の Proposition 1.2)に残していた。

★★2026-09-04: この準同型性は**すでに ℚ_p の場合だけ本リポジトリの別セクション
(`Found/IUTchIII/PadicLogMul.lean`)に存在した**——`logOneAdd_mul`(在庫確認、
`node tools/decl-index.mjs` の再生成後の grep で発見)。IUTchIII 側は `AbsTopIII`
の log-shell(ℚ_p 上)のために作られたもので、pGC 側の一般の `K` への一般化は
未着手のままだった。**本ファイルはその一般化を行う**——証明の筋は
`Found/IUTchIII/PadicLogMul.lean` と完全に並行(定理名・補題の対応は下記コメント参照)、
`ℚ_[p]` を `K.carrier` に、`ABC3.Found.IUTchIII.logOneAdd` を `ABC3.Found.PGC.padicLog`
に置き換えただけで、追加で要ったのは以下の2つの scoped instance のみ:

- `CharZero K.carrier`(`padicCharZero`、既存)
- そこから `DivisionRing.toRatAlgebra` で自動的に得られる `Algebra ℚ K.carrier`
- `IsAddTorsionFree K.carrier`(`IsAddTorsionFree.of_module_rat`)

これらを scoped instance として登録すると、mathlib 側の一般の道具
(`PowerSeries.log`/`PowerSeries.logOf`/`ABC3.Found.IUTchIII.logOf_mul`)が
そのまま `K.carrier` に対して動く。

## 結論(Proposition 1.2 の `.needs` への影響)

Proposition 1.2 の3つの `.needs` のうち「p進対数」の項目は**この定理で解消した**
(mathlib 不在 → 本ファイルで sorry 無しに構築)。残り2項目
(相互律の同型 `Γ_K^ab ≅ (K^×)^∧` そのもの・完全系列の典拠)は**依然として不在**——
本ファイルは Proposition 1.2 全体を閉じない。局所類体論の相互写像そのものが
absent なままである限り、`residueCard_and_degree_recoverable` は `sorry` を残す。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC PowerSeries

variable {p : ℕ} [Fact p.Prime]

/-! ### 部品0: `K.carrier` を `ℚ` 上の(捩れ無し)代数にする scoped instance -/

/-- `K.carrier` は特性 0(`padicCharZero`、`Found/PGC/PadicExp.lean` 既存)。
これがあると `DivisionRing.toRatAlgebra`(mathlib、`Algebra/Algebra/Rat.lean`)が
自動的に `Algebra ℚ K.carrier` を与える——`PowerSeries.log`/`logOf` が要求するもの。 -/
noncomputable scoped instance instCharZeroCarrier (K : PAdicLocalField p) :
    CharZero K.carrier := padicCharZero K

/-- `K.carrier` は(ℚ上の加群として)捩れ無し——`IsAddTorsionFree.of_module_rat` から。
`ABC3.Found.IUTchIII.logOf_mul` が要求する仮説。 -/
noncomputable scoped instance instIsAddTorsionFreeCarrier (K : PAdicLocalField p) :
    IsAddTorsionFree K.carrier :=
  IsAddTorsionFree.of_module_rat (M := K.carrier)

/-! ### 部品1: 形式的対数の係数のノルム(`ℚ_[p]` の場合から埋め込みで移す) -/

/-- `‖(n:K.carrier)‖ ≥ 1/n`。`K.carrier` へのノルムは `algebraMap ℚ_[p] K.carrier` で
等長(`norm_algebraMap`、`Found/PGC/LocalFieldNorm.lean`)なので、
`ABC3.Found.IUTchIII.one_div_le_norm_natCast`(ℚ_p の場合)からそのまま移せる。 -/
theorem one_div_le_norm_natCast (K : PAdicLocalField p) {n : ℕ} (hn : n ≠ 0) :
    (1 : ℝ) / n ≤ ‖((n : ℕ) : K.carrier)‖ := by
  have hnormK : ‖(n : K.carrier)‖ = ‖(n : ℚ_[p])‖ := by
    rw [show ((n : ℕ) : K.carrier) = algebraMap ℚ_[p] K.carrier ((n : ℕ) : ℚ_[p]) from
          (map_natCast _ _).symm, ABC3.Found.PGC.norm_algebraMap]
  rw [hnormK]
  exact ABC3.Found.IUTchIII.one_div_le_norm_natCast hn

/-- `‖coeff d (log K.carrier)‖ ≤ d`(`ABC3.Found.IUTchIII.norm_coeff_log_le` の一般化)。 -/
theorem norm_coeff_log_le (K : PAdicLocalField p) :
    ∀ d : ℕ, ‖coeff d (log K.carrier)‖ ≤ d := by
  intro d
  rcases Nat.eq_zero_or_pos d with rfl | hd
  · simp
  · rw [coeff_log, if_neg hd.ne']
    have hcast : (algebraMap ℚ K.carrier) ((-1 : ℚ) ^ (d + 1) / d)
        = (-1 : K.carrier) ^ (d + 1) / (d : K.carrier) := by
      rw [map_div₀, map_pow, map_neg, map_one, map_natCast]
    rw [hcast, norm_div, norm_pow, norm_neg, norm_one, one_pow]
    have h1 : (1 : ℝ) / d ≤ ‖((d : ℕ) : K.carrier)‖ := one_div_le_norm_natCast K hd.ne'
    have hdpos : (0 : ℝ) < d := by exact_mod_cast hd
    have hnpos : 0 < ‖((d : ℕ) : K.carrier)‖ := lt_of_lt_of_le (by positivity) h1
    rw [div_le_iff₀ hnpos]
    calc (1 : ℝ) = d * (1 / d) := by field_simp
      _ ≤ d * ‖((d : ℕ) : K.carrier)‖ := mul_le_mul_of_nonneg_left h1 hdpos.le

/-! ### 部品2: 多項式の冪の係数のノルム(`K.carrier` の超距離性のみを使う、汎用) -/

/-- 係数がすべてノルム `≤ r` なら、`P^d` の係数はノルム `≤ r^d`。 -/
theorem norm_coeff_pow_le {K : PAdicLocalField p} {P : Polynomial K.carrier} {r : ℝ}
    (hr0 : 0 ≤ r) (hr : ∀ n, ‖P.coeff n‖ ≤ r) (d : ℕ) :
    ∀ n, ‖(P ^ d).coeff n‖ ≤ r ^ d := by
  induction d with
  | zero =>
    intro n
    simp only [pow_zero, Polynomial.coeff_one]
    split_ifs <;> simp
  | succ k ih =>
    intro n
    rw [pow_succ, Polynomial.coeff_mul]
    refine IsUltrametricDist.norm_sum_le_of_forall_le_of_nonneg (by positivity) ?_
    rintro ⟨i, j⟩ -
    rw [norm_mul, pow_succ]
    exact mul_le_mul (ih i) (hr j) (norm_nonneg _) (by positivity)

/-! ### 部品3: 係数の総和 = `eval 1`(ノルム不要、汎用) -/

theorem tsum_coeff_eq_eval_one {K : PAdicLocalField p} (Q : Polynomial K.carrier) :
    ∑' n, Q.coeff n = Q.eval 1 := by
  rw [tsum_eq_sum (s := Finset.range (Q.natDegree + 1)) ?_, Polynomial.eval_eq_sum_range]
  · simp only [one_pow, mul_one]
  · intro n hn
    have h : Q.natDegree + 1 ≤ n := by simpa using hn
    exact Q.coeff_eq_zero_of_natDegree_lt h

/-! ### 部品4: `padicLog` を形式的対数の係数で書き直す

`Found/PGC/PadicLog.lean` の `logTerm` は添字 `n` が**そのまま次数**になる規約
(`n=0` はダミーの `0`)なので、`ABC3.Found.IUTchIII` 版のようなずらしが要らない
——`logTerm K x n` と `coeff n (log K.carrier) * x^n` は全ての `n` で一致する。 -/

theorem logTerm_eq_coeff_log_mul (K : PAdicLocalField p) (x : K.carrier) (n : ℕ) :
    logTerm K x n = coeff n (log K.carrier) * x ^ n := by
  rcases eq_or_ne n 0 with rfl | hn
  · simp [logTerm]
  · rw [logTerm, if_neg hn, coeff_log, if_neg hn, map_div₀, map_pow, map_neg, map_one,
      map_natCast]
    ring

/-- `padicLog K x = ∑' d, coeff d (log K.carrier) * x^d`。単なる添字の言い換え
(`‖x‖<1` の仮定すら不要——両辺とも同じ族の `tsum`)。 -/
theorem padicLog_eq_tsum_coeff_log (K : PAdicLocalField p) (x : K.carrier) :
    padicLog K x = ∑' d : ℕ, coeff d (log K.carrier) * x ^ d := by
  rw [padicLog]
  exact tsum_congr (logTerm_eq_coeff_log_mul K x)

/-! ### 部品5: `P^d` の係数の台(ノルム不要、汎用) -/

theorem coeff_pow_eq_zero_of_lt {K : PAdicLocalField p} {P : Polynomial K.carrier}
    (hP0 : P.coeff 0 = 0) : ∀ (d n : ℕ), n < d → (P ^ d).coeff n = 0 := by
  intro d
  induction d with
  | zero => intro n hn; exact absurd hn (Nat.not_lt_zero n)
  | succ k ih =>
    intro n hn
    rw [pow_succ, Polynomial.coeff_mul]
    refine Finset.sum_eq_zero ?_
    rintro ⟨i, j⟩ hij
    rw [Finset.mem_antidiagonal] at hij
    rcases Nat.eq_zero_or_pos j with rfl | hj
    · rw [hP0, mul_zero]
    · rw [ih i (by omega), zero_mul]

/-! ### 部品6: ★橋 —— 形式的対数の係数の総和 = `padicLog`

`ABC3.Found.IUTchIII.summable_and_tsum_coeff_logOf` の一般化。証明は一字一句
並行——`ℚ_[p]` を `K.carrier` に、`logOneAdd` を `padicLog K` に、
`NonarchimedeanAddGroup.summable_of_tendsto_cofinite_zero` を本プロジェクト独自の
`ABC3.Found.PGC.summable_of_tendsto_cofinite_zero`(`Found/PGC/PadicLog.lean`、
乗法版しか無かった mathlib の隙間を埋めた自前の一般形)に置き換えただけ。 -/
theorem summable_and_tsum_coeff_logOf {K : PAdicLocalField p} (P : Polynomial K.carrier)
    (hP0 : P.coeff 0 = 0) {r : ℝ} (hr0 : 0 ≤ r) (hr1 : r < 1) (hr : ∀ n, ‖P.coeff n‖ ≤ r) :
    (Summable fun n => PowerSeries.coeff n (logOf (1 + (P : K.carrier⟦X⟧)))) ∧
      ∑' n, PowerSeries.coeff n (logOf (1 + (P : K.carrier⟦X⟧))) = padicLog K (P.eval 1) := by
  classical
  set G : ℕ → ℕ → K.carrier := fun d n => coeff d (log K.carrier) * (P ^ d).coeff n with hG
  have hGzero_d : ∀ d n, n < d → G d n = 0 := fun d n h => by
    simp [hG, coeff_pow_eq_zero_of_lt hP0 d n h]
  have hGzero_n : ∀ d n, d * P.natDegree < n → G d n = 0 := fun d n h => by
    have h2 : (P ^ d).coeff n = 0 :=
      Polynomial.coeff_eq_zero_of_natDegree_lt (lt_of_le_of_lt Polynomial.natDegree_pow_le h)
    simp [hG, h2]
  have hGbound : ∀ d n, ‖G d n‖ ≤ d * r ^ d := by
    intro d n
    rw [hG, norm_mul]
    exact mul_le_mul (norm_coeff_log_le K d) (norm_coeff_pow_le hr0 hr d n)
      (norm_nonneg _) (Nat.cast_nonneg d)
  have h₁ : ∀ d, Summable (G d) := fun d =>
    (hasSum_sum_of_ne_finset_zero (s := Finset.range (d * P.natDegree + 1))
      (fun n hn => hGzero_n d n (by simpa using hn))).summable
  have h₂ : ∀ n, Summable fun d => G d n := fun n =>
    (hasSum_sum_of_ne_finset_zero (s := Finset.range (n + 1))
      (fun d hd => hGzero_d d n (by simpa using hd))).summable
  have hF : Summable (Function.uncurry G) := by
    refine summable_of_tendsto_cofinite_zero ?_
    rw [NormedAddGroup.tendsto_nhds_zero]
    intro ε hε
    obtain ⟨D, hD⟩ := (Metric.tendsto_atTop.mp
      (tendsto_self_mul_const_pow_of_lt_one hr0 hr1) ε hε)
    rw [Filter.eventually_cofinite]
    refine Set.Finite.subset
      ((Finset.range D ×ˢ Finset.range (D * P.natDegree + 1)).finite_toSet) ?_
    rintro ⟨d, n⟩ hq
    simp only [Set.mem_setOf_eq, not_lt] at hq
    have hne : G d n ≠ 0 := fun h => by simp [Function.uncurry, h] at hq; linarith
    have hdlt : d < D := by
      by_contra hcon
      have := hD d (not_lt.mp hcon)
      rw [Real.dist_eq, sub_zero] at this
      have hle : ‖G d n‖ ≤ d * r ^ d := hGbound d n
      have : (d : ℝ) * r ^ d < ε := lt_of_le_of_lt (le_abs_self _) this
      exact absurd hq (by simp only [Function.uncurry]; linarith)
    have hnle : n ≤ d * P.natDegree := by
      by_contra hcon
      exact hne (hGzero_n d n (not_le.mp hcon))
    simp only [Finset.coe_product, Set.mem_prod, Finset.mem_coe, Finset.mem_range]
    exact ⟨hdlt, by nlinarith [Nat.mul_le_mul_right P.natDegree (Nat.le_of_lt_succ
      (Nat.lt_succ_of_lt hdlt))]⟩
  have hsub : HasSubst ((P : K.carrier⟦X⟧)) :=
    HasSubst.of_constantCoeff_zero'
      (by rw [← PowerSeries.coeff_zero_eq_constantCoeff_apply, Polynomial.coeff_coe]; exact hP0)
  have hcoeff : ∀ n, PowerSeries.coeff n (logOf (1 + (P : K.carrier⟦X⟧))) = ∑' d, G d n := by
    intro n
    have hlog : logOf (1 + (P : K.carrier⟦X⟧)) = (log K.carrier).subst (P : K.carrier⟦X⟧) := by
      rw [logOf_eq, add_sub_cancel_left]
    rw [hlog, coeff_subst' hsub]
    have hterm : ∀ d,
        coeff d (log K.carrier) • PowerSeries.coeff n ((P : K.carrier⟦X⟧) ^ d) = G d n := by
      intro d
      rw [smul_eq_mul, hG, ← Polynomial.coe_pow, Polynomial.coeff_coe]
    simp only [hterm]
    rw [finsum_eq_sum_of_support_subset (s := Finset.range (n + 1)),
      tsum_eq_sum (s := Finset.range (n + 1)) (fun d hd => hGzero_d d n (by simpa using hd))]
    intro d hd
    simp only [Function.mem_support] at hd
    simp only [Finset.coe_range, Set.mem_Iio]
    by_contra hcon
    exact hd (hGzero_d d n (by omega))
  have hEval : ‖P.eval 1‖ < 1 := by
    rw [Polynomial.eval_eq_sum_range]
    refine lt_of_le_of_lt (IsUltrametricDist.norm_sum_le_of_forall_le_of_nonneg hr0 ?_) hr1
    intro i _
    simpa using hr i
  refine ⟨?_, ?_⟩
  · have hswap : Summable (Function.uncurry fun n d => G d n) := hF.prod_symm
    have := (hswap.hasSum.prod_fiberwise (fun n => (h₂ n).hasSum)).summable
    simpa only [← hcoeff] using this
  · calc ∑' n, PowerSeries.coeff n (logOf (1 + (P : K.carrier⟦X⟧)))
        = ∑' n, ∑' d, G d n := tsum_congr hcoeff
      _ = ∑' d, ∑' n, G d n := hF.tsum_comm' h₁ h₂
      _ = ∑' d, coeff d (log K.carrier) * (P.eval 1) ^ d := by
            refine tsum_congr fun d => ?_
            rw [hG]
            rw [tsum_mul_left, tsum_coeff_eq_eval_one, Polynomial.eval_pow]
      _ = padicLog K (P.eval 1) := (padicLog_eq_tsum_coeff_log K (P.eval 1)).symm

/-! ### 部品7: ★p進対数の乗法性(一般の局所体 `K`) -/

/-- `C a * X^k` の係数のノルムは `‖a‖` 以下(ノルム不要な汎用事実)。 -/
private theorem norm_coeff_CX_le {K : PAdicLocalField p} (a : K.carrier) (k n : ℕ) :
    ‖(Polynomial.C a * Polynomial.X ^ k).coeff n‖ ≤ ‖a‖ := by
  rw [Polynomial.coeff_C_mul, norm_mul, Polynomial.coeff_X_pow]
  split_ifs <;> simp

/-- ★★**一般の局所体 `K` 上の p進対数の乗法性**。`‖x‖, ‖y‖ < 1` のとき
`log((1+x)(1+y)) = log(1+x) + log(1+y)`。

`Skeleton/PGC/Section1.lean` の Proposition 1.2 の `.needs` に「absent」として
残っていた「p進対数」の入力を解消する——ただし相互律の同型そのもの
(`Γ_K^ab ≅ (K^×)^∧`)は依然として absent なので、Proposition 1.2 自体はこれだけでは
閉じない(module docstring 参照)。

証明は `ABC3.Found.IUTchIII.logOneAdd_mul`(ℚ_p の場合)と並行: 帳簿用の変数 `X` を
1つ入れた `K.carrier⟦X⟧` の中で `(1+Cx·X)(1+Cy·X) = 1+C(x+y)·X+C(xy)·X²` とし、
`ABC3.Found.IUTchIII.logOf_mul`(形式的加法公式、`[IsAddTorsionFree A]` のみを要求する
一般形)を当ててから、`summable_and_tsum_coeff_logOf`(上記、部品6)で係数の総和を
取る。`x`, `y` は**係数**であって代入する点ではないので、mathlib の評価 API
(`IsLinearTopology` を要求し体は満たさない)は一度も使わない。 -/
theorem padicLog_mul {K : PAdicLocalField p} {x y : K.carrier} (hx : ‖x‖ < 1) (hy : ‖y‖ < 1) :
    padicLog K (x + y + x * y) = padicLog K x + padicLog K y := by
  classical
  set r : ℝ := max ‖x‖ ‖y‖ with hrdef
  have hr0 : 0 ≤ r := le_trans (norm_nonneg x) (le_max_left _ _)
  have hr1 : r < 1 := max_lt hx hy
  set Px : Polynomial K.carrier := Polynomial.C x * Polynomial.X ^ 1 with hPx
  set Py : Polynomial K.carrier := Polynomial.C y * Polynomial.X ^ 1 with hPy
  set Pz : Polynomial K.carrier :=
    Polynomial.C (x + y) * Polynomial.X ^ 1 + Polynomial.C (x * y) * Polynomial.X ^ 2 with hPz
  have hPx0 : Px.coeff 0 = 0 := by simp [hPx]
  have hPy0 : Py.coeff 0 = 0 := by simp [hPy]
  have hPz0 : Pz.coeff 0 = 0 := by simp [hPz]
  have hbx : ∀ n, ‖Px.coeff n‖ ≤ r := fun n =>
    le_trans (by rw [hPx]; exact norm_coeff_CX_le x 1 n) (le_max_left _ _)
  have hby : ∀ n, ‖Py.coeff n‖ ≤ r := fun n =>
    le_trans (by rw [hPy]; exact norm_coeff_CX_le y 1 n) (le_max_right _ _)
  have hbz : ∀ n, ‖Pz.coeff n‖ ≤ r := by
    intro n
    rw [hPz, Polynomial.coeff_add]
    refine le_trans (IsUltrametricDist.norm_add_le_max _ _) (max_le ?_ ?_)
    · exact le_trans (norm_coeff_CX_le (x + y) 1 n)
        (le_trans (IsUltrametricDist.norm_add_le_max x y) (le_of_eq hrdef.symm))
    · refine le_trans (norm_coeff_CX_le (x * y) 2 n) ?_
      rw [norm_mul]
      calc ‖x‖ * ‖y‖ ≤ ‖x‖ * 1 := mul_le_mul_of_nonneg_left hy.le (norm_nonneg x)
        _ = ‖x‖ := mul_one _
        _ ≤ r := le_max_left _ _
  have hPoly : Px + Py + Px * Py = Pz := by
    simp only [hPx, hPy, hPz, Polynomial.C_add, Polynomial.C_mul]
    ring
  have hPS : (1 + (Px : K.carrier⟦X⟧)) * (1 + (Py : K.carrier⟦X⟧)) = 1 + (Pz : K.carrier⟦X⟧) := by
    rw [← hPoly]
    push_cast
    ring
  have hcx : PowerSeries.constantCoeff (1 + (Px : K.carrier⟦X⟧)) = 1 := by
    rw [map_add, map_one, ← PowerSeries.coeff_zero_eq_constantCoeff_apply,
      Polynomial.coeff_coe, hPx0, add_zero]
  have hcy : PowerSeries.constantCoeff (1 + (Py : K.carrier⟦X⟧)) = 1 := by
    rw [map_add, map_one, ← PowerSeries.coeff_zero_eq_constantCoeff_apply,
      Polynomial.coeff_coe, hPy0, add_zero]
  have hform : logOf (1 + (Pz : K.carrier⟦X⟧))
      = logOf (1 + (Px : K.carrier⟦X⟧)) + logOf (1 + (Py : K.carrier⟦X⟧)) := by
    rw [← hPS, ABC3.Found.IUTchIII.logOf_mul hcx hcy]
  obtain ⟨hsx, htx⟩ := summable_and_tsum_coeff_logOf Px hPx0 hr0 hr1 hbx
  obtain ⟨hsy, hty⟩ := summable_and_tsum_coeff_logOf Py hPy0 hr0 hr1 hby
  obtain ⟨_, htz⟩ := summable_and_tsum_coeff_logOf Pz hPz0 hr0 hr1 hbz
  have hze : Pz.eval 1 = x + y + x * y := by simp [hPz]
  have hxe : Px.eval 1 = x := by simp [hPx]
  have hye : Py.eval 1 = y := by simp [hPy]
  rw [← hze, ← htz, ← hxe, ← htx, ← hye, ← hty, hform]
  simp only [map_add]
  exact hsx.tsum_add hsy

end ABC3.Found.PGC
