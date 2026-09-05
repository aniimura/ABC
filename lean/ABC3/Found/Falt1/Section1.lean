import ABC3.Found.Falt1.Lemma11
import Mathlib.Analysis.SpecificLimits.Basic

/-!
# [Falt1] Chapter I §1 —— `Theorem 1.2` の最終段(`δₙ → 0`)(2026-09-05)

原典: G. Faltings, *p-Adic Hodge Theory*(1988)、Chapter I §1、
`Theorem 1.2`(物理 p.4 = 印字 p.257)。

> *1.2. Theorem. If `V ⊂ W` is any extension and `Wₙ` denotes the
> normalization of the tensor product `Vₙ ⊗_V W`, then the different
> `δ(Wₙ/Vₙ)` of `Wₙ` over `Vₙ` (or more precisely of each factor of the
> semilocal ring `Wₙ`) converges to 0 for `n → ∞`.*

## 原文の証明の骨格(2026-09-05 に読み直して判明)

`falt1-goal.md` はこれまで「`differentIdeal Wₙ Wₙ₊₁` の評価がブロッカー」と
記録していたが、Faltings の証明は **`differentIdeal` の tower diamond 公式
ではなく、`Lemma 1.1` を使った「長さ」の評価**である:

1. `Vₙ` を取り替えて `Wₙ` を整域としてよい。
2. 写像の列 `Ω_{Wₙ/Vₙ} ⊗ W_{n+1} → Ω_{W_{n+1}/Vₙ} → Ω_{W_{n+1}/V_{n+1}}` を見る。
3. 第 2 写像の核は `Ω_{V_{n+1}/Vₙ} ⊗ W_{n+1}` を含み、これは
   `(W_{n+1}/pW_{n+1})^{d+1}` を商に持つ。
4. `Ω_{V_{n+1}/Vₙ}` が `V_{n+1}/p^a V_{n+1}` 型の `d+1` 個の直和なので、
   合成は `Ω_{Wₙ/Vₙ} ⊗ W_{n+1}` 上の `p` 倍の核を零化する。
5. `δₙ` を差積とすると、その核の長さは
   `W_{n+1}/p^{δ'} W_{n+1}`(`δ' = min{1, δₙ/(d+1)}`)以上。
6. `p^{δₙ − δ_{n+1}}` が `W_{n+1}/(Wₙ ⊗_{Vₙ} V_{n+1})` を零化する
   (原文は *"it is clear that"* で畳んでいる)。
7. 以上から `δ_{n+1} ≤ δₙ − 1/(d+2)` または `δ_{n+1} ≤ (1 − 1/(d+2))·δₙ`
   が出て、*"In any case `δₙ → 0`"*。

`Lemma 1.1`(本トラックで証明済み——`Lemma11.lemma_1_1_falt1`)が
**差積 ↔ `Ω` の余核の長さ**の変換器になる:

    Module.length W Ω[W⁄V] = Module.length W (W ⧸ differentIdeal V W)

## 本ファイルの内容(第 7 段)

第 7 段——**再帰不等式から `δₙ → 0` を出す解析の部分**——を閉じる。
残るのは第 2〜6 段(`Module.length` の算術と `p^δ` の指数の評価)である。
-/

namespace ABC3.Found.Falt1

open Filter Topology

/-- **幾何減衰版**。`0 ≤ δₙ` かつ `δ_{n+1} ≤ r·δₙ`(`0 ≤ r < 1`)なら
`δₙ → 0`。 -/
theorem tendsto_zero_of_geometric (δ : ℕ → ℝ) (r : ℝ) (hr0 : 0 ≤ r) (hr1 : r < 1)
    (h0 : ∀ n, 0 ≤ δ n) (hstep : ∀ n, δ (n+1) ≤ r * δ n) :
    Filter.Tendsto δ Filter.atTop (nhds 0) := by
  have hbound : ∀ n, δ n ≤ r ^ n * δ 0 := by
    intro n
    induction n with
    | zero => simp
    | succ k ih =>
      calc δ (k+1) ≤ r * δ k := hstep k
        _ ≤ r * (r ^ k * δ 0) := mul_le_mul_of_nonneg_left ih hr0
        _ = r ^ (k+1) * δ 0 := by ring
  refine squeeze_zero h0 hbound ?_
  have h : Filter.Tendsto (fun n : ℕ => r ^ n) Filter.atTop (nhds 0) :=
    tendsto_pow_atTop_nhds_zero_of_lt_one hr0 hr1
  simpa using h.mul_const (δ 0)

/-- **`Theorem 1.2` の第 7 段**(原文の *"So if ... then
`δ_{n+1} ≤ δₙ − 1/(d+2)`, ... **In any case `δₙ → 0` for `n → ∞`**"*)。

各段で「一定量 `c` 減る」か「`r` 倍以下になる」かの**どちらかが起きれば
十分**である——原文が 2 つの場合分けを与えて *"In any case"* と結んで
いるのはこの形。証明は「単調減少 + 下に有界 ⟹ 収束」してから極限
`L` について `L ≤ max (L − c, r·L)` を取り、`L > 0` なら右辺 `< L` で
矛盾する、というもの。 -/
theorem tendsto_zero_of_two_regime (δ : ℕ → ℝ) (c r : ℝ)
    (hc : 0 < c) (hr1 : r < 1)
    (h0 : ∀ n, 0 ≤ δ n)
    (hstep : ∀ n, δ (n+1) ≤ δ n - c ∨ δ (n+1) ≤ r * δ n) :
    Filter.Tendsto δ Filter.atTop (nhds 0) := by
  have hmax : ∀ n, δ (n+1) ≤ max (δ n - c) (r * δ n) := by
    intro n
    rcases hstep n with h | h
    · exact h.trans (le_max_left _ _)
    · exact h.trans (le_max_right _ _)
  have hanti : Antitone δ := by
    refine antitone_nat_of_succ_le (fun n => ?_)
    refine (hmax n).trans (max_le ?_ ?_)
    · linarith
    · nlinarith [h0 n]
  have hbdd : BddBelow (Set.range δ) := ⟨0, by rintro _ ⟨n, rfl⟩; exact h0 n⟩
  set L := ⨅ i, δ i with hL
  have htend : Filter.Tendsto δ Filter.atTop (nhds L) := tendsto_atTop_ciInf hanti hbdd
  have hL0 : 0 ≤ L := le_ciInf (fun i => h0 i)
  have hsucc : Filter.Tendsto (fun n => δ (n+1)) Filter.atTop (nhds L) :=
    htend.comp (Filter.tendsto_add_atTop_nat 1)
  have hrhs : Filter.Tendsto (fun n => max (δ n - c) (r * δ n)) Filter.atTop
      (nhds (max (L - c) (r * L))) :=
    ((htend.sub_const c).max (htend.const_mul r))
  have hle : L ≤ max (L - c) (r * L) := le_of_tendsto_of_tendsto' hsucc hrhs hmax
  have hLzero : L = 0 := by
    by_contra hne
    have hpos : 0 < L := lt_of_le_of_ne hL0 (Ne.symm hne)
    have h1 : L - c < L := by linarith
    have h2 : r * L < L := by nlinarith
    exact absurd hle (not_le.mpr (max_lt h1 h2))
  rw [← hLzero]
  exact htend

/-- 非空虚性——`d = 0`(すなわち `c = 1/2`, `r = 1/2`)で
`δₙ := (1/2)^n` は仮定を満たし、実際に `0` へ収束する。 -/
example : Filter.Tendsto (fun n : ℕ => (1/2 : ℝ)^n) Filter.atTop (nhds 0) := by
  refine tendsto_zero_of_two_regime _ (1/2) (1/2) (by norm_num) (by norm_num)
    (fun n => by positivity) (fun n => Or.inr ?_)
  rw [pow_succ]
  ring_nf
  rfl

end ABC3.Found.Falt1
