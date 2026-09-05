import ABC3.Found.Falt1.Lemma11
import Mathlib.Analysis.SpecificLimits.Basic
import Mathlib.RingTheory.DiscreteValuationRing.Basic
import Mathlib.RingTheory.LocalRing.Length

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

/-! ## 原文の鍵の不等式(260dpi 目視確認、物理 p.5 = 印字 p.258)

★2026-09-05 に原典を 260dpi で描画して逐語確認した。OCR では潰れて
いた末尾の 3 行は正確には次のとおりである:

> *We derive that `δₙ − δ_{n+1} ≥ β − (d+1)(δₙ − δ_{n+1})`.
> So if `δₙ ≥ d+1`, then `δ_{n+1} ≤ δₙ − 1/(d+2)`, and otherwise
> `δ_{n+1} ≤ (1 − 1/((d+1)(d+2)))δₙ`. In any case `δₙ → 0` for `n → ∞`.*

ここで `β = min{1, δₙ/(d+1)}`(こちらも 260dpi で確認)。

★以前の記録では第 2 の場合を `(1 − 1/(d+2))δₙ` と書いていたが、
正しくは **`(1 − 1/((d+1)(d+2)))δₙ`** である(OCR の潰れによる誤り)。

導出は初等的:`Δ := δₙ − δ_{n+1}` と置くと鍵の不等式は `(d+2)Δ ≥ β`。
`δₙ ≥ d+1` なら `β = 1` で `Δ ≥ 1/(d+2)`、そうでなければ
`β = δₙ/(d+1)` で `Δ ≥ δₙ/((d+1)(d+2))`。 -/

/-- **原文の鍵の不等式から 2 領域の再帰へ**。 -/
theorem delta_two_regime_of_key (d : ℕ) (δ : ℕ → ℝ)
    (hkey : ∀ n, min 1 (δ n / ((d:ℝ)+1)) - ((d:ℝ)+1) * (δ n - δ (n+1)) ≤ δ n - δ (n+1))
    (n : ℕ) :
    δ (n+1) ≤ δ n - 1/((d:ℝ)+2)
      ∨ δ (n+1) ≤ (1 - 1/(((d:ℝ)+1)*((d:ℝ)+2))) * δ n := by
  have hd1 : (0:ℝ) < (d:ℝ)+1 := by positivity
  have hd2 : (0:ℝ) < (d:ℝ)+2 := by positivity
  have hk := hkey n
  have hΔ : min 1 (δ n / ((d:ℝ)+1)) ≤ ((d:ℝ)+2) * (δ n - δ (n+1)) := by linarith
  rcases le_or_gt ((d:ℝ)+1) (δ n) with h | h
  · left
    have hb : min 1 (δ n / ((d:ℝ)+1)) = 1 := by
      apply min_eq_left
      rw [le_div_iff₀ hd1]
      linarith
    rw [hb] at hΔ
    have hgoal : 1/((d:ℝ)+2) ≤ δ n - δ (n+1) := by
      rw [div_le_iff₀ hd2]
      nlinarith
    linarith
  · right
    have hb : min 1 (δ n / ((d:ℝ)+1)) = δ n / ((d:ℝ)+1) := by
      apply min_eq_right
      rw [div_le_one hd1]
      linarith
    rw [hb, div_le_iff₀ hd1] at hΔ
    have hexp : (1 - 1/(((d:ℝ)+1)*((d:ℝ)+2))) * δ n
        = δ n - δ n / (((d:ℝ)+1)*((d:ℝ)+2)) := by ring
    have h3 : δ n / (((d:ℝ)+1)*((d:ℝ)+2)) ≤ δ n - δ (n+1) := by
      rw [div_le_iff₀ (by positivity)]
      nlinarith
    linarith [hexp, h3]

/-- **`Theorem 1.2` の解析部分・完成形**——原文の鍵の不等式から
*"In any case `δₙ → 0` for `n → ∞`"* が出る。

★これで `Theorem 1.2` は**その 1 つの不等式**(長さの評価、原文の
第 2〜6 段)に完全に帰着した。残りは:

* `Ω_{V_{n+1}/Vₙ}` が `V_{n+1}/p^a V_{n+1}` 型の `d+1` 個の直和である
  という塔の仮定から、`p` 倍の核が第 2 写像の核に入ること
* `Lemma 1.1`(証明済み)で差積 `δₙ` を `Ω` の余核の長さに翻訳し、
  核の長さ `≥ length(W_{n+1}/p^β W_{n+1})`(`β = min{1, δₙ/(d+1)}`)
* 余核が `p^{δₙ−δ_{n+1}}` で零化される(原文の *"it is clear that"*)
-/
theorem thm_1_2_tendsto_zero (d : ℕ) (δ : ℕ → ℝ) (h0 : ∀ n, 0 ≤ δ n)
    (hkey : ∀ n, min 1 (δ n / ((d:ℝ)+1)) - ((d:ℝ)+1) * (δ n - δ (n+1)) ≤ δ n - δ (n+1)) :
    Filter.Tendsto δ Filter.atTop (nhds 0) := by
  refine tendsto_zero_of_two_regime δ (1/((d:ℝ)+2)) (1 - 1/(((d:ℝ)+1)*((d:ℝ)+2)))
    (by positivity) ?_ h0 (fun n => delta_two_regime_of_key d δ hkey n)
  have hp : (0:ℝ) < 1/(((d:ℝ)+1)*((d:ℝ)+2)) := by positivity
  linarith

/-- 非空虚性——`d = 0`、`δₙ := (1/2)^n` は原文の鍵の不等式を
**等号で**満たす(`β = min{1, δₙ} = δₙ`、`δₙ − δ_{n+1} = δₙ/2`)。
仮定が空虚に真になっていないことの対照。 -/
example : Filter.Tendsto (fun n : ℕ => (1/2 : ℝ)^n) Filter.atTop (nhds 0) := by
  refine thm_1_2_tendsto_zero 0 _ (fun n => by positivity) (fun n => ?_)
  have hle : ((1:ℝ)/2)^n ≤ 1 := pow_le_one₀ (by norm_num) (by norm_num)
  have hmin : min 1 (((1:ℝ)/2)^n / (((0:ℕ):ℝ) + 1)) = ((1:ℝ)/2)^n := by
    rw [Nat.cast_zero, zero_add, div_one]
    exact min_eq_right hle
  rw [hmin, Nat.cast_zero, zero_add, one_mul, pow_succ]
  ring_nf
  linarith

/-! ## 第 2〜6 段(長さの評価)の骨組み

残るのは原文の

> *... the kernel of the second map contains the kernel of multiplication
> by `p` on `Ω_{W_{n+1}/Vₙ}`, and hence **the composition of the two maps
> annihilates the kernel by `p`-multiplication** on `Ω_{Wₙ/Vₙ} ⊗_{Wₙ}
> W_{n+1}`. ... this kernel has **length at least** equal to that of
> `W_{n+1}/p^β W_{n+1}` ... the cokernel of the composition of the two maps
> is annihilated by `p^{δₙ−δ_{n+1}}`. So its **length is at most** ...*

の部分である。ここで効くのは「核と余核の長さの差が源と標的の長さの差に
等しい」という**指数の加法性**で、それを先に置いておく。 -/

/-- **長さの「指数」加法性**——任意の線形写像 `f : M → N` について

    length(ker f) + length N = length M + length(coker f)

`Theorem 1.2` の長さの評価はこの恒等式に核・余核の上下からの評価を
入れて `δₙ − δ_{n+1}` の不等式を出す、という形をしている。 -/
theorem length_ker_add_target {R M N : Type*} [Ring R] [AddCommGroup M] [Module R M]
    [AddCommGroup N] [Module R N] (f : M →ₗ[R] N) :
    Module.length R (LinearMap.ker f) + Module.length R N
      = Module.length R M + Module.length R (N ⧸ LinearMap.range f) := by
  have hM : Module.length R M
      = Module.length R (LinearMap.ker f) + Module.length R (LinearMap.range f) := by
    refine Module.length_eq_add_of_exact (LinearMap.ker f).subtype f.rangeRestrict
      Subtype.val_injective f.surjective_rangeRestrict ?_
    intro x
    constructor
    · intro hx
      exact ⟨⟨x, by simpa [LinearMap.mem_ker] using congrArg Subtype.val hx⟩, rfl⟩
    · rintro ⟨y, rfl⟩
      exact Subtype.ext y.2
  have hN : Module.length R N
      = Module.length R (LinearMap.range f) + Module.length R (N ⧸ LinearMap.range f) := by
    refine Module.length_eq_add_of_exact (LinearMap.range f).subtype (LinearMap.range f).mkQ
      Subtype.val_injective (Submodule.mkQ_surjective _) ?_
    intro x
    constructor
    · intro hx
      exact ⟨⟨x, (Submodule.Quotient.mk_eq_zero _).mp hx⟩, rfl⟩
    · rintro ⟨y, rfl⟩
      exact (Submodule.Quotient.mk_eq_zero _).mpr y.2
  rw [hM, hN]
  ring

/-- **原文の *"the composition of the two maps annihilates the kernel by
`p`-multiplication"***。第 2 写像の核が `p` 倍の核を含むなら、合成の核は
源の `p` 倍の核を含む。 -/
theorem ker_comp_contains_pTorsion {R M N P : Type*} [Ring R] [AddCommGroup M] [Module R M]
    [AddCommGroup N] [Module R N] [AddCommGroup P] [Module R P]
    (f : M →ₗ[R] N) (g : N →ₗ[R] P) (p : R)
    (hg : ∀ y : N, p • y = 0 → g y = 0)
    (x : M) (hx : p • x = 0) : g (f x) = 0 := by
  refine hg (f x) ?_
  rw [← map_smul, hx, map_zero]

/-! ### `length(W/p^α W)` の線型性

原文(物理 p.4 = 印字 p.257、260dpi 確認)は

> *The valuations `v : K* → Q` will be normalized such that `v(p) = 1`.
> If `α ∈ v(K*)` we denote by `p^α V` the ideal of `V` generated by an
> element of valuation `α`.*

と `p^α`(`α ∈ ℚ`)を定義する。長さの評価を `δₙ` の不等式に翻訳するのに
必要なのは `α ↦ length(W/p^α W)` の**線型性**である。その整数版は
mathlib の `IsDiscreteValuationRing.length_quotient_pow_maximalIdeal`
(`length(R/m^n) = n`)から直ちに出る。 -/

open IsLocalRing in
/-- **`length(W/p^k W) = k·e`**——`Ideal.span {p} = m^e`(`e` は絶対分岐指数)
のとき。`α ↦ length(W/p^α W)` の線型性の整数版。 -/
theorem length_quot_p_pow {W : Type*} [CommRing W] [IsDomain W] [IsDiscreteValuationRing W]
    (p : W) (e : ℕ) (hp : Ideal.span {p} = maximalIdeal W ^ e) (k : ℕ) :
    Module.length W (W ⧸ Ideal.span ({p^k} : Set W)) = (k : ℕ∞) * (e : ℕ∞) := by
  have hspan : Ideal.span ({p^k} : Set W) = maximalIdeal W ^ (e * k) := by
    rw [← Ideal.span_singleton_pow, hp, ← pow_mul]
  rw [hspan, IsDiscreteValuationRing.length_quotient_pow_maximalIdeal]
  push_cast
  ring

/-! ## `Theorem 1.2` を「長さの 3 事実」に完全に帰着させる

上の 3 つの道具(指数の加法性・`p` 倍の核の移送・長さの線型性)を使うと、
原文の証明は次の 3 つの数値的な事実に帰着する:

* **`hidx`**——`length(ker) − length(coker) = (δₙ − δ_{n+1})·e`
  (指数の加法性 `length_ker_add_target` と `Lemma 1.1` から)
* **`hk`**——*"this kernel has length at least equal to that of
  `W_{n+1}/p^β W_{n+1}`"*(`β = min{1, δₙ/(d+1)}`)
* **`hc`**——*"the cokernel ... is annihilated by `p^{δₙ−δ_{n+1}}`. So its
  length is at most equal to that of `W_{n+1}` divided by the `(d+1)`st
  power of this"*

この 3 つから原文の鍵の不等式が出て、あとは既に証明済みの
`delta_two_regime_of_key`・`tendsto_zero_of_two_regime` で `δₙ → 0`。 -/

/-- **3 つの長さの事実から原文の鍵の不等式へ**。 -/
theorem key_inequality_of_length_bounds (d : ℕ) (e β Δ Lker Lcoker : ℝ) (he : 0 < e)
    (hidx : Lker - Lcoker = Δ * e)
    (hk : β * e ≤ Lker)
    (hc : Lcoker ≤ ((d:ℝ)+1) * Δ * e) :
    β - ((d:ℝ)+1) * Δ ≤ Δ := by
  have h1 : (β - ((d:ℝ)+1) * Δ) * e ≤ Δ * e := by nlinarith
  exact le_of_mul_le_mul_right (by linarith [h1]) he

/-- **`Theorem 1.2`——長さの 3 事実に完全に帰着した形**。

★これが「解析部分は全部済んだ」という到達点である——残るのは上の
`hidx`・`hk`・`hc` を実際の `Ω` と差積について示すことだけ。 -/
theorem thm_1_2_of_length_bounds (d : ℕ) (δ : ℕ → ℝ) (h0 : ∀ n, 0 ≤ δ n)
    (e Lker Lcoker : ℕ → ℝ) (he : ∀ n, 0 < e n)
    (hidx : ∀ n, Lker n - Lcoker n = (δ n - δ (n+1)) * e n)
    (hk : ∀ n, min 1 (δ n / ((d:ℝ)+1)) * e n ≤ Lker n)
    (hc : ∀ n, Lcoker n ≤ ((d:ℝ)+1) * (δ n - δ (n+1)) * e n) :
    Filter.Tendsto δ Filter.atTop (nhds 0) := by
  refine thm_1_2_tendsto_zero d δ h0 (fun n => ?_)
  have h := key_inequality_of_length_bounds d (e n) (min 1 (δ n / ((d:ℝ)+1)))
    (δ n - δ (n+1)) (Lker n) (Lcoker n) (he n) (hidx n) (hk n) (hc n)
  linarith

/-- 非空虚性——`d = 0`、`δₙ = (1/2)^n`、`e ≡ 1`、`Lker = δₙ`、
`Lcoker = δₙ/2`。3 つの長さの事実が同時に成り立つ具体例。 -/
example : Filter.Tendsto (fun n : ℕ => (1/2 : ℝ)^n) Filter.atTop (nhds 0) := by
  refine thm_1_2_of_length_bounds 0 (fun n => (1/2 : ℝ)^n) (fun n => by positivity)
    (fun _ => 1) (fun n => (1/2 : ℝ)^n) (fun n => (1/2 : ℝ)^n / 2)
    (fun _ => one_pos) (fun n => ?_) (fun n => ?_) (fun n => ?_)
  · rw [pow_succ]; ring
  · have hle : ((1:ℝ)/2)^n ≤ 1 := pow_le_one₀ (by norm_num) (by norm_num)
    rw [Nat.cast_zero, zero_add, div_one, min_eq_right hle, mul_one]
  · rw [Nat.cast_zero, zero_add, pow_succ]; ring_nf; linarith

/-! ### 事実 (a) の材料——`Lemma 1.1` を基底変換する

事実 (a)(`length(ker) − length(coker) = (δₙ − δ_{n+1})·e`)は

* 指数の加法性 `length_ker_add_target`(上)
* 源 `M₁ = Ω_{Wₙ/Vₙ} ⊗_{Wₙ} W_{n+1}` の長さ = `δₙ·e`
* 標的 `M₃ = Ω_{W_{n+1}/V_{n+1}}` の長さ = `δ_{n+1}·e`

の 3 つで出る。真ん中が本節の `length_baseChange_kaehler` で、
`Lemma 1.1`(証明済み)と mathlib の `IsLocalRing.length_baseChange`
(**在庫にあった**)を繋ぐだけである。標的の方は `Lemma 1.1` を
そのまま `(V_{n+1}, W_{n+1})` に適用すればよい。 -/

open IsLocalRing in
/-- **`Lemma 1.1` を基底変換した形**。
`length_{W'}(W' ⊗_W Ω[W⁄V]) = length_W(W/𝔡) · length_{W'}(W'/𝔪_W W')`。
右の因子が `α ↦ length(W'/p^α W')` の比例定数 `e` に当たる。 -/
theorem length_baseChange_kaehler {Z V K L W W' : Type*} [CommRing Z] [CommRing V]
    [IsDedekindDomain V] [Field K] [Algebra V K] [IsFractionRing V K] [Field L] [Algebra K L]
    [FiniteDimensional K L] [Algebra.IsSeparable K L] [CommRing W] [Algebra W L] [Algebra V W]
    [Algebra V L] [IsScalarTower V K L] [IsScalarTower V W L] [IsIntegralClosure W V L]
    [IsDedekindDomain W] [Module.IsTorsionFree V W] [Algebra Z V] [Algebra Z W]
    [IsScalarTower Z V W]
    [IsLocalRing W] [CommRing W'] [IsLocalRing W'] [Algebra W W']
    [IsLocalHom (algebraMap W W')] [Module.Flat W W']
    (w : W) (hint : IsIntegral V w) (hadjoin : Algebra.adjoin V ({w} : Set W) = ⊤)
    (hw : Algebra.adjoin K ({(algebraMap W L) w} : Set L) = ⊤) :
    Module.length W' (TensorProduct W W' Ω[W⁄V])
      = Module.length W (W ⧸ differentIdeal V W)
        * Module.length W' (W' ⧸ Ideal.map (algebraMap W W') (maximalIdeal W)) := by
  rw [IsLocalRing.length_baseChange W W' Ω[W⁄V],
    (lemma_1_1_falt1 (Z := Z) w hint hadjoin hw).2]

/-! ### 事実 (c) の材料——「`d+1` 個の生成元」の効き方

原文の *"So its length is at most equal to that of `W_{n+1}` divided by
the `(d+1)`st power of this"* は、余核が `d+1` 個の生成元を持ち
`p^{δₙ−δ_{n+1}}` で消えることから

    length(coker) ≤ (d+1)·length(W_{n+1}/p^{δₙ−δ_{n+1}}W_{n+1})

となる、という形をしている(`d+1` は `Ω_V` が `≤ d+1` 個の元で生成される
——物理 p.4 = 印字 p.257、260dpi 確認——ことに由来する)。 -/

/-- 有限直積の長さ。 -/
theorem length_pi {R N : Type*} [Ring R] [AddCommGroup N] [Module R N] (r : ℕ) :
    Module.length R (Fin r → N) = (r : ℕ∞) * Module.length R N := by
  have e : (Fin r → N) ≃ₗ[R] (Fin r →₀ N) := (Finsupp.linearEquivFunOnFinite R N (Fin r)).symm
  rw [e.length_eq, Module.length_finsupp]
  congr 1
  simp

/-- **`r` 個の生成元で書ける加群の長さは `r·length(N)` 以下**
——事実 (c) の形(`N := W'/p^{δₙ−δ_{n+1}}W'`、`r := d+1`)。 -/
theorem length_le_of_surjective_pi {R M N : Type*} [Ring R] [AddCommGroup M] [Module R M]
    [AddCommGroup N] [Module R N] (r : ℕ) (f : (Fin r → N) →ₗ[R] M)
    (hf : Function.Surjective f) :
    Module.length R M ≤ (r : ℕ∞) * Module.length R N := by
  rw [← length_pi (R := R) (N := N) r]
  exact Module.length_le_of_surjective f hf

/-! ### 事実 (b) の材料——`β = min{1, δₙ/(d+1)}` の由来

原文は *"this kernel has length at least equal to that of
`W_{n+1}/p^β W_{n+1}`, where `β = min{1, δₙ/(d+1)}`"* と言う。

`Ω_{Wₙ/Vₙ} ⊗ W_{n+1}` が `r ≤ d+1` 個の巡回加群 `W'/π^{m_i}W'`
(`Σ m_i = δₙ·e`)の直和で、`p = (unit)·π^e` とすると:

* 巡回加群 `W'/π^m W'` の `p` 捩れ部分の長さは `min(m, e)`
* 直和なので全体では `Σ min(m_i, e)`
* **最大の `m_j` が平均 `(Σ m_i)/r` 以上**なので
  `Σ min(m_i, e) ≥ min((Σ m_i)/r, e) = min(δₙ·e/r, e) = min(δₙ/r, 1)·e`

——最後が `β = min{1, δₙ/(d+1)}` の正体である。その数値的な核を置く。 -/

/-- **`β = min{1, δₙ/(d+1)}` の由来**——最大の `m_j` が平均以上だから
`Σ min(m_i, e) ≥ min((Σ m_i)/r, e)`。 -/
theorem min_avg_le_sum_min {r : ℕ} (hr : 0 < r) (m : Fin r → ℝ) (hm : ∀ i, 0 ≤ m i)
    (e : ℝ) (he : 0 ≤ e) :
    min ((∑ i, m i) / (r:ℝ)) e ≤ ∑ i, min (m i) e := by
  haveI : NeZero r := ⟨hr.ne'⟩
  obtain ⟨j, -, hj⟩ := Finset.exists_max_image Finset.univ m ⟨0, Finset.mem_univ 0⟩
  have hsum : ∑ i, m i ≤ (r:ℝ) * m j := by
    calc ∑ i : Fin r, m i ≤ ∑ _i : Fin r, m j :=
          Finset.sum_le_sum (fun i _ => hj i (Finset.mem_univ i))
      _ = (r:ℝ) * m j := by simp [Finset.sum_const, mul_comm]
  have hdiv : (∑ i, m i) / (r:ℝ) ≤ m j := by
    rw [div_le_iff₀ (by exact_mod_cast hr)]
    linarith
  have h1 : min ((∑ i, m i) / (r:ℝ)) e ≤ min (m j) e := min_le_min hdiv (le_refl e)
  have h2 : min (m j) e ≤ ∑ i, min (m i) e := by
    refine Finset.single_le_sum (f := fun i => min (m i) e) (fun i _ => ?_) (Finset.mem_univ j)
    exact le_min (hm i) he
  linarith

end ABC3.Found.Falt1
