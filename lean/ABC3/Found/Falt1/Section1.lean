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
theorem length_ker_eq_length_coker {R M : Type*} [Ring R] [AddCommGroup M] [Module R M]
    (f : M →ₗ[R] M) (hfin : Module.length R M ≠ ⊤) :
    Module.length R (LinearMap.ker f) = Module.length R (M ⧸ LinearMap.range f) := by
  have h := length_ker_add_target f
  rw [add_comm (Module.length R M)] at h
  exact WithTop.add_right_cancel hfin h

theorem ker_comp_contains_pTorsion {R M N P : Type*} [Ring R] [AddCommGroup M] [Module R M]
    [AddCommGroup N] [Module R N] [AddCommGroup P] [Module R P]
    (f : M →ₗ[R] N) (g : N →ₗ[R] P) (p : R)
    (hg : ∀ y : N, p • y = 0 → g y = 0)
    (x : M) (hx : p • x = 0) : g (f x) = 0 := by
  refine hg (f x) ?_
  rw [← map_smul, hx, map_zero]

/-! ### ★節点 A —— 「第 2 写像の核は `p` 倍の核を含む」(2026-09-05)

原文(物理 p.4 = 印字 p.257、260dpi 確認):

> *The kernel of the second map contains `Ω_{V_{n+1}/Vₙ} ⊗_{V_{n+1}} W_{n+1}`,
> which has `(W_{n+1}/p W_{n+1})^{d+1}` as quotient. **As `Ω_{W_{n+1}/Vₙ}`
> is the direct sum of `d+1` modules of the form `W_{n+1}/p^α W_{n+1}`,
> the kernel of the second map contains the kernel of multiplication by
> `p` on `Ω_{W_{n+1}/Vₙ}`.***

Faltings は「どちらも `(W'/pW')^{d+1}` の大きさだから」という**大きさの
議論**で畳んでいる。ここではそれを、構造定理(巡回分解)を一切使わない
**長さの算術だけ**で完全に閉じる:

1. 任意の有限長 `L` について `length(L[p]) = length(L/pL)`
   ——`0 → L[p] → L →ᵖ L → L/pL → 0` の 2 本の完全列から
   (`length_ker_eq_length_coker` を `f := p·` に適用するだけ)。
2. `L` が `p` で消える `Q` を商に持てば `length Q ≤ length(L/pL) = length(L[p])`。
3. `L[p] ⊆ N[p]` で、`length(N[p]) ≤ length Q ≤ length(L[p])` なら
   有限長の部分加群の比較で `L[p] = N[p]`、ゆえに `N[p] ⊆ L`。

★原典は `Ω_{W_{n+1}/Vₙ}` の**巡回分解**(`d+1` 個の `W'/p^α W'`)を
使うが、上の論法はそれを使わない——必要なのは
`length(N[p]) ≤ length((W'/pW')^{d+1})` という**大きさの比較だけ**である。
この意味で原典より弱い仮定で足りる。 -/

/-- 部分加群 `A ≤ B` で `length B ≤ length A` なら `A = B`。 -/
theorem eq_of_le_of_length_le {R M : Type*} [Ring R] [AddCommGroup M] [Module R M]
    [IsArtinian R M] [IsNoetherian R M]
    {A B : Submodule R M} (h : A ≤ B)
    (hlen : Module.length R ↥B ≤ Module.length R ↥A) : A = B := by
  by_contra hne
  set A' : Submodule R ↥B := Submodule.comap B.subtype A with hA'
  have hA'ne : A' ≠ ⊤ := by
    intro htop
    refine hne (le_antisymm h ?_)
    intro x hxB
    have : (⟨x, hxB⟩ : ↥B) ∈ A' := htop ▸ Submodule.mem_top
    exact this
  have hlt := Submodule.length_lt (R := R) (M := ↥B) hA'ne
  have heq : Module.length R ↥A' = Module.length R ↥A :=
    (Submodule.comapSubtypeEquivOfLe h).length_eq
  rw [heq] at hlt
  exact absurd hlen (not_le.mpr hlt)

/-- `L` が `p` で消える加群 `Q` を商に持つなら `length Q ≤ length(L/pL)`。 -/
theorem length_le_length_quot_smul_of_surjective {R L Q : Type*} [CommRing R]
    [AddCommGroup L] [Module R L] [AddCommGroup Q] [Module R Q]
    (p : R) (hQ : ∀ y : Q, p • y = 0)
    (φ : L →ₗ[R] Q) (hφ : Function.Surjective φ) :
    Module.length R Q ≤ Module.length R (L ⧸ LinearMap.range (LinearMap.lsmul R L p)) := by
  have hker : LinearMap.range (LinearMap.lsmul R L p) ≤ LinearMap.ker φ := by
    rintro _ ⟨x, rfl⟩
    simp only [LinearMap.mem_ker, LinearMap.lsmul_apply, map_smul]
    exact hQ (φ x)
  refine Module.length_le_of_surjective (Submodule.liftQ _ φ hker) ?_
  intro y
  obtain ⟨x, rfl⟩ := hφ y
  exact ⟨Submodule.Quotient.mk x, rfl⟩

/-- **節点 A の代数的核**——`length(N[p]) ≤ length(L/pL)` なら
`N[p] ⊆ L`。構造定理は使わない。 -/
theorem pTorsion_le_of_length_le {R N : Type*} [CommRing R] [AddCommGroup N] [Module R N]
    [IsArtinian R N] [IsNoetherian R N] (p : R) (L : Submodule R N)
    (hge : Module.length R ↥(LinearMap.ker (LinearMap.lsmul R N p))
      ≤ Module.length R (↥L ⧸ LinearMap.range (LinearMap.lsmul R ↥L p))) :
    LinearMap.ker (LinearMap.lsmul R N p) ≤ L := by
  set S : Submodule R ↥L := LinearMap.ker (LinearMap.lsmul R ↥L p) with hS
  set A : Submodule R N := S.map L.subtype with hA
  have hAle : A ≤ LinearMap.ker (LinearMap.lsmul R N p) := by
    rintro _ ⟨x, hx, rfl⟩
    simp only [hS] at hx
    simp only [LinearMap.mem_ker, LinearMap.lsmul_apply, Submodule.subtype_apply]
    exact congrArg Subtype.val hx
  have hlenA : Module.length R ↥A = Module.length R ↥S :=
    ((Submodule.equivMapOfInjective L.subtype Subtype.val_injective S)).length_eq.symm
  have h1 : Module.length R ↥S
      = Module.length R (↥L ⧸ LinearMap.range (LinearMap.lsmul R ↥L p)) :=
    length_ker_eq_length_coker _ (Module.length_ne_top)
  have h2 : Module.length R ↥(LinearMap.ker (LinearMap.lsmul R N p)) ≤ Module.length R ↥A := by
    rw [hlenA, h1]; exact hge
  have h3 := eq_of_le_of_length_le hAle h2
  rw [← h3, hA]
  exact Submodule.map_subtype_le L S

/-- **節点 A 本体**——`L` が `p` で消える `Q` を商に持ち、`Q` が
`N` の `p`-捩れ以上の長さを持つなら、`N[p] ⊆ L`。

原文の `L := Ω_{V_{n+1}/Vₙ} ⊗ W_{n+1}` の像、`Q := (W'/pW')^{d+1}`、
`N := Ω_{W_{n+1}/Vₙ}` に対応する。 -/
theorem ker_contains_pTorsion_of_quotient {R N Q : Type*} [CommRing R]
    [AddCommGroup N] [Module R N] [IsArtinian R N] [IsNoetherian R N]
    [AddCommGroup Q] [Module R Q]
    (p : R) (L : Submodule R N)
    (hQ : ∀ y : Q, p • y = 0)
    (hlen : Module.length R ↥(LinearMap.ker (LinearMap.lsmul R N p)) ≤ Module.length R Q)
    (φ : ↥L →ₗ[R] Q) (hφ : Function.Surjective φ) :
    LinearMap.ker (LinearMap.lsmul R N p) ≤ L :=
  pTorsion_le_of_length_le p L
    (le_trans hlen (length_le_length_quot_smul_of_surjective p hQ φ hφ))

/-- **節点 A を `ker_comp_contains_pTorsion` の入力の形へ**——
`L := ker g` と取ると、原文の *"the composition of the two maps
annihilates the kernel by `p`-multiplication"* の前提
`hg : ∀ y, p·y = 0 → g y = 0` がそのまま出る。 -/
theorem hg_of_quotient {R N P Q : Type*} [CommRing R]
    [AddCommGroup N] [Module R N] [IsArtinian R N] [IsNoetherian R N]
    [AddCommGroup P] [Module R P] [AddCommGroup Q] [Module R Q]
    (p : R) (g : N →ₗ[R] P)
    (hQ : ∀ y : Q, p • y = 0)
    (hlen : Module.length R ↥(LinearMap.ker (LinearMap.lsmul R N p)) ≤ Module.length R Q)
    (φ : ↥(LinearMap.ker g) →ₗ[R] Q) (hφ : Function.Surjective φ) :
    ∀ y : N, p • y = 0 → g y = 0 := by
  intro y hy
  have := ker_contains_pTorsion_of_quotient p (LinearMap.ker g) hQ hlen φ hφ
  exact this (by simpa [LinearMap.mem_ker] using hy)

/-- **原文の literal な形**——`Q = (R/pR)^r`(原文の `(W'/pW')^{d+1}`)。 -/
theorem ker_contains_pTorsion_of_pi_quotient {R N : Type*} [CommRing R]
    [AddCommGroup N] [Module R N] [IsArtinian R N] [IsNoetherian R N]
    (p : R) (r : ℕ) (L : Submodule R N)
    (hlen : Module.length R ↥(LinearMap.ker (LinearMap.lsmul R N p))
      ≤ (r : ℕ∞) * Module.length R (R ⧸ Ideal.span ({p} : Set R)))
    (φ : ↥L →ₗ[R] (Fin r → R ⧸ Ideal.span ({p} : Set R))) (hφ : Function.Surjective φ) :
    LinearMap.ker (LinearMap.lsmul R N p) ≤ L := by
  refine ker_contains_pTorsion_of_quotient p L (fun y => ?_) ?_ φ hφ
  · funext i
    have : p • (y i) = (0 : R ⧸ Ideal.span ({p} : Set R)) := by
      obtain ⟨z, hz⟩ := Submodule.Quotient.mk_surjective (Ideal.span ({p} : Set R)) (y i)
      rw [← hz, ← Submodule.Quotient.mk_smul, Submodule.Quotient.mk_eq_zero, smul_eq_mul]
      exact Ideal.mul_mem_right _ _ (Ideal.subset_span rfl)
    simpa using this
  · have hpi : Module.length R (Fin r → R ⧸ Ideal.span ({p} : Set R))
        = (r : ℕ∞) * Module.length R (R ⧸ Ideal.span ({p} : Set R)) := by
      have e : (Fin r → (R ⧸ Ideal.span ({p} : Set R)))
          ≃ₗ[R] (Fin r →₀ (R ⧸ Ideal.span ({p} : Set R))) :=
        (Finsupp.linearEquivFunOnFinite R _ (Fin r)).symm
      rw [e.length_eq, Module.length_finsupp]
      congr 1
      simp
    rw [hpi]; exact hlen

/-- 非空虚性——`R = ℤ`、`N = ZMod 4`、`p = 2`。
`L := 2·N = {0,2}` は**真の非零部分加群**で、`Q := L` 自身(`2` で消える)を
商に持つ。節点 A の仮定がこの非退化な場合に実際に満たされ、結論
`N[2] ⊆ L` が得られる。 -/
example : LinearMap.ker (LinearMap.lsmul ℤ (ZMod 4) 2)
    ≤ LinearMap.range (LinearMap.lsmul ℤ (ZMod 4) 2) := by
  have hLK : LinearMap.range (LinearMap.lsmul ℤ (ZMod 4) 2)
      = LinearMap.ker (LinearMap.lsmul ℤ (ZMod 4) 2) := by
    ext x
    simp only [LinearMap.mem_range, LinearMap.mem_ker, LinearMap.lsmul_apply]
    revert x
    decide
  refine ker_contains_pTorsion_of_quotient (2 : ℤ)
    (LinearMap.range (LinearMap.lsmul ℤ (ZMod 4) 2))
    (Q := ↥(LinearMap.range (LinearMap.lsmul ℤ (ZMod 4) 2))) (fun y => ?_) ?_
    LinearMap.id Function.surjective_id
  · refine Subtype.ext ?_
    have hy : (y : ZMod 4) ∈ LinearMap.ker (LinearMap.lsmul ℤ (ZMod 4) 2) := hLK ▸ y.2
    simpa [LinearMap.mem_ker] using hy
  · rw [hLK]

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

/-- **`Theorem 1.2` の事実 (c) 本体**——`r` 個の元で生成され `a` で
零化される加群の長さは `r · length(R/(a))` 以下。原文の
*"So its length is at most equal to that of `W_{n+1}` divided by the
`(d+1)`st power of this"* そのもの。 -/
theorem length_le_of_span_and_annihilator {R M : Type*} [CommRing R] [AddCommGroup M]
    [Module R M] (r : ℕ) (g : Fin r → M)
    (hspan : Submodule.span R (Set.range g) = ⊤)
    (a : R) (ha : ∀ m : M, a • m = 0) :
    Module.length R M ≤ (r : ℕ∞) * Module.length R (R ⧸ Ideal.span ({a} : Set R)) := by
  have hker : ∀ i : Fin r, Ideal.span ({a} : Set R)
      ≤ LinearMap.ker (LinearMap.toSpanSingleton R M (g i)) := by
    intro i x hx
    obtain ⟨c, rfl⟩ := Ideal.mem_span_singleton'.mp hx
    show (c * a) • g i = 0
    rw [mul_comm, mul_smul, ha]
  set comp : Fin r → ((R ⧸ Ideal.span ({a} : Set R)) →ₗ[R] M) := fun i =>
    Submodule.liftQ _ (LinearMap.toSpanSingleton R M (g i)) (hker i) with hcomp
  set ψ : (Fin r → (R ⧸ Ideal.span ({a} : Set R))) →ₗ[R] M :=
    ∑ i, (comp i) ∘ₗ (LinearMap.proj i) with hψ
  refine length_le_of_surjective_pi r ψ ?_
  intro m
  have hm : m ∈ Submodule.span R (Set.range g) := hspan ▸ Submodule.mem_top
  obtain ⟨c, hc⟩ := (Submodule.mem_span_range_iff_exists_fun R).mp hm
  refine ⟨fun i => Ideal.Quotient.mk _ (c i), ?_⟩
  rw [hψ]
  simp only [LinearMap.coe_sum, Finset.sum_apply, LinearMap.coe_comp, Function.comp_apply,
    LinearMap.proj_apply, hcomp]
  rw [← hc]
  exact Finset.sum_congr rfl (fun i _ => rfl)

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

/-- 有限個の直積の長さは各因子の長さの和(`Module.length_prod` の
`Fin r` 版、`Fin.consLinearEquiv` で帰納)。構造定理
`Module.equiv_directSum_of_isTorsion` の出力(`⨁ i : ι, R ⧸ R∙p^e`)の
長さを計算するのに使う。 -/
theorem length_pi_fin {R : Type*} [Ring R] :
    ∀ (r : ℕ) (N : Fin r → Type*) (_ : ∀ i, AddCommGroup (N i)) (_ : ∀ i, Module R (N i)),
      Module.length R (∀ i, N i) = ∑ i, Module.length R (N i) := by
  intro r
  induction r with
  | zero =>
    intro N _ _
    haveI : Subsingleton (∀ i : Fin 0, N i) := ⟨fun a b => funext (fun i => absurd i.2 (by omega))⟩
    rw [Module.length_eq_zero]
    simp
  | succ k ih =>
    intro N hN hM
    have e : (∀ i : Fin (k+1), N i) ≃ₗ[R] (N 0 × ∀ i : Fin k, N i.succ) :=
      (Fin.consLinearEquiv (R := R) N).symm
    rw [e.length_eq, Module.length_prod, ih (fun i : Fin k => N i.succ)
      (fun i => hN i.succ) (fun i => hM i.succ), Fin.sum_univ_succ]

/-- **`min(Σkᵢ, e) ≤ Σ min(kᵢ, e)`**——★Faltings の
`β = min{1, δₙ/(d+1)}` は実は **`β = min{1, δₙ}` まで改善できる**
(`(d+1)` で割る必要が無い)。証明は 2 場合:全ての `kᵢ ≤ e` なら
`Σ min = Σ kᵢ`、そうでなければ或る `j` で `min(k_j,e) = e` が和の中に
現れる。

`min_avg_le_sum_min`(平均を使う弱い版、原文どおりの `β`)と
どちらでも `delta_two_regime_of_key` に渡せる。 -/
theorem min_sum_le_sum_min {r : ℕ} (k : Fin r → ℝ) (hk : ∀ i, 0 ≤ k i) (e : ℝ) (he : 0 ≤ e) :
    min (∑ i, k i) e ≤ ∑ i, min (k i) e := by
  by_cases hall : ∀ i, k i ≤ e
  · have hsum : ∑ i, min (k i) e = ∑ i, k i :=
      Finset.sum_congr rfl (fun i _ => min_eq_left (hall i))
    rw [hsum]
    exact min_le_left _ _
  · push Not at hall
    obtain ⟨j, hj⟩ := hall
    have hmj : min (k j) e = e := min_eq_right (le_of_lt hj)
    have hstep : min (k j) e ≤ ∑ i, min (k i) e :=
      Finset.single_le_sum (f := fun i => min (k i) e)
        (fun i _ => le_min (hk i) he) (Finset.mem_univ j)
    rw [hmj] at hstep
    exact le_trans (min_le_right _ _) hstep

/-! ### 事実 (b) を**分解なしで**閉じる

上の 2 つ(`min_avg_le_sum_min`・`min_sum_le_sum_min`)は巡回分解を
経由する道だが、実は**分解を全く使わずに**次が出る:

    min(length M, e) ≤ length(M / 𝔪^e·M)

証明は Nakayama だけ:`𝔪^e·M ≠ 0` なら鎖 `M ⊋ 𝔪M ⊋ ⋯ ⊋ 𝔪^eM` は
**真に**減るので `length(M/𝔪^eM) ≥ e`;`𝔪^e·M = 0` なら
`M/𝔪^eM = M` で `≥ length M`。

★これは Faltings の `β = min{1, δₙ/(d+1)}` より**強い**
(`(d+1)` で割る必要が無い)。構造定理も和因子の個数の評価も要らない。 -/

open IsLocalRing in
/-- **Nakayama**——`𝔪^e·M = 𝔪^{e+1}·M` なら `𝔪^e·M = 0`。 -/
theorem smul_pow_eq_bot_of_eq_succ {R M : Type*} [CommRing R] [IsLocalRing R] [IsNoetherianRing R]
    [AddCommGroup M] [Module R M] [Module.Finite R M] (e : ℕ)
    (h : (maximalIdeal R)^e • (⊤ : Submodule R M)
      = (maximalIdeal R)^(e+1) • (⊤ : Submodule R M)) :
    (maximalIdeal R)^e • (⊤ : Submodule R M) = ⊥ := by
  set N : Submodule R M := (maximalIdeal R)^e • ⊤ with hN
  have hfg : N.FG := Submodule.FG.of_finite
  have hle : N ≤ maximalIdeal R • N := by
    rw [hN] at h ⊢
    rw [pow_succ, mul_comm, mul_smul] at h
    exact le_of_eq h
  exact Submodule.eq_bot_of_le_smul_of_le_jacobson_bot (maximalIdeal R) N hfg hle
    (le_of_eq (IsLocalRing.jacobson_eq_maximalIdeal ⊥ bot_ne_top).symm)

open IsLocalRing in
/-- **`𝔪^e·M ≠ 0` なら `e ≤ length(M/𝔪^e·M)`**——鎖
`M ⊋ 𝔪M ⊋ ⋯ ⊋ 𝔪^eM` が真に減ることから。 -/
theorem le_length_quot_smul_pow {R M : Type*} [CommRing R] [IsLocalRing R] [IsNoetherianRing R]
    [AddCommGroup M] [Module R M] [Module.Finite R M] [IsArtinian R M] [IsNoetherian R M] :
    ∀ (e : ℕ), (maximalIdeal R)^e • (⊤ : Submodule R M) ≠ ⊥ →
      (e : ℕ∞) ≤ Module.length R (M ⧸ (maximalIdeal R)^e • (⊤ : Submodule R M)) := by
  intro e
  induction e with
  | zero => intro _; simp
  | succ k ih =>
    intro hne
    have hST : (maximalIdeal R)^(k+1) • (⊤ : Submodule R M)
        ≤ (maximalIdeal R)^k • (⊤ : Submodule R M) :=
      Submodule.smul_mono_left (Ideal.pow_le_pow_right (Nat.le_succ k))
    have hneT : (maximalIdeal R)^k • (⊤ : Submodule R M) ≠ ⊥ := by
      intro h
      exact hne (le_bot_iff.mp (h ▸ hST))
    have hmapne : Submodule.map ((maximalIdeal R)^(k+1) • (⊤ : Submodule R M)).mkQ
        ((maximalIdeal R)^k • (⊤ : Submodule R M)) ≠ ⊥ := by
      intro h
      have hle : (maximalIdeal R)^k • (⊤ : Submodule R M)
          ≤ (maximalIdeal R)^(k+1) • (⊤ : Submodule R M) := by
        intro x hx
        have hmem : ((maximalIdeal R)^(k+1) • (⊤ : Submodule R M)).mkQ x
            ∈ Submodule.map _ ((maximalIdeal R)^k • (⊤ : Submodule R M)) := ⟨x, hx, rfl⟩
        rw [h, Submodule.mem_bot] at hmem
        exact (Submodule.Quotient.mk_eq_zero _).mp hmem
      have heq : (maximalIdeal R)^k • (⊤ : Submodule R M)
          = (maximalIdeal R)^(k+1) • (⊤ : Submodule R M) := le_antisymm hle hST
      exact hneT (smul_pow_eq_bot_of_eq_succ k heq)
    have hlt := Submodule.length_quotient_lt
      (R := R) (M := M ⧸ (maximalIdeal R)^(k+1) • (⊤ : Submodule R M)) _ hmapne
    rw [(Submodule.quotientQuotientEquivQuotient _ _ hST).length_eq] at hlt
    have hih := ih hneT
    have hfin : ((k : ℕ∞) + 1) ≤ Module.length R
        (M ⧸ (maximalIdeal R)^(k+1) • (⊤ : Submodule R M)) :=
      Order.add_one_le_of_lt (lt_of_le_of_lt hih hlt)
    simpa using hfin

open IsLocalRing in
/-- **`Theorem 1.2` の事実 (b)**——`min(length M, e) ≤ length(M/𝔪^e·M)`。
★Faltings の `β = min{1, δₙ/(d+1)}` より強く、しかも巡回分解も
和因子の個数の評価も要らない。 -/
theorem min_le_length_quot_smul_pow {R M : Type*} [CommRing R] [IsLocalRing R]
    [IsNoetherianRing R] [AddCommGroup M] [Module R M] [Module.Finite R M]
    [IsArtinian R M] [IsNoetherian R M] (e : ℕ) :
    min (Module.length R M) (e : ℕ∞)
      ≤ Module.length R (M ⧸ (maximalIdeal R)^e • (⊤ : Submodule R M)) := by
  by_cases h : (maximalIdeal R)^e • (⊤ : Submodule R M) = ⊥
  · rw [h, (Submodule.quotEquivOfEqBot (⊥ : Submodule R M) rfl).length_eq]
    exact min_le_left _ _
  · exact le_trans (min_le_right _ _) (le_length_quot_smul_pow e h)

/-- DVR ではイデアルの冪は全順序なので、和は指数の `min`。 -/
theorem pow_sup_pow {R : Type*} [CommRing R] (I : Ideal R) (k e : ℕ) :
    I^k ⊔ I^e = I^(min k e) := by
  rcases le_total k e with h | h
  · rw [min_eq_left h, sup_eq_left]
    exact Ideal.pow_le_pow_right h
  · rw [min_eq_right h, sup_eq_right]
    exact Ideal.pow_le_pow_right h

/-- `(R/I)/(J·) ≅ R/(I ⊔ J)` を `R`-線型同型として
(mathlib の `DoubleQuot.quotQuotEquivQuotSup` は `RingEquiv`)。 -/
noncomputable def quotQuotLinearEquiv {R : Type*} [CommRing R] (I J : Ideal R) :
    ((R ⧸ I) ⧸ Ideal.map (Ideal.Quotient.mk I) J) ≃ₗ[R] (R ⧸ (I ⊔ J)) :=
  { DoubleQuot.quotQuotEquivQuotSup I J with
    map_smul' := by
      intro r x
      induction x using Submodule.Quotient.induction_on with
      | H y =>
        induction y using Submodule.Quotient.induction_on with
        | H z => rfl }

open IsLocalRing in
/-- **`length((R/𝔪^k)/(𝔪^e·)) = min(k,e)`**——事実 (b) の各巡回因子の計算。
`length_ker_eq_length_coker` と合わせると、巡回加群 `R/𝔪^k` 上の
`𝔪^e` 倍の**核**の長さがちょうど `min(k,e)` だと分かる。 -/
theorem length_quot_quot_pow {R : Type*} [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
    (k e : ℕ) :
    Module.length R ((R ⧸ maximalIdeal R ^ k) ⧸
      Ideal.map (Ideal.Quotient.mk (maximalIdeal R ^ k)) (maximalIdeal R ^ e))
      = (min k e : ℕ) := by
  rw [(quotQuotLinearEquiv (maximalIdeal R ^ k) (maximalIdeal R ^ e)).length_eq,
    pow_sup_pow, IsDiscreteValuationRing.length_quotient_pow_maximalIdeal]

/-! ## 翻訳層——`Module.length`(`ℕ∞`)から実数へ

`thm_1_2_of_length_bounds` は実数で述べてあるが(`δₙ − δ_{n+1}` や
`min{1, δₙ/(d+1)}` に引き算と割り算が要るため)、加群側の 3 事実は
`Module.length`(`ℕ∞`)で出る。その橋を架ける。 -/

/-- 有限な長さを実数として取り出す。 -/
noncomputable def lenR (R : Type*) [Ring R] (M : Type*) [AddCommGroup M] [Module R M] : ℝ :=
  ((Module.length R M).toNat : ℝ)

theorem lenR_nonneg (R : Type*) [Ring R] (M : Type*) [AddCommGroup M] [Module R M] :
    0 ≤ lenR R M := by
  unfold lenR; positivity

theorem lenR_mono {R M N : Type*} [Ring R] [AddCommGroup M] [Module R M]
    [AddCommGroup N] [Module R N] (h : Module.length R M ≤ Module.length R N)
    (hN : Module.length R N ≠ ⊤) : lenR R M ≤ lenR R N := by
  unfold lenR
  exact_mod_cast ENat.toNat_le_toNat h hN

/-- 事実 (a) の翻訳——`length M = length N + length P` を実数へ。 -/
theorem lenR_add_of_eq {R M N P : Type*} [Ring R] [AddCommGroup M] [Module R M]
    [AddCommGroup N] [Module R N] [AddCommGroup P] [Module R P]
    (h : Module.length R M = Module.length R N + Module.length R P)
    (hN : Module.length R N ≠ ⊤) (hP : Module.length R P ≠ ⊤) :
    lenR R M = lenR R N + lenR R P := by
  unfold lenR
  rw [h, ENat.toNat_add hN hP]
  push_cast
  ring

/-- 事実 (b) の翻訳——`min (length M) e ≤ length Q` を実数へ。 -/
theorem lenR_min_le {R M Q : Type*} [Ring R] [AddCommGroup M] [Module R M]
    [AddCommGroup Q] [Module R Q] (e : ℕ)
    (h : min (Module.length R M) (e : ℕ∞) ≤ Module.length R Q)
    (hQ : Module.length R Q ≠ ⊤) :
    min (lenR R M) (e : ℝ) ≤ lenR R Q := by
  by_cases hM : Module.length R M = ⊤
  · have h0 : lenR R M = 0 := by unfold lenR; rw [hM]; simp
    rw [h0, min_eq_left (by positivity : (0:ℝ) ≤ (e:ℝ))]
    exact lenR_nonneg R Q
  · have hmin : ((min (Module.length R M) (e : ℕ∞)).toNat : ℝ)
        = min (lenR R M) (e : ℝ) := by
      unfold lenR
      rcases le_total (Module.length R M) (e : ℕ∞) with hle | hle
      · rw [min_eq_left hle, min_eq_left]
        exact_mod_cast ENat.toNat_le_toNat hle (by simp)
      · rw [min_eq_right hle, min_eq_right, ENat.toNat_coe]
        exact_mod_cast ENat.toNat_le_toNat hle hM
    rw [← hmin]
    unfold lenR
    exact_mod_cast ENat.toNat_le_toNat h hQ

theorem enat_nsmul_ne_top (r : ℕ) {a : ℕ∞} (ha : a ≠ ⊤) : ((r : ℕ∞) * a) ≠ ⊤ := by
  lift a to ℕ using ha
  rw [← Nat.cast_mul]
  exact ENat.coe_ne_top _

theorem enat_toNat_nsmul (r : ℕ) {a : ℕ∞} (ha : a ≠ ⊤) :
    ((r : ℕ∞) * a).toNat = r * a.toNat := by
  lift a to ℕ using ha
  rw [← Nat.cast_mul, ENat.toNat_coe, ENat.toNat_coe]

/-- 事実 (c) の翻訳——`length M ≤ r · length N` を実数へ。 -/
theorem lenR_le_nsmul {R M N : Type*} [Ring R] [AddCommGroup M] [Module R M]
    [AddCommGroup N] [Module R N] (r : ℕ)
    (h : Module.length R M ≤ (r : ℕ∞) * Module.length R N)
    (hN : Module.length R N ≠ ⊤) :
    lenR R M ≤ (r : ℝ) * lenR R N := by
  have hle := ENat.toNat_le_toNat h (enat_nsmul_ne_top r hN)
  rw [enat_toNat_nsmul r hN] at hle
  unfold lenR
  exact_mod_cast hle

open IsLocalRing in
/-- **事実 (b) の使う形**——`a` が `𝔪^e` を生成し、部分加群 `K` が
`a` 倍の核を含むなら `min(length M, e) ≤ length K`。
原文の *"the composition of the two maps annihilates the kernel by
`p`-multiplication ... this kernel has length at least ..."* に対応。
`ker_comp_contains_pTorsion` が `K := ker(g∘f)` に対する仮定 `hK` を与える。 -/
theorem min_le_length_of_torsion_le {R M : Type*} [CommRing R] [IsLocalRing R]
    [IsNoetherianRing R] [AddCommGroup M] [Module R M] [Module.Finite R M]
    [IsArtinian R M] [IsNoetherian R M]
    (e : ℕ) (a : R) (ha : Ideal.span ({a} : Set R) = (maximalIdeal R)^e)
    (K : Submodule R M)
    (hK : LinearMap.ker (a • (LinearMap.id : M →ₗ[R] M)) ≤ K) :
    min (Module.length R M) (e : ℕ∞) ≤ Module.length R K := by
  have hrange : LinearMap.range (a • (LinearMap.id : M →ₗ[R] M))
      = (maximalIdeal R)^e • (⊤ : Submodule R M) := by
    rw [← ha, Submodule.ideal_span_singleton_smul, ← Submodule.map_top]
    rfl
  have hker : Module.length R (LinearMap.ker (a • (LinearMap.id : M →ₗ[R] M)))
      = Module.length R (M ⧸ (maximalIdeal R)^e • (⊤ : Submodule R M)) := by
    rw [length_ker_eq_length_coker _ (Module.length_ne_top), hrange]
  have h1 := min_le_length_quot_smul_pow (R := R) (M := M) e
  rw [← hker] at h1
  exact le_trans h1 (Module.length_le_of_injective
    (Submodule.inclusion hK) (Submodule.inclusion_injective hK))

/-! ## 最終配線——加群の 3 事実から `δₙ → 0` へ -/

/-- **加群の 3 事実(実数化済み)から原文の鍵の不等式へ**。
`δₙ·e = length(M₁)`、`δ_{n+1}·e = length(M₃)` と置くと、
(a) 指数の加法性・(b) 核の下界・(c) 余核の上界がそのまま
`min{1,δₙ} − (d+1)(δₙ−δ_{n+1}) ≤ δₙ−δ_{n+1}` を与える。 -/
theorem key_inequality_of_lenR (d : ℕ) (e δn δn1 Lker Lcoker : ℝ) (he : 0 < e)
    (hidx : Lker + δn1 * e = δn * e + Lcoker)
    (hb : min (δn * e) e ≤ Lker)
    (hc : Lcoker ≤ ((d:ℝ)+1) * (δn - δn1) * e) :
    min 1 δn - ((d:ℝ)+1) * (δn - δn1) ≤ δn - δn1 := by
  have hmin : min (δn * e) e = (min 1 δn) * e := by
    rcases le_total δn 1 with h | h
    · rw [min_eq_left (show δn * e ≤ e by nlinarith), min_eq_right h]
    · rw [min_eq_right (show e ≤ δn * e by nlinarith), min_eq_left h, one_mul]
  rw [hmin] at hb
  have hstep2 : (min 1 δn - ((d:ℝ)+1) * (δn - δn1)) * e ≤ (δn - δn1) * e := by
    nlinarith [hb, hidx, hc]
  exact le_of_mul_le_mul_right hstep2 he

open Filter Topology in
/-- **`Theorem 1.2`——加群の 3 事実(実数化済み)から `δₙ → 0`**。

`δₙ·e n = length(M₁ⁿ)`、`δ_{n+1}·e n = length(M₃ⁿ)` と置いたとき:

* `hidx`——指数の加法性
  `length(ker) + length(M₃) = length(M₁) + length(coker)`
  (`length_ker_add_target` を実数化したもの)
* `hb`——`min(length M₁, e) ≤ length(ker)`
  (`min_le_length_quot_smul_pow` 由来、★原文より強い形)
* `hc`——`length(coker) ≤ (d+1)(δₙ−δ_{n+1})·e`
  (`length_le_of_span_and_annihilator` 由来)

★これで `Theorem 1.2` の証明は、**実際の `Ω`・差積・塔をこの 3 つに
当てはめるだけ**になった。 -/
theorem thm_1_2_of_module_facts (d : ℕ) (δ e Lker Lcoker : ℕ → ℝ)
    (h0 : ∀ n, 0 ≤ δ n) (he : ∀ n, 0 < e n)
    (hidx : ∀ n, Lker n + δ (n+1) * e n = δ n * e n + Lcoker n)
    (hb : ∀ n, min (δ n * e n) (e n) ≤ Lker n)
    (hc : ∀ n, Lcoker n ≤ ((d:ℝ)+1) * (δ n - δ (n+1)) * e n) :
    Filter.Tendsto δ Filter.atTop (nhds 0) := by
  refine thm_1_2_tendsto_zero d δ h0 (fun n => ?_)
  have hkey := key_inequality_of_lenR d (e n) (δ n) (δ (n+1)) (Lker n) (Lcoker n)
    (he n) (hidx n) (hb n) (hc n)
  have hd1 : (0:ℝ) < (d:ℝ)+1 := by positivity
  have hdiv : min 1 (δ n / ((d:ℝ)+1)) ≤ min 1 (δ n) := by
    refine min_le_min (le_refl 1) ?_
    rw [div_le_iff₀ hd1]
    nlinarith [h0 n]
  linarith

/-- 非空虚性——`d = 0`、`δₙ = (1/2)^n`、`e ≡ 1`、`Lker = δₙ`、
`Lcoker = δₙ/2`。加群の 3 事実が同時に成り立つ具体例。 -/
example : Filter.Tendsto (fun n : ℕ => (1/2 : ℝ)^n) Filter.atTop (nhds 0) := by
  refine thm_1_2_of_module_facts 0 (fun n => (1/2 : ℝ)^n) (fun _ => 1)
    (fun n => (1/2 : ℝ)^n) (fun n => (1/2 : ℝ)^n / 2)
    (fun n => by positivity) (fun _ => one_pos) (fun n => ?_) (fun n => ?_) (fun n => ?_)
  · rw [pow_succ]; ring
  · have hle : ((1:ℝ)/2)^n ≤ 1 := pow_le_one₀ (by norm_num) (by norm_num)
    rw [mul_one, min_eq_left hle]
  · rw [Nat.cast_zero, zero_add, pow_succ]; ring_nf; linarith

open Filter Topology in
/-- **`Theorem 1.2`——Skeleton の `thm12`(ε-N 形)にそのまま繋いだ形**。

`Skeleton/Falt1/Section1.lean` の
`theorem_1_2 (E : RamificationSetup) := E.thm12` は
`∀ ε > 0, ∃ N, ∀ n ≥ N, δ n < ε` という形をしている。
加群の 3 事実からこれが**証明できる**——つまり `Theorem 1.2` を
`Found` にするには、実際の塔からこの 3 事実を出すだけでよい。 -/
theorem thm_1_2_eps_delta (d : ℕ) (δ e Lker Lcoker : ℕ → ℝ)
    (h0 : ∀ n, 0 ≤ δ n) (he : ∀ n, 0 < e n)
    (hidx : ∀ n, Lker n + δ (n+1) * e n = δ n * e n + Lcoker n)
    (hb : ∀ n, min (δ n * e n) (e n) ≤ Lker n)
    (hc : ∀ n, Lcoker n ≤ ((d:ℝ)+1) * (δ n - δ (n+1)) * e n) :
    ∀ ε : ℝ, 0 < ε → ∃ N : ℕ, ∀ n ≥ N, δ n < ε := by
  intro ε hε
  have htend := thm_1_2_of_module_facts d δ e Lker Lcoker h0 he hidx hb hc
  rw [Metric.tendsto_atTop] at htend
  obtain ⟨N, hN⟩ := htend ε hε
  refine ⟨N, fun n hn => ?_⟩
  have hd := hN n hn
  rw [Real.dist_eq, sub_zero, abs_of_nonneg (h0 n)] at hd
  exact hd

/-! ## 1 ステップ分の橋——加群レベルの入力から実数の 3 事実へ -/

theorem lenR_add_eq_add {R A B C D : Type*} [Ring R]
    [AddCommGroup A] [Module R A] [AddCommGroup B] [Module R B]
    [AddCommGroup C] [Module R C] [AddCommGroup D] [Module R D]
    (h : Module.length R A + Module.length R B = Module.length R C + Module.length R D)
    (hA : Module.length R A ≠ ⊤) (hB : Module.length R B ≠ ⊤)
    (hC : Module.length R C ≠ ⊤) (hD : Module.length R D ≠ ⊤) :
    lenR R A + lenR R B = lenR R C + lenR R D := by
  unfold lenR
  have h1 : (Module.length R A + Module.length R B).toNat
      = (Module.length R A).toNat + (Module.length R B).toNat := ENat.toNat_add hA hB
  have h2 : (Module.length R C + Module.length R D).toNat
      = (Module.length R C).toNat + (Module.length R D).toNat := ENat.toNat_add hC hD
  have h3 := congrArg ENat.toNat h
  rw [h1, h2] at h3
  exact_mod_cast h3

open IsLocalRing in
/-- **`Theorem 1.2` の 1 ステップ分の橋**——加群レベルの入力から
`thm_1_2_of_module_facts` が要求する実数の 3 事実を作る。

* `hker`——原文の *"the composition of the two maps annihilates the kernel
  by `p`-multiplication"*(`a` は `𝔪^e` の生成元)。
  `ker_comp_contains_pTorsion` がこの仮定を供給する。
* `g`/`hspan`/`hbann`——余核が `r` 個の元で生成され `b` で零化される
  (原文の *"it is clear that `p^{δₙ−δ_{n+1}}` annihilates
  `W_{n+1}/(Wₙ ⊗_{Vₙ} V_{n+1})`"*)。

★ここまでで `Theorem 1.2` は「実際の `Ω`・差積・塔をこの入力に
当てはめる」だけになる。 -/
theorem step_facts_of_modules {R M N : Type*} [CommRing R] [IsLocalRing R] [IsNoetherianRing R]
    [AddCommGroup M] [Module R M] [Module.Finite R M] [IsArtinian R M] [IsNoetherian R M]
    [AddCommGroup N] [Module R N] [IsArtinian R N] [IsNoetherian R N]
    (F : M →ₗ[R] N) (e : ℕ) (a : R) (ha : Ideal.span ({a} : Set R) = (maximalIdeal R)^e)
    (hker : LinearMap.ker (a • (LinearMap.id : M →ₗ[R] M)) ≤ LinearMap.ker F)
    (r : ℕ) (g : Fin r → (N ⧸ LinearMap.range F))
    (hspan : Submodule.span R (Set.range g) = ⊤)
    (b : R) (hbann : ∀ x : (N ⧸ LinearMap.range F), b • x = 0)
    (hbfin : Module.length R (R ⧸ Ideal.span ({b} : Set R)) ≠ ⊤) :
    (lenR R (LinearMap.ker F) + lenR R N
        = lenR R M + lenR R (N ⧸ LinearMap.range F))
    ∧ (min (lenR R M) (e : ℝ) ≤ lenR R (LinearMap.ker F))
    ∧ (lenR R (N ⧸ LinearMap.range F) ≤ (r : ℝ) * lenR R (R ⧸ Ideal.span ({b} : Set R))) := by
  refine ⟨?_, ?_, ?_⟩
  · exact lenR_add_eq_add (length_ker_add_target F)
      Module.length_ne_top Module.length_ne_top Module.length_ne_top Module.length_ne_top
  · exact lenR_min_le e (min_le_length_of_torsion_le e a ha _ hker) Module.length_ne_top
  · exact lenR_le_nsmul r (length_le_of_span_and_annihilator r g hspan b hbann) hbfin

/-! ## `Theorem 1.2` の残りの手順(2026-09-05 時点の確定版)

解析部分(`thm_1_2_of_length_bounds`)と 3 事実 (a)(b)(c) の**材料**は
本ファイルで揃った。残るのは次の 3 段だけである。

1. **`Ω_{Wₙ/Vₙ}` の巡回分解**。`Wₙ` は DVR(= PID)なので、有限生成
   捩れ加群の構造定理 mathlib `Module.equiv_directSum_of_isTorsion`
   (`M ≃ₗ ⨁ i : ι, R ⧸ R ∙ (p i)^(e i)`、`ι` は `Fintype`)が使える
   ——**在庫にあった**。DVR では既約元は一意化元の同伴なので
   各因子は `Wₙ/𝔪^{e i}`、長さは `e i`。
2. **和因子の個数 `≤ d+1`**。原文(物理 p.4 = 印字 p.257、260dpi 確認)の
   *"It is easy to see that such an `Ω_V` exists, and in fact it can be
   generated by `≤ d+1` elements. This follows from the familiar exact
   sequence `𝔪/𝔪² → Ω_V ⊗_V k → Ω_k → 0`."* から。
   DVR 上では和因子の個数 = `dim_k(M/𝔪M)` なので、生成元の個数で抑える。
3. **(a)(b)(c) を繋ぐ**。(a) は `length_ker_add_target` と
   `length_baseChange_kaehler`、(b) は各巡回因子の `p` 捩れの長さ
   `min(m_i, e)` と `min_avg_le_sum_min`、(c) は
   `length_le_of_surjective_pi` と `length_quot_p_pow`。
-/

end ABC3.Found.Falt1
