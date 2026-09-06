import ABC3.Found.PGC.ConjugateSumValuation

/-!
# π 進展開と `1 + v(σα_n − α_n) = n + i(σ^n)`(Yoshida 2008 Lemma 6.5)

典拠: T. Yoshida, *Local Class Field Theory via Lubin-Tate Theory*
(Ann. Fac. Sci. Toulouse 17-2, 2008; arXiv math/0606108) **Lemma 6.5**(物理 p.14)。
構造化済み原文は `ResearchPaper/1_Structured/Local Class Field Theory via Lubin-Tate Theory/`
`section-6.html` の `id="lemma-6-5"`(`data-pdf-page="14"`, `data-item="Lemma 6.5"`)。原文は

> Lemma 6.5. Let σ ∈ G_1. For each n ∈ Z, there exists α ∈ K'^× such that v(α) = n and
> v(σ(α) − α) = n + i(σ^n). Moreover, any x ∈ K'^× can be written as a sum
> x = Σ^∞_{n=v(x)} x_n (see Appendix I) where each x_n satisfies above two properties for n
> if x_n ≠ 0.

★**原文を `.txt` で読んではならない**: `pdftotext` は総和記号 `Σ` と `≠` を出力しない
(構造化 HTML が `data-txt=""` / `data-txt="="` で明示している)。

## ★★★★原典の字面は偽である —— `−1` のずれ(erratum)

`Def 6.1` の `i(σ) = v(σπ − π)` の下では、原文の `v(σα − α) = n + i(σ^n)` は成り立たない。
正しくは

```
v(σ α_n − α_n) = n + i(σ^n) − 1        (α_n := ∏_{i<n} σ^i π)
```

である。★`n = 1` が決定的で、そこでは `α_1 = π` だから左辺は `v(σπ − π) = i(σ)` であり、
原文の右辺 `1 + i(σ)` とは 1 ずれる。★2 度独立に検算した。

★★**本ファイルは `−1` を書かない。** 代わりに **環の恒等式**

```
π · (σ • α_n − α_n) = α_n · (σ^n • π − π)        (`uniformizer_mul_smul_sub`)
```

から `addVal` を取って **`1 + v(σ•α_n − α_n) = n + i(σ^n)`**
(`one_add_ramIndex_uniformizerProd`)として述べる。この形なら

* `ℕ∞` の**引き算が一切出ない**(切り詰め `⊤ − 1 = ⊤` による空虚な真を作らない)、
* `σ = 1` でも `n = 0` でも両辺 `⊤` / 恒等式として**込みで真**、
* erratum の追加コストが実質 0 行、

になる。★**逸脱**: 原典の等式を移項した形で述べている(数学的内容は同値)。

## 何を証明したか

**§1 抽象核(分岐・付値・Galois の語彙が 1 つも出てこない)**

* `uniformizerProd σ π n := ∏ i ∈ Finset.range n, σ ^ i • π`。
* `uniformizer_mul_smul_sub : π * (σ • α_n − α_n) = α_n * (σ^n • π − π)`。
  望遠鏡積 `Finset.prod_range_succ` / `Finset.prod_range_succ'` の 1 行。
  ★一般の可換環と一般の `MulSemiringAction` で成り立つ。`π` が素元である必要すらない。

**§2 付値の不変性**

* `maximalIdeal_eq_span_smul : 𝔪 = (π) → 𝔪 = (σ•π)`(環自己同型は素元を素元に写す)。
  ★`span_smul_uniformizer`(`RamificationJumpDivisibility.lean`)は `τ ∈ G_0` を要求するが、
  実は**任意の群元で成り立つ**(`σ⁻¹` を当てて割り算を写すだけ)。
* `addVal_smul : v(σ•x) = v(x)`。
* `addVal_uniformizerProd : v(α_n) = n`。

**§3 Lemma 6.5 の第 1 主張**

* `ramIndex α σ := addVal B (σ • α − α)`——原典 `Def 6.1` の `i(σ)`(★§6 も参照)。
* `one_add_ramIndex_uniformizerProd : 1 + i_{α_n}(σ) = n + i_π(σ^n)`。
* `exists_addVal_eq_and_one_add_ramIndex`——存在形(`α ≠ 0`、`v(α) = n` 込み)。

★**`σ ∈ G_1` は仮定しない。** 原文は Lemma 6.5 全体を `σ ∈ G_1` の下で述べているが、
第 1 主張は `σ` が任意でも成り立つ(仮定を落とした = 主張を強めた方向の逸脱)。

**§4 加法付値の有限和(抽象核、一般の `Γ`)**

* `lt_map_sum` / `map_sum_eq_of_lt`——付値が一意に最小の項があれば和の付値はその項の付値。
  ★mathlib に **`AddValuation` 版が無い**(乗法版 `Valuation.map_sum_eq_of_lt` はある)。
  `AddValuation` は手書きミラーなので真の欠落である。一般の
  `LinearOrderedAddCommMonoidWithTop` で書いてあるので `ℕ∞` に限らない。

**§5 打ち切り π 進展開(Lemma 6.5 の第 2 主張)**

* `exists_digits`——`x ∈ 𝔪^N`、`N ≤ M` なら桁 `a : ℕ → A` が取れて
  `x − Σ_{n ∈ [N,M)} a_n α_n ∈ 𝔪^M`、かつ各桁は **`0` か単元**。
* `smul_digitSum`——桁を `𝒪_K` から取ると `σ` で動かない(`smul_algebraMap`)ので、
  有限和への `σ` の作用は `α_n` への作用に落ちる(下流が使う項別作用)。
* `addVal_mul_uniformizerProd : v(u · α_n) = n`(`u` 単元)——桁の付値が相異なることの根拠。

**§6 `i` と `G_n` の橋**

* `mem_lowerRamificationGroup_iff_lt_ramIndex : σ ∈ G_n ↔ n < i(σ)`。
* `ramIndex_eq_iff_mem_sdiff : i(σ) = n+1 ↔ (σ ∈ G_n ∧ σ ∉ G_{n+1})`。
  ★この木は **Y1 が `addVal` 言語、Y3/Y4 が `ramCoeff` 言語**という 2 つの平行な言語に
  なっている。この橋があると下流が両方から型で引ける。

**§7 `PAdicLocalField` 具体層**

`B := adjoinIntegers K x`、`A := 𝒪[K.carrier]`、`G := Gal(K(x)/K)`。
剰余系 `hres` は完全分岐から
`exists_sub_algebraMap_mem_maximalIdeal_of_isTotallyRamifiedAdjoin` が供給する。

## ★π 進展開について落としたもの(意図的、下流で不要と測定済み)

| 原文が使うもの | 本ファイル |
|---|---|
| 展開の**存在** | ★入れた(`exists_digits`) |
| 展開の**一意性** | 落とした(`Lemma 5.2(i)` だけが使う) |
| **無限和・収束・完備性**(Appendix I) | 落とした(消費側は有限の付値を否定するだけで、`M := v(z)+1` で打ち切れる) |
| **`σ` の項別作用の連続性** | 落とした(有限和なら加法性だけ。`smul_digitSum`) |
| **`C = {0} ∪ μ_{q−1}`(Teichmüller 代表系)** | 落とした。代表系は `𝒪_K` から取る |

★代表系を `𝒪_K` から取ると **`σ` 不変が `smul_algebraMap` で無料**になる(完全分岐なので
`𝒪_K → 𝒪_{K'}/𝔭` は全射で、代表系として十分)。落とせないのは「**零または単元**」だけである。

## ★逸脱の記録

1. **原典の等式を掛け算形に移項した**(上の erratum 節)。原文 `v(σα−α) = n + i(σ^n)` は
   `Def 6.1` の下では偽で、正しくは `−1` ずれる。両辺に `1` を足した
   `1 + v(σα_n−α_n) = n + i(σ^n)` の形で述べる。
2. **`n ∈ Z` を `n : ℕ` に制限した**。原文は `n` が負でもよい(`K'^×` の元を取る)が、
   本ファイルの舞台は付値環 `B` であって商体ではない。負の `n` は `α_{-n}` の逆元を
   取るだけで、**落としたのは配線であって数学ではない**(`ConjugateSumValuation` の
   逸脱 1 と同じ理由)。
3. **`σ ∈ G_1` を仮定していない**(§3・§5 とも `σ` は任意の群元)。
4. **第 2 主張(展開)を「打ち切り」の形にした**。原文の `x = Σ^∞_{n=v(x)} x_n` は
   Appendix I の収束を要求するが、消費側(`Prop 6.6 (iii)` の背理法)は有限個の桁しか
   見ないので、`𝔪^M` を法とする有限和で足りる。★無限和・完備性は入れていない。
5. **展開の一意性を証明していない**。`Lemma 5.2(i)` 以外に消費者が無い。
6. `hres`(剰余系の存在)を仮定として持つ抽象形にし、具体層で完全分岐から供給した。

## ★退化の自己検査

* **(D1) `v(α_n) = n` を落とすと §5 は偽**。`α_n := 0` と置くと結論は「`x ∈ 𝔪^M`」に
  縮み、`x := 1`、`N := 0`、`M := 1` で反証される(`1 ∉ 𝔪`)。
  `addVal_uniformizerProd` が §5 の中で `exists_mul_eq_of_addVal_le` を起動する鍵である。
* **(D2) 桁の `0-or-unit` を落とすと項の付値がずれる**。`B := ℤ_p`、`α_n := p^n`、
  桁 `a_0 := p` を許すと `a_0 α_0 = p` と `a_1 α_1` の付値が衝突し、
  消費側の「相異なる付値」が壊れる。★単元性は **`B` 側**で言っている
  (`IsUnit (algebraMap A B (a n))`)ので `IsLocalHom` を要求しない。
* **(D3) `α_n := π^n` で §3 を主張すると偽**。`σ^n = 1` だが `σ(π^n) ≠ π^n` という状況で
  右辺は `n + ⊤ = ⊤`、左辺は有限になる。★望遠鏡積 `∏ σ^i π` でなければならない。
  ☆具体例(そういう `K'/K` と `σ`)の構成は高いので作っていない。反証の形だけ記す。
* **`σ` を `MulSemiringAction` でなく単なる加法射にすると §1 が偽**。望遠鏡は
  `σ • ∏ = ∏ σ •`(`Finset.smul_prod'`)に乗っている。
* **§4 で `j` の付値が一意最小でないと偽**。`v(f i) = v(f j)` を許すと
  `f i + f j` で付値が上がりうる(`p` と `−p` in `ℤ_p`)。
-/

namespace ABC3.Found.PGC

open IsLocalRing ABC3.Skeleton.PGC IsDiscreteValuationRing
open scoped NNReal Valued

def one_add_ramIndex_uniformizerProd.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 14, item := "Lemma 6.5", sectionId := "lemma-6-5" }

/-! ## §1 抽象核 —— 望遠鏡積

`α_n := ∏_{i<n} σ^i π` に対する恒等式 `π (σα_n − α_n) = α_n (σ^n π − π)`。
★一般の可換環・一般の `MulSemiringAction` で成り立ち、分岐・付値・Galois の語彙は
1 つも現れない。`π` が素元である必要すらない。 -/

/-- **共役の望遠鏡積** `α_n := ∏_{i<n} σ^i • π`。

`α_0 = 1`、`α_1 = π`、`α_{n+1} = α_n · σ^n π`。★`π^n` ではないことが本質(退化検査 D3)。 -/
def uniformizerProd {B : Type*} [CommRing B] {G : Type*} [Monoid G] [MulSemiringAction G B]
    (σ : G) (π : B) (n : ℕ) : B := ∏ i ∈ Finset.range n, σ ^ i • π

theorem uniformizerProd_zero {B : Type*} [CommRing B] {G : Type*} [Monoid G]
    [MulSemiringAction G B] (σ : G) (π : B) : uniformizerProd σ π 0 = 1 := by
  simp [uniformizerProd]

theorem uniformizerProd_one {B : Type*} [CommRing B] {G : Type*} [Monoid G]
    [MulSemiringAction G B] (σ : G) (π : B) : uniformizerProd σ π 1 = π := by
  simp [uniformizerProd]

theorem uniformizerProd_succ {B : Type*} [CommRing B] {G : Type*} [Monoid G]
    [MulSemiringAction G B] (σ : G) (π : B) (n : ℕ) :
    uniformizerProd σ π (n + 1) = uniformizerProd σ π n * σ ^ n • π := by
  simp [uniformizerProd, Finset.prod_range_succ]

/-- `σ • ∏ = ∏ σ •`(`Finset.smul_prod'`)で添字が 1 つずれる。 -/
theorem smul_uniformizerProd {B : Type*} [CommRing B] {G : Type*} [Monoid G]
    [MulSemiringAction G B] (σ : G) (π : B) (n : ℕ) :
    σ • uniformizerProd σ π n = ∏ i ∈ Finset.range n, σ ^ (i + 1) • π := by
  rw [uniformizerProd, Finset.smul_prod']
  exact Finset.prod_congr rfl fun i _ => by rw [smul_smul, ← pow_succ']

/-- **望遠鏡** `π · σα_n = α_n · σ^n π`。

左辺は `π · ∏_{i<n} σ^{i+1} π = ∏_{i<n+1} σ^i π`(`Finset.prod_range_succ'`)、
右辺は同じものを `Finset.prod_range_succ` で切ったもの。 -/
theorem uniformizer_mul_smul_uniformizerProd {B : Type*} [CommRing B] {G : Type*} [Monoid G]
    [MulSemiringAction G B] (σ : G) (π : B) (n : ℕ) :
    π * σ • uniformizerProd σ π n = uniformizerProd σ π n * σ ^ n • π := by
  rw [smul_uniformizerProd, ← uniformizerProd_succ, uniformizerProd,
    Finset.prod_range_succ' (fun i => σ ^ i • π) n]
  simp [mul_comm]

/-- ★★★**抽象核の心臓** `π · (σ•α_n − α_n) = α_n · (σ^n•π − π)`。

★これが Lemma 6.5 の `−1` を**引き算なしで**表す形である(冒頭の erratum 節)。
`σ = 1` でも `n = 0` でも(両辺 `0` になって)成り立つ。 -/
theorem uniformizer_mul_smul_sub {B : Type*} [CommRing B] {G : Type*} [Monoid G]
    [MulSemiringAction G B] (σ : G) (π : B) (n : ℕ) :
    π * (σ • uniformizerProd σ π n - uniformizerProd σ π n)
      = uniformizerProd σ π n * (σ ^ n • π - π) := by
  rw [mul_sub, mul_sub, uniformizer_mul_smul_uniformizerProd,
    mul_comm π (uniformizerProd σ π n)]

/-! ## §2 付値の不変性 -/

/-- **環自己同型は素元を素元に写す**——`𝔪 = (π)` なら `𝔪 = (σ•π)`。

★`span_smul_uniformizer`(`RamificationJumpDivisibility.lean`)は `τ ∈ G_0` を仮定するが、
実は**任意の群元**でよい: `y ∈ 𝔪` に対し `σ⁻¹•y = π c` を `σ` で写すだけ。 -/
theorem maximalIdeal_eq_span_smul {B : Type*} [CommRing B] [IsLocalRing B] {G : Type*} [Group G]
    [MulSemiringAction G B] {π : B} (hπ : maximalIdeal B = Ideal.span {π}) (σ : G) :
    maximalIdeal B = Ideal.span {σ • π} := by
  apply le_antisymm
  · intro y hy
    have h : σ⁻¹ • y ∈ Ideal.span ({π} : Set B) := by
      rw [← hπ]; exact smul_mem_maximalIdeal σ⁻¹ hy
    rw [Ideal.mem_span_singleton] at h ⊢
    obtain ⟨c, hc⟩ := h
    refine ⟨σ • c, ?_⟩
    have := congrArg (fun z : B => σ • z) hc
    simpa [smul_smul, smul_mul'] using this
  · rw [Ideal.span_le, Set.singleton_subset_iff]
    exact smul_mem_maximalIdeal σ (by rw [hπ]; exact Ideal.mem_span_singleton_self π)

/-- **付値は群作用で不変** `v(σ•x) = v(x)`。

`x = u π^k` と分解し、`σ•u` が単元・`σ•π` が素元であることを使う。 -/
theorem addVal_smul {B : Type*} [CommRing B] [IsDomain B] [IsDiscreteValuationRing B]
    {G : Type*} [Group G] [MulSemiringAction G B] (σ : G) (x : B) :
    addVal B (σ • x) = addVal B x := by
  rcases eq_or_ne x 0 with rfl | hx
  · simp
  obtain ⟨π, hπ⟩ := exists_irreducible B
  obtain ⟨k, u, rfl⟩ := eq_unit_mul_pow_irreducible hx hπ
  have hπ' : Irreducible (σ • π) :=
    (irreducible_iff_uniformizer _).mpr
      (maximalIdeal_eq_span_smul ((irreducible_iff_uniformizer π).mp hπ) σ)
  have hu : IsUnit (σ • (u : B)) := by
    simpa using (u.isUnit.map (MulSemiringAction.toRingHom G B σ))
  rw [smul_mul', smul_pow', addVal_mul, addVal_mul, addVal_pow, addVal_pow,
    addVal_uniformizer hπ, addVal_uniformizer hπ',
    addVal_eq_zero_iff.mpr hu, addVal_eq_zero_iff.mpr u.isUnit]

/-- ★**`v(α_n) = n`**——各因子 `σ^i π` の付値が `1` だから。
★退化検査 (D1): これを落とすと §5 の展開は空虚になる。 -/
theorem addVal_uniformizerProd {B : Type*} [CommRing B] [IsDomain B] [IsDiscreteValuationRing B]
    {G : Type*} [Group G] [MulSemiringAction G B] {π : B} (hπ : Irreducible π) (σ : G) (n : ℕ) :
    addVal B (uniformizerProd σ π n) = n := by
  induction n with
  | zero => simp [uniformizerProd_zero]
  | succ n ih =>
      rw [uniformizerProd_succ, addVal_mul, ih, addVal_smul, addVal_uniformizer hπ]
      push_cast
      ring

theorem uniformizerProd_ne_zero {B : Type*} [CommRing B] [IsDomain B] [IsDiscreteValuationRing B]
    {G : Type*} [Group G] [MulSemiringAction G B] {π : B} (hπ : Irreducible π) (σ : G) (n : ℕ) :
    uniformizerProd σ π n ≠ 0 := by
  intro h
  have hv := addVal_uniformizerProd hπ σ n
  rw [h, addVal_zero] at hv
  exact (ENat.top_ne_coe n) hv

/-- `y ∈ 𝔪^n ⟺ n ≤ v(y)`。 -/
theorem mem_maximalIdeal_pow_iff_le_addVal {B : Type*} [CommRing B] [IsDomain B]
    [IsDiscreteValuationRing B] {π : B} (hπ : Irreducible π) (n : ℕ) (x : B) :
    x ∈ (maximalIdeal B) ^ n ↔ (n : ℕ∞) ≤ addVal B x := by
  rw [(irreducible_iff_uniformizer π).1 hπ, Ideal.span_singleton_pow,
    Ideal.mem_span_singleton, ← addVal_le_iff_dvd, addVal_pow, addVal_uniformizer hπ]
  simp

theorem exists_mul_eq_of_addVal_le {B : Type*} [CommRing B] [IsDomain B]
    [IsDiscreteValuationRing B] {a q : B} (h : addVal B a ≤ addVal B q) : ∃ w : B, a * w = q := by
  obtain ⟨c, hc⟩ := addVal_le_iff_dvd.1 h
  exact ⟨c, hc.symm⟩

theorem uniformizerProd_mem_maximalIdeal_pow {B : Type*} [CommRing B] [IsDomain B]
    [IsDiscreteValuationRing B] {G : Type*} [Group G] [MulSemiringAction G B] {π : B}
    (hπ : Irreducible π) (σ : G) (n : ℕ) : uniformizerProd σ π n ∈ (maximalIdeal B) ^ n :=
  (mem_maximalIdeal_pow_iff_le_addVal hπ n _).mpr (le_of_eq (addVal_uniformizerProd hπ σ n).symm)

/-! ## §3 Yoshida Lemma 6.5 の第 1 主張

★掛け算形 `1 + i_{α_n}(σ) = n + i_π(σ^n)`。`ℕ∞` の引き算は現れない。 -/

/-- **Yoshida `Def 6.1` の `i(σ) := v(σα − α)`**(`i(1) = ⊤`)。

★`α` は素元であることを要求していない(素元のときが原典の `i`)。
★この木は Y1 が `addVal` 言語、Y3/Y4 が `ramCoeff` 言語という 2 つの平行な語彙に
なっていた。`ramIndex` は両者の橋の `addVal` 側の名前である(§6)。 -/
noncomputable def ramIndex {B : Type*} [CommRing B] [IsDomain B] [IsDiscreteValuationRing B]
    {G : Type*} [Monoid G] [MulSemiringAction G B] (α : B) (σ : G) : ℕ∞ :=
  addVal B (σ • α - α)

/-- `i(1) = ⊤`——原典が `i(id) = ∞` と置いている点。 -/
theorem ramIndex_one {B : Type*} [CommRing B] [IsDomain B] [IsDiscreteValuationRing B]
    {G : Type*} [Monoid G] [MulSemiringAction G B] (α : B) : ramIndex α (1 : G) = ⊤ := by
  simp [ramIndex]

theorem ramIndex_eq_top_iff {B : Type*} [CommRing B] [IsDomain B] [IsDiscreteValuationRing B]
    {G : Type*} [Monoid G] [MulSemiringAction G B] (α : B) (σ : G) :
    ramIndex α σ = ⊤ ↔ σ • α = α := by
  rw [ramIndex, addVal_eq_top_iff, sub_eq_zero]

/-- ★★★**Yoshida 2008 Lemma 6.5(第 1 主張・掛け算形)**——`α_n := ∏_{i<n} σ^i π` に対し

`1 + v(σ•α_n − α_n) = n + i(σ^n)`。

★原文は `v(σα − α) = n + i(σ^n)` と書くが、`Def 6.1` の下ではそれは偽で、
正しくは右辺が `1` 小さい(冒頭の erratum 節、`n = 1` で決定的)。ここでは
両辺に `1` を足した形で述べており、`ℕ∞` の引き算は出ない。
★`σ ∈ G_1` は仮定していない。`σ = 1` でも `n = 0` でも(両辺 `⊤` / `⊤` で)真。

証明は §1 の環の恒等式に `addVal` を当てるだけ(積は和になる)。 -/
theorem one_add_ramIndex_uniformizerProd {B : Type*} [CommRing B] [IsDomain B]
    [IsDiscreteValuationRing B] {G : Type*} [Group G] [MulSemiringAction G B] {π : B}
    (hπ : Irreducible π) (σ : G) (n : ℕ) :
    1 + ramIndex (uniformizerProd σ π n) σ = n + ramIndex π (σ ^ n) := by
  have h := congrArg (addVal B) (uniformizer_mul_smul_sub σ π n)
  rw [addVal_mul, addVal_mul, addVal_uniformizer hπ, addVal_uniformizerProd hπ] at h
  exact h

/-- **Lemma 6.5 の第 1 主張(存在形)**——各 `n : ℕ` に対し `v(α) = n` かつ
`1 + v(σα − α) = n + i(σ^n)` を満たす `α ≠ 0` が取れる。

★原文は `n ∈ Z` だが、本ファイルの舞台は付値環 `B` なので `n : ℕ` に制限した(逸脱 2)。 -/
theorem exists_addVal_eq_and_one_add_ramIndex {B : Type*} [CommRing B] [IsDomain B]
    [IsDiscreteValuationRing B] {G : Type*} [Group G] [MulSemiringAction G B] {π : B}
    (hπ : Irreducible π) (σ : G) (n : ℕ) :
    ∃ α : B, α ≠ 0 ∧ addVal B α = n ∧ 1 + ramIndex α σ = n + ramIndex π (σ ^ n) :=
  ⟨uniformizerProd σ π n, uniformizerProd_ne_zero hπ σ n, addVal_uniformizerProd hπ σ n,
    one_add_ramIndex_uniformizerProd hπ σ n⟩

/-! ## §4 抽象核 —— 加法付値の有限和

★mathlib には乗法版 `Valuation.map_sum_eq_of_lt` しかない(`AddValuation` は手書きの
ミラーなので真の欠落)。ここでは一般の `LinearOrderedAddCommMonoidWithTop` で書く。
★`⊤` の扱いを「`c < v 0`」という形の仮定に押し込めることで、`Γ` に `OrderTop` の
補題を要求せずに済ませている。 -/

/-- 各項の付値が `c` より真に大きければ、和の付値も `c` より真に大きい。
★空和 `v 0 = ⊤` を許すため仮定 `c < v 0`(= `c ≠ ⊤`)が要る。 -/
theorem lt_map_sum {R Γ ι : Type*} [Ring R] [LinearOrderedAddCommMonoidWithTop Γ]
    (v : AddValuation R Γ) {c : Γ} (hc : c < v (0 : R)) {s : Finset ι} {f : ι → R}
    (h : ∀ i ∈ s, c < v (f i)) : c < v (∑ i ∈ s, f i) := by
  classical
  induction s using Finset.induction with
  | empty => simpa using hc
  | insert a s ha ih =>
      rw [Finset.sum_insert ha]
      refine lt_of_lt_of_le ?_ (v.map_add _ _)
      exact lt_min (h a (Finset.mem_insert_self a s))
        (ih fun i hi => h i (Finset.mem_insert_of_mem hi))

/-- ★★**付値が一意に最小の項があれば、有限和の付値はその項の付値**。

★`v (f j) = ⊤` の退化(= `f j = 0`)も込みで真: そのとき仮定から `s = {j}` になる。 -/
theorem map_sum_eq_of_lt {R Γ ι : Type*} [Ring R] [LinearOrderedAddCommMonoidWithTop Γ]
    (v : AddValuation R Γ) {s : Finset ι} {f : ι → R} {j : ι} (hj : j ∈ s)
    (h : ∀ i ∈ s, i ≠ j → v (f j) < v (f i)) :
    v (∑ i ∈ s, f i) = v (f j) := by
  classical
  rcases lt_or_ge (v (f j)) (v (0 : R)) with hc | hc
  · rw [← Finset.add_sum_erase s f hj]
    have hlt : v (f j) < v (∑ i ∈ s.erase j, f i) :=
      lt_map_sum v hc fun i hi =>
        h i (Finset.mem_of_mem_erase hi) (Finset.ne_of_mem_erase hi)
    rw [v.map_add_of_distinct_val (ne_of_lt hlt), min_eq_left hlt.le]
  · have hst : v (f j) = ⊤ := by
      rw [v.map_zero] at hc
      exact le_antisymm le_top hc
    have hs : s = {j} := by
      refine Finset.eq_singleton_iff_unique_mem.mpr ⟨hj, fun i hi => ?_⟩
      by_contra hij
      exact absurd (hst ▸ h i hi hij) (not_lt_of_ge le_top)
    rw [hs, Finset.sum_singleton]

/-! ## §5 打ち切り π 進展開(Lemma 6.5 の第 2 主張)

原文の `x = Σ^∞_{n=v(x)} x_n` を **`𝔪^M` を法とする有限和**に切り詰めた形。
無限和・収束・完備性(Appendix I)も一意性も入れていない(逸脱 4・5)。 -/

/-- 局所環では極大イデアルの外は単元。 -/
theorem isUnit_of_notMem_maximalIdeal {B : Type*} [CommRing B] [IsLocalRing B] {x : B}
    (h : x ∉ maximalIdeal B) : IsUnit x := by
  rw [IsLocalRing.mem_maximalIdeal, mem_nonunits_iff, not_not] at h
  exact h

/-- ★**桁の付値は `n` ちょうど** `v(u · α_n) = n`(`u` 単元)。
★退化検査 (D2): これが消費側の「相異なる付値」を支える。 -/
theorem addVal_mul_uniformizerProd {B : Type*} [CommRing B] [IsDomain B]
    [IsDiscreteValuationRing B] {G : Type*} [Group G] [MulSemiringAction G B] {π : B}
    (hπ : Irreducible π) (σ : G) (n : ℕ) {u : B} (hu : IsUnit u) :
    addVal B (u * uniformizerProd σ π n) = n := by
  rw [addVal_mul, addVal_eq_zero_iff.mpr hu, zero_add, addVal_uniformizerProd hπ]

/-- ★★★**Yoshida 2008 Lemma 6.5(第 2 主張・打ち切り形)**——`x ∈ 𝔪^N` かつ `N ≤ M` なら
桁 `a : ℕ → A` が取れて

`x − Σ_{n ∈ [N,M)} a_n · α_n ∈ 𝔪^M`、かつ 各 `a_n` は **`0` か(`B` で)単元**。

`M` についての帰納。段では `w := x − (これまでの和)` に対し `v(α_M) = M ≤ v(w)` から
`w = α_M · b` と割り、`hres` で `b` の代表 `a_M ∈ A` を取る(`b ∈ 𝔪` なら `a_M := 0`)。

★仮定 `hres` は「剰余体が伸びない」= `A → B/𝔪_B` が全射。完全分岐から出る(§7)。
★**単元性は `B` 側**で言っている(`IsUnit (algebraMap A B (a n))`)ので
`IsLocalHom (algebraMap A B)` を要求しない。 -/
theorem exists_digits {A B : Type*} [CommRing A] [CommRing B] [Algebra A B] [IsDomain B]
    [IsDiscreteValuationRing B] {G : Type*} [Group G] [MulSemiringAction G B] {π : B}
    (hπ : Irreducible π) (σ : G)
    (hres : ∀ b : B, ∃ a : A, b - algebraMap A B a ∈ maximalIdeal B) (N : ℕ) :
    ∀ (M : ℕ), N ≤ M → ∀ x : B, x ∈ (maximalIdeal B) ^ N →
      ∃ a : ℕ → A, (∀ n, a n = 0 ∨ IsUnit (algebraMap A B (a n))) ∧
        x - ∑ n ∈ Finset.Ico N M, algebraMap A B (a n) * uniformizerProd σ π n
          ∈ (maximalIdeal B) ^ M := by
  classical
  intro M
  induction M with
  | zero =>
      intro _ x _
      exact ⟨0, fun n => Or.inl rfl, by rw [pow_zero, Ideal.one_eq_top]; exact Submodule.mem_top⟩
  | succ M ih =>
      intro hNM x hx
      rcases Nat.lt_or_ge M N with hMN | hMN
      · -- `N = M + 1`: 和は空で、結論は仮定そのもの
        have hN : N = M + 1 := le_antisymm hNM hMN
        subst hN
        exact ⟨0, fun n => Or.inl rfl, by simpa using hx⟩
      obtain ⟨a, ha, hrem⟩ := ih hMN x hx
      set S := ∑ n ∈ Finset.Ico N M, algebraMap A B (a n) * uniformizerProd σ π n with hS
      have hval : (M : ℕ∞) ≤ addVal B (x - S) :=
        (mem_maximalIdeal_pow_iff_le_addVal hπ M _).mp hrem
      obtain ⟨b, hb⟩ : ∃ b : B, uniformizerProd σ π M * b = x - S :=
        exists_mul_eq_of_addVal_le (by rw [addVal_uniformizerProd hπ]; exact hval)
      obtain ⟨c, hc⟩ := hres b
      have key : ∀ d : A, (d = 0 ∨ IsUnit (algebraMap A B d)) →
          b - algebraMap A B d ∈ maximalIdeal B →
          ∃ a' : ℕ → A, (∀ n, a' n = 0 ∨ IsUnit (algebraMap A B (a' n))) ∧
            x - ∑ n ∈ Finset.Ico N (M + 1), algebraMap A B (a' n) * uniformizerProd σ π n
              ∈ (maximalIdeal B) ^ (M + 1) := by
        intro d hdu hd
        have hupd : ∀ n ∈ Finset.Ico N M,
            algebraMap A B (Function.update a M d n) * uniformizerProd σ π n
              = algebraMap A B (a n) * uniformizerProd σ π n := fun n hn => by
          rw [Function.update_of_ne (Nat.ne_of_lt (Finset.mem_Ico.mp hn).2)]
        refine ⟨Function.update a M d, fun n => ?_, ?_⟩
        · rcases eq_or_ne n M with rfl | hn
          · rw [Function.update_self]; exact hdu
          · rw [Function.update_of_ne hn]; exact ha n
        · rw [Finset.sum_Ico_succ_top hMN, Finset.sum_congr rfl hupd, Function.update_self, ← hS]
          have hid : x - (S + algebraMap A B d * uniformizerProd σ π M)
              = uniformizerProd σ π M * (b - algebraMap A B d) := by linear_combination -hb
          rw [hid, pow_succ]
          exact Ideal.mul_mem_mul (uniformizerProd_mem_maximalIdeal_pow hπ σ M) hd
      by_cases hbm : b ∈ maximalIdeal B
      · exact key 0 (Or.inl rfl) (by simpa using hbm)
      · refine key c (Or.inr (isUnit_of_notMem_maximalIdeal fun hcm => hbm ?_)) hc
        have hbe : b = (b - algebraMap A B c) + algebraMap A B c := by ring
        rw [hbe]
        exact Ideal.add_mem _ hc hcm

/-- ★**桁は `σ` で動かない**——代表系を `A`(具体層では `𝒪_K`)から取ると
`σ • algebraMap A B a = algebraMap A B a`(`smul_algebraMap`)なので、
有限和への `σ` の作用は `α_n` への作用に落ちる。

★これが「原典が Appendix I の収束と項別作用で言っていること」の、有限和版である
(連続性は要らない。逸脱 4)。 -/
theorem smul_digitSum {A B : Type*} [CommRing A] [CommRing B] [Algebra A B] {G : Type*}
    [Monoid G] [MulSemiringAction G B] [SMulCommClass G A B] (σ : G) (π : B) (a : ℕ → A)
    (s : Finset ℕ) :
    σ • ∑ n ∈ s, algebraMap A B (a n) * uniformizerProd σ π n
      = ∑ n ∈ s, algebraMap A B (a n) * σ • uniformizerProd σ π n := by
  rw [show (σ • ∑ n ∈ s, algebraMap A B (a n) * uniformizerProd σ π n)
      = MulSemiringAction.toRingHom G B σ
          (∑ n ∈ s, algebraMap A B (a n) * uniformizerProd σ π n) from rfl,
    map_sum]
  exact Finset.sum_congr rfl fun n _ => by
    show σ • (algebraMap A B (a n) * uniformizerProd σ π n) = _
    rw [smul_mul', smul_algebraMap]

/-! ## §6 `i` と `G_n` の橋

★Y1 は `addVal` 言語、Y3/Y4 は `ramCoeff` 言語で書かれていて、`addVal` は Y3・Y4 に
1 度も現れない。ここで 2 つを型で往復できるようにしておく。 -/

/-- **原典 `Def 6.1` そのもの** `σ ∈ G_n ⟺ i(σ) > n`。 -/
theorem mem_lowerRamificationGroup_iff_lt_ramIndex {A B : Type*} [CommRing A] [CommRing B]
    [Algebra A B] [IsDomain B] [IsDiscreteValuationRing B] {G : Type*} [Group G]
    [MulSemiringAction G B] [SMulCommClass G A B] {α : B}
    (huni : maximalIdeal B = Ideal.span {α}) (hadj : Algebra.adjoin A ({α} : Set B) = ⊤)
    (n : ℕ) (σ : G) :
    σ ∈ lowerRamificationGroup B G n ↔ (n : ℕ∞) < ramIndex α σ :=
  mem_lowerRamificationGroup_iff_lt_addVal huni hadj n σ

/-- ★★**跳びの位置と `i` の橋** `i(σ) = n+1 ⟺ (σ ∈ G_n ∧ σ ∉ G_{n+1})`。

★Y3/Y4 は `ramCoeff` 側で `σ ∈ G_n ∖ G_{n+1}` を扱っている
(`residue_ramCoeff_eq_zero_iff`)。この補題でそちらから `i` の値が読める。 -/
theorem ramIndex_eq_iff_mem_sdiff {A B : Type*} [CommRing A] [CommRing B]
    [Algebra A B] [IsDomain B] [IsDiscreteValuationRing B] {G : Type*} [Group G]
    [MulSemiringAction G B] [SMulCommClass G A B] {α : B}
    (huni : maximalIdeal B = Ideal.span {α}) (hadj : Algebra.adjoin A ({α} : Set B) = ⊤)
    (n : ℕ) (σ : G) :
    ramIndex α σ = (n : ℕ∞) + 1 ↔
      (σ ∈ lowerRamificationGroup B G n ∧ σ ∉ lowerRamificationGroup B G (n + 1)) := by
  rw [mem_lowerRamificationGroup_iff_lt_ramIndex (A := A) huni hadj,
    mem_lowerRamificationGroup_iff_lt_ramIndex (A := A) huni hadj, not_lt,
    ← ENat.add_one_le_iff (by simp), Nat.cast_add, Nat.cast_one]
  exact ⟨fun h => ⟨h.ge, h.le⟩, fun h => le_antisymm h.2 h.1⟩

/-! ## §7 `PAdicLocalField` への具体化

`B := adjoinIntegers K x`(`= 𝒪_{K(x)}`)、`A := 𝒪[K.carrier]`、`G := Gal(K(x)/K)`。
★`IsDiscreteValuationRing.addVal` を statement に出すので DVR を
`attribute [local instance]` で入れる(`lean-idioms.md` #85)。 -/

variable {p : ℕ} [Fact p.Prime]

attribute [local instance] isDiscreteValuationRing_adjoinIntegers

/-- ★★★**Yoshida 2008 Lemma 6.5 第 1 主張(`PAdicLocalField` 版)**。 -/
theorem exists_addVal_eq_and_one_add_ramIndex_adjoin (K : PAdicLocalField p) (x : K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    {π : adjoinIntegers K x} (hπ : Irreducible π)
    (σ : (IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure))) (n : ℕ) :
    ∃ α : adjoinIntegers K x, α ≠ 0 ∧ addVal (adjoinIntegers K x) α = n ∧
      1 + ramIndex α σ = n + ramIndex π (σ ^ n) :=
  exists_addVal_eq_and_one_add_ramIndex hπ σ n

/-- ★★★**Yoshida 2008 Lemma 6.5 第 2 主張(`PAdicLocalField` 版・打ち切り π 進展開)**。

剰余系は完全分岐から出る
(`exists_sub_algebraMap_mem_maximalIdeal_of_isTotallyRamifiedAdjoin`)。
★桁は `𝒪_K` の元なので `σ` で動かない(`smul_digitSum`)。 -/
theorem exists_digits_adjoin (K : PAdicLocalField p) (x : K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (ht : IsTotallyRamifiedAdjoin K x) {π : adjoinIntegers K x} (hπ : Irreducible π)
    (σ : (IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure)))
    {N M : ℕ} (hNM : N ≤ M) {z : adjoinIntegers K x}
    (hz : z ∈ (maximalIdeal (adjoinIntegers K x)) ^ N) :
    ∃ a : ℕ → 𝒪[K.carrier],
      (∀ n, a n = 0 ∨ IsUnit (algebraMap 𝒪[K.carrier] (adjoinIntegers K x) (a n))) ∧
        z - ∑ n ∈ Finset.Ico N M,
              algebraMap 𝒪[K.carrier] (adjoinIntegers K x) (a n) * uniformizerProd σ π n
          ∈ (maximalIdeal (adjoinIntegers K x)) ^ M :=
  exists_digits hπ σ
    (exists_sub_algebraMap_mem_maximalIdeal_of_isTotallyRamifiedAdjoin K x ht) N M hNM z hz

end ABC3.Found.PGC
