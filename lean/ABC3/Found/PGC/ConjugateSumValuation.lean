import ABC3.Found.PGC.RamificationJumpDivisibility

/-!
# 共役和の付値が上がる —— `σ ∈ G_1` ⟹ `v(Σ_{i<p} σ^i α) > v(α)`(Yoshida 2008 Lemma 6.4)

典拠: T. Yoshida, *Local Class Field Theory via Lubin-Tate Theory*
(Ann. Fac. Sci. Toulouse 17-2, 2008; arXiv math/0606108) **Lemma 6.4**(物理 p.14)。
構造化済み原文は `ResearchPaper/1_Structured/Local Class Field Theory via Lubin-Tate Theory/`
`section-6.html` の `id="lemma-6-4"`(`data-pdf-page="14"`)。原文は

> Lemma 6.4. For σ ∈ G_1, we have v(Σ_{i=0}^{p−1} σ^i(α)) > v(α) for all α ∈ K'^×.

★★**原文を `.txt` で読んではならない**: `pdftotext` は大きな括弧と総和記号 `Σ` を
**一切出力しない**(構造化 HTML は `data-txt=""` でこの脱落を明示している)ので、
`.txt` だけでは `v(Σ σ^i(α))` が `v σ i(α)` のように潰れて読めない。
同じ論文の `Cor 6.3` では `≠` が、`Prop 6.6` では `⟨ ⟩` が同様に落ちる。

## 何を証明したか

**(A) 本ノードの心臓 —— `G_n` の元は「どの元でも」付値を `n` だけ上げる**

* `dvd_smul_sub_self`: `σ ∈ G_n`、`𝔪_B = (π)` ⟹ `α · π^n ∣ σ•α − α`(すべての `α : B`)。
* `addVal_add_le_addVal_smul_sub_self`: その付値版 `v(α) + n ≤ v(σ•α − α)`。

★これが実質の中身である。`G_n` の**定義**が与えるのは `σ•α − α ∈ 𝔪^{n+1}`、
すなわち `v(σ•α − α) ≥ n+1` という**絶対的な**評価にすぎない。`α` の付値が大きいときに
必要なのは `v(α) + n` という**相対的な**評価であり、両者は `v(α) = 0`(単数)のときだけ
一致する。橋渡しは「`α = u·π^k` と分解して、単数部分と素元冪部分に分けて積の法則を回す」
という一段で、`§1` がそれをやっている。

**(B) Lemma 6.4 そのもの**

* `dvd_sum_smul_pow`: `σ ∈ G_1`、`(m : B) ∈ 𝔪_B` ⟹ `α · π ∣ Σ_{i<m} σ^i • α`。
* `sum_smul_pow_mem_span_mul_maximalIdeal`: 同じことをイデアル `(α)·𝔪` の形で。
* `addVal_lt_addVal_sum_smul_pow` / `_char` / `_char'`: 付値の形
  `v(α) < v(Σ_{i<p} σ^i • α)`(`α ≠ 0`)。
* `addVal_lt_addVal_sum_smul_pow_adjoin`: `PAdicLocalField` への具体化。

## ★★原文の段取りとの差分(`(σ−1)^{p−1}` と `𝔽_p[X]` は使わなかった)

原典の証明は 2 行で、構造化 HTML の「読み」欄によれば

> 証明は `(σ−1)` の反復と `Σ_i X^i = (X−1)^{p−1}` in `𝔽_p[X]` の 2 行。

である。実際に形式化してみると、**`𝔽_p[X]` の恒等式は要らなかった**。理由は次のとおり:

* 原典の道は `Σ_{i<p} σ^i ≡ (σ−1)^{p−1} (mod p)` と書き換え、`(σ−1)` を `p−1` 回
  適用して `v` を `p−1` 上げ、法 `p` の誤差項を `v(p) ≥ 1` で吸収する。得られる評価は
  `v(α) + min(p−1, v(p))` で、必要な `v(α) + 1` より**強い**。
* こちらは代わりに **`Σ_{i<p} σ^i α = p·α + Σ_{i<p} (σ^i α − α)`** と分けた。
  `G_1` は部分群なので `σ^i ∈ G_1`、したがって (A) から各項が `v(α) + 1` 以上。
  `p·α` の側は `p ∈ 𝔪_B` から `v(α) + 1` 以上。よって和も `v(α) + 1` 以上。
  ★`p−1` 回の反復も、二項係数も、`𝔽_p[X]` も現れない。

すなわち **(A) を作れば (B) は 10 行**であり、原典が `(σ−1)^{p−1}` を持ち出すのは
(A) を暗黙に使わずに済ませるためだと読める。★ただし (A) 自身が「`σ` は `α` の付値に
比例して差を持ち上げる」という非自明な事実なので、**中身が消えたわけではなく移動した**。

なお、原典の道が与える強い評価 `v(α) + min(p−1, v(p))` は本ファイルでは出していない
(Lemma 6.4 の主張に不要で、下流の `Prop 6.6` も `> v(α)` しか使わない)。

## ★逸脱の記録

1. **`α` を `𝒪_{K'} ∖ {0}` に制限した**(原文は `α ∈ K'^×`)。
   本ファイルの舞台は DVR `B` であって商体 `K'` ではない。`K'` の側の主張は
   `B` の側から**割り算を払うだけ**で出る: `γ ∈ K'^×` に対し `G` で固定される
   `c ∈ B ∖ {0}`(たとえば `c ∈ 𝒪_K`)を `cγ ∈ B` となるように取れば
   `Σ σ^i(cγ) = c · Σ σ^i(γ)` かつ `v(cγ) = v(c) + v(γ)` なので、両辺から `v(c)` が消える。
   ★この木には `K'` 上の付値(`addVal` と `‖·‖` の橋)がまだ無いので、
   その一段は書いていない。**落としたのは配線であって数学ではない。**
2. **完全分岐(`G_0 = ⊤`)を仮定していない**。原典は §6.1 冒頭で `K'/K` を完全分岐と
   置いているが、本ファイルの評価は `G_1 = (𝔪^2).inertia G` という定義だけを使うので
   完全分岐は不要である(仮定を落とした = 主張を強めた方向の逸脱)。
   同じ理由で **Lemma 5.11(`𝒪_{K'} = 𝒪_K[π]`)も使っていない**。
3. **`p` の素数性を仮定していない**。抽象核が要求するのは `(p : B) ∈ 𝔪_B` だけで、
   これは `CharP (ResidueField B) p` から出る(`natCast_mem_maximalIdeal`)。
   ★`p` は**剰余体の標数**であって `B` の標数ではない(`B` は標数 `0` でよい)。
4. **`m` を一般の自然数のままにした版**(`dvd_sum_smul_pow` /
   `addVal_lt_addVal_sum_smul_pow`)を先に置き、`m := p` を `CharP` から供給する形に
   分けた。原典は最初から `p` と書いている。

## ★退化の自己検査

* **`σ ∈ G_1` を落とすと偽**。`e ∣ q−1` の順分岐拡大 `K' = K(π)` で
  `σ(π) = ζπ`(`ζ` は 1 の `e` 乗根、`ζ ≠ 1`、`ζ^p ≠ 1`)を取ると
  `Σ_{i<p} σ^i(π) = π·(ζ^p − 1)/(ζ − 1)` で、`(ζ^p−1)/(ζ−1)` は単数だから
  `v` は上がらない。★`σ ∈ G_0 ∖ G_1` では実際に等号が起きる。
* **`α ≠ 0` を落とすと偽**。`α = 0` なら両辺とも `v(0) = ⊤` で、`⊤ < ⊤` は偽。
  ★原典が `α ∈ K'^×` と書いているのはこの一点のためであり、`α = 0` を
  「込みで真になるように」書き直す余地は**無い**(`≤` に弱めるしかない)。
  よって本ファイルは原典どおり `α ≠ 0` を仮定する。
* ★★**結論を `≥` に弱めると自明**。`v(Σ σ^i α) ≥ min_i v(σ^i α) = v(α)` は
  「`σ` が付値を保つ」ことから無料で出る。**真に効いているのは `>` の一段**、
  すなわち (A) の `+n` の部分である。
* **`(p : B) ∈ 𝔪_B` を落とすと偽**。`p` を剰余体の標数以外の数にすると
  `p·α` の項が `v(α)` ちょうどになり、和の付値が上がらない
  (`σ = 1` に近い場合を考えればよい)。
* **`p = 0` は空虚**。`Finset.range 0` の和は `0` で `v(0) = ⊤`。局所体の剰余体は
  有限体なので `p = 0` は起きないが、抽象核は `p = 0` でも矛盾なく真である。

## 実測

前 6 ノードと同じく**抽象核と具体層を分ける**方針が効いた(7 回連続)。
抽象核 `§1`–`§2` の `lean_check` は 0.17–4.1 秒、`PAdicLocalField` 具体層 `§3` は
1.2–3.7 秒。★見立て 150–300 行に対して**証明本体は約 70 行**。
安く上がった理由は上の「原文の段取りとの差分」がすべてで、
`𝔽_p[X]` の恒等式(`Polynomial.map` の核、二項係数 `C(p−1,j) ≡ (−1)^j`)を
一切通らずに済んだこと。
-/

namespace ABC3.Found.PGC

open IsLocalRing ABC3.Skeleton.PGC IsDiscreteValuationRing
open scoped NNReal Valued

def addVal_lt_addVal_sum_smul_pow_char.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 14, item := "Lemma 6.4", sectionId := "lemma-6-4" }

/-! ## §1 抽象核 —— `G_n` の元は付値を `n` だけ「相対的に」上げる

分岐群の定義が直接与えるのは `σ•α − α ∈ 𝔪^{n+1}`(絶対評価)。ここで作るのは
`α·π^n ∣ σ•α − α`(相対評価)であり、両者の差が本ノードの中身である。 -/

/-- `σ ∈ G_n` かつ `𝔪_B = (π)` なら、すべての `x` で `π^{n+1} ∣ σ•x − x`。
`G_n` の定義(`𝔪^{n+1}` の惰性群)を単項生成に翻訳しただけ。 -/
theorem dvd_smul_sub_self_of_mem_lowerRamificationGroup {B : Type*} [CommRing B] [IsLocalRing B]
    {G : Type*} [Group G] [MulSemiringAction G B] {π : B}
    (hπ : maximalIdeal B = Ideal.span {π}) {n : ℕ} {σ : G}
    (hσ : σ ∈ lowerRamificationGroup B G n) (x : B) :
    π ^ (n + 1) ∣ σ • x - x := by
  have h := mem_lowerRamificationGroup_iff_forall.mp hσ x
  rwa [hπ, Ideal.span_singleton_pow, Ideal.mem_span_singleton] at h

/-- 環自己同型は `𝔪` を保つ(`smul_mem_maximalIdeal`)ので `π ∣ σ•π`。 -/
theorem dvd_smul_uniformizer {B : Type*} [CommRing B] [IsLocalRing B]
    {G : Type*} [Group G] [MulSemiringAction G B] {π : B}
    (hπ : maximalIdeal B = Ideal.span {π}) (σ : G) : π ∣ σ • π := by
  have hmem : π ∈ maximalIdeal B := by rw [hπ]; exact Ideal.mem_span_singleton_self π
  have := smul_mem_maximalIdeal σ hmem
  rwa [hπ, Ideal.mem_span_singleton] at this

/-- **素元の冪では差が `k + n` まで上がる**——`σ ∈ G_n` なら `π^{k+n} ∣ σ•(π^k) − π^k`。

`a := σ•π`、`b := π` として `a^{k+1} − b^{k+1} = (a^k − b^k)·a + b^k·(a − b)` で分ける。
第 1 項は帰納法の仮定 `π^{k+n} ∣ a^k − b^k` と `π ∣ a`、
第 2 項は `π^{n+1} ∣ a − b` から。 -/
theorem dvd_smul_pow_sub_pow {B : Type*} [CommRing B] [IsLocalRing B]
    {G : Type*} [Group G] [MulSemiringAction G B] {π : B}
    (hπ : maximalIdeal B = Ideal.span {π}) {n : ℕ} {σ : G}
    (hσ : σ ∈ lowerRamificationGroup B G n) (k : ℕ) :
    π ^ (k + n) ∣ σ • π ^ k - π ^ k := by
  induction k with
  | zero => simp
  | succ k ih =>
      have hpow : σ • π ^ (k + 1) = (σ • π) ^ (k + 1) := smul_pow' σ π (k + 1)
      have hpowk : σ • π ^ k = (σ • π) ^ k := smul_pow' σ π k
      rw [hpowk] at ih
      have hsplit : (σ • π) ^ (k + 1) - π ^ (k + 1)
          = ((σ • π) ^ k - π ^ k) * (σ • π) + π ^ k * (σ • π - π) := by ring
      rw [hpow, hsplit]
      refine dvd_add ?_ ?_
      · have : π ^ (k + n) * π ∣ ((σ • π) ^ k - π ^ k) * (σ • π) :=
          mul_dvd_mul ih (dvd_smul_uniformizer hπ σ)
        rwa [show π ^ (k + n) * π = π ^ (k + 1 + n) by ring] at this
      · have : π ^ k * π ^ (n + 1) ∣ π ^ k * (σ • π - π) :=
          mul_dvd_mul_left _ (dvd_smul_sub_self_of_mem_lowerRamificationGroup hπ hσ π)
        rwa [show π ^ k * π ^ (n + 1) = π ^ (k + 1 + n) by ring] at this

/-- ★★**本ノードの心臓**——`σ ∈ G_n` なら `σ` は**どの元でも**付値を `n` だけ上げる:
`α · π^n ∣ σ•α − α`。

`α = u·π^k`(`u` 単数)と分解し
`σ•(u π^k) − u π^k = (σ•u − u)·σ•(π^k) + u·(σ•(π^k) − π^k)` で分ける。
第 1 項は `π^{n+1} ∣ σ•u − u` と `π^k ∣ σ•(π^k)`、第 2 項は `dvd_smul_pow_sub_pow`。
★`α = 0` でも真(両辺 `0`)。★`u` が単数であることは `Units.mul_left_dvd` で消える。 -/
theorem dvd_smul_sub_self {B : Type*} [CommRing B] [IsDomain B]
    [IsDiscreteValuationRing B] {G : Type*} [Group G] [MulSemiringAction G B] {π : B}
    (hπ : maximalIdeal B = Ideal.span {π}) {n : ℕ} {σ : G}
    (hσ : σ ∈ lowerRamificationGroup B G n) (α : B) :
    α * π ^ n ∣ σ • α - α := by
  rcases eq_or_ne α 0 with rfl | hα
  · simp
  have hirr : Irreducible π := (irreducible_iff_uniformizer π).mpr hπ
  obtain ⟨k, u, rfl⟩ := eq_unit_mul_pow_irreducible hα hirr
  have hdvdk : π ^ k ∣ σ • π ^ k := by
    rw [smul_pow']
    exact pow_dvd_pow_of_dvd (dvd_smul_uniformizer hπ σ) k
  have hsplit : σ • ((u : B) * π ^ k) - (u : B) * π ^ k
      = (σ • (u : B) - (u : B)) * (σ • π ^ k) + (u : B) * (σ • π ^ k - π ^ k) := by
    rw [smul_mul']; ring
  rw [show (u : B) * π ^ k * π ^ n = (u : B) * π ^ (k + n) by ring, Units.mul_left_dvd, hsplit]
  refine dvd_add ?_ ?_
  · have h1 : π ^ (n + 1) * π ^ k ∣ (σ • (u : B) - (u : B)) * (σ • π ^ k) :=
      mul_dvd_mul (dvd_smul_sub_self_of_mem_lowerRamificationGroup hπ hσ _) hdvdk
    refine dvd_trans ?_ h1
    rw [← pow_add]
    exact pow_dvd_pow π (by omega)
  · exact Dvd.dvd.mul_left (dvd_smul_pow_sub_pow hπ hσ k) _

/-- `dvd_smul_sub_self` の付値版——`σ ∈ G_n` なら `v(α) + n ≤ v(σ•α − α)`。
★`n = 0` では「`σ` が付値を下げない」という自明な主張になる。 -/
theorem addVal_add_le_addVal_smul_sub_self {B : Type*} [CommRing B] [IsDomain B]
    [IsDiscreteValuationRing B] {G : Type*} [Group G] [MulSemiringAction G B] {π : B}
    (hπ : maximalIdeal B = Ideal.span {π}) {n : ℕ} {σ : G}
    (hσ : σ ∈ lowerRamificationGroup B G n) (α : B) :
    addVal B α + n ≤ addVal B (σ • α - α) := by
  have hirr : Irreducible π := (irreducible_iff_uniformizer π).mpr hπ
  have h := addVal_le_iff_dvd.mpr (dvd_smul_sub_self hπ hσ α)
  rwa [addVal_mul, addVal_pow, addVal_uniformizer hirr, nsmul_eq_mul, mul_one] at h

/-! ## §2 Yoshida Lemma 6.4

`Σ_{i<m} σ^i α = m·α + Σ_{i<m} (σ^i α − α)` と分けるだけ。`G_1` が部分群であることから
`σ^i ∈ G_1` なので §1 が各項に効き、`m·α` の側は `m ∈ 𝔪_B` で処理する。 -/

/-- 剰余体の標数が `p` なら `(p : B) ∈ 𝔪_B`。
★`B` 自身の標数は `0` でよい(混標数の局所環がまさにその場合)。 -/
theorem natCast_mem_maximalIdeal {B : Type*} [CommRing B] [IsLocalRing B] (p : ℕ)
    [CharP (ResidueField B) p] : (p : B) ∈ maximalIdeal B := by
  rw [← Ideal.Quotient.eq_zero_iff_mem]
  show residue B (p : B) = 0
  rw [map_natCast, CharP.cast_eq_zero]

/-- **Lemma 6.4 の割り切りの形**——`σ ∈ G_1` かつ `(m : B) ∈ 𝔪_B` なら
`α·π ∣ Σ_{i<m} σ^i • α`。★`α = 0` でも真。 -/
theorem dvd_sum_smul_pow {B : Type*} [CommRing B] [IsDomain B]
    [IsDiscreteValuationRing B] {G : Type*} [Group G] [MulSemiringAction G B] {π : B}
    (hπ : maximalIdeal B = Ideal.span {π}) {σ : G}
    (hσ : σ ∈ lowerRamificationGroup B G 1) {m : ℕ} (hm : (m : B) ∈ maximalIdeal B) (α : B) :
    α * π ∣ ∑ i ∈ Finset.range m, σ ^ i • α := by
  have hsplit : ∑ i ∈ Finset.range m, σ ^ i • α
      = (m : B) * α + ∑ i ∈ Finset.range m, (σ ^ i • α - α) := by
    rw [Finset.sum_sub_distrib, Finset.sum_const, Finset.card_range, nsmul_eq_mul]
    ring
  rw [hsplit]
  refine dvd_add ?_ (Finset.dvd_sum fun i _ => ?_)
  · have hpi : π ∣ (m : B) := by rwa [hπ, Ideal.mem_span_singleton] at hm
    rw [mul_comm (m : B) α]
    exact mul_dvd_mul_left α hpi
  · have := dvd_smul_sub_self hπ (pow_mem hσ i) α
    rwa [pow_one] at this

/-- 同じことをイデアル `(α)·𝔪_B` の言葉で。 -/
theorem sum_smul_pow_mem_span_mul_maximalIdeal {B : Type*} [CommRing B] [IsDomain B]
    [IsDiscreteValuationRing B] {G : Type*} [Group G] [MulSemiringAction G B] {π : B}
    (hπ : maximalIdeal B = Ideal.span {π}) {σ : G}
    (hσ : σ ∈ lowerRamificationGroup B G 1) {m : ℕ} (hm : (m : B) ∈ maximalIdeal B) (α : B) :
    (∑ i ∈ Finset.range m, σ ^ i • α) ∈ Ideal.span {α} * maximalIdeal B := by
  rw [hπ, Ideal.span_singleton_mul_span_singleton, Ideal.mem_span_singleton]
  exact dvd_sum_smul_pow hπ hσ hm α

/-- **Lemma 6.4 の付値の形(`m` 一般版)**——`v(α) < v(Σ_{i<m} σ^i • α)`。
★`α ≠ 0` は落とせない(`α = 0` なら両辺 `⊤` で偽)。 -/
theorem addVal_lt_addVal_sum_smul_pow {B : Type*} [CommRing B] [IsDomain B]
    [IsDiscreteValuationRing B] {G : Type*} [Group G] [MulSemiringAction G B] {π : B}
    (hπ : maximalIdeal B = Ideal.span {π}) {σ : G}
    (hσ : σ ∈ lowerRamificationGroup B G 1) {m : ℕ} (hm : (m : B) ∈ maximalIdeal B)
    {α : B} (hα : α ≠ 0) :
    addVal B α < addVal B (∑ i ∈ Finset.range m, σ ^ i • α) := by
  have hirr : Irreducible π := (irreducible_iff_uniformizer π).mpr hπ
  have h := addVal_le_iff_dvd.mpr (dvd_sum_smul_pow hπ hσ hm α)
  rw [addVal_mul, addVal_uniformizer hirr] at h
  obtain ⟨k, hk⟩ := ENat.ne_top_iff_exists.mp (fun hc => hα (addVal_eq_top_iff.mp hc))
  rw [← hk] at h ⊢
  refine lt_of_lt_of_le ?_ h
  exact_mod_cast Nat.lt_succ_self k

/-- ★★★**Yoshida 2008 Lemma 6.4(抽象核)**——`σ ∈ G_1` なら
`v(Σ_{i=0}^{p−1} σ^i(α)) > v(α)`(`α ≠ 0`、`p` は**剰余体**の標数)。

★完全分岐も `p` の素数性も要らない(逸脱の記録 2・3)。 -/
theorem addVal_lt_addVal_sum_smul_pow_char {B : Type*} [CommRing B] [IsDomain B]
    [IsDiscreteValuationRing B] {G : Type*} [Group G] [MulSemiringAction G B] (p : ℕ)
    [CharP (ResidueField B) p] {π : B} (hπ : maximalIdeal B = Ideal.span {π}) {σ : G}
    (hσ : σ ∈ lowerRamificationGroup B G 1) {α : B} (hα : α ≠ 0) :
    addVal B α < addVal B (∑ i ∈ Finset.range p, σ ^ i • α) :=
  addVal_lt_addVal_sum_smul_pow hπ hσ (natCast_mem_maximalIdeal p) hα

/-- **素元を明示しない形**——DVR には必ず素元があるので `hπ` は仮定から外せる。 -/
theorem addVal_lt_addVal_sum_smul_pow_char' {B : Type*} [CommRing B] [IsDomain B]
    [IsDiscreteValuationRing B] {G : Type*} [Group G] [MulSemiringAction G B] (p : ℕ)
    [CharP (ResidueField B) p] {σ : G} (hσ : σ ∈ lowerRamificationGroup B G 1)
    {α : B} (hα : α ≠ 0) :
    addVal B α < addVal B (∑ i ∈ Finset.range p, σ ^ i • α) := by
  obtain ⟨π, hπ⟩ := exists_irreducible B
  exact addVal_lt_addVal_sum_smul_pow_char p ((irreducible_iff_uniformizer π).mp hπ) hσ hα

/-! ## §3 `PAdicLocalField` への具体化

`B := adjoinIntegers K x`(`= 𝒪_{K(x)}`)、`G := Gal(K(x)/K)`。
★`IsDiscreteValuationRing.addVal` を statement に出すので DVR を
`attribute [local instance]` で入れる(`lean-idioms.md` #85)。 -/

variable {p : ℕ} [Fact p.Prime]

attribute [local instance] isDiscreteValuationRing_adjoinIntegers

/-- ★★★**Yoshida 2008 Lemma 6.4(`PAdicLocalField` 版)**——`σ ∈ G_1`、`α ≠ 0` なら
`v_{K(x)}(Σ_{i<p} σ^i • α) > v_{K(x)}(α)`。

★原典の `IsTotallyRamifiedAdjoin` は**不要**(逸脱の記録 2)。
★`α` は `𝒪_{K(x)} ∖ {0}` の元であって `K(x)^×` の元ではない(逸脱の記録 1)。 -/
theorem addVal_lt_addVal_sum_smul_pow_adjoin (K : PAdicLocalField p) (x : K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    {σ : (IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure))}
    (hσ : σ ∈ lowerRamificationGroupAdjoin K x 1)
    {α : adjoinIntegers K x} (hα : α ≠ 0) :
    addVal (adjoinIntegers K x) α
      < addVal (adjoinIntegers K x) (∑ i ∈ Finset.range p, σ ^ i • α) := by
  haveI := charP_residueField_adjoinIntegers K x
  exact addVal_lt_addVal_sum_smul_pow_char' p hσ hα

/-- 和を `K(x)` の元として書き下したもの——原典の `Σ_{i<p} σ^i(α)` の字面。
下流(`Prop 6.6 (Sen)`)が体の側で読むための橋。 -/
theorem coe_sum_smul_pow_adjoin (K : PAdicLocalField p) (x : K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (σ : (IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure)))
    (α : adjoinIntegers K x) (m : ℕ) :
    ((∑ i ∈ Finset.range m, σ ^ i • α : adjoinIntegers K x) :
        IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      = ∑ i ∈ Finset.range m,
          (σ ^ i) (α : IntermediateField.adjoin K.carrier ({x} : Set K.closure)) := by
  rw [AddSubmonoidClass.coe_finsetSum]
  exact Finset.sum_congr rfl fun i _ => coe_smul_adjoinIntegers K x (σ ^ i) α

end ABC3.Found.PGC
