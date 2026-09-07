import ABC3.Found.PGC.SenJumpFiltration

/-!
# Sen の定理 (iii) —— 跳び位置の合同 `i_{j−1} ≡ i_j (mod p^j)`(Yoshida 2008 Prop 6.6 (iii))

典拠: T. Yoshida, *Local Class Field Theory via Lubin-Tate Theory*
(Ann. Fac. Sci. Toulouse 17-2, 2008; arXiv math/0606108) **Proposition 6.6 (Sen)**(物理 p.15)。
構造化済み原文は `ResearchPaper/1_Structured/Local Class Field Theory via Lubin-Tate Theory/`
`section-6.html` の `id="prop-6-6"`(`data-pdf-page="15"`, `data-item="Proposition 6.6 (Sen)"`)。

原文 (Yoshida08 p.15):
> Proposition 6.6 (Sen [14]). Let σ ∈ G_1, and |σ| = p^m for m ≥ 1 (by Proposition 6.2). Let
> H_n := G_n ∩ σ for n ≥ 1 and i_j := i(σ^p^j) for j ≥ 0 (and i_j := ∞ for j ≥ m). Then:
> (i) i_j−1 < i_j if j ≤ m. Also, H_n = σ^p^j if and only if i_j−1 ≤ n < i_j.
> (ii) i(σ^a) = i_v_p(a) for a ≥ 1, where v_p := v_Q[bb]_p.
> (iii) i_j−1 ≡ i_j (mod p^j), where ∞ is understood to be congruent to any integer.

原文の (iii) の証明(同 p.15。上付き・下付きを `^` / `_` に開き、生成部分群の角括弧を補った):
> (iii): We can assume i_j < ∞, and use induction on j. The assertion is empty when j = 0.
> Let j = 1, and assume the Inductive Hypothesis (the assertion of (iii) for j − 1).
> We first prove the Claim: the i_{j−1} and n + i(σ^n) for n ∈ Z, v_p(n) < j are all distinct
> from each other. As v_p(n) ≤ j − 1, the Inductive Hypothesis shows i(σ^n) = i_{v_p(n)}
> ≡ i_{j−1} (mod p^{v_p(n)+1}), i.e. v_p(i_{j−1} − i(σ^n)) > v_p(n), hence i_{j−1} ≠ n + i(σ^n).
> Now assume n + i(σ^n) = n′ + i(σ^{n′}). If v_p(n) ≠ v_p(n′), then
> v_p(n − n′) = min{v_p(n), v_p(n′)}, but the Inductive Hypothesis shows
> v_p(i(σ^n) − i(σ^{n′})) > min{v_p(n), v_p(n′)}, which is impossible. Hence v_p(n) = v_p(n′),
> therefore i(σ^n) = i(σ^{n′}) and n = n′. Thus the Claim is proven. Now applying the Inductive
> Hypothesis to σ^p ∈ G_1, we have i_{j−1} ≡ i_j (mod p^{j−1}). Let s := i_{j−1} − i_j and
> assume v_p(s) = j − 1, to see it leads to contradiction. The first part of Lemma 6.5 for σ^p
> shows that there is x ∈ K′^× with v(x) = s and v(σ^p(x) − x) = s + i((σ^p)^s) = s + i_j
> = i_{j−1}. Letting y := Σ_{i=0}^{p−1} σ^i(x), we have v(y) > v(x) = s by Lemma 6.4 and
> v(σ(y) − y) = v(σ^p(x) − x) = i_{j−1}. Now expand y = Σ_{n≥v(y)} y_n as in Lemma 6.5:
> v(σ(y_n) − y_n) = n + i(σ^n) if y_n ≠ 0. Let z := σ(y) − y. Then v(z) = i_{j−1} and
> z = Σ_{n≥v(y)} z_n, where z_n := σ(y_n) − y_n, hence v(z_n) = n + i(σ^n) whenever z_n ≠ 0.
> The Claim shows v(z − Σ_{v_p(n)<j} z_n) ≤ i_{j−1}. If v_p(n) ≥ j and z_n ≠ 0, then
> v(z_n) = n + i(σ^n) ≥ n + i_j ≥ v(y) + i_j > i_{j−1}, hence v(Σ_{v_p(n)≥j} z_n) > i_{j−1},
> a contradiction.

## 本ファイルの範囲

**(iii) を丸ごと**埋めた(段 1〜12 すべて。落とした段は無い)。
(i)(ii) は `Found/PGC/SenJumpFiltration.lean`(Y7a)、Lemma 6.4 は
`Found/PGC/ConjugateSumValuation.lean`、Lemma 6.5 は `Found/PGC/UniformizerExpansion.lean`。

## 何を証明したか

**§1 抽象核 —— `p` 進整除性の初等算術(ℤ だけ。分岐・付値・Galois の語彙が 1 つも出ない)**

* `dvd_sub_of_forall_dvd_succ` : `∀ k < K, p^{k+1} ∣ f(k+1) − f k` ⟹ `k ≤ l ≤ K` で
  `p^{k+1} ∣ f l − f k`。★原文の「合同の鎖」(段 4)そのもの。
* `add_ne_add_of_not_dvd` : `p^{k+1} ∣ i_j − i_k`、`p^{k+1} ∣ E`、`¬ p^{k+1} ∣ n` ⟹
  `n + i_k ≠ i_j + E`。★原文 Claim の前半(段 5)。
* `eq_of_add_eq_add_of_chain` : 鎖の下で `n + f k = n′ + f k′` ⟹ `n = n′`。
  ★原文 Claim の後半(段 6)。`v_p` の三分法を `k < k′ / k = k′ / k > k′` に開いただけ。

**§2 抽象核 —— 付値の有限和と `ℕ∞` の `1 +`**

* `map_sum_ne_of_pairwise_ne` : 一般の環 `R` と一般の `AddValuation v` に対し、
  項の付値が**対どうし相異なり**すべて `c` と異なるなら `v(Σ) ≠ c`。
  ★原文の「The Claim shows `v(z − Σ_{v_p(n)<j} z_n) ≤ i_{j−1}`」を支える 1 本。
  `Finset.exists_min_image` で最小の項を選び Y6 の `map_sum_eq_of_lt` に渡す。
* `enat_one_add_inj` / `enat_one_add_lt_one_add_iff` : `ℕ∞` で `1 + ·` は単射かつ順序反射。
  ★Lemma 6.5 が掛け算形(`1 + i_α(σ) = n + i(σ^n)`)なので、比較は常に `1 +` を付けた
  ままで行う(引き算を出さないため)。
* `smul_sum_range_sub_self` : `σ•(Σ_{i<n} σ^i•x) − Σ_{i<n} σ^i•x = σ^n•x − x`。
  ★一般の加法群 + `DistribMulAction`。原文の「`v(σ(y) − y) = v(σ^p(x) − x)`」の中身。

**§3 Yoshida Prop 6.6 (iii)(抽象層)**

* `dvd_ramIndex_toNat_sub` : ★**強帰納法の外枠**。`σ` について**全称量化**してある
  (帰納法の仮説を `σ^p` に当てるため。原文が暗黙に使っている点)。
* `ramIndex_pow_pow_congr` : `i_j = c`、`i_{j+1} = d`(ともに有限)なら
  `p^{j+1} ∣ (d − c : ℤ)`。★これが (iii) の本体。
* `ramIndex_pow_pow_congr_or_top` : ★原文の規約「`∞` は任意の整数と合同」まで込めた形。

**§4 `PAdicLocalField` 具体層** —— Y7a §7 と同型。

## ★シフト法(原文の `s < 0` を `𝒪_{K′}` の中に収める)

原文は `s := i_{j−1} − i_j < 0` に Lemma 6.5 を当てて `x ∈ K′^×`(`v(x) = s < 0`)を取る。
本ファイルは `s′ := s + p^j·t`(`t := i_j + 1`)にずらし、`s′ ≥ 1` の**自然数**として扱う。

* `v_p(s′) = j − 1` は `p^j ∣ p^j·t` から保たれる。
* `i((σ^p)^{s′}) = i_j` も生成部分群が変わらないので不変((ii) が与える)。
* 関係式は `s′ + i_j = i_{j−1} + p^j·t` の形で書き、**引き算を一度も書かない**。

⇒ `x, y, z, 桁`がすべて `𝒪_{K′}` に収まり、`addVal : 𝒪 → ℕ∞` だけで足りる。
★分数体に `ℤ` 値付値を載せる必要は無い。★分岐指数 `e` も要らない
(必要なのは `s′ ≥ 1` と `p^{j+1} ∣ (s′ − s)` の 2 つだけ)。
★これで原文のもう 1 つの穴((ii) は `a ≥ 1` を要求するのに (iii) は `a < 0` で使う)も塞がる。

## ★逸脱の記録

1. **添字を 1 つずらした**。原文 `i_{j−1} ≡ i_j (mod p^j)` を
   `i_j ≡ i_{j+1} (mod p^{j+1})` と書いた。`ℕ` の切り詰め引き算 `j − 1` を出さないため
   (Y7a の逸脱 2 と同じ方針)。数学的内容は同じ。
2. ★★**`ℕ∞` の引き算を使わず、有限代表 `c, d : ℕ` を経由して `ℤ` で合同を書いた**。
   `ℕ∞` の引き算は切り詰めるので、(i) が `i_j < i_{j+1}` を与える以上
   `(p^{j+1} : ℕ∞) ∣ (i_j − i_{j+1})` は左辺が常に `0` になり **空虚に真**になる
   (下の退化の自己検査 D1)。原文の「`∞` は任意の整数と合同」は
   `ramIndex π (σ^{p^{j+1}}) = (d : ℕ∞)` という**有限性の仮定**として表した。
3. ★**原文 p.15 の `Let j = 1` は `Let j ≥ 1` の誤植と判断した**。`j = 1` では直後の
   `mod p^{j−1}` が `mod p^0` になり内容が消えるので、帰納法が動かない。
   本ファイルの(シフト後の)強帰納法は `j ≥ 0` のすべてで回っており、
   原文の `j = 1` に対応する `j = 0` も本体を通る(そこでは `σ^p` への帰納法仮説が
   `mod p^0` で自明になるだけである)。
4. **`v_p` を関数として使わず、整除性の述語 `p^k ∣ n ∧ ¬ p^{k+1} ∣ n` で回した**。
   `padicValNat p 0 = 0` というジャンク値を避けるため(Y7a の逸脱 4 と同じ)。
5. **Lemma 6.5 の第 1 主張は掛け算形で使った**。原文の `v(σ(α) − α) = n + i(σ^n)` は
   `Def 6.1` の下では偽で、正しくは `1 + v(σ(α) − α) = n + i(σ^n)`
   (Y6 の erratum)。したがって本ファイルの比較はすべて `1 + ·` を付けたまま行う。
6. **原文の無限和 `y = Σ_{n≥v(y)} y_n` を `𝔪^M` を法とする有限和に切り詰めた**
   (Y6 の `exists_digits`)。`M := i_j + E + 1` と取れば、剰余項 `r` は
   `v(σ•r − r) ≥ v(r) + 1 ≥ M + 1 > v(z)` を満たすので、無限和・収束(Appendix I)は要らない。
7. **原文の Claim は `i_{j−1}` と `n + i(σ^n)` の相異性だが、本ファイルは
   `i_j + E` と `n + i(σ^n)` の相異性になっている**(シフトのため)。
   `p^{j+1} ∣ E` なので `p` 進的な内容は変わらない(§1 `add_ne_add_of_not_dvd` の `E`)。
8. **帰納法を `σ` について全称量化した**。原文は「Now applying the Inductive Hypothesis to
   `σ^p ∈ G_1`」と書くだけだが、これは帰納法の主張が `σ` を動かせる形でないと使えない。
   ★仮定を強めたのではなく、原文が暗黙に使っている量化を明示しただけである。
9. **抽象層は「剰余体が伸びない」を仮定 `hres` として置いた**
   (`∀ b : B, ∃ a : A, b − algebraMap A B a ∈ 𝔪_B`)。具体層では完全分岐から出る
   (`exists_sub_algebraMap_mem_maximalIdeal_of_isTotallyRamifiedAdjoin`)。
10. **`m ≥ 1` を仮定していない**。`j + 1 < m` が主張の中に入っており、`m = 0` では
   仮定 `ramIndex π (σ^{p^{j+1}}) = (d : ℕ∞)` が偽になる(`i_j = ⊤`)ので不要。

## ★退化の自己検査

* **(D1) `ℕ∞` の引き算で書くと空虚に真**。(i)(`ramIndex_pow_pow_lt_succ`)より
  `i_j < i_{j+1}` だから `ℕ∞` の切り詰め引き算では `i_j − i_{j+1} = 0` で、
  `(p^{j+1} : ℕ∞) ∣ 0` は `simp` で通る ——★**内容ゼロの定理**。
  本ファイルは `c, d : ℕ` を経由して `ℤ` の `(d : ℤ) − (c : ℤ)` で書くのでこれを避けている。
* **(D2) `hd`(= `i_{j+1}` が有限)は落とせない**。落とすと原文の規約どおり
  「`∞` は任意の整数と合同」で主張が空虚になる。`ramIndex_pow_pow_congr_or_top` が
  その場合分けを明示している。
* **(D3) `hσ : σ ∈ G_1` は落とせない**。Lemma 6.4(`v(Σ_{i<p} σ^i α) > v(α)`)が
  `p` の剰余体での消滅を使う。これが `v(y) > v(x) = s′` を与え、
  「`v_p(n) ≥ j+1` の桁は付値が大きい」の評価を支えている。
* **(D4) `t := i_{j+1} + 1`(すなわち `E > i_{j+1}`)は落とせない**。`E ≤ i_{j+1}` だと
  `s′` が自然数として取れない(原文の `s < 0` に戻ってしまう)。
-/

namespace ABC3.Found.PGC

open IsLocalRing ABC3.Skeleton.PGC IsDiscreteValuationRing
open scoped NNReal Valued

def ramIndex_pow_pow_congr.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 15, item := "Proposition 6.6", sectionId := "prop-6-6" }

/-! ## §1 抽象核 —— `p` 進整除性の初等算術

★ℤ の整除性しか出てこない。分岐・付値・Galois の語彙は 1 つも現れない。 -/

/-- ★**合同の鎖** —— 隣接段の合同 `p^{k+1} ∣ f(k+1) − f k`(`k < K`)から
`k ≤ l ≤ K` の合同 `p^{k+1} ∣ f l − f k` が出る。

原文の「As `v_p(n) ≤ j − 1`, the Inductive Hypothesis shows
`i(σ^n) = i_{v_p(n)} ≡ i_{j−1} (mod p^{v_p(n)+1})`」がこれである。 -/
theorem dvd_sub_of_forall_dvd_succ {p : ℕ} {f : ℕ → ℤ} {K : ℕ}
    (hstep : ∀ k, k < K → (p : ℤ) ^ (k + 1) ∣ f (k + 1) - f k) {k l : ℕ}
    (hkl : k ≤ l) (hl : l ≤ K) : (p : ℤ) ^ (k + 1) ∣ f l - f k := by
  induction l, hkl using Nat.le_induction with
  | base => simp
  | succ l hkl ih =>
      have h1 : (p : ℤ) ^ (k + 1) ∣ f (l + 1) - f l :=
        dvd_trans (pow_dvd_pow _ (by omega)) (hstep l (by omega))
      have h2 := ih (by omega)
      have hsplit : f (l + 1) - f k = (f (l + 1) - f l) + (f l - f k) := by ring
      rw [hsplit]
      exact dvd_add h1 h2

/-- ★**原文 Claim の前半** —— `p^{k+1}` が `i_j − i_k` と `E` を割り、`n` を割らないなら
`n + i_k ≠ i_j + E`。

原文は `E = 0`(すなわち `i_{j−1} ≠ n + i(σ^n)`)の場合だが、本ファイルはシフト
`s′ = s + p^{j+1}·t` を使うので右辺に `E := p^{j+1}·t` が乗る(逸脱 7)。 -/
theorem add_ne_add_of_not_dvd {p : ℕ} {ik ij E n : ℤ} {k : ℕ}
    (h1 : (p : ℤ) ^ (k + 1) ∣ ij - ik) (h2 : (p : ℤ) ^ (k + 1) ∣ E)
    (h3 : ¬ (p : ℤ) ^ (k + 1) ∣ n) : n + ik ≠ ij + E := by
  intro h
  refine h3 ?_
  have hn : n = (ij - ik) + E := by linarith
  rw [hn]
  exact dvd_add h1 h2

/-- ★**原文 Claim の後半** —— 鎖 `hchain` の下で `n + f k = n′ + f k′` なら `n = n′`。

`p^k ∥ n`、`p^{k′} ∥ n′` として `k` と `k′` を三分する。`k < k′` なら
`p^{k+1} ∣ f k′ − f k = n − n′` と `p^{k+1} ∣ n′` から `p^{k+1} ∣ n` で矛盾。
`k = k′` なら `f k = f k′` を消して `n = n′`。★原文の `v_p` の言葉を整除性に開いただけ。 -/
theorem eq_of_add_eq_add_of_chain {p : ℕ} {f : ℕ → ℤ} {J : ℕ}
    (hchain : ∀ a b : ℕ, a ≤ b → b ≤ J → (p : ℤ) ^ (a + 1) ∣ f b - f a)
    {n n' : ℤ} {k k' : ℕ} (hkJ : k ≤ J) (hk'J : k' ≤ J)
    (hdk : (p : ℤ) ^ k ∣ n) (hnk : ¬ (p : ℤ) ^ (k + 1) ∣ n)
    (hdk' : (p : ℤ) ^ k' ∣ n') (hnk' : ¬ (p : ℤ) ^ (k' + 1) ∣ n')
    (h : n + f k = n' + f k') : n = n' := by
  rcases lt_trichotomy k k' with hlt | heq | hgt
  · exfalso
    have hd : (p : ℤ) ^ (k + 1) ∣ n - n' := by
      have he : n - n' = f k' - f k := by linarith
      rw [he]; exact hchain k k' (le_of_lt hlt) hk'J
    have hn' : (p : ℤ) ^ (k + 1) ∣ n' := dvd_trans (pow_dvd_pow _ (by omega)) hdk'
    exact hnk (by simpa using dvd_add hd hn')
  · subst heq; linarith
  · exfalso
    have hd : (p : ℤ) ^ (k' + 1) ∣ n' - n := by
      have he : n' - n = f k - f k' := by linarith
      rw [he]; exact hchain k' k (le_of_lt hgt) hkJ
    have hn : (p : ℤ) ^ (k' + 1) ∣ n := dvd_trans (pow_dvd_pow _ (by omega)) hdk
    exact hnk' (by simpa using dvd_add hd hn)

/-! ## §2 抽象核 —— 付値の有限和と `ℕ∞` の `1 +` -/

/-- ★★**付値が対どうし相異なる有限和は、どの項の付値とも `c` が違えば `c` にならない**。

一般の環 `R` と一般の `AddValuation v : AddValuation R Γ` で成り立つ。
空和は `v 0 = ⊤`(仮定 `h0`)で片付く。空でないときは `Finset.exists_min_image` で
最小の項を取り、Y6 の `map_sum_eq_of_lt` に渡す(相異性から最小は一意)。

★原文の「The Claim shows `v(z − Σ_{v_p(n)<j} z_n) ≤ i_{j−1}`」を支える 1 本。 -/
theorem map_sum_ne_of_pairwise_ne {R Γ ι : Type*} [Ring R]
    [LinearOrderedAddCommMonoidWithTop Γ] (v : AddValuation R Γ) {s : Finset ι} {f : ι → R}
    {c : Γ} (hpair : ∀ i ∈ s, ∀ i' ∈ s, i ≠ i' → v (f i) ≠ v (f i'))
    (hne : ∀ i ∈ s, v (f i) ≠ c) (h0 : v (0 : R) ≠ c) :
    v (∑ i ∈ s, f i) ≠ c := by
  rcases s.eq_empty_or_nonempty with rfl | hs
  · simpa using h0
  obtain ⟨i₀, hi₀, hmin⟩ := s.exists_min_image (fun i => v (f i)) hs
  rw [map_sum_eq_of_lt v hi₀ fun i hi hi' =>
    lt_of_le_of_ne (hmin i hi) (hpair i₀ hi₀ i hi (Ne.symm hi'))]
  exact hne i₀ hi₀

/-- `ℕ∞` で `1 + ·` は単射。★Lemma 6.5 の掛け算形を比較に使うための補題。 -/
theorem enat_one_add_inj {a b : ℕ∞} (h : 1 + a = 1 + b) : a = b :=
  ENat.add_left_injective_of_ne_top ENat.one_ne_top (by simpa [add_comm] using h)

/-- `ℕ∞` で `1 + ·` は順序を反射する。 -/
theorem enat_one_add_lt_one_add_iff {a b : ℕ∞} : 1 + a < 1 + b ↔ a < b :=
  WithTop.add_lt_add_iff_left ENat.one_ne_top

/-- **望遠鏡和(差の形)** `σ•(Σ_{i<n} σ^i•x) − Σ_{i<n} σ^i•x = σ^n•x − x`。

★原文の「`v(σ(y) − y) = v(σ^p(x) − x)`」(`y := Σ_{i<p} σ^i(x)`)の中身。
一般の加法群と一般の `DistribMulAction` で成り立つ。 -/
theorem smul_sum_range_sub_self {M : Type*} [AddCommGroup M] {G : Type*} [Monoid G]
    [DistribMulAction G M] (σ : G) (x : M) (n : ℕ) :
    σ • (∑ i ∈ Finset.range n, σ ^ i • x) - (∑ i ∈ Finset.range n, σ ^ i • x)
      = σ ^ n • x - x := by
  have h1 : σ • (∑ i ∈ Finset.range n, σ ^ i • x) - (∑ i ∈ Finset.range n, σ ^ i • x)
      = ∑ i ∈ Finset.range n, σ ^ i • (σ • x - x) := by
    rw [Finset.smul_sum, ← Finset.sum_sub_distrib]
    refine Finset.sum_congr rfl fun i _ => ?_
    rw [smul_sub, smul_smul, smul_smul, ← pow_succ, ← pow_succ']
  rw [h1, sum_range_smul_smul_sub]

/-! ## §3 Yoshida Prop 6.6 (iii)(抽象層)

★**強帰納法の外枠は `σ` について全称量化してある**(逸脱 8)。帰納法の仮説を
`σ^p` にも `j` 未満のすべての添字にも当てる必要があるからである。 -/

/-- ★★★**Yoshida 2008 Prop 6.6 (iii) の本体(帰納法の形)**。

`σ ∈ G_1`、`orderOf σ = p^m`、`j + 1 < m` のとき

`p^{j+1} ∣ (i_{j+1} − i_j)`  (`i_k := i(σ^{p^k})` の有限代表、ℤ の中の引き算)。

★添字は原文から 1 つずらしてある(原文 `i_{j−1} ≡ i_j (mod p^j)`。逸脱 1)。
★`toNat` を経由するのは `ℕ∞` の切り詰め引き算を出さないため(逸脱 2)。
仮定 `j + 1 < m` が `i_j`・`i_{j+1}` の有限性を保証するので、`toNat` はジャンクを返さない。

証明(原文の段取りをそのまま):
* 鎖(§1 `dvd_sub_of_forall_dvd_succ`)を帰納法の仮説から作る。
* `σ^p` に帰納法の仮説を当てて `p^j ∣ i_{j+1} − i_j` を得る(`j = 0` では自明)。
* 背理法。`s′` を `s′ + i_{j+1} = i_j + E`(`E := p^{j+1}(i_{j+1}+1)`)で定める(シフト法)。
* Lemma 6.5 で `x := α_{s′}(σ^p)` を取り、`y := Σ_{i<p} σ^i x`、`z := σ y − y`。
* `y` を打ち切り π 進展開し、桁を `p^{j+1} ∣ n` と否定に割る。
* 前者は付値が `v(z)` より真に大きく(単調性)、後者は Claim により互いに相異なる。
  よって `v(Σ_{後者}) = v(z)` と `v(Σ_{後者}) ≠ v(z)` が同時に出て矛盾。 -/
theorem dvd_ramIndex_toNat_sub {A B : Type*} [CommRing A] [CommRing B] [Algebra A B]
    [IsDomain B] [IsDiscreteValuationRing B] {G : Type*} [Group G] [MulSemiringAction G B]
    [SMulCommClass G A B] [FaithfulSMul G B] (p : ℕ) (hp : p.Prime)
    [CharP (ResidueField B) p] {π : B} (hπ : Irreducible π)
    (hadj : Algebra.adjoin A ({π} : Set B) = ⊤)
    (hres : ∀ b : B, ∃ a : A, b - algebraMap A B a ∈ maximalIdeal B) (j : ℕ) :
    ∀ σ : G, σ ∈ lowerRamificationGroup B G 1 → ∀ m : ℕ, orderOf σ = p ^ m → j + 1 < m →
      (p : ℤ) ^ (j + 1) ∣ ((ramIndex π (σ ^ p ^ (j + 1))).toNat : ℤ)
        - ((ramIndex π (σ ^ p ^ j)).toNat : ℤ) := by
  classical
  have huni : maximalIdeal B = Ideal.span {π} := (irreducible_iff_uniformizer π).mp hπ
  induction j using Nat.strong_induction_on with
  | _ j IH =>
  intro σ hσ m hord hjm
  obtain ⟨iN, hiNeq⟩ : ∃ iN : ℕ → ℕ, ∀ k, iN k = (ramIndex π (σ ^ p ^ k)).toNat :=
    ⟨_, fun _ => rfl⟩
  -- 段 3: 有限代表 `iN` と `ℕ∞` の橋
  have hiN : ∀ k, k < m → (iN k : ℕ∞) = ramIndex π (σ ^ p ^ k) := by
    intro k hk
    rw [hiNeq k]
    refine ENat.coe_toNat ?_
    intro hc
    have h2 := (ramIndex_pow_pow_eq_top_iff (A := A) hp.one_lt hadj hord k).1 hc
    omega
  rw [← hiNeq (j + 1), ← hiNeq j]
  -- 段 4: 合同の鎖
  have hchain : ∀ a b : ℕ, a ≤ b → b ≤ j → (p : ℤ) ^ (a + 1) ∣ (iN b : ℤ) - (iN a : ℤ) := by
    intro a b hab hbj
    refine dvd_sub_of_forall_dvd_succ (f := fun k => (iN k : ℤ)) (K := j) ?_ hab hbj
    intro k hk
    have h := IH k hk σ hσ m hord (by omega)
    rwa [← hiNeq (k + 1), ← hiNeq k] at h
  have hconj : ∀ k : ℕ, (σ ^ p) ^ p ^ k = σ ^ p ^ (k + 1) := fun k => by
    rw [← pow_mul, ← pow_succ']
  obtain ⟨m', hm'⟩ := exists_orderOf_pow_eq hp hord p
  -- 原文「Now applying the Inductive Hypothesis to σ^p ∈ G_1」
  have hstep : (p : ℤ) ^ j ∣ (iN (j + 1) : ℤ) - (iN j : ℤ) := by
    rcases Nat.eq_zero_or_pos j with rfl | hjpos
    · simp
    obtain ⟨j', rfl⟩ : ∃ j', j = j' + 1 := ⟨j - 1, by omega⟩
    have hj'm' : j' + 1 < m' := by
      by_contra hc
      have htop : ramIndex π ((σ ^ p) ^ p ^ (j' + 1)) = ⊤ :=
        (ramIndex_pow_pow_eq_top_iff (A := A) hp.one_lt hadj hm' (j' + 1)).2 (by omega)
      rw [hconj (j' + 1)] at htop
      have h2 := hiN (j' + 1 + 1) (by omega)
      rw [htop] at h2
      exact (ENat.coe_ne_top _) h2
    have h := IH j' (by omega) (σ ^ p) (pow_mem hσ p) m' hm' hj'm'
    rw [hconj (j' + 1), hconj j', ← hiNeq (j' + 1 + 1), ← hiNeq (j' + 1)] at h
    exact h
  by_contra hcon
  -- 段 7: シフト `s := (i_j − i_{j+1}) + p^{j+1}·t`(引き算を一度も書かない)
  set E : ℕ := p ^ (j + 1) * (iN (j + 1) + 1) with hEdef
  have hEgt : iN (j + 1) < E := by
    have h1 : iN (j + 1) + 1 ≤ p ^ (j + 1) * (iN (j + 1) + 1) :=
      Nat.le_mul_of_pos_left _ (pow_pos hp.pos (j + 1))
    omega
  obtain ⟨s, hs⟩ : ∃ s : ℕ, s + iN (j + 1) = iN j + E := ⟨iN j + E - iN (j + 1), by omega⟩
  have hEdvd : (p : ℤ) ^ (j + 1) ∣ (E : ℤ) :=
    ⟨(iN (j + 1) : ℤ) + 1, by rw [hEdef]; push_cast; ring⟩
  have hsZ : (s : ℤ) = ((iN j : ℤ) - (iN (j + 1) : ℤ)) + (E : ℤ) := by
    have h := congrArg (fun n : ℕ => (n : ℤ)) hs
    push_cast at h
    linarith
  have hsdvd : p ^ j ∣ s := by
    rw [← Int.natCast_dvd_natCast]
    push_cast
    rw [hsZ]
    exact dvd_add (by simpa using dvd_neg.mpr hstep)
      (dvd_trans (pow_dvd_pow _ (by omega)) hEdvd)
  have hnsdvd : ¬ p ^ (j + 1) ∣ s := by
    intro hd
    refine hcon ?_
    have hdZ : (p : ℤ) ^ (j + 1) ∣ (s : ℤ) := by
      rw [← Int.natCast_dvd_natCast] at hd; push_cast at hd; exact hd
    have he : (iN (j + 1) : ℤ) - (iN j : ℤ) = (E : ℤ) - (s : ℤ) := by linarith
    rw [he]
    exact dvd_sub hEdvd hdZ
  -- 段 8: `x = α_s(σ^p)`、`y = Σ σ^i x`、`z = σ y − y`
  have hram_s : ramIndex π ((σ ^ p) ^ s) = (iN (j + 1) : ℕ∞) := by
    rw [ramIndex_pow_eq_ramIndex_pow_pow (A := A) p hp huni hadj hm' hsdvd hnsdvd,
      hconj j, hiN (j + 1) (by omega)]
  have hα0 : uniformizerProd (σ ^ p) π s ≠ 0 := uniformizerProd_ne_zero hπ (σ ^ p) s
  have hαval : addVal B (uniformizerProd (σ ^ p) π s) = (s : ℕ∞) :=
    addVal_uniformizerProd hπ (σ ^ p) s
  set y : B := ∑ i ∈ Finset.range p, σ ^ i • uniformizerProd (σ ^ p) π s with hydef
  have hyval : (s : ℕ∞) < addVal B y := by
    rw [← hαval, hydef]
    exact addVal_lt_addVal_sum_smul_pow_char p huni hσ hα0
  set z : B := σ • y - y with hzdef
  have hz : z = (σ ^ p) • uniformizerProd (σ ^ p) π s - uniformizerProd (σ ^ p) π s := by
    rw [hzdef, hydef]
    exact smul_sum_range_sub_self σ (uniformizerProd (σ ^ p) π s) p
  have hzval : 1 + addVal B z = ((iN j + E : ℕ) : ℕ∞) := by
    rw [hz]
    have h1 : addVal B ((σ ^ p) • uniformizerProd (σ ^ p) π s - uniformizerProd (σ ^ p) π s)
        = ramIndex (uniformizerProd (σ ^ p) π s) (σ ^ p) := rfl
    rw [h1, one_add_ramIndex_uniformizerProd hπ (σ ^ p) s, hram_s, ← Nat.cast_add, hs]
  have hzne : addVal B z ≠ ⊤ := by
    intro h
    refine (ENat.coe_ne_top (iN j + E)) ?_
    rw [← hzval, h]
    simp
  -- 段 8 続き: 打ち切り π 進展開(原文の無限和の代わり。逸脱 6)
  have hyN : y ∈ (maximalIdeal B) ^ (s + 1) := by
    refine (mem_maximalIdeal_pow_iff_le_addVal hπ (s + 1) y).mpr ?_
    push_cast
    exact Order.add_one_le_of_lt hyval
  obtain ⟨a, ha, hrem⟩ :=
    exists_digits (A := A) hπ σ hres (s + 1) (iN j + E + 1) (by omega) y hyN
  set S : B := ∑ n ∈ Finset.Ico (s + 1) (iN j + E + 1),
    algebraMap A B (a n) * uniformizerProd σ π n with hSdef
  set r : B := y - S with hrdef
  obtain ⟨zz, hzz⟩ : ∃ zz : ℕ → B, ∀ n, zz n = algebraMap A B (a n) *
      (σ • uniformizerProd σ π n - uniformizerProd σ π n) := ⟨_, fun _ => rfl⟩
  have hSsub : σ • S - S = ∑ n ∈ Finset.Ico (s + 1) (iN j + E + 1), zz n := by
    rw [hSdef, smul_digitSum, ← Finset.sum_sub_distrib]
    exact Finset.sum_congr rfl fun n _ => by rw [hzz n]; ring
  have hzsplit : z = (∑ n ∈ Finset.Ico (s + 1) (iN j + E + 1), zz n) + (σ • r - r) := by
    rw [← hSsub, hzdef, hrdef, smul_sub]
    ring
  have hunit : ∀ n : ℕ, a n ≠ 0 → IsUnit (algebraMap A B (a n)) := fun n hn =>
    (ha n).resolve_left hn
  have hzzval : ∀ n : ℕ, a n ≠ 0 → 1 + addVal B (zz n) = (n : ℕ∞) + ramIndex π (σ ^ n) := by
    intro n hn
    rw [hzz n, addVal_mul, addVal_eq_zero_iff.mpr (hunit n hn), zero_add]
    exact one_add_ramIndex_uniformizerProd hπ σ n
  -- 段 9: 桁の分割と 3 つの評価
  have hkey : ∀ n : ℕ, ¬ p ^ (j + 1) ∣ n → a n ≠ 0 →
      ∃ k : ℕ, k ≤ j ∧ (p : ℤ) ^ k ∣ (n : ℤ) ∧ ¬ (p : ℤ) ^ (k + 1) ∣ (n : ℤ) ∧
        1 + addVal B (zz n) = ((n + iN k : ℕ) : ℕ∞) := by
    intro n hnd hn0
    have hne0 : n ≠ 0 := by rintro rfl; exact hnd (dvd_zero _)
    obtain ⟨k, u₀, hu₀, hnk⟩ := Nat.exists_eq_pow_mul_and_not_dvd hne0 p hp.ne_one
    have hpk : p ^ k ∣ n := ⟨u₀, hnk⟩
    have hnpk : ¬ p ^ (k + 1) ∣ n := by
      intro hc
      rw [hnk, pow_succ] at hc
      exact hu₀ ((mul_dvd_mul_iff_left (pow_ne_zero k hp.pos.ne')).1 hc)
    have hkj : k ≤ j := by
      by_contra hc2
      exact hnd (dvd_trans (pow_dvd_pow p (by omega)) hpk)
    refine ⟨k, hkj, ?_, ?_, ?_⟩
    · have h := Int.natCast_dvd_natCast.2 hpk
      push_cast at h
      exact h
    · intro hc
      refine hnpk ?_
      rw [← Int.natCast_dvd_natCast]
      push_cast
      exact hc
    · rw [hzzval n hn0, ramIndex_pow_eq_ramIndex_pow_pow (A := A) p hp huni hadj hord hpk hnpk,
        ← hiN k (by omega)]
      push_cast
      ring
  -- 原文「If v_p(n) ≥ j and z_n ≠ 0, then v(z_n) ≥ v(y) + i_j > i_{j−1}」
  have hbig : ∀ n ∈ Finset.Ico (s + 1) (iN j + E + 1),
      ¬ (¬ p ^ (j + 1) ∣ n ∧ a n ≠ 0) → addVal B z < addVal B (zz n) := by
    intro n hn hnot
    rcases eq_or_ne (a n) 0 with h0 | h0
    · have hz0 : zz n = 0 := by rw [hzz n, h0]; simp
      rw [hz0, addVal_zero]
      exact lt_of_le_of_ne le_top hzne
    · have hdvdn : p ^ (j + 1) ∣ n := by
        by_contra hc
        exact hnot ⟨hc, h0⟩
      have hnpos : s + 1 ≤ n := (Finset.mem_Ico.1 hn).1
      have hne0 : n ≠ 0 := by omega
      obtain ⟨k, u₀, hu₀, hnk⟩ := Nat.exists_eq_pow_mul_and_not_dvd hne0 p hp.ne_one
      have hpk : p ^ k ∣ n := ⟨u₀, hnk⟩
      have hnpk : ¬ p ^ (k + 1) ∣ n := by
        intro hc
        rw [hnk, pow_succ] at hc
        exact hu₀ ((mul_dvd_mul_iff_left (pow_ne_zero k hp.pos.ne')).1 hc)
      have hkj : j + 1 ≤ k := by
        by_contra hc2
        exact hnpk (dvd_trans (pow_dvd_pow p (by omega)) hdvdn)
      have hmono : ramIndex π (σ ^ p ^ (j + 1)) ≤ ramIndex π (σ ^ p ^ k) :=
        monotone_ramIndex_pow_pow (A := A) p hp huni hadj hσ hord hkj
      have hlow : ((iN j + E + 1 : ℕ) : ℕ∞) ≤ 1 + addVal B (zz n) := by
        rw [hzzval n h0, ramIndex_pow_eq_ramIndex_pow_pow (A := A) p hp huni hadj hord hpk hnpk]
        have hb : ((s + 1 : ℕ) : ℕ∞) ≤ (n : ℕ∞) := by exact_mod_cast hnpos
        have hc2 : (iN (j + 1) : ℕ∞) ≤ ramIndex π (σ ^ p ^ k) := by
          rw [hiN (j + 1) (by omega)]; exact hmono
        have heq : ((iN j + E + 1 : ℕ) : ℕ∞) = ((s + 1 : ℕ) : ℕ∞) + (iN (j + 1) : ℕ∞) := by
          rw [← Nat.cast_add]
          congr 1
          omega
        rw [heq]
        exact add_le_add hb hc2
      have hstrict : ((iN j + E : ℕ) : ℕ∞) < ((iN j + E + 1 : ℕ) : ℕ∞) := by
        exact_mod_cast Nat.lt_succ_self (iN j + E)
      rw [← enat_one_add_lt_one_add_iff, hzval]
      exact lt_of_lt_of_le hstrict hlow
  -- Claim(a): `i_j + E` は `n + i(σ^n)`(`v_p(n) < j+1`)のどれとも異なる
  have hS1ne : ∀ n : ℕ, (¬ p ^ (j + 1) ∣ n ∧ a n ≠ 0) → addVal B (zz n) ≠ addVal B z := by
    intro n hmem heq
    obtain ⟨k, hkj, hpkZ, hnpkZ, hval⟩ := hkey n hmem.1 hmem.2
    have h1 : ((n + iN k : ℕ) : ℕ∞) = ((iN j + E : ℕ) : ℕ∞) := by
      rw [← hval, ← hzval, heq]
    have h2 : n + iN k = iN j + E := by exact_mod_cast h1
    exact add_ne_add_of_not_dvd (hchain k j hkj (le_refl j))
      (dvd_trans (pow_dvd_pow _ (by omega)) hEdvd) hnpkZ (by exact_mod_cast h2)
  -- Claim(b): `n + i(σ^n)`(`v_p(n) < j+1`)は互いに相異なる
  have hpair : ∀ n : ℕ, (¬ p ^ (j + 1) ∣ n ∧ a n ≠ 0) → ∀ n' : ℕ,
      (¬ p ^ (j + 1) ∣ n' ∧ a n' ≠ 0) → n ≠ n' → addVal B (zz n) ≠ addVal B (zz n') := by
    intro n hm1 n' hm2 hne heq
    obtain ⟨k, hkj, hpkZ, hnpkZ, hval⟩ := hkey n hm1.1 hm1.2
    obtain ⟨k', hk'j, hpk'Z, hnpk'Z, hval'⟩ := hkey n' hm2.1 hm2.2
    have h1 : ((n + iN k : ℕ) : ℕ∞) = ((n' + iN k' : ℕ) : ℕ∞) := by
      rw [← hval, ← hval', heq]
    have h2 : n + iN k = n' + iN k' := by exact_mod_cast h1
    have h3 : (n : ℤ) + (iN k : ℤ) = (n' : ℤ) + (iN k' : ℤ) := by exact_mod_cast h2
    have h4 := eq_of_add_eq_add_of_chain (f := fun t => (iN t : ℤ)) (J := j) hchain hkj hk'j
      hpkZ hnpkZ hpk'Z hnpk'Z h3
    exact hne (by exact_mod_cast h4)
  -- 段 10: 矛盾
  have hrval : ((iN j + E + 1 : ℕ) : ℕ∞) ≤ addVal B r :=
    (mem_maximalIdeal_pow_iff_le_addVal hπ (iN j + E + 1) r).1 hrem
  have hru : addVal B z < addVal B (σ • r - r) := by
    have h1 : addVal B r + ((1 : ℕ) : ℕ∞) ≤ addVal B (σ • r - r) :=
      addVal_add_le_addVal_smul_sub_self huni hσ r
    have h2 : addVal B z ≤ ((iN j + E : ℕ) : ℕ∞) := by
      rw [← hzval]; exact le_add_self
    have h3 : ((iN j + E : ℕ) : ℕ∞) < ((iN j + E + 1 : ℕ) : ℕ∞) := by
      exact_mod_cast Nat.lt_succ_self (iN j + E)
    have h4 : addVal B r ≤ addVal B r + ((1 : ℕ) : ℕ∞) := le_self_add
    exact lt_of_le_of_lt h2 (lt_of_lt_of_le h3 (le_trans hrval (le_trans h4 h1)))
  have hfil := Finset.sum_filter_add_sum_filter_not
    (Finset.Ico (s + 1) (iN j + E + 1)) (fun n => ¬ p ^ (j + 1) ∣ n ∧ a n ≠ 0) zz
  set T₁ := Finset.filter (fun n => ¬ p ^ (j + 1) ∣ n ∧ a n ≠ 0)
    (Finset.Ico (s + 1) (iN j + E + 1)) with hT₁def
  set T₂ := Finset.filter (fun n => ¬ (¬ p ^ (j + 1) ∣ n ∧ a n ≠ 0))
    (Finset.Ico (s + 1) (iN j + E + 1)) with hT₂def
  set u : B := (∑ n ∈ T₂, zz n) + (σ • r - r) with hudef
  have hzu : z = (∑ n ∈ T₁, zz n) + u := by
    rw [hudef, ← add_assoc, hfil]
    exact hzsplit
  have hsum2 : addVal B z < addVal B (∑ n ∈ T₂, zz n) := by
    refine lt_map_sum (addVal B) (by rw [addVal_zero]; exact lt_of_le_of_ne le_top hzne) ?_
    intro n hn
    rw [hT₂def, Finset.mem_filter] at hn
    exact hbig n hn.1 hn.2
  have hu : addVal B z < addVal B u := by
    rw [hudef]
    exact lt_of_lt_of_le (lt_min hsum2 hru) (AddValuation.map_add _ _ _)
  have hne1 : addVal B (∑ n ∈ T₁, zz n) ≠ addVal B z := by
    refine map_sum_ne_of_pairwise_ne (addVal B) ?_ ?_ ?_
    · intro i hi i' hi' hii
      rw [hT₁def, Finset.mem_filter] at hi hi'
      exact hpair i hi.2 i' hi'.2 hii
    · intro i hi
      rw [hT₁def, Finset.mem_filter] at hi
      exact hS1ne i hi.2
    · rw [addVal_zero]
      exact fun h => hzne h.symm
  have hle1 : addVal B z ≤ addVal B (∑ n ∈ T₁, zz n) := by
    have heq : (∑ n ∈ T₁, zz n) = z - u := by rw [hzu]; ring
    rw [heq]
    refine le_trans ?_ (AddValuation.map_sub (addVal B) z u)
    exact le_min (le_refl _) (le_of_lt hu)
  have hle2 : addVal B (∑ n ∈ T₁, zz n) ≤ addVal B z := by
    by_contra hc
    rw [not_le] at hc
    have h := lt_of_lt_of_le (lt_min hc hu) (AddValuation.map_add (addVal B) (∑ n ∈ T₁, zz n) u)
    rw [← hzu] at h
    exact lt_irrefl _ h
  exact hne1 (le_antisymm hle2 hle1)

/-- ★★★**Yoshida 2008 Prop 6.6 (iii)** —— `i_j = c`、`i_{j+1} = d`(ともに有限)なら

`p^{j+1} ∣ (d − c)`  (ℤ の中の引き算)。

★`toNat` を statement から追い出した形。仮定 `hd` が `i_{j+1} ≠ ⊤` を含むので
`j + 1 < m` が従い、`hc` は `c` を一意に決める。
★★`ℕ∞` の引き算で書くと空虚に真になる(退化検査 D1)ため、有限代表 `c, d : ℕ` を
経由して `ℤ` で書いてある(逸脱 2)。 -/
theorem ramIndex_pow_pow_congr {A B : Type*} [CommRing A] [CommRing B] [Algebra A B]
    [IsDomain B] [IsDiscreteValuationRing B] {G : Type*} [Group G] [MulSemiringAction G B]
    [SMulCommClass G A B] [FaithfulSMul G B] (p : ℕ) (hp : p.Prime)
    [CharP (ResidueField B) p] {π : B} (hπ : Irreducible π)
    (hadj : Algebra.adjoin A ({π} : Set B) = ⊤)
    (hres : ∀ b : B, ∃ a : A, b - algebraMap A B a ∈ maximalIdeal B) {σ : G}
    (hσ : σ ∈ lowerRamificationGroup B G 1) {m : ℕ} (hord : orderOf σ = p ^ m)
    {j c d : ℕ} (hc : ramIndex π (σ ^ p ^ j) = (c : ℕ∞))
    (hd : ramIndex π (σ ^ p ^ (j + 1)) = (d : ℕ∞)) :
    (p : ℤ) ^ (j + 1) ∣ (d : ℤ) - (c : ℤ) := by
  have hjm : j + 1 < m := by
    by_contra hcon
    have htop : ramIndex π (σ ^ p ^ (j + 1)) = ⊤ :=
      (ramIndex_pow_pow_eq_top_iff (A := A) hp.one_lt hadj hord (j + 1)).2 (by omega)
    rw [hd] at htop
    exact (ENat.coe_ne_top d) htop
  have h := dvd_ramIndex_toNat_sub (A := A) p hp hπ hadj hres j σ hσ m hord hjm
  rw [hc, hd] at h
  simpa using h

/-- ★★★**Yoshida 2008 Prop 6.6 (iii)(原文の `∞` 規約まで込めた形)**。

`i_{j+1} = ∞` か、さもなくば `i_j`・`i_{j+1}` はともに有限で `p^{j+1} ∣ (i_{j+1} − i_j)`。

★原文の「where `∞` is understood to be congruent to any integer」を場合分けで表した
(逸脱 2 の後半)。★`i_{j+1} = ∞` の側は `⟺ m ≤ j+1`(Y7a `ramIndex_pow_pow_eq_top_iff`)。 -/
theorem ramIndex_pow_pow_congr_or_top {A B : Type*} [CommRing A] [CommRing B] [Algebra A B]
    [IsDomain B] [IsDiscreteValuationRing B] {G : Type*} [Group G] [MulSemiringAction G B]
    [SMulCommClass G A B] [FaithfulSMul G B] (p : ℕ) (hp : p.Prime)
    [CharP (ResidueField B) p] {π : B} (hπ : Irreducible π)
    (hadj : Algebra.adjoin A ({π} : Set B) = ⊤)
    (hres : ∀ b : B, ∃ a : A, b - algebraMap A B a ∈ maximalIdeal B) {σ : G}
    (hσ : σ ∈ lowerRamificationGroup B G 1) {m : ℕ} (hord : orderOf σ = p ^ m) (j : ℕ) :
    ramIndex π (σ ^ p ^ (j + 1)) = ⊤ ∨
      ∃ c d : ℕ, ramIndex π (σ ^ p ^ j) = (c : ℕ∞) ∧
        ramIndex π (σ ^ p ^ (j + 1)) = (d : ℕ∞) ∧ (p : ℤ) ^ (j + 1) ∣ (d : ℤ) - (c : ℤ) := by
  rcases eq_or_ne (ramIndex π (σ ^ p ^ (j + 1))) ⊤ with h | h
  · exact Or.inl h
  refine Or.inr ?_
  have hjm : j + 1 < m := by
    by_contra hcon
    exact h ((ramIndex_pow_pow_eq_top_iff (A := A) hp.one_lt hadj hord (j + 1)).2 (by omega))
  have hcfin : ramIndex π (σ ^ p ^ j) ≠ ⊤ := by
    intro hcc
    have := (ramIndex_pow_pow_eq_top_iff (A := A) hp.one_lt hadj hord j).1 hcc
    omega
  refine ⟨(ramIndex π (σ ^ p ^ j)).toNat, (ramIndex π (σ ^ p ^ (j + 1))).toNat,
    (ENat.coe_toNat hcfin).symm, (ENat.coe_toNat h).symm, ?_⟩
  exact dvd_ramIndex_toNat_sub (A := A) p hp hπ hadj hres j σ hσ m hord hjm

/-! ## §4 `PAdicLocalField` への具体化

`B := adjoinIntegers K x`、`A := 𝒪[K.carrier]`、`G := Gal(K(x)/K)`。★Y7a §7 と同型で、
仮定の束はそのまま流用できる。`hres`(剰余体が伸びない)は完全分岐から出る。 -/

variable {p : ℕ} [Fact p.Prime]

attribute [local instance] isDiscreteValuationRing_adjoinIntegers

/-- ★★★**Yoshida 2008 Prop 6.6 (iii)(`PAdicLocalField` 版)**。

`orderOf σ = p^m` は原文の "(by Proposition 6.2)"
(`Found/PGC/RamificationJumpDivisibility.lean` が供給する)。 -/
theorem ramIndex_pow_pow_congr_adjoin (K : PAdicLocalField p) (x : K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (ht : IsTotallyRamifiedAdjoin K x) {π : adjoinIntegers K x} (hπ : Irreducible π)
    {σ : (IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure))}
    (hσ : σ ∈ lowerRamificationGroupAdjoin K x 1) {j c d : ℕ}
    (hc : ramIndex π (σ ^ p ^ j) = (c : ℕ∞))
    (hd : ramIndex π (σ ^ p ^ (j + 1)) = (d : ℕ∞)) :
    (p : ℤ) ^ (j + 1) ∣ (d : ℤ) - (c : ℤ) := by
  haveI := charP_residueField_adjoinIntegers K x
  have huni := (irreducible_iff_uniformizer π).mp hπ
  have hadj := adjoin_uniformizer_eq_top_adjoinIntegers K x ht huni
  obtain ⟨m, hm⟩ := exists_orderOf_eq_pow_of_mem_lowerRamificationGroup_one
    (A := 𝒪[K.carrier]) p huni hπ.ne_zero hadj hσ
  exact ramIndex_pow_pow_congr (A := 𝒪[K.carrier]) p Fact.out hπ hadj
    (exists_sub_algebraMap_mem_maximalIdeal_of_isTotallyRamifiedAdjoin K x ht) hσ hm hc hd

/-- ★★★**Yoshida 2008 Prop 6.6 (iii)(`PAdicLocalField` 版・`∞` 規約込み)**。 -/
theorem ramIndex_pow_pow_congr_or_top_adjoin (K : PAdicLocalField p) (x : K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (ht : IsTotallyRamifiedAdjoin K x) {π : adjoinIntegers K x} (hπ : Irreducible π)
    {σ : (IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure))}
    (hσ : σ ∈ lowerRamificationGroupAdjoin K x 1) (j : ℕ) :
    ramIndex π (σ ^ p ^ (j + 1)) = ⊤ ∨
      ∃ c d : ℕ, ramIndex π (σ ^ p ^ j) = (c : ℕ∞) ∧
        ramIndex π (σ ^ p ^ (j + 1)) = (d : ℕ∞) ∧ (p : ℤ) ^ (j + 1) ∣ (d : ℤ) - (c : ℤ) := by
  haveI := charP_residueField_adjoinIntegers K x
  have huni := (irreducible_iff_uniformizer π).mp hπ
  have hadj := adjoin_uniformizer_eq_top_adjoinIntegers K x ht huni
  obtain ⟨m, hm⟩ := exists_orderOf_eq_pow_of_mem_lowerRamificationGroup_one
    (A := 𝒪[K.carrier]) p huni hπ.ne_zero hadj hσ
  exact ramIndex_pow_pow_congr_or_top (A := 𝒪[K.carrier]) p Fact.out hπ hadj
    (exists_sub_algebraMap_mem_maximalIdeal_of_isTotallyRamifiedAdjoin K x ht) hσ hm j

end ABC3.Found.PGC
