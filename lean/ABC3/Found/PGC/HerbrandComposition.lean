import ABC3.Found.PGC.HerbrandFunction

/-!
# Herbrand 関数の閉じた形と合成則(Yoshida 2008 Lemma 6.10)

典拠: T. Yoshida, *Local Class Field Theory via Lubin-Tate Theory*
(Ann. Fac. Sci. Toulouse 17-2, 2008; arXiv math/0606108) **Lemma 6.10**(物理 p.16)。
構造化済み原文は `ResearchPaper/1_Structured/Local Class Field Theory via Lubin-Tate Theory/`
`section-6.html` の `id="lemma-6-10"`(`data-pdf-page="16"`, `data-item="Lemma 6.10"`)。

原文 (Yoshida08 p.16):
> Lemma 6.10. Let φ_G(n) := −1 + 1|G| _τ∈G min{i(τ), n + 1} for n ∈ R[bb]_≥0. Then: (i) φ_G(0) = 0, φ_G(n) = 1|G| ^n_i=1 |G_i| for n ∈ Z[bb]_≥1. (ii) φ_G = φ_G/H ◦ φ_H on R[bb]_≥0.

★`pdftotext` は分数の横線と総和記号 `Σ` を出力しない(`data-txt` が `1|G|` /
空文字で明示している)。また合成記号は原典が `∘`、抽出結果は `◦`(U+25E6)である。読み下すと

  (i) `φ_G(0) = 0`、`φ_G(n) = (1/|G|) Σ_{i=1}^{n} |G_i|`(`n ∈ ℤ≥1`)、
  (ii) `φ_G = φ_{G/H} ∘ φ_H`(`ℝ≥0` の上で)。

設定は §6.2 (`#setup-6-2`):
> 6.2. The Hasse-Arf theorem. Let G = Gal(K′/K) with K′/K totally ramified as before,
> and let G ▷ H with G/H = Gal(K′′/K). For σ ∈ G, let σ[bar] = σH ∈ G/H be its image.

## 原典の Proof(`.txt` 1077–1101 行)と本ファイルの対応

> Proof. (i): Σ_{τ∈G} min{i(τ), n + 1} = Σ_{i=0}^{n−1} ( Σ_{τ∈G_i∖G_{i+1}} (i + 1) )
> + Σ_{τ∈G_n} (n + 1) = Σ_{i=0}^{n} |G_i|.
> (ii): As φ(0) = 0 and φ is continuous and piecewise linear, we only need to compare the
> derivatives of both sides at n ∈ (i − 1, i) for i ∈ Z>0. For LHS it is |G_n|/|G|, and for RHS
> it is (|(G/H)_{φ_H(n)}|/|G/H|) · (|H_n|/|H|) = |G_nH/H||H_n|/|G| = |G_n|/|G| by Proposition 6.9
> and G_n/H_n = G_n/(H ∩ G_n) ≅ G_nH/H.

### (i) —— 原文の層分けをそのまま「純組合せ」に切り出した

★★原文の第 1 式には分岐も付値も Galois も出てこない。出てくるのは
「有限集合 `S` と `f : S → ℕ∞`」だけである。実際

  `min{f(τ), n+1} = #{i ∈ [0, n] | f(τ) ≥ i+1}`

なので、`Σ_{τ∈S} min{f(τ), n+1}` は 2 重和の入れ替え(`Finset.sum_comm`)だけで
`Σ_{i=0}^{n} #{τ | f(τ) ≥ i+1}` になる。§2 の `sum_truncENat_natCast` がそれである。
★原文は `G_i ∖ G_{i+1}` で層に分けるが、**「`i` を数える」方向に読み替えると差集合が要らない**。

### (ii) —— ★★原文の「微分を比べる」を使わなかった(第 3 の道)

★**道 A(`deriv` + 区分線形関数の一致判定)も道 B(区間ごとの傾き)も採らなかった。**
理由は次の観察である。掛け算形で書くと (ii) は

  `Σ_{σ∈G} min{i(σ), n+1} = Σ_{σ∈G} min{i_ϖ(σ), φ_H(n)+1}`   (★両辺とも `G` 上の和)

と同値であり、**両辺を `H` の剰余類ごとに比べると、min の結合則だけで一致する**。
実際、剰余類 `σH` の中で `i` が最大の代表を `σρ`(`m := i(σρ)`)と取ると

* 左辺の寄与 `Σ_{τ∈H} min{i(στ), n+1} = Σ_{τ∈H} min{min{i(τ), m}, n+1}`
  (Proposition 6.9 の中核 `ultrametric_mul_eq_min`)
  `= Σ_{τ∈H} min{i(τ), min{m, n+1}}`,
* 右辺の寄与 `|H| · min{i_ϖ(σ), φ_H(n)+1} = |H| · (φ_H(min{m−1, n}) + 1)`
  (Lemma 6.8 と `φ_H` の単調性)
  `= Σ_{τ∈H} min{i(τ), min{m, n+1}}`(掛け算形 `card_mul_herbrandPhi_add_one`)

で、**同じ式**になる。★`min` の結合則 `truncENat_min` 1 本が原文の「微分の一致」を置き換えている。

★★**したがって原文が使う第 2 同型定理 `G_n/H_n = G_n/(H ∩ G_n) ≅ G_nH/H` は要らなかった。**
`|G_nH/H|` を Y11 の `(G_n : Set G) * (H : Set G)` から取り出す必要も生じていない
(取り出し方が難しいのではなく、**この道では現れない**)。

★剰余類への分割そのものも要らない。`Σ_{σ∈G} Σ_{τ∈H}` を `Finset.sum_comm` で入れ替え、
内側の `σ ↦ στ` が `G` の全単射であることを使うと `|H| · Σ_{σ∈G}` になる
(`herbrandPhiGroup_comp` の `hfub`)。★**商群 `G ⧸ H` を作らずに閉じる。**

### mathlib に「区分線形関数の一致判定」は在ったか

★**探す前に不要になった。** 上のとおり (ii) は `deriv` も `Continuous` も使わずに閉じたので、
`Monotone.map_min`(`φ_H(min a b) = min (φ_H a) (φ_H b)`)以外の解析的な道具を使っていない。

## `φ_G` と `φ_{G/H}` をどう書いたか

Y11 の `herbrandPhi α H n` は部分群 `H : Subgroup G` に対する `φ_H` である。
`Fintype ↥(⊤ : Subgroup G)` は `[Fintype G]` からは推論されない(実測)ので、
`φ_G` は `G` 全体の上の和として `herbrandPhiGroup G α n` を新たに定義した。
両者が同じものであることは `herbrandPhiGroup_eq_herbrandPhi_top` で確かめてある。

`φ_{G/H}` は **`herbrandPhiGroup G ϖ`**(`G` 上の和、分母 `|G|`)として書いた。
これが商群の上の `φ` と一致する根拠は §4 の

* `phiOf_comp_of_card_fiber` : `g : S → T` のファイバーの個数が一定なら `phiOf (f ∘ g) = phiOf f`
  (★純組合せ。群も分岐も出てこない)
* `phiOf_quotient` : その `g := QuotientGroup.mk : G → G ⧸ H` の場合
* `herbrandPhiGroup_eq_phiOf_quotient` : `i_ϖ` の `G ⧸ H` への降下 `f` を任意に取ると
  `herbrandPhiGroup G ϖ n = phiOf f n`(★降下の存在は `exists_quotient_ramIndex`)

である。★つまり `|G| = |H| · |G/H|` と「`i_ϖ` は `H` 剰余類上で定数」の 2 つで、
`G` 上の和と `G/H` 上の和は `φ` として同じ値を与える。

## 本ファイルの構成

**§1 抽象核 —— `truncENat` の代数**(`ℕ∞` と `ℝ` だけ)
* `truncENat_min` : ★`min{min{x,y}, r} = min{x, min{y, r}}`。(ii) の要。

**§2 抽象核 —— 有限和の層分け**(純組合せ。★分岐・付値・群すら出てこない)
* `card_filter_range_lt` / `truncENat_natCast_succ` : `min{x, n+1} = #{i ≤ n | x ≥ i+1}`。
* `sum_card_filter_comm` : 2 重和の入れ替え(原文の層分けの正体)。
* `phiOf f n := −1 + (1/|S|) Σ_{τ∈S} min{f(τ), n+1}` : ★`φ` の抽象な器。
* `phiOf_natCast` : `φ(n) = −1 + (1/|S|) Σ_{i=0}^{n} c_i`(★**添字は 0 から**)。
* `phiOf_natCast_of_pos` : `f > 0` なら `c_0 = |S|` で `−1` と相殺し、
  ★**`φ(n) = (1/|S|) Σ_{i=1}^{n} c_i`(添字は 1 から)**。
* `phiOf_zero` : `φ(0) = 0`。

**§3 抽象核 —— 剰余類上の和**(純群論。★分岐・付値・Galois の語彙が 1 つも出てこない)
* `sum_truncENat_coset` : `Σ_{τ∈H} min{f(στ), r} = Σ_{τ∈H} min{f(τ), min{f(σρ), r}}`。

**§4 抽象核 —— ファイバーが一定なら `phiOf` は不変**(純組合せ)
* `phiOf_comp_of_card_fiber` / `phiOf_quotient`。

**§5 具体層 —— `φ_H` / `φ_G` と `phiOf` の橋**
* `herbrandPhi_eq_phiOf` / `herbrandPhiGroup` / `herbrandPhiGroup_eq_herbrandPhi_top`。
* `pos_ramIndex` : ★`i(σ) ≥ 1`(退化検査。これが無いと (i) の添字が 1 ずれる)。

**§6 Lemma 6.10 (i)**
* `herbrandPhiGroup_zero` / `herbrandPhiGroup_natCast` : ★★**主結論 (i)**。
* `herbrandPhi_zero` / `herbrandPhi_natCast` : 一般の `H` 版((ii) の下流が使う)。

**§7 具体層 —— 剰余類ごとの等式**
* `ramIndex_mul_mem_eq` : `i_ϖ(στ) = i_ϖ(σ)`(`τ ∈ H`)。
* `coset_sum_truncENat` : ★★(ii) の心臓。Lemma 6.8 と Proposition 6.9 の中核を使う。

**§8 Lemma 6.10 (ii)**
* `herbrandPhiGroup_comp` : ★★★★**主結論 (ii)** `φ_G = φ_{G/H} ∘ φ_H`。
* `exists_quotient_ramIndex` / `herbrandPhiGroup_eq_phiOf_quotient` : `φ_{G/H}` の同定。

## ★★逸脱の記録

1. ★**原典は `n ∈ ℝ≥0` に制限しているが、本ファイルは `n : ℝ` 全体で述べている。**
   Y11 の `herbrandPhi` / `ramificationGroupReal` が既に `ℝ` 全体で定義されており、
   (ii) の証明は `n` の符号を一切使わない。★弱めても強めてもいない。
   (i) は `n : ℕ` についてなので、そもそも `n ≥ 0` である。
2. ★**原典の (i) は `n ∈ ℤ≥1` だが、本ファイルは `n : ℕ`(`n = 0` を含む)で述べている。**
   `n = 0` のとき `Finset.Icc 1 0 = ∅` で右辺は `0` になり、(i) の前半 `φ_G(0) = 0` に一致する。
   ★**原典の 2 つの主張が 1 本の式に収まった**(強めている。弱めてはいない)。
3. ★★**(ii) で `H` の正規性 `G ▷ H` を仮定していない。**
   本ファイルの `φ_{G/H}` は「`G` 上の和で書いた `φ`」なので商群を作らず、
   計算(剰余類ごとの比較)は正規性を使わない。
   ★正規性が要るのは `φ_{G/H}` を**商群の上の関数として**同定する補助定理
   `phiOf_quotient` / `exists_quotient_ramIndex` / `herbrandPhiGroup_eq_phiOf_quotient` だけで、
   そこには `[H.Normal]` を置いてある。★**主結論 (ii) は正規性なしで真である。**
4. ★**原典の (ii) の証明(微分の比較)は採らなかった。** 上の「第 3 の道」を参照。
   ★原文が使う第 2 同型定理 `G_n/H_n ≅ G_nH/H` は本ファイルに現れない。
5. (ii) の仮定 `hcomp` `hHtriv` `hπ'` `hinj` `hfixC` `hres` `hAC` `hϖ` `hadj` `hfix` は
   **すべて Y10(Lemma 6.8)・Y11(Proposition 6.9)のものをそのまま引き継いだ**。
   ★本ファイルが新たに足した仮定は `[Fintype G]` の 1 つだけである
   (原典は `G = Gal(K′/K)` が有限であることを前提にしている)。

## 退化の自己検査

* ★**`|G|` が 0(`G` 無限)だと除算が壊れる**。(ii) は `[Fintype G]` を要求する。
  (i) の一般 `H` 版は `[Fintype ↥H]` を要求する(Y11 と同じ)。
* ★**`i(1) = ⊤` を落とすと和が壊れる**(`ramIndex_one`)。`φ_H` の狭義単調性
  (`strictMono_herbrandPhi`、(ii) の `min` の交換で使う)がこれ 1 本に依っている。
* ★★**(i) の途中式の添字は `Σ_{i=0}^{n}`、結論の添字は `Σ_{i=1}^{n}` である。**
  `|G_0| = |G|` が `−1` と相殺する。相殺の根拠は `pos_ramIndex`(`i(σ) ≥ 1`、すなわち
  `σ•α − α ∈ 𝔪_B`)であり、これは `huni : 𝔪_B = (α)` から出る。
  ★**`pos_ramIndex` を落とすと結論が 1 ずれる**(`phiOf_natCast` と
  `phiOf_natCast_of_pos` の差がちょうどそれである)。
* ★**(ii) で `φ_H` の単調性を落とすと `min{i_ϖ(σ), φ_H(n)+1} = φ_H(min{m−1,n})+1`
  が出ない**(`Monotone.map_min`)。
* ★`ENat.mul_top` の穴: `i_ϖ(σ) ≠ ⊤` を出すところで `0 < |H|` を使う(Y11 と同じ)。
* ★`ℕ∞` の切り詰め引き算・除算は書いていない(#102)。除算は `phiOf` / `herbrandPhi` の
  定義式だけに現れる。
-/

namespace ABC3.Found.PGC

open IsLocalRing IsDiscreteValuationRing

/-! ## §1 抽象核 —— `truncENat` の代数

★`ℕ∞` と `ℝ` しか出てこない。 -/

/-- ★★**`min` の結合則**(切り詰めの合成) `min{min{x,y}, r} = min{x, min{y, r}}`。

★これが Lemma 6.10 (ii) で原文の「微分の比較」を置き換える 1 本である。 -/
theorem truncENat_min (x y : ℕ∞) (r : ℝ) :
    truncENat (min x y) r = truncENat x (truncENat y r) := by
  cases x with
  | top => simp
  | coe k =>
    cases y with
    | top => simp
    | coe l =>
      rcases le_total k l with h | h
      · rw [min_eq_left (by exact_mod_cast h : (k : ℕ∞) ≤ (l : ℕ∞))]
        simp only [truncENat_coe]
        rw [← min_assoc, min_eq_left (by exact_mod_cast h : (k : ℝ) ≤ (l : ℝ))]
      · rw [min_eq_right (by exact_mod_cast h : (l : ℕ∞) ≤ (k : ℕ∞))]
        simp only [truncENat_coe]
        rw [← min_assoc, min_eq_right (by exact_mod_cast h : (l : ℝ) ≤ (k : ℝ))]

/-! ## §2 抽象核 —— 有限和の層分け(Lemma 6.10 (i) の正体)

★★分岐も付値も Galois も、群さえも出てこない。有限集合 `S` と `f : S → ℕ∞` だけである。 -/

/-- **原文の層分けの核** —— `#{i ∈ [0, n] | f > i} = min{f, n+1}`。

★原文は `G_i ∖ G_{i+1}` で層に分けるが、「`i` の側を数える」と差集合が要らない。 -/
theorem card_filter_range_lt (x : ℕ∞) (n : ℕ) :
    ((Finset.range (n + 1)).filter (fun i : ℕ => (i : ℕ∞) < x)).card
      = (min x ((n + 1 : ℕ) : ℕ∞)).toNat := by
  cases x with
  | top =>
      rw [min_eq_right le_top, ENat.toNat_coe,
        Finset.filter_true_of_mem (fun i _ => (ENat.coe_lt_top i)), Finset.card_range]
  | coe k =>
      have h1 : (Finset.range (n + 1)).filter (fun i : ℕ => (i : ℕ∞) < (k : ℕ∞))
          = Finset.range (min (n + 1) k) := by
        ext i; simp
      have h2 : min ((k : ℕ∞)) ((n + 1 : ℕ) : ℕ∞) = ((min k (n + 1) : ℕ) : ℕ∞) := by
        rcases le_total k (n + 1) with h | h
        · rw [min_eq_left (by exact_mod_cast h : (k : ℕ∞) ≤ ((n + 1 : ℕ) : ℕ∞)), min_eq_left h]
        · rw [min_eq_right (by exact_mod_cast h : ((n + 1 : ℕ) : ℕ∞) ≤ (k : ℕ∞)), min_eq_right h]
      rw [h1, Finset.card_range, h2, ENat.toNat_coe, min_comm]

/-- `card_filter_range_lt` の実数値版(`truncENat` で書いた形)。 -/
theorem truncENat_natCast_succ (x : ℕ∞) (n : ℕ) :
    truncENat x ((n : ℝ) + 1)
      = ((((Finset.range (n + 1)).filter (fun i : ℕ => (i : ℕ∞) < x)).card : ℕ) : ℝ) := by
  rw [card_filter_range_lt, toNat_min_natCast]; norm_num

/-- **2 重和の入れ替え** —— 原文の「層に分けて数える」の正体。 -/
theorem sum_card_filter_comm {S : Type*} [Fintype S] (P : ℕ → S → Prop)
    [∀ i s, Decidable (P i s)] (t : Finset ℕ) :
    ∑ τ : S, (t.filter (fun i => P i τ)).card
      = ∑ i ∈ t, (Finset.univ.filter (fun τ : S => P i τ)).card := by
  simp only [Finset.card_filter]; exact Finset.sum_comm

/-- ★★**`φ` の抽象な器** `φ_f(n) := −1 + (1/|S|) Σ_{τ∈S} min{f(τ), n+1}`。

★原典の `φ_G` も `φ_H` も `φ_{G/H}` も、すべてこの 1 つの定義の実例である。
`S` は有限型、`f : S → ℕ∞`。分岐・付値・Galois の語彙は 1 つも要らない。 -/
noncomputable def phiOf {S : Type*} [Fintype S] (f : S → ℕ∞) (n : ℝ) : ℝ :=
  -1 + (∑ τ : S, truncENat (f τ) (n + 1)) / (Nat.card S : ℝ)

/-- ★★**原文 (i) の第 1 式**(層分け) —— `Σ_{τ∈S} min{f(τ), n+1} = Σ_{i=0}^{n} #{τ | f(τ) > i}`。

★★**添字は `0` から**である(原文の途中式も `Σ_{i=0}^{n} |G_i|`)。 -/
theorem sum_truncENat_natCast {S : Type*} [Fintype S] (f : S → ℕ∞) (n : ℕ) :
    ∑ τ : S, truncENat (f τ) ((n : ℝ) + 1)
      = ∑ i ∈ Finset.range (n + 1), (Nat.card {τ : S // (i : ℕ∞) < f τ} : ℝ) := by
  classical
  simp only [fun τ : S => truncENat_natCast_succ (f τ) n]
  rw [← Nat.cast_sum, sum_card_filter_comm (fun i (τ : S) => (i : ℕ∞) < f τ), Nat.cast_sum]
  exact Finset.sum_congr rfl fun i _ => by
    rw [Nat.card_eq_fintype_card, Fintype.card_subtype]

/-- `φ_f(n) = −1 + (1/|S|) Σ_{i=0}^{n} #{τ | f(τ) > i}`(★添字は `0` から)。 -/
theorem phiOf_natCast {S : Type*} [Fintype S] (f : S → ℕ∞) (n : ℕ) :
    phiOf f (n : ℝ)
      = -1 + (∑ i ∈ Finset.range (n + 1), (Nat.card {τ : S // (i : ℕ∞) < f τ} : ℝ))
          / (Nat.card S : ℝ) := by
  rw [phiOf, sum_truncENat_natCast]

/-- ★★★**Lemma 6.10 (i) の抽象核** —— `f > 0` なら
`φ_f(n) = (1/|S|) Σ_{i=1}^{n} #{τ | f(τ) > i}`(★**添字は `1` から**)。

★`i = 0` の項は `#{τ | f(τ) > 0} = |S|` で、定義の `−1` とちょうど相殺する。
★★**`f > 0` を落とすと結論が 1 ずれる。** 具体層ではこれが `i(σ) ≥ 1`
(`σ•α − α ∈ 𝔪_B`)にあたる(`pos_ramIndex`)。 -/
theorem phiOf_natCast_of_pos {S : Type*} [Fintype S] [Nonempty S] (f : S → ℕ∞)
    (hf : ∀ τ, 0 < f τ) (n : ℕ) :
    phiOf f (n : ℝ)
      = (∑ i ∈ Finset.Icc 1 n, (Nat.card {τ : S // (i : ℕ∞) < f τ} : ℝ)) / (Nat.card S : ℝ) := by
  have hcard : (Nat.card S : ℝ) ≠ 0 := Nat.cast_ne_zero.2 Nat.card_pos.ne'
  have h0 : Nat.card {τ : S // ((0 : ℕ) : ℕ∞) < f τ} = Nat.card S :=
    Nat.card_congr (Equiv.subtypeUnivEquiv (by simpa using hf))
  have hsplit : ∑ i ∈ Finset.range (n + 1), (Nat.card {τ : S // (i : ℕ∞) < f τ} : ℝ)
      = (∑ i ∈ Finset.Icc 1 n, (Nat.card {τ : S // (i : ℕ∞) < f τ} : ℝ)) + (Nat.card S : ℝ) := by
    rw [Finset.sum_range_succ' (fun i => (Nat.card {τ : S // (i : ℕ∞) < f τ} : ℝ)) n, h0]
    congr 1
    have hIcc : Finset.Icc 1 n = Finset.Ico 1 (n + 1) := by ext i; simp
    rw [hIcc, Finset.sum_Ico_eq_sum_range]
    simp [add_comm]
  rw [phiOf_natCast, hsplit]
  field_simp
  ring

/-- ★**原文 (i) の前半 `φ_G(0) = 0`** の抽象核。 -/
theorem phiOf_zero {S : Type*} [Fintype S] [Nonempty S] (f : S → ℕ∞) (hf : ∀ τ, 0 < f τ) :
    phiOf f 0 = 0 := by
  simpa using phiOf_natCast_of_pos f hf 0

/-! ## §3 抽象核 —— 剰余類上の和(Lemma 6.10 (ii) の骨)

★★純群論。分岐・付値・Galois の語彙が 1 つも出てこない。 -/

/-- ★★★**(ii) の骨** —— 超距離的な `f : G → ℕ∞` と有限部分群 `H` について、
剰余類 `σH` の中で `f` が最大の代表を `σρ` とすると

`Σ_{τ∈H} min{f(στ), r} = Σ_{τ∈H} min{f(τ), min{f(σρ), r}}`。

★段取り: (1) `τ ↦ ρτ` で添字を取り替え、(2) Proposition 6.9 の中核
`ultrametric_mul_eq_min` で `f(σρτ) = min{f(τ), f(σρ)}`、(3) `truncENat_min`。 -/
theorem sum_truncENat_coset {G : Type*} [Group G] {f : G → ℕ∞}
    (hmul : ∀ a b : G, min (f a) (f b) ≤ f (a * b)) (hinv : ∀ a : G, f a⁻¹ = f a)
    (H : Subgroup G) [Fintype H] {σ ρ : G} (hρ : ρ ∈ H)
    (hmax : ∀ τ ∈ H, f (σ * τ) ≤ f (σ * ρ)) (r : ℝ) :
    ∑ τ : H, truncENat (f (σ * (τ : G))) r
      = ∑ τ : H, truncENat (f (τ : G)) (truncENat (f (σ * ρ)) r) := by
  have hre : ∑ τ : H, truncENat (f (σ * (τ : G))) r
      = ∑ τ : H, truncENat (f (σ * ρ * (τ : G))) r := by
    refine (Fintype.sum_equiv (Equiv.mulLeft (⟨ρ, hρ⟩ : H))
      (fun τ : H => truncENat (f (σ * ρ * (τ : G))) r)
      (fun τ : H => truncENat (f (σ * (τ : G))) r) ?_).symm
    intro τ
    simp [Equiv.coe_mulLeft, mul_assoc]
  rw [hre]
  refine Finset.sum_congr rfl fun τ _ => ?_
  have hcos : f (σ * ρ * (τ : G)) = min (f (τ : G)) (f (σ * ρ)) := by
    refine ultrametric_mul_eq_min hmul hinv ?_
    rw [mul_assoc]
    exact hmax _ (H.mul_mem hρ τ.2)
  rw [hcos, truncENat_min]

/-! ## §4 抽象核 —— ファイバーが一定なら `phiOf` は不変

★★これが「`G` 上の和で書いた `φ` は `G/H` 上の `φ` に等しい」の根拠である(純組合せ)。 -/

/-- ★★`g : S → T` のファイバーの個数がどれも `k > 0` なら `φ_{f∘g} = φ_f`。

★分子は `k` 倍、分母 `|S| = k·|T|` も `k` 倍になって相殺する。 -/
theorem phiOf_comp_of_card_fiber {S T : Type*} [Fintype S] [Fintype T] [Nonempty T]
    (g : S → T) (k : ℕ) (hkpos : 0 < k) (hk : ∀ t : T, Nat.card {s : S // g s = t} = k)
    (f : T → ℕ∞) (n : ℝ) : phiOf (fun s => f (g s)) n = phiOf f n := by
  classical
  have hfib : ∀ t : T, (Finset.univ.filter (fun s : S => g s = t)).card = k := fun t => by
    rw [← hk t, Nat.card_eq_fintype_card, Fintype.card_subtype]
  have hsurj : Function.Surjective g := by
    intro t
    have hpos : 0 < (Finset.univ.filter (fun s : S => g s = t)).card := by
      rw [hfib t]; exact hkpos
    obtain ⟨s, hs⟩ := Finset.card_pos.1 hpos
    exact ⟨s, (Finset.mem_filter.1 hs).2⟩
  have himg : (Finset.univ : Finset S).image g = Finset.univ := by
    ext t; simpa using hsurj t
  have hsum : ∀ r : ℝ, ∑ s : S, truncENat (f (g s)) r = (k : ℝ) * ∑ t : T, truncENat (f t) r := by
    intro r
    rw [Finset.sum_comp (fun t => truncENat (f t) r) g, himg, Finset.mul_sum]
    exact Finset.sum_congr rfl fun t _ => by rw [hfib t, nsmul_eq_mul]
  have hcard : (Nat.card S : ℝ) = (k : ℝ) * (Nat.card T : ℝ) := by
    have h1 : Nat.card S = k * Nat.card T := by
      rw [Nat.card_eq_fintype_card, Nat.card_eq_fintype_card, ← Finset.card_univ,
        ← Finset.card_univ,
        Finset.card_eq_sum_card_fiberwise
          (f := g) (s := (Finset.univ : Finset S)) (t := (Finset.univ : Finset T))
          (fun a _ => Finset.mem_univ _)]
      simp [hfib, mul_comm]
    exact_mod_cast congrArg (fun m : ℕ => (m : ℝ)) h1
  have hk0 : (k : ℝ) ≠ 0 := Nat.cast_ne_zero.2 hkpos.ne'
  have hT : (Nat.card T : ℝ) ≠ 0 := Nat.cast_ne_zero.2 Nat.card_pos.ne'
  rw [phiOf, phiOf, hsum, hcard]
  field_simp

/-- ★★**`G` 上の和で書いた `φ` は `G ⧸ H` 上の `φ` に等しい**。

`QuotientGroup.mk` のファイバーはどれも剰余類で、個数は `|H|`。 -/
theorem phiOf_quotient {G : Type*} [Group G] [Fintype G] (H : Subgroup G) [H.Normal]
    [Fintype (G ⧸ H)] (f : G ⧸ H → ℕ∞) (n : ℝ) :
    phiOf (fun σ : G => f (QuotientGroup.mk σ)) n = phiOf f n := by
  refine phiOf_comp_of_card_fiber _ (Nat.card H) Nat.card_pos (fun b => ?_) f n
  refine Nat.card_congr ?_
  refine (Equiv.ofBijective (fun τ : H => (⟨(Quotient.out b) * (τ : G), ?_⟩ :
      {σ : G // (QuotientGroup.mk σ : G ⧸ H) = b})) ?_).symm
  · simp
  · constructor
    · intro a b' hab
      ext
      simpa using congrArg Subtype.val hab
    · rintro ⟨σ, hσ⟩
      refine ⟨⟨(Quotient.out b)⁻¹ * σ, ?_⟩, ?_⟩
      · rw [← QuotientGroup.eq]; simpa using hσ.symm
      · ext; simp

/-! ## §5 具体層 —— `φ_H` / `φ_G` と `phiOf` の橋 -/

/-- Y11 の `herbrandPhi α H` は `phiOf` の実例(`S := ↥H`, `f := i`)である。 -/
theorem herbrandPhi_eq_phiOf {B : Type*} [CommRing B] [IsDomain B] [IsDiscreteValuationRing B]
    {G : Type*} [Group G] [MulSemiringAction G B] (α : B) (H : Subgroup G) [Fintype H] (n : ℝ) :
    herbrandPhi α H n = phiOf (fun τ : H => ramIndex α (τ : G)) n := by
  rw [herbrandPhi, phiOf, herbrandSum, Nat.card_eq_fintype_card]

/-- ★★**原典の `φ_G`** —— 群 `G` 全体の上の Herbrand 関数
`φ_G(n) := −1 + (1/|G|) Σ_{τ∈G} min{i(τ), n+1}`。

★`Fintype ↥(⊤ : Subgroup G)` は `[Fintype G]` からは推論されない(実測)ので、
`herbrandPhi α (⊤ : Subgroup G)` ではなくこちらを定義した。
両者が一致することは `herbrandPhiGroup_eq_herbrandPhi_top` にある。 -/
noncomputable def herbrandPhiGroup {B : Type*} [CommRing B] [IsDomain B]
    [IsDiscreteValuationRing B] (G : Type*) [Group G] [MulSemiringAction G B] [Fintype G]
    (α : B) (n : ℝ) : ℝ :=
  phiOf (fun σ : G => ramIndex α σ) n

/-- `φ_G` は `H = G`(すなわち `⊤`)のときの `φ_H` に一致する。 -/
theorem herbrandPhiGroup_eq_herbrandPhi_top {B : Type*} [CommRing B] [IsDomain B]
    [IsDiscreteValuationRing B] {G : Type*} [Group G] [MulSemiringAction G B] [Fintype G]
    [Fintype ↥(⊤ : Subgroup G)] (α : B) (n : ℝ) :
    herbrandPhiGroup G α n = herbrandPhi α (⊤ : Subgroup G) n := by
  have hs : ∑ σ : G, truncENat (ramIndex α σ) (n + 1)
      = ∑ τ : ↥(⊤ : Subgroup G), truncENat (ramIndex α (τ : G)) (n + 1) :=
    (Fintype.sum_equiv Subgroup.topEquiv.toEquiv _ _ (fun _ => rfl)).symm
  have hc : (Nat.card G : ℝ) = (Nat.card ↥(⊤ : Subgroup G) : ℝ) := by
    norm_cast
    exact (Nat.card_congr Subgroup.topEquiv.toEquiv).symm
  rw [herbrandPhiGroup, herbrandPhi_eq_phiOf, phiOf, phiOf, hs, hc]

/-- ★★**退化検査の要** `i(σ) ≥ 1` —— `σ•α − α ∈ 𝔪_B` だから。

★これを落とすと Lemma 6.10 (i) の右辺の添字が `Σ_{i=0}` にずれる。 -/
theorem pos_ramIndex {B : Type*} [CommRing B] [IsDomain B] [IsDiscreteValuationRing B]
    {G : Type*} [Group G] [MulSemiringAction G B] {α : B}
    (huni : maximalIdeal B = Ideal.span {α}) (σ : G) : 0 < ramIndex α σ := by
  have hαm : α ∈ maximalIdeal B := by rw [huni]; exact Ideal.mem_span_singleton_self α
  have h := (mem_pow_maximalIdeal_iff_lt_addVal huni 0 (σ • α - α)).1
    (by simpa using Ideal.sub_mem _ (smul_mem_maximalIdeal σ hαm) hαm)
  simpa [ramIndex] using h

/-- `#{τ ∈ H | i(τ) > i}` は `H ∩ G_i = (G_i).subgroupOf H` の位数。 -/
theorem card_subtype_lt_ramIndex {A B : Type*} [CommRing A] [CommRing B] [Algebra A B]
    [IsDomain B] [IsDiscreteValuationRing B] {G : Type*} [Group G] [MulSemiringAction G B]
    [SMulCommClass G A B] {α : B} (huni : maximalIdeal B = Ideal.span {α})
    (hadj : Algebra.adjoin A ({α} : Set B) = ⊤) (H : Subgroup G) (i : ℕ) :
    Nat.card {τ : H // (i : ℕ∞) < ramIndex α (τ : G)}
      = Nat.card ((lowerRamificationGroup B G i).subgroupOf H) :=
  Nat.card_congr (Equiv.subtypeEquivRight fun τ => by
    rw [Subgroup.mem_subgroupOf,
      mem_lowerRamificationGroup_iff_lt_ramIndex huni hadj i (τ : G)]).symm

/-- `#{σ ∈ G | i(σ) > i}` は `G_i` の位数。 -/
theorem card_subtype_lt_ramIndex_top {A B : Type*} [CommRing A] [CommRing B] [Algebra A B]
    [IsDomain B] [IsDiscreteValuationRing B] {G : Type*} [Group G] [MulSemiringAction G B]
    [SMulCommClass G A B] {α : B} (huni : maximalIdeal B = Ideal.span {α})
    (hadj : Algebra.adjoin A ({α} : Set B) = ⊤) (i : ℕ) :
    Nat.card {σ : G // (i : ℕ∞) < ramIndex α σ} = Nat.card (lowerRamificationGroup B G i) :=
  Nat.card_congr (Equiv.subtypeEquivRight fun σ => by
    rw [mem_lowerRamificationGroup_iff_lt_ramIndex huni hadj i σ]).symm

/-! ## §6 Lemma 6.10 (i) -/

def herbrandPhiGroup_natCast.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 16, item := "Lemma 6.10", sectionId := "lemma-6-10" }

/-- ★★★★**Yoshida 2008 Lemma 6.10 (i)** —— `φ_G(n) = (1/|G|) Σ_{i=1}^{n} |G_i|`。

原文 (Yoshida08 p.16):
> (i) φ_G(0) = 0, φ_G(n) = 1|G| ^n_i=1 |G_i| for n ∈ Z[bb]_≥1.

★逸脱: 原典は `n ∈ ℤ≥1` だが、ここは `n : ℕ`(`n = 0` を含む)。
`n = 0` では `Finset.Icc 1 0 = ∅` で右辺 `= 0` となり、(i) の前半 `φ_G(0) = 0` に一致する
——★原典の 2 つの主張が 1 本に収まっている。

★段取り: 抽象核 `phiOf_natCast_of_pos`(純組合せ)に `f := i` を代入するだけ。
正値性 `i(σ) ≥ 1` は `pos_ramIndex`。 -/
theorem herbrandPhiGroup_natCast {A B : Type*} [CommRing A] [CommRing B] [Algebra A B]
    [IsDomain B] [IsDiscreteValuationRing B] {G : Type*} [Group G] [MulSemiringAction G B]
    [Fintype G] [SMulCommClass G A B] {α : B} (huni : maximalIdeal B = Ideal.span {α})
    (hadj : Algebra.adjoin A ({α} : Set B) = ⊤) (n : ℕ) :
    herbrandPhiGroup G α (n : ℝ)
      = (∑ i ∈ Finset.Icc 1 n, (Nat.card (lowerRamificationGroup B G i) : ℝ))
          / (Nat.card G : ℝ) := by
  rw [herbrandPhiGroup, phiOf_natCast_of_pos _ (fun σ => pos_ramIndex huni σ) n]
  congr 1
  exact Finset.sum_congr rfl fun i _ => by
    rw [card_subtype_lt_ramIndex_top (A := A) huni hadj i]

/-- ★★★**Yoshida 2008 Lemma 6.10 (i) の前半** `φ_G(0) = 0`。

原文 (Yoshida08 p.16):
> (i) φ_G(0) = 0
-/
theorem herbrandPhiGroup_zero {B : Type*} [CommRing B] [IsDomain B] [IsDiscreteValuationRing B]
    {G : Type*} [Group G] [MulSemiringAction G B] [Fintype G] {α : B}
    (huni : maximalIdeal B = Ideal.span {α}) : herbrandPhiGroup G α 0 = 0 :=
  phiOf_zero _ (fun σ => pos_ramIndex huni σ)

/-- ★**(i) の一般の部分群版** `φ_H(n) = (1/|H|) Σ_{i=1}^{n} |H ∩ G_i|`。

★下流(Hasse-Arf)は `φ_H` の側でもこの形を要求する。 -/
theorem herbrandPhi_natCast {A B : Type*} [CommRing A] [CommRing B] [Algebra A B]
    [IsDomain B] [IsDiscreteValuationRing B] {G : Type*} [Group G] [MulSemiringAction G B]
    [SMulCommClass G A B] {α : B} (huni : maximalIdeal B = Ideal.span {α})
    (hadj : Algebra.adjoin A ({α} : Set B) = ⊤) (H : Subgroup G) [Fintype H] (n : ℕ) :
    herbrandPhi α H (n : ℝ)
      = (∑ i ∈ Finset.Icc 1 n,
          (Nat.card ((lowerRamificationGroup B G i).subgroupOf H) : ℝ)) / (Nat.card H : ℝ) := by
  rw [herbrandPhi_eq_phiOf, phiOf_natCast_of_pos _ (fun τ : H => pos_ramIndex huni (τ : G)) n]
  congr 1
  exact Finset.sum_congr rfl fun i _ => by
    rw [card_subtype_lt_ramIndex (A := A) huni hadj H i]

/-- ★**(i) の前半の一般の部分群版** `φ_H(0) = 0`。 -/
theorem herbrandPhi_zero {B : Type*} [CommRing B] [IsDomain B] [IsDiscreteValuationRing B]
    {G : Type*} [Group G] [MulSemiringAction G B] {α : B}
    (huni : maximalIdeal B = Ideal.span {α}) (H : Subgroup G) [Fintype H] :
    herbrandPhi α H 0 = 0 := by
  rw [herbrandPhi_eq_phiOf]
  exact phiOf_zero _ (fun τ : H => pos_ramIndex huni (τ : G))

/-! ## §7 具体層 —— 剰余類ごとの等式((ii) の心臓) -/

/-- `H` が `C` に自明に作用するので `i_ϖ` は `H` 剰余類の上で定数。 -/
theorem ramIndex_mul_mem_eq {G : Type*} [Group G] {H : Subgroup G}
    {C : Type*} [CommRing C] [IsDomain C] [IsDiscreteValuationRing C] [MulSemiringAction G C]
    (hHtriv : ∀ ρ ∈ H, ∀ c : C, ρ • c = c) (ϖ : C) (σ : G) {τ : G} (hτ : τ ∈ H) :
    ramIndex ϖ (σ * τ) = ramIndex ϖ σ := by
  show addVal C ((σ * τ) • ϖ - ϖ) = addVal C (σ • ϖ - ϖ)
  rw [mul_smul, hHtriv τ hτ ϖ]

/-- ★★★★**Lemma 6.10 (ii) の心臓** —— 剰余類 `σH` ごとの等式

`Σ_{τ∈H} min{i(στ), n+1} = |H| · min{i_ϖ(σ), φ_H(n)+1}`。

★★これが原文の「両辺の微分が一致する」の代わりである。左辺・右辺とも
`Σ_{τ∈H} min{i(τ), min{m, n+1}}`(`m := i(σρ)`、`σρ` は剰余類の中の `i` 最大の代表)に等しい。

★段取り:
1. `exists_max_on_coset`(Y11)で最大代表 `σρ` を取る。
2. `sum_truncENat_coset`(§3, 純群論)で左辺 `= herbrandSum π' H (min{m, n+1} − 1)`。
3. `m = ⊤`(すなわち `σ ∈ H`)なら両辺とも `|H|·(φ_H(n)+1)`。
4. `m` が有限なら Lemma 6.8(`card_mul_ramIndex_eq_sum_ramIndex_of_fixedRing`)で
   `i_ϖ(σ) = φ_H(m−1)+1`、`φ_H` の単調性(`Monotone.map_min`)で
   `min{φ_H(m−1)+1, φ_H(n)+1} = φ_H(min{m−1, n})+1`。

★仮定はすべて Y10(Lemma 6.8)・Y11(Proposition 6.9)のものをそのまま引き継いだ。
本補題が新たに足した仮定は 1 つも無い。 -/
theorem coset_sum_truncENat {A B : Type*} [CommRing A] [CommRing B] [Algebra A B] [IsDomain B]
    [IsDiscreteValuationRing B] {G : Type*} [Group G] [MulSemiringAction G B]
    [SMulCommClass G A B] [FaithfulSMul G B] {H : Subgroup G} [Fintype H] {π' π'' : B}
    {C : Type*} [CommRing C] [IsDomain C] [IsDiscreteValuationRing C] [Algebra C B]
    [MulSemiringAction G C]
    (hcomp : ∀ (ρ : G) (c : C), algebraMap C B (ρ • c) = ρ • algebraMap C B c)
    (hHtriv : ∀ ρ ∈ H, ∀ c : C, ρ • c = c)
    (hπ' : Irreducible π')
    (hinj : Function.Injective (algebraMap C B))
    (hfixC : ∀ b : B, (∀ ρ ∈ H, ρ • b = b) → ∃ c : C, algebraMap C B c = b)
    (hres : ∀ b : B, ∃ c : C, b - algebraMap C B c ∈ maximalIdeal B)
    (hAC : ∀ a : A, ∃ c : C, algebraMap C B c = algebraMap A B a)
    {ϖ : C} (hϖ : algebraMap C B ϖ = π'')
    (hadj : Algebra.adjoin A ({π'} : Set B) = ⊤)
    (hfix : ∀ y : B, (∀ ρ ∈ H, ρ • y = y) → y ∈ Algebra.adjoin A ({π''} : Set B))
    (n : ℝ) (σ : G) :
    ∑ τ : H, truncENat (ramIndex π' (σ * (τ : G))) (n + 1)
      = (Nat.card H : ℝ) * truncENat (ramIndex ϖ σ) (herbrandPhi π' H n + 1) := by
  classical
  have hcardpos : 0 < Nat.card H := Nat.card_pos
  have hcardR : (0 : ℝ) < (Nat.card H : ℝ) := by exact_mod_cast hcardpos
  obtain ⟨ρ, hρ, hmax⟩ := exists_max_on_coset (fun s : G => ramIndex π' s) H σ
  have hcoset : ∀ τ : H, ramIndex π' (σ * ρ * (τ : G))
      = min (ramIndex π' (τ : G)) (ramIndex π' (σ * ρ)) := by
    intro τ
    refine ultrametric_mul_eq_min (min_ramIndex_le_ramIndex_mul π') (ramIndex_inv π') ?_
    rw [mul_assoc]
    exact hmax _ (H.mul_mem hρ τ.2)
  have hϖeq : ramIndex ϖ (σ * ρ) = ramIndex ϖ σ := ramIndex_mul_mem_eq hHtriv ϖ σ hρ
  have hL : ∑ τ : H, truncENat (ramIndex π' (σ * (τ : G))) (n + 1)
      = herbrandSum π' H (truncENat (ramIndex π' (σ * ρ)) (n + 1) - 1) := by
    rw [sum_truncENat_coset (min_ramIndex_le_ramIndex_mul π') (ramIndex_inv π') H hρ hmax (n + 1)]
    refine Finset.sum_congr rfl fun τ _ => ?_
    congr 1
    ring
  have hkey : truncENat (ramIndex ϖ σ) (herbrandPhi π' H n + 1)
      = herbrandPhi π' H (truncENat (ramIndex π' (σ * ρ)) (n + 1) - 1) + 1 := by
    rcases eq_or_ne (ramIndex π' (σ * ρ)) ⊤ with hm | hm
    · -- `m = ⊤`、すなわち `σρ = 1`(`σ ∈ H`)。
      have h1 : σ * ρ = 1 :=
        eq_one_of_smul_eq_of_adjoin_eq_top hadj ((ramIndex_eq_top_iff π' (σ * ρ)).1 hm)
      have h2 : ramIndex ϖ σ = ⊤ := by rw [← hϖeq, h1, ramIndex_one]
      rw [hm, h2, truncENat_top, truncENat_top]
      norm_num
    · -- `m = M : ℕ` が有限の場合。
      have hMeq : ramIndex π' (σ * ρ) = (((ramIndex π' (σ * ρ)).toNat : ℕ) : ℕ∞) :=
        (ENat.coe_toNat hm).symm
      set M : ℕ := (ramIndex π' (σ * ρ)).toNat with hMdef
      have h68 : (Nat.card H : ℕ∞) * ramIndex ϖ σ
          = ∑ τ : H, min (ramIndex π' (τ : G)) (ramIndex π' (σ * ρ)) := by
        rw [← hϖeq, card_mul_ramIndex_eq_sum_ramIndex_of_fixedRing (A := A) hcomp hHtriv hπ' hinj
          hfixC hres hAC hϖ hadj hfix (σ * ρ)]
        exact Finset.sum_congr rfl fun τ _ => hcoset τ
      have hsum : ∑ τ : H, min (ramIndex π' (τ : G)) (ramIndex π' (σ * ρ))
          = ((∑ τ : H, (min (ramIndex π' (τ : G)) (M : ℕ∞)).toNat : ℕ) : ℕ∞) := by
        rw [hMeq, Nat.cast_sum]
        exact Finset.sum_congr rfl fun τ _ => (ENat.coe_toNat (min_natCast_ne_top _ _)).symm
      rw [hsum] at h68
      have hJne : ramIndex ϖ σ ≠ ⊤ := by
        intro hc
        rw [hc, ENat.mul_top (by exact_mod_cast hcardpos.ne')] at h68
        exact (ENat.coe_ne_top _) h68.symm
      have hJeq : ramIndex ϖ σ = (((ramIndex ϖ σ).toNat : ℕ) : ℕ∞) := (ENat.coe_toNat hJne).symm
      set J : ℕ := (ramIndex ϖ σ).toNat with hJdef
      have hnat : Nat.card H * J = ∑ τ : H, (min (ramIndex π' (τ : G)) (M : ℕ∞)).toNat := by
        rw [hJeq] at h68
        exact_mod_cast h68
      have hreal : (Nat.card H : ℝ) * (J : ℝ) = herbrandSum π' H ((M : ℝ) - 1) := by
        have h1 : ((Nat.card H * J : ℕ) : ℝ)
            = ((∑ τ : H, (min (ramIndex π' (τ : G)) (M : ℕ∞)).toNat : ℕ) : ℝ) := by
          exact_mod_cast hnat
        rw [Nat.cast_mul, Nat.cast_sum] at h1
        rw [h1]
        unfold herbrandSum
        refine Finset.sum_congr rfl fun τ _ => ?_
        rw [show (M : ℝ) - 1 + 1 = (M : ℝ) by ring]
        exact toNat_min_natCast _ _
      -- ★原文「the Lemma 6.8 gives i(σ[bar]) = φH(m − 1) + 1」
      have hJval : (J : ℝ) = herbrandPhi π' H ((M : ℝ) - 1) + 1 := by
        rw [← card_mul_herbrandPhi_add_one π' H ((M : ℝ) - 1)] at hreal
        exact mul_left_cancel₀ hcardR.ne' hreal
      rw [hJeq, hMeq, truncENat_coe, truncENat_coe, hJval,
        show min (M : ℝ) (n + 1) - 1 = min ((M : ℝ) - 1) n by
          rw [← min_sub_sub_right]; norm_num,
        (strictMono_herbrandPhi π' H).monotone.map_min, min_add_add_right]
  rw [hL, hkey, ← card_mul_herbrandPhi_add_one]

/-! ## §8 Lemma 6.10 (ii) -/

def herbrandPhiGroup_comp.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 16, item := "Lemma 6.10", sectionId := "lemma-6-10" }

/-- ★★★★**Yoshida 2008 Lemma 6.10 (ii)** —— `φ_G = φ_{G/H} ∘ φ_H`。

原文 (Yoshida08 p.16):
> (ii) φ_G = φ_G/H ◦ φ_H on R[bb]_≥0.

`φ_{G/H}` は `herbrandPhiGroup G ϖ`(`G` 上の和、分母 `|G|`)として書いてある。
これが商群の上の `φ` と同じ値であることは `herbrandPhiGroup_eq_phiOf_quotient` にある。

★★**原文の証明(連続・区分線形・微分の比較)は使っていない。** 掛け算形にすると
両辺は `G` 上の和になり、`H` 剰余類ごとに `coset_sum_truncENat` で一致する。
★剰余類への分割も要らない: `Σ_{σ∈G} Σ_{τ∈H}` を入れ替え、`σ ↦ στ` が `G` の
全単射であることを使うと `|H| · Σ_{σ∈G}` になる(`hfub`)。

★★**`H` の正規性は仮定していない**(ファイル冒頭「逸脱の記録 (3)」)。 -/
theorem herbrandPhiGroup_comp {A B : Type*} [CommRing A] [CommRing B] [Algebra A B] [IsDomain B]
    [IsDiscreteValuationRing B] {G : Type*} [Group G] [MulSemiringAction G B] [Fintype G]
    [SMulCommClass G A B] [FaithfulSMul G B] {H : Subgroup G} [Fintype H] {π' π'' : B}
    {C : Type*} [CommRing C] [IsDomain C] [IsDiscreteValuationRing C] [Algebra C B]
    [MulSemiringAction G C]
    (hcomp : ∀ (ρ : G) (c : C), algebraMap C B (ρ • c) = ρ • algebraMap C B c)
    (hHtriv : ∀ ρ ∈ H, ∀ c : C, ρ • c = c)
    (hπ' : Irreducible π')
    (hinj : Function.Injective (algebraMap C B))
    (hfixC : ∀ b : B, (∀ ρ ∈ H, ρ • b = b) → ∃ c : C, algebraMap C B c = b)
    (hres : ∀ b : B, ∃ c : C, b - algebraMap C B c ∈ maximalIdeal B)
    (hAC : ∀ a : A, ∃ c : C, algebraMap C B c = algebraMap A B a)
    {ϖ : C} (hϖ : algebraMap C B ϖ = π'')
    (hadj : Algebra.adjoin A ({π'} : Set B) = ⊤)
    (hfix : ∀ y : B, (∀ ρ ∈ H, ρ • y = y) → y ∈ Algebra.adjoin A ({π''} : Set B))
    (n : ℝ) :
    herbrandPhiGroup G π' n = herbrandPhiGroup G ϖ (herbrandPhi π' H n) := by
  have hcardpos : 0 < Nat.card H := Nat.card_pos
  have hcardR : (0 : ℝ) < (Nat.card H : ℝ) := by exact_mod_cast hcardpos
  have hshift : ∀ τ : G, ∑ σ : G, truncENat (ramIndex π' (σ * τ)) (n + 1)
      = ∑ σ : G, truncENat (ramIndex π' σ) (n + 1) := fun τ =>
    Fintype.sum_equiv (Equiv.mulRight τ) _ _ (fun _ => rfl)
  have hfub : ∑ σ : G, ∑ τ : H, truncENat (ramIndex π' (σ * (τ : G))) (n + 1)
      = (Nat.card H : ℝ) * ∑ σ : G, truncENat (ramIndex π' σ) (n + 1) := by
    rw [Finset.sum_comm, Finset.sum_congr rfl (fun (τ : H) _ => hshift (τ : G)),
      Finset.sum_const, nsmul_eq_mul, Nat.card_eq_fintype_card, Finset.card_univ]
  have hcos : ∑ σ : G, ∑ τ : H, truncENat (ramIndex π' (σ * (τ : G))) (n + 1)
      = (Nat.card H : ℝ) * ∑ σ : G, truncENat (ramIndex ϖ σ) (herbrandPhi π' H n + 1) := by
    rw [Finset.mul_sum]
    exact Finset.sum_congr rfl fun σ _ =>
      coset_sum_truncENat (A := A) hcomp hHtriv hπ' hinj hfixC hres hAC hϖ hadj hfix n σ
  have hEq : ∑ σ : G, truncENat (ramIndex π' σ) (n + 1)
      = ∑ σ : G, truncENat (ramIndex ϖ σ) (herbrandPhi π' H n + 1) :=
    mul_left_cancel₀ hcardR.ne' (by rw [← hfub, hcos])
  rw [herbrandPhiGroup, herbrandPhiGroup, phiOf, phiOf, hEq]

/-- `i_ϖ` は `G ⧸ H` を経由する(`H` は `C` に自明に作用するから)。 -/
theorem exists_quotient_ramIndex {G : Type*} [Group G] {H : Subgroup G} [H.Normal]
    {C : Type*} [CommRing C] [IsDomain C] [IsDiscreteValuationRing C] [MulSemiringAction G C]
    (hHtriv : ∀ ρ ∈ H, ∀ c : C, ρ • c = c) (ϖ : C) :
    ∃ f : G ⧸ H → ℕ∞, ∀ σ : G, f (QuotientGroup.mk σ) = ramIndex ϖ σ := by
  refine ⟨Quotient.lift (fun σ : G => ramIndex ϖ σ) ?_, fun _ => rfl⟩
  intro a b hab
  have hmem : a⁻¹ * b ∈ H := (QuotientGroup.leftRel_apply).1 hab
  have h := ramIndex_mul_mem_eq hHtriv ϖ a hmem
  rw [mul_inv_cancel_left] at h
  exact h.symm

/-- ★★**`herbrandPhiGroup G ϖ` は本当に `φ_{G/H}` である**。

`i_ϖ` の `G ⧸ H` への降下 `f`(`exists_quotient_ramIndex` で存在する)に対して
`herbrandPhiGroup G ϖ n = phiOf f n`。★根拠は `|G| = |H|·|G/H|` と
「`i_ϖ` は `H` 剰余類上で定数」の 2 つだけ(§4 の `phiOf_quotient`)。 -/
theorem herbrandPhiGroup_eq_phiOf_quotient {G : Type*} [Group G] [Fintype G] {H : Subgroup G}
    [H.Normal] [Fintype (G ⧸ H)]
    {C : Type*} [CommRing C] [IsDomain C] [IsDiscreteValuationRing C] [MulSemiringAction G C]
    (ϖ : C) (f : G ⧸ H → ℕ∞) (hf : ∀ σ : G, f (QuotientGroup.mk σ) = ramIndex ϖ σ) (n : ℝ) :
    herbrandPhiGroup G ϖ n = phiOf f n := by
  rw [herbrandPhiGroup, ← phiOf_quotient H f n]
  simp only [hf]

end ABC3.Found.PGC
