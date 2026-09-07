import ABC3.Found.PGC.RamificationJumpDivisibility

/-!
# 可換な `G` の跳びの割り切り —— Corollary 6.3 の原文どおりの証明と、Hasse-Arf 最終段への受け渡し

典拠: T. Yoshida, *Local Class Field Theory via Lubin-Tate Theory*
(Ann. Fac. Sci. Toulouse 17-2, 2008; arXiv math/0606108) **Corollary 6.3**(物理 p.14)。
構造化済み原文は `ResearchPaper/1_Structured/Local Class Field Theory via Lubin-Tate Theory/`
`section-6.html` の `id="cor-6-3"`。

原文 (Yoshida08 p.14):
> Corollary 6.3. If G is abelian and G_n = G_n+1, then e_0 := |G_0/G_1| divides n.

★★**上の逐語は `pdftotext` に見えている形であり、原典ではない。**
原典は `G_n ≠ G_{n+1}` である。`≠` の斜線はベクター描画なので `pdftotext` が落とし、
出力は `=` だけになる。**テキストだけを読むと主張が反転する。**
構造化済み HTML は `data-txt="="` でこの脱落を明示している。
本ファイルの形式化はすべて `≠`(`hne : G_n ≠ G_{n+1}`)側で書いてある。

## ★★本ファイルの位置づけ(先行ノードとの関係)

**Corollary 6.3 そのものは既に `Found/PGC/RamificationJumpDivisibility.lean` で埋まっている**
(`nat_card_quot_zero_dvd_of_ne` / `nat_card_quot_top_dvd_of_ne` /
`nat_card_quot_lowerRamificationGroupAdjoin_one_dvd`)。本ファイルはその上に 2 つを足す。

**(A) 原文の段取りに忠実な別証明**(`nat_card_quot_zero_dvd_of_notMem_generator`)。
先行ノードは「`G_0/G_1` の全ての元が `n` 乗で `1`」→「指数 = 位数」(`Monoid.exponent`)
という経路を取った。**原文は違う経路を取っている**:

> choose σ ∈ G which generates G_0/G_1, i.e. θ_0(σ) = u mod p has order e_0
> in (O_{K′}/p_{K′})^×. Then a ≡ u^n a (mod p^{n+1}_{K′}) implies e_0 | n.

すなわち原文は **`u` の位数そのもの**を見ている。ここではその形
(`orderOf_thetaMul_dvd_of_notMem` + `exists_orderOf_thetaMul_eq_natCard_quot`)で書いた。
★2 つの経路は同じ結論を出すが、**中間の情報が違う**: 原文の経路は
「任意の `τ ∈ G_0` について `orderOf θ_0(τ) ∣ n`」という、生成元でない `τ` にも効く
形を経由する。

**(B) Hasse-Arf(Theorem 6.11)の最終段が Corollary 6.3 から取り出す形**。原文:

> Now when G ≠ G_1, set H = G_1 and |G/H| = e_0. ... it suffices to show e_0 | φ_H(n)
> when n ∈ Z≥0 and G_n ≠ G_{n+1}. ... For any i ∈ Z≥1 (where H_i = G_i) with H_i ≠ H_{i+1},
> we have e_0 | i by Corollary 6.3, hence e_0 | Σ_{i=1}^n |H_i|. As e_0 and |H| are coprime,
> we have e_0 | φ_H(n) by Lemma 6.10(i).

★原文が **"hence"** の 1 語で畳んでいるのは「跳びが `e_0` の倍数でしか起きない列の
`[1, n]` 上の和は `e_0` で割れる」という**純粋な算術**である(§1 の `dvd_sum_Icc_of_step`)。
本ファイルはこれを切り出し、`dvd_sum_natCard_lowerRamificationGroup` として
`e_0 ∣ Σ_{i=1}^n |G_i|` を出す。さらに "As e_0 and |H| are coprime" を
`coprime_natCard_quot_natCard_lowerRamificationGroup_one` として証明し、
最後の一歩 `dvd_of_sum_eq_natCard_mul` は **Lemma 6.10(i) を仮定の形で受け取る**
(`Σ_{i=1}^n |G_i| = |G_1| · k` なら `e_0 ∣ k`)。
★★**Lemma 6.10 を import しない**のはわざとである(同時並行で形式化中のため、
依存を作らずに「差し込むだけ」の形にしてある)。

## 抽象核と具体層(設計)

★**分岐・付値・Galois の語彙が 1 つも出てこない**ものを §1 に集めた。§2–§4 はそこへの代入である。

| 抽象核(§1) | 内容 |
|---|---|
| `orderOf_dvd_of_pow_mul_eq_self` | `u^n a = a`、`a ≠ 0` ⟹ `orderOf u ∣ n`(原文の最後の一歩) |
| `not_dvd_natCard_of_injective_units` | 標数 `p` の整域の単数群に埋め込める有限群の位数は `p` と素 |
| `exists_orderOf_map_eq_natCard` | 有限巡回群の生成元は像の位数が `|Q|`(原文「generates G_0/G_1」) |
| `exists_mem_notMem_of_antitone_ne` | 反単調な部分群の列で `F n ≠ F (n+1)` なら差に元がある |
| `eq_of_forall_eq_succ` | `f i = f (i+1)` が `[a,b)` で成り立てば `f b = f a` |
| `sum_Ico_block_eq` / `dvd_sum_Ico_mul_of_step` / `dvd_sum_Icc_of_step` | 原文の "hence" の中身(ブロック和) |

## ★逸脱の記録

1. **`G` の可換性は `∀ x y : G, x * y = y * x` という命題で渡す**(`CommGroup` 構造を
   要求しない)。先行ノード `RamificationJumpDivisibility.lean` の逸脱 2 と同じ理由:
   具体層の `G = Gal(K(x)/K)` は `Group` インスタンスしか持たない。数学的には同じ。
2. **`e_0` は `Nat.card` で表す**。原文の `e_0 := |G_0/G_1|` をそのまま
   `Nat.card (G_0 ⧸ (G_1).subgroupOf G_0)`(完全分岐版は `Nat.card (G ⧸ G_1)`)とした。
   ★「生成元 `σ` の `θ_0(σ)` の位数」という原文のもう 1 つの表し方も
   `exists_orderOf_thetaMul_eq_natCard_quot` で出してあり、**両者が一致することを証明済み**。
3. **`Σ_{i=1}^n |H_i|` は `Finset.Icc 1 n` 上の和で書く**。原文の `H_i = G_i`(`i ≥ 1`)は
   `subgroupOf_lowerRamificationGroup` により添字を `G` 側に寄せてよいので、
   `Nat.card (lowerRamificationGroup B G i)` をそのまま足している。
4. **最後の一歩(`e_0 ∣ φ_H(n)`)は Lemma 6.10(i) を仮定として受け取る**(上記 (B))。
   `φ_H` を持ち込むと Lemma 6.10 への依存が生じるため、`Σ = |G_1| · k` という
   等式を仮定に置き、`e_0 ∣ k` を結論とした。★これは弱化ではない: Lemma 6.10(i) が
   出来た時点で `k := φ_H(n)` を代入すれば原文の主張になる。
5. **`dvd_sum_Icc_of_step` に `0 < e` を仮定していない**。`e = 0` のときブロックが空で
   両辺 `0` になり、そのまま真になる(実測で不要と分かったので落とした)。

## ★退化の自己検査

* ★★**`G` 可換を落とすと偽**。使いどころは `orderOf_thetaMul_dvd_of_notMem` の
  `hconj : τ σ τ⁻¹ = σ` **ただ 1 箇所**である。落とすと残るのは共役の変換則
  `residue_ramCoeff_conj`(= `θ_n(τστ⁻¹) = θ_0(τ)^n θ_n(σ)`)だけで、
  `θ_0(τ)^n = 1` は出ない。
* ★★**`G_n ≠ G_{n+1}` を落とすと空虚**。`θ_n(τ) ≠ 0` なる `τ` が取れず、
  `u^n a = a` から `a` を約せない。★**原典のテキストは `=` にしか見えない**ので、
  ここを逆に書く危険が実際にある(ファイル冒頭の警告)。
* ★**`σ` が `G_0/G_1` を生成することを落とすと `e_0` ではなく `u` の位数しか出ない**。
  本ファイルはその区別を `orderOf_thetaMul_dvd_of_notMem`(任意の `τ`、`u` の位数)と
  `exists_orderOf_thetaMul_eq_natCard_quot`(生成元、`e_0`)の 2 宣言に分けてある。
* ★**`a ≠ 0`(`a ∈ 𝔭^n \ 𝔭^{n+1}`)を落とすと `u^n a = a` から何も出ない**。
  `orderOf_dvd_of_pow_mul_eq_self` の `ha : a ≠ 0` がそれで、
  整域でなければ約分できないので `[IsDomain R]` も落とせない。
* ★**巡回性を落とすと (A) が壊れる**。`(ℤ/2)²` と `n = 2` で `|Q| = 4 ∤ 2`。
  巡回性は `θ_0 : G_0/G_1 ↪ 𝓀ˣ` の単射性(Lemma 5.11)から来ている。
* ★**(B) の互いに素は `[CharP (ResidueField B) p]` と `Fact p.Prime` の両方が要る**。
  `e_0` が `p` と素なのは `𝓀ˣ` が標数 `p` の体の単数群だからであり(`x^p = 1 ⟹ x = 1`)、
  `|G_1|` が `p` 冪なのは `G_1` が `p` 群だから(`isPGroup_lowerRamificationGroup_one`)。
  どちらか一方でも落とすと互いに素は言えない。
* ★**`dvd_sum_Icc_of_step` は `e ∣ n` を落とすと偽**。`e = 2`、`f = id`、`n = 1` で
  和は `1`。原文が `e_0 | i` から `e_0 | Σ` へ渡るときに、**`n` 自身が `e_0` の倍数**
  であること(= Corollary 6.3 を `n` に適用したもの)を暗黙に使っている。

## 実測

抽象核(§1)は `lean_check` 0.01–0.30 秒、Corollary 6.3 の別証明(§2)は 0.24 秒、
Hasse-Arf 受け渡し(§3)は 0.12–0.13 秒、`PAdicLocalField` 具体層(§4)は 0.42–0.45 秒。
★見積 250–500 行に対し、証明本体は約 90 行。理由は
`RamificationJumpDivisibility.lean` の `residue_ramCoeff_conj` が原文の
`θ_n(στσ⁻¹) = u^n θ_n(τ)` をそのまま与えており、`≠` の側の仮定処理も
`nat_card_quot_top_dvd_of_ne` で済んだこと。
-/

namespace ABC3.Found.PGC

open IsLocalRing ABC3.Skeleton.PGC
open scoped NNReal Valued

def dvd_sum_natCard_lowerRamificationGroup.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 14, item := "Corollary 6.3", sectionId := "cor-6-3" }

/-! ## §1 抽象核

★**分岐も付値も Galois も出てこない**。§2 以降はここへの代入だけである。 -/

/-- ★★**原文の最後の一歩**「`a ≡ u^n a (mod 𝔭^{n+1})` implies `e_0 | n`」の核。

整域 `R` の単数 `u` と `a ≠ 0` について `u^n a = a` なら `orderOf u ∣ n`。
★`a ≠ 0` と `[IsDomain R]` の両方が要る(片方でも落とすと約分できない)。 -/
theorem orderOf_dvd_of_pow_mul_eq_self {R : Type*} [CommRing R] [IsDomain R] (u : Rˣ) {a : R}
    (ha : a ≠ 0) {n : ℕ} (h : ((u : R)) ^ n * a = a) : orderOf u ∣ n := by
  refine orderOf_dvd_of_pow_eq_one (Units.ext ?_)
  push_cast
  exact mul_right_cancel₀ ha (h.trans (one_mul a).symm)

/-- ★**標数 `p` の整域の単数群に単射に埋め込める有限群の位数は `p` と素**。

原文「As e_0 and |H| are coprime」の `e_0` 側。証明は Cauchy(`exists_prime_orderOf_dvd_card`)と
Frobenius(`sub_pow_char`): 位数 `p` の元 `q` があれば `(f q − 1)^p = (f q)^p − 1 = 0` で
`f q = 1`、単射性から `q = 1` となり位数 `1 ≠ p`。 -/
theorem not_dvd_natCard_of_injective_units {Q R : Type*} [Group Q] [Finite Q] [CommRing R]
    [IsDomain R] (p : ℕ) [Fact p.Prime] [CharP R p] (f : Q →* Rˣ) (hf : Function.Injective f) :
    ¬ p ∣ Nat.card Q := by
  intro hp
  haveI : Fintype Q := Fintype.ofFinite Q
  rw [Nat.card_eq_fintype_card] at hp
  obtain ⟨q, hq⟩ := exists_prime_orderOf_dvd_card (G := Q) p hp
  have h1 : (f q) ^ p = 1 := by
    rw [← map_pow, ← hq, pow_orderOf_eq_one, map_one]
  have h2 : ((f q : R) - 1) ^ p = 0 := by
    haveI : ExpChar R p := ExpChar.prime (Fact.out)
    rw [sub_pow_char, one_pow, ← Units.val_pow_eq_pow_val, h1, Units.val_one, sub_self]
  have h3 : (f q : R) = 1 := by
    have := (pow_eq_zero_iff (n := p) (Nat.Prime.pos (Fact.out)).ne').mp h2
    linear_combination this
  have h4 : q = 1 := hf (by rw [map_one]; exact Units.ext h3)
  rw [h4, orderOf_one] at hq
  exact (Nat.Prime.one_lt (Fact.out : p.Prime)).ne hq

/-- ★**原文「choose σ ∈ G which generates G_0/G_1, i.e. θ_0(σ) has order e_0」の核**。

有限巡回群 `Q` が `M` に単射に入るなら、**像の位数がちょうど `|Q|` になる元**が取れる。 -/
theorem exists_orderOf_map_eq_natCard {Q M : Type*} [Group Q] [Finite Q] [IsCyclic Q] [Monoid M]
    (f : Q →* M) (hf : Function.Injective f) : ∃ g : Q, orderOf (f g) = Nat.card Q := by
  obtain ⟨g, hg⟩ := IsCyclic.exists_generator (α := Q)
  exact ⟨g, by rw [orderOf_injective f hf g, orderOf_eq_card_of_forall_mem_zpowers hg]⟩

/-- ★**「商が非自明なら核の外に元がある」の純群論版**。
原文「If G_n ≠ G_{n+1}, we can choose τ ∈ G_n with θ_n(τ) ≠ 0」の第 1 段。 -/
theorem exists_mem_notMem_of_antitone_ne {G : Type*} [Group G] {F : ℕ → Subgroup G}
    (hF : Antitone F) {n : ℕ} (hne : F n ≠ F (n + 1)) : ∃ σ ∈ F n, σ ∉ F (n + 1) :=
  SetLike.exists_of_lt (lt_of_le_of_ne (hF (Nat.le_succ n)) (Ne.symm hne))

/-- `f i = f (i+1)` が `[a, b)` の全ての `i` で成り立てば `f b = f a`。 -/
theorem eq_of_forall_eq_succ {M : Type*} (f : ℕ → M) {a b : ℕ} (hab : a ≤ b)
    (h : ∀ i, a ≤ i → i < b → f i = f (i + 1)) : f b = f a := by
  induction b, hab using Nat.le_induction with
  | base => rfl
  | succ b hb ih =>
      rw [← h b hb (Nat.lt_succ_self b)]
      exact ih fun i hi hib => h i hi (hib.trans (Nat.lt_succ_self b))

/-- 跳びが `e` の倍数でしか起きないなら、長さ `e` の 1 ブロックの和は `e * f(始点)`。 -/
theorem sum_Ico_block_eq {e : ℕ} (f : ℕ → ℕ)
    (hstep : ∀ i, 1 ≤ i → ¬ (e ∣ i) → f i = f (i + 1)) (m : ℕ) :
    ∑ i ∈ Finset.Ico (e * m + 1) (e * m + e + 1), f i = e * f (e * m + 1) := by
  have hconst : ∀ j ∈ Finset.Ico (e * m + 1) (e * m + e + 1), f j = f (e * m + 1) := by
    intro j hj
    rw [Finset.mem_Ico] at hj
    refine eq_of_forall_eq_succ f hj.1 ?_
    intro i hi hij
    refine hstep i (by omega) ?_
    rintro ⟨k, rfl⟩
    have h1 : m < k := by
      by_contra hc
      have : e * k ≤ e * m := Nat.mul_le_mul_left e (Nat.not_lt.1 hc)
      omega
    have h2 : k < m + 1 := by
      by_contra hc
      have h3 : e * (m + 1) ≤ e * k := Nat.mul_le_mul_left e (Nat.not_lt.1 hc)
      rw [Nat.mul_succ] at h3
      omega
    omega
  rw [Finset.sum_congr rfl hconst, Finset.sum_const, Nat.card_Ico, smul_eq_mul]
  congr 1
  omega

/-- ブロックを `m` 個つないだ和は `e` で割れる。 -/
theorem dvd_sum_Ico_mul_of_step {e : ℕ} (f : ℕ → ℕ)
    (hstep : ∀ i, 1 ≤ i → ¬ (e ∣ i) → f i = f (i + 1)) (m : ℕ) :
    e ∣ ∑ i ∈ Finset.Ico 1 (e * m + 1), f i := by
  induction m with
  | zero => simp
  | succ m ih =>
      have hle : e * m + 1 ≤ e * (m + 1) + 1 := by rw [Nat.mul_succ]; omega
      rw [← Finset.sum_Ico_consecutive f (Nat.le_add_left 1 (e * m)) hle]
      refine Nat.dvd_add ih ?_
      have hb : e * (m + 1) + 1 = e * m + e + 1 := by rw [Nat.mul_succ]
      rw [hb, sum_Ico_block_eq f hstep m]
      exact Dvd.intro _ rfl

/-- ★★**原文の "hence e_0 | Σ_{i=1}^n |H_i|" の中身**——分岐も付値も出てこない純算術。

`f : ℕ → ℕ` の跳びが `e` の倍数の添字でしか起きないなら、`e ∣ n` のとき
`Σ_{i=1}^n f i` は `e` で割れる。
★`e ∣ n` を落とすと偽(`e = 2`、`f = id`、`n = 1` で和は `1`)。 -/
theorem dvd_sum_Icc_of_step {e : ℕ} (f : ℕ → ℕ)
    (hstep : ∀ i, 1 ≤ i → ¬ (e ∣ i) → f i = f (i + 1)) {n : ℕ} (hn : e ∣ n) :
    e ∣ ∑ i ∈ Finset.Icc 1 n, f i := by
  obtain ⟨m, rfl⟩ := hn
  have h : Finset.Icc 1 (e * m) = Finset.Ico 1 (e * m + 1) := by
    ext x; simp only [Finset.mem_Icc, Finset.mem_Ico]; omega
  rw [h]
  exact dvd_sum_Ico_mul_of_step f hstep m

/-! ## §2 Corollary 6.3 —— 原文の段取りに忠実な証明

先行ノード `RamificationJumpDivisibility.lean` は `Monoid.exponent` を経由したが、
原文は **`θ_0(σ) = u` の位数**を見ている。ここではその形で書く。 -/

/-- ★★**原文の計算の全部**——`G` 可換・`σ ∈ G_n \ G_{n+1}` なら、
**任意の** `τ ∈ G_0` について `orderOf θ_0(τ) ∣ n`。

原文:
> If G is abelian, then στσ⁻¹ = τ, hence a ≡ u^n a mod 𝔭^{n+1}.

の 1 行そのもの。`residue_ramCoeff_conj`(共役の変換則
`θ_n(τστ⁻¹) = θ_0(τ)^n θ_n(σ)`)に `τστ⁻¹ = σ` を入れ、
`orderOf_dvd_of_pow_mul_eq_self` で `a` を約す。

★★**`G` 可換を使うのは `hconj` ただ 1 行**である。
★`hσ' : σ ∉ G_{n+1}` を使うのは `hc : θ_n(σ) ≠ 0` ただ 1 行である。 -/
theorem orderOf_thetaMul_dvd_of_notMem {A B : Type*} [CommRing A] [CommRing B] [Algebra A B]
    [IsDomain B] [IsLocalRing B] {G : Type*} [Group G] [MulSemiringAction G B]
    [SMulCommClass G A B] {α : B} (hα : maximalIdeal B = Ideal.span {α}) (hα0 : α ≠ 0)
    (hadj : Algebra.adjoin A ({α} : Set B) = ⊤) (habel : ∀ x y : G, x * y = y * x)
    {n : ℕ} {σ : G} (hσ : σ ∈ lowerRamificationGroup B G n)
    (hσ' : σ ∉ lowerRamificationGroup B G (n + 1))
    (τ : lowerRamificationGroup B G 0) :
    orderOf (thetaMul hα hα0 τ) ∣ n := by
  have hc : residue B (ramCoeff α n σ) ≠ 0 := fun h =>
    hσ' ((residue_ramCoeff_eq_zero_iff hα hα0 hadj hσ).mp h)
  have hconj : (τ : G) * σ * (τ : G)⁻¹ = σ := by
    rw [habel (τ : G) σ, mul_assoc, mul_inv_cancel, mul_one]
  have h := residue_ramCoeff_conj hα hα0 hσ τ.2
  rw [hconj] at h
  refine orderOf_dvd_of_pow_mul_eq_self (thetaMul hα hα0 τ) hc ?_
  rw [thetaMul_apply_coe]
  exact h.symm

/-- ★**原文「choose σ ∈ G which generates G_0/G_1, i.e. θ_0(σ) = u mod p has order e_0」**。

`G_0/G_1` は巡回(`isCyclic_quot_lowerRamificationGroup_zero`)なので、
その生成元を `G_0` に持ち上げると `θ_0` の像の位数がちょうど `e_0 = |G_0/G_1|` になる。 -/
theorem exists_orderOf_thetaMul_eq_natCard_quot {A B : Type*} [CommRing A] [CommRing B]
    [Algebra A B] [IsDomain B] [IsLocalRing B] {G : Type*} [Group G] [MulSemiringAction G B]
    [SMulCommClass G A B] {α : B} (hα : maximalIdeal B = Ideal.span {α}) (hα0 : α ≠ 0)
    (hadj : Algebra.adjoin A ({α} : Set B) = ⊤)
    [Finite (lowerRamificationGroup B G 0 ⧸
      (lowerRamificationGroup B G 1).subgroupOf (lowerRamificationGroup B G 0))] :
    ∃ τ : lowerRamificationGroup B G 0, orderOf (thetaMul hα hα0 τ)
      = Nat.card (lowerRamificationGroup B G 0 ⧸
        (lowerRamificationGroup B G 1).subgroupOf (lowerRamificationGroup B G 0)) := by
  haveI := isCyclic_quot_lowerRamificationGroup_zero (A := A) (G := G) hα hα0 hadj
  obtain ⟨q, hq⟩ := exists_orderOf_map_eq_natCard (thetaMulQuot (G := G) hα hα0 hadj)
    (thetaMulQuot_injective hα hα0 hadj)
  obtain ⟨τ, rfl⟩ := QuotientGroup.mk_surjective q
  exact ⟨τ, hq⟩

/-- ★★**Yoshida Corollary 6.3(原文どおりの証明)**——`G` 可換で `σ ∈ G_n \ G_{n+1}` なら
`e_0 := |G_0/G_1|` は `n` を割る。

★結論は `RamificationJumpDivisibility.lean` の `nat_card_quot_zero_dvd_of_notMem` と同じだが、
**証明の経路が違う**: あちらは「全ての元が `n` 乗で `1`」→ 指数 = 位数、
こちらは原文どおり「生成元を取り、その像の位数が `e_0`」。 -/
theorem nat_card_quot_zero_dvd_of_notMem_generator {A B : Type*} [CommRing A] [CommRing B]
    [Algebra A B] [IsDomain B] [IsLocalRing B] {G : Type*} [Group G] [MulSemiringAction G B]
    [SMulCommClass G A B] {α : B} (hα : maximalIdeal B = Ideal.span {α}) (hα0 : α ≠ 0)
    (hadj : Algebra.adjoin A ({α} : Set B) = ⊤) (habel : ∀ x y : G, x * y = y * x)
    [Finite (lowerRamificationGroup B G 0 ⧸
      (lowerRamificationGroup B G 1).subgroupOf (lowerRamificationGroup B G 0))]
    {n : ℕ} {σ : G} (hσ : σ ∈ lowerRamificationGroup B G n)
    (hσ' : σ ∉ lowerRamificationGroup B G (n + 1)) :
    Nat.card (lowerRamificationGroup B G 0 ⧸
      (lowerRamificationGroup B G 1).subgroupOf (lowerRamificationGroup B G 0)) ∣ n := by
  obtain ⟨τ, hτ⟩ := exists_orderOf_thetaMul_eq_natCard_quot (A := A) (G := G) hα hα0 hadj
  rw [← hτ]
  exact orderOf_thetaMul_dvd_of_notMem hα hα0 hadj habel hσ hσ' τ

/-- 上の `≠` 版(`G_n ≠ G_{n+1}` から元を取り出す)。
★★仮定は `G_n ≠ G_{n+1}`(**`=` ではない**。ファイル冒頭の警告)。 -/
theorem nat_card_quot_zero_dvd_of_ne_generator {A B : Type*} [CommRing A] [CommRing B]
    [Algebra A B] [IsDomain B] [IsLocalRing B] {G : Type*} [Group G] [MulSemiringAction G B]
    [SMulCommClass G A B] {α : B} (hα : maximalIdeal B = Ideal.span {α}) (hα0 : α ≠ 0)
    (hadj : Algebra.adjoin A ({α} : Set B) = ⊤) (habel : ∀ x y : G, x * y = y * x)
    [Finite (lowerRamificationGroup B G 0 ⧸
      (lowerRamificationGroup B G 1).subgroupOf (lowerRamificationGroup B G 0))]
    {n : ℕ}
    (hne : lowerRamificationGroup B G n ≠ lowerRamificationGroup B G (n + 1)) :
    Nat.card (lowerRamificationGroup B G 0 ⧸
      (lowerRamificationGroup B G 1).subgroupOf (lowerRamificationGroup B G 0)) ∣ n := by
  obtain ⟨σ, hσ, hσ'⟩ :=
    exists_mem_notMem_of_antitone_ne (lowerRamificationGroup_antitone B G) hne
  exact nat_card_quot_zero_dvd_of_notMem_generator hα hα0 hadj habel hσ hσ'

/-! ## §3 Hasse-Arf(Theorem 6.11)最終段への受け渡し

原文:
> For any i ∈ Z≥1 (where H_i = G_i) with H_i ≠ H_{i+1}, we have e_0 | i by Corollary 6.3,
> hence e_0 | Σ_{i=1}^n |H_i|. As e_0 and |H| are coprime, we have e_0 | φ_H(n)
> by Lemma 6.10(i).

★完全分岐(`G_0 = ⊤`)を仮定する。原典 `§6.1` の設定「`K′/K` は完全分岐」がそれで、
このとき `e_0 = |G/G_1|` になる。 -/

/-- **Corollary 6.3 の対偶**——`e_0 ∤ i` なら `G_i = G_{i+1}`。 -/
theorem lowerRamificationGroup_eq_succ_of_not_dvd {A B : Type*} [CommRing A] [CommRing B]
    [Algebra A B] [IsDomain B] [IsLocalRing B] {G : Type*} [Group G] [MulSemiringAction G B]
    [SMulCommClass G A B] {α : B} (hα : maximalIdeal B = Ideal.span {α}) (hα0 : α ≠ 0)
    (hadj : Algebra.adjoin A ({α} : Set B) = ⊤) (h0 : lowerRamificationGroup B G 0 = ⊤)
    (habel : ∀ x y : G, x * y = y * x) [Finite (G ⧸ lowerRamificationGroup B G 1)] {i : ℕ}
    (hi : ¬ (Nat.card (G ⧸ lowerRamificationGroup B G 1) ∣ i)) :
    lowerRamificationGroup B G i = lowerRamificationGroup B G (i + 1) := by
  by_contra hne
  exact hi (nat_card_quot_top_dvd_of_ne hα hα0 hadj h0 habel hne)

/-- ★★**原文の "hence e_0 | Σ_{i=1}^n |H_i|"**——`G` 可換・完全分岐・`G_n ≠ G_{n+1}` なら
`e_0 = |G/G_1|` は `Σ_{i=1}^n |G_i|` を割る。

★2 つの入力を合成しただけ:
`n` 自身が `e_0` の倍数であること(Corollary 6.3)と、
跳びが `e_0` の倍数でしか起きないこと(Corollary 6.3 の対偶)。
どちらか一方だけでは出ない(`dvd_sum_Icc_of_step` の退化検査)。 -/
theorem dvd_sum_natCard_lowerRamificationGroup {A B : Type*} [CommRing A] [CommRing B]
    [Algebra A B] [IsDomain B] [IsLocalRing B] {G : Type*} [Group G] [MulSemiringAction G B]
    [SMulCommClass G A B] {α : B} (hα : maximalIdeal B = Ideal.span {α}) (hα0 : α ≠ 0)
    (hadj : Algebra.adjoin A ({α} : Set B) = ⊤) (h0 : lowerRamificationGroup B G 0 = ⊤)
    (habel : ∀ x y : G, x * y = y * x) [Finite (G ⧸ lowerRamificationGroup B G 1)] {n : ℕ}
    (hne : lowerRamificationGroup B G n ≠ lowerRamificationGroup B G (n + 1)) :
    Nat.card (G ⧸ lowerRamificationGroup B G 1) ∣
      ∑ i ∈ Finset.Icc 1 n, Nat.card (lowerRamificationGroup B G i) := by
  refine dvd_sum_Icc_of_step _ ?_ (nat_card_quot_top_dvd_of_ne hα hα0 hadj h0 habel hne)
  intro i _ hi
  rw [lowerRamificationGroup_eq_succ_of_not_dvd hα hα0 hadj h0 habel hi]

/-- `e_0 = |G/G_1|` は剰余体の標数 `p` で割れない
(`θ_0 : G/G_1 ↪ 𝓀ˣ` と `not_dvd_natCard_of_injective_units`)。 -/
theorem not_dvd_natCard_quot_lowerRamificationGroup_one {A B : Type*} [CommRing A] [CommRing B]
    [Algebra A B] [IsDomain B] [IsLocalRing B] {G : Type*} [Group G] [MulSemiringAction G B]
    [SMulCommClass G A B] (p : ℕ) [Fact p.Prime] [CharP (ResidueField B) p] {α : B}
    (hα : maximalIdeal B = Ideal.span {α}) (hα0 : α ≠ 0)
    (hadj : Algebra.adjoin A ({α} : Set B) = ⊤) (h0 : lowerRamificationGroup B G 0 = ⊤)
    [Finite (G ⧸ lowerRamificationGroup B G 1)] :
    ¬ p ∣ Nat.card (G ⧸ lowerRamificationGroup B G 1) :=
  not_dvd_natCard_of_injective_units p _ (thetaMulTop_injective hα hα0 hadj h0)

/-- ★★**原文の "As e_0 and |H| are coprime"**(`H = G_1`)。

`e_0` は `p` と素(上)、`|G_1|` は `p` 冪(`isPGroup_lowerRamificationGroup_one`)。 -/
theorem coprime_natCard_quot_natCard_lowerRamificationGroup_one {A B : Type*} [CommRing A]
    [CommRing B] [Algebra A B] [IsDomain B] [IsLocalRing B] [IsNoetherianRing B] {G : Type*}
    [Group G] [MulSemiringAction G B] [SMulCommClass G A B] [Finite G] [FaithfulSMul G B]
    (p : ℕ) [Fact p.Prime] [CharP (ResidueField B) p] {α : B}
    (hα : maximalIdeal B = Ideal.span {α}) (hα0 : α ≠ 0)
    (hadj : Algebra.adjoin A ({α} : Set B) = ⊤) (h0 : lowerRamificationGroup B G 0 = ⊤) :
    Nat.Coprime (Nat.card (G ⧸ lowerRamificationGroup B G 1))
      (Nat.card (lowerRamificationGroup B G 1)) := by
  obtain ⟨m, hm⟩ := (IsPGroup.iff_card (p := p)).mp
    (isPGroup_lowerRamificationGroup_one (A := A) (G := G) p hα hα0 hadj)
  rw [hm]
  exact Nat.Coprime.pow_right m
    (((Nat.Prime.coprime_iff_not_dvd Fact.out).mpr
      (not_dvd_natCard_quot_lowerRamificationGroup_one (A := A) p hα hα0 hadj h0)).symm)

/-- ★★★**Hasse-Arf(Theorem 6.11)最終段の結論**——
`Σ_{i=1}^n |G_i| = |G_1| · k` なら `e_0 ∣ k`。

原文の `k` は `φ_H(n)`(`H = G_1`)であり、`Σ_{i=1}^n |H_i| = |H| · φ_H(n)` は
**Lemma 6.10(i)** が与える。★本ファイルは Lemma 6.10 を import せず、
その等式を**仮定 `hk` として受け取る**(逸脱の記録 4)。
Lemma 6.10(i) が出来た時点で `k := φ_H(n)` を代入すれば原文の
「`e_0 | φ_H(n)`」がそのまま出る。 -/
theorem dvd_of_sum_eq_natCard_mul {A B : Type*} [CommRing A] [CommRing B] [Algebra A B]
    [IsDomain B] [IsLocalRing B] [IsNoetherianRing B] {G : Type*} [Group G]
    [MulSemiringAction G B] [SMulCommClass G A B] [Finite G] [FaithfulSMul G B]
    (p : ℕ) [Fact p.Prime] [CharP (ResidueField B) p] {α : B}
    (hα : maximalIdeal B = Ideal.span {α}) (hα0 : α ≠ 0)
    (hadj : Algebra.adjoin A ({α} : Set B) = ⊤) (h0 : lowerRamificationGroup B G 0 = ⊤)
    (habel : ∀ x y : G, x * y = y * x) {n : ℕ}
    (hne : lowerRamificationGroup B G n ≠ lowerRamificationGroup B G (n + 1)) {k : ℕ}
    (hk : ∑ i ∈ Finset.Icc 1 n, Nat.card (lowerRamificationGroup B G i)
      = Nat.card (lowerRamificationGroup B G 1) * k) :
    Nat.card (G ⧸ lowerRamificationGroup B G 1) ∣ k := by
  refine (coprime_natCard_quot_natCard_lowerRamificationGroup_one (A := A) p hα hα0 hadj
    h0).dvd_of_dvd_mul_left ?_
  rw [← hk]
  exact dvd_sum_natCard_lowerRamificationGroup hα hα0 hadj h0 habel hne

/-! ## §4 `PAdicLocalField` への具体化

`A := 𝒪[K.carrier]`、`B := adjoinIntegers K x`、`G := Gal(K(x)/K)`。
素元は `∃` の内側に閉じ込める(取り方に依らないことは Prop 6.2 の `§5` が保証する)。 -/

variable {p : ℕ} [Fact p.Prime]

/-- ★★**原文の "hence e_0 | Σ_{i=1}^n |H_i|"(`PAdicLocalField` 版)**。 -/
theorem dvd_sum_natCard_lowerRamificationGroupAdjoin (K : PAdicLocalField p) (x : K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (ht : IsTotallyRamifiedAdjoin K x)
    (habel : ∀ σ τ : ((IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure))),
      σ * τ = τ * σ)
    {n : ℕ}
    (hne : lowerRamificationGroupAdjoin K x n ≠ lowerRamificationGroupAdjoin K x (n + 1)) :
    Nat.card (((IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure))) ⧸
      lowerRamificationGroupAdjoin K x 1) ∣
      ∑ i ∈ Finset.Icc 1 n, Nat.card (lowerRamificationGroupAdjoin K x i) := by
  obtain ⟨α, hspan, hα0, hadj⟩ := exists_uniformizer_ne_zero_adjoin_eq K x ht
  exact dvd_sum_natCard_lowerRamificationGroup hspan hα0 hadj
    (lowerRamificationGroupAdjoin_zero_eq_top K x ht) habel hne

/-- ★**原文の "As e_0 and |H| are coprime"(`PAdicLocalField` 版)**。 -/
theorem coprime_natCard_quot_lowerRamificationGroupAdjoin_one (K : PAdicLocalField p)
    (x : K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (ht : IsTotallyRamifiedAdjoin K x) :
    Nat.Coprime (Nat.card (((IntermediateField.adjoin K.carrier ({x} : Set K.closure))
        ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure))) ⧸
      lowerRamificationGroupAdjoin K x 1))
      (Nat.card (lowerRamificationGroupAdjoin K x 1)) := by
  haveI := isDiscreteValuationRing_adjoinIntegers K x
  haveI := charP_residueField_adjoinIntegers K x
  obtain ⟨α, hspan, hα0, hadj⟩ := exists_uniformizer_ne_zero_adjoin_eq K x ht
  exact coprime_natCard_quot_natCard_lowerRamificationGroup_one p hspan hα0 hadj
    (lowerRamificationGroupAdjoin_zero_eq_top K x ht)

/-- ★★★**Hasse-Arf 最終段の結論(`PAdicLocalField` 版)**——
`Σ_{i=1}^n |G_i| = |G_1| · k` なら `e_0 ∣ k`。Lemma 6.10(i) を `hk` として差し込む形。 -/
theorem dvd_of_sum_eq_natCard_mul_adjoin (K : PAdicLocalField p) (x : K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (ht : IsTotallyRamifiedAdjoin K x)
    (habel : ∀ σ τ : ((IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure))),
      σ * τ = τ * σ)
    {n : ℕ}
    (hne : lowerRamificationGroupAdjoin K x n ≠ lowerRamificationGroupAdjoin K x (n + 1))
    {k : ℕ}
    (hk : ∑ i ∈ Finset.Icc 1 n, Nat.card (lowerRamificationGroupAdjoin K x i)
      = Nat.card (lowerRamificationGroupAdjoin K x 1) * k) :
    Nat.card (((IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure))) ⧸
      lowerRamificationGroupAdjoin K x 1) ∣ k := by
  refine (coprime_natCard_quot_lowerRamificationGroupAdjoin_one K x ht).dvd_of_dvd_mul_left ?_
  rw [← hk]
  exact dvd_sum_natCard_lowerRamificationGroupAdjoin K x ht habel hne

end ABC3.Found.PGC
