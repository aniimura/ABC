import ABC3.Found.PGC.UniformizerExpansion

/-!
# Sen の定理 (i)(ii) —— 跳び位置 `i_j` の狭義増加と `i(σ^a) = i_{v_p(a)}`(Yoshida 2008 Prop 6.6)

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

★`pdftotext` は生成部分群の角括弧 `⟨ ⟩` を落とす(`|⟨σ⟩| = p^m` が `|σ| = p^m` にしか見えない)。
上の逐語は `data-txt` 適用済み(= `pdftotext` に見えている形)である。

原文の証明 (i)(ii) の部分(同 p.15):
> Proof. (i): Lemma 6.4 for α = σ^p^j−1(π) −π shows i_j−1 < i_j. We have σ^p^j ⊂ H_n if and
> only if σ^p^j ∈ H_n, i.e. i_j > n. As all subgroups of σ are of the form σ^p^j, we have
> σ^p^j ⊃ H_n ⇔ σ^p^j−1 ̸⊂ H_n ⇔ i_j−1 ≤ n. (ii): This is ∞ = ∞ if p^m | a. If j := v_p(a) < m,
> then H_i_j−1 = σ^p^j and H_i_j = σ^p^j+1 by (i), therefore σ^a ∈ H_i_j−1 \ H_i_j,
> i.e. i(σ^a) = i_j.

## 本ファイルの範囲

**(i) 前半**・**(i) 後半**・**(ii)** を埋めた。**(iii) は入っていない**(別ノード)。

## 何を証明したか

**§1 抽象核 —— 純群論(分岐・付値・Galois の語彙が 1 つも出てこない)**

`orderOf σ = p^m`(`p` 素数)だけを仮定した巡回 `p` 群の束:

* `pow_pow_eq_one_iff` : `σ^{p^j} = 1 ↔ m ≤ j`。
* `exists_orderOf_pow_eq` : `orderOf (σ^b)` も `p` 冪。
* `exists_pow_eq_of_mem_zpowers` : `⟨σ⟩` の元は**自然数**冪で書ける(有限位数)。
* ★`zpowers_pow_eq_zpowers_pow` : `p^k ∥ a ⟹ ⟨σ^a⟩ = ⟨σ^{p^k}⟩`。**(ii) の全内容**。
* ★`pow_pow_notMem_zpowers_pow_succ` : `j < m ⟹ σ^{p^j} ∉ ⟨σ^{p^{j+1}}⟩`(狭義降下)。

★原典が「As all subgroups of `⟨σ⟩` are of the form `⟨σ^{p^j}⟩`」で畳んでいる
**部分群の分類は使っていない**。必要なのは上の 2 本(生成部分群の一致と狭義降下)だけで、
分類定理そのものは (i)(ii) の消費側に現れない。★見積(110〜180 行)に対して実測 40 行。

**§2 抽象核 —— `ℕ∞` の外延性と望遠鏡和**

* `enat_le_of_forall_natCast_lt` / `enat_eq_of_forall_natCast_lt_iff` :
  `ℕ∞` の元は「下にある自然数の集合」で決まる。★mathlib の
  `ENat.forall_natCast_le_iff_le` は `≤` 版で、`<` 版は無い。
* `sum_range_smul_smul_sub` : `Σ_{i<n} σ^i • (σ•x − x) = σ^n • x − x`。
  一般の加法群 + `DistribMulAction` で成り立つ(`Finset.sum_range_sub` の 1 行)。

**§3 Prop 6.6 (i) 前半**

* ★`ramIndex_lt_ramIndex_pow_char` : `σ ∈ G_1`、`σ•π ≠ π` ⟹ `i(σ) < i(σ^p)`。
  ★原文の「Lemma 6.4 for `α = σ^{p^{j−1}}(π) − π`」がこれである。`α := σ•π − π` に
  Lemma 6.4(`Found/PGC/ConjugateSumValuation.lean`)を当て、望遠鏡和 §2 で
  `Σ_{i<p} σ^i • α = σ^p • π − π` に潰す。
* `ramIndex_pow_pow_lt_succ` : `j < m ⟹ i_j < i_{j+1}`(原文の `i_{j−1} < i_j`, `j ≤ m`)。

**§4 `⊤` 判定と単調性**

* `ramIndex_pow_pow_eq_top_iff` : `i_j = ⊤ ↔ m ≤ j`。
  ★原文の規約「`i_j := ∞ for j ≥ m`」は**定義ではなく定理**になる(Lean では自動)。
* `monotone_ramIndex_pow_pow` : `j ↦ i_j` は単調(`m` 以降は `⊤` で一定)。

**§5 Prop 6.6 (ii)**

* `ramIndex_eq_of_zpowers_eq` : `⟨τ⟩ = ⟨τ'⟩ ⟹ i(τ) = i(τ')`。
  ★`σ ∈ G_n ⟺ ⟨σ⟩ ≤ G_n` と `ℕ∞` の外延性(§2)だけ。
* ★`ramIndex_pow_eq_ramIndex_pow_pow` : `p^k ∥ a ⟹ i(σ^a) = i_k`。
  ★★**原文は (ii) を (i) 経由で証明している**(`H_{i_j−1} = ⟨σ^{p^j}⟩` を使う)が、
  ここでは **(i) を一切使わない**。`⟨σ^a⟩ = ⟨σ^{p^k}⟩`(§1)から直ちに出る
  ——`i` は生成部分群だけで決まるからである。★原文より依存が浅い(逸脱ではなく別証明)。

**§6 Prop 6.6 (i) 後半**

* ★`lowerRamificationGroup_inf_zpowers_eq_iff` :
  `j + 1 ≤ m` のとき `G_n ⊓ ⟨σ⟩ = ⟨σ^{p^{j+1}}⟩ ↔ (i_j ≤ n ∧ n < i_{j+1})`。
  ★添字を 1 つずらして書いてある(原文 `H_n = ⟨σ^{p^j}⟩ ⟺ i_{j−1} ≤ n < i_j`、`1 ≤ j ≤ m`)。
  **`ℕ` の引き算を出さないため**である(落とし穴の回避)。
* `j ≤ m` は落とせない: `j > m` では左辺の `⟨σ^{p^j}⟩` が `⊥` に潰れ、
  `n` が大きいとき左辺は真・右辺は偽になる(下の退化の自己検査 D2)。

**§7 `PAdicLocalField` 具体層**

`B := adjoinIntegers K x`、`A := 𝒪[K.carrier]`、`G := Gal(K(x)/K)`。
`orderOf σ = p^m` は原文が "(by Proposition 6.2)" で畳んでいる部分で、
`exists_orderOf_eq_pow_of_mem_lowerRamificationGroup_one`(`RamificationJumpDivisibility.lean`)
がそのまま供給する。

## ★逸脱の記録

1. **`i_j := ∞ for j ≥ m` を規約ではなく定理にした**(§4)。原典は `j ≥ m` での `i_j` を
   記号的に `∞` と**定める**が、`ramIndex` の定義(`addVal (σ•π − π)`、`σ = 1` で `⊤`)の下では
   これは**証明できる命題**である。したがって規約を入れずに済ませた。
2. **(i) 後半の添字を 1 つずらした**。原文 `H_n = ⟨σ^{p^j}⟩ ⟺ i_{j−1} ≤ n < i_j`(`1 ≤ j ≤ m`)を
   `j + 1 ≤ m` の下で `H_n = ⟨σ^{p^{j+1}}⟩ ⟺ i_j ≤ n < i_{j+1}` と書いた。`ℕ` の切り詰め
   引き算 `j − 1` を出さないため(`j = 0` で空虚に真になる事故の回避)。数学的内容は同じ。
3. **(ii) の `a ≥ 1` を落とした**。原文は `a ≥ 1` を要求するが、仮定 `¬ p^{k+1} ∣ a` が
   すでに `a ≠ 0` を含む(`p^{k+1} ∣ 0` は常に真)ので冗長である。★仮定を減らした方向。
4. **(ii) を `v_p` の値ではなく整除性の述語で書いた**。原文は `i(σ^a) = i_{v_p(a)}` と
   `v_p : ℕ → ℕ` を関数として使うが、本ファイルは `p^k ∣ a ∧ ¬ p^{k+1} ∣ a` を仮定に置く。
   ★`padicValNat p 0 = 0` というジャンク値を避けるため(数学的内容は同じ)。
5. **(ii) の証明が原文と違う**(上の §5)。原文は (i) を経由するが、ここでは
   `⟨σ^a⟩ = ⟨σ^{p^k}⟩` から直接出す。★原文より依存が浅い。
6. **`m ≥ 1` を仮定していない**。原文は `σ ≠ 1` から `m ≥ 1` を出して置いているが、
   本ファイルの主張はどれも `m = 0`(= `σ = 1`)で自明に真か、仮定 `j < m` が空になる。

## ★退化の自己検査

* **(D1) (i) 前半で `σ ∈ G_1` を落とすと偽**。Lemma 6.4 は `σ ∈ G_1`(= `p` が剰余体で 0)を
  使う。落とすと `i(σ^p) = i(σ)` の例(`p` が可逆な tame な場合、`σ^p` が `σ` と同じ跳び)を
  排除できない。★`hσ` は §3 の中で本質的に 1 度だけ使われる。
* **(D2) (i) 後半で `j + 1 ≤ m` を落とすと偽**。`j + 1 > m` なら `σ^{p^{j+1}} = 1` で
  右辺の `⟨σ^{p^{j+1}}⟩ = ⊥`、また `i_{j+1} = ⊤` だから右辺の条件 `n < i_{j+1}` は真、
  `i_j ≤ n` も `j ≥ m` では `⊤ ≤ n` で偽。ところが `n` を十分大きく取れば
  `G_n = ⊥` なので左辺は真になる。★仮定 `j + 1 ≤ m` が消せない理由である。
* **(D3) (ii) で `¬ p^{k+1} ∣ a` を落とすと偽**。`k := 0`、`a := p` と取ると
  結論は `i(σ^p) = i(σ)` になり、(i) 前半(`m ≥ 2` のとき狭義増加)と矛盾する。
-/

namespace ABC3.Found.PGC

open IsLocalRing ABC3.Skeleton.PGC IsDiscreteValuationRing
open scoped NNReal Valued

def ramIndex_pow_pow_lt_succ.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 15, item := "Proposition 6.6", sectionId := "prop-6-6" }

/-! ## §1 抽象核 —— 巡回 `p` 群の束(純群論)

★分岐・付値・Galois の語彙は 1 つも現れない。仮定は `orderOf σ = p^m`(`p` 素数)だけ。 -/

/-- `σ^{p^j} = 1 ⟺ m ≤ j`。★原文の「`i_j := ∞ for j ≥ m`」の群論側の根拠。 -/
theorem pow_pow_eq_one_iff {G : Type*} [Group G] {σ : G} {p m : ℕ} (hp : 1 < p)
    (hord : orderOf σ = p ^ m) (j : ℕ) : σ ^ p ^ j = 1 ↔ m ≤ j := by
  rw [← orderOf_dvd_iff_pow_eq_one, hord, Nat.pow_dvd_pow_iff_le_right hp]

/-- `orderOf σ = p^m` なら `orderOf (σ^b)` も `p` 冪。 -/
theorem exists_orderOf_pow_eq {G : Type*} [Group G] {σ : G} {p m : ℕ} (hp : p.Prime)
    (hord : orderOf σ = p ^ m) (b : ℕ) : ∃ e, orderOf (σ ^ b) = p ^ e := by
  have h : orderOf (σ ^ b) ∣ p ^ m := by
    rw [← hord]
    exact orderOf_dvd_of_pow_eq_one
      (by rw [← pow_mul, mul_comm, pow_mul, pow_orderOf_eq_one, one_pow])
  obtain ⟨e, _, he⟩ := (Nat.dvd_prime_pow hp).1 h
  exact ⟨e, he⟩

/-- `⟨σ⟩` の元は**自然数**冪で書ける(`σ` は有限位数)。 -/
theorem exists_pow_eq_of_mem_zpowers {G : Type*} [Group G] {σ : G} {p m : ℕ} (hp : 0 < p)
    (hord : orderOf σ = p ^ m) {x : G} (hx : x ∈ Subgroup.zpowers σ) : ∃ a : ℕ, σ ^ a = x := by
  have hfin : IsOfFinOrder σ := orderOf_pos_iff.1 (by rw [hord]; exact Nat.pow_pos hp)
  have hmem : x ∈ (Submonoid.powers σ : Set G) := by rw [hfin.powers_eq_zpowers]; exact hx
  exact (Submonoid.mem_powers_iff x σ).1 hmem

/-- ★★**Prop 6.6 (ii) の群論的中身** —— `p^k ∥ a` なら `⟨σ^a⟩ = ⟨σ^{p^k}⟩`。

`a = p^k u`(`p ∤ u`)と書くと `σ^a = (σ^{p^k})^u` で、`orderOf (σ^{p^k})` は `p` 冪だから
`u` と互いに素。`mem_zpowers_pow_iff`(gcd 判定)で逆包含が出る。 -/
theorem zpowers_pow_eq_zpowers_pow {G : Type*} [Group G] {σ : G} {p m : ℕ} (hp : p.Prime)
    (hord : orderOf σ = p ^ m) {a k : ℕ} (hdvd : p ^ k ∣ a) (hndvd : ¬ p ^ (k + 1) ∣ a) :
    Subgroup.zpowers (σ ^ a) = Subgroup.zpowers (σ ^ p ^ k) := by
  obtain ⟨u, rfl⟩ := hdvd
  have hu : ¬ p ∣ u := fun ⟨c, hc⟩ => hndvd ⟨c, by rw [hc]; ring⟩
  rw [pow_mul]
  refine le_antisymm (Subgroup.zpowers_le.2 (pow_mem (Subgroup.mem_zpowers _) u))
    (Subgroup.zpowers_le.2 (mem_zpowers_pow_iff.2 ?_))
  obtain ⟨e, he⟩ := exists_orderOf_pow_eq hp hord (p ^ k)
  rw [he]
  exact Nat.Coprime.pow_right e ((Nat.Prime.coprime_iff_not_dvd hp).2 hu).symm

/-- ★**狭義降下** —— `j < m` なら `σ^{p^j} ∉ ⟨σ^{p^{j+1}}⟩`。

原文が「As all subgroups of `⟨σ⟩` are of the form `⟨σ^{p^j}⟩`」で畳んだ部分のうち、
(i) 後半が実際に必要とするのはこの 1 本だけである(部分群の分類は要らない)。 -/
theorem pow_pow_notMem_zpowers_pow_succ {G : Type*} [Group G] {σ : G} {p m : ℕ} (hp : p.Prime)
    (hord : orderOf σ = p ^ m) {j : ℕ} (hj : j < m) :
    σ ^ p ^ j ∉ Subgroup.zpowers (σ ^ p ^ (j + 1)) := by
  intro hmem
  rw [pow_succ, pow_mul] at hmem
  have hgcd := mem_zpowers_pow_iff.1 hmem
  obtain ⟨e, he⟩ := exists_orderOf_pow_eq hp hord (p ^ j)
  have hne : σ ^ p ^ j ≠ 1 := fun h =>
    absurd ((pow_pow_eq_one_iff hp.one_lt hord j).1 h) (by omega)
  have he1 : e ≠ 0 := by
    intro h
    rw [h, pow_zero, orderOf_eq_one_iff] at he
    exact hne he
  have hpd : p ∣ Nat.gcd p (orderOf (σ ^ p ^ j)) := Nat.dvd_gcd dvd_rfl (he ▸ dvd_pow_self p he1)
  rw [hgcd] at hpd
  exact absurd (Nat.le_of_dvd one_pos hpd) (not_le.2 hp.one_lt)

/-! ## §2 抽象核 —— `ℕ∞` の外延性と望遠鏡和 -/

/-- `ℕ∞` の元は「下にある自然数の集合」で決まる(`≤` 側)。
★mathlib の `ENat.forall_natCast_le_iff_le` は `≤` 版で、`<` 版は無い。 -/
theorem enat_le_of_forall_natCast_lt {x y : ℕ∞} (h : ∀ n : ℕ, (n : ℕ∞) < x → (n : ℕ∞) < y) :
    x ≤ y := by
  by_contra hxy
  rw [not_le] at hxy
  obtain ⟨n, rfl⟩ : ∃ n : ℕ, y = (n : ℕ∞) :=
    ⟨y.toNat, (ENat.coe_toNat (by rintro rfl; exact absurd hxy (by simp))).symm⟩
  exact absurd (h n hxy) (lt_irrefl _)

/-- `ℕ∞` の外延性(`<` 版)。 -/
theorem enat_eq_of_forall_natCast_lt_iff {x y : ℕ∞} (h : ∀ n : ℕ, (n : ℕ∞) < x ↔ (n : ℕ∞) < y) :
    x = y :=
  le_antisymm (enat_le_of_forall_natCast_lt fun n hn => (h n).1 hn)
    (enat_le_of_forall_natCast_lt fun n hn => (h n).2 hn)

/-- **望遠鏡和** `Σ_{i<n} σ^i • (σ•x − x) = σ^n • x − x`。
★一般の加法群と一般の `DistribMulAction` で成り立つ(`Finset.sum_range_sub` の 1 行)。 -/
theorem sum_range_smul_smul_sub {M : Type*} [AddCommGroup M] {G : Type*} [Monoid G]
    [DistribMulAction G M] (σ : G) (x : M) (n : ℕ) :
    ∑ i ∈ Finset.range n, σ ^ i • (σ • x - x) = σ ^ n • x - x := by
  simpa [smul_sub, pow_succ, mul_smul] using Finset.sum_range_sub (fun i => σ ^ i • x) n

/-! ## §3 Yoshida Prop 6.6 (i) 前半 —— `i_{j−1} < i_j`

原文は「Lemma 6.4 for `α = σ^{p^{j−1}}(π) − π` shows `i_{j−1} < i_j`」の 1 行。
その 1 行の中身が `ramIndex_lt_ramIndex_pow_char` である。 -/

/-- ★★★**Prop 6.6 (i) 前半の心臓** —— `σ ∈ G_1` かつ `σ•π ≠ π` なら `i(σ) < i(σ^p)`。

`α := σ•π − π ≠ 0` に Lemma 6.4(`addVal_lt_addVal_sum_smul_pow_char`)を当てると
`v(α) < v(Σ_{i<p} σ^i • α)` で、望遠鏡和(§2)より `Σ_{i<p} σ^i • α = σ^p•π − π`。
★`σ ∈ G_1` はここで 1 度だけ使う(退化検査 D1)。 -/
theorem ramIndex_lt_ramIndex_pow_char {B : Type*} [CommRing B] [IsDomain B]
    [IsDiscreteValuationRing B] {G : Type*} [Group G] [MulSemiringAction G B] (p : ℕ)
    [CharP (ResidueField B) p] {π : B} (hπ : maximalIdeal B = Ideal.span {π}) {σ : G}
    (hσ : σ ∈ lowerRamificationGroup B G 1) (hne : σ • π ≠ π) :
    ramIndex π σ < ramIndex π (σ ^ p) := by
  have h := addVal_lt_addVal_sum_smul_pow_char p hπ hσ (sub_ne_zero.2 hne)
  rw [sum_range_smul_smul_sub] at h
  exact h

/-- ★★★**Yoshida 2008 Prop 6.6 (i) 前半** —— `j < m` なら `i_j < i_{j+1}`。

原文は `i_{j−1} < i_j if j ≤ m`。添字を 1 つずらして `ℕ` の引き算を避けている(逸脱 2)。 -/
theorem ramIndex_pow_pow_lt_succ {A B : Type*} [CommRing A] [CommRing B] [Algebra A B]
    [IsDomain B] [IsDiscreteValuationRing B] {G : Type*} [Group G] [MulSemiringAction G B]
    [SMulCommClass G A B] [FaithfulSMul G B] (p : ℕ) (hp : p.Prime)
    [CharP (ResidueField B) p] {π : B} (huni : maximalIdeal B = Ideal.span {π})
    (hadj : Algebra.adjoin A ({π} : Set B) = ⊤) {σ : G}
    (hσ : σ ∈ lowerRamificationGroup B G 1) {m : ℕ} (hord : orderOf σ = p ^ m)
    {j : ℕ} (hj : j < m) :
    ramIndex π (σ ^ p ^ j) < ramIndex π (σ ^ p ^ (j + 1)) := by
  have hne : (σ ^ p ^ j) • π ≠ π := fun h =>
    absurd ((pow_pow_eq_one_iff hp.one_lt hord j).1
      (eq_one_of_smul_eq_of_adjoin_eq_top (A := A) hadj h)) (by omega)
  have h := ramIndex_lt_ramIndex_pow_char p huni (pow_mem hσ (p ^ j)) hne
  rwa [← pow_mul, ← pow_succ] at h

/-! ## §4 `⊤` 判定と単調性

★原文は「`i_j := ∞ for j ≥ m`」を**規約**として置くが、`ramIndex` の定義の下では
これは**定理**である(逸脱 1)。 -/

/-- **`i_j = ⊤ ⟺ m ≤ j`** —— 原文の規約「`i_j := ∞ for j ≥ m`」の中身。 -/
theorem ramIndex_pow_pow_eq_top_iff {A B : Type*} [CommRing A] [CommRing B] [Algebra A B]
    [IsDomain B] [IsDiscreteValuationRing B] {G : Type*} [Group G] [MulSemiringAction G B]
    [SMulCommClass G A B] [FaithfulSMul G B] {p : ℕ} (hp : 1 < p) {π : B}
    (hadj : Algebra.adjoin A ({π} : Set B) = ⊤) {σ : G} {m : ℕ} (hord : orderOf σ = p ^ m)
    (j : ℕ) : ramIndex π (σ ^ p ^ j) = ⊤ ↔ m ≤ j := by
  rw [ramIndex_eq_top_iff]
  refine ⟨fun h => (pow_pow_eq_one_iff hp hord j).1
    (eq_one_of_smul_eq_of_adjoin_eq_top (A := A) hadj h), fun h => ?_⟩
  rw [(pow_pow_eq_one_iff hp hord j).2 h, one_smul]

/-- **`j ↦ i_j` は単調**(`m` 以降は `⊤` で一定)。(i) 後半が使う。 -/
theorem monotone_ramIndex_pow_pow {A B : Type*} [CommRing A] [CommRing B] [Algebra A B]
    [IsDomain B] [IsDiscreteValuationRing B] {G : Type*} [Group G] [MulSemiringAction G B]
    [SMulCommClass G A B] [FaithfulSMul G B] (p : ℕ) (hp : p.Prime)
    [CharP (ResidueField B) p] {π : B} (huni : maximalIdeal B = Ideal.span {π})
    (hadj : Algebra.adjoin A ({π} : Set B) = ⊤) {σ : G}
    (hσ : σ ∈ lowerRamificationGroup B G 1) {m : ℕ} (hord : orderOf σ = p ^ m) :
    Monotone fun j => ramIndex π (σ ^ p ^ j) := by
  refine monotone_nat_of_le_succ fun j => ?_
  rcases lt_or_ge j m with hj | hj
  · exact le_of_lt (ramIndex_pow_pow_lt_succ (A := A) p hp huni hadj hσ hord hj)
  · rw [(ramIndex_pow_pow_eq_top_iff (A := A) hp.one_lt hadj hord j).2 hj,
      (ramIndex_pow_pow_eq_top_iff (A := A) hp.one_lt hadj hord (j + 1)).2 (by omega)]

/-! ## §5 Yoshida Prop 6.6 (ii) —— `i(σ^a) = i_{v_p(a)}`

★★**原文は (ii) を (i) 経由で証明する**が、ここでは (i) を一切使わない。
`i(τ)` は `τ` の**生成する部分群だけ**で決まる(`σ ∈ G_n ⟺ ⟨σ⟩ ≤ G_n`)ので、
§1 の `⟨σ^a⟩ = ⟨σ^{p^k}⟩` から直ちに出る。 -/

/-- ★**`i` は生成部分群だけで決まる** —— `⟨τ⟩ = ⟨τ'⟩ ⟹ i(τ) = i(τ')`。

`σ ∈ G_n ⟺ ⟨σ⟩ ≤ G_n`(`Subgroup.zpowers_le`)と `G_n` 判定(`Def 6.1`)、
`ℕ∞` の外延性(§2)の 3 つだけ。 -/
theorem ramIndex_eq_of_zpowers_eq {A B : Type*} [CommRing A] [CommRing B] [Algebra A B]
    [IsDomain B] [IsDiscreteValuationRing B] {G : Type*} [Group G] [MulSemiringAction G B]
    [SMulCommClass G A B] {π : B} (huni : maximalIdeal B = Ideal.span {π})
    (hadj : Algebra.adjoin A ({π} : Set B) = ⊤) {τ τ' : G}
    (h : Subgroup.zpowers τ = Subgroup.zpowers τ') : ramIndex π τ = ramIndex π τ' := by
  refine enat_eq_of_forall_natCast_lt_iff fun n => ?_
  rw [← mem_lowerRamificationGroup_iff_lt_ramIndex (A := A) huni hadj,
    ← mem_lowerRamificationGroup_iff_lt_ramIndex (A := A) huni hadj,
    ← Subgroup.zpowers_le, ← Subgroup.zpowers_le, h]

/-- ★★★**Yoshida 2008 Prop 6.6 (ii)** —— `p^k ∥ a` なら `i(σ^a) = i_k`。

★原文の `v_p` を関数として使わず、整除性の述語 `p^k ∣ a ∧ ¬ p^{k+1} ∣ a` で書いた(逸脱 4)。
★`a ≥ 1` は `¬ p^{k+1} ∣ a` に含まれるので落とした(逸脱 3)。
★(i) を使っていない(逸脱 5)。 -/
theorem ramIndex_pow_eq_ramIndex_pow_pow {A B : Type*} [CommRing A] [CommRing B] [Algebra A B]
    [IsDomain B] [IsDiscreteValuationRing B] {G : Type*} [Group G] [MulSemiringAction G B]
    [SMulCommClass G A B] (p : ℕ) (hp : p.Prime) {π : B}
    (huni : maximalIdeal B = Ideal.span {π}) (hadj : Algebra.adjoin A ({π} : Set B) = ⊤)
    {σ : G} {m : ℕ} (hord : orderOf σ = p ^ m) {a k : ℕ} (hdvd : p ^ k ∣ a)
    (hndvd : ¬ p ^ (k + 1) ∣ a) : ramIndex π (σ ^ a) = ramIndex π (σ ^ p ^ k) :=
  ramIndex_eq_of_zpowers_eq (A := A) huni hadj (zpowers_pow_eq_zpowers_pow hp hord hdvd hndvd)

/-! ## §6 Yoshida Prop 6.6 (i) 後半 —— `H_n = ⟨σ^{p^j}⟩ ⟺ i_{j−1} ≤ n < i_j` -/

/-- ★★★**Yoshida 2008 Prop 6.6 (i) 後半** —— `j + 1 ≤ m` のとき

`G_n ⊓ ⟨σ⟩ = ⟨σ^{p^{j+1}}⟩ ⟺ (i_j ≤ n ∧ n < i_{j+1})`。

★添字を 1 つずらしてある(原文 `H_n = ⟨σ^{p^j}⟩ ⟺ i_{j−1} ≤ n < i_j`、`1 ≤ j ≤ m`。逸脱 2)。
★`j + 1 ≤ m` は落とせない(退化検査 D2)。

証明: (⟸) `⊇` は `n < i_{j+1}` から `σ^{p^{j+1}} ∈ G_n`。`⊆` は `x = σ^a` と書いて
`p^k ∥ a` を取り、(ii) で `i(σ^a) = i_k`、単調性(§4)と `i_j ≤ n` から `k ≥ j+1`、
よって `p^{j+1} ∣ a`。(⟹) `n < i_{j+1}` は `σ^{p^{j+1}} ∈ H_n` から、
`i_j ≤ n` は狭義降下(§1)から。 -/
theorem lowerRamificationGroup_inf_zpowers_eq_iff {A B : Type*} [CommRing A] [CommRing B]
    [Algebra A B] [IsDomain B] [IsDiscreteValuationRing B] {G : Type*} [Group G]
    [MulSemiringAction G B] [SMulCommClass G A B] [FaithfulSMul G B] (p : ℕ) (hp : p.Prime)
    [CharP (ResidueField B) p] {π : B} (huni : maximalIdeal B = Ideal.span {π})
    (hadj : Algebra.adjoin A ({π} : Set B) = ⊤) {σ : G}
    (hσ : σ ∈ lowerRamificationGroup B G 1) {m : ℕ} (hord : orderOf σ = p ^ m)
    {j : ℕ} (hj : j + 1 ≤ m) (n : ℕ) :
    lowerRamificationGroup B G n ⊓ Subgroup.zpowers σ = Subgroup.zpowers (σ ^ p ^ (j + 1)) ↔
      ramIndex π (σ ^ p ^ j) ≤ (n : ℕ∞) ∧ (n : ℕ∞) < ramIndex π (σ ^ p ^ (j + 1)) := by
  constructor
  · intro h
    refine ⟨?_, ?_⟩
    · by_contra hlt
      rw [not_le] at hlt
      have hmem : σ ^ p ^ j ∈ lowerRamificationGroup B G n ⊓ Subgroup.zpowers σ :=
        Subgroup.mem_inf.2
          ⟨(mem_lowerRamificationGroup_iff_lt_ramIndex (A := A) huni hadj n _).2 hlt,
            pow_mem (Subgroup.mem_zpowers σ) _⟩
      rw [h] at hmem
      exact pow_pow_notMem_zpowers_pow_succ hp hord (by omega) hmem
    · have hmem : σ ^ p ^ (j + 1) ∈ lowerRamificationGroup B G n ⊓ Subgroup.zpowers σ := by
        rw [h]; exact Subgroup.mem_zpowers _
      exact (mem_lowerRamificationGroup_iff_lt_ramIndex (A := A) huni hadj n _).1
        (Subgroup.mem_inf.1 hmem).1
  · rintro ⟨h1, h2⟩
    apply le_antisymm
    · intro x hx
      obtain ⟨hxn, hxz⟩ := Subgroup.mem_inf.1 hx
      obtain ⟨a, rfl⟩ := exists_pow_eq_of_mem_zpowers hp.pos hord hxz
      rcases eq_or_ne a 0 with rfl | ha
      · simp
      obtain ⟨k, u, hu, hau⟩ := Nat.exists_eq_pow_mul_and_not_dvd ha p hp.ne_one
      subst hau
      have hndvd : ¬ p ^ (k + 1) ∣ p ^ k * u := by
        intro hc
        rw [pow_succ] at hc
        exact hu ((mul_dvd_mul_iff_left (pow_ne_zero k hp.pos.ne')).1 hc)
      have hik := ramIndex_pow_eq_ramIndex_pow_pow (A := A) p hp huni hadj hord
        (Dvd.intro u rfl) hndvd
      have hn : (n : ℕ∞) < ramIndex π (σ ^ p ^ k) := by
        rw [← hik]
        exact (mem_lowerRamificationGroup_iff_lt_ramIndex (A := A) huni hadj n _).1 hxn
      have hk : j + 1 ≤ k := by
        by_contra hk
        have hmono := monotone_ramIndex_pow_pow (A := A) p hp huni hadj hσ hord
          (show k ≤ j by omega)
        exact absurd (lt_of_lt_of_le hn hmono) (not_lt.2 h1)
      obtain ⟨c, hc⟩ : p ^ (j + 1) ∣ p ^ k * u := Dvd.dvd.mul_right (pow_dvd_pow p hk) u
      rw [hc, pow_mul]
      exact pow_mem (Subgroup.mem_zpowers _) c
    · exact Subgroup.zpowers_le.2 (Subgroup.mem_inf.2
        ⟨(mem_lowerRamificationGroup_iff_lt_ramIndex (A := A) huni hadj n _).2 h2,
          pow_mem (Subgroup.mem_zpowers σ) _⟩)

/-! ## §7 `PAdicLocalField` への具体化

`B := adjoinIntegers K x`、`A := 𝒪[K.carrier]`、`G := Gal(K(x)/K)`。
★`IsDiscreteValuationRing.addVal` を statement に出すので DVR を
`attribute [local instance]` で入れる(`lean-idioms.md` #85)。
★`orderOf σ = p^m` は原文の "(by Proposition 6.2)" で、Y4 が供給する。 -/

variable {p : ℕ} [Fact p.Prime]

attribute [local instance] isDiscreteValuationRing_adjoinIntegers

/-- ★★★**Yoshida 2008 Prop 6.6 (i) 前半(`PAdicLocalField` 版)** ——
`σ ∈ G_1` に対し `|⟨σ⟩| = p^m` となる `m` があり、`j < m` の範囲で `i_j < i_{j+1}`。

★`orderOf σ = p^m` は原文の "(by Proposition 6.2)"(`Found/PGC/RamificationJumpDivisibility`)。 -/
theorem exists_orderOf_eq_pow_and_ramIndex_pow_pow_lt_succ_adjoin (K : PAdicLocalField p)
    (x : K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (ht : IsTotallyRamifiedAdjoin K x) {π : adjoinIntegers K x} (hπ : Irreducible π)
    {σ : (IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure))}
    (hσ : σ ∈ lowerRamificationGroupAdjoin K x 1) :
    ∃ m : ℕ, orderOf σ = p ^ m ∧
      ∀ j : ℕ, j < m → ramIndex π (σ ^ p ^ j) < ramIndex π (σ ^ p ^ (j + 1)) := by
  haveI := charP_residueField_adjoinIntegers K x
  have huni := (irreducible_iff_uniformizer π).mp hπ
  have hadj := adjoin_uniformizer_eq_top_adjoinIntegers K x ht huni
  obtain ⟨m, hm⟩ := exists_orderOf_eq_pow_of_mem_lowerRamificationGroup_one
    (A := 𝒪[K.carrier]) p huni hπ.ne_zero hadj hσ
  exact ⟨m, hm, fun j hj => ramIndex_pow_pow_lt_succ (A := 𝒪[K.carrier]) p (Fact.out) huni hadj
    hσ hm hj⟩

/-- ★★★**Yoshida 2008 Prop 6.6 (ii)(`PAdicLocalField` 版)** —— `p^k ∥ a` なら
`i(σ^a) = i_k`。 -/
theorem ramIndex_pow_eq_ramIndex_pow_pow_adjoin (K : PAdicLocalField p) (x : K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (ht : IsTotallyRamifiedAdjoin K x) {π : adjoinIntegers K x} (hπ : Irreducible π)
    {σ : (IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure))}
    (hσ : σ ∈ lowerRamificationGroupAdjoin K x 1) {a k : ℕ} (hdvd : p ^ k ∣ a)
    (hndvd : ¬ p ^ (k + 1) ∣ a) : ramIndex π (σ ^ a) = ramIndex π (σ ^ p ^ k) := by
  haveI := charP_residueField_adjoinIntegers K x
  have huni := (irreducible_iff_uniformizer π).mp hπ
  have hadj := adjoin_uniformizer_eq_top_adjoinIntegers K x ht huni
  obtain ⟨m, hm⟩ := exists_orderOf_eq_pow_of_mem_lowerRamificationGroup_one
    (A := 𝒪[K.carrier]) p huni hπ.ne_zero hadj hσ
  exact ramIndex_pow_eq_ramIndex_pow_pow (A := 𝒪[K.carrier]) p (Fact.out) huni hadj hm hdvd hndvd

/-- ★★★**Yoshida 2008 Prop 6.6 (i) 後半(`PAdicLocalField` 版)** ——
`j + 1 ≤ m` のとき `H_n = G_n ⊓ ⟨σ⟩` が `⟨σ^{p^{j+1}}⟩` に一致するのは
`i_j ≤ n < i_{j+1}` のときに限る。 -/
theorem lowerRamificationGroupAdjoin_inf_zpowers_eq_iff (K : PAdicLocalField p) (x : K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (ht : IsTotallyRamifiedAdjoin K x) {π : adjoinIntegers K x} (hπ : Irreducible π)
    {σ : (IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure))}
    (hσ : σ ∈ lowerRamificationGroupAdjoin K x 1) {m : ℕ} (hord : orderOf σ = p ^ m)
    {j : ℕ} (hj : j + 1 ≤ m) (n : ℕ) :
    lowerRamificationGroupAdjoin K x n ⊓ Subgroup.zpowers σ =
        Subgroup.zpowers (σ ^ p ^ (j + 1)) ↔
      ramIndex π (σ ^ p ^ j) ≤ (n : ℕ∞) ∧ (n : ℕ∞) < ramIndex π (σ ^ p ^ (j + 1)) := by
  haveI := charP_residueField_adjoinIntegers K x
  have huni := (irreducible_iff_uniformizer π).mp hπ
  have hadj := adjoin_uniformizer_eq_top_adjoinIntegers K x ht huni
  exact lowerRamificationGroup_inf_zpowers_eq_iff (A := 𝒪[K.carrier]) p (Fact.out) huni hadj
    hσ hord hj n

end ABC3.Found.PGC
