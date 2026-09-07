import ABC3.Found.PGC.SubgroupActionBridge

/-!
# 順分岐商 `K′^{G_1}/K` の塔 —— Hasse-Arf(Theorem 6.11)の段 3 の最後の 2 本

★★★**本ノードは Yoshida 2008 Theorem 6.11 の「段 3 の塔」のみを担当する。**
段 1(`G` 巡回)は `Found/PGC/HasseArf.lean`
(`exists_natCast_herbrandPhiGroup_of_isCyclic`)、
段 2(`G = G_1`)は `Found/PGC/HasseArfStrongInduction.lean`
(`exists_natCast_herbrandPhiGroup_of_lowerRamificationGroup_one_eq_top`)、
3 段の合流は `Found/PGC/SubgroupActionBridge.lean`
(`exists_natCast_herbrandPhiGroup_of_abelian`)にある。
★本ファイルは**既存の `Found/PGC/*.lean` を 1 行も書き換えない**(import のみ)。

## ★★★到達点 —— Hasse-Arf が仮定ゼロで閉じた

先行ノード(Y22)の `exists_natCast_herbrandPhiGroup_of_abelian` は、`G ≠ G_1` の
分岐でだけ要求される仮定 `hC` を持っていた。その中身は順分岐商
`K″ := K′^{G_1}` に対応する DVR `C` と素元 `ϖ` についての 3 条件である。

| `hC` の成分 | 供給 |
|---|---|
| `htopC` : `σ ∈ G_1 ⟹ i_ϖ(σ) = ⊤` | Y22 `ramIndex_fixedRing_eq_top_of_mem`(供給済み) |
| `hcomp` : `φ_G(n) = φ_{G/G_1}(φ_{G_1}(n))` | ★**本ファイル §4** `herbrandPhiGroup_comp_fixedRing` |
| `honeC` : `σ ∉ G_1 ⟹ i_ϖ(σ) = 1` | ★**本ファイル §3** `ramIndex_fixedRing_eq_one_of_notMem` |

★★これで `C := ↥(fixedRing B G_1)` を代入した
`exists_natCast_herbrandPhiGroup_of_abelian_of_setup`(§5)からは **`hC` が消え**、
残る仮定は原典 §6.1–§6.2 の底の設定(`π′` は `𝒪_{K′}` の素元、`𝒪_{K′} = 𝒪_K[π′]`、
`𝒪_K → 𝒪_{K′}` 単射、`K′/K` 完全分岐、剰余標数 `p`)と `G` が可換であることだけになる。

## ★★設計 —— まず抽象核

`honeC` は原典が

> Now when `G ≠ G_1`, set `H = G_1` and `|G/H| = e_0`.
> As `φ_{G/H}(n) = n/e_0` for `n ∈ R≥0` by definition, …

と 1 語(`by definition`)で畳んだところである。`φ_{G/H}(n) = n/e_0` が「定義から」
言えるのは、`φ` の定義の中の `min{i(τ), n+1}` が `τ ≠ 1` でつねに `1` になるから、
すなわち**順分岐商の第 1 分岐群が自明**だからである。それを次の 4 段に切った。

| § | 宣言 | 語彙 |
|---|---|---|
| §1 | `eq_bot_of_isPGroup_of_not_dvd_natCard` | ★**純群論**。分岐・付値・Galois が 1 つも出てこない |
| §1 | `normal_of_comm` | ★**純群論**。可換群の部分群は正規 |
| §2 | `lowerRamificationGroup_one_eq_bot_of_not_dvd` | 一般の DVR + 一般の有限群作用。★`p ∤ |Q|` だけで `Q_1 = ⊥` |
| §2 | `ramIndex_eq_one_of_lowerRamificationGroup_one_eq_bot` | 同上。★`0 < i(τ)` と `¬(1 < i(τ))` を挟むだけ |
| §3–§5 | 具体層 | `C := ↥(fixedRing B G_1)`、`Q := G ⧸ G_1` を代入するだけ |

★**§2 の 2 本には `fixedRing` も `Gal` も出てこない。** 「順分岐なら第 1 分岐群は自明」は
原典の設定に依らない一般論であり、代入するだけで `honeC` になる。

★★**Y22 の見立て(「`htopC` と同じ道が `honeC` にも通るか」)は外れた。**
`htopC` は `ramIndex_eq_top_iff` + `mem_fixedRing` を 1 回ずつで、分岐の議論ゼロだった。
`honeC` はそうならない —— **第 1 分岐群が `p` 群であること**(Corollary 6.3 の中核)を
経由する必要がある。★逆に、`p` 群であることさえ使えば残りは 3 行である。

## ★★在庫は mathlib も先に引いた

★Y21 が「無い」と報告した 4 本を Y22 が mathlib で見つけた直後なので、本ノードでも
**新しい `instance` を書く前に mathlib を引いた**。使ったのは次のとおりで、
★**本ファイルが新たに定義する `instance` / `def` は 1 つも無い**。

| 要るもの | mathlib の宣言 |
|---|---|
| `IsPGroup p N → ∃ k, Nat.card N = p ^ k` | `IsPGroup.exists_card_eq` |
| `Nat.card ↥N ∣ Nat.card Q` | `Subgroup.card_subgroup_dvd_card` |
| `Nat.card ↥N = 1 ↔ N = ⊥` | `Subgroup.card_eq_one` |
| `1 ≤ x ↔ x ≠ 0`(`ℕ∞`) | `Order.one_le_iff_ne_zero` |
| `mk σ = 1 ↔ σ ∈ H` | `QuotientGroup.eq_one_iff` |
| DVR に素元がある | `IsDiscreteValuationRing.exists_irreducible` |
| `Irreducible ϖ ↔ 𝔪 = (ϖ)` | `IsDiscreteValuationRing.irreducible_iff_uniformizer` |

★`MulSemiringAction G ↥(fixedRing B H)` / `MulSemiringAction (G ⧸ H) ↥(fixedRing B H)` /
`SMulCommClass (G ⧸ H) A ↥(fixedRing B H)` / `FaithfulSMul (G ⧸ H) ↥(fixedRing B H)` /
`Algebra A ↥(fixedRing B H)` / `CharP (ResidueField ↥(fixedRing B H)) p` は
**すべて木(Y14・Y15b・Y16・Y17・Y20)にある**ので、本ファイルは `letI` / `haveI` で
呼び出すだけである。

## ★7 つの適合条件(Lemma 6.10(ii) の `hcomp`)

`herbrandPhiGroup_comp`(Y12)が要求する 10 個の仮定は、`C := ↥(fixedRing B H)` に
対して**そのまま在庫でそろった**(§4)。

| 仮定 | 供給 | 出所 |
|---|---|---|
| `hcomp`(包含の同変性) | `algebraMap_smul_fixedRing` | Y16(★`rfl`) |
| `hHtriv`(`H` は `C` に自明作用) | `smul_fixedRing_eq_self` | Y16 |
| `hinj` | `fixedRing_injective` | Y14 |
| `hfixC` | `exists_algebraMap_fixedRing` | Y14 |
| `hres`(`K′/K″` 完全分岐) | `exists_sub_mem_fixedRing` | Y14 |
| `hAC`(`𝒪_K ⊆ 𝒪_{K″}`) | `exists_algebraMap_fixedRing_eq` | Y14 |
| `hϖ`(`ι ϖ = π″`) | ★**`rfl`** | —— |
| `hfix`(`𝒪_{K′}^H ⊆ 𝒪_K[π″]`) | `fixedRing_mem_adjoin_uniformizer` | Y17 |
| `[MulSemiringAction G C]` | `fixedRingMulSemiringAction` | Y16(★`[H.Normal]`) |
| `[IsDiscreteValuationRing C]` | `fixedRing_isDiscreteValuationRing` | Y14(★`[Fintype ↥H]`) |

★**段取りより安かった。** §4 は `herbrandPhiGroup_comp` への 1 回の代入(項 1 本)である。

## ★★逸脱の記録

1. **`G` の可換性は `∀ x y : G, x * y = y * x` という命題で渡す**(`CommGroup` 構造を
   要求しない)。先行ノード(Y13・Y15・Y15b・Y21・Y22)の逸脱をそのまま引き継いだ。
   具体層の `G = Gal(K′/K)` は `Group` インスタンスしか持たないからである。
   ★§5 の中でだけ `haveI : (lowerRamificationGroup B G 1).Normal := normal_of_comm habel _`
   を立てている(`G` 自身のインスタンスは触らない)。
2. ★**原典は `e_0` と `|H|` が互いに素であることを使うが、本ファイルが実際に使うのは
   `p ∤ e_0` の方だけである**(`not_dvd_natCard_quot_lowerRamificationGroup_one`)。
   `|H| = |G_1|` が `p` 冪であること(`isPGroup_lowerRamificationGroup_one`)は
   §2 の抽象核の中で使う。★互除性そのものは経由していない(強めても弱めてもいない)。
3. ★**原典の `φ_{G/H}(n) = n/e_0` という等式そのものは述べていない。**
   段 3(`exists_natCast_herbrandPhiGroup_of_tame_quotient`、Y15)が要求するのは
   `htopC` / `honeC` の 2 条件であり、そこから `herbrandPhiGroup_eq_div_natCard_quot`
   が等式を作る。★本ファイルは**その 2 条件を供給する**という形で原典に対応する。
4. §2 の 2 本は `n = 1` に特化していない一般形ではなく、**`Q_1 = ⊥` の形**で述べてある。
   原典が必要としているのはこの形だけである。

## 退化の自己検査

* ★★**`G ≠ G_1` を落とすと段 3 は成り立たない。** §5 は
  `exists_natCast_herbrandPhiGroup_of_abelian`(Y22)の `by_cases` に乗っており、
  `G = G_1` の分岐は段 2 が処理する。★**本ファイルの §3・§4 は `G ≠ G_1` を仮定して
  いない**(`G = G_1` のときは `honeC` が空虚に真、`hcomp` は両辺とも `φ_G` になる)。
  それでも合流が正しいのは、`hC` が `G ≠ G_1` のときにだけ消費されるからである。
* ★★**`p ∤ e_0`(= `¬ p ∣ |G/G_1|`)を落とすと §2 が出ない。** これは
  `G_0/G_1 ↪ k^×`(Proposition 6.2)から来ており、`h0 : G_0 = ⊤`(`K′/K` 完全分岐)が
  無いと成り立たない。★`h0` を落とすと `honeC` は偽になりうる。
* ★★**`G_1` が `p` 群であること**(`isPGroup_lowerRamificationGroup_one`)を落とすと
  §2 の `eq_bot_of_isPGroup_of_not_dvd_natCard` に渡すものが無い。これは
  `[CharP (ResidueField C) p]`(Y20 `charP_residueField_fixedRing`)に依る。
* ★**`[H.Normal]`(`H = G_1`)は可換性から出る**(§1 `normal_of_comm`。実測で確認した)。
  落とすと `MulSemiringAction G ↥(fixedRing B H)` が立たず、`ramIndex ϖ σ` すら書けない
  (`lean-idioms.md` #125(iii))。
* ★**`[Fintype ↥H]` を落とすと `C = 𝒪_{K′}^H` が DVR にならない**(Y14 の実測)。
  `IsDiscreteValuationRing ↥(fixedRing B H)` が立たなければ `ramIndex ϖ` が定義できない。
* ★**`FaithfulSMul (G ⧸ H) C` を落とすと §2 が出ない** —— `i_ϖ(τ) = ⊤` が `τ = 1` を
  意味しなくなり、`(G/H)_1` が `p` 群であることを言う `exists_lowerRamificationGroup_eq_bot`
  が壊れる。供給は Y15b `faithfulSMul_quotient_fixedRing`(★`hfix` の下でのみ言える)。
* ★`ℕ∞` の切り詰め引き算・除算は 1 つも書いていない(`lean-idioms.md` #102)。
  除算は `phiOf` の定義の中にしか無い。§2 の `ℕ∞` の操作は `0 < x` と `¬(1 < x)` から
  `x = 1` を出す `le_antisymm` 1 回だけである。
-/

namespace ABC3.Found.PGC

open IsLocalRing IsDiscreteValuationRing

def ramIndex_fixedRing_eq_one_of_notMem.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 16, item := "Theorem 6.11", sectionId := "thm-6-11" }

def herbrandPhiGroup_comp_fixedRing.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 16, item := "Theorem 6.11", sectionId := "thm-6-11" }

def exists_natCast_herbrandPhiGroup_of_abelian_of_setup.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 16, item := "Theorem 6.11", sectionId := "thm-6-11" }

/-! ## §1 抽象核 —— 純群論

★分岐・付値・Galois の語彙が 1 つも出てこない。 -/

namespace AbstractCore

/-- ★★**抽象核** —— 位数が `p` で割れない有限群の中では、`p` 部分群は自明である。

`|N| = p^k`(`IsPGroup.exists_card_eq`)と `|N| ∣ |Q|`(Lagrange)から、
`k ≥ 1` なら `p ∣ |Q|` になってしまう。

★これが原典の "As `e_0` and `|H|` are coprime" の実質であり、
本ノードでは「順分岐商 `G/G_1` の中の `p` 部分群 `(G/G_1)_1` は自明」に使う。 -/
theorem eq_bot_of_isPGroup_of_not_dvd_natCard {Q : Type*} [Group Q] [Finite Q] {p : ℕ}
    [Fact p.Prime] {N : Subgroup Q} (hN : IsPGroup p N) (hp : ¬ p ∣ Nat.card Q) :
    N = ⊥ := by
  obtain ⟨k, hk⟩ := hN.exists_card_eq
  have hdvd : Nat.card N ∣ Nat.card Q := Subgroup.card_subgroup_dvd_card N
  rcases Nat.eq_zero_or_pos k with hk0 | hk0
  · rw [← Subgroup.card_eq_one, hk, hk0, pow_zero]
  · exact absurd ((hk ▸ dvd_pow_self p hk0.ne' : p ∣ Nat.card N).trans hdvd) hp

/-- ★**抽象核** —— 可換群の部分群は正規。

★`G` の可換性は `CommGroup` インスタンスではなく命題 `∀ x y, x * y = y * x` で
持ち回っている(逸脱 1)ので、mathlib の `Subgroup.normal_of_comm` は直接使えない。 -/
theorem normal_of_comm {Q : Type*} [Group Q] (habel : ∀ x y : Q, x * y = y * x)
    (N : Subgroup Q) : N.Normal :=
  ⟨fun a ha g => by rwa [habel g a, mul_assoc, mul_inv_cancel, mul_one]⟩

end AbstractCore

open AbstractCore

/-! ## §2 抽象核 —— 順分岐なら第 1 分岐群は自明

★**一般の DVR `C` と一般の有限群作用**についての主張である。
`fixedRing` も `Gal` も商群も出てこない。原典の設定に代入するのは §3 である。 -/

/-- ★★★**抽象核** —— `p ∤ |Q|`(順分岐)なら `Q_1 = ⊥`。

`Q_1` は `p` 群(`isPGroup_lowerRamificationGroup_one`、Corollary 6.3 の中核)であり、
`|Q|` は `p` で割れないので §1 の抽象核が効く。

★★これが原典の `φ_{G/H}(n) = n/e_0` の "by definition" の中身である。
`φ` の定義の `min{i(τ), n+1}` が `τ ≠ 1` でつねに `1` になるのは、
`i(τ) ≥ 2 ⟺ τ ∈ Q_1 = ⊥` だからである。

★**退化**: `[FaithfulSMul Q C]` を落とすと `Q_1` が `p` 群であることが言えない。 -/
theorem lowerRamificationGroup_one_eq_bot_of_not_dvd
    {A C : Type*} [CommRing A] [CommRing C] [Algebra A C] [IsDomain C]
    [IsDiscreteValuationRing C] {Q : Type*} [Group Q] [MulSemiringAction Q C]
    [SMulCommClass Q A C] [Finite Q] [FaithfulSMul Q C]
    (p : ℕ) [Fact p.Prime] [CharP (ResidueField C) p] {ϖ : C}
    (huni : maximalIdeal C = Ideal.span {ϖ}) (hϖ0 : ϖ ≠ 0)
    (hadj : Algebra.adjoin A ({ϖ} : Set C) = ⊤)
    (hnd : ¬ p ∣ Nat.card Q) :
    lowerRamificationGroup C Q 1 = ⊥ :=
  eq_bot_of_isPGroup_of_not_dvd_natCard
    (isPGroup_lowerRamificationGroup_one (A := A) p huni hϖ0 hadj) hnd

/-- ★★**抽象核** —— `Q_1 = ⊥` なら `τ ≠ 1` に対し `i_ϖ(τ) = 1`。

`Definition 6.1` の形(`τ ∈ Q_n ⟺ i(τ) > n`)を `n = 1` で使う:

* 上から: `τ ∉ Q_1` なので `¬ (1 < i(τ))`。
* 下から: `i(τ) ≥ 1`(`pos_ramIndex`。`ϖ` が素元なので `τϖ − ϖ ∈ 𝔪_C`)。

★`ℕ∞` の操作はこの `le_antisymm` 1 回だけである(切り詰め引き算・除算は書かない)。 -/
theorem ramIndex_eq_one_of_lowerRamificationGroup_one_eq_bot
    {A C : Type*} [CommRing A] [CommRing C] [Algebra A C] [IsDomain C]
    [IsDiscreteValuationRing C] {Q : Type*} [Group Q] [MulSemiringAction Q C]
    [SMulCommClass Q A C] {ϖ : C}
    (huni : maximalIdeal C = Ideal.span {ϖ})
    (hadj : Algebra.adjoin A ({ϖ} : Set C) = ⊤)
    (hbot : lowerRamificationGroup C Q 1 = ⊥) {τ : Q} (hτ : τ ≠ 1) :
    ramIndex ϖ τ = 1 := by
  have hnot : τ ∉ lowerRamificationGroup C Q 1 := by rw [hbot]; simpa using hτ
  rw [mem_lowerRamificationGroup_iff_lt_addVal (A := A) huni hadj 1 τ] at hnot
  exact le_antisymm (by exact_mod_cast not_lt.mp hnot)
    (Order.one_le_iff_ne_zero.mpr (pos_ramIndex huni τ).ne')

/-- ★★★**抽象核(合流)** —— 順分岐(`p ∤ |Q|`)な完全分岐拡大では
`τ ≠ 1 ⟹ i_ϖ(τ) = 1`。

★上の 2 本をつないだだけ。★**これが `honeC` の数学的中身である。** -/
theorem ramIndex_eq_one_of_ne_one_of_not_dvd
    {A C : Type*} [CommRing A] [CommRing C] [Algebra A C] [IsDomain C]
    [IsDiscreteValuationRing C] {Q : Type*} [Group Q] [MulSemiringAction Q C]
    [SMulCommClass Q A C] [Finite Q] [FaithfulSMul Q C]
    (p : ℕ) [Fact p.Prime] [CharP (ResidueField C) p] {ϖ : C}
    (huni : maximalIdeal C = Ideal.span {ϖ}) (hϖ0 : ϖ ≠ 0)
    (hadj : Algebra.adjoin A ({ϖ} : Set C) = ⊤)
    (hnd : ¬ p ∣ Nat.card Q) {τ : Q} (hτ : τ ≠ 1) :
    ramIndex ϖ τ = 1 :=
  ramIndex_eq_one_of_lowerRamificationGroup_one_eq_bot (A := A) huni hadj
    (lowerRamificationGroup_one_eq_bot_of_not_dvd (A := A) p huni hϖ0 hadj hnd) hτ

/-! ## §3 具体層 —— `honeC`(順分岐商の第 1 分岐群が自明)

★`C := ↥(fixedRing B G_1)`、`Q := G ⧸ G_1` を §2 に代入するだけ。 -/

/-- ★★★**順分岐商 `Gal(K″/K) = G/G_1` の第 1 分岐群は自明**(`K″ = K′^{G_1}`)。

原典の証明段落（`0_Source` の `.txt`、PyMuPDF 抽出。`̸=` は `≠`）が
`by definition` の 1 語で畳んだところ:
> Now when G ̸= G1, set H = G1 and |G/H| = e0. As φG/H(n) = n/e0 for n ∈R≥0 by definition,

★`φ_{G/H}(n) = n/e_0` が「定義から」出るのは、まさにこの
`(G/G_1)_1 = 1` があるからである。

`p ∤ |G/G_1|` は Corollary 6.3 の `not_dvd_natCard_quot_lowerRamificationGroup_one`
(`G_0/G_1 ↪ k^×` から)が与える。★そこに `h0 : G_0 = ⊤`(`K′/K` 完全分岐)が要る。 -/
theorem lowerRamificationGroup_quotient_fixedRing_one_eq_bot
    {A B : Type*} [CommRing A] [IsLocalRing A] [CommRing B] [Algebra A B]
    [IsDomain B] [IsDiscreteValuationRing B] [IsNoetherian A B]
    {G : Type*} [Group G] [MulSemiringAction G B] [Fintype G] [SMulCommClass G A B]
    [FaithfulSMul G B] [Fintype ↥(lowerRamificationGroup B G 1)]
    [(lowerRamificationGroup B G 1).Normal]
    [MulSemiringAction (G ⧸ lowerRamificationGroup B G 1)
      ↥(fixedRing B (lowerRamificationGroup B G 1))]
    (hq : ∀ (σ : G) (c : ↥(fixedRing B (lowerRamificationGroup B G 1))),
      (QuotientGroup.mk σ : G ⧸ lowerRamificationGroup B G 1) • c = σ • c)
    (p : ℕ) (hp : p.Prime) [CharP (ResidueField B) p]
    {π' : B} (hπ' : Irreducible π')
    (hresA : ∀ b : B, ∃ a : A, b - algebraMap A B a ∈ maximalIdeal B)
    (hadj : Algebra.adjoin A ({π'} : Set B) = ⊤)
    (hAne : ∃ a ∈ maximalIdeal A, algebraMap A B a ≠ 0)
    (h0 : lowerRamificationGroup B G 0 = ⊤)
    {ϖ : ↥(fixedRing B (lowerRamificationGroup B G 1))} (hϖ : Irreducible ϖ) :
    lowerRamificationGroup ↥(fixedRing B (lowerRamificationGroup B G 1))
      (G ⧸ lowerRamificationGroup B G 1) 1 = ⊥ := by
  haveI : Fact p.Prime := ⟨hp⟩
  letI := fixedRingAlgebra A B (lowerRamificationGroup B G 1)
  haveI := smulCommClass_quotient_fixedRing (A := A) (B := B) (G := G)
    (H := lowerRamificationGroup B G 1) hq
  haveI := charP_residueField_fixedRing (B := B) (H := lowerRamificationGroup B G 1) p hp
  haveI := faithfulSMul_quotient_fixedRing (A := A) hq hπ' hresA hadj
    (fixedRing_mem_adjoin_uniformizer hresA hAne hϖ)
  exact lowerRamificationGroup_one_eq_bot_of_not_dvd (A := A) p
    ((irreducible_iff_uniformizer ϖ).1 hϖ) hϖ.ne_zero
    (adjoin_fixedRing_uniformizer_eq_top hresA hAne hϖ)
    (not_dvd_natCard_quot_lowerRamificationGroup_one (A := A) p
      ((irreducible_iff_uniformizer π').1 hπ') hπ'.ne_zero hadj h0)

/-- ★★★★**段 3 の `honeC` の供給** —— `σ ∉ G_1 ⟹ i_ϖ(σ) = 1`。

`ϖ` は `𝒪_{K″} = 𝒪_{K′}^{G_1}` の素元、`i_ϖ` はその上での分岐指数である。

★`i_ϖ` は `G_1` を経由する(`ramIndex_quotient_mk`)ので、上の
`(G/G_1)_1 = ⊥` と §2 の抽象核から直ちに出る。

★★**Y22 の見立て(「`htopC` と同じく分岐の議論ゼロで出るか」)は外れた。**
`htopC` は固定環の定義だけで出たが、`honeC` は `(G/G_1)_1` が `p` 群であること
(Corollary 6.3 の中核)を経由する必要がある。 -/
theorem ramIndex_fixedRing_eq_one_of_notMem
    {A B : Type*} [CommRing A] [IsLocalRing A] [CommRing B] [Algebra A B]
    [IsDomain B] [IsDiscreteValuationRing B] [IsNoetherian A B]
    {G : Type*} [Group G] [MulSemiringAction G B] [Fintype G] [SMulCommClass G A B]
    [FaithfulSMul G B] [Fintype ↥(lowerRamificationGroup B G 1)]
    [(lowerRamificationGroup B G 1).Normal]
    (p : ℕ) (hp : p.Prime) [CharP (ResidueField B) p]
    {π' : B} (hπ' : Irreducible π')
    (hresA : ∀ b : B, ∃ a : A, b - algebraMap A B a ∈ maximalIdeal B)
    (hadj : Algebra.adjoin A ({π'} : Set B) = ⊤)
    (hAne : ∃ a ∈ maximalIdeal A, algebraMap A B a ≠ 0)
    (h0 : lowerRamificationGroup B G 0 = ⊤)
    {ϖ : ↥(fixedRing B (lowerRamificationGroup B G 1))} (hϖ : Irreducible ϖ)
    {σ : G} (hσ : σ ∉ lowerRamificationGroup B G 1) :
    ramIndex ϖ σ = 1 := by
  haveI : Fact p.Prime := ⟨hp⟩
  letI := fixedRingAlgebra A B (lowerRamificationGroup B G 1)
  letI := quotientMulSemiringActionOfTrivial (↥(fixedRing B (lowerRamificationGroup B G 1)))
    (smul_fixedRing_eq_self (B := B) (H := lowerRamificationGroup B G 1))
  haveI := smulCommClass_quotient_fixedRing (A := A) (B := B) (G := G)
    (H := lowerRamificationGroup B G 1) (fun σ c => quotientSMul_mk_fixedRing σ c)
  rw [← ramIndex_quotient_mk (H := lowerRamificationGroup B G 1)
    (fun σ c => quotientSMul_mk_fixedRing σ c) ϖ σ]
  exact ramIndex_eq_one_of_lowerRamificationGroup_one_eq_bot (A := A)
    ((irreducible_iff_uniformizer ϖ).1 hϖ)
    (adjoin_fixedRing_uniformizer_eq_top hresA hAne hϖ)
    (lowerRamificationGroup_quotient_fixedRing_one_eq_bot (A := A)
      (fun σ c => quotientSMul_mk_fixedRing σ c) p hp hπ' hresA hadj hAne h0 hϖ)
    (fun hc => hσ ((QuotientGroup.eq_one_iff σ).1 hc))

/-! ## §4 具体層 —— `hcomp`(Lemma 6.10(ii) の 7 つの適合条件)

★Y12 の `herbrandPhiGroup_comp` に `C := ↥(fixedRing B H)` を代入するだけ。
10 個の仮定はすべて Y14・Y16・Y17 の在庫でそろった(ファイル冒頭の表)。 -/

/-- ★★★★**段 3 の `hcomp` の供給** —— `φ_G = φ_{G/H} ∘ φ_H`(`C = 𝒪_{K′}^H`)。

原文 (Yoshida08 p.16) の Lemma 6.10(ii):
> (ii) φ_G = φ_G/H ◦ φ_H on R[bb]_≥0.

★★**7 つの適合条件はすべて在庫でそろった**(新しい補題は 1 本も要らなかった)。
`hϖ : algebraMap C B ϖ = π″` は `π″ := ι ϖ` と取って **`rfl`**、
`hcomp : ι (ρ • c) = ρ • ι c` も Y16 が **`rfl`** で出している。

★**`H` は `G_1` である必要が無い**(任意の有限正規部分群でよい)。 -/
theorem herbrandPhiGroup_comp_fixedRing
    {A B : Type*} [CommRing A] [IsLocalRing A] [CommRing B] [Algebra A B]
    [IsDomain B] [IsDiscreteValuationRing B] [IsNoetherian A B]
    {G : Type*} [Group G] [MulSemiringAction G B] [Fintype G] [SMulCommClass G A B]
    [FaithfulSMul G B] {H : Subgroup G} [H.Normal] [Fintype ↥H]
    {π' : B} (hπ' : Irreducible π')
    (hresA : ∀ b : B, ∃ a : A, b - algebraMap A B a ∈ maximalIdeal B)
    (hadj : Algebra.adjoin A ({π'} : Set B) = ⊤)
    (hAne : ∃ a ∈ maximalIdeal A, algebraMap A B a ≠ 0)
    {ϖ : ↥(fixedRing B H)} (hϖ : Irreducible ϖ) (n : ℝ) :
    herbrandPhiGroup G π' n = herbrandPhiGroup G ϖ (herbrandPhi π' H n) :=
  herbrandPhiGroup_comp (A := A) (fun ρ c => algebraMap_smul_fixedRing ρ c)
    smul_fixedRing_eq_self hπ' fixedRing_injective exists_algebraMap_fixedRing
    (exists_sub_mem_fixedRing hresA) exists_algebraMap_fixedRing_eq rfl hadj
    (fixedRing_mem_adjoin_uniformizer hresA hAne hϖ) n

/-! ## §5 ★★★★★ Hasse-Arf(Theorem 6.11)—— 仮定ゼロ -/

/-- ★★★★★**Yoshida 2008 Theorem 6.11 (Hasse-Arf) —— `hC` が消えた形**。

原文 (Yoshida08 p.16):
> Theorem 6.11 (Hasse-Arf). If G is abelian, n ∈ Z[bb]_≥0 and G_n = G_n+1, then φ_G(n) ∈ Z[bb]_≥0.

★逐語の `G_n = G_n+1` は `pdftotext` が `≠` の斜線を落とした形である
(原典は `G_n ≠ G_{n+1}`)。

Y22 の `exists_natCast_herbrandPhiGroup_of_abelian` は、`G ≠ G_1` の分岐でだけ
要求される仮定 `hC`(順分岐商の塔に関する 3 条件)を持っていた。
★★**本定理はそれを `C := ↥(fixedRing B G_1)` と `ϖ :=`(その素元)で埋めたものであり、
残る仮定は原典 §6.1–§6.2 の底の設定と `G` の可換性だけである**:

* `hπ'` : `π′` は `𝒪_{K′}` の素元。
* `hresA` : ★`K′/K` が完全分岐(剰余体が伸びない)。
* `hadj` : `𝒪_{K′} = 𝒪_K[π′]`(原典 Lemma 5.11)。
* `hAinj` : `𝒪_K → 𝒪_{K′}` が単射(体の拡大なので原典では自明)。
* `h0` : `G_0 = G`(`K′/K` が完全分岐)。
* `habel` : ★原文の `If G is abelian`(逸脱 1 の形)。
* `hne` : ★原文の `G_n ≠ G_{n+1}`。
* `[CharP (ResidueField B) p]` : 剰余標数 `p`。

★`ϖ` は `IsDiscreteValuationRing.exists_irreducible` が取る(選択の余地は無い ——
`ramIndex` は素元の取り方に依らない)。 -/
theorem exists_natCast_herbrandPhiGroup_of_abelian_of_setup
    {A B : Type*} [CommRing A] [IsDomain A] [IsDiscreteValuationRing A] [CommRing B] [Algebra A B]
    [IsDomain B] [IsDiscreteValuationRing B] [IsNoetherian A B]
    {G : Type*} [Group G] [MulSemiringAction G B] [Fintype G] [SMulCommClass G A B]
    [FaithfulSMul G B] [Fintype ↥(lowerRamificationGroup B G 1)]
    (p : ℕ) [hp : Fact p.Prime] [CharP (ResidueField B) p]
    {π' : B} (hπ' : Irreducible π')
    (hresA : ∀ b : B, ∃ a : A, b - algebraMap A B a ∈ maximalIdeal B)
    (hadj : Algebra.adjoin A ({π'} : Set B) = ⊤)
    (hAinj : Function.Injective (algebraMap A B))
    (h0 : lowerRamificationGroup B G 0 = ⊤) (habel : ∀ x y : G, x * y = y * x)
    {n : ℕ} (hne : lowerRamificationGroup B G n ≠ lowerRamificationGroup B G (n + 1)) :
    ∃ k : ℕ, herbrandPhiGroup G π' (n : ℝ) = (k : ℝ) := by
  haveI : (lowerRamificationGroup B G 1).Normal := normal_of_comm habel _
  obtain ⟨ϖ, hϖ⟩ :=
    IsDiscreteValuationRing.exists_irreducible ↥(fixedRing B (lowerRamificationGroup B G 1))
  have hAne : ∃ a ∈ maximalIdeal A, algebraMap A B a ≠ 0 :=
    exists_mem_maximalIdeal_map_ne_zero hAinj
  exact exists_natCast_herbrandPhiGroup_of_abelian (A := A)
    (C := ↥(fixedRing B (lowerRamificationGroup B G 1))) p hπ' hresA hadj hAinj h0 habel hne
    (ϖ := ϖ) (fun _ => ⟨herbrandPhiGroup_comp_fixedRing (A := A) hπ' hresA hadj hAne hϖ _,
      fun _ hσ => ramIndex_fixedRing_eq_top_of_mem ϖ hσ,
      fun _ hσ => ramIndex_fixedRing_eq_one_of_notMem (A := A) p hp.out hπ' hresA hadj hAne h0
        hϖ hσ⟩)

/-- ★★★★★**Theorem 6.11 の結論そのまま** —— `φ_G(n) ∈ ℤ_{≥0}` の `≥ 0` の側。

★上の `∃ k : ℕ, φ_G(n) = k` と合わせて、原文の `φ_G(n) ∈ Z[bb]_≥0` が
**仮定ゼロで**言えたことになる。 -/
theorem herbrandPhiGroup_nonneg_of_abelian_of_setup
    {A B : Type*} [CommRing A] [IsDomain A] [IsDiscreteValuationRing A] [CommRing B] [Algebra A B]
    [IsDomain B] [IsDiscreteValuationRing B] [IsNoetherian A B]
    {G : Type*} [Group G] [MulSemiringAction G B] [Fintype G] [SMulCommClass G A B]
    [FaithfulSMul G B] [Fintype ↥(lowerRamificationGroup B G 1)]
    (p : ℕ) [Fact p.Prime] [CharP (ResidueField B) p]
    {π' : B} (hπ' : Irreducible π')
    (hresA : ∀ b : B, ∃ a : A, b - algebraMap A B a ∈ maximalIdeal B)
    (hadj : Algebra.adjoin A ({π'} : Set B) = ⊤)
    (hAinj : Function.Injective (algebraMap A B))
    (h0 : lowerRamificationGroup B G 0 = ⊤) (habel : ∀ x y : G, x * y = y * x)
    {n : ℕ} (hne : lowerRamificationGroup B G n ≠ lowerRamificationGroup B G (n + 1)) :
    0 ≤ herbrandPhiGroup G π' (n : ℝ) :=
  nonneg_of_exists_natCast (exists_natCast_herbrandPhiGroup_of_abelian_of_setup (A := A) p hπ'
    hresA hadj hAinj h0 habel hne)

end ABC3.Found.PGC
