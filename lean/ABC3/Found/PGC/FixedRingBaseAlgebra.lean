import ABC3.Found.PGC.HasseArfInduction
import ABC3.Found.PGC.FixedRingMonogenic

/-!
# 固定環 `B^H` に底環 `A` の代数構造を降ろす(★原典が名前を付けていない基盤ノード)

★★**原典が名前を付けていない基盤ノードである。**
Y15b `Found/PGC/HasseArfInduction.lean` の
`exists_natCast_herbrandPhiGroup_of_cyclic_quotient_fixedRing`(Hasse-Arf 段 2)の
残り仮定 3 本を供給するために立てた。したがって本ファイルには
`.src`(`ABC3.Meta.Source`)を**置かない** —— 存在しない `sectionId` を書くことになるからである。

## ★何を埋めたか —— Y15b が名指しした 3 本

Y15b はファイル冒頭「新しく必要になったノード」でこう書いている:

1. `Algebra A ↥(fixedRing B H)`(＋ `SMulCommClass (G ⧸ H) A C`・`hresC`)
2. `CharP (ResidueField ↥(fixedRing B H)) p`
3. Y17 の `hfix` / その `K″` 側 `hadjC`

★★**3 本とも埋まった。** §3 の `exists_natCast_herbrandPhiGroup_of_cyclic_quotient_base` は
段 2 の主張を、これら 3 本と `[MulSemiringAction (G ⧸ H) C]`・`hq`・
`[FaithfulSMul (G ⧸ H) C]` を**すべて供給した形**で述べたものである。
★**残るのは `hind`(`j` に関する本当の再帰)だけ**になった。

## ★抽象核(§1)と具体層(§2)

★§1 は**分岐・付値・Galois・群作用の語彙が 1 つも出てこない**。

| §1 抽象核 | 内容 |
|---|---|
| `charP_residueField_of_algebraMap` | ★★**局所環の射 `C → B` があれば `B` の剰余体の標数が `C` の剰余体に降りる**。`(p : B) ∈ 𝔪_B` を Y17 の `mem_maximalIdeal_of_map_mem` で `(p : C) ∈ 𝔪_C` に引き戻し、`CharP.charP_iff_prime_eq_zero` を両側で使う |
| `smulCommClass_of_smul_algebraMap` | ★モノイド `M` が環 `S` に環作用し、`A` の像を各点固定するなら `SMulCommClass M A S`。`Algebra.smul_def` + `smul_mul'` の 1 行 |

| §2 具体層 | 内容 |
|---|---|
| `fixedRingAlgebra` | ★**Y17 の `Subring.algebraOfMapsTo` をそのまま再利用した**(新しく書いていない)。`algebraMap A B a ∈ fixedRing B H` は Y17 の `algebraMap_mem_fixedRing` |
| `coe_algebraMap_fixedRing` / `isScalarTower_fixedRing` | ★どちらも **`rfl`**(Y17 の実測どおり) |
| `exists_sub_mem_maximalIdeal_fixedRing` | `hresC`。Y17 の `exists_sub_mem_maximalIdeal_mid` に塔を渡すだけ |
| `charP_residueField_fixedRing` | ★`CharP (ResidueField ↥(fixedRing B H)) p`。§1 に代入するだけ |
| `smul_algebraMap_fixedRing_eq_self` | `G` は底の像を固定する(`Subtype.ext` 1 語) |
| `smulCommClass_quotient_fixedRing` | `SMulCommClass (G ⧸ H) A ↥(fixedRing B H)` |

## ★★1 を大域インスタンスにしなかった理由(実機で確かめた)

★**`letI` に留めた。** 理由は Y17(`lean-idioms.md` #119)と同じで、**測って確かめた**:

* Y15b の段 2 の**結論**は `∃ j : ℕ, herbrandPhiGroup G π' n = j` であって、
  `Algebra A ↥(fixedRing B H)` を**含まない**。含むのは仮定 `hadjC` / `hresC` の側だけで、
  それは本ファイルが**自分で供給する**。★したがって §3 の系の statement からは
  `Algebra A C` が完全に消え、`letI` で足りる。
* 大域 `instance` にすると、`A := ℤ` のとき mathlib の `Int.instAlgebra` と
  `Algebra ℤ ↥(fixedRing B H)` が二重になる(ダイヤモンド)。★避ける理由がある一方、
  上のとおり**大域にする理由が 1 つも無い**。
* ★`#117(i)`(素の `def` 包みは `instances` 透明度で展開されない)とは衝突しなかった。
  `fixedRing B H : Subring B` は**型としては `↥(Subring …)`** なので、
  `Subring.toAlgebra`(`Algebra ↥S B`)も `IsDomain ↥S` も普通に降りてくる。
  展開が要るのは `fixedRing` を**定義まで開く**ときだけで、本ファイルはそれをしない。

## ★★2 は剰余体の同型を作らなかった

★**`CharP` だけを移した。** 剰余体の同型 `ResidueField C ≅ ResidueField B` は作っていない。
理由は測って分かったことで、**同型どころか `IsLocalHom` すら要らない**からである:

* `CharP (ResidueField B) p` ⇒ `(p : B) ∈ 𝔪_B`。
* `algebraMap C B (p : C) = (p : B)` なので Y17 の `mem_maximalIdeal_of_map_mem`
  (★単元の像が単元、という**易しい向き**)で `(p : C) ∈ 𝔪_C`。
* `ResidueField C` は体(`Nontrivial`)なので `CharP.charP_iff_prime_eq_zero` で戻せる。

★★**したがって `hresA`(`K′/K` の完全分岐)は 2 には要らない。**
段取りの見立て(「完全分岐を落とすと 2 が偽」)は**外れていた** —— 不分岐でも剰余体は
伸びるだけで標数は `p` のままである。★本ファイルの `charP_residueField_fixedRing` は
`hresA` を**受け取っていない**(実機で確かめた)。

## ★3 は Y17 の 2 本がそのまま嵌った

* `hfix` = Y17 `fixedRing_mem_adjoin_uniformizer hresA hAne hϖ`。
* `hadjC` = Y17 `adjoin_fixedRing_uniformizer_eq_top hresA hAne hϖ`。
  ★**そのまま `hadjC` だった**(型を突き合わせて確かめた。`letI` の入れ方が
  Y17 の statement と同じ `Subring.algebraOfMapsTo (fixedRing B H) algebraMap_mem_fixedRing`
  なので、`exact` の既定透明度で一致する)。

## ★逸脱の記録

1. ★**既存の `Found/PGC/*.lean` は 1 行も書き換えていない**(import のみ)。
   Y15b の `HasseArfInduction.lean` から仮定を消して回っていない。仮定が消えた系は
   **本ファイルの中に新しい名前**(`exists_natCast_herbrandPhiGroup_of_cyclic_quotient_base`)
   で作った。
2. ★**`[IsLocalRing A]` / `[IsNoetherian A B]` / `hAne` を足した**。原典には無い。
   これは Y17 の逸脱 2・3 をそのまま引き継いだものである(Y17 の `hfix` / `hadjC` を
   呼ぶのに要る)。★`hAne` は `A` が DVR で `A → B` が単射なら自動であり、
   その形も §3 に置いた(`..._of_injective`)。
3. ★**`hind`(帰納法の再帰呼び出し)は Y15b と同じく仮定のまま**である。
   ★本ノードの守備範囲ではない(「部分群 `H ⊆ G` すべてについて `φ_H` の主張を
   強帰納法で回す」新ノードが要る。Y15b の報告どおり)。
4. ★§1 の 2 本は原典に対応物が無い(配管である)。★一般化であって逸脱ではない。

## 退化の自己検査

* ★**`[Fintype ↥H]` を落とすと `hresC` と 2 が出ない** —— `C = B^H` が DVR
  (したがって局所環)にならないからである(Y14 の実測)。★本ファイルの
  `exists_sub_mem_maximalIdeal_fixedRing` / `charP_residueField_fixedRing` は
  どちらも `[Fintype ↥H]` を要求している。
* ★★**`hresA`(`K′/K` の完全分岐)は 2 には要らない**(上記)。
  ★ただし **1 の `hresC` には要る** —— `hresC` は「`K″/K` も完全分岐」そのものだから、
  底が完全分岐でなければ偽である。★段 2 全体では `hresA` は必須のままである。
* ★**`H ⊴ G` は `fixedRingAlgebra` / `hresC` / 2 では要らない**(固定環は正規性なしに
  部分環である。Y17 の実測どおり)。★要るのは
  `smul_algebraMap_fixedRing_eq_self` 以降 —— そこで初めて `G` の**固定環への作用**
  (Y16 `fixedRingMulSemiringAction`、`[H.Normal]` が要る)を使うからである。
* ★**`hϖ : Irreducible ϖ` を落とすと 3 が偽**(Y17 の逸脱 4。`ϖ` が単元なら
  `𝒪_{K″} ≠ 𝒪[ϖ]`)。
* ★**`[SMulCommClass G A B]` を落とすと 1 が出ない** —— `algebraMap A B a` が
  `H` 不変である理由が消え、`fixedRing B H` が `A` の像を含まなくなる。
* ★`ℕ∞` の切り詰め引き算・除算は 1 箇所も書いていない(`lean-idioms.md` #102)。
-/

namespace ABC3.Found.PGC

open IsLocalRing IsDiscreteValuationRing

/-! ## §1 抽象核 —— 局所環の標数と、作用が可換になる十分条件

★分岐・付値・Galois・固定環の語彙が 1 つも出てこない。 -/

section AbstractCore

/-- ★★★**抽象核** —— 局所環の射に沿って剰余体の標数が**降りる**。

`C → B` が(局所環の間の)代数写像で `B` の剰余体が標数 `p` なら、`C` の剰余体も標数 `p`。

段取りは 3 行:
`(p : ResidueField B) = 0` ⇒ `(p : B) ∈ 𝔪_B` ⇒(Y17 `mem_maximalIdeal_of_map_mem`)
`(p : C) ∈ 𝔪_C` ⇒ `(p : ResidueField C) = 0` ⇒ `CharP`。

★★**`IsLocalHom` も剰余体の同型も要らない。** 使うのは「単元の像は単元」という
**易しい向き**だけである。★したがって `C → B` が**不分岐でも完全分岐でも**成り立つ
(標数は伸びない)。

★`p.Prime` は落とせない —— `CharP.charP_iff_prime_eq_zero` の仮定であり、
`(p : k) = 0` だけからは `CharP k p` は出ない(`p = 4`、`k` が標数 2 の体)。 -/
theorem charP_residueField_of_algebraMap {C B : Type*} [CommRing C] [IsLocalRing C] [CommRing B]
    [IsLocalRing B] [Algebra C B] (p : ℕ) (hp : p.Prime) [CharP (ResidueField B) p] :
    CharP (ResidueField C) p := by
  have hB : ((p : ℕ) : ResidueField B) = 0 := (CharP.charP_iff_prime_eq_zero hp).1 ‹_›
  have hmB : (p : B) ∈ maximalIdeal B := by
    rw [← residue_eq_zero_iff, map_natCast]; exact hB
  have hmC : (p : C) ∈ maximalIdeal C := by
    refine mem_maximalIdeal_of_map_mem (B := B) ?_
    rwa [map_natCast]
  refine (CharP.charP_iff_prime_eq_zero hp).2 ?_
  rw [← map_natCast (residue C) p, residue_eq_zero_iff]
  exact hmC

/-- ★★**抽象核** —— 環作用が底の像を各点固定するなら、その作用は係数と可換である。

`m • (a • s) = m • (algebraMap a * s) = algebraMap a * (m • s) = a • (m • s)`。
★`Algebra.smul_def` と `smul_mul'` の 1 行。**群も分岐も出てこない**(`Monoid M` でよい)。 -/
theorem smulCommClass_of_smul_algebraMap {M A S : Type*} [Monoid M] [CommSemiring A] [Semiring S]
    [Algebra A S] [MulSemiringAction M S]
    (h : ∀ (m : M) (a : A), m • algebraMap A S a = algebraMap A S a) : SMulCommClass M A S :=
  ⟨fun m a s => by rw [Algebra.smul_def, Algebra.smul_def, smul_mul', h]⟩

end AbstractCore

/-! ## §2 具体層 —— `C = B^H` に代入する -/

/-- ★★**Y15b の債務 1 を供給する** —— 底環 `A` の代数構造を固定環 `B^H` に降ろす。

★**Y17 の `Subring.algebraOfMapsTo` をそのまま再利用した**(新しく書いていない)。
必要な `∀ a : A, algebraMap A B a ∈ fixedRing B H` は Y17 の `algebraMap_mem_fixedRing`
(`[SMulCommClass G A B]` から Y14 の `smul_algebraMap_eq_self` 1 行)。

★★**大域インスタンスにしない。** ファイル冒頭「1 を大域インスタンスにしなかった理由」を
参照(消費側の statement に `Algebra A C` が現れないので `letI` で足りる。
`A := ℤ` のときの `Int.instAlgebra` とのダイヤモンドも避けられる)。 -/
@[reducible]
def fixedRingAlgebra (A : Type*) [CommRing A] (B : Type*) [CommRing B] [Algebra A B]
    {G : Type*} [Group G] [MulSemiringAction G B] [SMulCommClass G A B] (H : Subgroup G) :
    Algebra A ↥(fixedRing B H) :=
  Subring.algebraOfMapsTo (fixedRing B H) algebraMap_mem_fixedRing

/-- 固定環に降ろした代数写像は、台 `B` の上では底の代数写像そのものである。★**`rfl`**。 -/
theorem coe_algebraMap_fixedRing {A B : Type*} [CommRing A] [CommRing B] [Algebra A B]
    {G : Type*} [Group G] [MulSemiringAction G B] [SMulCommClass G A B] {H : Subgroup G}
    (a : A) :
    letI := fixedRingAlgebra A B H
    ((algebraMap A ↥(fixedRing B H) a : ↥(fixedRing B H)) : B) = algebraMap A B a := rfl

/-- 塔 `A → B^H → B`。★**`rfl` 一発**(Y17 の実測どおり)。 -/
theorem isScalarTower_fixedRing {A B : Type*} [CommRing A] [CommRing B] [Algebra A B]
    {G : Type*} [Group G] [MulSemiringAction G B] [SMulCommClass G A B] {H : Subgroup G} :
    letI := fixedRingAlgebra A B H
    IsScalarTower A ↥(fixedRing B H) B :=
  letI := fixedRingAlgebra A B H
  IsScalarTower.of_algebraMap_eq fun _ => rfl

/-- ★★**Y15b の債務 1 の `hresC` を供給する** —— `K″/K` も完全分岐である。

底 `K′/K` の完全分岐 `hresA` から、Y17 の抽象核 `exists_sub_mem_maximalIdeal_mid` 1 本で出る。

★**`[Fintype ↥H]` を落とすと使えない** —— `↥(fixedRing B H)` が DVR(したがって局所環)で
なくなるからである(Y14)。★**`H ⊴ G` は要らない。** -/
theorem exists_sub_mem_maximalIdeal_fixedRing {A B : Type*} [CommRing A] [CommRing B] [Algebra A B]
    [IsDomain B] [IsDiscreteValuationRing B] {G : Type*} [Group G] [MulSemiringAction G B]
    [SMulCommClass G A B] {H : Subgroup G} [Fintype H]
    (hresA : ∀ b : B, ∃ a : A, b - algebraMap A B a ∈ maximalIdeal B) :
    letI := fixedRingAlgebra A B H
    ∀ c : ↥(fixedRing B H), ∃ a : A,
      c - algebraMap A ↥(fixedRing B H) a ∈ maximalIdeal ↥(fixedRing B H) := by
  letI := fixedRingAlgebra A B H
  haveI := isScalarTower_fixedRing (A := A) (B := B) (H := H)
  exact fun c => exists_sub_mem_maximalIdeal_mid (B := B) hresA c

/-- ★★★**Y15b の債務 2 を供給する** —— `𝒪_{K″} = B^H` の剰余体も標数 `p` である。

§1 の抽象核に `C := ↥(fixedRing B H)` を代入するだけ(`Algebra ↥(fixedRing B H) B` は
mathlib の `Subring.toAlgebra`、`IsLocalRing ↥(fixedRing B H)` は Y14 の DVR インスタンス)。

★★**`hresA`(完全分岐)は要らない**(ファイル冒頭「2 は剰余体の同型を作らなかった」)。
★**`H ⊴ G` も要らない。**
★**`[Fintype ↥H]` は要る**(`C` が局所環であること自体が Y14 の DVR から来るため)。 -/
theorem charP_residueField_fixedRing {B : Type*} [CommRing B] [IsDomain B]
    [IsDiscreteValuationRing B] {G : Type*} [Group G] [MulSemiringAction G B] {H : Subgroup G}
    [Fintype H] (p : ℕ) (hp : p.Prime) [CharP (ResidueField B) p] :
    CharP (ResidueField ↥(fixedRing B H)) p :=
  charP_residueField_of_algebraMap (B := B) p hp

/-- `G` は固定環に降ろした底の像を各点固定する。★`Subtype.ext` 1 語。

★ここで初めて `[H.Normal]` が要る(Y16 `fixedRingMulSemiringAction` を使うから)。 -/
theorem smul_algebraMap_fixedRing_eq_self {A B : Type*} [CommRing A] [CommRing B] [Algebra A B]
    {G : Type*} [Group G] [MulSemiringAction G B] [SMulCommClass G A B] {H : Subgroup G}
    [H.Normal] (ρ : G) (a : A) :
    letI := fixedRingAlgebra A B H
    ρ • algebraMap A ↥(fixedRing B H) a = algebraMap A ↥(fixedRing B H) a :=
  letI := fixedRingAlgebra A B H
  Subtype.ext (smul_algebraMap_eq_self ρ a)

/-- ★★**Y15b の債務 1 の `SMulCommClass (G ⧸ H) A C` を供給する。**

商の作用は代表元の作用(`hq`、Y15b の §4 が `rfl` で供給する)なので、
`smul_algebraMap_fixedRing_eq_self` を §1 の抽象核に渡すだけ。

★`hq` を仮定で受けるのは Y15b の `faithfulSMul_quotient_fixedRing` と同じ形である
(`[MulSemiringAction (G ⧸ H) C]` は文脈ごとに `letI` で選ぶものだから)。 -/
theorem smulCommClass_quotient_fixedRing {A B : Type*} [CommRing A] [CommRing B] [Algebra A B]
    {G : Type*} [Group G] [MulSemiringAction G B] [SMulCommClass G A B] {H : Subgroup G}
    [H.Normal] [MulSemiringAction (G ⧸ H) ↥(fixedRing B H)]
    (hq : ∀ (σ : G) (c : ↥(fixedRing B H)), (QuotientGroup.mk σ : G ⧸ H) • c = σ • c) :
    letI := fixedRingAlgebra A B H
    SMulCommClass (G ⧸ H) A ↥(fixedRing B H) := by
  letI := fixedRingAlgebra A B H
  refine smulCommClass_of_smul_algebraMap (fun q a => ?_)
  induction q using QuotientGroup.induction_on with
  | _ σ =>
    rw [hq σ]
    exact smul_algebraMap_fixedRing_eq_self σ a

/-! ## §3 系 —— Hasse-Arf 段 2 から配管 3 本が消えた形 -/

/-- ★★★★★**Yoshida 2008 Theorem 6.11 の段 2 —— 配管 3 本が消えた形**。

Y15b `exists_natCast_herbrandPhiGroup_of_cyclic_quotient_fixedRing` から
次の**すべて**が消えた:

* `[Algebra A ↥(fixedRing B H)]`(§2 `fixedRingAlgebra`)
* `[SMulCommClass (G ⧸ H) A ↥(fixedRing B H)]`(§2 `smulCommClass_quotient_fixedRing`)
* `hresC`(§2 `exists_sub_mem_maximalIdeal_fixedRing`)
* `[CharP (ResidueField ↥(fixedRing B H)) p]`(§2 `charP_residueField_fixedRing`)
* `hfix`(Y17 `fixedRing_mem_adjoin_uniformizer`)
* `hadjC`(Y17 `adjoin_fixedRing_uniformizer_eq_top`)
* `[MulSemiringAction (G ⧸ H) ↥(fixedRing B H)]` と `hq`(Y15b §4)

★★**残る仮定は原典 §6.1–§6.2 の底の設定と `hind` だけ**である:

* `hπ'` / `hϖ` : `π′`・`ϖ` はそれぞれ `𝒪_{K′}`・`𝒪_{K″}` の素元。
* `hresA` : ★**`K′/K` が完全分岐**。★落とすと `hresC` が偽になる。
* `hadj` : `𝒪_{K′} = 𝒪[π′]`(原典 Lemma 5.11)。
* `hAne` : ★逸脱 2(`A` が DVR で `A → B` が単射なら自動。下の `..._of_injective`)。
* `h1` : 原文「First assume G = G_1」。
* `hgen` / `hxn` / `hxH` / `hle` : 原文「we can find H with G/H ≅ ℤ/p^{m_i}ℤ,
  and G_nH/H ≠ G_{n+1}H/H」が与えるデータ(Y15b §2 が構成する)。
* ★★`hind` : 原文「by inductive hypothesis」。★**これだけは残る**
  ——「部分群 `H ⊆ G` すべてについて強帰納法を回す」新ノードが要る(Y15b の報告どおり)。 -/
theorem exists_natCast_herbrandPhiGroup_of_cyclic_quotient_base
    {A B : Type*} [CommRing A] [IsLocalRing A] [CommRing B] [Algebra A B] [IsDomain B]
    [IsDiscreteValuationRing B] [IsNoetherian A B]
    {G : Type*} [Group G] [MulSemiringAction G B] [Fintype G] [SMulCommClass G A B]
    [FaithfulSMul G B] {H : Subgroup G} [H.Normal] [Fintype H] [Fintype (G ⧸ H)]
    (p : ℕ) (hp : p.Prime) [CharP (ResidueField B) p]
    {π' : B} (hπ' : Irreducible π') {ϖ : ↥(fixedRing B H)} (hϖ : Irreducible ϖ)
    (hresA : ∀ b : B, ∃ a : A, b - algebraMap A B a ∈ maximalIdeal B)
    (hadj : Algebra.adjoin A ({π'} : Set B) = ⊤)
    (hAne : ∃ a ∈ maximalIdeal A, algebraMap A B a ≠ 0)
    (h1 : lowerRamificationGroup B G 1 = ⊤)
    {σ : G} (hgen : ∀ g : G, ∃ t : ℤ, (σ ^ t)⁻¹ * g ∈ H)
    {n k : ℕ} (hind : herbrandPhi π' H (n : ℝ) = (k : ℝ))
    {x : G} (hxn : x ∈ lowerRamificationGroup B G n) (hxH : x ∉ H)
    (hle : lowerRamificationGroup B G (n + 1) ≤ H) :
    ∃ j : ℕ, herbrandPhiGroup G π' (n : ℝ) = (j : ℝ) := by
  letI := fixedRingAlgebra A B H
  letI := quotientMulSemiringActionOfTrivial (↥(fixedRing B H))
    (smul_fixedRing_eq_self (B := B) (H := H))
  haveI := smulCommClass_quotient_fixedRing (A := A) (B := B) (G := G) (H := H)
    (fun σ c => quotientSMul_mk_fixedRing σ c)
  haveI := charP_residueField_fixedRing (B := B) (H := H) p hp
  exact exists_natCast_herbrandPhiGroup_of_cyclic_quotient_fixedRing (A := A) p hp hπ' hϖ
    hresA hadj (fixedRing_mem_adjoin_uniformizer hresA hAne hϖ)
    (fun σ c => quotientSMul_mk_fixedRing σ c)
    (adjoin_fixedRing_uniformizer_eq_top hresA hAne hϖ)
    (exists_sub_mem_maximalIdeal_fixedRing hresA) h1 hgen hind hxn hxH hle

/-- ★**上の系の使いやすい形** —— `A` が DVR で `A → B` が単射なら `hAne` は自動。

★原典の設定(`K ⊆ K′` は体の拡大、`𝒪_K` は DVR)ではこれが常に成り立つので、
実質的に `hAne`(逸脱 2)は消えている。★Y17 の `..._of_injective` と同じ形である。 -/
theorem exists_natCast_herbrandPhiGroup_of_cyclic_quotient_base_of_injective
    {A B : Type*} [CommRing A] [IsDomain A] [IsDiscreteValuationRing A] [CommRing B] [Algebra A B]
    [IsDomain B] [IsDiscreteValuationRing B] [IsNoetherian A B]
    {G : Type*} [Group G] [MulSemiringAction G B] [Fintype G] [SMulCommClass G A B]
    [FaithfulSMul G B] {H : Subgroup G} [H.Normal] [Fintype H] [Fintype (G ⧸ H)]
    (p : ℕ) (hp : p.Prime) [CharP (ResidueField B) p]
    {π' : B} (hπ' : Irreducible π') {ϖ : ↥(fixedRing B H)} (hϖ : Irreducible ϖ)
    (hresA : ∀ b : B, ∃ a : A, b - algebraMap A B a ∈ maximalIdeal B)
    (hadj : Algebra.adjoin A ({π'} : Set B) = ⊤)
    (hAinj : Function.Injective (algebraMap A B))
    (h1 : lowerRamificationGroup B G 1 = ⊤)
    {σ : G} (hgen : ∀ g : G, ∃ t : ℤ, (σ ^ t)⁻¹ * g ∈ H)
    {n k : ℕ} (hind : herbrandPhi π' H (n : ℝ) = (k : ℝ))
    {x : G} (hxn : x ∈ lowerRamificationGroup B G n) (hxH : x ∉ H)
    (hle : lowerRamificationGroup B G (n + 1) ≤ H) :
    ∃ j : ℕ, herbrandPhiGroup G π' (n : ℝ) = (j : ℝ) :=
  exists_natCast_herbrandPhiGroup_of_cyclic_quotient_base p hp hπ' hϖ hresA hadj
    (exists_mem_maximalIdeal_map_ne_zero hAinj) h1 hgen hind hxn hxH hle

end ABC3.Found.PGC
