import ABC3.Found.PGC.FixedRingTower

/-!
# `H ⊴ G` なら `G` は固定環 `B^H` に環作用する(★原典が名前を付けていない基盤ノード)

★★**これは原典が独立の主張として立てていない基盤ノードである。**
Yoshida は §6.2 冒頭の設定 `G ▷ H with G/H = Gal(K″/K)` に畳んでおり、番号も名前も
付けていない。したがって本ファイルには `.src`(`ABC3.Meta.Source`)を**置かない**
—— 存在しない `sectionId` を書くことになるからである。

**なぜ立てたか**: Y10 `ABC3/Found/PGC/FixedRingRamificationIndex.lean` の Lemma 6.8 形
`card_mul_ramIndex_eq_sum_ramIndex_of_fixedRing` は、固定環 `C` について次の 3 つを
「供給されていない債務」として受け取っている:

```
[MulSemiringAction G C]                                  -- G が固定環 C に作用する
hcomp : ∀ (ρ : G) (c : C), algebraMap C B (ρ • c) = ρ • algebraMap C B c
hHtriv : ∀ ρ ∈ H, ∀ c : C, ρ • c = c
```

Y14 `ABC3/Found/PGC/FixedRingTower.lean` は `C := B^H` を実際に構成して
`hres` / `hAC` / `hinj` / `hinvC` / `hfix` / `hadj` / `[IsDiscreteValuationRing C]` を
全部供給したが、上の 3 つには手を付けなかった。理由は Y14 の逸脱記録 3 にある通り
**`G` が `B^H` に作用するのは `H ⊴ G` のときだけ**で、Y14 は `H` の正規性を
仮定していなかったからである。★**本ファイルがその 1 ノードである。**

## 何を示したか

**§1(抽象核・純群作用)** `smul_mem_fixedPoints_of_normal` ——
`H ⊴ G`、`a` が `H` で固定されるなら `ρ • a` も `H` で固定される。
★**環も付値も分岐も Galois も 1 語も出てこない。** `MulAction G α` だけで述べられる。
段取りは `τ • ρ • a = ρ • ((ρ⁻¹τρ) • a)` と書き替えて `Subgroup.Normal.conj_mem'` を
当てるだけ。★**`H ⊴ G` を使うのはこの 1 箇所だけである**(実機で確認した:
本ファイルで正規性を `.conj_mem'` の形で使っている場所は §1 の 1 行のみ)。

**★mathlib の在庫**: `Mathlib/RingTheory/Invariant/Basic.lean:98` に

```
instance (H : Subgroup G) [H.Normal] : MulSemiringAction (G ⧸ H) (FixedPoints.subring B H)
```

という**無名インスタンスが在る**。ただしこれは**商群 `G ⧸ H` の作用**であって、
Y10 が要求する `MulSemiringAction G C`(= `G` そのものの作用、`H` が自明に効く)ではない。
また `Mathlib.RingTheory.Invariant.Basic` は本木のどのファイルからも import されていない
(`Algebra.IsInvariant` が `Unknown constant` になる)。
★そこで本ファイルは `G` の作用を直接構成した —— 商の作用を
`MulSemiringAction.compHom (QuotientGroup.mk' H)` で引き戻す道もあるが、
mathlib import を 1 本増やすだけで得るものが無い(`hcomp` はどちらの道でも `rfl` である)。

**§2(具体層)** `fixedRingMulSemiringAction` : `MulSemiringAction G ↥(fixedRing B H)`
(`fixedRing` は Y14 の `FixedPoints.subring B ↥H`)。`ρ • c := ⟨ρ • (c : B), §1⟩`。
`one_smul` / `mul_smul` / `smul_zero` / `smul_add` / `smul_one` / `smul_mul` は
`Subtype.ext` で `B` から降りる。

**§2b(Y10 の債務の返済)** ★**`hcomp` は `rfl` 1 語、`hHtriv` は `Subtype.ext` 1 語で出た。**
`algebraMap ↥(fixedRing B H) B` が包含写像だからである。段取りの見立て
(「1–3 行で出るはず」)は当たった。

**§3(系)** `card_mul_ramIndex_eq_sum_ramIndex_fixedRing` —— Y10 の Lemma 6.8 形に
§1・§2 と Y14 の供給物を全部代入した形。★仮定は
「原典 §6.1 の設定 + `H ⊴ G` + `[Fintype ↥H]` + `hfix`」だけになった
(`[MulSemiringAction G C]` / `hcomp` / `hHtriv` / `hinj` / `hfixC` / `hres` / `hAC` /
`[IsDiscreteValuationRing C]` は**すべて消えた**)。

★`hfix`(`𝒪_{K′}^H ⊆ 𝒪[π″]`、原典 Lemma 5.11 を `K″/K` に当てた形)だけは残る。
これは本ノードの守備範囲ではない(`K″/K` 側の Lemma 5.11 であって、固定環への作用とは
無関係の債務である)。

## ★逸脱の記録

1. ★**既存の `Found/PGC/*.lean` は 1 行も書き換えていない。** Y10・Y11・Y12・Y14 から
   仮定を消して回るかどうかは上流の判断に委ねる。本ファイルは §3 の系を**新しい名前で**
   足しただけである。
2. ★**`H ⊴ G` は原典 §6.2 の設定そのもの**である(`G/H = Gal(K″/K)` と書くには
   `H ⊴ G` が要る)。弱めても強めてもいない。
3. ★§1 は `MulAction G α` で述べた(原典は環の話をしているが、正規性が効くのは
   作用の部分だけである)。★一般化であって逸脱ではない。

## 退化の自己検査

* ★★**`H ⊴ G` を落とすと §1・§2 は偽**である。`c ∈ B^H`、`ρ ∈ G` に対し
  `τ • (ρ • c) = ρ • ((ρ⁻¹τρ) • c)` であって、`ρ⁻¹τρ ∈ H` でなければ右辺は `ρ • c` に
  戻らない。★**これが本ノードの中身である。**
* ★**`[Fintype ↥H]` は §1・§2 では要らない**(作用の定義に有限性は 1 つも要らない)。
  実機で確認した —— §2 の `instance` は `[Fintype ↥H]` を持っていない。
  ★§3 では要る。理由は作用ではなく **Y14 の `IsDiscreteValuationRing ↥(fixedRing B H)`**
  (ノルム元 `∏_{τ∈H} τ•π` を作るのに有限性が要る)と、Y10 の `Nat.card H` の正値性である。
* ★**`[FaithfulSMul G B]` は §1・§2 では要らない**(`H` が自明に作用すれば
  `fixedRing B H = ⊤` で、`G` は `B` に作用しているのだからその上でも作用する)。
  ★§3 では要る(Y10 が発見した退化条件)。
* ★`fixedRing` の担い手 `B` は明示引数である(Y14 が `lean-idioms.md` #116(a) に
  記録した失敗形。`B` を暗黙にすると `MulSemiringAction G ?m` でメタ変数が残る)。
  本ファイルの `instance` も `(B : Type*)` を明示引数で受けている。
-/

namespace ABC3.Found.PGC

open IsLocalRing IsDiscreteValuationRing

/-! ## §1 抽象核 —— 正規部分群の固定点は全体の作用で保たれる -/

/-- ★★★**抽象核** —— `H ⊴ G` のとき、`H` の固定点全体は `G` の作用で閉じている。

`a` が `H` のすべての元で固定されるなら、`ρ • a`(`ρ : G` は任意)もそうである。

★**環・付値・分岐・Galois の語彙が 1 つも出てこない。** `MulAction G α` だけでよい。

★★**本ファイルで `H ⊴ G` を使うのはこの 1 行(`hH.conj_mem'`)だけである。**
正規性を落とすと主張は偽になる。 -/
theorem smul_mem_fixedPoints_of_normal {α G : Type*} [Group G] [MulAction G α]
    {H : Subgroup G} (hH : H.Normal) (ρ : G) {a : α} (ha : ∀ τ ∈ H, τ • a = a)
    {τ : G} (hτ : τ ∈ H) : τ • (ρ • a) = ρ • a := by
  have hc : ρ⁻¹ * τ * ρ ∈ H := hH.conj_mem' τ hτ ρ
  calc τ • (ρ • a) = (τ * ρ) • a := (mul_smul τ ρ a).symm
    _ = (ρ * (ρ⁻¹ * τ * ρ)) • a := by group
    _ = ρ • ((ρ⁻¹ * τ * ρ) • a) := mul_smul _ _ a
    _ = ρ • a := by rw [ha _ hc]

/-! ## §2 具体層 —— `G` は `B^H` に環作用する -/

/-- ★★★**`H ⊴ G` なら `G` は固定環 `B^H` に環作用する。**

`ρ • c` は `B` の中の `ρ • (c : B)` であり、それが再び `B^H` に入ることが §1 である。
環構造との両立(`smul_add` / `smul_mul` / `smul_one` / `smul_zero`)は
`Subtype.ext` で `B` から降りる。

★**`[Fintype ↥H]` は要らない** —— 作用の定義に有限性は 1 つも使わない。
★`B` は明示引数である(`lean-idioms.md` #116(a))。

★mathlib には `MulSemiringAction (G ⧸ H) (FixedPoints.subring B ↥H)`
(`Mathlib/RingTheory/Invariant/Basic.lean:98` の無名インスタンス)が在るが、
それは**商群の作用**であり、また当該ファイルは本木から import されていない。 -/
instance fixedRingMulSemiringAction (B : Type*) [CommRing B] {G : Type*} [Group G]
    [MulSemiringAction G B] (H : Subgroup G) [hH : H.Normal] :
    MulSemiringAction G ↥(fixedRing B H) where
  smul ρ c := ⟨ρ • (c : B), mem_fixedRing.2 fun _ hτ =>
    smul_mem_fixedPoints_of_normal hH ρ (mem_fixedRing.1 c.2) hτ⟩
  one_smul c := Subtype.ext (one_smul G (c : B))
  mul_smul ρ σ c := Subtype.ext (mul_smul ρ σ (c : B))
  smul_zero ρ := Subtype.ext (smul_zero ρ)
  smul_add ρ c d := Subtype.ext (smul_add ρ (c : B) (d : B))
  smul_one ρ := Subtype.ext (smul_one ρ)
  smul_mul ρ c d := Subtype.ext (MulSemiringAction.smul_mul ρ (c : B) (d : B))

/-- 作用は台 `B` の上では元の作用そのものである。★定義から `rfl`。 -/
@[simp]
theorem coe_smul_fixedRing {B : Type*} [CommRing B] {G : Type*} [Group G]
    [MulSemiringAction G B] {H : Subgroup G} [H.Normal] (ρ : G) (c : ↥(fixedRing B H)) :
    ((ρ • c : ↥(fixedRing B H)) : B) = ρ • (c : B) := rfl

/-! ### §2b Y10 の債務 `hcomp` / `hHtriv` の返済 -/

/-- ★★**Y10 の `hcomp` を供給する** —— 包含写像は `G` 同変である。

★`algebraMap ↥(fixedRing B H) B` は包含写像なので **`rfl` 1 語で出た**。 -/
theorem algebraMap_smul_fixedRing {B : Type*} [CommRing B] {G : Type*} [Group G]
    [MulSemiringAction G B] {H : Subgroup G} [H.Normal] (ρ : G) (c : ↥(fixedRing B H)) :
    algebraMap (↥(fixedRing B H)) B (ρ • c) = ρ • algebraMap (↥(fixedRing B H)) B c := rfl

/-- ★★**Y10 の `hHtriv` を供給する** —— `H` は固定環に自明に作用する。

★固定環の定義そのもの。**`Subtype.ext` 1 語で出た。** -/
theorem smul_fixedRing_eq_self {B : Type*} [CommRing B] {G : Type*} [Group G]
    [MulSemiringAction G B] {H : Subgroup G} [H.Normal] :
    ∀ ρ ∈ H, ∀ c : ↥(fixedRing B H), ρ • c = c :=
  fun ρ hρ c => Subtype.ext (mem_fixedRing.1 c.2 ρ hρ)

/-- ★`G` の作用は `H` を経由しない —— 商群 `G/H = Gal(K″/K)` が固定環に作用する、の実質。
`ρ` と `ρ * τ`(`τ ∈ H`)は固定環の上で同じ作用をする。 -/
theorem smul_mul_mem_fixedRing {B : Type*} [CommRing B] {G : Type*} [Group G]
    [MulSemiringAction G B] {H : Subgroup G} [H.Normal] {τ : G} (hτ : τ ∈ H) (ρ : G)
    (c : ↥(fixedRing B H)) : (ρ * τ) • c = ρ • c := by
  rw [mul_smul, smul_fixedRing_eq_self τ hτ c]

/-! ## §3 系 —— Y10 の Lemma 6.8 形を底の設定だけに落とす -/

/-- ★★★★**Y10 `card_mul_ramIndex_eq_sum_ramIndex_of_fixedRing` を `C = B^H` に代入した形**。

```
|H| · i_{K″}(σ ϖ) = Σ_{τ ∈ H} i_{K′}(σ τ)
```

★**仮定は原典 §6.1–§6.2 の底の設定だけになった**:
* `hπ'` : `π′` は `𝒪_{K′}` の素元。
* `hresA` : ★**`K′/K` が完全分岐**(剰余体が伸びない)。★落とすと偽。
* `hadjA` : `𝒪_{K′} = 𝒪[π′]`(原典 Lemma 5.11)。
* `hfix` : `𝒪_{K′}^H ⊆ 𝒪[π″]`(原典 Lemma 5.11 を `K″/K` に当てた形)。
  ★これだけは残る —— 固定環への作用とは無関係の債務である。
* `[SMulCommClass G A B]` : `G` は `K` 上の作用である。
* `[FaithfulSMul G B]` : ★落とすと偽(Y10 が発見した退化条件)。
* `[H.Normal]` : ★★**本ノードが供給した条件。落とすと `G` は `B^H` に作用しない。**
* `[Fintype ↥H]` : ★落とすと `Nat.card H = 0` で空虚に真。

★Y10 が持ち回っていた `[MulSemiringAction G C]` / `hcomp` / `hHtriv` /
`hinj` / `hfixC` / `hres` / `hAC` / `[IsDiscreteValuationRing C]` は**すべて消えた**。 -/
theorem card_mul_ramIndex_eq_sum_ramIndex_fixedRing
    {A B : Type*} [CommRing A] [CommRing B] [Algebra A B] [IsDomain B]
    [IsDiscreteValuationRing B] {G : Type*} [Group G] [MulSemiringAction G B]
    [SMulCommClass G A B] [FaithfulSMul G B] {H : Subgroup G} [H.Normal] [Fintype H]
    {π' : B} (hπ' : Irreducible π')
    (hresA : ∀ b : B, ∃ a : A, b - algebraMap A B a ∈ maximalIdeal B)
    (hadjA : Algebra.adjoin A ({π'} : Set B) = ⊤)
    (ϖ : ↥(fixedRing B H))
    (hfix : ∀ y : B, (∀ ρ ∈ H, ρ • y = y) →
      y ∈ Algebra.adjoin A ({algebraMap (↥(fixedRing B H)) B ϖ} : Set B))
    (σ : G) :
    (Nat.card H : ℕ∞) * ramIndex ϖ σ = ∑ τ : H, ramIndex π' (σ * (τ : G)) :=
  card_mul_ramIndex_eq_sum_ramIndex_of_fixedRing (C := ↥(fixedRing B H))
    algebraMap_smul_fixedRing smul_fixedRing_eq_self hπ' fixedRing_injective
    exists_algebraMap_fixedRing (exists_sub_mem_fixedRing hresA)
    (fun a => exists_algebraMap_fixedRing_eq a) rfl hadjA hfix σ

end ABC3.Found.PGC
