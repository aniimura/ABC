import ABC3.Found.PGC.HerbrandComposition

/-!
# 固定環の塔 `𝒪 ⊆ 𝒪_{K″} = 𝒪_{K′}^H ⊆ 𝒪_{K′}`(★原典が名前を付けていない基盤ノード)

★★**これは原典が独立の主張として立てていない基盤ノードである。**
Yoshida は §6.1–§6.2 の設定(`K′/K` は完全分岐、`𝒪_{K′} = 𝒪[π′]`、`K″ = (K′)^H`)の中に
畳んでおり、番号も名前も付けていない。したがって本ファイルには `.src`(`ABC3.Meta.Source`)を
**置かない** —— 存在しない `sectionId` を書くことになるからである。

**なぜ立てたか**: Y10 `ABC3/Found/PGC/FixedRingRamificationIndex.lean` の主結論
`addVal_map_eq_card_mul` は 2 つの仮定を「供給されていない債務」として持ち回っており、
それを Y11 `HerbrandFunction.lean`(Prop 6.9)と Y12 `HerbrandComposition.lean`(Lemma 6.10)が
そのまま継承している:

```
hres : ∀ b : B, ∃ c : C, b - algebraMap C B c ∈ maximalIdeal B   -- K′/K″ が完全分岐 (f = 1)
hAC  : ∀ a : A, ∃ c : C, algebraMap C B c = algebraMap A B a     -- 𝒪 ⊆ 𝒪_{K″}
```

さらに Y10 は `C`(= `𝒪_{K″}`)を部分環として**構成せず**、
`[Algebra C B] + hinj + hinvC + hfix + [IsDiscreteValuationRing C]` の形で受け取っている。
★**本ファイルは `C := B^H` を実際に構成し、Y10 が要求する引数を全部供給する。**

## 何を示したか

**§1(抽象核・純可換環)** `exists_sub_mem_of_exists_image` ——
`𝒪` からの剰余全射性(`hresA`)と `𝒪` の像が `C` に入ること(`hAC`)から、
`C` からの剰余全射性(`hres`)が出る。★分岐・付値・Galois の語彙が 1 つも出てこない。
イデアルは任意の `I : Ideal B` でよい。
★塔 `[IsScalarTower A C B]` があれば `hAC` は自動である
(`exists_image_of_isScalarTower`)。★★**本ノードの実質は `hres` ではなく
「`C` を作って `hAC` を出すこと」であった** —— `hres` は `hAC` から 3 行で出る。

**§2(具体層)** `fixedRing B H := FixedPoints.subring B ↥H`(mathlib の固定点部分環、**道 A**)。
`Subring B` なので `Algebra ↥(fixedRing B H) B` と `IsDomain ↥(fixedRing B H)` は
mathlib のインスタンスがそのまま付く。Y10 が要求する引数を
`fixedRing_injective` / `smul_algebraMap_fixedRing` / `exists_algebraMap_fixedRing` /
`exists_algebraMap_fixedRing_eq` / `exists_sub_mem_fixedRing` / `adjoin_fixedRing_eq_top`
として供給した。
★`hAC` は「`G` が `𝒪` の像を各点固定する」ことから出る。それは `[SMulCommClass G A B]` の
1 行(`smul_algebraMap_eq_self`)である。★これも抽象核であって、分岐を含まない。

**§3(抽象核・付値のみ)** `isDiscreteValuationRing_of_saturated` ——
★**離散付値環 `B` の部分環 `C` は、(i) 単元が降りる、(ii) 割り算で閉じている
(`c = d·x`、`d ≠ 0`、`c, d ∈ C` なら `x ∈ C`)、(iii) 非零非単元を持つ、の 3 条件で
それ自身 DVR である。** 群も Galois も出てこない。
段取りは「`v_B(C∖0) ⊆ ℕ` は差で閉じた部分モノイドだから最小正元 `n₀` の倍数全体」で、
`Nat.find` で `n₀` を取り、`k` に関する強帰納法で `c = u·ϖ₀^n` を出し、
`IsDiscreteValuationRing.ofHasUnitMulPowIrreducibleFactorization` に渡す。
具体層 `fixedRing_sat` / `fixedRing_exists_nonzero_nonunit` は
「`H` 不変な `c, d` の商は `H` 不変」「ノルム `∏_{τ∈H} τ•π` の付値は `|H| > 0`」だけ。
⇒ **インスタンス `IsDiscreteValuationRing ↥(fixedRing B H)`**。★§3 まで入った。

**§4(主結論)** `addVal_map_eq_card_mul_fixedRing` ——
Y10 の主結論を `C = B^H` に代入した形。★仮定は**底の設定だけ**になった:
`hπ`(`π` は `𝒪_{K′}` の素元)/ `hresA`(`K′/K` 完全分岐)/ `hadjA`(`𝒪_{K′} = 𝒪[π′]`)/
`[SMulCommClass G A B]` / `[FaithfulSMul G B]` / `[Fintype ↥H]`。
★`hres` も `hAC` も `hinj` も `hfix` も `hadj` も `[IsDiscreteValuationRing C]` も**消えた**。

## ★逸脱の記録

1. ★**`𝒪` を `A` として抽象的に受け、`[SMulCommClass G A B]` で「`G` が `K` 上の作用である」を
   表した。** 原典は `G = Gal(K′/K)` と書くので `K` の元を固定するのは定義から明らかだが、
   この木では `MulSemiringAction G B` + `SMulCommClass G A B` がその形式化である
   (Y9・Y10 の `card_mul_ramIndex_eq_sum_ramIndex_of_fixedRing` が同じ形を使っている)。
   ★弱めても強めてもいない。
2. ★**`𝒪` が体でないこと・`A → B` が単射であることは使っていない。** `A` は単に環でよい。
   実際に使うのは「`G` が `A` の像を固定する」ことだけである。
3. ★**Y10 の Lemma 6.8 形 `card_mul_ramIndex_eq_sum_ramIndex_of_fixedRing` は
   代入していない。** それは `[MulSemiringAction G C]` を要求するが、
   `C = B^H` に `G` が作用するのは `H ⊴ G` のときだけである。
   ★本ファイルは `H` の正規性を仮定していないので、そこには手を出さない
   (原典 §6.2 では `H ⊴ G` だが、必要になった時点で別ノードにする)。
4. ★**Y10・Y11・Y12 のファイルは書き換えていない。** 差し替えるかどうかは上流の判断に委ねる。

## 退化の自己検査

* ★★**`K′/K` の完全分岐(`hresA`)を落とすと §4 は偽**である。不分岐拡大では
  `e = 1`、`f = |H|` になり `addVal_map_eq_card_mul_fixedRing` の右辺が合わない。
  ★§1 の抽象核はこれを `hresA → hres` の形でしか触らないので、そこでは退化しない
  (仮定が空なら結論も空)。
* ★**`[Fintype ↥H]` を落とすと `Nat.card H = 0` になり §4 が空虚に真**になる
  (Y9 の `pos_natCard_subgroup` が塞いでいる穴と同じ)。
  ★§3 でも `[Fintype ↥H]` は落とせない —— ノルム元 `∏_{τ∈H} τ•π` が作れなくなり、
  「非零非単元が存在する」(`C` が体でない)が言えない。実際 `H` が無限だと
  `C` は体になりうる。
* ★**`[FaithfulSMul G B]` を落とすと §4 は偽**(Y10 が発見した退化条件)。
  `H` が `B` に自明に作用すれば `fixedRing B H = ⊤`、`C = B`、`e = 1 ≠ |H|`。
  ★§1・§2・§3 は `FaithfulSMul` を使わない(`C = B` でも真である)。
* ★**`algebraMap ↥(fixedRing B H) B` の単射性**は `Subtype.coe_injective` なので
  構成から自動である。★Y10 が `hinj` として仮定していたものが、ここでは定理になった。
* ★§3 で `hd : d ≠ 0` を落とすと飽和性が偽(`0 = 0 · x` は任意の `x` で成り立つ)。

★`FixedPoints.subring B ↥H` の `B` は `H : Subgroup G` から推論できないので、
`fixedRing` は `B` を**明示引数**にしてある(`fixedRing B H`)。
-/

namespace ABC3.Found.PGC

open IsLocalRing IsDiscreteValuationRing

/-! ## §1 抽象核 —— 塔に沿った剰余の移送 -/

/-- ★★★**抽象核** —— 剰余の全射性は塔の途中の環に移る。

`A` からの剰余全射性(`hresA`)と「`A` の像は `C` の像に含まれる」(`hAC`)から、
`C` からの剰余全射性が出る。

★**分岐・付値・Galois の語彙が 1 つも出てこない。** イデアル `I` は任意でよく、
`A`・`C`・`B` は単なる可換環である。

★これが Y10 `addVal_map_eq_card_mul` の仮定 `hres` の供給源である。 -/
theorem exists_sub_mem_of_exists_image {A C B : Type*} [CommRing A] [CommRing C] [CommRing B]
    [Algebra A B] [Algebra C B] {I : Ideal B}
    (hAC : ∀ a : A, ∃ c : C, algebraMap C B c = algebraMap A B a)
    (hresA : ∀ b : B, ∃ a : A, b - algebraMap A B a ∈ I) (b : B) :
    ∃ c : C, b - algebraMap C B c ∈ I := by
  obtain ⟨a, ha⟩ := hresA b
  obtain ⟨c, hc⟩ := hAC a
  exact ⟨c, by rw [hc]; exact ha⟩

/-- ★**抽象核** —— 塔 `A → C → B` が立っていれば `hAC` は自動である。

★本ファイルの具体層は `IsScalarTower` を使わず、`hAC` を固定環の性質から直接出す
(`exists_algebraMap_fixedRing_eq`)。それは `Algebra A ↥(fixedRing B H)` の
インスタンスを立てずに済ませるためである。★この補題は「塔があるなら不要」という
境界を記録するために置いてある。 -/
theorem exists_image_of_isScalarTower {A C B : Type*} [CommRing A] [CommRing C] [CommRing B]
    [Algebra A C] [Algebra A B] [Algebra C B] [IsScalarTower A C B] (a : A) :
    ∃ c : C, algebraMap C B c = algebraMap A B a :=
  ⟨algebraMap A C a, (IsScalarTower.algebraMap_apply A C B a).symm⟩

/-- ★**抽象核** —— `G` が `A` 上線形に作用するなら `G` は `A` の像を各点固定する。

`algebraMap A B a = a • 1` と `[SMulCommClass G A B]` から 1 行。
★これが `𝒪 ⊆ 𝒪_{K″} = B^H`(すなわち `hAC`)の実質である。 -/
theorem smul_algebraMap_eq_self {A B : Type*} [CommRing A] [CommRing B] [Algebra A B]
    {G : Type*} [Group G] [MulSemiringAction G B] [SMulCommClass G A B] (ρ : G) (a : A) :
    ρ • algebraMap A B a = algebraMap A B a := by
  rw [Algebra.algebraMap_eq_smul_one, smul_comm, smul_one]

/-! ## §2 具体層 —— 固定環 `B^H` を作る(道 A) -/

/-- **固定環 `𝒪_{K″} = 𝒪_{K′}^H`**。mathlib の `FixedPoints.subring` を `↥H` に対して使う
(★**道 A**)。`Subring B` なので `Algebra ↥(fixedRing B H) B` と
`IsDomain ↥(fixedRing B H)` は mathlib のインスタンスがそのまま付く。

★`B` は `H : Subgroup G` から推論できないので**明示引数**である。 -/
def fixedRing (B : Type*) [CommRing B] {G : Type*} [Group G] [MulSemiringAction G B]
    (H : Subgroup G) : Subring B :=
  FixedPoints.subring B ↥H

/-- 固定環の所属は「`H` の元がすべて固定する」ことに他ならない。 -/
theorem mem_fixedRing {B : Type*} [CommRing B] {G : Type*} [Group G] [MulSemiringAction G B]
    {H : Subgroup G} {b : B} : b ∈ fixedRing B H ↔ ∀ ρ ∈ H, ρ • b = b :=
  ⟨fun h ρ hρ => h ⟨ρ, hρ⟩, fun h ρ => h ρ ρ.2⟩

/-- `algebraMap ↥(fixedRing B H) B` は包含写像である。 -/
theorem algebraMap_fixedRing_apply {B : Type*} [CommRing B] {G : Type*} [Group G]
    [MulSemiringAction G B] {H : Subgroup G} (c : ↥(fixedRing B H)) :
    algebraMap (↥(fixedRing B H)) B c = (c : B) := rfl

/-- ★**Y10 の `hinj` を供給する** —— 包含写像は単射。 -/
theorem fixedRing_injective {B : Type*} [CommRing B] {G : Type*} [Group G]
    [MulSemiringAction G B] {H : Subgroup G} :
    Function.Injective (algebraMap (↥(fixedRing B H)) B) :=
  Subtype.coe_injective

/-- ★**Y10 の `hinvC` を供給する** —— `H` は `C` の像を各点固定する。 -/
theorem smul_algebraMap_fixedRing {B : Type*} [CommRing B] {G : Type*} [Group G]
    [MulSemiringAction G B] {H : Subgroup G} :
    ∀ ρ ∈ H, ∀ c : ↥(fixedRing B H),
      ρ • algebraMap (↥(fixedRing B H)) B c = algebraMap (↥(fixedRing B H)) B c :=
  fun ρ hρ c => mem_fixedRing.1 c.2 ρ hρ

/-- ★**Y10 の `hfix` を供給する** —— `H` 不変な `b` は `C` から来る。
★構成から自明である(Y10 では仮定だった)。 -/
theorem exists_algebraMap_fixedRing {B : Type*} [CommRing B] {G : Type*} [Group G]
    [MulSemiringAction G B] {H : Subgroup G} :
    ∀ b : B, (∀ ρ ∈ H, ρ • b = b) →
      ∃ c : ↥(fixedRing B H), algebraMap (↥(fixedRing B H)) B c = b :=
  fun b hb => ⟨⟨b, mem_fixedRing.2 hb⟩, rfl⟩

/-- ★★**Y10 の `hAC`(`𝒪 ⊆ 𝒪_{K″}`)を供給する。**

`G` が `A` 上線形に作用する(`[SMulCommClass G A B]`)なら `A` の像は `H` 不変なので
固定環に入る。★これが本ノードの実質である。 -/
theorem exists_algebraMap_fixedRing_eq {A B : Type*} [CommRing A] [CommRing B] [Algebra A B]
    {G : Type*} [Group G] [MulSemiringAction G B] [SMulCommClass G A B] {H : Subgroup G}
    (a : A) :
    ∃ c : ↥(fixedRing B H), algebraMap (↥(fixedRing B H)) B c = algebraMap A B a :=
  exists_algebraMap_fixedRing _ fun ρ _ => smul_algebraMap_eq_self ρ a

/-- ★★**Y10 の `hres`(`K′/K″` が完全分岐)を供給する。**

★底 `K′/K` の完全分岐(`hresA`)から、§1 の抽象核 1 本で出る。
★★**「塔が立てば `hres` はほぼ自明」という見立ては当たった** —— 実際に要ったのは
`hAC`(= `exists_algebraMap_fixedRing_eq`)だけである。 -/
theorem exists_sub_mem_fixedRing {A B : Type*} [CommRing A] [CommRing B] [Algebra A B]
    {G : Type*} [Group G] [MulSemiringAction G B] [SMulCommClass G A B] {H : Subgroup G}
    {I : Ideal B} (hresA : ∀ b : B, ∃ a : A, b - algebraMap A B a ∈ I) (b : B) :
    ∃ c : ↥(fixedRing B H), b - algebraMap (↥(fixedRing B H)) B c ∈ I :=
  exists_sub_mem_of_exists_image exists_algebraMap_fixedRing_eq hresA b

/-- ★**Y10 の `hadj`(`𝒪_{K′} = 𝒪_{K″}[π′]`)を供給する。**
Y10 の橋 `adjoin_eq_top_of_image` に `hAC` を流し込むだけ。 -/
theorem adjoin_fixedRing_eq_top {A B : Type*} [CommRing A] [CommRing B] [Algebra A B]
    {G : Type*} [Group G] [MulSemiringAction G B] [SMulCommClass G A B] {H : Subgroup G}
    {π : B} (hadjA : Algebra.adjoin A ({π} : Set B) = ⊤) :
    Algebra.adjoin (↥(fixedRing B H)) ({π} : Set B) = ⊤ :=
  adjoin_eq_top_of_image exists_algebraMap_fixedRing_eq hadjA

/-! ## §3 抽象核 —— 割り算で閉じた部分環は DVR -/

/-- ★★★**抽象核(§3)** —— 離散付値環 `B` の部分環 `C` が

1. `hunit` : 単元が降りる(`IsUnit (algebraMap C B c) → IsUnit c`)
2. `hsat`  : 割り算で閉じている(`c = d·x`、`d ≠ 0`、`c, d ∈ C` なら `x ∈ C`)
3. `hnu`   : 非零非単元を持つ(`C` が体でない)

を満たせば、`C` 自身が離散付値環である。

★**群作用・Galois・分岐は 1 つも出てこない。** 段取り:

* `v(C∖0) ⊆ ℕ` は加法で閉じ、さらに `hsat` により**差でも閉じる**。
  よって最小の正の値 `n₀`(`Nat.find`)を取れば `v(C∖0) = n₀ℕ`。
* `k` に関する強帰納法で `∀ c ≠ 0, ∃ n, Associated (ϖ₀^n) c`(`v(ϖ₀) = n₀`)。
  各段で `addVal_le_iff_dvd` が `B` での割り切りを与え、`hsat` が商を `C` に戻す。
* `ϖ₀` は既約(`v(a) + v(b) = n₀` で両方が非単元なら `2n₀ ≤ n₀`)。
* `IsDiscreteValuationRing.ofHasUnitMulPowIrreducibleFactorization`。

★退化検査: `hnu` を落とすと偽(`C` が体なら DVR でない)。
★`hd : d ≠ 0` を `hsat` から落とすと `hsat` 自体が偽になる(`0 = 0·x` は任意の `x` で真)。
★`ℕ∞` の切り詰め引き算・除算は 1 箇所も書いていない(`lean-idioms.md` #102)。 -/
theorem isDiscreteValuationRing_of_saturated {C B : Type*} [CommRing C] [IsDomain C]
    [CommRing B] [IsDomain B] [IsDiscreteValuationRing B] [Algebra C B]
    (hinj : Function.Injective (algebraMap C B))
    (hunit : ∀ c : C, IsUnit (algebraMap C B c) → IsUnit c)
    (hsat : ∀ (c d : C) (x : B), d ≠ 0 → algebraMap C B c = algebraMap C B d * x →
      ∃ q : C, algebraMap C B q = x)
    (hnu : ∃ c : C, c ≠ 0 ∧ ¬ IsUnit c) :
    IsDiscreteValuationRing C := by
  classical
  set f : C → ℕ∞ := fun c => addVal B (algebraMap C B c) with hf
  have hne0 : ∀ c : C, c ≠ 0 → algebraMap C B c ≠ 0 := fun c hc h =>
    hc (hinj (by rw [h, map_zero]))
  have hfin : ∀ c : C, c ≠ 0 → ∃ k : ℕ, f c = (k : ℕ∞) := by
    intro c hc
    obtain ⟨k, hk⟩ := ENat.ne_top_iff_exists.mp (fun h => hne0 c hc (addVal_eq_top_iff.mp h))
    exact ⟨k, hk.symm⟩
  have hzero : ∀ c : C, f c = 0 ↔ IsUnit c :=
    fun c => ⟨fun h => hunit c (addVal_eq_zero_iff.mp h),
      fun h => addVal_eq_zero_iff.mpr (h.map (algebraMap C B))⟩
  have hT : ∃ k : ℕ, 0 < k ∧ ∃ c : C, f c = (k : ℕ∞) := by
    obtain ⟨c, hc0, hcu⟩ := hnu
    obtain ⟨k, hk⟩ := hfin c hc0
    refine ⟨k, ?_, c, hk⟩
    rcases Nat.eq_zero_or_pos k with rfl | h
    · exact absurd ((hzero c).1 (by simpa using hk)) hcu
    · exact h
  set n₀ := Nat.find hT with hn₀
  obtain ⟨hn₀pos, ϖ₀, hϖ₀⟩ := Nat.find_spec hT
  have hϖ₀0 : ϖ₀ ≠ 0 := by
    intro h
    rw [h] at hϖ₀
    simp only [hf, map_zero, addVal_zero] at hϖ₀
    exact (ENat.coe_ne_top n₀) hϖ₀.symm
  have hmin : ∀ (k : ℕ) (c : C), f c = (k : ℕ∞) → ¬ IsUnit c → n₀ ≤ k := by
    intro k c hk hcu
    refine Nat.find_le ⟨?_, c, hk⟩
    rcases Nat.eq_zero_or_pos k with rfl | h
    · exact absurd ((hzero c).1 (by simpa using hk)) hcu
    · exact h
  have hv₀ : addVal B (algebraMap C B ϖ₀) = (n₀ : ℕ∞) := hϖ₀
  have key : ∀ k : ℕ, ∀ c : C, f c = (k : ℕ∞) → ∃ n : ℕ, Associated (ϖ₀ ^ n) c := by
    intro k
    induction k using Nat.strong_induction_on with
    | _ k ih =>
      intro c hk
      by_cases hcu : IsUnit c
      · exact ⟨0, by rw [pow_zero]; exact (associated_one_iff_isUnit.mpr hcu).symm⟩
      have hkc : addVal B (algebraMap C B c) = (k : ℕ∞) := hk
      have hc0 : c ≠ 0 := by
        intro h
        rw [h] at hkc
        simp only [map_zero, addVal_zero] at hkc
        exact (ENat.coe_ne_top k) hkc.symm
      have hle : n₀ ≤ k := hmin k c hk hcu
      have hdvd : algebraMap C B ϖ₀ ∣ algebraMap C B c := by
        rw [← addVal_le_iff_dvd, hv₀, hkc]
        exact_mod_cast hle
      obtain ⟨x, hx⟩ := hdvd
      obtain ⟨q, hq⟩ := hsat c ϖ₀ x hϖ₀0 hx
      have hcq : c = ϖ₀ * q := hinj (by rw [map_mul, hq, hx])
      have hq0 : q ≠ 0 := by rintro rfl; rw [mul_zero] at hcq; exact hc0 hcq
      obtain ⟨j, hj⟩ := hfin q hq0
      have hjq : addVal B (algebraMap C B q) = (j : ℕ∞) := hj
      have hsum : (k : ℕ∞) = (n₀ : ℕ∞) + (j : ℕ∞) := by
        rw [← hkc, ← hv₀, ← hjq, hcq, map_mul, addVal_mul]
      have hkj : k = n₀ + j := by exact_mod_cast hsum
      obtain ⟨n, hn⟩ := ih j (by omega) q hj
      refine ⟨n + 1, ?_⟩
      rw [pow_succ', hcq]
      exact hn.mul_left ϖ₀
  have hirr : Irreducible ϖ₀ := by
    constructor
    · intro hu
      have h0 : f ϖ₀ = 0 := (hzero ϖ₀).2 hu
      rw [hϖ₀] at h0
      have hz : n₀ = 0 := by exact_mod_cast h0
      omega
    · intro a b hab
      by_contra hcon
      have ha : ¬ IsUnit a := fun h => hcon (Or.inl h)
      have hb : ¬ IsUnit b := fun h => hcon (Or.inr h)
      have ha0 : a ≠ 0 := by rintro rfl; rw [zero_mul] at hab; exact hϖ₀0 hab
      have hb0 : b ≠ 0 := by rintro rfl; rw [mul_zero] at hab; exact hϖ₀0 hab
      obtain ⟨i, hi⟩ := hfin a ha0
      obtain ⟨j, hj⟩ := hfin b hb0
      have hia : addVal B (algebraMap C B a) = (i : ℕ∞) := hi
      have hjb : addVal B (algebraMap C B b) = (j : ℕ∞) := hj
      have h1 : n₀ ≤ i := hmin i a hi ha
      have h2 : n₀ ≤ j := hmin j b hj hb
      have hsum : (n₀ : ℕ∞) = (i : ℕ∞) + (j : ℕ∞) := by
        rw [← hv₀, ← hia, ← hjb, hab, map_mul, addVal_mul]
      have hij : n₀ = i + j := by exact_mod_cast hsum
      omega
  refine ofHasUnitMulPowIrreducibleFactorization ⟨ϖ₀, hirr, ?_⟩
  intro x hx
  obtain ⟨k, hk⟩ := hfin x hx
  exact key k x hk

/-! ## §3b 具体層 —— 固定環は割り算で閉じ、非零非単元を持つ -/

/-- **固定環は割り算で閉じている** —— `H` 不変な `c, d`(`d ≠ 0`)の商は `H` 不変。
`ρ` を `c = d·x` に当てて `d·(ρ•x) = d·x`、`B` は整域なので消去できる。 -/
theorem fixedRing_sat {B : Type*} [CommRing B] [IsDomain B] {G : Type*} [Group G]
    [MulSemiringAction G B] {H : Subgroup G} (c d : ↥(fixedRing B H)) (x : B) (hd : d ≠ 0)
    (h : algebraMap (↥(fixedRing B H)) B c = algebraMap (↥(fixedRing B H)) B d * x) :
    ∃ q : ↥(fixedRing B H), algebraMap (↥(fixedRing B H)) B q = x := by
  refine exists_algebraMap_fixedRing x fun ρ hρ => ?_
  have hd0 : algebraMap (↥(fixedRing B H)) B d ≠ 0 := fun hh =>
    hd (fixedRing_injective (by rw [hh, map_zero]))
  refine mul_left_cancel₀ hd0 ?_
  calc algebraMap (↥(fixedRing B H)) B d * (ρ • x)
      = (ρ • algebraMap (↥(fixedRing B H)) B d) * (ρ • x) := by
        rw [smul_algebraMap_fixedRing ρ hρ d]
    _ = ρ • (algebraMap (↥(fixedRing B H)) B d * x) := by rw [smul_mul']
    _ = ρ • algebraMap (↥(fixedRing B H)) B c := by rw [← h]
    _ = algebraMap (↥(fixedRing B H)) B d * x := by
        rw [smul_algebraMap_fixedRing ρ hρ c, h]

/-- **固定環は体でない** —— ノルム元 `ϖ₀ := ∏_{τ∈H} τ•π` は `H` 不変(Y10 `smul_prod_conj`)で
その付値は `|H| > 0`(Y10 `addVal_prod_conj`)だから、非零かつ非単元。

★`[Fintype ↥H]` を落とすとこの補題は使えない(積が作れない)。実際 `H` が無限だと
固定環は体になりうる。 -/
theorem fixedRing_exists_nonzero_nonunit {B : Type*} [CommRing B] [IsDomain B]
    [IsDiscreteValuationRing B] {G : Type*} [Group G] [MulSemiringAction G B]
    {H : Subgroup G} [Fintype H] :
    ∃ c : ↥(fixedRing B H), c ≠ 0 ∧ ¬ IsUnit c := by
  obtain ⟨π, hπ⟩ := exists_irreducible B
  obtain ⟨c, hc⟩ := exists_algebraMap_fixedRing (H := H) (∏ τ : H, (τ : G) • π)
    fun ρ hρ => smul_prod_conj π hρ
  have hval : addVal B (algebraMap (↥(fixedRing B H)) B c) = (Nat.card H : ℕ∞) := by
    rw [hc]; exact addVal_prod_conj hπ
  have hpos : (Nat.card H : ℕ∞) ≠ 0 := by
    exact_mod_cast (Nat.card_pos (α := H)).ne'
  refine ⟨c, ?_, ?_⟩
  · rintro rfl
    rw [map_zero, addVal_zero] at hval
    exact (ENat.coe_ne_top _) hval.symm
  · intro hu
    exact hpos (by
      rw [← hval, addVal_eq_zero_iff]
      exact hu.map (algebraMap (↥(fixedRing B H)) B))

/-- ★★★**固定環 `B^H` は離散付値環である**(`H` 有限)。

★Y10 が `[IsDiscreteValuationRing C]` として仮定していたものが、ここではインスタンスになった。
`hunit` は Y10 の `isUnit_of_isUnit_map`、`hsat` は `fixedRing_sat`、
`hnu` は `fixedRing_exists_nonzero_nonunit` が供給する。

★**完全分岐は要らない。** `H` が有限でありさえすれば `B^H` は DVR である。 -/
instance fixedRing_isDiscreteValuationRing {B : Type*} [CommRing B] [IsDomain B]
    [IsDiscreteValuationRing B] {G : Type*} [Group G] [MulSemiringAction G B]
    (H : Subgroup G) [Fintype H] : IsDiscreteValuationRing ↥(fixedRing B H) :=
  isDiscreteValuationRing_of_saturated fixedRing_injective
    (fun _ hu => isUnit_of_isUnit_map fixedRing_injective smul_algebraMap_fixedRing
      exists_algebraMap_fixedRing hu)
    fixedRing_sat fixedRing_exists_nonzero_nonunit

/-! ## §4 主結論 —— Y10 の債務を全部払った形 -/

/-- ★★★★**Y10 `addVal_map_eq_card_mul` を `C = B^H` に代入した形**。

`∀ c : B^H, addVal B c = |H| * addVal (B^H) c`。

★**仮定は原典 §6.1–§6.2 の底の設定だけになった**:
* `hπ` : `π′` は `𝒪_{K′}` の素元。
* `hresA` : ★**`K′/K` が完全分岐**(剰余体が伸びない)。★落とすと偽。
* `hadjA` : `𝒪_{K′} = 𝒪[π′]`(原典 Lemma 5.11)。
* `[SMulCommClass G A B]` : `G` は `K` 上の作用である。
* `[FaithfulSMul G B]` : ★落とすと偽(Y10 が発見した退化条件)。
* `[Fintype ↥H]` : ★落とすと `Nat.card H = 0` で空虚に真。

★Y10 が持ち回っていた `hres` / `hAC` / `hinj` / `hinvC` / `hfix` / `hadj` /
`[IsDiscreteValuationRing C]` は**すべて消えた**。 -/
theorem addVal_map_eq_card_mul_fixedRing {A B : Type*} [CommRing A] [CommRing B] [IsDomain B]
    [IsDiscreteValuationRing B] [Algebra A B] {G : Type*} [Group G] [MulSemiringAction G B]
    [SMulCommClass G A B] [FaithfulSMul G B] {H : Subgroup G} [Fintype H] {π : B}
    (hπ : Irreducible π)
    (hresA : ∀ b : B, ∃ a : A, b - algebraMap A B a ∈ maximalIdeal B)
    (hadjA : Algebra.adjoin A ({π} : Set B) = ⊤) (c : ↥(fixedRing B H)) :
    addVal B (algebraMap (↥(fixedRing B H)) B c)
      = (Nat.card H : ℕ∞) * addVal (↥(fixedRing B H)) c :=
  addVal_map_eq_card_mul hπ fixedRing_injective smul_algebraMap_fixedRing
    exists_algebraMap_fixedRing (exists_sub_mem_fixedRing hresA)
    (adjoin_fixedRing_eq_top hadjA) c

/-- ★★★**`e(K′/K″) = |H|`** —— 固定環の素元の `𝒪_{K′}` での付値はちょうど `|H|`。
底の設定だけから出る形。 -/
theorem addVal_map_uniformizer_eq_card_fixedRing {A B : Type*} [CommRing A] [CommRing B]
    [IsDomain B] [IsDiscreteValuationRing B] [Algebra A B] {G : Type*} [Group G]
    [MulSemiringAction G B] [SMulCommClass G A B] [FaithfulSMul G B] {H : Subgroup G} [Fintype H]
    {π : B} (hπ : Irreducible π)
    (hresA : ∀ b : B, ∃ a : A, b - algebraMap A B a ∈ maximalIdeal B)
    (hadjA : Algebra.adjoin A ({π} : Set B) = ⊤)
    {ϖ : ↥(fixedRing B H)} (hϖ : Irreducible ϖ) :
    addVal B (algebraMap (↥(fixedRing B H)) B ϖ) = (Nat.card H : ℕ∞) := by
  rw [addVal_map_eq_card_mul_fixedRing (A := A) hπ hresA hadjA ϖ, addVal_uniformizer hϖ, mul_one]

end ABC3.Found.PGC
