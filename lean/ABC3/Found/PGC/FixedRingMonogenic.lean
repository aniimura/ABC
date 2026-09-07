import ABC3.Found.PGC.FixedRingAction

/-!
# `𝒪_{K″} = 𝒪[ϖ]`(原典 Lemma 5.11 を固定環 `K″ = (K′)^H` に当てた形)

**原典**: Yoshida, *Local Class Field Theory via Lubin-Tate Theory* (2008), 物理 p.12,
Lemma 5.11 ——

> If `L′/L` is totally ramified and `α` is a uniformizer of `L′`, then `𝒪_{L′} = 𝒪_L[α]`.

★本ファイルはこれを **`L := K`(底)/ `L′ := K″ = (K′)^H`** に当てる。
一般形の `adjoin_uniformizer_eq_top`(`Found/PGC/LowerRamificationGroup.lean:200`、
中村補題による証明)は既に木にあり、原典 Lemma 5.11 の `.src` は
`exists_uniformizer_adjoin_eq`(同ファイル)に付いている。★本ファイルは**同じ原典項目を
固定環に当てた実例**であり、`.src` を重ねて置く(同一 item への複数 `.src` は
本木で普通に起きている —— `GenEll p17 Proposition 3.4` は 63 件)。

## なぜ立てたか —— Λ7 に残っていた**最後の 1 つ**の仮定

Y16 `Found/PGC/FixedRingAction.lean` の
`card_mul_ramIndex_eq_sum_ramIndex_fixedRing`(原典 Lemma 6.8 形)は、Y10 が持ち回っていた
8 つの債務のうち 7 つを返したが、次の 1 つだけを仮定として残していた:

```
hfix : ∀ y : B, (∀ ρ ∈ H, ρ • y = y) →
  y ∈ Algebra.adjoin A ({algebraMap (↥(fixedRing B H)) B ϖ} : Set B)
```

これは「`H` 不変な元は `𝒪[ϖ]` に入る」、すなわち `𝒪_{K″} ⊆ 𝒪[ϖ]` である。
★**本ファイルがそれを供給する**(`fixedRing_mem_adjoin_uniformizer`)。

## 何を示したか

**§1(抽象核・純可換環)** ★**群も Galois も分岐も 1 語も出てこない。**

* `mem_maximalIdeal_of_map_mem` —— 局所環の間では `𝔪_B` の逆像は `𝔪_C` に入る
  (非単元の逆像が非単元、というだけ)。
* `exists_sub_mem_maximalIdeal_mid` —— 剰余の全射性 `A → B` は塔の中間環 `A → C` に落ちる。
  ★**これが「`K″/K` も完全分岐」の実質**である(原典が Lemma 5.11 の仮定に置いた
  "totally ramified" は、ここでは `hresA` として底から降りてくる)。
* `map_maximalIdeal_ne_bot` —— `𝔪_A` の像が `B` の中で消えなければ `C` の中でも消えない。
* `module_finite_of_algebraMap_injective` —— `IsNoetherian A B` があれば
  `A`-部分加群は有限生成、よって `Module.Finite A C`。
* `adjoin_uniformizer_mid_eq_top` —— ★**Lemma 5.11 の中間環版**。
  `C` が DVR、`ϖ` がその素元、`A` 局所、`Module.Finite A C`、`hresA`、`hAne` から
  `Algebra.adjoin A {ϖ} = ⊤`。Y10 の `adjoin_uniformizer_eq_top` に代入するだけ。
* `mem_adjoin_map_of_adjoin_eq_top` —— ★**adjoin を像へ移す**。
  `AlgHom.map_adjoin` + `Set.image_singleton` の 1 行。
* `Subring.algebraOfMapsTo` —— ★★**`A` の像を含む部分環 `S ⊆ B` には `Algebra A ↥S` が
  `RingHom.codRestrict` で立つ。** 塔 `IsScalarTower A ↥S B` は **`rfl` で出る**。
* `adjoin_uniformizer_subring_eq_top` / `mem_adjoin_uniformizer_of_subring`
  —— ★★★**本ファイルの主抽象核**。「DVR `B` の部分環 `S` が DVR で、`A` の像を含み、
  `B/A` の剰余が伸びないなら `S = A[ϖ]`」。★固定環も群作用も出てこない。

**§2(具体層)** `S := fixedRing B H` を代入するだけ。`[IsDiscreteValuationRing ↥S]` は
Y14 の `fixedRing_isDiscreteValuationRing`(`[Fintype ↥H]` が要る)、
`hS`(`A` の像が固定環に入る)は Y14 の `smul_algebraMap_eq_self`(`[SMulCommClass G A B]`)。
★**`H ⊴ G` は §1・§2 では 1 度も要らない。**

**§3(系)** `card_mul_ramIndex_eq_sum_ramIndex_fixedRing_of_irreducible` ——
Y16 の Lemma 6.8 形から `hfix` が消えた形。★★**Λ7 の鎖
(Y10 → Y11 → Y12 → Y14 → Y16 → Y15)の仮定が、原典 §6.1–§6.2 の底の設定だけに載った。**

## ★逸脱の記録

1. ★**原典の証明経路とは別**である。原典は付値計算(`v(a_iα^i)` が互いに相異なる)で
   基底と整性を同時に出すが、本木は Y10 が既に中村補題で書き換えている
   (`LowerRamificationGroup.lean` 冒頭「逸脱の記録(2)」)。本ファイルはその一般形に
   代入するだけなので、同じ逸脱を引き継ぐ。★結論は同一である。
2. ★★**`[IsNoetherian A B]` を足した**(原典には無い)。原典は `L′/L` が有限次拡大の
   整数環という設定を暗黙に持っており、`𝒪_{L′}` が `𝒪_L` 上有限であることを
   「`{1,α,…,α^{n−1}}` が基底」の形で証明の中で使っている。中村補題経路では
   `Module.Finite A C` が要り、それを `A`-部分加群の有限生成性から出すのに
   `IsNoetherian A B` が要る。★§6.1 の設定(`𝒪_K` は DVR、`𝒪_{K′}` はその上有限)では
   真である。★**弱めていない(仮定を足しただけ)。**
3. ★★**`hAne : ∃ a ∈ maximalIdeal A, algebraMap A B a ≠ 0` を足した**(原典には無い)。
   原典は `L ⊆ L′` という体の拡大なので `𝒪_L → 𝒪_{L′}` は自明に単射である。
   本木は `A` を単なる環として受けている(Y14 の逸脱記録 2)ので、
   「`𝔪_A` が潰れない」を明示する必要がある。★`A` が DVR で `A → B` が単射なら
   `exists_mem_maximalIdeal_map_ne_zero` が自動で供給する(§1 末尾)。
4. ★**既存の `Found/PGC/*.lean` は 1 行も書き換えていない。** §3 の系は新しい名前で足した。
5. ★**`hadjA`(`𝒪_{K′} = 𝒪[π′]`)は §1・§2 では使っていない。**
   ★段取りより安い —— 「`K′/K` の完全分岐(`hresA`)だけで `K″/K` 版の Lemma 5.11 が出る」。
   `hadjA` が要るのは §3(Y16 の Lemma 6.8 形が要求する)だけである。

## 退化の自己検査

* ★★**`hresA`(`K′/K` の完全分岐)を落とすと偽**である。不分岐なら剰余体が伸び、
  `𝒪_{K″} ≠ 𝒪[ϖ]` になる(`𝒪[ϖ]` の剰余体は `k_K` のままで `k_{K″}` に届かない)。
  ★§1 の `exists_sub_mem_maximalIdeal_mid` がその唯一の入口である。
* ★**`ϖ` が `C` の素元(`Irreducible ϖ`)であることを落とすと偽**。
  `ϖ = 0` なら `A[ϖ] = A` の像、`ϖ` が単元なら同じく届かない。
  ★`hϖ` は `maximalIdeal C = span {ϖ}` として 2 箇所で使われる。
* ★**`[Fintype ↥H]` を落とすと `C` が DVR にならない**(Y14 の実測: ノルム元
  `∏_{τ∈H} τ•π` が作れず「`C` は体でない」が言えない)。★完全分岐は要らないが
  有限性は要る。
* ★**`hAne` を落とすと偽**。`𝔪_A` が `B` で潰れると `hpow` が作れず、実際
  `A` の像が体になって `A[ϖ]` は多項式環の像にしかならない。
* ★**`[IsNoetherian A B]` を落とすと中村補題が使えない**(有限生成性が要る)。
* ★**`H ⊴ G` は §1・§2 では要らない** —— 固定環は正規性なしに部分環である。
  §3 では Y16 の作用が要るので `[H.Normal]` が復活する。
* ★`ℕ∞` の切り詰め引き算・除算は 1 箇所も書いていない(`lean-idioms.md` #102)。
-/

namespace ABC3.Found.PGC

open IsLocalRing IsDiscreteValuationRing

/-! ## §1 抽象核 —— 純可換環(群・Galois・分岐の語彙は 1 つも出てこない) -/

section AbstractCore

variable {A C B : Type*}

/-- ★**抽象核** —— 局所環の間の代数写像では `𝔪_B` の逆像が `𝔪_C` に入る。

非単元の逆像は非単元、というだけである。★`C` が局所であることが本質
(そうでなければ「非単元 ⇒ 極大イデアルの元」が言えない)。 -/
theorem mem_maximalIdeal_of_map_mem [CommRing C] [IsLocalRing C] [CommRing B] [IsLocalRing B]
    [Algebra C B] {x : C} (h : algebraMap C B x ∈ maximalIdeal B) : x ∈ maximalIdeal C := by
  rw [IsLocalRing.mem_maximalIdeal] at h ⊢
  exact fun hu => h (hu.map (algebraMap C B))

/-- ★★**抽象核** —— 剰余の全射性は塔の**中間環**に落ちる。

`A → B` が剰余体を伸ばさない(`hresA`)なら `A → C` も伸ばさない。

★★**これが「`K′/K` が完全分岐なら `K″/K` も完全分岐」の実質**である。
原典 Lemma 5.11 は `L′/L` が完全分岐であることを仮定に置くが、本木では底の
`hresA` から中間環へ降ろす形になる。★分岐・付値・Galois の語彙は出てこない。 -/
theorem exists_sub_mem_maximalIdeal_mid [CommRing A] [CommRing C] [IsLocalRing C] [CommRing B]
    [IsLocalRing B] [Algebra A C] [Algebra A B] [Algebra C B] [IsScalarTower A C B]
    (hresA : ∀ b : B, ∃ a : A, b - algebraMap A B a ∈ maximalIdeal B) (c : C) :
    ∃ a : A, c - algebraMap A C a ∈ maximalIdeal C := by
  obtain ⟨a, ha⟩ := hresA (algebraMap C B c)
  refine ⟨a, mem_maximalIdeal_of_map_mem (B := B) ?_⟩
  rwa [map_sub, ← IsScalarTower.algebraMap_apply A C B]

/-- ★**抽象核** —— `𝔪_A` の像が `B` の中で消えないなら、中間環 `C` の中でも消えない。

★`adjoin_uniformizer_eq_top` の `hpow` を `exists_pow_maximalIdeal_le_map` 経由で
供給するための非退化条件である。 -/
theorem map_maximalIdeal_ne_bot [CommRing A] [IsLocalRing A] [CommRing C] [CommRing B]
    [Algebra A C] [Algebra A B] [Algebra C B] [IsScalarTower A C B]
    (hAne : ∃ a ∈ maximalIdeal A, algebraMap A B a ≠ 0) :
    Ideal.map (algebraMap A C) (maximalIdeal A) ≠ ⊥ := by
  obtain ⟨a, ha, ha0⟩ := hAne
  intro h
  refine ha0 ?_
  have hmem : algebraMap A C a ∈ Ideal.map (algebraMap A C) (maximalIdeal A) :=
    Ideal.mem_map_of_mem _ ha
  rw [h, Ideal.mem_bot] at hmem
  rw [IsScalarTower.algebraMap_apply A C B, hmem, map_zero]

/-- ★**抽象核** —— ネーター的な `B` の中の中間環は `A` 上有限である。

`C → B` が単射な `A`-線形写像だから `Module.Finite.of_injective` が効く。
★これが中村補題(`adjoin_uniformizer_eq_top` の `[Module.Finite A C]`)への入口である。 -/
theorem module_finite_of_algebraMap_injective [CommRing A] [CommRing C] [CommRing B]
    [Algebra A C] [Algebra A B] [Algebra C B] [IsScalarTower A C B] [IsNoetherian A B]
    (hinj : Function.Injective (algebraMap C B)) : Module.Finite A C :=
  Module.Finite.of_injective (IsScalarTower.toAlgHom A C B).toLinearMap hinj

/-- ★★★**抽象核 —— Yoshida Lemma 5.11(中間環版)**。`C = A[ϖ]`。

塔 `A → C → B` で `C` が DVR、`ϖ` がその素元、`A` が局所、`C` が `A` 上有限、
底 `A → B` の剰余が伸びない(`hresA`)、`𝔪_A` が潰れない(`hAne`)なら
`Algebra.adjoin A {ϖ} = ⊤`。

★Y10 の `adjoin_uniformizer_eq_top`(中村補題)に §1 の 3 本を代入するだけである。
★**分岐・Galois・群作用は 1 語も出てこない。** -/
theorem adjoin_uniformizer_mid_eq_top [CommRing A] [IsLocalRing A] [CommRing C] [IsDomain C]
    [IsDiscreteValuationRing C] [CommRing B] [IsLocalRing B] [Algebra A C] [Algebra A B]
    [Algebra C B] [IsScalarTower A C B] [Module.Finite A C] {ϖ : C} (hϖ : Irreducible ϖ)
    (hresA : ∀ b : B, ∃ a : A, b - algebraMap A B a ∈ maximalIdeal B)
    (hAne : ∃ a ∈ maximalIdeal A, algebraMap A B a ≠ 0) :
    Algebra.adjoin A ({ϖ} : Set C) = ⊤ :=
  adjoin_uniformizer_eq_top ((irreducible_iff_uniformizer ϖ).mp hϖ)
    (exists_sub_mem_maximalIdeal_mid (B := B) hresA)
    (exists_pow_maximalIdeal_le_map ((irreducible_iff_uniformizer ϖ).mp hϖ)
      (map_maximalIdeal_ne_bot (B := B) hAne))

/-- ★**抽象核** —— `Algebra.adjoin A S = ⊤` は環準同型の像へ移る。

`A[ϖ] = C` なら `ι c ∈ A[ι ϖ]`(`ι = algebraMap C B`)。
`AlgHom.map_adjoin` + `Set.image_singleton` の 1 行。★純可換環である。 -/
theorem mem_adjoin_map_of_adjoin_eq_top [CommRing A] [CommRing C] [CommRing B]
    [Algebra A C] [Algebra A B] [Algebra C B] [IsScalarTower A C B] {ϖ : C}
    (hadj : Algebra.adjoin A ({ϖ} : Set C) = ⊤) (c : C) :
    algebraMap C B c ∈ Algebra.adjoin A ({algebraMap C B ϖ} : Set B) := by
  have h : algebraMap C B c ∈
      (Algebra.adjoin A ({ϖ} : Set C)).map (IsScalarTower.toAlgHom A C B) :=
    Subalgebra.mem_map.2 ⟨c, hadj ▸ Algebra.mem_top, rfl⟩
  rwa [AlgHom.map_adjoin, Set.image_singleton, IsScalarTower.coe_toAlgHom'] at h

/-- ★★**抽象核** —— `A` の像を含む部分環 `S ⊆ B` には `Algebra A ↥S` が立つ。

`RingHom.codRestrict` で `A →+* ↥S` を作り `RingHom.toAlgebra` で束ねるだけ。
★★**塔 `IsScalarTower A ↥S B` は `rfl` で出る**(実測)。

★これを**大域インスタンスにはしない**(証明 `h` を引数に取るので不可能でもある)。
消費側では `letI` で局所的に入れる。★本ファイルの主結論の statement には
`Algebra A ↥S` が現れないので、それで困らない。 -/
@[reducible]
def Subring.algebraOfMapsTo [CommRing A] [CommRing B] [Algebra A B] (S : Subring B)
    (h : ∀ a : A, algebraMap A B a ∈ S) : Algebra A ↥S :=
  ((algebraMap A B).codRestrict S h).toAlgebra

/-- ★★★**抽象核 —— Yoshida Lemma 5.11(部分環版)**。

離散付値環 `B` の部分環 `S` が
(i) それ自身 DVR で、(ii) `A` の像を含み、(iii) 底 `A → B` の剰余が伸びず、
(iv) `𝔪_A` が潰れず、(v) `B` が `A` 上ネーターなら、`S = A[ϖ]`(`ϖ` は `S` の素元)。

★**固定環も群作用も分岐も出てこない。** これが原典 Lemma 5.11 の、この木での一般形である
(`L := Frac A`、`L′ := Frac S`)。 -/
theorem adjoin_uniformizer_subring_eq_top [CommRing A] [IsLocalRing A] [CommRing B] [IsDomain B]
    [IsDiscreteValuationRing B] [Algebra A B] [IsNoetherian A B] (S : Subring B)
    [IsDiscreteValuationRing ↥S] (hS : ∀ a : A, algebraMap A B a ∈ S)
    (hresA : ∀ b : B, ∃ a : A, b - algebraMap A B a ∈ maximalIdeal B)
    (hAne : ∃ a ∈ maximalIdeal A, algebraMap A B a ≠ 0)
    {ϖ : ↥S} (hϖ : Irreducible ϖ) :
    letI : Algebra A ↥S := Subring.algebraOfMapsTo S hS
    Algebra.adjoin A ({ϖ} : Set ↥S) = ⊤ := by
  letI : Algebra A ↥S := Subring.algebraOfMapsTo S hS
  haveI : IsScalarTower A ↥S B := IsScalarTower.of_algebraMap_eq fun _ => rfl
  haveI : Module.Finite A ↥S :=
    module_finite_of_algebraMap_injective (B := B) Subtype.coe_injective
  exact adjoin_uniformizer_mid_eq_top (B := B) hϖ hresA hAne

/-- **原典 Lemma 5.11**(Yoshida 2008, 物理 p.12)。

> If `L′/L` is totally ramified and `α` is a uniformizer of `L′`, then `𝒪_{L′} = 𝒪_L[α]`.

★同じ原典項目の一般形は `exists_uniformizer_adjoin_eq`
(`Found/PGC/LowerRamificationGroup.lean`)に既にある。本宣言はそれを
「`B` の部分環 `S`」の形に当てたもので、固定環 `K″ = (K′)^H` への適用がその用途である。 -/
def adjoin_uniformizer_subring_eq_top.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 12, item := "Lemma 5.11", sectionId := "lemma-5-11" }

/-- ★★★**抽象核(主)** —— 部分環 `S` の元は `A[ϖ]` に入る。

`adjoin_uniformizer_subring_eq_top` を像へ移した形。★これが `hfix` の抽象核である。 -/
theorem mem_adjoin_uniformizer_of_subring [CommRing A] [IsLocalRing A] [CommRing B] [IsDomain B]
    [IsDiscreteValuationRing B] [Algebra A B] [IsNoetherian A B] (S : Subring B)
    [IsDiscreteValuationRing ↥S] (hS : ∀ a : A, algebraMap A B a ∈ S)
    (hresA : ∀ b : B, ∃ a : A, b - algebraMap A B a ∈ maximalIdeal B)
    (hAne : ∃ a ∈ maximalIdeal A, algebraMap A B a ≠ 0)
    {ϖ : ↥S} (hϖ : Irreducible ϖ) (c : ↥S) :
    algebraMap (↥S) B c ∈ Algebra.adjoin A ({algebraMap (↥S) B ϖ} : Set B) := by
  letI : Algebra A ↥S := Subring.algebraOfMapsTo S hS
  haveI : IsScalarTower A ↥S B := IsScalarTower.of_algebraMap_eq fun _ => rfl
  haveI : Module.Finite A ↥S :=
    module_finite_of_algebraMap_injective (B := B) Subtype.coe_injective
  exact mem_adjoin_map_of_adjoin_eq_top
    (adjoin_uniformizer_mid_eq_top (B := B) hϖ hresA hAne) c

/-- ★**抽象核** —— `A` が DVR で `A → B` が単射なら `hAne` は自動である。

★原典の設定(`L ⊆ L′` は体の拡大)ではこれが常に成り立つ。 -/
theorem exists_mem_maximalIdeal_map_ne_zero [CommRing A] [IsDomain A]
    [IsDiscreteValuationRing A] [CommRing B] [Algebra A B]
    (hinj : Function.Injective (algebraMap A B)) :
    ∃ a ∈ maximalIdeal A, algebraMap A B a ≠ 0 := by
  obtain ⟨ϖ, hϖ⟩ := exists_irreducible A
  refine ⟨ϖ, (IsLocalRing.mem_maximalIdeal ϖ).2 hϖ.not_isUnit, fun h => hϖ.ne_zero ?_⟩
  exact hinj (by rw [h, map_zero])

end AbstractCore

/-! ## §2 具体層 —— 固定環 `C = B^H` に代入する -/

/-- ★**`A` の像は固定環に入る** —— `[SMulCommClass G A B]` から Y14 の
`smul_algebraMap_eq_self` 1 行。★これが §1 の `hS` である。 -/
theorem algebraMap_mem_fixedRing {A B : Type*} [CommRing A] [CommRing B] [Algebra A B]
    {G : Type*} [Group G] [MulSemiringAction G B] [SMulCommClass G A B] {H : Subgroup G}
    (a : A) : algebraMap A B a ∈ fixedRing B H :=
  mem_fixedRing.2 fun ρ _ => smul_algebraMap_eq_self ρ a

/-- ★★★**`𝒪_{K″} = 𝒪[ϖ]`** —— 原典 Lemma 5.11 を `K″ = (K′)^H` に当てた形。

★仮定は原典 §6.1 の底の設定だけである:
* `[IsLocalRing A]` : `𝒪 = 𝒪_K` は局所環。
* `[IsNoetherian A B]` : ★逸脱記録 2(原典には無いが §6.1 では真)。
* `hresA` : ★★**`K′/K` が完全分岐**(剰余体が伸びない)。★落とすと偽。
* `hAne` : ★逸脱記録 3(`𝔪_K` が `𝒪_{K′}` で潰れない)。
* `hϖ` : `ϖ` は `𝒪_{K″}` の素元。★落とすと偽。
* `[Fintype ↥H]` : ★落とすと `𝒪_{K″}` が DVR にならない(Y14)。

★**`H ⊴ G` は要らない**(固定環は正規性なしに部分環である)。
★**`hadjA`(`𝒪_{K′} = 𝒪[π′]`)も要らない** —— 段取りより安かった。 -/
theorem adjoin_fixedRing_uniformizer_eq_top {A B : Type*} [CommRing A] [IsLocalRing A]
    [CommRing B] [IsDomain B] [IsDiscreteValuationRing B] [Algebra A B] [IsNoetherian A B]
    {G : Type*} [Group G] [MulSemiringAction G B] [SMulCommClass G A B] {H : Subgroup G}
    [Fintype H]
    (hresA : ∀ b : B, ∃ a : A, b - algebraMap A B a ∈ maximalIdeal B)
    (hAne : ∃ a ∈ maximalIdeal A, algebraMap A B a ≠ 0)
    {ϖ : ↥(fixedRing B H)} (hϖ : Irreducible ϖ) :
    letI : Algebra A ↥(fixedRing B H) :=
      Subring.algebraOfMapsTo (fixedRing B H) algebraMap_mem_fixedRing
    Algebra.adjoin A ({ϖ} : Set ↥(fixedRing B H)) = ⊤ :=
  adjoin_uniformizer_subring_eq_top (fixedRing B H) algebraMap_mem_fixedRing hresA hAne hϖ

/-- ★★★★**Y16 の債務 `hfix` を供給する。**

```
hfix : ∀ y : B, (∀ ρ ∈ H, ρ • y = y) → y ∈ 𝒪[ϖ]
```

`H` 不変な `y` は Y14 の `exists_algebraMap_fixedRing` で `y = ι c`(`c : 𝒪_{K″}`)と書け、
`𝒪_{K″} = 𝒪[ϖ]`(上)を像へ移せばよい。

★★**これが Λ7 の鎖(Y10 → Y11 → Y12 → Y14 → Y16 → Y15)に残っていた最後の 1 つである。** -/
theorem fixedRing_mem_adjoin_uniformizer {A B : Type*} [CommRing A] [IsLocalRing A]
    [CommRing B] [IsDomain B] [IsDiscreteValuationRing B] [Algebra A B] [IsNoetherian A B]
    {G : Type*} [Group G] [MulSemiringAction G B] [SMulCommClass G A B] {H : Subgroup G}
    [Fintype H]
    (hresA : ∀ b : B, ∃ a : A, b - algebraMap A B a ∈ maximalIdeal B)
    (hAne : ∃ a ∈ maximalIdeal A, algebraMap A B a ≠ 0)
    {ϖ : ↥(fixedRing B H)} (hϖ : Irreducible ϖ) :
    ∀ y : B, (∀ ρ ∈ H, ρ • y = y) →
      y ∈ Algebra.adjoin A ({algebraMap (↥(fixedRing B H)) B ϖ} : Set B) := by
  intro y hy
  obtain ⟨c, rfl⟩ := exists_algebraMap_fixedRing y hy
  exact mem_adjoin_uniformizer_of_subring (fixedRing B H) algebraMap_mem_fixedRing hresA hAne hϖ c

/-! ## §3 系 —— Y16 の Lemma 6.8 形から `hfix` を消す -/

/-- ★★★★★**原典 Lemma 6.8(Yoshida §6.2)—— 仮定が底の設定だけになった形**。

```
|H| · i_{K″}(σ ϖ) = Σ_{τ ∈ H} i_{K′}(σ τ)
```

★Y16 `card_mul_ramIndex_eq_sum_ramIndex_fixedRing` から **`hfix` が消えた**。
残る仮定は原典 §6.1–§6.2 の設定そのものである:

* `hπ'` : `π′` は `𝒪_{K′}` の素元。
* `hresA` : ★★**`K′/K` が完全分岐**。★落とすと偽。
* `hadjA` : `𝒪_{K′} = 𝒪[π′]`(原典 Lemma 5.11 を `K′/K` に当てた形)。
* `hAne` : ★逸脱記録 3。
* `hϖ` : `ϖ` は `𝒪_{K″}` の素元。★★**Y16 では任意の `ϖ` でよかったが、
  `hfix` を供給するには素元であることが要る**(`ϖ` が単元なら `𝒪_{K″} ≠ 𝒪[ϖ]`)。
* `[IsNoetherian A B]` : ★逸脱記録 2。
* `[SMulCommClass G A B]` / `[FaithfulSMul G B]` / `[H.Normal]` / `[Fintype ↥H]`。 -/
theorem card_mul_ramIndex_eq_sum_ramIndex_fixedRing_of_irreducible {A B : Type*} [CommRing A]
    [IsLocalRing A] [CommRing B] [Algebra A B] [IsDomain B] [IsDiscreteValuationRing B]
    [IsNoetherian A B] {G : Type*} [Group G] [MulSemiringAction G B] [SMulCommClass G A B]
    [FaithfulSMul G B] {H : Subgroup G} [H.Normal] [Fintype H]
    {π' : B} (hπ' : Irreducible π')
    (hresA : ∀ b : B, ∃ a : A, b - algebraMap A B a ∈ maximalIdeal B)
    (hadjA : Algebra.adjoin A ({π'} : Set B) = ⊤)
    (hAne : ∃ a ∈ maximalIdeal A, algebraMap A B a ≠ 0)
    {ϖ : ↥(fixedRing B H)} (hϖ : Irreducible ϖ) (σ : G) :
    (Nat.card H : ℕ∞) * ramIndex ϖ σ = ∑ τ : H, ramIndex π' (σ * (τ : G)) :=
  card_mul_ramIndex_eq_sum_ramIndex_fixedRing hπ' hresA hadjA ϖ
    (fixedRing_mem_adjoin_uniformizer hresA hAne hϖ) σ

/-- ★**上の系の使いやすい形** —— `A` が DVR で `A → B` が単射なら `hAne` は自動。

★原典の設定(`K ⊆ K′` は体の拡大、`𝒪_K` は DVR)ではこれが常に成り立つので、
実質的に `hAne` は消えている。 -/
theorem card_mul_ramIndex_eq_sum_ramIndex_fixedRing_of_injective {A B : Type*} [CommRing A]
    [IsDomain A] [IsDiscreteValuationRing A] [CommRing B] [Algebra A B] [IsDomain B]
    [IsDiscreteValuationRing B] [IsNoetherian A B] {G : Type*} [Group G] [MulSemiringAction G B]
    [SMulCommClass G A B] [FaithfulSMul G B] {H : Subgroup G} [H.Normal] [Fintype H]
    {π' : B} (hπ' : Irreducible π')
    (hresA : ∀ b : B, ∃ a : A, b - algebraMap A B a ∈ maximalIdeal B)
    (hadjA : Algebra.adjoin A ({π'} : Set B) = ⊤)
    (hAinj : Function.Injective (algebraMap A B))
    {ϖ : ↥(fixedRing B H)} (hϖ : Irreducible ϖ) (σ : G) :
    (Nat.card H : ℕ∞) * ramIndex ϖ σ = ∑ τ : H, ramIndex π' (σ * (τ : G)) :=
  card_mul_ramIndex_eq_sum_ramIndex_fixedRing_of_irreducible hπ' hresA hadjA
    (exists_mem_maximalIdeal_map_ne_zero hAinj) hϖ σ

end ABC3.Found.PGC
