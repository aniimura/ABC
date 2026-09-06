import ABC3.Found.PGC.RamificationQuotientEmbedding

/-!
# 跳びの割り切り `e_0 ∣ n` と `G_1` の p 群性(Yoshida 2008 Corollary 6.3 とその周辺)

典拠: T. Yoshida, *Local Class Field Theory via Lubin-Tate Theory*
(Ann. Fac. Sci. Toulouse 17-2, 2008; arXiv math/0606108) **Corollary 6.3**(物理 p.14)。
構造化済み原文は `ResearchPaper/1_Structured/Local Class Field Theory via Lubin-Tate Theory/`
`section-6.html` の `id="cor-6-3"`。原文は

> Corollary 6.3. If G is abelian and G_n ≠ G_{n+1}, then e_0 := |G_0/G_1| divides n.

★★**原文を `.txt` で読んではならない**: `pdftotext` は `≠` の斜線(ベクター描画)を落とすので
`G_n = G_{n+1}` と出力し、**主張が反転して見える**。構造化済み HTML は
`data-txt="="` でこの脱落を明示している。同じ脱落は `Lemma 6.5` の `x_n ≠ 0`、
`Prop 6.6` の `⟨σ⟩`(角括弧が消える)にも出る。

## 何を証明したか

**(A) Corollary 6.3 そのもの**

* `nat_card_quot_zero_dvd_of_ne`: `G` 可換・`G_n ≠ G_{n+1}` ⟹ `|G_0/G_1| ∣ n`。
* `nat_card_quot_top_dvd_of_ne`: 完全分岐(`G_0 = ⊤`)版 `|G/G_1| ∣ n`。
  ★Hasse–Arf(Theorem 6.11)が実際に使うのはこちらの形。
* `nat_card_quot_lowerRamificationGroupAdjoin_one_dvd`: `PAdicLocalField` への具体化。

**(B) `Prop 6.6 (Sen)` が "by Proposition 6.2" で畳んでいる中身**

* `pow_char_mem_lowerRamificationGroup_succ`: `σ ∈ G_n`(`n ≥ 1`)⟹ `σ^p ∈ G_{n+1}`。
* `pow_char_eq_one_quot_lowerRamificationGroup`: `G_n/G_{n+1}` は**指数 `p`**(`n ≥ 1`)。
* `isPGroup_lowerRamificationGroup_one`: `G_1` は **`p` 群**。
* `exists_orderOf_eq_pow_of_mem_lowerRamificationGroup_one` /
  `exists_nat_card_zpowers_eq_pow`: `σ ∈ G_1` なら `|⟨σ⟩| = p^m`
  ——`Prop 6.6 (Sen)` の冒頭 "`|⟨σ⟩| = p^m` for `m ≥ 1` (by Proposition 6.2)" の
  `(by Proposition 6.2)` がこれである(`m ≥ 1` の部分は `σ ≠ 1` から出るので別問題、下記参照)。

## 証明の骨格(原文の 1 行を分解したもの)

原文は Prop 6.2 から Cor 6.3 を「immediately」で畳んでいる。中身は次の 3 段:

1. **共役の変換則**(`residue_ramCoeff_conj`)——`σ ∈ G_n`、`τ ∈ G_0` に対し
   `θ_n(τστ⁻¹) = θ_0(τ)^n · θ_n(σ)`。
   ★これが本ノードの心臓で、証明は **`Prop 6.2` の素元取り替えの捻れ**
   `residue_ramCoeff_change_uniformizer`(`c̄^α = w̄^n c̄^β`)**に完全に帰着する**:
   `β := τα` は再び素元(`span_smul_uniformizer`)で、割り算の一意性から
   `c^{τα}_{τστ⁻¹} = τ(c^α_σ)`(`ramCoeff_smul_conj`)。`τ ∈ G_0` は剰余体に自明に作用する
   ので `residue (τ c) = residue c`。捻れの因子 `w̄^n` がそのまま `θ_0(τ)^n` である。
2. **可換性を入れる**——`τστ⁻¹ = σ` なので `θ_n(σ) = θ_0(τ)^n θ_n(σ)`。
   `σ ∉ G_{n+1}` から `θ_n(σ) ≠ 0`(`Prop 6.2` の核の計算)、体で約分して `θ_0(τ)^n = 1`。
3. **`𝓀ˣ` の有限部分群は巡回**——`θ_0 : G_0/G_1 ↪ 𝓀ˣ` は単射だから `G_0/G_1` は巡回群
   (`isCyclic_of_injective_units`)で、指数 = 位数。全ての元が `n` 乗で 1 になるので
   `Monoid.exponent = |G_0/G_1| ∣ n`(`nat_card_dvd_of_injective_units`)。

★**`ramCoeff_mul` は本ノードでは直接は使っていない**。Y3 の見立ては
「`ramCoeff_mul` が準同型性の全情報を持つ」だったが、Cor 6.3 が要求するのは
**準同型性ではなく共変性**(`τ` で捻ったときの変換則)であり、その情報は
`ramCoeff_mul` ではなく `ramCoeff_independent`(素元取り替え)の側に入っていた。
`ramCoeff_mul` は `thetaMul`/`thetaAdd` が準同型であること(= (3) と (B) が使う)を
通じて間接的に効いている。

## ★逸脱の記録

1. **`n ≥ 1` を仮定していない**。原文も Cor 6.3 では `n ≥ 1` と書いていないし、
   上の証明は `n = 0` でもそのまま通る(そのとき `θ_0(τ)^0 = 1` で結論は
   `e_0 ∣ 0`、すなわち**空虚に真**)。仮定を足すと弱くなるだけなので足していない。
   ★退化検査としては「`n = 0` で主張が空虚になる」ことを**明示的に記録する**にとどめる。
2. **`G` の可換性は `∀ x y : G, x * y = y * x` という命題で渡す**(`CommGroup` 構造を
   要求しない)。具体層の `G = Gal(K(x)/K)` は `Group` インスタンスしか持たないので、
   構造で要求すると当てられない。数学的には同じ。
3. **`e_0` を `Nat.card` で書き、有限性は型クラス仮定 `[Finite _]` として置く**。
   抽象核は `G` の有限性を仮定しないので(`Prop 6.2` の埋め込み自体は不要)、
   商が有限であることだけを要求する。具体層では `Finite (L ≃ₐ[K] L)` から自動で出る。
4. **(B) の `m ≥ 1`(原文 `Prop 6.6` の "for m ≥ 1")は出していない**。
   `|⟨σ⟩| = p^m` は出したが `m ≥ 1` は `σ ≠ 1` と同値であり、原文が
   `σ ∈ G_1` に暗黙に課している「非自明」の条件に対応する。使う側で
   `σ ≠ 1 → orderOf σ ≠ 1` として足す方が安いので、ここでは `∃ m` のみを出す。
5. **剰余体の標数**は抽象核では `[CharP (ResidueField B) p]` として仮定する。
   具体層 `PAdicLocalField p` では `charP_residueField_adjoinIntegers` で供給する
   (`𝓀[K.carrier] → 𝓀_{K(x)}` が体の射なので単射、`charP_of_injective_algebraMap`)。

## ★退化の自己検査

* **可換性は落とせない**(原典が `G` abelian を明示的に仮定している)。
  形式化すると使いどころが 1 箇所に特定できる: (2) の `τστ⁻¹ = σ` だけである。
  可換性を外すと残るのは `residue_ramCoeff_conj`(= 共役の変換則)そのもので、
  そこから言えるのは「`G_0/G_1` が `G_n/G_{n+1}` に `θ_0(τ)^n` 倍で作用する」までであり、
  `θ_0(τ)^n = 1` は出ない。★本ファイルは可換性を**使わない**部分
  (`residue_ramCoeff_conj` まで)と**使う**部分(`residue_pow_eq_one_of_notMem` 以降)を
  宣言として分けてあるので、境界が読める。
* **`n = 0` では主張が空虚**(`e_0 ∣ 0` は常に真)。上の逸脱 1。
* **`σ ∉ G_{n+1}` を落とすと偽**。落とすと `θ_n(σ) = 0` が許され、
  (2) の約分ができない——`e_0 ∣ n` は任意の `n` に対して主張できてしまう。
  ここが `G_n ≠ G_{n+1}` という仮定の全部である。
* **単射性(= Lemma 5.11 `B = A[α]`)を落とすと (3) が壊れる**。`G_0/G_1` の巡回性も
  位数と指数の一致も、`θ_0` が単射でなければ言えない。`hadj` はそのために要る。
* **(B) `G_1` の p 群性は完全分岐が要る**。抽象核は「`𝔪_B = (α)` かつ `B = A[α]`」
  という形で完全分岐を吸収している(`hadj`)。不分岐部分があると `G_0 ≠ G` になり、
  `G_1` を測る `θ_n`(`n ≥ 1`)の列が `G` 全体を尽くさない。
* **(B) `[CharP (ResidueField B) p]` を落とすと偽**。`p` を勝手な素数にすると
  `𝓀⁺` が指数 `p` を持たなくなり、`σ^p ∈ G_{n+1}` が言えない。

## 実測

前ノードと同じく**抽象核と具体層を分ける**方針が効いた(6 回連続)。
抽象核 `§1`–`§4` の `lean_check` は 0.09–0.34 秒、`PAdicLocalField` 具体層 `§5` は
0.30–0.78 秒。★見立て 220–450 行に対して**証明本体は約 110 行**で、(A) と (B) の
両方が入った。理由は `Prop 6.2` の在庫(特に `residue_ramCoeff_change_uniformizer`)が
共役の変換則をそのまま与えたこと。
-/

namespace ABC3.Found.PGC

open IsLocalRing ABC3.Skeleton.PGC
open scoped NNReal Valued

def nat_card_quot_zero_dvd_of_ne.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 14, item := "Corollary 6.3", sectionId := "cor-6-3" }

/-! ## §1 抽象核(群論)——体の単数群に埋め込める有限群

`θ_0 : G_0/G_1 ↪ 𝓀ˣ` から「巡回群」と「指数 = 位数」を取り出す部分。分岐は一切現れない。 -/

/-- **整域の単数群に単射に埋め込める有限群は巡回群**。
`𝓀ˣ` の有限部分群が巡回であること(`isCyclic_subgroup_units`)を像に当て、
`MonoidHom.ofInjective` で戻す。 -/
theorem isCyclic_of_injective_units {Q R : Type*} [Group Q] [Finite Q] [CommRing R] [IsDomain R]
    (f : Q →* Rˣ) (hf : Function.Injective f) : IsCyclic Q := by
  haveI : Finite f.range := Finite.of_surjective _ f.rangeRestrict_surjective
  haveI : IsCyclic f.range := isCyclic_subgroup_units f.range
  exact isCyclic_of_surjective (MonoidHom.ofInjective hf).symm
    (MonoidHom.ofInjective hf).symm.surjective

/-- **Corollary 6.3 の群論部分**——`Q ↪ Rˣ` で全ての元が `n` 乗して `1` なら `|Q| ∣ n`。

★**巡回性は落とせない**: 一般の有限アーベル群なら「指数 ∣ n」しか出ず、
例えば `Q = (ℤ/2)²`、`n = 2` で `|Q| = 4 ∤ 2`。`Rˣ` への埋め込みが
`Q` を巡回にしているからこそ指数と位数が一致する。 -/
theorem nat_card_dvd_of_injective_units {Q R : Type*} [Group Q] [Finite Q] [CommRing R]
    [IsDomain R] (f : Q →* Rˣ) (hf : Function.Injective f) {n : ℕ} (hpow : ∀ q : Q, q ^ n = 1) :
    Nat.card Q ∣ n := by
  haveI : IsCyclic Q := isCyclic_of_injective_units f hf
  rw [← IsCyclic.exponent_eq_card]
  exact Monoid.exponent_dvd_iff_forall_pow_eq_one.mpr hpow

/-! ## §2 共役の変換則 `θ_n(τστ⁻¹) = θ_0(τ)^n · θ_n(σ)`

★本ノードの心臓。原文が `Prop 6.2` から `Cor 6.3` へ渡るときに畳んでいる 1 行。 -/

/-- `τ ∈ G_0` なら `τα` も素元。`τα = α(1 + c^0_τ)` で `1 + c^0_τ` は単元
(`isUnit_one_add_ramCoeff_zero`)。 -/
theorem span_smul_uniformizer {B : Type*} [CommRing B] [IsDomain B] [IsLocalRing B] {G : Type*}
    [Group G] [MulSemiringAction G B] {α : B} (hα : maximalIdeal B = Ideal.span {α}) (hα0 : α ≠ 0)
    {τ : G} (hτ : τ ∈ lowerRamificationGroup B G 0) :
    maximalIdeal B = Ideal.span {τ • α} := by
  rw [← mul_one_add_ramCoeff_zero hα hτ,
    Ideal.span_singleton_mul_right_unit (isUnit_one_add_ramCoeff_zero hα hα0 hτ) α]
  exact hα

/-- **共役と素元取り替えの整合**——`c^{τα}_{τστ⁻¹} = τ(c^α_σ)`。

`σα − α = α^{n+1} c` に `τ` を施すと `(τστ⁻¹)(τα) − τα = (τα)^{n+1} τ(c)` であり、
整域での割り算の一意性(`divOfDvd_mul_cancel`)からこれが `c^{τα}_{τστ⁻¹}` そのもの。
★`τ` は任意の群元でよい(`G_0` に入っている必要すらない)。 -/
theorem ramCoeff_smul_conj {B : Type*} [CommRing B] [IsDomain B] [IsLocalRing B] {G : Type*}
    [Group G] [MulSemiringAction G B] {α : B} (hα : maximalIdeal B = Ideal.span {α}) (hα0 : α ≠ 0)
    {n : ℕ} {σ τ : G} (hσ : σ ∈ lowerRamificationGroup B G n) :
    ramCoeff (τ • α) n (τ * σ * τ⁻¹) = τ • ramCoeff α n σ := by
  have h1 : (τ * σ * τ⁻¹) • (τ • α) = τ • (σ • α) := by
    rw [← mul_smul, ← mul_smul]
    congr 1
    group
  have key : (τ * σ * τ⁻¹) • (τ • α) - τ • α = (τ • α) ^ (n + 1) * (τ • ramCoeff α n σ) := by
    rw [h1, ← smul_sub, ← pow_mul_ramCoeff hα hσ, smul_mul', smul_pow']
  rw [ramCoeff, key]
  exact divOfDvd_mul_cancel (pow_ne_zero _ ((smul_ne_zero_iff_ne τ).mpr hα0)) _

/-- ★★**共役の変換則**——`θ_n(τστ⁻¹) = θ_0(τ)^n · θ_n(σ)`(`σ ∈ G_n`、`τ ∈ G_0`)。

証明は 3 本の在庫の合成だけ:
`span_smul_uniformizer`(`τα` も素元)+ `residue_ramCoeff_change_uniformizer`
(素元取り替えの捻れ `c̄^α = w̄^n c̄^β`、`w := 1 + c^0_τ`)+ `ramCoeff_smul_conj`。
★捻れの指数 `n` が、そのまま結論の `θ_0(τ)^n` の指数になる。 -/
theorem residue_ramCoeff_conj {B : Type*} [CommRing B] [IsDomain B] [IsLocalRing B] {G : Type*}
    [Group G] [MulSemiringAction G B] {α : B} (hα : maximalIdeal B = Ideal.span {α}) (hα0 : α ≠ 0)
    {n : ℕ} {σ τ : G} (hσ : σ ∈ lowerRamificationGroup B G n)
    (hτ : τ ∈ lowerRamificationGroup B G 0) :
    residue B (ramCoeff α n (τ * σ * τ⁻¹))
      = residue B (1 + ramCoeff α 0 τ) ^ n * residue B (ramCoeff α n σ) := by
  have hρ : τ * σ * τ⁻¹ ∈ lowerRamificationGroup B G n :=
    (lowerRamificationGroup_normal B G n).conj_mem σ hσ τ
  have hwv : ((isUnit_one_add_ramCoeff_zero hα hα0 hτ).unit : B) = 1 + ramCoeff α 0 τ :=
    IsUnit.unit_spec _
  have hw : α * ((isUnit_one_add_ramCoeff_zero hα hα0 hτ).unit : B) = τ • α := by
    rw [hwv]; exact mul_one_add_ramCoeff_zero hα hτ
  have hβ : maximalIdeal B = Ideal.span {τ • α} := span_smul_uniformizer hα hα0 hτ
  rw [residue_ramCoeff_change_uniformizer hα hα0 hβ hρ hw, hwv,
    ramCoeff_smul_conj hα hα0 hσ, residue_smul_of_mem_zero hτ]

/-! ## §3 Corollary 6.3

可換性を入れて `θ_0(τ)^n = 1` を出し、`G_0/G_1` の巡回性で `e_0 ∣ n` に変える。 -/

/-- **可換性を入れた段**——`G` 可換・`σ ∈ G_n \ G_{n+1}` なら
全ての `τ ∈ G_0` で `θ_0(τ)^n = 1`。

★`σ ∉ G_{n+1}` を使うのはただ 1 箇所、`residue (c_σ) ≠ 0` を出して体で約分するところ。
ここが Cor 6.3 の仮定 `G_n ≠ G_{n+1}` の使いどころの全部である。 -/
theorem residue_pow_eq_one_of_notMem {A B : Type*} [CommRing A] [CommRing B] [Algebra A B]
    [IsDomain B] [IsLocalRing B] {G : Type*} [Group G] [MulSemiringAction G B]
    [SMulCommClass G A B] {α : B} (hα : maximalIdeal B = Ideal.span {α}) (hα0 : α ≠ 0)
    (hadj : Algebra.adjoin A ({α} : Set B) = ⊤) (habel : ∀ x y : G, x * y = y * x)
    {n : ℕ} {σ : G} (hσ : σ ∈ lowerRamificationGroup B G n)
    (hσ' : σ ∉ lowerRamificationGroup B G (n + 1)) {τ : G}
    (hτ : τ ∈ lowerRamificationGroup B G 0) :
    residue B (1 + ramCoeff α 0 τ) ^ n = 1 := by
  have hc : residue B (ramCoeff α n σ) ≠ 0 := fun h =>
    hσ' ((residue_ramCoeff_eq_zero_iff hα hα0 hadj hσ).mp h)
  have hconj : τ * σ * τ⁻¹ = σ := by
    rw [habel τ σ, mul_assoc, mul_inv_cancel, mul_one]
  have h := residue_ramCoeff_conj hα hα0 hσ hτ
  rw [hconj] at h
  exact mul_right_cancel₀ hc (h.symm.trans (one_mul _).symm)

/-- 上を `θ_0 : G_0 →* 𝓀ˣ`(`thetaMul`)の言葉に直したもの。 -/
theorem thetaMul_pow_eq_one_of_notMem {A B : Type*} [CommRing A] [CommRing B] [Algebra A B]
    [IsDomain B] [IsLocalRing B] {G : Type*} [Group G] [MulSemiringAction G B]
    [SMulCommClass G A B] {α : B} (hα : maximalIdeal B = Ideal.span {α}) (hα0 : α ≠ 0)
    (hadj : Algebra.adjoin A ({α} : Set B) = ⊤) (habel : ∀ x y : G, x * y = y * x)
    {n : ℕ} {σ : G} (hσ : σ ∈ lowerRamificationGroup B G n)
    (hσ' : σ ∉ lowerRamificationGroup B G (n + 1))
    (τ : lowerRamificationGroup B G 0) :
    (thetaMul hα hα0 τ) ^ n = 1 := by
  apply Units.ext
  rw [Units.val_pow_eq_pow_val, thetaMul_apply_coe, Units.val_one]
  exact residue_pow_eq_one_of_notMem hα hα0 hadj habel hσ hσ' τ.2

/-- **`G_0/G_1` は巡回群**(`Prop 6.2` の `θ_0` の単射性の系)。
原文が `𝔽_q^×` と同一視して暗黙に使っている事実。 -/
theorem isCyclic_quot_lowerRamificationGroup_zero {A B : Type*} [CommRing A] [CommRing B]
    [Algebra A B] [IsDomain B] [IsLocalRing B] {G : Type*} [Group G] [MulSemiringAction G B]
    [SMulCommClass G A B] {α : B} (hα : maximalIdeal B = Ideal.span {α}) (hα0 : α ≠ 0)
    (hadj : Algebra.adjoin A ({α} : Set B) = ⊤)
    [Finite (lowerRamificationGroup B G 0 ⧸
      (lowerRamificationGroup B G 1).subgroupOf (lowerRamificationGroup B G 0))] :
    IsCyclic (lowerRamificationGroup B G 0 ⧸
      (lowerRamificationGroup B G 1).subgroupOf (lowerRamificationGroup B G 0)) :=
  isCyclic_of_injective_units _ (thetaMulQuot_injective hα hα0 hadj)

/-- **Yoshida Corollary 6.3(元を取り出した形)**——`G` 可換で `σ ∈ G_n \ G_{n+1}` なら
`e_0 := |G_0/G_1|` は `n` を割る。 -/
theorem nat_card_quot_zero_dvd_of_notMem {A B : Type*} [CommRing A] [CommRing B] [Algebra A B]
    [IsDomain B] [IsLocalRing B] {G : Type*} [Group G] [MulSemiringAction G B]
    [SMulCommClass G A B] {α : B} (hα : maximalIdeal B = Ideal.span {α}) (hα0 : α ≠ 0)
    (hadj : Algebra.adjoin A ({α} : Set B) = ⊤) (habel : ∀ x y : G, x * y = y * x)
    [Finite (lowerRamificationGroup B G 0 ⧸
      (lowerRamificationGroup B G 1).subgroupOf (lowerRamificationGroup B G 0))]
    {n : ℕ} {σ : G} (hσ : σ ∈ lowerRamificationGroup B G n)
    (hσ' : σ ∉ lowerRamificationGroup B G (n + 1)) :
    Nat.card (lowerRamificationGroup B G 0 ⧸
      (lowerRamificationGroup B G 1).subgroupOf (lowerRamificationGroup B G 0)) ∣ n := by
  refine nat_card_dvd_of_injective_units _ (thetaMulQuot_injective hα hα0 hadj) ?_
  intro q
  obtain ⟨τ, rfl⟩ := QuotientGroup.mk_surjective q
  apply thetaMulQuot_injective hα hα0 hadj
  rw [map_pow, map_one]
  exact thetaMul_pow_eq_one_of_notMem hα hα0 hadj habel hσ hσ' τ

/-- ★★**Yoshida Corollary 6.3**——`G` 可換かつ `G_n ≠ G_{n+1}` ⟹ `e_0 := |G_0/G_1|` が `n` を割る。

★`n = 0` のときは結論が `e_0 ∣ 0` で**空虚に真**(逸脱の記録 1)。
★可換性 `habel` を落とすと**偽**(退化の自己検査)。 -/
theorem nat_card_quot_zero_dvd_of_ne {A B : Type*} [CommRing A] [CommRing B] [Algebra A B]
    [IsDomain B] [IsLocalRing B] {G : Type*} [Group G] [MulSemiringAction G B]
    [SMulCommClass G A B] {α : B} (hα : maximalIdeal B = Ideal.span {α}) (hα0 : α ≠ 0)
    (hadj : Algebra.adjoin A ({α} : Set B) = ⊤) (habel : ∀ x y : G, x * y = y * x)
    [Finite (lowerRamificationGroup B G 0 ⧸
      (lowerRamificationGroup B G 1).subgroupOf (lowerRamificationGroup B G 0))]
    {n : ℕ}
    (hne : lowerRamificationGroup B G n ≠ lowerRamificationGroup B G (n + 1)) :
    Nat.card (lowerRamificationGroup B G 0 ⧸
      (lowerRamificationGroup B G 1).subgroupOf (lowerRamificationGroup B G 0)) ∣ n := by
  obtain ⟨σ, hσ, hσ'⟩ := SetLike.exists_of_lt
    (lt_of_le_of_ne (lowerRamificationGroup_antitone B G (Nat.le_succ n)) (Ne.symm hne))
  exact nat_card_quot_zero_dvd_of_notMem hα hα0 hadj habel hσ hσ'

/-- **完全分岐版の Corollary 6.3**(`G_0 = ⊤`、原文が実際に使う形)——
`|G/G_1| ∣ n`。原典 `Def 6.1` の直後の「`K'/K` は完全分岐だから `G = G_0`」を入れたもの。 -/
theorem nat_card_quot_top_dvd_of_ne {A B : Type*} [CommRing A] [CommRing B] [Algebra A B]
    [IsDomain B] [IsLocalRing B] {G : Type*} [Group G] [MulSemiringAction G B]
    [SMulCommClass G A B] {α : B} (hα : maximalIdeal B = Ideal.span {α}) (hα0 : α ≠ 0)
    (hadj : Algebra.adjoin A ({α} : Set B) = ⊤) (h0 : lowerRamificationGroup B G 0 = ⊤)
    (habel : ∀ x y : G, x * y = y * x)
    [Finite (G ⧸ lowerRamificationGroup B G 1)] {n : ℕ}
    (hne : lowerRamificationGroup B G n ≠ lowerRamificationGroup B G (n + 1)) :
    Nat.card (G ⧸ lowerRamificationGroup B G 1) ∣ n := by
  obtain ⟨σ, hσ, hσ'⟩ := SetLike.exists_of_lt
    (lt_of_le_of_ne (lowerRamificationGroup_antitone B G (Nat.le_succ n)) (Ne.symm hne))
  refine nat_card_dvd_of_injective_units _ (thetaMulTop_injective hα hα0 hadj h0) ?_
  intro q
  obtain ⟨τ, rfl⟩ := QuotientGroup.mk_surjective q
  apply thetaMulTop_injective hα hα0 hadj h0
  rw [map_pow, map_one]
  exact thetaMul_pow_eq_one_of_notMem hα hα0 hadj habel hσ hσ' (toLowerRamificationGroupZero h0 τ)

/-! ## §4 `G_1` は p 群 / `G_n/G_{n+1}` は指数 p

`Prop 6.6 (Sen)` が "`|⟨σ⟩| = p^m` for `m ≥ 1` (by Proposition 6.2)" で畳んでいる中身。
`θ_n : G_n/G_{n+1} ↪ 𝓀⁺`(`n ≥ 1`)の行き先が標数 `p` の体の加法群であることが全部。 -/

/-- 標数 `p` の環の加法群を乗法的に見ると指数 `p`。 -/
theorem pow_char_eq_one_multiplicative {R : Type*} [NonAssocSemiring R] (p : ℕ) [CharP R p]
    (y : Multiplicative R) : y ^ p = 1 := by
  rw [← ofAdd_toAdd y, ← ofAdd_nsmul, ofAdd_eq_one, nsmul_eq_mul, CharP.cast_eq_zero, zero_mul]

/-- **1 段上がる**——`σ ∈ G_n`(`n ≥ 1`)なら `σ^p ∈ G_{n+1}`。
`θ_n` が加法群 `𝓀⁺`(標数 `p`)への準同型で核が `G_{n+1}` であることから。 -/
theorem pow_char_mem_lowerRamificationGroup_succ {A B : Type*} [CommRing A] [CommRing B]
    [Algebra A B] [IsDomain B] [IsLocalRing B] {G : Type*} [Group G] [MulSemiringAction G B]
    [SMulCommClass G A B] (p : ℕ) [CharP (ResidueField B) p] {α : B}
    (hα : maximalIdeal B = Ideal.span {α}) (hα0 : α ≠ 0)
    (hadj : Algebra.adjoin A ({α} : Set B) = ⊤) {n : ℕ} (hn : 1 ≤ n) {σ : G}
    (hσ : σ ∈ lowerRamificationGroup B G n) :
    σ ^ p ∈ lowerRamificationGroup B G (n + 1) := by
  have h : (⟨σ, hσ⟩ : lowerRamificationGroup B G n) ^ p ∈ (thetaAdd (G := G) hα hα0 hn).ker := by
    rw [MonoidHom.mem_ker, map_pow]
    exact pow_char_eq_one_multiplicative p _
  rw [thetaAdd_ker hα hα0 hadj hn, Subgroup.mem_subgroupOf] at h
  simpa using h

/-- **`G_n/G_{n+1}` は指数 `p`**(`n ≥ 1`)。`Prop 6.2` の「初等 `p` アーベル」の内訳。 -/
theorem pow_char_eq_one_quot_lowerRamificationGroup {A B : Type*} [CommRing A] [CommRing B]
    [Algebra A B] [IsDomain B] [IsLocalRing B] {G : Type*} [Group G] [MulSemiringAction G B]
    [SMulCommClass G A B] (p : ℕ) [CharP (ResidueField B) p] {α : B}
    (hα : maximalIdeal B = Ideal.span {α}) (hα0 : α ≠ 0)
    (hadj : Algebra.adjoin A ({α} : Set B) = ⊤) {n : ℕ} (hn : 1 ≤ n)
    (x : lowerRamificationGroup B G n ⧸
      (lowerRamificationGroup B G (n + 1)).subgroupOf (lowerRamificationGroup B G n)) :
    x ^ p = 1 := by
  apply thetaAddQuot_injective hα hα0 hadj hn
  rw [map_pow, map_one]
  exact pow_char_eq_one_multiplicative p _

/-- `σ ∈ G_1` なら `σ^{p^k} ∈ G_{1+k}`(`pow_char_mem_lowerRamificationGroup_succ` の反復)。 -/
theorem pow_char_pow_mem_lowerRamificationGroup {A B : Type*} [CommRing A] [CommRing B]
    [Algebra A B] [IsDomain B] [IsLocalRing B] {G : Type*} [Group G] [MulSemiringAction G B]
    [SMulCommClass G A B] (p : ℕ) [CharP (ResidueField B) p] {α : B}
    (hα : maximalIdeal B = Ideal.span {α}) (hα0 : α ≠ 0)
    (hadj : Algebra.adjoin A ({α} : Set B) = ⊤) {σ : G}
    (hσ : σ ∈ lowerRamificationGroup B G 1) :
    ∀ k : ℕ, σ ^ p ^ k ∈ lowerRamificationGroup B G (1 + k) := by
  intro k
  induction k with
  | zero => simpa using hσ
  | succ k ih =>
      have := pow_char_mem_lowerRamificationGroup_succ p hα hα0 hadj (by omega) ih
      rw [← pow_mul, ← pow_succ] at this
      rwa [show 1 + (k + 1) = 1 + k + 1 from rfl]

/-- ★**`G_1` は `p` 群**(`p` = 剰余体の標数)。
`G_N = ⊥`(`exists_lowerRamificationGroup_eq_bot`)まで `p` 乗を繰り返すだけ。 -/
theorem isPGroup_lowerRamificationGroup_one {A B : Type*} [CommRing A] [CommRing B]
    [Algebra A B] [IsDomain B] [IsLocalRing B] [IsNoetherianRing B] {G : Type*} [Group G]
    [MulSemiringAction G B] [SMulCommClass G A B] [Finite G] [FaithfulSMul G B]
    (p : ℕ) [CharP (ResidueField B) p] {α : B}
    (hα : maximalIdeal B = Ideal.span {α}) (hα0 : α ≠ 0)
    (hadj : Algebra.adjoin A ({α} : Set B) = ⊤) :
    IsPGroup p (lowerRamificationGroup B G 1) := by
  obtain ⟨N, hN⟩ := exists_lowerRamificationGroup_eq_bot (A := A) (G := G) hadj
  intro g
  refine ⟨N, ?_⟩
  have h := pow_char_pow_mem_lowerRamificationGroup p hα hα0 hadj g.2 N
  rw [hN (1 + N) (by omega), Subgroup.mem_bot] at h
  exact Subtype.ext (by push_cast; exact h)

/-- **`σ ∈ G_1` の位数は `p` 冪**。 -/
theorem exists_orderOf_eq_pow_of_mem_lowerRamificationGroup_one {A B : Type*} [CommRing A]
    [CommRing B] [Algebra A B] [IsDomain B] [IsLocalRing B] [IsNoetherianRing B] {G : Type*}
    [Group G] [MulSemiringAction G B] [SMulCommClass G A B] [Finite G] [FaithfulSMul G B]
    (p : ℕ) [Fact p.Prime] [CharP (ResidueField B) p] {α : B}
    (hα : maximalIdeal B = Ideal.span {α}) (hα0 : α ≠ 0)
    (hadj : Algebra.adjoin A ({α} : Set B) = ⊤) {σ : G}
    (hσ : σ ∈ lowerRamificationGroup B G 1) :
    ∃ m : ℕ, orderOf σ = p ^ m := by
  obtain ⟨m, hm⟩ := (IsPGroup.iff_orderOf (p := p)).mp
    (isPGroup_lowerRamificationGroup_one (A := A) p hα hα0 hadj) ⟨σ, hσ⟩
  rw [← Subgroup.orderOf_coe] at hm
  exact ⟨m, hm⟩

/-- ★**`Prop 6.6 (Sen)` の "(by Proposition 6.2)" の中身**——`σ ∈ G_1` なら `|⟨σ⟩| = p^m`。
★`m ≥ 1`(原文)は `σ ≠ 1` と同値で、ここでは出していない(逸脱の記録 4)。 -/
theorem exists_nat_card_zpowers_eq_pow {A B : Type*} [CommRing A]
    [CommRing B] [Algebra A B] [IsDomain B] [IsLocalRing B] [IsNoetherianRing B] {G : Type*}
    [Group G] [MulSemiringAction G B] [SMulCommClass G A B] [Finite G] [FaithfulSMul G B]
    (p : ℕ) [Fact p.Prime] [CharP (ResidueField B) p] {α : B}
    (hα : maximalIdeal B = Ideal.span {α}) (hα0 : α ≠ 0)
    (hadj : Algebra.adjoin A ({α} : Set B) = ⊤) {σ : G}
    (hσ : σ ∈ lowerRamificationGroup B G 1) :
    ∃ m : ℕ, Nat.card (Subgroup.zpowers σ) = p ^ m := by
  obtain ⟨m, hm⟩ :=
    exists_orderOf_eq_pow_of_mem_lowerRamificationGroup_one (A := A) p hα hα0 hadj hσ
  exact ⟨m, by rw [Nat.card_zpowers, hm]⟩

/-! ## §5 `PAdicLocalField` への具体化

`A := 𝒪[K.carrier]`、`B := adjoinIntegers K x`、`G := Gal(K(x)/K)`。素元は `∃` の内側に閉じ込める
(取り方に依らないことは `Prop 6.2` の `§5` が保証する)。★有限性 `Finite G` は
`FiniteDimensional` から自動で出るので、`Cor 6.3` の `[Finite _]` 仮定は具体層では消える。 -/

variable {p : ℕ} [Fact p.Prime]

/-- **`𝓀_{K(x)}` の標数は `p`**。`𝓀[K.carrier] → 𝓀_{K(x)}` は体の射なので単射。 -/
theorem charP_residueField_adjoinIntegers (K : PAdicLocalField p) (x : K.closure) :
    CharP (IsLocalRing.ResidueField (adjoinIntegers K x)) p :=
  haveI := charP_residueField K
  charP_of_injective_algebraMap
    (algebraMap (IsLocalRing.ResidueField 𝒪[K.carrier])
      (IsLocalRing.ResidueField (adjoinIntegers K x))).injective p

/-- ★★**Yoshida Corollary 6.3(`PAdicLocalField` 版)**——`Gal(K(x)/K)` が可換で
完全分岐、`G_n ≠ G_{n+1}` なら `e_0 = |G/G_1|` が `n` を割る。 -/
theorem nat_card_quot_lowerRamificationGroupAdjoin_one_dvd (K : PAdicLocalField p) (x : K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (ht : IsTotallyRamifiedAdjoin K x)
    (habel : ∀ σ τ : ((IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure))),
      σ * τ = τ * σ)
    {n : ℕ}
    (hne : lowerRamificationGroupAdjoin K x n ≠ lowerRamificationGroupAdjoin K x (n + 1)) :
    Nat.card (((IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure))) ⧸
      lowerRamificationGroupAdjoin K x 1) ∣ n := by
  obtain ⟨α, hspan, hα0, hadj⟩ := exists_uniformizer_ne_zero_adjoin_eq K x ht
  exact nat_card_quot_top_dvd_of_ne hspan hα0 hadj
    (lowerRamificationGroupAdjoin_zero_eq_top K x ht) habel hne

/-- **`G_1` は `p` 群(`PAdicLocalField` 版)**。 -/
theorem isPGroup_lowerRamificationGroupAdjoin_one (K : PAdicLocalField p) (x : K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (ht : IsTotallyRamifiedAdjoin K x) :
    IsPGroup p (lowerRamificationGroupAdjoin K x 1) := by
  haveI := isDiscreteValuationRing_adjoinIntegers K x
  haveI := charP_residueField_adjoinIntegers K x
  obtain ⟨α, hspan, hα0, hadj⟩ := exists_uniformizer_ne_zero_adjoin_eq K x ht
  exact isPGroup_lowerRamificationGroup_one p hspan hα0 hadj

/-- **`G_n/G_{n+1}` は指数 `p`(`PAdicLocalField` 版、`n ≥ 1`)**。 -/
theorem pow_char_eq_one_quot_lowerRamificationGroupAdjoin (K : PAdicLocalField p) (x : K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (ht : IsTotallyRamifiedAdjoin K x) {n : ℕ} (hn : 1 ≤ n)
    (y : lowerRamificationGroupAdjoin K x n ⧸
      (lowerRamificationGroupAdjoin K x (n + 1)).subgroupOf
        (lowerRamificationGroupAdjoin K x n)) :
    y ^ p = 1 := by
  haveI := charP_residueField_adjoinIntegers K x
  obtain ⟨α, hspan, hα0, hadj⟩ := exists_uniformizer_ne_zero_adjoin_eq K x ht
  exact pow_char_eq_one_quot_lowerRamificationGroup p hspan hα0 hadj hn y

/-- ★**`Prop 6.6 (Sen)` の入力(`PAdicLocalField` 版)**——`σ ∈ G_1` なら `|⟨σ⟩| = p^m`。 -/
theorem exists_orderOf_eq_pow_of_mem_lowerRamificationGroupAdjoin_one (K : PAdicLocalField p)
    (x : K.closure)
    [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({x} : Set K.closure))]
    (ht : IsTotallyRamifiedAdjoin K x)
    {σ : (IntermediateField.adjoin K.carrier ({x} : Set K.closure))
      ≃ₐ[K.carrier] (IntermediateField.adjoin K.carrier ({x} : Set K.closure))}
    (hσ : σ ∈ lowerRamificationGroupAdjoin K x 1) :
    ∃ m : ℕ, Nat.card (Subgroup.zpowers σ) = p ^ m := by
  haveI := isDiscreteValuationRing_adjoinIntegers K x
  haveI := charP_residueField_adjoinIntegers K x
  obtain ⟨α, hspan, hα0, hadj⟩ := exists_uniformizer_ne_zero_adjoin_eq K x ht
  exact exists_nat_card_zpowers_eq_pow p hspan hα0 hadj hσ

end ABC3.Found.PGC
