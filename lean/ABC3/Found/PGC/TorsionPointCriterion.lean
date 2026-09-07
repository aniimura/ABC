import ABC3.Found.PGC.LubinTateCompletionNorm

/-!
# `α ∈ µ_{f,m} ⟺ [x](α) = 0 ⟺ [a](α) = 0 (∀a ∈ 𝔭^m)` —— Yoshida 2008 Lemma 4.3(ii)

典拠: T. Yoshida, *Local Class Field Theory via Lubin-Tate Theory*
(Ann. Fac. Sci. Toulouse 17-2, 2008; arXiv math/0606108) Lemma 4.3(物理 p.7)。
構造化済み原文は `ResearchPaper/1_Structured/Local Class Field Theory via Lubin-Tate Theory/`
の `section-4.html` の `#lemma-4-3-ii`。

原文 (Yoshida08 p.7):
> (ii) For x ∈ K^× with v(x) = m and α ∈ p[frak]_L^sep : α ∈ µ_f,m ⇐⇒ [x](α) = 0 ⇐⇒ [a](α) = 0 (∀a ∈ p[frak]^m).

原典の証明(全 3 行。★`>` を付けない ―― 逐語引用は上の 1 本だけである):

    (ii): By Definition 3.9, we have [x] = [x/π_m]_{f^{ϕ^m},f} ∘ f_m.
    As [x/π_m] is invertible, we see the first equivalence.
    The second one follows by 𝔭^m = (x).

## 何を示したか

`K` を p 進局所体、`f` を Lubin-Tate 冪級数(`coeff 0 f = 0`・`coeff 1 f = π`・
`f ≡ X^q mod 𝔭`)、`α ∈ K̄` を `‖α‖ < 1`(＝ `α ∈ 𝔭_{K̄}`)なる元、
`x ∈ 𝒪_K` を `(x) = 𝔭^m` なる元(＝ `v(x) = m`)とするとき、**3 つが同値**
(`torsionPointCriterion`、`List.TFAE`):

1. `α ∈ Λ_m`(`iteratedLubinTateTorsionPoints … m`、原典の `µ_{f,m}`)
2. `[x]_f(α) = 0`
3. `∀ a ∈ 𝔭^m, [a]_f(α) = 0`

個別の形でも出してある:
`mem_iteratedLubinTateTorsionPoints_iff_lubinTateActionAtPoint_eq_zero`(1 ⟺ 2)、
`lubinTateActionAtPoint_eq_zero_iff_forall_mem_maximalIdeal_pow`(2 ⟺ 3)、
`mem_iteratedLubinTateTorsionPoints_iff_forall_mem_maximalIdeal_pow`(1 ⟺ 3)。

## 抽象核(§0。分岐・付値・Lubin-Tate の語彙が 1 つも出てこない)

★原典の 2 行の証明は、次の**純代数**に落ちる ——
可換半環 `R` の部分集合 `S` が「吸収的」(`b ∈ S → a * b ∈ S`)であるとき:

| 宣言 | 内容 |
|---|---|
| `mem_iff_forall_mem_span_singleton_of_absorbing` | `x ∈ S ↔ ∀ a ∈ (x), a ∈ S` |
| `mem_iff_mem_of_span_singleton_eq_of_absorbing` | `(x) = (y)` ⟹ `x ∈ S ↔ y ∈ S` |

`S := {a | [a]_f(α) = 0}`(`lubinTateAnnihilator`)は吸収的である
(`lubinTateActionAtPoint_mul_eq_zero`、`[ab]_f = [a]_f ∘ [b]_f` から)。すると

* 第 2 の同値(`(x) = 𝔭^m` で `[a](α)=0 ∀a∈𝔭^m` と `[x](α)=0` が同値)は
  **`mem_iff_forall_mem_span_singleton_of_absorbing` そのもの**、
* 第 1 の同値(`[x](α)=0` と `[π^m](α)=0`、原典の「`[x/π_m]` が可逆」)は
  **`mem_iff_mem_of_span_singleton_eq_of_absorbing` そのもの**

になる。★**原典は「`[x/π_m]` の可逆性」と言うが、必要なのはそれより弱い
「`(x) = (π^m)`」(＝ `x` と `π^m` が互いに割り切る)だけである**
——`x = u·π^m` から `[x] = [u] ∘ [π^m]` と `[π^m] = [u⁻¹] ∘ [x]` の
両方向が出るので、逆写像の存在(可逆性)を言わずに済む。
★これが「原典より安い道」であり、`[x/π_m]_{f^{ϕ^m},f}` という
**2 つの形式群のあいだの準同型**(木にはまだ無い概念)を一切使わずに閉じる理由である。

## 具体層(§1–§4)

| 宣言 | 内容 |
|---|---|
| `mem_adjoinIntegers_self_of_spectralNorm_le_one` | `‖α‖ ≤ 1` なら `α ∈ 𝒪_{K(α)}` |
| `hasEval_self_of_spectralNorm_lt_one` | `‖α‖ < 1` なら `α` は位相的冪零(冪級数を代入できる) |
| ★`lubinTateActionAtPoint` | ★**`[a]_f(α)` —— `α` が捩れ点とは限らない一般の `α ∈ 𝔭_{K̄}` に対する評価** |
| ★`lubinTateActionAtPoint_eq_lubinTateActionAtTorsionPoint` | ★**`α ∈ Λ_n` のときは既存の `lubinTateActionAtTorsionPoint` と `rfl` で一致する** |
| `lubinTateAnnihilator` | `{a ∈ 𝒪_K | [a]_f(α) = 0}` |
| ★`lubinTateActionAtPoint_mul_eq_zero` | ★**吸収性 `[b](α)=0 ⟹ [ab](α)=0`**(抽象核への入力) |
| ★`mem_iteratedLubinTateTorsionPoints_iff_lubinTateActionAtPoint_pi_pow_eq_zero` | ★**`α ∈ Λ_m ⟺ [π^m]_f(α) = 0`** |

## 段取りとの差 —— 「思ったより安い」

★本体の見当では
`lubinTateActionAtTorsionPoint_eq_zero_iff_dvd_of_mem_iteratedLubinTatePsiTorsionPoints`
(`AdjoinIntegers.lean:1382`)が「(ii) の心臓」とされていたが、**それは使わなかった**。
その補題は「`x` が**原始的な** `π^n`-捩れ点(`ψ_n` の根)」を仮定して
`a·x = 0 ↔ π^n ∣ a` を言うもので、Lemma 4.3(ii) が扱う「`α` は `𝔭_{K̄}` の**任意**の元」
とは前提が逆である(あちらは点を固定して `a` を動かし、こちらは `a` を固定して点を動かす)。
★実際に効いたのは **`eq_zero_of_pi_pow_action_eq_zero`(同 `:1154`)** の方で、
これは「`adjoinIntegers K x` の**任意**の位相的冪零な元 `z`」について述べられており、
`z = α`・座標系も `α` 自身、と取るだけで `[π^m](α) = 0 ⟹ D_m(α) = 0` が出る。
★★**「作られた目的」ではなく「型」で在庫を引く**という指針がそのまま効いた例である。

## 原典との差(逸脱の記録)

1. ★**`L = K` に固定した。** 原典 §4.1 は「a complete unramified extension `L` of `K`」と
   一般に置き、`f ∈ 𝒪_L[[X]]`・`f_m = f^{ϕ^{m−1}} ∘ ⋯ ∘ f`・`π_m = ∏_{t<m} π^{ϕ^t}` とする。
   本ファイルは `L = K`(`ϕ = id`)の場合であり、そこでは
   **`π_m = π^m`・`f_m = [π^m]_f`** となる(木の `iteratedLubinTate f m` /
   `LubinTateAction_pi_pow`)。★`L^sep` は `L = K` のとき `K̄`(標数 0 なので分離閉＝代数閉)。
   §4-a・§4-b・§4-c と同じ逸脱である。
2. ★**`v(x) = m` を `Ideal.span {x} = 𝔭^m` と書いた。** 原典は `x ∈ K^×` かつ `v(x) = m` と
   書くが、`m ≥ 0` なら `x ∈ 𝒪_K` であり、離散付値環では `v(x) = m ⟺ (x) = 𝔭^m` である。
   ★この形にすると、原典が第 2 の同値の根拠として名指しする **`𝔭^m = (x)` が仮定そのもの**
   になり、証明が抽象核 1 本に落ちる。
3. ★**`α ∈ 𝔭_{L^sep}` を `spectralNorm K.carrier K.closure α < 1` と書いた。**
   木は `K̄` の整数環を `ValuationSubring` としても持つ(`absClosureInt`)が、
   `PowerSeries.aeval` を回すのに要るのは「`α` が位相的冪零」だけなので、
   ノルムの言葉で書く方が軽い(`LubinTateCompletionNorm` が「素元」をノルムで書いたのと同じ)。
4. ★**`[x/π_m]` の可逆性を使っていない**(上の「抽象核」を見よ)。原典より弱い
   `(x) = (π^m)` だけで両方向が出る。逆写像 `[π_m/x]` の構成は不要。
5. ★**`m ≥ 1` を課していない。** 原典は Lemma 4.3 全体を `m ≥ 1` の文脈で書くが、
   本ファイルの証明は `m = 0` でも通る(そのとき `Λ_0 = {0}`・`𝔭^0 = 𝒪_K`・
   `[1]_f(α) = α` で、主張は「`α = 0 ⟺ α = 0`」という真だが自明な形になる)。
   ★**仮定を落とした方向の逸脱なので、後続の消費(Prop 6.14)には影響しない。**

## 退化の自己検査

* ★★**`(x) = 𝔭^m` を落とすと偽。** `(x) ⊊ 𝔭^m`(例: `x = π^{m+1}`)なら
  `[x](α) = 0` は `α ∈ Λ_{m+1}` と同値で、`α ∈ Λ_m` より真に弱い。
  ★`𝔭^m = (x)` はまさに第 2 の同値の根拠であり、第 1 の同値の根拠でもある。
* ★★**`α ∈ 𝔭_{K̄}`(`‖α‖ < 1`)を落とすと冪級数が収束しない。**
  `lubinTateActionAtPoint` はそもそも定義できない(`PowerSeries.HasEval` が
  `hasEval_self_of_spectralNorm_lt_one` を通じて `‖α‖ < 1` を消費する)。
  ★`‖α‖ = 1` の元(単数)は `Λ_m` に入らないので、主張自体も `⟸` 向きが壊れる。
* ★**`FiniteDimensional K (K(α))` を落とすと `adjoinIntegers K α` が完備でない。**
  `K̄` は完備でないので、有限次部分体に降りて評価するのが要点である
  (`completeSpace_adjoinIntegers`)。★`α ∈ K̄` は代数的なのでこの仮定は無害だが、
  木ではインスタンスとして要求されている。
* ★**`coeff 1 f = π` を落とすと `LubinTateAction` が定義されない**(Lubin-Tate 冪級数の
  定義そのもの)。`[a]_f` の一意性は `f` の線形係数が `π` であることに依る。
* ★`ℕ∞` の切り詰め引き算・除算は 1 つも書いていない(`tools/lean-idioms.md` #102)。

## `tools/lean-idioms.md` #69(`adjoinField` と `adjoinIntegers` の境界)について

★本ファイルは `adjoinIntegers`(整数側)**だけ**を扱い、`adjoinField` を 1 度も参照しない。
`K.closure` へ降りるのは `Polynomial.hom_eval₂` を 1 回使う
`mem_iteratedLubinTateTorsionPoints_of_lubinTateActionAtPoint_pi_pow_eq_zero` のみで、
そこでも `adjoinIntegers K α →+* K.closure` という **1 本の環準同型**を明示的に組み立てている
(`AdjoinIntegers.lean` の既存の証明と同じ書き方)。★#69 の境界には当たらなかった。
-/

namespace ABC3.Found.PGC

open ABC3.Skeleton.PGC
open scoped NormedField Valued Classical

/-! ## 0. 抽象核 —— 「吸収的な部分集合」と単項イデアル

★分岐・付値・Lubin-Tate の語彙が 1 つも出てこない。原典の設定に依らない。 -/

section AbstractCore

/-- ★★★★★★★★★★★★★★★★**抽象核(1)**: 可換半環 `R` の部分集合 `S` が
**吸収的**(`b ∈ S` ならば任意の `a` について `a * b ∈ S`)であるとき、
`x ∈ S` と「単項イデアル `(x)` の元がすべて `S` に入る」は同値。

★これが原典 Lemma 4.3(ii) の**第 2 の同値**(「The second one follows
by `𝔭^m = (x)`」)の中身そのものである。`S = {a | [a]_f(α) = 0}` と
取れば、`⟸` は `x ∈ (x)`、`⟹` は `a = c·x` に対する
`[c·x](α) = [c]([x](α)) = [c](0) = 0` にあたる。 -/
theorem mem_iff_forall_mem_span_singleton_of_absorbing {R : Type*} [CommSemiring R] {S : Set R}
    (hS : ∀ a b : R, b ∈ S → a * b ∈ S) (x : R) :
    x ∈ S ↔ ∀ a ∈ Ideal.span ({x} : Set R), a ∈ S := by
  constructor
  · intro hx a ha
    obtain ⟨c, rfl⟩ := Ideal.mem_span_singleton'.mp ha
    exact hS c x hx
  · intro h
    exact h x (Ideal.mem_span_singleton_self x)

/-- ★★★★★★★★★★★★★★★★**抽象核(2)**: 吸収的な `S` について、
`(x) = (y)` ならば `x ∈ S ⟺ y ∈ S`。

★これが原典の**第 1 の同値**(「As `[x/π_m]` is invertible」)の中身である。
★★**原典は `[x/π_m]` の可逆性(逆写像の存在)を使うが、実際に要るのは
`(x) = (π_m)` だけ**である ——「互いに割り切る」から両方向が出るので、
形式群の準同型 `[x/π_m]_{f^{ϕ^m},f}` を作る必要がない。 -/
theorem mem_iff_mem_of_span_singleton_eq_of_absorbing {R : Type*} [CommSemiring R] {S : Set R}
    (hS : ∀ a b : R, b ∈ S → a * b ∈ S) {x y : R}
    (h : Ideal.span ({x} : Set R) = Ideal.span ({y} : Set R)) : x ∈ S ↔ y ∈ S := by
  rw [mem_iff_forall_mem_span_singleton_of_absorbing hS x,
    mem_iff_forall_mem_span_singleton_of_absorbing hS y, h]

end AbstractCore

/-! ## 1. `𝔭_{K̄}` の元を評価点にする

`AdjoinIntegers.lean` は評価点を `Λ_n` の元に限っていた(`Λ_n` の元は
`‖·‖ < 1` が自動だから)。Lemma 4.3(ii) は「`α ∈ 𝔭_{L^sep}` は任意」なので、
同じ 2 つの事実を `‖α‖ < 1` から直接立て直す。 -/

section Point

variable {p : ℕ} [Fact p.Prime]

/-- `‖α‖ ≤ 1` ならば `α` は `K(α)` の整数環 `adjoinIntegers K α` に入る
——`adjoinIntegers` は定義がそのまま「閉単位球」なので `show` だけ。 -/
theorem mem_adjoinIntegers_self_of_spectralNorm_le_one
    (K : PAdicLocalField p) (α : K.closure)
    (hmem : α ∈ IntermediateField.adjoin K.carrier ({α} : Set K.closure))
    (hα : spectralNorm K.carrier K.closure α ≤ 1) :
    (⟨α, hmem⟩ : IntermediateField.adjoin K.carrier ({α} : Set K.closure)) ∈
      adjoinIntegers K α := hα

/-- ★`‖α‖ < 1`(＝ `α ∈ 𝔭_{K̄}`)ならば `α` は `adjoinIntegers K α` の中で
位相的冪零、すなわち**冪級数を代入できる**
(`hasEval_mem_adjoinIntegers_of_mem_iteratedLubinTateTorsionPoints` の、
「捩れ点」を仮定しない版)。 -/
theorem hasEval_self_of_spectralNorm_lt_one
    (K : PAdicLocalField p) (α : K.closure)
    (hmem : α ∈ IntermediateField.adjoin K.carrier ({α} : Set K.closure))
    (hα : spectralNorm K.carrier K.closure α < 1) :
    PowerSeries.HasEval
      (⟨⟨α, hmem⟩, mem_adjoinIntegers_self_of_spectralNorm_le_one K α hmem hα.le⟩ :
        adjoinIntegers K α) := by
  apply tendsto_pow_atTop_nhds_zero_of_norm_lt_one
  show spectralNorm K.carrier K.closure α < 1
  exact hα

end Point

/-! ### `𝔭^m = (π^m)` -/

section Ideals

variable {p : ℕ} [Fact p.Prime]

/-- `𝔭^m = (π^m)`——`hπmax` と `Ideal.span_singleton_pow` だけ。 -/
theorem maximalIdeal_pow_eq_span_singleton_pi_pow (K : PAdicLocalField p)
    {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
    (m : ℕ) :
    IsLocalRing.maximalIdeal (𝒪[K.carrier]) ^ m =
      Ideal.span ({π ^ m} : Set (𝒪[K.carrier])) := by
  rw [hπmax, Ideal.span_singleton_pow]

end Ideals

/-! ## 2. 一般の `α ∈ 𝔭_{K̄}` に対する Lubin-Tate 作用 `[a]_f(α)` -/

section Action

variable {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
  [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])]
  {pp : ℕ} [ExpChar (IsLocalRing.ResidueField (𝒪[K.carrier])) pp]
  [Fintype (IsLocalRing.ResidueField (𝒪[K.carrier]))]
  {ff : ℕ} (hq : Fintype.card (IsLocalRing.ResidueField (𝒪[K.carrier])) = pp ^ ff)
  {π : 𝒪[K.carrier]} (hπmax : IsLocalRing.maximalIdeal (𝒪[K.carrier]) = Ideal.span {π})
  (hπne0 : π ≠ 0)
  (f : PowerSeries (𝒪[K.carrier])) (hf0 : PowerSeries.coeff 0 f = 0)
  (hf1 : PowerSeries.coeff 1 f = π)
  (hf : PowerSeries.map (IsLocalRing.residue (𝒪[K.carrier])) f = PowerSeries.X ^ (pp ^ ff))
  (α : K.closure)
  (hmem : α ∈ IntermediateField.adjoin K.carrier ({α} : Set K.closure))
  (hα : spectralNorm K.carrier K.closure α < 1)
  [FiniteDimensional K.carrier (IntermediateField.adjoin K.carrier ({α} : Set K.closure))]

/-- ★★★★★★★★★★★★★★★★**`[a]_f(α)` —— `α` が捩れ点とは限らない、
`𝔭_{K̄}` の任意の元に対する Lubin-Tate 作用**。既存の
`lubinTateActionAtTorsionPoint`(`AdjoinIntegers.lean:407`)は評価点に
`α ∈ Λ_n` を要求していたが、`PowerSeries.aeval` に本当に要るのは
`‖α‖ < 1` だけである。★`Λ_n` の元に対しては
`lubinTateActionAtPoint_eq_lubinTateActionAtTorsionPoint` により
既存の作用と `rfl` で一致するので、`AdjoinIntegers.lean` の在庫がそのまま使える。 -/
noncomputable def lubinTateActionAtPoint (a : 𝒪[K.carrier]) : adjoinIntegers K α :=
  lubinTateEvalAtPoint K α _ (hasEval_self_of_spectralNorm_lt_one K α hmem hα)
    (LubinTateAction hq hπmax f hf0 hf1 hf a)

/-- ★★★**橋**: `α ∈ Λ_n` のとき、`lubinTateActionAtPoint` は既存の
`lubinTateActionAtTorsionPoint` と**定義的に**一致する(`rfl`)。
評価点も `HasEval` の証明も `Prop` の証明部分が違うだけなので、
証明無関係で defeq になる。 -/
theorem lubinTateActionAtPoint_eq_lubinTateActionAtTorsionPoint
    (n : ℕ)
    (hx : α ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf n)
    (a : 𝒪[K.carrier]) :
    lubinTateActionAtPoint K hq hπmax f hf0 hf1 hf α hmem hα a =
      lubinTateActionAtTorsionPoint K hq hπmax hπne0 f hf0 hf1 hf n α hx hmem a := rfl

/-- `α` を消す `𝒪_K` の元全体 `{a | [a]_f(α) = 0}`。
★Lemma 4.3(ii) は「`α ∈ Λ_m` ⟺ この集合が `𝔭^m` を含む」と読める。 -/
noncomputable def lubinTateAnnihilator : Set (𝒪[K.carrier]) :=
  {a | lubinTateActionAtPoint K hq hπmax f hf0 hf1 hf α hmem hα a = 0}

omit [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])] in
theorem mem_lubinTateAnnihilator_iff (a : 𝒪[K.carrier]) :
    a ∈ lubinTateAnnihilator K hq hπmax f hf0 hf1 hf α hmem hα ↔
      lubinTateActionAtPoint K hq hπmax f hf0 hf1 hf α hmem hα a = 0 := Iff.rfl

/-! ### 吸収性 —— 抽象核への唯一の入力 -/

omit [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])] in
include hπne0 in
/-- ★★★★★★★★★★★★★★**吸収性**: `[b]_f(α) = 0` ならば `[a·b]_f(α) = 0`。
`LubinTateAction_comp`(`[ab]_f = subst([b]_f)([a]_f)`)と連鎖律
`aeval_subst_eq_aeval_aeval` で `[ab]_f(α) = [a]_f([b]_f(α)) = [a]_f(0)`
にし、`[a]_f` の定数項が `0` であること(`lubinTateEvalAtPoint_zero`)で閉じる。

★★これが**抽象核 `mem_iff_*_of_absorbing` への唯一の入力**であり、
ここから先に Lubin-Tate の語彙は出てこない。 -/
theorem lubinTateActionAtPoint_mul_eq_zero (a b : 𝒪[K.carrier])
    (hb : lubinTateActionAtPoint K hq hπmax f hf0 hf1 hf α hmem hα b = 0) :
    lubinTateActionAtPoint K hq hπmax f hf0 hf1 hf α hmem hα (a * b) = 0 := by
  haveI := completeSpace_adjoinIntegers K α
  haveI := isLinearTopology_adjoinIntegers K α
  haveI := continuousSMul_adjoinIntegers K α
  have hy : PowerSeries.HasEval (0 : adjoinIntegers K α) := IsTopologicallyNilpotent.zero
  show PowerSeries.aeval (hasEval_self_of_spectralNorm_lt_one K α hmem hα)
    (LubinTateAction hq hπmax f hf0 hf1 hf (a * b)) = 0
  rw [LubinTateAction_comp hq hπmax hπne0 f hf0 hf1 hf a b]
  rw [aeval_subst_eq_aeval_aeval (p := LubinTateAction hq hπmax f hf0 hf1 hf a)
    (show PowerSeries.HasSubst (LubinTateAction hq hπmax f hf0 hf1 hf b) by
      show IsNilpotent (PowerSeries.constantCoeff (LubinTateAction hq hπmax f hf0 hf1 hf b))
      rw [constantCoeff_LubinTateAction]; exact IsNilpotent.zero)
    (constantCoeff_LubinTateAction hq hπmax f hf0 hf1 hf b)
    (hasEval_self_of_spectralNorm_lt_one K α hmem hα) hy hb]
  exact lubinTateEvalAtPoint_zero K α hy _ (constantCoeff_LubinTateAction hq hπmax f hf0 hf1 hf a)

omit [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])] in
include hπne0 in
/-- ★吸収性を `lubinTateAnnihilator` の言葉で述べたもの(抽象核の仮定の形)。 -/
theorem absorbing_lubinTateAnnihilator :
    ∀ a b : 𝒪[K.carrier], b ∈ lubinTateAnnihilator K hq hπmax f hf0 hf1 hf α hmem hα →
      a * b ∈ lubinTateAnnihilator K hq hπmax f hf0 hf1 hf α hmem hα :=
  fun a b hb => lubinTateActionAtPoint_mul_eq_zero K hq hπmax hπne0 f hf0 hf1 hf α hmem hα a b hb

/-! ## 3. `α ∈ Λ_m ⟺ [π^m]_f(α) = 0`

★`L = K` では `f_m = [π^m]_f` である(`LubinTateAction_pi_pow`)。 -/

include hπne0 in
/-- `α ∈ Λ_m` ならば `[π^m]_f(α) = 0` —— 既存の `pi_pow_action_eq_zero` を
上の「橋」(`rfl`)を通して読み替えるだけ。 -/
theorem lubinTateActionAtPoint_pi_pow_eq_zero_of_mem (m : ℕ)
    (hx : α ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf m) :
    lubinTateActionAtPoint K hq hπmax f hf0 hf1 hf α hmem hα (π ^ m) = 0 :=
  pi_pow_action_eq_zero K hq hπmax hπne0 f hf0 hf1 hf m α hx hmem

include hπne0 in
/-- ★★★★★★★★★★`[π^m]_f(α) = 0` ならば `α ∈ Λ_m`。
`eq_zero_of_pi_pow_action_eq_zero`(`[π^m]_f = D_m · U_m` と `U_m(α)` が
単元であることから `D_m(α) = 0`)を `z = α` に適用し、`D_m(α) = 0` を
`Polynomial.hom_eval₂` で `K.closure` のレベルへ押し出して
`iteratedLubinTateTorsionPoints` の定義(`D_m` の根の `Finset`)に一致させる。

★**在庫 `eq_zero_of_pi_pow_action_eq_zero` は「`adjoinIntegers K x` の
任意の位相的冪零な元 `z`」について述べられている**ので、
`x = z = α` と取れるのがここでの要点(捩れ点であることを仮定しない)。 -/
theorem mem_iteratedLubinTateTorsionPoints_of_lubinTateActionAtPoint_pi_pow_eq_zero (m : ℕ)
    (h0 : lubinTateActionAtPoint K hq hπmax f hf0 hf1 hf α hmem hα (π ^ m) = 0) :
    α ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf m := by
  haveI := valuationRing_isDVR K
  have hD : Polynomial.aeval
      (⟨⟨α, hmem⟩, mem_adjoinIntegers_self_of_spectralNorm_le_one K α hmem hα.le⟩ :
        adjoinIntegers K α)
      (iteratedLubinTateDistinguished hq hπmax hπne0 f hf0 hf1 hf m) = 0 :=
    eq_zero_of_pi_pow_action_eq_zero K hq hπmax hπne0 f hf0 hf1 hf m α _
      (hasEval_self_of_spectralNorm_lt_one K α hmem hα) h0
  rw [iteratedLubinTateTorsionPoints, Multiset.mem_toFinset, Polynomial.mem_roots']
  refine ⟨(isDistinguishedAt_iteratedLubinTateDistinguished
    hq hπmax hπne0 f hf0 hf1 hf m).monic.map _ |>.ne_zero, ?_⟩
  set g : adjoinIntegers K α →+* K.closure :=
    (algebraMap (IntermediateField.adjoin K.carrier ({α} : Set K.closure)) K.closure).comp
      (algebraMap (adjoinIntegers K α)
        (IntermediateField.adjoin K.carrier ({α} : Set K.closure)))
    with hg_def
  have key := Polynomial.hom_eval₂ (iteratedLubinTateDistinguished hq hπmax hπne0 f hf0 hf1 hf m)
    (algebraMap (𝒪[K.carrier]) (adjoinIntegers K α)) g
    (⟨⟨α, hmem⟩, mem_adjoinIntegers_self_of_spectralNorm_le_one K α hmem hα.le⟩ :
      adjoinIntegers K α)
  rw [← Polynomial.aeval_def] at key
  have hgcomp : g.comp (algebraMap (𝒪[K.carrier]) (adjoinIntegers K α)) =
      algebraMap (𝒪[K.carrier]) K.closure := by
    apply RingHom.ext; intro y; rfl
  rw [hgcomp, hD, map_zero] at key
  show (Polynomial.map (algebraMap (𝒪[K.carrier]) K.closure)
      (iteratedLubinTateDistinguished hq hπmax hπne0 f hf0 hf1 hf m)).eval α = 0
  rw [Polynomial.eval_map]
  exact key.symm

include hπne0 in
/-- ★★★★★★★★★★★★★★**`α ∈ Λ_m ⟺ [π^m]_f(α) = 0`**
(`α ∈ 𝔭_{K̄}` は任意)。★これが Lemma 4.3(ii) の「`µ_{f,m}` は `f_m` の核」
という言い換えそのものであり、`L = K` では `f_m = [π^m]_f`。 -/
theorem mem_iteratedLubinTateTorsionPoints_iff_lubinTateActionAtPoint_pi_pow_eq_zero (m : ℕ) :
    α ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf m ↔
      lubinTateActionAtPoint K hq hπmax f hf0 hf1 hf α hmem hα (π ^ m) = 0 :=
  ⟨lubinTateActionAtPoint_pi_pow_eq_zero_of_mem K hq hπmax hπne0 f hf0 hf1 hf α hmem hα m,
    mem_iteratedLubinTateTorsionPoints_of_lubinTateActionAtPoint_pi_pow_eq_zero
      K hq hπmax hπne0 f hf0 hf1 hf α hmem hα m⟩

/-! ## 4. Lemma 4.3(ii) —— 2 つの同値 -/

include hπne0 in
/-- ★★★★★★★★★★★★★★★★★★**第 1 の同値**: `v(x) = m`(＝ `(x) = 𝔭^m`)の
とき、`α ∈ µ_{f,m} ⟺ [x]_f(α) = 0`。

原典は「`[x] = [x/π_m]_{f^{ϕ^m},f} ∘ f_m` で `[x/π_m]` は可逆」と言うが、
★**ここでは抽象核 `mem_iff_mem_of_span_singleton_eq_of_absorbing` に
`(x) = (π^m)` を渡すだけ**で済む(可逆性より弱い)。 -/
theorem mem_iteratedLubinTateTorsionPoints_iff_lubinTateActionAtPoint_eq_zero
    (m : ℕ) (x : 𝒪[K.carrier])
    (hxv : Ideal.span ({x} : Set (𝒪[K.carrier])) =
      IsLocalRing.maximalIdeal (𝒪[K.carrier]) ^ m) :
    α ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf m ↔
      lubinTateActionAtPoint K hq hπmax f hf0 hf1 hf α hmem hα x = 0 := by
  rw [mem_iteratedLubinTateTorsionPoints_iff_lubinTateActionAtPoint_pi_pow_eq_zero
    K hq hπmax hπne0 f hf0 hf1 hf α hmem hα m]
  exact (mem_iff_mem_of_span_singleton_eq_of_absorbing
    (absorbing_lubinTateAnnihilator K hq hπmax hπne0 f hf0 hf1 hf α hmem hα)
    (x := x) (y := π ^ m)
    (by rw [hxv, maximalIdeal_pow_eq_span_singleton_pi_pow K hπmax m])).symm

omit [IsAdicComplete (IsLocalRing.maximalIdeal (𝒪[K.carrier])) (𝒪[K.carrier])] in
include hπne0 in
/-- ★★★★★★★★★★★★★★★★★★**第 2 の同値**: `(x) = 𝔭^m` のとき、
`[x]_f(α) = 0 ⟺ [a]_f(α) = 0 (∀a ∈ 𝔭^m)`。

原典の「The second one follows by `𝔭^m = (x)`」の逐語訳であり、
★抽象核 `mem_iff_forall_mem_span_singleton_of_absorbing` そのものである。 -/
theorem lubinTateActionAtPoint_eq_zero_iff_forall_mem_maximalIdeal_pow
    (m : ℕ) (x : 𝒪[K.carrier])
    (hxv : Ideal.span ({x} : Set (𝒪[K.carrier])) =
      IsLocalRing.maximalIdeal (𝒪[K.carrier]) ^ m) :
    lubinTateActionAtPoint K hq hπmax f hf0 hf1 hf α hmem hα x = 0 ↔
      ∀ a ∈ IsLocalRing.maximalIdeal (𝒪[K.carrier]) ^ m,
        lubinTateActionAtPoint K hq hπmax f hf0 hf1 hf α hmem hα a = 0 := by
  rw [← hxv]
  exact mem_iff_forall_mem_span_singleton_of_absorbing
    (absorbing_lubinTateAnnihilator K hq hπmax hπne0 f hf0 hf1 hf α hmem hα) x

include hπne0 in
/-- ★`α ∈ µ_{f,m} ⟺ [a]_f(α) = 0 (∀a ∈ 𝔭^m)` —— 2 つの同値の合成。
★★**この形は `x` を含まない**ので、`v(x) = m` なる `x` を選ばずに使える。 -/
theorem mem_iteratedLubinTateTorsionPoints_iff_forall_mem_maximalIdeal_pow (m : ℕ) :
    α ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf m ↔
      ∀ a ∈ IsLocalRing.maximalIdeal (𝒪[K.carrier]) ^ m,
        lubinTateActionAtPoint K hq hπmax f hf0 hf1 hf α hmem hα a = 0 := by
  have hxv : Ideal.span ({π ^ m} : Set (𝒪[K.carrier])) =
      IsLocalRing.maximalIdeal (𝒪[K.carrier]) ^ m :=
    (maximalIdeal_pow_eq_span_singleton_pi_pow K hπmax m).symm
  exact (mem_iteratedLubinTateTorsionPoints_iff_lubinTateActionAtPoint_eq_zero
      K hq hπmax hπne0 f hf0 hf1 hf α hmem hα m (π ^ m) hxv).trans
    (lubinTateActionAtPoint_eq_zero_iff_forall_mem_maximalIdeal_pow
      K hq hπmax hπne0 f hf0 hf1 hf α hmem hα m (π ^ m) hxv)

include hπne0 in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★**Yoshida 2008 Lemma 4.3(ii)**
(`L = K` の場合):

    For `x ∈ K^×` with `v(x) = m` and `α ∈ 𝔭_{L^sep}`:
    `α ∈ µ_{f,m} ⟺ [x](α) = 0 ⟺ [a](α) = 0 (∀a ∈ 𝔭^m)`

`v(x) = m` は離散付値環では `(x) = 𝔭^m` と同値なので、そちらで書いてある
(逸脱 2)。`m = 0` でも真(逸脱 5)。 -/
theorem torsionPointCriterion (m : ℕ) (x : 𝒪[K.carrier])
    (hxv : Ideal.span ({x} : Set (𝒪[K.carrier])) =
      IsLocalRing.maximalIdeal (𝒪[K.carrier]) ^ m) :
    List.TFAE
      [α ∈ iteratedLubinTateTorsionPoints K hq hπmax hπne0 f hf0 hf1 hf m,
        lubinTateActionAtPoint K hq hπmax f hf0 hf1 hf α hmem hα x = 0,
        ∀ a ∈ IsLocalRing.maximalIdeal (𝒪[K.carrier]) ^ m,
          lubinTateActionAtPoint K hq hπmax f hf0 hf1 hf α hmem hα a = 0] := by
  tfae_have 1 ↔ 2 :=
    mem_iteratedLubinTateTorsionPoints_iff_lubinTateActionAtPoint_eq_zero
      K hq hπmax hπne0 f hf0 hf1 hf α hmem hα m x hxv
  tfae_have 2 ↔ 3 :=
    lubinTateActionAtPoint_eq_zero_iff_forall_mem_maximalIdeal_pow
      K hq hπmax hπne0 f hf0 hf1 hf α hmem hα m x hxv
  tfae_finish

def torsionPointCriterion.src : ABC3.Meta.Source :=
  { paper := "Yoshida08", pdfPage := 7, item := "Lemma 4.3", sectionId := "lemma-4-3-ii" }

end Action

end ABC3.Found.PGC
