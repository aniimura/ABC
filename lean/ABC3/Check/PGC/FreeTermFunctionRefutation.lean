import ABC3.Check.PGC.Prop12Degenerate
import ABC3.Check.PGC.Cor33Degenerate

/-!
# [pGC] §2・§3 の 3 定理は、いま `Skeleton` に在る形のままだと**偽**である

`Skeleton/PGC/Section2.lean::prop_2_2`(`:134`)・`Skeleton/PGC/Section3.lean::cor_3_1`
(`:111`)・同 `cor_3_3`(`:231`)は、いずれも原典に無い「**項** `K` の自由な関数」を
仮説に取っている:

| 定理 | 自由な項関数 | 型 |
|---|---|---|
| `prop_2_2` | `IntKbar` / `CompKbar` | `PAdicLocalField p → Type*` |
| `cor_3_1` | `isHodgeTate` | `∀ K V …, Prop` |
| `cor_3_3` | `toGal` | `∀ K, {x : K.carrier // ‖x‖ = 1} → K.absGal` |

`PAdicLocalField p` は**構造体の項**の型であって同型類の型ではない。したがって
「自由な項関数」は**同型な 2 つの項に別の値を割り当ててよい**。3 定理はどれも
「α で結ばれた 2 つの体で値(または真偽)が一致する」を結論に持つので、
値が食い違う関数を 1 つ代入すれば、そのまま反例になる。本ファイルはそれを
3 本とも `sorry` 無しで実行する(`not_prop_2_2_current_form`・
`not_cor_3_1_current_form`・`not_cor_3_3_current_form`)。

## 材料——相異なる 2 つの項とその間の同型

`Check/PGC/Prop12Degenerate.lean`(2026-09-05)が既に出している 3 点をそのまま使う:

* `twistedField p`——台の型は `ℚ_[p]` のまま、体構造だけを `x ↦ -x` に沿って移送した項。
  `twistedField_ne_selfField` により `selfField p` とは**項として異なる**。
* `twistedAlgEquiv p : TwistedQp p ≃ₐ[ℚ_[p]] ℚ_[p]`——2 体は ℚ_p-代数として**同型**。
* `twistedGalEquiv p`——したがって絶対 Galois 群も連続同型。
* `OneIsStandard`——2 つの項を区別する述語(`selfField` で真、`twistedField` で偽)。

★`Skeleton/PGC/Section3.lean:98-110` は「反例には `K ≠ K′` の witness が要るので
今のところ反証はできていない」と書いているが、この記述は 2026-09-05 に古くなっている
——`twistedField` が着地したので witness は在庫にある。本ファイルはその在庫を使う。

## ★これは原典の反証ではない

反例の 2 体は ℚ_p-代数として同型であり、原典の主張(Prop 2.2・Cor 3.1・Cor 3.3)は
この 2 体について何も食い違うことを言わない。落ちているのは原典の数学ではなく、
**我々の形式化**である:

> 自由な項関数 `IntKbar`・`isHodgeTate`・`toGal` には、同型不変性も、
> `K` の構造との関係も、一切課されていない。

すなわち `Check/PGC/Prop12Degenerate.lean` が旧 `ResidueCardinality.card` について
見つけたのと**同じ穴**である(あちらは `card_congr` を足して塞いだ)。

## 既存の反証との関係(重複していない)

* `Check/PGC/Prop22Degenerate.lean` は `prop_2_2` の**旧**形(作用が公理ゼロの `SMul`)を
  非可換性で倒す。現行形は `DistribMulAction` なのであの反例は塞がっている。
  本ファイルの反例は**作用を一切使わない**(自明な作用 `g • x = x` を入れる
  ——`trivialDistribMulAction`)ので、`DistribMulAction` への強化では塞がらない。
* `Check/PGC/Cor33Degenerate.lean` は `cor_3_3` の**旧**形(`ρ` と `ρ'` が無関係)を倒す。
  現行形は `_hρ : ∀ g, ρ' (α.equiv g) = ρ g` で塞いである。本ファイルの反例は
  その `_hρ` を**満たす** `ρ'`(`twistedRho`、`ρ` を α で運んだもの)を使い、
  代わりに `toGal` を退化させる。
* `cor_3_1` については既存の反証が無い(`Check/PGC/RefutationAttempts.lean` が
  「witness が作れない」で止まっていた)。本ファイルが初めて倒す。

## 修理の方向(本ファイルは `Skeleton` に触らない)

3 つとも、自由な項関数に**同型不変性**を課せば本ファイルの反例は消える。
`Interface/PGC/LocalFieldData.lean::ResidueCardinality.card_congr` と同じ形である:

* `prop_2_2`: 対象族に「ℚ_p-代数同型に沿った移送」を要求する
  (`no_transport_badFamily`——反例族はその移送を持てない)。
* `cor_3_1`・`cor_3_3`: 述語・辞書に `congr` 条件を要求する
  (`no_invariant_oneIsStandard`——判別述語は不変になれない)。

★どう直すかは人の判断待ち(`ResearchPaper/decisions-pending.md` の D24)。
本ファイルは判断の材料を出すだけで、`Skeleton`/`Interface` には 1 行も触っていない。

**これは原典の主張ではない**(我々のモデルと器具についての事実)ので `.src` を持たない。

## 逸脱の記録(CLAUDE.md「逸脱」)

* 仮説として書いた 3 つの `∀ …` は `Skeleton` の statement をそのまま全称量化した
  ものだが、`Skeleton` を import せずに済むよう**型を写して**いる
  (`prop_2_2` などは `sorry` なので、参照すると `sorryAx` が axiom 一覧に入る)。
  写しが原文と一致することは各定理の docstring のコードブロックで目視できる。
* 原典の `Type*`(宇宙多相)を `Type`(= `Type 0`)に固定した。`PAdicLocalField.carrier`
  自身が `Type` なので数学的な内容は変わらず、また Lean では宇宙多相な命題を
  `¬ (∀ …)` の内側で量化できない(宇宙変数を束縛できない)ため、この固定は不可避である。
  `Type 0` の場合を倒せば多相版も倒れる。`Check/PGC/Prop22Degenerate.lean` と
  `Check/PGC/Cor33Degenerate.lean` が既に同じ固定をしている。
* `cor_3_3` の反例で `E := (selfField p).carrier` と取っている(原典の `E` は自由)。
  これは `Cor33Degenerate.lean::isUniformizing_unitsHomField`(真になる側の材料)が
  `E := K.carrier` で作られているため。
-/

namespace ABC3.Check.PGC

open ABC3.Skeleton.PGC ABC3.Interface.PGC ABC3.Found.PGC
open scoped Classical

/-! ## 1. 抽象核——原典の設定に依らない部分

3 つの反例が共有しているのは「添字型の 2 点を分ける述語が 1 つでもあれば、
その添字型上の**自由な**関数を量化した主張は偽になる」という一点だけである。
分岐・付値・Galois の語彙は以下 3 宣言に一切現れない。 -/

/-- 任意のモノイドの任意の加法モノイドへの**自明な作用**(`g • x = x`)。

`DistribMulAction` の 4 公理はすべて `rfl` で従う。これは「作用を強めても
自由な型族の穴は塞がらない」ことを言うために要る——`Prop 2.2` の現行形は
`SMul` を `DistribMulAction` に強めて旧反例(`Check/PGC/Prop22Degenerate.lean`)を
塞いだが、自明な作用はどんな型族の上にも乗るので、型族が自由なままだと
何も塞がっていない。

★グローバルインスタンスにはしない(すべての `M`・`A` に自明な作用が付いてしまう)。 -/
@[implicit_reducible]
def trivialDistribMulAction (M A : Type*) [Monoid M] [AddMonoid A] : DistribMulAction M A where
  smul _ a := a
  one_smul _ := rfl
  mul_smul _ _ _ := rfl
  smul_zero _ := rfl
  smul_add _ _ _ := rfl

/-- **抽象核(型族版)**——添字 `a`・`b` を分ける述語 `P` があれば、
「任意の加法群の族 `Obj` について `Obj a ≃+ Obj b`」は偽。

`Obj i := ZMod (if P i then 2 else 3)` を代入すると `ZMod 3 ≃+ ZMod 2` を
要求することになり、濃度 3 = 2 で落ちる。★純粋に一般の添字型 `ι` の話であり、
p 進局所体も分岐も出てこない。 -/
theorem not_forall_family_addEquiv {ι : Type*} (P : ι → Prop) [DecidablePred P] {a b : ι}
    (ha : ¬ P a) (hb : P b)
    (h : ∀ (Obj : ι → Type) [∀ i, AddCommGroup (Obj i)], Nonempty (Obj a ≃+ Obj b)) :
    False := by
  obtain ⟨φ⟩ := h (fun i => ZMod (if P i then 2 else 3))
  have hc := Nat.card_congr φ.toEquiv
  rw [Nat.card_zmod, Nat.card_zmod, if_neg ha, if_pos hb] at hc
  exact absurd hc (by norm_num)

/-- **抽象核(述語版)**——添字 `a`・`b` を分ける述語 `P₀` があれば、
「任意の述語 `P` について `P a ↔ P b`」は偽。`P := P₀` を代入するだけ。 -/
theorem not_forall_pred_iff {ι : Sort*} (P₀ : ι → Prop) {a b : ι} (ha : ¬ P₀ a) (hb : P₀ b)
    (h : ∀ P : ι → Prop, P a ↔ P b) : False :=
  ha ((h P₀).mpr hb)

/-! ## 2. 具体層の共通部品

分ける述語は `OneIsStandard`(`Check/PGC/Prop12Degenerate.lean`)、
分けられる 2 点は `twistedField p`(偽)と `selfField p`(真)。
`RamificationFiltration` は `degenerateRF`(`Check/PGC/Theorem42Degenerate.lean`、
`Gv ≡ ⊤`)を使う——3 定理はどれも `RF` を結論で使っていないので、
どの `RF` を入れても反例は成り立つ。 -/

/-- 退化フィルトレーション(`Gv ≡ ⊤`)の下では、連続同型はそのまま
filtered group の同型になる——`map_Gv` は `Subgroup.map φ ⊤ = ⊤` に落ち、
`φ` の全射性だけで従う(`Check/PGC/Theorem42NaiveGC.lean::topIsoOfContinuousMulEquiv`
と同じ構成。あちらは `Found` の `filtOf`、こちらは `Skeleton/PGC/Section3` の
`filteredGroupOf` に対する版)。 -/
noncomputable def degenerateFilteredIso {p : ℕ} [Fact p.Prime] {K K' : PAdicLocalField p}
    (φ : ContinuousMulEquiv K.absGal K'.absGal) :
    FilteredGroup.Iso (filteredGroupOf (degenerateRF p) K)
      (filteredGroupOf (degenerateRF p) K') where
  equiv := φ
  map_Gv _ := Subgroup.map_top_of_surjective _ φ.surjective

/-! ## 3. Proposition 2.2 の現行形は偽 -/

/-- **★★★★★★★[pGC] Proposition 2.2 の現行形は偽**——`IntKbar`・`CompKbar` が
`PAdicLocalField p` の**項**の自由な関数だから。

原文 (pGC p.5):

> Suppose that we are given the following group-theoretic data: the topological group Γ_K,
> together with the indexed filtration Γ_K^v for all v > 0. Then the Γ_K-modules O[scr]_K[bar], and
> K[bar]∧ can be recovered group-theoretically from this group-theoretic data.

`Skeleton/PGC/Section2.lean:134` の現行形は

```
theorem prop_2_2 (_RF : RamificationFiltration p)
    (IntKbar CompKbar : PAdicLocalField p → Type*)
    [∀ K, AddCommGroup (IntKbar K)] [∀ K, DistribMulAction K.absGal (IntKbar K)]
    [∀ K, AddCommGroup (CompKbar K)] [∀ K, DistribMulAction K.absGal (CompKbar K)] :
    RecoverableAsAddModule IntKbar ∧ RecoverableAsAddModule CompKbar := sorry
```

である。`IntKbar := fun K => ZMod (if OneIsStandard p K then 2 else 3)`(抽象核の中で
代入される)と自明な作用を入れると、`twistedGalEquiv p` が
`ZMod 3 ≃+ ZMod 2` を要求するので落ちる。

★2026-09-05 の修理(`SMul` → `DistribMulAction`)はこの穴には効かない——
自明な作用はどんな型族の上にも乗る(`trivialDistribMulAction`)。 -/
theorem not_prop_2_2_current_form (p : ℕ) [Fact p.Prime] :
    ¬ (∀ (_RF : RamificationFiltration p)
        (IntKbar CompKbar : PAdicLocalField p → Type)
        [∀ K, AddCommGroup (IntKbar K)] [∀ K, DistribMulAction K.absGal (IntKbar K)]
        [∀ K, AddCommGroup (CompKbar K)] [∀ K, DistribMulAction K.absGal (CompKbar K)],
        RecoverableAsAddModule IntKbar ∧ RecoverableAsAddModule CompKbar) := by
  intro h
  refine not_forall_family_addEquiv (OneIsStandard p) (not_oneIsStandard_twisted p)
    (oneIsStandard_selfField p) ?_
  intro Obj inst
  letI : ∀ K : PAdicLocalField p, DistribMulAction K.absGal (Obj K) :=
    fun K => trivialDistribMulAction _ _
  obtain ⟨key, -⟩ := h (degenerateRF p) Obj Obj
  obtain ⟨φ, -⟩ := key (twistedGalEquiv p)
  exact ⟨φ⟩

#print axioms not_prop_2_2_current_form

/-! ## 4. Corollary 3.1 の現行形は偽 -/

/-- **★★★★★★★[pGC] Corollary 3.1 の現行形は偽**——`isHodgeTate` が
`PAdicLocalField p` の**項**の自由な述語だから。

原文 (pGC p.6):

> Given a continuous Q_p[Γ_K]-vector space of finite Q_p-dimension, the issue of whether or
> not V is Hodge-Tate (as well as the invariants d_V(i)) can be determined entirely
> group-theoretically from the filtered group Γ_K.

`Skeleton/PGC/Section3.lean:111` の現行形は

```
theorem cor_3_1 (RF : RamificationFiltration p)
    (isHodgeTate : ∀ (K : PAdicLocalField p) (V : Type)
      [AddCommGroup V] [Module ℚ_[p] V] [SMul K.absGal V], Prop) :
    ∀ {K K' : PAdicLocalField p}
      (α : FilteredGroup.Iso (filteredGroupOf RF K) (filteredGroupOf RF K'))
      (V : Type) [AddCommGroup V] [Module ℚ_[p] V] [FiniteDimensional ℚ_[p] V]
      (sK : SMul K.absGal V) (sK' : SMul K'.absGal V)
      (_hcompat : ∀ (g : K.absGal) (x : V), sK'.smul (α.equiv g) x = sK.smul g x),
      @isHodgeTate K V _ _ sK ↔ @isHodgeTate K' V _ _ sK' := sorry
```

である。`isHodgeTate K V := OneIsStandard p K`(`V` を無視する)と取り、
`V := ℚ_[p]`・両側とも自明な `SMul`(`_hcompat` は `rfl`)を入れると、結論は
`OneIsStandard p (twistedField p) ↔ OneIsStandard p (selfField p)`——すなわち
`False ↔ True` になる。

★2026-09-05 の修理(`isHodgeTate` を `V` にも依存させ、`_hcompat` で 2 つの作用を
結んだ)はこの穴には効かない——`V` を無視する述語を代入すればよいだけだから。 -/
theorem not_cor_3_1_current_form (p : ℕ) [Fact p.Prime] :
    ¬ (∀ (RF : RamificationFiltration p)
        (isHodgeTate : ∀ (K : PAdicLocalField p) (V : Type)
          [AddCommGroup V] [Module ℚ_[p] V] [SMul K.absGal V], Prop),
        ∀ {K K' : PAdicLocalField p}
          (α : FilteredGroup.Iso (filteredGroupOf RF K) (filteredGroupOf RF K'))
          (V : Type) [AddCommGroup V] [Module ℚ_[p] V] [FiniteDimensional ℚ_[p] V]
          (sK : SMul K.absGal V) (sK' : SMul K'.absGal V)
          (_hcompat : ∀ (g : K.absGal) (x : V), sK'.smul (α.equiv g) x = sK.smul g x),
          @isHodgeTate K V _ _ sK ↔ @isHodgeTate K' V _ _ sK') := by
  intro h
  refine not_forall_pred_iff (OneIsStandard p) (not_oneIsStandard_twisted p)
    (oneIsStandard_selfField p) ?_
  intro P
  exact h (degenerateRF p) (fun K _ _ _ _ => P K)
    (degenerateFilteredIso (twistedGalEquiv p)) ℚ_[p]
    ⟨fun _ x => x⟩ ⟨fun _ x => x⟩ (fun _ _ => rfl)

#print axioms not_cor_3_1_current_form

/-! ## 5. Corollary 3.3 の現行形は偽

`cor_3_3` は 2026-09-05 に `_hρ`(`ρ'` を α で `ρ` と結ぶ条件)を足して
`Check/PGC/Cor33Degenerate.lean` の反例を塞いだ。ところが `toGal`——
局所類体論の相互律 `U_K → Γ_K` の代役として置かれた辞書——は**項 `K` の
自由な関数のまま**である。そこで `ρ` の側は `_hρ` を満たすように作り、
`toGal` の側だけを退化させる。 -/

/-- `ρ ∘ toGal` が恒等的に `1` なら uniformizing でない。

`Check/PGC/Def32Degenerate.lean::not_isUniformizing_one`(`ρ = 1` の場合)の一般化で、
証明も同一——`I` は開で `1 ∈ I` だから `1` の近傍を含み、体準同型 `ι` は単射なので
`ι x = 1 = ι 1` は `x = 1` を強いるが、`K` の位相は非離散なのでそんな近傍は無い。
★`Def32Degenerate.lean` を書き換えないため、一般化はこちらに置く。 -/
theorem not_isUniformizing_of_comp_const {p : ℕ} [Fact p.Prime] (K : PAdicLocalField p)
    (E : Type*) [Field E] [Algebra ℚ_[p] E]
    (toGal : {x : K.carrier // ‖x‖ = (1 : ℝ)} → K.absGal)
    (ρ : K.absGal →* Eˣ) (hconst : ∀ x, ρ (toGal x) = 1) :
    ¬ IsUniformizing K E toGal ρ := by
  rintro ⟨I, hIU, hopen, hone, -, -, ι, hι⟩
  obtain ⟨ε, hε, hball⟩ := Metric.isOpen_iff.mp hopen 1 hone
  obtain ⟨δ, hδpos, hδlt⟩ := NormedField.exists_norm_lt K.carrier hε
  have hmem : (1 + δ) ∈ I := by
    refine hball ?_
    simp only [Metric.mem_ball, dist_eq_norm]
    simpa using hδlt
  have h1 := hι (1 + δ) hmem
  have h2 := hι 1 hone
  rw [hconst, Units.val_one] at h1 h2
  have h3 : ι (1 + δ) = ι 1 := by rw [← h1, ← h2]
  have h4 : (1 : K.carrier) + δ = 1 := ι.injective h3
  have h5 : δ = 0 := by linear_combination h4
  rw [h5, norm_zero] at hδpos
  exact lt_irrefl 0 hδpos

/-- 病的な `toGal`——標準側では `Check/PGC/Cor33Degenerate.lean::toGalChoice`
(`Γ_K ↠ 𝒪_K^×` の切断、`IsUniformizing` を**真**にする本物)、
捻り側では定数 `1`(`IsUniformizing` を**偽**にする)。

★2 体は `twistedAlgEquiv p` で ℚ_p-代数同型なのに、`toGal` にはそれを咎める条件が
1 つも無い——これが `cor_3_3` の穴である。 -/
noncomputable def badToGal (p : ℕ) [Fact p.Prime] (K : PAdicLocalField p) :
    {x : K.carrier // ‖x‖ = (1 : ℝ)} → K.absGal :=
  if OneIsStandard p K then toGalChoice K else fun _ => 1

theorem badToGal_selfField (p : ℕ) [Fact p.Prime] :
    badToGal p (selfField p) = toGalChoice (selfField p) :=
  if_pos (oneIsStandard_selfField p)

theorem badToGal_twisted (p : ℕ) [Fact p.Prime] :
    badToGal p (twistedField p) = fun _ => 1 :=
  if_neg (not_oneIsStandard_twisted p)

/-- 捻り側の表現——標準側の `unitsHomField` を `twistedGalEquiv` で運んだもの。
これが現行形の修理条件 `_hρ` を満たす(`twistedRho_symm`)。 -/
noncomputable def twistedRho (p : ℕ) [Fact p.Prime] :
    (twistedField p).absGal →* ((selfField p).carrier)ˣ :=
  (unitsHomField (selfField p)).comp (twistedGalEquiv p).toMulEquiv.toMonoidHom

/-- `_hρ` の中身——`α := (twistedGalEquiv p).symm` について
`ρ' (α g) = ρ g` が成り立つ。すなわち反例は 2026-09-05 の修理を**満たしている**。 -/
theorem twistedRho_symm (p : ℕ) [Fact p.Prime] (g : (selfField p).absGal) :
    twistedRho p ((twistedGalEquiv p).symm g) = unitsHomField (selfField p) g := by
  show unitsHomField (selfField p) ((twistedGalEquiv p) ((twistedGalEquiv p).symm g)) = _
  rw [ContinuousMulEquiv.apply_symm_apply]

/-- **★★★★★★★[pGC] Corollary 3.3 の現行形は偽**——`toGal` が
`PAdicLocalField p` の**項**の自由な関数だから。

原文 (pGC p.6):

> Given a continuous E[Γ_K]-module V of E-dimension 1, the issue of whether or not V is
> uniformizing can be determined entirely group-theoretically from the filtered group Γ_K.

`Skeleton/PGC/Section3.lean:231` の現行形は

```
theorem cor_3_3 (RF : RamificationFiltration p)
    (E : Type*) [Field E] [Algebra ℚ_[p] E]
    (toGal : ∀ K : PAdicLocalField p, {x : K.carrier // ‖x‖ = (1 : ℝ)} → K.absGal) :
    ∀ {K K' : PAdicLocalField p}
      (α : FilteredGroup.Iso (filteredGroupOf RF K) (filteredGroupOf RF K'))
      (ρ : K.absGal →* Eˣ) (ρ' : K'.absGal →* Eˣ)
      (_hρ : ∀ g : K.absGal, ρ' (α.equiv g) = ρ g),
      IsUniformizing K E (toGal K) ρ ↔ IsUniformizing K' E (toGal K') ρ' := sorry
```

である。`K := selfField p`・`K' := twistedField p`・`toGal := badToGal p` と取ると

* 左辺は真(`Cor33Degenerate.lean::isUniformizing_unitsHomField`)、
* 右辺は偽(`badToGal p (twistedField p) = fun _ => 1` なので `ρ' ∘ toGal ≡ 1`)

となって `↔` が破れる。★`_hρ` は満たしている(`twistedRho_symm`)ので、
2026-09-05 の修理では塞がらない。 -/
theorem not_cor_3_3_current_form (p : ℕ) [Fact p.Prime] :
    ¬ (∀ (RF : RamificationFiltration p) (E : Type) [Field E] [Algebra ℚ_[p] E]
        (toGal : ∀ K : PAdicLocalField p, {x : K.carrier // ‖x‖ = (1 : ℝ)} → K.absGal),
        ∀ {K K' : PAdicLocalField p}
          (α : FilteredGroup.Iso (filteredGroupOf RF K) (filteredGroupOf RF K'))
          (ρ : K.absGal →* Eˣ) (ρ' : K'.absGal →* Eˣ)
          (_hρ : ∀ g : K.absGal, ρ' (α.equiv g) = ρ g),
          IsUniformizing K E (toGal K) ρ ↔ IsUniformizing K' E (toGal K') ρ') := by
  intro h
  have key := (h (degenerateRF p) ((selfField p).carrier) (badToGal p)
      (degenerateFilteredIso (twistedGalEquiv p).symm)
      (unitsHomField (selfField p)) (twistedRho p) (fun g => twistedRho_symm p g)).mp
  refine not_isUniformizing_of_comp_const (twistedField p) ((selfField p).carrier)
      (badToGal p (twistedField p)) (twistedRho p) ?_ (key ?_)
  · intro x
    simp only [badToGal_twisted, map_one]
  · rw [badToGal_selfField]
    exact isUniformizing_unitsHomField (selfField p)

#print axioms not_cor_3_3_current_form

/-! ## 6. 修理の方向——同型不変性を課すと反例族は存在しなくなる

`Interface/PGC/LocalFieldData.lean::ResidueCardinality.card_congr`(2026-09-05)と
同じ形の条件を自由な項関数に課せばよい、ということを 2 つの補題で示す。
★どう直すかの判断そのものは人に委ねる(D24)。 -/

/-- 反例に使った型族の明示形(抽象核 `not_forall_family_addEquiv` が内部で作るものと同じ)。 -/
noncomputable def badFamily (p : ℕ) [Fact p.Prime] (K : PAdicLocalField p) : Type :=
  ZMod (if OneIsStandard p K then 2 else 3)

noncomputable instance instAddCommGroupBadFamily (p : ℕ) [Fact p.Prime]
    (K : PAdicLocalField p) : AddCommGroup (badFamily p K) :=
  inferInstanceAs (AddCommGroup (ZMod _))

theorem card_badFamily (p : ℕ) [Fact p.Prime] (K : PAdicLocalField p) :
    Nat.card (badFamily p K) = if OneIsStandard p K then 2 else 3 :=
  Nat.card_zmod _

/-- **★★★修理の確認(型族)**——対象族に「ℚ_p-代数同型に沿った移送」を要求すれば、
反例族 `badFamily` はもう項として許されない。 -/
theorem no_transport_badFamily (p : ℕ) [Fact p.Prime] :
    ¬ ∀ {K K' : PAdicLocalField p}, (K.carrier ≃ₐ[ℚ_[p]] K'.carrier) →
        Nonempty (badFamily p K ≃+ badFamily p K') := by
  intro h
  obtain ⟨φ⟩ := h (K := twistedField p) (K' := selfField p) (twistedAlgEquiv p)
  have hc := Nat.card_congr φ.toEquiv
  rw [card_badFamily, card_badFamily, if_neg (not_oneIsStandard_twisted p),
    if_pos (oneIsStandard_selfField p)] at hc
  exact absurd hc (by norm_num)

/-- **★★★修理の確認(述語・辞書)**——判別述語 `OneIsStandard` は
ℚ_p-代数同型で不変になれない。したがって `isHodgeTate`・`toGal` に
「同型で移り合う」条件を課せば、本ファイルの反例は 3 本とも消える。 -/
theorem no_invariant_oneIsStandard (p : ℕ) [Fact p.Prime] :
    ¬ ∀ {K K' : PAdicLocalField p}, (K.carrier ≃ₐ[ℚ_[p]] K'.carrier) →
        (OneIsStandard p K ↔ OneIsStandard p K') :=
  fun h => not_oneIsStandard_twisted p
    ((h (K := twistedField p) (K' := selfField p) (twistedAlgEquiv p)).mpr
      (oneIsStandard_selfField p))

end ABC3.Check.PGC
