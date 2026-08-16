import ABC3.Found.FrdI.Prop111

/-!
# [FrdI] Proposition 1.14 —— Irreducible Morphisms

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、
物理 p.41(**目視確認 2026-08-16**、親が p.41 と p.17–18 を描画して照合)。

原文 (FrdI p.41):
> on a connected, totally epimorphic category D; C →FΦ a Frobenioid of isotropic

原文 (FrdI p.41):
> type; φ ∈Arr(C). Suppose further that D is of FSMFF-type [cf. §0]. Then:

## ★この命題の規模(着手前の測定)

**5 条 (i)–(v)**。★検証役(私の文脈を持たない子)が原文を読んで数えたところ
**独立な主張は 12〜13**、★**原文が省いた段は 11 箇所**であった。
内訳は「may be assumed」型が 2、「follows immediately/formally」型が 5、
数え上げの圧縮が 2、`≼` への移行が 1、そして
★**引用符で名前だけ与える型が 1**((ii) の「adding the pull-backs of …」)。

## ★不足していた語彙(2026-08-16 に補った)

★検証役の調査で、当初把握していた 4 件のうち **2 件は既に実装済み**だった:
`mid-adjoint`(`IsMidAdjoint`)と `non-dilating`(`IsNonDilating`)。
★**足りなかったのは `FSMI` / `FSMFF-type`、そして把握していなかった
`subordinate`** である(`subordinate` は主張に現れず**証明にだけ**現れる)。
いずれも `CategoryVocabulary.lean` に置いた。

## ★依存

★検証役の機械集計と通読の両方で、★**`Proposition 1.6` には依存しない**。
必要なのは `Prop 1.10` の (i)(ii)(iv)(v) と `Prop 1.11` の (v)(vi)(vii) である。

## ★★(i) の構造

原文 (FrdI p.41):
> (i) φ is irreducible if and only if φ is one of the following: (a) a prime-

原文 (FrdI p.41):
> Frobenius morphism; (b) a step such that Div(φ) is irreducible; (c) a pull-back

原文 (FrdI p.41):
> morphism such that Base(φ) is an irreducible morphism of D.

★**原文の証明**(p.41、目視):
> The sufficiency of the condition of assertion (i) follows for morphisms as in (a)
> (respectively, (b); (c)) from Proposition 1.10, (iv) (respectively, the equivalences
> of categories of Definition 1.3, (iii), (d) [cf. also Propositions 1.4, (i); 1.7, (v)];
> Proposition 1.11, (vi)). To verify the necessity of the condition of assertion (i),
> observe that it follows formally from the factorization of Definition 1.3, (iv), (a),
> that φ is either a morphism of Frobenius type, a step, or a pull-back morphism.

★★**(a) と (c) は既存の定理の引用で終わる** —— `prop_1_10_iv_mp` と
`prop_1_11_vi_irred`。★**内容があるのは (b) だけ**である。

★(b) の両向きが、第1の圏同値(コスライス)の**2 つの使い方**に対応する:
- **十分性**: `Div(φ)` の分解 `Φ(Base β)(Div α) + Div β` に irreducible を当てる ——
  ★**圏同値は要らない**。`Remark 1.1.1` と isotropic 性だけで出る
- **必要性**: `Div(φ) = b + c` から `φ` の分解を作る ——
  ★**ここで圏同値が要る**(`coaPre_realize` ＋ `coaPre_factor_under_of_mle`)

★★**原文は両向きに同じ根拠(「the equivalences of categories」)を挙げているが、
実際に使うのは片方だけである** —— 「原文が挙げる根拠が実際の根拠の一部」の
また 1 例。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ)

/-! ## ★(i) の十分性 —— 3 つの場合 -/

include P in
/-- **(i)(a) の十分性** —— prime-Frobenius 射は irreducible。

★`Proposition 1.10, (iv)` そのもの。★**原文が `𝒞` を isotropic 型に取っているので、
`Prop 1.10, (iv)` が要求する「域が isotropic」は自動で満たされる。** -/
theorem prop_1_14_i_of_primeFrob (F : FrobenioidCore P) (hiso : ∀ X : C, IsIsotropic P X)
    {A B : C} (φ : A ⟶ B) (h : IsPrimeFrobenius P φ) : IsIrreducibleMor φ :=
  prop_1_10_iv_mp P F (hiso A) φ h

include P in
/-- **(i)(c) の十分性** —— 底が irreducible な pull-back 射は irreducible。

★`Proposition 1.11, (vi)` そのもの。 -/
theorem prop_1_14_i_of_pullBack (F : FrobenioidCore P) {A B : C} (φ : A ⟶ B)
    (hpb : IsPullBack P φ) (hb : IsIrreducibleMor (P.Base φ)) : IsIrreducibleMor φ :=
  prop_1_11_vi_irred_of_baseIrred P F φ hpb hb

include P in
/-- ★★**(i)(b) の十分性** —— `Div(φ)` が irreducible な step は irreducible。

★★**原文は「the equivalences of categories of Definition 1.3, (iii), (d)」を
根拠に挙げるが、この向きには圏同値は要らない。**
要るのは 2 つだけである:
1. `Proposition 1.7, (v)`(pre-step の分解は両因子とも pre-step)
2. `Remark 1.1.1`(`Div` の合成則)と、★**`𝒞` が isotropic 型であること**

★**isotropic 性の使いどころ**: 「`Div = 0` の pre-step は同型」。
`Div(φ) = Φ.map(Base β)(Div α) + Div β` に irreducible を当てると
どちらかの項が `0` になり、★**そこから片方の因子が isometric pre-step になって
同型に落ちる。**

★`Φ.map` の単射性(`Definition 1.1, (ii), (a)` の characteristic injectivity)が
`Φ.map(Base β)(Div α) = 0` から `Div α = 0` を出す —— ★**base が同型であることは
使っていない**(単射性はつねに成り立つ)。 -/
theorem prop_1_14_i_of_step (hiso : ∀ X : C, IsIsotropic P X)
    {A B : C} (φ : A ⟶ B) (hst : IsStep P φ) (hdiv : IsIrreducibleElt (P.Div φ)) :
    IsIrreducibleMor φ := by
  refine ⟨hst.2, fun X β α hfac => ?_⟩
  obtain ⟨hβs, hαs⟩ := prop_1_7_v_preStep P β α (hfac ▸ hst.1)
  have hlin : P.degFr α = 1 := hαs.1
  have hd : P.Div φ = Φ.map (P.Base β) (P.Div α) + P.Div β := by
    rw [hfac, P.Div_comp, hlin]
    simp
  rcases hdiv.2 _ _ hd with h1 | h2
  · right
    refine hiso X B α ?_ hαs
    refine Φ.map_injective (P.Base β) ?_
    rw [h1, map_zero]
  · exact Or.inl (hiso A X β h2 hβs)

/-! ## ★★(i) の必要性 —— `Div(φ)` の irreducibility

★**ここが (i) の内容のある部分**である。

★**組み立て**: `Div(φ) = b + c` が与えられたとき、
1. `coaPre_realize`(第1の圏同値の**本質的全射性**)が `Div ψ = b` なる
   co-angular pre-step `ψ : A ⟶ X` を作る
2. `coaPre_factor_under_of_mle`(第1の圏同値の**充満性**)が
   `ψ ≫ δ = φ` なる `δ` を作る
3. `φ` の irreducibility が `ψ` か `δ` を同型にする
4. 同型なら `Div = 0`(`Proposition 1.7, (i)`)なので、`b = 0` または `c = 0`

★★**圏同値の 2 つの性質(本質的全射性と充満性)が、それぞれ別の段で効いている。**
★原文の「the equivalences of categories」の一語には、この 2 段が入っている。
-/

include P in
/-- ★★★**(i)(b) の必要性** —— irreducible な step の `Div` は irreducible。 -/
theorem prop_1_14_i_div_irreducible (G : Frobenioid P) (hiso : ∀ X : C, IsIsotropic P X)
    {A B : C} (φ : A ⟶ B) (hst : IsStep P φ) (hirr : IsIrreducibleMor φ) :
    IsIrreducibleElt (P.Div φ) := by
  constructor
  · -- `Div φ = 0` なら `φ` は isometric pre-step、よって同型 —— step に反する
    intro h0
    exact hst.2 (hiso A B φ h0 hst.1)
  · intro b c hbc
    have hφc : IsCoAngular P φ := prop_1_4_i P φ (fun Y _ => hiso Y)
    obtain ⟨X, ψ, hψc, hψs, hψd⟩ := coaPre_realize P G A b
    have hle : MLe (P.Div ψ) (P.Div φ) := ⟨c, by rw [hψd, hbc]⟩
    obtain ⟨δ, _, hδs, hδ⟩ := coaPre_factor_under_of_mle P G ψ hψc hψs φ hφc hst.1 hle
    rcases hirr.2 X ψ δ hδ.symm with hiψ | hiδ
    · haveI := hiψ
      exact Or.inl (hψd ▸ isIsometric_of_isIso P ψ)
    · haveI := hiδ
      refine Or.inr ?_
      -- `δ` が同型なら `Div φ = Div ψ = b`、よって `b + c = b + 0`
      have hdd : P.Div φ = b := by
        rw [← hδ, P.Div_comp, isIsometric_of_isIso P δ,
          show P.degFr δ = 1 from hδs.1, hψd]
        simp
      letI := isCancelAdd_of_isIntegralMonoid _ (P.divisorial (P.toElem.obj A).base).1.1
      have : b + c = b + 0 := by rw [add_zero, ← hbc, hdd]
      exact add_left_cancel this

/-! ## ★★(i) の必要性 —— 3 分類そのもの

★**原文**(p.41、目視):
> observe that it follows formally from the factorization of Definition 1.3, (iv), (a),
> that φ is either a morphism of Frobenius type, a step, or a pull-back morphism.

★★**「follows formally」の中身は irreducibility の 2 回の適用である。**
`Definition 1.3, (iv), (a)` は `φ = γ ≫ β ≫ α`(Frobenius 型・pre-step・pull-back)を与える。
1. `φ = γ ≫ (β ≫ α)` に irreducibility ⟹ `γ` が同型 または `β ≫ α` が同型
2. 後者なら `φ` は Frobenius 型。前者なら `φ = (γ ≫ β) ≫ α` にもう一度当てて、
   `α` が同型(⟹ `φ` は pre-step)か `γ ≫ β` が同型(⟹ `φ` は pull-back)

★★**3 分類が「2 回の二分」から出る** —— 原文の「either … or … or」の 3 つは
対等ではなく、★**分解の順序が決めている。**
-/

include P in
/-- ★★**(i) の必要性の第1段** —— irreducible な射は
「Frobenius 型」「pre-step」「pull-back」のいずれかである。

★`Definition 1.3, (iv), (a)` の分解に irreducibility を 2 回当てるだけ。 -/
theorem prop_1_14_i_trichotomy (F : FrobenioidCore P) {A B : C} (φ : A ⟶ B)
    (hirr : IsIrreducibleMor φ) :
    IsFrobeniusType P φ ∨ IsPreStep P φ ∨ IsPullBack P φ := by
  obtain ⟨X, Y, γ, β, α, hfac, hγF, hβs, hαpb⟩ := F.arbFactor φ
  rcases hirr.2 X γ (β ≫ α) hfac with hγi | hri
  · -- `γ` が同型: もう一度当てる
    haveI := hγi
    rcases hirr.2 Y (γ ≫ β) α (by rw [hfac, Category.assoc]) with hgbi | hαi
    · -- `γ ≫ β` が同型 ⟹ `φ` は pull-back
      haveI := hgbi
      exact Or.inr (Or.inr (hfac ▸ (Category.assoc γ β α ▸
        IsPullBack.comp P (isPullBack_of_isIso P (γ ≫ β)) hαpb)))
    · -- `α` が同型 ⟹ `φ` は pre-step
      haveI := hαi
      refine Or.inr (Or.inl (hfac ▸ ?_))
      exact IsPreStep.comp P (isPreStep_of_isIso P γ)
        (IsPreStep.comp P hβs (isPreStep_of_isIso P α))
  · -- `β ≫ α` が同型 ⟹ `φ` は Frobenius 型
    haveI := hri
    exact Or.inl (hfac ▸ IsFrobeniusType.comp P F hγF (isFrobeniusType_of_isIso P (β ≫ α)))

include P in
/-- ★★★**[FrdI] Proposition 1.14, (i)** —— irreducible の 3 分類。

★**必要性の最後の一歩**は 3 つとも既存の定理である:
- Frobenius 型 ＋ irreducible ⟹ prime-Frobenius (`Proposition 1.10, (iv)`)
- pre-step ＋ ¬同型 ＝ step、その `Div` が irreducible (上の `prop_1_14_i_div_irreducible`)
- pull-back ＋ irreducible ⟹ 底も irreducible (`Proposition 1.11, (vi)`) -/
theorem prop_1_14_i (G : Frobenioid P) (hiso : ∀ X : C, IsIsotropic P X)
    {A B : C} (φ : A ⟶ B) :
    IsIrreducibleMor φ ↔
      (IsPrimeFrobenius P φ ∨
        (IsStep P φ ∧ IsIrreducibleElt (P.Div φ)) ∨
        (IsPullBack P φ ∧ IsIrreducibleMor (P.Base φ))) := by
  constructor
  · intro hirr
    rcases prop_1_14_i_trichotomy P G.core φ hirr with hF | hps | hpb
    · exact Or.inl (prop_1_10_iv_mpr P G.core (fun Y _ => hiso Y) φ hF hirr)
    · exact Or.inr (Or.inl ⟨⟨hps, hirr.1⟩,
        prop_1_14_i_div_irreducible P G hiso φ ⟨hps, hirr.1⟩ hirr⟩)
    · exact Or.inr (Or.inr ⟨hpb,
        prop_1_11_vi_baseIrred_of_irred P G.core φ hpb hirr⟩)
  · rintro (h | ⟨hst, hdiv⟩ | ⟨hpb, hb⟩)
    · exact prop_1_14_i_of_primeFrob P G.core hiso φ h
    · exact prop_1_14_i_of_step P hiso φ hst hdiv
    · exact prop_1_14_i_of_pullBack P G.core φ hpb hb

/-! ## ★(ii) —— pre-step ⟺ FSM ＋ mid-adjoint

原文 (FrdI p.41):
> (ii) φ is a pre-step if and only if it is an FSM-morphism that is mid-

原文 (FrdI p.41):
> adjoint [cf. §0] to the irreducible morphisms which are not pre-steps.

★**原文の証明の必要性**(p.42、目視):
> Thus, suppose that φ is a pre-step. By Proposition 1.11, (vii), φ is an FSM-morphism;
> by Proposition 1.7, (v), φ is mid-adjoint to the non-pre-steps.

★★**原文は主張より強いものを証明している** —— 主張は「**irreducible な**
非 pre-step に mid-adjoint」だが、証明は「**非 pre-step 全体**に mid-adjoint」を出す。
★`S` が大きいほど mid-adjoint は強い主張なので、これは余剰である。
★★**`Proposition 1.11, (vii)` でも同じことが起きた**(原文の証明が主張より強い `α` を作る)。
★**この命題群は「証明が主張より多くを語る」箇所を繰り返し持っている。**

★**必要性の中身は空虚性である** —— `Proposition 1.7, (v)` により
pre-step の 3 分解の真ん中は必ず pre-step だから、
★**「真ん中が非 pre-step」という仮定を満たす分解が存在しない。**
-/

/-- ★(ii) の `S` —— 「pre-step でない irreducible 射」。 -/
def irredNonPreStep : MorphismProperty C := fun _ _ f => IsIrreducibleMor f ∧ ¬ IsPreStep P f

/-- ★原文の証明が実際に使う `S` —— 「pre-step でない射」全体(★こちらが**大きい**)。 -/
def nonPreStep : MorphismProperty C := fun _ _ f => ¬ IsPreStep P f

include P in
/-- ★**大きい `S` に mid-adjoint なら、小さい `S` にも mid-adjoint**。 -/
theorem isMidAdjoint_irredNonPreStep_of_nonPreStep {A B : C} (φ : A ⟶ B)
    (h : IsMidAdjoint (nonPreStep P) φ) : IsMidAdjoint (irredNonPreStep P) φ :=
  fun X Y γ β α hfac hS => h X Y γ β α hfac hS.2

include P in
/-- ★★**(ii) の必要性(しかも原文の証明どおり、強い形で)** ——
pre-step は FSM 射であり、★**非 pre-step 全体に mid-adjoint** である。

★**FSM 性は `Proposition 1.11, (vii)`**(isotropic 型なので全射が co-angular)、
★**mid-adjoint 性は `Proposition 1.7, (v)` による空虚性**。 -/
theorem prop_1_14_ii_mp (F : FrobenioidCore P) (G : Frobenioid P)
    (hiso : ∀ X : C, IsIsotropic P X) {A B : C} (φ : A ⟶ B) (hφ : IsPreStep P φ) :
    IsFSMMorphism φ ∧ IsMidAdjoint (nonPreStep P) φ := by
  refine ⟨prop_1_11_vii_fsm_of_coaPre P F G hiso φ
    (prop_1_4_i P φ (fun Y _ => hiso Y)) hφ, ?_⟩
  intro X Y γ β α hfac hS
  exact absurd
    (prop_1_7_v_preStep P β α (prop_1_7_v_preStep P γ (β ≫ α) (hfac ▸ hφ)).2).1 hS

/-! ### ★(ii) の十分性 —— 未実装(2026-08-16 の測定)

★**原文の証明**(p.42、目視、要旨):
1. `Definition 1.3, (iv), (a)` で `φ = γ ≫ β ≫ α`(Frobenius 型・pre-step・pull-back)
2. (i)(＋`Proposition 1.10, (v)`)により `γ` は同型 ⟹ `φ = β ≫ α` としてよい
3. `φ` が FSM 射であることから `α` は fiberwise-surjective
4. ★**`α` が mono であることを示す**(ここが山) ——
   `α ≫ ϵ₁ = α ≫ ϵ₂` から `Remark 1.1.1` で `degFr`・`Div` の一致を出し、
   分解と全射性で `ϵ₁, ϵ₂` を pull-back としてよいとし、
   ★★**「`β∗(Div(β))` の引き戻しを `ϵ₁, ϵ₂` 経由で足す」**という構成で
   pre-step `ζ` と `γ₁, γ₂` を作って `φ ≫ γ₁ = φ ≫ γ₂` に持ち込む
5. よって `α` は FSM 射、`Proposition 1.11, (vi)` で `Base(α)` も FSM 射
6. `𝒟` が FSMFF 型なので、`α` が同型でなければ `Base(α)` は従属する FSMI 射を持ち、
   `Proposition 1.11, (vi)` で `α` も従属する FSMI 射(かつ pull-back)を持つ
7. それは mid-adjoint の仮定に矛盾。よって `α` は同型で `φ` は pre-step

★★**4 の「adding the pull-backs of β∗(Div(β)) via ϵ₁, ϵ₂」が
検証役の言う「引用符で名前だけ与える型」である** ——
★**構成の名前が引用符の中にあるだけで、構成そのものは書かれていない。**
★**ここが (ii) の山であり、まだ登っていない。**

★★**`sorry` は置かない**(`Found/` の規律)。
-/

/-! ## ★(iv) —— 四角形を渡る prime-Frobenius 性

★★**引用を選び直した記録(事故 #3 の 10 度目)**: (iv) の主張は 3 行だが、
★**3 行すべてが `′` を含むため引用できない**(`′` は layout 抽出で落ちる既知の文字で、
12/58 文字・12/44 文字で停止した)。★**(iv) を挟む前後の行を引く。**

原文 (FrdI p.41):
> — where α1, . . . , αn, ψ are FSMI-morphisms [cf. §0] — it holds that n ≤N.

原文 (FrdI p.41):
> (v) Suppose further that Φ is non-dilating, and that φ is a non-pre-step

★**私は最初、原文の行を手で打った**。★★**それが誤りだった** ——
規律は「内容は PDF 目視、文字列は抽出テキストから」であり、
★**手打ちは抽出と一致しない**。ゲートが捕まえた。

★★**原文は (iv) の証明を書いていない。** ★実際、(i) さえあれば 2 段で終わる:

1. **次数**: `Remark 1.1.1`(次数は合成で掛かる)から
   `degFr(α)·degFr(β) = degFr(β′)·degFr(α′)`。`degFr(β) = degFr(β′)` と
   ★**`ℕ≥1` が消約的**であることから `degFr(α) = degFr(α′)`
2. **prime-Frobenius 性**: (i) により irreducible な射は 3 種類しかない。
   ★**step も pull-back も次数 1** である(pull-back は `Proposition 1.4, (ii)` で linear)。
   ★**次数が素数なら 1 ではない**ので、残るのは prime-Frobenius だけ

★★**「次数が同型類の不変量になっている」という `Proposition 1.10, (iv)` の構図が、
ここでも効いている** —— (iv) は本質的に「**次数が型を決める**」ことの言い換えである。

★**原文が「moreover」で後置している `degFr(α) = degFr(α′)` のほうが先**であり、
iff はその系である(検証役の測定と一致する)。
-/

include P in
/-- ★★**(i) の系 —— 次数が 1 でない irreducible 射は prime-Frobenius**。

★(i) の 3 種類のうち step と pull-back はどちらも次数 1 だから。
★**これが (iv) の両向きを一手で与える。** -/
theorem prop_1_14_primeFrob_of_irred_of_degFr_ne_one (G : Frobenioid P)
    (hiso : ∀ X : C, IsIsotropic P X) {A B : C} (φ : A ⟶ B)
    (hirr : IsIrreducibleMor φ) (hd : P.degFr φ ≠ 1) : IsPrimeFrobenius P φ := by
  rcases (prop_1_14_i P G hiso φ).mp hirr with h | ⟨hst, -⟩ | ⟨hpb, -⟩
  · exact h
  · exact absurd (show P.degFr φ = 1 from hst.1.1) hd
  · exact absurd (prop_1_4_ii_mp P G.core φ hpb).2 hd

include P in
/-- ★**(iv) の次数の部分** —— `ℕ≥1` の消約性だけで出る。 -/
theorem prop_1_14_iv_degFr {X Y W Z : C} (β : X ⟶ Y) (α : Y ⟶ Z)
    (α' : X ⟶ W) (β' : W ⟶ Z) (hsq : β ≫ α = α' ≫ β')
    (hdeg : P.degFr β = P.degFr β') : P.degFr α = P.degFr α' := by
  have h := congrArg P.degFr hsq
  rw [P.degFr_comp, P.degFr_comp, hdeg] at h
  -- `degFr α * degFr β' = degFr α' * degFr β'`(★`ℕ≥1` の積は可換)
  have h' : ((P.degFr α : ℕ+) : ℕ) * ((P.degFr β' : ℕ+) : ℕ)
      = ((P.degFr α' : ℕ+) : ℕ) * ((P.degFr β' : ℕ+) : ℕ) := by
    have := congrArg (fun n : ℕ+ => (n : ℕ)) h
    simpa [mul_comm] using this
  refine PNat.coe_injective (Nat.eq_of_mul_eq_mul_right ?_ h')
  exact (P.degFr β').pos

include P in
/-- ★★★**[FrdI] Proposition 1.14, (iv)**。 -/
theorem prop_1_14_iv (G : Frobenioid P) (hiso : ∀ X : C, IsIsotropic P X)
    {X Y W Z : C} (β : X ⟶ Y) (α : Y ⟶ Z) (α' : X ⟶ W) (β' : W ⟶ Z)
    (hsq : β ≫ α = α' ≫ β') (hdeg : P.degFr β = P.degFr β')
    (hα : IsIrreducibleMor α) (hα' : IsIrreducibleMor α') :
    P.degFr α = P.degFr α' ∧ (IsPrimeFrobenius P α ↔ IsPrimeFrobenius P α') := by
  have hd : P.degFr α = P.degFr α' := prop_1_14_iv_degFr P β α α' β' hsq hdeg
  refine ⟨hd, ?_, ?_⟩
  · intro hp
    refine prop_1_14_primeFrob_of_irred_of_degFr_ne_one P G hiso α' hα' (fun h1 => ?_)
    have h2 : ((P.degFr α : ℕ+) : ℕ) = 1 := by rw [hd, h1]; rfl
    exact Nat.not_prime_one (h2 ▸ hp.2)
  · intro hp
    refine prop_1_14_primeFrob_of_irred_of_degFr_ne_one P G hiso α hα (fun h1 => ?_)
    have h2 : ((P.degFr α' : ℕ+) : ℕ) = 1 := by rw [← hd, h1]; rfl
    exact Nat.not_prime_one (h2 ▸ hp.2)

end ABC3.Found.FrdI
