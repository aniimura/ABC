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
  refine ⟨prop_1_11_vii_fsm_of_coaPre P F G φ
    (prop_1_4_i P φ (fun Y _ => hiso Y)) hφ, ?_⟩
  intro X Y γ β α hfac hS
  exact absurd
    (prop_1_7_v_preStep P β α (prop_1_7_v_preStep P γ (β ≫ α) (hfac ▸ hφ)).2).1 hS

/-! ### ★★(ii) 十分性の第 4 段 —— 「adding」の中身

★★**原文が引用符で名前だけ与えた構成**:
> Now by “adding the pull-backs of β∗(Div(β)) via ϵ1, ϵ2” …, it follows that there
> exists a pre-step ζ : E →D such that there exist γ1, γ2 ∈Arr(C) satisfying
> ϵ1 ◦ζ = β ◦γ1, ϵ2 ◦ζ = β ◦γ2.

★**まず記法**: `β : A ⟶ Cc` なので `β∗ = Φ.map (Base β) : Φ(Cc_𝒟) → Φ(A_𝒟)` であり、
★**`Div β ∈ Φ(A_𝒟)` だから `β∗(Div β)` は型が合わない。**
`Definition 1.3, (iii), (d)` が `(ψ∗)⁻¹(Div ψ)` と書いているとおり、
★**原文の意図は `(β∗)⁻¹(Div β)`、すなわち `β` の不変量**である
(検証役が発見した。私も型で確認した)。

★★**そして「adding」は本当に必要である。**
`ϵ1` 用の条件は `invζ = Φ.map (Base ϵ1) invβ`、`ϵ2` 用は `Φ.map (Base ϵ2) invβ` で、
★**2 つは一般に違う**(`Base ϵ1 = Base ϵ2` は出ない ——
`Base ϵ1 ≫ Base α = Base ϵ2 ≫ Base α` から `Base α` を右から消すには
**mono** が要るが、totally epimorphic が与えるのは **epi** である)。
★**だから両方を**上から押さえる**ものが要る。それが和である。**

★**組み立て**:
1. `zᵢ := Φ.map (Base ϵᵢ) invβ` を実現する co-angular pre-step `ζᵢ : Eᵢ ⟶ D`
2. `z := z₁ + z₂` を実現する `ζ : E ⟶ D`
3. `MLe zᵢ z` から `coaPre_factor_of_mle` が `κᵢ ≫ ζᵢ = ζ` を与える
4. `prop_1_11_v_exists_pullBack` が `γᵢ′ ≫ β = ζᵢ ≫ ϵᵢ` を与える
5. `γᵢ := κᵢ ≫ γᵢ′` と置くと `γᵢ ≫ β = ζ ≫ ϵᵢ`
6. `φ = β ≫ α` は mono なので `γ₁ = γ₂`、よって `ζ ≫ ϵ₁ = ζ ≫ ϵ₂`
7. ★**`ζ` は epi**(`𝒞` が totally epimorphic)—— こちらは**左から**消せる ⟹ `ϵ₁ = ϵ₂`

★★**epi と mono の使い分けがこの証明の要点である** ——
6 で mono を使い、7 で epi を使う。★**どちらも「消す側」が違う。**
-/

include P in
/-- ★★★**(ii) 十分性の第 4 段の核** ——
`ϵ₁, ϵ₂` が pull-back のとき、`α` は(その 2 本について)mono である。

★**`φ = β ≫ α` が mono であることだけを使う**(FSM 射の mono 部分)。 -/
theorem prop_1_14_ii_epsilon_eq (G : Frobenioid P) (hiso : ∀ X : C, IsIsotropic P X)
    {A Cc B : C} (β : A ⟶ Cc) (hβs : IsPreStep P β) (α : Cc ⟶ B)
    (hmono : Mono (β ≫ α)) {Dd : C} (ε₁ ε₂ : Dd ⟶ Cc)
    (hε₁ : IsPullBack P ε₁) (hε₂ : IsPullBack P ε₂)
    (hsq : ε₁ ≫ α = ε₂ ≫ α) : ε₁ = ε₂ := by
  haveI hbβ : IsIso (P.Base β) := hβs.2
  have hβc : IsCoAngular P β := prop_1_4_i P β (fun Y _ => hiso Y)
  set invβ := Φ.map (inv (P.Base β)) (P.Div β) with hinvβ
  -- 段 1: 各 `ϵᵢ` の条件を実現する `ζᵢ`
  obtain ⟨E₁, ζ₁, hζ₁c, hζ₁s, hζ₁inv⟩ :=
    coaPre_realize_over P G Dd (Φ.map (P.Base ε₁) invβ)
  obtain ⟨E₂, ζ₂, hζ₂c, hζ₂s, hζ₂inv⟩ :=
    coaPre_realize_over P G Dd (Φ.map (P.Base ε₂) invβ)
  -- 段 2: ★**和**を実現する `ζ`
  obtain ⟨E, ζ, hζc, hζs, hζinv⟩ :=
    coaPre_realize_over P G Dd (Φ.map (P.Base ε₁) invβ + Φ.map (P.Base ε₂) invβ)
  haveI hbζ₁ : IsIso (P.Base ζ₁) := hζ₁s.2
  haveI hbζ₂ : IsIso (P.Base ζ₂) := hζ₂s.2
  haveI hbζ : IsIso (P.Base ζ) := hζs.2
  -- 段 3: `MLe` から因子分解
  obtain ⟨κ₁, -, -, hκ₁⟩ := coaPre_factor_of_mle P G ζ₁ hζ₁c hζ₁s ζ hζc hζs
    (by rw [hζ₁inv, hζinv]; exact ⟨Φ.map (P.Base ε₂) invβ, rfl⟩)
  obtain ⟨κ₂, -, -, hκ₂⟩ := coaPre_factor_of_mle P G ζ₂ hζ₂c hζ₂s ζ hζc hζs
    (by rw [hζ₂inv, hζinv]; exact ⟨Φ.map (P.Base ε₁) invβ, by rw [add_comm]⟩)
  -- 段 4: `Proposition 1.11, (v)` の pull-back の場合
  obtain ⟨γ₁', hγ₁'⟩ :=
    prop_1_11_v_exists_pullBack P G ε₁ hε₁ ζ₁ hζ₁c hζ₁s β hβc hβs (by rw [hζ₁inv])
  obtain ⟨γ₂', hγ₂'⟩ :=
    prop_1_11_v_exists_pullBack P G ε₂ hε₂ ζ₂ hζ₂c hζ₂s β hβc hβs (by rw [hζ₂inv])
  -- 段 5: `γᵢ := κᵢ ≫ γᵢ′`
  have hγ₁ : (κ₁ ≫ γ₁') ≫ β = ζ ≫ ε₁ := by
    rw [Category.assoc, hγ₁', ← Category.assoc, hκ₁]
  have hγ₂ : (κ₂ ≫ γ₂') ≫ β = ζ ≫ ε₂ := by
    rw [Category.assoc, hγ₂', ← Category.assoc, hκ₂]
  -- 段 6: `φ = β ≫ α` は mono
  have hmm : (κ₁ ≫ γ₁') ≫ (β ≫ α) = (κ₂ ≫ γ₂') ≫ (β ≫ α) := by
    rw [← Category.assoc, hγ₁, ← Category.assoc, hγ₂, Category.assoc, Category.assoc, hsq]
  have hgeq : κ₁ ≫ γ₁' = κ₂ ≫ γ₂' := (cancel_mono (β ≫ α)).mp hmm
  -- 段 7: `ζ` は epi
  haveI : Epi ζ := P.totEpiC _ _ _
  refine (cancel_epi ζ).mp ?_
  rw [← hγ₁, ← hγ₂, hgeq]

include P in
/-- ★**(ii) 十分性の第 3 段** —— `φ = β ≫ α` が fiberwise-surjective なら
`α` もそうである。

★★**原文の「it follows formally」は本当に formal だった** ——
`φ` の証人 `δ_A` に `β` を後置するだけ(`δ_Cc := δ_A ≫ β`)。 -/
theorem prop_1_14_ii_alpha_fs {A Cc B : C} (β : A ⟶ Cc) (α : Cc ⟶ B)
    (hfs : IsFiberwiseSurjective (β ≫ α)) : IsFiberwiseSurjective α := by
  intro Z γ
  obtain ⟨Dd, δA, δZ, hsq⟩ := hfs γ
  exact ⟨Dd, δA ≫ β, δZ, by rw [Category.assoc, hsq]⟩

/-! ### ★★(ii) 十分性の第 4′ 段 —— 「WLOG `ϵ₁, ϵ₂` は pull-back」

★**原文**(p.42、目視):
> hence, by applying the factorization of Definition 1.3, (iv), (a) [and the total
> epimorphicity of C; cf. also Definition 1.3, (ii), and the equivalences of categories
> of Definition 1.3, (iii), (d)], we may assume without loss of generality [from the
> point of view of showing that α is a monomorphism] that ϵ1, ϵ2 are pull-back morphisms.

★★**原文より短い道があった。**
★原文は `ϵᵢ` から Frobenius 因子を剥がして pull-back に落とす道を取るが、
★**`α` が pull-back であること**を使えば、
「`α` が mono」は ★**「`Base α` が `𝒟` で mono」に帰着する**
(`Proposition 1.11, (vi)`)。

★**そして `Base α` の mono 性は、底の 2 射を `𝒞` へ持ち上げて
`prop_1_14_ii_epsilon_eq` に渡せば出る** ——
持ち上げは `Definition 1.3, (i), (c)` と **`α` の pull-back 全射性**の 2 つで、
★**どちらも「同じ対象 `U` の上で」できる**(これが要点。
`plBk_realize` を 2 回使うと対象が別々になってしまう)。
-/

include P in
/-- ★★★**(ii) 十分性の第 4 段(完全形)** —— `Base α` は mono。

★`Base α` の 2 射 `u`, `v` を、**同じ対象 `U`** の上の pull-back 射に持ち上げる:
`u` は `plBk_realize` で、★**`v` は `α` の pull-back 全射性で**(対象を変えずに)。 -/
theorem prop_1_14_ii_base_alpha_mono (F : FrobenioidCore P) (G : Frobenioid P)
    (hiso : ∀ X : C, IsIsotropic P X) {A Cc B : C} (β : A ⟶ Cc) (hβs : IsPreStep P β)
    (α : Cc ⟶ B) (hαpb : IsPullBack P α) (hmono : Mono (β ≫ α)) :
    Mono (P.Base α) := by
  refine ⟨fun {Y} u v huv => ?_⟩
  obtain ⟨U, pu, hpupb, θ, hθ⟩ := plBk_realize P F Cc u
  obtain ⟨-, hsurj⟩ := hαpb U
  obtain ⟨w, hw⟩ := hsurj ⟨(pu ≫ α, θ.hom ≫ v), by
    show P.Base (pu ≫ α) = (θ.hom ≫ v) ≫ P.Base α
    rw [P.Base_comp, hθ, Category.assoc, huv, ← Category.assoc]⟩
  have hw' := congrArg Subtype.val hw
  have hwα : w ≫ α = pu ≫ α := congrArg Prod.fst hw'
  have hbw : P.Base w = θ.hom ≫ v := congrArg Prod.snd hw'
  have hwpb : IsPullBack P w := by
    have hc : IsPullBack P (w ≫ α) := by rw [hwα]; exact IsPullBack.comp P hpupb hαpb
    obtain ⟨hlb, hlin⟩ := prop_1_4_ii_mp P F (w ≫ α) hc
    obtain ⟨⟨hwco, hwlin⟩, -⟩ := prop_1_7_v_coAngularLinear P F w α hlb.1 hlin
    exact (prop_1_4_ii P F w).mpr
      ⟨⟨hwco, (prop_1_7_v_isometric P w α hlb.2).1⟩, hwlin⟩
  have heq : w = pu :=
    prop_1_14_ii_epsilon_eq P G hiso β hβs α hmono w pu hwpb hpupb hwα
  have hbb : θ.hom ≫ v = θ.hom ≫ u := by rw [← hbw, heq, hθ]
  exact ((cancel_epi θ.hom).mp hbb).symm

/-! ### ★★(ii) 十分性の第 5–7 段 —— FSMFF 型から矛盾を出す

★**原文**(p.42、目視):
> Thus, it follows [cf. Proposition 1.11, (vi)] that Base(α) is an FSM-morphism
> of D. Since, however, we are operating under the assumption that D is of FSMFF-
> type, it follows that if α is not an isomorphism, then Base(α) admits a subordinate
> FSMI-morphism, which implies [cf. Proposition 1.11, (vi)] that α admits a subordinate
> FSMI-morphism [which is also a pull-back morphism].

★★**「which implies」が持ち上げの段である** —— 底の分解を `𝒞` へ持ち上げるのに
`Definition 1.3, (i), (c)`(`plBk_realize`)と **pull-back の全射性**を使う。

★★**持ち上げは同型を除いてしか底を再現しない** ——
`plBk_realize` が返す `θ : Base E ≅ Z` がそれである。
★**だから「irreducible は同型との合成で保たれる」(`IsIrreducibleMor.comp_isIso`)が要る。**

★**最後の一歩**: `a₁` は irreducible で pre-step でない(底が同型でないから)ので、
mid-adjoint の仮定が `a₁` を同型にする。★**すると底の `f₁` も同型になり、
`f₁` が FSMI(したがって irreducible、したがって非同型)であることに矛盾する。**
-/

include P in
/-- ★★★**(ii) 十分性の第 5–7 段** —— `α` は同型である。

★`φ = β ≫ α`(`α` は FSM な pull-back)で `φ` が非 pre-step irreducible 射に
mid-adjoint なら、`α` は同型。 -/
theorem prop_1_14_ii_alpha_isIso (F : FrobenioidCore P) (hFSMFF : IsOfFSMFFType D)
    {A Cc B : C} (β : A ⟶ Cc) (α : Cc ⟶ B) (hαpb : IsPullBack P α)
    (hαfsm : IsFSMMorphism α) (φ : A ⟶ B) (hfac : φ = β ≫ α)
    (hmid : IsMidAdjoint (irredNonPreStep P) φ) : IsIso α := by
  by_contra hni
  -- 段 5: `Base α` は FSM 射で、同型ではない
  have hbfsm : IsFSMMorphism (P.Base α) := (prop_1_11_vi_fsm P F α hαpb).mp hαfsm
  have hbni : ¬ IsIso (P.Base α) := fun h =>
    hni (isIso_of_isPullBack_of_isBaseIso P F α hαpb h)
  -- 段 6: FSMFF 型の条件 (a) から従属する FSMI 射
  obtain ⟨n, hnpos, hchain⟩ := hFSMFF.1 _ _ (P.Base α) hbfsm hbni
  obtain ⟨Z, f₁, g, hbfac, hf₁⟩ := hchain.exists_first hnpos
  -- 持ち上げ: `Definition 1.3, (i), (c)`
  obtain ⟨E, p, hppb, θ, hθ⟩ := plBk_realize P F B g
  obtain ⟨-, hsurj⟩ := hppb Cc
  obtain ⟨a₁, ha₁⟩ := hsurj ⟨(α, f₁ ≫ θ.inv), by
    show P.Base α = (f₁ ≫ θ.inv) ≫ P.Base p
    rw [hθ, ← Category.assoc, Category.assoc f₁ θ.inv θ.hom, θ.inv_hom_id,
      Category.comp_id, hbfac]⟩
  have ha₁' := congrArg Subtype.val ha₁
  have hsq : a₁ ≫ p = α := congrArg Prod.fst ha₁'
  have hba₁ : P.Base a₁ = f₁ ≫ θ.inv := congrArg Prod.snd ha₁'
  -- `a₁` は pull-back
  obtain ⟨hαlb, hαlin⟩ := prop_1_4_ii_mp P F α hαpb
  have hcoc : IsCoAngular P (a₁ ≫ p) := by rw [hsq]; exact hαlb.1
  have hlinc : IsLinear P (a₁ ≫ p) := by rw [hsq]; exact hαlin
  obtain ⟨⟨ha₁co, ha₁lin⟩, -⟩ := prop_1_7_v_coAngularLinear P F a₁ p hcoc hlinc
  have ha₁isom : IsIsometric P a₁ :=
    (prop_1_7_v_isometric P a₁ p (by rw [hsq]; exact hαlb.2)).1
  have ha₁pb : IsPullBack P a₁ := (prop_1_4_ii P F a₁).mpr ⟨⟨ha₁co, ha₁isom⟩, ha₁lin⟩
  -- `a₁` は irreducible
  have ha₁irr : IsIrreducibleMor a₁ := by
    refine (prop_1_11_vi_irred P F a₁ ha₁pb).mpr ?_
    rw [hba₁]
    exact hf₁.2.comp_isIso θ.inv
  -- `a₁` は pre-step でない
  have ha₁np : ¬ IsPreStep P a₁ := by
    intro hps
    refine hf₁.2.1 ?_
    haveI : IsIso (P.Base a₁) := hps.2
    have hf : f₁ = P.Base a₁ ≫ θ.hom := by rw [hba₁, Category.assoc, θ.inv_hom_id,
      Category.comp_id]
    rw [hf]
    infer_instance
  -- 段 7: mid-adjoint の仮定から矛盾
  haveI := hmid Cc E β a₁ p (by rw [hfac, ← hsq]) ⟨ha₁irr, ha₁np⟩
  refine hf₁.2.1 ?_
  have hf : f₁ = P.Base a₁ ≫ θ.hom := by
    rw [hba₁, Category.assoc, θ.inv_hom_id, Category.comp_id]
  haveI : IsIso (P.Base a₁) := by
    rw [show P.Base a₁ = P.Base a₁ from rfl]
    exact ⟨⟨P.Base (inv a₁), by rw [← P.Base_comp, IsIso.hom_inv_id, P.Base_id],
      by rw [← P.Base_comp, IsIso.inv_hom_id, P.Base_id]⟩⟩
  rw [hf]
  infer_instance

/-! ### ★★(ii) 十分性の第 1–2 段と組み立て -/

include P in
/-- ★prime-Frobenius の合成は「同型」か「prime-Frobenius 因子を持つ」かのいずれか。 -/
theorem IsPrimeFrobComposite.isIso_or_exists {A B : C} {γ : A ⟶ B}
    (h : IsPrimeFrobComposite P γ) :
    IsIso γ ∨ ∃ (X : C) (ζ : A ⟶ X) (rest : X ⟶ B),
      γ = ζ ≫ rest ∧ IsPrimeFrobenius P ζ := by
  cases h with
  | iso _ hi => exact Or.inl hi
  | cons hζ _ => exact Or.inr ⟨_, _, _, rfl, hζ⟩

include P in
/-- ★★**(ii) 十分性の第 2 段** —— Frobenius 因子は同型である。

★**原文**: 「By assertion (i) [cf. also Proposition 1.10, (v)], it follows that
γ is an isomorphism」。

★**中身**: `Proposition 1.10, (v)` が `γ` を prime-Frobenius の合成に分解する。
同型でなければ prime-Frobenius 因子 `ζ` があり、それは (i) により irreducible で、
次数が素数なので pre-step ではない。★**mid-adjoint の仮定が `ζ` を同型にし、
次数 1 になって素数性に矛盾する。** -/
theorem prop_1_14_ii_gamma_isIso (F : FrobenioidCore P) (hiso : ∀ X : C, IsIsotropic P X)
    {A X B : C} (γ : A ⟶ X) (hγF : IsFrobeniusType P γ) (ρ : X ⟶ B)
    (φ : A ⟶ B) (hfac : φ = γ ≫ ρ)
    (hmid : IsMidAdjoint (irredNonPreStep P) φ) : IsIso γ := by
  rcases IsPrimeFrobComposite.isIso_or_exists P ((prop_1_10_v P F γ).mp hγF) with
    h | ⟨Z, ζ, rest, hzfac, hζ⟩
  · exact h
  · exfalso
    have hirr : IsIrreducibleMor ζ := prop_1_10_iv_mp P F (hiso A) ζ hζ
    have hnp : ¬ IsPreStep P ζ := by
      intro hps
      have h1 : ((P.degFr ζ : ℕ+) : ℕ) = 1 := by rw [show P.degFr ζ = 1 from hps.1]; rfl
      exact Nat.not_prime_one (h1 ▸ hζ.2)
    haveI := hmid A Z (𝟙 A) ζ (rest ≫ ρ)
      (by rw [hfac, hzfac, Category.id_comp, Category.assoc]) ⟨hirr, hnp⟩
    have h1 : ((P.degFr ζ : ℕ+) : ℕ) = 1 := by rw [degFr_of_isIso P ζ]; rfl
    exact Nat.not_prime_one (h1 ▸ hζ.2)

include P in
/-- ★★★**[FrdI] Proposition 1.14, (ii) の十分性** ——
FSM 射で非 pre-step irreducible 射に mid-adjoint なら pre-step。 -/
theorem prop_1_14_ii_mpr (F : FrobenioidCore P) (G : Frobenioid P)
    (hiso : ∀ X : C, IsIsotropic P X) (hFSMFF : IsOfFSMFFType D)
    {A B : C} (φ : A ⟶ B) (hfsm : IsFSMMorphism φ)
    (hmid : IsMidAdjoint (irredNonPreStep P) φ) : IsPreStep P φ := by
  obtain ⟨X, Y, γ, β, α, hfac, hγF, hβs, hαpb⟩ := F.arbFactor φ
  haveI hγi : IsIso γ := prop_1_14_ii_gamma_isIso P F hiso γ hγF (β ≫ α) φ hfac hmid
  have hfac' : φ = (γ ≫ β) ≫ α := by rw [hfac, Category.assoc]
  have hβ's : IsPreStep P (γ ≫ β) := IsPreStep.comp P (isPreStep_of_isIso P γ) hβs
  have hmono : Mono ((γ ≫ β) ≫ α) := by rw [← hfac']; exact hfsm.2
  have hfs : IsFiberwiseSurjective α :=
    prop_1_14_ii_alpha_fs P (γ ≫ β) α (by rw [← hfac']; exact hfsm.1)
  have hαmono : Mono α := (prop_1_11_vi_mono P F α hαpb).mpr
    (prop_1_14_ii_base_alpha_mono P F G hiso (γ ≫ β) hβ's α hαpb hmono)
  haveI : IsIso α :=
    prop_1_14_ii_alpha_isIso P F hFSMFF (γ ≫ β) α hαpb ⟨hfs, hαmono⟩ φ hfac' hmid
  rw [hfac']
  exact IsPreStep.comp P hβ's (isPreStep_of_isIso P α)

include P in
/-- ★★★**[FrdI] Proposition 1.14, (ii)** —— pre-step ⟺ FSM ＋ mid-adjoint。 -/
theorem prop_1_14_ii (F : FrobenioidCore P) (G : Frobenioid P)
    (hiso : ∀ X : C, IsIsotropic P X) (hFSMFF : IsOfFSMFFType D)
    {A B : C} (φ : A ⟶ B) :
    IsPreStep P φ ↔ (IsFSMMorphism φ ∧ IsMidAdjoint (irredNonPreStep P) φ) := by
  constructor
  · intro h
    obtain ⟨h1, h2⟩ := prop_1_14_ii_mp P F G hiso φ h
    exact ⟨h1, isMidAdjoint_irredNonPreStep_of_nonPreStep P φ h2⟩
  · rintro ⟨h1, h2⟩
    exact prop_1_14_ii_mpr P F G hiso hFSMFF φ h1 h2

/-! ### ★(旧)(ii) の十分性 —— 段取りの記録(2026-08-16)

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

★★**その中身は上の `prop_1_14_ii_epsilon_eq` で実装した**(2026-08-16)。
★**残るのは 1・2・3・5・6・7 段**である:

| 段 | 内容 | 状態 |
|---|---|---|
| 1 | `Definition 1.3, (iv), (a)` で `φ = γ ≫ β ≫ α` | 道具あり(`arbFactor`) |
| 2 | (i) で `γ` は同型 ⟹ `φ = β ≫ α` としてよい | 道具あり |
| 3 | `φ` が FSM ⟹ `α` は fiberwise-surjective | ★**実装した**(`prop_1_14_ii_alpha_fs`) |
| 4 | ★**`α` が mono** | ★**実装した**(`prop_1_14_ii_epsilon_eq`、ただし `ϵᵢ` が pull-back の場合) |
| 4′ | 「WLOG `ϵ₁, ϵ₂` は pull-back」の正当化 | ★**未実装** |
| 5–7 | `Base(α)` も FSM ⟹ FSMFF 型 ⟹ 従属 FSMI 射 ⟹ mid-adjoint と矛盾 | ★**実装した**(`prop_1_14_ii_alpha_isIso`) |

★★**`sorry` は置かない**(`Found/` の規律)。
-/

/-! ## ★(iii) —— 非 pre-step ⟺ FSMI 分解の長さが有界

原文 (FrdI p.41):
> (iii) Suppose that φ is irreducible. Then φ is a non-pre-step if and only if

原文 (FrdI p.41):
> of composites in C

★**主張の条件を型で書いておく**(証明は未了。下の測定を見よ)。
-/

/-- ★(iii) の条件 —— 「`φ` の後ろに FSMI 射を 1 本足したものを FSMI 射の合成に
分解したとき、その長さに一様な上界がある」。

★原文の `αn ◦ … ◦ α1 = ψ ◦ φ` は Lean では
`IsFSMIChain n χ` かつ `χ = φ ≫ ψ` である(合成の向きが逆)。 -/
def BoundedFSMIFactor {A B : C} (φ : A ⟶ B) : Prop :=
  ∃ N : ℕ, ∀ (n : ℕ) (E : C) (ψ : B ⟶ E) (χ : A ⟶ E),
    IsFSMI ψ → IsFSMIChain n χ → χ = φ ≫ ψ → n ≤ N

/-! ### ★(iii) の測定(2026-08-16)

★**原文の証明**(p.42–43、目視)は (i) の 3 分類に沿って場合分けする。

**(a) `φ` が irreducible な pre-step のとき、条件は偽**:
> by taking ψ to be a prime-Frobenius morphism of increasingly large Frobenius degree
> [cf. Proposition 1.10, (ii)]

★**中身**: `φ` は step で `Div φ` は irreducible((i))。次数 `p` の
prime-Frobenius `ψ` を後置すると `Div (φ ≫ ψ) = p · Div φ` になり、
★**`p` 個の step の合成に分解できる**ので鎖の長さが `p` とともに伸びる。
★**この分解を作るのに第1の圏同値が要る。**

★★**さらに「prime-Frobenius 射が FSM である」ことが要る**(FSMI の M の部分)。
★**mono は `Definition 1.3, (v), (a)` が pre-step にしか与えない** ——
Frobenius 型射の mono 性がどこから来るかは未確認。

**(b) `φ` が非 pre-step のとき、条件は真**:
> since ψ◦φ and ψ are FSM-morphisms, it thus follows formally that φ is also an
> FSM-morphism … Div(ψ ◦φ) is either zero or irreducible; since, moreover,
> degFr(ψ ◦φ) always divides a product of two prime numbers …, it thus follows that
> in any factorization of ψ ◦φ by FSMI-morphisms, all but three … are pull-back
> morphisms … this implies that factorizations of arbitrarily large length determine
> chains of FSMI-morphisms … originating from the projection to D of the domain of φ
> which are also of arbitrarily large length, a contradiction

★★**「all but three」の数え上げが山である** ——
次数の素因数が高々 2 つ、零因子の既約因子が高々 1 つ、合わせて 3。
★**残りはすべて pull-back なので、`Proposition 1.11, (vi)` で底へ落ち、
`𝒟` の FSMFF 型の条件 (b)(鎖の長さの上界)に矛盾する。**

★★**`sorry` は置かない。**
-/

/-! ### ★★★(a) の側は原文の仮定では閉じない —— 障害を定理にした(2026-08-16)

★★**検証役の調査で確定した。**

★**(a) の側は「prime-Frobenius 射が FSM である」を避けて通れない** ——
条件文は鎖の各因子だけでなく ★**`ψ` 自身にも FSMI を要求している**
(原文 l.2136「where α1, . . . , αn, **ψ** are FSMI-morphisms」)。
そして合成の次数は乗法的なので、`degFr (φ ≫ ψ) = p > 1` である以上、
★**鎖のどこかに次数 > 1 の因子が必ずあり、(i) によりそれは prime-Frobenius である。**

★★**ところが `Definition 1.3` で mono を与える条項は `preStepMono` **1 つだけ**であり、
pre-step にしか mono を与えない。**

★★**そして「prime-Frobenius ⟹ mono」は偽になりうる。** 下の定理がその機構である ——
`Frobenius-normalized` の等式 `ζ ≫ α^(degFr ζ) = α ≫ ζ` に
★**`𝒪^×(A)` の `p`-捻れ元** `α`(`α^p = 1`, `α ≠ 1`)を当てると
`𝟙 ≫ ζ = α ≫ ζ` になり、`ζ` は mono でない。

★★**障害は `𝒪^×` の `p`-捻れに正確に局在している。**

★**「原文が誤り」とは断定しない** —— 上の条件を満たす Frobenioid を
実際に構成してはいないからである。★**言えるのは
「(iii) の (a) の側は、書かれている仮定だけでは閉じない」ことである。**
-/

include P in
/-- ★★★**prime-Frobenius 射が mono でなくなる機構**。

`A` が Frobenius-normalized で、`𝒪^×(A)`(実は `𝒪^▷(A)` で足りる)に
`α^p = 1` なる `α ≠ 1` があれば、次数 `p` の base-identity Frobenius 型自己射は
★**mono ではない**。

★★**`Prop 1.14, (iii)` の (a) の側が原文の仮定では閉じない理由がこれである** ——
そこは「prime-Frobenius 射が FSMI(したがって mono)」を要求する。 -/
theorem not_mono_of_frobNormalized_of_torsion {A : C} (hfn : IsFrobeniusNormalized P A)
    (ζ : End A) (hζb : IsBaseIdentity P ζ) (α : End A) (hα : α ∈ OTri P A)
    (hαp : α ^ ((P.degFr (ζ : A ⟶ A) : ℕ+) : ℕ) = 1) (hne : α ≠ 1) :
    ¬ Mono (ζ : A ⟶ A) := by
  intro hm
  refine hne ?_
  have h := hfn ζ hζb α hα
  rw [hαp] at h
  -- `ζ ≫ 𝟙 = α ≫ ζ`
  have h2 : (𝟙 A : A ⟶ A) ≫ (ζ : A ⟶ A) = (α : A ⟶ A) ≫ (ζ : A ⟶ A) := by
    rw [Category.id_comp, ← h]
    show (ζ : A ⟶ A) = (ζ : A ⟶ A) ≫ (𝟙 A : A ⟶ A)
    rw [Category.comp_id]
  exact ((cancel_mono (ζ : A ⟶ A)).mp h2).symm

/-! ### ★★(b) の側 —— 数え上げ

★**検証役が詰めた 4 段**:
1. `Div (φ ≫ ψ)` は **0 か既約**(`φ` が isometry で、`Φ.map (Base φ)` が
   `Definition 1.1, (ii), (b)` により同型だから既約性を保つ)
2. `degFr (φ ≫ ψ)` は**高々 2 つの素数の積**((i) より irreducible 射の次数は 1 か素数)
3. 鎖の各因子を (i) の 3 分類に当てる ——
   ★**(a) 型は高々 2 本**(次数)、★**(b) 型は高々 1 本**(`Div` の既約性)、
   残りはすべて **(c) 型 ＝ pull-back**
4. (c) 型は `Proposition 1.11, (vi)` で `𝒟` の FSMI 射に落ち、
   ★**(a)(b) 型は底が同型なので隣に吸収できる** ⟹
   `𝒟` に長さ `≥ n − 3` の FSMI 鎖ができ、FSMFF 型の条件 (b) に当たる

★★**求める `N` は `N_𝒟 + 3` で、証明は構成的である。**
★**「非有界」と仮定するのは `n` だけで、他はすべて有界。**
-/

include P in
/-- ★**FSMI 射の 3 分類**((i) の系、不変量の形で)。

★(a) 次数が素数・`Div = 0`・底は同型 /
(b) 次数 1・`Div` が既約・底は同型 /
(c) 次数 1・`Div = 0`・★**底が `𝒟` の FSMI 射**。 -/
theorem isFSMI_trichotomy (F : FrobenioidCore P) (G : Frobenioid P)
    (hiso : ∀ X : C, IsIsotropic P X) {A B : C} (χ : A ⟶ B) (h : IsFSMI χ) :
    (Nat.Prime ((P.degFr χ : ℕ+) : ℕ) ∧ P.Div χ = 0 ∧ IsIso (P.Base χ)) ∨
    (P.degFr χ = 1 ∧ IsIrreducibleElt (P.Div χ) ∧ IsIso (P.Base χ)) ∨
    (P.degFr χ = 1 ∧ P.Div χ = 0 ∧ IsFSMI (P.Base χ)) := by
  rcases (prop_1_14_i P G hiso χ).mp h.2 with hpf | ⟨hst, hdiv⟩ | ⟨hpb, hbirr⟩
  · exact Or.inl ⟨hpf.2, hpf.1.1.2, hpf.1.2⟩
  · exact Or.inr (Or.inl ⟨hst.1.1, hdiv, hst.1.2⟩)
  · obtain ⟨hlb, hlin⟩ := prop_1_4_ii_mp P F χ hpb
    exact Or.inr (Or.inr ⟨hlin, hlb.2,
      (prop_1_11_vi_fsm P F χ hpb).mp h.1, hbirr⟩)

include P in
/-- ★★★**鎖の射影と数え上げ**(2026-08-16)。

`𝒞` の長さ `n` の FSMI 鎖 `χ` について、3 つの数 `k`, `a`, `b` があって
`n = k + a + b` であり:
- **`k`** —— `𝒟` に `Base A` から出る長さ `k` の FSMI 鎖がある((c) 型の本数)
- **`a`** —— `degFr χ` の素因数の個数(重複込み)((a) 型の本数)
- **`b`** —— `Div ≠ 0` の因子の本数((b) 型の本数)。
  ★`b = 0` なら `Div χ = 0`、`b ≥ 1` なら `Div χ ≠ 0`、
  ★**`b ≥ 2` なら `Div χ` は非零の 2 項の和に分かれる**

★★**最後の条件が「(b) 型は高々 1 本」を与える** ——
`Div χ` が既約なら非零の 2 項の和には分かれないからである。 -/
theorem chain_project (F : FrobenioidCore P) (G : Frobenioid P)
    (hiso : ∀ X : C, IsIsotropic P X) :
    ∀ {n : ℕ} {A B : C} {χ : A ⟶ B}, IsFSMIChain n χ →
      ∃ k a b : ℕ, n = k + a + b ∧
        (∃ (Z : D) (g : (P.toElem.obj A).base ⟶ Z), IsFSMIChain k g) ∧
        (((P.degFr χ : ℕ+) : ℕ).primeFactorsList.length = a) ∧
        (b = 0 → P.Div χ = 0) ∧ (0 < b → P.Div χ ≠ 0) ∧
        (2 ≤ b → ∃ x y : Φ.val (P.toElem.obj A).base,
          P.Div χ = x + y ∧ x ≠ 0 ∧ y ≠ 0) := by
  intro n A B χ h
  induction h with
  | @nil A₀ =>
      refine ⟨0, 0, 0, rfl, ⟨_, 𝟙 _, IsFSMIChain.nil⟩, ?_, fun _ => P.Div_id A₀,
        by omega, by omega⟩
      rw [show ((P.degFr (𝟙 A₀) : ℕ+) : ℕ) = 1 by rw [P.degFr_id]; rfl,
        Nat.primeFactorsList_one]
      rfl
  | @cons m A₁ B₁ E₁ u v hu hv ih =>
      obtain ⟨k, a, b, hn, ⟨Z, g, hg⟩, hdeg, hb0, hb1, hb2⟩ := ih
      have hdv : P.Div (u ≫ v)
          = Φ.map (P.Base u) (P.Div v) + ((P.degFr v : ℕ+) : ℕ) • P.Div u := P.Div_comp u v
      have hdegc : ((P.degFr (u ≫ v) : ℕ+) : ℕ)
          = ((P.degFr v : ℕ+) : ℕ) * ((P.degFr u : ℕ+) : ℕ) := by
        rw [P.degFr_comp]; rfl
      -- 鎖の前置(同型の場合)
      have hchain_iso : ∀ (e : (P.toElem.obj A₁).base ⟶ (P.toElem.obj B₁).base), IsIso e →
          ∃ (Z' : D) (g' : (P.toElem.obj A₁).base ⟶ Z'), IsFSMIChain k g' := by
        intro e he
        rcases Nat.eq_zero_or_pos k with hk | hk
        · exact ⟨_, 𝟙 _, hk ▸ IsFSMIChain.nil⟩
        · haveI := he
          exact ⟨Z, e ≫ g, hg.isIso_comp e hk⟩
      -- 素因数の個数
      have hdeglen : ∀ c : ℕ, ((P.degFr u : ℕ+) : ℕ).primeFactorsList.length = c →
          ((P.degFr (u ≫ v) : ℕ+) : ℕ).primeFactorsList.length = a + c := by
        intro c hc
        have hperm := Nat.perm_primeFactorsList_mul
          (a := ((P.degFr v : ℕ+) : ℕ)) (b := ((P.degFr u : ℕ+) : ℕ))
          (P.degFr v).2.ne' (P.degFr u).2.ne'
        rw [hdegc, hperm.length_eq, List.length_append, hdeg, hc]
      rcases isFSMI_trichotomy P F G hiso u hu with
        ⟨hp, hz, hbi⟩ | ⟨h1, hirr, hbi⟩ | ⟨h1, hz, hfsmi⟩
      · -- (a) 型
        haveI := hbi
        refine ⟨k, a + 1, b, by omega, hchain_iso _ hbi, ?_, ?_, ?_, ?_⟩
        · exact hdeglen 1 (by rw [Nat.primeFactorsList_prime hp]; rfl)
        · intro hb
          rw [hdv, hz, smul_zero, add_zero, hb0 hb, map_zero]
        · intro hb hzero
          rw [hdv, hz, smul_zero, add_zero] at hzero
          exact hb1 hb (Φ.map_injective (P.Base u) (by rw [hzero, map_zero]))
        · intro hb
          obtain ⟨x, y, hxy, hx, hy⟩ := hb2 hb
          refine ⟨Φ.map (P.Base u) x, Φ.map (P.Base u) y, ?_, ?_, ?_⟩
          · rw [hdv, hz, smul_zero, add_zero, hxy, map_add]
          · exact fun hc => hx (Φ.map_injective (P.Base u) (by rw [hc, map_zero]))
          · exact fun hc => hy (Φ.map_injective (P.Base u) (by rw [hc, map_zero]))
      · -- (b) 型
        haveI := hbi
        have hdu : P.Div u ≠ 0 := hirr.1
        refine ⟨k, a, b + 1, by omega, hchain_iso _ hbi, ?_, by omega, ?_, ?_⟩
        · exact (hdeglen 0 (by rw [show ((P.degFr u : ℕ+) : ℕ) = 1 by rw [h1]; rfl,
            Nat.primeFactorsList_one]; rfl)).trans (by omega)
        · intro _ hzero
          rw [hdv] at hzero
          exact hdu (eq_zero_of_nsmul_eq_zero_of_isSharp
            (P.divisorial (P.toElem.obj A₁).base).2 (P.degFr v).2
            (eq_zero_of_add_eq_zero_of_isSharp
              (P.divisorial (P.toElem.obj A₁).base).2 hzero).2)
        · intro hb
          refine ⟨Φ.map (P.Base u) (P.Div v), ((P.degFr v : ℕ+) : ℕ) • P.Div u, hdv, ?_, ?_⟩
          · exact fun hc => hb1 (by omega)
              (Φ.map_injective (P.Base u) (by rw [hc, map_zero]))
          · exact fun hc => hdu (eq_zero_of_nsmul_eq_zero_of_isSharp
              (P.divisorial (P.toElem.obj A₁).base).2 (P.degFr v).2 hc)
      · -- (c) 型
        haveI := hfsmi
        refine ⟨k + 1, a, b, by omega, ⟨Z, P.Base u ≫ g, IsFSMIChain.cons hfsmi hg⟩,
          ?_, ?_, ?_, ?_⟩
        · exact (hdeglen 0 (by rw [show ((P.degFr u : ℕ+) : ℕ) = 1 by rw [h1]; rfl,
            Nat.primeFactorsList_one]; rfl)).trans (by omega)
        · intro hb
          rw [hdv, hz, smul_zero, add_zero, hb0 hb, map_zero]
        · intro hb hzero
          rw [hdv, hz, smul_zero, add_zero] at hzero
          exact hb1 hb (Φ.map_injective (P.Base u) (by rw [hzero, map_zero]))
        · intro hb
          obtain ⟨x, y, hxy, hx, hy⟩ := hb2 hb
          refine ⟨Φ.map (P.Base u) x, Φ.map (P.Base u) y, ?_, ?_, ?_⟩
          · rw [hdv, hz, smul_zero, add_zero, hxy, map_add]
          · exact fun hc => hx (Φ.map_injective (P.Base u) (by rw [hc, map_zero]))
          · exact fun hc => hy (Φ.map_injective (P.Base u) (by rw [hc, map_zero]))

include P in
/-- ★★★**[FrdI] Proposition 1.14, (iii) の `⟹`** ——
irreducible な非 pre-step は条件を満たす(`N := N_𝒟 + 3`)。

★**「非有界」と仮定するのは `n` だけで、他はすべて有界である**:
(a) 型は次数から高々 2 本、(b) 型は `Div` の既約性から高々 1 本、
(c) 型は `𝒟` の FSMFF 型の条件 (b) から高々 `N_𝒟` 本。 -/
theorem prop_1_14_iii_mp (F : FrobenioidCore P) (G : Frobenioid P)
    (hiso : ∀ X : C, IsIsotropic P X) (hFSMFF : IsOfFSMFFType D)
    {A B : C} (φ : A ⟶ B) (hirr : IsIrreducibleMor φ) (hnps : ¬ IsPreStep P φ) :
    BoundedFSMIFactor φ := by
  obtain ⟨N, hN⟩ := hFSMFF.2 (P.toElem.obj A).base
  refine ⟨N + 3, fun n E ψ χ hψ hchain hχ => ?_⟩
  obtain ⟨k, a, b, hn, ⟨Z, g, hg⟩, hdeg, -, -, hb2⟩ := chain_project P F G hiso hchain
  have hk : k ≤ N := hN k Z g hg
  -- `φ` の型: (i) で prime-Frobenius か pull-back(step は pre-step なので除外)
  have hφ : (Nat.Prime ((P.degFr φ : ℕ+) : ℕ) ∧ P.Div φ = 0 ∧ IsIso (P.Base φ)) ∨
      (P.Div φ = 0 ∧ IsPullBack P φ) := by
    rcases (prop_1_14_i P G hiso φ).mp hirr with hpf | ⟨hst, -⟩ | ⟨hpb, -⟩
    · exact Or.inl ⟨hpf.2, hpf.1.1.2, hpf.1.2⟩
    · exact absurd hst.1 hnps
    · exact Or.inr ⟨(prop_1_4_ii_mp P F φ hpb).1.2, hpb⟩
  -- ★`φ` は FSM(鎖が FSM で `ψ` が mono)
  have hφfsm : IsFSMMorphism φ := by
    refine isFSMMorphism_of_comp φ ψ ?_ hψ.1.2
    rw [← hχ]
    exact hchain.isFSMMorphism
  have hbfsm : IsFSMMorphism (P.Base φ) := by
    rcases hφ with ⟨-, -, hbi⟩ | ⟨-, hpb⟩
    · haveI := hbi; exact isFSMMorphism_of_isIso _
    · exact (prop_1_11_vi_fsm P F φ hpb).mp hφfsm
  have hbij : Function.Bijective (Φ.map (P.Base φ)) := Φ.fsmIso (P.Base φ) hbfsm
  -- ★`a ≤ 2` —— irreducible 射の次数は 1 か素数
  have hlen1 : ∀ {X Y : C} (ρ : X ⟶ Y), IsIrreducibleMor ρ →
      ((P.degFr ρ : ℕ+) : ℕ).primeFactorsList.length ≤ 1 := by
    intro X Y ρ hρ
    rcases (prop_1_14_i P G hiso ρ).mp hρ with hpf | ⟨hst, -⟩ | ⟨hpb, -⟩
    · rw [Nat.primeFactorsList_prime hpf.2]; simp
    · rw [show ((P.degFr ρ : ℕ+) : ℕ) = 1 by rw [show P.degFr ρ = 1 from hst.1.1]; rfl,
        Nat.primeFactorsList_one]; simp
    · rw [show ((P.degFr ρ : ℕ+) : ℕ) = 1 by
        rw [show P.degFr ρ = 1 from (prop_1_4_ii_mp P F ρ hpb).2]; rfl,
        Nat.primeFactorsList_one]; simp
  have ha : a ≤ 2 := by
    have h1 := hlen1 ψ hψ.2
    have h2 := hlen1 φ hirr
    have hperm := Nat.perm_primeFactorsList_mul
      (a := ((P.degFr ψ : ℕ+) : ℕ)) (b := ((P.degFr φ : ℕ+) : ℕ))
      (P.degFr ψ).2.ne' (P.degFr φ).2.ne'
    rw [← hdeg, hχ, show ((P.degFr (φ ≫ ψ) : ℕ+) : ℕ)
        = ((P.degFr ψ : ℕ+) : ℕ) * ((P.degFr φ : ℕ+) : ℕ) by rw [P.degFr_comp]; rfl,
      hperm.length_eq, List.length_append]
    omega
  -- ★`b ≤ 1` —— `Div χ` は 0 か既約
  have hb : b ≤ 1 := by
    by_contra hcon
    push_neg at hcon
    obtain ⟨x, y, hxy, hx, hy⟩ := hb2 (by omega)
    have hdφ : P.Div φ = 0 := by rcases hφ with ⟨-, h, -⟩ | ⟨h, -⟩ <;> exact h
    have hdχ : P.Div χ = Φ.map (P.Base φ) (P.Div ψ) := by
      rw [hχ, P.Div_comp, hdφ, smul_zero, add_zero]
    rcases isFSMI_trichotomy P F G hiso ψ hψ with ⟨-, hz, -⟩ | ⟨-, hirrψ, -⟩ | ⟨-, hz, -⟩
    · rw [hz, map_zero] at hdχ
      rw [hdχ] at hxy
      exact hx (eq_zero_of_add_eq_zero_of_isSharp
        (P.divisorial (P.toElem.obj A).base).2 hxy.symm).1
    · have hirrχ := isIrreducibleElt_of_bijective (Φ.map (P.Base φ)) hbij hirrψ
      rw [← hdχ] at hirrχ
      rcases hirrχ.2 x y hxy with h1 | h1
      · exact hx h1
      · exact hy h1
    · rw [hz, map_zero] at hdχ
      rw [hdχ] at hxy
      exact hx (eq_zero_of_add_eq_zero_of_isSharp
        (P.divisorial (P.toElem.obj A).base).2 hxy.symm).1
  omega

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

/-! ## ★(v) —— Div-identity prime-Frobenius 自己射の特徴づけ

原文 (FrdI p.41):
> (v) Suppose further that Φ is non-dilating, and that φ is a non-pre-step

原文 (FrdI p.41):
> condition holds: For every step α : A →B, there exists a non-pre-step irreducible

★**引用を選び直した記録(事故 #3 の 12 度目)**: 主張の最後の行
「morphism ψ : B →B′ and a step β : B →B′ such that ψ ◦α = β ◦α ◦φ.」は
★**`′` を含むため引用できない**(13/46 文字で停止)。1 行前を引く。
★**また手で打っていた** —— 抽出から差し込む規律を守れていない。

★**原文の必要性の証明**(p.43、目視):
> the necessity of the condition in the statement of assertion (v) [where we take ψ to
> be a prime-Frobenius morphism such that degFr(φ) = degFr(ψ)] follows immediately
> from Proposition 1.10, (i) [cf. also Definition 1.3, (ii); assertion (i); the first
> equivalence of categories of Definition 1.3, (iii), (d)].

★★**「follows immediately」の中身は 4 段である**:
1. `Definition 1.3, (ii)` が次数 `p` の Frobenius 型 `ψ : B ⟶ B′` を与える
   ((i) で irreducible、次数が素数なので非 pre-step)
2. `Proposition 1.10, (i)` が **一意の** `α′` で `α ≫ ψ = φ ≫ α′` を与える
3. ★**`φ` が Div-identity であることが効く** —— `Div α′ = p · Div α` になる
   (一般には `p · Φ(Base φ)⁻¹(Div α)`)
4. 第1の圏同値(`coaPre_factor_under_of_mle`)が `α ≫ β = α′` なる `β` を与え、
   ★**`Div β ≠ 0` は `(p−1)·Div α ≠ 0` から出る**(`p ≥ 2` と sharp)
-/

include P in
/-- ★★★**[FrdI] Proposition 1.14, (v) の必要性**。 -/
theorem prop_1_14_v_mp (F : FrobenioidCore P) (G : Frobenioid P)
    (hiso : ∀ X : C, IsIsotropic P X) {A : C} (φ : A ⟶ A)
    (hφF : IsFrobeniusType P φ) (hφp : Nat.Prime ((P.degFr φ : ℕ+) : ℕ))
    (hφd : IsDivIdentity P φ) {B : C} (α : A ⟶ B) (hαs : IsStep P α) :
    ∃ (B' : C) (ψ : B ⟶ B') (β : B ⟶ B'),
      IsIrreducibleMor ψ ∧ ¬ IsPreStep P ψ ∧ IsStep P β ∧
      α ≫ ψ = φ ≫ α ≫ β := by
  haveI hbφ : IsIso (P.Base φ) := hφF.2
  -- 段 1
  obtain ⟨B', ψ, hψF, hψd⟩ := F.frobDegSurj B (P.degFr φ)
  have hψp : IsPrimeFrobenius P ψ := ⟨hψF, by rw [hψd]; exact hφp⟩
  have hψirr : IsIrreducibleMor ψ := prop_1_10_iv_mp P F (hiso B) ψ hψp
  have hψnp : ¬ IsPreStep P ψ := by
    intro hps
    have h1 : ((P.degFr ψ : ℕ+) : ℕ) = 1 := by rw [show P.degFr ψ = 1 from hps.1]; rfl
    exact Nat.not_prime_one (h1 ▸ hψp.2)
  -- 段 2
  obtain ⟨α', hsq, -⟩ := prop_1_10_i_exists_given P F α φ hφF ψ hψF hψd.symm
  have hα's : IsPreStep P α' := prop_1_10_i_preStep_of P hφF hψF hψd.symm hsq hαs.1
  have hα'c : IsCoAngular P α' := prop_1_4_i P α' (fun Y _ => hiso Y)
  have hαc : IsCoAngular P α := prop_1_4_i P α (fun Y _ => hiso Y)
  -- 段 3: ★`Div-identity` が `Φ.map (inv (Base φ))` を恒等にする
  have hid : ∀ x : Φ.val (P.toElem.obj A).base, Φ.map (inv (P.Base φ)) x = x := by
    intro x
    have hb : Φ.map (P.Base φ) x = x := by
      rw [show Φ.map (P.Base φ) = Φ.map (P.Base (𝟙 A)) from hφd, P.Base_id, Φ.map_id]
    calc Φ.map (inv (P.Base φ)) x = Φ.map (inv (P.Base φ)) (Φ.map (P.Base φ) x) := by rw [hb]
      _ = x := by rw [← Φ.map_comp, IsIso.inv_hom_id, Φ.map_id]
  have hdα' : P.Div α' = ((P.degFr ψ : ℕ+) : ℕ) • P.Div α := by
    rw [prop_1_10_i_Div_formula' P α φ ψ α' hφF hψF.1.2 hsq, hid]
  -- `Div α ≠ 0`(step だから)
  have hdα : P.Div α ≠ 0 := fun h => hαs.2 (hiso A B α h hαs.1)
  -- 段 4
  obtain ⟨β, -, hβs, hβ⟩ := coaPre_factor_under_of_mle P G α hαc hαs.1 α' hα'c hα's
    (by rw [hdα']; exact mle_nsmul_self (P.degFr ψ).2 _)
  refine ⟨B', ψ, β, hψirr, hψnp, ⟨hβs, ?_⟩, by rw [hsq, hβ]⟩
  -- `β` は同型でない
  intro hβiso
  haveI := hβiso
  have hz : P.Div β = 0 := isIsometric_of_isIso P β
  have h2 : ((P.degFr ψ : ℕ+) : ℕ) • P.Div α = ((1 : ℕ)) • P.Div α := by
    rw [← hdα', ← hβ, P.Div_comp, hz, map_zero, zero_add,
      show P.degFr β = 1 from hβs.1]
    rfl
  have hp2 : 2 ≤ ((P.degFr ψ : ℕ+) : ℕ) := by rw [hψd]; exact hφp.two_le
  letI := isCancelAdd_of_isIntegralMonoid _ (P.divisorial (P.toElem.obj A).base).1.1
  obtain ⟨m, hm⟩ : ∃ m, ((P.degFr ψ : ℕ+) : ℕ) = m + 1 + 1 :=
    ⟨((P.degFr ψ : ℕ+) : ℕ) - 2, by omega⟩
  rw [hm, one_smul, succ_nsmul, succ_nsmul] at h2
  have h3 : (m • P.Div α + P.Div α) + P.Div α = 0 + P.Div α := by
    rw [zero_add]; exact h2
  exact hdα (eq_zero_of_add_eq_zero_of_isSharp
    (P.divisorial (P.toElem.obj A).base).2 (add_right_cancel h3)).2

end ABC3.Found.FrdI
