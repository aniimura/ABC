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

end ABC3.Found.FrdI
