import ABC3.Found.FrdI.Prop33v

/-!
# [FrdI] Proposition 3.11, (i) —— group-like・isotropic・unit-trivial 型なら `𝒞 ≃ 𝔽_Φ`

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.73–p.74。

原文 (FrdI p.73):
> (Frobenioids of Isotropic, Unit-trivial, and Group-

原文 (FrdI p.73):
> like and isotropic type, it follows that every pre-step of Ci is an isomorphism [cf.

## ★この命題の規模(測定)

**3 条、主張は 5**:

| 条 | # | 内容 | 状態 |
|---|---|---|---|
| (i) | 1 | `𝒞ᵢ → 𝔽_{Φᵢ}` は圏同値 | ★**ここで実装** |
| (ii) | 2 | `Ψ` は base-isomorphism・pull-back 射・linear 射・Frobenius 型射を保つ | 未 |
| (iii) | 3 | `Ψ_Base : 𝒟₁ → 𝒟₂` が 1-一意に存在し、図式が 1-可換 | 未 |
| (iii) | 4 | 両方が圏同値 | 未 |
| (iii) | 5 | `𝒟ᵢ` が slim なら合成関手は rigid | 未 |

★**(ii)(iii) は `Theorem 3.4` と同じ道具(FSMI 射の圏論的特徴づけ)を要する**ので別途。

## ★(i) の証明(原文どおり)

1. `Proposition 3.3, (iii), (iv)` で **本質的全射かつ忠実**
2. ★**group-like かつ isotropic なら、`𝒞` の pre-step はすべて同型**
   —— `Φ` は divisorial ゆえ **sharp** で、group-like と合わせると `Φ(A)` は**零**。
   よってどの射も **isometric** で、isotropic 性が pre-step を同型にする
3. ゆえに `Definition 1.3, (i), (b)`(`preStepSpan`)から
   **Aut-ample かつ base-trivial** 型
4. `Proposition 3.3, (v)` で圏同値
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ)

/-! ## ★手 2 —— group-like ＋ divisorial なら `Φ(A)` は零 -/

include P in
/-- ★**group-like ＋ sharp なら零モノイド** —— group-like は「すべてが可逆」、
sharp は「可逆なら零」。 -/
theorem eq_zero_of_groupLike_of_sharp {M : Type w} [AddCommMonoid M]
    (hgl : IsGroupLike M) (hsh : IsSharp M) (a : M) : a = 0 :=
  hsh a ((isGroupLike_iff M).mp hgl a)

include P in
/-- ★★**group-like 型ならどの射も isometric** —— `Φ` が divisorial(ゆえ sharp)
なので `Φ(A)` が零になる。 -/
theorem isIsometric_of_groupLikeType (hgl : IsOfGroupLikeType P)
    {A B : C} (φ : A ⟶ B) : IsIsometric P φ :=
  eq_zero_of_groupLike_of_sharp P (hgl A) (P.divisorial _).2 (P.Div φ)

include P in
/-- ★★★**group-like ＋ isotropic 型なら pre-step はすべて同型**。

★原文の「every pre-step of `𝒞` is an isomorphism [cf. Propositions 1.4, (i); 1.8, (iii)]」。 -/
theorem isIso_of_preStep_of_groupLike (hgl : IsOfGroupLikeType P)
    (hiso : IsOfIsotropicType P) {A B : C} (φ : A ⟶ B) (hφ : IsPreStep P φ) : IsIso φ :=
  hiso A B φ (isIsometric_of_groupLikeType P hgl φ) hφ

/-! ## ★手 3 —— Aut-ample と base-trivial -/

include P in
/-- ★★**base-trivial 型** —— 底の同型は pre-step の span で持ち上がり
(`Definition 1.3, (i), (b)`)、その pre-step が同型だから。 -/
theorem isOfBaseTrivialType_of_groupLike (Fc : FrobenioidCore P)
    (hgl : IsOfGroupLikeType P) (hiso : IsOfIsotropicType P) : IsOfBaseTrivialType P := by
  intro A Dd hbi
  obtain ⟨α⟩ := hbi
  obtain ⟨X, σ, τ, hσ, hτ, -⟩ := Fc.preStepSpan A Dd α.hom inferInstance
  haveI : IsIso σ := isIso_of_preStep_of_groupLike P hgl hiso σ hσ
  haveI : IsIso τ := isIso_of_preStep_of_groupLike P hgl hiso τ hτ
  exact ⟨(asIso τ).symm ≪≫ asIso σ⟩

include P in
/-- ★★**Aut-ample 型** —— 底の自己同型を `preStepSpan` で持ち上げ、
両辺の pre-step が同型なので `𝒞` の自己同型が得られる。 -/
theorem isOfAutAmpleType_of_groupLike (Fc : FrobenioidCore P)
    (hgl : IsOfGroupLikeType P) (hiso : IsOfIsotropicType P) : IsOfAutAmpleType P := by
  intro A g hg
  haveI := hg
  obtain ⟨X, σ, τ, hσ, hτ, hgeq⟩ := Fc.preStepSpan A A g hg
  haveI : IsIso σ := isIso_of_preStep_of_groupLike P hgl hiso σ hσ
  haveI : IsIso τ := isIso_of_preStep_of_groupLike P hgl hiso τ hτ
  haveI hbσ : IsIso (P.Base σ) := hσ.2
  refine ⟨inv σ ≫ τ, inferInstance, ?_⟩
  rw [P.Base_comp, hgeq]
  congr 1
  refine (cancel_epi (P.Base σ)).mp ?_
  rw [← P.Base_comp, IsIso.hom_inv_id, P.Base_id, IsIso.hom_inv_id]

/-! ## ★本体 -/

include P in
/-- ★★★**[FrdI] Proposition 3.11, (i)** —— `𝒞` が isotropic・unit-trivial・
group-like 型なら `𝒞 → 𝔽_Φ` は**圏同値**。

★原文どおり `Proposition 3.3, (v)` に帰着する —— 残る 2 型
(Aut-ample と base-trivial)を「pre-step はすべて同型」から出す。 -/
theorem prop_3_11_i (Fc : FrobenioidCore P) (G : Frobenioid P)
    (hiso : IsOfIsotropicType P) (hut : IsOfUnitTrivialType P)
    (hgl : IsOfGroupLikeType P) : P.toElem.IsEquivalence :=
  prop_3_3_v_mpr P Fc G
    (isOfAutAmpleType_of_groupLike P Fc hgl hiso) hut
    (isOfBaseTrivialType_of_groupLike P Fc hgl hiso)

end ABC3.Found.FrdI
