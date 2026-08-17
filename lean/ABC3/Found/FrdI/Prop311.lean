import ABC3.Found.FrdI.Prop33v
import ABC3.Found.FrdI.Thm34
import Mathlib.CategoryTheory.SingleObj
import Mathlib.CategoryTheory.Products.Basic

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
| (i) | 1 | `𝒞ᵢ → 𝔽_{Φᵢ}` は圏同値 | ★**実装** (`prop_3_11_i`) |
| (ii) | 2 | `Ψ` は base-isomorphism・pull-back 射・linear 射・Frobenius 型射を保つ | ★★**実装** (`prop_3_11_ii`) |
| (iii) | 3 | `Ψ_Base : 𝒟₁ → 𝒟₂` が 1-一意に存在し、図式が 1-可換 | 未 |
| (iii) | 4 | 両方が圏同値 | 未 |
| (iii) | 5 | `𝒟ᵢ` が slim なら合成関手は rigid | 未 |

## ★★(ii) の証明は **7 手**だった(原文は 6 行、2026-08-17 に実装)

| 手 | 内容 | 実装 |
|---|---|---|
| 1 | FSMI 自己射 ⟹ prime-Frobenius 自己射 | `isPrimeFrobenius_of_isFSMI_endo` |
| 2 | prime-Frobenius 自己射 ⟹ FSMI(★**`𝔽_Φ` で計算**) | `isFSMI_of_isPrimeFrobenius_endo` |
| 3 | `Ψ` は prime-Frobenius 射を保つ(★自己射に帰着) | `isPrimeFrobenius_map` |
| 4 | `Ψ` は Frobenius 型射を保つ(`Prop 1.10, (v)` の帰納法) | `isFrobeniusType_map_of_groupLike` |
| 5 | `Ψ` は底同型を保つ(`Prop 1.7, (ii)` ＋ pre-step は同型) | `isBaseIsomorphism_map_of_groupLike` |
| 6 | `Ψ` は linear 射を保つ(`Prop 1.7, (iii)` の minimal-adjoint) | `isLinear_map_of_groupLike` |
| 7 | `Ψ` は pull-back 射を保つ(`Prop 1.7, (ii)` の minimal-adjoint) | `isPullBack_map_of_groupLike` |

★★**手 2 が要**である。`𝔽_Φ` の射は `(base, div, deg)` の 3 つ組で、
`Φ` は group-like ＋ sharp なので `div` はつねに零 ——**実質 2 つ組**になり、
FSM 性が**明示的に計算できる**。★(i) の圏同値と FSMI の圏論性(`isFSMI_map_iff`)が
それを `𝒞` へ運ぶ。

★★**`Remark 3.4.1` は「次数を保つ」を仮定した。ここでは仮定しない** ——
(i) の圏同値があるので、FSMI の圏論性だけで次数の情報が回復する。

★**(iii) は別途**(`Ψ_Base` の構成・1-可換性・rigidity)。

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

/-! ## ★(ii) の第 1 歩 —— **FSMI 自己射は prime-Frobenius 自己射である**

原文 (FrdI p.73):
> follows that Di has no FSMI-endomorphisms [cf. §0], hence that a morphism of Ci

原文 (FrdI p.73):
> is an FSMI-endomorphism if and only if it is a prime-Frobenius endomorphism [cf.

★★**原文が挙げる 2 つの根拠がそのまま 2 つの場合を潰す**:

| `Proposition 1.14, (i)` の 3 分類 | 潰し方 |
|---|---|
| (a) prime-Frobenius | ★これが結論 |
| (b) `Div` が既約な step | ★**group-like なので `Div = 0`**、`0` は既約でない |
| (c) 底が既約な pull-back | ★★**底が `𝒟` の FSMI 自己射になる**(`Proposition 1.11, (vi)`)——`𝒟` は FSMFF 型なので矛盾 |

★★**`⟸` の側(prime-Frobenius 自己射は FSMI)は別問題**である ——
`Prop114.lean` の測定によれば「prime-Frobenius ⟹ mono」は一般には**偽になりうる**
(障害は `𝒪^×(A)` の `p`-捻れ)。★`Proposition 3.11` は **unit-trivial 型**を仮定するので
その障害は消えるが、**fiberwise-surjective の側は未測定**である。
-/

include P in
/-- ★★★**[FrdI] Proposition 3.11, (ii) の第 1 歩** ——
`𝒟` が FSMFF 型で `𝒞` が group-like・isotropic 型なら、
**`𝒞` の FSMI 自己射は prime-Frobenius 自己射である**。

★★`Proposition 1.14, (i)` の 3 分類のうち 2 つが潰れる:
group-like から `Div = 0`(既約でない)、
FSMFF 型から「底が既約な pull-back」が排除される。 -/
theorem isPrimeFrobenius_of_isFSMI_endo (Fc : FrobenioidCore P) (G : Frobenioid P)
    (hFSMFF : IsOfFSMFFType D) (hgl : IsOfGroupLikeType P)
    (hiso : ∀ X : C, IsIsotropic P X) {A : C} {φ : A ⟶ A} (h : IsFSMI φ) :
    IsPrimeFrobenius P φ := by
  rcases (prop_1_14_i P G hiso φ).mp h.2 with hpf | ⟨-, hd⟩ | ⟨hpb, hb⟩
  · exact hpf
  · -- ★(b): group-like なら `Div φ = 0`、しかし既約元は `0` でない
    exact absurd (isIsometric_of_groupLikeType P hgl φ) hd.1
  · -- ★(c): 底が `𝒟` の FSMI 自己射になる —— FSMFF 型に矛盾
    exact absurd
      (⟨(prop_1_11_vi_fsm P Fc φ hpb).mp h.1, hb⟩ : IsFSMI (P.Base φ))
      (not_isFSMI_endo_of_isOfFSMFFType hFSMFF (P.Base φ))

/-! ## ★(ii) の第 2 歩 —— **prime-Frobenius 自己射は FSMI である**

原文 (FrdI p.73):
> is an FSMI-endomorphism if and only if it is a prime-Frobenius endomorphism [cf.

★★**ここは `𝔽_Φ` で計算する。** (i) が `𝒞 ≃ 𝔽_Φ` を与えており、
`𝔽_Φ` の射は `(base, div, deg)` の 3 つ組なので **FSM 性が明示的に書ける**。
★`Φ` は group-like かつ sharp なので `Div` はつねに `0`——**3 つ組が 2 つ組になる**。

★★得た結果は `isFSMI_map_iff`(圏同値は FSMI を保ち、かつ**反射**する)で `𝒞` へ戻す。
-/

/-- ★★**`𝔽_Φ` では、底が fiberwise-surjective なら射も fiberwise-surjective**
(`Φ` の値がすべて零のとき)。

★**次数の側は `k₁ := γ.deg`, `k₂ := β.deg` と取れば `ℕ≥1` の可換性で揃う**。
★零因子の側は `Φ` が零なので自動。★**底の側だけが本質**である。 -/
theorem elemFrob_isFiberwiseSurjective_of_base
    (hzero : ∀ (X : D) (a : Φ.val X), a = 0)
    {B A : ElemFrobCat Φ} (β : B ⟶ A)
    (hb : IsFiberwiseSurjective (ElemFrobCat.Hom.base β)) :
    IsFiberwiseSurjective β := by
  intro Z γ
  obtain ⟨Dd, d₁, d₂, hd⟩ := hb (ElemFrobCat.Hom.base γ)
  refine ⟨⟨Dd⟩, ⟨d₁, 0, ElemFrobCat.Hom.deg γ⟩, ⟨d₂, 0, ElemFrobCat.Hom.deg β⟩, ?_⟩
  refine ElemFrobCat.Hom.ext ?_ ?_ ?_
  · exact hd
  · exact (hzero _ _).trans (hzero _ _).symm
  · exact mul_comm _ _

include P in
/-- ★★★**[FrdI] Proposition 3.11, (ii) の第 2 歩** —— group-like・isotropic・
unit-trivial 型なら、**prime-Frobenius 自己射は FSMI である**。

★★**3 成分を `𝔽_Φ` で潰す**:

| 成分 | `𝔽_Φ` での理由 |
|---|---|
| fiberwise-surjective | 底 `Base φ` は**同型**(Frobenius 型)ゆえ f.s.、あとは上の補題 |
| mono | 底が同型で次数は `ℕ≥1` で簡約的、零因子は零 |
| irreducible | `𝒞` 側で `Proposition 1.10, (iv)` を当ててから `Ψ` で送る |

★★**「prime-Frobenius ⟹ mono」は一般には偽になりうる**(`Prop114.lean` の測定、
障害は `𝒪^×(A)` の `p`-捻れ)。★**unit-trivial 型がその捻れをちょうど消す**——
ここでは (i) の圏同値を通して自動的に効いている。 -/
theorem isFSMI_of_isPrimeFrobenius_endo (Fc : FrobenioidCore P) (G : Frobenioid P)
    (hiso : IsOfIsotropicType P) (hut : IsOfUnitTrivialType P) (hgl : IsOfGroupLikeType P)
    {A : C} {φ : A ⟶ A} (h : IsPrimeFrobenius P φ) : IsFSMI φ := by
  haveI := prop_3_11_i P Fc G hiso hut hgl
  have hzero : ∀ (X : D) (a : Φ.val X), a = 0 := by
    intro X a
    obtain ⟨A₀, -, ⟨e⟩⟩ := Fc.baseSurj X
    have h0 : Φ.map e.hom a = 0 :=
      eq_zero_of_groupLike_of_sharp P (hgl A₀) (P.divisorial _).2 _
    exact Φ.map_injective e.hom (by rw [h0, map_zero])
  refine (isFSMI_map_iff P.toElem φ).mpr ?_
  haveI hbi : IsIso (P.Base φ) := h.1.2
  refine ⟨⟨?_, ?_⟩, isIrreducibleMor_map P.toElem
    ((prop_1_10_iv P Fc (hiso A) φ h.1).mp h)⟩
  · exact elemFrob_isFiberwiseSurjective_of_base hzero _
      (isFiberwiseSurjective_of_isIso (P.Base φ))
  · refine ⟨fun {X} f g hfg => ?_⟩
    have hb := congrArg ElemFrobCat.Hom.base hfg
    have hdg := congrArg ElemFrobCat.Hom.deg hfg
    refine ElemFrobCat.Hom.ext ?_ ((hzero _ _).trans (hzero _ _).symm) ?_
    · exact (cancel_mono (P.Base φ)).mp hb
    · exact mul_left_cancel hdg

include P in
/-- ★★★★**[FrdI] Proposition 3.11, (ii) の鍵** ——
**`𝒞` の自己射が FSMI であるのは prime-Frobenius であるとき、かつそのときに限る**。

原文 (FrdI p.73):
> is an FSMI-endomorphism if and only if it is a prime-Frobenius endomorphism [cf.

★★**これが「`Ψ` が prime-Frobenius 自己射を保つ」を与える** ——
FSMI は純粋に圏論的(`isFSMI_map_iff`)だからである。
★原文の「Thus, `Ψ` preserves the prime-Frobenius endomorphisms」の中身。 -/
theorem isFSMI_endo_iff_isPrimeFrobenius (Fc : FrobenioidCore P) (G : Frobenioid P)
    (hFSMFF : IsOfFSMFFType D) (hiso : IsOfIsotropicType P) (hut : IsOfUnitTrivialType P)
    (hgl : IsOfGroupLikeType P) {A : C} (φ : A ⟶ A) :
    IsFSMI φ ↔ IsPrimeFrobenius P φ :=
  ⟨isPrimeFrobenius_of_isFSMI_endo P Fc G hFSMFF hgl (fun X => hiso X),
   isFSMI_of_isPrimeFrobenius_endo P Fc G hiso hut hgl⟩

include P in
/-- ★★**prime-Frobenius 射は prime-Frobenius 自己射に abstractly equivalent** ——
group-like ＋ isotropic 型では **base-trivial** なので、終域を域に戻せる。

原文 (FrdI p.73):
> phisms [since every prime-Frobenius morphism is abstractly equivalent to a prime-

★Frobenius 型射は底同型なので域と終域は base-isomorphic、
base-trivial 型ならそれらは同型である。 -/
theorem exists_primeFrobenius_endo (Fc : FrobenioidCore P) (hgl : IsOfGroupLikeType P)
    (hiso : IsOfIsotropicType P) {A B : C} {φ : A ⟶ B} (h : IsPrimeFrobenius P φ) :
    ∃ k : B ≅ A, IsPrimeFrobenius P (φ ≫ k.hom) := by
  haveI : IsIso (P.Base φ) := h.1.2
  obtain ⟨k⟩ := isOfBaseTrivialType_of_groupLike P Fc hgl hiso A B ⟨asIso (P.Base φ)⟩
  refine ⟨k, IsFrobeniusType.comp P Fc h.1 (isFrobeniusType_of_isIso P k.hom), ?_⟩
  rw [P.degFr_comp, degFr_of_isIso, one_mul]
  exact h.2

/-! ## ★(ii) の第 3 歩 —— **圏同値は prime-Frobenius 射を保つ**

原文 (FrdI p.73):
> phisms [since every prime-Frobenius morphism is abstractly equivalent to a prime-

★★**自己射の場合に帰着してから FSMI で運ぶ** ——
FSMI は純粋に圏論的(`isFSMI_map_iff`)なので圏同値で移り、
第 1・第 2 歩の同値(`isFSMI_endo_iff_isPrimeFrobenius`)で両端を prime-Frobenius に翻訳する。
-/

section Preserve

universe v₃ u₃ v₄ u₄ w₃ w₄ v₅ u₅ v₆ u₆

variable {C₁ : Type u₃} [Category.{v₃} C₁] {C₂ : Type u₄} [Category.{v₄} C₂]
  (Ψ : C₁ ⥤ C₂) [Ψ.IsEquivalence]
  {D₁ : Type u₅} [Category.{v₅} D₁] {Φ₁ : MonoidOn.{v₅, u₅, w₃} D₁}
  (P₁ : PreFrobenioid C₁ Φ₁)
  {D₂ : Type u₆} [Category.{v₆} D₂] {Φ₂ : MonoidOn.{v₆, u₆, w₄} D₂}
  (P₂ : PreFrobenioid C₂ Φ₂)

/-- ★★★★**[FrdI] Proposition 3.11, (ii) の第 3 歩** ——
**圏同値は prime-Frobenius 射を保つ**。

★★**筋は 3 手**:
1. `exists_primeFrobenius_endo` で **prime-Frobenius 自己射に戻す**(base-trivial 型)
2. `isFSMI_endo_iff_isPrimeFrobenius` で **FSMI に翻訳**し、`isFSMI_map` で運ぶ
3. `𝒞₂` 側で同じ同値を逆向きに使い、最後に同型 `Ψ(k)` を剥がす

★**原文の「hence also prime-Frobenius morphisms」の 1 行が、この 3 手である。** -/
theorem isPrimeFrobenius_map
    (Fc₁ : FrobenioidCore P₁) (G₁ : Frobenioid P₁)
    (Fc₂ : FrobenioidCore P₂) (G₂ : Frobenioid P₂)
    (hF₁ : IsOfFSMFFType D₁) (hF₂ : IsOfFSMFFType D₂)
    (hiso₁ : IsOfIsotropicType P₁) (hut₁ : IsOfUnitTrivialType P₁)
    (hgl₁ : IsOfGroupLikeType P₁)
    (hiso₂ : IsOfIsotropicType P₂) (hut₂ : IsOfUnitTrivialType P₂)
    (hgl₂ : IsOfGroupLikeType P₂)
    {A B : C₁} {φ : A ⟶ B} (h : IsPrimeFrobenius P₁ φ) :
    IsPrimeFrobenius P₂ (Ψ.map φ) := by
  -- ★手 1: 自己射に戻す
  obtain ⟨k, hk⟩ := exists_primeFrobenius_endo P₁ Fc₁ hgl₁ hiso₁ h
  -- ★手 2: FSMI にして運ぶ
  have h2 : IsFSMI (Ψ.map (φ ≫ k.hom)) :=
    isFSMI_map Ψ
      ((isFSMI_endo_iff_isPrimeFrobenius P₁ Fc₁ G₁ hF₁ hiso₁ hut₁ hgl₁ _).mpr hk)
  -- ★手 3: `𝒞₂` 側で翻訳し、同型を剥がす
  have h3 : IsPrimeFrobenius P₂ (Ψ.map (φ ≫ k.hom)) :=
    (isFSMI_endo_iff_isPrimeFrobenius P₂ Fc₂ G₂ hF₂ hiso₂ hut₂ hgl₂ _).mp h2
  have hsplit : Ψ.map (φ ≫ k.hom) = Ψ.map φ ≫ Ψ.map k.hom := Ψ.map_comp _ _
  rw [hsplit] at h3
  haveI : IsIso (Ψ.map k.hom) := inferInstance
  have hback : Ψ.map φ = (Ψ.map φ ≫ Ψ.map k.hom) ≫ inv (Ψ.map k.hom) := by simp
  refine ⟨?_, ?_⟩
  · rw [hback]
    exact IsFrobeniusType.comp P₂ Fc₂ h3.1 (isFrobeniusType_of_isIso P₂ _)
  · have hd : P₂.degFr (Ψ.map φ ≫ Ψ.map k.hom) = P₂.degFr (Ψ.map φ) := by
      rw [P₂.degFr_comp, degFr_of_isIso, one_mul]
    rw [← hd]
    exact h3.2

/-- ★★★**[FrdI] Proposition 3.11, (ii) の第 4 歩** —— **圏同値は Frobenius 型射を保つ**。

原文 (FrdI p.73):
> Frobenius endomorphism]. But this implies that Ψ preserves the morphisms of Frobenius

★★鍵は `Proposition 1.10, (v)`(**Frobenius 型 ⟺ prime-Frobenius の合成**)。
★合成の**帰納法**をそのまま運ぶ: 基底(同型)は関手が保ち、
各因子は第 3 歩(`isPrimeFrobenius_map`)で保たれる。 -/
theorem isFrobeniusType_map_of_groupLike
    (Fc₁ : FrobenioidCore P₁) (G₁ : Frobenioid P₁)
    (Fc₂ : FrobenioidCore P₂) (G₂ : Frobenioid P₂)
    (hF₁ : IsOfFSMFFType D₁) (hF₂ : IsOfFSMFFType D₂)
    (hiso₁ : IsOfIsotropicType P₁) (hut₁ : IsOfUnitTrivialType P₁)
    (hgl₁ : IsOfGroupLikeType P₁)
    (hiso₂ : IsOfIsotropicType P₂) (hut₂ : IsOfUnitTrivialType P₂)
    (hgl₂ : IsOfGroupLikeType P₂)
    {A B : C₁} {φ : A ⟶ B} (h : IsFrobeniusType P₁ φ) :
    IsFrobeniusType P₂ (Ψ.map φ) := by
  have hcomp : IsPrimeFrobComposite P₁ φ :=
    (prop_1_10_v P₁ Fc₁ φ).mp h
  refine isFrobeniusType_of_isPrimeFrobComposite P₂ Fc₂ ?_
  clear h
  induction hcomp with
  | iso ψ hψ =>
    haveI := hψ
    exact IsPrimeFrobComposite.iso _ inferInstance
  | cons hφ _ ih =>
    rw [Ψ.map_comp]
    exact IsPrimeFrobComposite.cons
      (isPrimeFrobenius_map Ψ P₁ P₂ Fc₁ G₁ Fc₂ G₂ hF₁ hF₂ hiso₁ hut₁ hgl₁ hiso₂ hut₂ hgl₂ hφ)
      ih

/-- ★★★**[FrdI] Proposition 3.11, (ii) の第 5 歩** —— **圏同値は底同型を保つ**。

原文 (FrdI p.74):
> (iii)], it thus follows that Ψ preserves the pull-back morphisms [cf. Proposition 1.7,

★★`Proposition 1.7, (ii)` で底同型を `(Frobenius 型) ≫ (pre-step)` に分け、
★**group-like ＋ isotropic では pre-step はすべて同型**((i) の手 2)なので、
残るのは Frobenius 型の因子だけ。★第 4 歩がそれを運ぶ。

★★**`Remark 3.4.1` と同じ骨格**だが、あちらは「次数を保つ」を仮定した。
ここでは (i) の圏同値と FSMI の圏論性から**次数の仮定なしで**出る。 -/
theorem isBaseIsomorphism_map_of_groupLike
    (Fc₁ : FrobenioidCore P₁) (G₁ : Frobenioid P₁)
    (Fc₂ : FrobenioidCore P₂) (G₂ : Frobenioid P₂)
    (hF₁ : IsOfFSMFFType D₁) (hF₂ : IsOfFSMFFType D₂)
    (hiso₁ : IsOfIsotropicType P₁) (hut₁ : IsOfUnitTrivialType P₁)
    (hgl₁ : IsOfGroupLikeType P₁)
    (hiso₂ : IsOfIsotropicType P₂) (hut₂ : IsOfUnitTrivialType P₂)
    (hgl₂ : IsOfGroupLikeType P₂)
    {A B : C₁} {φ : A ⟶ B} (h : IsBaseIsomorphism P₁ φ) :
    IsBaseIsomorphism P₂ (Ψ.map φ) := by
  obtain ⟨X, γ, β, rfl, hγ, hβ⟩ := (prop_1_7_ii_baseIso_factor P₁ Fc₁ _).mp h
  haveI : IsIso β := isIso_of_preStep_of_groupLike P₁ hgl₁ hiso₁ β hβ
  rw [Ψ.map_comp]
  refine isBaseIsomorphism_comp P₂
    (isFrobeniusType_map_of_groupLike Ψ P₁ P₂ Fc₁ G₁ Fc₂ G₂ hF₁ hF₂
      hiso₁ hut₁ hgl₁ hiso₂ hut₂ hgl₂ hγ).2 ?_
  exact isBaseIsomorphism_of_isIso P₂ (Ψ.map β)

/-! ### ★保存の**反射** —— linear と pull-back に要る

★`Proposition 1.7` の特徴づけ(`linear ⟺ Frobenius 型に minimal-adjoint`、
`pull-back ⟺ 底同型に minimal-adjoint`)を運ぶには、分解を反射したあとで
**`𝒞₂` 側の因子が `𝒞₁` の因子から来ている**ことを言う必要がある。
★★そこで**保存を擬逆に当てて反射を作る**。 -/

variable (Fc₁ : FrobenioidCore P₁) (G₁ : Frobenioid P₁)
  (Fc₂ : FrobenioidCore P₂) (G₂ : Frobenioid P₂)
  (hF₁ : IsOfFSMFFType D₁) (hF₂ : IsOfFSMFFType D₂)
  (hiso₁ : IsOfIsotropicType P₁) (hut₁ : IsOfUnitTrivialType P₁)
  (hgl₁ : IsOfGroupLikeType P₁)
  (hiso₂ : IsOfIsotropicType P₂) (hut₂ : IsOfUnitTrivialType P₂)
  (hgl₂ : IsOfGroupLikeType P₂)

include Fc₁ G₁ Fc₂ G₂ hF₁ hF₂ hiso₁ hut₁ hgl₁ hiso₂ hut₂ hgl₂ in
/-- ★★**Frobenius 型は保存かつ反射される**。 -/
theorem isFrobeniusType_map_iff_of_groupLike {A B : C₁} (φ : A ⟶ B) :
    IsFrobeniusType P₁ φ ↔ IsFrobeniusType P₂ (Ψ.map φ) := by
  refine ⟨isFrobeniusType_map_of_groupLike Ψ P₁ P₂ Fc₁ G₁ Fc₂ G₂ hF₁ hF₂
    hiso₁ hut₁ hgl₁ hiso₂ hut₂ hgl₂, fun h => ?_⟩
  have h' : IsFrobeniusType P₁ (Ψ.inv.map (Ψ.map φ)) :=
    isFrobeniusType_map_of_groupLike Ψ.inv P₂ P₁ Fc₂ G₂ Fc₁ G₁ hF₂ hF₁
      hiso₂ hut₂ hgl₂ hiso₁ hut₁ hgl₁ h
  set eA := Ψ.asEquivalence.unitIso.app A with hEA
  set eB := Ψ.asEquivalence.unitIso.app B with hEB
  have hnat : φ ≫ eB.hom = eA.hom ≫ Ψ.inv.map (Ψ.map φ) :=
    Ψ.asEquivalence.unitIso.hom.naturality φ
  have h1 : IsFrobeniusType P₁ (φ ≫ eB.hom) := by
    rw [hnat]
    exact IsFrobeniusType.comp P₁ Fc₁ (isFrobeniusType_of_isIso P₁ eA.hom) h'
  have h2 := IsFrobeniusType.comp P₁ Fc₁ h1 (isFrobeniusType_of_isIso P₁ eB.inv)
  simpa using h2

include Fc₁ G₁ Fc₂ G₂ hF₁ hF₂ hiso₁ hut₁ hgl₁ hiso₂ hut₂ hgl₂ in
/-- ★★**底同型は保存かつ反射される**。 -/
theorem isBaseIsomorphism_map_iff_of_groupLike {A B : C₁} (φ : A ⟶ B) :
    IsBaseIsomorphism P₁ φ ↔ IsBaseIsomorphism P₂ (Ψ.map φ) := by
  refine ⟨isBaseIsomorphism_map_of_groupLike Ψ P₁ P₂ Fc₁ G₁ Fc₂ G₂ hF₁ hF₂
    hiso₁ hut₁ hgl₁ hiso₂ hut₂ hgl₂, fun h => ?_⟩
  have h' : IsBaseIsomorphism P₁ (Ψ.inv.map (Ψ.map φ)) :=
    isBaseIsomorphism_map_of_groupLike Ψ.inv P₂ P₁ Fc₂ G₂ Fc₁ G₁ hF₂ hF₁
      hiso₂ hut₂ hgl₂ hiso₁ hut₁ hgl₁ h
  set eA := Ψ.asEquivalence.unitIso.app A with hEA
  set eB := Ψ.asEquivalence.unitIso.app B with hEB
  have hnat : φ ≫ eB.hom = eA.hom ≫ Ψ.inv.map (Ψ.map φ) :=
    Ψ.asEquivalence.unitIso.hom.naturality φ
  have h1 : IsBaseIsomorphism P₁ (φ ≫ eB.hom) := by
    rw [hnat]
    exact isBaseIsomorphism_comp P₁ (isBaseIsomorphism_of_isIso P₁ eA.hom) h'
  have h2 := isBaseIsomorphism_comp P₁ h1 (isBaseIsomorphism_of_isIso P₁ eB.inv)
  simpa using h2

include Fc₁ G₁ Fc₂ G₂ hF₁ hF₂ hiso₁ hut₁ hgl₁ hiso₂ hut₂ hgl₂ in
/-- ★★★**[FrdI] Proposition 3.11, (ii) の第 6 歩** —— **圏同値は linear 射を保つ**。

原文 (FrdI p.73):
> type [cf. Proposition 1.10, (v)], hence also the linear morphisms [cf. Proposition

★★鍵は `Proposition 1.7, (iii)` の**圏論的な特徴づけ**
`IsLinear φ ↔ IsMinimalAdjoint (frobTypeProp) φ`。
★分解を反射し(`exists_factor_of_map_factor`)、因子の Frobenius 型性を
上の反射で `𝒞₁` へ戻してから、`𝒞₁` の minimal-adjoint 性を当てる。 -/
theorem isLinear_map_of_groupLike {A B : C₁} {φ : A ⟶ B} (h : IsLinear P₁ φ) :
    IsLinear P₂ (Ψ.map φ) := by
  rw [prop_1_7_iii_linear_minimalAdjoint P₂ Fc₂]
  rw [prop_1_7_iii_linear_minimalAdjoint P₁ Fc₁] at h
  intro X β α hfac hβ
  obtain ⟨X₀, β₀, α₀, e, hfac₀, hβe, hαe⟩ := exists_factor_of_map_factor Ψ φ β α hfac
  have hmapβ₀ : IsFrobeniusType P₂ (Ψ.map β₀) := by
    have : Ψ.map β₀ = β ≫ e.inv := by rw [hβe]; simp
    rw [this]
    exact IsFrobeniusType.comp P₂ Fc₂ hβ (isFrobeniusType_of_isIso P₂ e.inv)
  haveI : IsIso β₀ :=
    h X₀ β₀ α₀ hfac₀
      ((isFrobeniusType_map_iff_of_groupLike Ψ P₁ P₂ Fc₁ G₁ Fc₂ G₂ hF₁ hF₂
        hiso₁ hut₁ hgl₁ hiso₂ hut₂ hgl₂ β₀).mpr hmapβ₀)
  rw [hβe]
  infer_instance

include Fc₁ G₁ Fc₂ G₂ hF₁ hF₂ hiso₁ hut₁ hgl₁ hiso₂ hut₂ hgl₂ in
/-- ★★★**[FrdI] Proposition 3.11, (ii) の第 7 歩** —— **圏同値は pull-back 射を保つ**。

原文 (FrdI p.74):
> (iii)], it thus follows that Ψ preserves the pull-back morphisms [cf. Proposition 1.7,

★★鍵は `Proposition 1.7, (ii)` の**圏論的な特徴づけ**
`IsPullBack φ ↔ IsMinimalAdjoint (baseIsoProp) φ`。★第 6 歩と同じ骨格で、
`Frobenius 型`の代わりに`底同型`の反射を使う。 -/
theorem isPullBack_map_of_groupLike {A B : C₁} {φ : A ⟶ B} (h : IsPullBack P₁ φ) :
    IsPullBack P₂ (Ψ.map φ) := by
  rw [prop_1_7_ii_pullBack_minimalAdjoint P₂ Fc₂]
  rw [prop_1_7_ii_pullBack_minimalAdjoint P₁ Fc₁] at h
  intro X β α hfac hβ
  obtain ⟨X₀, β₀, α₀, e, hfac₀, hβe, hαe⟩ := exists_factor_of_map_factor Ψ φ β α hfac
  have hmapβ₀ : IsBaseIsomorphism P₂ (Ψ.map β₀) := by
    have : Ψ.map β₀ = β ≫ e.inv := by rw [hβe]; simp
    rw [this]
    exact isBaseIsomorphism_comp P₂ hβ (isBaseIsomorphism_of_isIso P₂ e.inv)
  haveI : IsIso β₀ :=
    h X₀ β₀ α₀ hfac₀
      ((isBaseIsomorphism_map_iff_of_groupLike Ψ P₁ P₂ Fc₁ G₁ Fc₂ G₂ hF₁ hF₂
        hiso₁ hut₁ hgl₁ hiso₂ hut₂ hgl₂ β₀).mpr hmapβ₀)
  rw [hβe]
  infer_instance

include Fc₁ G₁ Fc₂ G₂ hF₁ hF₂ hiso₁ hut₁ hgl₁ hiso₂ hut₂ hgl₂ in
/-- ★★★★★**[FrdI] Proposition 3.11, (ii)** —— `Ψ` は **base-isomorphism・
pull-back 射・linear 射・Frobenius 型射**を保つ。

原文 (FrdI p.73):
> (ii) Ψ preserves base-isomorphisms, pull-back morphisms, linear mor-

★★**4 主張が揃った。** 原文の証明の 7 手をすべて開いた:
FSMI 自己射の同定(第 1・2 歩)→ prime-Frobenius 射の保存(第 3 歩)→
Frobenius 型(第 4 歩)→ 底同型(第 5 歩)→ linear(第 6 歩)→ pull-back(第 7 歩)。 -/
theorem prop_3_11_ii {A B : C₁} (φ : A ⟶ B) :
    (IsBaseIsomorphism P₁ φ → IsBaseIsomorphism P₂ (Ψ.map φ)) ∧
      (IsPullBack P₁ φ → IsPullBack P₂ (Ψ.map φ)) ∧
      (IsLinear P₁ φ → IsLinear P₂ (Ψ.map φ)) ∧
      (IsFrobeniusType P₁ φ → IsFrobeniusType P₂ (Ψ.map φ)) :=
  ⟨isBaseIsomorphism_map_of_groupLike Ψ P₁ P₂ Fc₁ G₁ Fc₂ G₂ hF₁ hF₂
     hiso₁ hut₁ hgl₁ hiso₂ hut₂ hgl₂,
   isPullBack_map_of_groupLike Ψ P₁ P₂ Fc₁ G₁ Fc₂ G₂ hF₁ hF₂
     hiso₁ hut₁ hgl₁ hiso₂ hut₂ hgl₂,
   isLinear_map_of_groupLike Ψ P₁ P₂ Fc₁ G₁ Fc₂ G₂ hF₁ hF₂
     hiso₁ hut₁ hgl₁ hiso₂ hut₂ hgl₂,
   isFrobeniusType_map_of_groupLike Ψ P₁ P₂ Fc₁ G₁ Fc₂ G₂ hF₁ hF₂
     hiso₁ hut₁ hgl₁ hiso₂ hut₂ hgl₂⟩

end Preserve

/-! ## ★(iii) の第 1 歩 —— **`Φ` が零なら `𝔽_Φ ≌ 𝒟 × 𝒩`**

原文 (FrdI p.74):
> Finally, we consider assertion (iii). Write N for the one-object category whose

原文 (FrdI p.74):
> morphisms of Ci, where two morphisms of Ci are regarded as equivalent if they

★★**`Φ` が零モノイドなら `𝔽_Φ` の射 `(base, div, deg)` は `div` が消えて
`(base, deg)` の 2 つ組になる**——これはちょうど `𝒟 × 𝒩` の射である
(`𝒩` は自己射モノイドが `ℕ≥1` の 1 対象圏)。

★★★**合成則まで一致する**: `𝔽_Φ` は `(ψ ≫ φ).deg = φ.deg * ψ.deg`、
`SingleObj` は `comp x y := y * x`。★**原文が「≌」と書くところが `rfl` になる。**
-/

/-- ★★**`Φ` が零なら `𝔽_Φ → 𝒟 × 𝒩`** —— 射を `(base, deg)` に落とすだけ。

★`map_id` も `map_comp` も **`rfl`** である。 -/
def elemFrobToProd : ElemFrobCat Φ ⥤ D × SingleObj ℕ+ where
  obj A := (A.base, SingleObj.star ℕ+)
  map f := (ElemFrobCat.Hom.base f, ElemFrobCat.Hom.deg f)
  map_id _ := rfl
  map_comp _ _ := rfl

/-- ★★**忠実** —— `Φ` が零なので `div` 成分に情報が無い。 -/
theorem elemFrobToProd_faithful (hzero : ∀ (X : D) (a : Φ.val X), a = 0) :
    (elemFrobToProd (Φ := Φ)).Faithful where
  map_injective {A B f g} h := by
    refine ElemFrobCat.Hom.ext ?_ ((hzero _ _).trans (hzero _ _).symm) ?_
    · exact congrArg Prod.fst h
    · exact congrArg Prod.snd h

/-- ★★**充満** —— `(b, n)` は `⟨b, 0, n⟩` の像。 -/
theorem elemFrobToProd_full : (elemFrobToProd (Φ := Φ)).Full where
  map_surjective {A B} g := ⟨⟨g.1, 0, g.2⟩, rfl⟩

/-- ★★**本質的全射** —— 対象は `𝒟` の対象そのもの。 -/
theorem elemFrobToProd_essSurj : (elemFrobToProd (Φ := Φ)).EssSurj where
  mem_essImage Z := ⟨⟨Z.1⟩, ⟨Iso.refl _⟩⟩

/-- ★★★**[FrdI] Proposition 3.11, (iii) の第 1 歩** —— `Φ` が零モノイドなら
**`𝔽_Φ` は `𝒟 × 𝒩` と圏同値**である。

原文 (FrdI p.74):
> Finally, we consider assertion (iii). Write N for the one-object category whose

★★**(i) と合わせると `𝒞 ≌ 𝒟 × 𝒩`** になる。★原文はこの分解を使って
「`𝒟` は `𝒞` から圏論的に復元できる」(base-identity 自己射との合成が等しい射を同一視する)
と論じる。 -/
theorem elemFrobToProd_isEquivalence (hzero : ∀ (X : D) (a : Φ.val X), a = 0) :
    (elemFrobToProd (Φ := Φ)).IsEquivalence := by
  haveI := elemFrobToProd_faithful hzero
  haveI := elemFrobToProd_full (Φ := Φ)
  haveI := elemFrobToProd_essSurj (Φ := Φ)
  exact { }

include P in
/-- ★**`𝒞` の側へ運んだ形** —— group-like・isotropic・unit-trivial 型なら
`𝒞 ≌ 𝒟 × 𝒩`。

★(i)(`𝒞 ≌ 𝔽_Φ`)と上(`𝔽_Φ ≌ 𝒟 × 𝒩`)の合成である。 -/
theorem toElem_comp_prod_isEquivalence (Fc : FrobenioidCore P) (G : Frobenioid P)
    (hiso : IsOfIsotropicType P) (hut : IsOfUnitTrivialType P)
    (hgl : IsOfGroupLikeType P) :
    (P.toElem ⋙ elemFrobToProd (Φ := Φ)).IsEquivalence := by
  haveI := prop_3_11_i P Fc G hiso hut hgl
  haveI : (elemFrobToProd (Φ := Φ)).IsEquivalence := by
    refine elemFrobToProd_isEquivalence (fun X a => ?_)
    obtain ⟨A₀, -, ⟨e⟩⟩ := Fc.baseSurj X
    have h0 : Φ.map e.hom a = 0 :=
      eq_zero_of_groupLike_of_sharp P (hgl A₀) (P.divisorial _).2 _
    exact Φ.map_injective e.hom (by rw [h0, map_zero])
  infer_instance

/-! ### ★射影は `𝒟 × 𝒩` の第 1 成分そのもの

★★`ElemFrobCat.proj` は `obj A := A.base`、`map φ := φ.base`。
`elemFrobToProd` は `obj A := (A.base, ⋆)`、`map f := (f.base, f.deg)`。
★**したがって `elemFrobToProd ⋙ fst` は `proj` に等しい**——`rfl` である。 -/

theorem elemFrobToProd_comp_fst :
    elemFrobToProd (Φ := Φ) ⋙ CategoryTheory.Prod.fst D (SingleObj ℕ+)
      = ElemFrobCat.proj := rfl

/-- ★**`𝒟 → 𝒟 × 𝒩`**(第 2 成分は `⋆`、射は `1`)。

★★`SingleObj` では `𝟙 = 1` なので、**これが関手になるのは `1 * 1 = 1` だから**である。 -/
def toProdN : D ⥤ D × SingleObj ℕ+ :=
  CategoryTheory.Functor.prod' (𝟭 D) ((Functor.const D).obj (SingleObj.star ℕ+))

@[simp] theorem toProdN_obj (X : D) : (toProdN (D := D)).obj X = (X, SingleObj.star ℕ+) := rfl

/-- ★★**`𝔽_Φ` の対象は `(base, ⋆)` そのもの** —— `toProdN` は `elemFrobToProd` の
第 1 成分を取ってから戻すと元に戻る(対象について)。 -/
theorem toProdN_elemFrobToProd_obj (A : ElemFrobCat Φ) :
    (toProdN (D := D)).obj ((elemFrobToProd (Φ := Φ)).obj A).1
      = (elemFrobToProd (Φ := Φ)).obj A := rfl

/-! ## ★(iii) の第 2 歩 —— **`Ψ_Base : 𝒟₁ ⥤ 𝒟₂` の構成**

原文 (FrdI p.74):
> morphisms of Ci, where two morphisms of Ci are regarded as equivalent if they

★★**`𝒞ᵢ ≌ 𝒟ᵢ × 𝒩` を使えば `Ψ_Base` は合成で書ける**:

  `𝒟₁ --(X ↦ (X,⋆))--> 𝒟₁ × 𝒩 ≌ 𝒞₁ --Ψ--> 𝒞₂ ≌ 𝒟₂ × 𝒩 --fst--> 𝒟₂`

★★★**要点は「`𝒟₁ → 𝒟₁ × 𝒩` が射を次数 `1` に潰す」こと**である。
そのぶんの食い違い(`(Base φ, 1)` と `(Base φ, degFr φ)` の差)は
**base-identity 自己射**であり、★原文が (iii) で仮定する
「`Ψ` とその擬逆が base-identity 自己射を保つ」がちょうどそれを消す。
-/

section BaseFunctor

universe v₇ u₇ v₈ u₈ w₇ w₈ v₉ u₉ v₁₀ u₁₀

variable {C₁ : Type u₇} [Category.{v₇} C₁] {C₂ : Type u₈} [Category.{v₈} C₂]
  (Ψ : C₁ ⥤ C₂)
  {D₁ : Type u₉} [Category.{v₉} D₁] {Φ₁ : MonoidOn.{v₉, u₉, w₇} D₁}
  (P₁ : PreFrobenioid C₁ Φ₁)
  {D₂ : Type u₁₀} [Category.{v₁₀} D₂] {Φ₂ : MonoidOn.{v₁₀, u₁₀, w₈} D₂}
  (P₂ : PreFrobenioid C₂ Φ₂)

/-- ★**`𝒞 ⥤ 𝒟 × 𝒩`** —— `(i)` の圏同値の合成形。 -/
abbrev toProdCat : C₁ ⥤ D₁ × SingleObj ℕ+ := P₁.toElem ⋙ elemFrobToProd

/-- ★**射影は第 1 成分そのもの**(`rfl`)。 -/
theorem toProdCat_comp_fst :
    toProdCat P₁ ⋙ CategoryTheory.Prod.fst D₁ (SingleObj ℕ+) = P₁.proj := rfl

variable [hE₁ : (toProdCat P₁).IsEquivalence]

/-- ★★★**[FrdI] Proposition 3.11, (iii) の `Ψ_Base`**。

★`𝒟₁ → 𝒟₁ × 𝒩 ≌ 𝒞₁ --Ψ--> 𝒞₂ ≌ 𝒟₂ × 𝒩 → 𝒟₂` の合成。 -/
noncomputable def psiBase : D₁ ⥤ D₂ :=
  toProdN ⋙ (toProdCat P₁).inv ⋙ Ψ ⋙ toProdCat P₂ ⋙
    CategoryTheory.Prod.fst D₂ (SingleObj ℕ+)

/-- ★**対象の側の可換性** —— `(X, ⋆)` に戻す操作が `𝔽_Φ` の対象では**恒等**なので、
`Ψ_Base (Base A)` は `Base (Ψ (E₁⁻¹ (E₁ A)))` そのものである。 -/
theorem psiBase_obj (A : C₁) :
    (psiBase Ψ P₁ P₂).obj (P₁.proj.obj A)
      = ((toProdCat P₂).obj (Ψ.obj ((toProdCat P₁).inv.obj ((toProdCat P₁).obj A)))).1 := rfl

/-- ★★★**`𝒩` 成分だけの自己射は `E₁⁻¹` で base-identity 自己射に写る**。

★★これが `Ψ_Base` の 1-可換性の**要**である ——
`toProdN` が射を次数 `1` に潰すぶんの食い違いは、ちょうどこの形の自己射だから。

★`Functor.fun_inv_map`(counit で挟む)の**第 1 成分**を取ると
`𝟙` を同型で挟んだものになり、`𝟙` に戻る。 -/
theorem isBaseIdentity_inv_map_nComp (Z : D₁ × SingleObj ℕ+) (n : ℕ+) :
    IsBaseIdentity P₁ ((toProdCat P₁).inv.map ((𝟙 Z.1, n) : Z ⟶ Z)) := by
  have h1 := congrArg Prod.fst
    (Functor.fun_inv_map (toProdCat P₁) Z Z ((𝟙 Z.1, n) : Z ⟶ Z))
  have hc : ((toProdCat P₁).asEquivalence.counit.app Z).1 ≫
      ((toProdCat P₁).asEquivalence.counitInv.app Z).1
      = 𝟙 ((toProdCat P₁).obj ((toProdCat P₁).inv.obj Z)).1 :=
    congrArg Prod.fst ((toProdCat P₁).asEquivalence.counitIso.hom_inv_id_app Z)
  show P₁.Base _ = P₁.Base (𝟙 _)
  rw [P₁.Base_id]
  show ((toProdCat P₁).map ((toProdCat P₁).inv.map ((𝟙 Z.1, n) : Z ⟶ Z))).1 = 𝟙 _
  refine h1.trans (Eq.trans ?_ hc)
  exact congrArg (fun g => ((toProdCat P₁).asEquivalence.counit.app Z).1 ≫ g)
    (Category.id_comp _)

end BaseFunctor

end ABC3.Found.FrdI
