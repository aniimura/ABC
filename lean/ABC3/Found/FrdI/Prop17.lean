import ABC3.Found.FrdI.Prop14

/-!
# [FrdI] Proposition 1.7 —— Composites of Morphisms

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、
物理 p.28–p.30(**400 dpi 目視確認 2026-08-15**)。

原文 (FrdI p.28):
> Proposition 1.7. (Composites of Morphisms) Let Φ be a divisorial

原文 (FrdI p.28):
> monoid on a connected, totally epimorphic category D; C →FΦ a Frobenioid.

## ★この命題の規模(測定)

**5条 (i)–(v)**、主張の数は **25**:

| 条 | 主張の数 | 内容 |
|---|---|---|
| (i) | 9 | 合成で閉じる射のクラス 9 個 |
| (ii) | 3 | pull-back ⟺ minimal-adjoint / base-iso ⟺ minimal-coadjoint / base-iso の分解 |
| (iii) | 3 | Frobenius 型 ⟺ minimal-coadjoint / linear ⟺ minimal-adjoint / linear の分解 |
| (iv) | 1 | pre-step が co-angular ⟺ isometric pre-step に mid-adjoint |
| (v) | 9 | 合成が X なら両因子も X(8 クラス ＋ 始域が isotropic なときの Frobenius 型) |

★原文の証明は **3 段落**で、依存の向きは
`(i) → (ii)(iii) の十分性 → (ii)(iii) の必要性 → (v) → (iv) → (v) の残り`。
★**(iv) は (v) を使い、(v) の後半は (iv) を使う**ので、順序を守る必要がある。

## ★`Definition 1.3` のどの条項を引くか(段ごとの記録)

* **(i)** co-angular: `Definition 1.3, (iii), (a)`(= `FrobenioidCore.coAngularComp`)
* **(ii)(iii) 十分性**: `Definition 1.3, (iv), (a)` の**分解の存在**(`arbFactor`)
* **(ii)(iii) 必要性**: 同じ分解の**本質的一意性**(`arbFactorUniq`)と `𝒞` の totally epimorphic 性
* **(v) 同型**: `𝒞` が totally epimorphic
* **(v) base-iso**: `𝒟` が totally epimorphic
* **(v) linear**: `ℕ≥1` の乗法モノイドの構造(単元が 1 だけ)
* **(v) isometry**: `Φ` が **sharp**(`Definition 1.1, (i)`)＋ `Definition 1.1, (ii), (a)` の
  characteristic injectivity
-/

namespace ABC3.Found.FrdI

open CategoryTheory Opposite

universe v u w u2 v2

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ)

/-! ## ★(i) —— 合成で閉じる 9 クラス

原文 (FrdI p.28):
> (i) The following classes of morphisms are closed under composition: isome-

原文 (FrdI p.28):
> tries, base-isomorphisms, base-FSM-morphisms, pull-back morphisms,

原文 (FrdI p.28):
> linear morphisms, pre-steps, co-angular morphisms, LB-invertible mor-

★9 クラスのうち **5 つは既に持っていた**(`IsIsometric.comp` / `isBaseIsomorphism_comp` /
`IsPullBack.comp` / `IsLinear.comp` / `IsPreStep.comp`)。
残り 4 つをここで埋める。
-/

/-- **(i)** `base-FSM-morphism` は合成で閉じる。

★§0 の `IsFSMMorphism.comp`(原文 p.14 の「One verifies immediately」を
我々が証明したもの)を `Base` に沿って運ぶだけ。 -/
theorem IsBaseFSM.comp {A B E : C} {ψ : A ⟶ B} {φ : B ⟶ E}
    (hψ : IsBaseFSM P ψ) (hφ : IsBaseFSM P φ) : IsBaseFSM P (ψ ≫ φ) := by
  show IsFSMMorphism (P.Base (ψ ≫ φ))
  rw [P.Base_comp]
  exact IsFSMMorphism.comp hψ hφ

/-- **(i)** `LB-invertible` は合成で閉じる。

★co-angular の側は `Definition 1.3, (iii), (a)` を引く。**原文が明示する唯一の依存**。 -/
theorem IsLBInvertible.comp (F : FrobenioidCore P) {A B E : C} {ψ : A ⟶ B} {φ : B ⟶ E}
    (hψ : IsLBInvertible P ψ) (hφ : IsLBInvertible P φ) : IsLBInvertible P (ψ ≫ φ) :=
  ⟨F.coAngularComp ψ φ hψ.1 hφ.1, hψ.2.comp P hφ.2⟩

/-- **(i)** `morphism of Frobenius type` は合成で閉じる。 -/
theorem IsFrobeniusType.comp (F : FrobenioidCore P) {A B E : C} {ψ : A ⟶ B} {φ : B ⟶ E}
    (hψ : IsFrobeniusType P ψ) (hφ : IsFrobeniusType P φ) : IsFrobeniusType P (ψ ≫ φ) :=
  ⟨IsLBInvertible.comp P F hψ.1 hφ.1, isBaseIsomorphism_comp P hψ.2 hφ.2⟩

/-- **(i)** `co-angular` は合成で閉じる —— `Definition 1.3, (iii), (a)` そのもの。 -/
theorem IsCoAngular.comp (F : FrobenioidCore P) {A B E : C} {ψ : A ⟶ B} {φ : B ⟶ E}
    (hψ : IsCoAngular P ψ) (hφ : IsCoAngular P φ) : IsCoAngular P (ψ ≫ φ) :=
  F.coAngularComp ψ φ hψ hφ

/-! ## ★(v) —— 合成が X なら両因子も X

原文 (FrdI p.29):
> (v) If a composite morphism φ = α◦β of C is a(n) isomorphism (respectively,

原文 (FrdI p.29):
> base-isomorphism; linear morphism; pre-step; isometry; co-angular pre-

★原文の `φ = α ◦ β` は「先に `β`」なので、Lean では `φ = β ≫ α`。
★**(ii)(iii)(iv) より先に置く** —— (iv) が (v) を使うからである。
-/

/-- **(v)** 同型 —— `𝒞` が totally epimorphic であることから。

★これは §0 で既に証明した `isIso_of_isIso_comp`(原文 p.16 の「in passing」)である。
★`P` は `totEpiC` を取り出すためだけに要る(主張の型には現れない)ので、
`include` で明示的に引き入れる。 -/
theorem prop_1_7_v_isIso (hC : IsTotallyEpimorphic C) {A X B : C}
    (β : A ⟶ X) (α : X ⟶ B) [IsIso (β ≫ α)] :
    IsIso β ∧ IsIso α :=
  isIso_of_isIso_comp hC β α

/-- **(v)** base-isomorphism —— `𝒟` が totally epimorphic であることから。

★同じ補題を **`𝒟` の側で**使う。★(v) の 2 条が「同じ補題を別の圏で使う」形で
並んでいることは、原文の括弧書き
「from the fact that C is totally epimorphic (respectively, from the fact that
D is totally epimorphic; …)」がそのまま述べている。 -/
theorem prop_1_7_v_baseIso {A X B : C} (β : A ⟶ X) (α : X ⟶ B)
    (h : IsBaseIsomorphism P (β ≫ α)) :
    IsBaseIsomorphism P β ∧ IsBaseIsomorphism P α := by
  haveI : IsIso (P.Base β ≫ P.Base α) := by rw [← P.Base_comp]; exact h
  exact isIso_of_isIso_comp P.totEpiD (P.Base β) (P.Base α)

/-- **(v)** linear —— `ℕ≥1` の乗法モノイドの単元が `1` だけであることから。 -/
theorem prop_1_7_v_linear {A X B : C} (β : A ⟶ X) (α : X ⟶ B)
    (h : IsLinear P (β ≫ α)) : IsLinear P β ∧ IsLinear P α := by
  have hh : P.degFr α * P.degFr β = 1 := by rw [← P.degFr_comp]; exact h
  exact ⟨pnat_right_eq_one hh, pnat_left_eq_one hh⟩

/-- **(v)** pre-step —— base-isomorphism と linear の合わせ技。 -/
theorem prop_1_7_v_preStep {A X B : C} (β : A ⟶ X) (α : X ⟶ B)
    (h : IsPreStep P (β ≫ α)) : IsPreStep P β ∧ IsPreStep P α :=
  ⟨⟨(prop_1_7_v_linear P β α h.1).1, (prop_1_7_v_baseIso P β α h.2).1⟩,
   ⟨(prop_1_7_v_linear P β α h.1).2, (prop_1_7_v_baseIso P β α h.2).2⟩⟩

/-- ★sharp なモノイドで `x + y = 0` なら両方 `0`。 -/
private theorem eq_zero_of_add_eq_zero {M : Type w} [AddCommMonoid M] (h : IsSharp M)
    {x y : M} (hxy : x + y = 0) : x = 0 ∧ y = 0 :=
  ⟨h x ⟨⟨x, y, hxy, by rw [add_comm]; exact hxy⟩, rfl⟩,
   h y ⟨⟨y, x, by rw [add_comm]; exact hxy, hxy⟩, rfl⟩⟩

/-- **(v)** isometry。

★ここだけ **`Φ` が sharp**(= `divisorial`、`Definition 1.1, (i)`)を使う。
`Div(β ≫ α) = Φ(Base β)(Div α) + deg_Fr(α)·Div β = 0` から、sharp なので両項が `0`。

* 第1項が `0` ⟹ `Div α = 0` —— **`Definition 1.1, (ii), (a)` の単射性**
* 第2項が `0` ⟹ `Div β = 0` —— sharp から出る捩れの無さ(`isTorsionFreeNaive_of_isSharp`)

★`Proposition 1.5` は `Φ` が pre-divisorial でよかったが、
**ここは sharp が要る**。原文が 1.7 で `divisorial` に強めているのはこのためである。 -/
theorem prop_1_7_v_isometric {A X B : C} (β : A ⟶ X) (α : X ⟶ B)
    (h : IsIsometric P (β ≫ α)) : IsIsometric P β ∧ IsIsometric P α := by
  have hsharp : IsSharp (Φ.val (P.toElem.obj A).base) := (P.divisorial _).2
  have hd : Φ.map (P.Base β) (P.Div α) + ((P.degFr α : ℕ+) : ℕ) • P.Div β = 0 := by
    rw [← P.Div_comp]; exact h
  obtain ⟨h1, h2⟩ := eq_zero_of_add_eq_zero hsharp hd
  refine ⟨?_, ?_⟩
  · exact isTorsionFreeNaive_of_isSharp hsharp _ _ (P.degFr α).pos h2
  · refine Φ.map_injective (P.Base β) ?_
    rw [h1, map_zero]

/-! ### ★「同型ならそのクラスに属する」

★(ii)(iii) の必要性の議論で、`arbFactor` の因子が同型だと分かったあと
「では合成も元のクラスに入る」と言うために、毎回これが要る。
★**`isometric` と `co-angular` は自明ではない** —— 前者は (v) の isometry を、
後者は `𝒞` の totally epimorphic 性を使う。 -/

theorem isLinear_of_isIso {A B : C} (φ : A ⟶ B) [IsIso φ] : IsLinear P φ := by
  have h : P.degFr (inv φ) * P.degFr φ = 1 := by
    rw [← P.degFr_comp, IsIso.hom_inv_id, P.degFr_id]
  exact pnat_right_eq_one h

theorem isIsometric_of_isIso {A B : C} (φ : A ⟶ B) [IsIso φ] : IsIsometric P φ :=
  (prop_1_7_v_isometric P φ (inv φ) (by rw [IsIso.hom_inv_id]; exact P.Div_id A)).1

/-- ★**同型は co-angular** —— 3分解の真ん中は `isIso_of_isIso_comp` を2回で同型になる。 -/
theorem isCoAngular_of_isIso {A B : C} (φ : A ⟶ B) [IsIso φ] : IsCoAngular P φ := by
  rintro X Y γ β α rfl - - - -
  haveI : IsIso (β ≫ α) := (isIso_of_isIso_comp P.totEpiC γ (β ≫ α)).2
  exact (isIso_of_isIso_comp P.totEpiC β α).1

theorem isPreStep_of_isIso {A B : C} (φ : A ⟶ B) [IsIso φ] : IsPreStep P φ :=
  ⟨isLinear_of_isIso P φ, isBaseIsomorphism_of_isIso P φ⟩

theorem isLBInvertible_of_isIso {A B : C} (φ : A ⟶ B) [IsIso φ] : IsLBInvertible P φ :=
  ⟨isCoAngular_of_isIso P φ, isIsometric_of_isIso P φ⟩

theorem isFrobeniusType_of_isIso {A B : C} (φ : A ⟶ B) [IsIso φ] : IsFrobeniusType P φ :=
  ⟨isLBInvertible_of_isIso P φ, isBaseIsomorphism_of_isIso P φ⟩

/-! ## ★(iii) と (ii) の「分解」の形

原文 (FrdI p.29):
> a morphism of C is linear if and only if it is may be written as a composite α◦β,

原文 (FrdI p.29):
> where α is a pull-back morphism, and β is a pre-step.

★どちらも `Definition 1.3, (iv), (a)` の3分解 `φ = γ ≫ β ≫ α` から出る ——
**片方の端の因子が同型になる**ことを示せばよい。
-/

/-- **(iii)** `linear ⟺ pre-step ≫ pull-back` と書ける。

★`⇒` は3分解の **Frobenius 型因子 `γ` が同型になる**ことによる:
`φ` が linear なら (v) で `γ` も linear、Frobenius 型 ＋ linear は
**LB-invertible な pre-step** なので `Proposition 1.4, (iii)` で同型。 -/
theorem prop_1_7_iii_linear_factor (F : FrobenioidCore P) {A B : C} (φ : A ⟶ B) :
    IsLinear P φ ↔
      ∃ (X : C) (β : A ⟶ X) (α : X ⟶ B),
        φ = β ≫ α ∧ IsPreStep P β ∧ IsPullBack P α := by
  constructor
  · intro hlin
    obtain ⟨X, Y, γ, β, α, rfl, hγ, hβ, hα⟩ := F.arbFactor φ
    have hγlin : IsLinear P γ := (prop_1_7_v_linear P γ (β ≫ α) hlin).1
    haveI hγiso : IsIso γ :=
      prop_1_4_iii P F γ hγ.1 ⟨hγlin, hγ.2⟩
    exact ⟨Y, γ ≫ β, α, by rw [Category.assoc], IsPreStep.comp P (isPreStep_of_isIso P γ) hβ, hα⟩
  · rintro ⟨X, β, α, rfl, hβ, hα⟩
    exact IsLinear.comp P hβ.1 ((F.pullBackLB α hα).2)

/-- **(ii)** `base-isomorphism ⟺ Frobenius 型 ≫ pre-step` と書ける。

原文 (FrdI p.29):
> as a composite α ◦β, where α is a pre-step, and β is a morphism of Frobenius

★`⇒` は3分解の **pull-back 因子 `α` が同型になる**ことによる:
`φ` が base-iso なら (v) で `α` も base-iso、pull-back ＋ base-iso は
**LB-invertible な pre-step** なので `Proposition 1.4, (iii)` で同型。 -/
theorem prop_1_7_ii_baseIso_factor (F : FrobenioidCore P) {A B : C} (φ : A ⟶ B) :
    IsBaseIsomorphism P φ ↔
      ∃ (X : C) (γ : A ⟶ X) (β : X ⟶ B),
        φ = γ ≫ β ∧ IsFrobeniusType P γ ∧ IsPreStep P β := by
  constructor
  · intro hbi
    obtain ⟨X, Y, γ, β, α, rfl, hγ, hβ, hα⟩ := F.arbFactor φ
    have hαbi : IsBaseIsomorphism P α :=
      (prop_1_7_v_baseIso P β α (prop_1_7_v_baseIso P γ (β ≫ α) hbi).2).2
    haveI hαiso : IsIso α :=
      prop_1_4_iii P F α (F.pullBackLB α hα).1 ⟨(F.pullBackLB α hα).2, hαbi⟩
    exact ⟨X, γ, β ≫ α, rfl, hγ, IsPreStep.comp P hβ (isPreStep_of_isIso P α)⟩
  · rintro ⟨X, γ, β, rfl, hγ, hβ⟩
    exact isBaseIsomorphism_comp P hγ.2 hβ.2

/-! ### ★射のクラスを `MorphismProperty` として名付ける

`minimal-adjoint` などが `MorphismProperty` を取るので、各クラスを名付ける。
(`pullBackProp` は `Frobenioid.lean` に既にある。) -/

/-- `base-isomorphism` のクラス。 -/
def baseIsoProp : MorphismProperty C := fun _ _ φ => IsBaseIsomorphism P φ

/-- `linear morphism` のクラス。 -/
def linearProp : MorphismProperty C := fun _ _ φ => IsLinear P φ

/-- `morphism of Frobenius type` のクラス。 -/
def frobTypeProp : MorphismProperty C := fun _ _ φ => IsFrobeniusType P φ

/-- `isometric pre-step` のクラス。 -/
def isometricPreStepProp : MorphismProperty C :=
  fun _ _ φ => IsIsometric P φ ∧ IsPreStep P φ

/-! ### ★★co-angularity を「自明な因子」に当てる技

★(ii) の pull-back の必要性と (iii) の Frobenius 型の必要性は、
どちらも **`γ` か `α` を `𝟙` に取って co-angularity をそのまま使う**ことで出る。
★原文はここを「the essential uniqueness of this factorization」で説明するが、
**分解の一意性を使わずに済む**。★短い道があることを測定として記録する。 -/

/-- `φ = β ≫ α`、`β` が isometric pre-step、`α` が linear なら、
`φ` の co-angularity から `β` は同型(`γ := 𝟙`)。 -/
theorem isIso_of_isCoAngular_left {A X B : C} (β : A ⟶ X) (α : X ⟶ B)
    (hco : IsCoAngular P (β ≫ α)) (hβi : IsIsometric P β) (hβs : IsPreStep P β)
    (hα : IsLinear P α) : IsIso β :=
  hco A X (𝟙 A) β α (Category.id_comp _).symm hα hβi hβs
    (Or.inr (isBaseIsomorphism_of_isIso P (𝟙 A)))

/-- `φ = γ ≫ β`、`β` が isometric pre-step なら、
`φ` の co-angularity から `β` は同型(`α := 𝟙`)。 -/
theorem isIso_of_isCoAngular_right {A X B : C} (γ : A ⟶ X) (β : X ⟶ B)
    (hco : IsCoAngular P (γ ≫ β)) (hβi : IsIsometric P β) (hβs : IsPreStep P β) :
    IsIso β :=
  hco X B γ β (𝟙 B) (by rw [Category.comp_id]) (isLinear_id P B) hβi hβs
    (Or.inl (isBaseIsomorphism_of_isIso P (𝟙 B)))

/-! ## ★(ii)(iii) の adjoint 特徴づけ

原文 (FrdI p.28):
> (ii) A morphism of C is a pull-back morphism if and only if it is minimal-

原文 (FrdI p.28):
> adjoint to the base-isomorphisms of C. A morphism of C is a base-isomorphism

原文 (FrdI p.29):
> (iii) A morphism of C is of Frobenius type if and only if it is minimal-

原文 (FrdI p.29):
> coadjoint to the linear morphisms of C. A morphism of C is linear if and only if
-/

/-- **(iii)** `linear ⟺ Frobenius 型に minimal-adjoint`。

★`⇒` は `Proposition 1.4, (iii)`(LB-invertible な pre-step は同型)だけで済む。
★原文は「the well-known structure of the multiplicative monoid `N≥1` and the
essential uniqueness of morphisms of Frobenius type of a given Frobenius degree」
という別の道を挙げるが、**こちらのほうが短い**。 -/
theorem prop_1_7_iii_linear_minimalAdjoint (F : FrobenioidCore P) {A B : C} (φ : A ⟶ B) :
    IsLinear P φ ↔ IsMinimalAdjoint (frobTypeProp P) φ := by
  constructor
  · rintro hlin X β α rfl hβ
    exact prop_1_4_iii P F β hβ.1 ⟨(prop_1_7_v_linear P β α hlin).1, hβ.2⟩
  · intro hm
    obtain ⟨X, Y, γ, β, α, rfl, hγ, hβ, hα⟩ := F.arbFactor φ
    haveI : IsIso γ := hm X γ (β ≫ α) rfl hγ
    exact IsLinear.comp P (isLinear_of_isIso P γ)
      (IsLinear.comp P hβ.1 (F.pullBackLB α hα).2)

/-- **(iii)** `Frobenius 型 ⟺ linear に minimal-coadjoint`。

★`⇒` は **co-angularity を `α := 𝟙` に当てる技**で出る。 -/
theorem prop_1_7_iii_frobType_minimalCoadjoint (F : FrobenioidCore P) {A B : C}
    (φ : A ⟶ B) : IsFrobeniusType P φ ↔ IsMinimalCoadjoint (linearProp P) φ := by
  constructor
  · rintro hft X α β rfl hβ
    refine isIso_of_isCoAngular_right P α β hft.1.1 ?_ ⟨hβ, ?_⟩
    · exact (prop_1_7_v_isometric P α β hft.1.2).2
    · exact (prop_1_7_v_baseIso P α β hft.2).2
  · intro hm
    obtain ⟨X, Y, γ, β, α, rfl, hγ, hβ, hα⟩ := F.arbFactor φ
    haveI : IsIso (β ≫ α) :=
      hm X γ (β ≫ α) rfl (IsLinear.comp P hβ.1 (F.pullBackLB α hα).2)
    exact IsFrobeniusType.comp P F hγ (isFrobeniusType_of_isIso P (β ≫ α))

/-- **(ii)** `pull-back ⟺ base-isomorphism に minimal-adjoint`。

★`⇒` は **co-angularity を `γ := 𝟙` に当てる技**で出る。 -/
theorem prop_1_7_ii_pullBack_minimalAdjoint (F : FrobenioidCore P) {A B : C} (φ : A ⟶ B) :
    IsPullBack P φ ↔ IsMinimalAdjoint (baseIsoProp P) φ := by
  constructor
  · rintro hpb X β α rfl hβ
    obtain ⟨hlb, hlin⟩ := (prop_1_4_ii P F _).mp hpb
    refine isIso_of_isCoAngular_left P β α hlb.1 ?_ ⟨?_, hβ⟩ ?_
    · exact (prop_1_7_v_isometric P β α hlb.2).1
    · exact (prop_1_7_v_linear P β α hlin).1
    · exact (prop_1_7_v_linear P β α hlin).2
  · intro hm
    obtain ⟨X, Y, γ, β, α, rfl, hγ, hβ, hα⟩ := F.arbFactor φ
    haveI : IsIso (γ ≫ β) :=
      hm Y (γ ≫ β) α (Category.assoc _ _ _).symm (isBaseIsomorphism_comp P hγ.2 hβ.2)
    rw [← Category.assoc]
    exact IsPullBack.comp P (isPullBack_of_isIso P (γ ≫ β)) hα

/-- **(ii)** `base-isomorphism ⟺ pull-back に minimal-coadjoint`。 -/
theorem prop_1_7_ii_baseIso_minimalCoadjoint (F : FrobenioidCore P) {A B : C} (φ : A ⟶ B) :
    IsBaseIsomorphism P φ ↔ IsMinimalCoadjoint (pullBackProp P) φ := by
  constructor
  · rintro hbi X α β rfl hβ
    exact prop_1_4_iii P F β (F.pullBackLB β hβ).1
      ⟨(F.pullBackLB β hβ).2, (prop_1_7_v_baseIso P α β hbi).2⟩
  · intro hm
    obtain ⟨X, Y, γ, β, α, rfl, hγ, hβ, hα⟩ := F.arbFactor φ
    haveI : IsIso α := hm Y (γ ≫ β) α (Category.assoc _ _ _).symm hα
    exact isBaseIsomorphism_comp P hγ.2
      (isBaseIsomorphism_comp P hβ.2 (isBaseIsomorphism_of_isIso P α))

/-! ## ★(iv) —— pre-step が co-angular ⟺ isometric pre-step に mid-adjoint

原文 (FrdI p.29):
> (iv) A pre-step of C is co-angular if and only if it is mid-adjoint [cf. §0] to

★`⇒` が中身である —— co-angular の定義には「`α` が linear」「`α` か `γ` が base-iso」
という**余分な仮定**が付いているが、`φ` が pre-step ならこの2つは (v) から**自動で従う**。
★原文が「cf. the argument applied in the proof of Proposition 1.4, (iv)!」と
感嘆符つきで指す議論はこれである。 -/
theorem prop_1_7_iv {A B : C} (φ : A ⟶ B) (hφ : IsPreStep P φ) :
    IsCoAngular P φ ↔ IsMidAdjoint (isometricPreStepProp P) φ := by
  constructor
  · rintro hco X Y γ β α rfl ⟨hβi, hβs⟩
    obtain ⟨-, hrest⟩ := prop_1_7_v_preStep P γ (β ≫ α) hφ
    obtain ⟨-, hα⟩ := prop_1_7_v_preStep P β α hrest
    exact hco X Y γ β α rfl hα.1 hβi hβs
      (Or.inr (prop_1_7_v_baseIso P γ (β ≫ α) hφ.2).1)
  · rintro hmid X Y γ β α rfl hα hβi hβs
    exact fun _ => hmid X Y γ β α rfl ⟨hβi, hβs⟩

/-! ## ★(v) の後半 —— (iv) を使う部分

原文 (FrdI p.29):
> assertion (v) for co-angular pre-steps follows

原文 (FrdI p.29):
> from assertion (v) for pre-steps and assertion (iv). To prove assertion (v) for co-
-/

/-- **(v)** co-angular pre-step。

★(iv) で「co-angular」を「mid-adjoint」に**翻訳してから**分配する。
mid-adjoint は「真ん中の因子」についての条件なので、
`β` の3分解を `φ` の3分解に**そのまま埋め込める**。 -/
theorem prop_1_7_v_coAngularPreStep {A X B : C} (β : A ⟶ X) (α : X ⟶ B)
    (hco : IsCoAngular P (β ≫ α)) (hs : IsPreStep P (β ≫ α)) :
    (IsCoAngular P β ∧ IsPreStep P β) ∧ (IsCoAngular P α ∧ IsPreStep P α) := by
  obtain ⟨hβs, hαs⟩ := prop_1_7_v_preStep P β α hs
  have hmid := (prop_1_7_iv P (β ≫ α) hs).mp hco
  refine ⟨⟨?_, hβs⟩, ?_, hαs⟩
  · refine (prop_1_7_iv P β hβs).mpr ?_
    rintro X' Y' γ' β' α' hfac hβ'
    refine hmid X' Y' γ' β' (α' ≫ α) ?_ hβ'
    rw [hfac]; simp
  · refine (prop_1_7_iv P α hαs).mpr ?_
    rintro X' Y' γ' β' α' hfac hβ'
    refine hmid X' Y' (β ≫ γ') β' α' ?_ hβ'
    rw [hfac]; simp

/-- ★linear 射を「isometric pre-step ≫ co-angular linear」に分解する。

★原文が (v) の co-angular linear の証明で使う分解であり、
**(iii) の分解と `Definition 1.3, (v), (c)` の合成**である:

1. (iii): `φ = p ≫ q`(`p` pre-step、`q` pull-back)
2. `Definition 1.3, (v), (c)`: `p = p₁ ≫ p₂`(`p₁` isometric pre-step、`p₂` co-angular pre-step)
3. `Proposition 1.4, (ii)`: `q` は pull-back なので co-angular linear

したがって `φ = p₁ ≫ (p₂ ≫ q)`。 -/
theorem linear_factor_isoPreStep_coAngLinear (F : FrobenioidCore P) {A B : C}
    (φ : A ⟶ B) (hlin : IsLinear P φ) :
    ∃ (X : C) (β : A ⟶ X) (α : X ⟶ B),
      φ = β ≫ α ∧ (IsIsometric P β ∧ IsPreStep P β) ∧
        (IsCoAngular P α ∧ IsLinear P α) := by
  obtain ⟨X, p, q, rfl, hp, hq⟩ := (prop_1_7_iii_linear_factor P F φ).mp hlin
  obtain ⟨Y, p₁, p₂, rfl, hp₁i, hp₁s, hp₂c, hp₂s⟩ := F.preStepFactor' p hp
  obtain ⟨hqlb, hqlin⟩ := (prop_1_4_ii P F q).mp hq
  exact ⟨Y, p₁, p₂ ≫ q, by rw [Category.assoc], ⟨hp₁i, hp₁s⟩,
    F.coAngularComp p₂ q hp₂c hqlb.1, IsLinear.comp P hp₂s.1 hqlin⟩

/-- **(v)** co-angular linear。★原文の証明でいちばん長い段。

原文 (FrdI p.29):
> angular linear morphisms, suppose that φ is co-angular and linear. Then observe

★段取り(原文どおり):

1. `α`, `β` をそれぞれ「isometric pre-step ≫ co-angular linear」に分解: `α = α₂ ≫ α₁`、`β = β₂ ≫ β₁`
2. 真ん中に現れる `β₁ ≫ α₂` を同じ形に分解: `β₁ ≫ α₂ = γ₂ ≫ γ₁`
3. すると `φ = (β₂ ≫ γ₂) ≫ (γ₁ ≫ α₁)` で、`β₂ ≫ γ₂` は isometric pre-step
   ⟹ `φ` の co-angularity で **`β₂ ≫ γ₂` は同型**
4. (v) の同型で `β₂`, `γ₂` が同型
5. `β₁ ≫ α₂ = γ₂ ≫ γ₁` は co-angular なので、そこへ co-angularity を当てて **`α₂` が同型**
6. `α = α₂ ≫ α₁`、`β = β₂ ≫ β₁` の左因子が同型なので、両方 co-angular -/
theorem prop_1_7_v_coAngularLinear (F : FrobenioidCore P) {A X B : C}
    (β : A ⟶ X) (α : X ⟶ B)
    (hco : IsCoAngular P (β ≫ α)) (hlin : IsLinear P (β ≫ α)) :
    (IsCoAngular P β ∧ IsLinear P β) ∧ (IsCoAngular P α ∧ IsLinear P α) := by
  obtain ⟨hβlin, hαlin⟩ := prop_1_7_v_linear P β α hlin
  obtain ⟨Xα, α₂, α₁, hαeq, ⟨hα₂i, hα₂s⟩, hα₁c, hα₁lin⟩ :=
    linear_factor_isoPreStep_coAngLinear P F α hαlin
  obtain ⟨Xβ, β₂, β₁, hβeq, ⟨hβ₂i, hβ₂s⟩, hβ₁c, hβ₁lin⟩ :=
    linear_factor_isoPreStep_coAngLinear P F β hβlin
  obtain ⟨Xγ, γ₂, γ₁, hγeq, ⟨hγ₂i, hγ₂s⟩, hγ₁c, hγ₁lin⟩ :=
    linear_factor_isoPreStep_coAngLinear P F (β₁ ≫ α₂)
      (IsLinear.comp P hβ₁lin hα₂s.1)
  -- 3. `φ = (β₂ ≫ γ₂) ≫ (γ₁ ≫ α₁)`
  have hsplit : (β ≫ α : A ⟶ B) = (β₂ ≫ γ₂) ≫ (γ₁ ≫ α₁) := by
    rw [hβeq, hαeq, Category.assoc, ← Category.assoc β₁, hγeq]
    simp
  have hmidiso : IsIso (β₂ ≫ γ₂) := by
    refine isIso_of_isCoAngular_left P (β₂ ≫ γ₂) (γ₁ ≫ α₁) (by rw [← hsplit]; exact hco)
      (IsIsometric.comp P hβ₂i hγ₂i) (IsPreStep.comp P hβ₂s hγ₂s)
      (IsLinear.comp P hγ₁lin hα₁lin)
  -- 4. `β₂`, `γ₂` が同型
  obtain ⟨hβ₂iso, hγ₂iso⟩ := isIso_of_isIso_comp P.totEpiC β₂ γ₂
  -- 5. `α₂` が同型
  haveI := hγ₂iso
  have hco' : IsCoAngular P (β₁ ≫ α₂) := by
    rw [hγeq]
    exact F.coAngularComp γ₂ γ₁ (isCoAngular_of_isIso P γ₂) hγ₁c
  haveI hα₂iso : IsIso α₂ := isIso_of_isCoAngular_right P β₁ α₂ hco' hα₂i hα₂s
  -- 6. まとめ
  haveI := hβ₂iso
  refine ⟨⟨?_, hβlin⟩, ?_, hαlin⟩
  · rw [hβeq]
    exact F.coAngularComp β₂ β₁ (isCoAngular_of_isIso P β₂) hβ₁c
  · rw [hαeq]
    exact F.coAngularComp α₂ α₁ (isCoAngular_of_isIso P α₂) hα₁c

/-- **(v)** pull-back —— co-angular linear と isometry の合わせ技。

原文 (FrdI p.30):
> as desired. Now assertion (v) for pull-back morphisms follows from assertion (v) for
-/
theorem prop_1_7_v_pullBack (F : FrobenioidCore P) {A X B : C} (β : A ⟶ X) (α : X ⟶ B)
    (h : IsPullBack P (β ≫ α)) : IsPullBack P β ∧ IsPullBack P α := by
  obtain ⟨hlb, hlin⟩ := (prop_1_4_ii P F _).mp h
  obtain ⟨⟨hβco, hβlin⟩, hαco, hαlin⟩ := prop_1_7_v_coAngularLinear P F β α hlb.1 hlin
  obtain ⟨hβi, hαi⟩ := prop_1_7_v_isometric P β α hlb.2
  exact ⟨(prop_1_4_ii P F β).mpr ⟨⟨hβco, hβi⟩, hβlin⟩,
    (prop_1_4_ii P F α).mpr ⟨⟨hαco, hαi⟩, hαlin⟩⟩

/-- **(v)** Frobenius 型(始域が isotropic のとき)。

原文 (FrdI p.30):
> for morphisms of Frobenius type in Cistr [cf. Definition 1.3, (vii), (b)] follows from

★`Definition 1.3, (vii), (b)`(isotropic は射の先へ伝わる)で
`X` も isotropic になり、`Proposition 1.4, (i)` により
**`A` からも `X` からも出る射はすべて co-angular**。
残るのは isometric と base-isomorphism で、どちらも既に (v) にある。 -/
theorem prop_1_7_v_frobeniusType (F : FrobenioidCore P) {A X B : C} (β : A ⟶ X) (α : X ⟶ B)
    (hA : IsIsotropic P A) (h : IsFrobeniusType P (β ≫ α)) :
    IsFrobeniusType P β ∧ IsFrobeniusType P α := by
  obtain ⟨hβb, hαb⟩ := prop_1_7_v_baseIso P β α h.2
  obtain ⟨hβi, hαi⟩ := prop_1_7_v_isometric P β α h.1.2
  refine ⟨(prop_1_4_i_frobeniusType P β (fun Y f => F.isotropicClosed f hA)).mpr ⟨hβi, hβb⟩,
    (prop_1_4_i_frobeniusType P α (fun Y f => F.isotropicClosed f ?_)).mpr ⟨hαi, hαb⟩⟩
  exact F.isotropicClosed β hA

/-! ## ★負の対照 —— 原文が (i) の列挙から**外している**クラス

★原文の列挙は 9 個だが、我々が持っている射のクラスはもっと多い。
**外れているものが本当に閉じないのか**を確かめる。

| 外れているクラス | 閉じるか |
|---|---|
| `step`(同型でない pre-step) | ★**閉じる**(下の `IsStep.comp`) |
| `prime-Frobenius morphism` | ★**閉じない**(下の `not_isPrimeFrobenius_wDeg2_comp`) |

★**結論の2つの向き**:
* 原文の列挙は「合成で閉じるクラスの**全部**」ではない(`step` が漏れている)。
* しかし `prime-Frobenius` を外しているのは**正しい**。
-/

/-- ★**`step` は合成で閉じる** —— 原文の列挙には無いが、実際には閉じる。

`pre-step` であることは `IsPreStep.comp`、
同型でないことは **(v) の同型**(合成が同型なら両因子も同型)の対偶。 -/
theorem IsStep.comp {A B E : C} {ψ : A ⟶ B} {φ : B ⟶ E}
    (hψ : IsStep P ψ) (hφ : IsStep P φ) : IsStep P (ψ ≫ φ) := by
  refine ⟨IsPreStep.comp P hψ.1 hφ.1, fun hiso => ?_⟩
  haveI := hiso
  exact hψ.2 (isIso_of_isIso_comp P.totEpiC ψ φ).1

/-- ★模型 `wP` の `wDeg2 = ⟨𝟙, 0, 2⟩` は **prime-Frobenius morphism**。 -/
theorem wIsPrimeFrobenius_wDeg2 : IsPrimeFrobenius wP wDeg2 := by
  refine ⟨⟨⟨wIsCoAngular _, rfl⟩, ?_⟩, Nat.prime_two⟩
  show IsIso (𝟙 wTop.base)
  infer_instance

/-- ★★**prime-Frobenius morphism は合成で閉じない** —— 次数が `2 * 2 = 4` になる。

★したがって原文が (i) の列挙から `prime-Frobenius` を外しているのは**正しい**。
★対照的に `morphism of Frobenius type`(次数に条件を課さない)は
`IsFrobeniusType.comp` のとおり閉じる。**外した理由は「素数」の部分にある。** -/
theorem not_isPrimeFrobenius_wDeg2_comp :
    ¬ IsPrimeFrobenius wP (wDeg2 ≫ wDeg2) := by
  rintro ⟨-, hp⟩
  have h4 : ((wP.degFr (wDeg2 ≫ wDeg2) : ℕ+) : ℕ) = 4 := rfl
  rw [h4] at hp
  exact absurd hp (by decide)

/-- ★対照: 同じ `wDeg2 ≫ wDeg2` は **Frobenius 型としては閉じている**。 -/
theorem wIsFrobeniusType_wDeg2_comp : IsFrobeniusType wP (wDeg2 ≫ wDeg2) :=
  IsFrobeniusType.comp wP wIsFrobenioid.core wIsPrimeFrobenius_wDeg2.1
    wIsPrimeFrobenius_wDeg2.1

/-! ## ★(i) と (v) をそれぞれ1つの主張にまとめる

★`.src` を付ける単位を原文の条に揃えるため、9 クラスぶんを束ねた形も置く。 -/

/-- ★**(i) 全体** —— 原文が挙げる 9 クラスがすべて合成で閉じる。 -/
theorem prop_1_7_i (F : FrobenioidCore P) {A B E : C} (ψ : A ⟶ B) (φ : B ⟶ E) :
    (IsIsometric P ψ → IsIsometric P φ → IsIsometric P (ψ ≫ φ)) ∧
    (IsBaseIsomorphism P ψ → IsBaseIsomorphism P φ → IsBaseIsomorphism P (ψ ≫ φ)) ∧
    (IsBaseFSM P ψ → IsBaseFSM P φ → IsBaseFSM P (ψ ≫ φ)) ∧
    (IsPullBack P ψ → IsPullBack P φ → IsPullBack P (ψ ≫ φ)) ∧
    (IsLinear P ψ → IsLinear P φ → IsLinear P (ψ ≫ φ)) ∧
    (IsPreStep P ψ → IsPreStep P φ → IsPreStep P (ψ ≫ φ)) ∧
    (IsCoAngular P ψ → IsCoAngular P φ → IsCoAngular P (ψ ≫ φ)) ∧
    (IsLBInvertible P ψ → IsLBInvertible P φ → IsLBInvertible P (ψ ≫ φ)) ∧
    (IsFrobeniusType P ψ → IsFrobeniusType P φ → IsFrobeniusType P (ψ ≫ φ)) :=
  ⟨fun h1 h2 => IsIsometric.comp P h1 h2,
   fun h1 h2 => isBaseIsomorphism_comp P h1 h2,
   fun h1 h2 => IsBaseFSM.comp P h1 h2,
   fun h1 h2 => IsPullBack.comp P h1 h2,
   fun h1 h2 => IsLinear.comp P h1 h2,
   fun h1 h2 => IsPreStep.comp P h1 h2,
   fun h1 h2 => IsCoAngular.comp P F h1 h2,
   fun h1 h2 => IsLBInvertible.comp P F h1 h2,
   fun h1 h2 => IsFrobeniusType.comp P F h1 h2⟩

/-- ★**(v) 全体** —— 合成が属せば両因子も属する 8 クラス ＋ 始域が isotropic なときの
Frobenius 型。 -/
theorem prop_1_7_v (F : FrobenioidCore P) {A X B : C} (β : A ⟶ X) (α : X ⟶ B) :
    (IsIso (β ≫ α) → IsIso β ∧ IsIso α) ∧
    (IsBaseIsomorphism P (β ≫ α) → IsBaseIsomorphism P β ∧ IsBaseIsomorphism P α) ∧
    (IsLinear P (β ≫ α) → IsLinear P β ∧ IsLinear P α) ∧
    (IsPreStep P (β ≫ α) → IsPreStep P β ∧ IsPreStep P α) ∧
    (IsIsometric P (β ≫ α) → IsIsometric P β ∧ IsIsometric P α) ∧
    (IsCoAngular P (β ≫ α) → IsPreStep P (β ≫ α) →
      (IsCoAngular P β ∧ IsPreStep P β) ∧ (IsCoAngular P α ∧ IsPreStep P α)) ∧
    (IsCoAngular P (β ≫ α) → IsLinear P (β ≫ α) →
      (IsCoAngular P β ∧ IsLinear P β) ∧ (IsCoAngular P α ∧ IsLinear P α)) ∧
    (IsPullBack P (β ≫ α) → IsPullBack P β ∧ IsPullBack P α) ∧
    (IsIsotropic P A → IsFrobeniusType P (β ≫ α) →
      IsFrobeniusType P β ∧ IsFrobeniusType P α) :=
  ⟨fun h => by haveI := h; exact prop_1_7_v_isIso P.totEpiC β α,
   fun h => prop_1_7_v_baseIso P β α h,
   fun h => prop_1_7_v_linear P β α h,
   fun h => prop_1_7_v_preStep P β α h,
   fun h => prop_1_7_v_isometric P β α h,
   fun h1 h2 => prop_1_7_v_coAngularPreStep P β α h1 h2,
   fun h1 h2 => prop_1_7_v_coAngularLinear P F β α h1 h2,
   fun h => prop_1_7_v_pullBack P F β α h,
   fun h1 h2 => prop_1_7_v_frobeniusType P F β α h1 h2⟩

/-! ### ★出典の紐付け(`.src`) -/

def prop_1_7_i.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 28, item := "Proposition 1.7, (i)",
    sectionId := "frdi-prop-1-7-i" }

def prop_1_7_ii_pullBack_minimalAdjoint.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 28, item := "Proposition 1.7, (ii)",
    sectionId := "frdi-prop-1-7-ii" }

def prop_1_7_iii_linear_minimalAdjoint.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 29, item := "Proposition 1.7, (iii)",
    sectionId := "frdi-prop-1-7-iii" }

def prop_1_7_iv.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 29, item := "Proposition 1.7, (iv)",
    sectionId := "frdi-prop-1-7-iv" }

def prop_1_7_v.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 29, item := "Proposition 1.7, (v)",
    sectionId := "frdi-prop-1-7-v" }

end ABC3.Found.FrdI
