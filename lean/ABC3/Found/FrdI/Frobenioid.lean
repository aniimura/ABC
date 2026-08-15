import ABC3.Found.FrdI.MorphismTypes
import Mathlib.CategoryTheory.Widesubcategory
import Mathlib.CategoryTheory.Comma.Over.Basic

/-!
# [FrdI] Definition 1.3 —— Frobenioid 本体

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、
物理 p.24–p.25(**400 dpi 目視確認 2026-08-15、実装のための再読**)。

原文 (FrdI p.24):
> Definition 1.3. Let D, Φ, C →FΦ be as in Definition 1.2. Then we shall say

原文 (FrdI p.24):
> that the pre-Frobenioid C →FΦ [i.e., C equipped with this functor] is a Frobenioid

## ★添字の確定(目視でしか決まらない箇所)

§0(物理 p.13)の規約:

* `𝒞_A`(**後置**) = 対象が `B → A` の圏 = **スライス**(mathlib の `Over A`)
* `_A𝒞`(**前置**) = 対象が `A → B` の圏 = **コスライス**(mathlib の `Under A`)

これを Definition 1.3 に当てると:

* (i)(c) の `𝒞^pl-bk_A ≝ (𝒞^pl-bk)_A` は**後置** → **スライス**
* (iii)(d) は**両方**出る:
  `_A𝒞^coa-pre ≝ _A(𝒞^coa-pre)` は**前置** → **コスライス**、行き先は `Order(Φ(A))`
  `𝒞^coa-pre_A ≝ (𝒞^coa-pre)_A` は**後置** → **スライス**、行き先は `Order(Φ(A))^opp`

★`pdftotext` では前置の `_A` が落ちるか行頭に紛れるので、ここは目視でしか決まらない。

## ★mathlib の実測(2026-08-15)

| 必要なもの | mathlib | 判定 |
|---|---|---|
| 射のクラスが定める**広い**部分圏 | **`CategoryTheory.WideSubcategory (P) [IsMultiplicative P]`** | ★**使う** |
| スライス / コスライス | `CategoryTheory.Over` / `Under` | ★**使う** |
| 前順序を圏と見る | `Preorder.smallCategory` | ★**使う** |
| 圏同値 | `CategoryTheory.Functor.IsEquivalence` | ★**使う** |

★★**`MorphismProperty.overObj` に飛びつかなかった理由**を記す。
`overObj W : ObjectProperty (Over X)` は「構造射が `W` に属する `Over X` の対象」の
**充満**部分圏を作る。しかし `(𝒞^pl-bk)_A` は
「対象も射も pull-back morphism」であって、**充満ではない**
(`Over A` の射のうち、それ自身 pull-back morphism であるものだけを残す)。
★語が近いだけで別物である。正しいのは「先に `WideSubcategory` を作り、その `Over`」の順である。

## ★条件の順序の問題(原文には書かれていない)

`(iii)(d)` は `𝒞^coa-pre`(co-angular pre-step が定める部分圏)を使うが、
**それが圏になること**(合成で閉じること)は `(iii)(a)` そのものである。
すなわち **(iii)(d) は (iii)(a) を前提している**。

★本ファイルはこれを次のように扱う:

* `𝒞^pl-bk` の合成閉性は **Definition 1.3 の条件ではない**——
  `Definition 1.2, (ii)` の定義から従う。下の `IsPullBack.comp` で**証明する**。
  したがって (i)(c) は**圏同値としてそのまま**書ける。
* `𝒞^coa-pre` の合成閉性は (iii)(a) なので、`Frobenioid` 構造の**先の場に置き**、
  (iii)(d) はその場を使って**圏同値として**書く。
-/

namespace ABC3.Found.FrdI

open CategoryTheory Opposite

universe v u w

variable {D : Type u} [Category.{v} D] {C : Type u} [Category.{v} C]
  {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ)

/-! ### ★`Order(M)` —— §0 の `≤` が定める圏

原文 (FrdI p.12):
> observe that the relation “≤” on elements of M determines a category
-/

/-- **[FrdI] §0 `Order(M)`** —— `M` を `≤`(`a ≤ b ⟺ ∃ c, a + c = b`)で
前順序と見た型。`Preorder.smallCategory` により圏になる。

★型シノニムにするのは、`M` 自身に前順序を入れると他のインスタンスと衝突するためである。 -/
def OrderCat (M : Type w) [AddCommMonoid M] : Type w := M

instance (M : Type w) [AddCommMonoid M] : AddCommMonoid (OrderCat M) :=
  inferInstanceAs (AddCommMonoid M)

instance (M : Type w) [AddCommMonoid M] : Preorder (OrderCat M) where
  le a b := MLe (M := M) a b
  le_refl a := ⟨0, add_zero a⟩
  le_trans _ _ _ := fun ⟨x, hx⟩ ⟨y, hy⟩ => ⟨x + y, by rw [← add_assoc, hx, hy]⟩

/-- `M` の元を `Order(M)` の対象と見る。 -/
def toOrderCat {M : Type w} [AddCommMonoid M] (a : M) : OrderCat M := a

@[simp] theorem orderCat_le_iff {M : Type w} [AddCommMonoid M] {a b : M} :
    toOrderCat a ≤ toOrderCat b ↔ MLe a b := Iff.rfl

/-! ### ★`𝒞^pl-bk` —— pull-back morphism が定める広い部分圏

★これは `Definition 1.3` の条件ではなく、`Definition 1.2, (ii)` の定義から従う。
原文は (iv)(b) で「Every pull-back morphism of C is LB-invertible and linear」とだけ
述べ、合成閉性には触れていないが、**ファイバー積の定義から形式的に従う**。 -/

/-- ★**pull-back morphism は合成で閉じる**。

`φ : A ⟶ B`、`φ' : B ⟶ E` とする。単射性は `φ'` の単射性で
`f₁ ≫ φ = f₂ ≫ φ` を出し、次に `φ` の単射性を使う。全射性は逆順に、
`φ'` の全射性で `f' : X ⟶ B` を取り、`φ` の全射性で `f : X ⟶ A` を取る。 -/
theorem IsPullBack.comp {A B E : C} {φ : A ⟶ B} {φ' : B ⟶ E}
    (hφ : IsPullBack P φ) (hφ' : IsPullBack P φ') : IsPullBack P (φ ≫ φ') := by
  intro X
  constructor
  · intro f₁ f₂ h
    have h1 : f₁ ≫ φ ≫ φ' = f₂ ≫ φ ≫ φ' :=
      congrArg (fun p => (p : (X ⟶ E) × _).1) (congrArg Subtype.val h)
    have h2 : P.Base f₁ = P.Base f₂ :=
      congrArg (fun p => (p : (X ⟶ E) × _).2) (congrArg Subtype.val h)
    have hmid : f₁ ≫ φ = f₂ ≫ φ := by
      refine (hφ' X).1 (Subtype.ext ?_)
      have : (f₁ ≫ φ) ≫ φ' = (f₂ ≫ φ) ≫ φ' := by
        rw [Category.assoc, Category.assoc]; exact h1
      simp [this, P.Base_comp, h2]
    exact (hφ X).1 (Subtype.ext (by simp [hmid, h2]))
  · rintro ⟨⟨g, h'⟩, hgh⟩
    obtain ⟨f', hf'⟩ := (hφ' X).2 ⟨(g, h' ≫ P.Base φ), by
      simpa [Category.assoc, P.Base_comp] using hgh⟩
    have hf'1 : f' ≫ φ' = g := congrArg (fun p => (p : (X ⟶ E) × _).1) (congrArg Subtype.val hf')
    have hf'2 : P.Base f' = h' ≫ P.Base φ :=
      congrArg (fun p => (p : (X ⟶ E) × _).2) (congrArg Subtype.val hf')
    obtain ⟨f, hf⟩ := (hφ X).2 ⟨(f', h'), hf'2⟩
    have hf1 : f ≫ φ = f' := congrArg (fun p => (p : (X ⟶ B) × _).1) (congrArg Subtype.val hf)
    have hf2 : P.Base f = h' := congrArg (fun p => (p : (X ⟶ B) × _).2) (congrArg Subtype.val hf)
    exact ⟨f, Subtype.ext (by simp [← Category.assoc, hf1, hf'1, hf2])⟩

/-- pull-back morphism が定める射のクラス。 -/
def pullBackProp : MorphismProperty C := fun _ _ φ => IsPullBack P φ

instance : (pullBackProp P).ContainsIdentities :=
  ⟨fun A => isPullBack_of_isIso P (𝟙 A)⟩

instance : (pullBackProp P).IsStableUnderComposition :=
  ⟨fun _ _ hf hg => IsPullBack.comp P hf hg⟩

instance : MorphismProperty.IsMultiplicative (pullBackProp P) where

/-- **`𝒞^pl-bk`** —— pull-back morphism が定める広い部分圏。 -/
abbrev PlBk : Type u := WideSubcategory (pullBackProp P)

/-! ### ★`𝟙` は co-angular

`𝒞^coa-pre` が恒等射を含むことに要る。★`PreFrobenioid` が持つ
`totEpiC`(𝒞 が totally epimorphic)と、`CategoryVocabulary.lean` で証明した
`isIso_of_isIso_comp`(原文 p.16 の「in passing」の観察)から従う。 -/

/-- ★恒等射は co-angular。

`𝟙 = γ ≫ β ≫ α` なら、totally epimorphic な圏では各因子が同型になるので
特に `β` が同型。 -/
theorem isCoAngular_id (A : C) : IsCoAngular P (𝟙 A) := by
  intro X Y γ β α hfac _ _ _ _
  haveI : IsIso (γ ≫ β ≫ α) := hfac ▸ (inferInstance : IsIso (𝟙 A))
  have h1 := isIso_of_isIso_comp P.totEpiC γ (β ≫ α)
  haveI : IsIso (β ≫ α) := h1.2
  exact (isIso_of_isIso_comp P.totEpiC β α).1

/-- co-angular pre-step が定める射のクラス。 -/
def coaPreProp : MorphismProperty C := fun _ _ φ => IsCoAngular P φ ∧ IsPreStep P φ

instance : (coaPreProp P).ContainsIdentities :=
  ⟨fun A => ⟨isCoAngular_id P A, isPreStep_id P A⟩⟩

/-- base-isomorphism は合成で閉じる。 -/
theorem isBaseIsomorphism_comp {A B E : C} {ψ : A ⟶ B} {φ : B ⟶ E}
    (hψ : IsBaseIsomorphism P ψ) (hφ : IsBaseIsomorphism P φ) :
    IsBaseIsomorphism P (ψ ≫ φ) := by
  show IsIso (P.Base (ψ ≫ φ))
  haveI h1 : IsIso (P.Base ψ) := hψ
  haveI h2 : IsIso (P.Base φ) := hφ
  rw [P.Base_comp]
  infer_instance

/-- pre-step は合成で閉じる。 -/
theorem IsPreStep.comp {A B E : C} {ψ : A ⟶ B} {φ : B ⟶ E}
    (hψ : IsPreStep P ψ) (hφ : IsPreStep P φ) : IsPreStep P (ψ ≫ φ) :=
  ⟨IsLinear.comp P hψ.1 hφ.1, isBaseIsomorphism_comp P hψ.2 hφ.2⟩

/-- ★**(iii)(a) から `𝒞^coa-pre` が乗法的であることが出る。**

★これが「(iii)(d) は (iii)(a) を前提している」の中身である。恒等射の側は
`isCoAngular_id`(totally epimorphic から)と `isPreStep_id` で埋まる。 -/
theorem coaPreProp_isMultiplicative
    (hca : ∀ {A B E : C} (ψ : A ⟶ B) (φ : B ⟶ E), IsCoAngular P ψ → IsCoAngular P φ →
      IsCoAngular P (ψ ≫ φ)) :
    MorphismProperty.IsMultiplicative (coaPreProp P) where
  id_mem A := ⟨isCoAngular_id P A, isPreStep_id P A⟩
  comp_mem f g hf hg := ⟨hca f g hf.1 hg.1, IsPreStep.comp P hf.2 hg.2⟩

/-! ### ★(i)(c) の関手 `(𝒞^pl-bk)_A → 𝒟_{A_𝒟}`

原文 (FrdI p.24):
> [where A D def = Base(A)] determined by the natural projection functor C →D is an
-/

/-- **(i)(c)** の関手。対象 `(B → A)`(pull-back morphism)を `(B_𝒟 → A_𝒟)` に送る。 -/
def plBkOverFunctor (A : C) :
    Over (⟨A⟩ : PlBk P) ⥤ Over ((P.toElem.obj A).base) where
  obj Z := Over.mk (P.Base Z.hom.hom)
  map {Z W} f := Over.homMk (P.Base f.left.hom) (by
    have h : f.left.hom ≫ W.hom.hom = Z.hom.hom :=
      congrArg InducedWideCategory.Hom.hom (Over.w f)
    show P.Base f.left.hom ≫ P.Base W.hom.hom = P.Base Z.hom.hom
    rw [← P.Base_comp, h])
  map_id Z := by
    apply Over.OverMorphism.ext
    exact P.Base_id Z.left.obj
  map_comp {Z W V} f g := by
    apply Over.OverMorphism.ext
    show P.Base (f.left.hom ≫ g.left.hom) = P.Base f.left.hom ≫ P.Base g.left.hom
    rw [P.Base_comp]

/-! ### ★(iii)(d) の2つの関手

原文 (FrdI p.24):
> co-angular pre-steps. Then the natural functors

原文 (FrdI p.25):
> [obtained by assigning to an arrow φ : A →B the element Div(φ) ∈Φ(A) and to

原文 (FrdI p.25):
> an arrow ψ : B →A the element (ψ∗)−1(Div(ψ)) ∈Φ(A) [since ψ∗: Φ(A) ∼→Φ(B)

原文 (FrdI p.25):
> is a bijection — cf. the fact that ψ is a base-isomorphism!] are equivalences of

★`(ψ*)⁻¹` は `Φ.map (inv (Base ψ))` として書ける——`Φ` は関手なので
`Φ.map (inv h) ∘ Φ.map h = Φ.map (h ≫ inv h) = Φ.map 𝟙 = id`。
**逆写像を選ぶ必要はない。** -/

section
variable [MorphismProperty.IsMultiplicative (coaPreProp P)]

/-- **`_A(𝒞^coa-pre) → Order(Φ(A))`**(★**前置** = コスライス)。 -/
def coaPreUnderFunctor (A : C) :
    Under (⟨A⟩ : WideSubcategory (coaPreProp P)) ⥤
      OrderCat (Φ.val (P.toElem.obj A).base) where
  obj Z := toOrderCat (P.Div Z.hom.hom)
  map {Z W} f := homOfLE (by
    have h : Z.hom.hom ≫ f.right.hom = W.hom.hom :=
      congrArg InducedWideCategory.Hom.hom (Under.w f)
    refine ⟨Φ.map (P.Base Z.hom.hom) (P.Div f.right.hom), ?_⟩
    show P.Div Z.hom.hom + Φ.map (P.Base Z.hom.hom) (P.Div f.right.hom) = P.Div W.hom.hom
    rw [← h, P.Div_comp, show P.degFr f.right.hom = 1 from f.right.property.2.1]
    simp [add_comm])

/-- **`(𝒞^coa-pre)_A → Order(Φ(A))^opp`**(★**後置** = スライス)。 -/
noncomputable def coaPreOverFunctor (A : C) :
    Over (⟨A⟩ : WideSubcategory (coaPreProp P)) ⥤
      (OrderCat (Φ.val (P.toElem.obj A).base))ᵒᵖ where
  obj Z :=
    haveI : IsIso (P.Base Z.hom.hom) := Z.hom.property.2.2
    op (toOrderCat (Φ.map (inv (P.Base Z.hom.hom)) (P.Div Z.hom.hom)))
  map {Z W} f := by
    haveI hZ : IsIso (P.Base Z.hom.hom) := Z.hom.property.2.2
    haveI hW : IsIso (P.Base W.hom.hom) := W.hom.property.2.2
    haveI hf : IsIso (P.Base f.left.hom) := f.left.property.2.2
    refine (homOfLE (x := toOrderCat (Φ.map (inv (P.Base W.hom.hom)) (P.Div W.hom.hom)))
      (y := toOrderCat (Φ.map (inv (P.Base Z.hom.hom)) (P.Div Z.hom.hom))) ?_).op
    have h : f.left.hom ≫ W.hom.hom = Z.hom.hom :=
      congrArg InducedWideCategory.Hom.hom (Over.w f)
    -- ★witness は `(Base Z)⁻¹` を `Div f` に当てたもの。`(Base W)⁻¹` ではない。
    refine ⟨Φ.map (inv (P.Base Z.hom.hom)) (P.Div f.left.hom), ?_⟩
    show Φ.map (inv (P.Base W.hom.hom)) (P.Div W.hom.hom)
        + Φ.map (inv (P.Base Z.hom.hom)) (P.Div f.left.hom)
      = Φ.map (inv (P.Base Z.hom.hom)) (P.Div Z.hom.hom)
    -- `Div Z = Φ(Base f)(Div W) + Div f`(`deg_Fr W = 1`)
    have hdiv : P.Div Z.hom.hom
        = Φ.map (P.Base f.left.hom) (P.Div W.hom.hom) + P.Div f.left.hom := by
      rw [← h, P.Div_comp, show P.degFr W.hom.hom = 1 from W.hom.property.2.1]
      simp
    -- `(Base Z)⁻¹ ≫ … ` の第1項が `(Base W)⁻¹(Div W)` に落ちる
    have hbz : P.Base f.left.hom ≫ P.Base W.hom.hom = P.Base Z.hom.hom := by
      rw [← P.Base_comp, h]
    have hmor : inv (P.Base Z.hom.hom) ≫ P.Base f.left.hom = inv (P.Base W.hom.hom) := by
      rw [← cancel_mono (P.Base W.hom.hom), Category.assoc, hbz, IsIso.inv_hom_id,
        IsIso.inv_hom_id]
    have hcomp : Φ.map (inv (P.Base Z.hom.hom)) (Φ.map (P.Base f.left.hom) (P.Div W.hom.hom))
        = Φ.map (inv (P.Base W.hom.hom)) (P.Div W.hom.hom) := by
      rw [← Φ.map_comp, hmor]
    rw [hdiv, map_add, hcomp]

end

/-! ### ★[FrdI] Definition 1.3 —— Frobenioid

原文 (FrdI p.24):
> if the following conditions are satisfied:
-/

/-- **[FrdI] Definition 1.3** の (iii)(d) 以外の条件。

★(iii)(d) だけを分けるのは、それが (iii)(a)(= `coAngularComp`)を前提しており、
`𝒞^coa-pre` が圏になるのがその後だからである(ファイル冒頭の議論)。 -/
structure FrobenioidCore : Prop where
  /-- **(i)(a)** 底の圏への全射性。

  原文 (FrdI p.24):
  > (i) (Surjectivity to the Base Category via Pull-back Morphisms) (a) Every iso-

  原文 (FrdI p.24):
  > of an isomorphism class of a Frobenius-trivial object of C. (b) If A, B ∈Ob(C),
  -/
  baseSurj : ∀ Y : D, ∃ A : C, IsFrobeniusTrivial P A ∧
    Nonempty ((P.toElem.obj A).base ≅ Y)
  /-- **(i)(b)** 底の同型は pre-step の span で持ち上がる。

  原文 (FrdI p.24):
  > exist pre-steps φ : C →A, ψ : C →B such that α = Base(ψ) ◦Base(φ)−1. (c) For
  -/
  preStepSpan : ∀ (A B : C) (α : (P.toElem.obj A).base ⟶ (P.toElem.obj B).base), IsIso α →
    ∃ (X : C) (φ : X ⟶ A) (ψ : X ⟶ B) (hφ : IsPreStep P φ), IsPreStep P ψ ∧
      α = @inv _ _ _ _ (P.Base φ) hφ.2 ≫ P.Base ψ
  /-- **(i)(c)** `(𝒞^pl-bk)_A → 𝒟_{A_𝒟}` は圏同値。

  原文 (FrdI p.24):
  > equivalence of categories [cf. §0].
  -/
  plBkEquiv : ∀ A : C, (plBkOverFunctor P A).IsEquivalence
  /-- **(ii)** 各次数の Frobenius 型射が存在する。

  原文 (FrdI p.24):
  > Ob(C), n ∈N≥1, there exists a morphism of Frobenius type φ : A →B in C of
  -/
  frobDegSurj : ∀ (A : C) (n : ℕ+), ∃ (B : C) (φ : A ⟶ B),
    IsFrobeniusType P φ ∧ P.degFr φ = n
  /-- **(ii)** その本質的一意性。

  原文 (FrdI p.24):
  > epimorphic] isomorphism β : B ∼→C such that β ◦φ = ψ.
  -/
  frobDegUniq : ∀ (A B E : C) (φ : A ⟶ B) (ψ : A ⟶ E), IsFrobeniusType P φ →
    IsFrobeniusType P ψ → P.degFr φ = P.degFr ψ → ∃ β : B ⟶ E, IsIso β ∧ φ ≫ β = ψ
  /-- **(iii)(a)** co-angular 射は合成で閉じる。

  原文 (FrdI p.24):
  > (iii) (Surjectivity to the Divisor Monoid via Co-angular Morphisms) (a) The
  -/
  coAngularComp : ∀ {A B E : C} (ψ : A ⟶ B) (φ : B ⟶ E),
    IsCoAngular P ψ → IsCoAngular P φ → IsCoAngular P (ψ ≫ φ)
  /-- **(iii)(b)** co-angular pre-step が1本あれば、その始域終域の間の射は全部 co-angular。

  原文 (FrdI p.24):
  > co-angular pre-step φ : A →B, there exists a [uniquely determined] bijection of
  -/
  coAngularOfPreStep : ∀ {A' A : C} (α : A' ⟶ A), IsCoAngular P α → IsPreStep P α →
    ∀ φ : A' ⟶ A, IsCoAngular P φ
  /-- **(iii)(c)** co-angular pre-step は `𝒪^▷` の全単射を誘導する。

  原文 (FrdI p.24):
  > depends only [among the bijections induced by the various co-angular pre-steps
  -/
  otriFwd : ∀ {A B : C} (φ : A ⟶ B), IsCoAngular P φ → IsPreStep P φ →
    ∀ α ∈ OTri P A, ∃! β : End B, β ∈ OTri P B ∧ (φ ≫ β : A ⟶ B) = (α : A ⟶ A) ≫ φ
  /-- **(iii)(c)** その逆向き(全単射であること)。 -/
  otriBwd : ∀ {A B : C} (φ : A ⟶ B), IsCoAngular P φ → IsPreStep P φ →
    ∀ β ∈ OTri P B, ∃! α : End A, α ∈ OTri P A ∧ (φ ≫ β : A ⟶ B) = (α : A ⟶ A) ≫ φ
  /-- **(iii)(c)** 全単射は `Base(φ)` にしか依らない。

  原文 (FrdI p.24):
  > A →B] on Base(φ). (d) Denote by Ccoa-pre ⊆C the subcategory determined by the
  -/
  otriBase : ∀ {A B : C} (φ φ' : A ⟶ B), IsCoAngular P φ → IsPreStep P φ →
    IsCoAngular P φ' → IsPreStep P φ' → P.Base φ = P.Base φ' →
    ∀ (α : End A) (β : End B), (φ ≫ β : A ⟶ B) = (α : A ⟶ A) ≫ φ →
      (φ' ≫ β : A ⟶ B) = (α : A ⟶ A) ≫ φ'
  /-- **(iv)(a)** 任意射の分解 `φ = α ∘ β ∘ γ`(Lean では `φ = γ ≫ β ≫ α`)。

  原文 (FrdI p.25):
  > (iv) (Factorization of Arbitrary Morphisms) Let φ : A →B be a morphism of

  原文 (FrdI p.25):
  > where α is an pull-back morphism, β is a pre-step, and γ is a morphism of Frobenius
  -/
  arbFactor : ∀ {A B : C} (φ : A ⟶ B), ∃ (X Y : C) (γ : A ⟶ X) (β : X ⟶ Y) (α : Y ⟶ B),
    φ = γ ≫ β ≫ α ∧ IsFrobeniusType P γ ∧ IsPreStep P β ∧ IsPullBack P α
  /-- **(iv)(a)** 分解の一意性(同型 `δ`、`ε` を除く)。

  原文 (FrdI p.25):
  > type; this factorization is unique, up to replacing the triple (α, β, γ) by a triple of
  -/
  arbFactorUniq : ∀ {A B : C} (X Y X' Y' : C)
      (γ : A ⟶ X) (β : X ⟶ Y) (α : Y ⟶ B) (γ' : A ⟶ X') (β' : X' ⟶ Y') (α' : Y' ⟶ B),
    γ ≫ β ≫ α = γ' ≫ β' ≫ α' →
    IsFrobeniusType P γ → IsPreStep P β → IsPullBack P α →
    IsFrobeniusType P γ' → IsPreStep P β' → IsPullBack P α' →
    ∃ (δ : Y ≅ Y') (ε : X ≅ X'),
      α' = δ.inv ≫ α ∧ β' = ε.inv ≫ β ≫ δ.hom ∧ γ' = γ ≫ ε.hom
  /-- **(iv)(b)** pull-back morphism は LB-invertible かつ linear。

  原文 (FrdI p.25):
  > pull-back morphism of C is LB-invertible and linear.
  -/
  pullBackLB : ∀ {A B : C} (α : A ⟶ B), IsPullBack P α → IsLBInvertible P α ∧ IsLinear P α
  /-- **(v)(a)** pre-step は monomorphism。

  原文 (FrdI p.25):
  > φ is a monomorphism. (b) φ admits a factorization
  -/
  preStepMono : ∀ {A B : C} (φ : A ⟶ B), IsPreStep P φ → Mono φ
  /-- **(v)(b)** pre-step は「isometric ∘ co-angular」に分解する。

  原文 (FrdI p.25):
  > where α is an isometric pre-step, and β is a co-angular pre-step; this factorization
  -/
  preStepFactor : ∀ {A B : C} (φ : A ⟶ B), IsPreStep P φ →
    ∃ (X : C) (β : A ⟶ X) (α : X ⟶ B),
      φ = β ≫ α ∧ IsCoAngular P β ∧ IsPreStep P β ∧ IsIsometric P α ∧ IsPreStep P α
  /-- **(v)(c)** pre-step は「co-angular ∘ isometric」にも分解する。

  原文 (FrdI p.25):
  > (v) (Factorization of Pre-steps) Let φ : A →B be a pre-step of C. Then: (a)
  -/
  preStepFactor' : ∀ {A B : C} (φ : A ⟶ B), IsPreStep P φ →
    ∃ (X : C) (β : A ⟶ X) (α : X ⟶ B),
      φ = β ≫ α ∧ IsIsometric P β ∧ IsPreStep P β ∧ IsCoAngular P α ∧ IsPreStep P α
  /-- **(vi)** 単元を除く忠実性。

  原文 (FrdI p.25):
  > equivalent co-angular pre-steps of C. Then there exists a [necessarily unique] α ∈
  -/
  faithfulUpToUnits : ∀ {A B : C} (φ ψ : A ⟶ B), BaseEquivalent P φ ψ →
    MetricallyEquivalent P φ ψ → IsCoAngular P φ → IsPreStep P φ →
    IsCoAngular P ψ → IsPreStep P ψ →
    ∃! α : End B, α ∈ OTimes P B ∧ φ = ψ ≫ (α : B ⟶ B)
  /-- **(vii)(a)** isotropic hull が存在する。

  原文 (FrdI p.25):
  > (vii) (Isotropic Objects) (a) For every A ∈Ob(C), there exists a [necessarily
  -/
  isotropicHullExists : ∀ A : C, ∃ (B : C) (φ : A ⟶ B), IsIsotropicHull P φ
  /-- **(vii)(b)** isotropic な対象からの射の終域も isotropic。

  原文 (FrdI p.25):
  > isotropic, and A →C is a morphism of C, then C is also isotropic.
  -/
  isotropicClosed : ∀ {A B : C} (φ : A ⟶ B), IsIsotropic P A → IsIsotropic P B

/-- ★**[FrdI] Definition 1.3 —— Frobenioid**。

原文 (FrdI p.24):
> that the pre-Frobenioid C →FΦ [i.e., C equipped with this functor] is a Frobenioid

(iii)(d) は `core` の (iii)(a) から `𝒞^coa-pre` を圏にしてから述べる。 -/
structure Frobenioid : Prop where
  /-- (iii)(d) 以外の全条件 -/
  core : FrobenioidCore P
  /-- **(iii)(d)** `_A(𝒞^coa-pre) → Order(Φ(A))` は圏同値(★前置 = コスライス)。

  原文 (FrdI p.25):
  > categories.
  -/
  coaPreUnderEquiv :
    letI := coaPreProp_isMultiplicative P core.coAngularComp
    ∀ A : C, (coaPreUnderFunctor P A).IsEquivalence
  /-- **(iii)(d)** `(𝒞^coa-pre)_A → Order(Φ(A))^opp` は圏同値(★後置 = スライス)。 -/
  coaPreOverEquiv :
    letI := coaPreProp_isMultiplicative P core.coAngularComp
    ∀ A : C, (coaPreOverFunctor P A).IsEquivalence

end ABC3.Found.FrdI

/-! ### ★出典の紐付け(`.src`) -/

namespace ABC3.Found.FrdI

def OrderCat.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 12, item := "§0 Monoids — Order(M)",
    sectionId := "frdi-s0-order" }

def PlBk.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 24, item := "Definition 1.3, (i), (c) — 𝒞^pl-bk",
    sectionId := "frdi-def-1-3-ic" }

def coaPreProp.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 24, item := "Definition 1.3, (iii), (d) — 𝒞^coa-pre",
    sectionId := "frdi-def-1-3-iiid" }

def plBkOverFunctor.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 24, item := "Definition 1.3, (i), (c) — (𝒞^pl-bk)_A → 𝒟_{A_𝒟}",
    sectionId := "frdi-def-1-3-ic" }

def coaPreUnderFunctor.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 24, item := "Definition 1.3, (iii), (d) — _A(𝒞^coa-pre) → Order(Φ(A))",
    sectionId := "frdi-def-1-3-iiid" }

def coaPreOverFunctor.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 24, item := "Definition 1.3, (iii), (d) — (𝒞^coa-pre)_A → Order(Φ(A))^opp",
    sectionId := "frdi-def-1-3-iiid" }

def FrobenioidCore.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 24, item := "Definition 1.3, (i)(ii)(iii)(a)(b)(c)(iv)(v)(vi)(vii)",
    sectionId := "frdi-def-1-3-i" }

def Frobenioid.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 24, item := "Definition 1.3 — Frobenioid",
    sectionId := "frdi-def-1-3-i" }

end ABC3.Found.FrdI
