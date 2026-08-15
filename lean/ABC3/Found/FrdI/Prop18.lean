import ABC3.Found.FrdI.Prop15
import ABC3.Found.FrdI.Prop17

/-!
# [FrdI] Proposition 1.8 —— Pre-steps

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、
物理 p.30(**400 dpi 目視確認 2026-08-15**)。

原文 (FrdI p.30):
> Proposition 1.8. (Pre-steps) Let Φ be a divisorial monoid on a connected,

原文 (FrdI p.30):
> totally epimorphic category D; C →FΦ a Frobenioid. Then:

## ★この命題の規模(測定)

**3条 (i)–(iii)、主張は 6**:

| 条 | 主張 | 追加の仮定 |
|---|---|---|
| (i) | 2 | `𝒞 → 𝒟` が **full** / `𝒟` が **of Aut-type** |
| (ii) | 1 | `𝒞` が **metrically trivial** かつ **Aut-ample** type |
| (iii) | 3 | 追加仮定なし |

★**原文の証明は 2 段落**(私は着手前の測定で「3段落」と言ったが、
目視すると (iii) は第2段落の最後の1文である)。

## ★★`Definition 1.3, (iii), (d)` の**初使用**

原文 (FrdI p.30):
> we observe that the various equivalences of assertion (iii) follow formally from the

★我々は `Definition 1.3, (iii), (d)` の 2 本の圏同値を `Frobenioid` の
フィールドとして持っていたが、**下流で一度も使っていなかった**。
(iii) がその**初使用**である。どこでどう使ったかは各定理の docstring に書く。
-/

namespace ABC3.Found.FrdI

open CategoryTheory Opposite

universe v u w u2 v2

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ)

/-! ## ★(i) —— pre-step と linear End-equivalence

原文 (FrdI p.30):
> (i) If the natural projection functor C → D is full, then every pre-step of

原文 (FrdI p.30):
> C is a linear End-equivalence. If D is of Aut-type [cf. §0], then every linear
-/

/-- **(i) 前半** —— `𝒞 → 𝒟` が full なら、pre-step は linear End-equivalence。

★`linear` は pre-step の定義そのもの。`End-equivalence` は
「`Base φ` が同型なので逆射 `inv (Base φ)` が `𝒟` にあり、
**fullness でそれを `𝒞` へ持ち上げる**」で出る。

★**「follows formally」は本当に formally だった** —— 2 行。 -/
theorem prop_1_8_i_mp (hfull : (P.proj).Full) {A B : C} (φ : A ⟶ B)
    (hφ : IsPreStep P φ) : IsLinear P φ ∧ IsEndEquivalence φ := by
  refine ⟨hφ.1, ?_⟩
  haveI := hfull
  have hbi : IsIso (P.Base φ) := hφ.2
  obtain ⟨g, -, -⟩ := hbi.out
  obtain ⟨ψ, -⟩ := (P.proj).map_surjective (X := B) (Y := A) g
  exact ⟨ψ⟩

/-- **(i) 後半** —— `𝒟` が of Aut-type なら、linear End-equivalence は pre-step。

★`End-equivalence` から `ψ : B ⟶ A` を取ると、`Base φ ≫ Base ψ` も
`Base ψ ≫ Base φ` も `𝒟` の**自己射**なので、`of Aut-type` から同型。
両側から同型なら `Base φ` 自身が同型(`isIso_of_comp_isIso_both`)。

★**ここも formally** —— ただし「両側から同型なら同型」は
mathlib の `isIso_of_mono_of_isSplitEpi` を経由する**1つの補題**を要した。 -/
theorem prop_1_8_i_mpr (hAut : IsOfAutType D) {A B : C} (φ : A ⟶ B)
    (hlin : IsLinear P φ) (hee : IsEndEquivalence φ) : IsPreStep P φ := by
  refine ⟨hlin, ?_⟩
  obtain ⟨ψ⟩ := hee
  haveI : IsIso (P.Base φ ≫ P.Base ψ) := by
    rw [← P.Base_comp]; exact hAut _ (P.Base (φ ≫ ψ))
  haveI : IsIso (P.Base ψ ≫ P.Base φ) := by
    rw [← P.Base_comp]; exact hAut _ (P.Base (ψ ≫ φ))
  exact isIso_of_comp_isIso_both (P.Base φ) (P.Base ψ)

/-! ## ★(ii) —— co-angular pre-step の abstract equivalence による特徴づけ

原文 (FrdI p.30):
> (ii) Suppose further that C is of metrically trivial and Aut-ample type.

原文 (FrdI p.30):
> Then a morphism of C is a co-angular pre-step if and only if it is abstractly
-/

/-- ★**`𝒞` のすべての自己射は co-angular**。

`Definition 1.3, (iii), (b)`(`coAngularOfPreStep`)を `𝟙_A` に当てるだけ ——
`𝟙_A` は co-angular pre-step なので、`A ⟶ A` の射は**すべて** co-angular になる。

★原文が (ii) の括弧で「hence co-angular, by Definition 1.3, (iii), (b)」と
書いているのはこれである。 -/
theorem isCoAngular_of_endo (F : FrobenioidCore P) {A : C} (e : A ⟶ A) :
    IsCoAngular P e :=
  F.coAngularOfPreStep (𝟙 A) (isCoAngular_id P A) (isPreStep_id P A) e

/-- **(ii) ⇒** —— co-angular pre-step は base-identity pre-step 自己射に
abstractly equivalent。

★**2つの仮定がそれぞれ別の役割を持つ**:
* **metrically trivial** で `B ≅ A` を得て、`φ` を**自己射に変える**
* **Aut-ample** で `Base` を打ち消す自己同型 `u` を得て、**base-identity にする**

★どちらか一方では届かない(下の負の対照を見よ)。 -/
theorem prop_1_8_ii_mp (hmt : ∀ X : C, IsMetricallyTrivial P X)
    (hAA : ∀ X : C, IsAutAmple P X) {A B : C} (φ : A ⟶ B)
    (hco : IsCoAngular P φ) (hst : IsPreStep P φ) :
    ∃ (X : C) (e : X ⟶ X), IsBaseIdentity P e ∧ IsPreStep P e ∧
      AbstractlyEquivalent φ e := by
  obtain ⟨b⟩ := hmt A B φ hco hst
  haveI hbi : IsIso (P.Base (φ ≫ b.hom)) :=
    isBaseIsomorphism_comp P hst.2 (isBaseIsomorphism_of_isIso P b.hom)
  obtain ⟨u, huiso, hu⟩ := hAA A (inv (P.Base (φ ≫ b.hom))) inferInstance
  haveI := huiso
  refine ⟨A, (φ ≫ b.hom) ≫ (u : A ⟶ A), ?_, ?_, Iso.refl A, b ≪≫ asIso (u : A ⟶ A), ?_⟩
  · show P.Base ((φ ≫ b.hom) ≫ (u : A ⟶ A)) = P.Base (𝟙 A)
    rw [P.Base_comp, hu, IsIso.hom_inv_id, P.Base_id]
  · exact IsPreStep.comp P (IsPreStep.comp P hst (isPreStep_of_isIso P b.hom))
      (isPreStep_of_isIso P (u : A ⟶ A))
  · show φ ≫ (b.hom ≫ (u : A ⟶ A)) = 𝟙 A ≫ ((φ ≫ b.hom) ≫ (u : A ⟶ A))
    rw [Category.id_comp, Category.assoc]

/-- **(ii) ⇐** —— base-identity pre-step 自己射に abstractly equivalent なら
co-angular pre-step。

★**こちらは追加仮定を使わない** —— metrically trivial も Aut-ample も要らない。
`isCoAngular_of_endo` と「同型は co-angular pre-step」だけで出る。
★原文も「follows formally」としか書かないが、**仮定の非対称性は書いていない**。 -/
theorem prop_1_8_ii_mpr (F : FrobenioidCore P) {A B X : C} (φ : A ⟶ B) (e : X ⟶ X)
    (_hbe : IsBaseIdentity P e) (hes : IsPreStep P e)
    (hae : AbstractlyEquivalent φ e) : IsCoAngular P φ ∧ IsPreStep P φ := by
  obtain ⟨a, b, hab⟩ := hae
  have hφ : φ = a.hom ≫ e ≫ b.inv := by
    rw [← Category.assoc, ← hab, Category.assoc, b.hom_inv_id, Category.comp_id]
  rw [hφ]
  refine ⟨F.coAngularComp _ _ (isCoAngular_of_isIso P a.hom)
      (F.coAngularComp _ _ (isCoAngular_of_endo P F e) (isCoAngular_of_isIso P b.inv)), ?_⟩
  exact IsPreStep.comp P (isPreStep_of_isIso P a.hom)
    (IsPreStep.comp P hes (isPreStep_of_isIso P b.inv))

/-! ## ★★(iii) —— `Definition 1.3, (iii), (d)` の初使用

原文 (FrdI p.30):
> (iii) An object A ∈Ob(C) is non-group-like if and only if there exists a

原文 (FrdI p.30):
> co-angular step A →B; alternatively, an object A ∈Ob(C) is non-group-like

原文 (FrdI p.30):
> if and only if there exists a co-angular step B →A. Also, if A, B ∈Ob(C) are

★**圏同値のどの部分を使うか**:
* `⟸`(step から non-group-like)は **充満性**
* `⟹`(non-group-like から step)は **本質的全射性**

★**忠実性は使わない。** 原文は「the equivalences of categories」と一括で言うが、
実際に効くのは 2 つの部分だけである。 -/

section CoaPreUse

variable [MorphismProperty.IsMultiplicative (coaPreProp P)]

/-- ★`Φ(A)` の元が可逆であることと、対応する co-angular pre-step が同型であることは
`MLe` を通じて対応する —— 可逆なら `MLe _ 0` が成り立つ。 -/
private theorem mle_zero_of_isAddUnit {M : Type w} [AddCommMonoid M] {a : M}
    (h : IsAddUnit a) : MLe a (0 : M) := by
  obtain ⟨U, rfl⟩ := h
  exact ⟨U.neg, U.val_neg⟩

/-- ★`MLe` が両向きに成り立つとき、整域(= 簡約的)なら差は可逆元。 -/
private theorem isAddUnit_iff_of_mle_both {M : Type w} [AddCommMonoid M] [IsCancelAdd M]
    {a d : M} (h1 : MLe a d) (h2 : MLe d a) : IsAddUnit a ↔ IsAddUnit d := by
  obtain ⟨c, hc⟩ := h1
  obtain ⟨c', hc'⟩ := h2
  have hcc : c + c' = 0 := by
    refine add_left_cancel (a := a) ?_
    rw [← add_assoc, hc, hc', add_zero]
  have hcu : IsAddUnit c := ⟨⟨c, c', hcc, by rw [add_comm]; exact hcc⟩, rfl⟩
  constructor
  · intro ha; rw [← hc]; exact ha.add hcu
  · intro hd; rw [← hc']; exact hd.add ⟨⟨c', c, by rw [add_comm]; exact hcc, hcc⟩, rfl⟩

/-- **(iii) 第1主張** —— `A` が non-group-like ⟺ `A` から出る co-angular step がある。

★★**`Definition 1.3, (iii), (d)` のコスライス側の圏同値を実際に使う**:
* `⟸`: `Div φ` が可逆だと仮定すると `MLe (Div φ) 0` なので、
  **充満性**で `Z ⟶ Z₀`(`Z₀` は `𝟙_A`)が持ち上がり、`φ ≫ (その射) = 𝟙_A`。
  `𝒞` が totally epimorphic なので `φ` が同型になり、step であることに反する。
* `⟹`: 可逆でない `d ∈ Φ(A)` を取り、**本質的全射性**で `F.obj Z ≅ toOrderCat d`
  なる `Z` を得る。`Div Z.hom.hom` も可逆でないので、その co-angular pre-step は
  同型でない、すなわち step。 -/
theorem prop_1_8_iii_out (hequiv : ∀ X : C, (coaPreUnderFunctor P X).IsEquivalence) (A : C) :
    (¬ IsGroupLikeObj P A) ↔
      ∃ (B : C) (φ : A ⟶ B), IsCoAngular P φ ∧ IsStep P φ := by
  haveI hcancel := isCancelAdd_of_isIntegralMonoid _ (P.divisorial (P.toElem.obj A).base).1.1
  haveI := (hequiv A).full
  haveI := (hequiv A).essSurj
  constructor
  · intro hng
    obtain ⟨d, hd⟩ : ∃ d : Φ.val (P.toElem.obj A).base, ¬ IsAddUnit d := by
      by_contra hc
      exact hng ((isGroupLike_iff _).mpr (not_exists_not.mp hc))
    set Z := (coaPreUnderFunctor P A).objPreimage (toOrderCat d) with hZ
    have hiso := (coaPreUnderFunctor P A).objObjPreimageIso (toOrderCat d)
    have h1 : MLe (P.Div Z.hom.hom) d := leOfHom hiso.hom
    have h2 : MLe d (P.Div Z.hom.hom) := leOfHom hiso.inv
    have hnu : ¬ IsAddUnit (P.Div Z.hom.hom) := fun hu =>
      hd ((isAddUnit_iff_of_mle_both h1 h2).mp hu)
    refine ⟨Z.right.obj, Z.hom.hom, Z.hom.property.1, Z.hom.property.2, ?_⟩
    intro hi
    haveI := hi
    have hz : P.Div Z.hom.hom = 0 := isIsometric_of_isIso P _
    refine hnu ?_
    rw [hz]
    exact isAddUnit_zero
  · rintro ⟨B, φ, hco, hst, hnotiso⟩ hgl
    have hu : IsAddUnit (P.Div φ) := (isGroupLike_iff _).mp hgl _
    set Z : Under (⟨A⟩ : WideSubcategory (coaPreProp P)) :=
      Under.mk (show (⟨A⟩ : WideSubcategory (coaPreProp P)) ⟶ ⟨B⟩ from ⟨φ, hco, hst⟩) with hZ
    set Z₀ : Under (⟨A⟩ : WideSubcategory (coaPreProp P)) :=
      Under.mk (𝟙 (⟨A⟩ : WideSubcategory (coaPreProp P))) with hZ₀
    have hle : (coaPreUnderFunctor P A).obj Z ⟶ (coaPreUnderFunctor P A).obj Z₀ := by
      refine homOfLE ?_
      show MLe (P.Div φ) (P.Div (𝟙 A))
      rw [P.Div_id A]
      exact mle_zero_of_isAddUnit hu
    obtain ⟨f, -⟩ := (coaPreUnderFunctor P A).map_surjective hle
    have hw : φ ≫ f.right.hom = 𝟙 A :=
      congrArg InducedWideCategory.Hom.hom (Under.w f)
    haveI : IsIso (φ ≫ f.right.hom) := by rw [hw]; exact IsIso.id _
    exact hnotiso (isIso_of_isIso_comp P.totEpiC φ f.right.hom).1

/-- **(iii) 第2主張** —— `A` が non-group-like ⟺ `A` へ入る co-angular step がある。

★★**今度はスライス側の圏同値**(`coaPreOverFunctor`、`Order(Φ(A))^opp` へ行くほう)。
使う部分は第1主張と同じく **充満性**と**本質的全射性**だが、
★**射の向きが逆になるので、`𝟙_A` から `Z` への射を持ち上げる**ことになる。
そこから `g ≫ φ = 𝟙_A` が出て、`𝒞` の totally epimorphic 性で `φ` が同型。 -/
theorem prop_1_8_iii_in (hequiv : ∀ X : C, (coaPreOverFunctor P X).IsEquivalence) (A : C) :
    (¬ IsGroupLikeObj P A) ↔
      ∃ (B : C) (φ : B ⟶ A), IsCoAngular P φ ∧ IsStep P φ := by
  haveI hcancel := isCancelAdd_of_isIntegralMonoid _ (P.divisorial (P.toElem.obj A).base).1.1
  haveI := (hequiv A).full
  haveI := (hequiv A).essSurj
  constructor
  · intro hng
    obtain ⟨d, hd⟩ : ∃ d : Φ.val (P.toElem.obj A).base, ¬ IsAddUnit d := by
      by_contra hc
      exact hng ((isGroupLike_iff _).mpr (not_exists_not.mp hc))
    set Z := (coaPreOverFunctor P A).objPreimage (op (toOrderCat d)) with hZ
    have hiso := (coaPreOverFunctor P A).objObjPreimageIso (op (toOrderCat d))
    haveI hbz : IsIso (P.Base Z.hom.hom) := Z.hom.property.2.2
    have h1 : MLe d (Φ.map (inv (P.Base Z.hom.hom)) (P.Div Z.hom.hom)) := leOfHom hiso.hom.unop
    have h2 : MLe (Φ.map (inv (P.Base Z.hom.hom)) (P.Div Z.hom.hom)) d := leOfHom hiso.inv.unop
    have hnu : ¬ IsAddUnit (Φ.map (inv (P.Base Z.hom.hom)) (P.Div Z.hom.hom)) := fun hu =>
      hd ((isAddUnit_iff_of_mle_both h2 h1).mp hu)
    refine ⟨Z.left.obj, Z.hom.hom, Z.hom.property.1, Z.hom.property.2, ?_⟩
    intro hi
    haveI := hi
    have hz : P.Div Z.hom.hom = 0 := isIsometric_of_isIso P _
    refine hnu ?_
    rw [hz, map_zero]
    exact isAddUnit_zero
  · rintro ⟨B, φ, hco, hst, hnotiso⟩ hgl
    haveI hbz : IsIso (P.Base φ) := hst.2
    have hu : IsAddUnit (Φ.map (inv (P.Base φ)) (P.Div φ)) := (isGroupLike_iff _).mp hgl _
    set Z : Over (⟨A⟩ : WideSubcategory (coaPreProp P)) :=
      Over.mk (show (⟨B⟩ : WideSubcategory (coaPreProp P)) ⟶ ⟨A⟩ from ⟨φ, hco, hst⟩) with hZ
    set Z₀ : Over (⟨A⟩ : WideSubcategory (coaPreProp P)) :=
      Over.mk (𝟙 (⟨A⟩ : WideSubcategory (coaPreProp P))) with hZ₀
    have h0 : ((coaPreOverFunctor P A).obj Z₀).unop
        = toOrderCat (0 : Φ.val (P.toElem.obj A).base) := by
      show toOrderCat (Φ.map _ (P.Div (𝟙 A))) = toOrderCat 0
      rw [P.Div_id, map_zero]
    have hle : (coaPreOverFunctor P A).obj Z₀ ⟶ (coaPreOverFunctor P A).obj Z := by
      refine (homOfLE (x := ((coaPreOverFunctor P A).obj Z).unop)
        (y := ((coaPreOverFunctor P A).obj Z₀).unop) ?_).op
      rw [h0]
      exact mle_zero_of_isAddUnit hu
    obtain ⟨g, -⟩ := (coaPreOverFunctor P A).map_surjective hle
    have hw : g.left.hom ≫ φ = 𝟙 A :=
      congrArg InducedWideCategory.Hom.hom (Over.w g)
    haveI : IsIso (g.left.hom ≫ φ) := by rw [hw]; exact IsIso.id _
    exact hnotiso (isIso_of_isIso_comp P.totEpiC g.left.hom φ).2

end CoaPreUse

/-- **(iii) 第3主張** —— base-isomorphic な対象は group-like 性が一致する。

★**ここは `Definition 1.3, (iii), (d)` を使わない** —— `Φ` の関手性だけで出る。
原文は3主張を「follow formally from the definitions and the equivalences of
categories of Definition 1.3, (iii), (d)」と一括で言うが、
**3番目は圏同値に依らない**。 -/
theorem prop_1_8_iii_baseIso {A B : C} (h : BaseIsomorphic P A B) :
    IsGroupLikeObj P A ↔ IsGroupLikeObj P B := by
  obtain ⟨e⟩ := h
  have key : ∀ (X Y : C), ((P.toElem.obj X).base ≅ (P.toElem.obj Y).base) →
      IsGroupLikeObj P X → IsGroupLikeObj P Y := by
    intro X Y f hX
    refine (isGroupLike_iff _).mpr fun a => ?_
    have h1 : IsAddUnit (Φ.map f.hom a) := (isGroupLike_iff _).mp hX _
    have h2 : IsAddUnit (Φ.map f.inv (Φ.map f.hom a)) := h1.map _
    have h3 : Φ.map f.inv (Φ.map f.hom a) = a := by
      rw [← MonoidOn.map_comp, f.inv_hom_id, MonoidOn.map_id]
    rwa [h3] at h2
  exact ⟨key A B e, key B A e.symm⟩

/-! ## ★負の対照 —— (iii) の同値が**両側とも非自明**であること

★(iii) は「non-group-like ⟺ co-angular step が存在」という同値である。
これが空でないことを言うには、**両側が成り立つ模型**と
**両側が成り立たない模型**の 2 つが要る。

| 模型 | `Φ` | group-like か | co-angular step | 判定 |
|---|---|---|---|---|
| `wP` | `ℕ` | ★**でない** | ★**ある**(`wDiv1`) | 両側成立 |
| `zP` | `ℤ` | ★**である** | ★**ない**(pre-step はすべて同型) | 両側不成立 |
-/

/-- ★模型 `wP`(`Φ = ℕ`)の `wTop` は **group-like でない**。 -/
theorem wNot_groupLikeObj : ¬ IsGroupLikeObj wP wTop := not_isGroupLike_nat

/-- ★その `wTop` からは **co-angular step が出る**(`wDiv1`)。
★(iii) の同値の**両側が成り立つ**実例。 -/
theorem wCoAngularStep : IsCoAngular wP wDiv1 ∧ IsStep wP wDiv1 :=
  ⟨wIsCoAngular _, wStep_div1⟩

/-- ★`Φ = ℤ`(群)の初等 Frobenioid。`ℤ` は pre-divisorial なので
`Proposition 1.5` により本当に Frobenioid である。 -/
abbrev zΦ : MonoidOn.{0, 0, 0} Vee := MonoidOn.const Vee ℤ

/-- `𝔽_ℤ` の pre-Frobenioid 構造。 -/
abbrev zP : PreFrobenioid (ElemFrobCat zΦ) zΦ.charOn :=
  elemPreFrobenioid zΦ isTotallyEpimorphic_vee (fun _ => isPreDivisorial_int)

/-- ★`𝔽_ℤ` は Frobenioid(`Proposition 1.5` の帰結)。 -/
theorem zIsFrobenioid : Frobenioid zP :=
  elemFrob_isFrobenioid zΦ isTotallyEpimorphic_vee (fun _ => isPreDivisorial_int)

/-- ★`𝔽_ℤ` のすべての対象は **group-like**。 -/
theorem zGroupLikeObj (A : ElemFrobCat zΦ) : IsGroupLikeObj zP A :=
  elemFrob_groupLike zΦ isTotallyEpimorphic_vee (fun _ => isPreDivisorial_int)
    (fun _ => isGroupLike_int) A

/-- ★★`𝔽_ℤ` には **step が1つも無い** —— pre-step はすべて同型。

`ℤ` では零因子がつねに可逆なので、`ElemFrobCat.isIso_iff` の3条件が
pre-step の2条件から自動で揃う。
★(iii) の同値の**両側が成り立たない**実例。 -/
theorem zNo_step {A B : ElemFrobCat zΦ} (φ : A ⟶ B) (h : IsPreStep zP φ) : IsIso φ :=
  (ElemFrobCat.isIso_iff φ).mpr ⟨h.2, (isGroupLike_iff ℤ).mp isGroupLike_int _, h.1⟩

/-- ★したがって `𝔽_ℤ` には co-angular step が無い。 -/
theorem zNo_coAngularStep (A : ElemFrobCat zΦ) :
    ¬ ∃ (B : ElemFrobCat zΦ) (φ : A ⟶ B), IsCoAngular zP φ ∧ IsStep zP φ := by
  rintro ⟨B, φ, -, hst, hnot⟩
  exact hnot (zNo_step φ hst)

/-! ### ★出典の紐付け(`.src`) -/

def prop_1_8_i_mp.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 30, item := "Proposition 1.8, (i)",
    sectionId := "frdi-prop-1-8-i" }

def prop_1_8_ii_mp.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 30, item := "Proposition 1.8, (ii)",
    sectionId := "frdi-prop-1-8-ii" }

def prop_1_8_iii_out.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 30, item := "Proposition 1.8, (iii)",
    sectionId := "frdi-prop-1-8-iii" }

/-! ## ★★命題全体の `.src`(2026-08-16 に新設した段)

★条つきは **locator の記録**、条なしは ★★**命題全体を完全に実装したという主張**。

★**`Proposition 1.8` は 6 主張すべてが揃っている**:

| 条 | 主張 | 宣言 |
|---|---|---|
| (i) | `𝒞 → 𝒟` が充満なら pre-step は linear End-equivalence | `prop_1_8_i_mp` |
| (i) | ★**`𝒟` が Aut-type なら**逆 | `prop_1_8_i_mpr` |
| (ii) | ⟹ | `prop_1_8_ii_mp` |
| (ii) | ⟸ | `prop_1_8_ii_mpr` |
| (iii) | non-group-like ⟺ co-angular step `A → B` の存在 | `prop_1_8_iii_out` |
| (iii) | ★**alternatively**: `B → A` の側 | `prop_1_8_iii_in` |
| (iii) | ★**Also**: base-isomorphic なら group-like 性が一致 | `prop_1_8_iii_baseIso` |

★§0 の語彙の写しも照合済み ——
`IsEndEquivalence φ := Nonempty (B ⟶ A)` は原文 p.14 の
「an arrow A →B of C is an End-equivalence **if there exists an arrow B →A in C**」
そのものであり、★**射自身が条件に現れないのは原文どおり**である。
-/

def prop_1_8.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 30, item := "Proposition 1.8",
    sectionId := "frdi-prop-1-8-i" }

end ABC3.Found.FrdI
