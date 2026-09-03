/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop18
import ABC3.Found.FrdI.Prop19.Istr

/-!
# Prop19 —— (ii)(iii) の梱包

☆もとの 1 枚を**入れ子の切れ目**で割ったものである(第 1457)。
-/

namespace ABC3.Found.FrdI

open CategoryTheory Opposite

universe v u w u2 v2

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ)


section ImtrPreFunctor

variable (F : FrobenioidCore P)

include F in
/-- ★★**`𝒞^imtr-pre_A` は thin**(hom 集合が高々1元)。

★`Definition 1.3, (v), (a)`(pre-step は mono)だけから出る。
これにより `φ_*` / `φ^*` の**忠実性と関手則が無料**になる。

★**手順3**(移送・補助補題を書くときは最初に `include F in`)—— ここでも同じ。 -/
theorem imtrPreOver_hom_subsingleton {A : C} {Z W : Over (⟨A⟩ : ImtrPre P)} (f g : Z ⟶ W) :
    f = g :=
  Over.OverMorphism.ext (InducedWideCategory.Hom.ext
    (imtrPre_hom_uniq P F Z.hom.hom W.hom.hom W.hom.property.2 _ _
      (congrArg InducedWideCategory.Hom.hom (Over.w f))
      (congrArg InducedWideCategory.Hom.hom (Over.w g))))

/-! ### `φ_*`(base-isomorphism から) -/

/-- ★分解の左因子(co-angular base-isomorphism)。 -/
noncomputable def pushMid {A B : C} (φ : A ⟶ B) (hφ : IsBaseIsomorphism P φ)
    {Cc : C} (ε : Cc ⟶ A) (hεi : IsIsometric P ε) (hεs : IsPreStep P ε) :
    Cc ⟶ pushObj P F φ hφ ε hεi hεs :=
  (prop_1_9_ii_obj P F φ hφ ε hεi hεs).choose_spec.choose

theorem pushFac {A B : C} (φ : A ⟶ B) (hφ : IsBaseIsomorphism P φ)
    {Cc : C} (ε : Cc ⟶ A) (hεi : IsIsometric P ε) (hεs : IsPreStep P ε) :
    ε ≫ φ = pushMid P F φ hφ ε hεi hεs ≫ pushHom P F φ hφ ε hεi hεs :=
  (prop_1_9_ii_obj P F φ hφ ε hεi hεs).choose_spec.choose_spec.choose_spec.1

theorem pushMid_spec {A B : C} (φ : A ⟶ B) (hφ : IsBaseIsomorphism P φ)
    {Cc : C} (ε : Cc ⟶ A) (hεi : IsIsometric P ε) (hεs : IsPreStep P ε) :
    IsCoAngular P (pushMid P F φ hφ ε hεi hεs) ∧
      IsBaseIsomorphism P (pushMid P F φ hφ ε hεi hεs) :=
  (prop_1_9_ii_obj P F φ hφ ε hεi hεs).choose_spec.choose_spec.choose_spec.2.1

/-- ★`φ_*` の射への割り当て。★**選択は要らない**が、`∃` から取るので
`.choose` を使う —— 一意性は `imtrPre_hom_uniq` にある。 -/
noncomputable def pushMap {A B : C} (φ : A ⟶ B) (hφ : IsBaseIsomorphism P φ)
    {Z W : Over (⟨A⟩ : ImtrPre P)} (f : Z ⟶ W) :
    pushObj P F φ hφ Z.hom.hom Z.hom.property.1 Z.hom.property.2 ⟶
      pushObj P F φ hφ W.hom.hom W.hom.property.1 W.hom.property.2 :=
  (prop_1_9_ii_hom P F φ
    Z.hom.hom (pushMid P F φ hφ Z.hom.hom Z.hom.property.1 Z.hom.property.2)
      (pushHom P F φ hφ Z.hom.hom Z.hom.property.1 Z.hom.property.2)
    W.hom.hom (pushMid P F φ hφ W.hom.hom W.hom.property.1 W.hom.property.2)
      (pushHom P F φ hφ W.hom.hom W.hom.property.1 W.hom.property.2)
    (pushFac P F φ hφ Z.hom.hom Z.hom.property.1 Z.hom.property.2)
    (pushMid_spec P F φ hφ Z.hom.hom Z.hom.property.1 Z.hom.property.2).1
    (pushMid_spec P F φ hφ Z.hom.hom Z.hom.property.1 Z.hom.property.2).2
    (pushHom_spec P F φ hφ Z.hom.hom Z.hom.property.1 Z.hom.property.2).1
    (pushHom_spec P F φ hφ Z.hom.hom Z.hom.property.1 Z.hom.property.2).2
    (pushFac P F φ hφ W.hom.hom W.hom.property.1 W.hom.property.2)
    (pushMid_spec P F φ hφ W.hom.hom W.hom.property.1 W.hom.property.2).2
    (pushHom_spec P F φ hφ W.hom.hom W.hom.property.1 W.hom.property.2).1
    (pushHom_spec P F φ hφ W.hom.hom W.hom.property.1 W.hom.property.2).2
    f.left.hom f.left.property.2
    (congrArg InducedWideCategory.Hom.hom (Over.w f))).choose

theorem pushMap_spec {A B : C} (φ : A ⟶ B) (hφ : IsBaseIsomorphism P φ)
    {Z W : Over (⟨A⟩ : ImtrPre P)} (f : Z ⟶ W) :
    (IsIsometric P (pushMap P F φ hφ f) ∧ IsPreStep P (pushMap P F φ hφ f)) ∧
      pushMap P F φ hφ f ≫
          pushHom P F φ hφ W.hom.hom W.hom.property.1 W.hom.property.2 =
        pushHom P F φ hφ Z.hom.hom Z.hom.property.1 Z.hom.property.2 :=
  (prop_1_9_ii_hom P F φ
    Z.hom.hom (pushMid P F φ hφ Z.hom.hom Z.hom.property.1 Z.hom.property.2)
      (pushHom P F φ hφ Z.hom.hom Z.hom.property.1 Z.hom.property.2)
    W.hom.hom (pushMid P F φ hφ W.hom.hom W.hom.property.1 W.hom.property.2)
      (pushHom P F φ hφ W.hom.hom W.hom.property.1 W.hom.property.2)
    (pushFac P F φ hφ Z.hom.hom Z.hom.property.1 Z.hom.property.2)
    (pushMid_spec P F φ hφ Z.hom.hom Z.hom.property.1 Z.hom.property.2).1
    (pushMid_spec P F φ hφ Z.hom.hom Z.hom.property.1 Z.hom.property.2).2
    (pushHom_spec P F φ hφ Z.hom.hom Z.hom.property.1 Z.hom.property.2).1
    (pushHom_spec P F φ hφ Z.hom.hom Z.hom.property.1 Z.hom.property.2).2
    (pushFac P F φ hφ W.hom.hom W.hom.property.1 W.hom.property.2)
    (pushMid_spec P F φ hφ W.hom.hom W.hom.property.1 W.hom.property.2).2
    (pushHom_spec P F φ hφ W.hom.hom W.hom.property.1 W.hom.property.2).1
    (pushHom_spec P F φ hφ W.hom.hom W.hom.property.1 W.hom.property.2).2
    f.left.hom f.left.property.2
    (congrArg InducedWideCategory.Hom.hom (Over.w f))).choose_spec

/-- ★★**`φ_* : 𝒞^imtr-pre_A ⥤ 𝒞^imtr-pre_B`**(`φ` は base-isomorphism)。 -/
noncomputable def pushFunctor {A B : C} (φ : A ⟶ B) (hφ : IsBaseIsomorphism P φ) :
    Over (⟨A⟩ : ImtrPre P) ⥤ Over (⟨B⟩ : ImtrPre P) where
  obj Z := Over.mk (⟨pushHom P F φ hφ Z.hom.hom Z.hom.property.1 Z.hom.property.2,
      pushHom_spec P F φ hφ Z.hom.hom Z.hom.property.1 Z.hom.property.2⟩ :
    (⟨pushObj P F φ hφ Z.hom.hom Z.hom.property.1 Z.hom.property.2⟩ : ImtrPre P) ⟶
      (⟨B⟩ : ImtrPre P))
  map {Z W} f := Over.homMk (⟨pushMap P F φ hφ f, (pushMap_spec P F φ hφ f).1⟩ :
      (⟨pushObj P F φ hφ Z.hom.hom Z.hom.property.1 Z.hom.property.2⟩ : ImtrPre P) ⟶
        (⟨pushObj P F φ hφ W.hom.hom W.hom.property.1 W.hom.property.2⟩ : ImtrPre P))
    (InducedWideCategory.Hom.ext (pushMap_spec P F φ hφ f).2)
  map_id _ := imtrPreOver_hom_subsingleton P F _ _
  map_comp _ _ := imtrPreOver_hom_subsingleton P F _ _

/-! ### `φ^*`(pull-back morphism から) -/

noncomputable def pullObj {A B Dd : C} (φ : A ⟶ B) (hpb : IsPullBack P φ)
    (δ : Dd ⟶ B) (hδi : IsIsometric P δ) (hδs : IsPreStep P δ) : C :=
  (prop_1_9_iii_lift P F φ hpb δ hδi hδs).choose

noncomputable def pullHom {A B Dd : C} (φ : A ⟶ B) (hpb : IsPullBack P φ)
    (δ : Dd ⟶ B) (hδi : IsIsometric P δ) (hδs : IsPreStep P δ) :
    pullObj P F φ hpb δ hδi hδs ⟶ A :=
  (prop_1_9_iii_lift P F φ hpb δ hδi hδs).choose_spec.choose

/-- ★四角形の下辺(pull-back)。 -/
noncomputable def pullPsi {A B Dd : C} (φ : A ⟶ B) (hpb : IsPullBack P φ)
    (δ : Dd ⟶ B) (hδi : IsIsometric P δ) (hδs : IsPreStep P δ) :
    pullObj P F φ hpb δ hδi hδs ⟶ Dd :=
  (prop_1_9_iii_lift P F φ hpb δ hδi hδs).choose_spec.choose_spec.choose

theorem pullFac {A B Dd : C} (φ : A ⟶ B) (hpb : IsPullBack P φ)
    (δ : Dd ⟶ B) (hδi : IsIsometric P δ) (hδs : IsPreStep P δ) :
    pullHom P F φ hpb δ hδi hδs ≫ φ = pullPsi P F φ hpb δ hδi hδs ≫ δ :=
  (prop_1_9_iii_lift P F φ hpb δ hδi hδs).choose_spec.choose_spec.choose_spec.1

theorem pullPsi_spec {A B Dd : C} (φ : A ⟶ B) (hpb : IsPullBack P φ)
    (δ : Dd ⟶ B) (hδi : IsIsometric P δ) (hδs : IsPreStep P δ) :
    IsPullBack P (pullPsi P F φ hpb δ hδi hδs) :=
  (prop_1_9_iii_lift P F φ hpb δ hδi hδs).choose_spec.choose_spec.choose_spec.2.1

theorem pullHom_spec {A B Dd : C} (φ : A ⟶ B) (hpb : IsPullBack P φ)
    (δ : Dd ⟶ B) (hδi : IsIsometric P δ) (hδs : IsPreStep P δ) :
    IsIsometric P (pullHom P F φ hpb δ hδi hδs) ∧
      IsPreStep P (pullHom P F φ hpb δ hδi hδs) :=
  ⟨(prop_1_9_iii_lift P F φ hpb δ hδi hδs).choose_spec.choose_spec.choose_spec.2.2.1,
   (prop_1_9_iii_lift P F φ hpb δ hδi hδs).choose_spec.choose_spec.choose_spec.2.2.2⟩

/-- ★`φ^*` の射への割り当て。 -/
noncomputable def pullMap {A B : C} (φ : A ⟶ B) (hpb : IsPullBack P φ)
    {Z W : Over (⟨B⟩ : ImtrPre P)} (f : Z ⟶ W) :
    pullObj P F φ hpb Z.hom.hom Z.hom.property.1 Z.hom.property.2 ⟶
      pullObj P F φ hpb W.hom.hom W.hom.property.1 W.hom.property.2 :=
  (prop_1_9_iii_hom P φ hpb W.hom.hom W.hom.property.2
    (pullHom P F φ hpb Z.hom.hom Z.hom.property.1 Z.hom.property.2)
    (pullPsi P F φ hpb Z.hom.hom Z.hom.property.1 Z.hom.property.2 ≫ f.left.hom)
    (pullHom P F φ hpb W.hom.hom W.hom.property.1 W.hom.property.2)
    (pullPsi P F φ hpb W.hom.hom W.hom.property.1 W.hom.property.2)
    (by
      have hw : (Over.Hom.left f).hom ≫ W.hom.hom = Z.hom.hom :=
        congrArg InducedWideCategory.Hom.hom (Over.w f)
      rw [pullFac P F φ hpb Z.hom.hom Z.hom.property.1 Z.hom.property.2, Category.assoc, hw])
    (pullFac P F φ hpb W.hom.hom W.hom.property.1 W.hom.property.2)
    (pullHom_spec P F φ hpb W.hom.hom W.hom.property.1 W.hom.property.2).2
    (pullPsi_spec P F φ hpb W.hom.hom W.hom.property.1 W.hom.property.2)).choose

theorem pullMap_tri {A B : C} (φ : A ⟶ B) (hpb : IsPullBack P φ)
    {Z W : Over (⟨B⟩ : ImtrPre P)} (f : Z ⟶ W) :
    pullMap P F φ hpb f ≫
        pullHom P F φ hpb W.hom.hom W.hom.property.1 W.hom.property.2 =
      pullHom P F φ hpb Z.hom.hom Z.hom.property.1 Z.hom.property.2 :=
  (prop_1_9_iii_hom P φ hpb W.hom.hom W.hom.property.2
    (pullHom P F φ hpb Z.hom.hom Z.hom.property.1 Z.hom.property.2)
    (pullPsi P F φ hpb Z.hom.hom Z.hom.property.1 Z.hom.property.2 ≫ f.left.hom)
    (pullHom P F φ hpb W.hom.hom W.hom.property.1 W.hom.property.2)
    (pullPsi P F φ hpb W.hom.hom W.hom.property.1 W.hom.property.2)
    (by
      have hw : (Over.Hom.left f).hom ≫ W.hom.hom = Z.hom.hom :=
        congrArg InducedWideCategory.Hom.hom (Over.w f)
      rw [pullFac P F φ hpb Z.hom.hom Z.hom.property.1 Z.hom.property.2, Category.assoc, hw])
    (pullFac P F φ hpb W.hom.hom W.hom.property.1 W.hom.property.2)
    (pullHom_spec P F φ hpb W.hom.hom W.hom.property.1 W.hom.property.2).2
    (pullPsi_spec P F φ hpb W.hom.hom W.hom.property.1 W.hom.property.2)).choose_spec.1

/-- ★`pullMap` が本当に isometric pre-step であることは
**`Proposition 1.7, (v)`**(合成が属せば因子も属する)から出る ——
`prop_1_9_iii_hom` はそこまでは言わない。 -/
theorem pullMap_spec {A B : C} (φ : A ⟶ B) (hpb : IsPullBack P φ)
    {Z W : Over (⟨B⟩ : ImtrPre P)} (f : Z ⟶ W) :
    IsIsometric P (pullMap P F φ hpb f) ∧ IsPreStep P (pullMap P F φ hpb f) := by
  have ht := pullMap_tri P F φ hpb f
  have h1 := (pullHom_spec P F φ hpb Z.hom.hom Z.hom.property.1 Z.hom.property.2).1
  have h2 := (pullHom_spec P F φ hpb Z.hom.hom Z.hom.property.1 Z.hom.property.2).2
  rw [← ht] at h1 h2
  exact ⟨(prop_1_7_v_isometric P _ _ h1).1, (prop_1_7_v_preStep P _ _ h2).1⟩

/-- ★★**`φ^* : 𝒞^imtr-pre_B ⥤ 𝒞^imtr-pre_A`**(`φ` は pull-back morphism)。 -/
noncomputable def pullFunctor {A B : C} (φ : A ⟶ B) (hpb : IsPullBack P φ) :
    Over (⟨B⟩ : ImtrPre P) ⥤ Over (⟨A⟩ : ImtrPre P) where
  obj Z := Over.mk (⟨pullHom P F φ hpb Z.hom.hom Z.hom.property.1 Z.hom.property.2,
      pullHom_spec P F φ hpb Z.hom.hom Z.hom.property.1 Z.hom.property.2⟩ :
    (⟨pullObj P F φ hpb Z.hom.hom Z.hom.property.1 Z.hom.property.2⟩ : ImtrPre P) ⟶
      (⟨A⟩ : ImtrPre P))
  map {Z W} f := Over.homMk (⟨pullMap P F φ hpb f, pullMap_spec P F φ hpb f⟩ :
      (⟨pullObj P F φ hpb Z.hom.hom Z.hom.property.1 Z.hom.property.2⟩ : ImtrPre P) ⟶
        (⟨pullObj P F φ hpb W.hom.hom W.hom.property.1 W.hom.property.2⟩ : ImtrPre P))
    (InducedWideCategory.Hom.ext (pullMap_tri P F φ hpb f))
  map_id _ := imtrPreOver_hom_subsingleton P F _ _
  map_comp _ _ := imtrPreOver_hom_subsingleton P F _ _

/-- **`𝒪^×(A)^imtr-pre ⊆ 𝒪^×(A)`** —— `u_*` が恒等関手と同型になる `u` の全体。

原文 (FrdI p.31):
> the subgroup of v ∈O×(A) for which vimtr-pre is the identity.

★`u ∈ 𝒪^×(A)` は同型なので base-isomorphism、したがって `u_*` が定義できる。 -/
def OTimesImtrPre (A : C) : Set (End A) :=
  {u | u ∈ OTimes P A ∧ ∃ hb : IsBaseIsomorphism P (u : A ⟶ A),
    Nonempty (pushFunctor P F (u : A ⟶ A) hb ≅ 𝟭 (Over (⟨A⟩ : ImtrPre P)))}

theorem otimesImtrPre_subset (A : C) : OTimesImtrPre P F A ⊆ OTimes P A :=
  fun _ h => h.1

/-! ### ★★`𝒪^×(A)^{imtr-pre}` が**部分群**であること

原文 (FrdI p.31):
> the subgroup of v ∈O×(A) for which vimtr-pre is the identity.

★★**原文は「部分群」と述べている。** 上の `Set` だけでは主張を満たさない
(2026-08-16 の監査で判明)。単位元・積・逆元の閉性を示す。

★**楽になる事実**: `Over (⟨A⟩ : ImtrPre P)` の hom は subsingleton
(`imtrPreOver_hom_subsingleton`)。★**したがって自然性は自動**で、
対象ごとの同型さえ作れば関手の同型が得られる。 -/

/-- ★**`ImtrPre P` の同型を、`𝒞` の同型から作る**。

★同型は isometric(`isIsometric_of_isIso`)かつ pre-step(`isPreStep_of_isIso`)。 -/
def imtrPreIsoOfIso {X Y : C} (δ : X ≅ Y) : (⟨X⟩ : ImtrPre P) ≅ (⟨Y⟩ : ImtrPre P) where
  hom := ⟨δ.hom, ⟨isIsometric_of_isIso P δ.hom, isPreStep_of_isIso P δ.hom⟩⟩
  inv := ⟨δ.inv, ⟨isIsometric_of_isIso P δ.inv, isPreStep_of_isIso P δ.inv⟩⟩
  hom_inv_id := InducedWideCategory.Hom.ext (by simp)
  inv_hom_id := InducedWideCategory.Hom.ext (by simp)

include F in
/-- ★★**`(𝟙 A)_* ≅ 𝟭`** —— `𝒪^×(A)^{imtr-pre}` が**単位元を含む**ことの中身。

★`ε` 自身が isometric pre-step なので、`ε = 𝟙 ≫ ε` が第2の分解を与える。
★**分解の一意性**(`prop_1_9_i_uniq`)が対象ごとの同型を出し、
自然性は hom の subsingleton 性から自動。 -/
theorem pushId_uniq (A : C) (hb : IsBaseIsomorphism P (𝟙 A))
    (Z : Over (⟨A⟩ : ImtrPre P)) :
    ∃ δ : pushObj P F (𝟙 A) hb Z.hom.hom Z.hom.property.1 Z.hom.property.2 ≅ Z.left.obj,
      Z.hom.hom
        = δ.inv ≫ pushHom P F (𝟙 A) hb Z.hom.hom Z.hom.property.1 Z.hom.property.2 ∧
      (𝟙 Z.left.obj : Z.left.obj ⟶ Z.left.obj)
        = pushMid P F (𝟙 A) hb Z.hom.hom Z.hom.property.1 Z.hom.property.2 ≫ δ.hom :=
  prop_1_9_i_uniq P F _ _ _ _ (𝟙 Z.left.obj) Z.hom.hom
    (by rw [← pushFac P F (𝟙 A) hb Z.hom.hom Z.hom.property.1 Z.hom.property.2]; simp)
    (pushMid_spec P F (𝟙 A) hb Z.hom.hom Z.hom.property.1 Z.hom.property.2).1
    (pushMid_spec P F (𝟙 A) hb Z.hom.hom Z.hom.property.1 Z.hom.property.2).2
    (pushHom_spec P F (𝟙 A) hb Z.hom.hom Z.hom.property.1 Z.hom.property.2).1
    (pushHom_spec P F (𝟙 A) hb Z.hom.hom Z.hom.property.1 Z.hom.property.2).2
    (isCoAngular_id P _) (isBaseIsomorphism_of_isIso P _)
    Z.hom.property.1 Z.hom.property.2

include F in
/-- ★★**`(𝟙 A)_* ≅ 𝟭`**。★`∃` からデータを取るので `choose` を経由する。 -/
noncomputable def pushFunctorIdIso (A : C) (hb : IsBaseIsomorphism P (𝟙 A)) :
    pushFunctor P F (𝟙 A) hb ≅ 𝟭 (Over (⟨A⟩ : ImtrPre P)) :=
  NatIso.ofComponents
    (fun Z =>
      Over.isoMk (imtrPreIsoOfIso P (pushId_uniq P F A hb Z).choose)
        (InducedWideCategory.Hom.ext (by
          -- ★`rw` は motive が壊れる(`pushHom` の**証明引数**が `Z.hom.hom` に依存する)。
          -- ★`congrArg` なら関数が明示なので依存が起きない。
          have h2 := congrArg (fun t => (pushId_uniq P F A hb Z).choose.hom ≫ t)
            (pushId_uniq P F A hb Z).choose_spec.1
          simp only [Iso.hom_inv_id_assoc] at h2
          -- ★`simpa` ではなく `exact`(既定透明度でないと 2 つの綴りが同一視されない、表 #1)
          exact h2)))
    (fun _ => imtrPreOver_hom_subsingleton P F _ _)

include F in
/-- ★★**`𝒪^×(A)^{imtr-pre}` は単位元を含む**(部分群性の 3 条件のうち 1 つ目)。 -/
theorem one_mem_otimesImtrPre (A : C) : (1 : End A) ∈ OTimesImtrPre P F A :=
  ⟨(OTimes P A).one_mem, isBaseIsomorphism_of_isIso P (𝟙 A),
   ⟨pushFunctorIdIso P F A (isBaseIsomorphism_of_isIso P (𝟙 A))⟩⟩

include F in
/-- ★★**合成の分解** —— `ε ≫ (φ ≫ ψ)` の**第2の分解**。

★左因子は co-angular base-isomorphism の**合成**(`F.coAngularComp` と
`isBaseIsomorphism_comp`)、右因子は isometric pre-step。 -/
theorem pushComp_uniq {A B E : C} (φ : A ⟶ B) (hφ : IsBaseIsomorphism P φ)
    (ψ : B ⟶ E) (hψ : IsBaseIsomorphism P ψ) (Z : Over (⟨A⟩ : ImtrPre P)) :
    ∃ δ : pushObj P F (φ ≫ ψ) (isBaseIsomorphism_comp P hφ hψ)
            Z.hom.hom Z.hom.property.1 Z.hom.property.2
        ≅ pushObj P F ψ hψ
            (pushHom P F φ hφ Z.hom.hom Z.hom.property.1 Z.hom.property.2)
            (pushHom_spec P F φ hφ Z.hom.hom Z.hom.property.1 Z.hom.property.2).1
            (pushHom_spec P F φ hφ Z.hom.hom Z.hom.property.1 Z.hom.property.2).2,
      pushHom P F ψ hψ
          (pushHom P F φ hφ Z.hom.hom Z.hom.property.1 Z.hom.property.2)
          (pushHom_spec P F φ hφ Z.hom.hom Z.hom.property.1 Z.hom.property.2).1
          (pushHom_spec P F φ hφ Z.hom.hom Z.hom.property.1 Z.hom.property.2).2
        = δ.inv ≫ pushHom P F (φ ≫ ψ) (isBaseIsomorphism_comp P hφ hψ)
            Z.hom.hom Z.hom.property.1 Z.hom.property.2 ∧
      (pushMid P F φ hφ Z.hom.hom Z.hom.property.1 Z.hom.property.2
        ≫ pushMid P F ψ hψ
            (pushHom P F φ hφ Z.hom.hom Z.hom.property.1 Z.hom.property.2)
            (pushHom_spec P F φ hφ Z.hom.hom Z.hom.property.1 Z.hom.property.2).1
            (pushHom_spec P F φ hφ Z.hom.hom Z.hom.property.1 Z.hom.property.2).2)
        = pushMid P F (φ ≫ ψ) (isBaseIsomorphism_comp P hφ hψ)
            Z.hom.hom Z.hom.property.1 Z.hom.property.2 ≫ δ.hom := by
  set ε := Z.hom.hom with hε
  set hεi := Z.hom.property.1
  set hεs := Z.hom.property.2
  set h1 := (pushHom_spec P F φ hφ ε hεi hεs).1
  set h2 := (pushHom_spec P F φ hφ ε hεi hεs).2
  refine prop_1_9_i_uniq P F _ _ _ _
    (pushMid P F φ hφ ε hεi hεs ≫ pushMid P F ψ hψ (pushHom P F φ hφ ε hεi hεs) h1 h2)
    (pushHom P F ψ hψ (pushHom P F φ hφ ε hεi hεs) h1 h2) ?_ ?_ ?_ ?_ ?_ ?_ ?_ ?_ ?_
  · rw [← pushFac P F (φ ≫ ψ) (isBaseIsomorphism_comp P hφ hψ) ε hεi hεs,
      Category.assoc, ← pushFac P F ψ hψ (pushHom P F φ hφ ε hεi hεs) h1 h2,
      ← Category.assoc, pushFac P F φ hφ ε hεi hεs, Category.assoc]
  · exact (pushMid_spec P F (φ ≫ ψ) (isBaseIsomorphism_comp P hφ hψ) ε hεi hεs).1
  · exact (pushMid_spec P F (φ ≫ ψ) (isBaseIsomorphism_comp P hφ hψ) ε hεi hεs).2
  · exact (pushHom_spec P F (φ ≫ ψ) (isBaseIsomorphism_comp P hφ hψ) ε hεi hεs).1
  · exact (pushHom_spec P F (φ ≫ ψ) (isBaseIsomorphism_comp P hφ hψ) ε hεi hεs).2
  · exact F.coAngularComp _ _ (pushMid_spec P F φ hφ ε hεi hεs).1
      (pushMid_spec P F ψ hψ (pushHom P F φ hφ ε hεi hεs) h1 h2).1
  · exact isBaseIsomorphism_comp P (pushMid_spec P F φ hφ ε hεi hεs).2
      (pushMid_spec P F ψ hψ (pushHom P F φ hφ ε hεi hεs) h1 h2).2
  · exact (pushHom_spec P F ψ hψ (pushHom P F φ hφ ε hεi hεs) h1 h2).1
  · exact (pushHom_spec P F ψ hψ (pushHom P F φ hφ ε hεi hεs) h1 h2).2

include F in
/-- ★★**`(φ ≫ ψ)_* ≅ φ_* ⋙ ψ_*`**。★自然性は hom の subsingleton 性から自動。 -/
noncomputable def pushFunctorCompIso {A B E : C} (φ : A ⟶ B) (hφ : IsBaseIsomorphism P φ)
    (ψ : B ⟶ E) (hψ : IsBaseIsomorphism P ψ) :
    pushFunctor P F (φ ≫ ψ) (isBaseIsomorphism_comp P hφ hψ)
      ≅ pushFunctor P F φ hφ ⋙ pushFunctor P F ψ hψ :=
  NatIso.ofComponents
    (fun Z =>
      Over.isoMk (imtrPreIsoOfIso P (pushComp_uniq P F φ hφ ψ hψ Z).choose)
        (InducedWideCategory.Hom.ext (by
          have h2 := congrArg (fun t => (pushComp_uniq P F φ hφ ψ hψ Z).choose.hom ≫ t)
            (pushComp_uniq P F φ hφ ψ hψ Z).choose_spec.1
          simp only [Iso.hom_inv_id_assoc] at h2
          exact h2)))
    (fun _ => imtrPreOver_hom_subsingleton P F _ _)

include F in
/-- ★★**積で閉じる**(部分群性の 2 つ目)。

★**`End A` の積は `x * y = y ≫ x`** である(合成の向きに注意)。 -/
theorem mul_mem_otimesImtrPre (A : C) {u v : End A}
    (hu : u ∈ OTimesImtrPre P F A) (hv : v ∈ OTimesImtrPre P F A) :
    u * v ∈ OTimesImtrPre P F A := by
  obtain ⟨huo, hbu, ⟨eu⟩⟩ := hu
  obtain ⟨hvo, hbv, ⟨ev⟩⟩ := hv
  refine ⟨(OTimes P A).mul_mem huo hvo, isBaseIsomorphism_comp P hbv hbu, ⟨?_⟩⟩
  exact pushFunctorCompIso P F _ hbv _ hbu ≪≫
    Functor.isoWhiskerRight ev (pushFunctor P F (u : A ⟶ A) hbu) ≪≫
    Functor.leftUnitor _ ≪≫ eu

include F in
/-- ★★**逆元で閉じる**(部分群性の 3 つ目)。

★`u * w = 1`(すなわち `w ≫ u = 𝟙`)なる `w` が `𝒪^×(A)` にあれば、
`u ∈ ⟹ w ∈`。★**積と単位元の結果から出る。** -/
theorem inv_mem_otimesImtrPre (A : C) {u w : End A}
    (hu : u ∈ OTimesImtrPre P F A) (hwo : w ∈ OTimes P A) (h : u * w = 1) :
    w ∈ OTimesImtrPre P F A := by
  obtain ⟨-, hbu, ⟨eu⟩⟩ := hu
  -- `w` が base-isomorphism であること —— `𝒪^×` の元は同型
  haveI hwi : IsIso (w : A ⟶ A) := (CategoryTheory.isUnit_iff_isIso (w : End A)).mp hwo.2
  have hbw : IsBaseIsomorphism P (w : A ⟶ A) := isBaseIsomorphism_of_isIso P _
  refine ⟨hwo, hbw, ⟨?_⟩⟩
  -- `w ≫ u = 𝟙` から `pushFunctor w ⋙ pushFunctor u ≅ 𝟭`
  -- ★`End A` の積は `x * y = y ≫ x` なので、`u * w = 1` は `w ≫ u = 𝟙` である
  have hcomp : ((w : A ⟶ A) ≫ (u : A ⟶ A)) = 𝟙 A := h
  have hEq : pushFunctor P F ((w : A ⟶ A) ≫ (u : A ⟶ A)) (isBaseIsomorphism_comp P hbw hbu)
      = pushFunctor P F (𝟙 A) (isBaseIsomorphism_of_isIso P (𝟙 A)) := by
    -- ★`rw` は motive が壊れる(証明引数が射に依存)。`congr` なら証明無関係が効く。
    congr 1
  have e1 : pushFunctor P F (w : A ⟶ A) hbw ⋙ pushFunctor P F (u : A ⟶ A) hbu
      ≅ 𝟭 (Over (⟨A⟩ : ImtrPre P)) :=
    (pushFunctorCompIso P F _ hbw _ hbu).symm ≪≫ eqToIso hEq ≪≫
      pushFunctorIdIso P F A (isBaseIsomorphism_of_isIso P (𝟙 A))
  exact (Functor.rightUnitor _).symm ≪≫ Functor.isoWhiskerLeft _ eu.symm ≪≫ e1

/-! ### ★★`u^{imtr-pre} = id` の**十分条件**

★**`u` が isometric pre-step と「同型を除いて可換」なら `u_* ≅ 𝟭`。**

★これは `Proposition 1.13, (iii)` で使う ——
`𝒞_A → 𝒞` の自己同型 `α` は、自然性から
`γ ≫ α_X = α_Y ≫ γ`(`γ` は isometric pre-step)を満たし、
`α_Y` は同型だから、まさにこの形になる。 -/

include F in
/-- ★`u` と `γ` が同型 `v` を挟んで可換なら、`γ ≫ u` の第2の分解が `(v, γ)` である。 -/
theorem pushComm_uniq {A : C} (u : A ⟶ A) (hb : IsBaseIsomorphism P u)
    (Z : Over (⟨A⟩ : ImtrPre P)) (v : Z.left.obj ⟶ Z.left.obj) [IsIso v]
    (hcomm : Z.hom.hom ≫ u = v ≫ Z.hom.hom) :
    ∃ δ : pushObj P F u hb Z.hom.hom Z.hom.property.1 Z.hom.property.2 ≅ Z.left.obj,
      Z.hom.hom = δ.inv ≫ pushHom P F u hb Z.hom.hom Z.hom.property.1 Z.hom.property.2 ∧
      v = pushMid P F u hb Z.hom.hom Z.hom.property.1 Z.hom.property.2 ≫ δ.hom :=
  prop_1_9_i_uniq P F _ _ _ _ v Z.hom.hom
    (by rw [← pushFac P F u hb Z.hom.hom Z.hom.property.1 Z.hom.property.2]; exact hcomm)
    (pushMid_spec P F u hb Z.hom.hom Z.hom.property.1 Z.hom.property.2).1
    (pushMid_spec P F u hb Z.hom.hom Z.hom.property.1 Z.hom.property.2).2
    (pushHom_spec P F u hb Z.hom.hom Z.hom.property.1 Z.hom.property.2).1
    (pushHom_spec P F u hb Z.hom.hom Z.hom.property.1 Z.hom.property.2).2
    (isCoAngular_of_isIso P v) (isBaseIsomorphism_of_isIso P v)
    Z.hom.property.1 Z.hom.property.2

include F in
/-- ★★**`u` が isometric pre-step と同型を除いて可換なら `u_* ≅ 𝟭`**。

★`pushFunctorIdIso` の一般化(`v = 𝟙` がその場合)。自然性は
`imtrPreOver_hom_subsingleton` から自動。 -/
noncomputable def pushFunctorIsoOfCommutes {A : C} (u : A ⟶ A) (hb : IsBaseIsomorphism P u)
    (comm : ∀ Z : Over (⟨A⟩ : ImtrPre P),
      Σ' v : Z.left.obj ⟶ Z.left.obj, PLift (IsIso v) ×' Z.hom.hom ≫ u = v ≫ Z.hom.hom) :
    pushFunctor P F u hb ≅ 𝟭 (Over (⟨A⟩ : ImtrPre P)) :=
  NatIso.ofComponents
    (fun Z =>
      haveI : IsIso (comm Z).1 := (comm Z).2.1.down
      Over.isoMk (imtrPreIsoOfIso P (pushComm_uniq P F u hb Z (comm Z).1 (comm Z).2.2).choose)
        (InducedWideCategory.Hom.ext (by
          have h2 := congrArg
            (fun t => (pushComm_uniq P F u hb Z (comm Z).1 (comm Z).2.2).choose.hom ≫ t)
            (pushComm_uniq P F u hb Z (comm Z).1 (comm Z).2.2).choose_spec.1
          simp only [Iso.hom_inv_id_assoc] at h2
          exact h2)))
    (fun _ => imtrPreOver_hom_subsingleton P F _ _)

include F in
/-- ★★**`𝒪^×(A)^{imtr-pre}` を部分モノイドとして包む**。

原文 (FrdI p.31):
> the subgroup of v ∈O×(A) for which vimtr-pre is the identity.

★**`𝒪^×(A)` 自体が `End A` の部分モノイドとして実装されている**ので、
ここも部分モノイドとして構成し、★**逆元の閉性は `inv_mem_otimesImtrPre` で別に述べる**。
`𝒪^×` の元はすべて可逆なので、この 3 条件が原文の「subgroup」の内容である。 -/
def OTimesImtrPreSubmonoid (A : C) : Submonoid (End A) where
  carrier := OTimesImtrPre P F A
  one_mem' := one_mem_otimesImtrPre P F A
  mul_mem' hx hy := mul_mem_otimesImtrPre P F A hx hy

include F in
/-- ★`𝒪^×(A)^{imtr-pre} ≤ 𝒪^×(A)` —— 原文の包含。 -/
theorem otimesImtrPreSubmonoid_le (A : C) :
    OTimesImtrPreSubmonoid P F A ≤ OTimes P A :=
  fun _ h => h.1

end ImtrPreFunctor

/-! ## ★第5段 —— (ii) の圏同値

原文 (FrdI p.31):
> C →A with φ : A →B. Moreover, if φ is a co-angular pre-step, then φ∗is an

★中身は原文 p.32–33 の 6 段である。まず `(Φ(v))⁻¹(Div v)` の計算規則を用意する。
-/


def pullFunctor.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 31, item := "Proposition 1.9, (iii)",
    sectionId := "frdi-prop-1-9-iii" }

end ABC3.Found.FrdI
