/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop18
import ABC3.Found.FrdI.Prop19.ImtrPreFunctor

/-!
# Prop19 —— 押し出しの同値・負の対照・出典

☆もとの 1 枚を**入れ子の切れ目**で割ったものである(第 1457)。
-/

namespace ABC3.Found.FrdI

open CategoryTheory Opposite

universe v u w u2 v2

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ)


section PushEquiv

/-- ★**`(Φ(v))⁻¹(Div v) ∈ Φ(B)`** —— `Definition 1.3, (iii), (d)` のスライス側の不変量。

`coaPreOverFunctor` が対象に割り当てているものと**同じ式**である。 -/
noncomputable def divInv {X B : C} (v : X ⟶ B) (hv : IsBaseIsomorphism P v) :
    Φ.val (P.toElem.obj B).base :=
  haveI : IsIso (P.Base v) := hv
  Φ.map (inv (P.Base v)) (P.Div v)

/-- ★不変量は **isometric な射を左から合成しても変わらない**。

原文 (FrdI p.33) が「since `Div(β◦ψ) = Div(φ′◦α′)`, and `β`, `α′` are isometric」
と言って使うのがこれである。 -/
theorem divInv_comp_left {Y X B : C} (u : Y ⟶ X) (v : X ⟶ B)
    (hu : IsBaseIsomorphism P u) (hv : IsBaseIsomorphism P v) (hui : IsIsometric P u) :
    divInv P (u ≫ v) (isBaseIsomorphism_comp P hu hv) = divInv P v hv := by
  haveI hiu : IsIso (P.Base u) := hu
  haveI hiv : IsIso (P.Base v) := hv
  haveI : IsIso (P.Base (u ≫ v)) := isBaseIsomorphism_comp P hu hv
  show Φ.map (inv (P.Base (u ≫ v))) (P.Div (u ≫ v)) = Φ.map (inv (P.Base v)) (P.Div v)
  have hd : P.Div (u ≫ v) = Φ.map (P.Base u) (P.Div v) := by
    rw [P.Div_comp, show P.Div u = 0 from hui, smul_zero, add_zero]
  have hinv : inv (P.Base (u ≫ v)) = inv (P.Base v) ≫ inv (P.Base u) := by
    refine IsIso.inv_eq_of_hom_inv_id ?_
    rw [P.Base_comp]
    simp
  rw [hd, hinv, Φ.map_comp (inv (P.Base u)) (inv (P.Base v)),
    ← Φ.map_comp (P.Base u) (inv (P.Base u)), IsIso.inv_hom_id, Φ.map_id]

/-- ★不変量に **isometric な pre-step を右から合成する**と、`Φ` で引き戻される。 -/
theorem divInv_comp_right {Y X B : C} (u : Y ⟶ X) (v : X ⟶ B)
    (hu : IsBaseIsomorphism P u) (hv : IsBaseIsomorphism P v) (hvi : IsIsometric P v)
    (hvl : IsLinear P v) :
    divInv P (u ≫ v) (isBaseIsomorphism_comp P hu hv) =
      haveI : IsIso (P.Base v) := hv
      Φ.map (inv (P.Base v)) (divInv P u hu) := by
  haveI hiu : IsIso (P.Base u) := hu
  haveI hiv : IsIso (P.Base v) := hv
  haveI : IsIso (P.Base (u ≫ v)) := isBaseIsomorphism_comp P hu hv
  show Φ.map (inv (P.Base (u ≫ v))) (P.Div (u ≫ v))
    = Φ.map (inv (P.Base v)) (Φ.map (inv (P.Base u)) (P.Div u))
  have hd : P.Div (u ≫ v) = P.Div u := by
    rw [P.Div_comp, show P.Div v = 0 from hvi, map_zero, zero_add,
      show P.degFr v = 1 from hvl]
    simp
  have hinv : inv (P.Base (u ≫ v)) = inv (P.Base v) ≫ inv (P.Base u) := by
    refine IsIso.inv_eq_of_hom_inv_id ?_
    rw [P.Base_comp]
    simp
  rw [hd, hinv, Φ.map_comp (inv (P.Base u)) (inv (P.Base v))]

/-- ★`Φ` が divisorial なので、`Order(Φ(B))` の同型は**等式**に落ちる。

★`mle_antisymm` の 2 仮定(integral / sharp)は**どちらも `divisorial` に含まれる**。 -/
theorem eq_of_orderCat_iso {B : C} {a b : Φ.val (P.toElem.obj B).base}
    (e : toOrderCat a ≅ toOrderCat b) : a = b :=
  mle_antisymm (P.divisorial _).1.1 (P.divisorial _).2 (leOfHom e.hom) (leOfHom e.inv)

/-- ★射が等しければ不変量も等しい(証明の取り方には依らない)。 -/
theorem divInv_congr {X B : C} {v v' : X ⟶ B} (h : v = v')
    (hv : IsBaseIsomorphism P v) (hv' : IsBaseIsomorphism P v') :
    divInv P v hv = divInv P v' hv' := by
  subst h; rfl

end PushEquiv

section PushEquiv2

variable [(coaPreProp P).IsMultiplicative]

/-- ★★**同じ不変量を持つ2つの co-angular pre-step は同型で移り合う**(スライス側)。

★`coaPre_iso_of_div_eq`(コスライス側、`Proposition 1.8` で作った)の**スライス版**。
`Definition 1.3, (iii), (d)` の**第2**の圏同値の充満性と忠実性を使う。 -/
theorem coaPreOver_iso_of_divInv_eq
    (hover : ∀ X : C, (coaPreOverFunctor P X).IsEquivalence)
    {B : C} (Z W : Over (⟨B⟩ : WideSubcategory (coaPreProp P)))
    (h : divInv P Z.hom.hom Z.hom.property.2.2 = divInv P W.hom.hom W.hom.property.2.2) :
    Nonempty (Z ≅ W) := by
  haveI := (hover B).full
  haveI := (hover B).faithful
  have hobj : (coaPreOverFunctor P B).obj Z = (coaPreOverFunctor P B).obj W :=
    congrArg (fun x => op (toOrderCat x)) h
  exact ⟨(coaPreOverFunctor P B).preimageIso (eqToIso hobj)⟩

/-- ★★**原文 p.32–33 の構成**。

`φ : A ⟶ B` が co-angular pre-step、`β : D ⟶ B` が isometric pre-step のとき、
co-angular pre-step `ψ : Cc ⟶ D` と isometric pre-step `α : Cc ⟶ A` で
`ψ ≫ β = α ≫ φ`(原文の `β ◦ ψ = φ ◦ α`)となるものが取れる。

★★これが (ii) の**本質的全射性**であり、原文が言うとおり
**`φ` を `ψ` に取り替えれば充満性**にもなる。

★使う入力: (iii)(d) の**第2**の圏同値(本質的全射性・充満性・忠実性)、
`Definition 1.3, (v), (c)`(`preStepFactor'`)、そして `Φ` の **divisorial 性**
(`Order(Φ(B))` の反対称性)。 -/
theorem push_square (F : FrobenioidCore P)
    (hover : ∀ X : C, (coaPreOverFunctor P X).IsEquivalence)
    {A B Dd : C} (φ : A ⟶ B) (hφc : IsCoAngular P φ) (hφs : IsPreStep P φ)
    (β : Dd ⟶ B) (hβi : IsIsometric P β) (hβs : IsPreStep P β) :
    ∃ (Cc : C) (ψ : Cc ⟶ Dd) (α : Cc ⟶ A),
      ψ ≫ β = α ≫ φ ∧ (IsCoAngular P ψ ∧ IsPreStep P ψ) ∧
        (IsIsometric P α ∧ IsPreStep P α) := by
  haveI := (hover Dd).essSurj
  haveI hbβ : IsIso (P.Base β) := hβs.2
  -- 1. 目標の不変量 `e ∈ Φ(D)`
  obtain ⟨e, he⟩ : ∃ e : Φ.val (P.toElem.obj Dd).base,
      e = Φ.map (P.Base β) (divInv P φ hφs.2) := ⟨_, rfl⟩
  -- 2. (iii)(d) の**本質的全射性**で co-angular pre-step `ψ` を取る
  obtain ⟨Z, hZ⟩ : ∃ Z : Over (⟨Dd⟩ : WideSubcategory (coaPreProp P)),
      Z = (coaPreOverFunctor P Dd).objPreimage (op (toOrderCat e)) := ⟨_, rfl⟩
  have hiZ : (coaPreOverFunctor P Dd).obj Z ≅ op (toOrderCat e) := by
    rw [hZ]; exact (coaPreOverFunctor P Dd).objObjPreimageIso _
  have hobj : (coaPreOverFunctor P Dd).obj Z
      = op (toOrderCat (divInv P Z.hom.hom Z.hom.property.2.2)) := rfl
  rw [hobj] at hiZ
  have hdivψ : divInv P Z.hom.hom Z.hom.property.2.2 = e :=
    (eq_of_orderCat_iso P hiZ.unop).symm
  -- 3. `ψ ≫ β` の不変量は `φ` のそれに等しい
  have hψβs : IsPreStep P (Z.hom.hom ≫ β) := IsPreStep.comp P Z.hom.property.2 hβs
  have hkey : divInv P (Z.hom.hom ≫ β)
      (isBaseIsomorphism_comp P Z.hom.property.2.2 hβs.2) = divInv P φ hφs.2 := by
    rw [divInv_comp_right P Z.hom.hom β Z.hom.property.2.2 hβs.2 hβi hβs.1, hdivψ, he,
      ← Φ.map_comp (P.Base β) (inv (P.Base β)), IsIso.inv_hom_id, Φ.map_id]
  -- 4. `Definition 1.3, (v), (c)` で `ψ ≫ β = α' ≫ φ'`
  obtain ⟨Xx, α', φ', hfac, hα'i, hα's, hφ'c, hφ's⟩ := F.preStepFactor' (Z.hom.hom ≫ β) hψβs
  -- 5. `φ'` の不変量も `φ` のそれに等しい
  have hdiv' : divInv P φ' hφ's.2 = divInv P φ hφs.2 := by
    rw [← divInv_comp_left P α' φ' hα's.2 hφ's.2 hα'i,
      divInv_congr P hfac.symm (isBaseIsomorphism_comp P hα's.2 hφ's.2)
        (isBaseIsomorphism_comp P Z.hom.property.2.2 hβs.2)]
    exact hkey
  -- 6. (iii)(d) の**充満性・忠実性**で同型 `γ : Xx ≅ A`
  obtain ⟨iso⟩ := coaPreOver_iso_of_divInv_eq P hover
    (Over.mk (⟨φ', hφ'c, hφ's⟩ :
      (⟨Xx⟩ : WideSubcategory (coaPreProp P)) ⟶ (⟨B⟩ : WideSubcategory (coaPreProp P))))
    (Over.mk (⟨φ, hφc, hφs⟩ :
      (⟨A⟩ : WideSubcategory (coaPreProp P)) ⟶ (⟨B⟩ : WideSubcategory (coaPreProp P))))
    hdiv'
  obtain ⟨γ, hγ⟩ : ∃ γ : Xx ⟶ A, γ = (Over.Hom.left iso.hom).hom := ⟨_, rfl⟩
  obtain ⟨γ', hγ'⟩ : ∃ γ' : A ⟶ Xx, γ' = (Over.Hom.left iso.inv).hom := ⟨_, rfl⟩
  haveI hγiso : IsIso γ := by
    refine ⟨γ', ?_, ?_⟩
    · rw [hγ, hγ']
      exact congrArg (fun t => InducedWideCategory.Hom.hom (Over.Hom.left t)) iso.hom_inv_id
    · rw [hγ, hγ']
      exact congrArg (fun t => InducedWideCategory.Hom.hom (Over.Hom.left t)) iso.inv_hom_id
  have hγφ : γ ≫ φ = φ' := by
    rw [hγ]
    exact congrArg InducedWideCategory.Hom.hom (Over.w iso.hom)
  -- 7. `α := α' ≫ γ`
  refine ⟨Z.left.obj, Z.hom.hom, α' ≫ γ, ?_, Z.hom.property, ?_, ?_⟩
  · rw [Category.assoc, hγφ, hfac]
  · exact IsIsometric.comp P hα'i (isIsometric_of_isIso P γ)
  · exact IsPreStep.comp P hα's (isPreStep_of_isIso P γ)

/-- 補助: スライスの同型から `𝒞` の同型を取り出す。 -/
theorem coaPreOver_hom_of_iso {B : C} {Z W : Over (⟨B⟩ : WideSubcategory (coaPreProp P))}
    (iso : Z ≅ W) : ∃ γ : Z.left.obj ⟶ W.left.obj, IsIso γ ∧ γ ≫ W.hom.hom = Z.hom.hom := by
  refine ⟨(Over.Hom.left iso.hom).hom, ⟨(Over.Hom.left iso.inv).hom, ?_, ?_⟩, ?_⟩
  · exact congrArg (fun t => InducedWideCategory.Hom.hom (Over.Hom.left t)) iso.hom_inv_id
  · exact congrArg (fun t => InducedWideCategory.Hom.hom (Over.Hom.left t)) iso.inv_hom_id
  · exact congrArg InducedWideCategory.Hom.hom (Over.w iso.hom)

/-- ★(i) の分解の左因子は **pre-step でもある**(`Proposition 1.7, (v)`)。 -/
theorem pushMid_isPreStep (F : FrobenioidCore P) {A B : C} (φ : A ⟶ B)
    (hφ : IsBaseIsomorphism P φ) (hφl : IsLinear P φ) {Cc : C} (ε : Cc ⟶ A)
    (hεi : IsIsometric P ε) (hεs : IsPreStep P ε) :
    IsPreStep P (pushMid P F φ hφ ε hεi hεs) := by
  refine ⟨?_, (pushMid_spec P F φ hφ ε hεi hεs).2⟩
  have hlin : IsLinear P (ε ≫ φ) := IsLinear.comp P hεs.1 hφl
  rw [pushFac P F φ hφ ε hεi hεs] at hlin
  exact (prop_1_7_v_linear P _ _ hlin).1

/-- ★★★**`Proposition 1.9, (ii)` の圏同値** —— `φ` が co-angular pre-step なら
`φ_* : 𝒞^imtr-pre_A ⥤ 𝒞^imtr-pre_B` は圏同値。

★**忠実性は無料**(`𝒞^imtr-pre_A` が thin)。
★**本質的全射性**は `push_square` そのもの。
★**充満性**は `push_square` を **`φ` の代わりに `pushMid W`、`β` の代わりに `h`** に当てる ——
原文の「by possibly replacing φ by ψ」がこれである。 -/
theorem pushFunctor_isEquivalence (F : FrobenioidCore P)
    (hover : ∀ X : C, (coaPreOverFunctor P X).IsEquivalence)
    {A B : C} (φ : A ⟶ B) (hφc : IsCoAngular P φ) (hφs : IsPreStep P φ) :
    (pushFunctor P F φ hφs.2).IsEquivalence := by
  haveI hmono : Mono φ := F.preStepMono φ hφs
  haveI hfaith : (pushFunctor P F φ hφs.2).Faithful :=
    ⟨fun _ => imtrPreOver_hom_subsingleton P F _ _⟩
  -- ★★充満性
  haveI hfull : (pushFunctor P F φ hφs.2).Full := by
    constructor
    intro Z W h
    -- 記号
    obtain ⟨hh, hhe⟩ : ∃ hh : pushObj P F φ hφs.2 Z.hom.hom Z.hom.property.1 Z.hom.property.2 ⟶
        pushObj P F φ hφs.2 W.hom.hom W.hom.property.1 W.hom.property.2,
        hh = (Over.Hom.left h).hom := ⟨_, rfl⟩
    have hhp : IsIsometric P hh ∧ IsPreStep P hh := by
      rw [hhe]; exact (Over.Hom.left h).property
    have hhw : hh ≫ pushHom P F φ hφs.2 W.hom.hom W.hom.property.1 W.hom.property.2 =
        pushHom P F φ hφs.2 Z.hom.hom Z.hom.property.1 Z.hom.property.2 := by
      rw [hhe]; exact congrArg InducedWideCategory.Hom.hom (Over.w h)
    have hmW := pushMid_spec P F φ hφs.2 W.hom.hom W.hom.property.1 W.hom.property.2
    have hmZ := pushMid_spec P F φ hφs.2 Z.hom.hom Z.hom.property.1 Z.hom.property.2
    have hmWs := pushMid_isPreStep P F φ hφs.2 hφs.1 W.hom.hom W.hom.property.1 W.hom.property.2
    have hmZs := pushMid_isPreStep P F φ hφs.2 hφs.1 Z.hom.hom Z.hom.property.1 Z.hom.property.2
    have hHW := pushHom_spec P F φ hφs.2 W.hom.hom W.hom.property.1 W.hom.property.2
    -- ★`push_square` を `φ := pushMid W`, `β := hh` に当てる
    obtain ⟨Cc', ψ', a, heq', hψ', ha⟩ := push_square P F hover
      (pushMid P F φ hφs.2 W.hom.hom W.hom.property.1 W.hom.property.2) hmW.1 hmWs
      hh hhp.1 hhp.2
    -- ★不変量の一致
    have d1 : divInv P (ψ' ≫ hh) (isBaseIsomorphism_comp P hψ'.2.2 hhp.2.2)
        = divInv P (pushMid P F φ hφs.2 W.hom.hom W.hom.property.1 W.hom.property.2) hmW.2 := by
      rw [divInv_congr P heq' (isBaseIsomorphism_comp P hψ'.2.2 hhp.2.2)
        (isBaseIsomorphism_comp P ha.2.2 hmW.2)]
      exact divInv_comp_left P a _ ha.2.2 hmW.2 ha.1
    have d4 : divInv P (pushMid P F φ hφs.2 Z.hom.hom Z.hom.property.1 Z.hom.property.2 ≫ hh)
        (isBaseIsomorphism_comp P hmZ.2 hhp.2.2)
        = divInv P (pushMid P F φ hφs.2 W.hom.hom W.hom.property.1 W.hom.property.2) hmW.2 := by
      haveI : IsIso (P.Base (pushHom P F φ hφs.2 W.hom.hom
        W.hom.property.1 W.hom.property.2)) := hHW.2.2
      refine Φ.map_injective (inv (P.Base (pushHom P F φ hφs.2 W.hom.hom
        W.hom.property.1 W.hom.property.2))) ?_
      have e1 := divInv_comp_right P
        (pushMid P F φ hφs.2 Z.hom.hom Z.hom.property.1 Z.hom.property.2 ≫ hh)
        (pushHom P F φ hφs.2 W.hom.hom W.hom.property.1 W.hom.property.2)
        (isBaseIsomorphism_comp P hmZ.2 hhp.2.2) hHW.2.2 hHW.1 hHW.2.1
      have e2 := divInv_comp_right P
        (pushMid P F φ hφs.2 W.hom.hom W.hom.property.1 W.hom.property.2)
        (pushHom P F φ hφs.2 W.hom.hom W.hom.property.1 W.hom.property.2)
        hmW.2 hHW.2.2 hHW.1 hHW.2.1
      rw [← e1, ← e2]
      have f1 : (pushMid P F φ hφs.2 Z.hom.hom Z.hom.property.1 Z.hom.property.2 ≫ hh) ≫
          pushHom P F φ hφs.2 W.hom.hom W.hom.property.1 W.hom.property.2 = Z.hom.hom ≫ φ := by
        rw [Category.assoc, hhw, ← pushFac P F φ hφs.2 Z.hom.hom Z.hom.property.1 Z.hom.property.2]
      have f2 : pushMid P F φ hφs.2 W.hom.hom W.hom.property.1 W.hom.property.2 ≫
          pushHom P F φ hφs.2 W.hom.hom W.hom.property.1 W.hom.property.2 = W.hom.hom ≫ φ :=
        (pushFac P F φ hφs.2 W.hom.hom W.hom.property.1 W.hom.property.2).symm
      rw [divInv_congr P f1 _ (isBaseIsomorphism_comp P Z.hom.property.2.2 hφs.2),
        divInv_congr P f2 _ (isBaseIsomorphism_comp P W.hom.property.2.2 hφs.2),
        divInv_comp_left P Z.hom.hom φ Z.hom.property.2.2 hφs.2 Z.hom.property.1,
        divInv_comp_left P W.hom.hom φ W.hom.property.2.2 hφs.2 W.hom.property.1]
    have dkey : divInv P (pushMid P F φ hφs.2 Z.hom.hom Z.hom.property.1 Z.hom.property.2) hmZ.2
        = divInv P ψ' hψ'.2.2 := by
      haveI : IsIso (P.Base hh) := hhp.2.2
      refine Φ.map_injective (inv (P.Base hh)) ?_
      rw [← divInv_comp_right P _ hh hmZ.2 hhp.2.2 hhp.1 hhp.2.1,
        ← divInv_comp_right P _ hh hψ'.2.2 hhp.2.2 hhp.1 hhp.2.1, d4, ← d1]
    -- ★同型 `ee : Z.left.obj ≅ Cc'`
    obtain ⟨iso⟩ := coaPreOver_iso_of_divInv_eq P hover
      (Over.mk (⟨pushMid P F φ hφs.2 Z.hom.hom Z.hom.property.1 Z.hom.property.2, hmZ.1, hmZs⟩ :
        (⟨Z.left.obj⟩ : WideSubcategory (coaPreProp P)) ⟶ ⟨_⟩))
      (Over.mk (⟨ψ', hψ'.1, hψ'.2⟩ : (⟨Cc'⟩ : WideSubcategory (coaPreProp P)) ⟶ ⟨_⟩))
      dkey
    obtain ⟨ee, hee, heew0⟩ := coaPreOver_hom_of_iso P iso
    have heew : ee ≫ ψ' = pushMid P F φ hφs.2 Z.hom.hom Z.hom.property.1 Z.hom.property.2 := heew0
    -- ★求める射
    refine ⟨Over.homMk (⟨ee ≫ a, IsIsometric.comp P (isIsometric_of_isIso P ee) ha.1,
        IsPreStep.comp P (isPreStep_of_isIso P ee) ha.2⟩ : Z.left ⟶ W.left) ?_,
      imtrPreOver_hom_subsingleton P F _ _⟩
    refine InducedWideCategory.Hom.ext ?_
    show (ee ≫ a) ≫ W.hom.hom = Z.hom.hom
    refine (cancel_mono φ).mp ?_
    calc ((ee ≫ a) ≫ W.hom.hom) ≫ φ = (ee ≫ a) ≫ (W.hom.hom ≫ φ) := Category.assoc _ _ _
      _ = (ee ≫ a) ≫ (pushMid P F φ hφs.2 W.hom.hom W.hom.property.1 W.hom.property.2 ≫
            pushHom P F φ hφs.2 W.hom.hom W.hom.property.1 W.hom.property.2) := by
          rw [← pushFac P F φ hφs.2 W.hom.hom W.hom.property.1 W.hom.property.2]
      _ = ee ≫ ((a ≫ pushMid P F φ hφs.2 W.hom.hom W.hom.property.1 W.hom.property.2) ≫
            pushHom P F φ hφs.2 W.hom.hom W.hom.property.1 W.hom.property.2) := by
          simp only [Category.assoc]
      _ = ee ≫ ((ψ' ≫ hh) ≫
            pushHom P F φ hφs.2 W.hom.hom W.hom.property.1 W.hom.property.2) :=
          congrArg (fun t => ee ≫ (t ≫
            pushHom P F φ hφs.2 W.hom.hom W.hom.property.1 W.hom.property.2)) heq'.symm
      _ = (ee ≫ ψ') ≫ (hh ≫
            pushHom P F φ hφs.2 W.hom.hom W.hom.property.1 W.hom.property.2) := by
          simp only [Category.assoc]
      _ = pushMid P F φ hφs.2 Z.hom.hom Z.hom.property.1 Z.hom.property.2 ≫
            pushHom P F φ hφs.2 Z.hom.hom Z.hom.property.1 Z.hom.property.2 := by
          rw [heew, hhw]
          rfl
      _ = Z.hom.hom ≫ φ :=
          (pushFac P F φ hφs.2 Z.hom.hom Z.hom.property.1 Z.hom.property.2).symm
  -- ★★本質的全射性
  haveI hess : (pushFunctor P F φ hφs.2).EssSurj := by
    constructor
    intro Y
    obtain ⟨Cc, ψ, α, heq, hψ, hα⟩ := push_square P F hover φ hφc hφs
      Y.hom.hom Y.hom.property.1 Y.hom.property.2
    obtain ⟨δ, hδ1, -⟩ := prop_1_9_i_uniq P F
      (pushObj P F φ hφs.2 α hα.1 hα.2) Y.left.obj
      (pushMid P F φ hφs.2 α hα.1 hα.2) (pushHom P F φ hφs.2 α hα.1 hα.2)
      ψ Y.hom.hom
      ((pushFac P F φ hφs.2 α hα.1 hα.2).symm.trans heq.symm)
      (pushMid_spec P F φ hφs.2 α hα.1 hα.2).1 (pushMid_spec P F φ hφs.2 α hα.1 hα.2).2
      (pushHom_spec P F φ hφs.2 α hα.1 hα.2).1 (pushHom_spec P F φ hφs.2 α hα.1 hα.2).2
      hψ.1 hψ.2.2 Y.hom.property.1 Y.hom.property.2
    refine ⟨Over.mk (⟨α, hα⟩ : (⟨Cc⟩ : ImtrPre P) ⟶ (⟨A⟩ : ImtrPre P)), ⟨?_⟩⟩
    refine Over.isoMk ⟨⟨δ.hom, isIsometric_of_isIso P δ.hom, isPreStep_of_isIso P δ.hom⟩,
      ⟨δ.inv, isIsometric_of_isIso P δ.inv, isPreStep_of_isIso P δ.inv⟩,
      InducedWideCategory.Hom.ext δ.hom_inv_id,
      InducedWideCategory.Hom.ext δ.inv_hom_id⟩ ?_
    refine InducedWideCategory.Hom.ext ?_
    show δ.hom ≫ Y.hom.hom = pushHom P F φ hφs.2 α hα.1 hα.2
    rw [hδ1, ← Category.assoc, δ.hom_inv_id, Category.id_comp]
  exact ⟨hfaith, hfull, hess⟩

end PushEquiv2

/-! ## ★負の対照 —— (iv) の `co-angular` は落とせない

★(iv) は「co-angular **かつ** linear なら isotropic 性が両向きに移る」と言う。
`linear` だけでは `⇐` が壊れることを、`Proposition 1.4` で作った模型 `cP` で示す。

★`cP` は `Vee` 上の定数関手による pre-Frobenioid で、
**`vLeftTop : left ⟶ top` は linear だが co-angular でなく**、
**`top` は isotropic だが `left` は isotropic でない**。 -/

/-- `cP` の `vLeftTop` は **linear**。 -/
theorem cP_linear_vLeftTop : IsLinear cP vLeftTop := rfl

/-- `cP` の `top` は **isotropic** —— `top` から出る射は `𝟙` しかない。 -/
theorem cP_isotropic_top : IsIsotropic cP Vee.top := by
  intro Dd f _ _
  have hd : Dd = Vee.top := by
    rcases leOfHom f with h | h
    · exact h.symm
    · exact h
  subst hd
  haveI := Preorder.subsingleton_hom Vee.top Vee.top
  exact ⟨𝟙 _, Subsingleton.elim _ _, Subsingleton.elim _ _⟩

/-- ★★**(iv) の `co-angular` は落とせない** —— `vLeftTop` は linear で
終域 `top` は isotropic だが、始域 `left` は **isotropic でない**。

★`Proposition 1.4, (i)` の負の対照(`cP_not_coAngular`)がそのまま
`Proposition 1.9, (iv)` の負の対照になる。**同じ模型が2つの命題で効いた。** -/
theorem cP_prop_1_9_iv_coAngular_is_load_bearing :
    IsLinear cP vLeftTop ∧ IsIsotropic cP Vee.top ∧
      ¬ IsCoAngular cP vLeftTop ∧ ¬ IsIsotropic cP Vee.left :=
  ⟨cP_linear_vLeftTop, cP_isotropic_top, cP_not_coAngular, cP_not_isotropic_left⟩

/-! ### ★出典の紐付け(`.src`)

★`.src` は「その原典項目を**完全に**実装した」という主張である。
`Proposition 1.9` は (i)–(vii) がすべて揃ったので、ここで初めて付ける。
各条の内訳:

* **(i)**: `prop_1_9_i_factor`(分解)＋ `prop_1_9_i_uniq`(一意性)＋
  `prop_1_9_i_isometric_iff` / `prop_1_9_i_coAngular_iff` / `prop_1_9_i_pullBack_iff`
* **(ii)**: `pushFunctor`(関手)＋ `pushFunctor_isEquivalence`(圏同値)＋
  `OTimesImtrPre`(`𝒪^×(A)^imtr-pre`)
* **(iii)**: `pullFunctor`(関手)＋ `prop_1_9_iii_uniq`(同型を除く一意性)
* **(iv)**: `prop_1_9_iv`
* **(v)**: `istr_frobenioid`(`𝒞^istr` が Frobenioid)＋ `isotropificationAdj`(左随伴)＋
  `istr_compat_degFr` / `_Base` / `_Div`(`𝒞 → 𝔽_Φ` が経由すること)＋
  `isotropificationRestrictIso`(制限が恒等と同型)＋ 保存 11 種
  (`isotropification_frobType` / `_degFr` / `_preStep_iff` / `_pullBack` / `_baseIso_iff` /
  `_baseFSM` / `_baseIdentity` / `_divIdentity` / `_isometric_iff` / `istr_coAngular` /
  `_lbInvertible`)
* **(vi)**: `prop_1_9_vi`
* **(vii)**: `prop_1_9_vii`
-/


def pushFunctor_isEquivalence.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 31, item := "Proposition 1.9, (ii)",
    sectionId := "frdi-prop-1-9-ii" }






/-! ## ★★命題全体の `.src`(2026-08-16)

★上の `.src` はすべて**条つき**で、locator の記録である(器具の数には入らない)。
★下の `prop_1_9.src` が**命題全体を完全に実装したという主張**である。

★★**文脈を持たない検証役の監査(2026-08-16)で欠落が 2 件見つかり、どちらも同日に埋めた。**

原文 (FrdI p.32):
> forms a left adjoint to the inclusion functor Cistr →C, through which the functor

★★**(v) の 2 件**:
- 「through which the functor `𝒞 → 𝔽_Φ` factors」→ `isotropificationFactorIso`
  (`P.toElem ≅ isotropification P F ⋙ (istrPre P F).toElem`)。
  ★中身は「**isotropic hull は `𝔽_Φ` では同型になる**」(`toElem_map_hullMap_isIso`) ——
  isometric ＋ pre-step から `div = 0`・`deg = 1`・base 同型が揃うため。
- 保存 11 クラスのうち co-angular → `isotropification_coAngular`(仮定不要。原文より強い)。

★★**(ii) の 1 件**:

原文 (FrdI p.31):
> the subgroup of v ∈O×(A) for which vimtr-pre is the identity.

★**原文は `𝒪^×(A)^{imtr-pre}` が `𝒪^×(A)` の部分群であると述べている**が、
監査以前は `Set (End A)` と包含補題しか無かった。3 条件を埋めた:
- 単位元: `one_mem_otimesImtrPre`(`ε` 自身が isometric pre-step なので `ε = 𝟙 ≫ ε` が第2の分解)
- 積: `mul_mem_otimesImtrPre`(`ε ≫ (φ ≫ ψ)` の第2の分解 —— 左因子は
  co-angular base-isomorphism の**合成**、右因子は isometric pre-step)
- 逆元: `inv_mem_otimesImtrPre`(積と単位元から)
- 包装: `OTimesImtrPreSubmonoid`, `otimesImtrPreSubmonoid_le`

★★**共通の鍵**: `Over (⟨A⟩ : ImtrPre P)` の hom は subsingleton
(`imtrPreOver_hom_subsingleton`)。★**したがって自然性は自動**で、
対象ごとの同型さえ作れば関手の同型が得られる。

## ★条ごとの照合表

| 条 | 主張 | 宣言 |
|---|---|---|
| (i) | 分解の存在 | `prop_1_9_i_factor` |
| (i) | 一意性(`(α◦γ, γ⁻¹◦β)` を除く) | `prop_1_9_i_uniq` |
| (i) | φ isometric ⟺ β が Frobenius 型 | `prop_1_9_i_isometric_iff` |
| (i) | φ co-angular ⟺ α が同型 | `prop_1_9_i_coAngular_iff` |
| (i) | φ pull-back ⟺ φ が同型 | `prop_1_9_i_pullBack_iff` |
| (ii) | 関手 `φ_*` | `pushFunctor` |
| (ii) | **Moreover**: co-angular pre-step なら圏同値 | `pushFunctor_isEquivalence` |
| (ii) | `𝒪^×(A)^{imtr-pre}` と包含 | `OTimesImtrPre`, `otimesImtrPre_subset` |
| (ii) | ★**部分群であること**(監査で発見) | `OTimesImtrPreSubmonoid` ほか 3 本 |
| (iii) | 関手 `φ^*` | `pullFunctor` |
| (iii) | 四角形と、同型を除く一意性 | `pullFac`, `prop_1_9_iii_uniq`, `prop_1_9_iii_lift` |
| (iv) | iff | `prop_1_9_iv` |
| (v) | `𝒞^istr` が Frobenioid | `istr_frobenioid` |
| (v) | isotropification が左随伴 | `isotropificationAdj` |
| (v) | ★**`𝒞 → 𝔽_Φ` が経由する**(監査で発見) | `isotropificationFactorIso` |
| (v) | 制限が恒等と同型 | `isotropificationRestrictIso` |
| (v) | 保存 11 クラス(★`isotropification_coAngular` は監査で追加) | `isotropification_*` 11 本 |
| (v) | moreover: 包含関手との互換性 | `istr_compat_*`, `istr_frobType_iff`, `istr_isPullBack_*` |
| (vi) | iff ＋ moreover | `prop_1_9_vi` |
| (vii) | iff | `prop_1_9_vii` |

★**★印の 3 行が、監査以前には存在しなかった主張である。**
-/

def prop_1_9.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 30, item := "Proposition 1.9",
    sectionId := "frdi-prop-1-9-i" }

end ABC3.Found.FrdI
