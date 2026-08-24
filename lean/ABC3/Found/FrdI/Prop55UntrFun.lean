/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop55UntrCat

/-!
# [FrdI] Proposition 5.5, (ii) —— **関手**として組む

原文 (FrdI p.104):
> (ii) There is a natural equivalence of categories [compatible with the functors to the

★★`Prop55UntrCat.lean` は **射の集合の全単射**
`Hom_{(𝒞^pf)^un-tr}(A,B) ≃ Hom_{(𝒞^un-tr)^pf}(A,B)` まで作った。
★原文は「**圏同値**」と言うので、恒等と合成を保つことが要る。

## ★組み方

全単射は `HomColim.quotEquiv` ＋ `HasColimit.isoOfNatIso` ＋ `Functor.Final.colimitIso`
の合成で作られており、代表元の上での式が見えにくい。
★そこで本ファイルは**代表元の上で直に**写像を組み直す:

  `Hom^pf(A,B) ∋ mk Z φ ⟼ mk (idxToUnTr Z) (istrToUnTr φ) ∈ Hom^pf_{𝒞^un-tr}(A,B)`

★余錐の自然性は `istrToUnTr_frobTransport`(在庫)そのもの。
★恒等と合成は `toHomPf` と `compPf_mk`(在庫)に落ちる。
-/

namespace ABC3.Found.FrdI

open CategoryTheory CategoryTheory.Limits

universe v u w u2 v2

section UnTrFun

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ) (Fc : FrobenioidCore P)
  (F₁ : FrobenioidCore (istrPre P Fc)) (F₂ : FrobenioidCore (unTrPre P Fc))

/-- ★★`Hom^pf(A,B)` の上の余錐 —— 各段で `istrToUnTr` を当てる。

★自然性は `istrToUnTr_frobTransport`(遷移写像の対応)1 本。 -/
noncomputable def unTrPfCocone (A B : Istr P) :
    Cocone (homFunctorPf (istrPre P Fc) F₁ A B) where
  pt := HomPf (unTrPre P Fc) F₂ A B
  ι :=
    { app := fun Z => TypeCat.ofHom fun φ =>
        HomPf.mk ((idxToUnTr P Fc F₁ F₂ A B).obj Z) ((istrToUnTr P).map φ.down)
      naturality := by
        intro Z W u
        ext φ
        have hstep : (istrToUnTr P).map (idxTransport (istrPre P Fc) F₁ u φ.down)
            = idxTransport (unTrPre P Fc) F₂ ((idxToUnTr P Fc F₁ F₂ A B).map u)
              ((istrToUnTr P).map φ.down) :=
          istrToUnTr_frobTransport P Fc F₁ F₂
            u.right.hom.1 u.right.property.1 u.right.hom.2 u.right.property.2.1
            u.right.property.2.2 φ.down
            ((idxToUnTr P Fc F₁ F₂ A B).map u).right.property.1
            ((idxToUnTr P Fc F₁ F₂ A B).map u).right.property.2.1
            ((idxToUnTr P Fc F₁ F₂ A B).map u).right.property.2.2
        show HomPf.mk ((idxToUnTr P Fc F₁ F₂ A B).obj W)
            ((istrToUnTr P).map (idxTransport (istrPre P Fc) F₁ u φ.down))
          = HomPf.mk ((idxToUnTr P Fc F₁ F₂ A B).obj Z) ((istrToUnTr P).map φ.down)
        exact (congrArg (HomPf.mk ((idxToUnTr P Fc F₁ F₂ A B).obj W)) hstep).trans
          (HomPf.mk_map ((idxToUnTr P Fc F₁ F₂ A B).map u) _) }

/-- ★★**`Hom^pf_{𝒞^istr}(A,B) → Hom^pf_{𝒞^un-tr}(A,B)`**。 -/
noncomputable def homPfToUnTrPf (A B : Istr P) :
    HomPf (istrPre P Fc) F₁ A B → HomPf (unTrPre P Fc) F₂ A B :=
  colimit.desc (homFunctorPf (istrPre P Fc) F₁ A B) (unTrPfCocone P Fc F₁ F₂ A B)

variable {P Fc F₁ F₂} in
/-- ★代表元の上での式。 -/
@[simp] theorem homPfToUnTrPf_mk {A B : Istr P} (Z : IdxPf (istrPre P Fc) F₁ A B)
    (φ : Z.right.obj.1 ⟶ Z.right.obj.2) :
    homPfToUnTrPf P Fc F₁ F₂ A B (HomPf.mk Z φ)
      = HomPf.mk ((idxToUnTr P Fc F₁ F₂ A B).obj Z) ((istrToUnTr P).map φ) := by
  show colimit.desc (homFunctorPf (istrPre P Fc) F₁ A B) (unTrPfCocone P Fc F₁ F₂ A B)
      (colimit.ι (homFunctorPf (istrPre P Fc) F₁ A B) Z (ULift.up φ)) = _
  rw [← types_comp_apply (colimit.ι (homFunctorPf (istrPre P Fc) F₁ A B) Z)
    (colimit.desc (homFunctorPf (istrPre P Fc) F₁ A B) (unTrPfCocone P Fc F₁ F₂ A B)),
    colimit.ι_desc]
  rfl

/-! ## ★2. 3 つ組の添字圏でも同じこと -/

/-- ★**`𝒞^{tri-Fr}` のあいだの関手**(`biFrToUnTr` の 3 成分版)。 -/
def triFrToUnTr : TriFr (istrPre P Fc) F₁ ⥤ TriFr (unTrPre P Fc) F₂ where
  obj Z := ⟨((istrToUnTr P).obj Z.obj.1, (istrToUnTr P).obj Z.obj.2.1,
    (istrToUnTr P).obj Z.obj.2.2)⟩
  map {Z W} u :=
    ⟨((istrToUnTr P).map u.hom.1, (istrToUnTr P).map u.hom.2.1,
        (istrToUnTr P).map u.hom.2.2),
      (unTr_isFrobeniusType_iff' P Fc _).mpr u.property.1,
      (unTr_isFrobeniusType_iff' P Fc _).mpr u.property.2.1,
      (unTr_isFrobeniusType_iff' P Fc _).mpr u.property.2.2.1,
      u.property.2.2.2.1,
      u.property.2.2.2.2⟩
  map_id Z := by
    apply InducedWideCategory.Hom.ext
    exact Prod.ext ((istrToUnTr P).map_id _)
      (Prod.ext ((istrToUnTr P).map_id _) ((istrToUnTr P).map_id _))
  map_comp u u' := by
    apply InducedWideCategory.Hom.ext
    exact Prod.ext ((istrToUnTr P).map_comp _ _)
      (Prod.ext ((istrToUnTr P).map_comp _ _) ((istrToUnTr P).map_comp _ _))

/-- ★★**3 つ組の添字圏のあいだの関手**。 -/
def idx3ToUnTr (A B E : Istr P) :
    IdxPf3 (istrPre P Fc) F₁ A B E ⥤ IdxPf3 (unTrPre P Fc) F₂ A B E :=
  Under.post (triFrToUnTr P Fc F₁ F₂)

variable {P Fc F₁ F₂} in
/-- ★第 1・第 2 脚との両立(`rfl`)。 -/
theorem idx3ToUnTr_idx12 {A B E : Istr P} (V : IdxPf3 (istrPre P Fc) F₁ A B E) :
    (idxToUnTr P Fc F₁ F₂ A B).obj ((idx12 (istrPre P Fc) F₁ A B E).obj V)
      = (idx12 (unTrPre P Fc) F₂ A B E).obj ((idx3ToUnTr P Fc F₁ F₂ A B E).obj V) :=
  rfl

variable {P Fc F₁ F₂} in
/-- ★第 2・第 3 脚との両立(`rfl`)。 -/
theorem idx3ToUnTr_idx23 {A B E : Istr P} (V : IdxPf3 (istrPre P Fc) F₁ A B E) :
    (idxToUnTr P Fc F₁ F₂ B E).obj ((idx23 (istrPre P Fc) F₁ A B E).obj V)
      = (idx23 (unTrPre P Fc) F₂ A B E).obj ((idx3ToUnTr P Fc F₁ F₂ A B E).obj V) :=
  rfl

variable {P Fc F₁ F₂} in
/-- ★第 1・第 3 脚との両立(`rfl`)。 -/
theorem idx3ToUnTr_idx13 {A B E : Istr P} (V : IdxPf3 (istrPre P Fc) F₁ A B E) :
    (idxToUnTr P Fc F₁ F₂ A E).obj ((idx13 (istrPre P Fc) F₁ A B E).obj V)
      = (idx13 (unTrPre P Fc) F₂ A B E).obj ((idx3ToUnTr P Fc F₁ F₂ A B E).obj V) :=
  rfl

/-! ## ★3. 恒等と合成 -/

variable {P Fc F₁ F₂} in
/-- ★底の添字は動かない(`rfl`)。 -/
theorem idxToUnTr_idxOne (A B : Istr P) :
    (idxToUnTr P Fc F₁ F₂ A B).obj (idxOne (istrPre P Fc) F₁ A B)
      = idxOne (unTrPre P Fc) F₂ A B :=
  rfl

variable {P Fc F₁ F₂} in
/-- ★★`𝒞^istr` から来る射は `𝒞^un-tr` から来る射へ移る。 -/
theorem homPfToUnTrPf_toHomPf {A B : Istr P} (φ : A ⟶ B) :
    homPfToUnTrPf P Fc F₁ F₂ A B (toHomPf (F := F₁) φ)
      = toHomPf (F := F₂) ((istrToUnTr P).map φ) := by
  show homPfToUnTrPf P Fc F₁ F₂ A B (HomPf.mk (idxOne (istrPre P Fc) F₁ A B) φ) = _
  rw [homPfToUnTrPf_mk]
  rfl

variable {P Fc F₁ F₂} in
/-- ★★★**合成を保つ** —— 共通の 3 つ組の添字に揃えれば `≫` の保存 1 本。 -/
theorem homPfToUnTrPf_comp {A B E : Istr P}
    (f : HomPf (istrPre P Fc) F₁ A B) (g : HomPf (istrPre P Fc) F₁ B E) :
    homPfToUnTrPf P Fc F₁ F₂ A E (compPf (istrPre P Fc) F₁ f g)
      = compPf (unTrPre P Fc) F₂ (homPfToUnTrPf P Fc F₁ F₂ A B f)
          (homPfToUnTrPf P Fc F₁ F₂ B E g) := by
  obtain ⟨V, φ, ψ, rfl, rfl⟩ := exists_rep3 (P := istrPre P Fc) (F := F₁) f g
  rw [compPf_mk, homPfToUnTrPf_mk, homPfToUnTrPf_mk, homPfToUnTrPf_mk]
  refine Eq.trans (congrArg (HomPf.mk _) ((istrToUnTr P).map_comp φ ψ)) ?_
  exact (compPf_mk ((idx3ToUnTr P Fc F₁ F₂ A B E).obj V)
    ((istrToUnTr P).map φ) ((istrToUnTr P).map ψ)).symm

/-! ## ★4. unit-equivalence で割る -/

variable {P Fc F₁ F₂} in
set_option maxHeartbeats 800000 in
/-- ★★★**`𝔽_{Φ^pf}` で同じ像なら移り先も同じ** —— 商へ降りる根拠。

★共通の上界へ送れば `toHomUnTr_eq_iff`(`𝔽_Φ` の像の一致)そのものになる。 -/
theorem homPfToUnTrPf_congr {A B : Istr P} {f g : HomPf (istrPre P Fc) F₁ A B}
    (h : (pfPre (istrPre P Fc) F₁).toElem.map f
      = (pfPre (istrPre P Fc) F₁).toElem.map g) :
    homPfToUnTrPf P Fc F₁ F₂ A B f = homPfToUnTrPf P Fc F₁ F₂ A B g := by
  obtain ⟨Z, φ, rfl⟩ := HomPf.exists_rep f
  obtain ⟨W, ψ, rfl⟩ := HomPf.exists_rep g
  obtain ⟨V, u, u', hV⟩ :=
    (homPf_untr_eq_iff (P := istrPre P Fc) (F := F₁) (ULift.up φ) (ULift.up ψ)).mp
      ((toQuot_eq_iff_pfToElem (P := istrPre P Fc) (F := F₁) _ _).mpr h)
  rw [homPfToUnTrPf_mk, homPfToUnTrPf_mk]
  have hu : (istrToUnTr P).map (idxTransport (istrPre P Fc) F₁ u φ)
      = idxTransport (unTrPre P Fc) F₂ ((idxToUnTr P Fc F₁ F₂ A B).map u)
        ((istrToUnTr P).map φ) :=
    istrToUnTr_frobTransport P Fc F₁ F₂
      u.right.hom.1 u.right.property.1 u.right.hom.2 u.right.property.2.1
      u.right.property.2.2 φ
      ((idxToUnTr P Fc F₁ F₂ A B).map u).right.property.1
      ((idxToUnTr P Fc F₁ F₂ A B).map u).right.property.2.1
      ((idxToUnTr P Fc F₁ F₂ A B).map u).right.property.2.2
  have hu' : (istrToUnTr P).map (idxTransport (istrPre P Fc) F₁ u' ψ)
      = idxTransport (unTrPre P Fc) F₂ ((idxToUnTr P Fc F₁ F₂ A B).map u')
        ((istrToUnTr P).map ψ) :=
    istrToUnTr_frobTransport P Fc F₁ F₂
      u'.right.hom.1 u'.right.property.1 u'.right.hom.2 u'.right.property.2.1
      u'.right.property.2.2 ψ
      ((idxToUnTr P Fc F₁ F₂ A B).map u').right.property.1
      ((idxToUnTr P Fc F₁ F₂ A B).map u').right.property.2.1
      ((idxToUnTr P Fc F₁ F₂ A B).map u').right.property.2.2
  have hstep : (istrToUnTr P).map (idxTransport (istrPre P Fc) F₁ u φ)
      = (istrToUnTr P).map (idxTransport (istrPre P Fc) F₁ u' ψ) :=
    (toHomUnTr_eq_iff P _ _).mpr hV
  refine Eq.trans ?_ (Eq.trans (congrArg
      (HomPf.mk ((idxToUnTr P Fc F₁ F₂ A B).obj V)) (hu.symm.trans (hstep.trans hu'))) ?_)
  · exact (HomPf.mk_map ((idxToUnTr P Fc F₁ F₂ A B).map u) ((istrToUnTr P).map φ)).symm
  · exact HomPf.mk_map ((idxToUnTr P Fc F₁ F₂ A B).map u') ((istrToUnTr P).map ψ)

/-! ## ★5. 関手 -/

/-- ★★★★★★★**[FrdI] Proposition 5.5, (ii)** —— **`(𝒞^pf)^un-tr ⥤ (𝒞^un-tr)^pf`**。

原文 (FrdI p.104):
> (ii) There is a natural equivalence of categories [compatible with the functors to the

★対象は動かさず、射は `homPfToUnTrPf`(代表元に `istrToUnTr` を当てるだけ)。
★恒等・合成は `homPfToUnTrPf_toHomPf` / `homPfToUnTrPf_comp`。 -/
noncomputable def unTrPfFunctor :
    UnTr (pfPre (istrPre P Fc) F₁) ⥤ PfCat (unTrPre P Fc) F₂ where
  obj A := (show Istr P from (show Istr (pfPre (istrPre P Fc) F₁) from A).obj)
  map {A B} f :=
    Quotient.liftOn f
      (fun g => homPfToUnTrPf P Fc F₁ F₂
        (show Istr (pfPre (istrPre P Fc) F₁) from A).obj
        (show Istr (pfPre (istrPre P Fc) F₁) from B).obj g)
      (fun _ _ h => homPfToUnTrPf_congr h)
  map_id A := by
    show homPfToUnTrPf P Fc F₁ F₂ _ _ (toHomPf (F := F₁) (𝟙 _)) = toHomPf (F := F₂) (𝟙 _)
    rw [homPfToUnTrPf_toHomPf]
    exact congrArg (toHomPf (F := F₂)) ((istrToUnTr P).map_id _)
  map_comp {A B E} f g := by
    refine Quotient.inductionOn₂ f g (fun a b => ?_)
    exact homPfToUnTrPf_comp a b

/-! ## ★6. 既存の全単射との一致 -/

variable {P Fc F₁ F₂} in
/-- ★`prop_5_5_ii_untr` の代表元の上での式。 -/
theorem prop_5_5_ii_untr_mk {A B : Istr P} (Z : IdxPf (istrPre P Fc) F₁ A B)
    (φ : Z.right.obj.1 ⟶ Z.right.obj.2) :
    prop_5_5_ii_untr P Fc F₁ F₂ A B
        (Quotient.mk (HomColim.quotKer (unTrQuotPf (istrPre P Fc) F₁ A B)) (HomPf.mk Z φ))
      = HomPf.mk ((idxToUnTr P Fc F₁ F₂ A B).obj Z) ((istrToUnTr P).map φ) := by
  show (homPfUnTrColimIso P Fc F₁ F₂ A B).hom
      (HomColim.quotEquiv (unTrQuotPf (istrPre P Fc) F₁ A B)
        (Quotient.mk _ (HomPf.mk Z φ))) = _
  rw [HomColim.quotEquiv_apply]
  show (homPfUnTrColimIso P Fc F₁ F₂ A B).hom
      (HomColim.toQuot (unTrQuotPf (istrPre P Fc) F₁ A B)
        (HomColim.mk (homFunctorPf (istrPre P Fc) F₁ A B) Z (ULift.up φ))) = _
  rw [HomColim.toQuot_mk]
  show (Functor.Final.colimitIso (idxToUnTr P Fc F₁ F₂ A B)
        (homFunctorPf (unTrPre P Fc) F₂ A B)).hom
      ((HasColimit.isoOfNatIso (stageNatIso P Fc F₁ F₂ A B)).hom
        (colimit.ι (unTrQuotPf (istrPre P Fc) F₁ A B).functor Z
          (Quotient.mk _ (ULift.up φ)))) = _
  rw [← types_comp_apply (colimit.ι (unTrQuotPf (istrPre P Fc) F₁ A B).functor Z)
    (HasColimit.isoOfNatIso (stageNatIso P Fc F₁ F₂ A B)).hom,
    HasColimit.isoOfNatIso_ι_hom]
  show (Functor.Final.colimitIso (idxToUnTr P Fc F₁ F₂ A B)
        (homFunctorPf (unTrPre P Fc) F₂ A B)).hom
      (colimit.ι ((idxToUnTr P Fc F₁ F₂ A B) ⋙ homFunctorPf (unTrPre P Fc) F₂ A B) Z
        (ULift.up ((istrToUnTr P).map φ))) = _
  rw [← types_comp_apply
    (colimit.ι ((idxToUnTr P Fc F₁ F₂ A B) ⋙ homFunctorPf (unTrPre P Fc) F₂ A B) Z)
    (Functor.Final.colimitIso (idxToUnTr P Fc F₁ F₂ A B)
      (homFunctorPf (unTrPre P Fc) F₂ A B)).hom,
    Functor.Final.ι_colimitIso_hom]
  rfl

variable {P Fc F₁ F₂} in
/-- ★★★**関手の射の部分は、既存の全単射そのもの**。 -/
theorem unTrPfFunctor_map_eq {A B : UnTr (pfPre (istrPre P Fc) F₁)}
    (f : A ⟶ B) :
    (unTrPfFunctor P Fc F₁ F₂).map f
      = prop_5_5_ii_untr_cat P Fc F₁ F₂
          (show Istr (pfPre (istrPre P Fc) F₁) from A).obj
          (show Istr (pfPre (istrPre P Fc) F₁) from B).obj f := by
  refine Quotient.inductionOn f (fun g => ?_)
  obtain ⟨Z, φ, rfl⟩ := HomPf.exists_rep g
  show homPfToUnTrPf P Fc F₁ F₂ _ _ (HomPf.mk Z φ) = _
  rw [homPfToUnTrPf_mk]
  exact (prop_5_5_ii_untr_mk Z φ).symm

/-- ★★★★★★★**[FrdI] Proposition 5.5, (ii)** —— `unTrPfFunctor` は**充満忠実**。 -/
instance unTrPfFunctor_full : (unTrPfFunctor P Fc F₁ F₂).Full where
  map_surjective {A B} h := by
    obtain ⟨f, hf⟩ := (prop_5_5_ii_untr_cat P Fc F₁ F₂
      (show Istr (pfPre (istrPre P Fc) F₁) from A).obj
      (show Istr (pfPre (istrPre P Fc) F₁) from B).obj).surjective h
    exact ⟨f, (unTrPfFunctor_map_eq f).trans hf⟩

/-- ★★★★★★★**[FrdI] Proposition 5.5, (ii)** —— `unTrPfFunctor` は**忠実**。 -/
instance unTrPfFunctor_faithful : (unTrPfFunctor P Fc F₁ F₂).Faithful where
  map_injective {A B} {f g} h :=
    (prop_5_5_ii_untr_cat P Fc F₁ F₂
      (show Istr (pfPre (istrPre P Fc) F₁) from A).obj
      (show Istr (pfPre (istrPre P Fc) F₁) from B).obj).injective
      (((unTrPfFunctor_map_eq f).symm.trans h).trans (unTrPfFunctor_map_eq g))

/-! ## ★7. 圏同値 -/

/-- ★★★★★★★**[FrdI] Proposition 5.5, (ii)** ——
**`(𝒞^pf)^un-tr ≌ (𝒞^un-tr)^pf`**(圏同値)。

原文 (FrdI p.104):
> (ii) There is a natural equivalence of categories [compatible with the functors to the

★充満忠実は `prop_5_5_ii_untr_cat`(射の集合の全単射)から。
★本質的全射は「`𝒞^pf` の対象はすべて isotropic」(`hisoPf`)から ——
原文が「we may assume without loss of generality that `𝒞` is of isotropic type」と
置く条にあたる。 -/
noncomputable def unTrPfEquivalence
    (hisoPf : ∀ X : PfCat (istrPre P Fc) F₁, IsIsotropic (pfPre (istrPre P Fc) F₁) X) :
    UnTr (pfPre (istrPre P Fc) F₁) ≌ PfCat (unTrPre P Fc) F₂ := by
  haveI : (unTrPfFunctor P Fc F₁ F₂).EssSurj :=
    ⟨fun Y => ⟨(show UnTr (pfPre (istrPre P Fc) F₁) from
        (⟨(show PfCat (istrPre P Fc) F₁ from Y), hisoPf _⟩ :
          Istr (pfPre (istrPre P Fc) F₁))), ⟨Iso.refl _⟩⟩⟩
  haveI : (unTrPfFunctor P Fc F₁ F₂).IsEquivalence := { }
  exact (unTrPfFunctor P Fc F₁ F₂).asEquivalence

/-- ★★★★**`𝔽_{Φ^pf}` への射影と 1-可換(次数)** ——
原文の「compatible with the functors to the respective elementary Frobenioids」の一半。
★次数は代表元の次数そのものだから `rfl` に落ちる。 -/
theorem unTrPfFunctor_degFr (Fq : FrobenioidCore (pfPre (istrPre P Fc) F₁))
    {A B : UnTr (pfPre (istrPre P Fc) F₁)} (f : A ⟶ B) :
    (pfPre (unTrPre P Fc) F₂).degFr ((unTrPfFunctor P Fc F₁ F₂).map f)
      = (unTrPre (pfPre (istrPre P Fc) F₁) Fq).degFr f := by
  refine Quotient.inductionOn f (fun g => ?_)
  obtain ⟨Z, φ, rfl⟩ := HomPf.exists_rep g
  show pfDeg (homPfToUnTrPf P Fc F₁ F₂ _ _ (HomPf.mk Z φ)) = pfDeg (HomPf.mk Z φ)
  rw [homPfToUnTrPf_mk, pfDeg_mk, pfDeg_mk]
  rfl

/-- ★★★★**`𝔽_{Φ^pf}` への射影と 1-可換(底)**。 -/
theorem unTrPfFunctor_base (Fq : FrobenioidCore (pfPre (istrPre P Fc) F₁))
    {A B : UnTr (pfPre (istrPre P Fc) F₁)} (f : A ⟶ B) :
    (pfPre (unTrPre P Fc) F₂).Base ((unTrPfFunctor P Fc F₁ F₂).map f)
      = (unTrPre (pfPre (istrPre P Fc) F₁) Fq).Base f := by
  refine Quotient.inductionOn f (fun g => ?_)
  obtain ⟨Z, φ, rfl⟩ := HomPf.exists_rep g
  show pfBase (homPfToUnTrPf P Fc F₁ F₂ _ _ (HomPf.mk Z φ)) = pfBase (HomPf.mk Z φ)
  rw [homPfToUnTrPf_mk, pfBase_mk, pfBase_mk]
  rfl

end UnTrFun

/-! ### ★出典の紐付け -/

/-- ★★★★★★★locator —— `Proposition 5.5, (ii)` の圏同値(un-tr の側)。 -/
def unTrPfEquivalence.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 104,
    item := "Proposition 5.5, (ii) — (𝒞^pf)^un-tr ≌ (𝒞^un-tr)^pf(圏同値)",
    sectionId := "frdi-prop-5-5" }

/-- ★★★★locator —— その 1-可換性(底と次数)。 -/
def unTrPfFunctor_degFr.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 104,
    item := "Proposition 5.5, (ii) — 圏同値は 𝔽 への射影と 1-可換(底と次数)",
    sectionId := "frdi-prop-5-5" }

end ABC3.Found.FrdI
