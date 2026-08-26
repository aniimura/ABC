/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop55PfArb
import ABC3.Found.FrdI.Prop53PfRlf

/-!
# [FrdI] `Proposition 5.5, (iii)(iv)` / `Proposition 5.3` —— `hftr` を落とした版

原文 (FrdI p.104):
> the case of arbitrary A then follows by considering "pairs of pre-steps" as in Theorem

★★`Prop55PfArb.lean` で **`hftr : ∀ X : C, IsFrobeniusTrivial P X`** が落ちたので、
そこから作られていた

* `pfRoot_ratFnData` —— `𝒞^pf` の有理関数の単系の interface
* `pfRoot_modelFrobenioid` —— `𝒞^pf ≌ model Frobenioid`
* `pfRoot_ratFnData_bmon_val` —— その単系は `ℚ·Φ^birat`
* `pfRootToRlfFunctor` —— `Proposition 5.3` の図式の右の縦の矢印

を **`hftr` なし**で立て直す。★仮定は原文 `Proposition 5.5` の
「Frobenius-isotropic ＋ Frobenius-normalized 型」＋ unit-trivial 型だけになる。

★★★これは**仮定の弱化であって、新しい数学ではない**——
中身はすべて `Prop55PfArb.lean` の 1 段(pre-step の対)である。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

section ArbStd

variable {D : Type u} [Category.{v} D] [IsConnected D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ} {F : FrobenioidCore P} {G : Frobenioid P}

/-- ★★★★★**`𝒞^pf` の有理関数の単系の interface**(★`hftr` なしの版)。 -/
noncomputable def pfRoot_ratFnData_of_arb (hfi : IsOfFrobeniusIsotropicType P)
    (hiso : ∀ X : C, IsIsotropic P X)
    (hfn : ∀ X : C, IsFrobeniusNormalized P X)
    (hut : ∀ X : C, IsUnitTrivial P X)
    (Gpf : Frobenioid (pfRootPre P F)) (hfsmD : IsOfFSMType D) :
    RatFnData (pfRootPre P F) Gpf :=
  biratRatFnData Gpf (pfRoot_isOfIsotropicType (F := F) hfi)
    (fun Z => (pfRoot_isOfModelType_of_arb hfi hiso hfn hut Gpf).2 Z)
    (pfRoot_isOfUnitTrivialType_of_arb hfi hiso hfn hut Gpf)
    (fun A => (Pf.isDivisorial' (P.divisorial A)).1.1) hfsmD

set_option maxHeartbeats 800000 in
/-- ★★★★★★**[FrdI] Proposition 5.5, (iv) の圏の側**(★`hftr` なしの版)——
**`𝒞^pf` は model Frobenioid と圏同値**。 -/
noncomputable def pfRoot_modelFrobenioid_of_arb (hfi : IsOfFrobeniusIsotropicType P)
    (hiso : ∀ X : C, IsIsotropic P X)
    (hfn : ∀ X : C, IsFrobeniusNormalized P X)
    (hut : ∀ X : C, IsUnitTrivial P X)
    (Gpf : Frobenioid (pfRootPre P F)) (hfsmD : IsOfFSMType D) :
    PfRootObj P F ≌ ModelData.Obj (pfRoot_ratFnData_of_arb hfi hiso hfn hut Gpf hfsmD).model :=
  modelType_equiv _ (pfRoot_isOfIsotropicType (F := F) hfi)
    (pfRoot_isOfModelType_of_arb hfi hiso hfn hut Gpf)

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★**[FrdI] Proposition 5.5, (iv) の単系の同定**(★`hftr` なしの版)——
`𝒞^pf` の有理関数の単系は各 `d ∈ Ob(𝒟)` で **`ℚ·Φ^birat(d)`**。 -/
theorem pfRoot_ratFnData_of_arb_bmon_val (hfi : IsOfFrobeniusIsotropicType P)
    (hiso : ∀ X : C, IsIsotropic P X)
    (hfnC : ∀ X : C, IsFrobeniusNormalized P X)
    (hut : ∀ X : C, IsUnitTrivial P X)
    (Gpf : Frobenioid (pfRootPre P F)) (hfsmD : IsOfFSMType D)
    (F' : FrobenioidCore (biratPre P G))
    (hisoB : ∀ X : BiratCat P G, IsIsotropic (biratPre P G) X)
    (hfnBirat : ∀ X : BiratCat P G, IsFrobeniusNormalized (biratPre P G) X)
    (hAB : ∀ A : C, IsFrobeniusTrivial (biratPre P G) (biratUp P G A))
    (hfnB' : ∀ A : C, IsFrobeniusNormalized (biratPre P G)
      (rtObj (biratPre P G) F' (biratUp P G A) 1))
    (ζ : ∀ A : C, ℕ+ →* End (biratUp P G A))
    (hdeg : ∀ (A : C) (n : ℕ+), (biratPre P G).degFr
      ((ζ A n : End (biratUp P G A)) : biratUp P G A ⟶ biratUp P G A) = n)
    (hprop : ∀ (A : C) (n : ℕ+), IsBaseIdentity (biratPre P G) (ζ A n)
      ∧ IsFrobeniusType (biratPre P G) ((ζ A n : End (biratUp P G A)) : _ ⟶ _))
    (hfnPfBirat : ∀ X : BiratCat (pfRootPre P F) Gpf,
      IsFrobeniusNormalized (biratPre (pfRootPre P F) Gpf) X)
    (d : D) :
    (pfRoot_ratFnData_of_arb (F := F) hfi hiso hfnC hut Gpf hfsmD).bmon.val d
      = ↥(qPhiBiratOn P G d) := by
  show ↥(phiBiratOn Gpf d) = ↥(qPhiBiratOn P G d)
  rw [phiBiratOn_pf_eq_qPhiBiratOn_all hfi hiso Gpf F' hisoB hfnBirat hAB hfnB'
    ζ hdeg hprop hfnPfBirat d]
  rfl

set_option maxHeartbeats 800000 in
/-- ★★★★★★★**`Proposition 5.3` の図式の右の縦の矢印**(★`hftr` なしの版)
`𝒞^pf ⟶ 𝒞^rlf`。 -/
noncomputable def pfRootToRlfFunctor_of_arb (hfi : IsOfFrobeniusIsotropicType P)
    (hiso : ∀ X : C, IsIsotropic P X)
    (hfnC : ∀ X : C, IsFrobeniusNormalized P X)
    (hut : ∀ X : C, IsUnitTrivial P X)
    (Gpf : Frobenioid (pfRootPre P F)) (hfsmD : IsOfFSMType D)
    (hfnB : ∀ X : BiratCat P G, IsFrobeniusNormalized (biratPre P G) X)
    (hcharInj' : ∀ {A B : D} (α : B ⟶ A),
      IsCharacteristicallyInjective (scMap (S := NNReal) (Φ.map α)))
    (hint' : ∀ A : D, IsIntegralMonoid (RlfT (Φ.val A)))
    (heq : ∀ d : D, phiBiratOn Gpf d = qPhiBiratOn P G d) :
    ModelData.Obj (pfRoot_ratFnData_of_arb (F := F) hfi hiso hfnC hut Gpf hfsmD).model
      ⥤ Crlf G hiso hfnB hcharInj' hint' hfsmD :=
  ({ phiHom := fun d => pfToRlfHom (Φ.val d)
     phiNat := fun f x => pfToRlfHom_map (Φ.map f) x
     bmonHom := fun d =>
       AddMonoidHom.codRestrict
         ((gpMap _ (pfToRlfHom (Φ.val d))).comp (phiBiratOn Gpf d).subtype) _
         (fun x => qPhiBiratOn_map_rlf G (by rw [← heq d]; exact x.2))
     bmonNat := by
       intro A B f x
       exact Subtype.ext (gpMap_pfToRlfHom_map (Φ.map f) _)
     divCompat := fun _ _ => rfl } :
    ModelDataHom (pfRoot_ratFnData_of_arb (F := F) hfi hiso hfnC hut Gpf hfsmD).model
      (scModel NNReal G hiso hfnB hcharInj' hint' hfsmD)).functor

end ArbStd

/-! ### ★出典の紐付け -/

/-- ★★★★★★★locator —— `Proposition 5.5, (iv)`(`hftr` なしの版)。 -/
def pfRoot_ratFnData_of_arb_bmon_val.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 105,
    item := "Proposition 5.5, (iv) — 𝒞^pf の有理関数の単系は ℚ·Φ^birat(Frobenius-trivial 型の仮定なし)",
    sectionId := "frdi-prop-5-5" }

/-- ★★★★★★★locator —— `Proposition 5.3` の右の縦の矢印(`hftr` なしの版)。 -/
def pfRootToRlfFunctor_of_arb.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 103,
    item := "Proposition 5.3 — 図式の右の縦の矢印 𝒞^pf ⥤ 𝒞^rlf(Frobenius-trivial 型の仮定なし)",
    sectionId := "frdi-prop-5-3" }

end ABC3.Found.FrdI
