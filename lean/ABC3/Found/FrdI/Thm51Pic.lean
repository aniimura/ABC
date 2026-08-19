/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Thm51Cls

/-!
# [FrdI] Theorem 5.1, (i) —— `Pic_Φ(A) ≅ Pic_𝒞(A)`

`Thm51Span.lean` の判定条件と `Thm51Cls.lean` の well-defined 性を束ねて、
原文の**全単射**を作る。

原文 (FrdI p.96):
> determines a bijection PicΦ(A) →∼ PicC(A). Moreover, if (C, ζ : CD →∼ AD) ∈

★`Pic_𝒞(A)` は「`𝒞 ×_𝒟 𝒟^iso_{A_𝒟}` の同型類」だが、
**ファイバー積の圏を作らずに**、対 `(C, ζ : C_𝒟 ≅ A_𝒟)` の型と
「`A_𝒟` の上での同型」という同値関係の商として直接作る(圏の骨格は要らない)。

★後半の「Frobenius 型の射で類は次数倍される」も本ファイルで取る。

原文 (FrdI p.96):
> corresponds to the element degFr(κ) · γ ∈ PicΦ(A).
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

/-! ## ★0. 小道具 -/

/-- `M^gp` の元は `toGp a − toGp b` の形。 -/
theorem exists_toGp_sub {M : Type w} [AddCommMonoid M] (x : Gp M) :
    ∃ a b : M, x = toGp M a - toGp M b := by
  induction x using AddLocalization.induction_on with
  | _ p => exact ⟨p.1, (p.2 : M), by rw [← mk_add_toGp M p.1 p.2]; abel⟩

section PhiTools

variable {D : Type u} [Category.{v} D] {Φ : MonoidOn.{v, u, w} D}

theorem phi_map_inv_comp {X Y : D} (g : X ⟶ Y) [IsIso g] (z : Φ.val Y) :
    Φ.map (inv g) (Φ.map g z) = z := by
  have h : Φ.map (inv g) (Φ.map g z) = Φ.map (inv g ≫ g) z := (Φ.map_comp _ _ z).symm
  rw [h, IsIso.inv_hom_id]
  exact Φ.map_id _ z

end PhiTools

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D}

/-! ## ★1. `Pic_Φ(A)` と `Pic_𝒞(A)` -/

/-- **[FrdI] Theorem 5.1** —— `Pic_Φ(A) = Φ^gp(A)/Φ^birat(A)`。

原文 (FrdI p.96):
> PicΦ(A) d=ef Φgp(A)/Φbirat(A)
-/
abbrev PicPhi (P : PreFrobenioid C Φ) (G : Frobenioid P) (A : C) : Type w :=
  Gp (Φ.val (P.toElem.obj A).base) ⧸ phiBiratAt P G A

/-- **[FrdI] Theorem 5.1** —— `𝒞 ×_𝒟 𝒟^iso_{A_𝒟}` の対象 `(C, ζ)`。 -/
structure PicObj (P : PreFrobenioid C Φ) (A : C) where
  obj : C
  iso : (P.toElem.obj obj).base ≅ (P.toElem.obj A).base

variable {P : PreFrobenioid C Φ}

/-- `Pic_𝒞(A)` の同値関係 —— `A_𝒟` の上での同型。 -/
def PicObj.Rel {A : C} (Z W : PicObj P A) : Prop :=
  ∃ ι : Z.obj ≅ W.obj, P.Base ι.hom ≫ W.iso.hom = Z.iso.hom

theorem PicObj.Rel.rfl' {A : C} (Z : PicObj P A) : PicObj.Rel Z Z :=
  ⟨Iso.refl _, by rw [Iso.refl_hom, P.Base_id, Category.id_comp]⟩

theorem PicObj.Rel.symm' {A : C} {Z W : PicObj P A} (h : PicObj.Rel Z W) : PicObj.Rel W Z := by
  obtain ⟨ι, hι⟩ := h
  refine ⟨ι.symm, ?_⟩
  show P.Base ι.inv ≫ Z.iso.hom = W.iso.hom
  rw [← hι, ← Category.assoc, ← P.Base_comp, ι.inv_hom_id, P.Base_id, Category.id_comp]

theorem PicObj.Rel.trans' {A : C} {Z W V : PicObj P A} (h : PicObj.Rel Z W)
    (h' : PicObj.Rel W V) : PicObj.Rel Z V := by
  obtain ⟨ι, hι⟩ := h
  obtain ⟨ι', hι'⟩ := h'
  refine ⟨ι.trans ι', ?_⟩
  show P.Base (ι.hom ≫ ι'.hom) ≫ V.iso.hom = Z.iso.hom
  rw [P.Base_comp, Category.assoc, hι', hι]

instance picSetoid (A : C) : Setoid (PicObj P A) where
  r := PicObj.Rel
  iseqv := ⟨PicObj.Rel.rfl', PicObj.Rel.symm', PicObj.Rel.trans'⟩

/-- **[FrdI] Theorem 5.1** —— `Pic_𝒞(A)`。 -/
abbrev PicC (P : PreFrobenioid C Φ) (A : C) : Type (max u2 v) :=
  Quotient (picSetoid (P := P) A)

/-! ## ★2. `PicObj` の類 -/

/-- `Z : PicObj P A` の類の候補 —— span `(φ, ψ)` から作る。 -/
def PicObj.HasCls {A : C} (Z : PicObj P A) (c : Gp (Φ.val (P.toElem.obj A).base)) : Prop :=
  ∃ (X : C) (φ : X ⟶ A) (ψ : X ⟶ Z.obj) (hsφ : IsPreStep P φ), IsPreStep P ψ ∧
    P.Base ψ ≫ Z.iso.hom = P.Base φ ∧ spanCls φ hsφ.2 ψ = c

theorem PicObj.exists_cls (G : Frobenioid P) {A : C} (Z : PicObj P A) :
    ∃ c, PicObj.HasCls Z c := by
  haveI : IsIso Z.iso.hom := Z.iso.isIso_hom
  obtain ⟨X, ψ, φ, hsψ, hsφ, hspan⟩ := G.core.preStepSpan Z.obj A Z.iso.hom inferInstance
  haveI hbψ : IsIso (P.Base ψ) := hsψ.2
  refine ⟨spanCls φ hsφ.2 ψ, X, φ, ψ, hsφ, hsψ, ?_, rfl⟩
  rw [hspan, ← Category.assoc, IsIso.hom_inv_id, Category.id_comp]

/-- ★span から `ζ⁻¹` が復元される。 -/
theorem picSpan_alpha {A : C} {Z : PicObj P A} {X : C} {φ : X ⟶ A} {ψ : X ⟶ Z.obj}
    (hsφ : IsPreStep P φ) (hb : P.Base ψ ≫ Z.iso.hom = P.Base φ) :
    @inv _ _ _ _ (P.Base φ) hsφ.2 ≫ P.Base ψ = Z.iso.inv := by
  haveI hbφ : IsIso (P.Base φ) := hsφ.2
  have hψ : P.Base ψ = P.Base φ ≫ Z.iso.inv := by
    rw [← hb, Category.assoc, Z.iso.hom_inv_id, Category.comp_id]
  rw [hψ, ← Category.assoc, IsIso.inv_hom_id, Category.id_comp]

/-- ★★同じ `Z` の 2 つの類は `Φ^birat` を法として等しい。 -/
theorem PicObj.cls_sub_mem (G : Frobenioid P) (hiso : ∀ Y : C, IsIsotropic P Y)
    {A : C} {Z : PicObj P A} {c c' : Gp (Φ.val (P.toElem.obj A).base)}
    (h : Z.HasCls c) (h' : Z.HasCls c') : c - c' ∈ phiBiratAt P G A := by
  obtain ⟨X, φ, ψ, hsφ, hsψ, hb, rfl⟩ := h
  obtain ⟨X', φ', ψ', hsφ', hsψ', hb', rfl⟩ := h'
  exact span_class_sub_mem G hiso φ ψ φ' ψ' hsφ hsψ hsφ' hsψ'
    ((picSpan_alpha hsφ hb).trans (picSpan_alpha hsφ' hb').symm)

/-- ★同値な `PicObj` は同じ類を持つ。 -/
theorem PicObj.hasCls_of_rel {A : C} {Z W : PicObj P A} (hZW : PicObj.Rel Z W)
    {c : Gp (Φ.val (P.toElem.obj A).base)} (h : Z.HasCls c) : W.HasCls c := by
  obtain ⟨ι, hι⟩ := hZW
  obtain ⟨X, φ, ψ, hsφ, hsψ, hb, hc⟩ := h
  haveI hbι : IsIso (P.Base ι.hom) := base_isIso_of_iso ι
  refine ⟨X, φ, ψ ≫ ι.hom, hsφ, IsPreStep.comp P hsψ ⟨degFr_of_isIso P ι.hom, hbι⟩, ?_, ?_⟩
  · rw [P.Base_comp, Category.assoc, hι, hb]
  · rw [← hc, spanCls_eq, spanCls_eq,
      Div_comp_preStep _ _ (degFr_of_isIso P ι.hom),
      show P.Div ι.hom = 0 from isIsometric_of_isIso P ι.hom, map_zero, zero_add]

/-- ★★★**類が `Φ^birat` を法として等しければ `PicObj` は同値**。 -/
theorem PicObj.rel_of_cls_sub_mem (G : Frobenioid P) (hiso : ∀ Y : C, IsIsotropic P Y)
    {A : C} {Z W : PicObj P A} {c c' : Gp (Φ.val (P.toElem.obj A).base)}
    (h : Z.HasCls c) (h' : W.HasCls c') (hsub : c - c' ∈ phiBiratAt P G A) :
    PicObj.Rel Z W := by
  obtain ⟨X, φ, ψ, hsφ, hsψ, hb, rfl⟩ := h
  obtain ⟨X', φ', ψ', hsφ', hsψ', hb', rfl⟩ := h'
  haveI hbφ : IsIso (P.Base φ) := hsφ.2
  haveI hbφ' : IsIso (P.Base φ') := hsφ'.2
  obtain ⟨Z₀, u, u', hu, hu', hspan⟩ :=
    G.core.preStepSpan X X' (P.Base φ ≫ inv (P.Base φ')) inferInstance
  haveI hbu : IsIso (P.Base u) := hu.2
  haveI hbu' : IsIso (P.Base u') := hu'.2
  have hkey : P.Base u ≫ P.Base φ ≫ inv (P.Base φ') = P.Base u' := by
    rw [hspan, ← Category.assoc, IsIso.hom_inv_id, Category.id_comp]
  have hbφz : P.Base (u ≫ φ) = P.Base (u' ≫ φ') := by
    rw [P.Base_comp, P.Base_comp, ← hkey, Category.assoc, Category.assoc,
      IsIso.inv_hom_id, Category.comp_id]
  have hs1 : IsPreStep P (u ≫ φ) := IsPreStep.comp P hu hsφ
  have hE := span_pullback_mem G hiso u hu φ hsφ u' hu' φ' hsφ' hbφz
  have htr := mem_phiBiratAt_transport G (u ≫ φ) (prop_1_4_i P _ (fun Y _ => hiso Y)) hs1
    (spanCls φ hsφ.2 ψ - spanCls φ' hsφ'.2 ψ')
  have hA := htr.mpr hsub
  rw [map_sub, gpMap_base_comp_spanCls u φ ψ hsφ hsψ,
    gpMap_base_comp_spanCls' u u' φ φ' ψ' hsφ' hsψ' hu hspan] at hA
  have hd₀ : toGp _ (P.Div (u' ≫ ψ')) - toGp _ (P.Div (u ≫ ψ)) ∈ phiBiratAt P G Z₀ := by
    have hrw : toGp (Φ.val (P.toElem.obj Z₀).base) (P.Div (u' ≫ ψ'))
          - toGp _ (P.Div (u ≫ ψ))
        = (toGp _ (P.Div (u' ≫ φ')) - toGp _ (P.Div (u ≫ φ)))
          - ((toGp _ (P.Div (u ≫ ψ)) - toGp _ (P.Div (u ≫ φ)))
            - (toGp _ (P.Div (u' ≫ ψ')) - toGp _ (P.Div (u' ≫ φ')))) := by abel
    rw [hrw]
    exact AddSubgroup.sub_mem _ hE hA
  obtain ⟨θ, hθ⟩ := span_iso_of_mem_phiBirat G hiso (u ≫ ψ) (u' ≫ ψ')
    (IsPreStep.comp P hu hsψ) (IsPreStep.comp P hu' hsψ') hd₀
  refine ⟨θ, ?_⟩
  haveI hbuψ : IsIso (P.Base (u ≫ ψ)) := (IsPreStep.comp P hu hsψ).2
  have e1 : P.Base (u ≫ ψ) ≫ Z.iso.hom = P.Base (u ≫ φ) := by
    rw [P.Base_comp, P.Base_comp, Category.assoc, hb]
  have e2 : P.Base (u' ≫ ψ') ≫ W.iso.hom = P.Base (u' ≫ φ') := by
    rw [P.Base_comp, P.Base_comp, Category.assoc, hb']
  rw [hθ, Category.assoc, e2, ← hbφz, ← e1, ← Category.assoc, IsIso.inv_hom_id,
    Category.id_comp]

/-- ★★★**どの類も実現される**(全射性)。 -/
theorem PicObj.exists_hasCls (G : Frobenioid P) {A : C}
    (c : Gp (Φ.val (P.toElem.obj A).base)) : ∃ Z : PicObj P A, Z.HasCls c := by
  obtain ⟨y, x, rfl⟩ := exists_toGp_sub c
  obtain ⟨X, φ, hpφ, hxφ⟩ := exists_coaPreStep_into G A x
  haveI hbφ : IsIso (P.Base φ) := hpφ.2.2
  obtain ⟨Y, ψ, hcψ, hsψ, hdψ⟩ := exists_coaPreStep_div G X (Φ.map (P.Base φ) y)
  haveI hbψ : IsIso (P.Base ψ) := hsψ.2
  refine ⟨⟨Y, (asIso (P.Base ψ)).symm ≪≫ asIso (P.Base φ)⟩, X, φ, ψ, hpφ.2, hsψ, ?_, ?_⟩
  · show P.Base ψ ≫ (inv (P.Base ψ) ≫ P.Base φ) = P.Base φ
    rw [← Category.assoc, IsIso.hom_inv_id, Category.id_comp]
  · rw [spanCls_eq, map_sub, gpMap_toGp, gpMap_toGp, hdψ, phi_map_inv_comp, hxφ]

/-- ★★★★**[FrdI] Theorem 5.1, (i) の後半** —— Frobenius 型の射で類は次数倍される。

原文 (FrdI p.96):
> corresponds to the element degFr(κ) · γ ∈ PicΦ(A).
-/
theorem PicObj.hasCls_frobType (G : Frobenioid P) {A : C} (hA : IsFrobeniusTrivial P A)
    {Z : PicObj P A} {c : Gp (Φ.val (P.toElem.obj A).base)} (h : Z.HasCls c)
    {Y' : C} (κ : Z.obj ⟶ Y') (hκ : IsFrobeniusType P κ) :
    (⟨Y', (@asIso _ _ _ _ (P.Base κ) hκ.2).symm ≪≫ Z.iso⟩ : PicObj P A).HasCls
      (((P.degFr κ : ℕ+) : ℕ) • c) := by
  haveI hbκ : IsIso (P.Base κ) := hκ.2
  obtain ⟨X, φ, ψ, hsφ, hsψ, hb, rfl⟩ := h
  obtain ⟨ζ, hdegζ, hpropζ⟩ := hA
  set κA : A ⟶ A := ((ζ (P.degFr κ) : End A) : A ⟶ A) with hκAdef
  have hκAb : P.Base κA = P.Base (𝟙 A) := (hpropζ _).1
  have hκAd : P.Div κA = 0 := (hpropζ _).2.1.2
  have hκAdeg : P.degFr κA = P.degFr κ := hdegζ _
  obtain ⟨Y, μ, hμ, hμdeg⟩ := G.core.frobDegSurj X (P.degFr κ)
  have hμd : P.Div μ = 0 := hμ.1.2
  haveI hμb : IsIso (P.Base μ) := hμ.2
  haveI hbκA : IsIso (P.Base κA) := by rw [hκAb, P.Base_id]; infer_instance
  obtain ⟨φ₂, hsφ₂, heqφ⟩ := frobShift G φ hsφ κA hbκA μ hμ (by rw [hμdeg, hκAdeg])
  obtain ⟨ψ₂, hsψ₂, heqψ⟩ := frobShift G ψ hsψ κ hbκ μ hμ (by rw [hμdeg])
  obtain ⟨hbφ₂, hdφ₂⟩ := frobShift_data φ κA hκAd μ hμd φ₂ heqφ
  obtain ⟨hbψ₂, hdψ₂⟩ := frobShift_data ψ κ hκ.1.2 μ hμd ψ₂ heqψ
  rw [hκAdeg] at hdφ₂
  rw [hκAb, P.Base_id, Category.comp_id] at hbφ₂
  refine ⟨Y, φ₂, ψ₂, hsφ₂, hsψ₂, ?_, ?_⟩
  · haveI : Epi (P.Base μ) := IsIso.epi_of_iso _
    refine (cancel_epi (P.Base μ)).mp ?_
    show P.Base μ ≫ P.Base ψ₂ ≫ (inv (P.Base κ) ≫ Z.iso.hom) = P.Base μ ≫ P.Base φ₂
    rw [← Category.assoc, hbψ₂, hbφ₂, Category.assoc, ← Category.assoc (P.Base κ),
      IsIso.hom_inv_id, Category.id_comp, hb]
  · exact spanCls_frobShift_cls φ ψ hsφ μ hμb φ₂ ψ₂ hsφ₂ _ hbφ₂ hdφ₂ hdψ₂

/-! ## ★3. ★★★★全単射 `Pic_𝒞(A) ≃ Pic_Φ(A)` -/

noncomputable def PicObj.cls (G : Frobenioid P) {A : C} (Z : PicObj P A) :
    Gp (Φ.val (P.toElem.obj A).base) := (PicObj.exists_cls G Z).choose

theorem PicObj.cls_spec (G : Frobenioid P) {A : C} (Z : PicObj P A) : Z.HasCls (Z.cls G) :=
  (PicObj.exists_cls G Z).choose_spec

theorem picPhi_mk_eq {A : C} (G : Frobenioid P) {a b : Gp (Φ.val (P.toElem.obj A).base)} :
    (QuotientAddGroup.mk a : PicPhi P G A) = QuotientAddGroup.mk b
      ↔ a - b ∈ phiBiratAt P G A :=
  @QuotientAddGroup.eq_iff_sub_mem _ _ (phiBiratAt P G A)
    (AddSubgroup.normal_of_isAddCommutative _) a b

/-- ★★★★**[FrdI] Theorem 5.1, (i)** —— `Pic_𝒞(A) → Pic_Φ(A)`。 -/
noncomputable def picCls (G : Frobenioid P) (hiso : ∀ Y : C, IsIsotropic P Y) (A : C) :
    PicC P A → PicPhi P G A :=
  Quotient.lift (fun Z => (QuotientAddGroup.mk (Z.cls G) : PicPhi P G A)) (by
    intro Z W hZW
    refine (picPhi_mk_eq G).mpr ?_
    exact PicObj.cls_sub_mem G hiso (PicObj.hasCls_of_rel hZW (Z.cls_spec G)) (W.cls_spec G))

theorem picCls_mk (G : Frobenioid P) (hiso : ∀ Y : C, IsIsotropic P Y) {A : C}
    (Z : PicObj P A) :
    picCls G hiso A (Quotient.mk _ Z) = QuotientAddGroup.mk (Z.cls G) := rfl

theorem picCls_eq_of_hasCls (G : Frobenioid P) (hiso : ∀ Y : C, IsIsotropic P Y) {A : C}
    {Z : PicObj P A} {c : Gp (Φ.val (P.toElem.obj A).base)} (h : Z.HasCls c) :
    picCls G hiso A (Quotient.mk _ Z) = QuotientAddGroup.mk c :=
  (picPhi_mk_eq G).mpr (PicObj.cls_sub_mem G hiso (Z.cls_spec G) h)

/-- ★★★★**[FrdI] Theorem 5.1, (i)** —— 全単射性。

原文 (FrdI p.96):
> determines a bijection PicΦ(A) →∼ PicC(A). Moreover, if (C, ζ : CD →∼ AD) ∈
-/
theorem picCls_bijective (G : Frobenioid P) (hiso : ∀ Y : C, IsIsotropic P Y) (A : C) :
    Function.Bijective (picCls G hiso A) := by
  constructor
  · refine Quotient.ind fun Z => Quotient.ind fun W => ?_
    intro h
    rw [picCls_mk, picCls_mk] at h
    exact Quotient.sound (PicObj.rel_of_cls_sub_mem G hiso (Z.cls_spec G) (W.cls_spec G)
      ((picPhi_mk_eq G).mp h))
  · intro q
    induction q using QuotientAddGroup.induction_on with
    | _ c =>
      obtain ⟨Z, hZ⟩ := PicObj.exists_hasCls G c
      exact ⟨Quotient.mk _ Z, picCls_eq_of_hasCls G hiso hZ⟩

/-- ★★★★**[FrdI] Theorem 5.1, (i)** —— `Pic_𝒞(A) ≃ Pic_Φ(A)`。 -/
noncomputable def picEquiv (G : Frobenioid P) (hiso : ∀ Y : C, IsIsotropic P Y) (A : C) :
    PicC P A ≃ PicPhi P G A :=
  Equiv.ofBijective _ (picCls_bijective G hiso A)

/-! ### ★出典の紐付け(`.src`) -/

/-- ★locator —— `Theorem 5.1, (i)`。 -/
def picEquiv.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 96, item := "Theorem 5.1, (i)",
    sectionId := "frdi-thm-5-1" }

end ABC3.Found.FrdI
