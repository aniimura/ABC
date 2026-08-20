/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Divisor.NormPhiFunctor
import ABC3.Found.Divisor.SchemeDivFinite

/-!
# `Example 6.1` の model Frobenioid

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.109–110。

原文 (FrdI p.110):
> for the group of rational functions on Vi[Li] whose zeroes and poles are supported

## ★★`Example 6.3` と同じ骨組み

`Ex63Model.lean`(算術)と同じ 4 段である:

1. `Φ^gp ≃ Γ`(`effSubGpEquiv`)が引き戻しと両立すること
2. `B(L)` が `𝒟` 上の monoid であること
3. `div : B(L) → Φ(L)^gp` とその関手性
4. `Theorem 5.2` へ入れて model Frobenioid を得ること

## ★本ファイルで閉じること

| 定義 | 中身 |
|---|---|
| `phiGpHomC` | `Φ^gp → Γ`(Cartier 版) |
| `phiGpHomC_gpMapOn` | ★`Φ^gp ≃ Γ` は引き戻しと両立する |
-/

namespace ABC3.Found.Divisor

open CategoryTheory AlgebraicGeometry ABC3.Found.FrdI

universe v u w

/-! ## ★1. `Φ^gp ≃ Γ`(Cartier 版) -/

variable {D : Type u} [Category.{v} D] (Δ : CartierDatum.{v, u, w} D) (hD : IsOfFSMType D)

/-- ★`Φ^gp → Γ`(束ねた型の上で)。 -/
noncomputable def phiGpHomC {A : D} : Gp ((Δ.phi hD).val A) →+ Δ.grp A :=
  effSubGpHom (Δ.grp A)

@[simp] theorem phiGpHomC_toGp {A : D} (a : effSub (Δ.grp A)) :
    ((phiGpHomC Δ hD (toGp _ a) : Δ.grp A) : Δ.primes A →₀ ℤ)
      = ((a : Δ.primes A →₀ ℤ)) := by
  show ((effSubGpHom (Δ.grp A) (toGp _ a) : Δ.grp A) : Δ.primes A →₀ ℤ) = _
  rw [effSubGpHom, gpLift_toGp]
  rfl

theorem phiGpHomC_toGp' {A : D} (a : effSub (Δ.grp A)) :
    phiGpHomC Δ hD (toGp _ a) = effSubIncl (Δ.grp A) a :=
  Subtype.ext (phiGpHomC_toGp Δ hD a)

/-- ★★★**`Φ^gp ≃ Γ` は引き戻しと両立する**。 -/
theorem phiGpHomC_gpMapOn {A B : D} (α : B ⟶ A) (z : Gp ((Δ.phi hD).val A)) :
    phiGpHomC Δ hD ((Δ.phi hD).gpMapOn α z) = Δ.pull α (phiGpHomC Δ hD z) := by
  obtain ⟨a, b, rfl⟩ := gp_eq_sub z
  have key : ∀ c : (Δ.phi hD).val A,
      phiGpHomC Δ hD ((Δ.phi hD).gpMapOn α (toGp _ c))
        = Δ.pull α (phiGpHomC Δ hD (toGp _ c)) := by
    intro c
    have h : ((Δ.phi hD).gpMapOn α) (toGp ((Δ.phi hD).val A) c)
        = toGp ((Δ.phi hD).val B) ((Δ.phi hD).map α c) :=
      MonoidOn.gpMapOn_toGpHom (Δ.phi hD) α c
    rw [h, phiGpHomC_toGp', phiGpHomC_toGp']
    exact Subtype.ext (Δ.mapHom_coe α c)
  simp only [map_sub, key]

theorem phiGpHomC_injective {A : D} :
    Function.Injective (phiGpHomC Δ hD (A := A)) :=
  effSubGpHom_injective (Δ.grp A)


/-! ## ★2. `B(L)` は `𝒟` 上の monoid -/

section Geom

variable (V : Scheme.{u}) [IsIntegral V] {Kbar : Type u} [Field Kbar]
  [Algebra V.functionField Kbar] [IsGalois V.functionField Kbar]
  (DK : Set (PrimeDivisorPt V))
  [∀ L : FinSub V.functionField Kbar, IsLocallyNoetherian (normObj V L)]
  [∀ L : FinSub V.functionField Kbar, CompactSpace (normObj V L)]

/-- ★`B(L) → B(M)`(加法的に書いたもの)。 -/
noncomputable def bHom {L M : FinSub V.functionField Kbar} (f : L ⟶ M) :
    Additive (BSubgroup V DK L (normObj_isNormalScheme V L))
      →+ Additive (BSubgroup V DK M (normObj_isNormalScheme V M)) :=
  AddMonoidHom.mk' (fun x => Additive.ofMul (BSubgroupMap V DK f (Additive.toMul x)))
    (fun x y => by
      show Additive.ofMul (BSubgroupMap V DK f (Additive.toMul x * Additive.toMul y)) = _
      rw [map_mul]
      rfl)

omit [∀ L : FinSub V.functionField Kbar, CompactSpace (normObj V L)] in
omit [IsGalois V.functionField Kbar] in
@[simp] theorem bHom_val {L M : FinSub V.functionField Kbar} (f : L ⟶ M)
    (x : Additive (BSubgroup V DK L (normObj_isNormalScheme V L))) :
    ((Additive.toMul (bHom V DK f x) : BSubgroup V DK M _)
        : ((normObj V M).functionField)ˣ)
      = normFFUnits V f ((Additive.toMul x : BSubgroup V DK L _)
        : ((normObj V L).functionField)ˣ) := rfl

omit [∀ L : FinSub V.functionField Kbar, CompactSpace (normObj V L)] in
omit [IsGalois V.functionField Kbar] in
theorem bHom_injective {L M : FinSub V.functionField Kbar} (f : L ⟶ M) :
    Function.Injective (bHom V DK f) := by
  intro x y h
  have h' := congrArg (fun t => ((Additive.toMul t : BSubgroup V DK M _)
    : ((normObj V M).functionField)ˣ)) h
  rw [bHom_val, bHom_val] at h'
  refine Additive.toMul.injective (Subtype.ext (Units.ext ?_))
  exact (normFF V f).hom.injective (congrArg Units.val h')

omit [∀ L : FinSub V.functionField Kbar, CompactSpace (normObj V L)] in
omit [IsGalois V.functionField Kbar] in
theorem bHom_id (L : FinSub V.functionField Kbar) :
    bHom V DK (𝟙 L) = AddMonoidHom.id _ := by
  refine AddMonoidHom.ext fun x => ?_
  refine Additive.toMul.injective (Subtype.ext (Units.ext ?_))
  show (normFF V (𝟙 L)) (((Additive.toMul x : BSubgroup V DK L _)
    : ((normObj V L).functionField)ˣ) : (normObj V L).functionField) = _
  rw [normFF_id]
  rfl

omit [∀ L : FinSub V.functionField Kbar, CompactSpace (normObj V L)] in
omit [IsGalois V.functionField Kbar] in
theorem bHom_comp {L M N : FinSub V.functionField Kbar} (f : L ⟶ M) (g : M ⟶ N) :
    bHom V DK (f ≫ g) = (bHom V DK g).comp (bHom V DK f) := by
  refine AddMonoidHom.ext fun x => ?_
  refine Additive.toMul.injective (Subtype.ext (Units.ext ?_))
  show (normFF V (f ≫ g)) (((Additive.toMul x : BSubgroup V DK L _)
    : ((normObj V L).functionField)ˣ) : (normObj V L).functionField) = _
  rw [normFF_comp]
  rfl

/-- ★★**`L ↦ Additive (B(L))`** —— `𝒟ᵒᵖ ⥤ 𝔐𝔬𝔫`。 -/
noncomputable def bMonoidFunctor :
    ((FinSub V.functionField Kbar)ᵒᵖ)ᵒᵖ ⥤ AddCommMonCat.{u} where
  obj X := AddCommMonCat.of
    (Additive (BSubgroup V DK X.unop.unop (normObj_isNormalScheme V X.unop.unop)))
  map {_ _} f := AddCommMonCat.ofHom (bHom V DK f.unop.unop)
  map_id X := by
    show AddCommMonCat.ofHom (bHom V DK (𝟙 X.unop.unop)) = _
    rw [bHom_id]
    rfl
  map_comp {_ _ _} f g := by
    show AddCommMonCat.ofHom (bHom V DK (f.unop.unop ≫ g.unop.unop)) = _
    rw [bHom_comp]
    rfl

/-- ★★★★**`B(L)` は `𝒟` 上の monoid**。 -/
noncomputable def bmonGeom : MonoidOn.{u, u, u} ((FinSub V.functionField Kbar)ᵒᵖ) where
  functor := bMonoidFunctor V DK
  charInj {_ _} α := ⟨bHom_injective V DK α.unop,
    charMap_injective_of_addGroup (bHom V DK α.unop)⟩
  fsmIso {A B} α hα := by
    haveI : IsIso α := finSubOp_isOfFSMType B A α hα
    haveI : IsIso α.op := inferInstance
    exact ConcreteCategory.bijective_of_isIso _

omit [∀ L : FinSub V.functionField Kbar, CompactSpace (normObj V L)] in
/-- ★**`B(L)` は group-like**。 -/
theorem bmonGeom_isGroupLike (A : (FinSub V.functionField Kbar)ᵒᵖ) :
    IsGroupLike ((bmonGeom V DK).val A) := by
  show IsGroupLike (Additive (BSubgroup V DK A.unop (normObj_isNormalScheme V A.unop)))
  refine (isGroupLike_iff _).mpr (fun a => ?_)
  exact ⟨⟨a, -a, by simp, by simp⟩, rfl⟩

/-! ## ★3. `div : B(L) → Φ(L)^gp` -/

/-- ★**主因子は Cartier** —— 局所方程式は `f` そのもの。 -/
theorem isCartierDiv_weilDivOfFn {X : Scheme.{u}} [IsIntegral X] [IsLocallyNoetherian X]
    (hnorm : IsNormalScheme X) [CompactSpace X] {f : X.functionField} (hf : f ≠ 0) :
    IsCartierDiv hnorm (weilDivOfFn hnorm hf) :=
  fun _ => ⟨⊤, trivial, f, hf, fun _ _ => rfl⟩

/-- ★★**`div : B(L) → Φ(L)^gp`**(束ねた型の上で)。 -/
noncomputable def divBHom (L : FinSub V.functionField Kbar) :
    Additive (BSubgroup V DK L (normObj_isNormalScheme V L))
      →+ cartierOnDL V DK L (normObj_isNormalScheme V L) :=
  AddMonoidHom.mk'
    (fun x => ⟨Finsupp.subtypeDomain (· ∈ DLSet V DK L)
        (weilDivOfFn (normObj_isNormalScheme V L)
          (((Additive.toMul x : BSubgroup V DK L _)
            : ((normObj V L).functionField)ˣ)).ne_zero), by
      show toWeilOnDL V DK L _ ∈ cartierSubgroup (normObj_isNormalScheme V L)
      rw [embDomain_subtypeDomain V DK L _ (fun w hw => by
        rw [weilDivOfFn_apply]
        exact (Additive.toMul x).2 w hw)]
      exact isCartierDiv_weilDivOfFn (normObj_isNormalScheme V L) _⟩)
    (fun x y => by
      refine Subtype.ext (Finsupp.ext fun s => ?_)
      show ordPt (normObj V L) (normObj_isNormalScheme V L) s.1
          (((Additive.toMul (x + y) : BSubgroup V DK L _)
            : ((normObj V L).functionField)ˣ) : (normObj V L).functionField)
        = ordPt (normObj V L) (normObj_isNormalScheme V L) s.1
            (((Additive.toMul x : BSubgroup V DK L _)
              : ((normObj V L).functionField)ˣ) : (normObj V L).functionField)
          + ordPt (normObj V L) (normObj_isNormalScheme V L) s.1
            (((Additive.toMul y : BSubgroup V DK L _)
              : ((normObj V L).functionField)ˣ) : (normObj V L).functionField)
      exact ordPt_mul (normObj_isNormalScheme V L) s.1
        (Additive.toMul x : BSubgroup V DK L _).1.ne_zero
        (Additive.toMul y : BSubgroup V DK L _).1.ne_zero)

omit [IsGalois V.functionField Kbar] in
@[simp] theorem divBHom_apply (L : FinSub V.functionField Kbar)
    (x : Additive (BSubgroup V DK L (normObj_isNormalScheme V L)))
    (s : (DLSet V DK L : Type u)) :
    ((divBHom V DK L x : cartierOnDL V DK L _) : (DLSet V DK L) →₀ ℤ) s
      = ordPt (normObj V L) (normObj_isNormalScheme V L) s.1
        (((Additive.toMul x : BSubgroup V DK L _)
          : ((normObj V L).functionField)ˣ) : (normObj V L).functionField) := rfl

/-! ## ★4. `div` の関手性と model Frobenioid -/

variable (hkq : IsKQCartier V DK
  (fun (L : FinSub V.functionField Kbar) _ => normObj_isNormalScheme V L))

/-- ★★**`div : B(L) → Φ(L)^gp`**。 -/
noncomputable def divBGeom (A : (FinSub V.functionField Kbar)ᵒᵖ) :
    (bmonGeom V DK).val A →+
      Gp (((cartierDatumGeom V DK hkq).phi finSubOp_isOfFSMType).val A) :=
  ((effSubGpEquiv ((cartierDatumGeom V DK hkq).grp A)
      ((cartierDatumGeom V DK hkq).qc A)).symm.toAddMonoidHom).comp
    (divBHom V DK A.unop)

theorem phiGpHomC_divBGeom (A : (FinSub V.functionField Kbar)ᵒᵖ)
    (x : (bmonGeom V DK).val A) :
    phiGpHomC (cartierDatumGeom V DK hkq) finSubOp_isOfFSMType (divBGeom V DK hkq A x)
      = divBHom V DK A.unop x := by
  show (effSubGpEquiv ((cartierDatumGeom V DK hkq).grp A)
      ((cartierDatumGeom V DK hkq).qc A))
      ((effSubGpEquiv ((cartierDatumGeom V DK hkq).grp A)
        ((cartierDatumGeom V DK hkq).qc A)).symm (divBHom V DK A.unop x)) = _
  rw [AddEquiv.apply_symm_apply]


omit [IsGalois V.functionField Kbar] in
/-- ★★★**`div` は引き戻しと両立する** —— `u` 自身が大域の局所方程式だから。 -/
theorem divBHom_bHom {L M : FinSub V.functionField Kbar} (f : L ⟶ M)
    (x : Additive (BSubgroup V DK L (normObj_isNormalScheme V L))) :
    divBHom V DK M (bHom V DK f x) = phiPull V DK f (divBHom V DK L x) := by
  refine Subtype.ext (Finsupp.ext fun s => ?_)
  have hq : (((Additive.toMul x : BSubgroup V DK L _)
      : ((normObj V L).functionField)ˣ) : (normObj V L).functionField) ≠ 0 :=
    (Additive.toMul x : BSubgroup V DK L _).1.ne_zero
  have hDU : ∀ v : PrimeDivisorPt (normObj V L), (v : normObj V L) ∈ (⊤ : (normObj V L).Opens) →
      toWeilOnDL V DK L ((divBHom V DK L x : cartierOnDL V DK L _)
          : (DLSet V DK L) →₀ ℤ) v
        = ordPt (normObj V L) (normObj_isNormalScheme V L) v
          (((Additive.toMul x : BSubgroup V DK L _)
            : ((normObj V L).functionField)ˣ) : (normObj V L).functionField) := by
    intro v _
    rw [show ((divBHom V DK L x : cartierOnDL V DK L _) : (DLSet V DK L) →₀ ℤ)
        = Finsupp.subtypeDomain (· ∈ DLSet V DK L)
          (weilDivOfFn (normObj_isNormalScheme V L) hq) from rfl,
      embDomain_subtypeDomain V DK L _ (fun w hw => by
        rw [weilDivOfFn_apply]
        exact (Additive.toMul x).2 w hw)]
    rfl
  show ordPt (normObj V M) (normObj_isNormalScheme V M) s.1
      (((Additive.toMul (bHom V DK f x) : BSubgroup V DK M _)
        : ((normObj V M).functionField)ˣ) : (normObj V M).functionField)
    = pullCoeff (normMap V f) (normObj_isNormalScheme V M)
      (normObj_isNormalScheme V L) (divBHom V DK L x).2 s.1
  rw [pullCoeff_eq (normMap V f) (normObj_isNormalScheme V M) (normObj_isNormalScheme V L)
    (divBHom V DK L x).2 s.1 (normMap_hdim f s.1) (U := ⊤) hq (Set.mem_univ _) hDU]
  rw [bHom_val]
  rfl

/-- ★★★★**`div` の関手性**。 -/
theorem divBGeom_nat {A B : (FinSub V.functionField Kbar)ᵒᵖ} (f : A ⟶ B)
    (x : (bmonGeom V DK).val B) :
    divBGeom V DK hkq A ((bmonGeom V DK).map f x)
      = ((cartierDatumGeom V DK hkq).phi finSubOp_isOfFSMType).gpMapOn f
        (divBGeom V DK hkq B x) := by
  refine phiGpHomC_injective (cartierDatumGeom V DK hkq) finSubOp_isOfFSMType ?_
  rw [phiGpHomC_divBGeom, phiGpHomC_gpMapOn, phiGpHomC_divBGeom]
  exact divBHom_bHom V DK f.unop x

/-! ## ★5. model Frobenioid -/

/-- ★★★★**`Example 6.1` の `ModelData`** —— `Theorem 5.2` の入力そのもの。 -/
noncomputable def ex61ModelData : ModelData.{u, u, u} ((FinSub V.functionField Kbar)ᵒᵖ) :=
  (cartierDatumGeom V DK hkq).modelData finSubOp_isOfFSMType (bmonGeom V DK)
    (fun A => divBGeom V DK hkq A) (fun f x => divBGeom_nat V DK hkq f x)

@[simp] theorem ex61ModelData_phi :
    (ex61ModelData V DK hkq).phi
      = (cartierDatumGeom V DK hkq).phi finSubOp_isOfFSMType := rfl

@[simp] theorem ex61ModelData_bmon :
    (ex61ModelData V DK hkq).bmon = bmonGeom V DK := rfl

/-- ★★★★★★**[FrdI] `Example 6.1`** —— 幾何のデータから **Frobenioid** ができる。

原文 (FrdI p.109):
> If for every Spec(L) ∈Ob(B(G)0) [cf. §0], every prime divisor of DL is Q-Cartier, -/
theorem ex61Frobenioid :
    Frobenioid (ModelData.modelPre
      ((cartierDatumGeom V DK hkq).modelHyp finSubOp_isOfFSMType (bmonGeom V DK)
        (fun A => divBGeom V DK hkq A) (fun f x => divBGeom_nat V DK hkq f x)
        (bmonGeom_isGroupLike V DK) finSubOp_totallyEpimorphic finSubOp_isConnected)) :=
  (cartierDatumGeom V DK hkq).cartierFrobenioid finSubOp_isOfFSMType (bmonGeom V DK)
    (fun A => divBGeom V DK hkq A) (fun f x => divBGeom_nat V DK hkq f x)
    (bmonGeom_isGroupLike V DK) finSubOp_totallyEpimorphic finSubOp_isConnected

/-- ★★★★★**`Example 6.1`** —— できた Frobenioid は **isotropic 型**。 -/
theorem ex61Frobenioid_isotropicType :
    IsOfIsotropicType (ModelData.modelPre
      ((cartierDatumGeom V DK hkq).modelHyp finSubOp_isOfFSMType (bmonGeom V DK)
        (fun A => divBGeom V DK hkq A) (fun f x => divBGeom_nat V DK hkq f x)
        (bmonGeom_isGroupLike V DK) finSubOp_totallyEpimorphic finSubOp_isConnected)) :=
  (cartierDatumGeom V DK hkq).cartierFrobenioid_isotropicType finSubOp_isOfFSMType
    (bmonGeom V DK) (fun A => divBGeom V DK hkq A) (fun f x => divBGeom_nat V DK hkq f x)
    (bmonGeom_isGroupLike V DK) finSubOp_totallyEpimorphic finSubOp_isConnected

/-! ### ★出典の紐付け -/

/-- ★★★★★★locator —— `Example 6.1` の model Frobenioid の構成
(★**条つき**: 大域単数 `𝒪^×(A) = k_L^×` は節点 `ex61-units` に残る)。 -/
def ex61Frobenioid.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 109,
    item := "Example 6.1 — 幾何のデータから model Frobenioid ができる",
    sectionId := "frdi-example-6-1" }

end Geom

end ABC3.Found.Divisor
