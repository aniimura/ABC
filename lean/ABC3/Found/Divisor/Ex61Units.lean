/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Divisor.SchemeHartogs
import ABC3.Found.Divisor.Ex61Model
import ABC3.Found.FrdI.Thm52Otimes

/-!
# [FrdI] Example 6.1 —— `𝒪^×(A) = 𝒪^▷(A) = k_L^×`(鎖 `normalize` の `ex61-units`)

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.110。

原文 (FrdI p.110):
> that [since V [L] is a proper normal variety] for A ∈Ob(CV,

原文 (FrdI p.110):
> where kL denotes the algebraic closure of k in L [so kL is a finite separable extension

## ★★組み立て

| 段 | 実装 |
|---|---|
| model Frobenioid で `𝒪^×(A) ≃ {u ∈ B(L) \| Div_B(u) = 0}` | `Found/FrdI/Thm52Otimes.lean` |
| `Div_B(u) = 0` ＋ `u ∈ B(L)` ⟹ **すべての**余次元 1 の点で `ord = 0` | ★本ファイル |
| `ord = 0` がすべて ⟹ `Γ(V[L],⊤)` の単元 | `Found/Divisor/SchemeHartogs.lean` |
| `Γ(V[L],⊤)` は `k` 上有限次の体(＝ `k_L`) | `Found/Divisor/GlobalUnits.lean` |

## ★★★2 段目が本ファイルの実質

`B(L)` の定義は「`D_L` の**外**では `ord = 0`」(`BSubgroup`)であり、
`Div_B(u) = 0` は「`D_L` の**中**で `ord = 0`」である。
★合わせて **`V[L]` のすべての余次元 1 の点で `ord = 0`** になる。

## ★逆向き

`Γ(V[L],⊤)` の単元は各茎で単元なので、どの余次元 1 の点でも `ord = 0`。
★したがって `B(L)` に入り、`Div_B` も `0` になる。
-/

namespace ABC3.Found.Divisor

open AlgebraicGeometry CategoryTheory ABC3.Found.FrdI ABC3.Meta

universe u

/-! ## ★1. 大域単元は `ord = 0` -/

/-- ★**大域切断の芽を函数体へ送る 2 通りの道は一致する**。 -/
theorem germ_algebraMap_functionField {X : Scheme.{u}} [IsIntegral X] (U : X.Opens)
    (w : X) (hw : w ∈ U) [Nonempty U] (t : Γ(X, U)) :
    algebraMap (X.presheaf.stalk w) X.functionField
        ((ConcreteCategory.hom (X.presheaf.germ U w hw)) t)
      = (ConcreteCategory.hom (X.germToFunctionField U)) t := by
  show (ConcreteCategory.hom (X.presheaf.stalkSpecializes _))
    ((ConcreteCategory.hom (X.presheaf.germ U w hw)) t) = _
  rw [← CommRingCat.comp_apply, X.presheaf.germ_stalkSpecializes]

/-- ★★**大域切断の単元はどの余次元 1 の点でも `ord = 0`**。 -/
theorem ordPt_eq_zero_of_globalUnit {X : Scheme.{u}} [IsIntegral X] [IsLocallyNoetherian X]
    (hnorm : IsNormalScheme X) (t : Γ(X, ⊤)) (ht : IsUnit t) (w : PrimeDivisorPt X) :
    ordPt X hnorm w ((ConcreteCategory.hom (X.germToFunctionField ⊤)) t) = 0 := by
  haveI : Nonempty ((⊤ : X.Opens)) := ⟨⟨(IsIntegral.nonempty (X := X)).some, trivial⟩⟩
  have hu : IsUnit ((ConcreteCategory.hom (X.presheaf.germ ⊤ w.1 trivial)) t) :=
    ht.map (X.presheaf.germ ⊤ w.1 trivial).hom
  have hgm := germ_algebraMap_functionField (X := X) ⊤ w.1 trivial t
  rw [← hgm, ← hu.unit_spec]
  exact ordPt_eq_zero_of_isUnit hnorm w hu.unit

/-- ★★★**`ord ≥ 0` がすべてなら `ord = 0` がすべて**(proper normal のとき)。

★Hartogs で `Γ(X,⊤)` に入り、`Γ(X,⊤)` は体なので単元、単元なら `ord = 0`。
★★これが `𝒪^▷(A) = 𝒪^×(A)`(可逆性を課さなくても可逆になる)の中身である。 -/
theorem forall_ordPt_eq_zero_of_forall_nonneg {X : Scheme.{u}} [IsIntegral X]
    [IsLocallyNoetherian X] {k : Type u} [Field k]
    (g : X ⟶ Spec (CommRingCat.of k)) (hg : IsProper g)
    (hnorm : IsNormalScheme X)
    (hnoeth : ∀ U : X.Opens, IsAffineOpen U → Nonempty U → IsNoetherianRing Γ(X, U))
    (hic : ∀ U : X.Opens, IsAffineOpen U → Nonempty U → IsIntegrallyClosed Γ(X, U))
    (u : X.functionField) (hu : u ≠ 0)
    (h : ∀ x : PrimeDivisorPt X, 0 ≤ ordPt X hnorm x u) :
    ∀ x : PrimeDivisorPt X, ordPt X hnorm x u = 0 := by
  obtain ⟨t, ht⟩ := mem_range_germ_top_of_forall_ordPt_nonneg hnorm hnoeth hic u h
  letI := (globalSections_isField g hg).toField
  have htne : t ≠ 0 := fun h0 => hu (by rw [← ht, h0, map_zero])
  intro x
  rw [← ht]
  exact ordPt_eq_zero_of_globalUnit hnorm t (isUnit_iff_ne_zero.mpr htne) x

/-! ## ★2. `Div_B(u) = 0` は「すべての余次元 1 の点で `ord = 0`」 -/

variable (V : Scheme.{u}) [IsIntegral V] {Kbar : Type u} [Field Kbar]
  [Algebra V.functionField Kbar]
  (DK : Set (PrimeDivisorPt V))
  [∀ L : FinSub V.functionField Kbar, IsLocallyNoetherian (normObj V L)]
  [∀ L : FinSub V.functionField Kbar, CompactSpace (normObj V L)]

/-- ★★★**`B(L)` の元で `Div_B = 0` なら、`V[L]` のすべての余次元 1 の点で `ord = 0`**。

★`B(L)` は `D_L` の**外**、`Div_B = 0` は `D_L` の**中**を押さえる。 -/
theorem forall_ordPt_eq_zero_of_divBHom_eq_zero (L : FinSub V.functionField Kbar)
    (x : Additive (BSubgroup V DK L (normObj_isNormalScheme V L)))
    (h : divBHom V DK L x = 0) :
    ∀ w : PrimeDivisorPt (normObj V L),
      ordPt (normObj V L) (normObj_isNormalScheme V L) w
        (((Additive.toMul x : BSubgroup V DK L _) : ((normObj V L).functionField)ˣ)
          : (normObj V L).functionField) = 0 := by
  intro w
  by_cases hw : w ∈ DLSet V DK L
  · have hval := divBHom_apply V DK L x ⟨w, hw⟩
    rw [h] at hval
    simpa using hval.symm
  · exact (Additive.toMul x : BSubgroup V DK L _).2 w hw

/-! ## ★3. 出口 —— `𝒪^×(A) ⊆ k_L^×` -/

/-- ★★★★★★**[FrdI] Example 6.1** —— `Div_B(u) = 0` なる `u ∈ B(L)` は
`Γ(V[L],⊤)` の**単元**から来る。

★★`Γ(V[L],⊤)` は `k` 上有限次の体(`globalSections_isField` / `globalSections_finite`)
なので、これが原文の `𝒪^×(A) = 𝒪^▷(A) = k_L^×` の中身である。 -/
theorem exists_globalUnit_of_divBHom_eq_zero (L : FinSub V.functionField Kbar)
    {k : Type u} [Field k]
    (g : normObj V L ⟶ Spec (CommRingCat.of k)) (hg : IsProper g)
    (hnoeth : ∀ U : (normObj V L).Opens, IsAffineOpen U → Nonempty U →
      IsNoetherianRing Γ(normObj V L, U))
    (hic : ∀ U : (normObj V L).Opens, IsAffineOpen U → Nonempty U →
      IsIntegrallyClosed Γ(normObj V L, U))
    (x : Additive (BSubgroup V DK L (normObj_isNormalScheme V L)))
    (h : divBHom V DK L x = 0) :
    ∃ t : Γ(normObj V L, ⊤), IsUnit t ∧
      (ConcreteCategory.hom ((normObj V L).germToFunctionField ⊤)) t
        = (((Additive.toMul x : BSubgroup V DK L _) : ((normObj V L).functionField)ˣ)
            : (normObj V L).functionField) :=
  exists_unit_of_forall_ordPt_eq_zero g hg (normObj_isNormalScheme V L) hnoeth hic _
    ((Additive.toMul x : BSubgroup V DK L _) : ((normObj V L).functionField)ˣ).ne_zero
    (forall_ordPt_eq_zero_of_divBHom_eq_zero V DK L x h)

/-! ## ★4. 逆向き —— `k_L^× ⊆ 𝒪^×(A)` -/

omit [∀ L : FinSub V.functionField Kbar, CompactSpace (normObj V L)] in
/-- ★★**`Γ(V[L],⊤)` の単元は `B(L)` に入る**。 -/
theorem globalUnit_mem_BSubgroup (L : FinSub V.functionField Kbar)
    (t : Γ(normObj V L, ⊤)) (ht : IsUnit t)
    (hne : (ConcreteCategory.hom ((normObj V L).germToFunctionField ⊤)) t ≠ 0) :
    Units.mk0 _ hne ∈ BSubgroup V DK L (normObj_isNormalScheme V L) := by
  intro w _
  exact ordPt_eq_zero_of_globalUnit (normObj_isNormalScheme V L) t ht w

/-- ★★★**`Γ(V[L],⊤)` の単元は `Div_B = 0`**。 -/
theorem divBHom_eq_zero_of_globalUnit (L : FinSub V.functionField Kbar)
    (t : Γ(normObj V L, ⊤)) (ht : IsUnit t)
    (hne : (ConcreteCategory.hom ((normObj V L).germToFunctionField ⊤)) t ≠ 0) :
    divBHom V DK L
        (Additive.ofMul ⟨Units.mk0 _ hne, globalUnit_mem_BSubgroup V DK L t ht hne⟩) = 0 := by
  refine Subtype.ext (Finsupp.ext fun s => ?_)
  rw [divBHom_apply V DK L]
  show ordPt (normObj V L) (normObj_isNormalScheme V L) s.1
      ((ConcreteCategory.hom ((normObj V L).germToFunctionField ⊤)) t) = _
  rw [ordPt_eq_zero_of_globalUnit (normObj_isNormalScheme V L) t ht s.1]
  rfl

/-! ## ★5. `𝒪^▷(A) = 𝒪^×(A)`

原文 (FrdI p.110):
> where kL denotes the algebraic closure of k in L [so kL is a finite separable extension
-/

section Otri

variable [IsGalois V.functionField Kbar]
  (hkq : IsKQCartier V DK
    (fun (L : FinSub V.functionField Kbar) _ => normObj_isNormalScheme V L))

/-- ★`Example 6.1` の model Frobenioid の仮説束。 -/
theorem ex61Hyp : (ex61ModelData V DK hkq).Hyp :=
  (cartierDatumGeom V DK hkq).modelHyp finSubOp_isOfFSMType (bmonGeom V DK)
    (fun A => divBGeom V DK hkq A) (fun f x => divBGeom_nat V DK hkq f x)
    (bmonGeom_isGroupLike V DK) finSubOp_totallyEpimorphic finSubOp_isConnected

open ModelData in
/-- ★★★★★★**[FrdI] Example 6.1 —— `𝒪^▷(A) = 𝒪^×(A)`**。

★`𝒪^▷(A)` は可逆性を課さないので `Div` は有効因子でしかないが、
`Div = Div_B(u)` と `u ∈ B(L)` から**すべての**余次元 1 の点で `ord(u) ≥ 0` が出る。
`V[L]` が proper normal なので(Hartogs ＋ `Γ(V[L],⊤)` が体)`u` は単元、
したがって `ord(u) = 0`、すなわち `Div = 0` で `φ` は可逆になる。

★★これで原文の `𝒪^×(A) = 𝒪^▷(A) = k_L^×` の 3 つの等号が揃う。 -/
theorem ex61_otri_le_otimes (A : ModelData.Obj (ex61ModelData V DK hkq))
    {k : Type u} [Field k]
    (g : normObj V A.base.unop ⟶ Spec (CommRingCat.of k)) (hg : IsProper g)
    (hnoeth : ∀ U : (normObj V A.base.unop).Opens, IsAffineOpen U → Nonempty U →
      IsNoetherianRing Γ(normObj V A.base.unop, U))
    (hic : ∀ U : (normObj V A.base.unop).Opens, IsAffineOpen U → Nonempty U →
      IsIntegrallyClosed Γ(normObj V A.base.unop, U)) :
    OTri (modelPre (ex61Hyp V DK hkq)) A ≤ OTimes (modelPre (ex61Hyp V DK hkq)) A := by
  rintro ψ hψ
  -- ★段 1: `Div = Div_B(u)` を `Γ` の中で
  have hdiv : effSubIncl ((cartierDatumGeom V DK hkq).grp A.base) (ψ.div)
      = divBHom V DK A.base.unop (ψ.u) := by
    have hkey : toGp (((cartierDatumGeom V DK hkq).phi finSubOp_isOfFSMType).val A.base)
        (ψ.div) = divBGeom V DK hkq A.base (ψ.u) := by
      have h0 := ModelData.toGp_div_eq_divB_of_mem_otri (ex61Hyp V DK hkq) ⟨ψ, hψ⟩
      rw [toGpHom_apply] at h0
      exact h0
    have happ := congrArg (phiGpHomC (cartierDatumGeom V DK hkq) finSubOp_isOfFSMType) hkey
    rw [phiGpHomC_divBGeom, phiGpHomC_toGp'] at happ
    exact happ
  -- ★段 2: すべての余次元 1 の点で `ord ≥ 0`
  have hnn : ∀ w : PrimeDivisorPt (normObj V A.base.unop),
      0 ≤ ordPt (normObj V A.base.unop) (normObj_isNormalScheme V A.base.unop) w
        (((Additive.toMul ψ.u : BSubgroup V DK A.base.unop _)
          : ((normObj V A.base.unop).functionField)ˣ)
            : (normObj V A.base.unop).functionField) := by
    intro w
    by_cases hw : w ∈ DLSet V DK A.base.unop
    · have h1 := divBHom_apply V DK A.base.unop (ψ.u) ⟨w, hw⟩
      rw [← hdiv] at h1
      rw [← h1]
      exact ψ.div.2.2 ⟨w, hw⟩
    · exact le_of_eq ((Additive.toMul ψ.u : BSubgroup V DK A.base.unop _).2 w hw).symm
  -- ★段 3: proper normal なので `ord = 0`
  have hz := forall_ordPt_eq_zero_of_forall_nonneg g hg
    (normObj_isNormalScheme V A.base.unop) hnoeth hic _
    (((Additive.toMul ψ.u : BSubgroup V DK A.base.unop _)
      : ((normObj V A.base.unop).functionField)ˣ)).ne_zero hnn
  -- ★段 4: `Div = 0` ゆえ可逆
  have hdb : divBHom V DK A.base.unop (ψ.u) = 0 := by
    refine Subtype.ext (Finsupp.ext fun s => ?_)
    rw [divBHom_apply V DK A.base.unop, hz s.1]
    rfl
  have h2 : effSubIncl ((cartierDatumGeom V DK hkq).grp A.base) (ψ.div) = 0 := by
    rw [hdiv]; exact hdb
  rw [Subtype.ext_iff] at h2
  exact ModelData.mem_otimes_of_div_eq_zero _ ψ hψ.1 hψ.2 (Subtype.ext h2)

end Otri

/-! ### ★出典の紐付け(`.src`)と、証明が要求するもの(`.needs`) -/

def exists_globalUnit_of_divBHom_eq_zero.src : Source :=
  { paper := "FrdI", pdfPage := 110,
    item := "Example 6.1 — O^×(A) = O^▷(A) = k_L^×",
    sectionId := "frdi-example-6-1" }

def exists_globalUnit_of_divBHom_eq_zero.needs : List ProofObligation :=
  [ .citation "[ABC3]" "model Frobenioid の 𝒪^×(A) は Div_B の核"
      (.inProject "ABC3" "ABC3.Found.FrdI.ModelData.otimesModelEquiv") 110,
    .citation "[ABC3]" "代数的 Hartogs のスキーム版 ＋ Γ(X,⊤) は体"
      (.inProject "ABC3" "ABC3.Found.Divisor.exists_unit_of_forall_ordPt_eq_zero") 110,
    .citation "[ABC3]" "B(L) の定義(D_L の外では ord = 0)"
      (.inProject "ABC3" "ABC3.Found.Divisor.BSubgroup") 110,
    .derivation "B(L) が D_L の外を、Div_B = 0 が D_L の中を押さえるので、全ての余次元 1 の点で ord = 0" 110,
    .implicitStep
      "★原文は「we observe that [since V [L] is a proper normal variety]」の 1 文で畳む" 110 ]

def divBHom_eq_zero_of_globalUnit.src : Source :=
  { paper := "FrdI", pdfPage := 110,
    item := "Example 6.1 — k_L^× ⊆ O^×(A)",
    sectionId := "frdi-example-6-1" }

def divBHom_eq_zero_of_globalUnit.needs : List ProofObligation :=
  [ .citation "[ABC3]" "ordPt_eq_zero_of_isUnit(茎の単元は ord = 0)"
      (.inProject "ABC3" "ABC3.Found.Divisor.ordPt_eq_zero_of_isUnit") 110,
    .derivation "大域切断の単元は各茎で単元なので、どの余次元 1 の点でも ord = 0" 110 ]

def ordPt_eq_zero_of_globalUnit.src : Source :=
  { paper := "FrdI", pdfPage := 110,
    item := "Example 6.1 — 大域単元はどの余次元 1 の点でも ord = 0",
    sectionId := "frdi-example-6-1" }

def ordPt_eq_zero_of_globalUnit.needs : List ProofObligation :=
  [ .citation "[ABC3]" "ordPt_eq_zero_of_isUnit"
      (.inProject "ABC3" "ABC3.Found.Divisor.ordPt_eq_zero_of_isUnit") 110,
    .citation "[mathlib]" "TopCat.Presheaf.germ_stalkSpecializes"
      (.inMathlib "TopCat.Presheaf.germ_stalkSpecializes") 110 ]

def ex61_otri_le_otimes.src : Source :=
  { paper := "FrdI", pdfPage := 110,
    item := "Example 6.1 — O^▷(A) = O^×(A)",
    sectionId := "frdi-example-6-1" }

def ex61_otri_le_otimes.needs : List ProofObligation :=
  [ .citation "[ABC3]" "toGp_div_eq_divB_of_mem_otri(𝒪^▷ の元は Div = Div_B(u))"
      (.inProject "ABC3" "ABC3.Found.FrdI.ModelData.toGp_div_eq_divB_of_mem_otri") 110,
    .citation "[ABC3]" "forall_ordPt_eq_zero_of_forall_nonneg(Hartogs ＋ Γ(X,⊤) は体)"
      (.inProject "ABC3" "ABC3.Found.Divisor.forall_ordPt_eq_zero_of_forall_nonneg") 110,
    .citation "[ABC3]" "mem_otimes_of_div_eq_zero(Div = 0 なら可逆)"
      (.inProject "ABC3" "ABC3.Found.FrdI.ModelData.mem_otimes_of_div_eq_zero") 110,
    .derivation "Div が有効因子なので ord ≥ 0、proper normal なら単元、よって ord = 0 で Div = 0" 110 ]

def forall_ordPt_eq_zero_of_forall_nonneg.src : Source :=
  { paper := "FrdI", pdfPage := 110,
    item := "Example 6.1 — ord ≥ 0 がすべてなら ord = 0 がすべて",
    sectionId := "frdi-example-6-1" }

def forall_ordPt_eq_zero_of_forall_nonneg.needs : List ProofObligation :=
  [ .citation "[ABC3]" "mem_range_germ_top_of_forall_ordPt_nonneg(Hartogs のスキーム版)"
      (.inProject "ABC3" "ABC3.Found.Divisor.mem_range_germ_top_of_forall_ordPt_nonneg") 110,
    .citation "[ABC3]" "globalSections_isField"
      (.inProject "ABC3" "ABC3.Found.Divisor.globalSections_isField") 110 ]

def ex61Hyp.src : Source :=
  { paper := "FrdI", pdfPage := 109,
    item := "Example 6.1 — model Frobenioid の仮説束",
    sectionId := "frdi-example-6-1" }

def ex61Hyp.needs : List ProofObligation :=
  [ .citation "[ABC3]" "cartierDatumGeom(幾何から CartierDatum)"
      (.inProject "ABC3" "ABC3.Found.Divisor.cartierDatumGeom") 109 ]

end ABC3.Found.Divisor
