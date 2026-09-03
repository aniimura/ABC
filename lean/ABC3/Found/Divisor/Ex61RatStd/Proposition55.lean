/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Divisor.Ex61Units
import ABC3.Found.FrdI.Thm64RatStd
import ABC3.Found.Divisor.Ex61RatStd.Theorem62

/-!
# Ex61RatStd —— `[FrdI] Proposition 5.5` の分

☆もとの 1 枚を**条なしの項目ごと**に割ったものである（第 1458、案 a）。
★「1 ファイル = 1 ノード」を回復するための分割で、中身は動かしていない。
-/

namespace ABC3.Found.Divisor

open AlgebraicGeometry CategoryTheory ABC3.Found.FrdI ABC3.Meta
universe u
variable (V : Scheme.{u}) [IsIntegral V] {Kbar : Type u} [Field Kbar]
  [Algebra V.functionField Kbar]
  (DK : Set (PrimeDivisorPt V))
  [∀ L : FinSub V.functionField Kbar, IsLocallyNoetherian (normObj V L)]
  [∀ L : FinSub V.functionField Kbar, CompactSpace (normObj V L)]
  [IsGalois V.functionField Kbar]
  (hkq : IsKQCartier V DK
    (fun (L : FinSub V.functionField Kbar) _ => normObj_isNormalScheme V L))

/-- ★★★★★★**[FrdI] Proposition 5.5, (iii)** —— 幾何の **`𝒞^un-tr` も rationally standard**。

★一般形は `Thm64RatStd.lean` の `unTr_isOfRationallyStandardType`。 -/
theorem ex61_unTr_ratStd
    (hsupp : ∀ (A : (FinSub V.functionField Kbar)ᵒᵖ) (s : (DLSet V DK A.unop : Type u)),
      ∃ y : (bmonGeom V DK).val A,
        ((divBHom V DK A.unop y : cartierOnDL V DK A.unop _)
          : (DLSet V DK A.unop) →₀ ℤ) s ≠ 0)
    (s0 : (DLSet V DK (botSub V.functionField Kbar) : Type u)) :
    IsOfRationallyStandardType
      (unTrPre (ModelData.modelPre (ex61Hyp V DK hkq))
        (ModelData.model_frobenioid (ex61Hyp V DK hkq)).core)
      (unTr_frobenioid (ModelData.modelPre (ex61Hyp V DK hkq))
        (ModelData.model_frobenioid (ex61Hyp V DK hkq)).core
        (ModelData.model_frobenioid (ex61Hyp V DK hkq)))
      (ex61Iota V DK hkq) := by
  classical
  haveI := (ex61Hyp V DK hkq).connectedD
  set Y : (FinSub V.functionField Kbar)ᵒᵖ := Opposite.op (botSub V.functionField Kbar) with hY
  set M := ex61ModelData V DK hkq with hM
  set h := ex61Hyp V DK hkq with hh
  set X0 : ModelData.Obj M := ⟨Y, 0⟩ with hX0
  set Xi : Istr (ModelData.modelPre h) := ⟨X0, ModelData.model_isotropicType h X0⟩ with hXi
  set Xi2 : Istr (unTrPre (ModelData.modelPre h)
      (ModelData.model_frobenioid h).core) :=
    ⟨show UnTr (ModelData.modelPre h) from Xi,
      unTr_isotropic (ModelData.modelPre h) (ModelData.model_frobenioid h).core _⟩ with hXi2
  obtain ⟨y, hy⟩ := hsupp Y s0
  exact unTr_isOfRationallyStandardType h (ex61Iota V DK hkq) (ex61_hsp V DK hkq hsupp)
    (fun A => (isDivisorial_effSub ((cartierDatumGeom V DK hkq).grp A)).1.1)
    finSubOp_isOfFSMType
    (ex61_not_isOfGroupLikeType V DK hkq Y ⟨s0⟩)
    (ex61_standardType V DK hkq Y ⟨s0⟩)
    (show UnTr (unTrPre (ModelData.modelPre h)
      (ModelData.model_frobenioid h).core) from Xi2)
    (fun e => ex61_phi_map_bot_eq_id V DK hkq e)
    y
    (fun n hn => divBGeom_nsmul_ne_zero V DK hkq Y y s0 hy n hn)

/-! ### ★出典の紐付け(`.src`)と、証明が要求するもの(`.needs`) -/

def ex61_unTr_ratStd.src : Source :=
  { paper := "FrdI", pdfPage := 105,
    item := "Proposition 5.5, (iii) — 幾何の 𝒞^un-tr も rationally standard",
    sectionId := "frdi-prop-5-5" }

end ABC3.Found.Divisor
