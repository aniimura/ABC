/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Thm34Pre
import ABC3.Found.FrdI.Prop33UnTr
import Mathlib.Data.PNat.Factors

/-!
# [FrdI] Theorem 3.4, (ii) を quasi-isotropic 型へ

★現行の `thm_3_4_ii` は isotropic 型を仮定しているが、
原典は quasi-isotropic 型で十分だと言う。その帰着の材料を置く。
-/

universe v u w u2 v2

namespace ABC3.Found.FrdI

open CategoryTheory

/-! ## ★★★`Theorem 3.4, (ii)` を quasi-isotropic 型へ一般化する準備

原文 (FrdI p.62):
> (ii) Suppose that C1, C2 are of quasi-isotropic type, and that D1, D2 are of

★★現行の `thm_3_4_ii` は **isotropic 型**を仮定しているが、
原典は **quasi-isotropic 型**で十分だと言う。
★帰着の筋は `𝒞^istr` へ落とすことで、
\`psiIstr_isEquivalence\`(quasi-isotropic を仮定として受ける)が既にある。

★★**欠けていた環**は「`𝒞` での pre-step 性 ⇔ `𝒞^istr` での pre-step 性」であり、
★在庫の `isotropification_degFr`(次数の保存)と
`isotropification_baseIso_iff`(底同型の両向き保存)を組むだけで出る。 -/

/-- ★★★★**isotropification は pre-step 性を両向きに保つ**。

★`ii-quasi`(quasi-isotropic への一般化)の欠けていた環。 -/
theorem isotropification_isPreStep_iff {Dd : Type u} [Category.{v} Dd] {Cc : Type u2}
    [Category.{v2} Cc] {Φ₀ : MonoidOn.{v, u, w} Dd} (P : PreFrobenioid Cc Φ₀)
    (F : FrobenioidCore P) {A B : Cc} (f : A ⟶ B) :
    IsPreStep P f ↔ IsPreStep (istrPre P F) ((isotropification P F).map f) := by
  constructor
  · rintro ⟨hl, hb⟩
    refine ⟨?_, ?_⟩
    · show P.degFr (istrMap P F f) = 1
      rw [isotropification_degFr]
      exact hl
    · exact (isotropification_baseIso_iff P F f).mpr hb
  · rintro ⟨hl, hb⟩
    refine ⟨?_, ?_⟩
    · show P.degFr f = 1
      rw [← isotropification_degFr P F f]
      exact hl
    · exact (isotropification_baseIso_iff P F f).mp hb

def isotropification_isPreStep_iff.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 62,
    item := "Theorem 3.4, (ii) — isotropification は pre-step を保つ",
    sectionId := "frdi-thm-3-4" }

/-! ## ★同型四角形で pre-step 性は不変

★★`isotropificationCommute`(`Thm34.lean:774`)の自然同型を渡るのに要る。
★`Ψ^istr ∘ istr₁` と `istr₂ ∘ Ψ` は**同型だが同一ではない**ので、
射の性質を四角形を越えて運ぶ補題が必要である。 -/

/-- ★★★★**同型の四角形を越えて pre-step 性は保たれる**。

★次数は同型が 1 だから、底同型は同型の合成から出る。 -/
theorem isPreStep_congr_iso {Dd : Type u} [Category.{v} Dd] {Cc : Type u2}
    [Category.{v2} Cc] {Φ₀ : MonoidOn.{v, u, w} Dd} (P : PreFrobenioid Cc Φ₀)
    {A B A' B' : Cc} (f : A ⟶ B) (f' : A' ⟶ B')
    (α : A ≅ A') (β : B ≅ B') (h : f ≫ β.hom = α.hom ≫ f') :
    IsPreStep P f ↔ IsPreStep P f' := by
  have hda : P.degFr α.hom = 1 := degFr_of_isIso P α.hom
  have hdb : P.degFr β.hom = 1 := degFr_of_isIso P β.hom
  have hd : P.degFr f = P.degFr f' := by
    have := congrArg P.degFr h
    rwa [P.degFr_comp, P.degFr_comp, hda, hdb, mul_one, one_mul] at this
  have hb := congrArg P.Base h
  rw [P.Base_comp, P.Base_comp] at hb
  haveI : IsIso (P.Base α.hom) := ⟨P.Base α.inv, by
    rw [← P.Base_comp, α.hom_inv_id, P.Base_id], by
    rw [← P.Base_comp, α.inv_hom_id, P.Base_id]⟩
  haveI : IsIso (P.Base β.hom) := ⟨P.Base β.inv, by
    rw [← P.Base_comp, β.hom_inv_id, P.Base_id], by
    rw [← P.Base_comp, β.inv_hom_id, P.Base_id]⟩
  constructor
  · rintro ⟨hl, hbi⟩
    refine ⟨?_, ?_⟩
    · show P.degFr f' = 1
      rw [← hd]; exact hl
    · show IsIso (P.Base f')
      haveI : IsIso (P.Base f) := hbi
      have : P.Base f' = inv (P.Base α.hom) ≫ (P.Base f ≫ P.Base β.hom) := by
        rw [hb, ← Category.assoc, IsIso.inv_hom_id, Category.id_comp]
      rw [this]; infer_instance
  · rintro ⟨hl, hbi⟩
    refine ⟨?_, ?_⟩
    · show P.degFr f = 1
      rw [hd]; exact hl
    · show IsIso (P.Base f)
      haveI : IsIso (P.Base f') := hbi
      have : P.Base f = (P.Base α.hom ≫ P.Base f') ≫ inv (P.Base β.hom) := by
        rw [← hb, Category.assoc, IsIso.hom_inv_id, Category.comp_id]
      rw [this]; infer_instance

def isPreStep_congr_iso.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 62,
    item := "Theorem 3.4, (ii) — 同型四角形で pre-step 性は不変",
    sectionId := "frdi-thm-3-4" }

/-! ## ★★★★`Theorem 3.4, (ii)` の quasi-isotropic 版

原文 (FrdI p.62):
> (ii) Suppose that C1, C2 are of quasi-isotropic type, and that D1, D2 are of

★★**帰着で出す**: `𝒞^istr` へ落として既存の isotropic 版を使う。 -/

variable {D₁ : Type u} [Category.{v} D₁] {C₁ : Type u2} [Category.{v2} C₁]
  {Φ₁ : MonoidOn.{v, u, w} D₁} {P₁ : PreFrobenioid C₁ Φ₁}
  {D₂ : Type u} [Category.{v} D₂] {C₂ : Type u2} [Category.{v2} C₂]
  {Φ₂ : MonoidOn.{v, u, w} D₂} {P₂ : PreFrobenioid C₂ Φ₂}

/-- ★★★★★**[FrdI] Theorem 3.4, (ii)** の **quasi-isotropic 版**(pre-step の条)。

★手: `isotropification_isPreStep_iff` で `𝒞^istr` へ落とし、
`thm_3_4_ii` を `Ψ^istr` に適用し、
`isotropificationCommute` の自然同型を `isPreStep_congr_iso` で渡る。 -/
theorem thm_3_4_ii_quasi_preStep (Ψ : C₁ ⥤ C₂) [Ψ.IsEquivalence]
    (F₁ : FrobenioidCore P₁) (G₁ : Frobenioid P₁)
    (hq₁ : IsOfQuasiIsotropicType C₁ P₁) (hFSMFF₁ : IsOfFSMFFType D₁)
    (hFrobMono₁ : ∀ {X Y : Istr P₁} (ε : X ⟶ Y),
      IsFrobeniusType (istrPre P₁ F₁) ε → Mono ε)
    (hFrobFS₁ : ∀ {X Y : Istr P₁} (ε : X ⟶ Y),
      IsFrobeniusType (istrPre P₁ F₁) ε → IsFiberwiseSurjective ε)
    (F₂ : FrobenioidCore P₂) (G₂ : Frobenioid P₂)
    (hq₂ : IsOfQuasiIsotropicType C₂ P₂) (hFSMFF₂ : IsOfFSMFFType D₂)
    (hFrobMono₂ : ∀ {X Y : Istr P₂} (ε : X ⟶ Y),
      IsFrobeniusType (istrPre P₂ F₂) ε → Mono ε)
    (hFrobFS₂ : ∀ {X Y : Istr P₂} (ε : X ⟶ Y),
      IsFrobeniusType (istrPre P₂ F₂) ε → IsFiberwiseSurjective ε)
    {A B : C₁} (φ : A ⟶ B) :
    IsPreStep P₁ φ ↔ IsPreStep P₂ (Ψ.map φ) := by
  haveI := psiIstr_isEquivalence Ψ P₁ P₂ hq₁ hq₂
  have key : ∀ {X Y : Istr P₁} (f : X ⟶ Y),
      IsPreStep (istrPre P₁ F₁) f ↔
        IsPreStep (istrPre P₂ F₂) ((psiIstr Ψ P₁ P₂ hq₁ hq₂).map f) :=
    fun f => (thm_3_4_ii (istrPre P₁ F₁) (istrPre P₂ F₂)
      (psiIstr Ψ P₁ P₂ hq₁ hq₂).asEquivalence
      (istr_frobenioidCore P₁ F₁) (istr_frobenioid P₁ F₁ G₁)
      (istr_isotropic P₁ F₁) hFSMFF₁ hFrobMono₁ hFrobFS₁
      (istr_frobenioidCore P₂ F₂) (istr_frobenioid P₂ F₂ G₂)
      (istr_isotropic P₂ F₂) hFSMFF₂ hFrobMono₂ hFrobFS₂).1 f
  have hcom := isotropificationCommute Ψ P₁ P₂ F₁ F₂ hq₁ hq₂
  rw [isotropification_isPreStep_iff P₁ F₁ φ, isotropification_isPreStep_iff P₂ F₂ (Ψ.map φ)]
  refine (key ((isotropification P₁ F₁).map φ)).trans ?_
  exact isPreStep_congr_iso (istrPre P₂ F₂) _ _ (hcom.app A) (hcom.app B)
    (hcom.hom.naturality φ)

def thm_3_4_ii_quasi_preStep.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 62,
    item := "Theorem 3.4, (ii) — quasi-isotropic 型での pre-step の保存",
    sectionId := "frdi-thm-3-4" }

/-! ## ★★素数の全単射から `ℕ≥1 ≃* ℕ≥1` を作る

★★`Theorem 3.4, (iii)` の `Ψ_{N≥1}` の**拡張段**である。

原文 (FrdI p.62):
> such that Ψ maps morphisms of Frobenius degree d to morphisms of Frobenius

★`ℕ≥1` は素数の自由可換単系なので、素数の全単射がそのまま単系自己同型を与える。
★mathlib の `PNat.factorMultiset` / `PrimeMultiset` を使う。 -/

/-- ★★★★**素数の全単射 ⇒ `ℕ≥1` の単系自己同型**。 -/
def pnatMulEquivOfPrimeEquiv (σ : Nat.Primes ≃ Nat.Primes) : ℕ+ ≃* ℕ+ where
  toFun n := PrimeMultiset.prod (Multiset.map σ (PNat.factorMultiset n))
  invFun n := PrimeMultiset.prod (Multiset.map σ.symm (PNat.factorMultiset n))
  left_inv n := by
    show PrimeMultiset.prod (Multiset.map σ.symm (PNat.factorMultiset
      (PrimeMultiset.prod (Multiset.map σ (PNat.factorMultiset n))))) = n
    rw [PrimeMultiset.factorMultiset_prod, Multiset.map_map,
      show (σ.symm ∘ σ) = id from funext σ.symm_apply_apply, Multiset.map_id]
    exact PNat.prod_factorMultiset n
  right_inv n := by
    show PrimeMultiset.prod (Multiset.map σ (PNat.factorMultiset
      (PrimeMultiset.prod (Multiset.map σ.symm (PNat.factorMultiset n))))) = n
    rw [PrimeMultiset.factorMultiset_prod, Multiset.map_map,
      show (σ ∘ σ.symm) = id from funext σ.apply_symm_apply, Multiset.map_id]
    exact PNat.prod_factorMultiset n
  map_mul' m n := by
    show PrimeMultiset.prod (Multiset.map σ (PNat.factorMultiset (m * n)))
      = PrimeMultiset.prod (Multiset.map σ (PNat.factorMultiset m))
        * PrimeMultiset.prod (Multiset.map σ (PNat.factorMultiset n))
    rw [PNat.factorMultiset_mul]
    show PrimeMultiset.prod (Multiset.map σ
        (PNat.factorMultiset m + PNat.factorMultiset n)) = _
    have hma : Multiset.map σ (PNat.factorMultiset m + PNat.factorMultiset n)
        = Multiset.map σ (PNat.factorMultiset m) + Multiset.map σ (PNat.factorMultiset n) :=
      Multiset.map_add _ _ _
    rw [hma]
    exact PrimeMultiset.prod_add _ _

def pnatMulEquivOfPrimeEquiv.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 62,
    item := "Theorem 3.4, (iii) — Ψ_{N≥1} の拡張段",
    sectionId := "frdi-thm-3-4" }

end ABC3.Found.FrdI
