/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop48Cond3
import ABC3.Found.FrdI.Rmk451

/-!
# [FrdI] Proposition 4.8, (iii) —— `(𝒞^istr)^birat` は standard 型

原文 (FrdI p.88):
> Proposition 4.8. (Birationalization of a Frobenioid II)

★`Definition 3.1, (i)` の 5 条のうち、(a)(c)(d)(e) は在庫にある
(`birat_isOfQuasiIsotropicType` / `isOfFrobeniusIsotropicType_of_isotropic` /
`birat_isOfFrobeniusNormalizedType` / `trivialOn_isNonDilatingOn`)。
★★残っていた (b)(Frobenius-compact 対象)が `Prop48Cond3.lean` で埋まったので、
ここで **`Istr` へ持ち上げて**組み立てる。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

section IstrLift

variable {Dd : Type u} [Category.{v} Dd] {Cc : Type u2} [Category.{v2} Cc]
  {Φ₀ : MonoidOn.{v, u, w} Dd} (Q : PreFrobenioid Cc Φ₀) (Fq : FrobenioidCore Q)

/-- ★isotropic な対象を `Istr` の対象として見る。 -/
abbrev istrObj {A : Cc} (hA : IsIsotropic Q A) : Istr Q := ⟨A, hA⟩

/-- ★`Istr` の自己射モノイドは元の自己射モノイドと同型。 -/
def istrEndEquiv {A : Cc} (hA : IsIsotropic Q A) : End (istrObj Q hA) ≃* End A where
  toFun f := f.hom
  invFun g := ⟨g⟩
  left_inv _ := rfl
  right_inv _ := rfl
  map_mul' _ _ := rfl

/-- ★`OTimes` は `Istr` と行き来する。 -/
theorem istr_mem_otimes_iff {A : Cc} (hA : IsIsotropic Q A) (f : End (istrObj Q hA)) :
    f ∈ OTimes (istrPre Q Fq) (istrObj Q hA)
      ↔ (istrEndEquiv Q hA) f ∈ OTimes Q A := by
  constructor
  · rintro ⟨⟨hb, hl⟩, hu⟩
    refine ⟨⟨?_, ?_⟩, hu.map (istrEndEquiv Q hA).toMonoidHom⟩
    · show Q.Base (f.hom) = Q.Base (𝟙 A)
      have h1 : (istrPre Q Fq).Base f = (istrPre Q Fq).Base (𝟙 (istrObj Q hA)) := hb
      rw [istrPre_Base _ f, istrPre_Base _ (𝟙 (istrObj Q hA))] at h1
      exact h1
    · show Q.degFr (f.hom) = 1
      have h2 : (istrPre Q Fq).degFr f = 1 := hl
      rw [istrPre_degFr _ f] at h2
      exact h2
  · rintro ⟨⟨hb, hl⟩, hu⟩
    refine ⟨⟨?_, ?_⟩, hu.map (istrEndEquiv Q hA).symm.toMonoidHom⟩
    · show (istrPre Q Fq).Base f = (istrPre Q Fq).Base (𝟙 (istrObj Q hA))
      rw [istrPre_Base _ f, istrPre_Base _ (𝟙 (istrObj Q hA))]
      exact hb
    · show (istrPre Q Fq).degFr f = 1
      rw [istrPre_degFr _ f]
      exact hl

/-- ★`Istr` の自己同型は元の自己同型と 1 対 1。 -/
def istrIsoEquiv {A : Cc} (hA : IsIsotropic Q A) :
    ((istrObj Q hA) ≅ (istrObj Q hA)) ≃ (A ≅ A) where
  toFun θ := Iso.mk θ.hom.hom θ.inv.hom
    (congrArg InducedCategory.Hom.hom θ.hom_inv_id)
    (congrArg InducedCategory.Hom.hom θ.inv_hom_id)
  invFun θ := Iso.mk ⟨θ.hom⟩ ⟨θ.inv⟩
    (InducedCategory.Hom.ext θ.hom_inv_id)
    (InducedCategory.Hom.ext θ.inv_hom_id)
  left_inv _ := rfl
  right_inv _ := rfl

/-- ★共役は対応と交換する。 -/
theorem istr_endConj {A : Cc} (hA : IsIsotropic Q A)
    (θ : (istrObj Q hA) ≅ (istrObj Q hA)) (u : End (istrObj Q hA)) :
    (istrEndEquiv Q hA) (endConj θ u)
      = endConj (istrIsoEquiv Q hA θ) ((istrEndEquiv Q hA) u) := rfl

/-- ★★★★**`Frobenius-compact` は `Istr` へ持ち上がる**。 -/
theorem istr_frobeniusCompact {A : Cc} (hA : IsIsotropic Q A)
    (h : IsFrobeniusCompact Q A) :
    IsFrobeniusCompact (istrPre Q Fq) (istrObj Q hA) :=
  isFrobeniusCompact_transport Q (istrPre Q Fq)
    (istrEndEquiv Q hA).symm (istrIsoEquiv Q hA).symm
    (fun u => by
      rw [istr_mem_otimes_iff Q Fq hA ((istrEndEquiv Q hA).symm u),
        MulEquiv.apply_symm_apply])
    (fun θ u => by
      refine (istrEndEquiv Q hA).injective ?_
      rw [MulEquiv.apply_symm_apply, istr_endConj, MulEquiv.apply_symm_apply,
        Equiv.apply_symm_apply])
    h

def istr_frobeniusCompact.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 86,
    item := "Remark 4.5.1 — Frobenius-compact 対象は 𝒞^istr へ持ち上がる",
    sectionId := "frdi-rmk-4-5-1" }

end IstrLift

/-! ## ★★★★★★`Proposition 4.8, (iii)` の組み立て -/

section Assemble

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ} {G : Frobenioid P}

/-- ★★★★★★**[FrdI] Proposition 4.8, (iii)** ——
`𝒞` が isotropic 型で、birationally Frobenius-normalized 型で、
`Φ` に 0 でない値が 1 つあるなら、`𝒞^birat` は **standard 型**。

★★★(b) の Frobenius-compact 対象が `Prop48Cond3.lean` で埋まったので組み上がった。

★仮定 `hex`(0 でない `x₀` の存在)は「`Φ` が自明でない」ということであり、
`Φ` が全対象で自明な場合は `𝒞` 自身が group-like なので
`IsOfStandardType` の `groupLikeCompact` が別途効く(★その枝は未)。 -/
theorem prop_4_8_iii
    (hiso : IsOfIsotropicType P)
    (hfn : ∀ X : BiratCat P G, IsFrobeniusNormalized (biratPre P G) X)
    (hdivS : ∀ (Y : C) (a : Φ.val (P.toElem.obj Y).base),
      ∃ u : OTri P Y, P.Div (((u : End Y) : Y ⟶ Y)) = a)
    (hndOn : MonoidOn.IsNonDilatingOn Φ)
    (hFSMFF : IsOfFSMFFType D)
    (hex : ∃ (A₀ : BiratCat P G)
      (x₀ : Φ.val (P.toElem.obj (biratDown P G A₀)).base), x₀ ≠ 0) :
    IsOfStandardType D (BiratCat P G) (biratPre P G) (biratFrobenioid P G hfn).core where
  quasiIsotropic := birat_isOfQuasiIsotropicType P G hfn hiso
  frobIsotropic := isOfFrobeniusIsotropicType_of_isotropic (prop_4_8_i hiso)
  groupLikeCompact := fun _ => by
    obtain ⟨A₀, x₀, hx₀⟩ := hex
    exact ⟨istrObj (biratPre P G) (prop_4_8_i hiso A₀),
      istr_frobeniusCompact (biratPre P G) _ (prop_4_8_i hiso A₀)
        (birat_isFrobeniusCompact_of_ne_zero hdivS hfn hndOn A₀ x₀ hx₀)⟩
  frobNormalized := hfn
  baseFSMFF := hFSMFF
  phiNonDilating := trivialOn_isNonDilatingOn

def prop_4_8_iii.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 88, item := "Proposition 4.8, (iii)",
    sectionId := "frdi-prop-4-8" }

end Assemble

end ABC3.Found.FrdI
