/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop44Gl
import ABC3.Found.FrdI.Prop22Star
import ABC3.Found.FrdI.Def45
import ABC3.Found.FrdI.Remark311
import ABC3.Found.FrdI.Prop44Phi.Proposition44

/-!
# Prop44Phi —— `[FrdI] Proposition 4.8` の分

☆もとの 1 枚を**条なしの項目ごと**に割ったものである（第 1458、案 a）。
★「1 ファイル = 1 ノード」を回復するための分割で、中身は動かしていない。
-/

namespace ABC3.Found.FrdI

universe v u v2 u2 w
open CategoryTheory Opposite
variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ) (G : Frobenioid P)

/-! ## ★[FrdI] Proposition 4.8 —— 平方化 II

原文 (FrdI p.88):
> Proposition 4.8. (Birationalization of a Frobenioid II)

原文 (FrdI p.88):
> (i) If C is of isotropic type, then so is Cbirat.

★原典の証明は「Assertion (i) follows formally from Proposition 4.4, (iv)」であり、
★(iv) の辞書 `birat_isIsotropic_iff` がそのまま使える。 -/

variable {P G} in
/-- ★★★★**[FrdI] Proposition 4.8, (i)** —— isotropic 型は birat で保たれる。 -/
theorem prop_4_8_i (h : IsOfIsotropicType P) : IsOfIsotropicType (biratPre P G) :=
  fun X => (birat_isIsotropic_iff P G X).mpr (h X)

def prop_4_8_i.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 88, item := "Proposition 4.8, (i)",
    sectionId := "frdi-prop-4-8" }

/-! ## ★`Proposition 4.8, (iii)` の (d)(e) —— 底の圧と `Φ`

★(d) `𝒟` の FSMFF 性は **`𝒞^birat` でも底の圧が変わらない**のでそのまま。
★(e) `𝒞^birat` の `Φ` は `trivialOn D`(全値 `PUnit`)なので自明。 -/

/-- ★★**`trivialOn D` は non-dilating** —— `Proposition 4.8, (iii)` の (e)。

★値がすべて `PUnit` なので、`MChar` も subsingleton である。 -/
theorem trivialOn_isNonDilatingOn :
    MonoidOn.IsNonDilatingOn (trivialOn.{v, u, w} D) := by
  intro A e _
  haveI : Subsingleton (MChar ((trivialOn.{v, u, w} D).val A)) := ⟨fun a b => by
    induction a using AddCon.induction_on with | _ x =>
    induction b using AddCon.induction_on with | _ y =>
    exact congrArg _ (Subsingleton.elim (α := PUnit.{w + 1}) x y)⟩
  ext x
  exact Subsingleton.elim _ _

def trivialOn_isNonDilatingOn.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 88,
    item := "Proposition 4.8, (iii) — Φ が non-dilating",
    sectionId := "frdi-prop-4-8" }

/-! ## ★`Proposition 4.8, (iii)` の (a) と (c) -/

/-- ★★**(c)** birationally Frobenius-normalized 型 ⇒ `𝒞^birat` は Frobenius-normalized 型。

★`IsBirationallyFrobeniusNormalized` はそのまま `biratPre` 上の
`IsFrobeniusNormalized` として定義されており、
`BiratCat P G := C` なので**対象は同じ**。★したがって定義の展開だけで出る。 -/
theorem birat_isOfFrobeniusNormalizedType
    (h : IsOfBirationallyFrobeniusNormalizedType C P G) :
    IsOfFrobeniusNormalizedType (biratPre P G) := h

/-- ★★**(a) の後半** isotropic 型 ⇒ Frobenius-isotropic 型。

★`𝟙 A` を Frobenius 型射に取ればよい。 -/
theorem isOfFrobeniusIsotropicType_of_isotropic {C' : Type u2} [Category.{v2} C']
    {P' : PreFrobenioid C' Φ} (h : IsOfIsotropicType P') :
    IsOfFrobeniusIsotropicType P' :=
  fun A => ⟨A, 𝟙 A, ⟨⟨isCoAngular_id P' A, by
      show P'.Div (𝟙 A) = 0
      exact P'.Div_id A⟩, by
      show IsIso (P'.Base (𝟙 A))
      rw [P'.Base_id]; infer_instance⟩, h A⟩

def birat_isOfFrobeniusNormalizedType.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 88,
    item := "Proposition 4.8, (iii) — Frobenius-normalized 型",
    sectionId := "frdi-prop-4-8" }

/-- ★★**(a) の前半** —— `𝒞` が isotropic 型なら `𝒞^birat` は quasi-isotropic 型。

★原文の `Remark 3.1.1`(実装済 `isOfQuasiIsotropicType_of_isOfIsotropicType`)を
`𝒞^birat` に当てるだけ。 -/
theorem birat_isOfQuasiIsotropicType
    (hfn : ∀ X : BiratCat P G, IsFrobeniusNormalized (biratPre P G) X)
    (h : IsOfIsotropicType P) :
    IsOfQuasiIsotropicType (BiratCat P G) (biratPre P G) :=
  isOfQuasiIsotropicType_of_isOfIsotropicType (biratPre P G)
    (biratFrobenioid P G hfn).core (prop_4_8_i h)

def birat_isOfQuasiIsotropicType.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 88,
    item := "Proposition 4.8, (iii) — quasi-isotropic 型",
    sectionId := "frdi-prop-4-8" }

end ABC3.Found.FrdI
