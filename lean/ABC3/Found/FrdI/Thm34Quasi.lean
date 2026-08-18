/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Thm34Pre
import ABC3.Found.FrdI.Prop33UnTr

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

end ABC3.Found.FrdI
