/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop110
import ABC3.Found.ProL.ZetaPow

/-!
# [FrdI] Definition 2.8 —— unit-profinite 型・pro-`l` 部分・ζ 乗写像

原文 (FrdI p.52):
> We shall refer to the factor M [l] as the pro-l portion of M

## ★★★逸脱の記録(2026-08-19)

原文 (i) の角括弧「**uniquely determined**」——
「位相的有限生成な副有限群の位相は群構造だけで決まる」——は
**Nikolov–Segal の定理**(有限指数部分群は開)に依る深い事実で、mathlib に無い。
★★したがって (i) は**「そのような位相が存在する」**として形式化する。
★定義として下流(`Proposition 2.9` / `5.6` / `5.7`)が要求するのは
「位相を 1 つ選べること」なので、**消費側に影響しない逸脱**である。

## ★もう 1 つの逸脱(安全側)

原文 (ii)(iii) は `M` に**位相的有限生成**を仮定するが、
★分解 `M ≅ ∏_l M[l]`(`decompEquiv`)も ζ 乗写像の全単射性(`zetaPow_bijective`)も
**可換な副有限群なら成り立つ**。原文より**弱い仮定**で述べているので、
原文の主張はこれの系である。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

section UnitProfinite

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ)

/-- ★`𝒪^×(A)` の元の逆射も `𝒪^×(A)` に入る。 -/
theorem otimes_inv_mem {A : C} {δ : End A} (h : δ ∈ OTimes P A)
    (hi : IsIso ((δ : A ⟶ A))) :
    ((@inv _ _ _ _ ((δ : A ⟶ A)) hi) : End A) ∈ OTimes P A := by
  haveI := hi
  have hb : P.Base ((δ : A ⟶ A)) = 𝟙 _ := by
    rw [show P.Base ((δ : A ⟶ A)) = P.Base (𝟙 A) from h.1.1, P.Base_id]
  refine ⟨⟨?_, degFr_of_isIso P _⟩, ?_⟩
  · show P.Base (inv ((δ : A ⟶ A))) = P.Base (𝟙 A)
    have h1 : P.Base ((δ : A ⟶ A)) ≫ P.Base (inv ((δ : A ⟶ A))) = P.Base (𝟙 A) := by
      rw [← P.Base_comp, IsIso.hom_inv_id]
    rw [hb, Category.id_comp] at h1
    exact h1
  · exact (CategoryTheory.isUnit_iff_isIso _).mpr inferInstance

/-- ★★**`𝒪^×(A)` は群**。 -/
noncomputable instance otimesGroup (A : C) : Group (OTimes P A) where
  inv x :=
    haveI : IsIso (((x : End A) : A ⟶ A)) := (CategoryTheory.isUnit_iff_isIso _).mp x.2.2
    ⟨((inv (((x : End A)) : A ⟶ A)) : End A), otimes_inv_mem P x.2 inferInstance⟩
  inv_mul_cancel x := by
    haveI : IsIso (((x : End A) : A ⟶ A)) := (CategoryTheory.isUnit_iff_isIso _).mp x.2.2
    refine Subtype.ext ?_
    show ((x : End A) : A ⟶ A) ≫ (inv (((x : End A)) : A ⟶ A)) = 𝟙 A
    exact IsIso.hom_inv_id _

/-- ★★★**[FrdI] Definition 2.8, (i)** —— `𝒞` が **unit-profinite 型**。

★★角括弧「uniquely determined」は **Nikolov–Segal** に依るのでここでは主張しない
(ファイル冒頭の逸脱の記録を見よ)。 -/
def IsOfUnitProfiniteType : Prop :=
  ∀ A : C, ∃ t : TopologicalSpace (OTimes P A),
    (∀ a b : OTimes P A, a * b = b * a) ∧
    (letI := t; IsTopologicalGroup (OTimes P A)) ∧
    (letI := t; CompactSpace (OTimes P A)) ∧
    (letI := t; TotallyDisconnectedSpace (OTimes P A)) ∧
    ∃ S : Finset (OTimes P A),
      (letI := t; closure ((Subgroup.closure (S : Set (OTimes P A)) : Subgroup _) :
        Set (OTimes P A)) = Set.univ)

def IsOfUnitProfiniteType.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 52, item := "Definition 2.8, (i) — unit-profinite 型",
    sectionId := "frdi-def-2-8" }

end UnitProfinite

/-- ★★★★★★**[FrdI] Definition 2.8** —— 条なしの locator。

| 条 | 実装 |
|---|---|
| (i) unit-profinite 型 | `IsOfUnitProfiniteType` |
| (ii) pro-`l` 部分 `M[l]` と分解 | `ABC3.Found.ProL.lPart` / `lPartGrp` / `isProL_lPartGrp` / ★`decompEquiv` |
| (iii) ζ 乗写像と co-prime 型 | `ABC3.Found.ProL.zetaPow` / `IsCoPrimeType` / ★`zetaPow_bijective` |

★★逸脱 2 件(ファイル冒頭に記録):
(a) (i) の「uniquely determined」は主張しない(Nikolov–Segal)、
(b) (ii)(iii) は**位相的有限生成を仮定しない**(原文より弱い仮定＝安全側)。 -/
def def_2_8.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 52, item := "Definition 2.8",
    sectionId := "frdi-def-2-8" }

end ABC3.Found.FrdI
