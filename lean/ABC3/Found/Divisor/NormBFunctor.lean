/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Divisor.SchemeFFMap
import ABC3.Found.Divisor.NormB

/-!
# `B(L) → B(M)` の関手性(鎖 `normalize` の `B-functor` の最後の 1 段)

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.110。

原文 (FrdI p.110):
> for the group of rational functions on Vi[Li] whose zeroes and poles are supported

## ★★中身は 3 本

1. `x ∉ D_M` ⟹ `g(x) ∉ D_L`(`normMap ≫ normDown = normDown` から直ちに)
2. `ord_{g(x)}(u) = 0` ⟹ `u` は `𝒪_{V[L],g(x)}` の単元(`SchemeOrdUnit.lean`)
3. 単元は単元へ移る ⟹ `ord_x(g^*u) = 0`(`SchemeFFMap.lean` の四角形)

## ★★★残る 1 つの仮定 —— 余次元が上がらないこと

2 を使うには `g(x)` が**余次元 1 か生成点**でなければならない。
★これは「有限支配射は余次元を保つ」という一般論(整拡大 + 生成点の分離性 ⟹
going-down、したがって高さの保存)であり、**mathlib には高さの保存が無い**
(2026-08-21 実測: `Algebra.HasGoingDown` は在るが `height` の保存は 0 件)。
★★そこで本ファイルは**その言明を仮定として受け取る**形で書き、
節点 `codim-preserved` に切り出す。

## ★本ファイルで閉じること

| 定理 | 中身 |
|---|---|
| `normFF_eq_ffMap` | `normFF` は `ffMap` の特別な場合 |
| `not_mem_DLSet_of_not_mem` | ★`x ∉ D_M` ⟹ `g(x) ∉ D_L` |
| `exists_unit_of_eq_genericPoint` | 生成点の茎は関数体そのもの |
| `normFFUnits_mem_BSubgroup` | ★★★★**`B(L) → B(M)`**(余次元の仮定つき) |
-/

namespace ABC3.Found.Divisor

open AlgebraicGeometry CategoryTheory ABC3.Found.FrdI

universe u

/-- ★生成点では茎が関数体そのものなので、単元は茎から来る。 -/
theorem exists_unit_of_eq_genericPoint {X : Scheme.{u}} [IsIntegral X] {z : X}
    (hz : z = genericPoint (X : Type u)) (u : (X.functionField)ˣ) :
    ∃ t : (X.presheaf.stalk z)ˣ,
      algebraMap (X.presheaf.stalk z) X.functionField (t : X.presheaf.stalk z) = (u : _) := by
  subst hz
  refine ⟨u, ?_⟩
  show ((X.presheaf.stalkSpecializes _).hom) (u : X.functionField) = (u : X.functionField)
  rw [TopCat.Presheaf.stalkSpecializes_refl]
  rfl

/-- ★`WithBot ℕ∞` では `≤ 1` かつ `≠ 1` なら `≤ 0`(離散性)。 -/
theorem withBot_le_zero_of_le_one_of_ne {d : WithBot ℕ∞} (h : d ≤ 1) (hne : d ≠ 1) : d ≤ 0 := by
  induction d using WithBot.recBotCoe with
  | bot => exact bot_le
  | coe n =>
    induction n using ENat.recTopCoe with
    | top =>
      have htop : (1 : ℕ∞) < ⊤ := by simp
      have h1 : ((⊤ : ℕ∞) : WithBot ℕ∞) ≤ ((1 : ℕ∞) : WithBot ℕ∞) := h
      rw [WithBot.coe_le_coe] at h1
      exact absurd h1 (not_le_of_gt htop)
    | coe m =>
      have h1 : m ≤ 1 := by exact_mod_cast h
      have h2 : m ≠ 1 := fun hm => hne (by rw [hm]; rfl)
      have hm0 : m = 0 := by omega
      subst hm0
      rfl

/-- ★★**茎の次元が `0` なら茎は体** —— したがって関数体の元はすべて茎から来る。 -/
theorem exists_unit_of_ringKrullDim_le_zero {X : Scheme.{u}} [IsIntegral X] {z : X}
    (hz : ringKrullDim (X.presheaf.stalk z) ≤ 0) (u : (X.functionField)ˣ) :
    ∃ t : (X.presheaf.stalk z)ˣ,
      algebraMap (X.presheaf.stalk z) X.functionField (t : X.presheaf.stalk z)
        = (u : X.functionField) := by
  haveI : Ring.KrullDimLE 0 (X.presheaf.stalk z) := Ring.krullDimLE_iff.mpr (by simpa using hz)
  have hfield : IsField (X.presheaf.stalk z) := Ring.KrullDimLE.isField_of_isReduced
  have hsurj : Function.Surjective
      (algebraMap (X.presheaf.stalk z) X.functionField) :=
    IsFractionRing.surjective_iff_isField.mpr hfield
  obtain ⟨t0, ht0⟩ := hsurj (u : X.functionField)
  have ht0ne : t0 ≠ 0 := by
    intro h
    rw [h, map_zero] at ht0
    exact u.ne_zero ht0.symm
  obtain ⟨w, hw⟩ := hfield.mul_inv_cancel ht0ne
  have hunit : IsUnit t0 := ⟨⟨t0, w, hw, by rw [mul_comm]; exact hw⟩, rfl⟩
  exact ⟨hunit.unit, by rw [hunit.unit_spec]; exact ht0⟩

variable (V : Scheme.{u}) [IsIntegral V] {Kbar : Type u} [Field Kbar]
  [Algebra V.functionField Kbar]

/-- ★`normFF` は一般の `ffMap` の特別な場合。 -/
theorem normFF_eq_ffMap {L M : FinSub V.functionField Kbar} (f : L ⟶ M) :
    normFF V f = ffMap (normMap V f) := rfl

/-- ★★**`x ∉ D_M` なら `g(x) ∉ D_L`** —— `normMap ≫ normDown = normDown` から。 -/
theorem not_mem_DLSet_of_not_mem (DK : Set (PrimeDivisorPt V))
    {L M : FinSub V.functionField Kbar} (f : L ⟶ M)
    (x : PrimeDivisorPt (normObj V M)) (hx : x ∉ DLSet V DK M)
    (z : PrimeDivisorPt (normObj V L)) (hz : z.1 = (normMap V f).base x.1) :
    z ∉ DLSet V DK L := by
  rintro ⟨y, hy, hxy⟩
  refine hx ⟨y, hy, ?_⟩
  show (normDown V M).base x.1 ∈ closure ({y.1} : Set V)
  have hcomp : (normDown V M).base x.1
      = (normDown V L).base ((normMap V f).base x.1) :=
    congrArg (fun h : normObj V M ⟶ V => h.base x.1) (normMap_normDown V f).symm
  rw [hcomp, ← hz]
  exact hxy

/-- ★★★★**`B(L) → B(M)`** —— 零点・極が `D_L` にある有理函数は、
引き戻しても零点・極が `D_M` にある。

★仮定 `hdim` は「支配射で余次元が上がらない」ことで、節点 `codim-preserved` である。
★★`= 1` なら余次元 1 の点として `B(L)` の条件がそのまま効き、
`< 1` なら茎が**体**になるので `u` は自動的に茎の単元である。 -/
theorem normFFUnits_mem_BSubgroup (DK : Set (PrimeDivisorPt V))
    {L M : FinSub V.functionField Kbar} (f : L ⟶ M)
    [IsLocallyNoetherian (normObj V L)] [IsLocallyNoetherian (normObj V M)]
    (hdim : ∀ x : PrimeDivisorPt (normObj V M),
      ringKrullDim ((normObj V L).presheaf.stalk ((normMap V f).base x.1)) ≤ 1)
    (u : ((normObj V L).functionField)ˣ)
    (hu : u ∈ BSubgroup V DK L (normObj_isNormalScheme V L)) :
    normFFUnits V f u ∈ BSubgroup V DK M (normObj_isNormalScheme V M) := by
  intro x hx
  show ordPt (normObj V M) (normObj_isNormalScheme V M) x
      (normFF V f (u : (normObj V L).functionField)) = 0
  obtain ⟨t, ht⟩ : ∃ t : ((normObj V L).presheaf.stalk ((normMap V f).base x.1))ˣ,
      algebraMap ((normObj V L).presheaf.stalk ((normMap V f).base x.1))
        (normObj V L).functionField
        (t : (normObj V L).presheaf.stalk ((normMap V f).base x.1))
        = (u : (normObj V L).functionField) := by
    by_cases hcod : IsCodimOnePt (normObj V L) ((normMap V f).base x.1)
    · have hznot : (⟨(normMap V f).base x.1, hcod⟩ : PrimeDivisorPt (normObj V L))
          ∉ DLSet V DK L :=
        not_mem_DLSet_of_not_mem V DK f x hx _ rfl
      exact exists_unit_of_ordPt_eq_zero (normObj_isNormalScheme V L) _ u.ne_zero
        (hu _ hznot)
    · exact exists_unit_of_ringKrullDim_le_zero
        (withBot_le_zero_of_le_one_of_ne (hdim x) hcod) u
  rw [normFF_eq_ffMap, ← ht]
  exact ordPt_ffMap_eq_zero_of_isUnit (normMap V f) (normObj_isNormalScheme V M) x t

/-! ### ★出典の紐付け -/

/-- ★★★★locator —— `Example 6.1` の `B(L)` の関手性。 -/
def normFFUnits_mem_BSubgroup.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 110,
    item := "Example 6.1 — B(L) の関手性(L → M で B(L) → B(M))",
    sectionId := "frdi-example-6-1" }

end ABC3.Found.Divisor
