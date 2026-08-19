/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop48Std
import ABC3.Found.FrdI.Prop41Cat

/-!
# [FrdI] Proposition 4.8, (iii) —— `Φ` が自明な枝(＝ `𝒞` が group-like)

★`Prop48Std.lean` の `prop_4_8_iii` は「`Φ` に 0 でない値がある」ことを仮定した。
★★その補集合、すなわち **`Φ` が全対象で自明**な場合は、
`Definition 1.3` の意味で **step が 1 本も無い**ので
`𝒞 → 𝒞^birat` が**忠実充満**になり、`𝒞` 側の Frobenius-compact 対象がそのまま移る。

## ★手筋 —— `idxBiratOne` は添字圏の**始対象**

`IdxBirat P G A = (Over (coaPreObj P G A))ᵒᵖ` であり、
`Over.mk (𝟙 _)` は `Over` の**終対象**(`Over.mkIdTerminal`)。
★したがって `idxBiratOne` は `IdxBirat` の**始対象**であり、
- **全射性**: 任意の代表元 `(Z, φ)` は `Z.hom` が同型なので `𝟙` の側へ引き戻せる
- **単射性**: 終対象への 2 本の射は等しいので、`eq_iff` の `u` と `v` が一致する
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w u2 v2

section NoStep

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} {P : PreFrobenioid C Φ} {G : Frobenioid P}

/-- ★★**`Φ` が自明なら co-angular pre-step はすべて同型**。 -/
theorem isIso_of_coaPre_of_phiTrivial
    (hiso : ∀ X : C, IsIsotropic P X)
    (htriv : ∀ (Y : C) (x : Φ.val (P.toElem.obj Y).base), x = 0)
    {X Y : C} (f : X ⟶ Y) (hf : coaPreProp P f) : IsIso f := by
  by_contra hni
  have hstep : IsStep P f := ⟨hf.2, hni⟩
  exact ((isStep_iff_preStepVal_ne_zero hiso f hf.2).mp hstep) (htriv _ _)

variable (P G) in
/-- ★★`idxBiratOne` から任意の添字への射(終対象性の反対)。 -/
noncomputable def idxOneTo (A : C) (Z : IdxBirat P G A) : idxBiratOne P G A ⟶ Z :=
  Quiver.Hom.op ((Over.mkIdTerminal (X := coaPreObj P G A)).from Z.unop)

theorem idxOneTo_left (A : C) (Z : IdxBirat P G A) :
    (idxOneTo P G A Z).unop.left = Z.unop.hom := by
  have h := Over.w ((Over.mkIdTerminal (X := coaPreObj P G A)).from Z.unop)
  show ((Over.mkIdTerminal (X := coaPreObj P G A)).from Z.unop).left = Z.unop.hom
  rw [← h]
  exact (Category.comp_id _).symm

/-- ★★★**`Φ` が自明なら `toHomBirat` は全射**。 -/
theorem toHomBirat_surjective
    (hno : ∀ {X Y : C} (f : X ⟶ Y), coaPreProp P f → IsIso f)
    (A B : C) : Function.Surjective (toHomBirat (P := P) (G := G) (A := A) (B := B)) := by
  intro z
  obtain ⟨Z, φ, rfl⟩ := HomBirat.exists_rep z
  set u : idxBiratOne P G A ⟶ Z := idxOneTo P G A Z with hu
  have hleft : u.unop.left = Z.unop.hom := idxOneTo_left (P := P) (G := G) A Z
  haveI : IsIso (u.unop.left.hom) := by
    refine hno _ ?_
    have hh : u.unop.left.hom = Z.unop.hom.hom :=
      congrArg InducedWideCategory.Hom.hom hleft
    rw [hh]
    exact Z.unop.hom.property
  refine ⟨inv (u.unop.left.hom) ≫ φ, ?_⟩
  show HomBirat.mk (idxBiratOne P G A) (inv (u.unop.left.hom) ≫ φ) = HomBirat.mk Z φ
  rw [← HomBirat.mk_map u (inv (u.unop.left.hom) ≫ φ), ← Category.assoc,
    IsIso.hom_inv_id, Category.id_comp]

/-- ★★★**`Φ` が自明なら `toHomBirat` は単射**。 -/
theorem toHomBirat_injective
    (hno : ∀ {X Y : C} (f : X ⟶ Y), coaPreProp P f → IsIso f)
    (A B : C) : Function.Injective (toHomBirat (P := P) (G := G) (A := A) (B := B)) := by
  intro φ ψ hpq
  obtain ⟨V, u, v, huv⟩ := HomBirat.eq_iff.mp hpq
  have huv0 : u = v := by
    have h : u.unop = v.unop :=
      (Over.mkIdTerminal (X := coaPreObj P G A)).hom_ext _ _
    exact Quiver.Hom.unop_inj h
  rw [huv0] at huv
  haveI : IsIso (v.unop.left.hom) := hno _ v.unop.left.property
  exact (cancel_epi (v.unop.left.hom)).mp huv

/-- ★★★★**`Φ` が自明なら `𝒞 → 𝒞^birat` は忠実**。 -/
theorem toBiratCat_faithful_of_noStep
    (hno : ∀ {X Y : C} (f : X ⟶ Y), coaPreProp P f → IsIso f) :
    (toBiratCat P G).Faithful where
  map_injective {A B} := toHomBirat_injective hno A B

/-- ★★★★**`Φ` が自明なら `𝒞 → 𝒞^birat` は充満**。 -/
theorem toBiratCat_full_of_noStep
    (hno : ∀ {X Y : C} (f : X ⟶ Y), coaPreProp P f → IsIso f) :
    (toBiratCat P G).Full where
  map_surjective {A B} h := toHomBirat_surjective hno A B h

def toBiratCat_faithful_of_noStep.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 82,
    item := "Proposition 4.4, (i) — step が無ければ 𝒞 → 𝒞^birat は忠実充満",
    sectionId := "frdi-prop-4-4" }

end NoStep

end ABC3.Found.FrdI
