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
  have hh : (idxOneTo P G A Z).unop.left.hom = Z.unop.hom.hom :=
    congrArg InducedWideCategory.Hom.hom (idxOneTo_left (P := P) (G := G) A Z)
  haveI hIso : IsIso ((idxOneTo P G A Z).unop.left.hom) := by
    refine hno _ ?_
    rw [hh]
    exact Z.unop.hom.property
  refine ⟨(@inv _ _ _ _ ((idxOneTo P G A Z).unop.left.hom) hIso) ≫ φ, ?_⟩
  show HomBirat.mk (idxBiratOne P G A)
    ((@inv _ _ _ _ ((idxOneTo P G A Z).unop.left.hom) hIso) ≫ φ) = HomBirat.mk Z φ
  rw [← HomBirat.mk_map (idxOneTo P G A Z)
      ((@inv _ _ _ _ ((idxOneTo P G A Z).unop.left.hom) hIso) ≫ φ),
    ← Category.assoc, IsIso.hom_inv_id, Category.id_comp]

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

/-! ## ★★★★★★`Φ` 自明の枝での Frobenius-compact 対象 -/

/-- ★★`Φ` が自明なら `𝒞` は group-like 型。 -/
theorem isOfGroupLikeType_of_phiTrivial
    (htriv : ∀ (Y : C) (x : Φ.val (P.toElem.obj Y).base), x = 0) :
    IsOfGroupLikeType P := by
  intro A a
  obtain ⟨y, rfl⟩ := toChar_surjective _ a
  rw [htriv A y]
  exact map_zero _

/-- ★★★★★**`Φ` が自明なら `Frobenius-compact` は `𝒞^birat` へ移る**。

★`toBiratCat` が忠実充満なので、`End`・`OTimes`・自己同型の 3 つがそのまま対応する。 -/
theorem birat_isFrobeniusCompact_of_noStep
    (hno : ∀ {X Y : C} (f : X ⟶ Y), coaPreProp P f → IsIso f)
    (A : C) (h : IsFrobeniusCompact P A) :
    IsFrobeniusCompact (biratPre P G) ((toBiratCat P G).obj A) := by
  haveI := toBiratCat_faithful_of_noStep hno (G := G)
  haveI := toBiratCat_full_of_noStep hno (G := G)
  refine isFrobeniusCompact_transport P (biratPre P G)
    (MulEquiv.ofBijective (Functor.mapEnd A (toBiratCat P G))
      ⟨fun x y hxy => (toBiratCat P G).map_injective hxy,
       fun y => (toBiratCat P G).map_surjective y⟩)
    ((Functor.FullyFaithful.ofFullyFaithful (toBiratCat P G)).isoEquiv) ?_ ?_ h
  · intro u
    constructor
    · rintro ⟨⟨hb, hl⟩, hu⟩
      refine ⟨⟨?_, ?_⟩, ?_⟩
      · show (biratPre P G).Base ((toBiratCat P G).map (u : A ⟶ A))
          = (biratPre P G).Base (𝟙 ((toBiratCat P G).obj A))
        rw [biratPre_Base, biratPre_Base]
        show biratBase (toHomBirat (P := P) (G := G) (u : A ⟶ A))
          = biratBase (toHomBirat (P := P) (G := G) (𝟙 A))
        rw [biratBase_toHomBirat, biratBase_toHomBirat]
        exact hb
      · show (biratPre P G).degFr ((toBiratCat P G).map (u : A ⟶ A)) = 1
        rw [biratPre_degFr]
        show biratDeg (toHomBirat (P := P) (G := G) (u : A ⟶ A)) = 1
        rw [biratDeg_toHomBirat]
        exact hl
      · haveI : IsIso ((u : A ⟶ A)) := (CategoryTheory.isUnit_iff_isIso u).mp hu
        refine (CategoryTheory.isUnit_iff_isIso _).mpr ?_
        show IsIso ((toBiratCat P G).map (u : A ⟶ A))
        infer_instance
    · rintro ⟨⟨hb, hl⟩, hu⟩
      refine ⟨⟨?_, ?_⟩, ?_⟩
      · show P.Base (u : A ⟶ A) = P.Base (𝟙 A)
        have hb' : biratBase (toHomBirat (P := P) (G := G) (u : A ⟶ A))
            = biratBase (toHomBirat (P := P) (G := G) (𝟙 A)) := hb
        rw [biratBase_toHomBirat, biratBase_toHomBirat] at hb'
        exact hb'
      · show P.degFr (u : A ⟶ A) = 1
        have hl' : biratDeg (toHomBirat (P := P) (G := G) (u : A ⟶ A)) = 1 := hl
        rw [biratDeg_toHomBirat] at hl'
        exact hl'
      · haveI : IsIso ((toBiratCat P G).map (u : A ⟶ A)) := by
          have h1 : IsIso ((Functor.mapEnd A (toBiratCat P G)) u) :=
            (CategoryTheory.isUnit_iff_isIso _).mp hu
          exact h1
        refine (CategoryTheory.isUnit_iff_isIso u).mpr ?_
        exact isIso_of_reflects_iso (u : A ⟶ A) (toBiratCat P G)
  · intro θ u
    show (toBiratCat P G).map ((endConj θ u : End A) : A ⟶ A)
      = endConj _ ((toBiratCat P G).map ((u : End A) : A ⟶ A))
    show (toBiratCat P G).map (θ.inv ≫ (u : A ⟶ A) ≫ θ.hom)
      = ((toBiratCat P G).mapIso θ).inv ≫ (toBiratCat P G).map ((u : A ⟶ A))
        ≫ ((toBiratCat P G).mapIso θ).hom
    rw [CategoryTheory.Functor.map_comp, CategoryTheory.Functor.map_comp]
    rfl

def birat_isFrobeniusCompact_of_noStep.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 88,
    item := "Proposition 4.8, (iii) — Φ が自明な枝の Frobenius-compact 対象",
    sectionId := "frdi-prop-4-8" }

/-! ## ★★★★★★`Proposition 4.8, (iii)` —— 場合分けなしの形 -/

/-- ★★★★★★**[FrdI] Proposition 4.8, (iii)**(完全形) ——
`𝒞` が isotropic 型・birationally Frobenius-normalized 型・standard 型なら
`𝒞^birat` は **standard 型**。

★★`Φ` に 0 でない値があるかどうかで場合分けする:
- ある場合: `birat_isFrobeniusCompact_of_ne_zero`(`Prop48Cond3.lean`)
- 無い場合: `𝒞` は group-like なので `standard` の `groupLikeCompact` が
  compact 対象を与え、**step が無いので `𝒞 → 𝒞^birat` は忠実充満**、
  そのまま移送できる。 -/
theorem prop_4_8_iii_full
    (hiso : IsOfIsotropicType P)
    (hfn : ∀ X : BiratCat P G, IsFrobeniusNormalized (biratPre P G) X)
    (hdivS : ∀ (Y : C) (a : Φ.val (P.toElem.obj Y).base),
      ∃ u : OTri P Y, P.Div (((u : End Y) : Y ⟶ Y)) = a)
    (hndOn : MonoidOn.IsNonDilatingOn Φ)
    (hFSMFF : IsOfFSMFFType D)
    (hstd : IsOfStandardType D C P G.core) :
    IsOfStandardType D (BiratCat P G) (biratPre P G) (biratFrobenioid P G hfn).core := by
  classical
  by_cases hex : ∃ (A₀ : BiratCat P G)
      (x₀ : Φ.val (P.toElem.obj (biratDown P G A₀)).base), x₀ ≠ 0
  · exact prop_4_8_iii hiso hfn hdivS hndOn hFSMFF hex
  · -- ★`Φ` はすべての対象で自明
    have htriv : ∀ (Y : C) (x : Φ.val (P.toElem.obj Y).base), x = 0 := by
      intro Y x
      by_contra hx
      exact hex ⟨Y, x, hx⟩
    have hno : ∀ {X Y : C} (f : X ⟶ Y), coaPreProp P f → IsIso f :=
      fun f hf => isIso_of_coaPre_of_phiTrivial hiso htriv f hf
    obtain ⟨A, hA⟩ := hstd.groupLikeCompact (isOfGroupLikeType_of_phiTrivial htriv)
    refine
      { quasiIsotropic := birat_isOfQuasiIsotropicType P G hfn hiso
        frobIsotropic := isOfFrobeniusIsotropicType_of_isotropic (prop_4_8_i hiso)
        groupLikeCompact := fun _ => ?_
        frobNormalized := hfn
        baseFSMFF := hFSMFF
        phiNonDilating := trivialOn_isNonDilatingOn }
    refine ⟨istrObj (biratPre P G) (prop_4_8_i hiso A.obj), ?_⟩
    refine istr_frobeniusCompact (biratPre P G) _ (prop_4_8_i hiso A.obj) ?_
    exact birat_isFrobeniusCompact_of_noStep hno A.obj
      (istr_frobeniusCompact' P G.core A.property hA)

def prop_4_8_iii_full.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 88, item := "Proposition 4.8, (iii)",
    sectionId := "frdi-prop-4-8" }

/-- ★★★★★★**[FrdI] Proposition 4.8** —— (i)〜(iv) がすべて揃った。 -/
def prop_4_8.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 88, item := "Proposition 4.8",
    sectionId := "frdi-prop-4-8" }

end NoStep

end ABC3.Found.FrdI
