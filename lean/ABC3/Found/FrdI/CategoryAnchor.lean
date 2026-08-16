import ABC3.Found.FrdI.CategoryVocabulary
import Mathlib.CategoryTheory.Comma.Over.Basic
import Mathlib.Algebra.Group.Subgroup.Basic
import Mathlib.CategoryTheory.Endomorphism
import Mathlib.Data.Set.Finite.Basic

/-!
# [FrdI] §0 の圏の語彙(第 2 部)—— categorical quotient / anchor / iso-subanchor

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.18。

★**`Definition 3.1, (i)`(quasi-isotropic 型)へ向けた語彙**である。
`CategoryVocabulary.lean` の `IsTotallyEpimorphic` / `IsIrreducibleMor` の上に積む。

## ★この節で写すもの

| 原文 | ここでの名前 |
|---|---|
| `φ : A → B` が `A` の `G` による categorical quotient | `IsCategoricalQuotient` |
| その `mono-minimal` 性 | `IsMonoMinimalQuotient` |
| `anchor` | `IsAnchor` |
| `subanchor` | `IsSubanchor` |
| `iso-subanchor` | `IsIsoSubanchor` |

原文 (FrdI p.18):
> subgroup, then we shall say that an arrow

原文 (FrdI p.18):
> of A by G if the following conditions hold: (a)

## ★合成の向きについて

★**原文の `φ ◦ γ` は Lean では `γ ≫ φ`** である。
`categorical quotient` の条件 (a)「`φ ◦ γ = φ` for all `γ ∈ G`」は
`γ ≫ φ = φ` と書く。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u

variable {C : Type u} [Category.{v} C]

/-! ## ★categorical quotient

原文 (FrdI p.18):
> is mono-minimal if the following condition
-/

/-- **[FrdI] §0** —— `φ : A ⟶ B` が `A` の部分群 `G ⊆ Aut(A)` による
**categorical quotient** であること。

★(a) `G` の作用で不変、★(b) その普遍性。 -/
def IsCategoricalQuotient {A B : C} (G : Subgroup (Aut A)) (φ : A ⟶ B) : Prop :=
  (∀ γ ∈ G, (γ : Aut A).hom ≫ φ = φ) ∧
    ∀ (Cc : C) (ψ : A ⟶ Cc), (∀ γ ∈ G, (γ : Aut A).hom ≫ ψ = ψ) →
      ∃! ψ' : B ⟶ Cc, ψ = φ ≫ ψ'

/-- **[FrdI] §0** —— categorical quotient `φ` が **mono-minimal** であること。

★「`φ = φ' ∘ ζ`(`ζ : A ⟶ A'` が mono で、`G ≅ G' ⊆ Aut(A')` が `ζ` について
両者の作用と両立する)ならば `ζ` は同型」。 -/
def IsMonoMinimalQuotient {A B : C} (G : Subgroup (Aut A)) (φ : A ⟶ B) : Prop :=
  ∀ (A' : C) (ζ : A ⟶ A') (φ' : A' ⟶ B), Mono ζ → φ = ζ ≫ φ' →
    ∀ (G' : Subgroup (Aut A')) (e : G ≃* G'),
      (∀ γ : G, ((γ : Aut A).hom : A ⟶ A) ≫ ζ = ζ ≫ ((e γ : Aut A').hom : A' ⟶ A')) →
      IsIso ζ

/-- ★**非退化(上)** —— 同型はつねに自明群による mono-minimal categorical quotient。

原文 (FrdI p.18):
> follows that an isomorphism is always a mono-minimal categorical quotient of its

★**原文が「Thus, [by total epimorphicity] it follows that」と述べる事実**である。 -/
theorem isCategoricalQuotient_of_isIso {A B : C} (φ : A ⟶ B) [IsIso φ] :
    IsCategoricalQuotient (⊥ : Subgroup (Aut A)) φ := by
  constructor
  · intro γ hγ
    have hγ1 : γ = 1 := hγ
    rw [hγ1]
    exact Category.id_comp φ
  · intro Cc ψ _
    exact ⟨inv φ ≫ ψ, by simp, fun y hy => by
      rw [hy, ← Category.assoc, IsIso.inv_hom_id, Category.id_comp]⟩

/-- ★**非退化(上)** —— 同型は mono-minimal でもある。

★`ζ ≫ φ' = φ` で `φ` が同型なら `ζ` は split mono、mono と合わせて同型。 -/
theorem isMonoMinimalQuotient_of_isIso (hte : IsTotallyEpimorphic C) {A B : C}
    (φ : A ⟶ B) [IsIso φ] :
    IsMonoMinimalQuotient (⊥ : Subgroup (Aut A)) φ := by
  intro A' ζ φ' hmono hfac _ _ _
  haveI : Epi ζ := hte _ _ ζ
  refine ⟨φ' ≫ inv φ, ?_, ?_⟩
  · rw [← Category.assoc, ← hfac, IsIso.hom_inv_id]
  · refine (cancel_epi ζ).mp ?_
    rw [← Category.assoc, ← Category.assoc, ← hfac, IsIso.hom_inv_id,
      Category.id_comp, Category.comp_id]

/-! ## ★anchor / subanchor / iso-subanchor

原文 (FrdI p.18):
> exist finitely many isomorphism classes of objects of AC that arise from irreducible
-/

variable (C) in
/-- **[FrdI] §0** —— `A` が **anchor** であること。

★「`A` から出る **irreducible** な射がなす `_A𝒞` の対象の同型類が**有限個**」。
★有限集合 `s` を取り、どの irreducible 射も `s` の元と同型になる、と書く。 -/
def IsAnchor (A : C) : Prop :=
  ∃ s : Set (Under A), s.Finite ∧
    ∀ (B : C) (φ : A ⟶ B), IsIrreducibleMor φ → ∃ X ∈ s, Nonempty (Under.mk φ ≅ X)

variable (C) in
/-- **[FrdI] §0** —— `A` が **subanchor** であること(anchor への射がある)。 -/
def IsSubanchor (A : C) : Prop := ∃ (B : C) (_ : A ⟶ B), IsAnchor C B

variable (C) in
/-- **[FrdI] §0** —— `A` が **iso-subanchor** であること。

★「subanchor `B`、部分群 `G ⊆ Aut(B)`、および `B` の `G` による
**mono-minimal categorical quotient** `B ⟶ A` が存在する」。 -/
def IsIsoSubanchor (A : C) : Prop :=
  ∃ (B : C) (_ : IsSubanchor C B) (G : Subgroup (Aut B)) (φ : B ⟶ A),
    IsCategoricalQuotient G φ ∧ IsMonoMinimalQuotient G φ

/-- ★**非退化** —— subanchor はつねに iso-subanchor(恒等射を取ればよい)。 -/
theorem isIsoSubanchor_of_isSubanchor (hte : IsTotallyEpimorphic C) {A : C}
    (h : IsSubanchor C A) : IsIsoSubanchor C A :=
  ⟨A, h, ⊥, 𝟙 A, isCategoricalQuotient_of_isIso (𝟙 A),
    isMonoMinimalQuotient_of_isIso hte (𝟙 A)⟩

/-! ## ★★★出典の紐付け(`.src`) -/

/-- **[FrdI] §0** —— categorical quotient / mono-minimal。 -/
def IsCategoricalQuotient.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 18, item := "§0 Categories — categorical quotient",
    sectionId := "frdi-s0-cat-quotient" }

/-- **[FrdI] §0** —— anchor / subanchor / iso-subanchor。 -/
def IsAnchor.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 18, item := "§0 Categories — anchor / iso-subanchor",
    sectionId := "frdi-s0-anchor" }

end ABC3.Found.FrdI
