/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Prop44Ker

/-!
# [FrdI] `Proposition 4.4, (ii)` —— `𝒞^birat` の (iii)(d) 圏同値 2 本

## ★★測定の訂正(2026-08-18)—— 21 条だけでは `Frobenioid` にならない

★`Prop44Core.lean` / `Prop44Otri.lean` は `FrobenioidCore (biratPre P G)`(**21 条**)を
扱ってきたが、`Frobenioid` は

  `core` ＋ **`coaPreUnderEquiv`** ＋ **`coaPreOverEquiv`**((iii)(d) の 2 本)

の 3 フィールドである。★★**(iii)(d) の 2 本は未実装だった。**
`ResearchPaper/frdi-decomposition.json` の `otricomm` チェーンにも入っていなかった。

## ★★★しかし `𝒞^birat` では (iii)(d) は易しい

理由は 2 つ:

1. ★**`𝒞^birat` の単系は自明**(`biratPre P G : PreFrobenioid _ (trivialOn D)`)。
   よって `Order(Φ^birat(A))` は**1 対象・1 射**の圏である
   —— 目標側に条件が無いので、**充満性と本質的全射性は「射が在ればよい」**だけになる。
2. ★**co-angular pre-step はすべて同型**(`birat_isIso_of_coaPre_birat`)。
   よってスライス/コスライスの任意の 2 対象の間に射が作れる(`Z⁻¹ ≫ W`)。

★忠実性は一般論(全射性・単射性)から。

## ★これで `Prop 4.4, (ii)` に残るのは `otriBase` 1 条だけになった

`birat_frobenioid_of_comm`(本ファイル)は
**`𝒪^▷(A^birat)` の可換性さえあれば `Frobenioid (biratPre P G)` が出る**ことを与える。
★可換性が言えている場合(birat-Frobenius-normalized 型)の系も置く。
-/

namespace ABC3.Found.FrdI

open CategoryTheory

universe v u w v2 u2

/-! ## ★1. 自明単系の 1 元性 -/

section Trivial

variable {D : Type u} [Category.{v} D]

instance trivialOn_subsingleton (A : D) : Subsingleton ((trivialOn.{v, u, w} D).val A) :=
  inferInstanceAs (Subsingleton PUnit)

instance orderCat_trivialOn_subsingleton (A : D) :
    Subsingleton (OrderCat ((trivialOn.{v, u, w} D).val A)) :=
  inferInstanceAs (Subsingleton PUnit)

end Trivial

section Birat

variable {D : Type u} [Category.{v} D] {C : Type u2} [Category.{v2} C]
  {Φ : MonoidOn.{v, u, w} D} (P : PreFrobenioid C Φ) (G : Frobenioid P)

/-! ## ★2. 忠実性

★`Prop32Frob.lean` の `coaPreUnder_faithful` / `coaPreOver_faithful` と同じ議論だが、
そちらを import すると重いので、ここでは `𝒞^birat` について直接書く。 -/

/-- ★前置は忠実(全射性から)。 -/
theorem birat_coaPreUnder_faithful (A : BiratCat P G) :
    letI := coaPreProp_isMultiplicative (biratPre P G) (birat_coAngularComp P G)
    (coaPreUnderFunctor (biratPre P G) A).Faithful := by
  letI := coaPreProp_isMultiplicative (biratPre P G) (birat_coAngularComp P G)
  refine ⟨fun {Z W} {f g} _ => ?_⟩
  refine Under.UnderMorphism.ext (InducedWideCategory.Hom.ext ?_)
  haveI : Epi Z.hom.hom := (biratPre P G).totEpiC _ _ _
  exact (cancel_epi Z.hom.hom).mp
    ((congrArg InducedWideCategory.Hom.hom (Under.w f)).trans
      (congrArg InducedWideCategory.Hom.hom (Under.w g)).symm)

/-- ★後置は忠実(pre-step が mono から)。 -/
theorem birat_coaPreOver_faithful (A : BiratCat P G) :
    letI := coaPreProp_isMultiplicative (biratPre P G) (birat_coAngularComp P G)
    (coaPreOverFunctor (biratPre P G) A).Faithful := by
  letI := coaPreProp_isMultiplicative (biratPre P G) (birat_coAngularComp P G)
  refine ⟨fun {Z W} {f g} _ => ?_⟩
  refine Over.OverMorphism.ext (InducedWideCategory.Hom.ext ?_)
  haveI : Mono W.hom.hom := birat_preStepMono P G W.hom.hom W.hom.property.2
  exact (cancel_mono W.hom.hom).mp
    ((congrArg InducedWideCategory.Hom.hom (Over.w f)).trans
      (congrArg InducedWideCategory.Hom.hom (Over.w g)).symm)

/-! ## ★3. 充満性 —— co-angular pre-step が同型だから射が作れる -/

/-- ★★前置は充満 —— `Z⁻¹ ≫ W` を取ればよい。 -/
theorem birat_coaPreUnder_full (A : BiratCat P G) :
    letI := coaPreProp_isMultiplicative (biratPre P G) (birat_coAngularComp P G)
    (coaPreUnderFunctor (biratPre P G) A).Full := by
  letI := coaPreProp_isMultiplicative (biratPre P G) (birat_coAngularComp P G)
  refine ⟨fun {Z W} _ => ?_⟩
  haveI hz := birat_isIso_of_coaPre_birat Z.hom.hom Z.hom.property.1 Z.hom.property.2
  refine ⟨Under.homMk (show (⟨Z.right.obj⟩ : WideSubcategory (coaPreProp (biratPre P G)))
      ⟶ ⟨W.right.obj⟩ from
      ⟨inv Z.hom.hom ≫ W.hom.hom,
        birat_coAngularComp P G _ _ (isCoAngular_of_isIso (biratPre P G) (inv Z.hom.hom))
          W.hom.property.1,
        IsPreStep.comp (biratPre P G) (isPreStep_of_isIso (biratPre P G) (inv Z.hom.hom))
          W.hom.property.2⟩) ?_, ?_⟩
  · refine WideSubcategory.hom_ext _ ?_
    show Z.hom.hom ≫ (inv Z.hom.hom ≫ W.hom.hom) = W.hom.hom
    rw [← Category.assoc, IsIso.hom_inv_id, Category.id_comp]
  · exact Subsingleton.elim _ _

/-- ★★後置は充満 —— `Z ≫ W⁻¹` を取ればよい。 -/
theorem birat_coaPreOver_full (A : BiratCat P G) :
    letI := coaPreProp_isMultiplicative (biratPre P G) (birat_coAngularComp P G)
    (coaPreOverFunctor (biratPre P G) A).Full := by
  letI := coaPreProp_isMultiplicative (biratPre P G) (birat_coAngularComp P G)
  refine ⟨fun {Z W} _ => ?_⟩
  haveI hw := birat_isIso_of_coaPre_birat W.hom.hom W.hom.property.1 W.hom.property.2
  refine ⟨Over.homMk (show (⟨Z.left.obj⟩ : WideSubcategory (coaPreProp (biratPre P G)))
      ⟶ ⟨W.left.obj⟩ from
      ⟨Z.hom.hom ≫ inv W.hom.hom,
        birat_coAngularComp P G _ _ Z.hom.property.1
          (isCoAngular_of_isIso (biratPre P G) (inv W.hom.hom)),
        IsPreStep.comp (biratPre P G) Z.hom.property.2
          (isPreStep_of_isIso (biratPre P G) (inv W.hom.hom))⟩) ?_, ?_⟩
  · refine WideSubcategory.hom_ext _ ?_
    show (Z.hom.hom ≫ inv W.hom.hom) ≫ W.hom.hom = Z.hom.hom
    rw [Category.assoc, IsIso.inv_hom_id, Category.comp_id]
  · exact Subsingleton.elim _ _

/-! ## ★4. 本質的全射性 —— 目標が 1 対象なので `𝟙` でよい -/

theorem birat_coaPreUnder_essSurj (A : BiratCat P G) :
    letI := coaPreProp_isMultiplicative (biratPre P G) (birat_coAngularComp P G)
    (coaPreUnderFunctor (biratPre P G) A).EssSurj := by
  letI := coaPreProp_isMultiplicative (biratPre P G) (birat_coAngularComp P G)
  refine ⟨fun _ => ?_⟩
  exact ⟨Under.mk (𝟙 (⟨A⟩ : WideSubcategory (coaPreProp (biratPre P G)))),
    ⟨eqToIso (Subsingleton.elim _ _)⟩⟩

theorem birat_coaPreOver_essSurj (A : BiratCat P G) :
    letI := coaPreProp_isMultiplicative (biratPre P G) (birat_coAngularComp P G)
    (coaPreOverFunctor (biratPre P G) A).EssSurj := by
  letI := coaPreProp_isMultiplicative (biratPre P G) (birat_coAngularComp P G)
  refine ⟨fun _ => ?_⟩
  exact ⟨Over.mk (𝟙 (⟨A⟩ : WideSubcategory (coaPreProp (biratPre P G)))),
    ⟨eqToIso (congrArg Opposite.op (Subsingleton.elim _ _))⟩⟩

/-! ## ★5. (iii)(d) の 2 本 -/

/-- ★★★**(iii)(d) 前置** —— `_A(𝒞^birat,coa-pre) ≃ Order(0)`。 -/
theorem birat_coaPreUnderEquiv (A : BiratCat P G) :
    letI := coaPreProp_isMultiplicative (biratPre P G) (birat_coAngularComp P G)
    (coaPreUnderFunctor (biratPre P G) A).IsEquivalence :=
  ⟨birat_coaPreUnder_faithful P G A, birat_coaPreUnder_full P G A,
    birat_coaPreUnder_essSurj P G A⟩

/-- ★★★**(iii)(d) 後置** —— `(𝒞^birat,coa-pre)_A ≃ Order(0)^opp`。 -/
theorem birat_coaPreOverEquiv (A : BiratCat P G) :
    letI := coaPreProp_isMultiplicative (biratPre P G) (birat_coAngularComp P G)
    (coaPreOverFunctor (biratPre P G) A).IsEquivalence :=
  ⟨birat_coaPreOver_faithful P G A, birat_coaPreOver_full P G A,
    birat_coaPreOver_essSurj P G A⟩

/-! ## ★6. `Frobenioid (biratPre P G)` —— 残るのは可換性 1 条だけ -/

/-- ★★★★**`𝒪^▷(A^birat)` が可換なら `𝒞^birat` は Frobenioid**。

★★これで `Proposition 4.4, (ii)` に残るのは**可換性のみ**になった
(21 条は `birat_frobenioidCore_of_comm`、(iii)(d) は本ファイル)。 -/
theorem birat_frobenioid_of_comm
    (hcomm : ∀ (X : BiratCat P G) (x y : OTri (biratPre P G) X), x * y = y * x) :
    Frobenioid (biratPre P G) :=
  { core := birat_frobenioidCore_of_comm P G hcomm
    coaPreUnderEquiv := birat_coaPreUnderEquiv P G
    coaPreOverEquiv := birat_coaPreOverEquiv P G }

/-- ★★★**`𝒞` が birationally Frobenius-normalized 型なら `𝒞^birat` は Frobenioid**。 -/
theorem birat_frobenioid_of_frobNormalized
    (hfn : ∀ X : BiratCat P G, IsFrobeniusNormalized (biratPre P G) X) :
    Frobenioid (biratPre P G) :=
  { core := birat_frobenioidCore_of_frobNormalized P G hfn
    coaPreUnderEquiv := birat_coaPreUnderEquiv P G
    coaPreOverEquiv := birat_coaPreOverEquiv P G }

end Birat

/-! ## ★出典の紐付け(条つき) -/

def birat_coaPreUnderEquiv.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 85, item := "Proposition 4.4, (ii) — (iii)(d) の 2 本",
    sectionId := "frdi-prop-4-4" }

end ABC3.Found.FrdI
