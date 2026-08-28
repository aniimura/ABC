/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.RatioChartHom
import Mathlib.AlgebraicGeometry.ValuativeCriterion
import Mathlib.AlgebraicGeometry.ProjectiveSpectrum.Proper
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★★★★分離的なら生成点で決まる（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.6。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

## ★★★★★★★★★★★★★★★★★★これは何か —— 最後の紐を結ぶ道具

`§9-943`（有限素点では座標の 1 つが他を割る）と `§9-944`（比の組から `ℙᴺ` の点を作る）で

    `Spec (𝓞_F)_Q ⟶ D₊(x_j) ⊆ ℙᴺ`

が**構成できる**ようになった。★残るのは「これが局所化した点と一致する」ことである。

★★★★**本ファイルがその道具を与える**:

    `R` が付値環、`Z` が分離的なら、`Spec R ⟶ Z` は**生成点で決まる**

——分離性の付値判定法（`IsSeparated.valuativeCriterion`、mathlib）そのものである。

## ★★★機構 —— 分離性の付値判定法

`Scheme.IsSeparated Z` は `terminal.from Z` が分離的であることであり、
mathlib の `IsSeparated.valuativeCriterion` はそれから
**付値可換四角の持ち上げが高々 1 つ**（`ValuativeCriterion.Uniqueness`）を出す。
★`α`・`β` をその 2 つの持ち上げと見れば `Subsingleton` から `α = β` が出る。

## ★★これを `(𝓞_F)_Q` に当てる配線

* `(𝓞_F)_Q` は離散付値環（`§9-943` と同じ理由）
* `F` は `(𝓞_F)_Q` の商体
  （`IsFractionRing.isFractionRing_of_isDomain_of_isLocalization`）
* `ℙᴺ` は分離的（mathlib の `Scheme.IsSeparated (Proj 𝒜)`）

★★★これで `§9-940` の有限素点の整合を閉じる道具が揃った。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory Limits NumberField

/-! ## ★★★★★★★★★★★★★★★★★★付値環の点は生成点で決まる -/

/-- ★★★★★★★★★★★★★★★★★★**分離的なスキームへの射は、
付値環のスペクトルからなら生成点で決まる**。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

★機構は mathlib の `IsSeparated.valuativeCriterion`
（分離的 ⟹ 付値可換四角の持ち上げは高々 1 つ）だけである。
★★`α`・`β` を 2 つの持ち上げと見て `Subsingleton` を当てる。 -/
theorem eq_of_generic_eq (R K : Type) [CommRing R] [IsDomain R] [ValuationRing R]
    [Field K] [Algebra R K] [IsFractionRing R K]
    {Z : Scheme.{0}} [Z.IsSeparated]
    (α β : Spec (CommRingCat.of R) ⟶ Z)
    (h : Spec.map (CommRingCat.ofHom (algebraMap R K)) ≫ α
       = Spec.map (CommRingCat.ofHom (algebraMap R K)) ≫ β) : α = β := by
  have hsep : ValuativeCriterion.Uniqueness (terminal.from Z) :=
    IsSeparated.valuativeCriterion _
  let S : ValuativeCommSq (terminal.from Z) :=
    { R := R, K := K,
      i₁ := Spec.map (CommRingCat.ofHom (algebraMap R K)) ≫ α,
      i₂ := terminal.from _,
      commSq := ⟨terminal.hom_ext _ _⟩ }
  have hsub := hsep S
  exact congrArg CommSq.LiftStruct.l (Subsingleton.elim
    (⟨α, rfl, terminal.hom_ext _ _⟩ : S.commSq.LiftStruct)
    (⟨β, h.symm, terminal.hom_ext _ _⟩ : S.commSq.LiftStruct))

/-! ## ★★★`(𝓞_F)_Q` から `F` への配線 -/

/-- ★**`(𝓞_F)_Q` は `F` の部分環である**（`Q.primeCompl ≤ nonZeroDivisors`）。 -/
noncomputable instance atPrimeAlgebra (F : Type) [Field F] [NumberField F]
    (Q : Ideal (𝓞 F)) [Q.IsPrime] : Algebra (Localization.AtPrime Q) F :=
  IsLocalization.localizationAlgebraOfSubmonoidLe (Localization.AtPrime Q) F
    Q.primeCompl (nonZeroDivisors (𝓞 F)) Q.primeCompl_le_nonZeroDivisors

instance atPrimeIsScalarTower (F : Type) [Field F] [NumberField F]
    (Q : Ideal (𝓞 F)) [Q.IsPrime] : IsScalarTower (𝓞 F) (Localization.AtPrime Q) F :=
  IsLocalization.localization_isScalarTower_of_submonoid_le _ _ _ _
    Q.primeCompl_le_nonZeroDivisors

/-- ★★**`F` は `(𝓞_F)_Q` の商体である**。 -/
instance atPrimeIsFractionRing (F : Type) [Field F] [NumberField F]
    (Q : Ideal (𝓞 F)) [Q.IsPrime] : IsFractionRing (Localization.AtPrime Q) F :=
  IsFractionRing.isFractionRing_of_isDomain_of_isLocalization Q.primeCompl _ _

/-! ## ★★★★★★★★★★★★★★★★★★★数体の局所化での形 -/

/-- ★★★★★★★★★★★★★★★★★★★**`Spec (𝓞_F)_Q ⟶ Z`（`Z` 分離的）は
生成点で決まる**。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

★★★`Z = ℙᴺ` に当てるのが目的である（mathlib に `Scheme.IsSeparated (Proj 𝒜)` がある）。
★これで `§9-944` が構成した点と「局所化した点」を比べる道具が揃った。 -/
theorem eq_of_generic_eq_atPrime (F : Type) [Field F] [NumberField F]
    (Q : Ideal (𝓞 F)) (hQ : Q.IsMaximal)
    {Z : Scheme.{0}} [Z.IsSeparated]
    (α β : haveI := hQ.isPrime; Spec (CommRingCat.of (Localization.AtPrime Q)) ⟶ Z)
    (h : haveI := hQ.isPrime;
      Spec.map (CommRingCat.ofHom (algebraMap (Localization.AtPrime Q) F)) ≫ α
        = Spec.map (CommRingCat.ofHom (algebraMap (Localization.AtPrime Q) F)) ≫ β) :
    α = β := by
  haveI := hQ.isPrime
  haveI hQ0 : Q ≠ ⊥ := Ring.ne_bot_of_isMaximal_of_not_isField hQ (RingOfIntegers.not_isField F)
  haveI : IsDiscreteValuationRing (Localization.AtPrime Q) :=
    IsLocalization.AtPrime.isDiscreteValuationRing_of_dedekind_domain (𝓞 F) hQ0 _
  exact eq_of_generic_eq (Localization.AtPrime Q) F α β h

/-! ## ★出典の紐付け(`.src`) -/

def eq_of_generic_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(分離的なスキームへの射は付値環のスペクトルからなら生成点で決まる)",
    sectionId := "genell-prop-1-4" }

def atPrimeAlgebra.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)((𝓞_F)_Q は F の部分環である)",
    sectionId := "genell-prop-1-4" }

def atPrimeIsFractionRing.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(F は (𝓞_F)_Q の商体である)",
    sectionId := "genell-prop-1-4" }

def eq_of_generic_eq_atPrime.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(Spec (𝓞_F)_Q ⟶ Z(Z 分離的)は生成点で決まる)",
    sectionId := "genell-prop-1-4" }

def eq_of_generic_eq_atPrime.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]" "IsSeparated.valuativeCriterion(分離性の付値判定法)"
      (.inMathlib "AlgebraicGeometry.IsSeparated.valuativeCriterion") 3,
    .citation "[mathlib]" "Scheme.IsSeparated (Proj 𝒜)(ℙᴺ は分離的)"
      (.inMathlib "AlgebraicGeometry.Proj.isSeparated") 2,
    .citation "[mathlib]"
      "IsFractionRing.isFractionRing_of_isDomain_of_isLocalization(F は (𝓞_F)_Q の商体)"
      (.inMathlib "IsFractionRing.isFractionRing_of_isDomain_of_isLocalization") 2,
    .implicitStep
      ("★★★★測定(2026-08-29): §9-943(座標の 1 つが他を割る)と §9-944(比の組から点を作る)で " ++
       "Spec (𝓞_F)_Q ⟶ D₊(x_j) ⊆ ℙᴺ が**構成できる**ようになった。" ++
       "残るのは「これが局所化した点と一致する」ことであり、" ++
       "★分離性の付値判定法がその道具である" ++
       "——付値環のスペクトルからの射は生成点で決まる") 5,
    .implicitStep
      ("★配線は 3 つ: (𝓞_F)_Q は離散付値環、F はその商体、ℙᴺ は分離的。" ++
       "どれも mathlib にある") 3 ]

end ABC3.Found.GenEll
