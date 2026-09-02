/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import Mathlib.RingTheory.DedekindDomain.IntegralClosure
import Mathlib.RingTheory.DiscreteValuationRing.Basic
import Mathlib.RingTheory.Finiteness.Quotient
import Mathlib.LinearAlgebra.FreeModule.Finite.Basic
import ABC3.Found.GenEll.AdicComparable
import ABC3.Found.GenEll.AdicTransfer
import ABC3.Found.GenEll.ArtinianLocal
import ABC3.Found.GenEll.DvrFromDedekind
import ABC3.Found.GenEll.IdempotentLift
import ABC3.Found.GenEll.IntegralJacobson
import ABC3.Found.GenEll.LocalFromQuotient
import ABC3.Meta.Claim

/-!
# 第 1366 ブロック —— **完備 DVR の有限分離整閉包は完備 DVR**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★★★★★★★★★★★★★これは何か——第 1357-1365 の 9 段の**組み立て**

第 1352-1356 で「α 側（`exists_h2_h1_unipotent`）の欠けている原始命題」を測り、
**「完備なイデアルに沿う冪等元の持ち上げ」** に落ちることを突き止めた。
第 1357-1365 でその 9 段をすべて `Found` に積んだ:

| 段 | 内容 | 第 |
|---|---|---|
| 1 | 完備なイデアルに沿う冪等元の持ち上げ | 1357 |
| 2 | 完備な整域の商に非自明な冪等元なし | 1358 |
| 3 | `I`-進完備性は線型同型で移る | 1359 |
| 4 | 有限直積・係数環の移送 | 1360-1361 |
| 5 | 剰余体上有限なら Artin 環 | 1363 |
| 6 | Artin ＋ 冪等元なし ⟹ 局所 | 1362 |
| 7 | Dedekind ＋ 局所 ⟹ DVR | 1363 |
| 8 | 商が局所なら元も局所 | 1364 |
| 9 | 整拡大では `m·S ≤ jacobson ⊥` | 1365 |

★本ブロックはこれらを繋いで **`C = 整閉包(A, L)` は DVR** を出す。

☆道:

1. `C` は `A` 上**有限自由**（`IsIntegralClosure.finite` ＋ `module_free`、在庫）
2. ゆえに `C ≃ₗ[A] (ι → A)`、`A` が `m_A`-進完備なら `C` も（第 1359-1360）
3. `A`-加群としての完備性を `C`-環としての `m_A C`-進完備性に上げる（第 1361）
4. `C ⧸ m_A C` は剰余体 `A ⧸ m_A` 上有限（在庫）ゆえ Artin 環（第 1363）
5. `m_A C ≤ jacobson ⊥`（第 1365）・`m_A C ≠ ⊤`（第 1365）
6. よって `C` は局所（第 1364）
7. `C` は Dedekind（`IsIntegralClosure.isDedekindDomain`、在庫）で体でない
8. よって `C` は DVR（第 1363）
-/

namespace ABC3.Found.GenEll

open ABC3.Meta IsLocalRing

attribute [local instance] Ideal.Quotient.field

section IntegralClosure

set_option linter.unusedSectionVars false

variable {A : Type*} [CommRing A] [IsDomain A] [IsDiscreteValuationRing A]
variable (K : Type*) [Field K] [Algebra A K] [IsFractionRing A K]
variable (L : Type*) [Field L] [Algebra K L] [Algebra A L] [IsScalarTower A K L]
variable [FiniteDimensional K L] [Algebra.IsSeparable K L]
variable (C : Type*) [CommRing C] [IsDomain C] [Algebra A C] [Algebra C L]
variable [IsIntegralClosure C A L] [IsScalarTower A C L]

include A K L

/-- ★★★★★★**`A → L` は単射**（第 1366）。 -/
theorem algebraMap_injective_of_fractionRing :
    Function.Injective (algebraMap A L) := by
  rw [IsScalarTower.algebraMap_eq A K L]
  exact (algebraMap K L).injective.comp (IsFractionRing.injective A K)

/-- ★★★★★★**`L` は `A`-捻れ無し**（第 1366）。 -/
theorem isTorsionFree_of_fractionRing : Module.IsTorsionFree A L :=
  (Module.isTorsionFree_iff_algebraMap_injective (R := A) (A := L)).2
    (algebraMap_injective_of_fractionRing K L)

/-- ★★★★★★**`A → C` は単射**（第 1366）。 -/
theorem algebraMap_injective_integralClosure :
    Function.Injective (algebraMap A C) := by
  have h : Function.Injective ((algebraMap C L).comp (algebraMap A C)) := by
    rw [← IsScalarTower.algebraMap_eq A C L]
    exact algebraMap_injective_of_fractionRing K L
  exact Function.Injective.of_comp h

/-- ★★★★★★★★**`C` は `A` 上有限自由**（第 1366、道の 1）。 -/
theorem module_free_integralClosure : Module.Free A C := by
  haveI := isTorsionFree_of_fractionRing (A := A) K L
  exact IsIntegralClosure.module_free A K L C

/-- ★★★★★★★★★★★★
**完備 DVR `A` の有限分離整閉包 `C` は `m_A C`-進完備**（第 1366、道の 2-3）。 -/
theorem isAdicComplete_map_maximalIdeal [IsAdicComplete (maximalIdeal A) A] :
    IsAdicComplete ((maximalIdeal A).map (algebraMap A C)) C := by
  classical
  haveI : Module.Finite A C := IsIntegralClosure.finite A K L C
  haveI : Module.Free A C := module_free_integralClosure K L C
  haveI : Fintype (Module.Free.ChooseBasisIndex A C) := inferInstance
  have b : Module.Basis (Module.Free.ChooseBasisIndex A C) A C :=
    Module.Free.chooseBasis A C
  haveI : IsAdicComplete (maximalIdeal A) (Module.Free.ChooseBasisIndex A C → A) :=
    isAdicComplete_pi (maximalIdeal A)
  haveI : IsAdicComplete (maximalIdeal A) C :=
    isAdicComplete_of_linearEquiv (maximalIdeal A) b.equivFun.symm
  exact isAdicComplete_map (R := A) (A := C) (maximalIdeal A)

/-- ★★★★★★★★**`C ⧸ m_A C` は Artin 環**（第 1366、道の 4）。 -/
theorem isArtinianRing_quotient_map_maximalIdeal :
    IsArtinianRing (C ⧸ (maximalIdeal A).map (algebraMap A C)) := by
  haveI : IsNoetherian A C := IsIntegralClosure.isNoetherian A K L C
  haveI : IsNoetherian (A ⧸ maximalIdeal A)
      (C ⧸ (maximalIdeal A).map (algebraMap A C)) := inferInstance
  haveI : Module.Finite (A ⧸ maximalIdeal A)
      (C ⧸ (maximalIdeal A).map (algebraMap A C)) := Module.IsNoetherian.finite _ _
  exact isArtinianRing_of_module_finite (A ⧸ maximalIdeal A)

/-- ★★★★★★★★★★★★★★★★
**完備 DVR `A` の有限分離整閉包 `C` は局所環**（第 1366、道の 5-6）。 -/
theorem isLocalRing_integralClosure [IsAdicComplete (maximalIdeal A) A] : IsLocalRing C := by
  haveI : Algebra.IsIntegral A C := IsIntegralClosure.isIntegral_algebra A L
  haveI : IsAdicComplete ((maximalIdeal A).map (algebraMap A C)) C :=
    isAdicComplete_map_maximalIdeal K L C
  haveI : IsArtinianRing (C ⧸ (maximalIdeal A).map (algebraMap A C)) :=
    isArtinianRing_quotient_map_maximalIdeal K L C
  haveI : Nontrivial (C ⧸ (maximalIdeal A).map (algebraMap A C)) :=
    Ideal.Quotient.nontrivial_iff.mpr (map_maximalIdeal_ne_top (R := A) (S := C))
  exact isLocalRing_of_isAdicComplete_of_domain
    ((maximalIdeal A).map (algebraMap A C)) (map_maximalIdeal_le_jacobson (R := A) (S := C))

/-- ★★★★★★★★**`C` は体でない**（第 1366、道の 7）。 -/
theorem not_isField_integralClosure : ¬ IsField C := by
  haveI : Algebra.IsIntegral A C := IsIntegralClosure.isIntegral_algebra A L
  intro hfield
  obtain ⟨t, htmem, htne⟩ : ∃ t ∈ maximalIdeal A, t ≠ 0 :=
    Submodule.ne_bot_iff _ |>.mp (IsDiscreteValuationRing.not_a_field A)
  have hne : algebraMap A C t ≠ 0 := fun h =>
    htne (algebraMap_injective_integralClosure K L C (by simpa using h))
  have hunit : IsUnit (algebraMap A C t) := by
    obtain ⟨y, hy⟩ := (hfield.mul_inv_cancel hne)
    exact isUnit_iff_exists_inv.mpr ⟨y, hy⟩
  exact map_maximalIdeal_ne_top (R := A) (S := C)
    (Ideal.eq_top_of_isUnit_mem _ (Ideal.mem_map_of_mem _ htmem) hunit)

/-- ★★★★★★★★★★★★★★★★★★★★
**完備 DVR の有限分離整閉包は DVR**——★**無条件**（第 1366）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★★★これが第 1352-1356 で測った α 側の欠落を埋める**要石**である。
☆分岐した局所体の上でも「整数環は DVR」が言えるので、
`ζ_l ∈ L_p` の判定と分裂乗法的還元の議論がそのまま回る。 -/
theorem isDiscreteValuationRing_integralClosure [IsAdicComplete (maximalIdeal A) A] :
    IsDiscreteValuationRing C := by
  haveI : IsLocalRing C := isLocalRing_integralClosure (A := A) K L C
  haveI : IsDedekindDomain C := IsIntegralClosure.isDedekindDomain A K L C
  exact isDiscreteValuationRing_of_isDedekindDomain C
    (not_isField_integralClosure (A := A) K L C)

/-! ## ★★★★★★★★★★★★第 1368 —— 下流が使える形にする

★第 1320（`exists_h2_h1_of_bad_prime`）が要求するのは

* `IsDiscreteValuationRing C`（第 1366）
* `IsFractionRing C L`
* `IsAdicComplete (maximalIdeal C) C`

の 3 つである。☆残る 2 つをここで作る。 -/

/-- ★★★★★★**`m_A C ≠ ⊥`**（第 1368）。 -/
theorem map_maximalIdeal_ne_bot :
    (maximalIdeal A).map (algebraMap A C) ≠ ⊥ := by
  obtain ⟨t, htmem, htne⟩ : ∃ t ∈ maximalIdeal A, t ≠ 0 :=
    Submodule.ne_bot_iff _ |>.mp (IsDiscreteValuationRing.not_a_field A)
  intro hbot
  have hmem : algebraMap A C t ∈ (maximalIdeal A).map (algebraMap A C) :=
    Ideal.mem_map_of_mem _ htmem
  rw [hbot, Ideal.mem_bot] at hmem
  exact htne (algebraMap_injective_integralClosure K L C (by simpa using hmem))

/-- ★★★★★★★★★★★★★★★★★★★★
**完備 DVR の有限分離整閉包は `m_C`-進完備**——★**無条件**（第 1368）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★★★これで第 1320（`exists_h2_h1_of_bad_prime`）の
`[IsAdicComplete (IsLocalRing.maximalIdeal R) R]` が**拡大先でも満たせる**。 -/
theorem isAdicComplete_maximalIdeal_integralClosure [IsAdicComplete (maximalIdeal A) A]
    [IsLocalRing C] : IsAdicComplete (maximalIdeal C) C := by
  haveI : Algebra.IsIntegral A C := IsIntegralClosure.isIntegral_algebra A L
  haveI : IsAdicComplete ((maximalIdeal A).map (algebraMap A C)) C :=
    isAdicComplete_map_maximalIdeal K L C
  haveI : IsDiscreteValuationRing C := isDiscreteValuationRing_integralClosure (A := A) K L C
  exact isAdicComplete_maximalIdeal_of_isAdicComplete
    ((maximalIdeal A).map (algebraMap A C))
    (map_maximalIdeal_ne_bot K L C) (map_maximalIdeal_ne_top (R := A) (S := C))

/-- ★★★★★★★★**`C` の商体は `L`**（第 1368）。 -/
theorem isFractionRing_integralClosure : IsFractionRing C L :=
  IsIntegralClosure.isFractionRing_of_finite_extension A K L C

end IntegralClosure

/-! ## ★出典の紐付け(`.src`) -/

def isDiscreteValuationRing_integralClosure.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(完備 DVR の有限分離整閉包は DVR。★無条件)",
    sectionId := "genell-thm-3-8" }

def isLocalRing_integralClosure.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(完備 DVR の有限分離整閉包は局所環。★無条件)",
    sectionId := "genell-thm-3-8" }

def isAdicComplete_map_maximalIdeal.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(完備 DVR の有限分離整閉包は m_A C-進完備。★無条件)",
    sectionId := "genell-thm-3-8" }

def isAdicComplete_maximalIdeal_integralClosure.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(完備 DVR の有限分離整閉包は m_C-進完備。★無条件)",
    sectionId := "genell-thm-3-8" }

def isFractionRing_integralClosure.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(整閉包の商体はもとの体。★無条件)",
    sectionId := "genell-thm-3-8" }

def isDiscreteValuationRing_integralClosure.needs : List ProofObligation :=
  [ .citation "[ABC3]" "isAdicComplete_of_linearEquiv(第 1359、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.isAdicComplete_of_linearEquiv") 1,
    .citation "[ABC3]" "isAdicComplete_pi(第 1360、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.isAdicComplete_pi") 1,
    .citation "[ABC3]" "isAdicComplete_map(第 1361、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.isAdicComplete_map") 1,
    .citation "[ABC3]" "isLocalRing_of_isAdicComplete_of_domain(第 1364、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.isLocalRing_of_isAdicComplete_of_domain") 1,
    .citation "[ABC3]" "map_maximalIdeal_le_jacobson(第 1365、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.map_maximalIdeal_le_jacobson") 1,
    .citation "[ABC3]" "isDiscreteValuationRing_of_isDedekindDomain(第 1363、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.isDiscreteValuationRing_of_isDedekindDomain") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1366）**——第 1357-1365 の **9 段の組み立て**である。" ++
       "☆第 1352-1356 で測った α 側（`exists_h2_h1_unipotent`）の欠落" ++
       "「完備 DVR の有限整閉包は完備 DVR」がこれで埋まった。") 19 ]

end ABC3.Found.GenEll
