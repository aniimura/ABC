/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import Mathlib.RingTheory.Ideal.GoingUp
import Mathlib.RingTheory.LocalRing.MaximalIdeal.Basic
import Mathlib.RingTheory.Jacobson.Ideal
import ABC3.Meta.Claim

/-!
# 第 1365 ブロック —— **整拡大では `m·S` は Jacobson 根基に入る**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★これは何か——第 1364 が要る `I ≤ jacobson ⊥`

第 1364（`isLocalRing_of_isAdicComplete_of_domain`）は `I ≤ jacobson ⊥` を要求する。
★`I = m_R · S`（`S` は局所環 `R` の整拡大）ならこれは**自動**である
——`S` の極大イデアル `n` の縮約は整拡大なので極大（`isMaximal_comap_of_isIntegral_of_isMaximal`、
在庫）、`R` は局所なのでそれは `m_R` に他ならない。
-/

namespace ABC3.Found.GenEll

open ABC3.Meta

variable {R S : Type*} [CommRing R] [IsLocalRing R] [CommRing S] [Algebra R S]
  [Algebra.IsIntegral R S]

/-- ★★★★★★★★
**整拡大では `m_R · S` は Jacobson 根基に入る**——★**無条件**（第 1365）。 -/
theorem map_maximalIdeal_le_jacobson :
    (IsLocalRing.maximalIdeal R).map (algebraMap R S) ≤ Ideal.jacobson (⊥ : Ideal S) := by
  refine le_sInf ?_
  rintro n ⟨-, hn⟩
  haveI := hn
  have hcomap : (n.comap (algebraMap R S)).IsMaximal :=
    Ideal.isMaximal_comap_of_isIntegral_of_isMaximal n
  have heq : IsLocalRing.maximalIdeal R = n.comap (algebraMap R S) :=
    (IsLocalRing.eq_maximalIdeal hcomap).symm
  rw [heq]
  exact Ideal.map_comap_le

/-- ★★★★★★**整拡大では `m_R · S ≠ ⊤`**——★**無条件**（第 1365）。 -/
theorem map_maximalIdeal_ne_top [Nontrivial S] :
    (IsLocalRing.maximalIdeal R).map (algebraMap R S) ≠ ⊤ := by
  obtain ⟨n, hn, -⟩ := Ideal.exists_le_maximal (⊥ : Ideal S) bot_ne_top
  intro htop
  refine hn.ne_top (top_le_iff.mp ?_)
  have h1 : (IsLocalRing.maximalIdeal R).map (algebraMap R S) ≤ n := by
    refine le_trans map_maximalIdeal_le_jacobson ?_
    exact sInf_le ⟨bot_le, hn⟩
  rw [htop] at h1
  exact h1

/-! ## ★出典の紐付け(`.src`) -/

def map_maximalIdeal_le_jacobson.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(整拡大では m_R·S は Jacobson 根基に入る。★無条件)",
    sectionId := "genell-thm-3-8" }

def map_maximalIdeal_ne_top.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(整拡大では m_R·S ≠ ⊤。★無条件)",
    sectionId := "genell-thm-3-8" }

def map_maximalIdeal_le_jacobson.needs : List ProofObligation :=
  [ .implicitStep
      ("★★★★**2026-09-02（第 1365）**——第 1364 が要る `I ≤ jacobson ⊥` と " ++
       "`Nontrivial (A ⧸ I)` をまとめて与える段である。") 19 ]

end ABC3.Found.GenEll
