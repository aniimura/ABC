/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.ArtinianLocal
import ABC3.Found.GenEll.IdempotentLift
import ABC3.Meta.Claim

/-!
# 第 1364 ブロック —— **商が局所なら元も局所**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★★★★★これは何か——7 段を 1 本にまとめる

第 1357-1363 で 7 段が揃った。★本ブロックはそれらを繋いで

> `A` が整域で `I`-進完備、`I ≤ jacobson ⊥`、`A ⧸ I` が Artin 環なら `A` は局所

という**使う形**にする。

☆道は 2 段:

1. `A ⧸ I` は非自明な冪等元を持たない（第 1358）ので Artin より局所（第 1362）
2. `I ≤ jacobson ⊥` なので `A` の極大イデアルは `A ⧸ I` の極大イデアルと 1 対 1（★本ブロック）
-/

namespace ABC3.Found.GenEll

open ABC3.Meta

variable {A : Type*} [CommRing A]

/-- ★★★★★★★★
**商が局所で `I` が Jacobson 根基に入っていれば元も局所**——★**無条件**（第 1364）。 -/
theorem isLocalRing_of_quotient (I : Ideal A) (hI : I ≤ Ideal.jacobson (⊥ : Ideal A))
    [IsLocalRing (A ⧸ I)] : IsLocalRing A := by
  classical
  set f : A →+* A ⧸ I := Ideal.Quotient.mk I with hf
  have hsurj : Function.Surjective f := Ideal.Quotient.mk_surjective
  set m : Ideal A := Ideal.comap f (IsLocalRing.maximalIdeal (A ⧸ I)) with hm
  have hmmax : m.IsMaximal :=
    Ideal.comap_isMaximal_of_surjective f hsurj (K := IsLocalRing.maximalIdeal (A ⧸ I))
  refine IsLocalRing.of_unique_max_ideal ⟨m, hmmax, fun n hn => ?_⟩
  have hIn : I ≤ n := le_trans hI (sInf_le ⟨bot_le, hn⟩)
  have hcb : Ideal.comap f (⊥ : Ideal (A ⧸ I)) ≤ n := by
    intro x hx
    exact hIn (Ideal.Quotient.eq_zero_iff_mem.mp (by simpa using hx))
  have hcm : Ideal.comap f (Ideal.map f n) = n := by
    rw [Ideal.comap_map_of_surjective f hsurj n, sup_eq_left.mpr hcb]
  have hne : Ideal.map f n ≠ ⊤ := by
    intro htop
    rw [htop, Ideal.comap_top] at hcm
    exact hn.ne_top hcm.symm
  have hmax' : (Ideal.map f n).IsMaximal := by
    rcases Ideal.map_eq_top_or_isMaximal_of_surjective f hsurj hn with h | h
    · exact absurd h hne
    · exact h
  have heq : Ideal.map f n = IsLocalRing.maximalIdeal (A ⧸ I) :=
    IsLocalRing.eq_maximalIdeal hmax'
  rw [hm, ← heq, hcm]

/-- ★★★★★★★★★★★★★★★★★★★★
**完備な整域は商が Artin なら局所**——★**無条件**（第 1364）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★★★これが第 1357-1363 の 7 段を繋いだ**使う形**である。 -/
theorem isLocalRing_of_isAdicComplete_of_domain [IsDomain A] (I : Ideal A)
    [IsAdicComplete I A] (hI : I ≤ Ideal.jacobson (⊥ : Ideal A))
    [IsArtinianRing (A ⧸ I)] [Nontrivial (A ⧸ I)] : IsLocalRing A := by
  haveI : IsLocalRing (A ⧸ I) :=
    isLocalRing_of_isArtinian_of_no_idempotent
      (fun x hx => isIdempotentElem_quotient_eq_zero_or_one I hx)
  exact isLocalRing_of_quotient I hI

/-! ## ★出典の紐付け(`.src`) -/

def isLocalRing_of_quotient.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(商が局所で I が Jacobson 根基に入っていれば元も局所。★無条件)",
    sectionId := "genell-thm-3-8" }

def isLocalRing_of_isAdicComplete_of_domain.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(完備な整域は商が Artin なら局所。★無条件)",
    sectionId := "genell-thm-3-8" }

def isLocalRing_of_isAdicComplete_of_domain.needs : List ProofObligation :=
  [ .citation "[ABC3]" "isIdempotentElem_quotient_eq_zero_or_one(第 1358、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.isIdempotentElem_quotient_eq_zero_or_one") 1,
    .citation "[ABC3]" "isLocalRing_of_isArtinian_of_no_idempotent(第 1362、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.isLocalRing_of_isArtinian_of_no_idempotent") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1364）**——第 1357-1363 の 7 段を繋いだ**使う形**である。" ++
       "☆`R′ = integralClosure R L` に当てれば「完備 DVR の有限整閉包は局所」が出る。") 19 ]

end ABC3.Found.GenEll
