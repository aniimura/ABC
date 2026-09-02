/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import Mathlib.RingTheory.Artinian.Ring
import Mathlib.RingTheory.Artinian.Module
import Mathlib.RingTheory.Idempotents
import Mathlib.RingTheory.Jacobson.Semiprimary
import ABC3.Meta.Claim

/-!
# 第 1362 ブロック —— **Artin 環は非自明な冪等元が無ければ局所**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★★★★★これは何か——第 1357 の道の第 6 段

「完備 DVR `R` の有限整閉包 `R′` は局所」を出すには

* `R′/m_R R′` は Artin 環（剰余体上有限次元）
* `R′` は整域なので非自明な冪等元を持たない → 第 1358 より商も持たない
* **Artin ＋ 非自明な冪等元なし ⟹ 局所**（★本ブロック）

☆道は 3 段:

1. `J ≔ Ring.jacobson A` は冪零（`IsArtinianRing.isNilpotent_jacobson_bot`、在庫）
2. `B ≔ A ⧸ J` は**半単純**（`IsSemiprimaryRing.isSemisimpleRing`、在庫）で、
   冪等元は `J` に沿って `A` へ持ち上がる（`exists_isIdempotentElem_eq_of_ker_isNilpotent`、在庫）
3. 可換半単純環で非自明な冪等元が無ければ**体**——極大イデアル `m` の補元 `N` から
   `1 = e + f`（`e ∈ m`, `f ∈ N`）を取ると `e` は冪等元。`e = 1` は `m ≠ ⊤` に反し、
   `e = 0` なら `N = ⊤` より `m = ⊥`。

★したがって `J` は極大、`J` は全ての極大イデアルに含まれるので `A` は局所。
-/

namespace ABC3.Found.GenEll

open ABC3.Meta

variable {A : Type*} [CommRing A]

/-- ★★★★★★★★
**可換半単純環で非自明な冪等元が無ければ `⊥` は極大**——★**無条件**（第 1362）。 -/
theorem bot_isMaximal_of_isSemisimpleRing [IsSemisimpleRing A] [Nontrivial A]
    (h : ∀ e : A, e * e = e → e = 0 ∨ e = 1) {m : Ideal A} (hm : m.IsMaximal) :
    m = ⊥ := by
  obtain ⟨N, hN⟩ := ComplementedLattice.exists_isCompl (α := Submodule A A) m
  have h1 : (1 : A) ∈ m ⊔ N := by
    rw [hN.sup_eq_top]
    trivial
  obtain ⟨e, he, f, hf, hef⟩ := Submodule.mem_sup.1 h1
  have hzero : e * f = 0 := by
    have hmem : e * f ∈ m ⊓ N := by
      refine ⟨?_, ?_⟩
      · exact Ideal.mul_mem_right _ _ he
      · exact Submodule.smul_mem N e hf
    rw [hN.inf_eq_bot] at hmem
    simpa using hmem
  have hidem : e * e = e := by
    have : e * (e + f) = e * 1 := by rw [hef]
    rw [mul_add, hzero, add_zero, mul_one] at this
    exact this
  rcases h e hidem with h0 | h1'
  · -- ★ `e = 0` なら `f = 1`、よって `N = ⊤`、よって `m = ⊥`
    have hf1 : f = 1 := by rw [← hef, h0, zero_add]
    have hNtop : N = ⊤ := by
      refine Submodule.eq_top_iff'.2 fun x => ?_
      have : x • f ∈ N := Submodule.smul_mem N x hf
      simpa [hf1] using this
    have : m = m ⊓ N := by rw [hNtop, inf_top_eq]
    rw [this, hN.inf_eq_bot]
  · -- ★ `e = 1` は `m ≠ ⊤` に反する
    exact absurd (Ideal.eq_top_of_isUnit_mem m (h1' ▸ he) isUnit_one) hm.ne_top

/-- ★★★★★★★★★★★★★★★★★★★★
**Artin 環は非自明な冪等元が無ければ局所**——★**無条件**（第 1362）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★★★これが第 1357 の道の**第 6 段**である。 -/
theorem isLocalRing_of_isArtinian_of_no_idempotent [IsArtinianRing A] [Nontrivial A]
    (h : ∀ e : A, e * e = e → e = 0 ∨ e = 1) : IsLocalRing A := by
  classical
  set J : Ideal A := Ideal.jacobson (⊥ : Ideal A) with hJ
  -- ★ 段 1: `J` は冪零
  obtain ⟨n, hn⟩ := IsArtinianRing.isNilpotent_jacobson_bot (R := A)
  have hnil : ∀ x ∈ J, IsNilpotent x := by
    intro x hx
    refine ⟨n, ?_⟩
    have : x ^ n ∈ J ^ n := Ideal.pow_mem_pow hx n
    rw [hn] at this
    simpa using this
  -- ★ 段 2: `B = A ⧸ J` の冪等元は `A` へ持ち上がる
  have hjle : ∀ {m : Ideal A}, m.IsMaximal → J ≤ m := by
    intro m hm
    exact sInf_le ⟨bot_le, hm⟩
  have hJne : J ≠ ⊤ := by
    obtain ⟨m, hm, -⟩ := Ideal.exists_le_maximal (⊥ : Ideal A) bot_ne_top
    intro htop
    exact hm.ne_top (top_le_iff.mp (htop ▸ hjle hm))
  haveI : Nontrivial (A ⧸ J) := Ideal.Quotient.nontrivial_iff.mpr hJne
  have hidemB : ∀ x : A ⧸ J, x * x = x → x = 0 ∨ x = 1 := by
    intro x hx
    have hker : ∀ y ∈ RingHom.ker (Ideal.Quotient.mk J), IsNilpotent y := by
      intro y hy
      exact hnil y (Ideal.Quotient.eq_zero_iff_mem.mp hy)
    obtain ⟨e, he, hemk⟩ := exists_isIdempotentElem_eq_of_ker_isNilpotent
      (Ideal.Quotient.mk J) hker x (Ideal.Quotient.mk_surjective x) hx
    rcases h e he with h0 | h1
    · exact Or.inl (by rw [← hemk, h0, map_zero])
    · exact Or.inr (by rw [← hemk, h1, map_one])
  -- ★ 段 3: `B` は半単純で非自明な冪等元が無いので体
  haveI : IsSemisimpleRing (A ⧸ Ring.jacobson A) := IsSemiprimaryRing.isSemisimpleRing
  have hBsemi : IsSemisimpleRing (A ⧸ J) := by
    rw [hJ, Ideal.jacobson_bot]
    infer_instance
  haveI := hBsemi
  have hbot : ∀ {m : Ideal (A ⧸ J)}, m.IsMaximal → m = ⊥ := fun {m} hm =>
    bot_isMaximal_of_isSemisimpleRing hidemB hm
  have hJmax : J.IsMaximal := by
    obtain ⟨m, hm, -⟩ := Ideal.exists_le_maximal (⊥ : Ideal (A ⧸ J)) bot_ne_top
    have hmb : m = ⊥ := hbot hm
    have : (⊥ : Ideal (A ⧸ J)).IsMaximal := hmb ▸ hm
    exact (Ideal.Quotient.maximal_ideal_iff_isField_quotient J).2
      (Ring.isField_iff_maximal_bot.2 this)
  -- ★ `J` は全ての極大イデアルに含まれるので、`J` が唯一の極大イデアル
  refine IsLocalRing.of_unique_max_ideal ⟨J, hJmax, fun m hm => ?_⟩
  exact (hJmax.eq_of_le hm.ne_top (hjle hm)).symm

/-! ## ★出典の紐付け(`.src`) -/

def bot_isMaximal_of_isSemisimpleRing.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(可換半単純環で非自明な冪等元が無ければ ⊥ は極大。★無条件)",
    sectionId := "genell-thm-3-8" }

def isLocalRing_of_isArtinian_of_no_idempotent.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Artin 環は非自明な冪等元が無ければ局所。★無条件)",
    sectionId := "genell-thm-3-8" }

def isLocalRing_of_isArtinian_of_no_idempotent.needs : List ProofObligation :=
  [ .implicitStep
      ("★★★★**2026-09-02（第 1362）**——第 1357 の道の**第 6 段**である。" ++
       "☆`IsArtinianRing.isNilpotent_jacobson_bot`・`IsSemiprimaryRing.isSemisimpleRing`・" ++
       "`exists_isIdempotentElem_eq_of_ker_isNilpotent` はすべて mathlib の在庫。") 19 ]

end ABC3.Found.GenEll
