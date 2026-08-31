/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.GalRep
import ABC3.Found.GaloisRep.TranslateEquiv
import ABC3.Interface.GaloisRep.Torsion
import Mathlib.FieldTheory.KrullTopology
import ABC3.Meta.Claim

/-!
# `galRep` の連続性へ —— 座標を固定する `σ` は点を動かさない（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★これは何か

`Found/GenEll/GalRepClosed.lean`（`§9-1191`、第 765）で、
`Theorem 3.8` に残る位相の側は **`galRep` の連続性**だけになった。
その連続性の証明は

1. `σ` が `E[l^n]` の座標をすべて固定するなら `σ` は `E[l^n]` に自明に作用する ←**本ファイル**
2. `L(E[l^n])` は `L` の有限拡大であり、それを固定する部分群は Krull 位相で開
3. したがって `galMat σ ≡ 1 (mod l^n)` は開集合の上で成り立ち、`galRep` は連続

という 3 段である。本ファイルは**第 1 段**を取る。
-/

namespace ABC3.Found.GenEll

open ABC3.Found.GaloisRep ABC3.Interface.GaloisRep WeierstrassCurve WeierstrassCurve.Affine

variable {K L : Type} [Field K] [DecidableEq K] [Field L] [DecidableEq L] [Algebra K L]

/-- ★★★★★★**座標を固定する `σ` は点を動かさない**——★**無条件**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★`galPoint W σ = Point.map σ` なので、`some x y` の像は `some (σ x) (σ y)` である。 -/
theorem galPoint_eq_of_fixed (W : WeierstrassCurve K) (σ : L ≃ₐ[K] L)
    (P : (W.baseChange L).toAffine.Point)
    (h : ∀ (x y : L) (hns : (W.baseChange L).toAffine.Nonsingular x y),
      P = Point.some x y hns → σ x = x ∧ σ y = y) :
    galPoint W σ P = P := by
  cases P with
  | zero => exact Point.map_zero _
  | some x y hns =>
    obtain ⟨hx, hy⟩ := h x y hns rfl
    rw [galPoint, Point.map_some]
    exact point_some_congr hx hy _ _

/-- ★★★★★**すべての点で座標を固定するなら `galPoint` は恒等写像**。 -/
theorem galPoint_eq_id_of_fixed (W : WeierstrassCurve K) (σ : L ≃ₐ[K] L)
    (S : Set ((W.baseChange L).toAffine.Point))
    (h : ∀ P ∈ S, ∀ (x y : L) (hns : (W.baseChange L).toAffine.Nonsingular x y),
      P = Point.some x y hns → σ x = x ∧ σ y = y) :
    ∀ P ∈ S, galPoint W σ P = P :=
  fun P hP => galPoint_eq_of_fixed W σ P (h P hP)

/-! ## ★★★★★★★★捩れ点の座標が生成する有限拡大 -/

/-- ★点の座標の集合（`0` では空）。 -/
def ptCoordSet (W : WeierstrassCurve K) : (W.baseChange L).toAffine.Point → Set L
  | .zero => ∅
  | .some x y _ => {x, y}

theorem finite_ptCoordSet (W : WeierstrassCurve K) (P : (W.baseChange L).toAffine.Point) :
    (ptCoordSet W P).Finite := by
  cases P with
  | zero => exact Set.finite_empty
  | some x y _ => exact (Set.finite_singleton y).insert x

/-- ★★`σ` が座標集合を固定するなら点を動かさない。 -/
theorem galPoint_eq_of_fixed_coordSet (W : WeierstrassCurve K) (σ : L ≃ₐ[K] L)
    (P : (W.baseChange L).toAffine.Point) (h : ∀ z ∈ ptCoordSet W P, σ z = z) :
    galPoint W σ P = P := by
  refine galPoint_eq_of_fixed W σ P ?_
  intro x y hns hP
  subst hP
  refine ⟨h x ?_, h y ?_⟩
  · exact Or.inl rfl
  · exact Or.inr rfl

/-- ★★★**`n`-捩れ点の座標全体**。 -/
def torsionCoordSet (W : WeierstrassCurve K) (n : ℕ) : Set L :=
  ⋃ P ∈ (torsionPoints (W.baseChange L) n : Set ((W.baseChange L).toAffine.Point)),
    ptCoordSet W P

theorem finite_torsionCoordSet (W : WeierstrassCurve K) (n : ℕ)
    (hfin : (torsionPoints (W.baseChange L) n : Set ((W.baseChange L).toAffine.Point)).Finite) :
    (torsionCoordSet W (L := L) n).Finite :=
  hfin.biUnion (fun P _ => finite_ptCoordSet (L := L) W P)

/-- ★★★★★★★★★★**捩れ点の座標が生成する有限拡大**——その固定部分群は
`n`-捩れに自明に作用する。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★★これが `galRep` の連続性の**第 2 段**である
（Krull 位相の開部分群は有限次中間体の `fixingSubgroup` だから）。 -/
theorem exists_finiteDimensional_fixing_torsion [Algebra.IsAlgebraic K L]
    (W : WeierstrassCurve K) (n : ℕ)
    (hfin : (torsionPoints (W.baseChange L) n : Set ((W.baseChange L).toAffine.Point)).Finite) :
    ∃ F : IntermediateField K L, FiniteDimensional K F ∧
      ∀ σ ∈ F.fixingSubgroup, ∀ P ∈ torsionPoints (W.baseChange L) n,
        galPoint W σ P = P := by
  haveI : Finite (torsionCoordSet W (L := L) n) := (finite_torsionCoordSet W n hfin).to_subtype
  refine ⟨IntermediateField.adjoin K (torsionCoordSet W (L := L) n), ?_, ?_⟩
  · exact IntermediateField.finiteDimensional_adjoin
      (fun x _ => Algebra.IsIntegral.isIntegral x)
  · intro σ hσ P hP
    refine galPoint_eq_of_fixed_coordSet W σ P (fun z hz => ?_)
    have hmem : z ∈ IntermediateField.adjoin K (torsionCoordSet W (L := L) n) := by
      refine IntermediateField.subset_adjoin _ _ ?_
      exact Set.mem_biUnion hP hz
    exact (IntermediateField.mem_fixingSubgroup_iff _ _).1 hσ z hmem

/-! ## ★出典の紐付け(`.src`) -/

def exists_finiteDimensional_fixing_torsion.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(捩れ点の座標が生成する有限拡大——その固定部分群は捩れに自明に作用)",
    sectionId := "genell-thm-3-8" }

def galPoint_eq_of_fixed.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(座標を固定する σ は点を動かさない。★無条件)",
    sectionId := "genell-thm-3-8" }

def galPoint_eq_of_fixed.needs : List ABC3.Meta.ProofObligation :=
  [ .implicitStep
      ("☆連続性の残り 2 段: (2) L(E[l^n]) は L の有限拡大でありそれを固定する部分群は " ++
       "Krull 位相で開、(3) したがって galMat σ ≡ 1 (mod l^n) は開集合の上で成り立つ") 6 ]

end ABC3.Found.GenEll
