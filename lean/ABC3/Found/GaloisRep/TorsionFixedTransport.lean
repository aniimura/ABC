/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.TorsionTransport
import ABC3.Meta.Claim

/-!
# 第 1302 ブロック —— **固定点も埋め込みで降りる**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★これは何か——第 1271 の続き（`h2` の入力を運ぶ）

第 1271 では**非自明性**（動かされる点）を降ろした。
★本ブロックは**固定点**を降ろす——`h2`（第 1296）が要求するのはそちらである。

☆道具は同じ: `torsionMap_bijective`（個数の勘定）と `point_map_galPoint`（同変性）。
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine ABC3.Interface.GaloisRep ABC3.Meta

variable {K₀ F K : Type} [Field K₀] [DecidableEq K₀] [Field F] [DecidableEq F]
  [Field K] [DecidableEq K] [Algebra K₀ F] [Algebra K₀ K]

set_option linter.unusedSectionVars false in
/-- ★★★★★★★★★★★★★★★★
**固定点も埋め込みで降りる**——★**無条件**（第 1302）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

☆全単射（第 1271）で引き戻し、同変性で固定を移す。 -/
theorem exists_galPoint_fixed_of_map [IsAlgClosed F] [IsAlgClosed K] (W : WeierstrassCurve K₀)
    (hΔF : (W.baseChange F).Δ ≠ 0) (hΔK : (W.baseChange K).Δ ≠ 0)
    (f : F →ₐ[K₀] K) (σF : F ≃ₐ[K₀] F) (σK : K ≃ₐ[K₀] K)
    (hcomm : ∀ x : F, f (σF x) = σK (f x)) (n : ℕ) (hn : 1 ≤ n)
    (hcF : ∀ k : ℕ, 1 ≤ k → k ≤ n → ((k : F) ≠ 0))
    (hcK : ∀ k : ℕ, 1 ≤ k → k ≤ n → ((k : K) ≠ 0))
    (Q : (W.baseChange K).toAffine.Point) (hQn : addOrderOf Q = n)
    (hQfix : galPoint W σK Q = Q) :
    ∃ P : (W.baseChange F).toAffine.Point, addOrderOf P = n ∧ galPoint W σF P = P := by
  have hQ0 : n • Q = 0 := by
    rw [← hQn]
    exact addOrderOf_nsmul_eq_zero Q
  obtain ⟨P, hP⟩ := (torsionMap_bijective W hΔF hΔK f n hn hcF hcK).2 ⟨Q, hQ0⟩
  have hPQ : Point.map f ((P : (W.baseChange F).toAffine.Point)) = Q :=
    congrArg Subtype.val hP
  refine ⟨(P : (W.baseChange F).toAffine.Point), ?_, ?_⟩
  · have hord : addOrderOf (Point.map f ((P : (W.baseChange F).toAffine.Point)))
        = addOrderOf ((P : (W.baseChange F).toAffine.Point)) :=
      addOrderOf_injective (Point.map (W' := W) f) (Point.map_injective (f := f)) _
    rw [hPQ] at hord
    rw [← hord, hQn]
  · refine Point.map_injective (f := f) ?_
    rw [point_map_galPoint W f σF σK hcomm, hPQ, hQfix]

/-! ## ★出典の紐付け(`.src`) -/

def exists_galPoint_fixed_of_map.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(固定点も埋め込みで降りる。★無条件)",
    sectionId := "genell-thm-3-8" }

def exists_galPoint_fixed_of_map.needs : List ProofObligation :=
  [ .citation "[ABC3]" "torsionMap_bijective(第 1271、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.torsionMap_bijective") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1302）**——第 1271 の続きである。" ++
       "☆`h2`（第 1296）が要求する**固定点**を局所から大域へ運ぶ。") 2 ]

end ABC3.Found.GaloisRep
