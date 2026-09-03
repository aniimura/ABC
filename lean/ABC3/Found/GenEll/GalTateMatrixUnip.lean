/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.AlphaMatrix
import ABC3.Found.GenEll.UnipotentDet
import Mathlib.Tactic.NoncommRing
import ABC3.Meta.Claim

/-!
# 第 1293 ブロック —— **行列が幂単なら `T_l E` でも幂単**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★★★★★これは何か——第 1292 の受け口

第 1234（`glRed_unipotent_of_galTate`）は「`σ` が幂単 ⇒ 行列も幂単」だった。
★本ブロックは**逆向き**である——`(A − 1)² = 0` から `h2` を出す。

☆機構は `redVec_galTate`（第 1233、行列は `mod l` 座標に `mulVec` で作用）と
`redVec_eq_zero_iff`（第 1204、核はちょうど `l · ℤ_l²`）だけ。

★★★これで第 1292（`det = 1` ＋ 固定ベクトル ⇒ 幂単）が `h2` に繋がる。
-/

namespace ABC3.Found.GenEll

open WeierstrassCurve WeierstrassCurve.Affine ABC3.Interface.GaloisRep ABC3.Found.GaloisRep ABC3.Meta

open scoped Matrix

variable {K L : Type} [Field K] [DecidableEq K] [Field L] [DecidableEq L] [Algebra K L]

set_option maxHeartbeats 800000 in
/-- ★★★★★★★★★★★★★★★★
**行列が幂単なら `T_l E` でも `mod l` 幂単**——★**無条件**（第 1293）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

☆`redVec` で `mod l` に落とし、核が `l · T_l E` であることを使う。 -/
theorem galTate_unipotent_of_matrix (W : WeierstrassCurve K) (l : ℕ) [Fact l.Prime]
    (e : tateModule (W.baseChange L) l ≃+ (Fin 2 → ℤ_[l])) (σ : L ≃ₐ[K] L)
    (hA : (((glRedPadic l (galRep W l e σ) : GL (Fin 2) (ZMod l)) :
          Matrix (Fin 2) (Fin 2) (ZMod l)) - 1) *
        (((glRedPadic l (galRep W l e σ) : GL (Fin 2) (ZMod l)) :
          Matrix (Fin 2) (Fin 2) (ZMod l)) - 1) = 0) :
    ∀ x : tateModule (W.baseChange L) l, ∃ u : tateModule (W.baseChange L) l,
      galTate W l σ (galTate W l σ x) + x
        = galTate W l σ x + galTate W l σ x + l • u := by
  intro x
  have key : ∀ (B : Matrix (Fin 2) (Fin 2) (ZMod l)) (w : Fin 2 → ZMod l),
      B *ᵥ (B *ᵥ w) + w - (B *ᵥ w + B *ᵥ w) = ((B - 1) * (B - 1)) *ᵥ w := by
    intro B w
    have h1 : (B - 1) * (B - 1) = B * B - B - B + 1 := by noncomm_ring
    rw [h1]
    simp only [Matrix.add_mulVec, Matrix.sub_mulVec, Matrix.one_mulVec, Matrix.mulVec_mulVec]
    abel
  have hred : redVec l (e (galTate W l σ (galTate W l σ x) + x
      - (galTate W l σ x + galTate W l σ x))) = 0 := by
    rw [map_sub, map_add, map_add, redVec_sub, redVec_add, redVec_add,
      redVec_galTate, redVec_galTate, key, hA, Matrix.zero_mulVec]
  obtain ⟨u, hu⟩ := (redVec_eq_zero_iff l _).1 hred
  refine ⟨e.symm u, ?_⟩
  have hY : galTate W l σ (galTate W l σ x) + x - (galTate W l σ x + galTate W l σ x)
      = l • e.symm u := by
    have h2 := congrArg e.symm hu
    rw [e.symm_apply_apply, map_nsmul] at h2
    exact h2
  have h3 := sub_eq_iff_eq_add.mp hY
  rw [h3]
  abel

/-! ## ★出典の紐付け(`.src`) -/

def galTate_unipotent_of_matrix.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(行列が幂単なら T_l E でも mod l 幂単。★無条件)",
    sectionId := "genell-thm-3-8" }

def galTate_unipotent_of_matrix.needs : List ProofObligation :=
  [ .citation "[ABC3]" "redVec_galTate(第 1233、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.redVec_galTate") 1,
    .citation "[ABC3]" "redVec_eq_zero_iff(第 1204、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.redVec_eq_zero_iff") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1293）**——第 1234 の**逆向き**である。" ++
       "☆これで第 1292（`det = 1` ＋ 固定ベクトル ⇒ 幂単）が `h2` に繋がる。") 2 ]

end ABC3.Found.GenEll
