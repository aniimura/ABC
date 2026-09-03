/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.GalTateMatrixUnip
import ABC3.Found.GaloisRep.TateLevelOne
import ABC3.Found.GaloisRep.TateKerLevel
import ABC3.Found.GenEll.PointTransport
import ABC3.Meta.Claim

/-!
# 第 1294 ブロック —— **`σ` が固定する位数 `l` の点は固定ベクトルを与える**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★★★★★これは何か——第 1292 のもう一方の入力

第 1292 は「`det = 1` かつ**固定ベクトルがある**」から幂単性を出す。
★本ブロックはその**固定ベクトル**を、`σ` が固定する位数 `l` の点から作る。

☆機構は 3 つ:

| 段 | 内容 | 在庫 |
|---|---|---|
| 1 | 位数 `l` の点は `T_l E` の第 1 層に持ち上がる | 第 1202 |
| 2 | 持ち上げの `mod l` 座標は `0` でない（核は `l·T_l E`） | 第 1203・1204 |
| 3 | `σ` が点を固定するなら座標も固定される | 第 1233 |

★★★これで `μ_l ⊂ E(K₀)` から「行列は直線を点ごとに固定する」が出る。
-/

namespace ABC3.Found.GenEll

open WeierstrassCurve WeierstrassCurve.Affine ABC3.Interface.GaloisRep ABC3.Found.GaloisRep
open ABC3.Meta

open scoped Matrix

variable {K L : Type} [Field K] [DecidableEq K] [Field L] [DecidableEq L] [Algebra K L]

set_option maxHeartbeats 800000 in
/-- ★★★★★★★★★★★★★★★★
**`σ` が固定する位数 `l` の点は固定ベクトルを与える**——★**無条件**（第 1294）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

☆点を第 1 層に持ち上げ、その `mod l` 座標を取る。 -/
theorem exists_fixed_vec_of_galPoint_eq [IsAlgClosed L] [CharZero L]
    (W : WeierstrassCurve K) [((W.baseChange L).toAffine).IsElliptic]
    (l : ℕ) [Fact l.Prime]
    (e : tateModule (W.baseChange L) l ≃+ (Fin 2 → ℤ_[l])) (σ : L ≃ₐ[K] L)
    (P : (W.baseChange L).toAffine.Point) (hP : addOrderOf P = l)
    (hfix : galPoint W σ P = P) :
    ∃ v : Fin 2 → ZMod l, v ≠ 0 ∧
      ((glRedPadic l (galRep W l e σ) : GL (Fin 2) (ZMod l)) :
        Matrix (Fin 2) (Fin 2) (ZMod l)) *ᵥ v = v := by
  have hl : l.Prime := Fact.out
  have hPl : l • P = 0 := by
    rw [← hP]
    exact addOrderOf_nsmul_eq_zero P
  obtain ⟨f, hf⟩ := exists_tateProj_one_eq W l P hPl
  refine ⟨redVec l (e f), ?_, ?_⟩
  · intro hzero
    obtain ⟨u, hu⟩ := (redVec_eq_zero_iff l (e f)).1 hzero
    have hfl : ∃ g : tateModule (W.baseChange L) l, l • g = f := by
      refine ⟨e.symm u, ?_⟩
      have h2 := congrArg e.symm hu
      rw [e.symm_apply_apply, map_nsmul] at h2
      exact h2.symm
    have hP0 : P = 0 := by
      rw [← hf]
      exact (tateProj_one_eq_zero_iff (W.baseChange L) l f).2 hfl
    exact (ne_zero_of_addOrderOf_prime hl hP) hP0
  · have hproj : ((tateProj (W.baseChange L) l 1 (galTate W l σ f - f) :
        (W.baseChange L).toAffine.Point)) = 0 := by
      rw [map_sub]
      show ((tateProj (W.baseChange L) l 1 (galTate W l σ f)) :
        (W.baseChange L).toAffine.Point) - ((tateProj (W.baseChange L) l 1 f) :
        (W.baseChange L).toAffine.Point) = 0
      rw [tateProj_galTate, hf, hfix, sub_self]
    obtain ⟨g, hg⟩ := (tateProj_one_eq_zero_iff (W.baseChange L) l _).1 hproj
    have hdiff : e (galTate W l σ f) - e f = l • e g := by
      rw [← map_sub, ← hg, map_nsmul]
    have hred : redVec l (e (galTate W l σ f)) = redVec l (e f) := by
      have h3 : redVec l (e (galTate W l σ f) - e f) = 0 := by
        rw [hdiff, redVec_nsmul_self]
      rw [redVec_sub] at h3
      exact sub_eq_zero.mp h3
    rw [← redVec_galTate, hred]

/-! ## ★出典の紐付け(`.src`) -/

def exists_fixed_vec_of_galPoint_eq.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(σ が固定する位数 l の点は固定ベクトルを与える。★無条件)",
    sectionId := "genell-thm-3-8" }

def exists_fixed_vec_of_galPoint_eq.needs : List ProofObligation :=
  [ .citation "[ABC3]" "exists_tateProj_one_eq(第 1202、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.exists_tateProj_one_eq") 1,
    .citation "[ABC3]" "tateProj_one_eq_zero_iff(第 1203、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.tateProj_one_eq_zero_iff") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1294）**——第 1292 のもう一方の入力である。" ++
       "☆`μ_l ⊂ E(K₀)` から「行列は直線を点ごとに固定する」が出る。") 2 ]

end ABC3.Found.GenEll
