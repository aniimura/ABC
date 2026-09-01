/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.GalRepBasis
import ABC3.Meta.Claim

/-!
# 第 1231 ブロック —— **任意の `GL₂(ℤ_l)` の元は基底の取り替えで実現できる**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★★★★★これは何か

第 1229-1230 で「`σ` の `mod l` の行列は `α` に共役」「共役行列は `ℤ_l` に持ち上がる」
が取れた。☆その持ち上げ `B` を**実際に基底の取り替えとして実行する**には、
`B = basisChange e₀ e` なる `e` が要る。

★`e ≔ e₀` の後に `B` を掛けたものを取ればよい。
-/

namespace ABC3.Found.GaloisRep

open ABC3.Interface.GaloisRep WeierstrassCurve ABC3.Meta

variable {K L : Type} [Field K] [DecidableEq K] [Field L] [DecidableEq L] [Algebra K L]

/-- ☆`GL₂(ℤ_l)` の元は `ℤ_l²` の加法同型を与える。 -/
noncomputable def mulVecAddEquiv (l : ℕ) [Fact l.Prime] (B : GL (Fin 2) ℤ_[l]) :
    (Fin 2 → ℤ_[l]) ≃+ (Fin 2 → ℤ_[l]) where
  toFun v := (B : Matrix (Fin 2) (Fin 2) ℤ_[l]).mulVec v
  invFun v := ((B⁻¹ : GL (Fin 2) ℤ_[l]) : Matrix (Fin 2) (Fin 2) ℤ_[l]).mulVec v
  left_inv := fun v => by
    have h : ((B⁻¹ : GL (Fin 2) ℤ_[l]) : Matrix (Fin 2) (Fin 2) ℤ_[l])
        * (B : Matrix (Fin 2) (Fin 2) ℤ_[l]) = 1 := by
      simpa using congrArg (fun g : GL (Fin 2) ℤ_[l] =>
        (g : Matrix (Fin 2) (Fin 2) ℤ_[l])) (inv_mul_cancel B)
    show ((B⁻¹ : GL (Fin 2) ℤ_[l]) : Matrix (Fin 2) (Fin 2) ℤ_[l]).mulVec
      ((B : Matrix (Fin 2) (Fin 2) ℤ_[l]).mulVec v) = v
    rw [Matrix.mulVec_mulVec, h, Matrix.one_mulVec]
  right_inv := fun v => by
    have h : (B : Matrix (Fin 2) (Fin 2) ℤ_[l])
        * ((B⁻¹ : GL (Fin 2) ℤ_[l]) : Matrix (Fin 2) (Fin 2) ℤ_[l]) = 1 := by
      simpa using congrArg (fun g : GL (Fin 2) ℤ_[l] =>
        (g : Matrix (Fin 2) (Fin 2) ℤ_[l])) (mul_inv_cancel B)
    show (B : Matrix (Fin 2) (Fin 2) ℤ_[l]).mulVec
      (((B⁻¹ : GL (Fin 2) ℤ_[l]) : Matrix (Fin 2) (Fin 2) ℤ_[l]).mulVec v) = v
    rw [Matrix.mulVec_mulVec, h, Matrix.one_mulVec]
  map_add' := fun x y => Matrix.mulVec_add _ _ _

/-- ★★★★★★★★★★★★
**任意の `GL₂(ℤ_l)` の元は基底の取り替えで実現できる**——★**無条件**（第 1231）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

☆`e ≔ e₀` の後に `B` を掛けたものを取ればよい。

★★★これで第 1229-1230 の共役を実際の基底の取り替えとして実行できる。 -/
theorem basisChange_realize (W : WeierstrassCurve K) (l : ℕ) [Fact l.Prime]
    (e₀ : tateModule (W.baseChange L) l ≃+ (Fin 2 → ℤ_[l])) (B : GL (Fin 2) ℤ_[l]) :
    basisChange W l e₀ (e₀.trans (mulVecAddEquiv l B))
      = (B : Matrix (Fin 2) (Fin 2) ℤ_[l]) := by
  refine matrix_ext_of_mulVec (fun x => ?_)
  have h := basisChange_apply W l e₀ (e₀.trans (mulVecAddEquiv l B)) (e₀.symm x)
  simpa [mulVecAddEquiv] using h.symm

/-! ## ★出典の紐付け(`.src`) -/

def mulVecAddEquiv.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(GL₂(Z_l) の元が与える Z_l² の加法同型)",
    sectionId := "genell-thm-3-8" }

def basisChange_realize.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(任意の GL₂(Z_l) の元は基底の取り替えで実現できる。★無条件)",
    sectionId := "genell-thm-3-8" }

def basisChange_realize.needs : List ProofObligation :=
  [ .implicitStep
      ("★★★★**2026-09-02（第 1231）**——第 1229-1230 の共役を" ++
       "**実際の基底の取り替えとして実行する**ための段である。" ++
       "☆`e ≔ e₀` の後に `B` を掛けたものを取ればよい。") 2 ]

end ABC3.Found.GaloisRep
