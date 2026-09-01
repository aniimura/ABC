/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.PadicRedVec
import ABC3.Found.GaloisRep.GalRep
import ABC3.Meta.Claim

/-!
# 第 1228 ブロック —— **座標の作用から行列を読む**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★★★★★★★★★これは何か——α 橋の II 側の**最後の読み替え**

第 1212（`kummer_sigma_coord_alpha`）は `σ` が `mod l` の座標で
`α = (1 1 / 0 1)` として作用することを与える。

☆本ブロックはそれを **`glRedPadic l (galRep σ)` という行列の等式**に読み替える:

    ∀ x, redVec (e (galTate σ x)) = α ·ᵥ redVec (e x)
      ⟹  glRedPadic l (galRep σ) = α

★材料は `galMat_apply`（在庫）・`redVec_mulVec`（第 1204）・
`exists_redVec_eq`（第 1204、還元は全射）・
mathlib の `Matrix.ext_of_mulVec_single`。

★★★これで `α ∈ (galRep の像).map (glRedPadic l)` が言える。
-/

namespace ABC3.Found.GenEll

open ABC3.Found.GaloisRep ABC3.Interface.GaloisRep WeierstrassCurve ABC3.Meta

variable {K L : Type} [Field K] [DecidableEq K] [Field L] [DecidableEq L] [Algebra K L]

/-- ★★★★★★★★★★★★★★★★
**座標の作用から行列を読む**——★**無条件**（第 1228）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

☆`σ` が `mod l` の座標で `α` として作用するなら、
`glRedPadic l (galRep σ)` はちょうど `α` である。

★★★これが `α ∈ (galRep の像).map (glRedPadic l)` を出す段である。 -/
theorem glRed_galRep_eq_of_redVec (l : ℕ) [Fact l.Prime]
    (W : WeierstrassCurve K) (e : tateModule (W.baseChange L) l ≃+ (Fin 2 → ℤ_[l]))
    (σ : L ≃ₐ[K] L) (α : Matrix (Fin 2) (Fin 2) (ZMod l))
    (h : ∀ x : tateModule (W.baseChange L) l,
      redVec l (e (galTate W l σ x)) = α.mulVec (redVec l (e x))) :
    ((glRedPadic l (galRep W l e σ) : GL (Fin 2) (ZMod l)) :
      Matrix (Fin 2) (Fin 2) (ZMod l)) = α := by
  have hcoe : ((glRedPadic l (galRep W l e σ) : GL (Fin 2) (ZMod l)) :
      Matrix (Fin 2) (Fin 2) (ZMod l))
      = (galMat W l e σ).map PadicInt.toZMod := by
    rw [coe_glRedPadic, galRep_coe]
  rw [hcoe]
  refine Matrix.ext_of_mulVec_single ?_
  intro i
  obtain ⟨w, hw⟩ := exists_redVec_eq l (Pi.single i 1)
  have hex : e (e.symm w) = w := e.apply_symm_apply w
  have hgal : e (galTate W l σ (e.symm w)) = (galMat W l e σ).mulVec w := by
    rw [galMat_apply, hex]
  have hkey := h (e.symm w)
  rw [hgal, hex, redVec_mulVec, hw] at hkey
  exact hkey

/-! ## ★出典の紐付け(`.src`) -/

def glRed_galRep_eq_of_redVec.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(座標の作用から行列を読む。★無条件)",
    sectionId := "genell-thm-3-8" }

def glRed_galRep_eq_of_redVec.needs : List ProofObligation :=
  [ .citation "[ABC3]" "redVec_mulVec(第 1204、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.redVec_mulVec") 1,
    .citation "[ABC3]" "exists_redVec_eq(第 1204、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_redVec_eq") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1228）**——第 1212（`kummer_sigma_coord_alpha`）が" ++
       "与える「`σ` は `mod l` の座標で `α` として作用する」を、" ++
       "**`glRedPadic l (galRep σ)` という行列の等式**に読み替える段である。" ++
       "★★★これで `α ∈ (galRep の像).map (glRedPadic l)` が言える。") 3 ]

end ABC3.Found.GenEll
