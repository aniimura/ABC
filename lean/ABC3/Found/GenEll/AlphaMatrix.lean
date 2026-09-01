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

/-! ## ★★★★★★★★★★★★逆向き——行列は座標の作用を与える -/

/-- ★★★★★★★★★★★★
**`glRedPadic (galRep σ)` は `mod l` の座標に `mulVec` で作用する**——★**無条件**（第 1233）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

☆第 1228 の**逆向き**。`galMat_apply`（在庫）と `redVec_mulVec`（第 1204）で出る。

★★★これで `σ` の `T_l E` 上の作用（`galTate`）から
`mod l` の行列の性質（冪単性など）を読める。 -/
theorem redVec_galTate (l : ℕ) [Fact l.Prime]
    (W : WeierstrassCurve K) (e : tateModule (W.baseChange L) l ≃+ (Fin 2 → ℤ_[l]))
    (σ : L ≃ₐ[K] L) (x : tateModule (W.baseChange L) l) :
    redVec l (e (galTate W l σ x))
      = ((glRedPadic l (galRep W l e σ) : GL (Fin 2) (ZMod l)) :
          Matrix (Fin 2) (Fin 2) (ZMod l)).mulVec (redVec l (e x)) := by
  have hcoe : ((glRedPadic l (galRep W l e σ) : GL (Fin 2) (ZMod l)) :
      Matrix (Fin 2) (Fin 2) (ZMod l))
      = (galMat W l e σ).map PadicInt.toZMod := by
    rw [coe_glRedPadic, galRep_coe]
  rw [hcoe, ← redVec_mulVec, galMat_apply]

/-! ## ★★★★★★★★★★★★★★★★冪単性と非自明性 -/

/-- ☆`l` 倍は `mod l` で消える。 -/
theorem redVec_nsmul_self (l : ℕ) [Fact l.Prime] (w : Fin 2 → ℤ_[l]) :
    redVec l (l • w) = 0 := by
  rw [redVec_nsmul, ← Nat.cast_smul_eq_nsmul (ZMod l), ZMod.natCast_self, zero_smul]

/-- ★★★★★★★★★★★★★★★★
**`σ` が `mod l` で冪単に作用するなら行列も冪単**——★**無条件**（第 1234）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

☆Tate の描像では `σ` は `μ_l` を固定し `π ↦ ζπ` なので `(σ − 1)² = 0 (mod l)` である。
★★★これが第 1229（冪単は `α` に共役）の仮説を与える。 -/
theorem glRed_unipotent_of_galTate (l : ℕ) [Fact l.Prime]
    (W : WeierstrassCurve K) (e : tateModule (W.baseChange L) l ≃+ (Fin 2 → ℤ_[l]))
    (σ : L ≃ₐ[K] L)
    (h2 : ∀ x : tateModule (W.baseChange L) l, ∃ u : tateModule (W.baseChange L) l,
      galTate W l σ (galTate W l σ x) + x
        = galTate W l σ x + galTate W l σ x + l • u) :
    (((glRedPadic l (galRep W l e σ) : GL (Fin 2) (ZMod l)) :
        Matrix (Fin 2) (Fin 2) (ZMod l)) - 1)
      * (((glRedPadic l (galRep W l e σ) : GL (Fin 2) (ZMod l)) :
        Matrix (Fin 2) (Fin 2) (ZMod l)) - 1) = 0 := by
  set M := ((glRedPadic l (galRep W l e σ) : GL (Fin 2) (ZMod l)) :
    Matrix (Fin 2) (Fin 2) (ZMod l)) with hMdef
  refine Matrix.ext_of_mulVec_single ?_
  intro i
  obtain ⟨w, hw⟩ := exists_redVec_eq l (Pi.single i 1)
  set x : tateModule (W.baseChange L) l := e.symm w with hxdef
  have hex : e x = w := e.apply_symm_apply w
  obtain ⟨u, hu⟩ := h2 x
  have hkey := congrArg (fun y => redVec l (e y)) hu
  simp only [map_add, map_nsmul, redVec_add, redVec_nsmul_self, add_zero,
    redVec_galTate l W e σ, hex] at hkey
  have hkey2 : (M * M).mulVec (redVec l w) + redVec l w
      = M.mulVec (redVec l w) + M.mulVec (redVec l w) := by
    simpa [Matrix.mulVec_mulVec] using hkey
  have hy : (M - 1).mulVec (redVec l w) = M.mulVec (redVec l w) - redVec l w := by
    rw [Matrix.sub_mulVec, Matrix.one_mulVec]
  have hexp : ((M - 1) * (M - 1)).mulVec (redVec l w)
      = (M * M).mulVec (redVec l w) - M.mulVec (redVec l w)
        - (M.mulVec (redVec l w) - redVec l w) := by
    rw [← Matrix.mulVec_mulVec, hy, Matrix.sub_mulVec, Matrix.one_mulVec,
      Matrix.mulVec_sub, Matrix.mulVec_mulVec]
  rw [Matrix.zero_mulVec, ← hw, hexp, eq_sub_of_add_eq hkey2]
  abel

/-! ## ★出典の紐付け(`.src`) -/

def glRed_unipotent_of_galTate.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(σ が mod l で幂単に作用するなら行列も幂単。★無条件)",
    sectionId := "genell-thm-3-8" }

def redVec_galTate.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(glRedPadic (galRep σ) は mod l の座標に mulVec で作用する。★無条件)",
    sectionId := "genell-thm-3-8" }

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
