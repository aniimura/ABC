/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.GalRepBasis
import ABC3.Found.GenEll.Sl2Padic
import ABC3.Found.GenEll.UnipotentConj
import ABC3.Found.GenEll.GlLift
import ABC3.Found.GenEll.AlphaMatrix
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

open ABC3.Interface.GaloisRep ABC3.Found.GenEll WeierstrassCurve ABC3.Meta

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

/-! ## ★★★★★★★★★★★★★★★★基底を取り替えて共役にする -/

/-- ★★★★★★★★★★★★★★★★
**基底を取り替えれば `galRep` の像は任意の共役に取れる**——★**無条件**（第 1232）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

☆第 1231（基底の取り替えで任意の `C` が実現できる）と
`galRep_basisChange`（在庫）を繋いだもの。

★★★これで「`σ` の `mod l` の行列が `α` に共役」（第 1229）から
「ある基底で `α` そのもの」が出る——共役行列の持ち上げは第 1230。 -/
theorem exists_basis_glRed_conj (W : WeierstrassCurve K) (l : ℕ) [Fact l.Prime]
    (e : tateModule (W.baseChange L) l ≃+ (Fin 2 → ℤ_[l])) (σ : L ≃ₐ[K] L)
    (C : GL (Fin 2) ℤ_[l]) :
    ∃ e₀ : tateModule (W.baseChange L) l ≃+ (Fin 2 → ℤ_[l]),
      glRedPadic l (galRep W l e₀ σ)
        = glRedPadic l C * glRedPadic l (galRep W l e σ) * (glRedPadic l C)⁻¹ := by
  refine ⟨e.trans (mulVecAddEquiv l C), ?_⟩
  have hb : basisChangeGL W l e (e.trans (mulVecAddEquiv l C)) = C :=
    Units.ext (basisChange_realize W l e C)
  rw [galRep_basisChange W l e (e.trans (mulVecAddEquiv l C)) σ, hb, map_mul, map_mul, map_inv]

/-! ## ★★★★★★★★★★★★★★★★★★★★`α` を実現する基底 -/

set_option maxHeartbeats 800000 in
/-- ★★★★★★★★★★★★★★★★★★★★
**`σ` が `mod l` で冪単かつ非自明に作用するなら、ある基底で行列は `α`**
——★**無条件**（第 1236）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

☆繋いだのは 4 段: 冪単性（第 1234）・非自明性（第 1235）→
`α` への共役（第 1229）→ `ℤ_l` への持ち上げ（第 1230）→
基底の取り替え（第 1231-1232）。

★★★これで `α ∈ (galRep の像).map (glRedPadic l)` が出る
——`AlphaBridge` の II 側の**到達点**である。 -/
theorem exists_basis_glRed_eq_alpha (l : ℕ) [Fact l.Prime]
    (W : WeierstrassCurve K) (e : tateModule (W.baseChange L) l ≃+ (Fin 2 → ℤ_[l]))
    (σ : L ≃ₐ[K] L)
    (h2 : ∀ x : tateModule (W.baseChange L) l, ∃ u : tateModule (W.baseChange L) l,
      galTate W l σ (galTate W l σ x) + x
        = galTate W l σ x + galTate W l σ x + l • u)
    (h1 : ∃ x : tateModule (W.baseChange L) l,
      ∀ u : tateModule (W.baseChange L) l, galTate W l σ x ≠ x + l • u) :
    ∃ e₀ : tateModule (W.baseChange L) l ≃+ (Fin 2 → ℤ_[l]),
      ((glRedPadic l (galRep W l e₀ σ) : GL (Fin 2) (ZMod l)) :
        Matrix (Fin 2) (Fin 2) (ZMod l)) = !![1, 1; 0, 1] := by
  have hnil := glRed_unipotent_of_galTate l W e σ h2
  have hne := glRed_ne_one_of_galTate l W e σ h1
  obtain ⟨A, hAdet, hA⟩ := exists_conj_upperOne _ hnil hne
  obtain ⟨B, hBunit, hBmap⟩ := exists_gl_lift l A hAdet
  set Bu : (Matrix (Fin 2) (Fin 2) ℤ_[l])ˣ := Matrix.nonsingInvUnit B hBunit with hBu
  have hBuval : (Bu : Matrix (Fin 2) (Fin 2) ℤ_[l]) = B := rfl
  obtain ⟨e₀, he₀⟩ := exists_basis_glRed_conj W l e σ Bu⁻¹
  refine ⟨e₀, ?_⟩
  have hAu : IsUnit A.det := isUnit_iff_ne_zero.2 hAdet
  have hinvA : A⁻¹ * A = 1 := Matrix.nonsing_inv_mul A hAu
  have hgB : ((glRedPadic l Bu : GL (Fin 2) (ZMod l)) :
      Matrix (Fin 2) (Fin 2) (ZMod l)) = A := by
    rw [coe_glRedPadic, hBuval, hBmap]
  have hcoe := congrArg (fun g : GL (Fin 2) (ZMod l) =>
    (g : Matrix (Fin 2) (Fin 2) (ZMod l))) he₀
  simp only [map_inv, inv_inv, Units.val_mul, Matrix.coe_units_inv, hgB] at hcoe
  rw [hcoe, Matrix.mul_assoc, hA, ← Matrix.mul_assoc, hinvA, Matrix.one_mul]

/-! ## ★出典の紐付け(`.src`) -/

def exists_basis_glRed_eq_alpha.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(σ が mod l で幂単かつ非自明なら、ある基底で行列は α。★無条件)",
    sectionId := "genell-thm-3-8" }

def exists_basis_glRed_eq_alpha.needs : List ProofObligation :=
  [ .citation "[ABC3]" "exists_conj_upperOne(第 1229、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_conj_upperOne") 1,
    .citation "[ABC3]" "exists_gl_lift(第 1230、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_gl_lift") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1236）**——`AlphaBridge` の II 側の**到達点**である。" ++
       "☆繋いだのは 4 段: 幂単性（第 1234）・非自明性（第 1235）→" ++
       "`α` への共役（第 1229）→ `ℤ_l` への持ち上げ（第 1230）→" ++
       "基底の取り替え（第 1231-1232）。") 3 ]

def exists_basis_glRed_conj.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(基底を取り替えれば galRep の像は任意の共役に取れる。★無条件)",
    sectionId := "genell-thm-3-8" }

def exists_basis_glRed_conj.needs : List ProofObligation :=
  [ .citation "[ABC3]" "galRep_basisChange(証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.galRep_basisChange") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1232）**——「`σ` の `mod l` の行列が `α` に" ++
       "共役」（第 1229）から「ある基底で `α` そのもの」を出す段である" ++
       "——共役行列の持ち上げは第 1230。") 2 ]

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
