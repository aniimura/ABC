/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.Sl2Padic
import Mathlib.LinearAlgebra.Matrix.NonsingularInverse
import ABC3.Meta.Claim

/-!
# 第 1230 ブロック —— **`GL₂(𝔽_l)` の元は `GL₂(ℤ_l)` に持ち上がる**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★★★★★これは何か——基底の取り替えを `ℤ_l` の側で行うため

第 1229 は「`σ` の `mod l` の行列は `α` に共役」を与えるが、
共役を与える行列 `A` は `𝔽_l` の側にある。

☆基底の取り替え（`galRep_basisChange`、在庫）は `ℤ_l` の側の行列でしか行えないので、
`A` を `GL₂(ℤ_l)` に持ち上げる必要がある。

★成分ごとに代表元を取り（`(a.val : ℕ) : ℤ_l`）、
行列式が単元になることを `ℤ_l` が局所環であることから出す。
-/

namespace ABC3.Found.GenEll

open ABC3.Meta

/-- ★★★★★★★★★★★★
**`GL₂(𝔽_l)` の元は `GL₂(ℤ_l)` に持ち上がる**——★**無条件**（第 1230）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

☆成分ごとに代表元を取り、行列式が単元になることを
`ℤ_l` が局所環であること（`ker toZMod = maximalIdeal`）から出す。

★★★これで第 1229 の共役を `ℤ_l` の側の基底の取り替えとして実行できる。 -/
theorem exists_gl_lift (l : ℕ) [Fact l.Prime] (A : Matrix (Fin 2) (Fin 2) (ZMod l))
    (hA : A.det ≠ 0) :
    ∃ B : Matrix (Fin 2) (Fin 2) ℤ_[l], IsUnit B.det ∧ B.map PadicInt.toZMod = A := by
  haveI : NeZero l := ⟨(Fact.out : Nat.Prime l).ne_zero⟩
  refine ⟨A.map (fun a => ((a.val : ℕ) : ℤ_[l])), ?_, ?_⟩
  · rw [← IsLocalRing.notMem_maximalIdeal, ← PadicInt.ker_toZMod, RingHom.mem_ker]
    intro hz
    refine hA ?_
    have hmap : (A.map (fun a => ((a.val : ℕ) : ℤ_[l]))).map PadicInt.toZMod = A := by
      ext i j
      simp only [Matrix.map_apply]
      rw [map_natCast]
      simp
    have hdet := RingHom.map_det (PadicInt.toZMod (p := l))
      (A.map (fun a => ((a.val : ℕ) : ℤ_[l])))
    rw [hz, RingHom.mapMatrix_apply, hmap] at hdet
    exact hdet.symm
  · ext i j
    simp only [Matrix.map_apply]
    rw [map_natCast]
    simp

/-! ## ★出典の紐付け(`.src`) -/

def exists_gl_lift.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(GL₂(F_l) の元は GL₂(Z_l) に持ち上がる。★無条件)",
    sectionId := "genell-thm-3-8" }

def exists_gl_lift.needs : List ProofObligation :=
  [ .implicitStep
      ("★★★★**2026-09-02（第 1230）**——第 1229 の共役を与える行列は `𝔽_l` の側に" ++
       "あるが、基底の取り替え（`galRep_basisChange`、在庫）は `ℤ_l` の側でしか" ++
       "行えないので持ち上げが要る。☆行列式が単元になることは" ++
       "`ker toZMod = maximalIdeal`（`ℤ_l` は局所環）から出る。") 2 ]

end ABC3.Found.GenEll
