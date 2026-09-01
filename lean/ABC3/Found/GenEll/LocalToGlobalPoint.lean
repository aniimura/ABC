/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.Thm38Decomposition
import ABC3.Found.GaloisRep.TorsionTransport
import ABC3.Meta.Claim

/-!
# 第 1287 ブロック —— **局所の 2 条件は大域へ降りる**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★★★★★★★★★これは何か——節点 7 の完成形

第 1271 で「代数閉体の埋め込みで `E[n]` は全単射に移り、幂単性・非自明性は降りる」を取った。
★本ブロックはそれを **`L̄ ↪ M`（`M` は `L_v` の代数閉包）と
`restrictLocalHom`（第 1167）** に当てて、**局所の `σ` の 2 条件を大域の `σ` の 2 条件にする**。

☆埋め込みは `IsAlgClosed.lift`、可換性は `restrictLocalHom_commutes` である。

★★★これで「Tate 一意化で得た局所の `σ`」から
`alpha_mem_map_of_galTate`（第 1237）の `h2`・`h1` が出る道が繋がる
——残るのは第 1270（`E[l]` → `T_l E`）だけである。
-/

namespace ABC3.Found.GenEll

open WeierstrassCurve WeierstrassCurve.Affine ABC3.Found.GaloisRep ABC3.Interface.GaloisRep ABC3.Meta

open scoped Classical

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★★★★★★★★★★★★★★
**局所の幂単性・非自明性は大域へ降りる**——★**無条件**（第 1287）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

☆`M` は `L_v` の代数閉包、`L̄` は `L` の代数閉包。
`IsAlgClosed.lift` で `L̄ ↪ M` を取り、`restrictLocalHom` で `σ_M` を `L̄` に制限する。

★★★これが `AlphaBridge` の節点 7 の完成形である。 -/
theorem exists_galPoint_conditions_global
    (L Lv M : Type) [Field L] [CharZero L] [Field Lv] [Field M]
    [Algebra L Lv] [Algebra Lv M] [Algebra L M] [IsScalarTower L Lv M] [IsAlgClosed M]
    (W : WeierstrassCurve L) [W.IsElliptic] (n : ℕ) (hn : 1 ≤ n)
    (σM : M ≃ₐ[L] M)
    (hunip : ∀ Q : (W.baseChange M).toAffine.Point, n • Q = 0 →
      galPoint W σM (galPoint W σM Q) + Q = galPoint W σM Q + galPoint W σM Q)
    (hne : ∃ Q : (W.baseChange M).toAffine.Point, n • Q = 0 ∧ galPoint W σM Q ≠ Q)
    (σbar : AlgebraicClosure L ≃ₐ[L] AlgebraicClosure L)
    (ι : AlgebraicClosure L →ₐ[L] M)
    (hcomm : ∀ x : AlgebraicClosure L, ι (σbar x) = σM (ι x)) :
    (∀ P : (W.baseChange (AlgebraicClosure L)).toAffine.Point, n • P = 0 →
        galPoint W σbar (galPoint W σbar P) + P = galPoint W σbar P + galPoint W σbar P) ∧
      (∃ P : (W.baseChange (AlgebraicClosure L)).toAffine.Point,
        n • P = 0 ∧ galPoint W σbar P ≠ P) := by
  haveI : CharZero M := charZero_of_injective_algebraMap (algebraMap L M).injective
  haveI : CharZero (AlgebraicClosure L) :=
    charZero_of_injective_algebraMap (algebraMap L (AlgebraicClosure L)).injective
  haveI hEF : (W.baseChange (AlgebraicClosure L)).IsElliptic := by
    unfold WeierstrassCurve.baseChange
    infer_instance
  haveI hEM : (W.baseChange M).IsElliptic := by
    unfold WeierstrassCurve.baseChange
    infer_instance
  have hΔF : (W.baseChange (AlgebraicClosure L)).Δ ≠ 0 :=
    (W.baseChange (AlgebraicClosure L)).isUnit_Δ.ne_zero
  have hΔK : (W.baseChange M).Δ ≠ 0 := (W.baseChange M).isUnit_Δ.ne_zero
  have hcF : ∀ k : ℕ, 1 ≤ k → k ≤ n → ((k : AlgebraicClosure L) ≠ 0) :=
    fun k hk _ => Nat.cast_ne_zero.2 (by omega)
  have hcK : ∀ k : ℕ, 1 ≤ k → k ≤ n → ((k : M) ≠ 0) :=
    fun k hk _ => Nat.cast_ne_zero.2 (by omega)
  refine ⟨?_, ?_⟩
  · exact galPoint_unipotent_of_map W ι σbar σM hcomm n hunip
  · exact exists_galPoint_ne_of_map W hΔF hΔK ι σbar σM hcomm n hn hcF hcK hne

/-! ## ★出典の紐付け(`.src`) -/

def exists_galPoint_conditions_global.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(局所の幂単性・非自明性は大域へ降りる。★無条件)",
    sectionId := "genell-thm-3-8" }

def exists_galPoint_conditions_global.needs : List ProofObligation :=
  [ .citation "[ABC3]" "galPoint_unipotent_of_map・exists_galPoint_ne_of_map(第 1271、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.exists_galPoint_ne_of_map") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1287）**——`AlphaBridge` の節点 7 の完成形である。" ++
       "☆埋め込み `ι` と可換性 `hcomm` は `IsAlgClosed.lift` と" ++
       "`restrictLocalHom_commutes`（第 1167）が与える。" ++
       "★残るのは第 1270（`E[l]` → `T_l E`）だけである。") 2 ]

end ABC3.Found.GenEll
