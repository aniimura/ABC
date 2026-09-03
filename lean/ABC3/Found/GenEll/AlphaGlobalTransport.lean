/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.H2H1Package
import ABC3.Found.GaloisRep.TorsionFixedTransport
import ABC3.Found.GenEll.Thm38Decomposition
import ABC3.Meta.Claim

/-!
# 第 1308 ブロック —— **局所の `σ` から大域の `h2`・`h1` を出す**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★★★★★★★★★★★★★これは何か——II 側の**輸送の完成形**

☆局所（`M = L_v` の代数閉包）で得た

* `σ_M` が固定する位数 `l` の点
* `σ_M` が動かす `l`-捩れ点

を、`L̄` の側へ運んで `h2`・`h1` にする。

★道具は 3 つとも既にある:

| 段 | 内容 | 第 |
|---|---|---|
| 1 | `L̄ ↪ M` と制限準同型 | 1167（`restrictLocalHom`） |
| 2 | 固定点・動く点の輸送 | 1302・1271 |
| 3 | `h2`・`h1` の一括 | 1303 |
-/

namespace ABC3.Found.GenEll

open WeierstrassCurve WeierstrassCurve.Affine ABC3.Interface.GaloisRep ABC3.Found.GaloisRep
open ABC3.Meta

open scoped Classical

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★★★★★★★★★★★★★★
**局所の `σ` から大域の `h2`・`h1` を出す**——★**無条件**（第 1308）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

☆`L̄ ↪ M`（`IsAlgClosed.lift`）で運び、第 1303 に渡す。 -/
theorem exists_h2_h1_global_of_local
    (L Lv M : Type) [Field L] [CharZero L] [Field Lv] [Field M]
    [Algebra L Lv] [Algebra Lv M] [Algebra L M] [IsScalarTower L Lv M] [IsAlgClosed M]
    (W : WeierstrassCurve L) [W.IsElliptic] (l : ℕ) [Fact l.Prime]
    (e : tateModule (W.baseChange (AlgebraicClosure L)) l ≃+ (Fin 2 → ℤ_[l]))
    (ζ₀ : L) (hζ : IsPrimitiveRoot (algebraMap L (AlgebraicClosure L) ζ₀) l)
    (σM : M ≃ₐ[Lv] M)
    (Pfix : (W.baseChange M).toAffine.Point) (hPfixOrd : addOrderOf Pfix = l)
    (hfixed : galPoint W (σM.restrictScalars L) Pfix = Pfix)
    (Qmov : (W.baseChange M).toAffine.Point) (hQ : l • Qmov = 0)
    (hmoved : galPoint W (σM.restrictScalars L) Qmov ≠ Qmov) :
    ∃ σ : AlgebraicClosure L ≃ₐ[L] AlgebraicClosure L,
      (∀ x : tateModule (W.baseChange (AlgebraicClosure L)) l,
          ∃ u : tateModule (W.baseChange (AlgebraicClosure L)) l,
          galTate W l σ (galTate W l σ x) + x
            = galTate W l σ x + galTate W l σ x + l • u) ∧
        (∃ x : tateModule (W.baseChange (AlgebraicClosure L)) l,
          ∀ u : tateModule (W.baseChange (AlgebraicClosure L)) l,
          galTate W l σ x ≠ x + l • u) := by
  have hl : l.Prime := Fact.out
  haveI : CharZero M := charZero_of_injective_algebraMap (algebraMap L M).injective
  haveI : CharZero (AlgebraicClosure L) :=
    charZero_of_injective_algebraMap (algebraMap L (AlgebraicClosure L)).injective
  haveI hEF : (W.baseChange (AlgebraicClosure L)).IsElliptic := by
    unfold WeierstrassCurve.baseChange
    infer_instance
  haveI hEM : (W.baseChange M).IsElliptic := by
    unfold WeierstrassCurve.baseChange
    infer_instance
  haveI hEFa : WeierstrassCurve.IsElliptic (W.baseChange (AlgebraicClosure L)).toAffine :=
    inferInstanceAs (WeierstrassCurve.IsElliptic (W.baseChange (AlgebraicClosure L)).toAffine)
  -- `L̄ ↪ M`
  letI ι : AlgebraicClosure L →ₐ[L] M := IsAlgClosed.lift
  letI : Algebra (AlgebraicClosure L) M := ι.toAlgebra
  haveI : IsScalarTower L (AlgebraicClosure L) M :=
    IsScalarTower.of_algebraMap_eq (fun x => (ι.commutes x).symm)
  set σbar := restrictLocalHom L Lv M (AlgebraicClosure L) σM with hσbar
  have hcomm : ∀ x : AlgebraicClosure L, ι (σbar x) = (σM.restrictScalars L) (ι x) := by
    intro x
    exact restrictLocalHom_commutes L Lv M (AlgebraicClosure L) σM x
  have hΔF : (W.baseChange (AlgebraicClosure L)).Δ ≠ 0 :=
    (W.baseChange (AlgebraicClosure L)).isUnit_Δ.ne_zero
  have hΔK : (W.baseChange M).Δ ≠ 0 := (W.baseChange M).isUnit_Δ.ne_zero
  have hcF : ∀ k : ℕ, 1 ≤ k → k ≤ l → ((k : AlgebraicClosure L) ≠ 0) :=
    fun k hk _ => Nat.cast_ne_zero.2 (by omega)
  have hcK : ∀ k : ℕ, 1 ≤ k → k ≤ l → ((k : M) ≠ 0) :=
    fun k hk _ => Nat.cast_ne_zero.2 (by omega)
  -- 固定点を運ぶ
  obtain ⟨P, hPord, hPfix⟩ := exists_galPoint_fixed_of_map W hΔF hΔK ι σbar
    (σM.restrictScalars L) hcomm l hl.pos hcF hcK Pfix hPfixOrd hfixed
  -- 動く点を運ぶ
  obtain ⟨Q, hQ0, hQne⟩ := exists_galPoint_ne_of_map W hΔF hΔK ι σbar
    (σM.restrictScalars L) hcomm l hl.pos hcF hcK ⟨Qmov, hQ, hmoved⟩
  refine ⟨σbar, ?_⟩
  exact galTate_h2_h1_of_fixed_moved W l e σbar (algebraMap L (AlgebraicClosure L) ζ₀) hζ
    (σbar.commutes ζ₀) P hPord hPfix Q hQ0 hQne

/-! ## ★出典の紐付け(`.src`) -/

def exists_h2_h1_global_of_local.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(局所の σ から大域の h2・h1 を出す。★無条件)",
    sectionId := "genell-thm-3-8" }

def exists_h2_h1_global_of_local.needs : List ProofObligation :=
  [ .citation "[ABC3]" "restrictLocalHom・restrictLocalHom_commutes(第 1167、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.restrictLocalHom") 1,
    .citation "[ABC3]" "exists_galPoint_fixed_of_map(第 1302、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.exists_galPoint_fixed_of_map") 1,
    .citation "[ABC3]" "galTate_h2_h1_of_fixed_moved(第 1303、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.galTate_h2_h1_of_fixed_moved") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1308）**——II 側の**輸送の完成形**である。" ++
       "☆局所で得た「固定する位数 `l` の点」と「動かす `l`-捩れ点」を `L̄` へ運ぶだけ。") 3 ]

end ABC3.Found.GenEll
