/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.LocalFixedMoved
import ABC3.Found.GenEll.BaseChangeTower
import ABC3.Found.GenEll.AlphaGlobalTransport
import ABC3.Meta.Claim

/-!
# 第 1311 ブロック —— **局所の Tate データから大域の `h2`・`h1` を出す**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★★★★★★★★★★★★★★★★★これは何か——II 側の**橋の完成形**

第 1309（局所で固定点・動く点）＋第 1310（塔の突き合わせ）＋第 1308（輸送）を繋ぐ。

☆受け取るのは**基礎局所体 `L_v` についての 2 つ**と `ζ_l ∈ L` だけである:

| # | 入力 | Tate 曲線での出どころ |
|---|---|---|
| 1 | `ζ₀ ∈ L` が原始 `l` 乗根 | 大域で添加してよい |
| 2 | `P₀ ∈ E(L_v)` が位数 `l` | `μ_l ⊂ E(L_v)`（第 1297） |
| 3 | `L_v` の `l`-捩れはたかだか `l` 個 | `l ∤ v(q)`（第 1304） |

★★★これで `alpha_mem_map_of_galTate`（第 1237）の入力が**局所の Tate データだけ**から出る。
-/

namespace ABC3.Found.GenEll

open WeierstrassCurve WeierstrassCurve.Affine ABC3.Interface.GaloisRep ABC3.Found.GaloisRep
open ABC3.Meta

open scoped Classical

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★
**局所の Tate データから大域の `h2`・`h1` を出す**——★**無条件**（第 1311）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

☆第 1309・1310・1308 を繋いだだけである。 -/
theorem exists_h2_h1_global_of_localData
    (L Lv M : Type) [Field L] [CharZero L] [Field Lv] [Field M]
    [Algebra L Lv] [Algebra Lv M] [Algebra L M] [IsScalarTower L Lv M]
    [IsAlgClosed M] [IsGalois Lv M]
    (W : WeierstrassCurve L) [W.IsElliptic] (l : ℕ) [Fact l.Prime]
    (e : tateModule (W.baseChange (AlgebraicClosure L)) l ≃+ (Fin 2 → ℤ_[l]))
    (ζ₀ : L) (hζ : IsPrimitiveRoot (algebraMap L (AlgebraicClosure L) ζ₀) l)
    (P₀ : (W.baseChange Lv).toAffine.Point) (hP₀ : addOrderOf P₀ = l)
    (hcard : ∀ T : Finset ((W.baseChange Lv).toAffine.Point),
      (∀ p ∈ T, l • p = 0) → T.card ≤ l) :
    ∃ σ : AlgebraicClosure L ≃ₐ[L] AlgebraicClosure L,
      (∀ x : tateModule (W.baseChange (AlgebraicClosure L)) l,
          ∃ u : tateModule (W.baseChange (AlgebraicClosure L)) l,
          galTate W l σ (galTate W l σ x) + x
            = galTate W l σ x + galTate W l σ x + l • u) ∧
        (∃ x : tateModule (W.baseChange (AlgebraicClosure L)) l,
          ∀ u : tateModule (W.baseChange (AlgebraicClosure L)) l,
          galTate W l σ x ≠ x + l • u) := by
  have hl : l.Prime := Fact.out
  haveI : CharZero Lv := charZero_of_injective_algebraMap (algebraMap L Lv).injective
  haveI hEv : (W.baseChange Lv).IsElliptic := by
    unfold WeierstrassCurve.baseChange
    infer_instance
  haveI hEvM : ((W.baseChange Lv).baseChange M).IsElliptic := by
    unfold WeierstrassCurve.baseChange
    infer_instance
  haveI hEvMa : WeierstrassCurve.IsElliptic ((W.baseChange Lv).baseChange M).toAffine :=
    inferInstanceAs (WeierstrassCurve.IsElliptic ((W.baseChange Lv).baseChange M).toAffine)
  obtain ⟨σM, Pfix', Qmov', hord', hfix', hQ', hmov'⟩ :=
    exists_local_fixed_moved (M := M) (W.baseChange Lv) l hl P₀ hP₀ hcard
  set φ := pointEquivOfEq (baseChange_baseChange (Lv := Lv) (M := M) W) with hφ
  refine exists_h2_h1_global_of_local L Lv M W l e ζ₀ hζ σM (φ Pfix') ?_ ?_ (φ Qmov') ?_ ?_
  · rw [hφ]
    have h1 := addOrderOf_injective
      (pointEquivOfEq (baseChange_baseChange (Lv := Lv) (M := M) W)).toAddMonoidHom
      (pointEquivOfEq (baseChange_baseChange (Lv := Lv) (M := M) W)).injective Pfix'
    rw [hord'] at h1
    exact h1
  · rw [hφ, ← galPoint_pointEquivOfEq W σM Pfix', hfix']
  · rw [hφ, ← map_nsmul, hQ', map_zero]
  · rw [hφ, ← galPoint_pointEquivOfEq W σM Qmov']
    intro hcon
    exact hmov' ((pointEquivOfEq (baseChange_baseChange (Lv := Lv) (M := M) W)).injective hcon)

/-! ## ★出典の紐付け(`.src`) -/

def exists_h2_h1_global_of_localData.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(局所の Tate データから大域の h2・h1 を出す。★無条件)",
    sectionId := "genell-thm-3-8" }

def exists_h2_h1_global_of_localData.needs : List ProofObligation :=
  [ .citation "[ABC3]" "exists_local_fixed_moved(第 1309、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_local_fixed_moved") 1,
    .citation "[ABC3]" "exists_h2_h1_global_of_local(第 1308、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_h2_h1_global_of_local") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1311）**——II 側の**橋の完成形**である。" ++
       "☆受け取るのは `ζ_l ∈ L` と、基礎局所体 `L_v` についての 2 つ" ++
       "（位数 `l` の点・`l`-捩れの個数）だけ。" ++
       "★★★どちらも `L_v` の上の Tate 一意化から出る（第 1297・1304）。") 3 ]

end ABC3.Found.GenEll
