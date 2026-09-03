/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.AlphaBridgeFull
import ABC3.Found.GenEll.AlphaGlobalTransport2
import ABC3.Meta.Claim

/-!
# 第 1319 ブロック —— **`ζ` が局所体にある版の橋**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★★★★★★★★★これは何か——第 1311 の `ζ` を弱めた版

第 1311 は `ζ₀ ∈ L`（大域）を要求していた。
★第 1312 で「`ζ` の像が `L_v` にあれば足りる」ことが分かったので、
局所の点データ（第 1309）と合わせてその形に組み直す。

☆受け取るのは:

| # | 入力 |
|---|---|
| 1 | `ι`・`σbar`・`hcomm`（`L̄ ↪ M` と制限、第 1167） |
| 2 | `ζ ∈ L̄` が原始 `l` 乗根で、その像が `L_v` に入る |
| 3 | `P₀ ∈ E(L_v)` が位数 `l` |
| 4 | `L_v` の `l`-捩れはたかだか `l` 個 |

★★★ただし `σbar` は第 1309 が作る `σ_M` から来るので、
本ブロックでは `σ_M` も受け取る形にしている。
-/

namespace ABC3.Found.GenEll

open WeierstrassCurve WeierstrassCurve.Affine ABC3.Interface.GaloisRep ABC3.Found.GaloisRep
open ABC3.Meta

open scoped Classical

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★★★★★★★★★★★★★★
**`ζ` が局所体にある版の橋**——★**無条件**（第 1319）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

☆第 1309（局所の点データ）＋第 1310（塔）＋第 1312（`ζ` は局所体で足りる）。 -/
theorem exists_h2_h1_of_localData'
    (L Lv M : Type) [Field L] [CharZero L] [Field Lv] [Field M]
    [Algebra L Lv] [Algebra Lv M] [Algebra L M] [IsScalarTower L Lv M]
    [IsAlgClosed M] [IsGalois Lv M]
    (W : WeierstrassCurve L) [W.IsElliptic] (l : ℕ) [Fact l.Prime]
    (e : tateModule (W.baseChange (AlgebraicClosure L)) l ≃+ (Fin 2 → ℤ_[l]))
    (ι : AlgebraicClosure L →ₐ[L] M)
    (ζ : AlgebraicClosure L) (hζ : IsPrimitiveRoot ζ l)
    (z : Lv) (hz : algebraMap Lv M z = ι ζ)
    (P₀ : (W.baseChange Lv).toAffine.Point) (hP₀ : addOrderOf P₀ = l)
    (hcard : ∀ T : Finset ((W.baseChange Lv).toAffine.Point),
      (∀ p ∈ T, l • p = 0) → T.card ≤ l)
    (σbar : AlgebraicClosure L ≃ₐ[L] AlgebraicClosure L)
    (σM : M ≃ₐ[Lv] M)
    (hcomm : ∀ x : AlgebraicClosure L, ι (σbar x) = (σM.restrictScalars L) (ι x))
    (Pfix' Qmov' : ((W.baseChange Lv).baseChange M).toAffine.Point)
    (hord' : addOrderOf Pfix' = l) (hfix' : galPoint (W.baseChange Lv) σM Pfix' = Pfix')
    (hQ' : l • Qmov' = 0) (hmov' : galPoint (W.baseChange Lv) σM Qmov' ≠ Qmov') :
    (∀ x : tateModule (W.baseChange (AlgebraicClosure L)) l,
        ∃ u : tateModule (W.baseChange (AlgebraicClosure L)) l,
        galTate W l σbar (galTate W l σbar x) + x
          = galTate W l σbar x + galTate W l σbar x + l • u) ∧
      (∃ x : tateModule (W.baseChange (AlgebraicClosure L)) l,
        ∀ u : tateModule (W.baseChange (AlgebraicClosure L)) l,
        galTate W l σbar x ≠ x + l • u) := by
  set φ := pointEquivOfEq (baseChange_baseChange (Lv := Lv) (M := M) W) with hφ
  have hfixζ : σbar ζ = ζ :=
    galEquiv_fixes_of_mem_local L Lv M ι σM σbar hcomm ζ z hz
  refine exists_h2_h1_global_of_local' L Lv M W l e ι σM σbar hcomm ζ hζ hfixζ
    (φ Pfix') ?_ ?_ (φ Qmov') ?_ ?_
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

def exists_h2_h1_of_localData'.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(ζ が局所体にある版の橋。★無条件)",
    sectionId := "genell-thm-3-8" }

def exists_h2_h1_of_localData'.needs : List ProofObligation :=
  [ .citation "[ABC3]" "exists_h2_h1_global_of_local′(第 1312、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_h2_h1_global_of_local'") 1,
    .citation "[ABC3]" "galEquiv_fixes_of_mem_local(第 1312、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.galEquiv_fixes_of_mem_local") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1319）**——第 1311 の `ζ` の仮説を" ++
       "「像が `L_v` に入る」に弱めた版である。" ++
       "☆これで大域の底変換は完全に不要になった。") 2 ]

end ABC3.Found.GenEll
