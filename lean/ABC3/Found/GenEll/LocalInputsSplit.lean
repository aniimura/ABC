/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.TateModelVc
import ABC3.Found.GenEll.LocalInputVc2
import ABC3.Found.GenEll.TateLocalInputs
import ABC3.Found.GaloisRep.TateSetupDvr
import ABC3.Meta.Claim

/-!
# 第 1317 ブロック —— **分裂乗法還元から局所の 2 入力を出す**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★★★★★★★★★★★★★これは何か——残り (a) の**完成形**

☆完備 DVR `R`（分数体 `K`）の上で、`W` が極小かつ分裂乗法還元をもつとき:

| 段 | 在庫 |
|---|---|
| Tate 母数 `q` と変数変換 | `tateParamR_spec`（第 1315 で `K` へ） |
| `TateSetup`・`Φ` | `mkTateSetup`・`dvrTatePhiAddEquiv`（**無条件**） |
| 位数 `l` の点・`l`-捩れの個数 | 第 1314 |
| Tate モデルから `W` へ | 第 1316 |

★★★これで **`W` 自身についての 2 入力**が出る
——第 1311／1312（大域の `h2`・`h1`）が受け取る形である。
-/

namespace ABC3.Found.GenEll

open WeierstrassCurve WeierstrassCurve.Affine ABC3.Found.GaloisRep ABC3.Meta

open scoped Classical

variable {R : Type} [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
  [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
  {K : Type} [Field K] [Algebra R K] [IsFractionRing R K]

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★★★★★★★★★★★★★★
**分裂乗法還元から局所の 2 入力を出す**——★**無条件**（第 1317）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

☆第 1314（Tate モデルの 2 入力）を第 1316（変数変換）で `W` に戻すだけである。 -/
theorem local_inputs_of_split (W : WeierstrassCurve K) [hell : W.IsElliptic]
    [WeierstrassCurve.IsMinimal R W] (h : W.HasSplitMultiplicativeReduction R)
    {l : ℕ} (hl : l.Prime) {ζ : Kˣ} (hζ : IsPrimitiveRoot ((ζ : K)) l)
    (hnd : ¬ ((l : ℤ) ∣ vAdd (mkTateSetup (K := K) (tateParamR W h) (tateParamR_mem W h)
        (tateParamR_ne_zero W h)).v
      (mkTateSetup (K := K) (tateParamR W h) (tateParamR_mem W h)
        (tateParamR_ne_zero W h)).Q)) :
    (∃ P₀ : W.toAffine.Point, addOrderOf P₀ = l) ∧
      (∀ T : Finset (W.toAffine.Point), (∀ p ∈ T, l • p = 0) → T.card ≤ l) := by
  have hqmem := tateParamR_mem W h
  have hqne := tateParamR_ne_zero W h
  have hTΔ : ((tateCurveAt (mkTateSetup (K := K) (tateParamR W h) hqmem hqne).q
      (mkTateSetup (K := K) (tateParamR W h) hqmem hqne).hq).map
        (algebraMap R K)).toAffine.Δ ≠ 0 := by
    show ((tateCurveAt (tateParamR W h) hqmem).map (algebraMap R K)).Δ ≠ 0
    rw [WeierstrassCurve.map_Δ]
    exact (map_ne_zero_iff _ (IsFractionRing.injective R K)).2
      (tateCurveAt_Delta_ne_zero hqmem hqne)
  obtain ⟨hord, hcard⟩ :=
    tate_local_inputs (mkTateSetup (K := K) (tateParamR W h) hqmem hqne) hl hζ
      (dvrTatePhiAddEquiv (K := K) (tateParamR W h) hqmem hqne hTΔ) hnd
  obtain ⟨C', hC'⟩ := exists_vc_tateModel (R := R) W h
  haveI hCell : (C' • W).IsElliptic := by
    rw [hC']
    show (((tateCurveAt (tateParamR W h) hqmem).map (algebraMap R K))).IsElliptic
    exact ⟨isUnit_iff_ne_zero.2 (by
      rw [WeierstrassCurve.map_Δ]
      exact (map_ne_zero_iff _ (IsFractionRing.injective R K)).2
        (tateCurveAt_Delta_ne_zero hqmem hqne))⟩
  constructor
  · obtain ⟨P, hP⟩ := hord
    exact exists_point_order_of_vc' W C'
      ((tateCurveAt (mkTateSetup (K := K) (tateParamR W h) hqmem hqne).q
        (mkTateSetup (K := K) (tateParamR W h) hqmem hqne).hq).map (algebraMap R K))
      hC'.symm P hP
  · exact card_le_of_vc' W C'
      ((tateCurveAt (mkTateSetup (K := K) (tateParamR W h) hqmem hqne).q
        (mkTateSetup (K := K) (tateParamR W h) hqmem hqne).hq).map (algebraMap R K))
      hC'.symm hcard

/-! ## ★出典の紐付け(`.src`) -/

def local_inputs_of_split.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(分裂乗法還元から局所の 2 入力を出す。★無条件)",
    sectionId := "genell-thm-3-8" }

def local_inputs_of_split.needs : List ProofObligation :=
  [ .citation "[ABC3]" "mkTateSetup・dvrTatePhiAddEquiv(在庫、無条件)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.dvrTatePhiAddEquiv") 1,
    .citation "[ABC3]" "tate_local_inputs(第 1314、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.tate_local_inputs") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1317）**——残り (a) の**完成形**である。" ++
       "☆これで `W` 自身についての 2 入力が出るので、" ++
       "第 1311／1312（大域の `h2`・`h1`）に渡せる。") 3 ]

end ABC3.Found.GenEll
