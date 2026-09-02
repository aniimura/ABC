/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.GalPointCoord
import ABC3.Found.GenEll.GalActPoint
import ABC3.Found.GaloisRep.TateVeluPoints
import ABC3.Meta.Claim

/-!
# 第 1310 ブロック —— **塔での底変換と `galPoint` の突き合わせ**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★これは何か——第 1309 と第 1308 を繋ぐ

第 1309（局所で固定点・動く点を取る）は曲線を `K₀ = L_v` の上で見ている。
第 1308（大域へ運ぶ）は曲線を `L` の上で見ている。
★両者は `(W ⊗ L_v) ⊗ M = W ⊗ M` で繋がる——本ブロックがその突き合わせである。

☆座標は同じ（`pointCoords_pointEquivOfEq`、第 1285）で、
`galPoint` はどちらも座標に `σ` を当てるだけ（第 1298）なので、点の一意性で一致する。
-/

namespace ABC3.Found.GenEll

open WeierstrassCurve WeierstrassCurve.Affine ABC3.Interface.GaloisRep ABC3.Found.GaloisRep
open ABC3.Meta

open scoped Classical

variable {L Lv M : Type} [Field L] [Field Lv] [Field M]
  [Algebra L Lv] [Algebra Lv M] [Algebra L M] [IsScalarTower L Lv M]

/-- ★★★★★★**塔での底変換は 1 段で書ける**——★**無条件**（第 1310）。 -/
theorem baseChange_baseChange (W : WeierstrassCurve L) :
    (W.baseChange Lv).baseChange M = W.baseChange M := by
  show ((W.map (algebraMap L Lv)).map (algebraMap Lv M)) = W.map (algebraMap L M)
  rw [WeierstrassCurve.map_map, ← IsScalarTower.algebraMap_eq]

set_option maxHeartbeats 800000 in
/-- ★★★★★★★★★★★★
**`galPoint` は塔の底変換と可換**——★**無条件**（第 1310）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

☆どちらも座標に `σ` を当てるだけなので、点の一意性で一致する。 -/
theorem galPoint_pointEquivOfEq (W : WeierstrassCurve L) (σM : M ≃ₐ[Lv] M)
    (P : ((W.baseChange Lv).baseChange M).toAffine.Point) :
    pointEquivOfEq (baseChange_baseChange (Lv := Lv) W)
        (galPoint (W.baseChange Lv) σM P)
      = galPoint W (σM.restrictScalars L)
          (pointEquivOfEq (baseChange_baseChange (Lv := Lv) W) P) := by
  by_cases hP : P = 0
  · subst hP
    rw [map_zero, map_zero, map_zero]
  · have hgalP : galPoint (W.baseChange Lv) σM P ≠ 0 := by
      intro hcon
      refine hP ?_
      have h0 : galPoint (W.baseChange Lv) σM P = galPoint (W.baseChange Lv) σM 0 := by
        rw [hcon, map_zero]
      exact Point.map_injective (f := σM.toAlgHom) h0
    have hPe : pointEquivOfEq (baseChange_baseChange (Lv := Lv) W) P ≠ 0 := by
      intro hcon
      refine hP ?_
      have h0 := congrArg (pointEquivOfEq (baseChange_baseChange (Lv := Lv) W)).symm hcon
      simpa using h0
    have hne1 : pointEquivOfEq (baseChange_baseChange (Lv := Lv) W)
        (galPoint (W.baseChange Lv) σM P) ≠ 0 := by
      intro hcon
      refine hgalP ?_
      have h0 := congrArg (pointEquivOfEq (baseChange_baseChange (Lv := Lv) W)).symm hcon
      simpa using h0
    have hne2 : galPoint W (σM.restrictScalars L)
        (pointEquivOfEq (baseChange_baseChange (Lv := Lv) W) P) ≠ 0 := by
      intro hcon
      refine hPe ?_
      have h0 : galPoint W (σM.restrictScalars L)
          (pointEquivOfEq (baseChange_baseChange (Lv := Lv) W) P)
            = galPoint W (σM.restrictScalars L) 0 := by
        rw [hcon, map_zero]
      exact Point.map_injective (f := (σM.restrictScalars L).toAlgHom) h0
    refine pointCoords_injective hne1 hne2 ?_
    rw [pointCoords_pointEquivOfEq, pointCoords_galPoint, pointCoords_galPoint,
      pointCoords_pointEquivOfEq]
    rfl

/-! ## ★出典の紐付け(`.src`) -/

def baseChange_baseChange.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(塔での底変換は 1 段で書ける。★無条件)",
    sectionId := "genell-thm-3-8" }

def galPoint_pointEquivOfEq.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(galPoint は塔の底変換と可換。★無条件)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GenEll
