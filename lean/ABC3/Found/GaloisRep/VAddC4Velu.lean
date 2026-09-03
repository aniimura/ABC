/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.VeluTateC4Val
import ABC3.Meta.Claim

/-!
# 第 1326 ブロック —— **曲線の等式から `c₄` の付値へ**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★これは何か——在庫の `Δ` 版（第 1059）の `c₄` 版

在庫の `vAdd_Delta_of_veluCurve_eq`（第 1059）は

> `(C₀ • X) = veluCurve (tateCurveAt q) v w`（`K` の上で）＋ `v(C₀.u) = 0`
> ⟹ `v(Δ(X))` が計算できる

という形である。★本ブロックは同じ入力から **`v(c₄(X)) = 0`** を出す。

☆`c₄(veluCurve W v w) = c₄(W) + 240v`（在庫）で、それが単元（第 1323）だからである。

★★★これで **Vélu の商は悪い素点で `c₄` が単元**——
第 1322（`c₄` が単元の整モデルは半安定）に渡せば半安定性が出る。
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve ABC3.Found.GenEll ABC3.Meta

variable {R : Type} [CommRing R] [IsDomain R] [IsDiscreteValuationRing R]
  {K : Type} [Field K] [Algebra R K] [IsFractionRing R K]

/-- ★★★★★★★★★★★★★★★★
**曲線の等式から `c₄` の付値へ**——★**無条件**（第 1326）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆在庫の `vAdd_c4_variableChange`（`c₄` の付値の変換則）と
`tateDvrVal_eq_zero_of_isUnit` を繋ぐだけである。 -/
theorem vAdd_c4_of_veluCurve_eq [IsAdicComplete (IsLocalRing.maximalIdeal R) R]
    (q : R) (hq : q ∈ IsLocalRing.maximalIdeal R) (v w : R)
    (hunit : IsUnit ((tateCurveAt q hq).c₄ + 240 * v))
    (C₀ : WeierstrassCurve.VariableChange R) (X : WeierstrassCurve K)
    (hXc4 : X.c₄ ≠ 0)
    (hEq : (C₀.map (algebraMap R K)) • X
      = (veluCurve (tateCurveAt q hq) v w).map (algebraMap R K))
    (hu : vAdd (tateDvrVal R K) ((C₀.map (algebraMap R K)).u) = 0) :
    vAdd (tateDvrVal R K) (Units.mk0 X.c₄ hXc4) = 0 := by
  have hVne : algebraMap R K ((veluCurve (tateCurveAt q hq) v w).c₄) ≠ 0 := by
    rw [veluCurve_c₄]
    exact (hunit.map (algebraMap R K)).ne_zero
  have hc4' : ((C₀.map (algebraMap R K)) • X).c₄ ≠ 0 := by
    rw [hEq, WeierstrassCurve.map_c₄]
    exact hVne
  have hchg := vAdd_c4_variableChange (R := R) (K := K) X hXc4 (C₀.map (algebraMap R K)) hc4'
  rw [hu, mul_zero, sub_zero] at hchg
  have hEqc4 : (Units.mk0 (((C₀.map (algebraMap R K)) • X).c₄) hc4')
      = Units.mk0 (algebraMap R K ((veluCurve (tateCurveAt q hq) v w).c₄)) hVne := by
    refine Units.ext ?_
    rw [Units.val_mk0, Units.val_mk0, hEq, WeierstrassCurve.map_c₄]
  rw [← hchg, hEqc4]
  have hu2 : IsUnit ((veluCurve (tateCurveAt q hq) v w).c₄) := by
    rw [veluCurve_c₄]
    exact hunit
  exact tateDvrVal_eq_zero_of_isUnit _ hu2 hVne

/-! ## ★出典の紐付け(`.src`) -/

def vAdd_c4_of_veluCurve_eq.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(曲線の等式から c₄ の付値へ。★無条件)",
    sectionId := "genell-lemma-3-5" }

def vAdd_c4_of_veluCurve_eq.needs : List ProofObligation :=
  [ .citation "[ABC3]" "vAdd_c4_variableChange(在庫、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.vAdd_c4_variableChange") 1,
    .citation "[ABC3]" "isUnit_c4_velu_tate(第 1323、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.isUnit_c4_velu_tate") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1326）**——在庫の `Δ` 版（第 1059）の `c₄` 版である。" ++
       "☆これで **Vélu の商は悪い素点で `c₄` が単元**になり、" ++
       "第 1322（`c₄` が単元の整モデルは半安定）に渡せば半安定性が出る。") 2 ]

end ABC3.Found.GaloisRep
