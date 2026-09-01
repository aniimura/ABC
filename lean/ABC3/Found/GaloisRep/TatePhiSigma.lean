/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.TateGaloisStab
import ABC3.Found.GenEll.PointVariableChange
import ABC3.Meta.Claim

/-!
# 第 1283 ブロック —— **同変性を座標の言葉で**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★★★★★★★これは何か——第 1275 の訂正を受けた同変性

第 1275 で測ったとおり、在庫の `tatePhi_pointMap` は
`σA : K →ₐ[R] K` が恒等射しかないため**非自明な Galois 元に当たらない**。

☆正しい道具は `tatePhi_map`（在庫、`σR : R →+* R` と `σK : K →+* K` の**対**）である。
★本ブロックはそれを **`pointCoords` の言葉**に直す:

> `pointCoords (Φ (σ_* c)) = (σK x, σK y)`  ただし `(x, y) = pointCoords (Φ c)`

☆`Point.map` を経由しないので、曲線がどの環の上にあるかを気にしなくてよい
——`pointCoords` は座標を取るだけだからである。

★★★これが「`tateSigmaAct`（第 1276）＝ 実際の Galois 作用」を言うための形である。
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine ABC3.Found.GenEll ABC3.Meta

variable {R : Type} [CommRing R] [IsDomain R] {I : Ideal R} [IsAdicComplete I R]
  {K : Type} [Field K] [DecidableEq K] [Algebra R K]

set_option maxHeartbeats 800000 in
/-- ★★★★★★★★★★★★★★★★★★★★
**Tate 一意化の同変性——座標の言葉で**——★**無条件**（第 1283）。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

☆`σ_* c` の点の座標は、`c` の点の座標に `σK` を当てたものである。

★★★これが第 1275 で「書き直しが要る」と測った同変性の**正しい形**である。 -/
theorem pointCoords_tatePhi_sigma (S : TateSetup R I K)
    (hΔ : WeierstrassCurve.Δ ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine ≠ 0)
    (σR : R →+* R) (hσI : ∀ x ∈ I, σR x ∈ I) (σK : K →+* K)
    (hcompat : ∀ r : R, σK (algebraMap R K r) = algebraMap R K (σR r))
    (σU : Kˣ →* Kˣ) (hσU : ∀ u : Kˣ, ((σU u : Kˣ) : K) = σK (u : K))
    (hσq : σU S.Q = S.Q) (hσv : ∀ u, vAdd S.v (σU u) = vAdd S.v u)
    (hUinj : Function.Injective σU)
    (c : Kˣ ⧸ Subgroup.zpowers S.Q) :
    pointCoords (tatePhi S hΔ
        (QuotientGroup.map _ _ σU (zpowers_le_comap_self S.Q σU hσq) c))
      = (σK (pointCoords (tatePhi S hΔ c)).1, σK (pointCoords (tatePhi S hΔ c)).2) := by
  have hqq : σR S.q = S.q := tateSetup_q_map S σR σK hcompat σU hσU hσq
  by_cases hc : c = 1
  · subst hc
    rw [map_one, tatePhi_one]
    simp
  · have hc' := tateSetup_quotMap_ne_one S σU hσq hUinj hc
    have hA := tateAOf_map S σR σK hcompat σU hσU hσq hσv c
    have hW := tateWOf_map S σR σK hcompat σU hσU hσq hσv c
    rw [tatePhi_eq S hΔ hc', tatePhi_eq S hΔ hc, tatePtPair, tatePtPair,
      pointCoords_some, pointCoords_some, ← hA, ← hW]
    refine Prod.ext ?_ ?_
    · exact (tateXK_map_fixed σR hσI σK hcompat (tateAOf S c) (tateWOf S c) S.q S.hq hqq
        (tateWOf_mem S c)).symm
    · exact (tateYK_map_fixed σR hσI σK hcompat (tateAOf S c) (tateWOf S c) S.q S.hq hqq
        (tateWOf_mem S c)).symm

/-! ## ★出典の紐付け(`.src`) -/

def pointCoords_tatePhi_sigma.src : Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化の同変性——座標の言葉で。★無条件)",
    sectionId := "genell-def-3-3" }

def pointCoords_tatePhi_sigma.needs : List ProofObligation :=
  [ .citation "[ABC3]" "tateAOf_map・tateWOf_map(在庫、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.tateAOf_map") 1,
    .citation "[ABC3]" "tateXK_map_fixed・tateYK_map_fixed(在庫、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.tateXK_map_fixed") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1283）**——第 1275 で「書き直しが要る」と測った" ++
       "同変性の**正しい形**である。☆`Point.map` を経由しないので、" ++
       "曲線がどの環の上にあるかを気にしなくてよい。" ++
       "★★★これで「`tateSigmaAct`（第 1276）＝ 実際の Galois 作用」が言える。") 2 ]

end ABC3.Found.GaloisRep
