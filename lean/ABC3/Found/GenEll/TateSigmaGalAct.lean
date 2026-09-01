/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.GalActPoint
import ABC3.Found.GenEll.Thm38SigmaAct
import ABC3.Found.GaloisRep.TatePhiSigma
import ABC3.Found.GaloisRep.TateVeluPoints
import ABC3.Meta.Claim

/-!
# 第 1286 ブロック —— **`tateSigmaAct` は実際の Galois 作用である**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★★★★★★★★★★★★★これは何か——節点 6 が閉じる

第 1276 で `τ ≔ Φ ∘ σU ∘ Φ⁻¹` と定めた `tateSigmaAct` は、
定義から同変だが「本物の Galois 作用か」は別問題だった。

☆本ブロックはそれを**座標の比較 1 本**で片付ける:

| 側 | 座標 | 出典 |
|---|---|---|
| `tateSigmaAct` | `(σK x, σK y)` | 第 1283 |
| `galAct`（点に `σK` を当てる） | `(σK x, σK y)` | 第 1285 |

★点は座標で決まる（`pointCoords_injective`、原点以外）ので、両者は一致する。

★★★これで **`tateSigmaAct` の幂単性・非自明性（第 1276）は
本物の Galois 作用の幂単性・非自明性**になった。
-/

namespace ABC3.Found.GenEll

open WeierstrassCurve WeierstrassCurve.Affine ABC3.Found.GaloisRep ABC3.Meta

open scoped Classical

variable {R : Type} [CommRing R] [IsDomain R] {I : Ideal R} [IsAdicComplete I R]
  {K : Type} [Field K] [Algebra R K]

/-- ★★★★**`rhPoint` が原点に行くのは原点だけ**——★**無条件**（第 1286）。 -/
theorem rhPoint_eq_zero_iff {F K' : Type*} [Field F] [Field K'] (f : F →+* K')
    (W : WeierstrassCurve F) (P : W.toAffine.Point) :
    rhPoint f W P = 0 ↔ P = 0 := by
  cases P with
  | zero => exact ⟨fun _ => rfl, fun _ => rfl⟩
  | some x y h =>
      constructor
      · intro h'
        exact absurd h' (by simp [rhPoint])
      · intro h'
        exact absurd h' (by simp)

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★★★★★★★★★★★★★★
**`tateSigmaAct` は実際の Galois 作用である**——★**無条件**（第 1286）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

☆両者の座標がどちらも `(σK x, σK y)` になるので、点の一意性から一致する。

★★★これが第 1275 で見つけた欠陥の**最終的な埋め合わせ**である。 -/
theorem tateSigmaAct_eq_galAct (S : TateSetup R I K)
    (hΔ : WeierstrassCurve.Δ ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine ≠ 0)
    (Φ : Additive (Kˣ ⧸ Subgroup.zpowers S.Q) ≃+
      ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine.Point)
    (hΦ : ∀ c, Φ (Additive.ofMul c) = tatePhi S hΔ c)
    (σR : R →+* R) (hσI : ∀ x ∈ I, σR x ∈ I) (σK : K →+* K)
    (hcompat : ∀ r : R, σK (algebraMap R K r) = algebraMap R K (σR r))
    (σU : Kˣ →* Kˣ) (hσU : ∀ u : Kˣ, ((σU u : Kˣ) : K) = σK (u : K))
    (hσq : σU S.Q = S.Q) (hσv : ∀ u, vAdd S.v (σU u) = vAdd S.v u)
    (hUinj : Function.Injective σU)
    (hW : ((tateCurveAt S.q S.hq).map (algebraMap R K)).map σK
      = (tateCurveAt S.q S.hq).map (algebraMap R K))
    (P : ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine.Point) :
    tateSigmaAct S Φ σU hσq P = galAct σK _ hW P := by
  obtain ⟨c0, hc0⟩ := Φ.surjective P
  obtain ⟨x, hx⟩ := QuotientGroup.mk_surjective (Additive.toMul c0)
  have hPx : Φ (Additive.ofMul (QuotientGroup.mk x)) = P := by rw [hx]; exact hc0
  have hPt : tatePhi S hΔ (QuotientGroup.mk x) = P := by rw [← hΦ]; exact hPx
  by_cases hc : (QuotientGroup.mk x : Kˣ ⧸ Subgroup.zpowers S.Q) = 1
  · have hP0 : P = 0 := by rw [← hPt, hc, tatePhi_one]
    rw [hP0, map_zero, map_zero]
  · have hc' := tateSetup_quotMap_ne_one S σU hσq hUinj hc
    have hstep : tateSigmaAct S Φ σU hσq P
        = tatePhi S hΔ (QuotientGroup.map _ _ σU (zpowers_le_comap_self S.Q σU hσq)
            (QuotientGroup.mk x)) := by
      rw [← hPx, tateSigmaAct_phi, hΦ]
      rfl
    have hne1 : tateSigmaAct S Φ σU hσq P ≠ 0 := by
      rw [hstep]
      exact tatePhi_ne_zero S hΔ hc'
    have hP0 : P ≠ 0 := by
      rw [← hPt]
      exact tatePhi_ne_zero S hΔ hc
    have hg0 : galAct σK _ hW P ≠ 0 := by
      intro h
      refine hP0 ?_
      have h1 : rhPoint σK _ P = 0 := by
        have h2 := congrArg (pointEquivOfEq hW).symm h
        simpa [galAct, rhPointHom] using h2
      exact (rhPoint_eq_zero_iff σK _ P).1 h1
    refine pointCoords_injective hne1 hg0 ?_
    rw [hstep, pointCoords_galAct,
      pointCoords_tatePhi_sigma S hΔ σR hσI σK hcompat σU hσU hσq hσv hUinj, hPt]

/-! ## ★出典の紐付け(`.src`) -/

def rhPoint_eq_zero_iff.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(rhPoint が原点に行くのは原点だけ。★無条件)",
    sectionId := "genell-thm-3-8" }

def tateSigmaAct_eq_galAct.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(tateSigmaAct は実際の Galois 作用である。★無条件)",
    sectionId := "genell-thm-3-8" }

def tateSigmaAct_eq_galAct.needs : List ProofObligation :=
  [ .citation "[ABC3]" "pointCoords_tatePhi_sigma(第 1283、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.pointCoords_tatePhi_sigma") 1,
    .citation "[ABC3]" "pointCoords_galAct(第 1285、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.pointCoords_galAct") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1286）**——第 1275 で見つけた欠陥の**最終的な埋め合わせ**である。" ++
       "☆両者の座標がどちらも `(σK x, σK y)` になるので、点の一意性から一致する。" ++
       "★★★これで `tateSigmaAct` の幂単性・非自明性（第 1276）は" ++
       "**本物の Galois 作用の幂単性・非自明性**になった。") 3 ]

end ABC3.Found.GenEll
