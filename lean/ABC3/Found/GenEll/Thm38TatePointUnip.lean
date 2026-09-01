/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.Thm38TateUnip
import ABC3.Found.GaloisRep.TateGaloisStab
import ABC3.Meta.Claim

/-!
# 第 1273 ブロック —— **点の側で読んだ幂単性・非自明性**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★★★★★★★★★これは何か——配管 (A)

第 1272 は `Φ`・`τ`・同変性 `hτ` を抽象に受けていた。
★本ブロックは `τ ≔ Point.map σ` として**実際の Tate 曲線の点**に当てる
——同変性は `tatePhi_pointMap`（在庫、`Φ` は `G_K`-加群の同型）が与える。

☆これで「`σ` は Tate 曲線の `E[l]` に幂単かつ非自明に作用する」が
**点の言葉**で出る。★残る配管は 2 本（Tate 曲線と悪い素点の `E` を結ぶ段、
局所から大域へ運ぶ段＝第 1271）である。
-/

namespace ABC3.Found.GenEll

open WeierstrassCurve WeierstrassCurve.Affine ABC3.Found.GaloisRep ABC3.Meta

variable {R : Type} [CommRing R] [IsDomain R] {I : Ideal R} [IsAdicComplete I R]
  {K : Type} [Field K] [DecidableEq K] [Algebra R K]

/-- ★★★★**Tate 曲線の点の上の `σ` の作用**（第 1273）。

☆`Point.map` の値の型は `W⁄K` と書かれるので、
`W.map (algebraMap R K)` との間で `instances` 透明度の食い違いが出る。
★ここで一度名前を付けておくと以降が楽になる。 -/
noncomputable def tatePointMap (S : TateSetup R I K) (σA : K →ₐ[R] K) :
    ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine.Point →+
      ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine.Point :=
  Point.map (S := R) σA

/-- ★★★★★★★★**同変性を `Φ` の言葉に直した形**（第 1273）。 -/
theorem tatePhi_pointMap_mk (S : TateSetup R I K)
    (hΔ : WeierstrassCurve.Δ ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine ≠ 0)
    (Φ : Additive (Kˣ ⧸ Subgroup.zpowers S.Q)
      ≃+ ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine.Point)
    (hΦ : ∀ c, Φ (Additive.ofMul c) = tatePhi S hΔ c)
    (σA : K →ₐ[R] K) (σU : Kˣ →* Kˣ) (hσU : ∀ u : Kˣ, ((σU u : Kˣ) : K) = σA (u : K))
    (hσv : ∀ u, vAdd S.v (σU u) = vAdd S.v u) (x : Kˣ) :
    tatePointMap S σA (Φ (Additive.ofMul (QuotientGroup.mk x)))
      = Φ (Additive.ofMul (QuotientGroup.mk (σU x))) := by
  rw [hΦ, hΦ]
  show Point.map (S := R) σA (tatePhi S hΔ (QuotientGroup.mk x)) = _
  rw [tatePhi_pointMap S hΔ σA σU hσU hσv]
  rfl

set_option maxHeartbeats 800000 in
/-- ★★★★★★★★★★★★★★★★★★★★
**Tate 曲線の点の上で `σ` は幂単**——★**無条件**（第 1273）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

☆第 1272 に `τ ≔ Point.map σ` を代入しただけである。 -/
theorem tate_point_unipotent (S : TateSetup R I K)
    (hΔ : WeierstrassCurve.Δ ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine ≠ 0)
    (Φ : Additive (Kˣ ⧸ Subgroup.zpowers S.Q)
      ≃+ ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine.Point)
    (hΦ : ∀ c, Φ (Additive.ofMul c) = tatePhi S hΔ c)
    {l : ℕ} [NeZero l] {ζ π : Kˣ} (hζ : IsPrimitiveRoot (ζ : K) l) (hπl : π ^ l = S.Q)
    (σA : K →ₐ[R] K) (σU : Kˣ →* Kˣ) (hσU : ∀ u : Kˣ, ((σU u : Kˣ) : K) = σA (u : K))
    (hσv : ∀ u, vAdd S.v (σU u) = vAdd S.v u)
    (hσζ : σU ζ = ζ) (hσπ : σU π = ζ * π)
    (P : ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine.Point) (hP : l • P = 0) :
    tatePointMap S σA (tatePointMap S σA P) + P
      = tatePointMap S σA P + tatePointMap S σA P :=
  tate_unipotent_of_sigma S hζ hπl Φ (tatePointMap S σA) σU
    (tatePhi_pointMap_mk S hΔ Φ hΦ σA σU hσU hσv) hσζ hσπ P hP

/-- ★★★★★★★★★★★★★★★★★★★★
**Tate 曲線の点の上で `σ` は非自明**——★**無条件**（第 1273）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

☆証人は `Φ[π]` である。 -/
theorem tate_point_exists_ne (S : TateSetup R I K)
    (hΔ : WeierstrassCurve.Δ ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine ≠ 0)
    (Φ : Additive (Kˣ ⧸ Subgroup.zpowers S.Q)
      ≃+ ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine.Point)
    (hΦ : ∀ c, Φ (Additive.ofMul c) = tatePhi S hΔ c)
    {l : ℕ} [NeZero l] (hl : 1 < l) {ζ π : Kˣ}
    (hζ : IsPrimitiveRoot (ζ : K) l) (hπl : π ^ l = S.Q)
    (σA : K →ₐ[R] K) (σU : Kˣ →* Kˣ) (hσU : ∀ u : Kˣ, ((σU u : Kˣ) : K) = σA (u : K))
    (hσv : ∀ u, vAdd S.v (σU u) = vAdd S.v u) (hσπ : σU π = ζ * π) :
    ∃ P : ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine.Point,
      l • P = 0 ∧ tatePointMap S σA P ≠ P :=
  tate_exists_ne_of_sigma S hl hζ hπl Φ (tatePointMap S σA) σU
    (tatePhi_pointMap_mk S hΔ Φ hΦ σA σU hσU hσv) hσπ

/-! ## ★出典の紐付け(`.src`) -/

def tatePointMap.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Tate 曲線の点の上の σ の作用)",
    sectionId := "genell-thm-3-8" }

def tatePhi_pointMap_mk.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(同変性を Φ の言葉に直した形。★無条件)",
    sectionId := "genell-thm-3-8" }

def tate_point_unipotent.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Tate 曲線の点の上で σ は幂単。★無条件)",
    sectionId := "genell-thm-3-8" }

def tate_point_exists_ne.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Tate 曲線の点の上で σ は非自明。★無条件)",
    sectionId := "genell-thm-3-8" }

def tate_point_unipotent.needs : List ProofObligation :=
  [ .citation "[ABC3]" "tate_unipotent_of_sigma(第 1272、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.tate_unipotent_of_sigma") 1,
    .citation "[ABC3]" "tatePhi_pointMap(在庫、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.tatePhi_pointMap") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1273）**——第 1272 に `τ ≔ Point.map σ` を代入した形である。" ++
       "☆同変性は `tatePhi_pointMap`（在庫）が与える。" ++
       "★残る配管は 2 本——Tate 曲線と悪い素点の `E` を結ぶ段、" ++
       "局所から大域へ運ぶ段（第 1271）である。") 2 ]

end ABC3.Found.GenEll
