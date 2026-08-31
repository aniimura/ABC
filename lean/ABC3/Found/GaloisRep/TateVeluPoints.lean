/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.VeluImage
import ABC3.Found.GenEll.PointVariableChange
import ABC3.Found.GaloisRep.TateMuPoint

/-!
# Galois (G6) 第 886 ブロック —— **★★★★★★★★★★Tate の点の Vélu の商**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★これは何か

第 884 で「`Φ(ζⁱ)` は `tatePtPair ζⁱ (q(ζⁱ)^{l-1})`」を、
第 885 で「座標の対の像で取った Vélu の商は `veluCurve W v w` の底変換」を取った。
★本ブロックはその両者を繋ぐ釘である:

    `pointCoords (tatePtPair a w q hq …) = (φ (tateXpair a w q hq), φ (tateYpair a w q hq))`

☆`tatePtPair` は `Point.some _ _ (nonsingular_tateK …)` なので座標は `tateXK`・`tateYK`であり、
`1 - a` が単元なら `tateXK_eq`・`tateYK_eq` がそれを `R` の元の像に直す。

| 定理 | 内容 |
|---|---|
| `pointCoords_tatePtPair` | ★★★★Tate の点の座標 |
| `veluQuotientFull_points_eq` | ★★★★★★★★★★**点の族から商の一致へ** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine ABC3.Found.GenEll Finset
open scoped Classical

variable {R : Type} [CommRing R] {I : Ideal R} [IsAdicComplete I R]
  {K : Type} [Field K] [Algebra R K]

/-- ★★★★**Tate の点の座標**——`1 - a` が単元なら `R` の元の像である。 -/
theorem pointCoords_tatePtPair (a w q : R) (hq : q ∈ I) (haw : a * w = q)
    (hw : IsUnit (1 - w)) (hne : algebraMap R K (1 - a) ≠ 0)
    (hΔ : ((tateCurveAt q hq).map (algebraMap R K)).toAffine.Δ ≠ 0)
    (ha : IsUnit (1 - a)) :
    pointCoords (tatePtPair a w q hq haw hw hne hΔ)
      = ((algebraMap R K (tateXpair a w q hq), algebraMap R K (tateYpair a w q hq)) : K × K) := by
  rw [tatePtPair, pointCoords_some, tateXK_eq a w q hq ha, tateYK_eq a w q hq ha]

/-- ★★★★★★★★★★**点の族から Vélu の商の一致へ**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★座標が `R` の元の像であるような点の族を与えれば、
その Vélu の商は `veluCurve W v w` の底変換である。 -/
theorem veluQuotientFull_points_eq (W : WeierstrassCurve R) (s : Finset ℕ)
    (X Y : ℕ → R) (P : ℕ → (W.map (algebraMap R K)).toAffine.Point)
    (hP : ∀ i ∈ s, pointCoords (P i)
      = ((algebraMap R K (X i), algebraMap R K (Y i)) : K × K))
    (hinj : ∀ i ∈ s, ∀ j ∈ s,
      ((algebraMap R K (X i), algebraMap R K (Y i)) : K × K)
        = ((algebraMap R K (X j), algebraMap R K (Y j)) : K × K) → i = j)
    (v w : R) (h2 : (2 : K) ≠ 0)
    (hv : v = ∑ i ∈ s, veluV2 W (X i) (Y i))
    (hw : 2 * w = ∑ i ∈ s, (veluU W (X i) (Y i) + 2 * (veluV2 W (X i) (Y i) * X i))) :
    veluQuotientFull (W.map (algebraMap R K)) (s.image (fun i : ℕ => pointCoords (P i)))
      = (veluCurve W v w).map (algebraMap R K) := by
  have himg : s.image (fun i : ℕ => pointCoords (P i))
      = s.image (fun i : ℕ => ((algebraMap R K (X i), algebraMap R K (Y i)) : K × K)) :=
    Finset.image_congr (fun i hi => hP i hi)
  rw [himg]
  exact veluQuotientFull_image_eq W s X Y hinj v w h2 hv hw

/-! ## ★出典の紐付け(`.src`) -/

def pointCoords_tatePtPair.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化—点の座標は R の元の像である。★無条件)",
    sectionId := "genell-def-3-3" }

def veluQuotientFull_points_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(点の族から Vélu の商の一致へ。★無条件)",
    sectionId := "genell-lemma-3-5" }

end ABC3.Found.GaloisRep
