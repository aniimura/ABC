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

/-! ## ★★点の族の単射性 -/

/-- ★★**原点でない点は座標で決まる**。 -/
theorem pointCoords_injective {F : Type} [Field F] {W : WeierstrassCurve F}
    {P Q : W.toAffine.Point} (hP : P ≠ 0) (hQ : Q ≠ 0)
    (h : pointCoords P = pointCoords Q) : P = Q := by
  cases P with
  | zero => exact absurd rfl hP
  | some x y hx =>
    cases Q with
    | zero => exact absurd rfl hQ
    | some x' y' hx' =>
      rw [pointCoords_some, pointCoords_some, Prod.mk.injEq] at h
      obtain ⟨h1, h2⟩ := h
      subst h1
      subst h2
      rfl

/-- ★★★★★★★★★★**点の族が単射で原点を避けるなら Vélu の商は一致する**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★`hinj`（座標の単射性）を**点の単射性**に置き換えた形。
☆一意化 `Φ` は加法同型なので単射性はそこから出る。 -/
theorem veluQuotientFull_points_eq' (W : WeierstrassCurve R) (s : Finset ℕ)
    (X Y : ℕ → R) (P : ℕ → (W.map (algebraMap R K)).toAffine.Point)
    (hP : ∀ i ∈ s, pointCoords (P i)
      = ((algebraMap R K (X i), algebraMap R K (Y i)) : K × K))
    (hPne : ∀ i ∈ s, P i ≠ 0)
    (hPinj : ∀ i ∈ s, ∀ j ∈ s, P i = P j → i = j)
    (v w : R) (h2 : (2 : K) ≠ 0)
    (hv : v = ∑ i ∈ s, veluV2 W (X i) (Y i))
    (hw : 2 * w = ∑ i ∈ s, (veluU W (X i) (Y i) + 2 * (veluV2 W (X i) (Y i) * X i))) :
    veluQuotientFull (W.map (algebraMap R K)) (s.image (fun i : ℕ => pointCoords (P i)))
      = (veluCurve W v w).map (algebraMap R K) := by
  refine veluQuotientFull_points_eq W s X Y P hP (fun i hi j hj hij => ?_) v w h2 hv hw
  refine hPinj i hi j hj ?_
  refine pointCoords_injective (hPne i hi) (hPne j hj) ?_
  rw [hP i hi, hP j hj]
  exact hij

/-! ## ★★★★★★★★`Φ` の加法性から `k • Φ(ζ) = Φ(ζᵏ)` -/

section Phi

variable [IsDomain R] [DecidableEq K]

/-- ★★★★★★**`k • Φ(u) = Φ(uᵏ)`**——`Φ` が加法同型だから。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

☆これが「`⟨Q⟩ = μ_l` の各元は `ζᵏ` の行き先」を与える。 -/
theorem nsmul_tatePhi (S : TateSetup R I K)
    (hΔ : ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine.Δ ≠ 0)
    (Φ : Additive (Kˣ ⧸ Subgroup.zpowers S.Q)
      ≃+ ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine.Point)
    (hΦ : ∀ c, Φ (Additive.ofMul c) = tatePhi S hΔ c) (u : Kˣ) (k : ℕ) :
    k • tatePhi S hΔ (QuotientGroup.mk u)
      = tatePhi S hΔ (QuotientGroup.mk (u ^ k)) := by
  have hmk : (QuotientGroup.mk (u ^ k) : Kˣ ⧸ Subgroup.zpowers S.Q)
      = (QuotientGroup.mk u : Kˣ ⧸ Subgroup.zpowers S.Q) ^ k := by
    rw [← QuotientGroup.mk_pow]
  rw [hmk, ← hΦ, ← hΦ, ← map_nsmul]
  congr 1

end Phi

/-! ## ★出典の紐付け(`.src`) -/

def nsmul_tatePhi.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化—k • Φ(u) = Φ(uᵏ)。★無条件)",
    sectionId := "genell-def-3-3" }

def pointCoords_tatePtPair.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 一意化—点の座標は R の元の像である。★無条件)",
    sectionId := "genell-def-3-3" }

def veluQuotientFull_points_eq'.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(点の単射性で述べた形の Vélu の商の一致。★無条件)",
    sectionId := "genell-lemma-3-5" }

def veluQuotientFull_points_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(点の族から Vélu の商の一致へ。★無条件)",
    sectionId := "genell-lemma-3-5" }

end ABC3.Found.GaloisRep
