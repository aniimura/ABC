/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.TateVeluPoints

/-!
# Galois (G6) 第 890 ブロック —— **★★★★★★★★★★★★★★★★★★★★`E_q/μ_l` の組み立て**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★これは何か

第 884-889 で作った 5 本を並べて、

    **`⟨Φ(ζ)⟩` による Vélu の商は `veluCurve (E_q) v w` の底変換である**

を取る。ここの `v`・`w` は `c4_velu_tate`・`c6_velu_tate`（第 853・867）が
使っている和そのものである。

| 段 | 使うもの |
|---|---|
| `k • Φ(ζ) = Φ(ζᵏ)` | `nsmul_tatePhi`（第 888） |
| `Φ(ζᵏ) = tatePtPair ζᵏ (qζ⁻ᵏ)` | `tatePhi_of_pow_eq_one`（第 884） |
| 座標は `R` の元の像 | `pointCoords_tatePtPair`（第 886） |
| 点は別々 | `mk_pow_injOn`（第 889）＋ `Φ` の単射性 |
| 商の一致 | `veluQuotientFull_points_eq′`（第 887） |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine ABC3.Found.GenEll Finset QuotientGroup
open scoped Classical

variable {R : Type} [CommRing R] [IsDomain R] {I : Ideal R} [IsAdicComplete I R]
  {K : Type} [Field K] [Algebra R K]

/-- ★★★★★★★★★★★★★★★★★★★★**`⟨Φ(ζ)⟩` による Vélu の商は
`veluCurve (E_q) v w` の底変換**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★これが「`H` が `μ_l` に対応する」という仮説の**曲線の水準での中身**であり、
`j_velu_tate_mu_map`（第 881）の入力 `hquot` を与える。 -/
theorem veluQuotientFull_tate_mu (S : TateSetup R I K)
    (hΔ : ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine.Δ ≠ 0)
    (Φ : Additive (Kˣ ⧸ Subgroup.zpowers S.Q)
      ≃+ ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine.Point)
    (hΦ : ∀ c, Φ (Additive.ofMul c) = tatePhi S hΔ c)
    {l : ℕ} (hl : 0 < l) (ζ : R) (uζ : Kˣ)
    (hζu : algebraMap R K ζ = (uζ : K)) (hζl : uζ ^ l = 1)
    (hord : ∀ n : ℕ, 0 < n → n < l → uζ ^ n ≠ 1)
    (hu : ∀ i ∈ (range l).erase 0, IsUnit (1 - ζ ^ i))
    (v w : R) (h2 : (2 : K) ≠ 0)
    (hv : v = ∑ i ∈ (range l).erase 0,
      veluV2 (tateCurveAt S.q S.hq)
        (tateXpair (ζ ^ i) (S.q * (ζ ^ i) ^ (l - 1)) S.q S.hq)
        (tateYpair (ζ ^ i) (S.q * (ζ ^ i) ^ (l - 1)) S.q S.hq))
    (hw : 2 * w = ∑ i ∈ (range l).erase 0,
      (veluU (tateCurveAt S.q S.hq)
          (tateXpair (ζ ^ i) (S.q * (ζ ^ i) ^ (l - 1)) S.q S.hq)
          (tateYpair (ζ ^ i) (S.q * (ζ ^ i) ^ (l - 1)) S.q S.hq)
        + 2 * (veluV2 (tateCurveAt S.q S.hq)
                (tateXpair (ζ ^ i) (S.q * (ζ ^ i) ^ (l - 1)) S.q S.hq)
                (tateYpair (ζ ^ i) (S.q * (ζ ^ i) ^ (l - 1)) S.q S.hq)
              * tateXpair (ζ ^ i) (S.q * (ζ ^ i) ^ (l - 1)) S.q S.hq))) :
    veluQuotientFull ((tateCurveAt S.q S.hq).map (algebraMap R K))
        (((range l).erase 0).image
          (fun k : ℕ => pointCoords (k • tatePhi S hΔ (QuotientGroup.mk uζ))))
      = (veluCurve (tateCurveAt S.q S.hq) v w).map (algebraMap R K) := by
  -- ★準備: `ζ^l = 1` は `R` の中でも成り立つ
  have hζlR : ζ ^ l = 1 := pow_eq_one_of_map S ζ uζ hζu hζl
  have hmap : ∀ i : ℕ, algebraMap R K (ζ ^ i) = ((uζ ^ i : Kˣ) : K) := by
    intro i
    rw [map_pow, hζu, Units.val_pow_eq_pow_val]
  have hpowl : ∀ i : ℕ, (uζ ^ i) ^ l = 1 := by
    intro i
    rw [← pow_mul, mul_comm, pow_mul, hζl, one_pow]
  have hawR : ∀ i : ℕ, (ζ ^ i) * (S.q * (ζ ^ i) ^ (l - 1)) = S.q := by
    intro i
    have h1 : (ζ ^ i) * (ζ ^ i) ^ (l - 1) = 1 := by
      have hl' : l - 1 + 1 = l := Nat.succ_pred_eq_of_pos hl
      calc (ζ ^ i) * (ζ ^ i) ^ (l - 1) = (ζ ^ i) ^ (l - 1 + 1) := by ring
        _ = (ζ ^ i) ^ l := by rw [hl']
        _ = (ζ ^ l) ^ i := by rw [← pow_mul, mul_comm, pow_mul]
        _ = 1 := by rw [hζlR, one_pow]
    calc (ζ ^ i) * (S.q * (ζ ^ i) ^ (l - 1))
        = S.q * ((ζ ^ i) * (ζ ^ i) ^ (l - 1)) := by ring
      _ = S.q := by rw [h1, mul_one]
  have hwuR : ∀ i : ℕ, IsUnit (1 - S.q * (ζ ^ i) ^ (l - 1)) := fun i =>
    isUnit_one_sub (I := I) (Ideal.mul_mem_right _ _ S.hq)
  have hneR : ∀ i ∈ (range l).erase 0, algebraMap R K (1 - ζ ^ i) ≠ 0 := fun i hi =>
    ((hu i hi).map (algebraMap R K)).ne_zero
  -- ★類が単位でないこと
  have hmk0 : (QuotientGroup.mk (uζ ^ 0) : Kˣ ⧸ Subgroup.zpowers S.Q) = 1 := by
    rw [pow_zero]
    rfl
  have hcne : ∀ i ∈ (range l).erase 0,
      (QuotientGroup.mk (uζ ^ i) : Kˣ ⧸ Subgroup.zpowers S.Q) ≠ 1 := by
    intro i hi hc
    have hi0 : i ≠ 0 := (Finset.mem_erase.1 hi).1
    have hil : i < l := Finset.mem_range.1 (Finset.mem_erase.1 hi).2
    exact hi0 (mk_pow_injOn S.v S.Q S.hQ hl uζ hζl hord hil hl (by rw [hc, hmk0]))
  -- ★点の族
  refine veluQuotientFull_points_eq' (tateCurveAt S.q S.hq) ((range l).erase 0)
    (fun i => tateXpair (ζ ^ i) (S.q * (ζ ^ i) ^ (l - 1)) S.q S.hq)
    (fun i => tateYpair (ζ ^ i) (S.q * (ζ ^ i) ^ (l - 1)) S.q S.hq)
    (fun k => k • tatePhi S hΔ (QuotientGroup.mk uζ)) ?_ ?_ ?_ v w h2 hv hw
  · -- ★座標
    intro i hi
    rw [nsmul_tatePhi S hΔ Φ hΦ uζ i,
      tatePhi_of_pow_eq_one S hΔ hl (ζ ^ i) (uζ ^ i) (hmap i) (hpowl i) (hcne i hi)
        (hawR i) (hwuR i) (hneR i hi),
      pointCoords_tatePtPair _ _ _ S.hq _ _ _ hΔ (hu i hi)]
  · -- ★原点でない
    intro i hi
    rw [nsmul_tatePhi S hΔ Φ hΦ uζ i]
    exact tatePhi_ne_zero S hΔ (hcne i hi)
  · -- ★単射
    intro i hi j hj hij
    rw [nsmul_tatePhi S hΔ Φ hΦ uζ i, nsmul_tatePhi S hΔ Φ hΦ uζ j, ← hΦ, ← hΦ] at hij
    have hcl : (QuotientGroup.mk (uζ ^ i) : Kˣ ⧸ Subgroup.zpowers S.Q)
        = QuotientGroup.mk (uζ ^ j) := Φ.injective hij
    exact mk_pow_injOn S.v S.Q S.hQ hl uζ hζl hord
      (Finset.mem_range.1 (Finset.mem_erase.1 hi).2)
      (Finset.mem_range.1 (Finset.mem_erase.1 hj).2) hcl

/-! ## ★出典の紐付け(`.src`) -/

def veluQuotientFull_tate_mu.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(⟨Φ(ζ)⟩ による Vélu の商は veluCurve (E_q) v w の底変換。★無条件)",
    sectionId := "genell-lemma-3-5" }

end ABC3.Found.GaloisRep
