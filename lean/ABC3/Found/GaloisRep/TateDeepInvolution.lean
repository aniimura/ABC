/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.TateMuInvolution
import ABC3.Found.GaloisRep.TateVeluDeep
import ABC3.Meta.Claim

/-!
# 第 1415 ブロック —— **深い核でも Vélu の `w` は環の中で作れる**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★これは何か——第 959-961 の深い核版

第 959（`tateXpair_mu_inv`）は「添字 `i ↦ l−i` は点の反転」を `μ_l` の核で示し、
第 961（`exists_veluW_mu`）がそこから `2 ∣ Σ(u_Q + 2v_Q x_Q)` を出していた。

★★機構は**点の側**にある——`(l−i)•P = −(i•P)` と `pointCoords(−Q) = (x, negY x y)`。
☆したがって核が `μ_l` 型かどうかは**まったく関係ない**。
本ブロックは同じ議論を深い核（第 1410-1413 の二者択一の他方）で行う。

| 定理 | 内容 |
|---|---|
| `tateXpair_deep_inv` | ★★★★★★★★★★★★添字 `i ↦ l−i` は点の反転 |
| `exists_veluW_deep` | ★★★★★★★★★★★★★★★★**`w` が `R` の中で取れる** |

★`l = 2m+1`（奇）であることが対分けの本質である——不動点がないから対になる。
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine QuotientGroup ABC3.Found.GenEll Finset

open scoped Classical

variable {R : Type} [CommRing R] [IsDomain R] {I : Ideal R} [IsAdicComplete I R]
  {K : Type} [Field K] [Algebra R K]

variable (S : TateSetup R I K)
  (hΔ : ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine.Δ ≠ 0)

/-- ★★★★★★★★★★★★**深い核でも添字 `i ↦ l−i` は点の反転である**（第 1415）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆`X_{l−i} = X_i` と `Y_{l−i} = negY X_i Y_i`。
★第 959 と同じく**級数の側ではなく点の側**で取る——
`(l−i)•P = −(i•P)`（`l•P = 0` から）と `pointCoords(−Q) = (x, negY x y)`。 -/
theorem tateXpair_deep_inv
    (Φ : Additive (Kˣ ⧸ Subgroup.zpowers S.Q)
      ≃+ ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine.Point)
    (hΦ : ∀ c, Φ (Additive.ofMul c) = tatePhi S hΔ c)
    {l : ℕ} (y : Kˣ)
    (hyl : l • tatePhi S hΔ (QuotientGroup.mk y) = 0)
    (hdeep : ∀ i : ℕ, 0 < i → i < l → ¬ (vAdd S.v S.Q ∣ (i : ℤ) * vAdd S.v y))
    (i : ℕ) (hi : i ∈ (range l).erase 0) :
    tateXpair (tateAOf S (QuotientGroup.mk (y ^ (l - i))))
        (tateWOf S (QuotientGroup.mk (y ^ (l - i)))) S.q S.hq
      = tateXpair (tateAOf S (QuotientGroup.mk (y ^ i)))
        (tateWOf S (QuotientGroup.mk (y ^ i))) S.q S.hq
    ∧ tateYpair (tateAOf S (QuotientGroup.mk (y ^ (l - i))))
        (tateWOf S (QuotientGroup.mk (y ^ (l - i)))) S.q S.hq
      = (tateCurveAt S.q S.hq).toAffine.negY
          (tateXpair (tateAOf S (QuotientGroup.mk (y ^ i)))
            (tateWOf S (QuotientGroup.mk (y ^ i))) S.q S.hq)
          (tateYpair (tateAOf S (QuotientGroup.mk (y ^ i)))
            (tateWOf S (QuotientGroup.mk (y ^ i))) S.q S.hq) := by
  have hi0 : i ≠ 0 := (Finset.mem_erase.1 hi).1
  have hil : i < l := Finset.mem_range.1 (Finset.mem_erase.1 hi).2
  have hli : l - i ∈ (range l).erase 0 := by
    rw [Finset.mem_erase, Finset.mem_range]; omega
  -- ★段 1: `(l−i) • P = −(i • P)`
  have hneg : (l - i) • tatePhi S hΔ (QuotientGroup.mk y)
      = -(i • tatePhi S hΔ (QuotientGroup.mk y)) :=
    nsmul_eq_neg_nsmul_of_addOrderOf hyl hil.le
  -- ★段 2: どの `0 < k < l` でも `[yᵏ] ≠ 1`
  have hcne : ∀ k : ℕ, k ∈ (range l).erase 0 →
      (QuotientGroup.mk (y ^ k) : Kˣ ⧸ Subgroup.zpowers S.Q) ≠ 1 := by
    intro k hk
    exact mk_pow_ne_one_of_not_dvd S.v S.Q y k
      (hdeep k (Nat.pos_of_ne_zero (Finset.mem_erase.1 hk).1)
        (Finset.mem_range.1 (Finset.mem_erase.1 hk).2))
  -- ★段 3: 座標
  have hcoords : ∀ k : ℕ, k ∈ (range l).erase 0 →
      pointCoords (k • tatePhi S hΔ (QuotientGroup.mk y))
        = ((algebraMap R K (tateXpair (tateAOf S (QuotientGroup.mk (y ^ k)))
              (tateWOf S (QuotientGroup.mk (y ^ k))) S.q S.hq),
            algebraMap R K (tateYpair (tateAOf S (QuotientGroup.mk (y ^ k)))
              (tateWOf S (QuotientGroup.mk (y ^ k))) S.q S.hq)) : K × K) := by
    intro k hk
    rw [nsmul_tatePhi S hΔ Φ hΦ y k]
    exact pointCoords_tatePhi_of_mem S hΔ (hcne k hk)
      (tateAOf_pow_mem S y k (hdeep k (Nat.pos_of_ne_zero (Finset.mem_erase.1 hk).1)
        (Finset.mem_range.1 (Finset.mem_erase.1 hk).2)))
  have hine : i • tatePhi S hΔ (QuotientGroup.mk y) ≠ 0 := by
    rw [nsmul_tatePhi S hΔ Φ hΦ y i]
    exact tatePhi_ne_zero S hΔ (hcne i hi)
  -- ★段 4: 座標を比べる
  have hci := hcoords i hi
  have hcli := hcoords (l - i) hli
  rw [hneg, pointCoords_neg hine, hci] at hcli
  have hX := congrArg Prod.fst hcli
  have hY := congrArg Prod.snd hcli
  simp only at hX hY
  refine ⟨S.hinj hX.symm, S.hinj ?_⟩
  rw [← map_negY_algebraMap]
  exact hY.symm

/-- ★★★★★★★★★★★★★★★★**深い核でも Vélu の `w` は環の中で作れる**（第 1415）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆`exists_veluW_of_inv`（第 960、核の形に依らない）に `tateXpair_deep_inv` を当てるだけ。
★`l = 2m+1` の形で受ける——奇数であることが対分けの本質だから。 -/
theorem exists_veluW_deep
    (Φ : Additive (Kˣ ⧸ Subgroup.zpowers S.Q)
      ≃+ ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine.Point)
    (hΦ : ∀ c, Φ (Additive.ofMul c) = tatePhi S hΔ c)
    (m : ℕ) (y : Kˣ)
    (hyl : (2 * m + 1) • tatePhi S hΔ (QuotientGroup.mk y) = 0)
    (hdeep : ∀ i : ℕ, 0 < i → i < 2 * m + 1 →
      ¬ (vAdd S.v S.Q ∣ (i : ℤ) * vAdd S.v y)) :
    ∃ w : R, 2 * w = ∑ i ∈ (range (2 * m + 1)).erase 0,
      (veluU (tateCurveAt S.q S.hq)
          (tateXpair (tateAOf S (QuotientGroup.mk (y ^ i)))
            (tateWOf S (QuotientGroup.mk (y ^ i))) S.q S.hq)
          (tateYpair (tateAOf S (QuotientGroup.mk (y ^ i)))
            (tateWOf S (QuotientGroup.mk (y ^ i))) S.q S.hq)
        + 2 * (veluV2 (tateCurveAt S.q S.hq)
                (tateXpair (tateAOf S (QuotientGroup.mk (y ^ i)))
                  (tateWOf S (QuotientGroup.mk (y ^ i))) S.q S.hq)
                (tateYpair (tateAOf S (QuotientGroup.mk (y ^ i)))
                  (tateWOf S (QuotientGroup.mk (y ^ i))) S.q S.hq)
              * tateXpair (tateAOf S (QuotientGroup.mk (y ^ i)))
                  (tateWOf S (QuotientGroup.mk (y ^ i))) S.q S.hq)) := by
  have hmem : ∀ i ∈ Finset.Icc 1 m, i ∈ (range (2 * m + 1)).erase 0 := by
    intro i hi
    rw [Finset.mem_Icc] at hi
    rw [Finset.mem_erase, Finset.mem_range]
    omega
  refine ABC3.Found.GenEll.exists_veluW_of_inv (tateCurveAt S.q S.hq) m
    (fun i => tateXpair (tateAOf S (QuotientGroup.mk (y ^ i)))
      (tateWOf S (QuotientGroup.mk (y ^ i))) S.q S.hq)
    (fun i => tateYpair (tateAOf S (QuotientGroup.mk (y ^ i)))
      (tateWOf S (QuotientGroup.mk (y ^ i))) S.q S.hq) ?_ ?_
  · intro i hi
    exact (tateXpair_deep_inv S hΔ Φ hΦ y hyl hdeep i (hmem i hi)).1
  · intro i hi
    exact (tateXpair_deep_inv S hΔ Φ hΦ y hyl hdeep i (hmem i hi)).2

/-! ## ★出典の紐付け(`.src`) -/

def tateXpair_deep_inv.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(深い核でも添字 i ↦ l-i は点の反転である。★無条件)",
    sectionId := "genell-lemma-3-5" }

def exists_veluW_deep.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(深い核でも Vélu の w は環の中で作れる。★無条件)",
    sectionId := "genell-lemma-3-5" }

def exists_veluW_deep.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "exists_veluW_of_inv(第 960、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_veluW_of_inv") 1,
    .citation "[ABC3]" "pointCoords_tatePhi_of_mem(第 1411、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.pointCoords_tatePhi_of_mem") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1415）**——第 959-961 の議論は" ++
       "**点の側**（`(l−i)•P = −(i•P)`）にあるので、" ++
       "核が `μ_l` 型かどうかにまったく依らないと測った。" ++
       "☆深い核でそのまま繰り返して `w` を `R` の中で得た。" ++
       "★これで深い側でも `veluCurve (E_q) v w` の `v`・`w` が両方揃う。") 17 ]

end ABC3.Found.GaloisRep
