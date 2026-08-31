/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.TateVeluMu
import ABC3.Found.GenEll.VeluPointSet
import ABC3.Found.GenEll.SymmSum

/-!
# 第 959 ブロック —— **★★★★★★★★★★★★★★★★添字 `i ↦ l-i` は点の反転**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★★★★★これは何か——(D3) の (d2)

`hv`・`hw` は `tateXpair (ζ^i) (q(ζ^i)^{l-1})` という**添字の形**で書かれている。
一方 Vélu の `w` を環の中で作るには、**添字 `i ↦ l-i` が点の反転である**
ことが要る（第 957・958）。

★直接に級数をいじると難しい——`tateXpair_symm` は `X(q/u) = X(u)` であり、
欲しいのは `X(1/u) = X(u)` だから、週期性 `u ∼ qu` を級数の側で通す必要がある。

☆**点を経由すれば 1 行である**:

    `(l-i) • P = -(i • P)`（`l • P = 0` だから、第 949）

★座標は `pointCoords_tatePtPair`（第 912）が `tateXpair`・`tateYpair` で与え、
反転の座標は `pointCoords_neg`（第 949）が `(x, negY x y)` で与える。
☆`algebraMap R K` は単射（`S.hinj`）なので `R` の等式に降りる。

| 定理 | 内容 |
|---|---|
| `map_negY_algebraMap` | ★`negY` は底変換と可換 |
| `tateMu_pointCoords` | ★添字 `i` の点の座標 |
| `tateXpair_mu_inv` | ★★★★★★★★★★★★★★★★**`X_{l-i} = X_i`・`Y_{l-i} = negY`** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine ABC3.Found.GenEll Finset QuotientGroup
open scoped Classical

variable {R : Type} [CommRing R] [IsDomain R] {I : Ideal R} [IsAdicComplete I R]
  {K : Type} [Field K] [Algebra R K]

/-- ★**`negY` は底変換と可換する**。 -/
theorem map_negY_algebraMap (W : WeierstrassCurve R) (x y : R) :
    (W.map (algebraMap R K)).toAffine.negY (algebraMap R K x) (algebraMap R K y)
      = algebraMap R K (W.toAffine.negY x y) := by
  show -(algebraMap R K y) - (W.map (algebraMap R K)).a₁ * (algebraMap R K x)
      - (W.map (algebraMap R K)).a₃
    = algebraMap R K (-y - W.a₁ * x - W.a₃)
  show -(algebraMap R K y) - (algebraMap R K W.a₁) * (algebraMap R K x)
      - (algebraMap R K W.a₃)
    = algebraMap R K (-y - W.a₁ * x - W.a₃)
  rw [map_sub, map_sub, map_neg, map_mul]

section MuInv

variable (S : TateSetup R I K)
  (hΔ : ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine.Δ ≠ 0)

/-- ★**添字 `i` の `μ_l` の点の座標**——第 890 の中身を取り出したもの。 -/
theorem tateMu_pointCoords
    (Φ : Additive (Kˣ ⧸ Subgroup.zpowers S.Q)
      ≃+ ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine.Point)
    (hΦ : ∀ c, Φ (Additive.ofMul c) = tatePhi S hΔ c)
    {l : ℕ} (hl : 0 < l) (ζ : R) (uζ : Kˣ)
    (hζu : algebraMap R K ζ = (uζ : K)) (hζl : uζ ^ l = 1)
    (hord : ∀ n : ℕ, 0 < n → n < l → uζ ^ n ≠ 1)
    (hu : ∀ i ∈ (range l).erase 0, IsUnit (1 - ζ ^ i))
    (i : ℕ) (hi : i ∈ (range l).erase 0) :
    pointCoords (i • tatePhi S hΔ (QuotientGroup.mk uζ))
      = ((algebraMap R K (tateXpair (ζ ^ i) (S.q * (ζ ^ i) ^ (l - 1)) S.q S.hq),
          algebraMap R K (tateYpair (ζ ^ i) (S.q * (ζ ^ i) ^ (l - 1)) S.q S.hq)) : K × K) := by
  have hζlR : ζ ^ l = 1 := pow_eq_one_of_map S ζ uζ hζu hζl
  have hmap : algebraMap R K (ζ ^ i) = ((uζ ^ i : Kˣ) : K) := by
    rw [map_pow, hζu, Units.val_pow_eq_pow_val]
  have hpowl : (uζ ^ i) ^ l = 1 := by
    rw [← pow_mul, mul_comm, pow_mul, hζl, one_pow]
  have hawR : (ζ ^ i) * (S.q * (ζ ^ i) ^ (l - 1)) = S.q := by
    have h1 : (ζ ^ i) * (ζ ^ i) ^ (l - 1) = 1 := by
      have hl' : l - 1 + 1 = l := Nat.succ_pred_eq_of_pos hl
      calc (ζ ^ i) * (ζ ^ i) ^ (l - 1) = (ζ ^ i) ^ (l - 1 + 1) := by ring
        _ = (ζ ^ i) ^ l := by rw [hl']
        _ = (ζ ^ l) ^ i := by rw [← pow_mul, mul_comm, pow_mul]
        _ = 1 := by rw [hζlR, one_pow]
    calc (ζ ^ i) * (S.q * (ζ ^ i) ^ (l - 1))
        = S.q * ((ζ ^ i) * (ζ ^ i) ^ (l - 1)) := by ring
      _ = S.q := by rw [h1, mul_one]
  have hwuR : IsUnit (1 - S.q * (ζ ^ i) ^ (l - 1)) :=
    isUnit_one_sub (I := I) (Ideal.mul_mem_right _ _ S.hq)
  have hneR : algebraMap R K (1 - ζ ^ i) ≠ 0 := ((hu i hi).map (algebraMap R K)).ne_zero
  have hmk0 : (QuotientGroup.mk (uζ ^ 0) : Kˣ ⧸ Subgroup.zpowers S.Q) = 1 := by
    rw [pow_zero]; rfl
  have hcne : (QuotientGroup.mk (uζ ^ i) : Kˣ ⧸ Subgroup.zpowers S.Q) ≠ 1 := by
    intro hc
    have hi0 : i ≠ 0 := (Finset.mem_erase.1 hi).1
    have hil : i < l := Finset.mem_range.1 (Finset.mem_erase.1 hi).2
    exact hi0 (mk_pow_injOn S.v S.Q S.hQ hl uζ hζl hord hil hl (by rw [hc, hmk0]))
  rw [nsmul_tatePhi S hΔ Φ hΦ uζ i,
    tatePhi_of_pow_eq_one S hΔ hl (ζ ^ i) (uζ ^ i) hmap hpowl hcne hawR hwuR hneR,
    pointCoords_tatePtPair _ _ _ S.hq _ _ _ hΔ (hu i hi)]

/-- ★★★★★★★★★★★★★★★★**添字 `i ↦ l-i` は点の反転である**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★★**2026-09-01（第 959）**——これが (D3) の (d2) である。
☆`X_{l-i} = X_i` と `Y_{l-i} = negY X_i Y_i`。
★級数の側ではなく**点の側**で取るのが鍵である。 -/
theorem tateXpair_mu_inv
    (Φ : Additive (Kˣ ⧸ Subgroup.zpowers S.Q)
      ≃+ ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine.Point)
    (hΦ : ∀ c, Φ (Additive.ofMul c) = tatePhi S hΔ c)
    {l : ℕ} (hl : 0 < l) (ζ : R) (uζ : Kˣ)
    (hζu : algebraMap R K ζ = (uζ : K)) (hζl : uζ ^ l = 1)
    (hord : ∀ n : ℕ, 0 < n → n < l → uζ ^ n ≠ 1)
    (hu : ∀ i ∈ (range l).erase 0, IsUnit (1 - ζ ^ i))
    (i : ℕ) (hi : i ∈ (range l).erase 0) :
    tateXpair (ζ ^ (l - i)) (S.q * (ζ ^ (l - i)) ^ (l - 1)) S.q S.hq
        = tateXpair (ζ ^ i) (S.q * (ζ ^ i) ^ (l - 1)) S.q S.hq
      ∧ tateYpair (ζ ^ (l - i)) (S.q * (ζ ^ (l - i)) ^ (l - 1)) S.q S.hq
        = (tateCurveAt S.q S.hq).toAffine.negY
            (tateXpair (ζ ^ i) (S.q * (ζ ^ i) ^ (l - 1)) S.q S.hq)
            (tateYpair (ζ ^ i) (S.q * (ζ ^ i) ^ (l - 1)) S.q S.hq) := by
  have hi0 : i ≠ 0 := (Finset.mem_erase.1 hi).1
  have hil : i < l := Finset.mem_range.1 (Finset.mem_erase.1 hi).2
  have hli : l - i ∈ (range l).erase 0 := by
    rw [Finset.mem_erase, Finset.mem_range]
    omega
  -- ★段 1: `l • P = 0`
  have hl0 : l • tatePhi S hΔ (QuotientGroup.mk uζ) = 0 := by
    rw [nsmul_tatePhi S hΔ Φ hΦ uζ l, hζl]
    show tatePhi S hΔ (QuotientGroup.mk (1 : Kˣ)) = 0
    rw [show (QuotientGroup.mk (1 : Kˣ) : Kˣ ⧸ Subgroup.zpowers S.Q) = 1 from rfl,
      tatePhi_one]
  -- ★段 2: `(l-i) • P = -(i • P)`
  have hneg : (l - i) • tatePhi S hΔ (QuotientGroup.mk uζ)
      = -(i • tatePhi S hΔ (QuotientGroup.mk uζ)) :=
    nsmul_eq_neg_nsmul_of_addOrderOf hl0 hil.le
  -- ★段 3: `i • P ≠ 0`
  have hmk0 : (QuotientGroup.mk (uζ ^ 0) : Kˣ ⧸ Subgroup.zpowers S.Q) = 1 := by
    rw [pow_zero]; rfl
  have hcne : (QuotientGroup.mk (uζ ^ i) : Kˣ ⧸ Subgroup.zpowers S.Q) ≠ 1 := by
    intro hc
    exact hi0 (mk_pow_injOn S.v S.Q S.hQ hl uζ hζl hord hil hl (by rw [hc, hmk0]))
  have hine : i • tatePhi S hΔ (QuotientGroup.mk uζ) ≠ 0 := by
    rw [nsmul_tatePhi S hΔ Φ hΦ uζ i]
    exact tatePhi_ne_zero S hΔ hcne
  -- ★段 4: 座標を比べる
  have hci := tateMu_pointCoords S hΔ Φ hΦ hl ζ uζ hζu hζl hord hu i hi
  have hcli := tateMu_pointCoords S hΔ Φ hΦ hl ζ uζ hζu hζl hord hu (l - i) hli
  rw [hneg, pointCoords_neg hine, hci] at hcli
  -- `hcli : (X i, negY (X i) (Y i)) = (X (l-i), Y (l-i))` の向きに注意
  have hX := congrArg Prod.fst hcli
  have hY := congrArg Prod.snd hcli
  simp only at hX hY
  refine ⟨S.hinj hX.symm, S.hinj ?_⟩
  rw [← map_negY_algebraMap]
  exact hY.symm

/-- ★★★★★★★★★★★★★★★★★★★★**`μ_l` の Vélu の `w` は環の中で作れる**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★★**2026-09-01（第 961）**——これが (D3) の (d4) である。
`tateParam_quot_velu_of_torsion`（第 948）の `hvw` の `hw` の部分がこれで出る。

☆`exists_veluW_of_inv`（第 960）に `tateXpair_mu_inv`（第 959）を当てるだけである。
★`l` を `2m+1` の形で受ける——奇数であることが対分けの本質だから。 -/
theorem exists_veluW_mu
    (Φ : Additive (Kˣ ⧸ Subgroup.zpowers S.Q)
      ≃+ ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine.Point)
    (hΦ : ∀ c, Φ (Additive.ofMul c) = tatePhi S hΔ c)
    (m : ℕ) (ζ : R) (uζ : Kˣ)
    (hζu : algebraMap R K ζ = (uζ : K)) (hζl : uζ ^ (2 * m + 1) = 1)
    (hord : ∀ n : ℕ, 0 < n → n < 2 * m + 1 → uζ ^ n ≠ 1)
    (hu : ∀ i ∈ (range (2 * m + 1)).erase 0, IsUnit (1 - ζ ^ i)) :
    ∃ w : R, 2 * w = ∑ i ∈ (range (2 * m + 1)).erase 0,
      (veluU (tateCurveAt S.q S.hq)
          (tateXpair (ζ ^ i) (S.q * (ζ ^ i) ^ (2 * m + 1 - 1)) S.q S.hq)
          (tateYpair (ζ ^ i) (S.q * (ζ ^ i) ^ (2 * m + 1 - 1)) S.q S.hq)
        + 2 * (veluV2 (tateCurveAt S.q S.hq)
                (tateXpair (ζ ^ i) (S.q * (ζ ^ i) ^ (2 * m + 1 - 1)) S.q S.hq)
                (tateYpair (ζ ^ i) (S.q * (ζ ^ i) ^ (2 * m + 1 - 1)) S.q S.hq)
              * tateXpair (ζ ^ i) (S.q * (ζ ^ i) ^ (2 * m + 1 - 1)) S.q S.hq)) := by
  have hmem : ∀ i ∈ Finset.Icc 1 m, i ∈ (range (2 * m + 1)).erase 0 := by
    intro i hi
    rw [Finset.mem_Icc] at hi
    rw [Finset.mem_erase, Finset.mem_range]
    omega
  refine ABC3.Found.GenEll.exists_veluW_of_inv (tateCurveAt S.q S.hq) m
    (fun i => tateXpair (ζ ^ i) (S.q * (ζ ^ i) ^ (2 * m + 1 - 1)) S.q S.hq)
    (fun i => tateYpair (ζ ^ i) (S.q * (ζ ^ i) ^ (2 * m + 1 - 1)) S.q S.hq) ?_ ?_
  · intro i hi
    exact (tateXpair_mu_inv S hΔ Φ hΦ (Nat.succ_pos _) ζ uζ hζu hζl hord hu i (hmem i hi)).1
  · intro i hi
    exact (tateXpair_mu_inv S hΔ Φ hΦ (Nat.succ_pos _) ζ uζ hζu hζl hord hu i (hmem i hi)).2

end MuInv

/-! ## ★出典の紐付け(`.src`) -/

def tateMu_pointCoords.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(添字 i の μ_l の点の座標。★無条件)",
    sectionId := "genell-lemma-3-5" }

def tateXpair_mu_inv.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(添字 i ↦ l-i は点の反転である。★無条件)",
    sectionId := "genell-lemma-3-5" }


def exists_veluW_mu.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(μ_l の Vélu の w は環の中で作れる。★無条件)",
    sectionId := "genell-lemma-3-5" }

end ABC3.Found.GaloisRep
