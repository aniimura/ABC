/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.TateDeepPoints
import ABC3.Found.GaloisRep.TatePhiTwo
import ABC3.Found.GaloisRep.TateQUnique
import ABC3.Found.GaloisRep.TateVeluMu
import ABC3.Found.GaloisRep.TateLinearize
import ABC3.Found.GaloisRep.TateMultRed
import ABC3.Meta.Claim

/-!
# 第 1412 ブロック —— **★★★★★★★★★★★★★★★★深い核による Vélu の商**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★★★★★これは何か——`μ_l` でない核の側

第 890 の `veluQuotientFull_tate_mu` は核が **`μ_l` 型**のときの商を計算した。
★第 1410 の二者択一のもう一方——**深い代表を持つ核**——を本ブロックが扱う。

☆核の生成元 `y` が「どの `0 < i < l` でも `v(Q) ∤ i·v(y)`」を満たすとき:

| 段 | 使うもの |
|---|---|
| `[yⁱ] ≠ 1` | `mk_pow_ne_one_of_not_dvd`（第 1411） |
| 代表 `a([yⁱ]) ∈ 𝔪` | `normRep_vAdd_pos_of_not_dvd` ＋ `tateAOf_mem_of_pos` |
| 座標は `R` の元の像 | `pointCoords_tatePhi_of_mem`（第 1411） |
| 点は別々 | `mk_pow_injOn_of_ne_one`（第 1411）＋ `Φ` の単射性 |
| 商の一致 | `veluQuotientFull_points_eq'`（第 887） |

★★★**`μ_l` 型と違い、商が `E_{q'}` であることは示さない**——示さなくてよい。
☆必要なのは `c₄(veluCurve) = c₄(E_q) + 240v` が**単元**であることだけで、
それは `v ∈ 𝔪` から出る（`veluV_deep_mem`）。

| 定理 | 内容 |
|---|---|
| `tateAOf_pow_mem` | ★★深い類の代表は `𝔪` に入る |
| `veluQuotientFull_tate_deep` | ★★★★★★★★★★★★商 `= veluCurve (E_q) v w` |
| `veluV_deep_mem` | ★★★★★★★★**`v ∈ 𝔪`** |
| `isUnit_c4_add_240_deep` | ★★★★★★★★★★★★★★★★**`c₄ + 240v` は単元** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine QuotientGroup ABC3.Found.GenEll Finset

open scoped Classical

variable {R : Type} [CommRing R] [IsDomain R] {I : Ideal R} [IsAdicComplete I R]
  {K : Type} [Field K] [Algebra R K]

/-! ## ★★深い類の代表は `𝔪` に入る -/

/-- ★★**`v(Q) ∤ i·v(y)` なら `[yⁱ]` の代表は `I` に入る**。

☆`normRep_vAdd_pos_of_not_dvd`（第 1411）で付値が正になり、
`tateAOf_mem_of_pos`（在庫）で `I` の元になる。 -/
theorem tateAOf_pow_mem (S : TateSetup R I K) (y : Kˣ) (i : ℕ)
    (hi : ¬ (vAdd S.v S.Q ∣ (i : ℤ) * vAdd S.v y)) :
    tateAOf S (QuotientGroup.mk (y ^ i)) ∈ I := by
  refine tateAOf_mem_of_pos S _ ?_
  refine normRep_vAdd_pos_of_not_dvd S.v S.Q S.hQ (y ^ i) ?_
  rwa [show vAdd S.v (y ^ i) = (i : ℤ) * vAdd S.v y by rw [← zpow_natCast y i, vAdd_zpow]]

/-! ## ★★★★★★★★★★★★商の形 -/

/-- ★★★★★★★★★★★★**深い核による Vélu の商は `veluCurve (E_q) v w` の底変換**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★第 890 の `veluQuotientFull_tate_mu` の**深い核版**である。
☆`μ_l` 型では代表が `ζⁱ`（単元）だったが、ここでは `a([yⁱ]) ∈ 𝔪` である。 -/
theorem veluQuotientFull_tate_deep (S : TateSetup R I K)
    (hΔ : ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine.Δ ≠ 0)
    (Φ : Additive (Kˣ ⧸ Subgroup.zpowers S.Q)
      ≃+ ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine.Point)
    (hΦ : ∀ c, Φ (Additive.ofMul c) = tatePhi S hΔ c)
    {l : ℕ} (y : Kˣ)
    (hdeep : ∀ i : ℕ, 0 < i → i < l → ¬ (vAdd S.v S.Q ∣ (i : ℤ) * vAdd S.v y))
    (v w : R) (h2 : (2 : K) ≠ 0)
    (hv : v = ∑ i ∈ (range l).erase 0,
      veluV2 (tateCurveAt S.q S.hq)
        (tateXpair (tateAOf S (QuotientGroup.mk (y ^ i)))
          (tateWOf S (QuotientGroup.mk (y ^ i))) S.q S.hq)
        (tateYpair (tateAOf S (QuotientGroup.mk (y ^ i)))
          (tateWOf S (QuotientGroup.mk (y ^ i))) S.q S.hq))
    (hw : 2 * w = ∑ i ∈ (range l).erase 0,
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
                  (tateWOf S (QuotientGroup.mk (y ^ i))) S.q S.hq))) :
    veluQuotientFull ((tateCurveAt S.q S.hq).map (algebraMap R K))
        (((range l).erase 0).image
          (fun k : ℕ => pointCoords (k • tatePhi S hΔ (QuotientGroup.mk y))))
      = (veluCurve (tateCurveAt S.q S.hq) v w).map (algebraMap R K) := by
  have hne : ∀ k : ℕ, 0 < k → k < l →
      (QuotientGroup.mk (y ^ k) : Kˣ ⧸ Subgroup.zpowers S.Q) ≠ 1 := fun k hk hkl =>
    mk_pow_ne_one_of_not_dvd S.v S.Q y k (hdeep k hk hkl)
  refine veluQuotientFull_points_eq' (tateCurveAt S.q S.hq) ((range l).erase 0)
    (fun i => tateXpair (tateAOf S (QuotientGroup.mk (y ^ i)))
      (tateWOf S (QuotientGroup.mk (y ^ i))) S.q S.hq)
    (fun i => tateYpair (tateAOf S (QuotientGroup.mk (y ^ i)))
      (tateWOf S (QuotientGroup.mk (y ^ i))) S.q S.hq)
    (fun k => k • tatePhi S hΔ (QuotientGroup.mk y)) ?_ ?_ ?_ v w h2 hv hw
  · -- ★座標
    intro i hi
    have hi0 : i ≠ 0 := (Finset.mem_erase.1 hi).1
    have hil : i < l := Finset.mem_range.1 (Finset.mem_erase.1 hi).2
    rw [nsmul_tatePhi S hΔ Φ hΦ y i]
    exact pointCoords_tatePhi_of_mem S hΔ (hne i (Nat.pos_of_ne_zero hi0) hil)
      (tateAOf_pow_mem S y i (hdeep i (Nat.pos_of_ne_zero hi0) hil))
  · -- ★原点でない
    intro i hi
    have hi0 : i ≠ 0 := (Finset.mem_erase.1 hi).1
    have hil : i < l := Finset.mem_range.1 (Finset.mem_erase.1 hi).2
    rw [nsmul_tatePhi S hΔ Φ hΦ y i]
    exact tatePhi_ne_zero S hΔ (hne i (Nat.pos_of_ne_zero hi0) hil)
  · -- ★単射
    intro i hi j hj hij
    rw [nsmul_tatePhi S hΔ Φ hΦ y i, nsmul_tatePhi S hΔ Φ hΦ y j, ← hΦ, ← hΦ] at hij
    exact mk_pow_injOn_of_ne_one S.Q y hne
      (Finset.mem_range.1 (Finset.mem_erase.1 hi).2)
      (Finset.mem_range.1 (Finset.mem_erase.1 hj).2) (Φ.injective hij)

/-! ## ★★★★★★★★`v ∈ 𝔪` と `c₄` の単元性 -/

/-- ★★★★★★★★**深い核なら `v = Σ v_Q ∈ 𝔪`**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆各項は `veluV2_tateCurveAt_mem`（第 1411）——`a([yⁱ]) ∈ 𝔪`・`w([yⁱ]) ∈ 𝔪` なので
`tateXpair`・`tateYpair` がともに `𝔪` に入るからである。 -/
theorem veluV_deep_mem (S : TateSetup R I K) {l : ℕ} (y : Kˣ)
    (hdeep : ∀ i : ℕ, 0 < i → i < l → ¬ (vAdd S.v S.Q ∣ (i : ℤ) * vAdd S.v y))
    (v : R)
    (hv : v = ∑ i ∈ (range l).erase 0,
      veluV2 (tateCurveAt S.q S.hq)
        (tateXpair (tateAOf S (QuotientGroup.mk (y ^ i)))
          (tateWOf S (QuotientGroup.mk (y ^ i))) S.q S.hq)
        (tateYpair (tateAOf S (QuotientGroup.mk (y ^ i)))
          (tateWOf S (QuotientGroup.mk (y ^ i))) S.q S.hq)) :
    v ∈ I := by
  rw [hv]
  refine Ideal.sum_mem _ (fun i hi => ?_)
  have hi0 : i ≠ 0 := (Finset.mem_erase.1 hi).1
  have hil : i < l := Finset.mem_range.1 (Finset.mem_erase.1 hi).2
  have ha : tateAOf S (QuotientGroup.mk (y ^ i)) ∈ I :=
    tateAOf_pow_mem S y i (hdeep i (Nat.pos_of_ne_zero hi0) hil)
  have hwm : tateWOf S (QuotientGroup.mk (y ^ i)) ∈ I := tateWOf_mem S _
  exact veluV2_tateCurveAt_mem S.hq (tateXpair_mem _ _ _ S.hq ha hwm)
    (tateYpair_mem _ _ _ S.hq ha hwm)

/-- ★★★★★★★★★★★★★★★★**深い核なら `c₄(veluCurve) = c₄(E_q) + 240v` は単元**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★これが `semistableAt_velu_of_veluCurve_eq_ram`（第 1404）の入力 `hunit` である。
☆`μ_l` 型では `c₄ + 240v = l⁴·c₄(E_{q^l})`（第 1388 の `h4`）が単元性を与えていた。
★深い核では `v ∈ 𝔪` から直接出る——**`l` の情報を一切使わない**。 -/
theorem isUnit_c4_add_240_deep (S : TateSetup R I K) {l : ℕ} (y : Kˣ)
    (hdeep : ∀ i : ℕ, 0 < i → i < l → ¬ (vAdd S.v S.Q ∣ (i : ℤ) * vAdd S.v y))
    (v : R)
    (hv : v = ∑ i ∈ (range l).erase 0,
      veluV2 (tateCurveAt S.q S.hq)
        (tateXpair (tateAOf S (QuotientGroup.mk (y ^ i)))
          (tateWOf S (QuotientGroup.mk (y ^ i))) S.q S.hq)
        (tateYpair (tateAOf S (QuotientGroup.mk (y ^ i)))
          (tateWOf S (QuotientGroup.mk (y ^ i))) S.q S.hq)) :
    IsUnit ((tateCurveAt S.q S.hq).c₄ + 240 * v) :=
  isUnit_add_mem (tateCurveAt_c4_isUnit S.q S.hq)
    (Ideal.mul_mem_left _ _ (veluV_deep_mem S y hdeep v hv))

/-! ## ★出典の紐付け(`.src`) -/

def tateAOf_pow_mem.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(深い類の代表は 𝔪 に入る)",
    sectionId := "genell-def-3-3" }

def veluQuotientFull_tate_deep.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(深い核による Vélu の商は veluCurve (E_q) v w の底変換。★無条件)",
    sectionId := "genell-lemma-3-5" }

def veluV_deep_mem.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(深い核なら v = Σ v_Q ∈ 𝔪)",
    sectionId := "genell-lemma-3-5" }

def isUnit_c4_add_240_deep.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(深い核なら c₄ + 240v は単元。★l の情報を使わない)",
    sectionId := "genell-lemma-3-5" }

def isUnit_c4_add_240_deep.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "veluQuotientFull_points_eq'(第 887、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.veluQuotientFull_points_eq'") 1,
    .citation "[ABC3]" "tateAOf_mem_of_pos(第 292、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.tateAOf_mem_of_pos") 1,
    .citation "[ABC3]" "isUnit_add_mem(第 300 系、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.isUnit_add_mem") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1412）**——第 890 の `veluQuotientFull_tate_mu` の" ++
       "**深い核版**である。" ++
       "☆`μ_l` 型と違い、**商が `E_{q'}` であることは示さない**——示さなくてよい。" ++
       "★必要なのは `c₄(veluCurve) = c₄(E_q) + 240v` が単元であることだけで、" ++
       "それは `v ∈ 𝔪` から出る（`veluV_deep_mem`）。" ++
       "★★★これで `hcop`（`l ∤ jExp p E`）を仮定せずに悪い素点を閉じる道が通った" ++
       "——`hcopDoesNotDescend2026_09_02` の問題を根本から回避する。") 17 ]

end ABC3.Found.GaloisRep
