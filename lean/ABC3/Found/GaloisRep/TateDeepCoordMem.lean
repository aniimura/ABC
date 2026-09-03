/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.TateVeluDeep
import ABC3.Found.GaloisRep.TateCoordDescend
import ABC3.Meta.Claim

/-!
# 第 1422 ブロック —— **深い核では Tate の座標が `𝔪` に入る**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★これは何か——`p ∣ l` の深い核で整性を出す最後の 1 枚

第 1420（核の座標が整なら商も整）と第 1421（Tate の座標を `p` へ降ろす）を
繋ぐには「Tate の座標が `R` から来る」ことが要る。
★深い核（第 1410-1417）ではそれが**そのまま成り立つ**
——代表 `a([yᵏ])` も `w([yᵏ])` も `𝔪` に入るので、
`tateXpair`・`tateYpair` も `𝔪` に入る（第 1411 の在庫）。

☆第 1420 → 第 1421 → 本ブロックで、**`p ∣ l` の深い核でも
`IsIntegral (primeSubring p) (E/⟨Q⟩)` が出る**（`hlu` を使わずに）。

| 定理 | 内容 |
|---|---|
| `pointCoords_tatePhi_mem_of_deep` | ★★★★座標は `R` から来る（第 1421 の `hTmem`） |
| `pointCoords_tatePhi_mem_maximal_of_deep` | ★★★★★**座標は `𝔪` に入る** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine QuotientGroup ABC3.Found.GenEll Finset
open ABC3.Meta

open scoped Classical

variable {R : Type} [CommRing R] [IsDomain R] {I : Ideal R} [IsAdicComplete I R]
  {K : Type} [Field K] [Algebra R K]

/-- ★★★★**深い核では Tate の座標は `R` から来る**（第 1422）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆これが第 1421 の `pointCoords_mem_primeSubring_of_image_mem` の仮説 `hTmem` である。 -/
theorem pointCoords_tatePhi_mem_of_deep (S : TateSetup R I K)
    (hΔ : ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine.Δ ≠ 0)
    (Φ : Additive (Kˣ ⧸ Subgroup.zpowers S.Q)
      ≃+ ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine.Point)
    (hΦ : ∀ c, Φ (Additive.ofMul c) = tatePhi S hΔ c)
    {l : ℕ} (y : Kˣ)
    (hdeep : ∀ i : ℕ, 0 < i → i < l → ¬ (vAdd S.v S.Q ∣ (i : ℤ) * vAdd S.v y)) :
    ∀ z ∈ ((range l).erase 0).image
        (fun k : ℕ => pointCoords (k • tatePhi S hΔ (QuotientGroup.mk y))),
      (∃ a : R, algebraMap R K a = z.1) ∧ (∃ b : R, algebraMap R K b = z.2) := by
  intro z hz
  obtain ⟨k, hk, rfl⟩ := Finset.mem_image.1 hz
  have hk0 : 0 < k := Nat.pos_of_ne_zero (Finset.mem_erase.1 hk).1
  have hkl : k < l := Finset.mem_range.1 (Finset.mem_erase.1 hk).2
  have hne : (QuotientGroup.mk (y ^ k) : Kˣ ⧸ Subgroup.zpowers S.Q) ≠ 1 :=
    mk_pow_ne_one_of_not_dvd S.v S.Q y k (hdeep k hk0 hkl)
  rw [nsmul_tatePhi S hΔ Φ hΦ y k,
    pointCoords_tatePhi_of_mem S hΔ hne (tateAOf_pow_mem S y k (hdeep k hk0 hkl))]
  exact ⟨⟨_, rfl⟩, ⟨_, rfl⟩⟩

/-- ★★★★★**深い核では Tate の座標は `𝔪` に入る**（第 1422）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★これが「深い核なら `veluV2 ∈ 𝔪`」（第 1412 の `veluV_deep_mem`）の座標側の言い方であり、
同時に整性（第 1420-1421）の入口でもある。 -/
theorem pointCoords_tatePhi_mem_maximal_of_deep (S : TateSetup R I K)
    (hΔ : ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine.Δ ≠ 0)
    (Φ : Additive (Kˣ ⧸ Subgroup.zpowers S.Q)
      ≃+ ((tateCurveAt S.q S.hq).map (algebraMap R K)).toAffine.Point)
    (hΦ : ∀ c, Φ (Additive.ofMul c) = tatePhi S hΔ c)
    {l : ℕ} (y : Kˣ)
    (hdeep : ∀ i : ℕ, 0 < i → i < l → ¬ (vAdd S.v S.Q ∣ (i : ℤ) * vAdd S.v y))
    (k : ℕ) (hk0 : 0 < k) (hkl : k < l) :
    ∃ a ∈ I, ∃ b ∈ I,
      pointCoords (k • tatePhi S hΔ (QuotientGroup.mk y))
        = ((algebraMap R K a, algebraMap R K b) : K × K) := by
  have hne : (QuotientGroup.mk (y ^ k) : Kˣ ⧸ Subgroup.zpowers S.Q) ≠ 1 :=
    mk_pow_ne_one_of_not_dvd S.v S.Q y k (hdeep k hk0 hkl)
  have ha : tateAOf S (QuotientGroup.mk (y ^ k)) ∈ I :=
    tateAOf_pow_mem S y k (hdeep k hk0 hkl)
  have hwm : tateWOf S (QuotientGroup.mk (y ^ k)) ∈ I := tateWOf_mem S _
  refine ⟨_, tateXpair_mem _ _ _ S.hq ha hwm, _, tateYpair_mem _ _ _ S.hq ha hwm, ?_⟩
  rw [nsmul_tatePhi S hΔ Φ hΦ y k, pointCoords_tatePhi_of_mem S hΔ hne ha]

/-! ## ★出典の紐付け(`.src`) -/

def pointCoords_tatePhi_mem_of_deep.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(深い核では Tate の座標は R から来る。★p ∤ l 不要)",
    sectionId := "genell-lemma-3-5" }

def pointCoords_tatePhi_mem_maximal_of_deep.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(深い核では Tate の座標は 𝔪 に入る。★p ∤ l 不要)",
    sectionId := "genell-lemma-3-5" }

def pointCoords_tatePhi_mem_of_deep.needs : List ProofObligation :=
  [ .citation "[ABC3]" "pointCoords_tatePhi_of_mem(第 1411、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.pointCoords_tatePhi_of_mem") 1,
    .citation "[ABC3]" "tateAOf_pow_mem(第 1412、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.tateAOf_pow_mem") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1422）**——第 1420（核の座標が整なら商も整）と" ++
       "第 1421（Tate の座標を `p` へ降ろす）を繋ぐ最後の 1 枚である。" ++
       "☆深い核では代表 `a([yᵏ])`・`w([yᵏ])` がともに `𝔪` に入るので、" ++
       "`tateXpair`・`tateYpair` も `𝔪` に入る。" ++
       "★★★これで **`p ∣ l` の深い核でも商の整性が出る**（`hlu` を使わずに）。" ++
       "☆残るのは `p ∣ l` の `μ_l` 側だけである" ++
       "（`blocked-leaves.json` の `pDivLMeasured2026_09_02`）。") 17 ]

end ABC3.Found.GaloisRep
