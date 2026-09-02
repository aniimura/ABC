/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.TorsionIntegralGood
import ABC3.Meta.Claim

/-!
# 第 1420 ブロック —— **核の座標が整なら商も整**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★これは何か——第 1074 から `hlu` を外す

第 1074 の `isIntegral_veluQuotientFull_of_addOrderOf_prime` は
`hlu : IsUnit (l : primeSubring p)`（すなわち `p ∤ l`）を仮定していた。
★だが `hlu` を使うのは**最初の 1 段だけ**である
（`pointCoords_mem_of_addOrderOf_prime`——核の座標が整であること）。

☆本ブロックはその 1 段を**仮説として受ける**形に直す。
★★★これで `p ∣ l` でも、別の理由で核の座標が整なら商の整性が出る
——深い核（第 1410-1417）では `tateXpair(a, w, q)`（`a, w ∈ 𝔪`）なので
座標は `𝔪` に入る。
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve IsDedekindDomain NumberField Finset
open ABC3.Found.GenEll ABC3.Meta
open scoped Classical

variable {L : Type} [Field L] [NumberField L]

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★★★★★★**核の座標が `primeSubring p` に入れば Vélu の商も整**（第 1420）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆第 1074 の `hlu`（`p ∤ l`）を**仮説 `hmem` に置き換えた**形である。
★形式群は一度も使っていない。 -/
theorem isIntegral_veluQuotientFull_of_pointCoords_mem (p : HeightOneSpectrum (𝓞 L))
    (E : WeierstrassCurve L) [hE : WeierstrassCurve.IsIntegral (primeSubring p) E]
    {l : ℕ} (hl : l.Prime)
    (Q : E.toAffine.Point) (hQ : addOrderOf Q = l)
    (hmem : ∀ k ∈ (range l).erase 0,
      (pointCoords (k • Q)).1 ∈ primeSubring p ∧
        (pointCoords (k • Q)).2 ∈ primeSubring p) :
    (veluQuotientFull E (((range l).erase 0).image
      (fun k : ℕ => pointCoords (k • Q)))).IsIntegral (primeSubring p) := by
  classical
  obtain ⟨Wi, hWi⟩ := hE.integral
  have hlz : l • Q = 0 := by rw [← hQ]; exact addOrderOf_nsmul_eq_zero Q
  set X : ℕ → primeSubring p := fun i =>
    if h : (pointCoords (i • Q)).1 ∈ primeSubring p then ⟨(pointCoords (i • Q)).1, h⟩
    else 0 with hXdef
  set Y : ℕ → primeSubring p := fun i =>
    if h : (pointCoords (i • Q)).2 ∈ primeSubring p then ⟨(pointCoords (i • Q)).2, h⟩
    else 0 with hYdef
  have hXc : ∀ i ∈ (range l).erase 0,
      algebraMap (primeSubring p) L (X i) = (pointCoords (i • Q)).1 := by
    intro i hi
    simp only [hXdef, dif_pos (hmem i hi).1]
    rfl
  have hYc : ∀ i ∈ (range l).erase 0,
      algebraMap (primeSubring p) L (Y i) = (pointCoords (i • Q)).2 := by
    intro i hi
    simp only [hYdef, dif_pos (hmem i hi).2]
    rfl
  have hP : ∀ i ∈ (range l).erase 0, pointCoords (i • Q)
      = ((algebraMap (primeSubring p) L (X i),
          algebraMap (primeSubring p) L (Y i)) : L × L) := by
    intro i hi
    rw [hXc i hi, hYc i hi]
  -- ☆添字の反転は点の反転(★`l` の偶奇を問わない)
  have hneg : ∀ i ∈ (range l).erase 0,
      pointCoords ((l - i) • Q)
        = ((pointCoords (i • Q)).1,
           (Wi.map (algebraMap (primeSubring p) L)).toAffine.negY
             (pointCoords (i • Q)).1 (pointCoords (i • Q)).2) := by
    intro i hi
    rw [mem_erase, mem_range] at hi
    have hkne : i • Q ≠ 0 :=
      nsmul_ne_zero_of_lt_addOrderOf hQ (by omega) (by omega)
    have hns := nsmul_eq_neg_nsmul_of_addOrderOf hlz (by omega : i ≤ l)
    rw [hns]
    have hEq : E = Wi.map (algebraMap (primeSubring p) L) := hWi
    subst hEq
    exact pointCoords_neg hkne
  -- ★★★★★★★★`w` を作る——`l` が奇なら対で(第 960)、`l = 2` なら 2-捩れで(第 1149)
  obtain ⟨w, hw⟩ : ∃ w : primeSubring p, 2 * w = ∑ i ∈ (range l).erase 0,
      (ABC3.Found.GenEll.veluU Wi (X i) (Y i)
        + 2 * (ABC3.Found.GenEll.veluV2 Wi (X i) (Y i) * X i)) := by
    rcases eq_or_ne l 2 with rfl | hodd
    · -- ★`l = 2`: 唯一の点は 2-捩れなので `veluU = 0`、割り算が要らない
      refine ABC3.Found.GenEll.exists_veluW_two Wi X Y ?_
      have h1 : (1 : ℕ) ∈ (range 2).erase 0 := by decide
      have hn := hneg 1 h1
      have h21 : (2 : ℕ) - 1 = 1 := rfl
      rw [h21] at hn
      have hsnd := congrArg Prod.snd hn
      simp only at hsnd
      apply Subtype.ext
      show algebraMap (primeSubring p) L (Y 1)
        = algebraMap (primeSubring p) L (Wi.toAffine.negY (X 1) (Y 1))
      rw [hYc 1 h1, hsnd, ← hXc 1 h1, ← hYc 1 h1]
      exact (WeierstrassCurve.Affine.map_negY (W' := Wi)
        (algebraMap (primeSubring p) L) (X 1) (Y 1)).symm
    · -- ☆`l = 2m+1`: `i ↔ l−i` の対
      obtain ⟨m, rfl⟩ : ∃ m, l = 2 * m + 1 := hl.odd_of_ne_two hodd
      have hsub : ∀ i ∈ Icc 1 m, (2 * m + 1 - i) ∈ (range (2 * m + 1)).erase 0 := by
        intro i hi
        rw [mem_Icc] at hi
        rw [mem_erase, mem_range]
        omega
      have hin : ∀ i ∈ Icc 1 m, i ∈ (range (2 * m + 1)).erase 0 := by
        intro i hi
        rw [mem_Icc] at hi
        rw [mem_erase, mem_range]
        omega
      refine ABC3.Found.GenEll.exists_veluW_of_inv Wi m X Y ?_ ?_
      · intro i hi
        apply Subtype.ext
        show algebraMap (primeSubring p) L (X (2 * m + 1 - i))
          = algebraMap (primeSubring p) L (X i)
        rw [hXc _ (hsub i hi), hXc i (hin i hi), hneg i (hin i hi)]
      · intro i hi
        apply Subtype.ext
        show algebraMap (primeSubring p) L (Y (2 * m + 1 - i))
          = algebraMap (primeSubring p) L (Wi.toAffine.negY (X i) (Y i))
        rw [hYc _ (hsub i hi), hneg i (hin i hi)]
        rw [← hXc i (hin i hi), ← hYc i (hin i hi)]
        exact (WeierstrassCurve.Affine.map_negY (W' := Wi)
          (algebraMap (primeSubring p) L) (X i) (Y i)).symm
  -- ★単射性
  have hinj : ∀ i ∈ (range l).erase 0, ∀ j ∈ (range l).erase 0,
      ((algebraMap (primeSubring p) L (X i), algebraMap (primeSubring p) L (Y i)) : L × L)
        = ((algebraMap (primeSubring p) L (X j), algebraMap (primeSubring p) L (Y j)) : L × L)
      → i = j := by
    intro i hi j hj hij
    rw [mem_erase, mem_range] at hi hj
    have hne_i : i • Q ≠ 0 := nsmul_ne_zero_of_lt_addOrderOf hQ hi.1 hi.2
    have hne_j : j • Q ≠ 0 := nsmul_ne_zero_of_lt_addOrderOf hQ hj.1 hj.2
    have hEq : i • Q = j • Q := by
      refine pointCoords_injective hne_i hne_j ?_
      rw [hP i (by rw [mem_erase, mem_range]; exact hi),
        hP j (by rw [mem_erase, mem_range]; exact hj)]
      exact hij
    rcases le_total i j with hle | hle
    · have h1 : (j - i) • Q + i • Q = j • Q := by
        rw [← add_nsmul, Nat.sub_add_cancel hle]
      rw [← hEq] at h1
      have h2 : (j - i) • Q = 0 :=
        add_right_cancel (b := i • Q) (by rw [h1, zero_add])
      have h3 : addOrderOf Q ∣ (j - i) := addOrderOf_dvd_of_nsmul_eq_zero h2
      rw [hQ] at h3
      have := Nat.eq_zero_of_dvd_of_lt h3
      omega
    · have h1 : (i - j) • Q + j • Q = i • Q := by
        rw [← add_nsmul, Nat.sub_add_cancel hle]
      rw [hEq] at h1
      have h2 : (i - j) • Q = 0 :=
        add_right_cancel (b := j • Q) (by rw [h1, zero_add])
      have h3 : addOrderOf Q ∣ (i - j) := addOrderOf_dvd_of_nsmul_eq_zero h2
      rw [hQ] at h3
      have := Nat.eq_zero_of_dvd_of_lt h3
      omega
  -- ★組み立て
  have hS : ((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q))
      = ((range l).erase 0).image
          (fun i : ℕ => ((algebraMap (primeSubring p) L (X i),
                          algebraMap (primeSubring p) L (Y i)) : L × L)) :=
    Finset.image_congr (fun i hi => hP i hi)
  refine ⟨ABC3.Found.GenEll.veluCurve Wi
    (∑ i ∈ (range l).erase 0, ABC3.Found.GenEll.veluV2 Wi (X i) (Y i)) w, ?_⟩
  rw [hS, hWi]
  exact ABC3.Found.GenEll.veluQuotientFull_image_eq Wi ((range l).erase 0) X Y hinj
    _ w (two_ne_zero) rfl hw

/-! ## ★出典の紐付け(`.src`) -/

def isIntegral_veluQuotientFull_of_pointCoords_mem.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(核の座標が整なら Vélu の商も整。★p ∤ l 不要)",
    sectionId := "genell-lemma-3-5" }

def isIntegral_veluQuotientFull_of_pointCoords_mem.needs : List ProofObligation :=
  [ .citation "[ABC3]" "veluQuotientFull_image_eq(第 1074 系、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.veluQuotientFull_image_eq") 1,
    .citation "[ABC3]" "exists_veluW_of_inv(第 960、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_veluW_of_inv") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1420）**——第 1074 が `hlu`（`p ∤ l`）を使うのは" ++
       "**最初の 1 段だけ**（核の座標の整性）だと測った。☆その 1 段を仮説にした。" ++
       "★★★深い核（第 1410-1417）では座標が `𝔪` に入るので、" ++
       "`p ∣ l` でもこの補題が使える。") 17 ]

end ABC3.Found.GaloisRep
