/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.TateVeluPoints
import ABC3.Found.GenEll.VeluSetStable
import ABC3.Meta.Claim

/-!
# 第 1240 ブロック —— **Vélu の座標集合の個数は `l − 1`**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★これは何か——`IsQuotClassJ` が要る `S.card + 1 = l`

`IsQuotClassJ`（`Found/GenEll/EllModuliObjects.lean`）は
座標集合 `S` に **`S.card + 1 = l`** を要求する。

☆位数 `l` の点 `Q` から作る `S = {pointCoords (k•Q) : 0 < k < l}` は
`k ↦ k•Q` が単射（位数がちょうど `l`）で、`pointCoords` も原点でない点で単射
（`pointCoords_injective`、在庫）なので、**ちょうど `l − 1` 個**である。
-/

namespace ABC3.Found.GenEll

open ABC3.Found.GaloisRep ABC3.Meta

/-- ★★★★★★★★★★★★
**Vélu の座標集合の個数は `l − 1`**——★**無条件**（第 1240）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆`k ↦ k•Q` は `0 < k < l` で単射（位数がちょうど `l`）、
`pointCoords` は原点でない点で単射（在庫）。

★★★これが `IsQuotClassJ` が要求する `S.card + 1 = l` である。 -/
theorem card_image_pointCoords_nsmul {F : Type} [Field F] {W : WeierstrassCurve F}
    [W.IsElliptic] [DecidableEq F] {l : ℕ} (hl : l.Prime)
    {Q : W.toAffine.Point} (hQ : addOrderOf Q = l) :
    (((Finset.range l).erase 0).image (fun k : ℕ => pointCoords (k • Q))).card + 1 = l := by
  classical
  have hne : ∀ k ∈ (Finset.range l).erase 0, k • Q ≠ 0 := by
    intro k hk
    obtain ⟨hk0, hkr⟩ := Finset.mem_erase.1 hk
    exact nsmul_ne_zero_of_lt_addOrderOf hQ hk0 (Finset.mem_range.1 hkr)
  have hinj : ∀ k ∈ (Finset.range l).erase 0, ∀ k' ∈ (Finset.range l).erase 0,
      pointCoords (k • Q) = pointCoords (k' • Q) → k = k' := by
    intro k hk k' hk' h
    have hpt : k • Q = k' • Q := pointCoords_injective (hne k hk) (hne k' hk') h
    obtain ⟨hk0, hkr⟩ := Finset.mem_erase.1 hk
    obtain ⟨hk0', hkr'⟩ := Finset.mem_erase.1 hk'
    have hklt : k < l := Finset.mem_range.1 hkr
    have hklt' : k' < l := Finset.mem_range.1 hkr'
    have key : ∀ a b : ℕ, a ≤ b → b < l → a • Q = b • Q → b - a = 0 := by
      intro a b hab hbl hEq
      have h1 : (b - a) • Q + a • Q = b • Q := by
        rw [← add_nsmul]
        congr 1
        omega
      rw [← hEq] at h1
      have h2 : (b - a) • Q + a • Q = 0 + a • Q := by rw [h1, zero_add]
      have h3 : (b - a) • Q = 0 := add_right_cancel h2
      have h4 : addOrderOf Q ∣ (b - a) := addOrderOf_dvd_of_nsmul_eq_zero h3
      rw [hQ] at h4
      exact Nat.eq_zero_of_dvd_of_lt h4 (by omega)
    rcases le_total k k' with hle | hle
    · have := key k k' hle hklt' hpt
      omega
    · have := key k' k hle hklt hpt.symm
      omega
  rw [Finset.card_image_of_injOn
      (fun k hk k' hk' h => hinj k (Finset.mem_coe.1 hk) k' (Finset.mem_coe.1 hk') h),
    Finset.card_erase_of_mem (Finset.mem_range.2 hl.pos), Finset.card_range]
  have := hl.pos
  omega

/-! ## ★出典の紐付け(`.src`) -/

def card_image_pointCoords_nsmul.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Vélu の座標集合の個数は l − 1。★無条件)",
    sectionId := "genell-lemma-3-5" }

def card_image_pointCoords_nsmul.needs : List ProofObligation :=
  [ .citation "[ABC3]" "pointCoords_injective(証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.pointCoords_injective") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1240）**——`IsQuotClassJ`（`EllModuliObjects`）が" ++
       "要求する `S.card + 1 = l` である。☆`k ↦ k•Q` は `0 < k < l` で単射、" ++
       "`pointCoords` は原点でない点で単射なので、ちょうど `l − 1` 個になる。") 2 ]

end ABC3.Found.GenEll
