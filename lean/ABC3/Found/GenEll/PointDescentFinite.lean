/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.PointDescent
import ABC3.Found.GenEll.VeluDescent
import ABC3.Found.GenEll.PointTransport
import ABC3.Meta.Claim

/-!
# 第 1207 ブロック —— **`L̄` の位数 `l` の点は有限次拡大で有理になる**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★★★★★★★★★これは何か——`L(H)` の道の組み上げ

第 1205 が `HasLCyclicJ` から `L̄` 上の位数 `l` の点 `Q` を作り、
第 1199 が要るのは**有限拡大 `L''` の上で有理な** `Q` である。

本ブロックはその 2 つを繋ぐ:

| 段 | 材料 | 第 |
|---|---|---|
| `{k • Q}` の座標は有限次拡大 `M` に入る | `exists_finite_subextension` | 1195 |
| 座標が `M` に入る点は `M` の点から来る | `exists_rhPoint_eq` | 1206 |
| 曲線の等式に沿った輸送 | `castPoint`・`pointCoords_castPoint` | 1206-1207 |
| 位数は移る | `rhPoint_injective`・`rhPoint_nsmul` | 在庫 |

★★★結論は**`cast` を含まない**——`pointCoords` はどちらの曲線でも
`L̄ × L̄` に値を取るからである。☆これで Vélu の座標集合がそのまま一致する。
-/

namespace ABC3.Found.GenEll

open WeierstrassCurve ABC3.Meta

variable {L : Type} [Field L] [CharZero L]

local notation "Lbar" => AlgebraicClosure L

open scoped Classical in
set_option maxHeartbeats 1600000 in
/-- ★★★★★★★★★★★★★★★★★★★★
**`L̄` の位数 `l` の点は有限次拡大で有理になる**——★**無条件**（第 1207）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆`{k • Q : k < l}` の座標が生成する中間体 `M` は `L` 上有限次であり（第 1195）、
`Q` の座標は `M` に入るので `Q` は `W ⊗ M` の点 `Q'` から来る（第 1206）。
★位数は `rhPoint` が単射な加法写像だから移る。

★★★結論の座標集合の等式は、第 1197
（`veluQuotientFull_descends_to_intermediate`）が受け取る `hT` の形そのものである。

☆**`letI` が要る**——`↥M` は `Subtype` なので `DecidableEq` に
`Subtype.instDecidableEq` が当たってしまい、在庫の補題
（`Classical.propDecidable` で書かれている）と**群構造が別物になる**。 -/
theorem exists_finite_point_descent (W : WeierstrassCurve L) [W.IsElliptic]
    [(W.baseChange Lbar).IsElliptic]
    {l : ℕ} (hl : l.Prime) (Q : (W.baseChange Lbar).toAffine.Point)
    (hQ : addOrderOf Q = l) :
    ∃ M : IntermediateField L Lbar, FiniteDimensional L M ∧
      letI : DecidableEq (M : Type) := fun a b => Classical.propDecidable (a = b)
      ∃ Q' : (W.baseChange M).toAffine.Point, addOrderOf Q' = l ∧
        (((Finset.range l).erase 0).image (fun k : ℕ => pointCoords (k • Q'))).image
            (fun q : M × M => (algebraMap M Lbar q.1, algebraMap M Lbar q.2))
          = ((Finset.range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)) := by
  classical
  obtain ⟨M, hfin, hmem⟩ :=
    exists_finite_subextension (L := L)
      ((Finset.range l).image (fun k : ℕ => pointCoords (k • Q)))
  letI : DecidableEq (M : Type) := fun a b => Classical.propDecidable (a = b)
  have hC : (W.baseChange M).map (algebraMap M Lbar) = W.baseChange Lbar :=
    baseChange_map_intermediate W M
  haveI hellM : (W.baseChange M).IsElliptic := by
    show (W.map (algebraMap L M)).IsElliptic
    infer_instance
  haveI hell : ((W.baseChange M).map (algebraMap M Lbar)).IsElliptic := by
    rw [hC]; infer_instance
  have h1 : (1 : ℕ) ∈ Finset.range l := Finset.mem_range.2 hl.one_lt
  have hQmem : pointCoords Q
      ∈ (Finset.range l).image (fun k : ℕ => pointCoords (k • Q)) :=
    Finset.mem_image.2 ⟨1, h1, by rw [one_nsmul]⟩
  obtain ⟨hx, hy⟩ := hmem _ hQmem
  have hxr : (pointCoords (castPoint hC.symm Q)).1 ∈ Set.range (algebraMap M Lbar) := by
    rw [pointCoords_castPoint]
    exact ⟨⟨_, hx⟩, rfl⟩
  have hyr : (pointCoords (castPoint hC.symm Q)).2 ∈ Set.range (algebraMap M Lbar) := by
    rw [pointCoords_castPoint]
    exact ⟨⟨_, hy⟩, rfl⟩
  obtain ⟨Q', hQ'⟩ := exists_rhPoint_eq (algebraMap M Lbar) (W.baseChange M)
    (castPoint hC.symm Q) hxr hyr
  have hord : addOrderOf (castPoint hC.symm Q) = l := by
    rw [addOrderOf_castPoint]; exact hQ
  have hord' : addOrderOf Q' = l := by
    have hz : l • Q' = 0 := by
      refine rhPoint_injective (algebraMap M Lbar) (W.baseChange M) ?_
      rw [rhPoint_nsmul, hQ', rhPoint_zero]
      exact nsmul_eq_zero_of_addOrderOf hord
    have hne : Q' ≠ 0 := by
      intro h0
      rw [h0, rhPoint_zero] at hQ'
      exact (ne_zero_of_addOrderOf_prime hl hord) hQ'.symm
    have hdvd := addOrderOf_dvd_of_nsmul_eq_zero hz
    rcases hl.eq_one_or_self_of_dvd _ hdvd with h1' | h2'
    · exact absurd (AddMonoid.addOrderOf_eq_one_iff.mp h1') hne
    · exact h2'
  refine ⟨M, hfin, Q', hord', ?_⟩
  rw [← image_pointCoords_rhPoint_nsmul (algebraMap M Lbar) (W.baseChange M) hord']
  refine Finset.image_congr ?_
  intro k _
  dsimp only
  rw [hQ', ← castPoint_nsmul, pointCoords_castPoint]

/-! ## ★出典の紐付け(`.src`) -/

def exists_finite_point_descent.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(L̄ の位数 l の点は有限次拡大で有理になる。★無条件)",
    sectionId := "genell-lemma-3-5" }

def exists_finite_point_descent.needs : List ProofObligation :=
  [ .citation "[ABC3]" "exists_finite_subextension(第 1195、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_finite_subextension") 1,
    .citation "[ABC3]" "exists_rhPoint_eq(第 1206、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_rhPoint_eq") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1207）**——第 1205（`HasLCyclicJ` から `L̄` 上の" ++
       "位数 `l` の点）と第 1199（`Lemma 3.5` は有限拡大の上で確かめれば足りる）を" ++
       "繋ぐ段である。☆結論に `cast` が現れないのは、`pointCoords` が" ++
       "どちらの曲線でも `L̄ × L̄` に値を取るからで、" ++
       "**Vélu の座標集合がそのまま一致する**。") 2 ]

end ABC3.Found.GenEll
