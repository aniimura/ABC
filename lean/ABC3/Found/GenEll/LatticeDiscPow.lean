/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.VeluNorm
import ABC3.Meta.Claim

/-!
# 第 1401 ブロック —— **格子曲線の上の判別式の恒等式**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★これは何か

    latticeDisc Λ ^ l = latticeDisc Λ′ · ( ∏_{w ∈ T∖{0}} ℘′_Λ(w) )^4     （`l = 2m+1` 奇）

★★★これが `Skeleton/GenEll/VeluDiscIdentity.lean` の
`disc_pow_eq_veluQuot_mul_lattice` の中身である。

## ★★★★★★★★組み立て

| 部品 | 出どころ |
|---|---|
| `Δ = 16·((e₁−e₂)(e₁−e₃)(e₂−e₃))²` | 第 1397 |
| `∏_{i≠j}(℘(v_j+w) − e_i) = −D`（`w` に依らない） | 第 1398 |
| **同種のノルム** `∏_{w∈T}(℘(z+w) − e_i) = c_i(℘_{Λ′}(z) − e′_i)` | 第 1400 |

☆帳簿は 3 行である。`v_i` を `Λ` の半周期として

1. 同種のノルムを `z = v_j` で使う: `(e_j − e_i)·R_{ij} = c_i·(e′_j − e′_i)`（`i ≠ j`）
2. 6 本の順序対すべてを掛ける: `(−D)·∏R = (c₁c₂c₃)²·(−D′)`
3. `∏_{i≠j}R_{ij} = ∏_{w∈T∖0}(−D) = (−D)^{l−1} = D^{l−1}`（第 1398、`l` 奇）

★したがって `D^l = (c₁c₂c₃)²·D′`。☆あとは `N² = 4^{l−1}c₁c₂c₃`（`℘′² = 4∏(℘−e_i)`）と
`Δ = 16D`・`Δ′ = 16D′` を入れれば `Δ^l = Δ′·N⁴` が出る。

★★☆`Λ` の半周期は `Λ′` の 2-捩れでもある（`Λ′/Λ` の位数が奇だから）
——`hodd2` がその仮定である。
-/

namespace ABC3.Found.GenEll

open PeriodPair Finset ABC3.Meta

open scoped Classical

/-! ## ★★★★2-捩れ点の値は相異なる（一般形） -/

/-- ☆`u ± u′ ∉ M` なら `℘_M(u) ≠ ℘_M(u′)`。 -/
theorem weierstrassP_two_torsion_ne (M : PeriodPair) {u u' : ℂ}
    (hu : u ∉ M.lattice) (hu' : u' ∉ M.lattice)
    (hd : u - u' ∉ M.lattice) (hs : u + u' ∉ M.lattice) :
    M.weierstrassP u ≠ M.weierstrassP u' := by
  intro heq
  rcases sub_or_add_mem_of_weierstrassP_eq M hu hu' heq with h | h
  · exact hd h
  · exact hs h

/-- ★★★★★★★★★★★★
**判別式の差積表示（一般の 2-捩れ 3 点で）**——★**無条件**（第 1401）。

☆第 1397 の `latticeDisc_eq_prod_half` を `Λ` の半周期以外の 2-捩れ点にも使える形にした。 -/
theorem latticeDisc_eq_prod_two_torsion (M : PeriodPair) {u₁ u₂ u₃ : ℂ}
    (h1 : u₁ ∉ M.lattice) (h2 : u₂ ∉ M.lattice) (h3 : u₃ ∉ M.lattice)
    (t1 : 2 * u₁ ∈ M.lattice) (t2 : 2 * u₂ ∈ M.lattice) (t3 : 2 * u₃ ∈ M.lattice)
    (d12 : u₁ - u₂ ∉ M.lattice) (s12 : u₁ + u₂ ∉ M.lattice)
    (d13 : u₁ - u₃ ∉ M.lattice) (s13 : u₁ + u₃ ∉ M.lattice)
    (d23 : u₂ - u₃ ∉ M.lattice) (s23 : u₂ + u₃ ∉ M.lattice) :
    latticeDisc M
      = 16 * ((M.weierstrassP u₁ - M.weierstrassP u₂)
        * (M.weierstrassP u₁ - M.weierstrassP u₃)
        * (M.weierstrassP u₂ - M.weierstrassP u₃)) ^ 2 :=
  (vieta_of_three_roots _ _ _ M.g₂ M.g₃
    (weierstrassP_two_torsion_ne M h1 h2 d12 s12)
    (weierstrassP_two_torsion_ne M h1 h3 d13 s13)
    (weierstrassP_two_torsion_ne M h2 h3 d23 s23)
    (cubic_eq_zero_of_two_mem M u₁ h1 t1) (cubic_eq_zero_of_two_mem M u₂ h2 t2)
    (cubic_eq_zero_of_two_mem M u₃ h3 t3)).2.2.2

/-! ## ★★★★`Λ′/Λ` の位数が奇なら 2-捩れは動かない -/

/-- ☆`Λ′/Λ` に 2-捩れが無ければ、`Λ` の 2-捩れ点は `Λ′` にも入らない。 -/
theorem notMem_lattice'_of_two_mem (P P' : PeriodPair)
    (hodd2 : ∀ y ∈ P'.lattice, 2 * y ∈ P.lattice → y ∈ P.lattice)
    {x : ℂ} (hx : x ∉ P.lattice) (h2x : 2 * x ∈ P.lattice) : x ∉ P'.lattice :=
  fun hc => hx (hodd2 x hc h2x)

/-- ★★★★★★★★**`Λ` の半周期は `Λ′` の外にあり、差も和も外にある**（第 1401）。 -/
theorem half_period_data (P P' : PeriodPair)
    (hodd2 : ∀ y ∈ P'.lattice, 2 * y ∈ P.lattice → y ∈ P.lattice) :
    ((P.ω₁ / 2) ∉ P'.lattice ∧ (P.ω₂ / 2) ∉ P'.lattice ∧ ((P.ω₁ + P.ω₂) / 2) ∉ P'.lattice)
    ∧ (P.ω₁ / 2 - P.ω₂ / 2 ∉ P'.lattice ∧ P.ω₁ / 2 + P.ω₂ / 2 ∉ P'.lattice)
    ∧ (P.ω₁ / 2 - (P.ω₁ + P.ω₂) / 2 ∉ P'.lattice
        ∧ P.ω₁ / 2 + (P.ω₁ + P.ω₂) / 2 ∉ P'.lattice)
    ∧ (P.ω₂ / 2 - (P.ω₁ + P.ω₂) / 2 ∉ P'.lattice
        ∧ P.ω₂ / 2 + (P.ω₁ + P.ω₂) / 2 ∉ P'.lattice) := by
  have step : ∀ x : ℂ, x ∉ P.lattice → 2 * x ∈ P.lattice → x ∉ P'.lattice :=
    fun x hx h2x => notMem_lattice'_of_two_mem P P' hodd2 hx h2x
  have n1 : P.ω₁ / 2 ∉ P.lattice := P.ω₁_div_two_notMem_lattice
  have n2 : P.ω₂ / 2 ∉ P.lattice := P.ω₂_div_two_notMem_lattice
  have n3 : (P.ω₁ + P.ω₂) / 2 ∉ P.lattice := add_div_two_notMem_lattice P
  refine ⟨⟨step _ n1 ?_, step _ n2 ?_, step _ n3 ?_⟩, ⟨step _ ?_ ?_, step _ ?_ ?_⟩,
    ⟨step _ ?_ ?_, step _ ?_ ?_⟩, ⟨step _ ?_ ?_, step _ ?_ ?_⟩⟩
  · rw [show 2 * (P.ω₁ / 2) = P.ω₁ by ring]; exact P.ω₁_mem_lattice
  · rw [show 2 * (P.ω₂ / 2) = P.ω₂ by ring]; exact P.ω₂_mem_lattice
  · rw [show 2 * ((P.ω₁ + P.ω₂) / 2) = P.ω₁ + P.ω₂ by ring]
    exact P.lattice.add_mem P.ω₁_mem_lattice P.ω₂_mem_lattice
  · rw [show P.ω₁ / 2 - P.ω₂ / 2 = (P.ω₁ - P.ω₂) / 2 by ring]
    exact sub_div_two_notMem_lattice P
  · rw [show 2 * (P.ω₁ / 2 - P.ω₂ / 2) = P.ω₁ - P.ω₂ by ring]
    exact P.lattice.sub_mem P.ω₁_mem_lattice P.ω₂_mem_lattice
  · rw [show P.ω₁ / 2 + P.ω₂ / 2 = (P.ω₁ + P.ω₂) / 2 by ring]; exact n3
  · rw [show 2 * (P.ω₁ / 2 + P.ω₂ / 2) = P.ω₁ + P.ω₂ by ring]
    exact P.lattice.add_mem P.ω₁_mem_lattice P.ω₂_mem_lattice
  · intro h
    refine n2 ?_
    have he : P.ω₁ / 2 - (P.ω₁ + P.ω₂) / 2 = -(P.ω₂ / 2) := by ring
    rw [he] at h
    simpa using neg_mem h
  · rw [show 2 * (P.ω₁ / 2 - (P.ω₁ + P.ω₂) / 2) = -P.ω₂ by ring]
    exact neg_mem P.ω₂_mem_lattice
  · intro h
    refine n2 ?_
    have he : P.ω₁ / 2 + (P.ω₁ + P.ω₂) / 2 = P.ω₁ + P.ω₂ / 2 := by ring
    rw [he] at h
    simpa using P.lattice.sub_mem h P.ω₁_mem_lattice
  · rw [show 2 * (P.ω₁ / 2 + (P.ω₁ + P.ω₂) / 2) = 2 * P.ω₁ + P.ω₂ by ring]
    exact P.lattice.add_mem (by
      rw [show (2:ℂ) * P.ω₁ = P.ω₁ + P.ω₁ by ring]
      exact P.lattice.add_mem P.ω₁_mem_lattice P.ω₁_mem_lattice) P.ω₂_mem_lattice
  · intro h
    refine n1 ?_
    have he : P.ω₂ / 2 - (P.ω₁ + P.ω₂) / 2 = -(P.ω₁ / 2) := by ring
    rw [he] at h
    simpa using neg_mem h
  · rw [show 2 * (P.ω₂ / 2 - (P.ω₁ + P.ω₂) / 2) = -P.ω₁ by ring]
    exact neg_mem P.ω₁_mem_lattice
  · intro h
    refine n1 ?_
    have he : P.ω₂ / 2 + (P.ω₁ + P.ω₂) / 2 = P.ω₂ + P.ω₁ / 2 := by ring
    rw [he] at h
    simpa using P.lattice.sub_mem h P.ω₂_mem_lattice
  · rw [show 2 * (P.ω₂ / 2 + (P.ω₁ + P.ω₂) / 2) = P.ω₁ + 2 * P.ω₂ by ring]
    exact P.lattice.add_mem P.ω₁_mem_lattice (by
      rw [show (2:ℂ) * P.ω₂ = P.ω₂ + P.ω₂ by ring]
      exact P.lattice.add_mem P.ω₂_mem_lattice P.ω₂_mem_lattice)

/-- ☆代表系の非零元は 2-捩れでない。 -/
theorem two_notMem_of_mem_erase (P P' : PeriodPair) (T : Finset ℂ) (h0T : (0:ℂ) ∈ T)
    (hT : ∀ w ∈ T, w ∈ P'.lattice)
    (hrep : ∀ p ∈ P'.lattice, ∃ w₀ ∈ T, p + w₀ ∈ P.lattice
      ∧ ∀ w ∈ T, w ≠ w₀ → p + w ∉ P.lattice)
    (hodd2 : ∀ y ∈ P'.lattice, 2 * y ∈ P.lattice → y ∈ P.lattice)
    {w : ℂ} (hw : w ∈ T.erase 0) : 2 * w ∉ P.lattice := by
  intro hc
  exact rep_notMem_lattice P P' T h0T hT hrep hw
    (hodd2 w (hT w (Finset.mem_of_mem_erase hw)) hc)

/-! ## ★★★★★★★★6 本の積 -/

set_option maxHeartbeats 1000000 in
/-- ★★★★★★★★★★★★
**6 つの積の積は `(−D)^{l−1}`**——第 1398 を `w` について掛けた形（第 1401）。 -/
theorem prodA_eq (P P' : PeriodPair) (T : Finset ℂ) (h0T : (0:ℂ) ∈ T)
    (hT : ∀ w ∈ T, w ∈ P'.lattice)
    (hrep : ∀ p ∈ P'.lattice, ∃ w₀ ∈ T, p + w₀ ∈ P.lattice
      ∧ ∀ w ∈ T, w ≠ w₀ → p + w ∉ P.lattice)
    (hodd2 : ∀ y ∈ P'.lattice, 2 * y ∈ P.lattice → y ∈ P.lattice) :
    (∏ w ∈ T.erase 0, (P.weierstrassP (P.ω₁ / 2 + w) - P.weierstrassP (P.ω₂ / 2)))
      * (∏ w ∈ T.erase 0,
          (P.weierstrassP (P.ω₁ / 2 + w) - P.weierstrassP ((P.ω₁ + P.ω₂) / 2)))
      * (∏ w ∈ T.erase 0, (P.weierstrassP (P.ω₂ / 2 + w) - P.weierstrassP (P.ω₁ / 2)))
      * (∏ w ∈ T.erase 0,
          (P.weierstrassP (P.ω₂ / 2 + w) - P.weierstrassP ((P.ω₁ + P.ω₂) / 2)))
      * (∏ w ∈ T.erase 0,
          (P.weierstrassP ((P.ω₁ + P.ω₂) / 2 + w) - P.weierstrassP (P.ω₁ / 2)))
      * (∏ w ∈ T.erase 0,
          (P.weierstrassP ((P.ω₁ + P.ω₂) / 2 + w) - P.weierstrassP (P.ω₂ / 2)))
      = (-(((P.weierstrassP (P.ω₁ / 2) - P.weierstrassP (P.ω₂ / 2))
          * (P.weierstrassP (P.ω₁ / 2) - P.weierstrassP ((P.ω₁ + P.ω₂) / 2))
          * (P.weierstrassP (P.ω₂ / 2) - P.weierstrassP ((P.ω₁ + P.ω₂) / 2))) ^ 2))
        ^ (T.erase 0).card := by
  rw [← Finset.prod_mul_distrib, ← Finset.prod_mul_distrib, ← Finset.prod_mul_distrib,
    ← Finset.prod_mul_distrib, ← Finset.prod_mul_distrib]
  rw [Finset.prod_congr rfl (fun w hw =>
    prod_half_shift P (two_notMem_of_mem_erase P P' T h0T hT hrep hodd2 hw))]
  exact Finset.prod_const _

/-- ★★★★★★★★**同種のノルムを 2-捩れ点で使った形**（第 1401）。 -/
theorem key_eq (P P' : PeriodPair) (hle : P.lattice ≤ P'.lattice)
    (T : Finset ℂ) (h0T : (0:ℂ) ∈ T) (hT : ∀ w ∈ T, w ∈ P'.lattice)
    (hrep : ∀ p ∈ P'.lattice, ∃ w₀ ∈ T, p + w₀ ∈ P.lattice
      ∧ ∀ w ∈ T, w ≠ w₀ → p + w ∉ P.lattice)
    {u u' : ℂ} (hu : u ∉ P'.lattice) (hu' : u' ∉ P'.lattice) :
    (P.weierstrassP u' - P.weierstrassP u)
        * (∏ w ∈ T.erase 0, (P.weierstrassP (u' + w) - P.weierstrassP u))
      = (∏ w ∈ T.erase 0, (P.weierstrassP w - P.weierstrassP u))
        * (P'.weierstrassP u' - P'.weierstrassP u) := by
  have h := veluProd_eq P P' hle T h0T hT hrep hu hu'
  rw [veluProd_split P T h0T _ u'] at h
  exact h

set_option maxHeartbeats 1000000 in
/-- ★★★★★★★★**`N² = 4^{l−1}·c₁c₂c₃`**（第 1401）。 -/
theorem derivProd_sq (P P' : PeriodPair) (T : Finset ℂ) (h0T : (0:ℂ) ∈ T)
    (hT : ∀ w ∈ T, w ∈ P'.lattice)
    (hrep : ∀ p ∈ P'.lattice, ∃ w₀ ∈ T, p + w₀ ∈ P.lattice
      ∧ ∀ w ∈ T, w ≠ w₀ → p + w ∉ P.lattice) :
    (∏ w ∈ T.erase 0, P.derivWeierstrassP w) ^ 2
      = 4 ^ (T.erase 0).card
        * ((∏ w ∈ T.erase 0, (P.weierstrassP w - P.weierstrassP (P.ω₁ / 2)))
          * (∏ w ∈ T.erase 0, (P.weierstrassP w - P.weierstrassP (P.ω₂ / 2)))
          * (∏ w ∈ T.erase 0,
              (P.weierstrassP w - P.weierstrassP ((P.ω₁ + P.ω₂) / 2)))) := by
  rw [← Finset.prod_pow]
  rw [Finset.prod_congr rfl (fun w hw =>
    derivSq_factor P (rep_notMem_lattice P P' T h0T hT hrep hw))]
  rw [Finset.prod_mul_distrib, Finset.prod_mul_distrib, Finset.prod_mul_distrib,
    Finset.prod_const]
  ring

/-! ## ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★恒等式 -/

set_option maxHeartbeats 2000000 in
/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**[GenEll] 格子曲線の上の判別式の恒等式**——★**`Λ′/Λ` の位数が奇であることだけ**（第 1401）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

    latticeDisc Λ ^ l = latticeDisc Λ′ · ( ∏_{w ∈ T∖{0}} ℘′_Λ(w) )^4     （`l = 2m+1`）

★★★これが `Skeleton/GenEll/VeluDiscIdentity.lean` の
`disc_pow_eq_veluQuot_mul_lattice`（残っていた唯一の `sorry`）の中身である。 -/
theorem latticeDisc_pow_eq (P P' : PeriodPair) (hle : P.lattice ≤ P'.lattice)
    (T : Finset ℂ) (h0T : (0:ℂ) ∈ T) (hT : ∀ w ∈ T, w ∈ P'.lattice)
    (hrep : ∀ p ∈ P'.lattice, ∃ w₀ ∈ T, p + w₀ ∈ P.lattice
      ∧ ∀ w ∈ T, w ≠ w₀ → p + w ∉ P.lattice)
    (hodd2 : ∀ y ∈ P'.lattice, 2 * y ∈ P.lattice → y ∈ P.lattice)
    {m : ℕ} (hcard : T.card = 2 * m + 1) :
    latticeDisc P ^ (2 * m + 1)
      = latticeDisc P' * (∏ w ∈ T.erase 0, P.derivWeierstrassP w) ^ 4 := by
  classical
  obtain ⟨⟨q1, q2, q3⟩, ⟨qd12, qs12⟩, ⟨qd13, qs13⟩, ⟨qd23, qs23⟩⟩ :=
    half_period_data P P' hodd2
  have hcardS : (T.erase 0).card = 2 * m := by
    rw [Finset.card_erase_of_mem h0T, hcard]; omega
  set e₁ := P.weierstrassP (P.ω₁ / 2) with he₁
  set e₂ := P.weierstrassP (P.ω₂ / 2) with he₂
  set e₃ := P.weierstrassP ((P.ω₁ + P.ω₂) / 2) with he₃
  set E₁ := P'.weierstrassP (P.ω₁ / 2) with hE₁
  set E₂ := P'.weierstrassP (P.ω₂ / 2) with hE₂
  set E₃ := P'.weierstrassP ((P.ω₁ + P.ω₂) / 2) with hE₃
  set c₁ := ∏ w ∈ T.erase 0, (P.weierstrassP w - e₁) with hc₁
  set c₂ := ∏ w ∈ T.erase 0, (P.weierstrassP w - e₂) with hc₂
  set c₃ := ∏ w ∈ T.erase 0, (P.weierstrassP w - e₃) with hc₃
  have K₁ := key_eq P P' hle T h0T hT hrep q2 q1
  have K₂ := key_eq P P' hle T h0T hT hrep q3 q1
  have K₃ := key_eq P P' hle T h0T hT hrep q1 q2
  have K₄ := key_eq P P' hle T h0T hT hrep q3 q2
  have K₅ := key_eq P P' hle T h0T hT hrep q1 q3
  have K₆ := key_eq P P' hle T h0T hT hrep q2 q3
  have hA := prodA_eq P P' T h0T hT hrep hodd2
  rw [hcardS] at hA
  have hDD : latticeDisc P = 16 * ((e₁ - e₂) * (e₁ - e₃) * (e₂ - e₃)) ^ 2 :=
    latticeDisc_eq_prod_half P
  have hDD' : latticeDisc P' = 16 * ((E₁ - E₂) * (E₁ - E₃) * (E₂ - E₃)) ^ 2 := by
    refine latticeDisc_eq_prod_two_torsion P' q1 q2 q3 ?_ ?_ ?_ qd12 qs12 qd13 qs13 qd23 qs23
    · exact hle (by rw [show 2 * (P.ω₁ / 2) = P.ω₁ by ring]; exact P.ω₁_mem_lattice)
    · exact hle (by rw [show 2 * (P.ω₂ / 2) = P.ω₂ by ring]; exact P.ω₂_mem_lattice)
    · exact hle (by
        rw [show 2 * ((P.ω₁ + P.ω₂) / 2) = P.ω₁ + P.ω₂ by ring]
        exact P.lattice.add_mem P.ω₁_mem_lattice P.ω₂_mem_lattice)
  have hN2 := derivProd_sq P P' T h0T hT hrep
  rw [hcardS] at hN2
  have hprodK :
      ((e₁ - e₂) * (∏ w ∈ T.erase 0, (P.weierstrassP (P.ω₁ / 2 + w) - e₂)))
        * ((e₁ - e₃) * (∏ w ∈ T.erase 0, (P.weierstrassP (P.ω₁ / 2 + w) - e₃)))
        * ((e₂ - e₁) * (∏ w ∈ T.erase 0, (P.weierstrassP (P.ω₂ / 2 + w) - e₁)))
        * ((e₂ - e₃) * (∏ w ∈ T.erase 0, (P.weierstrassP (P.ω₂ / 2 + w) - e₃)))
        * ((e₃ - e₁) * (∏ w ∈ T.erase 0, (P.weierstrassP ((P.ω₁ + P.ω₂) / 2 + w) - e₁)))
        * ((e₃ - e₂) * (∏ w ∈ T.erase 0, (P.weierstrassP ((P.ω₁ + P.ω₂) / 2 + w) - e₂)))
      = (c₂ * (E₁ - E₂)) * (c₃ * (E₁ - E₃)) * (c₁ * (E₂ - E₁)) * (c₃ * (E₂ - E₃))
        * (c₁ * (E₃ - E₁)) * (c₂ * (E₃ - E₂)) := by
    rw [K₁, K₂, K₃, K₄, K₅, K₆]
  set DD := ((e₁ - e₂) * (e₁ - e₃) * (e₂ - e₃)) ^ 2 with hDDdef
  set DD' := ((E₁ - E₂) * (E₁ - E₃) * (E₂ - E₃)) ^ 2 with hDD'def
  have hneg : (-DD) ^ (2 * m) = DD ^ (2 * m) := Even.neg_pow ⟨m, by ring⟩ DD
  have hcombine : -DD * DD ^ (2 * m) = (c₁ * c₂ * c₃) ^ 2 * (-DD') := by
    rw [← hneg, ← hA, hDDdef, hDD'def]
    linear_combination hprodK
  have hDDpow : DD ^ (2 * m + 1) = (c₁ * c₂ * c₃) ^ 2 * DD' := by
    linear_combination -hcombine
  have hN4 : (∏ w ∈ T.erase 0, P.derivWeierstrassP w) ^ 4
      = ((4:ℂ) ^ (2 * m) * (c₁ * c₂ * c₃)) ^ 2 := by
    have hsq : (∏ w ∈ T.erase 0, P.derivWeierstrassP w) ^ 4
        = ((∏ w ∈ T.erase 0, P.derivWeierstrassP w) ^ 2) ^ 2 := by ring
    rw [hsq, hN2]
  have h16 : (16:ℂ) ^ (2 * m + 1) = 16 * ((4:ℂ) ^ (2 * m)) ^ 2 := by
    rw [← pow_mul, show (2 * m) * 2 = 2 * (2 * m) by ring, pow_mul]
    norm_num
    ring
  rw [hDD, hDD', hN4, mul_pow, hDDpow, h16]
  ring

/-! ## ★出典の紐付け(`.src`) -/

def weierstrassP_two_torsion_ne.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(u ± u′ ∉ M なら ℘_M(u) ≠ ℘_M(u′)。★無条件)",
    sectionId := "genell-lemma-3-5" }

def latticeDisc_eq_prod_two_torsion.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(判別式の差積表示——一般の 2-捩れ 3 点で。★無条件)",
    sectionId := "genell-lemma-3-5" }

def notMem_lattice'_of_two_mem.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Λ′/Λ に 2-捩れが無ければ Λ の 2-捩れ点は Λ′ の外。★無条件)",
    sectionId := "genell-lemma-3-5" }

def half_period_data.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Λ の半周期は Λ′ の外——差も和も。★Λ′/Λ の位数が奇)",
    sectionId := "genell-lemma-3-5" }

def two_notMem_of_mem_erase.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(代表系の非零元は 2-捩れでない。★Λ′/Λ の位数が奇)",
    sectionId := "genell-lemma-3-5" }

def prodA_eq.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(6 つの積の積は (−D)^{l−1}。★無条件)",
    sectionId := "genell-lemma-3-5" }

def key_eq.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(同種のノルムを 2-捩れ点で使った形。★無条件)",
    sectionId := "genell-lemma-3-5" }

def derivProd_sq.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(N² = 4^{l−1}·c₁c₂c₃。★無条件)",
    sectionId := "genell-lemma-3-5" }

def latticeDisc_pow_eq.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(格子曲線の上の判別式の恒等式。★Λ′/Λ の位数が奇)",
    sectionId := "genell-lemma-3-5" }

def latticeDisc_pow_eq.needs : List ProofObligation :=
  [ .citation "[ABC3]" "veluProd_eq(同種のノルム、第 1400、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.veluProd_eq") 1,
    .citation "[ABC3]" "prod_half_shift(第 1398、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.prod_half_shift") 1,
    .citation "[ABC3]" "latticeDisc_eq_prod_half(第 1397、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.latticeDisc_eq_prod_half") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1401）**——`Δ(E)^l = Δ(E/C)·N⁴` の**格子曲線版が閉じた**。" ++
       "☆帳簿は 3 行:(1) 同種のノルムを `z = v_j` で使う、" ++
       "(2) 6 本の順序対を掛ける、(3) `∏R = (−D)^{l−1} = D^{l−1}`（`l` 奇）。" ++
       "★これで `Skeleton/GenEll/VeluDiscIdentity.lean` の唯一の `sorry` が埋まる。") 17 ]

end ABC3.Found.GenEll
