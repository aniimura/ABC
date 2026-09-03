/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Skeleton.GenEll.TateIsogeny
import ABC3.Found.GaloisRep.VeluODEAny
import ABC3.Meta.Claim

/-!
# 第 1146 ブロック —— **`c₄`・`c₆` の恒等式を `hDX` なしで**（`Skeleton`、`l = 2` の枝の節点 2）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> parameter qE of E satisfies the relation qE = qEl ; in particular, we have

## ★★★★★★★★これは何か

`c4_velu_tate`・`c6_velu_tate` は `hDX : ∀ i, tateDXpair ≠ 0` を仮説に取る。
★`l = 2` では点が 2-捻れなので `DX = 0` になり、この仮説は**偽**である。

☆第 1144-1145 で `veluV2 = DY` と高階の ODE が `hDX` なしで出たので、
本ブロックで `c₄`・`c₆` も `hDX` なしになる。

★代わりに使うのは `IsUnit (ζ^i)` であり、これは `ζ^l = 1` から自動である。
-/

namespace ABC3.Skeleton.GenEll

open ABC3.Meta ABC3.Found.GaloisRep ABC3.Found.GenEll WeierstrassCurve Finset

open Finset in
/-- ★★★★★★★★**`∑ veluV2 = ∑ DY`——`hDX` を取らない**（第 1146）。 -/
theorem sum_veluV2_eq_sum_tateDYpair_any {R : Type} [CommRing R] [IsDomain R] {I : Ideal R}
    [IsAdicComplete I R] {l : ℕ} (hl : 0 < l) {ζ : R} (hζl : ζ ^ l = 1)
    (hu : ∀ i ∈ (range l).erase 0, IsUnit (1 - ζ ^ i))
    (q : R) (hq : q ∈ I)
    :
    (∑ i ∈ (range l).erase 0,
        veluV2 (tateCurveAt q hq) (tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
          (tateYpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq))
      = ∑ i ∈ (range l).erase 0, tateDYpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq := by
  refine Finset.sum_congr rfl (fun i hi => ?_)
  have hpow : (ζ ^ i) * (ζ ^ i) ^ (l - 1) = 1 := by
    rw [← pow_succ']
    rw [Nat.sub_add_cancel hl, ← pow_mul, mul_comm, pow_mul, hζl, one_pow]
  have haw : (ζ ^ i) * (q * (ζ ^ i) ^ (l - 1)) = q := by
    calc (ζ ^ i) * (q * (ζ ^ i) ^ (l - 1)) = q * ((ζ ^ i) * (ζ ^ i) ^ (l - 1)) := by ring
      _ = q := by rw [hpow, mul_one]
  have hwu : IsUnit (1 - q * (ζ ^ i) ^ (l - 1)) :=
    isUnit_one_sub (I := I) (Ideal.mul_mem_right _ _ hq)
  exact ABC3.Found.GaloisRep.veluV2_eq_tateDYpair_any _ _ q hq haw
    (IsUnit.of_mul_eq_one _ hpow) (hu i hi) hwu

def sum_veluV2_eq_sum_tateDYpair_any.src : Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(v = ∑ DY——★DX ≠ 0 を取らない)",
    sectionId := "genell-lemma-3-2" }

def sum_veluV2_eq_sum_tateDYpair_any.needs : List ProofObligation :=
  [ .citation "[ABC3]" "veluV2_eq_tateDYpair_any(第 1144、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.veluV2_eq_tateDYpair_any") 1 ]

open Finset in
/-- ★★★★★★★★★★**`c₄` の恒等式——`hDX` を取らない**（第 1146）。 -/
theorem c4_velu_tate_any {R : Type} [CommRing R] [IsDomain R] [CharZero R] {I : Ideal R}
    [IsAdicComplete I R] {l : ℕ} (hl : l.Prime) {ζ : R} (hζ : IsPrimitiveRoot ζ l)
    (hlu : IsUnit ((l : R))) (hu : ∀ i ∈ (range l).erase 0, IsUnit (1 - ζ ^ i))
    (q : R) (hq : q ∈ I) (hql : q ^ l ∈ I) (h2 : (2 : R) ≠ 0)
    :
    (tateCurveAt q hq).c₄
        + 240 * (∑ i ∈ (range l).erase 0,
            veluV2 (tateCurveAt q hq) (tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
              (tateYpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq))
      = (l : R) ^ 4 * (tateCurveAt (q ^ l) hql).c₄ := by
  have hζl : ζ ^ l = 1 := hζ.pow_eq_one
  have hsum1 := sum_veluV2_eq_sum_tateDYpair_any hl.pos hζl hu q hq
  have hstep : ∑ i ∈ (range l).erase 0, tateD2Xpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq
      = 2 * (∑ i ∈ (range l).erase 0, tateDYpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
        + ∑ i ∈ (range l).erase 0, tateDXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq := by
    rw [Finset.mul_sum, ← Finset.sum_add_distrib]
    refine Finset.sum_congr rfl (fun i hi => ?_)
    exact tateD2Xpair_eq _ _ q hq (hu i hi)
      (isUnit_one_sub (I := I) (Ideal.mul_mem_right _ _ hq))
  have hzero := sum_mu_dxpair_zero hl hζ hu q hq h2
  have hd2 := sum_mu_d2xpair hl hζ hu q hq hql
  rw [hstep, hzero, add_zero] at hd2
  rw [hsum1, tateCurveAt_c4_eq, tateCurveAt_c4_eq]
  linear_combination hd2

def c4_velu_tate_any.src : Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(c₄ の恒等式。★DX ≠ 0 を取らない)",
    sectionId := "genell-lemma-3-2" }

def c4_velu_tate_any.needs : List ProofObligation :=
  [ .citation "[ABC3]" "sum_veluV2_eq_sum_tateDYpair_any(第 1146、証明済み)"
      (.inProject "ABC3" "ABC3.Skeleton.GenEll.sum_veluV2_eq_sum_tateDYpair_any") 1 ]

set_option maxHeartbeats 2000000 in
open Finset in
/-- ★★★★★★★★★★**`c₆` の恒等式——`hDX` を取らない**（第 1146）。 -/
theorem c6_velu_tate_any {R : Type} [CommRing R] [IsDomain R] [CharZero R] {I : Ideal R}
    [IsAdicComplete I R] {l : ℕ} (hl : l.Prime) {ζ : R} (hζ : IsPrimitiveRoot ζ l)
    (hu : ∀ i ∈ (range l).erase 0, IsUnit (1 - ζ ^ i))
    (q : R) (hq : q ∈ I) (hql : q ^ l ∈ I) (h2 : (2 : R) ≠ 0)
    :
    (tateCurveAt q hq).c₆
        + 504 * (∑ i ∈ (range l).erase 0,
            veluV2 (tateCurveAt q hq) (tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
              (tateYpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq))
        + 3024 * (∑ i ∈ (range l).erase 0,
            (veluU (tateCurveAt q hq) (tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
                (tateYpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
              + 2 * (veluV2 (tateCurveAt q hq)
                      (tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
                      (tateYpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
                    * tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)))
      = (l : R) ^ 6 * (tateCurveAt (q ^ l) hql).c₆ := by
  have hζl : ζ ^ l = 1 := hζ.pow_eq_one
  have hawi : ∀ i ∈ (range l).erase 0,
      (ζ ^ i) * (q * (ζ ^ i) ^ (l - 1)) = q := by
    intro i hi
    have hpow : (ζ ^ i) * (ζ ^ i) ^ (l - 1) = 1 := by
      rw [← pow_succ']
      rw [Nat.sub_add_cancel hl.pos, ← pow_mul, mul_comm, pow_mul, hζl, one_pow]
    calc (ζ ^ i) * (q * (ζ ^ i) ^ (l - 1)) = q * ((ζ ^ i) * (ζ ^ i) ^ (l - 1)) := by ring
      _ = q := by rw [hpow, mul_one]
  have hζu : ∀ i ∈ (range l).erase 0, IsUnit (ζ ^ i) := fun i hi =>
    IsUnit.of_mul_eq_one _ (by
      have h1 : (ζ ^ i) * (ζ ^ i) ^ (l - 1) = 1 := by
        rw [← pow_succ']
        rw [Nat.sub_add_cancel hl.one_lt.le, ← pow_mul, mul_comm, pow_mul, hζl, one_pow]
      exact h1)
  have hwu : ∀ i : ℕ, IsUnit (1 - q * (ζ ^ i) ^ (l - 1)) := fun i =>
    isUnit_one_sub (I := I) (Ideal.mul_mem_right _ _ hq)
  -- 項ごとの恒等式
  have hterm : ∀ i ∈ (range l).erase 0,
      12 * (veluU (tateCurveAt q hq) (tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
              (tateYpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
            + 2 * (veluV2 (tateCurveAt q hq)
                    (tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
                    (tateYpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
                  * tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq))
        = tateD4Xpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq
          - tateD3Xpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq
          - tateD2Xpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq
          + tateDXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq := by
    intro i hi
    have hd4 := ABC3.Found.GaloisRep.tate_d4x_any (ζ ^ i)
      (q * (ζ ^ i) ^ (l - 1)) q hq (hawi i hi) (hζu i hi) (hu i hi)
    have hd3 := ABC3.Found.GaloisRep.tate_d3x_any (ζ ^ i)
      (q * (ζ ^ i) ^ (l - 1)) q hq (hawi i hi) (hζu i hi) (hu i hi)
    have hd2 := ABC3.Found.GaloisRep.tate_d2x_any (ζ ^ i)
      (q * (ζ ^ i) ^ (l - 1)) q hq (hawi i hi) (hζu i hi) (hu i hi)
    have hdx := tateDXpair_eq (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq (hu i hi) (hwu i)
    rw [veluU_tateCurveAt, veluV2_tateCurveAt, hd4, hd3, hd2, hdx]
    ring
  -- 和に直す
  have hsum12 : 12 * (∑ i ∈ (range l).erase 0,
      (veluU (tateCurveAt q hq) (tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
              (tateYpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
            + 2 * (veluV2 (tateCurveAt q hq)
                    (tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
                    (tateYpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
                  * tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)))
      = (∑ i ∈ (range l).erase 0, tateD4Xpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
        - (∑ i ∈ (range l).erase 0, tateD3Xpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
        - (∑ i ∈ (range l).erase 0, tateD2Xpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
        + (∑ i ∈ (range l).erase 0, tateDXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq) := by
    rw [Finset.mul_sum, Finset.sum_congr rfl hterm, Finset.sum_add_distrib,
      Finset.sum_sub_distrib, Finset.sum_sub_distrib]
  -- ∑DX = 0、∑D³X = 0
  have hz1 := sum_mu_dxpair_zero hl hζ hu q hq h2
  have hz3 := sum_mu_d3xpair_zero hl hζ hu q hq h2
  -- ∑v = ∑DY と ∑D²X = 2∑DY + ∑DX
  have hsum1 := sum_veluV2_eq_sum_tateDYpair_any hl.pos hζl hu q hq
  have hstep : ∑ i ∈ (range l).erase 0, tateD2Xpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq
      = 2 * (∑ i ∈ (range l).erase 0, tateDYpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
        + ∑ i ∈ (range l).erase 0, tateDXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq := by
    rw [Finset.mul_sum, ← Finset.sum_add_distrib]
    refine Finset.sum_congr rfl (fun i hi => ?_)
    exact tateD2Xpair_eq _ _ q hq (hu i hi) (hwu i)
  have hd4sum := sum_mu_d4xpair hl hζ hu q hq hql
  rw [hz1, add_zero] at hstep
  rw [hz1, hz3, sub_zero, add_zero] at hsum12
  rw [hsum1, tateCurveAt_c6_eq, tateCurveAt_c6_eq]
  linear_combination (252 : R) * hsum12 + hd4sum - (252 : R) * hstep

def c6_velu_tate_any.src : Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(c₆ の恒等式。★DX ≠ 0 を取らない)",
    sectionId := "genell-lemma-3-2" }

def c6_velu_tate_any.needs : List ProofObligation :=
  [ .citation "[ABC3]" "tate_d4x_any / tate_d3x_any / tate_d2x_any(第 1145、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.tate_d4x_any") 1,
    .implicitStep
      ("★★★★**2026-09-01（第 1146）**——`l = 2` では点が 2-捻れなので " ++
       "`hDX : tateDXpair ≠ 0` は**偽**である。" ++
       "☆第 1144-1145 で径数の道を通したので、`c₄`・`c₆` も `hDX` なしになった。") 8 ]

end ABC3.Skeleton.GenEll
