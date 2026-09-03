/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GaloisRep.PowerSeriesUniversal
import ABC3.Found.GaloisRep.MuPairFunctorial
import ABC3.Skeleton.GenEll.TateODE
import ABC3.Meta.Claim

/-!
# VeluDYDenomFree —— `[GenEll] Definition 3.3` の分

☆もとの 1 枚を**条なしの項目ごと**に割ったものである（第 1458、案 a）。
★「1 ファイル = 1 ノード」を回復するための分割で、中身は動かしていない。
-/

namespace ABC3.Found.GaloisRep

open PowerSeries Finset

/-! ## ★定数項 -/

/-- ☆`(X)` の元の定数項は `0`（`constantCoeff` の形）。 -/
theorem constantCoeff_of_mem_span_X {A : Type} [CommRing A] {f : PowerSeries A}
    (hf : f ∈ Ideal.span {(PowerSeries.X : PowerSeries A)}) :
    PowerSeries.constantCoeff f = 0 := by
  simpa using coeff_zero_of_mem_span_X hf

/-! ## ★★★★★★段 1 —— 体の上（`A₁`）では側条件が定数項で出る -/

/-- ★★★★★★★★**`DX ≠ 0` は定数項で出る**（第 1128）。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

☆`tateDXpair (C α) w X` の定数項は `α(1+α)·(1−α)⁻³` である——
尾（`tateDXtail`）も `w` 側（`w ∈ (X)`）もすべて `(X)` に入るので効かない。
★`TateDXNeZero.lean` の `tateDXpair_ne_zero_of_mu` は `TateSetup` と Tate 一意化を要するが、
`A₁` では**その重装備は要らない**。 -/
theorem tateDXpair_C_ne_zero {K : Type} [Field K] {α : K}
    (hα0 : α ≠ 0) (hαneg : α + 1 ≠ 0) (hα1 : α ≠ 1)
    {w : PowerSeries K} (hw : w ∈ Ideal.span {(PowerSeries.X : PowerSeries K)}) :
    tateDXpair (PowerSeries.C α) w PowerSeries.X (Ideal.mem_span_singleton_self _) ≠ 0 := by
  have hu : IsUnit (1 - α) := isUnit_iff_ne_zero.2 (sub_ne_zero.2 (Ne.symm hα1))
  intro h0
  have hcc := congrArg (PowerSeries.constantCoeff (R := K)) h0
  rw [tateDXpair, map_sub, map_add, map_add] at hcc
  have h1 : PowerSeries.constantCoeff (tateDXterm (PowerSeries.C α : PowerSeries K))
      = tateDXterm α := by
    rw [← map_tateDXterm (PowerSeries.C (R := K)) hu]; simp
  have h2 : PowerSeries.constantCoeff
      (tateDXtail (PowerSeries.C α : PowerSeries K) PowerSeries.X
        (Ideal.mem_span_singleton_self _)) = 0 :=
    constantCoeff_of_mem_span_X (tateDXtail_mem _ _ _)
  have h3mem : tateDXterm w ∈ Ideal.span {(PowerSeries.X : PowerSeries K)} := by
    have hh := tateDXterm_mem_pow (I := Ideal.span {(PowerSeries.X : PowerSeries K)})
      (k := 1) (t := w) (by simpa using hw)
    simpa using hh
  have h3 : PowerSeries.constantCoeff (tateDXterm w) = 0 := constantCoeff_of_mem_span_X h3mem
  have h4 : PowerSeries.constantCoeff
      (tateDXtail w PowerSeries.X (Ideal.mem_span_singleton_self _)) = 0 :=
    constantCoeff_of_mem_span_X (tateDXtail_mem _ _ _)
  rw [h1, h2, h3, h4] at hcc
  simp only [add_zero, sub_zero, map_zero] at hcc
  have hinv : Ring.inverse (1 - α) ≠ 0 := by
    rw [Ring.inverse_eq_inv]
    exact inv_ne_zero (sub_ne_zero.2 (Ne.symm hα1))
  rw [tateDXterm] at hcc
  have hne : α * (1 + α) * Ring.inverse (1 - α) ^ 3 ≠ 0 := by
    refine mul_ne_zero (mul_ne_zero hα0 ?_) (pow_ne_zero _ hinv)
    rw [add_comm]; exact hαneg
  exact hne hcc

def tateDXpair_C_ne_zero.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(PowerSeries の上では DX ≠ 0 が定数項で出ること)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
