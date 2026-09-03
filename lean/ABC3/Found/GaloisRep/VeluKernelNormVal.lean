/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.DeepDoubling
import ABC3.Found.GaloisRep.VeluKernelNorm
import ABC3.Meta.Claim

/-!
# 第 1396 ブロック —— **`3 ∣ v_p(N)`（積まで組み上げた形）**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★★★★★これは何か——点ごとの段を積に上げる

第 1395 は**核の点 1 個**について `3 ∣ v_p(t_P)` を与えた。
★`N = ∏_{P ∈ C∖{O}} t_P`（第 1386 `veluKernelNorm`）なので、
`v_p` は和になり **`3 ∣ v_p(N)`** が出る。

☆同時に **`N ≠ 0`** も出る——核の点は奇位数なので 2-捩れでなく、`t_P ≠ 0` である。

## ★★★★★★★★これで第 1393 に渡る形が揃った

    dvd_minDeltaExp_of_disc_pow_eq
      : v_p(Δ(E)) = 0 → Δ(E)^l = Δ(E′)·N⁴ → 3 ∣ v_p(N) → 12 ∣ minDeltaExp p E′

の**第 3 引数**が本ブロックである。★残るのは
**判別式の恒等式 `Δ(E)^l = Δ(E′)·N⁴`**（`Skeleton/GenEll/VeluDiscIdentity.lean`）と、
`12 ∣ minDeltaExp` から `= 0` への段だけになった。
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve Polynomial IsDedekindDomain NumberField Finset ABC3.Meta
open WeierstrassCurve.Affine ABC3.Found.GenEll

open scoped Classical

variable {L : Type} [Field L] [NumberField L]

/-- ★★★★★★**各因子の `valAdd` が 3 の倍数なら積の `valAdd` も**（第 1396）。 -/
theorem three_dvd_valAdd_prod {ι : Type} (p : HeightOneSpectrum (𝓞 L))
    (S : Finset ι) (f : ι → L)
    (h : ∀ q ∈ S, ∃ hq : f q ≠ 0, (3 : ℤ) ∣ valAdd p (Units.mk0 (f q) hq)) :
    ∃ hP : (∏ q ∈ S, f q) ≠ 0, (3 : ℤ) ∣ valAdd p (Units.mk0 (∏ q ∈ S, f q) hP) := by
  classical
  induction S using Finset.induction with
  | empty =>
      refine ⟨by simp, ?_⟩
      have hz : valAdd p (Units.mk0 (∏ q ∈ (∅ : Finset ι), f q) (by simp)) = 0 := by
        rw [show (Units.mk0 (∏ q ∈ (∅ : Finset ι), f q) (by simp)) = 1 from
          Units.ext (by simp)]
        exact valAdd_one p
      rw [hz]
      exact dvd_zero 3
  | insert a S ha ih =>
      obtain ⟨hane, hadvd⟩ := h a (Finset.mem_insert_self a S)
      obtain ⟨hPne, hPdvd⟩ := ih (fun q hq => h q (Finset.mem_insert_of_mem hq))
      have hprodeq : ∏ q ∈ insert a S, f q = f a * ∏ q ∈ S, f q :=
        Finset.prod_insert ha
      have hne : (∏ q ∈ insert a S, f q) ≠ 0 := by
        rw [hprodeq]; exact mul_ne_zero hane hPne
      refine ⟨hne, ?_⟩
      have hmulne : f a * ∏ q ∈ S, f q ≠ 0 := mul_ne_zero hane hPne
      rw [valAdd_congr p hne hmulne hprodeq, valAdd_mulL p hane hPne hmulne]
      exact dvd_add hadvd hPdvd

/-- ★★★★★★★★★★★★★★★★★★★★
**[GenEll] 良い素点では `N ≠ 0` かつ `3 ∣ v_p(N)`**——★（第 1396）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆点ごとの段（第 1395）を積に上げただけである。

★★★これが第 1393 の `dvd_minDeltaExp_of_disc_pow_eq` の第 3 引数である。 -/
theorem three_dvd_valAdd_veluKernelNorm (p : HeightOneSpectrum (𝓞 L))
    (E : WeierstrassCurve L) [WeierstrassCurve.IsIntegral (primeSubring p) E]
    (hΔ : E.Δ ≠ 0) (hΔ0 : valAdd p (Units.mk0 E.Δ hΔ) = 0)
    (h2 : valAdd p (Units.mk0 (2 : L) two_ne_zero) = 0)
    {l : ℕ} (hl : l.Prime) (hodd : l ≠ 2)
    (Q : E.toAffine.Point) (hQ : addOrderOf Q = l) :
    ∃ hN : veluKernelNorm E
        (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q))) ≠ 0,
      (3 : ℤ) ∣ valAdd p (Units.mk0 (veluKernelNorm E
        (((range l).erase 0).image (fun k : ℕ => pointCoords (k • Q)))) hN) := by
  rw [veluKernelNorm]
  refine three_dvd_valAdd_prod p _ _ ?_
  intro q hq
  obtain ⟨k, hk, rfl⟩ := Finset.mem_image.1 hq
  rw [mem_erase, mem_range] at hk
  have hkne : k • Q ≠ 0 := nsmul_ne_zero_of_lt_addOrderOf hQ hk.1 hk.2
  have hlz : l • Q = 0 := by rw [← hQ]; exact addOrderOf_nsmul_eq_zero Q
  have hdvd : addOrderOf (k • Q) ∣ l := by
    refine addOrderOf_dvd_of_nsmul_eq_zero ?_
    rw [smul_comm, hlz, smul_zero]
  have hord : addOrderOf (k • Q) = l := by
    rcases (Nat.Prime.eq_one_or_self_of_dvd hl _ hdvd) with h1 | h1
    · exact absurd (AddMonoid.addOrderOf_eq_one_iff.1 h1) hkne
    · exact h1
  rcases hkQ : k • Q with _ | ⟨x, y, h⟩
  · exact absurd hkQ hkne
  · have hord' : addOrderOf (WeierstrassCurve.Affine.Point.some x y h) = l := hkQ ▸ hord
    have hty : y ≠ E.toAffine.negY x y :=
      negYdiff_ne_zero_of_addOrderOf_prime E h hl hodd hord'
    have htne : y - E.toAffine.negY x y ≠ 0 := sub_ne_zero.mpr hty
    have hfac : 2 * y + E.a₁ * x + E.a₃ = y - E.toAffine.negY x y := by
      rw [WeierstrassCurve.Affine.negY]; ring
    have hne : 2 * y + E.a₁ * x + E.a₃ ≠ 0 := by rw [hfac]; exact htne
    rw [pointCoords_some]
    refine ⟨hne, ?_⟩
    rw [valAdd_congr p hne htne hfac]
    exact three_dvd_valAdd_negYdiff_of_addOrderOf_prime p E hΔ hΔ0 h2 h hl hodd hord' hty

/-! ## ★出典の紐付け(`.src`) -/

def three_dvd_valAdd_prod.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(各因子の valAdd が 3 の倍数なら積の valAdd も。★無条件)",
    sectionId := "genell-lemma-3-5" }

def three_dvd_valAdd_veluKernelNorm.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(良い素点では N ≠ 0 かつ 3 ∣ v_p(N)。★v_p(2) = 0)",
    sectionId := "genell-lemma-3-5" }

def three_dvd_valAdd_veluKernelNorm.needs : List ProofObligation :=
  [ .citation "[ABC3]" "three_dvd_valAdd_negYdiff_of_addOrderOf_prime(第 1395、証明済み)"
      (.inProject "ABC3"
        "ABC3.Found.GaloisRep.three_dvd_valAdd_negYdiff_of_addOrderOf_prime") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1396）**——第 1393 の " ++
       "`dvd_minDeltaExp_of_disc_pow_eq` の第 3 引数が揃った。" ++
       "☆残るのは判別式の恒等式 `Δ(E)^l = Δ(E′)·N⁴` と、" ++
       "`12 ∣ minDeltaExp` から `= 0` への段である。") 17 ]

end ABC3.Found.GaloisRep
