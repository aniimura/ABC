/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GaloisRep.MuGraded
import ABC3.Found.GaloisRep.TateVelu
import ABC3.Meta.Claim

/-!
# VeluMuSum —— `[GenEll] Lemma 3.2` の分

☆もとの 1 枚を**条なしの項目ごと**に割ったものである（第 1458、案 a）。
★「1 ファイル = 1 ノード」を回復するための分割で、中身は動かしていない。
-/

namespace ABC3.Found.GaloisRep

open Finset ABC3.Found.GenEll
variable {R : Type} [CommRing R] {I : Ideal R}

/-- ★★★★★★`a₄(q)` の ζ-free な μ-等級付き係数。 -/
noncomputable def a4C (l : ℕ) (q : R) (n a : ℕ) : R :=
  if a = 0 then ((PowerSeries.coeff n tateA4 : ℤ) : R) * q ^ n else 0

theorem a4C_mem {l : ℕ} {q : R} (hq : q ∈ I) (n a : ℕ) : a4C (R := R) l q n a ∈ I ^ n := by
  classical
  by_cases h : a = 0
  · simpa [a4C, h] using Ideal.mul_mem_left _ _ (Ideal.pow_mem_pow hq n)
  · simpa [a4C, h] using Submodule.zero_mem (I ^ n)

theorem a4_eq_muEval [IsAdicComplete I R] {l : ℕ} (hl : 0 < l) (q : R) (hq : q ∈ I) (z : R) :
    (tateCurveAt q hq).a₄ = muEval (I := I) l (a4C l q) (a4C_mem hq) z := by
  classical
  rw [tateCurveAt_a₄, evalAdic_eq_adicSum,
    adicSum_eq_muEval (l := l) hl _ (fun n => Ideal.mul_mem_left _ _ (Ideal.pow_mem_pow hq n)) z]
  exact muEval_congr _ _ _ _ (fun n a => by by_cases h : a = 0 <;> simp [a4C, h]) z

/-! ## ★★★★★★★★★★★★★★★★★★★★`∑_ζ X(ζ)` の閉じた形 -/

open Finset in
/-- ★★★★★★★★★★★★★★★★★★★★★★
**`∑_{ζ∈μ_l∖{1}} X(ζ, q)` の閉じた形**。

★古典形の作用で全係数和が `n ≥ 1` で消えるので、
`sum_mu_muEval'` の右辺は **`l·A n 0` だけ**になる。 -/
theorem sum_mu_X_closed [IsAdicComplete I R] [IsDomain R] {l : ℕ} (hl : l.Prime)
    (hlu : IsUnit ((l : R))) {ζ : R} (hζ : IsPrimitiveRoot ζ l)
    (hu : ∀ i ∈ (range l).erase 0, IsUnit (1 - ζ ^ i)) (q : R) (hq : q ∈ I) :
    ∑ i ∈ (range l).erase 0, tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq
      = adicSum (I := I)
          (fun n => (l : R) * tateXC l q n 0
            - (if n = 0 then (Ring.inverse ((l : R))) ^ 2 *
                ((∑ k ∈ range l, (k : R)) * ∑ m ∈ range l, (m : R)) else 0))
          (fun n => Submodule.sub_mem _
            (Ideal.mul_mem_left _ _ (tateXC_mem hq n 0))
            (by
              by_cases h : n = 0
              · simpa [h] using Submodule.mem_top
                  (x := (Ring.inverse ((l : R))) ^ 2 *
                    ((∑ k ∈ range l, (k : R)) * ∑ m ∈ range l, (m : R)))
              · simpa [h] using Submodule.zero_mem (I ^ n))) := by
  classical
  have hterm : ∀ i ∈ (range l).erase 0,
      tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq
        = muEval (I := I) l (tateXC l q) (tateXC_mem hq) (ζ ^ i) := by
    intro i hi
    have hi0 : i ≠ 0 := (Finset.mem_erase.1 hi).1
    have hil : i < l := Finset.mem_range.1 (Finset.mem_erase.1 hi).2
    have hnd : ¬ l ∣ i := fun h => hi0 (Nat.eq_zero_of_dvd_of_lt h hil)
    have hpow : (ζ ^ i) ^ l = 1 := by
      rw [← pow_mul, mul_comm, pow_mul, hζ.pow_eq_one, one_pow]
    have hsum : ∑ k ∈ range l, (ζ ^ i) ^ k = 0 :=
      (isPrimitiveRoot_pow_of_not_dvd (R := R) hl hζ hnd).geom_sum_eq_zero hl.one_lt
    exact tateXpair_eq_muEval hl.pos hlu (hu i hi) hpow hsum q hq
  rw [Finset.sum_congr rfl hterm, sum_mu_muEval' hl hζ (tateXC l q) (tateXC_mem hq)]
  refine adicSum_congr _ _ (fun n => ?_)
  rw [sum_tateXC hl.pos]

def sum_mu_X_closed.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(∑_ζ X(ζ,q) の閉じた形。★無条件)",
    sectionId := "genell-lemma-3-2" }

def a4C.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(a₄(q) の ζ-free な μ-等級付き係数)",
    sectionId := "genell-lemma-3-2" }

def a4C_mem.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(a4C の所属。★無条件)",
    sectionId := "genell-lemma-3-2" }

def a4_eq_muEval.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(a₄(q) の μ-等級付き形。★無条件)",
    sectionId := "genell-lemma-3-2" }

end ABC3.Found.GaloisRep
