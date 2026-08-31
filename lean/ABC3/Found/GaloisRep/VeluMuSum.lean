/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.MuGraded
import ABC3.Found.GaloisRep.TateVelu
import ABC3.Meta.Claim

/-!
# Vélu の和 `v` を μ-等級付きの枠で（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15–p.17。

原文 (GenEll p.15):
> parameter qE of E satisfies the relation qE = qEl ; in particular, we have

## ★★★★★★★★これは何か

`tateModel_of_quot_mu`（葉 1）の**手順 3**——Vélu の和

    `v = ∑_{ζ ∈ μ_l∖{1}} (3X(ζ)² + a₄ − Y(ζ))`

を、`Found/GaloisRep/MuGraded.lean` の枠の中で書く。

★`X(ζ)`・`Y(ζ)` はどちらも単一の `muEval`（第 807・809）なので、
`muEval_mul`（畳み込み）・`muEval_add`・`muEval_smul`・`muEval_sub` で
`v` の被加数も単一の `muEval` になり、`sum_mu_muEval'` により

    `∑_ζ muEval A ζ = adicSum (n ↦ l·A n 0 − ∑_{a<l} A n a)`

という**有限個の係数の計算**に落ちる。
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

/-- ★★★★★★★★★★★★★★★★**`v` の被加数の μ-等級付き係数**。 -/
noncomputable def veluVC (l : ℕ) (q : R) (n a : ℕ) : R :=
  3 * muConv l (tateXC l q) (tateXC l q) n a + a4C l q n a - tateYC l q n a

theorem veluVC_mem {l : ℕ} {q : R} (hq : q ∈ I) (n a : ℕ) : veluVC l q n a ∈ I ^ n :=
  Submodule.sub_mem _
    (Submodule.add_mem _
      (Ideal.mul_mem_left _ _ (muConv_mem (tateXC_mem hq) (tateXC_mem hq) n a))
      (a4C_mem hq n a))
    (tateYC_mem hq n a)

/-- ★★★★★★★★★★★★★★★★★★★★★★
**Vélu の被加数 `3X(ζ)² + a₄ − Y(ζ)` は単一の μ-等級付き級数**。 -/
theorem veluV2_tate_eq_muEval [IsAdicComplete I R] [IsDomain R] {l : ℕ} (hl : 0 < l)
    (hlu : IsUnit ((l : R))) {z : R} (hu : IsUnit (1 - z)) (hz : z ^ l = 1)
    (hsum : ∑ k ∈ range l, z ^ k = 0) (q : R) (hq : q ∈ I) :
    veluV2 (tateCurveAt q hq) (tateXpair z (q * z ^ (l - 1)) q hq)
        (tateYpair z (q * z ^ (l - 1)) q hq)
      = muEval (I := I) l (veluVC l q) (veluVC_mem hq) z := by
  classical
  rw [veluV2_tateCurveAt, tateXpair_eq_muEval hl hlu hu hz hsum q hq,
    tateYpair_eq_muEval hl hlu hu hz hsum q hq,
    a4_eq_muEval (I := I) hl q hq z, pow_two,
    muEval_mul hl _ _ (tateXC_mem hq) (tateXC_mem hq) hz,
    muEval_smul, muEval_add, muEval_sub]
  refine muEval_congr _ _ _ _ (fun n a => ?_) z
  simp only [veluVC]

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★
**`v = ∑_ζ (3X(ζ)² + a₄ − Y(ζ))` は有限個の係数の計算に落ちる**。

★★★これが `tateModel_of_quot_mu`（葉 1）の手順 3 の到達点である。 -/
theorem sum_mu_veluV2 [IsAdicComplete I R] [IsDomain R] {l : ℕ} (hl : l.Prime)
    (hlu : IsUnit ((l : R))) {ζ : R} (hζ : IsPrimitiveRoot ζ l)
    (hu : ∀ i ∈ (range l).erase 0, IsUnit (1 - ζ ^ i)) (q : R) (hq : q ∈ I) :
    ∑ i ∈ (range l).erase 0,
        veluV2 (tateCurveAt q hq) (tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
          (tateYpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
      = adicSum (I := I)
          (fun n => (l : R) * veluVC l q n 0 - ∑ a ∈ range l, veluVC l q n a)
          (fun n => Submodule.sub_mem _
            (Ideal.mul_mem_left _ _ (veluVC_mem hq n 0))
            (Submodule.sum_mem _ (fun a _ => veluVC_mem hq n a))) := by
  classical
  have hterm : ∀ i ∈ (range l).erase 0,
      veluV2 (tateCurveAt q hq) (tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
          (tateYpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
        = muEval (I := I) l (veluVC l q) (veluVC_mem hq) (ζ ^ i) := by
    intro i hi
    have hi0 : i ≠ 0 := (Finset.mem_erase.1 hi).1
    have hil : i < l := Finset.mem_range.1 (Finset.mem_erase.1 hi).2
    have hnd : ¬ l ∣ i := fun h => hi0 (Nat.eq_zero_of_dvd_of_lt h hil)
    have hpow : (ζ ^ i) ^ l = 1 := by
      rw [← pow_mul, mul_comm, pow_mul, hζ.pow_eq_one, one_pow]
    have hsum : ∑ k ∈ range l, (ζ ^ i) ^ k = 0 :=
      (isPrimitiveRoot_pow_of_not_dvd (R := R) hl hζ hnd).geom_sum_eq_zero hl.one_lt
    exact veluV2_tate_eq_muEval hl.pos hlu (hu i hi) hpow hsum q hq
  rw [Finset.sum_congr rfl hterm, sum_mu_muEval' hl hζ (veluVC l q) (veluVC_mem hq)]

/-! ## ★★★★★★★★★★★★★★★★★★★★`w` の側 -/

/-- ★`2Y + X` の μ-等級付き係数。 -/
noncomputable def twoYplusXC (l : ℕ) (q : R) (n a : ℕ) : R :=
  2 * tateYC l q n a + tateXC l q n a

theorem twoYplusXC_mem {l : ℕ} {q : R} (hq : q ∈ I) (n a : ℕ) :
    twoYplusXC l q n a ∈ I ^ n :=
  Submodule.add_mem _ (Ideal.mul_mem_left _ _ (tateYC_mem hq n a)) (tateXC_mem hq n a)

theorem twoYplusX_eq_muEval [IsAdicComplete I R] [IsDomain R] {l : ℕ} (hl : 0 < l)
    (hlu : IsUnit ((l : R))) {z : R} (hu : IsUnit (1 - z)) (hz : z ^ l = 1)
    (hsum : ∑ k ∈ range l, z ^ k = 0) (q : R) (hq : q ∈ I) :
    2 * tateYpair z (q * z ^ (l - 1)) q hq + tateXpair z (q * z ^ (l - 1)) q hq
      = muEval (I := I) l (twoYplusXC l q) (twoYplusXC_mem hq) z := by
  classical
  rw [tateXpair_eq_muEval hl hlu hu hz hsum q hq,
    tateYpair_eq_muEval hl hlu hu hz hsum q hq, muEval_smul, muEval_add]
  exact muEval_congr _ _ _ _ (fun n a => rfl) z

/-- ★★★★★★★★★★★★★★★★**`2w` の被加数の μ-等級付き係数**。

☆Vélu の `w_Q = u_Q/2 + g^x_Q·x_Q` は体の中の式だが、
**2 倍すれば環の中の式**になる:
`2w_Q = u_Q + 2·g^x_Q·x_Q = (2Y+X)² + 2(3X²+a₄−Y)X` -/
noncomputable def veluWC (l : ℕ) (q : R) (n a : ℕ) : R :=
  muConv l (twoYplusXC l q) (twoYplusXC l q) n a
  + 2 * muConv l (veluVC l q) (tateXC l q) n a

theorem veluWC_mem {l : ℕ} {q : R} (hq : q ∈ I) (n a : ℕ) : veluWC l q n a ∈ I ^ n :=
  Submodule.add_mem _
    (muConv_mem (twoYplusXC_mem hq) (twoYplusXC_mem hq) n a)
    (Ideal.mul_mem_left _ _ (muConv_mem (veluVC_mem hq) (tateXC_mem hq) n a))

/-- ★★★★★★★★★★★★★★★★★★★★
**`2w` の被加数も単一の μ-等級付き級数**。 -/
theorem veluW_tate_eq_muEval [IsAdicComplete I R] [IsDomain R] {l : ℕ} (hl : 0 < l)
    (hlu : IsUnit ((l : R))) {z : R} (hu : IsUnit (1 - z)) (hz : z ^ l = 1)
    (hsum : ∑ k ∈ range l, z ^ k = 0) (q : R) (hq : q ∈ I) :
    veluU (tateCurveAt q hq) (tateXpair z (q * z ^ (l - 1)) q hq)
        (tateYpair z (q * z ^ (l - 1)) q hq)
      + 2 * (veluV2 (tateCurveAt q hq) (tateXpair z (q * z ^ (l - 1)) q hq)
              (tateYpair z (q * z ^ (l - 1)) q hq)
            * tateXpair z (q * z ^ (l - 1)) q hq)
      = muEval (I := I) l (veluWC l q) (veluWC_mem hq) z := by
  classical
  rw [veluU_tateCurveAt, twoYplusX_eq_muEval hl hlu hu hz hsum q hq,
    veluV2_tate_eq_muEval hl hlu hu hz hsum q hq,
    tateXpair_eq_muEval hl hlu hu hz hsum q hq,
    pow_two, muEval_mul hl _ _ (twoYplusXC_mem hq) (twoYplusXC_mem hq) hz,
    muEval_mul hl _ _ (veluVC_mem hq) (tateXC_mem hq) hz,
    muEval_smul, muEval_add]
  exact muEval_congr _ _ _ _ (fun n a => rfl) z

/-- ★★★★★★★★★★★★★★★★★★★★★★★★
**`2w` も有限個の係数の計算に落ちる**。 -/
theorem sum_mu_veluW [IsAdicComplete I R] [IsDomain R] {l : ℕ} (hl : l.Prime)
    (hlu : IsUnit ((l : R))) {ζ : R} (hζ : IsPrimitiveRoot ζ l)
    (hu : ∀ i ∈ (range l).erase 0, IsUnit (1 - ζ ^ i)) (q : R) (hq : q ∈ I) :
    ∑ i ∈ (range l).erase 0,
        (veluU (tateCurveAt q hq) (tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
            (tateYpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
          + 2 * (veluV2 (tateCurveAt q hq) (tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
                  (tateYpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
                * tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq))
      = adicSum (I := I)
          (fun n => (l : R) * veluWC l q n 0 - ∑ a ∈ range l, veluWC l q n a)
          (fun n => Submodule.sub_mem _
            (Ideal.mul_mem_left _ _ (veluWC_mem hq n 0))
            (Submodule.sum_mem _ (fun a _ => veluWC_mem hq n a))) := by
  classical
  have hterm : ∀ i ∈ (range l).erase 0,
      (veluU (tateCurveAt q hq) (tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
          (tateYpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
        + 2 * (veluV2 (tateCurveAt q hq) (tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
                (tateYpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq)
              * tateXpair (ζ ^ i) (q * (ζ ^ i) ^ (l - 1)) q hq))
        = muEval (I := I) l (veluWC l q) (veluWC_mem hq) (ζ ^ i) := by
    intro i hi
    have hi0 : i ≠ 0 := (Finset.mem_erase.1 hi).1
    have hil : i < l := Finset.mem_range.1 (Finset.mem_erase.1 hi).2
    have hnd : ¬ l ∣ i := fun h => hi0 (Nat.eq_zero_of_dvd_of_lt h hil)
    have hpow : (ζ ^ i) ^ l = 1 := by
      rw [← pow_mul, mul_comm, pow_mul, hζ.pow_eq_one, one_pow]
    have hsum : ∑ k ∈ range l, (ζ ^ i) ^ k = 0 :=
      (isPrimitiveRoot_pow_of_not_dvd (R := R) hl hζ hnd).geom_sum_eq_zero hl.one_lt
    exact veluW_tate_eq_muEval hl.pos hlu (hu i hi) hpow hsum q hq
  rw [Finset.sum_congr rfl hterm, sum_mu_muEval' hl hζ (veluWC l q) (veluWC_mem hq)]

/-! ## ★★★★★★★★★★★★★★★★`veluVC`・`veluWC` の 2 つの量 -/

open Finset in
/-- ★★★★★★★★★★★★★★★★**`veluVC` の全係数和**。 -/
theorem sum_veluVC {l : ℕ} (hl : 0 < l) (q : R) (n : ℕ) :
    ∑ a ∈ range l, veluVC l q n a
      = 3 * (∑ k ∈ range (n + 1), (∑ a ∈ range l, tateXC l q k a)
              * ∑ b ∈ range l, tateXC l q (n - k) b)
        + ((PowerSeries.coeff n tateA4 : ℤ) : R) * q ^ n
        - ∑ a ∈ range l, tateYC l q n a := by
  classical
  simp only [veluVC]
  rw [Finset.sum_sub_distrib, Finset.sum_add_distrib]
  congr 1
  congr 1
  · rw [← Finset.mul_sum, sum_muConv hl]
  · simp only [a4C]
    rw [Finset.sum_eq_single 0]
    · simp
    · intro b _ hb
      simp [hb]
    · intro h
      exact absurd (Finset.mem_range.2 hl) h

open Finset in
/-- ★★★★★★★★★★★★★★★★**`veluVC` の 0 次係数**。 -/
theorem veluVC_zero {l : ℕ} (hl : 0 < l) (q : R) (n : ℕ) :
    veluVC l q n 0
      = 3 * (∑ k ∈ range (n + 1), ∑ a ∈ range l,
              tateXC l q k a * tateXC l q (n - k) ((l - a) % l))
        + ((PowerSeries.coeff n tateA4 : ℤ) : R) * q ^ n
        - tateYC l q n 0 := by
  classical
  simp only [veluVC, a4C, muConv_zero hl]
  norm_num

open Finset in
/-- ★★★★★★★★★★★★★★★★**`twoYplusXC` の 2 つの量**。 -/
theorem sum_twoYplusXC {l : ℕ} (q : R) (n : ℕ) :
    ∑ a ∈ range l, twoYplusXC l q n a
      = 2 * (∑ a ∈ range l, tateYC l q n a) + ∑ a ∈ range l, tateXC l q n a := by
  classical
  simp only [twoYplusXC]
  rw [Finset.sum_add_distrib, ← Finset.mul_sum]

theorem twoYplusXC_zero {l : ℕ} (q : R) (n : ℕ) :
    twoYplusXC l q n 0 = 2 * tateYC l q n 0 + tateXC l q n 0 := rfl

open Finset in
/-- ★★★★★★★★★★★★★★★★**`veluWC` の全係数和**。 -/
theorem sum_veluWC {l : ℕ} (hl : 0 < l) (q : R) (n : ℕ) :
    ∑ a ∈ range l, veluWC l q n a
      = (∑ k ∈ range (n + 1), (∑ a ∈ range l, twoYplusXC l q k a)
            * ∑ b ∈ range l, twoYplusXC l q (n - k) b)
        + 2 * (∑ k ∈ range (n + 1), (∑ a ∈ range l, veluVC l q k a)
            * ∑ b ∈ range l, tateXC l q (n - k) b) := by
  classical
  simp only [veluWC]
  rw [Finset.sum_add_distrib, sum_muConv hl, ← Finset.mul_sum, sum_muConv hl]

open Finset in
/-- ★★★★★★★★★★★★★★★★**`veluWC` の 0 次係数**。 -/
theorem veluWC_zero {l : ℕ} (hl : 0 < l) (q : R) (n : ℕ) :
    veluWC l q n 0
      = (∑ k ∈ range (n + 1), ∑ a ∈ range l,
            twoYplusXC l q k a * twoYplusXC l q (n - k) ((l - a) % l))
        + 2 * (∑ k ∈ range (n + 1), ∑ a ∈ range l,
            veluVC l q k a * tateXC l q (n - k) ((l - a) % l)) := by
  classical
  simp only [veluWC, muConv_zero hl]

/-! ## ★出典の紐付け(`.src`) -/

def twoYplusXC.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(2Y + X の μ-等級付き係数)",
    sectionId := "genell-lemma-3-5" }

def twoYplusXC_mem.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(twoYplusXC の所属。★無条件)",
    sectionId := "genell-lemma-3-5" }

def twoYplusX_eq_muEval.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(2Y + X の μ-等級付き形。★無条件)",
    sectionId := "genell-lemma-3-5" }

def veluWC.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(2w の被加数の μ-等級付き係数)",
    sectionId := "genell-lemma-3-5" }

def veluWC_mem.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(veluWC の所属。★無条件)",
    sectionId := "genell-lemma-3-5" }

def veluW_tate_eq_muEval.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(2w の被加数も単一の μ-等級付き級数。★無条件)",
    sectionId := "genell-lemma-3-5" }

def sum_mu_veluW.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(2w も有限個の係数の計算に落ちる。★無条件)",
    sectionId := "genell-lemma-3-5" }

def sum_veluVC.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(veluVC の全係数和。★無条件)",
    sectionId := "genell-lemma-3-5" }

def veluVC_zero.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(veluVC の 0 次係数。★無条件)",
    sectionId := "genell-lemma-3-5" }

def sum_twoYplusXC.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(twoYplusXC の全係数和。★無条件)",
    sectionId := "genell-lemma-3-5" }

def twoYplusXC_zero.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(twoYplusXC の 0 次係数。★無条件)",
    sectionId := "genell-lemma-3-5" }

def sum_veluWC.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(veluWC の全係数和。★無条件)",
    sectionId := "genell-lemma-3-5" }

def veluWC_zero.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(veluWC の 0 次係数。★無条件)",
    sectionId := "genell-lemma-3-5" }

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

def veluVC.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(v の被加数の μ-等級付き係数)",
    sectionId := "genell-lemma-3-5" }

def veluVC_mem.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(veluVC の所属。★無条件)",
    sectionId := "genell-lemma-3-5" }

def veluV2_tate_eq_muEval.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Vélu の被加数は単一の μ-等級付き級数。★無条件)",
    sectionId := "genell-lemma-3-5" }

def sum_mu_veluV2.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(v = ∑_ζ(3X² + a₄ − Y) は有限個の係数の計算に落ちる。★無条件)",
    sectionId := "genell-lemma-3-5" }

end ABC3.Found.GaloisRep
