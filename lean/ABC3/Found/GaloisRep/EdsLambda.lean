import ABC3.Found.GaloisRep.Wronskian

/-!
# Galois (G1) 第 56 ブロック —— **★★★★★★EDS の第 2 の 1 変数関係 (Λ)**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★Ward の定理へ向かう部品

`normEDS` が `IsEllSequence` を満たすこと(mathlib の TODO)を、
★**`n` についての帰納**で示す道を §9-396 で測った。帰納段は

    (S₁D₁)(S₋₁D₋₁) = (S₁D₋₁)(S₋₁D₁)

という**積の並べ替え**で回るが、そのとき残る差が

    B₀B₁B₋₁·Λ(A) − A₀A₁A₋₁·Λ(B),   Λ(A) := b²A₀³ + A₁²A₋₂ + A₂A₋₁²

であった。★★つまり **`Λ(m)/(W(m−1)W(m)W(m+1))` が `m` に依らない**ことが要る。

## ★★★実測した定数

    c·Λ(m) = (b⁴ + d)·W(m−1)W(m)W(m+1)

★`preNormEDS` の水準では(第 1 引数を `b` と書いて)

    c(ε(n)P(n)³ + P(n+1)²P(n−2) + P(n+2)P(n−1)²) = (b + d)P(n−1)P(n)P(n+1)

    ε(n) = b (n 偶), 1 (n 奇)

## ★★★★`normEDSRec` の窓とまた合った

Python(sympy)で測ったところ、(Λ) at `2k` / `2k+1` は
仮定「(Λ) と Somos-4 at `k−1, k, k+1`」に対して**4 ケースとも余り 0** ✅

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `LamAt` | ★(Λ) の主張 |
| `lam_even` / `lam_odd` | ★★帰納段(証明書は Python で算出) |
| `lam` | ★★★★★★**(Λ)——任意の可換環・任意の整数** |
-/

namespace ABC3.Found.GaloisRep

universe u

variable {R : Type u} [CommRing R] (b c d : R)

/-- ★★★★★★**(Λ)**——`Λ(n)/(P(n−1)P(n)P(n+1))` が `n` に依らないことの分母払い形。 -/
def LamAt (b c d : R) (n : ℤ) : Prop :=
  c * ((if Even n then b else 1) * preNormEDS b c d n ^ 3
      + preNormEDS b c d (n + 1) ^ 2 * preNormEDS b c d (n - 2)
      + preNormEDS b c d (n + 2) * preNormEDS b c d (n - 1) ^ 2)
    = (b + d) * preNormEDS b c d (n - 1) * preNormEDS b c d n * preNormEDS b c d (n + 1)

/-! ## ★基底 -/

theorem lam_zero : LamAt b c d 0 := by simp only [LamAt]; norm_num

theorem lam_one : LamAt b c d 1 := by
  simp only [LamAt]; rw [if_neg (by decide : ¬ Even (1 : ℤ))]; norm_num

theorem lam_two : LamAt b c d 2 := by
  simp only [LamAt]; rw [if_pos (by decide : Even (2 : ℤ))]; norm_num; ring

theorem lam_three : LamAt b c d 3 := by
  simp only [LamAt]; rw [if_neg (by decide : ¬ Even (3 : ℤ))]
  norm_num [pEDS5]; ring

theorem lam_four : LamAt b c d 4 := by
  simp only [LamAt]; rw [if_pos (by decide : Even (4 : ℤ))]
  norm_num [pEDS5, pEDS6]; ring

/-! ## ★★帰納段 -/

/-- ★★**偶数段**。証明書は Python(sympy)の多項式除算で算出した(§9-396)。 -/
theorem lam_even (k : ℤ)
    (l1 : LamAt b c d (k + 1)) (l0 : LamAt b c d k) (lm : LamAt b c d (k - 1)) :
    LamAt b c d (2 * k) := by
  have s1 := psomos b c d (k + 1)
  have s0 := psomos b c d k
  have sm := psomos b c d (k - 1)
  simp only [LamAt] at l1 l0 lm ⊢
  simp only [PSomosAt] at s1 s0 sm
  rw [show k + 1 + 1 = k + 2 by ring, show k + 1 - 2 = k - 1 by ring,
    show k + 1 + 2 = k + 3 by ring, show k + 1 - 1 = k by ring] at l1 s1
  rw [show k - 1 + 1 = k by ring, show k - 1 - 2 = k - 3 by ring,
    show k - 1 + 2 = k + 1 by ring, show k - 1 - 1 = k - 2 by ring] at lm sm
  rw [show 2 * k - 2 = 2 * (k - 1) by ring,
    show 2 * k + 2 = 2 * (k + 1) by ring, show 2 * k - 1 = 2 * (k - 1) + 1 by ring,
    if_pos (even_two_mul k)]
  have E1 := preNormEDS_even b c d (k + 1)
  have E0 := preNormEDS_even b c d k
  have Em := preNormEDS_even b c d (k - 1)
  have O0 := preNormEDS_odd b c d k
  have Om := preNormEDS_odd b c d (k - 1)
  rw [show k + 1 - 1 = k by ring, show k + 1 + 2 = k + 3 by ring,
    show k + 1 - 2 = k - 1 by ring, show k + 1 + 1 = k + 2 by ring] at E1
  rw [show k - 1 - 1 = k - 2 by ring, show k - 1 + 2 = k + 1 by ring,
    show k - 1 - 2 = k - 3 by ring, show k - 1 + 1 = k by ring] at Em
  rw [show k - 1 + 2 = k + 1 by ring, show k - 1 - 1 = k - 2 by ring,
    show k - 1 + 1 = k by ring] at Om
  rw [E1, E0, Em, O0, Om]
  set nm3 := preNormEDS b c d (k - 3)
  set nm2 := preNormEDS b c d (k - 2)
  set nm1 := preNormEDS b c d (k - 1)
  set n0 := preNormEDS b c d k
  set np1 := preNormEDS b c d (k + 1)
  set np2 := preNormEDS b c d (k + 2)
  set np3 := preNormEDS b c d (k + 3)
  by_cases hk : Even k
  · have hk1 : ¬ Even (k + 1) := by rw [Int.even_add_one]; exact not_not.mpr hk
    have hkm : ¬ Even (k - 1) := by rw [Int.even_sub_one]; exact not_not.mpr hk
    simp only [if_pos hk, if_neg hk1, if_neg hkm, one_mul, mul_one] at l1 l0 lm s1 s0 sm ⊢
    linear_combination (b^2*nm2^2*n0^6*np1 - 2*b*nm2*nm1^3*n0^3*np1^2 + nm1^6*np1^3) * l1
      + (-b^2*nm1*n0^6*np2^2 + 2*b*nm1^2*n0^3*np1^3*np2 - nm1^3*np1^6) * lm
      + (-b*nm2^2*n0^3*np1^4 + b*nm1^4*n0^3*np2^2 + 2*nm2*nm1^3*np1^5
          - 2*nm1^5*np1^3*np2) * l0
  · have hk1 : Even (k + 1) := Int.even_add_one.mpr hk
    have hkm : Even (k - 1) := Int.even_sub_one.mpr hk
    simp only [if_neg hk, if_pos hk1, if_pos hkm, one_mul, mul_one] at l1 l0 lm s1 s0 sm ⊢
    linear_combination (b^2*nm1^6*np1^3 - 2*b*nm2*nm1^3*n0^3*np1^2 + nm2^2*n0^6*np1) * l1
      + (-b^2*nm1^3*np1^6 + 2*b*nm1^2*n0^3*np1^3*np2 - nm1*n0^6*np2^2) * lm
      + (2*b^2*nm2*nm1^3*np1^5 - 2*b^2*nm1^5*np1^3*np2 - b*nm2^2*n0^3*np1^4
          + b*nm1^4*n0^3*np2^2) * l0

/-- ★★**奇数段**。 -/
theorem lam_odd (k : ℤ)
    (l1 : LamAt b c d (k + 1)) (l0 : LamAt b c d k) :
    LamAt b c d (2 * k + 1) := by
  have s1 := psomos b c d (k + 1)
  have s0 := psomos b c d k
  simp only [LamAt] at l1 l0 ⊢
  simp only [PSomosAt] at s1 s0
  rw [show k + 1 + 1 = k + 2 by ring, show k + 1 - 2 = k - 1 by ring,
    show k + 1 + 2 = k + 3 by ring, show k + 1 - 1 = k by ring] at l1 s1
  rw [show 2 * k + 1 + 1 = 2 * (k + 1) by ring, show 2 * k + 1 - 2 = 2 * (k - 1) + 1 by ring,
    show 2 * k + 1 + 2 = 2 * (k + 1) + 1 by ring, show 2 * k + 1 - 1 = 2 * k by ring,
    if_neg k.not_even_two_mul_add_one]
  have E1 := preNormEDS_even b c d (k + 1)
  have E0 := preNormEDS_even b c d k
  have O1 := preNormEDS_odd b c d (k + 1)
  have O0 := preNormEDS_odd b c d k
  have Om := preNormEDS_odd b c d (k - 1)
  rw [show k + 1 - 1 = k by ring, show k + 1 + 2 = k + 3 by ring,
    show k + 1 - 2 = k - 1 by ring, show k + 1 + 1 = k + 2 by ring] at E1
  rw [show k + 1 + 2 = k + 3 by ring, show k + 1 - 1 = k by ring,
    show k + 1 + 1 = k + 2 by ring] at O1
  rw [show k - 1 + 2 = k + 1 by ring, show k - 1 - 1 = k - 2 by ring,
    show k - 1 + 1 = k by ring] at Om
  rw [E1, E0, O1, O0, Om]
  set nm2 := preNormEDS b c d (k - 2)
  set nm1 := preNormEDS b c d (k - 1)
  set n0 := preNormEDS b c d k
  set np1 := preNormEDS b c d (k + 1)
  set np2 := preNormEDS b c d (k + 2)
  set np3 := preNormEDS b c d (k + 3)
  by_cases hk : Even k
  · have hk1 : ¬ Even (k + 1) := by rw [Int.even_add_one]; exact not_not.mpr hk
    have hkm : ¬ Even (k - 1) := by rw [Int.even_sub_one]; exact not_not.mpr hk
    simp only [if_pos hk, if_neg hk1, if_neg hkm, one_mul, mul_one] at l1 l0 s1 s0 O1 O0 Om ⊢
    linear_combination (3*b*nm2*nm1*n0^3*np1^2*np2^2 - b*nm2*n0^5*np1^2*np3 + b*nm2*n0^3*np1^5
        + nm2^2*np1^7 - 2*nm2*nm1^2*np1^5*np2 - 2*nm1^4*np1^3*np2^2 + nm1^3*n0^2*np1^3*np3
        - nm1^3*np1^6) * l1
      + (b^2*n0^6*np2^3 - b*nm2*n0^3*np1^2*np2^3 - b*nm1^2*n0^3*np2^4 - 3*b*nm1*n0^3*np1^3*np2^2
        - nm2*nm1*np1^5*np2^2 - nm2*np1^8 + 3*nm1^3*np1^3*np2^3 + 3*nm1^2*np1^6*np2) * l0
      + (-b^2*nm1*n0^6*np1*np2^2 - b*d*nm1*n0^6*np1*np2^2 - b*nm2*n0^3*np1^6
        + 2*b*nm1^2*n0^3*np1^4*np2 - d*nm2*n0^3*np1^6 + 2*d*nm1^2*n0^3*np1^4*np2) * s1
      + (b^2*nm1*n0^4*np1^3*np2^2 + b*d*nm1*n0^4*np1^3*np2^2 + b*nm2*n0*np1^8
        - 2*b*nm1^2*n0*np1^6*np2 + d*nm2*n0*np1^8 - 2*d*nm1^2*n0*np1^6*np2) * s0
  · have hk1 : Even (k + 1) := Int.even_add_one.mpr hk
    have hkm : Even (k - 1) := Int.even_sub_one.mpr hk
    simp only [if_neg hk, if_pos hk1, if_pos hkm, one_mul, mul_one] at l1 l0 s1 s0 O1 O0 Om ⊢
    linear_combination (-b^2*nm1^3*np1^6 + b*nm2^2*np1^7 - 2*b*nm2*nm1^2*np1^5*np2
        + b*nm2*n0^3*np1^5 - 2*b*nm1^4*np1^3*np2^2 + b*nm1^3*n0^2*np1^3*np3
        + 3*nm2*nm1*n0^3*np1^2*np2^2 - nm2*n0^5*np1^2*np3) * l1
      + (-b^2*nm2*np1^8 + 3*b^2*nm1^2*np1^6*np2 - b*nm2*nm1*np1^5*np2^2 + 3*b*nm1^3*np1^3*np2^3
        - 3*b*nm1*n0^3*np1^3*np2^2 - nm2*n0^3*np1^2*np2^3 - nm1^2*n0^3*np2^4 + n0^6*np2^3) * l0
      + (-b^2*nm2*n0^3*np1^6 + 2*b^2*nm1^2*n0^3*np1^4*np2 - b*d*nm2*n0^3*np1^6
        + 2*b*d*nm1^2*n0^3*np1^4*np2 - b*nm1*n0^6*np1*np2^2 - d*nm1*n0^6*np1*np2^2) * s1
      + (b^2*nm2*n0*np1^8 - 2*b^2*nm1^2*n0*np1^6*np2 + b*d*nm2*n0*np1^8
        - 2*b*d*nm1^2*n0*np1^6*np2 + b*nm1*n0^4*np1^3*np2^2 + d*nm1*n0^4*np1^3*np2^2) * s0

/-! ## ★★★★★★(Λ) 本体 -/

theorem lam_nat (n : ℕ) : LamAt b c d (n : ℤ) := by
  induction n using normEDSRec with
  | zero => exact_mod_cast lam_zero b c d
  | one => exact_mod_cast lam_one b c d
  | two => exact_mod_cast lam_two b c d
  | three => exact_mod_cast lam_three b c d
  | four => exact_mod_cast lam_four b c d
  | even m h1 h2 h3 h4 h5 =>
      have key := lam_even b c d ((m : ℤ) + 3)
        (by rw [show (m : ℤ) + 3 + 1 = ((m + 4 : ℕ) : ℤ) by push_cast; ring]; exact h4)
        (by rw [show (m : ℤ) + 3 = ((m + 3 : ℕ) : ℤ) by push_cast; ring]; exact h3)
        (by rw [show (m : ℤ) + 3 - 1 = ((m + 2 : ℕ) : ℤ) by push_cast; ring]; exact h2)
      rw [show ((2 * (m + 3) : ℕ) : ℤ) = 2 * ((m : ℤ) + 3) by push_cast; ring]
      exact key
  | odd m h1 h2 h3 h4 =>
      have key := lam_odd b c d ((m : ℤ) + 2)
        (by rw [show (m : ℤ) + 2 + 1 = ((m + 3 : ℕ) : ℤ) by push_cast; ring]; exact h3)
        (by rw [show (m : ℤ) + 2 = ((m + 2 : ℕ) : ℤ) by push_cast; ring]; exact h2)
      rw [show ((2 * (m + 2) + 1 : ℕ) : ℤ) = 2 * ((m : ℤ) + 2) + 1 by push_cast; ring]
      exact key

/-- ★★★★★★**(Λ)——任意の可換環・任意の整数**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★Ward の定理(`normEDS` が `IsEllSequence` を満たすこと)の帰納段で使う。 -/
theorem lam (n : ℤ) : LamAt b c d n := by
  rcases le_or_gt 0 n with h | h
  · lift n to ℕ using h; exact lam_nat b c d n
  · have key := lam_nat b c d (-n).toNat
    rw [Int.toNat_of_nonneg (by omega)] at key
    simp only [LamAt] at key ⊢
    rw [show -n + 1 = -(n - 1) by ring, show -n - 2 = -(n + 2) by ring,
      show -n + 2 = -(n - 2) by ring, show -n - 1 = -(n + 1) by ring] at key
    simp only [preNormEDS_neg, even_neg] at key
    linear_combination -key

/-! ## ★出典の紐付け(`.src`) -/

def lam.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(EDS の第 2 の 1 変数関係)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
