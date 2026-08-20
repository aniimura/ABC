import ABC3.Found.GaloisRep.MulFormula

/-!
# Galois (G1) 第 47 ブロック —— **★★★★★★★EDS の Somos-4 関係**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★§9-392 で「律速」と判定したものが 1 本に絞れた

§9-393 で乗法公式の帰納から退化を消したあと、残ったのは **y 側**であった。
★§9-394 で測ると、y 側に要るのは 3 項関係の**ただ 1 つの実例**

    ψ_{n+3}ψ_{n−1}ψ_n² − ψ_{n+2}ψ_{n−2}ψ_{n+1}² = ψ_{2n+1}ψ₂²     …(★)

であり、★★さらに (★) は **Somos-4 関係**

    W(m+2)W(m−2) = W(m+1)W(m−1)W(2)² − W(3)W(1)W(m)²             …(R2)

**ただ 1 本**から出る(§9-394 で実測)。

## ★★★★★★そして (R2) は `normEDSRec` の窓とぴったり合う

Python(sympy)で測ったところ:

| 段 | 帰納の仮定 | 余り |
|---|---|---|
| `m = 2k` | (R2) at `k−1, k, k+1` | ★**0** |
| `m = 2k+1` | (R2) at `k−1, k, k+1` | ★**0** |

★`normEDSRec` は偶数段で `P(m+1)…P(m+5)`、奇数段で `P(m+1)…P(m+4)` をくれる
——**どちらも `k−1, k, k+1` を含む** ✅

## ★★★★★`preNormEDS` の水準で書くと割り算が要らない

`normEDS` で書くと偶数段で `b²` を割る必要があり、**整域が要る**。
★★**`preNormEDS` で書けば `b` の因子が最初から出ない**——
`preNormEDS_even` に `b` が現れないからである。

    PSomosAt n : P(n+2)P(n−2) = P(n+1)P(n−1)·ε(n) − c·P(n)²,  ε(n) = 1 (n 偶), b (n 奇)

★★★したがって**任意の可換環**でそのまま通る。整域も普遍環も要らない。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `PSomosAt` | ★Somos-4 の主張 |
| `psomos_even` / `psomos_odd` | ★★帰納段(証明書は Python で算出) |
| `psomos` | ★★★★**任意の可換環・任意の整数で Somos-4** |
| `eds_star` | ★★★★★**(★)——y 側に要る 3 項関係の実例** |
| `normEDS_somos` | ★★`normEDS` 版(mathlib の TODO の一部) |
-/

namespace ABC3.Found.GaloisRep

universe u

variable {R : Type u} [CommRing R] (b c d : R)

/-! ## ★Somos-4 の主張 -/

/-- ★**Somos-4 (R2) の主張**——`preNormEDS` の水準で書く。

`ε(n) = 1` (`n` 偶)、`b` (`n` 奇) が付くのは、`normEDS` から `b` の因子を
外した分の帳尻である。 -/
def PSomosAt (b c d : R) (n : ℤ) : Prop :=
  preNormEDS b c d (n + 2) * preNormEDS b c d (n - 2)
    = preNormEDS b c d (n + 1) * preNormEDS b c d (n - 1) * (if Even n then 1 else b)
      - c * preNormEDS b c d n ^ 2

/-! ## ★基底 -/

theorem pEDS5 : preNormEDS b c d 5 = d * b - c ^ 3 := by
  rw [show (5 : ℤ) = 2 * 2 + 1 by ring, preNormEDS_odd]; norm_num

theorem pEDS6 : preNormEDS b c d 6 = c * (d * b - c ^ 3) - c * d ^ 2 := by
  rw [show (6 : ℤ) = 2 * 3 by ring, preNormEDS_even, show (3 : ℤ) + 2 = 5 by ring,
    show (3 : ℤ) - 1 = 2 by ring, show (3 : ℤ) - 2 = 1 by ring, show (3 : ℤ) + 1 = 4 by ring,
    pEDS5]
  norm_num

theorem psomos_zero : PSomosAt b c d 0 := by simp only [PSomosAt]; norm_num

theorem psomos_one : PSomosAt b c d 1 := by simp only [PSomosAt]; norm_num

theorem psomos_two : PSomosAt b c d 2 := by simp only [PSomosAt]; norm_num

theorem psomos_three : PSomosAt b c d 3 := by
  simp only [PSomosAt]
  rw [if_neg (by decide : ¬ Even (3 : ℤ))]
  norm_num [pEDS5]; ring

theorem psomos_four : PSomosAt b c d 4 := by
  simp only [PSomosAt]
  rw [if_pos (by decide : Even (4 : ℤ))]
  norm_num [pEDS5, pEDS6]; ring

/-! ## ★★帰納段 -/

/-- ★★**偶数段**——`P(k−1), P(k), P(k+1)` から `P(2k)`。

証明書は Python(sympy)の多項式除算で算出した(§9-394)。 -/
theorem psomos_even (k : ℤ)
    (h1 : PSomosAt b c d (k + 1)) (h0 : PSomosAt b c d k) (hm : PSomosAt b c d (k - 1)) :
    PSomosAt b c d (2 * k) := by
  simp only [PSomosAt] at h1 h0 hm ⊢
  rw [show k + 1 + 2 = k + 3 by ring, show k + 1 - 2 = k - 1 by ring,
    show k + 1 + 1 = k + 2 by ring, show k + 1 - 1 = k by ring] at h1
  rw [show k - 1 + 2 = k + 1 by ring, show k - 1 - 2 = k - 3 by ring,
    show k - 1 + 1 = k by ring, show k - 1 - 1 = k - 2 by ring] at hm
  rw [show 2 * k + 2 = 2 * (k + 1) by ring, show 2 * k - 2 = 2 * (k - 1) by ring,
    show 2 * k - 1 = 2 * (k - 1) + 1 by ring, if_pos (even_two_mul k), mul_one]
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
  set nm3 := preNormEDS b c d (k - 3)
  set nm2 := preNormEDS b c d (k - 2)
  set nm1 := preNormEDS b c d (k - 1)
  set n0 := preNormEDS b c d k
  set np1 := preNormEDS b c d (k + 1)
  set np2 := preNormEDS b c d (k + 2)
  set np3 := preNormEDS b c d (k + 3)
  set Q := preNormEDS b c d (2 * (k - 1))
  set S := preNormEDS b c d (2 * k)
  set Vm := preNormEDS b c d (2 * (k - 1) + 1)
  by_cases hk : Even k
  · have hk1 : ¬ Even (k + 1) := by rw [Int.even_add_one]; exact not_not.mpr hk
    have hkm : ¬ Even (k - 1) := by rw [Int.even_sub_one]; exact not_not.mpr hk
    simp only [if_pos hk, if_neg hk1, if_neg hkm, mul_one] at h1 h0 hm O0 Om
    linear_combination Q * E1 + (n0^2*np1*np3 - nm1*np1*np2^2) * Em
      - Vm * O0 - (np2*n0^3*b - nm1*np1^3) * Om
      + c * (S + (nm1^2*n0*np2 - nm2*n0*np1^2)) * E0
      + (nm2^2*n0^2*np1^2 - nm3*n0^4*np1) * h1
      + (nm1^2*n0^2*np2^2 - n0^5*np2*b + n0^4*np1^2*c) * hm
      + (-nm1^3*np1^3 - nm1^2*nm2*np1^2*np2 + nm1^2*n0^3*np2*b - nm1^2*n0^2*np1^2*c
          + nm2*n0^3*np1^2*b) * h0
  · have hk1 : Even (k + 1) := Int.even_add_one.mpr hk
    have hkm : Even (k - 1) := Int.even_sub_one.mpr hk
    simp only [if_neg hk, if_pos hk1, if_pos hkm, mul_one] at h1 h0 hm O0 Om
    linear_combination Q * E1 + (n0^2*np1*np3 - nm1*np1*np2^2) * Em
      - Vm * O0 - (np2*n0^3 - nm1*np1^3*b) * Om
      + c * (S + (nm1^2*n0*np2 - nm2*n0*np1^2)) * E0
      + (nm2^2*n0^2*np1^2 - nm3*n0^4*np1) * h1
      + (nm1^2*n0^2*np2^2 - n0^5*np2 + n0^4*np1^2*c) * hm
      + (-nm1^3*np1^3*b - nm1^2*nm2*np1^2*np2 + nm1^2*n0^3*np2 - nm1^2*n0^2*np1^2*c
          + nm2*n0^3*np1^2) * h0

/-- ★★**奇数段**——`P(k), P(k+1)` から `P(2k+1)`。 -/
theorem psomos_odd (k : ℤ)
    (h1 : PSomosAt b c d (k + 1)) (h0 : PSomosAt b c d k) :
    PSomosAt b c d (2 * k + 1) := by
  simp only [PSomosAt] at h1 h0 ⊢
  rw [show k + 1 + 2 = k + 3 by ring, show k + 1 - 2 = k - 1 by ring,
    show k + 1 + 1 = k + 2 by ring, show k + 1 - 1 = k by ring] at h1
  rw [show 2 * k + 1 + 2 = 2 * (k + 1) + 1 by ring,
    show 2 * k + 1 - 2 = 2 * (k - 1) + 1 by ring,
    show 2 * k + 1 + 1 = 2 * (k + 1) by ring, show 2 * k + 1 - 1 = 2 * k by ring,
    if_neg k.not_even_two_mul_add_one]
  have Op := preNormEDS_odd b c d (k + 1)
  have Om := preNormEDS_odd b c d (k - 1)
  have O0 := preNormEDS_odd b c d k
  have E1 := preNormEDS_even b c d (k + 1)
  have E0 := preNormEDS_even b c d k
  rw [show k + 1 + 2 = k + 3 by ring, show k + 1 - 1 = k by ring,
    show k + 1 + 1 = k + 2 by ring] at Op
  rw [show k - 1 + 2 = k + 1 by ring, show k - 1 - 1 = k - 2 by ring,
    show k - 1 + 1 = k by ring] at Om
  rw [show k + 1 - 1 = k by ring, show k + 1 + 2 = k + 3 by ring,
    show k + 1 - 2 = k - 1 by ring, show k + 1 + 1 = k + 2 by ring] at E1
  set nm2 := preNormEDS b c d (k - 2)
  set nm1 := preNormEDS b c d (k - 1)
  set n0 := preNormEDS b c d k
  set np1 := preNormEDS b c d (k + 1)
  set np2 := preNormEDS b c d (k + 2)
  set np3 := preNormEDS b c d (k + 3)
  set S := preNormEDS b c d (2 * k)
  set T := preNormEDS b c d (2 * k + 1)
  set Vm := preNormEDS b c d (2 * (k - 1) + 1)
  by_cases hk : Even k
  · have hk1 : ¬ Even (k + 1) := by rw [Int.even_add_one]; exact not_not.mpr hk
    have hkm : ¬ Even (k - 1) := by rw [Int.even_sub_one]; exact not_not.mpr hk
    simp only [if_pos hk, if_neg hk1, if_neg hkm, mul_one] at h1 h0 Op Om O0
    linear_combination Vm * Op + (np3*np1^3 - n0*np2^3*b) * Om
      + c * (T + (np2*n0^3*b - nm1*np1^3)) * O0
      - b * S * E1 - b * (n0^2*np1*np3 - nm1*np1*np2^2) * E0
      + (nm1^2*np1^4 - nm1*n0^3*np1*np2*b) * h1
      + (-nm1*n0*np1^3*np2*b + n0^4*np2^2*b^2) * h0
  · have hk1 : Even (k + 1) := Int.even_add_one.mpr hk
    have hkm : Even (k - 1) := Int.even_sub_one.mpr hk
    simp only [if_neg hk, if_pos hk1, if_pos hkm, mul_one] at h1 h0 Op Om O0
    linear_combination Vm * Op + (np3*np1^3*b - n0*np2^3) * Om
      + c * (T + (np2*n0^3 - nm1*np1^3*b)) * O0
      - b * S * E1 - b * (n0^2*np1*np3 - nm1*np1*np2^2) * E0
      + (nm1^2*np1^4*b^2 - nm1*n0^3*np1*np2*b) * h1
      + (-nm1*n0*np1^3*np2*b + n0^4*np2^2) * h0

/-! ## ★★★★Somos-4 関係(本体) -/

theorem psomos_nat (n : ℕ) : PSomosAt b c d (n : ℤ) := by
  induction n using normEDSRec with
  | zero => exact_mod_cast psomos_zero b c d
  | one => exact_mod_cast psomos_one b c d
  | two => exact_mod_cast psomos_two b c d
  | three => exact_mod_cast psomos_three b c d
  | four => exact_mod_cast psomos_four b c d
  | even m h1 h2 h3 h4 h5 =>
      have key := psomos_even b c d ((m : ℤ) + 3)
        (by rw [show (m : ℤ) + 3 + 1 = ((m + 4 : ℕ) : ℤ) by push_cast; ring]; exact h4)
        (by rw [show (m : ℤ) + 3 = ((m + 3 : ℕ) : ℤ) by push_cast; ring]; exact h3)
        (by rw [show (m : ℤ) + 3 - 1 = ((m + 2 : ℕ) : ℤ) by push_cast; ring]; exact h2)
      rw [show ((2 * (m + 3) : ℕ) : ℤ) = 2 * ((m : ℤ) + 3) by push_cast; ring]
      exact key
  | odd m h1 h2 h3 h4 =>
      have key := psomos_odd b c d ((m : ℤ) + 2)
        (by rw [show (m : ℤ) + 2 + 1 = ((m + 3 : ℕ) : ℤ) by push_cast; ring]; exact h3)
        (by rw [show (m : ℤ) + 2 = ((m + 2 : ℕ) : ℤ) by push_cast; ring]; exact h2)
      rw [show ((2 * (m + 2) + 1 : ℕ) : ℤ) = 2 * ((m : ℤ) + 2) + 1 by push_cast; ring]
      exact key

/-- ★★★★**Somos-4 関係——任意の可換環、任意の整数で成り立つ**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★★★★★これが §9-392 で「G1 の律速」と判定したものである。 -/
theorem psomos (n : ℤ) : PSomosAt b c d n := by
  rcases le_or_gt 0 n with h | h
  · lift n to ℕ using h; exact psomos_nat b c d n
  · have key := psomos_nat b c d (-n).toNat
    rw [Int.toNat_of_nonneg (by omega)] at key
    simp only [PSomosAt] at key ⊢
    rw [show -n + 2 = -(n - 2) by ring, show -n - 2 = -(n + 2) by ring,
      show -n + 1 = -(n - 1) by ring, show -n - 1 = -(n + 1) by ring] at key
    simp only [preNormEDS_neg, even_neg] at key
    linear_combination key

/-! ## ★★★★★(★)——y 側に要る 3 項関係の実例 -/

/-- ★★★★★**3 項関係の実例 (★)**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★`ψ` に直すと `ψ_{n+3}ψ_{n−1}ψ_n² − ψ_{n+2}ψ_{n−2}ψ_{n+1}² = ψ_{2n+1}ψ₂²`
——`IsEllSequence` の `(m,n,r) = (n+1, 2, n)` の場合である。

★★**パリティに依らない**——`preNormEDS` の水準では `if` が消える。 -/
theorem eds_star (n : ℤ) :
    preNormEDS b c d (n + 3) * preNormEDS b c d (n - 1) * preNormEDS b c d n ^ 2
      - preNormEDS b c d (n + 2) * preNormEDS b c d (n - 2) * preNormEDS b c d (n + 1) ^ 2
      = preNormEDS b c d (2 * n + 1) := by
  have s1 := psomos b c d (n + 1)
  have s0 := psomos b c d n
  have o := preNormEDS_odd b c d n
  simp only [PSomosAt] at s1 s0
  rw [show n + 1 + 2 = n + 3 by ring, show n + 1 - 2 = n - 1 by ring,
    show n + 1 + 1 = n + 2 by ring, show n + 1 - 1 = n by ring] at s1
  by_cases hn : Even n
  · have hn1 : ¬ Even (n + 1) := by rw [Int.even_add_one]; exact not_not.mpr hn
    simp only [if_pos hn, if_neg hn1, mul_one] at s1 s0 o
    linear_combination (preNormEDS b c d n ^ 2) * s1
      - (preNormEDS b c d (n + 1) ^ 2) * s0 - o
  · have hn1 : Even (n + 1) := Int.even_add_one.mpr hn
    simp only [if_neg hn, if_pos hn1, mul_one] at s1 s0 o
    linear_combination (preNormEDS b c d n ^ 2) * s1
      - (preNormEDS b c d (n + 1) ^ 2) * s0 - o

/-! ## ★★`normEDS` 版(mathlib の TODO の一部) -/

/-- ★★**`normEDS` の Somos-4 関係**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★mathlib の `IsEllSequence` は `(m,n,r)` の 3 変数だが、
これはその `(m, 2, 1)` の場合にあたる。 -/
theorem normEDS_somos (n : ℤ) :
    normEDS b c d (n + 2) * normEDS b c d (n - 2)
      = normEDS b c d (n + 1) * normEDS b c d (n - 1) * b ^ 2
        - c * normEDS b c d n ^ 2 := by
  have h := psomos (b ^ 4) c d n
  simp only [PSomosAt] at h
  simp only [normEDS]
  have e2 : Even (n + 2) ↔ Even n := by simp
  have em2 : Even (n - 2) ↔ Even n := by simp
  have e1 : Even (n + 1) ↔ ¬ Even n := Int.even_add_one
  have em1 : Even (n - 1) ↔ ¬ Even n := Int.even_sub_one
  by_cases hn : Even n
  · rw [if_pos hn, if_pos (e2.mpr hn), if_pos (em2.mpr hn),
      if_neg (fun hh => (e1.mp hh) hn), if_neg (fun hh => (em1.mp hh) hn)] at *
    linear_combination b ^ 2 * h
  · rw [if_neg hn, if_neg (fun hh => hn (e2.mp hh)), if_neg (fun hh => hn (em2.mp hh)),
      if_pos (e1.mpr hn), if_pos (em1.mpr hn)] at *
    linear_combination h

/-! ## ★出典の紐付け(`.src`) -/

def eds_star.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(EDS の Somos-4 関係と 3 項関係の実例)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
