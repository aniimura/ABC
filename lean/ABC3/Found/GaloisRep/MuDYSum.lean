/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.MuPowerSum
import ABC3.Found.GaloisRep.TateDSeries

/-!
# Galois (G6) 第 849 ブロック —— **★★★★★★★★★★`∑_ζ Dg(ζ) = (l⁴−1)/240`**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> parameter qE of E satisfies the relation qE = qEl ; in particular, we have

## ★★★★★★★★★★★★★★★★★★★★葉が一つ閉じた

第 846 の葉 `sum_mu_dyterm`

    `240 · ∑_{ζ≠1} ζ²(2+ζ)/(1−ζ)⁴ = l⁴ − 1`

が**証明済みになった**（体・標数 0 の版）。★道は第 847-848 の機械 `muPow`/`muM` だけである。

`x = 1/(1−ζ)` と置くと

    `ζ²(2+ζ)/(1−ζ)⁴ = 3x⁴ − 7x³ + 5x² − x`

であり（`tateDYterm_eq_powers`）、`p_k = ∑_{ζ≠1}x^k` は第 848 で全部出ている:

    `p₁ = (l−1)/2`、`p₂ = (l−1)(5−l)/12`、`p₃ = (l−1)(3−l)/8`
    `p₄ = (l−1)(l³+l²−109l+251)/720`

よって `240(3p₄ − 7p₃ + 5p₂ − p₁) = (l−1)·3(l³+l²+l+1) = (l²−1)(l²+1) = l⁴ − 1` ✓。

☆`p₄` に「新しい種」が要ると思われていたが、`∑_{a<l}ζ^a = 0` だけで足りた（第 847）。
-/

namespace ABC3.Found.GaloisRep

open Finset

variable {F : Type} [Field F]

/-- ★★★★**`Dg(t) = 3x⁴ − 7x³ + 5x² − x`**（`x = 1/(1−t)`）。

`t = (x−1)/x` を入れると `t²(2+t)x⁴ = (x−1)²(3x−1)x = 3x⁴ − 7x³ + 5x² − x`。 -/
theorem tateDYterm_eq_powers {t : F} (hne : (1 : F) - t ≠ 0) :
    tateDYterm t
      = 3 * ((1 - t)⁻¹) ^ 4 - 7 * ((1 - t)⁻¹) ^ 3 + 5 * ((1 - t)⁻¹) ^ 2 - (1 - t)⁻¹ := by
  rw [tateDYterm, Ring.inverse_eq_inv']
  field_simp
  ring

/-- ★★★★★★**`∑_{ζ≠1} Dg(ζ)` は `p₁…p₄` の一次結合**。 -/
theorem sum_mu_dyterm_eq_muPow {l : ℕ} {ζ : F} (hζ : IsPrimitiveRoot ζ l) :
    ∑ i ∈ (range l).erase 0, tateDYterm (ζ ^ i)
      = 3 * muPow l ζ 4 0 - 7 * muPow l ζ 3 0 + 5 * muPow l ζ 2 0 - muPow l ζ 1 0 := by
  simp only [muPow, pow_zero, one_mul, Finset.mul_sum]
  rw [← Finset.sum_sub_distrib, ← Finset.sum_add_distrib, ← Finset.sum_sub_distrib]
  refine Finset.sum_congr rfl (fun i hi => ?_)
  rw [tateDYterm_eq_powers (one_sub_zeta_ne_zero hζ hi)]
  ring

/-- ★★★★★★★★★★**`240·∑_{ζ≠1} ζ²(2+ζ)/(1−ζ)⁴ = l⁴ − 1`**——第 846 の葉。 -/
theorem sum_mu_dyterm_field [CharZero F] {l : ℕ} (hl : l.Prime) {ζ : F}
    (hζ : IsPrimitiveRoot ζ l) :
    240 * (∑ i ∈ (range l).erase 0, tateDYterm (ζ ^ i)) = (l : F) ^ 4 - 1 := by
  rw [sum_mu_dyterm_eq_muPow hζ, muPow_four_zero hl hζ, muPow_three_zero hl hζ,
    muPow_two_zero hl hζ, muPow_one_zero hl hζ]
  ring

/-! ## ★★★★★★`D²f` と `Df` の指標和 -/

/-- ★★**`D²f(t) = 6x⁴ − 12x³ + 7x² − x`**（`x = 1/(1−t)`）。 -/
theorem tateD2Xterm_eq_powers {t : F} (hne : (1 : F) - t ≠ 0) :
    tateD2Xterm t
      = 6 * ((1 - t)⁻¹) ^ 4 - 12 * ((1 - t)⁻¹) ^ 3 + 7 * ((1 - t)⁻¹) ^ 2
        - (1 - t)⁻¹ := by
  rw [tateD2Xterm, Ring.inverse_eq_inv']
  field_simp
  ring

/-- ★★**`Df(t) = 2x³ − 3x² + x`**（`x = 1/(1−t)`）。 -/
theorem tateDXterm_eq_powers {t : F} (hne : (1 : F) - t ≠ 0) :
    tateDXterm t = 2 * ((1 - t)⁻¹) ^ 3 - 3 * ((1 - t)⁻¹) ^ 2 + (1 - t)⁻¹ := by
  rw [tateDXterm, Ring.inverse_eq_inv']
  field_simp
  ring

/-- ★★★★★★★★**`120·∑_{ζ≠1} D²f(ζ) = l⁴ − 1`**。 -/
theorem sum_mu_d2xterm_field [CharZero F] {l : ℕ} (hl : l.Prime) {ζ : F}
    (hζ : IsPrimitiveRoot ζ l) :
    120 * (∑ i ∈ (range l).erase 0, tateD2Xterm (ζ ^ i)) = (l : F) ^ 4 - 1 := by
  have hsplit : ∑ i ∈ (range l).erase 0, tateD2Xterm (ζ ^ i)
      = 6 * muPow l ζ 4 0 - 12 * muPow l ζ 3 0 + 7 * muPow l ζ 2 0 - muPow l ζ 1 0 := by
    simp only [muPow, pow_zero, one_mul, Finset.mul_sum]
    rw [← Finset.sum_sub_distrib, ← Finset.sum_add_distrib, ← Finset.sum_sub_distrib]
    refine Finset.sum_congr rfl (fun i hi => ?_)
    rw [tateD2Xterm_eq_powers (one_sub_zeta_ne_zero hζ hi)]
    ring
  rw [hsplit, muPow_four_zero hl hζ, muPow_three_zero hl hζ,
    muPow_two_zero hl hζ, muPow_one_zero hl hζ]
  ring

/-- ★★★★★★**`∑_{ζ≠1} Df(ζ) = 0`**——`∑_ζ X(ζ)` は `q` だけの関数なので
微分で消える、という事実の定数項。 -/
theorem sum_mu_dxterm_field [CharZero F] {l : ℕ} (hl : l.Prime) {ζ : F}
    (hζ : IsPrimitiveRoot ζ l) :
    ∑ i ∈ (range l).erase 0, tateDXterm (ζ ^ i) = 0 := by
  have hsplit : ∑ i ∈ (range l).erase 0, tateDXterm (ζ ^ i)
      = 2 * muPow l ζ 3 0 - 3 * muPow l ζ 2 0 + muPow l ζ 1 0 := by
    simp only [muPow, pow_zero, one_mul, Finset.mul_sum]
    rw [← Finset.sum_sub_distrib, ← Finset.sum_add_distrib]
    refine Finset.sum_congr rfl (fun i hi => ?_)
    rw [tateDXterm_eq_powers (one_sub_zeta_ne_zero hζ hi)]
    ring
  rw [hsplit, muPow_three_zero hl hζ, muPow_two_zero hl hζ, muPow_one_zero hl hζ]
  ring

/-! ## ★★★★★★★★整域への転送（分数体に埋め込む） -/

section Domain

variable {A : Type} [CommRing A] [IsDomain A] [CharZero A]

theorem map_ringInverse {B : Type} [CommRing B] (f : A →+* B) {x : A} (hx : IsUnit x) :
    f (Ring.inverse x) = Ring.inverse (f x) := by
  refine (ring_inverse_eq_of_mul_eq_one (hx.map f) ?_).symm
  rw [← map_mul, Ring.mul_inverse_cancel x hx, map_one]

theorem map_tateD2Xterm {B : Type} [CommRing B] (f : A →+* B) {t : A}
    (ht : IsUnit (1 - t)) : f (tateD2Xterm t) = tateD2Xterm (f t) := by
  rw [tateD2Xterm, tateD2Xterm, map_mul, map_mul, map_pow, map_ringInverse f ht]
  simp [map_ofNat]

/-- ★★★★★★★★**整域版**: `120·∑_{ζ≠1} D²f(ζ) = l⁴ − 1`。

★分数体に埋め込んで体版（`sum_mu_d2xterm_field`）を使い、単射性で戻す。 -/
theorem sum_mu_d2xterm {l : ℕ} (hl : l.Prime) {ζ : A} (hζ : IsPrimitiveRoot ζ l)
    (hu : ∀ i ∈ (range l).erase 0, IsUnit (1 - ζ ^ i)) :
    120 * (∑ i ∈ (range l).erase 0, tateD2Xterm (ζ ^ i)) = (l : A) ^ 4 - 1 := by
  have hinj : Function.Injective (algebraMap A (FractionRing A)) :=
    IsFractionRing.injective A (FractionRing A)
  haveI : CharZero (FractionRing A) := charZero_of_injective_algebraMap hinj
  refine hinj ?_
  have hmap : (algebraMap A (FractionRing A))
      (120 * ∑ i ∈ (range l).erase 0, tateD2Xterm (ζ ^ i))
      = 120 * ∑ i ∈ (range l).erase 0,
          tateD2Xterm ((algebraMap A (FractionRing A)) ζ ^ i) := by
    rw [map_mul, map_sum]
    congr 1
    · exact map_ofNat _ 120
    · refine Finset.sum_congr rfl (fun i hi => ?_)
      rw [map_tateD2Xterm _ (hu i hi), map_pow]
  rw [hmap, map_sub, map_pow, map_natCast, map_one]
  exact sum_mu_d2xterm_field hl (hζ.map_of_injective hinj)

end Domain

/-! ## ★出典の紐付け(`.src`) -/

def tateDYterm_eq_powers.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(Dg(t) = 3x⁴ − 7x³ + 5x² − x。★1−t ≠ 0)",
    sectionId := "genell-lemma-3-2" }

def sum_mu_dyterm_field.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(240·∑ ζ²(2+ζ)/(1−ζ)⁴ = l⁴ − 1。★体・標数 0)",
    sectionId := "genell-lemma-3-2" }

end ABC3.Found.GaloisRep
