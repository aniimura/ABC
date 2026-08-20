import ABC3.Found.GaloisRep.CountPow
import Mathlib.RingTheory.DedekindDomain.AdicValuation

/-!
# Galois (G5) 第 143 ブロック —— **★★★★★★D1' 場合 A**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★場合 A が閉じた

D1'(`∀ v, n ∣ count_v((μ f_P))`)は 2 つの場合に分かれる。
★本ブロックは**アフィンの場合(`[n]Q` がアフィン点、すなわち `w(μ ·) ≤ 1`)**を閉じる:

    (a, b)^n = (f)   かつ   ∀ r, w(μ r) ≤ 1
      ⟹  w(μ f) = max(w(μ a), w(μ b))^n
      ⟹  n ∣ count_v((μ f))

★★**分岐指数も素点の対応も使わない。**「`[n]` の下にある素点は何か」を
一切問わずに済むのが要点である。

## ★★★★機構——2 つの不等式

`M := max(w(μ a), w(μ b))` と置く。

* `≤`: `{r | w(μ r) ≤ γ}` は**イデアル**である(`w(μ ·) ≤ 1` だから `R` 倍で閉じる)。
  `(a,b) ≤ levelIdeal M` から帰納法で `(a,b)^n ≤ levelIdeal (M^n)`、
  `f ∈ (a,b)^n` なので `w(μ f) ≤ M^n`。
* `≥`: `a^n ∈ (f)` だから `a^n = s·f` で `w(μ a)^n = w(μ s)·w(μ f) ≤ w(μ f)`。
  `b` も同様なので `M^n ≤ w(μ f)`。

## ★★★★★付値と `count` の橋

mathlib は `FractionalIdeal.count`(整数値)と
`HeightOneSpectrum.valuation`(乗法的)を**別々に**持っており、
両者を結ぶ補題は無かった(2026-08-20 実測)。★本ブロックで橋を架ける:

    v(x) = exp(−count_v((x)))          (x ≠ 0)

★★これにより**超距離不等式が `Valuation` API から無料で使える**ようになり、
結論は `count` の言葉(整数、加法的)で述べられる。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `intValuation_eq_exp_neg_count` | ★★★整元での橋 |
| `valuation_eq_exp_neg_count` | ★★★★★**付値と `count` の橋** |
| `levelIdeal` | ★★`{r | w(μ r) ≤ γ}` はイデアル |
| `valuation_pow_of_span_pair_pow` | ★★★★★**`w(μ f) = max(w(μa), w(μb))^n`** |
| `dvd_count_of_span_pair_pow` | ★★★★★★**場合 A —— 指数は `n` で割れる** |
-/

namespace ABC3.Found.GaloisRep

open IsDedekindDomain nonZeroDivisors

/-! ## ★★★★★付値と `count` の橋 -/

variable {R K : Type} [CommRing R] [IsDedekindDomain R] [Field K] [Algebra R K]
  [IsFractionRing R K]

/-- ★★★整元での橋——`v.intValuation r = exp(−count_v((r)))`。 -/
theorem intValuation_eq_exp_neg_count (v : HeightOneSpectrum R) {r : R} (hr : r ≠ 0) :
    v.intValuation r
      = WithZero.exp (- FractionalIdeal.count K v
          (FractionalIdeal.spanSingleton R⁰ (algebraMap R K r))) := by
  classical
  rw [← FractionalIdeal.coeIdeal_span_singleton,
    FractionalIdeal.count_coe K v (by simpa [Ideal.zero_eq_bot, Ideal.span_singleton_eq_bot]
      using hr)]
  exact v.intValuationDef_if_neg hr

/-- ★★★★★**付値と `count` の橋**——`v(x) = exp(−count_v((x)))`。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★mathlib は両者を別々に持っており、繋ぐ補題は無かった(2026-08-20 実測)。
★★これで超距離不等式が `Valuation` API から無料になる。 -/
theorem valuation_eq_exp_neg_count (v : HeightOneSpectrum R) {x : K} (hx : x ≠ 0) :
    v.valuation K x
      = WithZero.exp (- FractionalIdeal.count K v (FractionalIdeal.spanSingleton R⁰ x)) := by
  have hinj : Function.Injective (algebraMap R K) := IsFractionRing.injective R K
  obtain ⟨a, b, hb, rfl⟩ := IsFractionRing.div_surjective (A := R) x
  have hb0 : b ≠ 0 := nonZeroDivisors.ne_zero hb
  have hbK : algebraMap R K b ≠ 0 := fun h => hb0 (hinj (by rw [h, map_zero]))
  have ha0 : a ≠ 0 := by
    intro h0
    rw [h0, map_zero, zero_div] at hx
    exact hx rfl
  have haK : algebraMap R K a ≠ 0 := fun h => ha0 (hinj (by rw [h, map_zero]))
  have hsplit : FractionalIdeal.spanSingleton R⁰ (algebraMap R K a / algebraMap R K b)
      = FractionalIdeal.spanSingleton R⁰ (algebraMap R K a)
        * (FractionalIdeal.spanSingleton R⁰ (algebraMap R K b))⁻¹ := by
    rw [FractionalIdeal.spanSingleton_inv, FractionalIdeal.spanSingleton_mul_spanSingleton,
      div_eq_mul_inv]
  rw [hsplit, FractionalIdeal.count_mul K v (FractionalIdeal.spanSingleton_ne_zero_iff.2 haK)
      (by rw [Ne, inv_eq_zero]; exact FractionalIdeal.spanSingleton_ne_zero_iff.2 hbK),
    FractionalIdeal.count_inv, map_div₀, v.valuation_of_algebraMap, v.valuation_of_algebraMap,
    intValuation_eq_exp_neg_count (K := K) v ha0, intValuation_eq_exp_neg_count (K := K) v hb0,
    ← WithZero.exp_sub]
  congr 1
  ring

/-! ## ★★水準イデアル -/

section Level

variable {A L Γ₀ : Type} [CommRing A] [Field L] [LinearOrderedCommGroupWithZero Γ₀]

/-- ★★`{r | w(μ r) ≤ γ}` はイデアルである(`w(μ ·) ≤ 1` のとき)。 -/
def levelIdeal (w : Valuation L Γ₀) (μ : A →+* L) (hnn : ∀ r : A, w (μ r) ≤ 1) (γ : Γ₀) :
    Ideal A where
  carrier := {r | w (μ r) ≤ γ}
  zero_mem' := by simp
  add_mem' := by
    intro x y hx hy
    simp only [Set.mem_setOf_eq] at *
    have hμ : μ (x + y) = μ x + μ y := map_add μ x y
    rw [hμ]
    exact le_trans (Valuation.map_add w _ _) (max_le hx hy)
  smul_mem' := by
    intro c x hx
    simp only [Set.mem_setOf_eq] at *
    have hμ : μ (c • x) = μ c * μ x := by rw [smul_eq_mul, map_mul]
    rw [hμ, Valuation.map_mul]
    calc w (μ c) * w (μ x) ≤ 1 * γ := mul_le_mul' (hnn c) hx
      _ = γ := one_mul γ

theorem mem_levelIdeal {w : Valuation L Γ₀} {μ : A →+* L} {hnn : ∀ r : A, w (μ r) ≤ 1} {γ : Γ₀}
    {r : A} : r ∈ levelIdeal w μ hnn γ ↔ w (μ r) ≤ γ := Iff.rfl

/-- ★★★★★**`(a,b)^n = (f)` なら `w(μ f) = max(w(μa), w(μb))^n`**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★`≤` は水準イデアルへの帰納法、`≥` は `a^n, b^n ∈ (f)` から。
★★**分岐指数も素点の対応も使わない。** -/
theorem valuation_pow_of_span_pair_pow (w : Valuation L Γ₀) (μ : A →+* L)
    (hnn : ∀ r : A, w (μ r) ≤ 1) {a b f : A} {n : ℕ}
    (hf : (Ideal.span ({a, b} : Set A)) ^ n = Ideal.span {f}) :
    w (μ f) = (max (w (μ a)) (w (μ b))) ^ n := by
  set M := max (w (μ a)) (w (μ b)) with hM
  set I : Ideal A := Ideal.span ({a, b} : Set A) with hI
  have hIlev : I ≤ levelIdeal w μ hnn M := by
    rw [hI, Ideal.span_le]
    rintro z (rfl | rfl)
    · exact le_max_left _ _
    · exact le_max_right _ _
  have hpow : ∀ k : ℕ, I ^ k ≤ levelIdeal w μ hnn (M ^ k) := by
    intro k
    induction k with
    | zero =>
      intro z _
      show w (μ z) ≤ M ^ 0
      rw [pow_zero]
      exact hnn z
    | succ k ih =>
      rw [pow_succ]
      refine Ideal.mul_le.2 fun r hr s hs => ?_
      show w (μ (r * s)) ≤ M ^ (k + 1)
      rw [map_mul, Valuation.map_mul, pow_succ]
      exact mul_le_mul' (ih hr) (hIlev hs)
  have hle : w (μ f) ≤ M ^ n := by
    refine hpow n ?_
    rw [hf]
    exact Ideal.subset_span rfl
  have hmem : ∀ z : A, z ∈ I → w (μ z) ^ n ≤ w (μ f) := by
    intro z hz
    have hzn : z ^ n ∈ Ideal.span ({f} : Set A) := by
      rw [← hf]
      exact Ideal.pow_mem_pow hz n
    obtain ⟨s, hs⟩ := Ideal.mem_span_singleton'.1 hzn
    have hws : w (μ (z ^ n)) = w (μ s) * w (μ f) := by rw [← hs, map_mul, Valuation.map_mul]
    rw [map_pow, Valuation.map_pow] at hws
    rw [hws]
    calc w (μ s) * w (μ f) ≤ 1 * w (μ f) := mul_le_mul' (hnn s) le_rfl
      _ = w (μ f) := one_mul _
  have hge : M ^ n ≤ w (μ f) := by
    rcases max_cases (w (μ a)) (w (μ b)) with ⟨hmax, _⟩ | ⟨hmax, _⟩
    · rw [hM, hmax]; exact hmem a (Ideal.subset_span (by simp))
    · rw [hM, hmax]; exact hmem b (Ideal.subset_span (by simp))
  exact le_antisymm hle hge

end Level

/-! ## ★★★★★★場合 A の結論 -/

/-- ★★★★★★**場合 A —— 指数は `n` で割れる**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★`w(μ ·) ≤ 1`(引き戻しが `v` で極を持たない)ときの結論である。
★★橋(`valuation_eq_exp_neg_count`)で `count` の言葉に移すだけ。 -/
theorem dvd_count_of_span_pair_pow (v : HeightOneSpectrum R) (μ : R →+* K)
    (hμinj : Function.Injective μ)
    (hnn : ∀ r : R, v.valuation K (μ r) ≤ 1)
    {a b f : R} (ha : a ≠ 0) (hb : b ≠ 0) (hf0 : f ≠ 0) {n : ℕ}
    (hf : (Ideal.span ({a, b} : Set R)) ^ n = Ideal.span {f}) :
    (n : ℤ) ∣ FractionalIdeal.count K v (FractionalIdeal.spanSingleton R⁰ (μ f)) := by
  have hne : ∀ {r : R}, r ≠ 0 → μ r ≠ 0 := fun {r} hr h => hr (hμinj (by rw [h, map_zero]))
  set ca := FractionalIdeal.count K v (FractionalIdeal.spanSingleton R⁰ (μ a)) with hca
  set cb := FractionalIdeal.count K v (FractionalIdeal.spanSingleton R⁰ (μ b)) with hcb
  set cf := FractionalIdeal.count K v (FractionalIdeal.spanSingleton R⁰ (μ f)) with hcf
  have hva := valuation_eq_exp_neg_count (K := K) v (hne ha)
  have hvb := valuation_eq_exp_neg_count (K := K) v (hne hb)
  have hvf := valuation_eq_exp_neg_count (K := K) v (hne hf0)
  have hkey := valuation_pow_of_span_pair_pow (v.valuation K) μ hnn hf
  rw [hva, hvb, hvf] at hkey
  have hmax : max (WithZero.exp (-ca)) (WithZero.exp (-cb))
      = WithZero.exp (-(min ca cb)) := by
    rcases le_total ca cb with hle | hle
    · rw [min_eq_left hle, max_eq_left (WithZero.exp_le_exp.2 (by omega))]
    · rw [min_eq_right hle, max_eq_right (WithZero.exp_le_exp.2 (by omega))]
  rw [hmax, ← WithZero.exp_nsmul, WithZero.exp_inj] at hkey
  refine ⟨min ca cb, ?_⟩
  have h2 : -cf = (n : ℤ) * (-(min ca cb)) := by rw [← nsmul_eq_mul]; exact hkey
  rw [mul_neg] at h2
  exact neg_inj.1 h2

/-! ## ★出典の紐付け(`.src`) -/

def dvd_count_of_span_pair_pow.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——引き戻しが極を持たない素点での指数の n 可除性)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
