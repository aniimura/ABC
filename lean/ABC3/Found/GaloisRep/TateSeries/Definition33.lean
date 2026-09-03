/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Meta.Claim
import Mathlib.NumberTheory.ArithmeticFunction.Misc
import Mathlib.RingTheory.PowerSeries.Basic
import Mathlib.Data.ZMod.Basic
import Mathlib.AlgebraicGeometry.EllipticCurve.Weierstrass

/-!
# TateSeries —— `[GenEll] Definition 3.3` の分

☆もとの 1 枚を**条なしの項目ごと**に割ったものである（第 1458、案 a）。
★「1 ファイル = 1 ノード」を回復するための分割で、中身は動かしていない。
-/

namespace ABC3.Found.GaloisRep



/-! ## ★★整数性 -/

/-- ★★**`12 ∣ 5σ₃(n) + 7σ₅(n)`**——`a₆(q)` が整数係数であることの中身。 -/
theorem twelve_dvd_sigma_comb (n : ℕ) :
    (12 : ℤ) ∣ 5 * (ArithmeticFunction.sigma 3 n : ℤ)
      + 7 * (ArithmeticFunction.sigma 5 n : ℤ) := by
  have hterm : ∀ d : ℕ, (12 : ℤ) ∣ 5 * (d : ℤ) ^ 3 + 7 * (d : ℤ) ^ 5 := by
    intro d
    have hz : ((5 * (d : ℤ) ^ 3 + 7 * (d : ℤ) ^ 5 : ℤ) : ZMod 12) = 0 := by
      push_cast
      have hall : ∀ x : ZMod 12, 5 * x ^ 3 + 7 * x ^ 5 = 0 := by decide
      exact hall _
    exact_mod_cast (ZMod.intCast_zmod_eq_zero_iff_dvd _ 12).1 hz
  have hsum : 5 * (ArithmeticFunction.sigma 3 n : ℤ) + 7 * (ArithmeticFunction.sigma 5 n : ℤ)
      = ∑ d ∈ n.divisors, (5 * (d : ℤ) ^ 3 + 7 * (d : ℤ) ^ 5) := by
    rw [ArithmeticFunction.sigma_apply, ArithmeticFunction.sigma_apply]
    push_cast
    rw [Finset.mul_sum, Finset.mul_sum, ← Finset.sum_add_distrib]
  rw [hsum]
  exact Finset.dvd_sum (fun d _ => hterm d)

/-! ## ★q 展開 -/

/-- ★**`s_k(q) = ∑_{N≥1} σ_k(N) q^N`**。 -/
def sigmaSeries (k : ℕ) : PowerSeries ℤ :=
  PowerSeries.mk (fun N => if N = 0 then 0 else (ArithmeticFunction.sigma k N : ℤ))

@[simp] theorem coeff_sigmaSeries (k N : ℕ) :
    PowerSeries.coeff N (sigmaSeries k)
      = if N = 0 then 0 else (ArithmeticFunction.sigma k N : ℤ) := by
  rw [sigmaSeries, PowerSeries.coeff_mk]

/-- ★★Tate 曲線の係数 `a₄(q) = −5 s₃(q)`。 -/
noncomputable def tateA4 : PowerSeries ℤ := PowerSeries.C (-5 : ℤ) * sigmaSeries 3

/-- ★★Tate 曲線の係数 `a₆(q) = −(5 s₃(q) + 7 s₅(q))/12`(整数係数)。 -/
noncomputable def tateA6 : PowerSeries ℤ :=
  PowerSeries.mk (fun N => if N = 0 then 0 else
    -(5 * (ArithmeticFunction.sigma 3 N : ℤ) + 7 * (ArithmeticFunction.sigma 5 N : ℤ)) / 12)

/-- ★★★**`12 a₆(q) = −(5 s₃(q) + 7 s₅(q))`**——整数性の帰結。 -/
theorem twelve_mul_tateA6 :
    PowerSeries.C (12 : ℤ) * tateA6
      = -(PowerSeries.C (5 : ℤ) * sigmaSeries 3 + PowerSeries.C (7 : ℤ) * sigmaSeries 5) := by
  ext N
  rw [PowerSeries.coeff_C_mul, map_neg, map_add, PowerSeries.coeff_C_mul,
    PowerSeries.coeff_C_mul, coeff_sigmaSeries, coeff_sigmaSeries, tateA6, PowerSeries.coeff_mk]
  by_cases hN : N = 0
  · subst hN
    simp
  · rw [if_neg hN, if_neg hN, if_neg hN]
    have hdvd := twelve_dvd_sigma_comb N
    have h12 : (12 : ℤ) * (-(5 * (ArithmeticFunction.sigma 3 N : ℤ)
        + 7 * (ArithmeticFunction.sigma 5 N : ℤ)) / 12)
        = -(5 * (ArithmeticFunction.sigma 3 N : ℤ) + 7 * (ArithmeticFunction.sigma 5 N : ℤ)) :=
      Int.mul_ediv_cancel' (dvd_neg.2 hdvd)
    rw [h12]

/-! ## ★★★Tate 曲線 -/

/-- ★★★**Tate 曲線 `E_q : y² + xy = x³ + a₄(q)x + a₆(q)`**(`ℤ⟦q⟧` 上)。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
noncomputable def tateCurve : WeierstrassCurve (PowerSeries ℤ) where
  a₁ := 1
  a₂ := 0
  a₃ := 0
  a₄ := tateA4
  a₆ := tateA6

def tateCurve.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 曲線の q 展開による構成)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
