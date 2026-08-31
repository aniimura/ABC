/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import Mathlib.NumberTheory.Primorial
import Mathlib.NumberTheory.NumberField.Discriminant.Basic
import ABC3.Meta.Claim

/-!
# `θ(d) ≤ 2·log|disc L|`（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.23。

原文 (GenEll p.22):
> Corollary 4.3. (Full Galois Actions for Degenerating Elliptic Curves)

## ★★★★★★★★これは何か

`EllModuliData` の `sum_log_ramPrimes_le` 欄

    ∑_{p ∈ ramPrimes} log p ≤ ∑_{p ∈ badPrimes} log p + 3·d·log-diff([E_L])

を埋めるのに要る**数論的な核**である。`ramPrimes` には

* `badPrimes`
* `L` で分岐する素数（＝ `disc L` を割る素数）
* `L/ℚ` の分岐指数を割る素数（分岐指数は `≤ d` なので `d` 以下の素数すべて）

が入る。★`log-diff = log|disc L| / d` と取ると、右辺の余りは `3·log|disc L|` であり、

* 分岐する素数の分は `∑_{p ∣ disc} log p ≤ log|disc|`（根基は本体以下）で **1 単位**
* `d` 以下の素数の分は本ファイルの `θ(d) ≤ 2·log|disc|` で **2 単位**

でちょうど収まる。★★係数 `3` が原文の値と一致するのは偶然ではない。

## ★★★機構——2 つの mathlib の定理

* `primorial_le_four_pow`（Chebyshev）: `d# ≤ 4^d`
* `NumberField.abs_discr_ge`（Hermite–Minkowski）: `(4/9)(3π/4)^d ≤ |disc|`（`d ≥ 2`）

★`4^d ≤ |disc|²` を示せばよく、それは `2^d ≤ |disc|`、すなわち `9/4 ≤ (3π/8)^d` に帰着する。
`3π/8 > 1.178` なので `d ≥ 5` で成り立つ（`1.178⁵ = 2.2636 > 2.25`）。
☆`d ≤ 4` は `d# ≤ 6` と `|disc| ≥ 3` で処理する。
-/

namespace ABC3.Found.GaloisRep

open NumberField

/-! ## ★★★数値の核 -/

/-- ★★★**`d ≥ 5` なら `9/4 ≤ (3π/8)^d`**。 -/
theorem nine_div_four_le_pow (d : ℕ) (hd : 5 ≤ d) : (9:ℝ)/4 ≤ (3 * Real.pi / 8) ^ d := by
  have hpi : (3.141592 : ℝ) < Real.pi := Real.pi_gt_d6
  have hb : (1.178 : ℝ) ≤ 3 * Real.pi / 8 := by linarith
  have h5 : (9:ℝ)/4 ≤ (1.178 : ℝ) ^ 5 := by norm_num
  have hmono : ((1.178 : ℝ)) ^ 5 ≤ (3 * Real.pi / 8) ^ 5 :=
    pow_le_pow_left₀ (by norm_num) hb 5
  have hgrow : (3 * Real.pi / 8) ^ 5 ≤ (3 * Real.pi / 8) ^ d :=
    pow_le_pow_right₀ (by linarith) hd
  linarith

/-! ## ★★★★★★`d# ≤ |disc|²` -/

variable (L : Type) [Field L] [NumberField L]

/-- ★★★★★★**`d# ≤ |disc L|²`**——Chebyshev と Hermite–Minkowski。 -/
theorem primorial_le_discr_sq :
    (primorial (Module.finrank ℚ L) : ℝ) ≤ ((|NumberField.discr L| : ℤ) : ℝ) ^ 2 := by
  set d := Module.finrank ℚ L with hd
  have hdpos : 0 < d := Module.finrank_pos
  have hD0 : (1:ℝ) ≤ ((|NumberField.discr L| : ℤ) : ℝ) := by
    have : (1:ℤ) ≤ |NumberField.discr L| :=
      Int.one_le_abs (NumberField.discr_ne_zero (K := L))
    exact_mod_cast this
  rcases Nat.lt_or_ge d 5 with hgt | hd5
  · -- `d ≤ 4`: `d# ≤ 6` と `|disc| ≥ 3`（`d ≥ 2`）／`d = 1` は `d# = 1`
    have hle : d ≤ 4 := by omega
    rcases Nat.lt_or_ge d 2 with h1 | h2
    · -- `d = 1`
      have hd1 : d = 1 := by omega
      rw [hd1]
      have hp1 : (primorial 1 : ℝ) = 1 := by rw [primorial_one]; norm_num
      rw [hp1]
      nlinarith
    · -- `2 ≤ d ≤ 4`
      have hd2 : 1 < d := by omega
      have hD3 : (3:ℝ) ≤ ((|NumberField.discr L| : ℤ) : ℝ) := by
        have : (2:ℤ) < |NumberField.discr L| := NumberField.abs_discr_gt_two (K := L) (by omega)
        have : (3:ℤ) ≤ |NumberField.discr L| := by omega
        exact_mod_cast this
      have hp6 : (primorial d : ℝ) ≤ 6 := by
        have hnat : primorial d ≤ 6 := by interval_cases d <;> decide
        exact_mod_cast hnat
      nlinarith
  · -- `d ≥ 5`
    have hpi : (0:ℝ) < 3 * Real.pi / 8 := by positivity
    have hkey : (9:ℝ)/4 ≤ (3 * Real.pi / 8) ^ d := nine_div_four_le_pow d hd5
    have hpow : ((2:ℝ)) ^ d ≤ (4/9 : ℝ) * (3 * Real.pi / 4) ^ d := by
      have hsplit : (3 * Real.pi / 4 : ℝ) ^ d = (2:ℝ) ^ d * (3 * Real.pi / 8) ^ d := by
        rw [← mul_pow]
        ring_nf
      rw [hsplit]
      have h2d : (0:ℝ) < (2:ℝ) ^ d := by positivity
      nlinarith
    have hmink : (4/9 : ℝ) * (3 * Real.pi / 4) ^ d ≤ ((|NumberField.discr L| : ℤ) : ℝ) := by
      have := NumberField.abs_discr_ge (K := L) (by omega)
      exact_mod_cast this
    have h2D : ((2:ℝ)) ^ d ≤ ((|NumberField.discr L| : ℤ) : ℝ) := le_trans hpow hmink
    have hcheb : (primorial d : ℝ) ≤ (4:ℝ) ^ d := by
      have := primorial_le_four_pow d
      exact_mod_cast this
    have h4 : (4:ℝ) ^ d = ((2:ℝ) ^ d) ^ 2 := by
      rw [← pow_mul, mul_comm d 2, pow_mul]
      norm_num
    have h2pos : (0:ℝ) ≤ (2:ℝ) ^ d := by positivity
    calc (primorial d : ℝ) ≤ (4:ℝ) ^ d := hcheb
      _ = ((2:ℝ) ^ d) ^ 2 := h4
      _ ≤ ((|NumberField.discr L| : ℤ) : ℝ) ^ 2 := by
          exact pow_le_pow_left₀ h2pos h2D 2

/-- ★★★★★★★★**`θ(d) ≤ 2·log|disc L|`**——★**無条件**。 -/
theorem log_primorial_le_two_log_discr :
    Real.log (primorial (Module.finrank ℚ L))
      ≤ 2 * Real.log ((|NumberField.discr L| : ℤ) : ℝ) := by
  have hD0 : (0:ℝ) < ((|NumberField.discr L| : ℤ) : ℝ) := by
    have : (1:ℤ) ≤ |NumberField.discr L| :=
      Int.one_le_abs (NumberField.discr_ne_zero (K := L))
    have : (1:ℝ) ≤ ((|NumberField.discr L| : ℤ) : ℝ) := by exact_mod_cast this
    linarith
  have hp0 : (0:ℝ) < (primorial (Module.finrank ℚ L) : ℝ) := by
    exact_mod_cast primorial_pos _
  have h := primorial_le_discr_sq L
  calc Real.log (primorial (Module.finrank ℚ L))
      ≤ Real.log (((|NumberField.discr L| : ℤ) : ℝ) ^ 2) := Real.log_le_log hp0 h
    _ = 2 * Real.log ((|NumberField.discr L| : ℤ) : ℝ) := by
        rw [Real.log_pow]; push_cast; ring

/-! ## ★出典の紐付け(`.src`) -/

def log_primorial_le_two_log_discr.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 22,
    item := "Corollary 4.3(θ(d) ≤ 2·log|disc L|——分岐指数を割る素数の分。★無条件)",
    sectionId := "genell-cor-4-3" }

def log_primorial_le_two_log_discr.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]" "primorial_le_four_pow(Chebyshev: d# ≤ 4^d)"
      (.inMathlib "primorial_le_four_pow") 1,
    .citation "[mathlib]" "NumberField.abs_discr_ge(Hermite–Minkowski)"
      (.inMathlib "NumberField.abs_discr_ge") 1 ]

end ABC3.Found.GaloisRep
