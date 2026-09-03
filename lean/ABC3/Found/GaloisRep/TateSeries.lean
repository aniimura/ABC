import ABC3.Meta.Claim
import Mathlib.NumberTheory.ArithmeticFunction.Misc
import Mathlib.RingTheory.PowerSeries.Basic
import Mathlib.Data.ZMod.Basic
import Mathlib.AlgebraicGeometry.EllipticCurve.Weierstrass

/-!
# Galois (G6) 第 94 ブロック —— **★★★Tate 曲線の q 展開**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★Tate 曲線は**形式冪級数**として作れる——収束は要らない

古典的な `s_k(q) = ∑_{n≥1} n^k q^n/(1−q^n)` は、展開すると

    s_k(q) = ∑_{N≥1} σ_k(N) q^N        (σ_k(N) = ∑_{d ∣ N} d^k)

★**係数が約数和で書ける**ので、`ℤ⟦q⟧` の元として**収束の議論なしに**定義できる。

★★Tate 曲線は

    E_q : y² + xy = x³ + a₄(q) x + a₆(q),
    a₄(q) = −5 s₃(q),   a₆(q) = −(5 s₃(q) + 7 s₅(q))/12

★★★`a₆` が**整数係数**であることは `12 ∣ 5σ₃(N) + 7σ₅(N)` から出る
——これは `∀ d : ZMod 12, 5d³ + 7d⁵ = 0` を `decide` で確かめれば十分である。

## ★★これは (G6) の臨界路の第 1 層である

(G6) が要求するのは **Tate 一意化** `E(K) ≅ Kˣ/q^ℤ` である。
★その前に `E_q` そのものを作らねばならない——本ブロックがそれである。
★★残るのは (i) 完備 DVR への特殊化、(ii) `j(q)` の反転、(iii) 一意化定理。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `twelve_dvd_sigma_comb` | ★★`12 ∣ 5σ₃(n) + 7σ₅(n)` |
| `sigmaSeries` | ★`s_k(q) ∈ ℤ⟦q⟧` |
| `tateA4` / `tateA6` | ★★Tate 曲線の係数(整数係数) |
| `tateCurve` | ★★★**Tate 曲線 `E_q` over `ℤ⟦q⟧`** |
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

/-! ## ★★★★★★★★★★★★`c₄ = E₄`・`c₆ = −E₆` -/

/-- ★★★★★★★★★★★★**Tate 曲線の `c₄` は `E₄`**——★**無条件**（第 1254）。

原文 (GenEll p.15):
> parameter qE of E satisfies the relation qE = qEl ; in particular, we have

☆`a₁ = 1`・`a₂ = a₃ = 0`・`a₄ = −5 σ₃` なので
`c₄ = b₂² − 24 b₄ = 1 − 48·(−5 σ₃) = 1 + 240 σ₃`。 -/
theorem c₄_tateCurve :
    tateCurve.c₄ = 1 + 240 * sigmaSeries 3 := by
  simp only [WeierstrassCurve.c₄, WeierstrassCurve.b₂, WeierstrassCurve.b₄,
    tateCurve, tateA4, map_neg, map_ofNat]
  ring

/-- ★★★★★★★★★★★★**Tate 曲線の `c₆` は `−E₆`**——★**無条件**（第 1254）。

☆`c₆ = −b₂³ + 36 b₂ b₄ − 216 b₆ = −1 + 72 a₄ − 864 a₆`
であり、`12 a₆ = −(5σ₃ + 7σ₅)` から `−1 + 504 σ₅` になる。 -/
theorem c₆_tateCurve :
    tateCurve.c₆ = -1 + 504 * sigmaSeries 5 := by
  have h12 : (12 : PowerSeries ℤ) * tateA6
      = -(5 * sigmaSeries 3 + 7 * sigmaSeries 5) := by
    have h := twelve_mul_tateA6
    simpa [map_ofNat] using h
  have hc : tateCurve.c₆ = -1 + 72 * tateA4 - 72 * ((12 : PowerSeries ℤ) * tateA6) := by
    simp only [WeierstrassCurve.c₆, WeierstrassCurve.b₂, WeierstrassCurve.b₄,
      WeierstrassCurve.b₆, tateCurve]
    ring
  rw [hc, h12, tateA4]
  simp only [map_neg, map_ofNat]
  ring

/-! ## ★出典の紐付け(`.src`) -/

def c₄_tateCurve.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(Tate 曲線の c₄ は E₄ = 1 + 240σ₃。★無条件)",
    sectionId := "genell-lemma-3-2" }

def c₆_tateCurve.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(Tate 曲線の c₆ は −E₆ = −1 + 504σ₅。★無条件)",
    sectionId := "genell-lemma-3-2" }

def tateCurve.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 曲線の q 展開による構成)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
