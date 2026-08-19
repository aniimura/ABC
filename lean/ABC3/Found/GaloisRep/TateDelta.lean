import ABC3.Found.GaloisRep.TateSpecialize

/-!
# Galois (G6) 第 97 ブロック —— **★★★`Δ(E_q) = q + O(q²)`**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★これが Tate 曲線の非退化性である

Tate 曲線 `E_q : y² + xy = x³ + a₄(q)x + a₆(q)` の判別式は、古典的には

    Δ(q) = q ∏_{n≥1} (1 − qⁿ)^24

★★本ブロックはその**低次 2 項**だけを取る:

    coeff₀ Δ = 0,   coeff₁ Δ = 1        すなわち `Δ = q + O(q²)`

★★★これで十分である——`Δ` の**位数がちょうど 1** であることが
(G6) の「局所高さ `v(q)` が正」の根拠であり、
★また `j = c₄³/Δ` の**極**がここから出る。

## ★★機構

`a₄ = −5s₃`、`a₆ = −(5s₃+7s₅)/12` の低次係数は `σ_k(1) = 1` から出る:

    coeff₀ a₄ = 0,  coeff₁ a₄ = −5,   coeff₀ a₆ = 0,  coeff₁ a₆ = −1

★`Δ = −a₆ + a₄² − 64a₄³ − 432a₆² + 72a₄a₆` に入れると、
定数項が 0 の項の積はすべて 2 次以上に落ちるので `coeff₁ Δ = −coeff₁ a₆ = 1`。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `coeff_zero_sigmaSeries` / `coeff_one_sigmaSeries` | ★`s_k = q + O(q²)` |
| `coeff_*_tateA4` / `coeff_*_tateA6` | ★★係数の低次 |
| `tateCurve_Delta_eq` | ★★`Δ` の明示式 |
| `coeff_zero_tateDelta` / `coeff_one_tateDelta` | ★★★**`Δ = q + O(q²)`** |
| `tateCurveAt_Delta` | ★特殊化との両立 |
-/

namespace ABC3.Found.GaloisRep

open PowerSeries

/-! ## ★低次係数 -/

theorem coeff_zero_sigmaSeries (k : ℕ) : PowerSeries.coeff 0 (sigmaSeries k) = 0 := by
  rw [coeff_sigmaSeries]
  simp

theorem coeff_one_sigmaSeries (k : ℕ) : PowerSeries.coeff 1 (sigmaSeries k) = 1 := by
  rw [coeff_sigmaSeries]
  simp [ArithmeticFunction.sigma_apply]

theorem coeff_zero_tateA4 : PowerSeries.coeff 0 tateA4 = 0 := by
  rw [tateA4, PowerSeries.coeff_C_mul, coeff_zero_sigmaSeries]
  ring

theorem coeff_one_tateA4 : PowerSeries.coeff 1 tateA4 = -5 := by
  rw [tateA4, PowerSeries.coeff_C_mul, coeff_one_sigmaSeries]
  ring

theorem coeff_zero_tateA6 : PowerSeries.coeff 0 tateA6 = 0 := by
  rw [tateA6, PowerSeries.coeff_mk]
  simp

theorem coeff_one_tateA6 : PowerSeries.coeff 1 tateA6 = -1 := by
  rw [tateA6, PowerSeries.coeff_mk]
  simp [ArithmeticFunction.sigma_apply]

/-! ## ★★判別式 -/

theorem tateCurve_Delta_eq :
    tateCurve.Δ = -tateA6 + tateA4 ^ 2 - 64 * tateA4 ^ 3 - 432 * tateA6 ^ 2
      + 72 * tateA4 * tateA6 := by
  simp only [WeierstrassCurve.Δ, WeierstrassCurve.b₂, WeierstrassCurve.b₄, WeierstrassCurve.b₆,
    WeierstrassCurve.b₈, tateCurve]
  ring

theorem constantCoeff_tateA4 : PowerSeries.constantCoeff tateA4 = 0 := by
  rw [← PowerSeries.coeff_zero_eq_constantCoeff_apply]
  exact coeff_zero_tateA4

theorem constantCoeff_tateA6 : PowerSeries.constantCoeff tateA6 = 0 := by
  rw [← PowerSeries.coeff_zero_eq_constantCoeff_apply]
  exact coeff_zero_tateA6

/-- ★★★`Δ(E_q)` の定数項は 0。 -/
theorem coeff_zero_tateDelta : PowerSeries.coeff 0 tateCurve.Δ = 0 := by
  rw [tateCurve_Delta_eq]
  simp [PowerSeries.coeff_zero_eq_constantCoeff_apply, constantCoeff_tateA4, constantCoeff_tateA6]

/-- ★★★**`Δ(E_q) = q + O(q²)`**——1 次係数は 1。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

★これが局所高さ `v_K(q_E)` が正であることの根拠である。 -/
theorem coeff_one_tateDelta : PowerSeries.coeff 1 tateCurve.Δ = 1 := by
  have hexp : tateCurve.Δ = -tateA6 + tateA4 * tateA4 - 64 * (tateA4 * tateA4 * tateA4)
      - 432 * (tateA6 * tateA6) + 72 * tateA4 * tateA6 := by
    rw [tateCurve_Delta_eq]; ring
  rw [hexp]
  simp [PowerSeries.coeff_one_mul, constantCoeff_tateA4, constantCoeff_tateA6,
    coeff_one_tateA6, coeff_one_tateA4]

/-- ★`Δ` は 0 でない(1 次係数が 1 だから)。 -/
theorem tateDelta_ne_zero : tateCurve.Δ ≠ 0 := by
  intro h
  have h1 := coeff_one_tateDelta
  rw [h] at h1
  simp at h1

/-! ## ★特殊化との両立 -/

variable {R : Type} [CommRing R] {I : Ideal R}

/-- ★特殊化した Tate 曲線の判別式は、判別式の特殊化である。 -/
theorem tateCurveAt_Delta [IsAdicComplete I R] (q : R) (hq : q ∈ I) :
    (tateCurveAt q hq).Δ = evalAdicHom q hq tateCurve.Δ := by
  rw [tateCurveAt, WeierstrassCurve.map_Δ]

/-! ## ★出典の紐付け(`.src`) -/

def coeff_one_tateDelta.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Tate 曲線の判別式の位数が 1 であること)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
