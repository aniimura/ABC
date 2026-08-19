import ABC3.Found.GaloisRep.UniversalCurve

/-!
# Galois (G1) 第 4 ブロック —— **`ω_n` の初期値**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★`2 ∣ omegaNum` の帰納法の**土台**

`omegaNum n := psiComp n - ψ_n ⬝ (a₁φ_n + a₃ψ_n²)` について
`2 ∣ omegaNum n` を示すには、初期値から始める帰納法が要る。

| `n` | `omegaNum n` | `2 ∣` |
|---|---|---|
| 0 | ★`2` | ★**明らか** |
| 1 | `ψ₂ - (a₁φ₁ + a₃)` | ★`ψ₂ = 2Y + a₁X + a₃`、`φ₁ = X` から |

## ★★具体形は `rfl` と `simp` で出る

| 定理 | 証明 |
|---|---|
| `psi2_val`(`ψ₂ = 2Y + (a₁X + a₃)`) | ★**`rfl`**(`polynomialY` の定義そのもの) |
| `phi_one`(`φ₁ = X`) | ★`simp`(`ψ₁ = 1`、`ψ₀ = 0`) |

★★★`ψ₂` が **`2Y + …`** の形をしているのが、
`ω_1` の分子が偶数になる理由である——**`2Y` の項が `a₁X + a₃` を相殺する**。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `two_dvd_zero` | ★`2 ∣ omegaNum 0` |
| `psi2_val` | ★★`ψ₂` の具体形 |
| `phi_one` | ★`φ₁ = X` |
-/

universe u

namespace ABC3.Found.GaloisRep

open WeierstrassCurve Polynomial Polynomial.Bivariate

variable {R : Type u} [CommRing R] (W : WeierstrassCurve R)

/-- ★**`2 ∣ omegaNum 0`**——分子は `2` そのもの。 -/
theorem two_dvd_zero : (2 : R[X][Y]) ∣ omegaNum W 0 := by
  rw [omegaNum_zero]

/-- ★★**`ψ₂ = 2Y + (a₁X + a₃)`**(`rfl`)。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★★★`2Y` の項があるのが `ω_n` の分子が偶数になる根拠である。 -/
theorem psi2_val : W.ψ₂ = C (C 2) * Y + C (C W.a₁ * X + C W.a₃) := rfl

/-- ★**`φ₁ = X`**。 -/
theorem phi_one : W.φ 1 = Polynomial.C Polynomial.X := by
  rw [WeierstrassCurve.φ]
  simp

/-! ## ★出典の紐付け(`.src`) -/

def psi2_val.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(G1——ω_n の初期値)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
