import Mathlib.AlgebraicGeometry.EllipticCurve.DivisionPolynomial.Basic
import ABC3.Meta.Claim

/-!
# Galois (G1) 第 1 ブロック —— **`ω_n` の分子**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★Galois 側の**最初のブロック**

§9-224 / §9-234 / §9-235 で測った結果、G1(`E[n] ≅ (ℤ/n)²`)の律速は
`#E[n] = n²` であり、そのまた律速が**乗法公式の第 3 成分 `ω_n`** である。

    ωₙ := (ψ₂ₙ / ψₙ - ψₙ ⬝ (a₁φₙ + a₃ψₙ²)) / 2

mathlib は `DivisionPolynomial/Basic.lean` でこれを **TODO** と書き、
FLT も `n_torsion_card` を sorry のまま「David Angdinata が作業中」と注記している。

## ★★★★★測って分かった —— **障害の半分は既に済んでいた**

docstring は 2 つの障害を挙げる:

| 障害 | ★実測(2026-08-19) |
|---|---|
| (1) `ψₙ ∣ ψ₂ₙ`(「帰納法で示せる」) | ★★**`NumberTheory/EllipticDivisibilitySequence.lean` に在る** |
| (2) `2 ∣ (分子)`(「普遍環で示して降ろす」) | ★無い |

★★(1) は `complEDS₂`(2-補完列)と `normEDS_mul_complEDS₂` として在る——
**`ψ₂ₙ / ψₙ` は `complEDS₂` そのもの**である。
★★★docstring の TODO を額面で受け取らず、**別ファイルまで grep した**ので分かった。

## ★★本ブロックで取れるもの

| 定義・定理 | 内容 |
|---|---|
| `psiComp` | ★★`ψ₂ₙ / ψₙ`(= `complEDS₂`) |
| `psi_mul_psiComp` | ★★★**`ψₙ · psiComp n = ψ₂ₙ`**(割り算が定義できる根拠) |
| `omegaNum` | ★★★★**`ω_n` の分子** |
| `omegaNum_zero` / `omegaNum_one` | ★初期値 |
-/

universe u

namespace ABC3.Found.GaloisRep

open WeierstrassCurve Polynomial Polynomial.Bivariate

variable {R : Type u} [CommRing R] (W : WeierstrassCurve R)

/-- ★★**`ψ₂ₙ / ψₙ`** —— mathlib の 2-補完列そのもの。 -/
noncomputable def psiComp (n : ℤ) : R[X][Y] :=
  complEDS₂ W.ψ₂ (C W.Ψ₃) (C W.preΨ₄) n

/-- ★★★**`ψₙ · psiComp n = ψ₂ₙ`**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★★これが「`ψ₂ₙ / ψₙ` が多項式として定義できる」ことの根拠であり、
mathlib の docstring が「帰納法で示せる」と書いた部分である
——★**別ファイルに既に在った**。 -/
theorem psi_mul_psiComp (n : ℤ) : W.ψ n * psiComp W n = W.ψ (2 * n) :=
  normEDS_mul_complEDS₂ _ _ _ n

/-- ★★★★**`ω_n` の分子** —— これを 2 で割ると `ω_n` になる。 -/
noncomputable def omegaNum (n : ℤ) : R[X][Y] :=
  psiComp W n - W.ψ n * (C (C W.a₁) * W.φ n + C (C W.a₃) * W.ψ n ^ 2)

theorem omegaNum_zero : omegaNum W 0 = 2 := by
  simp [omegaNum, psiComp]

theorem omegaNum_one : omegaNum W 1 = W.ψ₂ - (C (C W.a₁) * W.φ 1 + C (C W.a₃)) := by
  simp [omegaNum, psiComp]

/-! ## ★出典の紐付け(`.src`) -/

def omegaNum.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(G1——分多項式の乗法公式の第 3 成分 ω_n の分子)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
