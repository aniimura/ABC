import ABC3.Found.GaloisRep.OmegaMulPsi

/-!
# Galois (G1) 第 9 ブロック —— **`ψ` の漸化式を揃える**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★帰納に使う 4 本を 1 箇所に

| 定理 | 内容 | 出どころ |
|---|---|---|
| `psi_even` | `ψ(2m) ψ₂ = ψ(m-1)²ψ(m)ψ(m+2) − ψ(m-2)ψ(m)ψ(m+1)²` | ★mathlib `normEDS_even` |
| `psi_odd` | `ψ(2m+1) = ψ(m+2)ψ(m)³ − ψ(m-1)ψ(m+1)³` | ★mathlib `normEDS_odd` |
| `first_term` | `omegaNum` の第 1 項は `psiComp × ψ₂` | ★第 7 |
| `psi_two_mul` | `ψ_n · psiComp n = ψ(2n)` | ★第 1 |

## ★★★これで `psi_even` と `psi_two_mul` が**同じことを言っている**と見える

    psi_even:      ψ(2m) ψ₂ = ψ(m-1)²ψ(m)ψ(m+2) − ψ(m-2)ψ(m)ψ(m+1)²
    psi_two_mul:   ψ(2m)    = ψ(m) · psiComp m
    first_term:    ψ(m-1)²ψ(m+2) − ψ(m-2)ψ(m+1)² = psiComp m · ψ₂

★★3 本を並べると **`ψ(m) × first_term = psi_even`** である——整合している。
★★★つまり `psiComp` の定義は `ψ` の漸化式と**同じ情報**であり、
帰納は「`ψ` の漸化式を標数 2 で読む」ことに尽きる。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `psi_even` / `psi_odd` | ★`ψ` の漸化式(mathlib の再掲) |
| `first_term` | ★★第 1 項の同定 |
| `psi_two_mul` | ★`ψ(2n)` の分解 |
-/

universe u

namespace ABC3.Found.GaloisRep

open WeierstrassCurve Polynomial Polynomial.Bivariate

variable {R : Type u} [CommRing R] (W : WeierstrassCurve R)

/-- ★**`ψ` の偶数漸化式**(mathlib `normEDS_even`)。 -/
theorem psi_even (m : ℤ) :
    W.ψ (2 * m) * W.ψ₂
      = W.ψ (m - 1) ^ 2 * W.ψ m * W.ψ (m + 2) - W.ψ (m - 2) * W.ψ m * W.ψ (m + 1) ^ 2 :=
  normEDS_even _ _ _ m

/-- ★**`ψ` の奇数漸化式**(mathlib `normEDS_odd`)。 -/
theorem psi_odd (m : ℤ) :
    W.ψ (2 * m + 1) = W.ψ (m + 2) * W.ψ m ^ 3 - W.ψ (m - 1) * W.ψ (m + 1) ^ 3 :=
  normEDS_odd _ _ _ m

/-- ★★**`omegaNum × ψ₂` の第 1 項は `psiComp × ψ₂`**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★★★これと `psi_even` を並べると、`psiComp` の定義が `ψ` の漸化式と
**同じ情報**であることが見える。 -/
theorem first_term (n : ℤ) :
    W.ψ (n - 1) ^ 2 * W.ψ (n + 2) - W.ψ (n - 2) * W.ψ (n + 1) ^ 2
      = psiComp W n * W.ψ₂ :=
  (psiComp_mul_psi2 W n).symm

/-- ★**`ψ(2n) = ψ_n · psiComp n`**(第 1 ブロックの再掲)。 -/
theorem psi_two_mul (n : ℤ) : W.ψ n * psiComp W n = W.ψ (2 * n) :=
  psi_mul_psiComp W n

/-! ## ★出典の紐付け(`.src`) -/

def first_term.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(G1——ψ の漸化式を揃える)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
