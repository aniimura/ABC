import ABC3.Found.GaloisRep.PsiCompMul

/-!
# Galois (G1) 第 8 ブロック —— **`omegaNum × ψ₂` を `ψ` だけで書く**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★帰納の舞台が整った

    omegaNum n × ψ₂ = (ψ(n-1)² ψ(n+2) − ψ(n-2) ψ(n+1)²) − ψ_n(a₁φ_n + a₃ψ_n²) ψ₂

★第 7(`psiComp × ψ₂`)を代入しただけである。★★これで右辺は
**`ψ` と `φ` と係数だけ**——`complEDS₂` が消えた。

## ★★標数 2 ではさらに `ψ₂ = a₁X + a₃`

    omegaNum n × ψ₂ = (…) − ψ_n(a₁φ_n + a₃ψ_n²)(a₁X + a₃)      （標数 2)

★★★これで**`ψ` の漸化式**(`normEDS_even` / `normEDS_odd`)で帰納が回せる。

## ★★残るのは `ψ` の漸化式の計算

`ψ` は EDS なので

    ψ(2m) ψ₂ = ψ(m-1)² ψ(m) ψ(m+2) − ψ(m-2) ψ(m) ψ(m+1)²
    ψ(2m+1) = ψ(m+2) ψ(m)³ − ψ(m-1) ψ(m+1)³

★これを標数 2 で使い、`omegaNum × ψ₂ = 0` を示す。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `omegaNum_mul_psi2` | ★★★`omegaNum × ψ₂` を `ψ` だけで |
| `omegaNum_mul_psi2_char_two` | ★★標数 2 版 |
-/

universe u

namespace ABC3.Found.GaloisRep

open WeierstrassCurve Polynomial Polynomial.Bivariate

section General

variable {R : Type u} [CommRing R] (W : WeierstrassCurve R)

/-- ★★★**`omegaNum × ψ₂` を `ψ` だけで書く**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★★第 7 を代入するだけで `complEDS₂` が消える。 -/
theorem omegaNum_mul_psi2 (n : ℤ) :
    omegaNum W n * W.ψ₂
      = (W.ψ (n - 1) ^ 2 * W.ψ (n + 2) - W.ψ (n - 2) * W.ψ (n + 1) ^ 2)
        - W.ψ n * (C (C W.a₁) * W.φ n + C (C W.a₃) * W.ψ n ^ 2) * W.ψ₂ := by
  rw [omegaNum, sub_mul, psiComp_mul_psi2]

end General

section CharTwo

variable {R : Type u} [CommRing R] [CharP R 2] (W : WeierstrassCurve R)

/-- ★★**標数 2 版**——`ψ₂ = a₁X + a₃` を代入した形。 -/
theorem omegaNum_mul_psi2_char_two (n : ℤ) :
    omegaNum W n * W.ψ₂
      = (W.ψ (n - 1) ^ 2 * W.ψ (n + 2) - W.ψ (n - 2) * W.ψ (n + 1) ^ 2)
        - W.ψ n * (C (C W.a₁) * W.φ n + C (C W.a₃) * W.ψ n ^ 2)
          * C (C W.a₁ * X + C W.a₃) := by
  rw [omegaNum_mul_psi2, psi2_char_two]

end CharTwo

/-! ## ★出典の紐付け(`.src`) -/

def omegaNum_mul_psi2.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(G1——omegaNum × ψ₂ を ψ だけで書く)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
