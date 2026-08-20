import ABC3.Found.GaloisRep.CharTwo

/-!
# Galois (G1) 第 7 ブロック —— **`psiComp × ψ₂` の公式**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★mathlib が `complEDS₂` の**明示公式**を持っていた

    complEDS₂ b c d k * b = W(k-1)² W(k+2) − W(k-2) W(k+1)²

★これを `psiComp` に翻訳すると

    psiComp n * ψ₂ = ψ(n-1)² ψ(n+2) − ψ(n-2) ψ(n+1)²

★★**`psiComp` が `ψ` だけで書ける**(`ψ₂` を掛ければ)。
★★★これで標数 2 の計算が **`ψ` の漸化式**に帰着する。

## ★★これが帰納段の入口である

標数 2 で `omegaNum n = psiComp n − ψ_n(a₁φ_n + a₃ψ_n²) = 0` を示したい。
★両辺に `ψ₂` を掛ければ、本ブロックで左辺が `ψ` だけになる
——`ψ₂` は零因子ではない(普遍環は整域)ので**割り戻せる**。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `psiComp_def` | ★`psiComp` の定義(`rfl`) |
| `psiComp_mul_psi2` | ★★★**`psiComp × ψ₂` の明示公式** |
-/

universe u

namespace ABC3.Found.GaloisRep

open WeierstrassCurve Polynomial Polynomial.Bivariate

variable {R : Type u} [CommRing R] (W : WeierstrassCurve R)

/-- ★**`psiComp` の定義**(`rfl`)。 -/
theorem psiComp_def (n : ℤ) :
    psiComp W n = complEDS₂ W.ψ₂ (C W.Ψ₃) (C W.preΨ₄) n := rfl

/-- ★★★**`psiComp × ψ₂` の明示公式**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★★mathlib の `complEDS₂_mul_b` そのものである——
**`psiComp` が `ψ` だけで書ける**ので、標数 2 の計算が `ψ` の漸化式に帰着する。 -/
theorem psiComp_mul_psi2 (n : ℤ) :
    psiComp W n * W.ψ₂ = W.ψ (n - 1) ^ 2 * W.ψ (n + 2) - W.ψ (n - 2) * W.ψ (n + 1) ^ 2 :=
  complEDS₂_mul_b _ _ _ n

/-! ## ★出典の紐付け(`.src`) -/

def psiComp_mul_psi2.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(G1——psiComp × ψ₂ の公式)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
