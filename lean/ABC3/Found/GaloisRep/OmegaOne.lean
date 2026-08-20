import ABC3.Found.GaloisRep.PsiRec

/-!
# Galois (G1) 第 10 ブロック —— ★★★★★**標数 2 で `omegaNum 1 = 0`**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★帰納の**基底が 2 つ揃った**

| `n` | 標数 2 での `omegaNum n` | ブロック |
|---|---|---|
| 0 | ★`0`(分子は `2`) | 第 6 |
| ★1 | ★**`0`** | ★本ブロック |

## ★★`n = 1` が消える理由

    omegaNum 1 = ψ₂ − (a₁φ₁ + a₃)
               = (a₁X + a₃) − (a₁·X + a₃)      （標数 2、φ₁ = X)
               = 0

★★★**`ψ₂` の `2Y` が消えると、残りが `a₁X + a₃` でちょうど第 2 項と一致する**
——これが `ω_n` が整数係数で定義できる根拠の**最小の場合**である。

## ★★`simp` 一発で閉じた

第 4(`omegaNum_one`)+ 第 6(`psi2_char_two`)+ 第 4(`phi_one`)を
`rw` してから `simp` するだけであった。★`ring` は要らなかった。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `psiComp_one` | ★`psiComp 1 = ψ₂` |
| `psi_one` | ★`ψ 1 = 1` |
| `omegaNum_one_char_two` | ★★★★★**標数 2 で `omegaNum 1 = 0`** |
-/

universe u

namespace ABC3.Found.GaloisRep

open WeierstrassCurve Polynomial Polynomial.Bivariate

section General

variable {R : Type u} [CommRing R] (W : WeierstrassCurve R)

/-- ★**`psiComp 1 = ψ₂`**。 -/
theorem psiComp_one : psiComp W 1 = W.ψ₂ := by
  rw [psiComp]
  exact complEDS₂_one _ _ _

/-- ★**`ψ 1 = 1`**。 -/
theorem psi_one : W.ψ 1 = 1 := by
  rw [WeierstrassCurve.ψ]
  exact normEDS_one _ _ _

end General

section CharTwo

variable {R : Type u} [CommRing R] [CharP R 2] (W : WeierstrassCurve R)

/-- ★★★★★**標数 2 では `omegaNum 1 = 0`**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★★★`ψ₂` の `2Y` が消えると残りが `a₁X + a₃` で、第 2 項とちょうど一致する
——これが `ω_n` が整数係数で定義できる根拠の**最小の場合**である。 -/
theorem omegaNum_one_char_two : omegaNum W 1 = 0 := by
  rw [omegaNum_one, psi2_char_two, phi_one]
  simp

end CharTwo

/-! ## ★出典の紐付け(`.src`) -/

def omegaNum_one_char_two.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(G1——標数 2 で omegaNum 1 = 0)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
