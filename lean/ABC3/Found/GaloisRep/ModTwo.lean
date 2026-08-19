import ABC3.Found.GaloisRep.OmegaBase
import Mathlib.Data.ZMod.Basic

/-!
# Galois (G1) 第 5 ブロック —— **`ZMod 2` への還元**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★`2 ∣ omegaNum` を**標数 2 の計算**に落とす

普遍環 `URing = ℤ[A₁,…,A₆]` は**整数係数の多項式環**なので、

    2 ∣ p   ⟺   p の係数を `ZMod 2` に落とすと 0

★★これで「割り切れる」という**存在命題**が、
「標数 2 で消える」という**等式**に変わる——★★★帰納法が回しやすくなる。

## ★★第 2 ブロック(底変換と可換)がここで効く

`map_omegaNum : omegaNum (W.map f) n = (omegaNum W n).map (mapRingHom f)` により、
★`omegaNum uCurve n` を `ZMod 2` に落としたものは
**`omegaNum (uCurve.map modTwo) n`** である——
つまり**標数 2 の Weierstrass 曲線の `omegaNum`** を計算すればよい。

★★★標数 2 では `ψ₂ = a₁X + a₃`(`2Y` が消える)ので、
`omegaNum` の 2 項が相殺する筋が見える。

## ★★本ブロックで取れるもの

| 定義・定理 | 内容 |
|---|---|
| `modTwo` | ★`URing → ℤ/2[A₁,…,A₆]` |
| `modTwo_two` | ★`modTwo 2 = 0` |
-/

universe u

namespace ABC3.Found.GaloisRep

open WeierstrassCurve Polynomial Polynomial.Bivariate

/-- ★**`ZMod 2` への係数写像**。 -/
noncomputable abbrev modTwo : URing →+* MvPolynomial (Fin 5) (ZMod 2) :=
  MvPolynomial.map (Int.castRingHom (ZMod 2))

/-- ★**`modTwo 2 = 0`**。 -/
theorem modTwo_two : modTwo (2 : URing) = 0 := by
  rw [show (2 : URing) = MvPolynomial.C (2 : ℤ) by simp, MvPolynomial.map_C]
  norm_num
  decide

/-! ## ★出典の紐付け(`.src`) -/

def modTwo.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(G1——ZMod 2 への還元)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
