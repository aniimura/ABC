import ABC3.Found.GaloisRep.OmegaFour
import ABC3.Found.GaloisRep.OmegaMap

/-!
# Galois (G1) 第 21 ブロック —— **普遍曲線から降ろす**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★「普遍環で示して降ろす」の**降ろす**段

mathlib の docstring が `ω_n` の TODO に書いた道の後半である。

★第 19 の `uCurveF2_map`(任意の標数 2 の曲線は普遍曲線の像)と
第 9 の `map_omegaNum`(`ω_n` の分子は底変換と可換)を繋ぐだけ。

## ★★これで残りは帰納段だけになった

    帰納段(普遍曲線の上、`normEDSRec`) ⟹ 本ブロック ⟹ 任意の標数 2 の曲線で `omegaNum = 0`
      ⟹ 普遍環 ℤ[a₁,…,a₆] で `2 ∣ omegaNum` ⟹ `ω_n` が多項式として定義できる

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `omegaNum_char_two_of_univ` | ★★★★**普遍曲線から任意の標数 2 の曲線へ降ろす** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve Polynomial Polynomial.Bivariate

universe u

/-- ★★★★**普遍曲線で消えるなら、任意の標数 2 の曲線で消える**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★★★これが「普遍環で示して降ろす」の**降ろす**段である。 -/
theorem omegaNum_char_two_of_univ (h : ∀ n : ℤ, omegaNum uCurveF2 n = 0)
    {R : Type u} [CommRing R] [CharP R 2] (W : WeierstrassCurve R) (n : ℤ) :
    omegaNum W n = 0 := by
  have hmap : omegaNum (uCurveF2.map (uCurveF2Hom W)) n
      = (omegaNum uCurveF2 n).map (mapRingHom (uCurveF2Hom W)) :=
    map_omegaNum uCurveF2 (uCurveF2Hom W) n
  rw [uCurveF2_map W] at hmap
  rw [hmap, h n, Polynomial.map_zero]

/-! ## ★出典の紐付け(`.src`) -/

def omegaNum_char_two_of_univ.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(普遍曲線から任意の標数 2 の曲線へ降ろす段)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
