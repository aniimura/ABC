import Mathlib.Algebra.CharP.Algebra
import Mathlib.Algebra.Field.ZMod
import ABC3.Found.GaloisRep.UniversalCurve
import ABC3.Found.GaloisRep.CharTwo

/-!
# Galois (G1) 第 19 ブロック —— **標数 2 の普遍曲線と `ψ₂ ≠ 0`**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★§9-349 の罠を避ける場所

標数 2 の一般の環では `ψ₂ = a₁X + a₃` が **0 になりうる**(`a₁ = a₃ = 0` の曲線)。
★したがって帰納の途中で `ψ₂` を約せない。

★★**普遍環 `𝔽₂[A₁,…,A₆]` の上では `a₁ = A₁`、`a₃ = A₃` が不定元なので `ψ₂ ≠ 0`**であり、
`R[X][Y]` は整域だから**約せる**。

## ★★降ろす道は在庫にある

第 9 ブロック `map_omegaNum` が `omegaNum (W.map f) n = (omegaNum W n).map …` を与えるので、
★普遍曲線で `omegaNum = 0` を示せば、**任意の標数 2 の曲線へ降りる**。

## ★★本ブロックで取れるもの

| 定義・定理 | 内容 |
|---|---|
| `UF2Ring` | ★標数 2 の普遍環 `𝔽₂[A₁,…,A₆]` |
| `uCurveF2` | ★普遍曲線 |
| `uCurveF2_psi2_ne_zero` | ★★★★**`ψ₂ ≠ 0`**——約分が使える |
| `uCurveF2_map` | ★★★任意の標数 2 の曲線は普遍曲線の像である |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve Polynomial Polynomial.Bivariate MvPolynomial

universe u

/-- ★**標数 2 の普遍環** `𝔽₂[A₁,A₂,A₃,A₄,A₆]`。 -/
abbrev UF2Ring : Type := MvPolynomial (Fin 5) (ZMod 2)

/-- ★**標数 2 の普遍 Weierstrass 曲線**。 -/
noncomputable def uCurveF2 : WeierstrassCurve UF2Ring where
  a₁ := X 0
  a₂ := X 1
  a₃ := X 2
  a₄ := X 3
  a₆ := X 4

instance : CharP UF2Ring 2 :=
  charP_of_injective_ringHom (MvPolynomial.C_injective (Fin 5) (ZMod 2)) 2

instance : Fact (Nat.Prime 2) := ⟨Nat.prime_two⟩

/-- ★**普遍環は整域である**——ここでだけ約分が使える。 -/
instance : IsDomain UF2Ring := inferInstance

instance : IsDomain (Polynomial UF2Ring) := inferInstance

instance : IsDomain (Polynomial (Polynomial UF2Ring)) := inferInstance

/-- ★★★★**普遍曲線では `ψ₂ ≠ 0`**——ここでだけ約分が使える。 -/
theorem uCurveF2_psi2_ne_zero : uCurveF2.ψ₂ ≠ 0 := by
  rw [psi2_char_two]
  intro h
  rw [← map_zero Polynomial.C] at h
  have h1 : (Polynomial.C (uCurveF2.a₁) * Polynomial.X
      + Polynomial.C (uCurveF2.a₃) : UF2Ring[X]) = 0 :=
    Polynomial.C_injective h
  have h2 := congrArg (fun p : UF2Ring[X] => p.coeff 1) h1
  simp [uCurveF2] at h2

/-- ★★★**任意の標数 2 の曲線は普遍曲線の像である**。 -/
noncomputable def uCurveF2Hom {R : Type u} [CommRing R] [CharP R 2] (W : WeierstrassCurve R) :
    UF2Ring →+* R :=
  MvPolynomial.eval₂Hom (ZMod.castHom (dvd_refl 2) R)
    ![W.a₁, W.a₂, W.a₃, W.a₄, W.a₆]

theorem uCurveF2_map {R : Type u} [CommRing R] [CharP R 2] (W : WeierstrassCurve R) :
    uCurveF2.map (uCurveF2Hom W) = W := by
  cases W
  simp [WeierstrassCurve.map, uCurveF2, uCurveF2Hom]

/-! ## ★出典の紐付け(`.src`) -/

def uCurveF2_psi2_ne_zero.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(標数 2 の普遍曲線——ψ₂ が非零であること)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
