import ABC3.Found.GaloisRep.OmegaMap
import Mathlib.Algebra.MvPolynomial.CommRing

/-!
# Galois (G1) 第 3 ブロック —— **普遍 Weierstrass 曲線**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★`2 ∣ omegaNum` を示す舞台

mathlib の docstring が示す筋は

    標数 0 の普遍環 ℤ[A₁,A₂,A₃,A₄,A₆][X,Y] で 2 ∣ omegaNum を示し、普遍射で降ろす

★本ブロックはその**普遍環と普遍曲線**を用意する。降ろす段は第 2 ブロック済み。

## ★★測って分かった —— `CharZero (MvPolynomial …)` の instance が**無い**

`Polynomial.charZero`(1 変数)は mathlib に在るが、
★**`MvPolynomial` 版は無い**(2026-08-19 実測)。

★★`constantCoeff` で `ℤ` に落として `Nat.cast_injective` に帰着すれば
**3 行で作れる**ので、自前で置く。

★★★これも「mathlib に無い」の実例だが、
**1 変数版が在るので探せば見つかる形**であった——
`Polynomial` と `MvPolynomial` で整備の粗密があるのは自然である。

## ★★本ブロックで取れるもの

| 定義・定理 | 内容 |
|---|---|
| `URing` | ★普遍環 `ℤ[A₁,…,A₆]` |
| `uCurve` | ★★普遍 Weierstrass 曲線 |
| `charZero_URing` | ★標数 0(自前) |
-/

universe u

namespace ABC3.Found.GaloisRep

open WeierstrassCurve Polynomial Polynomial.Bivariate MvPolynomial

/-- ★**普遍環** `ℤ[A₁,A₂,A₃,A₄,A₆]`。 -/
abbrev URing : Type := MvPolynomial (Fin 5) ℤ

/-- ★★**普遍 Weierstrass 曲線**——係数が不定元。 -/
noncomputable def uCurve : WeierstrassCurve URing where
  a₁ := X 0
  a₂ := X 1
  a₃ := X 2
  a₄ := X 3
  a₆ := X 4

/-- ★**普遍環は標数 0 である**(mathlib に `MvPolynomial` 版が無いので自前)。 -/
instance charZero_URing : CharZero URing := ⟨fun a b h => by
  have := congrArg (MvPolynomial.constantCoeff) h
  simpa using this⟩

example : IsDomain URing := by infer_instance

/-! ## ★出典の紐付け(`.src`) -/

def uCurve.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(G1——普遍 Weierstrass 曲線)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
