import ABC3.Found.GaloisRep.OmegaDef

/-!
# Galois (G1) 第 28 ブロック —— **★★★★★★任意の曲線の `ω_n`**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★普遍曲線の `ω_n` を特殊化する

第 27 で `uOmega n : ℤ[A₁,…,A₆][X][Y]` を得た。
★任意の Weierstrass 曲線 `W / R` は普遍曲線の像なので、そのまま移せる。

    ω_n(W) := (uOmega n).map (uHom W)
    2 · ω_n(W) = omegaNum W n

## ★★これで乗法公式の第 3 成分が揃った

原文が使う `[n](x, y) = (φₙ/ψₙ², ωₙ/ψₙ³)` の `ωₙ` である。

## ★★本ブロックで取れるもの

| 定義・定理 | 内容 |
|---|---|
| `uHom` / `uCurve_map_uHom` | ★任意の曲線は普遍曲線の像 |
| `omega` | ★★★★★**任意の曲線の `ω_n`** |
| `two_mul_omega` | ★★★**`2·ω_n = 分子`** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve Polynomial Polynomial.Bivariate MvPolynomial

universe u

variable {R : Type u} [CommRing R] (W : WeierstrassCurve R)

/-- ★**普遍環から `R` への特殊化**。 -/
noncomputable def uHom : URing →+* R :=
  MvPolynomial.eval₂Hom (Int.castRingHom R) ![W.a₁, W.a₂, W.a₃, W.a₄, W.a₆]

/-- ★**任意の Weierstrass 曲線は普遍曲線の像である**。 -/
theorem uCurve_map_uHom : uCurve.map (uHom W) = W := by
  cases W
  simp [WeierstrassCurve.map, uCurve, uHom]

/-- ★★★★★**任意の曲線の `ω_n`**。 -/
noncomputable def omega (n : ℤ) : R[X][Y] :=
  (uOmega n).map (mapRingHom (uHom W))

/-- ★★★**`2·ω_n = 分子`**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★★これが mathlib の `ω_n` の TODO の結論そのものである。 -/
theorem two_mul_omega (n : ℤ) : 2 * omega W n = omegaNum W n := by
  have h : ((2 : URing[X][Y]) * uOmega n).map (mapRingHom (uHom W))
      = (omegaNum uCurve n).map (mapRingHom (uHom W)) := by
    rw [two_mul_uOmega]
  rw [Polynomial.map_mul, Polynomial.map_ofNat] at h
  rw [omega, h, ← map_omegaNum uCurve (uHom W) n, uCurve_map_uHom]

/-! ## ★出典の紐付け(`.src`) -/

def two_mul_omega.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(任意の曲線の ω_n)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
