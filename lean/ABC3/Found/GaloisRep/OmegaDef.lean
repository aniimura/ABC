import ABC3.Found.GaloisRep.OmegaAll
import ABC3.Found.GaloisRep.UniversalCurve
import ABC3.Found.GaloisRep.UniversalF2
import ABC3.Found.GaloisRep.OmegaMap

/-!
# Galois (G1) 第 27 ブロック —— **★★★★★★★★★`ω_n` が定義できる**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★★mathlib の TODO が閉じた

`Mathlib/AlgebraicGeometry/EllipticCurve/DivisionPolynomial/Basic.lean` は

> TODO: define `ω_n`

と書き、その障害を 2 つ挙げていた。★どちらも本区間で消えた:

| 障害 | 解決 |
|---|---|
| `ψₙ ∣ ψ₂ₙ` | ★在庫に在った(`normEDS_mul_complEDS₂`、第 6 ブロックで確認) |
| **`2 ∣ 分子`** | ★★★**第 22–26 ブロック**(標数 2 で `omegaNum = 0`) |

★★普遍環 `ℤ[A₁,…,A₆]` へ落として係数ごとに見れば、
標数 2 での消滅がそのまま **2 での割り切れ**になる。

## ★★本ブロックで取れるもの

| 定義・定理 | 内容 |
|---|---|
| `red` / `uCurve_map_red` | ★普遍環の標数 2 還元 |
| `two_dvd_omegaNum` | ★★★★**`2 ∣ omegaNum uCurve n`** |
| `uOmega` | ★★★★★**`ω_n`(普遍曲線)** |
| `two_mul_uOmega` | ★★`2·ω_n = 分子` |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve Polynomial Polynomial.Bivariate MvPolynomial

/-- 普遍環から標数 2 の普遍環への還元。 -/
noncomputable def red : URing →+* UF2Ring := MvPolynomial.map (Int.castRingHom (ZMod 2))

theorem uCurve_map_red : uCurve.map red = uCurveF2 := by
  simp [WeierstrassCurve.map, uCurve, uCurveF2, red]

theorem map_zero_of_char_two (n : ℤ) :
    (omegaNum uCurve n).map (mapRingHom red) = 0 := by
  rw [← map_omegaNum uCurve red n, uCurve_map_red]
  exact omegaNum_char_two_all uCurveF2 n

theorem two_dvd_omegaNum (n : ℤ) : (2 : URing[X][Y]) ∣ omegaNum uCurve n := by
  have h0 := map_zero_of_char_two n
  rw [show (2 : URing[X][Y]) = Polynomial.C (Polynomial.C (2 : URing)) by
    rw [map_ofNat, map_ofNat]]
  rw [Polynomial.C_dvd_iff_dvd_coeff]
  intro i
  rw [Polynomial.C_dvd_iff_dvd_coeff]
  intro j
  rw [show (2 : URing) = MvPolynomial.C (2 : ℤ) by rw [map_ofNat]]
  rw [MvPolynomial.C_dvd_iff_dvd_coeff]
  intro m
  have h1 : ((omegaNum uCurve n).map (mapRingHom red)).coeff i = 0 := by rw [h0]; simp
  rw [Polynomial.coeff_map] at h1
  have h2 : ((((omegaNum uCurve n).coeff i).map red).coeff j) = 0 := by
    rw [show ((omegaNum uCurve n).coeff i).map red
      = (mapRingHom red) ((omegaNum uCurve n).coeff i) from rfl, h1]
    simp
  rw [Polynomial.coeff_map] at h2
  have h3 : (MvPolynomial.map (Int.castRingHom (ZMod 2))
      (((omegaNum uCurve n).coeff i).coeff j)).coeff m = 0 := by
    rw [show MvPolynomial.map (Int.castRingHom (ZMod 2))
      (((omegaNum uCurve n).coeff i).coeff j) = red (((omegaNum uCurve n).coeff i).coeff j)
        from rfl, h2]
    simp
  rw [MvPolynomial.coeff_map] at h3
  exact (ZMod.intCast_zmod_eq_zero_iff_dvd _ 2).1 h3


/-- ★★★★★**`ω_n`(普遍曲線)**——mathlib が TODO と書いたもの。 -/
noncomputable def uOmega (n : ℤ) : URing[X][Y] := (two_dvd_omegaNum n).choose

/-- ★★**`2·ω_n = 分子`**。 -/
theorem two_mul_uOmega (n : ℤ) : 2 * uOmega n = omegaNum uCurve n :=
  (two_dvd_omegaNum n).choose_spec.symm

/-! ## ★出典の紐付け(`.src`) -/

def two_dvd_omegaNum.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(ω_n の分子が 2 で割り切れること——mathlib の TODO)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
