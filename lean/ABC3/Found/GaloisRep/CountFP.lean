import ABC3.Found.GaloisRep.WeilMul

/-!
# Galois (G5) 第 185 ブロック —— **★★★★★`f_P` の指数と素点の輸送(一般形)**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★交代性に向けた足場

`e_n(P,P) = 1` の古典的な証明は、

    ∏_{i=0}^{n-1} τ_{iP}^*(f_P)  は定数

を使う。★これを**指数の計算**で言うために、次の 2 つを用意する。

### ★★★★★`f_P` の指数

`(f_P) = XYIdeal(P)^n` なので、mathlib の `count_pow` と `count_maximal` で

    count_v(f_P) = if v = P の素点 then n else 0

が直ちに出る(`count_fP`)。

### ★★★★★★素点の輸送(一般形)

第 170 の `count_translate_eq` は `τ z = z` を仮定していた(`z = μ f_P` に使うため)。
★交代性では `z = f_P` で、これは `τ` で動く。★★そこで**仮定を外した形**を切り出す:

    count_v(τ z) = count_{v'}(z)        (Q_{v'} = Q_v + T)

★★★中身は第 169 の `count_comp_eq` に第 170 の仮説をそのまま渡すだけである。

## ★★これで交代性の指数計算ができる

    count_v(∏_{i} τ_{iP}^* f_P) = Σ_i [ n·[Q_v + iP = P] − n·[Q_v + iP = O] ]

となり、`Q_v ∈ ⟨P⟩` なら各和が 1 項ずつで打ち消し、そうでなければ両方 0。
★どちらでも 0 なので `∏` は単元、すなわち**定数**(第 128)。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `count_fP` | ★★★★★**`f_P` の指数** |
| `count_translate_gen` | ★★★★★★**素点の輸送(一般形)** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine Polynomial IsDedekindDomain nonZeroDivisors

variable {F : Type} [Field F] [DecidableEq F] (W : WeierstrassCurve.Affine F)
  [inst : IsDedekindDomain W.CoordinateRing]

/-- ★★★★★**`f_P` の指数**——`P` の素点で `n`、他では 0。 -/
theorem count_fP [DecidableEq (HeightOneSpectrum W.CoordinateRing)]
    (v : HeightOneSpectrum W.CoordinateRing) {x y : F} (hP : W.Nonsingular x y)
    (n : ℕ) (f : W.CoordinateRing)
    (hf : (CoordinateRing.XYIdeal W x (Polynomial.C y)) ^ n = Ideal.span {f}) :
    FractionalIdeal.count W.FunctionField v
      (FractionalIdeal.spanSingleton W.CoordinateRing⁰
        (algebraMap W.CoordinateRing W.FunctionField f))
      = if pointSpectrum W hP.1 = v then (n : ℤ) else 0 := by
  rw [← FractionalIdeal.coeIdeal_span_singleton, ← hf, FractionalIdeal.coeIdeal_pow,
    FractionalIdeal.count_pow]
  have hid : (CoordinateRing.XYIdeal W x (Polynomial.C y))
      = (pointSpectrum W hP.1).asIdeal := rfl
  rw [hid, FractionalIdeal.count_maximal]
  by_cases hc : pointSpectrum W hP.1 = v
  · rw [if_pos hc, if_pos hc, mul_one]
  · rw [if_neg hc, if_neg hc, mul_zero]

variable [W.IsElliptic]

/-- ★★★★★★**素点の輸送(一般形)**——`τ z = z` を要求しない。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★第 170 の `count_translate_eq` から不変性の仮定を外したもの。 -/
theorem count_translate_gen (h2 : IsUnit (2 : F))
    (v v' : HeightOneSpectrum W.CoordinateRing)
    {c y₀ c' y₀' : F} (h : W.Equation c y₀)
    (hv : v.asIdeal = CoordinateRing.XYIdeal W c (Polynomial.C y₀))
    (h' : W.Equation c' y₀')
    (hv' : v'.asIdeal = CoordinateRing.XYIdeal W c' (Polynomial.C y₀'))
    (τ : W.FunctionField ≃ₐ[F] W.FunctionField)
    {x₀ y₀T : F} (hQ : W.Nonsingular x₀ y₀T)
    (hx : τ (coordX W) = translateX W x₀ y₀T) (hy : τ (coordY W) = translateY W x₀ y₀T)
    (hsum : Point.some c y₀ (equation_iff_nonsingular.mp h) + Point.some x₀ y₀T hQ
      = Point.some c' y₀' (equation_iff_nonsingular.mp h'))
    {z : W.FunctionField} (hz : z ≠ 0) :
    FractionalIdeal.count W.FunctionField v
        (FractionalIdeal.spanSingleton W.CoordinateRing⁰ (τ z))
      = FractionalIdeal.count W.FunctionField v'
        (FractionalIdeal.spanSingleton W.CoordinateRing⁰ z) := by
  obtain ⟨hle, hlt⟩ := aut_transport_hyps W v h hv h2 τ hQ hx hy
    (equation_iff_nonsingular.mp h') hsum
  have hQ' : W.Nonsingular x₀ (W.negY x₀ y₀T) := (nonsingular_neg x₀ y₀T).mpr hQ
  obtain ⟨hx', hy'⟩ := autFF_symm_coord W τ hQ hx hy
  have hsum' : Point.some c' y₀' (equation_iff_nonsingular.mp h')
      + Point.some x₀ (W.negY x₀ y₀T) hQ'
      = Point.some c y₀ (equation_iff_nonsingular.mp h) := by
    rw [← hsum, ← Point.neg_some hQ, add_assoc, add_neg_cancel, add_zero]
  obtain ⟨hle', hlt'⟩ := aut_transport_hyps W v' h' hv' h2 τ.symm hQ' hx' hy'
    (equation_iff_nonsingular.mp h) hsum'
  rw [← hv'] at hlt
  rw [← hv] at hlt'
  exact count_comp_eq τ.toRingEquiv v v' hle hlt hle' hlt' hz

/-! ## ★出典の紐付け(`.src`) -/

def count_fP.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——f_P の指数)",
    sectionId := "genell-thm-3-8" }

def count_translate_gen.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——素点の輸送、一般形)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
