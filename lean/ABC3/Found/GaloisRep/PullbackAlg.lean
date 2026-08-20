import ABC3.Found.GaloisRep.DivisorTools

/-!
# Galois (G5) 第 141 ブロック —— **★★★★★`f_P · f_{−P}` は `n` 乗**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★mathlib が持っていた関係式

mathlib の `CoordinateRing.XYIdeal_neg_mul` は

    XYIdeal(−P) · XYIdeal(P) = XIdeal(x_P) = (x − x_P)

を与える(2026-08-20 実測)。★これを `n` 乗すると `(f_P) · (f_{−P}) = ((x − x_P)^n)` となり、
第 128 ブロック(単元は定数)から

    f_P · f_{−P} = c · (x − x_P)^n        (c ∈ F^×)

が出る。★★**幾何的には `div(f_P) + div(f_{−P}) = n(P) + n(−P) − 2n(O)` そのもの**である。

## ★★★これが D1 で効く

`ord_v(μ f_P) + ord_v(μ f_{−P}) = n · ord_v(μ x − x_P)` となり、
**右辺は `n` で割れる**。★因子の `n` 可除性を調べる際の拘束条件になる。

## ★★★★`μ` は `F`-代数射である

第 119 ブロック(`pointHom_algebraMap`)により、`exists_mulByNHom` が返す `μ` は
`F` 上の代数射である。★スケルトンの仮定にこれを含められるよう、
`exists_mulByNHom` を**その情報も返す形に強めた**。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `exists_mulByNHom_alg` | ★★★`μ` が `F`-代数射であることも返す版 |
| `exists_mulByNHom_alg_charZero` | ★★★その標数 0 版 |
| `fP_mul_fNegP` | ★★★★★**`f_P · f_{−P} = c · (x − x_P)^n`** |
| `mu_fP_mul_mu_fNegP` | ★★★★★**その引き戻し版** |
-/

namespace ABC3.Found.GaloisRep

open WeierstrassCurve WeierstrassCurve.Affine Polynomial

variable {F : Type} [Field F]

/-! ## ★★★`μ` が `F`-代数射であることも返す版 -/

/-- ★★★**`[n]` の引き戻し `μ`——`F`-代数射であることも返す**。

★第 119 ブロック `pointHom_algebraMap` を添えただけ。 -/
theorem exists_mulByNHom_alg (W : WeierstrassCurve.Affine F) [W.IsElliptic] (n : ℕ)
    (hn : n • genericPoint W ≠ 0) :
    ∃ (x y : W.FunctionField) (h : (W.map (algebraMap F W.FunctionField)).Nonsingular x y),
      n • genericPoint W = Point.some x y h ∧
      ∃ μ : W.CoordinateRing →+* W.FunctionField,
        μ (genX W) = x ∧ μ (genY W) = y ∧
        ∀ c : F, μ (algebraMap F W.CoordinateRing c) = algebraMap F W.FunctionField c := by
  obtain ⟨x, y, h, hP⟩ := exists_coords_of_ne_zero (n • genericPoint W) hn
  exact ⟨x, y, h, hP, pointHom W h.1, pointHom_genX W h.1, pointHom_genY W h.1,
    pointHom_algebraMap W h.1⟩

/-- ★★★その標数 0 版——生成点が捩れ点でないこと(第 125)を使う。 -/
theorem exists_mulByNHom_alg_charZero [DecidableEq F] [CharZero F]
    (W : WeierstrassCurve.Affine F) [W.IsElliptic] (n : ℕ) (hn : 1 ≤ n) :
    ∃ (x y : W.FunctionField) (h : (W.map (algebraMap F W.FunctionField)).Nonsingular x y),
      n • genericPoint W = Point.some x y h ∧
      ∃ μ : W.CoordinateRing →+* W.FunctionField,
        μ (genX W) = x ∧ μ (genY W) = y ∧
        ∀ c : F, μ (algebraMap F W.CoordinateRing c) = algebraMap F W.FunctionField c :=
  exists_mulByNHom_alg W n (genericPoint_not_torsion_charZero W n hn)

/-! ## ★★★★★`f_P · f_{−P}` は `n` 乗 -/

/-- ★★★★★**`f_P · f_{−P} = c · (x − x_P)^n`**(`c ∈ F^×`)。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★mathlib の `XYIdeal_neg_mul` を `n` 乗し、第 128(単元は定数)で単元を落とす。
★★幾何的には `div(f_P) + div(f_{−P}) = n(P) + n(−P) − 2n(O)` である。 -/
theorem fP_mul_fNegP (W : WeierstrassCurve.Affine F) [W.IsElliptic] {x y : F}
    (h : W.Nonsingular x y) (n : ℕ) (fP fN : W.CoordinateRing)
    (hfP : (CoordinateRing.XYIdeal W x (Polynomial.C y)) ^ n = Ideal.span {fP})
    (hfN : (CoordinateRing.XYIdeal W x (Polynomial.C (W.negY x y))) ^ n = Ideal.span {fN}) :
    ∃ c : F, c ≠ 0 ∧
      fP * fN = algebraMap F W.CoordinateRing c
        * (genX W - algebraMap F W.CoordinateRing x) ^ n := by
  have hmul : Ideal.span {fP * fN}
      = Ideal.span {(genX W - algebraMap F W.CoordinateRing x) ^ n} := by
    rw [← Ideal.span_singleton_mul_span_singleton, ← hfP, ← hfN, ← mul_pow,
      mul_comm (CoordinateRing.XYIdeal W x (Polynomial.C y)),
      CoordinateRing.XYIdeal_neg_mul h, CoordinateRing.XIdeal, XClass_eq,
      ← Ideal.span_singleton_pow]
  obtain ⟨u, hu⟩ := Ideal.span_singleton_eq_span_singleton.1 hmul.symm
  obtain ⟨c, hc0, hc⟩ := isUnit_coordinateRing (u := (u : W.CoordinateRing)) u.isUnit
  exact ⟨c, hc0, by rw [← hu, hc, mul_comm]⟩

/-- ★★★★★**その引き戻し版**——`μ f_P · μ f_{−P} = c · (μ x − x_P)^n`。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★これにより `ord_v(μ f_P) + ord_v(μ f_{−P})` は `n` で割れる。 -/
theorem mu_fP_mul_mu_fNegP (W : WeierstrassCurve.Affine F) [W.IsElliptic] {x y : F}
    (h : W.Nonsingular x y) (n : ℕ) (fP fN : W.CoordinateRing)
    (hfP : (CoordinateRing.XYIdeal W x (Polynomial.C y)) ^ n = Ideal.span {fP})
    (hfN : (CoordinateRing.XYIdeal W x (Polynomial.C (W.negY x y))) ^ n = Ideal.span {fN})
    (μ : W.CoordinateRing →+* W.FunctionField)
    (hμF : ∀ c : F, μ (algebraMap F W.CoordinateRing c) = algebraMap F W.FunctionField c) :
    ∃ c : F, c ≠ 0 ∧
      μ fP * μ fN = algebraMap F W.FunctionField c
        * (μ (genX W) - algebraMap F W.FunctionField x) ^ n := by
  obtain ⟨c, hc0, hc⟩ := fP_mul_fNegP W h n fP fN hfP hfN
  refine ⟨c, hc0, ?_⟩
  have hq := congrArg μ hc
  rwa [map_mul, map_mul, map_pow, map_sub, hμF, hμF] at hq

/-! ## ★出典の紐付け(`.src`) -/

def fP_mul_fNegP.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(Weil 対の構成——f_P·f_{−P} が (x − x_P)^n の定数倍であること)",
    sectionId := "genell-thm-3-8" }

end ABC3.Found.GaloisRep
