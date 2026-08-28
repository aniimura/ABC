/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.HomogRatio
import ABC3.Found.GenEll.ProjSpaceCover
import ABC3.Meta.Claim

/-!
# ★★★★`A⁰_{x_i}` の分母は `x_i^k`、次数もちょうど `k` —— 段 E1c の包装（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.7。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

## ★★★★★★これは何か

環準同型 `A⁰_{x_i} → Γ(X, X_{s_i})` を作る（段 E1c）には、`Away` の元

    `c : NumDenSameDeg 𝒜 (powers (x_i))`   （`deg`・`num ∈ 𝒜 deg`・`den ∈ 𝒜 deg`・`den ∈ powers x_i`）

から **`den = x_i^k` と `deg = k`** を取り出す必要がある。★本ファイルはそれである。

## ★★★機構 —— 斉次な非零元の次数は一意

`den_mem` が `den = x_i^k` を与える。★`x_i^k` は次数 `k` の斉次式であり、
**`0` ではない**（単項式の係数が `1`）。
★★したがって `MvPolynomial.IsHomogeneous.inj_right`（斉次な非零元の次数は一意）で
`deg = k` が出る。★★★分子 `num` も `𝒜 deg = 𝒜 k` の元なので次数 `k` の斉次式である。

## ★測定の記録

`x_i^k ≠ 0` に `pow_ne_zero` は使えない——`R` 一般では `MvPolynomial` に
零因子がありうるからである（`IsReduced` のインスタンスが要求された）。
★`X_pow_eq_monomial` で単項式に直すと `simp` が落とす（係数が `1`）。

## ★残っている段（明示）

★★これで段 E1c に要る材料は揃った:

| 材料 | 場所 |
|---|---|
| 同値関係を潰す（`x_i^m·a = x_i^n·b` ⟹ 比が等しい） | `homogRatio_congr_of_mul_eq`（`§9-812`） |
| 分母と次数の取り出し | ★本ファイル |
| `Quotient.lift` で環準同型に仕上げる | ★★残る |

★★★残るのは **`Quotient.lift` の包装**（加法・乗法・`0`・`1` の確認）だけである。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory MvPolynomial

/-- ★★★★**`Away (x_i)` の分母は `x_i^k` であり、次数もちょうど `k` である**。

原文 (GenEll p.7):
> that [some positive tensor power of] the ample line bundle LQ yields an embedding

★機構は「斉次な非零元の次数は一意」（`IsHomogeneous.inj_right`）である。 -/
theorem exists_pow_of_numDenSameDeg {N : ℕ} {R : Type} [CommRing R] [Nontrivial R]
    (i : Fin (N + 1))
    (c : HomogeneousLocalization.NumDenSameDeg
      (MvPolynomial.homogeneousSubmodule (Fin (N + 1)) R)
      (Submonoid.powers (MvPolynomial.X i : MvPolynomial (Fin (N + 1)) R))) :
    ∃ k : ℕ, (c.den : MvPolynomial (Fin (N + 1)) R) = MvPolynomial.X i ^ k ∧ c.deg = k := by
  obtain ⟨k, hk⟩ := c.den_mem
  refine ⟨k, hk.symm, ?_⟩
  have hpow : ((MvPolynomial.X i : MvPolynomial (Fin (N + 1)) R) ^ k).IsHomogeneous k := by
    simpa using (MvPolynomial.isHomogeneous_X R i).pow k
  have hne : (MvPolynomial.X i : MvPolynomial (Fin (N + 1)) R) ^ k ≠ 0 := by
    rw [MvPolynomial.X_pow_eq_monomial]; simp
  have hden : ((c.den : MvPolynomial (Fin (N + 1)) R)).IsHomogeneous c.deg :=
    (MvPolynomial.mem_homogeneousSubmodule _ _).1 c.den.2
  rw [← hk] at hden
  exact MvPolynomial.IsHomogeneous.inj_right hden hpow hne

/-- ★★★**分子も同じ次数 `k` の斉次式である**。 -/
theorem num_isHomogeneous_of_numDenSameDeg {N : ℕ} {R : Type} [CommRing R] [Nontrivial R]
    (i : Fin (N + 1))
    (c : HomogeneousLocalization.NumDenSameDeg
      (MvPolynomial.homogeneousSubmodule (Fin (N + 1)) R)
      (Submonoid.powers (MvPolynomial.X i : MvPolynomial (Fin (N + 1)) R)))
    {k : ℕ} (hk : c.deg = k) :
    ((c.num : MvPolynomial (Fin (N + 1)) R)).IsHomogeneous k :=
  hk ▸ (MvPolynomial.mem_homogeneousSubmodule _ _).1 c.num.2

/-! ## ★出典の紐付け(`.src`) -/

def exists_pow_of_numDenSameDeg.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(A⁰_{x_i} の分母は x_i^k で次数もちょうど k)",
    sectionId := "genell-prop-1-4" }

def num_isHomogeneous_of_numDenSameDeg.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 7,
    item := "Proposition 1.4, (iv)(分子も同じ次数 k の斉次式である)",
    sectionId := "genell-prop-1-4" }

def exists_pow_of_numDenSameDeg.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]" "MvPolynomial.IsHomogeneous.inj_right(斉次な非零元の次数は一意)"
      (.inMathlib "MvPolynomial.IsHomogeneous.inj_right") 7,
    .citation "[ABC3]" "homogRatio_congr_of_mul_eq(同値関係を潰す、§9-812)"
      (.inProject "ABC3" "ABC3.Found.GenEll.homogRatio_congr_of_mul_eq") 7,
    .implicitStep
      ("★x_i^k ≠ 0 に pow_ne_zero は使えない——R 一般では MvPolynomial に" ++
       "零因子がありうる(IsReduced のインスタンスが要求された、2026-08-28 実測)。" ++
       "★★X_pow_eq_monomial で単項式に直すと simp が落とす") 7,
    .implicitStep
      ("★★★これで段 E1c に要る材料は揃った。残るのは Quotient.lift の包装" ++
       "(加法・乗法・0・1 の確認)だけである") 7 ]

end ABC3.Found.GenEll
