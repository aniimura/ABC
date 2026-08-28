/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.DifferentKummer
import ABC3.Found.GenEll.TameRamification
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★★★★★★★★馴分岐では different はちょうど `𝔪^{e−1}`（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.9。

原文 (GenEll p.9):
> Proposition 1.7. (Conductors and Log Differents) Let

## ★★★★★★★★★★★★★★★★★★★★★★これは何か —— 残っていた最後の 1 本

`§9-958` で測ったとおり、`Proposition 1.7` の elementary claim に残っていたのは

    **`P^{e} ∤ 𝔡`（馴分岐）**——上からの側

だけであった。★mathlib には**不分岐（`e = 1`）の場合しか無く**
（`not_dvd_differentIdeal_of_isCoprime`）、本プロジェクトの馴分岐 6/6 が持つのも
`λ^m ∈ 𝔡`（`m ≥ e−1`）——**下からの側**であった。

★★★★**本ファイルはそれを埋める**。しかも**等式**が出る:

    馴分岐（`e` が `B` で単元）＋ Eisenstein ＋ 単項生成  ⟹  **`𝔡 = (λ^{e−1})`**

★★したがって `(λ)^e ∤ 𝔡` である——`λ` が単元でないからである。

## ★★★機構 —— すでにある 2 本を掛け合わせるだけだった

| 部品 | 出どころ |
|---|---|
| `𝔡 = ( f'(λ) )`（単項生成なら different は微分の生成） | `§9-397`（`DifferentKummer.lean`） |
| `f'(λ) = λ^{e−1}·(単元)`（Eisenstein ＋ 馴） | `§9-393`（`TameRamification.lean`） |

★**単元倍は `span` で消える**（`Ideal.span_singleton_mul_right_unit`）ので、
`𝔡 = (λ^{e−1})` が**そのまま**出る。
★★★2026-08-29 の実測で「上からの側は誰も持っていない」と書いたが、
**部品は両方とも手元にあり、掛け合わせるだけだった**。

## ★これで `Proposition 1.7` の elementary claim は

`hup`（下から）は `§9-958` で分岐指数だけになり、
`hlow`（上から）に要る `P^e ∤ 𝔡` は本ファイルで出た。
★★残るのは**局所（`IsLocalRing B`）から大域（数体の整数環）への持ち上げ**である
——各素点で局所化して本定理を当てる段。
-/

namespace ABC3.Found.GenEll

open Polynomial Finset IsLocalRing

/-! ## ★★★★★★★★★★★★★★★★★★★★★★馴分岐の different は `(λ^{e−1})` である -/

/-- ★★★★★★★★★★★★★★★★★★★★★★**馴分岐では `𝔡 = (λ^{e−1})` ちょうど**。

原文 (GenEll p.9):
> Proposition 1.7. (Conductors and Log Differents) Let

★機構は 2 本の掛け合わせ:
* `𝔡 = ( f'(λ) )`（単項生成、`§9-397`）
* `f'(λ) = λ^{e−1}·(単元)`（Eisenstein ＋ 馴、`§9-393`）

★★単元倍は `span` で消えるので等式がそのまま出る。 -/
theorem differentIdeal_eq_span_pow_of_eisenstein_tame
    (A : Type*) (K : Type*) (L : Type*) {B : Type*} [CommRing A] [Field K] [CommRing B] [Field L]
    [Algebra A K] [Algebra B L] [Algebra A B] [Algebra K L] [Algebra A L]
    [IsScalarTower A K L] [IsScalarTower A B L] [IsDomain A] [IsFractionRing A K]
    [FiniteDimensional K L] [Algebra.IsSeparable K L] [IsIntegralClosure B A L]
    [IsIntegrallyClosed A] [IsDedekindDomain B] [Module.IsTorsionFree A B]
    [IsLocalRing B]
    (e : ℕ) (he : 0 < e) (a : ℕ → B) (lam : B)
    (hdvd : ∀ i ∈ range e, lam ^ e ∣ a i)
    (hlam : lam ∈ maximalIdeal B) (hunit : IsUnit (e : B))
    (hmono : Algebra.adjoin A {lam} = ⊤)
    (hgen : Algebra.adjoin K {(algebraMap B L) lam} = ⊤)
    (hmin : (minpoly A lam).map (algebraMap A B)
      = X ^ e + ∑ i ∈ range e, C (a i) * X ^ i) :
    differentIdeal A B = Ideal.span {lam ^ (e - 1)} := by
  obtain ⟨v, hv, hveq⟩ := aeval_derivative_eisenstein_tame e he a lam hdvd hlam hunit
  have hkey : (aeval lam) (derivative (minpoly A lam)) = lam ^ (e - 1) * v := by
    rw [← hveq, ← hmin, derivative_map, Polynomial.aeval_map_algebraMap]
  rw [differentIdeal_eq_span_of_adjoin_eq_top A K L lam hmono hgen, hkey]
  exact Ideal.span_singleton_mul_right_unit hv _

/-! ## ★★★★★★★★★★★上からの非割り切り -/

/-- ★★★★★★★★★★★**`(λ)^e` は `(λ^{e−1})` を割らない**——`λ` が単元でないから。 -/
theorem not_span_pow_dvd_span_pow_sub_one {B : Type*} [CommRing B]
    [IsDedekindDomain B] (lam : B) (e : ℕ) (he : 0 < e)
    (hlam0 : lam ≠ 0) (hlamnu : ¬ IsUnit lam) :
    ¬ ((Ideal.span {lam} : Ideal B) ^ e ∣ Ideal.span {lam ^ (e - 1)}) := by
  rw [Ideal.span_singleton_pow]
  intro hdvd
  have hle : (Ideal.span {lam ^ (e - 1)} : Ideal B) ≤ Ideal.span {lam ^ e} :=
    Ideal.le_of_dvd hdvd
  rw [Ideal.span_singleton_le_span_singleton] at hle
  obtain ⟨c, hc⟩ := hle
  apply hlamnu
  have hsplit : lam ^ e = lam ^ (e - 1) * lam := by
    conv_lhs => rw [show e = (e - 1) + 1 by omega]
    rw [pow_succ]
  rw [hsplit, mul_assoc] at hc
  have hpow : lam ^ (e - 1) ≠ 0 := pow_ne_zero _ hlam0
  have h1 : lam * c = 1 := (mul_left_cancel₀ hpow (by rw [mul_one]; exact hc)).symm
  exact IsUnit.of_mul_eq_one c h1

/-! ## ★★★★★★★★★★★★★★★★★★★★★★★到達 —— `hlow` に要る非割り切り -/

/-- ★★★★★★★★★★★★★★★★★★★★★★★**馴分岐なら `(λ)^e ∤ 𝔡`**。

原文 (GenEll p.9):
> Proposition 1.7. (Conductors and Log Differents) Let

★★★これが `§9-956`／`§9-957` が受けていた `hbound` の**局所版**であり、
`Proposition 1.7` の elementary claim に残っていた**最後の 1 本**である。
★★mathlib は不分岐（`e = 1`）しか持っていなかった。 -/
theorem not_pow_dvd_differentIdeal_of_eisenstein_tame
    (A : Type*) (K : Type*) (L : Type*) {B : Type*} [CommRing A] [Field K] [CommRing B] [Field L]
    [Algebra A K] [Algebra B L] [Algebra A B] [Algebra K L] [Algebra A L]
    [IsScalarTower A K L] [IsScalarTower A B L] [IsDomain A] [IsFractionRing A K]
    [FiniteDimensional K L] [Algebra.IsSeparable K L] [IsIntegralClosure B A L]
    [IsIntegrallyClosed A] [IsDedekindDomain B] [Module.IsTorsionFree A B]
    [IsLocalRing B]
    (e : ℕ) (he : 0 < e) (a : ℕ → B) (lam : B)
    (hdvd : ∀ i ∈ range e, lam ^ e ∣ a i)
    (hlam : lam ∈ maximalIdeal B) (hunit : IsUnit (e : B))
    (hmono : Algebra.adjoin A {lam} = ⊤)
    (hgen : Algebra.adjoin K {(algebraMap B L) lam} = ⊤)
    (hmin : (minpoly A lam).map (algebraMap A B)
      = X ^ e + ∑ i ∈ range e, C (a i) * X ^ i)
    (hlam0 : lam ≠ 0) :
    ¬ ((Ideal.span {lam} : Ideal B) ^ e ∣ differentIdeal A B) := by
  rw [differentIdeal_eq_span_pow_of_eisenstein_tame A K L e he a lam hdvd hlam hunit
    hmono hgen hmin]
  exact not_span_pow_dvd_span_pow_sub_one lam e he hlam0 (mem_nonunits_iff.mp hlam)

/-! ## ★出典の紐付け(`.src`) -/

def differentIdeal_eq_span_pow_of_eisenstein_tame.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 9,
    item := "Proposition 1.7(馴分岐では different はちょうど (λ^{e−1}))",
    sectionId := "genell-prop-1-7" }

def not_span_pow_dvd_span_pow_sub_one.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 9,
    item := "Proposition 1.7((λ)^e は (λ^{e−1}) を割らない)",
    sectionId := "genell-prop-1-7" }

def not_pow_dvd_differentIdeal_of_eisenstein_tame.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 9,
    item := "Proposition 1.7(馴分岐なら (λ)^e ∤ 𝔡——hlow の局所版)",
    sectionId := "genell-prop-1-7" }

def not_pow_dvd_differentIdeal_of_eisenstein_tame.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "differentIdeal_eq_span_of_adjoin_eq_top(単項生成なら 𝔡 = (f'(λ))、§9-397)"
      (.inProject "ABC3" "ABC3.Found.GenEll.differentIdeal_eq_span_of_adjoin_eq_top") 3,
    .citation "[ABC3]" "aeval_derivative_eisenstein_tame(f'(λ) = λ^{e−1}·単元、§9-393)"
      (.inProject "ABC3" "ABC3.Found.GenEll.aeval_derivative_eisenstein_tame") 3,
    .implicitStep
      ("★★★★★測定(2026-08-29): Proposition 1.7 の elementary claim に残っていた" ++
       "最後の 1 本(馴分岐での P^e ∤ 𝔡)は、**部品が両方とも手元にあった**" ++
       "——𝔡 = (f'(λ))(§9-397)と f'(λ) = λ^{e−1}·単元(§9-393)。" ++
       "★単元倍は span で消えるので **𝔡 = (λ^{e−1}) が等式で出る**") 6,
    .implicitStep
      ("★★残るのは局所(IsLocalRing B)から大域(数体の整数環)への持ち上げである" ++
       "——各素点で局所化して本定理を当てる段。" ++
       "★hup の側は §9-958 で分岐指数だけになっている") 4 ]

end ABC3.Found.GenEll
