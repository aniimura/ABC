/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.DifferentGlobalBound
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★★★分岐指数から `log-diff` を上から押さえる（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.9。

原文 (GenEll p.9):
> Proposition 1.7. (Conductors and Log Differents) Let

## ★★★★★★★★★★★★★★★★★これは何か —— `hlow` の数値の側

`§9-954` は `Proposition 1.7` の `hup`（**下からの評価**）の
局所から大域への組み立てを取った。★本ファイルは**上からの評価**の側を取る:

    `𝔡_{K/F} ∣ ∏_{P∈S} P^{a_P}`  ⟹
      `log-diff(K) − log-diff(F)  ≤  (∑_P a_P · log N(P)) / [K:ℚ]`

★★これが `Skeleton` の `prop_1_7` が受けている `hlow` の**数値の側**である。

## ★★★機構 —— `§9-954` の鏡像

1. `Ideal.absNorm` は割り切りを保つので `N(𝔡) ≤ N(∏ P^{a_P}) = ∏ N(P)^{a_P}`
2. 対数を取って `log N(𝔡) ≤ ∑ a_P·log N(P)`
3. `logDiffOfField_tower` で読み替える

## ★残っている段（明示）

★★★**割り切り `𝔡 ∣ ∏_{P∈S} P^{a_P}` そのもの**が残る。
これには 2 つが要る:

* `a_P ≔ e_P − 1` で **`P^{e_P} ∤ 𝔡`**（馴分岐）
  ——mathlib の `not_dvd_differentIdeal_of_isCoprime`（`e_P` が剰余標数と素）にある
* **`𝔡` の台が `S` に含まれる**こと（`S` の外では `P ∤ 𝔡`）

★そこから `𝔡 ∣ ∏ P^{a_P}` を出す段は
`dvd_iff_normalizedFactors_le_normalizedFactors` ＋
`count_normalizedFactors_eq`（どちらも mathlib）による**重複度の組み立て**である。
★★★★**野分岐の素点をどう集めるか**が残る——原文が `Σ` の上の食い違いとして
分離している部分にあたる。
-/

namespace ABC3.Found.GenEll

open NumberField

/-! ## ★★★★上からの評価（ノルムの水準） -/

/-- ★★★★**割り切りからノルムの対数の上からの不等式へ**。

★`§9-954` の `sum_log_absNorm_le_of_prod_dvd` の鏡像である
——割り切りの向きが逆になっただけである。 -/
theorem log_absNorm_le_sum_of_dvd_prod (K : Type) [Field K] [NumberField K]
    (S : Finset (Ideal (𝓞 K))) (a : Ideal (𝓞 K) → ℕ) (I : Ideal (𝓞 K))
    (hne : ∀ P ∈ S, P ≠ ⊥) (hI : I ≠ 0)
    (hdvd : I ∣ (∏ P ∈ S, P ^ (a P))) :
    Real.log (Ideal.absNorm I) ≤ ∑ P ∈ S, (a P : ℝ) * Real.log (Ideal.absNorm P) := by
  have hPne : ∀ P ∈ S, ((Ideal.absNorm P : ℝ)) ≠ 0 := by
    intro P hP
    have h : Ideal.absNorm P ≠ 0 := by rw [Ne, Ideal.absNorm_eq_zero_iff]; exact hne P hP
    exact_mod_cast h
  have hIne : Ideal.absNorm I ≠ 0 := by
    rw [Ne, Ideal.absNorm_eq_zero_iff]
    simpa [Ideal.zero_eq_bot] using hI
  have hcast : ((Ideal.absNorm (∏ P ∈ S, P ^ (a P)) : ℕ) : ℝ)
      = ∏ P ∈ S, ((Ideal.absNorm P : ℝ) ^ (a P)) := by
    rw [map_prod]
    push_cast
    exact Finset.prod_congr rfl (fun P _ => by rw [map_pow]; push_cast; ring)
  have hprodne : Ideal.absNorm (∏ P ∈ S, P ^ (a P)) ≠ 0 := by
    intro h0
    have hz : ((Ideal.absNorm (∏ P ∈ S, P ^ (a P)) : ℕ) : ℝ) = 0 := by exact_mod_cast h0
    rw [hcast, Finset.prod_eq_zero_iff] at hz
    obtain ⟨P, hP, hP0⟩ := hz
    exact hPne P hP (pow_eq_zero_iff'.mp hP0).1
  have hexp : Real.log (∏ P ∈ S, ((Ideal.absNorm P : ℝ) ^ (a P)))
      = ∑ P ∈ S, (a P : ℝ) * Real.log (Ideal.absNorm P) := by
    rw [Real.log_prod]
    · exact Finset.sum_congr rfl (fun P _ => Real.log_pow _ _)
    · intro P hP
      exact pow_ne_zero _ (hPne P hP)
  rw [← hexp, ← hcast]
  refine Real.log_le_log ?_ ?_
  · exact_mod_cast Nat.pos_of_ne_zero hIne
  · exact_mod_cast Nat.le_of_dvd (Nat.pos_of_ne_zero hprodne) (map_dvd Ideal.absNorm hdvd)

/-! ## ★★★★★★★★★★★★★★★★★`log-diff` の上からの評価 -/

/-- ★★★★★★★★★★★★★★★★★**分岐指数から `log-diff` を上から押さえる**
—— `Proposition 1.7` の `hlow` の数値の側。

原文 (GenEll p.9):
> Proposition 1.7. (Conductors and Log Differents) Let

    `𝔡_{K/F} ∣ ∏_{P∈S} P^{a_P}`  ⟹
      `log-diff(K) − log-diff(F)  ≤  (∑_P a_P · log N(P)) / [K:ℚ]`

★`§9-954`（下から）と合わせて、`log-diff` の差は**両側から**分岐データで挟まれた。 -/
theorem logDiff_tower_le_sum (F K : Type) [Field F] [NumberField F]
    [Field K] [NumberField K] [Algebra F K]
    (S : Finset (Ideal (𝓞 K))) (hne : ∀ P ∈ S, P ≠ ⊥)
    (a : Ideal (𝓞 K) → ℕ)
    (hI : differentIdeal (𝓞 F) (𝓞 K) ≠ 0)
    (hdvd : differentIdeal (𝓞 F) (𝓞 K) ∣ (∏ P ∈ S, P ^ (a P))) :
    logDiffOfField K - logDiffOfField F
      ≤ (∑ P ∈ S, (a P : ℝ) * Real.log (Ideal.absNorm P)) / (Module.finrank ℚ K : ℝ) := by
  have hfr : (0:ℝ) < (Module.finrank ℚ K : ℝ) := by exact_mod_cast Module.finrank_pos
  have hsum := log_absNorm_le_sum_of_dvd_prod K S a _ hne hI hdvd
  rw [logDiffOfField_tower F K]
  have hdiv : Real.log ((Ideal.absNorm (differentIdeal (𝓞 F) (𝓞 K)) : ℝ)) /
        (Module.finrank ℚ K : ℝ)
      ≤ (∑ P ∈ S, (a P : ℝ) * Real.log (Ideal.absNorm P)) / (Module.finrank ℚ K : ℝ) :=
    (div_le_div_iff_of_pos_right hfr).mpr hsum
  linarith

/-! ## ★出典の紐付け(`.src`) -/

def log_absNorm_le_sum_of_dvd_prod.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 9,
    item := "Proposition 1.7(割り切りからノルムの対数の上からの不等式へ)",
    sectionId := "genell-prop-1-7" }

def logDiff_tower_le_sum.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 9,
    item := "Proposition 1.7(分岐指数から log-diff を上から押さえる)",
    sectionId := "genell-prop-1-7" }

def logDiff_tower_le_sum.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "logDiffOfField_tower(log-diff(K) − log-diff(F) = log N(𝔡)/[K:ℚ])"
      (.inProject "ABC3" "ABC3.Found.GenEll.logDiffOfField_tower") 2,
    .citation "[mathlib]" "not_dvd_differentIdeal_of_isCoprime(馴分岐なら P^{e_P} ∤ 𝔡)"
      (.inMathlib "IsDedekindDomain.not_dvd_differentIdeal_of_isCoprime") 3,
    .citation "[mathlib]"
      "dvd_iff_normalizedFactors_le_normalizedFactors / count_normalizedFactors_eq"
      (.inMathlib "UniqueFactorizationMonoid.count_normalizedFactors_eq") 3,
    .implicitStep
      ("★★★★測定(2026-08-29): §9-954(下から)と本ファイル(上から)で、" ++
       "log-diff の差は**両側から分岐データで挟まれた**。" ++
       "★残るのは割り切り 𝔡 ∣ ∏ P^{a_P} そのものであり、" ++
       "馴分岐の P^{e_P} ∤ 𝔡(mathlib にある)と「𝔡 の台が S に含まれる」" ++
       "から重複度の組み立てで出る") 5,
    .implicitStep
      ("★★野分岐の素点をどう集めるかが残る" ++
       "——原文が Σ の上の食い違いとして分離している部分にあたる") 4 ]

end ABC3.Found.GenEll
