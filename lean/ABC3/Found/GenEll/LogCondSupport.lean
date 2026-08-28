/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.LogCondSigma
import ABC3.Found.GenEll.DifferentDivides
import ABC3.Found.GenEll.DifferentGlobalBound
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★★★★★★★★★`log-cond` は台の上の `∑ log N(P)` である（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.9。

原文 (GenEll p.9):
> Proposition 1.7. (Conductors and Log Differents) Let

## ★★★★★★★★★★★★★★★★★★★★★★★これは何か —— 導手の側の等式

`§9-965` で `Proposition 1.7` の `hup` は

    `(1 − 1/e)·(∑_{P∈S} log N(P)) / [K:ℚ]  ≤  log-diff(K) − log-diff(F)`

の形になった。★残っていたのは左辺の和が原文の **`log-cond`** と一致することである。

★★★★本ファイルはその等式を取る:

    **`log-cond_D(x) = (∑_{P ∈ 台} log N(P)) / [F:ℚ]`**

（`台` ＝ 引き戻したイデアルを割る素イデアルの全体）。

## ★★★機構 —— 根基は台の素イデアルの積である

1. `radical I` は**平方因子をもたない**（`LogCondSigma.lean` の `isRadical_ideal_radical`）
2. だから各素イデアルの重複度は `≤ 1`——`§9-956` の `dvd_prod_pow_of_bounded`（`a_P = 1`）で
   `radical I ∣ ∏_{P∈S} P`
3. 逆向きは `I ≤ P ⟹ radical I ≤ radical P = P` と `§9-954` の
   `prod_pow_dvd_of_pairwise_maximal`
4. ★したがって **`radical I = ∏_{P∈S} P`**——あとは `absNorm` の乗法性と `log` である

## ★残っている段（明示）

★★★`Proposition 1.7` に残るのは**「被覆 `Y → Z` が `E` の上でだけ分岐する」**
——すなわち `§9-965` の `S`（分岐する素点）と本ファイルの `S`（導手の台）が
**同じ集合である**という幾何の仮定だけになった。
-/

namespace ABC3.Found.GenEll

open NumberField AlgebraicGeometry

/-! ## ★★★★★★★★★★★★根基は台の素イデアルの積 -/

/-- ★★★★★★★★★★★★**根基イデアルは台の素イデアルの積である**。

★`radical I` が平方因子をもたないこと（`isRadical_ideal_radical`）と、
`§9-954`／`§9-956` の割り切りの道具だけである。 -/
theorem radical_eq_prod_of_support (F : Type) [Field F] [NumberField F]
    [DecidableEq (Ideal (𝓞 F))]
    (I : Ideal (𝓞 F)) (hI : I ≠ 0)
    (S : Finset (Ideal (𝓞 F))) (hmax : ∀ P ∈ S, P.IsMaximal)
    (hmem : ∀ P ∈ S, P ∣ I)
    (hsupp : ∀ Q : Ideal (𝓞 F), Irreducible Q → Q ∣ I → Q ∈ S) :
    Ideal.radical I = ∏ P ∈ S, P := by
  have hrad0 : I.radical ≠ 0 := radical_ne_zero hI
  have hsq : Squarefree I.radical := (isRadical_ideal_radical I).squarefree hrad0
  have hradI : I.radical ∣ I := Ideal.dvd_iff_le.mpr Ideal.le_radical
  have hPne : ∀ P ∈ S, P ≠ ⊥ := by
    intro P hP h0
    apply hI
    have hd : (0 : Ideal (𝓞 F)) ∣ I := by
      rw [Ideal.zero_eq_bot, ← h0]
      exact hmem P hP
    exact zero_dvd_iff.mp hd
  have hprodne : (∏ P ∈ S, P ^ 1) ≠ 0 := by
    intro h0
    rw [Finset.prod_eq_zero_iff] at h0
    obtain ⟨P, hP, hP0⟩ := h0
    exact hPne P hP (by simpa using hP0)
  have h1 : I.radical ∣ ∏ P ∈ S, P ^ 1 := by
    refine dvd_prod_pow_of_bounded S (fun _ => 1) I.radical hrad0 hprodne ?_ ?_
    · intro Q hQirr hQdvd
      exact hsupp Q hQirr (hQdvd.trans hradI)
    · intro P hP hcon
      have h2 : P * P ∣ I.radical := by
        have hp2 : P ^ (1 + 1) = P * P := by ring
        rwa [hp2] at hcon
      exact (hmax P hP).ne_top (Ideal.isUnit_iff.1 (hsq P h2))
  have h2 : (∏ P ∈ S, P ^ 1) ∣ I.radical := by
    refine prod_pow_dvd_of_pairwise_maximal S hmax (fun _ => 1) I.radical ?_
    intro P hP
    rw [pow_one]
    refine Ideal.dvd_iff_le.mpr ?_
    have hle : I ≤ P := Ideal.dvd_iff_le.mp (hmem P hP)
    calc I.radical ≤ P.radical := Ideal.radical_mono hle
      _ = P := (hmax P hP).isPrime.radical
  have heq : I.radical = ∏ P ∈ S, P ^ 1 :=
    le_antisymm (Ideal.dvd_iff_le.mp h2) (Ideal.dvd_iff_le.mp h1)
  simpa [pow_one] using heq

/-! ## ★★★★★★★★★★★★★★★★★★★★★★★導手の等式 -/

/-- ★★★★★★★★★★★★★★★★★★★★★★★**`log-cond` は台の上の `∑ log N(P)`**。

原文 (GenEll p.9):
> Proposition 1.7. (Conductors and Log Differents) Let

    `log-cond_D(x) = (∑_{P ∈ 台} log N(P)) / [F:ℚ]`

★★これで `Proposition 1.7` の**導手の側の等式**が取れた
——`§9-965`（`hup`）の左辺と同じ形である。
★★★残るのは「その `S` が**分岐する素点の集合と同じ**」という幾何の仮定だけである。 -/
theorem logCond_eq_sum_log (F : Type) [Field F] [NumberField F]
    [DecidableEq (Ideal (𝓞 F))] {X : Scheme.{0}}
    (D : X.IdealSheafData) (xF : specRingOfIntegers F ⟶ X)
    (hI : pullbackIdeal F D xF ≠ 0)
    (S : Finset (Ideal (𝓞 F))) (hmax : ∀ P ∈ S, P.IsMaximal)
    (hne : ∀ P ∈ S, P ≠ ⊥)
    (hmem : ∀ P ∈ S, P ∣ pullbackIdeal F D xF)
    (hsupp : ∀ Q : Ideal (𝓞 F), Irreducible Q → Q ∣ pullbackIdeal F D xF → Q ∈ S) :
    logCond F D xF
      = (∑ P ∈ S, Real.log (Ideal.absNorm P)) / (Module.finrank ℚ F : ℝ) := by
  have hrad0 : (pullbackIdeal F D xF).radical ≠ 0 := radical_ne_zero hI
  rw [logCond, conductorADiv, degNormalized, deg_idealADiv F _ hrad0,
    radical_eq_prod_of_support F _ hI S hmax hmem hsupp]
  congr 1
  rw [map_prod]
  push_cast
  rw [Real.log_prod]
  intro P hP
  have h : Ideal.absNorm P ≠ 0 := by rw [Ne, Ideal.absNorm_eq_zero_iff]; exact hne P hP
  exact_mod_cast h

/-! ## ★出典の紐付け(`.src`) -/

def radical_eq_prod_of_support.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 9,
    item := "Proposition 1.7(根基イデアルは台の素イデアルの積である)",
    sectionId := "genell-prop-1-7" }

def logCond_eq_sum_log.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 9,
    item := "Proposition 1.7(log-cond は台の上の ∑ log N(P) である)",
    sectionId := "genell-prop-1-7" }

def logCond_eq_sum_log.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "isRadical_ideal_radical / radical_ne_zero(根基は平方因子をもたない)"
      (.inProject "ABC3" "ABC3.Found.GenEll.isRadical_ideal_radical") 2,
    .citation "[ABC3]" "dvd_prod_pow_of_bounded(重複度から割り切り、§9-956)"
      (.inProject "ABC3" "ABC3.Found.GenEll.dvd_prod_pow_of_bounded") 3,
    .citation "[ABC3]" "prod_pow_dvd_of_pairwise_maximal(互いに素な積、§9-954)"
      (.inProject "ABC3" "ABC3.Found.GenEll.prod_pow_dvd_of_pairwise_maximal") 2,
    .implicitStep
      ("★★★★測定(2026-08-29): 導手の側の等式は" ++
       "**根基が平方因子をもたない**ことと本日の割り切りの道具" ++
       "(§9-954・§9-956)だけで出る——`radical I = ∏_{P∈台} P` である") 5,
    .implicitStep
      ("★★★これで Proposition 1.7 に残るのは" ++
       "「被覆 Y → Z が E の上でだけ分岐する」" ++
       "——すなわち §9-965 の S(分岐する素点)と本ファイルの S(導手の台)が" ++
       "**同じ集合である**という幾何の仮定だけになった") 5 ]

end ABC3.Found.GenEll
