/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.LogDiffTower
import Mathlib.RingTheory.DedekindDomain.Different
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★★分岐指数から `log-diff` を下から押さえる（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.9。

原文 (GenEll p.9):
> Proposition 1.7. (Conductors and Log Differents) Let

## ★★★★★★★★★★★★★★★★これは何か —— `Proposition 1.7` の**局所から大域へ**

`Skeleton/GenEll/Section1.lean` の `prop_1_7` は `hlow`／`hup` を**仮定として受けて**おり、
`.needs` に「★局所の核は馴分岐 6/6 で実装した。★★残っているのは
**局所から大域への組み立て**である」と記録していた。

★★★★本ファイルはその組み立ての**下からの評価**（`hup` の側）を取る:

    各素点 `P` で `P^{a_P} ∣ 𝔡_{K/F}` なら
      `(∑_P a_P · log N(P)) / [K:ℚ]  ≤  log-diff(K) − log-diff(F)`

★`a_P ≔ e_P − 1` と取れるのは mathlib の `pow_sub_one_dvd_differentIdeal`
（`P^e ∣ p·𝓞_K` ⟹ `P^{e−1} ∣ 𝔡`）である——**局所の核は mathlib にあった**。

## ★★★機構 —— 3 段

1. **相異なる極大イデアルの冪は互いに素**なので、個々に割り切れば**積も割り切る**
   （`Finset.prod_dvd_of_coprime`）
2. `Ideal.absNorm` は乗法的かつ**割り切りを保つ**ので、ノルムの不等式になる
3. `logDiffOfField_tower`（`log-diff(K) − log-diff(F) = log N(𝔡)/[K:ℚ]`）で読み替える

## ★残っている段（明示）

★★`hlow`（上からの評価）は**馴分岐の仮定が要る**——`P^{e_P} ∤ 𝔡` は
mathlib の `not_dvd_differentIdeal_of_isCoprime`（`e_P` が剰余標数と素）で出るが、
**野分岐の素点をどう集めるか**が残る。
★これは原文が `Σ` の上の食い違いとして分離している部分にあたる。
-/

namespace ABC3.Found.GenEll

open NumberField

/-! ## ★★★段 1 —— 互いに素な冪の積は割り切る -/

/-- ★★★**相異なる極大イデアルの冪は、個々に割り切れば積も割り切る**。

★Dedekind 環では相異なる極大イデアルが互いに素（`P ⊔ Q = ⊤`）だから、
`Finset.prod_dvd_of_coprime` がそのまま当たる。 -/
theorem prod_pow_dvd_of_pairwise_maximal {B : Type} [CommRing B] [IsDedekindDomain B]
    (S : Finset (Ideal B)) (hmax : ∀ P ∈ S, P.IsMaximal)
    (a : Ideal B → ℕ) (I : Ideal B) (h : ∀ P ∈ S, P ^ (a P) ∣ I) :
    (∏ P ∈ S, P ^ (a P)) ∣ I := by
  refine Finset.prod_dvd_of_coprime ?_ h
  intro P hP Q hQ hPQ
  exact (Ideal.isCoprime_iff_sup_eq.mpr
    (Ideal.IsMaximal.coprime_of_ne (hmax P hP) (hmax Q hQ) hPQ)).pow

/-! ## ★★★★段 2 —— ノルムを取って対数にする -/

/-- ★★★★**割り切りからノルムの対数の不等式へ**。

★`Ideal.absNorm` は乗法的で割り切りを保つので、
`∏ P^{a_P} ∣ I` から `∑ a_P·log N(P) ≤ log N(I)` が出る。 -/
theorem sum_log_absNorm_le_of_prod_dvd (K : Type) [Field K] [NumberField K]
    (S : Finset (Ideal (𝓞 K))) (a : Ideal (𝓞 K) → ℕ) (I : Ideal (𝓞 K)) (hI : I ≠ 0)
    (hne : ∀ P ∈ S, P ≠ ⊥)
    (hdvd : (∏ P ∈ S, P ^ (a P)) ∣ I) :
    ∑ P ∈ S, (a P : ℝ) * Real.log (Ideal.absNorm P) ≤ Real.log (Ideal.absNorm I) := by
  have hPne : ∀ P ∈ S, ((Ideal.absNorm P : ℝ)) ≠ 0 := by
    intro P hP
    have h : Ideal.absNorm P ≠ 0 := by rw [Ne, Ideal.absNorm_eq_zero_iff]; exact hne P hP
    exact_mod_cast h
  have hIne : Ideal.absNorm I ≠ 0 := by
    rw [Ne, Ideal.absNorm_eq_zero_iff]
    simpa [Ideal.zero_eq_bot] using hI
  have hnormdvd : Ideal.absNorm (∏ P ∈ S, P ^ (a P)) ∣ Ideal.absNorm I :=
    map_dvd Ideal.absNorm hdvd
  have hcast : ((Ideal.absNorm (∏ P ∈ S, P ^ (a P)) : ℕ) : ℝ)
      = ∏ P ∈ S, ((Ideal.absNorm P : ℝ) ^ (a P)) := by
    rw [map_prod]
    push_cast
    exact Finset.prod_congr rfl (fun P _ => by rw [map_pow]; push_cast; ring)
  have hprodpos : (0:ℝ) < ∏ P ∈ S, ((Ideal.absNorm P : ℝ) ^ (a P)) := by
    refine Finset.prod_pos (fun P hP => ?_)
    have h1 : (0:ℝ) ≤ (Ideal.absNorm P : ℝ) := by positivity
    exact pow_pos (lt_of_le_of_ne h1 (Ne.symm (hPne P hP))) _
  have hexp : Real.log (∏ P ∈ S, ((Ideal.absNorm P : ℝ) ^ (a P)))
      = ∑ P ∈ S, (a P : ℝ) * Real.log (Ideal.absNorm P) := by
    rw [Real.log_prod]
    · exact Finset.sum_congr rfl (fun P _ => Real.log_pow _ _)
    · intro P hP
      exact pow_ne_zero _ (hPne P hP)
  rw [← hexp, ← hcast]
  refine Real.log_le_log ?_ ?_
  · rw [hcast]; exact hprodpos
  · exact_mod_cast Nat.le_of_dvd (Nat.pos_of_ne_zero hIne) hnormdvd

/-! ## ★★★★★★★★★★★★★★★★段 3 —— `log-diff` の下からの評価 -/

/-- ★★★★★★★★★★★★★★★★**分岐指数から `log-diff` を下から押さえる**
—— `Proposition 1.7` の `hup` の局所から大域への組み立て。

原文 (GenEll p.9):
> Proposition 1.7. (Conductors and Log Differents) Let

    各素点 `P` で `P^{a_P} ∣ 𝔡_{K/F}` なら
      `(∑_P a_P · log N(P)) / [K:ℚ]  ≤  log-diff(K) − log-diff(F)`

★`a_P ≔ e_P − 1` と取れるのは mathlib の `pow_sub_one_dvd_differentIdeal` である。
★★これが原文の『the elementary theory of differents』の**下からの側**にあたる。 -/
theorem sum_log_le_logDiff_tower (F K : Type) [Field F] [NumberField F]
    [Field K] [NumberField K] [Algebra F K]
    (S : Finset (Ideal (𝓞 K))) (hmax : ∀ P ∈ S, P.IsMaximal) (hne : ∀ P ∈ S, P ≠ ⊥)
    (a : Ideal (𝓞 K) → ℕ)
    (hdvd : ∀ P ∈ S, P ^ (a P) ∣ differentIdeal (𝓞 F) (𝓞 K))
    (hI : differentIdeal (𝓞 F) (𝓞 K) ≠ 0) :
    (∑ P ∈ S, (a P : ℝ) * Real.log (Ideal.absNorm P)) / (Module.finrank ℚ K : ℝ)
      ≤ logDiffOfField K - logDiffOfField F := by
  have hfr : (0:ℝ) < (Module.finrank ℚ K : ℝ) := by exact_mod_cast Module.finrank_pos
  have hsum := sum_log_absNorm_le_of_prod_dvd K S a _ hI hne
    (prod_pow_dvd_of_pairwise_maximal S hmax a _ hdvd)
  rw [logDiffOfField_tower F K]
  have hdiv : (∑ P ∈ S, (a P : ℝ) * Real.log (Ideal.absNorm P)) / (Module.finrank ℚ K : ℝ)
      ≤ Real.log ((Ideal.absNorm (differentIdeal (𝓞 F) (𝓞 K)) : ℝ)) /
        (Module.finrank ℚ K : ℝ) :=
    (div_le_div_iff_of_pos_right hfr).mpr hsum
  linarith

/-! ## ★出典の紐付け(`.src`) -/

def prod_pow_dvd_of_pairwise_maximal.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 9,
    item := "Proposition 1.7(相異なる極大イデアルの冪は個々に割り切れば積も割り切る)",
    sectionId := "genell-prop-1-7" }

def sum_log_absNorm_le_of_prod_dvd.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 9,
    item := "Proposition 1.7(割り切りからノルムの対数の不等式へ)",
    sectionId := "genell-prop-1-7" }

def sum_log_le_logDiff_tower.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 9,
    item := "Proposition 1.7(分岐指数から log-diff を下から押さえる——局所から大域へ)",
    sectionId := "genell-prop-1-7" }

def sum_log_le_logDiff_tower.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]" "pow_sub_one_dvd_differentIdeal(P^e ∣ p·𝓞_K ⟹ P^{e−1} ∣ 𝔡)"
      (.inMathlib "IsDedekindDomain.pow_sub_one_dvd_differentIdeal") 3,
    .citation "[ABC3]" "logDiffOfField_tower(log-diff(K) − log-diff(F) = log N(𝔡)/[K:ℚ])"
      (.inProject "ABC3" "ABC3.Found.GenEll.logDiffOfField_tower") 2,
    .citation "[mathlib]" "Finset.prod_dvd_of_coprime(互いに素なら積も割り切る)"
      (.inMathlib "Finset.prod_dvd_of_coprime") 1,
    .implicitStep
      ("★★★★測定(2026-08-29): Skeleton の prop_1_7 が仮定として受けている hup の" ++
       "**局所から大域への組み立て**のうち、下からの評価はこれで取れた。" ++
       "★局所の核(P^{e−1} ∣ 𝔡)は**mathlib にあった**" ++
       "——pow_sub_one_dvd_differentIdeal である") 5,
    .implicitStep
      ("★★残るのは hlow(上からの評価)であり、**馴分岐の仮定が要る**" ++
       "——P^{e_P} ∤ 𝔡 は mathlib の not_dvd_differentIdeal_of_isCoprime" ++
       "(e_P が剰余標数と素)で出るが、野分岐の素点をどう集めるかが残る。" ++
       "★これは原文が Σ の上の食い違いとして分離している部分にあたる") 4 ]

end ABC3.Found.GenEll
