/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.DifferentGlobalBoundUpper
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★★★★★重複度から割り切りを出す（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.9。

原文 (GenEll p.9):
> Proposition 1.7. (Conductors and Log Differents) Let

## ★★★★★★★★★★★★★★★★★★★これは何か —— `hlow` の最後の一段

`§9-955` は「`𝔡 ∣ ∏_{P∈S} P^{a_P}` なら `log-diff` の差が上から押さえられる」を取った。
★その**割り切りそのもの**が残っていた。

★★★★本ファイルはそれを**重複度から**出す:

    `𝔡` の既約因子はすべて `S` に入り（`hsupp`）、
    各 `P ∈ S` で `P^{a_P+1} ∤ 𝔡`（`hbound`）
      ⟹ `𝔡 ∣ ∏_{P∈S} P^{a_P}`

★★これで `Proposition 1.7` の elementary claim は
**両側とも「局所の分岐データ ⟹ 大域の不等式」の形で閉じた**。

## ★★★機構 —— 数えるだけ

* `P^n ∣ I ↔ n ≤ count P (normalizedFactors I)`（`pow_dvd_iff_le_emultiplicity` ＋
  `emultiplicity_eq_count_normalizedFactors`、どちらも mathlib）
* `I ∣ J ↔ normalizedFactors I ≤ normalizedFactors J`（mathlib）
* ★あとは `Multiset.le_iff_count` で `Q` ごとに場合分けするだけ:
  * `Q` が既約でない → `count Q (nf I) = 0`
  * `Q ∉ S` → `Q ∤ I` だから `count = 0`
  * `Q ∈ S` → `hbound` から `count ≤ a_Q`、`Q^{a_Q} ∣ ∏` から `a_Q ≤ count Q (nf ∏)`

★★**積の因子分解を計算する必要はない**——`Finset.dvd_prod_of_mem` で足りる。

## ★残っている段（明示）

★★`hsupp`（`𝔡` の台が `S` に入る）と `hbound`（`P^{e_P} ∤ 𝔡`）を
**分岐の言葉から供給する**段が残る:

* `hbound` は馴分岐なら mathlib の `not_dvd_differentIdeal_of_isCoprime` にある
* `hsupp` は「不分岐なら `P ∤ 𝔡`」——★これも `not_dvd_differentIdeal_iff` の系である

★★★したがって残るのは**野分岐の素点の扱い**だけであり、
それが原文の `Σ` の上の食い違いにあたる。
-/

namespace ABC3.Found.GenEll

open NumberField UniqueFactorizationMonoid Multiset

/-! ## ★★★冪の割り切りは重複度で測れる -/

/-- ★★★**`P^n ∣ I ↔ n ≤ (`I` の正規化因子に `P` が現れる回数`)`**。

★mathlib の `pow_dvd_iff_le_emultiplicity` と
`emultiplicity_eq_count_normalizedFactors` を繋いだだけである
（イデアルでは `normalize = id` なので `normalize P` が消える）。 -/
theorem pow_dvd_iff_le_count {B : Type} [CommRing B] [IsDedekindDomain B]
    [DecidableEq (Ideal B)] {P I : Ideal B} (hP : Irreducible P) (hI : I ≠ 0) (n : ℕ) :
    P ^ n ∣ I ↔ n ≤ Multiset.count P (normalizedFactors I) := by
  rw [pow_dvd_iff_le_emultiplicity, le_emultiplicity_iff_replicate_le_normalizedFactors hP hI,
    normalize_eq, Multiset.le_count_iff_replicate_le]

/-! ## ★★★★★★★★★★★★★★★★★★★重複度の上界から割り切りへ -/

/-- ★★★★★★★★★★★★★★★★★★★**台と重複度の上界から割り切りが出る**。

原文 (GenEll p.9):
> Proposition 1.7. (Conductors and Log Differents) Let

    `I` の既約因子がすべて `S` に入り、各 `P ∈ S` で `P^{a_P+1} ∤ I`
      ⟹ `I ∣ ∏_{P∈S} P^{a_P}`

★★これが `§9-955` が受けていた割り切りの供給元である。
★★★**積の因子分解を計算する必要はない**——`Q ∈ S` の側は
`Finset.dvd_prod_of_mem` で `Q^{a_Q} ∣ ∏` を出せば足りる。 -/
theorem dvd_prod_pow_of_bounded {B : Type} [CommRing B] [IsDedekindDomain B]
    [DecidableEq (Ideal B)]
    (S : Finset (Ideal B)) (a : Ideal B → ℕ) (I : Ideal B) (hI : I ≠ 0)
    (hprodne : (∏ P ∈ S, P ^ (a P)) ≠ 0)
    (hsupp : ∀ Q : Ideal B, Irreducible Q → Q ∣ I → Q ∈ S)
    (hbound : ∀ P ∈ S, ¬ (P ^ (a P + 1) ∣ I)) :
    I ∣ ∏ P ∈ S, P ^ (a P) := by
  rw [dvd_iff_normalizedFactors_le_normalizedFactors hI hprodne, Multiset.le_iff_count]
  intro Q
  by_cases hQirr : Irreducible Q
  · by_cases hQS : Q ∈ S
    · have hle : Multiset.count Q (normalizedFactors I) ≤ a Q := by
        rcases Nat.lt_or_ge (a Q) (Multiset.count Q (normalizedFactors I)) with h | h
        · exact absurd ((pow_dvd_iff_le_count hQirr hI _).mpr h) (hbound Q hQS)
        · exact h
      have hge : a Q ≤ Multiset.count Q (normalizedFactors (∏ P ∈ S, P ^ (a P))) :=
        (pow_dvd_iff_le_count hQirr hprodne _).mp (Finset.dvd_prod_of_mem _ hQS)
      exact le_trans hle hge
    · have h0 : Multiset.count Q (normalizedFactors I) = 0 := by
        rw [Multiset.count_eq_zero]
        intro hmem
        exact hQS (hsupp Q hQirr (dvd_of_mem_normalizedFactors hmem))
      rw [h0]
      exact Nat.zero_le _
  · have h0 : Multiset.count Q (normalizedFactors I) = 0 := by
      rw [Multiset.count_eq_zero]
      intro hmem
      exact hQirr (irreducible_of_normalized_factor Q hmem)
    rw [h0]
    exact Nat.zero_le _

/-! ## ★★★★★★★★★★★★★★★★★★★★`hlow` の到達形 -/

/-- ★★★★★★★★★★★★★★★★★★★★**分岐データだけから `log-diff` の上からの評価**
—— `Proposition 1.7` の `hlow`。

原文 (GenEll p.9):
> Proposition 1.7. (Conductors and Log Differents) Let

    `𝔡` の既約因子がすべて `S` に入り、各 `P ∈ S` で `P^{a_P+1} ∤ 𝔡`
      ⟹ `log-diff(K) − log-diff(F) ≤ (∑_P a_P·log N(P)) / [K:ℚ]`

★★`§9-954`（下から）と合わせて、`log-diff` の差は
**局所の分岐データだけから両側で挟まれた**。 -/
theorem logDiff_tower_le_sum_of_bounded (F K : Type) [Field F] [NumberField F]
    [Field K] [NumberField K] [Algebra F K] [DecidableEq (Ideal (𝓞 K))]
    (S : Finset (Ideal (𝓞 K))) (hne : ∀ P ∈ S, P ≠ ⊥) (a : Ideal (𝓞 K) → ℕ)
    (hI : differentIdeal (𝓞 F) (𝓞 K) ≠ 0)
    (hprodne : (∏ P ∈ S, P ^ (a P)) ≠ 0)
    (hsupp : ∀ Q : Ideal (𝓞 K), Irreducible Q → Q ∣ differentIdeal (𝓞 F) (𝓞 K) → Q ∈ S)
    (hbound : ∀ P ∈ S, ¬ (P ^ (a P + 1) ∣ differentIdeal (𝓞 F) (𝓞 K))) :
    logDiffOfField K - logDiffOfField F
      ≤ (∑ P ∈ S, (a P : ℝ) * Real.log (Ideal.absNorm P)) / (Module.finrank ℚ K : ℝ) :=
  logDiff_tower_le_sum F K S hne a hI
    (dvd_prod_pow_of_bounded S a _ hI hprodne hsupp hbound)

/-! ## ★出典の紐付け(`.src`) -/

def pow_dvd_iff_le_count.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 9,
    item := "Proposition 1.7(冪の割り切りは重複度で測れる)",
    sectionId := "genell-prop-1-7" }

def dvd_prod_pow_of_bounded.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 9,
    item := "Proposition 1.7(台と重複度の上界から割り切りが出る)",
    sectionId := "genell-prop-1-7" }

def logDiff_tower_le_sum_of_bounded.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 9,
    item := "Proposition 1.7(分岐データだけから log-diff の上からの評価)",
    sectionId := "genell-prop-1-7" }

def logDiff_tower_le_sum_of_bounded.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]" "pow_dvd_iff_le_emultiplicity / emultiplicity_eq_count_normalizedFactors"
      (.inMathlib "UniqueFactorizationMonoid.emultiplicity_eq_count_normalizedFactors") 2,
    .citation "[mathlib]" "dvd_iff_normalizedFactors_le_normalizedFactors"
      (.inMathlib "UniqueFactorizationMonoid.dvd_iff_normalizedFactors_le_normalizedFactors") 2,
    .citation "[ABC3]" "logDiff_tower_le_sum(数値の側、§9-955)"
      (.inProject "ABC3" "ABC3.Found.GenEll.logDiff_tower_le_sum") 2,
    .implicitStep
      ("★★★★測定(2026-08-29): §9-955 が受けていた割り切り 𝔡 ∣ ∏ P^{a_P} は" ++
       "**重複度を数えるだけ**で出る。積の因子分解を計算する必要はない" ++
       "——Q ∈ S の側は Finset.dvd_prod_of_mem で Q^{a_Q} ∣ ∏ を出せば足りる") 5,
    .implicitStep
      ("★★残るのは hsupp(𝔡 の台が S に入る)と hbound(P^{e_P} ∤ 𝔡)を" ++
       "**分岐の言葉から供給する**段である。" ++
       "hbound は馴分岐なら mathlib の not_dvd_differentIdeal_of_isCoprime、" ++
       "hsupp は不分岐なら P ∤ 𝔡(not_dvd_differentIdeal_iff の系)。" ++
       "★したがって残るのは**野分岐の素点の扱い**だけであり、" ++
       "それが原文の Σ の上の食い違いにあたる") 4 ]

end ABC3.Found.GenEll
