/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.DifferentDivides
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★★★★★★分岐の言葉で `log-diff` を押さえる（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.9。

原文 (GenEll p.9):
> Proposition 1.7. (Conductors and Log Differents) Let

## ★★★★★★★★★★★★★★★★★★★★これは何か —— `hsupp` を分岐に翻訳する

`§9-956` は「`𝔡` の既約因子がすべて `S` に入る」（`hsupp`）を仮定として受けていた。
★mathlib の **`dvd_differentIdeal_iff`**（`P ∣ 𝔡 ↔ P は分岐している`）で、
それは**「分岐する素点はすべて `S` に入る」**に翻訳される。

★★★★到達:

    分岐する素点がすべて `S` に入り、各 `P ∈ S` で `P^{a_P+1} ∤ 𝔡`
      ⟹ `log-diff(K) − log-diff(F) ≤ (∑_P a_P·log N(P)) / [K:ℚ]`

★★これで `Proposition 1.7` の `hlow` は**分岐の言葉だけ**になった。

## ★★★系 —— 不分岐なら `log-diff` は増えない

`S = ∅` と取れば

    **すべての素点で不分岐 ⟹ `log-diff(K) ≤ log-diff(F)`**

★これは仮定の要らない**完全な定理**であり、
`§9-954`〜`§9-956` の鎖が空虚でないことの witness でもある。

## ★残っている段（明示）

★★★`hbound`（`P^{a_P+1} ∤ 𝔡`）だけが残る。
`a_P ≔ e_P − 1` と取ったときの `P^{e_P} ∤ 𝔡` は**馴分岐の主張**であり、
mathlib には**不分岐の場合（`e_P = 1`）しか無い**
（`not_dvd_differentIdeal_of_isCoprime`、2026-08-29 実測）。
★★本プロジェクトの `TameRamification.lean` / `DifferentKummer.lean` /
`TotallyRamified.lean`（馴分岐 6/6）が持っている**局所の主張を大域へ上げる**のが次である。
-/

namespace ABC3.Found.GenEll

open NumberField

/-! ## ★★★★★★★★★★★★★★★★★★★★分岐の言葉での上からの評価 -/

/-- ★★★★★★★★★★★★★★★★★★★★**分岐する素点を集めれば `log-diff` は上から押さえられる**。

原文 (GenEll p.9):
> Proposition 1.7. (Conductors and Log Differents) Let

★mathlib の `dvd_differentIdeal_iff`（`P ∣ 𝔡 ↔ P は分岐`）で
`§9-956` の `hsupp` を翻訳しただけである。
★★イデアルの既約性から素性への移行は Dedekind 環の一意分解による。 -/
theorem logDiff_tower_le_sum_of_ramified (F K : Type) [Field F] [NumberField F]
    [Field K] [NumberField K] [Algebra F K] [DecidableEq (Ideal (𝓞 K))]
    (S : Finset (Ideal (𝓞 K))) (hne : ∀ P ∈ S, P ≠ ⊥) (a : Ideal (𝓞 K) → ℕ)
    (hI : differentIdeal (𝓞 F) (𝓞 K) ≠ 0)
    (hprodne : (∏ P ∈ S, P ^ (a P)) ≠ 0)
    (hram : ∀ Q : Ideal (𝓞 K), ∀ _ : Q.IsPrime,
      ¬ Algebra.IsUnramifiedAt (𝓞 F) Q → Q ∈ S)
    (hbound : ∀ P ∈ S, ¬ (P ^ (a P + 1) ∣ differentIdeal (𝓞 F) (𝓞 K))) :
    logDiffOfField K - logDiffOfField F
      ≤ (∑ P ∈ S, (a P : ℝ) * Real.log (Ideal.absNorm P)) / (Module.finrank ℚ K : ℝ) := by
  refine logDiff_tower_le_sum_of_bounded F K S hne a hI hprodne ?_ hbound
  intro Q hQirr hQdvd
  have hp : Prime Q := UniqueFactorizationMonoid.irreducible_iff_prime.mp hQirr
  haveI hQP : Q.IsPrime := (Ideal.prime_iff_isPrime hp.ne_zero).mp hp
  exact hram Q hQP (dvd_differentIdeal_iff.mp hQdvd)

/-! ## ★★★★★★★★★系 —— 不分岐なら `log-diff` は増えない -/

/-- ★★★★★★★★★**すべての素点で不分岐なら `log-diff(K) ≤ log-diff(F)`**。

原文 (GenEll p.9):
> Proposition 1.7. (Conductors and Log Differents) Let

★`S = ∅` と取っただけである——仮定は「不分岐」だけで、**完全な定理**である。
★★`§9-954`〜`§9-956` の鎖が空虚でないことの witness でもある
（`Check/` の非空虚性検査と同じ役割）。 -/
theorem logDiff_le_of_unramified (F K : Type) [Field F] [NumberField F]
    [Field K] [NumberField K] [Algebra F K] [DecidableEq (Ideal (𝓞 K))]
    (hI : differentIdeal (𝓞 F) (𝓞 K) ≠ 0)
    (hunram : ∀ Q : Ideal (𝓞 K), ∀ _ : Q.IsPrime, Algebra.IsUnramifiedAt (𝓞 F) Q) :
    logDiffOfField K ≤ logDiffOfField F := by
  have h := logDiff_tower_le_sum_of_ramified F K ∅ (by simp) (fun _ => 0) hI (by simp)
    (fun Q hQ hnu => absurd (hunram Q hQ) hnu) (by simp)
  simpa using h

/-! ## ★出典の紐付け(`.src`) -/

def logDiff_tower_le_sum_of_ramified.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 9,
    item := "Proposition 1.7(分岐する素点を集めれば log-diff は上から押さえられる)",
    sectionId := "genell-prop-1-7" }

def logDiff_le_of_unramified.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 9,
    item := "Proposition 1.7(不分岐なら log-diff は増えない)",
    sectionId := "genell-prop-1-7" }

def logDiff_tower_le_sum_of_ramified.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]" "dvd_differentIdeal_iff(P ∣ 𝔡 ↔ P は分岐している)"
      (.inMathlib "IsDedekindDomain.dvd_differentIdeal_iff") 2,
    .citation "[ABC3]" "logDiff_tower_le_sum_of_bounded(重複度から割り切り、§9-956)"
      (.inProject "ABC3" "ABC3.Found.GenEll.logDiff_tower_le_sum_of_bounded") 3,
    .implicitStep
      ("★★★★測定(2026-08-29): §9-956 の hsupp(𝔡 の既約因子が S に入る)は" ++
       "mathlib の dvd_differentIdeal_iff で**分岐の言葉に翻訳できる**。" ++
       "★これで Proposition 1.7 の hlow は分岐の言葉だけになった") 5,
    .implicitStep
      ("★★★残るのは hbound(P^{a_P+1} ∤ 𝔡)だけである。" ++
       "a_P ≔ e_P − 1 のときの P^{e_P} ∤ 𝔡 は**馴分岐の主張**であり、" ++
       "mathlib には**不分岐の場合(e_P = 1)しか無い**" ++
       "(not_dvd_differentIdeal_of_isCoprime、2026-08-29 実測)。" ++
       "★本プロジェクトの TameRamification / DifferentKummer / TotallyRamified" ++
       "(馴分岐 6/6)が持つ局所の主張を大域へ上げるのが次である") 5 ]

end ABC3.Found.GenEll
