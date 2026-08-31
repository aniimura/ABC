/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.DifferentRamified
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★★★★★★★`hup` は分岐指数だけで出る（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.9。

原文 (GenEll p.9):
> Proposition 1.7. (Conductors and Log Differents) Let

## ★★★★★★★★★★★★★★★★★★★★★これは何か —— `hup` が仮定なしになった

`§9-954` は「各 `P` で `P^{a_P} ∣ 𝔡`」を仮定として受けていた。
★その仮定は **mathlib の `pow_sub_one_dvd_differentIdeal`**
（`P^e ∣ p·𝓞_K` ⟹ `P^{e−1} ∣ 𝔡`）で**そのまま供給できる**。

★★★★到達:

    各 `P ∈ S` で `P^{e_P} ∣ p_P·𝓞_K`（＝分岐指数の定義そのもの）
      ⟹ `(∑_P (e_P − 1)·log N(P)) / [K:ℚ]  ≤  log-diff(K) − log-diff(F)`

★★これで `Proposition 1.7` の **`hup` は分岐指数だけで出る**——
仮定は「`e_P` が分岐指数である」以外に何も要らない。

## ★★★`hlow` との非対称（実測）

| 向き | 必要なもの | 状態 |
|---|---|---|
| `hup`（下から） | `P^{e−1} ∣ 𝔡` | ★**mathlib にある**（本ファイルで供給） |
| `hlow`（上から） | `P^{e} ∤ 𝔡`（馴分岐） | ☆**mathlib には不分岐（`e = 1`）しか無い** |

★★★★★**2026-08-29 実測**: 本プロジェクトの馴分岐 6/6
（`TameRamification` / `DifferentKummer` / `TotallyRamified`）が持つのも
`λ^m ∈ 𝔡`（`m ≥ e−1`）——**下からの側**である。
★上からの側（`P^e ∤ 𝔡`）は**まだ誰も持っていない**。
それが `Proposition 1.7` の elementary claim に残る最後の 1 本である。
-/

namespace ABC3.Found.GenEll

open NumberField

/-! ## ★★★★★★★★★★★★★★★★★★★★★分岐指数からの下からの評価 -/

/-- ★★★★★★★★★★★★★★★★★★★★★**`hup` は分岐指数だけで出る**。

原文 (GenEll p.9):
> Proposition 1.7. (Conductors and Log Differents) Let

    各 `P ∈ S` で `P^{e_P} ∣ p_P·𝓞_K`
      ⟹ `(∑_P (e_P − 1)·log N(P)) / [K:ℚ]  ≤  log-diff(K) − log-diff(F)`

★仮定は「`e_P` が分岐指数である」（`P^{e_P} ∣ p_P·𝓞_K`）以外に何も要らない。
★★`§9-954` の `hdvd` を mathlib の `pow_sub_one_dvd_differentIdeal` で埋めただけである。 -/
theorem sum_log_le_logDiff_of_ramificationIdx (F K : Type) [Field F] [NumberField F]
    [Field K] [NumberField K] [Algebra F K]
    (S : Finset (Ideal (𝓞 K))) (hmax : ∀ P ∈ S, P.IsMaximal) (hne : ∀ P ∈ S, P ≠ ⊥)
    (p : Ideal (𝓞 K) → Ideal (𝓞 F)) (e : Ideal (𝓞 K) → ℕ)
    (hpmax : ∀ P ∈ S, (p P).IsMaximal) (hpne : ∀ P ∈ S, p P ≠ ⊥)
    (hdvd : ∀ P ∈ S, P ^ (e P) ∣ (p P).map (algebraMap (𝓞 F) (𝓞 K)))
    (hI : differentIdeal (𝓞 F) (𝓞 K) ≠ 0) :
    (∑ P ∈ S, ((e P - 1 : ℕ) : ℝ) * Real.log (Ideal.absNorm P))
        / (Module.finrank ℚ K : ℝ)
      ≤ logDiffOfField K - logDiffOfField F := by
  refine sum_log_le_logDiff_tower F K S hmax hne (fun P => e P - 1) ?_ hI
  intro P hP
  haveI := hpmax P hP
  exact pow_sub_one_dvd_differentIdeal (𝓞 F) P (e P) (hpne P hP) (hdvd P hP)

/-! ## ★出典の紐付け(`.src`) -/

def sum_log_le_logDiff_of_ramificationIdx.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 9,
    item := "Proposition 1.7(hup は分岐指数だけで出る)",
    sectionId := "genell-prop-1-7" }

def sum_log_le_logDiff_of_ramificationIdx.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]" "pow_sub_one_dvd_differentIdeal(P^e ∣ p·𝓞_K ⟹ P^{e−1} ∣ 𝔡)"
      (.inMathlib "pow_sub_one_dvd_differentIdeal") 3,
    .citation "[ABC3]" "sum_log_le_logDiff_tower(局所から大域へ、§9-954)"
      (.inProject "ABC3" "ABC3.Found.GenEll.sum_log_le_logDiff_tower") 3,
    .implicitStep
      ("★★★★★測定(2026-08-29): Proposition 1.7 の hup(下から)は" ++
       "**分岐指数だけで出る**——mathlib の pow_sub_one_dvd_differentIdeal が" ++
       "そのまま仮定を埋める。仮定は「e_P が分岐指数である」以外に何も要らない") 5,
    .implicitStep
      ("★★★★★非対称の実測: hlow(上から)に要る P^e ∤ 𝔡(馴分岐)は" ++
       "**mathlib には不分岐(e = 1)しか無い**。" ++
       "本プロジェクトの馴分岐 6/6 が持つのも λ^m ∈ 𝔡(m ≥ e−1)——下からの側である。" ++
       "★上からの側は**まだ誰も持っていない**" ++
       "——それが Proposition 1.7 の elementary claim に残る最後の 1 本である") 6 ]

end ABC3.Found.GenEll
