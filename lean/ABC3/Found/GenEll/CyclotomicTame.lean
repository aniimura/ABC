/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import Mathlib.NumberTheory.Cyclotomic.Gal
import ABC3.Found.GenEll.PGroupChain
import ABC3.Meta.Claim

/-!
# ★★★★★★★段 EC5 —— `ζ_p` を添加する拡大は馴である（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.10。

原文 (GenEll p.10):
> integer n such that for any finite Galois extension L/K of finite extensions

## ★★★★★★★これは何か —— elementary claim の段 EC5

原文の証明は

> By separating the extension L[ζ]/K, where ζ is a primitive p-th root of unity,
> into a composite of wildly ramified and tamely ramified extensions

と書く。★その「馴の部分」の中身は **`K(ζ_p)/K` の次数が `p−1` を割る**ことである
——`p` 進体の剰余標数は `p` なので、`p−1` を割る次数は `p` で割れず、したがって**馴**である。

## ★★★機構 —— Galois 群が `(ℤ/p)ˣ` に埋まる

★mathlib の `IsPrimitiveRoot.autToPow`（`Gal(L/K) →* (ZMod n)ˣ`）と
`IsPrimitiveRoot.autToPow_injective`（円分拡大なら単射）で

    `[K(ζ_n) : K] = |Gal| ∣ |(ℤ/n)ˣ|`

★★`n = p` なら右辺は `p − 1` である。

## ★★★★測定の記録

★★★★★`IsPrimitiveRoot.autToPow_injective` は
`Mathlib/NumberTheory/Cyclotomic/Gal.lean` にあるが、
**本プロジェクトはこれまで `Mathlib.NumberTheory.Cyclotomic.*` を一度も import していなかった**
（2026-08-28 実測——REPL の基準環境で `IsCyclotomicExtension` が名前空間ごと未知だった）。
★本ファイルが最初の import である。
-/

namespace ABC3.Found.GenEll

/-! ## ★★★★★★円分拡大の次数は `(ℤ/n)ˣ` の位数を割る -/

/-- ★★★★★★**円分拡大の次数は `|(ℤ/n)ˣ|` を割る**。

原文 (GenEll p.10):
> integer n such that for any finite Galois extension L/K of finite extensions

★`IsPrimitiveRoot.autToPow_injective`（円分拡大なら `Gal → (ZMod n)ˣ` は単射）と
`IsGalois.card_aut_eq_finrank` の合成である。 -/
theorem finrank_cyclotomic_dvd (n : ℕ) [NeZero n] (K : Type) [Field K] (L : Type) [Field L]
    {μ : L} (hμ : IsPrimitiveRoot μ n) [Algebra K L] [IsCyclotomicExtension {n} K L]
    [FiniteDimensional K L] [IsGalois K L] :
    Module.finrank K L ∣ Nat.card (ZMod n)ˣ := by
  rw [← IsGalois.card_aut_eq_finrank K L]
  exact Subgroup.card_dvd_of_injective _ (hμ.autToPow_injective K)

/-- ★★**素数の場合** —— `[K(ζ_p) : K]` は `p − 1` を割る。

★これが「`ζ_p` を添加しても馴である」ことの中身である
——剰余標数は `p` だが、次数は `p−1` を割るので `p` では割れない。 -/
theorem finrank_cyclotomic_prime_dvd (p : ℕ) [Fact p.Prime] (K : Type) [Field K]
    (L : Type) [Field L] {μ : L} (hμ : IsPrimitiveRoot μ p) [Algebra K L]
    [IsCyclotomicExtension {p} K L] [FiniteDimensional K L] [IsGalois K L] :
    Module.finrank K L ∣ p - 1 := by
  haveI : NeZero p := ⟨(Fact.out : p.Prime).ne_zero⟩
  have h := finrank_cyclotomic_dvd p K L hμ
  rwa [Nat.card_eq_fintype_card, ZMod.card_units_eq_totient p,
    Nat.totient_prime (Fact.out : p.Prime)] at h

/-- ★★★**したがって次数は `p` で割れない**（＝馴分岐の次数条件）。 -/
theorem not_dvd_finrank_cyclotomic_prime (p : ℕ) [Fact p.Prime] (K : Type) [Field K]
    (L : Type) [Field L] {μ : L} (hμ : IsPrimitiveRoot μ p) [Algebra K L]
    [IsCyclotomicExtension {p} K L] [FiniteDimensional K L] [IsGalois K L]
    (hpos : 0 < Module.finrank K L) :
    ¬ (p ∣ Module.finrank K L) := by
  intro hdvd
  have hle : p ≤ Module.finrank K L := Nat.le_of_dvd hpos hdvd
  have hdvd2 := finrank_cyclotomic_prime_dvd p K L hμ
  have hp1 : 0 < p - 1 := Nat.sub_pos_of_lt (Fact.out : p.Prime).one_lt
  have hle2 : Module.finrank K L ≤ p - 1 := Nat.le_of_dvd hp1 hdvd2
  omega

/-! ## ★出典の紐付け(`.src`) -/

def finrank_cyclotomic_dvd.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 10,
    item := "Proposition 1.7(elementary claim——円分拡大の次数は |(ℤ/n)ˣ| を割る)",
    sectionId := "genell-prop-1-7" }

def finrank_cyclotomic_prime_dvd.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 10,
    item := "Proposition 1.7(elementary claim——[K(ζ_p):K] は p−1 を割る)",
    sectionId := "genell-prop-1-7" }

def not_dvd_finrank_cyclotomic_prime.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 10,
    item := "Proposition 1.7(elementary claim——ζ_p の添加は馴である)",
    sectionId := "genell-prop-1-7" }

def not_dvd_finrank_cyclotomic_prime.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]" "IsPrimitiveRoot.autToPow_injective(円分拡大なら Gal → (ZMod n)ˣ は単射)"
      (.inMathlib "IsPrimitiveRoot.autToPow_injective") 1,
    .citation "[mathlib]" "IsGalois.card_aut_eq_finrank"
      (.inMathlib "IsGalois.card_aut_eq_finrank") 1,
    .implicitStep
      ("★★★★★測定: IsPrimitiveRoot.autToPow_injective は " ++
       "Mathlib/NumberTheory/Cyclotomic/Gal.lean にあるが、" ++
       "**本プロジェクトはこれまで Mathlib.NumberTheory.Cyclotomic.* を" ++
       "一度も import していなかった**(2026-08-28 実測)。本ファイルが最初の import である") 2,
    .implicitStep
      ("★★これが原文の「By separating the extension L[ζ]/K … into a composite of " ++
       "wildly ramified and tamely ramified extensions」の**馴の側**である" ++
       "——p 進体の剰余標数は p だが、次数は p−1 を割るので p では割れない") 2,
    .implicitStep
      ("★★★残るのは EC6(局所から大域へ)と、" ++
       "本補題を IsTameDegree(§9-TameRamification)に繋ぐ配管である") 3 ]

end ABC3.Found.GenEll
