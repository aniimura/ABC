/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import Mathlib.NumberTheory.NumberField.Completion.FinitePlace
import ABC3.Found.GaloisRep.DegInfTateParam

/-!
# 界面の測定 —— **完備化の橋に欠けているのは `IsPrecomplete` 1 つだけ**（`Check`）

**これは原典の主張ではない**（我々のモデルについての事実）ので `.src` を持たない。

## ★★★★★★★★★★ 2026-08-31 の測定（第 896）

`Lemma 3.5` の残る道は、大域の曲線 `E/L` を各悪い素点 `p` で**完備化**に落とし、
`minDeltaExp_eq_mul_of_tateParamR`（第 892）に渡すことである。
そのためには具体的な `R`・`Lv` を選ばなければならない。

★自然な選択は

    `R := v.adicCompletionIntegers K`   `Lv := v.adicCompletion K`

である。**何が揃っていて何が欠けているかを `#synth` で測った**:

| インスタンス | 状態 |
|---|---|
| `IsDiscreteValuationRing R` | ★**ある**（`NumberField/Completion/FinitePlace.lean`） |
| `IsDomain R`・`IsNoetherianRing R`・`IsLocalRing R` | ★ある |
| `IsFractionRing R Lv` | ★ある |
| `Algebra K Lv`・`CompleteSpace Lv` | ★ある |
| `IsHausdorff (maximalIdeal R) R` | ★**ある** |
| **`IsPrecomplete (maximalIdeal R) R`** | ☆☆☆**無い** |
| `IsAdicComplete (maximalIdeal R) R` | ☆無い（上の 1 つが欠けているためだけ） |
| `CompactSpace R`・`CompleteSpace R` | ☆無い |
| `LocallyCompactSpace Lv`・`ValuativeRel Lv` | ☆無い（`IsNonarchimedeanLocalField` 経由は使えない） |

★★★★★★★★**つまり、完備化の橋はインスタンス 1 つに落ちている**。

    `IsPrecomplete (IsLocalRing.maximalIdeal (v.adicCompletionIntegers K))`
    `              (v.adicCompletionIntegers K)`

## ☆道筋

mathlib には `IsNonarchimedeanLocalField K` の下で `IsAdicComplete 𝓂[K] 𝒪[K]` がある
（`NumberTheory/LocalField/Basic.lean`）が、`v.adicCompletion K` をそのクラスに
登録するインスタンスは無い（`ValuativeRel` も `LocallyCompactSpace` も合成できない）。

☆したがって直接作るのが早い。段は 3 つ:

1. `R` は `Lv` の**閉集合**（`Valued.v x ≤ 1`）なので `CompleteSpace R`
2. `SMOD I^n` を付値の不等式に翻訳し、仮説の列が Cauchy であることを見る
3. 極限 `L` を取り、`I^n` が閉じていることから `f n ≡ L [SMOD I^n]`

## ★以下は「ある」側の実測（このファイルが通ることが証拠）
-/

namespace ABC3.Check.GenEll

open IsDedekindDomain NumberField IsLocalRing

/-- ★完備化の整数環は離散付値環である（mathlib にある）。 -/
example (K : Type) [Field K] [NumberField K] (v : HeightOneSpectrum (𝓞 K)) :
    IsDiscreteValuationRing (v.adicCompletionIntegers K) := inferInstance

/-- ★分数体の関係もある。 -/
example (K : Type) [Field K] [NumberField K] (v : HeightOneSpectrum (𝓞 K)) :
    IsFractionRing (v.adicCompletionIntegers K) (v.adicCompletion K) := inferInstance

/-- ★Noether で局所である。 -/
example (K : Type) [Field K] [NumberField K] (v : HeightOneSpectrum (𝓞 K)) :
    IsNoetherianRing (v.adicCompletionIntegers K) := inferInstance

example (K : Type) [Field K] [NumberField K] (v : HeightOneSpectrum (𝓞 K)) :
    IsLocalRing (v.adicCompletionIntegers K) := inferInstance

/-- ★★**Hausdorff の側はある**——欠けているのは `IsPrecomplete` だけである。 -/
example (K : Type) [Field K] [NumberField K] (v : HeightOneSpectrum (𝓞 K)) :
    IsHausdorff (maximalIdeal (v.adicCompletionIntegers K))
      (v.adicCompletionIntegers K) := inferInstance

/-- ★体の側は完備である（部分環の側が無いだけ）。 -/
example (K : Type) [Field K] [NumberField K] (v : HeightOneSpectrum (𝓞 K)) :
    CompleteSpace (v.adicCompletion K) := inferInstance

end ABC3.Check.GenEll
