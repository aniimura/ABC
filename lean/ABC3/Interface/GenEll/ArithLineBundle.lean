import ABC3.Meta.Claim
import Mathlib.NumberTheory.NumberField.InfinitePlace.Basic

/-!
# [GenEll] §1 —— 算術直線束の `Interface`(層 B・C の境界)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、
物理 p.3–p.4。**260 dpi 目視確認 2026-08-16**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★この `Interface` が受けているもの

原文の高さの定義(p.4、目視)は

```
ht_M̄(x) ≝ deg_F(x_F^* M̄) ∈ ℝ
```

である。★**一般の `X` 上の算術直線束が要るのは、
`x_F^* M̄ ∈ APic(Spec 𝒪_F)` を作るためだけ**であり、
その `APic(Spec 𝒪_F)` は原文自身の同型

```
ADiv(F)/APrc(F) ⥲ APic(Spec(𝒪_F))
```

によって **`ADiv(F)/APrc(F)`(数体の素点だけで書ける)** と同一視できる。

★**したがって境界をここに引く**——`Interface` が受けるのは
「点 `x` と算術直線束 `M̄` から `ADiv(F)` の類を作る操作」1 つだけであり、
下流(`deg`・高さ・BD-class)は `Found/GenEll/ArithDiv.lean`(層 A)だけで書ける。

## ★何が無いのか(2026-08-16 実測、探索範囲つき)

| 層 | 要るもの | mathlib |
|---|---|---|
| B | スキーム上の直線束・可逆層・`APic(X)` | ★**0 件**(`LineBundle` / `invertible sheaf` / `Scheme…Picard`) |
| C | `X^arc`(スキームの解析化)・hermitian 計量・`ι_X` | ★**0 件**(`analytification` / `GAGA` / `complex analytic space`) |

★**層 C が「もう一段下」である**——hermitian 計量を載せる先の `X^arc` そのものが無い。
`ResearchPaper/foundations.json` の `arakelov` 節点はこの 2 層をまとめて指している。

## ★`waiting` を置く理由(退化 witness を作らない)

`PilotObjects.lean` と同じ判断である。**退化 witness(すべての類を `0` に送る操作)なら
今すぐ作れる**が、それで G2 を満たしても**作業キューから消えるだけで何も進まない**
(`tools/check.mjs` 冒頭 B5 の穴)。ゆえに `nonvacuous` ではなく `waiting` を置く。
-/

namespace ABC3.Interface.GenEll

open ABC3.Meta NumberField

/-- **算術直線束から算術因子の類を作る操作**を受ける `Interface`。

★これは「算術直線束そのもの」ではなく、**高さの定義が実際に使う 1 操作**である。
原文の `x_F^* M̄` を `ADiv(F)` の類として受け取る。

★**意図的に弱く取っている**——`X` も `M̄` も型に出さず、
「数体 `F` ごとに `ADiv(F)` の元が定まり、有限次拡大で正規化次数が保たれる」
という**下流が実際に使う性質だけ**を要求する。
`X` や `M̄` を型に出すには層 B・C が要り、それは無い。 -/
structure PulledBackClassData where
  /-- 数体 `F` と点 `x ∈ X(F)` に対する `x_F^* M̄` の類(の代表)。

  ★`ADiv(F)` そのものではなく `ℝ` 値の正規化次数を受けるのは、
  `APrc(F)` による商を型に出すと層 D(同型 `ADiv/APrc ≅ APic`)まで
  この `Interface` が抱え込むからである。**境界は薄く取る。** -/
  normalizedDegree : (F : Type) → [Field F] → [NumberField F] → (Type) → ℝ
  /-- 原文 p.4 の「`deg_K(L̄|_{Spec 𝒪_K}) = deg_F(L̄)`」。

  ★**この不変性が高さを `X(ℚ̄)` 上で well-defined にする**。飾りではない——
  これが破れたら `ht` は定義できない。 -/
  base_change_invariant :
    ∀ (F K : Type) [Field F] [NumberField F] [Field K] [NumberField K] (P : Type),
      normalizedDegree K P = normalizedDegree F P

/-- Track B は何を作らねばならないか。

★**退化 witness(すべてを `0` に送る)なら今すぐ作れる**が、それで G2 を満たしても
作業キューから消えるだけで何も進まない。ゆえに `waiting` を置く。 -/
def PulledBackClassData.waiting : WaitingFor :=
  { what := "算術直線束(Definition 1.1, (i))とその引き戻し(同 (ii))——すなわち スキーム上の直線束(層 B)と、X^arc(スキームの解析化)上の hermitian 計量(層 C)"
    trackB := "Found/GenEll — ただし層 C(解析化・GAGA)は mathlib に 0 件で、これは Arakelov 理論の**さらに下**の層である。層 A(ADiv・deg_F)は Found/GenEll/ArithDiv.lean に実装済み" }

/-! ## ★出典の紐付け(`.src`) -/

def PulledBackClassData.src : Source :=
  { paper := "GenEll", pdfPage := 3, item := "Definition 1.1, (i)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Interface.GenEll
