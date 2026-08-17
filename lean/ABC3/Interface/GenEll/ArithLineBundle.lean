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
その `APic(Spec 𝒪_F)` は原文自身の同型 `ADiv(F)/APrc(F) ⥲ APic(Spec(𝒪_F))` によって
**`ADiv(F)/APrc(F)`(数体の素点だけで書ける)** と同一視できる。

★**したがって境界をここに引く**——`Interface` が受けるのは
「点 `x` と算術直線束 `M̄` から、`x` を定義する数体 `F` ごとに正規化次数を作る操作」であり、
下流(`deg`・高さ・BD-class)は `Found/GenEll/ArithDiv.lean`(層 A)だけで書ける。

## ★何が無いのか(2026-08-16 実測、探索範囲つき)

| 層 | 要るもの | mathlib |
|---|---|---|
| B | スキーム上の直線束・可逆層・`APic(X)` | ★**0 件**(`LineBundle` / `invertible sheaf` / `Scheme…Picard`) |
| C | `X^arc`(スキームの解析化)・hermitian 計量・`ι_X` | ★**0 件**(`analytification` / `GAGA` / `complex analytic space`) |

★**層 C が「もう一段下」である**——hermitian 計量を載せる先の `X^arc` そのものが無い。

## ★★仮説を強めない工夫(`PLAN.md` A6)

原文の不変性は「`x` が `F` でも `K` でも定義され、`K ⊇ F`」という状況の主張である。
**それを `∀ F K` の無条件な等式として posit すると、原文より強い仮説になる**
(`check.mjs` 冒頭 A6 が警告する向き)。ゆえに

- `DefinedOver F x`(`x` が `F` 上で定義される)という述語を持たせ、
- `[Algebra F K]`(`K` が `F` の拡大)を要求し、
- 不変性も高さの一致も **`DefinedOver` の下でだけ**主張する。

★これで `Interface` は「原文が言っている以上のこと」を言わない。

## ★`waiting` を置く理由(退化 witness を作らない)

`PilotObjects.lean` と同じ判断である。**退化 witness(すべてを `0` に送る操作)なら
今すぐ作れる**が、それで G2 を満たしても**作業キューから消えるだけで何も進まない**
(`tools/check.mjs` 冒頭 B5 の穴)。ゆえに `nonvacuous` ではなく `waiting` を置く。
-/

namespace ABC3.Interface.GenEll

open ABC3.Meta NumberField

/-- **算術直線束から高さを作る操作**を受ける `Interface`。

★これは「算術直線束そのもの」ではなく、**高さの定義が実際に使うものだけ**である。
`X` や `M̄` の中身(層 B・C)は型に出さない——出すには解析化と可逆層が要り、それは無い。 -/
structure PulledBackClassData.{u, v} where
  /-- `X(ℚ̄)` の点を表す型。

  ★**宇宙は多相にしてある。**`Found/` の witness では点の型が
  「定義体 `F : Type` を含む構造の商」なので `Type 1` に居る
  (`Found/GenEll/UPoint.lean` の `UPoint`)。 -/
  Point : Type u
  /-- `X` 上の算術直線束 `M̄` を表す型。 -/
  Bundle : Type v
  /-- 点 `x` が数体 `F` の上で定義されること(原文 `x ∈ X(F)`)。 -/
  DefinedOver : (F : Type) → [Field F] → [NumberField F] → Point → Prop
  /-- `x` を定義する数体 `F` の上での正規化次数 `deg_F(x_F^* M̄)`。 -/
  degOver : (F : Type) → [Field F] → [NumberField F] → Bundle → Point → ℝ
  /-- 原文 p.4 の「`deg_K(L̄|_{Spec 𝒪_K}) = deg_F(L̄)`」。

  ★**この不変性が高さを `X(ℚ̄)` 上で well-defined にする**。飾りではない——
  これが破れたら `ht` は定義できない。
  ★`DefinedOver` と `[Algebra F K]` で条件を付けているのは、
  無条件にすると原文より強い仮説になるからである(上の docstring)。 -/
  base_change_invariant :
    ∀ (F K : Type) [Field F] [NumberField F] [Field K] [NumberField K] [Algebra F K]
      (M : Bundle) (x : Point),
      DefinedOver F x → DefinedOver K x → degOver K M x = degOver F M x
  /-- 高さ `ht_M̄(x)`——`F` の取り方に依らない値。 -/
  height : Bundle → Point → ℝ
  /-- `height` が、`x` を定義するどの `F` の上での `degOver` とも一致すること。

  ★**これが「`ht` が定義できる」ことの中身**である。 -/
  height_eq :
    ∀ (F : Type) [Field F] [NumberField F] (M : Bundle) (x : Point),
      DefinedOver F x → height M x = degOver F M x

/-- Track B は何を作らねばならないか。

★**退化 witness(すべてを `0` に送る)なら今すぐ作れる**が、それで G2 を満たしても
作業キューから消えるだけで何も進まない。ゆえに `waiting` を置く。 -/
def PulledBackClassData.waiting : WaitingFor :=
  { what := "算術直線束(Definition 1.1, (i))とその引き戻し(同 (ii))——すなわち スキーム上の直線束(層 B)と、X^arc(スキームの解析化)上の hermitian 計量(層 C)"
    trackB := "Found/GenEll — ただし層 C(解析化・GAGA)は mathlib に 0 件で、これは Arakelov 理論の**さらに下**の層である。層 A(ADiv・deg_F・ord_v・主因子)は Found/GenEll/ArithDiv.lean に実装済み(sorry 無し)" }

/-! ## ★出典の紐付け(`.src`) -/

def PulledBackClassData.src : Source :=
  { paper := "GenEll", pdfPage := 3, item := "Definition 1.1, (i)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Interface.GenEll
