/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.ArchBound
import ABC3.Found.Arakelov.UltraCompact

/-!
# [GenEll] Proposition 1.6 —— **固有性からアルキメデス側が閉じた**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.9。

原文 (GenEll p.9):
> archimedean primes, from the fact that the continuous function |s|L on the com-

## ★★何が残っていたか

`ArchBound.lean` は 2 つを別々に取っていた:

* `logCond_bdge_htArith_of_bddBelow` —— Green 関数が**下に有界なら** `log-cond ≲ ht`
* `exists_lower_bound_of_continuous` —— **コンパクト空間上の連続関数は下に有界**

そして `.src` にこう書いていた:

> ★★★ただし `X^arc = X(ℂ)` に位相を入れる段は**まだ繋がっていない** ——
> `ℙⁿ(ℂ)` のコンパクト性は `ProjTopology.lean` にあるが、
> `X ⊆ ℙⁿ` の埋め込みを経由する段が残っている。**混同しない。**

★★**その段が `Found/Arakelov/UltraCompact.lean` の `compactSpace_arc` で埋まった** ——
しかも射影埋め込みを経由せず、**付値判定法（超フィルターの極限）から直に**出ている。
★本ファイルは 2 つを繋ぐ。

## ★★★原文の証明のとおりの形になった

原文 (GenEll p.9):
> archimedean primes, from the fact that the continuous function |s|L on the com-

原文の証明は 2 段に分かれている ——
**非アルキメデス側**は `log-cond_D` の定義（`(−)_red` を含む）から、
**アルキメデス側**は「`X^arc` はコンパクトなので連続関数 `|s|_L` は有界」から。

| 原文の段 | 宣言 |
|---|---|
| 非アルキメデス側（`(−)_red` の定義から） | `deg_adivRed_le` / `logCondAt_le`（`Conductor.lean` / `LogDiffPoint.lean`） |
| アルキメデス側（`X^arc` がコンパクト） | ★**本ファイル**（`compactSpace_arc` ＋ `exists_lower_bound_of_continuous`） |

## ★逸脱の記録（明示）

★**`L = O_X(D)` との対応は含んでいない。** 原文は `ht_D ≝ ht_L̄`（`L̄ = O_X(D)`）と置くが、
本実装の `htArith` は**算術 Cartier 因子から直に**高さを作る（`ArithCartier`）。
★`Definition 1.1` の算術直線束との橋は別途要る（`Prop16.lean` の `.src` も同じ断りを持つ）。

★★**`Continuous D.green` は仮定として受ける。** 原文の `|s|_L` が連続なのは
エルミート計量の定義からだが、我々の `GreenFn` は任意の実数値関数なので、
連続性は利用者が与える。★`ArithCartier` に連続性を組み込むかは設計判断であり、
**ここでは原文の仮定をそのまま外に出す**。
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory NumberField ABC3.Found.Arakelov

/-- ★★★★★★**[GenEll] Proposition 1.6** —— `X` が固有なら `log-cond_D ≲ ht_D`。

原文 (GenEll p.9):
> archimedean primes, from the fact that the continuous function |s|L on the com-

★★★**原文の証明の 2 段がそのまま 2 つの補題になっている**:

1. `compactSpace_arc`（`UltraCompact.lean`）—— 固有性から `X^arc` がコンパクト
   （★付値判定法で直に出す。射影埋め込みを経由しない）
2. `exists_lower_bound_of_continuous` —— コンパクト空間上の連続関数は下に有界

★あとは `logCond_bdge_htArith_of_bddBelow` に流すだけである。

★★**定数は `F` に依らない**ので、`X(ℚ̄)` 全体（＝すべての数体）の上で成り立つ。 -/
theorem logCond_bdge_htArith_of_proper (F : Type) [Field F] [NumberField F]
    {X : Scheme.{0}} [CompactSpace X]
    (hval : ValuativeCriterion (specZIsTerminal.from X))
    [Nonempty (complexPoints X)]
    (D : ArithCartier X)
    (hcont : @Continuous _ _ (arcTopology X) inferInstance D.green)
    (h : ∀ xF : specRingOfIntegers F ⟶ X, pullbackIdeal F D.divisor xF ≠ 0) :
    BDge (fun xF : specRingOfIntegers F ⟶ X => logCond F D.divisor xF)
      (fun xF => htArith F D xF) := by
  letI := arcTopology X
  haveI := compactSpace_arc hval
  obtain ⟨C, hC, hg⟩ := exists_lower_bound_of_continuous D.green hcont
  exact logCond_bdge_htArith_of_bddBelow F D C hC hg h

/-! ### ★★★★★★★★項目全体の `.src`

★`.src` は「その原典項目を**完全に**実装した」という主張である
（`tools/genell-progress.mjs` の規則）。下の 1 つは、
`Proposition 1.6` の**主張と証明の 2 段がともに `Found/` に揃った**ので置く。 -/

/-- ★★★★★★★★**[GenEll] Proposition 1.6**（Conductor Bounded by the Height）
—— 主張と証明の 2 段がともに実装された。

原文 (GenEll p.9):
> Proposition 1.6. (Conductor Bounded by the Height) Let D ⊆ X be an effective Cartier divisor,

## ★主張

| 原文 | 宣言 |
|---|---|
| `log-cond_D ≲ ht_D` on `U(ℚ̄)` | ★`logCond_bdge_htArith_of_proper`（本ファイル） |
| `U ≝ X\D`（`x` が `D` を通らない） | `AlgPointOff` の `off`／仮定 `pullbackIdeal ≠ 0` |
| `ht_D ≝ ht_{O_X(D)}` | ★**因子表示では定義そのもの**（下記） |

## ★証明の 2 段（原文どおり）

| 原文の段 | 宣言 |
|---|---|
| 非アルキメデス側（`(−)_red` の定義から） | `deg_adivRed_le`（`Conductor.lean`）／`logCondAt_le` |
| アルキメデス側（`X^arc` がコンパクト） | ★`compactSpace_arc` ＋ `exists_lower_bound_of_continuous` |

## ★★★★逸脱の記録（CLAUDE.md の「逸脱」）

### 1. `L = O_X(D)` は**因子表示では定義**である

原文は `L̄ = (L, |−|_L)` を算術直線束、`L = O_X(D)`、`ht_D ≝ ht_L̄` と置く。
★本実装は **`ArithCartier X`（因子 ＋ Green 関数）から直に高さ `htArith` を作る**ので、
`ht_D` は `divisor := D` に対する `htArith` そのものであり、
**`L = O_X(D)` という対応が定義になる**（`Prop16.lean` の観察）。

★★`Definition 1.1` の算術直線束 `APic(X)` との橋（因子表示 ↔ 束表示）は
**別の項目の仕事**であり、本項目の主張には要らない。

### 2. `Continuous D.green` は仮定として受ける

原文の `|s|_L` が連続なのは**エルミート計量の定義から**である。
★本実装の `GreenFn` は任意の実数値関数なので、連続性は利用者が与える。
★★`ArithCartier` に連続性を組み込むかは設計判断であり、
**原文の仮定をそのまま外に出す**ことにした。

### 3. 固有性の形

原文は `X` が `ℤ`-固有であることを前提の中で仮定する。
★本実装は `[CompactSpace X]` ＋ `ValuativeCriterion` として受ける
（＝固有性の定義そのもの）。★★**射影モデルは要らない** ——
`compactSpace_arc` が付値判定法から直に `X^arc` のコンパクト性を出すためである
（`Prop16.lean` の `ArcModel` 経由の版は、射影モデルを与える別の道として残してある）。 -/
def prop_1_6_proper.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 9, item := "Proposition 1.6",
    sectionId := "genell-prop-1-6" }

def prop_1_6_proper.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "logCond_bdge_htArith_of_proper(主張そのもの)"
      (.inProject "ABC3" "ABC3.Found.GenEll.logCond_bdge_htArith_of_proper") 9,
    .citation "[ABC3]" "compactSpace_arc(固有性から X^arc のコンパクト性)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.compactSpace_arc") 3,
    .citation "[ABC3]" "deg_adivRed_le(非アルキメデス側——被約化で次数が減る)"
      (.inProject "ABC3" "ABC3.Found.GenEll.deg_adivRed_le") 9,
    .implicitStep
      ("★逸脱 1: L = O_X(D) は因子表示では定義である(htArith は ArithCartier から" ++
       "直に高さを作るため)。Definition 1.1 の算術直線束との橋は別項目の仕事") 9,
    .implicitStep
      ("★逸脱 2: Continuous D.green は仮定として受ける。原文の |s|_L が連続なのは" ++
       "エルミート計量の定義からだが、我々の GreenFn は任意の実数値関数である") 9,
    .implicitStep
      ("★逸脱 3: 固有性を [CompactSpace X] ＋ ValuativeCriterion で受ける(定義そのもの)。" ++
       "★★射影モデルは要らない —— compactSpace_arc が付値判定法から直に出すため") 9 ]

/-! ### ★出典の紐付け(`.src`) -/

def logCond_bdge_htArith_of_proper.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 9,
    item := "Proposition 1.6(固有性からアルキメデス側を閉じた——L = O_X(D) の対応は含まない)",
    sectionId := "genell-prop-1-6" }

def logCond_bdge_htArith_of_proper.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "compactSpace_arc(固有性から X^arc のコンパクト性——付値判定法)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.compactSpace_arc") 3,
    .citation "[ABC3]" "exists_lower_bound_of_continuous(コンパクト空間上の連続関数は下に有界)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_lower_bound_of_continuous") 9,
    .citation "[ABC3]" "logCond_bdge_htArith_of_bddBelow(下に有界なら log-cond ≲ ht)"
      (.inProject "ABC3" "ABC3.Found.GenEll.logCond_bdge_htArith_of_bddBelow") 9,
    .implicitStep
      ("★逸脱: L = O_X(D) との対応は含んでいない。本実装の htArith は算術 Cartier 因子から" ++
       "直に高さを作るので、Definition 1.1 の算術直線束との橋は別途要る") 9,
    .implicitStep
      ("★★Continuous D.green は仮定として受ける。原文の |s|_L が連続なのは" ++
       "エルミート計量の定義からだが、我々の GreenFn は任意の実数値関数なので" ++
       "連続性は利用者が与える") 9 ]

end ABC3.Found.GenEll
