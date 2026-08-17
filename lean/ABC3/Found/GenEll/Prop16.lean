import ABC3.Found.GenEll.ArcModel
import ABC3.Found.GenEll.UPoint

/-!
# [GenEll] Proposition 1.6 —— **導手は高さで抑えられる**(組み上げ)(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.9。

原文 (GenEll p.9):
> Proposition 1.6. (Conductor Bounded by the Height) Let D ⊆ X be an effective Cartier divisor,

## ★★★本日の部品を 1 つに組む

原文の主張は

    `log-cond_D ≲ ht_D`   (`U_X(ℚ̄)` の上、`ht_D ≝ ht_{O_X(D)}`)

★因子表示では **`ht_D` は `ArithCartier` の `divisor := D` に対する `htArith`** である
——`L = O_X(D)` という対応が**定義そのもの**になる。

## ★★組む部品(すべて本日取得)

| 部品 | 場所 |
|---|---|
| 非アルキメデス側(被約化で次数が減る) | `ConductorHeight.lean` |
| アルキメデス側の下界(正規化で定数が一様) | `ArchBound.lean` |
| `X^arc` のコンパクト性 | `ArcModel.lean` |
| 連続関数は下に有界 | `ArchBound.lean` |

★★★**射影モデルと Green 関数の連続性さえ与えれば、定数は自動で出る。**
-/

namespace ABC3.Found.GenEll

open AlgebraicGeometry CategoryTheory NumberField

open scoped LinearAlgebra.Projectivization

variable (F : Type) [Field F] [NumberField F]

variable {X : Scheme.{0}} {V : Type} [NormedAddCommGroup V] [NormedSpace ℂ V]
  [FiniteDimensional ℂ V]

/-! ## ★★★`Proposition 1.6` -/

/-- ★★★**[GenEll] Proposition 1.6**(Conductor Bounded by the Height)。

原文 (GenEll p.9):
> Proposition 1.6. (Conductor Bounded by the Height) Let D ⊆ X be an effective Cartier divisor,

    `log-cond_D(x) ≤ ht_D(x) + C`   (`C` は `F` にも `x` にも依らない)

★★★**射影モデルと Green 関数の連続性から、定数が自動で出る。**

★機構は 3 段:
1. `X^arc` はコンパクト(`ArcModel.compactSpace`)
2. 連続な Green 関数は下に有界(`ArcModel.exists_bound`)
3. 導手は高さ + 定数で抑えられる(`ArchBound.logCond_le_htArith_add`)

★★定数が `F` に依らないのは**正規化のおかげ**である
(`Σ_v mult v = [F:ℚ]`。`ArchBound.lean` の測定)。 -/
theorem prop_1_6 (M : ArcModel X V) [Nonempty (complexPoints X)]
    (D : ArithCartier X) (hg : @Continuous _ _ M.topology _ D.green) :
    ∃ C : ℝ, 0 ≤ C ∧
      ∀ (F : Type) [Field F] [NumberField F] (xF : specRingOfIntegers F ⟶ X),
        pullbackIdeal F D.divisor xF ≠ 0 →
        logCond F D.divisor xF ≤ htArith F D xF + C := by
  obtain ⟨C, hC, hlo, -⟩ := M.exists_bound D.green hg
  exact ⟨C, hC, fun F _ _ xF h => logCond_le_htArith_add F D xF C hC h hlo⟩

/-- ★★**BD-class の水準**(原文の `≲` そのもの)。

★★★定数が `F` に依らないので、**すべての数体の上で同時に**成り立つ。 -/
theorem prop_1_6_bdge (M : ArcModel X V) [Nonempty (complexPoints X)]
    (D : ArithCartier X) (hg : @Continuous _ _ M.topology _ D.green)
    (F : Type) [Field F] [NumberField F]
    (h : ∀ xF : specRingOfIntegers F ⟶ X, pullbackIdeal F D.divisor xF ≠ 0) :
    BDge (fun xF : specRingOfIntegers F ⟶ X => logCond F D.divisor xF)
      (fun xF => htArith F D xF) := by
  obtain ⟨C, hC, hlo, -⟩ := M.exists_bound D.green hg
  exact logCond_bdge_htArith_of_bddBelow F D C hC hlo h

/-! ## ★出典の紐付け(`.src`)

★条つきである。原文は `X` が ℤ-固有であることから射影埋め込みを得るが、
本ファイルは `ArcModel` として**与えられたものとして受けている**。 -/

def prop_1_6.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 9,
    item := "Proposition 1.6(射影モデルを与えられたものとして——ℤ-固有からの構成は含まない)",
    sectionId := "genell-prop-1-6" }

end ABC3.Found.GenEll
