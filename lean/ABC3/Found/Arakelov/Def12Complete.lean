/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Arakelov.Def11Complete
import ABC3.Found.GenEll.BDClass
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★[GenEll] Definition 1.2 —— **項目全体が揃った**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.5。

原文 (GenEll p.5):
> as the height function associated to the arithmetic line bundle

## ★★★★★★★★★★これは何か

原文 Definition 1.2 は 2 条からなる:

* **(i)** `ht_M̄ : X(ℚ̄) → ℝ` を「`M̄` に伴う**高さ関数**」と呼ぶ（**命名**である）。
* **(ii)** `F ⊆ X(ℚ̄)` を固定し、`≲_F` / `≳_F` / `≈_F` と **BD-類**を定める。

★(ii) は `BDClass.lean` で閉じていた（`BDle` / `BDge` / `BDeq` / `bdSetoid` /
`BDClass` / `BDClass.le` / `BDClass.ge` / `BDClass.eq_of_le_of_ge`）。
★★(i) は「高さ関数がある」ことが中身であり、**そこが空いていた**。

## ★★★★★★★★★2026-08-27 の取り下げと、それが解けた理由

`Found/GenEll/Def12Height.lean` は 2026-08-27 に本項目全体の `.src` を**下げた**。
理由は `Check/GenEll/Def12Degenerate.lean` が機械にかけたとおり:

    古い degFOf F x = -(∑_σ x.2.fn (embPoint F σ)) / [F:ℚ]

には **`x.1`（`Pic` 類）が現れない**——`ht_indep_pic` が `rfl` で落ちた。
★アルキメデス側しか見ない写像は `deg_F` ではない。

★★**いまはそこが埋まっている**。`§9-776`〜`§9-782` の

    degArithPre L s = degFinPre L s + archDeg L s
    degFinPre L s   = log #(Γ_pre(L) / Γ(X,⊤)·s)

は有限素点側を **`Γ_pre(L)` の商の位数として実際に測る**。
★★★`Check/GenEll/Def12NonDegenerate.lean` がその非退化性を機械にかけた
（`degFinPre_two_ne` / `exists_degFinPre_ne_zero`）——古い失敗形の**反証**である。

## ★★★★★★★主張 —— 原文の条と宣言の対応

| 原文 | 宣言 |
|---|---|
| (i) `ht_M̄ : X(ℚ̄) → ℝ`（高さ関数） | ★`heightFunction` ＝ `htMetricAlg`（`§9-803`） |
| (i) `x_F` の取り方に依らない | `htMetricU_baseChange`（`§9-799`）／`htMetricAny_respects` |
| (i) 加法性・等長不変性 | `htMetricAlg_mul` ／ `htMetricAlg_congr` |
| (ii) `≲_F` / `≳_F` / `≈_F` | `BDle` / `BDge` / `BDeq` |
| (ii) `≈` は同値関係 | `bdeq_equivalence` |
| (ii) BD-類 `[α]_F` | `bdSetoid` / `BDClass` / `BDClass.mk` |
| (ii) `α ≈ β ↔ α ≲ β ∧ α ≳ β` | `bdeq_iff` / `BDClass.eq_of_le_of_ge` |
| (ii) 記号を BD-類に適用する | `BDClass.le` / `BDClass.ge` |
| (i)＋(ii) `F ⊆ X(ℚ̄)` の上の `[ht_M̄]_F` | ★`heightOn` / `heightBDClass`（**本ファイル**） |

## ★★★★★逸脱の記録

### 1. `X(ℚ̄)` は「定義体つきの点を底変換で同一視した商」である

`AlgPointAnyClass X`（`AlgPointClass.lean`）。★原文は `ℚ̄` 上の点の集合として扱うが、
「`x` を生じさせる**任意の**射がすべて同じ値を与える」という原文の要求は
`htMetricAny_respects` が**定理として**満たしている。
★★`X(ℚ̄) ≃ AlgPointAnyClass X` の同一視は本項目には要らない。

### 2. `F ⊆ X(ℚ̄)` は `Set (AlgPointAnyClass X)` として取る

`heightOn M S : ↥S → ℝ`。★`BDClass` は台の型について多相なのでそのまま乗る。

### 3. 高さは `Definition 1.1` の設計に依る

`AInv X`（前層の水準）＋ `degAPicM`。`Definition 1.1` の逸脱の記録
（`Def11Complete.lean`）がそのまま引き継がれる。
-/

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory NumberField
open ABC3.Found.GenEll

/-! ## ★★★★★★★★★★(i) 高さ関数 -/

/-- ★★★★★★★★★★**[GenEll] Definition 1.2, (i)** —— `M̄` に伴う**高さ関数**。

原文 (GenEll p.5):
> as the height function associated to the arithmetic line bundle

    `ht_M̄ : X(ℚ̄) → ℝ`

★中身は `§9-803` の `htMetricAlg` である（[Szp] Prop 1.1 を**仮定しない**）。 -/
noncomputable abbrev heightFunction {X : Scheme.{0}} (M : AInv X) :
    AlgPointAnyClass X → ℝ := htMetricAlg M

/-! ## ★★★★★★★★(ii) と繋ぐ —— `F ⊆ X(ℚ̄)` の上の BD-類 -/

/-- ★★★★★★**原文の `F ⊆ X(ℚ̄)` へ制限した高さ関数**。

原文 (GenEll p.5):
> we shall refer to an equivalence class relative to this equivalence relation -/
noncomputable def heightOn {X : Scheme.{0}} (M : AInv X) (S : Set (AlgPointAnyClass X)) :
    S → ℝ := fun x => htMetricAlg M (x : AlgPointAnyClass X)

/-- ★★★★★★★★**原文の `[ht_M̄]_F`** —— 高さ関数の BD-類。

原文 (GenEll p.5):
> we shall refer to an equivalence class relative to this equivalence relation

★★これが (i) と (ii) を繋ぐ——§1 の後続（`Proposition 1.4` 以降）が消費する形である。 -/
noncomputable def heightBDClass {X : Scheme.{0}} (M : AInv X)
    (S : Set (AlgPointAnyClass X)) : BDClass S :=
  BDClass.mk (heightOn M S)

/-- ★★★★**高さは加法的である**（部分集合の上でも）。 -/
theorem heightOn_mul {X : Scheme.{0}} (M N : AInv X) (S : Set (AlgPointAnyClass X)) :
    heightOn (M.mul N) S = heightOn M S + heightOn N S := by
  funext x
  exact htMetricAlg_mul M N _

/-- ★★★**高さは等長同型類だけで決まる**（部分集合の上でも）。 -/
theorem heightOn_congr {X : Scheme.{0}} {M N : AInv X}
    (h : Isometric M.carrier N.carrier) (S : Set (AlgPointAnyClass X)) :
    heightOn M S = heightOn N S := by
  funext x
  exact htMetricAlg_congr h _

/-- ★★★★★**BD-類も等長同型類だけで決まる**。 -/
theorem heightBDClass_congr {X : Scheme.{0}} {M N : AInv X}
    (h : Isometric M.carrier N.carrier) (S : Set (AlgPointAnyClass X)) :
    heightBDClass M S = heightBDClass N S := by
  rw [heightBDClass, heightBDClass, heightOn_congr h S]

/-! ## ★出典の紐付け(`.src`) -/

def heightFunction.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 5,
    item := "Definition 1.2, (i)(ht_M̄ : X(ℚ̄) → ℝ を高さ関数と呼ぶ——[Szp] 不使用の形)",
    sectionId := "genell-def-1-2" }

def heightBDClass.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 5,
    item := "Definition 1.2(高さ関数の BD-類 [ht_M̄]_F——(i) と (ii) を繋ぐ形)",
    sectionId := "genell-def-1-2" }

/-- ★★★★★★★★★★**[GenEll] Definition 1.2**（高さ関数と BD-類）—— 実装された。

原文 (GenEll p.5):
> as the height function associated to the arithmetic line bundle

★対応表と逸脱の記録は本ファイル冒頭の docstring にある。

★★**2026-08-27 の取り下げとの違い**: 当時の `ht` は
`Check/GenEll/Def12Degenerate.lean` が示したとおり `Pic` 類を見ていなかった
（`ht_indep_pic` が `rfl` で落ちる）。★★★いまの `htMetricAlg` は
有限素点側を `log #(Γ_pre(L)/Γ(X,⊤)·s)` として実際に測っており、
その非退化性は `Check/GenEll/Def12NonDegenerate.lean` で機械にかけた。 -/
def definition_1_2_metric.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 5, item := "Definition 1.2",
    sectionId := "genell-def-1-2" }

def definition_1_2_metric.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "htMetricAlg(ht_M̄ : X(ℚ̄) → ℝ、§9-803)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.htMetricAlg") 5,
    .citation "[ABC3]" "htMetricU_baseChange(x_F の取り方に依らない、§9-799)"
      (.inProject "ABC3" "ABC3.Found.Arakelov.htMetricU_baseChange") 4,
    .citation "[ABC3]" "BDle / BDge / BDeq / bdeq_equivalence / BDClass((ii) の全体)"
      (.inProject "ABC3" "ABC3.Found.GenEll.BDClass") 5,
    .citation "[ABC3]" "BDClass.le / BDClass.ge(記号を BD-類に適用する段)"
      (.inProject "ABC3" "ABC3.Found.GenEll.BDClass.le") 5,
    .citation "[ABC3]" "degFinPre_two_ne(有限素点側が退化していないことの機械検査)"
      (.inProject "ABC3" "ABC3.Check.GenEll.degFinPre_two_ne") 5,
    .otherPaper "GenEll" "Definition 1.1" 3,
    .implicitStep
      ("★逸脱 1: X(ℚ̄) は定義体つきの点を底変換で同一視した商(AlgPointAnyClass)である。" ++
       "原文の「x を生じさせる任意の射がすべて同じ値を与える」は " ++
       "htMetricAny_respects が定理として満たしている") 5 ]

end ABC3.Found.Arakelov
