/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.ProjCompact
import ABC3.Meta.Claim

/-!
# ★★★★★★★超平面はチャートの上で `x_0/x_i` で切られる（片側）（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.6。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

## ★★★★★★★これは何か —— 段 C2c-1 の本体（片側）

`§9-857` で「チャートの上で超平面イデアルを読む」道具が入った。
★あとは `ker (Away.map (hyperplaneHom) (x_i))` が **`x_0/x_i` で生成される**ことを見ればよい。

★★本ファイルはその**片側**（`⊇`、すなわち `span {x_0/x_i} ≤ ker`）を取る:

| 補題 | 内容 |
|---|---|
| `hyperplaneHom_X_zero` | `hyperplaneHom (x₀) = 0`（定義そのもの） |
| `awayMap_projCoord_zero` | ★**`x_0/x_i ↦ 0`** |
| `span_projCoord_le_ker` | ★★`span {x_0/x_i} ≤ ker` |

## ★★★機構 —— mathlib の 2 本

| 道具 | 役割 |
|---|---|
| `Proj.awayι_comp_map`（mathlib） | チャートの上で `Proj.map` は `Away.map` である |
| `HomogeneousLocalization.Away.map_mk`（mathlib） | `Away.map g f (mk n x) = mk n (g x)` |

★あとは `hyperGen N R 0 = 0`（`Fin.cases 0 X` の第 0 成分）だけである。

## ★残っている段（明示）

★★★**逆向き（`ker ≤ span {x_0/x_i}`）が残る**。それには

> `A⁰_{x_i}` は `R` 上 `x_j/x_i`（`j ≠ i`）を変数とする**多項式環である**

ことが要る——`§9-850` の `ext_of_projCoord` は**生成**（全射性）の半分しか取っていない。
★★もう半分（**関係が無いこと**、すなわち自由性）が必要である。
★★★これは `Proj` の標準チャート `D₊(x_i) ≅ 𝔸ᴺ` の**残り半分**であり、mathlib に無い
（2026-08-28 実測）。
-/

namespace ABC3.Found.GenEll

open MvPolynomial AlgebraicGeometry CategoryTheory HomogeneousLocalization

attribute [local instance] MvPolynomial.gradedAlgebra

/-! ## ★超平面の定義そのもの -/

/-- ★**`hyperplaneHom (x₀) = 0`** —— 超平面 `{x₀ = 0}` の定義そのもの。 -/
theorem hyperplaneHom_X_zero (N : ℕ) (R : Type) [CommRing R] :
    hyperplaneHom N R (MvPolynomial.X (0 : Fin (N + 1))) = 0 := by
  rw [hyperplaneHom_apply, MvPolynomial.aeval_X]
  show hyperGen N R 0 = 0
  rw [hyperGen]
  simp

/-! ## ★★★★★★チャートの上で `x_0/x_i` は消える -/

/-- ★★★★★★**チャートの上で `x_0/x_i` は消える**。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

★`Proj.awayι_comp_map`（mathlib）により、チャート `D₊(x_i)` の上で
`Proj.map` は `Away.map (hyperplaneHom) (x_i)` である。
★★`Away.map_mk` で分子だけを見れば `hyperplaneHom (x₀) = 0` に落ちる。 -/
theorem awayMap_projCoord_zero (N : ℕ) (R : Type) [CommRing R] [Nontrivial R]
    (i : Fin (N + 1)) :
    Away.map (hyperplaneHom N R) (MvPolynomial.X i) (projCoord N R i 0) = 0 := by
  rw [projCoord, Away.map_mk]
  refine HomogeneousLocalization.val_injective _ ?_
  rw [Away.val_mk, HomogeneousLocalization.val_zero]
  show Localization.mk (hyperplaneHom N R (MvPolynomial.X (0 : Fin (N+1)))) _ = 0
  rw [hyperplaneHom_X_zero, Localization.mk_zero]

/-- ★★**したがって `span {x_0/x_i} ≤ ker`** —— 段 C2c-1 の片側。

★逆向きには「`A⁰_{x_i}` は `x_j/x_i` を変数とする**多項式環である**」ことが要る
——`§9-850` は**生成**の半分しか取っていない（残りは自由性）。 -/
theorem span_projCoord_le_ker (N : ℕ) (R : Type) [CommRing R] [Nontrivial R]
    (i : Fin (N + 1)) :
    Ideal.span {projCoord N R i 0}
      ≤ RingHom.ker (Away.map (hyperplaneHom N R) (MvPolynomial.X i)) := by
  rw [Ideal.span_le]
  rintro z rfl
  exact awayMap_projCoord_zero N R i

/-! ## ★出典の紐付け(`.src`) -/

def hyperplaneHom_X_zero.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(hyperplaneHom (x₀) = 0)",
    sectionId := "genell-prop-1-4" }

def awayMap_projCoord_zero.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(チャートの上で x_0/x_i は消える)",
    sectionId := "genell-prop-1-4" }

def span_projCoord_le_ker.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(段 C2c-1 の片側——span {x_0/x_i} ≤ ker)",
    sectionId := "genell-prop-1-4" }

def span_projCoord_le_ker.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]" "Proj.awayι_comp_map(チャートの上で Proj.map は Away.map)"
      (.inMathlib "AlgebraicGeometry.Proj.awayι_comp_map") 2,
    .citation "[mathlib]" "HomogeneousLocalization.Away.map_mk"
      (.inMathlib "HomogeneousLocalization.Away.map_mk") 2,
    .citation "[ABC3]" "hyperplaneIdeal_apply(チャートの上で超平面イデアルを読む、§9-857)"
      (.inProject "ABC3" "ABC3.Found.GenEll.hyperplaneIdeal_apply") 2,
    .implicitStep
      ("★★★**逆向き(ker ≤ span {x_0/x_i})が残る**。それには " ++
       "「A⁰_{x_i} は R 上 x_j/x_i(j ≠ i)を変数とする**多項式環である**」ことが要る" ++
       "——§9-850 の ext_of_projCoord は**生成**(全射性)の半分しか取っていない。" ++
       "★もう半分(**関係が無いこと**、すなわち自由性)が必要である。" ++
       "★★これは Proj の標準チャート D₊(x_i) ≅ 𝔸ᴺ の**残り半分**であり mathlib に無い" ++
       "(2026-08-28 実測)") 6 ]

end ABC3.Found.GenEll
