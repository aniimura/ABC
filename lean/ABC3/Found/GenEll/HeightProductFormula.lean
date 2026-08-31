/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.HeightIntegralRepr
import ABC3.Meta.Claim

/-!
# ★★★★★★積公式の帰結と正規化 —— 段 C2c-1 が消費する 2 本（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.6。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

## ★★★★★★これは何か —— 段 C2c-1 の材料

段 C2c-1（超平面因子の引き戻しが座標の生成するイデアルであること）を
高さの計算へ繋ぐには、次の 2 本が要る:

| 補題 | 内容 |
|---|---|
| `prod_infinitePlace_eq_absNorm` | ★**積公式**: `∏_{v 無限素点} v(x)^{mult} = N((x))` |
| `mulHeight_normalized` | ★★**正規化しても素朴高さは変わらない**: `H(x_j/x_{i₀}) = H(x)` |

★1 本目は「超平面の切断 `x_{i₀}` が定める主因子の次数が消える」ことの中身である
——`degFin((x_{i₀})) − degArch((x_{i₀})) = 0`、すなわち**積公式そのもの**。
★★2 本目は消費側（`projPointCoord`、`§9-C2b`）が**正規化座標**を出すことに対応する。

## ★★★これで高さの側の材料が揃った

`§9-853`（`H(x)·N(I) = ∏_arch (sup_i v(x_i))^mult`）と本ファイルの 2 本を合わせると、

    `log H(x) = [log N((x_{i₀})) − log N(I)] + Σ_v mult·log(sup_i v(x_i)/v(x_{i₀}))`

——★**右辺の第 1 項が `degFin`、第 2 項が `degArch`（Fubini–Study 計量）**である。
★★あとは「左の括弧が超平面因子の引き戻しの `degFin` である」ことを言えばよく、
それが段 C2c-1 の中身である。

## ★測定の記録

★★★ここでも在庫が効いた——`NumberField.InfinitePlace.prod_eq_abs_norm` と
`Ideal.absNorm_span_singleton` と `Algebra.coe_norm_int` を繋ぐだけであった。
★`|(n:ℝ)| = (n.natAbs : ℝ)` は `Int.cast_natAbs` では**なく**
`← Int.cast_abs` ＋ `Int.abs_eq_natAbs` で通る（2026-08-28 実測）。
-/

namespace ABC3.Found.GenEll

open NumberField Ideal Height

/-! ## ★★★★★★積公式の帰結 -/

/-- ★★★★★★**積公式の帰結** —— アルキメデス積は主イデアルのノルムに等しい。

    `∏_{v 無限素点} v(x)^{mult v} = N((x))`

★これは「主因子の算術次数が消える」ことの中身である
——`degFin` と `degArch` がちょうど打ち消し合う。
★★機構は `InfinitePlace.prod_eq_abs_norm` ＋ `absNorm_span_singleton` ＋ `Algebra.coe_norm_int`。 -/
theorem prod_infinitePlace_eq_absNorm (K : Type) [Field K] [NumberField K] (x : 𝓞 K) :
    ∏ w : InfinitePlace K, w ((x : K)) ^ w.mult
      = ((Ideal.span {x}).absNorm : ℝ) := by
  rw [NumberField.InfinitePlace.prod_eq_abs_norm, Ideal.absNorm_span_singleton,
    ← Algebra.coe_norm_int]
  push_cast
  rw [← Int.cast_abs, Int.abs_eq_natAbs]
  simp

/-! ## ★★正規化しても素朴高さは変わらない -/

/-- ★★**正規化しても素朴高さは変わらない**: `H(y_j/c) = H(y)`（`c ≠ 0`）。

★消費側（`projPointCoord`、`§9-C2b`）が出すのは正規化座標 `x_j/x_{i₀}` なので、
この形が要る。★★機構は `Height.mulHeight_smul_eq_mulHeight`（射影不変）だけである。 -/
theorem mulHeight_normalized (K : Type) [Field K] [NumberField K] {ι : Type} [Finite ι]
    (y : ι → K) (c : K) (hc : c ≠ 0) :
    Height.mulHeight (fun i => y i / c) = Height.mulHeight y := by
  have h : (fun i => y i / c) = c⁻¹ • y := by
    funext i
    simp [div_eq_inv_mul]
  rw [h, Height.mulHeight_smul_eq_mulHeight y (inv_ne_zero hc)]

/-! ## ★出典の紐付け(`.src`) -/

def prod_infinitePlace_eq_absNorm.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(積公式の帰結——アルキメデス積は主イデアルのノルム)",
    sectionId := "genell-prop-1-4" }

def mulHeight_normalized.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(正規化しても素朴高さは変わらない)",
    sectionId := "genell-prop-1-4" }

def prod_infinitePlace_eq_absNorm.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]" "NumberField.InfinitePlace.prod_eq_abs_norm(積公式)"
      (.inMathlib "NumberField.InfinitePlace.prod_eq_abs_norm") 2,
    .citation "[mathlib]" "Ideal.absNorm_span_singleton / Algebra.coe_norm_int"
      (.inMathlib "Ideal.absNorm_span_singleton") 2,
    .citation "[ABC3]" "mulHeight_mul_absNorm(段 C2c-2、§9-853)"
      (.inProject "ABC3" "ABC3.Found.GenEll.mulHeight_mul_absNorm") 2,
    .implicitStep
      ("★§9-853 と本ファイルの 2 本を合わせると " ++
       "log H(x) = [log N((x_{i₀})) − log N(I)] + Σ_v mult·log(sup_i v(x_i)/v(x_{i₀})) となり、" ++
       "**第 1 項が degFin、第 2 項が degArch(Fubini–Study 計量)** である。" ++
       "★★あとは「左の括弧が超平面因子の引き戻しの degFin である」ことを言えばよく、" ++
       "それが段 C2c-1 の中身である") 5,
    .implicitStep
      ("★|(n:ℝ)| = (n.natAbs : ℝ) は Int.cast_natAbs では**なく** " ++
       "← Int.cast_abs ＋ Int.abs_eq_natAbs で通る(2026-08-28 実測)") 2 ]

end ABC3.Found.GenEll
