/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.NorthcottCoord
import ABC3.Found.GenEll.DenominatorBound
import ABC3.Meta.Claim

/-!
# ★★★★★★★★段 C2c-2 —— 素朴高さ ＝ アルキメデス積 ÷ イデアルのノルム（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.6。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

## ★★★★★★★★これは何か —— 段 C2c-2

`§9-852` で段 C2c 前半を 3 つに分けた。本ファイルはその**第 2 段**である:

> 整数座標 `x : ι → 𝓞_K`（全部 0 でない）について
>
>     `H(x) · N(span{x_i}) = ∏_{v 無限素点} (sup_i v(x_i))^{mult v}`

★★すなわち **素朴高さ ＝ アルキメデス積 ÷ 座標イデアルのノルム**。
★これは `degArith = degArch + degFin` と**同じ形**である
——`log` を取れば `log H(x) = (アルキメデス側) − log N(I)` で、
`degFin = log N(I)` が符号違いで対応する。

## ★★★機構 —— mathlib に**すでにあった**

    `NumberField.absNorm_mul_finprod_finitePlace_eq_one`
      : `N(span{x_i}) · ∏ᶠ_{v 有限素点} sup_i v(x_i) = 1`

★これは「素朴高さの有限素点部分は座標イデアルのノルムの逆数である」ことの
mathlib 側の言い方であり（`Mathlib/NumberTheory/Height/NumberField.lean:174`）、
★★`NumberField.mulHeight_eq`（`mulHeight` を無限素点の積と有限素点の積に分ける）と
合わせるだけで本段が出る。

★★★**在庫を引いて分かったこと**——段 C2c-2 は**書く必要がほとんど無かった**。
`ResearchPaper/mathlib-gap.json` に「在庫にある」と記録する。

## ★残っている段（明示）

★★段 C2c-1（超平面因子の引き戻しが座標の生成するイデアルであること）と
段 C2c-3（与えられた計量と Fubini–Study 型の差が一様に有界であること）が残る。
-/

namespace ABC3.Found.GenEll

open NumberField Ideal Height

/-! ## ★イデアルのノルムは 0 でない -/

/-- ★**座標が全部 0 でなければ、生成イデアルのノルムは 0 でない**。 -/
theorem absNorm_span_range_ne_zero (K : Type) [Field K] [NumberField K] {ι : Type} [Finite ι]
    (x : ι → 𝓞 K) (hx : x ≠ 0) :
    ((Ideal.span (Set.range x)).absNorm : ℝ) ≠ 0 := by
  have h := NumberField.absNorm_mul_finprod_finitePlace_eq_one (K := K) (x := x) hx
  intro h0
  rw [h0, zero_mul] at h
  exact zero_ne_one h

/-! ## ★★★有限素点側 —— `log N(I)` と非アルキメデス積 -/

/-- ★★★**`log N(I)` は非アルキメデス積の対数の符号違いである**。

★これが「有限部分の高さ `degFin = log #(𝓞_K/I)` が素朴高さの有限素点部分に対応する」
という段の中身である。 -/
theorem log_absNorm_eq_neg_log_finprod (K : Type) [Field K] [NumberField K]
    {ι : Type} [Finite ι] (x : ι → 𝓞 K) (hx : x ≠ 0) :
    Real.log ((Ideal.span (Set.range x)).absNorm : ℝ)
      = - Real.log (∏ᶠ v : NumberField.FinitePlace K, ⨆ i, v ((x i : K))) := by
  have h := NumberField.absNorm_mul_finprod_finitePlace_eq_one (K := K) (x := x) hx
  have hA := absNorm_span_range_ne_zero K x hx
  have hB : (∏ᶠ v : NumberField.FinitePlace K, ⨆ i, v ((x i : K))) ≠ 0 := by
    intro h0
    rw [h0, mul_zero] at h
    exact zero_ne_one h
  have h2 : Real.log (((Ideal.span (Set.range x)).absNorm : ℝ)
      * ∏ᶠ v : NumberField.FinitePlace K, ⨆ i, v ((x i : K))) = 0 := by rw [h]; simp
  rw [Real.log_mul hA hB] at h2
  linarith

/-! ## ★★★★★★★★段 C2c-2 の本体 -/

/-- ★★★★★★★★**素朴高さ ＝ アルキメデス積 ÷ 座標イデアルのノルム** —— 段 C2c-2。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

    `H(x) · N(span{x_i}) = ∏_{v 無限素点} (sup_i v(x_i))^{mult v}`

★★これは `degArith = degArch + degFin` と**同じ形**である。
★★★機構は mathlib の `absNorm_mul_finprod_finitePlace_eq_one` と `mulHeight_eq` だけである
——**在庫を引いたら書く必要がほとんど無かった**。 -/
theorem mulHeight_mul_absNorm (K : Type) [Field K] [NumberField K] {ι : Type} [Finite ι]
    (x : ι → 𝓞 K) (hx : x ≠ 0) :
    Height.mulHeight (fun i => (x i : K)) * ((Ideal.span (Set.range x)).absNorm : ℝ)
      = ∏ v : InfinitePlace K, (⨆ i, v ((x i : K))) ^ v.mult := by
  have hx0 : (fun i => (x i : K)) ≠ 0 := by
    obtain ⟨i, hi⟩ := Function.ne_iff.mp hx
    refine Function.ne_iff.mpr ⟨i, ?_⟩
    show ((x i : K)) ≠ (0 : ι → K) i
    simpa using (by exact_mod_cast hi : ((x i : K)) ≠ 0)
  have h := NumberField.absNorm_mul_finprod_finitePlace_eq_one (K := K) (x := x) hx
  rw [NumberField.mulHeight_eq hx0]
  calc (∏ v : InfinitePlace K, (⨆ i, v ((x i : K))) ^ v.mult)
        * (∏ᶠ v : NumberField.FinitePlace K, ⨆ i, v ((x i : K)))
        * ((Ideal.span (Set.range x)).absNorm : ℝ)
      = (∏ v : InfinitePlace K, (⨆ i, v ((x i : K))) ^ v.mult)
        * (((Ideal.span (Set.range x)).absNorm : ℝ)
          * ∏ᶠ v : NumberField.FinitePlace K, ⨆ i, v ((x i : K))) := by ring
    _ = _ := by rw [h, mul_one]

/-- ★★★★★**対数の形** —— `degArith = degArch + degFin` と直に見比べられる。 -/
theorem log_mulHeight_add_log_absNorm (K : Type) [Field K] [NumberField K]
    {ι : Type} [Finite ι] (x : ι → 𝓞 K) (hx : x ≠ 0) :
    Real.log (Height.mulHeight (fun i => (x i : K)))
        + Real.log (((Ideal.span (Set.range x)).absNorm : ℝ))
      = Real.log (∏ v : InfinitePlace K, (⨆ i, v ((x i : K))) ^ v.mult) := by
  have h1 : (0:ℝ) < Height.mulHeight (fun i => (x i : K)) :=
    lt_of_lt_of_le zero_lt_one (Height.one_le_mulHeight _)
  rw [← Real.log_mul h1.ne' (absNorm_span_range_ne_zero K x hx),
    mulHeight_mul_absNorm K x hx]

/-! ## ★出典の紐付け(`.src`) -/

def log_absNorm_eq_neg_log_finprod.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(log N(I) は非アルキメデス積の対数の符号違い)",
    sectionId := "genell-prop-1-4" }

def mulHeight_mul_absNorm.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(段 C2c-2——素朴高さ ＝ アルキメデス積 ÷ イデアルのノルム)",
    sectionId := "genell-prop-1-4" }

def log_mulHeight_add_log_absNorm.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(対数の形——degArith = degArch + degFin と見比べる)",
    sectionId := "genell-prop-1-4" }

def mulHeight_mul_absNorm.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]" "NumberField.absNorm_mul_finprod_finitePlace_eq_one"
      (.inMathlib "NumberField.absNorm_mul_finprod_finitePlace_eq_one") 2,
    .citation "[mathlib]" "NumberField.mulHeight_eq(無限素点の積と有限素点の積に分ける)"
      (.inMathlib "NumberField.mulHeight_eq") 2,
    .implicitStep
      ("★★★**在庫を引いて分かったこと**——段 C2c-2 は**書く必要がほとんど無かった**。" ++
       "mathlib の absNorm_mul_finprod_finitePlace_eq_one が" ++
       "「素朴高さの有限素点部分は座標イデアルのノルムの逆数である」を既に持っていた" ++
       "(Mathlib/NumberTheory/Height/NumberField.lean:174、2026-08-28 実測)") 3,
    .implicitStep
      ("★★段 C2c-1(超平面因子の引き戻しが座標の生成するイデアルであること)と " ++
       "段 C2c-3(与えられた計量と Fubini–Study 型の差が一様に有界であること)が残る") 5 ]

end ABC3.Found.GenEll
