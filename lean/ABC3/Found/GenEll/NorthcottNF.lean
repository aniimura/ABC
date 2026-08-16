import ABC3.Meta.Claim
import Mathlib.NumberTheory.Height.NumberField
import Mathlib.NumberTheory.Height.Northcott
import Mathlib.NumberTheory.NumberField.InfinitePlace.Embeddings

/-!
# [GenEll] Proposition 1.4, (iv) の基底 —— 数体上の Northcott 性(整元の場合)(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.6。

原文 (GenEll p.6):
> (iv) Let d be a positive integer, C ∈ R. Suppose further that the line bundle LQ is ample on XQ. Then the set of points x ∈ X(Q)≤d [cf. Example 1.3, (i)] such that htL(x) ≤ C is ﬁnite.

## ★★mathlib の状況(2026-08-17 実測)

`Mathlib/NumberTheory/Height/Northcott.lean` は

> A field that satisfies the Northcott property for `mulHeight₁` also does for `logHeight₁`.

という**条件つきの instance を 1 つ持つだけ**で、
★**`Northcott (mulHeight₁ (K := K))` の基底 instance を 1 つも持たない**。

`Found/GenEll/NorthcottRat.lean` で `K = ℚ` の場合は与えた。
★本ファイルは**任意の数体 `K` の整元**について与える。

## ★★何が取れて、何が取れていないか

- ★**取れた**: `{x : K | x は代数的整数 ∧ H(x) ≤ B}` は有限
- ★★**取れていない**: `{x : K | H(x) ≤ B}`(整でない元を含む全体)

★残っているのは**分母の評価**である——`H(x) ≤ B` から
「分母イデアルのノルムが `B` 以下」を出し、`Ideal.absNorm I ∈ I` で
有界な自然数 `d` を取って `d·x` を整にする、という段。
そこには**分母イデアルの構成**が要り、本ファイルの外である。

★**半分取れたことを全部取れたと読まない。** `.src` は条つきである。

## ★機構

**`v(x) ≤ H(x)`(`v` は無限素点)** ——
高さの定義に現れる因子が**すべて 1 以上**なので、1 つの因子だけ残せばよい。
★これは整性を一切使わない。

あとは mathlib の **`NumberField.Embeddings.finite_of_norm_le`**
(「共役がすべて有界な代数的整数は有限個」)にそのまま渡す。
-/

namespace ABC3.Found.GenEll

open NumberField

variable {K : Type*} [Field K] [NumberField K]

/-- ★**無限素点での絶対値は高さで抑えられる**: `v(x) ≤ H(x)`。

★高さ `H(x) = (∏_{v|∞} max(v x,1)^{mult}) · (∏ᶠ_{v∤∞} max(v x,1))` の
因子はすべて 1 以上だから、1 つだけ残せる。
★★**整性は一切使わない。** -/
theorem infinitePlace_le_mulHeight₁ (v : InfinitePlace K) (x : K) :
    v x ≤ Height.mulHeight₁ x := by
  classical
  rw [NumberField.mulHeight₁_eq]
  -- 有限素点の積は 1 以上
  have hP : (1 : ℝ) ≤ ∏ᶠ w : FinitePlace K, max (w x) 1 :=
    one_le_finprod (fun w => le_max_right _ _)
  -- 無限素点の積は自分の因子以上
  have hone : ∀ w : InfinitePlace K, (1 : ℝ) ≤ max (w x) 1 ^ w.mult := by
    intro w
    exact one_le_pow₀ (le_max_right _ _)
  have hA : max (v x) 1 ^ v.mult
      ≤ ∏ w : InfinitePlace K, max (w x) 1 ^ w.mult := by
    rw [← Finset.mul_prod_erase _ _ (Finset.mem_univ v)]
    nth_rewrite 1 [← mul_one (max (v x) 1 ^ v.mult)]
    exact mul_le_mul_of_nonneg_left
      (Finset.one_le_prod (fun w _ => hone w)) (le_trans zero_le_one (hone v))
  have hApos : (0 : ℝ) ≤ ∏ w : InfinitePlace K, max (w x) 1 ^ w.mult :=
    le_trans zero_le_one (le_trans (hone v) hA)
  have hv1 : v x ≤ max (v x) 1 := le_max_left _ _
  have hpow : max (v x) 1 ≤ max (v x) 1 ^ v.mult := by
    refine le_self_pow₀ (le_max_right _ _) ?_
    exact InfinitePlace.mult_ne_zero
  calc v x ≤ max (v x) 1 := hv1
    _ ≤ max (v x) 1 ^ v.mult := hpow
    _ ≤ ∏ w : InfinitePlace K, max (w x) 1 ^ w.mult := hA
    _ = (∏ w : InfinitePlace K, max (w x) 1 ^ w.mult) * 1 := (mul_one _).symm
    _ ≤ (∏ w : InfinitePlace K, max (w x) 1 ^ w.mult)
          * ∏ᶠ w : FinitePlace K, max (w x) 1 := by
        exact mul_le_mul_of_nonneg_left hP hApos

/-- ★★**数体の整元についての Northcott 性**。

> `{x : K | x は代数的整数 ∧ H(x) ≤ B}` は**有限**

★mathlib の `NumberField.Embeddings.finite_of_norm_le`
(「共役がすべて有界な代数的整数は有限個」)に上の評価を渡すだけ。

★★**整でない元まで含めた `{x : K | H(x) ≤ B}` は本ファイルでは取れていない。**
残るのは分母の評価である(上の docstring)。 -/
theorem finite_isIntegral_mulHeight₁_le (B : ℝ) :
    {x : K | IsIntegral ℤ x ∧ Height.mulHeight₁ x ≤ B}.Finite := by
  refine Set.Finite.subset (NumberField.Embeddings.finite_of_norm_le K ℂ B) ?_
  rintro x ⟨hint, hB⟩
  refine ⟨hint, fun φ => ?_⟩
  calc ‖φ x‖ = (InfinitePlace.mk φ) x := (InfinitePlace.apply φ x).symm
    _ ≤ Height.mulHeight₁ x := infinitePlace_le_mulHeight₁ _ x
    _ ≤ B := hB

/-! ## ★出典の紐付け(`.src`) -/

def finite_isIntegral_mulHeight₁_le.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 6,
    item := "Proposition 1.4, (iv)(数体の整元の場合のみ)",
    sectionId := "genell-prop-1-4" }

end ABC3.Found.GenEll
