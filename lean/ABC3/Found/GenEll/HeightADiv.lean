import ABC3.Found.GenEll.Conductor
import Mathlib.NumberTheory.Height.NumberField

/-!
# [GenEll] §1 —— **高さは算術因子の次数である**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.4。

原文 (GenEll p.4):
> [cf. [Szp], §1.1] determines a homomorphism APic(Spec(OF )) → R, which we shall

## ★★何を取るか —— mathlib の高さと本プロジェクトの Arakelov 層を繋ぐ

mathlib は数体上の高さ `Height.logHeight₁` を持つ:

> `logHeight₁(x) = Σ_{v|∞} mult_v·log⁺|x|_v + Σᶠ_{v∤∞} log⁺|x|_v`

本プロジェクトは算術因子 `ADiv(F)` と次数 `deg_F` を持つ。★その 2 つを繋ぐ:

> **`logHeight₁(x) = deg_F(polarADiv x)`**

ここで `polarADiv x` は `x` の**極因子**——有限素点では `max(−ord_v(x), 0)`、
無限素点では `mult_v·log⁺|x|_v` を係数とする有効算術因子である。

## ★★これが `Definition 1.1`/`1.2` に対して持つ意味

原文の高さ `ht_L̄(x)` は「点 `x` に沿って**算術直線束を引き戻し、その次数を取る**」
ものである。★本ファイルはその**最も単純な場合**——
`X = ℙ¹`、`L = O(1)`、`x` を有理点とする場合——を、
**`X^arc` を一切使わずに**取ったものにあたる。

★★**ただし一般の `X` には届かない。** 一般には `L` の計量が `X^arc` 上の
エルミート計量として与えられねばならず、そこは依然として空いている
(`blocked-leaves.json` の `[GenEll] Definition 1.1`)。
**単純な場合が取れたことを一般の場合が取れたと読まない。**

## ★機構

`ProductFormula.lean` の **`w(x) = q_v^{−ord_v(x)}`** をそのまま使う。

- `log⁺|x|_w = max(0, −ord_v(x)·log q_v) = max(−ord_v(x), 0)·log q_v`
  (`log q_v > 0` なので `max` の外へ出せる)
- 無限素点側は定義がそのまま一致する

★**`polarADiv` は有効である**——極因子だから当然だが、
`deg` が `logHeight₁` と一致するには符号が合っていなければならず、
`ordv` の符号は 2026-08-17 に一度直している。★ここが**その符号の独立な検算**になる。
-/

namespace ABC3.Found.GenEll

open NumberField IsDedekindDomain

variable {F : Type*} [Field F] [NumberField F]

/-! ## ★極因子 -/

/-- ★**`x` の極因子** —— 有限素点では `max(−ord_v(x), 0)`、
無限素点では `mult_v·log⁺|x|_v`。

★台が有限なのは `ordv_finite_support` から。 -/
noncomputable def polarADiv (x : Fˣ) : ADiv F :=
  (Finsupp.onFinset (ordv_finite_support x).toFinset
     (fun v => max (-(ordv v x)) 0)
     (by
       intro v hv
       simp only [Set.Finite.mem_toFinset, Set.mem_setOf_eq]
       intro hcon
       exact hv (by rw [hcon]; simp)),
   Finsupp.onFinset Finset.univ
     (fun v : InfinitePlace F => (v.mult : ℝ) * Real.posLog (v (x : F)))
     (fun v _ => Finset.mem_univ v))

/-- ★極因子は**有効**である。 -/
theorem polarADiv_isEffective (x : Fˣ) : (polarADiv x).IsEffective := by
  constructor
  · intro v
    exact le_max_right _ _
  · intro v
    exact mul_nonneg (Nat.cast_nonneg _) (Real.posLog_nonneg)

/-! ## ★有限側 -/

/-- ★**`log⁺|x|_w = max(−ord_v(x), 0)·log q_v`**。

★`log q_v > 0` なので `max` を掛け算の外へ出せる。 -/
theorem posLog_finitePlace_apply (w : NumberField.FinitePlace F) (x : Fˣ) :
    Real.posLog (w ((x : F)))
      = ((max (-(ordv (NumberField.FinitePlace.maximalIdeal w) x)) 0 : ℤ) : ℝ)
          * Real.log (residueCard (NumberField.FinitePlace.maximalIdeal w)) := by
  set v := NumberField.FinitePlace.maximalIdeal w with hv
  have hq : (0 : ℝ) ≤ Real.log (residueCard v) := (log_residueCard_pos v).le
  rw [Real.posLog, log_finitePlace_apply w x]
  push_cast
  rw [max_comm, max_mul_of_nonneg _ _ hq, zero_mul]

/-- ★有限素点の和が、極因子の有限側の次数に一致する。 -/
theorem finsum_posLog_finitePlace (x : Fˣ) :
    ∑ᶠ w : NumberField.FinitePlace F, Real.posLog (w ((x : F)))
      = (polarADiv x).fin.sum (fun v n => (n : ℝ) * Real.log (residueCard v)) := by
  classical
  rw [finsum_congr (fun w => posLog_finitePlace_apply w x)]
  -- 添字を `HeightOneSpectrum` に付け替える
  have h2 : ∑ᶠ w : NumberField.FinitePlace F,
      ((max (-(ordv (NumberField.FinitePlace.maximalIdeal w) x)) 0 : ℤ) : ℝ)
        * Real.log (residueCard (NumberField.FinitePlace.maximalIdeal w))
      = ∑ᶠ v : FinitePlace F,
          ((max (-(ordv v x)) 0 : ℤ) : ℝ) * Real.log (residueCard v) :=
    finsum_comp_equiv (M := ℝ)
      (e := NumberField.FinitePlace.equivHeightOneSpectrum (K := F))
      (f := fun v : FinitePlace F =>
        ((max (-(ordv v x)) 0 : ℤ) : ℝ) * Real.log (residueCard v))
  rw [h2]
  simp only [polarADiv, ADiv.fin]
  rw [Finsupp.onFinset_sum _ (fun a => by simp)]
  refine finsum_eq_finsetSum_of_support_subset _ ?_
  intro v hv
  rw [Function.mem_support] at hv
  have h0 : ordv v x ≠ 0 := by
    intro hcon
    apply hv
    rw [hcon]
    simp
  simpa using h0

/-! ## ★★高さ = 次数 -/

/-- ★★★**`logHeight₁(x) = deg_F(polarADiv x)`**。

原文 (GenEll p.4):
> [cf. [Szp], §1.1] determines a homomorphism APic(Spec(OF )) → R, which we shall

★mathlib の高さ(`Height.logHeight₁`)と、本プロジェクトの算術因子の次数
(`deg`)が**同じ数**であることを示す。

★★これは「高さ = 引き戻した算術直線束の次数」という原文の見方の、
**`X^arc` を要しない最小の場合**である。一般の `X` には届かない。 -/
theorem logHeight₁_eq_deg_polarADiv (x : Fˣ) :
    Height.logHeight₁ ((x : F)) = deg (polarADiv x) := by
  classical
  rw [NumberField.logHeight₁_eq, deg, add_comm]
  congr 1
  · exact finsum_posLog_finitePlace x
  · simp only [polarADiv, ADiv.arc]
    rw [Finsupp.onFinset_sum _ (fun _ => rfl)]

/-- ★**高さは非負である** —— 上の等式と極因子の有効性からの帰結。

★mathlib にも同種の補題はあるが、ここでは**我々の `deg` を経由して**出す。
これが符号の独立な検算になる。 -/
theorem logHeight₁_nonneg (x : Fˣ) : 0 ≤ Height.logHeight₁ ((x : F)) := by
  classical
  rw [logHeight₁_eq_deg_polarADiv x, deg]
  have hfin : (0 : ℝ) ≤ (polarADiv x).fin.sum
      (fun v n => (n : ℝ) * Real.log (residueCard v)) := by
    refine Finset.sum_nonneg fun v _ => ?_
    exact mul_nonneg (by exact_mod_cast (polarADiv_isEffective x).1 v)
      (log_residueCard_pos v).le
  have harc : (0 : ℝ) ≤ (polarADiv x).arc.sum (fun _ r => r) :=
    Finset.sum_nonneg fun v _ => (polarADiv_isEffective x).2 v
  linarith

/-! ## ★出典の紐付け(`.src`) -/

def logHeight₁_eq_deg_polarADiv.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 4, item := "§1 地の文(次数写像 deg_F)",
    sectionId := "genell-deg" }

end ABC3.Found.GenEll
