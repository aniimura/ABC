/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Meta.Claim
import Mathlib.NumberTheory.NumberField.ProductFormula

/-!
# 算術因子 —— `deg^arith` が `B(L)` の像で消えること(`Example 6.3`)

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.113–114。

原文 (FrdI p.113):
> as an effective arithmetic divisor on F, and to an element of the group

原文 (FrdI p.114):
> vanishes on the

## ★★何を閉じたか

依存グラフ(`ResearchPaper/frdi-decomposition.json`)の鎖 `arith` のうち
**`deg-vanish`** ——「`deg^arith` が `B(L)` の像で消える」—— を閉じた。
★原文は `Example 6.3` の末尾で「one verifies immediately that `deg^arith_L`
vanishes on the image of `B(L)`」と一言で畳むが、その中身は**積公式**である。

★★mathlib に `NumberField.prod_abs_eq_one` が在る:

  `(∏ w : InfinitePlace, w x ^ w.mult) * ∏ᶠ (w : FinitePlace), w x = 1`

両辺の**対数**を取れば `Σ_v log|x|_v = 0` になり、
`deg^arith(div(x)) = −Σ_v log|x|_v = 0` が出る。

## ★設計 —— 対数を先に取る

原文は `ord(F_v)` を `F_v^×/O_v^×` として定め、非アルキメデスなら `≅ ℤ`、
アルキメデスなら `≅ ℝ` と述べ、`deg^arith` を素点ごとに
「`−[F_v:ℝ]·log λ`」「`log #(O_v/(λ))`」で定める。

★★**本ファイルは `deg^arith` の値の側を直接扱う** —— すなわち各素点の成分を
`−log|x|_v`(無限素点は重複度 `mult` つき)と置く。
非アルキメデスでは `log|x|_v = −ord_v(x)·log(Nv)` なので、
これは原文の `ord_v(x)·log(Nv)` と**同じ数**である。

★★★**逸脱(記録)**: したがって本ファイルの `ArithPlace L →₀ ℝ` は
原文の `Φ(L)^gp = ⊕_v ord(F_v)`(有限素点成分が `ℤ`)を
**`deg^arith` で見た像**に相当する。`ord(F_v) ≅ ℤ` の同一視そのものは
鎖 `arith` の節点 `ord-mon` であり、まだ実装していない
(`Skeleton/Divisor/ArithDivisor.lean` に statement がある)。
-/

namespace ABC3.Found.Divisor

open NumberField
open scoped Classical

variable {L : Type*} [Field L] [NumberField L]

/-! ## ★1. 素点と `−log|x|_v` -/

/-- ★**`V(L)`** —— 有限素点と無限素点の合併。 -/
abbrev ArithPlace (L : Type*) [Field L] [NumberField L] : Type _ :=
  FinitePlace L ⊕ InfinitePlace L

/-- ★各素点での `−log|x|_v`(無限素点は重複度 `mult` つき)。

★非アルキメデスでは `= ord_v(x)·log(Nv)`、
アルキメデスでは `= −[F_v:ℝ]·log|x|_v` である。 -/
noncomputable def arithPlaceLog (x : L) : ArithPlace L → ℝ :=
  Sum.elim (fun w : FinitePlace L => -Real.log (w x))
    (fun w : InfinitePlace L => -(w.mult : ℝ) * Real.log (w x))

/-- ★`|x|_v ≠ 1` となる有限素点は有限個(mathlib)。 -/
noncomputable def arithFinSupport {x : L} (hx : x ≠ 0) : Finset (FinitePlace L) :=
  (FinitePlace.hasFiniteMulSupport hx).toFinset

theorem prod_arithFinSupport {x : L} (hx : x ≠ 0) :
    ∏ᶠ (w : FinitePlace L), w x = ∏ w ∈ arithFinSupport hx, w x :=
  finprod_eq_prod_of_mulSupport_toFinset_subset _ (FinitePlace.hasFiniteMulSupport hx)
    (Finset.Subset.refl _)

/-- ★★★**積公式の対数版** —— `Σ_v log|x|_v = 0`。 -/
theorem sum_log_place_eq_zero {x : L} (hx : x ≠ 0) :
    (∑ w : InfinitePlace L, (w.mult : ℝ) * Real.log (w x))
      + ∑ w ∈ arithFinSupport hx, Real.log (w x) = 0 := by
  have hpf := NumberField.prod_abs_eq_one hx
  rw [prod_arithFinSupport hx] at hpf
  have h1 : ∀ w : InfinitePlace L, w x ≠ 0 := fun w => (InfinitePlace.pos_iff.mpr hx).ne'
  have h2 : ∀ w : FinitePlace L, w x ≠ 0 := fun w => (FinitePlace.pos_iff.mpr hx).ne'
  have hA : (∏ w : InfinitePlace L, (w x) ^ w.mult) ≠ 0 :=
    Finset.prod_ne_zero_iff.mpr (fun w _ => pow_ne_zero _ (h1 w))
  have hB : (∏ w ∈ arithFinSupport hx, w x) ≠ 0 :=
    Finset.prod_ne_zero_iff.mpr (fun w _ => h2 w)
  have hlog := congrArg Real.log hpf
  rw [Real.log_mul hA hB, Real.log_one] at hlog
  rw [Real.log_prod (f := fun w : InfinitePlace L => (w x) ^ w.mult)
      (fun w _ => pow_ne_zero _ (h1 w)),
    Real.log_prod (f := fun w : FinitePlace L => w x) (fun w _ => h2 w)] at hlog
  simpa [Real.log_pow] using hlog

/-! ## ★2. 算術因子 `div(x)` -/

theorem finite_support_arithPlaceLog {x : L} (hx : x ≠ 0) :
    (Function.support (arithPlaceLog x)).Finite := by
  refine Set.Finite.subset
    ((Set.Finite.image Sum.inl (FinitePlace.hasFiniteMulSupport hx)).union
      (Set.finite_range (Sum.inr : InfinitePlace L → ArithPlace L))) ?_
  rintro (w | w) hw
  · refine Or.inl ⟨w, ?_, rfl⟩
    simp only [Function.mem_mulSupport]
    intro h
    exact hw (by simp [arithPlaceLog, h])
  · exact Or.inr ⟨w, rfl⟩

/-- ★★**`x ∈ L^×` の算術因子** —— 各素点の成分が `−log|x|_v`。 -/
noncomputable def arithDiv (x : Lˣ) : ArithPlace L →₀ ℝ :=
  Finsupp.ofSupportFinite (arithPlaceLog (x : L)) (finite_support_arithPlaceLog x.ne_zero)

@[simp] theorem arithDiv_apply (x : Lˣ) (v : ArithPlace L) :
    arithDiv x v = arithPlaceLog (x : L) v := rfl

theorem arithPlaceLog_mul {x y : L} (hx : x ≠ 0) (hy : y ≠ 0) (v : ArithPlace L) :
    arithPlaceLog (x * y) v = arithPlaceLog x v + arithPlaceLog y v := by
  rcases v with w | w
  · show -Real.log (w (x * y)) = -Real.log (w x) + -Real.log (w y)
    rw [map_mul, Real.log_mul (FinitePlace.pos_iff.mpr hx).ne' (FinitePlace.pos_iff.mpr hy).ne']
    ring
  · show -(w.mult : ℝ) * Real.log (w (x * y))
      = -(w.mult : ℝ) * Real.log (w x) + -(w.mult : ℝ) * Real.log (w y)
    rw [map_mul, Real.log_mul (InfinitePlace.pos_iff.mpr hx).ne'
      (InfinitePlace.pos_iff.mpr hy).ne']
    ring

/-- ★★**`B(L) → Φ(L)^gp` は準同型**。 -/
theorem arithDiv_mul (x y : Lˣ) : arithDiv (x * y) = arithDiv x + arithDiv y := by
  refine Finsupp.ext (fun v => ?_)
  show arithPlaceLog ((x : L) * (y : L)) v = arithPlaceLog (x : L) v + arithPlaceLog (y : L) v
  exact arithPlaceLog_mul x.ne_zero y.ne_zero v

/-- ★★`B(L) = L^× → Φ(L)^gp` を群準同型として。 -/
noncomputable def arithDivHom : Lˣ →* Multiplicative (ArithPlace L →₀ ℝ) where
  toFun x := Multiplicative.ofAdd (arithDiv x)
  map_one' := by
    refine congrArg Multiplicative.ofAdd (Finsupp.ext (fun v => ?_))
    show arithPlaceLog ((1 : Lˣ) : L) v = 0
    rcases v with w | w <;> simp [arithPlaceLog]
  map_mul' x y := congrArg Multiplicative.ofAdd (arithDiv_mul x y)

/-! ## ★3. 算術次数と、その `B(L)` 上での消滅 -/

/-- ★**算術次数** `deg^arith : Φ(L)^gp → ℝ` —— 成分の総和。 -/
noncomputable def arithDegree : (ArithPlace L →₀ ℝ) →+ ℝ :=
  Finsupp.liftAddHom (fun _ : ArithPlace L => (AddMonoidHom.id ℝ))

@[simp] theorem arithDegree_apply (d : ArithPlace L →₀ ℝ) :
    arithDegree d = d.sum (fun _ r => r) := rfl

/-- ★台を含む有限集合。 -/
noncomputable def arithBigSupport {x : L} (hx : x ≠ 0) : Finset (ArithPlace L) :=
  (arithFinSupport hx).image Sum.inl ∪ (Finset.univ : Finset (InfinitePlace L)).image Sum.inr

theorem support_subset_arithBigSupport (x : Lˣ) :
    (arithDiv x).support ⊆ arithBigSupport x.ne_zero := by
  intro v hv
  rw [Finsupp.mem_support_iff, arithDiv_apply] at hv
  rcases v with w | w
  · refine Finset.mem_union_left _ (Finset.mem_image.mpr ⟨w, ?_, rfl⟩)
    simp only [arithFinSupport, Set.Finite.mem_toFinset, Function.mem_mulSupport]
    intro h
    exact hv (by simp [arithPlaceLog, h])
  · exact Finset.mem_union_right _ (Finset.mem_image.mpr ⟨w, Finset.mem_univ w, rfl⟩)

/-- ★★★★★**[FrdI] Example 6.3** —— **`deg^arith` は `B(L)` の像で消える**。

原文 (FrdI p.114):
> image of B(L).

★★これは**積公式**そのものである(`NumberField.prod_abs_eq_one` の対数版)。 -/
theorem arithDegree_arithDiv (x : Lˣ) : arithDegree (arithDiv x) = 0 := by
  have hx : (x : L) ≠ 0 := x.ne_zero
  have hsum : arithDegree (arithDiv x)
      = ∑ v ∈ arithBigSupport hx, arithPlaceLog (x : L) v := by
    rw [arithDegree_apply]
    exact Finsupp.sum_of_support_subset _ (support_subset_arithBigSupport x) _ (fun _ _ => rfl)
  rw [hsum, arithBigSupport]
  have hdisj : Disjoint ((arithFinSupport hx).image (Sum.inl : FinitePlace L → ArithPlace L))
      ((Finset.univ : Finset (InfinitePlace L)).image Sum.inr) := by
    rw [Finset.disjoint_left]
    rintro v hv1 hv2
    obtain ⟨a, _, rfl⟩ := Finset.mem_image.mp hv1
    obtain ⟨b, _, hb⟩ := Finset.mem_image.mp hv2
    cases hb
  rw [Finset.sum_union hdisj,
    Finset.sum_image (fun a _ b _ h => Sum.inl_injective h),
    Finset.sum_image (fun a _ b _ h => Sum.inr_injective h)]
  have h := sum_log_place_eq_zero hx
  show (∑ w ∈ arithFinSupport hx, -Real.log (w (x : L)))
      + ∑ w : InfinitePlace L, -(w.mult : ℝ) * Real.log (w (x : L)) = 0
  have h1 : (∑ w ∈ arithFinSupport hx, -Real.log (w (x : L)))
      = -∑ w ∈ arithFinSupport hx, Real.log (w (x : L)) := by
    rw [← Finset.sum_neg_distrib]
  have h2 : ∑ w : InfinitePlace L, -(w.mult : ℝ) * Real.log (w (x : L))
      = -∑ w : InfinitePlace L, (w.mult : ℝ) * Real.log (w (x : L)) := by
    rw [← Finset.sum_neg_distrib]
    exact Finset.sum_congr rfl (fun w _ => by ring)
  rw [h1, h2]
  linarith [h]

/-! ### ★出典の紐付け -/

/-- ★locator —— `Example 6.3` の素点の集合 `V(F)`。 -/
def ArithPlace.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 112, item := "Example 6.3 — V(F)(素点の集合)",
    sectionId := "frdi-example-6-3" }

/-- ★locator —— `Example 6.3` の `B(F) = F^× → Φ(F)^gp`。 -/
def arithDivHom.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 113, item := "Example 6.3 — B(F) = F^× → Φ(F)^gp",
    sectionId := "frdi-example-6-3" }

/-- ★locator —— `Example 6.3` の `deg^arith_L`。 -/
def arithDegree.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 113, item := "Example 6.3 — deg^arith_L : Φ(L)^gp → ℝ",
    sectionId := "frdi-example-6-3" }

/-- ★★locator —— `Example 6.3` の「`deg^arith` は `B(L)` の像で消える」(=積公式)。 -/
def arithDegree_arithDiv.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 114, item := "Example 6.3 — deg^arith は B(L) の像で消える",
    sectionId := "frdi-example-6-3" }

end ABC3.Found.Divisor
