import ABC3.Found.NCBelyi.NestedInduction
import ABC3.Found.NCBelyi.CoeffBound
import Mathlib.FieldTheory.Minpoly.IsIntegrallyClosed

/-!
# [NCBelyi] Lemma 2.4 —— 共役で閉じた集合(`Found`)

原典: S. Mochizuki, *Noncritical Belyi Maps* [NCBelyi]、物理 p.5。

原文 (NCBelyi p.5):
> Also, we may replace S by the union of all Gal(Q/Q)-conjugates of S and assume,

## ★★★★なぜ要るか —— ノルムの評価が要求している

`Lemma 2.4` は正規化で `|α| ≤ 1`(`∀ α ∈ S∖{∞}`)を達成する。
`CoeffBound.lean` の評価が要求するのは
**`f₀` の根がすべて単位閉円板の中にあること**である。

★`f₀` は `α₀` の最小多項式なので、その根は `α₀` の**共役**である。
★★したがって「`S` が共役で閉じている」ことが無ければ、
`|α| ≤ 1` は `f₀` の根の位置を何も言わない。
★★★**原文の『without loss of generality』の 1 行は、そこを支えている。**

## ★★★★★`Gal(ℚ̄/ℚ)` を持ち込まない

原文は `Gal(ℚ̄/ℚ)`-共役と書く。★ここでは**最小多項式の根**として定義する:

    IsConjStable S  ⟺  x ∈ S の最小多項式の(ℂ での)根はすべて S に入る

★★`ℂ` の体自己同型は選択公理でしか作れない(連続でないものが大量にある)ので、
**根の側で書くほうが素直**である。★★★有限次拡大の範囲では両者は一致する。
-/

namespace ABC3.Found.NCBelyi

open Polynomial IntermediateField

/-! ## ★共役の集まり -/

/-- `x` の(`ℂ` での)共役の集まり —— `minpoly ℚ x` の根。 -/
noncomputable def conjSet (x : ℂ) : Finset ℂ :=
  (((minpoly ℚ x).map (algebraMap ℚ ℂ)).roots).toFinset

theorem mem_conjSet_iff {x z : ℂ} (hx : IsIntegral ℚ x) :
    z ∈ conjSet x ↔ Polynomial.aeval z (minpoly ℚ x) = 0 := by
  have hne : ((minpoly ℚ x).map (algebraMap ℚ ℂ)) ≠ 0 :=
    (Polynomial.map_ne_zero_iff (algebraMap ℚ ℂ).injective).2 (minpoly.ne_zero hx)
  rw [conjSet, Multiset.mem_toFinset, Polynomial.mem_roots hne]
  simp [Polynomial.IsRoot, Polynomial.eval_map, ← Polynomial.aeval_def]

/-- ★`x` は自分の共役である。 -/
theorem self_mem_conjSet {x : ℂ} (hx : IsIntegral ℚ x) : x ∈ conjSet x :=
  (mem_conjSet_iff hx).2 (minpoly.aeval ℚ x)

/-- ★★**共役は同じ最小多項式を持つ** —— `minpoly` が既約かつモニックだから。 -/
theorem minpoly_eq_of_mem_conjSet {x z : ℂ} (hx : IsIntegral ℚ x) (hz : z ∈ conjSet x) :
    minpoly ℚ z = minpoly ℚ x :=
  (minpoly.eq_of_irreducible_of_monic (minpoly.irreducible hx)
    ((mem_conjSet_iff hx).1 hz) (minpoly.monic hx)).symm

theorem isIntegral_of_mem_conjSet {x z : ℂ} (hx : IsIntegral ℚ x) (hz : z ∈ conjSet x) :
    IsIntegral ℚ z :=
  isIntegral_of_aeval_eq_zero (minpoly.ne_zero hx) ((mem_conjSet_iff hx).1 hz)

/-- ★★★**共役は同じ reduced degree を持つ** —— `[ℚ(z):ℚ] = deg (minpoly ℚ z)` だから。 -/
theorem redDeg_eq_of_mem_conjSet {x z : ℂ} (hx : IsIntegral ℚ x) (hz : z ∈ conjSet x) :
    redDeg z = redDeg x := by
  have hzint : IsIntegral ℚ z := isIntegral_of_mem_conjSet hx hz
  rw [redDeg, redDeg, IntermediateField.adjoin.finrank hzint,
    IntermediateField.adjoin.finrank hx, minpoly_eq_of_mem_conjSet hx hz]

/-- ★★**共役の共役はもとの共役** —— 最小多項式が一致するから。 -/
theorem conjSet_eq_of_mem {x z : ℂ} (hx : IsIntegral ℚ x) (hz : z ∈ conjSet x) :
    conjSet z = conjSet x := by
  rw [conjSet, conjSet, minpoly_eq_of_mem_conjSet hx hz]

/-! ## ★★★★★★共役で閉じていること -/

/-- **共役で閉じた集合** —— 原文の `Gal(ℚ̄/ℚ)`-stable。

原文 (NCBelyi p.5):
> without loss of generality, that S is Gal(Q/Q)-stable. -/
def IsConjStable (S : Finset ℂ) : Prop := ∀ x ∈ S, conjSet x ⊆ S

/-- **共役閉包** —— 原文の「`S` を `Gal(ℚ̄/ℚ)`-共役の合併で置き換える」。 -/
noncomputable def conjClosure (S : Finset ℂ) : Finset ℂ := S.biUnion conjSet

theorem subset_conjClosure {S : Finset ℂ} (hint : ∀ x ∈ S, IsIntegral ℚ x) :
    S ⊆ conjClosure S := by
  intro x hx
  exact Finset.mem_biUnion.2 ⟨x, hx, self_mem_conjSet (hint x hx)⟩

theorem isIntegral_of_mem_conjClosure {S : Finset ℂ} (hint : ∀ x ∈ S, IsIntegral ℚ x)
    {z : ℂ} (hz : z ∈ conjClosure S) : IsIntegral ℚ z := by
  obtain ⟨x, hxS, hzx⟩ := Finset.mem_biUnion.1 hz
  exact isIntegral_of_mem_conjSet (hint x hxS) hzx

/-- ★★★★**共役閉包は共役で閉じている**。 -/
theorem isConjStable_conjClosure {S : Finset ℂ} (hint : ∀ x ∈ S, IsIntegral ℚ x) :
    IsConjStable (conjClosure S) := by
  intro z hz w hw
  obtain ⟨x, hxS, hzx⟩ := Finset.mem_biUnion.1 hz
  have hxint : IsIntegral ℚ x := hint x hxS
  rw [conjSet_eq_of_mem hxint hzx] at hw
  exact Finset.mem_biUnion.2 ⟨x, hxS, hw⟩

/-- ★★★★★**共役閉包は `m(S)` を変えない** —— 共役の reduced degree は等しいから。

★これが「without loss of generality」の中身である
——測度を変えないので帰納法に影響しない。 -/
theorem maxRedDeg_conjClosure {S : Finset ℂ} (hint : ∀ x ∈ S, IsIntegral ℚ x) :
    maxRedDeg (conjClosure S) = maxRedDeg S := by
  refine le_antisymm ?_ ?_
  · refine Finset.sup_le (fun z hz => ?_)
    obtain ⟨x, hxS, hzx⟩ := Finset.mem_biUnion.1 hz
    rw [redDeg_eq_of_mem_conjSet (hint x hxS) hzx]
    exact redDeg_le_maxRedDeg hxS
  · exact Finset.sup_mono (subset_conjClosure hint)

/-! ## ★★★★★★★ノルムの評価への橋 -/

/-- ★★★★★★★★**共役で閉じていれば `f₀` の根はすべて単位閉円板の中**。

原文 (NCBelyi p.5):
> Also, we may replace S by the union of all Gal(Q/Q)-conjugates of S and assume,

★★これが `CoeffBound.lean` の仮定 `∀ z ∈ f.roots, ‖z‖ ≤ 1` を供給する。
★★★**原文の『without loss of generality』はここへ効く**
——共役で閉じていなければ、`|α| ≤ 1`(`α ∈ S`)は `f₀` の根について何も言わない。 -/
theorem norm_roots_minpoly_le_one {S : Finset ℂ} (hstab : IsConjStable S)
    (hnorm : ∀ x ∈ S, ‖x‖ ≤ 1) {a : ℂ} (haS : a ∈ S) (hint : IsIntegral ℚ a) :
    ∀ z ∈ ((minpoly ℚ a).map (algebraMap ℚ ℂ)).roots, ‖z‖ ≤ 1 := by
  intro z hz
  have hzc : z ∈ conjSet a := by rw [conjSet, Multiset.mem_toFinset]; exact hz
  exact hnorm z (hstab a haS hzc)

/-! ## ★出典の紐付け(`.src`) -/

def IsConjStable.src : ABC3.Meta.Source :=
  { paper := "NCBelyi", pdfPage := 5,
    item := "Lemma 2.4(S が Gal(ℚ̄/ℚ)-stable であるとしてよいこと)",
    sectionId := "ncbelyi-lemma-2-4" }

def conjClosure.src : ABC3.Meta.Source :=
  { paper := "NCBelyi", pdfPage := 5,
    item := "Lemma 2.4(S を Gal(ℚ̄/ℚ)-共役の合併で置き換える)",
    sectionId := "ncbelyi-lemma-2-4" }

def norm_roots_minpoly_le_one.src : ABC3.Meta.Source :=
  { paper := "NCBelyi", pdfPage := 5,
    item := "Lemma 2.4(共役で閉じていれば f₀ の根はすべて |z| ≤ 1)",
    sectionId := "ncbelyi-lemma-2-4" }

end ABC3.Found.NCBelyi
