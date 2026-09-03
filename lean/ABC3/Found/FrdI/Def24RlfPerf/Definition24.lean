/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Def24RlfCone
import ABC3.Found.FrdI.Def24
import ABC3.Found.FrdI.Thm42

/-!
# Def24RlfPerf —— `[FrdI] Definition 2.4` の分

☆もとの 1 枚を**条なしの項目ごと**に割ったものである（第 1458、案 a）。
★「1 ファイル = 1 ノード」を回復するための分割で、中身は動かしていない。
-/

namespace ABC3.Found.FrdI

open CategoryTheory
open scoped TensorProduct NNReal
universe v u w
variable {M : Type w} [AddCommMonoid M]

/-! ## ★1. `M^gp` の普遍性(群行き)

★`Found/Divisor/CartierMonoid.lean` の `gpLift` と同じ構成である。
★**そちらは `Found/Divisor` にあり本ファイル(`Found/FrdI`)からは読めない**ので、
ここに置き直す(名前を変えて衝突を避ける)。 -/

/-- ★★**`M^gp` の普遍性** —— 群 `G` への準同型は一意に持ち上がる。 -/
noncomputable def gpLiftTo {G : Type*} [AddCommGroup G] (f : M →+ G) : Gp M →+ G :=
  Algebra.GrothendieckAddGroup.lift f

@[simp] theorem gpLiftTo_toGp {G : Type*} [AddCommGroup G] (f : M →+ G) (m : M) :
    gpLiftTo f (toGp M m) = f m := by
  have hb : toGp M m = (AddLocalization.addMonoidOf (⊤ : AddSubmonoid M)) m := rfl
  rw [hb]
  simp only [gpLiftTo, Algebra.GrothendieckAddGroup.lift, Equiv.coe_fn_mk]
  rw [AddSubmonoid.LocalizationMap.lift_eq]

/-! ## ★2. 素点ごとの `ord_p` -/

/-- ★★**素点 `p` での `ord_p : M →+ ℝ`** —— `m ↦ (因子分解写像の `p` 成分)`。

★加法性は `Pf.of` の加法性と `factorAdd` の 2 本だけ。 -/
noncomputable def perfOrd {ι : Prime M → Pf M → ℝ≥0} (H : IsPerfFactorialWith M ι)
    (p : Prime M) : M →+ ℝ where
  toFun m := (factorMap ι (Pf.of m) p : ℝ)
  map_zero' := by
    show ((factorMap ι (Pf.of (0 : M)) p : ℝ≥0) : ℝ) = 0
    rw [map_zero, factorMap_zero H]
    rfl
  map_add' a b := by
    show ((factorMap ι (Pf.of (a + b)) p : ℝ≥0) : ℝ) = _
    rw [map_add, H.factorAdd]
    rfl

@[simp] theorem perfOrd_apply {ι : Prime M → Pf M → ℝ≥0} (H : IsPerfFactorialWith M ι)
    (p : Prime M) (m : M) : perfOrd H p m = (factorMap ι (Pf.of m) p : ℝ) := rfl

/-- ★`ord_p` は**非負**(`ℝ≥0` に値を取るから)。 -/
theorem perfOrd_nonneg {ι : Prime M → Pf M → ℝ≥0} (H : IsPerfFactorialWith M ι)
    (p : Prime M) (m : M) : 0 ≤ perfOrd H p m :=
  (factorMap ι (Pf.of m) p).coe_nonneg

/-! ## ★3. 族が `M` の点を分離すること -/

/-- ★★**`M` が sharp なら `Pf.of m = 0 ⟹ m = 0`**。

★`Pf.of m = 0` は `∃ k ≥ 1, k • m = 0` であり、`k ≥ 1` なので `m` は加法的単元。 -/
theorem eq_zero_of_pf_of_eq_zero (hsharp : IsSharp M) {m : M} (h : Pf.of m = 0) : m = 0 := by
  obtain ⟨k, e⟩ := Quotient.exact h
  have e' : ((k : ℕ) * 1) • m = ((k : ℕ) * 1) • (0 : M) := e
  rw [mul_one, smul_zero] at e'
  refine hsharp m (isAddUnit_iff_exists_neg.mpr ⟨((k : ℕ) - 1) • m, ?_⟩)
  rw [← succ_nsmul', show ((k : ℕ) - 1 + 1) = (k : ℕ) from Nat.sub_add_cancel k.2]
  exact e'

/-- ★★★**族 `{ord_p}` は `M` の点を分離する** —— `factorInj` そのもの。 -/
theorem perfOrd_separates {ι : Prime M → Pf M → ℝ≥0} (H : IsPerfFactorialWith M ι)
    (m : M) (h : ∀ p : Prime M, perfOrd H p m = 0) : m = 0 := by
  refine eq_zero_of_pf_of_eq_zero H.divisorial.2 (H.factorInj ?_)
  rw [factorMap_zero H]
  refine funext fun p => ?_
  have hp : ((factorMap ι (Pf.of m) p : ℝ≥0) : ℝ) = 0 := h p
  exact NNReal.coe_injective hp

/-! ## ★4. perf-factorial ⟹ 錐は sharp -/

/-- ★★★★★★**perf-factorial なら実化の錐は sharp**。

★★これが鎖 `rlf` の最後の穴(`rlf-flat`)である。
**有理性の議論も半環上の平坦性も要らなかった** —— 錐の生成元が 1 点しか
含まないので、素点ごとの `ord_p` を**束ねずに**使えばよい。 -/
theorem isSharp_rlfCone_of_perfFactorial (h : IsPerfFactorial M) : IsSharp (rlfCone M) := by
  obtain ⟨ι, H⟩ := h
  refine isSharp_rlfCone_of_family (fun p : Prime M => gpLiftTo (perfOrd H p)) ?_ ?_
  · intro p m
    rw [gpLiftTo_toGp]
    exact perfOrd_nonneg H p m
  · intro m hm
    refine perfOrd_separates H m fun p => ?_
    have := hm p
    rwa [gpLiftTo_toGp] at this

/-- ★locator —— perf-factorial なら実化の錐は sharp。 -/
def isSharp_rlfCone_of_perfFactorial.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 51,
    item := "Definition 2.4, (i) — realification(perf-factorial なら錐は sharp)",
    sectionId := "frdi-def-2-4" }

end ABC3.Found.FrdI
