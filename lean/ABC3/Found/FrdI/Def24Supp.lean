/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.Def24

/-!
# [FrdI] `Definition 2.4, (i)` —— **台(`Supp`)の基本性質**

原文 (FrdI p.48):
> Now suppose that M is perf-factorial. Then we shall refer to the [subset which is

★`Def24.lean` は `perf-factorial` の**定義**(`IsPerfFactorialWith` の 11 フィールド)を
置いたが、**台の理論**は無かった。★★このファイルがそれを作る。

## ★★なぜ要るか —— `Proposition 4.1, (ii)–(v)` が全部これで動く

原文 (FrdI p.76):
> Now the necessity of this condition follows immediately from the structure of the

★原文が「the structure of the `Φ(A)_p`」「primary factorizations」と一言で書く所が、
**素点ごとの分解**である。★具体的には

  ★★★**`Supp` の任意の部分集合 `S` について `x = x_S + x_{Sᶜ}` と分解できる**
  (`exists_split_supp`)

が要る。★これが `Definition 2.4, (i), (d)` の条件のはたらきである ——
(d) は「台が `M^pf` の元の台に収まる `M^rlf_factor` の元は `M^pf` に戻る」と言うので、
**`x` の因子を `S` の外で `0` に落としたもの**が `M^pf` に戻る。

## ★道具

| 宣言 | 内容 |
|---|---|
| `iota_zero` | `ι p 0 = 0`(`p` 成分での加法性から) |
| `exists_restrict` | ★**`S` への制限は `M^pf` に戻る**((d) がここで効く) |
| `exists_split_supp` | ★★**`x = x_S + x_{Sᶜ}`** |
| `supp_indicator` | `Supp(1_S · f) = S ∩ Supp f` |
| `supp_mono_of_mle` | `x ≼ y ⟹ Supp x ⊆ Supp y` |
| `supp_factorMap_add` | `Supp(x+y) = Supp x ∪ Supp y` |
| `eq_zero_of_supp_empty` | 台が空なら `0` |
-/

namespace ABC3.Found.FrdI

open CategoryTheory
open scoped NNReal

universe w

variable {M : Type w} [AddCommMonoid M] {ι : Prime M → Pf M → ℝ≥0}

/-! ## ★1. `ι p 0 = 0` -/

/-- ★`ι p` は `0` を `0` に送る(`p` 成分での加法性から)。 -/
theorem iota_zero (H : IsPerfFactorialWith M ι) (p : Prime M) : ι p 0 = 0 := by
  have h := H.embedAdd p (zero_mem_pCarrierPf p) (zero_mem_pCarrierPf p)
  rw [add_zero] at h
  exact left_eq_add.mp h

/-! ## ★2. ★★★素点の集合への制限 —— `Definition 2.4, (i), (d)` のはたらき -/

/-- ★★**素点の集合 `S` への制限は `M^pf` に戻る**((d) の条件がここで効く)。 -/
theorem exists_restrict (H : IsPerfFactorialWith M ι) (x : Pf M) (S : Set (Prime M)) :
    ∃ y : Pf M, factorMap ι y = Set.indicator S (factorMap ι x) := by
  classical
  have hmem : ∀ p, Set.indicator S (factorMap ι x) p ∈ ι p '' pCarrierPf M p := by
    intro p
    by_cases hp : p ∈ S
    · rw [Set.indicator_of_mem hp]
      exact H.factorMem x p
    · rw [Set.indicator_of_notMem hp]
      exact ⟨0, zero_mem_pCarrierPf p, iota_zero H p⟩
  have hsub : Supp (Set.indicator S (factorMap ι x)) ⊆ Supp (factorMap ι x) := by
    intro p hp
    by_cases h : p ∈ S
    · rw [Supp, Set.mem_setOf_eq, Set.indicator_of_mem h] at hp
      exact hp
    · rw [Supp, Set.mem_setOf_eq, Set.indicator_of_notMem h] at hp
      exact absurd rfl hp
  obtain ⟨y, hy⟩ := H.supp _ hmem x hsub
  exact ⟨y, hy⟩

/-- ★指示関数の台。 -/
theorem supp_indicator (S : Set (Prime M)) (f : RlfFactor M) :
    Supp (Set.indicator S f) = S ∩ Supp f := by
  classical
  ext p
  by_cases h : p ∈ S
  · simp [Supp, h]
  · simp [Supp, h]

/-- ★★★**素点の集合 `S` で `x` を 2 つに割る**。

★★原文の「primary factorizations」の中身である。 -/
theorem exists_split_supp (H : IsPerfFactorialWith M ι) (x : Pf M) (S : Set (Prime M)) :
    ∃ y z : Pf M, x = y + z ∧ factorMap ι y = Set.indicator S (factorMap ι x)
      ∧ factorMap ι z = Set.indicator Sᶜ (factorMap ι x) := by
  classical
  obtain ⟨y, hy⟩ := exists_restrict H x S
  obtain ⟨z, hz⟩ := exists_restrict H x Sᶜ
  refine ⟨y, z, H.factorInj ?_, hy, hz⟩
  rw [H.factorAdd, hy, hz]
  funext p
  by_cases h : p ∈ S
  · rw [Pi.add_apply, Set.indicator_of_mem h, Set.indicator_of_notMem (by simpa using h),
      add_zero]
  · rw [Pi.add_apply, Set.indicator_of_notMem h, Set.indicator_of_mem (by simpa using h),
      zero_add]

/-! ## ★3. 台の単調性と加法性 -/

/-- ★`MLe` は台の包含を与える。 -/
theorem supp_mono_of_mle (H : IsPerfFactorialWith M ι) {x y : Pf M} (h : MLe x y) :
    Supp (factorMap ι x) ⊆ Supp (factorMap ι y) := by
  obtain ⟨c, rfl⟩ := h
  intro p hp
  rw [Supp, Set.mem_setOf_eq, H.factorAdd, Pi.add_apply]
  intro h0
  exact hp (by simpa using (add_eq_zero.mp h0).1)

/-- ★和の台は台の和。 -/
theorem supp_factorMap_add (H : IsPerfFactorialWith M ι) (x y : Pf M) :
    Supp (factorMap ι (x + y)) = Supp (factorMap ι x) ∪ Supp (factorMap ι y) := by
  ext p
  simp only [Supp, Set.mem_union, Set.mem_setOf_eq, H.factorAdd, Pi.add_apply]
  constructor
  · intro h
    by_contra hc
    rw [not_or, not_not, not_not] at hc
    exact h (by rw [hc.1, hc.2, add_zero])
  · rintro (h | h) h0
    · exact h (add_eq_zero.mp h0).1
    · exact h (add_eq_zero.mp h0).2

/-- ★台が空なら `0`。 -/
theorem eq_zero_of_supp_empty (H : IsPerfFactorialWith M ι) {x : Pf M}
    (h : Supp (factorMap ι x) = ∅) : x = 0 := by
  refine H.factorInj ?_
  rw [factorMap_zero H]
  funext p
  by_contra hp
  exact absurd (Set.eq_empty_iff_forall_notMem.mp h p) (by simpa [Supp] using hp)

/-- ★★**台が交わらない 2 元に共通の非零下界は無い**。

★★これが `Proposition 4.1, (ii)(iii)` の「disjoint supports」の条の骨である。 -/
theorem eq_zero_of_mle_of_supp_disjoint (H : IsPerfFactorialWith M ι) {x y d : Pf M}
    (hdx : MLe d x) (hdy : MLe d y)
    (hdisj : Supp (factorMap ι x) ∩ Supp (factorMap ι y) = ∅) : d = 0 := by
  refine eq_zero_of_supp_empty H (Set.eq_empty_iff_forall_notMem.mpr fun p hp => ?_)
  exact absurd (Set.eq_empty_iff_forall_notMem.mp hdisj p)
    (by simp only [not_not]; exact ⟨supp_mono_of_mle H hdx hp, supp_mono_of_mle H hdy hp⟩)

/-! ## ★出典の紐付け(条つき) -/

def exists_split_supp.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 48, item := "Definition 2.4, (i), (d) — 台による分解",
    sectionId := "frdi-def-2-4" }

end ABC3.Found.FrdI
