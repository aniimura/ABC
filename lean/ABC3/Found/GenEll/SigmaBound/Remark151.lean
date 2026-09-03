/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.Conductor
import Mathlib.NumberTheory.RamificationInertia.Basic
import ABC3.Found.GenEll.SigmaBound.Proposition17

/-!
# SigmaBound —— `[GenEll] Remark 1.5.1` の分

☆もとの 1 枚を**条なしの項目ごと**に割ったものである（第 1458、案 a）。
★「1 ファイル = 1 ノード」を回復するための分割で、中身は動かしていない。
-/

namespace ABC3.Found.GenEll

open scoped BigOperators
variable {F : Type*} [Field F] [NumberField F]

/-! ## ★★★★★★★★★★★2 つの因子の差 —— BD-class の中身 -/

/-- ★★**台が `Σ` の上にある集まりの `log N v` の総和の限界**(上の証明の中核を切り出したもの)。 -/
theorem sum_log_residueCard_le_of_over (V : Finset (FinitePlace F)) (Sig : Finset ℕ)
    (hprime : ∀ q ∈ Sig, q.Prime) (ch : FinitePlace F → ℕ)
    (hmem : ∀ v ∈ V, ch v ∈ Sig)
    (hover : ∀ v ∈ V, (v.asIdeal).LiesOver (Ideal.span {((ch v : ℕ) : ℤ)})) :
    ∑ v ∈ V, Real.log (residueCard v)
      ≤ (Module.finrank ℚ F : ℝ) * ∑ q ∈ Sig, Real.log q := by
  classical
  have hfib : ∑ q ∈ Sig, ∑ v ∈ V with ch v = q, Real.log (residueCard v)
      = ∑ v ∈ V, Real.log (residueCard v) :=
    Finset.sum_fiberwise_of_maps_to hmem _
  have hq : ∀ q ∈ Sig, ∑ v ∈ V with ch v = q,
      Real.log (residueCard v) ≤ (Module.finrank ℚ F : ℝ) * Real.log q := by
    intro q hqS
    refine sum_log_residueCard_le q (hprime q hqS) _ (fun v hv => ?_)
    obtain ⟨hv1, hv2⟩ := Finset.mem_filter.1 hv
    have hh := hover v hv1
    rwa [hv2] at hh
  calc ∑ v ∈ V, Real.log (residueCard v)
      = ∑ q ∈ Sig, ∑ v ∈ V with ch v = q, Real.log (residueCard v) := hfib.symm
    _ ≤ ∑ q ∈ Sig, (Module.finrank ℚ F : ℝ) * Real.log q := Finset.sum_le_sum hq
    _ = (Module.finrank ℚ F : ℝ) * ∑ q ∈ Sig, Real.log q := by rw [Finset.mul_sum]

/-- ★★★★★★★★★★★**`Σ` の外で一致する 2 つの因子の次数の差は `Σ` だけで決まる定数以下**。

原文 (GenEll p.9):
> immediately that the BD-class of log-condD on UX(Q) depends only on the pair

★★これが `Remark 1.5.1` の『BD-class は `(X_ℚ, D_ℚ)` だけに依る』の**計算部分**である
——`Σ` の外で導手が一致すれば、正規化次数の差は `Σ_{q∈Σ} log q` を超えない。
★★★**右辺は点にも定義体にも依らない**。だから差が有界(= BD-同値)になる。 -/
theorem degNormalized_sub_le_sum_log (c c' : ADiv F)
    (harc : c.arc = 0) (harc' : c'.arc = 0)
    (hc1 : ∀ v, c.fin v ≤ 1) (hc'0 : ∀ v, 0 ≤ c'.fin v)
    (Sig : Finset ℕ) (hprime : ∀ q ∈ Sig, q.Prime)
    (ch : FinitePlace F → ℕ)
    (hagree : ∀ v, ch v ∉ Sig → c.fin v = c'.fin v)
    (hover : ∀ v : FinitePlace F, (v.asIdeal).LiesOver (Ideal.span {((ch v : ℕ) : ℤ)})) :
    degNormalized c - degNormalized c' ≤ ∑ q ∈ Sig, Real.log q := by
  classical
  have hn : 0 < (Module.finrank ℚ F : ℝ) := by
    exact_mod_cast Module.finrank_pos (R := ℚ) (M := F)
  set T : Finset (FinitePlace F) := c.fin.support ∪ c'.fin.support with hT
  have hgz : ∀ (v : FinitePlace F), ((0 : ℤ) : ℝ) * Real.log (residueCard v) = 0 := by
    intro v; simp
  have hdegc : deg c = ∑ v ∈ T, ((c.fin v : ℤ) : ℝ) * Real.log (residueCard v) := by
    rw [deg, harc]
    simp only [Finsupp.sum_zero_index, add_zero]
    exact Finsupp.sum_of_support_subset _ Finset.subset_union_left _ (fun v _ => hgz v)
  have hdegc' : deg c' = ∑ v ∈ T, ((c'.fin v : ℤ) : ℝ) * Real.log (residueCard v) := by
    rw [deg, harc']
    simp only [Finsupp.sum_zero_index, add_zero]
    exact Finsupp.sum_of_support_subset _ Finset.subset_union_right _ (fun v _ => hgz v)
  have hdiff : deg c - deg c'
      = ∑ v ∈ T, (((c.fin v : ℤ) : ℝ) - ((c'.fin v : ℤ) : ℝ)) * Real.log (residueCard v) := by
    rw [hdegc, hdegc', ← Finset.sum_sub_distrib]
    exact Finset.sum_congr rfl (fun v _ => by ring)
  -- `Σ` の外の項は消える
  have hrestrict : ∑ v ∈ T,
        (((c.fin v : ℤ) : ℝ) - ((c'.fin v : ℤ) : ℝ)) * Real.log (residueCard v)
      = ∑ v ∈ T with ch v ∈ Sig,
        (((c.fin v : ℤ) : ℝ) - ((c'.fin v : ℤ) : ℝ)) * Real.log (residueCard v) := by
    refine (Finset.sum_subset (Finset.filter_subset _ _) (fun v hv hnot => ?_)).symm
    have hns : ch v ∉ Sig := by
      intro hc
      exact hnot (Finset.mem_filter.2 ⟨hv, hc⟩)
    rw [hagree v hns]
    ring
  -- 各項を `log N v` で抑える
  have hbound : ∑ v ∈ T with ch v ∈ Sig,
        (((c.fin v : ℤ) : ℝ) - ((c'.fin v : ℤ) : ℝ)) * Real.log (residueCard v)
      ≤ ∑ v ∈ T with ch v ∈ Sig, Real.log (residueCard v) := by
    refine Finset.sum_le_sum (fun v _ => ?_)
    have h1 : ((c.fin v : ℤ) : ℝ) ≤ 1 := by exact_mod_cast hc1 v
    have h2 : (0 : ℝ) ≤ ((c'.fin v : ℤ) : ℝ) := by exact_mod_cast hc'0 v
    nlinarith [(log_residueCard_pos v).le]
  have hSig : ∑ v ∈ T with ch v ∈ Sig, Real.log (residueCard v)
      ≤ (Module.finrank ℚ F : ℝ) * ∑ q ∈ Sig, Real.log q := by
    refine sum_log_residueCard_le_of_over _ Sig hprime ch (fun v hv => (Finset.mem_filter.1 hv).2)
      (fun v _ => hover v)
  have hkey : deg c - deg c' ≤ (Module.finrank ℚ F : ℝ) * ∑ q ∈ Sig, Real.log q := by
    rw [hdiff, hrestrict]
    exact le_trans hbound hSig
  rw [degNormalized, degNormalized, div_sub_div_same, div_le_iff₀ hn]
  linarith [hkey]

def degNormalized_sub_le_sum_log.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 9,
    item := "Remark 1.5.1(BD-class が (X_ℚ, D_ℚ) だけに依ることの計算部分)",
    sectionId := "genell-rem-1-5-1" }

end ABC3.Found.GenEll
