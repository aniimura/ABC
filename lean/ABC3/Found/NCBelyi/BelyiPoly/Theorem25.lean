/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.NCBelyi.Lemma22
import ABC3.Found.NCBelyi.DescendData
import ABC3.Found.NCBelyi.BelyiPoly.Lemma22

/-!
# BelyiPoly —— `[NCBelyi] Theorem 2.5` の分

☆もとの 1 枚を**条なしの項目ごと**に割ったものである（第 1458、案 a）。
★「1 ファイル = 1 ノード」を回復するための分割で、中身は動かしていない。
-/

namespace ABC3.Found.NCBelyi

open Polynomial

/-! ## ★★★★★★★★★★★★`ℚ` の有限集合に対する Belyi 多項式 -/

/-- ★★★★★★★★★★★★**有理点の有限集合には Belyi 多項式がある**。

原文 (NCBelyi p.6):
> by applying Lemma 2.4, one reduces to the case X = P1

`S ⊆ ℚ` が有限なら、**非定数 `f ∈ ℚ[x]` が存在して**
`f(S) ⊆ {0,1}` かつ **臨界値が `{0,1}` に入る**。

★これは古典的な Belyi の定理の `ℙ¹` 版((a) と (c))である。
★★条件 (i)(ii)(iii) はアフィン変換で作れる(`exists_shift_nonneg` / `exists_beta_large`)ので、
`Lemma 2.2` をそのまま適用できる。
★★★**1 次式は臨界点を持たない**ので、前処理は不分岐性を壊さない。 -/
theorem exists_belyi_poly_rat (S : Finset ℚ) :
    ∃ f : ℚ[X], 0 < f.natDegree
      ∧ (∀ α ∈ S, f.eval α = 0 ∨ f.eval α = 1)
      ∧ (∀ x : ℂ, (derivative (f.map (algebraMap ℚ ℂ))).eval x = 0 →
          (f.map (algebraMap ℚ ℂ)).eval x = 0 ∨ (f.map (algebraMap ℚ ℂ)).eval x = 1) := by
  classical
  rcases S.eq_empty_or_nonempty with rfl | hne
  · refine ⟨X, by simp, fun α hα => absurd hα (by simp), fun x hx => ?_⟩
    exfalso
    rw [Polynomial.map_X, derivative_X, eval_one] at hx
    exact one_ne_zero hx
  obtain ⟨m, h0, hpos⟩ := exists_shift_nonneg S hne
  have hne' : (S.image (fun x => x - m)).Nonempty := ⟨0, h0⟩
  obtain ⟨β, hβ, hβ0, hratio⟩ := exists_beta_large (S.image (fun x => x - m)) hne' h0
  obtain ⟨g, hgdeg, hgS, hgβ0, hgβ1, hgcrit⟩ :=
    lemma_2_2 (S.image (fun x => x - m)) β h0 hpos hβ hβ0 hratio
  -- ★`f ≝ g ∘ (x − m)`
  have hmc : (X - C m : ℚ[X]).natDegree = 1 := natDegree_X_sub_C m
  have hmapc : ((X - C m : ℚ[X]).map (algebraMap ℚ ℂ))
      = (X - C (algebraMap ℚ ℂ m) : ℂ[X]) := by
    simp [Polynomial.map_sub]
  refine ⟨g.comp (X - C m), ?_, ?_, ?_⟩
  · rw [natDegree_comp, hmc, mul_one]
    exact hgdeg
  · intro α hα
    rw [eval_comp, eval_sub, eval_X, eval_C]
    exact hgS _ (Finset.mem_image.2 ⟨α, hα, rfl⟩)
  · intro x hx
    rw [Polynomial.map_comp, hmapc] at hx ⊢
    rw [derivative_comp, eval_mul, eval_comp, derivative_sub, derivative_X, derivative_C,
      sub_zero, eval_one, one_mul, eval_sub, eval_X, eval_C] at hx
    rw [eval_comp, eval_sub, eval_X, eval_C]
    exact hgcrit _ hx

/-! ## ★★★★★★★★★★★★★代数的数の有限集合に対する Belyi 多項式 -/

/-- ★★★★★★★★★★★★★**代数的数の有限集合には Belyi 多項式がある**。

原文 (NCBelyi p.6):
> by applying Lemma 2.4, one reduces to the case X = P1

`S ⊆ ℚ̄` が有限なら、**非定数 `F ∈ ℚ[x]` が存在して**
`F(S) ⊆ {0,1}` かつ **臨界値が `{0,1}` に入る**。

★これが `Theorem 2.5` の `X = ℙ¹` の場合の (a)(c) である。
★★2 段の合成である:

| 段 | 何をする | どこ |
|---|---|---|
| `Lemma 2.4` の多項式の段 | 代数的数を `ℚ` の中へ(臨界値も込み) | `DescendData.lean`(第 417–418) |
| `Lemma 2.2` | `ℚ` の有限集合を `{0,1}` へ | `exists_belyi_poly_rat`(上) |

★★★合成の臨界点は 2 種類に分かれ、どちらも吸収される
——`h′(w) = 0` なら `h(w)` は `h` の臨界値なので `S₁` に入れてあり、
`f′(h(w)) = 0` なら `f` の条件が効く。 -/
theorem exists_belyi_poly (S : Finset ℂ) (hint : ∀ x ∈ S, IsIntegral ℚ x) :
    ∃ F : ℚ[X], 0 < F.natDegree
      ∧ (∀ x ∈ S, Polynomial.aeval x F = 0 ∨ Polynomial.aeval x F = 1)
      ∧ (∀ w : ℂ, Polynomial.aeval w (derivative F) = 0 →
          Polynomial.aeval w F = 0 ∨ Polynomial.aeval w F = 1) := by
  classical
  obtain ⟨h, hdeg, hval, hcrit⟩ := exists_poly_image_rat_crit S hint
  set T : Finset ℂ := S.image (fun x => Polynomial.aeval x h)
      ∪ (((derivative h).map (algebraMap ℚ ℂ)).roots.toFinset.image
          (fun w => Polynomial.aeval w h)) with hT
  set S₁ : Finset ℚ := T.preimage (algebraMap ℚ ℂ)
      ((algebraMap ℚ ℂ).injective.injOn) with hS₁
  have hcastq : ∀ q : ℚ, (algebraMap ℚ ℂ) q = (q : ℂ) := by intro q; simp
  have hmemS₁ : ∀ (q : ℚ), ((q : ℂ)) ∈ T → q ∈ S₁ := by
    intro q hq
    rw [hS₁, Finset.mem_preimage, hcastq]
    exact hq
  obtain ⟨f, hfdeg, hfS, hfcrit⟩ := exists_belyi_poly_rat S₁
  have hbridge : ∀ (p : ℚ[X]) (w : ℂ),
      Polynomial.aeval w p = (p.map (algebraMap ℚ ℂ)).eval w := by
    intro p w
    rw [Polynomial.eval_map, Polynomial.aeval_def]
  have hratval : ∀ q : ℚ, Polynomial.aeval ((q : ℂ)) f = ((f.eval q : ℚ) : ℂ) := by
    intro q
    have h1 : ((q : ℂ)) = (algebraMap ℚ ℂ) q := by simp
    rw [hbridge, h1, eval_map_ratCast]
    simp
  refine ⟨f.comp h, ?_, ?_, ?_⟩
  · rw [natDegree_comp]
    exact Nat.mul_pos hfdeg hdeg
  · -- ★(a) `F(S) ⊆ {0,1}`
    intro x hx
    obtain ⟨q, hq⟩ := hval x hx
    have hqT : ((q : ℂ)) ∈ T := by
      rw [hT, Finset.mem_union]
      exact Or.inl (Finset.mem_image.2 ⟨x, hx, hq⟩)
    have hqS₁ : q ∈ S₁ := hmemS₁ q hqT
    rw [Polynomial.aeval_comp, hq, hratval]
    rcases hfS q hqS₁ with hv | hv <;> rw [hv] <;> simp
  · -- ★(c) 臨界値が `{0,1}` に入る
    intro w hw
    rw [Polynomial.derivative_comp, map_mul, Polynomial.aeval_comp] at hw
    rw [Polynomial.aeval_comp]
    rcases mul_eq_zero.1 hw with h1 | h2
    · -- `h′(w) = 0`
      obtain ⟨q, hq⟩ := hcrit w h1
      have hd0 : (derivative h) ≠ 0 := by
        intro hc
        have hz := (Polynomial.derivative_eq_zero (p := h)).1 hc
        omega
      have hmapne : ((derivative h).map (algebraMap ℚ ℂ)) ≠ 0 :=
        (Polynomial.map_ne_zero_iff (algebraMap ℚ ℂ).injective).2 hd0
      have hwroot : w ∈ ((derivative h).map (algebraMap ℚ ℂ)).roots := by
        rw [Polynomial.mem_roots hmapne]
        simpa [Polynomial.IsRoot, Polynomial.eval_map, ← Polynomial.aeval_def] using h1
      have hqT : ((q : ℂ)) ∈ T := by
        rw [hT, Finset.mem_union]
        exact Or.inr (Finset.mem_image.2 ⟨w, Multiset.mem_toFinset.2 hwroot, hq⟩)
      have hqS₁ : q ∈ S₁ := hmemS₁ q hqT
      rw [hq, hratval]
      rcases hfS q hqS₁ with hv | hv <;> rw [hv] <;> simp
    · -- `f′(h(w)) = 0`
      have hc := hfcrit (Polynomial.aeval w h) (by
        rw [Polynomial.derivative_map, ← hbridge]
        exact h2)
      rw [← hbridge] at hc
      exact hc


def exists_belyi_poly_rat.src : ABC3.Meta.Source :=
  { paper := "NCBelyi", pdfPage := 6,
    item := "Theorem 2.5(X = ℙ¹ の場合——有理点の有限集合に対する Belyi 多項式)",
    sectionId := "ncbelyi-thm-2-5" }

def exists_belyi_poly.src : ABC3.Meta.Source :=
  { paper := "NCBelyi", pdfPage := 6,
    item := "Theorem 2.5(X = ℙ¹ の場合——代数的数の有限集合に対する Belyi 多項式)",
    sectionId := "ncbelyi-thm-2-5" }

end ABC3.Found.NCBelyi
