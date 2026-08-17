/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import Mathlib.Analysis.Complex.AbsMax
import Mathlib.Analysis.Complex.RemovableSingularity
import Mathlib.Analysis.Calculus.DSlope

/-!
# 多点で消える正則関数の小ささ —— チェーン `sixexp` の葉 `schwarz-many-zeros`

★`ResearchPaper/frdi-decomposition.json` の `sixexp` チェーンの葉。
最終目標は [FrdI] `Lemma 6.5, (ii)`(six exponentials theorem)。

## ★何を示すか

`f` が `|z| ≤ R` で正則、`|z| ≤ r` の中の **相異なる `N` 点で消える**、
`|z| = R` で `|f| ≤ M` とすると、`|z| ≤ r` で

  ★★★**`|f(z)| ≤ M · (2r/(R−r))^N`**   (`norm_le_mul_pow_of_forall_eq_zero`)

★これが超越性証明で「解析側の小ささ」を作る段である。零点を増やすほど、
また `R` を大きくするほど小さくなる。`5r ≤ R` なら係数は `1/2` 以下
(`norm_le_half_pow`)。

## ★★測定(2026-08-18)—— mathlib にあったもの・無かったもの

| 要るもの | mathlib | 判定 |
|---|---|---|
| 最大値原理 | `Complex.norm_le_of_forall_mem_frontier_norm_le`(`Analysis/Complex/AbsMax.lean:400`) | ★**ある** |
| 除去可能特異点(`f/(z−a)` が正則) | `Complex.differentiableOn_dslope`(`Analysis/Complex/RemovableSingularity.lean:60`) | ★**ある** |
| `(b−a)·dslope f a b = f b`(`f a = 0` のとき) | `sub_smul_dslope_of_zero` | ★**ある** |
| **多零点版の小ささ評価** | `Analysis/Complex/` の宣言一覧を目視、1 点の高位版はあるが多点版は無し | ★**無い**(本ファイル) |

★★**すなわち、部品はすべて在庫にあった。**組み方(`dslope` による零点の 1 個ずつの
消去と、`Finset` についての帰納)が我々の仕事である。

## ★設計の要 —— なぜ `Finset ℂ` で添字づけるか

★零点を `Fin n → ℂ` で持つと**相異なることを別に仮定**しなければならず、
帰納の各段で「残りの零点が新しい関数でも零点である」ことを言うのに
その仮定を使い回すことになる。★`Finset ℂ` にすると**相異なることが型から出る**ので、
`Finset.induction_on` の `a ∉ s` がそのまま使える。
-/

namespace ABC3.Found.SixExp

open Metric Complex

/-- 最大値原理(閉円板版)——球面での上界が内部でも成り立つ。 -/
theorem norm_le_of_sphere_bound {f : ℂ → ℂ} {R M : ℝ} (hR : 0 < R)
    (hf : DifferentiableOn ℂ f (closedBall 0 R))
    (hM : ∀ ζ : ℂ, ‖ζ‖ = R → ‖f ζ‖ ≤ M)
    {z : ℂ} (hz : ‖z‖ ≤ R) : ‖f z‖ ≤ M := by
  have hcl : closure (ball (0:ℂ) R) = closedBall 0 R := closure_ball 0 hR.ne'
  refine Complex.norm_le_of_forall_mem_frontier_norm_le (U := ball (0:ℂ) R) isBounded_ball ?_ ?_ ?_
  · exact ⟨hf.mono ball_subset_closedBall, by rw [hcl]; exact hf.continuousOn⟩
  · intro ζ hζ
    rw [frontier_ball 0 hR.ne'] at hζ
    exact hM ζ (by simpa [mem_sphere_iff_norm] using hζ)
  · rw [hcl]; simpa [mem_closedBall_zero_iff] using hz

/-- ★★★**多零点版の Schwarz の補題**。

`f` が `|z| ≤ R` で正則、`s ⊆ {|z| ≤ r}` の各点で消え、`|z| = R` で `|f| ≤ M` なら、
`|z| ≤ r` で `|f(z)| ≤ M · (2r/(R−r))^{|s|}`。

★証明は `s` についての帰納。各段で `dslope f a`(= 除去可能特異点をならした `f/(z−a)`)
に移ると、零点が 1 個減り、球面での上界は `M/(R−r)` になり、
`|z| ≤ r` では `|z−a| ≤ 2r` を掛けて戻る。 -/
theorem norm_le_mul_pow_of_forall_eq_zero {R r : ℝ} (hr : 0 ≤ r) (hrR : r < R) :
    ∀ (s : Finset ℂ) (f : ℂ → ℂ) (M : ℝ),
      DifferentiableOn ℂ f (closedBall 0 R) →
      (∀ w ∈ s, ‖w‖ ≤ r) → (∀ w ∈ s, f w = 0) →
      (∀ ζ : ℂ, ‖ζ‖ = R → ‖f ζ‖ ≤ M) →
      ∀ z : ℂ, ‖z‖ ≤ r → ‖f z‖ ≤ M * (2 * r / (R - r)) ^ s.card := by
  have hR : 0 < R := lt_of_le_of_lt hr hrR
  have hRr : (0:ℝ) < R - r := by linarith
  intro s
  induction s using Finset.induction_on with
  | empty =>
    intro f M hf _ _ hM z hz
    simpa using norm_le_of_sphere_bound hR hf hM (le_trans hz hrR.le)
  | insert a s ha ih =>
    intro f M hf hs hzero hM z hzr
    have hM0 : (0:ℝ) ≤ M := by
      refine le_trans (norm_nonneg (f (R : ℂ))) (hM (R : ℂ) ?_)
      simp [abs_of_pos hR]
    have haf : ‖a‖ ≤ r := hs a (Finset.mem_insert_self a s)
    have hfa : f a = 0 := hzero a (Finset.mem_insert_self a s)
    -- ★除去可能特異点: `dslope f a` は閉円板全体で正則
    have hnhds : closedBall (0:ℂ) R ∈ nhds a :=
      mem_nhds_iff.mpr ⟨ball 0 R, ball_subset_closedBall, isOpen_ball,
        mem_ball_zero_iff.mpr (lt_of_le_of_lt haf hrR)⟩
    have hg : DifferentiableOn ℂ (dslope f a) (closedBall 0 R) :=
      (Complex.differentiableOn_dslope hnhds).mpr hf
    have hgs : ∀ w ∈ s, ‖w‖ ≤ r := fun w hw => hs w (Finset.mem_insert_of_mem hw)
    -- ★残りの零点は `dslope f a` の零点でもある(`a ∉ s` が効く)
    have hgzero : ∀ w ∈ s, dslope f a w = 0 := by
      intro w hw
      have hne : w ≠ a := fun h => ha (h ▸ hw)
      rw [dslope_of_ne _ hne, slope_def_field, hfa, sub_zero,
        hzero w (Finset.mem_insert_of_mem hw), zero_div]
    -- ★球面での上界は `M/(R−r)` に落ちる
    have hgM : ∀ ζ : ℂ, ‖ζ‖ = R → ‖dslope f a ζ‖ ≤ M / (R - r) := by
      intro ζ hζ
      have hne : ζ ≠ a := by
        intro h
        rw [h] at hζ
        linarith [hζ ▸ haf]
      rw [dslope_of_ne _ hne, slope_def_field, hfa, sub_zero, norm_div]
      have hden : R - r ≤ ‖ζ - a‖ := by
        have h1 := norm_sub_norm_le ζ a
        rw [hζ] at h1
        linarith [h1, haf]
      gcongr
      exact hM ζ hζ
    have hih := ih (dslope f a) (M / (R - r)) hg hgs hgzero hgM z hzr
    have hfz : (z - a) • dslope f a z = f z := sub_smul_dslope_of_zero hfa z
    have hnorm : ‖f z‖ = ‖z - a‖ * ‖dslope f a z‖ := by rw [← hfz, norm_smul]
    have hza : ‖z - a‖ ≤ 2 * r := by
      calc ‖z - a‖ ≤ ‖z‖ + ‖a‖ := norm_sub_le z a
        _ ≤ r + r := by linarith
        _ = 2 * r := by ring
    rw [hnorm, Finset.card_insert_of_notMem ha]
    calc ‖z - a‖ * ‖dslope f a z‖
        ≤ (2 * r) * (M / (R - r) * (2 * r / (R - r)) ^ s.card) := by gcongr
      _ = M * (2 * r / (R - r)) ^ (s.card + 1) := by rw [pow_succ]; field_simp

/-- ★`5r ≤ R` なら係数は `1/2` 以下 —— 零点の個数だけ指数的に小さくなる。 -/
theorem norm_le_half_pow {R r : ℝ} (hr : 0 ≤ r) (h5 : 5 * r ≤ R) (hrR : r < R)
    (s : Finset ℂ) (f : ℂ → ℂ) (M : ℝ)
    (hf : DifferentiableOn ℂ f (closedBall 0 R))
    (hs : ∀ w ∈ s, ‖w‖ ≤ r) (hzero : ∀ w ∈ s, f w = 0)
    (hM : ∀ ζ : ℂ, ‖ζ‖ = R → ‖f ζ‖ ≤ M)
    (z : ℂ) (hz : ‖z‖ ≤ r) : ‖f z‖ ≤ M * (1/2 : ℝ) ^ s.card := by
  have hM0 : (0:ℝ) ≤ M := by
    have hR : 0 < R := lt_of_le_of_lt hr hrR
    refine le_trans (norm_nonneg (f (R : ℂ))) (hM (R : ℂ) ?_)
    simp [abs_of_pos hR]
  have hRr : (0:ℝ) < R - r := by linarith
  have hfac : 2 * r / (R - r) ≤ 1/2 := by
    rw [div_le_iff₀ hRr]
    linarith
  refine le_trans (norm_le_mul_pow_of_forall_eq_zero hr hrR s f M hf hs hzero hM z hz) ?_
  gcongr

/-- ★**負の対照** —— 「各零点で消える」を落とすと結論は破れる。

`f ≡ 1`・`R = 5`・`r = 1`・`s = {0}` は、`f 0 = 0` 以外の仮定をすべて満たすが、
`‖f 0‖ = 1 > 1/2 = M · (2r/(R−r))^{|s|}` である。 -/
theorem not_norm_le_of_not_zero :
    ¬ (‖(1 : ℂ)‖ ≤ (1:ℝ) * (2 * 1 / ((5:ℝ) - 1)) ^ ({0} : Finset ℂ).card) := by
  norm_num

end ABC3.Found.SixExp
