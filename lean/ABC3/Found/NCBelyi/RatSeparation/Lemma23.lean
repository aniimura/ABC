/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.NCBelyi.Separation

/-!
# RatSeparation —— `[NCBelyi] Lemma 2.3` の分

☆もとの 1 枚を**条なしの項目ごと**に割ったものである（第 1458、案 a）。
★「1 ファイル = 1 ノード」を回復するための分割で、中身は動かしていない。
-/

namespace ABC3.Found.NCBelyi



/-! ## ★`absInvShift` の基本性質 -/

theorem absInvShift_nonneg (lam : ℂ) (p : P1C) : 0 ≤ absInvShift lam p := by
  cases p with
  | none => simp
  | some x => simp

theorem absInvShift_pos_of_ne {lam b : ℂ} (h : lam ≠ b) :
    0 < absInvShift lam (some b) := by
  rw [absInvShift_some]
  refine norm_pos_iff.2 ?_
  refine div_ne_zero one_ne_zero ?_
  exact sub_ne_zero.2 (Ne.symm h)

/-! ## ★★`S` の点と `β` の距離の下界 -/

/-- ★**`S` の有限点はすべて `β` から `δ` 以上離れている**ような `δ > 0` が取れる。

★`Separation.lean` の `lemma_2_3` の証明の前半をそのまま切り出したものである
(`∞` は距離 `1` に潰して有限集合の `min'` を取る)。 -/
theorem exists_dist_lower_bound (S : Finset P1C) (b : ℂ) (hb : (some b) ∉ S) :
    ∃ delta : ℝ, 0 < delta ∧ ∀ p ∈ S, ∀ x : ℂ, p = some x → delta ≤ ‖x - b‖ := by
  classical
  set g : P1C → ℝ := fun p => match p with | none => 1 | some x => ‖x - b‖ with hg
  set D : Finset ℝ := S.image g with hD
  have hgpos : ∀ p ∈ S, 0 < g p := by
    intro p hp
    cases p with
    | none => simp [hg]
    | some x =>
      have hxb : x ≠ b := by
        intro hcon; rw [hcon] at hp; exact hb hp
      simp only [hg]
      exact norm_pos_iff.2 (sub_ne_zero.2 hxb)
  set delta : ℝ := if h : D.Nonempty then D.min' h else 1 with hdel
  have hdelta : 0 < delta := by
    simp only [hdel]
    split_ifs with h
    · obtain ⟨p, hp, hpe⟩ := Finset.mem_image.1 (D.min'_mem h)
      rw [← hpe]
      exact hgpos p hp
    · norm_num
  refine ⟨delta, hdelta, ?_⟩
  intro p hp x hpx
  have hmem : g p ∈ D := Finset.mem_image_of_mem g hp
  have hne : D.Nonempty := ⟨g p, hmem⟩
  have hle : delta ≤ g p := by simp only [hdel, dif_pos hne]; exact D.min'_le _ hmem
  rw [hpx] at hle
  simpa [hg] using hle

/-! ## ★★★★★`Lemma 2.3` の有理版 -/

/-- ★★★★★**[NCBelyi] Lemma 2.3 の追記** —— `β ∈ ℚ` なら `λ ∈ ℚ` に取れる。

原文 (NCBelyi p.4):
> if β ∈Q, then one may take λ ∈Q.

★`separation_core` は `λ ≝ β + ε` を返し、`ε` は
`0 < ε ≤ min(δ/4, δ/(4C))` を満たす**任意の実数**でよい。
★★だから **`ε` を有理数に取ればそのまま `λ ∈ ℚ`** である
——原文の「one may take」はこの自由度のことである。 -/
theorem lemma_2_3_rat (S : Finset P1C) (C : ℝ) (hC : 0 < C) (b : ℚ)
    (hb : (some (b : ℂ)) ∉ S) :
    ∃ lam : ℚ, (lam : ℂ) ≠ (b : ℂ) ∧ (some (lam : ℂ) ∉ S)
      ∧ ∀ p ∈ S, C * absInvShift (lam : ℂ) p
          ≤ absInvShift (lam : ℂ) (some (b : ℂ)) := by
  obtain ⟨delta, hdelta, hSle⟩ := exists_dist_lower_bound S (b : ℂ) hb
  have hmin : (0 : ℝ) < min (delta / 4) (delta / (4 * C)) :=
    lt_min (by positivity) (by positivity)
  obtain ⟨q, hq0, hqlt⟩ := exists_rat_btwn hmin
  refine ⟨b + q, ?_⟩
  have heps : (0 : ℝ) < (q : ℝ) := hq0
  have h1 : (q : ℝ) ≤ delta / 4 := le_trans hqlt.le (min_le_left _ _)
  have h2 : C * (q : ℝ) ≤ delta / 4 := by
    have hq2 : (q : ℝ) ≤ delta / (4 * C) := le_trans hqlt.le (min_le_right _ _)
    calc C * (q : ℝ) ≤ C * (delta / (4 * C)) := mul_le_mul_of_nonneg_left hq2 hC.le
      _ = delta / 4 := by field_simp
  have hcast : (((b + q : ℚ)) : ℂ) = (b : ℂ) + ((q : ℝ) : ℂ) := by push_cast; ring
  rw [hcast]
  exact separation_core S C hC (b : ℂ) delta (q : ℝ) hdelta heps hSle h1 h2

/-! ## ★出典の紐付け(`.src`) -/

def lemma_2_3_rat.src : ABC3.Meta.Source :=
  { paper := "NCBelyi", pdfPage := 4,
    item := "Lemma 2.3(Moreover, if β ∈ ℚ, then one may take λ ∈ ℚ)",
    sectionId := "ncbelyi-lemma-2-3" }

end ABC3.Found.NCBelyi
