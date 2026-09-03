/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.NCBelyi.Separation
import ABC3.Found.NCBelyi.RatSeparation.Lemma23

/-!
# RatSeparation —— `[NCBelyi] Lemma 2.4` の分

☆もとの 1 枚を**条なしの項目ごと**に割ったものである（第 1458、案 a）。
★「1 ファイル = 1 ノード」を回復するための分割で、中身は動かしていない。
-/

namespace ABC3.Found.NCBelyi



/-! ## ★★★★★★★`Lemma 2.4` の正規化 -/

/-- ★★★★★★★**正規化** —— `|f(α)| ≤ 1`(`∀ α ∈ S`)かつ `|f(β)| ≥ C`。

原文 (NCBelyi p.5):
> multiplying by some positive rational number, we may assume that |α| ≤1, for all

`f(x) = c/(x − λ)`(`λ, c ∈ ℚ`、`c > 0`)が取れる。

★★★★**`2C` で `Lemma 2.3` を引くのが要点**である。
`C` で引くと倍率は `C/M` ちょうどしか許されず有理数に取れないが、
`2C` で引けば `c ∈ (C/M, 2C/M)` の開区間が空でなくなり `ℚ` の稠密性が効く。 -/
theorem exists_rat_normalization (S : Finset P1C) (C : ℝ) (hC : 0 < C) (b : ℚ)
    (hb : (some (b : ℂ)) ∉ S) :
    ∃ lam c : ℚ, 0 < c ∧ (lam : ℂ) ≠ (b : ℂ) ∧ (some (lam : ℂ) ∉ S)
      ∧ (∀ p ∈ S, (c : ℝ) * absInvShift (lam : ℂ) p ≤ 1)
      ∧ C ≤ (c : ℝ) * absInvShift (lam : ℂ) (some (b : ℂ)) := by
  obtain ⟨lam, hne, hnotmem, hsep⟩ := lemma_2_3_rat S (2 * C) (by linarith) b hb
  have hM : 0 < absInvShift (lam : ℂ) (some (b : ℂ)) := absInvShift_pos_of_ne hne
  have hlt : C / absInvShift (lam : ℂ) (some (b : ℂ))
      < 2 * C / absInvShift (lam : ℂ) (some (b : ℂ)) := by
    rw [div_lt_div_iff_of_pos_right hM]
    linarith
  obtain ⟨c, hc1, hc2⟩ := exists_rat_btwn hlt
  have hcpos : (0 : ℝ) < (c : ℝ) := lt_trans (by positivity) hc1
  refine ⟨lam, c, by exact_mod_cast hcpos, hne, hnotmem, ?_, ?_⟩
  · intro p hp
    have hle : absInvShift (lam : ℂ) p
        ≤ absInvShift (lam : ℂ) (some (b : ℂ)) / (2 * C) := by
      rw [le_div_iff₀ (by linarith)]
      have := hsep p hp
      linarith
    have hcle : (c : ℝ) ≤ 2 * C / absInvShift (lam : ℂ) (some (b : ℂ)) := hc2.le
    calc (c : ℝ) * absInvShift (lam : ℂ) p
        ≤ (2 * C / absInvShift (lam : ℂ) (some (b : ℂ)))
            * (absInvShift (lam : ℂ) (some (b : ℂ)) / (2 * C)) :=
          mul_le_mul hcle hle (absInvShift_nonneg _ _) (by positivity)
      _ = 1 := by field_simp
  · have hlt2 := (div_lt_iff₀ hM).1 hc1
    linarith

def exists_rat_normalization.src : ABC3.Meta.Source :=
  { paper := "NCBelyi", pdfPage := 5,
    item := "Lemma 2.4(正規化——有理係数の自己同型と正の有理数倍)",
    sectionId := "ncbelyi-lemma-2-4" }

end ABC3.Found.NCBelyi
