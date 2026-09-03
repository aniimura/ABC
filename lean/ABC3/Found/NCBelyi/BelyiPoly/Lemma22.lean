/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.NCBelyi.Lemma22
import ABC3.Found.NCBelyi.DescendData

/-!
# BelyiPoly —— `[NCBelyi] Lemma 2.2` の分

☆もとの 1 枚を**条なしの項目ごと**に割ったものである（第 1458、案 a）。
★「1 ファイル = 1 ノード」を回復するための分割で、中身は動かしていない。
-/

namespace ABC3.Found.NCBelyi

open Polynomial

/-! ## ★条件 (i)(ii) はアフィン変換で作れる -/

/-- ★**`x ↦ x − min S` で `0 ∈ S′` と非負性が出る**(原文の条件 (i)(ii))。 -/
theorem exists_shift_nonneg (S : Finset ℚ) (hne : S.Nonempty) :
    ∃ m : ℚ, (0 : ℚ) ∈ S.image (fun x => x - m)
      ∧ ∀ α ∈ S.image (fun x => x - m), α ≠ 0 → 0 < α := by
  classical
  refine ⟨S.min' hne, Finset.mem_image.2 ⟨S.min' hne, S.min'_mem hne, by ring⟩, ?_⟩
  intro α hα hα0
  obtain ⟨x, hx, rfl⟩ := Finset.mem_image.1 hα
  have hge : S.min' hne ≤ x := S.min'_le x hx
  rcases lt_or_eq_of_le hge with h | h
  · linarith
  · exact absurd (by rw [← h]; ring) hα0

/-- ★★**条件 (iii) を満たす `β` は常に取れる** —— `β ≝ 2·max S + 1`。

原文 (NCBelyi p.4):
> (iii) β/α ≥ 2, for all α ∈ S\{0, ∞}.

★`0 ∈ S` なら `max S ≥ 0` なので `β ≥ 1 > max S` であり、`β ∉ S` も出る。 -/
theorem exists_beta_large (S : Finset ℚ) (hne : S.Nonempty) (h0 : (0 : ℚ) ∈ S) :
    ∃ β : ℚ, β ∉ S ∧ β ≠ 0 ∧ ∀ α ∈ S, α ≠ 0 → 2 * α ≤ β := by
  classical
  have hmax0 : (0 : ℚ) ≤ S.max' hne := S.le_max' 0 h0
  refine ⟨2 * S.max' hne + 1, ?_, by linarith, ?_⟩
  · intro hmem
    have := S.le_max' _ hmem
    linarith
  · intro α hα _
    have := S.le_max' α hα
    linarith

/-! ## ★出典の紐付け(`.src`) -/

def exists_shift_nonneg.src : ABC3.Meta.Source :=
  { paper := "NCBelyi", pdfPage := 4,
    item := "Lemma 2.2 の条件 (i)(ii)(アフィン変換で作れること)",
    sectionId := "ncbelyi-lemma-2-2" }

def exists_beta_large.src : ABC3.Meta.Source :=
  { paper := "NCBelyi", pdfPage := 4,
    item := "Lemma 2.2 の条件 (iii)(β ≝ 2·max S + 1 で満たせること)",
    sectionId := "ncbelyi-lemma-2-2" }

end ABC3.Found.NCBelyi
