import ABC3.Found.GenEll.LogDiffValue
import Mathlib.NumberTheory.NumberField.Discriminant.Basic

/-!
# [GenEll] Definition 1.5, (iii) の帰結 —— log-diff の Northcott 性(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.8。

原文 (GenEll p.8):
> determines a well-deﬁned log-diﬀerent function log-diﬀX on X(Q).

## ★★何を取るか

> **`log-diff(F) ≤ C` かつ `[F:ℚ] ≤ d` を満たす数体 `F ⊆ ℚ̄` は有限個である。**

★[GenEll] の `Corollary 4.4` 型の主張(`X(ℚ̄)^{≤d}` 上で
`log-diff + log-cond` を高さで抑える)は、**次数を有界にしたうえで**
log-diff を有界にする。その裏にあるのが本定理である。

★★**次数の有界性は落とせない。** `log-diff` だけを有界にしても体は無限にある
(`[F:ℚ]` がいくらでも大きくなれば `log|disc F|` も大きくてよい)。
本ファイルの `d` はそのために要る——原文の `X(ℚ̄)^{≤d}` の `d` と同じものである。

## ★機構

`LogDiffValue.lean` の **`log-diff(F) = log|disc F| / [F:ℚ]`** に
`[F:ℚ] ≤ d` を掛けると `log|disc F| ≤ C·d`、すなわち `|disc F| ≤ exp(C·d)`。
あとは **Hermite の定理**(`NumberField.finite_of_discr_bdd`、mathlib にある)である。

★mathlib は Hermite を**固定した拡大 `A` の中間体**の言葉で持っている。
原文の `ℚ̄` の中で考えるという設定とそのまま合う。
-/

namespace ABC3.Found.GenEll

open NumberField

section LogDiffFinite

variable (A : Type*) [Field A] [CharZero A]

/-- ★★**log-diff の Northcott 性**。

> `log-diff(F) ≤ C` かつ `[F:ℚ] ≤ d` なる `F ⊆ A` は**有限個**。

★`C` が負でも成り立つ(そのときは集合が空である)——
`|disc F| ≥ 1` より `log-diff ≥ 0` だからである。
証明では `max C 0` を使ってその場合分けを回避している。 -/
theorem finite_logDiffOfField_le (C : ℝ) (d : ℕ) :
    {K : {F : IntermediateField ℚ A // FiniteDimensional ℚ F} |
      haveI : NumberField K := @NumberField.mk _ _ inferInstance K.prop
      logDiffOfField K ≤ C ∧ Module.finrank ℚ K ≤ d}.Finite := by
  classical
  set N : ℕ := ⌈Real.exp (max C 0 * d)⌉₊ with hN
  refine Set.Finite.subset (NumberField.finite_of_discr_bdd A N) ?_
  rintro ⟨K, hK₀⟩ ⟨hC, hd⟩
  haveI : NumberField K := @NumberField.mk _ _ inferInstance hK₀
  -- `|disc K|` を実数の不等式へ
  have hdisc : (NumberField.discr K) ≠ 0 := NumberField.discr_ne_zero K
  have hnat : (NumberField.discr K).natAbs ≠ 0 := Int.natAbs_ne_zero.2 hdisc
  have hpos : (0 : ℝ) < ((NumberField.discr K).natAbs : ℝ) := by
    exact_mod_cast Nat.pos_of_ne_zero hnat
  have hrk : (0 : ℝ) < (Module.finrank ℚ K : ℝ) := by
    exact_mod_cast Module.finrank_pos (R := ℚ) (M := K)
  -- `log|disc K| = log-diff(K)·[K:ℚ] ≤ (max C 0)·d`
  have heq : Real.log ((NumberField.discr K).natAbs)
      = logDiffOfField K * (Module.finrank ℚ K : ℝ) := by
    rw [logDiffOfField_eq]
    field_simp
  have hCmax : logDiffOfField K ≤ max C 0 := le_trans hC (le_max_left _ _)
  have hstep1 : logDiffOfField K * (Module.finrank ℚ K : ℝ)
      ≤ max C 0 * (Module.finrank ℚ K : ℝ) :=
    mul_le_mul_of_nonneg_right hCmax hrk.le
  have hstep2 : max C 0 * (Module.finrank ℚ K : ℝ) ≤ max C 0 * (d : ℝ) := by
    refine mul_le_mul_of_nonneg_left ?_ (le_max_right _ _)
    exact_mod_cast hd
  have hlog : Real.log ((NumberField.discr K).natAbs) ≤ max C 0 * (d : ℝ) := by
    rw [heq]; linarith
  -- 対数を外す
  have hle : ((NumberField.discr K).natAbs : ℝ) ≤ Real.exp (max C 0 * (d : ℝ)) :=
    (Real.log_le_iff_le_exp hpos).1 hlog
  have hleN : (NumberField.discr K).natAbs ≤ N := by
    have hcast : ((NumberField.discr K).natAbs : ℝ) ≤ (N : ℝ) :=
      le_trans hle (Nat.le_ceil _)
    exact_mod_cast hcast
  have : |NumberField.discr K| ≤ (N : ℤ) := by
    rw [Int.abs_eq_natAbs]
    exact_mod_cast hleN
  exact this

end LogDiffFinite

/-! ## ★出典の紐付け(`.src`) -/

def finite_logDiffOfField_le.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 8,
    item := "Definition 1.5, (iii)(log-diff の Northcott 性)",
    sectionId := "genell-def-1-5" }

end ABC3.Found.GenEll
