/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Divisor.ArithTower

/-!
# 素点の制限は全射、引き戻しは単射(`Example 6.3` の `pull_inj` / `pull_nonneg`)

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.113。

原文 (FrdI p.113):
> finite subsets of V(L). Thus, by Theorem 5.2, (ii), this data determines a [model]

## ★★`ArithDatum` の残り 2 フィールド

`ArithTower.lean` が `pull_id` / `pull_comp` を与えたので、残るのは
`pull_nonneg`(有効性を保つ)と `pull_inj`(単射)である。

★`pull_nonneg` は式 `(arithExtend d) V = localDeg V · d(V|_L)` から直ちに
(`localDeg > 0`)。
★★`pull_inj` は **`resPlace` が全射**であることに帰着する ——
「`L` のどの素点も `M` へ延びる」:

| 種別 | 根拠 |
|---|---|
| 非アルキメデス | `Ideal.exists_ideal_over_prime_of_isIntegral`(整拡大で素イデアルは持ち上がる) |
| アルキメデス | `InfinitePlace.comap_surjective`(代数拡大) |

## ★本ファイルで閉じること

| 定理 | 中身 |
|---|---|
| `resFin_surjective` / `resInf_surjective` / `resPlace_surjective` | 素点の制限は全射 |
| `arithExtend_nonneg` | ★`pull_nonneg` |
| `arithExtend_injective` | ★★**`pull_inj`** |
-/

namespace ABC3.Found.Divisor

open NumberField IsDedekindDomain Ideal

universe u

variable {L M : Type u} [Field L] [Field M] [NumberField L] [NumberField M] [Algebra L M]

/-! ## ★1. 素点の制限は全射 -/

/-- ★**非アルキメデス素点は延びる** —— 整拡大で素イデアルが持ち上がることから。 -/
theorem resFin_surjective : Function.Surjective (resFin (L := L) (M := M)) := by
  intro v
  obtain ⟨P, -, hPp, hPo⟩ :=
    Ideal.exists_ideal_over_prime_of_isIntegral (R := 𝓞 L) (S := 𝓞 M)
      v.maximalIdeal.asIdeal ⊥ (by simp)
  haveI := hPp
  have hPne : P ≠ ⊥ := by
    intro h
    rw [h] at hPo
    exact v.maximalIdeal.ne_bot (by simpa using hPo.symm)
  refine ⟨FinitePlace.mk ⟨P, hPp, hPne⟩, ?_⟩
  rw [resFin, FinitePlace.maximalIdeal_mk]
  have : resHOS (L := L) (⟨P, hPp, hPne⟩ : HeightOneSpectrum (𝓞 M)) = v.maximalIdeal := by
    refine HeightOneSpectrum.ext ?_
    exact hPo
  rw [this, FinitePlace.mk_maximalIdeal]

/-- ★**アルキメデス素点は延びる**(`InfinitePlace.comap_surjective`)。 -/
theorem resInf_surjective : Function.Surjective (resInf (L := L) (M := M)) := by
  intro v
  obtain ⟨w, hw⟩ := InfinitePlace.comap_surjective (k := L) (K := M) v
  exact ⟨w, hw⟩

/-- ★★**素点の制限は全射**。 -/
theorem resPlace_surjective : Function.Surjective (resPlace (L := L) (M := M)) := by
  intro v
  cases v with
  | inl w =>
      obtain ⟨W, hW⟩ := resFin_surjective (L := L) (M := M) w
      exact ⟨Sum.inl W, by simp [hW]⟩
  | inr w =>
      obtain ⟨W, hW⟩ := resInf_surjective (L := L) (M := M) w
      exact ⟨Sum.inr W, by simp [hW]⟩

/-! ## ★2. `pull_nonneg` と `pull_inj` -/

theorem arithExtend_apply_eq (d : ArithPlace L →₀ ℝ) (V : ArithPlace M) :
    arithExtend (L := L) d V = (localDeg (L := L) V : ℝ) * d (resPlace (L := L) V) := rfl

/-- ★**引き戻しは有効性を保つ**(`ArithDatum.pull_nonneg`)。 -/
theorem arithExtend_nonneg {d : ArithPlace L →₀ ℝ} (hd : 0 ≤ d) :
    0 ≤ arithExtend (L := L) (M := M) d := by
  refine Finsupp.le_def.mpr fun V => ?_
  rw [Finsupp.coe_zero, Pi.zero_apply, arithExtend_apply_eq]
  exact mul_nonneg (Nat.cast_nonneg _) (Finsupp.le_def.mp hd _)

/-- ★★★**引き戻しは単射**(`ArithDatum.pull_inj`)——
`localDeg > 0` と `resPlace` の全射性から。 -/
theorem arithExtend_injective :
    Function.Injective (arithExtend (L := L) (M := M)) := by
  intro d e h
  refine Finsupp.ext fun v => ?_
  obtain ⟨V, rfl⟩ := resPlace_surjective (L := L) (M := M) v
  have hV : arithExtend (L := L) (M := M) d V = arithExtend (L := L) (M := M) e V := by rw [h]
  rw [arithExtend_apply_eq, arithExtend_apply_eq] at hV
  have hpos : (0 : ℝ) < (localDeg (L := L) V : ℝ) :=
    Nat.cast_pos.mpr (localDeg_pos (L := L) V)
  exact mul_left_cancel₀ (ne_of_gt hpos) hV

/-! ### ★出典の紐付け -/

/-- ★★★locator —— `Example 6.3` の引き戻しの単射性(素点の制限の全射性)。 -/
def arithExtend_injective.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 113,
    item := "Example 6.3 — 因子の引き戻しは単射(素点の制限が全射だから)",
    sectionId := "frdi-example-6-3" }

end ABC3.Found.Divisor
