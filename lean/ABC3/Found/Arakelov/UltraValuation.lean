import ABC3.Found.Arakelov.ArcSpaceInterface
import Mathlib.Order.Filter.FilterProduct
import Mathlib.RingTheory.Valuation.ValuationRing

/-!
# Arakelov (C2) 第 80 ブロック —— **★★★超積の付値環と標準部分**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★★Chow の補題を迂回する道

(C2) は「`X` が ℤ 上固有・平坦 ⟹ `X^arc` はコンパクト」である。
★教科書の道は **Chow の補題**(固有 → 射影的な変更)を経るが、mathlib に無い。

★★**もう一つの道がある**——mathlib は
`AlgebraicGeometry/ValuativeCriterion.lean` に**付値判定法**を持っている
(2026-08-20 実測)。★★★超フィルターの極限を**付値判定法で作る**:

| 段 | 内容 |
|---|---|
| ★ | `Arc X` 上の超フィルター `𝒰` を取る |
| ★★ | 超積体 `*ℂ = ℂ^{Arc X}/𝒰` の点を作る |
| ★★★ | **有限元のなす付値環 `O ⊆ *ℂ`**(本ブロック) |
| ★★★★ | 付値判定法で `Spec O → X` へ持ち上げる |
| ★★★★★ | 閉点の像が `𝒰` の極限になる |

★★★★★★`O` の剰余体は `ℂ` である——**標準部分**(超フィルター極限)がそれを与える。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `IsFinGerm` | ★有限な超積元 |
| `finGermSub` | ★★**有限元のなす付値部分環** |
| `stOf` / `stG` | ★標準部分(超フィルター極限) |
| `stHom` | ★★★**`O →+* ℂ`** |
-/

namespace ABC3.Found.Arakelov

open Filter

variable {α : Type} (𝒰 : Ultrafilter α)

/-! ## ★有限な超積元 -/

/-- ★**有限な超積元**——ノルムが 𝒰-ほとんど至る所で有界。 -/
def IsFinGerm (x : Germ (𝒰 : Filter α) ℂ) : Prop :=
  ∃ C : ℝ, Germ.LiftPred (fun z => ‖z‖ ≤ C) x

theorem isFinGerm_coe (f : α → ℂ) :
    IsFinGerm 𝒰 (f : Germ (𝒰 : Filter α) ℂ) ↔ ∃ C : ℝ, ∀ᶠ a in (𝒰 : Filter α), ‖f a‖ ≤ C := by
  simp [IsFinGerm, Germ.liftPred_coe]

/-- ★★**有限な超積元のなす付値部分環**。

★`x` が有界でなければ `x⁻¹` は無限小——超フィルターだからどちらか一方が起きる。 -/
def finGermSub : ValuationSubring (Germ (𝒰 : Filter α) ℂ) where
  carrier := {x | IsFinGerm 𝒰 x}
  zero_mem' := ⟨0, by
    show Germ.LiftPred _ (((fun _ => (0 : ℂ)) : α → ℂ) : Germ (𝒰 : Filter α) ℂ)
    rw [Germ.liftPred_coe]
    filter_upwards with a
    simp⟩
  one_mem' := ⟨1, by
    show Germ.LiftPred _ (((fun _ => (1 : ℂ)) : α → ℂ) : Germ (𝒰 : Filter α) ℂ)
    rw [Germ.liftPred_coe]
    filter_upwards with a
    simp⟩
  add_mem' := by
    rintro x y ⟨C, hC⟩ ⟨D, hD⟩
    induction x using Germ.inductionOn with
    | _ f =>
      induction y using Germ.inductionOn with
      | _ g =>
        rw [Germ.liftPred_coe] at hC hD
        refine ⟨C + D, ?_⟩
        show Germ.LiftPred _ ((f + g : α → ℂ) : Germ (𝒰 : Filter α) ℂ)
        rw [Germ.liftPred_coe]
        filter_upwards [hC, hD] with a ha hb
        exact le_trans (norm_add_le _ _) (add_le_add ha hb)
  mul_mem' := by
    rintro x y ⟨C, hC⟩ ⟨D, hD⟩
    induction x using Germ.inductionOn with
    | _ f =>
      induction y using Germ.inductionOn with
      | _ g =>
        rw [Germ.liftPred_coe] at hC hD
        refine ⟨C * D, ?_⟩
        show Germ.LiftPred _ ((f * g : α → ℂ) : Germ (𝒰 : Filter α) ℂ)
        rw [Germ.liftPred_coe]
        filter_upwards [hC, hD] with a ha hb
        calc ‖f a * g a‖ = ‖f a‖ * ‖g a‖ := norm_mul _ _
          _ ≤ C * D := mul_le_mul ha hb (norm_nonneg _) (le_trans (norm_nonneg (f a)) ha)
  neg_mem' := by
    rintro x ⟨C, hC⟩
    induction x using Germ.inductionOn with
    | _ f =>
      rw [Germ.liftPred_coe] at hC
      refine ⟨C, ?_⟩
      show Germ.LiftPred _ ((-f : α → ℂ) : Germ (𝒰 : Filter α) ℂ)
      rw [Germ.liftPred_coe]
      filter_upwards [hC] with a ha
      simpa using ha
  mem_or_inv_mem' := by
    intro x
    induction x using Germ.inductionOn with
    | _ f =>
      by_cases h : ∃ C : ℝ, ∀ᶠ a in (𝒰 : Filter α), ‖f a‖ ≤ C
      · left
        obtain ⟨C, hC⟩ := h
        exact ⟨C, by rw [Germ.liftPred_coe]; exact hC⟩
      · right
        push_neg at h
        refine ⟨1, ?_⟩
        show Germ.LiftPred _ (((fun a => (f a)⁻¹) : α → ℂ) : Germ (𝒰 : Filter α) ℂ)
        rw [Germ.liftPred_coe]
        filter_upwards [h 1] with a ha
        rw [norm_inv]
        exact inv_le_one_of_one_le₀ (le_of_lt ha)

/-! ## ★標準部分 -/

/-- ★**標準部分**——有界な族の 𝒰-極限。 -/
noncomputable def stOf (f : α → ℂ) : ℂ := (Filter.map f (𝒰 : Filter α)).lim

theorem stOf_congr {f g : α → ℂ} (h : f =ᶠ[(𝒰 : Filter α)] g) : stOf 𝒰 f = stOf 𝒰 g := by
  unfold stOf
  rw [Filter.map_congr h]

/-- ★★有界なら実際に収束する。 -/
theorem tendsto_stOf {f : α → ℂ} (h : ∃ C : ℝ, ∀ᶠ a in (𝒰 : Filter α), ‖f a‖ ≤ C) :
    Filter.Tendsto f (𝒰 : Filter α) (nhds (stOf 𝒰 f)) := by
  obtain ⟨C, hC⟩ := h
  have hle : (Ultrafilter.map f 𝒰 : Filter ℂ) ≤ Filter.principal (Metric.closedBall (0 : ℂ) C) := by
    rw [Filter.le_principal_iff]
    show {z : ℂ | z ∈ Metric.closedBall (0 : ℂ) C} ∈ Filter.map f (𝒰 : Filter α)
    rw [Filter.mem_map]
    filter_upwards [hC] with a ha
    simpa [Metric.mem_closedBall, dist_zero_right] using ha
  obtain ⟨z, _, hz⟩ :=
    (isCompact_closedBall (0 : ℂ) C).ultrafilter_le_nhds (Ultrafilter.map f 𝒰) hle
  exact le_nhds_lim ⟨z, hz⟩

/-- ★芽の標準部分。 -/
noncomputable def stG (x : Germ (𝒰 : Filter α) ℂ) : ℂ :=
  Quotient.liftOn x (stOf 𝒰) (fun _ _ h => stOf_congr 𝒰 h)

theorem stG_coe (f : α → ℂ) : stG 𝒰 (f : Germ (𝒰 : Filter α) ℂ) = stOf 𝒰 f := rfl

theorem stG_add {x y : Germ (𝒰 : Filter α) ℂ} (hx : IsFinGerm 𝒰 x) (hy : IsFinGerm 𝒰 y) :
    stG 𝒰 (x + y) = stG 𝒰 x + stG 𝒰 y := by
  revert hx hy
  induction x using Germ.inductionOn with
  | _ f =>
    induction y using Germ.inductionOn with
    | _ g =>
      intro hx hy
      rw [isFinGerm_coe] at hx hy
      obtain ⟨C, hC⟩ := hx
      obtain ⟨D, hD⟩ := hy
      have h1 := tendsto_stOf 𝒰 ⟨C, hC⟩
      have h2 := tendsto_stOf 𝒰 ⟨D, hD⟩
      have hb : ∃ E : ℝ, ∀ᶠ a in (𝒰 : Filter α), ‖(f + g) a‖ ≤ E := by
        refine ⟨C + D, ?_⟩
        filter_upwards [hC, hD] with a ha hb
        exact le_trans (norm_add_le _ _) (add_le_add ha hb)
      have h4 := tendsto_stOf 𝒰 hb
      show stOf 𝒰 (f + g) = stOf 𝒰 f + stOf 𝒰 g
      exact tendsto_nhds_unique h4 (h1.add h2)

theorem stG_mul {x y : Germ (𝒰 : Filter α) ℂ} (hx : IsFinGerm 𝒰 x) (hy : IsFinGerm 𝒰 y) :
    stG 𝒰 (x * y) = stG 𝒰 x * stG 𝒰 y := by
  revert hx hy
  induction x using Germ.inductionOn with
  | _ f =>
    induction y using Germ.inductionOn with
    | _ g =>
      intro hx hy
      rw [isFinGerm_coe] at hx hy
      obtain ⟨C, hC⟩ := hx
      obtain ⟨D, hD⟩ := hy
      have h1 := tendsto_stOf 𝒰 ⟨C, hC⟩
      have h2 := tendsto_stOf 𝒰 ⟨D, hD⟩
      have hb : ∃ E : ℝ, ∀ᶠ a in (𝒰 : Filter α), ‖(f * g) a‖ ≤ E := by
        refine ⟨C * D, ?_⟩
        filter_upwards [hC, hD] with a ha hb
        calc ‖f a * g a‖ = ‖f a‖ * ‖g a‖ := norm_mul _ _
          _ ≤ C * D := mul_le_mul ha hb (norm_nonneg _) (le_trans (norm_nonneg (f a)) ha)
      have h4 := tendsto_stOf 𝒰 hb
      show stOf 𝒰 (f * g) = stOf 𝒰 f * stOf 𝒰 g
      exact tendsto_nhds_unique h4 (h1.mul h2)

theorem stG_const (c : ℂ) : stG 𝒰 (((fun _ => c) : α → ℂ) : Germ (𝒰 : Filter α) ℂ) = c := by
  show stOf 𝒰 (fun _ => c) = c
  have h1 := tendsto_stOf 𝒰 (f := fun _ => c) ⟨‖c‖, by filter_upwards with a; exact le_rfl⟩
  exact tendsto_nhds_unique h1 tendsto_const_nhds

/-- ★★★**標準部分の環準同型** `O → ℂ`。 -/
noncomputable def stHom : finGermSub 𝒰 →+* ℂ where
  toFun := fun x => stG 𝒰 (x : Germ (𝒰 : Filter α) ℂ)
  map_one' := stG_const 𝒰 1
  map_mul' := fun x y => stG_mul 𝒰 x.2 y.2
  map_zero' := stG_const 𝒰 0
  map_add' := fun x y => stG_add 𝒰 x.2 y.2

theorem stHom_coe (f : α → ℂ) (hf : IsFinGerm 𝒰 (f : Germ (𝒰 : Filter α) ℂ)) :
    stHom 𝒰 ⟨(f : Germ (𝒰 : Filter α) ℂ), hf⟩ = stOf 𝒰 f := rfl

/-! ## ★出典の紐付け(`.src`) -/

def finGermSub.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 C——固有性からコンパクト性を出す道具:超積の付値環)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
