/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.Divisor.ArithPerfFactorial

/-!
# 実係数の `Φ` —— 台(`SuppElt`)と「`p` で正・`p` を含まない」分割

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.114。

原文 (FrdI p.115):
> C is of [strictly] rational type, and that every object of (Cun-tr)birat is Frobenius-

## ★★幾何版(`CartierPerfFactorial.lean`)の 2 本を実係数へ

| 幾何(`ℤ` 係数) | 算術(`ℝ` 係数、本ファイル) |
|---|---|
| `mem_suppElt_iff` | `mem_suppElt_iff_R` |
| `exists_split_of_qc` | `exists_split_R` |
| `exists_split_suppElt_of_qc` | `exists_split_suppElt_R` |

## ★★★幾何より**易しい**理由

幾何では「Cartier 因子の正部分は Cartier とは限らない」ため、
`K`-`Q`-Cartier(`n_t·[t] ∈ Γ`)で負の座標を埋める必要があった。
★算術の `Γ = arithDivGroup L` は**座標ごとの条件**(`IsCoordwiseR`)なので、
`single t (y t) ∈ Γ` がそのまま使え、**`|y t|` で埋めるだけで済む**。
-/

namespace ABC3.Found.Divisor

open Finsupp ABC3.Found.FrdI

open scoped NNReal

variable {S : Type*} {Γ : AddSubgroup (S →₀ ℝ)}

/-! ## ★1. `SuppElt` は「その素点の係数が 0 でない」ことに他ならない -/

/-- ★★**`SuppElt` は係数が 0 でないことに他ならない**(実係数版)。 -/
theorem mem_suppElt_iff_R (hC : IsCoordwiseR Γ) (hG : IsGenSubgroupR Γ)
    (a : effR Γ) (p : Prime (effR Γ)) :
    p ∈ SuppElt (iotaR hG) a ↔ ((a : effR Γ) : S →₀ ℝ) (effRPrimeEquiv hG p) ≠ 0 := by
  have hkey : realCoeffAt hG p (Pf.mk a 1)
      = ((a : effR Γ) : S →₀ ℝ) (effRPrimeEquiv hG p) := by
    show (pfCoeffRHom (Pf.mk a 1) : S →₀ ℝ) (effRPrimeEquiv hG p) = _
    exact congrArg (fun t : S →₀ ℝ => t (effRPrimeEquiv hG p)) (pfCoeffR_of a)
  show factorMap (iotaR hG) (Pf.mk a 1) p ≠ 0 ↔ _
  rw [factorMap_iotaR hC hG (Pf.mk a 1) p, Ne, iotaR_eq_zero_iff hG p (Pf.mk a 1), hkey]

/-- ★locator —— `Definition 2.4, (i), (d)` の台と係数の対応(実係数)。 -/
def mem_suppElt_iff_R.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 114,
    item := "Theorem 6.4, (i) — SuppElt は係数が 0 でないこと(実係数)",
    sectionId := "frdi-thm-6-4" }

/-! ## ★2. `Γ` の元を「`s` で正・`s` を含まない」差に割る -/

section Split

variable [DecidableEq S]

/-- ★★★★**座標ごとに閉じていれば、`Γ` の元は「`s` で正・`s` を含まない」差に書ける**。

★埋め草は `b := ∑_{t ∈ supp y \ {s}} [t]·|y t|` である。
`|y t| + y t ≥ 0` なので `a := b + y` は有効、`b s = 0`、`a s = y s > 0`。 -/
theorem exists_split_R (hC : IsCoordwiseR Γ) {y : S →₀ ℝ} (hy : y ∈ Γ)
    {s : S} (hs : 0 < y s) :
    ∃ a b : S →₀ ℝ, a ∈ effR Γ ∧ b ∈ effR Γ ∧ a = b + y ∧ 0 < a s ∧ b s = 0 := by
  classical
  set b : S →₀ ℝ := ∑ t ∈ y.support.erase s, single t |y t| with hb
  have hsingle : ∀ t : S, single t |y t| ∈ Γ := by
    intro t
    rcases abs_choice (y t) with hch | hch
    · rw [hch]; exact hC y hy t
    · rw [hch, Finsupp.single_neg]; exact Γ.neg_mem (hC y hy t)
  have hbmem : b ∈ Γ := AddSubgroup.sum_mem _ fun t _ => hsingle t
  have hbapp : ∀ u : S, b u = if u ∈ y.support.erase s then |y u| else 0 := by
    intro u
    rw [hb, Finsupp.finsetSum_apply]
    by_cases hu : u ∈ y.support.erase s
    · rw [if_pos hu, Finset.sum_eq_single u]
      · rw [Finsupp.single_eq_same]
      · intro t _ htu
        exact Finsupp.single_eq_of_ne (Ne.symm htu)
      · intro hcon; exact absurd hu hcon
    · rw [if_neg hu, Finset.sum_eq_zero]
      intro t ht
      refine Finsupp.single_eq_of_ne ?_
      intro hc
      exact hu (hc ▸ ht)
  have hbnn : ∀ u, 0 ≤ b u := by
    intro u
    rw [hbapp u]
    by_cases hu : u ∈ y.support.erase s
    · rw [if_pos hu]; exact abs_nonneg _
    · rw [if_neg hu]
  have hbs : b s = 0 := by
    rw [hbapp s, if_neg (by simp)]
  refine ⟨b + y, b, mem_effR.mpr ⟨Γ.add_mem hbmem hy, ?_⟩,
    mem_effR.mpr ⟨hbmem, hbnn⟩, rfl, ?_, hbs⟩
  · intro u
    rw [Finsupp.add_apply]
    by_cases hu : u ∈ y.support.erase s
    · rw [hbapp u, if_pos hu]
      have h1 := neg_abs_le (y u)
      linarith
    · rw [hbapp u, if_neg hu]
      rcases eq_or_ne u s with rfl | hus
      · linarith
      · have hy0 : y u = 0 := by
          by_contra hc
          exact hu (Finset.mem_erase.mpr ⟨hus, Finsupp.mem_support_iff.mpr hc⟩)
        rw [hy0]
        linarith
  · rw [Finsupp.add_apply, hbs]
    linarith

/-- ★★★locator —— `Theorem 6.4, (i)` の rationally standard に要る分割(実係数)。 -/
def exists_split_R.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 114,
    item := "Theorem 6.4, (i) — 座標ごとに閉じていれば Γ の元は p で正・p を含まない差に割れる",
    sectionId := "frdi-thm-6-4" }

/-- ★★★★**`IsStrictlyRational` が要求する分割**(台の言葉で、実係数版)。 -/
theorem exists_split_suppElt_R (hC : IsCoordwiseR Γ) (hG : IsGenSubgroupR Γ)
    {y : S →₀ ℝ} (hy : y ∈ Γ) (p : Prime (effR Γ))
    (hyp : y (effRPrimeEquiv hG p) ≠ 0) :
    ∃ (a b : effR Γ) (z : S →₀ ℝ), z ∈ Γ ∧ (z = y ∨ z = -y) ∧
      ((a : effR Γ) : S →₀ ℝ) = ((b : effR Γ) : S →₀ ℝ) + z ∧
      p ∈ SuppElt (iotaR hG) a ∧ p ∉ SuppElt (iotaR hG) b := by
  classical
  set s := effRPrimeEquiv hG p with hsdef
  rcases lt_or_gt_of_ne hyp with hneg | hpos
  · have hnegy : 0 < (-y) s := by
      show 0 < -(y s)
      linarith
    obtain ⟨a, b, ha, hb, hab, hap, hbp⟩ := exists_split_R hC (Γ.neg_mem hy) hnegy
    refine ⟨⟨a, ha⟩, ⟨b, hb⟩, -y, Γ.neg_mem hy, Or.inr rfl, hab, ?_, ?_⟩
    · exact (mem_suppElt_iff_R hC hG ⟨a, ha⟩ p).mpr (ne_of_gt hap)
    · intro hc
      exact ((mem_suppElt_iff_R hC hG ⟨b, hb⟩ p).mp hc) hbp
  · obtain ⟨a, b, ha, hb, hab, hap, hbp⟩ := exists_split_R hC hy hpos
    refine ⟨⟨a, ha⟩, ⟨b, hb⟩, y, hy, Or.inl rfl, hab, ?_, ?_⟩
    · exact (mem_suppElt_iff_R hC hG ⟨a, ha⟩ p).mpr (ne_of_gt hap)
    · intro hc
      exact ((mem_suppElt_iff_R hC hG ⟨b, hb⟩ p).mp hc) hbp

/-- ★★★★locator —— `Theorem 6.4, (i)` の分割(台の言葉、実係数)。 -/
def exists_split_suppElt_R.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 114,
    item := "Theorem 6.4, (i) — Φ^birat の元は p で正・p を含まない差に割れる(実係数)",
    sectionId := "frdi-thm-6-4" }

end Split

end ABC3.Found.Divisor
