/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Meta.Claim
import Mathlib.RingTheory.Ideal.Norm.AbsNorm
import Mathlib.RingTheory.DedekindDomain.Ideal.Lemmas
import Mathlib.NumberTheory.NumberField.Basic
import Mathlib.Analysis.SpecialFunctions.Log.Basic
import Mathlib.Analysis.SpecialFunctions.Pow.Continuity
import Mathlib.Topology.Algebra.InfiniteSum.Basic

/-!
# Dirichlet 密度(鎖 `cheb` の `cheb-density`)

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.116。

原文 (FrdI p.116):
> over a prime pi ∈V(Li), then pi splits completely in Li if and only if deg(Li, vi) =

## ★★何のためか

`Theorem 6.4, (iv)` が使う Tchebotarev 密度定理の**語彙**である。
`ResearchPaper/frdi-decomposition.json` の鎖 `cheb` の節点 `cheb-density`。

  `i(T) = lim_{s→1+} (Σ_{𝔭∈T} N𝔭^{-s}) / log(1/(s-1))`

★mathlib に `DirichletDensity` は **0 件**(2026-08-20 実測)。

## ★本ファイルで閉じること

| 定理 | 中身 |
|---|---|
| `HasDirichletDensity.unique` | 密度は一意(極限の一意性) |
| `tendsto_logInv_atTop` | 分母 `log(1/(s-1))` は `s → 1+` で `+∞` |
| `hasDirichletDensity_zero_of_finite` | **有限集合の密度は 0** |
| `infinite_of_hasDirichletDensity_ne_zero` | ★**密度が 0 でなければ無限集合** |

★★最後の 1 本が `Theorem 6.4, (iv)` の使い方そのものである ——
「完全分解する素点の密度は `1/[L:K] ≠ 0`、ゆえにそういう素点は無限個ある」。
-/

namespace ABC3.Found.NF

open NumberField IsDedekindDomain Filter Topology

variable {K : Type*} [Field K] [NumberField K]

/-! ## ★1. 定義 -/

/-- ★★**Dirichlet 密度** —— `i(T) = lim_{s→1+} (Σ_{𝔭∈T} N𝔭^{-s}) / log(1/(s-1))`。 -/
def HasDirichletDensity (S : Set (HeightOneSpectrum (𝓞 K))) (d : ℝ) : Prop :=
  Tendsto
    (fun s : ℝ =>
      (∑' 𝔭 : S, ((Ideal.absNorm (𝔭 : HeightOneSpectrum (𝓞 K)).asIdeal : ℝ)) ^ (-s))
        / Real.log (1 / (s - 1)))
    (𝓝[>] (1 : ℝ)) (𝓝 d)

/-- ★密度は一意 —— `𝓝[>] 1` は自明でないので極限が一意。 -/
theorem HasDirichletDensity.unique {S : Set (HeightOneSpectrum (𝓞 K))} {d d' : ℝ}
    (h : HasDirichletDensity S d) (h' : HasDirichletDensity S d') : d = d' :=
  tendsto_nhds_unique h h'

/-! ## ★2. 分母は `+∞` へ行く -/

/-- ★`s → 1+` のとき `s - 1 → 0+`。 -/
theorem tendsto_sub_one_nhdsGT :
    Tendsto (fun s : ℝ => s - 1) (𝓝[>] (1 : ℝ)) (𝓝[>] (0 : ℝ)) := by
  refine tendsto_nhdsWithin_of_tendsto_nhds_of_eventually_within _ ?_ ?_
  · have : Tendsto (fun s : ℝ => s - 1) (𝓝 (1 : ℝ)) (𝓝 (1 - 1)) :=
      (continuous_id.sub continuous_const).tendsto 1
    simpa using this.mono_left nhdsWithin_le_nhds
  · filter_upwards [self_mem_nhdsWithin] with s hs
    exact sub_pos.mpr (Set.mem_Ioi.mp hs)

/-- ★★**分母は `s → 1+` で `+∞`**。 -/
theorem tendsto_logInv_atTop :
    Tendsto (fun s : ℝ => Real.log (1 / (s - 1))) (𝓝[>] (1 : ℝ)) atTop := by
  have h1 : Tendsto (fun x : ℝ => 1 / x) (𝓝[>] (0 : ℝ)) atTop := by
    simpa [one_div] using tendsto_inv_nhdsGT_zero
  exact Real.tendsto_log_atTop.comp (h1.comp tendsto_sub_one_nhdsGT)

/-! ## ★3. 有限集合の密度は 0 -/

/-- ★有限集合上の和は `Finset` の和で書ける。 -/
theorem tsum_eq_finsetSum_of_finite {S : Set (HeightOneSpectrum (𝓞 K))} (hS : S.Finite)
    (f : HeightOneSpectrum (𝓞 K) → ℝ) :
    ∑' 𝔭 : S, f (𝔭 : HeightOneSpectrum (𝓞 K)) = ∑ 𝔭 ∈ hS.toFinset, f 𝔭 := by
  classical
  haveI : Fintype S := hS.fintype
  rw [tsum_fintype]
  exact (Finset.sum_subtype hS.toFinset (fun x => hS.mem_toFinset) f).symm

/-- ★★**有限集合の Dirichlet 密度は 0**。

★分子は `s → 1` で収束する有限和、分母は `+∞` へ行くから。 -/
theorem hasDirichletDensity_zero_of_finite {S : Set (HeightOneSpectrum (𝓞 K))}
    (hS : S.Finite) : HasDirichletDensity S 0 := by
  classical
  set f : ℝ → ℝ := fun s =>
    ∑ 𝔭 ∈ hS.toFinset, ((Ideal.absNorm 𝔭.asIdeal : ℝ)) ^ (-s) with hf
  have hnum : ∀ s : ℝ,
      (∑' 𝔭 : S, ((Ideal.absNorm (𝔭 : HeightOneSpectrum (𝓞 K)).asIdeal : ℝ)) ^ (-s)) = f s := by
    intro s
    exact tsum_eq_finsetSum_of_finite hS
      (fun 𝔭 => ((Ideal.absNorm 𝔭.asIdeal : ℝ)) ^ (-s))
  have hcont : Tendsto f (𝓝[>] (1 : ℝ)) (𝓝 (f 1)) := by
    refine Tendsto.mono_left ?_ nhdsWithin_le_nhds
    refine tendsto_finsetSum _ (fun 𝔭 _ => ?_)
    have hne : (Ideal.absNorm 𝔭.asIdeal) ≠ 0 := by
      simpa [Ideal.absNorm_eq_zero_iff] using 𝔭.ne_bot
    have hpos : (0 : ℝ) < (Ideal.absNorm 𝔭.asIdeal : ℝ) :=
      Nat.cast_pos.mpr (Nat.pos_of_ne_zero hne)
    exact ((Real.continuousAt_const_rpow (ne_of_gt hpos)).comp
      (continuous_neg.continuousAt)).tendsto
  have : Tendsto (fun s : ℝ => f s / Real.log (1 / (s - 1))) (𝓝[>] (1 : ℝ)) (𝓝 0) :=
    hcont.div_atTop tendsto_logInv_atTop
  refine this.congr (fun s => ?_)
  rw [hnum s]

/-! ## ★4. 密度が 0 でなければ無限集合 -/

/-- ★★★**密度が 0 でない素点の集合は無限**。

★★これが `Theorem 6.4, (iv)` の使い方そのものである ——
「完全分解する素点の密度は `1/[L:K] ≠ 0`、ゆえにそういう素点は無限個ある」。 -/
theorem infinite_of_hasDirichletDensity_ne_zero {S : Set (HeightOneSpectrum (𝓞 K))} {d : ℝ}
    (h : HasDirichletDensity S d) (hd : d ≠ 0) : S.Infinite := by
  intro hfin
  exact hd (h.unique (hasDirichletDensity_zero_of_finite hfin))

/-! ### ★出典の紐付け -/

/-- ★★locator —— `Theorem 6.4, (iv)` が使う Tchebotarev の語彙(Dirichlet 密度)。 -/
def HasDirichletDensity.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 116,
    item := "Theorem 6.4, (iv) — Dirichlet 密度の定義と基本性質",
    sectionId := "frdi-thm-6-4" }

end ABC3.Found.NF
