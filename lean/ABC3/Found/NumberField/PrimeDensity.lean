/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.NumberField.LogZeta
import ABC3.Found.NumberField.DirichletDensity

/-!
# `Σ_p a_K(p)·p^{-s} ~ log(1/(s−1))`(鎖 `cheb` の `cheb-split-density` の解析側)

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.116。

原文 (FrdI p.116):
> over a prime pi ∈V(Li), then pi splits completely in Li if and only if deg(Li, vi) =

## ★★何をするか

`LogZeta.lean` の `|log ζ_L(s) − Σ_p a_L(p)p^{-s}| ≤ C` と、
mathlib の**類数公式**(`tendsto_sub_one_mul_dedekindZeta_nhdsGT`)を合わせて

  `Σ_p a_L(p)·p^{-s} / log(1/(s−1)) → 1`  (`s → 1+`)

を出す。★これが「完全分解する素点の密度は `1/[L:K]`」の**解析側の全部**である。
残るのは `a_L(p)` を分解の様子で読む**代数側**だけになる。

## ★本ファイルで閉じること

| 定理 | 中身 |
|---|---|
| `tendsto_sub_one_mul_zetaR` | 実の側での単純極(mathlib の類数公式の実部) |
| `tendsto_log_zetaR_sub` | ★`log ζ_L(s) − log(1/(s−1)) → log(留数)` |
| `tendsto_tsum_primes_div_log` | ★★★★★★**`Σ_p a_L(p)p^{-s}/log(1/(s−1)) → 1`** |

★★中身は `A/D = 1 + (Z−D)/D + (A−Z)/D` の 3 項分解である ——
第 2 項は分子が収束して分母が `+∞`、第 3 項は分子が**一様有界**(`logConst`)。
-/

namespace ABC3.Found.NF

open NumberField Filter Topology

variable (K : Type*) [Field K] [NumberField K]

/-- ★実の側での単純極 —— mathlib の類数公式の実部を取るだけ。 -/
theorem tendsto_sub_one_mul_zetaR :
    Tendsto (fun s : ℝ ↦ (s - 1) * zetaR K s) (𝓝[>] 1)
      (𝓝 (dedekindZeta_residue K)) := by
  have h : Tendsto (fun s : ℝ ↦ (((s : ℂ) - 1) * dedekindZeta K (s : ℂ)).re) (𝓝[>] 1)
      (𝓝 ((((dedekindZeta_residue K : ℝ)) : ℂ)).re) :=
    (Complex.continuous_re.tendsto _).comp (tendsto_sub_one_mul_dedekindZeta_nhdsGT K)
  simp only [Complex.ofReal_re] at h
  refine h.congr' ?_
  filter_upwards [self_mem_nhdsWithin] with s hs
  have hs1 : (1 : ℝ) < s := hs
  have heq : ((s : ℂ) - 1) * dedekindZeta K (s : ℂ)
      = (((s - 1) * zetaR K s : ℝ) : ℂ) := by
    rw [Complex.ofReal_mul, ofReal_zetaR K hs1]
    push_cast
    ring
  rw [heq, Complex.ofReal_re]

/-- ★★**`log ζ_L(s) − log(1/(s−1)) → log(留数)`**。 -/
theorem tendsto_log_zetaR_sub :
    Tendsto (fun s : ℝ ↦ Real.log (zetaR K s) - Real.log (1 / (s - 1))) (𝓝[>] 1)
      (𝓝 (Real.log (dedekindZeta_residue K))) := by
  have hres := dedekindZeta_residue_pos K
  have h := (Real.continuousAt_log (ne_of_gt hres)).tendsto.comp
    (tendsto_sub_one_mul_zetaR K)
  refine h.congr' ?_
  filter_upwards [self_mem_nhdsWithin] with s hs
  have hs1 : (1 : ℝ) < s := hs
  have hsub : (0 : ℝ) < s - 1 := by linarith
  have hz : (0 : ℝ) < zetaR K s := lt_of_lt_of_le zero_lt_one (one_le_zetaR K hs1)
  show Real.log ((s - 1) * zetaR K s) = _
  rw [Real.log_mul (ne_of_gt hsub) (ne_of_gt hz), one_div, Real.log_inv]
  ring

/-- ★★★★★★**`Σ_p a_K(p)·p^{-s} / log(1/(s−1)) → 1`**。

原文 (FrdI p.116):
> over a prime pi ∈V(Li), then pi splits completely in Li if and only if deg(Li, vi) =

★★中身は `A/D = 1 + (Z−D)/D + (A−Z)/D` の 3 項分解 ——
第 2 項は `tendsto_log_zetaR_sub`(分子が収束・分母が `+∞`)、
第 3 項は `abs_log_zetaR_sub_tsum_le`(分子が**一様有界**)。 -/
theorem tendsto_tsum_primes_div_log :
    Tendsto (fun s : ℝ ↦
        (∑' p : Nat.Primes, zetaSummandR K s p) / Real.log (1 / (s - 1)))
      (𝓝[>] 1) (𝓝 1) := by
  set D : ℝ → ℝ := fun s => Real.log (1 / (s - 1)) with hD
  set Z : ℝ → ℝ := fun s => Real.log (zetaR K s) with hZ
  set A : ℝ → ℝ := fun s => ∑' p : Nat.Primes, zetaSummandR K s p with hA
  have hDtop : Tendsto D (𝓝[>] (1:ℝ)) atTop := tendsto_logInv_atTop
  have h1 : Tendsto (fun s => (Z s - D s) / D s) (𝓝[>] (1:ℝ)) (𝓝 0) :=
    (tendsto_log_zetaR_sub K).div_atTop hDtop
  have hconst : Tendsto (fun _ : ℝ => logConst K) (𝓝[>] (1:ℝ)) (𝓝 (logConst K)) :=
    tendsto_const_nhds
  have h2 : Tendsto (fun s => (A s - Z s) / D s) (𝓝[>] (1:ℝ)) (𝓝 0) := by
    refine squeeze_zero_norm' ?_ (hconst.div_atTop hDtop)
    filter_upwards [self_mem_nhdsWithin, hDtop.eventually_gt_atTop 0] with s hs hDs
    have hs1 : (1 : ℝ) < s := hs
    have hbound := abs_log_zetaR_sub_tsum_le K hs1
    rw [Real.norm_eq_abs, abs_div, abs_of_pos hDs]
    refine div_le_div_of_nonneg_right ?_ (le_of_lt hDs)
    rw [abs_sub_comm]
    exact hbound
  have hcomb : Tendsto (fun s => 1 + ((Z s - D s) / D s + (A s - Z s) / D s))
      (𝓝[>] (1:ℝ)) (𝓝 (1 + (0 + 0))) := tendsto_const_nhds.add (h1.add h2)
  rw [show (1 : ℝ) + (0 + 0) = 1 by ring] at hcomb
  refine hcomb.congr' ?_
  filter_upwards [hDtop.eventually_gt_atTop 0] with s hDs
  have hne : D s ≠ 0 := ne_of_gt hDs
  show 1 + ((Z s - D s) / D s + (A s - Z s) / D s) = A s / D s
  field_simp
  ring

/-! ### ★出典の紐付け -/

/-- ★★★★★★locator —— `Theorem 6.4, (iv)` が使う Tchebotarev の解析側の全部。 -/
def tendsto_tsum_primes_div_log.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 116,
    item := "Theorem 6.4, (iv) — Σ_p a(p)p^{-s} ~ log(1/(s−1))",
    sectionId := "frdi-thm-6-4" }

def tendsto_tsum_primes_div_log.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "abs_log_zetaR_sub_tsum_le(log ζ の一様評価)"
      (.inProject "ABC3" "ABC3.Found.NF.abs_log_zetaR_sub_tsum_le") 116,
    .citation "[mathlib]" "tendsto_sub_one_mul_dedekindZeta_nhdsGT(類数公式)"
      (.inMathlib "NumberField.tendsto_sub_one_mul_dedekindZeta_nhdsGT") 116 ]

end ABC3.Found.NF
