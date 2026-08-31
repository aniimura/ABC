/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GaloisRep.HtJWeil
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★対数版の積公式（`Found`、無条件）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.18。

原文 (GenEll p.18):
> First, observe that if v is any local height of EL, then d · deg∞([EL]) ≥

## ★★★★★★★★★★★★これは何か

    **`Σ_{σ : L →+* ℂ} log‖σ(x)‖ = Σᶠ_p v_p(x)·log N(p)`**

★数体の積公式（mathlib の `NumberField.prod_abs_eq_one`）を、
**本プロジェクトの `valAdd`（`p` での付値）の言葉で**対数の形に書いたものである。

## ★なぜ要るのか —— `deg∞` と `archSum` を繋ぐ唯一の道具

`12·htFaltOf = deg∞ − archSum/d`（`§9-670`）で、

* `deg∞` は**有限素点**の `minDeltaExp = v_p(Δ) − 12·neronExp` から、
* `archSum` は**アルキメデス素点**の `curveArchInv = |Δ|·covol⁶` から

作られている。★★両者を突き合わせるには `Σ_σ log|σΔ|` と `Σ_p v_p(Δ)·log N(p)` を
繋がねばならず、それがこの公式である。

★★★`§9-1020`（第 576）の測定はこれを紙の上で使って

    `12·htFaltOf = −12log(2π) − (12/d)Σ_p neronExp_p·log N(p) − (6/d)Σ_σ log covol(Λ_{E^σ})`

を得た。★本ファイルはその**道具の側を Lean に載せる**。

## ★★機構

`NumberField.prod_abs_eq_one`:「`(∏_w (w x)^{mult w})·(∏ᶠ_w w x) = 1`」の対数を取る。

* 無限素点側 `Σ_w mult(w)·log(w x)` は、素点の上の埋め込みが `mult w` 個あること
  （`InfinitePlace.card_filter_mk_eq`）から `Σ_σ log‖σ x‖` になる。
* 有限素点側 `Σᶠ_w log(w x)` は、`FinitePlace.mk p x = N(p)^{−v_p(x)}` から
  `Σᶠ_p (−v_p(x))·log N(p)` になる（`FinitePlace.equivHeightOneSpectrum` で添字を移す）。

★★どちらも `Found/GaloisRep/HtJWeil.lean`（`§9-1004`）で `log⁺` 版を作ったときの
機構をそのまま `log` に取り替えたものである。
-/

namespace ABC3.Found.GaloisRep

open NumberField IsDedekindDomain WeierstrassCurve ABC3.Found.GenEll
open scoped Classical

variable {L : Type} [Field L] [NumberField L]

/-! ## ★★★★★アルキメデス側 -/

/-- ★★★★★**埋め込みの和 ＝ 無限素点の重みつき和**（`log` 版）。

★`§9-1004` の `sum_arch_eq`（`log⁺` 版）の `log` 版である。 -/
theorem sum_arch_log_eq (x : L) :
    ∑ σ : (L →+* ℂ), Real.log ‖σ x‖
      = ∑ w : InfinitePlace L, (w.mult : ℝ) * Real.log (w x) := by
  rw [← Finset.sum_fiberwise (g := fun σ : L →+* ℂ => InfinitePlace.mk σ)]
  refine Finset.sum_congr rfl (fun w _ => ?_)
  have h : ∀ σ ∈ Finset.univ.filter (fun σ : L →+* ℂ => InfinitePlace.mk σ = w),
      Real.log ‖σ x‖ = Real.log (w x) := by
    intro σ hσ
    simp only [Finset.mem_filter] at hσ
    rw [← hσ.2, InfinitePlace.apply]
  rw [Finset.sum_congr rfl h, Finset.sum_const, nsmul_eq_mul]
  congr 1
  exact_mod_cast InfinitePlace.card_filter_mk_eq w

/-! ## ★★★★★有限素点側 -/

/-- ★★★★★**`log|x|_p = −v_p(x)·log N(p)`**。

★`§9-1004` の `posLog_finitePlace_eq` の `log` 版（`log⁺` を外したもの）。 -/
theorem log_finitePlace_eq (p : HeightOneSpectrum (𝓞 L)) (x : L) (hx : x ≠ 0) :
    Real.log (FinitePlace.mk p x : ℝ)
      = -(valAdd p (Units.mk0 x hx) : ℝ) * Real.log (Ideal.absNorm p.asIdeal) := by
  have hne : (p.valuation L) x ≠ 0 := valuationP_ne_zero p (Units.mk0 x hx)
  have hval : (FinitePlace.mk p x : ℝ)
      = (Ideal.absNorm p.asIdeal : ℝ) ^ (-(valAdd p (Units.mk0 x hx))) := by
    rw [FinitePlace.mk_apply, FinitePlace.norm_embedding',
      WithZeroMulInt.toNNReal_neg_apply _ hne, valAdd, neg_neg]
    push_cast
    rfl
  rw [hval, Real.log_zpow]
  push_cast
  ring

/-- ★★★★★**有限素点の和を `HeightOneSpectrum` の和に移す**。 -/
theorem finsum_log_finite_eq (x : L) (hx : x ≠ 0) :
    (∑ᶠ w : FinitePlace L, Real.log (w x))
      = ∑ᶠ p : HeightOneSpectrum (𝓞 L),
          -(valAdd p (Units.mk0 x hx) : ℝ) * Real.log (Ideal.absNorm p.asIdeal) := by
  rw [← finsum_comp_equiv (FinitePlace.equivHeightOneSpectrum (K := L))
    (f := fun p : HeightOneSpectrum (𝓞 L) =>
      -(valAdd p (Units.mk0 x hx) : ℝ) * Real.log (Ideal.absNorm p.asIdeal))]
  refine finsum_congr (fun v => ?_)
  rw [← log_finitePlace_eq (FinitePlace.equivHeightOneSpectrum v) x hx,
    FinitePlace.equivHeightOneSpectrum_apply, FinitePlace.mk_maximalIdeal]

/-! ## ★★★★★★★★★★積公式 -/

/-- ★★★★★★★**積公式の対数版**（素点の言葉のまま）。

★mathlib の `NumberField.prod_abs_eq_one` の対数を取っただけである。 -/
theorem log_product_formula (x : L) (hx : x ≠ 0) :
    (∑ w : InfinitePlace L, (w.mult : ℝ) * Real.log (w x))
      + (∑ᶠ w : FinitePlace L, Real.log (w x)) = 0 := by
  have hprod := NumberField.prod_abs_eq_one hx
  have hinf : ∀ w : InfinitePlace L, (0:ℝ) < w x := fun w => w.pos_iff.2 hx
  have hfin : ∀ w : FinitePlace L, (0:ℝ) < w x := fun w => FinitePlace.pos_iff.2 hx
  have hA : (0:ℝ) < ∏ w : InfinitePlace L, (w x) ^ w.mult :=
    Finset.prod_pos (fun w _ => pow_pos (hinf w) _)
  have hB : (0:ℝ) < ∏ᶠ w : FinitePlace L, (w x) := by
    rw [NumberField.FinitePlace.prod_eq_inv_abs_norm hx]
    have hn : (Algebra.norm ℚ) x ≠ 0 := by
      simpa using (Algebra.norm_ne_zero_iff (R := ℚ) (S := L)).2 hx
    positivity
  have hlog := congrArg Real.log hprod
  rw [Real.log_mul (ne_of_gt hA) (ne_of_gt hB), Real.log_one] at hlog
  rw [Real.log_prod (fun w _ => ne_of_gt (pow_pos (hinf w) _))] at hlog
  simp only [Real.log_pow] at hlog
  rw [Real.log_finprod hfin] at hlog
  exact hlog

/-- ★★★★★★★★★★★★**`Σ_σ log‖σ(x)‖ = Σᶠ_p v_p(x)·log N(p)`**——★**無条件**。

原文 (GenEll p.18):
> First, observe that if v is any local height of EL, then d · deg∞([EL]) ≥

★★これが `deg∞`（有限素点）と `archSum`（アルキメデス素点）を繋ぐ道具である。
★★★`§9-1020`（第 576）が紙の上で使った関係の、Lean 上の実体。 -/
theorem sum_arch_log_eq_finsum_valAdd (x : L) (hx : x ≠ 0) :
    ∑ σ : (L →+* ℂ), Real.log ‖σ x‖
      = ∑ᶠ p : HeightOneSpectrum (𝓞 L),
          (valAdd p (Units.mk0 x hx) : ℝ) * Real.log (Ideal.absNorm p.asIdeal) := by
  have h1 := log_product_formula x hx
  rw [finsum_log_finite_eq x hx] at h1
  rw [sum_arch_log_eq x]
  have h2 : (∑ᶠ p : HeightOneSpectrum (𝓞 L),
      -(valAdd p (Units.mk0 x hx) : ℝ) * Real.log (Ideal.absNorm p.asIdeal))
      = -∑ᶠ p : HeightOneSpectrum (𝓞 L),
          (valAdd p (Units.mk0 x hx) : ℝ) * Real.log (Ideal.absNorm p.asIdeal) := by
    rw [← finsum_neg_distrib]
    refine finsum_congr (fun p => ?_)
    ring
  rw [h2] at h1
  linarith

/-! ## ★出典の紐付け(`.src`) -/

def sum_arch_log_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 18,
    item := "Proposition 3.4(埋め込みの和 = 無限素点の重みつき和。log 版)",
    sectionId := "genell-prop-3-4" }

def sum_arch_log_eq_finsum_valAdd.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 18,
    item := "Proposition 3.4(対数版の積公式——deg∞ と archSum を繋ぐ道具。★無条件)",
    sectionId := "genell-prop-3-4" }

def sum_arch_log_eq_finsum_valAdd.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]" "NumberField.prod_abs_eq_one(数体の積公式)"
      (.inMathlib "NumberField.prod_abs_eq_one") 2,
    .citation "[mathlib]"
      "NumberField.FinitePlace.prod_eq_inv_abs_norm(有限素点の積はノルムの逆数)"
      (.inMathlib "NumberField.FinitePlace.prod_eq_inv_abs_norm") 2,
    .citation "[mathlib]"
      "NumberField.InfinitePlace.card_filter_mk_eq(素点の上の埋め込みは mult 個)"
      (.inMathlib "NumberField.InfinitePlace.card_filter_mk_eq") 2,
    .implicitStep
      ("★★★★★★★★到達点(2026-08-29、第 577): deg∞(有限素点)と archSum" ++
       "(アルキメデス素点)を繋ぐ道具を Lean に載せた。" ++
       "★§9-1020(第 576)の測定" ++
       "『12·htFaltOf = −12log(2π) − (12/d)Σ_p neronExp_p·log N(p) " ++
       "− (6/d)Σ_σ log covol(Λ_{E^σ})』は紙の上でこの公式を使っている。" ++
       "★★機構は §9-1004(HtJWeil.lean)の log⁺ 版をそのまま log に取り替えたもの") 8,
    .implicitStep
      ("☆次の一手: 本公式と curveArchInv = |Δ|·covol⁶(§9-353)を合わせて、" ++
       "上の htFaltOf の共体積表示を Lean で証明すること。" ++
       "★そうすれば同種写像評価は『各素点で neronExp_p(E/H) − neronExp_p(E) が " ++
       "log(l) 程度』に完全に局所化される") 9 ]

end ABC3.Found.GaloisRep
