/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GaloisRep.HtJBound
import Mathlib.NumberTheory.Height.NumberField
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★★★★`htJ` は mathlib の Weil 高さである（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

## ★★★★★★★★★★★★★★★★★★これは何か

`Found/GaloisRep/HtJBound.lean`（`§9-1002`）の `htJ L E` は

    `htJ L E = (1/d)·( Σ_p log⁺|j|_p·log N(p) + Σ_σ log⁺‖σ(j)‖ )`

と**素点ごとに手で書いた**ものである。★★本ファイルはそれが

    **`htJ L E = Height.logHeight₁ (E.j) / [L:ℚ]`**

——すなわち **mathlib の（絶対対数）Weil 高さそのもの**であることを示す。

## ★なぜ要るのか —— Northcott への接続

`Interface/GenEll/EllModuli.lean` の `northcott`（原文 `Proposition 1.4, (iv)`）は
「`ht^Falt` が有界な点は有限個」である。★`Found/GenEll/FiniteFromNorthcott.lean` の
`finite_of_le_of_northcott` はそれを `ht∞` の Northcott 性に帰す。

★★★`ht∞` を `htJ` と同定できれば、必要なのは
**`Height.logHeight₁` の Northcott 性**だけになる——
これは mathlib の `Mathlib/NumberTheory/Height/Northcott.lean` が
`Northcott (mulHeight₁)` から出す形で用意している枠組みそのものである
（☆数体でのインスタンスは 2026-08-29 時点で mathlib に無い。
本プロジェクトの `Found/GenEll/NorthcottImage.lean` から供給できる）。

## ★機構——素点ごとの照合

| 側 | 本プロジェクト | mathlib |
|---|---|---|
| アルキメデス | `Σ_{σ : L →+* ℂ} max(log‖σ j‖, 0)` | `Σ_{v : InfinitePlace} mult(v)·log⁺(v j)` |
| 有限 | `Σᶠ_p max(0, −v_p(j))·log N(p)` | `Σᶠ_{v : FinitePlace} log⁺(v j)` |

★アルキメデス側は `InfinitePlace.card_filter_mk_eq`（素点 `w` の上の埋め込みの個数が
`mult w`）でファイバーごとに足し合わせる。
★★有限側は `FinitePlace.equivHeightOneSpectrum` で添字を移し、
`FinitePlace.mk p x = N(p)^{−v_p(x)}` を対数に通す。
-/

namespace ABC3.Found.GaloisRep

open NumberField IsDedekindDomain Real WeierstrassCurve
open scoped Classical

variable {L : Type} [Field L] [NumberField L]

/-! ## ★★★★★有限素点側の照合 -/

/-- ★★★★★**`log⁺(|j|_p) = max(0, −v_p(j))·log N(p)`**。

★`FinitePlace.mk p x = N(p)^{−v_p(x)}`（`FinitePlace.norm_embedding'` と
`WithZeroMulInt.toNNReal_neg_apply`）を対数に通しただけである。
★★`j = 0` の場合は両辺 `0`（`Real.log 0 = 0`)。 -/
theorem posLog_finitePlace_eq (p : HeightOneSpectrum (𝓞 L))
    (W : WeierstrassCurve L) [W.IsElliptic] :
    (FinitePlace.mk p W.j : ℝ).posLog
      = ((max 0 (-jExp p W) : ℤ) : ℝ) * Real.log (Ideal.absNorm p.asIdeal) := by
  have hlog : (0:ℝ) ≤ Real.log (Ideal.absNorm p.asIdeal) := log_absNorm_nonneg p
  by_cases hj : W.j = 0
  · rw [jExp, dif_pos hj, hj]
    simp [Real.posLog_apply]
  · rw [jExp, dif_neg hj]
    have hne : (p.valuation L) W.j ≠ 0 := valuationP_ne_zero p (Units.mk0 W.j hj)
    have hval : (FinitePlace.mk p W.j : ℝ)
        = (Ideal.absNorm p.asIdeal : ℝ) ^ (-(valAdd p (Units.mk0 W.j hj))) := by
      rw [FinitePlace.mk_apply, FinitePlace.norm_embedding',
        WithZeroMulInt.toNNReal_neg_apply _ hne, valAdd, neg_neg]
      push_cast
      rfl
    rw [Real.posLog_apply, hval, Real.log_zpow]
    rcases le_or_gt (0:ℤ) (-(valAdd p (Units.mk0 W.j hj))) with h | h
    · rw [max_eq_right h, max_eq_right (by exact_mod_cast mul_nonneg (by exact_mod_cast h) hlog)]
    · rw [max_eq_left h.le, max_eq_left]
      · push_cast; ring
      · exact mul_nonpos_of_nonpos_of_nonneg (by exact_mod_cast h.le) hlog

/-- ★★★★★★**有限素点の和が一致する**——添字を `FinitePlace ≃ HeightOneSpectrum` で移す。 -/
theorem finsum_finite_eq (W : WeierstrassCurve L) [W.IsElliptic] :
    (∑ᶠ p : HeightOneSpectrum (𝓞 L),
        ((max 0 (-jExp p W) : ℤ) : ℝ) * Real.log (Ideal.absNorm p.asIdeal))
      = ∑ᶠ v : FinitePlace L, (v W.j).posLog := by
  rw [← finsum_comp_equiv (FinitePlace.equivHeightOneSpectrum (K := L))
    (f := fun p : HeightOneSpectrum (𝓞 L) =>
      ((max 0 (-jExp p W) : ℤ) : ℝ) * Real.log (Ideal.absNorm p.asIdeal))]
  refine finsum_congr (fun v => ?_)
  rw [← posLog_finitePlace_eq (FinitePlace.equivHeightOneSpectrum v) W,
    FinitePlace.equivHeightOneSpectrum_apply, FinitePlace.mk_maximalIdeal]

/-! ## ★★★★★★アルキメデス素点側の照合 -/

/-- ★★★★★★**埋め込みの和 ＝ 無限素点の重みつき和**。

★素点 `w` の上の埋め込みはちょうど `mult w` 個ある（`InfinitePlace.card_filter_mk_eq`）
——実素点で 1 個、複素素点で 2 個（互いに共役）。 -/
theorem sum_arch_eq (x : L) :
    ∑ σ : (L →+* ℂ), max (Real.log ‖σ x‖) 0
      = ∑ v : InfinitePlace L, (v.mult : ℝ) * (v x).posLog := by
  rw [← Finset.sum_fiberwise (g := fun σ : L →+* ℂ => InfinitePlace.mk σ)]
  refine Finset.sum_congr rfl (fun w _ => ?_)
  have h : ∀ σ ∈ Finset.univ.filter (fun σ : L →+* ℂ => InfinitePlace.mk σ = w),
      max (Real.log ‖σ x‖) 0 = (w x).posLog := by
    intro σ hσ
    simp only [Finset.mem_filter] at hσ
    rw [← hσ.2, Real.posLog_apply, InfinitePlace.apply, max_comm]
  rw [Finset.sum_congr rfl h, Finset.sum_const, nsmul_eq_mul]
  congr 1
  exact_mod_cast InfinitePlace.card_filter_mk_eq w

/-! ## ★★★★★★★★★★★★★★★★★★同定 -/

/-- ★★★★★★★★★★★★★★★★★★**`htJ` は `j` の絶対対数 Weil 高さである**。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

★`Height.logHeight₁` は mathlib の（相対）対数 Weil 高さ
——`[L:ℚ]` で割って絶対高さになる。 -/
theorem htJ_eq_logHeight (E : WeierstrassCurve L) [E.IsElliptic] :
    htJ L E = Height.logHeight₁ E.j / (Module.finrank ℚ L : ℝ) := by
  rw [htJ, htFinJ, htArchJ, NumberField.logHeight₁_eq, finsum_finite_eq E, sum_arch_eq]
  ring

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★**mathlib の言葉での到達点**:

    `Height.logHeight₁ (E.j) / [L:ℚ] ≤ 12(1+ϵ)·ht^Falt(E) + C`

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

★`C` は `L` にも `E` にも依らない普遍定数である（`§9-1000〜1003`）。
★★★これが [Silv2] `Proposition 2.1` の内容であり、
`Proposition 3.4` の第 3 の `≲` である。 -/
theorem exists_logHeight_j_le_htFalt (eps : ℝ) (heps : 0 < eps) :
    ∃ C : ℝ, ∀ (L : Type) [Field L] [NumberField L] (E : WeierstrassCurve L) [E.IsElliptic],
      Height.logHeight₁ E.j / (Module.finrank ℚ L : ℝ)
        ≤ 12 * (1 + eps) * htFaltOf L E + C := by
  obtain ⟨C, hC⟩ := exists_htJ_le_htFalt' eps heps
  refine ⟨C, fun L _ _ E _ => ?_⟩
  rw [← htJ_eq_logHeight E]
  exact hC L E

/-! ## ★出典の紐付け(`.src`) -/

def htJ_eq_logHeight.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(htJ は mathlib の Weil 高さ Height.logHeight₁ である)",
    sectionId := "genell-prop-3-4" }

def exists_logHeight_j_le_htFalt.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(第 3 の ≲ を mathlib の Height.logHeight₁ の言葉で——★無条件)",
    sectionId := "genell-prop-3-4" }

def htJ_eq_logHeight.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]"
      "NumberField.logHeight₁_eq(数体上の対数 Weil 高さの素点分解)"
      (.inMathlib "NumberField.logHeight₁_eq") 2,
    .citation "[mathlib]"
      "NumberField.InfinitePlace.card_filter_mk_eq(無限素点の上の埋め込みの個数 = mult)"
      (.inMathlib "NumberField.InfinitePlace.card_filter_mk_eq") 2,
    .citation "[mathlib]"
      "NumberField.FinitePlace.equivHeightOneSpectrum(有限素点と極大イデアルの対応)"
      (.inMathlib "NumberField.FinitePlace.equivHeightOneSpectrum") 2,
    .folklore
      ("☆**Northcott の定理**: Height.logHeight₁ が有界な点は次数を止めれば有限個。" ++
       "★mathlib は Northcott (mulHeight₁) → Northcott (logHeight₁) の枠組みだけ持ち、" ++
       "**数体でのインスタンスは無い**(2026-08-29 に infer_instance で確認)。" ++
       "★★本プロジェクトの Found/GenEll/NorthcottImage.lean の " ++
       "northcott_of_log_mulHeight_image から供給できる(未接続)") 7,
    .implicitStep
      ("★★★★★★★★到達点(2026-08-29、第 555): §9-1002 の htJ が" ++
       "**mathlib の Weil 高さそのもの**であることを示した。" ++
       "★これで Prop 3.4 の第 3 の ≲ は mathlib native な形 " ++
       "Height.logHeight₁ (E.j) / [L:ℚ] ≤ 12(1+ϵ)·ht^Falt + C で書ける。" ++
       "☆残るのは ht∞ の同定(M̄_ell の無限遠因子の高さが h(j) であること)と Northcott") 8 ]

end ABC3.Found.GaloisRep
