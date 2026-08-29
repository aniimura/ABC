/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.JPeterssonBound
import ABC3.Found.GenEll.ArchNormTotal
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★★★★★★★★アルキメデス素点での `log⁺|j|`（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

## ★★★★★★★★★★★★★★★★★★★★★★これは何か

`§9-1000`（第 551）で **`‖j(τ)‖·(‖Δ‖_Pet(τ))^{1+ϵ}` は有界**が取れた。
★★本ファイルはそれを**曲線の言葉に翻訳する**——`τ` を消して

    **`log⁺‖j(W)‖ + (1+ϵ)·log((2π)¹²·curveArchInv W) ≤ C`**

にする。★★★左の第 2 項は **`archSum` の 1 項そのもの**であり
（`Found/GaloisRep/FaltingsWitness.lean` の `archSum L E = ∑_σ log((2π)¹²·archNorm E σ)`）、
`σ` について足して `[L:ℚ]` で割れば

    `h∞(j) ≤ C − (1+ϵ)·archSum/[L:ℚ]`

になる。★★★★`12·ht^Falt = deg∞ − archSum/[L:ℚ]` を代入すると

    `h∞(j) ≤ C + 12(1+ϵ)·ht^Falt − (1+ϵ)·deg∞`

——★★★★★これが `Prop 3.4` の**第 3 の `≲`** の中身である。

## ★機構

`curveArchInv W = 4096π¹²·peterssonDelta τ`（`jFun τ = W.j` なる `τ`、`§9-354`）なので

    `(2π)¹²·curveArchInv W = A·peterssonDelta τ`,  `A ≔ (2π)¹²·4096π¹²`

★あとは `‖jFun τ‖·pD^{1+ϵ} ≤ M` の対数を取るだけである。
★★`log⁺` の `0` 側は `pD` の一様上界（`§9-355`）で処理する。
-/

namespace ABC3.Found.GenEll

open Real Complex UpperHalfPlane

/-! ## ★★★★★★★★★★★★★★★★★★★★★★曲線の言葉に翻訳する -/

/-- ★★★★★★★★★★★★★★★★★★★★★★**`log⁺‖j(W)‖ + (1+ϵ)·log((2π)¹²·curveArchInv W)` は有界**。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

★`§9-1000` の `exists_bound_jFun_peterssonDelta` を曲線の言葉に翻訳したもの。
★★左の第 2 項は `archSum` の 1 項そのものである。 -/
theorem exists_logPlus_j_bound (eps : ℝ) (heps : 0 < eps) :
    ∃ C : ℝ, ∀ (W : WeierstrassCurve ℂ) [W.IsElliptic],
      max (Real.log ‖W.j‖) 0
        + (1 + eps) * Real.log ((2 * Real.pi) ^ 12 * curveArchInv W) ≤ C := by
  obtain ⟨M, hMpos, hM⟩ := exists_bound_jFun_peterssonDelta eps heps
  obtain ⟨M₀, hM₀⟩ := exists_bound_peterssonDelta
  set M₁ : ℝ := max M₀ 1 with hM₁def
  have hM₁pos : (0:ℝ) < M₁ := lt_of_lt_of_le zero_lt_one (le_max_right _ _)
  have hM₁ : ∀ τ : UpperHalfPlane, peterssonDelta τ ≤ M₁ :=
    fun τ => le_trans (hM₀ τ) (le_max_left _ _)
  set A : ℝ := (2 * Real.pi) ^ 12 * (4096 * Real.pi ^ 12) with hAdef
  have hApos : (0:ℝ) < A := by
    have := Real.pi_pos
    rw [hAdef]; positivity
  have h1e : (0:ℝ) < 1 + eps := by linarith
  refine ⟨max (Real.log M + (1 + eps) * Real.log A) ((1 + eps) * Real.log (A * M₁)),
    fun W hW => ?_⟩
  obtain ⟨τ, hjτ, hci⟩ := exists_curveArchInv_eq_petersson W hW
  have hpd : (0:ℝ) < peterssonDelta τ := peterssonDelta_pos τ
  have hkey : (2 * Real.pi) ^ 12 * curveArchInv W = A * peterssonDelta τ := by
    rw [hci, hAdef]; ring
  rw [hkey]
  have hlogA : Real.log (A * peterssonDelta τ) = Real.log A + Real.log (peterssonDelta τ) :=
    Real.log_mul (ne_of_gt hApos) (ne_of_gt hpd)
  rcases le_or_gt (Real.log ‖W.j‖) 0 with hj | hj
  · rw [max_eq_right hj, zero_add]
    refine le_trans ?_ (le_max_right _ _)
    refine mul_le_mul_of_nonneg_left ?_ h1e.le
    exact Real.log_le_log (by positivity) (by nlinarith [hM₁ τ, hApos.le])
  · rw [max_eq_left hj.le]
    refine le_trans ?_ (le_max_left _ _)
    have hjn : ‖W.j‖ = ‖jFun τ‖ := by rw [hjτ]
    have hjpos : (0:ℝ) < ‖jFun τ‖ := by
      rcases (norm_nonneg (jFun τ)).lt_or_eq with h | h
      · exact h
      · rw [hjn, ← h] at hj; simp [Real.log_zero] at hj
    have hlog := Real.log_le_log (by positivity) (hM τ)
    rw [Real.log_mul (ne_of_gt hjpos) (by positivity), Real.log_rpow hpd] at hlog
    rw [hlogA, hjn]
    linarith

/-- ★★★★★★**埋め込みごとの形**——`archNorm` の言葉で書き直したもの。

★`(E.map σ).j = σ (E.j)`（`WeierstrassCurve.map_j`）と
`archNorm E σ = curveArchInv (E.map σ)` を使うだけである。 -/
theorem logPlus_archNorm_le {L : Type*} [Field L] (E : WeierstrassCurve L) [E.IsElliptic]
    (σ : L →+* ℂ) (eps C : ℝ)
    (hC : ∀ (W : WeierstrassCurve ℂ) [W.IsElliptic],
      max (Real.log ‖W.j‖) 0
        + (1 + eps) * Real.log ((2 * Real.pi) ^ 12 * curveArchInv W) ≤ C) :
    max (Real.log ‖σ E.j‖) 0
      + (1 + eps) * Real.log ((2 * Real.pi) ^ 12 * archNorm E σ) ≤ C := by
  simpa [archNorm, WeierstrassCurve.map_j] using hC (E.map σ)

/-! ## ★出典の紐付け(`.src`) -/

def exists_logPlus_j_bound.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(第 3 の ≲ ——log⁺|j| ＋ (1+ϵ)log((2π)¹²‖Δ‖_arch) は有界)",
    sectionId := "genell-prop-3-4" }

def logPlus_archNorm_le.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(第 3 の ≲ を archNorm の言葉で)",
    sectionId := "genell-prop-3-4" }

def exists_logPlus_j_bound.needs : List ABC3.Meta.ProofObligation :=
  [ .otherPaper "[Silv2]"
      ("Proposition 2.1——★原文が Prop 3.4 の根拠として引く。" ++
       "★★**そのアルキメデス側は本ファイルで閉じた**") 4,
    .implicitStep
      ("★★★★★★★★到達点(2026-08-29、第 552): §9-1000 の " ++
       "exists_bound_jFun_peterssonDelta を**曲線の言葉に翻訳した**。" ++
       "★curveArchInv W = 4096π¹²·peterssonDelta τ(jFun τ = W.j なる τ、§9-354)により " ++
       "τ が消える。★★左の第 2 項は archSum の 1 項そのものなので、" ++
       "σ について足して [L:ℚ] で割れば h∞(j) ≤ C − (1+ϵ)·archSum/[L:ℚ] になる") 8,
    .implicitStep
      ("☆残るのは**有限素点側**: 半安定のとき h(j) の有限素点の寄与 ≤ deg∞ " ++
       "(v(j) = v(c₄³) − v(Δ_min) で、乗法還元なら v(c₄) = 0 なので " ++
       "v(j) = −v(Δ_min) ≤ 0)。★これは代数の段である") 7 ]

end ABC3.Found.GenEll
