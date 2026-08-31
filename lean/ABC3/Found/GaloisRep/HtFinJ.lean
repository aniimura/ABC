/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GaloisRep.DegInf
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★★`h(j)` の有限素点側は `deg∞` で抑えられる（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.18。

原文 (GenEll p.18):
> First, observe that if v is any local height of EL, then d · deg∞([EL]) ≥

## ★★★★★★★★★★★★★★★★これは何か

`Found/GaloisRep/HtJBound.lean`（`§9-1002`、第 553）は

    `h(j) ≤ 12(1+ϵ)·ht^Falt + C`

を**有限素点側 1 本を仮定として**証明していた。★★本ファイルはその仮定

    **`h_fin(j) ≤ deg∞`**

を**証明する**——★★★しかも**半安定を仮定せずに**。

## ★機構——極小モデルの `c₄` が整だから

素点 `p` で極小なモデル `W'` を取ると `j` は不変で

    `j·Δ' = c₄'³`  すなわち  `v_p(j) = 3·v_p(c₄') − v_p(Δ_min)`

★極小モデルは**整**なので `v_p(c₄') ≥ 0`、したがって

    `−v_p(j) = v_p(Δ_min) − 3·v_p(c₄') ≤ v_p(Δ_min)`

★★`v_p(Δ_min) ≥ 0` と合わせて `log⁺|j|_p ≔ max(0, −v_p(j)) ≤ v_p(Δ_min)`。

☆半安定（乗法還元）のときは `v_p(c₄') = 0` なので**等号**になる
——原文が念頭に置く場合である。★★★不等号だけなら**半安定は要らない**。

## ★★得られるもの

| 定理 | 内容 |
|---|---|
| `jExp` | ★★★`j` の `p` での付値（`j = 0` では `0`） |
| `neg_jExp_le_minDeltaExp` | ★★★★★★**`−v_p(j) ≤ v_p(Δ_min)`** |
| `htFinJ` | ★★★★`h(j)` の有限素点側 |
| `htFinJ_le_degInfOf` | ★★★★★★★★**`h_fin(j) ≤ deg∞`（無条件）** |
-/

namespace ABC3.Found.GaloisRep

open IsDedekindDomain NumberField WeierstrassCurve
open scoped Classical

variable {L : Type} [Field L] [NumberField L]

/-! ## ★★★`j` の付値 -/

/-- ★★★**`j` の `p` での付値**——`j = 0` のときは `0` と定める。

★`j = 0`（`c₄ = 0`）は `ρ = e^{2πi/3}` の曲線であり、`log⁺|j|_p = 0` で正しい。 -/
noncomputable def jExp (p : HeightOneSpectrum (𝓞 L)) (W : WeierstrassCurve L) [W.IsElliptic] : ℤ :=
  if h : W.j = 0 then 0 else valAdd p (Units.mk0 W.j h)

/-- ★★★★★★**`−v_p(j) ≤ v_p(Δ_min)`**——極小モデルの `c₄` が整だから。

原文 (GenEll p.18):
> First, observe that if v is any local height of EL, then d · deg∞([EL]) ≥

★`j·Δ' = c₄'³` で `v_p(c₄') ≥ 0`。★★半安定は要らない
（半安定なら `v_p(c₄') = 0` で等号になる）。 -/
theorem neg_jExp_le_minDeltaExp (p : HeightOneSpectrum (𝓞 L)) (W : WeierstrassCurve L)
    [W.IsElliptic] : -jExp p W ≤ minDeltaExp p W := by
  by_cases hj : W.j = 0
  · rw [jExp, dif_pos hj, neg_zero]
    exact minDeltaExp_nonneg p W
  · rw [jExp, dif_neg hj]
    have hΔ : W.Δ ≠ 0 := W.isUnit_Δ.ne_zero
    obtain ⟨C, hC⟩ := WeierstrassCurve.exists_isMinimal (primeSubring p) W
    haveI : (C • W).IsElliptic := inferInstance
    have hjC : (C • W).j = W.j := WeierstrassCurve.variableChange_j W C
    have hΔC : (C • W).Δ ≠ 0 := variableChange_Delta_ne_zero W hΔ C
    have hmd := minDeltaExp_eq p W hΔ C hC
    have hmul : (C • W).j * (C • W).Δ = (C • W).c₄ ^ 3 := by
      rw [WeierstrassCurve.j, ← WeierstrassCurve.coe_Δ', Units.val_inv_eq_inv_val]
      have : ((C • W).Δ' : L) ≠ 0 := ((C • W).Δ').ne_zero
      field_simp
    have hc4 : (C • W).c₄ ≠ 0 := by
      intro h
      rw [h] at hmul
      simp only [zero_pow (by norm_num : (3:ℕ) ≠ 0)] at hmul
      rcases mul_eq_zero.1 hmul with h1 | h1
      · rw [hjC] at h1; exact hj h1
      · exact hΔC h1
    have hunits : Units.mk0 ((C • W).j) (by rw [hjC]; exact hj) * Units.mk0 ((C • W).Δ) hΔC
        = (Units.mk0 ((C • W).c₄) hc4) ^ 3 := by
      refine Units.ext ?_
      push_cast
      exact hmul
    have hval := congrArg (valAdd p) hunits
    rw [valAdd_mul, valAdd_pow] at hval
    have hc4nonneg : 0 ≤ valAdd p (Units.mk0 ((C • W).c₄) hc4) := by
      rw [valAdd_nonneg_iff]
      have hmem : ((C • W).c₄) ∈ primeSubring p := by
        have h1 := WeierstrassCurve.integralModel_c₄_eq (primeSubring p) (C • W)
        rw [← h1]
        exact ((integralModel (primeSubring p) (C • W)).c₄).2
      exact (mem_primeSubring_iff p _).1 hmem
    have hjeq : Units.mk0 W.j hj = Units.mk0 ((C • W).j) (by rw [hjC]; exact hj) := by
      refine Units.ext ?_
      simp [hjC]
    rw [hjeq, hmd]
    omega

/-- ★★★★★**`log⁺|j|_p ≤ v_p(Δ_min)`**。 -/
theorem maxJ_le_minDeltaExp (p : HeightOneSpectrum (𝓞 L)) (W : WeierstrassCurve L)
    [W.IsElliptic] : max 0 (-jExp p W) ≤ minDeltaExp p W :=
  max_le (minDeltaExp_nonneg p W) (neg_jExp_le_minDeltaExp p W)

/-! ## ★★★★有限素点側の高さ -/

/-- ★★★★**`j` の高さの有限素点側** `h_fin(j) = (1/d)·Σ_p log⁺|j|_p·log N(p)`。 -/
noncomputable def htFinJ (L : Type) [Field L] [NumberField L] (W : WeierstrassCurve L)
    [W.IsElliptic] : ℝ :=
  (∑ᶠ p : HeightOneSpectrum (𝓞 L),
      ((max 0 (-jExp p W) : ℤ) : ℝ) * Real.log (Ideal.absNorm p.asIdeal))
    / (Module.finrank ℚ L)

/-- ★★★★★★★★**`h_fin(j) ≤ deg∞`**——★**無条件**（半安定を仮定しない）。

原文 (GenEll p.18):
> First, observe that if v is any local height of EL, then d · deg∞([EL]) ≥

★項ごとに `max(0, −v_p(j)) ≤ v_p(Δ_min)`。
★★台の有限性は `minDeltaExp_finite` から従う
（`v_p(Δ_min) = 0` なら `log⁺|j|_p = 0` だから）。 -/
theorem htFinJ_le_degInfOf (W : WeierstrassCurve L) [W.IsElliptic] :
    htFinJ L W ≤ degInfOf L W := by
  have hΔ : W.Δ ≠ 0 := W.isUnit_Δ.ne_zero
  have hfin := minDeltaExp_finite W hΔ
  have hg : (Function.support (fun p : HeightOneSpectrum (𝓞 L) =>
      (minDeltaExp p W : ℝ) * Real.log (Ideal.absNorm p.asIdeal))).Finite := by
    refine hfin.subset (fun p hp => ?_)
    simp only [Function.mem_support, ne_eq, mul_eq_zero, not_or] at hp
    simp only [Set.mem_setOf_eq]
    exact_mod_cast hp.1
  have hf : (Function.support (fun p : HeightOneSpectrum (𝓞 L) =>
      ((max 0 (-jExp p W) : ℤ) : ℝ) * Real.log (Ideal.absNorm p.asIdeal))).Finite := by
    refine hfin.subset (fun p hp => ?_)
    simp only [Function.mem_support, ne_eq, mul_eq_zero, not_or] at hp
    simp only [Set.mem_setOf_eq]
    intro hz
    have h1 := maxJ_le_minDeltaExp p W
    rw [hz] at h1
    have h2 : max 0 (-jExp p W) = 0 := le_antisymm h1 (le_max_left _ _)
    exact hp.1 (by exact_mod_cast h2)
  have hS : (∑ᶠ p : HeightOneSpectrum (𝓞 L),
        ((max 0 (-jExp p W) : ℤ) : ℝ) * Real.log (Ideal.absNorm p.asIdeal))
      ≤ ∑ᶠ p : HeightOneSpectrum (𝓞 L),
        (minDeltaExp p W : ℝ) * Real.log (Ideal.absNorm p.asIdeal) := by
    refine finsum_le_finsum' hf hg (fun p => ?_)
    refine mul_le_mul_of_nonneg_right ?_ (log_absNorm_nonneg p)
    exact Int.cast_le.2 (maxJ_le_minDeltaExp p W)
  rw [htFinJ, degInfOf]
  have hd : (0:ℝ) ≤ (Module.finrank ℚ L : ℝ) := by positivity
  gcongr

/-- ★★★★**`h_fin(j) ≥ 0`**。 -/
theorem htFinJ_nonneg (W : WeierstrassCurve L) [W.IsElliptic] : 0 ≤ htFinJ L W := by
  rw [htFinJ]
  refine div_nonneg ?_ (by positivity)
  refine finsum_nonneg (fun p => ?_)
  exact mul_nonneg (by exact_mod_cast le_max_left 0 (-jExp p W)) (log_absNorm_nonneg p)

/-! ## ★出典の紐付け(`.src`) -/

def jExp.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 18,
    item := "Proposition 3.4(j の有限素点での付値)",
    sectionId := "genell-prop-3-4" }

def htFinJ_le_degInfOf.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 18,
    item := "Proposition 3.4(h(j) の有限素点側 ≤ deg∞——★無条件)",
    sectionId := "genell-prop-3-4" }

def htFinJ_le_degInfOf.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[mathlib]"
      "WeierstrassCurve.integralModel_c₄_eq(極小モデルの c₄ は整)"
      (.inMathlib "WeierstrassCurve.integralModel_c₄_eq") 2,
    .citation "[mathlib]"
      "WeierstrassCurve.variableChange_j(j は変数変換で不変)"
      (.inMathlib "WeierstrassCurve.variableChange_j") 1,
    .implicitStep
      ("★★★★★★★★到達点(2026-08-29、第 554): " ++
       "Found/GaloisRep/HtJBound.lean が仮定として受けていた h_fin(j) ≤ deg∞ を" ++
       "**証明した**。★しかも**半安定を仮定せずに**——" ++
       "極小モデルの c₄ が整であること(v_p(c₄') ≥ 0)だけで不等号は出る。" ++
       "☆半安定なら v_p(c₄') = 0 で等号になる(原文が念頭に置く場合)。" ++
       "★★これで h(j) ≤ 12(1+ϵ)ht^Falt + C は**無条件**になる") 9 ]

end ABC3.Found.GaloisRep
