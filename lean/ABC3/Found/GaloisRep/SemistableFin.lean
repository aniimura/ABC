/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GaloisRep.HtJLower
import ABC3.Meta.Claim

/-!
# ★★★★★★★★★★★★★★★★★★★★★★★★★★半安定なら `deg∞ = h_fin(j)`（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.18。

原文 (GenEll p.18):
> First, observe that if v is any local height of EL, then d · deg∞([EL]) ≥

## ★★★★★★★★★★★★★★★★★★★★★★★★★★これは何か

`Found/GaloisRep/HtJLower.lean` の `prop_3_4_chain` は
第 1・第 3 の `≲` を仮定 `hss : deg∞ ≤ h_fin(j)` のもとで出していた。
★★本ファイルはその仮定を、**素点ごとの半安定性**から**証明する**。

    素点 `p` で半安定  ⟹  `v_p(Δ_min) ≤ log⁺|j|_p`

★逆向き（`log⁺|j|_p ≤ v_p(Δ_min)`）は `Found/GaloisRep/HtFinJ.lean` で
**無条件**に取れているので、半安定なら**等号**である。

## ★機構

極小モデル `W'` で `j·Δ' = c₄'³`、すなわち `v_p(j) = 3v_p(c₄') − v_p(Δ_min)`。

* **良還元**（`v_p(Δ_min) = 0`）: 左辺が `0` なので不等式は自明。
* **乗法還元**（`v_p(c₄') = 0`）: `v_p(j) = −v_p(Δ_min) ≤ 0` で、
  `log⁺|j|_p = v_p(Δ_min)` ——★**等号**。
* ☆**加法還元**では `v_p(c₄') > 0` になり、破れる（`3v_p(c₄')` のぶんだけ `deg∞` が大きい）。

## ★★`SemistableAt` の書き方 —— 本プロジェクトの付値の言葉で

mathlib の `WeierstrassCurve.HasMultiplicativeReduction`
（`AlgebraicGeometry/EllipticCurve/Reduction.lean`）は
`valuation K (maximalIdeal R) W.c₄ = 1` で乗法還元を定める。
★これは本ファイルの第 2 の選言肢と**同じ内容**だが、
`primeSubring p` の DVR 付値と `p.valuation L` を繋ぐ橋がまだ無い。
★★そこで本ファイルは**本プロジェクトの `valAdd` の言葉で**書く
——☆その橋渡しは `.needs` に残す。
-/

namespace ABC3.Found.GaloisRep

open IsDedekindDomain NumberField WeierstrassCurve
open scoped Classical

variable {L : Type} [Field L] [NumberField L]

/-! ## ★★★★★★素点ごとの半安定性 -/

/-- ★★★★★★**素点 `p` で半安定**——良還元（`v_p(Δ_min) = 0`）または
乗法還元（極小モデルの `c₄` が単元）。

原文 (GenEll p.18):
> First, observe that if v is any local height of EL, then d · deg∞([EL]) ≥

★mathlib の `HasGoodReduction`／`HasMultiplicativeReduction` と同じ内容を、
本プロジェクトの `valAdd` の言葉で書いたものである。 -/
def SemistableAt (p : HeightOneSpectrum (𝓞 L)) (W : WeierstrassCurve L) : Prop :=
  minDeltaExp p W = 0 ∨
    ∃ C : WeierstrassCurve.VariableChange L, IsMinimal (primeSubring p) (C • W) ∧
      ∃ h : (C • W).c₄ ≠ 0, valAdd p (Units.mk0 ((C • W).c₄) h) = 0

/-- ★★★★★★★★**半安定なら `v_p(Δ_min) ≤ log⁺|j|_p`**。

原文 (GenEll p.18):
> First, observe that if v is any local height of EL, then d · deg∞([EL]) ≥

★乗法還元では `v_p(c₄') = 0` なので `v_p(j) = −v_p(Δ_min)`、
したがって `log⁺|j|_p = v_p(Δ_min)` で**等号**になる。 -/
theorem minDeltaExp_le_maxJ (p : HeightOneSpectrum (𝓞 L)) (W : WeierstrassCurve L)
    [W.IsElliptic] (hss : SemistableAt p W) :
    minDeltaExp p W ≤ max 0 (-jExp p W) := by
  rcases hss with h0 | ⟨C, hC, hc4, hv⟩
  · rw [h0]; exact le_max_left _ _
  · have hΔ : W.Δ ≠ 0 := W.isUnit_Δ.ne_zero
    haveI : (C • W).IsElliptic := inferInstance
    have hjC : (C • W).j = W.j := WeierstrassCurve.variableChange_j W C
    have hΔC : (C • W).Δ ≠ 0 := variableChange_Delta_ne_zero W hΔ C
    have hmd := minDeltaExp_eq p W hΔ C hC
    have hmul : (C • W).j * (C • W).Δ = (C • W).c₄ ^ 3 := by
      rw [WeierstrassCurve.j, ← WeierstrassCurve.coe_Δ', Units.val_inv_eq_inv_val]
      have : ((C • W).Δ' : L) ≠ 0 := ((C • W).Δ').ne_zero
      field_simp
    have hj : W.j ≠ 0 := by
      intro h
      rw [hjC, h, zero_mul] at hmul
      exact hc4 (pow_eq_zero_iff (by norm_num) |>.1 hmul.symm)
    have hunits : Units.mk0 ((C • W).j) (by rw [hjC]; exact hj) * Units.mk0 ((C • W).Δ) hΔC
        = (Units.mk0 ((C • W).c₄) hc4) ^ 3 := by
      refine Units.ext ?_
      push_cast
      exact hmul
    have hval := congrArg (valAdd p) hunits
    rw [valAdd_mul, valAdd_pow, hv] at hval
    have hjeq : Units.mk0 W.j hj = Units.mk0 ((C • W).j) (by rw [hjC]; exact hj) := by
      refine Units.ext ?_
      simp [hjC]
    have hjE : jExp p W = -minDeltaExp p W := by
      rw [jExp, dif_neg hj, hjeq, hmd]
      omega
    rw [hjE, neg_neg]
    exact le_max_right _ _

/-! ## ★★★★★★★★★★大域の形 -/

/-- ★★★★★★★★★★**すべての素点で半安定なら `deg∞ ≤ h_fin(j)`**。

原文 (GenEll p.18):
> First, observe that if v is any local height of EL, then d · deg∞([EL]) ≥

★項ごとの `minDeltaExp_le_maxJ` を足すだけである。 -/
theorem degInfOf_le_htFinJ (W : WeierstrassCurve L) [W.IsElliptic]
    (hss : ∀ p, SemistableAt p W) : degInfOf L W ≤ htFinJ L W := by
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
        (minDeltaExp p W : ℝ) * Real.log (Ideal.absNorm p.asIdeal))
      ≤ ∑ᶠ p : HeightOneSpectrum (𝓞 L),
        ((max 0 (-jExp p W) : ℤ) : ℝ) * Real.log (Ideal.absNorm p.asIdeal) := by
    refine finsum_le_finsum' hg hf (fun p => ?_)
    refine mul_le_mul_of_nonneg_right ?_ (log_absNorm_nonneg p)
    exact Int.cast_le.2 (minDeltaExp_le_maxJ p W (hss p))
  rw [htFinJ, degInfOf]
  have hd : (0:ℝ) ≤ (Module.finrank ℚ L : ℝ) := by positivity
  gcongr

/-- ★★★★★★★★★★★★**半安定なら `deg∞ = h_fin(j)`**——両向きが取れた。

★`≥` は `Found/GaloisRep/HtFinJ.lean` の `htFinJ_le_degInfOf`（**無条件**）。 -/
theorem htFinJ_eq_degInfOf (W : WeierstrassCurve L) [W.IsElliptic]
    (hss : ∀ p, SemistableAt p W) : htFinJ L W = degInfOf L W :=
  le_antisymm (htFinJ_le_degInfOf W) (degInfOf_le_htFinJ W hss)

/-! ## ★★★★★★★★★★★★★★★★★★★★★★★★★★★★半安定曲線での `Prop 3.4` の鎖 -/

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★**半安定曲線での原文の鎖**

    `deg_∞ ≲ ht_∞ ≲ 12(1+ϵ)·ht^Falt ≲ (1+ϵ)·ht_∞`

を `ht_∞ ≔ htJ`（`j` の Weil 高さ）で。

原文 (GenEll p.18):
> First, observe that if v is any local height of EL, then d · deg∞([EL]) ≥

★★★仮定は**素点ごとの半安定性だけ**である
——原文が `Proposition 3.4` の文脈で置いている仮定そのもの。 -/
theorem prop_3_4_chain_semistable (eps : ℝ) (heps : 0 < eps) :
    ∃ C : ℝ, ∀ (L : Type) [Field L] [NumberField L] (E : WeierstrassCurve L) [E.IsElliptic],
      (∀ p, SemistableAt p E) →
        degInfOf L E ≤ htJ L E
        ∧ htJ L E ≤ 12 * (1 + eps) * htFaltOf L E + C
        ∧ 12 * (1 + eps) * htFaltOf L E ≤ (1 + eps) * htJ L E + C := by
  obtain ⟨C, hC⟩ := prop_3_4_chain eps heps
  exact ⟨C, fun L _ _ E _ hss => hC L E (degInfOf_le_htFinJ E hss)⟩

/-! ## ★出典の紐付け(`.src`) -/

def SemistableAt.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 18,
    item := "Proposition 3.4(素点ごとの半安定性——valAdd の言葉で)",
    sectionId := "genell-prop-3-4" }

def htFinJ_eq_degInfOf.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 18,
    item := "Proposition 3.4(半安定なら deg∞ = h_fin(j))",
    sectionId := "genell-prop-3-4" }

def prop_3_4_chain_semistable.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 18,
    item := "Proposition 3.4(3 本の ≲ ——半安定曲線で、ht∞ = htJ として)",
    sectionId := "genell-prop-3-4" }

def prop_3_4_chain_semistable.needs : List ABC3.Meta.ProofObligation :=
  [ .implicitStep
      ("☆**mathlib の還元の言葉との橋**: " ++
       "WeierstrassCurve.HasMultiplicativeReduction は " ++
       "`valuation K (maximalIdeal R) W.c₄ = 1` で乗法還元を定める。" ++
       "★本ファイルの SemistableAt の第 2 の選言肢と同じ内容だが、" ++
       "primeSubring p の DVR 付値と p.valuation L を繋ぐ橋がまだ無い。" ++
       "★★橋が架かれば SemistableAt ⟺ ∀ p, HasSemistableReduction が言える" ++
       "(Found/GaloisRep/Semistable.lean が受け皿)") 6,
    .implicitStep
      ("★★★★★★★★★★★★★★到達点(2026-08-29、第 559): " ++
       "HtJLower.lean の prop_3_4_chain が受けていた hss : deg∞ ≤ h_fin(j) を" ++
       "**素点ごとの半安定性から証明した**。" ++
       "★逆向きは §9-1003 で無条件に取れているので、半安定なら **deg∞ = h_fin(j)**。" ++
       "★★これで原文の鎖 deg∞ ≲ ht∞ ≲ 12(1+ϵ)ht^Falt ≲ (1+ϵ)ht∞ は、" ++
       "ht∞ = htJ として**半安定性だけを仮定に**証明された") 9,
    .implicitStep
      ("☆残るのは ht∞ の**同定**である: Interface/GenEll/EllModuli.lean の " ++
       "EllClass(M_ell(ℚ̄) の点)を構成し、htInf = htJ と置くこと。" ++
       "★northcott 欄の中身は Found/GaloisRep/NorthcottHtJ.lean、" ++
       "degInf_le_htInf・htInf_bdeq_faltings の中身は本ファイルと HtJLower.lean にある") 8 ]

end ABC3.Found.GaloisRep
