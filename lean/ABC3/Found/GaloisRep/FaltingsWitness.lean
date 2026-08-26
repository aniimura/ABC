import ABC3.Found.GaloisRep.DegInfLocal
import ABC3.Found.GaloisRep.NeronWitness
import ABC3.Found.GenEll.ArchNormTotal

/-!
# Galois (G8) 第 357 ブロック —— **★★★★★★★★★★G8 の witness(アルキメデス項つき)**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

## ★★★★★★★★★★到達点——**欠陥 #6 が塞がった**

> **`FaltingsHeightData.nonvacuous`**——(G8) の欄がすべて埋まり、
> **`htFalt` は Faltings 高さの式そのもの**になった

    12·ht^Falt(E) = deg∞(E) − (1/d)·Σ_{σ:L↪ℂ} log( (2π)¹²·‖Δ‖_arch(E^σ) )

★`d = [L:ℚ]`。★★有限素点側が `deg∞`(極小判別式の指数の和)、
アルキメデス側が第 2 項である。

## ★★★★★★★★★2026-08-26 —— 以前の witness は `deg∞/12` だった

第 329 の witness は `htFalt := deg∞/12` で、界面が課す 2 条件
(`htFalt_variableChange`・`prop_3_4`)を満たすだけのものだった。
★★`Proposition 3.4` の内容が**恒等的に成り立つ形**で埋まってしまう
——界面の欠陥 #6(初めての「弱すぎる」型)である。

★★★★塞ぐのに要ったのは**アルキメデス素点での計量**であり、
そこに至る鎖(一意化 → 共体積 → `(g₂,g₃)` が束を決める → `archNorm`)は
第 332-356 の 25 ブロックで閉じた。

★★★★★★界面には 4 つの欄・条件を足した:
`archNorm`・`archNorm_pos`・**`archNorm_eq`**(モジュラー判別式で同定)・
**`htFalt_eq`**(Faltings 高さの式)。
★`archNorm_eq` は `j` の全射性(第 348)により `archNorm` を**一意に決める**。

## ★★`prop_3_4` を保てた理由

`htFalt = deg∞/12 − S/(12d)` で `S := Σ_σ log((2π)¹²·archNorm)` である。
★`archNorm` は**一様に有界**(第 355)なので `S ≤ d·log((2π)¹²M)`、
したがって `S/(12d)` は普遍定数で抑えられ、`C := log((2π)¹²M)/12` で `prop_3_4` が出る。
★★`deg∞ ≥ 0` と `1+ε > 1` から `deg∞/(12(1+ε)) ≤ deg∞/12` は従来どおり。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `minDeltaExp_variableChange` | ★★★★★極小判別式の指数は変数変換で不変 |
| `degInfOf_variableChange` | ★★★★★★`deg∞` は変数変換で不変 |
| `archSum`・`archSum_variableChange` | ★★★★アルキメデス和とその不変性 |
| `archSum_le` | ★★★★★**一様な上界** |
| `htFaltOf` | ★★★★★★★**Faltings 高さの式** |
| `faltingsHeightDataWitness` | ★★★★★★★★**`FaltingsHeightData` の実装** |
| `FaltingsHeightData.nonvacuous` | ★★★★★★★★★★**G8 の欄が埋まる** |
-/

namespace ABC3.Found.GaloisRep

open IsDedekindDomain NumberField WeierstrassCurve ABC3.Interface.GaloisRep ABC3.Found.GenEll
open scoped Classical

variable {L : Type} [Field L] [NumberField L]

/-! ## ★★★★★★`deg∞` は変数変換で不変(真の定理) -/

set_option maxHeartbeats 1600000 in
/-- ★★★★★**極小判別式の指数は変数変換で不変である**。

★`Δ(C•W) = u⁻¹²Δ` で付値は `−12v(u)`、Néron 指数は `−v(u)` だけ動くので、
`v(Δ) − 12·(Néron 指数)` は打ち消し合う。 -/
theorem minDeltaExp_variableChange (p : HeightOneSpectrum (𝓞 L)) (W : WeierstrassCurve L)
    (C : VariableChange L) : minDeltaExp p (C • W) = minDeltaExp p W := by
  by_cases hΔ : W.Δ = 0
  · have h2 : (C • W).Δ = 0 := by
      rw [WeierstrassCurve.variableChange_Δ, hΔ, mul_zero]
    rw [minDeltaExp, dif_pos h2, minDeltaExp, dif_pos hΔ]
  · have h2 : (C • W).Δ ≠ 0 := variableChange_Delta_ne_zero W hΔ C
    rw [minDeltaExp, dif_neg h2, minDeltaExp, dif_neg hΔ, neronExp_variableChange p W hΔ C]
    have hu : (Units.mk0 ((C • W).Δ) h2) = C.u⁻¹ ^ 12 * Units.mk0 W.Δ hΔ := by
      refine Units.ext ?_
      show (C • W).Δ = _
      rw [WeierstrassCurve.variableChange_Δ]
      push_cast
      simp
    rw [hu, valAdd_mul, valAdd_pow, valAdd_inv]
    omega

/-- ★★★★★★**`deg∞` は変数変換で不変である**。 -/
theorem degInfOf_variableChange (E : WeierstrassCurve L) (C : VariableChange L) :
    degInfOf L (C • E) = degInfOf L E := by
  rw [degInfOf, degInfOf]
  congr 1
  refine finsum_congr (fun q => ?_)
  rw [minDeltaExp_variableChange]

/-! ## ★★★★アルキメデス和 -/

/-- ★★★★**アルキメデス和** `Σ_σ log((2π)¹²·‖Δ‖_arch(E^σ))`。 -/
noncomputable def archSum (L : Type) [Field L] [NumberField L] (E : WeierstrassCurve L) : ℝ :=
  ∑ σ : (L →+* ℂ), Real.log ((2 * Real.pi) ^ 12 * ABC3.Found.GenEll.archNorm E σ)

theorem archSum_variableChange (L : Type) [Field L] [NumberField L] (E : WeierstrassCurve L)
    (C : VariableChange L) : archSum L (C • E) = archSum L E := by
  rw [archSum, archSum]
  refine Finset.sum_congr rfl (fun σ _ => ?_)
  rw [ABC3.Found.GenEll.archNorm_variableChange']

/-- ★★★★★**アルキメデス和の一様な上界**——`prop_3_4` を保つのに要る。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

★埋め込みの個数は `[L:ℚ]` なので、項ごとの上界を `d` 倍すればよい。 -/
theorem archSum_le (L : Type) [Field L] [NumberField L] (E : WeierstrassCurve L)
    (M : ℝ) (hM : 1 ≤ M) (hMb : ∀ W : WeierstrassCurve ℂ, curveArchInv W ≤ M) :
    archSum L E ≤ (Module.finrank ℚ L : ℝ) * Real.log ((2 * Real.pi) ^ 12 * M) := by
  have hpi : (1:ℝ) ≤ (2 * Real.pi) ^ 12 := by
    have h2 : (1:ℝ) ≤ 2 * Real.pi := by nlinarith [Real.two_le_pi]
    exact one_le_pow₀ h2
  have hbig : (1:ℝ) ≤ (2 * Real.pi) ^ 12 * M := by nlinarith
  have hterm : ∀ σ : (L →+* ℂ),
      Real.log ((2 * Real.pi) ^ 12 * ABC3.Found.GenEll.archNorm E σ)
        ≤ Real.log ((2 * Real.pi) ^ 12 * M) := by
    intro σ
    rcases eq_or_lt_of_le (ABC3.Found.GenEll.archNorm_nonneg E σ) with h0 | h0
    · rw [← h0, mul_zero, Real.log_zero]
      exact Real.log_nonneg hbig
    · refine Real.log_le_log (by positivity) ?_
      have hb := hMb (E.map σ)
      exact mul_le_mul_of_nonneg_left hb (by positivity)
  calc archSum L E ≤ ∑ _σ : (L →+* ℂ), Real.log ((2 * Real.pi) ^ 12 * M) :=
        Finset.sum_le_sum (fun σ _ => hterm σ)
    _ = (Fintype.card (L →+* ℂ) : ℝ) * Real.log ((2 * Real.pi) ^ 12 * M) := by
        rw [Finset.sum_const, Finset.card_univ, nsmul_eq_mul]
    _ = (Module.finrank ℚ L : ℝ) * Real.log ((2 * Real.pi) ^ 12 * M) := by
        rw [NumberField.Embeddings.card L ℂ]

/-! ## ★★★★★★★Faltings 高さ -/

/-- ★★★★★★★**Faltings 高さ** `12·ht^Falt = deg∞ − (1/d)·Σ_σ log((2π)¹²‖Δ‖_arch)`。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

★有限素点側が `deg∞`、アルキメデス側が第 2 項である。 -/
noncomputable def htFaltOf (L : Type) [Field L] [NumberField L] (E : WeierstrassCurve L) : ℝ :=
  degInfOf L E / 12 - archSum L E / (12 * (Module.finrank ℚ L : ℝ))

/-! ## ★★★★★★★★G8 の witness -/

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★★**`FaltingsHeightData` の実装**(アルキメデス項つき)。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

★★★★★★2026-08-26: `htFalt` は **Faltings 高さの式**になった
——`deg∞/12` だった第 329 の witness を置き換える(界面の欠陥 #6 を塞いだ)。 -/
noncomputable def faltingsHeightDataWitness : FaltingsHeightData where
  toSemistableModelData := semistableModelDataWitness
  htFalt := fun L _ _ E => htFaltOf L E
  htFalt_variableChange := by
    intro L _ _ E C
    rw [htFaltOf, htFaltOf, degInfOf_variableChange, archSum_variableChange]
  degInf := fun L _ _ E => degInfOf L E
  degInf_nonneg := by
    intro L _ _ E
    exact degInfOf_nonneg E
  degInf_ge_localHeight := by
    intro L _ _ E Lv _ _ R _ _ _ _ _ _ _ _ h p hp
    exact degInfOf_ge_localHeight E h p hp
  archNorm := fun L _ _ E σ => ABC3.Found.GenEll.archNorm E σ
  archNorm_pos := by
    intro L _ _ E _ σ
    exact ABC3.Found.GenEll.archNorm_pos E σ
  archNorm_eq := by
    intro L _ _ E _ σ τ h
    haveI : (latticeCurve (tauPair τ)).IsElliptic := isElliptic_latticeCurve' _
    have hj : (E.map σ).j = (latticeCurve (tauPair τ)).j := by
      rw [← h]; exact jFun_eq_latticeCurve_j τ
    show ABC3.Found.GenEll.archNorm E σ = _
    rw [ABC3.Found.GenEll.archNorm,
      curveArchInv_congr_j (inferInstance : (E.map σ).IsElliptic) inferInstance hj,
      curveArchInv_tauPair, peterssonDelta]
  htFalt_eq := by
    intro L _ _ E _
    show 12 * htFaltOf L E = degInfOf L E - archSum L E / (Module.finrank ℚ L : ℝ)
    rw [htFaltOf]
    ring
  prop_3_4 := by
    intro ε hε
    obtain ⟨M, hM1, hMb⟩ := ABC3.Found.GenEll.exists_bound_curveArchInv'
    refine ⟨Real.log ((2 * Real.pi) ^ 12 * M) / 12, ?_⟩
    intro L _ _ E _
    have hd : (0:ℝ) < (Module.finrank ℚ L : ℝ) := by exact_mod_cast Module.finrank_pos
    have h0 : 0 ≤ degInfOf L E := degInfOf_nonneg E
    have hS := archSum_le L E M hM1 hMb
    show degInfOf L E / (12 * (1 + ε))
      ≤ htFaltOf L E + Real.log ((2 * Real.pi) ^ 12 * M) / 12
    rw [htFaltOf]
    have hle1 : degInfOf L E / (12 * (1 + ε)) ≤ degInfOf L E / 12 := by
      gcongr
      nlinarith
    have hle2 : archSum L E / (12 * (Module.finrank ℚ L : ℝ))
        ≤ Real.log ((2 * Real.pi) ^ 12 * M) / 12 := by
      rw [show Real.log ((2 * Real.pi) ^ 12 * M) / 12
          = ((Module.finrank ℚ L : ℝ) * Real.log ((2 * Real.pi) ^ 12 * M))
            / (12 * (Module.finrank ℚ L : ℝ)) from by field_simp]
      gcongr
    linarith

/-- ★★★★★★★★★★**`FaltingsHeightData` は非空虚である**——G8 の欄が埋まる。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

★★★★★★★2026-08-26: `htFalt` は **Faltings 高さの式**であり、
界面の `htFalt_eq` により `degInf` と `archNorm` から一意に決まる
(`Check/GaloisRep/HtFaltPinned.lean`)。 -/
theorem FaltingsHeightData.nonvacuous : Nonempty FaltingsHeightData :=
  ⟨faltingsHeightDataWitness⟩

/-! ## ★出典の紐付け(`.src`) -/

def degInfOf_variableChange.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 18,
    item := "Proposition 3.4(無限遠因子の次数 deg∞)",
    sectionId := "genell-prop-3-4" }

def htFaltOf.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(Faltings Heights and the Divisor at Infinity)",
    sectionId := "genell-prop-3-4" }

def faltingsHeightDataWitness.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(Faltings Heights and the Divisor at Infinity)",
    sectionId := "genell-prop-3-4" }

end ABC3.Found.GaloisRep
